#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

u64 pow_u64(u64 p, u64 e) {
    u128 result = 1;

    for (u64 i = 0; i < e; ++i) {
        result *= p;
    }

    return static_cast<u64>(result);
}

std::vector<u64> base_p_digits(
    u64 p,
    u64 n
) {
    std::vector<u64> digits;

    while (n > 0) {
        digits.push_back(n % p);
        n /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

u64 digit_product(
    u64 p,
    u64 n
) {
    const auto digits =
        base_p_digits(p, n);

    u64 product = 1;

    for (const u64 d : digits) {
        product *= d + 1;
    }

    return product;
}

std::vector<u64> interval_counts(
    u64 p,
    u64 q
) {
    const auto digits =
        base_p_digits(p, q);

    std::vector<u64> counts(
        digits.size(),
        0
    );

    u64 product_higher = 1;

    for (std::size_t r = digits.size();
         r-- > 0;) {

        counts[r] =
            digits[r] *
            product_higher;

        product_higher *=
            digits[r] + 1;
    }

    return counts;
}

/*
 * Compute

 *   S(q)
 *
 * = sum C_r *
 *   (p^r - (q mod p^r)).
 */
u64 inner_sum(
    u64 p,
    u64 q
) {
    const auto counts =
        interval_counts(p, q);

    u128 total = 0;

    for (std::size_t r = 0;
         r < counts.size();
         ++r) {

        const u64 pr =
            pow_u64(
                p,
                static_cast<u64>(r)
            );

        const u64 qlow =
            r == 0
            ? 0
            : q % pr;

        total +=
            static_cast<u128>(
                counts[r]
            ) *
            static_cast<u128>(
                pr - qlow
            );
    }

    return static_cast<u64>(total);
}

/*
 * Full weighted interval sum:

 * H =
 * sum C_r *
 * [
 *   (p^r - q_low)*p^e - s0
 * ]
 */
u64 full_sum(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const auto counts =
        interval_counts(
            p,
            q
        );

    const u64 pe =
        pow_u64(
            p,
            e
        );

    u128 total = 0;

    for (std::size_t r = 0;
         r < counts.size();
         ++r) {

        const u64 pr =
            pow_u64(
                p,
                static_cast<u64>(r)
            );

        const u64 qlow =
            r == 0
            ? 0
            : q % pr;

        const u64 length =
            static_cast<u64>(
                static_cast<u128>(
                    pr - qlow
                ) *
                static_cast<u128>(pe) -
                static_cast<u128>(s0)
            );

        total +=
            static_cast<u128>(
                counts[r]
            ) *
            static_cast<u128>(length);
    }

    return static_cast<u64>(total);
}

/*
 * Simplified formula:

 *   H =
 *   q*p^e
 *   - s0*(prod(q_i+1)-1)
 */
u64 simplified_sum(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const u64 pe =
        pow_u64(p, e);

    const u64 count =
        digit_product(p, q) - 1;

    const u128 result =
        static_cast<u128>(q) *
        static_cast<u128>(pe) -
        static_cast<u128>(s0) *
        static_cast<u128>(count);

    return static_cast<u64>(result);
}

/*
 * Split

 *   q = Q + a*p^R

 * where a is the most significant base-p digit.
 */
bool split_q(
    u64 p,
    u64 q,
    u64 &Q,
    u64 &a,
    u64 &pR
) {
    if (q == 0) {
        return false;
    }

    u64 temp = q;
    u64 R = 0;

    while (temp >= p) {
        temp /= p;
        ++R;
    }

    pR =
        pow_u64(
            p,
            R
        );

    a = q / pR;
    Q = q % pR;

    return true;
}

/*
 * Candidate recurrence for the full weighted sum.

 * Since

 *   H(q) = p^e*S(q) - s0*N(q)

 * with

 *   S(q)=q
 *   N(q)=prod(q_i+1)-1,

 * and

 *   N(Q+a*p^R)
 *       = (a+1)(N(Q)+1)-1,

 * the recurrence becomes

 *   H(q)
 *     =
 *     (a+1)H(Q)
 *     +
 *     a*p^e*(p^R-Q)
 *     -
 *     a*s0.
 */
u64 recurrence_value(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    if (q == 0) {
        return 0;
    }

    u64 Q = 0;
    u64 a = 0;
    u64 pR = 0;

    split_q(
        p,
        q,
        Q,
        a,
        pR
    );

    const u64 lower =
        full_sum(
            p,
            e,
            s0,
            Q
        );

    const u64 pe =
        pow_u64(
            p,
            e
        );

    const u128 result =
        static_cast<u128>(a + 1) *
        static_cast<u128>(lower) +

        static_cast<u128>(a) *
        static_cast<u128>(pe) *
        static_cast<u128>(pR - Q) -

        static_cast<u128>(a) *
        static_cast<u128>(s0);

    return static_cast<u64>(result);
}

void run_small_exhaustive(
    u64 max_q,
    u64 &cases,
    u64 &recurrence_pass,
    u64 &simplified_pass,
    u64 &agreement_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (const u64 p : primes) {
        for (u64 q = 1;
             q <= max_q;
             ++q) {

            /*
             * Pick several e,s0 values. The recurrence itself
             * does not depend on how q was originally embedded
             * into m.
             */
            const std::vector<u64> e_values = {
                1, 2, 3
            };

            const std::vector<u64> s0_values = {
                1,
                p,
                p * p
            };

            for (const u64 e : e_values) {
                for (const u64 s0 : s0_values) {

                    ++cases;

                    const u64 direct =
                        full_sum(
                            p,
                            e,
                            s0,
                            q
                        );

                    const u64 recurrence =
                        recurrence_value(
                            p,
                            e,
                            s0,
                            q
                        );

                    const u64 simplified =
                        simplified_sum(
                            p,
                            e,
                            s0,
                            q
                        );

                    if (direct ==
                        recurrence) {
                        ++recurrence_pass;
                    }

                    if (direct ==
                        simplified) {
                        ++simplified_pass;
                    }

                    if (recurrence ==
                            simplified &&
                        direct ==
                            simplified) {
                        ++agreement_pass;
                    }
                }
            }
        }
    }
}

void run_large_random(
    u64 trials,
    u64 &recurrence_pass,
    u64 &simplified_pass,
    u64 &agreement_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    std::mt19937_64 rng(
        0x19320260913ULL
    );

    std::uniform_int_distribution<u64> q_dist(
        1ULL,
        1000000000000000000ULL
    );

    std::uniform_int_distribution<u64> e_dist(
        1,
        10
    );

    std::uniform_int_distribution<u64> s_dist(
        1ULL,
        1000000ULL
    );

    for (u64 i = 0;
         i < trials;
         ++i) {

        const u64 p =
            primes[i % primes.size()];

        const u64 q =
            q_dist(rng);

        const u64 e =
            e_dist(rng);

        const u64 s0 =
            s_dist(rng);

        const u64 direct =
            full_sum(
                p,
                e,
                s0,
                q
            );

        const u64 recurrence =
            recurrence_value(
                p,
                e,
                s0,
                q
            );

        const u64 simplified =
            simplified_sum(
                p,
                e,
                s0,
                q
            );

        if (direct ==
            recurrence) {
            ++recurrence_pass;
        }

        if (direct ==
            simplified) {
            ++simplified_pass;
        }

        if (direct ==
                simplified &&
            recurrence ==
                simplified) {
            ++agreement_pass;
        }
    }
}

void print_example(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    u64 Q = 0;
    u64 a = 0;
    u64 pR = 0;

    split_q(
        p,
        q,
        Q,
        a,
        pR
    );

    const u64 direct =
        full_sum(
            p,
            e,
            s0,
            q
        );

    const u64 recurrence =
        recurrence_value(
            p,
            e,
            s0,
            q
        );

    const u64 simplified =
        simplified_sum(
            p,
            e,
            s0,
            q
        );

    std::cout
        << "p=" << p
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
        << "\n";

    std::cout
        << "  Q=" << Q
        << " a=" << a
        << " pR=" << pR
        << "\n";

    std::cout
        << "  direct="
        << direct
        << "\n";

    std::cout
        << "  recurrence="
        << recurrence
        << "\n";

    std::cout
        << "  simplified="
        << simplified
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 193\n";

    u64 cases = 0;
    u64 recurrence_pass = 0;
    u64 simplified_pass = 0;
    u64 agreement_pass = 0;

    run_small_exhaustive(
        5000,
        cases,
        recurrence_pass,
        simplified_pass,
        agreement_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "cases="
        << cases
        << "\n";

    std::cout
        << "full_recurrence_pass="
        << recurrence_pass
        << "/" << cases
        << "\n";

    std::cout
        << "simplified_formula_pass="
        << simplified_pass
        << "/" << cases
        << "\n";

    std::cout
        << "full_agreement_pass="
        << agreement_pass
        << "/" << cases
        << "\n";

    u64 large_recurrence_pass = 0;
    u64 large_simplified_pass = 0;
    u64 large_agreement_pass = 0;

    run_large_random(
        200000,
        large_recurrence_pass,
        large_simplified_pass,
        large_agreement_pass
    );

    std::cout
        << "\nLARGE RANDOM\n";

    std::cout
        << "full_recurrence_pass="
        << large_recurrence_pass
        << "/200000\n";

    std::cout
        << "simplified_formula_pass="
        << large_simplified_pass
        << "/200000\n";

    std::cout
        << "full_agreement_pass="
        << large_agreement_pass
        << "/200000\n";

    std::cout
        << "\nEXAMPLES\n";

    print_example(
        2,
        2,
        1,
        5
    );

    print_example(
        2,
        3,
        1,
        5
    );

    print_example(
        3,
        2,
        6,
        4
    );

    print_example(
        5,
        2,
        20,
        7
    );

    std::cout
        << "\nRECURRENCE UNDER TEST\n";

    std::cout
        << "H(Q+a*p^R) = "
        << "(a+1)H(Q)"
        << " + a*p^e*(p^R-Q)"
        << " - a*s0\n";

    std::cout
        << "H(q) = "
        << "q*p^e "
        << "- s0*(product(q_i+1)-1)\n";

    std::cout
        << "FINISHED EXPERIMENT 193\n";

    return 0;
}
