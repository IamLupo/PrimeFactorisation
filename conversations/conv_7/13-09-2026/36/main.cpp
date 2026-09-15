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

/*
 * Return the base-p digits of n, least-significant first.
 */
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

/*
 * C_r =
 *
 *   q_r * product_{i>r}(q_i+1)
 */
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
 * q mod p^r.
 */
u64 q_low_part(
    u64 p,
    u64 q,
    u64 r
) {
    if (r == 0) {
        return 0;
    }

    return q % pow_u64(p, r);
}

/*
 * The experimentally observed inner identity:

 *
 *   sum_r C_r *
 *       (p^r - (q mod p^r))
 *
 *   = q
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

        if (counts[r] == 0) {
            continue;
        }

        const u64 pr =
            pow_u64(
                p,
                static_cast<u64>(r)
            );

        const u64 qlow =
            q_low_part(
                p,
                q,
                static_cast<u64>(r)
            );

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
 * Sum C_r:

 *   product_i(q_i+1) - 1
 */
u64 interval_count_sum(
    u64 p,
    u64 q
) {
    const auto counts =
        interval_counts(
            p,
            q
        );

    u64 total = 0;

    for (const u64 c : counts) {
        total += c;
    }

    return total;
}

u64 digit_product(
    u64 p,
    u64 n
) {
    const auto digits =
        base_p_digits(
            p,
            n
        );

    u64 product = 1;

    for (const u64 d : digits) {
        product *= d + 1;
    }

    return product;
}

/*
 * Derive the structural parameters from m.

 * m = s0 + q*p^e - 1.
 */
bool derive_parameters(
    u64 p,
    u64 m,
    u64 &e,
    u64 &s0,
    u64 &q
) {
    u64 x = m;
    u64 a = 0;

    while (x > 0 && x % p == p - 1) {
        ++a;
        x /= p;
    }

    if (x == 0) {
        return false;
    }

    const u64 b = x % p;

    x /= p;

    if (b >= p - 1) {
        return false;
    }

    u64 z = 0;

    while (x > 0 && x % p == 0) {
        ++z;
        x /= p;
    }

    e = a + 1 + z;

    s0 =
        static_cast<u64>(
            static_cast<u128>(b + 1) *
            static_cast<u128>(
                pow_u64(p, a)
            )
        );

    const u64 pe =
        pow_u64(p, e);

    const u128 numerator =
        static_cast<u128>(m) +
        1 -
        static_cast<u128>(s0);

    if (numerator == 0 ||
        numerator % pe != 0) {
        return false;
    }

    q =
        static_cast<u64>(
            numerator / pe
        );

    return q > 0;
}

/*
 * Complete weighted interval sum.

 * H =
 *
 *   sum_r C_r *
 *       [(p^r - q_low_r)*p^e - s0]
 */
u64 weighted_interval_sum(
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
        pow_u64(p, e);

    u128 total = 0;

    for (std::size_t r = 0;
         r < counts.size();
         ++r) {

        if (counts[r] == 0) {
            continue;
        }

        const u64 pr =
            pow_u64(
                p,
                static_cast<u64>(r)
            );

        const u64 qlow =
            q_low_part(
                p,
                q,
                static_cast<u64>(r)
            );

        const u64 length =
            static_cast<u64>(
                static_cast<u128>(
                    pr - qlow
                ) *
                static_cast<u128>(
                    pe
                ) -
                static_cast<u128>(
                    s0
                )
            );

        total +=
            static_cast<u128>(
                counts[r]
            ) *
            static_cast<u128>(
                length
            );
    }

    return static_cast<u64>(total);
}

/*
 * Simplified symbolic expression:

 *   q*p^e
 *   - s0*(product(q_i+1)-1)
 */
u64 simplified_formula(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const u64 pe =
        pow_u64(p, e);

    const u64 product_q =
        digit_product(
            p,
            q
        );

    const u128 value =
        static_cast<u128>(q) *
        static_cast<u128>(pe) -
        static_cast<u128>(s0) *
        static_cast<u128>(
            product_q - 1
        );

    return static_cast<u64>(value);
}

/*
 * Direct Lucas digit count:

 *   H = m+1 - product_i(m_i+1)
 */
u64 direct_digit_formula(
    u64 p,
    u64 m
) {
    return
        m + 1 -
        digit_product(
            p,
            m
        );
}

/*
 * Test the inner q-identity exhaustively.
 */
void run_q_exhaustive(
    u64 max_q,
    u64 &cases,
    u64 &inner_pass,
    u64 &count_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (const u64 p : primes) {
        for (u64 q = 1;
             q <= max_q;
             ++q) {

            ++cases;

            if (inner_sum(
                    p,
                    q
                ) == q) {
                ++inner_pass;
            }

            const u64 expected_count =
                digit_product(
                    p,
                    q
                ) - 1;

            if (interval_count_sum(
                    p,
                    q
                ) == expected_count) {
                ++count_pass;
            }
        }
    }
}

/*
 * Test complete structural collapse for all admissible
 * small m.
 */
void run_m_exhaustive(
    u64 max_m,
    u64 &cases,
    u64 &simplification_pass,
    u64 &digit_formula_pass,
    u64 &full_agreement_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (const u64 p : primes) {
        for (u64 m = 2;
             m <= max_m;
             ++m) {

            u64 e = 0;
            u64 s0 = 0;
            u64 q = 0;

            if (!derive_parameters(
                    p,
                    m,
                    e,
                    s0,
                    q
                )) {
                continue;
            }

            ++cases;

            const u64 weighted =
                weighted_interval_sum(
                    p,
                    e,
                    s0,
                    q
                );

            const u64 simplified =
                simplified_formula(
                    p,
                    e,
                    s0,
                    q
                );

            const u64 digit_formula =
                direct_digit_formula(
                    p,
                    m
                );

            if (weighted ==
                simplified) {
                ++simplification_pass;
            }

            if (weighted ==
                digit_formula) {
                ++digit_formula_pass;
            }

            if (simplified ==
                    digit_formula &&
                weighted ==
                    digit_formula) {
                ++full_agreement_pass;
            }
        }
    }
}

/*
 * Large random q test.

 * This does not need m at all. It attacks the
 * inner digit identity directly.
 */
void run_q_random(
    u64 trials,
    u64 &inner_pass,
    u64 &count_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31
    };

    std::mt19937_64 rng(
        0x19120260913ULL
    );

    std::uniform_int_distribution<u64> dist(
        1ULL,
        1000000000000000000ULL
    );

    for (u64 i = 0;
         i < trials;
         ++i) {

        const u64 p =
            primes[i % primes.size()];

        const u64 q =
            dist(rng);

        if (inner_sum(
                p,
                q
            ) == q) {
            ++inner_pass;
        }

        const u64 expected =
            digit_product(
                p,
                q
            ) - 1;

        if (interval_count_sum(
                p,
                q
            ) == expected) {
            ++count_pass;
        }
    }
}

/*
 * Large random complete m test.
 */
void run_m_random(
    u64 trials,
    u64 &simplification_pass,
    u64 &digit_formula_pass,
    u64 &agreement_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31
    };

    std::mt19937_64 rng(
        0x19120260914ULL
    );

    std::uniform_int_distribution<u64> dist(
        1000000000000ULL,
        1000000000000000000ULL
    );

    for (u64 i = 0;
         i < trials;
         ++i) {

        const u64 p =
            primes[i % primes.size()];

        const u64 m =
            dist(rng);

        u64 e = 0;
        u64 s0 = 0;
        u64 q = 0;

        if (!derive_parameters(
                p,
                m,
                e,
                s0,
                q
            )) {
            --i;
            continue;
        }

        const u64 weighted =
            weighted_interval_sum(
                p,
                e,
                s0,
                q
            );

        const u64 simplified =
            simplified_formula(
                p,
                e,
                s0,
                q
            );

        const u64 digit_formula =
            direct_digit_formula(
                p,
                m
            );

        if (weighted ==
            simplified) {
            ++simplification_pass;
        }

        if (weighted ==
            digit_formula) {
            ++digit_formula_pass;
        }

        if (simplified ==
                digit_formula &&
            weighted ==
                digit_formula) {
            ++agreement_pass;
        }
    }
}

void print_example(
    u64 p,
    u64 m
) {
    u64 e = 0;
    u64 s0 = 0;
    u64 q = 0;

    if (!derive_parameters(
            p,
            m,
            e,
            s0,
            q
        )) {
        return;
    }

    const auto counts =
        interval_counts(
            p,
            q
        );

    const u64 inner =
        inner_sum(
            p,
            q
        );

    const u64 weighted =
        weighted_interval_sum(
            p,
            e,
            s0,
            q
        );

    const u64 simplified =
        simplified_formula(
            p,
            e,
            s0,
            q
        );

    const u64 digit =
        direct_digit_formula(
            p,
            m
        );

    std::cout
        << "p=" << p
        << " m=" << m
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
        << "\n";

    std::cout
        << "  inner_sum="
        << inner
        << "\n";

    std::cout
        << "  q="
        << q
        << "\n";

    std::cout
        << "  sum_C="
        << interval_count_sum(
            p,
            q
        )
        << "\n";

    std::cout
        << "  product(q_i+1)-1="
        << digit_product(
            p,
            q
        ) - 1
        << "\n";

    std::cout
        << "  weighted="
        << weighted
        << "\n";

    std::cout
        << "  simplified="
        << simplified
        << "\n";

    std::cout
        << "  digit_formula="
        << digit
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 191\n";

    u64 q_cases = 0;
    u64 q_inner_pass = 0;
    u64 q_count_pass = 0;

    run_q_exhaustive(
        5000,
        q_cases,
        q_inner_pass,
        q_count_pass
    );

    std::cout
        << "\nQ-EXHAUSTIVE\n";

    std::cout
        << "cases="
        << q_cases
        << "\n";

    std::cout
        << "inner_identity_pass="
        << q_inner_pass
        << "/" << q_cases
        << "\n";

    std::cout
        << "count_identity_pass="
        << q_count_pass
        << "/" << q_cases
        << "\n";

    u64 m_cases = 0;
    u64 simplification_pass = 0;
    u64 digit_formula_pass = 0;
    u64 full_agreement_pass = 0;

    run_m_exhaustive(
        5000,
        m_cases,
        simplification_pass,
        digit_formula_pass,
        full_agreement_pass
    );

    std::cout
        << "\nM-EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << m_cases
        << "\n";

    std::cout
        << "weighted_equals_simplified="
        << simplification_pass
        << "/" << m_cases
        << "\n";

    std::cout
        << "weighted_equals_digit_formula="
        << digit_formula_pass
        << "/" << m_cases
        << "\n";

    std::cout
        << "full_three_way_agreement="
        << full_agreement_pass
        << "/" << m_cases
        << "\n";

    u64 random_q_inner = 0;
    u64 random_q_count = 0;

    run_q_random(
        200000,
        random_q_inner,
        random_q_count
    );

    std::cout
        << "\nQ-LARGE RANDOM\n";

    std::cout
        << "inner_identity_pass="
        << random_q_inner
        << "/200000\n";

    std::cout
        << "count_identity_pass="
        << random_q_count
        << "/200000\n";

    u64 random_m_simplification = 0;
    u64 random_m_digit = 0;
    u64 random_m_agreement = 0;

    run_m_random(
        200000,
        random_m_simplification,
        random_m_digit,
        random_m_agreement
    );

    std::cout
        << "\nM-LARGE RANDOM\n";

    std::cout
        << "weighted_equals_simplified="
        << random_m_simplification
        << "/200000\n";

    std::cout
        << "weighted_equals_digit_formula="
        << random_m_digit
        << "/200000\n";

    std::cout
        << "full_three_way_agreement="
        << random_m_agreement
        << "/200000\n";

    std::cout
        << "\nEXAMPLES\n";

    print_example(2, 20);
    print_example(2, 40);
    print_example(3, 41);
    print_example(5, 194);

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "sum_r C_r * "
        << "(p^r - (q mod p^r)) = q\n";

    std::cout
        << "sum_r C_r = "
        << "product_i(q_i+1) - 1\n";

    std::cout
        << "weighted interval sum = "
        << "q*p^e - "
        << "s0*(product_i(q_i+1)-1)\n";

    std::cout
        << "FINISHED EXPERIMENT 191\n";

    return 0;
}
