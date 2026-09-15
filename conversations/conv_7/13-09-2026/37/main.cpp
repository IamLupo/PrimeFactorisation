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

u64 inner_sum(
    u64 p,
    u64 q
) {
    const auto counts =
        interval_counts(
            p,
            q
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
 * Split

 *   q = Q + a*p^R

 * where a is the most significant base-p digit.
 */
bool split_highest_digit(
    u64 p,
    u64 q,
    u64 &Q,
    u64 &a,
    u64 &pR
) {
    if (q == 0) {
        Q = 0;
        a = 0;
        pR = 1;
        return false;
    }

    std::size_t highest = 0;

    u64 temp = q;

    while (temp >= p) {
        temp /= p;
        ++highest;
    }

    pR =
        pow_u64(
            p,
            static_cast<u64>(highest)
        );

    a = q / pR;
    Q = q % pR;

    return true;
}

/*
 * Recurrence proposed for Experiment 192:

 *   S(q)
 *     = (a+1) S(Q)
 *       + a(p^R-Q)

 * where

 *   q = Q + a*p^R.
 */
u64 recurrence_value(
    u64 p,
    u64 q
) {
    if (q == 0) {
        return 0;
    }

    u64 Q = 0;
    u64 a = 0;
    u64 pR = 0;

    split_highest_digit(
        p,
        q,
        Q,
        a,
        pR
    );

    const u64 lower =
        inner_sum(
            p,
            Q
        );

    const u128 value =
        static_cast<u128>(a + 1) *
        static_cast<u128>(lower) +
        static_cast<u128>(a) *
        static_cast<u128>(
            pR - Q
        );

    return static_cast<u64>(value);
}

/*
 * Recursive version of the same identity.

 * This is included as a structural check:
 *
 *   S(0)=0
 *   S(q)=q
 *
 * follows recursively if the recurrence holds.
 */
u64 recursive_identity_value(
    u64 p,
    u64 q
) {
    if (q == 0) {
        return 0;
    }

    u64 Q = 0;
    u64 a = 0;
    u64 pR = 0;

    split_highest_digit(
        p,
        q,
        Q,
        a,
        pR
    );

    const u64 lower =
        recursive_identity_value(
            p,
            Q
        );

    const u128 value =
        static_cast<u128>(a + 1) *
        static_cast<u128>(lower) +
        static_cast<u128>(a) *
        static_cast<u128>(
            pR - Q
        );

    return static_cast<u64>(value);
}

/*
 * Also test the interval-count recurrence.

 * Let

 *   N(q)=sum_r C_r
 *       = product_i(q_i+1)-1.

 * If q=Q+a*p^R then

 *   N(q)=(a+1)(N(Q)+1)-1.
 */
u64 interval_count(
    u64 p,
    u64 q
) {
    const auto counts =
        interval_counts(
            p,
            q
        );

    u64 result = 0;

    for (const u64 c : counts) {
        result += c;
    }

    return result;
}

u64 interval_count_recurrence(
    u64 p,
    u64 q
) {
    if (q == 0) {
        return 0;
    }

    u64 Q = 0;
    u64 a = 0;
    u64 pR = 0;

    split_highest_digit(
        p,
        q,
        Q,
        a,
        pR
    );

    const u64 lower =
        interval_count(
            p,
            Q
        );

    const u128 value =
        static_cast<u128>(a + 1) *
        static_cast<u128>(lower + 1) -
        1;

    return static_cast<u64>(value);
}

void run_small_exhaustive(
    u64 max_q,
    u64 &cases,
    u64 &recurrence_pass,
    u64 &direct_pass,
    u64 &recursive_pass,
    u64 &count_recurrence_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (const u64 p : primes) {
        for (u64 q = 1;
             q <= max_q;
             ++q) {

            ++cases;

            const u64 direct =
                inner_sum(
                    p,
                    q
                );

            const u64 recurrence =
                recurrence_value(
                    p,
                    q
                );

            const u64 recursive =
                recursive_identity_value(
                    p,
                    q
                );

            const u64 count =
                interval_count(
                    p,
                    q
                );

            const u64 count_recurrence =
                interval_count_recurrence(
                    p,
                    q
                );

            if (direct == recurrence) {
                ++recurrence_pass;
            }

            if (direct == q) {
                ++direct_pass;
            }

            if (recursive == q) {
                ++recursive_pass;
            }

            if (count == count_recurrence) {
                ++count_recurrence_pass;
            }
        }
    }
}

void run_large_random(
    u64 trials,
    u64 &recurrence_pass,
    u64 &direct_pass,
    u64 &count_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    std::mt19937_64 rng(
        0x19220260913ULL
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

        const u64 direct =
            inner_sum(
                p,
                q
            );

        const u64 recurrence =
            recurrence_value(
                p,
                q
            );

        const u64 count =
            interval_count(
                p,
                q
            );

        const u64 count_recurrence =
            interval_count_recurrence(
                p,
                q
            );

        if (direct == recurrence) {
            ++recurrence_pass;
        }

        if (direct == q) {
            ++direct_pass;
        }

        if (count == count_recurrence) {
            ++count_pass;
        }
    }
}

void print_example(
    u64 p,
    u64 q
) {
    u64 Q = 0;
    u64 a = 0;
    u64 pR = 0;

    split_highest_digit(
        p,
        q,
        Q,
        a,
        pR
    );

    const u64 lower =
        inner_sum(
            p,
            Q
        );

    const u64 direct =
        inner_sum(
            p,
            q
        );

    const u64 recurrence =
        recurrence_value(
            p,
            q
        );

    std::cout
        << "p=" << p
        << " q=" << q
        << "\n";

    std::cout
        << "  Q=" << Q
        << " a=" << a
        << " pR=" << pR
        << "\n";

    std::cout
        << "  S(Q)=" << lower
        << "\n";

    std::cout
        << "  (a+1)S(Q)+a(pR-Q)="
        << recurrence
        << "\n";

    std::cout
        << "  S(q)=" << direct
        << "\n";

    std::cout
        << "  q=" << q
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 192\n";

    u64 cases = 0;
    u64 recurrence_pass = 0;
    u64 direct_pass = 0;
    u64 recursive_pass = 0;
    u64 count_recurrence_pass = 0;

    run_small_exhaustive(
        5000,
        cases,
        recurrence_pass,
        direct_pass,
        recursive_pass,
        count_recurrence_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "cases="
        << cases
        << "\n";

    std::cout
        << "inner_recurrence_pass="
        << recurrence_pass
        << "/" << cases
        << "\n";

    std::cout
        << "direct_inner_identity_pass="
        << direct_pass
        << "/" << cases
        << "\n";

    std::cout
        << "recursive_identity_pass="
        << recursive_pass
        << "/" << cases
        << "\n";

    std::cout
        << "interval_count_recurrence_pass="
        << count_recurrence_pass
        << "/" << cases
        << "\n";

    u64 large_recurrence_pass = 0;
    u64 large_direct_pass = 0;
    u64 large_count_pass = 0;

    run_large_random(
        200000,
        large_recurrence_pass,
        large_direct_pass,
        large_count_pass
    );

    std::cout
        << "\nLARGE RANDOM\n";

    std::cout
        << "inner_recurrence_pass="
        << large_recurrence_pass
        << "/200000\n";

    std::cout
        << "direct_inner_identity_pass="
        << large_direct_pass
        << "/200000\n";

    std::cout
        << "interval_count_recurrence_pass="
        << large_count_pass
        << "/200000\n";

    std::cout
        << "\nEXAMPLES\n";

    print_example(
        2,
        5
    );

    print_example(
        2,
        17
    );

    print_example(
        3,
        14
    );

    print_example(
        5,
        194
    );

    std::cout
        << "\nRECURRENCE UNDER TEST\n";

    std::cout
        << "q = Q + a*p^R\n";

    std::cout
        << "S(q) = (a+1)S(Q) + a(p^R-Q)\n";

    std::cout
        << "S(0)=0 => S(q)=q\n";

    std::cout
        << "FINISHED EXPERIMENT 192\n";

    return 0;
}
