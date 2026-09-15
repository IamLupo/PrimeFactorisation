#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Parameters {
    u64 p;
    u64 m;
    u64 e;
    u64 s0;
    u64 q;
    u64 pe;
};

u64 pow_u64(u64 p, u64 e) {
    u128 result = 1;

    for (u64 i = 0; i < e; ++i) {
        result *= p;
    }

    return static_cast<u64>(result);
}

bool is_hit(u64 p, u64 m, u64 x) {
    u64 a = x;
    u64 b = m;

    while (a > 0 || b > 0) {
        const u64 ad = a % p;
        const u64 bd = b % p;

        if (ad > bd) {
            return true;
        }

        a /= p;
        b /= p;
    }

    return false;
}

bool derive_parameters(
    u64 p,
    u64 m,
    Parameters &out
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

    const u64 e = a + 1 + z;

    const u64 pa = pow_u64(p, a);
    const u64 pe = pow_u64(p, e);

    const u128 numerator =
        static_cast<u128>(m) +
        1 -
        static_cast<u128>((b + 1) * pa);

    if (numerator == 0 ||
        numerator % pe != 0) {
        return false;
    }

    const u64 q =
        static_cast<u64>(numerator / pe);

    if (q == 0) {
        return false;
    }

    out = {
        p,
        m,
        e,
        (b + 1) * pa,
        q,
        pe
    };

    return true;
}

/*
 * Direct count of HIT integers by scanning t.
 *
 * This is only used for small exhaustive validation.
 */
u64 actual_hit_count(
    u64 p,
    u64 m
) {
    u64 count = 0;

    for (u64 t = 1; t <= m; ++t) {
        if (is_hit(p, m, t)) {
            ++count;
        }
    }

    return count;
}

/*
 * Predicted number of intervals of type r:

 *   C_r = q_r * product_{i>r}(q_i+1)
 */
std::vector<u64> interval_counts_by_r(
    u64 p,
    u64 q
) {
    std::vector<u64> digits;

    u64 x = q;

    while (x > 0) {
        digits.push_back(x % p);
        x /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    std::vector<u64> result(
        digits.size(),
        0
    );

    u64 product_higher = 1;

    for (std::size_t i = digits.size();
         i-- > 0;) {

        result[i] =
            digits[i] *
            product_higher;

        product_higher *=
            digits[i] + 1;
    }

    return result;
}

/*
 * Length of every interval of type r:

 *   L_r =
 *     (p^r - (q mod p^r))*p^e - s0
 */
u64 interval_length(
    u64 p,
    u64 e,
    u64 s0,
    u64 q,
    u64 r
) {
    const u64 pr =
        pow_u64(p, r);

    const u64 qlow =
        q % pr;

    const u64 pe =
        pow_u64(p, e);

    const u128 result =
        (static_cast<u128>(pr) -
         static_cast<u128>(qlow)) *
        static_cast<u128>(pe) -
        static_cast<u128>(s0);

    return static_cast<u64>(result);
}

/*
 * Complete predicted HIT count:

 *   H =
 *     sum_r C_r * L_r
 */
u64 predicted_hit_count(
    const Parameters &par
) {
    const auto counts =
        interval_counts_by_r(
            par.p,
            par.q
        );

    u128 total = 0;

    for (std::size_t r = 0;
         r < counts.size();
         ++r) {

        if (counts[r] == 0) {
            continue;
        }

        const u64 length =
            interval_length(
                par.p,
                par.e,
                par.s0,
                par.q,
                static_cast<u64>(r)
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
 * Candidate simplification.

 * We test whether

 *   H = m - (number of MISS values)
 *
 * can itself be written in a simple digit product form.

 * First candidate:

 *   MISS = product_i(q_i+1)

 * This is NOT assumed to be true; it is tested.
 */
u64 candidate_miss_count_1(
    u64 p,
    u64 q
) {
    u64 product = 1;

    while (q > 0) {
        product *=
            (q % p) + 1;

        q /= p;
    }

    return product;
}

void run_small_exhaustive(
    u64 max_m,
    u64 &cases,
    u64 &predicted_pass,
    u64 &candidate_pass,
    u64 &total_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (const u64 p : primes) {
        for (u64 m = 2;
             m <= max_m;
             ++m) {

            Parameters par{};

            if (!derive_parameters(
                    p,
                    m,
                    par
                )) {
                continue;
            }

            ++cases;

            const u64 actual =
                actual_hit_count(
                    p,
                    m
                );

            const u64 predicted =
                predicted_hit_count(
                    par
                );

            if (actual == predicted) {
                ++predicted_pass;
            }

            /*
             * Total number of valid start indices:
             *
             * product(q_i+1)-1.
             *
             * This is the number of intervals, NOT the
             * number of HIT integers.
             */
            const u64 candidate_miss =
                candidate_miss_count_1(
                    p,
                    par.q
                );

            if (actual ==
                par.m -
                candidate_miss) {
                ++candidate_pass;
            }

            /*
             * Total-length consistency.
             */
            u128 sum = 0;

            const auto counts =
                interval_counts_by_r(
                    p,
                    par.q
                );

            for (std::size_t r = 0;
                 r < counts.size();
                 ++r) {

                const u64 length =
                    interval_length(
                        p,
                        par.e,
                        par.s0,
                        par.q,
                        static_cast<u64>(r)
                    );

                sum +=
                    static_cast<u128>(
                        counts[r]
                    ) *
                    static_cast<u128>(
                        length
                    );
            }

            if (sum ==
                static_cast<u128>(actual)) {
                ++total_pass;
            }
        }
    }
}

void print_example(
    u64 p,
    u64 m
) {
    Parameters par{};

    if (!derive_parameters(
            p,
            m,
            par
        )) {
        return;
    }

    const auto counts =
        interval_counts_by_r(
            p,
            par.q
        );

    std::cout
        << "p=" << p
        << " m=" << m
        << " e=" << par.e
        << " s0=" << par.s0
        << " q=" << par.q
        << "\n";

    u64 total = 0;

    for (std::size_t r = 0;
         r < counts.size();
         ++r) {

        if (counts[r] == 0) {
            continue;
        }

        const u64 length =
            interval_length(
                p,
                par.e,
                par.s0,
                par.q,
                static_cast<u64>(r)
            );

        const u64 contribution =
            counts[r] * length;

        total += contribution;

        std::cout
            << "  r=" << r
            << " count="
            << counts[r]
            << " length="
            << length
            << " contribution="
            << contribution
            << "\n";
    }

    std::cout
        << "  predicted_total="
        << total
        << "\n";

    std::cout
        << "  actual_total="
        << actual_hit_count(
            p,
            m
        )
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 187\n";

    u64 cases = 0;
    u64 predicted_pass = 0;
    u64 candidate_pass = 0;
    u64 total_pass = 0;

    run_small_exhaustive(
        5000,
        cases,
        predicted_pass,
        candidate_pass,
        total_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << cases
        << "\n";

    std::cout
        << "weighted_interval_sum_pass="
        << predicted_pass
        << "/" << cases
        << "\n";

    std::cout
        << "candidate_miss_formula_pass="
        << candidate_pass
        << "/" << cases
        << "\n";

    std::cout
        << "total_length_pass="
        << total_pass
        << "/" << cases
        << "\n";

    std::cout
        << "\nEXAMPLES\n";

    print_example(2, 20);
    print_example(2, 40);
    print_example(3, 41);
    print_example(5, 194);

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "H = sum_r C_r * L_r\n";

    std::cout
        << "C_r = q_r * product_{i>r}(q_i+1)\n";

    std::cout
        << "L_r = (p^r - (q mod p^r))*p^e - s0\n";

    std::cout
        << "FINISHED EXPERIMENT 187\n";

    return 0;
}
