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

bool digitwise_leq(u64 p, u64 j, u64 q) {
    while (j > 0 || q > 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd > qd) {
            return false;
        }

        j /= p;
        q /= p;
    }

    return true;
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

u64 lowest_r(
    u64 p,
    u64 j,
    u64 q
) {
    u64 r = 0;

    while (j > 0 || q > 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd < qd) {
            return r;
        }

        j /= p;
        q /= p;
        ++r;
    }

    return UINT64_MAX;
}

/*
 * Count valid start indices by their first strict
 * digit position r.

 * A valid start index satisfies:
 *
 *   j < q
 *   j_i <= q_i
 *
 * and r is the first position with j_r < q_r.
 */
std::vector<u64> observed_r_counts(
    u64 p,
    u64 q
) {
    std::vector<u64> counts;

    u64 temp = q;

    while (temp > 0) {
        counts.push_back(0);
        temp /= p;
    }

    if (counts.empty()) {
        counts.push_back(0);
    }

    for (u64 j = 0; j < q; ++j) {
        if (!digitwise_leq(
                p,
                j,
                q
            )) {
            continue;
        }

        const u64 r =
            lowest_r(
                p,
                j,
                q
            );

        if (r == UINT64_MAX) {
            continue;
        }

        if (r >= counts.size()) {
            counts.resize(r + 1, 0);
        }

        ++counts[r];
    }

    return counts;
}

/*
 * Closed formula:

 *   C_r =
 *       q_r * product_{i>r}(q_i+1)
 */
std::vector<u64> predicted_r_counts(
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

    u64 suffix_product = 1;

    /*
     * Work from MSB toward LSB.

     * At position r:
     *
     *   q_r choices for j_r
     *   product(q_i+1) above r
     */
    for (std::size_t i = digits.size();
         i-- > 0;) {

        const u64 qr = digits[i];

        result[i] =
            qr * suffix_product;

        suffix_product *= qr + 1;
    }

    return result;
}

void run_small_exhaustive(
    u64 max_m,
    u64 &cases,
    u64 &r_cases,
    u64 &count_pass,
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

            const auto observed =
                observed_r_counts(
                    p,
                    par.q
                );

            const auto predicted =
                predicted_r_counts(
                    p,
                    par.q
                );

            ++r_cases;

            if (observed.size() ==
                    predicted.size()) {

                bool ok = true;

                for (std::size_t r = 0;
                     r < predicted.size();
                     ++r) {

                    if (observed[r] !=
                        predicted[r]) {

                        ok = false;
                        break;
                    }
                }

                if (ok) {
                    ++count_pass;
                }
            }

            u64 observed_total = 0;
            u64 predicted_total = 0;

            for (u64 value : observed) {
                observed_total += value;
            }

            for (u64 value : predicted) {
                predicted_total += value;
            }

            /*
             * Number of valid starts should be
             *
             *   product(q_i+1)-1.
             */
            u64 expected_total = 1;

            u64 q = par.q;

            while (q > 0) {
                expected_total *=
                    (q % p) + 1;

                q /= p;
            }

            if (expected_total > 0) {
                --expected_total;
            }

            if (observed_total ==
                    expected_total &&
                predicted_total ==
                    expected_total) {
                ++total_pass;
            }
        }
    }
}

void print_example(
    u64 p,
    u64 q
) {
    const auto observed =
        observed_r_counts(
            p,
            q
        );

    const auto predicted =
        predicted_r_counts(
            p,
            q
        );

    std::cout
        << "p=" << p
        << " q=" << q
        << "\n";

    for (std::size_t r = 0;
         r < predicted.size();
         ++r) {

        std::cout
            << "  r=" << r
            << " observed="
            << (
                r < observed.size()
                ? observed[r]
                : 0
            )
            << " predicted="
            << predicted[r]
            << "\n";
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 186\n";

    u64 cases = 0;
    u64 r_cases = 0;
    u64 count_pass = 0;
    u64 total_pass = 0;

    run_small_exhaustive(
        5000,
        cases,
        r_cases,
        count_pass,
        total_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << cases
        << "\n";

    std::cout
        << "r_distribution_cases="
        << r_cases
        << "\n";

    std::cout
        << "r_count_formula_pass="
        << count_pass
        << "/" << r_cases
        << "\n";

    std::cout
        << "total_start_count_pass="
        << total_pass
        << "/" << cases
        << "\n";

    std::cout
        << "\nEXAMPLES\n";

    print_example(2, 5);
    print_example(2, 17);
    print_example(3, 14);
    print_example(5, 7);

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "C_r = q_r * product_{i>r}(q_i+1)\n";

    std::cout
        << "sum_r C_r = product_i(q_i+1) - 1\n";

    std::cout
        << "FINISHED EXPERIMENT 186\n";

    return 0;
}
