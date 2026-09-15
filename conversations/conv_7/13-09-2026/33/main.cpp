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
 * Product of (m_i + 1) over the base-p digits of m.
 */
u64 digit_product(
    u64 p,
    u64 m
) {
    u64 product = 1;

    while (m > 0) {
        product *=
            (m % p) + 1;

        m /= p;
    }

    return product;
}

/*
 * Direct Lucas count:

 *   MISS values among 0..m:
 *
 *       product_i (m_i + 1)
 *
 * Therefore:
 *
 *   HIT values among 1..m:
 *
 *       m + 1 - product_i(m_i + 1)
 */
u64 predicted_hit_count(
    u64 p,
    u64 m
) {
    const u64 misses_including_zero =
        digit_product(
            p,
            m
        );

    return m + 1 -
           misses_including_zero;
}

/*
 * Keep the interval decomposition from the previous
 * experiments as an independent calculation.
 */
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

std::vector<u64> interval_counts_by_r(
    u64 p,
    u64 q
) {
    std::vector<u64> digits;

    u64 x = q;

    while (x > 0) {
        digits.push_back(
            x % p
        );

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

u64 weighted_interval_sum(
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

void run_small_exhaustive(
    u64 max_m,
    u64 &cases,
    u64 &direct_digit_pass,
    u64 &interval_pass,
    u64 &agreement_pass
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

            const u64 digit_formula =
                predicted_hit_count(
                    p,
                    m
                );

            const u64 interval_formula =
                weighted_interval_sum(
                    par
                );

            if (actual ==
                digit_formula) {
                ++direct_digit_pass;
            }

            if (actual ==
                interval_formula) {
                ++interval_pass;
            }

            if (digit_formula ==
                interval_formula) {
                ++agreement_pass;
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

    const u64 digits_product =
        digit_product(
            p,
            m
        );

    const u64 digit_formula =
        predicted_hit_count(
            p,
            m
        );

    const u64 interval_formula =
        weighted_interval_sum(
            par
        );

    const u64 actual =
        actual_hit_count(
            p,
            m
        );

    std::cout
        << "p=" << p
        << " m=" << m
        << "\n"
        << "  product(m_i+1)="
        << digits_product
        << "\n"
        << "  digit_formula="
        << digit_formula
        << "\n"
        << "  interval_formula="
        << interval_formula
        << "\n"
        << "  actual="
        << actual
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 188\n";

    u64 cases = 0;
    u64 direct_digit_pass = 0;
    u64 interval_pass = 0;
    u64 agreement_pass = 0;

    run_small_exhaustive(
        5000,
        cases,
        direct_digit_pass,
        interval_pass,
        agreement_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << cases
        << "\n";

    std::cout
        << "direct_digit_formula_pass="
        << direct_digit_pass
        << "/" << cases
        << "\n";

    std::cout
        << "interval_formula_pass="
        << interval_pass
        << "/" << cases
        << "\n";

    std::cout
        << "digit_vs_interval_agreement="
        << agreement_pass
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
        << "MISS(0..m) = product_i(m_i+1)\n";

    std::cout
        << "HIT(1..m) = m + 1 - product_i(m_i+1)\n";

    std::cout
        << "This must equal the weighted interval sum.\n";

    std::cout
        << "FINISHED EXPERIMENT 188\n";

    return 0;
}
