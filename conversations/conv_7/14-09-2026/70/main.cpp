#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

struct CoeffSet {
    i64 A;
    i64 B;
    i64 C;
    i64 D;
};

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        __builtin_sqrtl(
            static_cast<long double>(n)
        )
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(
        static_cast<std::size_t>(limit + 1),
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<long long>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 50000) {
            usable.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        usable.size() - 1
    );

    std::vector<CaseData> result;
    result.reserve(
        static_cast<std::size_t>(count)
    );

    while (
        static_cast<int>(result.size()) < count
    ) {
        const std::size_t i = dist(rng);
        const std::size_t j = dist(rng);

        if (i == j) {
            continue;
        }

        const u64 p =
            static_cast<u64>(
                std::min(
                    usable[i],
                    usable[j]
                )
            );

        const u64 q =
            static_cast<u64>(
                std::max(
                    usable[i],
                    usable[j]
                )
            );

        const u64 N = p * q;

        if (!used_N.insert(N).second) {
            continue;
        }

        result.push_back({
            p,
            q,
            N,
            isqrt_u64(N)
        });
    }

    return result;
}

static i128 evaluate_family(
    const CoeffSet& coeff,
    u64 x,
    const CaseData& c
) {
    const i128 xx =
        static_cast<i128>(x);

    const i128 s =
        static_cast<i128>(c.s);

    const i128 N =
        static_cast<i128>(c.N);

    const i128 R =
        static_cast<i128>(coeff.A) * xx
        +
        static_cast<i128>(coeff.B) * s
        +
        static_cast<i128>(coeff.D);

    return
        static_cast<i128>(coeff.C) * N
        +
        (xx - (s + 1)) * R;
}

static i128 quotient_formula(
    const CoeffSet& coeff,
    const CaseData& c
) {
    return
        static_cast<i128>(coeff.A) *
            static_cast<i128>(c.p)
        +
        static_cast<i128>(coeff.C) *
            static_cast<i128>(c.q)
        -
        static_cast<i128>(
            coeff.A + coeff.B
        ) *
            static_cast<i128>(c.s)
        -
        static_cast<i128>(
            coeff.A + coeff.D
        );
}

static u64 abs_i128_mod(
    i128 value,
    u64 mod
) {
    if (value < 0) {
        value = -value;
    }

    return static_cast<u64>(
        static_cast<u128>(value) %
        static_cast<u128>(mod)
    );
}

static u64 gcd_abs(
    i128 value,
    u64 N
) {
    return std::gcd(
        abs_i128_mod(value, N),
        N
    );
}

static void print_i128(i128 value) {
    if (value == 0) {
        std::cout << "0";
        return;
    }

    if (value < 0) {
        std::cout << '-';
        value = -value;
    }

    std::string out;

    while (value > 0) {
        const unsigned digit =
            static_cast<unsigned>(value % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        value /= 10;
    }

    std::reverse(
        out.begin(),
        out.end()
    );

    std::cout << out;
}

int main() {
    std::cout << "START EXPERIMENT 359\n";

    constexpr int CASE_COUNT = 1000;
    constexpr int PRIME_LIMIT = 50000;

    std::mt19937_64 rng(
        0x359359359ULL
    );

    const std::vector<int> primes =
        sieve_primes(
            PRIME_LIMIT
        );

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    /*
     * ------------------------------------------------------------
     * Three distinguished polynomial families.
     *
     * H:
     *   (A,B,C,D)=(4,3,3,0)
     *
     * Difference polynomial:
     *   quotient = q-p
     *
     * Center polynomial:
     *   quotient = p+q-2s
     * ------------------------------------------------------------
     */

    const CoeffSet H{
        4, 3, 3, 0
    };

    const CoeffSet P_DIFF{
        -1, 1, 1, 1
    };

    const CoeffSet P_CENTER{
        1, 1, 1, -1
    };

    u64 H_failures = 0;
    u64 diff_failures = 0;
    u64 center_failures = 0;

    u64 H_gcd_failures = 0;
    u64 diff_gcd_failures = 0;
    u64 center_gcd_failures = 0;

    u64 reconstructed_difference_failures = 0;
    u64 reconstructed_center_failures = 0;
    u64 reconstructed_H_failures = 0;

    u64 linear_combination_failures = 0;

    u64 H_gcd_exact_p = 0;
    u64 diff_gcd_exact_p = 0;
    u64 center_gcd_exact_p = 0;

    u64 H_gcd_exact_q = 0;
    u64 diff_gcd_exact_q = 0;
    u64 center_gcd_exact_q = 0;

    u64 diff_zero = 0;
    u64 center_zero = 0;
    u64 H_zero = 0;

    u64 diff_gcd_one = 0;
    u64 center_gcd_one = 0;
    u64 H_gcd_one = 0;

    u64 diff_abs_lt_s = 0;
    u64 center_abs_lt_s = 0;
    u64 H_abs_lt_s = 0;

    /*
     * Track gcd relationships between the three quotients.
     */
    u64 gcd_quotients_all_one = 0;
    u64 gcd_diff_center_one = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c =
            cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 xp =
            s + 1 - p;

        /*
         * --------------------------------------------------------
         * Evaluate all three polynomials at the hidden coordinate.
         * --------------------------------------------------------
         */

        const i128 H_value =
            evaluate_family(
                H,
                xp,
                c
            );

        const i128 diff_value =
            evaluate_family(
                P_DIFF,
                xp,
                c
            );

        const i128 center_value =
            evaluate_family(
                P_CENTER,
                xp,
                c
            );

        /*
         * Exact quotients.
         */
        const i128 H_q =
            quotient_formula(
                H,
                c
            );

        const i128 diff_q =
            quotient_formula(
                P_DIFF,
                c
            );

        const i128 center_q =
            quotient_formula(
                P_CENTER,
                c
            );

        /*
         * --------------------------------------------------------
         * Expected expressions.
         * --------------------------------------------------------
         */

        const i128 expected_diff =
            static_cast<i128>(q)
            -
            static_cast<i128>(p);

        const i128 expected_center =
            static_cast<i128>(p)
            +
            static_cast<i128>(q)
            -
            static_cast<i128>(2) *
            static_cast<i128>(s);

        const i128 expected_H =
            static_cast<i128>(4) *
                static_cast<i128>(p)
            +
            static_cast<i128>(3) *
                static_cast<i128>(q)
            -
            static_cast<i128>(7) *
                static_cast<i128>(s)
            -
            static_cast<i128>(4);

        if (diff_q != expected_diff) {
            ++reconstructed_difference_failures;
        }

        if (center_q != expected_center) {
            ++reconstructed_center_failures;
        }

        if (H_q != expected_H) {
            ++reconstructed_H_failures;
        }

        /*
         * Check evaluation = p * quotient.
         */
        if (
            H_value !=
            static_cast<i128>(p) * H_q
        ) {
            ++H_failures;
        }

        if (
            diff_value !=
            static_cast<i128>(p) * diff_q
        ) {
            ++diff_failures;
        }

        if (
            center_value !=
            static_cast<i128>(p) * center_q
        ) {
            ++center_failures;
        }

        /*
         * All three should therefore be divisible by p.
         */
        if (
            abs_i128_mod(
                H_value,
                p
            ) != 0
        ) {
            ++H_gcd_failures;
        }

        if (
            abs_i128_mod(
                diff_value,
                p
            ) != 0
        ) {
            ++diff_gcd_failures;
        }

        if (
            abs_i128_mod(
                center_value,
                p
            ) != 0
        ) {
            ++center_gcd_failures;
        }

        /*
         * --------------------------------------------------------
         * GCD with N.
         * --------------------------------------------------------
         */

        const u64 H_g =
            gcd_abs(
                H_value,
                N
            );

        const u64 diff_g =
            gcd_abs(
                diff_value,
                N
            );

        const u64 center_g =
            gcd_abs(
                center_value,
                N
            );

        if (H_g == p) {
            ++H_gcd_exact_p;
        } else if (H_g == q) {
            ++H_gcd_exact_q;
        } else if (H_g == 1) {
            ++H_gcd_one;
        }

        if (diff_g == p) {
            ++diff_gcd_exact_p;
        } else if (diff_g == q) {
            ++diff_gcd_exact_q;
        } else if (diff_g == 1) {
            ++diff_gcd_one;
        }

        if (center_g == p) {
            ++center_gcd_exact_p;
        } else if (center_g == q) {
            ++center_gcd_exact_q;
        } else if (center_g == 1) {
            ++center_gcd_one;
        }

        /*
         * p should be the gcd unless the quotient happens
         * to contain q as well.
         */
        if (
            H_g != p
        ) {
            ++H_gcd_failures;
        }

        if (
            diff_g != p
        ) {
            ++diff_gcd_failures;
        }

        if (
            center_g != p
        ) {
            ++center_gcd_failures;
        }

        /*
         * --------------------------------------------------------
         * Linear relationship:
         *
         * 2H_q = 7 center_q - diff_q - 8.
         * --------------------------------------------------------
         */

        const i128 reconstructed_H =
            static_cast<i128>(7) *
            center_q
            -
            diff_q
            -
            static_cast<i128>(8);

        if (
            static_cast<i128>(2) * H_q
            != reconstructed_H
        ) {
            ++linear_combination_failures;
        }

        /*
         * --------------------------------------------------------
         * Direct hidden-coordinate quantities:
         *
         * q-p
         * p+q-2s
         * H quotient.
         * --------------------------------------------------------
         */

        const i128 abs_diff =
            diff_q < 0
                ? -diff_q
                : diff_q;

        const i128 abs_center =
            center_q < 0
                ? -center_q
                : center_q;

        const i128 abs_H =
            H_q < 0
                ? -H_q
                : H_q;

        if (diff_q == 0) {
            ++diff_zero;
        }

        if (center_q == 0) {
            ++center_zero;
        }

        if (H_q == 0) {
            ++H_zero;
        }

        if (
            abs_diff < static_cast<i128>(s)
        ) {
            ++diff_abs_lt_s;
        }

        if (
            abs_center < static_cast<i128>(s)
        ) {
            ++center_abs_lt_s;
        }

        if (
            abs_H < static_cast<i128>(s)
        ) {
            ++H_abs_lt_s;
        }

        /*
         * The quotient of the two constructed objects:
         *
         * gcd(q-p, p+q-2s)
         */
        const u64 quotient_gcd =
            std::gcd(
                static_cast<u64>(abs_diff),
                static_cast<u64>(abs_center)
            );

        if (quotient_gcd == 1) {
            ++gcd_diff_center_one;
        }

        const u64 all_gcd =
            std::gcd(
                quotient_gcd,
                static_cast<u64>(abs_H)
            );

        if (all_gcd == 1) {
            ++gcd_quotients_all_one;
        }

        /*
         * --------------------------------------------------------
         * Print a few representative cases.
         * --------------------------------------------------------
         */

        if (case_id < 10) {
            std::cout
                << "\nCASE="
                << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << "\n";

            std::cout
                << "q-p=";

            print_i128(
                expected_diff
            );

            std::cout
                << "\n";

            std::cout
                << "p+q-2s=";

            print_i128(
                expected_center
            );

            std::cout
                << "\n";

            std::cout
                << "H_quotient=";

            print_i128(
                H_q
            );

            std::cout
                << "\n";

            std::cout
                << "2H_quotient=";

            print_i128(
                static_cast<i128>(2) * H_q
            );

            std::cout
                << "\n";

            std::cout
                << "7center-diff-8=";

            print_i128(
                reconstructed_H
            );

            std::cout
                << "\n";

            std::cout
                << "GCD_DIFF="
                << diff_g
                << "\n";

            std::cout
                << "GCD_CENTER="
                << center_g
                << "\n";

            std::cout
                << "GCD_H="
                << H_g
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "H_EVALUATION_FAILURES="
        << H_failures
        << "\n";

    std::cout
        << "DIFF_EVALUATION_FAILURES="
        << diff_failures
        << "\n";

    std::cout
        << "CENTER_EVALUATION_FAILURES="
        << center_failures
        << "\n";

    std::cout
        << "H_DIVISIBILITY_FAILURES="
        << H_gcd_failures
        << "\n";

    std::cout
        << "DIFF_DIVISIBILITY_FAILURES="
        << diff_gcd_failures
        << "\n";

    std::cout
        << "CENTER_DIVISIBILITY_FAILURES="
        << center_gcd_failures
        << "\n";

    std::cout
        << "RECONSTRUCTED_DIFFERENCE_FAILURES="
        << reconstructed_difference_failures
        << "\n";

    std::cout
        << "RECONSTRUCTED_CENTER_FAILURES="
        << reconstructed_center_failures
        << "\n";

    std::cout
        << "RECONSTRUCTED_H_FAILURES="
        << reconstructed_H_failures
        << "\n";

    std::cout
        << "LINEAR_COMBINATION_FAILURES="
        << linear_combination_failures
        << "\n";

    std::cout
        << "H_GCD_EXACT_P="
        << H_gcd_exact_p
        << "\n";

    std::cout
        << "DIFF_GCD_EXACT_P="
        << diff_gcd_exact_p
        << "\n";

    std::cout
        << "CENTER_GCD_EXACT_P="
        << center_gcd_exact_p
        << "\n";

    std::cout
        << "H_GCD_EXACT_Q="
        << H_gcd_exact_q
        << "\n";

    std::cout
        << "DIFF_GCD_EXACT_Q="
        << diff_gcd_exact_q
        << "\n";

    std::cout
        << "CENTER_GCD_EXACT_Q="
        << center_gcd_exact_q
        << "\n";

    std::cout
        << "H_GCD_ONE="
        << H_gcd_one
        << "\n";

    std::cout
        << "DIFF_GCD_ONE="
        << diff_gcd_one
        << "\n";

    std::cout
        << "CENTER_GCD_ONE="
        << center_gcd_one
        << "\n";

    std::cout
        << "DIFF_ZERO="
        << diff_zero
        << "\n";

    std::cout
        << "CENTER_ZERO="
        << center_zero
        << "\n";

    std::cout
        << "H_ZERO="
        << H_zero
        << "\n";

    std::cout
        << "DIFF_ABS_LT_S="
        << diff_abs_lt_s
        << "\n";

    std::cout
        << "CENTER_ABS_LT_S="
        << center_abs_lt_s
        << "\n";

    std::cout
        << "H_ABS_LT_S="
        << H_abs_lt_s
        << "\n";

    std::cout
        << "GCD_DIFF_CENTER_ONE="
        << gcd_diff_center_one
        << "\n";

    std::cout
        << "GCD_ALL_THREE_QUOTIENTS_ONE="
        << gcd_quotients_all_one
        << "\n";

    std::cout
        << "\n"
        << "P_DIFF_COEFFICIENTS="
        << "(-1,1,1,1)"
        << "\n";

    std::cout
        << "P_CENTER_COEFFICIENTS="
        << "(1,1,1,-1)"
        << "\n";

    std::cout
        << "H_COEFFICIENTS="
        << "(4,3,3,0)"
        << "\n";

    if (
        H_failures == 0 &&
        diff_failures == 0 &&
        center_failures == 0 &&
        reconstructed_difference_failures == 0 &&
        reconstructed_center_failures == 0 &&
        reconstructed_H_failures == 0 &&
        linear_combination_failures == 0 &&
        H_gcd_failures == 0 &&
        diff_gcd_failures == 0 &&
        center_gcd_failures == 0
    ) {
        std::cout
            << "STRUCTURE_STATUS="
            << "QUOTIENT_BASIS_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 359\n";

    return 0;
}
