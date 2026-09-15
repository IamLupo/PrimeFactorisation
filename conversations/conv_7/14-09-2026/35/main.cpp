#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Case {
    u64 p;
    u64 q;
};

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

/* =========================================================
 * Integer square root
 * ========================================================= */

static u64 isqrt_u64(u64 n) {
    u64 lo = 0;
    u64 hi = std::min<u64>(
        n,
        1ULL << 32
    );

    while (lo <= hi) {
        const u64 mid =
            lo + (hi - lo) / 2;

        if (mid == 0) {
            lo = 1;
            continue;
        }

        if (mid <= n / mid) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return hi;
}

/* =========================================================
 * Correct prefix MISS oracle
 * ========================================================= */

static u64 prefix_miss(
    u64 m,
    u64 y,
    int base
) {
    if (y >= m) {
        u64 result = 1;
        u64 x = m;

        if (x == 0) {
            return 1;
        }

        while (x > 0) {
            result *=
                (x % base) + 1;

            x /= base;
        }

        return result;
    }

    std::vector<u64> mdigits;
    std::vector<u64> ydigits;

    u64 mx = m;
    u64 yx = y;

    while (
        mx > 0 ||
        yx > 0
    ) {
        mdigits.push_back(
            mx % base
        );

        ydigits.push_back(
            yx % base
        );

        mx /= base;
        yx /= base;
    }

    const std::size_t L =
        std::max(
            mdigits.size(),
            ydigits.size()
        );

    while (mdigits.size() < L) {
        mdigits.push_back(0);
    }

    while (ydigits.size() < L) {
        ydigits.push_back(0);
    }

    std::vector<u64> W(
        L + 1,
        1
    );

    for (
        std::size_t i = 0;
        i < L;
        ++i
    ) {
        W[i + 1] =
            W[i] *
            (mdigits[i] + 1);
    }

    int h =
        static_cast<int>(L) - 1;

    while (
        h >= 0 &&
        ydigits[h] ==
            mdigits[h]
    ) {
        --h;
    }

    if (h < 0) {
        return W[L];
    }

    u64 result = 0;

    for (
        int i =
            static_cast<int>(L) - 1;
        i > h;
        --i
    ) {
        result +=
            ydigits[i] *
            W[i];
    }

    result +=
        ydigits[h] *
        W[h];

    if (h == 0) {
        return result;
    }

    u64 low_m = 0;
    u64 low_y = 0;
    u64 power = 1;

    for (int i = 0; i < h; ++i) {
        low_m +=
            mdigits[i] *
            power;

        low_y +=
            ydigits[i] *
            power;

        power *=
            static_cast<u64>(base);
    }

    result +=
        prefix_miss(
            low_m,
            low_y,
            base
        );

    return result;
}

/* =========================================================
 * Exact coefficient formula
 * ========================================================= */

static std::vector<i128> coefficient_formula(
    u64 m,
    u64 y,
    int base
) {
    std::vector<u64> digits;

    u64 x = m;

    while (x > 0) {
        digits.push_back(
            x % base
        );

        x /= base;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    const int L =
        static_cast<int>(
            digits.size()
        );

    std::vector<u64> W(
        L + 1,
        1
    );

    for (int i = 0; i < L; ++i) {
        W[i + 1] =
            W[i] *
            (digits[i] + 1);
    }

    std::vector<i128> coeff(
        L,
        0
    );

    u64 power_r = 1;

    for (int r = 0; r < L; ++r) {
        const u64 Q =
            y / power_r;

        const u64 R =
            y % power_r;

        if (Q == 0) {
            coeff[r] = 0;
        } else if (
            Q > digits[r]
        ) {
            coeff[r] =
                static_cast<i128>(
                    digits[r]
                ) *
                static_cast<i128>(
                    W[r]
                );
        } else {
            const u64 low_m =
                m % power_r;

            const u64 A =
                prefix_miss(
                    low_m,
                    R,
                    base
                );

            coeff[r] =
                static_cast<i128>(
                    Q - 1
                ) *
                static_cast<i128>(
                    W[r]
                )
                +
                static_cast<i128>(
                    A
                );
        }

        if (r + 1 < L) {
            power_r *=
                static_cast<u64>(base);
        }
    }

    return coeff;
}

/* =========================================================
 * b-adic valuation
 * ========================================================= */

static int valuation_base(
    u64 m,
    int base
) {
    if (m == 0) {
        return -1;
    }

    int v = 0;

    while (
        m % static_cast<u64>(base)
        == 0
    ) {
        m /=
            static_cast<u64>(base);

        ++v;
    }

    return v;
}

/* =========================================================
 * Actual multiplicity of T=0
 *
 * P(T) = c_0 + c_1 T + ...
 *
 * multiplicity is first nonzero coefficient index.
 * ========================================================= */

static int zero_multiplicity(
    const std::vector<i128>& coeff
) {
    for (
        std::size_t i = 0;
        i < coeff.size();
        ++i
    ) {
        if (coeff[i] != 0) {
            return static_cast<int>(i);
        }
    }

    /*
     * Zero polynomial is represented as -1.
     */
    return -1;
}

/* =========================================================
 * Power divisibility helper
 * ========================================================= */

static bool divisible_by_power(
    u64 m,
    int base,
    int k
) {
    if (k <= 0) {
        return true;
    }

    for (int i = 0; i < k; ++i) {
        if (
            m % static_cast<u64>(base)
            != 0
        ) {
            return false;
        }

        m /=
            static_cast<u64>(base);
    }

    return true;
}

/* =========================================================
 * Test one construction
 * ========================================================= */

static bool test_case(
    u64 m,
    u64 y,
    u64 D,
    int base,
    bool print_failure
) {
    const auto P =
        coefficient_formula(
            m,
            y,
            base
        );

    const int expected =
        valuation_base(
            m,
            base
        );

    const int actual =
        zero_multiplicity(P);

    if (actual != expected) {
        if (print_failure) {
            std::cout
                << "FAIL_RAW"
                << " m=" << m
                << " y=" << y
                << " D=" << D
                << " base=" << base
                << " expected="
                << expected
                << " actual="
                << actual
                << "\n";
        }

        return false;
    }

    /*
     * The affine polynomial D+P(T).
     *
     * Its constant coefficient is D+c_0.
     *
     * For D>0 in our semiprime construction,
     * this must be nonzero.
     */
    if (D > 0) {
        auto A = P;

        if (A.empty()) {
            A.push_back(
                static_cast<i128>(D)
            );
        } else {
            A[0] +=
                static_cast<i128>(D);
        }

        const int affine_mult =
            zero_multiplicity(A);

        if (affine_mult != 0) {
            if (print_failure) {
                std::cout
                    << "FAIL_AFFINE"
                    << " m=" << m
                    << " y=" << y
                    << " D=" << D
                    << " base=" << base
                    << " affine_mult="
                    << affine_mult
                    << "\n";
            }

            return false;
        }
    }

    /*
     * Explicitly verify:
     *
     * c_0 = ... = c_(v-1) = 0
     * c_v != 0
     */
    if (expected > 0) {
        for (int r = 0; r < expected; ++r) {
            if (P[r] != 0) {
                if (print_failure) {
                    std::cout
                        << "FAIL_LOWER"
                        << " m=" << m
                        << " base=" << base
                        << " r=" << r
                        << "\n";
                }

                return false;
            }
        }

        if (
            expected >=
            static_cast<int>(
                P.size()
            )
            ||
            P[expected] == 0
        ) {
            if (print_failure) {
                std::cout
                    << "FAIL_FIRST_NONZERO"
                    << " m=" << m
                    << " y=" << y
                    << " base=" << base
                    << " v=" << expected
                    << "\n";
            }

            return false;
        }
    }

    return true;
}

int main() {
    std::cout
        << "START EXPERIMENT 325\n"
        << "T-ADIC MULTIPLICITY THEOREM\n"
        << "IS mult_0(P_b)=v_b(m) EXACTLY?\n\n";

    const std::array<Case, 30> cases = {{
        {81077, 162749},
        {125017, 174259},
        {57107, 88261},
        {52757, 148457},
        {19483, 34123},
        {123001, 188291},
        {97987, 112583},
        {76129, 192113},
        {41257, 73643},
        {137239, 140419},
        {64187, 115319},
        {162527, 184087},
        {66179, 124753},
        {87943, 158047},
        {42683, 52529},
        {107473, 149711},
        {75853, 82759},
        {150131, 175267},
        {99859, 124769},
        {75337, 162229},
        {124471, 168043},
        {51563, 87683},
        {52237, 142123},
        {19069, 28751},
        {117043, 187637},
        {92203, 112067},
        {75617, 185869},
        {40823, 67843},
        {134369, 136709},
        {63667, 109211}
    }};

    u64 semiprime_tests = 0;
    u64 semiprime_failures = 0;

    u64 random_tests = 0;
    u64 random_failures = 0;

    u64 edge_tests = 0;
    u64 edge_failures = 0;

    std::array<u64, 6> valuation_hist{};
    u64 max_v = 0;

    /*
     * --------------------------------------------------
     * 1. Existing semiprime constructions
     * --------------------------------------------------
     */

    for (
        std::size_t case_idx = 0;
        case_idx < cases.size();
        ++case_idx
    ) {
        const u64 p =
            cases[case_idx].p;

        const u64 q =
            cases[case_idx].q;

        const u64 N =
            p * q;

        const u64 s =
            isqrt_u64(N);

        for (
            int offset = 1;
            offset <= 5;
            ++offset
        ) {
            const u64 off =
                static_cast<u64>(
                    offset
                );

            if (s < off) {
                continue;
            }

            const u64 m =
                s - off + 1;

            const u64 D =
                N - m * m;

            for (int base : BASES) {
                ++semiprime_tests;

                const int v =
                    valuation_base(
                        m,
                        base
                    );

                if (v >= 0 &&
                    v <
                    static_cast<int>(
                        valuation_hist.size()
                    )) {

                    ++valuation_hist[v];
                }

                max_v =
                    std::max<u64>(
                        max_v,
                        v < 0
                            ? 0
                            : static_cast<u64>(v)
                    );

                if (
                    !test_case(
                        m,
                        N,
                        D,
                        base,
                        semiprime_failures < 10
                    )
                ) {
                    ++semiprime_failures;
                }
            }
        }
    }

    /*
     * --------------------------------------------------
     * 2. Exhaustive small m/y tests
     *
     * This is independent of the semiprime structure.
     * --------------------------------------------------
     */

    for (
        int base : BASES
    ) {
        for (
            u64 m = 1;
            m <= 500;
            ++m
        ) {
            /*
             * Test several y regimes:
             * y=m, y>m, y<m.
             */
            const std::array<u64, 4> ys = {
                m,
                m + 1,
                m > 1 ? m - 1 : 0,
                2 * m + 3
            };

            for (u64 y : ys) {
                ++random_tests;

                if (
                    !test_case(
                        m,
                        y,
                        1,
                        base,
                        random_failures < 10
                    )
                ) {
                    ++random_failures;
                }
            }
        }
    }

    /*
     * --------------------------------------------------
     * 3. Random larger tests
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x325325325ULL
    );

    for (
        u64 i = 0;
        i < 100000;
        ++i
    ) {
        const int base =
            BASES[
                rng() % BASES.size()
            ];

        const u64 m =
            1 +
            (
                rng()
                % 1000000000000ULL
            );

        /*
         * Deliberately use several y relations to m.
         */
        const u64 mode =
            rng() % 4;

        u64 y;

        if (mode == 0) {
            y = m;
        } else if (mode == 1) {
            y = m - 1;
        } else if (mode == 2) {
            y = m + 1;
        } else {
            y =
                rng()
                % (
                    m +
                    1000000ULL
                );
        }

        ++edge_tests;

        if (
            !test_case(
                m,
                y,
                1,
                base,
                edge_failures < 10
            )
        ) {
            ++edge_failures;
        }
    }

    /*
     * --------------------------------------------------
     * Results
     * --------------------------------------------------
     */

    std::cout
        << "============================\n"
        << "SEMIPRIME TESTS\n"
        << "tests="
        << semiprime_tests
        << "\n"
        << "failures="
        << semiprime_failures
        << "\n\n";

    std::cout
        << "SMALL EXHAUSTIVE TESTS\n"
        << "tests="
        << random_tests
        << "\n"
        << "failures="
        << random_failures
        << "\n\n";

    std::cout
        << "RANDOM TESTS\n"
        << "tests="
        << edge_tests
        << "\n"
        << "failures="
        << edge_failures
        << "\n\n";

    std::cout
        << "VALUATION HISTOGRAM\n";

    for (
        std::size_t i = 0;
        i < valuation_hist.size();
        ++i
    ) {
        std::cout
            << "v="
            << i
            << " count="
            << valuation_hist[i]
            << "\n";
    }

    std::cout
        << "max_v="
        << max_v
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 325\n";

    return 0;
}
