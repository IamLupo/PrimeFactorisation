
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 local_formula_failures = 0;
    u64 global_formula_failures = 0;
    u64 exact_failures = 0;
    u64 total_failures = 0;
    u64 zero_layer_failures = 0;

    u64 targeted = 0;
};

/* =========================================================
 * v_p(n!)
 * ========================================================= */

static u64 factorial_valuation(
    u64 n,
    int p
) {
    u64 result = 0;
    const u64 prime =
        static_cast<u64>(p);

    while (n > 0) {
        n /= prime;
        result += n;
    }

    return result;
}

/* =========================================================
 * Exact v_p(C(n,k))
 * ========================================================= */

static u64 binomial_valuation(
    u64 n,
    u64 k,
    int p
) {
    return
        factorial_valuation(n, p)
        - factorial_valuation(k, p)
        - factorial_valuation(n - k, p);
}

/* =========================================================
 * p^r
 * ========================================================= */

static u64 prime_power(
    int p,
    int r
) {
    u64 result = 1;

    for (int i = 0;
         i < r;
         ++i) {
        result *=
            static_cast<u64>(p);
    }

    return result;
}

/* =========================================================
 * Base-p digits
 * ========================================================= */

static std::vector<u64>
digits_of(
    u64 n,
    int p,
    int length
) {
    std::vector<u64> digits(
        static_cast<std::size_t>(length),
        0
    );

    const u64 prime =
        static_cast<u64>(p);

    for (int i = 0;
         i < length;
         ++i) {
        digits[i] =
            n % prime;

        n /= prime;
    }

    return digits;
}

/* =========================================================
 * MISS count
 * ========================================================= */

static u64 miss_count(
    u64 m,
    int p
) {
    const auto digits =
        digits_of(
            m,
            p,
            64
        );

    u64 result = 1;

    bool still_nonzero = true;

    for (u64 d : digits) {
        if (still_nonzero) {
            result *= d + 1;
        }

        if (d == 0) {
            /*
             * This does not mean higher digits vanish, so we
             * cannot stop. The full fixed-length digit vector
             * is not used here as a stopping condition.
             */
        }
    }

    /*
     * Recompute properly using m itself.
     */
    result = 1;

    u64 n = m;
    const u64 prime =
        static_cast<u64>(p);

    while (n > 0) {
        result *=
            n % prime + 1;

        n /= prime;
    }

    /*
     * m=0 has one MISS value.
     */
    if (m == 0) {
        return 1;
    }

    return result;
}

/* =========================================================
 * Exact global valuation histogram
 *
 * Small m only.
 * ========================================================= */

static std::vector<u64>
exact_global_histogram(
    u64 m,
    int p
) {
    const auto digits =
        digits_of(
            m,
            p,
            64
        );

    std::vector<u64> result(
        digits.size() + 1,
        0
    );

    for (u64 t = 0;
         t <= m;
         ++t) {

        const u64 v =
            binomial_valuation(
                m,
                t,
                p
            );

        if (
            result.size() <= v
        ) {
            result.resize(
                static_cast<std::size_t>(v + 1),
                0
            );
        }

        ++result[v];
    }

    return result;
}

/* =========================================================
 * Global Kummer digit DP
 *
 * Independent reference implementation.
 * ========================================================= */

static std::vector<u64>
global_kummer_dp(
    u64 m,
    int p
) {
    std::vector<u64> digits;

    u64 n = m;
    const u64 prime =
        static_cast<u64>(p);

    if (n == 0) {
        digits.push_back(0);
    } else {
        while (n > 0) {
            digits.push_back(
                n % prime
            );

            n /= prime;
        }
    }

    const int length =
        static_cast<int>(
            digits.size()
        );

    std::vector<u64> dp0(
        static_cast<std::size_t>(length + 1),
        0
    );

    std::vector<u64> dp1(
        static_cast<std::size_t>(length + 1),
        0
    );

    dp0[0] = 1;

    for (int i = 0;
         i < length;
         ++i) {

        std::vector<u64> next0(
            static_cast<std::size_t>(length + 1),
            0
        );

        std::vector<u64> next1(
            static_cast<std::size_t>(length + 1),
            0
        );

        const u64 md =
            digits[i];

        for (int k = 0;
             k <= i;
             ++k) {

            if (dp0[k] != 0) {
                for (
                    u64 td = 0;
                    td < static_cast<u64>(p);
                    ++td
                ) {
                    if (md >= td) {
                        next0[k] +=
                            dp0[k];
                    } else {
                        next1[k + 1] +=
                            dp0[k];
                    }
                }
            }

            if (dp1[k] != 0) {
                for (
                    u64 td = 0;
                    td < static_cast<u64>(p);
                    ++td
                ) {
                    if (md >= td + 1) {
                        next0[k] +=
                            dp1[k];
                    } else {
                        next1[k + 1] +=
                            dp1[k];
                    }
                }
            }
        }

        dp0.swap(next0);
        dp1.swap(next1);
    }

    return dp0;
}

/* =========================================================
 * Lower-borrow distribution B_h(j)
 *
 * We process the first h digits of L.
 *
 * state 0:
 *   no borrow entering current digit
 *
 * state 1:
 *   borrow entering current digit
 *
 * The arrays count choices of lower y digits.
 * ========================================================= */

static std::vector<u64>
lower_borrow_distribution(
    const std::vector<u64>& Ldigits,
    int p,
    int h
) {
    std::vector<u64> dp0(
        static_cast<std::size_t>(h + 1),
        0
    );

    std::vector<u64> dp1(
        static_cast<std::size_t>(h + 1),
        0
    );

    dp0[0] = 1;

    for (int i = 0;
         i < h;
         ++i) {

        const u64 Li =
            Ldigits[
                static_cast<std::size_t>(i)
            ];

        std::vector<u64> next0(
            static_cast<std::size_t>(h + 1),
            0
        );

        std::vector<u64> next1(
            static_cast<std::size_t>(h + 1),
            0
        );

        for (int j = 0;
             j <= i;
             ++j) {

            /*
             * Incoming borrow = 0.
             *
             * y_i <= L_i:
             *   L_i+1 choices, no new borrow.
             *
             * y_i > L_i:
             *   p-L_i-1 choices, new borrow.
             */
            if (dp0[j] != 0) {
                next0[j] +=
                    dp0[j] *
                    (Li + 1);

                next1[j + 1] +=
                    dp0[j] *
                    (
                        static_cast<u64>(p) -
                        Li -
                        1
                    );
            }

            /*
             * Incoming borrow = 1.
             *
             * y_i <= L_i-1:
             *   L_i choices, borrow disappears.
             *
             * y_i >= L_i:
             *   p-L_i choices, borrow continues.
             */
            if (dp1[j] != 0) {
                next0[j] +=
                    dp1[j] *
                    Li;

                next1[j + 1] +=
                    dp1[j] *
                    (
                        static_cast<u64>(p) -
                        Li
                    );
            }
        }

        dp0.swap(next0);
        dp1.swap(next1);
    }

    /*
     * B_h(j) counts both possible final borrow states.
     */
    std::vector<u64> result(
        static_cast<std::size_t>(h + 1),
        0
    );

    for (int j = 0;
         j <= h;
         ++j) {
        result[j] =
            dp0[j] +
            dp1[j];
    }

    return result;
}

/* =========================================================
 * Direct local coefficient formula
 *
 * For:
 *
 *   L = sum_i L_i p^i
 *
 * and 1 <= k <= r:
 *
 *   A_{r,L}(k)
 *
 * is obtained by summing over the highest digit h where
 * y_h > L_h.
 *
 * valuation = (r-h) + B
 *
 * where B is the lower borrow count.
 * ========================================================= */

static std::vector<u64>
local_coefficient_formula(
    u64 L,
    int p,
    int r
) {
    std::vector<u64> result(
        static_cast<std::size_t>(r + 1),
        0
    );

    const auto Ldigits =
        digits_of(
            L,
            p,
            r
        );

    for (int h = 0;
         h < r;
         ++h) {

        const u64 Lh =
            Ldigits[
                static_cast<std::size_t>(h)
            ];

        /*
         * y_h may be any digit strictly above L_h.
         */
        const u64 high_choices =
            static_cast<u64>(p) -
            Lh -
            1;

        if (high_choices == 0) {
            continue;
        }

        const auto lower =
            lower_borrow_distribution(
                Ldigits,
                p,
                h
            );

        /*
         * The borrow at digit h itself is guaranteed
         * because y_h > L_h.
         *
         * Every higher digit  h+1,...,r-1 also carries.
         *
         * Therefore the fixed contribution is r-h.
         */
        const int fixed =
            r - h;

        for (int j = 0;
             j <= h;
             ++j) {

            const u64 count =
                lower[
                    static_cast<
                        std::size_t
                    >(j)
                ];

            if (count == 0) {
                continue;
            }

            const int valuation =
                fixed + j;

            result[
                static_cast<
                    std::size_t
                >(valuation)
            ] +=
                high_choices *
                count;
        }
    }

    while (
        result.size() > 1 &&
        result.back() == 0
    ) {
        result.pop_back();
    }

    return result;
}

/* =========================================================
 * Exact local histogram
 * ========================================================= */

static std::vector<u64>
exact_local_histogram(
    u64 L,
    int p,
    int r
) {
    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 n =
        power + L;

    std::vector<u64> result(
        static_cast<std::size_t>(r + 1),
        0
    );

    for (
        u64 y = L + 1;
        y < power;
        ++y
    ) {
        const u64 v =
            binomial_valuation(
                n,
                y,
                p
            );

        ++result[
            static_cast<std::size_t>(v)
        ];
    }

    while (
        result.size() > 1 &&
        result.back() == 0
    ) {
        result.pop_back();
    }

    return result;
}

/* =========================================================
 * Compare vectors
 * ========================================================= */

static bool equal_histogram(
    const std::vector<u64>& a,
    const std::vector<u64>& b
) {
    const std::size_t n =
        std::max(
            a.size(),
            b.size()
        );

    for (
        std::size_t i = 0;
        i < n;
        ++i
    ) {
        const u64 x =
            i < a.size()
                ? a[i]
                : 0;

        const u64 y =
            i < b.size()
                ? b[i]
                : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
}

/* =========================================================
 * Add scaled histogram
 * ========================================================= */

static void add_scaled_histogram(
    std::vector<u64>& dst,
    const std::vector<u64>& src,
    u64 scale
) {
    if (scale == 0) {
        return;
    }

    if (
        dst.size() <
        src.size()
    ) {
        dst.resize(
            src.size(),
            0
        );
    }

    for (
        std::size_t k = 0;
        k < src.size();
        ++k
    ) {
        dst[k] +=
            src[k] *
            scale;
    }

    while (
        dst.size() > 1 &&
        dst.back() == 0
    ) {
        dst.pop_back();
    }
}

/* =========================================================
 * Global coefficient formula
 *
 * A_0 = product_i(m_i+1)
 *
 * For k>=1:
 *
 * A_k =
 *   sum_r C_r * A_{r,L_r}(k)
 *
 * where
 *
 * C_r = m_r * product_{i>r}(m_i+1)
 * ========================================================= */

static std::vector<u64>
global_coefficient_formula(
    u64 m,
    int p
) {
    /*
     * Extract actual digits.
     */
    std::vector<u64> digits;

    u64 n = m;
    const u64 prime =
        static_cast<u64>(p);

    if (n == 0) {
        digits.push_back(0);
    } else {
        while (n > 0) {
            digits.push_back(
                n % prime
            );

            n /= prime;
        }
    }

    std::vector<u64> result(
        digits.size() + 1,
        0
    );

    /*
     * Constant term = MISS count.
     */
    u64 miss = 1;

    for (u64 d : digits) {
        miss *= d + 1;
    }

    result[0] = miss;

    /*
     * Suffix product:
     *
     * product_{i>r}(m_i+1)
     */
    std::vector<u64> suffix(
        digits.size(),
        1
    );

    u64 product = 1;

    for (
        std::size_t i =
            digits.size();
        i > 0;
        --i
    ) {
        const std::size_t r =
            i - 1;

        suffix[r] =
            product;

        product *=
            digits[r] + 1;
    }

    /*
     * Every digit position r>=1 with m_r>0 contributes
     * C_r identical local profiles.
     */
    for (
        int r = 1;
        r < static_cast<int>(
            digits.size()
        );
        ++r
    ) {
        const u64 mr =
            digits[
                static_cast<std::size_t>(r)
            ];

        if (mr == 0) {
            continue;
        }

        const u64 C =
            mr *
            suffix[
                static_cast<std::size_t>(r)
            ];

        const u64 L =
            m %
            prime_power(
                p,
                r
            );

        const auto local =
            local_coefficient_formula(
                L,
                p,
                r
            );

        add_scaled_histogram(
            result,
            local,
            C
        );
    }

    while (
        result.size() > 1 &&
        result.back() == 0
    ) {
        result.pop_back();
    }

    return result;
}

/* =========================================================
 * Histogram sum
 * ========================================================= */

static u64 histogram_sum(
    const std::vector<u64>& h
) {
    u64 result = 0;

    for (u64 x : h) {
        result += x;
    }

    return result;
}

/* =========================================================
 * Verify one case
 * ========================================================= */

static bool verify_case(
    u64 m,
    int p,
    bool exact_check,
    bool print_failure,
    Stats& stats
) {
    ++stats.tested;

    bool ok = true;

    const auto global_dp =
        global_kummer_dp(
            m,
            p
        );

    const auto formula =
        global_coefficient_formula(
            m,
            p
        );

    /*
     * --------------------------------------------------
     * Global coefficient formula vs global Kummer DP
     * --------------------------------------------------
     */

    if (
        !equal_histogram(
            global_dp,
            formula
        )
    ) {
        ok = false;

        ++stats.global_formula_failures;

        if (print_failure) {
            std::cout
                << "FAIL_GLOBAL_FORMULA "
                << "m=" << m
                << " p=" << p
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Total number of t values
     * --------------------------------------------------
     */

    if (
        histogram_sum(global_dp) !=
        m + 1
    ) {
        ok = false;

        ++stats.total_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TOTAL "
                << "m=" << m
                << " p=" << p
                << " actual="
                << histogram_sum(global_dp)
                << " expected="
                << (m + 1)
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Lucas layer
     * --------------------------------------------------
     */

    const u64 expected_zero =
        miss_count(
            m,
            p
        );

    const u64 dp_zero =
        global_dp.empty()
            ? 0
            : global_dp[0];

    const u64 formula_zero =
        formula.empty()
            ? 0
            : formula[0];

    if (
        dp_zero != expected_zero ||
        formula_zero != expected_zero
    ) {
        ok = false;

        ++stats.zero_layer_failures;

        if (print_failure) {
            std::cout
                << "FAIL_ZERO_LAYER "
                << "m=" << m
                << " p=" << p
                << " dp="
                << dp_zero
                << " formula="
                << formula_zero
                << " expected="
                << expected_zero
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Exact global enumeration for small m.
     * --------------------------------------------------
     */

    if (exact_check) {
        const auto exact =
            exact_global_histogram(
                m,
                p
            );

        if (
            !equal_histogram(
                exact,
                global_dp
            )
        ) {
            ok = false;

            ++stats.exact_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_EXACT "
                    << "m=" << m
                    << " p=" << p
                    << "\n";
            }
        }

        /*
         * ------------------------------------------------
         * Independently test every local coefficient
         * formula for every r present in m.
         * ------------------------------------------------
         */

        std::vector<u64> digits;

        u64 n = m;
        const u64 prime =
            static_cast<u64>(p);

        if (n == 0) {
            digits.push_back(0);
        } else {
            while (n > 0) {
                digits.push_back(
                    n % prime
                );

                n /= prime;
            }
        }

        for (
            int r = 1;
            r < static_cast<int>(
                digits.size()
            );
            ++r
        ) {
            const u64 L =
                m %
                prime_power(
                    p,
                    r
                );

            const auto local_formula =
                local_coefficient_formula(
                    L,
                    p,
                    r
                );

            const auto local_exact =
                exact_local_histogram(
                    L,
                    p,
                    r
                );

            if (
                !equal_histogram(
                    local_formula,
                    local_exact
                )
            ) {
                ok = false;

                ++stats.local_formula_failures;

                if (print_failure) {
                    std::cout
                        << "FAIL_LOCAL_FORMULA "
                        << "m=" << m
                        << " p=" << p
                        << " r=" << r
                        << " L=" << L
                        << "\n";
                }

                break;
            }
        }
    }

    if (!ok) {
        ++stats.failed;
    }

    return ok;
}

/* =========================================================
 * Build m from digits
 * ========================================================= */

static bool build_from_digits(
    const std::vector<u64>& digits,
    int p,
    u64& m
) {
    m = 0;

    u64 power = 1;

    for (u64 d : digits) {
        if (
            d != 0 &&
            power >
                UINT64_MAX / d
        ) {
            return false;
        }

        const u64 term =
            d * power;

        if (
            m >
            UINT64_MAX - term
        ) {
            return false;
        }

        m += term;

        if (
            power >
            UINT64_MAX /
                static_cast<u64>(p)
        ) {
            return false;
        }

        power *=
            static_cast<u64>(p);
    }

    return true;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 337\n"
        << "DIRECT COEFFICIENT FORMULA FOR KUMMER LAYERS\n"
        << "DO LOWER-BORROW DISTRIBUTIONS GIVE EVERY GLOBAL VALUATION COEFFICIENT?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small m.
     *
     * Exact individual binomial valuations are checked.
     * --------------------------------------------------
     */

    constexpr u64 SMALL_M = 500;

    for (int p : PRIMES) {
        for (u64 m = 0;
             m <= SMALL_M;
             ++m) {

            verify_case(
                m,
                p,
                true,
                stats.failed < 10,
                stats
            );
        }
    }

    /*
     * --------------------------------------------------
     * Phase 2:
     * Random arbitrary m.
     *
     * No t-enumeration.
     *
     * Only the digit-level formula and independent
     * global Kummer DP are evaluated.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x337337337ULL
    );

    constexpr u64 RANDOM_CASES = 15000;

    for (
        u64 i = 0;
        i < RANDOM_CASES;
        ++i
    ) {
        const int p =
            PRIMES[
                rng() % PRIMES.size()
            ];

        const u64 m =
            rng() %
            1000000000000ULL;

        verify_case(
            m,
            p,
            false,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 3:
     * Targeted digit structures.
     * --------------------------------------------------
     */

    for (int p : PRIMES) {
        for (int length = 1;
             length <= 14;
             ++length) {

            std::vector<
                std::vector<u64>
            > patterns;

            /*
             * all zero
             */
            patterns.emplace_back(
                length,
                0
            );

            /*
             * all maximum
             */
            patterns.emplace_back(
                length,
                static_cast<u64>(
                    p - 1
                )
            );

            /*
             * alternating 0/max
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    if (i & 1) {
                        d[i] =
                            static_cast<u64>(
                                p - 1
                            );
                    }
                }

                patterns.push_back(d);
            }

            /*
             * alternating max/0
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    if (!(i & 1)) {
                        d[i] =
                            static_cast<u64>(
                                p - 1
                            );
                    }
                }

                patterns.push_back(d);
            }

            /*
             * all ones
             */
            patterns.emplace_back(
                length,
                1
            );

            /*
             * increasing digits
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    d[i] =
                        static_cast<u64>(
                            i % p
                        );
                }

                patterns.push_back(d);
            }

            /*
             * Sparse low digits.
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = 0;
                     i < length / 2;
                     ++i) {
                    d[i] = 1;
                }

                patterns.push_back(d);
            }

            /*
             * Sparse high digits.
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = length / 2;
                     i < length;
                     ++i) {
                    d[i] = 1;
                }

                patterns.push_back(d);
            }

            for (const auto& digits :
                 patterns) {

                u64 m = 0;

                if (
                    !build_from_digits(
                        digits,
                        p,
                        m
                    )
                ) {
                    continue;
                }

                ++stats.targeted;

                verify_case(
                    m,
                    p,
                    m <= SMALL_M,
                    stats.failed < 10,
                    stats
                );
            }
        }
    }

    /*
     * --------------------------------------------------
     * Results
     * --------------------------------------------------
     */

    std::cout
        << "============================\n"
        << "tested="
        << stats.tested
        << "\n"
        << "failed="
        << stats.failed
        << "\n"
        << "local_formula_failures="
        << stats.local_formula_failures
        << "\n"
        << "global_formula_failures="
        << stats.global_formula_failures
        << "\n"
        << "exact_failures="
        << stats.exact_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "zero_layer_failures="
        << stats.zero_layer_failures
        << "\n"
        << "targeted="
        << stats.targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 337\n";

    return 0;
}
