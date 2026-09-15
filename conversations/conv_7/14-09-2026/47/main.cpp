#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

struct Poly {
    std::vector<u64> c;
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 exact_failures = 0;
    u64 global_dp_failures = 0;
    u64 decomposition_failures = 0;
    u64 coefficient_sum_failures = 0;
    u64 zero_layer_failures = 0;
    u64 gap_count_failures = 0;
    u64 local_failures = 0;

    u64 targeted = 0;
};

/* =========================================================
 * Polynomial helpers
 * ========================================================= */

static void trim(
    Poly& p
) {
    while (
        p.c.size() > 1 &&
        p.c.back() == 0
    ) {
        p.c.pop_back();
    }

    if (p.c.empty()) {
        p.c.push_back(0);
    }
}

static bool equal_poly(
    const Poly& a,
    const Poly& b
) {
    const std::size_t n =
        std::max(
            a.c.size(),
            b.c.size()
        );

    for (std::size_t i = 0;
         i < n;
         ++i) {

        const u64 x =
            i < a.c.size()
                ? a.c[i]
                : 0;

        const u64 y =
            i < b.c.size()
                ? b.c[i]
                : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
}

static void add_scaled(
    Poly& dst,
    const Poly& src,
    u64 scale
) {
    if (scale == 0) {
        return;
    }

    if (
        dst.c.size() <
        src.c.size()
    ) {
        dst.c.resize(
            src.c.size(),
            0
        );
    }

    for (
        std::size_t i = 0;
        i < src.c.size();
        ++i
    ) {
        dst.c[i] +=
            src.c[i] * scale;
    }

    trim(dst);
}

static u64 coefficient_sum(
    const Poly& p
) {
    u64 result = 0;

    for (u64 x : p.c) {
        result += x;
    }

    return result;
}

static void print_poly(
    const Poly& p
) {
    std::cout << "[";

    for (
        std::size_t i = 0;
        i < p.c.size();
        ++i
    ) {
        if (i != 0) {
            std::cout << ",";
        }

        std::cout << p.c[i];
    }

    std::cout << "]";
}

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
    int p
) {
    std::vector<u64> digits;

    if (n == 0) {
        digits.push_back(0);
        return digits;
    }

    const u64 prime =
        static_cast<u64>(p);

    while (n > 0) {
        digits.push_back(
            n % prime
        );

        n /= prime;
    }

    return digits;
}

/* =========================================================
 * MISS count
 *
 * product_i (m_i+1)
 * ========================================================= */

static u64 miss_count(
    u64 m,
    int p
) {
    const auto digits =
        digits_of(
            m,
            p
        );

    u64 result = 1;

    for (u64 d : digits) {
        result *= d + 1;
    }

    return result;
}

/* =========================================================
 * Global exact polynomial
 *
 * Used only for small m.
 * ========================================================= */

static Poly exact_global_polynomial(
    u64 m,
    int p
) {
    Poly result;

    result.c.assign(
        1,
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
            result.c.size() <= v
        ) {
            result.c.resize(
                static_cast<std::size_t>(v + 1),
                0
            );
        }

        ++result.c[v];
    }

    trim(result);
    return result;
}

/* =========================================================
 * Global Kummer digit DP
 *
 * This computes the complete distribution of
 * v_p(C(m,t)) for all 0 <= t <= m.
 *
 * State:
 *
 *   borrow = current borrow in m-t
 *
 * valuation is the polynomial degree.
 *
 * We don't need an explicit y>m relation here because
 * digitwise subtraction is performed against the actual
 * upper bound m and only final borrow=0 is accepted.
 * ========================================================= */

static Poly global_kummer_dp(
    u64 m,
    int p
) {
    const auto mdigits =
        digits_of(
            m,
            p
        );

    const int length =
        static_cast<int>(
            mdigits.size()
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
            mdigits[
                static_cast<std::size_t>(i)
            ];

        for (int v = 0;
             v <= i;
             ++v) {

            if (dp0[v] != 0) {
                for (
                    u64 td = 0;
                    td <
                        static_cast<u64>(p);
                    ++td
                ) {
                    if (md >= td) {
                        next0[v] +=
                            dp0[v];
                    } else {
                        next1[v + 1] +=
                            dp0[v];
                    }
                }
            }

            if (dp1[v] != 0) {
                for (
                    u64 td = 0;
                    td <
                        static_cast<u64>(p);
                    ++td
                ) {
                    if (md >= td + 1) {
                        next0[v] +=
                            dp1[v];
                    } else {
                        next1[v + 1] +=
                            dp1[v];
                    }
                }
            }
        }

        dp0.swap(next0);
        dp1.swap(next1);
    }

    Poly result;

    result.c =
        dp0;

    trim(result);

    return result;
}

/* =========================================================
 * Local Kummer transfer DP
 *
 * This is the same local polynomial construction validated
 * in Experiment 335, implemented directly here.
 *
 * We count:
 *
 *   L < y < p^r
 *
 * and the valuation of
 *
 *   C(p^r+L, y).
 * ========================================================= */

static Poly local_kummer_polynomial(
    u64 L,
    int p,
    int r
) {
    /*
     * State:
     *
     * borrow:
     *   0,1
     *
     * relation:
     *   0 = equal
     *   1 = y > L
     *   2 = y < L
     *
     * index = borrow*3 + relation
     */

    constexpr int STATES = 6;

    std::vector<Poly> dp(
        STATES
    );

    dp[0].c = {1};

    const auto digits =
        digits_of(
            L,
            p
        );

    for (int i = 0;
         i < r;
         ++i) {

        std::vector<Poly> next(
            STATES
        );

        const u64 Li =
            i < static_cast<int>(
                digits.size()
            )
                ? digits[
                    static_cast<std::size_t>(i)
                ]
                : 0;

        for (int state = 0;
             state < STATES;
             ++state) {

            if (
                dp[state].c.size() == 1 &&
                dp[state].c[0] == 0
            ) {
                continue;
            }

            const int borrow =
                state / 3;

            const int relation =
                state % 3;

            for (
                u64 yi = 0;
                yi <
                    static_cast<u64>(p);
                ++yi
            ) {
                int new_relation =
                    relation;

                /*
                 * Higher digits dominate lower digits.
                 * Since we process low -> high, replace the
                 * relation whenever this digit differs.
                 */
                if (yi > Li) {
                    new_relation = 1;
                } else if (yi < Li) {
                    new_relation = 2;
                }

                int new_borrow = 0;
                int valuation_add = 0;

                if (
                    Li <
                    yi +
                    static_cast<u64>(borrow)
                ) {
                    new_borrow = 1;
                    valuation_add = 1;
                }

                const int new_state =
                    new_borrow * 3 +
                    new_relation;

                if (
                    next[new_state].c.size() <
                    dp[state].c.size() +
                    static_cast<std::size_t>(
                        valuation_add
                    )
                ) {
                    next[new_state].c.resize(
                        dp[state].c.size() +
                        static_cast<std::size_t>(
                            valuation_add
                        ),
                        0
                    );
                }

                for (
                    std::size_t v = 0;
                    v < dp[state].c.size();
                    ++v
                ) {
                    next[
                        new_state
                    ].c[
                        v +
                        static_cast<std::size_t>(
                            valuation_add
                        )
                    ] +=
                        dp[state].c[v];
                }
            }
        }

        dp.swap(next);

        for (Poly& poly : dp) {
            trim(poly);
        }
    }

    Poly result;

    result.c = {0};

    /*
     * relation=1 means y>L.
     *
     * Both borrow states are accepted because the leading
     * digit 1 of p^r+L absorbs a possible final borrow.
     */
    add_scaled(
        result,
        dp[1],
        1
    );

    add_scaled(
        result,
        dp[4],
        1
    );

    trim(result);

    return result;
}

/* =========================================================
 * p^r gap-count coefficient
 *
 * C_r = m_r * product_{i>r}(m_i+1)
 * ========================================================= */

static u64 gap_count_for_digit(
    const std::vector<u64>& digits,
    int r
) {
    if (
        r < 0 ||
        r >= static_cast<int>(
            digits.size()
        )
    ) {
        return 0;
    }

    u64 result =
        digits[
            static_cast<std::size_t>(r)
        ];

    for (
        int i = r + 1;
        i < static_cast<int>(
            digits.size()
        );
        ++i
    ) {
        result *=
            digits[
                static_cast<std::size_t>(i)
            ] + 1;
    }

    return result;
}

/* =========================================================
 * Global gap decomposition
 *
 * F(m,z) =
 *
 *   MISS
 *   +
 *   sum_r C_r K_{r,L_r}(z)
 * ========================================================= */

static Poly decomposed_global_polynomial(
    u64 m,
    int p,
    bool& local_ok
) {
    local_ok = true;

    const auto digits =
        digits_of(
            m,
            p
        );

    Poly result;

    result.c.assign(
        1,
        miss_count(
            m,
            p
        )
    );

    for (
        int r = 1;
        r < static_cast<int>(
            digits.size()
        );
        ++r
    ) {
        const u64 C =
            gap_count_for_digit(
                digits,
                r
            );

        if (C == 0) {
            continue;
        }

        const u64 L =
            m %
            prime_power(
                p,
                r
            );

        const Poly local =
            local_kummer_polynomial(
                L,
                p,
                r
            );

        /*
         * Local polynomial coefficient sum must equal
         * the length of a type-r gap.
         */
        const u64 expected_local =
            prime_power(
                p,
                r
            ) -
            L -
            1;

        if (
            coefficient_sum(local) !=
            expected_local
        ) {
            local_ok = false;
        }

        add_scaled(
            result,
            local,
            C
        );
    }

    trim(result);

    return result;
}

/* =========================================================
 * Verify one m,p
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

    /*
     * --------------------------------------------------
     * Global Kummer DP
     * --------------------------------------------------
     */

    const Poly global_dp =
        global_kummer_dp(
            m,
            p
        );

    /*
     * --------------------------------------------------
     * Decomposition
     * --------------------------------------------------
     */

    bool local_ok = true;

    const Poly decomposed =
        decomposed_global_polynomial(
            m,
            p,
            local_ok
        );

    if (!local_ok) {
        ok = false;

        ++stats.local_failures;

        if (print_failure) {
            std::cout
                << "FAIL_LOCAL_SUM "
                << "m=" << m
                << " p=" << p
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Compare decomposition with global Kummer DP
     * --------------------------------------------------
     */

    if (
        !equal_poly(
            global_dp,
            decomposed
        )
    ) {
        ok = false;

        ++stats.decomposition_failures;

        if (print_failure) {
            std::cout
                << "FAIL_DECOMPOSITION "
                << "m=" << m
                << " p=" << p
                << "\n";

            std::cout
                << "global_dp=";

            print_poly(
                global_dp
            );

            std::cout
                << "\ndecomposed=";

            print_poly(
                decomposed
            );

            std::cout
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Total coefficient count
     *
     * Must be m+1.
     * --------------------------------------------------
     */

    if (
        coefficient_sum(global_dp) !=
        m + 1
    ) {
        ok = false;

        ++stats.coefficient_sum_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TOTAL "
                << "m=" << m
                << " p=" << p
                << " total="
                << coefficient_sum(global_dp)
                << " expected="
                << (m + 1)
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Constant term must be Lucas MISS count.
     * --------------------------------------------------
     */

    const u64 misses =
        miss_count(
            m,
            p
        );

    const u64 global_zero =
        global_dp.c.empty()
            ? 0
            : global_dp.c[0];

    const u64 decomposed_zero =
        decomposed.c.empty()
            ? 0
            : decomposed.c[0];

    if (
        global_zero != misses ||
        decomposed_zero != misses
    ) {
        ok = false;

        ++stats.zero_layer_failures;

        if (print_failure) {
            std::cout
                << "FAIL_ZERO_LAYER "
                << "m=" << m
                << " p=" << p
                << " global="
                << global_zero
                << " decomposed="
                << decomposed_zero
                << " expected="
                << misses
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Exact enumeration for small cases.
     * --------------------------------------------------
     */

    if (exact_check) {
        const Poly exact =
            exact_global_polynomial(
                m,
                p
            );

        if (
            !equal_poly(
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

                std::cout
                    << "exact=";

                print_poly(
                    exact
                );

                std::cout
                    << "\nglobal_dp=";

                print_poly(
                    global_dp
                );

                std::cout
                    << "\n";
            }
        }
    }

    if (!ok) {
        ++stats.failed;
    }

    return ok;
}

/* =========================================================
 * Construct m from digits
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
        << "START EXPERIMENT 336\n"
        << "GLOBAL KUMMER POLYNOMIAL FROM LUCAS GAPS\n"
        << "DO LOCAL GAP POLYNOMIALS RECONSTRUCT THE FULL VALUATION DISTRIBUTION?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small global m.
     *
     * Exact global polynomial is independently computed.
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
     * No enumeration over t.
     *
     * Global Kummer DP is compared against the gap
     * decomposition.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x336336336ULL
    );

    constexpr u64 RANDOM_CASES = 10000;

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
             length <= 12;
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
        << "exact_failures="
        << stats.exact_failures
        << "\n"
        << "global_dp_failures="
        << stats.global_dp_failures
        << "\n"
        << "decomposition_failures="
        << stats.decomposition_failures
        << "\n"
        << "coefficient_sum_failures="
        << stats.coefficient_sum_failures
        << "\n"
        << "zero_layer_failures="
        << stats.zero_layer_failures
        << "\n"
        << "local_failures="
        << stats.local_failures
        << "\n"
        << "gap_count_failures="
        << stats.gap_count_failures
        << "\n"
        << "targeted="
        << stats.targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 336\n";

    return 0;
}
