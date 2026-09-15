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

    u64 histogram_failures = 0;
    u64 individual_failures = 0;
    u64 total_failures = 0;
    u64 maxval_failures = 0;
    u64 recurrence_failures = 0;
};

/* =========================================================
 * p-adic valuation of n!
 *
 * v_p(n!) = floor(n/p) + floor(n/p^2) + ...
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
 * Exact Legendre valuation
 *
 * v_p(C(m,t))
 * ========================================================= */

static u64 binomial_valuation(
    u64 m,
    u64 t,
    int p
) {
    return
        factorial_valuation(m, p)
        - factorial_valuation(t, p)
        - factorial_valuation(m - t, p);
}

/* =========================================================
 * Base-p digits of m
 * ========================================================= */

static std::vector<u64>
digits_of(
    u64 m,
    int p
) {
    std::vector<u64> digits;

    if (m == 0) {
        digits.push_back(0);
        return digits;
    }

    const u64 prime =
        static_cast<u64>(p);

    while (m > 0) {
        digits.push_back(
            m % prime
        );

        m /= prime;
    }

    return digits;
}

/* =========================================================
 * Borrow count for subtracting t from m in base p
 *
 * This is the Kummer carry number.
 *
 * Returns:
 *   - number of borrows
 *   - UINT64_MAX if t > m
 * ========================================================= */

static u64
subtraction_borrow_count(
    u64 m,
    u64 t,
    int p
) {
    if (t > m) {
        return UINT64_MAX;
    }

    const u64 prime =
        static_cast<u64>(p);

    u64 borrow = 0;
    u64 count = 0;

    while (m > 0 || t > 0 || borrow > 0) {
        const u64 md =
            m % prime;

        const u64 td =
            t % prime;

        const u64 value =
            md - td - borrow;

        /*
         * Unsigned underflow means the subtraction
         * requires a borrow.
         */
        if (
            md < td + borrow
        ) {
            borrow = 1;
            ++count;
        } else {
            borrow = 0;
        }

        (void)value;

        m /= prime;
        t /= prime;
    }

    /*
     * For t <= m this must finish with no borrow.
     */
    if (borrow != 0) {
        return UINT64_MAX;
    }

    return count;
}

/* =========================================================
 * Carry/Borrow DP
 *
 * dp[borrow][k] =
 * number of partial digit choices producing
 * k borrows and current borrow state.
 *
 * At each digit:
 *
 *   m_i - t_i - borrow_in >= 0
 *
 *       => borrow_out = 0
 *
 * otherwise
 *
 *   m_i - t_i - borrow_in + p >= 0
 *
 *       => borrow_out = 1
 *
 * Every outgoing borrow contributes one to the
 * Kummer valuation.
 *
 * Final borrow must be zero.
 * ========================================================= */

static std::vector<u64>
kummer_histogram_dp(
    u64 m,
    int p
) {
    const auto digits =
        digits_of(
            m,
            p
        );

    /*
     * The maximum number of borrows cannot exceed
     * the number of base-p digits.
     */
    const std::size_t L =
        digits.size();

    std::vector<u64> dp0(
        L + 1,
        0
    );

    std::vector<u64> dp1(
        L + 1,
        0
    );

    dp0[0] = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        std::vector<u64> next0(
            L + 1,
            0
        );

        std::vector<u64> next1(
            L + 1,
            0
        );

        const u64 md =
            digits[i];

        for (std::size_t k = 0;
             k <= i;
             ++k) {

            /*
             * Current borrow = 0.
             */
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

            /*
             * Current borrow = 1.
             */
            if (dp1[k] != 0) {
                for (
                    u64 td = 0;
                    td < static_cast<u64>(p);
                    ++td
                ) {
                    if (
                        md >= td + 1
                    ) {
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

    /*
     * Only final borrow = 0 represents 0 <= t <= m.
     */
    return dp0;
}

/* =========================================================
 * Exhaustive histogram using Legendre
 *
 * Used only for small m.
 * ========================================================= */

static std::vector<u64>
legendre_histogram(
    u64 m,
    int p
) {
    const auto digits =
        digits_of(
            m,
            p
        );

    const std::size_t L =
        digits.size();

    std::vector<u64> histogram(
        L + 1,
        0
    );

    for (u64 t = 0; t <= m; ++t) {
        const u64 value =
            binomial_valuation(
                m,
                t,
                p
            );

        if (
            value >= histogram.size()
        ) {
            /*
             * Should never happen.
             */
            return {};
        }

        ++histogram[value];
    }

    return histogram;
}

/* =========================================================
 * Remove trailing zero buckets
 *
 * Makes histogram equality cleaner.
 * ========================================================= */

static void
trim_histogram(
    std::vector<u64>& histogram
) {
    while (
        !histogram.empty() &&
        histogram.back() == 0
    ) {
        histogram.pop_back();
    }
}

/* =========================================================
 * Histogram equality
 * ========================================================= */

static bool
same_histogram(
    const std::vector<u64>& a,
    const std::vector<u64>& b
) {
    std::vector<u64> x = a;
    std::vector<u64> y = b;

    trim_histogram(x);
    trim_histogram(y);

    return x == y;
}

/* =========================================================
 * Sum histogram
 * ========================================================= */

static u64
histogram_total(
    const std::vector<u64>& histogram
) {
    u64 total = 0;

    for (u64 x : histogram) {
        total += x;
    }

    return total;
}

/* =========================================================
 * Highest nonzero valuation
 * ========================================================= */

static u64
histogram_max(
    const std::vector<u64>& histogram
) {
    for (
        std::size_t i =
            histogram.size();
        i > 0;
        --i
    ) {
        if (histogram[i - 1] != 0) {
            return static_cast<u64>(
                i - 1
            );
        }
    }

    return 0;
}

/* =========================================================
 * Verify individual Kummer correspondence
 *
 * For every t:
 *
 *   Legendre valuation
 *       ==
 *   subtraction borrow count.
 * ========================================================= */

static bool
verify_individual_borrows(
    u64 m,
    int p,
    bool print_failure,
    Stats& stats
) {
    for (u64 t = 0; t <= m; ++t) {
        const u64 exact =
            binomial_valuation(
                m,
                t,
                p
            );

        const u64 carries =
            subtraction_borrow_count(
                m,
                t,
                p
            );

        if (
            exact != carries
        ) {
            ++stats.individual_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_INDIVIDUAL "
                    << "m=" << m
                    << " p=" << p
                    << " t=" << t
                    << " exact=" << exact
                    << " borrows="
                    << carries
                    << "\n";
            }

            return false;
        }
    }

    return true;
}

/* =========================================================
 * Verify one case
 * ========================================================= */

static bool
verify_case(
    u64 m,
    int p,
    bool brute_check,
    bool individual_check,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    /*
     * --------------------------------------------------
     * Carry DP histogram
     * --------------------------------------------------
     */

    const auto dp =
        kummer_histogram_dp(
            m,
            p
        );

    /*
     * --------------------------------------------------
     * Total coefficient count
     *
     * Every t in [0,m] occurs exactly once.
     * --------------------------------------------------
     */

    const u64 total =
        histogram_total(
            dp
        );

    if (
        total != m + 1
    ) {
        ok = false;

        ++stats.total_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TOTAL "
                << "m=" << m
                << " p=" << p
                << " total=" << total
                << " expected="
                << (m + 1)
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Maximum valuation sanity bound
     * --------------------------------------------------
     */

    const auto digits =
        digits_of(
            m,
            p
        );

    const u64 maximum =
        histogram_max(
            dp
        );

    if (
        maximum >
        static_cast<u64>(
            digits.size()
        )
    ) {
        ok = false;

        ++stats.maxval_failures;

        if (print_failure) {
            std::cout
                << "FAIL_MAX "
                << "m=" << m
                << " p=" << p
                << " max="
                << maximum
                << " digits="
                << digits.size()
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Exhaustive Legendre comparison
     * --------------------------------------------------
     */

    if (brute_check) {
        const auto exact =
            legendre_histogram(
                m,
                p
            );

        if (
            !same_histogram(
                dp,
                exact
            )
        ) {
            ok = false;

            ++stats.histogram_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_HISTOGRAM "
                    << "m=" << m
                    << " p=" << p
                    << "\n";

                const std::size_t n =
                    std::max(
                        dp.size(),
                        exact.size()
                    );

                for (std::size_t k = 0;
                     k < n;
                     ++k) {

                    const u64 a =
                        k < dp.size()
                            ? dp[k]
                            : 0;

                    const u64 b =
                        k < exact.size()
                            ? exact[k]
                            : 0;

                    if (a != b) {
                        std::cout
                            << "  k="
                            << k
                            << " dp="
                            << a
                            << " exact="
                            << b
                            << "\n";
                    }
                }
            }
        }

        /*
         * Individual valuation check.
         */
        if (individual_check) {
            if (
                !verify_individual_borrows(
                    m,
                    p,
                    print_failure,
                    stats
                )
            ) {
                ok = false;
            }
        }
    }

    /*
     * --------------------------------------------------
     * Direct Lucas layer check:
     *
     * valuation 0 count must equal
     *
     * product_i (m_i+1).
     * --------------------------------------------------
     */

    u64 lucas_count = 1;

    for (u64 d : digits) {
        lucas_count *= d + 1;
    }

    const u64 dp_zero =
        dp.empty()
            ? 0
            : dp[0];

    if (
        dp_zero != lucas_count
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_LUCAS_LAYER "
                << "m=" << m
                << " p=" << p
                << " dp0="
                << dp_zero
                << " expected="
                << lucas_count
                << "\n";
        }
    }

    if (!ok) {
        ++stats.failed;
    }

    return ok;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 331\n"
        << "KUMMER CARRY FILTRATION\n"
        << "DO BORROW-COUNT DIGIT DP REPRODUCE THE FULL p-ADIC VALUATION HISTOGRAM?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small m.
     *
     * Full comparison against Legendre valuation,
     * plus individual borrow verification.
     * --------------------------------------------------
     */

    for (int p : PRIMES) {
        for (u64 m = 0; m <= 1000; ++m) {
            ++stats.tested;

            verify_case(
                m,
                p,
                true,
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
     * Only digit DP is used.
     *
     * No enumeration over t.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x331331331ULL
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
            rng() % 1000000000000ULL;

        ++stats.tested;

        verify_case(
            m,
            p,
            false,
            false,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 3:
     * Targeted digit structures.
     *
     * Build m directly from base-p digits.
     * --------------------------------------------------
     */

    u64 targeted = 0;

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
             * alternating 0 / max
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
             * alternating max / 0
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
             * Sparse digits:
             * only high positions nonzero.
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    if (
                        i >= length / 2
                    ) {
                        d[i] = 1;
                    }
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

            for (const auto& ds :
                 patterns) {

                u64 m = 0;
                u64 power = 1;
                bool overflow = false;

                for (u64 digit : ds) {
                    if (
                        digit != 0 &&
                        power >
                            UINT64_MAX / digit
                    ) {
                        overflow = true;
                        break;
                    }

                    const u64 term =
                        digit * power;

                    if (
                        m >
                        UINT64_MAX - term
                    ) {
                        overflow = true;
                        break;
                    }

                    m += term;

                    if (
                        power >
                        UINT64_MAX /
                            static_cast<u64>(p)
                    ) {
                        overflow = true;
                        break;
                    }

                    power *=
                        static_cast<u64>(p);
                }

                if (overflow) {
                    continue;
                }

                ++targeted;
                ++stats.tested;

                /*
                 * Small targeted cases get the full
                 * individual check as well.
                 */
                verify_case(
                    m,
                    p,
                    m <= 5000,
                    m <= 5000,
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
        << "histogram_failures="
        << stats.histogram_failures
        << "\n"
        << "individual_failures="
        << stats.individual_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "maxval_failures="
        << stats.maxval_failures
        << "\n"
        << "recurrence_failures="
        << stats.recurrence_failures
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 331\n";

    return 0;
}
