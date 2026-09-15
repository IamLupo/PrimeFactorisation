#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <vector>

using u64 = std::uint64_t;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

struct Gap {
    u64 left;
    u64 right;
    int type;
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 point_failures = 0;
    u64 histogram_failures = 0;
    u64 digit_formula_failures = 0;
    u64 length_failures = 0;
    u64 total_failures = 0;
    u64 maximum_failures = 0;
};

/* =========================================================
 * Base-p digits
 * ========================================================= */

static std::vector<u64> digits_of(
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
 * Lucas MISS predicate
 * ========================================================= */

static bool is_miss(
    u64 x,
    u64 m,
    int p
) {
    const u64 prime =
        static_cast<u64>(p);

    while (x > 0 || m > 0) {
        const u64 xd =
            x % prime;

        const u64 md =
            m % prime;

        if (xd > md) {
            return false;
        }

        x /= prime;
        m /= prime;
    }

    return true;
}

/* =========================================================
 * p^r
 * ========================================================= */

static u64 prime_power(
    int p,
    int r
) {
    u64 result = 1;

    for (int i = 0; i < r; ++i) {
        result *=
            static_cast<u64>(p);
    }

    return result;
}

/* =========================================================
 * m mod p^r
 * ========================================================= */

static u64 lower_block(
    u64 m,
    int p,
    int r
) {
    if (r == 0) {
        return 0;
    }

    return m %
        prime_power(
            p,
            r
        );
}

/* =========================================================
 * Gap type from the preceding MISS value
 * ========================================================= */

static int gap_type_from_left_miss(
    u64 x,
    u64 m,
    int p
) {
    const auto digits =
        digits_of(
            m,
            p
        );

    const u64 prime =
        static_cast<u64>(p);

    for (
        std::size_t r = 0;
        r < digits.size();
        ++r
    ) {
        const u64 xd =
            x % prime;

        if (xd < digits[r]) {
            return static_cast<int>(r);
        }

        x /= prime;
    }

    return -1;
}

/* =========================================================
 * Find all HIT gaps
 * ========================================================= */

static std::vector<Gap>
find_typed_gaps(
    u64 m,
    int p
) {
    std::vector<Gap> result;

    bool have_previous_miss = false;
    u64 previous_miss = 0;

    for (u64 x = 0; x <= m; ++x) {
        if (!is_miss(x, m, p)) {
            continue;
        }

        if (have_previous_miss) {
            if (x > previous_miss + 1) {
                const int type =
                    gap_type_from_left_miss(
                        previous_miss,
                        m,
                        p
                    );

                result.push_back(
                    Gap{
                        previous_miss + 1,
                        x - 1,
                        type
                    }
                );
            }
        }

        previous_miss = x;
        have_previous_miss = true;
    }

    return result;
}

/* =========================================================
 * Actual histogram inside a gap
 *
 * histogram[k] =
 * number of t in the gap with
 * v_p(C(m,t)) = k
 * ========================================================= */

static std::map<u64, u64>
actual_gap_histogram(
    const Gap& gap,
    u64 m,
    int p
) {
    std::map<u64, u64> histogram;

    for (
        u64 t = gap.left;
        t <= gap.right;
        ++t
    ) {
        const u64 v =
            binomial_valuation(
                m,
                t,
                p
            );

        ++histogram[v];
    }

    return histogram;
}

/* =========================================================
 * Closed digit-only histogram
 *
 * For:
 *
 *   L = m mod p^r
 *
 * and k = 1,...,r:
 *
 *   A_k =
 *     (p - L_{r-k} - 1) * p^(r-k)
 *
 * ========================================================= */

static std::map<u64, u64>
closed_gap_histogram(
    u64 m,
    int p,
    int r
) {
    std::map<u64, u64> histogram;

    if (r <= 0) {
        return histogram;
    }

    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 L =
        lower_block(
            m,
            p,
            r
        );

    const auto digits =
        digits_of(
            L,
            p
        );

    /*
     * Ensure there are r low digits.
     */
    const std::size_t needed =
        static_cast<std::size_t>(r);

    /*
     * For k=1:
     *   h=r-1
     *
     * For k=r:
     *   h=0
     */
    for (int k = 1;
         k <= r;
         ++k) {

        const int h =
            r - k;

        u64 Lh = 0;

        if (
            static_cast<std::size_t>(h) <
            digits.size()
        ) {
            Lh = digits[h];
        }

        const u64 choices =
            static_cast<u64>(p) -
            Lh -
            1;

        if (choices == 0) {
            continue;
        }

        const u64 free_lower =
            prime_power(
                p,
                h
            );

        const u64 count =
            choices *
            free_lower;

        histogram[
            static_cast<u64>(k)
        ] += count;
    }

    (void)power;
    (void)needed;

    return histogram;
}

/* =========================================================
 * Sum a histogram
 * ========================================================= */

static u64
histogram_total(
    const std::map<u64, u64>& histogram
) {
    u64 total = 0;

    for (const auto& entry :
         histogram) {
        total += entry.second;
    }

    return total;
}

/* =========================================================
 * Weighted valuation sum
 *
 * sum k*A_k
 * ========================================================= */

static u64
histogram_weighted(
    const std::map<u64, u64>& histogram
) {
    u64 total = 0;

    for (const auto& entry :
         histogram) {
        total +=
            entry.first *
            entry.second;
    }

    return total;
}

/* =========================================================
 * Histogram printer
 * ========================================================= */

static void print_histogram(
    const std::map<u64, u64>& histogram
) {
    std::cout << "{";

    bool first = true;

    for (const auto& entry :
         histogram) {

        if (!first) {
            std::cout << ", ";
        }

        first = false;

        std::cout
            << entry.first
            << ":"
            << entry.second;
    }

    std::cout << "}";
}

/* =========================================================
 * Verify one gap
 * ========================================================= */

static bool verify_gap(
    const Gap& gap,
    u64 m,
    int p,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    const int r =
        gap.type;

    if (r <= 0) {
        ok = false;

        ++stats.length_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TYPE "
                << "m=" << m
                << " p=" << p
                << " type=" << r
                << "\n";
        }

        return false;
    }

    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 L =
        lower_block(
            m,
            p,
            r
        );

    /*
     * --------------------------------------------------
     * Actual gap length
     *
     * G = p^r - L - 1
     * --------------------------------------------------
     */

    const u64 actual_length =
        gap.right -
        gap.left +
        1;

    const u64 expected_length =
        power -
        L -
        1;

    if (
        actual_length !=
        expected_length
    ) {
        ok = false;

        ++stats.length_failures;

        if (print_failure) {
            std::cout
                << "FAIL_LENGTH "
                << "m=" << m
                << " p=" << p
                << " r=" << r
                << " L=" << L
                << " actual="
                << actual_length
                << " expected="
                << expected_length
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Actual vs closed histogram
     * --------------------------------------------------
     */

    const auto actual =
        actual_gap_histogram(
            gap,
            m,
            p
        );

    const auto closed =
        closed_gap_histogram(
            m,
            p,
            r
        );

    if (actual != closed) {
        ok = false;

        ++stats.histogram_failures;

        if (print_failure) {
            std::cout
                << "FAIL_HISTOGRAM "
                << "m=" << m
                << " p=" << p
                << " r=" << r
                << " L=" << L
                << "\n";

            std::cout
                << "actual=";

            print_histogram(
                actual
            );

            std::cout
                << "\nclosed=";

            print_histogram(
                closed
            );

            std::cout
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Histogram total must equal gap length
     * --------------------------------------------------
     */

    const u64 actual_total =
        histogram_total(
            actual
        );

    const u64 closed_total =
        histogram_total(
            closed
        );

    if (
        actual_total !=
        actual_length ||
        closed_total !=
        expected_length
    ) {
        ok = false;

        ++stats.total_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TOTAL "
                << "m=" << m
                << " p=" << p
                << " r=" << r
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Maximum valuation must be <= r
     * --------------------------------------------------
     */

    if (!actual.empty()) {
        const u64 max_value =
            actual.rbegin()->first;

        if (
            max_value >
            static_cast<u64>(r)
        ) {
            ok = false;

            ++stats.maximum_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_MAX "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " max="
                    << max_value
                    << "\n";
            }
        }
    }

    /*
     * --------------------------------------------------
     * Pointwise digit formula
     *
     * For each t:
     *
     *   k = v_p(C(m,t))
     *
     * and h=r-k.
     *
     * We directly verify that the most significant
     * differing digit between y=t mod p^r and L is h.
     * --------------------------------------------------
     */

    const auto Ldigits =
        digits_of(
            L,
            p
        );

    for (
        u64 t = gap.left;
        t <= gap.right;
        ++t
    ) {
        const u64 y =
            t % power;

        const u64 valuation =
            binomial_valuation(
                m,
                t,
                p
            );

        if (
            valuation == 0 ||
            valuation >
                static_cast<u64>(r)
        ) {
            ok = false;

            ++stats.point_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_POINT_RANGE "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " t=" << t
                    << " v="
                    << valuation
                    << "\n";
            }

            break;
        }

        const int h =
            r -
            static_cast<int>(
                valuation
            );

        const auto Ydigits =
            digits_of(
                y,
                p
            );

        /*
         * Find most significant digit where y != L.
         */
        int observed_h = -1;

        for (int i = r - 1;
             i >= 0;
             --i) {

            const u64 yd =
                static_cast<
                    u64
                >(
                    static_cast<
                        std::size_t
                    >(i) < Ydigits.size()
                        ? Ydigits[i]
                        : 0
                );

            const u64 ld =
                static_cast<
                    u64
                >(
                    static_cast<
                        std::size_t
                    >(i) < Ldigits.size()
                        ? Ldigits[i]
                        : 0
                );

            if (yd != ld) {
                observed_h = i;
                break;
            }
        }

        if (
            observed_h != h
        ) {
            ok = false;

            ++stats.point_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_DIGIT_POSITION "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " t=" << t
                    << " y=" << y
                    << " v="
                    << valuation
                    << " expected_h="
                    << h
                    << " observed_h="
                    << observed_h
                    << "\n";
            }

            break;
        }

        /*
         * Verify the differing digit is actually above L.
         */
        const u64 yd =
            static_cast<
                std::size_t
            >(h) < Ydigits.size()
                ? Ydigits[h]
                : 0;

        const u64 ld =
            static_cast<
                std::size_t
            >(h) < Ldigits.size()
                ? Ldigits[h]
                : 0;

        if (yd <= ld) {
            ok = false;

            ++stats.digit_formula_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_DIGIT_DIRECTION "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " t=" << t
                    << " h=" << h
                    << " y_h=" << yd
                    << " L_h=" << ld
                    << "\n";
            }

            break;
        }
    }

    if (!ok) {
        return false;
    }

    return true;
}

/* =========================================================
 * Verify one m,p
 * ========================================================= */

static bool verify_case(
    u64 m,
    int p,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    const auto gaps =
        find_typed_gaps(
            m,
            p
        );

    for (const Gap& gap : gaps) {
        if (
            !verify_gap(
                gap,
                m,
                p,
                print_failure,
                stats
            )
        ) {
            ok = false;
            break;
        }
    }

    if (!ok) {
        ++stats.failed;
    }

    return ok;
}

/* =========================================================
 * Build m from base-p digits
 * ========================================================= */

static bool build_from_digits(
    const std::vector<u64>& digits,
    int p,
    u64& m
) {
    m = 0;

    u64 power = 1;

    for (u64 digit : digits) {
        if (
            digit != 0 &&
            power >
                UINT64_MAX / digit
        ) {
            return false;
        }

        const u64 term =
            digit * power;

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
        << "START EXPERIMENT 333\n"
        << "CLOSED DIGIT HISTOGRAM OF LOCAL KUMMER GAPS\n"
        << "DO GAP VALUATIONS REDUCE TO A SINGLE DIGIT FORMULA?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small m.
     *
     * Every gap and every t in every gap.
     * --------------------------------------------------
     */

    for (int p : PRIMES) {
        for (u64 m = 0; m <= 750; ++m) {
            ++stats.tested;

            verify_case(
                m,
                p,
                stats.failed < 10,
                stats
            );
        }
    }

    /*
     * --------------------------------------------------
     * Phase 2:
     * Targeted digit patterns.
     * --------------------------------------------------
     */

    u64 targeted = 0;

    for (int p : PRIMES) {
        for (int length = 1;
             length <= 10;
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
             * A single high nonzero digit.
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                d[length - 1] =
                    static_cast<u64>(
                        p - 1
                    );

                patterns.push_back(d);
            }

            /*
             * Sparse mixed pattern.
             */
            {
                std::vector<u64> d(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    if (i % 3 == 0) {
                        d[i] = 1;
                    }

                    if (i % 4 == 1) {
                        d[i] =
                            static_cast<u64>(
                                p - 1
                            );
                    }
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

                /*
                 * Keep pointwise verification bounded.
                 */
                if (m > 200000) {
                    continue;
                }

                ++targeted;
                ++stats.tested;

                verify_case(
                    m,
                    p,
                    stats.failed < 10,
                    stats
                );
            }
        }
    }

    /*
     * --------------------------------------------------
     * Phase 3:
     * Random arbitrary m.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x333333333ULL
    );

    constexpr u64 RANDOM_CASES = 1500;

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
            rng() % 100000ULL;

        ++stats.tested;

        verify_case(
            m,
            p,
            stats.failed < 10,
            stats
        );
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
        << "point_failures="
        << stats.point_failures
        << "\n"
        << "histogram_failures="
        << stats.histogram_failures
        << "\n"
        << "digit_formula_failures="
        << stats.digit_formula_failures
        << "\n"
        << "length_failures="
        << stats.length_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "maximum_failures="
        << stats.maximum_failures
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 333\n";

    return 0;
}
