#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <vector>

using u64 = std::uint64_t;

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

struct Interval {
    u64 left;
    u64 right;
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 histogram_failures = 0;
    u64 count_failures = 0;
    u64 weighted_failures = 0;
    u64 brute_failures = 0;
};

/* =========================================================
 * Base-b digits
 * ========================================================= */

static std::vector<u64> digits_of(
    u64 m,
    int base
) {
    std::vector<u64> digits;

    if (m == 0) {
        digits.push_back(0);
        return digits;
    }

    while (m > 0) {
        digits.push_back(
            m % static_cast<u64>(base)
        );

        m /= static_cast<u64>(base);
    }

    return digits;
}

/* =========================================================
 * MISS predicate
 * ========================================================= */

static bool is_miss(
    u64 x,
    u64 m,
    int base
) {
    while (x > 0 || m > 0) {
        const u64 xd =
            x % static_cast<u64>(base);

        const u64 md =
            m % static_cast<u64>(base);

        if (xd > md) {
            return false;
        }

        x /= static_cast<u64>(base);
        m /= static_cast<u64>(base);
    }

    return true;
}

/* =========================================================
 * Direct HIT intervals
 *
 * Only used for small cases.
 * ========================================================= */

static std::vector<Interval>
brute_intervals(
    u64 m,
    int base
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 x = 0; x <= m; ++x) {
        const bool hit =
            !is_miss(
                x,
                m,
                base
            );

        if (hit && !inside) {
            inside = true;
            start = x;
        }

        if (!hit && inside) {
            result.push_back(
                Interval{
                    start,
                    x - 1
                }
            );

            inside = false;
        }
    }

    if (inside) {
        result.push_back(
            Interval{
                start,
                m
            }
        );
    }

    return result;
}

/* =========================================================
 * Convert intervals to length histogram
 *
 * histogram[length] = number of intervals of that length.
 * ========================================================= */

static std::map<u64, u64>
histogram_from_intervals(
    const std::vector<Interval>& intervals
) {
    std::map<u64, u64> result;

    for (const auto& in : intervals) {
        const u64 length =
            in.right - in.left + 1;

        ++result[length];
    }

    return result;
}

/* =========================================================
 * MISS count
 * ========================================================= */

static u64 miss_count(
    u64 m,
    int base
) {
    u64 result = 1;

    const auto digits =
        digits_of(m, base);

    for (u64 d : digits) {
        result *= d + 1;
    }

    return result;
}

/* =========================================================
 * HIT count
 * ========================================================= */

static u64 direct_hit_count(
    u64 m,
    int base
) {
    return m + 1 -
           miss_count(
               m,
               base
           );
}

/* =========================================================
 * Highest power b^k <= m
 * ========================================================= */

static u64 highest_power(
    u64 m,
    int base
) {
    u64 power = 1;

    while (
        power <=
        m / static_cast<u64>(base)
    ) {
        power *=
            static_cast<u64>(base);
    }

    return power;
}

/* =========================================================
 * Recursive interval-length histogram
 *
 * If
 *
 *     m = d*b^k + r
 *
 * then
 *
 *     A_m(l)
 *       = (d+1) A_r(l)
 *         + d * [l = b^k-r-1]
 *
 * for positive bridge length.
 * ========================================================= */

static std::map<u64, u64>
recursive_histogram(
    u64 m,
    int base
) {
    /*
     * A one-digit m has no HIT values.
     */
    if (m < static_cast<u64>(base)) {
        return {};
    }

    const u64 power =
        highest_power(
            m,
            base
        );

    const u64 d =
        m / power;

    const u64 r =
        m % power;

    const auto lower =
        recursive_histogram(
            r,
            base
        );

    std::map<u64, u64> result;

    /*
     * Replicate the lower spectrum in all d+1 blocks.
     */
    for (const auto& entry : lower) {
        result[entry.first] =
            entry.second * (d + 1);
    }

    /*
     * Add d bridge intervals if the bridge has
     * positive length.
     */
    const u64 bridge_length =
        power - r - 1;

    if (
        d > 0 &&
        bridge_length > 0
    ) {
        result[bridge_length] += d;
    }

    return result;
}

/* =========================================================
 * Histogram equality
 * ========================================================= */

static bool same_histogram(
    const std::map<u64, u64>& a,
    const std::map<u64, u64>& b
) {
    return a == b;
}

/* =========================================================
 * Sum number of intervals
 *
 * sum_l A(l)
 * ========================================================= */

static u64 histogram_interval_count(
    const std::map<u64, u64>& histogram
) {
    u64 total = 0;

    for (const auto& entry : histogram) {
        total += entry.second;
    }

    return total;
}

/* =========================================================
 * Weighted sum
 *
 * sum_l l*A(l)
 *
 * This must equal the total number of HIT integers.
 * ========================================================= */

static u64 histogram_weighted_count(
    const std::map<u64, u64>& histogram
) {
    u64 total = 0;

    for (const auto& entry : histogram) {
        total +=
            entry.first *
            entry.second;
    }

    return total;
}

/* =========================================================
 * Maximum interval length in histogram
 * ========================================================= */

static u64 histogram_max_length(
    const std::map<u64, u64>& histogram
) {
    if (histogram.empty()) {
        return 0;
    }

    return histogram.rbegin()->first;
}

/* =========================================================
 * Print a histogram
 *
 * Only used for failures.
 * ========================================================= */

static void print_histogram(
    const std::map<u64, u64>& histogram
) {
    std::cout << "{";

    bool first = true;

    for (const auto& entry : histogram) {
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
 * Verify one case
 * ========================================================= */

static bool verify_case(
    u64 m,
    int base,
    bool brute_check,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    /*
     * --------------------------------------------------
     * Recursive histogram
     * --------------------------------------------------
     */

    const auto recursive =
        recursive_histogram(
            m,
            base
        );

    /*
     * --------------------------------------------------
     * Direct HIT count
     * --------------------------------------------------
     */

    const u64 expected_hits =
        direct_hit_count(
            m,
            base
        );

    /*
     * --------------------------------------------------
     * Histogram count
     * --------------------------------------------------
     */

    const u64 interval_count =
        histogram_interval_count(
            recursive
        );

    if (interval_count == 0 && expected_hits != 0) {
        ok = false;
    }

    /*
     * Weighted histogram must equal total HIT count.
     */
    const u64 weighted =
        histogram_weighted_count(
            recursive
        );

    if (weighted != expected_hits) {
        ok = false;

        ++stats.weighted_failures;

        if (print_failure) {
            std::cout
                << "FAIL_WEIGHTED "
                << "m=" << m
                << " base=" << base
                << " weighted="
                << weighted
                << " expected="
                << expected_hits
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Independent brute-force histogram
     * --------------------------------------------------
     */

    if (brute_check) {
        const auto brute_intervals_result =
            brute_intervals(
                m,
                base
            );

        const auto brute =
            histogram_from_intervals(
                brute_intervals_result
            );

        /*
         * Exact histogram equality.
         */
        if (
            !same_histogram(
                recursive,
                brute
            )
        ) {
            ok = false;

            ++stats.histogram_failures;
            ++stats.brute_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_HISTOGRAM "
                    << "m=" << m
                    << " base=" << base
                    << "\n";

                std::cout
                    << "recursive=";

                print_histogram(
                    recursive
                );

                std::cout
                    << "\nbrute=";

                print_histogram(
                    brute
                );

                std::cout
                    << "\n";
            }
        }

        /*
         * Direct interval count.
         */
        if (
            interval_count !=
            brute_intervals_result.size()
        ) {
            ok = false;

            ++stats.count_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_INTERVAL_COUNT "
                    << "m=" << m
                    << " base=" << base
                    << " recursive="
                    << interval_count
                    << " brute="
                    << brute_intervals_result.size()
                    << "\n";
            }
        }

        /*
         * Independent direct sum of all lengths.
         */
        u64 brute_weighted = 0;

        for (const auto& in :
             brute_intervals_result) {
            brute_weighted +=
                in.right -
                in.left +
                1;
        }

        if (brute_weighted != weighted) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_BRUTE_WEIGHT "
                    << "m=" << m
                    << " base=" << base
                    << "\n";
            }
        }
    }

    /*
     * --------------------------------------------------
     * Basic internal checks
     * --------------------------------------------------
     */

    /*
     * Every recorded length must be positive.
     */
    for (const auto& entry : recursive) {
        if (
            entry.first == 0 ||
            entry.second == 0
        ) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_ZERO_BUCKET "
                    << "m=" << m
                    << " base=" << base
                    << "\n";
            }

            break;
        }
    }

    /*
     * Maximum possible interval length cannot exceed m.
     */
    if (
        histogram_max_length(
            recursive
        ) > m
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_MAX_LENGTH "
                << "m=" << m
                << " base=" << base
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
        << "START EXPERIMENT 328\n"
        << "ARBITRARY-m HIT INTERVAL LENGTH SPECTRUM\n"
        << "DOES THE RECURSIVE HISTOGRAM FORMULA HOLD?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small cases.
     *
     * Exact comparison of the complete histogram.
     * --------------------------------------------------
     */

    for (int base : BASES) {
        for (u64 m = 0; m <= 750; ++m) {
            ++stats.tested;

            verify_case(
                m,
                base,
                true,
                stats.failed < 10,
                stats
            );
        }
    }

    /*
     * --------------------------------------------------
     * Phase 2:
     * Random large arbitrary m.
     *
     * No integer scanning and no interval generation.
     * Only the recursive histogram is constructed.
     *
     * The histogram contains at most O(number of digits)
     * distinct lengths.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x328328328ULL
    );

    constexpr u64 RANDOM_CASES = 10000;

    for (
        u64 i = 0;
        i < RANDOM_CASES;
        ++i
    ) {
        const int base =
            BASES[
                rng() % BASES.size()
            ];

        const u64 m =
            rng() %
            1000000000000000000ULL;

        ++stats.tested;

        verify_case(
            m,
            base,
            false,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 3:
     * Targeted digit patterns.
     * --------------------------------------------------
     */

    u64 targeted = 0;

    for (int base : BASES) {
        for (int length = 1;
             length <= 12;
             ++length) {

            std::vector<
                std::vector<u64>
            > patterns;

            /*
             * All zero.
             */
            patterns.emplace_back(
                length,
                0
            );

            /*
             * All maximum.
             */
            patterns.emplace_back(
                length,
                static_cast<u64>(
                    base - 1
                )
            );

            /*
             * Alternating 0 / max.
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
                                base - 1
                            );
                    }
                }

                patterns.push_back(d);
            }

            /*
             * Alternating max / 0.
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
                                base - 1
                            );
                    }
                }

                patterns.push_back(d);
            }

            /*
             * All ones.
             *
             * Skip if base=2? It is still valid.
             */
            patterns.emplace_back(
                length,
                1
            );

            /*
             * Increasing digits.
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
                            i % base
                        );
                }

                patterns.push_back(d);
            }

            /*
             * Low digits maximum,
             * high digits one.
             */
            {
                std::vector<u64> d(
                    length,
                    1
                );

                for (int i = 0;
                     i < length / 2;
                     ++i) {
                    d[i] =
                        static_cast<u64>(
                            base - 1
                        );
                }

                patterns.push_back(d);
            }

            /*
             * High digits maximum,
             * low digits one.
             */
            {
                std::vector<u64> d(
                    length,
                    1
                );

                for (int i = length / 2;
                     i < length;
                     ++i) {
                    d[i] =
                        static_cast<u64>(
                            base - 1
                        );
                }

                patterns.push_back(d);
            }

            for (const auto& digits :
                 patterns) {

                u64 m = 0;
                u64 power = 1;
                bool overflow = false;

                for (u64 d : digits) {
                    if (
                        d != 0 &&
                        power >
                            UINT64_MAX / d
                    ) {
                        overflow = true;
                        break;
                    }

                    const u64 term =
                        d * power;

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
                            static_cast<u64>(base)
                    ) {
                        overflow = true;
                        break;
                    }

                    power *=
                        static_cast<u64>(base);
                }

                if (overflow) {
                    continue;
                }

                ++targeted;
                ++stats.tested;

                verify_case(
                    m,
                    base,
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
        << "count_failures="
        << stats.count_failures
        << "\n"
        << "weighted_failures="
        << stats.weighted_failures
        << "\n"
        << "brute_failures="
        << stats.brute_failures
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 328\n";

    return 0;
}
