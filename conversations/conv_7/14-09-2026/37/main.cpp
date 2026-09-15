#include <array>
#include <cstdint>
#include <iostream>
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
    u64 brute_failures = 0;
    u64 count_failures = 0;
    u64 interval_count_failures = 0;
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
 * HIT count directly from MISS count
 * ========================================================= */

static u64 hit_count(
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
 * Direct brute-force HIT intervals
 *
 * Used only for small m.
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
 * Highest power b^k <= m
 *
 * Returns b^k.
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
 * Recursive interval COUNT
 *
 * If
 *
 *     m = d*b^k + r
 *
 * then
 *
 *     R(m) = (d+1)R(r) + d
 *
 * because:
 *
 *   - each of the d+1 blocks contains R(r)
 *     recursive intervals;
 *   - each of the d boundaries contributes one
 *     bridge interval whenever the bridge is nonempty.
 *
 * We calculate this without constructing intervals.
 * ========================================================= */

static u64 recursive_interval_count(
    u64 m,
    int base
) {
    if (m < static_cast<u64>(base)) {
        return 0;
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

    const u64 lower =
        recursive_interval_count(
            r,
            base
        );

    const u64 bridge_count =
        (r + 1 < power)
            ? d
            : 0;

    return
        (d + 1) * lower +
        bridge_count;
}

/* =========================================================
 * Recursive HIT cardinality
 *
 * H(m) =
 *
 *   (d+1) H(r)
 *   +
 *   d*(b^k-r-1)
 *
 * The second term counts the bridge intervals.
 * ========================================================= */

static u64 recursive_hit_count(
    u64 m,
    int base
) {
    if (m < static_cast<u64>(base)) {
        return 0;
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

    const u64 lower =
        recursive_hit_count(
            r,
            base
        );

    const u64 bridge_length =
        power - r - 1;

    return
        (d + 1) * lower +
        d * bridge_length;
}

/* =========================================================
 * Compare interval vectors
 * ========================================================= */

static bool same_intervals(
    const std::vector<Interval>& a,
    const std::vector<Interval>& b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (
        std::size_t i = 0;
        i < a.size();
        ++i
    ) {
        if (
            a[i].left != b[i].left ||
            a[i].right != b[i].right
        ) {
            return false;
        }
    }

    return true;
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
     * Direct cardinality formula
     * --------------------------------------------------
     */

    const u64 direct_hits =
        hit_count(
            m,
            base
        );

    /*
     * --------------------------------------------------
     * Recursive cardinality
     * --------------------------------------------------
     */

    const u64 recursive_hits =
        recursive_hit_count(
            m,
            base
        );

    if (
        direct_hits !=
        recursive_hits
    ) {
        ok = false;
        ++stats.count_failures;

        if (print_failure) {
            std::cout
                << "FAIL_HIT_COUNT "
                << "m=" << m
                << " base=" << base
                << " direct="
                << direct_hits
                << " recursive="
                << recursive_hits
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Recursive interval count
     * --------------------------------------------------
     */

    const u64 direct_miss =
        miss_count(
            m,
            base
        );

    /*
     * The direct interval count can be obtained from
     * the MISS ordering:
     *
     * every positive gap between consecutive MISS values
     * produces exactly one HIT interval.
     *
     * For small m we obtain this directly.
     */

    if (brute_check) {
        const auto brute =
            brute_intervals(
                m,
                base
            );

        const u64 recursive_count =
            recursive_interval_count(
                m,
                base
            );

        if (
            recursive_count !=
            brute.size()
        ) {
            ok = false;
            ++stats.interval_count_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_INTERVAL_COUNT "
                    << "m=" << m
                    << " base=" << base
                    << " recursive="
                    << recursive_count
                    << " brute="
                    << brute.size()
                    << "\n";
            }
        }

        /*
         * Full independent interval comparison.
         *
         * This is the expensive part, but only used
         * for small m.
         *
         * Construct the recursive grammar explicitly.
         */

        std::vector<Interval> recursive;

        /*
         * Small recursive constructor.
         */
        struct Builder {
            static void build(
                u64 m,
                int base,
                std::vector<Interval>& out
            ) {
                if (
                    m <
                    static_cast<u64>(base)
                ) {
                    return;
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

                std::vector<Interval> lower;

                build(
                    r,
                    base,
                    lower
                );

                for (
                    u64 q = 0;
                    q <= d;
                    ++q
                ) {
                    const u64 offset =
                        q * power;

                    for (
                        const auto& in :
                        lower
                    ) {
                        out.push_back(
                            Interval{
                                offset + in.left,
                                offset + in.right
                            }
                        );
                    }

                    /*
                     * Bridge between q and q+1.
                     */
                    if (
                        q < d &&
                        r + 1 < power
                    ) {
                        out.push_back(
                            Interval{
                                offset + r + 1,
                                offset + power - 1
                            }
                        );
                    }
                }
            }
        };

        Builder::build(
            m,
            base,
            recursive
        );

        if (
            !same_intervals(
                recursive,
                brute
            )
        ) {
            ok = false;
            ++stats.brute_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_BRUTE_INTERVALS "
                    << "m=" << m
                    << " base=" << base
                    << " recursive="
                    << recursive.size()
                    << " brute="
                    << brute.size()
                    << "\n";
            }
        }

        /*
         * Also verify no unexpected overlap.
         */
        for (
            std::size_t i = 0;
            i + 1 < recursive.size();
            ++i
        ) {
            if (
                recursive[i].right >=
                recursive[i + 1].left
            ) {
                ok = false;

                if (print_failure) {
                    std::cout
                        << "FAIL_OVERLAP "
                        << "m=" << m
                        << " base=" << base
                        << " index=" << i
                        << "\n";
                }

                break;
            }
        }

        /*
         * Sanity:
         * total interval length must equal HIT count.
         */
        u64 total_length = 0;

        for (const auto& in :
             recursive) {
            total_length +=
                in.right -
                in.left +
                1;
        }

        if (
            total_length !=
            direct_hits
        ) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_LENGTH "
                    << "m=" << m
                    << " base=" << base
                    << "\n";
            }
        }
    }

    /*
     * The variable is deliberately evaluated here so
     * the MISS-count side of the theorem is exercised
     * for every test.
     */
    (void)direct_miss;

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
        << "START EXPERIMENT 327\n"
        << "RECURSIVE ARBITRARY-m INTERVAL GRAMMAR\n"
        << "DO DIGIT-BLOCK RECURSION REPRODUCE THE HIT INTERVAL SET?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small cases.
     *
     * Full interval-by-interval comparison.
     * --------------------------------------------------
     */

    for (int base : BASES) {
        for (u64 m = 0; m <= 500; ++m) {
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
     * Random large cases.
     *
     * Only scalar recursive identities are evaluated.
     *
     * No interval vectors are created.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x327327327ULL
    );

    constexpr u64 RANDOM_CASES = 5000;

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
     * Hand-selected digit structures.
     *
     * These stay scalar, so large values are safe.
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
                    base - 1
                )
            );

            /*
             * 0,max,0,max,...
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
             * max,0,max,0,...
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
             * all 1
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
                            i % base
                        );
                }

                patterns.push_back(d);
            }

            for (const auto& d :
                 patterns) {

                u64 m = 0;
                u64 power = 1;
                bool overflow = false;

                for (u64 x : d) {
                    if (
                        x != 0 &&
                        power >
                            UINT64_MAX / x
                    ) {
                        overflow = true;
                        break;
                    }

                    const u64 term =
                        x * power;

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
        << "brute_failures="
        << stats.brute_failures
        << "\n"
        << "count_failures="
        << stats.count_failures
        << "\n"
        << "interval_count_failures="
        << stats.interval_count_failures
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 327\n";

    return 0;
}