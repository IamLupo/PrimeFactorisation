#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

struct Interval {
    u64 left;
    u64 right;
};

struct Moments {
    u128 s0 = 0;
    u128 s1 = 0;
    u128 s2 = 0;
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 s0_failures = 0;
    u64 s1_failures = 0;
    u64 s2_failures = 0;
    u64 brute_failures = 0;
    u64 recurrence_failures = 0;
};

/* =========================================================
 * Printing unsigned __int128
 * ========================================================= */

static void print_u128(u128 value) {
    if (value == 0) {
        std::cout << '0';
        return;
    }

    char buffer[64];
    int pos = 0;

    while (value > 0) {
        const unsigned digit =
            static_cast<unsigned>(value % 10);

        buffer[pos++] =
            static_cast<char>('0' + digit);

        value /= 10;
    }

    while (pos > 0) {
        std::cout << buffer[--pos];
    }
}

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
 * Moments from explicit intervals
 *
 * S0 = number of intervals
 * S1 = sum(length)
 * S2 = sum(length^2)
 * ========================================================= */

static Moments
moments_from_intervals(
    const std::vector<Interval>& intervals
) {
    Moments result;

    for (const auto& in : intervals) {
        const u64 length =
            in.right - in.left + 1;

        const u128 L =
            static_cast<u128>(length);

        ++result.s0;
        result.s1 += L;
        result.s2 += L * L;
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
 * Direct HIT moments from the known HIT count
 *
 * Only S1 can be obtained this way.
 * S0/S2 are checked through the recursive spectrum.
 * ========================================================= */

static u128 direct_s1(
    u64 m,
    int base
) {
    const u64 misses =
        miss_count(
            m,
            base
        );

    return static_cast<u128>(m + 1 - misses);
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
 * Recursive interval moments
 *
 * m = d*b^k + r
 *
 * There are d+1 copies of the lower interval spectrum
 * and d bridge intervals of length
 *
 *     L = b^k-r-1
 *
 * when L > 0.
 *
 * Therefore:
 *
 *     S0(m) = (d+1)S0(r) + d
 *
 * when L>0, and
 *
 *     Sj(m) = (d+1)Sj(r) + d*L^j
 *
 * for j=1,2.
 * ========================================================= */

static Moments
recursive_moments(
    u64 m,
    int base
) {
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

    const Moments lower =
        recursive_moments(
            r,
            base
        );

    Moments result;

    result.s0 =
        static_cast<u128>(d + 1) *
        lower.s0;

    result.s1 =
        static_cast<u128>(d + 1) *
        lower.s1;

    result.s2 =
        static_cast<u128>(d + 1) *
        lower.s2;

    const u64 bridge_length =
        power - r - 1;

    if (
        d > 0 &&
        bridge_length > 0
    ) {
        const u128 L =
            static_cast<u128>(
                bridge_length
            );

        result.s0 +=
            static_cast<u128>(d);

        result.s1 +=
            static_cast<u128>(d) * L;

        result.s2 +=
            static_cast<u128>(d) * L * L;
    }

    return result;
}

/* =========================================================
 * Independently calculate recursive moments by repeatedly
 * peeling the highest digit.
 *
 * This is deliberately iterative, so we have a second
 * implementation of the recurrence.
 * ========================================================= */

static Moments
iterative_moments(
    u64 m,
    int base
) {
    Moments result;

    u64 current = m;

    while (
        current >=
        static_cast<u64>(base)
    ) {
        const u64 power =
            highest_power(
                current,
                base
            );

        const u64 d =
            current / power;

        const u64 r =
            current % power;

        /*
         * We currently have the contribution from all
         * lower-level structures. Starting from zero,
         * scale it when moving upward.
         */
        result.s0 =
            static_cast<u128>(d + 1) *
            result.s0;

        result.s1 =
            static_cast<u128>(d + 1) *
            result.s1;

        result.s2 =
            static_cast<u128>(d + 1) *
            result.s2;

        const u64 bridge_length =
            power - r - 1;

        if (
            d > 0 &&
            bridge_length > 0
        ) {
            const u128 L =
                static_cast<u128>(
                    bridge_length
                );

            result.s0 +=
                static_cast<u128>(d);

            result.s1 +=
                static_cast<u128>(d) * L;

            result.s2 +=
                static_cast<u128>(d) * L * L;
        }

        current = r;
    }

    return result;
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

    const Moments recursive =
        recursive_moments(
            m,
            base
        );

    const Moments iterative =
        iterative_moments(
            m,
            base
        );

    /*
     * --------------------------------------------------
     * Internal recurrence implementation agreement
     * --------------------------------------------------
     */

    if (
        recursive.s0 != iterative.s0 ||
        recursive.s1 != iterative.s1 ||
        recursive.s2 != iterative.s2
    ) {
        ok = false;
        ++stats.recurrence_failures;

        if (print_failure) {
            std::cout
                << "FAIL_IMPLEMENTATIONS "
                << "m=" << m
                << " base=" << base
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * S1 against independent MISS cardinality
     * --------------------------------------------------
     */

    const u128 expected_s1 =
        direct_s1(
            m,
            base
        );

    if (recursive.s1 != expected_s1) {
        ok = false;
        ++stats.s1_failures;

        if (print_failure) {
            std::cout
                << "FAIL_S1 "
                << "m=" << m
                << " base=" << base
                << " actual=";

            print_u128(recursive.s1);

            std::cout
                << " expected=";

            print_u128(expected_s1);

            std::cout
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Exact brute-force moment comparison for small m
     * --------------------------------------------------
     */

    if (brute_check) {
        const auto intervals =
            brute_intervals(
                m,
                base
            );

        const Moments brute =
            moments_from_intervals(
                intervals
            );

        if (recursive.s0 != brute.s0) {
            ok = false;
            ++stats.s0_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_S0 "
                    << "m=" << m
                    << " base=" << base
                    << " actual="
                    << static_cast<u64>(recursive.s0)
                    << " expected="
                    << static_cast<u64>(brute.s0)
                    << "\n";
            }
        }

        if (recursive.s1 != brute.s1) {
            ok = false;
            ++stats.s1_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_BRUTE_S1 "
                    << "m=" << m
                    << " base=" << base
                    << "\n";
            }
        }

        if (recursive.s2 != brute.s2) {
            ok = false;
            ++stats.s2_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_S2 "
                    << "m=" << m
                    << " base=" << base
                    << " actual=";

                print_u128(recursive.s2);

                std::cout
                    << " expected=";

                print_u128(brute.s2);

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
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 329\n"
        << "ARBITRARY-m INTERVAL MOMENTS\n"
        << "DO THE FIRST AND SECOND MOMENTS OBEY THE RECURSIVE GRAMMAR?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small cases.
     *
     * Full interval construction allows exact S0,S1,S2
     * comparison.
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
     * Random arbitrary m.
     *
     * m <= 10^12 keeps S2 safely inside uint128:
     *
     *     S2 <= m^3 < 10^36.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x329329329ULL
    );

    constexpr u64 RANDOM_CASES = 15000;

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
            rng() % 1000000000000ULL;

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
             * Alternating 0 / maximum.
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
             * Alternating maximum / 0.
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
             */
            patterns.emplace_back(
                length,
                1
            );

            /*
             * Increasing modulo-base digits.
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

                for (u64 digit : d) {
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
        << "s0_failures="
        << stats.s0_failures
        << "\n"
        << "s1_failures="
        << stats.s1_failures
        << "\n"
        << "s2_failures="
        << stats.s2_failures
        << "\n"
        << "recurrence_failures="
        << stats.recurrence_failures
        << "\n"
        << "brute_failures="
        << stats.brute_failures
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 329\n";

    return 0;
}
