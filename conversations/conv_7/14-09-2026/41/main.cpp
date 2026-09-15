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

constexpr int MAX_MOMENT = 3;

struct Interval {
    u64 left;
    u64 right;
};

struct Moments {
    u128 value[MAX_MOMENT + 1] = {};
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 closed_failures = 0;
    u64 recurrence_failures = 0;
    u64 brute_failures = 0;
};

/* =========================================================
 * Print uint128
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

static std::vector<u64>
digits_of(
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

        result.value[0] += 1;
        result.value[1] += L;
        result.value[2] += L * L;
        result.value[3] += L * L * L;
    }

    return result;
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
 * Recursive moment hierarchy
 *
 * m = d*b^k + r
 *
 * S_j(m)
 *   = (d+1) S_j(r) + d G^j
 *
 * for j >= 1,
 * with S_0 counting intervals.
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

    for (int j = 0;
         j <= MAX_MOMENT;
         ++j) {
        result.value[j] =
            static_cast<u128>(d + 1) *
            lower.value[j];
    }

    const u64 G =
        power - r - 1;

    if (
        d > 0 &&
        G > 0
    ) {
        result.value[0] +=
            static_cast<u128>(d);

        const u128 Gu =
            static_cast<u128>(G);

        u128 power_G = 1;

        for (int j = 1;
             j <= MAX_MOMENT;
             ++j) {
            power_G *= Gu;

            result.value[j] +=
                static_cast<u128>(d) *
                power_G;
        }
    }

    return result;
}

/* =========================================================
 * Closed digit-sum formula
 *
 * digits[k] = m_k
 *
 * R_k = sum_{i<k} m_i b^i
 *
 * G_k = b^k - R_k - 1
 *
 * C_k = m_k * product_{i>k}(m_i+1)
 *
 * S_j = sum_k C_k G_k^j
 *
 * for j >= 1.
 *
 * For j=0, only G_k > 0 contributes.
 * ========================================================= */

static Moments
closed_digit_moments(
    u64 m,
    int base
) {
    const auto digits =
        digits_of(
            m,
            base
        );

    const std::size_t L =
        digits.size();

    Moments result;

    /*
     * Suffix products:
     *
     * suffix[k] =
     * product_{i>k}(m_i+1)
     *
     * Construct from high to low.
     */
    std::vector<u128> suffix(
        L,
        1
    );

    u128 product = 1;

    for (
        std::size_t k = L;
        k > 0;
        --k
    ) {
        const std::size_t i =
            k - 1;

        suffix[i] = product;

        product *=
            static_cast<u128>(
                digits[i] + 1
            );
    }

    /*
     * lower = R_k
     *
     * power = b^k
     *
     * We walk k from low to high.
     */
    u128 lower = 0;
    u128 power = 1;

    for (
        std::size_t k = 0;
        k < L;
        ++k
    ) {
        const u64 digit =
            digits[k];

        const u128 G128 =
            power >
                    static_cast<u128>(
                        lower
                    ) +
                        static_cast<u128>(1)
                ? power -
                      lower -
                      static_cast<u128>(1)
                : 0;

        const bool positive_gap =
            (G128 > 0);

        const u128 coefficient =
            static_cast<u128>(digit) *
            suffix[k];

        /*
         * j = 0
         */
        if (
            positive_gap &&
            coefficient > 0
        ) {
            result.value[0] +=
                coefficient;
        }

        /*
         * j = 1,2,3
         */
        if (
            positive_gap &&
            coefficient > 0
        ) {
            u128 Gpow = 1;

            for (
                int j = 1;
                j <= MAX_MOMENT;
                ++j
            ) {
                Gpow *= G128;

                result.value[j] +=
                    coefficient *
                    Gpow;
            }
        }

        /*
         * Advance to R_{k+1}.
         */
        lower +=
            static_cast<u128>(digit) *
            power;

        power *=
            static_cast<u128>(base);
    }

    return result;
}

/* =========================================================
 * Compare moments
 * ========================================================= */

static bool same_moments(
    const Moments& a,
    const Moments& b
) {
    for (int j = 0;
         j <= MAX_MOMENT;
         ++j) {
        if (
            a.value[j] !=
            b.value[j]
        ) {
            return false;
        }
    }

    return true;
}

/* =========================================================
 * Print first mismatch
 * ========================================================= */

static void print_moment_mismatch(
    const char* name,
    const Moments& a,
    const Moments& b
) {
    for (int j = 0;
         j <= MAX_MOMENT;
         ++j) {
        if (
            a.value[j] !=
            b.value[j]
        ) {
            std::cout
                << "  "
                << name
                << " S"
                << j
                << ": actual=";

            print_u128(
                a.value[j]
            );

            std::cout
                << " expected=";

            print_u128(
                b.value[j]
            );

            std::cout
                << "\n";
        }
    }
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

    const Moments closed =
        closed_digit_moments(
            m,
            base
        );

    /*
     * --------------------------------------------------
     * Closed digit formula vs recursive formula
     * --------------------------------------------------
     */

    if (
        !same_moments(
            recursive,
            closed
        )
    ) {
        ok = false;

        ++stats.closed_failures;

        if (print_failure) {
            std::cout
                << "FAIL_CLOSED_FORM "
                << "m=" << m
                << " base=" << base
                << "\n";

            print_moment_mismatch(
                "recursive-vs-closed",
                recursive,
                closed
            );
        }
    }

    /*
     * --------------------------------------------------
     * Exact brute-force comparison
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

        if (
            !same_moments(
                recursive,
                brute
            )
        ) {
            ok = false;

            ++stats.brute_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_BRUTE "
                    << "m=" << m
                    << " base=" << base
                    << "\n";

                print_moment_mismatch(
                    "recursive-vs-brute",
                    recursive,
                    brute
                );
            }
        }
    }

    /*
     * --------------------------------------------------
     * Basic consistency:
     *
     * S1 must equal total HIT count.
     * --------------------------------------------------
     */

    u64 miss = 1;

    const auto digits =
        digits_of(
            m,
            base
        );

    for (u64 d : digits) {
        miss *= d + 1;
    }

    const u128 expected_s1 =
        static_cast<u128>(
            m + 1 - miss
        );

    if (
        recursive.value[1] !=
        expected_s1
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_HIT_COUNT "
                << "m=" << m
                << " base=" << base
                << " actual=";

            print_u128(
                recursive.value[1]
            );

            std::cout
                << " expected=";

            print_u128(
                expected_s1
            );

            std::cout
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
        << "START EXPERIMENT 330\n"
        << "CLOSED DIGIT EXPANSION OF HIT INTERVAL MOMENTS\n"
        << "DO THE MOMENTS REDUCE TO A DIRECT SUM OVER DIGIT LEVELS?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small cases.
     *
     * Full S0..S3 brute-force comparison.
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
     * The closed formula uses only O(number of digits).
     * No HIT/MISS enumeration occurs.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x330330330ULL
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
                                base - 1
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
                                base - 1
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

            for (const auto& ds :
                 patterns) {

                u64 m = 0;
                u64 power = 1;
                bool overflow = false;

                for (u64 d : ds) {
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
        << "closed_failures="
        << stats.closed_failures
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
        << "FINISHED EXPERIMENT 330\n";

    return 0;
}
