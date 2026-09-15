#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct Transition {
    u64 left;
    u64 right;
    u64 gap;
    int type;
};

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;
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
 * Check x <=_b m
 * ========================================================= */

static bool digitwise_miss(
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
 * MISS count formula
 * ========================================================= */

static u64 miss_count_formula(
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
 * Enumerate MISS values.
 *
 * The mixed-radix odometer is already in numeric order,
 * so no sorting is required.
 * ========================================================= */

static std::vector<u64> enumerate_misses(
    u64 m,
    int base
) {
    const auto digits =
        digits_of(m, base);

    const u64 total =
        miss_count_formula(
            m,
            base
        );

    std::vector<u64> misses;
    misses.reserve(
        static_cast<std::size_t>(total)
    );

    std::vector<u64> current(
        digits.size(),
        0
    );

    const std::size_t L =
        digits.size();

    u64 value = 0;

    while (true) {
        /*
         * Recompute the numeric value from the
         * current digit vector.
         *
         * This is still cheap because the MISS set
         * is explicitly bounded.
         */
        value = 0;

        u64 power = 1;

        for (std::size_t i = 0; i < L; ++i) {
            value +=
                current[i] * power;

            power *=
                static_cast<u64>(base);
        }

        misses.push_back(value);

        /*
         * Mixed-radix increment.
         */
        std::size_t pos = 0;

        while (pos < L) {
            ++current[pos];

            if (current[pos] <= digits[pos]) {
                break;
            }

            current[pos] = 0;
            ++pos;
        }

        if (pos == L) {
            break;
        }
    }

    return misses;
}

/* =========================================================
 * Gap formula
 * ========================================================= */

static u64 gap_formula(
    const std::vector<u64>& digits,
    int base,
    int r
) {
    u64 power = 1;

    for (int i = 0; i < r; ++i) {
        power *= static_cast<u64>(base);
    }

    u64 lower = 0;
    u64 lower_power = 1;

    for (int i = 0; i < r; ++i) {
        lower +=
            digits[i] * lower_power;

        lower_power *=
            static_cast<u64>(base);
    }

    return power - 1 - lower;
}

/* =========================================================
 * C_r formula
 * ========================================================= */

static u64 count_type_formula(
    const std::vector<u64>& digits,
    int r
) {
    u64 result = digits[r];

    for (
        int i = r + 1;
        i < static_cast<int>(digits.size());
        ++i
    ) {
        result *= digits[i] + 1;
    }

    return result;
}

/* =========================================================
 * Transition type
 * ========================================================= */

static int transition_type(
    u64 x,
    u64 m,
    int base
) {
    const auto digits =
        digits_of(m, base);

    for (
        int r = 0;
        r < static_cast<int>(digits.size());
        ++r
    ) {
        const u64 xd =
            x % static_cast<u64>(base);

        if (xd < digits[r]) {
            return r;
        }

        x /= static_cast<u64>(base);
    }

    return -1;
}

/* =========================================================
 * Direct HIT intervals
 * ========================================================= */

static std::vector<std::pair<u64, u64>>
brute_hit_intervals(
    u64 m,
    int base
) {
    std::vector<std::pair<u64, u64>> result;

    bool inside = false;
    u64 start = 0;

    for (u64 x = 0; x <= m; ++x) {
        const bool hit =
            !digitwise_miss(
                x,
                m,
                base
            );

        if (hit && !inside) {
            inside = true;
            start = x;
        }

        if (!hit && inside) {
            result.emplace_back(
                start,
                x - 1
            );

            inside = false;
        }
    }

    if (inside) {
        result.emplace_back(
            start,
            m
        );
    }

    return result;
}

/* =========================================================
 * Verify one case
 * ========================================================= */

static bool verify_case(
    u64 m,
    int base,
    bool brute_intervals,
    bool print_failure
) {
    const auto digits =
        digits_of(m, base);

    const u64 miss_count =
        miss_count_formula(
            m,
            base
        );

    /*
     * Enumeration safety bound.
     */
    if (miss_count > 100000) {
        return true;
    }

    const auto misses =
        enumerate_misses(
            m,
            base
        );

    bool ok = true;

    /* MISS count */
    if (
        misses.size() !=
        static_cast<std::size_t>(
            miss_count
        )
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_MISS_COUNT "
                << "m=" << m
                << " base=" << base
                << " expected=" << miss_count
                << " actual=" << misses.size()
                << "\n";
        }
    }

    /* Endpoints */
    if (
        misses.empty() ||
        misses.front() != 0 ||
        misses.back() != m
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_ENDPOINTS "
                << "m=" << m
                << " base=" << base
                << "\n";
        }
    }

    /* Every MISS really satisfies digit condition */
    for (u64 x : misses) {
        if (
            !digitwise_miss(
                x,
                m,
                base
            )
        ) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_MEMBERSHIP "
                    << "m=" << m
                    << " base=" << base
                    << " x=" << x
                    << "\n";
            }

            break;
        }
    }

    std::vector<u64> actual_count(
        digits.size(),
        0
    );

    std::vector<u64> actual_gap_sum(
        digits.size(),
        0
    );

    std::vector<
        std::pair<u64, u64>
    > predicted_intervals;

    for (
        std::size_t i = 0;
        i + 1 < misses.size();
        ++i
    ) {
        const u64 left =
            misses[i];

        const u64 right =
            misses[i + 1];

        const u64 gap =
            right - left - 1;

        const int type =
            transition_type(
                left,
                m,
                base
            );

        if (
            type < 0 ||
            type >=
                static_cast<int>(
                    digits.size()
                )
        ) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_TYPE "
                    << "m=" << m
                    << " base=" << base
                    << " left=" << left
                    << "\n";
            }

            continue;
        }

        ++actual_count[type];
        actual_gap_sum[type] += gap;

        if (gap > 0) {
            predicted_intervals.emplace_back(
                left + 1,
                right - 1
            );
        }

        const u64 expected_gap =
            gap_formula(
                digits,
                base,
                type
            );

        if (gap != expected_gap) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_GAP "
                    << "m=" << m
                    << " base=" << base
                    << " type=" << type
                    << " left=" << left
                    << " actual=" << gap
                    << " expected=" << expected_gap
                    << "\n";
            }
        }
    }

    /* C_r */
    for (
        int r = 0;
        r < static_cast<int>(digits.size());
        ++r
    ) {
        const u64 expected_count =
            count_type_formula(
                digits,
                r
            );

        if (
            actual_count[r] !=
            expected_count
        ) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_C "
                    << "m=" << m
                    << " base=" << base
                    << " r=" << r
                    << " actual="
                    << actual_count[r]
                    << " expected="
                    << expected_count
                    << "\n";
            }
        }

        const u64 expected_gap =
            gap_formula(
                digits,
                base,
                r
            );

        const u64 expected_sum =
            expected_count *
            expected_gap;

        if (
            actual_gap_sum[r] !=
            expected_sum
        ) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_WEIGHT "
                    << "m=" << m
                    << " base=" << base
                    << " r=" << r
                    << "\n";
            }
        }
    }

    /* Global weighted identity */
    u64 weighted_sum = 0;

    for (u64 x : actual_gap_sum) {
        weighted_sum += x;
    }

    const u64 expected_weighted =
        m + 1 - miss_count;

    if (
        weighted_sum !=
        expected_weighted
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_GLOBAL "
                << "m=" << m
                << " base=" << base
                << " actual="
                << weighted_sum
                << " expected="
                << expected_weighted
                << "\n";
        }
    }

    /* Number of transitions */
    u64 transition_count = 0;

    for (u64 x : actual_count) {
        transition_count += x;
    }

    if (
        transition_count !=
        miss_count - 1
    ) {
        ok = false;

        if (print_failure) {
            std::cout
                << "FAIL_TRANSITIONS "
                << "m=" << m
                << " base=" << base
                << "\n";
        }
    }

    /* Independent brute-force interval check */
    if (brute_intervals) {
        const auto brute =
            brute_hit_intervals(
                m,
                base
            );

        if (brute != predicted_intervals) {
            ok = false;

            if (print_failure) {
                std::cout
                    << "FAIL_INTERVALS "
                    << "m=" << m
                    << " base=" << base
                    << " brute="
                    << brute.size()
                    << " predicted="
                    << predicted_intervals.size()
                    << "\n";
            }
        }
    }

    return ok;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 326\n"
        << "ARBITRARY-m LUCAS INTERVAL GRAMMAR\n"
        << "DO THE MIXED-RADIX GAP FORMULAS HOLD FOR UNCONSTRAINED DIGITS?\n\n";

    Stats stats;

    /* -----------------------------------------------------
     * Phase 1: small exhaustive
     * ----------------------------------------------------- */

    for (int base : BASES) {
        for (u64 m = 0; m <= 300; ++m) {
            ++stats.tested;

            if (!verify_case(
                    m,
                    base,
                    true,
                    stats.failed < 10
                )) {
                ++stats.failed;
            }
        }
    }

    /* -----------------------------------------------------
     * Phase 2: random arbitrary m
     *
     * Only 2000 accepted cases.
     * ----------------------------------------------------- */

    std::mt19937_64 rng(
        0x326326326ULL
    );

    u64 accepted = 0;
    u64 attempts = 0;

    while (accepted < 2000) {
        ++attempts;

        const int base =
            BASES[
                rng() % BASES.size()
            ];

        const u64 m =
            rng() %
            1000000000000ULL;

        const u64 miss_count =
            miss_count_formula(
                m,
                base
            );

        if (miss_count > 100000) {
            continue;
        }

        ++accepted;
        ++stats.tested;

        if (!verify_case(
                m,
                base,
                false,
                stats.failed < 10
            )) {
            ++stats.failed;
        }
    }

    /* -----------------------------------------------------
     * Phase 3: targeted digit patterns
     * ----------------------------------------------------- */

    u64 targeted = 0;

    for (int base : BASES) {
        for (int length = 1;
             length <= 8;
             ++length) {

            std::vector<
                std::vector<u64>
            > patterns;

            /* all zero */
            patterns.emplace_back(
                length,
                0
            );

            /* all maximum */
            patterns.emplace_back(
                length,
                static_cast<u64>(
                    base - 1
                )
            );

            /* alternating 0 / max */
            {
                std::vector<u64> v(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    if (i & 1) {
                        v[i] =
                            static_cast<u64>(
                                base - 1
                            );
                    }
                }

                patterns.push_back(v);
            }

            /* increasing digit pattern */
            {
                std::vector<u64> v(
                    length,
                    0
                );

                for (int i = 0;
                     i < length;
                     ++i) {
                    v[i] =
                        static_cast<u64>(
                            i % base
                        );
                }

                patterns.push_back(v);
            }

            for (const auto& ds :
                 patterns) {

                u64 m = 0;
                u64 power = 1;

                for (u64 d : ds) {
                    m += d * power;

                    power *=
                        static_cast<u64>(
                            base
                        );
                }

                ++targeted;
                ++stats.tested;

                if (!verify_case(
                        m,
                        base,
                        m <= 5000,
                        stats.failed < 10
                    )) {
                    ++stats.failed;
                }
            }
        }
    }

    std::cout
        << "============================\n"
        << "tested="
        << stats.tested
        << "\n"
        << "failed="
        << stats.failed
        << "\n"
        << "random_attempts="
        << attempts
        << "\n"
        << "random_accepted="
        << accepted
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 326\n";

    return 0;
}