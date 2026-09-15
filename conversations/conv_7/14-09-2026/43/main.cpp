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
    u64 repeated_profile_failures = 0;
    u64 type_failures = 0;
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
 * Kummer valuation via borrows
 *
 * Counts borrows in n-k in base p.
 * ========================================================= */

static u64 kummer_borrow_count(
    u64 n,
    u64 k,
    int p
) {
    if (k > n) {
        return UINT64_MAX;
    }

    const u64 prime =
        static_cast<u64>(p);

    u64 borrow = 0;
    u64 count = 0;

    while (
        n > 0 ||
        k > 0 ||
        borrow != 0
    ) {
        const u64 nd =
            n % prime;

        const u64 kd =
            k % prime;

        if (nd < kd + borrow) {
            borrow = 1;
            ++count;
        } else {
            borrow = 0;
        }

        n /= prime;
        k /= prime;
    }

    return count;
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
 * Highest power p^r
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
 * Lower block m mod p^r
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
 * Determine gap type from preceding MISS value
 *
 * r is the lowest digit position whose digit can be
 * incremented in the mixed-radix MISS successor.
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
 * Enumerate all maximal HIT gaps
 *
 * Only used for manageable m.
 * ========================================================= */

static std::vector<Gap> find_typed_gaps(
    u64 m,
    int p
) {
    std::vector<Gap> result;

    bool have_previous_miss = false;
    u64 previous_miss = 0;

    for (u64 x = 0; x <= m; ++x) {
        if (is_miss(x, m, p)) {
            if (have_previous_miss) {
                if (
                    x >
                    previous_miss + 1
                ) {
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
    }

    return result;
}

/* =========================================================
 * Predicted local profile
 *
 * For type r:
 *
 *   lower = m mod p^r
 *   local_n = p^r + lower
 *
 * and y = 1,...,p^r-lower-1
 *
 * because the actual gap starts immediately after the
 * lower MISS boundary.
 *
 * This is the corrected arbitrary-m form.
 * ========================================================= */

static std::map<u64, u64> predicted_local_profile(
    u64 m,
    int p,
    int r
) {
    std::map<u64, u64> result;

    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 lower =
        lower_block(
            m,
            p,
            r
        );

    const u64 gap_length =
        power - lower - 1;

    if (gap_length == 0) {
        return result;
    }

    const u64 local_n =
        power + lower;

    for (
        u64 y = lower + 1;
        y < power;
        ++y
    ) {
        const u64 value =
            binomial_valuation(
                local_n,
                y,
                p
            );

        ++result[value];
    }

    return result;
}

/* =========================================================
 * Actual profile in a gap
 * ========================================================= */

static std::map<u64, u64> actual_gap_profile(
    const Gap& gap,
    u64 m,
    int p
) {
    std::map<u64, u64> result;

    for (
        u64 t = gap.left;
        t <= gap.right;
        ++t
    ) {
        const u64 value =
            binomial_valuation(
                m,
                t,
                p
            );

        ++result[value];
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
 * Print histogram
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

    /*
     * A repeated (r, lower) pair must produce exactly
     * the same valuation histogram.
     */
    std::map<
        std::pair<int, u64>,
        std::map<u64, u64>
    > reference_profiles;

    for (const Gap& gap : gaps) {
        const int r =
            gap.type;

        if (r < 0) {
            ok = false;
            ++stats.type_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_TYPE "
                    << "m=" << m
                    << " p=" << p
                    << " left="
                    << gap.left
                    << "\n";
            }

            continue;
        }

        /*
         * ------------------------------------------------
         * 1. Correct arbitrary-m gap length
         *
         * G_r = p^r - 1 - (m mod p^r)
         * ------------------------------------------------
         */

        const u64 power =
            prime_power(
                p,
                r
            );

        const u64 lower =
            lower_block(
                m,
                p,
                r
            );

        const u64 expected_length =
            power - lower - 1;

        const u64 actual_length =
            gap.right -
            gap.left +
            1;

        if (
            actual_length !=
            expected_length
        ) {
            ok = false;
            ++stats.type_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_LENGTH "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " lower=" << lower
                    << " actual="
                    << actual_length
                    << " expected="
                    << expected_length
                    << "\n";
            }

            continue;
        }

        /*
         * ------------------------------------------------
         * 2. Pointwise local Kummer identity
         *
         * t in the gap corresponds to
         *
         *   y = t mod p^r
         *
         * and the proposed identity is
         *
         *   v_p(C(m,t))
         *   =
         *   v_p(C(p^r+lower,y)).
         * ------------------------------------------------
         */

        const u64 local_n =
            power + lower;

        for (
            u64 t = gap.left;
            t <= gap.right;
            ++t
        ) {
            const u64 y =
                t % power;

            /*
             * Correct gap coordinates:
             *
             * lower < y < p^r
             */
            if (
                y <= lower ||
                y >= power
            ) {
                ok = false;
                ++stats.point_failures;

                if (print_failure) {
                    std::cout
                        << "FAIL_COORDINATE "
                        << "m=" << m
                        << " p=" << p
                        << " r=" << r
                        << " lower=" << lower
                        << " t=" << t
                        << " y=" << y
                        << "\n";
                }

                break;
            }

            const u64 exact =
                binomial_valuation(
                    m,
                    t,
                    p
                );

            const u64 local =
                binomial_valuation(
                    local_n,
                    y,
                    p
                );

            /*
             * Also independently verify Kummer on the
             * actual (m,t) pair.
             */
            const u64 borrow =
                kummer_borrow_count(
                    m,
                    t,
                    p
                );

            if (
                exact != local ||
                exact != borrow
            ) {
                ok = false;
                ++stats.point_failures;

                if (print_failure) {
                    std::cout
                        << "FAIL_POINT "
                        << "m=" << m
                        << " p=" << p
                        << " r=" << r
                        << " t=" << t
                        << " y=" << y
                        << " exact="
                        << exact
                        << " local="
                        << local
                        << " borrow="
                        << borrow
                        << "\n";
                }

                break;
            }
        }

        /*
         * ------------------------------------------------
         * 3. Complete histogram comparison
         * ------------------------------------------------
         */

        const auto actual =
            actual_gap_profile(
                gap,
                m,
                p
            );

        const auto expected =
            predicted_local_profile(
                m,
                p,
                r
            );

        if (
            !same_histogram(
                actual,
                expected
            )
        ) {
            ok = false;
            ++stats.histogram_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_HISTOGRAM "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " lower=" << lower
                    << "\n";

                std::cout
                    << "actual=";

                print_histogram(
                    actual
                );

                std::cout
                    << "\nexpected=";

                print_histogram(
                    expected
                );

                std::cout
                    << "\n";
            }
        }

        /*
         * ------------------------------------------------
         * 4. Same (r, lower) => same profile
         * ------------------------------------------------
         */

        const auto key =
            std::make_pair(
                r,
                lower
            );

        const auto it =
            reference_profiles.find(
                key
            );

        if (
            it == reference_profiles.end()
        ) {
            reference_profiles.emplace(
                key,
                actual
            );
        } else if (
            !same_histogram(
                it->second,
                actual
            )
        ) {
            ok = false;
            ++stats.repeated_profile_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_REPEATED_PROFILE "
                    << "m=" << m
                    << " p=" << p
                    << " r=" << r
                    << " lower=" << lower
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
        << "START EXPERIMENT 332\n"
        << "LOCAL KUMMER PROFILE OF HIT GAPS\n"
        << "DO GAP VALUATIONS DEPEND ONLY ON THE LOWER DIGIT BLOCK?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small cases.
     *
     * Every HIT gap and every point inside each gap.
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
             * Sparse high half.
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

            /*
             * Sparse low half.
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

                if (
                    overflow ||
                    m > 200000
                ) {
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
     * Random unstructured medium-size m.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x332332332ULL
    );

    constexpr u64 RANDOM_CASES = 1000;

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
            rng() % 50000ULL;

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
        << "repeated_profile_failures="
        << stats.repeated_profile_failures
        << "\n"
        << "type_failures="
        << stats.type_failures
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 332\n";

    return 0;
}