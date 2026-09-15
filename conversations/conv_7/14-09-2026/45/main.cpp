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

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 histogram_failures = 0;
    u64 point_failures = 0;
    u64 total_failures = 0;
    u64 max_failures = 0;

    u64 full_histogram_cases = 0;
    u64 sampled_point_cases = 0;
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
 * Direct Kummer borrow count
 *
 * Independent implementation of the valuation.
 * ========================================================= */

static u64 kummer_borrow_count(
    u64 n,
    u64 k,
    int p
) {
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

        if (
            nd <
            kd + borrow
        ) {
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

    for (int i = 0; i < length; ++i) {
        digits[i] =
            n % prime;

        n /= prime;
    }

    return digits;
}

/* =========================================================
 * Exact local histogram
 *
 * This is used only for modest p^r.
 * ========================================================= */

static std::map<u64, u64>
exact_local_histogram(
    u64 L,
    int p,
    int r
) {
    std::map<u64, u64> histogram;

    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 n =
        power + L;

    for (
        u64 y = L + 1;
        y < power;
        ++y
    ) {
        const u64 value =
            binomial_valuation(
                n,
                y,
                p
            );

        ++histogram[value];
    }

    return histogram;
}

/* =========================================================
 * Local digit DP
 *
 * State:
 *
 *   borrow  = current subtraction borrow
 *   relation:
 *
 *     0 = y == L on processed higher significance
 *     1 = y >  L
 *     2 = y <  L
 *
 * We process low -> high because borrow propagates
 * low -> high. The latest differing digit dominates
 * the numerical relation because it is the higher digit.
 * ========================================================= */

static std::map<u64, u64>
local_digit_dp(
    u64 L,
    int p,
    int r
) {
    std::vector<
        std::vector<
            std::vector<u64>
        >
    > dp(
        2,
        std::vector<
            std::vector<u64>
        >(
            3,
            std::vector<u64>(
                static_cast<std::size_t>(r + 1),
                0
            )
        )
    );

    dp[0][0][0] = 1;

    const auto Ldigits =
        digits_of(
            L,
            p,
            r
        );

    for (int i = 0; i < r; ++i) {
        std::vector<
            std::vector<
                std::vector<u64>
            >
        > next(
            2,
            std::vector<
                std::vector<u64>
            >(
                3,
                std::vector<u64>(
                    static_cast<std::size_t>(r + 1),
                    0
                )
            )
        );

        const u64 Li =
            Ldigits[i];

        for (int borrow = 0;
             borrow <= 1;
             ++borrow) {

            for (int relation = 0;
                 relation < 3;
                 ++relation) {

                for (int value = 0;
                     value <= i;
                     ++value) {

                    const u64 count =
                        dp[borrow][relation][value];

                    if (count == 0) {
                        continue;
                    }

                    for (
                        u64 yi = 0;
                        yi < static_cast<u64>(p);
                        ++yi
                    ) {
                        int new_relation =
                            relation;

                        if (yi > Li) {
                            new_relation = 1;
                        } else if (yi < Li) {
                            new_relation = 2;
                        }

                        int new_borrow = 0;
                        int added = 0;

                        if (
                            Li <
                            yi +
                            static_cast<u64>(borrow)
                        ) {
                            new_borrow = 1;
                            added = 1;
                        }

                        next
                            [new_borrow]
                            [new_relation]
                            [value + added]
                            += count;
                    }
                }
            }
        }

        dp.swap(next);
    }

    /*
     * y < p^r, so the top digit of y is zero while
     * p^r+L has top digit one. The top digit absorbs
     * any outstanding borrow.
     *
     * Only relation y>L belongs to the gap.
     */
    std::map<u64, u64> result;

    for (int borrow = 0;
         borrow <= 1;
         ++borrow) {

        for (int relation = 0;
             relation < 3;
             ++relation) {

            for (int value = 0;
                 value <= r;
                 ++value) {

                if (
                    relation == 1 &&
                    dp[borrow][relation][value] != 0
                ) {
                    result[
                        static_cast<u64>(value)
                    ] +=
                        dp[borrow][relation][value];
                }
            }
        }
    }

    return result;
}

/* =========================================================
 * Histogram total
 * ========================================================= */

static u64 histogram_total(
    const std::map<u64, u64>& histogram
) {
    u64 result = 0;

    for (const auto& entry :
         histogram) {
        result += entry.second;
    }

    return result;
}

/* =========================================================
 * Print histogram
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
 * Full histogram verification
 * ========================================================= */

static bool verify_full_histogram(
    u64 L,
    int p,
    int r,
    bool print_failure,
    Stats& stats
) {
    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 expected_length =
        power - L - 1;

    const auto exact =
        exact_local_histogram(
            L,
            p,
            r
        );

    const auto dp =
        local_digit_dp(
            L,
            p,
            r
        );

    bool ok = true;

    if (exact != dp) {
        ok = false;
        ++stats.histogram_failures;

        if (print_failure) {
            std::cout
                << "FAIL_HISTOGRAM "
                << "L=" << L
                << " p=" << p
                << " r=" << r
                << "\n";

            std::cout
                << "exact=";

            print_histogram(exact);

            std::cout
                << "\ndp=";

            print_histogram(dp);

            std::cout
                << "\n";
        }
    }

    const u64 exact_total =
        histogram_total(exact);

    const u64 dp_total =
        histogram_total(dp);

    if (
        exact_total != expected_length ||
        dp_total != expected_length
    ) {
        ok = false;
        ++stats.total_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TOTAL "
                << "L=" << L
                << " p=" << p
                << " r=" << r
                << " exact="
                << exact_total
                << " dp="
                << dp_total
                << " expected="
                << expected_length
                << "\n";
        }
    }

    if (!exact.empty()) {
        const u64 maximum =
            exact.rbegin()->first;

        if (
            maximum >
            static_cast<u64>(r)
        ) {
            ok = false;
            ++stats.max_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_MAX "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << " max="
                    << maximum
                    << "\n";
            }
        }
    }

    return ok;
}

/* =========================================================
 * Sampled point verification for a large local block
 * ========================================================= */

static bool verify_sampled_points(
    u64 L,
    int p,
    int r,
    std::mt19937_64& rng,
    bool print_failure,
    Stats& stats
) {
    const u64 power =
        prime_power(
            p,
            r
        );

    if (L + 1 >= power) {
        return true;
    }

    const u64 n =
        power + L;

    /*
     * Include deterministic boundary points.
     */
    std::vector<u64> samples;

    samples.push_back(
        L + 1
    );

    samples.push_back(
        power - 1
    );

    if (
        power - L - 1 > 2
    ) {
        samples.push_back(
            L + 2
        );

        samples.push_back(
            L + 3
        );

        samples.push_back(
            power - 2
        );
    }

    /*
     * Add random points.
     */
    constexpr int RANDOM_SAMPLES = 100;

    for (int i = 0;
         i < RANDOM_SAMPLES;
         ++i) {

        const u64 span =
            power - L - 1;

        samples.push_back(
            L + 1 +
            rng() % span
        );
    }

    bool ok = true;

    for (u64 y : samples) {
        if (
            y <= L ||
            y >= power
        ) {
            continue;
        }

        const u64 exact =
            binomial_valuation(
                n,
                y,
                p
            );

        const u64 borrow =
            kummer_borrow_count(
                n,
                y,
                p
            );

        if (
            exact != borrow
        ) {
            ok = false;
            ++stats.point_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_POINT "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << " y=" << y
                    << " exact="
                    << exact
                    << " borrow="
                    << borrow
                    << "\n";
            }

            break;
        }
    }

    return ok;
}

/* =========================================================
 * Full case
 * ========================================================= */

static bool verify_case(
    u64 L,
    int p,
    int r,
    bool full_histogram,
    std::mt19937_64& rng,
    bool print_failure,
    Stats& stats
) {
    ++stats.tested;

    bool ok = true;

    if (full_histogram) {
        ++stats.full_histogram_cases;

        if (
            !verify_full_histogram(
                L,
                p,
                r,
                print_failure,
                stats
            )
        ) {
            ok = false;
        }
    } else {
        ++stats.sampled_point_cases;

        if (
            !verify_sampled_points(
                L,
                p,
                r,
                rng,
                print_failure,
                stats
            )
        ) {
            ok = false;
        }

        /*
         * For large cases we also verify the DP's total
         * number of y values.
         */
        const u64 power =
            prime_power(
                p,
                r
            );

        const u64 expected =
            power - L - 1;

        const auto dp =
            local_digit_dp(
                L,
                p,
                r
            );

        const u64 actual =
            histogram_total(dp);

        if (actual != expected) {
            ok = false;
            ++stats.total_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_DP_TOTAL "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << " actual="
                    << actual
                    << " expected="
                    << expected
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
        << "START EXPERIMENT 334\n"
        << "LOCAL KUMMER TRANSFER-MATRIX HISTOGRAM\n"
        << "CAN THE GAP PROFILE BE COMPUTED FROM L DIGITS BY A BORROW DP?\n\n";

    Stats stats;

    std::mt19937_64 rng(
        0x334334334ULL
    );

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive full histograms for small blocks.
     *
     * Every L is tested.
     * --------------------------------------------------
     */

    constexpr u64 FULL_LIMIT = 2000;

    for (int p : PRIMES) {
        for (int r = 1; r <= 12; ++r) {
            const u64 power =
                prime_power(
                    p,
                    r
                );

            if (
                power > FULL_LIMIT
            ) {
                break;
            }

            for (
                u64 L = 0;
                L < power;
                ++L
            ) {
                verify_case(
                    L,
                    p,
                    r,
                    true,
                    rng,
                    stats.failed < 10,
                    stats
                );
            }
        }
    }

    /*
     * --------------------------------------------------
     * Phase 2:
     * Random full histograms.
     *
     * Full exact enumeration is only done where
     * p^r is reasonably small.
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_FULL_CASES = 1000;
    constexpr u64 RANDOM_FULL_LIMIT = 50000;

    u64 random_full = 0;

    while (
        random_full <
        RANDOM_FULL_CASES
    ) {
        const int p =
            PRIMES[
                rng() % PRIMES.size()
            ];

        const int r =
            1 +
            static_cast<int>(
                rng() % 12
            );

        const u64 power =
            prime_power(
                p,
                r
            );

        if (
            power > RANDOM_FULL_LIMIT
        ) {
            continue;
        }

        const u64 L =
            rng() % power;

        ++random_full;

        verify_case(
            L,
            p,
            r,
            true,
            rng,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 3:
     * Large local blocks.
     *
     * We do not enumerate the whole interval.
     * Instead we:
     *
     *   - evaluate the digit DP
     *   - verify its total
     *   - compare exact Legendre and Kummer values
     *     on 100+ selected points.
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_LARGE_CASES = 3000;

    for (
        u64 i = 0;
        i < RANDOM_LARGE_CASES;
        ++i
    ) {
        const int p =
            PRIMES[
                rng() % PRIMES.size()
            ];

        const int r =
            1 +
            static_cast<int>(
                rng() % 20
            );

        const u64 power =
            prime_power(
                p,
                r
            );

        /*
         * Stay inside uint64.
         */
        if (
            power == 0
        ) {
            continue;
        }

        const u64 L =
            rng() % power;

        verify_case(
            L,
            p,
            r,
            false,
            rng,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 4:
     * Targeted L digit structures.
     * --------------------------------------------------
     */

    u64 targeted = 0;

    for (int p : PRIMES) {
        for (int r = 1; r <= 16; ++r) {
            const u64 power =
                prime_power(
                    p,
                    r
                );

            if (
                power == 0
            ) {
                break;
            }

            std::vector<u64> candidates;

            candidates.push_back(0);

            if (power > 1) {
                candidates.push_back(
                    power - 1
                );
            }

            if (power > 3) {
                candidates.push_back(1);
                candidates.push_back(2);
                candidates.push_back(power / 2);
            }

            /*
             * Build structured digit patterns.
             */
            for (int pattern = 0;
                 pattern < 4;
                 ++pattern) {

                std::vector<u64> digits(
                    static_cast<std::size_t>(r),
                    0
                );

                for (int i = 0;
                     i < r;
                     ++i) {

                    if (pattern == 0) {
                        /*
                         * alternating 0 / max
                         */
                        if (i & 1) {
                            digits[i] =
                                static_cast<u64>(
                                    p - 1
                                );
                        }
                    }

                    if (pattern == 1) {
                        /*
                         * alternating max / 0
                         */
                        if (!(i & 1)) {
                            digits[i] =
                                static_cast<u64>(
                                    p - 1
                                );
                        }
                    }

                    if (pattern == 2) {
                        /*
                         * all ones
                         */
                        digits[i] = 1;
                    }

                    if (pattern == 3) {
                        /*
                         * digit equals i mod p
                         */
                        digits[i] =
                            static_cast<u64>(
                                i % p
                            );
                    }
                }

                u64 value = 0;
                u64 place = 1;
                bool overflow = false;

                for (u64 d : digits) {
                    if (
                        d != 0 &&
                        place >
                            UINT64_MAX / d
                    ) {
                        overflow = true;
                        break;
                    }

                    value += d * place;

                    if (
                        place >
                        UINT64_MAX /
                            static_cast<u64>(p)
                    ) {
                        overflow = true;
                        break;
                    }

                    place *=
                        static_cast<u64>(p);
                }

                if (
                    !overflow &&
                    value < power
                ) {
                    candidates.push_back(value);
                }
            }

            for (u64 L : candidates) {
                ++targeted;

                /*
                 * Full exact comparison only when practical.
                 */
                const bool full =
                    power <= FULL_LIMIT;

                verify_case(
                    L,
                    p,
                    r,
                    full,
                    rng,
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
        << "point_failures="
        << stats.point_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "max_failures="
        << stats.max_failures
        << "\n"
        << "full_histogram_cases="
        << stats.full_histogram_cases
        << "\n"
        << "sampled_point_cases="
        << stats.sampled_point_cases
        << "\n"
        << "targeted="
        << targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 334\n";

    return 0;
}