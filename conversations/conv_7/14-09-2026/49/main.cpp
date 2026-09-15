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

    u64 path_failures = 0;
    u64 local_failures = 0;
    u64 exact_failures = 0;
    u64 total_failures = 0;

    u64 targeted = 0;
};

/* =========================================================
 * Utility
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

    for (int i = 0;
         i < length;
         ++i) {
        digits[i] =
            n % prime;

        n /= prime;
    }

    return digits;
}

/* =========================================================
 * Histogram equality
 * ========================================================= */

static bool equal_histogram(
    const std::vector<u64>& a,
    const std::vector<u64>& b
) {
    const std::size_t n =
        std::max(
            a.size(),
            b.size()
        );

    for (
        std::size_t i = 0;
        i < n;
        ++i
    ) {
        const u64 x =
            i < a.size()
                ? a[i]
                : 0;

        const u64 y =
            i < b.size()
                ? b[i]
                : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
}

static u64 histogram_sum(
    const std::vector<u64>& h
) {
    u64 result = 0;

    for (u64 x : h) {
        result += x;
    }

    return result;
}

/* =========================================================
 * Compute ALL lower-borrow distributions B_h(j)
 * in one pass.
 *
 * all_dp0[h][j] and all_dp1[h][j]
 * are not both necessary; we immediately combine them
 * after each level because the next transition still needs
 * the separated states.
 *
 * The returned vector contains:
 *
 * distributions[h][j] = B_h(j)
 *
 * for every h=0..r.
 * ========================================================= */

static std::vector<
    std::vector<u64>
>
all_borrow_distributions(
    u64 L,
    int p,
    int r
) {
    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    std::vector<
        std::vector<u64>
    > result(
        static_cast<std::size_t>(
            r + 1
        )
    );

    /*
     * Current states.
     */
    std::vector<u64> dp0(
        static_cast<std::size_t>(r + 1),
        0
    );

    std::vector<u64> dp1(
        static_cast<std::size_t>(r + 1),
        0
    );

    dp0[0] = 1;

    /*
     * h=0.
     */
    result[0] = {1};

    for (int i = 0;
         i < r;
         ++i) {

        const u64 Li =
            digits[
                static_cast<
                    std::size_t
                >(i)
            ];

        std::vector<u64> next0(
            static_cast<std::size_t>(r + 1),
            0
        );

        std::vector<u64> next1(
            static_cast<std::size_t>(r + 1),
            0
        );

        /*
         * Only j<=i can be nonzero.
         */
        for (int j = 0;
             j <= i;
             ++j) {

            if (dp0[j] != 0) {
                next0[j] +=
                    dp0[j] *
                    (Li + 1);

                next1[j + 1] +=
                    dp0[j] *
                    (
                        static_cast<u64>(p)
                        - Li
                        - 1
                    );
            }

            if (dp1[j] != 0) {
                next0[j] +=
                    dp1[j] * Li;

                next1[j + 1] +=
                    dp1[j] *
                    (
                        static_cast<u64>(p)
                        - Li
                    );
            }
        }

        dp0.swap(next0);
        dp1.swap(next1);

        result[
            static_cast<
                std::size_t
            >(i + 1)
        ].assign(
            static_cast<
                std::size_t
            >(i + 2),
            0
        );

        for (int j = 0;
             j <= i + 1;
             ++j) {

            result[
                static_cast<
                    std::size_t
                >(i + 1)
            ][j] =
                dp0[j] +
                dp1[j];
        }
    }

    return result;
}

/* =========================================================
 * Explicit path expansion for ONE h
 *
 * This is deliberately restricted to small h.
 * ========================================================= */

static std::vector<u64>
borrow_distribution_paths(
    const std::vector<u64>& digits,
    int p,
    int h
) {
    std::vector<u64> result(
        static_cast<
            std::size_t
        >(h + 1),
        0
    );

    if (h == 0) {
        result[0] = 1;
        return result;
    }

    /*
     * Enumerate 2^h paths.
     */
    const u64 paths =
        1ULL << h;

    for (
        u64 mask = 0;
        mask < paths;
        ++mask
    ) {
        u64 weight = 1;

        int borrow = 0;
        int count = 0;

        bool valid = true;

        for (int i = 0;
             i < h;
             ++i) {

            const int next_borrow =
                static_cast<int>(
                    (mask >> i) & 1ULL
                );

            const u64 Li =
                digits[
                    static_cast<
                        std::size_t
                    >(i)
                ];

            u64 transition = 0;

            if (
                borrow == 0 &&
                next_borrow == 0
            ) {
                transition =
                    Li + 1;
            } else if (
                borrow == 0 &&
                next_borrow == 1
            ) {
                transition =
                    static_cast<u64>(p)
                    - Li
                    - 1;
            } else if (
                borrow == 1 &&
                next_borrow == 0
            ) {
                transition = Li;
            } else {
                transition =
                    static_cast<u64>(p)
                    - Li;
            }

            if (transition == 0) {
                valid = false;
                break;
            }

            weight *= transition;

            count +=
                next_borrow;

            borrow =
                next_borrow;
        }

        if (valid) {
            result[
                static_cast<
                    std::size_t
                >(count)
            ] += weight;
        }
    }

    return result;
}

/* =========================================================
 * Local formula using precomputed B_h
 * ========================================================= */

static std::vector<u64>
local_from_all_B(
    u64 L,
    int p,
    int r,
    const std::vector<
        std::vector<u64>
    >& B
) {
    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    std::vector<u64> result(
        static_cast<
            std::size_t
        >(r + 1),
        0
    );

    for (int h = 0;
         h < r;
         ++h) {

        const u64 Lh =
            digits[
                static_cast<
                    std::size_t
                >(h)
            ];

        const u64 high_choices =
            static_cast<u64>(p)
            - Lh
            - 1;

        if (high_choices == 0) {
            continue;
        }

        const int fixed =
            r - h;

        for (
            int j = 0;
            j < static_cast<int>(
                B[h].size()
            );
            ++j
        ) {
            if (
                B[h][
                    static_cast<
                        std::size_t
                    >(j)
                ] == 0
            ) {
                continue;
            }

            const int valuation =
                fixed + j;

            result[
                static_cast<
                    std::size_t
                >(valuation)
            ] +=
                high_choices *
                B[h][
                    static_cast<
                        std::size_t
                    >(j)
                ];
        }
    }

    while (
        result.size() > 1 &&
        result.back() == 0
    ) {
        result.pop_back();
    }

    return result;
}

/* =========================================================
 * Exact local histogram
 *
 * Only used for modest p^r.
 * ========================================================= */

static std::vector<u64>
exact_local_histogram(
    u64 L,
    int p,
    int r
) {
    const u64 power =
        prime_power(
            p,
            r
        );

    const u64 n =
        power + L;

    std::vector<u64> result(
        static_cast<
            std::size_t
        >(r + 1),
        0
    );

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

        ++result[
            static_cast<
                std::size_t
            >(value)
        ];
    }

    while (
        result.size() > 1 &&
        result.back() == 0
    ) {
        result.pop_back();
    }

    return result;
}

/* =========================================================
 * Verify one local block
 * ========================================================= */

static bool verify_local(
    u64 L,
    int p,
    int r,
    bool exact_check,
    bool path_check,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    /*
     * ONE DP pass computes all B_h.
     */
    const auto B =
        all_borrow_distributions(
            L,
            p,
            r
        );

    /*
     * --------------------------------------------------
     * Explicit path validation.
     *
     * Only a few small h are needed. We do not enumerate
     * all 2^h paths for every h.
     * --------------------------------------------------
     */

    if (path_check) {
        const int max_h =
            std::min(
                r,
                8
            );

        for (
            int h = 0;
            h <= max_h;
            ++h
        ) {
            const auto paths =
                borrow_distribution_paths(
                    digits,
                    p,
                    h
                );

            if (
                !equal_histogram(
                    B[h],
                    paths
                )
            ) {
                ok = false;

                ++stats.path_failures;

                if (print_failure) {
                    std::cout
                        << "FAIL_PATH "
                        << "L=" << L
                        << " p=" << p
                        << " r=" << r
                        << " h=" << h
                        << "\n";
                }

                break;
            }
        }
    }

    /*
     * --------------------------------------------------
     * Local coefficient formula
     * --------------------------------------------------
     */

    const auto local =
        local_from_all_B(
            L,
            p,
            r,
            B
        );

    const u64 expected_length =
        prime_power(
            p,
            r
        ) -
        L -
        1;

    if (
        histogram_sum(local) !=
        expected_length
    ) {
        ok = false;

        ++stats.total_failures;

        if (print_failure) {
            std::cout
                << "FAIL_TOTAL "
                << "L=" << L
                << " p=" << p
                << " r=" << r
                << " actual="
                << histogram_sum(local)
                << " expected="
                << expected_length
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Exact local comparison
     * --------------------------------------------------
     */

    if (exact_check) {
        const auto exact =
            exact_local_histogram(
                L,
                p,
                r
            );

        if (
            !equal_histogram(
                local,
                exact
            )
        ) {
            ok = false;

            ++stats.exact_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_EXACT "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
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
        << "START EXPERIMENT 338\n"
        << "EXPLICIT BORROW-PATH EXPANSION\n"
        << "DO BINARY BORROW HISTORIES REPRODUCE THE KUMMER LAYERS?\n\n";

    Stats stats;

    std::mt19937_64 rng(
        0x338338338ULL
    );

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive local blocks.
     *
     * All explicit exact checks have p^r <= 1000.
     * Path enumeration only goes up to h<=8.
     * --------------------------------------------------
     */

    constexpr u64 EXHAUSTIVE_LIMIT = 1000;

    for (int p : PRIMES) {
        for (int r = 1;
             r <= 16;
             ++r) {

            const u64 power =
                prime_power(
                    p,
                    r
                );

            if (
                power > EXHAUSTIVE_LIMIT
            ) {
                break;
            }

            for (
                u64 L = 0;
                L < power;
                ++L
            ) {
                ++stats.tested;

                verify_local(
                    L,
                    p,
                    r,
                    true,
                    true,
                    stats.failed < 10,
                    stats
                );
            }
        }
    }

    /*
     * --------------------------------------------------
     * Phase 2:
     * Random local blocks.
     *
     * Much larger r is allowed, but the explicit path
     * check still stops at h<=8.
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_CASES = 5000;

    for (
        u64 i = 0;
        i < RANDOM_CASES;
        ++i
    ) {
        const int p =
            PRIMES[
                rng() % PRIMES.size()
            ];

        const int r =
            1 +
            static_cast<int>(
                rng() % 24
            );

        const u64 power =
            prime_power(
                p,
                r
            );

        if (power == 0) {
            continue;
        }

        const u64 L =
            rng() % power;

        ++stats.tested;

        /*
         * Exact enumeration is only used for small blocks.
         */
        const bool exact =
            power <= 50000;

        verify_local(
            L,
            p,
            r,
            exact,
            true,
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
        for (int r = 1;
             r <= 20;
             ++r) {

            const u64 power =
                prime_power(
                    p,
                    r
                );

            if (power == 0) {
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
             * Structured digit patterns.
             */
            for (int pattern = 0;
                 pattern < 5;
                 ++pattern) {

                std::vector<u64> digits(
                    static_cast<
                        std::size_t
                    >(r),
                    0
                );

                for (int i = 0;
                     i < r;
                     ++i) {

                    if (pattern == 0) {
                        /*
                         * 0,max,0,max,...
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
                         * max,0,max,0,...
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
                         * increasing digits
                         */
                        digits[i] =
                            static_cast<u64>(
                                i % p
                            );
                    }

                    if (pattern == 4) {
                        /*
                         * sparse high half
                         */
                        if (
                            i >= r / 2
                        ) {
                            digits[i] = 1;
                        }
                    }
                }

                u64 value = 0;

                if (
                    build_from_digits(
                        digits,
                        p,
                        value
                    ) &&
                    value < power
                ) {
                    candidates.push_back(
                        value
                    );
                }
            }

            for (u64 L : candidates) {
                ++stats.targeted;
                ++stats.tested;

                const bool exact =
                    power <= 1000;

                verify_local(
                    L,
                    p,
                    r,
                    exact,
                    true,
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
        << "path_failures="
        << stats.path_failures
        << "\n"
        << "local_failures="
        << stats.local_failures
        << "\n"
        << "exact_failures="
        << stats.exact_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "targeted="
        << stats.targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 338\n";

    return 0;
}