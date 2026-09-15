#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

using Poly = std::vector<u128>;
using SPoly = std::vector<i128>;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 coefficient_formula_failures = 0;
    u64 matching_degree_failures = 0;
    u64 exact_failures = 0;
    u64 total_failures = 0;

    u64 matching_cases = 0;
    u64 matching_terms = 0;

    u64 targeted = 0;
};

/* =========================================================
 * Printing
 * ========================================================= */

static void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    char buffer[64];
    int pos = 0;

    while (x > 0) {
        buffer[pos++] =
            static_cast<char>(
                '0' + static_cast<unsigned>(x % 10)
            );

        x /= 10;
    }

    while (pos > 0) {
        std::cout << buffer[--pos];
    }
}

/* =========================================================
 * Unsigned polynomial helpers
 * ========================================================= */

static void trim_poly(Poly& p) {
    while (p.size() > 1 && p.back() == 0) {
        p.pop_back();
    }

    if (p.empty()) {
        p.push_back(0);
    }
}

static bool equal_poly(
    const Poly& a,
    const Poly& b
) {
    const std::size_t n =
        std::max(a.size(), b.size());

    for (std::size_t i = 0; i < n; ++i) {
        const u128 x =
            i < a.size() ? a[i] : 0;

        const u128 y =
            i < b.size() ? b[i] : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
}

static Poly poly_add(
    const Poly& a,
    const Poly& b
) {
    Poly result(
        std::max(a.size(), b.size()),
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] += a[i];
    }

    for (std::size_t i = 0; i < b.size(); ++i) {
        result[i] += b[i];
    }

    trim_poly(result);
    return result;
}

static Poly scalar_mul(
    const Poly& a,
    u128 scale
) {
    Poly result = a;

    for (u128& x : result) {
        x *= scale;
    }

    trim_poly(result);
    return result;
}

/* =========================================================
 * Signed polynomial helpers
 * ========================================================= */

static void trim_signed(SPoly& p) {
    while (p.size() > 1 && p.back() == 0) {
        p.pop_back();
    }

    if (p.empty()) {
        p.push_back(0);
    }
}

static bool signed_equals_unsigned(
    const SPoly& a,
    const Poly& b
) {
    const std::size_t n =
        std::max(a.size(), b.size());

    for (std::size_t i = 0; i < n; ++i) {
        const i128 x =
            i < a.size()
                ? a[i]
                : static_cast<i128>(0);

        const i128 y =
            i < b.size()
                ? static_cast<i128>(b[i])
                : static_cast<i128>(0);

        if (x != y) {
            return false;
        }
    }

    return true;
}

static void add_signed_term(
    SPoly& result,
    const Poly& coefficients,
    std::size_t degree_shift,
    u128 multiplier,
    bool negative
) {
    const std::size_t required =
        coefficients.size() + degree_shift;

    if (result.size() < required) {
        result.resize(required, 0);
    }

    for (
        std::size_t k = 0;
        k < coefficients.size();
        ++k
    ) {
        const i128 value =
            static_cast<i128>(
                coefficients[k] * multiplier
            );

        const std::size_t degree =
            k + degree_shift;

        if (negative) {
            result[degree] -= value;
        } else {
            result[degree] += value;
        }
    }
}

/* =========================================================
 * Base-p helpers
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

static std::vector<u64> digits_of(
    u64 n,
    int p,
    int r
) {
    std::vector<u64> digits(
        static_cast<std::size_t>(r),
        0
    );

    const u64 P =
        static_cast<u64>(p);

    for (int i = 0; i < r; ++i) {
        digits[
            static_cast<std::size_t>(i)
        ] = n % P;

        n /= P;
    }

    return digits;
}

/* =========================================================
 * Two-state Kummer polynomial system
 *
 * A_{i+1} =
 *     (d_i+1) A_i + d_i C_i
 *
 * C_{i+1} =
 *     z(p-d_i-1) A_i
 *     + z(p-d_i) C_i
 * ========================================================= */

static void compute_history(
    u64 L,
    int p,
    int r,
    std::vector<Poly>& A,
    std::vector<Poly>& C
) {
    const auto digits =
        digits_of(L, p, r);

    A.resize(
        static_cast<std::size_t>(r + 1)
    );

    C.resize(
        static_cast<std::size_t>(r + 1)
    );

    A[0] = {1};
    C[0] = {0};

    for (int i = 0; i < r; ++i) {
        const u128 d =
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(i)
                ]
            );

        A[i + 1] =
            poly_add(
                scalar_mul(
                    A[i],
                    d + 1
                ),
                scalar_mul(
                    C[i],
                    d
                )
            );

        Poly zA(
            A[i].size() + 1,
            0
        );

        for (
            std::size_t k = 0;
            k < A[i].size();
            ++k
        ) {
            zA[k + 1] = A[i][k];
        }

        Poly zC(
            C[i].size() + 1,
            0
        );

        for (
            std::size_t k = 0;
            k < C[i].size();
            ++k
        ) {
            zC[k + 1] = C[i][k];
        }

        C[i + 1] =
            poly_add(
                scalar_mul(
                    zA,
                    static_cast<u128>(p)
                    - d
                    - 1
                ),
                scalar_mul(
                    zC,
                    static_cast<u128>(p)
                    - d
                )
            );
    }
}

/* =========================================================
 * D_r A_r
 * ========================================================= */

static Poly scaled_A(
    const std::vector<Poly>& A,
    const std::vector<u64>& digits,
    int r
) {
    u128 scale = 1;

    for (int i = 0; i <= r - 2; ++i) {
        scale *=
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(i)
                ]
            );
    }

    return scalar_mul(
        A[
            static_cast<std::size_t>(r)
        ],
        scale
    );
}

/* =========================================================
 * Matching factors
 *
 * a_0 = alpha_0
 *
 * a_i = alpha_i + beta_i z, i >= 1
 *
 * edge_i = gamma_i z
 * ========================================================= */

static u128 alpha(
    const std::vector<u64>& digits,
    int i
) {
    if (i == 0) {
        return static_cast<u128>(
            digits[0] + 1
        );
    }

    const u128 d_prev =
        static_cast<u128>(
            digits[
                static_cast<std::size_t>(i - 1)
            ]
        );

    const u128 d_curr =
        static_cast<u128>(
            digits[
                static_cast<std::size_t>(i)
            ]
        );

    return
        d_prev *
        (d_curr + 1);
}

static u128 beta(
    const std::vector<u64>& digits,
    int p,
    int i
) {
    if (i == 0) {
        return 0;
    }

    const u128 d_prev =
        static_cast<u128>(
            digits[
                static_cast<std::size_t>(i - 1)
            ]
        );

    const u128 d_curr =
        static_cast<u128>(
            digits[
                static_cast<std::size_t>(i)
            ]
        );

    return
        d_curr *
        (
            static_cast<u128>(p)
            - d_prev
        );
}

static u128 gamma(
    const std::vector<u64>& digits,
    int p,
    int edge
) {
    u128 result =
        static_cast<u128>(
            digits[
                static_cast<std::size_t>(edge)
            ]
        ) *
        static_cast<u128>(p);

    if (edge >= 2) {
        result *=
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(
                        edge - 2
                    )
                ]
            );
    }

    return result;
}

/* =========================================================
 * Product of unmatched diagonal factors
 *
 * Computes
 *
 * prod_{i in U} (alpha_i + beta_i z)
 *
 * directly at coefficient level.
 * ========================================================= */

static Poly unmatched_distribution(
    const std::vector<int>& unmatched,
    const std::vector<u64>& digits,
    int p
) {
    Poly dp{1};

    for (int vertex : unmatched) {
        const u128 a =
            alpha(
                digits,
                vertex
            );

        const u128 b =
            beta(
                digits,
                p,
                vertex
            );

        Poly next(
            dp.size() + 1,
            0
        );

        for (
            std::size_t k = 0;
            k < dp.size();
            ++k
        ) {
            next[k] +=
                dp[k] * a;

            next[k + 1] +=
                dp[k] * b;
        }

        dp =
            std::move(next);
    }

    trim_poly(dp);
    return dp;
}

/* =========================================================
 * Explicit path matching enumeration
 * ========================================================= */

static void enumerate_matchings(
    int r,
    int p,
    const std::vector<u64>& digits,
    int next_edge,
    std::vector<int>& selected,
    SPoly& result,
    u64& matching_count
) {
    /*
     * A matching is complete once there are no edges
     * left to consider.
     */
    if (next_edge >= r) {
        ++matching_count;

        std::vector<bool> covered(
            static_cast<std::size_t>(r),
            false
        );

        u128 edge_weight = 1;

        for (int edge : selected) {
            covered[
                static_cast<std::size_t>(
                    edge - 1
                )
            ] = true;

            covered[
                static_cast<std::size_t>(
                    edge
                )
            ] = true;

            edge_weight *=
                gamma(
                    digits,
                    p,
                    edge
                );
        }

        /*
         * Find unmatched vertices.
         */
        std::vector<int> unmatched;

        for (int vertex = 0;
             vertex < r;
             ++vertex) {

            if (
                !covered[
                    static_cast<
                        std::size_t
                    >(vertex)
                ]
            ) {
                unmatched.push_back(
                    vertex
                );
            }
        }

        /*
         * Coefficients from the unmatched diagonal
         * factors.
         */
        const Poly diagonal =
            unmatched_distribution(
                unmatched,
                digits,
                p
            );

        /*
         * Every selected edge contributes one z.
         */
        const std::size_t shift =
            selected.size();

        /*
         * Sign = (-1)^|M|.
         */
        const bool negative =
            (selected.size() & 1U) != 0;

        add_signed_term(
            result,
            diagonal,
            shift,
            edge_weight,
            negative
        );

        return;
    }

    /*
     * Do not select this edge.
     */
    enumerate_matchings(
        r,
        p,
        digits,
        next_edge + 1,
        selected,
        result,
        matching_count
    );

    /*
     * Select this edge.
     */
    selected.push_back(next_edge);

    enumerate_matchings(
        r,
        p,
        digits,
        next_edge + 2,
        selected,
        result,
        matching_count
    );

    selected.pop_back();
}

static SPoly matching_coefficient_formula(
    u64 L,
    int p,
    int r,
    u64& matching_count
) {
    const auto digits =
        digits_of(L, p, r);

    SPoly result{0};

    std::vector<int> selected;

    matching_count = 0;

    enumerate_matchings(
        r,
        p,
        digits,
        1,
        selected,
        result,
        matching_count
    );

    trim_signed(result);

    return result;
}

/* =========================================================
 * Exact valuation distribution
 * ========================================================= */

static u64 factorial_vp(
    u64 n,
    int p
) {
    u64 result = 0;

    const u64 P =
        static_cast<u64>(p);

    while (n > 0) {
        n /= P;
        result += n;
    }

    return result;
}

static u64 binomial_vp(
    u64 n,
    u64 k,
    int p
) {
    return
        factorial_vp(n, p)
        - factorial_vp(k, p)
        - factorial_vp(n - k, p);
}

static Poly exact_local(
    u64 L,
    int p,
    int r
) {
    const u64 power =
        prime_power(p, r);

    const u64 n =
        power + L;

    Poly result{0};

    for (
        u64 y = L + 1;
        y < power;
        ++y
    ) {
        const u64 v =
            binomial_vp(
                n,
                y,
                p
            );

        if (
            result.size() <= v
        ) {
            result.resize(
                static_cast<
                    std::size_t
                >(v + 1),
                0
            );
        }

        ++result[
            static_cast<
                std::size_t
            >(v)
        ];
    }

    trim_poly(result);
    return result;
}

/* =========================================================
 * Local polynomial from A_h + C_h
 * ========================================================= */

static Poly local_from_history(
    const std::vector<Poly>& A,
    const std::vector<Poly>& C,
    const std::vector<u64>& digits,
    int p,
    int r
) {
    Poly result{0};

    for (int h = 0; h < r; ++h) {
        Poly B =
            poly_add(
                A[
                    static_cast<
                        std::size_t
                    >(h)
                ],
                C[
                    static_cast<
                        std::size_t
                    >(h)
                ]
            );

        const u128 choices =
            static_cast<u128>(p)
            -
            static_cast<u128>(
                digits[
                    static_cast<
                        std::size_t
                    >(h)
                ]
            )
            - 1;

        if (choices == 0) {
            continue;
        }

        const std::size_t shift =
            static_cast<std::size_t>(
                r - h
            );

        Poly shifted(
            B.size() + shift,
            0
        );

        for (
            std::size_t k = 0;
            k < B.size();
            ++k
        ) {
            shifted[k + shift] =
                B[k];
        }

        result =
            poly_add(
                result,
                scalar_mul(
                    shifted,
                    choices
                )
            );
    }

    trim_poly(result);
    return result;
}

/* =========================================================
 * Test one case
 * ========================================================= */

static bool verify_case(
    u64 L,
    int p,
    int r,
    bool exact_check,
    bool matching_check,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    const auto digits =
        digits_of(L, p, r);

    std::vector<Poly> A;
    std::vector<Poly> C;

    compute_history(
        L,
        p,
        r,
        A,
        C
    );

    /*
     * Reference polynomial:
     *
     * D_r A_r
     */
    const Poly reference =
        scaled_A(
            A,
            digits,
            r
        );

    /*
     * --------------------------------------------------
     * Coefficient-level matching formula
     * --------------------------------------------------
     */

    if (matching_check) {
        ++stats.matching_cases;

        u64 matching_count = 0;

        const SPoly formula =
            matching_coefficient_formula(
                L,
                p,
                r,
                matching_count
            );

        stats.matching_terms +=
            matching_count;

        if (
            !signed_equals_unsigned(
                formula,
                reference
            )
        ) {
            ok = false;
            ++stats.coefficient_formula_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_COEFFICIENT_FORMULA "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << " matchings="
                    << matching_count
                    << "\n";
            }
        }

        /*
         * Explicit coefficient-by-coefficient check.
         */
        const std::size_t n =
            std::max(
                formula.size(),
                reference.size()
            );

        for (
            std::size_t k = 0;
            k < n;
            ++k
        ) {
            const i128 lhs =
                k < formula.size()
                    ? formula[k]
                    : static_cast<i128>(0);

            const i128 rhs =
                k < reference.size()
                    ? static_cast<i128>(
                        reference[k]
                    )
                    : static_cast<i128>(0);

            if (lhs != rhs) {
                ok = false;
                ++stats.matching_degree_failures;

                if (print_failure) {
                    std::cout
                        << "FAIL_COEFFICIENT "
                        << "L=" << L
                        << " p=" << p
                        << " r=" << r
                        << " k=" << k
                        << "\n";

                    std::cout
                        << "expected=";

                    print_u128(
                        k < reference.size()
                            ? reference[k]
                            : 0
                    );

                    std::cout
                        << "\n";

                    break;
                }
            }
        }
    }

    /*
     * --------------------------------------------------
     * Exact local distribution
     * --------------------------------------------------
     */

    if (exact_check) {
        const Poly local =
            local_from_history(
                A,
                C,
                digits,
                p,
                r
            );

        const Poly exact =
            exact_local(
                L,
                p,
                r
            );

        if (!equal_poly(local, exact)) {
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

        u64 actual = 0;

        for (u128 x : local) {
            actual +=
                static_cast<u64>(x);
        }

        const u64 expected_total =
            prime_power(p, r)
            - L
            - 1;

        if (actual != expected_total) {
            ok = false;
            ++stats.total_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_TOTAL "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << " actual="
                    << actual
                    << " expected="
                    << expected_total
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
        << "START EXPERIMENT 342\n"
        << "COEFFICIENT-LEVEL MATCHING FORMULA\n"
        << "CAN EACH KUMMER COEFFICIENT BE WRITTEN "
           "DIRECTLY AS A MATCHING SUM?\n\n";

    Stats stats;

    std::mt19937_64 rng(
        0x342342342ULL
    );

    /*
     * --------------------------------------------------
     * Exhaustive small cases
     * --------------------------------------------------
     */

    constexpr u64 SMALL_LIMIT = 1200;

    for (int p : PRIMES) {
        for (int r = 1; r <= 14; ++r) {
            const u64 power =
                prime_power(p, r);

            if (power > SMALL_LIMIT) {
                break;
            }

            for (u64 L = 0; L < power; ++L) {
                ++stats.tested;

                verify_case(
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
     * Random matching cases
     *
     * Explicit matching enumeration grows
     * approximately like Fibonacci(r).
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_CASES = 5000;

    for (
        u64 n = 0;
        n < RANDOM_CASES;
        ++n
    ) {
        const int p =
            PRIMES[
                rng() % PRIMES.size()
            ];

        const int r =
            3 +
            static_cast<int>(
                rng() % 13
            );

        const u64 power =
            prime_power(p, r);

        const u64 L =
            rng() % power;

        ++stats.tested;

        verify_case(
            L,
            p,
            r,
            power <= 50000,
            true,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Larger cases
     * --------------------------------------------------
     */

    constexpr u64 LARGE_CASES = 1800;

    for (
        u64 n = 0;
        n < LARGE_CASES;
        ++n
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
            prime_power(p, r);

        const u64 L =
            rng() % power;

        ++stats.tested;

        const bool matching =
            (r <= 16);

        verify_case(
            L,
            p,
            r,
            power <= 50000,
            matching,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Targeted digit structures
     * --------------------------------------------------
     */

    for (int p : PRIMES) {
        for (int r = 1; r <= 16; ++r) {
            const u64 power =
                prime_power(p, r);

            std::vector<u64> cases;

            cases.push_back(0);

            if (power > 1) {
                cases.push_back(1);
                cases.push_back(power - 1);
            }

            if (power > 4) {
                cases.push_back(2);
                cases.push_back(power / 2);
                cases.push_back(power / 2 + 1);
            }

            /*
             * Equal digits.
             */
            for (
                u64 digit = 0;
                digit < static_cast<u64>(p);
                ++digit
            ) {
                u64 value = 0;
                u64 place = 1;

                for (int i = 0; i < r; ++i) {
                    value +=
                        digit * place;

                    place *=
                        static_cast<u64>(p);
                }

                cases.push_back(value);
            }

            /*
             * Alternating 0 / (p-1).
             */
            for (int mode = 0; mode < 2; ++mode) {
                u64 value = 0;
                u64 place = 1;

                for (int i = 0; i < r; ++i) {
                    const u64 digit =
                        ((i & 1) == mode)
                            ? static_cast<u64>(p - 1)
                            : 0;

                    value +=
                        digit * place;

                    place *=
                        static_cast<u64>(p);
                }

                cases.push_back(value);
            }

            /*
             * Single nonzero digit.
             */
            for (int pos = 0; pos < r; ++pos) {
                u64 value = 0;
                u64 place = 1;

                for (int i = 0; i < r; ++i) {
                    if (i == pos) {
                        value +=
                            static_cast<u64>(p - 1)
                            * place;
                    }

                    place *=
                        static_cast<u64>(p);
                }

                cases.push_back(value);
            }

            for (u64 L : cases) {
                ++stats.tested;
                ++stats.targeted;

                verify_case(
                    L,
                    p,
                    r,
                    power <= 50000,
                    r <= 16,
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
        << "coefficient_formula_failures="
        << stats.coefficient_formula_failures
        << "\n"
        << "matching_degree_failures="
        << stats.matching_degree_failures
        << "\n"
        << "exact_failures="
        << stats.exact_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "matching_cases="
        << stats.matching_cases
        << "\n"
        << "matching_terms="
        << stats.matching_terms
        << "\n"
        << "targeted="
        << stats.targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 342\n";

    return 0;
}