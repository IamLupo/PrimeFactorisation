#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

using UPoly = std::vector<u128>;
using SPoly = std::vector<i128>;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 matching_failures = 0;
    u64 coefficient_failures = 0;
    u64 scaled_a_failures = 0;
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

static void trim_u(UPoly& p) {
    while (p.size() > 1 && p.back() == 0) {
        p.pop_back();
    }

    if (p.empty()) {
        p.push_back(0);
    }
}

static bool equal_u(
    const UPoly& a,
    const UPoly& b
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

static UPoly add_u(
    const UPoly& a,
    const UPoly& b
) {
    UPoly r(
        std::max(a.size(), b.size()),
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        r[i] += a[i];
    }

    for (std::size_t i = 0; i < b.size(); ++i) {
        r[i] += b[i];
    }

    trim_u(r);
    return r;
}

static UPoly scalar_u(
    const UPoly& a,
    u128 s
) {
    UPoly r = a;

    for (u128& x : r) {
        x *= s;
    }

    trim_u(r);
    return r;
}

static UPoly mul_u(
    const UPoly& a,
    const UPoly& b
) {
    if (
        (a.size() == 1 && a[0] == 0) ||
        (b.size() == 1 && b[0] == 0)
    ) {
        return UPoly{0};
    }

    UPoly r(
        a.size() + b.size() - 1,
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] +=
                a[i] * b[j];
        }
    }

    trim_u(r);
    return r;
}

/* =========================================================
 * Signed polynomial helpers
 * ========================================================= */

static void trim_s(SPoly& p) {
    while (p.size() > 1 && p.back() == 0) {
        p.pop_back();
    }

    if (p.empty()) {
        p.push_back(0);
    }
}

static SPoly signed_from_unsigned(
    const UPoly& p
) {
    SPoly r(p.size());

    for (std::size_t i = 0; i < p.size(); ++i) {
        r[i] =
            static_cast<i128>(p[i]);
    }

    trim_s(r);
    return r;
}

static SPoly add_s(
    const SPoly& a,
    const SPoly& b
) {
    SPoly r(
        std::max(a.size(), b.size()),
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        r[i] += a[i];
    }

    for (std::size_t i = 0; i < b.size(); ++i) {
        r[i] += b[i];
    }

    trim_s(r);
    return r;
}

static SPoly neg_s(
    const SPoly& a
) {
    SPoly r = a;

    for (i128& x : r) {
        x = -x;
    }

    trim_s(r);
    return r;
}

static SPoly mul_s(
    const SPoly& a,
    const SPoly& b
) {
    if (
        (a.size() == 1 && a[0] == 0) ||
        (b.size() == 1 && b[0] == 0)
    ) {
        return SPoly{0};
    }

    SPoly r(
        a.size() + b.size() - 1,
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] +=
                a[i] * b[j];
        }
    }

    trim_s(r);
    return r;
}

static bool equal_signed_unsigned(
    const SPoly& a,
    const UPoly& b
) {
    const std::size_t n =
        std::max(a.size(), b.size());

    for (std::size_t i = 0; i < n; ++i) {
        const i128 x =
            i < a.size() ? a[i] : 0;

        const i128 y =
            i < b.size()
                ? static_cast<i128>(b[i])
                : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
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
    std::vector<u64> d(
        static_cast<std::size_t>(r),
        0
    );

    const u64 P =
        static_cast<u64>(p);

    for (int i = 0; i < r; ++i) {
        d[
            static_cast<std::size_t>(i)
        ] = n % P;

        n /= P;
    }

    return d;
}

/* =========================================================
 * Two-state Kummer system
 * ========================================================= */

static void compute_history(
    u64 L,
    int p,
    int r,
    std::vector<UPoly>& A,
    std::vector<UPoly>& C
) {
    const auto d =
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
        const u128 digit =
            static_cast<u128>(
                d[
                    static_cast<std::size_t>(i)
                ]
            );

        /*
         * A_{i+1}
         */
        A[i + 1] =
            scalar_u(
                A[i],
                digit + 1
            );

        A[i + 1] =
            add_u(
                A[i + 1],
                scalar_u(
                    C[i],
                    digit
                )
            );

        /*
         * C_{i+1}
         */
        UPoly zA(
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

        UPoly zC(
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
            scalar_u(
                zA,
                static_cast<u128>(p)
                - digit
                - 1
            );

        C[i + 1] =
            add_u(
                C[i + 1],
                scalar_u(
                    zC,
                    static_cast<u128>(p)
                    - digit
                )
            );
    }
}

/* =========================================================
 * Exact local polynomial
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

static UPoly exact_local(
    u64 L,
    int p,
    int r
) {
    const u64 power =
        prime_power(p, r);

    const u64 n =
        power + L;

    UPoly result{0};

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
                static_cast<std::size_t>(v + 1),
                0
            );
        }

        ++result[
            static_cast<std::size_t>(v)
        ];
    }

    trim_u(result);
    return result;
}

/* =========================================================
 * Digit-dependent factors
 *
 * a_0 = d_0 + 1
 *
 * a_i =
 * d_{i-1}(d_i+1)
 * + d_i(p-d_{i-1}) z
 *
 * s_1 = d_1 p z
 *
 * s_i =
 * d_i p z d_{i-2}
 *
 * for i >= 2.
 * ========================================================= */

static UPoly diagonal_factor(
    const std::vector<u64>& d,
    int p,
    int i
) {
    if (i == 0) {
        return UPoly{
            static_cast<u128>(
                d[0] + 1
            )
        };
    }

    UPoly result(2, 0);

    const u128 d_prev =
        static_cast<u128>(
            d[
                static_cast<std::size_t>(i - 1)
            ]
        );

    const u128 d_curr =
        static_cast<u128>(
            d[
                static_cast<std::size_t>(i)
            ]
        );

    result[0] =
        d_prev *
        (d_curr + 1);

    result[1] =
        d_curr *
        (
            static_cast<u128>(p)
            - d_prev
        );

    trim_u(result);
    return result;
}

static UPoly edge_factor(
    const std::vector<u64>& d,
    int p,
    int edge
) {
    /*
     * edge is between vertices edge-1 and edge.
     */
    const u128 d_curr =
        static_cast<u128>(
            d[
                static_cast<std::size_t>(edge)
            ]
        );

    u128 coefficient =
        d_curr *
        static_cast<u128>(p);

    if (edge >= 2) {
        coefficient *=
            static_cast<u128>(
                d[
                    static_cast<std::size_t>(edge - 2)
                ]
            );
    }

    return UPoly{
        0,
        coefficient
    };
}

/* =========================================================
 * Explicit matching enumeration
 *
 * A matching of a path is represented by the set of
 * selected edges.
 *
 * Vertices incident to a selected edge are removed from
 * the diagonal product.
 *
 * Each selected edge contributes a NEGATIVE edge factor.
 * ========================================================= */

static void enumerate_matchings(
    int r,
    int next_edge,
    std::vector<int>& selected_edges,
    SPoly& total,
    u64& term_count,
    const std::vector<u64>& digits,
    int p
) {
    if (next_edge >= r) {
        /*
         * Build the term for this matching.
         */
        std::vector<bool> covered(
            static_cast<std::size_t>(r),
            false
        );

        for (int edge : selected_edges) {
            covered[
                static_cast<std::size_t>(edge - 1)
            ] = true;

            covered[
                static_cast<std::size_t>(edge)
            ] = true;
        }

        SPoly term{1};

        /*
         * Unmatched vertices contribute their diagonal
         * factors.
         */
        for (int vertex = 0; vertex < r; ++vertex) {
            if (covered[
                static_cast<std::size_t>(vertex)
            ]) {
                continue;
            }

            const UPoly factor =
                diagonal_factor(
                    digits,
                    p,
                    vertex
                );

            term =
                mul_s(
                    term,
                    signed_from_unsigned(factor)
                );
        }

        /*
         * Matched edges contribute -s_i.
         */
        for (int edge : selected_edges) {
            const UPoly factor =
                edge_factor(
                    digits,
                    p,
                    edge
                );

            SPoly signed_factor =
                signed_from_unsigned(
                    factor
                );

            signed_factor =
                neg_s(
                    signed_factor
                );

            term =
                mul_s(
                    term,
                    signed_factor
                );
        }

        total =
            add_s(
                total,
                term
            );

        ++term_count;
        return;
    }

    /*
     * Case 1:
     * Do not select this edge.
     */
    enumerate_matchings(
        r,
        next_edge + 1,
        selected_edges,
        total,
        term_count,
        digits,
        p
    );

    /*
     * Case 2:
     * Select this edge, so the next adjacent edge
     * cannot be selected.
     */
    selected_edges.push_back(
        next_edge
    );

    enumerate_matchings(
        r,
        next_edge + 2,
        selected_edges,
        total,
        term_count,
        digits,
        p
    );

    selected_edges.pop_back();
}

static SPoly matching_expansion(
    u64 L,
    int p,
    int r,
    u64& term_count
) {
    const auto digits =
        digits_of(L, p, r);

    SPoly result{0};

    std::vector<int> selected;

    term_count = 0;

    enumerate_matchings(
        r,
        1,
        selected,
        result,
        term_count,
        digits,
        p
    );

    trim_s(result);
    return result;
}

/* =========================================================
 * D_r A_r
 * ========================================================= */

static UPoly scaled_A(
    const std::vector<UPoly>& A,
    const std::vector<u64>& digits,
    int r
) {
    u128 scale = 1;

    if (r >= 2) {
        for (int i = 0; i <= r - 2; ++i) {
            scale *=
                static_cast<u128>(
                    digits[
                        static_cast<std::size_t>(i)
                    ]
                );
        }
    }

    return scalar_u(
        A[
            static_cast<std::size_t>(r)
        ],
        scale
    );
}

/* =========================================================
 * Local polynomial from history
 * ========================================================= */

static UPoly local_from_history(
    const std::vector<UPoly>& A,
    const std::vector<UPoly>& C,
    const std::vector<u64>& digits,
    int p,
    int r
) {
    UPoly result{0};

    for (int h = 0; h < r; ++h) {
        UPoly B =
            add_u(
                A[
                    static_cast<std::size_t>(h)
                ],
                C[
                    static_cast<std::size_t>(h)
                ]
            );

        const u128 choices =
            static_cast<u128>(p)
            -
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(h)
                ]
            )
            - 1;

        if (choices == 0) {
            continue;
        }

        const int shift =
            r - h;

        UPoly shifted(
            B.size() +
            static_cast<std::size_t>(shift),
            0
        );

        for (
            std::size_t k = 0;
            k < B.size();
            ++k
        ) {
            shifted[
                k +
                static_cast<std::size_t>(shift)
            ] = B[k];
        }

        result =
            add_u(
                result,
                scalar_u(
                    shifted,
                    choices
                )
            );
    }

    trim_u(result);
    return result;
}

/* =========================================================
 * One test case
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

    std::vector<UPoly> A;
    std::vector<UPoly> C;

    compute_history(
        L,
        p,
        r,
        A,
        C
    );

    const UPoly expected =
        scaled_A(
            A,
            digits,
            r
        );

    /*
     * --------------------------------------------------
     * Explicit matching expansion
     * --------------------------------------------------
     */

    if (matching_check) {
        ++stats.matching_cases;

        u64 term_count = 0;

        const SPoly expansion =
            matching_expansion(
                L,
                p,
                r,
                term_count
            );

        stats.matching_terms +=
            term_count;

        if (
            !equal_signed_unsigned(
                expansion,
                expected
            )
        ) {
            ok = false;
            ++stats.matching_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_MATCHING "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << " terms="
                    << term_count
                    << "\n";
            }
        }
    }

    /*
     * --------------------------------------------------
     * Coefficient check
     *
     * The previous comparison is polynomial-wide.
     * This explicitly checks every coefficient.
     * --------------------------------------------------
     */

    if (matching_check) {
        u64 dummy = 0;

        const SPoly expansion =
            matching_expansion(
                L,
                p,
                r,
                dummy
            );

        const std::size_t n =
            std::max(
                expansion.size(),
                expected.size()
            );

        for (
            std::size_t k = 0;
            k < n;
            ++k
        ) {
            const i128 x =
                k < expansion.size()
                    ? expansion[k]
                    : 0;

            const i128 y =
                k < expected.size()
                    ? static_cast<i128>(
                        expected[k]
                    )
                    : 0;

            if (x != y) {
                ok = false;
                ++stats.coefficient_failures;

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
                        k < expected.size()
                            ? expected[k]
                            : 0
                    );

                    std::cout << "\n";

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
        const UPoly local =
            local_from_history(
                A,
                C,
                digits,
                p,
                r
            );

        const UPoly exact =
            exact_local(
                L,
                p,
                r
            );

        if (!equal_u(local, exact)) {
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
        << "START EXPERIMENT 341\n"
        << "MATCHING EXPANSION OF THE KUMMER CONTINUANT\n"
        << "CAN THE DETERMINANT BE EXPANDED DIRECTLY OVER "
           "PATH MATCHINGS?\n\n";

    Stats stats;

    std::mt19937_64 rng(
        0x341341341ULL
    );

    /*
     * --------------------------------------------------
     * Exhaustive small cases
     *
     * Explicit matching enumeration is used here.
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
     * Random explicit-matching cases.
     *
     * Keep r modest because explicit matching enumeration
     * grows like Fibonacci(r).
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_MATCHING_CASES = 2200;

    for (
        u64 n = 0;
        n < RANDOM_MATCHING_CASES;
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
     * Larger cases:
     *
     * No explicit matching enumeration beyond r=16.
     * We still verify the original local structure.
     * --------------------------------------------------
     */

    constexpr u64 LARGE_CASES = 3500;

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

        /*
         * Matching expansion disabled here when r > 16.
         */
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
     * Targeted digit patterns
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
             * Equal-digit patterns.
             */
            for (
                u64 d = 0;
                d < static_cast<u64>(p);
                ++d
            ) {
                u64 value = 0;
                u64 place = 1;

                for (int i = 0; i < r; ++i) {
                    value +=
                        d * place;

                    place *=
                        static_cast<u64>(p);
                }

                cases.push_back(value);
            }

            /*
             * Alternating extremes.
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
            for (
                int pos = 0;
                pos < r;
                ++pos
            ) {
                if (p <= 2) {
                    continue;
                }

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

    std::cout
        << "============================\n"
        << "tested="
        << stats.tested
        << "\n"
        << "failed="
        << stats.failed
        << "\n"
        << "matching_failures="
        << stats.matching_failures
        << "\n"
        << "coefficient_failures="
        << stats.coefficient_failures
        << "\n"
        << "scaled_a_failures="
        << stats.scaled_a_failures
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
        << "FINISHED EXPERIMENT 341\n";

    return 0;
}
