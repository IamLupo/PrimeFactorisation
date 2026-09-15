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

    u64 continuant_failures = 0;
    u64 determinant_failures = 0;
    u64 exact_failures = 0;
    u64 total_failures = 0;

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

    char buf[64];
    int n = 0;

    while (x > 0) {
        buf[n++] =
            static_cast<char>('0' + (x % 10));
        x /= 10;
    }

    while (n > 0) {
        std::cout << buf[--n];
    }
}

/* =========================================================
 * Polynomial helpers
 * ========================================================= */

static void trim(Poly& p) {
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

static Poly scalar_mul(
    const Poly& a,
    u128 s
) {
    Poly r = a;

    for (u128& x : r) {
        x *= s;
    }

    trim(r);
    return r;
}

static Poly poly_add(
    const Poly& a,
    const Poly& b
) {
    Poly r(
        std::max(a.size(), b.size()),
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        r[i] += a[i];
    }

    for (std::size_t i = 0; i < b.size(); ++i) {
        r[i] += b[i];
    }

    trim(r);
    return r;
}

static Poly poly_mul(
    const Poly& a,
    const Poly& b
) {
    if (
        (a.size() == 1 && a[0] == 0) ||
        (b.size() == 1 && b[0] == 0)
    ) {
        return Poly{0};
    }

    Poly r(
        a.size() + b.size() - 1,
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            r[i + j] +=
                a[i] * b[j];
        }
    }

    trim(r);
    return r;
}

/*
 * Polynomial subtraction with a coefficient-wise
 * nonnegative result.
 */
static Poly poly_sub_nonnegative(
    const Poly& a,
    const Poly& b,
    bool& ok
) {
    const std::size_t n =
        std::max(a.size(), b.size());

    Poly r(n, 0);

    for (std::size_t i = 0; i < n; ++i) {
        const u128 x =
            i < a.size() ? a[i] : 0;

        const u128 y =
            i < b.size() ? b[i] : 0;

        if (x < y) {
            ok = false;
            return Poly{0};
        }

        r[i] = x - y;
    }

    trim(r);
    return r;
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

static SPoly signed_from_poly(
    const Poly& p
) {
    SPoly r(p.size());

    for (std::size_t i = 0; i < p.size(); ++i) {
        r[i] =
            static_cast<i128>(p[i]);
    }

    trim_signed(r);
    return r;
}

static SPoly signed_mul(
    const SPoly& a,
    const SPoly& b
) {
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

    trim_signed(r);
    return r;
}

static void signed_add(
    SPoly& dst,
    const SPoly& src,
    int sign
) {
    if (dst.size() < src.size()) {
        dst.resize(src.size(), 0);
    }

    for (std::size_t i = 0; i < src.size(); ++i) {
        if (sign > 0) {
            dst[i] += src[i];
        } else {
            dst[i] -= src[i];
        }
    }

    trim_signed(dst);
}

static bool signed_equals_unsigned(
    const SPoly& a,
    const Poly& b
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
    u64 x = 1;

    for (int i = 0; i < r; ++i) {
        x *= static_cast<u64>(p);
    }

    return x;
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
 * Kummer two-state polynomial system
 *
 * A_i = no-borrow state
 * C_i = borrow state
 *
 * A_{i+1} =
 *     (d_i+1) A_i
 *     + d_i C_i
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
                    static_cast<
                        std::size_t
                    >(i)
                ]
            );

        /*
         * A_{i+1}
         */
        A[i + 1] =
            scalar_mul(
                A[i],
                digit + 1
            );

        A[i + 1] =
            poly_add(
                A[i + 1],
                scalar_mul(
                    C[i],
                    digit
                )
            );

        /*
         * C_{i+1}
         */
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
            scalar_mul(
                zA,
                static_cast<u128>(p)
                - digit
                - 1
            );

        C[i + 1] =
            poly_add(
                C[i + 1],
                scalar_mul(
                    zC,
                    static_cast<u128>(p)
                    - digit
                )
            );
    }
}

/* =========================================================
 * Exact valuation
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

    trim(result);
    return result;
}

/* =========================================================
 * Scalar recurrence coefficients
 *
 * d_{i-1} A_{i+1}
 *
 * =
 * c_i A_i - e_i A_{i-1}
 *
 * where
 *
 * c_i =
 * d_{i-1}(d_i+1)
 * +d_i z(p-d_{i-1})
 *
 * e_i = d_i p z
 * ========================================================= */

static Poly make_c(
    u64 d_prev,
    u64 d_curr,
    int p
) {
    Poly result(2, 0);

    result[0] =
        static_cast<u128>(
            d_prev *
            (d_curr + 1)
        );

    result[1] =
        static_cast<u128>(d_curr) *
        (
            static_cast<u128>(p)
            - static_cast<u128>(d_prev)
        );

    trim(result);
    return result;
}

static Poly make_e(
    u64 d_curr,
    int p
) {
    Poly result{
        0,
        static_cast<u128>(d_curr) *
        static_cast<u128>(p)
    };

    trim(result);
    return result;
}

/* =========================================================
 * D_r A_r
 *
 * D_r =
 * product_{j=0}^{r-2} d_j
 *
 * Hence:
 *
 * K_r = D_r A_r
 * ========================================================= */

static Poly scaled_A(
    const std::vector<Poly>& A,
    const std::vector<u64>& d,
    int r
) {
    u128 scale = 1;

    if (r >= 2) {
        for (int j = 0; j <= r - 2; ++j) {
            scale *=
                static_cast<u128>(
                    d[
                        static_cast<
                            std::size_t
                        >(j)
                    ]
                );
        }
    }

    return scalar_mul(
        A[
            static_cast<
                std::size_t
            >(r)
        ],
        scale
    );
}

/* =========================================================
 * Continuant
 *
 * K_0 = 1
 * K_1 = d_0 + 1
 *
 * For i >= 1:
 *
 * K_{i+1}
 * =
 * c_i K_i
 * -
 * e_i d_{i-2} K_{i-1}
 *
 * For i=1, the second term is simply e_1.
 * ========================================================= */

static Poly continuant(
    u64 L,
    int p,
    int r
) {
    const auto d =
        digits_of(L, p, r);

    if (r == 0) {
        return Poly{1};
    }

    /*
     * K_1 = A_1 = d_0 + 1.
     */
    Poly K_prev{
        static_cast<u128>(
            d[0] + 1
        )
    };

    if (r == 1) {
        return K_prev;
    }

    /*
     * K_2 =
     * c_1 K_1 - e_1
     */
    const Poly c1 =
        make_c(
            d[0],
            d[1],
            p
        );

    const Poly e1 =
        make_e(
            d[1],
            p
        );

    const Poly first =
        poly_mul(
            c1,
            K_prev
        );

    bool ok = true;

    Poly K_curr =
        poly_sub_nonnegative(
            first,
            e1,
            ok
        );

    if (!ok) {
        return Poly{0};
    }

    /*
     * K_{i+1}, i >= 2
     */
    for (int i = 2; i < r; ++i) {
        const Poly ci =
            make_c(
                d[
                    static_cast<
                        std::size_t
                    >(i - 1)
                ],
                d[
                    static_cast<
                        std::size_t
                    >(i)
                ],
                p
            );

        const Poly ei =
            make_e(
                d[
                    static_cast<
                        std::size_t
                    >(i)
                ],
                p
            );

        const Poly second_coefficient =
            scalar_mul(
                ei,
                static_cast<u128>(
                    d[
                        static_cast<
                            std::size_t
                        >(i - 2)
                    ]
                )
            );

        const Poly first_term =
            poly_mul(
                ci,
                K_curr
            );

        const Poly second_term =
            poly_mul(
                second_coefficient,
                K_prev
            );

        bool local_ok = true;

        Poly next =
            poly_sub_nonnegative(
                first_term,
                second_term,
                local_ok
            );

        if (!local_ok) {
            return Poly{0};
        }

        K_prev =
            std::move(K_curr);

        K_curr =
            std::move(next);
    }

    return K_curr;
}

/* =========================================================
 * Build the tridiagonal matrix
 *
 * T_r =
 *
 * [ d0+1      e1          0        ... ]
 * [ 1         c1          e2*d0    ... ]
 * [ 0         1           c2       ... ]
 * [ ...                           ... ]
 *
 * Its determinant is K_r.
 * ========================================================= */

static std::vector<std::vector<Poly>>
build_matrix(
    u64 L,
    int p,
    int r
) {
    const auto d =
        digits_of(L, p, r);

    std::vector<std::vector<Poly>> M(
        static_cast<std::size_t>(r),
        std::vector<Poly>(
            static_cast<std::size_t>(r),
            Poly{0}
        )
    );

    if (r == 0) {
        return M;
    }

    M[0][0] = {
        static_cast<u128>(
            d[0] + 1
        )
    };

    for (int i = 1; i < r; ++i) {
        /*
         * Diagonal: c_i.
         */
        M[i][i] =
            make_c(
                d[
                    static_cast<
                        std::size_t
                    >(i - 1)
                ],
                d[
                    static_cast<
                        std::size_t
                    >(i)
                ],
                p
            );

        /*
         * Subdiagonal = 1.
         */
        M[i][i - 1] = {1};

        /*
         * Superdiagonal:
         *
         * i=1:
         *   e_1
         *
         * i>=2:
         *   e_i d_{i-1?}
         *
         * Careful: the continuant recurrence is
         *
         * K_{i+1}
         * =
         * c_i K_i
         * -
         * e_i d_{i-2} K_{i-1}.
         *
         * Therefore the superdiagonal for row i-1
         * is e_i d_{i-2}.
         */
        const Poly e =
            make_e(
                d[
                    static_cast<
                        std::size_t
                    >(i)
                ],
                p
            );

        if (i == 1) {
            M[0][1] = e;
        } else {
            M[i - 1][i] =
                scalar_mul(
                    e,
                    static_cast<u128>(
                        d[
                            static_cast<
                                std::size_t
                            >(i - 2)
                        ]
                    )
                );
        }
    }

    return M;
}

/* =========================================================
 * Independent determinant by permutations
 *
 * Only used for r <= 7.
 * ========================================================= */

static SPoly direct_determinant(
    const std::vector<std::vector<Poly>>& M
) {
    const int n =
        static_cast<int>(M.size());

    if (n == 0) {
        return SPoly{1};
    }

    std::vector<int> permutation(
        static_cast<std::size_t>(n)
    );

    for (int i = 0; i < n; ++i) {
        permutation[
            static_cast<
                std::size_t
            >(i)
        ] = i;
    }

    SPoly result{0};

    do {
        int inversions = 0;

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (
                    permutation[
                        static_cast<
                            std::size_t
                        >(i)
                    ] >
                    permutation[
                        static_cast<
                            std::size_t
                        >(j)
                    ]
                ) {
                    ++inversions;
                }
            }
        }

        const int sign =
            (inversions & 1)
                ? -1
                : 1;

        SPoly term{1};

        bool zero = false;

        for (int i = 0; i < n; ++i) {
            const Poly& entry =
                M[i][
                    permutation[
                        static_cast<
                            std::size_t
                        >(i)
                    ]
                ];

            if (
                entry.size() == 1 &&
                entry[0] == 0
            ) {
                zero = true;
                break;
            }

            term =
                signed_mul(
                    term,
                    signed_from_poly(entry)
                );
        }

        if (!zero) {
            signed_add(
                result,
                term,
                sign
            );
        }

    } while (
        std::next_permutation(
            permutation.begin(),
            permutation.end()
        )
    );

    trim_signed(result);
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
            A[
                static_cast<
                    std::size_t
                >(h)
            ];

        B =
            poly_add(
                B,
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

        const int shift =
            r - h;

        Poly shifted(
            B.size() +
            static_cast<
                std::size_t
            >(shift),
            0
        );

        for (
            std::size_t k = 0;
            k < B.size();
            ++k
        ) {
            shifted[
                k +
                static_cast<
                    std::size_t
                >(shift)
            ] = B[k];
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

    trim(result);
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
    bool determinant_check,
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
     * --------------------------------------------------
     * Continuant identity
     *
     * K_r = D_r A_r
     * --------------------------------------------------
     */

    const Poly expected_continuant =
        scaled_A(
            A,
            digits,
            r
        );

    const Poly actual_continuant =
        continuant(
            L,
            p,
            r
        );

    if (
        !equal_poly(
            expected_continuant,
            actual_continuant
        )
    ) {
        ok = false;
        ++stats.continuant_failures;

        if (print_failure) {
            std::cout
                << "FAIL_CONTINUANT "
                << "L=" << L
                << " p=" << p
                << " r=" << r
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Independent determinant
     * --------------------------------------------------
     */

    if (determinant_check) {
        const auto matrix =
            build_matrix(
                L,
                p,
                r
            );

        const SPoly determinant =
            direct_determinant(matrix);

        if (
            !signed_equals_unsigned(
                determinant,
                expected_continuant
            )
        ) {
            ok = false;
            ++stats.determinant_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_DETERMINANT "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << "\n";
            }
        }
    }

    /*
     * --------------------------------------------------
     * Exact local construction
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

        const u64 expected =
            prime_power(p, r)
            - L
            - 1;

        if (actual != expected) {
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
        << "START EXPERIMENT 340\n"
        << "TRIDIAGONAL CONTINUANT / DETERMINANT FORM\n"
        << "CAN THE SCALAR RECURRENCE BE WRITTEN AS A "
           "DIVISION-FREE DETERMINANT?\n\n";

    Stats stats;

    std::mt19937_64 rng(
        0x340340340ULL
    );

    /*
     * --------------------------------------------------
     * Exhaustive small cases
     * --------------------------------------------------
     */

    constexpr u64 SMALL_LIMIT = 1000;

    for (int p : PRIMES) {
        for (int r = 1; r <= 12; ++r) {
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
                    r <= 7,
                    stats.failed < 10,
                    stats
                );
            }
        }
    }

    /*
     * --------------------------------------------------
     * Random larger cases
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_CASES = 7000;

    for (u64 n = 0; n < RANDOM_CASES; ++n) {
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

        verify_case(
            L,
            p,
            r,
            power <= 50000,
            r <= 7,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Targeted cases
     * --------------------------------------------------
     */

    for (int p : PRIMES) {
        for (int r = 1; r <= 20; ++r) {
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
             * Equal-digit values.
             */
            for (
                u64 d = 0;
                d < static_cast<u64>(p);
                ++d
            ) {
                u64 value = 0;
                u64 place = 1;

                for (int i = 0; i < r; ++i) {
                    value += d * place;
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
                    const u64 d =
                        ((i & 1) == mode)
                            ? static_cast<u64>(p - 1)
                            : 0;

                    value += d * place;

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
                    r <= 7,
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
        << "continuant_failures="
        << stats.continuant_failures
        << "\n"
        << "determinant_failures="
        << stats.determinant_failures
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
        << "FINISHED EXPERIMENT 340\n";

    return 0;
}