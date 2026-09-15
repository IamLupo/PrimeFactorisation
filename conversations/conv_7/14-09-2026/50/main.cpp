#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

using Poly = std::vector<u128>;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 path_failures = 0;
    u64 recurrence_failures = 0;
    u64 constant_digit_failures = 0;
    u64 exact_failures = 0;
    u64 total_failures = 0;

    u64 targeted = 0;
};

/* =========================================================
 * u128 printing
 * ========================================================= */

static void print_u128(
    u128 value
) {
    if (value == 0) {
        std::cout << '0';
        return;
    }

    char buffer[64];
    int pos = 0;

    while (value > 0) {
        const unsigned digit =
            static_cast<unsigned>(
                value % 10
            );

        buffer[pos++] =
            static_cast<char>(
                '0' + digit
            );

        value /= 10;
    }

    while (pos > 0) {
        std::cout
            << buffer[--pos];
    }
}

/* =========================================================
 * Polynomial helpers
 * ========================================================= */

static void trim(
    Poly& p
) {
    while (
        p.size() > 1 &&
        p.back() == 0
    ) {
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
        std::max(
            a.size(),
            b.size()
        );

    for (
        std::size_t i = 0;
        i < n;
        ++i
    ) {
        const u128 x =
            i < a.size()
                ? a[i]
                : 0;

        const u128 y =
            i < b.size()
                ? b[i]
                : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
}

static Poly scalar_mul(
    const Poly& p,
    u128 scale
) {
    Poly result = p;

    for (u128& x : result) {
        x *= scale;
    }

    trim(result);
    return result;
}

static Poly shift_z(
    const Poly& p
) {
    Poly result(
        p.size() + 1,
        0
    );

    for (
        std::size_t i = 0;
        i < p.size();
        ++i
    ) {
        result[i + 1] = p[i];
    }

    trim(result);
    return result;
}

static void add_poly(
    Poly& dst,
    const Poly& src
) {
    if (
        dst.size() <
        src.size()
    ) {
        dst.resize(
            src.size(),
            0
        );
    }

    for (
        std::size_t i = 0;
        i < src.size();
        ++i
    ) {
        dst[i] += src[i];
    }

    trim(dst);
}

static void add_scaled(
    Poly& dst,
    const Poly& src,
    u128 scale
) {
    if (scale == 0) {
        return;
    }

    if (
        dst.size() <
        src.size()
    ) {
        dst.resize(
            src.size(),
            0
        );
    }

    for (
        std::size_t i = 0;
        i < src.size();
        ++i
    ) {
        dst[i] +=
            src[i] * scale;
    }

    trim(dst);
}

/* =========================================================
 * p-adic helpers
 * ========================================================= */

static u64 prime_power(
    int p,
    int r
) {
    u64 result = 1;

    for (
        int i = 0;
        i < r;
        ++i
    ) {
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

/* =========================================================
 * Digits
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

    for (
        int i = 0;
        i < length;
        ++i
    ) {
        digits[
            static_cast<std::size_t>(i)
        ] =
            n % prime;

        n /= prime;
    }

    return digits;
}

/* =========================================================
 * Two-state history
 *
 * A_i = no-borrow polynomial
 * C_i = borrow polynomial
 * ========================================================= */

static void compute_history(
    u64 L,
    int p,
    int r,
    std::vector<Poly>& A,
    std::vector<Poly>& C
) {
    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    A.resize(
        static_cast<std::size_t>(r + 1)
    );

    C.resize(
        static_cast<std::size_t>(r + 1)
    );

    A[0] = {1};
    C[0] = {0};

    for (
        int i = 0;
        i < r;
        ++i
    ) {
        const u128 d =
            static_cast<u128>(
                digits[
                    static_cast<
                        std::size_t
                    >(i)
                ]
            );

        /*
         * A_{i+1}
         *
         * = (d+1) A_i
         *   + d C_i
         */
        A[i + 1] =
            scalar_mul(
                A[i],
                d + 1
            );

        add_scaled(
            A[i + 1],
            C[i],
            d
        );

        /*
         * C_{i+1}
         *
         * = z(p-d-1) A_i
         *   + z(p-d) C_i
         */
        C[i + 1] =
            scalar_mul(
                shift_z(
                    A[i]
                ),
                static_cast<u128>(p)
                - d
                - 1
            );

        add_scaled(
            C[i + 1],
            shift_z(
                C[i]
            ),
            static_cast<u128>(p)
            - d
        );
    }
}

/* =========================================================
 * Explicit path expansion
 *
 * Returns A_h and C_h separately.
 *
 * Only used for small h.
 * ========================================================= */

static void explicit_paths(
    const std::vector<u64>& digits,
    int p,
    int h,
    Poly& A,
    Poly& C
) {
    A.assign(
        static_cast<std::size_t>(h + 1),
        0
    );

    C.assign(
        static_cast<std::size_t>(h + 1),
        0
    );

    if (h == 0) {
        A[0] = 1;
        return;
    }

    const u64 path_count =
        1ULL << h;

    for (
        u64 mask = 0;
        mask < path_count;
        ++mask
    ) {
        u128 weight = 1;

        int borrow = 0;
        int valuation = 0;

        bool valid = true;

        for (
            int i = 0;
            i < h;
            ++i
        ) {
            const int next_borrow =
                static_cast<int>(
                    (mask >> i) & 1ULL
                );

            const u128 d =
                static_cast<u128>(
                    digits[
                        static_cast<
                            std::size_t
                        >(i)
                    ]
                );

            u128 ways = 0;

            if (
                borrow == 0 &&
                next_borrow == 0
            ) {
                ways = d + 1;
            } else if (
                borrow == 0 &&
                next_borrow == 1
            ) {
                ways =
                    static_cast<u128>(p)
                    - d
                    - 1;
            } else if (
                borrow == 1 &&
                next_borrow == 0
            ) {
                ways = d;
            } else {
                ways =
                    static_cast<u128>(p)
                    - d;
            }

            if (ways == 0) {
                valid = false;
                break;
            }

            weight *= ways;

            valuation +=
                next_borrow;

            borrow =
                next_borrow;
        }

        if (!valid) {
            continue;
        }

        if (borrow == 0) {
            A[
                static_cast<
                    std::size_t
                >(valuation)
            ] += weight;
        } else {
            C[
                static_cast<
                    std::size_t
                >(valuation)
            ] += weight;
        }
    }

    trim(A);
    trim(C);
}

/* =========================================================
 * Local polynomial reconstructed from A_h + C_h
 * ========================================================= */

static Poly local_from_history(
    u64 L,
    int p,
    int r,
    const std::vector<Poly>& A,
    const std::vector<Poly>& C
) {
    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    Poly result{0};

    for (
        int h = 0;
        h < r;
        ++h
    ) {
        const u128 d =
            static_cast<u128>(
                digits[
                    static_cast<
                        std::size_t
                    >(h)
                ]
            );

        const u128 choices =
            static_cast<u128>(p)
            - d
            - 1;

        if (choices == 0) {
            continue;
        }

        /*
         * B_h = A_h + C_h.
         */
        Poly B =
            A[
                static_cast<
                    std::size_t
                >(h)
            ];

        add_poly(
            B,
            C[
                static_cast<
                    std::size_t
                >(h)
            ]
        );

        /*
         * First strict digit creates one borrow and
         * all higher positions retain that borrow.
         *
         * Fixed valuation contribution = r-h.
         */
        const int fixed =
            r - h;

        Poly shifted(
            B.size() +
                static_cast<
                    std::size_t
                >(fixed),
            0
        );

        for (
            std::size_t j = 0;
            j < B.size();
            ++j
        ) {
            shifted[
                j +
                static_cast<
                    std::size_t
                >(fixed)
            ] =
                B[j];
        }

        add_scaled(
            result,
            shifted,
            choices
        );
    }

    trim(result);
    return result;
}

/* =========================================================
 * Exact local polynomial
 * ========================================================= */

static Poly exact_local(
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

    Poly result{0};

    for (
        u64 y = L + 1;
        y < power;
        ++y
    ) {
        const u64 v =
            binomial_valuation(
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
 * Correct scalar recurrence
 *
 * d_{i-1} A_{i+1}
 *
 * =
 * [d_{i-1}(d_i+1)
 *  +d_i(p-d_{i-1})z] A_i
 *
 * -d_i p z A_{i-1}.
 *
 * Use signed __int128 during subtraction.
 * ========================================================= */

static bool verify_scalar_recurrence(
    u64 L,
    int p,
    int r,
    const std::vector<Poly>& A,
    bool print_failure
) {
    if (r < 2) {
        return true;
    }

    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    for (
        int i = 1;
        i < r;
        ++i
    ) {
        const u128 d_prev =
            static_cast<u128>(
                digits[
                    static_cast<
                        std::size_t
                    >(i - 1)
                ]
            );

        const u128 d_curr =
            static_cast<u128>(
                digits[
                    static_cast<
                        std::size_t
                    >(i)
                ]
            );

        /*
         * Left:
         *
         * d_{i-1} A_{i+1}
         */
        Poly left =
            scalar_mul(
                A[
                    static_cast<
                        std::size_t
                    >(i + 1)
                ],
                d_prev
            );

        /*
         * Right:
         *
         * d_prev(d_curr+1) A_i
         *
         * + d_curr(p-d_prev) z A_i
         *
         * - d_curr*p*z*A_{i-1}
         */
        Poly positive{
            0
        };

        add_scaled(
            positive,
            A[
                static_cast<
                    std::size_t
                >(i)
            ],
            d_prev *
            (d_curr + 1)
        );

        add_scaled(
            positive,
            shift_z(
                A[
                    static_cast<
                        std::size_t
                    >(i)
                ]
            ),
            d_curr *
            (
                static_cast<u128>(p)
                - d_prev
            )
        );

        const Poly negative =
            scalar_mul(
                shift_z(
                    A[
                        static_cast<
                            std::size_t
                        >(i - 1)
                    ]
                ),
                d_curr *
                static_cast<u128>(p)
            );

        const std::size_t n =
            std::max(
                positive.size(),
                negative.size()
            );

        std::vector<i128> right(
            n,
            0
        );

        std::vector<i128> lhs(
            n,
            0
        );

        for (
            std::size_t k = 0;
            k < left.size();
            ++k
        ) {
            lhs[k] =
                static_cast<i128>(
                    left[k]
                );
        }

        for (
            std::size_t k = 0;
            k < positive.size();
            ++k
        ) {
            right[k] +=
                static_cast<i128>(
                    positive[k]
                );
        }

        for (
            std::size_t k = 0;
            k < negative.size();
            ++k
        ) {
            right[k] -=
                static_cast<i128>(
                    negative[k]
                );
        }

        for (
            std::size_t k = 0;
            k < n;
            ++k
        ) {
            if (
                lhs[k] != right[k]
            ) {
                if (print_failure) {
                    std::cout
                        << "FAIL_SCALAR "
                        << "L=" << L
                        << " p=" << p
                        << " r=" << r
                        << " i=" << i
                        << " k=" << k
                        << "\n"
                        << "lhs=";

                    print_u128(
                        left.size() > k
                            ? left[k]
                            : 0
                    );

                    std::cout
                        << "\n";
                }

                return false;
            }
        }
    }

    return true;
}

/* =========================================================
 * Constant-digit recurrence
 *
 * For every digit d:
 *
 * d A_{i+1}
 * =
 * d[(d+1)+(p-d)z]A_i
 * -d p z A_{i-1}.
 * ========================================================= */

static bool verify_constant_digit_recurrence(
    int d,
    int p,
    int r,
    bool print_failure
) {
    if (r < 2) {
        return true;
    }

    u64 L = 0;
    u64 power = 1;

    for (
        int i = 0;
        i < r;
        ++i
    ) {
        L +=
            static_cast<u64>(d) *
            power;

        power *=
            static_cast<u64>(p);
    }

    std::vector<Poly> A;
    std::vector<Poly> C;

    compute_history(
        L,
        p,
        r,
        A,
        C
    );

    for (
        int i = 1;
        i < r;
        ++i
    ) {
        const u128 D =
            static_cast<u128>(d);

        Poly left =
            scalar_mul(
                A[
                    static_cast<
                        std::size_t
                    >(i + 1)
                ],
                D
            );

        Poly positive{
            0
        };

        add_scaled(
            positive,
            A[
                static_cast<
                    std::size_t
                >(i)
            ],
            D * (D + 1)
        );

        add_scaled(
            positive,
            shift_z(
                A[
                    static_cast<
                        std::size_t
                    >(i)
                ]
            ),
            D *
            (
                static_cast<u128>(p)
                - D
            )
        );

        const Poly negative =
            scalar_mul(
                shift_z(
                    A[
                        static_cast<
                            std::size_t
                        >(i - 1)
                    ]
                ),
                D *
                static_cast<u128>(p)
            );

        const std::size_t n =
            std::max(
                positive.size(),
                negative.size()
            );

        for (
            std::size_t k = 0;
            k < n;
            ++k
        ) {
            const i128 lhs =
                k < left.size()
                    ? static_cast<i128>(
                        left[k]
                    )
                    : 0;

            const i128 pos =
                k < positive.size()
                    ? static_cast<i128>(
                        positive[k]
                    )
                    : 0;

            const i128 neg =
                k < negative.size()
                    ? static_cast<i128>(
                        negative[k]
                    )
                    : 0;

            if (
                lhs !=
                pos - neg
            ) {
                if (print_failure) {
                    std::cout
                        << "FAIL_CONSTANT "
                        << "p=" << p
                        << " d=" << d
                        << " r=" << r
                        << " i=" << i
                        << "\n";
                }

                return false;
            }
        }
    }

    return true;
}

/* =========================================================
 * Build L from digits
 * ========================================================= */

static bool build_from_digits(
    const std::vector<u64>& digits,
    int p,
    u64& L
) {
    L = 0;

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
            L >
            UINT64_MAX - term
        ) {
            return false;
        }

        L += term;

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
 * Verify one local case
 * ========================================================= */

static bool verify_case(
    u64 L,
    int p,
    int r,
    bool exact_check,
    bool path_check,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

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
     * Explicit paths
     * --------------------------------------------------
     */

    if (path_check) {
        const auto digits =
            digits_of(
                L,
                p,
                r
            );

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
            Poly path_A;
            Poly path_C;

            explicit_paths(
                digits,
                p,
                h,
                path_A,
                path_C
            );

            if (
                !equal_poly(
                    A[
                        static_cast<
                            std::size_t
                        >(h)
                    ],
                    path_A
                ) ||
                !equal_poly(
                    C[
                        static_cast<
                            std::size_t
                        >(h)
                    ],
                    path_C
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
     * Scalar recurrence
     * --------------------------------------------------
     */

    if (
        !verify_scalar_recurrence(
            L,
            p,
            r,
            A,
            print_failure
        )
    ) {
        ok = false;
        ++stats.recurrence_failures;
    }

    /*
     * --------------------------------------------------
     * Local polynomial
     * --------------------------------------------------
     */

    const Poly local =
        local_from_history(
            L,
            p,
            r,
            A,
            C
        );

    /*
     * Gap length.
     */
    const u64 expected_length =
        prime_power(
            p,
            r
        ) -
        L -
        1;

    u64 local_total = 0;

    for (u128 x : local) {
        /*
         * Current parameter sizes keep this within u64
         * for the tested cases.
         */
        local_total +=
            static_cast<u64>(x);
    }

    if (
        local_total !=
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
                << local_total
                << " expected="
                << expected_length
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Exact small-block validation
     * --------------------------------------------------
     */

    if (exact_check) {
        const Poly exact =
            exact_local(
                L,
                p,
                r
            );

        if (
            !equal_poly(
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

    /*
     * --------------------------------------------------
     * Constant-digit recurrence
     * --------------------------------------------------
     */

    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    bool constant =
        true;

    for (
        int i = 1;
        i < r;
        ++i
    ) {
        if (
            digits[
                static_cast<
                    std::size_t
                >(i)
            ] !=
            digits[0]
        ) {
            constant = false;
            break;
        }
    }

    if (constant) {
        if (
            !verify_constant_digit_recurrence(
                static_cast<int>(
                    digits[0]
                ),
                p,
                r,
                print_failure
            )
        ) {
            ok = false;
            ++stats.constant_digit_failures;
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
        << "START EXPERIMENT 339\n"
        << "SCALAR SECOND-ORDER KUMMER RECURRENCE\n"
        << "CAN THE TWO-STATE BORROW SYSTEM BE ELIMINATED?\n\n";

    Stats stats;

    std::mt19937_64 rng(
        0x339339339ULL
    );

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small blocks.
     * --------------------------------------------------
     */

    constexpr u64 EXHAUSTIVE_LIMIT = 1200;

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
     * Phase 2:
     * Random larger cases.
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_CASES = 6000;

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
     * Phase 3:
     * Targeted equal-digit and structured cases.
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

            /*
             * Equal-digit values.
             */
            for (
                u64 d = 0;
                d < static_cast<u64>(p);
                ++d
            ) {
                std::vector<u64> digits(
                    static_cast<
                        std::size_t
                    >(r),
                    d
                );

                u64 value = 0;

                if (
                    build_from_digits(
                        digits,
                        p,
                        value
                    )
                ) {
                    candidates.push_back(
                        value
                    );
                }
            }

            /*
             * Alternating / increasing patterns.
             */
            for (int pattern = 0;
                 pattern < 3;
                 ++pattern) {

                std::vector<u64> digits(
                    static_cast<
                        std::size_t
                    >(r),
                    0
                );

                for (
                    int i = 0;
                    i < r;
                    ++i
                ) {
                    if (pattern == 0) {
                        if (i & 1) {
                            digits[i] =
                                static_cast<u64>(
                                    p - 1
                                );
                        }
                    } else if (
                        pattern == 1
                    ) {
                        if (!(i & 1)) {
                            digits[i] =
                                static_cast<u64>(
                                    p - 1
                                );
                        }
                    } else {
                        digits[i] =
                            static_cast<u64>(
                                i % p
                            );
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

            /*
             * Deterministic boundaries.
             */
            candidates.push_back(0);

            if (power > 1) {
                candidates.push_back(
                    power - 1
                );
            }

            if (power > 3) {
                candidates.push_back(1);
                candidates.push_back(2);
                candidates.push_back(
                    power / 2
                );
            }

            for (u64 L : candidates) {
                ++stats.targeted;
                ++stats.tested;

                verify_case(
                    L,
                    p,
                    r,
                    power <= EXHAUSTIVE_LIMIT,
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
        << "recurrence_failures="
        << stats.recurrence_failures
        << "\n"
        << "constant_digit_failures="
        << stats.constant_digit_failures
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
        << "FINISHED EXPERIMENT 339\n";

    return 0;
}