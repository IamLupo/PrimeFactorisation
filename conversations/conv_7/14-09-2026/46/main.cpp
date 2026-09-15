#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <vector>

using u64 = std::uint64_t;

/*
 * Polynomial:
 *
 * coeff[k] = coefficient of z^k
 */
using Poly = std::vector<u64>;

static const std::array<int, 6> PRIMES = {
    2, 3, 5, 7, 11, 13
};

/*
 * State:
 *
 * borrow = 0/1
 *
 * relation:
 *   0 = y == L
 *   1 = y >  L
 *   2 = y <  L
 *
 * state = borrow*3 + relation
 *
 * 0..2  -> borrow 0
 * 3..5  -> borrow 1
 */

constexpr int STATES = 6;

struct Matrix {
    Poly cell[STATES][STATES];
};

struct Stats {
    u64 tested = 0;
    u64 failed = 0;

    u64 exact_failures = 0;
    u64 dp_failures = 0;
    u64 matrix_failures = 0;
    u64 split_failures = 0;
    u64 total_failures = 0;
    u64 degree_failures = 0;

    u64 targeted = 0;
};

/* =========================================================
 * Polynomial helpers
 * ========================================================= */

static Poly poly_zero() {
    return Poly{0};
}

static Poly poly_one() {
    return Poly{1};
}

static void trim_poly(
    Poly& p
) {
    while (
        p.size() > 1 &&
        p.back() == 0
    ) {
        p.pop_back();
    }
}

static bool poly_equal(
    const Poly& a,
    const Poly& b
) {
    std::size_t n =
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

static void poly_add_inplace(
    Poly& dst,
    const Poly& src
) {
    if (dst.size() < src.size()) {
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

    trim_poly(dst);
}

static Poly poly_shift(
    const Poly& p,
    int shift
) {
    if (
        p.empty() ||
        (p.size() == 1 && p[0] == 0)
    ) {
        return Poly{0};
    }

    Poly result(
        p.size() +
            static_cast<std::size_t>(shift),
        0
    );

    for (
        std::size_t i = 0;
        i < p.size();
        ++i
    ) {
        result[
            i +
            static_cast<std::size_t>(shift)
        ] = p[i];
    }

    return result;
}

static Poly poly_mul(
    const Poly& a,
    const Poly& b
) {
    if (
        a.empty() ||
        b.empty()
    ) {
        return Poly{0};
    }

    if (
        (a.size() == 1 && a[0] == 0) ||
        (b.size() == 1 && b[0] == 0)
    ) {
        return Poly{0};
    }

    Poly result(
        a.size() +
            b.size() -
            1,
        0
    );

    for (
        std::size_t i = 0;
        i < a.size();
        ++i
    ) {
        if (a[i] == 0) {
            continue;
        }

        for (
            std::size_t j = 0;
            j < b.size();
            ++j
        ) {
            if (b[j] == 0) {
                continue;
            }

            result[i + j] +=
                a[i] * b[j];
        }
    }

    trim_poly(result);
    return result;
}

static std::ostream&
operator<<(
    std::ostream& out,
    const Poly& p
) {
    out << "[";

    for (
        std::size_t i = 0;
        i < p.size();
        ++i
    ) {
        if (i != 0) {
            out << ",";
        }

        out << p[i];
    }

    out << "]";
    return out;
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
 * p^r
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
 * Exact local generating polynomial
 *
 * K_{r,L}(z)
 *
 * = sum_{L<y<p^r}
 *     z^{v_p(C(p^r+L,y))}
 * ========================================================= */

static Poly exact_local_polynomial(
    u64 L,
    int p,
    int r
) {
    Poly result{0};

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
                static_cast<std::size_t>(v + 1),
                0
            );
        }

        ++result[v];
    }

    trim_poly(result);
    return result;
}

/* =========================================================
 * Add one digit to a state
 * ========================================================= */

static void add_digit_transition(
    Poly& destination,
    const Poly& source,
    int source_borrow,
    int source_relation,
    int new_borrow,
    int new_relation,
    int borrow_added
) {
    (void)new_borrow;
    (void)new_relation;

    Poly shifted =
        poly_shift(
            source,
            borrow_added
        );

    poly_add_inplace(
        destination,
        shifted
    );

    (void)source_borrow;
    (void)source_relation;
}

/* =========================================================
 * Construct the single-digit transfer matrix
 *
 * For a fixed L_i = digit:
 *
 * each possible y_i in [0,p-1] creates one transition.
 *
 * The polynomial weight is z if the transition creates
 * a borrow, otherwise 1.
 * ========================================================= */

static Matrix digit_matrix(
    u64 Ldigit,
    int p
) {
    Matrix M{};

    for (int s = 0;
         s < STATES;
         ++s) {

        const int borrow =
            s / 3;

        const int relation =
            s % 3;

        for (
            u64 ydigit = 0;
            ydigit < static_cast<u64>(p);
            ++ydigit
        ) {
            /*
             * Numerical relation:
             *
             * Higher digits dominate lower digits.
             */
            int new_relation =
                relation;

            if (ydigit > Ldigit) {
                new_relation = 1;
            } else if (ydigit < Ldigit) {
                new_relation = 2;
            }

            /*
             * Borrow in subtraction L-y.
             */
            int new_borrow = 0;
            int borrow_added = 0;

            if (
                Ldigit <
                ydigit +
                static_cast<u64>(borrow)
            ) {
                new_borrow = 1;
                borrow_added = 1;
            }

            const int ns =
                new_borrow * 3 +
                new_relation;

            if (
                M.cell[s][ns].empty()
            ) {
                M.cell[s][ns] =
                    poly_zero();
            }

            if (
                borrow_added == 0
            ) {
                poly_add_inplace(
                    M.cell[s][ns],
                    Poly{1}
                );
            } else {
                poly_add_inplace(
                    M.cell[s][ns],
                    Poly{0, 1}
                );
            }
        }
    }

    return M;
}

/* =========================================================
 * Identity matrix
 * ========================================================= */

static Matrix matrix_identity() {
    Matrix I{};

    for (int i = 0;
         i < STATES;
         ++i) {
        I.cell[i][i] =
            Poly{1};
    }

    return I;
}

/* =========================================================
 * Matrix multiplication
 * ========================================================= */

static Matrix matrix_multiply(
    const Matrix& A,
    const Matrix& B
) {
    Matrix C{};

    for (int i = 0;
         i < STATES;
         ++i) {

        for (int k = 0;
             k < STATES;
             ++k) {

            if (
                A.cell[i][k].empty()
            ) {
                continue;
            }

            for (int j = 0;
                 j < STATES;
                 ++j) {

                if (
                    B.cell[k][j].empty()
                ) {
                    continue;
                }

                const Poly product =
                    poly_mul(
                        A.cell[i][k],
                        B.cell[k][j]
                    );

                poly_add_inplace(
                    C.cell[i][j],
                    product
                );
            }
        }
    }

    return C;
}

/* =========================================================
 * Extract local polynomial from matrix product
 *
 * Initial state:
 *
 *   borrow=0
 *   relation=equal
 *
 * => state 0.
 *
 * Accept relation y>L.
 *
 * Both final borrow states are accepted because the leading
 * digit 1 of p^r+L absorbs an outstanding borrow.
 * ========================================================= */

static Poly extract_local_polynomial(
    const Matrix& M
) {
    Poly result{0};

    const int start_state =
        0;

    for (
        int borrow = 0;
        borrow <= 1;
        ++borrow
    ) {
        const int state =
            borrow * 3 + 1;

        poly_add_inplace(
            result,
            M.cell[
                start_state
            ][state]
        );
    }

    trim_poly(result);
    return result;
}

/* =========================================================
 * Transfer polynomial from all L digits
 *
 * Digits are processed low -> high.
 * ========================================================= */

static Poly transfer_polynomial(
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

    Matrix product =
        matrix_identity();

    /*
     * The state transition for digit 0 applies first,
     * then digit 1, etc.
     */
    for (int i = 0;
         i < r;
         ++i) {

        const Matrix D =
            digit_matrix(
                digits[i],
                p
            );

        product =
            matrix_multiply(
                product,
                D
            );
    }

    return extract_local_polynomial(
        product
    );
}

/* =========================================================
 * Direct DP polynomial
 *
 * This is an independent non-matrix implementation of the
 * same finite-state process.
 * ========================================================= */

static Poly direct_digit_dp_polynomial(
    u64 L,
    int p,
    int r
) {
    std::vector<
        std::vector<Poly>
    > dp(
        STATES,
        std::vector<Poly>(
            1,
            Poly{0}
        )
    );

    dp[0][0] =
        Poly{1};

    const auto digits =
        digits_of(
            L,
            p,
            r
        );

    for (int i = 0;
         i < r;
         ++i) {

        std::vector<
            std::vector<Poly>
        > next(
            STATES
        );

        const u64 Ldigit =
            digits[i];

        for (int s = 0;
             s < STATES;
             ++s) {

            if (
                dp[s].empty()
            ) {
                continue;
            }

            const Poly source =
                dp[s][0];

            if (
                source.size() == 1 &&
                source[0] == 0
            ) {
                continue;
            }

            const int borrow =
                s / 3;

            const int relation =
                s % 3;

            for (
                u64 ydigit = 0;
                ydigit <
                    static_cast<u64>(p);
                ++ydigit
            ) {
                int new_relation =
                    relation;

                if (
                    ydigit >
                    Ldigit
                ) {
                    new_relation = 1;
                } else if (
                    ydigit <
                    Ldigit
                ) {
                    new_relation = 2;
                }

                int new_borrow = 0;
                int added = 0;

                if (
                    Ldigit <
                    ydigit +
                    static_cast<u64>(borrow)
                ) {
                    new_borrow = 1;
                    added = 1;
                }

                const int ns =
                    new_borrow * 3 +
                    new_relation;

                Poly shifted =
                    poly_shift(
                        source,
                        added
                    );

                if (
                    next[ns].empty()
                ) {
                    next[ns].push_back(
                        Poly{0}
                    );
                }

                poly_add_inplace(
                    next[ns][0],
                    shifted
                );
            }
        }

        /*
         * Normalize representation so every state has
         * exactly one polynomial.
         */
        dp.swap(next);

        for (int s = 0;
             s < STATES;
             ++s) {
            if (
                dp[s].empty()
            ) {
                dp[s].push_back(
                    Poly{0}
                );
            }
        }
    }

    Poly result{0};

    for (int borrow = 0;
         borrow <= 1;
         ++borrow) {

        const int state =
            borrow * 3 + 1;

        if (
            dp[state].empty()
        ) {
            continue;
        }

        poly_add_inplace(
            result,
            dp[state][0]
        );
    }

    trim_poly(result);
    return result;
}

/* =========================================================
 * Compute one digit-block matrix independently
 *
 * Used to test block composition.
 * ========================================================= */

static Matrix digit_product(
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

    Matrix product =
        matrix_identity();

    for (int i = 0;
         i < r;
         ++i) {

        product =
            matrix_multiply(
                product,
                digit_matrix(
                    digits[i],
                    p
                )
            );
    }

    return product;
}

/* =========================================================
 * Split L's digits into low/high blocks
 *
 * T(L) must equal:
 *
 *   T(low) * T(high)
 *
 * because the transfer matrices compose digit by digit.
 * ========================================================= */

static bool verify_split(
    u64 L,
    int p,
    int r,
    bool print_failure
) {
    if (r <= 1) {
        return true;
    }

    const int split =
        r / 2;

    const u64 low_power =
        prime_power(
            p,
            split
        );

    const u64 low =
        L % low_power;

    const u64 high =
        L / low_power;

    const Matrix low_matrix =
        digit_product(
            low,
            p,
            split
        );

    const Matrix high_matrix =
        digit_product(
            high,
            p,
            r - split
        );

    const Matrix composed =
        matrix_multiply(
            low_matrix,
            high_matrix
        );

    const Matrix whole =
        digit_product(
            L,
            p,
            r
        );

    for (int i = 0;
         i < STATES;
         ++i) {

        for (int j = 0;
             j < STATES;
             ++j) {

            if (
                !poly_equal(
                    composed.cell[i][j],
                    whole.cell[i][j]
                )
            ) {
                if (print_failure) {
                    std::cout
                        << "FAIL_SPLIT "
                        << "L=" << L
                        << " p=" << p
                        << " r=" << r
                        << " split="
                        << split
                        << " state="
                        << i
                        << "->"
                        << j
                        << "\n";
                }

                return false;
            }
        }
    }

    return true;
}

/* =========================================================
 * Verify one local polynomial
 * ========================================================= */

static bool verify_case(
    u64 L,
    int p,
    int r,
    bool brute_exact,
    bool print_failure,
    Stats& stats
) {
    bool ok = true;

    /*
     * --------------------------------------------------
     * Exact polynomial
     * --------------------------------------------------
     */

    if (brute_exact) {
        const Poly exact =
            exact_local_polynomial(
                L,
                p,
                r
            );

        const Poly matrix =
            transfer_polynomial(
                L,
                p,
                r
            );

        if (
            !poly_equal(
                exact,
                matrix
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
                    << "\n"
                    << "exact="
                    << exact
                    << "\nmatrix="
                    << matrix
                    << "\n";
            }
        }

        const Poly dp =
            direct_digit_dp_polynomial(
                L,
                p,
                r
            );

        if (
            !poly_equal(
                exact,
                dp
            )
        ) {
            ok = false;
            ++stats.dp_failures;

            if (print_failure) {
                std::cout
                    << "FAIL_DP "
                    << "L=" << L
                    << " p=" << p
                    << " r=" << r
                    << "\n"
                    << "exact="
                    << exact
                    << "\ndp="
                    << dp
                    << "\n";
            }
        }
    }

    /*
     * --------------------------------------------------
     * Matrix product vs direct DP
     * --------------------------------------------------
     */

    const Poly matrix =
        transfer_polynomial(
            L,
            p,
            r
        );

    const Poly dp =
        direct_digit_dp_polynomial(
            L,
            p,
            r
        );

    if (
        !poly_equal(
            matrix,
            dp
        )
    ) {
        ok = false;
        ++stats.matrix_failures;

        if (print_failure) {
            std::cout
                << "FAIL_MATRIX_DP "
                << "L=" << L
                << " p=" << p
                << " r=" << r
                << "\n";
        }
    }

    /*
     * --------------------------------------------------
     * Split composition
     * --------------------------------------------------
     */

    if (
        !verify_split(
            L,
            p,
            r,
            print_failure
        )
    ) {
        ok = false;
        ++stats.split_failures;
    }

    /*
     * --------------------------------------------------
     * Degree bound
     * --------------------------------------------------
     */

    const Poly result =
        matrix;

    for (
        std::size_t k = result.size();
        k > 0;
        --k
    ) {
        if (
            result[k - 1] != 0
        ) {
            const std::size_t degree =
                k - 1;

            if (
                degree >
                static_cast<std::size_t>(
                    r
                )
            ) {
                ok = false;
                ++stats.degree_failures;

                if (print_failure) {
                    std::cout
                        << "FAIL_DEGREE "
                        << "L=" << L
                        << " p=" << p
                        << " r=" << r
                        << " degree="
                        << degree
                        << "\n";
                }
            }

            break;
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
        << "START EXPERIMENT 335\n"
        << "KUMMER LOCAL GENERATING POLYNOMIAL\n"
        << "DO DIGIT TRANSFER MATRICES REPRODUCE THE FULL GAP POLYNOMIAL?\n\n";

    Stats stats;

    /*
     * --------------------------------------------------
     * Phase 1:
     * Exhaustive small local blocks.
     *
     * Full exact polynomial comparison.
     * --------------------------------------------------
     */

    constexpr u64 FULL_LIMIT = 1500;

    for (int p : PRIMES) {
        for (int r = 1;
             r <= 12;
             ++r) {

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
                ++stats.tested;

                verify_case(
                    L,
                    p,
                    r,
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
     * Random full polynomial checks.
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x335335335ULL
    );

    constexpr u64 RANDOM_FULL = 2500;

    for (
        u64 i = 0;
        i < RANDOM_FULL;
        ++i
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
            power > 50000
        ) {
            --i;
            continue;
        }

        const u64 L =
            rng() % power;

        ++stats.tested;

        verify_case(
            L,
            p,
            r,
            true,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 3:
     * Larger blocks.
     *
     * No full enumeration.
     *
     * The matrix and direct digit-DP constructions must
     * agree, and block composition must hold.
     * --------------------------------------------------
     */

    constexpr u64 RANDOM_LARGE = 4000;

    for (
        u64 i = 0;
        i < RANDOM_LARGE;
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
            false,
            stats.failed < 10,
            stats
        );
    }

    /*
     * --------------------------------------------------
     * Phase 4:
     * Targeted digit patterns.
     * --------------------------------------------------
     */

    for (int p : PRIMES) {
        for (int r = 1;
             r <= 16;
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

            for (int pattern = 0;
                 pattern < 5;
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
                         * alternating 0/max
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
                         * alternating max/0
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
                         * sparse high digits
                         */
                        if (
                            i >= r / 2
                        ) {
                            digits[i] = 1;
                        }
                    }
                }

                u64 L = 0;
                u64 place = 1;
                bool overflow = false;

                for (u64 digit : digits) {
                    if (
                        digit != 0 &&
                        place >
                            UINT64_MAX / digit
                    ) {
                        overflow = true;
                        break;
                    }

                    L +=
                        digit * place;

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
                    L < power
                ) {
                    candidates.push_back(L);
                }
            }

            for (u64 L : candidates) {
                ++stats.targeted;
                ++stats.tested;

                const bool full =
                    power <= FULL_LIMIT;

                verify_case(
                    L,
                    p,
                    r,
                    full,
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
        << "exact_failures="
        << stats.exact_failures
        << "\n"
        << "dp_failures="
        << stats.dp_failures
        << "\n"
        << "matrix_failures="
        << stats.matrix_failures
        << "\n"
        << "split_failures="
        << stats.split_failures
        << "\n"
        << "total_failures="
        << stats.total_failures
        << "\n"
        << "degree_failures="
        << stats.degree_failures
        << "\n"
        << "targeted="
        << stats.targeted
        << "\n"
        << "============================\n";

    std::cout
        << "FINISHED EXPERIMENT 335\n";

    return 0;
}
