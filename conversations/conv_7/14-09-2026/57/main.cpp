#include <array>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

using SPoly = std::vector<i128>;

struct SemiprimeCase {
    u64 N;
    u64 p;
    u64 q;
    u64 sum;
};

struct Features {
    /*
     * 10 observable quantities:
     *
     * 0  N
     *
     * 1  B20
     * 2  B21
     * 3  B22
     *
     * 4  B30
     * 5  B31
     * 6  B32
     *
     * 7  B50
     * 8  B51
     * 9  B52
     */
    std::array<u64, 10> x{};
};

struct VariableExponent {
    /*
     * Variables:
     *
     * 0..9 = observable features
     * 10   = S = p+q
     */
    std::array<int, 11> exponent{};
};

struct Relation {
    bool found = false;

    int degree = -1;

    std::vector<
        VariableExponent
    > monomials;

    /*
     * Integer coefficients obtained from the exact
     * rational nullspace normalization.
     */
    std::vector<mpz_class> coefficients;

    u64 training_failures = 0;
    u64 validation_failures = 0;
};

static const std::array<int, 3> BASES = {
    2, 3, 5
};

/* =========================================================
 * GMP helpers
 * ========================================================= */

static mpq_class mpq_from_u64(
    u64 value
) {
    return mpq_class(
        mpz_class(value)
    );
}

/* =========================================================
 * Prime generation
 * ========================================================= */

static std::vector<int> sieve_primes(
    int limit
) {
    std::vector<bool> composite(
        static_cast<std::size_t>(
            limit + 1
        ),
        false
    );

    std::vector<int> primes;

    for (
        int i = 2;
        i <= limit;
        ++i
    ) {
        if (
            composite[
                static_cast<
                    std::size_t
                >(i)
            ]
        ) {
            continue;
        }

        primes.push_back(i);

        if (
            static_cast<long long>(i) *
            i <= limit
        ) {
            for (
                int j = i * i;
                j <= limit;
                j += i
            ) {
                composite[
                    static_cast<
                        std::size_t
                    >(j)
                ] = true;
            }
        }
    }

    return primes;
}

/* =========================================================
 * Kummer valuation histogram
 *
 * H_k(N,p) =
 *
 * #{ t in [0,N] :
 *    v_p(C(N,t)) = k }.
 *
 * Input is only N.
 * ========================================================= */

static std::vector<u64>
valuation_histogram(
    u64 N,
    int p
) {
    std::vector<u64> digits;

    u64 x = N;

    const u64 P =
        static_cast<u64>(p);

    while (x > 0) {
        digits.push_back(
            x % P
        );

        x /= P;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    const int r =
        static_cast<int>(
            digits.size()
        );

    /*
     * dp[borrow][valuation]
     */
    std::vector<
        std::vector<u64>
    > dp(
        2,
        std::vector<u64>(
            static_cast<
                std::size_t
            >(r + 1),
            0
        )
    );

    std::vector<
        std::vector<u64>
    > next(
        2,
        std::vector<u64>(
            static_cast<
                std::size_t
            >(r + 1),
            0
        )
    );

    dp[0][0] = 1;

    for (
        int i = 0;
        i < r;
        ++i
    ) {
        for (auto& row : next) {
            std::fill(
                row.begin(),
                row.end(),
                0
            );
        }

        const u64 d =
            digits[
                static_cast<
                    std::size_t
                >(i)
            ];

        /*
         * Borrow = 0.
         */
        for (
            int k = 0;
            k <= i;
            ++k
        ) {
            const u64 count =
                dp[0][
                    static_cast<
                        std::size_t
                    >(k)
                ];

            if (count == 0) {
                continue;
            }

            /*
             * No new borrow.
             */
            next[0][
                static_cast<
                    std::size_t
                >(k)
            ] +=
                count * (d + 1);

            /*
             * New borrow.
             */
            next[1][
                static_cast<
                    std::size_t
                >(k + 1)
            ] +=
                count *
                (
                    static_cast<u64>(p)
                    - d
                    - 1
                );
        }

        /*
         * Borrow = 1.
         */
        for (
            int k = 0;
            k <= i;
            ++k
        ) {
            const u64 count =
                dp[1][
                    static_cast<
                        std::size_t
                    >(k)
                ];

            if (count == 0) {
                continue;
            }

            /*
             * Borrow resolves.
             */
            next[0][
                static_cast<
                    std::size_t
                >(k)
            ] +=
                count * d;

            /*
             * Borrow continues.
             */
            next[1][
                static_cast<
                    std::size_t
                >(k + 1)
            ] +=
                count *
                (
                    static_cast<u64>(p)
                    - d
                );
        }

        dp.swap(next);
    }

    return dp[0];
}

/* =========================================================
 * Feature extraction
 * ========================================================= */

static Features extract_features(
    u64 N
) {
    Features result{};

    result.x[0] = N;

    int index = 1;

    for (int base : BASES) {
        const auto histogram =
            valuation_histogram(
                N,
                base
            );

        for (int k = 0; k <= 2; ++k) {
            result.x[
                static_cast<
                    std::size_t
                >(index++)
            ] =
                histogram.size() >
                static_cast<
                    std::size_t
                >(k)
                    ? histogram[
                        static_cast<
                            std::size_t
                        >(k)
                    ]
                    : 0;
        }
    }

    return result;
}

/* =========================================================
 * Generate all monomials of total degree <= D
 *
 * Variables:
 *
 * 0..9  = observable features
 * 10    = S
 * ========================================================= */

static void generate_monomials_recursive(
    int variable,
    int remaining_degree,
    std::array<int, 11>& current,
    std::vector<
        VariableExponent
    >& result
) {
    if (
        variable == 11
    ) {
        result.push_back(
            VariableExponent{
                current
            }
        );

        return;
    }

    for (
        int exponent = 0;
        exponent <= remaining_degree;
        ++exponent
    ) {
        current[
            static_cast<
                std::size_t
            >(variable)
        ] = exponent;

        generate_monomials_recursive(
            variable + 1,
            remaining_degree - exponent,
            current,
            result
        );
    }

    current[
        static_cast<
            std::size_t
        >(variable)
    ] = 0;
}

static std::vector<
    VariableExponent
>
generate_monomials(
    int degree
) {
    std::vector<
        VariableExponent
    > result;

    std::array<int, 11> current{};

    generate_monomials_recursive(
        0,
        degree,
        current,
        result
    );

    return result;
}

/* =========================================================
 * Modular arithmetic
 * ========================================================= */

static u64 mod_mul(
    u64 a,
    u64 b,
    u64 mod
) {
    return static_cast<u64>(
        (
            static_cast<u128>(a) *
            static_cast<u128>(b)
        ) % mod
    );
}

static u64 mod_pow(
    u64 base,
    u64 exponent,
    u64 mod
) {
    u64 result = 1;

    while (exponent > 0) {
        if (exponent & 1ULL) {
            result =
                mod_mul(
                    result,
                    base,
                    mod
                );
        }

        base =
            mod_mul(
                base,
                base,
                mod
            );

        exponent >>= 1;
    }

    return result;
}

static u64 mod_inverse(
    u64 x,
    u64 mod
) {
    return mod_pow(
        x,
        mod - 2,
        mod
    );
}

/* =========================================================
 * Monomial evaluation modulo p
 * ========================================================= */

static u64 evaluate_monomial_mod(
    const Features& features,
    u64 sum,
    const VariableExponent& monomial,
    u64 mod
) {
    u64 result = 1;

    for (
        int i = 0;
        i < 11;
        ++i
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<
                    std::size_t
                >(i)
            ];

        if (exponent == 0) {
            continue;
        }

        const u64 value =
            i < 10
                ? features.x[
                    static_cast<
                        std::size_t
                    >(i)
                ]
                : sum;

        result =
            mod_mul(
                result,
                mod_pow(
                    value % mod,
                    static_cast<
                        u64
                    >(exponent),
                    mod
                ),
                mod
            );
    }

    return result;
}

/* =========================================================
 * Exact monomial evaluation
 *
 * GMP is used because feature powers can become large.
 * ========================================================= */

static mpz_class evaluate_monomial_exact(
    const Features& features,
    u64 sum,
    const VariableExponent& monomial
) {
    mpz_class result(1);

    for (
        int i = 0;
        i < 11;
        ++i
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<
                    std::size_t
                >(i)
            ];

        if (exponent == 0) {
            continue;
        }

        const u64 value =
            i < 10
                ? features.x[
                    static_cast<
                        std::size_t
                    >(i)
                ]
                : sum;

        mpz_class base(value);

        for (
            int j = 0;
            j < exponent;
            ++j
        ) {
            result *= base;
        }
    }

    return result;
}

/* =========================================================
 * Modular nullspace
 *
 * Returns:
 *
 * 0 = no nonzero nullspace vector
 * 1 = at least one nullspace vector
 *
 * The returned vector is one nonzero null vector.
 * ========================================================= */

static bool find_modular_nullvector(
    const std::vector<Features>& data,
    const std::vector<u64>& sums,
    const std::vector<
        VariableExponent
    >& monomials,
    u64 mod,
    std::vector<u64>& nullvector
) {
    const std::size_t rows =
        data.size();

    const std::size_t cols =
        monomials.size();

    /*
     * Augmented variable matrix is not needed:
     * this is a homogeneous system.
     */
    std::vector<
        std::vector<u64>
    > matrix(
        rows,
        std::vector<u64>(
            cols,
            0
        )
    );

    for (
        std::size_t r = 0;
        r < rows;
        ++r
    ) {
        for (
            std::size_t c = 0;
            c < cols;
            ++c
        ) {
            matrix[r][c] =
                evaluate_monomial_mod(
                    data[r],
                    sums[r],
                    monomials[c],
                    mod
                );
        }
    }

    std::size_t pivot_row = 0;

    std::vector<int> pivot_column_for_row;

    for (
        std::size_t col = 0;
        col < cols &&
        pivot_row < rows;
        ++col
    ) {
        std::size_t pivot =
            pivot_row;

        while (
            pivot < rows &&
            matrix[pivot][col] == 0
        ) {
            ++pivot;
        }

        if (pivot == rows) {
            continue;
        }

        std::swap(
            matrix[pivot],
            matrix[pivot_row]
        );

        const u64 inverse =
            mod_inverse(
                matrix[pivot_row][col],
                mod
            );

        for (
            std::size_t j = col;
            j < cols;
            ++j
        ) {
            matrix[pivot_row][j] =
                mod_mul(
                    matrix[pivot_row][j],
                    inverse,
                    mod
                );
        }

        /*
         * Eliminate this pivot from every other row.
         */
        for (
            std::size_t r = 0;
            r < rows;
            ++r
        ) {
            if (r == pivot_row) {
                continue;
            }

            const u64 factor =
                matrix[r][col];

            if (factor == 0) {
                continue;
            }

            for (
                std::size_t j = col;
                j < cols;
                ++j
            ) {
                const u64 product =
                    mod_mul(
                        factor,
                        matrix[pivot_row][j],
                        mod
                    );

                if (
                    matrix[r][j] >=
                    product
                ) {
                    matrix[r][j] -=
                        product;
                } else {
                    matrix[r][j] =
                        mod -
                        (
                            product -
                            matrix[r][j]
                        );
                }
            }
        }

        ++pivot_row;
    }

    /*
     * If rank < number of columns, a nullspace exists.
     */
    if (
        pivot_row >= cols
    ) {
        return false;
    }

    /*
     * Record pivot columns.
     */
    std::vector<bool> pivot(
        cols,
        false
    );

    std::size_t discovered_pivots = 0;

    for (
        std::size_t r = 0;
        r < rows &&
        discovered_pivots < pivot_row;
        ++r
    ) {
        for (
            std::size_t c = 0;
            c < cols;
            ++c
        ) {
            if (
                matrix[r][c] != 0
            ) {
                pivot[c] = true;
                ++discovered_pivots;
                break;
            }
        }
    }

    /*
     * Pick the first free variable.
     */
    std::size_t free_column =
        cols;

    for (
        std::size_t c = 0;
        c < cols;
        ++c
    ) {
        if (!pivot[c]) {
            free_column = c;
            break;
        }
    }

    if (
        free_column == cols
    ) {
        return false;
    }

    nullvector.assign(
        cols,
        0
    );

    nullvector[
        free_column
    ] = 1;

    /*
     * Because the matrix is in reduced row-echelon form,
     * each pivot variable is determined directly by the
     * free variables.
     */
    for (
        std::size_t r = 0;
        r < pivot_row;
        ++r
    ) {
        std::size_t pivot_col =
            cols;

        for (
            std::size_t c = 0;
            c < cols;
            ++c
        ) {
            if (
                matrix[r][c] != 0
            ) {
                pivot_col = c;
                break;
            }
        }

        if (
            pivot_col == cols
        ) {
            continue;
        }

        const u64 coefficient =
            matrix[r][
                free_column
            ];

        if (coefficient == 0) {
            nullvector[pivot_col] = 0;
        } else {
            nullvector[pivot_col] =
                (
                    mod -
                    coefficient
                ) % mod;
        }
    }

    return true;
}

/* =========================================================
 * Exact rational nullspace
 *
 * We only call this after the modular tests show that a
 * nullspace exists in BOTH primes.
 * ========================================================= */

static bool exact_nullvector(
    const std::vector<Features>& data,
    const std::vector<u64>& sums,
    const std::vector<
        VariableExponent
    >& monomials,
    std::vector<mpq_class>& vector_out
) {
    const std::size_t rows =
        data.size();

    const std::size_t cols =
        monomials.size();

    std::vector<
        std::vector<mpq_class>
    > matrix(
        rows,
        std::vector<mpq_class>(
            cols,
            mpq_class(0)
        )
    );

    for (
        std::size_t r = 0;
        r < rows;
        ++r
    ) {
        for (
            std::size_t c = 0;
            c < cols;
            ++c
        ) {
            matrix[r][c] =
                mpq_class(
                    evaluate_monomial_exact(
                        data[r],
                        sums[r],
                        monomials[c]
                    )
                );
        }
    }

    std::size_t pivot_row = 0;

    std::vector<int> pivot_column;

    for (
        std::size_t col = 0;
        col < cols &&
        pivot_row < rows;
        ++col
    ) {
        std::size_t pivot =
            pivot_row;

        while (
            pivot < rows &&
            matrix[pivot][col] == 0
        ) {
            ++pivot;
        }

        if (pivot == rows) {
            continue;
        }

        std::swap(
            matrix[pivot],
            matrix[pivot_row]
        );

        const mpq_class pivot_value =
            matrix[pivot_row][col];

        for (
            std::size_t j = col;
            j < cols;
            ++j
        ) {
            matrix[pivot_row][j] /=
                pivot_value;
        }

        /*
         * Gauss-Jordan.
         */
        for (
            std::size_t r = 0;
            r < rows;
            ++r
        ) {
            if (r == pivot_row) {
                continue;
            }

            const mpq_class factor =
                matrix[r][col];

            if (factor == 0) {
                continue;
            }

            for (
                std::size_t j = col;
                j < cols;
                ++j
            ) {
                matrix[r][j] -=
                    factor *
                    matrix[pivot_row][j];
            }
        }

        pivot_column.push_back(
            static_cast<int>(col)
        );

        ++pivot_row;
    }

    if (
        pivot_row >= cols
    ) {
        /*
         * No nullspace.
         */
        return false;
    }

    /*
     * Mark pivot columns.
     */
    std::vector<bool> is_pivot(
        cols,
        false
    );

    for (
        int col :
        pivot_column
    ) {
        is_pivot[
            static_cast<
                std::size_t
            >(col)
        ] = true;
    }

    /*
     * Select first free variable = 1.
     */
    std::size_t free_column =
        cols;

    for (
        std::size_t c = 0;
        c < cols;
        ++c
    ) {
        if (!is_pivot[c]) {
            free_column = c;
            break;
        }
    }

    if (
        free_column == cols
    ) {
        return false;
    }

    vector_out.assign(
        cols,
        mpq_class(0)
    );

    vector_out[
        free_column
    ] = mpq_class(1);

    /*
     * For each pivot row:
     *
     * x_pivot +
     * sum free a_j x_j = 0
     *
     */
    for (
        std::size_t r = 0;
        r < pivot_column.size();
        ++r
    ) {
        const int pivot_col =
            pivot_column[r];

        vector_out[
            static_cast<
                std::size_t
            >(pivot_col)
        ] =
            -matrix[r][free_column];
    }

    return true;
}

/* =========================================================
 * Normalize rational relation to primitive integers
 * ========================================================= */

static bool normalize_integer_relation(
    const std::vector<mpq_class>& rational,
    std::vector<mpz_class>& integer
) {
    if (rational.empty()) {
        return false;
    }

    mpz_class common_denominator = 1;

    for (
        const mpq_class& x :
        rational
    ) {
        common_denominator =
            lcm(
                common_denominator,
                x.get_den()
            );
    }

    integer.resize(
        rational.size()
    );

    for (
        std::size_t i = 0;
        i < rational.size();
        ++i
    ) {
        integer[i] =
            rational[i].get_num() *
            (
                common_denominator /
                rational[i].get_den()
            );
    }

    /*
     * Divide by the gcd.
     */
    mpz_class g = 0;

    for (
        const mpz_class& x :
        integer
    ) {
        if (x == 0) {
            continue;
        }

        if (g == 0) {
            g = abs(x);
        } else {
            g = gcd(
                g,
                abs(x)
            );
        }
    }

    if (g == 0) {
        return false;
    }

    for (
        mpz_class& x :
        integer
    ) {
        x /= g;
    }

    /*
     * Canonical sign:
     * first nonzero coefficient positive.
     */
    for (
        const mpz_class& x :
        integer
    ) {
        if (x == 0) {
            continue;
        }

        if (x < 0) {
            for (
                mpz_class& y :
                integer
            ) {
                y = -y;
            }
        }

        break;
    }

    return true;
}

/* =========================================================
 * Evaluate integer relation exactly
 * ========================================================= */

static mpz_class evaluate_relation(
    const Features& features,
    u64 sum,
    const std::vector<
        VariableExponent
    >& monomials,
    const std::vector<
        mpz_class
    >& coefficients
) {
    mpz_class result = 0;

    for (
        std::size_t i = 0;
        i < monomials.size();
        ++i
    ) {
        if (
            coefficients[i] == 0
        ) {
            continue;
        }

        const mpz_class monomial =
            evaluate_monomial_exact(
                features,
                sum,
                monomials[i]
            );

        result +=
            coefficients[i] *
            monomial;
    }

    return result;
}

/* =========================================================
 * Pretty-print relation
 * ========================================================= */

static std::string variable_name(
    int index
) {
    static const std::array<
        const char*,
        11
    > names = {
        "N",
        "B20",
        "B21",
        "B22",
        "B30",
        "B31",
        "B32",
        "B50",
        "B51",
        "B52",
        "S"
    };

    return names[
        static_cast<
            std::size_t
        >(index)
    ];
}

static std::string monomial_name(
    const VariableExponent& monomial
) {
    std::string result;

    for (
        int i = 0;
        i < 11;
        ++i
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<
                    std::size_t
                >(i)
            ];

        if (exponent == 0) {
            continue;
        }

        if (!result.empty()) {
            result += "*";
        }

        result +=
            variable_name(i);

        if (exponent != 1) {
            result += "^";
            result +=
                std::to_string(
                    exponent
                );
        }
    }

    if (result.empty()) {
        return "1";
    }

    return result;
}

/* =========================================================
 * Search relation of fixed degree
 * ========================================================= */

static Relation search_relation_degree(
    int degree,
    const std::vector<Features>& training,
    const std::vector<u64>& training_sums,
    const std::vector<Features>& validation,
    const std::vector<u64>& validation_sums
) {
    Relation relation;

    relation.degree =
        degree;

    const auto monomials =
        generate_monomials(
            degree
        );

    std::cout
        << "DEGREE="
        << degree
        << " MONOMIALS="
        << monomials.size()
        << "\n";

    constexpr u64 MOD1 =
        1000000007ULL;

    constexpr u64 MOD2 =
        1000000009ULL;

    std::vector<u64> null1;

    const bool exists1 =
        find_modular_nullvector(
            training,
            training_sums,
            monomials,
            MOD1,
            null1
        );

    std::vector<u64> null2;

    const bool exists2 =
        find_modular_nullvector(
            training,
            training_sums,
            monomials,
            MOD2,
            null2
        );

    std::cout
        << "MOD1_NULLSPACE="
        << (exists1 ? 1 : 0)
        << " MOD2_NULLSPACE="
        << (exists2 ? 1 : 0)
        << "\n";

    if (
        !exists1 ||
        !exists2
    ) {
        std::cout
            << "NO_MODULAR_RELATION\n";

        return relation;
    }

    /*
     * Exact rational nullspace.
     */
    std::vector<mpq_class> rational_relation;

    if (
        !exact_nullvector(
            training,
            training_sums,
            monomials,
            rational_relation
        )
    ) {
        std::cout
            << "EXACT_NULLSPACE_FAILED\n";

        return relation;
    }

    std::vector<mpz_class> coefficients;

    if (
        !normalize_integer_relation(
            rational_relation,
            coefficients
        )
    ) {
        std::cout
            << "NORMALIZATION_FAILED\n";

        return relation;
    }

    relation.found = true;

    relation.monomials =
        monomials;

    relation.coefficients =
        coefficients;

    /*
     * Training verification.
     */
    for (
        std::size_t i = 0;
        i < training.size();
        ++i
    ) {
        const mpz_class value =
            evaluate_relation(
                training[i],
                training_sums[i],
                monomials,
                coefficients
            );

        if (value != 0) {
            ++relation.training_failures;
        }
    }

    /*
     * Completely unseen validation data.
     */
    for (
        std::size_t i = 0;
        i < validation.size();
        ++i
    ) {
        const mpz_class value =
            evaluate_relation(
                validation[i],
                validation_sums[i],
                monomials,
                coefficients
            );

        if (
            value != 0
        ) {
            ++relation.validation_failures;

            if (
                relation.validation_failures <= 5
            ) {
                std::cout
                    << "VALIDATION_FAILURE "
                    << "index="
                    << i
                    << "\n"
                    << "value="
                    << value.get_str()
                    << "\n";
            }
        }
    }

    std::cout
        << "TRAINING_FAILURES="
        << relation.training_failures
        << "\n";

    std::cout
        << "VALIDATION_FAILURES="
        << relation.validation_failures
        << "\n";

    /*
     * Print relation only when it actually survives
     * unseen data.
     */
    if (
        relation.validation_failures == 0
    ) {
        std::cout
            << "VALIDATED_RELATION=YES\n";

        std::size_t nonzero = 0;

        for (
            const mpz_class& coefficient :
            coefficients
        ) {
            if (coefficient != 0) {
                ++nonzero;
            }
        }

        std::cout
            << "NONZERO_TERMS="
            << nonzero
            << "\n";

        std::cout
            << "RELATION:\n";

        bool first = true;

        for (
            std::size_t i = 0;
            i < coefficients.size();
            ++i
        ) {
            if (
                coefficients[i] == 0
            ) {
                continue;
            }

            if (!first) {
                std::cout
                    << "\n";
            }

            first = false;

            std::cout
                << "("
                << coefficients[i].get_str()
                << ")*"
                << monomial_name(
                    monomials[i]
                );
        }

        std::cout
            << "\n";
    }

    return relation;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 346\n"
        << "ALGEBRAIC RELATION SEARCH FOR p+q\n"
        << "CAN OBSERVABLE KUMMER DATA AND S=p+q "
           "SATISFY A LOW-DEGREE POLYNOMIAL RELATION?\n\n";

    /*
     * --------------------------------------------------
     * Prime pool
     * --------------------------------------------------
     */

    const auto all_primes =
        sieve_primes(
            997
        );

    std::vector<int> usable_primes;

    for (int p : all_primes) {
        if (p >= 101) {
            usable_primes.push_back(
                p
            );
        }
    }

    /*
     * --------------------------------------------------
     * Dataset
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x346346346ULL
    );

    constexpr std::size_t TOTAL_CASES =
        1000;

    constexpr std::size_t TRAINING_CASES =
        700;

    constexpr std::size_t VALIDATION_CASES =
        300;

    std::set<u64> used_N;

    std::vector<SemiprimeCase> cases;

    cases.reserve(
        TOTAL_CASES
    );

    while (
        cases.size() <
        TOTAL_CASES
    ) {
        const int p =
            usable_primes[
                rng() %
                usable_primes.size()
            ];

        const int q =
            usable_primes[
                rng() %
                usable_primes.size()
            ];

        if (p == q) {
            continue;
        }

        const int small =
            std::min(
                p,
                q
            );

        const int large =
            std::max(
                p,
                q
            );

        const u64 N =
            static_cast<u64>(
                small
            ) *
            static_cast<u64>(
                large
            );

        if (
            !used_N.insert(N).second
        ) {
            continue;
        }

        cases.push_back(
            SemiprimeCase{
                N,
                static_cast<u64>(
                    small
                ),
                static_cast<u64>(
                    large
                ),
                static_cast<u64>(
                    small + large
                )
            }
        );
    }

    std::shuffle(
        cases.begin(),
        cases.end(),
        rng
    );

    /*
     * --------------------------------------------------
     * Extract observable features.
     * --------------------------------------------------
     */

    std::vector<Features> training_features;
    std::vector<u64> training_sums;

    std::vector<Features> validation_features;
    std::vector<u64> validation_sums;

    training_features.reserve(
        TRAINING_CASES
    );

    training_sums.reserve(
        TRAINING_CASES
    );

    validation_features.reserve(
        VALIDATION_CASES
    );

    validation_sums.reserve(
        VALIDATION_CASES
    );

    for (
        std::size_t i = 0;
        i < cases.size();
        ++i
    ) {
        const Features features =
            extract_features(
                cases[i].N
            );

        if (
            i < TRAINING_CASES
        ) {
            training_features.push_back(
                features
            );

            training_sums.push_back(
                cases[i].sum
            );
        } else {
            validation_features.push_back(
                features
            );

            validation_sums.push_back(
                cases[i].sum
            );
        }
    }

    std::cout
        << "TRAINING_CASES="
        << training_features.size()
        << "\n";

    std::cout
        << "VALIDATION_CASES="
        << validation_features.size()
        << "\n";

    std::cout
        << "OBSERVABLE_FEATURES=10\n";

    std::cout
        << "VARIABLES="
        << "N,"
        << "B20,B21,B22,"
        << "B30,B31,B32,"
        << "B50,B51,B52,"
        << "S"
        << "\n\n";

    /*
     * --------------------------------------------------
     * Search fixed algebraic degrees.
     *
     * Do NOT increase degree because a validation
     * failure occurs.
     * --------------------------------------------------
     */

    Relation validated_relation;

    for (
        int degree = 1;
        degree <= 3;
        ++degree
    ) {
        std::cout
            << "============================\n";

        const Relation relation =
            search_relation_degree(
                degree,
                training_features,
                training_sums,
                validation_features,
                validation_sums
            );

        if (
            relation.found &&
            relation.training_failures == 0 &&
            relation.validation_failures == 0
        ) {
            validated_relation =
                relation;

            break;
        }
    }

    /*
     * --------------------------------------------------
     * Result
     * --------------------------------------------------
     */

    std::cout
        << "\n"
        << "============================\n";

    if (
        validated_relation.found &&
        validated_relation.validation_failures == 0
    ) {
        std::cout
            << "DISCOVERY_STATUS=VALIDATED_RELATION\n";

        std::cout
            << "DISCOVERY_DEGREE="
            << validated_relation.degree
            << "\n";
    } else {
        std::cout
            << "DISCOVERY_STATUS="
               "NO_VALIDATED_RELATION_DEGREE_1_TO_3\n";
    }

    std::cout
        << "============================\n"
        << "FINISHED EXPERIMENT 346\n";

    return 0;
}
