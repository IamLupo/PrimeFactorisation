#include <array>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>
#include <numeric>

#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

/*
 * Observable variables:
 *
 * 0  = N
 *
 * For every base b:
 *
 * C_b = sum_{i>=1} d_{i-1}(d_i+1)
 * B_b = sum_{i>=1} d_i(b-d_{i-1})
 * G_b = d_1*b + sum_{i>=2} d_i*b*d_{i-2}
 *
 * 6 bases * 3 = 18
 * + N = 19 observable variables.
 *
 * Variable 19 = S = p+q.
 */
static constexpr int OBSERVABLE_COUNT = 19;
static constexpr int TOTAL_VARIABLES = 20;

struct SemiprimeCase {
    u64 N;
    u64 p;
    u64 q;
    u64 S;
};

struct Features {
    std::array<u64, OBSERVABLE_COUNT> x{};
};

struct Monomial {
    std::array<int, TOTAL_VARIABLES> exponent{};
};

struct RREF {
    std::vector<std::vector<u64>> matrix;
    std::vector<int> pivot_columns;
    std::vector<bool> is_pivot;
    std::size_t rank = 0;
};

struct Relation {
    bool found = false;
    bool exact_training = false;
    bool exact_validation = false;

    int degree = -1;

    std::vector<Monomial> monomials;
    std::vector<mpz_class> coefficients;

    u64 training_failures = 0;
    u64 validation_failures = 0;
};

/* =========================================================
 * GMP conversion helpers
 *
 * IMPORTANT:
 * Do not construct mpz_class directly from long long.
 * Use strings for signed 64-bit values.
 * ========================================================= */

static mpz_class mpz_from_signed64(
    std::int64_t value
) {
    return mpz_class(
        std::to_string(value)
    );
}

static mpz_class mpz_from_u64(
    u64 value
) {
    return mpz_class(
        std::to_string(value)
    );
}

/* =========================================================
 * Prime sieve
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
                static_cast<std::size_t>(i)
            ]
        ) {
            continue;
        }

        primes.push_back(i);

        if (
            static_cast<long long>(i) * i
            <= limit
        ) {
            for (
                int j = i * i;
                j <= limit;
                j += i
            ) {
                composite[
                    static_cast<std::size_t>(j)
                ] = true;
            }
        }
    }

    return primes;
}

/* =========================================================
 * Base-b digits
 * ========================================================= */

static std::vector<u64> digits_of(
    u64 N,
    int base
) {
    std::vector<u64> digits;

    u64 x = N;
    const u64 b =
        static_cast<u64>(base);

    while (x > 0) {
        digits.push_back(
            x % b
        );

        x /= b;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

/* =========================================================
 * Local recurrence invariants
 *
 * d_{i-1} A_{i+1}
 *
 * =
 * [
 *   d_{i-1}(d_i+1)
 *   +
 *   d_i(b-d_{i-1}) z
 * ] A_i
 *
 * -
 * d_i*b*d_{i-2}*z*A_{i-1}
 * ========================================================= */

static std::array<u64, 3>
local_invariants(
    u64 N,
    int base
) {
    const auto digits =
        digits_of(
            N,
            base
        );

    const int r =
        static_cast<int>(
            digits.size()
        );

    u128 C = 0;
    u128 B = 0;
    u128 G = 0;

    /*
     * Constant and z coefficients.
     */
    for (
        int i = 1;
        i < r;
        ++i
    ) {
        const u128 d_prev =
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(
                        i - 1
                    )
                ]
            );

        const u128 d_curr =
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(
                        i
                    )
                ]
            );

        C +=
            d_prev *
            (d_curr + 1);

        B +=
            d_curr *
            (
                static_cast<u128>(base)
                - d_prev
            );
    }

    /*
     * Previous-layer coefficient.
     */
    if (r >= 2) {
        G +=
            static_cast<u128>(
                digits[1]
            ) *
            static_cast<u128>(base);
    }

    for (
        int i = 2;
        i < r;
        ++i
    ) {
        const u128 d_im2 =
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(
                        i - 2
                    )
                ]
            );

        const u128 d_i =
            static_cast<u128>(
                digits[
                    static_cast<std::size_t>(
                        i
                    )
                ]
            );

        G +=
            d_i *
            static_cast<u128>(base) *
            d_im2;
    }

    if (
        C > static_cast<u128>(UINT64_MAX) ||
        B > static_cast<u128>(UINT64_MAX) ||
        G > static_cast<u128>(UINT64_MAX)
    ) {
        std::cerr
            << "FATAL: local invariant overflow\n";

        std::exit(1);
    }

    return {
        static_cast<u64>(C),
        static_cast<u64>(B),
        static_cast<u64>(G)
    };
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
        const auto inv =
            local_invariants(
                N,
                base
            );

        result.x[
            static_cast<std::size_t>(
                index++
            )
        ] = inv[0];

        result.x[
            static_cast<std::size_t>(
                index++
            )
        ] = inv[1];

        result.x[
            static_cast<std::size_t>(
                index++
            )
        ] = inv[2];
    }

    return result;
}

/* =========================================================
 * Monomial generation
 * ========================================================= */

static void generate_monomials_recursive(
    int variable,
    int remaining_degree,
    std::array<int, TOTAL_VARIABLES>& current,
    std::vector<Monomial>& result
) {
    if (
        variable ==
        TOTAL_VARIABLES
    ) {
        result.push_back(
            Monomial{
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
            static_cast<std::size_t>(
                variable
            )
        ] = exponent;

        generate_monomials_recursive(
            variable + 1,
            remaining_degree - exponent,
            current,
            result
        );
    }

    current[
        static_cast<std::size_t>(
            variable
        )
    ] = 0;
}

static std::vector<Monomial>
generate_monomials(
    int degree
) {
    std::vector<Monomial> result;

    std::array<int, TOTAL_VARIABLES> current{};

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
    u64 value,
    u64 mod
) {
    return mod_pow(
        value,
        mod - 2,
        mod
    );
}

/* =========================================================
 * Monomial evaluation modulo prime
 * ========================================================= */

static u64 evaluate_monomial_mod(
    const Features& features,
    u64 S,
    const Monomial& monomial,
    u64 mod
) {
    u64 result = 1;

    for (
        int variable = 0;
        variable < TOTAL_VARIABLES;
        ++variable
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<std::size_t>(
                    variable
                )
            ];

        if (exponent == 0) {
            continue;
        }

        const u64 value =
            variable < OBSERVABLE_COUNT
                ? features.x[
                    static_cast<std::size_t>(
                        variable
                    )
                ]
                : S;

        result =
            mod_mul(
                result,
                mod_pow(
                    value % mod,
                    static_cast<u64>(
                        exponent
                    ),
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
 * GMP deliberately used here.
 * ========================================================= */

static mpz_class evaluate_monomial_exact(
    const Features& features,
    u64 S,
    const Monomial& monomial
) {
    mpz_class result(
        1
    );

    for (
        int variable = 0;
        variable < TOTAL_VARIABLES;
        ++variable
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<std::size_t>(
                    variable
                )
            ];

        if (exponent == 0) {
            continue;
        }

        const u64 value =
            variable < OBSERVABLE_COUNT
                ? features.x[
                    static_cast<std::size_t>(
                        variable
                    )
                ]
                : S;

        const mpz_class base =
            mpz_from_u64(
                value
            );

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
 * RREF modulo prime
 * ========================================================= */

static RREF build_rref(
    const std::vector<Features>& data,
    const std::vector<u64>& sums,
    const std::vector<Monomial>& monomials,
    u64 mod
) {
    RREF system;

    const std::size_t rows =
        data.size();

    const std::size_t cols =
        monomials.size();

    system.matrix.assign(
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
            system.matrix[r][c] =
                evaluate_monomial_mod(
                    data[r],
                    sums[r],
                    monomials[c],
                    mod
                );
        }
    }

    std::size_t pivot_row = 0;

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
            system.matrix[pivot][col] == 0
        ) {
            ++pivot;
        }

        if (pivot == rows) {
            continue;
        }

        std::swap(
            system.matrix[pivot],
            system.matrix[pivot_row]
        );

        const u64 inverse =
            mod_inverse(
                system.matrix[pivot_row][col],
                mod
            );

        for (
            std::size_t j = col;
            j < cols;
            ++j
        ) {
            system.matrix[pivot_row][j] =
                mod_mul(
                    system.matrix[pivot_row][j],
                    inverse,
                    mod
                );
        }

        for (
            std::size_t r = 0;
            r < rows;
            ++r
        ) {
            if (r == pivot_row) {
                continue;
            }

            const u64 factor =
                system.matrix[r][col];

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
                        system.matrix[pivot_row][j],
                        mod
                    );

                if (
                    system.matrix[r][j] >=
                    product
                ) {
                    system.matrix[r][j] -=
                        product;
                } else {
                    system.matrix[r][j] =
                        mod -
                        (
                            product -
                            system.matrix[r][j]
                        );
                }
            }
        }

        system.pivot_columns.push_back(
            static_cast<int>(col)
        );

        ++pivot_row;
    }

    system.rank =
        pivot_row;

    system.is_pivot.assign(
        cols,
        false
    );

    for (
        int column :
        system.pivot_columns
    ) {
        system.is_pivot[
            static_cast<std::size_t>(
                column
            )
        ] = true;
    }

    return system;
}

/* =========================================================
 * Observable-only monomials
 * ========================================================= */

static std::vector<Monomial>
remove_S_monomials(
    const std::vector<Monomial>& all
) {
    std::vector<Monomial> result;

    for (
        const Monomial& monomial :
        all
    ) {
        if (
            monomial.exponent[
                OBSERVABLE_COUNT
            ] == 0
        ) {
            result.push_back(
                monomial
            );
        }
    }

    return result;
}

/* =========================================================
 * Compare nullities
 * ========================================================= */

static bool relation_involves_S(
    const RREF& full,
    const RREF& observable,
    std::size_t full_columns,
    std::size_t observable_columns
) {
    const std::size_t full_nullity =
        full_columns -
        full.rank;

    const std::size_t observable_nullity =
        observable_columns -
        observable.rank;

    std::cout
        << "FULL_RANK="
        << full.rank
        << " FULL_NULLITY="
        << full_nullity
        << "\n";

    std::cout
        << "OBS_RANK="
        << observable.rank
        << " OBS_NULLITY="
        << observable_nullity
        << "\n";

    return
        full_nullity >
        observable_nullity;
}

/* =========================================================
 * Extract null relation involving S
 * ========================================================= */

static bool find_relation_vector(
    const RREF& system,
    const std::vector<Monomial>& monomials,
    u64 mod,
    std::vector<u64>& relation
) {
    const std::size_t cols =
        monomials.size();

    std::vector<std::size_t> free_columns;

    for (
        std::size_t c = 0;
        c < cols;
        ++c
    ) {
        if (
            !system.is_pivot[c]
        ) {
            free_columns.push_back(c);
        }
    }

    /*
     * Try every free variable containing S.
     */
    for (
        const std::size_t free_column :
        free_columns
    ) {
        if (
            monomials[free_column].exponent[
                OBSERVABLE_COUNT
            ] == 0
        ) {
            continue;
        }

        std::vector<u64> candidate(
            cols,
            0
        );

        candidate[
            free_column
        ] = 1;

        /*
         * x_pivot =
         * -sum(a_j x_j)
         *
         * with one free variable set to 1.
         */
        for (
            std::size_t row = 0;
            row <
            system.pivot_columns.size();
            ++row
        ) {
            const std::size_t pivot_column =
                static_cast<std::size_t>(
                    system.pivot_columns[row]
                );

            const u64 coefficient =
                system.matrix[row][
                    free_column
                ];

            candidate[
                pivot_column
            ] =
                coefficient == 0
                    ? 0
                    : (
                        mod -
                        coefficient
                    ) % mod;
        }

        relation =
            std::move(candidate);

        return true;
    }

    /*
     * It is possible for every free-variable basis
     * vector to have no S coefficient even though the
     * nullspace difference was detected. In that case
     * try a few combinations of free variables.
     */
    if (
        free_columns.size() >= 2
    ) {
        const std::size_t a =
            free_columns[0];

        const std::size_t b =
            free_columns[1];

        std::vector<u64> candidate(
            cols,
            0
        );

        candidate[a] = 1;
        candidate[b] = 1;

        for (
            std::size_t row = 0;
            row <
            system.pivot_columns.size();
            ++row
        ) {
            const std::size_t pivot_column =
                static_cast<std::size_t>(
                    system.pivot_columns[row]
                );

            u64 value = 0;

            const u64 ca =
                system.matrix[row][a];

            const u64 cb =
                system.matrix[row][b];

            value =
                (
                    mod -
                    ca
                ) % mod;

            const u64 add_b =
                (
                    mod -
                    cb
                ) % mod;

            value += add_b;

            if (value >= mod) {
                value -= mod;
            }

            candidate[pivot_column] =
                value;
        }

        for (
            std::size_t i = 0;
            i < cols;
            ++i
        ) {
            if (
                candidate[i] != 0 &&
                monomials[i].exponent[
                    OBSERVABLE_COUNT
                ] > 0
            ) {
                relation =
                    std::move(candidate);

                return true;
            }
        }
    }

    return false;
}

/* =========================================================
 * Lift modular coefficients to mpz_class
 *
 * The conversion is performed through decimal strings,
 * avoiding all GMP signed-long-long constructor issues.
 * ========================================================= */

static std::vector<mpz_class>
lift_relation(
    const std::vector<u64>& coefficients,
    u64 mod
) {
    std::vector<mpz_class> result;

    result.reserve(
        coefficients.size()
    );

    const u64 half =
        mod / 2;

    for (u64 value : coefficients) {
        if (value == 0) {
            result.emplace_back(
                mpz_class(0)
            );
        } else if (value <= half) {
            result.emplace_back(
                mpz_from_u64(
                    value
                )
            );
        } else {
            const u64 magnitude =
                mod - value;

            const mpz_class positive =
                mpz_from_u64(
                    magnitude
                );

            result.emplace_back(
                -positive
            );
        }
    }

    /*
     * Divide by the gcd.
     */
    mpz_class g(0);

    for (
        const mpz_class& value :
        result
    ) {
        if (value == 0) {
            continue;
        }

        if (g == 0) {
            g = abs(value);
        } else {
            g = gcd(
                g,
                abs(value)
            );
        }
    }

    if (g != 0) {
        for (
            mpz_class& value :
            result
        ) {
            value /= g;
        }
    }

    /*
     * Canonical sign.
     */
    for (
        const mpz_class& value :
        result
    ) {
        if (value == 0) {
            continue;
        }

        if (value < 0) {
            for (
                mpz_class& x :
                result
            ) {
                x = -x;
            }
        }

        break;
    }

    return result;
}

/* =========================================================
 * Exact relation evaluation
 * ========================================================= */

static mpz_class evaluate_exact_relation(
    const Features& features,
    u64 S,
    const std::vector<Monomial>& monomials,
    const std::vector<mpz_class>& coefficients
) {
    mpz_class result(0);

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
                S,
                monomials[i]
            );

        result +=
            coefficients[i] *
            monomial;
    }

    return result;
}

/* =========================================================
 * Variable names
 * ========================================================= */

static std::string variable_name(
    int index
) {
    if (index == 0) {
        return "N";
    }

    if (
        index == OBSERVABLE_COUNT
    ) {
        return "S";
    }

    const int local =
        index - 1;

    const int base_index =
        local / 3;

    const int component =
        local % 3;

    const std::array<
        const char*,
        3
    > names = {
        "C",
        "B",
        "G"
    };

    return
        std::string(
            names[
                static_cast<std::size_t>(
                    component
                )
            ]
        ) +
        std::to_string(
            BASES[
                static_cast<std::size_t>(
                    base_index
                )
            ]
        );
}

/* =========================================================
 * Monomial pretty printing
 * ========================================================= */

static std::string monomial_name(
    const Monomial& monomial
) {
    std::string result;

    for (
        int i = 0;
        i < TOTAL_VARIABLES;
        ++i
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<std::size_t>(
                    i
                )
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
 * Print relation
 * ========================================================= */

static void print_relation(
    const std::vector<Monomial>& monomials,
    const std::vector<mpz_class>& coefficients
) {
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

/* =========================================================
 * Search one degree / modulus
 * ========================================================= */

static Relation search_degree_modulus(
    int degree,
    u64 mod,
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

    const auto observable_monomials =
        remove_S_monomials(
            monomials
        );

    std::cout
        << "MOD="
        << mod
        << "\n";

    std::cout
        << "ALL_MONOMIALS="
        << monomials.size()
        << "\n";

    std::cout
        << "OBS_MONOMIALS="
        << observable_monomials.size()
        << "\n";

    const RREF full =
        build_rref(
            training,
            training_sums,
            monomials,
            mod
        );

    const RREF observable =
        build_rref(
            training,
            training_sums,
            observable_monomials,
            mod
        );

    if (
        !relation_involves_S(
            full,
            observable,
            monomials.size(),
            observable_monomials.size()
        )
    ) {
        std::cout
            << "NO_S_RELATION_MODULAR\n";

        return relation;
    }

    std::cout
        << "S_RELATION_MODULAR=YES\n";

    std::vector<u64> modular_relation;

    if (
        !find_relation_vector(
            full,
            monomials,
            mod,
            modular_relation
        )
    ) {
        std::cout
            << "RELATION_EXTRACTION_FAILED\n";

        return relation;
    }

    const auto lifted =
        lift_relation(
            modular_relation,
            mod
        );

    /*
     * Exact training verification.
     */
    for (
        std::size_t i = 0;
        i < training.size();
        ++i
    ) {
        const mpz_class value =
            evaluate_exact_relation(
                training[i],
                training_sums[i],
                monomials,
                lifted
            );

        if (value != 0) {
            ++relation.training_failures;

            if (
                relation.training_failures <= 3
            ) {
                std::cout
                    << "TRAINING_FAILURE "
                    << "index="
                    << i
                    << " value="
                    << value.get_str()
                    << "\n";
            }
        }
    }

    if (
        relation.training_failures != 0
    ) {
        std::cout
            << "MODULAR_RELATION_NOT_EXACT\n";

        return relation;
    }

    relation.exact_training =
        true;

    /*
     * Exact validation.
     */
    for (
        std::size_t i = 0;
        i < validation.size();
        ++i
    ) {
        const mpz_class value =
            evaluate_exact_relation(
                validation[i],
                validation_sums[i],
                monomials,
                lifted
            );

        if (value != 0) {
            ++relation.validation_failures;

            if (
                relation.validation_failures <= 5
            ) {
                std::cout
                    << "VALIDATION_FAILURE "
                    << "index="
                    << i
                    << " value="
                    << value.get_str()
                    << "\n";
            }
        }
    }

    if (
        relation.validation_failures == 0
    ) {
        relation.found = true;
        relation.exact_validation =
            true;

        relation.monomials =
            monomials;

        relation.coefficients =
            lifted;

        std::cout
            << "VALIDATED_RELATION=YES\n";

        print_relation(
            monomials,
            lifted
        );
    } else {
        std::cout
            << "VALIDATED_RELATION=NO\n";
    }

    return relation;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 348\n"
        << "LOCAL RECURRENCE INVARIANT SEARCH\n"
        << "CAN p+q SATISFY A LOW-DEGREE RELATION WITH "
           "CROSS-BASE RECURRENCE COEFFICIENT DATA?\n\n";

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
        0x348348348ULL
    );

    constexpr std::size_t TOTAL_CASES =
        1200;

    constexpr std::size_t TRAINING_CASES =
        900;

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
                static_cast<u64>(small),
                static_cast<u64>(large),
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
     * Extract observables from N only.
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
                cases[i].S
            );
        } else {
            validation_features.push_back(
                features
            );

            validation_sums.push_back(
                cases[i].S
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
        << "OBSERVABLE_FEATURES="
        << OBSERVABLE_COUNT
        << "\n";

    std::cout
        << "FEATURES=N";

    for (int base : BASES) {
        std::cout
            << ",C"
            << base
            << ",B"
            << base
            << ",G"
            << base;
    }

    std::cout
        << ",S\n\n";

    /*
     * --------------------------------------------------
     * Degree 1 and 2.
     * --------------------------------------------------
     */

    Relation validated;

    for (
        int degree = 1;
        degree <= 2;
        ++degree
    ) {
        std::cout
            << "============================\n"
            << "DEGREE="
            << degree
            << "\n";

        /*
         * First modulus.
         */
        const Relation first =
            search_degree_modulus(
                degree,
                1000000007ULL,
                training_features,
                training_sums,
                validation_features,
                validation_sums
            );

        if (
            first.found
        ) {
            validated =
                first;

            break;
        }

        /*
         * Second independent modulus.
         */
        const Relation second =
            search_degree_modulus(
                degree,
                1000000009ULL,
                training_features,
                training_sums,
                validation_features,
                validation_sums
            );

        if (
            second.found
        ) {
            validated =
                second;

            break;
        }
    }

    /*
     * --------------------------------------------------
     * Final result
     * --------------------------------------------------
     */

    std::cout
        << "\n"
        << "============================\n";

    if (
        validated.found &&
        validated.exact_training &&
        validated.exact_validation
    ) {
        std::cout
            << "DISCOVERY_STATUS="
               "VALIDATED_RELATION\n";

        std::cout
            << "DISCOVERY_DEGREE="
            << validated.degree
            << "\n";
    } else {
        std::cout
            << "DISCOVERY_STATUS="
               "NO_VALIDATED_RELATION_DEGREE_1_TO_2\n";
    }

    std::cout
        << "============================\n"
        << "FINISHED EXPERIMENT 348\n";

    return 0;
}