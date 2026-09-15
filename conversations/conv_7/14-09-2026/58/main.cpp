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

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

static constexpr int OBSERVABLE_COUNT = 25;
// 0 = N
// 1..4   = A0..A3 base 2
// 5..8   = A0..A3 base 3
// 9..12  = A0..A3 base 5
// 13..16 = A0..A3 base 7
// 17..20 = A0..A3 base 11
// 21..24 = A0..A3 base 13
//
// Variable 25 = S = p+q

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
    std::array<int, OBSERVABLE_COUNT + 1> exponent{};
};

struct ModularSystem {
    std::vector<std::vector<u64>> matrix;
    std::vector<int> pivot_column_for_row;
    std::vector<bool> is_pivot;
    std::size_t rank = 0;
};

struct RelationCandidate {
    bool found = false;
    bool exact_training = false;
    bool exact_validation = false;

    int degree = -1;

    std::vector<Monomial> monomials;
    std::vector<mpz_class> coefficients;

    u64 training_failures = 0;
    u64 validation_failures = 0;
};

struct Stats {
    u64 tested = 0;
    u64 validated = 0;
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
 * Prime sieve
 * ========================================================= */

static std::vector<int> sieve_primes(
    int limit
) {
    std::vector<bool> composite(
        static_cast<std::size_t>(limit + 1),
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (
            composite[
                static_cast<std::size_t>(i)
            ]
        ) {
            continue;
        }

        primes.push_back(i);

        if (
            static_cast<long long>(i) * i <= limit
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
 * Kummer histogram
 *
 * Returns A_k =
 *
 * #{ t in [0,N] :
 *    v_p(C(N,t)) = k }.
 *
 * Only N is used.
 * ========================================================= */

static std::vector<u64> valuation_histogram(
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
    std::vector<std::vector<u64>> dp(
        2,
        std::vector<u64>(
            static_cast<std::size_t>(r + 1),
            0
        )
    );

    std::vector<std::vector<u64>> next(
        2,
        std::vector<u64>(
            static_cast<std::size_t>(r + 1),
            0
        )
    );

    dp[0][0] = 1;

    for (int i = 0; i < r; ++i) {
        for (auto& row : next) {
            std::fill(
                row.begin(),
                row.end(),
                0
            );
        }

        const u64 d =
            digits[
                static_cast<std::size_t>(i)
            ];

        /*
         * Previous borrow = 0.
         */
        for (
            int k = 0;
            k <= i;
            ++k
        ) {
            const u64 count =
                dp[0][
                    static_cast<std::size_t>(k)
                ];

            if (count == 0) {
                continue;
            }

            /*
             * No borrow.
             */
            next[0][
                static_cast<std::size_t>(k)
            ] +=
                count * (d + 1);

            /*
             * Borrow begins.
             */
            next[1][
                static_cast<std::size_t>(k + 1)
            ] +=
                count *
                (
                    static_cast<u64>(p)
                    - d
                    - 1
                );
        }

        /*
         * Previous borrow = 1.
         */
        for (
            int k = 0;
            k <= i;
            ++k
        ) {
            const u64 count =
                dp[1][
                    static_cast<std::size_t>(k)
                ];

            if (count == 0) {
                continue;
            }

            /*
             * Borrow resolves.
             */
            next[0][
                static_cast<std::size_t>(k)
            ] +=
                count * d;

            /*
             * Borrow continues.
             */
            next[1][
                static_cast<std::size_t>(k + 1)
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

    /*
     * Feature 0 = N.
     */
    result.x[0] = N;

    int index = 1;

    /*
     * Four Kummer coefficients for each base.
     */
    for (int base : BASES) {
        const auto histogram =
            valuation_histogram(
                N,
                base
            );

        for (int k = 0; k <= 3; ++k) {
            result.x[
                static_cast<std::size_t>(index++)
            ] =
                histogram.size() >
                static_cast<std::size_t>(k)
                    ? histogram[
                        static_cast<std::size_t>(k)
                    ]
                    : 0;
        }
    }

    return result;
}

/* =========================================================
 * Generate monomials total degree <= D
 *
 * Variables 0..24 = observables
 * Variable 25     = S
 * ========================================================= */

static void generate_monomials_recursive(
    int variable,
    int remaining_degree,
    std::array<int, OBSERVABLE_COUNT + 1>& current,
    std::vector<Monomial>& result
) {
    if (
        variable ==
        OBSERVABLE_COUNT + 1
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
            static_cast<std::size_t>(variable)
        ] = exponent;

        generate_monomials_recursive(
            variable + 1,
            remaining_degree - exponent,
            current,
            result
        );
    }

    current[
        static_cast<std::size_t>(variable)
    ] = 0;
}

static std::vector<Monomial> generate_monomials(
    int degree
) {
    std::vector<Monomial> result;

    std::array<int, OBSERVABLE_COUNT + 1> current{};

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
        variable < OBSERVABLE_COUNT + 1;
        ++variable
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<
                    std::size_t
                >(variable)
            ];

        if (exponent == 0) {
            continue;
        }

        const u64 value =
            variable < OBSERVABLE_COUNT
                ? features.x[
                    static_cast<
                        std::size_t
                    >(variable)
                ]
                : S;

        result =
            mod_mul(
                result,
                mod_pow(
                    value % mod,
                    static_cast<u64>(exponent),
                    mod
                ),
                mod
            );
    }

    return result;
}

/* =========================================================
 * Exact monomial evaluation
 * ========================================================= */

static mpz_class evaluate_monomial_exact(
    const Features& features,
    u64 S,
    const Monomial& monomial
) {
    mpz_class result(1);

    for (
        int variable = 0;
        variable < OBSERVABLE_COUNT + 1;
        ++variable
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<
                    std::size_t
                >(variable)
            ];

        if (exponent == 0) {
            continue;
        }

        const u64 value =
            variable < OBSERVABLE_COUNT
                ? features.x[
                    static_cast<
                        std::size_t
                    >(variable)
                ]
                : S;

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
 * Build modular RREF
 * ========================================================= */

static ModularSystem build_rref(
    const std::vector<Features>& data,
    const std::vector<u64>& sums,
    const std::vector<Monomial>& monomials,
    u64 mod
) {
    const std::size_t rows =
        data.size();

    const std::size_t cols =
        monomials.size();

    ModularSystem system;

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

    system.pivot_column_for_row.clear();

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
                const u64 value =
                    mod_mul(
                        factor,
                        system.matrix[pivot_row][j],
                        mod
                    );

                if (
                    system.matrix[r][j] >= value
                ) {
                    system.matrix[r][j] -= value;
                } else {
                    system.matrix[r][j] =
                        mod -
                        (
                            value -
                            system.matrix[r][j]
                        );
                }
            }
        }

        system.pivot_column_for_row.push_back(
            static_cast<int>(col)
        );

        ++pivot_row;
    }

    system.rank = pivot_row;

    system.is_pivot.assign(
        cols,
        false
    );

    for (
        int column :
        system.pivot_column_for_row
    ) {
        system.is_pivot[
            static_cast<
                std::size_t
            >(column)
        ] = true;
    }

    return system;
}

/* =========================================================
 * Determine whether there is a relation involving S
 *
 * A relation among observables only is not interesting.
 *
 * Let:
 *
 *   R_full = rank(all monomials)
 *   R_obs  = rank(monoms with S exponent = 0)
 *
 * A genuine relation involving S exists precisely when
 *
 *   nullity(full) > nullity(observable-only)
 *
 * modulo the chosen prime.
 * ========================================================= */

static bool relation_involves_S(
    const std::vector<Features>& data,
    const std::vector<u64>& sums,
    const std::vector<Monomial>& all_monomials,
    const std::vector<Monomial>& observable_monomials,
    u64 mod
) {
    const ModularSystem full =
        build_rref(
            data,
            sums,
            all_monomials,
            mod
        );

    const ModularSystem observable =
        build_rref(
            data,
            sums,
            observable_monomials,
            mod
        );

    const std::size_t full_columns =
        all_monomials.size();

    const std::size_t observable_columns =
        observable_monomials.size();

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
 * Find one modular null relation with S
 *
 * The system must have a free variable corresponding to
 * a monomial containing S.
 * ========================================================= */

static bool find_modular_relation_with_S(
    const std::vector<Features>& data,
    const std::vector<u64>& sums,
    const std::vector<Monomial>& monomials,
    u64 mod,
    std::vector<u64>& relation
) {
    const ModularSystem system =
        build_rref(
            data,
            sums,
            monomials,
            mod
        );

    const std::size_t cols =
        monomials.size();

    /*
     * Locate free columns that contain S.
     */
    std::vector<std::size_t> free_S_columns;

    for (
        std::size_t c = 0;
        c < cols;
        ++c
    ) {
        if (
            system.is_pivot[c]
        ) {
            continue;
        }

        if (
            monomials[c].exponent[
                OBSERVABLE_COUNT
            ] > 0
        ) {
            free_S_columns.push_back(c);
        }
    }

    if (
        free_S_columns.empty()
    ) {
        return false;
    }

    /*
     * Select the first S-containing free column.
     */
    const std::size_t free_column =
        free_S_columns.front();

    relation.assign(
        cols,
        0
    );

    relation[
        free_column
    ] = 1;

    /*
     * RREF rows correspond to pivot columns.
     *
     * x_pivot +
     * sum_j a_j x_j = 0
     *
     * therefore
     *
     * x_pivot = -sum_j a_j x_j.
     */
    for (
        std::size_t row = 0;
        row <
        system.pivot_column_for_row.size();
        ++row
    ) {
        const std::size_t pivot_column =
            static_cast<
                std::size_t
            >(
                system.pivot_column_for_row[row]
            );

        const u64 coefficient =
            system.matrix[row][
                free_column
            ];

        relation[pivot_column] =
            coefficient == 0
                ? 0
                : (
                    mod -
                    coefficient
                ) % mod;
    }

    return true;
}

/* =========================================================
 * Convert modular relation to signed coefficients
 *
 * Each residue is represented in
 * [-mod/2, mod/2].
 * ========================================================= */

static std::vector<mpz_class>
lift_modular_relation(
    const std::vector<u64>& modular,
    u64 mod
) {
    std::vector<mpz_class> result;

    result.reserve(
        modular.size()
    );

    const u64 half =
        mod / 2;

    for (u64 value : modular) {
        if (value == 0) {
            result.push_back(
                mpz_class(0)
            );
        } else if (value <= half) {
            result.push_back(
                mpz_class(value)
            );
        } else {
            result.push_back(
                -mpz_class(
                    mod - value
                )
            );
        }
    }

    /*
     * Normalize by gcd.
     */
    mpz_class g = 0;

    for (
        const mpz_class& x :
        result
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

    if (g != 0) {
        for (mpz_class& x : result) {
            x /= g;
        }
    }

    /*
     * Canonical sign.
     */
    for (
        const mpz_class& x :
        result
    ) {
        if (x == 0) {
            continue;
        }

        if (x < 0) {
            for (mpz_class& y : result) {
                y = -y;
            }
        }

        break;
    }

    return result;
}

/* =========================================================
 * Exact relation evaluation
 * ========================================================= */

static mpz_class evaluate_relation(
    const Features& features,
    u64 S,
    const std::vector<Monomial>& monomials,
    const std::vector<mpz_class>& coefficients
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

    if (index == OBSERVABLE_COUNT) {
        return "S";
    }

    const int local =
        index - 1;

    const int base_index =
        local / 4;

    const int k =
        local % 4;

    return
        "B" +
        std::to_string(
            BASES[
                static_cast<
                    std::size_t
                >(base_index)
            ]
        ) +
        std::to_string(k);
}

/* =========================================================
 * Monomial string
 * ========================================================= */

static std::string monomial_name(
    const Monomial& monomial
) {
    std::string result;

    for (
        int i = 0;
        i < OBSERVABLE_COUNT + 1;
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
 * Main search for one degree/modulus
 * ========================================================= */

static RelationCandidate search_degree_modulus(
    int degree,
    u64 mod,
    const std::vector<Features>& training,
    const std::vector<u64>& training_sums,
    const std::vector<Features>& validation,
    const std::vector<u64>& validation_sums
) {
    RelationCandidate candidate;

    candidate.degree =
        degree;

    const auto all_monomials =
        generate_monomials(
            degree
        );

    /*
     * Extract observable-only monomials.
     */
    std::vector<Monomial> observable_monomials;

    for (
        const Monomial& monomial :
        all_monomials
    ) {
        if (
            monomial.exponent[
                OBSERVABLE_COUNT
            ] == 0
        ) {
            observable_monomials.push_back(
                monomial
            );
        }
    }

    std::cout
        << "MOD="
        << mod
        << "\n";

    std::cout
        << "ALL_MONOMIALS="
        << all_monomials.size()
        << "\n";

    std::cout
        << "OBS_MONOMIALS="
        << observable_monomials.size()
        << "\n";

    if (
        !relation_involves_S(
            training,
            training_sums,
            all_monomials,
            observable_monomials,
            mod
        )
    ) {
        std::cout
            << "NO_S_RELATION_MODULAR\n";

        return candidate;
    }

    std::cout
        << "S_RELATION_MODULAR=YES\n";

    std::vector<u64> modular_relation;

    if (
        !find_modular_relation_with_S(
            training,
            training_sums,
            all_monomials,
            mod,
            modular_relation
        )
    ) {
        std::cout
            << "FAILED_TO_EXTRACT_RELATION\n";

        return candidate;
    }

    const auto integer_coefficients =
        lift_modular_relation(
            modular_relation,
            mod
        );

    /*
     * Exact training verification.
     */
    u64 training_failures = 0;

    for (
        std::size_t i = 0;
        i < training.size();
        ++i
    ) {
        const mpz_class value =
            evaluate_relation(
                training[i],
                training_sums[i],
                all_monomials,
                integer_coefficients
            );

        if (
            value != 0
        ) {
            ++training_failures;

            if (training_failures <= 3) {
                std::cout
                    << "EXACT_TRAINING_FAILURE "
                    << "index="
                    << i
                    << " value="
                    << value.get_str()
                    << "\n";
            }
        }
    }

    candidate.monomials =
        all_monomials;

    candidate.coefficients =
        integer_coefficients;

    candidate.training_failures =
        training_failures;

    if (
        training_failures != 0
    ) {
        std::cout
            << "MODULAR_RELATION_BUT_NOT_EXACT\n";

        return candidate;
    }

    candidate.exact_training = true;

    /*
     * Exact validation.
     */
    u64 validation_failures = 0;

    for (
        std::size_t i = 0;
        i < validation.size();
        ++i
    ) {
        const mpz_class value =
            evaluate_relation(
                validation[i],
                validation_sums[i],
                all_monomials,
                integer_coefficients
            );

        if (
            value != 0
        ) {
            ++validation_failures;

            if (validation_failures <= 5) {
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

    candidate.validation_failures =
        validation_failures;

    if (
        validation_failures == 0
    ) {
        candidate.found = true;
        candidate.exact_validation = true;

        std::cout
            << "VALIDATED_RELATION=YES\n";

        std::size_t nonzero = 0;

        for (
            const mpz_class& coefficient :
            integer_coefficients
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
            i < integer_coefficients.size();
            ++i
        ) {
            if (
                integer_coefficients[i] == 0
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
                << integer_coefficients[i].get_str()
                << ")*"
                << monomial_name(
                    all_monomials[i]
                );
        }

        std::cout
            << "\n";
    } else {
        std::cout
            << "VALIDATED_RELATION=NO\n";
    }

    return candidate;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 347\n"
        << "FULL KUMMER PROFILE ALGEBRAIC SEARCH\n"
        << "CAN p+q SATISFY A LOW-DEGREE RELATION "
           "WITH MULTI-BASE KUMMER COEFFICIENTS?\n\n";

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
            usable_primes.push_back(p);
        }
    }

    /*
     * --------------------------------------------------
     * Dataset
     * --------------------------------------------------
     */

    std::mt19937_64 rng(
        0x347347347ULL
    );

    constexpr std::size_t TOTAL_CASES =
        1400;

    constexpr std::size_t TRAINING_CASES =
        1000;

    constexpr std::size_t VALIDATION_CASES =
        400;

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
            static_cast<u64>(small) *
            static_cast<u64>(large);

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
     * Extract features from N only.
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
        << "FEATURES="
        << "N";

    for (int base : BASES) {
        std::cout
            << ",B"
            << base
            << "0"
            << ",B"
            << base
            << "1"
            << ",B"
            << base
            << "2"
            << ",B"
            << base
            << "3";
    }

    std::cout
        << ",S\n\n";

    /*
     * --------------------------------------------------
     * Search degree 1 and 2.
     *
     * Degree 3 is intentionally not attempted here:
     * with 26 variables it has 3654 monomials and would
     * turn a clean exploratory experiment into a much
     * heavier linear-algebra calculation.
     * --------------------------------------------------
     */

    RelationCandidate validated;

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
         * First modular test.
         */
        const RelationCandidate first =
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
         *
         * We only proceed to accept the result if a
         * relation survives both exact checks.
         */
        const RelationCandidate second =
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
     * Result
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
            << "DISCOVERY_STATUS=VALIDATED_RELATION\n";

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
        << "FINISHED EXPERIMENT 347\n";

    return 0;
}
