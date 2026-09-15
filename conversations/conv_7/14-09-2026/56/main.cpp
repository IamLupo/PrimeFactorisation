#include <array>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

using Poly = std::vector<u64>;

struct SemiprimeCase {
    u64 N;
    u64 p;
    u64 q;
    u64 sum;
};

struct Features {
    /*
     * [0] = N
     *
     * Base 2:
     * [1] = A0
     * [2] = A1
     * [3] = A2
     *
     * Base 3:
     * [4] = A0
     * [5] = A1
     * [6] = A2
     *
     * Base 5:
     * [7] = A0
     * [8] = A1
     * [9] = A2
     */
    std::array<u64, 10> x{};
};

struct Monomial {
    std::array<int, 10> exponent{};
};

struct Candidate {
    bool found = false;
    bool unique = false;

    int degree = -1;

    std::vector<Monomial> monomials;
    std::vector<mpq_class> coefficients;

    u64 training_failures = 0;
    u64 validation_failures = 0;
};

struct Stats {
    u64 total_cases = 0;
    u64 training_cases = 0;
    u64 validation_cases = 0;
};

/* =========================================================
 * Constants
 * ========================================================= */

static const std::array<int, 3> BASES = {
    2, 3, 5
};

/* =========================================================
 * GMP conversion
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
        static_cast<std::size_t>(limit + 1),
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
 * Kummer valuation histogram
 *
 * Returns
 *
 * A_k =
 * #{ t in [0,N] :
 *    v_p(C(N,t)) = k }.
 *
 * The input is N only.
 * ========================================================= */

static Poly valuation_histogram(
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
    std::vector<Poly> dp(
        2,
        Poly(
            static_cast<std::size_t>(r + 1),
            0
        )
    );

    std::vector<Poly> next(
        2,
        Poly(
            static_cast<std::size_t>(r + 1),
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
         * Previous borrow = 0.
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
             * No borrow.
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
         * Previous borrow = 1.
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
        const Poly histogram =
            valuation_histogram(
                N,
                base
            );

        /*
         * A_0.
         */
        result.x[
            static_cast<
                std::size_t
            >(index++)
        ] =
            histogram.size() > 0
                ? histogram[0]
                : 0;

        /*
         * A_1.
         */
        result.x[
            static_cast<
                std::size_t
            >(index++)
        ] =
            histogram.size() > 1
                ? histogram[1]
                : 0;

        /*
         * A_2.
         */
        result.x[
            static_cast<
                std::size_t
            >(index++)
        ] =
            histogram.size() > 2
                ? histogram[2]
                : 0;
    }

    return result;
}

/* =========================================================
 * Monomial generation
 * ========================================================= */

static void generate_monomials_recursive(
    int variable,
    int remaining_degree,
    std::array<int, 10>& current,
    std::vector<Monomial>& result
) {
    if (
        variable == 10
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

static std::vector<Monomial>
generate_monomials(
    int degree
) {
    std::vector<Monomial> result;

    std::array<int, 10> current{};

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
    u64 a,
    u64 mod
) {
    return mod_pow(
        a,
        mod - 2,
        mod
    );
}

/* =========================================================
 * Monomial evaluation
 *
 * u128 is used internally to avoid accidental intermediate
 * overflow.
 * ========================================================= */

static u64 evaluate_monomial(
    const Features& features,
    const Monomial& monomial
) {
    u128 result = 1;

    for (
        int i = 0;
        i < 10;
        ++i
    ) {
        const int exponent =
            monomial.exponent[
                static_cast<
                    std::size_t
                >(i)
            ];

        for (
            int j = 0;
            j < exponent;
            ++j
        ) {
            result *=
                static_cast<u128>(
                    features.x[
                        static_cast<
                            std::size_t
                        >(i)
                    ]
                );
        }
    }

    if (
        result >
        static_cast<u128>(
            UINT64_MAX
        )
    ) {
        std::cerr
            << "FATAL: monomial overflow\n";

        std::exit(1);
    }

    return static_cast<u64>(
        result
    );
}

static u64 evaluate_monomial_mod(
    const Features& features,
    const Monomial& monomial,
    u64 mod
) {
    u64 result = 1;

    for (
        int i = 0;
        i < 10;
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

        result =
            mod_mul(
                result,
                mod_pow(
                    features.x[
                        static_cast<
                            std::size_t
                        >(i)
                    ] % mod,
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
 * Modular linear system
 *
 * 0 = inconsistent
 * 1 = consistent but underdetermined
 * 2 = consistent and unique
 * ========================================================= */

static int modular_consistency(
    const std::vector<Features>& data,
    const std::vector<u64>& labels,
    const std::vector<Monomial>& monomials,
    u64 mod
) {
    const std::size_t rows =
        data.size();

    const std::size_t cols =
        monomials.size();

    std::vector<std::vector<u64>> matrix(
        rows,
        std::vector<u64>(
            cols + 1,
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
                    monomials[c],
                    mod
                );
        }

        matrix[r][cols] =
            labels[r] % mod;
    }

    std::size_t pivot_row = 0;
    std::size_t rank = 0;

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
            j <= cols;
            ++j
        ) {
            matrix[pivot_row][j] =
                mod_mul(
                    matrix[pivot_row][j],
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
                matrix[r][col];

            if (factor == 0) {
                continue;
            }

            for (
                std::size_t j = col;
                j <= cols;
                ++j
            ) {
                const u64 product =
                    mod_mul(
                        factor,
                        matrix[pivot_row][j],
                        mod
                    );

                if (
                    matrix[r][j] >= product
                ) {
                    matrix[r][j] -= product;
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
        ++rank;
    }

    /*
     * Inconsistency check.
     */
    for (
        std::size_t r = 0;
        r < rows;
        ++r
    ) {
        bool all_zero = true;

        for (
            std::size_t c = 0;
            c < cols;
            ++c
        ) {
            if (matrix[r][c] != 0) {
                all_zero = false;
                break;
            }
        }

        if (
            all_zero &&
            matrix[r][cols] != 0
        ) {
            return 0;
        }
    }

    if (rank < cols) {
        return 1;
    }

    return 2;
}

/* =========================================================
 * Exact rational solve
 * ========================================================= */

static bool exact_solve(
    const std::vector<Features>& data,
    const std::vector<u64>& labels,
    const std::vector<Monomial>& monomials,
    std::vector<mpq_class>& solution
) {
    const std::size_t rows =
        data.size();

    const std::size_t cols =
        monomials.size();

    std::vector<std::vector<mpq_class>> matrix(
        rows,
        std::vector<mpq_class>(
            cols + 1,
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
                mpq_from_u64(
                    evaluate_monomial(
                        data[r],
                        monomials[c]
                    )
                );
        }

        matrix[r][cols] =
            mpq_from_u64(
                labels[r]
            );
    }

    std::size_t pivot_row = 0;

    std::vector<int> pivot_for_column(
        cols,
        -1
    );

    for (
        std::size_t col = 0;
        col < cols;
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
            j <= cols;
            ++j
        ) {
            matrix[pivot_row][j] /=
                pivot_value;
        }

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
                j <= cols;
                ++j
            ) {
                matrix[r][j] -=
                    factor *
                    matrix[pivot_row][j];
            }
        }

        pivot_for_column[col] =
            static_cast<int>(
                pivot_row
            );

        ++pivot_row;

        if (
            pivot_row == rows
        ) {
            break;
        }
    }

    /*
     * Require unique solution.
     */
    if (
        pivot_row < cols
    ) {
        return false;
    }

    /*
     * Consistency.
     */
    for (
        std::size_t r = 0;
        r < rows;
        ++r
    ) {
        bool all_zero = true;

        for (
            std::size_t c = 0;
            c < cols;
            ++c
        ) {
            if (matrix[r][c] != 0) {
                all_zero = false;
                break;
            }
        }

        if (
            all_zero &&
            matrix[r][cols] != 0
        ) {
            return false;
        }
    }

    solution.resize(
        cols
    );

    for (
        std::size_t c = 0;
        c < cols;
        ++c
    ) {
        const int row =
            pivot_for_column[c];

        if (row < 0) {
            return false;
        }

        solution[c] =
            matrix[
                static_cast<
                    std::size_t
                >(row)
            ][cols];
    }

    return true;
}

/* =========================================================
 * Exact polynomial evaluation
 * ========================================================= */

static mpq_class evaluate_polynomial(
    const Features& features,
    const std::vector<Monomial>& monomials,
    const std::vector<mpq_class>& coefficients
) {
    mpq_class result(0);

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

        result +=
            coefficients[i] *
            mpq_from_u64(
                evaluate_monomial(
                    features,
                    monomials[i]
                )
            );
    }

    return result;
}

/* =========================================================
 * Polynomial pretty-printing
 * ========================================================= */

static std::string monomial_name(
    const Monomial& monomial
) {
    static const std::array<
        const char*,
        10
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
        "B52"
    };

    std::string result;

    for (
        int i = 0;
        i < 10;
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
            names[
                static_cast<
                    std::size_t
                >(i)
            ];

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

static void print_candidate(
    const Candidate& candidate
) {
    if (!candidate.found) {
        return;
    }

    std::size_t nonzero = 0;

    for (
        const auto& coefficient :
        candidate.coefficients
    ) {
        if (
            coefficient != 0
        ) {
            ++nonzero;
        }
    }

    std::cout
        << "CANDIDATE_DEGREE="
        << candidate.degree
        << "\n";

    std::cout
        << "NONZERO_TERMS="
        << nonzero
        << "\n";

    std::cout
        << "POLYNOMIAL:\n";

    bool first = true;

    for (
        std::size_t i = 0;
        i < candidate.coefficients.size();
        ++i
    ) {
        const auto& coefficient =
            candidate.coefficients[i];

        if (coefficient == 0) {
            continue;
        }

        if (!first) {
            std::cout << "\n";
        }

        first = false;

        std::cout
            << "("
            << coefficient.get_str()
            << ")*"
            << monomial_name(
                candidate.monomials[i]
            );
    }

    std::cout << "\n";
}

/* =========================================================
 * Search one degree
 * ========================================================= */

static Candidate search_degree(
    int degree,
    const std::vector<Features>& training,
    const std::vector<u64>& training_labels,
    const std::vector<Features>& validation,
    const std::vector<u64>& validation_labels
) {
    Candidate candidate;

    candidate.degree =
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

    /*
     * Two independent modular checks.
     */
    constexpr u64 MOD1 =
        1000000007ULL;

    constexpr u64 MOD2 =
        1000000009ULL;

    const int status1 =
        modular_consistency(
            training,
            training_labels,
            monomials,
            MOD1
        );

    const int status2 =
        modular_consistency(
            training,
            training_labels,
            monomials,
            MOD2
        );

    std::cout
        << "MOD1_STATUS="
        << status1
        << " MOD2_STATUS="
        << status2
        << "\n";

    /*
     * 0 = no solution
     * 1 = nonunique
     * 2 = unique
     */
    if (
        status1 != 2 ||
        status2 != 2
    ) {
        std::cout
            << "EXACT_SOLVE_SKIPPED\n";

        return candidate;
    }

    std::vector<mpq_class> coefficients;

    if (
        !exact_solve(
            training,
            training_labels,
            monomials,
            coefficients
        )
    ) {
        std::cout
            << "EXACT_SOLVE_FAILED\n";

        return candidate;
    }

    candidate.found = true;
    candidate.unique = true;

    candidate.monomials =
        monomials;

    candidate.coefficients =
        coefficients;

    /*
     * Training verification.
     */
    for (
        std::size_t i = 0;
        i < training.size();
        ++i
    ) {
        const mpq_class predicted =
            evaluate_polynomial(
                training[i],
                monomials,
                coefficients
            );

        const mpq_class actual =
            mpq_from_u64(
                training_labels[i]
            );

        if (
            predicted != actual
        ) {
            ++candidate.training_failures;
        }
    }

    /*
     * Completely held-out validation.
     */
    for (
        std::size_t i = 0;
        i < validation.size();
        ++i
    ) {
        const mpq_class predicted =
            evaluate_polynomial(
                validation[i],
                monomials,
                coefficients
            );

        const mpq_class actual =
            mpq_from_u64(
                validation_labels[i]
            );

        if (
            predicted != actual
        ) {
            ++candidate.validation_failures;

            if (
                candidate.validation_failures <= 5
            ) {
                std::cout
                    << "VALIDATION_FAILURE "
                    << "index="
                    << i
                    << " predicted="
                    << predicted.get_str()
                    << " actual="
                    << actual.get_str()
                    << "\n";
            }
        }
    }

    std::cout
        << "TRAINING_FAILURES="
        << candidate.training_failures
        << "\n";

    std::cout
        << "VALIDATION_FAILURES="
        << candidate.validation_failures
        << "\n";

    if (
        candidate.validation_failures == 0
    ) {
        std::cout
            << "VALIDATED_CANDIDATE=YES\n";

        print_candidate(
            candidate
        );
    } else {
        std::cout
            << "VALIDATED_CANDIDATE=NO\n";
    }

    return candidate;
}

/* =========================================================
 * Main
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 345\n"
        << "RICH KUMMER FEATURE POLYNOMIAL SEARCH\n"
        << "CAN LOW-DEGREE POLYNOMIAL DATA DETERMINE p+q?\n\n";

    /*
     * --------------------------------------------------
     * Candidate prime pool
     * --------------------------------------------------
     */

    const auto all_primes =
        sieve_primes(997);

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
        0x345345345ULL
    );

    constexpr std::size_t TOTAL_CASES =
        900;

    constexpr std::size_t TRAINING_CASES =
        600;

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
            std::min(p, q);

        const int large =
            std::max(p, q);

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
     * Feature extraction
     *
     * Only N is given to extract_features().
     * --------------------------------------------------
     */

    std::vector<Features> training_features;
    std::vector<u64> training_labels;

    std::vector<Features> validation_features;
    std::vector<u64> validation_labels;

    training_features.reserve(
        TRAINING_CASES
    );

    training_labels.reserve(
        TRAINING_CASES
    );

    validation_features.reserve(
        VALIDATION_CASES
    );

    validation_labels.reserve(
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

            training_labels.push_back(
                cases[i].sum
            );
        } else {
            validation_features.push_back(
                features
            );

            validation_labels.push_back(
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
        << "FEATURE_COUNT=10\n";

    std::cout
        << "FEATURES="
        << "N,"
        << "B20,B21,B22,"
        << "B30,B31,B32,"
        << "B50,B51,B52"
        << "\n\n";

    /*
     * --------------------------------------------------
     * Polynomial search
     *
     * Degree is fixed independently.
     * No validation failure causes automatic degree
     * escalation within a candidate.
     * --------------------------------------------------
     */

    Candidate validated;

    for (
        int degree = 1;
        degree <= 3;
        ++degree
    ) {
        std::cout
            << "============================\n";

        const Candidate candidate =
            search_degree(
                degree,
                training_features,
                training_labels,
                validation_features,
                validation_labels
            );

        if (
            candidate.found &&
            candidate.validation_failures == 0
        ) {
            validated =
                candidate;

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
        validated.validation_failures == 0
    ) {
        std::cout
            << "DISCOVERY_STATUS=VALIDATED\n";

        std::cout
            << "DISCOVERY_DEGREE="
            << validated.degree
            << "\n";
    } else {
        std::cout
            << "DISCOVERY_STATUS="
               "NO_VALIDATED_POLYNOMIAL_DEGREE_1_TO_3\n";
    }

    std::cout
        << "============================\n"
        << "FINISHED EXPERIMENT 345\n";

    return 0;
}
