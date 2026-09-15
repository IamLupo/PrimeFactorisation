#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 S;
};

struct Candidate {
    std::string name;
    std::vector<u64> train_x;
    std::vector<u64> valid_x;
};

struct FitResult {
    bool possible = false;
    bool consistent = false;
    int rank = 0;
    std::vector<mpq_class> coeffs;
};

static mpz_class mpz_from_u64(u64 value) {
    return mpz_class(std::to_string(value));
}

static mpq_class mpq_from_u64(u64 value) {
    return mpq_class(mpz_from_u64(value));
}

static u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(limit + 1, false);
    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (!composite[i]) {
            primes.push_back(i);

            if (static_cast<long long>(i) * i <= limit) {
                for (int j = i * i; j <= limit; j += i) {
                    composite[j] = true;
                }
            }
        }
    }

    return primes;
}

static std::vector<int> digits_base(u64 n, int base) {
    std::vector<int> d;

    if (n == 0) {
        d.push_back(0);
        return d;
    }

    while (n > 0) {
        d.push_back(static_cast<int>(n % static_cast<u64>(base)));
        n /= static_cast<u64>(base);
    }

    return d;
}

static u64 digit_sum(u64 n, int base) {
    u64 result = 0;

    while (n > 0) {
        result += n % static_cast<u64>(base);
        n /= static_cast<u64>(base);
    }

    return result;
}

static u64 digit_zero_count(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);
    u64 result = 0;

    for (int x : d) {
        if (x == 0) {
            ++result;
        }
    }

    return result;
}

static u64 digit_nonzero_count(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);
    u64 result = 0;

    for (int x : d) {
        if (x != 0) {
            ++result;
        }
    }

    return result;
}

static u64 digit_square_sum(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);
    u64 result = 0;

    for (int x : d) {
        result += static_cast<u64>(x * x);
    }

    return result;
}

static u64 digit_product_plus_one(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);

    u64 result = 1;

    for (int x : d) {
        result *= static_cast<u64>(x + 1);
    }

    return result;
}

static u64 recurrence_C(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);

    u64 result = 0;

    for (std::size_t i = 1; i < d.size(); ++i) {
        result += static_cast<u64>(d[i - 1]) *
                  static_cast<u64>(d[i] + 1);
    }

    return result;
}

static u64 recurrence_B(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);

    u64 result = 0;
    const u64 b = static_cast<u64>(base);

    for (std::size_t i = 1; i < d.size(); ++i) {
        result += static_cast<u64>(d[i]) *
                  (b - static_cast<u64>(d[i - 1]));
    }

    return result;
}

static u64 recurrence_G(u64 n, int base) {
    std::vector<int> d = digits_base(n, base);

    if (d.size() < 2) {
        return 0;
    }

    const u64 b = static_cast<u64>(base);
    u64 result = static_cast<u64>(d[1]) * b;

    for (std::size_t i = 2; i < d.size(); ++i) {
        result += static_cast<u64>(d[i]) *
                  b *
                  static_cast<u64>(d[i - 2]);
    }

    return result;
}

static std::vector<u64> kummer_histogram_prefix(
    u64 m,
    int base,
    int max_k
) {
    std::vector<int> md = digits_base(m, base);

    struct State {
        mpz_class count[2][64];
    };

    std::vector<std::vector<mpz_class>> dp(
        2,
        std::vector<mpz_class>(static_cast<std::size_t>(max_k + 2))
    );

    dp[0][0] = 1;

    for (int digit : md) {
        std::vector<std::vector<mpz_class>> next(
            2,
            std::vector<mpz_class>(
                static_cast<std::size_t>(max_k + 2)
            )
        );

        const int d = digit;

        /*
         * Subtract t_i from m_i with incoming borrow b.
         *
         * From borrow 0:
         *   t_i <= d     -> no new borrow
         *   t_i >  d     -> new borrow
         *
         * From borrow 1:
         *   t_i <= d-1   -> no new borrow
         *   t_i >= d     -> new borrow
         */
        for (int k = 0; k <= max_k; ++k) {
            if (dp[0][k] != 0) {
                next[0][k] += dp[0][k] * (d + 1);

                if (k + 1 <= max_k) {
                    next[1][k + 1] +=
                        dp[0][k] * (base - d - 1);
                }
            }

            if (dp[1][k] != 0) {
                next[0][k] += dp[1][k] * d;

                if (k + 1 <= max_k) {
                    next[1][k + 1] +=
                        dp[1][k] * (base - d);
                }
            }
        }

        dp.swap(next);
    }

    std::vector<u64> result(static_cast<std::size_t>(max_k + 1), 0);

    /*
     * The valid t range is 0..m, so only final borrow=0 is allowed.
     */
    for (int k = 0; k <= max_k; ++k) {
        result[k] = dp[0][k].get_ui();
    }

    return result;
}

static u64 feature_sqrtN(const CaseData& c) {
    return isqrt_u64(c.N);
}

static u64 feature_D(const CaseData& c) {
    const u64 s = isqrt_u64(c.N);
    return c.N - s * s;
}

static void add_candidate(
    std::vector<Candidate>& candidates,
    const std::string& name,
    const std::vector<CaseData>& train_cases,
    const std::vector<CaseData>& valid_cases,
    u64 (*fn)(const CaseData&)
) {
    Candidate c;
    c.name = name;

    c.train_x.reserve(train_cases.size());
    c.valid_x.reserve(valid_cases.size());

    for (const CaseData& x : train_cases) {
        c.train_x.push_back(fn(x));
    }

    for (const CaseData& x : valid_cases) {
        c.valid_x.push_back(fn(x));
    }

    candidates.push_back(std::move(c));
}

static u64 make_feature_digit_sum_2(const CaseData& c) {
    return digit_sum(c.N, 2);
}

static u64 make_feature_digit_sum_3(const CaseData& c) {
    return digit_sum(c.N, 3);
}

static u64 make_feature_digit_sum_5(const CaseData& c) {
    return digit_sum(c.N, 5);
}

static u64 make_feature_digit_sum_7(const CaseData& c) {
    return digit_sum(c.N, 7);
}

static u64 make_feature_digit_sum_11(const CaseData& c) {
    return digit_sum(c.N, 11);
}

static u64 make_feature_digit_sum_13(const CaseData& c) {
    return digit_sum(c.N, 13);
}

static u64 make_feature_zero_2(const CaseData& c) {
    return digit_zero_count(c.N, 2);
}

static u64 make_feature_zero_3(const CaseData& c) {
    return digit_zero_count(c.N, 3);
}

static u64 make_feature_zero_5(const CaseData& c) {
    return digit_zero_count(c.N, 5);
}

static u64 make_feature_zero_7(const CaseData& c) {
    return digit_zero_count(c.N, 7);
}

static u64 make_feature_square_2(const CaseData& c) {
    return digit_square_sum(c.N, 2);
}

static u64 make_feature_square_3(const CaseData& c) {
    return digit_square_sum(c.N, 3);
}

static u64 make_feature_square_5(const CaseData& c) {
    return digit_square_sum(c.N, 5);
}

static u64 make_feature_prod_2(const CaseData& c) {
    return digit_product_plus_one(c.N, 2);
}

static u64 make_feature_prod_3(const CaseData& c) {
    return digit_product_plus_one(c.N, 3);
}

static u64 make_feature_prod_5(const CaseData& c) {
    return digit_product_plus_one(c.N, 5);
}

static u64 make_feature_C2(const CaseData& c) {
    return recurrence_C(c.N, 2);
}

static u64 make_feature_B2(const CaseData& c) {
    return recurrence_B(c.N, 2);
}

static u64 make_feature_G2(const CaseData& c) {
    return recurrence_G(c.N, 2);
}

static u64 make_feature_C3(const CaseData& c) {
    return recurrence_C(c.N, 3);
}

static u64 make_feature_B3(const CaseData& c) {
    return recurrence_B(c.N, 3);
}

static u64 make_feature_G3(const CaseData& c) {
    return recurrence_G(c.N, 3);
}

static u64 make_feature_C5(const CaseData& c) {
    return recurrence_C(c.N, 5);
}

static u64 make_feature_B5(const CaseData& c) {
    return recurrence_B(c.N, 5);
}

static u64 make_feature_G5(const CaseData& c) {
    return recurrence_G(c.N, 5);
}

static u64 make_feature_C7(const CaseData& c) {
    return recurrence_C(c.N, 7);
}

static u64 make_feature_B7(const CaseData& c) {
    return recurrence_B(c.N, 7);
}

static u64 make_feature_G7(const CaseData& c) {
    return recurrence_G(c.N, 7);
}

static u64 make_feature_C11(const CaseData& c) {
    return recurrence_C(c.N, 11);
}

static u64 make_feature_B11(const CaseData& c) {
    return recurrence_B(c.N, 11);
}

static u64 make_feature_G11(const CaseData& c) {
    return recurrence_G(c.N, 11);
}

static u64 make_feature_C13(const CaseData& c) {
    return recurrence_C(c.N, 13);
}

static u64 make_feature_B13(const CaseData& c) {
    return recurrence_B(c.N, 13);
}

static u64 make_feature_G13(const CaseData& c) {
    return recurrence_G(c.N, 13);
}

static mpq_class evaluate_polynomial(
    const std::vector<mpq_class>& coeffs,
    u64 x
) {
    mpq_class xx = mpq_from_u64(x);
    mpq_class result = 0;

    for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) {
        result *= xx;
        result += *it;
    }

    return result;
}

static FitResult fit_exact_polynomial(
    const std::vector<u64>& x,
    const std::vector<u64>& y,
    int degree
) {
    FitResult result;

    const int rows = static_cast<int>(x.size());
    const int cols = degree + 1;

    if (rows == 0 || static_cast<int>(y.size()) != rows) {
        return result;
    }

    std::vector<std::vector<mpq_class>> a(
        rows,
        std::vector<mpq_class>(cols + 1, 0)
    );

    for (int i = 0; i < rows; ++i) {
        mpq_class power = 1;
        mpq_class xx = mpq_from_u64(x[i]);

        for (int j = 0; j < cols; ++j) {
            a[i][j] = power;
            power *= xx;
        }

        a[i][cols] = mpq_from_u64(y[i]);
    }

    int pivot_row = 0;
    std::vector<int> pivot_col;

    for (int col = 0; col < cols && pivot_row < rows; ++col) {
        int found = -1;

        for (int row = pivot_row; row < rows; ++row) {
            if (a[row][col] != 0) {
                found = row;
                break;
            }
        }

        if (found == -1) {
            continue;
        }

        std::swap(a[pivot_row], a[found]);

        const mpq_class pivot = a[pivot_row][col];

        for (int j = col; j <= cols; ++j) {
            a[pivot_row][j] /= pivot;
        }

        for (int row = 0; row < rows; ++row) {
            if (row == pivot_row) {
                continue;
            }

            if (a[row][col] == 0) {
                continue;
            }

            const mpq_class factor = a[row][col];

            for (int j = col; j <= cols; ++j) {
                a[row][j] -= factor * a[pivot_row][j];
            }
        }

        pivot_col.push_back(col);
        ++pivot_row;
    }

    result.rank = pivot_row;

    /*
     * Check consistency:
     *   0 ... 0 | nonzero
     */
    for (int row = 0; row < rows; ++row) {
        bool all_zero = true;

        for (int col = 0; col < cols; ++col) {
            if (a[row][col] != 0) {
                all_zero = false;
                break;
            }
        }

        if (all_zero && a[row][cols] != 0) {
            result.consistent = false;
            return result;
        }
    }

    result.consistent = true;

    if (result.rank != cols) {
        result.possible = false;
        return result;
    }

    result.coeffs.assign(cols, 0);

    for (int row = 0; row < result.rank; ++row) {
        const int col = pivot_col[row];
        result.coeffs[col] = a[row][cols];
    }

    result.possible = true;
    return result;
}

static int count_collisions_with_different_target(
    const std::vector<u64>& x,
    const std::vector<u64>& y
) {
    std::map<u64, u64> first_target;
    int collisions = 0;

    for (std::size_t i = 0; i < x.size(); ++i) {
        auto it = first_target.find(x[i]);

        if (it == first_target.end()) {
            first_target[x[i]] = y[i];
        } else if (it->second != y[i]) {
            ++collisions;
        }
    }

    return collisions;
}

static int count_exact_failures(
    const std::vector<u64>& x,
    const std::vector<u64>& y,
    const std::vector<mpq_class>& coeffs
) {
    int failures = 0;

    for (std::size_t i = 0; i < x.size(); ++i) {
        const mpq_class predicted = evaluate_polynomial(coeffs, x[i]);

        if (predicted != mpq_from_u64(y[i])) {
            ++failures;
        }
    }

    return failures;
}

static void print_polynomial(
    const std::vector<mpq_class>& coeffs
) {
    std::cout << "POLYNOMIAL=";

    bool first = true;

    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        if (coeffs[i] == 0) {
            continue;
        }

        if (!first) {
            std::cout << " + ";
        }

        first = false;

        std::cout << "(" << coeffs[i].get_str() << ")";

        if (i >= 1) {
            std::cout << "*X";
        }

        if (i >= 2) {
            std::cout << "^" << i;
        }
    }

    if (first) {
        std::cout << "0";
    }

    std::cout << "\n";
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<CaseData> result;
    result.reserve(count);

    std::uniform_int_distribution<std::size_t> dist(
        0,
        primes.size() - 1
    );

    while (static_cast<int>(result.size()) < count) {
        const std::size_t i = dist(rng);
        const std::size_t j = dist(rng);

        if (i == j) {
            continue;
        }

        const u64 p = static_cast<u64>(
            std::min(primes[i], primes[j])
        );
        const u64 q = static_cast<u64>(
            std::max(primes[i], primes[j])
        );

        const u64 N = p * q;

        if (!used_N.insert(N).second) {
            continue;
        }

        result.push_back({
            p,
            q,
            N,
            p + q
        });
    }

    return result;
}

int main() {
    std::cout << "START EXPERIMENT 349\n";

    constexpr int TRAIN_COUNT = 600;
    constexpr int VALID_COUNT = 300;
    constexpr int MAX_DEGREE = 6;

    constexpr int PRIME_LIMIT = 100000;

    std::mt19937_64 rng(
        0x349349349ULL
    );

    const std::vector<int> primes = sieve_primes(PRIME_LIMIT);

    /*
     * Avoid tiny primes so sqrt(N), base expansions, and digit structures
     * have a reasonably broad range.
     */
    std::vector<int> usable_primes;

    for (int p : primes) {
        if (p >= 1009) {
            usable_primes.push_back(p);
        }
    }

    std::set<u64> used_N;

    const std::vector<CaseData> train_cases =
        generate_cases(
            usable_primes,
            TRAIN_COUNT,
            rng,
            used_N
        );

    const std::vector<CaseData> valid_cases =
        generate_cases(
            usable_primes,
            VALID_COUNT,
            rng,
            used_N
        );

    std::vector<Candidate> candidates;

    add_candidate(
        candidates,
        "sqrtN",
        train_cases,
        valid_cases,
        feature_sqrtN
    );

    add_candidate(
        candidates,
        "D",
        train_cases,
        valid_cases,
        feature_D
    );

    add_candidate(
        candidates,
        "digit_sum_b2",
        train_cases,
        valid_cases,
        make_feature_digit_sum_2
    );

    add_candidate(
        candidates,
        "digit_sum_b3",
        train_cases,
        valid_cases,
        make_feature_digit_sum_3
    );

    add_candidate(
        candidates,
        "digit_sum_b5",
        train_cases,
        valid_cases,
        make_feature_digit_sum_5
    );

    add_candidate(
        candidates,
        "digit_sum_b7",
        train_cases,
        valid_cases,
        make_feature_digit_sum_7
    );

    add_candidate(
        candidates,
        "digit_sum_b11",
        train_cases,
        valid_cases,
        make_feature_digit_sum_11
    );

    add_candidate(
        candidates,
        "digit_sum_b13",
        train_cases,
        valid_cases,
        make_feature_digit_sum_13
    );

    add_candidate(
        candidates,
        "zero_count_b2",
        train_cases,
        valid_cases,
        make_feature_zero_2
    );

    add_candidate(
        candidates,
        "zero_count_b3",
        train_cases,
        valid_cases,
        make_feature_zero_3
    );

    add_candidate(
        candidates,
        "zero_count_b5",
        train_cases,
        valid_cases,
        make_feature_zero_5
    );

    add_candidate(
        candidates,
        "zero_count_b7",
        train_cases,
        valid_cases,
        make_feature_zero_7
    );

    add_candidate(
        candidates,
        "digit_square_sum_b2",
        train_cases,
        valid_cases,
        make_feature_square_2
    );

    add_candidate(
        candidates,
        "digit_square_sum_b3",
        train_cases,
        valid_cases,
        make_feature_square_3
    );

    add_candidate(
        candidates,
        "digit_square_sum_b5",
        train_cases,
        valid_cases,
        make_feature_square_5
    );

    add_candidate(
        candidates,
        "digit_product_plus_one_b2",
        train_cases,
        valid_cases,
        make_feature_prod_2
    );

    add_candidate(
        candidates,
        "digit_product_plus_one_b3",
        train_cases,
        valid_cases,
        make_feature_prod_3
    );

    add_candidate(
        candidates,
        "digit_product_plus_one_b5",
        train_cases,
        valid_cases,
        make_feature_prod_5
    );

    add_candidate(
        candidates,
        "C2",
        train_cases,
        valid_cases,
        make_feature_C2
    );

    add_candidate(
        candidates,
        "B2",
        train_cases,
        valid_cases,
        make_feature_B2
    );

    add_candidate(
        candidates,
        "G2",
        train_cases,
        valid_cases,
        make_feature_G2
    );

    add_candidate(
        candidates,
        "C3",
        train_cases,
        valid_cases,
        make_feature_C3
    );

    add_candidate(
        candidates,
        "B3",
        train_cases,
        valid_cases,
        make_feature_B3
    );

    add_candidate(
        candidates,
        "G3",
        train_cases,
        valid_cases,
        make_feature_G3
    );

    add_candidate(
        candidates,
        "C5",
        train_cases,
        valid_cases,
        make_feature_C5
    );

    add_candidate(
        candidates,
        "B5",
        train_cases,
        valid_cases,
        make_feature_B5
    );

    add_candidate(
        candidates,
        "G5",
        train_cases,
        valid_cases,
        make_feature_G5
    );

    add_candidate(
        candidates,
        "C7",
        train_cases,
        valid_cases,
        make_feature_C7
    );

    add_candidate(
        candidates,
        "B7",
        train_cases,
        valid_cases,
        make_feature_B7
    );

    add_candidate(
        candidates,
        "G7",
        train_cases,
        valid_cases,
        make_feature_G7
    );

    add_candidate(
        candidates,
        "C11",
        train_cases,
        valid_cases,
        make_feature_C11
    );

    add_candidate(
        candidates,
        "B11",
        train_cases,
        valid_cases,
        make_feature_B11
    );

    add_candidate(
        candidates,
        "G11",
        train_cases,
        valid_cases,
        make_feature_G11
    );

    add_candidate(
        candidates,
        "C13",
        train_cases,
        valid_cases,
        make_feature_C13
    );

    add_candidate(
        candidates,
        "B13",
        train_cases,
        valid_cases,
        make_feature_B13
    );

    add_candidate(
        candidates,
        "G13",
        train_cases,
        valid_cases,
        make_feature_G13
    );

    std::vector<u64> train_y;
    std::vector<u64> valid_y;

    train_y.reserve(train_cases.size());
    valid_y.reserve(valid_cases.size());

    for (const CaseData& c : train_cases) {
        train_y.push_back(c.S);
    }

    for (const CaseData& c : valid_cases) {
        valid_y.push_back(c.S);
    }

    int validated_candidates = 0;

    std::cout << "TRAIN_CASES=" << train_cases.size() << "\n";
    std::cout << "VALID_CASES=" << valid_cases.size() << "\n";
    std::cout << "CANDIDATES=" << candidates.size() << "\n";
    std::cout << "MAX_DEGREE=" << MAX_DEGREE << "\n";

    for (const Candidate& candidate : candidates) {
        const int train_collisions =
            count_collisions_with_different_target(
                candidate.train_x,
                train_y
            );

        const int valid_collisions =
            count_collisions_with_different_target(
                candidate.valid_x,
                valid_y
            );

        std::cout
            << "\n"
            << "CANDIDATE=" << candidate.name
            << "\n"
            << "TRAIN_COLLISIONS=" << train_collisions
            << "\n"
            << "VALID_COLLISIONS=" << valid_collisions
            << "\n";

        for (int degree = 1; degree <= MAX_DEGREE; ++degree) {
            std::cout
                << "DEGREE=" << degree
                << " ";

            if (train_collisions != 0) {
                std::cout
                    << "FIT_STATUS=IMPOSSIBLE_COLLISION"
                    << "\n";
                continue;
            }

            FitResult fit = fit_exact_polynomial(
                candidate.train_x,
                train_y,
                degree
            );

            if (!fit.consistent) {
                std::cout
                    << "FIT_STATUS=INCONSISTENT"
                    << " RANK=" << fit.rank
                    << "\n";
                continue;
            }

            if (!fit.possible) {
                std::cout
                    << "FIT_STATUS=UNDERDETERMINED"
                    << " RANK=" << fit.rank
                    << " REQUIRED_RANK=" << (degree + 1)
                    << "\n";
                continue;
            }

            const int train_failures =
                count_exact_failures(
                    candidate.train_x,
                    train_y,
                    fit.coeffs
                );

            const int valid_failures =
                count_exact_failures(
                    candidate.valid_x,
                    valid_y,
                    fit.coeffs
                );

            std::cout
                << "FIT_STATUS=FULL_RANK"
                << " TRAIN_FAILURES=" << train_failures
                << " VALID_FAILURES=" << valid_failures
                << "\n";

            if (train_failures == 0 &&
                valid_failures == 0) {

                ++validated_candidates;

                std::cout
                    << "VALIDATED_POLYNOMIAL=YES\n";

                print_polynomial(fit.coeffs);
            } else {
                std::cout
                    << "VALIDATED_POLYNOMIAL=NO\n";
            }
        }
    }

    /*
     * A final deterministic sanity check on the generated semiprimes.
     */
    bool generation_ok = true;

    for (const CaseData& c : train_cases) {
        if (c.p * c.q != c.N ||
            c.p + c.q != c.S ||
            c.p == c.q) {
            generation_ok = false;
            break;
        }
    }

    for (const CaseData& c : valid_cases) {
        if (c.p * c.q != c.N ||
            c.p + c.q != c.S ||
            c.p == c.q) {
            generation_ok = false;
            break;
        }
    }

    std::cout
        << "\n"
        << "GENERATION_SANITY="
        << (generation_ok ? "PASS" : "FAIL")
        << "\n";

    std::cout
        << "VALIDATED_CANDIDATES="
        << validated_candidates
        << "\n";

    if (validated_candidates == 0) {
        std::cout
            << "DISCOVERY_STATUS="
            << "NO_VALIDATED_UNIVARIATE_POLYNOMIAL\n";
    } else {
        std::cout
            << "DISCOVERY_STATUS="
            << "VALIDATED_UNIVARIATE_POLYNOMIAL_FOUND\n";
    }

    std::cout << "FINISHED EXPERIMENT 349\n";

    return 0;
}
