#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct Pair {
    int k1;
    int k2;
    u64 a;
    u64 M;
};

struct CaseData {
    u64 N;
    u64 p;
    u64 q;
    u64 s;
};

struct Stats {
    u64 success = 0;
    u64 matching = 0;

    u64 total_integer_candidates = 0;
    u64 max_integer_candidates = 0;

    u64 total_unique_candidates = 0;
    u64 max_unique_candidates = 0;

    u64 total_prime_candidates = 0;
    u64 max_prime_candidates = 0;

    u64 total_factor_tests = 0;
    u64 max_factor_tests = 0;
};

std::vector<bool> build_sieve(int limit) {
    std::vector<bool> prime(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    prime[0] = false;

    if (limit >= 1) {
        prime[1] = false;
    }

    for (int i = 2; 1LL * i * i <= limit; ++i) {
        if (!prime[i]) {
            continue;
        }

        for (int j = i * i; j <= limit; j += i) {
            prime[j] = false;
        }
    }

    return prime;
}

std::vector<int> generate_primes(
    int low,
    int high,
    const std::vector<bool>& prime) {

    std::vector<int> result;

    for (int x = low; x <= high; ++x) {
        if (prime[x]) {
            result.push_back(x);
        }
    }

    return result;
}

/*
    Build the determinant-one family

        m2*k1 - m1*k2 = 1

    with

        1 <= m1,m2 <= 7
        2 <= k1 <= k2 <= 1000
        gcd(k1,k2)=1.

    We exclude k=1 so M=1 cannot trivialize the search.
*/
std::vector<Pair> build_determinant_one_pairs(
    int k_limit,
    int m_limit) {

    std::vector<Pair> pairs;

    for (int m1 = 1; m1 <= m_limit; ++m1) {
        for (int m2 = 1; m2 <= m_limit; ++m2) {

            for (int k1 = 2; k1 <= k_limit; ++k1) {

                const long long numerator =
                    1LL * m2 * k1 - 1;

                if (numerator <= 0) {
                    continue;
                }

                if (numerator % m1 != 0) {
                    continue;
                }

                const int k2 =
                    static_cast<int>(numerator / m1);

                if (k2 < k1 || k2 > k_limit) {
                    continue;
                }

                if (k2 < 2) {
                    continue;
                }

                if (std::gcd(k1, k2) != 1) {
                    continue;
                }

                const u64 a =
                    static_cast<u64>(k1) +
                    static_cast<u64>(k2);

                const u64 M =
                    static_cast<u64>(k1) *
                    static_cast<u64>(k2);

                pairs.push_back({
                    k1,
                    k2,
                    a,
                    M
                });
            }
        }
    }

    std::sort(
        pairs.begin(),
        pairs.end(),
        [](const Pair& lhs, const Pair& rhs) {
            if (lhs.M != rhs.M) {
                return lhs.M < rhs.M;
            }

            if (lhs.a != rhs.a) {
                return lhs.a < rhs.a;
            }

            if (lhs.k1 != rhs.k1) {
                return lhs.k1 < rhs.k1;
            }

            return lhs.k2 < rhs.k2;
        }
    );

    pairs.erase(
        std::unique(
            pairs.begin(),
            pairs.end(),
            [](const Pair& lhs, const Pair& rhs) {
                return
                    lhs.k1 == rhs.k1 &&
                    lhs.k2 == rhs.k2;
            }
        ),
        pairs.end()
    );

    return pairs;
}

CaseData generate_semiprime_case(
    const std::vector<int>& primes,
    std::mt19937_64& rng) {

    std::uniform_int_distribution<std::size_t> dist(
        0,
        primes.size() - 1
    );

    while (true) {
        u64 p =
            static_cast<u64>(
                primes[dist(rng)]
            );

        u64 q =
            static_cast<u64>(
                primes[dist(rng)]
            );

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        const u64 N = p * q;

        u64 s =
            static_cast<u64>(
                std::sqrt(
                    static_cast<long double>(N)
                )
            );

        while ((s + 1) <= N / (s + 1)) {
            ++s;
        }

        while (s > N / s) {
            --s;
        }

        return {
            N,
            p,
            q,
            s
        };
    }
}

u64 ordinary_factor_tests(
    u64 N,
    u64 s,
    u64& found_factor) {

    u64 tests = 0;

    for (u64 candidate = s; candidate >= 2; --candidate) {
        ++tests;

        if (N % candidate == 0) {
            found_factor = candidate;
            return tests;
        }
    }

    found_factor = 0;
    return tests;
}

/*
    Generate the blind candidate set from the congruence

        r == a (mod M)

    for every determinant-one pair with

        M <= M_LIMIT.

    Only values <= s are generated.

    We retain all integers first, then sort/unique them.
*/
std::vector<u64> generate_congruence_candidates(
    const std::vector<Pair>& pairs,
    u64 s,
    u64 M_LIMIT,
    u64& integer_candidate_count) {

    std::vector<u64> values;

    integer_candidate_count = 0;

    for (const Pair& pair : pairs) {

        if (pair.M > M_LIMIT) {
            break;
        }

        if (pair.a > s) {
            continue;
        }

        /*
            p = a + jM, j >= 0.
        */
        for (u64 value = pair.a;
             value <= s;
             ) {

            values.push_back(value);
            ++integer_candidate_count;

            /*
                Overflow-safe step.
            */
            if (value > s - pair.M) {
                break;
            }

            value += pair.M;
        }
    }

    std::sort(
        values.begin(),
        values.end()
    );

    values.erase(
        std::unique(
            values.begin(),
            values.end()
        ),
        values.end()
    );

    return values;
}

/*
    Test the candidate values in descending numerical order.

    Returns the number of actual modular divisibility tests.
*/
u64 test_prime_candidates(
    u64 N,
    const std::vector<u64>& candidates,
    const std::vector<bool>& prime,
    u64& found_factor,
    u64& prime_candidate_count) {

    u64 factor_tests = 0;
    prime_candidate_count = 0;

    for (auto it = candidates.rbegin();
         it != candidates.rend();
         ++it) {

        const u64 candidate = *it;

        if (candidate < 2) {
            continue;
        }

        if (candidate >= prime.size()) {
            continue;
        }

        if (!prime[
                static_cast<std::size_t>(candidate)
            ]) {
            continue;
        }

        ++prime_candidate_count;
        ++factor_tests;

        if (N % candidate == 0) {
            found_factor = candidate;
            return factor_tests;
        }
    }

    found_factor = 0;
    return factor_tests;
}

double average(
    u64 total,
    int count) {

    if (count == 0) {
        return 0.0;
    }

    return
        static_cast<double>(total) /
        static_cast<double>(count);
}

int main() {
    constexpr int CASE_COUNT = 500;

    constexpr int PRIME_LOW = 10000;
    constexpr int PRIME_HIGH = 100000;
    constexpr int PRIME_LIMIT = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_FAMILY_LIMIT = 7;

    const std::vector<u64> M_LIMITS = {
        6,
        10,
        20,
        50,
        100,
        250,
        500,
        1000
    };

    std::cout
        << "START EXPERIMENT 484\n";

    const auto sieve =
        build_sieve(PRIME_LIMIT);

    const auto primes =
        generate_primes(
            PRIME_LOW,
            PRIME_HIGH,
            sieve
        );

    const auto pairs =
        build_determinant_one_pairs(
            K_LIMIT,
            M_FAMILY_LIMIT
        );

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << "\n";

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << "\n";

    std::mt19937_64 rng(
        4841234567ULL
    );

    std::vector<Stats> stats(
        M_LIMITS.size()
    );

    u64 ordinary_total_tests = 0;
    u64 ordinary_total_distance = 0;

    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;
    u64 first_s = 0;
    u64 first_ordinary_tests = 0;

    std::vector<u64> first_integer_candidates(
        M_LIMITS.size(),
        0
    );

    std::vector<u64> first_unique_candidates(
        M_LIMITS.size(),
        0
    );

    std::vector<u64> first_prime_candidates(
        M_LIMITS.size(),
        0
    );

    std::vector<u64> first_factor_tests(
        M_LIMITS.size(),
        0
    );

    std::vector<u64> first_found(
        M_LIMITS.size(),
        0
    );

    const auto start =
        std::chrono::steady_clock::now();

    for (int case_id = 0;
         case_id < CASE_COUNT;
         ++case_id) {

        const CaseData data =
            generate_semiprime_case(
                primes,
                rng
            );

        u64 ordinary_factor = 0;

        const u64 ordinary_tests =
            ordinary_factor_tests(
                data.N,
                data.s,
                ordinary_factor
            );

        ordinary_total_tests += ordinary_tests;

        ordinary_total_distance +=
            data.s - ordinary_factor;

        if (case_id == 0) {
            first_N = data.N;
            first_p = data.p;
            first_q = data.q;
            first_s = data.s;
            first_ordinary_tests = ordinary_tests;
        }

        for (std::size_t i = 0;
             i < M_LIMITS.size();
             ++i) {

            const u64 M_LIMIT =
                M_LIMITS[i];

            u64 integer_candidate_count = 0;

            const auto candidates =
                generate_congruence_candidates(
                    pairs,
                    data.s,
                    M_LIMIT,
                    integer_candidate_count
                );

            const u64 unique_candidate_count =
                static_cast<u64>(
                    candidates.size()
                );

            u64 found_factor = 0;
            u64 prime_candidate_count = 0;

            const u64 factor_tests =
                test_prime_candidates(
                    data.N,
                    candidates,
                    sieve,
                    found_factor,
                    prime_candidate_count
                );

            Stats& current = stats[i];

            if (found_factor != 0) {
                ++current.success;

                if (found_factor == data.p) {
                    ++current.matching;
                }
            }

            current.total_integer_candidates +=
                integer_candidate_count;

            current.max_integer_candidates =
                std::max(
                    current.max_integer_candidates,
                    integer_candidate_count
                );

            current.total_unique_candidates +=
                unique_candidate_count;

            current.max_unique_candidates =
                std::max(
                    current.max_unique_candidates,
                    unique_candidate_count
                );

            current.total_prime_candidates +=
                prime_candidate_count;

            current.max_prime_candidates =
                std::max(
                    current.max_prime_candidates,
                    prime_candidate_count
                );

            current.total_factor_tests +=
                factor_tests;

            current.max_factor_tests =
                std::max(
                    current.max_factor_tests,
                    factor_tests
                );

            if (case_id == 0) {
                first_integer_candidates[i] =
                    integer_candidate_count;

                first_unique_candidates[i] =
                    unique_candidate_count;

                first_prime_candidates[i] =
                    prime_candidate_count;

                first_factor_tests[i] =
                    factor_tests;

                first_found[i] =
                    found_factor;
            }
        }

        if ((case_id + 1) % 100 == 0) {
            std::cout
                << "PROGRESS="
                << (case_id + 1)
                << "/"
                << CASE_COUNT
                << "\n";
        }
    }

    const auto end =
        std::chrono::steady_clock::now();

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    std::cout
        << "FIRST_N="
        << first_N
        << "\n";

    std::cout
        << "FIRST_P="
        << first_p
        << "\n";

    std::cout
        << "FIRST_Q="
        << first_q
        << "\n";

    std::cout
        << "FIRST_S="
        << first_s
        << "\n";

    std::cout
        << "FIRST_ORDINARY_GCD_TESTS="
        << first_ordinary_tests
        << "\n";

    for (std::size_t i = 0;
         i < M_LIMITS.size();
         ++i) {

        std::cout
            << "FIRST_M_LIMIT="
            << M_LIMITS[i]
            << " INTEGER_CANDIDATES="
            << first_integer_candidates[i]
            << " UNIQUE_CANDIDATES="
            << first_unique_candidates[i]
            << " PRIME_CANDIDATES="
            << first_prime_candidates[i]
            << " FACTOR_TESTS="
            << first_factor_tests[i]
            << " FOUND="
            << first_found[i]
            << "\n";
    }

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << "\n";

    std::cout
        << "ORDINARY_TOTAL_GCD_TESTS="
        << ordinary_total_tests
        << "\n";

    std::cout
        << "ORDINARY_AVG_GCD_TESTS="
        << average(
               ordinary_total_tests,
               CASE_COUNT
           )
        << "\n";

    std::cout
        << "ORDINARY_AVG_DISTANCE="
        << average(
               ordinary_total_distance,
               CASE_COUNT
           )
        << "\n";

    for (std::size_t i = 0;
         i < M_LIMITS.size();
         ++i) {

        const Stats& current = stats[i];

        std::cout
            << "M_LIMIT="
            << M_LIMITS[i]
            << " SUCCESS="
            << current.success
            << " MATCHING="
            << current.matching
            << " AVG_INTEGER_CANDIDATES="
            << average(
                   current.total_integer_candidates,
                   CASE_COUNT
               )
            << " MAX_INTEGER_CANDIDATES="
            << current.max_integer_candidates
            << " AVG_UNIQUE_CANDIDATES="
            << average(
                   current.total_unique_candidates,
                   CASE_COUNT
               )
            << " MAX_UNIQUE_CANDIDATES="
            << current.max_unique_candidates
            << " AVG_PRIME_CANDIDATES="
            << average(
                   current.total_prime_candidates,
                   CASE_COUNT
               )
            << " MAX_PRIME_CANDIDATES="
            << current.max_prime_candidates
            << " AVG_FACTOR_TESTS="
            << average(
                   current.total_factor_tests,
                   CASE_COUNT
               )
            << " MAX_FACTOR_TESTS="
            << current.max_factor_tests
            << "\n";
    }

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 484\n";

    return 0;
}
