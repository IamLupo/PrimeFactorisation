#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_set>
#include <vector>

using u64 = std::uint64_t;

struct Pair {
    int m1;
    int m2;
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

struct CutStats {
    u64 success = 0;
    u64 matching = 0;
    u64 total_unique = 0;
    u64 max_unique = 0;
    u64 total_prime_tests = 0;
    u64 max_prime_tests = 0;
};

std::vector<bool> build_sieve(int limit) {
    std::vector<bool> prime(limit + 1, true);

    if (limit >= 0) {
        prime[0] = false;
    }

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
    Generate the determinant-one pairs

        m2*k1 - m1*k2 = 1

    with

        1 <= m1,m2 <= M_LIMIT
        1 <= k1 <= k2 <= K_LIMIT
        gcd(k1,k2)=1

    The earlier experiments found 3647 such pairs
    for M_LIMIT=7 and K_LIMIT=1000.
*/
std::vector<Pair> build_determinant_one_pairs(
    int k_limit,
    int m_limit) {

    std::vector<Pair> pairs;

    for (int m1 = 1; m1 <= m_limit; ++m1) {
        for (int m2 = 1; m2 <= m_limit; ++m2) {

            for (int k1 = 1; k1 <= k_limit; ++k1) {

                long long numerator =
                    1LL * m2 * k1 - 1;

                if (numerator <= 0) {
                    continue;
                }

                if (numerator % m1 != 0) {
                    continue;
                }

                int k2 =
                    static_cast<int>(numerator / m1);

                if (k2 < k1 || k2 > k_limit) {
                    continue;
                }

                if (std::gcd(k1, k2) != 1) {
                    continue;
                }

                u64 a =
                    static_cast<u64>(k1) +
                    static_cast<u64>(k2);

                u64 M =
                    static_cast<u64>(k1) *
                    static_cast<u64>(k2);

                pairs.push_back({
                    m1,
                    m2,
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

            if (lhs.m1 != rhs.m1) {
                return lhs.m1 < rhs.m1;
            }

            if (lhs.m2 != rhs.m2) {
                return lhs.m2 < rhs.m2;
            }

            if (lhs.k1 != rhs.k1) {
                return lhs.k1 < rhs.k1;
            }

            return lhs.k2 < rhs.k2;
        }
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
        u64 p = static_cast<u64>(primes[dist(rng)]);
        u64 q = static_cast<u64>(primes[dist(rng)]);

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        u64 N = p * q;

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

        return {N, p, q, s};
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

std::vector<u64> generate_candidates(
    const std::vector<Pair>& pairs,
    u64 s,
    unsigned int j_limit) {

    std::vector<u64> values;

    for (const Pair& pair : pairs) {

        if (pair.a > s) {
            continue;
        }

        for (unsigned int j = 0; j <= j_limit; ++j) {

            u64 value =
                pair.a +
                static_cast<u64>(j) * pair.M;

            if (value > s) {
                break;
            }

            if (value < 2) {
                continue;
            }

            values.push_back(value);
        }
    }

    std::sort(values.begin(), values.end());

    values.erase(
        std::unique(values.begin(), values.end()),
        values.end()
    );

    return values;
}

u64 test_candidates(
    u64 N,
    u64 p,
    const std::vector<u64>& candidates,
    const std::vector<bool>& prime,
    u64& found_factor) {

    u64 tests = 0;

    /*
        Test numerically descending, matching the natural
        factor-search direction.
    */
    for (auto it = candidates.rbegin();
         it != candidates.rend();
         ++it) {

        u64 candidate = *it;

        if (candidate >= prime.size()) {
            continue;
        }

        if (!prime[
                static_cast<std::size_t>(candidate)
            ]) {
            continue;
        }

        ++tests;

        if (N % candidate == 0) {
            found_factor = candidate;
            return tests;
        }
    }

    found_factor = 0;
    return tests;
}

u64 find_factor_in_candidate_set(
    u64 N,
    const std::vector<u64>& candidates) {

    for (u64 candidate : candidates) {
        if (candidate != 0 && N % candidate == 0) {
            return candidate;
        }
    }

    return 0;
}

double average_u64(u64 value, int count) {
    if (count == 0) {
        return 0.0;
    }

    return static_cast<double>(value) /
           static_cast<double>(count);
}

int main() {
    constexpr int CASE_COUNT = 500;
    constexpr int PRIME_LOW = 10000;
    constexpr int PRIME_HIGH = 100000;
    constexpr int PRIME_LIMIT = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_LIMIT = 7;

    const std::vector<unsigned int> J_CUTS = {
        1,
        2,
        5,
        10,
        20,
        50,
        100,
        250,
        500,
        1000,
        2000,
        5000
    };

    std::cout << "START EXPERIMENT 480\n";

    const auto sieve = build_sieve(PRIME_LIMIT);

    const auto primes =
        generate_primes(
            PRIME_LOW,
            PRIME_HIGH,
            sieve
        );

    const auto pairs =
        build_determinant_one_pairs(
            K_LIMIT,
            M_LIMIT
        );

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << "\n";

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << "\n";

    /*
        Statistics for each J cutoff.
    */
    std::vector<CutStats> stats(J_CUTS.size());

    std::mt19937_64 rng(4801234567ULL);

    u64 ordinary_total_tests = 0;
    u64 ordinary_total_distance = 0;

    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;
    u64 first_s = 0;
    u64 first_ordinary_tests = 0;

    std::vector<u64> first_candidate_counts(
        J_CUTS.size(),
        0
    );

    std::vector<u64> first_prime_tests(
        J_CUTS.size(),
        0
    );

    std::vector<u64> first_found(
        J_CUTS.size(),
        0
    );

    const auto start =
        std::chrono::steady_clock::now();

    for (int case_id = 0;
         case_id < CASE_COUNT;
         ++case_id) {

        CaseData data =
            generate_semiprime_case(
                primes,
                rng
            );

        u64 ordinary_factor = 0;

        u64 ordinary_tests =
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
            first_ordinary_tests =
                ordinary_tests;
        }

        for (std::size_t cut_index = 0;
             cut_index < J_CUTS.size();
             ++cut_index) {

            unsigned int j_limit =
                J_CUTS[cut_index];

            std::vector<u64> candidates =
                generate_candidates(
                    pairs,
                    data.s,
                    j_limit
                );

            u64 candidate_count =
                static_cast<u64>(
                    candidates.size()
                );

            stats[cut_index].total_unique +=
                candidate_count;

            stats[cut_index].max_unique =
                std::max(
                    stats[cut_index].max_unique,
                    candidate_count
                );

            u64 found_factor = 0;

            u64 prime_tests =
                test_candidates(
                    data.N,
                    data.p,
                    candidates,
                    sieve,
                    found_factor
                );

            stats[cut_index].total_prime_tests +=
                prime_tests;

            stats[cut_index].max_prime_tests =
                std::max(
                    stats[cut_index].max_prime_tests,
                    prime_tests
                );

            if (found_factor != 0) {
                ++stats[cut_index].success;

                if (found_factor == data.p) {
                    ++stats[cut_index].matching;
                }
            }

            if (case_id == 0) {
                first_candidate_counts[cut_index] =
                    candidate_count;

                first_prime_tests[cut_index] =
                    prime_tests;

                first_found[cut_index] =
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
         i < J_CUTS.size();
         ++i) {

        std::cout
            << "FIRST_J="
            << J_CUTS[i]
            << " CANDIDATES="
            << first_candidate_counts[i]
            << " PRIME_TESTS="
            << first_prime_tests[i]
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
        << average_u64(
               ordinary_total_tests,
               CASE_COUNT
           )
        << "\n";

    std::cout
        << "ORDINARY_AVG_DISTANCE="
        << average_u64(
               ordinary_total_distance,
               CASE_COUNT
           )
        << "\n";

    for (std::size_t i = 0;
         i < J_CUTS.size();
         ++i) {

        std::cout
            << "J="
            << J_CUTS[i]
            << " SUCCESS="
            << stats[i].success
            << " MATCHING="
            << stats[i].matching
            << " AVG_UNIQUE_CANDIDATES="
            << average_u64(
                   stats[i].total_unique,
                   CASE_COUNT
               )
            << " MAX_UNIQUE_CANDIDATES="
            << stats[i].max_unique
            << " AVG_PRIME_TESTS="
            << average_u64(
                   stats[i].total_prime_tests,
                   CASE_COUNT
               )
            << " MAX_PRIME_TESTS="
            << stats[i].max_prime_tests
            << "\n";
    }

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 480\n";

    return 0;
}
