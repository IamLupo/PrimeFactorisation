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

struct StageStats {
    u64 cases_entered = 0;
    u64 cases_found = 0;
    u64 factor_tests = 0;
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
                    static_cast<int>(
                        numerator / m1
                    );

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
    u64& factor) {

    u64 tests = 0;

    for (u64 candidate = s; candidate >= 2; --candidate) {
        ++tests;

        if (N % candidate == 0) {
            factor = candidate;
            return tests;
        }
    }

    factor = 0;
    return tests;
}

/*
    Generate all prime candidates belonging to the newly
    added determinant-one modulus range:

        previous_M < k1*k2 <= new_M

    Candidate values satisfy

        r = (k1+k2) + j*(k1*k2)

    and are restricted to r <= s.

    'seen' prevents duplicate prime candidates from being
    tested twice across stages.
*/
std::vector<u64> generate_new_candidates(
    const std::vector<Pair>& pairs,
    u64 s,
    u64 previous_M,
    u64 new_M,
    const std::vector<bool>& prime,
    std::vector<bool>& seen) {

    std::vector<u64> candidates;

    for (const Pair& pair : pairs) {

        if (pair.M <= previous_M) {
            continue;
        }

        if (pair.M > new_M) {
            break;
        }

        if (pair.a > s) {
            continue;
        }

        u64 value = pair.a;

        while (value <= s) {

            if (value >= 2 &&
                value < prime.size()) {

                const std::size_t index =
                    static_cast<std::size_t>(
                        value
                    );

                if (prime[index] &&
                    !seen[index]) {

                    seen[index] = true;
                    candidates.push_back(value);
                }
            }

            if (value > s - pair.M) {
                break;
            }

            value += pair.M;
        }
    }

    /*
        The search order is descending numerical candidate.
    */
    std::sort(
        candidates.begin(),
        candidates.end(),
        std::greater<u64>()
    );

    return candidates;
}

double average(
    u64 total,
    u64 count) {

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

    const std::vector<u64> STAGES = {
        10,
        20,
        50,
        250,
        1000
    };

    std::cout
        << "START EXPERIMENT 486\n";

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
        4861234567ULL
    );

    std::vector<StageStats> stage_stats(
        STAGES.size()
    );

    u64 ordinary_total_tests = 0;
    u64 ordinary_total_distance = 0;

    u64 adaptive_total_tests = 0;
    u64 adaptive_max_tests = 0;

    u64 adaptive_success = 0;
    u64 adaptive_matching = 0;
    u64 fallback_count = 0;

    /*
        First-case diagnostics.
    */
    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;
    u64 first_s = 0;
    u64 first_ordinary_tests = 0;

    std::vector<u64> first_entered(
        STAGES.size(),
        0
    );

    std::vector<u64> first_new_candidates(
        STAGES.size(),
        0
    );

    std::vector<u64> first_stage_tests(
        STAGES.size(),
        0
    );

    std::vector<u64> first_stage_found(
        STAGES.size(),
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

        ordinary_total_tests +=
            ordinary_tests;

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

        /*
            Every stage uses the same seen set.
            Therefore the stage sets are genuinely nested.
        */
        std::vector<bool> seen(
            static_cast<std::size_t>(
                PRIME_LIMIT
            ) + 1,
            false
        );

        u64 case_total_tests = 0;
        u64 found_factor = 0;

        u64 previous_M = 0;
        int found_stage = -1;

        for (std::size_t stage = 0;
             stage < STAGES.size();
             ++stage) {

            const u64 current_M =
                STAGES[stage];

            StageStats& stats =
                stage_stats[stage];

            ++stats.cases_entered;

            const auto candidates =
                generate_new_candidates(
                    pairs,
                    data.s,
                    previous_M,
                    current_M,
                    sieve,
                    seen
                );

            if (case_id == 0) {
                first_entered[stage] =
                    stats.cases_entered;

                first_new_candidates[stage] =
                    candidates.size();
            }

            u64 stage_tests = 0;
            u64 stage_factor = 0;

            for (u64 candidate : candidates) {

                ++stage_tests;
                ++case_total_tests;

                if (data.N % candidate == 0) {
                    stage_factor = candidate;
                    found_factor = candidate;
                    break;
                }
            }

            stats.factor_tests += stage_tests;

            stats.max_factor_tests =
                std::max(
                    stats.max_factor_tests,
                    stage_tests
                );

            if (case_id == 0) {
                first_stage_tests[stage] =
                    stage_tests;

                /*
                    This value is recorded at the exact
                    stage where the factor is found.
                */
                first_stage_found[stage] =
                    stage_factor;
            }

            if (stage_factor != 0) {

                ++stats.cases_found;

                found_stage =
                    static_cast<int>(stage);

                if (stage_factor == data.p) {
                    ++adaptive_matching;
                }

                break;
            }

            previous_M = current_M;
        }

        adaptive_total_tests +=
            case_total_tests;

        adaptive_max_tests =
            std::max(
                adaptive_max_tests,
                case_total_tests
            );

        if (found_factor != 0) {
            ++adaptive_success;
        } else {
            ++fallback_count;
        }

        if (case_id == 0) {
            (void)found_stage;
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

    for (std::size_t stage = 0;
         stage < STAGES.size();
         ++stage) {

        std::cout
            << "FIRST_M="
            << STAGES[stage]
            << " NEW_CANDIDATES="
            << first_new_candidates[stage]
            << " STAGE_FACTOR_TESTS="
            << first_stage_tests[stage]
            << " FOUND="
            << first_stage_found[stage]
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

    std::cout
        << "ADAPTIVE_SUCCESS="
        << adaptive_success
        << "\n";

    std::cout
        << "ADAPTIVE_MATCHING="
        << adaptive_matching
        << "\n";

    std::cout
        << "ADAPTIVE_FALLBACK_COUNT="
        << fallback_count
        << "\n";

    std::cout
        << "ADAPTIVE_TOTAL_FACTOR_TESTS="
        << adaptive_total_tests
        << "\n";

    std::cout
        << "ADAPTIVE_AVG_FACTOR_TESTS="
        << average(
               adaptive_total_tests,
               CASE_COUNT
           )
        << "\n";

    std::cout
        << "ADAPTIVE_MAX_FACTOR_TESTS="
        << adaptive_max_tests
        << "\n";

    if (adaptive_success > 0) {
        std::cout
            << "ADAPTIVE_SUCCESS_PERCENT="
            << 100.0 *
               static_cast<double>(
                   adaptive_success
               ) /
               static_cast<double>(
                   CASE_COUNT
               )
            << "\n";
    } else {
        std::cout
            << "ADAPTIVE_SUCCESS_PERCENT=0\n";
    }

    for (std::size_t stage = 0;
         stage < STAGES.size();
         ++stage) {

        const StageStats& stats =
            stage_stats[stage];

        std::cout
            << "STAGE_M="
            << STAGES[stage]
            << " CASES_ENTERED="
            << stats.cases_entered
            << " CASES_FOUND="
            << stats.cases_found
            << " AVG_STAGE_TESTS="
            << average(
                   stats.factor_tests,
                   stats.cases_entered
               )
            << " MAX_STAGE_TESTS="
            << stats.max_factor_tests
            << "\n";
    }

    if (ordinary_total_tests > 0) {

        const double ratio =
            static_cast<double>(
                adaptive_total_tests
            ) /
            static_cast<double>(
                ordinary_total_tests
            );

        std::cout
            << "ADAPTIVE_OVER_ORDINARY="
            << ratio
            << "\n";

        std::cout
            << "ADAPTIVE_TEST_REDUCTION_PERCENT="
            << 100.0 *
               (1.0 - ratio)
            << "\n";
    } else {
        std::cout
            << "ADAPTIVE_OVER_ORDINARY=0\n";

        std::cout
            << "ADAPTIVE_TEST_REDUCTION_PERCENT=100\n";
    }

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 486\n";

    return 0;
}
