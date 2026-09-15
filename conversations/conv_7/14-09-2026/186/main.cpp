#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

struct Pair {
    int k1;
    int k2;
    u64 a;
    u64 M;
};

struct Candidate {
    u64 value;
    int pair_index;
    u64 j;
};

struct CaseData {
    u64 N;
    u64 p;
    u64 q;
    u64 s;
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

std::vector<int> generate_primes(int low,
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

std::vector<Pair> build_determinant_one_pairs(int k_limit, int m_limit) {
    std::vector<Pair> pairs;

    for (int k1 = 1; k1 <= k_limit; ++k1) {
        for (int k2 = k1; k2 <= k_limit; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            u64 M = static_cast<u64>(k1) * static_cast<u64>(k2);

            if (M > static_cast<u64>(m_limit) *
                        static_cast<u64>(k_limit)) {
                // Conservative pruning only for very large products.
                // Kept intentionally weak so valid pairs are not lost.
            }

            u64 a = static_cast<u64>(k1) + static_cast<u64>(k2);

            pairs.push_back({k1, k2, a, M});
        }
    }

    return pairs;
}

std::vector<Candidate> build_best_candidates(const std::vector<Pair>& pairs,
                                             u64 s) {
    std::vector<Candidate> candidates;
    candidates.reserve(pairs.size());

    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const Pair& pair = pairs[i];

        if (pair.a > s) {
            continue;
        }

        u64 j = (s - pair.a) / pair.M;
        u64 value = pair.a + j * pair.M;

        if (value < 2) {
            continue;
        }

        candidates.push_back({
            value,
            static_cast<int>(i),
            j
        });
    }

    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const Candidate& lhs, const Candidate& rhs) {
            if (lhs.value != rhs.value) {
                return lhs.value > rhs.value;
            }

            if (lhs.pair_index != rhs.pair_index) {
                return lhs.pair_index < rhs.pair_index;
            }

            return lhs.j > rhs.j;
        }
    );

    return candidates;
}

std::vector<Candidate> deduplicate_candidates(
    const std::vector<Candidate>& candidates) {

    std::vector<Candidate> unique_candidates;
    unique_candidates.reserve(candidates.size());

    u64 previous = 0;
    bool first = true;

    for (const Candidate& candidate : candidates) {
        if (first || candidate.value != previous) {
            unique_candidates.push_back(candidate);
            previous = candidate.value;
            first = false;
        }
    }

    return unique_candidates;
}

CaseData generate_semiprime_case(const std::vector<int>& primes,
                                 std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> dist(0, primes.size() - 1);

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

        u64 s = static_cast<u64>(
            std::sqrt(static_cast<long double>(N))
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

u64 ordinary_factor_tests(u64 N, u64 s, u64& found_factor) {
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

u64 find_factor_from_candidates(
    u64 N,
    const std::vector<Candidate>& candidates,
    const std::vector<bool>& prime,
    std::size_t limit,
    u64& found_factor) {

    u64 tests = 0;

    std::size_t count = std::min(limit, candidates.size());

    for (std::size_t i = 0; i < count; ++i) {
        u64 candidate = candidates[i].value;

        if (candidate > static_cast<u64>(prime.size() - 1)) {
            continue;
        }

        if (!prime[static_cast<std::size_t>(candidate)]) {
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

u64 candidate_rank_of_factor(
    const std::vector<Candidate>& candidates,
    u64 factor) {

    for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i].value == factor) {
            return static_cast<u64>(i + 1);
        }
    }

    return 0;
}

double average_u64(u64 total, int count) {
    if (count == 0) {
        return 0.0;
    }

    return static_cast<double>(total) /
           static_cast<double>(count);
}

void main_experiment() {
    constexpr int CASE_COUNT = 500;
    constexpr int PRIME_LOW = 10000;
    constexpr int PRIME_HIGH = 100000;
    constexpr int PRIME_LIMIT = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_LIMIT = 7;

    const std::vector<std::size_t> CUTS = {
        10,
        25,
        50,
        100,
        250,
        500,
        1000,
        2000,
        4000
    };

    std::cout << "START EXPERIMENT 479\n";

    const auto sieve = build_sieve(PRIME_LIMIT);
    const auto primes = generate_primes(PRIME_LOW, PRIME_HIGH, sieve);

    const auto pairs = build_determinant_one_pairs(
        K_LIMIT,
        M_LIMIT
    );

    std::cout << "PAIR_COUNT=" << pairs.size() << "\n";
    std::cout << "PRIME_POOL_SIZE=" << primes.size() << "\n";

    std::mt19937_64 rng(4791234567ULL);

    u64 ordinary_total_tests = 0;
    u64 ordinary_total_distance = 0;

    u64 all_unique_total = 0;
    u64 all_unique_max = 0;

    u64 rank_sum = 0;
    u64 rank_max = 0;
    u64 rank_missing = 0;

    struct CutStats {
        u64 success = 0;
        u64 matching = 0;
        u64 total_prime_tests = 0;
        u64 max_prime_tests = 0;
    };

    std::vector<CutStats> stats(CUTS.size());

    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;
    u64 first_s = 0;
    u64 first_ordinary = 0;
    u64 first_rank = 0;
    u64 first_unique = 0;
    u64 first_top10 = 0;
    u64 first_top100 = 0;
    u64 first_top500 = 0;

    const auto start = std::chrono::steady_clock::now();

    for (int case_id = 0; case_id < CASE_COUNT; ++case_id) {
        CaseData data = generate_semiprime_case(primes, rng);

        if (case_id == 0) {
            first_N = data.N;
            first_p = data.p;
            first_q = data.q;
            first_s = data.s;
        }

        u64 ordinary_factor = 0;
        u64 ordinary_tests = ordinary_factor_tests(
            data.N,
            data.s,
            ordinary_factor
        );

        ordinary_total_tests += ordinary_tests;

        if (data.p >= ordinary_factor) {
            ordinary_total_distance +=
                data.s - ordinary_factor;
        } else {
            ordinary_total_distance +=
                data.s - data.p;
        }

        std::vector<Candidate> candidates =
            build_best_candidates(pairs, data.s);

        std::vector<Candidate> unique_candidates =
            deduplicate_candidates(candidates);

        all_unique_total += unique_candidates.size();

        if (unique_candidates.size() > all_unique_max) {
            all_unique_max = unique_candidates.size();
        }

        u64 rank = candidate_rank_of_factor(
            unique_candidates,
            data.p
        );

        if (rank == 0) {
            ++rank_missing;
        } else {
            rank_sum += rank;
            rank_max = std::max(rank_max, rank);
        }

        if (case_id == 0) {
            first_ordinary = ordinary_tests;
            first_rank = rank;
            first_unique = unique_candidates.size();

            u64 found = 0;
            first_top10 = find_factor_from_candidates(
                data.N,
                unique_candidates,
                sieve,
                10,
                found
            );

            found = 0;
            first_top100 = find_factor_from_candidates(
                data.N,
                unique_candidates,
                sieve,
                100,
                found
            );

            found = 0;
            first_top500 = find_factor_from_candidates(
                data.N,
                unique_candidates,
                sieve,
                500,
                found
            );
        }

        for (std::size_t cut_index = 0;
             cut_index < CUTS.size();
             ++cut_index) {

            std::size_t cut = CUTS[cut_index];

            u64 found_factor = 0;

            u64 prime_tests = find_factor_from_candidates(
                data.N,
                unique_candidates,
                sieve,
                cut,
                found_factor
            );

            stats[cut_index].total_prime_tests += prime_tests;

            if (prime_tests > stats[cut_index].max_prime_tests) {
                stats[cut_index].max_prime_tests = prime_tests;
            }

            if (found_factor != 0) {
                ++stats[cut_index].success;

                if (found_factor == data.p) {
                    ++stats[cut_index].matching;
                }
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

    const auto end = std::chrono::steady_clock::now();

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    std::cout << "FIRST_N=" << first_N << "\n";
    std::cout << "FIRST_P=" << first_p << "\n";
    std::cout << "FIRST_Q=" << first_q << "\n";
    std::cout << "FIRST_S=" << first_s << "\n";
    std::cout << "FIRST_ORDINARY_GCD_TESTS=" << first_ordinary << "\n";
    std::cout << "FIRST_FACTOR_RANK=" << first_rank << "\n";
    std::cout << "FIRST_UNIQUE_CANDIDATES=" << first_unique << "\n";
    std::cout << "FIRST_TOP10_PRIME_TESTS=" << first_top10 << "\n";
    std::cout << "FIRST_TOP100_PRIME_TESTS=" << first_top100 << "\n";
    std::cout << "FIRST_TOP500_PRIME_TESTS=" << first_top500 << "\n";

    std::cout << "CASE_COUNT=" << CASE_COUNT << "\n";

    std::cout << "ORDINARY_TOTAL_GCD_TESTS="
              << ordinary_total_tests << "\n";

    std::cout << "ORDINARY_AVG_GCD_TESTS="
              << average_u64(ordinary_total_tests, CASE_COUNT)
              << "\n";

    std::cout << "ORDINARY_AVG_DISTANCE="
              << average_u64(ordinary_total_distance, CASE_COUNT)
              << "\n";

    std::cout << "TOTAL_UNIQUE_CANDIDATES="
              << all_unique_total << "\n";

    std::cout << "AVG_UNIQUE_CANDIDATES="
              << average_u64(all_unique_total, CASE_COUNT)
              << "\n";

    std::cout << "MAX_UNIQUE_CANDIDATES="
              << all_unique_max << "\n";

    std::cout << "FACTOR_RANK_MISSING="
              << rank_missing << "\n";

    if (CASE_COUNT - rank_missing > 0) {
        std::cout << "AVG_FACTOR_RANK="
                  << static_cast<double>(rank_sum) /
                         static_cast<double>(CASE_COUNT - rank_missing)
                  << "\n";
    } else {
        std::cout << "AVG_FACTOR_RANK=0\n";
    }

    std::cout << "MAX_FACTOR_RANK="
              << rank_max << "\n";

    for (std::size_t i = 0; i < CUTS.size(); ++i) {
        std::cout
            << "CUT="
            << CUTS[i]
            << " SUCCESS="
            << stats[i].success
            << " MATCHING="
            << stats[i].matching
            << " AVG_PRIME_TESTS="
            << average_u64(
                   stats[i].total_prime_tests,
                   CASE_COUNT
               )
            << " MAX_PRIME_TESTS="
            << stats[i].max_prime_tests
            << "\n";
    }

    std::cout << "ELAPSED_TIME_MS="
              << elapsed_ms
              << "\n";

    std::cout << "FINISHED EXPERIMENT 479\n";
}

int main() {
    main_experiment();
    return 0;
}
