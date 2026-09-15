#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_set>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 464;

constexpr int PRIME_LIMIT = 100000;
constexpr int PRIME_MIN = 10000;
constexpr int CASE_COUNT = 500;

constexpr int K_LIMIT = 1000;
constexpr int M_LIMIT = 7;

struct PrimeCase {
    u64 p;
    u64 q;
    u64 n;
    u64 s;
};

struct CRTPair {
    int k1;
    int k2;
    u64 base;
    u64 modulus;
};

struct CRTNode {
    u64 value;
    int pair_index;

    bool operator<(const CRTNode& other) const {
        return value < other.value;
    }
};

struct ScanResult {
    u64 factor = 0;
    u64 gcd_tests = 0;
    u64 distance = 0;
    bool found = false;
};

struct CRTScanResult {
    u64 factor = 0;
    u64 gcd_tests = 0;
    u64 candidate_count = 0;
    bool found = false;
};

std::vector<int> generate_primes(int limit) {
    std::vector<bool> composite(limit + 1, false);
    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<std::int64_t>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

std::vector<CRTPair> build_crt_pairs() {
    std::vector<CRTPair> pairs;

    for (int k1 = 1; k1 <= K_LIMIT; ++k1) {
        for (int k2 = 1; k2 <= K_LIMIT; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            /*
             * We need a determinant-one pair
             *
             *     m2*k1 - m1*k2 = 1
             *
             * with 1 <= m1,m2 <= M_LIMIT.
             *
             * For the CRT progression itself, only (k1,k2) matters:
             *
             *     p = k1 + k2 + j*k1*k2
             */
            bool valid = false;

            for (int m1 = 1; m1 <= M_LIMIT; ++m1) {
                const int numerator = 1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const int m2 = numerator / k1;

                if (m2 >= 1 && m2 <= M_LIMIT) {
                    valid = true;
                    break;
                }
            }

            if (!valid) {
                continue;
            }

            CRTPair pair;
            pair.k1 = k1;
            pair.k2 = k2;
            pair.base = static_cast<u64>(k1) + static_cast<u64>(k2);
            pair.modulus =
                static_cast<u64>(k1) * static_cast<u64>(k2);

            pairs.push_back(pair);
        }
    }

    std::sort(
        pairs.begin(),
        pairs.end(),
        [](const CRTPair& a, const CRTPair& b) {
            if (a.base != b.base) {
                return a.base < b.base;
            }
            return a.modulus < b.modulus;
        }
    );

    pairs.erase(
        std::unique(
            pairs.begin(),
            pairs.end(),
            [](const CRTPair& a, const CRTPair& b) {
                return a.base == b.base &&
                       a.modulus == b.modulus;
            }
        ),
        pairs.end()
    );

    return pairs;
}

std::vector<PrimeCase> generate_cases(
    const std::vector<int>& primes,
    std::mt19937_64& rng
) {
    std::vector<int> candidates;

    for (int p : primes) {
        if (p >= PRIME_MIN && p <= PRIME_LIMIT) {
            candidates.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        candidates.size() - 1
    );

    std::vector<PrimeCase> cases;
    cases.reserve(CASE_COUNT);

    for (int i = 0; i < CASE_COUNT; ++i) {
        int p = candidates[dist(rng)];
        int q = candidates[dist(rng)];

        while (q == p) {
            q = candidates[dist(rng)];
        }

        if (p > q) {
            std::swap(p, q);
        }

        PrimeCase c;
        c.p = static_cast<u64>(p);
        c.q = static_cast<u64>(q);
        c.n = c.p * c.q;
        c.s = static_cast<u64>(
            std::sqrt(static_cast<long double>(c.n))
        );

        while ((c.s + 1) * (c.s + 1) <= c.n) {
            ++c.s;
        }

        while (c.s * c.s > c.n) {
            --c.s;
        }

        cases.push_back(c);
    }

    return cases;
}

ScanResult ordinary_scan(const PrimeCase& c) {
    ScanResult result;

    for (u64 p_candidate = c.s; p_candidate >= 2; --p_candidate) {
        ++result.gcd_tests;

        const u64 g = std::gcd(c.n, p_candidate);

        if (g > 1 && g < c.n) {
            result.factor = g;
            result.distance = c.s - p_candidate;
            result.found = true;
            return result;
        }

        if (p_candidate == 2) {
            break;
        }
    }

    return result;
}

CRTScanResult crt_filtered_scan(
    const PrimeCase& c,
    const std::vector<CRTPair>& pairs
) {
    CRTScanResult result;

    std::priority_queue<CRTNode> heap;

    /*
     * For every CRT progression
     *
     *     p_j = base + j*modulus
     *
     * find the largest member <= s.
     */
    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        const CRTPair& pair = pairs[i];

        if (pair.base > c.s) {
            continue;
        }

        const u64 j =
            (c.s - pair.base) / pair.modulus;

        const u64 value =
            pair.base + j * pair.modulus;

        if (value < 2) {
            continue;
        }

        heap.push({value, i});
    }

    u64 last_value = 0;
    bool have_last_value = false;

    while (!heap.empty()) {
        CRTNode node = heap.top();
        heap.pop();

        const CRTPair& pair = pairs[node.pair_index];

        /*
         * Several CRT classes can produce the same candidate.
         * Only test each integer once.
         */
        if (!have_last_value || node.value != last_value) {
            last_value = node.value;
            have_last_value = true;

            ++result.candidate_count;
            ++result.gcd_tests;

            const u64 g = std::gcd(c.n, node.value);

            if (g > 1 && g < c.n) {
                result.factor = g;
                result.found = true;
                return result;
            }
        }

        /*
         * Move this arithmetic progression down by one step.
         */
        if (node.value > pair.base) {
            const u64 next_value =
                node.value - pair.modulus;

            if (next_value >= 2) {
                heap.push({
                    next_value,
                    node.pair_index
                });
            }
        }
    }

    return result;
}

void print_pair_summary(
    const std::vector<CRTPair>& pairs
) {
    std::cout
        << "CRT_PAIR_COUNT="
        << pairs.size()
        << '\n';

    int max_k1 = 0;
    int max_k2 = 0;

    for (const CRTPair& pair : pairs) {
        max_k1 = std::max(max_k1, pair.k1);
        max_k2 = std::max(max_k2, pair.k2);
    }

    std::cout
        << "CRT_MAX_K1="
        << max_k1
        << '\n';

    std::cout
        << "CRT_MAX_K2="
        << max_k2
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x46420260915ULL);

    const std::vector<int> primes =
        generate_primes(PRIME_LIMIT);

    const std::vector<PrimeCase> cases =
        generate_cases(primes, rng);

    const std::vector<CRTPair> crt_pairs =
        build_crt_pairs();

    print_pair_summary(crt_pairs);

    u64 ordinary_total_gcd = 0;
    u64 ordinary_total_distance = 0;

    u64 crt_total_gcd = 0;
    u64 crt_total_candidates = 0;

    u64 ordinary_hits = 0;
    u64 crt_hits = 0;
    u64 matching_factors = 0;
    u64 crt_misses = 0;

    double ordinary_time_ms = 0.0;
    double crt_time_ms = 0.0;

    u64 ordinary_max_distance = 0;
    u64 crt_max_candidates = 0;

    for (int i = 0; i < CASE_COUNT; ++i) {
        const PrimeCase& c = cases[i];

        const auto ordinary_start =
            std::chrono::steady_clock::now();

        const ScanResult ordinary =
            ordinary_scan(c);

        const auto ordinary_end =
            std::chrono::steady_clock::now();

        const auto crt_start =
            std::chrono::steady_clock::now();

        const CRTScanResult crt =
            crt_filtered_scan(c, crt_pairs);

        const auto crt_end =
            std::chrono::steady_clock::now();

        ordinary_time_ms +=
            std::chrono::duration<double, std::milli>(
                ordinary_end - ordinary_start
            ).count();

        crt_time_ms +=
            std::chrono::duration<double, std::milli>(
                crt_end - crt_start
            ).count();

        if (ordinary.found) {
            ++ordinary_hits;
            ordinary_total_gcd += ordinary.gcd_tests;
            ordinary_total_distance += ordinary.distance;
            ordinary_max_distance =
                std::max(
                    ordinary_max_distance,
                    ordinary.distance
                );
        }

        if (crt.found) {
            ++crt_hits;
            crt_total_gcd += crt.gcd_tests;
            crt_total_candidates += crt.candidate_count;
            crt_max_candidates =
                std::max(
                    crt_max_candidates,
                    crt.candidate_count
                );
        } else {
            ++crt_misses;
        }

        if (ordinary.found && crt.found &&
            ordinary.factor == crt.factor) {
            ++matching_factors;
        }

        if ((i + 1) % 100 == 0) {
            std::cout
                << "PROGRESS="
                << (i + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const double ordinary_avg_gcd =
        ordinary_hits == 0
            ? 0.0
            : static_cast<double>(ordinary_total_gcd) /
              static_cast<double>(ordinary_hits);

    const double ordinary_avg_distance =
        ordinary_hits == 0
            ? 0.0
            : static_cast<double>(ordinary_total_distance) /
              static_cast<double>(ordinary_hits);

    const double crt_avg_gcd =
        crt_hits == 0
            ? 0.0
            : static_cast<double>(crt_total_gcd) /
              static_cast<double>(crt_hits);

    const double crt_avg_candidates =
        crt_hits == 0
            ? 0.0
            : static_cast<double>(crt_total_candidates) /
              static_cast<double>(crt_hits);

    const double gcd_reduction =
        ordinary_total_gcd == 0
            ? 0.0
            : 100.0 *
              (1.0 -
               static_cast<double>(crt_total_gcd) /
               static_cast<double>(ordinary_total_gcd));

    const double distance_reduction =
        ordinary_total_distance == 0
            ? 0.0
            : 100.0 *
              (1.0 -
               static_cast<double>(crt_total_candidates) /
               static_cast<double>(ordinary_total_distance));

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "ORDINARY_HITS="
        << ordinary_hits
        << '\n';

    std::cout
        << "CRT_HITS="
        << crt_hits
        << '\n';

    std::cout
        << "CRT_MISSES="
        << crt_misses
        << '\n';

    std::cout
        << "MATCHING_FACTORS="
        << matching_factors
        << '\n';

    std::cout
        << "ORDINARY_TOTAL_GCD_TESTS="
        << ordinary_total_gcd
        << '\n';

    std::cout
        << "CRT_TOTAL_GCD_TESTS="
        << crt_total_gcd
        << '\n';

    std::cout
        << "ORDINARY_AVG_GCD_TESTS="
        << ordinary_avg_gcd
        << '\n';

    std::cout
        << "CRT_AVG_GCD_TESTS="
        << crt_avg_gcd
        << '\n';

    std::cout
        << "ORDINARY_TOTAL_DISTANCE="
        << ordinary_total_distance
        << '\n';

    std::cout
        << "CRT_TOTAL_UNIQUE_CANDIDATES="
        << crt_total_candidates
        << '\n';

    std::cout
        << "ORDINARY_AVG_DISTANCE="
        << ordinary_avg_distance
        << '\n';

    std::cout
        << "CRT_AVG_CANDIDATES="
        << crt_avg_candidates
        << '\n';

    std::cout
        << "ORDINARY_MAX_DISTANCE="
        << ordinary_max_distance
        << '\n';

    std::cout
        << "CRT_MAX_CANDIDATES="
        << crt_max_candidates
        << '\n';

    std::cout
        << "GCD_TEST_REDUCTION_PERCENT="
        << gcd_reduction
        << '\n';

    std::cout
        << "CANDIDATE_REDUCTION_PERCENT="
        << distance_reduction
        << '\n';

    std::cout
        << "ORDINARY_TIME_MS="
        << ordinary_time_ms
        << '\n';

    std::cout
        << "CRT_TIME_MS="
        << crt_time_ms
        << '\n';

    std::cout
        << "TIME_RATIO_CRT_OVER_ORDINARY="
        << (
            ordinary_time_ms == 0.0
                ? 0.0
                : crt_time_ms / ordinary_time_ms
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';
}

int main() {
    main_experiment();
    return 0;
}
