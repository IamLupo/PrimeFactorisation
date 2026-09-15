#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 478;

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
    u64 k1;
    u64 k2;
    u64 base;
    u64 modulus;
};

struct HeapNode {
    u64 value;
    std::size_t pair_index;

    bool operator<(const HeapNode& other) const {
        return value < other.value;
    }
};

struct OrdinaryResult {
    bool found = false;
    u64 factor = 0;
    u64 gcd_tests = 0;
    u64 distance = 0;
};

struct StructuralResult {
    bool found = false;
    u64 factor = 0;

    u64 progression_values = 0;
    u64 unique_values = 0;
    u64 prime_candidates = 0;
    u64 factor_tests = 0;
    u64 duplicate_values = 0;
};

std::vector<int> generate_primes(int limit) {
    std::vector<bool> composite(
        limit + 1,
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (
            static_cast<std::int64_t>(i) * i <=
            limit
        ) {
            for (
                int j = i * i;
                j <= limit;
                j += i
            ) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

std::vector<PrimeCase> generate_cases(
    const std::vector<int>& primes,
    std::mt19937_64& rng
) {
    std::vector<int> candidates;

    for (int prime : primes) {
        if (
            prime >= PRIME_MIN &&
            prime <= PRIME_LIMIT
        ) {
            candidates.push_back(prime);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        candidates.size() - 1
    );

    std::vector<PrimeCase> cases;
    cases.reserve(CASE_COUNT);

    for (int i = 0; i < CASE_COUNT; ++i) {
        u64 p =
            static_cast<u64>(
                candidates[dist(rng)]
            );

        u64 q =
            static_cast<u64>(
                candidates[dist(rng)]
            );

        while (q == p) {
            q =
                static_cast<u64>(
                    candidates[dist(rng)]
                );
        }

        if (p > q) {
            std::swap(p, q);
        }

        PrimeCase c;

        c.p = p;
        c.q = q;
        c.n = p * q;

        c.s = static_cast<u64>(
            std::sqrt(
                static_cast<long double>(c.n)
            )
        );

        while (
            (c.s + 1) *
            (c.s + 1) <=
            c.n
        ) {
            ++c.s;
        }

        while (c.s * c.s > c.n) {
            --c.s;
        }

        cases.push_back(c);
    }

    return cases;
}

std::vector<CRTPair> build_determinant_one_pairs() {
    std::vector<CRTPair> pairs;

    for (u64 k1 = 1; k1 <= K_LIMIT; ++k1) {
        for (u64 k2 = k1; k2 <= K_LIMIT; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            bool valid = false;

            for (u64 m1 = 1; m1 <= M_LIMIT; ++m1) {
                const u64 numerator =
                    1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const u64 m2 =
                    numerator / k1;

                if (
                    m2 >= 1 &&
                    m2 <= M_LIMIT
                ) {
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
            pair.base = k1 + k2;
            pair.modulus = k1 * k2;

            pairs.push_back(pair);
        }
    }

    return pairs;
}

OrdinaryResult ordinary_scan(
    const PrimeCase& c
) {
    OrdinaryResult result;

    for (
        u64 candidate = c.s;
        candidate >= 2;
        --candidate
    ) {
        ++result.gcd_tests;

        const u64 g =
            std::gcd(
                c.n,
                candidate
            );

        if (
            g > 1 &&
            g < c.n
        ) {
            result.found = true;
            result.factor = g;
            result.distance =
                c.s - candidate;

            return result;
        }

        if (candidate == 2) {
            break;
        }
    }

    return result;
}

void initialize_heap(
    const PrimeCase& c,
    const std::vector<CRTPair>& pairs,
    std::priority_queue<HeapNode>& heap
) {
    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        const CRTPair& pair = pairs[i];

        if (pair.base > c.s) {
            continue;
        }

        const u64 j =
            (c.s - pair.base) /
            pair.modulus;

        const u64 value =
            pair.base +
            j * pair.modulus;

        if (value < 2) {
            continue;
        }

        heap.push({
            value,
            i
        });
    }
}

StructuralResult structural_scan(
    const PrimeCase& c,
    const std::vector<CRTPair>& pairs,
    const std::vector<bool>& is_prime
) {
    StructuralResult result;

    std::priority_queue<HeapNode> heap;

    initialize_heap(
        c,
        pairs,
        heap
    );

    u64 last_value = 0;
    bool have_last_value = false;

    while (!heap.empty()) {
        const HeapNode node =
            heap.top();

        heap.pop();

        ++result.progression_values;

        const CRTPair& pair =
            pairs[node.pair_index];

        const bool duplicate =
            have_last_value &&
            node.value == last_value;

        if (duplicate) {
            ++result.duplicate_values;
        } else {
            last_value = node.value;
            have_last_value = true;

            ++result.unique_values;

            if (
                node.value <=
                static_cast<u64>(PRIME_LIMIT) &&
                is_prime[
                    static_cast<std::size_t>(
                        node.value
                    )
                ]
            ) {
                ++result.prime_candidates;
                ++result.factor_tests;

                if (
                    c.n % node.value == 0
                ) {
                    result.found = true;
                    result.factor =
                        node.value;

                    return result;
                }
            }
        }

        if (node.value > pair.base) {
            const u64 next =
                node.value -
                pair.modulus;

            if (next >= 2) {
                heap.push({
                    next,
                    node.pair_index
                });
            }
        }
    }

    return result;
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47820260915ULL
    );

    const std::vector<int> primes =
        generate_primes(
            PRIME_LIMIT
        );

    std::vector<bool> is_prime(
        PRIME_LIMIT + 1,
        false
    );

    for (int prime : primes) {
        is_prime[
            static_cast<std::size_t>(prime)
        ] = true;
    }

    const std::vector<PrimeCase> cases =
        generate_cases(
            primes,
            rng
        );

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << '\n';

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << '\n';

    u64 ordinary_success = 0;
    u64 structural_success = 0;
    u64 matching_factors = 0;

    u64 ordinary_gcd_tests = 0;
    u64 ordinary_distance = 0;

    u64 structural_progression_values = 0;
    u64 structural_unique_values = 0;
    u64 structural_prime_candidates = 0;
    u64 structural_factor_tests = 0;
    u64 structural_duplicates = 0;

    u64 max_progression_values = 0;
    u64 max_unique_values = 0;
    u64 max_prime_candidates = 0;
    u64 max_factor_tests = 0;

    double ordinary_time_ms = 0.0;
    double structural_time_ms = 0.0;

    bool first_example = false;

    for (
        int case_index = 0;
        case_index < CASE_COUNT;
        ++case_index
    ) {
        const PrimeCase& c =
            cases[case_index];

        const auto ordinary_start =
            std::chrono::steady_clock::now();

        const OrdinaryResult ordinary =
            ordinary_scan(c);

        const auto ordinary_end =
            std::chrono::steady_clock::now();

        const auto structural_start =
            std::chrono::steady_clock::now();

        const StructuralResult structural =
            structural_scan(
                c,
                pairs,
                is_prime
            );

        const auto structural_end =
            std::chrono::steady_clock::now();

        ordinary_time_ms +=
            std::chrono::duration<
                double,
                std::milli
            >(
                ordinary_end -
                ordinary_start
            ).count();

        structural_time_ms +=
            std::chrono::duration<
                double,
                std::milli
            >(
                structural_end -
                structural_start
            ).count();

        if (ordinary.found) {
            ++ordinary_success;
        }

        if (structural.found) {
            ++structural_success;
        }

        if (
            ordinary.found &&
            structural.found &&
            ordinary.factor ==
                structural.factor
        ) {
            ++matching_factors;
        }

        ordinary_gcd_tests +=
            ordinary.gcd_tests;

        ordinary_distance +=
            ordinary.distance;

        structural_progression_values +=
            structural.progression_values;

        structural_unique_values +=
            structural.unique_values;

        structural_prime_candidates +=
            structural.prime_candidates;

        structural_factor_tests +=
            structural.factor_tests;

        structural_duplicates +=
            structural.duplicate_values;

        max_progression_values =
            std::max(
                max_progression_values,
                structural.progression_values
            );

        max_unique_values =
            std::max(
                max_unique_values,
                structural.unique_values
            );

        max_prime_candidates =
            std::max(
                max_prime_candidates,
                structural.prime_candidates
            );

        max_factor_tests =
            std::max(
                max_factor_tests,
                structural.factor_tests
            );

        if (!first_example) {
            std::cout
                << "FIRST_N="
                << c.n
                << '\n';

            std::cout
                << "FIRST_P="
                << c.p
                << '\n';

            std::cout
                << "FIRST_Q="
                << c.q
                << '\n';

            std::cout
                << "FIRST_S="
                << c.s
                << '\n';

            std::cout
                << "FIRST_ORDINARY_GCD_TESTS="
                << ordinary.gcd_tests
                << '\n';

            std::cout
                << "FIRST_ORDINARY_DISTANCE="
                << ordinary.distance
                << '\n';

            std::cout
                << "FIRST_STRUCTURAL_PROGRESSION_VALUES="
                << structural.progression_values
                << '\n';

            std::cout
                << "FIRST_STRUCTURAL_UNIQUE_VALUES="
                << structural.unique_values
                << '\n';

            std::cout
                << "FIRST_STRUCTURAL_PRIME_CANDIDATES="
                << structural.prime_candidates
                << '\n';

            std::cout
                << "FIRST_STRUCTURAL_FACTOR_TESTS="
                << structural.factor_tests
                << '\n';

            std::cout
                << "FIRST_STRUCTURAL_FACTOR="
                << structural.factor
                << '\n';

            first_example = true;
        }

        if (
            (case_index + 1) % 100 ==
            0
        ) {
            std::cout
                << "PROGRESS="
                << (case_index + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const double avg_ordinary_gcd =
        static_cast<double>(
            static_cast<long double>(
                ordinary_gcd_tests
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double avg_ordinary_distance =
        static_cast<double>(
            static_cast<long double>(
                ordinary_distance
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double avg_progression =
        static_cast<double>(
            static_cast<long double>(
                structural_progression_values
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double avg_unique =
        static_cast<double>(
            static_cast<long double>(
                structural_unique_values
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double avg_prime_candidates =
        static_cast<double>(
            static_cast<long double>(
                structural_prime_candidates
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double avg_factor_tests =
        static_cast<double>(
            static_cast<long double>(
                structural_factor_tests
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double avg_duplicates =
        static_cast<double>(
            static_cast<long double>(
                structural_duplicates
            ) /
            static_cast<long double>(
                CASE_COUNT
            )
        );

    const double progression_reduction =
        100.0 *
        (
            1.0 -
            static_cast<double>(
                structural_progression_values
            ) /
            static_cast<double>(
                ordinary_gcd_tests
            )
        );

    const double prime_candidate_reduction =
        100.0 *
        (
            1.0 -
            static_cast<double>(
                structural_prime_candidates
            ) /
            static_cast<double>(
                ordinary_gcd_tests
            )
        );

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "ORDINARY_SUCCESS="
        << ordinary_success
        << '\n';

    std::cout
        << "STRUCTURAL_SUCCESS="
        << structural_success
        << '\n';

    std::cout
        << "MATCHING_FACTORS="
        << matching_factors
        << '\n';

    std::cout
        << "ORDINARY_TOTAL_GCD_TESTS="
        << ordinary_gcd_tests
        << '\n';

    std::cout
        << "ORDINARY_AVG_GCD_TESTS="
        << avg_ordinary_gcd
        << '\n';

    std::cout
        << "ORDINARY_AVG_DISTANCE="
        << avg_ordinary_distance
        << '\n';

    std::cout
        << "STRUCTURAL_TOTAL_PROGRESSION_VALUES="
        << structural_progression_values
        << '\n';

    std::cout
        << "STRUCTURAL_AVG_PROGRESSION_VALUES="
        << avg_progression
        << '\n';

    std::cout
        << "STRUCTURAL_MAX_PROGRESSION_VALUES="
        << max_progression_values
        << '\n';

    std::cout
        << "STRUCTURAL_TOTAL_UNIQUE_VALUES="
        << structural_unique_values
        << '\n';

    std::cout
        << "STRUCTURAL_AVG_UNIQUE_VALUES="
        << avg_unique
        << '\n';

    std::cout
        << "STRUCTURAL_MAX_UNIQUE_VALUES="
        << max_unique_values
        << '\n';

    std::cout
        << "STRUCTURAL_TOTAL_PRIME_CANDIDATES="
        << structural_prime_candidates
        << '\n';

    std::cout
        << "STRUCTURAL_AVG_PRIME_CANDIDATES="
        << avg_prime_candidates
        << '\n';

    std::cout
        << "STRUCTURAL_MAX_PRIME_CANDIDATES="
        << max_prime_candidates
        << '\n';

    std::cout
        << "STRUCTURAL_TOTAL_FACTOR_TESTS="
        << structural_factor_tests
        << '\n';

    std::cout
        << "STRUCTURAL_AVG_FACTOR_TESTS="
        << avg_factor_tests
        << '\n';

    std::cout
        << "STRUCTURAL_MAX_FACTOR_TESTS="
        << max_factor_tests
        << '\n';

    std::cout
        << "STRUCTURAL_AVG_DUPLICATE_VALUES="
        << avg_duplicates
        << '\n';

    std::cout
        << "PROGRESSION_VALUE_REDUCTION_PERCENT="
        << progression_reduction
        << '\n';

    std::cout
        << "PRIME_CANDIDATE_REDUCTION_PERCENT="
        << prime_candidate_reduction
        << '\n';

    std::cout
        << "ORDINARY_TIME_MS="
        << ordinary_time_ms
        << '\n';

    std::cout
        << "STRUCTURAL_TIME_MS="
        << structural_time_ms
        << '\n';

    std::cout
        << "STRUCTURAL_OVER_ORDINARY="
        << (
            ordinary_time_ms == 0.0
                ? 0.0
                : structural_time_ms /
                  ordinary_time_ms
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