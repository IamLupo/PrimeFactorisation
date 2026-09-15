#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 477;

constexpr int PRIME_LIMIT = 100000;
constexpr int PRIME_MIN = 10000;

constexpr int CASE_COUNT = 500;

constexpr int K_LIMIT = 1000;
constexpr int M_LIMIT = 7;

struct PrimeCase {
    u64 p;
    u64 q;
    u64 n;
};

struct CRTPair {
    int k1;
    int k2;
};

struct Representation {
    bool found = false;

    u64 prime = 0;
    u64 k1 = 0;
    u64 k2 = 0;

    u64 a = 0;
    u64 m = 0;
    u64 j = 0;
};

struct StructuredResult {
    bool found = false;

    u64 factor = 0;

    u64 pool_tests = 0;
    u64 progression_tests = 0;

    u64 one_sided_candidates = 0;
    u64 two_sided_candidates = 0;
    u64 certificate_tests = 0;

    double time_ms = 0.0;
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

        cases.push_back(c);
    }

    return cases;
}

std::vector<CRTPair> build_determinant_one_pairs() {
    std::vector<CRTPair> pairs;

    for (int k1 = 1; k1 <= K_LIMIT; ++k1) {
        for (int k2 = k1; k2 <= K_LIMIT; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            bool valid = false;

            for (int m1 = 1; m1 <= M_LIMIT; ++m1) {
                const int numerator =
                    1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const int m2 =
                    numerator / k1;

                if (
                    m2 >= 1 &&
                    m2 <= M_LIMIT
                ) {
                    valid = true;
                    break;
                }
            }

            if (valid) {
                pairs.push_back({k1, k2});
            }
        }
    }

    return pairs;
}

Representation find_min_j_representation(
    u64 prime,
    const std::vector<CRTPair>& pairs,
    u64& pair_tests
) {
    Representation best;
    best.prime = prime;

    for (const CRTPair& pair : pairs) {
        ++pair_tests;

        const u64 k1 =
            static_cast<u64>(pair.k1);

        const u64 k2 =
            static_cast<u64>(pair.k2);

        const u64 a =
            k1 + k2;

        const u64 m =
            k1 * k2;

        if (a > prime) {
            continue;
        }

        const u64 remainder =
            prime - a;

        if (remainder % m != 0) {
            continue;
        }

        const u64 j =
            remainder / m;

        bool take = false;

        if (!best.found) {
            take = true;
        } else if (j < best.j) {
            take = true;
        } else if (
            j == best.j &&
            m < best.m
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) <
                std::max(
                    best.k1,
                    best.k2
                )
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) ==
                std::max(
                    best.k1,
                    best.k2
                ) &&
            k1 < best.k1
        ) {
            take = true;
        }

        if (!take) {
            continue;
        }

        best.found = true;
        best.k1 = k1;
        best.k2 = k2;
        best.a = a;
        best.m = m;
        best.j = j;
    }

    return best;
}

bool one_sided_match(
    u64 n,
    const Representation& p,
    u64 candidate
) {
    if (p.m == 0) {
        return false;
    }

    return
        (
            (p.a % p.m) *
            (candidate % p.m)
        ) %
        p.m
        ==
        n % p.m;
}

bool two_sided_match(
    u64 n,
    const Representation& p,
    const Representation& r
) {
    if (!one_sided_match(
            n,
            p,
            r.prime
        )) {
        return false;
    }

    if (r.m == 0) {
        return false;
    }

    return
        (
            (r.a % r.m) *
            (p.prime % r.m)
        ) %
        r.m
        ==
        n % r.m;
}

bool forced_zero(
    const PrimeCase& c,
    const Representation& p,
    const Representation& r
) {
    const u64 g =
        std::gcd(
            p.m,
            r.m
        );

    if (g == 0) {
        return false;
    }

    const u64 lcm_m =
        (p.m / g) * r.m;

    const u64 product =
        p.prime * r.prime;

    const u64 difference =
        c.n >= product
            ? c.n - product
            : product - c.n;

    return
        difference < lcm_m &&
        difference % lcm_m == 0;
}

StructuredResult structured_pool_scan(
    const PrimeCase& c,
    const Representation& p,
    const std::vector<int>& primes,
    const std::vector<Representation>& cache
) {
    StructuredResult result;

    const auto start =
        std::chrono::steady_clock::now();

    for (int prime : primes) {
        if (
            prime < PRIME_MIN ||
            prime > PRIME_LIMIT
        ) {
            continue;
        }

        const u64 r =
            static_cast<u64>(prime);

        if (r == c.p) {
            continue;
        }

        ++result.pool_tests;

        if (
            !one_sided_match(
                c.n,
                p,
                r
            )
        ) {
            continue;
        }

        ++result.one_sided_candidates;

        const Representation& rr =
            cache[prime];

        if (!rr.found) {
            continue;
        }

        if (
            !two_sided_match(
                c.n,
                p,
                rr
            )
        ) {
            continue;
        }

        ++result.two_sided_candidates;
        ++result.certificate_tests;

        if (
            forced_zero(
                c,
                p,
                rr
            )
        ) {
            result.found = true;
            result.factor = r;
            break;
        }
    }

    const auto end =
        std::chrono::steady_clock::now();

    result.time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(
            end - start
        ).count();

    return result;
}

StructuredResult structured_progression_scan(
    const PrimeCase& c,
    const Representation& p,
    const std::vector<int>& is_prime,
    const std::vector<Representation>& cache
) {
    StructuredResult result;

    const auto start =
        std::chrono::steady_clock::now();

    const u64 modulus =
        p.m;

    if (modulus == 0) {
        return result;
    }

    /*
     * Solve:
     *
     *     a_p*r == N (mod M_p)
     *
     * We know gcd(a_p,M_p)=1.
     *
     * Find the inverse using a small brute-force
     * search here. M_p <= 1,000,000.
     */
    u64 inverse = 0;
    bool found_inverse = false;

    for (u64 x = 1; x < modulus; ++x) {
        if (
            (p.a * x) %
            modulus == 1
        ) {
            inverse = x;
            found_inverse = true;
            break;
        }
    }

    if (!found_inverse) {
        return result;
    }

    const u64 residue =
        (
            (c.n % modulus) *
            inverse
        ) %
        modulus;

    /*
     * First positive integer in range:
     *
     *     r = residue + t*M_p
     */
    u64 first = residue;

    if (first < PRIME_MIN) {
        const u64 delta =
            static_cast<u64>(
                PRIME_MIN
            ) - first;

        const u64 steps =
            (
                delta +
                modulus -
                1
            ) /
            modulus;

        first +=
            steps * modulus;
    }

    for (
        u64 r = first;
        r <=
            static_cast<u64>(
                PRIME_LIMIT
            );
        r += modulus
    ) {
        ++result.progression_tests;

        if (r == c.p) {
            continue;
        }

        if (
            !is_prime[
                static_cast<std::size_t>(r)
            ]
        ) {
            continue;
        }

        ++result.one_sided_candidates;

        const Representation& rr =
            cache[
                static_cast<std::size_t>(r)
            ];

        if (!rr.found) {
            continue;
        }

        if (
            !two_sided_match(
                c.n,
                p,
                rr
            )
        ) {
            continue;
        }

        ++result.two_sided_candidates;
        ++result.certificate_tests;

        if (
            forced_zero(
                c,
                p,
                rr
            )
        ) {
            result.found = true;
            result.factor = r;
            break;
        }
    }

    const auto end =
        std::chrono::steady_clock::now();

    result.time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(
            end - start
        ).count();

    return result;
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47720260915ULL
    );

    const auto total_start =
        std::chrono::steady_clock::now();

    const std::vector<int> primes =
        generate_primes(
            PRIME_LIMIT
        );

    std::vector<int> is_prime(
        PRIME_LIMIT + 1,
        0
    );

    for (int p : primes) {
        if (p <= PRIME_LIMIT) {
            is_prime[p] = 1;
        }
    }

    const std::vector<PrimeCase> cases =
        generate_cases(
            primes,
            rng
        );

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    /*
     * Cache all prime representations.
     */
    std::vector<Representation> cache(
        PRIME_LIMIT + 1
    );

    u64 preprocessing_pair_tests = 0;

    const auto prep_start =
        std::chrono::steady_clock::now();

    for (int prime : primes) {
        if (prime < PRIME_MIN) {
            continue;
        }

        cache[prime] =
            find_min_j_representation(
                static_cast<u64>(prime),
                pairs,
                preprocessing_pair_tests
            );
    }

    const auto prep_end =
        std::chrono::steady_clock::now();

    const double preprocessing_time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(
            prep_end -
            prep_start
        ).count();

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << '\n';

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << '\n';

    /*
     * Pool scan statistics.
     */
    u64 pool_success = 0;
    u64 pool_candidates = 0;
    u64 pool_two_sided = 0;
    u64 pool_certificate_tests = 0;

    double pool_time_ms = 0.0;

    /*
     * Progression scan statistics.
     */
    u64 progression_success = 0;
    u64 progression_integer_tests = 0;
    u64 progression_prime_candidates = 0;
    u64 progression_two_sided = 0;
    u64 progression_certificate_tests = 0;

    double progression_time_ms = 0.0;

    /*
     * Correctness comparison.
     */
    u64 matching_factors = 0;

    /*
     * Reduction statistics.
     */
    long double integer_scan_reduction_sum = 0.0L;
    long double candidate_test_reduction_sum = 0.0L;

    u64 max_progression_tests = 0;
    u64 max_progression_candidates = 0;

    bool first_example = false;

    for (
        int i = 0;
        i < CASE_COUNT;
        ++i
    ) {
        const PrimeCase& c =
            cases[i];

        const Representation& p =
            cache[c.p];

        if (!p.found) {
            continue;
        }

        const StructuredResult pool =
            structured_pool_scan(
                c,
                p,
                primes,
                cache
            );

        const StructuredResult progression =
            structured_progression_scan(
                c,
                p,
                is_prime,
                cache
            );

        if (pool.found) {
            ++pool_success;
        }

        if (progression.found) {
            ++progression_success;
        }

        if (
            pool.found &&
            progression.found &&
            pool.factor ==
                progression.factor
        ) {
            ++matching_factors;
        }

        pool_candidates +=
            pool.one_sided_candidates;

        pool_two_sided +=
            pool.two_sided_candidates;

        pool_certificate_tests +=
            pool.certificate_tests;

        progression_integer_tests +=
            progression.progression_tests;

        progression_prime_candidates +=
            progression.one_sided_candidates;

        progression_two_sided +=
            progression.two_sided_candidates;

        progression_certificate_tests +=
            progression.certificate_tests;

        pool_time_ms +=
            pool.time_ms;

        progression_time_ms +=
            progression.time_ms;

        max_progression_tests =
            std::max(
                max_progression_tests,
                progression.progression_tests
            );

        max_progression_candidates =
            std::max(
                max_progression_candidates,
                progression.one_sided_candidates
            );

        if (
            pool.pool_tests != 0
        ) {
            integer_scan_reduction_sum +=
                1.0L -
                static_cast<long double>(
                    progression.progression_tests
                ) /
                static_cast<long double>(
                    pool.pool_tests
                );
        }

        if (
            pool.one_sided_candidates != 0
        ) {
            candidate_test_reduction_sum +=
                1.0L -
                static_cast<long double>(
                    progression.one_sided_candidates
                ) /
                static_cast<long double>(
                    pool.one_sided_candidates
                );
        }

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
                << "FIRST_P_M="
                << p.m
                << '\n';

            std::cout
                << "FIRST_POOL_TESTS="
                << pool.pool_tests
                << '\n';

            std::cout
                << "FIRST_POOL_ONE_SIDED="
                << pool.one_sided_candidates
                << '\n';

            std::cout
                << "FIRST_POOL_TWO_SIDED="
                << pool.two_sided_candidates
                << '\n';

            std::cout
                << "FIRST_PROGRESSION_TESTS="
                << progression.progression_tests
                << '\n';

            std::cout
                << "FIRST_PROGRESSION_PRIME_CANDIDATES="
                << progression.one_sided_candidates
                << '\n';

            std::cout
                << "FIRST_PROGRESSION_TWO_SIDED="
                << progression.two_sided_candidates
                << '\n';

            std::cout
                << "FIRST_PROGRESSION_FACTOR="
                << progression.factor
                << '\n';

            first_example = true;
        }

        if (
            (i + 1) % 100 ==
            0
        ) {
            std::cout
                << "PROGRESS="
                << (i + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const auto total_end =
        std::chrono::steady_clock::now();

    const double total_time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(
            total_end -
            total_start
        ).count();

    const double avg_pool_candidates =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    pool_candidates
                ) /
                CASE_COUNT
            );

    const double avg_pool_two_sided =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    pool_two_sided
                ) /
                CASE_COUNT
            );

    const double avg_progression_tests =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    progression_integer_tests
                ) /
                CASE_COUNT
            );

    const double avg_progression_candidates =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    progression_prime_candidates
                ) /
                CASE_COUNT
            );

    const double avg_progression_two_sided =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    progression_two_sided
                ) /
                CASE_COUNT
            );

    const double avg_integer_reduction =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                100.0L *
                integer_scan_reduction_sum /
                CASE_COUNT
            );

    const double avg_candidate_reduction =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                100.0L *
                candidate_test_reduction_sum /
                CASE_COUNT
            );

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "POOL_SUCCESS="
        << pool_success
        << '\n';

    std::cout
        << "PROGRESSION_SUCCESS="
        << progression_success
        << '\n';

    std::cout
        << "MATCHING_FACTORS="
        << matching_factors
        << '\n';

    std::cout
        << "POOL_TOTAL_CANDIDATES="
        << pool_candidates
        << '\n';

    std::cout
        << "POOL_AVG_CANDIDATES="
        << avg_pool_candidates
        << '\n';

    std::cout
        << "POOL_TOTAL_TWO_SIDED="
        << pool_two_sided
        << '\n';

    std::cout
        << "POOL_AVG_TWO_SIDED="
        << avg_pool_two_sided
        << '\n';

    std::cout
        << "PROGRESSION_TOTAL_INTEGER_TESTS="
        << progression_integer_tests
        << '\n';

    std::cout
        << "PROGRESSION_AVG_INTEGER_TESTS="
        << avg_progression_tests
        << '\n';

    std::cout
        << "PROGRESSION_MAX_INTEGER_TESTS="
        << max_progression_tests
        << '\n';

    std::cout
        << "PROGRESSION_TOTAL_PRIME_CANDIDATES="
        << progression_prime_candidates
        << '\n';

    std::cout
        << "PROGRESSION_AVG_PRIME_CANDIDATES="
        << avg_progression_candidates
        << '\n';

    std::cout
        << "PROGRESSION_MAX_PRIME_CANDIDATES="
        << max_progression_candidates
        << '\n';

    std::cout
        << "PROGRESSION_TOTAL_TWO_SIDED="
        << progression_two_sided
        << '\n';

    std::cout
        << "PROGRESSION_AVG_TWO_SIDED="
        << avg_progression_two_sided
        << '\n';

    std::cout
        << "PROGRESSION_CERTIFICATE_TESTS="
        << progression_certificate_tests
        << '\n';

    std::cout
        << "AVG_INTEGER_SCAN_REDUCTION_PERCENT="
        << avg_integer_reduction
        << '\n';

    std::cout
        << "AVG_PRIME_CANDIDATE_REDUCTION_PERCENT="
        << avg_candidate_reduction
        << '\n';

    std::cout
        << "PREPROCESSING_PAIR_TESTS="
        << preprocessing_pair_tests
        << '\n';

    std::cout
        << "PREPROCESSING_TIME_MS="
        << preprocessing_time_ms
        << '\n';

    std::cout
        << "POOL_QUERY_TIME_MS="
        << pool_time_ms
        << '\n';

    std::cout
        << "PROGRESSION_QUERY_TIME_MS="
        << progression_time_ms
        << '\n';

    std::cout
        << "PROGRESSION_TOTAL_TIME_MS="
        << (
            preprocessing_time_ms +
            progression_time_ms
        )
        << '\n';

    std::cout
        << "POOL_TOTAL_TIME_MS="
        << (
            preprocessing_time_ms +
            pool_time_ms
        )
        << '\n';

    std::cout
        << "PROGRESSION_OVER_POOL_QUERY="
        << (
            pool_time_ms == 0.0
                ? 0.0
                : progression_time_ms /
                  pool_time_ms
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
