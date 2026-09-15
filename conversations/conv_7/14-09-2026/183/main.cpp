#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 476;

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

    u64 prime_pool_tests = 0;
    u64 one_sided_candidates = 0;
    u64 two_sided_candidates = 0;
    u64 certificate_tests = 0;

    u64 representation_pair_tests = 0;

    double time_ms = 0.0;
};

struct GcdScanResult {
    bool found = false;

    u64 factor = 0;

    u64 gcd_tests = 0;
    u64 distance = 0;

    double time_ms = 0.0;
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

std::vector<PrimeCase> generate_cases(
    const std::vector<int>& primes,
    std::mt19937_64& rng
) {
    std::vector<int> candidates;

    for (int prime : primes) {
        if (prime >= PRIME_MIN &&
            prime <= PRIME_LIMIT) {
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
                static_cast<long double>(
                    c.n
                )
            )
        );

        while (
            (c.s + 1) * (c.s + 1) <=
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

                if (m2 >= 1 &&
                    m2 <= M_LIMIT) {
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
    if (!one_sided_match(n, p, r.prime)) {
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
        std::gcd(p.m, r.m);

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

StructuredResult structured_factor(
    const PrimeCase& c,
    const Representation& p_rep,
    const std::vector<int>& primes,
    const std::vector<Representation>& cache
) {
    StructuredResult result;

    const auto start =
        std::chrono::steady_clock::now();

    std::vector<int> one_sided;
    one_sided.reserve(512);

    /*
     * Stage 1:
     *
     *   M_p | N - a_p*r
     */
    for (int prime : primes) {
        if (prime < PRIME_MIN ||
            prime > PRIME_LIMIT) {
            continue;
        }

        const u64 r =
            static_cast<u64>(prime);

        if (r == c.p) {
            continue;
        }

        ++result.prime_pool_tests;

        if (
            one_sided_match(
                c.n,
                p_rep,
                r
            )
        ) {
            one_sided.push_back(prime);
        }
    }

    result.one_sided_candidates =
        one_sided.size();

    /*
     * Stage 2:
     *
     * Each candidate supplies its own M_r,a_r.
     */
    for (int prime : one_sided) {
        const Representation& r =
            cache[prime];

        if (!r.found) {
            continue;
        }

        if (
            two_sided_match(
                c.n,
                p_rep,
                r
            )
        ) {
            /*
             * Stage 3:
             *
             * Exact forced-zero certificate.
             */
            ++result.two_sided_candidates;
            ++result.certificate_tests;

            if (
                forced_zero(
                    c,
                    p_rep,
                    r
                )
            ) {
                result.found = true;
                result.factor = r.prime;
                break;
            }
        }
    }

    const auto end =
        std::chrono::steady_clock::now();

    result.time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(end - start).count();

    return result;
}

GcdScanResult ordinary_gcd_scan(
    const PrimeCase& c
) {
    GcdScanResult result;

    const auto start =
        std::chrono::steady_clock::now();

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
            break;
        }

        if (candidate == 2) {
            break;
        }
    }

    const auto end =
        std::chrono::steady_clock::now();

    result.time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(end - start).count();

    return result;
}

void print_representation(
    const char* label,
    const Representation& rep
) {
    std::cout
        << label
        << "_PRIME="
        << rep.prime
        << '\n';

    std::cout
        << label
        << "_K1="
        << rep.k1
        << '\n';

    std::cout
        << label
        << "_K2="
        << rep.k2
        << '\n';

    std::cout
        << label
        << "_A="
        << rep.a
        << '\n';

    std::cout
        << label
        << "_M="
        << rep.m
        << '\n';

    std::cout
        << label
        << "_J="
        << rep.j
        << '\n';

    std::cout
        << label
        << "_CHECK="
        << rep.a +
           rep.j * rep.m
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47620260915ULL
    );

    const auto total_start =
        std::chrono::steady_clock::now();

    const std::vector<int> primes =
        generate_primes(PRIME_LIMIT);

    const std::vector<PrimeCase> cases =
        generate_cases(primes, rng);

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    /*
     * Build representation cache.
     *
     * We count this separately because it is a
     * preprocessing cost.
     */
    std::vector<Representation> cache(
        PRIME_LIMIT + 1
    );

    u64 preprocessing_pair_tests = 0;

    const auto preprocessing_start =
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

    const auto preprocessing_end =
        std::chrono::steady_clock::now();

    const double preprocessing_time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(
            preprocessing_end -
            preprocessing_start
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
     * Correctness.
     */
    u64 ordinary_success = 0;
    u64 structured_success = 0;
    u64 matching_factors = 0;

    /*
     * Ordinary scan totals.
     */
    u64 ordinary_gcd_tests = 0;
    u64 ordinary_distance = 0;

    /*
     * Structured totals.
     */
    u64 structured_pool_tests = 0;
    u64 structured_one_sided = 0;
    u64 structured_two_sided = 0;
    u64 structured_certificate_tests = 0;
    u64 structured_pair_tests = 0;

    /*
     * Maximums.
     */
    u64 max_ordinary_gcd_tests = 0;
    u64 max_one_sided = 0;
    u64 max_two_sided = 0;
    u64 max_certificate_tests = 0;

    /*
     * Timing.
     */
    double ordinary_time_ms = 0.0;
    double structured_time_ms = 0.0;

    /*
     * Speed statistics.
     */
    long double gcd_reduction_sum = 0.0L;
    long double pool_reduction_sum = 0.0L;
    long double candidate_reduction_sum = 0.0L;

    bool first_example = false;

    for (
        int i = 0;
        i < CASE_COUNT;
        ++i
    ) {
        const PrimeCase& c =
            cases[i];

        const Representation& p_rep =
            cache[c.p];

        if (!p_rep.found) {
            continue;
        }

        const GcdScanResult ordinary =
            ordinary_gcd_scan(c);

        const StructuredResult structured =
            structured_factor(
                c,
                p_rep,
                primes,
                cache
            );

        if (ordinary.found) {
            ++ordinary_success;
        }

        if (structured.found) {
            ++structured_success;
        }

        if (
            ordinary.found &&
            structured.found &&
            ordinary.factor ==
                structured.factor
        ) {
            ++matching_factors;
        }

        ordinary_gcd_tests +=
            ordinary.gcd_tests;

        ordinary_distance +=
            ordinary.distance;

        structured_pool_tests +=
            structured.prime_pool_tests;

        structured_one_sided +=
            structured.one_sided_candidates;

        structured_two_sided +=
            structured.two_sided_candidates;

        structured_certificate_tests +=
            structured.certificate_tests;

        structured_pair_tests +=
            structured.representation_pair_tests;

        max_ordinary_gcd_tests =
            std::max(
                max_ordinary_gcd_tests,
                ordinary.gcd_tests
            );

        max_one_sided =
            std::max(
                max_one_sided,
                structured.one_sided_candidates
            );

        max_two_sided =
            std::max(
                max_two_sided,
                structured.two_sided_candidates
            );

        max_certificate_tests =
            std::max(
                max_certificate_tests,
                structured.certificate_tests
            );

        ordinary_time_ms +=
            ordinary.time_ms;

        structured_time_ms +=
            structured.time_ms;

        if (
            ordinary.gcd_tests != 0
        ) {
            gcd_reduction_sum +=
                1.0L -
                static_cast<long double>(
                    structured.certificate_tests
                ) /
                static_cast<long double>(
                    ordinary.gcd_tests
                );
        }

        if (
            structured.prime_pool_tests != 0
        ) {
            pool_reduction_sum +=
                1.0L -
                static_cast<long double>(
                    structured.one_sided_candidates
                ) /
                static_cast<long double>(
                    structured.prime_pool_tests
                );
        }

        if (
            structured.one_sided_candidates != 0
        ) {
            candidate_reduction_sum +=
                1.0L -
                static_cast<long double>(
                    structured.two_sided_candidates
                ) /
                static_cast<long double>(
                    structured.one_sided_candidates
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

            print_representation(
                "FIRST_P_REP",
                p_rep
            );

            std::cout
                << "FIRST_ORDINARY_GCD_TESTS="
                << ordinary.gcd_tests
                << '\n';

            std::cout
                << "FIRST_ORDINARY_DISTANCE="
                << ordinary.distance
                << '\n';

            std::cout
                << "FIRST_STRUCTURED_POOL_TESTS="
                << structured.prime_pool_tests
                << '\n';

            std::cout
                << "FIRST_STRUCTURED_ONE_SIDED="
                << structured.one_sided_candidates
                << '\n';

            std::cout
                << "FIRST_STRUCTURED_TWO_SIDED="
                << structured.two_sided_candidates
                << '\n';

            std::cout
                << "FIRST_STRUCTURED_CERTIFICATE_TESTS="
                << structured.certificate_tests
                << '\n';

            std::cout
                << "FIRST_STRUCTURED_FACTOR="
                << structured.factor
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

    const double avg_ordinary_gcd =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    ordinary_gcd_tests
                ) /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_ordinary_distance =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    ordinary_distance
                ) /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_structured_pool =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    structured_pool_tests
                ) /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_one_sided =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    structured_one_sided
                ) /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_two_sided =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    structured_two_sided
                ) /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_certificate =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    structured_certificate_tests
                ) /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_gcd_reduction =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                100.0L *
                gcd_reduction_sum /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_pool_reduction =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                100.0L *
                pool_reduction_sum /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double avg_candidate_reduction =
        CASE_COUNT == 0
            ? 0.0
            : static_cast<double>(
                100.0L *
                candidate_reduction_sum /
                static_cast<long double>(
                    CASE_COUNT
                )
            );

    const double structured_total_time =
        preprocessing_time_ms +
        structured_time_ms;

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "ORDINARY_SUCCESS="
        << ordinary_success
        << '\n';

    std::cout
        << "STRUCTURED_SUCCESS="
        << structured_success
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
        << "ORDINARY_MAX_GCD_TESTS="
        << max_ordinary_gcd_tests
        << '\n';

    std::cout
        << "ORDINARY_TOTAL_DISTANCE="
        << ordinary_distance
        << '\n';

    std::cout
        << "ORDINARY_AVG_DISTANCE="
        << avg_ordinary_distance
        << '\n';

    std::cout
        << "STRUCTURED_POOL_TESTS="
        << structured_pool_tests
        << '\n';

    std::cout
        << "STRUCTURED_AVG_POOL_TESTS="
        << avg_structured_pool
        << '\n';

    std::cout
        << "STRUCTURED_ONE_SIDED="
        << structured_one_sided
        << '\n';

    std::cout
        << "STRUCTURED_AVG_ONE_SIDED="
        << avg_one_sided
        << '\n';

    std::cout
        << "STRUCTURED_TWO_SIDED="
        << structured_two_sided
        << '\n';

    std::cout
        << "STRUCTURED_AVG_TWO_SIDED="
        << avg_two_sided
        << '\n';

    std::cout
        << "STRUCTURED_CERTIFICATE_TESTS="
        << structured_certificate_tests
        << '\n';

    std::cout
        << "STRUCTURED_AVG_CERTIFICATE_TESTS="
        << avg_certificate
        << '\n';

    std::cout
        << "STRUCTURED_MAX_ONE_SIDED="
        << max_one_sided
        << '\n';

    std::cout
        << "STRUCTURED_MAX_TWO_SIDED="
        << max_two_sided
        << '\n';

    std::cout
        << "STRUCTURED_MAX_CERTIFICATE_TESTS="
        << max_certificate_tests
        << '\n';

    std::cout
        << "PAIR_PREPROCESSING_TESTS="
        << preprocessing_pair_tests
        << '\n';

    std::cout
        << "PAIR_PREPROCESSING_TIME_MS="
        << preprocessing_time_ms
        << '\n';

    std::cout
        << "AVG_GCD_TEST_REDUCTION_PERCENT="
        << avg_gcd_reduction
        << '\n';

    std::cout
        << "AVG_ONE_SIDED_POOL_REDUCTION_PERCENT="
        << avg_pool_reduction
        << '\n';

    std::cout
        << "AVG_TWO_SIDED_CANDIDATE_REDUCTION_PERCENT="
        << avg_candidate_reduction
        << '\n';

    std::cout
        << "ORDINARY_TIME_MS="
        << ordinary_time_ms
        << '\n';

    std::cout
        << "STRUCTURED_QUERY_TIME_MS="
        << structured_time_ms
        << '\n';

    std::cout
        << "STRUCTURED_TOTAL_TIME_MS="
        << structured_total_time
        << '\n';

    std::cout
        << "STRUCTURED_TOTAL_OVER_ORDINARY="
        << (
            ordinary_time_ms == 0.0
                ? 0.0
                : structured_total_time /
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
