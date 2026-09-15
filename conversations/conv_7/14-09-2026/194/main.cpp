#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
#include <limits>

using u64 = std::uint64_t;

struct Pair {
    int k1;
    int k2;
    u64 a;
    u64 M;
};

struct PrimeRecord {
    bool found = false;

    u64 representation_count = 0;

    u64 min_M = 0;
    u64 min_J = 0;

    int min_k1 = 0;
    int min_k2 = 0;

    int min_max_k = 0;

    u64 max_M = 0;
    u64 max_J = 0;
};

struct CoverageStat {
    int limit;
    u64 covered = 0;
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

std::vector<int> extract_primes(
    int low,
    int high,
    const std::vector<bool>& prime) {

    std::vector<int> result;

    for (int p = low; p <= high; ++p) {
        if (prime[p]) {
            result.push_back(p);
        }
    }

    return result;
}

/*
    Generate all unique (k1,k2) pairs satisfying

        m2*k1 - m1*k2 = 1

    with

        2 <= k1 <= k2 <= k_limit
        1 <= m1,m2 <= m_limit
        gcd(k1,k2)=1.
*/
std::vector<Pair> build_pairs(
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

                if (k2 < k1 ||
                    k2 > k_limit ||
                    k2 < 2) {
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
            if (lhs.k1 != rhs.k1) {
                return lhs.k1 < rhs.k1;
            }

            if (lhs.k2 != rhs.k2) {
                return lhs.k2 < rhs.k2;
            }

            if (lhs.M != rhs.M) {
                return lhs.M < rhs.M;
            }

            return lhs.a < rhs.a;
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

/*
    Examine one prime against all determinant-one pairs.

    Representation condition:

        p = a + j*M

    with

        a = k1+k2
        M = k1*k2

    and j >= 0.

    Equivalent congruence:

        p mod M = a mod M.
*/
PrimeRecord analyze_prime(
    int p,
    const std::vector<Pair>& pairs) {

    PrimeRecord record;

    for (const Pair& pair : pairs) {

        if (pair.a > static_cast<u64>(p)) {
            continue;
        }

        const u64 pp =
            static_cast<u64>(p);

        if (pp % pair.M != pair.a % pair.M) {
            continue;
        }

        const u64 j =
            (pp - pair.a) / pair.M;

        ++record.representation_count;

        if (!record.found) {
            record.found = true;

            record.min_M = pair.M;
            record.min_J = j;

            record.min_k1 = pair.k1;
            record.min_k2 = pair.k2;

            record.min_max_k =
                std::max(
                    pair.k1,
                    pair.k2
                );

            record.max_M = pair.M;
            record.max_J = j;

            continue;
        }

        if (pair.M < record.min_M ||
            (pair.M == record.min_M &&
             j < record.min_J)) {

            record.min_M = pair.M;
            record.min_J = j;

            record.min_k1 = pair.k1;
            record.min_k2 = pair.k2;

            record.min_max_k =
                std::max(
                    pair.k1,
                    pair.k2
                );
        }

        if (pair.M > record.max_M) {
            record.max_M = pair.M;
        }

        if (j > record.max_J) {
            record.max_J = j;
        }
    }

    return record;
}

double average_u64(
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
    constexpr int PRIME_LOW = 2;
    constexpr int PRIME_HIGH = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_LIMIT = 7;

    /*
        Coverage as k-space expands.
    */
    const std::vector<int> K_CUTS = {
        2,
        3,
        5,
        10,
        20,
        50,
        100,
        250,
        500,
        1000
    };

    /*
        Coverage as the determinant coefficients expand.
    */
    const std::vector<int> M_CUTS = {
        1,
        2,
        3,
        4,
        5,
        6,
        7
    };

    std::cout
        << "START EXPERIMENT 487\n";

    const auto sieve =
        build_sieve(PRIME_HIGH);

    const auto primes =
        extract_primes(
            PRIME_LOW,
            PRIME_HIGH,
            sieve
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << "\n";

    std::cout
        << "M_LIMIT="
        << M_LIMIT
        << "\n";

    const auto start =
        std::chrono::steady_clock::now();

    /*
        Full representation family.
    */
    const auto pairs =
        build_pairs(
            K_LIMIT,
            M_LIMIT
        );

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << "\n";

    /*
        Full census.
    */
    u64 total_found = 0;
    u64 total_representations = 0;

    u64 min_M_total = 0;
    u64 max_M_total = 0;

    u64 min_J_total = 0;
    u64 max_J_total = 0;

    u64 min_k1_total = 0;
    u64 min_k2_total = 0;

    u64 max_representation_count = 0;

    u64 max_representation_prime = 0;

    u64 min_representation_count =
        std::numeric_limits<u64>::max();

    u64 min_representation_prime = 0;

    int global_max_min_M = 0;
    u64 global_max_min_M_prime = 0;

    int global_max_min_J = 0;
    u64 global_max_min_J_prime = 0;

    u64 representation_histogram[11] = {};

    PrimeRecord largest_prime_record;

    /*
        Keep the records for the first few primes and a few
        selected primes near the top of the range.
    */
    std::vector<int> selected_primes;

    for (int p : primes) {
        if (p <= 50 ||
            p >= 99900 ||
            p == 99991 ||
            p == 99989) {

            selected_primes.push_back(p);
        }
    }

    for (int p : primes) {

        const PrimeRecord record =
            analyze_prime(
                p,
                pairs
            );

        if (record.found) {
            ++total_found;

            total_representations +=
                record.representation_count;

            min_M_total +=
                record.min_M;

            max_M_total +=
                record.max_M;

            min_J_total +=
                record.min_J;

            max_J_total +=
                record.max_J;

            min_k1_total +=
                static_cast<u64>(
                    record.min_k1
                );

            min_k2_total +=
                static_cast<u64>(
                    record.min_k2
                );

            if (record.representation_count >
                max_representation_count) {

                max_representation_count =
                    record.representation_count;

                max_representation_prime =
                    static_cast<u64>(p);
            }

            if (record.representation_count <
                min_representation_count) {

                min_representation_count =
                    record.representation_count;

                min_representation_prime =
                    static_cast<u64>(p);
            }

            if (record.min_M >
                static_cast<u64>(
                    global_max_min_M
                )) {

                global_max_min_M =
                    static_cast<int>(
                        record.min_M
                    );

                global_max_min_M_prime =
                    static_cast<u64>(p);
            }

            if (record.min_J >
                static_cast<u64>(
                    global_max_min_J
                )) {

                global_max_min_J =
                    static_cast<int>(
                        record.min_J
                    );

                global_max_min_J_prime =
                    static_cast<u64>(p);
            }

            const u64 bucket =
                std::min<u64>(
                    10,
                    record.representation_count
                );

            ++representation_histogram[bucket];

            if (p >= 99900) {
                largest_prime_record =
                    record;
            }
        }
    }

    const u64 total_primes =
        static_cast<u64>(
            primes.size()
        );

    /*
        Print overall census.
    */
    std::cout
        << "TOTAL_PRIMES="
        << total_primes
        << "\n";

    std::cout
        << "REPRESENTABLE_PRIMES="
        << total_found
        << "\n";

    std::cout
        << "NONREPRESENTABLE_PRIMES="
        << total_primes - total_found
        << "\n";

    std::cout
        << "COVERAGE_PERCENT="
        << 100.0 *
           static_cast<double>(
               total_found
           ) /
           static_cast<double>(
               total_primes
           )
        << "\n";

    std::cout
        << "TOTAL_REPRESENTATIONS="
        << total_representations
        << "\n";

    if (total_found > 0) {

        std::cout
            << "AVG_REPRESENTATIONS_PER_PRIME="
            << average_u64(
                   total_representations,
                   total_found
               )
            << "\n";

        std::cout
            << "AVG_MIN_M="
            << average_u64(
                   min_M_total,
                   total_found
               )
            << "\n";

        std::cout
            << "AVG_MIN_J="
            << average_u64(
                   min_J_total,
                   total_found
               )
            << "\n";

        std::cout
            << "AVG_MIN_K1="
            << average_u64(
                   min_k1_total,
                   total_found
               )
            << "\n";

        std::cout
            << "AVG_MIN_K2="
            << average_u64(
                   min_k2_total,
                   total_found
               )
            << "\n";

        std::cout
            << "MAX_OF_MIN_M="
            << global_max_min_M
            << "\n";

        std::cout
            << "PRIME_WITH_MAX_OF_MIN_M="
            << global_max_min_M_prime
            << "\n";

        std::cout
            << "MAX_OF_MIN_J="
            << global_max_min_J
            << "\n";

        std::cout
            << "PRIME_WITH_MAX_OF_MIN_J="
            << global_max_min_J_prime
            << "\n";

        std::cout
            << "MAX_REPRESENTATION_COUNT="
            << max_representation_count
            << "\n";

        std::cout
            << "PRIME_WITH_MAX_REPRESENTATION_COUNT="
            << max_representation_prime
            << "\n";

        std::cout
            << "MIN_REPRESENTATION_COUNT="
            << min_representation_count
            << "\n";

        std::cout
            << "PRIME_WITH_MIN_REPRESENTATION_COUNT="
            << min_representation_prime
            << "\n";
    }

    /*
        Histogram of representation multiplicity.

        Bucket 10 means "10 or more".
    */
    for (int bucket = 1;
         bucket <= 10;
         ++bucket) {

        std::cout
            << "REPRESENTATION_COUNT_BUCKET="
            << bucket
            << " COUNT="
            << representation_histogram[bucket]
            << "\n";
    }

    /*
        Representation counts grouped by minimum M.
    */
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

    for (u64 limit : M_LIMITS) {

        u64 covered = 0;

        for (int p : primes) {

            bool found = false;

            for (const Pair& pair : pairs) {

                if (pair.M > limit) {
                    break;
                }

                if (pair.a > static_cast<u64>(p)) {
                    continue;
                }

                const u64 pp =
                    static_cast<u64>(p);

                if (pp % pair.M ==
                    pair.a % pair.M) {

                    found = true;
                    break;
                }
            }

            if (found) {
                ++covered;
            }
        }

        std::cout
            << "M_LIMIT="
            << limit
            << " COVERAGE="
            << covered
            << " COVERAGE_PERCENT="
            << 100.0 *
               static_cast<double>(
                   covered
               ) /
               static_cast<double>(
                   total_primes
               )
            << "\n";
    }

    /*
        Coverage grouped by k-limit.

        Rebuild the pair family for each cutoff.
    */
    for (int k_cut : K_CUTS) {

        const auto cut_pairs =
            build_pairs(
                k_cut,
                M_LIMIT
            );

        u64 covered = 0;

        u64 representation_total = 0;

        for (int p : primes) {

            u64 count = 0;

            for (const Pair& pair : cut_pairs) {

                if (pair.a >
                    static_cast<u64>(p)) {

                    continue;
                }

                const u64 pp =
                    static_cast<u64>(p);

                if (pp % pair.M ==
                    pair.a % pair.M) {

                    ++count;
                }
            }

            if (count > 0) {
                ++covered;
            }

            representation_total += count;
        }

        std::cout
            << "K_LIMIT_CUT="
            << k_cut
            << " PAIR_COUNT="
            << cut_pairs.size()
            << " COVERAGE="
            << covered
            << " COVERAGE_PERCENT="
            << 100.0 *
               static_cast<double>(
                   covered
               ) /
               static_cast<double>(
                   total_primes
               )
            << " TOTAL_REPRESENTATIONS="
            << representation_total
            << "\n";
    }

    /*
        Selected-prime detailed records.
    */
    for (int p : selected_primes) {

        const PrimeRecord record =
            analyze_prime(
                p,
                pairs
            );

        std::cout
            << "PRIME="
            << p
            << " FOUND="
            << (record.found ? 1 : 0)
            << " REPRESENTATIONS="
            << record.representation_count
            << " MIN_M="
            << record.min_M
            << " MIN_J="
            << record.min_J
            << " MIN_K1="
            << record.min_k1
            << " MIN_K2="
            << record.min_k2
            << " MAX_M="
            << record.max_M
            << " MAX_J="
            << record.max_J
            << "\n";
    }

    /*
        Final selected high-prime record.
    */
    if (!selected_primes.empty()) {

        /*
            This section is intentionally only a summary,
            not an additional search.
        */
        std::cout
            << "LAST_SELECTED_PRIME="
            << selected_primes.back()
            << "\n";
    }

    const auto end =
        std::chrono::steady_clock::now();

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 487\n";

    return 0;
}
