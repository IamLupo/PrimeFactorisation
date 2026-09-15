#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
#include <limits>

using u64 = std::uint64_t;

struct Factor {
    u64 prime;
    int exponent;
};

struct Representation {
    u64 j;
    u64 p;
    u64 n;

    u64 d1;
    u64 d2;

    u64 k1;
    u64 k2;

    u64 gcd_k;
};

struct PrimeRecord {
    u64 count = 0;
    u64 coprime_count = 0;

    u64 min_j = 0;
    u64 max_j = 0;

    u64 first_k1 = 0;
    u64 first_k2 = 0;
};

struct JStats {
    u64 j = 0;

    u64 primes_with_representation = 0;
    u64 total_representations = 0;
    u64 coprime_representations = 0;
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

    for (int n = low; n <= high; ++n) {
        if (prime[n]) {
            result.push_back(n);
        }
    }

    return result;
}

/*
    Smallest-prime-factor sieve.

    This is deliberately used here because we want to inspect
    many values

        n = j*p + 1

    and factor them repeatedly.
*/
std::vector<std::uint32_t> build_spf(int limit) {
    std::vector<std::uint32_t> spf(
        static_cast<std::size_t>(limit) + 1,
        0
    );

    std::vector<int> primes;

    primes.reserve(
        static_cast<std::size_t>(
            limit / 10
        )
    );

    for (int i = 2; i <= limit; ++i) {

        if (spf[i] == 0) {
            spf[i] =
                static_cast<std::uint32_t>(i);

            primes.push_back(i);
        }

        for (int p : primes) {

            if (p > static_cast<int>(spf[i])) {
                break;
            }

            const long long value =
                1LL * p * i;

            if (value > limit) {
                break;
            }

            spf[
                static_cast<std::size_t>(value)
            ] =
                static_cast<std::uint32_t>(p);
        }
    }

    return spf;
}

std::vector<Factor> factorize(
    u64 n,
    const std::vector<std::uint32_t>& spf) {

    std::vector<Factor> factors;

    while (n > 1) {

        const u64 p =
            spf[
                static_cast<std::size_t>(n)
            ];

        int exponent = 0;

        while (n % p == 0) {
            n /= p;
            ++exponent;
        }

        factors.push_back({
            p,
            exponent
        });
    }

    return factors;
}

/*
    Generate all positive divisors of n.

    n <= 25000001 in this experiment, so the divisor count
    stays small.
*/
std::vector<u64> generate_divisors(
    const std::vector<Factor>& factors) {

    std::vector<u64> divisors = {1};

    for (const Factor& factor : factors) {

        const std::size_t old_size =
            divisors.size();

        u64 power = 1;

        for (int e = 1;
             e <= factor.exponent;
             ++e) {

            power *= factor.prime;

            for (std::size_t i = 0;
                 i < old_size;
                 ++i) {

                divisors.push_back(
                    divisors[i] * power
                );
            }
        }
    }

    std::sort(
        divisors.begin(),
        divisors.end()
    );

    return divisors;
}

/*
    For

        n = j*p + 1

    a representation

        p = k1 + k2 + j*k1*k2

    is equivalent to

        n = (j*k1+1)(j*k2+1).

    Therefore we only need a nontrivial divisor d satisfying

        d == 1 (mod j).

    Because n == 1 (mod j), the complementary divisor
    is automatically 1 (mod j) as well.
*/
std::vector<Representation> find_representations(
    u64 p,
    u64 j,
    const std::vector<std::uint32_t>& spf) {

    std::vector<Representation> result;

    const u64 n =
        j * p + 1;

    const auto factors =
        factorize(
            n,
            spf
        );

    const auto divisors =
        generate_divisors(
            factors
        );

    for (u64 d1 : divisors) {

        if (d1 <= 1) {
            continue;
        }

        if (d1 * d1 > n) {
            break;
        }

        if (d1 >= n) {
            continue;
        }

        if (d1 % j != 1) {
            continue;
        }

        const u64 d2 =
            n / d1;

        if (d2 <= d1) {
            continue;
        }

        if (d2 % j != 1) {
            continue;
        }

        const u64 k1 =
            (d1 - 1) / j;

        const u64 k2 =
            (d2 - 1) / j;

        if (k1 == 0 || k2 == 0) {
            continue;
        }

        const u64 gcd_k =
            std::gcd(
                k1,
                k2
            );

        /*
            Explicit reconstruction check.
        */
        const u64 reconstructed =
            k1 +
            k2 +
            j * k1 * k2;

        if (reconstructed != p) {
            continue;
        }

        result.push_back({
            j,
            p,
            n,
            d1,
            d2,
            k1,
            k2,
            gcd_k
        });
    }

    return result;
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
    constexpr int PRIME_LOW = 2;
    constexpr int PRIME_HIGH = 100000;

    /*
        This is intentionally independent of the earlier
        determinant-one m1,m2 restriction.
    */
    constexpr int J_LIMIT = 250;

    const std::vector<u64> J_COVERAGE = {
        1,
        2,
        5,
        10,
        20,
        50,
        100,
        150,
        200,
        250
    };

    const std::vector<int> SELECTED_PRIMES = {
        5,
        7,
        11,
        13,
        17,
        19,
        23,
        29,
        31,
        37,
        41,
        43,
        47,
        97,
        101,
        997,
        99991
    };

    std::cout
        << "START EXPERIMENT 489\n";

    const auto prime_sieve =
        build_sieve(
            PRIME_HIGH
        );

    const auto primes =
        extract_primes(
            PRIME_LOW,
            PRIME_HIGH,
            prime_sieve
        );

    const u64 max_n =
        static_cast<u64>(J_LIMIT) *
        static_cast<u64>(PRIME_HIGH) +
        1;

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    std::cout
        << "J_LIMIT="
        << J_LIMIT
        << "\n";

    std::cout
        << "MAX_JP_PLUS_1="
        << max_n
        << "\n";

    const auto start =
        std::chrono::steady_clock::now();

    /*
        SPF is the main acceleration structure.
    */
    const auto spf =
        build_spf(
            static_cast<int>(max_n)
        );

    std::cout
        << "SPF_READY=1\n";

    /*
        One JStats object for every j.
    */
    std::vector<JStats> j_stats(
        static_cast<std::size_t>(J_LIMIT)
    );

    for (int j = 1; j <= J_LIMIT; ++j) {
        j_stats[
            static_cast<std::size_t>(j - 1)
        ].j = static_cast<u64>(j);
    }

    /*
        Per-prime totals.
    */
    u64 representable_primes = 0;
    u64 nonrepresentable_primes = 0;

    u64 total_representations = 0;
    u64 total_coprime_representations = 0;
    u64 total_non_coprime_representations = 0;

    u64 min_j_total = 0;
    u64 max_j_total = 0;

    u64 max_representation_count = 0;
    u64 max_representation_prime = 0;

    u64 min_representation_count =
        std::numeric_limits<u64>::max();

    u64 min_representation_prime = 0;

    /*
        Cumulative coverage for the selected J limits.
    */
    std::vector<u64> cumulative_coverage(
        J_COVERAGE.size(),
        0
    );

    std::vector<u64> cumulative_coprime_coverage(
        J_COVERAGE.size(),
        0
    );

    /*
        Distribution of smallest j.
    */
    std::vector<u64> min_j_histogram(
        static_cast<std::size_t>(J_LIMIT) + 1,
        0
    );

    /*
        Selected-prime detailed output is recorded after
        the main census as well.
    */
    for (int p : primes) {

        bool prime_has_representation = false;
        bool prime_has_coprime_representation = false;

        u64 prime_representation_count = 0;
        u64 prime_coprime_count = 0;

        u64 prime_min_j =
            std::numeric_limits<u64>::max();

        u64 prime_max_j = 0;

        u64 first_k1 = 0;
        u64 first_k2 = 0;

        for (int j = 1;
             j <= J_LIMIT;
             ++j) {

            const auto representations =
                find_representations(
                    static_cast<u64>(p),
                    static_cast<u64>(j),
                    spf
                );

            if (representations.empty()) {
                continue;
            }

            prime_has_representation = true;

            const u64 j_value =
                static_cast<u64>(j);

            auto& current =
                j_stats[
                    static_cast<std::size_t>(j - 1)
                ];

            ++current.primes_with_representation;

            current.total_representations +=
                representations.size();

            for (const auto& rep :
                 representations) {

                ++prime_representation_count;
                ++total_representations;

                if (rep.gcd_k == 1) {

                    ++prime_coprime_count;
                    ++total_coprime_representations;

                    ++current.coprime_representations;

                    prime_has_coprime_representation =
                        true;
                } else {
                    ++total_non_coprime_representations;
                }

                prime_min_j =
                    std::min(
                        prime_min_j,
                        rep.j
                    );

                prime_max_j =
                    std::max(
                        prime_max_j,
                        rep.j
                    );

                if (first_k1 == 0) {
                    first_k1 = rep.k1;
                    first_k2 = rep.k2;
                }
            }
        }

        if (prime_has_representation) {

            ++representable_primes;

            min_j_total +=
                prime_min_j;

            max_j_total +=
                prime_max_j;

            ++min_j_histogram[
                static_cast<std::size_t>(
                    prime_min_j
                )
            ];

            if (prime_representation_count >
                max_representation_count) {

                max_representation_count =
                    prime_representation_count;

                max_representation_prime =
                    static_cast<u64>(p);
            }

            if (prime_representation_count <
                min_representation_count) {

                min_representation_count =
                    prime_representation_count;

                min_representation_prime =
                    static_cast<u64>(p);
            }
        } else {
            ++nonrepresentable_primes;
        }

        /*
            Cumulative J coverage.
        */
        for (std::size_t i = 0;
             i < J_COVERAGE.size();
             ++i) {

            if (prime_min_j !=
                std::numeric_limits<u64>::max() &&
                prime_min_j <= J_COVERAGE[i]) {

                ++cumulative_coverage[i];
            }

            /*
                A separate pass for coprime
                representation coverage.
            */
            bool coprime_within_limit = false;

            for (int jj = 1;
                 jj <=
                    static_cast<int>(
                        J_COVERAGE[i]
                    );
                 ++jj) {

                const auto reps =
                    find_representations(
                        static_cast<u64>(p),
                        static_cast<u64>(jj),
                        spf
                    );

                for (const auto& rep : reps) {
                    if (rep.gcd_k == 1) {
                        coprime_within_limit = true;
                        break;
                    }
                }

                if (coprime_within_limit) {
                    break;
                }
            }

            if (coprime_within_limit) {
                ++cumulative_coprime_coverage[i];
            }
        }
    }

    /*
        Basic census.
    */
    std::cout
        << "TOTAL_PRIMES="
        << primes.size()
        << "\n";

    std::cout
        << "REPRESENTABLE_PRIMES="
        << representable_primes
        << "\n";

    std::cout
        << "NONREPRESENTABLE_PRIMES="
        << nonrepresentable_primes
        << "\n";

    std::cout
        << "COVERAGE_PERCENT="
        << 100.0 *
           static_cast<double>(
               representable_primes
           ) /
           static_cast<double>(
               primes.size()
           )
        << "\n";

    std::cout
        << "TOTAL_REPRESENTATIONS="
        << total_representations
        << "\n";

    std::cout
        << "TOTAL_COPRIME_REPRESENTATIONS="
        << total_coprime_representations
        << "\n";

    std::cout
        << "TOTAL_NONCOPRIME_REPRESENTATIONS="
        << total_non_coprime_representations
        << "\n";

    if (representable_primes > 0) {

        std::cout
            << "AVG_REPRESENTATIONS_PER_REPRESENTABLE_PRIME="
            << average(
                   total_representations,
                   representable_primes
               )
            << "\n";

        std::cout
            << "AVG_MIN_J="
            << average(
                   min_j_total,
                   representable_primes
               )
            << "\n";

        std::cout
            << "AVG_MAX_J="
            << average(
                   max_j_total,
                   representable_primes
               )
            << "\n";
    } else {
        std::cout
            << "AVG_REPRESENTATIONS_PER_REPRESENTABLE_PRIME=0\n";

        std::cout
            << "AVG_MIN_J=0\n";

        std::cout
            << "AVG_MAX_J=0\n";
    }

    std::cout
        << "MAX_REPRESENTATION_COUNT="
        << max_representation_count
        << "\n";

    std::cout
        << "MAX_REPRESENTATION_PRIME="
        << max_representation_prime
        << "\n";

    std::cout
        << "MIN_REPRESENTATION_COUNT="
        << min_representation_count
        << "\n";

    std::cout
        << "MIN_REPRESENTATION_PRIME="
        << min_representation_prime
        << "\n";

    /*
        Coverage as J grows.
    */
    for (std::size_t i = 0;
         i < J_COVERAGE.size();
         ++i) {

        std::cout
            << "J_LIMIT="
            << J_COVERAGE[i]
            << " COVERAGE="
            << cumulative_coverage[i]
            << " COVERAGE_PERCENT="
            << 100.0 *
               static_cast<double>(
                   cumulative_coverage[i]
               ) /
               static_cast<double>(
                   primes.size()
               )
            << " COPRIME_COVERAGE="
            << cumulative_coprime_coverage[i]
            << " COPRIME_COVERAGE_PERCENT="
            << 100.0 *
               static_cast<double>(
                   cumulative_coprime_coverage[i]
               ) /
               static_cast<double>(
                   primes.size()
               )
            << "\n";
    }

    /*
        Per-j statistics.
    */
    for (const auto& stat : j_stats) {

        if (stat.total_representations == 0) {
            continue;
        }

        std::cout
            << "J="
            << stat.j
            << " PRIMES_WITH_REPRESENTATION="
            << stat.primes_with_representation
            << " TOTAL_REPRESENTATIONS="
            << stat.total_representations
            << " COPRIME_REPRESENTATIONS="
            << stat.coprime_representations
            << "\n";
    }

    /*
        Distribution of minimum j.
    */
    for (std::size_t j = 1;
         j < min_j_histogram.size();
         ++j) {

        if (min_j_histogram[j] == 0) {
            continue;
        }

        std::cout
            << "MIN_J="
            << j
            << " PRIME_COUNT="
            << min_j_histogram[j]
            << "\n";
    }

    /*
        Detailed selected primes.
    */
    for (int p : SELECTED_PRIMES) {

        u64 total = 0;
        u64 coprime = 0;

        u64 min_j =
            std::numeric_limits<u64>::max();

        Representation first_rep{};
        bool first_found = false;

        for (int j = 1;
             j <= J_LIMIT;
             ++j) {

            const auto reps =
                find_representations(
                    static_cast<u64>(p),
                    static_cast<u64>(j),
                    spf
                );

            for (const auto& rep :
                 reps) {

                ++total;

                if (rep.gcd_k == 1) {
                    ++coprime;
                }

                if (rep.j < min_j) {
                    min_j = rep.j;

                    first_rep = rep;
                    first_found = true;
                }
            }
        }

        if (!first_found) {

            std::cout
                << "SELECTED_PRIME="
                << p
                << " REPRESENTATIONS=0\n";

            continue;
        }

        std::cout
            << "SELECTED_PRIME="
            << p
            << " REPRESENTATIONS="
            << total
            << " COPRIME="
            << coprime
            << " MIN_J="
            << min_j
            << " K1="
            << first_rep.k1
            << " K2="
            << first_rep.k2
            << " D1="
            << first_rep.d1
            << " D2="
            << first_rep.d2
            << " JP_PLUS_1="
            << first_rep.n
            << " RECONSTRUCTED_P="
            << first_rep.k1 +
                   first_rep.k2 +
                   first_rep.j *
                       first_rep.k1 *
                       first_rep.k2
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
        << "FINISHED EXPERIMENT 489\n";

    return 0;
}
