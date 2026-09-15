#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

using u64 = std::uint64_t;

struct Factor {
    u64 prime;
    int exponent;
};

struct JStats {
    u64 j = 0;

    u64 total = 0;
    u64 composite = 0;
    u64 prime = 0;

    u64 qualifying = 0;
    u64 composite_without_qualifying = 0;

    u64 total_nontrivial_factor_pairs = 0;
    u64 total_qualifying_factor_pairs = 0;

    u64 max_qualifying_factor_pairs = 0;
    u64 prime_with_max_qualifying_pairs = 0;

    u64 first_composite_without_qualifying = 0;
    u64 first_qualifying_number = 0;
};

struct Example {
    bool found = false;

    u64 p = 0;
    u64 j = 0;
    u64 n = 0;

    u64 d1 = 0;
    u64 d2 = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 factor_pair_count = 0;
    u64 qualifying_pair_count = 0;
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

std::vector<std::uint32_t> build_spf(int limit) {
    std::vector<std::uint32_t> spf(
        static_cast<std::size_t>(limit) + 1,
        0
    );

    std::vector<int> primes;

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
    Analyze a single number

        n = j*p + 1.

    Since n == 1 (mod j), if a divisor d satisfies

        d == 1 (mod j),

    then its complementary divisor n/d automatically also
    satisfies

        n/d == 1 (mod j).

    Therefore only one side of each factor pair needs to
    be checked.
*/
void analyze_number(
    u64 p,
    u64 j,
    const std::vector<std::uint32_t>& spf,
    JStats& stats,
    Example& first_missing,
    Example& first_qualifying) {

    const u64 n =
        j * p + 1;

    ++stats.total;

    if (n < 4) {
        ++stats.prime;
        return;
    }

    const auto factors =
        factorize(
            n,
            spf
        );

    const bool is_prime_number =
        factors.size() == 1 &&
        factors[0].exponent == 1;

    if (is_prime_number) {
        ++stats.prime;
        return;
    }

    ++stats.composite;

    const auto divisors =
        generate_divisors(
            factors
        );

    u64 factor_pairs = 0;
    u64 qualifying_pairs = 0;

    u64 first_d1 = 0;
    u64 first_d2 = 0;
    u64 first_k1 = 0;
    u64 first_k2 = 0;

    for (u64 d1 : divisors) {

        if (d1 <= 1) {
            continue;
        }

        if (d1 * d1 >= n) {
            break;
        }

        if (n % d1 != 0) {
            continue;
        }

        ++factor_pairs;

        if (d1 % j != 1) {
            continue;
        }

        const u64 d2 =
            n / d1;

        if (d2 % j != 1) {
            continue;
        }

        const u64 k1 =
            (d1 - 1) / j;

        const u64 k2 =
            (d2 - 1) / j;

        /*
            Explicit identity verification:

                p = k1+k2+j*k1*k2
        */
        const u64 reconstructed =
            k1 +
            k2 +
            j * k1 * k2;

        if (reconstructed != p) {
            continue;
        }

        ++qualifying_pairs;

        if (first_d1 == 0) {
            first_d1 = d1;
            first_d2 = d2;
            first_k1 = k1;
            first_k2 = k2;
        }
    }

    stats.total_nontrivial_factor_pairs +=
        factor_pairs;

    stats.total_qualifying_factor_pairs +=
        qualifying_pairs;

    stats.max_qualifying_factor_pairs =
        std::max(
            stats.max_qualifying_factor_pairs,
            qualifying_pairs
        );

    if (qualifying_pairs >
        stats.max_qualifying_factor_pairs) {

        stats.prime_with_max_qualifying_pairs =
            p;
    }

    if (qualifying_pairs > 0) {

        ++stats.qualifying;

        if (stats.first_qualifying_number == 0) {
            stats.first_qualifying_number = n;
        }

        if (!first_qualifying.found) {

            first_qualifying.found = true;
            first_qualifying.p = p;
            first_qualifying.j = j;
            first_qualifying.n = n;
            first_qualifying.d1 = first_d1;
            first_qualifying.d2 = first_d2;
            first_qualifying.k1 = first_k1;
            first_qualifying.k2 = first_k2;
            first_qualifying.factor_pair_count =
                factor_pairs;
            first_qualifying.qualifying_pair_count =
                qualifying_pairs;
        }

    } else {

        ++stats.composite_without_qualifying;

        if (!first_missing.found) {

            first_missing.found = true;
            first_missing.p = p;
            first_missing.j = j;
            first_missing.n = n;
            first_missing.factor_pair_count =
                factor_pairs;
            first_missing.qualifying_pair_count =
                qualifying_pairs;
        }
    }
}

double percentage(
    u64 numerator,
    u64 denominator) {

    if (denominator == 0) {
        return 0.0;
    }

    return
        100.0 *
        static_cast<double>(numerator) /
        static_cast<double>(denominator);
}

double average(
    u64 numerator,
    u64 denominator) {

    if (denominator == 0) {
        return 0.0;
    }

    return
        static_cast<double>(numerator) /
        static_cast<double>(denominator);
}

int main() {
    constexpr int PRIME_LOW = 2;
    constexpr int PRIME_HIGH = 100000;

    constexpr int J_LIMIT = 100;

    const std::vector<int> SELECTED_J = {
        2,
        3,
        4,
        5,
        6,
        7,
        10,
        20,
        50,
        100
    };

    std::cout
        << "START EXPERIMENT 491\n";

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

    /*
        Maximum number we factor:

            J_LIMIT * PRIME_HIGH + 1.
    */
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

    const auto spf =
        build_spf(
            static_cast<int>(max_n)
        );

    std::cout
        << "SPF_READY=1\n";

    std::vector<JStats> stats(
        static_cast<std::size_t>(J_LIMIT) + 1
    );

    for (int j = 1;
         j <= J_LIMIT;
         ++j) {

        stats[
            static_cast<std::size_t>(j)
        ].j = static_cast<u64>(j);
    }

    std::vector<Example> first_missing(
        static_cast<std::size_t>(J_LIMIT) + 1
    );

    std::vector<Example> first_qualifying(
        static_cast<std::size_t>(J_LIMIT) + 1
    );

    /*
        Global totals.
    */
    u64 total_composite = 0;
    u64 total_prime = 0;
    u64 total_qualifying = 0;
    u64 total_composite_without_qualifying = 0;

    u64 total_factor_pairs = 0;
    u64 total_qualifying_pairs = 0;

    for (int j = 1;
         j <= J_LIMIT;
         ++j) {

        JStats& current =
            stats[
                static_cast<std::size_t>(j)
            ];

        for (int p : primes) {

            analyze_number(
                static_cast<u64>(p),
                static_cast<u64>(j),
                spf,
                current,
                first_missing[
                    static_cast<std::size_t>(j)
                ],
                first_qualifying[
                    static_cast<std::size_t>(j)
                ]
            );
        }

        total_composite += current.composite;
        total_prime += current.prime;
        total_qualifying += current.qualifying;
        total_composite_without_qualifying +=
            current.composite_without_qualifying;

        total_factor_pairs +=
            current.total_nontrivial_factor_pairs;

        total_qualifying_pairs +=
            current.total_qualifying_factor_pairs;
    }

    /*
        Selected-j summary.
    */
    for (int j : SELECTED_J) {

        const JStats& current =
            stats[
                static_cast<std::size_t>(j)
            ];

        std::cout
            << "J="
            << j
            << " TOTAL="
            << current.total
            << " COMPOSITE="
            << current.composite
            << " PRIME="
            << current.prime
            << " QUALIFYING="
            << current.qualifying
            << " COMPOSITE_WITHOUT_QUALIFYING="
            << current.composite_without_qualifying
            << " COMPOSITE_TO_QUALIFYING_PERCENT="
            << percentage(
                   current.qualifying,
                   current.composite
               )
            << " TOTAL_FACTOR_PAIRS="
            << current.total_nontrivial_factor_pairs
            << " QUALIFYING_FACTOR_PAIRS="
            << current.total_qualifying_factor_pairs
            << "\n";
    }

    /*
        Full j table.
    */
    for (int j = 1;
         j <= J_LIMIT;
         ++j) {

        const JStats& current =
            stats[
                static_cast<std::size_t>(j)
            ];

        std::cout
            << "FULL_J="
            << j
            << " COMPOSITE="
            << current.composite
            << " QUALIFYING="
            << current.qualifying
            << " WITHOUT="
            << current.composite_without_qualifying
            << " QUALIFY_PERCENT="
            << percentage(
                   current.qualifying,
                   current.composite
               )
            << "\n";
    }

    /*
        First composite counterexample for each selected j.
    */
    for (int j : SELECTED_J) {

        const Example& ex =
            first_missing[
                static_cast<std::size_t>(j)
            ];

        if (!ex.found) {

            std::cout
                << "FIRST_MISSING_J="
                << j
                << " FOUND=0\n";

            continue;
        }

        std::cout
            << "FIRST_MISSING_J="
            << j
            << " FOUND=1"
            << " P="
            << ex.p
            << " N="
            << ex.n
            << " FACTOR_PAIRS="
            << ex.factor_pair_count
            << " QUALIFYING_PAIRS="
            << ex.qualifying_pair_count
            << "\n";
    }

    /*
        First successful example for each selected j.
    */
    for (int j : SELECTED_J) {

        const Example& ex =
            first_qualifying[
                static_cast<std::size_t>(j)
            ];

        if (!ex.found) {

            std::cout
                << "FIRST_QUALIFYING_J="
                << j
                << " FOUND=0\n";

            continue;
        }

        std::cout
            << "FIRST_QUALIFYING_J="
            << j
            << " FOUND=1"
            << " P="
            << ex.p
            << " N="
            << ex.n
            << " D1="
            << ex.d1
            << " D2="
            << ex.d2
            << " K1="
            << ex.k1
            << " K2="
            << ex.k2
            << " FACTOR_PAIRS="
            << ex.factor_pair_count
            << " QUALIFYING_PAIRS="
            << ex.qualifying_pair_count
            << "\n";
    }

    /*
        Global summary.
    */
    std::cout
        << "TOTAL_COMPOSITE="
        << total_composite
        << "\n";

    std::cout
        << "TOTAL_PRIME="
        << total_prime
        << "\n";

    std::cout
        << "TOTAL_QUALIFYING="
        << total_qualifying
        << "\n";

    std::cout
        << "TOTAL_COMPOSITE_WITHOUT_QUALIFYING="
        << total_composite_without_qualifying
        << "\n";

    std::cout
        << "TOTAL_FACTOR_PAIRS="
        << total_factor_pairs
        << "\n";

    std::cout
        << "TOTAL_QUALIFYING_FACTOR_PAIRS="
        << total_qualifying_pairs
        << "\n";

    std::cout
        << "GLOBAL_COMPOSITE_TO_QUALIFYING_PERCENT="
        << percentage(
               total_qualifying,
               total_composite
           )
        << "\n";

    std::cout
        << "GLOBAL_COMPOSITE_WITHOUT_QUALIFYING_PERCENT="
        << percentage(
               total_composite_without_qualifying,
               total_composite
           )
        << "\n";

    std::cout
        << "AVG_FACTOR_PAIRS_PER_COMPOSITE="
        << average(
               total_factor_pairs,
               total_composite
           )
        << "\n";

    std::cout
        << "AVG_QUALIFYING_PAIRS_PER_COMPOSITE="
        << average(
               total_qualifying_pairs,
               total_composite
           )
        << "\n";

    /*
        Important special case j=2.
    */
    const JStats& j2 =
        stats[2];

    std::cout
        << "J2_COMPOSITE="
        << j2.composite
        << "\n";

    std::cout
        << "J2_QUALIFYING="
        << j2.qualifying
        << "\n";

    std::cout
        << "J2_COMPOSITE_WITHOUT_QUALIFYING="
        << j2.composite_without_qualifying
        << "\n";

    std::cout
        << "J2_COMPOSITE_TO_QUALIFYING_PERCENT="
        << percentage(
               j2.qualifying,
               j2.composite
           )
        << "\n";

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
        << "FINISHED EXPERIMENT 491\n";

    return 0;
}
