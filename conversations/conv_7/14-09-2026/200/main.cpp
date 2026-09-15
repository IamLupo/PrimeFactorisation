#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using u64 = std::uint64_t;

struct Factor {
    u64 prime;
    int exponent;
};

struct JStats {
    u64 j = 0;

    u64 composite = 0;
    u64 qualifying = 0;
    u64 failures = 0;

    u64 total_factor_pairs = 0;
    u64 total_qualifying_pairs = 0;

    u64 max_factor_pairs = 0;
    u64 max_qualifying_pairs = 0;

    u64 total_omega = 0;
    u64 max_omega = 0;

    u64 first_failure_p = 0;
    u64 first_failure_n = 0;
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
    Direct mathematical criterion.

    For

        n = j*p + 1,

    a representation

        p = k1 + k2 + j*k1*k2

    exists exactly when n has a proper factor d satisfying

        d == 1 (mod j).

    If

        d = j*k1 + 1

    and

        n/d = j*k2 + 1,

    then the representation follows immediately.
*/
u64 count_qualifying_factor_pairs(
    u64 p,
    u64 j,
    const std::vector<std::uint32_t>& spf) {

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

    u64 count = 0;

    for (u64 d : divisors) {

        /*
            We only count one representative from each
            pair d, n/d by restricting to d < sqrt(n).
        */
        if (d <= 1) {
            continue;
        }

        if (d * d >= n) {
            break;
        }

        if (n % d != 0) {
            continue;
        }

        if (d % j != 1) {
            continue;
        }

        const u64 complementary =
            n / d;

        if (complementary % j != 1) {
            continue;
        }

        /*
            Explicitly reconstruct k1,k2 and p.
        */
        const u64 k1 =
            (d - 1) / j;

        const u64 k2 =
            (complementary - 1) / j;

        const u64 reconstructed =
            k1 +
            k2 +
            j * k1 * k2;

        if (reconstructed != p) {
            continue;
        }

        ++count;
    }

    return count;
}

u64 count_factor_pairs(
    u64 n,
    const std::vector<std::uint32_t>& spf) {

    const auto factors =
        factorize(
            n,
            spf
        );

    const auto divisors =
        generate_divisors(
            factors
        );

    u64 count = 0;

    for (u64 d : divisors) {

        if (d <= 1) {
            continue;
        }

        if (d * d >= n) {
            break;
        }

        if (n % d == 0) {
            ++count;
        }
    }

    return count;
}

int total_prime_factor_multiplicity(
    const std::vector<Factor>& factors) {

    int result = 0;

    for (const auto& factor : factors) {
        result += factor.exponent;
    }

    return result;
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
        8,
        10,
        20,
        50,
        100
    };

    std::cout
        << "START EXPERIMENT 494\n";

    const auto prime_sieve =
        build_sieve(
            PRIME_HIGH
        );

    const auto primes =
        generate_primes(
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

    const auto spf =
        build_spf(
            static_cast<int>(max_n)
        );

    std::cout
        << "SPF_READY=1\n";

    std::vector<JStats> stats(
        static_cast<std::size_t>(
            J_LIMIT + 1
        )
    );

    for (int j = 2;
         j <= J_LIMIT;
         ++j) {

        JStats& current =
            stats[
                static_cast<std::size_t>(j)
            ];

        current.j =
            static_cast<u64>(j);

        for (int p : primes) {

            const u64 n =
                static_cast<u64>(j) *
                static_cast<u64>(p) +
                1;

            const auto factors =
                factorize(
                    n,
                    spf
                );

            const bool is_prime_number =
                factors.size() == 1 &&
                factors[0].exponent == 1;

            if (is_prime_number) {
                continue;
            }

            ++current.composite;

            const int omega =
                total_prime_factor_multiplicity(
                    factors
                );

            current.total_omega +=
                static_cast<u64>(omega);

            current.max_omega =
                std::max(
                    current.max_omega,
                    static_cast<u64>(omega)
                );

            const u64 factor_pairs =
                count_factor_pairs(
                    n,
                    spf
                );

            const u64 qualifying_pairs =
                count_qualifying_factor_pairs(
                    static_cast<u64>(p),
                    static_cast<u64>(j),
                    spf
                );

            current.total_factor_pairs +=
                factor_pairs;

            current.total_qualifying_pairs +=
                qualifying_pairs;

            current.max_factor_pairs =
                std::max(
                    current.max_factor_pairs,
                    factor_pairs
                );

            current.max_qualifying_pairs =
                std::max(
                    current.max_qualifying_pairs,
                    qualifying_pairs
                );

            if (qualifying_pairs > 0) {

                ++current.qualifying;

            } else {

                ++current.failures;

                if (current.first_failure_p == 0) {
                    current.first_failure_p =
                        static_cast<u64>(p);

                    current.first_failure_n =
                        n;
                }
            }
        }
    }

    /*
        Sanity check explicitly required by parity.

        For j=2, every composite odd n has an odd proper
        divisor, so every composite 2p+1 must qualify.
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
        << "J2_FAILURES="
        << j2.failures
        << "\n";

    std::cout
        << "J2_SANITY_PASS="
        << (
            j2.composite == j2.qualifying
                ? 1
                : 0
        )
        << "\n";

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
            << " COMPOSITE="
            << current.composite
            << " QUALIFYING="
            << current.qualifying
            << " FAILURES="
            << current.failures
            << " QUALIFY_PERCENT="
            << percentage(
                   current.qualifying,
                   current.composite
               )
            << " AVG_FACTOR_PAIRS="
            << average(
                   current.total_factor_pairs,
                   current.composite
               )
            << " AVG_QUALIFYING_PAIRS="
            << average(
                   current.total_qualifying_pairs,
                   current.composite
               )
            << " AVG_OMEGA="
            << average(
                   current.total_omega,
                   current.composite
               )
            << " MAX_OMEGA="
            << current.max_omega
            << "\n";
    }

    /*
        Full j table.
    */
    for (int j = 2;
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
            << " FAILURES="
            << current.failures
            << " QUALIFY_PERCENT="
            << percentage(
                   current.qualifying,
                   current.composite
               )
            << "\n";
    }

    /*
        First failure examples.
    */
    for (int j : SELECTED_J) {

        const JStats& current =
            stats[
                static_cast<std::size_t>(j)
            ];

        if (current.first_failure_p == 0) {

            std::cout
                << "FIRST_FAILURE_J="
                << j
                << " FOUND=0\n";

            continue;
        }

        const u64 p =
            current.first_failure_p;

        const u64 n =
            current.first_failure_n;

        const auto factors =
            factorize(
                n,
                spf
            );

        std::cout
            << "FIRST_FAILURE_J="
            << j
            << " P="
            << p
            << " N="
            << n
            << " FACTORS=";

        for (std::size_t i = 0;
             i < factors.size();
             ++i) {

            if (i != 0) {
                std::cout << "*";
            }

            std::cout
                << factors[i].prime;

            if (factors[i].exponent > 1) {
                std::cout
                    << "^"
                    << factors[i].exponent;
            }
        }

        std::cout
            << "\n";
    }

    /*
        Global totals.
    */
    u64 total_composite = 0;
    u64 total_qualifying = 0;
    u64 total_failures = 0;

    for (int j = 2;
         j <= J_LIMIT;
         ++j) {

        const JStats& current =
            stats[
                static_cast<std::size_t>(j)
            ];

        total_composite +=
            current.composite;

        total_qualifying +=
            current.qualifying;

        total_failures +=
            current.failures;
    }

    std::cout
        << "TOTAL_COMPOSITE="
        << total_composite
        << "\n";

    std::cout
        << "TOTAL_QUALIFYING="
        << total_qualifying
        << "\n";

    std::cout
        << "TOTAL_FAILURES="
        << total_failures
        << "\n";

    std::cout
        << "GLOBAL_QUALIFY_PERCENT="
        << percentage(
               total_qualifying,
               total_composite
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
        << "FINISHED EXPERIMENT 494\n";

    return 0;
}
