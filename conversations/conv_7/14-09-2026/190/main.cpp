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

struct Representation {
    bool found = false;
    u64 M = 0;
    u64 a = 0;
    u64 j = 0;
    int k1 = 0;
    int k2 = 0;
};

struct Coverage {
    u64 limit;
    u64 count = 0;
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
                    static_cast<int>(numerator / m1);

                if (k2 < 2 || k2 < k1 || k2 > k_limit) {
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
            static_cast<u64>(primes[dist(rng)]);

        u64 q =
            static_cast<u64>(primes[dist(rng)]);

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

Representation find_minimum_M(
    u64 p,
    const std::vector<Pair>& pairs) {

    Representation best;

    for (const Pair& pair : pairs) {

        if (pair.a > p) {
            continue;
        }

        if (p % pair.M != pair.a % pair.M) {
            continue;
        }

        const u64 j =
            (p - pair.a) / pair.M;

        if (!best.found ||
            pair.M < best.M ||
            (pair.M == best.M &&
             j < best.j)) {

            best.found = true;
            best.M = pair.M;
            best.a = pair.a;
            best.j = j;
            best.k1 = pair.k1;
            best.k2 = pair.k2;
        }
    }

    return best;
}

bool verify_representation(
    u64 p,
    const Representation& rep) {

    if (!rep.found) {
        return false;
    }

    const u64 reconstructed =
        rep.a + rep.j * rep.M;

    const u64 lhs =
        p % rep.M;

    const u64 rhs =
        rep.a % rep.M;

    return reconstructed == p && lhs == rhs;
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
    constexpr int CASE_COUNT = 500;

    constexpr int PRIME_LOW = 10000;
    constexpr int PRIME_HIGH = 100000;
    constexpr int PRIME_LIMIT = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_LIMIT = 7;

    const std::vector<u64> M_LIMITS = {
        6,
        10,
        20,
        50,
        100,
        250,
        500,
        1000,
        2000,
        5000,
        10000,
        20000,
        50000,
        100000,
        250000,
        500000,
        1000000
    };

    std::cout << "START EXPERIMENT 483\n";

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
            M_LIMIT
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
        4831234567ULL
    );

    std::vector<Coverage> coverage;

    for (u64 limit : M_LIMITS) {
        coverage.push_back({
            limit,
            0
        });
    }

    u64 found = 0;
    u64 missing = 0;

    u64 total_M = 0;
    u64 max_M = 0;

    u64 total_a = 0;
    u64 total_j = 0;

    u64 total_k1 = 0;
    u64 total_k2 = 0;

    u64 total_p = 0;
    u64 total_distance = 0;

    u64 total_Mj = 0;

    u64 verify_ok = 0;
    u64 verify_fail = 0;

    u64 M_le_10 = 0;
    u64 M_le_100 = 0;
    u64 M_le_1000 = 0;
    u64 M_le_10000 = 0;
    u64 M_le_100000 = 0;

    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;
    u64 first_s = 0;
    Representation first_rep;

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

        const Representation rep =
            find_minimum_M(
                data.p,
                pairs
            );

        if (case_id == 0) {
            first_N = data.N;
            first_p = data.p;
            first_q = data.q;
            first_s = data.s;
            first_rep = rep;
        }

        if (!rep.found) {
            ++missing;
        } else {
            ++found;

            total_M += rep.M;
            max_M =
                std::max(
                    max_M,
                    rep.M
                );

            total_a += rep.a;
            total_j += rep.j;

            total_k1 +=
                static_cast<u64>(
                    rep.k1
                );

            total_k2 +=
                static_cast<u64>(
                    rep.k2
                );

            total_p += data.p;

            total_distance +=
                data.s - data.p;

            total_Mj +=
                rep.M * rep.j;

            if (verify_representation(
                    data.p,
                    rep)) {

                ++verify_ok;
            } else {
                ++verify_fail;
            }

            if (rep.M <= 10) {
                ++M_le_10;
            }

            if (rep.M <= 100) {
                ++M_le_100;
            }

            if (rep.M <= 1000) {
                ++M_le_1000;
            }

            if (rep.M <= 10000) {
                ++M_le_10000;
            }

            if (rep.M <= 100000) {
                ++M_le_100000;
            }

            for (auto& entry : coverage) {
                if (rep.M <= entry.limit) {
                    ++entry.count;
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
        << "FIRST_FOUND="
        << (first_rep.found ? 1 : 0)
        << "\n";

    if (first_rep.found) {
        std::cout
            << "FIRST_MIN_M="
            << first_rep.M
            << "\n";

        std::cout
            << "FIRST_K1="
            << first_rep.k1
            << "\n";

        std::cout
            << "FIRST_K2="
            << first_rep.k2
            << "\n";

        std::cout
            << "FIRST_A="
            << first_rep.a
            << "\n";

        std::cout
            << "FIRST_J="
            << first_rep.j
            << "\n";

        std::cout
            << "FIRST_P_MOD_M="
            << first_p % first_rep.M
            << "\n";

        std::cout
            << "FIRST_A_MOD_M="
            << first_rep.a % first_rep.M
            << "\n";

        std::cout
            << "FIRST_RECONSTRUCTED_P="
            << first_rep.a +
                   first_rep.j *
                   first_rep.M
            << "\n";
    }

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << "\n";

    std::cout
        << "REPRESENTATION_FOUND="
        << found
        << "\n";

    std::cout
        << "REPRESENTATION_MISSING="
        << missing
        << "\n";

    std::cout
        << "VERIFY_OK="
        << verify_ok
        << "\n";

    std::cout
        << "VERIFY_FAIL="
        << verify_fail
        << "\n";

    std::cout
        << "AVG_MIN_M="
        << average_u64(
               total_M,
               found
           )
        << "\n";

    std::cout
        << "MAX_MIN_M="
        << max_M
        << "\n";

    std::cout
        << "AVG_A="
        << average_u64(
               total_a,
               found
           )
        << "\n";

    std::cout
        << "AVG_J="
        << average_u64(
               total_j,
               found
           )
        << "\n";

    std::cout
        << "AVG_K1="
        << average_u64(
               total_k1,
               found
           )
        << "\n";

    std::cout
        << "AVG_K2="
        << average_u64(
               total_k2,
               found
           )
        << "\n";

    std::cout
        << "AVG_P="
        << average_u64(
               total_p,
               found
           )
        << "\n";

    std::cout
        << "AVG_DISTANCE="
        << average_u64(
               total_distance,
               found
           )
        << "\n";

    std::cout
        << "AVG_M_TIMES_J="
        << average_u64(
               total_Mj,
               found
           )
        << "\n";

    std::cout
        << "MIN_M_LE_10="
        << M_le_10
        << "\n";

    std::cout
        << "MIN_M_LE_100="
        << M_le_100
        << "\n";

    std::cout
        << "MIN_M_LE_1000="
        << M_le_1000
        << "\n";

    std::cout
        << "MIN_M_LE_10000="
        << M_le_10000
        << "\n";

    std::cout
        << "MIN_M_LE_100000="
        << M_le_100000
        << "\n";

    for (const auto& entry : coverage) {
        std::cout
            << "M_LIMIT="
            << entry.limit
            << " COVERAGE="
            << entry.count
            << "\n";
    }

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 483\n";

    return 0;
}
