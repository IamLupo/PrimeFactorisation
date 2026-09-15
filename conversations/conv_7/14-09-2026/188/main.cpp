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
    u64 j = 0;
    int k1 = 0;
    int k2 = 0;
    u64 a = 0;
    u64 M = 0;
};

struct HistogramEntry {
    u64 limit;
    u64 count;
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

/*
    Determinant-one structure:

        m2*k1 - m1*k2 = 1

    with

        1 <= m1,m2 <= 7
        1 <= k1 <= k2 <= 1000
        gcd(k1,k2)=1

    This is the same family used in the earlier CRT experiments.
*/
std::vector<Pair> build_determinant_one_pairs(
    int k_limit,
    int m_limit) {

    std::vector<Pair> pairs;

    for (int m1 = 1; m1 <= m_limit; ++m1) {
        for (int m2 = 1; m2 <= m_limit; ++m2) {

            for (int k1 = 1; k1 <= k_limit; ++k1) {

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

                if (k2 < k1 || k2 > k_limit) {
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

/*
    Find the minimum j such that

        p = k1 + k2 + j*k1*k2

    for one of the determinant-one pairs.
*/
Representation find_minimum_representation(
    u64 p,
    const std::vector<Pair>& pairs) {

    Representation best;

    for (const Pair& pair : pairs) {

        if (pair.a > p) {
            continue;
        }

        const u64 remainder =
            p - pair.a;

        if (remainder % pair.M != 0) {
            continue;
        }

        const u64 j =
            remainder / pair.M;

        if (!best.found ||
            j < best.j ||
            (j == best.j &&
             pair.M < best.M) ||
            (j == best.j &&
             pair.M == best.M &&
             pair.a < best.a)) {

            best.found = true;
            best.j = j;
            best.k1 = pair.k1;
            best.k2 = pair.k2;
            best.a = pair.a;
            best.M = pair.M;
        }
    }

    return best;
}

u64 ordinary_factor_distance(
    u64 s,
    u64 p) {

    return s - p;
}

double average_u64(
    u64 value,
    int count) {

    if (count == 0) {
        return 0.0;
    }

    return
        static_cast<double>(value) /
        static_cast<double>(count);
}

int main() {
    constexpr int CASE_COUNT = 500;

    constexpr int PRIME_LOW = 10000;
    constexpr int PRIME_HIGH = 100000;
    constexpr int PRIME_LIMIT = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_LIMIT = 7;

    /*
        Coverage thresholds for the minimum j.
    */
    const std::vector<u64> J_LIMITS = {
        0,
        1,
        2,
        5,
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
        20000
    };

    std::cout << "START EXPERIMENT 481\n";

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
        4811234567ULL
    );

    /*
        Aggregate minimum-j statistics.
    */
    u64 representation_found = 0;
    u64 representation_missing = 0;

    u64 j_total = 0;
    u64 j_max = 0;

    u64 k1_total = 0;
    u64 k2_total = 0;
    u64 a_total = 0;
    u64 M_total = 0;

    u64 p_total = 0;
    u64 s_total = 0;
    u64 distance_total = 0;

    u64 j_times_distance_total = 0;
    u64 M_times_j_total = 0;

    u64 smallest_j_zero = 0;
    u64 smallest_j_one = 0;

    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;
    u64 first_s = 0;
    u64 first_distance = 0;

    Representation first_representation;

    std::vector<HistogramEntry> histogram;

    histogram.reserve(J_LIMITS.size());

    for (u64 limit : J_LIMITS) {
        histogram.push_back({
            limit,
            0
        });
    }

    u64 min_j_0_10 = 0;
    u64 min_j_11_100 = 0;
    u64 min_j_101_1000 = 0;
    u64 min_j_1001_5000 = 0;
    u64 min_j_5001_10000 = 0;
    u64 min_j_above_10000 = 0;

    u64 p_le_s = 0;

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

        if (data.p <= data.s) {
            ++p_le_s;
        }

        const Representation representation =
            find_minimum_representation(
                data.p,
                pairs
            );

        if (case_id == 0) {
            first_N = data.N;
            first_p = data.p;
            first_q = data.q;
            first_s = data.s;

            first_distance =
                ordinary_factor_distance(
                    data.s,
                    data.p
                );

            first_representation =
                representation;
        }

        if (!representation.found) {
            ++representation_missing;
        } else {
            ++representation_found;

            j_total += representation.j;
            j_max =
                std::max(
                    j_max,
                    representation.j
                );

            k1_total +=
                static_cast<u64>(
                    representation.k1
                );

            k2_total +=
                static_cast<u64>(
                    representation.k2
                );

            a_total += representation.a;
            M_total += representation.M;

            p_total += data.p;
            s_total += data.s;

            const u64 distance =
                data.s - data.p;

            distance_total += distance;

            j_times_distance_total +=
                representation.j * distance;

            M_times_j_total +=
                representation.M *
                representation.j;

            if (representation.j == 0) {
                ++smallest_j_zero;
            }

            if (representation.j <= 1) {
                ++smallest_j_one;
            }

            if (representation.j <= 10) {
                ++min_j_0_10;
            } else if (representation.j <= 100) {
                ++min_j_11_100;
            } else if (representation.j <= 1000) {
                ++min_j_101_1000;
            } else if (representation.j <= 5000) {
                ++min_j_1001_5000;
            } else if (representation.j <= 10000) {
                ++min_j_5001_10000;
            } else {
                ++min_j_above_10000;
            }

            for (auto& entry : histogram) {
                if (representation.j <= entry.limit) {
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

    /*
        First case.
    */
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
        << "FIRST_DISTANCE="
        << first_distance
        << "\n";

    std::cout
        << "FIRST_REPRESENTATION_FOUND="
        << (first_representation.found ? 1 : 0)
        << "\n";

    if (first_representation.found) {
        std::cout
            << "FIRST_MIN_J="
            << first_representation.j
            << "\n";

        std::cout
            << "FIRST_K1="
            << first_representation.k1
            << "\n";

        std::cout
            << "FIRST_K2="
            << first_representation.k2
            << "\n";

        std::cout
            << "FIRST_A="
            << first_representation.a
            << "\n";

        std::cout
            << "FIRST_M="
            << first_representation.M
            << "\n";

        std::cout
            << "FIRST_RECONSTRUCTED_P="
            << first_representation.a +
                   first_representation.j *
                   first_representation.M
            << "\n";
    }

    /*
        Overall coverage.
    */
    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << "\n";

    std::cout
        << "REPRESENTATION_FOUND="
        << representation_found
        << "\n";

    std::cout
        << "REPRESENTATION_MISSING="
        << representation_missing
        << "\n";

    std::cout
        << "P_LE_S="
        << p_le_s
        << "\n";

    if (representation_found > 0) {

        std::cout
            << "AVG_MIN_J="
            << average_u64(
                   j_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "MAX_MIN_J="
            << j_max
            << "\n";

        std::cout
            << "AVG_K1="
            << average_u64(
                   k1_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_K2="
            << average_u64(
                   k2_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_A="
            << average_u64(
                   a_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_M="
            << average_u64(
                   M_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_P="
            << average_u64(
                   p_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_S="
            << average_u64(
                   s_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_DISTANCE="
            << average_u64(
                   distance_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_J_TIMES_DISTANCE="
            << average_u64(
                   j_times_distance_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";

        std::cout
            << "AVG_M_TIMES_J="
            << average_u64(
                   M_times_j_total,
                   static_cast<int>(
                       representation_found
                   )
               )
            << "\n";
    } else {
        std::cout << "AVG_MIN_J=0\n";
        std::cout << "MAX_MIN_J=0\n";
        std::cout << "AVG_K1=0\n";
        std::cout << "AVG_K2=0\n";
        std::cout << "AVG_A=0\n";
        std::cout << "AVG_M=0\n";
        std::cout << "AVG_P=0\n";
        std::cout << "AVG_S=0\n";
        std::cout << "AVG_DISTANCE=0\n";
        std::cout << "AVG_J_TIMES_DISTANCE=0\n";
        std::cout << "AVG_M_TIMES_J=0\n";
    }

    std::cout
        << "MIN_J_0_COUNT="
        << smallest_j_zero
        << "\n";

    std::cout
        << "MIN_J_0_OR_1_COUNT="
        << smallest_j_one
        << "\n";

    std::cout
        << "MIN_J_0_10="
        << min_j_0_10
        << "\n";

    std::cout
        << "MIN_J_11_100="
        << min_j_11_100
        << "\n";

    std::cout
        << "MIN_J_101_1000="
        << min_j_101_1000
        << "\n";

    std::cout
        << "MIN_J_1001_5000="
        << min_j_1001_5000
        << "\n";

    std::cout
        << "MIN_J_5001_10000="
        << min_j_5001_10000
        << "\n";

    std::cout
        << "MIN_J_ABOVE_10000="
        << min_j_above_10000
        << "\n";

    /*
        Cumulative minimum-j coverage.
    */
    for (const auto& entry : histogram) {
        std::cout
            << "J_LIMIT="
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
        << "FINISHED EXPERIMENT 481\n";

    return 0;
}
