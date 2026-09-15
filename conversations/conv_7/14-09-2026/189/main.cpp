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

/*
    Build the determinant-one family

        m2*k1 - m1*k2 = 1

    with

        1 <= m1,m2 <= 7
        1 <= k1 <= k2 <= 1000
        gcd(k1,k2)=1.

    Only (k1,k2) matter for the representation

        p = k1 + k2 + j*k1*k2.
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
    Find the representation of p with the smallest

        M = k1*k2.

    The defining congruence is

        p mod M = k1 + k2.

    Since p = a + j*M, this is equivalent to an exact
    non-negative integer j.
*/
Representation find_minimum_M_representation(
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
            pair.M < best.M ||
            (pair.M == best.M &&
             j < best.j) ||
            (pair.M == best.M &&
             j == best.j &&
             pair.a < best.a)) {

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

/*
    Find the minimum M among representations satisfying
    the congruence directly:

        p mod M = k1+k2.

    This function is deliberately separate from the exact
    reconstruction test so the congruence itself is explicit.
*/
Representation find_minimum_M_by_congruence(
    u64 p,
    const std::vector<Pair>& pairs) {

    Representation best;

    for (const Pair& pair : pairs) {

        if (p < pair.a) {
            continue;
        }

        if (p % pair.M != pair.a) {
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

double average_u64(u64 total, u64 count) {
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

    /*
        M thresholds are deliberately concentrated around
        the low range because the question is whether a
        small-M search can localize the factor.
    */
    const std::vector<u64> M_LIMITS = {
        1,
        2,
        3,
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
        20000,
        50000,
        100000,
        250000,
        500000,
        1000000
    };

    std::cout << "START EXPERIMENT 482\n";

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
        4821234567ULL
    );

    std::vector<Coverage> coverage;

    for (u64 limit : M_LIMITS) {
        coverage.push_back({limit, 0});
    }

    u64 representation_found = 0;
    u64 representation_missing = 0;

    u64 M_total = 0;
    u64 M_max = 0;

    u64 a_total = 0;
    u64 j_total = 0;

    u64 k1_total = 0;
    u64 k2_total = 0;

    u64 p_total = 0;
    u64 distance_total = 0;

    u64 M_times_j_total = 0;

    u64 min_M_eq_0 = 0;
    u64 min_M_le_10 = 0;
    u64 min_M_le_100 = 0;
    u64 min_M_le_1000 = 0;
    u64 min_M_le_10000 = 0;
    u64 min_M_le_100000 = 0;

    u64 congruence_matches = 0;
    u64 congruence_mismatches = 0;

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

        const Representation exact_rep =
            find_minimum_M_representation(
                data.p,
                pairs
            );

        const Representation congruence_rep =
            find_minimum_M_by_congruence(
                data.p,
                pairs
            );

        if (case_id == 0) {
            first_N = data.N;
            first_p = data.p;
            first_q = data.q;
            first_s = data.s;
            first_rep = exact_rep;
        }

        if (!exact_rep.found) {
            ++representation_missing;
        } else {
            ++representation_found;

            M_total += exact_rep.M;
            M_max =
                std::max(
                    M_max,
                    exact_rep.M
                );

            a_total += exact_rep.a;
            j_total += exact_rep.j;

            k1_total +=
                static_cast<u64>(
                    exact_rep.k1
                );

            k2_total +=
                static_cast<u64>(
                    exact_rep.k2
                );

            p_total += data.p;

            distance_total +=
                data.s - data.p;

            M_times_j_total +=
                exact_rep.M *
                exact_rep.j;

            if (exact_rep.M == 1) {
                ++min_M_eq_0;
            }

            if (exact_rep.M <= 10) {
                ++min_M_le_10;
            }

            if (exact_rep.M <= 100) {
                ++min_M_le_100;
            }

            if (exact_rep.M <= 1000) {
                ++min_M_le_1000;
            }

            if (exact_rep.M <= 10000) {
                ++min_M_le_10000;
            }

            if (exact_rep.M <= 100000) {
                ++min_M_le_100000;
            }

            for (auto& entry : coverage) {
                if (exact_rep.M <= entry.limit) {
                    ++entry.count;
                }
            }
        }

        if (!exact_rep.found &&
            !congruence_rep.found) {
            /*
                Both methods missing is consistent.
            */
        } else if (exact_rep.found &&
                   congruence_rep.found &&
                   exact_rep.M == congruence_rep.M &&
                   exact_rep.k1 == congruence_rep.k1 &&
                   exact_rep.k2 == congruence_rep.k2 &&
                   exact_rep.a == congruence_rep.a &&
                   exact_rep.j == congruence_rep.j) {
            ++congruence_matches;
        } else {
            ++congruence_mismatches;
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
        << "FIRST_REPRESENTATION_FOUND="
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
            << "FIRST_CONGRUENCE_RESIDUE="
            << first_rep.a
            << "\n";

        std::cout
            << "FIRST_P_MOD_M="
            << first_p % first_rep.M
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
        << representation_found
        << "\n";

    std::cout
        << "REPRESENTATION_MISSING="
        << representation_missing
        << "\n";

    std::cout
        << "CONGRUENCE_MATCHES="
        << congruence_matches
        << "\n";

    std::cout
        << "CONGRUENCE_MISMATCHES="
        << congruence_mismatches
        << "\n";

    if (representation_found > 0) {

        std::cout
            << "AVG_MIN_M="
            << average_u64(
                   M_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "MAX_MIN_M="
            << M_max
            << "\n";

        std::cout
            << "AVG_A="
            << average_u64(
                   a_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "AVG_J="
            << average_u64(
                   j_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "AVG_K1="
            << average_u64(
                   k1_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "AVG_K2="
            << average_u64(
                   k2_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "AVG_P="
            << average_u64(
                   p_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "AVG_DISTANCE="
            << average_u64(
                   distance_total,
                   representation_found
               )
            << "\n";

        std::cout
            << "AVG_M_TIMES_J="
            << average_u64(
                   M_times_j_total,
                   representation_found
               )
            << "\n";
    } else {
        std::cout << "AVG_MIN_M=0\n";
        std::cout << "MAX_MIN_M=0\n";
        std::cout << "AVG_A=0\n";
        std::cout << "AVG_J=0\n";
        std::cout << "AVG_K1=0\n";
        std::cout << "AVG_K2=0\n";
        std::cout << "AVG_P=0\n";
        std::cout << "AVG_DISTANCE=0\n";
        std::cout << "AVG_M_TIMES_J=0\n";
    }

    std::cout
        << "MIN_M_EQ_1_COUNT="
        << min_M_eq_0
        << "\n";

    std::cout
        << "MIN_M_LE_10="
        << min_M_le_10
        << "\n";

    std::cout
        << "MIN_M_LE_100="
        << min_M_le_100
        << "\n";

    std::cout
        << "MIN_M_LE_1000="
        << min_M_le_1000
        << "\n";

    std::cout
        << "MIN_M_LE_10000="
        << min_M_le_10000
        << "\n";

    std::cout
        << "MIN_M_LE_100000="
        << min_M_le_100000
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
        << "FINISHED EXPERIMENT 482\n";

    return 0;
}
