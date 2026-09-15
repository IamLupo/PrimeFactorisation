#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 k = 0;

    u64 residue = 0;
    u64 modulus = 0;

    u64 candidate_index = 0;
    u64 candidate_tests = 0;
};

struct StructuralInfo {
    u64 min_k_p = 0;
    u64 min_k_q = 0;

    bool p_has_small_k = false;
    bool q_has_small_k = false;
};

static std::vector<u64> generate_primes(
    int limit
) {
    std::vector<bool> sieve(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    sieve[0] = false;
    sieve[1] = false;

    for (
        int p = 2;
        static_cast<long long>(p) * p <= limit;
        ++p
    ) {
        if (!sieve[p]) {
            continue;
        }

        for (
            int x = p * p;
            x <= limit;
            x += p
        ) {
            sieve[x] = false;
        }
    }

    std::vector<u64> primes;

    for (int x = 2; x <= limit; ++x) {
        if (sieve[x]) {
            primes.push_back(
                static_cast<u64>(x)
            );
        }
    }

    return primes;
}

static u64 integer_sqrt(
    u64 n
) {
    u64 lo = 0;
    u64 hi = n;
    u64 answer = 0;

    while (lo <= hi) {
        const u64 mid =
            lo + (hi - lo) / 2;

        const i128 sq =
            static_cast<i128>(mid) *
            static_cast<i128>(mid);

        if (
            sq <=
            static_cast<i128>(n)
        ) {
            answer = mid;
            lo = mid + 1;
        } else {
            if (mid == 0) {
                break;
            }

            hi = mid - 1;
        }
    }

    return answer;
}

/*
    m = 1 CRT family:

        r ≡ 1  (mod k)
        r ≡ -1 (mod k+1)

    Therefore

        r ≡ 2k+1 (mod k(k+1)).
*/
static Hit blind_search(
    u64 N,
    u64 s,
    u64 K_LIMIT
) {
    Hit hit;

    for (
        u64 k = 2;
        k <= K_LIMIT;
        ++k
    ) {
        const u64 residue =
            2 * k + 1;

        const u64 modulus =
            k * (k + 1);

        if (residue > s) {
            continue;
        }

        /*
            Candidates are

                residue + j*modulus

            up to sqrt(N).
        */

        const u64 count =
            1 +
            (s - residue) /
            modulus;

        for (
            u64 j = 0;
            j < count;
            ++j
        ) {
            const u64 candidate =
                residue +
                j * modulus;

            if (candidate < 2) {
                continue;
            }

            ++hit.candidate_tests;

            const u64 g =
                std::gcd(
                    N,
                    candidate
                );

            if (
                g != 1 &&
                g != N
            ) {
                hit.valid = true;
                hit.factor = g;
                hit.k = k;
                hit.residue = residue;
                hit.modulus = modulus;
                hit.candidate_index = j;

                return hit;
            }
        }
    }

    return hit;
}

/*
    This is analysis only.

    For a prime r, the m=1 CRT condition is

        r ≡ 1  (mod k)
        r ≡ -1 (mod k+1).

    Equivalently:

        r % k     == 1
        r % (k+1) == k.
*/
static u64 structural_min_k(
    u64 r,
    u64 K_LIMIT
) {
    for (
        u64 k = 2;
        k <= K_LIMIT;
        ++k
    ) {
        if (
            r % k == 1 &&
            r % (k + 1) == k
        ) {
            return k;
        }
    }

    return 0;
}

static StructuralInfo analyse_structure(
    u64 p,
    u64 q,
    u64 K_LIMIT
) {
    StructuralInfo info;

    info.min_k_p =
        structural_min_k(
            p,
            K_LIMIT
        );

    info.min_k_q =
        structural_min_k(
            q,
            K_LIMIT
        );

    info.p_has_small_k =
        info.min_k_p != 0;

    info.q_has_small_k =
        info.min_k_q != 0;

    return info;
}

static std::vector<
    std::pair<u64, u64>
> generate_cases(
    const std::vector<u64>& primes,
    int prime_min,
    std::size_t count,
    std::uint64_t seed
) {
    std::vector<u64> usable;

    for (u64 p : primes) {
        if (
            p >=
            static_cast<u64>(prime_min)
        ) {
            usable.push_back(p);
        }
    }

    std::mt19937_64 rng(seed);

    std::vector<
        std::pair<u64, u64>
    > cases;

    cases.reserve(count);

    for (
        std::size_t i = 0;
        i < count;
        ++i
    ) {
        std::size_t a =
            static_cast<std::size_t>(
                rng() % usable.size()
            );

        std::size_t b =
            static_cast<std::size_t>(
                rng() % usable.size()
            );

        while (a == b) {
            b =
                static_cast<std::size_t>(
                    rng() % usable.size()
                );
        }

        if (usable[a] > usable[b]) {
            std::swap(a, b);
        }

        cases.emplace_back(
            usable[a],
            usable[b]
        );
    }

    return cases;
}

static void print_progress(
    std::size_t index,
    std::size_t total
) {
    if (
        index == 0 ||
        index % 100 == 0 ||
        index + 1 == total
    ) {
        std::cout
            << "PROGRESS "
            << index + 1
            << "/"
            << total
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 451;

    constexpr int PRIME_LIMIT = 100000;

    /*
        Exclude small primes so that the trivial

            k = (p-1)/2

        solution is generally far outside K_LIMIT.
    */
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 1000;

    constexpr u64 K_LIMIT = 500;

    constexpr std::uint64_t SEED =
        0x451451451ULL;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT
        << "\n";

    std::cout
        << "PRIME_MIN="
        << PRIME_MIN
        << "\n";

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << "\n";

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << "\n";

    std::cout
        << "FAMILY=m=1"
        << "\n";

    std::cout
        << "SEED="
        << SEED
        << "\n";

    const std::vector<u64> primes =
        generate_primes(
            PRIME_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    const std::vector<
        std::pair<u64, u64>
    > cases =
        generate_cases(
            primes,
            PRIME_MIN,
            CASE_COUNT,
            SEED
        );

    std::size_t hits = 0;
    std::size_t misses = 0;

    std::size_t structural_p_hits = 0;
    std::size_t structural_q_hits = 0;

    std::size_t structural_either_hits = 0;

    u64 min_blind_k = UINT64_MAX;
    u64 max_blind_k = 0;

    u64 min_structural_k = UINT64_MAX;
    u64 max_structural_k = 0;

    u64 total_candidate_tests = 0;
    u64 max_candidate_tests = 0;
    u64 min_candidate_tests = UINT64_MAX;

    u64 total_s = 0;

    bool printed_first_hit = false;
    bool printed_first_miss = false;
    bool printed_first_large_k = false;

    for (
        std::size_t i = 0;
        i < cases.size();
        ++i
    ) {
        const u64 p =
            cases[i].first;

        const u64 q =
            cases[i].second;

        const u64 N =
            p * q;

        const u64 s =
            integer_sqrt(N);

        const Hit hit =
            blind_search(
                N,
                s,
                K_LIMIT
            );

        const StructuralInfo structure =
            analyse_structure(
                p,
                q,
                K_LIMIT
            );

        total_s += s;

        if (structure.p_has_small_k) {
            ++structural_p_hits;

            min_structural_k =
                std::min(
                    min_structural_k,
                    structure.min_k_p
                );

            max_structural_k =
                std::max(
                    max_structural_k,
                    structure.min_k_p
                );
        }

        if (structure.q_has_small_k) {
            ++structural_q_hits;

            min_structural_k =
                std::min(
                    min_structural_k,
                    structure.min_k_q
                );

            max_structural_k =
                std::max(
                    max_structural_k,
                    structure.min_k_q
                );
        }

        if (
            structure.p_has_small_k ||
            structure.q_has_small_k
        ) {
            ++structural_either_hits;
        }

        total_candidate_tests +=
            hit.candidate_tests;

        min_candidate_tests =
            std::min(
                min_candidate_tests,
                hit.candidate_tests
            );

        max_candidate_tests =
            std::max(
                max_candidate_tests,
                hit.candidate_tests
            );

        if (hit.valid) {
            ++hits;

            min_blind_k =
                std::min(
                    min_blind_k,
                    hit.k
                );

            max_blind_k =
                std::max(
                    max_blind_k,
                    hit.k
                );

            if (!printed_first_hit) {
                printed_first_hit = true;

                std::cout
                    << "\nFIRST_HIT\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << N
                    << " S="
                    << s
                    << "\n";

                std::cout
                    << "BLIND_K="
                    << hit.k
                    << "\n";

                std::cout
                    << "BLIND_FACTOR="
                    << hit.factor
                    << "\n";

                std::cout
                    << "BLIND_RESIDUE="
                    << hit.residue
                    << "\n";

                std::cout
                    << "BLIND_MODULUS="
                    << hit.modulus
                    << "\n";

                std::cout
                    << "BLIND_CANDIDATE_INDEX="
                    << hit.candidate_index
                    << "\n";

                std::cout
                    << "BLIND_CANDIDATE_TESTS="
                    << hit.candidate_tests
                    << "\n";

                std::cout
                    << "STRUCTURAL_MIN_K_P="
                    << structure.min_k_p
                    << "\n";

                std::cout
                    << "STRUCTURAL_MIN_K_Q="
                    << structure.min_k_q
                    << "\n";
            }

            if (
                hit.k >= 100 &&
                !printed_first_large_k
            ) {
                printed_first_large_k = true;

                std::cout
                    << "\nFIRST_BLIND_HIT_K_GE_100\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << N
                    << " S="
                    << s
                    << "\n";

                std::cout
                    << "BLIND_K="
                    << hit.k
                    << "\n";

                std::cout
                    << "BLIND_FACTOR="
                    << hit.factor
                    << "\n";

                std::cout
                    << "STRUCTURAL_MIN_K_P="
                    << structure.min_k_p
                    << "\n";

                std::cout
                    << "STRUCTURAL_MIN_K_Q="
                    << structure.min_k_q
                    << "\n";
            }
        } else {
            ++misses;

            if (!printed_first_miss) {
                printed_first_miss = true;

                std::cout
                    << "\nFIRST_MISS\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << N
                    << " S="
                    << s
                    << "\n";

                std::cout
                    << "STRUCTURAL_MIN_K_P="
                    << structure.min_k_p
                    << "\n";

                std::cout
                    << "STRUCTURAL_MIN_K_Q="
                    << structure.min_k_q
                    << "\n";

                std::cout
                    << "CANDIDATE_TESTS="
                    << hit.candidate_tests
                    << "\n";
            }
        }

        print_progress(
            i,
            cases.size()
        );
    }

    if (
        min_blind_k == UINT64_MAX
    ) {
        min_blind_k = 0;
    }

    if (
        min_structural_k == UINT64_MAX
    ) {
        min_structural_k = 0;
    }

    if (
        min_candidate_tests == UINT64_MAX
    ) {
        min_candidate_tests = 0;
    }

    const double average_candidate_tests =
        cases.empty()
            ? 0.0
            : static_cast<double>(
                  total_candidate_tests
              ) /
              static_cast<double>(
                  cases.size()
              );

    const double average_s =
        cases.empty()
            ? 0.0
            : static_cast<double>(
                  total_s
              ) /
              static_cast<double>(
                  cases.size()
              );

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "HITS="
        << hits
        << "\n";

    std::cout
        << "MISSES="
        << misses
        << "\n";

    std::cout
        << "STRUCTURAL_P_HITS="
        << structural_p_hits
        << "\n";

    std::cout
        << "STRUCTURAL_Q_HITS="
        << structural_q_hits
        << "\n";

    std::cout
        << "STRUCTURAL_EITHER_HITS="
        << structural_either_hits
        << "\n";

    std::cout
        << "MIN_BLIND_K="
        << min_blind_k
        << "\n";

    std::cout
        << "MAX_BLIND_K="
        << max_blind_k
        << "\n";

    std::cout
        << "MIN_STRUCTURAL_K="
        << min_structural_k
        << "\n";

    std::cout
        << "MAX_STRUCTURAL_K="
        << max_structural_k
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_candidate_tests
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_candidate_tests
        << "\n";

    std::cout
        << "TOTAL_CANDIDATE_TESTS="
        << total_candidate_tests
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << average_candidate_tests
        << "\n";

    std::cout
        << "AVERAGE_S="
        << average_s
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
