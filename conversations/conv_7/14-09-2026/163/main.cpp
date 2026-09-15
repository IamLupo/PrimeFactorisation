#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct CRTPair {
    int m1 = 0;
    int m2 = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 modulus = 0;
    u64 residue = 0;
};

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 candidate = 0;

    u64 j = 0;

    u64 pair_rank = 0;
    u64 candidate_rank = 0;

    CRTPair pair;
};

struct CaseResult {
    u64 N = 0;
    u64 s = 0;

    Hit hit;
};

static std::vector<u64> generate_primes(
    int limit
) {
    std::vector<bool> sieve(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    if (limit >= 0) {
        sieve[0] = false;
    }

    if (limit >= 1) {
        sieve[1] = false;
    }

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
            sq <= static_cast<i128>(n)
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

static bool determinant_one(
    int m1,
    u64 k1,
    int m2,
    u64 k2
) {
    const i128 value =
        static_cast<i128>(m2) *
        static_cast<i128>(k1) -
        static_cast<i128>(m1) *
        static_cast<i128>(k2);

    return value == 1;
}

static std::vector<CRTPair>
generate_crt_pairs(
    const std::vector<int>& multipliers,
    u64 K_LIMIT
) {
    std::vector<CRTPair> pairs;

    for (int m1 : multipliers) {
        for (int m2 : multipliers) {
            for (
                u64 k1 = 2;
                k1 <= K_LIMIT;
                ++k1
            ) {
                const i128 numerator =
                    static_cast<i128>(m2) *
                    static_cast<i128>(k1) -
                    static_cast<i128>(1);

                const i128 denominator =
                    static_cast<i128>(m1);

                if (numerator <= 0) {
                    continue;
                }

                if (
                    numerator %
                    denominator != 0
                ) {
                    continue;
                }

                const i128 k2_i =
                    numerator /
                    denominator;

                if (
                    k2_i < 2 ||
                    k2_i >
                        static_cast<i128>(
                            K_LIMIT
                        )
                ) {
                    continue;
                }

                const u64 k2 =
                    static_cast<u64>(k2_i);

                if (
                    !determinant_one(
                        m1,
                        k1,
                        m2,
                        k2
                    )
                ) {
                    continue;
                }

                if (
                    std::gcd(k1, k2) != 1
                ) {
                    continue;
                }

                const u64 modulus =
                    k1 * k2;

                if (modulus == 0) {
                    continue;
                }

                CRTPair pair;

                pair.m1 = m1;
                pair.m2 = m2;

                pair.k1 = k1;
                pair.k2 = k2;

                pair.modulus = modulus;
                pair.residue =
                    (k1 + k2) %
                    modulus;

                pairs.push_back(pair);
            }
        }
    }

    /*
        Deduplicate identical (residue, modulus)
        classes. Keep the first structural pair.
    */

    std::set<
        std::pair<u64, u64>
    > seen;

    std::vector<CRTPair> unique_pairs;

    unique_pairs.reserve(
        pairs.size()
    );

    for (const CRTPair& pair : pairs) {
        const auto key =
            std::make_pair(
                pair.residue,
                pair.modulus
            );

        if (
            seen.insert(key).second
        ) {
            unique_pairs.push_back(pair);
        }
    }

    return unique_pairs;
}

static Hit blind_search(
    u64 N,
    const std::vector<CRTPair>& pairs,
    u64 J_LIMIT
) {
    Hit hit;

    u64 candidate_rank = 0;

    /*
        j is searched first.

        Therefore the first hit has the smallest
        j among all generated determinant-one
        structures, subject to duplicate removal.
    */

    for (
        u64 j = 0;
        j <= J_LIMIT;
        ++j
    ) {
        std::set<u64> candidates_this_j;

        for (
            std::size_t pair_index = 0;
            pair_index < pairs.size();
            ++pair_index
        ) {
            const CRTPair& pair =
                pairs[pair_index];

            const i128 candidate_i =
                static_cast<i128>(pair.k1) +
                static_cast<i128>(pair.k2) +
                static_cast<i128>(j) *
                static_cast<i128>(pair.k1) *
                static_cast<i128>(pair.k2);

            if (
                candidate_i < 2 ||
                candidate_i >
                    static_cast<i128>(
                        N
                    )
            ) {
                continue;
            }

            const u64 candidate =
                static_cast<u64>(
                    candidate_i
                );

            candidates_this_j.insert(
                candidate
            );
        }

        for (u64 candidate :
             candidates_this_j) {

            ++candidate_rank;

            const u64 g =
                std::gcd(
                    N,
                    candidate
                );

            if (
                g != 1 &&
                g != N
            ) {
                /*
                    Find one structural representation
                    of this candidate for reporting.
                */
                for (
                    std::size_t pair_index = 0;
                    pair_index < pairs.size();
                    ++pair_index
                ) {
                    const CRTPair& pair =
                        pairs[pair_index];

                    const i128 candidate_i =
                        static_cast<i128>(
                            pair.k1
                        ) +
                        static_cast<i128>(
                            pair.k2
                        ) +
                        static_cast<i128>(j) *
                        static_cast<i128>(
                            pair.k1
                        ) *
                        static_cast<i128>(
                            pair.k2
                        );

                    if (
                        candidate_i ==
                        static_cast<i128>(
                            candidate
                        )
                    ) {
                        hit.valid = true;
                        hit.factor = g;
                        hit.candidate = candidate;
                        hit.j = j;

                        hit.pair_rank =
                            static_cast<u64>(
                                pair_index + 1
                            );

                        hit.candidate_rank =
                            candidate_rank;

                        hit.pair = pair;

                        return hit;
                    }
                }
            }
        }
    }

    return hit;
}

static std::string pair_to_string(
    const CRTPair& pair
) {
    std::string out;

    out +=
        "m1=" +
        std::to_string(pair.m1);

    out +=
        " m2=" +
        std::to_string(pair.m2);

    out +=
        " k1=" +
        std::to_string(pair.k1);

    out +=
        " k2=" +
        std::to_string(pair.k2);

    out +=
        " M=" +
        std::to_string(pair.modulus);

    out +=
        " residue=" +
        std::to_string(pair.residue);

    return out;
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
    constexpr int EXPERIMENT = 455;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 2000;

    constexpr u64 K_LIMIT = 1000;

    constexpr u64 J_LIMIT = 2000;

    constexpr std::uint64_t SEED =
        0x455455455ULL;

    const std::vector<int> MULTIPLIERS = {
        1, 2, 3, 4, 5, 6, 7
    };

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
        << "J_LIMIT="
        << J_LIMIT
        << "\n";

    std::cout
        << "MULTIPLIERS=1,2,3,4,5,6,7"
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

    const std::vector<CRTPair> pairs =
        generate_crt_pairs(
            MULTIPLIERS,
            K_LIMIT
        );

    std::cout
        << "UNIQUE_CRT_PAIRS="
        << pairs.size()
        << "\n";

    const auto cases =
        generate_cases(
            primes,
            PRIME_MIN,
            CASE_COUNT,
            SEED
        );

    std::size_t hits = 0;
    std::size_t misses = 0;

    std::size_t j_zero = 0;
    std::size_t j_one = 0;
    std::size_t j_le_2 = 0;
    std::size_t j_le_5 = 0;
    std::size_t j_le_10 = 0;
    std::size_t j_le_20 = 0;
    std::size_t j_le_50 = 0;
    std::size_t j_le_100 = 0;
    std::size_t j_le_500 = 0;
    std::size_t j_le_1000 = 0;
    std::size_t j_gt_1000 = 0;

    u64 min_j = UINT64_MAX;
    u64 max_j = 0;
    u64 total_j = 0;

    u64 min_candidate_rank =
        UINT64_MAX;

    u64 max_candidate_rank = 0;
    u64 total_candidate_rank = 0;

    u64 first_factor_hits = 0;

    bool printed_first_hit = false;
    bool printed_first_miss = false;
    bool printed_first_large_j = false;

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

        const CaseResult result{
            N,
            s,
            blind_search(
                N,
                pairs,
                J_LIMIT
            )
        };

        print_progress(
            i,
            cases.size()
        );

        if (!result.hit.valid) {
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
            }

            continue;
        }

        ++hits;

        const Hit& hit =
            result.hit;

        min_j =
            std::min(
                min_j,
                hit.j
            );

        max_j =
            std::max(
                max_j,
                hit.j
            );

        total_j +=
            hit.j;

        min_candidate_rank =
            std::min(
                min_candidate_rank,
                hit.candidate_rank
            );

        max_candidate_rank =
            std::max(
                max_candidate_rank,
                hit.candidate_rank
            );

        total_candidate_rank +=
            hit.candidate_rank;

        if (
            hit.factor == p ||
            hit.factor == q
        ) {
            ++first_factor_hits;
        }

        if (hit.j == 0) {
            ++j_zero;
        }

        if (hit.j == 1) {
            ++j_one;
        }

        if (hit.j <= 2) {
            ++j_le_2;
        }

        if (hit.j <= 5) {
            ++j_le_5;
        }

        if (hit.j <= 10) {
            ++j_le_10;
        }

        if (hit.j <= 20) {
            ++j_le_20;
        }

        if (hit.j <= 50) {
            ++j_le_50;
        }

        if (hit.j <= 100) {
            ++j_le_100;
        }

        if (hit.j <= 500) {
            ++j_le_500;
        }

        if (hit.j <= 1000) {
            ++j_le_1000;
        }

        if (hit.j > 1000) {
            ++j_gt_1000;

            if (!printed_first_large_j) {
                printed_first_large_j = true;

                std::cout
                    << "\nFIRST_J_GT_1000\n";

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
                    << "J="
                    << hit.j
                    << "\n";

                std::cout
                    << "CANDIDATE="
                    << hit.candidate
                    << "\n";

                std::cout
                    << "FACTOR="
                    << hit.factor
                    << "\n";

                std::cout
                    << "PAIR="
                    << pair_to_string(
                        hit.pair
                    )
                    << "\n";
            }
        }

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
                << "J="
                << hit.j
                << "\n";

            std::cout
                << "CANDIDATE="
                << hit.candidate
                << "\n";

            std::cout
                << "FACTOR="
                << hit.factor
                << "\n";

            std::cout
                << "PAIR="
                << pair_to_string(
                    hit.pair
                )
                << "\n";

            std::cout
                << "PAIR_RANK="
                << hit.pair_rank
                << "\n";

            std::cout
                << "CANDIDATE_RANK="
                << hit.candidate_rank
                << "\n";

            std::cout
                << "IDENTITY_CHECK=";

            const i128 lhs =
                (
                    static_cast<i128>(
                        hit.j
                    ) *
                    static_cast<i128>(
                        hit.pair.k1
                    ) +
                    static_cast<i128>(1)
                ) *
                (
                    static_cast<i128>(
                        hit.j
                    ) *
                    static_cast<i128>(
                        hit.pair.k2
                    ) +
                    static_cast<i128>(1)
                );

            const i128 rhs =
                static_cast<i128>(
                    hit.j
                ) *
                static_cast<i128>(
                    hit.candidate
                ) +
                static_cast<i128>(1);

            std::cout
                << (
                    lhs == rhs
                        ? "PASS"
                        : "FAIL"
                )
                << "\n";
        }
    }

    if (min_j == UINT64_MAX) {
        min_j = 0;
    }

    if (
        min_candidate_rank ==
        UINT64_MAX
    ) {
        min_candidate_rank = 0;
    }

    const double average_j =
        hits == 0
            ? 0.0
            : static_cast<double>(
                  total_j
              ) /
              static_cast<double>(
                  hits
              );

    const double average_candidate_rank =
        hits == 0
            ? 0.0
            : static_cast<double>(
                  total_candidate_rank
              ) /
              static_cast<double>(
                  hits
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
        << "TRUE_FACTOR_HITS="
        << first_factor_hits
        << "\n";

    std::cout
        << "J_ZERO="
        << j_zero
        << "\n";

    std::cout
        << "J_ONE="
        << j_one
        << "\n";

    std::cout
        << "J_LE_2="
        << j_le_2
        << "\n";

    std::cout
        << "J_LE_5="
        << j_le_5
        << "\n";

    std::cout
        << "J_LE_10="
        << j_le_10
        << "\n";

    std::cout
        << "J_LE_20="
        << j_le_20
        << "\n";

    std::cout
        << "J_LE_50="
        << j_le_50
        << "\n";

    std::cout
        << "J_LE_100="
        << j_le_100
        << "\n";

    std::cout
        << "J_LE_500="
        << j_le_500
        << "\n";

    std::cout
        << "J_LE_1000="
        << j_le_1000
        << "\n";

    std::cout
        << "J_GT_1000="
        << j_gt_1000
        << "\n";

    std::cout
        << "MIN_J="
        << min_j
        << "\n";

    std::cout
        << "MAX_J="
        << max_j
        << "\n";

    std::cout
        << "AVERAGE_J="
        << average_j
        << "\n";

    std::cout
        << "MIN_CANDIDATE_RANK="
        << min_candidate_rank
        << "\n";

    std::cout
        << "MAX_CANDIDATE_RANK="
        << max_candidate_rank
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_RANK="
        << average_candidate_rank
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
