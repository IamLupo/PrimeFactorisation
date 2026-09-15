#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct CRTPair {
    int m1 = 0;
    int m2 = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 residue = 0;
    u64 modulus = 0;
};

struct BestResult {
    bool valid = false;

    u64 modulus = 0;
    u64 residue = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    int m1 = 0;
    int m2 = 0;

    u64 candidate_count = 0;
    u64 candidates_below_factor = 0;
    u64 candidates_above_factor = 0;
};

struct CaseResult {
    u64 p = 0;
    u64 q = 0;
    u64 N = 0;
    u64 s = 0;

    BestResult best;

    bool modulus_gt_p = false;
    bool modulus_gt_sqrt_p = false;
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

                const i128 modulus_i =
                    static_cast<i128>(k1) *
                    static_cast<i128>(k2);

                if (
                    modulus_i >
                    static_cast<i128>(
                        UINT64_MAX
                    )
                ) {
                    continue;
                }

                const u64 modulus =
                    static_cast<u64>(
                        modulus_i
                    );

                const u64 residue =
                    (k1 + k2) %
                    modulus;

                CRTPair pair;

                pair.m1 = m1;
                pair.m2 = m2;
                pair.k1 = k1;
                pair.k2 = k2;
                pair.modulus = modulus;
                pair.residue = residue;

                pairs.push_back(pair);
            }
        }
    }

    return pairs;
}

static bool crt_matches_factor(
    u64 factor,
    const CRTPair& pair
) {
    return
        factor % pair.modulus ==
        pair.residue;
}

static BestResult analyze_factor(
    u64 factor,
    u64 s,
    const std::vector<CRTPair>& pairs
) {
    BestResult best;

    for (const CRTPair& pair : pairs) {
        if (
            !crt_matches_factor(
                factor,
                pair
            )
        ) {
            continue;
        }

        const u64 M =
            pair.modulus;

        /*
            All positive integers x <= s with

                x ≡ residue (mod M)

            are candidate numbers.

            residue may be zero, in which case
            the first positive candidate is M.
        */

        u64 first = pair.residue;

        if (first == 0) {
            first = M;
        }

        if (first > s) {
            continue;
        }

        const u64 total =
            1 +
            (s - first) / M;

        /*
            Count candidates strictly below factor.
        */

        u64 below = 0;

        if (factor > first) {
            below =
                (factor - 1 - first) /
                M + 1;
        }

        /*
            Number at or above factor.
            The factor itself is always included.
        */

        const u64 above =
            total - below - 1;

        if (
            !best.valid ||
            M > best.modulus ||
            (
                M == best.modulus &&
                total < best.candidate_count
            )
        ) {
            best.valid = true;
            best.modulus = M;
            best.residue = pair.residue;

            best.k1 = pair.k1;
            best.k2 = pair.k2;

            best.m1 = pair.m1;
            best.m2 = pair.m2;

            best.candidate_count = total;
            best.candidates_below_factor =
                below;
            best.candidates_above_factor =
                above;
        }
    }

    return best;
}

static CaseResult run_case(
    u64 p,
    u64 q,
    const std::vector<CRTPair>& pairs
) {
    CaseResult result;

    result.p = p;
    result.q = q;

    result.N = p * q;

    result.s =
        integer_sqrt(
            result.N
        );

    result.best =
        analyze_factor(
            p,
            result.s,
            pairs
        );

    if (result.best.valid) {
        result.modulus_gt_p =
            result.best.modulus > p;

        result.modulus_gt_sqrt_p =
            result.best.modulus >
            integer_sqrt(p);
    }

    return result;
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

static std::string pair_to_string(
    const BestResult& best
) {
    std::string result;

    result +=
        "m1=" +
        std::to_string(best.m1);

    result +=
        " m2=" +
        std::to_string(best.m2);

    result +=
        " k1=" +
        std::to_string(best.k1);

    result +=
        " k2=" +
        std::to_string(best.k2);

    result +=
        " residue=" +
        std::to_string(best.residue);

    result +=
        " M=" +
        std::to_string(best.modulus);

    return result;
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
    constexpr int EXPERIMENT = 452;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 2000;

    constexpr u64 K_LIMIT = 1000;

    constexpr std::uint64_t SEED =
        0x452452452ULL;

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
        << "DETERMINANT_ONE_PAIRS="
        << pairs.size()
        << "\n";

    std::size_t hits = 0;

    std::size_t modulus_gt_p = 0;
    std::size_t modulus_gt_sqrt_p = 0;

    std::size_t candidate_count_one = 0;
    std::size_t candidate_count_le_2 = 0;
    std::size_t candidate_count_le_10 = 0;
    std::size_t candidate_count_le_100 = 0;
    std::size_t candidate_count_le_1000 = 0;

    u64 min_M = UINT64_MAX;
    u64 max_M = 0;

    u64 min_candidates = UINT64_MAX;
    u64 max_candidates = 0;

    u64 min_M_over_p = UINT64_MAX;
    u64 max_M_over_p = 0;

    u64 total_candidates = 0;
    u64 total_M = 0;

    std::size_t printed_large_M = 0;

    bool printed_first =
        false;

    bool printed_first_one_candidate =
        false;

    bool printed_first_M_gt_p =
        false;

    bool printed_first_M_gt_sqrt_p =
        false;

    const auto cases =
        generate_cases(
            primes,
            PRIME_MIN,
            CASE_COUNT,
            SEED
        );

    for (
        std::size_t i = 0;
        i < cases.size();
        ++i
    ) {
        const u64 p =
            cases[i].first;

        const u64 q =
            cases[i].second;

        const CaseResult result =
            run_case(
                p,
                q,
                pairs
            );

        print_progress(
            i,
            cases.size()
        );

        if (!result.best.valid) {
            continue;
        }

        ++hits;

        const u64 M =
            result.best.modulus;

        const u64 candidates =
            result.best.candidate_count;

        min_M =
            std::min(
                min_M,
                M
            );

        max_M =
            std::max(
                max_M,
                M
            );

        min_candidates =
            std::min(
                min_candidates,
                candidates
            );

        max_candidates =
            std::max(
                max_candidates,
                candidates
            );

        total_candidates +=
            candidates;

        total_M += M;

        const u64 M_over_p =
            M / p;

        min_M_over_p =
            std::min(
                min_M_over_p,
                M_over_p
            );

        max_M_over_p =
            std::max(
                max_M_over_p,
                M_over_p
            );

        if (
            result.modulus_gt_p
        ) {
            ++modulus_gt_p;

            if (!printed_first_M_gt_p) {
                printed_first_M_gt_p = true;

                std::cout
                    << "\nFIRST_M_GT_P\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << result.N
                    << " S="
                    << result.s
                    << "\n";

                std::cout
                    << pair_to_string(
                        result.best
                    )
                    << "\n";

                std::cout
                    << "CANDIDATE_COUNT="
                    << candidates
                    << "\n";
            }
        }

        if (
            result.modulus_gt_sqrt_p
        ) {
            ++modulus_gt_sqrt_p;

            if (!printed_first_M_gt_sqrt_p) {
                printed_first_M_gt_sqrt_p = true;

                std::cout
                    << "\nFIRST_M_GT_SQRT_P\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << result.N
                    << " S="
                    << result.s
                    << "\n";

                std::cout
                    << pair_to_string(
                        result.best
                    )
                    << "\n";

                std::cout
                    << "SQRT_P="
                    << integer_sqrt(p)
                    << "\n";

                std::cout
                    << "CANDIDATE_COUNT="
                    << candidates
                    << "\n";
            }
        }

        if (candidates == 1) {
            ++candidate_count_one;

            if (!printed_first_one_candidate) {
                printed_first_one_candidate = true;

                std::cout
                    << "\nFIRST_ONE_CANDIDATE\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << result.N
                    << " S="
                    << result.s
                    << "\n";

                std::cout
                    << pair_to_string(
                        result.best
                    )
                    << "\n";

                std::cout
                    << "CANDIDATE_COUNT=1\n";
            }
        }

        if (candidates <= 2) {
            ++candidate_count_le_2;
        }

        if (candidates <= 10) {
            ++candidate_count_le_10;
        }

        if (candidates <= 100) {
            ++candidate_count_le_100;
        }

        if (candidates <= 1000) {
            ++candidate_count_le_1000;
        }

        if (
            !printed_first
        ) {
            printed_first = true;

            std::cout
                << "\nFIRST_VALID_CASE\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " N="
                << result.N
                << " S="
                << result.s
                << "\n";

            std::cout
                << pair_to_string(
                    result.best
                )
                << "\n";

            std::cout
                << "CANDIDATE_COUNT="
                << candidates
                << "\n";

            std::cout
                << "BELOW_FACTOR="
                << result.best
                    .candidates_below_factor
                << "\n";

            std::cout
                << "ABOVE_FACTOR="
                << result.best
                    .candidates_above_factor
                << "\n";
        }

        if (
            M >= 100000 &&
            printed_large_M < 5
        ) {
            ++printed_large_M;

            std::cout
                << "\nLARGE_M_CASE\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << "\n";

            std::cout
                << pair_to_string(
                    result.best
                )
                << "\n";

            std::cout
                << "CANDIDATE_COUNT="
                << candidates
                << "\n";
        }
    }

    if (
        min_M == UINT64_MAX
    ) {
        min_M = 0;
    }

    if (
        min_candidates == UINT64_MAX
    ) {
        min_candidates = 0;
    }

    if (
        min_M_over_p == UINT64_MAX
    ) {
        min_M_over_p = 0;
    }

    const double average_M =
        hits == 0
            ? 0.0
            : static_cast<double>(
                  total_M
              ) /
              static_cast<double>(
                  hits
              );

    const double average_candidates =
        hits == 0
            ? 0.0
            : static_cast<double>(
                  total_candidates
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
        << "VALID_CRT_CASES="
        << hits
        << "\n";

    std::cout
        << "NO_CRT_CASES="
        << (
            cases.size() - hits
        )
        << "\n";

    std::cout
        << "MODULUS_GT_P="
        << modulus_gt_p
        << "\n";

    std::cout
        << "MODULUS_GT_SQRT_P="
        << modulus_gt_sqrt_p
        << "\n";

    std::cout
        << "CANDIDATE_COUNT_1="
        << candidate_count_one
        << "\n";

    std::cout
        << "CANDIDATE_COUNT_LE_2="
        << candidate_count_le_2
        << "\n";

    std::cout
        << "CANDIDATE_COUNT_LE_10="
        << candidate_count_le_10
        << "\n";

    std::cout
        << "CANDIDATE_COUNT_LE_100="
        << candidate_count_le_100
        << "\n";

    std::cout
        << "CANDIDATE_COUNT_LE_1000="
        << candidate_count_le_1000
        << "\n";

    std::cout
        << "MIN_M="
        << min_M
        << "\n";

    std::cout
        << "MAX_M="
        << max_M
        << "\n";

    std::cout
        << "AVERAGE_M="
        << average_M
        << "\n";

    std::cout
        << "MIN_CANDIDATE_COUNT="
        << min_candidates
        << "\n";

    std::cout
        << "MAX_CANDIDATE_COUNT="
        << max_candidates
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_COUNT="
        << average_candidates
        << "\n";

    std::cout
        << "MIN_M_OVER_P_FLOOR="
        << min_M_over_p
        << "\n";

    std::cout
        << "MAX_M_OVER_P_FLOOR="
        << max_M_over_p
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
