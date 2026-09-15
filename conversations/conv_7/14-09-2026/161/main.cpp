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

    u64 residue = 0;
    u64 modulus = 0;
};

struct SearchHit {
    bool valid = false;

    u64 factor = 0;

    u64 candidate_tests = 0;
    u64 crt_classes_tested = 0;

    u64 pair_index = 0;
    u64 candidate_index = 0;

    CRTPair pair;
};

struct CaseResult {
    u64 N = 0;
    u64 s = 0;

    SearchHit hit;
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

                const u64 residue =
                    (k1 + k2) %
                    modulus;

                CRTPair pair;

                pair.m1 = m1;
                pair.m2 = m2;
                pair.k1 = k1;
                pair.k2 = k2;
                pair.residue = residue;
                pair.modulus = modulus;

                pairs.push_back(pair);
            }
        }
    }

    /*
        Sort strongest CRT classes first.

        Larger M means fewer integers remain in
        the residue class below sqrt(N).
    */
    std::sort(
        pairs.begin(),
        pairs.end(),
        [](
            const CRTPair& a,
            const CRTPair& b
        ) {
            if (a.modulus != b.modulus) {
                return a.modulus >
                       b.modulus;
            }

            if (a.k1 != b.k1) {
                return a.k1 >
                       b.k1;
            }

            if (a.k2 != b.k2) {
                return a.k2 >
                       b.k2;
            }

            if (a.m1 != b.m1) {
                return a.m1 <
                       b.m1;
            }

            return a.m2 <
                   b.m2;
        }
    );

    /*
        Different determinant-one pairs can produce
        exactly the same CRT class. Keep one copy.
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
            unique_pairs.push_back(
                pair
            );
        }
    }

    return unique_pairs;
}

static SearchHit blind_search(
    u64 N,
    u64 s,
    const std::vector<CRTPair>& pairs
) {
    SearchHit hit;

    for (
        std::size_t pair_index = 0;
        pair_index < pairs.size();
        ++pair_index
    ) {
        const CRTPair& pair =
            pairs[pair_index];

        ++hit.crt_classes_tested;

        const u64 M =
            pair.modulus;

        if (M == 0) {
            continue;
        }

        /*
            Positive solutions:

                x = residue + j*M

            with x <= sqrt(N).

            residue = 0 means the first positive
            representative is M.
        */
        u64 first =
            pair.residue;

        if (first == 0) {
            first = M;
        }

        if (first > s) {
            continue;
        }

        const u64 count =
            1 +
            (s - first) / M;

        for (
            u64 j = 0;
            j < count;
            ++j
        ) {
            const u64 candidate =
                first +
                j * M;

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

                hit.pair_index =
                    static_cast<u64>(
                        pair_index
                    );

                hit.candidate_index = j;

                hit.pair = pair;

                return hit;
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
        " residue=" +
        std::to_string(pair.residue);

    out +=
        " M=" +
        std::to_string(pair.modulus);

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
    constexpr int EXPERIMENT = 453;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 2000;

    constexpr u64 K_LIMIT = 1000;

    constexpr std::uint64_t SEED =
        0x453453453ULL;

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

    std::vector<CRTPair> pairs =
        generate_crt_pairs(
            MULTIPLIERS,
            K_LIMIT
        );

    std::cout
        << "UNIQUE_CRT_CLASSES="
        << pairs.size()
        << "\n";

    for (
        std::size_t i = 0;
        i < std::min<std::size_t>(
            20,
            pairs.size()
        );
        ++i
    ) {
        std::cout
            << "CRT_RANK "
            << i
            << " "
            << pair_to_string(
                pairs[i]
            )
            << "\n";
    }

    const auto cases =
        generate_cases(
            primes,
            PRIME_MIN,
            CASE_COUNT,
            SEED
        );

    std::size_t hits = 0;
    std::size_t misses = 0;

    u64 total_candidate_tests = 0;
    u64 total_classes_tested = 0;

    u64 min_candidate_tests =
        UINT64_MAX;

    u64 max_candidate_tests = 0;

    u64 min_classes_tested =
        UINT64_MAX;

    u64 max_classes_tested = 0;

    u64 min_hit_rank =
        UINT64_MAX;

    u64 max_hit_rank = 0;

    u64 total_hit_rank = 0;

    u64 hits_with_M_gt_s = 0;
    u64 hits_with_M_gt_factor = 0;
    u64 hits_with_one_candidate_class = 0;
    u64 hits_with_le_10_candidates = 0;
    u64 hits_with_le_100_candidates = 0;
    u64 hits_with_le_1000_candidates = 0;

    u64 first_factor_overshoot_examples = 0;

    bool printed_first_hit = false;
    bool printed_first_miss = false;
    bool printed_first_M_gt_s = false;
    bool printed_first_large_candidate_count = false;

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

        const SearchHit hit =
            blind_search(
                N,
                s,
                pairs
            );

        print_progress(
            i,
            cases.size()
        );

        total_candidate_tests +=
            hit.candidate_tests;

        total_classes_tested +=
            hit.crt_classes_tested;

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

        min_classes_tested =
            std::min(
                min_classes_tested,
                hit.crt_classes_tested
            );

        max_classes_tested =
            std::max(
                max_classes_tested,
                hit.crt_classes_tested
            );

        if (hit.valid) {
            ++hits;

            const u64 rank =
                hit.pair_index + 1;

            min_hit_rank =
                std::min(
                    min_hit_rank,
                    rank
                );

            max_hit_rank =
                std::max(
                    max_hit_rank,
                    rank
                );

            total_hit_rank += rank;

            if (
                hit.pair.modulus > s
            ) {
                ++hits_with_M_gt_s;

                if (!printed_first_M_gt_s) {
                    printed_first_M_gt_s = true;

                    std::cout
                        << "\nFIRST_HIT_M_GT_S\n";

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
                        << pair_to_string(
                            hit.pair
                        )
                        << "\n";

                    std::cout
                        << "HIT_FACTOR="
                        << hit.factor
                        << "\n";

                    std::cout
                        << "PAIR_RANK="
                        << rank
                        << "\n";

                    std::cout
                        << "CANDIDATE_INDEX="
                        << hit.candidate_index
                        << "\n";

                    std::cout
                        << "CANDIDATE_TESTS="
                        << hit.candidate_tests
                        << "\n";

                    std::cout
                        << "CLASSES_TESTED="
                        << hit.crt_classes_tested
                        << "\n";
                }
            }

            /*
                This is mainly diagnostic. For the intended
                lower factor p <= sqrt(N), M > p means the
                residue class is so sparse that at most one
                relevant integer lies below p, though M > p
                is not necessary for a successful blind hit.
            */

            if (
                hit.pair.modulus > p
            ) {
                ++hits_with_M_gt_factor;
            }

            if (
                hit.candidate_index == 0
            ) {
                ++hits_with_one_candidate_class;
            }

            if (
                hit.candidate_tests <= 10
            ) {
                ++hits_with_le_10_candidates;
            }

            if (
                hit.candidate_tests <= 100
            ) {
                ++hits_with_le_100_candidates;
            }

            if (
                hit.candidate_tests <= 1000
            ) {
                ++hits_with_le_1000_candidates;
            }

            if (
                !printed_first_hit
            ) {
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
                    << pair_to_string(
                        hit.pair
                    )
                    << "\n";

                std::cout
                    << "HIT_FACTOR="
                    << hit.factor
                    << "\n";

                std::cout
                    << "PAIR_RANK="
                    << rank
                    << "\n";

                std::cout
                    << "CANDIDATE_INDEX="
                    << hit.candidate_index
                    << "\n";

                std::cout
                    << "CANDIDATE_TESTS="
                    << hit.candidate_tests
                    << "\n";

                std::cout
                    << "CLASSES_TESTED="
                    << hit.crt_classes_tested
                    << "\n";
            }

            if (
                hit.candidate_tests > 10000 &&
                !printed_first_large_candidate_count
            ) {
                printed_first_large_candidate_count = true;

                std::cout
                    << "\nFIRST_LARGE_SEARCH\n";

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
                    << pair_to_string(
                        hit.pair
                    )
                    << "\n";

                std::cout
                    << "HIT_FACTOR="
                    << hit.factor
                    << "\n";

                std::cout
                    << "PAIR_RANK="
                    << rank
                    << "\n";

                std::cout
                    << "CANDIDATE_TESTS="
                    << hit.candidate_tests
                    << "\n";

                std::cout
                    << "CLASSES_TESTED="
                    << hit.crt_classes_tested
                    << "\n";
            }

            if (
                hit.pair.modulus > s &&
                first_factor_overshoot_examples < 5
            ) {
                ++first_factor_overshoot_examples;

                std::cout
                    << "\nM_GT_S_CASE\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << "\n";

                std::cout
                    << "S="
                    << s
                    << "\n";

                std::cout
                    << pair_to_string(
                        hit.pair
                    )
                    << "\n";

                std::cout
                    << "HIT_FACTOR="
                    << hit.factor
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
                    << "CLASSES_TESTED="
                    << hit.crt_classes_tested
                    << "\n";

                std::cout
                    << "CANDIDATE_TESTS="
                    << hit.candidate_tests
                    << "\n";
            }
        }
    }

    if (
        min_candidate_tests ==
        UINT64_MAX
    ) {
        min_candidate_tests = 0;
    }

    if (
        min_classes_tested ==
        UINT64_MAX
    ) {
        min_classes_tested = 0;
    }

    if (
        min_hit_rank ==
        UINT64_MAX
    ) {
        min_hit_rank = 0;
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

    const double average_classes_tested =
        cases.empty()
            ? 0.0
            : static_cast<double>(
                  total_classes_tested
              ) /
              static_cast<double>(
                  cases.size()
              );

    const double average_hit_rank =
        hits == 0
            ? 0.0
            : static_cast<double>(
                  total_hit_rank
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
        << "TOTAL_UNIQUE_CRT_CLASSES="
        << pairs.size()
        << "\n";

    std::cout
        << "TOTAL_CRT_CLASSES_TESTED="
        << total_classes_tested
        << "\n";

    std::cout
        << "AVERAGE_CRT_CLASSES_TESTED="
        << average_classes_tested
        << "\n";

    std::cout
        << "MIN_CRT_CLASSES_TESTED="
        << min_classes_tested
        << "\n";

    std::cout
        << "MAX_CRT_CLASSES_TESTED="
        << max_classes_tested
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
        << "MIN_CANDIDATE_TESTS="
        << min_candidate_tests
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_candidate_tests
        << "\n";

    std::cout
        << "MIN_HIT_RANK="
        << min_hit_rank
        << "\n";

    std::cout
        << "MAX_HIT_RANK="
        << max_hit_rank
        << "\n";

    std::cout
        << "AVERAGE_HIT_RANK="
        << average_hit_rank
        << "\n";

    std::cout
        << "HITS_M_GT_S="
        << hits_with_M_gt_s
        << "\n";

    std::cout
        << "HITS_M_GT_FACTOR="
        << hits_with_M_gt_factor
        << "\n";

    std::cout
        << "HITS_FIRST_CANDIDATE="
        << hits_with_one_candidate_class
        << "\n";

    std::cout
        << "HITS_CANDIDATE_TESTS_LE_10="
        << hits_with_le_10_candidates
        << "\n";

    std::cout
        << "HITS_CANDIDATE_TESTS_LE_100="
        << hits_with_le_100_candidates
        << "\n";

    std::cout
        << "HITS_CANDIDATE_TESTS_LE_1000="
        << hits_with_le_1000_candidates
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
