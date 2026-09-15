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

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 residue = 0;
    u64 modulus = 0;

    int m1 = 0;
    int m2 = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 candidate_index = 0;
};

struct CaseResult {
    u64 N = 0;
    u64 p = 0;
    u64 q = 0;
    u64 s = 0;

    bool hit = false;

    Hit first_hit;

    u64 crt_classes = 0;
    u64 candidate_tests = 0;
};

static std::vector<u64> generate_primes(int limit) {
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

static u64 integer_sqrt(u64 n) {
    u64 x = 0;

    {
        u64 lo = 0;
        u64 hi = std::min<u64>(
            n,
            1ULL << 32
        );

        while (lo <= hi) {
            const u64 mid =
                lo + (hi - lo) / 2;

            const i128 sq =
                static_cast<i128>(mid) *
                static_cast<i128>(mid);

            if (sq <= static_cast<i128>(n)) {
                x = mid;
                lo = mid + 1;
            } else {
                if (mid == 0) {
                    break;
                }

                hi = mid - 1;
            }
        }
    }

    return x;
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

static CRTPair make_crt_pair(
    int m1,
    u64 k1,
    int m2,
    u64 k2
) {
    CRTPair result;

    result.m1 = m1;
    result.m2 = m2;
    result.k1 = k1;
    result.k2 = k2;

    result.modulus = k1 * k2;

    /*
        r ≡ k2 (mod k1)
        r ≡ k1 (mod k2)

        Since determinant = 1 implies
        gcd(k1,k2)=1, CRT gives

        r ≡ k1 + k2 (mod k1*k2).
    */

    result.residue =
        (k1 + k2) %
        result.modulus;

    return result;
}

static std::vector<CRTPair>
generate_crt_pairs(
    const std::vector<int>& multipliers,
    u64 K_LIMIT
) {
    std::vector<CRTPair> pairs;

    for (int m1 : multipliers) {
        for (int m2 : multipliers) {
            for (u64 k1 = 2; k1 <= K_LIMIT; ++k1) {
                /*
                    Solve

                    m2*k1 - m1*k2 = 1

                    for k2.

                    This avoids a full k1*k2 scan.
                */

                const i128 numerator =
                    static_cast<i128>(m2) *
                    static_cast<i128>(k1) -
                    static_cast<i128>(1);

                if (numerator <= 0) {
                    continue;
                }

                const i128 denominator =
                    static_cast<i128>(m1);

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
                        static_cast<i128>(K_LIMIT)
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

                pairs.push_back(
                    make_crt_pair(
                        m1,
                        k1,
                        m2,
                        k2
                    )
                );
            }
        }
    }

    return pairs;
}

static std::string crt_pair_string(
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
        " modulus=" +
        std::to_string(pair.modulus);

    return out;
}

static Hit search_crt_classes(
    u64 N,
    u64 s,
    const std::vector<CRTPair>& pairs,
    u64& candidate_tests,
    u64& class_count
) {
    Hit first_hit;

    /*
        Deduplicate by (residue, modulus).
        Different multiplier pairs can produce
        the same CRT class.
    */

    std::set<std::pair<u64, u64>> seen;

    for (const CRTPair& pair : pairs) {
        const auto key =
            std::make_pair(
                pair.residue,
                pair.modulus
            );

        if (!seen.insert(key).second) {
            continue;
        }

        ++class_count;

        const u64 residue =
            pair.residue;

        const u64 modulus =
            pair.modulus;

        if (modulus == 0) {
            continue;
        }

        if (residue > s) {
            continue;
        }

        /*
            All positive candidates <= s:

                residue + j*modulus

            Special case residue = 0.
        */

        u64 first = residue;

        if (first == 0) {
            first = modulus;
        }

        if (first > s) {
            continue;
        }

        const u64 count =
            1 +
            (s - first) / modulus;

        for (u64 j = 0; j < count; ++j) {
            const u64 candidate =
                first + j * modulus;

            if (candidate < 2) {
                continue;
            }

            ++candidate_tests;

            const u64 g =
                std::gcd(
                    N,
                    candidate
                );

            if (
                g != 1 &&
                g != N
            ) {
                first_hit.valid = true;
                first_hit.factor = g;
                first_hit.residue = residue;
                first_hit.modulus = modulus;

                first_hit.m1 = pair.m1;
                first_hit.m2 = pair.m2;

                first_hit.k1 = pair.k1;
                first_hit.k2 = pair.k2;

                first_hit.candidate_index =
                    j;

                return first_hit;
            }
        }
    }

    return first_hit;
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

    result.first_hit =
        search_crt_classes(
            result.N,
            result.s,
            pairs,
            result.candidate_tests,
            result.crt_classes
        );

    result.hit =
        result.first_hit.valid;

    return result;
}

static std::vector<std::pair<u64, u64>>
generate_cases(
    const std::vector<u64>& primes,
    std::size_t count,
    std::uint64_t seed
) {
    std::mt19937_64 rng(seed);

    std::vector<
        std::pair<u64, u64>
    > cases;

    if (primes.size() < 2) {
        return cases;
    }

    cases.reserve(count);

    for (
        std::size_t i = 0;
        i < count;
        ++i
    ) {
        std::size_t a =
            static_cast<std::size_t>(
                rng() % primes.size()
            );

        std::size_t b =
            static_cast<std::size_t>(
                rng() % primes.size()
            );

        while (a == b) {
            b =
                static_cast<std::size_t>(
                    rng() % primes.size()
                );
        }

        if (primes[a] > primes[b]) {
            std::swap(a, b);
        }

        cases.emplace_back(
            primes[a],
            primes[b]
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
            << (index + 1)
            << "/"
            << total
            << "\n";
    }
}

static bool hit_is_true_factor(
    const Hit& hit,
    u64 p,
    u64 q
) {
    if (!hit.valid) {
        return false;
    }

    return
        hit.factor == p ||
        hit.factor == q;
}

int main() {
    constexpr int EXPERIMENT = 450;

    constexpr int PRIME_LIMIT = 100000;

    constexpr std::size_t CASE_COUNT = 5000;

    constexpr u64 K_LIMIT = 500;

    constexpr std::uint64_t SEED =
        0x450450450ULL;

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

    for (
        std::size_t i = 0;
        i < std::min<std::size_t>(
            pairs.size(),
            20
        );
        ++i
    ) {
        std::cout
            << "CRT_PAIR "
            << i
            << " "
            << crt_pair_string(
                pairs[i]
            )
            << "\n";
    }

    const auto cases =
        generate_cases(
            primes,
            CASE_COUNT,
            SEED
        );

    std::size_t hits = 0;
    std::size_t misses = 0;
    std::size_t true_factor_hits = 0;

    u64 total_candidate_tests = 0;
    u64 max_candidate_tests = 0;

    u64 total_classes = 0;
    u64 max_classes = 0;

    u64 min_candidate_tests = UINT64_MAX;

    bool printed_first_hit = false;
    bool printed_first_miss = false;

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

        total_candidate_tests +=
            result.candidate_tests;

        total_classes +=
            result.crt_classes;

        max_candidate_tests =
            std::max(
                max_candidate_tests,
                result.candidate_tests
            );

        max_classes =
            std::max(
                max_classes,
                result.crt_classes
            );

        min_candidate_tests =
            std::min(
                min_candidate_tests,
                result.candidate_tests
            );

        if (result.hit) {
            ++hits;

            if (
                hit_is_true_factor(
                    result.first_hit,
                    p,
                    q
                )
            ) {
                ++true_factor_hits;
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
                    << result.N
                    << " S="
                    << result.s
                    << "\n";

                std::cout
                    << "HIT_FACTOR="
                    << result.first_hit.factor
                    << "\n";

                std::cout
                    << "HIT_RESIDUE="
                    << result.first_hit.residue
                    << "\n";

                std::cout
                    << "HIT_MODULUS="
                    << result.first_hit.modulus
                    << "\n";

                std::cout
                    << "HIT_M1="
                    << result.first_hit.m1
                    << "\n";

                std::cout
                    << "HIT_M2="
                    << result.first_hit.m2
                    << "\n";

                std::cout
                    << "HIT_K1="
                    << result.first_hit.k1
                    << "\n";

                std::cout
                    << "HIT_K2="
                    << result.first_hit.k2
                    << "\n";

                std::cout
                    << "HIT_CANDIDATE_INDEX="
                    << result.first_hit.candidate_index
                    << "\n";

                std::cout
                    << "CANDIDATE_TESTS="
                    << result.candidate_tests
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
                    << result.N
                    << " S="
                    << result.s
                    << "\n";

                std::cout
                    << "CRT_CLASSES="
                    << result.crt_classes
                    << "\n";

                std::cout
                    << "CANDIDATE_TESTS="
                    << result.candidate_tests
                    << "\n";
            }
        }
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

    const double average_classes =
        cases.empty()
            ? 0.0
            : static_cast<double>(
                  total_classes
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
        << "TRUE_FACTOR_HITS="
        << true_factor_hits
        << "\n";

    std::cout
        << "DETERMINANT_ONE_PAIRS="
        << pairs.size()
        << "\n";

    std::cout
        << "TOTAL_CRT_CLASSES_TESTED="
        << total_classes
        << "\n";

    std::cout
        << "AVERAGE_CRT_CLASSES_PER_CASE="
        << average_classes
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
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
