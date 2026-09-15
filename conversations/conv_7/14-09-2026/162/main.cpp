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

struct Witness {
    bool valid = false;

    u64 factor = 0;

    u64 j = 0;
    u64 modulus = 0;
    u64 residue = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    int m1 = 0;
    int m2 = 0;

    u64 candidate_index = 0;
};

struct CaseResult {
    u64 p = 0;
    u64 q = 0;
    u64 N = 0;
    u64 s = 0;

    Witness smallest_j;
    Witness smallest_modulus;
    Witness smallest_candidate_index;
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

    return pairs;
}

static bool better_j(
    const Witness& a,
    const Witness& b
) {
    if (!a.valid) {
        return false;
    }

    if (!b.valid) {
        return true;
    }

    if (a.j != b.j) {
        return a.j < b.j;
    }

    if (a.modulus != b.modulus) {
        return a.modulus > b.modulus;
    }

    if (a.k1 != b.k1) {
        return a.k1 > b.k1;
    }

    if (a.k2 != b.k2) {
        return a.k2 > b.k2;
    }

    return a.m1 < b.m1;
}

static bool better_modulus(
    const Witness& a,
    const Witness& b
) {
    if (!a.valid) {
        return false;
    }

    if (!b.valid) {
        return true;
    }

    if (a.modulus != b.modulus) {
        return a.modulus > b.modulus;
    }

    if (a.j != b.j) {
        return a.j < b.j;
    }

    return a.k1 > b.k1;
}

static bool better_candidate_index(
    const Witness& a,
    const Witness& b
) {
    if (!a.valid) {
        return false;
    }

    if (!b.valid) {
        return true;
    }

    if (a.candidate_index != b.candidate_index) {
        return a.candidate_index <
               b.candidate_index;
    }

    if (a.j != b.j) {
        return a.j < b.j;
    }

    return a.modulus > b.modulus;
}

static Witness make_witness(
    u64 factor,
    const CRTPair& pair
) {
    Witness result;

    const u64 M =
        pair.modulus;

    if (M == 0) {
        return result;
    }

    if (
        factor < pair.residue
    ) {
        return result;
    }

    const u64 difference =
        factor -
        pair.residue;

    if (
        difference % M != 0
    ) {
        return result;
    }

    const u64 j =
        difference / M;

    /*
        factor = k1 + k2 + j*k1*k2
    */
    const i128 reconstructed =
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
        reconstructed !=
        static_cast<i128>(factor)
    ) {
        return result;
    }

    /*
        The CRT class is

            factor = residue + j*M

        and residue = k1+k2 here because
        k1+k2 < k1*k2 for k1,k2 >= 2,
        so candidate_index equals j.
    */

    result.valid = true;

    result.factor = factor;

    result.j = j;

    result.modulus =
        pair.modulus;

    result.residue =
        pair.residue;

    result.k1 =
        pair.k1;

    result.k2 =
        pair.k2;

    result.m1 =
        pair.m1;

    result.m2 =
        pair.m2;

    result.candidate_index =
        j;

    return result;
}

static bool verify_factorization_identity(
    const Witness& w
) {
    if (!w.valid) {
        return false;
    }

    const i128 left =
        (
            static_cast<i128>(w.j) *
            static_cast<i128>(w.k1) +
            static_cast<i128>(1)
        ) *
        (
            static_cast<i128>(w.j) *
            static_cast<i128>(w.k2) +
            static_cast<i128>(1)
        );

    const i128 right =
        static_cast<i128>(w.j) *
        static_cast<i128>(w.factor) +
        static_cast<i128>(1);

    return left == right;
}

static std::string witness_to_string(
    const Witness& w
) {
    if (!w.valid) {
        return "NONE";
    }

    std::string out;

    out +=
        "m1=" +
        std::to_string(w.m1);

    out +=
        " m2=" +
        std::to_string(w.m2);

    out +=
        " k1=" +
        std::to_string(w.k1);

    out +=
        " k2=" +
        std::to_string(w.k2);

    out +=
        " j=" +
        std::to_string(w.j);

    out +=
        " residue=" +
        std::to_string(w.residue);

    out +=
        " M=" +
        std::to_string(w.modulus);

    out +=
        " candidate_index=" +
        std::to_string(
            w.candidate_index
        );

    return out;
}

static CaseResult analyze_case(
    u64 p,
    u64 q,
    const std::vector<CRTPair>& pairs
) {
    CaseResult result;

    result.p = p;
    result.q = q;

    result.N =
        p * q;

    result.s =
        integer_sqrt(
            result.N
        );

    for (const CRTPair& pair : pairs) {
        const Witness w =
            make_witness(
                p,
                pair
            );

        if (!w.valid) {
            continue;
        }

        if (
            better_j(
                w,
                result.smallest_j
            )
        ) {
            result.smallest_j = w;
        }

        if (
            better_modulus(
                w,
                result.smallest_modulus
            )
        ) {
            result.smallest_modulus = w;
        }

        if (
            better_candidate_index(
                w,
                result.smallest_candidate_index
            )
        ) {
            result.smallest_candidate_index = w;
        }
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
    constexpr int EXPERIMENT = 454;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 2000;

    constexpr u64 K_LIMIT = 1000;

    constexpr std::uint64_t SEED =
        0x454454454ULL;

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

    const auto cases =
        generate_cases(
            primes,
            PRIME_MIN,
            CASE_COUNT,
            SEED
        );

    std::size_t valid_cases = 0;

    std::size_t identity_failures = 0;

    std::size_t j_zero = 0;
    std::size_t j_one = 0;
    std::size_t j_le_2 = 0;
    std::size_t j_le_5 = 0;
    std::size_t j_le_10 = 0;
    std::size_t j_le_20 = 0;
    std::size_t j_le_50 = 0;
    std::size_t j_le_100 = 0;
    std::size_t j_gt_100 = 0;
    std::size_t j_gt_1000 = 0;

    u64 min_j = UINT64_MAX;
    u64 max_j = 0;

    u64 total_j = 0;

    u64 min_M_at_min_j =
        UINT64_MAX;

    u64 max_M_at_min_j = 0;

    u64 min_candidate_index =
        UINT64_MAX;

    u64 max_candidate_index = 0;

    u64 total_candidate_index = 0;

    u64 min_modulus =
        UINT64_MAX;

    u64 max_modulus = 0;

    u64 total_modulus = 0;

    std::size_t
        printed_first_small_j = 0;

    bool printed_first_case = false;

    bool printed_first_j_gt_100 =
        false;

    bool printed_first_j_one =
        false;

    bool printed_first_identity =
        false;

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
            analyze_case(
                p,
                q,
                pairs
            );

        print_progress(
            i,
            cases.size()
        );

        if (
            !result.smallest_j.valid
        ) {
            continue;
        }

        ++valid_cases;

        const Witness& w =
            result.smallest_j;

        if (
            !verify_factorization_identity(w)
        ) {
            ++identity_failures;
        }

        min_j =
            std::min(
                min_j,
                w.j
            );

        max_j =
            std::max(
                max_j,
                w.j
            );

        total_j += w.j;

        min_M_at_min_j =
            std::min(
                min_M_at_min_j,
                w.modulus
            );

        max_M_at_min_j =
            std::max(
                max_M_at_min_j,
                w.modulus
            );

        min_candidate_index =
            std::min(
                min_candidate_index,
                w.candidate_index
            );

        max_candidate_index =
            std::max(
                max_candidate_index,
                w.candidate_index
            );

        total_candidate_index +=
            w.candidate_index;

        min_modulus =
            std::min(
                min_modulus,
                result.smallest_modulus.modulus
            );

        max_modulus =
            std::max(
                max_modulus,
                result.smallest_modulus.modulus
            );

        total_modulus +=
            result.smallest_modulus.modulus;

        if (w.j == 0) {
            ++j_zero;
        }

        if (w.j == 1) {
            ++j_one;

            if (!printed_first_j_one) {
                printed_first_j_one = true;

                std::cout
                    << "\nFIRST_J_ONE\n";

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
                    << "SMALLEST_J="
                    << witness_to_string(w)
                    << "\n";

                std::cout
                    << "IDENTITY="
                    << (
                        verify_factorization_identity(
                            w
                        )
                            ? "PASS"
                            : "FAIL"
                    )
                    << "\n";
            }
        }

        if (w.j <= 2) {
            ++j_le_2;
        }

        if (w.j <= 5) {
            ++j_le_5;
        }

        if (w.j <= 10) {
            ++j_le_10;
        }

        if (w.j <= 20) {
            ++j_le_20;
        }

        if (w.j <= 50) {
            ++j_le_50;
        }

        if (w.j <= 100) {
            ++j_le_100;
        }

        if (w.j > 100) {
            ++j_gt_100;
        }

        if (w.j > 1000) {
            ++j_gt_1000;

            if (!printed_first_j_gt_100) {
                printed_first_j_gt_100 = true;

                std::cout
                    << "\nFIRST_J_GT_1000\n";

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
                    << "SMALLEST_J="
                    << witness_to_string(w)
                    << "\n";
            }
        }

        if (
            w.j <= 10 &&
            printed_first_small_j < 5
        ) {
            ++printed_first_small_j;

            std::cout
                << "\nSMALL_J_CASE\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << "\n";

            std::cout
                << "SMALLEST_J="
                << witness_to_string(w)
                << "\n";

            std::cout
                << "IDENTITY="
                << (
                    verify_factorization_identity(w)
                        ? "PASS"
                        : "FAIL"
                )
                << "\n";
        }

        if (!printed_first_case) {
            printed_first_case = true;

            std::cout
                << "\nFIRST_VALID_CASE\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " N="
                << result.N
                << "\n";

            std::cout
                << "SMALLEST_J="
                << witness_to_string(
                    result.smallest_j
                )
                << "\n";

            std::cout
                << "SMALLEST_MODULUS="
                << witness_to_string(
                    result.smallest_modulus
                )
                << "\n";

            std::cout
                << "SMALLEST_CANDIDATE_INDEX="
                << witness_to_string(
                    result.smallest_candidate_index
                )
                << "\n";
        }
    }

    if (
        valid_cases == 0
    ) {
        min_j = 0;
        min_M_at_min_j = 0;
        min_candidate_index = 0;
        min_modulus = 0;
    }

    if (
        min_j == UINT64_MAX
    ) {
        min_j = 0;
    }

    if (
        min_M_at_min_j == UINT64_MAX
    ) {
        min_M_at_min_j = 0;
    }

    if (
        min_candidate_index == UINT64_MAX
    ) {
        min_candidate_index = 0;
    }

    if (
        min_modulus == UINT64_MAX
    ) {
        min_modulus = 0;
    }

    const double average_j =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                  total_j
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double average_candidate_index =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                  total_candidate_index
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double average_max_modulus =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                  total_modulus
              ) /
              static_cast<double>(
                  valid_cases
              );

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "VALID_CASES="
        << valid_cases
        << "\n";

    std::cout
        << "INVALID_CASES="
        << (
            cases.size() -
            valid_cases
        )
        << "\n";

    std::cout
        << "IDENTITY_FAILURES="
        << identity_failures
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
        << "J_GT_100="
        << j_gt_100
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
        << "MIN_CANDIDATE_INDEX="
        << min_candidate_index
        << "\n";

    std::cout
        << "MAX_CANDIDATE_INDEX="
        << max_candidate_index
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_INDEX="
        << average_candidate_index
        << "\n";

    std::cout
        << "MIN_MODULUS_AT_MIN_J="
        << min_M_at_min_j
        << "\n";

    std::cout
        << "MAX_MODULUS_AT_MIN_J="
        << max_M_at_min_j
        << "\n";

    std::cout
        << "MIN_MODULUS_OVERALL="
        << min_modulus
        << "\n";

    std::cout
        << "MAX_MODULUS_OVERALL="
        << max_modulus
        << "\n";

    std::cout
        << "AVERAGE_MAX_MODULUS="
        << average_max_modulus
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
