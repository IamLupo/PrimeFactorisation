#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
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

struct Witness {
    bool valid = false;

    u64 j = 0;

    int m1 = 0;
    int m2 = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 modulus = 0;
    u64 residue = 0;
};

struct CaseResult {
    u64 p = 0;
    u64 q = 0;
    u64 N = 0;

    Witness best;

    u64 pair_count = 0;
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

static std::vector<CRTPair>
generate_pairs(
    int M_LIMIT,
    u64 K_LIMIT
) {
    std::vector<CRTPair> pairs;

    for (
        int m1 = 1;
        m1 <= M_LIMIT;
        ++m1
    ) {
        for (
            int m2 = 1;
            m2 <= M_LIMIT;
            ++m2
        ) {
            /*
                m2*k1 - m1*k2 = 1

                Therefore

                    k2 = (m2*k1 - 1)/m1.

                We only need to scan k1.
            */

            for (
                u64 k1 = 2;
                k1 <= K_LIMIT;
                ++k1
            ) {
                const i128 numerator =
                    static_cast<i128>(m2) *
                    static_cast<i128>(k1) -
                    static_cast<i128>(1);

                if (numerator <= 0) {
                    continue;
                }

                if (
                    numerator %
                        static_cast<i128>(m1) !=
                    0
                ) {
                    continue;
                }

                const i128 k2_i =
                    numerator /
                    static_cast<i128>(m1);

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
                    std::gcd(k1, k2) != 1
                ) {
                    continue;
                }

                CRTPair pair;

                pair.m1 = m1;
                pair.m2 = m2;

                pair.k1 = k1;
                pair.k2 = k2;

                pair.modulus =
                    k1 * k2;

                pair.residue =
                    k1 + k2;

                pairs.push_back(pair);
            }
        }
    }

    return pairs;
}

static bool better_witness(
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

    return a.k2 > b.k2;
}

static Witness find_smallest_j(
    u64 factor,
    const std::vector<CRTPair>& pairs
) {
    Witness best;

    for (
        const CRTPair& pair :
        pairs
    ) {
        /*
            factor =
                k1 + k2 + j*k1*k2

            Hence

                j =
                    (factor-k1-k2)/(k1*k2).
        */

        const i128 difference =
            static_cast<i128>(factor) -
            static_cast<i128>(pair.k1) -
            static_cast<i128>(pair.k2);

        if (difference < 0) {
            continue;
        }

        const i128 M =
            static_cast<i128>(
                pair.k1
            ) *
            static_cast<i128>(
                pair.k2
            );

        if (M <= 0) {
            continue;
        }

        if (
            difference % M != 0
        ) {
            continue;
        }

        const i128 j_i =
            difference / M;

        if (
            j_i < 0 ||
            j_i >
                static_cast<i128>(
                    UINT64_MAX
                )
        ) {
            continue;
        }

        const u64 j =
            static_cast<u64>(j_i);

        Witness witness;

        witness.valid = true;
        witness.j = j;

        witness.m1 = pair.m1;
        witness.m2 = pair.m2;

        witness.k1 = pair.k1;
        witness.k2 = pair.k2;

        witness.modulus =
            pair.modulus;

        witness.residue =
            pair.residue;

        if (
            better_witness(
                witness,
                best
            )
        ) {
            best = witness;
        }
    }

    return best;
}

static bool verify_identity(
    u64 factor,
    const Witness& witness
) {
    if (!witness.valid) {
        return false;
    }

    const i128 lhs =
        static_cast<i128>(
            witness.k1
        ) +
        static_cast<i128>(
            witness.k2
        ) +
        static_cast<i128>(
            witness.j
        ) *
        static_cast<i128>(
            witness.k1
        ) *
        static_cast<i128>(
            witness.k2
        );

    if (
        lhs !=
        static_cast<i128>(factor)
    ) {
        return false;
    }

    const i128 left2 =
        (
            static_cast<i128>(
                witness.j
            ) *
            static_cast<i128>(
                witness.k1
            ) +
            static_cast<i128>(1)
        ) *
        (
            static_cast<i128>(
                witness.j
            ) *
            static_cast<i128>(
                witness.k2
            ) +
            static_cast<i128>(1)
        );

    const i128 right2 =
        static_cast<i128>(
            witness.j
        ) *
        static_cast<i128>(
            factor
        ) +
        static_cast<i128>(1);

    return left2 == right2;
}

static CaseResult analyze_case(
    u64 p,
    u64 q,
    const std::vector<CRTPair>& pairs
) {
    CaseResult result;

    result.p = p;
    result.q = q;
    result.N = p * q;

    result.pair_count =
        pairs.size();

    /*
        We deliberately analyze the smaller
        factor because that is the factor a
        sqrt(N)-bounded blind search can reach.
    */

    result.best =
        find_smallest_j(
            p,
            pairs
        );

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

static std::string witness_to_string(
    const Witness& witness
) {
    if (!witness.valid) {
        return "NONE";
    }

    std::string out;

    out +=
        "j=" +
        std::to_string(
            witness.j
        );

    out +=
        " m1=" +
        std::to_string(
            witness.m1
        );

    out +=
        " m2=" +
        std::to_string(
            witness.m2
        );

    out +=
        " k1=" +
        std::to_string(
            witness.k1
        );

    out +=
        " k2=" +
        std::to_string(
            witness.k2
        );

    out +=
        " M=" +
        std::to_string(
            witness.modulus
        );

    return out;
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
    constexpr int EXPERIMENT = 457;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 5000;

    constexpr int M_LIMIT = 7;

    /*
        This is deliberately much larger than
        the previously observed maximum j=28.
    */
    constexpr u64 K_LIMIT = 1000;

    constexpr std::uint64_t SEED =
        0x457457457ULL;

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
        << "M_LIMIT="
        << M_LIMIT
        << "\n";

    std::cout
        << "K_LIMIT="
        << K_LIMIT
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
        generate_pairs(
            M_LIMIT,
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
    std::size_t invalid_cases = 0;

    std::size_t identity_failures = 0;

    std::size_t j_one = 0;
    std::size_t j_le_2 = 0;
    std::size_t j_le_3 = 0;
    std::size_t j_le_5 = 0;
    std::size_t j_le_10 = 0;
    std::size_t j_le_20 = 0;
    std::size_t j_le_28 = 0;
    std::size_t j_le_50 = 0;
    std::size_t j_le_100 = 0;
    std::size_t j_le_500 = 0;
    std::size_t j_le_1000 = 0;
    std::size_t j_gt_1000 = 0;

    u64 min_j = UINT64_MAX;
    u64 max_j = 0;
    u64 total_j = 0;

    u64 max_j_p = 0;
    u64 max_j_q = 0;

    std::size_t printed_small_examples = 0;

    bool printed_first_case = false;
    bool printed_first_large_j = false;
    bool printed_first_invalid = false;

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

        if (!result.best.valid) {
            ++invalid_cases;

            if (!printed_first_invalid) {
                printed_first_invalid = true;

                std::cout
                    << "\nFIRST_INVALID_CASE\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << result.N
                    << "\n";
            }

            continue;
        }

        ++valid_cases;

        const Witness& w =
            result.best;

        if (
            !verify_identity(
                p,
                w
            )
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

        max_j_p =
            std::max(
                max_j_p,
                w.j
            );

        if (w.j == 1) {
            ++j_one;
        }

        if (w.j <= 2) {
            ++j_le_2;
        }

        if (w.j <= 3) {
            ++j_le_3;
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

        if (w.j <= 28) {
            ++j_le_28;
        }

        if (w.j <= 50) {
            ++j_le_50;
        }

        if (w.j <= 100) {
            ++j_le_100;
        }

        if (w.j <= 500) {
            ++j_le_500;
        }

        if (w.j <= 1000) {
            ++j_le_1000;
        }

        if (w.j > 1000) {
            ++j_gt_1000;
        }

        if (
            !printed_first_case
        ) {
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
                    w
                )
                << "\n";

            std::cout
                << "IDENTITY="
                << (
                    verify_identity(
                        p,
                        w
                    )
                        ? "PASS"
                        : "FAIL"
                )
                << "\n";
        }

        if (
            w.j <= 10 &&
            printed_small_examples < 5
        ) {
            ++printed_small_examples;

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
                << witness_to_string(
                    w
                )
                << "\n";

            std::cout
                << "IDENTITY="
                << (
                    verify_identity(
                        p,
                        w
                    )
                        ? "PASS"
                        : "FAIL"
                )
                << "\n";
        }

        if (
            w.j > 28 &&
            !printed_first_large_j
        ) {
            printed_first_large_j = true;

            std::cout
                << "\nFIRST_J_GT_28\n";

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
                    w
                )
                << "\n";
        }
    }

    if (
        min_j == UINT64_MAX
    ) {
        min_j = 0;
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
        << invalid_cases
        << "\n";

    std::cout
        << "IDENTITY_FAILURES="
        << identity_failures
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
        << "J_LE_3="
        << j_le_3
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
        << "J_LE_28="
        << j_le_28
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
        << "MAX_J_P="
        << max_j_p
        << "\n";

    std::cout
        << "MAX_J_Q="
        << max_j_q
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
