#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Pair {
    int m1 = 0;
    int m2 = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 M = 0;
    u64 base = 0;
};

struct Witness {
    bool valid = false;

    u64 d = 0;
    u64 j = 0;

    u64 M = 0;
    u64 base = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    int m1 = 0;
    int m2 = 0;
};

struct CaseResult {
    u64 p = 0;
    u64 q = 0;
    u64 N = 0;
    u64 s = 0;

    Witness min_d;
    Witness min_d_over_M;
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

static std::vector<Pair> generate_pairs(
    int M_LIMIT,
    u64 K_LIMIT
) {
    std::vector<Pair> pairs;

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
                    std::gcd(
                        k1,
                        k2
                    ) != 1
                ) {
                    continue;
                }

                Pair pair;

                pair.m1 = m1;
                pair.m2 = m2;
                pair.k1 = k1;
                pair.k2 = k2;

                pair.M =
                    k1 * k2;

                pair.base =
                    k1 + k2;

                pairs.push_back(pair);
            }
        }
    }

    /*
        Remove identical arithmetic progressions.
    */
    std::sort(
        pairs.begin(),
        pairs.end(),
        [](
            const Pair& a,
            const Pair& b
        ) {
            if (a.base != b.base) {
                return a.base <
                       b.base;
            }

            if (a.M != b.M) {
                return a.M <
                       b.M;
            }

            if (a.m1 != b.m1) {
                return a.m1 <
                       b.m1;
            }

            if (a.m2 != b.m2) {
                return a.m2 <
                       b.m2;
            }

            if (a.k1 != b.k1) {
                return a.k1 <
                       b.k1;
            }

            return a.k2 <
                   b.k2;
        }
    );

    std::vector<Pair> unique_pairs;

    for (const Pair& pair : pairs) {
        if (
            unique_pairs.empty() ||
            unique_pairs.back().base !=
                pair.base ||
            unique_pairs.back().M !=
                pair.M
        ) {
            unique_pairs.push_back(pair);
        }
    }

    return unique_pairs;
}

static Witness make_witness(
    u64 p,
    u64 s,
    const Pair& pair
) {
    Witness witness;

    /*
        p = base + j*M
    */
    if (p < pair.base) {
        return witness;
    }

    const u64 difference =
        p - pair.base;

    if (
        pair.M == 0 ||
        difference % pair.M != 0
    ) {
        return witness;
    }

    const u64 j =
        difference / pair.M;

    const u64 d =
        s - p;

    witness.valid = true;
    witness.d = d;
    witness.j = j;
    witness.M = pair.M;
    witness.base = pair.base;
    witness.k1 = pair.k1;
    witness.k2 = pair.k2;
    witness.m1 = pair.m1;
    witness.m2 = pair.m2;

    return witness;
}

static bool better_d(
    const Witness& a,
    const Witness& b
) {
    if (!a.valid) {
        return false;
    }

    if (!b.valid) {
        return true;
    }

    if (a.d != b.d) {
        return a.d < b.d;
    }

    if (a.M != b.M) {
        return a.M > b.M;
    }

    if (a.j != b.j) {
        return a.j < b.j;
    }

    return a.k1 > b.k1;
}

static bool better_d_over_M(
    const Witness& a,
    const Witness& b
) {
    if (!a.valid) {
        return false;
    }

    if (!b.valid) {
        return true;
    }

    const i128 lhs =
        static_cast<i128>(a.d) *
        static_cast<i128>(b.M);

    const i128 rhs =
        static_cast<i128>(b.d) *
        static_cast<i128>(a.M);

    if (lhs != rhs) {
        return lhs < rhs;
    }

    if (a.d != b.d) {
        return a.d < b.d;
    }

    return a.M > b.M;
}

static bool verify_witness(
    u64 p,
    const Witness& witness
) {
    if (!witness.valid) {
        return false;
    }

    const i128 reconstructed =
        static_cast<i128>(
            witness.base
        ) +
        static_cast<i128>(
            witness.j
        ) *
        static_cast<i128>(
            witness.M
        );

    return
        reconstructed ==
        static_cast<i128>(p);
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

    std::mt19937_64 rng(
        seed
    );

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

static CaseResult analyze_case(
    u64 p,
    u64 q,
    const std::vector<Pair>& pairs
) {
    CaseResult result;

    result.p = p;
    result.q = q;
    result.N = p * q;
    result.s = integer_sqrt(
        result.N
    );

    for (const Pair& pair : pairs) {
        const Witness witness =
            make_witness(
                p,
                result.s,
                pair
            );

        if (!witness.valid) {
            continue;
        }

        if (
            better_d(
                witness,
                result.min_d
            )
        ) {
            result.min_d =
                witness;
        }

        if (
            better_d_over_M(
                witness,
                result.min_d_over_M
            )
        ) {
            result.min_d_over_M =
                witness;
        }
    }

    return result;
}

static std::string witness_to_string(
    const Witness& witness
) {
    if (!witness.valid) {
        return "NONE";
    }

    std::string result;

    result +=
        "m1=" +
        std::to_string(
            witness.m1
        );

    result +=
        " m2=" +
        std::to_string(
            witness.m2
        );

    result +=
        " k1=" +
        std::to_string(
            witness.k1
        );

    result +=
        " k2=" +
        std::to_string(
            witness.k2
        );

    result +=
        " j=" +
        std::to_string(
            witness.j
        );

    result +=
        " M=" +
        std::to_string(
            witness.M
        );

    result +=
        " d=" +
        std::to_string(
            witness.d
        );

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
    constexpr int EXPERIMENT = 462;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 2000;

    constexpr int M_LIMIT = 7;
    constexpr u64 K_LIMIT = 1000;

    constexpr std::uint64_t SEED =
        0x462462462ULL;

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

    const auto primes =
        generate_primes(
            PRIME_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    const auto pairs =
        generate_pairs(
            M_LIMIT,
            K_LIMIT
        );

    std::cout
        << "UNIQUE_PAIRS="
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

    std::size_t d_lt_M = 0;
    std::size_t d_le_M = 0;
    std::size_t d_le_2M = 0;
    std::size_t d_le_5M = 0;
    std::size_t d_le_10M = 0;
    std::size_t d_gt_10M = 0;

    std::size_t d_zero = 0;

    u64 min_d = UINT64_MAX;
    u64 max_d = 0;
    u64 total_d = 0;

    u64 min_M = UINT64_MAX;
    u64 max_M = 0;
    u64 total_M = 0;

    u64 min_j = UINT64_MAX;
    u64 max_j = 0;
    u64 total_j = 0;

    u64 min_d_over_M_d = UINT64_MAX;
    u64 max_d_over_M_d = 0;

    u64 min_distance_case_p = UINT64_MAX;
    u64 max_distance_case_p = 0;

    std::size_t printed_small_ratio =
        0;

    bool printed_first = false;
    bool printed_first_large_ratio =
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
            !result.min_d.valid
        ) {
            continue;
        }

        ++valid_cases;

        const Witness& w =
            result.min_d;

        if (
            !verify_witness(
                p,
                w
            )
        ) {
            ++identity_failures;
        }

        min_d =
            std::min(
                min_d,
                w.d
            );

        max_d =
            std::max(
                max_d,
                w.d
            );

        total_d +=
            w.d;

        min_M =
            std::min(
                min_M,
                w.M
            );

        max_M =
            std::max(
                max_M,
                w.M
            );

        total_M +=
            w.M;

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

        total_j +=
            w.j;

        if (w.d == 0) {
            ++d_zero;
        }

        if (w.d < w.M) {
            ++d_lt_M;
        }

        if (w.d <= w.M) {
            ++d_le_M;
        }

        if (
            static_cast<i128>(w.d) <=
            static_cast<i128>(2) *
            static_cast<i128>(w.M)
        ) {
            ++d_le_2M;
        }

        if (
            static_cast<i128>(w.d) <=
            static_cast<i128>(5) *
            static_cast<i128>(w.M)
        ) {
            ++d_le_5M;
        }

        if (
            static_cast<i128>(w.d) <=
            static_cast<i128>(10) *
            static_cast<i128>(w.M)
        ) {
            ++d_le_10M;
        } else {
            ++d_gt_10M;
        }

        const u64 q_distance =
            w.M == 0
                ? 0
                : w.d / w.M;

        min_d_over_M_d =
            std::min(
                min_d_over_M_d,
                q_distance
            );

        max_d_over_M_d =
            std::max(
                max_d_over_M_d,
                q_distance
            );

        if (
            !printed_first
        ) {
            printed_first = true;

            std::cout
                << "\nFIRST_CASE\n";

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
                << "MIN_D="
                << witness_to_string(
                    result.min_d
                )
                << "\n";

            std::cout
                << "MIN_D_OVER_M="
                << witness_to_string(
                    result.min_d_over_M
                )
                << "\n";
        }

        /*
            Print a few especially strong cases
            where d < M.
        */
        if (
            w.d < w.M &&
            printed_small_ratio < 5
        ) {
            ++printed_small_ratio;

            std::cout
                << "\nD_LT_M_CASE\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " S="
                << result.s
                << "\n";

            std::cout
                << witness_to_string(
                    w
                )
                << "\n";

            std::cout
                << "D_OVER_M="
                << std::setprecision(10)
                << (
                    static_cast<double>(w.d) /
                    static_cast<double>(w.M)
                )
                << "\n";
        }

        if (
            w.d >
            10 * w.M &&
            !printed_first_large_ratio
        ) {
            printed_first_large_ratio =
                true;

            std::cout
                << "\nFIRST_D_GT_10M\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " S="
                << result.s
                << "\n";

            std::cout
                << witness_to_string(
                    w
                )
                << "\n";

            std::cout
                << "D_OVER_M="
                << std::setprecision(10)
                << (
                    static_cast<double>(w.d) /
                    static_cast<double>(w.M)
                )
                << "\n";
        }
    }

    if (
        min_d == UINT64_MAX
    ) {
        min_d = 0;
    }

    if (
        min_M == UINT64_MAX
    ) {
        min_M = 0;
    }

    if (
        min_j == UINT64_MAX
    ) {
        min_j = 0;
    }

    if (
        min_d_over_M_d == UINT64_MAX
    ) {
        min_d_over_M_d = 0;
    }

    const double average_d =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                  total_d
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double average_M =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                  total_M
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double average_j =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                  total_j
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double average_d_over_M =
        valid_cases == 0
            ? 0.0
            : (
                average_d /
                average_M
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
        << "IDENTITY_FAILURES="
        << identity_failures
        << "\n";

    std::cout
        << "D_ZERO="
        << d_zero
        << "\n";

    std::cout
        << "D_LT_M="
        << d_lt_M
        << "\n";

    std::cout
        << "D_LE_M="
        << d_le_M
        << "\n";

    std::cout
        << "D_LE_2M="
        << d_le_2M
        << "\n";

    std::cout
        << "D_LE_5M="
        << d_le_5M
        << "\n";

    std::cout
        << "D_LE_10M="
        << d_le_10M
        << "\n";

    std::cout
        << "D_GT_10M="
        << d_gt_10M
        << "\n";

    std::cout
        << "MIN_D="
        << min_d
        << "\n";

    std::cout
        << "MAX_D="
        << max_d
        << "\n";

    std::cout
        << "AVERAGE_D="
        << average_d
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
        << "MIN_FLOOR_D_OVER_M="
        << min_d_over_M_d
        << "\n";

    std::cout
        << "MAX_FLOOR_D_OVER_M="
        << max_d_over_M_d
        << "\n";

    std::cout
        << "AVERAGE_D_OVER_M="
        << std::setprecision(10)
        << average_d_over_M
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
