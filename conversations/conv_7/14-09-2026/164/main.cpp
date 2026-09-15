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

struct RootClass {
    u64 modulus = 0;
    u64 residue = 0;
};

struct ParameterSet {
    u64 j = 0;
    int m1 = 0;
    int m2 = 0;

    u64 a = 0;
    u64 b = 0;
    u64 A = 0;
    u64 B = 0;

    std::vector<RootClass> roots;
};

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 candidate = 0;

    u64 j = 0;

    int m1 = 0;
    int m2 = 0;

    u64 X = 0;
    u64 A = 0;
    u64 B = 0;

    u64 root_residue = 0;
    u64 root_index = 0;

    u64 candidate_rank = 0;
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

static std::vector<RootClass>
find_square_roots_mod(
    u64 B,
    u64 modulus
) {
    std::vector<RootClass> roots;

    if (modulus == 0) {
        return roots;
    }

    /*
        Current parameter bounds give very small
        moduli, so direct enumeration is sufficient.
    */

    for (
        u64 x = 0;
        x < modulus;
        ++x
    ) {
        const u64 xx =
            static_cast<u64>(
                (
                    static_cast<i128>(x) *
                    static_cast<i128>(x)
                ) %
                static_cast<i128>(modulus)
            );

        if (
            xx ==
            B % modulus
        ) {
            RootClass root;

            root.modulus = modulus;
            root.residue = x;

            roots.push_back(root);
        }
    }

    return roots;
}

static std::vector<ParameterSet>
generate_parameters(
    u64 J_LIMIT,
    int M_LIMIT
) {
    std::vector<ParameterSet> parameters;

    for (
        u64 j = 1;
        j <= J_LIMIT;
        ++j
    ) {
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
                    Starting from

                    p =
                        k1 + k2 + j*k1*k2

                    and

                    m2*k1 - m1*k2 = 1,

                    the resulting quadratic is

                    j*m2*k1^2
                    + (m1+m2-j)*k1
                    - (m1*p+1) = 0.

                    Its discriminant is

                    X^2 =
                        4*j*m1*m2*p
                        + (m1+m2-j)^2
                        + 4*j*m2.
                */

                const i128 A_i =
                    static_cast<i128>(4) *
                    static_cast<i128>(j) *
                    static_cast<i128>(m1) *
                    static_cast<i128>(m2);

                const i128 linear_i =
                    static_cast<i128>(m1) +
                    static_cast<i128>(m2) -
                    static_cast<i128>(j);

                const i128 B_i =
                    linear_i * linear_i +
                    static_cast<i128>(4) *
                    static_cast<i128>(j) *
                    static_cast<i128>(m2);

                if (
                    A_i <= 0 ||
                    A_i >
                        static_cast<i128>(
                            UINT64_MAX
                        )
                ) {
                    continue;
                }

                if (
                    B_i < 0 ||
                    B_i >
                        static_cast<i128>(
                            UINT64_MAX
                        )
                ) {
                    continue;
                }

                ParameterSet parameter;

                parameter.j = j;
                parameter.m1 = m1;
                parameter.m2 = m2;

                parameter.a =
                    j *
                    static_cast<u64>(m2);

                parameter.b =
                    linear_i >= 0
                        ? static_cast<u64>(
                              linear_i
                          )
                        : static_cast<u64>(
                              -linear_i
                          );

                parameter.A =
                    static_cast<u64>(
                        A_i
                    );

                parameter.B =
                    static_cast<u64>(
                        B_i
                    );

                parameter.roots =
                    find_square_roots_mod(
                        parameter.B,
                        parameter.A
                    );

                if (
                    parameter.roots.empty()
                ) {
                    continue;
                }

                parameters.push_back(
                    parameter
                );
            }
        }
    }

    std::sort(
        parameters.begin(),
        parameters.end(),
        [](
            const ParameterSet& a,
            const ParameterSet& b
        ) {
            if (a.j != b.j) {
                return a.j < b.j;
            }

            if (a.m1 != b.m1) {
                return a.m1 < b.m1;
            }

            return a.m2 < b.m2;
        }
    );

    return parameters;
}

static Hit blind_search(
    u64 N,
    u64 s,
    const std::vector<ParameterSet>& parameters
) {
    Hit hit;

    u64 candidate_rank = 0;

    for (
        const ParameterSet& parameter :
        parameters
    ) {
        /*
            X^2 = A*p + B

            and p <= sqrt(N).
        */

        const i128 max_square_i =
            static_cast<i128>(
                parameter.A
            ) *
            static_cast<i128>(s) +
            static_cast<i128>(
                parameter.B
            );

        if (
            max_square_i <= 0 ||
            max_square_i >
                static_cast<i128>(
                    UINT64_MAX
                )
        ) {
            continue;
        }

        const u64 max_square =
            static_cast<u64>(
                max_square_i
            );

        const u64 X_max =
            integer_sqrt(
                max_square
            );

        for (
            const RootClass& root :
            parameter.roots
        ) {
            u64 X =
                root.residue;

            if (X == 0) {
                X = root.modulus;
            }

            const u64 two_a =
                static_cast<u64>(2) *
                parameter.a;

            if (two_a == 0) {
                continue;
            }

            while (X <= X_max) {
                /*
                    X =
                        2*j*m2*k1
                        + (m1+m2-j)

                    hence k1 must be integral.
                */

                const i128 linear_i =
                    static_cast<i128>(
                        parameter.m1
                    ) +
                    static_cast<i128>(
                        parameter.m2
                    ) -
                    static_cast<i128>(
                        parameter.j
                    );

                const i128 numerator =
                    static_cast<i128>(X) -
                    linear_i;

                if (
                    numerator >= 0 &&
                    numerator %
                        static_cast<i128>(
                            two_a
                        ) ==
                        0
                ) {
                    const i128 k1_i =
                        numerator /
                        static_cast<i128>(
                            two_a
                        );

                    if (
                        k1_i >= 2 &&
                        k1_i <=
                            static_cast<i128>(
                                UINT64_MAX
                            )
                    ) {
                        const u64 k1 =
                            static_cast<u64>(
                                k1_i
                            );

                        /*
                            m2*k1 - m1*k2 = 1
                        */

                        const i128 k2_num =
                            static_cast<i128>(
                                parameter.m2
                            ) *
                            static_cast<i128>(
                                k1
                            ) -
                            static_cast<i128>(1);

                        if (
                            k2_num > 0 &&
                            k2_num %
                                static_cast<i128>(
                                    parameter.m1
                                ) ==
                                0
                        ) {
                            const i128 k2_i =
                                k2_num /
                                static_cast<i128>(
                                    parameter.m1
                                );

                            if (
                                k2_i >= 2 &&
                                k2_i <=
                                    static_cast<i128>(
                                        UINT64_MAX
                                    )
                            ) {
                                const u64 k2 =
                                    static_cast<u64>(
                                        k2_i
                                    );

                                if (
                                    determinant_one(
                                        parameter.m1,
                                        k1,
                                        parameter.m2,
                                        k2
                                    )
                                ) {
                                    /*
                                        Candidate:

                                        p =
                                            k1 + k2
                                            + j*k1*k2
                                    */

                                    const i128 candidate_i =
                                        static_cast<i128>(
                                            k1
                                        ) +
                                        static_cast<i128>(
                                            k2
                                        ) +
                                        static_cast<i128>(
                                            parameter.j
                                        ) *
                                        static_cast<i128>(
                                            k1
                                        ) *
                                        static_cast<i128>(
                                            k2
                                        );

                                    if (
                                        candidate_i >= 2 &&
                                        candidate_i <=
                                            static_cast<i128>(
                                                s
                                            )
                                    ) {
                                        const u64 candidate =
                                            static_cast<u64>(
                                                candidate_i
                                            );

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
                                            const i128 lhs =
                                                static_cast<i128>(
                                                    X
                                                ) *
                                                static_cast<i128>(
                                                    X
                                                );

                                            const i128 rhs =
                                                static_cast<i128>(
                                                    parameter.A
                                                ) *
                                                static_cast<i128>(
                                                    candidate
                                                ) +
                                                static_cast<i128>(
                                                    parameter.B
                                                );

                                            if (
                                                lhs == rhs
                                            ) {
                                                hit.valid = true;

                                                hit.factor = g;
                                                hit.candidate =
                                                    candidate;

                                                hit.j =
                                                    parameter.j;

                                                hit.m1 =
                                                    parameter.m1;

                                                hit.m2 =
                                                    parameter.m2;

                                                hit.X = X;
                                                hit.A =
                                                    parameter.A;
                                                hit.B =
                                                    parameter.B;

                                                hit.root_residue =
                                                    root.residue;

                                                hit.root_index =
                                                    (
                                                        X -
                                                        root.residue
                                                    ) /
                                                    root.modulus;

                                                hit.candidate_rank =
                                                    candidate_rank;

                                                return hit;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (
                    X >
                    X_max -
                    root.modulus
                ) {
                    break;
                }

                X += root.modulus;
            }
        }
    }

    return hit;
}

static std::string hit_to_string(
    const Hit& hit
) {
    if (!hit.valid) {
        return "NONE";
    }

    std::string out;

    out +=
        "j=" +
        std::to_string(hit.j);

    out +=
        " m1=" +
        std::to_string(hit.m1);

    out +=
        " m2=" +
        std::to_string(hit.m2);

    out +=
        " factor=" +
        std::to_string(hit.factor);

    out +=
        " candidate=" +
        std::to_string(hit.candidate);

    out +=
        " X=" +
        std::to_string(hit.X);

    out +=
        " A=" +
        std::to_string(hit.A);

    out +=
        " B=" +
        std::to_string(hit.B);

    out +=
        " candidate_rank=" +
        std::to_string(
            hit.candidate_rank
        );

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
            static_cast<u64>(
                prime_min
            )
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
    constexpr int EXPERIMENT = 456;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 2000;

    constexpr u64 J_LIMIT = 28;
    constexpr int M_LIMIT = 7;

    constexpr std::uint64_t SEED =
        0x456456456ULL;

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
        << "J_LIMIT="
        << J_LIMIT
        << "\n";

    std::cout
        << "M_LIMIT="
        << M_LIMIT
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

    const std::vector<ParameterSet>
        parameters =
            generate_parameters(
                J_LIMIT,
                M_LIMIT
            );

    std::size_t total_root_classes = 0;

    for (
        const ParameterSet& parameter :
        parameters
    ) {
        total_root_classes +=
            parameter.roots.size();
    }

    std::cout
        << "PARAMETER_SETS="
        << parameters.size()
        << "\n";

    std::cout
        << "TOTAL_ROOT_CLASSES="
        << total_root_classes
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
    std::size_t true_factor_hits = 0;

    std::size_t j_le_1 = 0;
    std::size_t j_le_2 = 0;
    std::size_t j_le_5 = 0;
    std::size_t j_le_10 = 0;
    std::size_t j_le_20 = 0;
    std::size_t j_le_28 = 0;

    u64 min_j = UINT64_MAX;
    u64 max_j = 0;
    u64 total_j = 0;

    u64 min_candidate_rank =
        UINT64_MAX;

    u64 max_candidate_rank = 0;
    u64 total_candidate_rank = 0;

    std::set<
        std::tuple<
            u64,
            int,
            int
        >
    > hit_parameter_types;

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

        const u64 N =
            p * q;

        const u64 s =
            integer_sqrt(N);

        const Hit hit =
            blind_search(
                N,
                s,
                parameters
            );

        print_progress(
            i,
            cases.size()
        );

        if (!hit.valid) {
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

        if (
            hit.factor == p ||
            hit.factor == q
        ) {
            ++true_factor_hits;
        }

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

        if (hit.j <= 1) {
            ++j_le_1;
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

        if (hit.j <= 28) {
            ++j_le_28;
        }

        hit_parameter_types.insert(
            std::make_tuple(
                hit.j,
                hit.m1,
                hit.m2
            )
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
                << hit_to_string(hit)
                << "\n";

            std::cout
                << "TRUE_FACTOR="
                << (
                    hit.factor == p ||
                    hit.factor == q
                        ? "YES"
                        : "NO"
                )
                << "\n";
        }
    }

    if (
        min_j == UINT64_MAX
    ) {
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
        << true_factor_hits
        << "\n";

    std::cout
        << "PARAMETER_SETS="
        << parameters.size()
        << "\n";

    std::cout
        << "TOTAL_ROOT_CLASSES="
        << total_root_classes
        << "\n";

    std::cout
        << "HIT_PARAMETER_TYPES="
        << hit_parameter_types.size()
        << "\n";

    std::cout
        << "J_LE_1="
        << j_le_1
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
        << "J_LE_28="
        << j_le_28
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