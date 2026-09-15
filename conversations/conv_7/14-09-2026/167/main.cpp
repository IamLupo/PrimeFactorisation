#include <algorithm>
#include <cstdint>
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

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 candidate = 0;

    u64 j = 0;
    u64 M = 0;

    u64 candidate_tests = 0;
    u64 pair_rank = 0;

    std::size_t pair_index = 0;
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
        Remove duplicate arithmetic progressions.

            p = base + j*M

        has the same sequence whenever (base,M)
        is identical.
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

            return
                std::tie(
                    a.m1,
                    a.m2,
                    a.k1,
                    a.k2
                ) <
                std::tie(
                    b.m1,
                    b.m2,
                    b.k1,
                    b.k2
                );
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
            unique_pairs.push_back(
                pair
            );
        }
    }

    return unique_pairs;
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

/*
    Strategy 1:
        M descending.

    Large M means fewer p candidates per progression.
*/
static Hit search_M_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::vector<std::size_t> order(
        pairs.size()
    );

    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        order[i] = i;
    }

    std::sort(
        order.begin(),
        order.end(),
        [&pairs](
            std::size_t a,
            std::size_t b
        ) {
            if (
                pairs[a].M !=
                pairs[b].M
            ) {
                return pairs[a].M >
                       pairs[b].M;
            }

            if (
                pairs[a].base !=
                pairs[b].base
            ) {
                return pairs[a].base <
                       pairs[b].base;
            }

            return a < b;
        }
    );

    for (
        std::size_t rank = 0;
        rank < order.size();
        ++rank
    ) {
        const std::size_t index =
            order[rank];

        const Pair& pair =
            pairs[index];

        const u64 max_j =
            pair.base > s
                ? 0
                : (
                    s - pair.base
                ) /
                  pair.M;

        if (max_j == 0) {
            continue;
        }

        for (
            u64 j = 1;
            j <= max_j;
            ++j
        ) {
            ++hit.candidate_tests;

            const u64 candidate =
                pair.base +
                j * pair.M;

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
                hit.candidate = candidate;

                hit.j = j;
                hit.M = pair.M;

                hit.pair_rank =
                    static_cast<u64>(
                        rank + 1
                    );

                hit.pair_index =
                    index;

                return hit;
            }
        }
    }

    return hit;
}

/*
    Strategy 2:
        Work-first.

    Approximate amount of work for a pair is

        floor((sqrt(N)-base)/M).

    Search pairs requiring the least work first.
*/
static Hit search_work_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    struct RankedPair {
        u64 work = 0;
        std::size_t index = 0;
    };

    std::vector<RankedPair> order;

    order.reserve(
        pairs.size()
    );

    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        const Pair& pair =
            pairs[i];

        if (
            pair.base >= s
        ) {
            continue;
        }

        const u64 work =
            (
                s -
                pair.base
            ) /
            pair.M;

        if (work == 0) {
            continue;
        }

        RankedPair ranked;

        ranked.work = work;
        ranked.index = i;

        order.push_back(
            ranked
        );
    }

    std::sort(
        order.begin(),
        order.end(),
        [&pairs](
            const RankedPair& a,
            const RankedPair& b
        ) {
            if (a.work != b.work) {
                return a.work <
                       b.work;
            }

            if (
                pairs[a.index].M !=
                pairs[b.index].M
            ) {
                return pairs[a.index].M >
                       pairs[b.index].M;
            }

            return a.index <
                   b.index;
        }
    );

    for (
        std::size_t rank = 0;
        rank < order.size();
        ++rank
    ) {
        const std::size_t index =
            order[rank].index;

        const Pair& pair =
            pairs[index];

        const u64 max_j =
            (
                s -
                pair.base
            ) /
            pair.M;

        for (
            u64 j = 1;
            j <= max_j;
            ++j
        ) {
            ++hit.candidate_tests;

            const u64 candidate =
                pair.base +
                j * pair.M;

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
                hit.candidate = candidate;

                hit.j = j;
                hit.M = pair.M;

                hit.pair_rank =
                    static_cast<u64>(
                        rank + 1
                    );

                hit.pair_index =
                    index;

                return hit;
            }
        }
    }

    return hit;
}

/*
    Strategy 3:
        smallest base/M ratio.

    Since

        p = base + j*M,

    this favors progressions whose first
    terms are relatively large compared with
    their step.
*/
static Hit search_density_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::vector<std::size_t> order(
        pairs.size()
    );

    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        order[i] = i;
    }

    std::sort(
        order.begin(),
        order.end(),
        [&pairs](
            std::size_t a,
            std::size_t b
        ) {
            const i128 lhs =
                static_cast<i128>(
                    pairs[a].base
                ) *
                static_cast<i128>(
                    pairs[b].M
                );

            const i128 rhs =
                static_cast<i128>(
                    pairs[b].base
                ) *
                static_cast<i128>(
                    pairs[a].M
                );

            if (lhs != rhs) {
                return lhs > rhs;
            }

            if (
                pairs[a].M !=
                pairs[b].M
            ) {
                return pairs[a].M >
                       pairs[b].M;
            }

            return a < b;
        }
    );

    for (
        std::size_t rank = 0;
        rank < order.size();
        ++rank
    ) {
        const std::size_t index =
            order[rank];

        const Pair& pair =
            pairs[index];

        if (
            pair.base >= s
        ) {
            continue;
        }

        const u64 max_j =
            (
                s -
                pair.base
            ) /
            pair.M;

        for (
            u64 j = 1;
            j <= max_j;
            ++j
        ) {
            ++hit.candidate_tests;

            const u64 candidate =
                pair.base +
                j * pair.M;

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
                hit.candidate = candidate;

                hit.j = j;
                hit.M = pair.M;

                hit.pair_rank =
                    static_cast<u64>(
                        rank + 1
                    );

                hit.pair_index =
                    index;

                return hit;
            }
        }
    }

    return hit;
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
            usable.push_back(
                p
            );
        }
    }

    std::mt19937_64 rng(
        seed
    );

    std::vector<
        std::pair<u64, u64>
    > cases;

    cases.reserve(
        count
    );

    for (
        std::size_t i = 0;
        i < count;
        ++i
    ) {
        std::size_t a =
            static_cast<std::size_t>(
                rng() %
                usable.size()
            );

        std::size_t b =
            static_cast<std::size_t>(
                rng() %
                usable.size()
            );

        while (
            a == b
        ) {
            b =
                static_cast<std::size_t>(
                    rng() %
                    usable.size()
                );
        }

        if (
            usable[a] >
            usable[b]
        ) {
            std::swap(
                a,
                b
            );
        }

        cases.emplace_back(
            usable[a],
            usable[b]
        );
    }

    return cases;
}

static void print_hit(
    const char* name,
    const Hit& hit,
    u64 p,
    u64 q
) {
    std::cout
        << "\n"
        << name
        << "\n";

    std::cout
        << "P="
        << p
        << " Q="
        << q
        << "\n";

    if (!hit.valid) {
        std::cout
            << "HIT=NONE\n";

        return;
    }

    std::cout
        << "FACTOR="
        << hit.factor
        << "\n";

    std::cout
        << "CANDIDATE="
        << hit.candidate
        << "\n";

    std::cout
        << "J="
        << hit.j
        << "\n";

    std::cout
        << "M="
        << hit.M
        << "\n";

    std::cout
        << "PAIR_RANK="
        << hit.pair_rank
        << "\n";

    std::cout
        << "CANDIDATE_TESTS="
        << hit.candidate_tests
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 460;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 500;

    constexpr int M_LIMIT = 7;
    constexpr u64 K_LIMIT = 1000;

    constexpr std::uint64_t SEED =
        0x460460460ULL;

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

    u64 hits_M = 0;
    u64 hits_work = 0;
    u64 hits_density = 0;

    u64 total_M = 0;
    u64 total_work = 0;
    u64 total_density = 0;

    u64 min_M = UINT64_MAX;
    u64 min_work = UINT64_MAX;
    u64 min_density = UINT64_MAX;

    u64 max_M = 0;
    u64 max_work = 0;
    u64 max_density = 0;

    u64 total_j_M = 0;
    u64 total_j_work = 0;
    u64 total_j_density = 0;

    u64 min_j_M = UINT64_MAX;
    u64 min_j_work = UINT64_MAX;
    u64 min_j_density = UINT64_MAX;

    u64 max_j_M = 0;
    u64 max_j_work = 0;
    u64 max_j_density = 0;

    bool printed_first = false;

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

        const Hit hit_M =
            search_M_first(
                N,
                s,
                pairs
            );

        const Hit hit_work =
            search_work_first(
                N,
                s,
                pairs
            );

        const Hit hit_density =
            search_density_first(
                N,
                s,
                pairs
            );

        if (
            hit_M.valid
        ) {
            ++hits_M;

            total_M +=
                hit_M.candidate_tests;

            min_M =
                std::min(
                    min_M,
                    hit_M.candidate_tests
                );

            max_M =
                std::max(
                    max_M,
                    hit_M.candidate_tests
                );

            total_j_M +=
                hit_M.j;

            min_j_M =
                std::min(
                    min_j_M,
                    hit_M.j
                );

            max_j_M =
                std::max(
                    max_j_M,
                    hit_M.j
                );
        }

        if (
            hit_work.valid
        ) {
            ++hits_work;

            total_work +=
                hit_work.candidate_tests;

            min_work =
                std::min(
                    min_work,
                    hit_work.candidate_tests
                );

            max_work =
                std::max(
                    max_work,
                    hit_work.candidate_tests
                );

            total_j_work +=
                hit_work.j;

            min_j_work =
                std::min(
                    min_j_work,
                    hit_work.j
                );

            max_j_work =
                std::max(
                    max_j_work,
                    hit_work.j
                );
        }

        if (
            hit_density.valid
        ) {
            ++hits_density;

            total_density +=
                hit_density.candidate_tests;

            min_density =
                std::min(
                    min_density,
                    hit_density.candidate_tests
                );

            max_density =
                std::max(
                    max_density,
                    hit_density.candidate_tests
                );

            total_j_density +=
                hit_density.j;

            min_j_density =
                std::min(
                    min_j_density,
                    hit_density.j
                );

            max_j_density =
                std::max(
                    max_j_density,
                    hit_density.j
                );
        }

        if (
            !printed_first
        ) {
            printed_first = true;

            print_hit(
                "FIRST_M_HIT",
                hit_M,
                p,
                q
            );

            print_hit(
                "FIRST_WORK_HIT",
                hit_work,
                p,
                q
            );

            print_hit(
                "FIRST_DENSITY_HIT",
                hit_density,
                p,
                q
            );
        }

        if (
            i == 0 ||
            i % 100 == 0 ||
            i + 1 == cases.size()
        ) {
            std::cout
                << "PROGRESS "
                << i + 1
                << "/"
                << cases.size()
                << "\n";
        }
    }

    if (min_M == UINT64_MAX) {
        min_M = 0;
    }

    if (min_work == UINT64_MAX) {
        min_work = 0;
    }

    if (min_density == UINT64_MAX) {
        min_density = 0;
    }

    if (min_j_M == UINT64_MAX) {
        min_j_M = 0;
    }

    if (min_j_work == UINT64_MAX) {
        min_j_work = 0;
    }

    if (min_j_density == UINT64_MAX) {
        min_j_density = 0;
    }

    const double avg_M =
        hits_M == 0
            ? 0.0
            : static_cast<double>(
                  total_M
              ) /
              static_cast<double>(
                  hits_M
              );

    const double avg_work =
        hits_work == 0
            ? 0.0
            : static_cast<double>(
                  total_work
              ) /
              static_cast<double>(
                  hits_work
              );

    const double avg_density =
        hits_density == 0
            ? 0.0
            : static_cast<double>(
                  total_density
              ) /
              static_cast<double>(
                  hits_density
              );

    const double avg_j_M =
        hits_M == 0
            ? 0.0
            : static_cast<double>(
                  total_j_M
              ) /
              static_cast<double>(
                  hits_M
              );

    const double avg_j_work =
        hits_work == 0
            ? 0.0
            : static_cast<double>(
                  total_j_work
              ) /
              static_cast<double>(
                  hits_work
              );

    const double avg_j_density =
        hits_density == 0
            ? 0.0
            : static_cast<double>(
                  total_j_density
              ) /
              static_cast<double>(
                  hits_density
              );

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "\nM_FIRST\n";

    std::cout
        << "HITS="
        << hits_M
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << avg_M
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_M
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_M
        << "\n";

    std::cout
        << "AVERAGE_J="
        << avg_j_M
        << "\n";

    std::cout
        << "MIN_J="
        << min_j_M
        << "\n";

    std::cout
        << "MAX_J="
        << max_j_M
        << "\n";

    std::cout
        << "\nWORK_FIRST\n";

    std::cout
        << "HITS="
        << hits_work
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << avg_work
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_work
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_work
        << "\n";

    std::cout
        << "AVERAGE_J="
        << avg_j_work
        << "\n";

    std::cout
        << "MIN_J="
        << min_j_work
        << "\n";

    std::cout
        << "MAX_J="
        << max_j_work
        << "\n";

    std::cout
        << "\nDENSITY_FIRST\n";

    std::cout
        << "HITS="
        << hits_density
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << avg_density
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_density
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_density
        << "\n";

    std::cout
        << "AVERAGE_J="
        << avg_j_density
        << "\n";

    std::cout
        << "MIN_J="
        << min_j_density
        << "\n";

    std::cout
        << "MAX_J="
        << max_j_density
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}