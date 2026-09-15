#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <tuple>
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
    u64 pair_states_examined = 0;

    std::size_t pair_index = 0;
};

struct CaseResult {
    Hit j_first;
    Hit M_first;
    Hit jM_first;
};

struct HeapNode {
    u64 key = 0;
    u64 j = 0;
    std::size_t pair_index = 0;
};

struct HeapCompare {
    bool operator()(
        const HeapNode& a,
        const HeapNode& b
    ) const {
        if (a.key != b.key) {
            return a.key > b.key;
        }

        if (a.j != b.j) {
            return a.j > b.j;
        }

        return a.pair_index >
               b.pair_index;
    }
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

static std::vector<Pair>
generate_pairs(
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
                    std::gcd(k1, k2) != 1
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
        Deduplicate identical (base,M) sequences.
        Structural copies of the same arithmetic
        progression need not be tested repeatedly.
    */

    std::sort(
        pairs.begin(),
        pairs.end(),
        [](
            const Pair& a,
            const Pair& b
        ) {
            if (a.base != b.base) {
                return a.base < b.base;
            }

            if (a.M != b.M) {
                return a.M < b.M;
            }

            if (a.m1 != b.m1) {
                return a.m1 < b.m1;
            }

            if (a.m2 != b.m2) {
                return a.m2 < b.m2;
            }

            if (a.k1 != b.k1) {
                return a.k1 < b.k1;
            }

            return a.k2 < b.k2;
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

static Hit search_j_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs,
    u64 J_LIMIT
) {
    Hit hit;

    /*
        At fixed j:

            candidate = base + j*M.

        Search j=1,2,... and all structural
        progressions at that j.
    */

    for (
        u64 j = 1;
        j <= J_LIMIT;
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < pairs.size();
            ++i
        ) {
            const Pair& pair =
                pairs[i];

            const i128 candidate_i =
                static_cast<i128>(
                    pair.base
                ) +
                static_cast<i128>(j) *
                static_cast<i128>(pair.M);

            if (
                candidate_i >
                static_cast<i128>(s)
            ) {
                continue;
            }

            if (
                candidate_i < 2
            ) {
                continue;
            }

            ++hit.candidate_tests;

            const u64 candidate =
                static_cast<u64>(
                    candidate_i
                );

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
                hit.pair_index = i;

                return hit;
            }
        }
    }

    return hit;
}

static Hit search_M_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs,
    u64 J_LIMIT
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

    for (std::size_t rank = 0;
         rank < order.size();
         ++rank) {

        const std::size_t i =
            order[rank];

        const Pair& pair =
            pairs[i];

        ++hit.pair_states_examined;

        if (
            pair.base >
            s
        ) {
            continue;
        }

        u64 max_j =
            (
                s -
                pair.base
            ) /
            pair.M;

        if (max_j > J_LIMIT) {
            max_j = J_LIMIT;
        }

        for (
            u64 j = 1;
            j <= max_j;
            ++j
        ) {
            const u64 candidate =
                pair.base +
                j * pair.M;

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
                hit.candidate = candidate;
                hit.j = j;
                hit.M = pair.M;
                hit.pair_index = i;

                return hit;
            }
        }
    }

    return hit;
}

static Hit search_jM_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs,
    u64 J_LIMIT
) {
    Hit hit;

    std::priority_queue<
        HeapNode,
        std::vector<HeapNode>,
        HeapCompare
    > heap;

    /*
        Each pair generates

            cost = j*M

        for j=1,2,... .

        Merge all sequences with a min-heap.
    */

    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        const Pair& pair =
            pairs[i];

        if (
            pair.base +
            pair.M >
            s
        ) {
            continue;
        }

        HeapNode node;

        node.key = pair.M;
        node.j = 1;
        node.pair_index = i;

        heap.push(node);
    }

    while (!heap.empty()) {
        const HeapNode node =
            heap.top();

        heap.pop();

        const Pair& pair =
            pairs[node.pair_index];

        const u64 j =
            node.j;

        const i128 candidate_i =
            static_cast<i128>(
                pair.base
            ) +
            static_cast<i128>(j) *
            static_cast<i128>(pair.M);

        if (
            candidate_i >= 2 &&
            candidate_i <=
                static_cast<i128>(s)
        ) {
            ++hit.candidate_tests;

            const u64 candidate =
                static_cast<u64>(
                    candidate_i
                );

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
                hit.pair_index =
                    node.pair_index;

                return hit;
            }
        }

        if (
            j < J_LIMIT
        ) {
            const u64 next_j =
                j + 1;

            const i128 next_key_i =
                static_cast<i128>(
                    next_j
                ) *
                static_cast<i128>(
                    pair.M
                );

            if (
                pair.base +
                next_key_i <=
                static_cast<i128>(s)
            ) {
                if (
                    next_key_i <=
                    static_cast<i128>(
                        UINT64_MAX
                    )
                ) {
                    HeapNode next;

                    next.key =
                        static_cast<u64>(
                            next_key_i
                        );

                    next.j =
                        next_j;

                    next.pair_index =
                        node.pair_index;

                    heap.push(next);
                }
            }
        }
    }

    return hit;
}

static std::string
hit_to_string(
    const Hit& hit
) {
    if (!hit.valid) {
        return "NONE";
    }

    std::string out;

    out +=
        "factor=" +
        std::to_string(
            hit.factor
        );

    out +=
        " candidate=" +
        std::to_string(
            hit.candidate
        );

    out +=
        " j=" +
        std::to_string(
            hit.j
        );

    out +=
        " M=" +
        std::to_string(
            hit.M
        );

    out +=
        " candidate_tests=" +
        std::to_string(
            hit.candidate_tests
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

static void update_stats(
    const Hit& hit,
    u64& hits,
    u64& total_candidates,
    u64& min_candidates,
    u64& max_candidates,
    u64& total_j,
    u64& min_j,
    u64& max_j,
    u64& total_M,
    u64& min_M,
    u64& max_M
) {
    if (!hit.valid) {
        return;
    }

    ++hits;

    total_candidates +=
        hit.candidate_tests;

    min_candidates =
        std::min(
            min_candidates,
            hit.candidate_tests
        );

    max_candidates =
        std::max(
            max_candidates,
            hit.candidate_tests
        );

    total_j += hit.j;

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

    total_M += hit.M;

    min_M =
        std::min(
            min_M,
            hit.M
        );

    max_M =
        std::max(
            max_M,
            hit.M
        );
}

int main() {
    constexpr int EXPERIMENT = 458;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 500;

    constexpr int M_LIMIT = 7;
    constexpr u64 K_LIMIT = 1000;
    constexpr u64 J_LIMIT = 5000;

    constexpr std::uint64_t SEED =
        0x458458458ULL;

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
        << "J_LIMIT="
        << J_LIMIT
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

    const std::vector<Pair> pairs =
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

    u64 j_hits = 0;
    u64 m_hits = 0;
    u64 jm_hits = 0;

    u64 j_failures = 0;
    u64 m_failures = 0;
    u64 jm_failures = 0;

    u64 total_j_candidates = 0;
    u64 total_m_candidates = 0;
    u64 total_jm_candidates = 0;

    u64 min_j_candidates = UINT64_MAX;
    u64 min_m_candidates = UINT64_MAX;
    u64 min_jm_candidates = UINT64_MAX;

    u64 max_j_candidates = 0;
    u64 max_m_candidates = 0;
    u64 max_jm_candidates = 0;

    u64 j_total_j = 0;
    u64 m_total_j = 0;
    u64 jm_total_j = 0;

    u64 j_min_j = UINT64_MAX;
    u64 m_min_j = UINT64_MAX;
    u64 jm_min_j = UINT64_MAX;

    u64 j_max_j = 0;
    u64 m_max_j = 0;
    u64 jm_max_j = 0;

    u64 j_total_M = 0;
    u64 m_total_M = 0;
    u64 jm_total_M = 0;

    u64 j_min_M = UINT64_MAX;
    u64 m_min_M = UINT64_MAX;
    u64 jm_min_M = UINT64_MAX;

    u64 j_max_M = 0;
    u64 m_max_M = 0;
    u64 jm_max_M = 0;

    bool printed_first_j = false;
    bool printed_first_m = false;
    bool printed_first_jm = false;

    bool printed_miss = false;

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
            [&]() {
                u64 lo = 0;
                u64 hi = N;
                u64 ans = 0;

                while (lo <= hi) {
                    const u64 mid =
                        lo +
                        (hi - lo) / 2;

                    const i128 sq =
                        static_cast<i128>(
                            mid
                        ) *
                        static_cast<i128>(
                            mid
                        );

                    if (
                        sq <=
                        static_cast<i128>(
                            N
                        )
                    ) {
                        ans = mid;
                        lo = mid + 1;
                    } else {
                        if (mid == 0) {
                            break;
                        }

                        hi = mid - 1;
                    }
                }

                return ans;
            }();

        const Hit hit_j =
            search_j_first(
                N,
                s,
                pairs,
                J_LIMIT
            );

        const Hit hit_m =
            search_M_first(
                N,
                s,
                pairs,
                J_LIMIT
            );

        const Hit hit_jm =
            search_jM_first(
                N,
                s,
                pairs,
                J_LIMIT
            );

        print_progress(
            i,
            cases.size()
        );

        if (!hit_j.valid) {
            ++j_failures;
        }

        if (!hit_m.valid) {
            ++m_failures;
        }

        if (!hit_jm.valid) {
            ++jm_failures;
        }

        update_stats(
            hit_j,
            j_hits,
            total_j_candidates,
            min_j_candidates,
            max_j_candidates,
            j_total_j,
            j_min_j,
            j_max_j,
            j_total_M,
            j_min_M,
            j_max_M
        );

        update_stats(
            hit_m,
            m_hits,
            total_m_candidates,
            min_m_candidates,
            max_m_candidates,
            m_total_j,
            m_min_j,
            m_max_j,
            m_total_M,
            m_min_M,
            m_max_M
        );

        update_stats(
            hit_jm,
            jm_hits,
            total_jm_candidates,
            min_jm_candidates,
            max_jm_candidates,
            jm_total_j,
            jm_min_j,
            jm_max_j,
            jm_total_M,
            jm_min_M,
            jm_max_M
        );

        if (
            !hit_j.valid ||
            !hit_m.valid ||
            !hit_jm.valid
        ) {
            if (!printed_miss) {
                printed_miss = true;

                std::cout
                    << "\nFIRST_METHOD_MISS\n";

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
                    << "J_FIRST="
                    << hit_to_string(
                        hit_j
                    )
                    << "\n";

                std::cout
                    << "M_FIRST="
                    << hit_to_string(
                        hit_m
                    )
                    << "\n";

                std::cout
                    << "JM_FIRST="
                    << hit_to_string(
                        hit_jm
                    )
                    << "\n";
            }
        }

        if (
            !printed_first_j &&
            hit_j.valid
        ) {
            printed_first_j = true;

            std::cout
                << "\nFIRST_J_FIRST_HIT\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " N="
                << N
                << "\n";

            std::cout
                << hit_to_string(
                    hit_j
                )
                << "\n";
        }

        if (
            !printed_first_m &&
            hit_m.valid
        ) {
            printed_first_m = true;

            std::cout
                << "\nFIRST_M_FIRST_HIT\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " N="
                << N
                << "\n";

            std::cout
                << hit_to_string(
                    hit_m
                )
                << "\n";
        }

        if (
            !printed_first_jm &&
            hit_jm.valid
        ) {
            printed_first_jm = true;

            std::cout
                << "\nFIRST_JM_FIRST_HIT\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " N="
                << N
                << "\n";

            std::cout
                << hit_to_string(
                    hit_jm
                )
                << "\n";
        }
    }

    if (min_j_candidates == UINT64_MAX) {
        min_j_candidates = 0;
    }

    if (min_m_candidates == UINT64_MAX) {
        min_m_candidates = 0;
    }

    if (min_jm_candidates == UINT64_MAX) {
        min_jm_candidates = 0;
    }

    if (j_min_j == UINT64_MAX) {
        j_min_j = 0;
    }

    if (m_min_j == UINT64_MAX) {
        m_min_j = 0;
    }

    if (jm_min_j == UINT64_MAX) {
        jm_min_j = 0;
    }

    if (j_min_M == UINT64_MAX) {
        j_min_M = 0;
    }

    if (m_min_M == UINT64_MAX) {
        m_min_M = 0;
    }

    if (jm_min_M == UINT64_MAX) {
        jm_min_M = 0;
    }

    const double
        j_avg_candidates =
            j_hits == 0
                ? 0.0
                : static_cast<double>(
                      total_j_candidates
                  ) /
                  static_cast<double>(
                      j_hits
                  );

    const double
        m_avg_candidates =
            m_hits == 0
                ? 0.0
                : static_cast<double>(
                      total_m_candidates
                  ) /
                  static_cast<double>(
                      m_hits
                  );

    const double
        jm_avg_candidates =
            jm_hits == 0
                ? 0.0
                : static_cast<double>(
                      total_jm_candidates
                  ) /
                  static_cast<double>(
                      jm_hits
                  );

    const double
        j_avg_j =
            j_hits == 0
                ? 0.0
                : static_cast<double>(
                      j_total_j
                  ) /
                  static_cast<double>(
                      j_hits
                  );

    const double
        m_avg_j =
            m_hits == 0
                ? 0.0
                : static_cast<double>(
                      m_total_j
                  ) /
                  static_cast<double>(
                      m_hits
                  );

    const double
        jm_avg_j =
            jm_hits == 0
                ? 0.0
                : static_cast<double>(
                      jm_total_j
                  ) /
                  static_cast<double>(
                      jm_hits
                  );

    const double
        j_avg_M =
            j_hits == 0
                ? 0.0
                : static_cast<double>(
                      j_total_M
                  ) /
                  static_cast<double>(
                      j_hits
                  );

    const double
        m_avg_M =
            m_hits == 0
                ? 0.0
                : static_cast<double>(
                      m_total_M
                  ) /
                  static_cast<double>(
                      m_hits
                  );

    const double
        jm_avg_M =
            jm_hits == 0
                ? 0.0
                : static_cast<double>(
                      jm_total_M
                  ) /
                  static_cast<double>(
                      jm_hits
                  );

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "J_FIRST_HITS="
        << j_hits
        << "\n";

    std::cout
        << "M_FIRST_HITS="
        << m_hits
        << "\n";

    std::cout
        << "JM_FIRST_HITS="
        << jm_hits
        << "\n";

    std::cout
        << "J_FIRST_MISSES="
        << j_failures
        << "\n";

    std::cout
        << "M_FIRST_MISSES="
        << m_failures
        << "\n";

    std::cout
        << "JM_FIRST_MISSES="
        << jm_failures
        << "\n";

    std::cout
        << "\nJ_FIRST\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << j_avg_candidates
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_j_candidates
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_j_candidates
        << "\n";

    std::cout
        << "AVERAGE_J="
        << j_avg_j
        << "\n";

    std::cout
        << "MIN_J="
        << j_min_j
        << "\n";

    std::cout
        << "MAX_J="
        << j_max_j
        << "\n";

    std::cout
        << "AVERAGE_M="
        << j_avg_M
        << "\n";

    std::cout
        << "\nM_FIRST\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << m_avg_candidates
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_m_candidates
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_m_candidates
        << "\n";

    std::cout
        << "AVERAGE_J="
        << m_avg_j
        << "\n";

    std::cout
        << "MIN_J="
        << m_min_j
        << "\n";

    std::cout
        << "MAX_J="
        << m_max_j
        << "\n";

    std::cout
        << "AVERAGE_M="
        << m_avg_M
        << "\n";

    std::cout
        << "\nJM_FIRST\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << jm_avg_candidates
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << min_jm_candidates
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << max_jm_candidates
        << "\n";

    std::cout
        << "AVERAGE_J="
        << jm_avg_j
        << "\n";

    std::cout
        << "MIN_J="
        << jm_min_j
        << "\n";

    std::cout
        << "MAX_J="
        << jm_max_j
        << "\n";

    std::cout
        << "AVERAGE_M="
        << jm_avg_M
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
