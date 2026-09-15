#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>
#include <queue>
#include <unordered_set>

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

struct RankedPair {
    std::size_t index = 0;

    u64 d0 = 0;
    u64 M = 0;
};

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 candidate = 0;

    u64 d = 0;
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

    sieve[0] = false;
    sieve[1] = false;

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
            sq <=
            static_cast<i128>(n)
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
        Remove duplicate arithmetic progressions.
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

            return std::tie(
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

static RankedPair make_ranked_pair(
    const Pair& pair,
    std::size_t index,
    u64 s
) {
    RankedPair result;

    result.index = index;
    result.M = pair.M;

    const u64 base_mod =
        pair.base % pair.M;

    const u64 s_mod =
        s % pair.M;

    result.d0 =
        (s_mod + pair.M - base_mod) %
        pair.M;

    return result;
}

/*
    Strategy A:
    globally descending candidate.

    This is the best baseline from Experiment 461.
*/
static Hit search_candidate_descending(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    struct Node {
        u64 candidate = 0;
        u64 j = 0;
        std::size_t index = 0;
    };

    struct Compare {
        bool operator()(
            const Node& a,
            const Node& b
        ) const {
            if (a.candidate != b.candidate) {
                return a.candidate <
                       b.candidate;
            }

            return a.index >
                   b.index;
        }
    };

    Hit hit;

    std::priority_queue<
        Node,
        std::vector<Node>,
        Compare
    > heap;

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

        const u64 max_j =
            (
                s -
                pair.base
            ) /
            pair.M;

        if (max_j == 0) {
            continue;
        }

        Node node;

        node.j = max_j;
        node.index = i;

        node.candidate =
            pair.base +
            max_j * pair.M;

        heap.push(node);
    }

    std::unordered_set<u64> tested;

    tested.reserve(65536);

    while (!heap.empty()) {
        const Node node =
            heap.top();

        heap.pop();

        const Pair& pair =
            pairs[node.index];

        if (
            tested.insert(
                node.candidate
            ).second
        ) {
            ++hit.candidate_tests;

            const u64 g =
                std::gcd(
                    N,
                    node.candidate
                );

            if (
                g != 1 &&
                g != N
            ) {
                hit.valid = true;
                hit.factor = g;
                hit.candidate =
                    node.candidate;

                hit.d =
                    s -
                    node.candidate;

                hit.j = node.j;
                hit.M = pair.M;
                hit.pair_index =
                    node.index;

                return hit;
            }
        }

        if (node.j > 1) {
            const u64 next_j =
                node.j - 1;

            Node next;

            next.j = next_j;
            next.index =
                node.index;

            next.candidate =
                pair.base +
                next_j * pair.M;

            heap.push(next);
        }
    }

    return hit;
}

/*
    Strategy B:
    rank pairs by d0.

    For a pair,

        d = d0 + t*M.

    Search smaller d0 pairs first.
*/
static Hit search_d0_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::vector<RankedPair> ranked;

    ranked.reserve(
        pairs.size()
    );

    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        ranked.push_back(
            make_ranked_pair(
                pairs[i],
                i,
                s
            )
        );
    }

    std::sort(
        ranked.begin(),
        ranked.end(),
        [](
            const RankedPair& a,
            const RankedPair& b
        ) {
            if (a.d0 != b.d0) {
                return a.d0 <
                       b.d0;
            }

            if (a.M != b.M) {
                return a.M >
                       b.M;
            }

            return a.index <
                   b.index;
        }
    );

    std::unordered_set<u64> tested;

    tested.reserve(65536);

    for (
        std::size_t rank = 0;
        rank < ranked.size();
        ++rank
    ) {
        const RankedPair& rp =
            ranked[rank];

        const Pair& pair =
            pairs[rp.index];

        if (
            rp.d0 > s -
                pair.base
        ) {
            continue;
        }

        /*
            d values:

                d0,
                d0+M,
                d0+2M,...

            Candidate:

                p=s-d.
        */

        u64 d =
            rp.d0;

        while (
            d < s
        ) {
            const u64 candidate =
                s - d;

            if (
                candidate >= 2
            ) {
                if (
                    tested.insert(
                        candidate
                    ).second
                ) {
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
                        const u64 j =
                            candidate >=
                            pair.base
                                ? (
                                    candidate -
                                    pair.base
                                ) /
                                  pair.M
                                : 0;

                        hit.valid = true;
                        hit.factor = g;
                        hit.candidate =
                            candidate;

                        hit.d = d;
                        hit.j = j;
                        hit.M = pair.M;
                        hit.pair_index =
                            rp.index;

                        return hit;
                    }
                }
            }

            if (
                d >
                s -
                    pair.M
            ) {
                break;
            }

            d +=
                pair.M;
        }
    }

    return hit;
}

/*
    Strategy C:
    rank by d0/M.

    Small normalized first offset.
*/
static Hit search_d0_ratio_first(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::vector<RankedPair> ranked;

    ranked.reserve(
        pairs.size()
    );

    for (
        std::size_t i = 0;
        i < pairs.size();
        ++i
    ) {
        ranked.push_back(
            make_ranked_pair(
                pairs[i],
                i,
                s
            )
        );
    }

    std::sort(
        ranked.begin(),
        ranked.end(),
        [](
            const RankedPair& a,
            const RankedPair& b
        ) {
            const i128 lhs =
                static_cast<i128>(
                    a.d0
                ) *
                static_cast<i128>(
                    b.M
                );

            const i128 rhs =
                static_cast<i128>(
                    b.d0
                ) *
                static_cast<i128>(
                    a.M
                );

            if (lhs != rhs) {
                return lhs < rhs;
            }

            if (a.d0 != b.d0) {
                return a.d0 <
                       b.d0;
            }

            if (a.M != b.M) {
                return a.M >
                       b.M;
            }

            return a.index <
                   b.index;
        }
    );

    std::unordered_set<u64> tested;

    tested.reserve(65536);

    for (
        std::size_t rank = 0;
        rank < ranked.size();
        ++rank
    ) {
        const RankedPair& rp =
            ranked[rank];

        const Pair& pair =
            pairs[rp.index];

        u64 d =
            rp.d0;

        while (
            d < s
        ) {
            const u64 candidate =
                s - d;

            if (
                candidate >= 2
            ) {
                if (
                    tested.insert(
                        candidate
                    ).second
                ) {
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
                        const u64 j =
                            candidate >=
                            pair.base
                                ? (
                                    candidate -
                                    pair.base
                                ) /
                                  pair.M
                                : 0;

                        hit.valid = true;
                        hit.factor = g;
                        hit.candidate =
                            candidate;

                        hit.d = d;
                        hit.j = j;
                        hit.M = pair.M;
                        hit.pair_index =
                            rp.index;

                        return hit;
                    }
                }
            }

            if (
                d >
                s -
                    pair.M
            ) {
                break;
            }

            d +=
                pair.M;
        }
    }

    return hit;
}

static void print_hit(
    const char* label,
    const Hit& hit,
    u64 p,
    u64 q
) {
    std::cout
        << "\n"
        << label
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
        << "D="
        << hit.d
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
        << "CANDIDATE_TESTS="
        << hit.candidate_tests
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 463;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 500;

    constexpr int M_LIMIT = 7;
    constexpr u64 K_LIMIT = 1000;

    /*
        Same seed as Experiment 461.
    */
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
        [&]() {
            std::vector<u64> usable;

            for (u64 p : primes) {
                if (
                    p >=
                    static_cast<u64>(
                        PRIME_MIN
                    )
                ) {
                    usable.push_back(p);
                }
            }

            std::mt19937_64 rng(
                SEED
            );

            std::vector<
                std::pair<u64, u64>
            > result;

            result.reserve(
                CASE_COUNT
            );

            for (
                std::size_t i = 0;
                i < CASE_COUNT;
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

                while (a == b) {
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

                result.emplace_back(
                    usable[a],
                    usable[b]
                );
            }

            return result;
        }();

    u64 desc_hits = 0;
    u64 d0_hits = 0;
    u64 ratio_hits = 0;

    u64 desc_total = 0;
    u64 d0_total = 0;
    u64 ratio_total = 0;

    u64 desc_min = UINT64_MAX;
    u64 d0_min = UINT64_MAX;
    u64 ratio_min = UINT64_MAX;

    u64 desc_max = 0;
    u64 d0_max = 0;
    u64 ratio_max = 0;

    u64 desc_d_total = 0;
    u64 d0_d_total = 0;
    u64 ratio_d_total = 0;

    u64 desc_M_total = 0;
    u64 d0_M_total = 0;
    u64 ratio_M_total = 0;

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

        const Hit desc =
            search_candidate_descending(
                N,
                s,
                pairs
            );

        const Hit d0 =
            search_d0_first(
                N,
                s,
                pairs
            );

        const Hit ratio =
            search_d0_ratio_first(
                N,
                s,
                pairs
            );

        if (desc.valid) {
            ++desc_hits;
            desc_total +=
                desc.candidate_tests;
            desc_min =
                std::min(
                    desc_min,
                    desc.candidate_tests
                );
            desc_max =
                std::max(
                    desc_max,
                    desc.candidate_tests
                );
            desc_d_total += desc.d;
            desc_M_total += desc.M;
        }

        if (d0.valid) {
            ++d0_hits;
            d0_total +=
                d0.candidate_tests;
            d0_min =
                std::min(
                    d0_min,
                    d0.candidate_tests
                );
            d0_max =
                std::max(
                    d0_max,
                    d0.candidate_tests
                );
            d0_d_total += d0.d;
            d0_M_total += d0.M;
        }

        if (ratio.valid) {
            ++ratio_hits;
            ratio_total +=
                ratio.candidate_tests;
            ratio_min =
                std::min(
                    ratio_min,
                    ratio.candidate_tests
                );
            ratio_max =
                std::max(
                    ratio_max,
                    ratio.candidate_tests
                );
            ratio_d_total += ratio.d;
            ratio_M_total += ratio.M;
        }

        if (!printed_first) {
            printed_first = true;

            print_hit(
                "FIRST_DESCENDING",
                desc,
                p,
                q
            );

            print_hit(
                "FIRST_D0_FIRST",
                d0,
                p,
                q
            );

            print_hit(
                "FIRST_D0_RATIO_FIRST",
                ratio,
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

    if (desc_min == UINT64_MAX) {
        desc_min = 0;
    }

    if (d0_min == UINT64_MAX) {
        d0_min = 0;
    }

    if (ratio_min == UINT64_MAX) {
        ratio_min = 0;
    }

    const double desc_avg =
        desc_hits == 0
            ? 0.0
            : static_cast<double>(
                  desc_total
              ) /
              static_cast<double>(
                  desc_hits
              );

    const double d0_avg =
        d0_hits == 0
            ? 0.0
            : static_cast<double>(
                  d0_total
              ) /
              static_cast<double>(
                  d0_hits
              );

    const double ratio_avg =
        ratio_hits == 0
            ? 0.0
            : static_cast<double>(
                  ratio_total
              ) /
              static_cast<double>(
                  ratio_hits
              );

    const double desc_d_avg =
        desc_hits == 0
            ? 0.0
            : static_cast<double>(
                  desc_d_total
              ) /
              static_cast<double>(
                  desc_hits
              );

    const double d0_d_avg =
        d0_hits == 0
            ? 0.0
            : static_cast<double>(
                  d0_d_total
              ) /
              static_cast<double>(
                  d0_hits
              );

    const double ratio_d_avg =
        ratio_hits == 0
            ? 0.0
            : static_cast<double>(
                  ratio_d_total
              ) /
              static_cast<double>(
                  ratio_hits
              );

    const double desc_M_avg =
        desc_hits == 0
            ? 0.0
            : static_cast<double>(
                  desc_M_total
              ) /
              static_cast<double>(
                  desc_hits
              );

    const double d0_M_avg =
        d0_hits == 0
            ? 0.0
            : static_cast<double>(
                  d0_M_total
              ) /
              static_cast<double>(
                  d0_hits
              );

    const double ratio_M_avg =
        ratio_hits == 0
            ? 0.0
            : static_cast<double>(
                  ratio_M_total
              ) /
              static_cast<double>(
                  ratio_hits
              );

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "\nDESCENDING_BASELINE\n";

    std::cout
        << "HITS="
        << desc_hits
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << desc_avg
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << desc_min
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << desc_max
        << "\n";

    std::cout
        << "AVERAGE_D="
        << desc_d_avg
        << "\n";

    std::cout
        << "AVERAGE_M="
        << desc_M_avg
        << "\n";

    std::cout
        << "\nD0_FIRST\n";

    std::cout
        << "HITS="
        << d0_hits
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << d0_avg
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << d0_min
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << d0_max
        << "\n";

    std::cout
        << "AVERAGE_D="
        << d0_d_avg
        << "\n";

    std::cout
        << "AVERAGE_M="
        << d0_M_avg
        << "\n";

    std::cout
        << "\nD0_RATIO_FIRST\n";

    std::cout
        << "HITS="
        << ratio_hits
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << ratio_avg
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << ratio_min
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << ratio_max
        << "\n";

    std::cout
        << "AVERAGE_D="
        << ratio_d_avg
        << "\n";

    std::cout
        << "AVERAGE_M="
        << ratio_M_avg
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
