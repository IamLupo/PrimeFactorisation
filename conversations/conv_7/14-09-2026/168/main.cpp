#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <string>
#include <unordered_set>
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

struct Node {
    u64 candidate = 0;
    u64 j = 0;
    std::size_t pair_index = 0;
};

struct MinCompare {
    bool operator()(
        const Node& a,
        const Node& b
    ) const {
        if (a.candidate != b.candidate) {
            return a.candidate > b.candidate;
        }

        if (a.pair_index != b.pair_index) {
            return a.pair_index > b.pair_index;
        }

        return a.j > b.j;
    }
};

struct MaxCompare {
    bool operator()(
        const Node& a,
        const Node& b
    ) const {
        if (a.candidate != b.candidate) {
            return a.candidate < b.candidate;
        }

        if (a.pair_index != b.pair_index) {
            return a.pair_index > b.pair_index;
        }

        return a.j < b.j;
    }
};

struct Hit {
    bool valid = false;

    u64 factor = 0;
    u64 candidate = 0;

    u64 j = 0;
    u64 M = 0;

    u64 candidate_tests = 0;
    u64 pair_states = 0;

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
                    static_cast<i128>(m1) != 0
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
                return a.base < b.base;
            }

            if (a.M != b.M) {
                return a.M < b.M;
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
            unique_pairs.push_back(pair);
        }
    }

    return unique_pairs;
}

static Hit search_ascending(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::priority_queue<
        Node,
        std::vector<Node>,
        MinCompare
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

        Node node;

        node.j = 1;
        node.pair_index = i;
        node.candidate =
            pair.base +
            pair.M;

        heap.push(node);
    }

    u64 last_candidate = 0;
    bool have_last = false;

    while (!heap.empty()) {
        const Node node =
            heap.top();

        heap.pop();

        const Pair& pair =
            pairs[node.pair_index];

        if (
            !have_last ||
            node.candidate !=
                last_candidate
        ) {
            have_last = true;
            last_candidate =
                node.candidate;

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

                hit.j = node.j;
                hit.M = pair.M;
                hit.pair_index =
                    node.pair_index;

                return hit;
            }
        }

        const u64 next_j =
            node.j + 1;

        const i128 next_candidate_i =
            static_cast<i128>(
                pair.base
            ) +
            static_cast<i128>(
                next_j
            ) *
            static_cast<i128>(
                pair.M
            );

        if (
            next_candidate_i <=
            static_cast<i128>(s)
        ) {
            Node next;

            next.j = next_j;
            next.pair_index =
                node.pair_index;

            next.candidate =
                static_cast<u64>(
                    next_candidate_i
                );

            heap.push(next);
        }
    }

    return hit;
}

static Hit search_descending(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::priority_queue<
        Node,
        std::vector<Node>,
        MaxCompare
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
        node.pair_index = i;

        node.candidate =
            pair.base +
            max_j * pair.M;

        heap.push(node);
    }

    u64 last_candidate = 0;
    bool have_last = false;

    while (!heap.empty()) {
        const Node node =
            heap.top();

        heap.pop();

        const Pair& pair =
            pairs[node.pair_index];

        if (
            !have_last ||
            node.candidate !=
                last_candidate
        ) {
            have_last = true;
            last_candidate =
                node.candidate;

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

                hit.j = node.j;
                hit.M = pair.M;
                hit.pair_index =
                    node.pair_index;

                return hit;
            }
        }

        if (node.j > 1) {
            const u64 next_j =
                node.j - 1;

            const u64 next_candidate =
                pair.base +
                next_j * pair.M;

            Node next;

            next.j = next_j;
            next.pair_index =
                node.pair_index;

            next.candidate =
                next_candidate;

            heap.push(next);
        }
    }

    return hit;
}

static Hit search_alternating(
    u64 N,
    u64 s,
    const std::vector<Pair>& pairs
) {
    Hit hit;

    std::priority_queue<
        Node,
        std::vector<Node>,
        MinCompare
    > min_heap;

    std::priority_queue<
        Node,
        std::vector<Node>,
        MaxCompare
    > max_heap;

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

        Node low;

        low.j = 1;
        low.pair_index = i;
        low.candidate =
            pair.base +
            pair.M;

        min_heap.push(low);

        const u64 max_j =
            (
                s -
                pair.base
            ) /
            pair.M;

        Node high;

        high.j = max_j;
        high.pair_index = i;
        high.candidate =
            pair.base +
            max_j * pair.M;

        max_heap.push(high);
    }

    std::unordered_set<u64> tested;

    tested.reserve(
        65536
    );

    bool choose_low = true;

    while (
        !min_heap.empty() ||
        !max_heap.empty()
    ) {
        Node node;
        bool got_node = false;

        if (
            choose_low &&
            !min_heap.empty()
        ) {
            node =
                min_heap.top();

            min_heap.pop();

            got_node = true;
        } else if (
            !choose_low &&
            !max_heap.empty()
        ) {
            node =
                max_heap.top();

            max_heap.pop();

            got_node = true;
        } else if (
            !min_heap.empty()
        ) {
            node =
                min_heap.top();

            min_heap.pop();

            got_node = true;
        } else {
            node =
                max_heap.top();

            max_heap.pop();

            got_node = true;
        }

        if (!got_node) {
            break;
        }

        const Pair& pair =
            pairs[node.pair_index];

        const bool new_candidate =
            tested.insert(
                node.candidate
            ).second;

        if (new_candidate) {
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

                hit.j = node.j;
                hit.M = pair.M;
                hit.pair_index =
                    node.pair_index;

                return hit;
            }
        }

        /*
            Advance the sequence that produced
            this node.

            The two heaps contain independent
            traversals, so duplicate values are
            harmless and are removed by tested.
        */

        if (choose_low) {
            const u64 next_j =
                node.j + 1;

            const i128 next_candidate_i =
                static_cast<i128>(
                    pair.base
                ) +
                static_cast<i128>(
                    next_j
                ) *
                static_cast<i128>(
                    pair.M
                );

            if (
                next_candidate_i <=
                static_cast<i128>(s)
            ) {
                Node next;

                next.j = next_j;
                next.pair_index =
                    node.pair_index;

                next.candidate =
                    static_cast<u64>(
                        next_candidate_i
                    );

                min_heap.push(next);
            }
        } else {
            if (node.j > 1) {
                const u64 next_j =
                    node.j - 1;

                Node next;

                next.j = next_j;
                next.pair_index =
                    node.pair_index;

                next.candidate =
                    pair.base +
                    next_j * pair.M;

                max_heap.push(next);
            }
        }

        choose_low =
            !choose_low;
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
    constexpr int EXPERIMENT = 461;

    constexpr int PRIME_LIMIT = 100000;
    constexpr int PRIME_MIN = 10000;

    constexpr std::size_t CASE_COUNT = 500;

    constexpr int M_LIMIT = 7;
    constexpr u64 K_LIMIT = 1000;

    /*
        Same seed as Experiment 460 so the
        semiprime population is directly comparable.
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
        generate_cases(
            primes,
            PRIME_MIN,
            CASE_COUNT,
            SEED
        );

    u64 ascending_hits = 0;
    u64 descending_hits = 0;
    u64 alternating_hits = 0;

    u64 ascending_total = 0;
    u64 descending_total = 0;
    u64 alternating_total = 0;

    u64 ascending_min = UINT64_MAX;
    u64 descending_min = UINT64_MAX;
    u64 alternating_min = UINT64_MAX;

    u64 ascending_max = 0;
    u64 descending_max = 0;
    u64 alternating_max = 0;

    u64 ascending_j_total = 0;
    u64 descending_j_total = 0;
    u64 alternating_j_total = 0;

    u64 ascending_j_min = UINT64_MAX;
    u64 descending_j_min = UINT64_MAX;
    u64 alternating_j_min = UINT64_MAX;

    u64 ascending_j_max = 0;
    u64 descending_j_max = 0;
    u64 alternating_j_max = 0;

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

        const Hit ascending =
            search_ascending(
                N,
                s,
                pairs
            );

        const Hit descending =
            search_descending(
                N,
                s,
                pairs
            );

        const Hit alternating =
            search_alternating(
                N,
                s,
                pairs
            );

        if (ascending.valid) {
            ++ascending_hits;

            ascending_total +=
                ascending.candidate_tests;

            ascending_min =
                std::min(
                    ascending_min,
                    ascending.candidate_tests
                );

            ascending_max =
                std::max(
                    ascending_max,
                    ascending.candidate_tests
                );

            ascending_j_total +=
                ascending.j;

            ascending_j_min =
                std::min(
                    ascending_j_min,
                    ascending.j
                );

            ascending_j_max =
                std::max(
                    ascending_j_max,
                    ascending.j
                );
        }

        if (descending.valid) {
            ++descending_hits;

            descending_total +=
                descending.candidate_tests;

            descending_min =
                std::min(
                    descending_min,
                    descending.candidate_tests
                );

            descending_max =
                std::max(
                    descending_max,
                    descending.candidate_tests
                );

            descending_j_total +=
                descending.j;

            descending_j_min =
                std::min(
                    descending_j_min,
                    descending.j
                );

            descending_j_max =
                std::max(
                    descending_j_max,
                    descending.j
                );
        }

        if (alternating.valid) {
            ++alternating_hits;

            alternating_total +=
                alternating.candidate_tests;

            alternating_min =
                std::min(
                    alternating_min,
                    alternating.candidate_tests
                );

            alternating_max =
                std::max(
                    alternating_max,
                    alternating.candidate_tests
                );

            alternating_j_total +=
                alternating.j;

            alternating_j_min =
                std::min(
                    alternating_j_min,
                    alternating.j
                );

            alternating_j_max =
                std::max(
                    alternating_j_max,
                    alternating.j
                );
        }

        if (!printed_first) {
            printed_first = true;

            print_hit(
                "FIRST_ASCENDING_HIT",
                ascending,
                p,
                q
            );

            print_hit(
                "FIRST_DESCENDING_HIT",
                descending,
                p,
                q
            );

            print_hit(
                "FIRST_ALTERNATING_HIT",
                alternating,
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

    if (ascending_min == UINT64_MAX) {
        ascending_min = 0;
    }

    if (descending_min == UINT64_MAX) {
        descending_min = 0;
    }

    if (alternating_min == UINT64_MAX) {
        alternating_min = 0;
    }

    if (ascending_j_min == UINT64_MAX) {
        ascending_j_min = 0;
    }

    if (descending_j_min == UINT64_MAX) {
        descending_j_min = 0;
    }

    if (alternating_j_min == UINT64_MAX) {
        alternating_j_min = 0;
    }

    const double ascending_avg =
        ascending_hits == 0
            ? 0.0
            : static_cast<double>(
                  ascending_total
              ) /
              static_cast<double>(
                  ascending_hits
              );

    const double descending_avg =
        descending_hits == 0
            ? 0.0
            : static_cast<double>(
                  descending_total
              ) /
              static_cast<double>(
                  descending_hits
              );

    const double alternating_avg =
        alternating_hits == 0
            ? 0.0
            : static_cast<double>(
                  alternating_total
              ) /
              static_cast<double>(
                  alternating_hits
              );

    const double ascending_j_avg =
        ascending_hits == 0
            ? 0.0
            : static_cast<double>(
                  ascending_j_total
              ) /
              static_cast<double>(
                  ascending_hits
              );

    const double descending_j_avg =
        descending_hits == 0
            ? 0.0
            : static_cast<double>(
                  descending_j_total
              ) /
              static_cast<double>(
                  descending_hits
              );

    const double alternating_j_avg =
        alternating_hits == 0
            ? 0.0
            : static_cast<double>(
                  alternating_j_total
              ) /
              static_cast<double>(
                  alternating_hits
              );

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "\nASCENDING_CANDIDATE\n";

    std::cout
        << "HITS="
        << ascending_hits
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << ascending_avg
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << ascending_min
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << ascending_max
        << "\n";

    std::cout
        << "AVERAGE_J="
        << ascending_j_avg
        << "\n";

    std::cout
        << "MIN_J="
        << ascending_j_min
        << "\n";

    std::cout
        << "MAX_J="
        << ascending_j_max
        << "\n";

    std::cout
        << "\nDESCENDING_CANDIDATE\n";

    std::cout
        << "HITS="
        << descending_hits
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << descending_avg
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << descending_min
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << descending_max
        << "\n";

    std::cout
        << "AVERAGE_J="
        << descending_j_avg
        << "\n";

    std::cout
        << "MIN_J="
        << descending_j_min
        << "\n";

    std::cout
        << "MAX_J="
        << descending_j_max
        << "\n";

    std::cout
        << "\nALTERNATING_CANDIDATE\n";

    std::cout
        << "HITS="
        << alternating_hits
        << "\n";

    std::cout
        << "AVERAGE_CANDIDATE_TESTS="
        << alternating_avg
        << "\n";

    std::cout
        << "MIN_CANDIDATE_TESTS="
        << alternating_min
        << "\n";

    std::cout
        << "MAX_CANDIDATE_TESTS="
        << alternating_max
        << "\n";

    std::cout
        << "AVERAGE_J="
        << alternating_j_avg
        << "\n";

    std::cout
        << "MIN_J="
        << alternating_j_min
        << "\n";

    std::cout
        << "MAX_J="
        << alternating_j_max
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
