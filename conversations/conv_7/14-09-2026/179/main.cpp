#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;

constexpr int EXPERIMENT = 472;

constexpr int PRIME_LIMIT = 100000;
constexpr int PRIME_MIN = 10000;

constexpr int CASE_COUNT = 500;
constexpr int CONTROL_COUNT = 100;

constexpr int K_LIMIT = 1000;
constexpr int M_LIMIT = 7;

struct PrimeCase {
    u64 p;
    u64 q;
    u64 n;
};

struct CRTPair {
    int k1;
    int k2;
};

struct Representation {
    bool found = false;

    u64 prime = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 a = 0;
    u64 m = 0;
    u64 j = 0;
};

struct ResidueFilter {
    bool valid = false;

    u64 inverse_a = 0;
    u64 target_residue = 0;

    u64 modulus = 0;
};

std::vector<int> generate_primes(int limit) {
    std::vector<bool> composite(
        limit + 1,
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (
            static_cast<std::int64_t>(i) * i <=
            limit
        ) {
            for (
                int j = i * i;
                j <= limit;
                j += i
            ) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

std::vector<PrimeCase> generate_cases(
    const std::vector<int>& primes,
    std::mt19937_64& rng
) {
    std::vector<int> candidates;

    for (int prime : primes) {
        if (
            prime >= PRIME_MIN &&
            prime <= PRIME_LIMIT
        ) {
            candidates.push_back(prime);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        candidates.size() - 1
    );

    std::vector<PrimeCase> cases;
    cases.reserve(CASE_COUNT);

    for (int i = 0; i < CASE_COUNT; ++i) {
        u64 p =
            static_cast<u64>(
                candidates[dist(rng)]
            );

        u64 q =
            static_cast<u64>(
                candidates[dist(rng)]
            );

        while (q == p) {
            q =
                static_cast<u64>(
                    candidates[dist(rng)]
                );
        }

        if (p > q) {
            std::swap(p, q);
        }

        PrimeCase c;
        c.p = p;
        c.q = q;
        c.n = p * q;

        cases.push_back(c);
    }

    return cases;
}

std::vector<CRTPair> build_determinant_one_pairs() {
    std::vector<CRTPair> pairs;

    for (
        int k1 = 1;
        k1 <= K_LIMIT;
        ++k1
    ) {
        for (
            int k2 = k1;
            k2 <= K_LIMIT;
            ++k2
        ) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            bool valid = false;

            /*
             * m2*k1 - m1*k2 = 1
             *
             * for 1 <= m1,m2 <= M_LIMIT.
             */
            for (
                int m1 = 1;
                m1 <= M_LIMIT;
                ++m1
            ) {
                const int numerator =
                    1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const int m2 =
                    numerator / k1;

                if (
                    m2 >= 1 &&
                    m2 <= M_LIMIT
                ) {
                    valid = true;
                    break;
                }
            }

            if (valid) {
                pairs.push_back({k1, k2});
            }
        }
    }

    return pairs;
}

Representation find_min_j_representation(
    u64 prime,
    const std::vector<CRTPair>& pairs
) {
    Representation best;
    best.prime = prime;

    for (const CRTPair& pair : pairs) {
        const u64 k1 =
            static_cast<u64>(pair.k1);

        const u64 k2 =
            static_cast<u64>(pair.k2);

        const u64 a =
            k1 + k2;

        const u64 m =
            k1 * k2;

        if (a > prime) {
            continue;
        }

        const u64 remainder =
            prime - a;

        if (remainder % m != 0) {
            continue;
        }

        const u64 j =
            remainder / m;

        bool take = false;

        if (!best.found) {
            take = true;
        } else if (j < best.j) {
            take = true;
        } else if (
            j == best.j &&
            m < best.m
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) <
                std::max(best.k1, best.k2)
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) ==
                std::max(best.k1, best.k2) &&
            k1 < best.k1
        ) {
            take = true;
        }

        if (!take) {
            continue;
        }

        best.found = true;
        best.k1 = k1;
        best.k2 = k2;
        best.a = a;
        best.m = m;
        best.j = j;
    }

    return best;
}

i64 extended_gcd(
    i64 a,
    i64 b,
    i64& x,
    i64& y
) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }

    i64 x1 = 0;
    i64 y1 = 0;

    const i64 g =
        extended_gcd(
            b,
            a % b,
            x1,
            y1
        );

    x = y1;
    y = x1 - (a / b) * y1;

    return g;
}

u64 modular_inverse(
    u64 a,
    u64 modulus,
    bool& valid
) {
    valid = false;

    if (modulus == 1) {
        valid = true;
        return 0;
    }

    i64 x = 0;
    i64 y = 0;

    const i64 g =
        extended_gcd(
            static_cast<i64>(a % modulus),
            static_cast<i64>(modulus),
            x,
            y
        );

    if (g != 1) {
        return 0;
    }

    i64 result =
        x % static_cast<i64>(modulus);

    if (result < 0) {
        result += static_cast<i64>(modulus);
    }

    valid = true;

    return static_cast<u64>(result);
}

ResidueFilter build_residue_filter(
    const PrimeCase& c,
    const Representation& rp
) {
    ResidueFilter filter;

    filter.modulus = rp.m;

    if (rp.m == 0) {
        return filter;
    }

    /*
     * gcd(a,M) must be 1 because
     *
     *   a = k1+k2
     *   M = k1*k2
     *   gcd(k1,k2)=1.
     */
    if (std::gcd(rp.a, rp.m) != 1) {
        return filter;
    }

    bool inverse_valid = false;

    filter.inverse_a =
        modular_inverse(
            rp.a,
            rp.m,
            inverse_valid
        );

    if (!inverse_valid) {
        return filter;
    }

    filter.target_residue =
        (
            (c.n % rp.m) *
            filter.inverse_a
        ) % rp.m;

    filter.valid = true;

    return filter;
}

bool matches_residue(
    u64 value,
    const ResidueFilter& filter
) {
    if (!filter.valid ||
        filter.modulus == 0) {
        return false;
    }

    return
        value % filter.modulus ==
        filter.target_residue;
}

std::vector<int> build_nearest_controls(
    const std::vector<int>& primes,
    u64 q,
    u64 p
) {
    std::vector<int> controls;
    controls.reserve(CONTROL_COUNT);

    auto lower =
        std::lower_bound(
            primes.begin(),
            primes.end(),
            static_cast<int>(q)
        );

    int right =
        static_cast<int>(
            lower - primes.begin()
        );

    int left =
        right - 1;

    while (
        static_cast<int>(
            controls.size()
        ) < CONTROL_COUNT &&
        (
            left >= 0 ||
            right <
                static_cast<int>(
                    primes.size()
                )
        )
    ) {
        bool take_left;

        if (left < 0) {
            take_left = false;
        } else if (
            right >=
            static_cast<int>(
                primes.size()
            )
        ) {
            take_left = true;
        } else {
            const u64 left_value =
                static_cast<u64>(
                    primes[left]
                );

            const u64 right_value =
                static_cast<u64>(
                    primes[right]
                );

            const u64 left_distance =
                q > left_value
                    ? q - left_value
                    : left_value - q;

            const u64 right_distance =
                q > right_value
                    ? q - right_value
                    : right_value - q;

            take_left =
                left_distance <=
                right_distance;
        }

        int candidate;

        if (take_left) {
            candidate = primes[left];
            --left;
        } else {
            candidate = primes[right];
            ++right;
        }

        if (
            static_cast<u64>(candidate) == p ||
            static_cast<u64>(candidate) == q
        ) {
            continue;
        }

        if (candidate < PRIME_MIN) {
            continue;
        }

        controls.push_back(candidate);
    }

    return controls;
}

void print_representation(
    const char* label,
    const Representation& rep
) {
    std::cout
        << label
        << "_PRIME="
        << rep.prime
        << '\n';

    std::cout
        << label
        << "_K1="
        << rep.k1
        << '\n';

    std::cout
        << label
        << "_K2="
        << rep.k2
        << '\n';

    std::cout
        << label
        << "_A="
        << rep.a
        << '\n';

    std::cout
        << label
        << "_M="
        << rep.m
        << '\n';

    std::cout
        << label
        << "_J="
        << rep.j
        << '\n';

    std::cout
        << label
        << "_CHECK="
        << rep.a +
           rep.j * rep.m
        << '\n';
}

void print_filter(
    const ResidueFilter& filter
) {
    std::cout
        << "FIRST_RESIDUE_VALID="
        << (filter.valid ? 1 : 0)
        << '\n';

    std::cout
        << "FIRST_MODULUS="
        << filter.modulus
        << '\n';

    std::cout
        << "FIRST_A_INVERSE="
        << filter.inverse_a
        << '\n';

    std::cout
        << "FIRST_TARGET_RESIDUE="
        << filter.target_residue
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47220260915ULL
    );

    const std::vector<int> primes =
        generate_primes(
            PRIME_LIMIT
        );

    const std::vector<PrimeCase> cases =
        generate_cases(
            primes,
            rng
        );

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    std::vector<Representation> cache(
        PRIME_LIMIT + 1
    );

    for (int prime : primes) {
        if (prime < PRIME_MIN) {
            continue;
        }

        cache[prime] =
            find_min_j_representation(
                static_cast<u64>(prime),
                pairs
            );
    }

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << '\n';

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << '\n';

    std::cout
        << "CONTROL_COUNT="
        << CONTROL_COUNT
        << '\n';

    u64 valid_filters = 0;
    u64 gcd_a_m_failures = 0;
    u64 inverse_failures = 0;

    /*
     * True q must always pass.
     */
    u64 true_hits = 0;

    /*
     * Nearby controls.
     */
    u64 control_pairs = 0;
    u64 control_hits = 0;

    /*
     * All primes in PRIME_MIN..PRIME_LIMIT.
     */
    u64 all_pool_tests = 0;
    u64 all_pool_hits = 0;

    /*
     * How many primes remain after the residue filter.
     */
    u64 minimum_candidates = UINT64_MAX;
    u64 maximum_candidates = 0;
    long double total_candidates = 0.0L;

    /*
     * Rank of the true q among all matching primes,
     * ordered by numerical distance to q.
     */
    u64 all_rank_sum = 0;
    u64 all_matching_ties_sum = 0;

    /*
     * Is q the only prime satisfying the residue
     * among the entire configured prime pool?
     */
    u64 true_unique = 0;
    u64 true_not_unique = 0;

    /*
     * Local control statistics.
     */
    u64 nearby_zero_hits = 0;
    u64 nearby_one_hit = 0;
    u64 nearby_two_plus_hits = 0;

    /*
     * Modulus size statistics.
     */
    u64 min_modulus = UINT64_MAX;
    u64 max_modulus = 0;

    long double sum_modulus = 0.0L;

    /*
     * First example.
     */
    bool first_example = false;

    for (
        int case_index = 0;
        case_index < CASE_COUNT;
        ++case_index
    ) {
        const PrimeCase& c =
            cases[case_index];

        const Representation& rp =
            cache[c.p];

        const Representation& rq =
            cache[c.q];

        if (!rp.found || !rq.found) {
            continue;
        }

        const u64 gcd_a_m =
            std::gcd(
                rp.a,
                rp.m
            );

        if (gcd_a_m != 1) {
            ++gcd_a_m_failures;
            continue;
        }

        const ResidueFilter filter =
            build_residue_filter(
                c,
                rp
            );

        if (!filter.valid) {
            ++inverse_failures;
            continue;
        }

        ++valid_filters;

        if (
            matches_residue(
                c.q,
                filter
            )
        ) {
            ++true_hits;
        }

        min_modulus =
            std::min(
                min_modulus,
                filter.modulus
            );

        max_modulus =
            std::max(
                max_modulus,
                filter.modulus
            );

        sum_modulus +=
            static_cast<long double>(
                filter.modulus
            );

        /*
         * Count every prime in the complete
         * search pool that satisfies the residue.
         */
        std::vector<int> matching_primes;
        matching_primes.reserve(64);

        for (int prime : primes) {
            if (
                prime < PRIME_MIN ||
                prime > PRIME_LIMIT
            ) {
                continue;
            }

            ++all_pool_tests;

            if (
                matches_residue(
                    static_cast<u64>(prime),
                    filter
                )
            ) {
                ++all_pool_hits;
                matching_primes.push_back(
                    prime
                );
            }
        }

        const u64 candidate_count =
            matching_primes.size();

        minimum_candidates =
            std::min(
                minimum_candidates,
                candidate_count
            );

        maximum_candidates =
            std::max(
                maximum_candidates,
                candidate_count
            );

        total_candidates +=
            static_cast<long double>(
                candidate_count
            );

        /*
         * Sort matching primes by absolute
         * distance from the real q.
         */
        std::sort(
            matching_primes.begin(),
            matching_primes.end(),
            [q = c.q](
                int x,
                int y
            ) {
                const u64 dx =
                    q > static_cast<u64>(x)
                        ? q -
                          static_cast<u64>(x)
                        : static_cast<u64>(x) -
                          q;

                const u64 dy =
                    q > static_cast<u64>(y)
                        ? q -
                          static_cast<u64>(y)
                        : static_cast<u64>(y) -
                          q;

                if (dx != dy) {
                    return dx < dy;
                }

                return x < y;
            }
        );

        auto q_position =
            std::find(
                matching_primes.begin(),
                matching_primes.end(),
                static_cast<int>(c.q)
            );

        if (
            q_position !=
            matching_primes.end()
        ) {
            const u64 rank =
                1 +
                static_cast<u64>(
                    q_position -
                    matching_primes.begin()
                );

            all_rank_sum += rank;

            u64 same_distance = 0;

            const u64 q_distance =
                0;

            for (int prime :
                 matching_primes) {

                const u64 distance =
                    prime > static_cast<int>(c.q)
                        ? static_cast<u64>(prime) -
                          c.q
                        : c.q -
                          static_cast<u64>(prime);

                if (distance == q_distance) {
                    ++same_distance;
                }
            }

            all_matching_ties_sum +=
                same_distance - 1;
        }

        if (candidate_count == 2) {
            /*
             * In the pool there are exactly two
             * matching primes. Since q is one of
             * them, this is especially interesting.
             */
        }

        if (candidate_count == 1) {
            ++true_unique;
        } else {
            ++true_not_unique;
        }

        /*
         * Nearby controls.
         */
        const std::vector<int> controls =
            build_nearest_controls(
                primes,
                c.q,
                c.p
            );

        u64 nearby_hits_for_case = 0;

        for (int candidate : controls) {
            ++control_pairs;

            if (
                matches_residue(
                    static_cast<u64>(candidate),
                    filter
                )
            ) {
                ++control_hits;
                ++nearby_hits_for_case;
            }
        }

        if (nearby_hits_for_case == 0) {
            ++nearby_zero_hits;
        } else if (
            nearby_hits_for_case == 1
        ) {
            ++nearby_one_hit;
        } else {
            ++nearby_two_plus_hits;
        }

        if (!first_example) {
            std::cout
                << "FIRST_N="
                << c.n
                << '\n';

            print_representation(
                "FIRST_P",
                rp
            );

            print_representation(
                "FIRST_Q",
                rq
            );

            print_filter(filter);

            std::cout
                << "FIRST_TRUE_Q_MATCH="
                << (
                    matches_residue(
                        c.q,
                        filter
                    )
                        ? 1
                        : 0
                )
                << '\n';

            std::cout
                << "FIRST_ALL_POOL_MATCHING_COUNT="
                << candidate_count
                << '\n';

            std::cout
                << "FIRST_NEARBY_MATCHING_COUNT="
                << nearby_hits_for_case
                << '\n';

            first_example = true;
        }

        if (
            (case_index + 1) % 100 ==
            0
        ) {
            std::cout
                << "PROGRESS="
                << (case_index + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const double avg_modulus =
        valid_filters == 0
            ? 0.0
            : static_cast<double>(
                sum_modulus /
                static_cast<long double>(
                    valid_filters
                )
            );

    const double avg_candidates =
        valid_filters == 0
            ? 0.0
            : static_cast<double>(
                total_candidates /
                static_cast<long double>(
                    valid_filters
                )
            );

    const double true_pool_hit_rate =
        valid_filters == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  true_hits
              ) /
              static_cast<double>(
                  valid_filters
              );

    const double nearby_hit_rate =
        control_pairs == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  control_hits
              ) /
              static_cast<double>(
                  control_pairs
              );

    const double pool_density =
        all_pool_tests == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  all_pool_hits
              ) /
              static_cast<double>(
                  all_pool_tests
              );

    const double avg_true_rank =
        valid_filters == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    all_rank_sum
                ) /
                static_cast<long double>(
                    valid_filters
                )
            );

    const double avg_matching_ties =
        valid_filters == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    all_matching_ties_sum
                ) /
                static_cast<long double>(
                    valid_filters
                )
            );

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "CONTROL_COUNT="
        << CONTROL_COUNT
        << '\n';

    std::cout
        << "VALID_FILTERS="
        << valid_filters
        << '\n';

    std::cout
        << "GCD_A_M_FAILURES="
        << gcd_a_m_failures
        << '\n';

    std::cout
        << "INVERSE_FAILURES="
        << inverse_failures
        << '\n';

    std::cout
        << "TRUE_Q_MATCHES="
        << true_hits
        << '\n';

    std::cout
        << "TRUE_Q_MATCH_RATE="
        << true_pool_hit_rate
        << '\n';

    std::cout
        << "NEARBY_CONTROL_PAIRS="
        << control_pairs
        << '\n';

    std::cout
        << "NEARBY_CONTROL_HITS="
        << control_hits
        << '\n';

    std::cout
        << "NEARBY_CONTROL_HIT_RATE="
        << nearby_hit_rate
        << '\n';

    std::cout
        << "NEARBY_ZERO_HITS="
        << nearby_zero_hits
        << '\n';

    std::cout
        << "NEARBY_ONE_HIT="
        << nearby_one_hit
        << '\n';

    std::cout
        << "NEARBY_TWO_PLUS_HITS="
        << nearby_two_plus_hits
        << '\n';

    std::cout
        << "ALL_POOL_TESTS="
        << all_pool_tests
        << '\n';

    std::cout
        << "ALL_POOL_HITS="
        << all_pool_hits
        << '\n';

    std::cout
        << "ALL_POOL_HIT_DENSITY_PERCENT="
        << pool_density
        << '\n';

    std::cout
        << "MIN_MODULUS="
        << (
            min_modulus == UINT64_MAX
                ? 0
                : min_modulus
        )
        << '\n';

    std::cout
        << "MAX_MODULUS="
        << max_modulus
        << '\n';

    std::cout
        << "AVG_MODULUS="
        << avg_modulus
        << '\n';

    std::cout
        << "MIN_MATCHING_PRIMES="
        << (
            minimum_candidates ==
                UINT64_MAX
                ? 0
                : minimum_candidates
        )
        << '\n';

    std::cout
        << "MAX_MATCHING_PRIMES="
        << maximum_candidates
        << '\n';

    std::cout
        << "AVG_MATCHING_PRIMES="
        << avg_candidates
        << '\n';

    std::cout
        << "TRUE_Q_UNIQUE_IN_POOL="
        << true_unique
        << '\n';

    std::cout
        << "TRUE_Q_NOT_UNIQUE_IN_POOL="
        << true_not_unique
        << '\n';

    std::cout
        << "AVG_TRUE_Q_DISTANCE_RANK="
        << avg_true_rank
        << '\n';

    std::cout
        << "AVG_MATCHING_TIES_AT_Q="
        << avg_matching_ties
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';
}

int main() {
    main_experiment();
    return 0;
}
