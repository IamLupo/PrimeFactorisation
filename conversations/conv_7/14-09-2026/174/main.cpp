#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 467;

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

struct PairScore {
    int j_score = 0;
    int m_score = 0;
    int total_score = 0;
};

std::vector<int> generate_primes(int limit) {
    std::vector<bool> composite(limit + 1, false);
    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<std::int64_t>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
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
        if (prime >= PRIME_MIN &&
            prime <= PRIME_LIMIT) {
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
        u64 p = static_cast<u64>(candidates[dist(rng)]);
        u64 q = static_cast<u64>(candidates[dist(rng)]);

        while (q == p) {
            q = static_cast<u64>(candidates[dist(rng)]);
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

    for (int k1 = 1; k1 <= K_LIMIT; ++k1) {
        for (int k2 = k1; k2 <= K_LIMIT; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            bool valid = false;

            /*
             * m2*k1 - m1*k2 = 1
             *
             * with 1 <= m1,m2 <= M_LIMIT.
             */
            for (int m1 = 1; m1 <= M_LIMIT; ++m1) {
                const int numerator =
                    1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const int m2 =
                    numerator / k1;

                if (m2 >= 1 &&
                    m2 <= M_LIMIT) {
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

        const u64 a = k1 + k2;
        const u64 m = k1 * k2;

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

PairScore score_pair(
    const Representation& fixed,
    const Representation& other
) {
    PairScore score;

    if (!fixed.found || !other.found) {
        return score;
    }

    /*
     * j-coupling:
     *
     * jp | jr*Mr
     * jr | jp*Mp
     */
    if ((other.j * other.m) % fixed.j == 0) {
        ++score.j_score;
    }

    if ((fixed.j * fixed.m) % other.j == 0) {
        ++score.j_score;
    }

    /*
     * M-coupling:
     *
     * Mp | (r-ar)
     * Mr | (p-ap)
     */
    if ((other.prime - other.a) % fixed.m == 0) {
        ++score.m_score;
    }

    if ((fixed.prime - fixed.a) % other.m == 0) {
        ++score.m_score;
    }

    score.total_score =
        score.j_score + score.m_score;

    return score;
}

void print_representation(
    const char* name,
    const Representation& rep
) {
    std::cout
        << name
        << "_PRIME="
        << rep.prime
        << '\n';

    std::cout
        << name
        << "_K1="
        << rep.k1
        << '\n';

    std::cout
        << name
        << "_K2="
        << rep.k2
        << '\n';

    std::cout
        << name
        << "_A="
        << rep.a
        << '\n';

    std::cout
        << name
        << "_M="
        << rep.m
        << '\n';

    std::cout
        << name
        << "_J="
        << rep.j
        << '\n';

    std::cout
        << name
        << "_CHECK="
        << rep.a + rep.j * rep.m
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x46720260915ULL
    );

    const std::vector<int> primes =
        generate_primes(PRIME_LIMIT);

    const std::vector<PrimeCase> cases =
        generate_cases(primes, rng);

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    /*
     * Cache the representation of every prime.
     * This avoids recomputing the 3647-pair search
     * for every control comparison.
     */
    std::vector<Representation> representation_cache(
        PRIME_LIMIT + 1
    );

    for (int prime : primes) {
        if (prime < PRIME_MIN) {
            continue;
        }

        representation_cache[prime] =
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
        << "CONTROL_COUNT="
        << CONTROL_COUNT
        << '\n';

    /*
     * True-factor statistics.
     */
    u64 true_found = 0;

    u64 true_j_score_0 = 0;
    u64 true_j_score_1 = 0;
    u64 true_j_score_2 = 0;

    u64 true_m_score_0 = 0;
    u64 true_m_score_1 = 0;
    u64 true_m_score_2 = 0;

    u64 true_total_score_0 = 0;
    u64 true_total_score_1 = 0;
    u64 true_total_score_2 = 0;
    u64 true_total_score_3 = 0;
    u64 true_total_score_4 = 0;

    /*
     * For the true q:
     *
     * number of controls with score >
     * number of controls with score =
     * number of controls with score <
     */
    u64 true_total_rank_sum = 0;
    u64 true_total_ties_sum = 0;
    u64 true_best_count = 0;
    u64 true_strict_best_count = 0;

    u64 true_highest_possible = 0;

    /*
     * Aggregate control score distribution.
     */
    u64 control_pairs = 0;

    u64 control_j_score_0 = 0;
    u64 control_j_score_1 = 0;
    u64 control_j_score_2 = 0;

    u64 control_m_score_0 = 0;
    u64 control_m_score_1 = 0;
    u64 control_m_score_2 = 0;

    u64 control_total_score_0 = 0;
    u64 control_total_score_1 = 0;
    u64 control_total_score_2 = 0;
    u64 control_total_score_3 = 0;
    u64 control_total_score_4 = 0;

    u64 controls_scored_above_true = 0;
    u64 controls_scored_equal_true = 0;

    /*
     * Track whether the true q beats all controls.
     */
    u64 true_beats_all = 0;
    u64 true_tied_best = 0;

    /*
     * Strong 4/4 structural match counts.
     */
    u64 true_score_4 = 0;
    u64 control_score_4 = 0;

    /*
     * First example.
     */
    bool first_example_printed = false;

    for (int case_index = 0;
         case_index < CASE_COUNT;
         ++case_index) {

        const PrimeCase& c =
            cases[case_index];

        const Representation& rp =
            representation_cache[c.p];

        const Representation& rq =
            representation_cache[c.q];

        if (!rp.found || !rq.found) {
            continue;
        }

        ++true_found;

        const PairScore true_score =
            score_pair(rp, rq);

        if (true_score.j_score == 0) {
            ++true_j_score_0;
        } else if (true_score.j_score == 1) {
            ++true_j_score_1;
        } else {
            ++true_j_score_2;
        }

        if (true_score.m_score == 0) {
            ++true_m_score_0;
        } else if (true_score.m_score == 1) {
            ++true_m_score_1;
        } else {
            ++true_m_score_2;
        }

        switch (true_score.total_score) {
            case 0:
                ++true_total_score_0;
                break;

            case 1:
                ++true_total_score_1;
                break;

            case 2:
                ++true_total_score_2;
                break;

            case 3:
                ++true_total_score_3;
                break;

            case 4:
                ++true_total_score_4;
                break;
        }

        if (true_score.total_score == 4) {
            ++true_highest_possible;
        }

        int controls_above = 0;
        int controls_equal = 0;
        int controls_below = 0;

        int controls_best = -1;

        std::vector<int> used_indices;
        used_indices.reserve(CONTROL_COUNT);

        std::uniform_int_distribution<std::size_t> dist(
            0,
            primes.size() - 1
        );

        while (
            static_cast<int>(used_indices.size()) <
            CONTROL_COUNT
        ) {
            const int candidate =
                primes[dist(rng)];

            if (candidate < PRIME_MIN) {
                continue;
            }

            if (
                static_cast<u64>(candidate) ==
                    c.p ||
                static_cast<u64>(candidate) ==
                    c.q
            ) {
                continue;
            }

            bool duplicate = false;

            for (int value : used_indices) {
                if (value == candidate) {
                    duplicate = true;
                    break;
                }
            }

            if (duplicate) {
                continue;
            }

            used_indices.push_back(candidate);
        }

        for (int candidate : used_indices) {
            const Representation& rr =
                representation_cache[candidate];

            if (!rr.found) {
                continue;
            }

            const PairScore score =
                score_pair(rp, rr);

            ++control_pairs;

            if (score.j_score == 0) {
                ++control_j_score_0;
            } else if (score.j_score == 1) {
                ++control_j_score_1;
            } else {
                ++control_j_score_2;
            }

            if (score.m_score == 0) {
                ++control_m_score_0;
            } else if (score.m_score == 1) {
                ++control_m_score_1;
            } else {
                ++control_m_score_2;
            }

            switch (score.total_score) {
                case 0:
                    ++control_total_score_0;
                    break;

                case 1:
                    ++control_total_score_1;
                    break;

                case 2:
                    ++control_total_score_2;
                    break;

                case 3:
                    ++control_total_score_3;
                    break;

                case 4:
                    ++control_total_score_4;
                    break;
            }

            if (score.total_score == 4) {
                ++control_score_4;
            }

            controls_best =
                std::max(
                    controls_best,
                    score.total_score
                );

            if (
                score.total_score >
                true_score.total_score
            ) {
                ++controls_above;
            } else if (
                score.total_score ==
                true_score.total_score
            ) {
                ++controls_equal;
            } else {
                ++controls_below;
            }
        }

        if (controls_above == 0) {
            ++true_beats_all;
        }

        if (
            controls_above == 0 &&
            controls_equal > 0
        ) {
            ++true_tied_best;
        }

        if (
            controls_above == 0 &&
            controls_equal == 0
        ) {
            ++true_strict_best_count;
        }

        /*
         * Rank convention:
         *
         * rank = 1 + number of controls
         *        with strictly larger score.
         */
        const u64 rank =
            1 + static_cast<u64>(controls_above);

        true_total_rank_sum += rank;
        true_total_ties_sum +=
            static_cast<u64>(controls_equal);

        if (
            controls_above == 0 &&
            controls_equal > 0
        ) {
            ++true_best_count;
        }

        if (!first_example_printed) {
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

            std::cout
                << "FIRST_TRUE_J_SCORE="
                << true_score.j_score
                << '\n';

            std::cout
                << "FIRST_TRUE_M_SCORE="
                << true_score.m_score
                << '\n';

            std::cout
                << "FIRST_TRUE_TOTAL_SCORE="
                << true_score.total_score
                << '\n';

            first_example_printed = true;
        }

        if ((case_index + 1) % 100 == 0) {
            std::cout
                << "PROGRESS="
                << (case_index + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const double avg_true_rank =
        true_found == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    true_total_rank_sum
                ) /
                static_cast<long double>(
                    true_found
                )
            );

    const double avg_true_ties =
        true_found == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    true_total_ties_sum
                ) /
                static_cast<long double>(
                    true_found
                )
            );

    const double average_percentile =
        true_found == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      true_total_rank_sum
                  ) /
                  static_cast<double>(
                      true_found * CONTROL_COUNT
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
        << "TRUE_FOUND="
        << true_found
        << '\n';

    std::cout
        << "TRUE_J_SCORE_0="
        << true_j_score_0
        << '\n';

    std::cout
        << "TRUE_J_SCORE_1="
        << true_j_score_1
        << '\n';

    std::cout
        << "TRUE_J_SCORE_2="
        << true_j_score_2
        << '\n';

    std::cout
        << "TRUE_M_SCORE_0="
        << true_m_score_0
        << '\n';

    std::cout
        << "TRUE_M_SCORE_1="
        << true_m_score_1
        << '\n';

    std::cout
        << "TRUE_M_SCORE_2="
        << true_m_score_2
        << '\n';

    std::cout
        << "TRUE_TOTAL_SCORE_0="
        << true_total_score_0
        << '\n';

    std::cout
        << "TRUE_TOTAL_SCORE_1="
        << true_total_score_1
        << '\n';

    std::cout
        << "TRUE_TOTAL_SCORE_2="
        << true_total_score_2
        << '\n';

    std::cout
        << "TRUE_TOTAL_SCORE_3="
        << true_total_score_3
        << '\n';

    std::cout
        << "TRUE_TOTAL_SCORE_4="
        << true_total_score_4
        << '\n';

    std::cout
        << "CONTROL_PAIRS_SCORED="
        << control_pairs
        << '\n';

    std::cout
        << "CONTROL_J_SCORE_0="
        << control_j_score_0
        << '\n';

    std::cout
        << "CONTROL_J_SCORE_1="
        << control_j_score_1
        << '\n';

    std::cout
        << "CONTROL_J_SCORE_2="
        << control_j_score_2
        << '\n';

    std::cout
        << "CONTROL_M_SCORE_0="
        << control_m_score_0
        << '\n';

    std::cout
        << "CONTROL_M_SCORE_1="
        << control_m_score_1
        << '\n';

    std::cout
        << "CONTROL_M_SCORE_2="
        << control_m_score_2
        << '\n';

    std::cout
        << "CONTROL_TOTAL_SCORE_0="
        << control_total_score_0
        << '\n';

    std::cout
        << "CONTROL_TOTAL_SCORE_1="
        << control_total_score_1
        << '\n';

    std::cout
        << "CONTROL_TOTAL_SCORE_2="
        << control_total_score_2
        << '\n';

    std::cout
        << "CONTROL_TOTAL_SCORE_3="
        << control_total_score_3
        << '\n';

    std::cout
        << "CONTROL_TOTAL_SCORE_4="
        << control_total_score_4
        << '\n';

    std::cout
        << "TRUE_STRICT_BEST="
        << true_strict_best_count
        << '\n';

    std::cout
        << "TRUE_BEATS_ALL_OR_TIES_BEST="
        << true_best_count
        << '\n';

    std::cout
        << "TRUE_TIED_BEST="
        << true_tied_best
        << '\n';

    std::cout
        << "AVG_TRUE_RANK="
        << avg_true_rank
        << '\n';

    std::cout
        << "AVG_TRUE_TIES="
        << avg_true_ties
        << '\n';

    std::cout
        << "AVG_TRUE_PERCENTILE="
        << average_percentile
        << '\n';

    std::cout
        << "TRUE_SCORE_4="
        << true_score_4
        << '\n';

    std::cout
        << "CONTROL_SCORE_4="
        << control_score_4
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
