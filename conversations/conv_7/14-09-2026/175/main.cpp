#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 468;

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
                std::max(
                    best.k1,
                    best.k2
                )
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) ==
                std::max(
                    best.k1,
                    best.k2
                ) &&
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

    if (
        !fixed.found ||
        !other.found
    ) {
        return score;
    }

    /*
     * J-score:
     *
     *   jp | jr*Mr
     *   jr | jp*Mp
     */
    if (
        other.j != 0 &&
        (other.j * other.m) %
            fixed.j == 0
    ) {
        ++score.j_score;
    }

    if (
        fixed.j != 0 &&
        (fixed.j * fixed.m) %
            other.j == 0
    ) {
        ++score.j_score;
    }

    /*
     * M-score:
     *
     *   Mp | (r-ar)
     *   Mr | (p-ap)
     */
    if (
        (other.prime - other.a) %
            fixed.m == 0
    ) {
        ++score.m_score;
    }

    if (
        (fixed.prime - fixed.a) %
            other.m == 0
    ) {
        ++score.m_score;
    }

    score.total_score =
        score.j_score +
        score.m_score;

    return score;
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

    int left = right - 1;

    while (
        static_cast<int>(
            controls.size()
        ) < CONTROL_COUNT &&
        (left >= 0 ||
         right < static_cast<int>(primes.size()))
    ) {
        bool take_left = false;

        if (left < 0) {
            take_left = false;
        } else if (
            right >=
            static_cast<int>(primes.size())
        ) {
            take_left = true;
        } else {
            const u64 left_distance =
                q >
                    static_cast<u64>(
                        primes[left]
                    )
                    ? q -
                      static_cast<u64>(
                          primes[left]
                      )
                    : static_cast<u64>(
                          primes[left]
                      ) - q;

            const u64 right_distance =
                q >
                    static_cast<u64>(
                        primes[right]
                    )
                    ? q -
                      static_cast<u64>(
                          primes[right]
                      )
                    : static_cast<u64>(
                          primes[right]
                      ) - q;

            take_left =
                left_distance <=
                right_distance;
        }

        int candidate;

        if (take_left) {
            candidate =
                primes[left];
            --left;
        } else {
            candidate =
                primes[right];
            ++right;
        }

        if (
            static_cast<u64>(candidate) == q ||
            static_cast<u64>(candidate) == p
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

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x46820260915ULL
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

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << '\n';

    std::cout
        << "CONTROL_COUNT="
        << CONTROL_COUNT
        << '\n';

    /*
     * Cache every representation once.
     */
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

    u64 true_found = 0;
    u64 control_pairs = 0;

    /*
     * True J-score distribution.
     */
    u64 true_j0 = 0;
    u64 true_j1 = 0;
    u64 true_j2 = 0;

    /*
     * Control J-score distribution.
     */
    u64 control_j0 = 0;
    u64 control_j1 = 0;
    u64 control_j2 = 0;

    /*
     * True total-score distribution.
     */
    u64 true_total0 = 0;
    u64 true_total1 = 0;
    u64 true_total2 = 0;
    u64 true_total3 = 0;
    u64 true_total4 = 0;

    /*
     * Control total-score distribution.
     */
    u64 control_total0 = 0;
    u64 control_total1 = 0;
    u64 control_total2 = 0;
    u64 control_total3 = 0;
    u64 control_total4 = 0;

    /*
     * For every case:
     *
     * rank = 1 + controls with
     *        strictly larger J-score.
     */
    u64 j_rank_sum = 0;
    u64 total_rank_sum = 0;

    u64 j_tie_sum = 0;
    u64 total_tie_sum = 0;

    u64 j_strict_best = 0;
    u64 j_tied_best = 0;

    u64 total_strict_best = 0;
    u64 total_tied_best = 0;

    /*
     * How often does true q beat at least
     * 75%, 90%, 95%, 99% of matched controls?
     */
    u64 j_percentile75 = 0;
    u64 j_percentile90 = 0;
    u64 j_percentile95 = 0;
    u64 j_percentile99 = 0;

    u64 total_percentile75 = 0;
    u64 total_percentile90 = 0;
    u64 total_percentile95 = 0;
    u64 total_percentile99 = 0;

    /*
     * Average absolute distance of controls
     * from q.
     */
    long double sum_control_distance = 0.0L;
    u64 max_control_distance = 0;

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

        if (
            !rp.found ||
            !rq.found
        ) {
            continue;
        }

        ++true_found;

        const PairScore true_score =
            score_pair(
                rp,
                rq
            );

        switch (true_score.j_score) {
            case 0:
                ++true_j0;
                break;

            case 1:
                ++true_j1;
                break;

            case 2:
                ++true_j2;
                break;
        }

        switch (true_score.total_score) {
            case 0:
                ++true_total0;
                break;

            case 1:
                ++true_total1;
                break;

            case 2:
                ++true_total2;
                break;

            case 3:
                ++true_total3;
                break;

            case 4:
                ++true_total4;
                break;
        }

        const std::vector<int> controls =
            build_nearest_controls(
                primes,
                c.q,
                c.p
            );

        int j_above = 0;
        int j_equal = 0;

        int total_above = 0;
        int total_equal = 0;

        for (int candidate : controls) {
            const Representation& rr =
                cache[candidate];

            if (!rr.found) {
                continue;
            }

            const PairScore score =
                score_pair(
                    rp,
                    rr
                );

            ++control_pairs;

            switch (score.j_score) {
                case 0:
                    ++control_j0;
                    break;

                case 1:
                    ++control_j1;
                    break;

                case 2:
                    ++control_j2;
                    break;
            }

            switch (score.total_score) {
                case 0:
                    ++control_total0;
                    break;

                case 1:
                    ++control_total1;
                    break;

                case 2:
                    ++control_total2;
                    break;

                case 3:
                    ++control_total3;
                    break;

                case 4:
                    ++control_total4;
                    break;
            }

            if (
                score.j_score >
                true_score.j_score
            ) {
                ++j_above;
            } else if (
                score.j_score ==
                true_score.j_score
            ) {
                ++j_equal;
            }

            if (
                score.total_score >
                true_score.total_score
            ) {
                ++total_above;
            } else if (
                score.total_score ==
                true_score.total_score
            ) {
                ++total_equal;
            }

            const u64 distance =
                candidate > static_cast<int>(c.q)
                    ? static_cast<u64>(
                        candidate
                    ) - c.q
                    : c.q -
                      static_cast<u64>(
                          candidate
                      );

            sum_control_distance +=
                static_cast<long double>(
                    distance
                );

            max_control_distance =
                std::max(
                    max_control_distance,
                    distance
                );
        }

        /*
         * Rank 1 means no control scored higher.
         */
        const u64 j_rank =
            1 +
            static_cast<u64>(j_above);

        const u64 total_rank =
            1 +
            static_cast<u64>(total_above);

        j_rank_sum += j_rank;
        total_rank_sum += total_rank;

        j_tie_sum +=
            static_cast<u64>(j_equal);

        total_tie_sum +=
            static_cast<u64>(total_equal);

        if (j_above == 0) {
            if (j_equal == 0) {
                ++j_strict_best;
            } else {
                ++j_tied_best;
            }
        }

        if (total_above == 0) {
            if (total_equal == 0) {
                ++total_strict_best;
            } else {
                ++total_tied_best;
            }
        }

        /*
         * Percentile convention:
         *
         * Controls beaten = controls with
         * strictly LOWER score.
         *
         * This ignores ties.
         */
        const int j_lower =
            static_cast<int>(
                controls.size()
            ) -
            j_above -
            j_equal;

        const int total_lower =
            static_cast<int>(
                controls.size()
            ) -
            total_above -
            total_equal;

        if (j_lower >= 75) {
            ++j_percentile75;
        }

        if (j_lower >= 90) {
            ++j_percentile90;
        }

        if (j_lower >= 95) {
            ++j_percentile95;
        }

        if (j_lower >= 99) {
            ++j_percentile99;
        }

        if (total_lower >= 75) {
            ++total_percentile75;
        }

        if (total_lower >= 90) {
            ++total_percentile90;
        }

        if (total_lower >= 95) {
            ++total_percentile95;
        }

        if (total_lower >= 99) {
            ++total_percentile99;
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

            std::cout
                << "FIRST_TRUE_J_SCORE="
                << true_score.j_score
                << '\n';

            std::cout
                << "FIRST_TRUE_TOTAL_SCORE="
                << true_score.total_score
                << '\n';

            std::cout
                << "FIRST_CONTROL_MIN_PRIME="
                << controls.front()
                << '\n';

            std::cout
                << "FIRST_CONTROL_MAX_PRIME="
                << controls.back()
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

    const double avg_j_rank =
        true_found == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    j_rank_sum
                ) /
                static_cast<long double>(
                    true_found
                )
            );

    const double avg_total_rank =
        true_found == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    total_rank_sum
                ) /
                static_cast<long double>(
                    true_found
                )
            );

    const double avg_j_ties =
        true_found == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    j_tie_sum
                ) /
                static_cast<long double>(
                    true_found
                )
            );

    const double avg_total_ties =
        true_found == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    total_tie_sum
                ) /
                static_cast<long double>(
                    true_found
                )
            );

    const double avg_control_distance =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_control_distance /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double j_percentile =
        true_found == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      j_rank_sum
                  ) /
                  static_cast<double>(
                      true_found *
                      CONTROL_COUNT
                  )
              );

    const double total_percentile =
        true_found == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      total_rank_sum
                  ) /
                  static_cast<double>(
                      true_found *
                      CONTROL_COUNT
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
        << "CONTROL_PAIRS="
        << control_pairs
        << '\n';

    std::cout
        << "TRUE_J0="
        << true_j0
        << '\n';

    std::cout
        << "TRUE_J1="
        << true_j1
        << '\n';

    std::cout
        << "TRUE_J2="
        << true_j2
        << '\n';

    std::cout
        << "CONTROL_J0="
        << control_j0
        << '\n';

    std::cout
        << "CONTROL_J1="
        << control_j1
        << '\n';

    std::cout
        << "CONTROL_J2="
        << control_j2
        << '\n';

    std::cout
        << "TRUE_TOTAL0="
        << true_total0
        << '\n';

    std::cout
        << "TRUE_TOTAL1="
        << true_total1
        << '\n';

    std::cout
        << "TRUE_TOTAL2="
        << true_total2
        << '\n';

    std::cout
        << "TRUE_TOTAL3="
        << true_total3
        << '\n';

    std::cout
        << "TRUE_TOTAL4="
        << true_total4
        << '\n';

    std::cout
        << "CONTROL_TOTAL0="
        << control_total0
        << '\n';

    std::cout
        << "CONTROL_TOTAL1="
        << control_total1
        << '\n';

    std::cout
        << "CONTROL_TOTAL2="
        << control_total2
        << '\n';

    std::cout
        << "CONTROL_TOTAL3="
        << control_total3
        << '\n';

    std::cout
        << "CONTROL_TOTAL4="
        << control_total4
        << '\n';

    std::cout
        << "J_STRICT_BEST="
        << j_strict_best
        << '\n';

    std::cout
        << "J_TIED_BEST="
        << j_tied_best
        << '\n';

    std::cout
        << "TOTAL_STRICT_BEST="
        << total_strict_best
        << '\n';

    std::cout
        << "TOTAL_TIED_BEST="
        << total_tied_best
        << '\n';

    std::cout
        << "AVG_J_RANK="
        << avg_j_rank
        << '\n';

    std::cout
        << "AVG_TOTAL_RANK="
        << avg_total_rank
        << '\n';

    std::cout
        << "AVG_J_TIES="
        << avg_j_ties
        << '\n';

    std::cout
        << "AVG_TOTAL_TIES="
        << avg_total_ties
        << '\n';

    std::cout
        << "AVG_J_PERCENTILE="
        << j_percentile
        << '\n';

    std::cout
        << "AVG_TOTAL_PERCENTILE="
        << total_percentile
        << '\n';

    std::cout
        << "J_PERCENTILE_75="
        << j_percentile75
        << '\n';

    std::cout
        << "J_PERCENTILE_90="
        << j_percentile90
        << '\n';

    std::cout
        << "J_PERCENTILE_95="
        << j_percentile95
        << '\n';

    std::cout
        << "J_PERCENTILE_99="
        << j_percentile99
        << '\n';

    std::cout
        << "TOTAL_PERCENTILE_75="
        << total_percentile75
        << '\n';

    std::cout
        << "TOTAL_PERCENTILE_90="
        << total_percentile90
        << '\n';

    std::cout
        << "TOTAL_PERCENTILE_95="
        << total_percentile95
        << '\n';

    std::cout
        << "TOTAL_PERCENTILE_99="
        << total_percentile99
        << '\n';

    std::cout
        << "AVG_CONTROL_DISTANCE="
        << avg_control_distance
        << '\n';

    std::cout
        << "MAX_CONTROL_DISTANCE="
        << max_control_distance
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
