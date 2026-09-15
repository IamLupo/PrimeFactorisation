#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 470;

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

struct RelationData {
    u64 j_fixed = 0;
    u64 j_other = 0;
    u64 m_other = 0;

    u64 gcd_j_m = 0;
    u64 reduced_j = 0;

    bool exact_division = false;
    bool reduced_division = false;

    u64 quotient = 0;
    u64 reduced_quotient = 0;

    u64 remainder = 0;
    u64 reduced_remainder = 0;
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

RelationData analyze_relation(
    const Representation& fixed,
    const Representation& other
) {
    RelationData result;

    result.j_fixed = fixed.j;
    result.j_other = other.j;
    result.m_other = other.m;

    result.gcd_j_m =
        std::gcd(
            fixed.j,
            other.m
        );

    result.reduced_j =
        fixed.j /
        result.gcd_j_m;

    result.exact_division =
        (other.j * other.m) %
        fixed.j == 0;

    result.reduced_division =
        other.j %
        result.reduced_j == 0;

    if (result.exact_division) {
        result.quotient =
            (other.j * other.m) /
            fixed.j;
    }

    if (result.reduced_division) {
        result.reduced_quotient =
            other.j /
            result.reduced_j;
    }

    result.remainder =
        (other.j * other.m) %
        fixed.j;

    result.reduced_remainder =
        other.j %
        result.reduced_j;

    return result;
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
        bool take_left;

        if (left < 0) {
            take_left = false;
        } else if (
            right >=
            static_cast<int>(primes.size())
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

void print_relation(
    const char* label,
    const RelationData& relation
) {
    std::cout
        << label
        << "_J_FIXED="
        << relation.j_fixed
        << '\n';

    std::cout
        << label
        << "_J_OTHER="
        << relation.j_other
        << '\n';

    std::cout
        << label
        << "_M_OTHER="
        << relation.m_other
        << '\n';

    std::cout
        << label
        << "_GCD_J_M="
        << relation.gcd_j_m
        << '\n';

    std::cout
        << label
        << "_REDUCED_J="
        << relation.reduced_j
        << '\n';

    std::cout
        << label
        << "_EXACT_DIVISION="
        << (relation.exact_division ? 1 : 0)
        << '\n';

    std::cout
        << label
        << "_REDUCED_DIVISION="
        << (relation.reduced_division ? 1 : 0)
        << '\n';

    std::cout
        << label
        << "_QUOTIENT="
        << relation.quotient
        << '\n';

    std::cout
        << label
        << "_REDUCED_QUOTIENT="
        << relation.reduced_quotient
        << '\n';

    std::cout
        << label
        << "_REMAINDER="
        << relation.remainder
        << '\n';

    std::cout
        << label
        << "_REDUCED_REMAINDER="
        << relation.reduced_remainder
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47020260915ULL
    );

    const std::vector<int> primes =
        generate_primes(PRIME_LIMIT);

    const std::vector<PrimeCase> cases =
        generate_cases(primes, rng);

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
        << "CONTROL_COUNT="
        << CONTROL_COUNT
        << '\n';

    u64 true_cases = 0;
    u64 control_pairs = 0;

    /*
     * Exact-relation counts.
     */
    u64 true_exact = 0;
    u64 control_exact = 0;

    u64 true_reduced = 0;
    u64 control_reduced = 0;

    /*
     * Distribution of reduced_j.
     *
     * We group it into useful size ranges.
     */
    u64 true_reduced_j_1 = 0;
    u64 true_reduced_j_2_10 = 0;
    u64 true_reduced_j_11_100 = 0;
    u64 true_reduced_j_gt100 = 0;

    u64 control_reduced_j_1 = 0;
    u64 control_reduced_j_2_10 = 0;
    u64 control_reduced_j_11_100 = 0;
    u64 control_reduced_j_gt100 = 0;

    /*
     * gcd(j_fixed,M_other) / j_fixed.
     *
     * This is another direct measure of how much
     * of j_fixed is already supplied by M_other.
     */
    long double true_gcd_fraction_sum = 0.0L;
    long double control_gcd_fraction_sum = 0.0L;

    /*
     * Quotient when exact divisibility holds.
     */
    long double true_quotient_sum = 0.0L;
    long double control_quotient_sum = 0.0L;

    u64 true_quotient_count = 0;
    u64 control_quotient_count = 0;

    /*
     * Reduced quotient when reduced divisibility holds.
     */
    long double true_reduced_quotient_sum = 0.0L;
    long double control_reduced_quotient_sum = 0.0L;

    u64 true_reduced_quotient_count = 0;
    u64 control_reduced_quotient_count = 0;

    /*
     * Remainders.
     */
    long double true_remainder_sum = 0.0L;
    long double control_remainder_sum = 0.0L;

    /*
     * Compare true q against controls using:
     *
     *   reduced_j
     *   gcd fraction
     *   normalized remainder
     */
    u64 reduced_j_rank_sum = 0;
    u64 gcd_fraction_rank_sum = 0;
    u64 remainder_rank_sum = 0;

    u64 reduced_j_strict_best = 0;
    u64 gcd_fraction_strict_best = 0;
    u64 remainder_strict_best = 0;

    u64 reduced_j_tied_best = 0;
    u64 gcd_fraction_tied_best = 0;
    u64 remainder_tied_best = 0;

    u64 reduced_j_tie_sum = 0;
    u64 gcd_fraction_tie_sum = 0;
    u64 remainder_tie_sum = 0;

    /*
     * Threshold tests.
     */
    u64 true_gcd_fraction_ge_25 = 0;
    u64 true_gcd_fraction_ge_50 = 0;
    u64 true_gcd_fraction_eq_1 = 0;

    u64 control_gcd_fraction_ge_25 = 0;
    u64 control_gcd_fraction_ge_50 = 0;
    u64 control_gcd_fraction_eq_1 = 0;

    /*
     * Stratify by j_fixed.
     */
    u64 jp_small_cases = 0;
    u64 jp_medium_cases = 0;
    u64 jp_large_cases = 0;

    u64 jp_small_true_exact = 0;
    u64 jp_medium_true_exact = 0;
    u64 jp_large_true_exact = 0;

    u64 jp_small_control_exact = 0;
    u64 jp_medium_control_exact = 0;
    u64 jp_large_control_exact = 0;

    u64 first_example_printed = 0;

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

        ++true_cases;

        const RelationData true_relation =
            analyze_relation(
                rp,
                rq
            );

        if (true_relation.exact_division) {
            ++true_exact;
        }

        if (true_relation.reduced_division) {
            ++true_reduced;
        }

        if (true_relation.reduced_j == 1) {
            ++true_reduced_j_1;
        } else if (
            true_relation.reduced_j <= 10
        ) {
            ++true_reduced_j_2_10;
        } else if (
            true_relation.reduced_j <= 100
        ) {
            ++true_reduced_j_11_100;
        } else {
            ++true_reduced_j_gt100;
        }

        const long double true_gcd_fraction =
            static_cast<long double>(
                true_relation.gcd_j_m
            ) /
            static_cast<long double>(
                true_relation.j_fixed
            );

        true_gcd_fraction_sum +=
            true_gcd_fraction;

        if (true_gcd_fraction >= 0.25L) {
            ++true_gcd_fraction_ge_25;
        }

        if (true_gcd_fraction >= 0.50L) {
            ++true_gcd_fraction_ge_50;
        }

        if (true_gcd_fraction == 1.0L) {
            ++true_gcd_fraction_eq_1;
        }

        if (true_relation.exact_division) {
            true_quotient_sum +=
                static_cast<long double>(
                    true_relation.quotient
                );

            ++true_quotient_count;
        }

        if (true_relation.reduced_division) {
            true_reduced_quotient_sum +=
                static_cast<long double>(
                    true_relation.reduced_quotient
                );

            ++true_reduced_quotient_count;
        }

        true_remainder_sum +=
            static_cast<long double>(
                true_relation.reduced_remainder
            );

        /*
         * Stratification:
         *
         * small  <= 10
         * medium 11..100
         * large > 100
         */
        if (rp.j <= 10) {
            ++jp_small_cases;

            if (true_relation.exact_division) {
                ++jp_small_true_exact;
            }
        } else if (rp.j <= 100) {
            ++jp_medium_cases;

            if (true_relation.exact_division) {
                ++jp_medium_true_exact;
            }
        } else {
            ++jp_large_cases;

            if (true_relation.exact_division) {
                ++jp_large_true_exact;
            }
        }

        const std::vector<int> controls =
            build_nearest_controls(
                primes,
                c.q,
                c.p
            );

        int reduced_j_above = 0;
        int reduced_j_equal = 0;

        int gcd_fraction_above = 0;
        int gcd_fraction_equal = 0;

        int remainder_above = 0;
        int remainder_equal = 0;

        for (int candidate : controls) {
            const Representation& rr =
                cache[candidate];

            if (!rr.found) {
                continue;
            }

            const RelationData relation =
                analyze_relation(
                    rp,
                    rr
                );

            ++control_pairs;

            if (relation.exact_division) {
                ++control_exact;
            }

            if (relation.reduced_division) {
                ++control_reduced;
            }

            if (relation.reduced_j == 1) {
                ++control_reduced_j_1;
            } else if (
                relation.reduced_j <= 10
            ) {
                ++control_reduced_j_2_10;
            } else if (
                relation.reduced_j <= 100
            ) {
                ++control_reduced_j_11_100;
            } else {
                ++control_reduced_j_gt100;
            }

            const long double gcd_fraction =
                static_cast<long double>(
                    relation.gcd_j_m
                ) /
                static_cast<long double>(
                    relation.j_fixed
                );

            control_gcd_fraction_sum +=
                gcd_fraction;

            if (gcd_fraction >= 0.25L) {
                ++control_gcd_fraction_ge_25;
            }

            if (gcd_fraction >= 0.50L) {
                ++control_gcd_fraction_ge_50;
            }

            if (gcd_fraction == 1.0L) {
                ++control_gcd_fraction_eq_1;
            }

            if (relation.exact_division) {
                control_quotient_sum +=
                    static_cast<long double>(
                        relation.quotient
                    );

                ++control_quotient_count;
            }

            if (relation.reduced_division) {
                control_reduced_quotient_sum +=
                    static_cast<long double>(
                        relation.reduced_quotient
                    );

                ++control_reduced_quotient_count;
            }

            control_remainder_sum +=
                static_cast<long double>(
                    relation.reduced_remainder
                );

            if (rp.j <= 10 &&
                relation.exact_division) {
                ++jp_small_control_exact;
            }

            if (rp.j > 10 &&
                rp.j <= 100 &&
                relation.exact_division) {
                ++jp_medium_control_exact;
            }

            if (rp.j > 100 &&
                relation.exact_division) {
                ++jp_large_control_exact;
            }

            /*
             * Smaller reduced_j is better.
             */
            if (
                relation.reduced_j <
                true_relation.reduced_j
            ) {
                ++reduced_j_above;
            } else if (
                relation.reduced_j ==
                true_relation.reduced_j
            ) {
                ++reduced_j_equal;
            }

            /*
             * Larger gcd fraction is better.
             */
            if (
                gcd_fraction >
                true_gcd_fraction
            ) {
                ++gcd_fraction_above;
            } else if (
                gcd_fraction ==
                true_gcd_fraction
            ) {
                ++gcd_fraction_equal;
            }

            /*
             * Smaller reduced remainder is better.
             */
            if (
                relation.reduced_remainder <
                true_relation.reduced_remainder
            ) {
                ++remainder_above;
            } else if (
                relation.reduced_remainder ==
                true_relation.reduced_remainder
            ) {
                ++remainder_equal;
            }
        }

        /*
         * Rank 1 is best.
         */
        reduced_j_rank_sum +=
            1 +
            static_cast<u64>(
                reduced_j_above
            );

        gcd_fraction_rank_sum +=
            1 +
            static_cast<u64>(
                gcd_fraction_above
            );

        remainder_rank_sum +=
            1 +
            static_cast<u64>(
                remainder_above
            );

        reduced_j_tie_sum +=
            static_cast<u64>(
                reduced_j_equal
            );

        gcd_fraction_tie_sum +=
            static_cast<u64>(
                gcd_fraction_equal
            );

        remainder_tie_sum +=
            static_cast<u64>(
                remainder_equal
            );

        if (reduced_j_above == 0) {
            if (reduced_j_equal == 0) {
                ++reduced_j_strict_best;
            } else {
                ++reduced_j_tied_best;
            }
        }

        if (gcd_fraction_above == 0) {
            if (gcd_fraction_equal == 0) {
                ++gcd_fraction_strict_best;
            } else {
                ++gcd_fraction_tied_best;
            }
        }

        if (remainder_above == 0) {
            if (remainder_equal == 0) {
                ++remainder_strict_best;
            } else {
                ++remainder_tied_best;
            }
        }

        if (!first_example_printed) {
            std::cout
                << "FIRST_N="
                << c.n
                << '\n';

            std::cout
                << "FIRST_P_PRIME="
                << rp.prime
                << '\n';

            std::cout
                << "FIRST_P_A="
                << rp.a
                << '\n';

            std::cout
                << "FIRST_P_M="
                << rp.m
                << '\n';

            std::cout
                << "FIRST_P_J="
                << rp.j
                << '\n';

            std::cout
                << "FIRST_Q_PRIME="
                << rq.prime
                << '\n';

            std::cout
                << "FIRST_Q_A="
                << rq.a
                << '\n';

            std::cout
                << "FIRST_Q_M="
                << rq.m
                << '\n';

            std::cout
                << "FIRST_Q_J="
                << rq.j
                << '\n';

            print_relation(
                "FIRST_TRUE",
                true_relation
            );

            first_example_printed = 1;
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

    const double avg_true_gcd_fraction =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                true_gcd_fraction_sum /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_control_gcd_fraction =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                control_gcd_fraction_sum /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_true_remainder =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                true_remainder_sum /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_control_remainder =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                control_remainder_sum /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_true_quotient =
        true_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                true_quotient_sum /
                static_cast<long double>(
                    true_quotient_count
                )
            );

    const double avg_control_quotient =
        control_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                control_quotient_sum /
                static_cast<long double>(
                    control_quotient_count
                )
            );

    const double avg_true_reduced_quotient =
        true_reduced_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                true_reduced_quotient_sum /
                static_cast<long double>(
                    true_reduced_quotient_count
                )
            );

    const double avg_control_reduced_quotient =
        control_reduced_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                control_reduced_quotient_sum /
                static_cast<long double>(
                    control_reduced_quotient_count
                )
            );

    const double avg_reduced_j_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    reduced_j_rank_sum
                ) /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_gcd_fraction_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    gcd_fraction_rank_sum
                ) /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_remainder_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    remainder_rank_sum
                ) /
                static_cast<long double>(
                    true_cases
                )
            );

    const double reduced_j_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      reduced_j_rank_sum
                  ) /
                  static_cast<double>(
                      true_cases *
                      CONTROL_COUNT
                  )
              );

    const double gcd_fraction_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      gcd_fraction_rank_sum
                  ) /
                  static_cast<double>(
                      true_cases *
                      CONTROL_COUNT
                  )
              );

    const double remainder_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      remainder_rank_sum
                  ) /
                  static_cast<double>(
                      true_cases *
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
        << "TRUE_CASES="
        << true_cases
        << '\n';

    std::cout
        << "CONTROL_PAIRS="
        << control_pairs
        << '\n';

    std::cout
        << "TRUE_EXACT_DIVISION="
        << true_exact
        << '\n';

    std::cout
        << "CONTROL_EXACT_DIVISION="
        << control_exact
        << '\n';

    std::cout
        << "TRUE_REDUCED_DIVISION="
        << true_reduced
        << '\n';

    std::cout
        << "CONTROL_REDUCED_DIVISION="
        << control_reduced
        << '\n';

    std::cout
        << "TRUE_REDUCED_J_1="
        << true_reduced_j_1
        << '\n';

    std::cout
        << "TRUE_REDUCED_J_2_10="
        << true_reduced_j_2_10
        << '\n';

    std::cout
        << "TRUE_REDUCED_J_11_100="
        << true_reduced_j_11_100
        << '\n';

    std::cout
        << "TRUE_REDUCED_J_GT100="
        << true_reduced_j_gt100
        << '\n';

    std::cout
        << "CONTROL_REDUCED_J_1="
        << control_reduced_j_1
        << '\n';

    std::cout
        << "CONTROL_REDUCED_J_2_10="
        << control_reduced_j_2_10
        << '\n';

    std::cout
        << "CONTROL_REDUCED_J_11_100="
        << control_reduced_j_11_100
        << '\n';

    std::cout
        << "CONTROL_REDUCED_J_GT100="
        << control_reduced_j_gt100
        << '\n';

    std::cout
        << "AVG_TRUE_GCD_FRACTION="
        << avg_true_gcd_fraction
        << '\n';

    std::cout
        << "AVG_CONTROL_GCD_FRACTION="
        << avg_control_gcd_fraction
        << '\n';

    std::cout
        << "TRUE_GCD_FRACTION_GE_25="
        << true_gcd_fraction_ge_25
        << '\n';

    std::cout
        << "TRUE_GCD_FRACTION_GE_50="
        << true_gcd_fraction_ge_50
        << '\n';

    std::cout
        << "TRUE_GCD_FRACTION_EQ_1="
        << true_gcd_fraction_eq_1
        << '\n';

    std::cout
        << "CONTROL_GCD_FRACTION_GE_25="
        << control_gcd_fraction_ge_25
        << '\n';

    std::cout
        << "CONTROL_GCD_FRACTION_GE_50="
        << control_gcd_fraction_ge_50
        << '\n';

    std::cout
        << "CONTROL_GCD_FRACTION_EQ_1="
        << control_gcd_fraction_eq_1
        << '\n';

    std::cout
        << "AVG_TRUE_REDUCED_REMAINDER="
        << avg_true_remainder
        << '\n';

    std::cout
        << "AVG_CONTROL_REDUCED_REMAINDER="
        << avg_control_remainder
        << '\n';

    std::cout
        << "TRUE_EXACT_QUOTIENT_COUNT="
        << true_quotient_count
        << '\n';

    std::cout
        << "CONTROL_EXACT_QUOTIENT_COUNT="
        << control_quotient_count
        << '\n';

    std::cout
        << "AVG_TRUE_EXACT_QUOTIENT="
        << avg_true_quotient
        << '\n';

    std::cout
        << "AVG_CONTROL_EXACT_QUOTIENT="
        << avg_control_quotient
        << '\n';

    std::cout
        << "TRUE_REDUCED_QUOTIENT_COUNT="
        << true_reduced_quotient_count
        << '\n';

    std::cout
        << "CONTROL_REDUCED_QUOTIENT_COUNT="
        << control_reduced_quotient_count
        << '\n';

    std::cout
        << "AVG_TRUE_REDUCED_QUOTIENT="
        << avg_true_reduced_quotient
        << '\n';

    std::cout
        << "AVG_CONTROL_REDUCED_QUOTIENT="
        << avg_control_reduced_quotient
        << '\n';

    std::cout
        << "AVG_REDUCED_J_RANK="
        << avg_reduced_j_rank
        << '\n';

    std::cout
        << "AVG_GCD_FRACTION_RANK="
        << avg_gcd_fraction_rank
        << '\n';

    std::cout
        << "AVG_REMAINDER_RANK="
        << avg_remainder_rank
        << '\n';

    std::cout
        << "REDUCED_J_PERCENTILE="
        << reduced_j_percentile
        << '\n';

    std::cout
        << "GCD_FRACTION_PERCENTILE="
        << gcd_fraction_percentile
        << '\n';

    std::cout
        << "REMAINDER_PERCENTILE="
        << remainder_percentile
        << '\n';

    std::cout
        << "REDUCED_J_STRICT_BEST="
        << reduced_j_strict_best
        << '\n';

    std::cout
        << "REDUCED_J_TIED_BEST="
        << reduced_j_tied_best
        << '\n';

    std::cout
        << "GCD_FRACTION_STRICT_BEST="
        << gcd_fraction_strict_best
        << '\n';

    std::cout
        << "GCD_FRACTION_TIED_BEST="
        << gcd_fraction_tied_best
        << '\n';

    std::cout
        << "REMAINDER_STRICT_BEST="
        << remainder_strict_best
        << '\n';

    std::cout
        << "REMAINDER_TIED_BEST="
        << remainder_tied_best
        << '\n';

    std::cout
        << "JP_SMALL_CASES="
        << jp_small_cases
        << '\n';

    std::cout
        << "JP_MEDIUM_CASES="
        << jp_medium_cases
        << '\n';

    std::cout
        << "JP_LARGE_CASES="
        << jp_large_cases
        << '\n';

    std::cout
        << "JP_SMALL_TRUE_EXACT="
        << jp_small_true_exact
        << '\n';

    std::cout
        << "JP_MEDIUM_TRUE_EXACT="
        << jp_medium_true_exact
        << '\n';

    std::cout
        << "JP_LARGE_TRUE_EXACT="
        << jp_large_true_exact
        << '\n';

    std::cout
        << "JP_SMALL_CONTROL_EXACT="
        << jp_small_control_exact
        << '\n';

    std::cout
        << "JP_MEDIUM_CONTROL_EXACT="
        << jp_medium_control_exact
        << '\n';

    std::cout
        << "JP_LARGE_CONTROL_EXACT="
        << jp_large_control_exact
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
