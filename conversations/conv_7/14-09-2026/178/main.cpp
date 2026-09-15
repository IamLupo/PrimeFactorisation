#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 471;

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

struct ProductRelation {
    u64 gcd_m = 0;
    u64 residual = 0;
    u64 gcd_match = 0;

    bool exact = false;

    long double strength = 0.0L;

    u64 quotient = 0;
    u64 remainder = 0;
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
             * Determinant-one condition:
             *
             *     m2*k1 - m1*k2 = 1
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

ProductRelation analyze_product_relation(
    u64 n,
    const Representation& p,
    const Representation& other
) {
    ProductRelation result;

    result.gcd_m =
        std::gcd(
            p.m,
            other.m
        );

    const u64 base_product =
        p.a * other.a;

    if (n >= base_product) {
        result.residual =
            n - base_product;
    } else {
        result.residual =
            base_product - n;
    }

    result.gcd_match =
        std::gcd(
            result.gcd_m,
            result.residual
        );

    result.exact =
        result.gcd_m != 0 &&
        result.residual %
            result.gcd_m == 0;

    if (result.gcd_m != 0) {
        result.strength =
            static_cast<long double>(
                result.gcd_match
            ) /
            static_cast<long double>(
                result.gcd_m
            );

        result.quotient =
            result.residual /
            result.gcd_m;

        result.remainder =
            result.residual %
            result.gcd_m;
    }

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

void print_relation(
    const char* label,
    const ProductRelation& relation
) {
    std::cout
        << label
        << "_GCD_M="
        << relation.gcd_m
        << '\n';

    std::cout
        << label
        << "_RESIDUAL="
        << relation.residual
        << '\n';

    std::cout
        << label
        << "_GCD_MATCH="
        << relation.gcd_match
        << '\n';

    std::cout
        << label
        << "_EXACT="
        << (relation.exact ? 1 : 0)
        << '\n';

    std::cout
        << label
        << "_STRENGTH="
        << relation.strength
        << '\n';

    std::cout
        << label
        << "_QUOTIENT="
        << relation.quotient
        << '\n';

    std::cout
        << label
        << "_REMAINDER="
        << relation.remainder
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47120260915ULL
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
        << "CONTROL_COUNT="
        << CONTROL_COUNT
        << '\n';

    u64 true_cases = 0;
    u64 control_pairs = 0;

    /*
     * Exact necessary congruence:
     *
     * gcd(Mp,Mq) | N-ap*aq
     */
    u64 true_exact = 0;
    u64 control_exact = 0;

    /*
     * Strength:
     *
     * gcd(g, residual) / g
     *
     * where g = gcd(Mp,Mr).
     */
    long double true_strength_sum = 0.0L;
    long double control_strength_sum = 0.0L;

    u64 true_strength_ge_25 = 0;
    u64 true_strength_ge_50 = 0;
    u64 true_strength_eq_1 = 0;

    u64 control_strength_ge_25 = 0;
    u64 control_strength_ge_50 = 0;
    u64 control_strength_eq_1 = 0;

    /*
     * GCD sizes.
     */
    long double true_gcd_m_sum = 0.0L;
    long double control_gcd_m_sum = 0.0L;

    long double true_residual_sum = 0.0L;
    long double control_residual_sum = 0.0L;

    /*
     * Residual normalized by gcd(Mp,Mr).
     */
    long double true_quotient_sum = 0.0L;
    long double control_quotient_sum = 0.0L;

    /*
     * Remainders.
     */
    long double true_remainder_sum = 0.0L;
    long double control_remainder_sum = 0.0L;

    /*
     * Rank:
     *
     * larger strength = better
     * smaller remainder = better
     * larger gcd match = better
     */
    u64 strength_rank_sum = 0;
    u64 remainder_rank_sum = 0;
    u64 gcd_match_rank_sum = 0;

    u64 strength_strict_best = 0;
    u64 strength_tied_best = 0;

    u64 remainder_strict_best = 0;
    u64 remainder_tied_best = 0;

    u64 gcd_match_strict_best = 0;
    u64 gcd_match_tied_best = 0;

    u64 strength_tie_sum = 0;
    u64 remainder_tie_sum = 0;
    u64 gcd_match_tie_sum = 0;

    /*
     * Check the exact factor-specific identity directly:
     *
     * N = (ap+jpMp)(aq+jqMq)
     */
    u64 factor_identity_failures = 0;

    /*
     * Is the true gcd condition much stronger than
     * the random-control gcd condition?
     */
    u64 true_gcd_eq_1 = 0;
    u64 control_gcd_eq_1 = 0;

    u64 first_example = 0;

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

        /*
         * Verify exact product identity.
         */
        const u64 reconstructed =
            (
                rp.a +
                rp.j * rp.m
            ) *
            (
                rq.a +
                rq.j * rq.m
            );

        if (reconstructed != c.n) {
            ++factor_identity_failures;
        }

        const ProductRelation true_relation =
            analyze_product_relation(
                c.n,
                rp,
                rq
            );

        if (true_relation.exact) {
            ++true_exact;
        }

        if (true_relation.gcd_m == 1) {
            ++true_gcd_eq_1;
        }

        true_strength_sum +=
            true_relation.strength;

        true_gcd_m_sum +=
            static_cast<long double>(
                true_relation.gcd_m
            );

        true_residual_sum +=
            static_cast<long double>(
                true_relation.residual
            );

        true_quotient_sum +=
            static_cast<long double>(
                true_relation.quotient
            );

        true_remainder_sum +=
            static_cast<long double>(
                true_relation.remainder
            );

        if (true_relation.strength >= 0.25L) {
            ++true_strength_ge_25;
        }

        if (true_relation.strength >= 0.50L) {
            ++true_strength_ge_50;
        }

        if (true_relation.strength == 1.0L) {
            ++true_strength_eq_1;
        }

        const std::vector<int> controls =
            build_nearest_controls(
                primes,
                c.q,
                c.p
            );

        int strength_above = 0;
        int strength_equal = 0;

        int remainder_above = 0;
        int remainder_equal = 0;

        int gcd_match_above = 0;
        int gcd_match_equal = 0;

        for (int candidate : controls) {
            const Representation& rr =
                cache[candidate];

            if (!rr.found) {
                continue;
            }

            const ProductRelation relation =
                analyze_product_relation(
                    c.n,
                    rp,
                    rr
                );

            ++control_pairs;

            if (relation.exact) {
                ++control_exact;
            }

            if (relation.gcd_m == 1) {
                ++control_gcd_eq_1;
            }

            control_strength_sum +=
                relation.strength;

            control_gcd_m_sum +=
                static_cast<long double>(
                    relation.gcd_m
                );

            control_residual_sum +=
                static_cast<long double>(
                    relation.residual
                );

            control_quotient_sum +=
                static_cast<long double>(
                    relation.quotient
                );

            control_remainder_sum +=
                static_cast<long double>(
                    relation.remainder
                );

            if (relation.strength >= 0.25L) {
                ++control_strength_ge_25;
            }

            if (relation.strength >= 0.50L) {
                ++control_strength_ge_50;
            }

            if (relation.strength == 1.0L) {
                ++control_strength_eq_1;
            }

            /*
             * Strength: larger is better.
             */
            if (
                relation.strength >
                true_relation.strength
            ) {
                ++strength_above;
            } else if (
                relation.strength ==
                true_relation.strength
            ) {
                ++strength_equal;
            }

            /*
             * Remainder: smaller is better.
             */
            if (
                relation.remainder <
                true_relation.remainder
            ) {
                ++remainder_above;
            } else if (
                relation.remainder ==
                true_relation.remainder
            ) {
                ++remainder_equal;
            }

            /*
             * gcd match: larger is better.
             */
            if (
                relation.gcd_match >
                true_relation.gcd_match
            ) {
                ++gcd_match_above;
            } else if (
                relation.gcd_match ==
                true_relation.gcd_match
            ) {
                ++gcd_match_equal;
            }
        }

        strength_rank_sum +=
            1 +
            static_cast<u64>(
                strength_above
            );

        remainder_rank_sum +=
            1 +
            static_cast<u64>(
                remainder_above
            );

        gcd_match_rank_sum +=
            1 +
            static_cast<u64>(
                gcd_match_above
            );

        strength_tie_sum +=
            static_cast<u64>(
                strength_equal
            );

        remainder_tie_sum +=
            static_cast<u64>(
                remainder_equal
            );

        gcd_match_tie_sum +=
            static_cast<u64>(
                gcd_match_equal
            );

        if (strength_above == 0) {
            if (strength_equal == 0) {
                ++strength_strict_best;
            } else {
                ++strength_tied_best;
            }
        }

        if (remainder_above == 0) {
            if (remainder_equal == 0) {
                ++remainder_strict_best;
            } else {
                ++remainder_tied_best;
            }
        }

        if (gcd_match_above == 0) {
            if (gcd_match_equal == 0) {
                ++gcd_match_strict_best;
            } else {
                ++gcd_match_tied_best;
            }
        }

        if (first_example == 0) {
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

            print_relation(
                "FIRST_TRUE",
                true_relation
            );

            first_example = 1;
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

    const double avg_true_strength =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                true_strength_sum /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_control_strength =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                control_strength_sum /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_true_gcd_m =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                true_gcd_m_sum /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_control_gcd_m =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                control_gcd_m_sum /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_true_residual =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                true_residual_sum /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_control_residual =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                control_residual_sum /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_true_quotient =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                true_quotient_sum /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_control_quotient =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                control_quotient_sum /
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

    const double avg_strength_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    strength_rank_sum
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

    const double avg_gcd_match_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    gcd_match_rank_sum
                ) /
                static_cast<long double>(
                    true_cases
                )
            );

    const double strength_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      strength_rank_sum
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

    const double gcd_match_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      gcd_match_rank_sum
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
        << "FACTOR_IDENTITY_FAILURES="
        << factor_identity_failures
        << '\n';

    std::cout
        << "TRUE_EXACT_GCD_CONGRUENCE="
        << true_exact
        << '\n';

    std::cout
        << "CONTROL_EXACT_GCD_CONGRUENCE="
        << control_exact
        << '\n';

    std::cout
        << "TRUE_GCD_M_EQUALS_1="
        << true_gcd_eq_1
        << '\n';

    std::cout
        << "CONTROL_GCD_M_EQUALS_1="
        << control_gcd_eq_1
        << '\n';

    std::cout
        << "AVG_TRUE_GCD_M="
        << avg_true_gcd_m
        << '\n';

    std::cout
        << "AVG_CONTROL_GCD_M="
        << avg_control_gcd_m
        << '\n';

    std::cout
        << "AVG_TRUE_STRENGTH="
        << avg_true_strength
        << '\n';

    std::cout
        << "AVG_CONTROL_STRENGTH="
        << avg_control_strength
        << '\n';

    std::cout
        << "TRUE_STRENGTH_GE_25="
        << true_strength_ge_25
        << '\n';

    std::cout
        << "TRUE_STRENGTH_GE_50="
        << true_strength_ge_50
        << '\n';

    std::cout
        << "TRUE_STRENGTH_EQ_1="
        << true_strength_eq_1
        << '\n';

    std::cout
        << "CONTROL_STRENGTH_GE_25="
        << control_strength_ge_25
        << '\n';

    std::cout
        << "CONTROL_STRENGTH_GE_50="
        << control_strength_ge_50
        << '\n';

    std::cout
        << "CONTROL_STRENGTH_EQ_1="
        << control_strength_eq_1
        << '\n';

    std::cout
        << "AVG_TRUE_RESIDUAL="
        << avg_true_residual
        << '\n';

    std::cout
        << "AVG_CONTROL_RESIDUAL="
        << avg_control_residual
        << '\n';

    std::cout
        << "AVG_TRUE_QUOTIENT="
        << avg_true_quotient
        << '\n';

    std::cout
        << "AVG_CONTROL_QUOTIENT="
        << avg_control_quotient
        << '\n';

    std::cout
        << "AVG_TRUE_REMAINDER="
        << avg_true_remainder
        << '\n';

    std::cout
        << "AVG_CONTROL_REMAINDER="
        << avg_control_remainder
        << '\n';

    std::cout
        << "AVG_STRENGTH_RANK="
        << avg_strength_rank
        << '\n';

    std::cout
        << "AVG_REMAINDER_RANK="
        << avg_remainder_rank
        << '\n';

    std::cout
        << "AVG_GCD_MATCH_RANK="
        << avg_gcd_match_rank
        << '\n';

    std::cout
        << "STRENGTH_PERCENTILE="
        << strength_percentile
        << '\n';

    std::cout
        << "REMAINDER_PERCENTILE="
        << remainder_percentile
        << '\n';

    std::cout
        << "GCD_MATCH_PERCENTILE="
        << gcd_match_percentile
        << '\n';

    std::cout
        << "STRENGTH_STRICT_BEST="
        << strength_strict_best
        << '\n';

    std::cout
        << "STRENGTH_TIED_BEST="
        << strength_tied_best
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
        << "GCD_MATCH_STRICT_BEST="
        << gcd_match_strict_best
        << '\n';

    std::cout
        << "GCD_MATCH_TIED_BEST="
        << gcd_match_tied_best
        << '\n';

    std::cout
        << "AVG_STRENGTH_TIES="
        << (
            true_cases == 0
                ? 0.0
                : static_cast<double>(
                    strength_tie_sum /
                    static_cast<long double>(
                        true_cases
                    )
                )
        )
        << '\n';

    std::cout
        << "AVG_REMAINDER_TIES="
        << (
            true_cases == 0
                ? 0.0
                : static_cast<double>(
                    remainder_tie_sum /
                    static_cast<long double>(
                        true_cases
                    )
                )
        )
        << '\n';

    std::cout
        << "AVG_GCD_MATCH_TIES="
        << (
            true_cases == 0
                ? 0.0
                : static_cast<double>(
                    gcd_match_tie_sum /
                    static_cast<long double>(
                        true_cases
                    )
                )
        )
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
