#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 469;

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

struct Strength {
    u64 gcd_forward = 0;
    u64 gcd_reverse = 0;

    u64 quotient_forward = 0;
    u64 quotient_reverse = 0;

    long double rho_forward = 0.0L;
    long double rho_reverse = 0.0L;

    u64 lcm_residual_forward = 0;
    u64 lcm_residual_reverse = 0;

    u64 remainder_forward = 0;
    u64 remainder_reverse = 0;
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
            static_cast<u64>(candidates[dist(rng)]);

        u64 q =
            static_cast<u64>(candidates[dist(rng)]);

        while (q == p) {
            q =
                static_cast<u64>(candidates[dist(rng)]);
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
             * for 1 <= m1,m2 <= M_LIMIT.
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
        } else if (j == best.j && m < best.m) {
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

Strength calculate_strength(
    const Representation& fixed,
    const Representation& other
) {
    Strength result;

    if (!fixed.found || !other.found) {
        return result;
    }

    /*
     * Forward:
     *
     *   jp | jr*Mr
     */
    const u64 forward_value =
        other.j * other.m;

    result.gcd_forward =
        std::gcd(
            fixed.j,
            forward_value
        );

    result.rho_forward =
        static_cast<long double>(
            result.gcd_forward
        ) /
        static_cast<long double>(
            fixed.j
        );

    if (result.gcd_forward != 0 &&
        forward_value %
            fixed.j == 0) {
        result.quotient_forward =
            forward_value /
            fixed.j;
    }

    result.remainder_forward =
        forward_value %
        fixed.j;

    /*
     * Reverse:
     *
     *   jr | jp*Mp
     */
    const u64 reverse_value =
        fixed.j * fixed.m;

    result.gcd_reverse =
        std::gcd(
            other.j,
            reverse_value
        );

    result.rho_reverse =
        static_cast<long double>(
            result.gcd_reverse
        ) /
        static_cast<long double>(
            other.j
        );

    if (result.gcd_reverse != 0 &&
        reverse_value %
            other.j == 0) {
        result.quotient_reverse =
            reverse_value /
            other.j;
    }

    result.remainder_reverse =
        reverse_value %
        other.j;

    /*
     * How much of the target j is covered
     * after combining the other j and M.
     *
     * This gives a useful alternative to the
     * normalized gcd alone.
     */
    const u64 gcd_fm =
        std::gcd(
            fixed.j,
            other.m
        );

    result.lcm_residual_forward =
        fixed.j / gcd_fm;

    const u64 gcd_rm =
        std::gcd(
            other.j,
            fixed.m
        );

    result.lcm_residual_reverse =
        other.j / gcd_rm;

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
         right <
             static_cast<int>(
                 primes.size()
             ))
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
            const u64 left_prime =
                static_cast<u64>(
                    primes[left]
                );

            const u64 right_prime =
                static_cast<u64>(
                    primes[right]
                );

            const u64 left_distance =
                q > left_prime
                    ? q - left_prime
                    : left_prime - q;

            const u64 right_distance =
                q > right_prime
                    ? q - right_prime
                    : right_prime - q;

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

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x46920260915ULL
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

    /*
     * True sums.
     */
    long double sum_rho_f_true = 0.0L;
    long double sum_rho_r_true = 0.0L;

    long double sum_gcd_f_ratio_true = 0.0L;
    long double sum_gcd_r_ratio_true = 0.0L;

    long double sum_residual_f_true = 0.0L;
    long double sum_residual_r_true = 0.0L;

    long double sum_lcm_residual_f_true = 0.0L;
    long double sum_lcm_residual_r_true = 0.0L;

    /*
     * Control sums.
     */
    long double sum_rho_f_control = 0.0L;
    long double sum_rho_r_control = 0.0L;

    long double sum_gcd_f_ratio_control = 0.0L;
    long double sum_gcd_r_ratio_control = 0.0L;

    long double sum_residual_f_control = 0.0L;
    long double sum_residual_r_control = 0.0L;

    long double sum_lcm_residual_f_control = 0.0L;
    long double sum_lcm_residual_r_control = 0.0L;

    u64 true_cases = 0;
    u64 control_pairs = 0;

    /*
     * Exact divisibility counts.
     */
    u64 true_forward_divides = 0;
    u64 true_reverse_divides = 0;
    u64 true_both_divide = 0;

    u64 control_forward_divides = 0;
    u64 control_reverse_divides = 0;
    u64 control_both_divide = 0;

    /*
     * Strong thresholds for rho.
     */
    u64 true_forward_rho_ge_50 = 0;
    u64 true_reverse_rho_ge_50 = 0;

    u64 true_forward_rho_ge_25 = 0;
    u64 true_reverse_rho_ge_25 = 0;

    u64 control_forward_rho_ge_50 = 0;
    u64 control_reverse_rho_ge_50 = 0;

    u64 control_forward_rho_ge_25 = 0;
    u64 control_reverse_rho_ge_25 = 0;

    /*
     * Rank statistics for rho_sum and rho_product.
     */
    u64 rho_sum_rank_sum = 0;
    u64 rho_product_rank_sum = 0;

    u64 rho_sum_true_strict_best = 0;
    u64 rho_product_true_strict_best = 0;

    u64 rho_sum_true_tied_best = 0;
    u64 rho_product_true_tied_best = 0;

    u64 rho_sum_tie_sum = 0;
    u64 rho_product_tie_sum = 0;

    /*
     * Average quotient when exact divisibility occurs.
     */
    long double sum_true_forward_quotient = 0.0L;
    long double sum_true_reverse_quotient = 0.0L;

    long double sum_control_forward_quotient = 0.0L;
    long double sum_control_reverse_quotient = 0.0L;

    u64 true_forward_quotient_count = 0;
    u64 true_reverse_quotient_count = 0;

    u64 control_forward_quotient_count = 0;
    u64 control_reverse_quotient_count = 0;

    bool first_example = false;

    for (int case_index = 0;
         case_index < CASE_COUNT;
         ++case_index) {

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

        const Strength true_strength =
            calculate_strength(
                rp,
                rq
            );

        const long double true_rho_sum =
            true_strength.rho_forward +
            true_strength.rho_reverse;

        const long double true_rho_product =
            true_strength.rho_forward *
            true_strength.rho_reverse;

        sum_rho_f_true +=
            true_strength.rho_forward;

        sum_rho_r_true +=
            true_strength.rho_reverse;

        sum_gcd_f_ratio_true +=
            true_strength.rho_forward;

        sum_gcd_r_ratio_true +=
            true_strength.rho_reverse;

        sum_residual_f_true +=
            static_cast<long double>(
                true_strength.remainder_forward
            );

        sum_residual_r_true +=
            static_cast<long double>(
                true_strength.remainder_reverse
            );

        sum_lcm_residual_f_true +=
            static_cast<long double>(
                true_strength
                    .lcm_residual_forward
            );

        sum_lcm_residual_r_true +=
            static_cast<long double>(
                true_strength
                    .lcm_residual_reverse
            );

        if (
            true_strength.remainder_forward ==
            0
        ) {
            ++true_forward_divides;
        }

        if (
            true_strength.remainder_reverse ==
            0
        ) {
            ++true_reverse_divides;
        }

        if (
            true_strength.remainder_forward == 0 &&
            true_strength.remainder_reverse == 0
        ) {
            ++true_both_divide;
        }

        if (
            true_strength.rho_forward >=
            0.5L
        ) {
            ++true_forward_rho_ge_50;
        }

        if (
            true_strength.rho_reverse >=
            0.5L
        ) {
            ++true_reverse_rho_ge_50;
        }

        if (
            true_strength.rho_forward >=
            0.25L
        ) {
            ++true_forward_rho_ge_25;
        }

        if (
            true_strength.rho_reverse >=
            0.25L
        ) {
            ++true_reverse_rho_ge_25;
        }

        if (
            true_strength.quotient_forward !=
            0
        ) {
            sum_true_forward_quotient +=
                static_cast<long double>(
                    true_strength
                        .quotient_forward
                );

            ++true_forward_quotient_count;
        }

        if (
            true_strength.quotient_reverse !=
            0
        ) {
            sum_true_reverse_quotient +=
                static_cast<long double>(
                    true_strength
                        .quotient_reverse
                );

            ++true_reverse_quotient_count;
        }

        const std::vector<int> controls =
            build_nearest_controls(
                primes,
                c.q,
                c.p
            );

        int sum_above = 0;
        int sum_equal = 0;

        int product_above = 0;
        int product_equal = 0;

        for (int candidate : controls) {
            const Representation& rr =
                cache[candidate];

            if (!rr.found) {
                continue;
            }

            const Strength strength =
                calculate_strength(
                    rp,
                    rr
                );

            ++control_pairs;

            sum_rho_f_control +=
                strength.rho_forward;

            sum_rho_r_control +=
                strength.rho_reverse;

            sum_gcd_f_ratio_control +=
                strength.rho_forward;

            sum_gcd_r_ratio_control +=
                strength.rho_reverse;

            sum_residual_f_control +=
                static_cast<long double>(
                    strength.remainder_forward
                );

            sum_residual_r_control +=
                static_cast<long double>(
                    strength.remainder_reverse
                );

            sum_lcm_residual_f_control +=
                static_cast<long double>(
                    strength
                        .lcm_residual_forward
                );

            sum_lcm_residual_r_control +=
                static_cast<long double>(
                    strength
                        .lcm_residual_reverse
                );

            if (
                strength.remainder_forward ==
                0
            ) {
                ++control_forward_divides;
            }

            if (
                strength.remainder_reverse ==
                0
            ) {
                ++control_reverse_divides;
            }

            if (
                strength.remainder_forward == 0 &&
                strength.remainder_reverse == 0
            ) {
                ++control_both_divide;
            }

            if (
                strength.rho_forward >=
                0.5L
            ) {
                ++control_forward_rho_ge_50;
            }

            if (
                strength.rho_reverse >=
                0.5L
            ) {
                ++control_reverse_rho_ge_50;
            }

            if (
                strength.rho_forward >=
                0.25L
            ) {
                ++control_forward_rho_ge_25;
            }

            if (
                strength.rho_reverse >=
                0.25L
            ) {
                ++control_reverse_rho_ge_25;
            }

            if (
                strength.quotient_forward !=
                0
            ) {
                sum_control_forward_quotient +=
                    static_cast<long double>(
                        strength
                            .quotient_forward
                    );

                ++control_forward_quotient_count;
            }

            if (
                strength.quotient_reverse !=
                0
            ) {
                sum_control_reverse_quotient +=
                    static_cast<long double>(
                        strength
                            .quotient_reverse
                    );

                ++control_reverse_quotient_count;
            }

            const long double rho_sum =
                strength.rho_forward +
                strength.rho_reverse;

            const long double rho_product =
                strength.rho_forward *
                strength.rho_reverse;

            if (rho_sum > true_rho_sum) {
                ++sum_above;
            } else if (
                rho_sum == true_rho_sum
            ) {
                ++sum_equal;
            }

            if (
                rho_product >
                true_rho_product
            ) {
                ++product_above;
            } else if (
                rho_product ==
                true_rho_product
            ) {
                ++product_equal;
            }
        }

        rho_sum_rank_sum +=
            1 +
            static_cast<u64>(sum_above);

        rho_product_rank_sum +=
            1 +
            static_cast<u64>(product_above);

        rho_sum_tie_sum +=
            static_cast<u64>(sum_equal);

        rho_product_tie_sum +=
            static_cast<u64>(product_equal);

        if (sum_above == 0) {
            if (sum_equal == 0) {
                ++rho_sum_true_strict_best;
            } else {
                ++rho_sum_true_tied_best;
            }
        }

        if (product_above == 0) {
            if (product_equal == 0) {
                ++rho_product_true_strict_best;
            } else {
                ++rho_product_true_tied_best;
            }
        }

        if (!first_example) {
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

            std::cout
                << "FIRST_TRUE_GCD_FORWARD="
                << true_strength.gcd_forward
                << '\n';

            std::cout
                << "FIRST_TRUE_GCD_REVERSE="
                << true_strength.gcd_reverse
                << '\n';

            std::cout
                << "FIRST_TRUE_RHO_FORWARD="
                << true_strength.rho_forward
                << '\n';

            std::cout
                << "FIRST_TRUE_RHO_REVERSE="
                << true_strength.rho_reverse
                << '\n';

            std::cout
                << "FIRST_TRUE_RHO_SUM="
                << true_rho_sum
                << '\n';

            std::cout
                << "FIRST_TRUE_RHO_PRODUCT="
                << true_rho_product
                << '\n';

            first_example = true;
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

    const double avg_rho_f_true =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_rho_f_true /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_rho_r_true =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_rho_r_true /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_rho_f_control =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_rho_f_control /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_rho_r_control =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_rho_r_control /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_remainder_f_true =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_residual_f_true /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_remainder_r_true =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_residual_r_true /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_remainder_f_control =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_residual_f_control /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_remainder_r_control =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_residual_r_control /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_lcm_f_true =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_lcm_residual_f_true /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_lcm_r_true =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_lcm_residual_r_true /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_lcm_f_control =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_lcm_residual_f_control /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_lcm_r_control =
        control_pairs == 0
            ? 0.0
            : static_cast<double>(
                sum_lcm_residual_r_control /
                static_cast<long double>(
                    control_pairs
                )
            );

    const double avg_true_forward_quotient =
        true_forward_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                sum_true_forward_quotient /
                static_cast<long double>(
                    true_forward_quotient_count
                )
            );

    const double avg_true_reverse_quotient =
        true_reverse_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                sum_true_reverse_quotient /
                static_cast<long double>(
                    true_reverse_quotient_count
                )
            );

    const double avg_control_forward_quotient =
        control_forward_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                sum_control_forward_quotient /
                static_cast<long double>(
                    control_forward_quotient_count
                )
            );

    const double avg_control_reverse_quotient =
        control_reverse_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                sum_control_reverse_quotient /
                static_cast<long double>(
                    control_reverse_quotient_count
                )
            );

    const double avg_rho_sum_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    rho_sum_rank_sum
                ) /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_rho_product_rank =
        true_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    rho_product_rank_sum
                ) /
                static_cast<long double>(
                    true_cases
                )
            );

    const double avg_rho_sum =
        avg_rho_f_true +
        avg_rho_r_true;

    const double avg_rho_product =
        avg_rho_f_true *
        avg_rho_r_true;

    const double avg_rho_sum_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      rho_sum_rank_sum
                  ) /
                  static_cast<double>(
                      true_cases *
                      CONTROL_COUNT
                  )
              );

    const double avg_rho_product_percentile =
        true_cases == 0
            ? 0.0
            : 100.0 *
              (
                  1.0 -
                  static_cast<double>(
                      rho_product_rank_sum
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
        << "AVG_TRUE_RHO_FORWARD="
        << avg_rho_f_true
        << '\n';

    std::cout
        << "AVG_TRUE_RHO_REVERSE="
        << avg_rho_r_true
        << '\n';

    std::cout
        << "AVG_CONTROL_RHO_FORWARD="
        << avg_rho_f_control
        << '\n';

    std::cout
        << "AVG_CONTROL_RHO_REVERSE="
        << avg_rho_r_control
        << '\n';

    std::cout
        << "TRUE_FORWARD_DIVIDES="
        << true_forward_divides
        << '\n';

    std::cout
        << "TRUE_REVERSE_DIVIDES="
        << true_reverse_divides
        << '\n';

    std::cout
        << "TRUE_BOTH_DIVIDE="
        << true_both_divide
        << '\n';

    std::cout
        << "CONTROL_FORWARD_DIVIDES="
        << control_forward_divides
        << '\n';

    std::cout
        << "CONTROL_REVERSE_DIVIDES="
        << control_reverse_divides
        << '\n';

    std::cout
        << "CONTROL_BOTH_DIVIDE="
        << control_both_divide
        << '\n';

    std::cout
        << "TRUE_FORWARD_RHO_GE_50="
        << true_forward_rho_ge_50
        << '\n';

    std::cout
        << "TRUE_REVERSE_RHO_GE_50="
        << true_reverse_rho_ge_50
        << '\n';

    std::cout
        << "TRUE_FORWARD_RHO_GE_25="
        << true_forward_rho_ge_25
        << '\n';

    std::cout
        << "TRUE_REVERSE_RHO_GE_25="
        << true_reverse_rho_ge_25
        << '\n';

    std::cout
        << "CONTROL_FORWARD_RHO_GE_50="
        << control_forward_rho_ge_50
        << '\n';

    std::cout
        << "CONTROL_REVERSE_RHO_GE_50="
        << control_reverse_rho_ge_50
        << '\n';

    std::cout
        << "CONTROL_FORWARD_RHO_GE_25="
        << control_forward_rho_ge_25
        << '\n';

    std::cout
        << "CONTROL_REVERSE_RHO_GE_25="
        << control_reverse_rho_ge_25
        << '\n';

    std::cout
        << "AVG_TRUE_REMAINDER_FORWARD="
        << avg_remainder_f_true
        << '\n';

    std::cout
        << "AVG_TRUE_REMAINDER_REVERSE="
        << avg_remainder_r_true
        << '\n';

    std::cout
        << "AVG_CONTROL_REMAINDER_FORWARD="
        << avg_remainder_f_control
        << '\n';

    std::cout
        << "AVG_CONTROL_REMAINDER_REVERSE="
        << avg_remainder_r_control
        << '\n';

    std::cout
        << "AVG_TRUE_LCM_RESIDUAL_FORWARD="
        << avg_lcm_f_true
        << '\n';

    std::cout
        << "AVG_TRUE_LCM_RESIDUAL_REVERSE="
        << avg_lcm_r_true
        << '\n';

    std::cout
        << "AVG_CONTROL_LCM_RESIDUAL_FORWARD="
        << avg_lcm_f_control
        << '\n';

    std::cout
        << "AVG_CONTROL_LCM_RESIDUAL_REVERSE="
        << avg_lcm_r_control
        << '\n';

    std::cout
        << "TRUE_FORWARD_QUOTIENT_COUNT="
        << true_forward_quotient_count
        << '\n';

    std::cout
        << "TRUE_REVERSE_QUOTIENT_COUNT="
        << true_reverse_quotient_count
        << '\n';

    std::cout
        << "CONTROL_FORWARD_QUOTIENT_COUNT="
        << control_forward_quotient_count
        << '\n';

    std::cout
        << "CONTROL_REVERSE_QUOTIENT_COUNT="
        << control_reverse_quotient_count
        << '\n';

    std::cout
        << "AVG_TRUE_FORWARD_QUOTIENT="
        << avg_true_forward_quotient
        << '\n';

    std::cout
        << "AVG_TRUE_REVERSE_QUOTIENT="
        << avg_true_reverse_quotient
        << '\n';

    std::cout
        << "AVG_CONTROL_FORWARD_QUOTIENT="
        << avg_control_forward_quotient
        << '\n';

    std::cout
        << "AVG_CONTROL_REVERSE_QUOTIENT="
        << avg_control_reverse_quotient
        << '\n';

    std::cout
        << "AVG_RHO_SUM="
        << avg_rho_sum
        << '\n';

    std::cout
        << "AVG_RHO_PRODUCT="
        << avg_rho_product
        << '\n';

    std::cout
        << "AVG_RHO_SUM_RANK="
        << avg_rho_sum_rank
        << '\n';

    std::cout
        << "AVG_RHO_PRODUCT_RANK="
        << avg_rho_product_rank
        << '\n';

    std::cout
        << "AVG_RHO_SUM_PERCENTILE="
        << avg_rho_sum_percentile
        << '\n';

    std::cout
        << "AVG_RHO_PRODUCT_PERCENTILE="
        << avg_rho_product_percentile
        << '\n';

    std::cout
        << "RHO_SUM_STRICT_BEST="
        << rho_sum_true_strict_best
        << '\n';

    std::cout
        << "RHO_SUM_TIED_BEST="
        << rho_sum_true_tied_best
        << '\n';

    std::cout
        << "RHO_PRODUCT_STRICT_BEST="
        << rho_product_true_strict_best
        << '\n';

    std::cout
        << "RHO_PRODUCT_TIED_BEST="
        << rho_product_true_tied_best
        << '\n';

    std::cout
        << "AVG_RHO_SUM_TIES="
        << (
            true_cases == 0
                ? 0.0
                : static_cast<double>(
                    rho_sum_tie_sum
                    /
                    static_cast<long double>(
                        true_cases
                    )
                )
        )
        << '\n';

    std::cout
        << "AVG_RHO_PRODUCT_TIES="
        << (
            true_cases == 0
                ? 0.0
                : static_cast<double>(
                    rho_product_tie_sum
                    /
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
