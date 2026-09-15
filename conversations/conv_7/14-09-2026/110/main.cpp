#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct PrimePair {
    u64 p;
    u64 q;
    u64 n;
    u64 s;
};

struct Witness {
    bool exists = false;

    u64 k = 0;
    u64 t = 0;
    u64 m = 0;

    int sign = 0; // -1 => kt = mr-1, +1 => kt = mr+1
};

struct TransitionRecord {
    bool exists = false;

    u64 K = 0;

    Witness m1;
    Witness m2;

    bool simple_ratio_predicts = false;
    bool exact_condition = false;

    u64 d1 = 0;
    u64 d2 = 0;

    u128 lhs = 0;
    u128 rhs = 0;
};

struct Stats {
    u64 cases = 0;

    u64 transition_cases = 0;
    u64 no_transition_cases = 0;

    u64 simple_ratio_true = 0;
    u64 simple_ratio_false = 0;

    u64 simple_ratio_correct = 0;
    u64 simple_ratio_wrong = 0;

    u64 exact_condition_failures = 0;

    u64 min_transition_K =
        std::numeric_limits<u64>::max();

    u64 max_transition_K = 0;

    u64 sum_transition_K = 0;

    u64 min_d1 =
        std::numeric_limits<u64>::max();

    u64 max_d1 = 0;

    u64 min_d2 =
        std::numeric_limits<u64>::max();

    u64 max_d2 = 0;

    u64 min_d2_over_d1_scaled =
        std::numeric_limits<u64>::max();

    u64 max_d2_over_d1_scaled = 0;

    u64 sum_d2_over_d1_scaled = 0;

    u64 sign_m1_minus = 0;
    u64 sign_m1_plus = 0;

    u64 sign_m2_minus = 0;
    u64 sign_m2_plus = 0;

    u64 verification_failures = 0;

    u64 printed_examples = 0;
};

static mpz_class mpz_from_u64(u64 x) {
    return mpz_class(std::to_string(x));
}

static bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    const mpz_class z =
        mpz_from_u64(n);

    return mpz_probab_prime_p(
        z.get_mpz_t(),
        25
    ) > 0;
}

static u64 integer_sqrt(u64 n) {
    u64 lo = 0;
    u64 hi = 1000000000ULL;

    while (lo <= hi) {
        const u64 mid =
            lo + (hi - lo) / 2;

        const u128 sq =
            static_cast<u128>(mid) * mid;

        if (sq <= n) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return hi;
}

static u64 random_prime(
    std::mt19937_64& rng,
    u64 lo,
    u64 hi
) {
    std::uniform_int_distribution<u64> dist(
        lo,
        hi
    );

    while (true) {
        const u64 x =
            dist(rng) | 1ULL;

        if (x < lo || x > hi) {
            continue;
        }

        if (is_prime(x)) {
            return x;
        }
    }
}

static PrimePair generate_case(
    std::mt19937_64& rng
) {
    const u64 a =
        random_prime(
            rng,
            10000,
            1000000
        );

    u64 b =
        random_prime(
            rng,
            10000,
            1000000
        );

    while (b == a) {
        b =
            random_prime(
                rng,
                10000,
                1000000
            );
    }

    const u64 p =
        std::min(a, b);

    const u64 q =
        std::max(a, b);

    const u64 n = p * q;

    const u64 s =
        integer_sqrt(n);

    return {p, q, n, s};
}

/*
 * Best witness for fixed m and K.

 * We seek
 *
 *     k*t = m*r +/- 1
 *
 * with
 *
 *     2 <= k <= K.
 *
 * For a fixed value m*r +/- 1,
 * scanning k downward finds the largest
 * admissible divisor and therefore the
 * smallest t.
 */
static Witness best_for_m(
    u64 r,
    u64 K,
    u64 m
) {
    Witness best;

    const u128 mr =
        static_cast<u128>(m) * r;

    for (int sign_index = 0;
         sign_index < 2;
         ++sign_index) {

        const int sign =
            sign_index == 0
                ? -1
                : +1;

        u128 value = 0;

        if (sign == -1) {
            if (mr <= 1) {
                continue;
            }

            value = mr - 1;
        } else {
            value = mr + 1;
        }

        if (
            value >
            static_cast<u128>(
                std::numeric_limits<u64>::max()
            )
        ) {
            continue;
        }

        const u64 value_u64 =
            static_cast<u64>(value);

        for (u64 k = K;
             k >= 2;
             --k) {

            if (value_u64 % k != 0) {
                continue;
            }

            const u64 t =
                value_u64 / k;

            if (t == 0) {
                continue;
            }

            if (t > r / 2) {
                continue;
            }

            best.exists = true;
            best.k = k;
            best.t = t;
            best.m = m;
            best.sign = sign;

            break;
        }
    }

    return best;
}

static bool verify_witness(
    u64 r,
    const Witness& w
) {
    if (!w.exists) {
        return false;
    }

    const u128 lhs =
        static_cast<u128>(w.k) *
        w.t;

    const u128 rhs =
        static_cast<u128>(w.m) *
        r;

    if (w.sign == +1) {
        return lhs == rhs + 1;
    }

    if (w.sign == -1) {
        return lhs + 1 == rhs;
    }

    return false;
}

/*
 * Find the first K <= K_LIMIT where m=2
 * beats m=1.
 */
static TransitionRecord find_transition(
    u64 r,
    u64 K_LIMIT
) {
    TransitionRecord result;

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        const Witness w1 =
            best_for_m(r, K, 1);

        const Witness w2 =
            best_for_m(r, K, 2);

        if (!w1.exists || !w2.exists) {
            continue;
        }

        if (w2.t >= w1.t) {
            continue;
        }

        result.exists = true;
        result.K = K;
        result.m1 = w1;
        result.m2 = w2;

        result.d1 = w1.k;
        result.d2 = w2.k;

        /*
         * Exact comparison:
         *
         *     (2r +/- 1)/d2 < (r +/- 1)/d1
         *
         * Cross-multiply.
         *
         * These are exact values because d1,d2
         * actually divide the corresponding numbers.
         */
        const u128 value1 =
            w1.sign == -1
                ? static_cast<u128>(r) - 1
                : static_cast<u128>(r) + 1;

        const u128 value2 =
            w2.sign == -1
                ? static_cast<u128>(2) * r - 1
                : static_cast<u128>(2) * r + 1;

        result.lhs =
            value2 *
            result.d1;

        result.rhs =
            value1 *
            result.d2;

        result.exact_condition =
            result.lhs < result.rhs;

        /*
         * Simple heuristic:
         *
         *     d2 > 2*d1
         */
        result.simple_ratio_predicts =
            result.d2 > 2 * result.d1;

        return result;
    }

    return result;
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    Stats& stats,
    bool print_examples
) {
    ++stats.cases;

    const TransitionRecord tr =
        find_transition(
            r,
            K_LIMIT
        );

    if (!tr.exists) {
        ++stats.no_transition_cases;
        return;
    }

    ++stats.transition_cases;

    stats.min_transition_K =
        std::min(
            stats.min_transition_K,
            tr.K
        );

    stats.max_transition_K =
        std::max(
            stats.max_transition_K,
            tr.K
        );

    stats.sum_transition_K +=
        tr.K;

    stats.min_d1 =
        std::min(
            stats.min_d1,
            tr.d1
        );

    stats.max_d1 =
        std::max(
            stats.max_d1,
            tr.d1
        );

    stats.min_d2 =
        std::min(
            stats.min_d2,
            tr.d2
        );

    stats.max_d2 =
        std::max(
            stats.max_d2,
            tr.d2
        );

    const u64 ratio_scaled =
        static_cast<u64>(
            (
                static_cast<u128>(tr.d2) *
                1000000ULL
            ) /
            tr.d1
        );

    stats.min_d2_over_d1_scaled =
        std::min(
            stats.min_d2_over_d1_scaled,
            ratio_scaled
        );

    stats.max_d2_over_d1_scaled =
        std::max(
            stats.max_d2_over_d1_scaled,
            ratio_scaled
        );

    stats.sum_d2_over_d1_scaled +=
        ratio_scaled;

    if (tr.simple_ratio_predicts) {
        ++stats.simple_ratio_true;
    } else {
        ++stats.simple_ratio_false;
    }

    /*
     * Since this is an actual observed takeover,
     * exact_condition should ALWAYS be true.
     */
    if (!tr.exact_condition) {
        ++stats.exact_condition_failures;
    }

    /*
     * The simple d2 > 2*d1 heuristic is only
     * approximately sufficient because of the +/-1
     * terms.
     */
    if (
        tr.simple_ratio_predicts ==
        tr.exact_condition
    ) {
        ++stats.simple_ratio_correct;
    } else {
        ++stats.simple_ratio_wrong;
    }

    if (tr.m1.sign == -1) {
        ++stats.sign_m1_minus;
    } else {
        ++stats.sign_m1_plus;
    }

    if (tr.m2.sign == -1) {
        ++stats.sign_m2_minus;
    } else {
        ++stats.sign_m2_plus;
    }

    if (!verify_witness(r, tr.m1)) {
        ++stats.verification_failures;
    }

    if (!verify_witness(r, tr.m2)) {
        ++stats.verification_failures;
    }

    if (
        print_examples &&
        stats.printed_examples < 15
    ) {
        ++stats.printed_examples;

        std::cout
            << "EXAMPLE"
            << " r=" << r
            << " K=" << tr.K
            << " d1=" << tr.d1
            << " t1=" << tr.m1.t
            << " s1=" << tr.m1.sign
            << " d2=" << tr.d2
            << " t2=" << tr.m2.t
            << " s2=" << tr.m2.sign
            << " D2_OVER_D1="
            << static_cast<double>(ratio_scaled) /
               1000000.0
            << " SIMPLE="
            << (tr.simple_ratio_predicts ? 1 : 0)
            << '\n';
    }
}

static void print_stats(
    const Stats& stats
) {
    std::cout
        << "  CASES="
        << stats.cases
        << '\n';

    std::cout
        << "  M2_TRANSITION_CASES="
        << stats.transition_cases
        << '\n';

    std::cout
        << "  NO_M2_TRANSITION="
        << stats.no_transition_cases
        << '\n';

    std::cout
        << "  MIN_TRANSITION_K="
        << stats.min_transition_K
        << '\n';

    std::cout
        << "  MAX_TRANSITION_K="
        << stats.max_transition_K
        << '\n';

    if (stats.transition_cases != 0) {
        std::cout
            << "  AVERAGE_TRANSITION_K="
            << static_cast<double>(
                   stats.sum_transition_K
               ) /
               static_cast<double>(
                   stats.transition_cases
               )
            << '\n';

        std::cout
            << "  MIN_D1="
            << stats.min_d1
            << '\n';

        std::cout
            << "  MAX_D1="
            << stats.max_d1
            << '\n';

        std::cout
            << "  MIN_D2="
            << stats.min_d2
            << '\n';

        std::cout
            << "  MAX_D2="
            << stats.max_d2
            << '\n';

        std::cout
            << "  MIN_D2_OVER_D1="
            << static_cast<double>(
                   stats.min_d2_over_d1_scaled
               ) /
               1000000.0
            << '\n';

        std::cout
            << "  MAX_D2_OVER_D1="
            << static_cast<double>(
                   stats.max_d2_over_d1_scaled
               ) /
               1000000.0
            << '\n';

        std::cout
            << "  AVERAGE_D2_OVER_D1="
            << static_cast<double>(
                   stats.sum_d2_over_d1_scaled
               ) /
               static_cast<double>(
                   stats.transition_cases
               ) /
               1000000.0
            << '\n';
    }

    std::cout
        << "  SIMPLE_RATIO_TRUE="
        << stats.simple_ratio_true
        << '\n';

    std::cout
        << "  SIMPLE_RATIO_FALSE="
        << stats.simple_ratio_false
        << '\n';

    std::cout
        << "  SIMPLE_RATIO_CORRECT="
        << stats.simple_ratio_correct
        << '\n';

    std::cout
        << "  SIMPLE_RATIO_WRONG="
        << stats.simple_ratio_wrong
        << '\n';

    std::cout
        << "  M1_SIGN_MINUS="
        << stats.sign_m1_minus
        << '\n';

    std::cout
        << "  M1_SIGN_PLUS="
        << stats.sign_m1_plus
        << '\n';

    std::cout
        << "  M2_SIGN_MINUS="
        << stats.sign_m2_minus
        << '\n';

    std::cout
        << "  M2_SIGN_PLUS="
        << stats.sign_m2_plus
        << '\n';

    std::cout
        << "  EXACT_CONDITION_FAILURES="
        << stats.exact_condition_failures
        << '\n';

    std::cout
        << "  VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';
}

int main() {
    constexpr u64 EXPERIMENT = 401;
    constexpr u64 CASES = 3000;
    constexpr u64 K_LIMIT = 40;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << '\n';

    std::cout
        << "PRIME_MIN=10000\n";

    std::cout
        << "PRIME_MAX=1000000\n\n";

    std::mt19937_64 rng(
        401123456789ULL
    );

    std::vector<PrimePair> cases;
    cases.reserve(CASES);

    for (u64 i = 0; i < CASES; ++i) {
        cases.push_back(
            generate_case(rng)
        );
    }

    Stats p_stats;
    Stats q_stats;

    for (const PrimePair& c : cases) {
        analyze_prime(
            c.p,
            K_LIMIT,
            p_stats,
            true
        );

        analyze_prime(
            c.q,
            K_LIMIT,
            q_stats,
            false
        );
    }

    std::cout
        << "\nP_STATS\n";

    print_stats(p_stats);

    std::cout
        << "\nQ_STATS\n";

    print_stats(q_stats);

    const bool pass =
        p_stats.exact_condition_failures == 0 &&
        q_stats.exact_condition_failures == 0 &&
        p_stats.verification_failures == 0 &&
        q_stats.verification_failures == 0;

    std::cout
        << "\nSTATUS="
        << (pass ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
