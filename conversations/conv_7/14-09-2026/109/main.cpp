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

struct Transition {
    bool exists = false;
    u64 k = 0;
    u64 t = 0;
    u64 m = 0;
    int sign = 0;
};

struct Stats {
    u64 cases = 0;

    u64 m1_exists = 0;
    u64 m2_exists = 0;
    u64 m3_exists = 0;

    u64 m2_beats_m1 = 0;
    u64 m3_beats_m1 = 0;
    u64 m3_beats_m2 = 0;

    u64 transition_m2 = 0;
    u64 transition_m3 = 0;
    u64 no_transition = 0;

    u64 min_k_transition_m2 =
        std::numeric_limits<u64>::max();

    u64 max_k_transition_m2 = 0;

    u64 min_k_transition_m3 =
        std::numeric_limits<u64>::max();

    u64 max_k_transition_m3 = 0;

    u64 sum_transition_m2 = 0;
    u64 sum_transition_m3 = 0;

    u64 min_t_m1 =
        std::numeric_limits<u64>::max();

    u64 max_t_m1 = 0;

    u64 min_t_m2 =
        std::numeric_limits<u64>::max();

    u64 max_t_m2 = 0;

    u64 min_t_m3 =
        std::numeric_limits<u64>::max();

    u64 max_t_m3 = 0;

    u64 exact_examples_printed = 0;

    u64 verification_failures = 0;
};

static mpz_class mpz_from_u64(u64 x) {
    return mpz_class(std::to_string(x));
}

static bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    const mpz_class z = mpz_from_u64(n);

    return mpz_probab_prime_p(
        z.get_mpz_t(),
        25
    ) > 0;
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        const u64 r = a % b;
        a = b;
        b = r;
    }

    return a;
}

static u64 integer_sqrt(u64 n) {
    u64 lo = 0;
    u64 hi = 1000000000ULL;

    while (lo <= hi) {
        const u64 mid = lo + (hi - lo) / 2;

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
    std::uniform_int_distribution<u64> dist(lo, hi);

    while (true) {
        const u64 x = dist(rng) | 1ULL;

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
        random_prime(rng, 10000, 1000000);

    u64 b =
        random_prime(rng, 10000, 1000000);

    while (b == a) {
        b =
            random_prime(rng, 10000, 1000000);
    }

    const u64 p = std::min(a, b);
    const u64 q = std::max(a, b);

    const u64 n = p * q;
    const u64 s = integer_sqrt(n);

    return {p, q, n, s};
}

/*
 * Find the best t for a fixed m.

 * We need

 *     k*t = m*r +/- 1

 * with 2 <= k <= K.
 *
 * Scanning k downward finds the largest divisor <= K,
 * hence the smallest t for this m.
 */
static Witness best_for_m(
    u64 r,
    u64 k_max,
    u64 m
) {
    Witness best;

    const u128 mr =
        static_cast<u128>(m) * r;

    for (int sign_index = 0;
         sign_index < 2;
         ++sign_index) {

        const int sign =
            sign_index == 0 ? -1 : +1;

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

        for (u64 k = k_max; k >= 2; --k) {
            if (value_u64 % k != 0) {
                continue;
            }

            const u64 t =
                value_u64 / k;

            if (t == 0) {
                continue;
            }

            /*
             * The modular inverse representative is
             * always <= r/2.
             */
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
        static_cast<u128>(w.k) * w.t;

    const u128 rhs =
        static_cast<u128>(w.m) * r;

    if (w.sign == +1) {
        return lhs == rhs + 1;
    }

    if (w.sign == -1) {
        return lhs + 1 == rhs;
    }

    return false;
}

/*
 * Find the first K for which m=2 can beat m=1.

 * We increment K and compare:
 *
 *     best_{m=1}(K)
 *
 * against
 *
 *     best_{m=2}(K).
 */
static u64 find_m2_transition(
    u64 r,
    u64 k_limit,
    Witness& final_m1,
    Witness& final_m2
) {
    for (u64 k = 2; k <= k_limit; ++k) {
        const Witness w1 =
            best_for_m(r, k, 1);

        const Witness w2 =
            best_for_m(r, k, 2);

        if (!w1.exists || !w2.exists) {
            continue;
        }

        if (w2.t < w1.t) {
            final_m1 = w1;
            final_m2 = w2;
            return k;
        }
    }

    return 0;
}

static u64 find_m3_transition(
    u64 r,
    u64 k_limit,
    Witness& final_m1,
    Witness& final_m2,
    Witness& final_m3
) {
    for (u64 k = 2; k <= k_limit; ++k) {
        const Witness w1 =
            best_for_m(r, k, 1);

        const Witness w2 =
            best_for_m(r, k, 2);

        const Witness w3 =
            best_for_m(r, k, 3);

        if (!w1.exists || !w2.exists || !w3.exists) {
            continue;
        }

        if (
            w3.t < w1.t &&
            w3.t < w2.t
        ) {
            final_m1 = w1;
            final_m2 = w2;
            final_m3 = w3;
            return k;
        }
    }

    return 0;
}

static void analyze_prime(
    u64 r,
    u64 k_limit,
    Stats& stats,
    bool print_examples
) {
    ++stats.cases;

    const Witness m1 =
        best_for_m(r, k_limit, 1);

    const Witness m2 =
        best_for_m(r, k_limit, 2);

    const Witness m3 =
        best_for_m(r, k_limit, 3);

    if (m1.exists) {
        ++stats.m1_exists;

        stats.min_t_m1 =
            std::min(stats.min_t_m1, m1.t);

        stats.max_t_m1 =
            std::max(stats.max_t_m1, m1.t);

        if (!verify_witness(r, m1)) {
            ++stats.verification_failures;
        }
    }

    if (m2.exists) {
        ++stats.m2_exists;

        stats.min_t_m2 =
            std::min(stats.min_t_m2, m2.t);

        stats.max_t_m2 =
            std::max(stats.max_t_m2, m2.t);

        if (!verify_witness(r, m2)) {
            ++stats.verification_failures;
        }

        if (m2.t < m1.t) {
            ++stats.m2_beats_m1;
        }
    }

    if (m3.exists) {
        ++stats.m3_exists;

        stats.min_t_m3 =
            std::min(stats.min_t_m3, m3.t);

        stats.max_t_m3 =
            std::max(stats.max_t_m3, m3.t);

        if (!verify_witness(r, m3)) {
            ++stats.verification_failures;
        }

        if (m3.t < m1.t) {
            ++stats.m3_beats_m1;
        }

        if (m3.t < m2.t) {
            ++stats.m3_beats_m2;
        }
    }

    Witness transition_m1;
    Witness transition_m2;

    const u64 k2 =
        find_m2_transition(
            r,
            k_limit,
            transition_m1,
            transition_m2
        );

    if (k2 != 0) {
        ++stats.transition_m2;

        stats.min_k_transition_m2 =
            std::min(
                stats.min_k_transition_m2,
                k2
            );

        stats.max_k_transition_m2 =
            std::max(
                stats.max_k_transition_m2,
                k2
            );

        stats.sum_transition_m2 += k2;

        if (
            print_examples &&
            stats.exact_examples_printed < 10
        ) {
            ++stats.exact_examples_printed;

            std::cout
                << "EXAMPLE_M2_TRANSITION"
                << " r=" << r
                << " K=" << k2
                << " M1_K=" << transition_m1.k
                << " M1_T=" << transition_m1.t
                << " M1_SIGN=" << transition_m1.sign
                << " M2_K=" << transition_m2.k
                << " M2_T=" << transition_m2.t
                << " M2_SIGN=" << transition_m2.sign
                << '\n';
        }
    } else {
        ++stats.no_transition;
    }

    Witness m3_m1;
    Witness m3_m2;
    Witness m3_m3;

    const u64 k3 =
        find_m3_transition(
            r,
            k_limit,
            m3_m1,
            m3_m2,
            m3_m3
        );

    if (k3 != 0) {
        ++stats.transition_m3;

        stats.min_k_transition_m3 =
            std::min(
                stats.min_k_transition_m3,
                k3
            );

        stats.max_k_transition_m3 =
            std::max(
                stats.max_k_transition_m3,
                k3
            );

        stats.sum_transition_m3 += k3;
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
        << "  M1_EXISTS="
        << stats.m1_exists
        << '\n';

    std::cout
        << "  M2_EXISTS="
        << stats.m2_exists
        << '\n';

    std::cout
        << "  M3_EXISTS="
        << stats.m3_exists
        << '\n';

    std::cout
        << "  M2_BEATS_M1="
        << stats.m2_beats_m1
        << '\n';

    std::cout
        << "  M3_BEATS_M1="
        << stats.m3_beats_m1
        << '\n';

    std::cout
        << "  M3_BEATS_M2="
        << stats.m3_beats_m2
        << '\n';

    std::cout
        << "  M2_TRANSITION_CASES="
        << stats.transition_m2
        << '\n';

    std::cout
        << "  M3_TRANSITION_CASES="
        << stats.transition_m3
        << '\n';

    std::cout
        << "  NO_M2_TRANSITION="
        << stats.no_transition
        << '\n';

    std::cout
        << "  MIN_M2_TRANSITION_K="
        << stats.min_k_transition_m2
        << '\n';

    std::cout
        << "  MAX_M2_TRANSITION_K="
        << stats.max_k_transition_m2
        << '\n';

    if (stats.transition_m2 != 0) {
        std::cout
            << "  AVG_M2_TRANSITION_K="
            << static_cast<double>(
                   stats.sum_transition_m2
               ) /
               static_cast<double>(
                   stats.transition_m2
               )
            << '\n';
    }

    std::cout
        << "  MIN_M3_TRANSITION_K="
        << stats.min_k_transition_m3
        << '\n';

    std::cout
        << "  MAX_M3_TRANSITION_K="
        << stats.max_k_transition_m3
        << '\n';

    if (stats.transition_m3 != 0) {
        std::cout
            << "  AVG_M3_TRANSITION_K="
            << static_cast<double>(
                   stats.sum_transition_m3
               ) /
               static_cast<double>(
                   stats.transition_m3
               )
            << '\n';
    }

    std::cout
        << "  M1_MIN_T="
        << stats.min_t_m1
        << '\n';

    std::cout
        << "  M1_MAX_T="
        << stats.max_t_m1
        << '\n';

    std::cout
        << "  M2_MIN_T="
        << stats.min_t_m2
        << '\n';

    std::cout
        << "  M2_MAX_T="
        << stats.max_t_m2
        << '\n';

    std::cout
        << "  M3_MIN_T="
        << stats.min_t_m3
        << '\n';

    std::cout
        << "  M3_MAX_T="
        << stats.max_t_m3
        << '\n';

    std::cout
        << "  VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';
}

int main() {
    constexpr u64 EXPERIMENT = 400;
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
        400123456789ULL
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

    std::cout
        << "\nSTATUS="
        << (
            p_stats.verification_failures == 0 &&
            q_stats.verification_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
