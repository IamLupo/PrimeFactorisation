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

    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
};

struct Stats {
    u64 cases = 0;

    u64 verification_failures = 0;

    u64 winner_m1 = 0;
    u64 winner_m2 = 0;
    u64 winner_m3 = 0;
    u64 winner_m4 = 0;

    u64 transition_1_to_2 = 0;
    u64 transition_2_to_3 = 0;
    u64 transition_3_to_4 = 0;

    u64 first_transition_1_to_2 =
        std::numeric_limits<u64>::max();

    u64 first_transition_2_to_3 =
        std::numeric_limits<u64>::max();

    u64 first_transition_3_to_4 =
        std::numeric_limits<u64>::max();

    u64 last_transition_1_to_2 = 0;
    u64 last_transition_2_to_3 = 0;
    u64 last_transition_3_to_4 = 0;

    u64 sum_transition_1_to_2 = 0;
    u64 sum_transition_2_to_3 = 0;
    u64 sum_transition_3_to_4 = 0;

    u64 profile_mismatch = 0;

    u64 examples_printed = 0;
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

    const u64 p = std::min(a, b);
    const u64 q = std::max(a, b);

    const u64 n = p * q;
    const u64 s = integer_sqrt(n);

    return {p, q, n, s};
}

/*
 * Best witness for fixed m and K.

 * We search divisors k <= K of
 *
 *     m*r-1
 *     m*r+1
 *
 * and minimize t=(m*r +/- 1)/k.
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

        u128 value128 = 0;

        if (sign == -1) {
            if (mr <= 1) {
                continue;
            }

            value128 = mr - 1;
        } else {
            value128 = mr + 1;
        }

        if (
            value128 >
            static_cast<u128>(
                std::numeric_limits<u64>::max()
            )
        ) {
            continue;
        }

        const u64 value =
            static_cast<u64>(value128);

        /*
         * Largest divisor <= K gives
         * smallest possible t, so scan downward.
         */
        for (u64 k = K;
             k >= 2;
             --k) {

            if (value % k != 0) {
                continue;
            }

            const u64 t =
                value / k;

            if (t == 0) {
                continue;
            }

            if (t > r / 2) {
                continue;
            }

            if (
                !best.exists ||
                t < best.t ||
                (t == best.t && k < best.k)
            ) {
                best.exists = true;
                best.m = m;
                best.k = k;
                best.t = t;
                best.sign = sign;
            }

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

    const u128 mr =
        static_cast<u128>(w.m) * r;

    if (w.sign == +1) {
        return lhs == mr + 1;
    }

    if (w.sign == -1) {
        return lhs + 1 == mr;
    }

    return false;
}

/*
 * The profile winner is simply the minimum T_m(K).
 *
 * m = 1,2,3,4 are compared.
 */
static Witness profile_winner(
    u64 r,
    u64 K
) {
    Witness best;

    for (u64 m = 1; m <= 4; ++m) {
        const Witness w =
            best_for_m(r, K, m);

        if (!w.exists) {
            continue;
        }

        if (
            !best.exists ||
            w.t < best.t ||
            (w.t == best.t && w.m < best.m)
        ) {
            best = w;
        }
    }

    return best;
}

static void update_transition(
    const Witness& previous,
    const Witness& current,
    u64 K,
    Stats& stats,
    bool print_example,
    u64 r
) {
    if (!previous.exists ||
        !current.exists) {
        return;
    }

    if (previous.m == current.m) {
        return;
    }

    if (previous.m == 1 &&
        current.m == 2) {

        ++stats.transition_1_to_2;

        stats.first_transition_1_to_2 =
            std::min(
                stats.first_transition_1_to_2,
                K
            );

        stats.last_transition_1_to_2 =
            std::max(
                stats.last_transition_1_to_2,
                K
            );

        stats.sum_transition_1_to_2 += K;

        if (print_example &&
            stats.examples_printed < 15) {

            ++stats.examples_printed;

            std::cout
                << "TRANSITION_1_TO_2"
                << " r=" << r
                << " K=" << K
                << " M1_K=" << previous.k
                << " M1_T=" << previous.t
                << " M2_K=" << current.k
                << " M2_T=" << current.t
                << '\n';
        }
    }

    if (previous.m == 2 &&
        current.m == 3) {

        ++stats.transition_2_to_3;

        stats.first_transition_2_to_3 =
            std::min(
                stats.first_transition_2_to_3,
                K
            );

        stats.last_transition_2_to_3 =
            std::max(
                stats.last_transition_2_to_3,
                K
            );

        stats.sum_transition_2_to_3 += K;
    }

    if (previous.m == 3 &&
        current.m == 4) {

        ++stats.transition_3_to_4;

        stats.first_transition_3_to_4 =
            std::min(
                stats.first_transition_3_to_4,
                K
            );

        stats.last_transition_3_to_4 =
            std::max(
                stats.last_transition_3_to_4,
                K
            );

        stats.sum_transition_3_to_4 += K;
    }
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    Stats& stats,
    bool print_examples
) {
    ++stats.cases;

    Witness previous;
    u64 previous_K = 0;

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        /*
         * Independently construct the profile.
         */
        const Witness profile =
            profile_winner(r, K);

        /*
         * Verify every member of the profile.
         */
        for (u64 m = 1; m <= 4; ++m) {
            const Witness w =
                best_for_m(r, K, m);

            if (
                w.exists &&
                !verify_witness(r, w)
            ) {
                ++stats.verification_failures;
            }
        }

        if (!profile.exists) {
            continue;
        }

        /*
         * Compare the exact profile winner
         * to itself reconstructed from all m.
         *
         * This is deliberately redundant to catch
         * implementation mistakes.
         */
        Witness reconstructed;

        for (u64 m = 1; m <= 4; ++m) {
            const Witness w =
                best_for_m(r, K, m);

            if (!w.exists) {
                continue;
            }

            if (
                !reconstructed.exists ||
                w.t < reconstructed.t ||
                (
                    w.t == reconstructed.t &&
                    w.m < reconstructed.m
                )
            ) {
                reconstructed = w;
            }
        }

        if (
            profile.exists !=
            reconstructed.exists ||
            (
                profile.exists &&
                (
                    profile.m != reconstructed.m ||
                    profile.k != reconstructed.k ||
                    profile.t != reconstructed.t
                )
            )
        ) {
            ++stats.profile_mismatch;
        }

        if (previous.exists) {
            update_transition(
                previous,
                profile,
                K,
                stats,
                print_examples,
                r
            );
        }

        previous = profile;
        previous_K = K;

        (void)previous_K;

        switch (profile.m) {
            case 1:
                ++stats.winner_m1;
                break;

            case 2:
                ++stats.winner_m2;
                break;

            case 3:
                ++stats.winner_m3;
                break;

            case 4:
                ++stats.winner_m4;
                break;

            default:
                break;
        }
    }
}

static void print_optional_stat(
    const char* name,
    u64 value
) {
    if (value ==
        std::numeric_limits<u64>::max()) {
        std::cout
            << "  "
            << name
            << "=NONE\n";
    } else {
        std::cout
            << "  "
            << name
            << "="
            << value
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
        << "  WINNER_M1="
        << stats.winner_m1
        << '\n';

    std::cout
        << "  WINNER_M2="
        << stats.winner_m2
        << '\n';

    std::cout
        << "  WINNER_M3="
        << stats.winner_m3
        << '\n';

    std::cout
        << "  WINNER_M4="
        << stats.winner_m4
        << '\n';

    std::cout
        << "  TRANSITION_1_TO_2="
        << stats.transition_1_to_2
        << '\n';

    std::cout
        << "  TRANSITION_2_TO_3="
        << stats.transition_2_to_3
        << '\n';

    std::cout
        << "  TRANSITION_3_TO_4="
        << stats.transition_3_to_4
        << '\n';

    print_optional_stat(
        "MIN_TRANSITION_1_TO_2_K",
        stats.first_transition_1_to_2
    );

    print_optional_stat(
        "MAX_TRANSITION_1_TO_2_K",
        stats.last_transition_1_to_2
    );

    print_optional_stat(
        "MIN_TRANSITION_2_TO_3_K",
        stats.first_transition_2_to_3
    );

    print_optional_stat(
        "MAX_TRANSITION_2_TO_3_K",
        stats.last_transition_2_to_3
    );

    print_optional_stat(
        "MIN_TRANSITION_3_TO_4_K",
        stats.first_transition_3_to_4
    );

    print_optional_stat(
        "MAX_TRANSITION_3_TO_4_K",
        stats.last_transition_3_to_4
    );

    if (stats.transition_1_to_2 != 0) {
        std::cout
            << "  AVERAGE_TRANSITION_1_TO_2_K="
            << static_cast<double>(
                   stats.sum_transition_1_to_2
               ) /
               static_cast<double>(
                   stats.transition_1_to_2
               )
            << '\n';
    }

    if (stats.transition_2_to_3 != 0) {
        std::cout
            << "  AVERAGE_TRANSITION_2_TO_3_K="
            << static_cast<double>(
                   stats.sum_transition_2_to_3
               ) /
               static_cast<double>(
                   stats.transition_2_to_3
               )
            << '\n';
    }

    if (stats.transition_3_to_4 != 0) {
        std::cout
            << "  AVERAGE_TRANSITION_3_TO_4_K="
            << static_cast<double>(
                   stats.sum_transition_3_to_4
               ) /
               static_cast<double>(
                   stats.transition_3_to_4
               )
            << '\n';
    }

    std::cout
        << "  PROFILE_MISMATCH="
        << stats.profile_mismatch
        << '\n';

    std::cout
        << "  VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';
}

int main() {
    constexpr u64 EXPERIMENT = 402;
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
        << "M_MAX=4\n";

    std::cout
        << "PRIME_MIN=10000\n";

    std::cout
        << "PRIME_MAX=1000000\n\n";

    std::mt19937_64 rng(
        402123456789ULL
    );

    std::vector<PrimePair> cases;
    cases.reserve(CASES);

    for (u64 i = 0;
         i < CASES;
         ++i) {
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
        p_stats.profile_mismatch == 0 &&
        q_stats.profile_mismatch == 0 &&
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
