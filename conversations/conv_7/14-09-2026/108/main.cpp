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

    int sign = 0; // +1 => kt = mr + 1, -1 => kt = mr - 1
};

struct Stats {
    u64 valid = 0;

    u64 reconstruction_failures = 0;
    u64 verification_failures = 0;

    u64 m1_winner = 0;
    u64 m2_winner = 0;
    u64 m3_winner = 0;
    u64 m4plus_winner = 0;

    u64 m_gt1_beats_m1 = 0;

    u64 sum_t = 0;
    u64 max_t = 0;

    u64 min_m = std::numeric_limits<u64>::max();
    u64 max_m = 0;

    u64 min_k = std::numeric_limits<u64>::max();
    u64 max_k = 0;

    u64 m_frequency[32] = {};
};

static mpz_class mpz_from_u64(u64 x) {
    return mpz_class(std::to_string(x));
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        const u64 r = a % b;
        a = b;
        b = r;
    }

    return a;
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

static bool inverse_mod_gmp(
    u64 a,
    u64 r,
    u64& inverse
) {
    if (r <= 1) {
        return false;
    }

    if (gcd_u64(a, r) != 1) {
        return false;
    }

    const mpz_class aa =
        mpz_from_u64(a);

    const mpz_class rr =
        mpz_from_u64(r);

    mpz_class inv;

    if (
        mpz_invert(
            inv.get_mpz_t(),
            aa.get_mpz_t(),
            rr.get_mpz_t()
        ) == 0
    ) {
        return false;
    }

    inverse = inv.get_ui();

    return inverse != 0;
}

static Witness best_inverse_witness(
    u64 r,
    u64 k_max
) {
    Witness best;

    for (u64 k = 2; k <= k_max; ++k) {
        u64 inv = 0;

        if (!inverse_mod_gmp(k, r, inv)) {
            continue;
        }

        const u64 t =
            std::min(inv, r - inv);

        if (t == 0) {
            continue;
        }

        if (
            !best.exists ||
            t < best.t ||
            (t == best.t && k < best.k)
        ) {
            best.exists = true;
            best.k = k;
            best.t = t;
        }
    }

    if (!best.exists) {
        return best;
    }

    const u128 kt =
        static_cast<u128>(best.k) * best.t;

    const u64 rem =
        static_cast<u64>(kt % r);

    if (rem == 1) {
        best.sign = +1;
        best.m =
            static_cast<u64>((kt - 1) / r);
    } else if (rem == r - 1) {
        best.sign = -1;
        best.m =
            static_cast<u64>((kt + 1) / r);
    } else {
        best.exists = false;
    }

    return best;
}

/*
 * For fixed m and sign, value = m*r +/- 1.
 *
 * Any divisor k <= K gives
 *
 *     t = value / k.
 *
 * For the best t we only need the LARGEST divisor
 * of value which does not exceed K.
 */
static Witness best_divisor_witness(
    u64 r,
    u64 k_max,
    u64 m_max
) {
    Witness best;

    for (u64 m = 1; m <= m_max; ++m) {
        for (int sign_index = 0;
             sign_index < 2;
             ++sign_index) {

            const int sign =
                sign_index == 0 ? -1 : +1;

            const u128 mr =
                static_cast<u128>(m) * r;

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
             * Find the largest divisor <= k_max.
             *
             * Scanning downward is cheap because k_max <= 20.
             */
            for (u64 k = k_max; k >= 2; --k) {
                if (value % k != 0) {
                    continue;
                }

                const u64 t =
                    value / k;

                /*
                 * The modular-inverse representative always
                 * has t <= r/2.
                 */
                if (t > r / 2) {
                    continue;
                }

                if (
                    !best.exists ||
                    t < best.t ||
                    (t == best.t && k < best.k)
                ) {
                    best.exists = true;
                    best.k = k;
                    best.t = t;
                    best.m = m;
                    best.sign = sign;
                }

                /*
                 * Since k is scanned downward, this is the
                 * best divisor for this particular value.
                 */
                break;
            }
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

static void analyze_prime(
    u64 r,
    u64 k_max,
    Stats& stats
) {
    /*
     * Only m <= k_max/2 + 2 can possibly beat the
     * modular-inverse optimum.
     */
    const u64 m_max =
        k_max / 2 + 2;

    const Witness inverse =
        best_inverse_witness(
            r,
            k_max
        );

    const Witness divisor =
        best_divisor_witness(
            r,
            k_max,
            m_max
        );

    if (!inverse.exists ||
        !divisor.exists) {
        ++stats.reconstruction_failures;
        return;
    }

    ++stats.valid;

    if (
        inverse.k != divisor.k ||
        inverse.t != divisor.t ||
        inverse.m != divisor.m ||
        inverse.sign != divisor.sign
    ) {
        ++stats.reconstruction_failures;
    }

    if (!verify_witness(r, divisor)) {
        ++stats.verification_failures;
    }

    if (divisor.m < 32) {
        ++stats.m_frequency[divisor.m];
    }

    stats.min_m =
        std::min(stats.min_m, divisor.m);

    stats.max_m =
        std::max(stats.max_m, divisor.m);

    stats.min_k =
        std::min(stats.min_k, divisor.k);

    stats.max_k =
        std::max(stats.max_k, divisor.k);

    stats.sum_t += divisor.t;

    stats.max_t =
        std::max(stats.max_t, divisor.t);

    if (divisor.m == 1) {
        ++stats.m1_winner;
    } else if (divisor.m == 2) {
        ++stats.m2_winner;
        ++stats.m_gt1_beats_m1;
    } else if (divisor.m == 3) {
        ++stats.m3_winner;
        ++stats.m_gt1_beats_m1;
    } else {
        ++stats.m4plus_winner;
        ++stats.m_gt1_beats_m1;
    }
}

static void print_stats(
    const Stats& stats
) {
    std::cout
        << "    VALID="
        << stats.valid
        << '\n';

    std::cout
        << "    RECONSTRUCTION_FAILURES="
        << stats.reconstruction_failures
        << '\n';

    std::cout
        << "    VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';

    std::cout
        << "    M1_WINNER="
        << stats.m1_winner
        << '\n';

    std::cout
        << "    M2_WINNER="
        << stats.m2_winner
        << '\n';

    std::cout
        << "    M3_WINNER="
        << stats.m3_winner
        << '\n';

    std::cout
        << "    M4PLUS_WINNER="
        << stats.m4plus_winner
        << '\n';

    std::cout
        << "    M_GT1_BEATS_M1="
        << stats.m_gt1_beats_m1
        << '\n';

    std::cout
        << "    MIN_M="
        << stats.min_m
        << '\n';

    std::cout
        << "    MAX_M="
        << stats.max_m
        << '\n';

    std::cout
        << "    MIN_K="
        << stats.min_k
        << '\n';

    std::cout
        << "    MAX_K="
        << stats.max_k
        << '\n';

    std::cout
        << "    AVERAGE_T="
        << static_cast<double>(
               stats.sum_t
           ) /
           static_cast<double>(
               stats.valid
           )
        << '\n';

    std::cout
        << "    MAX_T="
        << stats.max_t
        << '\n';

    std::cout
        << "    M_FREQUENCY:\n";

    for (u64 m = 1; m < 32; ++m) {
        if (stats.m_frequency[m] == 0) {
            continue;
        }

        std::cout
            << "      M="
            << m
            << " COUNT="
            << stats.m_frequency[m]
            << '\n';
    }
}

int main() {
    constexpr u64 EXPERIMENT = 399;
    constexpr u64 CASES = 3000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "K_MIN=2\n";

    std::cout
        << "K_MAX=20\n";

    std::cout
        << "PRIME_MIN=10000\n";

    std::cout
        << "PRIME_MAX=1000000\n\n";

    std::mt19937_64 rng(
        399123456789ULL
    );

    std::vector<PrimePair> cases;
    cases.reserve(CASES);

    for (u64 i = 0; i < CASES; ++i) {
        cases.push_back(
            generate_case(rng)
        );
    }

    for (u64 k_max = 2;
         k_max <= 20;
         ++k_max) {

        Stats p_stats;
        Stats q_stats;

        for (const PrimePair& c : cases) {
            analyze_prime(
                c.p,
                k_max,
                p_stats
            );

            analyze_prime(
                c.q,
                k_max,
                q_stats
            );
        }

        std::cout
            << "K_MAX="
            << k_max
            << '\n';

        std::cout
            << "  P\n";

        print_stats(p_stats);

        std::cout
            << "  Q\n";

        print_stats(q_stats);

        std::cout
            << '\n';
    }

    std::cout
        << "STATUS=PASS_IF_RECONSTRUCTION_FAILURES_AND_VERIFICATION_FAILURES_ARE_ZERO"
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
