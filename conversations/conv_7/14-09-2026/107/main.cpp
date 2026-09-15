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

    u64 multiplier_frequency[32] = {};

    u64 exact_m1 = 0;
    u64 exact_m2 = 0;
    u64 exact_m3 = 0;
    u64 exact_m4_or_more = 0;

    u64 min_m = std::numeric_limits<u64>::max();
    u64 max_m = 0;

    u64 max_k = 0;

    u64 sum_t = 0;
    u64 max_t = 0;
};

static mpz_class mpz_from_u64(u64 x) {
    return mpz_class(std::to_string(x));
}

static bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    mpz_class z = mpz_from_u64(n);

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

/*
 * Exact optimum from the modular-inverse formulation.

 * For each k:
 *
 *     kt == +/-1 (mod r)
 *
 * and the smallest positive t is
 *
 *     min(k^{-1} mod r, r - (k^{-1} mod r)).
 */
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

    const u64 remainder =
        static_cast<u64>(kt % r);

    if (remainder == 1) {
        best.sign = +1;
        best.m =
            static_cast<u64>((kt - 1) / r);
    } else if (remainder == r - 1) {
        best.sign = -1;
        best.m =
            static_cast<u64>((kt + 1) / r);
    } else {
        best.exists = false;
    }

    return best;
}

/*
 * Independent reconstruction.

 * We seek
 *
 *     k*t = m*r +/- 1
 *
 * with
 *
 *     2 <= k <= k_max.
 *
 * The exact winning solution always has t <= r/2.
 *
 * We search m up to k_max/2 + 2, which safely covers
 * the possible winning multiplier.
 */
static Witness best_multiplier_witness(
    u64 r,
    u64 k_max
) {
    Witness best;

    const u64 m_max =
        k_max / 2 + 2;

    for (u64 m = 1; m <= m_max; ++m) {

        for (int sign_index = 0;
             sign_index < 2;
             ++sign_index) {

            const int sign =
                sign_index == 0 ? -1 : +1;

            /*
             * m*r +/- 1.
             *
             * m <= ~11 and r <= 1e6 here,
             * so uint64 is more than sufficient.
             */
            const u128 mr =
                static_cast<u128>(m) * r;

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

            for (u64 k = 2;
                 k <= k_max;
                 ++k) {

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

    const u128 kt =
        static_cast<u128>(w.k) * w.t;

    const u128 mr =
        static_cast<u128>(w.m) * r;

    if (w.sign == +1) {
        return kt == mr + 1;
    }

    if (w.sign == -1) {
        return kt + 1 == mr;
    }

    return false;
}

static void analyze_prime(
    u64 r,
    u64 k_max,
    Stats& stats
) {
    const Witness exact =
        best_inverse_witness(r, k_max);

    const Witness reconstructed =
        best_multiplier_witness(r, k_max);

    if (!exact.exists ||
        !reconstructed.exists) {
        ++stats.reconstruction_failures;
        return;
    }

    ++stats.valid;

    if (
        exact.k != reconstructed.k ||
        exact.t != reconstructed.t ||
        exact.m != reconstructed.m ||
        exact.sign != reconstructed.sign
    ) {
        ++stats.reconstruction_failures;
    }

    if (!verify_witness(r, reconstructed)) {
        ++stats.verification_failures;
    }

    if (exact.m < 32) {
        ++stats.multiplier_frequency[exact.m];
    }

    if (exact.m == 1) {
        ++stats.exact_m1;
    } else if (exact.m == 2) {
        ++stats.exact_m2;
    } else if (exact.m == 3) {
        ++stats.exact_m3;
    } else {
        ++stats.exact_m4_or_more;
    }

    stats.min_m =
        std::min(stats.min_m, exact.m);

    stats.max_m =
        std::max(stats.max_m, exact.m);

    stats.max_k =
        std::max(stats.max_k, exact.k);

    stats.sum_t += exact.t;

    stats.max_t =
        std::max(stats.max_t, exact.t);
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
        << "    MIN_M="
        << stats.min_m
        << '\n';

    std::cout
        << "    MAX_M="
        << stats.max_m
        << '\n';

    std::cout
        << "    M1="
        << stats.exact_m1
        << '\n';

    std::cout
        << "    M2="
        << stats.exact_m2
        << '\n';

    std::cout
        << "    M3="
        << stats.exact_m3
        << '\n';

    std::cout
        << "    M4_OR_MORE="
        << stats.exact_m4_or_more
        << '\n';

    std::cout
        << "    MAX_WINNING_K="
        << stats.max_k
        << '\n';

    std::cout
        << "    AVERAGE_T="
        << static_cast<double>(stats.sum_t) /
           static_cast<double>(stats.valid)
        << '\n';

    std::cout
        << "    MAX_T="
        << stats.max_t
        << '\n';

    std::cout
        << "    M_FREQUENCY:\n";

    for (u64 m = 1; m < 32; ++m) {
        if (stats.multiplier_frequency[m] == 0) {
            continue;
        }

        std::cout
            << "      M="
            << m
            << " COUNT="
            << stats.multiplier_frequency[m]
            << '\n';
    }
}

int main() {
    constexpr u64 EXPERIMENT = 398;
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
        398123456789ULL
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