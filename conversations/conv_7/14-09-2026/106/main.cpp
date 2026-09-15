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
    u64 kt = 0;
    int sign = 0;          // -1: kt = mr-1, +1: kt = mr+1
    u64 multiplier = 0;
};

struct Stats {
    u64 valid = 0;
    u64 verification_failures = 0;

    u64 lower_bound_hits = 0;
    u64 lower_bound_misses = 0;

    u64 multiplier_1 = 0;
    u64 multiplier_gt1 = 0;

    u64 exact_r_minus_1 = 0;
    u64 exact_r_plus_1 = 0;

    u64 sum_t = 0;
    u64 min_t = std::numeric_limits<u64>::max();
    u64 max_t = 0;

    u64 sum_gap = 0;
    u64 min_gap = std::numeric_limits<u64>::max();
    u64 max_gap = 0;

    u64 sum_ratio_scaled = 0;
    u64 max_ratio_scaled = 0;

    u64 sum_t_over_lower_scaled = 0;
    u64 max_t_over_lower_scaled = 0;

    u64 k_frequency[21] = {};
};

static mpz_class mpz_from_u64(u64 x) {
    return mpz_class(std::to_string(x));
}

static u64 mpz_to_u64(const mpz_class& x) {
    return x.get_ui();
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
    const u64 a = random_prime(
        rng,
        10000,
        1000000
    );

    u64 b = random_prime(
        rng,
        10000,
        1000000
    );

    while (b == a) {
        b = random_prime(
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

    const mpz_class aa = mpz_from_u64(a);
    const mpz_class rr = mpz_from_u64(r);

    mpz_class inv;

    if (mpz_invert(
            inv.get_mpz_t(),
            aa.get_mpz_t(),
            rr.get_mpz_t()
        ) == 0) {
        return false;
    }

    inverse = mpz_to_u64(inv);

    return inverse != 0;
}

static Witness best_witness(
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

        const u128 kt128 =
            static_cast<u128>(k) * t;

        if (
            kt128 >
            static_cast<u128>(
                std::numeric_limits<u64>::max()
            )
        ) {
            continue;
        }

        const u64 kt =
            static_cast<u64>(kt128);

        int sign = 0;
        u64 multiplier = 0;

        if (kt % r == r - 1) {
            sign = -1;
            multiplier = (kt + 1) / r;
        } else if (kt % r == 1) {
            sign = +1;
            multiplier = (kt - 1) / r;
        } else {
            std::cerr
                << "INTERNAL_INVERSE_ERROR r="
                << r
                << " k="
                << k
                << " inv="
                << inv
                << " t="
                << t
                << " kt="
                << kt
                << '\n';

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
            best.kt = kt;
            best.sign = sign;
            best.multiplier = multiplier;
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

    if (w.kt % r == 1) {
        return true;
    }

    if (w.kt % r == r - 1) {
        return true;
    }

    return false;
}

static u64 ceil_div(
    u64 a,
    u64 b
) {
    return (a + b - 1) / b;
}

static void analyze_prime(
    u64 r,
    u64 k_max,
    Stats& stats
) {
    const Witness w =
        best_witness(r, k_max);

    if (!w.exists) {
        return;
    }

    ++stats.valid;

    if (!verify_witness(r, w)) {
        ++stats.verification_failures;
    }

    ++stats.k_frequency[w.k];

    stats.sum_t += w.t;

    stats.min_t =
        std::min(stats.min_t, w.t);

    stats.max_t =
        std::max(stats.max_t, w.t);

    const u64 lower =
        ceil_div(r - 1, k_max);

    if (w.t == lower) {
        ++stats.lower_bound_hits;
    } else {
        ++stats.lower_bound_misses;
    }

    if (w.multiplier == 1) {
        ++stats.multiplier_1;
    } else {
        ++stats.multiplier_gt1;
    }

    if (
        w.sign == -1 &&
        w.kt == r - 1
    ) {
        ++stats.exact_r_minus_1;
    }

    if (
        w.sign == +1 &&
        w.kt == r + 1
    ) {
        ++stats.exact_r_plus_1;
    }

    const u64 gap =
        w.kt - (r - 1);

    stats.sum_gap += gap;

    stats.min_gap =
        std::min(stats.min_gap, gap);

    stats.max_gap =
        std::max(stats.max_gap, gap);

    const u64 ratio_scaled =
        static_cast<u64>(
            (
                static_cast<u128>(w.kt) *
                1000000ULL
            ) / r
        );

    stats.sum_ratio_scaled +=
        ratio_scaled;

    stats.max_ratio_scaled =
        std::max(
            stats.max_ratio_scaled,
            ratio_scaled
        );

    const u64 t_over_lower_scaled =
        static_cast<u64>(
            (
                static_cast<u128>(w.t) *
                1000000ULL
            ) / lower
        );

    stats.sum_t_over_lower_scaled +=
        t_over_lower_scaled;

    stats.max_t_over_lower_scaled =
        std::max(
            stats.max_t_over_lower_scaled,
            t_over_lower_scaled
        );
}

static void print_stats(
    const char* label,
    const Stats& stats
) {
    std::cout << "  " << label << '\n';

    std::cout
        << "    VALID="
        << stats.valid
        << '\n';

    std::cout
        << "    ROOT_VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';

    std::cout
        << "    LOWER_BOUND_HITS="
        << stats.lower_bound_hits
        << '\n';

    std::cout
        << "    LOWER_BOUND_MISSES="
        << stats.lower_bound_misses
        << '\n';

    std::cout
        << "    MULTIPLIER_1="
        << stats.multiplier_1
        << '\n';

    std::cout
        << "    MULTIPLIER_GT1="
        << stats.multiplier_gt1
        << '\n';

    std::cout
        << "    EXACT_R_MINUS_1="
        << stats.exact_r_minus_1
        << '\n';

    std::cout
        << "    EXACT_R_PLUS_1="
        << stats.exact_r_plus_1
        << '\n';

    std::cout
        << "    MIN_T="
        << stats.min_t
        << '\n';

    std::cout
        << "    MAX_T="
        << stats.max_t
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
        << "    MIN_GAP="
        << stats.min_gap
        << '\n';

    std::cout
        << "    MAX_GAP="
        << stats.max_gap
        << '\n';

    std::cout
        << "    AVERAGE_GAP="
        << static_cast<double>(
               stats.sum_gap
           ) /
           static_cast<double>(
               stats.valid
           )
        << '\n';

    std::cout
        << "    MAX_KT_OVER_R="
        << static_cast<double>(
               stats.max_ratio_scaled
           ) /
           1000000.0
        << '\n';

    std::cout
        << "    AVERAGE_KT_OVER_R="
        << static_cast<double>(
               stats.sum_ratio_scaled
           ) /
           static_cast<double>(
               stats.valid
           ) /
           1000000.0
        << '\n';

    std::cout
        << "    MAX_T_OVER_LOWER="
        << static_cast<double>(
               stats.max_t_over_lower_scaled
           ) /
           1000000.0
        << '\n';

    std::cout
        << "    AVERAGE_T_OVER_LOWER="
        << static_cast<double>(
               stats.sum_t_over_lower_scaled
           ) /
           static_cast<double>(
               stats.valid
           ) /
           1000000.0
        << '\n';

    std::cout
        << "    K_FREQUENCY:\n";

    for (u64 k = 2; k <= 20; ++k) {
        if (stats.k_frequency[k] == 0) {
            continue;
        }

        std::cout
            << "      K="
            << k
            << " COUNT="
            << stats.k_frequency[k]
            << '\n';
    }
}

int main() {
    constexpr u64 EXPERIMENT = 397;
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
        397123456789ULL
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

        print_stats(
            "P",
            p_stats
        );

        print_stats(
            "Q",
            q_stats
        );

        std::cout << '\n';
    }

    std::cout
        << "STATUS=PASS_IF_ROOT_VERIFICATION_FAILURES_ARE_ZERO"
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}