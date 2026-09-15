#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Witness {
    bool exists = false;
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
};

struct PrimeCase {
    u64 p;
    u64 q;
};

struct Stats {
    u64 cases = 0;
    u64 total_points = 0;

    u64 exact_matches = 0;
    u64 exact_mismatches = 0;

    u64 normalized_ties = 0;

    u64 verification_failures = 0;

    bool first_mismatch_found = false;

    u64 mismatch_r = 0;
    u64 mismatch_k = 0;

    Witness mismatch_exact;
    u64 mismatch_normalized_m = 0;

    u64 exact_m_frequency[33] = {};
    u64 normalized_m_frequency[33] = {};

    u64 mismatch_exact_to_normalized[33][33] = {};
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
        30
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

static Witness make_candidate(
    u64 r,
    u64 m,
    u64 k,
    int sign
) {
    Witness w;

    const u128 mr =
        static_cast<u128>(m) * r;

    u128 value = 0;

    if (sign == -1) {
        if (mr <= 1) {
            return w;
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
        return w;
    }

    const u64 value_u64 =
        static_cast<u64>(value);

    if (k < 2 ||
        value_u64 % k != 0) {
        return w;
    }

    const u64 t =
        value_u64 / k;

    if (t == 0) {
        return w;
    }

    w.exists = true;
    w.m = m;
    w.k = k;
    w.t = t;
    w.sign = sign;

    return w;
}

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

        /*
         * Largest divisor <= K gives
         * smallest t.
         */
        for (u64 k = K;
             k >= 2;
             --k) {

            if (value_u64 % k != 0) {
                continue;
            }

            const u64 t =
                value_u64 / k;

            if (
                !best.exists ||
                t < best.t ||
                (
                    t == best.t &&
                    (
                        k < best.k ||
                        (
                            k == best.k &&
                            sign < best.sign
                        )
                    )
                )
            ) {
                best = make_candidate(
                    r,
                    m,
                    k,
                    sign
                );
            }

            break;
        }
    }

    return best;
}

static Witness exact_winner(
    u64 r,
    u64 K,
    u64 m_max
) {
    Witness best;

    for (u64 m = 1;
         m <= m_max;
         ++m) {

        const Witness w =
            best_for_m(
                r,
                K,
                m
            );

        if (!w.exists) {
            continue;
        }

        if (
            !best.exists ||
            w.t < best.t ||
            (
                w.t == best.t &&
                (
                    w.m < best.m ||
                    (
                        w.m == best.m &&
                        w.k < best.k
                    )
                )
            )
        ) {
            best = w;
        }
    }

    return best;
}

/*
 * Incremental divisor profile.

 * d[m] is the largest k<=K dividing either
 * mr-1 or mr+1.
 *
 * The normalized score is d[m]/m.
 */
static u64 normalized_winner(
    const u64 best_divisor[],
    u64 m_max,
    bool& unique
) {
    unique = true;

    u64 best_m = 0;
    u64 best_d = 0;

    for (u64 m = 1;
         m <= m_max;
         ++m) {

        const u64 d =
            best_divisor[m];

        if (d < 2) {
            continue;
        }

        if (best_m == 0) {
            best_m = m;
            best_d = d;
            continue;
        }

        const u128 lhs =
            static_cast<u128>(d) *
            best_m;

        const u128 rhs =
            static_cast<u128>(best_d) *
            m;

        if (lhs > rhs) {
            best_m = m;
            best_d = d;
            unique = true;
        } else if (lhs == rhs) {
            unique = false;
        }
    }

    return best_m;
}

static PrimeCase generate_case(
    std::mt19937_64& rng
) {
    const u64 p =
        random_prime(
            rng,
            1000,
            1000000
        );

    u64 q =
        random_prime(
            rng,
            1000,
            1000000
        );

    while (q == p) {
        q =
            random_prime(
                rng,
                1000,
                1000000
            );
    }

    return {
        std::min(p, q),
        std::max(p, q)
    };
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    u64 m_max,
    Stats& stats
) {
    ++stats.cases;

    u64 best_divisor[33] = {};

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        /*
         * Update the divisor profile.
         */
        for (u64 m = 1;
             m <= m_max;
             ++m) {

            const u128 mr =
                static_cast<u128>(m) * r;

            const u128 a =
                mr - 1;

            const u128 b =
                mr + 1;

            if (a % K == 0) {
                best_divisor[m] = K;
            }

            if (b % K == 0) {
                best_divisor[m] = K;
            }
        }

        const Witness exact =
            exact_winner(
                r,
                K,
                m_max
            );

        if (!exact.exists) {
            continue;
        }

        if (!verify_witness(r, exact)) {
            ++stats.verification_failures;
        }

        ++stats.total_points;

        ++stats.exact_m_frequency[exact.m];

        bool normalized_unique = true;

        const u64 normalized_m =
            normalized_winner(
                best_divisor,
                m_max,
                normalized_unique
            );

        if (normalized_m == 0) {
            ++stats.exact_mismatches;
            continue;
        }

        ++stats.normalized_m_frequency[
            normalized_m
        ];

        if (!normalized_unique) {
            ++stats.normalized_ties;
        }

        if (
            normalized_unique &&
            normalized_m == exact.m
        ) {
            ++stats.exact_matches;
        } else {
            ++stats.exact_mismatches;

            if (
                exact.m <= m_max &&
                normalized_m <= m_max
            ) {
                ++stats.mismatch_exact_to_normalized[
                    exact.m
                ][
                    normalized_m
                ];
            }

            if (!stats.first_mismatch_found) {
                stats.first_mismatch_found = true;

                stats.mismatch_r = r;
                stats.mismatch_k = K;

                stats.mismatch_exact = exact;
                stats.mismatch_normalized_m =
                    normalized_m;
            }
        }
    }
}

static void print_frequency(
    const char* name,
    const u64 values[],
    u64 max_m
) {
    std::cout
        << "  "
        << name
        << ":\n";

    for (u64 m = 1;
         m <= max_m;
         ++m) {

        if (values[m] == 0) {
            continue;
        }

        std::cout
            << "    M="
            << m
            << " COUNT="
            << values[m]
            << '\n';
    }
}

static void print_stats(
    const Stats& stats,
    u64 m_max
) {
    std::cout
        << "  CASES="
        << stats.cases
        << '\n';

    std::cout
        << "  TOTAL_POINTS="
        << stats.total_points
        << '\n';

    std::cout
        << "  EXACT_WINNER_MATCHES="
        << stats.exact_matches
        << '\n';

    std::cout
        << "  EXACT_WINNER_MISMATCHES="
        << stats.exact_mismatches
        << '\n';

    std::cout
        << "  NORMALIZED_TIES="
        << stats.normalized_ties
        << '\n';

    std::cout
        << "  VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';

    print_frequency(
        "EXACT_WINNER_FREQUENCY",
        stats.exact_m_frequency,
        m_max
    );

    print_frequency(
        "NORMALIZED_WINNER_FREQUENCY",
        stats.normalized_m_frequency,
        m_max
    );

    if (stats.first_mismatch_found) {
        std::cout
            << "  FIRST_MISMATCH:\n";

        std::cout
            << "    R="
            << stats.mismatch_r
            << '\n';

        std::cout
            << "    K="
            << stats.mismatch_k
            << '\n';

        std::cout
            << "    EXACT_M="
            << stats.mismatch_exact.m
            << '\n';

        std::cout
            << "    EXACT_K="
            << stats.mismatch_exact.k
            << '\n';

        std::cout
            << "    EXACT_T="
            << stats.mismatch_exact.t
            << '\n';

        std::cout
            << "    EXACT_SIGN="
            << stats.mismatch_exact.sign
            << '\n';

        std::cout
            << "    NORMALIZED_M="
            << stats.mismatch_normalized_m
            << '\n';
    } else {
        std::cout
            << "  FIRST_MISMATCH=NONE\n";
    }

    bool any_matrix = false;

    for (u64 a = 1;
         a <= m_max;
         ++a) {

        for (u64 b = 1;
             b <= m_max;
             ++b) {

            if (
                stats.mismatch_exact_to_normalized[a][b]
                != 0
            ) {
                any_matrix = true;
            }
        }
    }

    if (any_matrix) {
        std::cout
            << "  MISMATCH_MATRIX:\n";

        for (u64 a = 1;
             a <= m_max;
             ++a) {

            for (u64 b = 1;
                 b <= m_max;
                 ++b) {

                const u64 count =
                    stats.mismatch_exact_to_normalized[a][b];

                if (count == 0) {
                    continue;
                }

                std::cout
                    << "    EXACT_M="
                    << a
                    << " NORMALIZED_M="
                    << b
                    << " COUNT="
                    << count
                    << '\n';
            }
        }
    }
}

int main() {
    constexpr u64 EXPERIMENT = 405;
    constexpr u64 CASES = 1000;
    constexpr u64 K_LIMIT = 200;
    constexpr u64 M_MAX = 32;

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
        << "M_MAX="
        << M_MAX
        << '\n';

    std::cout
        << "PRIME_MIN=1000\n";

    std::cout
        << "PRIME_MAX=1000000\n\n";

    std::mt19937_64 rng(
        405123456789ULL
    );

    std::vector<PrimeCase> cases;
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

    for (const PrimeCase& c : cases) {
        analyze_prime(
            c.p,
            K_LIMIT,
            M_MAX,
            p_stats
        );

        analyze_prime(
            c.q,
            K_LIMIT,
            M_MAX,
            q_stats
        );
    }

    std::cout
        << "\nP_STATS\n";

    print_stats(
        p_stats,
        M_MAX
    );

    std::cout
        << "\nQ_STATS\n";

    print_stats(
        q_stats,
        M_MAX
    );

    const bool pass =
        p_stats.exact_mismatches == 0 &&
        q_stats.exact_mismatches == 0 &&
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
