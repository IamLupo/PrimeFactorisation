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

struct FamilyState {
    u64 best_k_minus = 0;
    u64 best_k_plus = 0;
};

struct Stats {
    u64 cases = 0;
    u64 total_points = 0;

    u64 exact_winner_matches = 0;
    u64 exact_winner_mismatches = 0;

    u64 first_mismatch_points = 0;

    u64 winner_m1 = 0;
    u64 winner_m2 = 0;
    u64 winner_m3 = 0;
    u64 winner_m4 = 0;
    u64 winner_m5 = 0;
    u64 winner_m6 = 0;
    u64 winner_m7 = 0;
    u64 winner_m8 = 0;

    u64 normalized_winner_m1 = 0;
    u64 normalized_winner_m2 = 0;
    u64 normalized_winner_m3 = 0;
    u64 normalized_winner_m4 = 0;
    u64 normalized_winner_m5 = 0;
    u64 normalized_winner_m6 = 0;
    u64 normalized_winner_m7 = 0;
    u64 normalized_winner_m8 = 0;

    u64 mismatch_by_k[81] = {};

    u64 mismatch_exact_m[9][9] = {};

    u64 normalized_tie_points = 0;

    u64 normalized_tie_exact_winner = 0;
    u64 normalized_tie_exact_winner_not_unique = 0;

    u64 max_normalized_score_scaled = 0;

    u64 verification_failures = 0;
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

    const u64 n =
        p * q;

    const u64 s =
        integer_sqrt(n);

    return {p, q, n, s};
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

    if (t == 0 ||
        t > r / 2) {
        return w;
    }

    w.exists = true;
    w.m = m;
    w.k = k;
    w.t = t;
    w.sign = sign;

    return w;
}

static Witness exact_family_witness(
    u64 r,
    u64 m,
    const FamilyState& state
) {
    Witness best;

    if (state.best_k_minus >= 2) {
        best =
            make_candidate(
                r,
                m,
                state.best_k_minus,
                -1
            );
    }

    if (state.best_k_plus >= 2) {
        const Witness candidate =
            make_candidate(
                r,
                m,
                state.best_k_plus,
                +1
            );

        if (
            candidate.exists &&
            (
                !best.exists ||
                candidate.t < best.t ||
                (
                    candidate.t == best.t &&
                    candidate.k < best.k
                )
            )
        ) {
            best = candidate;
        }
    }

    return best;
}

static Witness exact_global_winner(
    u64 r,
    const FamilyState states[],
    u64 m_max
) {
    Witness best;

    for (u64 m = 1;
         m <= m_max;
         ++m) {

        const Witness w =
            exact_family_witness(
                r,
                m,
                states[m]
            );

        if (!w.exists) {
            continue;
        }

        if (
            !best.exists ||
            w.t < best.t ||
            (
                w.t == best.t &&
                w.m < best.m
            )
        ) {
            best = w;
        }
    }

    return best;
}

/*
 * Normalized divisor score:
 *
 *     S_m = max(d_minus, d_plus) / m
 *
 * The comparison is performed exactly through
 *
 *     d_a / m_a > d_b / m_b
 *
 * using cross multiplication.
 */
static u64 normalized_winner(
    const FamilyState states[],
    u64 m_max,
    bool& unique
) {
    unique = true;

    u64 best_m = 0;
    u64 best_k = 0;

    for (u64 m = 1;
         m <= m_max;
         ++m) {

        const u64 d =
            std::max(
                states[m].best_k_minus,
                states[m].best_k_plus
            );

        if (d < 2) {
            continue;
        }

        if (best_m == 0) {
            best_m = m;
            best_k = d;
            continue;
        }

        const u128 lhs =
            static_cast<u128>(d) *
            best_m;

        const u128 rhs =
            static_cast<u128>(best_k) *
            m;

        if (lhs > rhs) {
            best_m = m;
            best_k = d;
            unique = true;
        } else if (lhs == rhs) {
            unique = false;

            /*
             * Keep the smaller m as deterministic
             * representative.
             */
        }
    }

    return best_m;
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

static void count_exact_winner(
    Stats& stats,
    u64 m
) {
    switch (m) {
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

        case 5:
            ++stats.winner_m5;
            break;

        case 6:
            ++stats.winner_m6;
            break;

        case 7:
            ++stats.winner_m7;
            break;

        case 8:
            ++stats.winner_m8;
            break;

        default:
            break;
    }
}

static void count_normalized_winner(
    Stats& stats,
    u64 m
) {
    switch (m) {
        case 1:
            ++stats.normalized_winner_m1;
            break;

        case 2:
            ++stats.normalized_winner_m2;
            break;

        case 3:
            ++stats.normalized_winner_m3;
            break;

        case 4:
            ++stats.normalized_winner_m4;
            break;

        case 5:
            ++stats.normalized_winner_m5;
            break;

        case 6:
            ++stats.normalized_winner_m6;
            break;

        case 7:
            ++stats.normalized_winner_m7;
            break;

        case 8:
            ++stats.normalized_winner_m8;
            break;

        default:
            break;
    }
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    u64 m_max,
    Stats& stats
) {
    ++stats.cases;

    FamilyState states[9] = {};

    bool saw_mismatch =
        false;

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        /*
         * A new divisor can only appear when
         * K itself divides mr-1 or mr+1.
         */
        for (u64 m = 1;
             m <= m_max;
             ++m) {

            const u128 mr =
                static_cast<u128>(m) * r;

            const u128 minus_value =
                mr - 1;

            const u128 plus_value =
                mr + 1;

            if (
                minus_value % K == 0 &&
                K >= 2
            ) {
                states[m].best_k_minus = K;
            }

            if (
                plus_value % K == 0 &&
                K >= 2
            ) {
                states[m].best_k_plus = K;
            }
        }

        const Witness exact =
            exact_global_winner(
                r,
                states,
                m_max
            );

        bool unique_normalized = true;

        const u64 normalized_m =
            normalized_winner(
                states,
                m_max,
                unique_normalized
            );

        if (!exact.exists ||
            normalized_m == 0) {
            continue;
        }

        ++stats.total_points;

        count_exact_winner(
            stats,
            exact.m
        );

        count_normalized_winner(
            stats,
            normalized_m
        );

        if (!unique_normalized) {
            ++stats.normalized_tie_points;
        }

        if (
            unique_normalized &&
            normalized_m == exact.m
        ) {
            ++stats.exact_winner_matches;
        } else {
            ++stats.exact_winner_mismatches;

            if (!saw_mismatch) {
                saw_mismatch = true;
                ++stats.first_mismatch_points;
            }

            if (
                exact.m <= 8 &&
                normalized_m <= 8
            ) {
                ++stats.mismatch_exact_m
                    [exact.m]
                    [normalized_m];
            }
        }

        /*
         * Score = d/m.
         *
         * Use the largest divisor currently available.
         */
        u64 best_score_m = 0;
        u64 best_score_d = 0;

        for (u64 m = 1;
             m <= m_max;
             ++m) {

            const u64 d =
                std::max(
                    states[m].best_k_minus,
                    states[m].best_k_plus
                );

            if (d < 2) {
                continue;
            }

            if (best_score_m == 0) {
                best_score_m = m;
                best_score_d = d;
                continue;
            }

            const u128 lhs =
                static_cast<u128>(d) *
                best_score_m;

            const u128 rhs =
                static_cast<u128>(best_score_d) *
                m;

            if (lhs > rhs) {
                best_score_m = m;
                best_score_d = d;
            }
        }

        if (best_score_m != 0) {
            const u64 score_scaled =
                static_cast<u64>(
                    (
                        static_cast<u128>(
                            best_score_d
                        ) *
                        1000000ULL
                    ) /
                    best_score_m
                );

            stats.max_normalized_score_scaled =
                std::max(
                    stats.max_normalized_score_scaled,
                    score_scaled
                );
        }

        if (!verify_witness(r, exact)) {
            ++stats.verification_failures;
        }
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
        << "  TOTAL_POINTS="
        << stats.total_points
        << '\n';

    std::cout
        << "  EXACT_WINNER_MATCHES="
        << stats.exact_winner_matches
        << '\n';

    std::cout
        << "  EXACT_WINNER_MISMATCHES="
        << stats.exact_winner_mismatches
        << '\n';

    std::cout
        << "  FIRST_MISMATCH_CASES="
        << stats.first_mismatch_points
        << '\n';

    std::cout
        << "  NORMALIZED_TIE_POINTS="
        << stats.normalized_tie_points
        << '\n';

    std::cout
        << "  VERIFICATION_FAILURES="
        << stats.verification_failures
        << '\n';

    std::cout
        << "  EXACT_WINNER_FREQUENCY:\n";

    std::cout
        << "    M1="
        << stats.winner_m1
        << '\n';

    std::cout
        << "    M2="
        << stats.winner_m2
        << '\n';

    std::cout
        << "    M3="
        << stats.winner_m3
        << '\n';

    std::cout
        << "    M4="
        << stats.winner_m4
        << '\n';

    std::cout
        << "    M5="
        << stats.winner_m5
        << '\n';

    std::cout
        << "    M6="
        << stats.winner_m6
        << '\n';

    std::cout
        << "    M7="
        << stats.winner_m7
        << '\n';

    std::cout
        << "    M8="
        << stats.winner_m8
        << '\n';

    std::cout
        << "  NORMALIZED_WINNER_FREQUENCY:\n";

    std::cout
        << "    M1="
        << stats.normalized_winner_m1
        << '\n';

    std::cout
        << "    M2="
        << stats.normalized_winner_m2
        << '\n';

    std::cout
        << "    M3="
        << stats.normalized_winner_m3
        << '\n';

    std::cout
        << "    M4="
        << stats.normalized_winner_m4
        << '\n';

    std::cout
        << "    M5="
        << stats.normalized_winner_m5
        << '\n';

    std::cout
        << "    M6="
        << stats.normalized_winner_m6
        << '\n';

    std::cout
        << "    M7="
        << stats.normalized_winner_m7
        << '\n';

    std::cout
        << "    M8="
        << stats.normalized_winner_m8
        << '\n';

    std::cout
        << "  MISMATCH_MATRIX_EXACT_TO_NORMALIZED:\n";

    for (u64 a = 1;
         a <= 8;
         ++a) {

        bool row_nonzero = false;

        for (u64 b = 1;
             b <= 8;
             ++b) {

            if (
                stats.mismatch_exact_m[a][b]
                != 0
            ) {
                row_nonzero = true;
            }
        }

        if (!row_nonzero) {
            continue;
        }

        std::cout
            << "    EXACT_M="
            << a;

        for (u64 b = 1;
             b <= 8;
             ++b) {

            if (
                stats.mismatch_exact_m[a][b]
                != 0
            ) {
                std::cout
                    << " NORMALIZED_M"
                    << b
                    << "="
                    << stats.mismatch_exact_m[a][b];
            }
        }

        std::cout << '\n';
    }
}

int main() {
    constexpr u64 EXPERIMENT = 404;
    constexpr u64 CASES = 3000;
    constexpr u64 K_LIMIT = 80;
    constexpr u64 M_MAX = 8;

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
        << "PRIME_MIN=10000\n";

    std::cout
        << "PRIME_MAX=1000000\n\n";

    std::mt19937_64 rng(
        404123456789ULL
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

    print_stats(p_stats);

    std::cout
        << "\nQ_STATS\n";

    print_stats(q_stats);

    const bool pass =
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
