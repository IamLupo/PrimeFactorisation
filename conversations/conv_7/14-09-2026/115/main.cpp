#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Witness {
    bool exists = false;
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0; // -1 => kt = mr-1, +1 => kt = mr+1
};

struct PrimeStats {
    u64 prime = 0;

    u64 points = 0;
    u64 matches = 0;
    u64 mismatches = 0;

    u64 pair_tests = 0;

    u64 sufficient_condition_tests = 0;
    u64 sufficient_condition_failures = 0;

    bool first_mismatch_found = false;
    u64 mismatch_K = 0;

    Witness mismatch_exact;
    u64 mismatch_normalized_m = 0;
};

struct GlobalStats {
    u64 primes = 0;
    u64 prime_limit = 0;

    u64 total_points = 0;
    u64 exact_matches = 0;
    u64 exact_mismatches = 0;

    u64 pair_tests = 0;
    u64 sufficient_condition_tests = 0;
    u64 sufficient_condition_failures = 0;

    u64 first_counterexample_r = 0;
    u64 counterexample_K = 0;

    Witness counterexample_exact;
    u64 counterexample_normalized_m = 0;
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
        30
    ) > 0;
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

static Witness make_candidate(
    u64 r,
    u64 m,
    u64 k,
    int sign
) {
    Witness w;

    const u128 mr =
        static_cast<u128>(m) * r;

    u128 value;

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
            sign_index == 0
                ? -1
                : +1;

        u128 value;

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
                best =
                    make_candidate(
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

static u64 divisor_score(
    u64 r,
    u64 K,
    u64 m
) {
    const u128 mr =
        static_cast<u128>(m) * r;

    u64 best_d = 0;

    /*
     * mr - 1
     */
    if (mr > 1) {
        const u64 value =
            static_cast<u64>(mr - 1);

        for (u64 k = K;
             k >= 2;
             --k) {

            if (value % k == 0) {
                best_d = std::max(best_d, k);
                break;
            }
        }
    }

    /*
     * mr + 1
     */
    {
        const u64 value =
            static_cast<u64>(mr + 1);

        for (u64 k = K;
             k >= 2;
             --k) {

            if (value % k == 0) {
                best_d = std::max(best_d, k);
                break;
            }
        }
    }

    return best_d;
}

static u64 normalized_winner(
    u64 r,
    u64 K,
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
            divisor_score(
                r,
                K,
                m
            );

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

/*
 * Compare exact values:

 *   A = (m1*r + sign1) / k1
 *   B = (m2*r + sign2) / k2
 *
 * No signed arithmetic is needed.
 */
static bool exact_prefers_first(
    u64 r,
    const Witness& a,
    const Witness& b
) {
    u128 value_a;

    if (a.sign == +1) {
        value_a =
            static_cast<u128>(a.m) * r + 1;
    } else {
        value_a =
            static_cast<u128>(a.m) * r - 1;
    }

    u128 value_b;

    if (b.sign == +1) {
        value_b =
            static_cast<u128>(b.m) * r + 1;
    } else {
        value_b =
            static_cast<u128>(b.m) * r - 1;
    }

    const u128 lhs =
        value_a * b.k;

    const u128 rhs =
        value_b * a.k;

    return lhs < rhs;
}

static bool normalized_prefers_first(
    const Witness& a,
    const Witness& b
) {
    const u128 lhs =
        static_cast<u128>(a.k) * b.m;

    const u128 rhs =
        static_cast<u128>(b.k) * a.m;

    return lhs > rhs;
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    u64 m_max,
    GlobalStats& global
) {
    PrimeStats local;
    local.prime = r;

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        const Witness exact =
            exact_winner(
                r,
                K,
                m_max
            );

        if (!exact.exists) {
            continue;
        }

        ++local.points;
        ++global.total_points;

        if (!verify_witness(r, exact)) {
            ++local.sufficient_condition_failures;
            ++global.sufficient_condition_failures;
        }

        bool normalized_unique = true;

        const u64 normalized_m =
            normalized_winner(
                r,
                K,
                m_max,
                normalized_unique
            );

        if (normalized_m == 0) {
            continue;
        }

        if (
            normalized_unique &&
            normalized_m == exact.m
        ) {
            ++local.matches;
            ++global.exact_matches;
        } else {
            ++local.mismatches;
            ++global.exact_mismatches;

            if (!local.first_mismatch_found) {
                local.first_mismatch_found = true;
                local.mismatch_K = K;
                local.mismatch_exact = exact;
                local.mismatch_normalized_m =
                    normalized_m;
            }

            if (global.first_counterexample_r == 0) {
                global.first_counterexample_r = r;
                global.counterexample_K = K;
                global.counterexample_exact = exact;
                global.counterexample_normalized_m =
                    normalized_m;
            }
        }

        /*
         * Pairwise test.
         */
        for (u64 m1 = 1;
             m1 <= m_max;
             ++m1) {

            const Witness w1 =
                best_for_m(
                    r,
                    K,
                    m1
                );

            if (!w1.exists) {
                continue;
            }

            for (u64 m2 = m1 + 1;
                 m2 <= m_max;
                 ++m2) {

                const Witness w2 =
                    best_for_m(
                        r,
                        K,
                        m2
                    );

                if (!w2.exists) {
                    continue;
                }

                ++local.pair_tests;
                ++global.pair_tests;

                const bool normalized_first =
                    normalized_prefers_first(
                        w1,
                        w2
                    );

                const bool exact_first =
                    exact_prefers_first(
                        r,
                        w1,
                        w2
                    );

                if (!normalized_first) {
                    continue;
                }

                /*
                 * Sufficient condition:
                 *
                 * r > k1 + k2
                 *
                 * makes the r-term dominate the
                 * +/-1 correction.
                 */
                if (
                    r >
                    static_cast<u64>(w1.k + w2.k)
                ) {
                    ++local.sufficient_condition_tests;
                    ++global.sufficient_condition_tests;

                    if (!exact_first) {
                        ++local.sufficient_condition_failures;
                        ++global.sufficient_condition_failures;
                    }
                }
            }
        }
    }

    if (local.first_mismatch_found) {
        std::cout
            << "PRIME_COUNTEREXAMPLE"
            << " r="
            << r
            << " K="
            << local.mismatch_K
            << " EXACT_M="
            << local.mismatch_exact.m
            << " EXACT_K="
            << local.mismatch_exact.k
            << " EXACT_T="
            << local.mismatch_exact.t
            << " NORMALIZED_M="
            << local.mismatch_normalized_m
            << '\n';
    }
}

int main() {
    constexpr u64 EXPERIMENT = 406;
    constexpr u64 PRIME_LIMIT = 5000;
    constexpr u64 K_LIMIT = 200;
    constexpr u64 M_MAX = 32;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT
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
        << "GMP_PRIMALITY=1\n\n";

    GlobalStats global;
    global.prime_limit = PRIME_LIMIT;

    for (u64 r = 2;
         r <= PRIME_LIMIT;
         ++r) {

        if (!is_prime(r)) {
            continue;
        }

        ++global.primes;

        analyze_prime(
            r,
            K_LIMIT,
            M_MAX,
            global
        );
    }

    std::cout
        << "\nGLOBAL_STATS\n";

    std::cout
        << "  PRIMES="
        << global.primes
        << '\n';

    std::cout
        << "  TOTAL_POINTS="
        << global.total_points
        << '\n';

    std::cout
        << "  EXACT_MATCHES="
        << global.exact_matches
        << '\n';

    std::cout
        << "  EXACT_MISMATCHES="
        << global.exact_mismatches
        << '\n';

    std::cout
        << "  PAIR_TESTS="
        << global.pair_tests
        << '\n';

    std::cout
        << "  SUFFICIENT_CONDITION_TESTS="
        << global.sufficient_condition_tests
        << '\n';

    std::cout
        << "  SUFFICIENT_CONDITION_FAILURES="
        << global.sufficient_condition_failures
        << '\n';

    if (global.first_counterexample_r != 0) {
        std::cout
            << "  FIRST_COUNTEREXAMPLE_R="
            << global.first_counterexample_r
            << '\n';

        std::cout
            << "  COUNTEREXAMPLE_K="
            << global.counterexample_K
            << '\n';

        std::cout
            << "  COUNTEREXAMPLE_EXACT_M="
            << global.counterexample_exact.m
            << '\n';

        std::cout
            << "  COUNTEREXAMPLE_EXACT_K="
            << global.counterexample_exact.k
            << '\n';

        std::cout
            << "  COUNTEREXAMPLE_EXACT_T="
            << global.counterexample_exact.t
            << '\n';

        std::cout
            << "  COUNTEREXAMPLE_NORMALIZED_M="
            << global.counterexample_normalized_m
            << '\n';
    } else {
        std::cout
            << "  FIRST_COUNTEREXAMPLE=NONE\n";
    }

    const bool pass =
        global.exact_mismatches == 0 &&
        global.sufficient_condition_failures == 0;

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