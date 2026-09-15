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

struct MismatchStats {
    u64 total = 0;

    u64 delta_one = 0;
    u64 delta_gt_one = 0;

    u64 sign_pp = 0;
    u64 sign_pm = 0;
    u64 sign_mp = 0;
    u64 sign_mm = 0;

    u64 correction_positive = 0;
    u64 correction_zero = 0;
    u64 correction_negative = 0;

    u64 exact_condition_true = 0;
    u64 exact_condition_false = 0;

    u64 max_required_r = 0;

    u64 max_k1 = 0;
    u64 max_k2 = 0;
    u64 max_delta = 0;

    u64 r_below_threshold = 0;
    u64 r_at_or_above_threshold = 0;
};

struct GlobalStats {
    u64 primes = 0;

    u64 pair_tests = 0;

    u64 normalized_first = 0;
    u64 normalized_second = 0;

    u64 exact_agree = 0;
    u64 exact_disagree = 0;

    u64 correction_prediction_failures = 0;

    MismatchStats mismatch;

    u64 counterexample_r =
        std::numeric_limits<u64>::max();

    u64 counterexample_k = 0;

    Witness counterexample_w1;
    Witness counterexample_w2;
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

        /*
         * Largest divisor k <= K gives
         * the smallest possible t.
         */
        for (u64 k = K;
             k >= 2;
             --k) {

            if (value_u64 % k != 0) {
                continue;
            }

            const Witness candidate =
                make_candidate(
                    r,
                    m,
                    k,
                    sign
                );

            if (!candidate.exists) {
                continue;
            }

            if (
                !best.exists ||
                candidate.t < best.t ||
                (
                    candidate.t == best.t &&
                    (
                        candidate.k < best.k ||
                        (
                            candidate.k == best.k &&
                            candidate.sign < best.sign
                        )
                    )
                )
            ) {
                best = candidate;
            }

            break;
        }
    }

    return best;
}

static u128 numerator(
    u64 r,
    const Witness& w
) {
    const u128 mr =
        static_cast<u128>(w.m) * r;

    if (w.sign == +1) {
        return mr + 1;
    }

    return mr - 1;
}

static bool normalized_first(
    const Witness& a,
    const Witness& b
) {
    /*
     * k1/m1 > k2/m2
     *
     * iff
     *
     * k1*m2 > k2*m1.
     */
    const u128 lhs =
        static_cast<u128>(a.k) * b.m;

    const u128 rhs =
        static_cast<u128>(b.k) * a.m;

    return lhs > rhs;
}

static bool exact_first(
    u64 r,
    const Witness& a,
    const Witness& b
) {
    /*
     * (m1*r+e1)/k1 <
     * (m2*r+e2)/k2
     */
    const u128 lhs =
        numerator(r, a) * b.k;

    const u128 rhs =
        numerator(r, b) * a.k;

    return lhs < rhs;
}

/*
 * Called only when normalized_first(a,b) == true.
 *
 * Define
 *
 *   Delta = m2*k1 - m1*k2 > 0
 *
 * and
 *
 *   C = e1*k2 - e2*k1.
 *
 * Exact a < b iff
 *
 *   Delta*r > C.
 */
static void classify_mismatch(
    u64 r,
    const Witness& a,
    const Witness& b,
    GlobalStats& stats
) {
    const u128 delta128 =
        static_cast<u128>(b.m) * a.k -
        static_cast<u128>(a.m) * b.k;

    const u64 delta =
        static_cast<u64>(delta128);

    ++stats.mismatch.total;

    if (delta == 1) {
        ++stats.mismatch.delta_one;
    } else {
        ++stats.mismatch.delta_gt_one;
    }

    stats.mismatch.max_delta =
        std::max(
            stats.mismatch.max_delta,
            delta
        );

    stats.mismatch.max_k1 =
        std::max(
            stats.mismatch.max_k1,
            a.k
        );

    stats.mismatch.max_k2 =
        std::max(
            stats.mismatch.max_k2,
            b.k
        );

    if (a.sign == +1 && b.sign == +1) {
        ++stats.mismatch.sign_pp;
    } else if (a.sign == +1 && b.sign == -1) {
        ++stats.mismatch.sign_pm;
    } else if (a.sign == -1 && b.sign == +1) {
        ++stats.mismatch.sign_mp;
    } else {
        ++stats.mismatch.sign_mm;
    }

    /*
     * C = e1*k2 - e2*k1.
     *
     * Compute positive and negative contributions
     * without signed arithmetic.
     */
    u64 positive_part = 0;
    u64 negative_part = 0;

    if (a.sign == +1) {
        positive_part += b.k;
    } else {
        negative_part += b.k;
    }

    if (b.sign == +1) {
        negative_part += a.k;
    } else {
        positive_part += a.k;
    }

    int correction_sign = 0;
    u64 correction_abs = 0;

    if (positive_part > negative_part) {
        correction_sign = +1;
        correction_abs =
            positive_part - negative_part;
    } else if (negative_part > positive_part) {
        correction_sign = -1;
        correction_abs =
            negative_part - positive_part;
    } else {
        correction_sign = 0;
        correction_abs = 0;
    }

    if (correction_sign > 0) {
        ++stats.mismatch.correction_positive;
    } else if (correction_sign < 0) {
        ++stats.mismatch.correction_negative;
    } else {
        ++stats.mismatch.correction_zero;
    }

    /*
     * Exact condition:
     *
     *   Delta*r > C.
     *
     * If C <= 0, this is automatically true.
     */
    bool predicted_exact_first = false;

    if (correction_sign <= 0) {
        predicted_exact_first = true;
    } else {
        predicted_exact_first =
            delta128 * static_cast<u128>(r)
            >
            static_cast<u128>(correction_abs);
    }

    if (predicted_exact_first) {
        ++stats.mismatch.exact_condition_true;
        ++stats.correction_prediction_failures;
    } else {
        ++stats.mismatch.exact_condition_false;
    }

    if (correction_sign > 0) {
        /*
         * Exact threshold:
         *
         *   r >= floor(C/Delta) + 1
         */
        const u64 threshold =
            correction_abs / delta + 1;

        stats.mismatch.max_required_r =
            std::max(
                stats.mismatch.max_required_r,
                threshold
            );

        if (r < threshold) {
            ++stats.mismatch.r_below_threshold;
        } else {
            ++stats.mismatch.r_at_or_above_threshold;
        }
    }
}

static void analyze_prime(
    u64 r,
    u64 K_LIMIT,
    u64 M_MAX,
    GlobalStats& stats
) {
    ++stats.primes;

    for (u64 K = 2;
         K <= K_LIMIT;
         ++K) {

        std::vector<Witness> witnesses(
            M_MAX + 1
        );

        for (u64 m = 1;
             m <= M_MAX;
             ++m) {

            witnesses[m] =
                best_for_m(
                    r,
                    K,
                    m
                );
        }

        for (u64 m1 = 1;
             m1 <= M_MAX;
             ++m1) {

            if (!witnesses[m1].exists) {
                continue;
            }

            for (u64 m2 = m1 + 1;
                 m2 <= M_MAX;
                 ++m2) {

                if (!witnesses[m2].exists) {
                    continue;
                }

                const Witness& a =
                    witnesses[m1];

                const Witness& b =
                    witnesses[m2];

                ++stats.pair_tests;

                const bool nf =
                    normalized_first(
                        a,
                        b
                    );

                const bool ef =
                    exact_first(
                        r,
                        a,
                        b
                    );

                if (nf) {
                    ++stats.normalized_first;
                } else {
                    ++stats.normalized_second;
                }

                if (nf == ef) {
                    ++stats.exact_agree;
                    continue;
                }

                ++stats.exact_disagree;

                classify_mismatch(
                    r,
                    a,
                    b,
                    stats
                );

                if (
                    r <
                    stats.counterexample_r
                ) {
                    stats.counterexample_r = r;
                    stats.counterexample_k = K;
                    stats.counterexample_w1 = a;
                    stats.counterexample_w2 = b;
                }
            }
        }
    }
}

static void print_witness(
    const char* label,
    const Witness& w
) {
    std::cout
        << "  "
        << label
        << "_M="
        << w.m
        << '\n';

    std::cout
        << "  "
        << label
        << "_K="
        << w.k
        << '\n';

    std::cout
        << "  "
        << label
        << "_T="
        << w.t
        << '\n';

    std::cout
        << "  "
        << label
        << "_SIGN="
        << w.sign
        << '\n';
}

int main() {
    constexpr u64 EXPERIMENT = 408;
    constexpr u64 PRIME_LIMIT = 5000;
    constexpr u64 K_LIMIT = 500;
    constexpr u64 M_MAX = 64;

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
        << "EXACT_CORRECTION_CLASSIFICATION=1\n\n";

    GlobalStats stats;

    for (u64 r = 2;
         r <= PRIME_LIMIT;
         ++r) {

        if (!is_prime(r)) {
            continue;
        }

        analyze_prime(
            r,
            K_LIMIT,
            M_MAX,
            stats
        );
    }

    std::cout
        << "\nGLOBAL_STATS\n";

    std::cout
        << "  PRIMES="
        << stats.primes
        << '\n';

    std::cout
        << "  PAIR_TESTS="
        << stats.pair_tests
        << '\n';

    std::cout
        << "  NORMALIZED_FIRST="
        << stats.normalized_first
        << '\n';

    std::cout
        << "  NORMALIZED_SECOND="
        << stats.normalized_second
        << '\n';

    std::cout
        << "  EXACT_AGREE="
        << stats.exact_agree
        << '\n';

    std::cout
        << "  EXACT_DISAGREE="
        << stats.exact_disagree
        << '\n';

    std::cout
        << "  MISMATCH_DELTA_1="
        << stats.mismatch.delta_one
        << '\n';

    std::cout
        << "  MISMATCH_DELTA_GT1="
        << stats.mismatch.delta_gt_one
        << '\n';

    std::cout
        << "  SIGN_PP="
        << stats.mismatch.sign_pp
        << '\n';

    std::cout
        << "  SIGN_PM="
        << stats.mismatch.sign_pm
        << '\n';

    std::cout
        << "  SIGN_MP="
        << stats.mismatch.sign_mp
        << '\n';

    std::cout
        << "  SIGN_MM="
        << stats.mismatch.sign_mm
        << '\n';

    std::cout
        << "  CORRECTION_POSITIVE="
        << stats.mismatch.correction_positive
        << '\n';

    std::cout
        << "  CORRECTION_ZERO="
        << stats.mismatch.correction_zero
        << '\n';

    std::cout
        << "  CORRECTION_NEGATIVE="
        << stats.mismatch.correction_negative
        << '\n';

    std::cout
        << "  MAX_DELTA="
        << stats.mismatch.max_delta
        << '\n';

    std::cout
        << "  MAX_K1="
        << stats.mismatch.max_k1
        << '\n';

    std::cout
        << "  MAX_K2="
        << stats.mismatch.max_k2
        << '\n';

    std::cout
        << "  EXACT_CORRECTION_TRUE="
        << stats.mismatch.exact_condition_true
        << '\n';

    std::cout
        << "  EXACT_CORRECTION_FALSE="
        << stats.mismatch.exact_condition_false
        << '\n';

    std::cout
        << "  CORRECTION_PREDICTION_FAILURES="
        << stats.correction_prediction_failures
        << '\n';

    std::cout
        << "  MAX_REQUIRED_R="
        << stats.mismatch.max_required_r
        << '\n';

    std::cout
        << "  MISMATCHES_R_BELOW_THRESHOLD="
        << stats.mismatch.r_below_threshold
        << '\n';

    std::cout
        << "  MISMATCHES_R_AT_OR_ABOVE_THRESHOLD="
        << stats.mismatch.r_at_or_above_threshold
        << '\n';

    if (
        stats.counterexample_r !=
        std::numeric_limits<u64>::max()
    ) {
        std::cout
            << "\nFIRST_MISMATCH\n";

        std::cout
            << "  R="
            << stats.counterexample_r
            << '\n';

        std::cout
            << "  K="
            << stats.counterexample_k
            << '\n';

        print_witness(
            "A",
            stats.counterexample_w1
        );

        print_witness(
            "B",
            stats.counterexample_w2
        );
    } else {
        std::cout
            << "\nFIRST_MISMATCH=NONE\n";
    }

    const bool pass =
        stats.correction_prediction_failures == 0;

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