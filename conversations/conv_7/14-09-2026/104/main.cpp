#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

struct RootData {
    u64 root = 0;
    bool exists = false;
};

static mpz_class mpz_from_u64(u64 value) {
    return mpz_class(std::to_string(value));
}

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        std::sqrt(static_cast<long double>(n))
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static bool is_prime_u64(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d <= n / d; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

static u64 random_prime(
    std::mt19937_64& rng,
    u64 lo,
    u64 hi
) {
    std::uniform_int_distribution<u64> dist(lo, hi);

    while (true) {
        u64 x = dist(rng);

        if (x < 2) {
            x = 2;
        }

        if (x > 2 && x % 2 == 0) {
            ++x;
        }

        while (x <= hi && !is_prime_u64(x)) {
            x += 2;
        }

        if (x >= lo && x <= hi) {
            return x;
        }
    }
}

static CaseData make_case(
    std::mt19937_64& rng
) {
    while (true) {
        u64 p = random_prime(rng, 5, 5000);
        u64 q = random_prime(rng, 5001, 30000);

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        const u64 N = p * q;
        const u64 s = isqrt_u64(N);

        if (s * s > N) {
            continue;
        }

        if ((s + 1) * (s + 1) <= N) {
            continue;
        }

        return {p, q, N, s};
    }
}

/*
    Smallest t>0 satisfying

        k*t == +1 (mod r)

    or

        k*t == -1 (mod r).

    We use a direct scan here because the experiment's primes are
    small enough. This gives an unambiguous ground truth.
*/
static RootData inverse_root(
    u64 r,
    u64 k
) {
    const u64 km = k % r;

    if (km == 0) {
        return {};
    }

    const u64 limit = r;

    for (u64 t = 1; t < limit; ++t) {
        const u64 v =
            (km * (t % r)) % r;

        if (v == 1 ||
            v == r - 1) {
            return {t, true};
        }
    }

    return {};
}

static u64 gcd_candidate(
    const CaseData& c,
    u64 k,
    u64 t
) {
    const mpz_class N =
        mpz_from_u64(c.N);

    const mpz_class K =
        mpz_from_u64(k);

    const mpz_class T =
        mpz_from_u64(t);

    const mpz_class value =
        K * K * T * T - 1;

    const mpz_class g =
        gcd(N, value);

    return g.get_ui();
}

int main() {
    constexpr int EXPERIMENT = 393;
    constexpr int CASES = 3000;

    constexpr u64 K_MIN = 2;
    constexpr u64 K_MAX = 20;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x39320260914ULL);

    struct KStats {
        u64 cases = 0;

        u64 p_first = 0;
        u64 q_first = 0;
        u64 neither = 0;

        u64 inverse_prediction_failures = 0;
        u64 gcd_verification_failures = 0;

        u64 total_rho_p = 0;
        u64 total_rho_q = 0;

        u64 min_rho_p = UINT64_MAX;
        u64 max_rho_p = 0;

        u64 min_rho_q = UINT64_MAX;
        u64 max_rho_q = 0;

        u64 q_first_rho_cases = 0;

        u64 total_ratio_qp_num = 0;
        u64 total_ratio_qp_den = 0;
    };

    std::vector<KStats> stats(
        K_MAX + 1
    );

    u64 total_inverse_checks = 0;

    u64 cases_all_k_explained = 0;
    u64 cases_any_prediction_failure = 0;

    /*
        How often does q-first occur because its inverse root is
        genuinely smaller than the p inverse root?
    */
    u64 total_q_first_inverse_wins = 0;

    /*
        Track the smallest observed rho values.
    */
    u64 global_min_rho_p = UINT64_MAX;
    u64 global_max_rho_p = 0;

    u64 global_min_rho_q = UINT64_MAX;
    u64 global_max_rho_q = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c =
            make_case(rng);

        bool case_all_correct = true;

        for (u64 k = K_MIN;
             k <= K_MAX;
             ++k) {

            ++total_inverse_checks;

            KStats& st =
                stats[k];

            ++st.cases;

            const RootData rp =
                inverse_root(
                    c.p,
                    k
                );

            const RootData rq =
                inverse_root(
                    c.q,
                    k
                );

            const u64 limit =
                c.s / 2;

            const bool p_exists =
                rp.exists &&
                rp.root <= limit;

            const bool q_exists =
                rq.exists &&
                rq.root <= limit;

            if (rp.exists) {
                st.total_rho_p +=
                    rp.root;

                st.min_rho_p =
                    std::min(
                        st.min_rho_p,
                        rp.root
                    );

                st.max_rho_p =
                    std::max(
                        st.max_rho_p,
                        rp.root
                    );

                global_min_rho_p =
                    std::min(
                        global_min_rho_p,
                        rp.root
                    );

                global_max_rho_p =
                    std::max(
                        global_max_rho_p,
                        rp.root
                    );
            }

            if (rq.exists) {
                st.total_rho_q +=
                    rq.root;

                st.min_rho_q =
                    std::min(
                        st.min_rho_q,
                        rq.root
                    );

                st.max_rho_q =
                    std::max(
                        st.max_rho_q,
                        rq.root
                    );

                global_min_rho_q =
                    std::min(
                        global_min_rho_q,
                        rq.root
                    );

                global_max_rho_q =
                    std::max(
                        global_max_rho_q,
                        rq.root
                    );
            }

            /*
                Prediction from inverse roots.
            */
            if (!p_exists &&
                !q_exists) {

                ++st.neither;

            } else if (
                q_exists &&
                (
                    !p_exists ||
                    rq.root < rp.root
                )
            ) {

                ++st.q_first;

            } else {

                ++st.p_first;
            }

            /*
                Verify that the predicted winner really produces
                the expected factor using GMP.
            */
            if (q_exists &&
                (!p_exists ||
                 rq.root < rp.root)) {

                const u64 g =
                    gcd_candidate(
                        c,
                        k,
                        rq.root
                    );

                if (g != c.q &&
                    g != c.N) {

                    ++st.gcd_verification_failures;
                    case_all_correct = false;
                }

                ++st.q_first_rho_cases;
                ++total_q_first_inverse_wins;

            } else if (p_exists) {

                const u64 g =
                    gcd_candidate(
                        c,
                        k,
                        rp.root
                    );

                if (g != c.p &&
                    g != c.N) {

                    ++st.gcd_verification_failures;
                    case_all_correct = false;
                }
            }

            /*
                Compare the inverse-root prediction against an actual
                direct search over t <= s/2.
            */
            u64 actual_t = 0;
            u64 actual_factor = 0;

            for (u64 t = 1;
                 t <= limit;
                 ++t) {

                const u64 g =
                    gcd_candidate(
                        c,
                        k,
                        t
                    );

                if (g > 1 &&
                    g < c.N) {

                    actual_t = t;
                    actual_factor = g;
                    break;
                }
            }

            u64 predicted_t = 0;
            u64 predicted_factor = 0;

            if (q_exists &&
                (!p_exists ||
                 rq.root < rp.root)) {

                predicted_t = rq.root;
                predicted_factor = c.q;

            } else if (p_exists) {

                predicted_t = rp.root;
                predicted_factor = c.p;
            }

            if (actual_t != predicted_t) {
                ++st.inverse_prediction_failures;
                case_all_correct = false;
            }

            if (actual_t != 0 &&
                predicted_t != 0 &&
                actual_factor != predicted_factor &&
                actual_factor != c.N) {

                ++st.inverse_prediction_failures;
                case_all_correct = false;
            }
        }

        if (!case_all_correct) {
            ++cases_any_prediction_failure;
        } else {
            ++cases_all_k_explained;
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "K_MIN="
        << K_MIN
        << '\n';

    std::cout
        << "K_MAX="
        << K_MAX
        << '\n';

    std::cout
        << "TOTAL_INVERSE_CHECKS="
        << total_inverse_checks
        << '\n';

    std::cout
        << "CASES_ALL_K_EXPLAINED="
        << cases_all_k_explained
        << '\n';

    std::cout
        << "CASES_ANY_PREDICTION_FAILURE="
        << cases_any_prediction_failure
        << '\n';

    std::cout
        << "TOTAL_Q_FIRST_INVERSE_WINS="
        << total_q_first_inverse_wins
        << '\n';

    std::cout
        << "GLOBAL_MIN_RHO_P="
        << global_min_rho_p
        << '\n';

    std::cout
        << "GLOBAL_MAX_RHO_P="
        << global_max_rho_p
        << '\n';

    std::cout
        << "GLOBAL_MIN_RHO_Q="
        << global_min_rho_q
        << '\n';

    std::cout
        << "GLOBAL_MAX_RHO_Q="
        << global_max_rho_q
        << '\n';

    for (u64 k = K_MIN;
         k <= K_MAX;
         ++k) {

        const KStats& st =
            stats[k];

        std::cout
            << "K="
            << k
            << '\n';

        std::cout
            << "  P_FIRST="
            << st.p_first
            << '\n';

        std::cout
            << "  Q_FIRST="
            << st.q_first
            << '\n';

        std::cout
            << "  NEITHER="
            << st.neither
            << '\n';

        std::cout
            << "  Q_FIRST_INVERSE_WINS="
            << st.q_first_rho_cases
            << '\n';

        std::cout
            << "  TOTAL_RHO_P="
            << st.total_rho_p
            << '\n';

        std::cout
            << "  TOTAL_RHO_Q="
            << st.total_rho_q
            << '\n';

        std::cout
            << "  MIN_RHO_P="
            << st.min_rho_p
            << '\n';

        std::cout
            << "  MAX_RHO_P="
            << st.max_rho_p
            << '\n';

        std::cout
            << "  MIN_RHO_Q="
            << st.min_rho_q
            << '\n';

        std::cout
            << "  MAX_RHO_Q="
            << st.max_rho_q
            << '\n';

        std::cout
            << "  INVERSE_PREDICTION_FAILURES="
            << st.inverse_prediction_failures
            << '\n';

        std::cout
            << "  GCD_VERIFICATION_FAILURES="
            << st.gcd_verification_failures
            << '\n';

        std::cout
            << "  AVERAGE_RHO_P="
            << (
                st.cases > 0
                    ? static_cast<double>(
                        st.total_rho_p
                      ) /
                      static_cast<double>(
                        st.cases
                      )
                    : 0.0
            )
            << '\n';

        std::cout
            << "  AVERAGE_RHO_Q="
            << (
                st.cases > 0
                    ? static_cast<double>(
                        st.total_rho_q
                      ) /
                      static_cast<double>(
                        st.cases
                      )
                    : 0.0
            )
            << '\n';
    }

    std::cout
        << "STATUS="
        << (
            cases_any_prediction_failure == 0
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
