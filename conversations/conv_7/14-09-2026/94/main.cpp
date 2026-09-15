#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
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
        u64 p = random_prime(rng, 5, 3000);
        u64 q = random_prime(rng, 3001, 10000);

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
    F(t)=t(N-t(s-1))
*/
static mpz_class F(
    const CaseData& c,
    u64 t
) {
    const mpz_class T = mpz_from_u64(t);
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);

    return T * (N - T * A);
}

/*
    C(t)=F(t+1)F(t-1)-F(t)^2
*/
static mpz_class C(
    const CaseData& c,
    u64 t
) {
    const mpz_class fm = F(c, t - 1);
    const mpz_class f0 = F(c, t);
    const mpz_class fp = F(c, t + 1);

    return fp * fm - f0 * f0;
}

static bool boolean_hit(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class value = C(c, t);

    const mpz_class g = gcd(N, value);

    return g > 1 && g < N;
}

/*
    Find all roots x in [1,r-1] satisfying

        2x^2 == 1 (mod r).

    r <= 10000 here, so a direct scan is completely adequate.
*/
static std::vector<u64> roots_mod_r(u64 r) {
    std::vector<u64> roots;

    for (u64 x = 1; x < r; ++x) {
        if (((2ULL * x * x) % r) == 1) {
            roots.push_back(x);
        }
    }

    return roots;
}

static bool root_exists_below_s(
    u64 r,
    u64 s,
    std::vector<u64>* roots_out = nullptr
) {
    const std::vector<u64> roots =
        roots_mod_r(r);

    if (roots_out != nullptr) {
        *roots_out = roots;
    }

    for (u64 root : roots) {
        if (root < s) {
            return true;
        }
    }

    return false;
}

static bool predicts_factor_hit_in_scan(
    u64 r,
    u64 s
) {
    if ((s - 1) % r == 0) {
        return true;
    }

    return root_exists_below_s(r, s);
}

int main() {
    constexpr int EXPERIMENT = 383;
    constexpr int CASES = 2000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38320260914ULL);

    u64 total_scan_points = 0;

    u64 cases_hit = 0;
    u64 cases_no_hit = 0;

    u64 cases_p_hit_predicted = 0;
    u64 cases_q_hit_predicted = 0;
    u64 cases_any_hit_predicted = 0;

    u64 prediction_failures = 0;

    u64 q_residue_but_no_root_below_s = 0;
    u64 q_residue_with_root_below_s = 0;

    u64 q_root1_below_s = 0;
    u64 q_root2_below_s = 0;

    u64 cases_p_exception = 0;
    u64 cases_p_exception_hit = 0;

    u64 total_first_hit_t = 0;
    u64 min_first_hit_t = UINT64_MAX;
    u64 max_first_hit_t = 0;

    u64 first_hit_p = 0;
    u64 first_hit_q = 0;
    u64 first_hit_one = 0;
    u64 first_hit_other = 0;

    /*
        Extra diagnostic:
        for q-residue cases, compare the smaller root to s.
    */
    u64 q_smaller_root_below_s = 0;
    u64 q_smaller_root_above_s = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        const bool p_exception =
            ((s - 1) % p == 0);

        if (p_exception) {
            ++cases_p_exception;
        }

        /*
            Exact finite-interval prediction for p.
        */
        const bool p_predicted =
            predicts_factor_hit_in_scan(
                p,
                s
            );

        /*
            Exact finite-interval prediction for q.
        */
        const bool q_predicted =
            predicts_factor_hit_in_scan(
                q,
                s
            );

        if (p_predicted) {
            ++cases_p_hit_predicted;
        }

        if (q_predicted) {
            ++cases_q_hit_predicted;
        }

        if (p_predicted || q_predicted) {
            ++cases_any_hit_predicted;
        }

        /*
            Examine q roots explicitly.
        */
        std::vector<u64> q_roots;

        const bool q_has_scan_root =
            root_exists_below_s(
                q,
                s,
                &q_roots
            );

        if (q_roots.size() == 2) {
            /*
                q>s, so both roots are distinct.
            */
            if (q_roots[0] < s) {
                ++q_root1_below_s;
            }

            if (q_roots[1] < s) {
                ++q_root2_below_s;
            }

            const u64 smaller =
                std::min(
                    q_roots[0],
                    q_roots[1]
                );

            if (smaller < s) {
                ++q_smaller_root_below_s;
            } else {
                ++q_smaller_root_above_s;
            }
        }

        /*
            Detect q being a quadratic residue without requiring
            a root in the scanned interval.
        */
        if (q % 8 == 1 ||
            q % 8 == 7) {

            if (q_has_scan_root) {
                ++q_residue_with_root_below_s;
            } else {
                ++q_residue_but_no_root_below_s;
            }
        }

        /*
            Scan the actual interval.
        */
        bool found = false;
        u64 first_t = 0;
        u64 first_factor = 0;

        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;

            if (!boolean_hit(c, t)) {
                continue;
            }

            found = true;
            first_t = t;

            const mpz_class N =
                mpz_from_u64(c.N);

            const mpz_class T =
                mpz_from_u64(t);

            const mpz_class poly =
                2 * T * T - 1;

            const mpz_class g =
                gcd(N, poly);

            first_factor =
                g.get_ui();

            break;
        }

        if (found) {
            ++cases_hit;

            total_first_hit_t += first_t;

            min_first_hit_t =
                std::min(
                    min_first_hit_t,
                    first_t
                );

            max_first_hit_t =
                std::max(
                    max_first_hit_t,
                    first_t
                );

            if (first_factor == p) {
                ++first_hit_p;
            } else if (first_factor == q) {
                ++first_hit_q;
            } else if (first_factor == 1) {
                ++first_hit_one;
            } else {
                ++first_hit_other;
            }

            if (p_exception) {
                ++cases_p_exception_hit;
            }
        } else {
            ++cases_no_hit;
        }

        /*
            Compare actual and finite-interval prediction.
        */
        if (found !=
            (p_predicted || q_predicted)) {

            ++prediction_failures;
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "TOTAL_SCAN_POINTS="
        << total_scan_points
        << '\n';

    std::cout
        << "CASES_HIT="
        << cases_hit
        << '\n';

    std::cout
        << "CASES_NO_HIT="
        << cases_no_hit
        << '\n';

    std::cout
        << "CASES_P_HIT_PREDICTED="
        << cases_p_hit_predicted
        << '\n';

    std::cout
        << "CASES_Q_HIT_PREDICTED="
        << cases_q_hit_predicted
        << '\n';

    std::cout
        << "CASES_ANY_HIT_PREDICTED="
        << cases_any_hit_predicted
        << '\n';

    std::cout
        << "PREDICTION_FAILURES="
        << prediction_failures
        << '\n';

    std::cout
        << "CASES_Q_RESIDUE_BUT_NO_ROOT_BELOW_S="
        << q_residue_but_no_root_below_s
        << '\n';

    std::cout
        << "CASES_Q_RESIDUE_WITH_ROOT_BELOW_S="
        << q_residue_with_root_below_s
        << '\n';

    std::cout
        << "Q_ROOT1_BELOW_S="
        << q_root1_below_s
        << '\n';

    std::cout
        << "Q_ROOT2_BELOW_S="
        << q_root2_below_s
        << '\n';

    std::cout
        << "Q_SMALLER_ROOT_BELOW_S="
        << q_smaller_root_below_s
        << '\n';

    std::cout
        << "Q_SMALLER_ROOT_ABOVE_S="
        << q_smaller_root_above_s
        << '\n';

    std::cout
        << "CASES_P_EXCEPTION="
        << cases_p_exception
        << '\n';

    std::cout
        << "CASES_P_EXCEPTION_HIT="
        << cases_p_exception_hit
        << '\n';

    std::cout
        << "FIRST_HIT_P="
        << first_hit_p
        << '\n';

    std::cout
        << "FIRST_HIT_Q="
        << first_hit_q
        << '\n';

    std::cout
        << "FIRST_HIT_ONE="
        << first_hit_one
        << '\n';

    std::cout
        << "FIRST_HIT_OTHER="
        << first_hit_other
        << '\n';

    std::cout
        << "TOTAL_FIRST_HIT_T="
        << total_first_hit_t
        << '\n';

    std::cout
        << "MIN_FIRST_HIT_T="
        << min_first_hit_t
        << '\n';

    std::cout
        << "MAX_FIRST_HIT_T="
        << max_first_hit_t
        << '\n';

    std::cout
        << "AVERAGE_FIRST_HIT_T="
        << (
            cases_hit > 0
                ? static_cast<double>(
                    total_first_hit_t
                  ) /
                  static_cast<double>(
                    cases_hit
                  )
                : 0.0
        )
        << '\n';

    const bool status =
        prediction_failures == 0;

    std::cout
        << "FINITE_INTERVAL_CLASSIFICATION_STATUS="
        << (
            status
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
