#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

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
    F(t) = t(N - t(s-1))
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
    C(t) = F(t+1)F(t-1) - F(t)^2
*/
static mpz_class cross_difference(
    const CaseData& c,
    u64 t
) {
    const mpz_class fp = F(c, t + 1);
    const mpz_class fm = F(c, t - 1);
    const mpz_class f0 = F(c, t);

    return fp * fm - f0 * f0;
}

/*
    Theoretical formula:

        gcd(N,C(t))
        =
        gcd(N,(s-1)^2(1-t^2))
*/
static u64 theoretical_gcd(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        A * A * (1 - T * T);

    const mpz_class g = gcd(N, value);

    return g.get_ui();
}

static bool residue_pm_one(
    u64 t,
    u64 r
) {
    return
        (t % r == 1) ||
        (t % r == r - 1);
}

int main() {
    constexpr int EXPERIMENT = 369;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x36920260914ULL);

    u64 total_points = 0;

    u64 algebra_failures = 0;
    u64 gcd_formula_failures = 0;

    u64 p_classification_failures = 0;
    u64 q_classification_failures = 0;

    u64 cases_p_divides_s_minus_1 = 0;
    u64 cases_q_divides_s_minus_1 = 0;

    u64 total_p_hits = 0;
    u64 total_q_hits = 0;

    u64 total_p_hits_predicted = 0;
    u64 total_q_hits_predicted = 0;

    u64 cases_p_hit = 0;
    u64 cases_q_hit = 0;

    u64 cases_p_only = 0;
    u64 cases_q_only = 0;
    u64 cases_both = 0;

    u64 cases_other_gcd = 0;
    u64 total_other_gcds = 0;

    u64 max_gcd = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const bool p_divides_A =
            ((s - 1) % p == 0);

        const bool q_divides_A =
            ((s - 1) % q == 0);

        if (p_divides_A) {
            ++cases_p_divides_s_minus_1;
        }

        if (q_divides_A) {
            ++cases_q_divides_s_minus_1;
        }

        bool saw_p = false;
        bool saw_q = false;

        for (u64 t = 1; t < s; ++t) {
            ++total_points;

            /*
                Closed form:

                    C(t)
                    =
                    (s-1)^2(1-t^2)
                    -
                    (N-(s-1)t)^2
            */
            {
                const mpz_class T = mpz_from_u64(t);
                const mpz_class NN = mpz_from_u64(N);
                const mpz_class A = mpz_from_u64(s - 1);

                const mpz_class direct =
                    cross_difference(c, t);

                const mpz_class closed =
                    A * A * (1 - T * T)
                    -
                    (NN - A * T) *
                    (NN - A * T);

                if (direct != closed) {
                    ++algebra_failures;
                }
            }

            /*
                Actual gcd.
            */
            const mpz_class N_mpz =
                mpz_from_u64(N);

            const mpz_class C_mpz =
                cross_difference(c, t);

            const mpz_class g_mpz =
                gcd(N_mpz, C_mpz);

            const u64 actual =
                g_mpz.get_ui();

            /*
                Theoretical gcd.
            */
            const u64 theory =
                theoretical_gcd(c, t);

            if (actual != theory) {
                ++gcd_formula_failures;
            }

            max_gcd =
                std::max(max_gcd, actual);

            /*
                p prediction.
            */
            bool p_predicted;

            if (p_divides_A) {
                p_predicted = true;
            } else {
                p_predicted =
                    residue_pm_one(t, p);
            }

            /*
                q prediction.
            */
            bool q_predicted;

            if (q_divides_A) {
                q_predicted = true;
            } else {
                q_predicted =
                    residue_pm_one(t, q);
            }

            const bool p_actual =
                (actual % p == 0);

            const bool q_actual =
                (actual % q == 0);

            if (p_actual != p_predicted) {
                ++p_classification_failures;
            }

            if (q_actual != q_predicted) {
                ++q_classification_failures;
            }

            if (p_predicted) {
                ++total_p_hits_predicted;
            }

            if (q_predicted) {
                ++total_q_hits_predicted;
            }

            if (p_actual) {
                ++total_p_hits;
                saw_p = true;
            }

            if (q_actual) {
                ++total_q_hits;
                saw_q = true;
            }

            /*
                Any proper divisor of N other than p or q
                would indicate an unexpected structure.
            */
            if (actual > 1 && actual < N) {
                if (actual != p && actual != q) {
                    ++total_other_gcds;
                }
            }
        }

        if (saw_p) {
            ++cases_p_hit;
        }

        if (saw_q) {
            ++cases_q_hit;
        }

        /*
            Determine whether this case ever produced an
            unexpected proper gcd.
        */
        {
            bool other = false;

            for (u64 t = 1; t < s; ++t) {
                const mpz_class N_mpz =
                    mpz_from_u64(N);

                const mpz_class C_mpz =
                    cross_difference(c, t);

                const mpz_class g_mpz =
                    gcd(N_mpz, C_mpz);

                const u64 g =
                    g_mpz.get_ui();

                if (g > 1 &&
                    g < N &&
                    g != p &&
                    g != q) {
                    other = true;
                    break;
                }
            }

            if (other) {
                ++cases_other_gcd;
            }
        }

        if (saw_p && !saw_q) {
            ++cases_p_only;
        } else if (!saw_p && saw_q) {
            ++cases_q_only;
        } else if (saw_p && saw_q) {
            ++cases_both;
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "TOTAL_POINTS="
        << total_points
        << '\n';

    std::cout
        << "ALGEBRA_FAILURES="
        << algebra_failures
        << '\n';

    std::cout
        << "GCD_FORMULA_FAILURES="
        << gcd_formula_failures
        << '\n';

    std::cout
        << "P_CLASSIFICATION_FAILURES="
        << p_classification_failures
        << '\n';

    std::cout
        << "Q_CLASSIFICATION_FAILURES="
        << q_classification_failures
        << '\n';

    std::cout
        << "CASES_P_DIVIDES_S_MINUS_1="
        << cases_p_divides_s_minus_1
        << '\n';

    std::cout
        << "CASES_Q_DIVIDES_S_MINUS_1="
        << cases_q_divides_s_minus_1
        << '\n';

    std::cout
        << "TOTAL_P_HITS="
        << total_p_hits
        << '\n';

    std::cout
        << "TOTAL_Q_HITS="
        << total_q_hits
        << '\n';

    std::cout
        << "TOTAL_P_HITS_PREDICTED="
        << total_p_hits_predicted
        << '\n';

    std::cout
        << "TOTAL_Q_HITS_PREDICTED="
        << total_q_hits_predicted
        << '\n';

    std::cout
        << "CASES_P_HIT="
        << cases_p_hit
        << '\n';

    std::cout
        << "CASES_Q_HIT="
        << cases_q_hit
        << '\n';

    std::cout
        << "CASES_P_ONLY="
        << cases_p_only
        << '\n';

    std::cout
        << "CASES_Q_ONLY="
        << cases_q_only
        << '\n';

    std::cout
        << "CASES_BOTH="
        << cases_both
        << '\n';

    std::cout
        << "CASES_OTHER_GCD="
        << cases_other_gcd
        << '\n';

    std::cout
        << "TOTAL_OTHER_GCDS="
        << total_other_gcds
        << '\n';

    std::cout
        << "MAX_GCD="
        << max_gcd
        << '\n';

    const bool status =
        algebra_failures == 0 &&
        gcd_formula_failures == 0 &&
        p_classification_failures == 0 &&
        q_classification_failures == 0 &&
        total_other_gcds == 0;

    std::cout
        << "COMPLETE_ROOT_CLASSIFICATION_STATUS="
        << (status ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}