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
static mpz_class C(
    const CaseData& c,
    u64 t
) {
    const mpz_class fm = F(c, t - 1);
    const mpz_class f0 = F(c, t);
    const mpz_class fp = F(c, t + 1);

    return fp * fm - f0 * f0;
}

/*
    Exact closed form:

        C(t)
        =
        A^2(1-2t^2)
        + 2ANt
        - N^2

    where A=s-1.
*/
static mpz_class C_closed(
    const CaseData& c,
    u64 t
) {
    const mpz_class A = mpz_from_u64(c.s - 1);
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);

    return
        A * A * (1 - 2 * T * T)
        + 2 * A * N * T
        - N * N;
}

/*
    Theoretical gcd from the congruence:

        C(t) = A^2(1-2t^2) mod N
*/
static u64 theoretical_gcd(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        A * A * (1 - 2 * T * T);

    const mpz_class g = gcd(N, value);

    return g.get_ui();
}

/*
    Check whether

        2t^2 == 1 mod r.
*/
static bool satisfies_root(
    u64 t,
    u64 r
) {
    const u64 t_mod = t % r;

    /*
        Avoid multiplication overflow.
        r <= 10000 in this experiment.
    */
    const u64 value =
        (2ULL * t_mod * t_mod) % r;

    return value == 1 % r;
}

int main() {
    constexpr int EXPERIMENT = 370;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37020260914ULL);

    u64 total_points = 0;

    u64 algebra_failures = 0;
    u64 gcd_formula_failures = 0;

    u64 p_classification_failures = 0;
    u64 q_classification_failures = 0;

    u64 total_p_hits = 0;
    u64 total_q_hits = 0;

    u64 predicted_p_hits = 0;
    u64 predicted_q_hits = 0;

    u64 cases_p_hit = 0;
    u64 cases_q_hit = 0;

    u64 cases_both = 0;
    u64 cases_p_only = 0;
    u64 cases_q_only = 0;

    u64 cases_other = 0;
    u64 total_other = 0;

    u64 cases_p_not_residue = 0;
    u64 cases_q_not_residue = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;

        bool saw_p = false;
        bool saw_q = false;

        for (u64 t = 1; t < c.s; ++t) {
            ++total_points;

            /*
                Verify exact algebra.
            */
            if (C(c, t) != C_closed(c, t)) {
                ++algebra_failures;
            }

            /*
                Verify gcd formula.
            */
            const mpz_class N_mpz =
                mpz_from_u64(N);

            const mpz_class C_mpz =
                C(c, t);

            const mpz_class actual_g_mpz =
                gcd(N_mpz, C_mpz);

            const u64 actual =
                actual_g_mpz.get_ui();

            const u64 theory =
                theoretical_gcd(c, t);

            if (actual != theory) {
                ++gcd_formula_failures;
            }

            const bool predicted_p =
                satisfies_root(t, p);

            const bool predicted_q =
                satisfies_root(t, q);

            const bool actual_p =
                (actual % p == 0);

            const bool actual_q =
                (actual % q == 0);

            /*
                If A=s-1 is divisible by a factor, that factor
                divides C for every t. Record this separately.
            */
            const bool p_divides_A =
                ((c.s - 1) % p == 0);

            const bool q_divides_A =
                ((c.s - 1) % q == 0);

            const bool predicted_p_full =
                p_divides_A || predicted_p;

            const bool predicted_q_full =
                q_divides_A || predicted_q;

            if (actual_p != predicted_p_full) {
                ++p_classification_failures;
            }

            if (actual_q != predicted_q_full) {
                ++q_classification_failures;
            }

            /*
                Count only the generic quadratic-residue roots.
            */
            if (predicted_p) {
                ++predicted_p_hits;
            }

            if (predicted_q) {
                ++predicted_q_hits;
            }

            if (actual_p) {
                ++total_p_hits;
                saw_p = true;
            }

            if (actual_q) {
                ++total_q_hits;
                saw_q = true;
            }

            if (actual > 1 &&
                actual < N &&
                actual != p &&
                actual != q) {
                ++total_other;
            }

            /*
                A useful diagnostic: if the prime does not admit
                a square root of 1/2, there can be no generic hit.
            */
            if (!predicted_p &&
                !((c.s - 1) % p == 0)) {
                ++cases_p_not_residue;
            }

            if (!predicted_q &&
                !((c.s - 1) % q == 0)) {
                ++cases_q_not_residue;
            }
        }

        if (saw_p) {
            ++cases_p_hit;
        }

        if (saw_q) {
            ++cases_q_hit;
        }

        if (saw_p && saw_q) {
            ++cases_both;
        } else if (saw_p) {
            ++cases_p_only;
        } else if (saw_q) {
            ++cases_q_only;
        }

        if (total_other > 0) {
            /*
                This is intentionally cumulative only for the
                final diagnostic below.
            */
        }
    }

    /*
        Recompute whether any unexpected proper gcd was observed.
        The global total is sufficient because N=pq.
    */
    if (total_other > 0) {
        cases_other = 1;
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
        << "TOTAL_P_HITS="
        << total_p_hits
        << '\n';

    std::cout
        << "TOTAL_Q_HITS="
        << total_q_hits
        << '\n';

    std::cout
        << "P_ROOTS_PREDICTED="
        << predicted_p_hits
        << '\n';

    std::cout
        << "Q_ROOTS_PREDICTED="
        << predicted_q_hits
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
        << "CASES_BOTH="
        << cases_both
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
        << "TOTAL_OTHER_GCDS="
        << total_other
        << '\n';

    std::cout
        << "CASES_WITH_NO_GENERIC_P_ROOT="
        << cases_p_not_residue
        << '\n';

    std::cout
        << "CASES_WITH_NO_GENERIC_Q_ROOT="
        << cases_q_not_residue
        << '\n';

    const bool status =
        algebra_failures == 0 &&
        gcd_formula_failures == 0 &&
        p_classification_failures == 0 &&
        q_classification_failures == 0 &&
        total_other == 0;

    std::cout
        << "STATUS="
        << (status ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
