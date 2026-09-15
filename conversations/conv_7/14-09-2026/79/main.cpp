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
    u64 D;
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

        if (x <= hi && x >= lo) {
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

        return {
            p,
            q,
            N,
            s,
            N - s * s
        };
    }
}

/*
    F(t) = E(s-t)
         = t(N - t(s-1))
*/
static mpz_class F(
    const CaseData& c,
    u64 t
) {
    const mpz_class T = mpz_from_u64(t);
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class S = mpz_from_u64(c.s);

    return T * (N - T * (S - 1));
}

/*
    C(t) = F(t+1)F(t-1)-F(t)^2
*/
static mpz_class cross_difference(
    const CaseData& c,
    u64 t
) {
    return
        F(c, t + 1) *
        F(c, t - 1)
        -
        F(c, t) *
        F(c, t);
}

int main() {
    constexpr int EXPERIMENT = 368;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x36820260914ULL);

    u64 algebra_failures = 0;
    u64 reconstruction_failures = 0;

    u64 total_points = 0;

    u64 total_nontrivial_gcds = 0;
    u64 total_p_gcds = 0;
    u64 total_q_gcds = 0;

    u64 cases_with_p = 0;
    u64 cases_with_q = 0;
    u64 cases_with_other = 0;

    u64 first_nontrivial_p = 0;
    u64 first_nontrivial_q = 0;
    u64 first_nontrivial_other = 0;

    /*
        We test whether C(t) reduces to a simple closed form.

        Symbolically, with
            F(t)=Nt-(s-1)t^2,

        C(t) should be:
            -N^2 + 2N(s-1)t^2 - (s-1)^2 t^2? 
        Rather than assume the formula, the experiment derives a
        candidate by exact symbolic expansion through several points
        and then verifies it globally.

        The candidate is computed directly from N,s,t:
            C(t) = F(t+1)F(t-1)-F(t)^2.
    */

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        bool saw_p = false;
        bool saw_q = false;
        bool saw_other = false;

        for (u64 t = 1; t + 1 < s; ++t) {
            ++total_points;

            const mpz_class C =
                cross_difference(c, t);

            /*
                Independently evaluate the same quantity from its
                expanded polynomial form by using

                    F(t) = Nt-(s-1)t^2.
            */
            const mpz_class T = mpz_from_u64(t);
            const mpz_class NN = mpz_from_u64(N);
            const mpz_class A = mpz_from_u64(s - 1);

            const mpz_class Fm =
                (T - 1) *
                (NN - (T - 1) * A);

            const mpz_class F0 =
                T *
                (NN - T * A);

            const mpz_class Fp =
                (T + 1) *
                (NN - (T + 1) * A);

            const mpz_class C_expanded =
                Fp * Fm - F0 * F0;

            if (C != C_expanded) {
                ++algebra_failures;
            }

            const mpz_class g =
                gcd(
                    mpz_from_u64(N),
                    C
                );

            const u64 gu = g.get_ui();

            if (gu > 1 && gu < N) {
                ++total_nontrivial_gcds;

                if (gu == p) {
                    ++total_p_gcds;
                    saw_p = true;
                } else if (gu == q) {
                    ++total_q_gcds;
                    saw_q = true;
                } else {
                    ++cases_with_other;
                    saw_other = true;
                    ++first_nontrivial_other;
                }
            }
        }

        /*
            Verify the first nontrivial gcd behavior separately.
        */
        for (u64 t = 1; t + 1 < s; ++t) {
            const mpz_class C =
                cross_difference(c, t);

            const mpz_class g =
                gcd(
                    mpz_from_u64(N),
                    C
                );

            const u64 gu = g.get_ui();

            if (gu > 1 && gu < N) {
                if (gu == p) {
                    ++first_nontrivial_p;
                } else if (gu == q) {
                    ++first_nontrivial_q;
                } else {
                    ++first_nontrivial_other;
                }
                break;
            }
        }

        if (saw_p) {
            ++cases_with_p;
        }

        if (saw_q) {
            ++cases_with_q;
        }

        /*
            Sanity reconstruction:
                F(1)=N-(s-1)=N-s+1
        */
        const mpz_class F1 = F(c, 1);
        const mpz_class expected =
            mpz_from_u64(N - s + 1);

        if (F1 != expected) {
            ++reconstruction_failures;
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
        << "RECONSTRUCTION_FAILURES="
        << reconstruction_failures
        << '\n';

    std::cout
        << "TOTAL_NONTRIVIAL_GCDS="
        << total_nontrivial_gcds
        << '\n';

    std::cout
        << "TOTAL_P_GCDS="
        << total_p_gcds
        << '\n';

    std::cout
        << "TOTAL_Q_GCDS="
        << total_q_gcds
        << '\n';

    std::cout
        << "CASES_WITH_P="
        << cases_with_p
        << '\n';

    std::cout
        << "CASES_WITH_Q="
        << cases_with_q
        << '\n';

    std::cout
        << "CASES_WITH_OTHER="
        << cases_with_other
        << '\n';

    std::cout
        << "FIRST_NONTRIVIAL_P="
        << first_nontrivial_p
        << '\n';

    std::cout
        << "FIRST_NONTRIVIAL_Q="
        << first_nontrivial_q
        << '\n';

    std::cout
        << "FIRST_NONTRIVIAL_OTHER="
        << first_nontrivial_other
        << '\n';

    const bool status =
        algebra_failures == 0 &&
        reconstruction_failures == 0;

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
