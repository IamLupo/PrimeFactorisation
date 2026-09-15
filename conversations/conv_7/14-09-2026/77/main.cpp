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
    u64 D;
};

static mpz_class mpz_from_u64(u64 value) {
    return mpz_class(std::to_string(value));
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        const u64 r = a % b;
        a = b;
        b = r;
    }
    return a;
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

static CaseData make_case(std::mt19937_64& rng) {
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

        const u64 D = N - s * s;

        return {p, q, N, s, D};
    }
}

/*
    Original E(z):

        E(z) = (s-z)(D+s+z(s-1))
*/
static mpz_class E_from_z(
    const CaseData& c,
    u64 z
) {
    const mpz_class S = mpz_from_u64(c.s);
    const mpz_class Z = mpz_from_u64(z);
    const mpz_class D = mpz_from_u64(c.D);

    return (S - Z) *
           (D + S + Z * (S - 1));
}

/*
    Substitute z=s-t:

        E(s-t)
          = t(D+s+(s-t)(s-1))
          = t(N-t(s-1)).
*/
static mpz_class E_from_t(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);
    const mpz_class S = mpz_from_u64(c.s);

    return T * (N - T * (S - 1));
}

static u64 gcd_N_E_from_t(
    const CaseData& c,
    u64 t
) {
    const mpz_class value = E_from_t(c, t);
    const mpz_class N = mpz_from_u64(c.N);

    const mpz_class g = gcd(N, value);

    return g.get_ui();
}

static u64 gcd_N_theoretical(
    const CaseData& c,
    u64 t
) {
    const u64 t2 = t * t;
    const u64 s_minus_1 = c.s - 1;

    const u64 product = t2 * s_minus_1;

    return gcd_u64(c.N, product);
}

int main() {
    constexpr int EXPERIMENT = 366;
    constexpr int CASES = 300;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x36620260914ULL);

    u64 total_identity_checks = 0;
    u64 identity_failures = 0;

    u64 total_gcd_checks = 0;
    u64 gcd_formula_failures = 0;

    u64 cases_with_nontrivial_s_minus_1 = 0;
    u64 cases_p_divides_s_minus_1 = 0;

    u64 cases_found_factor = 0;
    u64 cases_first_hit_equals_p = 0;
    u64 cases_first_hit_equals_other = 0;
    u64 cases_no_nontrivial_hit = 0;

    u64 total_nontrivial_hits = 0;
    u64 total_p_hits = 0;
    u64 total_q_hits = 0;

    u64 total_pre_p_nontrivial_hits = 0;
    u64 cases_pre_p_hit = 0;

    u64 max_first_hit_t = 0;
    u64 max_nontrivial_hits = 0;

    u64 cases_gcd_s_minus_1_nontrivial = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 g_s1 = gcd_u64(N, s - 1);

        if (g_s1 > 1) {
            ++cases_with_nontrivial_s_minus_1;
        }

        if ((s - 1) % p == 0) {
            ++cases_p_divides_s_minus_1;
        }

        /*
            Verify the exact t-substitution identity at every t.
        */
        for (u64 t = 1; t <= s; ++t) {
            const u64 z = s - t;

            const mpz_class from_z =
                E_from_z(c, z);

            const mpz_class from_t =
                E_from_t(c, t);

            ++total_identity_checks;

            if (from_z != from_t) {
                ++identity_failures;
            }
        }

        /*
            Verify

                gcd(N,E(s-t))
                  =
                gcd(N,t^2(s-1))
        */
        u64 first_nontrivial_t = 0;
        u64 first_nontrivial_g = 1;

        u64 nontrivial_hits_this_case = 0;
        u64 p_hits_this_case = 0;
        u64 q_hits_this_case = 0;

        for (u64 t = 1; t <= s; ++t) {
            ++total_gcd_checks;

            const u64 g_actual =
                gcd_N_E_from_t(c, t);

            const u64 g_theory =
                gcd_N_theoretical(c, t);

            if (g_actual != g_theory) {
                ++gcd_formula_failures;
            }

            if (g_actual > 1 && g_actual < N) {
                ++nontrivial_hits_this_case;
                ++total_nontrivial_hits;

                if (first_nontrivial_t == 0) {
                    first_nontrivial_t = t;
                    first_nontrivial_g = g_actual;

                    max_first_hit_t =
                        std::max(
                            max_first_hit_t,
                            t
                        );
                }

                if (g_actual == p) {
                    ++p_hits_this_case;
                    ++total_p_hits;
                }

                if (g_actual == q) {
                    ++q_hits_this_case;
                    ++total_q_hits;
                }
            }

            /*
                Count nontrivial hits before the hidden p-coordinate.
                The hidden coordinate is t=p.
            */
            if (t < p &&
                g_actual > 1 &&
                g_actual < N) {
                ++total_pre_p_nontrivial_hits;
                ++cases_pre_p_hit;
            }
        }

        max_nontrivial_hits =
            std::max(
                max_nontrivial_hits,
                nontrivial_hits_this_case
            );

        if (first_nontrivial_t == 0) {
            ++cases_no_nontrivial_hit;
        } else {
            ++cases_found_factor;

            if (first_nontrivial_g == p) {
                ++cases_first_hit_equals_p;
            } else {
                ++cases_first_hit_equals_other;
            }
        }

        /*
            For p < q and t <= s < q, a q-hit should only be possible
            through the s-1 factor. This explicitly records whether it
            occurs.
        */
        if (q_hits_this_case > 0) {
            /*
                Kept as information; no assertion is made here because
                q | (s-1) is possible in principle.
            */
        }

        (void)p_hits_this_case;
        (void)q_hits_this_case;
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "TOTAL_IDENTITY_CHECKS="
        << total_identity_checks
        << '\n';

    std::cout
        << "IDENTITY_FAILURES="
        << identity_failures
        << '\n';

    std::cout
        << "TOTAL_GCD_CHECKS="
        << total_gcd_checks
        << '\n';

    std::cout
        << "GCD_FORMULA_FAILURES="
        << gcd_formula_failures
        << '\n';

    std::cout
        << "CASES_GCD_N_S_MINUS_1_GT_1="
        << cases_with_nontrivial_s_minus_1
        << '\n';

    std::cout
        << "CASES_P_DIVIDES_S_MINUS_1="
        << cases_p_divides_s_minus_1
        << '\n';

    std::cout
        << "CASES_FOUND_NONTRIVIAL_FACTOR="
        << cases_found_factor
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_EQUALS_P="
        << cases_first_hit_equals_p
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_EQUALS_OTHER="
        << cases_first_hit_equals_other
        << '\n';

    std::cout
        << "CASES_NO_NONTRIVIAL_HIT="
        << cases_no_nontrivial_hit
        << '\n';

    std::cout
        << "TOTAL_NONTRIVIAL_HITS="
        << total_nontrivial_hits
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
        << "TOTAL_PRE_P_NONTRIVIAL_HITS="
        << total_pre_p_nontrivial_hits
        << '\n';

    std::cout
        << "CASES_PRE_P_HIT="
        << cases_pre_p_hit
        << '\n';

    std::cout
        << "MAX_FIRST_HIT_T="
        << max_first_hit_t
        << '\n';

    std::cout
        << "MAX_NONTRIVIAL_HITS="
        << max_nontrivial_hits
        << '\n';

    const bool algebra_ok =
        identity_failures == 0 &&
        gcd_formula_failures == 0;

    const bool first_hit_structure_ok =
        cases_first_hit_equals_other == 0;

    std::cout
        << "ALGEBRA_STATUS="
        << (algebra_ok ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FIRST_HIT_STRUCTURE_STATUS="
        << (first_hit_structure_ok
                ? "PASS"
                : "OTHER_FACTOR_DETECTED")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
