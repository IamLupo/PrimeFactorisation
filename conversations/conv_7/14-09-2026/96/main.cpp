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

static u64 gcd_mpz_u64(
    const mpz_class& a,
    const mpz_class& b
) {
    const mpz_class g = gcd(a, b);
    return g.get_ui();
}

/*
    gcd(N,C(t))
*/
static u64 original_gcd(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class value = C(c, t);

    return gcd_mpz_u64(N, value);
}

/*
    gcd(N,(s-1)^2)
*/
static u64 s_square_gcd(
    const CaseData& c
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);

    const mpz_class value = A * A;

    return gcd_mpz_u64(N, value);
}

/*
    gcd(N,2t^2-1)
*/
static u64 quadratic_gcd(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        2 * T * T - 1;

    return gcd_mpz_u64(N, value);
}

/*
    gcd(N,(s-1)^2(2t^2-1))
*/
static u64 combined_gcd(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        A * A *
        (2 * T * T - 1);

    return gcd_mpz_u64(N, value);
}

int main() {
    constexpr int EXPERIMENT = 385;
    constexpr int CASES = 1000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38520260914ULL);

    u64 total_points = 0;

    u64 original_hits = 0;
    u64 quadratic_hits = 0;
    u64 combined_hits = 0;

    u64 combined_identity_failures = 0;

    u64 original_quadratic_same = 0;
    u64 original_s_factor_only = 0;
    u64 original_quadratic_only = 0;
    u64 original_both_sources = 0;

    u64 cases_with_s_factor = 0;

    u64 total_s_factor_gcd = 0;

    u64 p_from_s_factor = 0;
    u64 q_from_s_factor = 0;
    u64 n_from_s_factor = 0;

    u64 p_from_quadratic = 0;
    u64 q_from_quadratic = 0;

    u64 cases_original_no_quadratic =
        0;

    u64 cases_original_no_s_factor =
        0;

    u64 cases_original_explained =
        0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;

        const u64 s_factor_gcd =
            s_square_gcd(c);

        if (s_factor_gcd > 1 &&
            s_factor_gcd < N) {
            ++cases_with_s_factor;
        }

        total_s_factor_gcd +=
            s_factor_gcd;

        if (s_factor_gcd == p) {
            ++p_from_s_factor;
        } else if (s_factor_gcd == q) {
            ++q_from_s_factor;
        } else if (s_factor_gcd == N) {
            ++n_from_s_factor;
        }

        for (u64 t = 1;
             t < c.s;
             ++t) {

            ++total_points;

            const u64 g_original =
                original_gcd(c, t);

            const u64 g_quadratic =
                quadratic_gcd(c, t);

            const u64 g_combined =
                combined_gcd(c, t);

            /*
                Fundamental identity test.
            */
            if (g_original != g_combined) {
                ++combined_identity_failures;
            }

            const bool original_hit =
                g_original > 1 &&
                g_original < N;

            const bool quadratic_hit =
                g_quadratic > 1 &&
                g_quadratic < N;

            const bool combined_hit =
                g_combined > 1 &&
                g_combined < N;

            if (original_hit) {
                ++original_hits;
            }

            if (quadratic_hit) {
                ++quadratic_hits;
            }

            if (combined_hit) {
                ++combined_hits;
            }

            if (!original_hit) {
                continue;
            }

            const bool s_source =
                (s_factor_gcd > 1 &&
                 s_factor_gcd < N);

            const bool q_source =
                quadratic_hit;

            if (s_source && q_source) {
                ++original_both_sources;
            } else if (s_source) {
                ++original_s_factor_only;
            } else if (q_source) {
                ++original_quadratic_only;
            } else {
                /*
                    This is the interesting impossible-looking
                    remainder.
                */
            }

            if (q_source) {
                if (g_quadratic == p) {
                    ++p_from_quadratic;
                } else if (g_quadratic == q) {
                    ++q_from_quadratic;
                }
            }
        }

        /*
            Case-level diagnostics.
        */
        bool has_original =
            false;

        bool has_quadratic =
            false;

        for (u64 t = 1;
             t < c.s;
             ++t) {

            const u64 go =
                original_gcd(c, t);

            const u64 gq =
                quadratic_gcd(c, t);

            if (go > 1 && go < N) {
                has_original = true;
            }

            if (gq > 1 && gq < N) {
                has_quadratic = true;
            }
        }

        if (has_original &&
            !has_quadratic) {
            ++cases_original_no_quadratic;
        }

        if (has_original &&
            s_factor_gcd == 1) {
            ++cases_original_no_s_factor;
        }

        if (!has_original ||
            has_quadratic ||
            s_factor_gcd > 1) {
            ++cases_original_explained;
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
        << "ORIGINAL_HITS="
        << original_hits
        << '\n';

    std::cout
        << "QUADRATIC_HITS="
        << quadratic_hits
        << '\n';

    std::cout
        << "COMBINED_HITS="
        << combined_hits
        << '\n';

    std::cout
        << "COMBINED_IDENTITY_FAILURES="
        << combined_identity_failures
        << '\n';

    std::cout
        << "ORIGINAL_QUADRATIC_ONLY="
        << original_quadratic_only
        << '\n';

    std::cout
        << "ORIGINAL_S_FACTOR_ONLY="
        << original_s_factor_only
        << '\n';

    std::cout
        << "ORIGINAL_BOTH_SOURCES="
        << original_both_sources
        << '\n';

    std::cout
        << "CASES_WITH_S_FACTOR="
        << cases_with_s_factor
        << '\n';

    std::cout
        << "TOTAL_S_FACTOR_GCD="
        << total_s_factor_gcd
        << '\n';

    std::cout
        << "P_FROM_S_FACTOR="
        << p_from_s_factor
        << '\n';

    std::cout
        << "Q_FROM_S_FACTOR="
        << q_from_s_factor
        << '\n';

    std::cout
        << "N_FROM_S_FACTOR="
        << n_from_s_factor
        << '\n';

    std::cout
        << "P_FROM_QUADRATIC="
        << p_from_quadratic
        << '\n';

    std::cout
        << "Q_FROM_QUADRATIC="
        << q_from_quadratic
        << '\n';

    std::cout
        << "CASES_ORIGINAL_WITHOUT_QUADRATIC="
        << cases_original_no_quadratic
        << '\n';

    std::cout
        << "CASES_ORIGINAL_WITHOUT_S_FACTOR="
        << cases_original_no_s_factor
        << '\n';

    std::cout
        << "CASES_ORIGINAL_EXPLAINED="
        << cases_original_explained
        << '\n';

    const bool status =
        combined_identity_failures == 0;

    std::cout
        << "EXACT_GCD_DECOMPOSITION_STATUS="
        << (status ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
