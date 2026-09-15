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
    Recover factor from:

        gcd(N, 2t^2-1)
*/
static u64 factor_from_hit(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        2 * T * T - 1;

    const mpz_class g =
        gcd(N, value);

    return g.get_ui();
}

/*
    For odd prime r:

        2 is a quadratic residue mod r
        iff r == 1 or 7 mod 8.
*/
static bool predicts_generic_root(
    u64 r
) {
    const u64 m = r % 8;

    return m == 1 || m == 7;
}

int main() {
    constexpr int EXPERIMENT = 382;
    constexpr int CASES = 2000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38220260914ULL);

    u64 total_scan_points = 0;

    u64 cases_hit = 0;
    u64 cases_no_hit = 0;

    u64 cases_p_exception = 0;

    u64 cases_predicted_hit = 0;
    u64 cases_predicted_no_hit = 0;

    u64 prediction_failures = 0;

    u64 first_hit_p = 0;
    u64 first_hit_q = 0;
    u64 first_hit_one = 0;
    u64 first_hit_other = 0;

    u64 cases_both_factors_residue = 0;
    u64 cases_only_p_residue = 0;
    u64 cases_only_q_residue = 0;
    u64 cases_neither_residue = 0;

    u64 no_hit_with_exception = 0;
    u64 hit_without_predicted_residue = 0;
    u64 predicted_hit_without_actual = 0;

    u64 total_first_hit_t = 0;
    u64 min_first_hit_t = UINT64_MAX;
    u64 max_first_hit_t = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        const bool p_residue =
            predicts_generic_root(p);

        const bool q_residue =
            predicts_generic_root(q);

        const bool p_exception =
            ((s - 1) % p == 0);

        if (p_exception) {
            ++cases_p_exception;
        }

        if (p_residue && q_residue) {
            ++cases_both_factors_residue;
        } else if (p_residue) {
            ++cases_only_p_residue;
        } else if (q_residue) {
            ++cases_only_q_residue;
        } else {
            ++cases_neither_residue;
        }

        /*
            Generic prediction of at least one hit.

            The exceptional p|s-1 branch also guarantees hits.
        */
        const bool predicted_hit =
            p_exception ||
            p_residue ||
            q_residue;

        if (predicted_hit) {
            ++cases_predicted_hit;
        } else {
            ++cases_predicted_no_hit;
        }

        bool found = false;
        u64 first_t = 0;
        u64 recovered = 0;

        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;

            if (!boolean_hit(c, t)) {
                continue;
            }

            found = true;
            first_t = t;

            recovered =
                factor_from_hit(c, t);

            break;
        }

        if (found) {
            ++cases_hit;

            total_first_hit_t +=
                first_t;

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

            if (recovered == p) {
                ++first_hit_p;
            } else if (recovered == q) {
                ++first_hit_q;
            } else if (recovered == 1) {
                ++first_hit_one;
            } else {
                ++first_hit_other;
            }

            if (!predicted_hit) {
                ++hit_without_predicted_residue;
            }
        } else {
            ++cases_no_hit;

            if (p_exception) {
                ++no_hit_with_exception;
            }

            if (predicted_hit) {
                ++predicted_hit_without_actual;
            }
        }

        if (found != predicted_hit) {
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
        << "CASES_P_EXCEPTION="
        << cases_p_exception
        << '\n';

    std::cout
        << "CASES_PREDICTED_HIT="
        << cases_predicted_hit
        << '\n';

    std::cout
        << "CASES_PREDICTED_NO_HIT="
        << cases_predicted_no_hit
        << '\n';

    std::cout
        << "PREDICTION_FAILURES="
        << prediction_failures
        << '\n';

    std::cout
        << "CASES_BOTH_FACTORS_RESIDUE="
        << cases_both_factors_residue
        << '\n';

    std::cout
        << "CASES_ONLY_P_RESIDUE="
        << cases_only_p_residue
        << '\n';

    std::cout
        << "CASES_ONLY_Q_RESIDUE="
        << cases_only_q_residue
        << '\n';

    std::cout
        << "CASES_NEITHER_RESIDUE="
        << cases_neither_residue
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
        << "HIT_WITHOUT_PREDICTED_RESIDUE="
        << hit_without_predicted_residue
        << '\n';

    std::cout
        << "PREDICTED_HIT_WITHOUT_ACTUAL="
        << predicted_hit_without_actual
        << '\n';

    std::cout
        << "NO_HIT_WITH_EXCEPTION="
        << no_hit_with_exception
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
        prediction_failures == 0 &&
        hit_without_predicted_residue == 0 &&
        predicted_hit_without_actual == 0 &&
        no_hit_with_exception == 0;

    std::cout
        << "MOD8_HIT_CLASSIFICATION_STATUS="
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
