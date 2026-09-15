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
    C(t) = F(t+1)F(t-1)-F(t)^2
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
    Exact algebraic form:

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
    Generic root test:

        2t^2 == 1 (mod r)
*/
static bool root_mod(
    u64 t,
    u64 r
) {
    const u64 x = t % r;
    return ((2ULL * x * x) % r) == 1;
}

/*
    Complete factor-hit prediction.

    If r | (s-1), then r divides C(t) for EVERY t.

    Otherwise:
        r | C(t) iff 2t^2 == 1 mod r.
*/
static bool factor_hit_predicted(
    u64 t,
    u64 r,
    u64 s
) {
    if ((s - 1) % r == 0) {
        return true;
    }

    return root_mod(t, r);
}

int main() {
    constexpr int EXPERIMENT = 373;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37320260914ULL);

    u64 total_points = 0;

    u64 algebra_failures = 0;
    u64 gcd_classification_failures = 0;

    u64 total_p_only = 0;
    u64 total_q_only = 0;
    u64 total_both = 0;
    u64 total_gcd_one = 0;

    u64 total_predicted_p = 0;
    u64 total_predicted_q = 0;
    u64 total_predicted_both = 0;
    u64 total_predicted_union = 0;

    u64 total_actual_union = 0;

    u64 cases_p_divides_s_minus_1 = 0;
    u64 cases_q_divides_s_minus_1 = 0;

    u64 cases_with_p_exception = 0;
    u64 cases_with_q_exception = 0;

    u64 p_exception_points = 0;
    u64 q_exception_points = 0;

    u64 inclusion_exclusion_failures = 0;

    u64 min_union_hits = UINT64_MAX;
    u64 max_union_hits = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const bool p_exception =
            ((s - 1) % p == 0);

        const bool q_exception =
            ((s - 1) % q == 0);

        if (p_exception) {
            ++cases_p_divides_s_minus_1;
        }

        if (q_exception) {
            ++cases_q_divides_s_minus_1;
        }

        u64 case_p_only = 0;
        u64 case_q_only = 0;
        u64 case_both = 0;
        u64 case_union = 0;

        u64 case_predicted_p = 0;
        u64 case_predicted_q = 0;
        u64 case_predicted_both = 0;

        for (u64 t = 1; t < s; ++t) {
            ++total_points;

            /*
                Verify exact algebra.
            */
            if (C(c, t) != C_closed(c, t)) {
                ++algebra_failures;
            }

            /*
                Actual gcd.
            */
            const mpz_class N_mpz =
                mpz_from_u64(N);

            const mpz_class C_mpz =
                C(c, t);

            const mpz_class G =
                gcd(N_mpz, C_mpz);

            const u64 g =
                G.get_ui();

            /*
                Complete prediction.
            */
            const bool predicted_p =
                factor_hit_predicted(
                    t,
                    p,
                    s
                );

            const bool predicted_q =
                factor_hit_predicted(
                    t,
                    q,
                    s
                );

            if (predicted_p) {
                ++case_predicted_p;
            }

            if (predicted_q) {
                ++case_predicted_q;
            }

            if (predicted_p && predicted_q) {
                ++case_predicted_both;
            }

            /*
                Actual gcd class.
            */
            if (g == 1) {
                ++total_gcd_one;

                if (predicted_p || predicted_q) {
                    ++gcd_classification_failures;
                }
            } else if (g == p) {
                ++case_p_only;

                if (!predicted_p || predicted_q) {
                    ++gcd_classification_failures;
                }
            } else if (g == q) {
                ++case_q_only;

                if (!predicted_q || predicted_p) {
                    ++gcd_classification_failures;
                }
            } else if (g == N) {
                ++case_both;

                if (!predicted_p || !predicted_q) {
                    ++gcd_classification_failures;
                }
            } else {
                ++gcd_classification_failures;
            }

            /*
                Count specifically the exceptional branch.
            */
            if (p_exception) {
                ++p_exception_points;

                if (!predicted_p) {
                    ++gcd_classification_failures;
                }
            }

            if (q_exception) {
                ++q_exception_points;

                if (!predicted_q) {
                    ++gcd_classification_failures;
                }
            }
        }

        /*
            Actual and predicted union counts.
        */
        case_union =
            case_p_only +
            case_q_only +
            case_both;

        const u64 predicted_union =
            case_predicted_p +
            case_predicted_q -
            case_predicted_both;

        if (case_union != predicted_union) {
            ++inclusion_exclusion_failures;
        }

        total_p_only += case_p_only;
        total_q_only += case_q_only;
        total_both += case_both;

        total_predicted_p +=
            case_predicted_p;

        total_predicted_q +=
            case_predicted_q;

        total_predicted_both +=
            case_predicted_both;

        total_predicted_union +=
            predicted_union;

        total_actual_union +=
            case_union;

        min_union_hits =
            std::min(
                min_union_hits,
                case_union
            );

        max_union_hits =
            std::max(
                max_union_hits,
                case_union
            );

        if (p_exception) {
            ++cases_with_p_exception;
        }

        if (q_exception) {
            ++cases_with_q_exception;
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
        << "TOTAL_GCD_ONE="
        << total_gcd_one
        << '\n';

    std::cout
        << "TOTAL_P_ONLY="
        << total_p_only
        << '\n';

    std::cout
        << "TOTAL_Q_ONLY="
        << total_q_only
        << '\n';

    std::cout
        << "TOTAL_BOTH="
        << total_both
        << '\n';

    std::cout
        << "TOTAL_ACTUAL_UNION="
        << total_actual_union
        << '\n';

    std::cout
        << "TOTAL_PREDICTED_P="
        << total_predicted_p
        << '\n';

    std::cout
        << "TOTAL_PREDICTED_Q="
        << total_predicted_q
        << '\n';

    std::cout
        << "TOTAL_PREDICTED_BOTH="
        << total_predicted_both
        << '\n';

    std::cout
        << "TOTAL_PREDICTED_UNION="
        << total_predicted_union
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
        << "CASES_WITH_P_EXCEPTION="
        << cases_with_p_exception
        << '\n';

    std::cout
        << "CASES_WITH_Q_EXCEPTION="
        << cases_with_q_exception
        << '\n';

    std::cout
        << "P_EXCEPTION_POINTS="
        << p_exception_points
        << '\n';

    std::cout
        << "Q_EXCEPTION_POINTS="
        << q_exception_points
        << '\n';

    std::cout
        << "INCLUSION_EXCLUSION_FAILURES="
        << inclusion_exclusion_failures
        << '\n';

    std::cout
        << "GCD_CLASSIFICATION_FAILURES="
        << gcd_classification_failures
        << '\n';

    std::cout
        << "ALGEBRA_FAILURES="
        << algebra_failures
        << '\n';

    std::cout
        << "MIN_UNION_HITS="
        << min_union_hits
        << '\n';

    std::cout
        << "MAX_UNION_HITS="
        << max_union_hits
        << '\n';

    const bool status =
        algebra_failures == 0 &&
        gcd_classification_failures == 0 &&
        inclusion_exclusion_failures == 0 &&
        total_actual_union == total_predicted_union;

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
