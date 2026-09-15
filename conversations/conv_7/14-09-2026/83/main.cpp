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
    const mpz_class fm = F(c, t - 1);
    const mpz_class f0 = F(c, t);
    const mpz_class fp = F(c, t + 1);

    return fp * fm - f0 * f0;
}

/*
    Exact closed form:

        C(t)
        =
        (s-1)^2(1-2t^2)
        + 2N(s-1)t
        - N^2
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
    Generic root condition:

        2t^2 == 1 (mod r)
*/
static bool root_mod(
    u64 t,
    u64 r
) {
    const u64 x = t % r;
    return ((2ULL * x * x) % r) == 1;
}

int main() {
    constexpr int EXPERIMENT = 372;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37220260914ULL);

    u64 total_points = 0;

    u64 algebra_failures = 0;
    u64 gcd_classification_failures = 0;
    u64 inclusion_exclusion_failures = 0;

    u64 total_boolean_hits = 0;

    u64 total_actual_p_only = 0;
    u64 total_actual_q_only = 0;
    u64 total_actual_both = 0;

    u64 total_predicted_p = 0;
    u64 total_predicted_q = 0;
    u64 total_predicted_both = 0;
    u64 total_predicted_union = 0;

    u64 total_actual_union = 0;

    u64 cases_with_overlap = 0;
    u64 cases_without_overlap = 0;

    u64 max_overlap = 0;
    u64 max_union_hits = 0;
    u64 min_union_hits = UINT64_MAX;

    u64 estimate_10pct = 0;
    u64 estimate_25pct = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        u64 case_actual_p_only = 0;
        u64 case_actual_q_only = 0;
        u64 case_actual_both = 0;

        u64 case_predicted_p = 0;
        u64 case_predicted_q = 0;
        u64 case_predicted_both = 0;

        for (u64 t = 1; t < s; ++t) {
            ++total_points;

            /*
                Algebra verification.
            */
            if (cross_difference(c, t) != C_closed(c, t)) {
                ++algebra_failures;
            }

            /*
                Actual gcd class.
            */
            const mpz_class N_mpz =
                mpz_from_u64(N);

            const mpz_class C_mpz =
                cross_difference(c, t);

            const mpz_class G =
                gcd(N_mpz, C_mpz);

            const u64 g =
                G.get_ui();

            /*
                Predicted modular roots.
            */
            const bool p_hit =
                root_mod(t, p);

            const bool q_hit =
                root_mod(t, q);

            if (p_hit) {
                ++case_predicted_p;
            }

            if (q_hit) {
                ++case_predicted_q;
            }

            if (p_hit && q_hit) {
                ++case_predicted_both;
            }

            /*
                Actual gcd classification.
            */
            if (g == 1) {
                if (p_hit || q_hit) {
                    ++gcd_classification_failures;
                }
            } else if (g == p) {
                ++case_actual_p_only;

                if (!p_hit || q_hit) {
                    ++gcd_classification_failures;
                }
            } else if (g == q) {
                ++case_actual_q_only;

                if (!q_hit || p_hit) {
                    ++gcd_classification_failures;
                }
            } else if (g == N) {
                ++case_actual_both;

                if (!p_hit || !q_hit) {
                    ++gcd_classification_failures;
                }
            } else {
                /*
                    N=pq, so there should be no other gcd.
                */
                ++gcd_classification_failures;
            }
        }

        /*
            Inclusion-exclusion for actual hit positions.
        */
        const u64 case_actual_union =
            case_actual_p_only +
            case_actual_q_only +
            case_actual_both;

        const u64 case_predicted_union =
            case_predicted_p +
            case_predicted_q -
            case_predicted_both;

        if (case_actual_union != case_predicted_union) {
            ++inclusion_exclusion_failures;
        }

        /*
            Global counters.
        */
        total_actual_p_only += case_actual_p_only;
        total_actual_q_only += case_actual_q_only;
        total_actual_both += case_actual_both;

        total_predicted_p += case_predicted_p;
        total_predicted_q += case_predicted_q;
        total_predicted_both += case_predicted_both;
        total_predicted_union += case_predicted_union;

        total_actual_union += case_actual_union;
        total_boolean_hits += case_actual_union;

        if (case_predicted_both > 0) {
            ++cases_with_overlap;
        } else {
            ++cases_without_overlap;
        }

        max_overlap =
            std::max(
                max_overlap,
                case_predicted_both
            );

        max_union_hits =
            std::max(
                max_union_hits,
                case_actual_union
            );

        min_union_hits =
            std::min(
                min_union_hits,
                case_actual_union
            );

        /*
            Factor-blind count estimate.
        */
        if (case_actual_union > 0) {
            const u64 estimate =
                (2 * s) / case_actual_union;

            const u64 diff =
                estimate > p
                    ? estimate - p
                    : p - estimate;

            if (diff * 10 <= p) {
                ++estimate_10pct;
            }

            if (diff * 4 <= p) {
                ++estimate_25pct;
            }
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
        << "TOTAL_BOOLEAN_HITS="
        << total_boolean_hits
        << '\n';

    std::cout
        << "TOTAL_ACTUAL_P_ONLY="
        << total_actual_p_only
        << '\n';

    std::cout
        << "TOTAL_ACTUAL_Q_ONLY="
        << total_actual_q_only
        << '\n';

    std::cout
        << "TOTAL_ACTUAL_BOTH="
        << total_actual_both
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
        << "TOTAL_ACTUAL_UNION="
        << total_actual_union
        << '\n';

    std::cout
        << "ALGEBRA_FAILURES="
        << algebra_failures
        << '\n';

    std::cout
        << "GCD_CLASSIFICATION_FAILURES="
        << gcd_classification_failures
        << '\n';

    std::cout
        << "INCLUSION_EXCLUSION_FAILURES="
        << inclusion_exclusion_failures
        << '\n';

    std::cout
        << "CASES_WITH_OVERLAP="
        << cases_with_overlap
        << '\n';

    std::cout
        << "CASES_WITHOUT_OVERLAP="
        << cases_without_overlap
        << '\n';

    std::cout
        << "MAX_OVERLAP="
        << max_overlap
        << '\n';

    std::cout
        << "MIN_UNION_HITS="
        << min_union_hits
        << '\n';

    std::cout
        << "MAX_UNION_HITS="
        << max_union_hits
        << '\n';

    std::cout
        << "ESTIMATE_WITHIN_10_PERCENT="
        << estimate_10pct
        << '\n';

    std::cout
        << "ESTIMATE_WITHIN_25_PERCENT="
        << estimate_25pct
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