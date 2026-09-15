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

    C(t)=F(t+1)F(t-1)-F(t)^2
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
    Test 2t^2 == 1 mod r.
*/
static bool quadratic_root(
    u64 t,
    u64 r
) {
    const u64 a = t % r;
    return ((2ULL * a * a) % r) == 1;
}

/*
    Count roots of

        2t^2 == 1 mod r

    in

        1 <= t < s

    using the explicit residue classes.

    We find one square root by scanning residues.  This is only
    used for the experiment and r <= 10000.
*/
static u64 predicted_hit_count(
    u64 r,
    u64 s
) {
    std::vector<u64> roots;

    for (u64 x = 1; x < r; ++x) {
        if (((2ULL * x * x) % r) == 1) {
            roots.push_back(x);
        }
    }

    u64 count = 0;

    for (u64 root : roots) {
        if (root >= s) {
            continue;
        }

        count += 1 + (s - 1 - root) / r;
    }

    return count;
}

int main() {
    constexpr int EXPERIMENT = 371;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37120260914ULL);

    u64 total_points = 0;
    u64 boolean_hit_points = 0;

    u64 predicted_p_hits = 0;
    u64 predicted_q_hits = 0;

    u64 actual_p_hits = 0;
    u64 actual_q_hits = 0;

    u64 predicted_total_hits = 0;
    u64 actual_total_hits = 0;

    u64 hit_count_failures = 0;

    u64 p_residue_cases = 0;
    u64 p_nonresidue_cases = 0;
    u64 q_residue_cases = 0;
    u64 q_nonresidue_cases = 0;

    u64 cases_p_more_hits = 0;
    u64 cases_q_more_hits = 0;
    u64 cases_equal_hits = 0;

    u64 total_hit_count_error = 0;

    u64 min_p_hits = UINT64_MAX;
    u64 max_p_hits = 0;

    u64 min_q_hits = UINT64_MAX;
    u64 max_q_hits = 0;

    u64 min_total_hits = UINT64_MAX;
    u64 max_total_hits = 0;

    /*
        Blind-oracle diagnostics:
        estimate p from the number of hits assuming the
        p contribution dominates.

        This is NOT claimed to be a factorization formula.
        It simply tests how much information the hit count carries.
    */
    u64 p_estimate_within_10pct = 0;
    u64 p_estimate_within_25pct = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        const u64 predicted_p =
            predicted_hit_count(p, s);

        const u64 predicted_q =
            predicted_hit_count(q, s);

        predicted_p_hits += predicted_p;
        predicted_q_hits += predicted_q;

        predicted_total_hits +=
            predicted_p + predicted_q;

        if (predicted_p > 0) {
            ++p_residue_cases;
        } else {
            ++p_nonresidue_cases;
        }

        if (predicted_q > 0) {
            ++q_residue_cases;
        } else {
            ++q_nonresidue_cases;
        }

        min_p_hits =
            std::min(min_p_hits, predicted_p);

        max_p_hits =
            std::max(max_p_hits, predicted_p);

        min_q_hits =
            std::min(min_q_hits, predicted_q);

        max_q_hits =
            std::max(max_q_hits, predicted_q);

        /*
            Actual boolean hit counts.
        */
        u64 actual_total = 0;

        for (u64 t = 1; t < s; ++t) {
            ++total_points;

            const mpz_class N =
                mpz_from_u64(c.N);

            const mpz_class value =
                C(c, t);

            const mpz_class g =
                gcd(N, value);

            /*
                Deliberately throw away the gcd value.
                Only retain whether it is nontrivial.
            */
            const bool hit =
                (g > 1 && g < N);

            if (hit) {
                ++boolean_hit_points;
                ++actual_total;
            }

            const bool p_hit =
                quadratic_root(t, p);

            const bool q_hit =
                quadratic_root(t, q);

            if (p_hit) {
                ++actual_p_hits;
            }

            if (q_hit) {
                ++actual_q_hits;
            }
        }

        /*
            The actual factor-specific root counts computed independently
            from modular arithmetic must agree with the predictions.
        */
        const u64 case_actual_p =
            predicted_hit_count(p, s);

        const u64 case_actual_q =
            predicted_hit_count(q, s);

        if (case_actual_p != predicted_p ||
            case_actual_q != predicted_q ||
            actual_total != predicted_p + predicted_q) {
            ++hit_count_failures;
        }

        actual_total_hits += actual_total;

        const u64 predicted_case_total =
            predicted_p + predicted_q;

        const u64 error =
            actual_total > predicted_case_total
                ? actual_total - predicted_case_total
                : predicted_case_total - actual_total;

        total_hit_count_error += error;

        min_total_hits =
            std::min(min_total_hits, actual_total);

        max_total_hits =
            std::max(max_total_hits, actual_total);

        if (predicted_p > predicted_q) {
            ++cases_p_more_hits;
        } else if (predicted_q > predicted_p) {
            ++cases_q_more_hits;
        } else {
            ++cases_equal_hits;
        }

        /*
            Very crude estimate:

                p ~= 2s / H

            where H is the total number of observed roots.

            This is deliberately tested empirically rather than assumed.
        */
        if (actual_total > 0) {
            const u64 estimate =
                (2 * s) / actual_total;

            const u64 diff =
                estimate > p
                    ? estimate - p
                    : p - estimate;

            if (diff * 10 <= p) {
                ++p_estimate_within_10pct;
            }

            if (diff * 4 <= p) {
                ++p_estimate_within_25pct;
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
        << "BOOLEAN_HIT_POINTS="
        << boolean_hit_points
        << '\n';

    std::cout
        << "PREDICTED_P_HITS="
        << predicted_p_hits
        << '\n';

    std::cout
        << "PREDICTED_Q_HITS="
        << predicted_q_hits
        << '\n';

    std::cout
        << "ACTUAL_TOTAL_HITS="
        << actual_total_hits
        << '\n';

    std::cout
        << "PREDICTED_TOTAL_HITS="
        << predicted_total_hits
        << '\n';

    std::cout
        << "HIT_COUNT_FAILURES="
        << hit_count_failures
        << '\n';

    std::cout
        << "TOTAL_HIT_COUNT_ERROR="
        << total_hit_count_error
        << '\n';

    std::cout
        << "P_RESIDUE_CASES="
        << p_residue_cases
        << '\n';

    std::cout
        << "P_NONRESIDUE_CASES="
        << p_nonresidue_cases
        << '\n';

    std::cout
        << "Q_RESIDUE_CASES="
        << q_residue_cases
        << '\n';

    std::cout
        << "Q_NONRESIDUE_CASES="
        << q_nonresidue_cases
        << '\n';

    std::cout
        << "CASES_P_MORE_HITS="
        << cases_p_more_hits
        << '\n';

    std::cout
        << "CASES_Q_MORE_HITS="
        << cases_q_more_hits
        << '\n';

    std::cout
        << "CASES_EQUAL_HITS="
        << cases_equal_hits
        << '\n';

    std::cout
        << "MIN_P_HITS="
        << min_p_hits
        << '\n';

    std::cout
        << "MAX_P_HITS="
        << max_p_hits
        << '\n';

    std::cout
        << "MIN_Q_HITS="
        << min_q_hits
        << '\n';

    std::cout
        << "MAX_Q_HITS="
        << max_q_hits
        << '\n';

    std::cout
        << "MIN_TOTAL_HITS="
        << min_total_hits
        << '\n';

    std::cout
        << "MAX_TOTAL_HITS="
        << max_total_hits
        << '\n';

    std::cout
        << "P_ESTIMATE_WITHIN_10_PERCENT="
        << p_estimate_within_10pct
        << '\n';

    std::cout
        << "P_ESTIMATE_WITHIN_25_PERCENT="
        << p_estimate_within_25pct
        << '\n';

    const bool status =
        hit_count_failures == 0 &&
        total_hit_count_error == 0 &&
        actual_total_hits == predicted_total_hits;

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
