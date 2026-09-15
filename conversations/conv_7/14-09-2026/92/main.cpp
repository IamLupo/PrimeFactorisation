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

/*
    Factor-blind Boolean oracle.
*/
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
    Recover a candidate factor from a hit position:

        gcd(N, 2t^2-1)
*/
static u64 factor_from_hit(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class polynomial =
        2 * T * T - 1;

    const mpz_class g =
        gcd(N, polynomial);

    return g.get_ui();
}

int main() {
    constexpr int EXPERIMENT = 381;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38120260914ULL);

    u64 total_scan_points = 0;

    u64 cases_with_hit = 0;
    u64 cases_no_hit = 0;

    u64 total_first_hit_t = 0;
    u64 min_first_hit_t = UINT64_MAX;
    u64 max_first_hit_t = 0;

    u64 cases_first_hit_gives_p = 0;
    u64 cases_first_hit_gives_q = 0;
    u64 cases_first_hit_gives_n = 0;
    u64 cases_first_hit_gives_one = 0;
    u64 cases_first_hit_gives_other = 0;

    u64 cases_exception_p = 0;
    u64 exception_first_hit_p_recovered = 0;
    u64 exception_gcd_s_minus_1_recovers_p = 0;

    u64 total_hit_tests = 0;

    /*
        Also test every Boolean hit position, not just the first.
    */
    u64 total_hit_positions = 0;
    u64 hit_position_factor_success = 0;
    u64 hit_position_factor_failure = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const bool p_exception =
            ((s - 1) % p == 0);

        if (p_exception) {
            ++cases_exception_p;
        }

        /*
            First check the exceptional branch independently:

                gcd(N,s-1)=p
        */
        if (p_exception) {
            const mpz_class N_mpz =
                mpz_from_u64(N);

            const mpz_class A_mpz =
                mpz_from_u64(s - 1);

            const mpz_class g =
                gcd(N_mpz, A_mpz);

            if (g.get_ui() == p) {
                ++exception_gcd_s_minus_1_recovers_p;
            }
        }

        bool found = false;
        u64 first_hit_t = 0;
        u64 first_hit_factor = 0;

        /*
            Scan t=1,...,s-1 until first Boolean hit.
        */
        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;
            ++total_hit_tests;

            if (!boolean_hit(c, t)) {
                continue;
            }

            ++cases_with_hit;

            first_hit_t = t;
            found = true;

            first_hit_factor =
                factor_from_hit(c, t);

            if (first_hit_factor == p) {
                ++cases_first_hit_gives_p;
            } else if (first_hit_factor == q) {
                ++cases_first_hit_gives_q;
            } else if (first_hit_factor == N) {
                ++cases_first_hit_gives_n;
            } else if (first_hit_factor == 1) {
                ++cases_first_hit_gives_one;
            } else {
                ++cases_first_hit_gives_other;
            }

            if (p_exception &&
                first_hit_factor == p) {
                ++exception_first_hit_p_recovered;
            }

            break;
        }

        if (!found) {
            ++cases_no_hit;
        } else {
            total_first_hit_t += first_hit_t;

            min_first_hit_t =
                std::min(
                    min_first_hit_t,
                    first_hit_t
                );

            max_first_hit_t =
                std::max(
                    max_first_hit_t,
                    first_hit_t
                );
        }

        /*
            Now test every hit position independently.
        */
        for (u64 t = 1; t < s; ++t) {
            if (!boolean_hit(c, t)) {
                continue;
            }

            ++total_hit_positions;

            const u64 recovered =
                factor_from_hit(c, t);

            /*
                In the generic branch we expect the recovered factor
                to be one of p or q.

                In the exceptional p | (s-1) branch, the hit itself
                does not necessarily imply p | 2t^2-1.
            */
            if (recovered == p ||
                recovered == q ||
                recovered == N) {

                ++hit_position_factor_success;
            } else {
                ++hit_position_factor_failure;
            }
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
        << "CASES_WITH_HIT="
        << cases_with_hit
        << '\n';

    std::cout
        << "CASES_NO_HIT="
        << cases_no_hit
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
            cases_with_hit > 0
                ? static_cast<double>(
                    total_first_hit_t
                  ) /
                  static_cast<double>(
                    cases_with_hit
                  )
                : 0.0
        )
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_GIVES_P="
        << cases_first_hit_gives_p
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_GIVES_Q="
        << cases_first_hit_gives_q
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_GIVES_N="
        << cases_first_hit_gives_n
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_GIVES_ONE="
        << cases_first_hit_gives_one
        << '\n';

    std::cout
        << "CASES_FIRST_HIT_GIVES_OTHER="
        << cases_first_hit_gives_other
        << '\n';

    std::cout
        << "CASES_P_EXCEPTION="
        << cases_exception_p
        << '\n';

    std::cout
        << "EXCEPTION_GCD_S_MINUS_1_RECOVERS_P="
        << exception_gcd_s_minus_1_recovers_p
        << '\n';

    std::cout
        << "EXCEPTION_FIRST_HIT_P_RECOVERED="
        << exception_first_hit_p_recovered
        << '\n';

    std::cout
        << "TOTAL_HIT_POSITIONS="
        << total_hit_positions
        << '\n';

    std::cout
        << "HIT_POSITION_FACTOR_SUCCESS="
        << hit_position_factor_success
        << '\n';

    std::cout
        << "HIT_POSITION_FACTOR_FAILURE="
        << hit_position_factor_failure
        << '\n';

    const bool first_hit_success =
        cases_first_hit_gives_p > 0 ||
        cases_first_hit_gives_q > 0 ||
        cases_first_hit_gives_n > 0;

    std::cout
        << "FIRST_HIT_RECOVERY_STATUS="
        << (
            first_hit_success
                ? "OBSERVED"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
