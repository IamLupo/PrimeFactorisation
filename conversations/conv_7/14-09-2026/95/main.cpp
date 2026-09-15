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
    Original method:

        gcd(N,C(t))

    Return the full gcd.
*/
static u64 gcd_original(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class value = C(c, t);

    const mpz_class g = gcd(N, value);

    return g.get_ui();
}

/*
    Direct method:

        gcd(N, 2t^2-1)
*/
static u64 gcd_direct(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        2 * T * T - 1;

    const mpz_class g = gcd(N, value);

    return g.get_ui();
}

/*
    Original boolean hit.
*/
static bool original_hit(
    u64 g,
    u64 N
) {
    return g > 1 && g < N;
}

/*
    Direct boolean hit.
*/
static bool direct_hit(
    u64 g,
    u64 N
) {
    return g > 1 && g < N;
}

int main() {
    constexpr int EXPERIMENT = 384;
    constexpr int CASES = 1000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38420260914ULL);

    u64 total_points = 0;

    u64 original_hits = 0;
    u64 direct_hits = 0;

    u64 hit_position_mismatches = 0;
    u64 gcd_value_mismatches = 0;

    u64 cases_hit_set_mismatch = 0;
    u64 cases_recovered_factor_mismatch = 0;

    u64 cases_no_hit_original = 0;
    u64 cases_no_hit_direct = 0;

    u64 original_first_hit_p = 0;
    u64 original_first_hit_q = 0;
    u64 original_first_hit_other = 0;
    u64 original_first_hit_none = 0;

    u64 direct_first_hit_p = 0;
    u64 direct_first_hit_q = 0;
    u64 direct_first_hit_other = 0;
    u64 direct_first_hit_none = 0;

    u64 cases_p_exception = 0;
    u64 exception_original_hits = 0;
    u64 exception_direct_hits = 0;

    u64 total_original_first_t = 0;
    u64 total_direct_first_t = 0;

    u64 min_original_first_t = UINT64_MAX;
    u64 max_original_first_t = 0;

    u64 min_direct_first_t = UINT64_MAX;
    u64 max_direct_first_t = 0;

    /*
        Directly compare the full hit sets in every case.
    */
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
            ++cases_p_exception;
        }

        std::vector<u64> original_positions;
        std::vector<u64> direct_positions;

        original_positions.reserve(s);
        direct_positions.reserve(s);

        for (u64 t = 1; t < s; ++t) {
            ++total_points;

            const u64 g_original =
                gcd_original(c, t);

            const u64 g_direct =
                gcd_direct(c, t);

            const bool hit_original =
                original_hit(g_original, N);

            const bool hit_direct =
                direct_hit(g_direct, N);

            if (hit_original) {
                ++original_hits;
                original_positions.push_back(t);
            }

            if (hit_direct) {
                ++direct_hits;
                direct_positions.push_back(t);
            }

            /*
                In generic cases the two gcd values need not be
                numerically equal if C contains an extra factor,
                but the factor-relevant part should agree.

                We therefore separately test whether their gcds
                have the same nontrivial factor content.
            */
            if (hit_original != hit_direct) {
                ++hit_position_mismatches;
            }

            if (hit_original && hit_direct) {
                /*
                    If both are proper divisors of N, they must be
                    p or q. Compare them directly.
                */
                if (g_original != g_direct) {
                    ++gcd_value_mismatches;
                }
            }
        }

        /*
            Full hit-set equality.
        */
        if (original_positions != direct_positions) {
            ++cases_hit_set_mismatch;
        }

        /*
            First-hit comparison.
        */
        if (original_positions.empty()) {
            ++original_first_hit_none;
            ++cases_no_hit_original;
        } else {
            const u64 t =
                original_positions.front();

            total_original_first_t += t;

            min_original_first_t =
                std::min(
                    min_original_first_t,
                    t
                );

            max_original_first_t =
                std::max(
                    max_original_first_t,
                    t
                );

            const u64 g =
                gcd_original(c, t);

            if (g == p) {
                ++original_first_hit_p;
            } else if (g == q) {
                ++original_first_hit_q;
            } else {
                ++original_first_hit_other;
            }
        }

        if (direct_positions.empty()) {
            ++direct_first_hit_none;
            ++cases_no_hit_direct;
        } else {
            const u64 t =
                direct_positions.front();

            total_direct_first_t += t;

            min_direct_first_t =
                std::min(
                    min_direct_first_t,
                    t
                );

            max_direct_first_t =
                std::max(
                    max_direct_first_t,
                    t
                );

            const u64 g =
                gcd_direct(c, t);

            if (g == p) {
                ++direct_first_hit_p;
            } else if (g == q) {
                ++direct_first_hit_q;
            } else {
                ++direct_first_hit_other;
            }
        }

        /*
            Exceptional branch diagnostic.
        */
        if (p_exception) {
            if (!original_positions.empty()) {
                ++exception_original_hits;
            }

            if (!direct_positions.empty()) {
                ++exception_direct_hits;
            }
        }

        /*
            The first recovered factor should agree whenever both
            methods have a hit and the gcd is nontrivial.
        */
        if (!original_positions.empty() &&
            !direct_positions.empty()) {

            const u64 t_original =
                original_positions.front();

            const u64 t_direct =
                direct_positions.front();

            const u64 g_original =
                gcd_original(c, t_original);

            const u64 g_direct =
                gcd_direct(c, t_direct);

            if (g_original != g_direct) {
                ++cases_recovered_factor_mismatch;
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
        << "ORIGINAL_HITS="
        << original_hits
        << '\n';

    std::cout
        << "DIRECT_HITS="
        << direct_hits
        << '\n';

    std::cout
        << "HIT_POSITION_MISMATCHES="
        << hit_position_mismatches
        << '\n';

    std::cout
        << "GCD_VALUE_MISMATCHES="
        << gcd_value_mismatches
        << '\n';

    std::cout
        << "CASES_HIT_SET_MISMATCH="
        << cases_hit_set_mismatch
        << '\n';

    std::cout
        << "CASES_RECOVERED_FACTOR_MISMATCH="
        << cases_recovered_factor_mismatch
        << '\n';

    std::cout
        << "CASES_NO_HIT_ORIGINAL="
        << cases_no_hit_original
        << '\n';

    std::cout
        << "CASES_NO_HIT_DIRECT="
        << cases_no_hit_direct
        << '\n';

    std::cout
        << "ORIGINAL_FIRST_HIT_P="
        << original_first_hit_p
        << '\n';

    std::cout
        << "ORIGINAL_FIRST_HIT_Q="
        << original_first_hit_q
        << '\n';

    std::cout
        << "ORIGINAL_FIRST_HIT_OTHER="
        << original_first_hit_other
        << '\n';

    std::cout
        << "ORIGINAL_FIRST_HIT_NONE="
        << original_first_hit_none
        << '\n';

    std::cout
        << "DIRECT_FIRST_HIT_P="
        << direct_first_hit_p
        << '\n';

    std::cout
        << "DIRECT_FIRST_HIT_Q="
        << direct_first_hit_q
        << '\n';

    std::cout
        << "DIRECT_FIRST_HIT_OTHER="
        << direct_first_hit_other
        << '\n';

    std::cout
        << "DIRECT_FIRST_HIT_NONE="
        << direct_first_hit_none
        << '\n';

    std::cout
        << "CASES_P_EXCEPTION="
        << cases_p_exception
        << '\n';

    std::cout
        << "EXCEPTION_ORIGINAL_HITS="
        << exception_original_hits
        << '\n';

    std::cout
        << "EXCEPTION_DIRECT_HITS="
        << exception_direct_hits
        << '\n';

    std::cout
        << "TOTAL_ORIGINAL_FIRST_T="
        << total_original_first_t
        << '\n';

    std::cout
        << "TOTAL_DIRECT_FIRST_T="
        << total_direct_first_t
        << '\n';

    std::cout
        << "MIN_ORIGINAL_FIRST_T="
        << min_original_first_t
        << '\n';

    std::cout
        << "MAX_ORIGINAL_FIRST_T="
        << max_original_first_t
        << '\n';

    std::cout
        << "MIN_DIRECT_FIRST_T="
        << min_direct_first_t
        << '\n';

    std::cout
        << "MAX_DIRECT_FIRST_T="
        << max_direct_first_t
        << '\n';

    const bool generic_equivalence =
        hit_position_mismatches == 0 &&
        cases_hit_set_mismatch == 0;

    std::cout
        << "GENERIC_SIGNAL_EQUIVALENCE_STATUS="
        << (
            generic_equivalence
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
