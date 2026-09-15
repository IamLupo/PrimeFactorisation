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
    Return the smallest positive root t < r of

        2t^2 == 1 (mod r).

    r <= 10000, so direct enumeration is fine here.
*/
static u64 smallest_root(u64 r) {
    for (u64 t = 1; t < r; ++t) {
        if (((2ULL * t * t) % r) == 1) {
            return t;
        }
    }

    return 0;
}

int main() {
    constexpr int EXPERIMENT = 386;
    constexpr int CASES = 2000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38620260914ULL);

    u64 generic_cases = 0;
    u64 no_p_root_cases = 0;

    u64 total_smallest_p_root = 0;
    u64 total_smallest_q_root = 0;

    u64 min_p_root = UINT64_MAX;
    u64 max_p_root = 0;

    u64 min_q_root = UINT64_MAX;
    u64 max_q_root = 0;

    u64 cases_p_root_lt_s = 0;
    u64 cases_q_root_lt_s = 0;

    u64 cases_p_root_lt_p_over_2 = 0;
    u64 cases_q_root_lt_q_over_2 = 0;

    /*
        Compare the smallest root with simple geometric quantities.
    */
    u64 cases_p_root_eq_s_minus_p = 0;
    u64 cases_p_root_eq_s_mod_p = 0;

    /*
        Test whether the smallest root correlates with the
        factor deviations from sqrt(N).
    */
    u64 p_root_lt_s_minus_p = 0;
    u64 p_root_ge_s_minus_p = 0;

    u64 p_root_lt_q_minus_s = 0;
    u64 p_root_ge_q_minus_s = 0;

    /*
        Distribution bins for root/p and root/s.
    */
    u64 p_root_lt_p_10pct = 0;
    u64 p_root_lt_p_25pct = 0;
    u64 p_root_lt_p_50pct = 0;

    u64 p_root_lt_s_10pct = 0;
    u64 p_root_lt_s_25pct = 0;
    u64 p_root_lt_s_50pct = 0;

    /*
        Actual direct scan length needed to find p.
    */
    u64 total_scan_to_p_root = 0;
    u64 min_scan_to_p_root = UINT64_MAX;
    u64 max_scan_to_p_root = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        const bool p_residue =
            (p % 8 == 1 || p % 8 == 7);

        const bool q_residue =
            (q % 8 == 1 || q % 8 == 7);

        /*
            This experiment is specifically about the generic
            quadratic-root branch.
        */
        if (!p_residue) {
            ++no_p_root_cases;
            continue;
        }

        ++generic_cases;

        const u64 rp =
            smallest_root(p);

        const u64 rq =
            q_residue
                ? smallest_root(q)
                : 0;

        if (rp == 0) {
            ++no_p_root_cases;
            continue;
        }

        total_smallest_p_root += rp;

        min_p_root =
            std::min(min_p_root, rp);

        max_p_root =
            std::max(max_p_root, rp);

        if (rq != 0) {
            total_smallest_q_root += rq;

            min_q_root =
                std::min(min_q_root, rq);

            max_q_root =
                std::max(max_q_root, rq);

            if (rq < s) {
                ++cases_q_root_lt_s;
            }

            if (2 * rq < q) {
                ++cases_q_root_lt_q_over_2;
            }
        }

        if (rp < s) {
            ++cases_p_root_lt_s;
        }

        if (2 * rp < p) {
            ++cases_p_root_lt_p_over_2;
        }

        /*
            Compare to simple factor geometry.
        */
        const u64 gap_p =
            s > p
                ? s - p
                : 0;

        const u64 gap_q =
            q > s
                ? q - s
                : 0;

        if (rp == gap_p) {
            ++cases_p_root_eq_s_minus_p;
        }

        if (rp % p == s % p) {
            ++cases_p_root_eq_s_mod_p;
        }

        if (rp < gap_p) {
            ++p_root_lt_s_minus_p;
        } else {
            ++p_root_ge_s_minus_p;
        }

        if (rp < gap_q) {
            ++p_root_lt_q_minus_s;
        } else {
            ++p_root_ge_q_minus_s;
        }

        /*
            Root relative to p.
        */
        if (10 * rp < p) {
            ++p_root_lt_p_10pct;
        }

        if (4 * rp < p) {
            ++p_root_lt_p_25pct;
        }

        if (2 * rp < p) {
            ++p_root_lt_p_50pct;
        }

        /*
            Root relative to s.
        */
        if (10 * rp < s) {
            ++p_root_lt_s_10pct;
        }

        if (4 * rp < s) {
            ++p_root_lt_s_25pct;
        }

        if (2 * rp < s) {
            ++p_root_lt_s_50pct;
        }

        /*
            The actual direct scan reaches the first p-root at rp.
        */
        total_scan_to_p_root += rp;

        min_scan_to_p_root =
            std::min(
                min_scan_to_p_root,
                rp
            );

        max_scan_to_p_root =
            std::max(
                max_scan_to_p_root,
                rp
            );
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "GENERIC_P_ROOT_CASES="
        << generic_cases
        << '\n';

    std::cout
        << "NO_P_ROOT_CASES="
        << no_p_root_cases
        << '\n';

    std::cout
        << "P_ROOT_LT_S="
        << cases_p_root_lt_s
        << '\n';

    std::cout
        << "P_ROOT_LT_P_OVER_2="
        << cases_p_root_lt_p_over_2
        << '\n';

    std::cout
        << "Q_ROOT_LT_S="
        << cases_q_root_lt_s
        << '\n';

    std::cout
        << "Q_ROOT_LT_Q_OVER_2="
        << cases_q_root_lt_q_over_2
        << '\n';

    std::cout
        << "TOTAL_SMALLEST_P_ROOT="
        << total_smallest_p_root
        << '\n';

    std::cout
        << "MIN_P_ROOT="
        << min_p_root
        << '\n';

    std::cout
        << "MAX_P_ROOT="
        << max_p_root
        << '\n';

    std::cout
        << "TOTAL_SMALLEST_Q_ROOT="
        << total_smallest_q_root
        << '\n';

    std::cout
        << "MIN_Q_ROOT="
        << min_q_root
        << '\n';

    std::cout
        << "MAX_Q_ROOT="
        << max_q_root
        << '\n';

    std::cout
        << "P_ROOT_EQ_S_MINUS_P="
        << cases_p_root_eq_s_minus_p
        << '\n';

    std::cout
        << "P_ROOT_EQ_S_MOD_P="
        << cases_p_root_eq_s_mod_p
        << '\n';

    std::cout
        << "P_ROOT_LT_S_MINUS_P="
        << p_root_lt_s_minus_p
        << '\n';

    std::cout
        << "P_ROOT_GE_S_MINUS_P="
        << p_root_ge_s_minus_p
        << '\n';

    std::cout
        << "P_ROOT_LT_Q_MINUS_S="
        << p_root_lt_q_minus_s
        << '\n';

    std::cout
        << "P_ROOT_GE_Q_MINUS_S="
        << p_root_ge_q_minus_s
        << '\n';

    std::cout
        << "P_ROOT_LT_10_PERCENT_OF_P="
        << p_root_lt_p_10pct
        << '\n';

    std::cout
        << "P_ROOT_LT_25_PERCENT_OF_P="
        << p_root_lt_p_25pct
        << '\n';

    std::cout
        << "P_ROOT_LT_50_PERCENT_OF_P="
        << p_root_lt_p_50pct
        << '\n';

    std::cout
        << "P_ROOT_LT_10_PERCENT_OF_S="
        << p_root_lt_s_10pct
        << '\n';

    std::cout
        << "P_ROOT_LT_25_PERCENT_OF_S="
        << p_root_lt_s_25pct
        << '\n';

    std::cout
        << "P_ROOT_LT_50_PERCENT_OF_S="
        << p_root_lt_s_50pct
        << '\n';

    std::cout
        << "TOTAL_SCAN_TO_P_ROOT="
        << total_scan_to_p_root
        << '\n';

    std::cout
        << "MIN_SCAN_TO_P_ROOT="
        << min_scan_to_p_root
        << '\n';

    std::cout
        << "MAX_SCAN_TO_P_ROOT="
        << max_scan_to_p_root
        << '\n';

    std::cout
        << "AVERAGE_SCAN_TO_P_ROOT="
        << (
            generic_cases > 0
                ? static_cast<double>(
                    total_scan_to_p_root
                  ) /
                  static_cast<double>(
                    generic_cases
                  )
                : 0.0
        )
        << '\n';

    const bool root_status =
        generic_cases == 0 ||
        cases_p_root_lt_s == generic_cases;

    std::cout
        << "P_ROOT_FINITE_SCAN_STATUS="
        << (
            root_status
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
