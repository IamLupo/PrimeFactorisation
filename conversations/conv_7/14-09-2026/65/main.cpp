#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        __builtin_sqrtl(static_cast<long double>(n))
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(
        static_cast<std::size_t>(limit + 1),
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<long long>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 50000) {
            usable.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        usable.size() - 1
    );

    std::vector<CaseData> result;
    result.reserve(static_cast<std::size_t>(count));

    while (static_cast<int>(result.size()) < count) {
        const std::size_t i = dist(rng);
        const std::size_t j = dist(rng);

        if (i == j) {
            continue;
        }

        const u64 p =
            static_cast<u64>(
                std::min(usable[i], usable[j])
            );

        const u64 q =
            static_cast<u64>(
                std::max(usable[i], usable[j])
            );

        const u64 N = p * q;

        if (!used_N.insert(N).second) {
            continue;
        }

        result.push_back({
            p,
            q,
            N,
            isqrt_u64(N)
        });
    }

    return result;
}

static i128 linear_exact(
    u64 r,
    u64 s
) {
    /*
     * L(r) = 4r - (7s+4)
     */
    return
        static_cast<i128>(4) *
        static_cast<i128>(r)
        -
        (
            static_cast<i128>(7) *
            static_cast<i128>(s)
            +
            static_cast<i128>(4)
        );
}

static u64 abs_i128_to_u64(i128 x) {
    if (x < 0) {
        x = -x;
    }

    return static_cast<u64>(x);
}

static u64 linear_gcd(
    u64 r,
    u64 s,
    u64 N
) {
    const i128 L =
        linear_exact(r, s);

    const u64 value =
        abs_i128_to_u64(L);

    return std::gcd(value, N);
}

static u64 mod_pow(
    u64 a,
    u64 e,
    u64 mod
) {
    u64 result = 1 % mod;
    a %= mod;

    while (e != 0) {
        if (e & 1ULL) {
            result = static_cast<u64>(
                (
                    static_cast<u128>(result) *
                    static_cast<u128>(a)
                ) %
                static_cast<u128>(mod)
            );
        }

        a = static_cast<u64>(
            (
                static_cast<u128>(a) *
                static_cast<u128>(a)
            ) %
            static_cast<u128>(mod)
        );

        e >>= 1ULL;
    }

    return result;
}

static u64 second_root(
    u64 s,
    u64 p
) {
    /*
     * 4r = 7s+4 (mod p)
     */
    const u64 numerator =
        (
            static_cast<u64>(
                (
                    static_cast<u128>(7) *
                    static_cast<u128>(s % p)
                ) %
                static_cast<u128>(p)
            )
            +
            4
        ) % p;

    const u64 inv4 =
        mod_pow(4, p - 2, p);

    return static_cast<u64>(
        (
            static_cast<u128>(numerator) *
            static_cast<u128>(inv4)
        ) %
        static_cast<u128>(p)
    );
}

static bool is_linear_root_for_factor(
    u64 r,
    u64 factor,
    u64 second_root_value
) {
    return
        (r % factor) == second_root_value;
}

int main() {
    std::cout << "START EXPERIMENT 354\n";

    constexpr int CASE_COUNT = 1000;
    constexpr int PRIME_LIMIT = 50000;

    std::mt19937_64 rng(
        0x354354354ULL
    );

    const std::vector<int> primes =
        sieve_primes(PRIME_LIMIT);

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    u64 total_scan_points = 0;

    u64 gcd1_count = 0;
    u64 gcdp_count = 0;
    u64 gcdq_count = 0;
    u64 gcdN_count = 0;
    u64 gcd_other_count = 0;

    u64 unexpected_p_hits = 0;
    u64 unexpected_q_hits = 0;
    u64 unexpected_other_hits = 0;

    u64 cases_with_p_hit = 0;
    u64 cases_with_q_hit = 0;
    u64 cases_with_N_hit = 0;

    u64 cases_with_exactly_one_nontrivial_hit = 0;
    u64 cases_with_exactly_two_nontrivial_hits = 0;
    u64 cases_with_exactly_three_nontrivial_hits = 0;

    u64 min_nontrivial_hits = UINT64_MAX;
    u64 max_nontrivial_hits = 0;
    u64 total_nontrivial_hits = 0;

    u64 p_root_found_exactly = 0;
    u64 q_root_found_exactly = 0;

    u64 p_root_mismatch = 0;
    u64 q_root_mismatch = 0;

    u64 linear_vs_quadratic_same_factor_hits = 0;

    u64 smallest_nontrivial_is_p = 0;
    u64 smallest_nontrivial_is_q = 0;
    u64 smallest_nontrivial_is_N = 0;

    u64 all_nontrivial_hits_are_predicted = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 rp =
            second_root(s, p);

        const u64 rq =
            second_root(s, q);

        u64 nontrivial_hits = 0;

        u64 first_nontrivial_r = 0;
        u64 first_nontrivial_g = 0;

        bool found_p = false;
        bool found_q = false;
        bool found_N = false;

        bool all_predicted = true;

        /*
         * ---------------------------------------------------------
         * Linear scan:
         *
         *     gcd(4r-(7s+4), N)
         *
         * ---------------------------------------------------------
         */

        for (u64 r = 1; r <= s; ++r) {
            ++total_scan_points;

            const u64 g =
                linear_gcd(r, s, N);

            if (g == 1) {
                ++gcd1_count;
                continue;
            }

            ++nontrivial_hits;

            if (first_nontrivial_r == 0) {
                first_nontrivial_r = r;
                first_nontrivial_g = g;
            }

            if (g == p) {
                ++gcdp_count;
                found_p = true;

                /*
                 * Since the linear expression vanishes modulo p,
                 * r must be the second root modulo p.
                 */
                if (!is_linear_root_for_factor(
                        r,
                        p,
                        rp)) {

                    ++unexpected_p_hits;
                    all_predicted = false;
                }

                continue;
            }

            if (g == q) {
                ++gcdq_count;
                found_q = true;

                if (!is_linear_root_for_factor(
                        r,
                        q,
                        rq)) {

                    ++unexpected_q_hits;
                    all_predicted = false;
                }

                continue;
            }

            if (g == N) {
                ++gcdN_count;
                found_N = true;

                /*
                 * N means the same r is simultaneously a
                 * second root modulo p and q.
                 */
                const bool p_ok =
                    is_linear_root_for_factor(
                        r,
                        p,
                        rp
                    );

                const bool q_ok =
                    is_linear_root_for_factor(
                        r,
                        q,
                        rq
                    );

                if (!p_ok || !q_ok) {
                    ++unexpected_other_hits;
                    all_predicted = false;
                }

                continue;
            }

            ++gcd_other_count;
            ++unexpected_other_hits;
            all_predicted = false;
        }

        if (found_p) {
            ++cases_with_p_hit;
        }

        if (found_q) {
            ++cases_with_q_hit;
        }

        if (found_N) {
            ++cases_with_N_hit;
        }

        total_nontrivial_hits +=
            nontrivial_hits;

        min_nontrivial_hits =
            std::min(
                min_nontrivial_hits,
                nontrivial_hits
            );

        max_nontrivial_hits =
            std::max(
                max_nontrivial_hits,
                nontrivial_hits
            );

        if (nontrivial_hits == 1) {
            ++cases_with_exactly_one_nontrivial_hit;
        }

        if (nontrivial_hits == 2) {
            ++cases_with_exactly_two_nontrivial_hits;
        }

        if (nontrivial_hits == 3) {
            ++cases_with_exactly_three_nontrivial_hits;
        }

        if (all_predicted) {
            ++all_nontrivial_hits_are_predicted;
        }

        /*
         * ---------------------------------------------------------
         * Verify that r_p itself is found.
         * ---------------------------------------------------------
         */

        if (rp >= 1 && rp <= s) {
            const u64 g =
                linear_gcd(rp, s, N);

            if (g == p) {
                ++p_root_found_exactly;
            } else if (g == N) {
                ++p_root_found_exactly;
            } else {
                ++p_root_mismatch;
            }
        } else {
            ++p_root_mismatch;
        }

        /*
         * q's second root may or may not be in the scan interval.
         */
        if (rq >= 1 && rq <= s) {
            const u64 g =
                linear_gcd(rq, s, N);

            if (g == q || g == N) {
                ++q_root_found_exactly;
            } else {
                ++q_root_mismatch;
            }
        }

        /*
         * ---------------------------------------------------------
         * Direct comparison with quadratic Q.
         *
         * At the nonzero root class, the quadratic and linear
         * constructions detect the same factor.
         * ---------------------------------------------------------
         */

        if (rp >= 1 && rp <= s) {
            const u64 gl =
                linear_gcd(rp, s, N);

            /*
             * Q(rp) has the same modular root because:
             *
             * Q(k) = k L(k) + 3N.
             *
             * Therefore gcd(Q(rp),N) should equal the linear
             * gcd at rp.
             */
            i128 q_value =
                static_cast<i128>(4) *
                static_cast<i128>(rp) *
                static_cast<i128>(rp)
                -
                (
                    static_cast<i128>(7) *
                    static_cast<i128>(s)
                    +
                    static_cast<i128>(4)
                ) *
                static_cast<i128>(rp)
                +
                static_cast<i128>(3) *
                static_cast<i128>(N);

            const u64 q_abs =
                abs_i128_to_u64(q_value);

            const u64 gq =
                std::gcd(q_abs, N);

            if (gl == gq) {
                ++linear_vs_quadratic_same_factor_hits;
            }
        }

        /*
         * Determine which factor appears at the first
         * nontrivial linear hit.
         */
        if (first_nontrivial_r != 0) {
            if (first_nontrivial_g == p) {
                ++smallest_nontrivial_is_p;
            } else if (first_nontrivial_g == q) {
                ++smallest_nontrivial_is_q;
            } else if (first_nontrivial_g == N) {
                ++smallest_nontrivial_is_N;
            }
        }

        if (case_id < 10) {
            std::cout
                << "\nCASE="
                << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << "\n";

            std::cout
                << "r_p="
                << rp
                << "\n";

            std::cout
                << "r_q="
                << rq
                << "\n";

            std::cout
                << "r_p_in_scan="
                << (
                    rp >= 1 && rp <= s
                        ? "YES"
                        : "NO"
                )
                << "\n";

            std::cout
                << "r_q_in_scan="
                << (
                    rq >= 1 && rq <= s
                        ? "YES"
                        : "NO"
                )
                << "\n";

            std::cout
                << "NONTRIVIAL_HITS="
                << nontrivial_hits
                << "\n";

            std::cout
                << "FIRST_HIT_R="
                << first_nontrivial_r
                << "\n";

            std::cout
                << "FIRST_HIT_GCD="
                << first_nontrivial_g
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "TOTAL_SCAN_POINTS="
        << total_scan_points
        << "\n";

    std::cout
        << "GCD_1_COUNT="
        << gcd1_count
        << "\n";

    std::cout
        << "GCD_P_COUNT="
        << gcdp_count
        << "\n";

    std::cout
        << "GCD_Q_COUNT="
        << gcdq_count
        << "\n";

    std::cout
        << "GCD_N_COUNT="
        << gcdN_count
        << "\n";

    std::cout
        << "GCD_OTHER_COUNT="
        << gcd_other_count
        << "\n";

    std::cout
        << "CASES_WITH_P_HIT="
        << cases_with_p_hit
        << "\n";

    std::cout
        << "CASES_WITH_Q_HIT="
        << cases_with_q_hit
        << "\n";

    std::cout
        << "CASES_WITH_N_HIT="
        << cases_with_N_hit
        << "\n";

    std::cout
        << "P_ROOT_FOUND_EXACTLY="
        << p_root_found_exactly
        << "\n";

    std::cout
        << "P_ROOT_MISMATCH="
        << p_root_mismatch
        << "\n";

    std::cout
        << "Q_ROOT_FOUND_EXACTLY="
        << q_root_found_exactly
        << "\n";

    std::cout
        << "Q_ROOT_MISMATCH="
        << q_root_mismatch
        << "\n";

    std::cout
        << "TOTAL_NONTRIVIAL_HITS="
        << total_nontrivial_hits
        << "\n";

    std::cout
        << "AVG_NONTRIVIAL_HITS="
        << total_nontrivial_hits /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "MIN_NONTRIVIAL_HITS="
        << min_nontrivial_hits
        << "\n";

    std::cout
        << "MAX_NONTRIVIAL_HITS="
        << max_nontrivial_hits
        << "\n";

    std::cout
        << "CASES_EXACTLY_1_HIT="
        << cases_with_exactly_one_nontrivial_hit
        << "\n";

    std::cout
        << "CASES_EXACTLY_2_HITS="
        << cases_with_exactly_two_nontrivial_hits
        << "\n";

    std::cout
        << "CASES_EXACTLY_3_HITS="
        << cases_with_exactly_three_nontrivial_hits
        << "\n";

    std::cout
        << "SMALLEST_HIT_IS_P="
        << smallest_nontrivial_is_p
        << "\n";

    std::cout
        << "SMALLEST_HIT_IS_Q="
        << smallest_nontrivial_is_q
        << "\n";

    std::cout
        << "SMALLEST_HIT_IS_N="
        << smallest_nontrivial_is_N
        << "\n";

    std::cout
        << "LINEAR_QUADRATIC_SAME_FACTOR="
        << linear_vs_quadratic_same_factor_hits
        << "\n";

    std::cout
        << "UNEXPECTED_P_HITS="
        << unexpected_p_hits
        << "\n";

    std::cout
        << "UNEXPECTED_Q_HITS="
        << unexpected_q_hits
        << "\n";

    std::cout
        << "UNEXPECTED_OTHER_HITS="
        << unexpected_other_hits
        << "\n";

    std::cout
        << "ALL_HITS_PREDICTED="
        << all_nontrivial_hits_are_predicted
        << "\n";

    const bool structure_ok =
        p_root_mismatch == 0 &&
        q_root_mismatch == 0 &&
        gcd_other_count == 0 &&
        unexpected_p_hits == 0 &&
        unexpected_q_hits == 0 &&
        unexpected_other_hits == 0 &&
        all_nontrivial_hits_are_predicted ==
            cases.size();

    if (structure_ok) {
        std::cout
            << "STRUCTURE_STATUS="
            << "LINEAR_ROOT_STRUCTURE_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS="
            << "UNEXPLAINED_BEHAVIOR_FOUND"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 354\n";

    return 0;
}
