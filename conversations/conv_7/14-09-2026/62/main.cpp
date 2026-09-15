#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

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

static u64 mod_mul(u64 a, u64 b, u64 mod) {
    return static_cast<u64>(
        (static_cast<u128>(a) * static_cast<u128>(b)) %
        static_cast<u128>(mod)
    );
}

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 result = 1 % mod;

    while (e != 0) {
        if (e & 1ULL) {
            result = mod_mul(result, a, mod);
        }

        a = mod_mul(a, a, mod);
        e >>= 1ULL;
    }

    return result;
}

static u64 mod_inverse_prime(u64 a, u64 p) {
    return mod_pow(a % p, p - 2, p);
}

/*
 * Exact polynomial:
 *
 * Q(k) = 4k^2 - (7s+4)k + 3N
 */
static i128 Q_exact(u64 k, u64 s, u64 N) {
    return
        static_cast<i128>(4) *
            static_cast<i128>(k) *
            static_cast<i128>(k)
        -
        static_cast<i128>(7) *
            static_cast<i128>(s) *
            static_cast<i128>(k)
        -
        static_cast<i128>(4) *
            static_cast<i128>(k)
        +
        static_cast<i128>(3) *
            static_cast<i128>(N);
}

static u64 abs_i128_to_u64(i128 x) {
    if (x < 0) {
        x = -x;
    }

    return static_cast<u64>(x);
}

/*
 * Compute Q(k) modulo prime without ever constructing a negative
 * unsigned integer.
 *
 * Q(k) = 4k^2 - (7s+4)k + 3N
 */
static u64 Q_mod_prime(
    u64 k,
    u64 s,
    u64 N,
    u64 prime
) {
    const u64 kk = k % prime;
    const u64 ss = s % prime;
    const u64 nn = N % prime;

    const u64 term1 =
        mod_mul(
            4 % prime,
            mod_mul(kk, kk, prime),
            prime
        );

    const u64 coeff =
        (
            mod_mul(7 % prime, ss, prime) +
            4 % prime
        ) % prime;

    const u64 term2 =
        mod_mul(coeff, kk, prime);

    const u64 term3 =
        mod_mul(3 % prime, nn, prime);

    /*
     * term1 - term2 + term3 mod prime
     */
    u64 result = term1;

    if (result >= term2) {
        result -= term2;
    } else {
        result = prime - (term2 - result);
    }

    result += term3;

    if (result >= prime) {
        result -= prime;
    }

    return result;
}

static u64 Q_mod_N(
    u64 k,
    u64 s,
    u64 N
) {
    const u64 a =
        static_cast<u64>(
            static_cast<u128>(4) *
            static_cast<u128>(k % N) *
            static_cast<u128>(k % N) %
            static_cast<u128>(N)
        );

    const u64 coeff =
        static_cast<u64>(
            (
                static_cast<u128>(7) *
                static_cast<u128>(s % N)
                +
                4
            ) % N
        );

    const u64 b =
        static_cast<u64>(
            static_cast<u128>(coeff) *
            static_cast<u128>(k % N) %
            static_cast<u128>(N)
        );

    const u64 c =
        static_cast<u64>(
            static_cast<u128>(3) *
            static_cast<u128>(N % N)
        );

    (void)c;

    /*
     * Since 3N == 0 (mod N),
     * Q(k) mod N = 4k^2 - (7s+4)k mod N.
     */
    u64 result;

    if (a >= b) {
        result = a - b;
    } else {
        result = N - (b - a);
    }

    return result;
}

static u64 predicted_second_root(
    u64 s,
    u64 prime
) {
    /*
     * Q(k) mod prime =
     *
     * k(4k-(7s+4)).
     *
     * Therefore:
     *
     * k = (7s+4) / 4 mod prime.
     */
    const u64 seven_s =
        mod_mul(
            7 % prime,
            s % prime,
            prime
        );

    const u64 c =
        (seven_s + 4 % prime) % prime;

    const u64 inv4 =
        mod_inverse_prime(4, prime);

    return mod_mul(c, inv4, prime);
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 20000) {
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

static bool second_root_is_valid(
    u64 root,
    u64 prime
) {
    return root < prime;
}

static u64 gcd_Q_N(
    u64 k,
    u64 s,
    u64 N
) {
    /*
     * We only need Q(k) modulo N.
     *
     * gcd(Q(k),N) = gcd(Q(k) mod N,N).
     */
    const u64 residue =
        Q_mod_N(k, s, N);

    return std::gcd(residue, N);
}

static bool is_predicted_root_mod_factor(
    u64 k,
    u64 factor,
    u64 second_root
) {
    return
        (k % factor == 0) ||
        (k % factor == second_root);
}

int main() {
    std::cout << "START EXPERIMENT 351\n";

    constexpr int CASE_COUNT = 300;
    constexpr int PRIME_LIMIT = 20000;

    std::mt19937_64 rng(
        0x351351351ULL
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

    int exact_identity_failures = 0;
    int p_divisibility_failures = 0;
    int q_divisibility_failures = 0;

    int modular_root_failures = 0;

    int p_direct_hits = 0;
    int q_direct_hits = 0;

    int second_p_root_in_scan = 0;
    int second_q_root_in_scan = 0;

    int gcd1_count = 0;
    int gcdp_count = 0;
    int gcdq_count = 0;
    int gcdN_count = 0;
    int gcd_other_count = 0;

    int unexpected_p_hits = 0;
    int unexpected_q_hits = 0;
    int unexpected_other_hits = 0;

    long long total_scan_points = 0;

    for (int case_id = 0;
         case_id < static_cast<int>(cases.size());
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        /*
         * ---------------------------------------------------------
         * Exact polynomial identities
         * ---------------------------------------------------------
         */

        const i128 Qp =
            Q_exact(p, s, N);

        const i128 expected_p =
            static_cast<i128>(p) *
            (
                static_cast<i128>(4) *
                    static_cast<i128>(p)
                +
                static_cast<i128>(3) *
                    static_cast<i128>(q)
                -
                static_cast<i128>(7) *
                    static_cast<i128>(s)
                -
                static_cast<i128>(4)
            );

        if (Qp != expected_p) {
            ++exact_identity_failures;
        }

        const i128 Qq =
            Q_exact(q, s, N);

        const i128 expected_q =
            static_cast<i128>(q) *
            (
                static_cast<i128>(3) *
                    static_cast<i128>(p)
                +
                static_cast<i128>(4) *
                    static_cast<i128>(q)
                -
                static_cast<i128>(7) *
                    static_cast<i128>(s)
                -
                static_cast<i128>(4)
            );

        if (Qq != expected_q) {
            ++exact_identity_failures;
        }

        if (abs_i128_to_u64(Qp) % p != 0) {
            ++p_divisibility_failures;
        }

        if (abs_i128_to_u64(Qq) % q != 0) {
            ++q_divisibility_failures;
        }

        /*
         * ---------------------------------------------------------
         * Modular roots
         * ---------------------------------------------------------
         */

        const u64 rp =
            predicted_second_root(s, p);

        const u64 rq =
            predicted_second_root(s, q);

        if (!second_root_is_valid(rp, p) ||
            !second_root_is_valid(rq, q)) {
            ++modular_root_failures;
        }

        if (Q_mod_prime(rp, s, N, p) != 0) {
            ++modular_root_failures;
        }

        if (Q_mod_prime(rq, s, N, q) != 0) {
            ++modular_root_failures;
        }

        /*
         * k=p is always in the search interval because p<sqrt(N)
         * for p<q.
         */
        if (p <= s) {
            if (gcd_Q_N(p, s, N) == p) {
                ++p_direct_hits;
            }
        }

        /*
         * q is outside [1,s].
         */
        if (q <= s) {
            if (gcd_Q_N(q, s, N) == q) {
                ++q_direct_hits;
            }
        }

        if (rp >= 1 && rp <= s) {
            ++second_p_root_in_scan;
        }

        if (rq >= 1 && rq <= s) {
            ++second_q_root_in_scan;
        }

        /*
         * ---------------------------------------------------------
         * Exhaustive scan
         * ---------------------------------------------------------
         */

        for (u64 k = 1; k <= s; ++k) {
            ++total_scan_points;

            const u64 g =
                gcd_Q_N(k, s, N);

            if (g == 1) {
                ++gcd1_count;
                continue;
            }

            if (g == p) {
                ++gcdp_count;

                if (!is_predicted_root_mod_factor(
                        k, p, rp)) {
                    ++unexpected_p_hits;
                }

                continue;
            }

            if (g == q) {
                ++gcdq_count;

                if (!is_predicted_root_mod_factor(
                        k, q, rq)) {
                    ++unexpected_q_hits;
                }

                continue;
            }

            if (g == N) {
                ++gcdN_count;

                const bool p_root =
                    is_predicted_root_mod_factor(
                        k, p, rp
                    );

                const bool q_root =
                    is_predicted_root_mod_factor(
                        k, q, rq
                    );

                if (!p_root || !q_root) {
                    ++unexpected_other_hits;
                }

                continue;
            }

            ++gcd_other_count;
            ++unexpected_other_hits;
        }

        if (case_id < 5) {
            std::cout
                << "\nCASE=" << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << "\n";

            std::cout
                << "r_p=" << rp
                << " r_q=" << rq
                << "\n";

            std::cout
                << "Q(p)="
                << static_cast<long long>(Qp)
                << "\n";

            std::cout
                << "Q(q)="
                << static_cast<long long>(Qq)
                << "\n";

            std::cout
                << "gcd(Q(p),N)="
                << gcd_Q_N(p, s, N)
                << "\n";

            std::cout
                << "p_root_in_scan="
                << (p <= s ? "YES" : "NO")
                << "\n";

            std::cout
                << "second_p_root_in_scan="
                << ((rp >= 1 && rp <= s) ? "YES" : "NO")
                << "\n";

            std::cout
                << "second_q_root_in_scan="
                << ((rq >= 1 && rq <= s) ? "YES" : "NO")
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
        << "EXACT_IDENTITY_FAILURES="
        << exact_identity_failures
        << "\n";

    std::cout
        << "P_DIVISIBILITY_FAILURES="
        << p_divisibility_failures
        << "\n";

    std::cout
        << "Q_DIVISIBILITY_FAILURES="
        << q_divisibility_failures
        << "\n";

    std::cout
        << "MODULAR_ROOT_FAILURES="
        << modular_root_failures
        << "\n";

    std::cout
        << "P_DIRECT_HITS="
        << p_direct_hits
        << "\n";

    std::cout
        << "Q_DIRECT_HITS="
        << q_direct_hits
        << "\n";

    std::cout
        << "SECOND_P_ROOT_IN_SCAN="
        << second_p_root_in_scan
        << "\n";

    std::cout
        << "SECOND_Q_ROOT_IN_SCAN="
        << second_q_root_in_scan
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

    const bool algebra_status =
        exact_identity_failures == 0 &&
        p_divisibility_failures == 0 &&
        q_divisibility_failures == 0 &&
        modular_root_failures == 0;

    const bool scan_status =
        unexpected_p_hits == 0 &&
        unexpected_q_hits == 0 &&
        unexpected_other_hits == 0 &&
        gcd_other_count == 0;

    std::cout
        << "ALGEBRA_STATUS="
        << (algebra_status ? "PASS" : "FAIL")
        << "\n";

    std::cout
        << "SCAN_STATUS="
        << (scan_status ? "PASS" : "FAIL")
        << "\n";

    if (algebra_status && scan_status) {
        std::cout
            << "STRUCTURE_STATUS="
            << "ALL_OBSERVED_GCD_HITS_EXPLAINED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS="
            << "UNEXPLAINED_BEHAVIOR_FOUND"
            << "\n";
    }

    std::cout << "FINISHED EXPERIMENT 351\n";

    return 0;
}
