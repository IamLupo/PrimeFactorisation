#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;

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
        (static_cast<__uint128_t>(a) *
         static_cast<__uint128_t>(b)) %
        static_cast<__uint128_t>(mod)
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

static u64 mod_inverse(u64 a, u64 mod) {
    /*
     * Fermat because p and q are prime.
     */
    return mod_pow(a % mod, mod - 2, mod);
}

static u64 q_polynomial(
    u64 k,
    u64 s,
    u64 N
) {
    /*
     * Q(k) = 4 k^2 - (7s+4) k + 3N.
     *
     * The chosen parameter range is small enough that this fits in u64,
     * but the intermediate product is computed in __uint128_t.
     */
    const __uint128_t value =
        static_cast<__uint128_t>(4) *
            static_cast<__uint128_t>(k) *
            static_cast<__uint128_t>(k)
        -
        static_cast<__uint128_t>(7 * s + 4) *
            static_cast<__uint128_t>(k)
        +
        static_cast<__uint128_t>(3) *
            static_cast<__uint128_t>(N);

    return static_cast<u64>(value);
}

static u64 H_polynomial(
    u64 x,
    u64 s,
    u64 D
) {
    /*
     * H(x)=4x^2-(s+4)x+3D-3s
     */
    const __int128_t value =
        static_cast<__int128_t>(4) *
            static_cast<__int128_t>(x) *
            static_cast<__int128_t>(x)
        -
        static_cast<__int128_t>(s + 4) *
            static_cast<__int128_t>(x)
        +
        static_cast<__int128_t>(3) *
            static_cast<__int128_t>(D)
        -
        static_cast<__int128_t>(3) *
            static_cast<__int128_t>(s);

    return static_cast<u64>(value);
}

static u64 predicted_second_root(
    u64 s,
    u64 prime
) {
    /*
     * Q(k) mod prime =
     * k(4k-(7s+4)).
     *
     * The second root is
     *
     * k = (7s+4) / 4 mod prime.
     */
    const u64 c = (7 * (s % prime) + 4) % prime;
    const u64 inv4 = mod_inverse(4, prime);

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

        u64 p = static_cast<u64>(
            std::min(usable[i], usable[j])
        );

        u64 q = static_cast<u64>(
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

int main() {
    std::cout << "START EXPERIMENT 350\n";

    constexpr int CASE_COUNT = 300;
    constexpr int PRIME_LIMIT = 20000;

    std::mt19937_64 rng(
        0x350350350ULL
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

    int identity_failures = 0;
    int p_divisibility_failures = 0;
    int q_divisibility_failures = 0;

    int p_gcd_failures = 0;
    int q_gcd_hits = 0;
    int q_root_in_scan = 0;

    int predicted_root_failures = 0;
    int unexpected_gcd_hits = 0;

    int p_root_hits = 0;
    int second_p_root_hits = 0;
    int second_q_root_hits = 0;

    int gcd_one_count = 0;
    int gcd_p_count = 0;
    int gcd_q_count = 0;
    int gcd_N_count = 0;

    int gcd_other_count = 0;

    long long total_scan_points = 0;

    for (int case_id = 0;
         case_id < static_cast<int>(cases.size());
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 D = N - s * s;

        /*
         * ------------------------------------------------------------
         * 1. Verify the change of variable:
         *
         * x = s+1-k
         *
         * H(x) = Q(k)
         * ------------------------------------------------------------
         */

        const u64 xp = s + 1 - p;

        const u64 Hxp =
            H_polynomial(xp, s, D);

        const u64 Qp =
            q_polynomial(p, s, N);

        if (Hxp != Qp) {
            ++identity_failures;
        }

        /*
         * ------------------------------------------------------------
         * 2. Verify the exact hidden-factor identities:
         *
         * Q(p)=p(4p+3q-7s-4)
         * Q(q)=q(3p+4q-7s-4)
         * ------------------------------------------------------------
         */

        const __uint128_t expected_p =
            static_cast<__uint128_t>(p) *
            static_cast<__uint128_t>(
                4 * p + 3 * q - 7 * s - 4
            );

        const __uint128_t expected_q =
            static_cast<__uint128_t>(q) *
            static_cast<__uint128_t>(
                3 * p + 4 * q - 7 * s - 4
            );

        if (static_cast<__uint128_t>(Qp) != expected_p) {
            ++identity_failures;
        }

        const u64 Qq =
            q_polynomial(q, s, N);

        if (static_cast<__uint128_t>(Qq) != expected_q) {
            ++identity_failures;
        }

        if (Qp % p != 0) {
            ++p_divisibility_failures;
        }

        if (Qq % q != 0) {
            ++q_divisibility_failures;
        }

        /*
         * ------------------------------------------------------------
         * 3. Predicted second modular roots.
         *
         * Q(k)=0 mod p at
         *
         * k=0
         * k=r_p
         *
         * and analogously modulo q.
         * ------------------------------------------------------------
         */

        const u64 rp =
            predicted_second_root(s, p);

        const u64 rq =
            predicted_second_root(s, q);

        const u64 rp_mod =
            rp % p;

        const u64 rq_mod =
            rq % q;

        /*
         * Verify the modular prediction directly.
         */
        if (q_polynomial(rp_mod, s, N) % p != 0) {
            ++predicted_root_failures;
        }

        if (q_polynomial(rq_mod, s, N) % q != 0) {
            ++predicted_root_failures;
        }

        /*
         * p itself is in [1,s] for distinct semiprimes.
         */
        if (p <= s) {
            if (q_polynomial(p, s, N) % p == 0) {
                ++p_root_hits;
            } else {
                ++p_gcd_failures;
            }
        }

        /*
         * Is the second p-root also visible in the search interval?
         */
        if (rp_mod >= 1 && rp_mod <= s) {
            ++second_p_root_hits;

            if (q_polynomial(rp_mod, s, N) % p != 0) {
                ++predicted_root_failures;
            }
        }

        /*
         * q > sqrt(N), so q itself is not in the scan.
         *
         * But its second modular root may have a representative
         * inside [1,s].
         */
        if (rq_mod >= 1 && rq_mod <= s) {
            ++q_root_in_scan;
            ++second_q_root_hits;

            if (q_polynomial(rq_mod, s, N) % q != 0) {
                ++predicted_root_failures;
            }
        }

        /*
         * ------------------------------------------------------------
         * 4. Exhaustive localization scan:
         *
         *     1 <= k <= s
         *
         * and classify gcd(Q(k),N).
         *
         * This is deliberately small:
         * 300 cases x at most 20000 scan points.
         * ------------------------------------------------------------
         */

        for (u64 k = 1; k <= s; ++k) {
            ++total_scan_points;

            const u64 Qk =
                q_polynomial(k, s, N);

            const u64 g =
                std::gcd(Qk, N);

            if (g == 1) {
                ++gcd_one_count;
                continue;
            }

            if (g == p) {
                ++gcd_p_count;

                /*
                 * Any p-hit must satisfy one of the two
                 * modular root equations.
                 */
                const bool root_zero =
                    (k % p) == 0;

                const bool root_second =
                    (k % p) == rp_mod;

                if (!root_zero && !root_second) {
                    ++unexpected_gcd_hits;
                }

                continue;
            }

            if (g == q) {
                ++gcd_q_count;

                const bool root_zero =
                    (k % q) == 0;

                const bool root_second =
                    (k % q) == rq_mod;

                if (!root_zero && !root_second) {
                    ++unexpected_gcd_hits;
                }

                continue;
            }

            if (g == N) {
                ++gcd_N_count;

                /*
                 * For 1 <= k <= s < q, k cannot be 0 mod q.
                 * Therefore this is a CRT coincidence where Q(k)
                 * vanishes modulo both primes.
                 */
                const bool p_root =
                    ((k % p) == 0) ||
                    ((k % p) == rp_mod);

                const bool q_root =
                    ((k % q) == 0) ||
                    ((k % q) == rq_mod);

                if (!p_root || !q_root) {
                    ++unexpected_gcd_hits;
                }

                continue;
            }

            ++gcd_other_count;
            ++unexpected_gcd_hits;
        }

        /*
         * Print a few concrete cases so the modular coordinates
         * can be inspected without producing a huge result file.
         */
        if (case_id < 5) {
            std::cout
                << "\nCASE=" << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << "\n";

            std::cout
                << "r_p=" << rp_mod
                << " r_q=" << rq_mod
                << "\n";

            std::cout
                << "p_in_scan="
                << (p <= s ? "YES" : "NO")
                << "\n";

            std::cout
                << "r_p_in_scan="
                << ((rp_mod >= 1 && rp_mod <= s) ? "YES" : "NO")
                << "\n";

            std::cout
                << "r_q_in_scan="
                << ((rq_mod >= 1 && rq_mod <= s) ? "YES" : "NO")
                << "\n";

            std::cout
                << "Q(p)=" << Qp
                << "\n";

            std::cout
                << "Q(q)=" << Qq
                << "\n";
        }
    }

    /*
     * ------------------------------------------------------------
     * Final results.
     * ------------------------------------------------------------
     */

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
        << "IDENTITY_FAILURES="
        << identity_failures
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
        << "P_GCD_FAILURES="
        << p_gcd_failures
        << "\n";

    std::cout
        << "P_ROOT_HITS="
        << p_root_hits
        << "\n";

    std::cout
        << "SECOND_P_ROOT_HITS="
        << second_p_root_hits
        << "\n";

    std::cout
        << "Q_ROOT_IN_SCAN="
        << q_root_in_scan
        << "\n";

    std::cout
        << "SECOND_Q_ROOT_HITS="
        << second_q_root_hits
        << "\n";

    std::cout
        << "PREDICTED_ROOT_FAILURES="
        << predicted_root_failures
        << "\n";

    std::cout
        << "GCD_1_COUNT="
        << gcd_one_count
        << "\n";

    std::cout
        << "GCD_P_COUNT="
        << gcd_p_count
        << "\n";

    std::cout
        << "GCD_Q_COUNT="
        << gcd_q_count
        << "\n";

    std::cout
        << "GCD_N_COUNT="
        << gcd_N_count
        << "\n";

    std::cout
        << "GCD_OTHER_COUNT="
        << gcd_other_count
        << "\n";

    std::cout
        << "UNEXPECTED_GCD_HITS="
        << unexpected_gcd_hits
        << "\n";

    const bool algebra_ok =
        identity_failures == 0 &&
        p_divisibility_failures == 0 &&
        q_divisibility_failures == 0 &&
        p_gcd_failures == 0 &&
        predicted_root_failures == 0 &&
        unexpected_gcd_hits == 0;

    std::cout
        << "ALGEBRA_STATUS="
        << (algebra_ok ? "PASS" : "FAIL")
        << "\n";

    /*
     * The important structural statement is not whether
     * factor recovery succeeds here. It is whether every observed
     * nontrivial gcd is explained by the two predicted roots modulo
     * p and q.
     */
    if (algebra_ok) {
        std::cout
            << "STRUCTURE_STATUS="
            << "ALL_GCD_HITS_EXPLAINED_BY_MODULAR_ROOTS"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS="
            << "UNEXPLAINED_BEHAVIOR_FOUND"
            << "\n";
    }

    std::cout << "FINISHED EXPERIMENT 350\n";

    return 0;
}
