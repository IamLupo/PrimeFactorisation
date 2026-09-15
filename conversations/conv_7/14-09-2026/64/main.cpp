#include <algorithm>
#include <cstdint>
#include <iostream>
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

static std::vector<int> usable_primes(
    const std::vector<int>& primes,
    int min_p,
    int max_p
) {
    std::vector<int> result;

    for (int p : primes) {
        if (p >= min_p && p <= max_p) {
            result.push_back(p);
        }
    }

    return result;
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::uniform_int_distribution<std::size_t> dist(
        0,
        primes.size() - 1
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
                std::min(primes[i], primes[j])
            );

        const u64 q =
            static_cast<u64>(
                std::max(primes[i], primes[j])
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

static u64 mod_mul(
    u64 a,
    u64 b,
    u64 mod
) {
    return static_cast<u64>(
        (
            static_cast<u128>(a) *
            static_cast<u128>(b)
        ) %
        static_cast<u128>(mod)
    );
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
            result = mod_mul(
                result,
                a,
                mod
            );
        }

        a = mod_mul(
            a,
            a,
            mod
        );

        e >>= 1ULL;
    }

    return result;
}

static u64 inverse_mod_prime(
    u64 a,
    u64 p
) {
    /*
     * Fermat:
     *
     * a^(p-2) = a^(-1) mod p
     */
    return mod_pow(
        a % p,
        p - 2,
        p
    );
}

static u64 second_root(
    u64 s,
    u64 p
) {
    /*
     * 4r = 7s+4 (mod p)
     *
     * r = (7s+4) * 4^(-1) (mod p)
     */
    const u64 numerator =
        (
            mod_mul(
                7 % p,
                s % p,
                p
            )
            +
            4
        ) % p;

    const u64 inv4 =
        inverse_mod_prime(
            4,
            p
        );

    return mod_mul(
        numerator,
        inv4,
        p
    );
}

static u64 compute_t(
    u64 p,
    u64 s,
    u64 r
) {
    /*
     * 4r + tp = 7s+4
     *
     * Since r<p, the numerator is positive
     * for all generated cases.
     */
    const u64 numerator =
        7 * s + 4 - 4 * r;

    return numerator / p;
}

static u64 ceil_div(
    u64 a,
    u64 b
) {
    return a / b + (a % b != 0);
}

static bool prime_residue_allowed(
    int p,
    u64 t,
    u64 A
) {
    /*
     * From:
     *
     *     tp = A (mod 4)
     *
     * with p an odd prime.
     */
    return
        (
            (
                t *
                static_cast<u64>(p % 4)
            ) % 4
        )
        ==
        (A % 4);
}

static void print_u128(
    u128 x
) {
    if (x == 0) {
        std::cout << "0";
        return;
    }

    std::string out;

    while (x > 0) {
        const unsigned digit =
            static_cast<unsigned>(x % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(
        out.begin(),
        out.end()
    );

    std::cout << out;
}

int main() {
    std::cout << "START EXPERIMENT 353\n";

    constexpr int CASE_COUNT = 1000;
    constexpr int PRIME_LIMIT = 100000;

    constexpr int MIN_FACTOR = 1009;
    constexpr int MAX_FACTOR = 100000;

    std::mt19937_64 rng(
        0x353353353ULL
    );

    const std::vector<int> all_primes =
        sieve_primes(
            PRIME_LIMIT
        );

    const std::vector<int> primes =
        usable_primes(
            all_primes,
            MIN_FACTOR,
            MAX_FACTOR
        );

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    u64 total_integer_candidates = 0;
    u64 total_prime_candidates = 0;
    u64 total_congruence_candidates = 0;

    u64 min_integer_candidates = UINT64_MAX;
    u64 max_integer_candidates = 0;

    u64 min_prime_candidates = UINT64_MAX;
    u64 max_prime_candidates = 0;

    u64 min_congruence_candidates = UINT64_MAX;
    u64 max_congruence_candidates = 0;

    u64 factor_not_in_interval = 0;
    u64 factor_not_prime_candidate = 0;
    u64 congruence_failures = 0;

    u64 interval_width_lt_10 = 0;
    u64 interval_width_lt_100 = 0;
    u64 interval_width_lt_1000 = 0;
    u64 interval_width_lt_10000 = 0;

    u64 prime_candidates_lt_10 = 0;
    u64 prime_candidates_lt_100 = 0;
    u64 prime_candidates_lt_1000 = 0;
    u64 prime_candidates_lt_10000 = 0;

    u64 congruence_candidates_lt_10 = 0;
    u64 congruence_candidates_lt_100 = 0;
    u64 congruence_candidates_lt_1000 = 0;
    u64 congruence_candidates_lt_10000 = 0;

    u64 congruence_exactly_1 = 0;
    u64 congruence_exactly_2 = 0;
    u64 congruence_exactly_3 = 0;
    u64 congruence_exactly_4 = 0;
    u64 congruence_exactly_5 = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 A =
            7 * s + 4;

        const u64 r =
            second_root(
                s,
                p
            );

        const u64 t =
            compute_t(
                p,
                s,
                r
            );

        /*
         * From
         *
         *     4r + tp = 7s+4
         *
         * and 1 <= r <= p-1:
         *
         *     4 + tp <= 7s+4
         *
         *     7s+4 < 4p+tp
         *
         * hence
         *
         *     (7s+8)/(t+4) <= p
         *
         *     p < (7s+4)/t.
         */
        const u64 lower =
            ceil_div(
                A + 4,
                t + 4
            );

        const u64 upper =
            (A - 1) / t;

        if (lower > upper ||
            p < lower ||
            p > upper) {

            ++factor_not_in_interval;
        }

        const u64 integer_candidates =
            lower <= upper
                ? upper - lower + 1
                : 0;

        /*
         * Prime candidates.
         */
        auto first_it =
            std::lower_bound(
                primes.begin(),
                primes.end(),
                static_cast<int>(lower)
            );

        auto last_it =
            std::upper_bound(
                primes.begin(),
                primes.end(),
                static_cast<int>(upper)
            );

        const u64 prime_candidates =
            static_cast<u64>(
                last_it - first_it
            );

        /*
         * Apply:
         *
         *     tp = 7s+4 (mod 4)
         */
        u64 congruence_candidates = 0;

        for (auto it = first_it;
             it != last_it;
             ++it) {

            if (prime_residue_allowed(
                    *it,
                    t,
                    A
                )) {

                ++congruence_candidates;
            }
        }

        if (!prime_residue_allowed(
                static_cast<int>(p),
                t,
                A
            )) {

            ++congruence_failures;
        }

        /*
         * Check that p is actually in the prime list.
         */
        bool found_p = false;

        for (auto it = first_it;
             it != last_it;
             ++it) {

            if (static_cast<u64>(*it) == p) {
                found_p = true;
                break;
            }
        }

        if (!found_p) {
            ++factor_not_prime_candidate;
        }

        total_integer_candidates +=
            integer_candidates;

        total_prime_candidates +=
            prime_candidates;

        total_congruence_candidates +=
            congruence_candidates;

        min_integer_candidates =
            std::min(
                min_integer_candidates,
                integer_candidates
            );

        max_integer_candidates =
            std::max(
                max_integer_candidates,
                integer_candidates
            );

        min_prime_candidates =
            std::min(
                min_prime_candidates,
                prime_candidates
            );

        max_prime_candidates =
            std::max(
                max_prime_candidates,
                prime_candidates
            );

        min_congruence_candidates =
            std::min(
                min_congruence_candidates,
                congruence_candidates
            );

        max_congruence_candidates =
            std::max(
                max_congruence_candidates,
                congruence_candidates
            );

        if (integer_candidates < 10) {
            ++interval_width_lt_10;
        }

        if (integer_candidates < 100) {
            ++interval_width_lt_100;
        }

        if (integer_candidates < 1000) {
            ++interval_width_lt_1000;
        }

        if (integer_candidates < 10000) {
            ++interval_width_lt_10000;
        }

        if (prime_candidates < 10) {
            ++prime_candidates_lt_10;
        }

        if (prime_candidates < 100) {
            ++prime_candidates_lt_100;
        }

        if (prime_candidates < 1000) {
            ++prime_candidates_lt_1000;
        }

        if (prime_candidates < 10000) {
            ++prime_candidates_lt_10000;
        }

        if (congruence_candidates < 10) {
            ++congruence_candidates_lt_10;
        }

        if (congruence_candidates < 100) {
            ++congruence_candidates_lt_100;
        }

        if (congruence_candidates < 1000) {
            ++congruence_candidates_lt_1000;
        }

        if (congruence_candidates < 10000) {
            ++congruence_candidates_lt_10000;
        }

        if (congruence_candidates == 1) {
            ++congruence_exactly_1;
        }

        if (congruence_candidates == 2) {
            ++congruence_exactly_2;
        }

        if (congruence_candidates == 3) {
            ++congruence_exactly_3;
        }

        if (congruence_candidates == 4) {
            ++congruence_exactly_4;
        }

        if (congruence_candidates == 5) {
            ++congruence_exactly_5;
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
                << "r="
                << r
                << "\n";

            std::cout
                << "t="
                << t
                << "\n";

            std::cout
                << "LOWER="
                << lower
                << "\n";

            std::cout
                << "UPPER="
                << upper
                << "\n";

            std::cout
                << "INTEGER_CANDIDATES="
                << integer_candidates
                << "\n";

            std::cout
                << "PRIME_CANDIDATES="
                << prime_candidates
                << "\n";

            std::cout
                << "CONGRUENCE_CANDIDATES="
                << congruence_candidates
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "TOTAL_INTEGER_CANDIDATES="
        << total_integer_candidates
        << "\n";

    std::cout
        << "TOTAL_PRIME_CANDIDATES="
        << total_prime_candidates
        << "\n";

    std::cout
        << "TOTAL_CONGRUENCE_CANDIDATES="
        << total_congruence_candidates
        << "\n";

    std::cout
        << "AVG_INTEGER_CANDIDATES="
        << total_integer_candidates /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "AVG_PRIME_CANDIDATES="
        << total_prime_candidates /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "AVG_CONGRUENCE_CANDIDATES="
        << total_congruence_candidates /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "MIN_INTEGER_CANDIDATES="
        << min_integer_candidates
        << "\n";

    std::cout
        << "MAX_INTEGER_CANDIDATES="
        << max_integer_candidates
        << "\n";

    std::cout
        << "MIN_PRIME_CANDIDATES="
        << min_prime_candidates
        << "\n";

    std::cout
        << "MAX_PRIME_CANDIDATES="
        << max_prime_candidates
        << "\n";

    std::cout
        << "MIN_CONGRUENCE_CANDIDATES="
        << min_congruence_candidates
        << "\n";

    std::cout
        << "MAX_CONGRUENCE_CANDIDATES="
        << max_congruence_candidates
        << "\n";

    std::cout
        << "INTERVAL_WIDTH_LT_10="
        << interval_width_lt_10
        << "\n";

    std::cout
        << "INTERVAL_WIDTH_LT_100="
        << interval_width_lt_100
        << "\n";

    std::cout
        << "INTERVAL_WIDTH_LT_1000="
        << interval_width_lt_1000
        << "\n";

    std::cout
        << "INTERVAL_WIDTH_LT_10000="
        << interval_width_lt_10000
        << "\n";

    std::cout
        << "PRIME_CANDIDATES_LT_10="
        << prime_candidates_lt_10
        << "\n";

    std::cout
        << "PRIME_CANDIDATES_LT_100="
        << prime_candidates_lt_100
        << "\n";

    std::cout
        << "PRIME_CANDIDATES_LT_1000="
        << prime_candidates_lt_1000
        << "\n";

    std::cout
        << "PRIME_CANDIDATES_LT_10000="
        << prime_candidates_lt_10000
        << "\n";

    std::cout
        << "CONGRUENCE_CANDIDATES_LT_10="
        << congruence_candidates_lt_10
        << "\n";

    std::cout
        << "CONGRUENCE_CANDIDATES_LT_100="
        << congruence_candidates_lt_100
        << "\n";

    std::cout
        << "CONGRUENCE_CANDIDATES_LT_1000="
        << congruence_candidates_lt_1000
        << "\n";

    std::cout
        << "CONGRUENCE_CANDIDATES_LT_10000="
        << congruence_candidates_lt_10000
        << "\n";

    std::cout
        << "CONGRUENCE_EXACTLY_1="
        << congruence_exactly_1
        << "\n";

    std::cout
        << "CONGRUENCE_EXACTLY_2="
        << congruence_exactly_2
        << "\n";

    std::cout
        << "CONGRUENCE_EXACTLY_3="
        << congruence_exactly_3
        << "\n";

    std::cout
        << "CONGRUENCE_EXACTLY_4="
        << congruence_exactly_4
        << "\n";

    std::cout
        << "CONGRUENCE_EXACTLY_5="
        << congruence_exactly_5
        << "\n";

    std::cout
        << "FACTOR_NOT_IN_INTERVAL="
        << factor_not_in_interval
        << "\n";

    std::cout
        << "FACTOR_NOT_PRIME_CANDIDATE="
        << factor_not_prime_candidate
        << "\n";

    std::cout
        << "CONGRUENCE_FAILURES="
        << congruence_failures
        << "\n";

    if (factor_not_in_interval == 0 &&
        factor_not_prime_candidate == 0 &&
        congruence_failures == 0) {

        std::cout
            << "STRUCTURE_STATUS="
            << "FACTOR_ALWAYS_LOCALIZED_BY_T_INTERVAL"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 353\n";

    return 0;
}