#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
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
        (static_cast<u128>(a) *
         static_cast<u128>(b)) %
        static_cast<u128>(mod)
    );
}

static u64 mod_pow(u64 a, u64 e, u64 mod) {
    u64 result = 1 % mod;

    while (e != 0) {
        if (e & 1ULL) {
            result = mod_mul(result, a, result == 0 ? 1 : mod);
        }

        a = mod_mul(a, a, mod);
        e >>= 1ULL;
    }

    return result;
}

static u64 mod_inverse_prime(u64 a, u64 p) {
    /*
     * p is prime and p != 2 here because our generated primes
     * are all >= 1009.
     */
    return mod_pow(a % p, p - 2, p);
}

static u64 second_root(
    u64 s,
    u64 p
) {
    const u64 c =
        (
            mod_mul(7 % p, s % p, p) +
            4
        ) % p;

    return mod_mul(
        c,
        mod_inverse_prime(4, p),
        p
    );
}

static u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
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

static void print_i128(i128 x) {
    if (x == 0) {
        std::cout << "0";
        return;
    }

    if (x < 0) {
        std::cout << "-";
        x = -x;
    }

    std::string out;

    while (x > 0) {
        const int digit =
            static_cast<int>(x % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(out.begin(), out.end());

    std::cout << out;
}

int main() {
    std::cout << "START EXPERIMENT 352\n";

    constexpr int CASE_COUNT = 1000;
    constexpr int PRIME_LIMIT = 50000;

    std::mt19937_64 rng(
        0x352352352ULL
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

    u64 min_t = UINT64_MAX;
    u64 max_t = 0;

    u64 min_abs_t = UINT64_MAX;
    u64 max_abs_t = 0;

    u64 min_difference = UINT64_MAX;
    u64 max_difference = 0;

    u64 second_root_less_than_p = 0;
    u64 second_root_equal_p = 0;
    u64 second_root_greater_than_p = 0;

    u64 t_equal_1 = 0;
    u64 t_equal_2 = 0;
    u64 t_equal_3 = 0;
    u64 t_equal_4 = 0;
    u64 t_equal_5 = 0;
    u64 t_equal_6 = 0;
    u64 t_equal_7 = 0;
    u64 t_equal_8 = 0;
    u64 t_equal_9 = 0;
    u64 t_equal_10 = 0;

    u64 t_less_than_20 = 0;
    u64 t_less_than_50 = 0;
    u64 t_less_than_100 = 0;
    u64 t_less_than_1000 = 0;

    u64 t_zero = 0;

    u64 equation_failures = 0;
    u64 range_failures = 0;
    u64 modular_failures = 0;

    /*
     * Test whether the multiplier t has some simple relationship
     * with the floor/square-gap quantities.
     */
    u64 t_equal_q_minus_p = 0;
    u64 t_equal_q_div_p_floor = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 D =
            N - s * s;

        const u64 r =
            second_root(s, p);

        /*
         * The defining exact integer equation:
         *
         *     4r + t p = 7s+4.
         *
         * Since r is the canonical residue 0<=r<p,
         * t should be a positive integer.
         */
        const i128 numerator =
            static_cast<i128>(7) *
            static_cast<i128>(s)
            +
            static_cast<i128>(4)
            -
            static_cast<i128>(4) *
            static_cast<i128>(r);

        const bool divisible =
            (numerator % static_cast<i128>(p)) == 0;

        if (!divisible) {
            ++equation_failures;
            continue;
        }

        const i128 t128 =
            numerator / static_cast<i128>(p);

        if (t128 < 0 ||
            t128 > static_cast<i128>(UINT64_MAX)) {
            ++range_failures;
            continue;
        }

        const u64 t =
            static_cast<u64>(t128);

        /*
         * Reconstruct the equation exactly.
         */
        const i128 reconstructed =
            static_cast<i128>(4) *
            static_cast<i128>(r)
            +
            static_cast<i128>(t) *
            static_cast<i128>(p);

        if (reconstructed !=
            static_cast<i128>(7) *
            static_cast<i128>(s)
            +
            static_cast<i128>(4)) {

            ++equation_failures;
        }

        /*
         * Check the root directly.
         *
         * 4r-(7s+4) = -tp,
         * therefore it must be 0 modulo p.
         */
        const i128 root_expression =
            static_cast<i128>(4) *
            static_cast<i128>(r)
            -
            static_cast<i128>(7) *
            static_cast<i128>(s)
            -
            static_cast<i128>(4);

        if (root_expression %
            static_cast<i128>(p) != 0) {

            ++modular_failures;
        }

        /*
         * r is always a residue modulo p.
         */
        if (!(r < p)) {
            ++range_failures;
        }

        /*
         * r is empirically always inside the factor-search
         * interval [1,s] in experiment 351.
         */
        if (!(r >= 1 && r <= s)) {
            ++range_failures;
        }

        if (r < p) {
            ++second_root_less_than_p;
        } else if (r == p) {
            ++second_root_equal_p;
        } else {
            ++second_root_greater_than_p;
        }

        /*
         * Statistics on t.
         */
        if (t < min_t) {
            min_t = t;
        }

        if (t > max_t) {
            max_t = t;
        }

        const u64 abs_t = t;

        if (abs_t < min_abs_t) {
            min_abs_t = abs_t;
        }

        if (abs_t > max_abs_t) {
            max_abs_t = abs_t;
        }

        const u64 difference =
            (p >= r) ? (p - r) : (r - p);

        if (difference < min_difference) {
            min_difference = difference;
        }

        if (difference > max_difference) {
            max_difference = difference;
        }

        if (t == 0) {
            ++t_zero;
        }

        if (t == 1) ++t_equal_1;
        if (t == 2) ++t_equal_2;
        if (t == 3) ++t_equal_3;
        if (t == 4) ++t_equal_4;
        if (t == 5) ++t_equal_5;
        if (t == 6) ++t_equal_6;
        if (t == 7) ++t_equal_7;
        if (t == 8) ++t_equal_8;
        if (t == 9) ++t_equal_9;
        if (t == 10) ++t_equal_10;

        if (t < 20) ++t_less_than_20;
        if (t < 50) ++t_less_than_50;
        if (t < 100) ++t_less_than_100;
        if (t < 1000) ++t_less_than_1000;

        /*
         * A few possible simple relationships.
         */
        if (t == q - p) {
            ++t_equal_q_minus_p;
        }

        if (t ==
            static_cast<u64>(
                static_cast<i128>(q) /
                static_cast<i128>(p)
            )) {
            ++t_equal_q_div_p_floor;
        }

        if (case_id < 10) {
            std::cout
                << "\nCASE=" << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << " D=" << D
                << "\n";

            std::cout
                << "r_p=" << r
                << "\n";

            std::cout
                << "p-r=";

            if (p >= r) {
                std::cout << (p - r);
            } else {
                std::cout << "-";
                print_i128(
                    static_cast<i128>(r) -
                    static_cast<i128>(p)
                );
            }

            std::cout << "\n";

            std::cout
                << "t=" << t
                << "\n";

            std::cout
                << "4r+tp=";

            print_i128(reconstructed);

            std::cout
                << "\n7s+4="
                << (7 * s + 4)
                << "\n";

            std::cout
                << "q_over_p_floor="
                << (q / p)
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "EQUATION_FAILURES="
        << equation_failures
        << "\n";

    std::cout
        << "RANGE_FAILURES="
        << range_failures
        << "\n";

    std::cout
        << "MODULAR_FAILURES="
        << modular_failures
        << "\n";

    std::cout
        << "MIN_T="
        << min_t
        << "\n";

    std::cout
        << "MAX_T="
        << max_t
        << "\n";

    std::cout
        << "SECOND_ROOT_LT_P="
        << second_root_less_than_p
        << "\n";

    std::cout
        << "SECOND_ROOT_EQ_P="
        << second_root_equal_p
        << "\n";

    std::cout
        << "SECOND_ROOT_GT_P="
        << second_root_greater_than_p
        << "\n";

    std::cout
        << "T_ZERO="
        << t_zero
        << "\n";

    std::cout
        << "T_EQ_1="
        << t_equal_1
        << "\n";

    std::cout
        << "T_EQ_2="
        << t_equal_2
        << "\n";

    std::cout
        << "T_EQ_3="
        << t_equal_3
        << "\n";

    std::cout
        << "T_EQ_4="
        << t_equal_4
        << "\n";

    std::cout
        << "T_EQ_5="
        << t_equal_5
        << "\n";

    std::cout
        << "T_EQ_6="
        << t_equal_6
        << "\n";

    std::cout
        << "T_EQ_7="
        << t_equal_7
        << "\n";

    std::cout
        << "T_EQ_8="
        << t_equal_8
        << "\n";

    std::cout
        << "T_EQ_9="
        << t_equal_9
        << "\n";

    std::cout
        << "T_EQ_10="
        << t_equal_10
        << "\n";

    std::cout
        << "T_LT_20="
        << t_less_than_20
        << "\n";

    std::cout
        << "T_LT_50="
        << t_less_than_50
        << "\n";

    std::cout
        << "T_LT_100="
        << t_less_than_100
        << "\n";

    std::cout
        << "T_LT_1000="
        << t_less_than_1000
        << "\n";

    std::cout
        << "T_EQ_Q_MINUS_P="
        << t_equal_q_minus_p
        << "\n";

    std::cout
        << "T_EQ_FLOOR_Q_OVER_P="
        << t_equal_q_div_p_floor
        << "\n";

    if (equation_failures == 0 &&
        range_failures == 0 &&
        modular_failures == 0) {

        std::cout
            << "STRUCTURE_STATUS="
            << "SECOND_ROOT_EXACT_INTEGER_RELATION_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS="
            << "FAIL"
            << "\n";
    }

    std::cout << "FINISHED EXPERIMENT 352\n";

    return 0;
}
