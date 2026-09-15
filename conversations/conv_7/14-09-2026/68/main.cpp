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
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

struct CoeffSet {
    i64 a;
    i64 b;
    i64 c;
    i64 d;
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

/*
 * Generic polynomial:
 *
 * P(x) =
 *     cN + (x-(s+1))(ax+bs+d)
 */
static i128 evaluate_family(
    const CoeffSet& coeff,
    u64 x,
    const CaseData& c
) {
    const i128 xx = static_cast<i128>(x);
    const i128 s = static_cast<i128>(c.s);
    const i128 N = static_cast<i128>(c.N);

    const i128 R =
        static_cast<i128>(coeff.a) * xx
        +
        static_cast<i128>(coeff.b) * s
        +
        static_cast<i128>(coeff.d);

    return
        static_cast<i128>(coeff.c) * N
        +
        (xx - (s + 1)) * R;
}

/*
 * Exact quotient after substituting x=s+1-p:
 *
 * P(x_p)/p =
 * ap + cq - (a+b)s - a-d
 */
static i128 quotient_formula(
    const CoeffSet& coeff,
    const CaseData& c
) {
    return
        static_cast<i128>(coeff.a) *
            static_cast<i128>(c.p)
        +
        static_cast<i128>(coeff.c) *
            static_cast<i128>(c.q)
        -
        static_cast<i128>(
            coeff.a + coeff.b
        ) *
            static_cast<i128>(c.s)
        -
        static_cast<i128>(
            coeff.a + coeff.d
        );
}

static u64 abs_i128_mod(
    i128 value,
    u64 mod
) {
    if (value < 0) {
        value = -value;
    }

    return static_cast<u64>(
        static_cast<u128>(value) %
        static_cast<u128>(mod)
    );
}

static u64 abs_i128_to_u64(
    i128 value
) {
    if (value < 0) {
        value = -value;
    }

    return static_cast<u64>(value);
}

static void print_i128(i128 value) {
    if (value == 0) {
        std::cout << "0";
        return;
    }

    if (value < 0) {
        std::cout << '-';
        value = -value;
    }

    std::string out;

    while (value > 0) {
        const unsigned digit =
            static_cast<unsigned>(value % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        value /= 10;
    }

    std::reverse(
        out.begin(),
        out.end()
    );

    std::cout << out;
}

static u64 gcd_quotient_N(
    i128 quotient,
    const CaseData& c
) {
    return std::gcd(
        abs_i128_mod(quotient, c.N),
        c.N
    );
}

static bool is_original_H(
    const CoeffSet& c
) {
    return
        c.a == 4 &&
        c.b == 3 &&
        c.c == 3 &&
        c.d == 0;
}

int main() {
    std::cout << "START EXPERIMENT 357\n";

    constexpr int CASE_COUNT = 500;
    constexpr int PRIME_LIMIT = 50000;

    constexpr int COEFF_MIN = -5;
    constexpr int COEFF_MAX = 5;

    std::mt19937_64 rng(
        0x357357357ULL
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

    /*
     * Enumerate all small coefficient combinations.
     *
     * There are 11^4 = 14641 combinations.
     */
    std::vector<CoeffSet> coefficient_sets;

    for (i64 a = COEFF_MIN;
         a <= COEFF_MAX;
         ++a) {

        for (i64 b = COEFF_MIN;
             b <= COEFF_MAX;
             ++b) {

            for (i64 c = COEFF_MIN;
                 c <= COEFF_MAX;
                 ++c) {

                for (i64 d = COEFF_MIN;
                     d <= COEFF_MAX;
                     ++d) {

                    if (a == 0 &&
                        b == 0 &&
                        c == 0 &&
                        d == 0) {
                        continue;
                    }

                    coefficient_sets.push_back({
                        a, b, c, d
                    });
                }
            }
        }
    }

    /*
     * Statistics specifically for H.
     */
    u64 H_identity_failures = 0;
    u64 H_exact_p_gcd = 0;
    u64 H_exact_q_gcd = 0;

    u64 H_quotient_positive = 0;
    u64 H_quotient_zero = 0;
    u64 H_quotient_negative = 0;

    u64 H_abs_quotient_lt_s = 0;
    u64 H_abs_quotient_lt_p = 0;
    u64 H_abs_quotient_lt_q = 0;
    u64 H_abs_quotient_lt_N = 0;

    u64 H_quotient_gcd_N_one = 0;
    u64 H_quotient_gcd_N_p = 0;
    u64 H_quotient_gcd_N_q = 0;
    u64 H_quotient_gcd_N_other = 0;

    /*
     * Count how many coefficient families share properties
     * with H across the entire dataset.
     */
    u64 family_divisibility_failures = 0;

    u64 family_exact_factor_all_cases = 0;

    u64 total_family_tests = 0;

    /*
     * Rather than storing all families' entire statistics,
     * test several exact properties:
     *
     * 1. quotient always nonzero;
     * 2. gcd(quotient,N)=1 for all cases;
     * 3. gcd(quotient,N)=p for all cases;
     * 4. gcd(quotient,N)=q for all cases.
     */

    u64 families_all_gcd_one = 0;
    u64 families_all_gcd_p = 0;
    u64 families_all_gcd_q = 0;

    u64 smallest_abs_quotient_sum =
        UINT64_MAX;

    CoeffSet smallest_family{
        0, 0, 0, 0
    };

    for (const CaseData& c : cases) {
        const CoeffSet H{
            4, 3, 3, 0
        };

        const u64 xp =
            c.s + 1 - c.p;

        const i128 P =
            evaluate_family(
                H,
                xp,
                c
            );

        const i128 Q =
            quotient_formula(
                H,
                c
            );

        /*
         * P(x_p) must equal p*Q.
         */
        const i128 expected =
            static_cast<i128>(c.p) * Q;

        if (P != expected) {
            ++H_identity_failures;
        }

        const u64 gp =
            std::gcd(
                abs_i128_mod(P, c.N),
                c.N
            );

        if (gp == c.p) {
            ++H_exact_p_gcd;
        } else if (gp == c.q) {
            ++H_exact_q_gcd;
        }

        if (Q > 0) {
            ++H_quotient_positive;
        } else if (Q == 0) {
            ++H_quotient_zero;
        } else {
            ++H_quotient_negative;
        }

        const u64 absQ =
            abs_i128_to_u64(Q);

        if (absQ < c.s) {
            ++H_abs_quotient_lt_s;
        }

        if (absQ < c.p) {
            ++H_abs_quotient_lt_p;
        }

        if (absQ < c.q) {
            ++H_abs_quotient_lt_q;
        }

        if (absQ < c.N) {
            ++H_abs_quotient_lt_N;
        }

        const u64 gQ =
            std::gcd(
                absQ,
                c.N
            );

        if (gQ == 1) {
            ++H_quotient_gcd_N_one;
        } else if (gQ == c.p) {
            ++H_quotient_gcd_N_p;
        } else if (gQ == c.q) {
            ++H_quotient_gcd_N_q;
        } else {
            ++H_quotient_gcd_N_other;
        }
    }

    /*
     * Test all coefficient families.
     */
    for (const CoeffSet& coeff :
         coefficient_sets) {

        bool all_divisibility = true;

        bool all_gcd_one = true;
        bool all_gcd_p = true;
        bool all_gcd_q = true;

        for (const CaseData& c : cases) {
            ++total_family_tests;

            const u64 xp =
                c.s + 1 - c.p;

            const i128 value =
                evaluate_family(
                    coeff,
                    xp,
                    c
                );

            const u64 remainder =
                abs_i128_mod(
                    value,
                    c.p
                );

            if (remainder != 0) {
                all_divisibility = false;
            }

            const u64 g =
                std::gcd(
                    abs_i128_mod(
                        value,
                        c.N
                    ),
                    c.N
                );

            if (g != 1) {
                all_gcd_one = false;
            }

            if (g != c.p) {
                all_gcd_p = false;
            }

            if (g != c.q) {
                all_gcd_q = false;
            }
        }

        if (!all_divisibility) {
            ++family_divisibility_failures;
        }

        if (all_gcd_one) {
            ++families_all_gcd_one;
        }

        if (all_gcd_p) {
            ++families_all_gcd_p;
        }

        if (all_gcd_q) {
            ++families_all_gcd_q;
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "COEFFICIENT_FAMILIES="
        << coefficient_sets.size()
        << "\n";

    std::cout
        << "TOTAL_FAMILY_TESTS="
        << total_family_tests
        << "\n";

    std::cout
        << "H_IDENTITY_FAILURES="
        << H_identity_failures
        << "\n";

    std::cout
        << "H_GCD_P_EXACT="
        << H_exact_p_gcd
        << "\n";

    std::cout
        << "H_GCD_Q_EXACT="
        << H_exact_q_gcd
        << "\n";

    std::cout
        << "H_QUOTIENT_POSITIVE="
        << H_quotient_positive
        << "\n";

    std::cout
        << "H_QUOTIENT_ZERO="
        << H_quotient_zero
        << "\n";

    std::cout
        << "H_QUOTIENT_NEGATIVE="
        << H_quotient_negative
        << "\n";

    std::cout
        << "H_ABS_QUOTIENT_LT_S="
        << H_abs_quotient_lt_s
        << "\n";

    std::cout
        << "H_ABS_QUOTIENT_LT_P="
        << H_abs_quotient_lt_p
        << "\n";

    std::cout
        << "H_ABS_QUOTIENT_LT_Q="
        << H_abs_quotient_lt_q
        << "\n";

    std::cout
        << "H_ABS_QUOTIENT_LT_N="
        << H_abs_quotient_lt_N
        << "\n";

    std::cout
        << "H_QUOTIENT_GCD_N_1="
        << H_quotient_gcd_N_one
        << "\n";

    std::cout
        << "H_QUOTIENT_GCD_N_P="
        << H_quotient_gcd_N_p
        << "\n";

    std::cout
        << "H_QUOTIENT_GCD_N_Q="
        << H_quotient_gcd_N_q
        << "\n";

    std::cout
        << "H_QUOTIENT_GCD_N_OTHER="
        << H_quotient_gcd_N_other
        << "\n";

    std::cout
        << "FAMILIES_WITH_DIVISIBILITY_FAILURE="
        << family_divisibility_failures
        << "\n";

    std::cout
        << "FAMILIES_ALL_GCD_1="
        << families_all_gcd_one
        << "\n";

    std::cout
        << "FAMILIES_ALL_GCD_P="
        << families_all_gcd_p
        << "\n";

    std::cout
        << "FAMILIES_ALL_GCD_Q="
        << families_all_gcd_q
        << "\n";

    /*
     * Directly print the symbolic form for H.
     */
    std::cout
        << "\nH_COEFFICIENTS="
        << "(a=4,b=3,c=3,d=0)"
        << "\n";

    std::cout
        << "H_QUOTIENT_FORM="
        << "4p+3q-7s-4"
        << "\n";

    if (H_identity_failures == 0 &&
        family_divisibility_failures == 0) {

        std::cout
            << "STRUCTURE_STATUS="
            << "GENERIC_QUOTIENT_FORM_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 357\n";

    return 0;
}
