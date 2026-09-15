#include <algorithm>
#include <cstdint>
#include <iostream>
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

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 100000) {
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

static u64 second_root(
    u64 s,
    u64 prime
) {
    /*
     * 4r = 7s+4 (mod prime).
     */
    const u64 numerator =
        (
            mod_mul(
                7 % prime,
                s % prime,
                prime
            )
            +
            4
        ) % prime;

    const u64 inv4 =
        mod_pow(
            4,
            prime - 2,
            prime
        );

    return mod_mul(
        numerator,
        inv4,
        prime
    );
}

/*
 * Signed quotient:
 *
 *     t = (7s+4 - 4r) / prime
 *
 * IMPORTANT:
 * t can be negative.
 */
static i64 compute_t(
    u64 s,
    u64 r,
    u64 prime
) {
    const i128 numerator =
        static_cast<i128>(7) *
        static_cast<i128>(s)
        +
        static_cast<i128>(4)
        -
        static_cast<i128>(4) *
        static_cast<i128>(r);

    const i128 denominator =
        static_cast<i128>(prime);

    const i128 quotient =
        numerator / denominator;

    const i128 remainder =
        numerator % denominator;

    if (remainder != 0) {
        std::cerr
            << "INTERNAL ERROR: non-integral t\n";
        std::exit(1);
    }

    return static_cast<i64>(quotient);
}

static i128 exact_root_expression(
    i64 t,
    u64 factor
) {
    return
        static_cast<i128>(t) *
        static_cast<i128>(factor);
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

    std::string result;

    while (value > 0) {
        const unsigned digit =
            static_cast<unsigned>(value % 10);

        result.push_back(
            static_cast<char>('0' + digit)
        );

        value /= 10;
    }

    std::reverse(
        result.begin(),
        result.end()
    );

    std::cout << result;
}

static u128 abs_i128(i128 value) {
    if (value < 0) {
        value = -value;
    }

    return static_cast<u128>(value);
}

int main() {
    std::cout << "START EXPERIMENT 355\n";

    constexpr int CASE_COUNT = 2000;
    constexpr int PRIME_LIMIT = 100000;

    std::mt19937_64 rng(
        0x355355355ULL
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

    u64 equation_failures = 0;
    u64 product_identity_failures = 0;

    i64 min_tp = INT64_MAX;
    i64 max_tp = INT64_MIN;

    i64 min_tq = INT64_MAX;
    i64 max_tq = INT64_MIN;

    i64 min_product = INT64_MAX;
    i64 max_product = INT64_MIN;

    u64 negative_tp = 0;
    u64 zero_tp = 0;
    u64 positive_tp = 0;

    u64 negative_tq = 0;
    u64 zero_tq = 0;
    u64 positive_tq = 0;

    u64 product_negative = 0;
    u64 product_zero = 0;
    u64 product_positive = 0;

    u64 product_abs_lt_10 = 0;
    u64 product_abs_lt_20 = 0;
    u64 product_abs_lt_50 = 0;
    u64 product_abs_lt_100 = 0;

    u64 tp_abs_lt_10 = 0;
    u64 tq_abs_lt_10 = 0;

    u64 same_sign = 0;
    u64 opposite_sign = 0;

    u64 tp_eq_tq = 0;
    u64 tp_eq_minus_tq = 0;

    u64 product_eq_16 = 0;
    u64 product_eq_24 = 0;
    u64 product_eq_25 = 0;
    u64 product_eq_32 = 0;
    u64 product_eq_36 = 0;
    u64 product_eq_42 = 0;
    u64 product_eq_48 = 0;
    u64 product_eq_49 = 0;

    std::set<i64> distinct_tp;
    std::set<i64> distinct_tq;
    std::set<i64> distinct_products;

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

        const u64 rp =
            second_root(
                s,
                p
            );

        const u64 rq =
            second_root(
                s,
                q
            );

        const i64 tp =
            compute_t(
                s,
                rp,
                p
            );

        const i64 tq =
            compute_t(
                s,
                rq,
                q
            );

        /*
         * Exact identities:
         *
         * 4rp + tp*p = A
         * 4rq + tq*q = A
         */
        const i128 lhs_p =
            static_cast<i128>(4) *
            static_cast<i128>(rp)
            +
            exact_root_expression(
                tp,
                p
            );

        const i128 lhs_q =
            static_cast<i128>(4) *
            static_cast<i128>(rq)
            +
            exact_root_expression(
                tq,
                q
            );

        if (lhs_p != static_cast<i128>(A) ||
            lhs_q != static_cast<i128>(A)) {

            ++equation_failures;
        }

        /*
         * Signed product.
         */
        const i128 product128 =
            static_cast<i128>(tp) *
            static_cast<i128>(tq);

        const i64 product =
            static_cast<i64>(product128);

        distinct_tp.insert(tp);
        distinct_tq.insert(tq);
        distinct_products.insert(product);

        min_tp = std::min(min_tp, tp);
        max_tp = std::max(max_tp, tp);

        min_tq = std::min(min_tq, tq);
        max_tq = std::max(max_tq, tq);

        min_product =
            std::min(min_product, product);

        max_product =
            std::max(max_product, product);

        if (tp < 0) {
            ++negative_tp;
        } else if (tp == 0) {
            ++zero_tp;
        } else {
            ++positive_tp;
        }

        if (tq < 0) {
            ++negative_tq;
        } else if (tq == 0) {
            ++zero_tq;
        } else {
            ++positive_tq;
        }

        if (product < 0) {
            ++product_negative;
        } else if (product == 0) {
            ++product_zero;
        } else {
            ++product_positive;
        }

        const u64 abs_product =
            product < 0
                ? static_cast<u64>(-product)
                : static_cast<u64>(product);

        if (abs_product < 10) {
            ++product_abs_lt_10;
        }

        if (abs_product < 20) {
            ++product_abs_lt_20;
        }

        if (abs_product < 50) {
            ++product_abs_lt_50;
        }

        if (abs_product < 100) {
            ++product_abs_lt_100;
        }

        if (
            (tp >= -9 && tp <= 9)
        ) {
            ++tp_abs_lt_10;
        }

        if (
            (tq >= -9 && tq <= 9)
        ) {
            ++tq_abs_lt_10;
        }

        if (
            (tp < 0 && tq < 0) ||
            (tp > 0 && tq > 0)
        ) {
            ++same_sign;
        }

        if (
            (tp < 0 && tq > 0) ||
            (tp > 0 && tq < 0)
        ) {
            ++opposite_sign;
        }

        if (tp == tq) {
            ++tp_eq_tq;
        }

        if (tp == -tq) {
            ++tp_eq_minus_tq;
        }

        if (product == 16) ++product_eq_16;
        if (product == 24) ++product_eq_24;
        if (product == 25) ++product_eq_25;
        if (product == 32) ++product_eq_32;
        if (product == 36) ++product_eq_36;
        if (product == 42) ++product_eq_42;
        if (product == 48) ++product_eq_48;
        if (product == 49) ++product_eq_49;

        /*
         * ---------------------------------------------------------
         * Exact product identity:
         *
         * tp*tq*N
         *
         * =
         *
         * (A-4rp)(A-4rq)
         *
         * ---------------------------------------------------------
         */
        const i128 left =
            product128 *
            static_cast<i128>(N);

        const i128 right =
            (
                static_cast<i128>(A)
                -
                static_cast<i128>(4) *
                static_cast<i128>(rp)
            )
            *
            (
                static_cast<i128>(A)
                -
                static_cast<i128>(4) *
                static_cast<i128>(rq)
            );

        if (left != right) {
            ++product_identity_failures;
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
                << "t_p="
                << tp
                << "\n";

            std::cout
                << "t_q="
                << tq
                << "\n";

            std::cout
                << "T=t_p*t_q="
                << product
                << "\n";

            std::cout
                << "4r_p+t_p*p="
                << static_cast<u64>(lhs_p)
                << "\n";

            std::cout
                << "4r_q+t_q*q="
                << static_cast<u64>(lhs_q)
                << "\n";

            std::cout
                << "7s+4="
                << A
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
        << "PRODUCT_IDENTITY_FAILURES="
        << product_identity_failures
        << "\n";

    std::cout
        << "MIN_TP="
        << min_tp
        << "\n";

    std::cout
        << "MAX_TP="
        << max_tp
        << "\n";

    std::cout
        << "MIN_TQ="
        << min_tq
        << "\n";

    std::cout
        << "MAX_TQ="
        << max_tq
        << "\n";

    std::cout
        << "MIN_TP_TQ="
        << min_product
        << "\n";

    std::cout
        << "MAX_TP_TQ="
        << max_product
        << "\n";

    std::cout
        << "DISTINCT_TP="
        << distinct_tp.size()
        << "\n";

    std::cout
        << "DISTINCT_TQ="
        << distinct_tq.size()
        << "\n";

    std::cout
        << "DISTINCT_TP_TQ="
        << distinct_products.size()
        << "\n";

    std::cout
        << "TP_NEGATIVE="
        << negative_tp
        << "\n";

    std::cout
        << "TP_ZERO="
        << zero_tp
        << "\n";

    std::cout
        << "TP_POSITIVE="
        << positive_tp
        << "\n";

    std::cout
        << "TQ_NEGATIVE="
        << negative_tq
        << "\n";

    std::cout
        << "TQ_ZERO="
        << zero_tq
        << "\n";

    std::cout
        << "TQ_POSITIVE="
        << positive_tq
        << "\n";

    std::cout
        << "TP_TQ_NEGATIVE="
        << product_negative
        << "\n";

    std::cout
        << "TP_TQ_ZERO="
        << product_zero
        << "\n";

    std::cout
        << "TP_TQ_POSITIVE="
        << product_positive
        << "\n";

    std::cout
        << "ABS_PRODUCT_LT_10="
        << product_abs_lt_10
        << "\n";

    std::cout
        << "ABS_PRODUCT_LT_20="
        << product_abs_lt_20
        << "\n";

    std::cout
        << "ABS_PRODUCT_LT_50="
        << product_abs_lt_50
        << "\n";

    std::cout
        << "ABS_PRODUCT_LT_100="
        << product_abs_lt_100
        << "\n";

    std::cout
        << "ABS_TP_LT_10="
        << tp_abs_lt_10
        << "\n";

    std::cout
        << "ABS_TQ_LT_10="
        << tq_abs_lt_10
        << "\n";

    std::cout
        << "TP_TQ_SAME_SIGN="
        << same_sign
        << "\n";

    std::cout
        << "TP_TQ_OPPOSITE_SIGN="
        << opposite_sign
        << "\n";

    std::cout
        << "TP_EQ_TQ="
        << tp_eq_tq
        << "\n";

    std::cout
        << "TP_EQ_MINUS_TQ="
        << tp_eq_minus_tq
        << "\n";

    std::cout
        << "PRODUCT_EQ_16="
        << product_eq_16
        << "\n";

    std::cout
        << "PRODUCT_EQ_24="
        << product_eq_24
        << "\n";

    std::cout
        << "PRODUCT_EQ_25="
        << product_eq_25
        << "\n";

    std::cout
        << "PRODUCT_EQ_32="
        << product_eq_32
        << "\n";

    std::cout
        << "PRODUCT_EQ_36="
        << product_eq_36
        << "\n";

    std::cout
        << "PRODUCT_EQ_42="
        << product_eq_42
        << "\n";

    std::cout
        << "PRODUCT_EQ_48="
        << product_eq_48
        << "\n";

    std::cout
        << "PRODUCT_EQ_49="
        << product_eq_49
        << "\n";

    if (equation_failures == 0 &&
        product_identity_failures == 0) {

        std::cout
            << "STRUCTURE_STATUS="
            << "SIGNED_TWO_ROOT_RELATION_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 355\n";

    return 0;
}