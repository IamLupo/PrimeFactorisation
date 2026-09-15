#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
};

struct Polynomial {
    std::vector<u64> coeff;
};

struct Stats {
    u64 monomial = 0;
    u64 non_monomial = 0;
    u64 zero = 0;

    u64 min_t_valuation = UINT64_MAX;
    u64 max_t_valuation = 0;

    u64 support_sum = 0;
    u64 max_support = 0;
};

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool negative = x < 0;

    u128 v =
        negative
            ? static_cast<u128>(-x)
            : static_cast<u128>(x);

    std::string s;

    while (v > 0) {
        int d =
            static_cast<int>(v % 10);

        s.push_back(
            static_cast<char>('0' + d)
        );

        v /= 10;
    }

    if (negative) {
        s.push_back('-');
    }

    std::reverse(
        s.begin(),
        s.end()
    );

    return s;
}

u64 gcd_u64(
    u64 a,
    u64 b
) {
    return std::gcd(a, b);
}

u64 integer_sqrt(
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    u64 lo = 0;

    u64 hi =
        std::min<u64>(
            n,
            1ULL << 32
        );

    while (lo + 1 < hi) {
        u64 mid =
            lo + (hi - lo) / 2;

        if (mid <= n / mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    if (hi <= n / hi) {
        return hi;
    }

    return lo;
}

std::vector<u64> base_digits(
    u64 x,
    u64 base
) {
    std::vector<u64> d;

    if (x == 0) {
        d.push_back(0);
        return d;
    }

    while (x > 0) {
        d.push_back(
            x % base
        );

        x /= base;
    }

    return d;
}

/*
    CORRECT implementation of

        M_m(y)
        =
        #{ x : 0 <= x <= y and x <=_p m }

    Digits are scanned MSB -> LSB.

    If the first differing digit is i, then
    all LOWER digits are free subject to x_j <= m_j.

    Therefore the multiplier is

        prod_{j<i} (m_j + 1)

    and NOT a higher-digit suffix product.
*/
u64 miss_prefix(
    u64 m,
    i128 y_signed,
    u64 base
) {
    if (y_signed < 0) {
        return 0;
    }

    u64 y =
        static_cast<u64>(
            y_signed
        );

    std::vector<u64> md =
        base_digits(
            m,
            base
        );

    std::vector<u64> yd =
        base_digits(
            y,
            base
        );

    std::size_t L =
        std::max(
            md.size(),
            yd.size()
        );

    md.resize(
        L,
        0
    );

    yd.resize(
        L,
        0
    );

    /*
        lower_product[i] =
            prod_{j<i} (m_j + 1)
    */
    std::vector<u128> lower_product(
        L,
        1
    );

    u128 running = 1;

    for (
        std::size_t i = 0;
        i < L;
        ++i
    ) {
        lower_product[i] =
            running;

        running *=
            static_cast<u128>(
                md[i] + 1
            );
    }

    u128 answer = 0;

    bool tight = true;

    for (
        std::size_t pos = L;
        pos-- > 0;
    ) {
        if (!tight) {
            break;
        }

        u64 yi = yd[pos];
        u64 mi = md[pos];

        /*
            Digits x_i < y_i which are also <= m_i.
        */
        u64 choices =
            std::min(
                yi,
                mi + 1
            );

        answer +=
            static_cast<u128>(
                choices
            ) *
            lower_product[pos];

        /*
            Continue with x_i = y_i only if allowed.
        */
        if (yi <= mi) {
            tight = true;
        } else {
            tight = false;
            break;
        }
    }

    /*
        Every digit matched y, hence y itself is admissible.
    */
    if (tight) {
        ++answer;
    }

    return static_cast<u64>(
        answer
    );
}

/*
    Independent brute-force definition.
*/
u64 brute_miss_prefix(
    u64 m,
    u64 y,
    u64 base
) {
    u64 count = 0;

    for (
        u64 x = 0;
        x <= y;
        ++x
    ) {
        u64 a = x;
        u64 b = m;

        bool allowed = true;

        while (
            a != 0 ||
            b != 0
        ) {
            u64 xd = a % base;
            u64 md = b % base;

            if (xd > md) {
                allowed = false;
                break;
            }

            a /= base;
            b /= base;
        }

        if (allowed) {
            ++count;
        }
    }

    return count;
}

u64 local_delta(
    u64 old_m,
    u64 new_m,
    i128 y,
    u64 base
) {
    u64 a =
        miss_prefix(
            old_m,
            y,
            base
        );

    u64 b =
        miss_prefix(
            new_m,
            y,
            base
        );

    if (b < a) {
        std::cerr
            << "ERROR: negative local delta\n";

        std::exit(1);
    }

    return b - a;
}

u64 base_power(
    u64 base,
    std::size_t exponent
) {
    u64 result = 1;

    for (
        std::size_t i = 0;
        i < exponent;
        ++i
    ) {
        result *= base;
    }

    return result;
}

Polynomial build_polynomial(
    u64 target_m,
    i128 y,
    u64 base
) {
    std::vector<u64> digits =
        base_digits(
            target_m,
            base
        );

    Polynomial P;

    P.coeff.assign(
        digits.size(),
        0
    );

    u64 current = 0;

    for (
        std::size_t r = 0;
        r < digits.size();
        ++r
    ) {
        u64 step =
            base_power(
                base,
                r
            );

        for (
            u64 k = 0;
            k < digits[r];
            ++k
        ) {
            u64 next =
                current + step;

            u64 delta =
                local_delta(
                    current,
                    next,
                    y,
                    base
                );

            P.coeff[r] +=
                delta;

            current = next;
        }
    }

    if (current != target_m) {
        std::cerr
            << "ERROR: path did not reach target\n"
            << "target=" << target_m
            << " current=" << current
            << "\n";

        std::exit(1);
    }

    return P;
}

bool is_zero(
    const Polynomial &P
) {
    for (u64 c : P.coeff) {
        if (c != 0) {
            return false;
        }
    }

    return true;
}

u64 support_size(
    const Polynomial &P
) {
    u64 count = 0;

    for (u64 c : P.coeff) {
        if (c != 0) {
            ++count;
        }
    }

    return count;
}

std::size_t t_valuation(
    const Polynomial &P
) {
    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (P.coeff[i] != 0) {
            return i;
        }
    }

    return P.coeff.size();
}

void print_support(
    const Polynomial &P
) {
    std::cout
        << "support=[";

    bool first = true;

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (P.coeff[i] == 0) {
            continue;
        }

        if (!first) {
            std::cout << ",";
        }

        first = false;

        std::cout
            << "("
            << i
            << ","
            << P.coeff[i]
            << ")";
    }

    std::cout
        << "]\n";
}

void accumulate(
    Stats &dst,
    const Stats &src
) {
    dst.monomial +=
        src.monomial;

    dst.non_monomial +=
        src.non_monomial;

    dst.zero +=
        src.zero;

    dst.min_t_valuation =
        std::min(
            dst.min_t_valuation,
            src.min_t_valuation
        );

    dst.max_t_valuation =
        std::max(
            dst.max_t_valuation,
            src.max_t_valuation
        );

    dst.support_sum +=
        src.support_sum;

    dst.max_support =
        std::max(
            dst.max_support,
            src.max_support
        );
}

int main() {
    std::cout
        << "START EXPERIMENT 308\n"
        << "CORRECTED PREFIX-ORACLE REGRESSION\n"
        << "REBUILD PATH POLYNOMIAL STRUCTURE AFTER FIXING MSB DIGIT DP\n"
        << "\n";

    /*
        ============================================================
        PHASE 1
        Exhaustive small regression against brute force.
        ============================================================
    */

    u64 regression_cases = 0;
    u64 regression_fail = 0;

    const std::vector<u64> regression_bases = {
        2, 3, 5, 7
    };

    for (u64 base : regression_bases) {
        for (u64 m = 0; m <= 200; ++m) {
            for (u64 y = 0; y <= 250; ++y) {

                u64 fast =
                    miss_prefix(
                        m,
                        static_cast<i128>(y),
                        base
                    );

                u64 brute =
                    brute_miss_prefix(
                        m,
                        y,
                        base
                    );

                ++regression_cases;

                if (fast != brute) {
                    ++regression_fail;

                    if (regression_fail <= 20) {
                        std::cout
                            << "REGRESSION_FAIL\n"
                            << "base="
                            << base
                            << " m="
                            << m
                            << " y="
                            << y
                            << "\n"
                            << "fast="
                            << fast
                            << "\n"
                            << "brute="
                            << brute
                            << "\n";
                    }
                }
            }
        }
    }

    std::cout
        << "REGRESSION\n"
        << "cases="
        << regression_cases
        << "\n"
        << "fail="
        << regression_fail
        << "\n\n";

    if (regression_fail != 0) {
        std::cout
            << "FATAL: corrected prefix oracle failed regression\n"
            << "FINISHED EXPERIMENT 308\n";

        return 1;
    }

    /*
        ============================================================
        PHASE 2
        Genuine prime*prime semiprimes.
        ============================================================
    */

    const std::vector<CaseData> cases = {
        {81077, 162749},
        {125017, 174259},
        {57107, 88261},
        {52757, 148457},
        {19483, 34123},
        {123001, 188291},
        {97987, 112583},
        {76129, 192113},
        {41257, 73643},
        {137239, 140419},
        {64187, 115319},
        {162527, 184087},
        {66179, 124753},
        {87943, 158047},
        {42683, 52529},
        {107473, 149711},
        {75853, 82759},
        {150131, 175267},
        {99859, 124769},
        {75337, 162229},
        {124471, 168043},
        {51563, 87683},
        {52237, 142123},
        {19069, 28751},
        {117043, 187637},
        {92203, 112067},
        {75617, 185869},
        {40823, 67843},
        {134369, 136709},
        {63667, 109211}
    };

    const std::vector<u64> bases = {
        2, 3, 5, 7, 11, 13
    };

    const std::vector<int> offsets = {
        -4, -3, -2, -1,
         0,
         1, 2, 3, 4
    };

    Stats total;

    u64 total_polynomials = 0;
    u64 total_zero = 0;
    u64 total_nonzero = 0;

    for (
        std::size_t ci = 0;
        ci < cases.size();
        ++ci
    ) {
        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        u64 case_polynomials = 0;
        u64 case_zero = 0;
        u64 case_nonzero = 0;

        Stats case_stats;

        std::cout
            << "CASE "
            << ci
            << " p="
            << p
            << " q="
            << q
            << "\n";

        u64 root =
            integer_sqrt(N);

        for (u64 base : bases) {
            for (int offset : offsets) {

                i128 m_signed =
                    static_cast<i128>(root) +
                    static_cast<i128>(offset);

                if (m_signed <= 0) {
                    continue;
                }

                u64 m =
                    static_cast<u64>(
                        m_signed
                    );

                std::vector<u64> ns;

                ns.push_back(N);
                ns.push_back(m);
                ns.push_back(m + 1);

                if (
                    m != 0 &&
                    m <= N / m
                ) {
                    u64 D =
                        N - m * m;

                    if (D != 0) {
                        ns.push_back(D);
                    }
                }

                for (u64 n : ns) {
                    if (n == 0) {
                        continue;
                    }

                    i128 y =
                        static_cast<i128>(n) - 1;

                    /*
                        Independent local sanity check for
                        the first few path states.
                    */
                    u64 expected_final =
                        miss_prefix(
                            m,
                            y,
                            base
                        );

                    Polynomial P =
                        build_polynomial(
                            m,
                            y,
                            base
                        );

                    /*
                        Since the path starts at m=0 and
                        M_0(y)=1 for y>=0:

                            sum_r P_r = M_m(y)-1
                    */
                    u64 sum_coeff = 0;

                    for (u64 c : P.coeff) {
                        sum_coeff += c;
                    }

                    if (
                        y >= 0 &&
                        sum_coeff !=
                        expected_final - 1
                    ) {
                        std::cerr
                            << "ERROR: path telescoping failure\n"
                            << "case=" << ci
                            << " base=" << base
                            << " m=" << m
                            << " n=" << n
                            << "\n"
                            << "sum_coeff="
                            << sum_coeff
                            << "\n"
                            << "expected="
                            << expected_final - 1
                            << "\n";

                        return 1;
                    }

                    ++case_polynomials;
                    ++total_polynomials;

                    if (is_zero(P)) {
                        ++case_zero;
                        ++total_zero;

                        Stats current;
                        current.zero = 1;

                        accumulate(
                            case_stats,
                            current
                        );

                        accumulate(
                            total,
                            current
                        );

                        continue;
                    }

                    ++case_nonzero;
                    ++total_nonzero;

                    u64 s =
                        support_size(P);

                    std::size_t k =
                        t_valuation(P);

                    Stats current;

                    if (s == 1) {
                        current.monomial = 1;
                    } else {
                        current.non_monomial = 1;

                        std::cout
                            << "NON_MONOMIAL\n"
                            << "base="
                            << base
                            << " m="
                            << m
                            << " n="
                            << n
                            << "\n"
                            << "N="
                            << N
                            << "\n"
                            << "t_valuation="
                            << k
                            << "\n"
                            << "support_size="
                            << s
                            << "\n"
                            << "m_mod_base="
                            << (m % base)
                            << "\n";

                        print_support(P);

                        if (current.non_monomial) {
                            std::cout << "\n";
                        }
                    }

                    current.min_t_valuation =
                        static_cast<u64>(k);

                    current.max_t_valuation =
                        static_cast<u64>(k);

                    current.support_sum = s;
                    current.max_support = s;

                    accumulate(
                        case_stats,
                        current
                    );

                    accumulate(
                        total,
                        current
                    );
                }
            }
        }

        std::cout
            << "polynomials="
            << case_polynomials
            << "\n";

        std::cout
            << "zero_polynomials="
            << case_zero
            << "\n";

        std::cout
            << "nonzero_polynomials="
            << case_nonzero
            << "\n";

        std::cout
            << "monomial="
            << case_stats.monomial
            << "\n";

        std::cout
            << "non_monomial="
            << case_stats.non_monomial
            << "\n";

        std::cout
            << "min_t_valuation="
            << case_stats.min_t_valuation
            << "\n";

        std::cout
            << "max_t_valuation="
            << case_stats.max_t_valuation
            << "\n";

        std::cout
            << "max_support="
            << case_stats.max_support
            << "\n\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "regression_fail="
        << regression_fail
        << "\n";

    std::cout
        << "polynomials="
        << total_polynomials
        << "\n";

    std::cout
        << "zero_polynomials="
        << total_zero
        << "\n";

    std::cout
        << "nonzero_polynomials="
        << total_nonzero
        << "\n";

    std::cout
        << "monomial="
        << total.monomial
        << "\n";

    std::cout
        << "non_monomial="
        << total.non_monomial
        << "\n";

    std::cout
        << "min_t_valuation="
        << total.min_t_valuation
        << "\n";

    std::cout
        << "max_t_valuation="
        << total.max_t_valuation
        << "\n";

    std::cout
        << "support_size_sum="
        << total.support_sum
        << "\n";

    std::cout
        << "max_support="
        << total.max_support
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 308\n";

    return 0;
}
