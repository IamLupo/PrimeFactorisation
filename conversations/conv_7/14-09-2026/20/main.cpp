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
    u64 expressions = 0;
    u64 formula_fail = 0;

    u64 coeff_gcd1 = 0;
    u64 coeff_gcdN = 0;
    u64 coeff_nontrivial = 0;

    u64 coeff_p_only = 0;
    u64 coeff_q_only = 0;
};

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool neg = (x < 0);

    u128 v =
        neg
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

    if (neg) {
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
    Correct prefix count:

        M_m(y)
        =
        #{x : 0 <= x <= y and x <=_b m}
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
        lower_product[i]
        =
        product_{j<i}(m_j+1)
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

        if (yi <= mi) {
            tight = true;
        } else {
            tight = false;
            break;
        }
    }

    if (tight) {
        ++answer;
    }

    return static_cast<u64>(
        answer
    );
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
    std::size_t r
) {
    u64 v = 1;

    for (
        std::size_t i = 0;
        i < r;
        ++i
    ) {
        v *= base;
    }

    return v;
}

/*
    Build the coefficient vector from the original
    canonical path definition.

    This is the independent reference for the closed formula.
*/
Polynomial build_path_polynomial(
    u64 m,
    i128 y,
    u64 base
) {
    std::vector<u64> md =
        base_digits(
            m,
            base
        );

    Polynomial P;

    P.coeff.assign(
        md.size(),
        0
    );

    u64 current = 0;

    for (
        std::size_t r = 0;
        r < md.size();
        ++r
    ) {
        u64 step =
            base_power(
                base,
                r
            );

        for (
            u64 k = 0;
            k < md[r];
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

    if (current != m) {
        std::cerr
            << "ERROR: path did not reach m\n";

        std::exit(1);
    }

    return P;
}

/*
    Closed formula for coefficient c_r.

    Let

        m_r       = digit r of m
        low       = m mod b^r
        W_r       = product_{j<r}(m_j+1)

    and

        y = Q*b^r + R.

    A newly-created point at digit r has

        x = a*b^r + L

    with

        1 <= a <= m_r
        L <=_b low.

    For L <= R, a can reach Q.
    For L > R, a can reach Q-1.
*/
u64 closed_coefficient(
    u64 m,
    u64 y,
    u64 base,
    std::size_t r
) {
    u64 pow =
        base_power(
            base,
            r
        );

    std::vector<u64> md =
        base_digits(
            m,
            base
        );

    u64 mr =
        r < md.size()
            ? md[r]
            : 0;

    if (mr == 0) {
        return 0;
    }

    u64 low =
        m % pow;

    /*
        W_r = number of admissible lower-digit
        combinations.
    */
    u64 W =
        miss_prefix(
            low,
            static_cast<i128>(low),
            base
        );

    /*
        y = Q*b^r + R
    */
    u64 Q = y / pow;
    u64 R = y % pow;

    /*
        Number of admissible lower patterns L <= R.
    */
    u64 A =
        miss_prefix(
            low,
            static_cast<i128>(R),
            base
        );

    u64 a1 =
        std::min(
            mr,
            Q
        );

    u64 q_minus_one =
        Q == 0
            ? 0
            : Q - 1;

    u64 a2 =
        std::min(
            mr,
            q_minus_one
        );

    u128 result =
        static_cast<u128>(a1) * A +
        static_cast<u128>(a2) *
        static_cast<u128>(W - A);

    return static_cast<u64>(
        result
    );
}

u64 coefficient_gcd(
    const Polynomial &P
) {
    u64 g = 0;

    for (u64 c : P.coeff) {
        g =
            gcd_u64(
                g,
                c
            );
    }

    return g;
}

bool has_factor(
    const Polynomial &P,
    u64 factor
) {
    if (factor == 0) {
        return false;
    }

    for (u64 c : P.coeff) {
        if (c % factor != 0) {
            return false;
        }
    }

    return true;
}

void print_coefficients(
    const Polynomial &P
) {
    std::cout
        << "[";

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (i != 0) {
            std::cout
                << ",";
        }

        std::cout
            << P.coeff[i];
    }

    std::cout
        << "]\n";
}

void accumulate(
    Stats &dst,
    const Stats &src
) {
    dst.expressions +=
        src.expressions;

    dst.formula_fail +=
        src.formula_fail;

    dst.coeff_gcd1 +=
        src.coeff_gcd1;

    dst.coeff_gcdN +=
        src.coeff_gcdN;

    dst.coeff_nontrivial +=
        src.coeff_nontrivial;

    dst.coeff_p_only +=
        src.coeff_p_only;

    dst.coeff_q_only +=
        src.coeff_q_only;
}

int main() {
    std::cout
        << "START EXPERIMENT 309\n"
        << "CLOSED FORM FOR CANONICAL PATH COEFFICIENTS\n"
        << "CAN THE ENTIRE PATH POLYNOMIAL BE WRITTEN DIRECTLY FROM DIGITS?\n"
        << "\n";

    /*
        ------------------------------------------------------------
        Phase I: independent formula regression.
        ------------------------------------------------------------
    */

    u64 regression_cases = 0;
    u64 regression_fail = 0;

    const std::vector<u64> regression_bases = {
        2, 3, 5, 7, 11
    };

    for (u64 base : regression_bases) {
        for (u64 m = 0; m <= 300; ++m) {
            for (u64 y = 0; y <= 500; ++y) {

                Polynomial P =
                    build_path_polynomial(
                        m,
                        static_cast<i128>(y),
                        base
                    );

                for (
                    std::size_t r = 0;
                    r < P.coeff.size();
                    ++r
                ) {
                    u64 expected =
                        closed_coefficient(
                            m,
                            y,
                            base,
                            r
                        );

                    ++regression_cases;

                    if (
                        P.coeff[r] !=
                        expected
                    ) {
                        ++regression_fail;

                        if (regression_fail <= 20) {
                            std::cout
                                << "FORMULA_FAIL\n"
                                << "base="
                                << base
                                << " m="
                                << m
                                << " y="
                                << y
                                << " r="
                                << r
                                << "\n"
                                << "path="
                                << P.coeff[r]
                                << "\n"
                                << "closed="
                                << expected
                                << "\n";
                        }
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
            << "FATAL: closed coefficient formula failed\n"
            << "FINISHED EXPERIMENT 309\n";

        return 1;
    }

    /*
        ------------------------------------------------------------
        Phase II: semiprime scan.
        ------------------------------------------------------------
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

    for (
        std::size_t ci = 0;
        ci < cases.size();
        ++ci
    ) {
        const u64 p =
            cases[ci].p;

        const u64 q =
            cases[ci].q;

        const u64 N =
            p * q;

        Stats local;

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

                    ++local.expressions;

                    i128 y_signed =
                        static_cast<i128>(n) - 1;

                    Polynomial P =
                        build_path_polynomial(
                            m,
                            y_signed,
                            base
                        );

                    /*
                        Validate every coefficient directly
                        against the closed formula.
                    */
                    for (
                        std::size_t r = 0;
                        r < P.coeff.size();
                        ++r
                    ) {
                        u64 expected =
                            closed_coefficient(
                                m,
                                static_cast<u64>(y_signed),
                                base,
                                r
                            );

                        if (
                            P.coeff[r] !=
                            expected
                        ) {
                            ++local.formula_fail;

                            if (
                                local.formula_fail <=
                                5
                            ) {
                                std::cout
                                    << "SEMIPRIME_FORMULA_FAIL\n"
                                    << "base="
                                    << base
                                    << " m="
                                    << m
                                    << " n="
                                    << n
                                    << " r="
                                    << r
                                    << "\n"
                                    << "path="
                                    << P.coeff[r]
                                    << "\n"
                                    << "closed="
                                    << expected
                                    << "\n";
                            }
                        }
                    }

                    u64 g =
                        coefficient_gcd(P);

                    if (g == 1) {
                        ++local.coeff_gcd1;
                    } else if (g == N) {
                        ++local.coeff_gcdN;
                    } else {
                        ++local.coeff_nontrivial;

                        if (g == p) {
                            ++local.coeff_p_only;
                        }

                        if (g == q) {
                            ++local.coeff_q_only;
                        }
                    }
                }
            }
        }

        accumulate(
            total,
            local
        );

        std::cout
            << "expressions="
            << local.expressions
            << "\n";

        std::cout
            << "formula_fail="
            << local.formula_fail
            << "\n";

        std::cout
            << "coeff_gcd1="
            << local.coeff_gcd1
            << "\n";

        std::cout
            << "coeff_gcdN="
            << local.coeff_gcdN
            << "\n";

        std::cout
            << "coeff_nontrivial="
            << local.coeff_nontrivial
            << "\n";

        std::cout
            << "coeff_p_only="
            << local.coeff_p_only
            << "\n";

        std::cout
            << "coeff_q_only="
            << local.coeff_q_only
            << "\n\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "expressions="
        << total.expressions
        << "\n";

    std::cout
        << "formula_fail="
        << total.formula_fail
        << "\n";

    std::cout
        << "coeff_gcd1="
        << total.coeff_gcd1
        << "\n";

    std::cout
        << "coeff_gcdN="
        << total.coeff_gcdN
        << "\n";

    std::cout
        << "coeff_nontrivial="
        << total.coeff_nontrivial
        << "\n";

    std::cout
        << "coeff_p_only="
        << total.coeff_p_only
        << "\n";

    std::cout
        << "coeff_q_only="
        << total.coeff_q_only
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 309\n";

    return 0;
}
