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

    u64 coeff_gcd1 = 0;
    u64 coeff_gcdN = 0;
    u64 coeff_nontrivial = 0;
    u64 coeff_p_only = 0;
    u64 coeff_q_only = 0;

    u64 difference_gcd1 = 0;
    u64 difference_gcdN = 0;
    u64 difference_nontrivial = 0;
    u64 difference_p_only = 0;
    u64 difference_q_only = 0;

    u64 pair_gcd1 = 0;
    u64 pair_gcdN = 0;
    u64 pair_nontrivial = 0;
    u64 pair_p_only = 0;
    u64 pair_q_only = 0;
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

u64 abs_i128_to_u64(
    i128 x
) {
    u128 v =
        x < 0
            ? static_cast<u128>(-x)
            : static_cast<u128>(x);

    return static_cast<u64>(v);
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

u64 base_power(
    u64 base,
    std::size_t r
) {
    u64 result = 1;

    for (
        std::size_t i = 0;
        i < r;
        ++i
    ) {
        result *= base;
    }

    return result;
}

/*
    Direct closed coefficient formula.

    c_r =
        0                                  if Q=0
        m_r W                              if Q>m_r
        (Q-1)W + A                         otherwise

    with

        Q = floor(y / b^r)
        R = y mod b^r
        W = prod_{j<r}(m_j+1)
        A = M_{m mod b^r}(R)
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

    u64 W =
        miss_prefix(
            low,
            static_cast<i128>(low),
            base
        );

    u64 Q =
        y / pow;

    u64 R =
        y % pow;

    if (Q == 0) {
        return 0;
    }

    if (Q > mr) {
        return mr * W;
    }

    u64 A =
        miss_prefix(
            low,
            static_cast<i128>(R),
            base
        );

    u128 result =
        static_cast<u128>(
            Q - 1
        ) *
        static_cast<u128>(W)
        +
        static_cast<u128>(A);

    return static_cast<u64>(
        result
    );
}

Polynomial build_closed_polynomial(
    u64 m,
    u64 y,
    u64 base
) {
    std::vector<u64> digits =
        base_digits(
            m,
            base
        );

    Polynomial P;

    P.coeff.assign(
        digits.size(),
        0
    );

    for (
        std::size_t r = 0;
        r < digits.size();
        ++r
    ) {
        P.coeff[r] =
            closed_coefficient(
                m,
                y,
                base,
                r
            );
    }

    return P;
}

u64 polynomial_gcd(
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

void classify_gcd(
    u64 g,
    u64 N,
    u64 p,
    u64 q,
    u64 &gcd1,
    u64 &gcdN,
    u64 &nontrivial,
    u64 &p_only,
    u64 &q_only
) {
    if (g == 1) {
        ++gcd1;
    } else if (g == N) {
        ++gcdN;
    } else {
        ++nontrivial;

        if (g == p) {
            ++p_only;
        }

        if (g == q) {
            ++q_only;
        }
    }
}

void accumulate(
    Stats &dst,
    const Stats &src
) {
    dst.expressions +=
        src.expressions;

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

    dst.difference_gcd1 +=
        src.difference_gcd1;

    dst.difference_gcdN +=
        src.difference_gcdN;

    dst.difference_nontrivial +=
        src.difference_nontrivial;

    dst.difference_p_only +=
        src.difference_p_only;

    dst.difference_q_only +=
        src.difference_q_only;

    dst.pair_gcd1 +=
        src.pair_gcd1;

    dst.pair_gcdN +=
        src.pair_gcdN;

    dst.pair_nontrivial +=
        src.pair_nontrivial;

    dst.pair_p_only +=
        src.pair_p_only;

    dst.pair_q_only +=
        src.pair_q_only;
}

void print_polynomial(
    const std::string &name,
    const Polynomial &P
) {
    std::cout
        << name
        << "=[";

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (i != 0) {
            std::cout << ",";
        }

        std::cout
            << P.coeff[i];
    }

    std::cout
        << "]\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 310\n"
        << "FACTOR-SENSITIVE N-m^2 COEFFICIENT STRUCTURE\n"
        << "DOES THE D=N-m^2 BOUNDARY CONTAIN CRT FACTOR INFORMATION?\n"
        << "\n";

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

                /*
                    We only care about the N-m^2
                    boundary here.

                    D is guaranteed positive for
                    m <= floor(sqrt(N)).
                */
                if (m == 0 || m > N / m) {
                    continue;
                }

                u64 D =
                    N - m * m;

                if (D == 0) {
                    continue;
                }

                ++local.expressions;

                /*
                    D-boundary polynomial:
                        y = D-1
                */
                u64 yD =
                    D - 1;

                Polynomial PD =
                    build_closed_polynomial(
                        m,
                        yD,
                        base
                    );

                u64 gD =
                    polynomial_gcd(PD);

                classify_gcd(
                    gD,
                    N,
                    p,
                    q,
                    local.coeff_gcd1,
                    local.coeff_gcdN,
                    local.coeff_nontrivial,
                    local.coeff_p_only,
                    local.coeff_q_only
                );

                /*
                    Comparison polynomials.
                */
                u64 yM =
                    m - 1;

                Polynomial PM =
                    build_closed_polynomial(
                        m,
                        yM,
                        base
                    );

                Polynomial PN =
                    build_closed_polynomial(
                        m,
                        N - 1,
                        base
                    );

                Polynomial PM1 =
                    build_closed_polynomial(
                        m,
                        m,
                        base
                    );

                /*
                    Difference:
                        D-polynomial - m-polynomial
                */
                const std::size_t L =
                    std::max(
                        PD.coeff.size(),
                        PM.coeff.size()
                    );

                std::vector<i128> diff(
                    L,
                    0
                );

                for (
                    std::size_t r = 0;
                    r < L;
                    ++r
                ) {
                    u64 a =
                        r < PD.coeff.size()
                            ? PD.coeff[r]
                            : 0;

                    u64 b =
                        r < PM.coeff.size()
                            ? PM.coeff[r]
                            : 0;

                    diff[r] =
                        static_cast<i128>(a) -
                        static_cast<i128>(b);
                }

                u64 diff_gcd = 0;

                for (i128 v : diff) {
                    diff_gcd =
                        gcd_u64(
                            diff_gcd,
                            abs_i128_to_u64(v)
                        );
                }

                classify_gcd(
                    diff_gcd,
                    N,
                    p,
                    q,
                    local.difference_gcd1,
                    local.difference_gcdN,
                    local.difference_nontrivial,
                    local.difference_p_only,
                    local.difference_q_only
                );

                /*
                    Pairwise coefficient gcds.

                    For every pair r<s with at least one
                    nonzero coefficient, inspect

                        gcd(c_r - c_s, N).

                    This searches for a hidden CRT separation
                    that would not appear in the gcd of all
                    coefficients.
                */
                for (
                    std::size_t r = 0;
                    r < PD.coeff.size();
                    ++r
                ) {
                    for (
                        std::size_t s = r + 1;
                        s < PD.coeff.size();
                        ++s
                    ) {
                        i128 v =
                            static_cast<i128>(
                                PD.coeff[r]
                            ) -
                            static_cast<i128>(
                                PD.coeff[s]
                            );

                        u64 g =
                            gcd_u64(
                                abs_i128_to_u64(v),
                                N
                            );

                        if (g == 1) {
                            ++local.pair_gcd1;
                        } else if (g == N) {
                            ++local.pair_gcdN;
                        } else {
                            ++local.pair_nontrivial;

                            if (g == p) {
                                ++local.pair_p_only;
                            }

                            if (g == q) {
                                ++local.pair_q_only;
                            }

                            /*
                                First diagnostic witness.
                            */
                            static u64 printed = 0;

                            if (printed < 20) {
                                ++printed;

                                std::cout
                                    << "PAIR_HIT\n"
                                    << "base="
                                    << base
                                    << " m="
                                    << m
                                    << " D="
                                    << D
                                    << "\n"
                                    << "r="
                                    << r
                                    << " s="
                                    << s
                                    << "\n"
                                    << "c_r="
                                    << PD.coeff[r]
                                    << "\n"
                                    << "c_s="
                                    << PD.coeff[s]
                                    << "\n"
                                    << "difference="
                                    << to_string_i128(v)
                                    << "\n"
                                    << "gcd="
                                    << g
                                    << "\n";
                            }
                        }
                    }
                }

                /*
                    Print the D polynomial only for small
                    potentially interesting gcd cases.
                */
                if (
                    gD != 1 &&
                    gD != N
                ) {
                    static u64 printed_coeff = 0;

                    if (printed_coeff < 10) {
                        ++printed_coeff;

                        std::cout
                            << "D_COEFFICIENT_GCD\n"
                            << "base="
                            << base
                            << " m="
                            << m
                            << " D="
                            << D
                            << "\n"
                            << "gcd="
                            << gD
                            << "\n";

                        print_polynomial(
                            "PD",
                            PD
                        );
                    }
                }

                /*
                    Avoid warnings about the comparison
                    polynomials being intentionally materialized.
                */
                (void)PN;
                (void)PM1;
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
            << "coefficient_gcd1="
            << local.coeff_gcd1
            << "\n";

        std::cout
            << "coefficient_gcdN="
            << local.coeff_gcdN
            << "\n";

        std::cout
            << "coefficient_nontrivial="
            << local.coeff_nontrivial
            << "\n";

        std::cout
            << "coefficient_p_only="
            << local.coeff_p_only
            << "\n";

        std::cout
            << "coefficient_q_only="
            << local.coeff_q_only
            << "\n";

        std::cout
            << "difference_gcd1="
            << local.difference_gcd1
            << "\n";

        std::cout
            << "difference_gcdN="
            << local.difference_gcdN
            << "\n";

        std::cout
            << "difference_nontrivial="
            << local.difference_nontrivial
            << "\n";

        std::cout
            << "difference_p_only="
            << local.difference_p_only
            << "\n";

        std::cout
            << "difference_q_only="
            << local.difference_q_only
            << "\n";

        std::cout
            << "\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "expressions="
        << total.expressions
        << "\n";

    std::cout
        << "coefficient_gcd1="
        << total.coeff_gcd1
        << "\n";

    std::cout
        << "coefficient_gcdN="
        << total.coeff_gcdN
        << "\n";

    std::cout
        << "coefficient_nontrivial="
        << total.coeff_nontrivial
        << "\n";

    std::cout
        << "coefficient_p_only="
        << total.coeff_p_only
        << "\n";

    std::cout
        << "coefficient_q_only="
        << total.coeff_q_only
        << "\n";

    std::cout
        << "difference_gcd1="
        << total.difference_gcd1
        << "\n";

    std::cout
        << "difference_gcdN="
        << total.difference_gcdN
        << "\n";

    std::cout
        << "difference_nontrivial="
        << total.difference_nontrivial
        << "\n";

    std::cout
        << "difference_p_only="
        << total.difference_p_only
        << "\n";

    std::cout
        << "difference_q_only="
        << total.difference_q_only
        << "\n";

    std::cout
        << "pair_gcd1="
        << total.pair_gcd1
        << "\n";

    std::cout
        << "pair_gcdN="
        << total.pair_gcdN
        << "\n";

    std::cout
        << "pair_nontrivial="
        << total.pair_nontrivial
        << "\n";

    std::cout
        << "pair_p_only="
        << total.pair_p_only
        << "\n";

    std::cout
        << "pair_q_only="
        << total.pair_q_only
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 310\n";

    return 0;
}
