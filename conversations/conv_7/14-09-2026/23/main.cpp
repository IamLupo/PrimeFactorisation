#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

struct CaseData {
    u64 p;
    u64 q;
};

struct Polynomial {
    std::vector<u64> coeff;
};

struct Stats {
    u64 expressions = 0;
    u64 invariant_tests = 0;

    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;

    u64 p_only = 0;
    u64 q_only = 0;

    u64 zero = 0;

    u64 control_base_equal = 0;
};

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }

    std::string s;

    while (x > 0) {
        int d =
            static_cast<int>(x % 10);

        s.push_back(
            static_cast<char>('0' + d)
        );

        x /= 10;
    }

    std::reverse(
        s.begin(),
        s.end()
    );

    return s;
}

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool negative =
        (x < 0);

    u128 v =
        negative
            ? static_cast<u128>(-x)
            : static_cast<u128>(x);

    std::string s =
        to_string_u128(v);

    if (negative) {
        s.insert(
            s.begin(),
            '-'
        );
    }

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

/*
    Correct prefix oracle.
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
    Closed coefficient formula.
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

    if (Q == 0) {
        return 0;
    }

    if (Q > mr) {
        return mr * W;
    }

    u64 R =
        y % pow;

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
        static_cast<u128>(
            W
        )
        +
        static_cast<u128>(
            A
        );

    return static_cast<u64>(
        result
    );
}

Polynomial build_polynomial(
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

    P.coeff.resize(
        digits.size()
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

/*
    Exact polynomial evaluation in signed 128-bit.

    For the tested sizes, the evaluated values fit comfortably
    inside 128 bits.
*/
i128 evaluate(
    const Polynomial &P,
    i128 t
) {
    i128 result = 0;

    for (
        std::size_t i = P.coeff.size();
        i-- > 0;
    ) {
        result =
            result * t +
            static_cast<i128>(
                P.coeff[i]
            );
    }

    return result;
}

void classify(
    i128 value,
    u64 N,
    u64 p,
    u64 q,
    Stats &stats
) {
    u64 a =
        abs_i128_to_u64(
            value
        );

    if (a == 0) {
        ++stats.zero;
        return;
    }

    u64 g =
        gcd_u64(
            a,
            N
        );

    if (g == 1) {
        ++stats.gcd1;
    } else if (g == N) {
        ++stats.gcdN;
    } else {
        ++stats.nontrivial;

        if (g == p) {
            ++stats.p_only;
        }

        if (g == q) {
            ++stats.q_only;
        }
    }
}

void accumulate(
    Stats &dst,
    const Stats &src
) {
    dst.expressions +=
        src.expressions;

    dst.invariant_tests +=
        src.invariant_tests;

    dst.gcd1 +=
        src.gcd1;

    dst.gcdN +=
        src.gcdN;

    dst.nontrivial +=
        src.nontrivial;

    dst.p_only +=
        src.p_only;

    dst.q_only +=
        src.q_only;

    dst.zero +=
        src.zero;

    dst.control_base_equal +=
        src.control_base_equal;
}

int main() {
    std::cout
        << "START EXPERIMENT 312\n"
        << "CROSS-BASE PATH POLYNOMIAL RELATIONS\n"
        << "CAN DIFFERENT RADIX REPRESENTATIONS EXPOSE A HIDDEN FACTOR?\n"
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

    /*
        For each base b we evaluate at:

            1
            b-1
            b
            b+1

        These remain small enough for exact i128 evaluation.
    */
    const int POINT_COUNT = 4;

    Stats total;

    u64 printed_hit = 0;

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
                We isolate the N-dependent boundary.
            */
            if (
                m == 0 ||
                m > N / m
            ) {
                continue;
            }

            u64 D =
                N - m * m;

            if (D == 0) {
                continue;
            }

            u64 y =
                D - 1;

            /*
                Build all base polynomials first.
            */
            std::vector<Polynomial> polys;
            polys.reserve(
                bases.size()
            );

            for (u64 base : bases) {
                polys.push_back(
                    build_polynomial(
                        m,
                        y,
                        base
                    )
                );
            }

            ++local.expressions;

            /*
                Base-value table:

                  values[b][point]
                where point index is

                    0 -> 1
                    1 -> b-1
                    2 -> b
                    3 -> b+1
            */
            std::vector<
                std::vector<i128>
            > values(
                bases.size(),
                std::vector<i128>(
                    POINT_COUNT,
                    0
                )
            );

            for (
                std::size_t bi = 0;
                bi < bases.size();
                ++bi
            ) {
                u64 b =
                    bases[bi];

                i128 points[
                    POINT_COUNT
                ] = {
                    1,
                    static_cast<i128>(b) - 1,
                    static_cast<i128>(b),
                    static_cast<i128>(b) + 1
                };

                for (
                    int pi = 0;
                    pi < POINT_COUNT;
                    ++pi
                ) {
                    values[bi][pi] =
                        evaluate(
                            polys[bi],
                            points[pi]
                        );
                }
            }

            /*
                CONTROL:
                P_b(1) should be the same for every base.
            */
            bool control_equal = true;

            for (
                std::size_t bi = 1;
                bi < bases.size();
                ++bi
            ) {
                if (
                    values[bi][0] !=
                    values[0][0]
                ) {
                    control_equal = false;
                    break;
                }
            }

            if (control_equal) {
                ++local.control_base_equal;
            }

            /*
                Compare every pair of bases.
            */
            for (
                std::size_t bi = 0;
                bi < bases.size();
                ++bi
            ) {
                for (
                    std::size_t bj = bi + 1;
                    bj < bases.size();
                    ++bj
                ) {
                    for (
                        int pi = 0;
                        pi < POINT_COUNT;
                        ++pi
                    ) {
                        i128 difference =
                            values[bi][pi] -
                            values[bj][pi];

                        ++local.invariant_tests;

                        if (
                            difference != 0
                        ) {
                            classify(
                                difference,
                                N,
                                p,
                                q,
                                local
                            );

                            /*
                                classify() increments
                                gcd counters but does not
                                tell us explicitly whether
                                the result was a factor.
                                Recalculate for witness output.
                            */
                            u64 g =
                                gcd_u64(
                                    abs_i128_to_u64(
                                        difference
                                    ),
                                    N
                                );

                            if (
                                (
                                    g == p ||
                                    g == q ||
                                    (
                                        g != 1 &&
                                        g != N
                                    )
                                ) &&
                                printed_hit < 20
                            ) {
                                ++printed_hit;

                                std::cout
                                    << "CROSS_BASE_HIT\n"
                                    << "m="
                                    << m
                                    << " D="
                                    << D
                                    << "\n"
                                    << "base_a="
                                    << bases[bi]
                                    << "\n"
                                    << "base_b="
                                    << bases[bj]
                                    << "\n"
                                    << "point_index="
                                    << pi
                                    << "\n"
                                    << "value_a="
                                    << to_string_i128(
                                        values[bi][pi]
                                    )
                                    << "\n"
                                    << "value_b="
                                    << to_string_i128(
                                        values[bj][pi]
                                    )
                                    << "\n"
                                    << "difference="
                                    << to_string_i128(
                                        difference
                                    )
                                    << "\n"
                                    << "gcd="
                                    << g
                                    << "\n";
                            }
                        }
                    }
                }
            }

            /*
                Additional normalization:
                compare ratios at the same t by
                cross multiplication:

                    A_i(1)*A_j(t)
                    -
                    A_j(1)*A_i(t)

                This asks whether the different base
                polynomials have the same projective
                scaling modulo a factor.
            */
            for (
                std::size_t bi = 0;
                bi < bases.size();
                ++bi
            ) {
                for (
                    std::size_t bj = bi + 1;
                    bj < bases.size();
                    ++bj
                ) {
                    for (
                        int pi = 1;
                        pi < POINT_COUNT;
                        ++pi
                    ) {
                        i128 relation =
                            values[bi][0] *
                            values[bj][pi]
                            -
                            values[bj][0] *
                            values[bi][pi];

                        ++local.invariant_tests;

                        if (
                            relation != 0
                        ) {
                            classify(
                                relation,
                                N,
                                p,
                                q,
                                local
                            );

                            u64 g =
                                gcd_u64(
                                    abs_i128_to_u64(
                                        relation
                                    ),
                                    N
                                );

                            if (
                                (
                                    g == p ||
                                    g == q ||
                                    (
                                        g != 1 &&
                                        g != N
                                    )
                                ) &&
                                printed_hit < 20
                            ) {
                                ++printed_hit;

                                std::cout
                                    << "CROSS_BASE_PROJECTIVE_HIT\n"
                                    << "m="
                                    << m
                                    << " D="
                                    << D
                                    << "\n"
                                    << "base_a="
                                    << bases[bi]
                                    << "\n"
                                    << "base_b="
                                    << bases[bj]
                                    << "\n"
                                    << "point_index="
                                    << pi
                                    << "\n"
                                    << "relation="
                                    << to_string_i128(
                                        relation
                                    )
                                    << "\n"
                                    << "gcd="
                                    << g
                                    << "\n";
                            }
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
            << "invariant_tests="
            << local.invariant_tests
            << "\n";

        std::cout
            << "gcd1="
            << local.gcd1
            << "\n";

        std::cout
            << "gcdN="
            << local.gcdN
            << "\n";

        std::cout
            << "nontrivial="
            << local.nontrivial
            << "\n";

        std::cout
            << "p_only="
            << local.p_only
            << "\n";

        std::cout
            << "q_only="
            << local.q_only
            << "\n";

        std::cout
            << "zero="
            << local.zero
            << "\n";

        std::cout
            << "control_base_equal="
            << local.control_base_equal
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
        << "invariant_tests="
        << total.invariant_tests
        << "\n";

    std::cout
        << "gcd1="
        << total.gcd1
        << "\n";

    std::cout
        << "gcdN="
        << total.gcdN
        << "\n";

    std::cout
        << "nontrivial="
        << total.nontrivial
        << "\n";

    std::cout
        << "p_only="
        << total.p_only
        << "\n";

    std::cout
        << "q_only="
        << total.q_only
        << "\n";

    std::cout
        << "zero="
        << total.zero
        << "\n";

    std::cout
        << "control_base_equal="
        << total.control_base_equal
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 312\n";

    return 0;
}
