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

    u64 determinant_tests = 0;

    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;

    u64 p_only = 0;
    u64 q_only = 0;

    u64 zero_det = 0;

    u64 adjacent_gcd1 = 0;
    u64 adjacent_gcdN = 0;
    u64 adjacent_nontrivial = 0;
    u64 adjacent_p_only = 0;
    u64 adjacent_q_only = 0;
};

std::string to_string_u128(
    u128 x
) {
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

std::string to_string_i128(
    i128 x
) {
    if (x == 0) {
        return "0";
    }

    bool negative = x < 0;

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
    c_r*c_{s+1} - c_{r+1}*c_s

    Computed using signed 128-bit arithmetic.
*/
i128 determinant(
    const Polynomial &P,
    std::size_t r,
    std::size_t s
) {
    i128 a =
        static_cast<i128>(
            P.coeff[r]
        );

    i128 b =
        static_cast<i128>(
            P.coeff[r + 1]
        );

    i128 c =
        static_cast<i128>(
            P.coeff[s]
        );

    i128 d =
        static_cast<i128>(
            P.coeff[s + 1]
        );

    return
        a * d -
        b * c;
}

void classify(
    i128 value,
    u64 N,
    u64 p,
    u64 q,
    u64 &gcd1,
    u64 &gcdN,
    u64 &nontrivial,
    u64 &p_only,
    u64 &q_only,
    u64 &zero
) {
    u64 a =
        abs_i128_to_u64(
            value
        );

    if (a == 0) {
        ++zero;
        return;
    }

    u64 g =
        gcd_u64(
            a,
            N
        );

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

    dst.determinant_tests +=
        src.determinant_tests;

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

    dst.zero_det +=
        src.zero_det;

    dst.adjacent_gcd1 +=
        src.adjacent_gcd1;

    dst.adjacent_gcdN +=
        src.adjacent_gcdN;

    dst.adjacent_nontrivial +=
        src.adjacent_nontrivial;

    dst.adjacent_p_only +=
        src.adjacent_p_only;

    dst.adjacent_q_only +=
        src.adjacent_q_only;
}

void print_polynomial(
    const Polynomial &P
) {
    std::cout
        << "coeff=[";
    
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

int main() {
    std::cout
        << "START EXPERIMENT 311\n"
        << "MULTIPLICATIVE RELATIONS BETWEEN PATH COEFFICIENTS\n"
        << "DO 2x2 COEFFICIENT DETERMINANTS EXPOSE A HIDDEN PRIME?\n"
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

    u64 printed_hits = 0;

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

                /*
                    We test only the genuinely
                    N-dependent boundary.
                */
                u64 y =
                    D - 1;

                Polynomial P =
                    build_polynomial(
                        m,
                        y,
                        base
                    );

                ++local.expressions;

                std::size_t L =
                    P.coeff.size();

                /*
                    Adjacent determinant:

                        c_r*c_{r+2}
                        -
                        c_{r+1}^2
                */
                for (
                    std::size_t r = 0;
                    r + 2 < L;
                    ++r
                ) {
                    i128 det =
                        determinant(
                            P,
                            r,
                            r + 1
                        );

                    u64 dummy1 = 0;
                    u64 dummyN = 0;
                    u64 dummyNontrivial = 0;
                    u64 dummyP = 0;
                    u64 dummyQ = 0;
                    u64 dummyZero = 0;

                    classify(
                        det,
                        N,
                        p,
                        q,
                        dummy1,
                        dummyN,
                        dummyNontrivial,
                        dummyP,
                        dummyQ,
                        dummyZero
                    );

                    if (dummy1) {
                        ++local.adjacent_gcd1;
                    }

                    if (dummyN) {
                        ++local.adjacent_gcdN;
                    }

                    if (dummyNontrivial) {
                        ++local.adjacent_nontrivial;
                    }

                    if (dummyP) {
                        ++local.adjacent_p_only;
                    }

                    if (dummyQ) {
                        ++local.adjacent_q_only;
                    }
                }

                /*
                    All pairwise 2x2 determinants.

                    For r < s:

                        c_r c_{s+1}
                        -
                        c_{r+1} c_s
                */
                for (
                    std::size_t r = 0;
                    r + 1 < L;
                    ++r
                ) {
                    for (
                        std::size_t s = r + 1;
                        s + 1 < L;
                        ++s
                    ) {
                        i128 det =
                            determinant(
                                P,
                                r,
                                s
                            );

                        ++local.determinant_tests;

                        u64 before1 =
                            local.gcd1;

                        u64 beforeN =
                            local.gcdN;

                        u64 beforeNontrivial =
                            local.nontrivial;

                        u64 beforeP =
                            local.p_only;

                        u64 beforeQ =
                            local.q_only;

                        u64 beforeZero =
                            local.zero_det;

                        classify(
                            det,
                            N,
                            p,
                            q,
                            local.gcd1,
                            local.gcdN,
                            local.nontrivial,
                            local.p_only,
                            local.q_only,
                            local.zero_det
                        );

                        bool isP =
                            local.p_only >
                            beforeP;

                        bool isQ =
                            local.q_only >
                            beforeQ;

                        bool isNontrivial =
                            local.nontrivial >
                            beforeNontrivial;

                        if (
                            (
                                isP ||
                                isQ ||
                                isNontrivial
                            ) &&
                            printed_hits < 20
                        ) {
                            ++printed_hits;

                            u64 g =
                                gcd_u64(
                                    abs_i128_to_u64(
                                        det
                                    ),
                                    N
                                );

                            std::cout
                                << "DETERMINANT_HIT\n"
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
                                << P.coeff[r]
                                << "\n"
                                << "c_r1="
                                << P.coeff[r + 1]
                                << "\n"
                                << "c_s="
                                << P.coeff[s]
                                << "\n"
                                << "c_s1="
                                << P.coeff[s + 1]
                                << "\n"
                                << "det="
                                << to_string_i128(
                                    det
                                )
                                << "\n"
                                << "gcd="
                                << g
                                << "\n";

                            print_polynomial(
                                P
                            );
                        }

                        (void)before1;
                        (void)beforeN;
                        (void)beforeZero;
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
            << "determinant_tests="
            << local.determinant_tests
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
            << "zero_det="
            << local.zero_det
            << "\n";

        std::cout
            << "adjacent_gcd1="
            << local.adjacent_gcd1
            << "\n";

        std::cout
            << "adjacent_gcdN="
            << local.adjacent_gcdN
            << "\n";

        std::cout
            << "adjacent_nontrivial="
            << local.adjacent_nontrivial
            << "\n";

        std::cout
            << "adjacent_p_only="
            << local.adjacent_p_only
            << "\n";

        std::cout
            << "adjacent_q_only="
            << local.adjacent_q_only
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
        << "determinant_tests="
        << total.determinant_tests
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
        << "zero_det="
        << total.zero_det
        << "\n";

    std::cout
        << "adjacent_gcd1="
        << total.adjacent_gcd1
        << "\n";

    std::cout
        << "adjacent_gcdN="
        << total.adjacent_gcdN
        << "\n";

    std::cout
        << "adjacent_nontrivial="
        << total.adjacent_nontrivial
        << "\n";

    std::cout
        << "adjacent_p_only="
        << total.adjacent_p_only
        << "\n";

    std::cout
        << "adjacent_q_only="
        << total.adjacent_q_only
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 311\n";

    return 0;
}
