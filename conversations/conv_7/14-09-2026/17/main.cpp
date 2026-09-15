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

struct PrimitivePolynomial {
    std::vector<i128> coeff;
};

struct ReducedPolynomial {
    std::vector<i128> coeff;
};

struct Stats {
    u64 q_zero = 0;
    u64 q_one_zero = 0;
    u64 q_minus_one_zero = 0;
    u64 q_two_zero = 0;
    u64 q_minus_two_zero = 0;

    u64 q_mod_pattern_equal = 0;

    u64 monomial = 0;
    u64 binomial = 0;

    u64 degree1 = 0;
    u64 degree2 = 0;
    u64 degree3plus = 0;
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
        int digit =
            static_cast<int>(v % 10);

        s.push_back(
            static_cast<char>('0' + digit)
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
    std::vector<u64> digits;

    if (x == 0) {
        digits.push_back(0);
        return digits;
    }

    while (x > 0) {
        digits.push_back(
            x % base
        );

        x /= base;
    }

    return digits;
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

    std::vector<u128> suffix(
        L + 1,
        1
    );

    for (
        std::size_t i = L;
        i-- > 0;
    ) {
        suffix[i] =
            suffix[i + 1] *
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

        u64 upper =
            std::min(
                yi,
                mi + 1
            );

        for (
            u64 xi = 0;
            xi < upper;
            ++xi
        ) {
            answer +=
                suffix[pos + 1];
        }

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
            << "ERROR: negative local derivative\n";

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
            << "ERROR: path did not reach target\n";

        std::exit(1);
    }

    return P;
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

PrimitivePolynomial primitive_part(
    const Polynomial &P,
    u64 g
) {
    PrimitivePolynomial R;

    if (g == 0) {
        return R;
    }

    R.coeff.resize(
        P.coeff.size()
    );

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        R.coeff[i] =
            static_cast<i128>(
                P.coeff[i] / g
            );
    }

    while (
        R.coeff.size() > 1 &&
        R.coeff.back() == 0
    ) {
        R.coeff.pop_back();
    }

    return R;
}

/*
    Experiment 305 established that every nonzero
    primitive polynomial has coefficient 0 equal to 0.

        P(t) = t Q(t)

    Remove the forced t factor.
*/
ReducedPolynomial divide_forced_t(
    const PrimitivePolynomial &P
) {
    ReducedPolynomial Q;

    if (P.coeff.empty()) {
        return Q;
    }

    if (P.coeff[0] != 0) {
        std::cerr
            << "ERROR: primitive polynomial does not "
            << "have forced t factor\n";

        std::exit(1);
    }

    Q.coeff.resize(
        P.coeff.size() - 1
    );

    for (
        std::size_t i = 1;
        i < P.coeff.size();
        ++i
    ) {
        Q.coeff[i - 1] =
            P.coeff[i];
    }

    return Q;
}

i128 evaluate(
    const ReducedPolynomial &Q,
    i128 x
) {
    i128 result = 0;

    for (
        std::size_t i = Q.coeff.size();
        i-- > 0;
    ) {
        result =
            result * x +
            Q.coeff[i];
    }

    return result;
}

std::size_t degree(
    const ReducedPolynomial &Q
) {
    if (Q.coeff.empty()) {
        return 0;
    }

    return Q.coeff.size() - 1;
}

bool is_zero(
    const ReducedPolynomial &Q
) {
    for (i128 c : Q.coeff) {
        if (c != 0) {
            return false;
        }
    }

    return true;
}

bool is_monomial(
    const ReducedPolynomial &Q
) {
    bool seen = false;

    for (i128 c : Q.coeff) {
        if (c != 0) {
            if (seen) {
                return false;
            }

            seen = true;
        }
    }

    return seen;
}

bool is_binomial(
    const ReducedPolynomial &Q
) {
    u64 count = 0;

    for (i128 c : Q.coeff) {
        if (c != 0) {
            ++count;
        }
    }

    return count == 2;
}

u64 coefficient_gcd_i128(
    const ReducedPolynomial &Q
) {
    u64 g = 0;

    for (i128 c : Q.coeff) {
        if (c < 0) {
            c = -c;
        }

        u64 v =
            static_cast<u64>(
                c
            );

        g =
            gcd_u64(
                g,
                v
            );
    }

    return g;
}

void print_reduced(
    const ReducedPolynomial &Q
) {
    std::cout
        << "Q_coefficients=[";

    for (
        std::size_t i = 0;
        i < Q.coeff.size();
        ++i
    ) {
        if (i != 0) {
            std::cout << ",";
        }

        std::cout
            << to_string_i128(
                Q.coeff[i]
            );
    }

    std::cout
        << "]\n";
}

void accumulate_stats(
    Stats &dst,
    const Stats &src
) {
    dst.q_zero += src.q_zero;
    dst.q_one_zero += src.q_one_zero;
    dst.q_minus_one_zero +=
        src.q_minus_one_zero;

    dst.q_two_zero +=
        src.q_two_zero;

    dst.q_minus_two_zero +=
        src.q_minus_two_zero;

    dst.q_mod_pattern_equal +=
        src.q_mod_pattern_equal;

    dst.monomial +=
        src.monomial;

    dst.binomial +=
        src.binomial;

    dst.degree1 +=
        src.degree1;

    dst.degree2 +=
        src.degree2;

    dst.degree3plus +=
        src.degree3plus;
}

int main() {
    std::cout
        << "START EXPERIMENT 306\n"
        << "REMOVE THE FORCED t FACTOR\n"
        << "WHAT STRUCTURE REMAINS IN Q(t)=P(t)/t?\n"
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

    u64 total_polynomials = 0;
    u64 total_zero_polynomials = 0;
    u64 total_nonzero_polynomials = 0;

    u64 total_q_degree_sum = 0;
    u64 min_q_degree = UINT64_MAX;
    u64 max_q_degree = 0;

    u64 printed_interesting = 0;

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
            << " p=" << p
            << " q=" << q
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

                    Polynomial P =
                        build_polynomial(
                            m,
                            y,
                            base
                        );

                    ++case_polynomials;
                    ++total_polynomials;

                    u64 g =
                        coefficient_gcd(P);

                    if (g == 0) {
                        ++case_zero;
                        ++total_zero_polynomials;
                        continue;
                    }

                    ++case_nonzero;
                    ++total_nonzero_polynomials;

                    PrimitivePolynomial R =
                        primitive_part(
                            P,
                            g
                        );

                    ReducedPolynomial Q =
                        divide_forced_t(R);

                    std::size_t d =
                        degree(Q);

                    total_q_degree_sum +=
                        static_cast<u64>(d);

                    min_q_degree =
                        std::min(
                            min_q_degree,
                            static_cast<u64>(d)
                        );

                    max_q_degree =
                        std::max(
                            max_q_degree,
                            static_cast<u64>(d)
                        );

                    bool q_zero =
                        is_zero(Q);

                    i128 q0 =
                        evaluate(Q, 0);

                    i128 q1 =
                        evaluate(Q, 1);

                    i128 qm1 =
                        evaluate(Q, -1);

                    i128 q2 =
                        evaluate(Q, 2);

                    i128 qm2 =
                        evaluate(Q, -2);

                    Stats current;

                    if (q_zero) {
                        current.q_zero = 1;
                    }

                    if (q1 == 0) {
                        current.q_one_zero = 1;
                    }

                    if (qm1 == 0) {
                        current.q_minus_one_zero = 1;
                    }

                    if (q2 == 0) {
                        current.q_two_zero = 1;
                    }

                    if (qm2 == 0) {
                        current.q_minus_two_zero = 1;
                    }

                    /*
                        Compare Q(1) with the known scalar
                        potential. Since P=tQ, we have:

                            P(1)=Q(1).

                        This is a basic consistency relation.
                    */
                    u64 miss =
                        miss_prefix(
                            m,
                            y,
                            base
                        );

                    i128 expected =
                        static_cast<i128>(
                            miss
                        ) - 1;

                    if (
                        q1 ==
                        expected
                    ) {
                        current.q_mod_pattern_equal = 1;
                    }

                    if (is_monomial(Q)) {
                        current.monomial = 1;
                    }

                    if (is_binomial(Q)) {
                        current.binomial = 1;
                    }

                    if (d == 1) {
                        current.degree1 = 1;
                    } else if (d == 2) {
                        current.degree2 = 1;
                    } else if (d >= 3) {
                        current.degree3plus = 1;
                    }

                    accumulate_stats(
                        case_stats,
                        current
                    );

                    accumulate_stats(
                        total,
                        current
                    );

                    /*
                        Print unusual nonzero Q structures.
                    */
                    bool interesting =
                        current.q_one_zero ||
                        current.q_minus_one_zero ||
                        current.q_two_zero ||
                        current.q_minus_two_zero ||
                        current.binomial;

                    if (
                        interesting &&
                        printed_interesting < 30
                    ) {
                        ++printed_interesting;

                        std::cout
                            << "INTERESTING_Q\n"
                            << "base="
                            << base
                            << " m="
                            << m
                            << " n="
                            << n
                            << "\n"
                            << "content_gcd="
                            << g
                            << "\n"
                            << "degree_Q="
                            << d
                            << "\n"
                            << "Q(0)="
                            << to_string_i128(q0)
                            << "\n"
                            << "Q(1)="
                            << to_string_i128(q1)
                            << "\n"
                            << "Q(-1)="
                            << to_string_i128(qm1)
                            << "\n"
                            << "Q(2)="
                            << to_string_i128(q2)
                            << "\n"
                            << "Q(-2)="
                            << to_string_i128(qm2)
                            << "\n"
                            << "monomial="
                            << (current.monomial ? 1 : 0)
                            << "\n"
                            << "binomial="
                            << (current.binomial ? 1 : 0)
                            << "\n";

                        print_reduced(Q);
                    }
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
            << "Q(0)=0="
            << case_stats.q_zero
            << "\n";

        std::cout
            << "Q(1)=0="
            << case_stats.q_one_zero
            << "\n";

        std::cout
            << "Q(-1)=0="
            << case_stats.q_minus_one_zero
            << "\n";

        std::cout
            << "Q(2)=0="
            << case_stats.q_two_zero
            << "\n";

        std::cout
            << "Q(-2)=0="
            << case_stats.q_minus_two_zero
            << "\n";

        std::cout
            << "Q(1)_matches_potential="
            << case_stats.q_mod_pattern_equal
            << "\n";

        std::cout
            << "monomial="
            << case_stats.monomial
            << "\n";

        std::cout
            << "binomial="
            << case_stats.binomial
            << "\n";

        std::cout
            << "degree1="
            << case_stats.degree1
            << "\n";

        std::cout
            << "degree2="
            << case_stats.degree2
            << "\n";

        std::cout
            << "degree3plus="
            << case_stats.degree3plus
            << "\n\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "polynomials="
        << total_polynomials
        << "\n";

    std::cout
        << "zero_polynomials="
        << total_zero_polynomials
        << "\n";

    std::cout
        << "nonzero_polynomials="
        << total_nonzero_polynomials
        << "\n";

    if (total_nonzero_polynomials > 0) {
        std::cout
            << "average_Q_degree="
            << static_cast<double>(
                total_q_degree_sum
            ) /
               static_cast<double>(
                   total_nonzero_polynomials
               )
            << "\n";
    }

    std::cout
        << "min_Q_degree="
        << min_q_degree
        << "\n";

    std::cout
        << "max_Q_degree="
        << max_q_degree
        << "\n";

    std::cout
        << "Q(0)=0="
        << total.q_zero
        << "\n";

    std::cout
        << "Q(1)=0="
        << total.q_one_zero
        << "\n";

    std::cout
        << "Q(-1)=0="
        << total.q_minus_one_zero
        << "\n";

    std::cout
        << "Q(2)=0="
        << total.q_two_zero
        << "\n";

    std::cout
        << "Q(-2)=0="
        << total.q_minus_two_zero
        << "\n";

    std::cout
        << "Q(1)_matches_potential="
        << total.q_mod_pattern_equal
        << "\n";

    std::cout
        << "monomial="
        << total.monomial
        << "\n";

    std::cout
        << "binomial="
        << total.binomial
        << "\n";

    std::cout
        << "degree1="
        << total.degree1
        << "\n";

    std::cout
        << "degree2="
        << total.degree2
        << "\n";

    std::cout
        << "degree3plus="
        << total.degree3plus
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 306\n";

    return 0;
}
