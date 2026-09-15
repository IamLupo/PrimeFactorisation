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

struct FactorStats {
    u64 factor_t = 0;
    u64 factor_t_minus_1 = 0;
    u64 factor_t_plus_1 = 0;
    u64 factor_t_minus_2 = 0;
    u64 factor_t_plus_2 = 0;
    u64 factor_t2_minus_1 = 0;
    u64 factor_t2_plus_1 = 0;
};

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool negative = (x < 0);

    u128 v = negative
        ? static_cast<u128>(-x)
        : static_cast<u128>(x);

    std::string s;

    while (v > 0) {
        int digit = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + digit));
        v /= 10;
    }

    if (negative) {
        s.push_back('-');
    }

    std::reverse(s.begin(), s.end());

    return s;
}

u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

u64 integer_sqrt(u64 n) {
    if (n == 0) {
        return 0;
    }

    u64 lo = 0;
    u64 hi = std::min<u64>(n, 1ULL << 32);

    while (lo + 1 < hi) {
        u64 mid = lo + (hi - lo) / 2;

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
        digits.push_back(x % base);
        x /= base;
    }

    return digits;
}

/*
    M_m(y) =
        #{ x : 0 <= x <= y and x <=_p m }
*/
u64 miss_prefix(
    u64 m,
    i128 y_signed,
    u64 base
) {
    if (y_signed < 0) {
        return 0;
    }

    u64 y = static_cast<u64>(y_signed);

    std::vector<u64> md =
        base_digits(m, base);

    std::vector<u64> yd =
        base_digits(y, base);

    std::size_t L =
        std::max(md.size(), yd.size());

    md.resize(L, 0);
    yd.resize(L, 0);

    std::vector<u128> suffix(
        L + 1,
        1
    );

    for (std::size_t i = L; i-- > 0;) {
        suffix[i] =
            suffix[i + 1] *
            static_cast<u128>(md[i] + 1);
    }

    u128 answer = 0;
    bool tight = true;

    for (std::size_t pos = L; pos-- > 0;) {
        if (!tight) {
            break;
        }

        u64 yi = yd[pos];
        u64 mi = md[pos];

        u64 upper =
            std::min(yi, mi + 1);

        for (u64 xi = 0; xi < upper; ++xi) {
            answer += suffix[pos + 1];
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

    return static_cast<u64>(answer);
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
            << "ERROR: negative local derivative\n"
            << "old_m=" << old_m
            << " new_m=" << new_m
            << " y=" << to_string_i128(y)
            << " base=" << base
            << "\n";

        std::exit(1);
    }

    return b - a;
}

u64 base_power(
    u64 base,
    std::size_t exponent
) {
    u64 result = 1;

    for (std::size_t i = 0; i < exponent; ++i) {
        result *= base;
    }

    return result;
}

/*
    P(t) = sum_r Delta_r t^r
*/
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

            P.coeff[r] += delta;

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

u64 coefficient_gcd(
    const Polynomial &P
) {
    u64 g = 0;

    for (u64 c : P.coeff) {
        g = gcd_u64(g, c);
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

i128 evaluate(
    const PrimitivePolynomial &P,
    i128 x
) {
    i128 result = 0;

    for (
        std::size_t i = P.coeff.size();
        i-- > 0;
    ) {
        result =
            result * x +
            P.coeff[i];
    }

    return result;
}

std::size_t degree(
    const PrimitivePolynomial &P
) {
    if (P.coeff.empty()) {
        return 0;
    }

    return P.coeff.size() - 1;
}

bool divisible_by_t(
    const PrimitivePolynomial &P
) {
    return !P.coeff.empty() &&
           P.coeff[0] == 0;
}

bool divisible_by_t_minus_1(
    const PrimitivePolynomial &P
) {
    return evaluate(P, 1) == 0;
}

bool divisible_by_t_plus_1(
    const PrimitivePolynomial &P
) {
    return evaluate(P, -1) == 0;
}

bool divisible_by_t_minus_2(
    const PrimitivePolynomial &P
) {
    return evaluate(P, 2) == 0;
}

bool divisible_by_t_plus_2(
    const PrimitivePolynomial &P
) {
    return evaluate(P, -2) == 0;
}

/*
    t^2 - 1 = (t-1)(t+1)
*/
bool divisible_by_t2_minus_1(
    const PrimitivePolynomial &P
) {
    return
        divisible_by_t_minus_1(P) &&
        divisible_by_t_plus_1(P);
}

/*
    P(i) = 0 iff its real and imaginary parts vanish.

    For

        P(t) = a_0 + a_1 t + a_2 t^2 + ...

    we have

        P(i) = E(-1) + i O(-1).
*/
bool divisible_by_t2_plus_1(
    const PrimitivePolynomial &P
) {
    i128 real = 0;
    i128 imag = 0;

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if ((i & 1U) == 0) {
            std::size_t k = i / 2;

            i128 sign =
                (k & 1U)
                    ? -1
                    : 1;

            real +=
                sign *
                P.coeff[i];
        } else {
            std::size_t k = (i - 1) / 2;

            i128 sign =
                (k & 1U)
                    ? -1
                    : 1;

            imag +=
                sign *
                P.coeff[i];
        }
    }

    return real == 0 &&
           imag == 0;
}

void print_coefficients(
    const PrimitivePolynomial &P
) {
    std::cout
        << "primitive_coefficients=[";

    for (
        std::size_t i = 0;
        i < P.coeff.size();
        ++i
    ) {
        if (i != 0) {
            std::cout << ",";
        }

        std::cout
            << to_string_i128(
                P.coeff[i]
            );
    }

    std::cout
        << "]\n";
}

void accumulate_stats(
    FactorStats &dst,
    const FactorStats &src
) {
    dst.factor_t += src.factor_t;
    dst.factor_t_minus_1 += src.factor_t_minus_1;
    dst.factor_t_plus_1 += src.factor_t_plus_1;
    dst.factor_t_minus_2 += src.factor_t_minus_2;
    dst.factor_t_plus_2 += src.factor_t_plus_2;
    dst.factor_t2_minus_1 += src.factor_t2_minus_1;
    dst.factor_t2_plus_1 += src.factor_t2_plus_1;
}

int main() {
    std::cout
        << "START EXPERIMENT 305\n"
        << "PRIMITIVE PATH POLYNOMIAL STRUCTURE\n"
        << "DO THE NORMALIZED POLYNOMIALS HAVE UNIVERSAL FACTORS?\n"
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

    FactorStats total_stats;

    u64 total_polynomials = 0;
    u64 total_zero = 0;
    u64 total_nonzero = 0;

    u64 degree_sum = 0;
    u64 min_degree = UINT64_MAX;
    u64 max_degree = 0;

    u64 polynomials_with_simple_factor = 0;

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

        FactorStats case_stats;

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
                    static_cast<u64>(m_signed);

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
                        ++total_zero;
                        continue;
                    }

                    ++case_nonzero;
                    ++total_nonzero;

                    PrimitivePolynomial R =
                        primitive_part(
                            P,
                            g
                        );

                    u64 d =
                        static_cast<u64>(
                            degree(R)
                        );

                    degree_sum += d;

                    min_degree =
                        std::min(
                            min_degree,
                            d
                        );

                    max_degree =
                        std::max(
                            max_degree,
                            d
                        );

                    bool f_t =
                        divisible_by_t(R);

                    bool f_t_minus_1 =
                        divisible_by_t_minus_1(R);

                    bool f_t_plus_1 =
                        divisible_by_t_plus_1(R);

                    bool f_t_minus_2 =
                        divisible_by_t_minus_2(R);

                    bool f_t_plus_2 =
                        divisible_by_t_plus_2(R);

                    bool f_t2_minus_1 =
                        divisible_by_t2_minus_1(R);

                    bool f_t2_plus_1 =
                        divisible_by_t2_plus_1(R);

                    FactorStats current;

                    if (f_t) {
                        current.factor_t = 1;
                    }

                    if (f_t_minus_1) {
                        current.factor_t_minus_1 = 1;
                    }

                    if (f_t_plus_1) {
                        current.factor_t_plus_1 = 1;
                    }

                    if (f_t_minus_2) {
                        current.factor_t_minus_2 = 1;
                    }

                    if (f_t_plus_2) {
                        current.factor_t_plus_2 = 1;
                    }

                    if (f_t2_minus_1) {
                        current.factor_t2_minus_1 = 1;
                    }

                    if (f_t2_plus_1) {
                        current.factor_t2_plus_1 = 1;
                    }

                    accumulate_stats(
                        case_stats,
                        current
                    );

                    accumulate_stats(
                        total_stats,
                        current
                    );

                    bool simple_factor =
                        f_t ||
                        f_t_minus_1 ||
                        f_t_plus_1 ||
                        f_t_minus_2 ||
                        f_t_plus_2 ||
                        f_t2_plus_1;

                    if (simple_factor) {
                        ++polynomials_with_simple_factor;

                        if (printed_interesting < 20) {
                            ++printed_interesting;

                            std::cout
                                << "INTERESTING_POLYNOMIAL\n"
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
                                << "degree="
                                << d
                                << "\n"
                                << "P(0)="
                                << to_string_i128(
                                    evaluate(R, 0)
                                )
                                << "\n"
                                << "P(1)="
                                << to_string_i128(
                                    evaluate(R, 1)
                                )
                                << "\n"
                                << "P(-1)="
                                << to_string_i128(
                                    evaluate(R, -1)
                                )
                                << "\n"
                                << "P(2)="
                                << to_string_i128(
                                    evaluate(R, 2)
                                )
                                << "\n"
                                << "P(-2)="
                                << to_string_i128(
                                    evaluate(R, -2)
                                )
                                << "\n"
                                << "factor_t="
                                << (f_t ? 1 : 0)
                                << "\n"
                                << "factor_t_minus_1="
                                << (f_t_minus_1 ? 1 : 0)
                                << "\n"
                                << "factor_t_plus_1="
                                << (f_t_plus_1 ? 1 : 0)
                                << "\n"
                                << "factor_t_minus_2="
                                << (f_t_minus_2 ? 1 : 0)
                                << "\n"
                                << "factor_t_plus_2="
                                << (f_t_plus_2 ? 1 : 0)
                                << "\n"
                                << "factor_t2_minus_1="
                                << (f_t2_minus_1 ? 1 : 0)
                                << "\n"
                                << "factor_t2_plus_1="
                                << (f_t2_plus_1 ? 1 : 0)
                                << "\n";

                            print_coefficients(R);
                        }
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
            << "factor_t="
            << case_stats.factor_t
            << "\n";

        std::cout
            << "factor_t_minus_1="
            << case_stats.factor_t_minus_1
            << "\n";

        std::cout
            << "factor_t_plus_1="
            << case_stats.factor_t_plus_1
            << "\n";

        std::cout
            << "factor_t_minus_2="
            << case_stats.factor_t_minus_2
            << "\n";

        std::cout
            << "factor_t_plus_2="
            << case_stats.factor_t_plus_2
            << "\n";

        std::cout
            << "factor_t2_minus_1="
            << case_stats.factor_t2_minus_1
            << "\n";

        std::cout
            << "factor_t2_plus_1="
            << case_stats.factor_t2_plus_1
            << "\n";

        std::cout
            << "\n";
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
        << total_zero
        << "\n";

    std::cout
        << "nonzero_polynomials="
        << total_nonzero
        << "\n";

    if (total_nonzero > 0) {
        std::cout
            << "average_degree="
            << static_cast<double>(degree_sum) /
               static_cast<double>(total_nonzero)
            << "\n";
    }

    std::cout
        << "min_degree="
        << min_degree
        << "\n";

    std::cout
        << "max_degree="
        << max_degree
        << "\n";

    std::cout
        << "factor_t="
        << total_stats.factor_t
        << "\n";

    std::cout
        << "factor_t_minus_1="
        << total_stats.factor_t_minus_1
        << "\n";

    std::cout
        << "factor_t_plus_1="
        << total_stats.factor_t_plus_1
        << "\n";

    std::cout
        << "factor_t_minus_2="
        << total_stats.factor_t_minus_2
        << "\n";

    std::cout
        << "factor_t_plus_2="
        << total_stats.factor_t_plus_2
        << "\n";

    std::cout
        << "factor_t2_minus_1="
        << total_stats.factor_t2_minus_1
        << "\n";

    std::cout
        << "factor_t2_plus_1="
        << total_stats.factor_t2_plus_1
        << "\n";

    std::cout
        << "polynomials_with_simple_factor="
        << polynomials_with_simple_factor
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 305\n";

    return 0;
}