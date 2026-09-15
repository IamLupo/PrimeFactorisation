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
    std::vector<i128> coeff;
    u64 steps = 0;
};

struct Stats {
    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;
    u64 p_only = 0;
    u64 q_only = 0;
    u64 both = 0;
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
        int d = static_cast<int>(v % 10);
        s.push_back(static_cast<char>('0' + d));
        v /= 10;
    }

    if (neg) {
        s.push_back('-');
    }

    std::reverse(s.begin(), s.end());

    return s;
}

u64 abs_u64(i128 x) {
    u128 v =
        x < 0
            ? static_cast<u128>(-x)
            : static_cast<u128>(x);

    return static_cast<u64>(v);
}

u64 integer_sqrt(u64 n) {
    if (n == 0) {
        return 0;
    }

    u64 lo = 0;
    u64 hi = std::min<u64>(n, 1ULL << 32);

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
        d.push_back(x % base);
        x /= base;
    }

    return d;
}

/*
    M_m(y) =
        #{x : 0 <= x <= y and x <=_p m}

    Independent MSB-first digit DP.
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
        static_cast<u64>(y_signed);

    std::vector<u64> md =
        base_digits(m, base);

    std::vector<u64> yd =
        base_digits(y, base);

    std::size_t L =
        std::max(
            md.size(),
            yd.size()
        );

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

    u128 ans = 0;
    bool tight = true;

    for (std::size_t pos = L; pos-- > 0;) {
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

        for (u64 xi = 0;
             xi < upper;
             ++xi) {

            ans += suffix[pos + 1];
        }

        if (yi <= mi) {
            tight = true;
        } else {
            tight = false;
            break;
        }
    }

    if (tight) {
        ++ans;
    }

    return static_cast<u64>(ans);
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
    std::size_t r
) {
    u64 v = 1;

    for (std::size_t i = 0;
         i < r;
         ++i) {
        v *= base;
    }

    return v;
}

/*
    Build

        P(t) = sum Delta_r * t^r

    where Delta_r is the aggregate derivative
    of the canonical path at digit position r.

    We store one coefficient per digit position.
*/
Polynomial canonical_polynomial(
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

    for (std::size_t r = 0;
         r < digits.size();
         ++r) {

        u64 step =
            base_power(
                base,
                r
            );

        for (u64 t = 0;
             t < digits[r];
             ++t) {

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
                static_cast<i128>(delta);

            ++P.steps;

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

i128 evaluate_polynomial(
    const Polynomial &P,
    i128 t
) {
    /*
        Horner evaluation.
    */
    i128 result = 0;

    for (std::size_t i = P.coeff.size();
         i-- > 0;) {

        result =
            result * t +
            P.coeff[i];
    }

    return result;
}

u64 gcd_u64(
    u64 a,
    u64 b
) {
    return std::gcd(a, b);
}

Stats classify(
    i128 value,
    u64 N,
    u64 p,
    u64 q
) {
    Stats s;

    u64 g =
        gcd_u64(
            abs_u64(value),
            N
        );

    if (g == 1) {
        ++s.gcd1;
        return s;
    }

    if (g == N) {
        ++s.gcdN;
        return s;
    }

    ++s.nontrivial;

    if (g == p && g == q) {
        ++s.both;
    } else if (g == p) {
        ++s.p_only;
    } else if (g == q) {
        ++s.q_only;
    }

    return s;
}

void accumulate(
    Stats &dst,
    const Stats &src
) {
    dst.gcd1 += src.gcd1;
    dst.gcdN += src.gcdN;
    dst.nontrivial += src.nontrivial;
    dst.p_only += src.p_only;
    dst.q_only += src.q_only;
    dst.both += src.both;
}

void print_stats(
    const std::string &name,
    const Stats &s
) {
    std::cout
        << name
        << " gcd1=" << s.gcd1
        << " gcdN=" << s.gcdN
        << " nontrivial=" << s.nontrivial
        << " p_only=" << s.p_only
        << " q_only=" << s.q_only
        << " both=" << s.both
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 301\n"
        << "PATH DERIVATIVE GENERATING POLYNOMIAL\n"
        << "CAN EVALUATION AWAY FROM t=1 EXPOSE A HIDDEN FACTOR?\n"
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
        Evaluation points.

        t = 1 reproduces S0.

        t = 0 isolates the lowest digit-position
        contribution.

        Negative t tests cancellation patterns.

        Larger positive t gives progressively stronger
        weighting toward high digit positions.
    */
    const std::vector<int> eval_points = {
        -8, -7, -6, -5, -4,
        -3, -2, -1,
         0,
         1,
         2, 3, 4, 5, 6, 7, 8,
        10, 12, 16, 20
    };

    Stats total;

    u64 total_polynomials = 0;
    u64 total_evaluations = 0;
    u64 total_nontrivial = 0;

    u64 max_abs_value = 0;

    bool printed_first_hit = false;

    for (std::size_t ci = 0;
         ci < cases.size();
         ++ci) {

        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        u64 case_polynomials = 0;
        u64 case_evaluations = 0;
        u64 case_nontrivial = 0;

        std::cout
            << "CASE "
            << ci
            << " p=" << p
            << " q=" << q
            << "\n";

        u64 root =
            integer_sqrt(N);

        for (u64 base : bases) {
            for (int off : offsets) {

                i128 m128 =
                    static_cast<i128>(root) +
                    static_cast<i128>(off);

                if (m128 <= 0) {
                    continue;
                }

                u64 m =
                    static_cast<u64>(m128);

                std::vector<u64> ns;

                ns.push_back(N);
                ns.push_back(m);
                ns.push_back(m + 1);

                if (m != 0 &&
                    m <= N / m) {

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
                        canonical_polynomial(
                            m,
                            y,
                            base
                        );

                    ++case_polynomials;
                    ++total_polynomials;

                    for (int t_int :
                         eval_points) {

                        i128 t =
                            static_cast<i128>(t_int);

                        i128 value =
                            evaluate_polynomial(
                                P,
                                t
                            );

                        ++case_evaluations;
                        ++total_evaluations;

                        max_abs_value =
                            std::max(
                                max_abs_value,
                                abs_u64(value)
                            );

                        Stats s =
                            classify(
                                value,
                                N,
                                p,
                                q
                            );

                        accumulate(
                            total,
                            s
                        );

                        if (s.nontrivial != 0) {
                            ++case_nontrivial;
                            ++total_nontrivial;

                            if (!printed_first_hit) {
                                printed_first_hit =
                                    true;

                                u64 g =
                                    gcd_u64(
                                        abs_u64(value),
                                        N
                                    );

                                std::cout
                                    << "FIRST_GLOBAL_HIT\n"
                                    << "case=" << ci
                                    << " p=" << p
                                    << " q=" << q
                                    << "\n"
                                    << "base=" << base
                                    << " m=" << m
                                    << " n=" << n
                                    << "\n"
                                    << "t=" << t_int
                                    << "\n"
                                    << "value="
                                    << to_string_i128(
                                        value
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

        std::cout
            << "polynomials="
            << case_polynomials
            << "\n";

        std::cout
            << "evaluations="
            << case_evaluations
            << "\n";

        std::cout
            << "nontrivial="
            << case_nontrivial
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
        << "evaluations="
        << total_evaluations
        << "\n";

    std::cout
        << "nontrivial="
        << total_nontrivial
        << "\n";

    print_stats(
        "ALL",
        total
    );

    std::cout
        << "max_abs_value="
        << max_abs_value
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 301\n";

    return 0;
}
