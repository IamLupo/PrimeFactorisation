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

struct Moments {
    i128 S[4] = {0, 0, 0, 0};
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

struct Coeff {
    int a[4];
};

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool neg = (x < 0);

    u128 v = neg
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

std::vector<u64> base_digits(u64 x, u64 base) {
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
        #{ x : 0 <= x <= y and x <=_p m }

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

    u128 ans = 0;
    bool tight = true;

    for (std::size_t pos = L; pos-- > 0;) {
        if (!tight) {
            break;
        }

        u64 yi = yd[pos];
        u64 mi = md[pos];

        u64 upper =
            std::min(yi, mi + 1);

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

    for (std::size_t i = 0; i < r; ++i) {
        v *= base;
    }

    return v;
}

/*
    S_k = sum (r+1)^k * Delta_r
    for k = 0,1,2,3.
*/
Moments canonical_path_moments(
    u64 target_m,
    i128 y,
    u64 base
) {
    std::vector<u64> digits =
        base_digits(target_m, base);

    Moments out;

    u64 current = 0;

    for (std::size_t r = 0;
         r < digits.size();
         ++r) {

        u64 step =
            base_power(base, r);

        i128 w1 =
            static_cast<i128>(r + 1);

        i128 w2 = w1 * w1;
        i128 w3 = w2 * w1;

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

            i128 d =
                static_cast<i128>(delta);

            out.S[0] += d;
            out.S[1] += w1 * d;
            out.S[2] += w2 * d;
            out.S[3] += w3 * d;

            ++out.steps;

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

    return out;
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

bool is_zero_coeff(
    const Coeff &c
) {
    return
        c.a[0] == 0 &&
        c.a[1] == 0 &&
        c.a[2] == 0 &&
        c.a[3] == 0;
}

i128 evaluate(
    const Moments &M,
    const Coeff &c
) {
    i128 v = 0;

    for (int k = 0; k < 4; ++k) {
        v +=
            static_cast<i128>(c.a[k]) *
            M.S[k];
    }

    return v;
}

void print_coeff(
    const Coeff &c
) {
    std::cout
        << "("
        << c.a[0] << ","
        << c.a[1] << ","
        << c.a[2] << ","
        << c.a[3]
        << ")";
}

int main() {
    std::cout
        << "START EXPERIMENT 300\n"
        << "LINEAR SPAN OF PATH-INTEGRATED DIGIT DERIVATIVES\n"
        << "CAN ANY SMALL INTEGER COMBINATION EXPOSE A HIDDEN FACTOR?\n"
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
         1,  2,  3,  4
    };

    /*
        All coefficient vectors in [-2,2]^4.
    */
    std::vector<Coeff> coefficients;

    for (int a0 = -2; a0 <= 2; ++a0) {
        for (int a1 = -2; a1 <= 2; ++a1) {
            for (int a2 = -2; a2 <= 2; ++a2) {
                for (int a3 = -2; a3 <= 2; ++a3) {
                    Coeff c = {
                        {a0, a1, a2, a3}
                    };

                    if (!is_zero_coeff(c)) {
                        coefficients.push_back(c);
                    }
                }
            }
        }
    }

    std::cout
        << "coefficient_vectors="
        << coefficients.size()
        << "\n\n";

    Stats total;
    u64 total_moment_cases = 0;
    u64 total_linear_tests = 0;
    u64 total_nontrivial = 0;

    u64 first_hit_count = 0;

    bool printed_first_hit = false;

    for (std::size_t ci = 0;
         ci < cases.size();
         ++ci) {

        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        u64 case_moment_cases = 0;
        u64 case_linear_tests = 0;
        u64 case_nontrivial = 0;

        bool case_hit = false;

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

                if (m != 0 && m <= N / m) {
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

                    Moments M =
                        canonical_path_moments(
                            m,
                            y,
                            base
                        );

                    ++case_moment_cases;
                    ++total_moment_cases;

                    for (const Coeff &c :
                         coefficients) {

                        i128 value =
                            evaluate(M, c);

                        Stats s =
                            classify(
                                value,
                                N,
                                p,
                                q
                            );

                        accumulate(total, s);

                        ++case_linear_tests;
                        ++total_linear_tests;

                        if (s.nontrivial != 0) {
                            ++case_nontrivial;
                            ++total_nontrivial;

                            if (!case_hit) {
                                case_hit = true;
                                ++first_hit_count;

                                std::cout
                                    << "FIRST_CASE_HIT\n"
                                    << "base=" << base
                                    << " m=" << m
                                    << " n=" << n
                                    << "\n"
                                    << "coeff=";

                                print_coeff(c);

                                u64 g =
                                    gcd_u64(
                                        abs_u64(value),
                                        N
                                    );

                                std::cout
                                    << "\nvalue="
                                    << to_string_i128(value)
                                    << "\ngcd="
                                    << g
                                    << "\n";
                            }

                            if (!printed_first_hit) {
                                printed_first_hit = true;

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
                                    << "coeff=";

                                print_coeff(c);

                                u64 g =
                                    gcd_u64(
                                        abs_u64(value),
                                        N
                                    );

                                std::cout
                                    << "\nvalue="
                                    << to_string_i128(value)
                                    << "\ngcd="
                                    << g
                                    << "\n";
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "moment_cases="
            << case_moment_cases
            << "\n";

        std::cout
            << "linear_tests="
            << case_linear_tests
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
        << "moment_cases="
        << total_moment_cases
        << "\n";

    std::cout
        << "linear_tests="
        << total_linear_tests
        << "\n";

    std::cout
        << "nontrivial="
        << total_nontrivial
        << "\n";

    std::cout
        << "cases_with_hit="
        << first_hit_count
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
        << "both="
        << total.both
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 300\n";

    return 0;
}
