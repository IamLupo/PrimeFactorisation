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

struct RootInfo {
    u64 roots_p = 0;
    u64 roots_q = 0;
    u64 common_roots = 0;
    u64 p_only = 0;
    u64 q_only = 0;
};

struct DivClass {
    u64 neither = 0;
    u64 p_only = 0;
    u64 q_only = 0;
    u64 both = 0;
};

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }

    std::string s;

    while (x > 0) {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    return s;
}

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool negative = x < 0;

    u128 v = negative
        ? static_cast<u128>(-x)
        : static_cast<u128>(x);

    std::string s = to_string_u128(v);

    if (negative) {
        s.insert(s.begin(), '-');
    }

    return s;
}

u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

u64 abs_i128_to_u64(i128 x) {
    u128 v = x < 0
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
    Independent reference:

        M_m(y) =
        #{x : 0 <= x <= y and x <=_p m}
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
            << "ERROR: negative derivative\n"
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

    for (std::size_t i = 0;
         i < exponent;
         ++i) {

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
        base_digits(target_m, base);

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
            base_power(base, r);

        for (u64 k = 0;
             k < digits[r];
             ++k) {

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
        g = gcd_u64(g, c);
    }

    return g;
}

bool all_coeff_divisible(
    const Polynomial &P,
    u64 mod
) {
    for (u64 c : P.coeff) {
        if (c % mod != 0) {
            return false;
        }
    }

    return true;
}

std::size_t polynomial_degree(
    const Polynomial &P
) {
    for (std::size_t i = P.coeff.size();
         i-- > 0;) {

        if (P.coeff[i] != 0) {
            return i;
        }
    }

    return 0;
}

u64 mod_add(
    u64 a,
    u64 b,
    u64 mod
) {
    u64 r = a + b;

    if (r >= mod || r < a) {
        r %= mod;
    }

    return r;
}

u64 mod_mul(
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

u64 evaluate_mod(
    const Polynomial &P,
    u64 x,
    u64 mod
) {
    u64 result = 0;

    for (std::size_t i = P.coeff.size();
         i-- > 0;) {

        result =
            mod_mul(
                result,
                x,
                mod
            );

        result =
            mod_add(
                result,
                P.coeff[i] % mod,
                mod
            );
    }

    return result;
}

RootInfo scan_roots(
    const Polynomial &P,
    u64 p,
    u64 q
) {
    RootInfo R;

    std::vector<u64> p_roots;
    std::vector<u64> q_roots;

    /*
        Only scan a field if the polynomial is not
        identically zero modulo that field.
    */

    if (!all_coeff_divisible(P, p)) {
        for (u64 x = 0; x < p; ++x) {
            if (evaluate_mod(P, x, p) == 0) {
                p_roots.push_back(x);
            }
        }
    } else {
        R.roots_p = p;
    }

    if (!all_coeff_divisible(P, q)) {
        for (u64 x = 0; x < q; ++x) {
            if (evaluate_mod(P, x, q) == 0) {
                q_roots.push_back(x);
            }
        }
    } else {
        R.roots_q = q;
    }

    if (!p_roots.empty()) {
        R.roots_p =
            static_cast<u64>(
                p_roots.size()
            );
    }

    if (!q_roots.empty()) {
        R.roots_q =
            static_cast<u64>(
                q_roots.size()
            );
    }

    std::size_t i = 0;
    std::size_t j = 0;

    while (
        i < p_roots.size() &&
        j < q_roots.size()
    ) {
        if (p_roots[i] == q_roots[j]) {
            ++R.common_roots;
            ++i;
            ++j;
        } else if (
            p_roots[i] <
            q_roots[j]
        ) {
            ++i;
        } else {
            ++j;
        }
    }

    if (
        R.roots_p == p &&
        R.roots_q == q
    ) {
        R.common_roots = p;
    } else if (
        R.roots_p == p
    ) {
        /*
            If P is identically zero mod p,
            every x in 0..p-1 is a root modulo p.
            Common roots are the roots found modulo q
            among those integers.
        */

        if (!q_roots.empty()) {
            R.common_roots =
                static_cast<u64>(
                    std::min<std::size_t>(
                        p_roots.size(),
                        q_roots.size()
                    )
                );
        }
    } else if (
        R.roots_q == q
    ) {
        /*
            Symmetric situation.
        */
        if (!p_roots.empty()) {
            R.common_roots =
                static_cast<u64>(
                    std::min<std::size_t>(
                        p_roots.size(),
                        q_roots.size()
                    )
                );
        }
    }

    R.p_only =
        R.roots_p -
        R.common_roots;

    R.q_only =
        R.roots_q -
        R.common_roots;

    return R;
}

void accumulate_class(
    DivClass &dst,
    bool div_p,
    bool div_q
) {
    if (div_p && div_q) {
        ++dst.both;
    } else if (div_p) {
        ++dst.p_only;
    } else if (div_q) {
        ++dst.q_only;
    } else {
        ++dst.neither;
    }
}

void print_polynomial(
    const Polynomial &P
) {
    std::cout << "coefficients=[";

    for (std::size_t i = 0;
         i < P.coeff.size();
         ++i) {

        if (i != 0) {
            std::cout << ",";
        }

        std::cout << P.coeff[i];
    }

    std::cout << "]\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 303\n"
        << "COEFFICIENT DIVISIBILITY OF PATH POLYNOMIALS\n"
        << "ARE THE ROOTS CAUSED BY WHOLE-POLYNOMIAL FACTOR DIVISIBILITY?\n"
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

    u64 total_polynomials = 0;

    u64 total_identical_p = 0;
    u64 total_identical_q = 0;
    u64 total_identical_N = 0;

    u64 total_p_only_coeff = 0;
    u64 total_q_only_coeff = 0;
    u64 total_neither_coeff = 0;

    u64 total_nontrivial_coeff_gcd = 0;

    u64 total_root_p = 0;
    u64 total_root_q = 0;
    u64 total_root_common = 0;
    u64 total_root_p_only = 0;
    u64 total_root_q_only = 0;

    u64 max_coefficient_gcd = 0;

    bool printed_p_only_coeff = false;
    bool printed_q_only_coeff = false;
    bool printed_nontrivial_root = false;

    for (
        std::size_t ci = 0;
        ci < cases.size();
        ++ci
    ) {
        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        u64 case_polynomials = 0;
        u64 case_identical_p = 0;
        u64 case_identical_q = 0;
        u64 case_identical_N = 0;

        u64 case_p_only_coeff = 0;
        u64 case_q_only_coeff = 0;
        u64 case_neither_coeff = 0;

        u64 case_root_p_only = 0;
        u64 case_root_q_only = 0;

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

                    u64 cg =
                        coefficient_gcd(P);

                    max_coefficient_gcd =
                        std::max(
                            max_coefficient_gcd,
                            cg
                        );

                    bool div_p =
                        all_coeff_divisible(
                            P,
                            p
                        );

                    bool div_q =
                        all_coeff_divisible(
                            P,
                            q
                        );

                    bool div_N =
                        all_coeff_divisible(
                            P,
                            N
                        );

                    if (div_p) {
                        ++case_identical_p;
                        ++total_identical_p;
                    }

                    if (div_q) {
                        ++case_identical_q;
                        ++total_identical_q;
                    }

                    if (div_N) {
                        ++case_identical_N;
                        ++total_identical_N;
                    }

                    if (div_p && !div_q) {
                        ++case_p_only_coeff;
                        ++total_p_only_coeff;

                        if (!printed_p_only_coeff) {
                            printed_p_only_coeff = true;

                            std::cout
                                << "FIRST_P_ONLY_COEFFICIENT_DIVISIBILITY\n"
                                << "base=" << base
                                << " m=" << m
                                << " n=" << n
                                << "\n"
                                << "degree="
                                << polynomial_degree(P)
                                << "\n"
                                << "coefficient_gcd="
                                << cg
                                << "\n";

                            print_polynomial(P);
                        }
                    } else if (
                        !div_p && div_q
                    ) {
                        ++case_q_only_coeff;
                        ++total_q_only_coeff;

                        if (!printed_q_only_coeff) {
                            printed_q_only_coeff = true;

                            std::cout
                                << "FIRST_Q_ONLY_COEFFICIENT_DIVISIBILITY\n"
                                << "base=" << base
                                << " m=" << m
                                << " n=" << n
                                << "\n"
                                << "degree="
                                << polynomial_degree(P)
                                << "\n"
                                << "coefficient_gcd="
                                << cg
                                << "\n";

                            print_polynomial(P);
                        }
                    } else if (!div_p && !div_q) {
                        ++case_neither_coeff;
                        ++total_neither_coeff;
                    }

                    if (cg != 1 &&
                        cg != p &&
                        cg != q &&
                        cg != N) {

                        ++total_nontrivial_coeff_gcd;

                        if (!printed_nontrivial_root) {
                            printed_nontrivial_root = true;

                            std::cout
                                << "FIRST_NONTRIVIAL_COEFFICIENT_GCD\n"
                                << "case=" << ci
                                << "\n"
                                << "base=" << base
                                << " m=" << m
                                << " n=" << n
                                << "\n"
                                << "degree="
                                << polynomial_degree(P)
                                << "\n"
                                << "coefficient_gcd="
                                << cg
                                << "\n";

                            print_polynomial(P);
                        }
                    }

                    /*
                        Root search only when at least one
                        field is not identically zero.
                    */
                    if (!div_p || !div_q) {
                        RootInfo R =
                            scan_roots(
                                P,
                                p,
                                q
                            );

                        total_root_p += R.roots_p;
                        total_root_q += R.roots_q;
                        total_root_common +=
                            R.common_roots;

                        total_root_p_only +=
                            R.p_only;

                        total_root_q_only +=
                            R.q_only;

                        if (R.p_only > 0) {
                            case_root_p_only +=
                                R.p_only;
                        }

                        if (R.q_only > 0) {
                            case_root_q_only +=
                                R.q_only;
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
            << "identical_mod_p="
            << case_identical_p
            << "\n";

        std::cout
            << "identical_mod_q="
            << case_identical_q
            << "\n";

        std::cout
            << "identical_mod_N="
            << case_identical_N
            << "\n";

        std::cout
            << "p_only_coefficient_divisibility="
            << case_p_only_coeff
            << "\n";

        std::cout
            << "q_only_coefficient_divisibility="
            << case_q_only_coeff
            << "\n";

        std::cout
            << "neither_coefficient_divisibility="
            << case_neither_coeff
            << "\n";

        std::cout
            << "root_p_only="
            << case_root_p_only
            << "\n";

        std::cout
            << "root_q_only="
            << case_root_q_only
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
        << "identical_mod_p="
        << total_identical_p
        << "\n";

    std::cout
        << "identical_mod_q="
        << total_identical_q
        << "\n";

    std::cout
        << "identical_mod_N="
        << total_identical_N
        << "\n";

    std::cout
        << "p_only_coefficient_divisibility="
        << total_p_only_coeff
        << "\n";

    std::cout
        << "q_only_coefficient_divisibility="
        << total_q_only_coeff
        << "\n";

    std::cout
        << "neither_coefficient_divisibility="
        << total_neither_coeff
        << "\n";

    std::cout
        << "nontrivial_coefficient_gcd="
        << total_nontrivial_coeff_gcd
        << "\n";

    std::cout
        << "root_p="
        << total_root_p
        << "\n";

    std::cout
        << "root_q="
        << total_root_q
        << "\n";

    std::cout
        << "root_common="
        << total_root_common
        << "\n";

    std::cout
        << "root_p_only="
        << total_root_p_only
        << "\n";

    std::cout
        << "root_q_only="
        << total_root_q_only
        << "\n";

    std::cout
        << "max_coefficient_gcd="
        << max_coefficient_gcd
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 303\n";

    return 0;
}
