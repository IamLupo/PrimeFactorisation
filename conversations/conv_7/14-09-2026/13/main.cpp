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

struct RootStats {
    u64 roots_p = 0;
    u64 roots_q = 0;
    u64 common_roots = 0;
    u64 p_only = 0;
    u64 q_only = 0;

    bool has_p_root = false;
    bool has_q_root = false;
    bool has_p_only = false;
    bool has_q_only = false;
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

u64 integer_sqrt(u64 n) {
    if (n == 0) {
        return 0;
    }

    u64 lo = 0;
    u64 hi = std::min<u64>(
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
        digits.push_back(x % base);
        x /= base;
    }

    return digits;
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

    u128 answer = 0;
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
    u64 old_count =
        miss_prefix(
            old_m,
            y,
            base
        );

    u64 new_count =
        miss_prefix(
            new_m,
            y,
            base
        );

    if (new_count < old_count) {
        std::cerr
            << "ERROR: negative local derivative\n"
            << "old_m=" << old_m
            << " new_m=" << new_m
            << " y=" << to_string_i128(y)
            << " base=" << base
            << "\n";

        std::exit(1);
    }

    return new_count - old_count;
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

/*
    Construct the exact path polynomial

        P(t) = sum_r Delta_r t^r

    where Delta_r is the total derivative contribution
    from digit position r along the canonical low-to-high path.
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

    for (std::size_t r = 0;
         r < digits.size();
         ++r) {

        u64 step =
            base_power(
                base,
                r
            );

        for (u64 count = 0;
             count < digits[r];
             ++count) {

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
            << "ERROR: canonical path failed\n"
            << "target_m=" << target_m
            << " current=" << current
            << "\n";

        std::exit(1);
    }

    return P;
}

/*
    Horner evaluation modulo mod.
*/
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

/*
    Exhaustive finite-field root search.

    For p:
        x = 0,...,p-1

    For q:
        x = 0,...,q-1
*/
RootStats scan_roots(
    const Polynomial &P,
    u64 p,
    u64 q,
    std::vector<u64> &p_roots,
    std::vector<u64> &q_roots
) {
    RootStats result;

    for (u64 x = 0;
         x < p;
         ++x) {

        if (
            evaluate_mod(
                P,
                x,
                p
            ) == 0
        ) {
            p_roots.push_back(x);
        }
    }

    for (u64 x = 0;
         x < q;
         ++x) {

        if (
            evaluate_mod(
                P,
                x,
                q
            ) == 0
        ) {
            q_roots.push_back(x);
        }
    }

    result.roots_p =
        static_cast<u64>(
            p_roots.size()
        );

    result.roots_q =
        static_cast<u64>(
            q_roots.size()
        );

    result.has_p_root =
        (result.roots_p > 0);

    result.has_q_root =
        (result.roots_q > 0);

    /*
        Since p and q are both < 2^32 here, integer x values
        are directly comparable for identifying common roots.
    */
    std::size_t i = 0;
    std::size_t j = 0;

    while (
        i < p_roots.size() &&
        j < q_roots.size()
    ) {
        if (p_roots[i] == q_roots[j]) {
            ++result.common_roots;
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

    result.p_only =
        result.roots_p -
        result.common_roots;

    result.q_only =
        result.roots_q -
        result.common_roots;

    result.has_p_only =
        (result.p_only > 0);

    result.has_q_only =
        (result.q_only > 0);

    return result;
}

int main() {
    std::cout
        << "START EXPERIMENT 302\n"
        << "FULL FINITE-FIELD ROOT SEARCH\n"
        << "DO THE PATH POLYNOMIALS HAVE FACTOR-SPECIFIC ROOTS?\n"
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

    u64 total_p_roots = 0;
    u64 total_q_roots = 0;
    u64 total_common_roots = 0;

    u64 polynomials_with_p_root = 0;
    u64 polynomials_with_q_root = 0;

    u64 polynomials_with_p_only_root = 0;
    u64 polynomials_with_q_only_root = 0;

    bool printed_p_only = false;
    bool printed_q_only = false;

    u64 first_p_only_case = 0;
    u64 first_q_only_case = 0;

    for (
        std::size_t ci = 0;
        ci < cases.size();
        ++ci
    ) {
        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        u64 case_polynomials = 0;
        u64 case_p_roots = 0;
        u64 case_q_roots = 0;
        u64 case_common_roots = 0;
        u64 case_p_only = 0;
        u64 case_q_only = 0;

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

                    std::vector<u64> p_roots;
                    std::vector<u64> q_roots;

                    RootStats S =
                        scan_roots(
                            P,
                            p,
                            q,
                            p_roots,
                            q_roots
                        );

                    ++case_polynomials;
                    ++total_polynomials;

                    case_p_roots +=
                        S.roots_p;

                    case_q_roots +=
                        S.roots_q;

                    case_common_roots +=
                        S.common_roots;

                    case_p_only +=
                        S.p_only;

                    case_q_only +=
                        S.q_only;

                    total_p_roots +=
                        S.roots_p;

                    total_q_roots +=
                        S.roots_q;

                    total_common_roots +=
                        S.common_roots;

                    if (S.has_p_root) {
                        ++polynomials_with_p_root;
                    }

                    if (S.has_q_root) {
                        ++polynomials_with_q_root;
                    }

                    if (S.has_p_only) {
                        ++polynomials_with_p_only_root;

                        if (!printed_p_only) {
                            printed_p_only = true;
                            first_p_only_case = ci;

                            std::cout
                                << "FIRST_P_ONLY_ROOT\n"
                                << "base=" << base
                                << " m=" << m
                                << " n=" << n
                                << "\n"
                                << "p_roots="
                                << S.roots_p
                                << "\n"
                                << "q_roots="
                                << S.roots_q
                                << "\n"
                                << "common_roots="
                                << S.common_roots
                                << "\n"
                                << "p_only="
                                << S.p_only
                                << "\n"
                                << "first_p_root="
                                << p_roots.front()
                                << "\n";
                        }
                    }

                    if (S.has_q_only) {
                        ++polynomials_with_q_only_root;

                        if (!printed_q_only) {
                            printed_q_only = true;
                            first_q_only_case = ci;

                            std::cout
                                << "FIRST_Q_ONLY_ROOT\n"
                                << "base=" << base
                                << " m=" << m
                                << " n=" << n
                                << "\n"
                                << "p_roots="
                                << S.roots_p
                                << "\n"
                                << "q_roots="
                                << S.roots_q
                                << "\n"
                                << "common_roots="
                                << S.common_roots
                                << "\n"
                                << "q_only="
                                << S.q_only
                                << "\n"
                                << "first_q_root="
                                << q_roots.front()
                                << "\n";
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
            << "total_p_roots="
            << case_p_roots
            << "\n";

        std::cout
            << "total_q_roots="
            << case_q_roots
            << "\n";

        std::cout
            << "common_roots="
            << case_common_roots
            << "\n";

        std::cout
            << "p_only="
            << case_p_only
            << "\n";

        std::cout
            << "q_only="
            << case_q_only
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
        << "total_p_roots="
        << total_p_roots
        << "\n";

    std::cout
        << "total_q_roots="
        << total_q_roots
        << "\n";

    std::cout
        << "common_roots="
        << total_common_roots
        << "\n";

    std::cout
        << "polynomials_with_p_root="
        << polynomials_with_p_root
        << "\n";

    std::cout
        << "polynomials_with_q_root="
        << polynomials_with_q_root
        << "\n";

    std::cout
        << "polynomials_with_p_only_root="
        << polynomials_with_p_only_root
        << "\n";

    std::cout
        << "polynomials_with_q_only_root="
        << polynomials_with_q_only_root
        << "\n";

    std::cout
        << "p_only_case_seen="
        << (printed_p_only ? 1 : 0)
        << "\n";

    std::cout
        << "q_only_case_seen="
        << (printed_q_only ? 1 : 0)
        << "\n";

    if (printed_p_only) {
        std::cout
            << "first_p_only_case="
            << first_p_only_case
            << "\n";
    }

    if (printed_q_only) {
        std::cout
            << "first_q_only_case="
            << first_q_only_case
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT 302\n";

    return 0;
}