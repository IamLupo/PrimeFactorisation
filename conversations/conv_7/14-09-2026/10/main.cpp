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

struct Aggregate {
    i128 S0 = 0;
    i128 S1 = 0;
    i128 S2 = 0;
    u64 count_steps = 0;
};

struct GcdStats {
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
    u128 v = neg ? static_cast<u128>(-x)
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
    u128 v = x < 0
        ? static_cast<u128>(-x)
        : static_cast<u128>(x);

    return static_cast<u64>(v);
}

u64 integer_sqrt(u64 n) {
    u64 lo = 0;
    u64 hi = 1;

    while (hi <= n / hi) {
        if (hi > n / 2) {
            hi = n;
            break;
        }

        hi *= 2;
    }

    if (hi > n / hi) {
        // hi is an upper bound with hi^2 > n.
    }

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

    MSB-first digit DP.

    For a digit position:
      x_i < y_i:
          contribution is determined by the allowed m digits.

      x_i = y_i:
          continue only if y_i <= m_i.
*/
u64 miss_prefix(u64 m, i128 y_signed, u64 base) {
    if (y_signed < 0) {
        return 0;
    }

    u64 y = static_cast<u64>(y_signed);

    std::vector<u64> md = base_digits(m, base);
    std::vector<u64> yd = base_digits(y, base);

    std::size_t L = std::max(md.size(), yd.size());

    md.resize(L, 0);
    yd.resize(L, 0);

    std::vector<u128> suffix(L + 1, 1);

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

        /*
            Digits x_i < y_i.

            Since x_i must also satisfy x_i <= m_i,
            there are min(y_i, m_i+1) choices.
        */
        u64 upper = std::min(yi, mi + 1);

        for (u64 xi = 0; xi < upper; ++xi) {
            ans += suffix[pos + 1];
        }

        /*
            Continue with x_i = y_i iff y_i <= m_i.
        */
        if (yi <= mi) {
            tight = true;
        } else {
            tight = false;
            break;
        }
    }

    /*
        Every digit matched y, so y itself is admissible.
    */
    if (tight) {
        ans += 1;
    }

    return static_cast<u64>(ans);
}

u64 local_delta(
    u64 old_m,
    u64 new_m,
    i128 y,
    u64 base
) {
    u64 a = miss_prefix(old_m, y, base);
    u64 b = miss_prefix(new_m, y, base);

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

u64 base_power(u64 base, std::size_t exponent) {
    u64 result = 1;

    for (std::size_t i = 0; i < exponent; ++i) {
        result *= base;
    }

    return result;
}

/*
    Canonical path:

        0 -> target_m

    with digits changed from low to high.

    For every individual no-carry increment:

        Delta = M_new(y) - M_old(y)

    and collect:

        S0 = sum Delta

        S1 = sum (r+1) Delta

        S2 = sum (r+1)^2 Delta
*/
Aggregate canonical_path_aggregate(
    u64 target_m,
    i128 y,
    u64 base
) {
    std::vector<u64> digits =
        base_digits(target_m, base);

    Aggregate out;

    u64 current = 0;

    for (std::size_t r = 0; r < digits.size(); ++r) {
        u64 step = base_power(base, r);

        for (u64 t = 0; t < digits[r]; ++t) {
            u64 next = current + step;

            u64 delta =
                local_delta(
                    current,
                    next,
                    y,
                    base
                );

            i128 w1 =
                static_cast<i128>(r + 1);

            i128 w2 =
                w1 * w1;

            i128 d =
                static_cast<i128>(delta);

            out.S0 += d;
            out.S1 += w1 * d;
            out.S2 += w2 * d;

            ++out.count_steps;

            current = next;
        }
    }

    if (current != target_m) {
        std::cerr
            << "ERROR: canonical path did not reach target\n"
            << "target=" << target_m
            << " current=" << current
            << "\n";

        std::exit(1);
    }

    return out;
}

u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

GcdStats classify(
    i128 value,
    u64 N,
    u64 p,
    u64 q
) {
    GcdStats out;

    u64 g =
        gcd_u64(
            abs_u64(value),
            N
        );

    if (g == 1) {
        ++out.gcd1;
        return out;
    }

    if (g == N) {
        ++out.gcdN;
        return out;
    }

    ++out.nontrivial;

    if (g == p && g == q) {
        ++out.both;
    } else if (g == p) {
        ++out.p_only;
    } else if (g == q) {
        ++out.q_only;
    }

    return out;
}

void accumulate(
    GcdStats &a,
    const GcdStats &b
) {
    a.gcd1 += b.gcd1;
    a.gcdN += b.gcdN;
    a.nontrivial += b.nontrivial;
    a.p_only += b.p_only;
    a.q_only += b.q_only;
    a.both += b.both;
}

void print_stats(
    const std::string &name,
    const GcdStats &s
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
        << "START EXPERIMENT 299\n"
        << "PATH-INTEGRATED DIGIT DERIVATIVES\n"
        << "CAN AN AGGREGATE OF SMALL D1/D2 VALUES EXPOSE A HIDDEN FACTOR?\n"
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
        -4, -3, -2, -1, 0,
         1,  2,  3,  4
    };

    u64 total_expressions = 0;
    u64 total_steps = 0;

    u64 path_identity_fail = 0;

    GcdStats total_S0;
    GcdStats total_S1;
    GcdStats total_S2;
    GcdStats total_C10;
    GcdStats total_C21;
    GcdStats total_C20;

    u64 max_abs_S0 = 0;
    u64 max_abs_S1 = 0;
    u64 max_abs_S2 = 0;

    for (std::size_t ci = 0; ci < cases.size(); ++ci) {
        const u64 p = cases[ci].p;
        const u64 q = cases[ci].q;
        const u64 N = p * q;

        std::cout
            << "CASE " << ci
            << " p=" << p
            << " q=" << q
            << "\n";

        u64 case_expressions = 0;
        u64 case_steps = 0;

        GcdStats case_S0;
        GcdStats case_S1;
        GcdStats case_S2;

        u64 root = integer_sqrt(N);

        for (u64 base : bases) {
            for (int off : offsets) {
                i128 m128 =
                    static_cast<i128>(root) +
                    static_cast<i128>(off);

                if (m128 <= 0) {
                    continue;
                }

                u64 m = static_cast<u64>(m128);

                /*
                    Keep exactly the same n family as the
                    intended experiment.
                */
                std::vector<u64> ns;

                ns.push_back(N);
                ns.push_back(m);
                ns.push_back(m + 1);

                if (m != 0 && m <= N / m) {
                    u64 D = N - m * m;

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

                    Aggregate A =
                        canonical_path_aggregate(
                            m,
                            y,
                            base
                        );

                    ++case_expressions;
                    ++total_expressions;

                    case_steps += A.count_steps;
                    total_steps += A.count_steps;

                    max_abs_S0 =
                        std::max(
                            max_abs_S0,
                            abs_u64(A.S0)
                        );

                    max_abs_S1 =
                        std::max(
                            max_abs_S1,
                            abs_u64(A.S1)
                        );

                    max_abs_S2 =
                        std::max(
                            max_abs_S2,
                            abs_u64(A.S2)
                        );

                    GcdStats a0 =
                        classify(
                            A.S0,
                            N,
                            p,
                            q
                        );

                    GcdStats a1 =
                        classify(
                            A.S1,
                            N,
                            p,
                            q
                        );

                    GcdStats a2 =
                        classify(
                            A.S2,
                            N,
                            p,
                            q
                        );

                    accumulate(case_S0, a0);
                    accumulate(case_S1, a1);
                    accumulate(case_S2, a2);

                    accumulate(total_S0, a0);
                    accumulate(total_S1, a1);
                    accumulate(total_S2, a2);

                    /*
                        Since the path starts at m=0,

                            M_0(y) = 1

                        for y >= 0.

                        Therefore:

                            S0 = M_m(y) - 1.
                    */
                    u64 final_M =
                        miss_prefix(
                            m,
                            y,
                            base
                        );

                    i128 expected =
                        static_cast<i128>(final_M) - 1;

                    if (A.S0 != expected) {
                        ++path_identity_fail;

                        if (path_identity_fail <= 10) {
                            std::cout
                                << "PATH_IDENTITY_FAIL"
                                << " base=" << base
                                << " m=" << m
                                << " n=" << n
                                << " S0="
                                << to_string_i128(A.S0)
                                << " expected="
                                << to_string_i128(expected)
                                << "\n";
                        }
                    }

                    /*
                        Remove trivial low-order weighting.
                    */
                    i128 C10 =
                        A.S1 - A.S0;

                    i128 C21 =
                        A.S2 - A.S1;

                    i128 C20 =
                        A.S2 -
                        2 * A.S1 +
                        A.S0;

                    accumulate(
                        total_C10,
                        classify(C10, N, p, q)
                    );

                    accumulate(
                        total_C21,
                        classify(C21, N, p, q)
                    );

                    accumulate(
                        total_C20,
                        classify(C20, N, p, q)
                    );
                }
            }
        }

        std::cout
            << "expressions="
            << case_expressions
            << "\n";

        std::cout
            << "path_steps="
            << case_steps
            << "\n";

        std::cout
            << "path_identity_failures="
            << path_identity_fail
            << "\n";

        print_stats("S0", case_S0);
        print_stats("S1", case_S1);
        print_stats("S2", case_S2);

        std::cout << "\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "expressions="
        << total_expressions
        << "\n";

    std::cout
        << "path_steps="
        << total_steps
        << "\n";

    std::cout
        << "path_identity_failures="
        << path_identity_fail
        << "\n";

    print_stats("S0", total_S0);
    print_stats("S1", total_S1);
    print_stats("S2", total_S2);
    print_stats("C10=S1-S0", total_C10);
    print_stats("C21=S2-S1", total_C21);
    print_stats(
        "C20=S2-2S1+S0",
        total_C20
    );

    std::cout
        << "max_abs_S0="
        << max_abs_S0
        << "\n";

    std::cout
        << "max_abs_S1="
        << max_abs_S1
        << "\n";

    std::cout
        << "max_abs_S2="
        << max_abs_S2
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 299\n";

    return 0;
}