#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Digits {
    std::vector<u64> d;
};

struct Stats {
    std::size_t expressions = 0;
    std::size_t nontrivial = 0;
    std::size_t gcd_one = 0;
    std::size_t gcd_n = 0;
};

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }

    std::string s;

    while (x > 0) {
        const unsigned digit =
            static_cast<unsigned>(x % 10);

        s.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    return s;
}

u128 isqrt_u128(u128 n) {
    if (n == 0) {
        return 0;
    }

    u128 lo = 0;
    u128 hi = n;

    while (lo + 1 < hi) {
        const u128 mid =
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

u128 power_u128(
    u128 base,
    unsigned exp
) {
    u128 result = 1;

    while (exp > 0) {
        if (exp & 1u) {
            result *= base;
        }

        base *= base;
        exp >>= 1u;
    }

    return result;
}

Digits to_base_p(
    u128 x,
    u64 p
) {
    Digits out;

    if (x == 0) {
        out.d.push_back(0);
        return out;
    }

    while (x > 0) {
        out.d.push_back(
            static_cast<u64>(
                x % static_cast<u128>(p)
            )
        );

        x /= static_cast<u128>(p);
    }

    return out;
}

u64 digit_at(
    const Digits& d,
    unsigned r
) {
    if (r >= d.d.size()) {
        return 0;
    }

    return d.d[r];
}

u128 total_miss(
    const Digits& md
) {
    u128 result = 1;

    for (u64 d : md.d) {
        result *=
            static_cast<u128>(d + 1);
    }

    return result;
}

/*
    Correct closed MISS-prefix formula.

        M_m(y)
        =
        #{x <= y : x <=_p m}
*/
u128 miss_prefix(
    u128 m,
    u128 y,
    u64 p
) {
    const Digits md =
        to_base_p(m, p);

    if (y >= m) {
        return total_miss(md);
    }

    const Digits yd =
        to_base_p(y, p);

    const unsigned len =
        static_cast<unsigned>(
            std::max(
                md.d.size(),
                yd.d.size()
            )
        );

    std::vector<u128> weight(
        len + 1,
        1
    );

    for (unsigned i = 0;
         i < len;
         ++i) {

        const u64 mi =
            (i < md.d.size())
                ? md.d[i]
                : 0;

        weight[i + 1] =
            weight[i] *
            static_cast<u128>(
                mi + 1
            );
    }

    int h = -1;

    for (int i =
             static_cast<int>(len) - 1;
         i >= 0;
         --i) {

        const u64 mi =
            (i < md.d.size())
                ? md.d[i]
                : 0;

        const u64 yi =
            (i < yd.d.size())
                ? yd.d[i]
                : 0;

        if (mi != yi) {
            h = i;
            break;
        }
    }

    if (h < 0) {
        return 1;
    }

    u128 result = 0;

    for (int i =
             static_cast<int>(len) - 1;
         i > h;
         --i) {

        const u64 yi =
            (i < yd.d.size())
                ? yd.d[i]
                : 0;

        result +=
            static_cast<u128>(yi) *
            weight[i];
    }

    const u64 yh =
        (h < static_cast<int>(yd.d.size()))
            ? yd.d[h]
            : 0;

    result +=
        static_cast<u128>(yh) *
        weight[h];

    const u128 ph =
        power_u128(
            static_cast<u128>(p),
            static_cast<unsigned>(h)
        );

    const u128 low_m =
        m % ph;

    const u128 low_y =
        y % ph;

    result +=
        miss_prefix(
            low_m,
            low_y,
            p
        );

    return result;
}

u128 H(
    u128 m,
    u128 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    return n -
        miss_prefix(
            m,
            n - 1,
            p
        );
}

bool legal_update(
    u128 m,
    unsigned r,
    u64 p
) {
    const u128 pr =
        power_u128(
            static_cast<u128>(p),
            r
        );

    const u64 mr =
        static_cast<u64>(
            (m / pr) %
            static_cast<u128>(p)
        );

    return mr + 1 < p;
}

u128 delta1(
    u128 m,
    unsigned r,
    u128 n,
    u64 p
) {
    if (!legal_update(m, r, p)) {
        return 0;
    }

    const u128 step =
        power_u128(
            static_cast<u128>(p),
            r
        );

    return H(
               m + step,
               n,
               p
           )
           -
           H(
               m,
               n,
               p
           );
}

bool legal_pair(
    u128 m,
    unsigned r,
    unsigned s,
    u64 p
) {
    if (r == s) {
        return false;
    }

    if (!legal_update(m, r, p)) {
        return false;
    }

    const u128 mr_step =
        power_u128(
            static_cast<u128>(p),
            r
        );

    const u128 m1 =
        m + mr_step;

    return legal_update(
        m1,
        s,
        p
    );
}

u128 delta2(
    u128 m,
    unsigned r,
    unsigned s,
    u128 n,
    u64 p
) {
    if (!legal_pair(
            m,
            r,
            s,
            p
        )) {
        return 0;
    }

    const u128 pr =
        power_u128(
            static_cast<u128>(p),
            r
        );

    const u128 ps =
        power_u128(
            static_cast<u128>(p),
            s
        );

    const u128 a = H(m, n, p);

    const u128 b =
        H(m + pr, n, p);

    const u128 c =
        H(m + ps, n, p);

    const u128 d =
        H(m + pr + ps, n, p);

    /*
        Mixed second difference.
    */
    return d - b - c + a;
}

void test_expression(
    const char* label,
    u128 value,
    u128 N,
    Stats& stats,
    bool print
) {
    ++stats.expressions;

    const u64 g =
        std::gcd(
            static_cast<u64>(value % N),
            static_cast<u64>(N)
        );

    if (g == 1) {
        ++stats.gcd_one;
        return;
    }

    if (g == N) {
        ++stats.gcd_n;
        return;
    }

    ++stats.nontrivial;

    if (print) {
        std::cout
            << "FACTOR HIT\n";

        std::cout
            << "expression="
            << label
            << '\n';

        std::cout
            << "gcd="
            << g
            << '\n';

        std::cout
            << "N="
            << to_string_u128(N)
            << '\n';
    }
}

void test_semiprime(
    u64 p_factor,
    u64 q_factor,
    Stats& stats,
    std::size_t index
) {
    const u128 N =
        static_cast<u128>(p_factor) *
        static_cast<u128>(q_factor);

    const u128 s =
        isqrt_u128(N);

    /*
        The algorithm only sees N.
        p_factor and q_factor are used only to
        classify a discovered gcd afterward.
    */

    const u64 bases[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    bool printed_factor = false;

    /*
        Test a small neighborhood around sqrt(N).
    */
    for (long long offset = -4;
         offset <= 4;
         ++offset) {

        u128 m = s;

        if (offset < 0) {
            const u128 d =
                static_cast<u128>(
                    -offset
                );

            if (m < d) {
                continue;
            }

            m -= d;
        } else {
            m +=
                static_cast<u128>(
                    offset
                );
        }

        /*
            Multiple natural n choices.
        */
        const u128 ns[] = {
            N,
            s,
            s + 1,
            N - s * s
        };

        for (u128 n : ns) {
            if (n == 0) {
                continue;
            }

            for (u64 base : bases) {
                /*
                    First differences.
                */
                for (unsigned r = 0;
                     r < 10;
                     ++r) {

                    if (!legal_update(
                            m,
                            r,
                            base
                        )) {
                        continue;
                    }

                    const u128 value =
                        delta1(
                            m,
                            r,
                            n,
                            base
                        );

                    const std::string label =
                        "D1";

                    const bool before =
                        stats.nontrivial > 0;

                    test_expression(
                        label.c_str(),
                        value,
                        N,
                        stats,
                        !printed_factor
                    );

                    if (!before &&
                        stats.nontrivial > 0) {

                        printed_factor = true;

                        std::cout
                            << "case="
                            << index
                            << '\n';

                        std::cout
                            << "base="
                            << base
                            << '\n';

                        std::cout
                            << "m="
                            << to_string_u128(m)
                            << '\n';

                        std::cout
                            << "r="
                            << r
                            << '\n';

                        std::cout
                            << "n="
                            << to_string_u128(n)
                            << '\n';

                        std::cout
                            << "known_factor="
                            << p_factor
                            << '\n';

                        std::cout
                            << "other_factor="
                            << q_factor
                            << '\n';
                    }
                }

                /*
                    Second differences.
                */
                for (unsigned r = 0;
                     r < 7;
                     ++r) {

                    for (unsigned ss = r + 1;
                         ss < 8;
                         ++ss) {

                        if (!legal_pair(
                                m,
                                r,
                                ss,
                                base
                            )) {
                            continue;
                        }

                        const u128 value =
                            delta2(
                                m,
                                r,
                                ss,
                                n,
                                base
                            );

                        test_expression(
                            "D2",
                            value,
                            N,
                            stats,
                            !printed_factor
                        );

                        if (!printed_factor &&
                            stats.nontrivial > 0) {

                            printed_factor = true;

                            std::cout
                                << "case="
                                << index
                                << '\n';

                            std::cout
                                << "base="
                                << base
                                << '\n';

                            std::cout
                                << "m="
                                << to_string_u128(m)
                                << '\n';

                            std::cout
                                << "r="
                                << r
                                << '\n';

                            std::cout
                                << "s="
                                << ss
                                << '\n';

                            std::cout
                                << "n="
                                << to_string_u128(n)
                                << '\n';

                            std::cout
                                << "known_factor="
                                << p_factor
                                << '\n';

                            std::cout
                                << "other_factor="
                                << q_factor
                                << '\n';
                        }
                    }
                }
            }
        }
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 290\n";

    std::cout
        << "FACTOR-SENSITIVE DIGIT-CALCULUS SCAN\n";

    std::cout
        << "FIRST AND SECOND LOCAL DIGIT DIFFERENCES\n";

    std::cout
        << "GCD WITH HIDDEN SEMIPRIMES\n\n";

    /*
        Moderate semiprimes.

        The actual algorithm only receives N.
        The factors are retained solely for validation.
    */
    const std::vector<std::pair<u64, u64>> cases = {
        {1009, 1013},
        {2003, 2011},
        {3001, 3011},
        {4001, 4007},
        {5003, 5009},
        {7001, 7013},
        {10007, 10009},
        {12011, 12037},
        {15013, 15017},
        {20011, 20021},
        {30011, 30013},
        {40009, 40013},
        {50021, 50023},
        {70001, 70009},
        {100003, 100019},
        {120011, 120017},
        {150001, 150013}
    };

    Stats stats;

    for (std::size_t i = 0;
         i < cases.size();
         ++i) {

        test_semiprime(
            cases[i].first,
            cases[i].second,
            stats,
            i
        );
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "semiprimes="
        << cases.size()
        << '\n';

    std::cout
        << "expressions="
        << stats.expressions
        << '\n';

    std::cout
        << "gcd_one="
        << stats.gcd_one
        << '\n';

    std::cout
        << "gcd_N="
        << stats.gcd_n
        << '\n';

    std::cout
        << "nontrivial_factor_hits="
        << stats.nontrivial
        << '\n';

    const bool found =
        stats.nontrivial > 0;

    std::cout
        << "\nFACTORING_SIGNAL="
        << (found ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 290\n";

    return 0;
}
