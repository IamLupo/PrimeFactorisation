#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

struct Digits {
    std::vector<u64> d;
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

std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    const bool negative = x < 0;

    u128 magnitude;

    if (negative) {
        magnitude =
            static_cast<u128>(-(x + 1)) + 1;
    } else {
        magnitude =
            static_cast<u128>(x);
    }

    std::string s =
        to_string_u128(magnitude);

    if (negative) {
        s.insert(s.begin(), '-');
    }

    return s;
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
    Correct MISS prefix formula validated by 288.
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
            i < md.d.size()
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
            i < static_cast<int>(md.d.size())
                ? md.d[i]
                : 0;

        const u64 yi =
            i < static_cast<int>(yd.d.size())
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
            i < static_cast<int>(yd.d.size())
                ? yd.d[i]
                : 0;

        result +=
            static_cast<u128>(yi) *
            weight[i];
    }

    const u64 yh =
        h < static_cast<int>(yd.d.size())
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

i128 H(
    u128 m,
    u128 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    return static_cast<i128>(
        n -
        miss_prefix(
            m,
            n - 1,
            p
        )
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

i128 delta1_signed(
    u128 m,
    unsigned r,
    u128 n,
    u64 p
) {
    if (!legal_update(
            m,
            r,
            p
        )) {
        return 0;
    }

    const u128 step =
        power_u128(
            static_cast<u128>(p),
            r
        );

    return
        H(
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

u128 delta1_unsigned_wrap(
    u128 m,
    unsigned r,
    u128 n,
    u64 p
) {
    if (!legal_update(
            m,
            r,
            p
        )) {
        return 0;
    }

    const u128 step =
        power_u128(
            static_cast<u128>(p),
            r
        );

    /*
        Reproduce the dangerous unsigned calculation
        from the earlier experiments.
    */
    const u128 a =
        static_cast<u128>(
            H(
                m,
                n,
                p
            )
        );

    const u128 b =
        static_cast<u128>(
            H(
                m + step,
                n,
                p
            )
        );

    return b - a;
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

    if (!legal_update(
            m,
            r,
            p
        )) {
        return false;
    }

    const u128 pr =
        power_u128(
            static_cast<u128>(p),
            r
        );

    return legal_update(
        m + pr,
        s,
        p
    );
}

i128 delta2_signed(
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

    const i128 a =
        H(m, n, p);

    const i128 b =
        H(m + pr, n, p);

    const i128 c =
        H(m + ps, n, p);

    const i128 d =
        H(m + pr + ps, n, p);

    return d - b - c + a;
}

u128 abs_i128(
    i128 x
) {
    if (x >= 0) {
        return static_cast<u128>(x);
    }

    return
        static_cast<u128>(
            -(x + 1)
        ) + 1;
}

u64 gcd_value(
    u128 value,
    u128 N
) {
    return std::gcd(
        static_cast<u64>(
            value % N
        ),
        static_cast<u64>(N)
    );
}

u64 gcd_signed(
    i128 value,
    u128 N
) {
    return gcd_value(
        abs_i128(value),
        N
    );
}

void print_diagnostic(
    const char* type,
    u64 p,
    u64 q,
    u64 base,
    u128 m,
    u128 n,
    unsigned r,
    unsigned s,
    i128 signed_value,
    u128 unsigned_value
) {
    const u128 N =
        static_cast<u128>(p) *
        static_cast<u128>(q);

    const u64 gcd_signed_value =
        gcd_signed(
            signed_value,
            N
        );

    const u64 gcd_unsigned_value =
        gcd_value(
            unsigned_value,
            N
        );

    const u64 signed_mod_p =
        static_cast<u64>(
            abs_i128(signed_value) %
            static_cast<u128>(p)
        );

    const u64 signed_mod_q =
        static_cast<u64>(
            abs_i128(signed_value) %
            static_cast<u128>(q)
        );

    const u64 unsigned_mod_p =
        static_cast<u64>(
            unsigned_value %
            static_cast<u128>(p)
        );

    const u64 unsigned_mod_q =
        static_cast<u64>(
            unsigned_value %
            static_cast<u128>(q)
        );

    std::cout
        << "\nREGRESSION HIT\n";

    std::cout
        << "type="
        << type
        << '\n';

    std::cout
        << "N="
        << to_string_u128(N)
        << '\n';

    std::cout
        << "p="
        << p
        << '\n';

    std::cout
        << "q="
        << q
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
        << "n="
        << to_string_u128(n)
        << '\n';

    std::cout
        << "r="
        << r
        << '\n';

    if (s != 999) {
        std::cout
            << "s="
            << s
            << '\n';
    }

    std::cout
        << "signed="
        << to_string_i128(signed_value)
        << '\n';

    std::cout
        << "unsigned_wrap="
        << to_string_u128(unsigned_value)
        << '\n';

    std::cout
        << "signed_mod_p="
        << signed_mod_p
        << '\n';

    std::cout
        << "signed_mod_q="
        << signed_mod_q
        << '\n';

    std::cout
        << "unsigned_mod_p="
        << unsigned_mod_p
        << '\n';

    std::cout
        << "unsigned_mod_q="
        << unsigned_mod_q
        << '\n';

    std::cout
        << "signed_gcd="
        << gcd_signed_value
        << '\n';

    std::cout
        << "unsigned_gcd="
        << gcd_unsigned_value
        << '\n';
}

void run_known_case(
    u64 p,
    u64 q,
    u64 base,
    u128 m,
    u128 n,
    unsigned r,
    unsigned s
) {
    const u128 N =
        static_cast<u128>(p) *
        static_cast<u128>(q);

    if (s == 999) {
        const i128 signed_value =
            delta1_signed(
                m,
                r,
                n,
                base
            );

        const u128 unsigned_value =
            delta1_unsigned_wrap(
                m,
                r,
                n,
                base
            );

        const u64 gs =
            gcd_signed(
                signed_value,
                N
            );

        const u64 gu =
            gcd_value(
                unsigned_value,
                N
            );

        /*
            Print every case where either path produces
            a nontrivial gcd, plus selected known case 0.
        */
        if ((gs != 1 && gs != N) ||
            (gu != 1 && gu != N)) {

            print_diagnostic(
                "D1",
                p,
                q,
                base,
                m,
                n,
                r,
                999,
                signed_value,
                unsigned_value
            );
        }

        return;
    }

    if (!legal_pair(
            m,
            r,
            s,
            base
        )) {
        return;
    }

    const i128 signed_value =
        delta2_signed(
            m,
            r,
            s,
            n,
            base
        );

    /*
        Reproduce old unsigned D2 behavior.
    */
    const u128 pr =
        power_u128(
            static_cast<u128>(base),
            r
        );

    const u128 ps =
        power_u128(
            static_cast<u128>(base),
            s
        );

    const u128 a =
        static_cast<u128>(
            H(
                m,
                n,
                base
            )
        );

    const u128 b =
        static_cast<u128>(
            H(
                m + pr,
                n,
                base
            )
        );

    const u128 c =
        static_cast<u128>(
            H(
                m + ps,
                n,
                base
            )
        );

    const u128 d =
        static_cast<u128>(
            H(
                m + pr + ps,
                n,
                base
            )
        );

    const u128 unsigned_value =
        d - b - c + a;

    const u64 gs =
        gcd_signed(
            signed_value,
            N
        );

    const u64 gu =
        gcd_value(
            unsigned_value,
            N
        );

    if ((gs != 1 && gs != N) ||
        (gu != 1 && gu != N)) {

        print_diagnostic(
            "D2",
            p,
            q,
            base,
            m,
            n,
            r,
            s,
            signed_value,
            unsigned_value
        );
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 294\n";

    std::cout
        << "290/292/293 REGRESSION FORENSICS\n";

    std::cout
        << "SIGNED VS UNSIGNED D1/D2\n";

    std::cout
        << "REPRODUCE NONTRIVIAL GCD CASES\n\n";

    /*
        First reproduce the explicitly printed
        Experiment 290 D1 hit.
    */
    std::cout
        << "KNOWN EXPERIMENT 290 CASE\n";

    run_known_case(
        1009,
        1013,
        7,
        1007,
        1010,
        1,
        999
    );

    /*
        Reproduce several explicitly printed D2 cases
        from Experiment 290.
    */
    std::cout
        << "\nKNOWN EXPERIMENT 290 D2 CASES\n";

    run_known_case(
        2003,
        2011,
        2,
        2002,
        4028033,
        0,
        2
    );

    run_known_case(
        3001,
        3011,
        2,
        3001,
        9036011,
        1,
        2
    );

    run_known_case(
        4001,
        4007,
        2,
        3999,
        16032007,
        5,
        6
    );

    run_known_case(
        5003,
        5009,
        2,
        5001,
        25060027,
        1,
        2
    );

    run_known_case(
        7001,
        7013,
        2,
        7002,
        49098013,
        0,
        2
    );

    run_known_case(
        10007,
        10009,
        2,
        10003,
        100160063,
        2,
        3
    );

    /*
        Full original expression generator, but only print
        disagreements between signed and unsigned gcds.
    */
    std::cout
        << "\nFULL REGRESSION SCAN\n";

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

    std::size_t expressions = 0;
    std::size_t signed_factor = 0;
    std::size_t unsigned_factor = 0;
    std::size_t disagreement = 0;

    const u64 bases[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    for (const auto& [p, q] : cases) {
        const u128 N =
            static_cast<u128>(p) *
            static_cast<u128>(q);

        const u128 sqrtN =
            isqrt_u128(N);

        const u128 n_values[] = {
            N,
            sqrtN,
            sqrtN + 1,
            N - sqrtN * sqrtN
        };

        for (long long offset = -4;
             offset <= 4;
             ++offset) {

            u128 m = sqrtN;

            if (offset < 0) {
                m -=
                    static_cast<u128>(
                        -offset
                    );
            } else {
                m +=
                    static_cast<u128>(
                        offset
                    );
            }

            for (u128 n : n_values) {
                if (n == 0) {
                    continue;
                }

                for (u64 base : bases) {

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

                        const i128 signed_value =
                            delta1_signed(
                                m,
                                r,
                                n,
                                base
                            );

                        const u128 unsigned_value =
                            delta1_unsigned_wrap(
                                m,
                                r,
                                n,
                                base
                            );

                        const u64 gs =
                            gcd_signed(
                                signed_value,
                                N
                            );

                        const u64 gu =
                            gcd_value(
                                unsigned_value,
                                N
                            );

                        ++expressions;

                        const bool sf =
                            gs != 1 &&
                            gs != N;

                        const bool uf =
                            gu != 1 &&
                            gu != N;

                        if (sf) {
                            ++signed_factor;
                        }

                        if (uf) {
                            ++unsigned_factor;
                        }

                        if (sf != uf) {
                            ++disagreement;

                            if (disagreement <= 10) {
                                print_diagnostic(
                                    "D1-DISAGREEMENT",
                                    p,
                                    q,
                                    base,
                                    m,
                                    n,
                                    r,
                                    999,
                                    signed_value,
                                    unsigned_value
                                );
                            }
                        }
                    }

                    for (unsigned r = 0;
                         r < 7;
                         ++r) {

                        for (unsigned s = r + 1;
                             s < 8;
                             ++s) {

                            if (!legal_pair(
                                    m,
                                    r,
                                    s,
                                    base
                                )) {
                                continue;
                            }

                            const i128 signed_value =
                                delta2_signed(
                                    m,
                                    r,
                                    s,
                                    n,
                                    base
                                );

                            const u128 pr =
                                power_u128(
                                    static_cast<u128>(base),
                                    r
                                );

                            const u128 ps =
                                power_u128(
                                    static_cast<u128>(base),
                                    s
                                );

                            const u128 a =
                                static_cast<u128>(
                                    H(m, n, base)
                                );

                            const u128 b =
                                static_cast<u128>(
                                    H(
                                        m + pr,
                                        n,
                                        base
                                    )
                                );

                            const u128 c =
                                static_cast<u128>(
                                    H(
                                        m + ps,
                                        n,
                                        base
                                    )
                                );

                            const u128 d =
                                static_cast<u128>(
                                    H(
                                        m + pr + ps,
                                        n,
                                        base
                                    )
                                );

                            const u128 unsigned_value =
                                d - b - c + a;

                            const u64 gs =
                                gcd_signed(
                                    signed_value,
                                    N
                                );

                            const u64 gu =
                                gcd_value(
                                    unsigned_value,
                                    N
                                );

                            ++expressions;

                            const bool sf =
                                gs != 1 &&
                                gs != N;

                            const bool uf =
                                gu != 1 &&
                                gu != N;

                            if (sf) {
                                ++signed_factor;
                            }

                            if (uf) {
                                ++unsigned_factor;
                            }

                            if (sf != uf) {
                                ++disagreement;

                                if (disagreement <= 10) {
                                    print_diagnostic(
                                        "D2-DISAGREEMENT",
                                        p,
                                        q,
                                        base,
                                        m,
                                        n,
                                        r,
                                        s,
                                        signed_value,
                                        unsigned_value
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "expressions="
        << expressions
        << '\n';

    std::cout
        << "signed_factor_hits="
        << signed_factor
        << '\n';

    std::cout
        << "unsigned_factor_hits="
        << unsigned_factor
        << '\n';

    std::cout
        << "signed_unsigned_disagreements="
        << disagreement
        << '\n';

    std::cout
        << "\nFINISHED EXPERIMENT 294\n";

    return 0;
}
