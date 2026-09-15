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

struct Stats {
    std::size_t total_hits = 0;
    std::size_t d1_hits = 0;
    std::size_t d2_hits = 0;

    std::size_t factor_p_hits = 0;
    std::size_t factor_q_hits = 0;

    std::size_t exact_p = 0;
    std::size_t exact_q = 0;

    std::size_t multiple_p = 0;
    std::size_t multiple_q = 0;
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

    return (hi <= n / hi) ? hi : lo;
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
    Correct MISS-prefix formula:

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
            (i < static_cast<int>(md.d.size()))
                ? md.d[i]
                : 0;

        const u64 yi =
            (i < static_cast<int>(yd.d.size()))
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
            (i < static_cast<int>(yd.d.size()))
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

i128 H(
    u128 m,
    u128 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    const u128 misses =
        miss_prefix(
            m,
            n - 1,
            p
        );

    return static_cast<i128>(
        n - misses
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

i128 delta1(
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

    return
        H(m + step, n, p)
        -
        H(m, n, p);
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

i128 delta2(
    u128 m,
    unsigned r,
    unsigned s,
    u128 n,
    u64 p
) {
    if (!legal_pair(m, r, s, p)) {
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
        H(
            m + pr + ps,
            n,
            p
        );

    return d - b - c + a;
}

u128 abs_i128(i128 x) {
    if (x >= 0) {
        return static_cast<u128>(x);
    }

    return
        static_cast<u128>(-(x + 1))
        + 1;
}

u64 gcd_abs(
    i128 value,
    u128 N
) {
    const u128 magnitude =
        abs_i128(value);

    const u64 reduced =
        static_cast<u64>(
            magnitude % N
        );

    return std::gcd(
        reduced,
        static_cast<u64>(N)
    );
}

void analyze_factor_hit(
    const char* type,
    i128 value,
    u64 p,
    u64 q,
    u128 N,
    u64 base,
    u128 m,
    u128 n,
    unsigned r,
    unsigned s,
    Stats& stats
) {
    const u64 g =
        gcd_abs(
            value,
            N
        );

    if (g != p && g != q) {
        return;
    }

    ++stats.total_hits;

    if (type[1] == '1') {
        ++stats.d1_hits;
    } else {
        ++stats.d2_hits;
    }

    const u64 factor =
        (g == p) ? p : q;

    if (factor == p) {
        ++stats.factor_p_hits;
    } else {
        ++stats.factor_q_hits;
    }

    const u128 magnitude =
        abs_i128(value);

    const bool exact =
        magnitude ==
        static_cast<u128>(factor);

    if (factor == p && exact) {
        ++stats.exact_p;
    }

    if (factor == q && exact) {
        ++stats.exact_q;
    }

    if (magnitude % factor == 0) {
        if (factor == p) {
            ++stats.multiple_p;
        } else {
            ++stats.multiple_q;
        }
    }

    std::cout
        << "\nFACTOR HIT "
        << type
        << '\n';

    std::cout
        << "p="
        << p
        << " q="
        << q
        << '\n';

    std::cout
        << "N="
        << to_string_u128(N)
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
        << "value="
        << to_string_i128(value)
        << '\n';

    std::cout
        << "abs="
        << to_string_u128(magnitude)
        << '\n';

    std::cout
        << "gcd="
        << g
        << '\n';

    std::cout
        << "abs_div_factor="
        << to_string_u128(
               magnitude /
               static_cast<u128>(factor)
           )
        << '\n';

    std::cout
        << "abs_mod_p="
        << static_cast<u64>(
               magnitude %
               static_cast<u128>(p)
           )
        << '\n';

    std::cout
        << "abs_mod_q="
        << static_cast<u64>(
               magnitude %
               static_cast<u128>(q)
           )
        << '\n';

    std::cout
        << "abs_equals_factor="
        << (exact ? 1 : 0)
        << '\n';
}

void scan_semiprime(
    u64 p,
    u64 q,
    Stats& stats
) {
    const u128 N =
        static_cast<u128>(p) *
        static_cast<u128>(q);

    const u128 sqrtN =
        isqrt_u128(N);

    const u64 bases[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

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

                    const i128 value =
                        delta1(
                            m,
                            r,
                            n,
                            base
                        );

                    analyze_factor_hit(
                        "D1",
                        value,
                        p,
                        q,
                        N,
                        base,
                        m,
                        n,
                        r,
                        999,
                        stats
                    );
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

                        const i128 value =
                            delta2(
                                m,
                                r,
                                s,
                                n,
                                base
                            );

                        analyze_factor_hit(
                            "D2",
                            value,
                            p,
                            q,
                            N,
                            base,
                            m,
                            n,
                            r,
                            s,
                            stats
                        );
                    }
                }
            }
        }
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 295\n";

    std::cout
        << "ANALYSIS OF GENUINE SIGNED FACTOR HITS\n";

    std::cout
        << "D1/D2 MAGNITUDE AND FACTOR MULTIPLICITY\n\n";

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

    for (const auto& pair : cases) {
        scan_semiprime(
            pair.first,
            pair.second,
            stats
        );
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "factor_hits="
        << stats.total_hits
        << '\n';

    std::cout
        << "D1_hits="
        << stats.d1_hits
        << '\n';

    std::cout
        << "D2_hits="
        << stats.d2_hits
        << '\n';

    std::cout
        << "factor_p_hits="
        << stats.factor_p_hits
        << '\n';

    std::cout
        << "factor_q_hits="
        << stats.factor_q_hits
        << '\n';

    std::cout
        << "abs_equals_p="
        << stats.exact_p
        << '\n';

    std::cout
        << "abs_equals_q="
        << stats.exact_q
        << '\n';

    std::cout
        << "multiples_of_p="
        << stats.multiple_p
        << '\n';

    std::cout
        << "multiples_of_q="
        << stats.multiple_q
        << '\n';

    std::cout
        << "\nFINISHED EXPERIMENT 295\n";

    return 0;
}