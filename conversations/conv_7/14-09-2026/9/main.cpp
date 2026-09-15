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
    std::size_t expressions = 0;

    std::size_t gcd1 = 0;
    std::size_t gcdN = 0;
    std::size_t nontrivial = 0;

    std::size_t abs_lt_p = 0;
    std::size_t abs_p_to_q = 0;
    std::size_t abs_ge_q = 0;

    u128 max_abs_d1 = 0;
    u128 max_abs_d2 = 0;

    u128 max_abs_all = 0;

    std::size_t d1_abs_lt_p = 0;
    std::size_t d1_abs_p_to_q = 0;
    std::size_t d1_abs_ge_q = 0;

    std::size_t d2_abs_lt_p = 0;
    std::size_t d2_abs_p_to_q = 0;
    std::size_t d2_abs_ge_q = 0;
};

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }

    std::string s;

    while (x > 0) {
        const unsigned digit =
            static_cast<unsigned>(
                x % 10
            );

        s.push_back(
            static_cast<char>(
                '0' + digit
            )
        );

        x /= 10;
    }

    std::reverse(
        s.begin(),
        s.end()
    );

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

    return hi <= n / hi ? hi : lo;
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
                x %
                static_cast<u128>(p)
            )
        );

        x /=
            static_cast<u128>(p);
    }

    return out;
}

u128 total_miss(
    const Digits& md
) {
    u128 result = 1;

    for (u64 d : md.d) {
        result *=
            static_cast<u128>(
                d + 1
            );
    }

    return result;
}

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
            i < static_cast<int>(
                md.d.size()
            )
                ? md.d[i]
                : 0;

        const u64 yi =
            i < static_cast<int>(
                yd.d.size()
            )
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
            i < static_cast<int>(
                yd.d.size()
            )
                ? yd.d[i]
                : 0;

        result +=
            static_cast<u128>(yi) *
            weight[i];
    }

    const u64 yh =
        h < static_cast<int>(
            yd.d.size()
        )
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
    u64 base
) {
    if (n == 0) {
        return 0;
    }

    return static_cast<i128>(
        n -
        miss_prefix(
            m,
            n - 1,
            base
        )
    );
}

bool legal_update(
    u128 m,
    unsigned r,
    u64 base
) {
    const u128 pr =
        power_u128(
            static_cast<u128>(base),
            r
        );

    const u64 mr =
        static_cast<u64>(
            (m / pr) %
            static_cast<u128>(base)
        );

    return mr + 1 < base;
}

i128 delta1(
    u128 m,
    unsigned r,
    u128 n,
    u64 base
) {
    if (!legal_update(
            m,
            r,
            base
        )) {
        return 0;
    }

    const u128 step =
        power_u128(
            static_cast<u128>(base),
            r
        );

    return
        H(
            m + step,
            n,
            base
        )
        -
        H(
            m,
            n,
            base
        );
}

bool legal_pair(
    u128 m,
    unsigned r,
    unsigned s,
    u64 base
) {
    if (r == s) {
        return false;
    }

    if (!legal_update(
            m,
            r,
            base
        )) {
        return false;
    }

    const u128 pr =
        power_u128(
            static_cast<u128>(base),
            r
        );

    return legal_update(
        m + pr,
        s,
        base
    );
}

i128 delta2(
    u128 m,
    unsigned r,
    unsigned s,
    u128 n,
    u64 base
) {
    if (!legal_pair(
            m,
            r,
            s,
            base
        )) {
        return 0;
    }

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

    const i128 a =
        H(
            m,
            n,
            base
        );

    const i128 b =
        H(
            m + pr,
            n,
            base
        );

    const i128 c =
        H(
            m + ps,
            n,
            base
        );

    const i128 d =
        H(
            m + pr + ps,
            n,
            base
        );

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

void classify_magnitude(
    u128 magnitude,
    u64 p,
    u64 q,
    std::size_t& less_p,
    std::size_t& p_to_q,
    std::size_t& ge_q,
    u128& maximum
) {
    if (magnitude < p) {
        ++less_p;
    } else if (magnitude < q) {
        ++p_to_q;
    } else {
        ++ge_q;
    }

    maximum =
        std::max(
            maximum,
            magnitude
        );
}

void classify_gcd(
    i128 value,
    u64 p,
    u64 q,
    u128 N,
    std::size_t& gcd1,
    std::size_t& gcdN,
    std::size_t& nontrivial
) {
    const u128 magnitude =
        abs_i128(value);

    const u64 g =
        std::gcd(
            static_cast<u64>(
                magnitude % N
            ),
            static_cast<u64>(N)
        );

    if (g == 1) {
        ++gcd1;
    } else if (g ==
               static_cast<u64>(N)) {
        ++gcdN;
    } else {
        ++nontrivial;
    }
}

void process_d1(
    u128 m,
    u128 n,
    unsigned r,
    u64 base,
    u64 p,
    u64 q,
    u128 N,
    Stats& stats
) {
    const i128 value =
        delta1(
            m,
            r,
            n,
            base
        );

    const u128 magnitude =
        abs_i128(value);

    ++stats.expressions;

    classify_magnitude(
        magnitude,
        p,
        q,
        stats.abs_lt_p,
        stats.abs_p_to_q,
        stats.abs_ge_q,
        stats.max_abs_all
    );

    classify_magnitude(
        magnitude,
        p,
        q,
        stats.d1_abs_lt_p,
        stats.d1_abs_p_to_q,
        stats.d1_abs_ge_q,
        stats.max_abs_d1
    );

    classify_gcd(
        value,
        p,
        q,
        N,
        stats.gcd1,
        stats.gcdN,
        stats.nontrivial
    );
}

void process_d2(
    u128 m,
    u128 n,
    unsigned r,
    unsigned s,
    u64 base,
    u64 p,
    u64 q,
    u128 N,
    Stats& stats
) {
    const i128 value =
        delta2(
            m,
            r,
            s,
            n,
            base
        );

    const u128 magnitude =
        abs_i128(value);

    ++stats.expressions;

    classify_magnitude(
        magnitude,
        p,
        q,
        stats.abs_lt_p,
        stats.abs_p_to_q,
        stats.abs_ge_q,
        stats.max_abs_all
    );

    classify_magnitude(
        magnitude,
        p,
        q,
        stats.d2_abs_lt_p,
        stats.d2_abs_p_to_q,
        stats.d2_abs_ge_q,
        stats.max_abs_d2
    );

    classify_gcd(
        value,
        p,
        q,
        N,
        stats.gcd1,
        stats.gcdN,
        stats.nontrivial
    );
}

void run_semiprime(
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
        2, 3, 5, 7,
        11, 13
    };

    const long long offsets[] = {
        -4, -3, -2, -1,
        0, 1, 2, 3, 4
    };

    for (long long offset : offsets) {
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

        const u128 n_values[] = {
            N,
            m,
            m + 1,
            N - (
                m * m <= N
                    ? m * m
                    : N
            )
        };

        for (u128 n : n_values) {
            if (n == 0) {
                continue;
            }

            for (u64 base : bases) {

                for (unsigned r = 0;
                     r < 8;
                     ++r) {

                    if (!legal_update(
                            m,
                            r,
                            base
                        )) {
                        continue;
                    }

                    process_d1(
                        m,
                        n,
                        r,
                        base,
                        p,
                        q,
                        N,
                        stats
                    );
                }

                for (unsigned r = 0;
                     r < 6;
                     ++r) {

                    for (unsigned s = r + 1;
                         s < 7;
                         ++s) {

                        if (!legal_pair(
                                m,
                                r,
                                s,
                                base
                            )) {
                            continue;
                        }

                        process_d2(
                            m,
                            n,
                            r,
                            s,
                            base,
                            p,
                            q,
                            N,
                            stats
                        );
                    }
                }
            }
        }
    }
}

void print_stats(
    const char* name,
    const Stats& s
) {
    std::cout
        << '\n'
        << name
        << '\n';

    std::cout
        << "expressions="
        << s.expressions
        << '\n';

    std::cout
        << "gcd1="
        << s.gcd1
        << '\n';

    std::cout
        << "gcdN="
        << s.gcdN
        << '\n';

    std::cout
        << "nontrivial="
        << s.nontrivial
        << '\n';

    std::cout
        << "abs_lt_p="
        << s.abs_lt_p
        << '\n';

    std::cout
        << "abs_p_to_q="
        << s.abs_p_to_q
        << '\n';

    std::cout
        << "abs_ge_q="
        << s.abs_ge_q
        << '\n';

    std::cout
        << "max_abs_D1="
        << to_string_u128(
               s.max_abs_d1
           )
        << '\n';

    std::cout
        << "max_abs_D2="
        << to_string_u128(
               s.max_abs_d2
           )
        << '\n';

    std::cout
        << "max_abs_all="
        << to_string_u128(
               s.max_abs_all
           )
        << '\n';

    std::cout
        << "D1_abs_lt_p="
        << s.d1_abs_lt_p
        << '\n';

    std::cout
        << "D1_abs_p_to_q="
        << s.d1_abs_p_to_q
        << '\n';

    std::cout
        << "D1_abs_ge_q="
        << s.d1_abs_ge_q
        << '\n';

    std::cout
        << "D2_abs_lt_p="
        << s.d2_abs_lt_p
        << '\n';

    std::cout
        << "D2_abs_p_to_q="
        << s.d2_abs_p_to_q
        << '\n';

    std::cout
        << "D2_abs_ge_q="
        << s.d2_abs_ge_q
        << '\n';
}

int main() {
    std::cout
        << "START EXPERIMENT 298\n";

    std::cout
        << "MAGNITUDE OF DIGIT DERIVATIVES\n";

    std::cout
        << "CAN D1/D2 EVEN REACH A HIDDEN PRIME FACTOR?\n\n";

    const std::vector<std::pair<u64, u64>> cases = {
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

    Stats total;

    for (std::size_t i = 0;
         i < cases.size();
         ++i) {

        Stats local;

        run_semiprime(
            cases[i].first,
            cases[i].second,
            local
        );

        std::cout
            << "\nCASE "
            << i
            << " p="
            << cases[i].first
            << " q="
            << cases[i].second
            << '\n';

        print_stats(
            "LOCAL",
            local
        );

        total.expressions +=
            local.expressions;

        total.gcd1 +=
            local.gcd1;

        total.gcdN +=
            local.gcdN;

        total.nontrivial +=
            local.nontrivial;

        total.abs_lt_p +=
            local.abs_lt_p;

        total.abs_p_to_q +=
            local.abs_p_to_q;

        total.abs_ge_q +=
            local.abs_ge_q;

        total.d1_abs_lt_p +=
            local.d1_abs_lt_p;

        total.d1_abs_p_to_q +=
            local.d1_abs_p_to_q;

        total.d1_abs_ge_q +=
            local.d1_abs_ge_q;

        total.d2_abs_lt_p +=
            local.d2_abs_lt_p;

        total.d2_abs_p_to_q +=
            local.d2_abs_p_to_q;

        total.d2_abs_ge_q +=
            local.d2_abs_ge_q;

        total.max_abs_d1 =
            std::max(
                total.max_abs_d1,
                local.max_abs_d1
            );

        total.max_abs_d2 =
            std::max(
                total.max_abs_d2,
                local.max_abs_d2
            );

        total.max_abs_all =
            std::max(
                total.max_abs_all,
                local.max_abs_all
            );
    }

    std::cout
        << "\n============================\n";

    print_stats(
        "TOTAL",
        total
    );

    std::cout
        << "\nFINISHED EXPERIMENT 298\n";

    return 0;
}
