#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

struct Stats {
    std::size_t expressions = 0;

    std::size_t structured_one = 0;
    std::size_t structured_n = 0;
    std::size_t structured_factor = 0;

    std::size_t random_one = 0;
    std::size_t random_n = 0;
    std::size_t random_factor = 0;
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

    /*
        Avoid signed overflow when negating the minimum
        representable i128.
    */
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
        s.insert(
            s.begin(),
            '-'
        );
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

struct Digits {
    std::vector<u64> d;
};

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
    Correct MISS prefix formula:

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

u128 H_unsigned(
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

i128 H_signed(
    u128 m,
    u128 n,
    u64 p
) {
    return static_cast<i128>(
        H_unsigned(
            m,
            n,
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

/*
    Signed first difference.
*/
i128 delta1(
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

    const i128 a =
        H_signed(
            m,
            n,
            p
        );

    const i128 b =
        H_signed(
            m + step,
            n,
            p
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

    const u128 m1 =
        m + pr;

    return legal_update(
        m1,
        s,
        p
    );
}

/*
    Signed mixed second difference.
*/
i128 delta2(
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
        H_signed(
            m,
            n,
            p
        );

    const i128 b =
        H_signed(
            m + pr,
            n,
            p
        );

    const i128 c =
        H_signed(
            m + ps,
            n,
            p
        );

    const i128 d =
        H_signed(
            m + pr + ps,
            n,
            p
        );

    return d - b - c + a;
}

u128 abs_i128(
    i128 x
) {
    if (x >= 0) {
        return static_cast<u128>(x);
    }

    return static_cast<u128>(
        -(x + 1)
    ) + 1;
}

u64 gcd_value_with_N(
    i128 value,
    u128 N
) {
    const u128 magnitude =
        abs_i128(value);

    const u64 reduced =
        static_cast<u64>(
            magnitude %
            N
        );

    return std::gcd(
        reduced,
        static_cast<u64>(N)
    );
}

void record_signed_value(
    i128 value,
    u128 N,
    bool structured,
    Stats& stats
) {
    const u64 g =
        gcd_value_with_N(
            value,
            N
        );

    if (structured) {
        if (g == 1) {
            ++stats.structured_one;
        } else if (g == N) {
            ++stats.structured_n;
        } else {
            ++stats.structured_factor;
        }
    } else {
        if (g == 1) {
            ++stats.random_one;
        } else if (g == N) {
            ++stats.random_n;
        } else {
            ++stats.random_factor;
        }
    }
}

void run_case(
    u64 p_factor,
    u64 q_factor,
    Stats& stats,
    std::mt19937_64& rng
) {
    const u128 N =
        static_cast<u128>(p_factor) *
        static_cast<u128>(q_factor);

    const u128 s =
        isqrt_u128(N);

    const u64 bases[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    const u128 n_values[] = {
        N,
        s,
        s + 1,
        N - s * s
    };

    /*
        Preserve exactly the Experiment 290/291
        expression-generation strategy.
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

        for (u128 n : n_values) {
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

                    const i128 structured =
                        delta1(
                            m,
                            r,
                            n,
                            base
                        );

                    record_signed_value(
                        structured,
                        N,
                        true,
                        stats
                    );

                    const i128 random_value =
                        static_cast<i128>(
                            static_cast<u128>(
                                rng()
                            ) %
                            N
                        );

                    record_signed_value(
                        random_value,
                        N,
                        false,
                        stats
                    );

                    ++stats.expressions;
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

                        const i128 structured =
                            delta2(
                                m,
                                r,
                                ss,
                                n,
                                base
                            );

                        record_signed_value(
                            structured,
                            N,
                            true,
                            stats
                        );

                        const i128 random_value =
                            static_cast<i128>(
                                static_cast<u128>(
                                    rng()
                                ) %
                                N
                            );

                        record_signed_value(
                            random_value,
                            N,
                            false,
                            stats
                        );

                        ++stats.expressions;
                    }
                }
            }
        }
    }
}

void print_case_result(
    std::size_t index,
    u64 p,
    u64 q,
    const Stats& stats
) {
    const double structured_rate =
        stats.expressions == 0
            ? 0.0
            : static_cast<double>(
                  stats.structured_factor
              ) /
              static_cast<double>(
                  stats.expressions
              );

    const double random_rate =
        stats.expressions == 0
            ? 0.0
            : static_cast<double>(
                  stats.random_factor
              ) /
              static_cast<double>(
                  stats.expressions
              );

    const double enrichment =
        random_rate == 0.0
            ? 0.0
            : structured_rate /
              random_rate;

    std::cout
        << "CASE "
        << index
        << " N="
        << p
        << "*"
        << q
        << '\n';

    std::cout
        << "  expressions="
        << stats.expressions
        << '\n';

    std::cout
        << "  structured_factor_hits="
        << stats.structured_factor
        << '\n';

    std::cout
        << "  random_factor_hits="
        << stats.random_factor
        << '\n';

    std::cout
        << "  structured_hit_rate="
        << structured_rate
        << '\n';

    std::cout
        << "  random_hit_rate="
        << random_rate
        << '\n';

    std::cout
        << "  enrichment="
        << enrichment
        << '\n';

    std::cout
        << "  structured_gcd_N="
        << stats.structured_n
        << '\n';

    std::cout
        << "  random_gcd_N="
        << stats.random_n
        << '\n';

    std::cout
        << "  structured_gcd_1="
        << stats.structured_one
        << '\n';

    std::cout
        << "  random_gcd_1="
        << stats.random_one
        << "\n\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 292\n";

    std::cout
        << "SIGNED DIGIT-CALCULUS GCD BASELINE\n";

    std::cout
        << "FIXED D1/D2 UNSIGNED WRAP\n";

    std::cout
        << "MATCHED RANDOM CONTROL\n\n";

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

    std::mt19937_64 rng(
        0x292292ULL
    );

    Stats total;

    for (std::size_t i = 0;
         i < cases.size();
         ++i) {

        Stats stats;

        run_case(
            cases[i].first,
            cases[i].second,
            stats,
            rng
        );

        print_case_result(
            i,
            cases[i].first,
            cases[i].second,
            stats
        );

        total.expressions +=
            stats.expressions;

        total.structured_one +=
            stats.structured_one;

        total.structured_n +=
            stats.structured_n;

        total.structured_factor +=
            stats.structured_factor;

        total.random_one +=
            stats.random_one;

        total.random_n +=
            stats.random_n;

        total.random_factor +=
            stats.random_factor;
    }

    const double structured_rate =
        static_cast<double>(
            total.structured_factor
        ) /
        static_cast<double>(
            total.expressions
        );

    const double random_rate =
        static_cast<double>(
            total.random_factor
        ) /
        static_cast<double>(
            total.expressions
        );

    const double enrichment =
        random_rate == 0.0
            ? 0.0
            : structured_rate /
              random_rate;

    std::cout
        << "SUMMARY\n";

    std::cout
        << "semiprimes="
        << cases.size()
        << '\n';

    std::cout
        << "expression_pairs="
        << total.expressions
        << '\n';

    std::cout
        << "structured_factor_hits="
        << total.structured_factor
        << '\n';

    std::cout
        << "random_factor_hits="
        << total.random_factor
        << '\n';

    std::cout
        << "structured_hit_rate="
        << structured_rate
        << '\n';

    std::cout
        << "random_hit_rate="
        << random_rate
        << '\n';

    std::cout
        << "enrichment="
        << enrichment
        << '\n';

    std::cout
        << "structured_gcd_N="
        << total.structured_n
        << '\n';

    std::cout
        << "random_gcd_N="
        << total.random_n
        << '\n';

    std::cout
        << "structured_gcd_1="
        << total.structured_one
        << '\n';

    std::cout
        << "random_gcd_1="
        << total.random_one
        << '\n';

    std::cout
        << "\nFINISHED EXPERIMENT 292\n";

    return 0;
}
