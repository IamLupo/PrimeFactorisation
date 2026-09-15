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

struct Digits {
    std::vector<u64> d;
};

struct HitExample {
    u64 p;
    u64 q;
    u64 base;

    u128 N;
    u128 m;
    u128 n;

    unsigned r;
    unsigned s;

    char type;

    i128 value;
    i128 second_value;
};

struct Aggregate {
    std::size_t total = 0;

    std::size_t zeroN = 0;
    std::size_t factorP = 0;
    std::size_t factorQ = 0;
    std::size_t unit = 0;

    std::vector<std::pair<u64, std::size_t>> base_hits;
    std::vector<std::pair<unsigned, std::size_t>> r_hits;
    std::vector<std::pair<unsigned, std::size_t>> s_hits;

    std::size_t d1_hits = 0;
    std::size_t d2_hits = 0;
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
        H(
            m,
            n,
            p
        );

    const i128 b =
        H(
            m + pr,
            n,
            p
        );

    const i128 c =
        H(
            m + ps,
            n,
            p
        );

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
        static_cast<u128>(
            -(x + 1)
        ) + 1;
}

u64 gcd_abs_with_N(
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

char classify_gcd(
    u64 g,
    u64 p,
    u64 q,
    u64 N
) {
    if (g == 1) {
        return 'U';
    }

    if (g == N) {
        return 'N';
    }

    if (g == p) {
        return 'P';
    }

    if (g == q) {
        return 'Q';
    }

    return 'X';
}

void add_sorted_count(
    std::vector<std::pair<u64, std::size_t>>& v,
    u64 key
) {
    for (auto& item : v) {
        if (item.first == key) {
            ++item.second;
            return;
        }
    }

    v.push_back({key, 1});
}

void add_sorted_count_u(
    std::vector<std::pair<unsigned, std::size_t>>& v,
    unsigned key
) {
    for (auto& item : v) {
        if (item.first == key) {
            ++item.second;
            return;
        }
    }

    v.push_back({key, 1});
}

void process_hit(
    Aggregate& agg,
    char kind,
    u64 base,
    unsigned r,
    unsigned s
) {
    if (kind == 'N') {
        ++agg.zeroN;
        return;
    }

    if (kind == 'U') {
        ++agg.unit;
        return;
    }

    if (kind == 'P') {
        ++agg.factorP;
    } else if (kind == 'Q') {
        ++agg.factorQ;
    } else {
        return;
    }

    if (kind == 'P' || kind == 'Q') {
        add_sorted_count(
            agg.base_hits,
            base
        );

        add_sorted_count_u(
            agg.r_hits,
            r
        );

        if (s != 999) {
            add_sorted_count_u(
                agg.s_hits,
                s
            );
        }

        if (s == 999) {
            ++agg.d1_hits;
        } else {
            ++agg.d2_hits;
        }
    }
}

void print_example(
    const HitExample& h
) {
    std::cout
        << "\nFACTOR HIT\n";

    std::cout
        << "type="
        << h.type
        << '\n';

    std::cout
        << "p="
        << h.p
        << " q="
        << h.q
        << '\n';

    std::cout
        << "N="
        << to_string_u128(h.N)
        << '\n';

    std::cout
        << "base="
        << h.base
        << '\n';

    std::cout
        << "m="
        << to_string_u128(h.m)
        << '\n';

    std::cout
        << "n="
        << to_string_u128(h.n)
        << '\n';

    std::cout
        << "r="
        << h.r
        << '\n';

    if (h.s != 999) {
        std::cout
            << "s="
            << h.s
            << '\n';
    }

    std::cout
        << "value="
        << to_string_i128(h.value)
        << '\n';

    if (h.s != 999) {
        std::cout
            << "second_value="
            << to_string_i128(
                h.second_value
            )
            << '\n';
    }

    std::cout
        << "actual_factor="
        << (h.type == 'P'
                ? h.p
                : h.q)
        << '\n';
}

void run_semiprime(
    u64 p,
    u64 q,
    std::size_t case_index,
    Aggregate& agg,
    std::vector<HitExample>& examples
) {
    const u128 N =
        static_cast<u128>(p) *
        static_cast<u128>(q);

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
                    D1
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

                    const i128 value =
                        delta1(
                            m,
                            r,
                            n,
                            base
                        );

                    ++agg.total;

                    const u64 g =
                        gcd_abs_with_N(
                            value,
                            N
                        );

                    const char type =
                        classify_gcd(
                            g,
                            p,
                            q,
                            static_cast<u64>(N)
                        );

                    process_hit(
                        agg,
                        type,
                        base,
                        r,
                        999
                    );

                    if ((type == 'P' ||
                         type == 'Q') &&
                        examples.size() < 20) {

                        examples.push_back({
                            p,
                            q,
                            base,
                            N,
                            m,
                            n,
                            r,
                            999,
                            type,
                            value,
                            0
                        });
                    }
                }

                /*
                    D2
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

                        const i128 value =
                            delta2(
                                m,
                                r,
                                ss,
                                n,
                                base
                            );

                        ++agg.total;

                        const u64 g =
                            gcd_abs_with_N(
                                value,
                                N
                            );

                        const char type =
                            classify_gcd(
                                g,
                                p,
                                q,
                                static_cast<u64>(N)
                            );

                        process_hit(
                            agg,
                            type,
                            base,
                            r,
                            ss
                        );

                        if ((type == 'P' ||
                             type == 'Q') &&
                            examples.size() < 20) {

                            examples.push_back({
                                p,
                                q,
                                base,
                                N,
                                m,
                                n,
                                r,
                                ss,
                                type,
                                value,
                                value
                            });
                        }
                    }
                }
            }
        }
    }
}

void print_vector_u64(
    const std::vector<std::pair<u64, std::size_t>>& v
) {
    for (const auto& item : v) {
        std::cout
            << "  "
            << item.first
            << ": "
            << item.second
            << '\n';
    }
}

void print_vector_unsigned(
    const std::vector<std::pair<unsigned, std::size_t>>& v
) {
    for (const auto& item : v) {
        std::cout
            << "  "
            << item.first
            << ": "
            << item.second
            << '\n';
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 293\n";

    std::cout
        << "CLASSIFYING RARE FACTOR HITS\n";

    std::cout
        << "D1/D2 STRUCTURAL PARAMETER ANALYSIS\n\n";

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

    Aggregate total;

    std::vector<HitExample> examples;

    for (std::size_t i = 0;
         i < cases.size();
         ++i) {

        Aggregate local;

        run_semiprime(
            cases[i].first,
            cases[i].second,
            i,
            local,
            examples
        );

        std::cout
            << "CASE "
            << i
            << " N="
            << cases[i].first
            << "*"
            << cases[i].second
            << '\n';

        std::cout
            << "  expressions="
            << local.total
            << '\n';

        std::cout
            << "  factor_p="
            << local.factorP
            << '\n';

        std::cout
            << "  factor_q="
            << local.factorQ
            << '\n';

        std::cout
            << "  factor_total="
            << (
                local.factorP +
                local.factorQ
            )
            << '\n';

        /*
            Merge totals.
        */
        total.total += local.total;
        total.zeroN += local.zeroN;
        total.factorP += local.factorP;
        total.factorQ += local.factorQ;
        total.unit += local.unit;

        total.d1_hits += local.d1_hits;
        total.d2_hits += local.d2_hits;

        for (const auto& item :
             local.base_hits) {

            for (std::size_t k = 0;
                 k < item.second;
                 ++k) {

                add_sorted_count(
                    total.base_hits,
                    item.first
                );
            }
        }

        for (const auto& item :
             local.r_hits) {

            for (std::size_t k = 0;
                 k < item.second;
                 ++k) {

                add_sorted_count_u(
                    total.r_hits,
                    item.first
                );
            }
        }

        for (const auto& item :
             local.s_hits) {

            for (std::size_t k = 0;
                 k < item.second;
                 ++k) {

                add_sorted_count_u(
                    total.s_hits,
                    item.first
                );
            }
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "expressions="
        << total.total
        << '\n';

    std::cout
        << "gcd_N="
        << total.zeroN
        << '\n';

    std::cout
        << "factor_p="
        << total.factorP
        << '\n';

    std::cout
        << "factor_q="
        << total.factorQ
        << '\n';

    std::cout
        << "factor_total="
        << (
            total.factorP +
            total.factorQ
        )
        << '\n';

    std::cout
        << "gcd_1="
        << total.unit
        << '\n';

    std::cout
        << "factor_hits_D1="
        << total.d1_hits
        << '\n';

    std::cout
        << "factor_hits_D2="
        << total.d2_hits
        << '\n';

    std::cout
        << "\nFACTOR HITS BY BASE\n";

    print_vector_u64(
        total.base_hits
    );

    std::cout
        << "\nFACTOR HITS BY r\n";

    print_vector_unsigned(
        total.r_hits
    );

    std::cout
        << "\nFACTOR HITS BY s\n";

    print_vector_unsigned(
        total.s_hits
    );

    std::cout
        << "\nEXAMPLES\n";

    for (const auto& example :
         examples) {

        print_example(example);
    }

    std::cout
        << "\nFINISHED EXPERIMENT 293\n";

    return 0;
}
