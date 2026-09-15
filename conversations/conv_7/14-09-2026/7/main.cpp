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

struct Factor {
    u64 p;
    unsigned exponent;
};

struct Stats {
    std::size_t expressions = 0;
    std::size_t nontrivial = 0;

    std::size_t matched_original_p = 0;
    std::size_t matched_original_q = 0;
    std::size_t matched_other_prime = 0;

    std::size_t gcd_N = 0;
    std::size_t gcd_1 = 0;

    std::vector<std::pair<u64, std::size_t>>
        divisor_hits;
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

/*
    Complete trial-division factorization.

    The current test N values are only around 10^10,
    so this is easily fast enough.
*/
std::vector<Factor> factorize(u64 n) {
    std::vector<Factor> factors;

    if (n % 2 == 0) {
        unsigned e = 0;

        while (n % 2 == 0) {
            n /= 2;
            ++e;
        }

        factors.push_back({2, e});
    }

    for (u64 p = 3;
         p <= n / p;
         p += 2) {

        if (n % p != 0) {
            continue;
        }

        unsigned e = 0;

        while (n % p == 0) {
            n /= p;
            ++e;
        }

        factors.push_back({p, e});
    }

    if (n > 1) {
        factors.push_back({n, 1});
    }

    return factors;
}

bool is_prime_u64(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 p = 3;
         p <= n / p;
         p += 2) {

        if (n % p == 0) {
            return false;
        }
    }

    return true;
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
        H(m, n, p);

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
    i128 value,
    u128 N
) {
    const u128 magnitude =
        abs_i128(value);

    return std::gcd(
        static_cast<u64>(
            magnitude % N
        ),
        static_cast<u64>(N)
    );
}

void add_divisor_hit(
    Stats& stats,
    u64 divisor
) {
    for (auto& item :
         stats.divisor_hits) {

        if (item.first == divisor) {
            ++item.second;
            return;
        }
    }

    stats.divisor_hits.push_back(
        {divisor, 1}
    );
}

void classify_expression(
    const char* type,
    i128 value,
    u64 original_p,
    u64 original_q,
    const std::vector<Factor>& factors,
    u128 N,
    u64 base,
    u128 m,
    u128 n,
    unsigned r,
    unsigned s,
    Stats& stats
) {
    ++stats.expressions;

    const u64 g =
        gcd_value(
            value,
            N
        );

    if (g == 1) {
        ++stats.gcd_1;
        return;
    }

    if (g == static_cast<u64>(N)) {
        ++stats.gcd_N;
        return;
    }

    ++stats.nontrivial;

    bool matches_original_p =
        g == original_p;

    bool matches_original_q =
        g == original_q;

    if (matches_original_p) {
        ++stats.matched_original_p;
    }

    if (matches_original_q) {
        ++stats.matched_original_q;
    }

    bool is_actual_prime_factor =
        false;

    for (const Factor& f : factors) {
        if (g == f.p) {
            is_actual_prime_factor = true;
            break;
        }
    }

    if (!is_actual_prime_factor) {
        /*
            This is a proper divisor that is not itself a
            prime factor.
        */
        add_divisor_hit(
            stats,
            g
        );
    }

    std::cout
        << "\nNONTRIVIAL GCD\n";

    std::cout
        << "type="
        << type
        << '\n';

    std::cout
        << "original_p="
        << original_p
        << '\n';

    std::cout
        << "original_q="
        << original_q
        << '\n';

    std::cout
        << "actual_gcd="
        << g
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
        << "is_prime_factor="
        << (is_actual_prime_factor ? 1 : 0)
        << '\n';

    std::cout
        << "matches_original_p="
        << (matches_original_p ? 1 : 0)
        << '\n';

    std::cout
        << "matches_original_q="
        << (matches_original_q ? 1 : 0)
        << '\n';

    /*
        Do not print endlessly.
        This still records every classification in stats.
    */
    static std::size_t printed = 0;

    ++printed;

    if (printed >= 50) {
        /*
            Suppress detailed printing after the first 50.
        */
        return;
    }
}

void scan_case(
    u64 original_p,
    u64 original_q,
    Stats& stats
) {
    const u128 N =
        static_cast<u128>(original_p) *
        static_cast<u128>(original_q);

    const auto factors =
        factorize(
            static_cast<u64>(N)
        );

    std::cout
        << "\n============================\n";

    std::cout
        << "N="
        << to_string_u128(N)
        << '\n';

    std::cout
        << "original_p="
        << original_p
        << " prime="
        << (
            is_prime_u64(
                original_p
            )
                ? 1
                : 0
        )
        << '\n';

    std::cout
        << "original_q="
        << original_q
        << " prime="
        << (
            is_prime_u64(
                original_q
            )
                ? 1
                : 0
        )
        << '\n';

    std::cout
        << "factorization=";

    for (std::size_t i = 0;
         i < factors.size();
         ++i) {

        if (i > 0) {
            std::cout
                << " * ";
        }

        std::cout
            << factors[i].p;

        if (factors[i].exponent > 1) {
            std::cout
                << "^"
                << factors[i].exponent;
        }
    }

    std::cout
        << '\n';

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

                /*
                    D1.
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

                    classify_expression(
                        "D1",
                        value,
                        original_p,
                        original_q,
                        factors,
                        N,
                        base,
                        m,
                        n,
                        r,
                        999,
                        stats
                    );
                }

                /*
                    D2.
                */
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

                        classify_expression(
                            "D2",
                            value,
                            original_p,
                            original_q,
                            factors,
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
        << "START EXPERIMENT 296\n";

    std::cout
        << "COMPLETE FACTORIZATION OF ALL TEST N\n";

    std::cout
        << "RECLASSIFYING ALL NONTRIVIAL GCDS\n";

    std::cout
        << "REPRODUCING EXPERIMENT 294 WITH FULL FACTOR DATA\n\n";

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

    Stats total;

    for (const auto& pair : cases) {
        Stats local;

        scan_case(
            pair.first,
            pair.second,
            local
        );

        std::cout
            << "\nCASE SUMMARY\n";

        std::cout
            << "original_p="
            << pair.first
            << '\n';

        std::cout
            << "original_q="
            << pair.second
            << '\n';

        std::cout
            << "expressions="
            << local.expressions
            << '\n';

        std::cout
            << "gcd_1="
            << local.gcd_1
            << '\n';

        std::cout
            << "gcd_N="
            << local.gcd_N
            << '\n';

        std::cout
            << "nontrivial="
            << local.nontrivial
            << '\n';

        std::cout
            << "matched_original_p="
            << local.matched_original_p
            << '\n';

        std::cout
            << "matched_original_q="
            << local.matched_original_q
            << '\n';

        total.expressions +=
            local.expressions;

        total.gcd_1 +=
            local.gcd_1;

        total.gcd_N +=
            local.gcd_N;

        total.nontrivial +=
            local.nontrivial;

        total.matched_original_p +=
            local.matched_original_p;

        total.matched_original_q +=
            local.matched_original_q;

        for (const auto& item :
             local.divisor_hits) {

            add_divisor_hit(
                total,
                item.first
            );

            for (std::size_t k = 1;
                 k < item.second;
                 ++k) {

                add_divisor_hit(
                    total,
                    item.first
                );
            }
        }
    }

    std::cout
        << "\n============================\n";

    std::cout
        << "TOTAL\n";

    std::cout
        << "expressions="
        << total.expressions
        << '\n';

    std::cout
        << "gcd_1="
        << total.gcd_1
        << '\n';

    std::cout
        << "gcd_N="
        << total.gcd_N
        << '\n';

    std::cout
        << "nontrivial="
        << total.nontrivial
        << '\n';

    std::cout
        << "matched_original_p="
        << total.matched_original_p
        << '\n';

    std::cout
        << "matched_original_q="
        << total.matched_original_q
        << '\n';

    std::cout
        << "\nNON-ORIGINAL PROPER DIVISORS\n";

    for (const auto& item :
         total.divisor_hits) {

        std::cout
            << "divisor="
            << item.first
            << " hits="
            << item.second
            << '\n';
    }

    std::cout
        << "\nFINISHED EXPERIMENT 296\n";

    return 0;
}
