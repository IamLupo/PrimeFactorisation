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

struct Stats {
    std::size_t expressions = 0;

    std::size_t gcd1 = 0;
    std::size_t gcdN = 0;
    std::size_t nontrivial = 0;

    std::size_t p_only = 0;
    std::size_t q_only = 0;
    std::size_t both = 0;

    std::size_t d1_p_only = 0;
    std::size_t d1_q_only = 0;

    std::size_t d2_p_only = 0;
    std::size_t d2_q_only = 0;

    std::size_t printed = 0;
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

    return hi <= n / hi ? hi : lo;
}

/*
    Simple primality test.

    The generated primes are only <= 200000,
    so trial division is more than sufficient.
*/
bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3;
         d <= n / d;
         d += 2) {

        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

std::vector<u64> generate_primes(
    u64 limit
) {
    std::vector<bool> composite(
        limit + 1,
        false
    );

    std::vector<u64> primes;

    for (u64 i = 2;
         i <= limit;
         ++i) {

        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (i <= limit / i) {
            for (u64 j = i * i;
                 j <= limit;
                 j += i) {

                composite[j] = true;
            }
        }
    }

    return primes;
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
    Correct closed MISS prefix formula.
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
    u64 base
) {
    if (n == 0) {
        return 0;
    }

    const u128 misses =
        miss_prefix(
            m,
            n - 1,
            base
        );

    return static_cast<i128>(
        n - misses
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
        H(m, n, base);

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

u64 gcd_signed(
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

void classify(
    i128 value,
    char type,
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
    ++stats.expressions;

    const u64 g =
        gcd_signed(
            value,
            N
        );

    if (g == 1) {
        ++stats.gcd1;
        return;
    }

    if (g == static_cast<u64>(N)) {
        ++stats.gcdN;
        return;
    }

    ++stats.nontrivial;

    const bool divisible_p =
        g == p;

    const bool divisible_q =
        g == q;

    if (divisible_p && divisible_q) {
        ++stats.both;
    } else if (divisible_p) {
        ++stats.p_only;
    } else if (divisible_q) {
        ++stats.q_only;
    }

    if (type == '1') {
        if (divisible_p) {
            ++stats.d1_p_only;
        }

        if (divisible_q) {
            ++stats.d1_q_only;
        }
    } else {
        if (divisible_p) {
            ++stats.d2_p_only;
        }

        if (divisible_q) {
            ++stats.d2_q_only;
        }
    }

    /*
        Print only the first ten genuine factor-specific hits.
    */
    if ((divisible_p || divisible_q) &&
        stats.printed < 10) {

        ++stats.printed;

        std::cout
            << "\nFACTOR-SPECIFIC HIT\n";

        std::cout
            << "type=D"
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

        if (type == '2') {
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
            << "gcd="
            << g
            << '\n';

        std::cout
            << "divisible_p="
            << (divisible_p ? 1 : 0)
            << '\n';

        std::cout
            << "divisible_q="
            << (divisible_q ? 1 : 0)
            << '\n';
    }
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

    /*
        We use several m values around sqrt(N),
        plus N itself as n.
    */
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
            N - (m * m <= N
                     ? m * m
                     : N)
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

                    const i128 value =
                        delta1(
                            m,
                            r,
                            n,
                            base
                        );

                    classify(
                        value,
                        '1',
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

                        const i128 value =
                            delta2(
                                m,
                                r,
                                s,
                                n,
                                base
                            );

                        classify(
                            value,
                            '2',
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
        << "START EXPERIMENT 297\n";

    std::cout
        << "VERIFIED PRIME-SEMIPRIME FACTOR SIGNAL SEARCH\n";

    std::cout
        << "GENUINE p*q ONLY\n";

    std::cout
        << "D1/D2 CRT ASYMMETRY TEST\n\n";

    /*
        Prime pool.
    */
    const std::vector<u64> primes =
        generate_primes(
            200000
        );

    /*
        Deterministic selection of genuinely prime
        balanced semiprimes.
    */
    std::vector<std::pair<u64, u64>> cases;

    const std::size_t wanted_cases = 30;

    std::size_t index = 0;

    while (cases.size() < wanted_cases) {
        const std::size_t i =
            (index * 7919 + 17) %
            primes.size();

        const std::size_t j =
            (index * 104729 + 101) %
            primes.size();

        ++index;

        const u64 p =
            primes[i];

        const u64 q =
            primes[j];

        if (p == q) {
            continue;
        }

        const u64 lo =
            std::min(p, q);

        const u64 hi =
            std::max(p, q);

        /*
            Keep the factors reasonably balanced.
        */
        if (hi > lo * 3) {
            continue;
        }

        bool duplicate = false;

        for (const auto& item : cases) {
            if (item.first == lo &&
                item.second == hi) {

                duplicate = true;
                break;
            }
        }

        if (duplicate) {
            continue;
        }

        cases.push_back({
            lo,
            hi
        });
    }

    std::cout
        << "verified_semiprimes="
        << cases.size()
        << '\n';

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
            << '\n';

        std::cout
            << "p="
            << cases[i].first
            << " q="
            << cases[i].second
            << '\n';

        std::cout
            << "N="
            << to_string_u128(
                   static_cast<u128>(
                       cases[i].first
                   ) *
                   static_cast<u128>(
                       cases[i].second
                   )
               )
            << '\n';

        std::cout
            << "expressions="
            << local.expressions
            << '\n';

        std::cout
            << "nontrivial="
            << local.nontrivial
            << '\n';

        std::cout
            << "p_only="
            << local.p_only
            << '\n';

        std::cout
            << "q_only="
            << local.q_only
            << '\n';

        std::cout
            << "both="
            << local.both
            << '\n';

        std::cout
            << "D1_p_only="
            << local.d1_p_only
            << '\n';

        std::cout
            << "D1_q_only="
            << local.d1_q_only
            << '\n';

        std::cout
            << "D2_p_only="
            << local.d2_p_only
            << '\n';

        std::cout
            << "D2_q_only="
            << local.d2_q_only
            << '\n';

        total.expressions +=
            local.expressions;

        total.gcd1 +=
            local.gcd1;

        total.gcdN +=
            local.gcdN;

        total.nontrivial +=
            local.nontrivial;

        total.p_only +=
            local.p_only;

        total.q_only +=
            local.q_only;

        total.both +=
            local.both;

        total.d1_p_only +=
            local.d1_p_only;

        total.d1_q_only +=
            local.d1_q_only;

        total.d2_p_only +=
            local.d2_p_only;

        total.d2_q_only +=
            local.d2_q_only;
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "expressions="
        << total.expressions
        << '\n';

    std::cout
        << "gcd1="
        << total.gcd1
        << '\n';

    std::cout
        << "gcdN="
        << total.gcdN
        << '\n';

    std::cout
        << "nontrivial="
        << total.nontrivial
        << '\n';

    std::cout
        << "p_only="
        << total.p_only
        << '\n';

    std::cout
        << "q_only="
        << total.q_only
        << '\n';

    std::cout
        << "both="
        << total.both
        << '\n';

    std::cout
        << "D1_p_only="
        << total.d1_p_only
        << '\n';

    std::cout
        << "D1_q_only="
        << total.d1_q_only
        << '\n';

    std::cout
        << "D2_p_only="
        << total.d2_p_only
        << '\n';

    std::cout
        << "D2_q_only="
        << total.d2_q_only
        << '\n';

    std::cout
        << "\nFINISHED EXPERIMENT 297\n";

    return 0;
}
