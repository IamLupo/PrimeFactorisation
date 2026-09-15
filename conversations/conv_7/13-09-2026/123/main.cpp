#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u128 U64_MAX_U128 =
    static_cast<u128>(~static_cast<u64>(0));

constexpr int MAX_DIGITS = 128;

struct Digits {
    std::array<u64, MAX_DIGITS> d{};
    int len = 1;
};

struct GroupStats {
    std::size_t cases = 0;
    std::size_t failures = 0;

    std::size_t safe_cases = 0;
    std::size_t safe_failures = 0;

    std::size_t overflow_cases = 0;
    std::size_t overflow_failures = 0;

    std::size_t printed_failures = 0;
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

u128 ipow(u128 base, unsigned exp) {
    u128 result = 1;

    while (exp > 0) {
        if (exp & 1u) {
            result *= base;
        }

        exp >>= 1u;

        if (exp != 0) {
            base *= base;
        }
    }

    return result;
}

Digits to_base_p(
    u128 x,
    u64 p
) {
    Digits out;

    if (x == 0) {
        out.d[0] = 0;
        out.len = 1;
        return out;
    }

    int i = 0;

    while (x > 0) {
        out.d[i++] =
            static_cast<u64>(
                x % static_cast<u128>(p)
            );

        x /= static_cast<u128>(p);
    }

    out.len = i;
    return out;
}

u64 digit_at(
    const Digits& d,
    unsigned r
) {
    if (r >= static_cast<unsigned>(d.len)) {
        return 0;
    }

    return d.d[r];
}

u128 total_miss_count(
    const Digits& md
) {
    u128 result = 1;

    for (int i = 0; i < md.len; ++i) {
        result *=
            static_cast<u128>(md.d[i] + 1);
    }

    return result;
}

/*
    Exact prefix count:

        M_m(y)
        =
        #{x : 0 <= x <= y and x <=_p m}
*/
u128 miss_prefix(
    u128 m,
    u128 y,
    u64 p
) {
    const Digits md =
        to_base_p(m, p);

    if (y >= m) {
        return total_miss_count(md);
    }

    const Digits yd =
        to_base_p(y, p);

    const int len =
        std::max(md.len, yd.len);

    std::array<u128, MAX_DIGITS + 1> weight{};
    weight[0] = 1;

    for (int i = 0; i < len; ++i) {
        const u64 mi =
            i < md.len ? md.d[i] : 0;

        weight[i + 1] =
            weight[i] *
            static_cast<u128>(mi + 1);
    }

    int h = -1;

    for (int i = len - 1; i >= 0; --i) {
        const u64 mi =
            i < md.len ? md.d[i] : 0;

        const u64 yi =
            i < yd.len ? yd.d[i] : 0;

        if (mi != yi) {
            h = i;
            break;
        }
    }

    if (h < 0) {
        return 1;
    }

    u128 result = 0;

    for (int i = len - 1; i > h; --i) {
        const u64 yi =
            i < yd.len ? yd.d[i] : 0;

        result +=
            static_cast<u128>(yi) *
            weight[i];
    }

    const u64 yh =
        h < yd.len ? yd.d[h] : 0;

    result +=
        static_cast<u128>(yh + 1) *
        weight[h];

    return result;
}

/*
    Exact reference for the local MISS delta.
*/
u128 exact_delta(
    u128 current,
    u128 next,
    u128 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    const u128 y = n - 1;

    const u128 before =
        miss_prefix(
            current,
            y,
            p
        );

    const u128 after =
        miss_prefix(
            next,
            y,
            p
        );

    if (after < before) {
        std::cerr
            << "INTERNAL MONOTONICITY ERROR\n";
        std::abort();
    }

    return after - before;
}

/*
    Candidate local-slice formula under investigation.
*/
u128 local_slice_formula(
    u128 current,
    unsigned r,
    u128 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    const u128 y = n - 1;

    const u128 pr =
        ipow(
            static_cast<u128>(p),
            r
        );

    const u128 pr1 =
        pr *
        static_cast<u128>(p);

    const Digits md =
        to_base_p(current, p);

    const u64 mr =
        digit_at(md, r);

    if (mr + 1 >= p) {
        return 0;
    }

    const u128 high_m =
        current / pr1;

    const u128 high_y =
        y / pr1;

    const u64 digit_y =
        static_cast<u64>(
            (y / pr) %
            static_cast<u128>(p)
        );

    const u128 low_m =
        current % pr;

    const u128 low_y =
        y % pr;

    const u128 lower_total =
        total_miss_count(
            to_base_p(low_m, p)
        );

    u128 result = 0;

    /*
        Higher prefixes strictly below high_y.
    */
    if (high_y > 0) {
        const u128 high_count =
            miss_prefix(
                high_m,
                high_y - 1,
                p
            );

        result +=
            high_count *
            lower_total;
    }

    /*
        Is high_y itself admissible?
    */
    const Digits hm =
        to_base_p(high_m, p);

    const Digits hy =
        to_base_p(high_y, p);

    const int high_len =
        std::max(hm.len, hy.len);

    for (int i = 0; i < high_len; ++i) {
        const u64 a =
            i < hm.len ? hm.d[i] : 0;

        const u64 b =
            i < hy.len ? hy.d[i] : 0;

        if (b > a) {
            return result;
        }
    }

    const u64 new_digit =
        mr + 1;

    if (digit_y < new_digit) {
        return result;
    }

    if (digit_y > new_digit) {
        return result + lower_total;
    }

    result +=
        miss_prefix(
            low_m,
            low_y,
            p
        );

    return result;
}

const char* group_name(
    int group
) {
    switch (group) {
        case 0:
            return "LOW";
        case 1:
            return "MID";
        default:
            return "HIGH";
    }
}

void record_case(
    int group,
    u64 p,
    u128 current,
    unsigned r,
    u128 n,
    GroupStats& stats
) {
    const Digits md =
        to_base_p(current, p);

    const u64 mr =
        digit_at(md, r);

    /*
        Skip carry cases.
    */
    if (mr + 1 >= p) {
        return;
    }

    const u128 step =
        ipow(
            static_cast<u128>(p),
            r
        );

    const u128 next =
        current + step;

    const bool overflow64 =
        next > U64_MAX_U128;

    ++stats.cases;

    if (overflow64) {
        ++stats.overflow_cases;
    } else {
        ++stats.safe_cases;
    }

    const u128 exact =
        exact_delta(
            current,
            next,
            n,
            p
        );

    const u128 predicted =
        local_slice_formula(
            current,
            r,
            n,
            p
        );

    if (exact == predicted) {
        return;
    }

    ++stats.failures;

    if (overflow64) {
        ++stats.overflow_failures;
    } else {
        ++stats.safe_failures;
    }

    /*
        Only print the first three failures
        for this group.
    */
    if (stats.printed_failures >= 3) {
        return;
    }

    ++stats.printed_failures;

    std::cout
        << "\nFAIL "
        << group_name(group)
        << '\n';

    std::cout
        << "p="
        << p
        << '\n';

    std::cout
        << "current="
        << to_string_u128(current)
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
        << "step="
        << to_string_u128(step)
        << '\n';

    std::cout
        << "next="
        << to_string_u128(next)
        << '\n';

    std::cout
        << "next_over_u64="
        << (overflow64 ? 1 : 0)
        << '\n';

    std::cout
        << "m_r="
        << mr
        << '\n';

    std::cout
        << "exact="
        << to_string_u128(exact)
        << '\n';

    std::cout
        << "predicted="
        << to_string_u128(predicted)
        << '\n';
}

void run_group(
    int group,
    const u64* primes,
    int prime_count,
    std::size_t target_cases,
    GroupStats& stats,
    std::uint64_t seed
) {
    std::mt19937_64 rng(seed);

    while (stats.cases < target_cases) {
        const u64 p =
            primes[
                rng() %
                static_cast<std::size_t>(prime_count)
            ];

        /*
            r range intentionally extends into
            the region where p^r may exceed 64 bits.
        */
        const unsigned r =
            static_cast<unsigned>(
                rng() % 14
            );

        const u128 current =
            static_cast<u128>(
                rng() %
                1000000000000000000ULL
            );

        const u128 n =
            static_cast<u128>(
                rng() %
                1000000000000000000ULL
            );

        record_case(
            group,
            p,
            current,
            r,
            n,
            stats
        );
    }
}

void print_group(
    const char* name,
    const GroupStats& s
) {
    std::cout
        << '\n'
        << name
        << '\n';

    std::cout
        << "cases="
        << s.cases
        << " fail="
        << s.failures
        << '\n';

    std::cout
        << "safe64_cases="
        << s.safe_cases
        << " safe64_fail="
        << s.safe_failures
        << '\n';

    std::cout
        << "overflow64_cases="
        << s.overflow_cases
        << " overflow64_fail="
        << s.overflow_failures
        << '\n';
}

int main() {
    std::cout
        << "START EXPERIMENT 283\n";

    std::cout
        << "LOW / MID / HIGH PRIME OVERFLOW DIAGNOSTIC\n";

    std::cout
        << "SAFE-64 VS NEXT-OVER-64 LOCAL DIGIT SLICE\n\n";

    /*
        LOW:
            powers usually remain relatively small.

        MID:
            starts stressing 64-bit arithmetic.

        HIGH:
            p^r crosses 2^64 much earlier.
    */
    const u64 low_primes[] = {
        2, 3, 5, 7, 11, 13
    };

    const u64 mid_primes[] = {
        17, 19, 23, 29, 31
    };

    const u64 high_primes[] = {
        37, 41, 43
    };

    GroupStats low;
    GroupStats mid;
    GroupStats high;

    run_group(
        0,
        low_primes,
        static_cast<int>(
            sizeof(low_primes) /
            sizeof(low_primes[0])
        ),
        10000,
        low,
        0x283001ULL
    );

    run_group(
        1,
        mid_primes,
        static_cast<int>(
            sizeof(mid_primes) /
            sizeof(mid_primes[0])
        ),
        10000,
        mid,
        0x283002ULL
    );

    run_group(
        2,
        high_primes,
        static_cast<int>(
            sizeof(high_primes) /
            sizeof(high_primes[0])
        ),
        10000,
        high,
        0x283003ULL
    );

    print_group("LOW", low);
    print_group("MID", mid);
    print_group("HIGH", high);

    const std::size_t total_cases =
        low.cases +
        mid.cases +
        high.cases;

    const std::size_t total_failures =
        low.failures +
        mid.failures +
        high.failures;

    const std::size_t total_safe_failures =
        low.safe_failures +
        mid.safe_failures +
        high.safe_failures;

    const std::size_t total_overflow_failures =
        low.overflow_failures +
        mid.overflow_failures +
        high.overflow_failures;

    std::cout
        << "\nTOTAL\n";

    std::cout
        << "cases="
        << total_cases
        << " fail="
        << total_failures
        << '\n';

    std::cout
        << "safe64_fail="
        << total_safe_failures
        << '\n';

    std::cout
        << "overflow64_fail="
        << total_overflow_failures
        << '\n';

    /*
        Diagnostic interpretation.

        If:
            safe64_fail == 0
            overflow64_fail > 0

        then the failures strongly correlate with the
        64-bit transition boundary.

        If safe64_fail > 0, the problem is not solely
        caused by current+p^r overflowing uint64.
    */
    const bool pass =
        total_failures == 0;

    std::cout
        << "\nOVERALL PASS="
        << (pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 283\n";

    return pass ? 0 : 1;
}