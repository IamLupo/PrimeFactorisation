#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int MAX_DIGITS = 128;

struct Digits {
    std::array<u64, MAX_DIGITS> d{};
    int len = 1;
};

struct Stats {
    std::size_t cases = 0;
    std::size_t failures = 0;

    std::size_t safe64_cases = 0;
    std::size_t safe64_failures = 0;

    std::size_t overflow64_cases = 0;
    std::size_t overflow64_failures = 0;

    std::size_t new_digit_cases = 0;
    std::size_t new_digit_failures = 0;

    std::size_t printed = 0;
};

constexpr u128 U64_MAX_U128 =
    static_cast<u128>(~static_cast<u64>(0));

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

u64 get_digit(
    const Digits& d,
    int i
) {
    if (i < 0 || i >= d.len) {
        return 0;
    }

    return d.d[i];
}

u64 digit_at(
    const Digits& d,
    unsigned r
) {
    return get_digit(
        d,
        static_cast<int>(r)
    );
}

u128 total_miss(
    const Digits& md
) {
    u128 result = 1;

    for (int i = 0; i < md.len; ++i) {
        result *=
            static_cast<u128>(
                md.d[i] + 1
            );
    }

    return result;
}

/*
    Independent MSB-first digit DP reference.

        M_m(y)
        =
        #{x : 0 <= x <= y and x <=_p m}
*/
u128 miss_prefix_dp(
    u128 m,
    u128 y,
    u64 p
) {
    if (y >= m) {
        return total_miss(
            to_base_p(m, p)
        );
    }

    const Digits md =
        to_base_p(m, p);

    const Digits yd =
        to_base_p(y, p);

    const int len =
        std::max(
            md.len,
            yd.len
        );

    u128 less = 0;
    u128 tight = 1;

    for (int i = len - 1;
         i >= 0;
         --i) {

        const u64 mi =
            get_digit(md, i);

        const u64 yi =
            get_digit(yd, i);

        /*
            Already below y:
            any digit 0..mi is allowed.
        */
        const u128 from_less =
            less *
            static_cast<u128>(
                mi + 1
            );

        /*
            Tight prefix becomes less:

                x_i < y_i

            with x_i <= m_i.

            Number of choices:

                min(m_i+1, y_i)
        */
        const u128 from_tight =
            tight *
            static_cast<u128>(
                std::min<u64>(
                    mi + 1,
                    yi
                )
            );

        /*
            Remain tight with x_i = y_i.
        */
        const u128 next_tight =
            (yi <= mi)
                ? tight
                : 0;

        less =
            from_less +
            from_tight;

        tight =
            next_tight;
    }

    return less + tight;
}

/*
    Correct closed MISS-prefix formula validated
    in Experiment 288.

    For y < m, let h be the highest digit where
    y_h != m_h.

        W_i = product_{j<i}(m_j+1)

        M_m(y)
          =
          sum_{i>h} y_i W_i
          +
          y_h W_h
          +
          M_{m_<h}(y_<h)
*/
u128 miss_prefix_closed(
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

    const int len =
        std::max(
            md.len,
            yd.len
        );

    int h = -1;

    for (int i = len - 1;
         i >= 0;
         --i) {

        if (get_digit(md, i) !=
            get_digit(yd, i)) {

            h = i;
            break;
        }
    }

    if (h < 0) {
        return 1;
    }

    std::array<u128, MAX_DIGITS + 1> weight{};

    weight[0] = 1;

    for (int i = 0;
         i < len;
         ++i) {

        weight[i + 1] =
            weight[i] *
            static_cast<u128>(
                get_digit(md, i) + 1
            );
    }

    u128 result = 0;

    /*
        Complete lower blocks from higher digits.
    */
    for (int i = len - 1;
         i > h;
         --i) {

        result +=
            static_cast<u128>(
                get_digit(yd, i)
            ) *
            weight[i];
    }

    /*
        At h:

            x_h = 0,...,y_h-1

        gives y_h complete lower blocks.
    */
    result +=
        static_cast<u128>(
            get_digit(yd, h)
        ) *
        weight[h];

    /*
        Equal h digit:
        lower digits must satisfy both bounds.
    */
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
        miss_prefix_dp(
            low_m,
            low_y,
            p
        );

    return result;
}

/*
    Exact local delta from independent digit DP:

        M_{m+p^r}(n-1)
        -
        M_m(n-1)
*/
u128 exact_local_delta(
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
        miss_prefix_dp(
            current,
            y,
            p
        );

    const u128 after =
        miss_prefix_dp(
            next,
            y,
            p
        );

    /*
        Diagnostic guard.
    */
    if (after < before) {
        return static_cast<u128>(-1);
    }

    return after - before;
}

/*
    Closed local-slice formula.

    Update:

        m' = m + p^r

    with m_r < p-1.

    Newly created MISS points have:

        x_r = m_r + 1

    and every other digit is bounded by m.

    Split

        y = high_y * p^(r+1)
            + digit_y * p^r
            + low_y

    and similarly for m.
*/
u128 local_slice_closed(
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
        power_u128(
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

    /*
        No-carry condition.
    */
    if (mr + 1 >= p) {
        return 0;
    }

    const u64 new_digit =
        mr + 1;

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

    /*
        Number of admissible lower MISS values.
    */
    const u128 lower_total =
        miss_prefix_closed(
            low_m,
            low_m,
            p
        );

    u128 result = 0;

    /*
        high_x < high_y.

        Every such high MISS prefix gets a complete
        lower block.
    */
    if (high_y > 0) {
        const u128 high_count =
            miss_prefix_closed(
                high_m,
                high_y - 1,
                p
            );

        result +=
            high_count *
            lower_total;
    }

    /*
        Determine whether high_y itself is admissible.
    */
    const u128 high_through =
        miss_prefix_closed(
            high_m,
            high_y,
            p
        );

    const u128 high_before =
        (high_y == 0)
            ? 0
            : miss_prefix_closed(
                high_m,
                high_y - 1,
                p
            );

    const bool high_admissible =
        high_through > high_before;

    if (!high_admissible) {
        return result;
    }

    /*
        Compare the fixed new digit with y_r.
    */
    if (digit_y < new_digit) {
        return result;
    }

    if (digit_y > new_digit) {
        return result + lower_total;
    }

    /*
        digit_y == new_digit:

        partial lower block.
    */
    result +=
        miss_prefix_closed(
            low_m,
            low_y,
            p
        );

    return result;
}

bool legal_update(
    u128 current,
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
            (current / pr) %
            static_cast<u128>(p)
        );

    return mr + 1 < p;
}

void record_case(
    const char* group_name,
    u64 p,
    u128 current,
    unsigned r,
    u128 n,
    Stats& stats
) {
    if (!legal_update(
            current,
            r,
            p
        )) {
        return;
    }

    const u128 step =
        power_u128(
            static_cast<u128>(p),
            r
        );

    const u128 next =
        current + step;

    const Digits md =
        to_base_p(current, p);

    const bool new_digit =
        r >= static_cast<unsigned>(md.len);

    const bool overflow64 =
        next > U64_MAX_U128;

    ++stats.cases;

    if (overflow64) {
        ++stats.overflow64_cases;
    } else {
        ++stats.safe64_cases;
    }

    if (new_digit) {
        ++stats.new_digit_cases;
    }

    const u128 exact =
        exact_local_delta(
            current,
            next,
            n,
            p
        );

    const u128 predicted =
        local_slice_closed(
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
        ++stats.overflow64_failures;
    } else {
        ++stats.safe64_failures;
    }

    if (new_digit) {
        ++stats.new_digit_failures;
    }

    if (stats.printed >= 3) {
        return;
    }

    ++stats.printed;

    const u64 mr =
        digit_at(
            md,
            r
        );

    std::cout
        << "\nFAIL "
        << group_name
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
        << "m_r="
        << mr
        << '\n';

    std::cout
        << "new_digit_case="
        << (new_digit ? 1 : 0)
        << '\n';

    std::cout
        << "next_over_u64="
        << (overflow64 ? 1 : 0)
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
    const char* group_name,
    const u64* primes,
    int prime_count,
    std::size_t target,
    std::uint64_t seed,
    Stats& stats
) {
    std::mt19937_64 rng(seed);

    while (stats.cases < target) {
        const u64 p =
            primes[
                rng() %
                static_cast<std::size_t>(
                    prime_count
                )
            ];

        /*
            r extends into and beyond the
            64-bit transition range.
        */
        const unsigned r =
            static_cast<unsigned>(
                rng() % 14
            );

        const u128 current =
            static_cast<u128>(
                rng()
            ) %
            static_cast<u128>(
                1000000000000000000ULL
            );

        const u128 n =
            static_cast<u128>(
                rng()
            ) %
            static_cast<u128>(
                1000000000000000000ULL
            );

        record_case(
            group_name,
            p,
            current,
            r,
            n,
            stats
        );
    }
}

void print_stats(
    const char* name,
    const Stats& stats
) {
    std::cout
        << '\n'
        << name
        << '\n';

    std::cout
        << "cases="
        << stats.cases
        << " fail="
        << stats.failures
        << '\n';

    std::cout
        << "safe64_cases="
        << stats.safe64_cases
        << " safe64_fail="
        << stats.safe64_failures
        << '\n';

    std::cout
        << "overflow64_cases="
        << stats.overflow64_cases
        << " overflow64_fail="
        << stats.overflow64_failures
        << '\n';

    std::cout
        << "new_digit_cases="
        << stats.new_digit_cases
        << " new_digit_fail="
        << stats.new_digit_failures
        << '\n';
}

void explicit_edge_cases(
    Stats& stats
) {
    const u64 primes[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31, 37, 41, 43
    };

    for (u64 p : primes) {
        for (unsigned r = 0;
             r <= 13;
             ++r) {

            const u128 pr =
                power_u128(
                    static_cast<u128>(p),
                    r
                );

            const u128 currents[] = {
                0,
                1,
                pr,
                pr + 1,
                2 * pr,
                3 * pr + 7
            };

            for (u128 current : currents) {
                if (!legal_update(
                        current,
                        r,
                        p
                    )) {
                    continue;
                }

                const u128 step = pr;
                const u128 next =
                    current + step;

                const u128 ns[] = {
                    0,
                    1,
                    current,
                    current + 1,
                    next,
                    next + 1,
                    step,
                    step + 1
                };

                for (u128 n : ns) {
                    if (stats.cases >=
                        10000) {
                        return;
                    }

                    record_case(
                        "EDGE",
                        p,
                        current,
                        r,
                        n,
                        stats
                    );
                }
            }
        }
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 289\n";

    std::cout
        << "LARGE ARBITRARY-CURRENT LOCAL DIGIT SLICE\n";

    std::cout
        << "CORRECTED PREFIX FORMULA VS INDEPENDENT DIGIT DP\n";

    std::cout
        << "SAFE64 / OVERFLOW64 / NEW-DIGIT DIAGNOSTICS\n\n";

    const u64 low_primes[] = {
        2, 3, 5, 7, 11, 13
    };

    const u64 mid_primes[] = {
        17, 19, 23, 29, 31
    };

    const u64 high_primes[] = {
        37, 41, 43
    };

    Stats low;
    Stats mid;
    Stats high;
    Stats edge;

    explicit_edge_cases(edge);

    std::cout
        << "edge_cases="
        << edge.cases
        << " fail="
        << edge.failures
        << '\n';

    run_group(
        "LOW",
        low_primes,
        static_cast<int>(
            sizeof(low_primes) /
            sizeof(low_primes[0])
        ),
        100000,
        0x289001ULL,
        low
    );

    run_group(
        "MID",
        mid_primes,
        static_cast<int>(
            sizeof(mid_primes) /
            sizeof(mid_primes[0])
        ),
        100000,
        0x289002ULL,
        mid
    );

    run_group(
        "HIGH",
        high_primes,
        static_cast<int>(
            sizeof(high_primes) /
            sizeof(high_primes[0])
        ),
        100000,
        0x289003ULL,
        high
    );

    print_stats(
        "EDGE",
        edge
    );

    print_stats(
        "LOW",
        low
    );

    print_stats(
        "MID",
        mid
    );

    print_stats(
        "HIGH",
        high
    );

    const std::size_t total_cases =
        edge.cases +
        low.cases +
        mid.cases +
        high.cases;

    const std::size_t total_failures =
        edge.failures +
        low.failures +
        mid.failures +
        high.failures;

    const std::size_t total_safe_failures =
        edge.safe64_failures +
        low.safe64_failures +
        mid.safe64_failures +
        high.safe64_failures;

    const std::size_t total_overflow_failures =
        edge.overflow64_failures +
        low.overflow64_failures +
        mid.overflow64_failures +
        high.overflow64_failures;

    const std::size_t total_new_digit_failures =
        edge.new_digit_failures +
        low.new_digit_failures +
        mid.new_digit_failures +
        high.new_digit_failures;

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

    std::cout
        << "new_digit_fail="
        << total_new_digit_failures
        << '\n';

    const bool pass =
        total_failures == 0;

    std::cout
        << "\nOVERALL PASS="
        << (pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 289\n";

    return pass ? 0 : 1;
}
