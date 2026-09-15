#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr int MAX_DIGITS = 128;

struct Digits {
    std::array<u64, MAX_DIGITS> d{};
    int len = 1;
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
    Independent MSB-first digit DP.

    Counts

        #{x : 0 <= x <= y and x <=_p m}

    State meanings:

        less  = number of prefixes strictly below y
        tight = number of prefixes equal to y

    For a tight state:

        d < y_i:
            number of choices =
            min(m_i + 1, y_i)

        d = y_i:
            one tight choice if y_i <= m_i
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
        const u128 next_less_from_less =
            less *
            static_cast<u128>(
                mi + 1
            );

        /*
            Tight prefix.

            Choose d < yi with d <= mi.

            Number of such d:

                min(mi+1, yi)
        */
        const u128 next_less_from_tight =
            tight *
            static_cast<u128>(
                std::min<u64>(
                    mi + 1,
                    yi
                )
            );

        /*
            Choose d == yi.

            This remains tight iff yi <= mi.
        */
        const u128 next_tight =
            (yi <= mi)
                ? tight
                : 0;

        less =
            next_less_from_less +
            next_less_from_tight;

        tight =
            next_tight;
    }

    return less + tight;
}

/*
    Closed mixed-radix prefix formula.

    For y < m, let h be the most significant
    digit where y_h != m_h.

        M_m(y)
        =
        sum_{i>h} y_i W_i
        +
        (y_h+1) W_h

    W_i = product_{j<i}(m_j+1)
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

    std::array<u128, MAX_DIGITS + 1> weight{};

    weight[0] = 1;

    for (int i = 0;
         i < len;
         ++i) {

        const u64 mi =
            get_digit(md, i);

        weight[i + 1] =
            weight[i] *
            static_cast<u128>(
                mi + 1
            );
    }

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

    u128 result = 0;

    for (int i = len - 1;
         i > h;
         --i) {

        result +=
            static_cast<u128>(
                get_digit(yd, i)
            ) *
            weight[i];
    }

    result +=
        static_cast<u128>(
            get_digit(yd, h) + 1
        ) *
        weight[h];

    return result;
}

u128 exact_local_delta_dp(
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
        Increasing m cannot remove MISS points.
    */
    if (after < before) {
        return static_cast<u128>(-1);
    }

    return after - before;
}

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
        get_digit(
            md,
            static_cast<int>(r)
        );

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
        miss_prefix_dp(
            low_m,
            low_m,
            p
        );

    u128 result = 0;

    /*
        High prefix strictly below high_y.
    */
    if (high_y > 0) {
        result +=
            miss_prefix_dp(
                high_m,
                high_y - 1,
                p
            ) *
            lower_total;
    }

    /*
        Determine whether high_y itself is a MISS
        prefix of high_m.
    */
    const u128 high_before =
        (high_y == 0)
            ? 0
            : miss_prefix_dp(
                high_m,
                high_y - 1,
                p
            );

    const u128 high_through =
        miss_prefix_dp(
            high_m,
            high_y,
            p
        );

    const bool high_admissible =
        high_through > high_before;

    if (!high_admissible) {
        return result;
    }

    const u64 new_digit =
        mr + 1;

    if (digit_y < new_digit) {
        return result;
    }

    if (digit_y > new_digit) {
        return result + lower_total;
    }

    /*
        Equal changed digit:
        only lower MISS points <= low_y.
    */
    result +=
        miss_prefix_dp(
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

struct Stats {
    std::size_t cases = 0;
    std::size_t failures = 0;
    std::size_t printed = 0;
};

void print_prefix_failure(
    u64 p,
    u128 m,
    u128 y,
    u128 dp,
    u128 closed
) {
    std::cout
        << "\nPREFIX FAIL\n";

    std::cout
        << "p="
        << p
        << '\n';

    std::cout
        << "m="
        << to_string_u128(m)
        << '\n';

    std::cout
        << "y="
        << to_string_u128(y)
        << '\n';

    std::cout
        << "dp="
        << to_string_u128(dp)
        << '\n';

    std::cout
        << "closed="
        << to_string_u128(closed)
        << '\n';
}

void print_local_failure(
    u64 p,
    u128 current,
    unsigned r,
    u128 n,
    u128 exact,
    u128 closed
) {
    std::cout
        << "\nLOCAL FAIL\n";

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
        << "exact="
        << to_string_u128(exact)
        << '\n';

    std::cout
        << "closed="
        << to_string_u128(closed)
        << '\n';
}

int main() {
    std::cout
        << "START EXPERIMENT 286\n";

    std::cout
        << "CORRECTED INDEPENDENT DIGIT-DP\n";

    std::cout
        << "PREFIX FORMULA VS DIGIT-DP\n";

    std::cout
        << "LOCAL SLICE VS DIGIT-DP\n\n";

    const u64 primes[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31, 37, 41, 43
    };

    /*
        Small exhaustive prefix comparison.
    */
    Stats small_prefix;

    for (u64 p : {
        u64(2),
        u64(3),
        u64(5),
        u64(7),
        u64(11)
    }) {
        for (u128 m = 0;
             m <= 300;
             ++m) {

            for (u128 y = 0;
                 y <= 300;
                 ++y) {

                ++small_prefix.cases;

                const u128 dp =
                    miss_prefix_dp(
                        m,
                        y,
                        p
                    );

                const u128 closed =
                    miss_prefix_closed(
                        m,
                        y,
                        p
                    );

                if (dp != closed) {
                    ++small_prefix.failures;

                    if (small_prefix.printed < 5) {
                        ++small_prefix.printed;

                        print_prefix_failure(
                            p,
                            m,
                            y,
                            dp,
                            closed
                        );
                    }
                }
            }
        }
    }

    std::cout
        << "small_prefix_cases="
        << small_prefix.cases
        << " fail="
        << small_prefix.failures
        << '\n';

    /*
        Large arbitrary prefix comparison.
    */
    Stats random_prefix;

    std::mt19937_64 rng(
        0x286286ULL
    );

    for (std::size_t i = 0;
         i < 50000;
         ++i) {

        const u64 p =
            primes[
                rng() %
                (sizeof(primes) /
                 sizeof(primes[0]))
            ];

        const u128 m =
            static_cast<u128>(
                rng()
            ) %
            static_cast<u128>(
                1000000000000000000ULL
            );

        const u128 y =
            static_cast<u128>(
                rng()
            ) %
            static_cast<u128>(
                1000000000000000000ULL
            );

        ++random_prefix.cases;

        const u128 dp =
            miss_prefix_dp(
                m,
                y,
                p
            );

        const u128 closed =
            miss_prefix_closed(
                m,
                y,
                p
            );

        if (dp != closed) {
            ++random_prefix.failures;

            if (random_prefix.printed < 5) {
                ++random_prefix.printed;

                print_prefix_failure(
                    p,
                    m,
                    y,
                    dp,
                    closed
                );
            }
        }
    }

    std::cout
        << "random_prefix_cases="
        << random_prefix.cases
        << " fail="
        << random_prefix.failures
        << '\n';

    /*
        Large local-update comparison.
    */
    Stats local;

    for (std::size_t i = 0;
         i < 50000;
         ++i) {

        const u64 p =
            primes[
                rng() %
                (sizeof(primes) /
                 sizeof(primes[0]))
            ];

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

        if (!legal_update(
                current,
                r,
                p
            )) {
            continue;
        }

        const u128 step =
            power_u128(
                static_cast<u128>(p),
                r
            );

        const u128 next =
            current + step;

        const u128 n =
            static_cast<u128>(
                rng()
            ) %
            static_cast<u128>(
                1000000000000000000ULL
            );

        ++local.cases;

        const u128 exact =
            exact_local_delta_dp(
                current,
                next,
                n,
                p
            );

        const u128 closed =
            local_slice_closed(
                current,
                r,
                n,
                p
            );

        if (exact == static_cast<u128>(-1) ||
            exact != closed) {

            ++local.failures;

            if (local.printed < 5) {
                ++local.printed;

                print_local_failure(
                    p,
                    current,
                    r,
                    n,
                    exact,
                    closed
                );
            }
        }
    }

    std::cout
        << "local_cases="
        << local.cases
        << " fail="
        << local.failures
        << '\n';

    const bool pass =
        small_prefix.failures == 0 &&
        random_prefix.failures == 0 &&
        local.failures == 0;

    std::cout
        << "\nOVERALL PASS="
        << (pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 286\n";

    return pass ? 0 : 1;
}