#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using Clock = std::chrono::steady_clock;

// ============================================================
// Printing
// ============================================================

void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    std::string s;

    while (x > 0) {
        s.push_back(char('0' + (x % 10)));
        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    std::cout << s;
}

// ============================================================
// Base-p digits
// ============================================================

std::vector<u64> base_digits(
    u64 p,
    u64 n
) {
    std::vector<u64> d;

    if (n == 0) {
        d.push_back(0);
        return d;
    }

    while (n > 0) {
        d.push_back(n % p);
        n /= p;
    }

    return d;
}

// ============================================================
// p^r
// ============================================================

u128 p_power(
    u64 p,
    std::size_t r
) {
    u128 result = 1;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        result *= (u128)p;
    }

    return result;
}

// ============================================================
// Mixed-radix weights
//
// Missing m digits are zero, hence radix = 1.
// ============================================================

std::vector<u128> mixed_weights(
    const std::vector<u64>& md,
    std::size_t L
) {
    std::vector<u128> W(L + 1);

    W[0] = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        const u64 mi =
            i < md.size()
                ? md[i]
                : 0;

        W[i + 1] =
            W[i] * (u128)(mi + 1);
    }

    return W;
}

// ============================================================
// Digitwise MISS
// ============================================================

bool is_miss(
    u64 p,
    u64 m,
    u64 x
) {
    const auto xd =
        base_digits(p, x);

    const auto md =
        base_digits(p, m);

    const std::size_t L =
        std::max(
            xd.size(),
            md.size()
        );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        const u64 xi =
            i < xd.size()
                ? xd[i]
                : 0;

        const u64 mi =
            i < md.size()
                ? md[i]
                : 0;

        if (xi > mi) {
            return false;
        }
    }

    return true;
}

// ============================================================
// MISS prefix
//
// Number of MISS x with 0 <= x <= y.
// ============================================================

u128 miss_prefix(
    u64 p,
    u64 m,
    u64 y
) {
    const auto yd =
        base_digits(p, y);

    const auto md =
        base_digits(p, m);

    const std::size_t L =
        std::max(
            yd.size(),
            md.size()
        );

    const auto W =
        mixed_weights(
            md,
            L
        );

    u128 result = 0;

    for (std::size_t pos = L;
         pos-- > 0;) {

        const u64 yi =
            pos < yd.size()
                ? yd[pos]
                : 0;

        const u64 mi =
            pos < md.size()
                ? md[pos]
                : 0;

        const u64 choices =
            std::min<u64>(
                yi,
                mi + 1
            );

        result +=
            (u128)choices * W[pos];

        if (yi > mi) {
            return result;
        }
    }

    return result + 1;
}

// ============================================================
// Closed HIT prefix
// ============================================================

u128 closed_hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    const auto nd =
        base_digits(p, n);

    const auto md =
        base_digits(p, m);

    const std::size_t L =
        std::max(
            nd.size(),
            md.size()
        );

    const auto W =
        mixed_weights(
            md,
            L
        );

    bool offending = false;
    std::size_t h = 0;

    for (std::size_t pos = L;
         pos-- > 0;) {

        const u64 ni =
            pos < nd.size()
                ? nd[pos]
                : 0;

        const u64 mi =
            pos < md.size()
                ? md[pos]
                : 0;

        if (ni > mi) {
            offending = true;
            h = pos;
            break;
        }
    }

    u128 miss = 0;

    if (!offending) {

        for (std::size_t i = 0;
             i < L;
             ++i) {

            const u64 ni =
                i < nd.size()
                    ? nd[i]
                    : 0;

            miss +=
                (u128)ni * W[i];
        }

    } else {

        for (std::size_t i = h + 1;
             i < L;
             ++i) {

            const u64 ni =
                i < nd.size()
                    ? nd[i]
                    : 0;

            miss +=
                (u128)ni * W[i];
        }

        const u64 mh =
            h < md.size()
                ? md[h]
                : 0;

        miss +=
            (u128)(mh + 1) * W[h];
    }

    return
        (u128)n - miss;
}

// ============================================================
// Correct local new-MISS slice
//
// current -> current + p^r
//
// IMPORTANT:
// If r >= md.size(), then current_r = 0.
// This is a legal update and can create a new digit.
// ============================================================

u128 new_miss_slice(
    u64 p,
    u64 current,
    std::size_t r,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    const auto md =
        base_digits(
            p,
            current
        );

    const u64 mr =
        r < md.size()
            ? md[r]
            : 0;

    if (mr >= p - 1) {
        return 0;
    }

    const u128 pr =
        p_power(
            p,
            r
        );

    const u128 pr1 =
        pr * (u128)p;

    const u128 y =
        (u128)n - 1;

    const u128 higher_y =
        y / pr1;

    const u128 digit_y =
        (y / pr) % p;

    const u64 higher_m =
        (u64)(
            (u128)current / pr1
        );

    const u64 new_digit =
        mr + 1;

    // Number of admissible lower-digit combinations.
    u128 lower_total = 1;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        const u64 mi =
            i < md.size()
                ? md[i]
                : 0;

        lower_total *=
            (u128)(mi + 1);
    }

    // --------------------------------------------------------
    // Higher prefixes strictly below higher_y.
    // --------------------------------------------------------

    u128 result = 0;

    if (higher_y > 0) {

        const u128 count =
            miss_prefix(
                p,
                higher_m,
                (u64)(higher_y - 1)
            );

        result +=
            count * lower_total;
    }

    // --------------------------------------------------------
    // Equal higher prefix must be admissible.
    // --------------------------------------------------------

    if (!is_miss(
            p,
            higher_m,
            (u64)higher_y
        )) {

        return result;
    }

    // --------------------------------------------------------
    // r digit below new digit.
    // --------------------------------------------------------

    if (digit_y < new_digit) {
        return result;
    }

    // --------------------------------------------------------
    // r digit above new digit.
    // --------------------------------------------------------

    if (digit_y > new_digit) {
        return result + lower_total;
    }

    // --------------------------------------------------------
    // Equal r digit:
    // count admissible lower values <= lower_y.
    // --------------------------------------------------------

    const u128 lower_y =
        y % pr;

    if (lower_y == 0) {
        return result + 1;
    }

    const u128 lower_m =
        (u128)current % pr;

    const u128 partial =
        miss_prefix(
            p,
            (u64)lower_m,
            (u64)lower_y
        );

    return result + partial;
}

// ============================================================
// Update multiset
// ============================================================

std::vector<std::size_t> update_multiset(
    u64 p,
    u64 m
) {
    const auto d =
        base_digits(
            p,
            m
        );

    std::vector<std::size_t> updates;

    for (std::size_t r = 0;
         r < d.size();
         ++r) {

        for (u64 k = 0;
             k < d[r];
             ++k) {

            updates.push_back(r);
        }
    }

    return updates;
}

// ============================================================
// Fast reconstruction
//
// H_0(n) - H_m(n)
// ============================================================

u128 fast_path_gain(
    u64 p,
    u64 n,
    const std::vector<std::size_t>& order
) {
    u64 current = 0;
    u128 total = 0;

    for (std::size_t r : order) {

        total +=
            new_miss_slice(
                p,
                current,
                r,
                n
            );

        current +=
            (u64)p_power(
                p,
                r
            );
    }

    return total;
}

// ============================================================
// Direct endpoint
// ============================================================

u128 direct_endpoint_gain(
    u64 p,
    u64 m,
    u64 n
) {
    const u128 H0 =
        n == 0
            ? (u128)0
            : (u128)n - 1;

    const u128 Hm =
        closed_hit_prefix(
            p,
            m,
            n
        );

    return H0 - Hm;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout
        << "START EXPERIMENT 277\n";

    std::cout
        << "FAST GLOBAL POTENTIAL RECONSTRUCTION\n";

    std::cout
        << "NEW HIGHER ZERO DIGIT UPDATES INCLUDED\n";

    std::cout
        << "m=0 -> m THROUGH ELEMENTARY DIGIT UPDATES\n\n";

    const auto start =
        Clock::now();

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 0},
        {2, 1},
        {2, 3},
        {2, 7},
        {3, 5},
        {3, 10},
        {3, 26},
        {5, 12},
        {5, 124},
        {7, 100},
        {11, 12345},
        {13, 987654321},
        {17, 1000000},
        {19, 123456789},
        {31, 987654321012ULL},
        {13, 987654321012345678ULL}
    };

    std::size_t deterministic_fail = 0;

    for (const auto& [p, m] : cases) {

        bool ok = true;

        const auto updates =
            update_multiset(
                p,
                m
            );

        const std::vector<u64> tests = {
            0,
            1,
            2,
            m / 3,
            m / 2,
            m,
            m + 1
        };

        for (u64 n : tests) {

            const u128 path =
                fast_path_gain(
                    p,
                    n,
                    updates
                );

            const u128 direct =
                direct_endpoint_gain(
                    p,
                    m,
                    n
                );

            if (path != direct) {
                ok = false;
            }
        }

        if (!ok) {
            ++deterministic_fail;
        }

        std::cout
            << "p=" << p
            << " m=" << m
            << " updates="
            << updates.size()
            << " pass="
            << (ok ? 1 : 0)
            << '\n';
    }

    std::cout
        << "deterministic_cases="
        << cases.size()
        << " fail="
        << deterministic_fail
        << "\n\n";

    // --------------------------------------------------------
    // Exhaustive small
    // --------------------------------------------------------

    std::size_t exhaustive_cases = 0;
    std::size_t exhaustive_fail = 0;

    for (u64 p : {
        2ULL,
        3ULL,
        5ULL,
        7ULL
    }) {

        for (u64 m = 0;
             m <= 1000;
             ++m) {

            const auto updates =
                update_multiset(
                    p,
                    m
                );

            for (u64 n = 0;
                 n <= m + 2;
                 ++n) {

                ++exhaustive_cases;

                const u128 path =
                    fast_path_gain(
                        p,
                        n,
                        updates
                    );

                const u128 direct =
                    direct_endpoint_gain(
                        p,
                        m,
                        n
                    );

                if (path != direct) {
                    ++exhaustive_fail;
                }
            }
        }
    }

    std::cout
        << "exhaustive_cases="
        << exhaustive_cases
        << " fail="
        << exhaustive_fail
        << "\n\n";

    // --------------------------------------------------------
    // Local slice verification
    //
    // Includes r >= digits(current.size()) deliberately.
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x277277277ULL
    );

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        41, 43, 47
    };

    const std::size_t local_cases =
        30000;

    std::size_t local_fail = 0;

    for (std::size_t tc = 0;
         tc < local_cases;
         ++tc) {

        const u64 p =
            primes[
                rng() % primes.size()
            ];

        const u64 current =
            rng() %
            1000000000000000000ULL;

        const auto md =
            base_digits(
                p,
                current
            );

        // Deliberately allow one position beyond
        // the current digit vector.
        const std::size_t max_r =
            md.size();

        const std::size_t r =
            rng() %
            (max_r + 1);

        const u64 mr =
            r < md.size()
                ? md[r]
                : 0;

        if (mr >= p - 1) {
            continue;
        }

        const u64 n =
            rng() %
            1000000000000000001ULL;

        const u64 next =
            current +
            (u64)p_power(
                p,
                r
            );

        const u128 old_hit =
            closed_hit_prefix(
                p,
                current,
                n
            );

        const u128 new_hit =
            closed_hit_prefix(
                p,
                next,
                n
            );

        const u128 actual =
            old_hit - new_hit;

        const u128 predicted =
            new_miss_slice(
                p,
                current,
                r,
                n
            );

        if (actual != predicted) {
            ++local_fail;
        }
    }

    std::cout
        << "local_slice_cases="
        << local_cases
        << " fail="
        << local_fail
        << "\n\n";

    // --------------------------------------------------------
    // Random global cases
    // --------------------------------------------------------

    const std::size_t random_cases =
        20000;

    std::size_t random_fail = 0;

    for (std::size_t tc = 0;
         tc < random_cases;
         ++tc) {

        const u64 p =
            primes[
                rng() % primes.size()
            ];

        const u64 m =
            rng() %
            1000000000000000000ULL;

        const u64 n =
            rng() %
            1000000000000000001ULL;

        const auto updates =
            update_multiset(
                p,
                m
            );

        const u128 path =
            fast_path_gain(
                p,
                n,
                updates
            );

        const u128 direct =
            direct_endpoint_gain(
                p,
                m,
                n
            );

        if (path != direct) {
            ++random_fail;
        }
    }

    std::cout
        << "random_cases="
        << random_cases
        << " fail="
        << random_fail
        << "\n\n";

    // --------------------------------------------------------
    // Explicit tiny carry-chain examples
    //
    // These specifically test creation of new digits.
    // --------------------------------------------------------

    {
        bool tiny_ok = true;

        const std::vector<
            std::pair<u64, u64>
        > tiny = {
            {2, 3},
            {2, 7},
            {2, 15},
            {3, 8},
            {3, 26},
            {5, 24},
            {5, 124}
        };

        std::cout
            << "NEW_DIGIT_CASES\n";

        for (const auto& [p, m] : tiny) {

            for (u64 n = 0;
                 n <= m + 2;
                 ++n) {

                const auto updates =
                    update_multiset(
                        p,
                        m
                    );

                const u128 path =
                    fast_path_gain(
                        p,
                        n,
                        updates
                    );

                const u128 direct =
                    direct_endpoint_gain(
                        p,
                        m,
                        n
                    );

                if (path != direct) {
                    tiny_ok = false;
                }
            }

            std::cout
                << "p=" << p
                << " m=" << m
                << " pass="
                << (
                    tiny_ok
                        ? 1
                        : 0
                )
                << '\n';
        }

        std::cout
            << "new_digit_cases_pass="
            << (tiny_ok ? 1 : 0)
            << "\n\n";
    }

    // --------------------------------------------------------
    // Large case
    // --------------------------------------------------------

    {
        const u64 p = 13;

        const u64 m =
            987654321012345678ULL;

        const u64 n =
            876543210123456789ULL;

        const auto updates =
            update_multiset(
                p,
                m
            );

        const u128 path =
            fast_path_gain(
                p,
                n,
                updates
            );

        const u128 direct =
            direct_endpoint_gain(
                p,
                m,
                n
            );

        const u128 H0 =
            (u128)n - 1;

        const u128 Hm =
            closed_hit_prefix(
                p,
                m,
                n
            );

        std::cout
            << "LARGE_CASE\n";

        std::cout
            << "p=" << p << '\n';

        std::cout
            << "m=" << m << '\n';

        std::cout
            << "n=" << n << '\n';

        std::cout
            << "updates="
            << updates.size()
            << '\n';

        std::cout
            << "H0=";

        print_u128(
            H0
        );

        std::cout
            << "\nHm=";

        print_u128(
            Hm
        );

        std::cout
            << "\npath_gain=";

        print_u128(
            path
        );

        std::cout
            << "\ndirect_gain=";

        print_u128(
            direct
        );

        std::cout
            << "\nlarge_case_pass="
            << (path == direct ? 1 : 0)
            << '\n';
    }

    // --------------------------------------------------------
    // Runtime
    // --------------------------------------------------------

    const auto end =
        Clock::now();

    const double seconds =
        std::chrono::duration<double>(
            end - start
        ).count();

    std::cout
        << "\nruntime_seconds="
        << seconds
        << '\n';

    // --------------------------------------------------------
    // Overall
    // --------------------------------------------------------

    const bool overall =
        deterministic_fail == 0 &&
        exhaustive_fail == 0 &&
        local_fail == 0 &&
        random_fail == 0;

    std::cout
        << "OVERALL PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 277\n";

    return overall ? 0 : 1;
}