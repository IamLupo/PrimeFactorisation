#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

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
// Base-p digits, least-significant first
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
// Digitwise MISS predicate
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
// HIT prefix
//
// Number of HIT x with 0 <= x < n.
// ============================================================

u128 hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    return
        (u128)n -
        miss_prefix(
            p,
            m,
            n - 1
        );
}

// ============================================================
// Signed 128
// ============================================================

struct Signed128 {
    bool negative;
    u128 magnitude;
};

// ============================================================
// Signed difference
//
// a-b
// ============================================================

Signed128 signed_difference(
    u128 a,
    u128 b
) {
    if (a >= b) {
        return {
            false,
            a - b
        };
    }

    return {
        true,
        b - a
    };
}

// ============================================================
// Signed addition
// ============================================================

Signed128 add_signed(
    Signed128 a,
    Signed128 b
) {
    if (a.negative == b.negative) {
        return {
            a.negative,
            a.magnitude + b.magnitude
        };
    }

    if (a.magnitude >= b.magnitude) {
        return {
            a.negative,
            a.magnitude - b.magnitude
        };
    }

    return {
        b.negative,
        b.magnitude - a.magnitude
    };
}

// ============================================================
// Signed equality
// ============================================================

bool signed_equal(
    Signed128 a,
    Signed128 b
) {
    return
        a.negative == b.negative &&
        a.magnitude == b.magnitude;
}

// ============================================================
// One legal local digit increment
//
// current -> current + p^r
// ============================================================

u64 digit_increment(
    u64 p,
    u64 current,
    std::size_t r
) {
    return current + (u64)p_power(p, r);
}

// ============================================================
// Local potential increment
//
// Delta_r(current,n)
//     = H_current(n) - H_{current+p^r}(n)
// ============================================================

Signed128 local_delta(
    u64 p,
    u64 current,
    std::size_t r,
    u64 n
) {
    const u64 next =
        digit_increment(
            p,
            current,
            r
        );

    return signed_difference(
        hit_prefix(
            p,
            current,
            n
        ),
        hit_prefix(
            p,
            next,
            n
        )
    );
}

// ============================================================
// Construct target from base-p digits.
//
// Returns m itself; the function makes the intended
// digit-vector construction explicit for the experiment.
// ============================================================

u64 reconstruct_from_digits(
    u64 p,
    const std::vector<u64>& digits
) {
    u128 result = 0;
    u128 power = 1;

    for (u64 d : digits) {
        result +=
            (u128)d * power;

        power *=
            (u128)p;
    }

    return (u64)result;
}

// ============================================================
// Build the multiset of elementary updates.
//
// Digit r appears exactly m_r times.
// ============================================================

std::vector<std::size_t> update_multiset(
    const std::vector<u64>& digits
) {
    std::vector<std::size_t> updates;

    for (std::size_t r = 0;
         r < digits.size();
         ++r) {

        for (u64 k = 0;
             k < digits[r];
             ++k) {

            updates.push_back(r);
        }
    }

    return updates;
}

// ============================================================
// Follow one complete update path.
//
// Starts at m=0 and applies every elementary update.
//
// Returns accumulated
//
//   H_0(n) - H_m(n).
// ============================================================

Signed128 path_sum(
    u64 p,
    u64 n,
    const std::vector<std::size_t>& order
) {
    u64 current = 0;

    Signed128 total{
        false,
        0
    };

    for (std::size_t r : order) {

        const Signed128 delta =
            local_delta(
                p,
                current,
                r,
                n
            );

        total =
            add_signed(
                total,
                delta
            );

        current =
            digit_increment(
                p,
                current,
                r
            );
    }

    return total;
}

// ============================================================
// Closed-form endpoint difference
//
// H_0(n) - H_m(n).
//
// H_0(n) has MISS set {0}, so:
//
//   H_0(n) = n-1 for n>0,
//   H_0(0) = 0.
// ============================================================

Signed128 direct_endpoint_difference(
    u64 p,
    u64 m,
    u64 n
) {
    const u128 H0 =
        n == 0
            ? (u128)0
            : (u128)n - 1;

    const u128 Hm =
        hit_prefix(
            p,
            m,
            n
        );

    return signed_difference(
        H0,
        Hm
    );
}

// ============================================================
// Closed digit formula for H_m(n)
//
// This is the formula validated previously.
//
// If n <=_p m:
//
//   H_m(n) = n - sum n_i W_i
//
// Otherwise, h is the most-significant digit with n_h > m_h:
//
//   H_m(n)
//     = n - [
//         sum_{i>h} n_i W_i
//         + (m_h+1)W_h
//       ].
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
        base_digits(
            p,
            n
        );

    const auto md =
        base_digits(
            p,
            m
        );

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
// Verify one complete endpoint
//
// Checks:
//
// 1. digit reconstruction m
// 2. direct H_m
// 3. closed H_m
// 4. H_0-H_m endpoint difference
// ============================================================

bool verify_endpoint(
    u64 p,
    u64 m,
    u64 n
) {
    const auto md =
        base_digits(
            p,
            m
        );

    const u64 reconstructed =
        reconstruct_from_digits(
            p,
            md
        );

    if (reconstructed != m) {
        return false;
    }

    const u128 direct =
        hit_prefix(
            p,
            m,
            n
        );

    const u128 closed =
        closed_hit_prefix(
            p,
            m,
            n
        );

    if (direct != closed) {
        return false;
    }

    const Signed128 expected =
        direct_endpoint_difference(
            p,
            m,
            n
        );

    const Signed128 zero_path =
        signed_difference(
            n == 0
                ? (u128)0
                : (u128)n - 1,
            direct
        );

    return signed_equal(
        expected,
        zero_path
    );
}

// ============================================================
// Enumerate all unique multiset permutations.
//
// Only used when total update count is small.
// ============================================================

bool exhaustive_paths_recursive(
    u64 p,
    u64 n,
    const std::vector<std::size_t>& updates,
    std::vector<bool>& used,
    std::vector<std::size_t>& path,
    Signed128 expected,
    std::size_t& tested,
    std::size_t& failures
) {
    if (path.size() == updates.size()) {

        ++tested;

        const Signed128 actual =
            path_sum(
                p,
                n,
                path
            );

        if (!signed_equal(
                actual,
                expected
            )) {

            ++failures;
            return false;
        }

        return true;
    }

    bool ok = true;

    std::size_t previous =
        static_cast<std::size_t>(-1);

    for (std::size_t i = 0;
         i < updates.size();
         ++i) {

        if (used[i]) {
            continue;
        }

        const std::size_t digit =
            updates[i];

        // Skip identical multiset choices at this level.
        if (previous !=
            static_cast<std::size_t>(-1) &&
            previous == digit) {

            continue;
        }

        previous = digit;

        used[i] = true;
        path.push_back(digit);

        if (!exhaustive_paths_recursive(
                p,
                n,
                updates,
                used,
                path,
                expected,
                tested,
                failures)) {

            ok = false;
        }

        path.pop_back();
        used[i] = false;
    }

    return ok;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout
        << "START EXPERIMENT 273\n";

    std::cout
        << "GLOBAL POTENTIAL RECONSTRUCTION FROM m=0\n";

    std::cout
        << "LOCAL DIGIT UPDATES -> H_0 - H_m\n";

    std::cout
        << "EXHAUSTIVE SMALL PATHS + RANDOM LARGE PATHS\n\n";

    // --------------------------------------------------------
    // Deterministic cases
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 0},
        {2, 1},
        {2, 2},
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
        {31, 987654321012ULL}
    };

    std::size_t deterministic_fail = 0;
    std::size_t deterministic_tests = 0;

    for (const auto& [p, m] : cases) {

        bool ok = true;

        const auto digits =
            base_digits(
                p,
                m
            );

        const auto updates =
            update_multiset(
                digits
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

            ++deterministic_tests;

            if (!verify_endpoint(
                    p,
                    m,
                    n)) {

                ok = false;
            }
        }

        // Full path enumeration only when small.
        if (updates.size() <= 8) {

            for (u64 n : tests) {

                const Signed128 expected =
                    direct_endpoint_difference(
                        p,
                        m,
                        n
                    );

                std::vector<bool> used(
                    updates.size(),
                    false
                );

                std::vector<std::size_t> path;

                std::size_t path_tests = 0;
                std::size_t path_fail = 0;

                exhaustive_paths_recursive(
                    p,
                    n,
                    updates,
                    used,
                    path,
                    expected,
                    path_tests,
                    path_fail
                );

                if (path_fail != 0) {
                    ok = false;
                }
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
        << " tests="
        << deterministic_tests
        << " fail="
        << deterministic_fail
        << "\n\n";

    // --------------------------------------------------------
    // Exhaustive small m
    //
    // Enumerate all unique paths whenever the total number
    // of elementary updates is <= 7.
    // --------------------------------------------------------

    std::size_t exhaustive_endpoints = 0;
    std::size_t exhaustive_paths = 0;
    std::size_t exhaustive_fail = 0;

    for (u64 p : {
        2ULL,
        3ULL,
        5ULL,
        7ULL
    }) {

        for (u64 m = 0;
             m <= 250;
             ++m) {

            const auto digits =
                base_digits(
                    p,
                    m
                );

            const auto updates =
                update_multiset(
                    digits
                );

            if (updates.size() > 7) {
                continue;
            }

            for (u64 n = 0;
                 n <= m + 10;
                 ++n) {

                ++exhaustive_endpoints;

                if (!verify_endpoint(
                        p,
                        m,
                        n)) {

                    ++exhaustive_fail;
                    continue;
                }

                const Signed128 expected =
                    direct_endpoint_difference(
                        p,
                        m,
                        n
                    );

                std::vector<bool> used(
                    updates.size(),
                    false
                );

                std::vector<std::size_t> path;

                exhaustive_paths_recursive(
                    p,
                    n,
                    updates,
                    used,
                    path,
                    expected,
                    exhaustive_paths,
                    exhaustive_fail
                );
            }
        }
    }

    std::cout
        << "exhaustive_endpoints="
        << exhaustive_endpoints
        << " paths="
        << exhaustive_paths
        << " fail="
        << exhaustive_fail
        << "\n\n";

    // --------------------------------------------------------
    // Random large cases
    //
    // Randomize m, n, and many random paths.
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x273273273ULL
    );

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        41, 43, 47
    };

    const std::size_t random_cases =
        100000;

    std::size_t random_fail = 0;
    std::size_t random_tested = 0;
    std::size_t random_paths = 0;

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

        const auto digits =
            base_digits(
                p,
                m
            );

        const auto updates =
            update_multiset(
                digits
            );

        ++random_tested;

        if (!verify_endpoint(
                p,
                m,
                n)) {

            ++random_fail;
            continue;
        }

        const Signed128 expected =
            direct_endpoint_difference(
                p,
                m,
                n
            );

        if (updates.empty()) {
            continue;
        }

        std::vector<std::size_t> order =
            updates;

        const std::size_t trials =
            std::min<std::size_t>(
                50,
                std::max<std::size_t>(
                    5,
                    updates.size() * 2
                )
            );

        for (std::size_t trial = 0;
             trial < trials;
             ++trial) {

            std::shuffle(
                order.begin(),
                order.end(),
                rng
            );

            ++random_paths;

            const Signed128 actual =
                path_sum(
                    p,
                    n,
                    order
                );

            if (!signed_equal(
                    actual,
                    expected
                )) {

                ++random_fail;
                break;
            }
        }
    }

    std::cout
        << "random_cases="
        << random_cases
        << " tested="
        << random_tested
        << " paths="
        << random_paths
        << " fail="
        << random_fail
        << "\n\n";

    // --------------------------------------------------------
    // Large explicit cases
    // --------------------------------------------------------

    {
        const u64 p = 13;

        const u64 m =
            987654321012345678ULL;

        const u64 n =
            876543210123456789ULL;

        const auto digits =
            base_digits(
                p,
                m
            );

        const auto updates =
            update_multiset(
                digits
            );

        const Signed128 expected =
            direct_endpoint_difference(
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
            << "update_count="
            << updates.size()
            << '\n';

        bool large_ok = true;

        // Sorted/canonical path.
        {
            const Signed128 actual =
                path_sum(
                    p,
                    n,
                    updates
                );

            const bool pass =
                signed_equal(
                    actual,
                    expected
                );

            std::cout
                << "canonical_path_pass="
                << (pass ? 1 : 0)
                << '\n';

            if (!pass) {
                large_ok = false;
            }
        }

        // 500 random paths.
        std::vector<std::size_t> order =
            updates;

        for (int trial = 0;
             trial < 500;
             ++trial) {

            std::shuffle(
                order.begin(),
                order.end(),
                rng
            );

            const Signed128 actual =
                path_sum(
                    p,
                    n,
                    order
                );

            if (!signed_equal(
                    actual,
                    expected
                )) {

                large_ok = false;
                break;
            }
        }

        std::cout
            << "random_500_paths_pass="
            << (large_ok ? 1 : 0)
            << '\n';

        std::cout
            << "direct_H0_minus_Hm=";

        print_u128(
            expected.magnitude
        );

        std::cout << '\n';

        // Explicit direct endpoint values.
        const u128 H0 =
            n == 0
                ? (u128)0
                : (u128)n - 1;

        const u128 Hm =
            hit_prefix(
                p,
                m,
                n
            );

        std::cout
            << "H0=";

        print_u128(H0);

        std::cout
            << "\nHm=";

        print_u128(Hm);

        std::cout
            << "\nlarge_case_pass="
            << (large_ok ? 1 : 0)
            << "\n\n";
    }

    // --------------------------------------------------------
    // Overall
    // --------------------------------------------------------

    const bool overall =
        deterministic_fail == 0 &&
        exhaustive_fail == 0 &&
        random_fail == 0;

    std::cout
        << "OVERALL PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 273\n";

    return overall ? 0 : 1;
}