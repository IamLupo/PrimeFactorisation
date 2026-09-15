#include <algorithm>
#include <array>
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
// Legal no-carry digit increment
// ============================================================

bool legal_digit_increment(
    u64 p,
    u64 m,
    std::size_t r
) {
    const auto md =
        base_digits(p, m);

    const u64 mr =
        r < md.size()
            ? md[r]
            : 0;

    return mr < p - 1;
}

// ============================================================
// Apply m -> m+p^r
// ============================================================

u64 digit_increment(
    u64 p,
    u64 m,
    std::size_t r
) {
    return
        m +
        (u64)p_power(
            p,
            r
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
// One local delta
//
// Delta_r(m,n) = H_m(n) - H_{m+p^r}(n)
// ============================================================

Signed128 actual_delta(
    u64 p,
    u64 m,
    std::size_t r,
    u64 n
) {
    const u64 m2 =
        digit_increment(
            p,
            m,
            r
        );

    const u128 old_hit =
        hit_prefix(
            p,
            m,
            n
        );

    const u128 new_hit =
        hit_prefix(
            p,
            m2,
            n
        );

    return signed_difference(
        old_hit,
        new_hit
    );
}

// ============================================================
// Sum the local deltas in a specified order
//
// Returns H_m - H_target for that path.
// ============================================================

Signed128 permutation_path(
    u64 p,
    u64 m,
    u64 n,
    const std::vector<std::size_t>& order
) {
    u64 current = m;

    Signed128 total{
        false,
        0
    };

    for (std::size_t r : order) {

        const Signed128 delta =
            actual_delta(
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
// Target after all selected updates
// ============================================================

u64 target_after_updates(
    u64 p,
    u64 m,
    const std::vector<std::size_t>& digits
) {
    u128 result = (u128)m;

    for (std::size_t r : digits) {
        result += p_power(p, r);
    }

    return (u64)result;
}

// ============================================================
// Verify one update set and one order
// ============================================================

bool verify_order(
    u64 p,
    u64 m,
    u64 n,
    const std::vector<std::size_t>& digits,
    const std::vector<std::size_t>& order
) {
    const u64 target =
        target_after_updates(
            p,
            m,
            digits
        );

    const Signed128 path =
        permutation_path(
            p,
            m,
            n,
            order
        );

    const Signed128 direct =
        signed_difference(
            hit_prefix(
                p,
                m,
                n
            ),
            hit_prefix(
                p,
                target,
                n
            )
        );

    return signed_equal(
        path,
        direct
    );
}

// ============================================================
// Build legal digit positions
// ============================================================

std::vector<std::size_t> legal_digits(
    u64 p,
    u64 m
) {
    const auto md =
        base_digits(
            p,
            m
        );

    std::vector<std::size_t> result;

    for (std::size_t i = 0;
         i < md.size();
         ++i) {

        if (md[i] < p - 1) {
            result.push_back(i);
        }
    }

    return result;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout
        << "START EXPERIMENT 272\n";

    std::cout
        << "K-DIGIT PATH INDEPENDENCE\n";

    std::cout
        << "ARBITRARY NUMBER OF INDEPENDENT NO-CARRY UPDATES\n";

    std::cout
        << "EXHAUSTIVE SMALL ORDERS + RANDOM LARGE ORDERS\n\n";

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {3, 0},
        {3, 5},
        {3, 10},
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

        const auto valid =
            legal_digits(
                p,
                m
            );

        const std::size_t max_k =
            std::min<std::size_t>(
                valid.size(),
                6
            );

        for (std::size_t k = 1;
             k <= max_k;
             ++k) {

            std::vector<std::size_t>
                selected(
                    valid.begin(),
                    valid.begin() + k
                );

            std::vector<std::size_t>
                order = selected;

            do {
                const std::vector<u64> tests = {
                    0,
                    1,
                    m / 2,
                    m,
                    m + 1,
                    m + 10
                };

                for (u64 n : tests) {

                    ++deterministic_tests;

                    if (!verify_order(
                            p,
                            m,
                            n,
                            selected,
                            order)) {

                        ok = false;
                    }
                }

            } while (
                std::next_permutation(
                    order.begin(),
                    order.end()
                )
            );
        }

        if (!ok) {
            ++deterministic_fail;
        }

        std::cout
            << "p=" << p
            << " m=" << m
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
    // Exhaustive small
    //
    // Enumerate every subset of size k <= 5 and every
    // permutation for very small m.
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
             m <= 250;
             ++m) {

            const auto valid =
                legal_digits(
                    p,
                    m
                );

            const std::size_t V =
                valid.size();

            for (std::size_t k = 1;
                 k <= std::min<std::size_t>(
                        V,
                        5
                    );
                 ++k) {

                std::vector<bool>
                    choose(V, false);

                std::fill(
                    choose.begin(),
                    choose.begin() + k,
                    true
                );

                do {
                    std::vector<std::size_t>
                        selected;

                    for (std::size_t i = 0;
                         i < V;
                         ++i) {

                        if (choose[i]) {
                            selected.push_back(
                                valid[i]
                            );
                        }
                    }

                    std::vector<std::size_t>
                        order = selected;

                    do {
                        for (u64 n = 0;
                             n <= m + 20;
                             ++n) {

                            ++exhaustive_cases;

                            if (!verify_order(
                                    p,
                                    m,
                                    n,
                                    selected,
                                    order)) {

                                ++exhaustive_fail;
                            }
                        }

                    } while (
                        std::next_permutation(
                            order.begin(),
                            order.end()
                        )
                    );

                } while (
                    std::prev_permutation(
                        choose.begin(),
                        choose.end()
                    )
                );
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
    // Random large cases
    //
    // Test k from 1 up to 12.
    // For each case, test 20 random permutations.
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x272272272ULL
    );

    const std::vector<u64> primes = {
        3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31, 37, 41,
        43, 47
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

        const auto valid =
            legal_digits(
                p,
                m
            );

        if (valid.empty()) {
            continue;
        }

        const std::size_t max_k =
            std::min<std::size_t>(
                valid.size(),
                12
            );

        const std::size_t k =
            1 +
            (rng() % max_k);

        // Random distinct subset.
        std::vector<std::size_t>
            pool = valid;

        std::shuffle(
            pool.begin(),
            pool.end(),
            rng
        );

        std::vector<std::size_t>
            selected(
                pool.begin(),
                pool.begin() + k
            );

        // Ensure endpoint fits uint64.
        const u128 target128 =
            (u128)m;

        u128 target =
            target128;

        for (std::size_t r :
             selected) {

            target +=
                p_power(
                    p,
                    r
                );
        }

        if (target >
            (u128)UINT64_MAX) {

            continue;
        }

        const u64 n =
            rng() %
            1000000000000000001ULL;

        // ----------------------------------------------------
        // Canonical sorted order.
        // ----------------------------------------------------

        std::vector<std::size_t>
            sorted =
                selected;

        std::sort(
            sorted.begin(),
            sorted.end()
        );

        ++random_tested;

        if (!verify_order(
                p,
                m,
                n,
                selected,
                sorted)) {

            ++random_fail;
            continue;
        }

        // ----------------------------------------------------
        // Random permutations.
        // ----------------------------------------------------

        std::vector<std::size_t>
            order =
                selected;

        for (int path = 0;
             path < 20;
             ++path) {

            std::shuffle(
                order.begin(),
                order.end(),
                rng
            );

            ++random_paths;

            if (!verify_order(
                    p,
                    m,
                    n,
                    selected,
                    order)) {

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
    // Large explicit cases for k = 1..10
    // --------------------------------------------------------

    {
        const u64 p = 13;

        const u64 m =
            987654321012345678ULL;

        const auto valid =
            legal_digits(
                p,
                m
            );

        const u64 n =
            876543210123456789ULL;

        std::cout
            << "LARGE_CASE\n";

        std::cout
            << "p=" << p << '\n';

        std::cout
            << "m=" << m << '\n';

        std::cout
            << "n=" << n << '\n';

        bool large_ok = true;

        const std::size_t max_k =
            std::min<std::size_t>(
                valid.size(),
                10
            );

        for (std::size_t k = 1;
             k <= max_k;
             ++k) {

            std::vector<std::size_t>
                selected(
                    valid.begin(),
                    valid.begin() + k
                );

            const u128 target128 =
                (u128)m;

            u128 target =
                target128;

            for (std::size_t r :
                 selected) {

                target +=
                    p_power(
                        p,
                        r
                    );
            }

            if (target >
                (u128)UINT64_MAX) {

                continue;
            }

            const u64 target64 =
                (u64)target;

            const Signed128 direct =
                signed_difference(
                    hit_prefix(
                        p,
                        m,
                        n
                    ),
                    hit_prefix(
                        p,
                        target64,
                        n
                    )
                );

            std::vector<std::size_t>
                order =
                selected;

            const Signed128
                sorted_path =
                permutation_path(
                    p,
                    m,
                    n,
                    order
                );

            bool k_ok =
                signed_equal(
                    sorted_path,
                    direct
                );

            // Also test 100 random orders for
            // this particular large case.
            for (int trial = 0;
                 trial < 100;
                 ++trial) {

                std::shuffle(
                    order.begin(),
                    order.end(),
                    rng
                );

                const Signed128 path =
                    permutation_path(
                        p,
                        m,
                        n,
                        order
                    );

                if (!signed_equal(
                        path,
                        direct
                    )) {

                    k_ok = false;
                }
            }

            std::cout
                << "k=" << k
                << " pass="
                << (k_ok ? 1 : 0)
                << " target="
                << target64
                << '\n';

            if (!k_ok) {
                large_ok = false;
            }
        }

        std::cout
            << "large_case_pass="
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
        << "FINISHED EXPERIMENT 272\n";

    return overall ? 0 : 1;
}
