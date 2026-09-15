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
// MISS predicate
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
// Apply m -> m + p^r
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
// Signed subtraction
//
// a - b
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
// Apply one permutation and sum the three local deltas
// ============================================================

Signed128 permutation_path(
    u64 p,
    u64 m,
    u64 n,
    const std::array<std::size_t, 3>& order
) {
    u64 current = m;

    Signed128 total{
        false,
        0
    };

    for (std::size_t k = 0;
         k < 3;
         ++k) {

        const std::size_t r =
            order[k];

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
// Verify all six permutations
// ============================================================

bool verify_triple(
    u64 p,
    u64 m,
    std::size_t r,
    std::size_t s,
    std::size_t t,
    u64 n
) {
    if (r == s ||
        r == t ||
        s == t) {

        return false;
    }

    // All three initial digit updates must be legal.
    if (!legal_digit_increment(p, m, r) ||
        !legal_digit_increment(p, m, s) ||
        !legal_digit_increment(p, m, t)) {

        return false;
    }

    // Because r,s,t are distinct and each update is
    // no-carry, every later update remains legal.
    const u64 target =
        m +
        (u64)p_power(p, r) +
        (u64)p_power(p, s) +
        (u64)p_power(p, t);

    const std::array<std::array<std::size_t, 3>, 6>
        permutations = {{
            {{r, s, t}},
            {{r, t, s}},
            {{s, r, t}},
            {{s, t, r}},
            {{t, r, s}},
            {{t, s, r}}
        }};

    const Signed128 reference =
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

    for (const auto& order : permutations) {
        const Signed128 path =
            permutation_path(
                p,
                m,
                n,
                order
            );

        if (!signed_equal(
                path,
                reference)) {

            return false;
        }
    }

    // Explicitly verify all six paths are equal as well.
    const Signed128 path0 =
        permutation_path(
            p,
            m,
            n,
            permutations[0]
        );

    for (std::size_t i = 1;
         i < permutations.size();
         ++i) {

        const Signed128 path =
            permutation_path(
                p,
                m,
                n,
                permutations[i]
            );

        if (!signed_equal(
                path,
                path0)) {

            return false;
        }
    }

    return true;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout
        << "START EXPERIMENT 271\n";

    std::cout
        << "THREE-DIGIT PATH INDEPENDENCE\n";

    std::cout
        << "ALL 6 ORDERS OF THREE NO-CARRY UPDATES\n\n";

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 0},
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

        const auto md =
            base_digits(
                p,
                m
            );

        std::vector<std::size_t> valid;

        for (std::size_t i = 0;
             i < md.size();
             ++i) {

            if (md[i] < p - 1) {
                valid.push_back(i);
            }
        }

        for (std::size_t a = 0;
             a < valid.size();
             ++a) {

            for (std::size_t b = a + 1;
                 b < valid.size();
                 ++b) {

                for (std::size_t c = b + 1;
                     c < valid.size();
                     ++c) {

                    const std::size_t r =
                        valid[a];

                    const std::size_t s =
                        valid[b];

                    const std::size_t t =
                        valid[c];

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

                        if (!verify_triple(
                                p,
                                m,
                                r,
                                s,
                                t,
                                n)) {

                            ok = false;
                        }
                    }
                }
            }
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
             m <= 500;
             ++m) {

            const auto md =
                base_digits(
                    p,
                    m
                );

            std::vector<std::size_t> valid;

            for (std::size_t i = 0;
                 i < md.size();
                 ++i) {

                if (md[i] < p - 1) {
                    valid.push_back(i);
                }
            }

            if (valid.size() < 3) {
                continue;
            }

            for (std::size_t a = 0;
                 a < valid.size();
                 ++a) {

                for (std::size_t b = a + 1;
                     b < valid.size();
                     ++b) {

                    for (std::size_t c = b + 1;
                         c < valid.size();
                         ++c) {

                        const std::size_t r =
                            valid[a];

                        const std::size_t s =
                            valid[b];

                        const std::size_t t =
                            valid[c];

                        for (u64 n = 0;
                             n <= m + 20;
                             ++n) {

                            ++exhaustive_cases;

                            if (!verify_triple(
                                    p,
                                    m,
                                    r,
                                    s,
                                    t,
                                    n)) {

                                ++exhaustive_fail;
                            }
                        }
                    }
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
    // Random large cases
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x271271271ULL
    );

    const std::vector<u64> primes = {
        3, 5, 7, 11, 13,
        17, 19, 23, 29,
        31, 37, 41, 43, 47
    };

    const std::size_t random_cases =
        100000;

    std::size_t random_fail = 0;
    std::size_t random_tested = 0;

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

        const auto md =
            base_digits(
                p,
                m
            );

        std::vector<std::size_t> valid;

        for (std::size_t i = 0;
             i < md.size();
             ++i) {

            if (md[i] < p - 1) {
                valid.push_back(i);
            }
        }

        if (valid.size() < 3) {
            continue;
        }

        // Pick three distinct valid digit positions.
        std::size_t a =
            rng() % valid.size();

        std::size_t b =
            rng() % valid.size();

        while (b == a) {
            b =
                rng() % valid.size();
        }

        std::size_t c =
            rng() % valid.size();

        while (c == a || c == b) {
            c =
                rng() % valid.size();
        }

        const std::size_t r =
            valid[a];

        const std::size_t s =
            valid[b];

        const std::size_t t =
            valid[c];

        const u128 target128 =
            (u128)m +
            p_power(p, r) +
            p_power(p, s) +
            p_power(p, t);

        // Keep target inside uint64_t.
        if (target128 >
            (u128)UINT64_MAX) {

            continue;
        }

        const u64 n =
            rng() %
            1000000000000000001ULL;

        ++random_tested;

        if (!verify_triple(
                p,
                m,
                r,
                s,
                t,
                n)) {

            ++random_fail;
        }
    }

    std::cout
        << "random_cases="
        << random_cases
        << " tested="
        << random_tested
        << " fail="
        << random_fail
        << "\n\n";

    // --------------------------------------------------------
    // Large explicit case
    // --------------------------------------------------------

    {
        const u64 p = 13;

        const u64 m =
            987654321012345678ULL;

        const auto md =
            base_digits(
                p,
                m
            );

        std::vector<std::size_t> valid;

        for (std::size_t i = 0;
             i < md.size();
             ++i) {

            if (md[i] < p - 1) {
                valid.push_back(i);
            }
        }

        const std::size_t r =
            valid[0];

        const std::size_t s =
            valid[1];

        const std::size_t t =
            valid[2];

        const u64 target =
            (u64)(
                (u128)m +
                p_power(p, r) +
                p_power(p, s) +
                p_power(p, t)
            );

        const u64 n =
            876543210123456789ULL;

        const std::array<
            std::array<std::size_t, 3>,
            6
        > permutations = {{
            {{r, s, t}},
            {{r, t, s}},
            {{s, r, t}},
            {{s, t, r}},
            {{t, r, s}},
            {{t, s, r}}
        }};

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

        std::cout
            << "LARGE_CASE\n";

        std::cout
            << "p=" << p << '\n';

        std::cout
            << "m=" << m << '\n';

        std::cout
            << "r=" << r << '\n';

        std::cout
            << "s=" << s << '\n';

        std::cout
            << "t=" << t << '\n';

        std::cout
            << "target=" << target << '\n';

        std::cout
            << "n=" << n << '\n';

        for (std::size_t i = 0;
             i < permutations.size();
             ++i) {

            const Signed128 path =
                permutation_path(
                    p,
                    m,
                    n,
                    permutations[i]
                );

            std::cout
                << "path_"
                << (i + 1)
                << "="
                << (
                    path.negative
                        ? "-"
                        : "+"
                );

            print_u128(
                path.magnitude
            );

            std::cout << '\n';
        }

        std::cout
            << "direct="
            << (
                direct.negative
                    ? "-"
                    : "+"
            );

        print_u128(
            direct.magnitude
        );

        std::cout
            << "\nlarge_case_pass="
            << (
                verify_triple(
                    p,
                    m,
                    r,
                    s,
                    t,
                    n
                )
                    ? 1
                    : 0
            )
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
        << "FINISHED EXPERIMENT 271\n";

    return overall ? 0 : 1;
}
