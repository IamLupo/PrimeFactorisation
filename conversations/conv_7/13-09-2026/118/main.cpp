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
// Local update
//
// Delta_r = H_m - H_{m+p^r}
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
// Verify commuting pair
// ============================================================

bool verify_pair(
    u64 p,
    u64 m,
    std::size_t r,
    std::size_t s,
    u64 n
) {
    if (r == s) {
        return false;
    }

    // Both initial updates must be legal.
    if (!legal_digit_increment(
            p,
            m,
            r) ||
        !legal_digit_increment(
            p,
            m,
            s)) {

        return false;
    }

    const u64 m_r =
        digit_increment(
            p,
            m,
            r
        );

    const u64 m_s =
        digit_increment(
            p,
            m,
            s
        );

    // Both cross-updates must remain legal.
    if (!legal_digit_increment(
            p,
            m_r,
            s) ||
        !legal_digit_increment(
            p,
            m_s,
            r)) {

        return false;
    }

    const u64 target =
        digit_increment(
            p,
            m_r,
            s
        );

    // --------------------------------------------------------
    // Path A
    // --------------------------------------------------------

    const Signed128 path_a =
        add_signed(
            actual_delta(
                p,
                m,
                r,
                n
            ),
            actual_delta(
                p,
                m_r,
                s,
                n
            )
        );

    // --------------------------------------------------------
    // Path B
    // --------------------------------------------------------

    const Signed128 path_b =
        add_signed(
            actual_delta(
                p,
                m,
                s,
                n
            ),
            actual_delta(
                p,
                m_s,
                r,
                n
            )
        );

    // --------------------------------------------------------
    // Direct endpoint difference
    // --------------------------------------------------------

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

    return
        signed_equal(path_a, path_b) &&
        signed_equal(path_a, direct);
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout
        << "START EXPERIMENT 270\n";

    std::cout
        << "DIGIT-UPDATE PATH INDEPENDENCE\n";

    std::cout
        << "COMMUTING NO-CARRY DIGIT UPDATES\n";

    std::cout
        << "VALID PAIRS ONLY\n\n";

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 0},
        {2, 1},
        {2, 2},
        {2, 4},
        {2, 8},
        {3, 5},
        {3, 10},
        {3, 20},
        {5, 12},
        {5, 124},
        {7, 100},
        {11, 12345},
        {13, 987654321},
        {17, 1000000}
    };

    std::size_t deterministic_fail = 0;
    std::size_t deterministic_tests = 0;

    for (const auto& [p, m] : cases) {

        bool ok = true;

        const auto md =
            base_digits(p, m);

        for (std::size_t r = 0;
             r < md.size();
             ++r) {

            if (!legal_digit_increment(
                    p,
                    m,
                    r)) {

                continue;
            }

            for (std::size_t s = r + 1;
                 s < md.size();
                 ++s) {

                if (!legal_digit_increment(
                        p,
                        m,
                        s)) {

                    continue;
                }

                const u64 m_r =
                    digit_increment(
                        p,
                        m,
                        r
                    );

                if (!legal_digit_increment(
                        p,
                        m_r,
                        s)) {

                    continue;
                }

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

                    if (!verify_pair(
                            p,
                            m,
                            r,
                            s,
                            n)) {

                        ok = false;
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
    // Exhaustive small verification
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

            for (std::size_t r = 0;
                 r < md.size();
                 ++r) {

                if (!legal_digit_increment(
                        p,
                        m,
                        r)) {

                    continue;
                }

                for (std::size_t s = r + 1;
                     s < md.size();
                     ++s) {

                    if (!legal_digit_increment(
                            p,
                            m,
                            s)) {

                        continue;
                    }

                    const u64 m_r =
                        digit_increment(
                            p,
                            m,
                            r
                        );

                    if (!legal_digit_increment(
                            p,
                            m_r,
                            s)) {

                        continue;
                    }

                    for (u64 n = 0;
                         n <= m + 20;
                         ++n) {

                        ++exhaustive_cases;

                        if (!verify_pair(
                                p,
                                m,
                                r,
                                s,
                                n)) {

                            ++exhaustive_fail;
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
    // Random
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x270270270ULL
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

        if (valid.size() < 2) {
            continue;
        }

        const std::size_t index_r =
            rng() % valid.size();

        std::size_t index_s =
            rng() % valid.size();

        if (index_r == index_s) {
            index_s =
                (index_s + 1) %
                valid.size();
        }

        std::size_t r =
            valid[index_r];

        std::size_t s =
            valid[index_s];

        // Normalize mutable indices.
        if (r > s) {
            const std::size_t tmp = r;
            r = s;
            s = tmp;
        }

        const u64 m_r =
            digit_increment(
                p,
                m,
                r
            );

        if (!legal_digit_increment(
                p,
                m_r,
                s)) {

            continue;
        }

        const u64 n =
            rng() %
            1000000000000000001ULL;

        ++random_tested;

        if (!verify_pair(
                p,
                m,
                r,
                s,
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
    // Large case
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

        const u64 m_r =
            digit_increment(
                p,
                m,
                r
            );

        const u64 m_s =
            digit_increment(
                p,
                m,
                s
            );

        const u64 target =
            digit_increment(
                p,
                m_r,
                s
            );

        const u64 n =
            876543210123456789ULL;

        const Signed128 path_a =
            add_signed(
                actual_delta(
                    p,
                    m,
                    r,
                    n
                ),
                actual_delta(
                    p,
                    m_r,
                    s,
                    n
                )
            );

        const Signed128 path_b =
            add_signed(
                actual_delta(
                    p,
                    m,
                    s,
                    n
                ),
                actual_delta(
                    p,
                    m_s,
                    r,
                    n
                )
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
            << "target=" << target << '\n';

        std::cout
            << "n=" << n << '\n';

        std::cout
            << "path_a="
            << (path_a.negative ? "-" : "+");

        print_u128(
            path_a.magnitude
        );

        std::cout
            << "\npath_b="
            << (path_b.negative ? "-" : "+");

        print_u128(
            path_b.magnitude
        );

        std::cout
            << "\ndirect="
            << (direct.negative ? "-" : "+");

        print_u128(
            direct.magnitude
        );

        std::cout
            << "\nlarge_case_pass="
            << (
                signed_equal(
                    path_a,
                    path_b
                ) &&
                signed_equal(
                    path_a,
                    direct
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
        << "FINISHED EXPERIMENT 270\n";

    return overall ? 0 : 1;
}
