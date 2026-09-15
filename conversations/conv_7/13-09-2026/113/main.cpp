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
// Mixed-radix weights, padded to arbitrary length
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
            i < md.size() ? md[i] : 0;

        W[i + 1] =
            W[i] * (u128)(mi + 1);
    }

    return W;
}

// ============================================================
// Digitwise MISS predicate
//
// x <=_p m.
// ============================================================

bool is_miss(
    u64 p,
    u64 m,
    u64 x
) {
    const auto xd = base_digits(p, x);
    const auto md = base_digits(p, m);

    const std::size_t L =
        std::max(xd.size(), md.size());

    for (std::size_t i = 0;
         i < L;
         ++i) {

        const u64 xi =
            i < xd.size() ? xd[i] : 0;

        const u64 mi =
            i < md.size() ? md[i] : 0;

        if (xi > mi) {
            return false;
        }
    }

    return true;
}

// ============================================================
// Number of MISS values x with 0 <= x <= y.
//
// Supports y having more digits than m.
// ============================================================

u128 direct_miss_prefix(
    u64 p,
    u64 m,
    u64 y
) {
    const auto yd = base_digits(p, y);
    const auto md = base_digits(p, m);

    const std::size_t L =
        std::max(yd.size(), md.size());

    const auto W =
        mixed_weights(md, L);

    u128 result = 0;

    for (std::size_t pos = L;
         pos-- > 0;) {

        const u64 yi =
            pos < yd.size() ? yd[pos] : 0;

        const u64 mi =
            pos < md.size() ? md[pos] : 0;

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
// Number of HIT values x with 0 <= x < n.
// ============================================================

u128 direct_hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    return
        (u128)n -
        direct_miss_prefix(
            p,
            m,
            n - 1
        );
}

// ============================================================
// Find carry position h for m -> m+1.
//
// Digits 0..h-1 are p-1.
// Digit h is < p-1.
//
// If all digits are p-1, h = digits.size().
// ============================================================

std::size_t carry_position(
    u64 p,
    u64 m
) {
    const auto md =
        base_digits(p, m);

    for (std::size_t i = 0;
         i < md.size();
         ++i) {

        if (md[i] < p - 1) {
            return i;
        }
    }

    return md.size();
}

// ============================================================
// Predicted NEW MISS points gained when m -> m+1.
//
// Let h be the first non-(p-1) digit.
//
// New MISS points have:
//
//   x_i = 0                 for i < h
//   x_h = m_h + 1
//   x_i <= m_i              for i > h
//
// Count those with x < n.
// ============================================================

u128 predicted_gain(
    u64 p,
    u64 m,
    u64 n
) {
    const auto md =
        base_digits(p, m);

    const std::size_t h =
        carry_position(p, m);

    // m = p^L - 1.
    if (h == md.size()) {
        // m+1 is the first number outside the old
        // base-p digit length. It is a single new MISS point.
        //
        // For x < n this contributes exactly when
        // m+1 < n.
        const u128 mp1 =
            (u128)m + 1;

        return mp1 < (u128)n ? 1 : 0;
    }

    const u128 ph =
        p_power(p, h);

    const u128 ph1 =
        ph * (u128)p;

    const u128 T =
        (u128)n - 1;

    // New digit at h.
    const u64 new_digit =
        md[h] + 1;

    // Higher part of m.
    const u64 higher_m =
        (u64)((u128)m / ph1);

    // Higher part of T.
    const u128 higher_n =
        T / ph1;

    const u128 digit_n =
        (T / ph) % p;

    // Number of admissible higher prefixes strictly below
    // higher_n.
    u128 result = 0;

    if (higher_n > 0) {
        result +=
            direct_miss_prefix(
                p,
                higher_m,
                (u64)(higher_n - 1)
            );
    }

    // Equal higher prefix.
    if (digit_n < new_digit) {
        return result;
    }

    // If the higher prefix is not admissible, nothing more.
    if (!is_miss(
            p,
            higher_m,
            (u64)higher_n)) {

        return result;
    }

    // digit_n == or > new_digit gives one complete
    // lower-zero slice because all lower digits must be zero.
    return result + 1;
}

// ============================================================
// Predicted OLD MISS points LOST when m -> m+1.
//
// Old MISS points that disappear satisfy:
//
//   x_i arbitrary          for i < h
//   x_h <= m_h
//   x_i <= m_i             for i > h
//
// Since old lower digits had all p choices, the new system
// only retains the lower-zero point.
//
// Therefore:
//
//   loss = (p^h - 1) * A
//
// where A counts admissible prefixes from digit h upward
// below n.
// ============================================================

u128 predicted_loss(
    u64 p,
    u64 m,
    u64 n
) {
    const auto md =
        base_digits(p, m);

    const std::size_t h =
        carry_position(p, m);

    // No digit can absorb the carry:
    // m = p^L - 1.
    //
    // Nothing is lost; all old MISS values remain MISS.
    if (h == md.size()) {
        return 0;
    }

    const u128 ph =
        p_power(p, h);

    const u128 T =
        (u128)n - 1;

    // Truncate both n-1 and m at digit h.
    const u64 n_high =
        (u64)(T / ph);

    const u64 m_high =
        (u64)((u128)m / ph);

    // Number of admissible high-prefix values <= n_high.
    const u128 A =
        direct_miss_prefix(
            p,
            m_high,
            n_high
        );

    // Old lower digits had p^h possibilities,
    // but only the all-zero lower suffix survives.
    const u128 lower_nonzero =
        ph - 1;

    return A * lower_nonzero;
}

// ============================================================
// Predicted HIT-prefix delta.
//
// HIT = total - MISS.
//
// Therefore:
//
// C_m(n) - C_{m+1}(n)
//     = MISS_{m+1}(n) - MISS_m(n)
//     = gain - loss.
// ============================================================

u128 predicted_hit_delta(
    u64 p,
    u64 m,
    u64 n
) {
    const u128 gain =
        predicted_gain(
            p,
            m,
            n
        );

    const u128 loss =
        predicted_loss(
            p,
            m,
            n
        );

    if (gain >= loss) {
        return gain - loss;
    }

    // The HIT delta can be negative.
    // This function is therefore only used for equality
    // after comparing signed differences elsewhere.
    return 0;
}

// ============================================================
// Brute-force gain/loss
//
// Small cases only.
// ============================================================

u64 brute_gain(
    u64 p,
    u64 m,
    u64 n
) {
    const u64 mp1 =
        m + 1;

    u64 count = 0;

    for (u64 x = 0;
         x < n;
         ++x) {

        if (is_miss(p, mp1, x) &&
            !is_miss(p, m, x)) {

            ++count;
        }
    }

    return count;
}

u64 brute_loss(
    u64 p,
    u64 m,
    u64 n
) {
    const u64 mp1 =
        m + 1;

    u64 count = 0;

    for (u64 x = 0;
         x < n;
         ++x) {

        if (is_miss(p, m, x) &&
            !is_miss(p, mp1, x)) {

            ++count;
        }
    }

    return count;
}

// ============================================================
// Signed comparison of HIT delta
// ============================================================

bool verify_point(
    u64 p,
    u64 m,
    u64 n
) {
    const u64 mp1 =
        m + 1;

    const u128 old_hit =
        direct_hit_prefix(
            p,
            m,
            n
        );

    const u128 new_hit =
        direct_hit_prefix(
            p,
            mp1,
            n
        );

    const u128 gain =
        predicted_gain(
            p,
            m,
            n
        );

    const u128 loss =
        predicted_loss(
            p,
            m,
            n
        );

    if (gain >= loss) {
        return
            old_hit - new_hit
            ==
            gain - loss;
    }

    return
        new_hit - old_hit
        ==
        loss - gain;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout << "START EXPERIMENT 265\n";
    std::cout << "ARBITRARY m -> m+1 CARRY RECURRENCE\n";
    std::cout << "GAINED MISS SLICE VS LOST MISS SLICES\n";
    std::cout << "AND RESULTING HIT-PREFIX DELTA\n\n";

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 0},
        {2, 1},
        {2, 2},
        {2, 3},
        {2, 4},
        {2, 6},
        {2, 7},
        {2, 10},
        {3, 8},
        {3, 10},
        {3, 26},
        {5, 24},
        {5, 124},
        {7, 999},
        {11, 12345},
        {13, 987654321},
        {17, 1000000}
    };

    std::size_t deterministic_fail = 0;

    for (const auto& [p, m] : cases) {

        bool ok = true;

        const std::vector<u64> tests = {
            0,
            1,
            2,
            m / 4,
            m / 2,
            m,
            m + 1,
            m + 2
        };

        for (u64 n : tests) {

            if (!verify_point(
                    p,
                    m,
                    n)) {

                ok = false;
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
        << " fail="
        << deterministic_fail
        << "\n\n";

    // --------------------------------------------------------
    // Exhaustive small recurrence
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

            // Avoid overflow and enormous work.
            if (m == UINT64_MAX) {
                continue;
            }

            for (u64 n = 0;
                 n <= m + 2;
                 ++n) {

                ++exhaustive_cases;

                if (!verify_point(
                        p,
                        m,
                        n)) {

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
    // Independent brute-force decomposition
    // --------------------------------------------------------

    std::size_t brute_cases = 0;
    std::size_t brute_fail = 0;

    for (u64 p : {
        2ULL,
        3ULL,
        5ULL
    }) {

        for (u64 m = 0;
             m <= 150;
             ++m) {

            for (u64 n = 0;
                 n <= m + 2;
                 ++n) {

                ++brute_cases;

                const u128 gain =
                    predicted_gain(
                        p,
                        m,
                        n
                    );

                const u128 loss =
                    predicted_loss(
                        p,
                        m,
                        n
                    );

                const u64 brute_g =
                    brute_gain(
                        p,
                        m,
                        n
                    );

                const u64 brute_l =
                    brute_loss(
                        p,
                        m,
                        n
                    );

                if (gain != (u128)brute_g ||
                    loss != (u128)brute_l) {

                    ++brute_fail;
                }
            }
        }
    }

    std::cout
        << "brute_gain_loss_cases="
        << brute_cases
        << " fail="
        << brute_fail
        << "\n\n";

    // --------------------------------------------------------
    // Random large cases
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x265265265ULL
    );

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        41, 43, 47
    };

    const std::size_t random_cases = 100000;
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

        if (!verify_point(
                p,
                m,
                n)) {

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
    // Large case
    // --------------------------------------------------------

    {
        const u64 p = 13;

        const u64 m =
            987654321012345678ULL;

        const u64 n =
            876543210123456789ULL;

        const u64 mp1 =
            m + 1;

        const u128 old_hit =
            direct_hit_prefix(
                p,
                m,
                n
            );

        const u128 new_hit =
            direct_hit_prefix(
                p,
                mp1,
                n
            );

        const u128 gain =
            predicted_gain(
                p,
                m,
                n
            );

        const u128 loss =
            predicted_loss(
                p,
                m,
                n
            );

        std::cout << "LARGE_CASE\n";
        std::cout << "p=" << p << '\n';
        std::cout << "m=" << m << '\n';
        std::cout << "m_plus_1=" << mp1 << '\n';
        std::cout << "n=" << n << '\n';

        std::cout << "old_hit=";
        print_u128(old_hit);

        std::cout << "\nnew_hit=";
        print_u128(new_hit);

        std::cout << "\ngained_miss=";
        print_u128(gain);

        std::cout << "\nlost_miss=";
        print_u128(loss);

        std::cout
            << "\nactual_hit_delta=";

        if (old_hit >= new_hit) {
            print_u128(
                old_hit - new_hit
            );
        } else {
            std::cout << '-';
            print_u128(
                new_hit - old_hit
            );
        }

        std::cout
            << "\npredicted_hit_delta=";

        if (gain >= loss) {
            print_u128(
                gain - loss
            );
        } else {
            std::cout << '-';
            print_u128(
                loss - gain
            );
        }

        std::cout
            << "\nlarge_case_pass="
            << (
                verify_point(
                    p,
                    m,
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
        brute_fail == 0 &&
        random_fail == 0;

    std::cout
        << "OVERALL PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 265\n";

    return overall ? 0 : 1;
}
