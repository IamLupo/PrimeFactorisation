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

std::vector<u64> base_digits(u64 p, u64 n) {
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
// Digits -> integer
// ============================================================

u128 digits_to_value(
    u64 p,
    const std::vector<u64>& d
) {
    u128 result = 0;
    u128 power = 1;

    for (u64 x : d) {
        result += (u128)x * power;
        power *= (u128)p;
    }

    return result;
}

// ============================================================
// Mixed-radix weights
//
// W_i = product_{j<i}(m_j+1)
// ============================================================

std::vector<u128> mixed_weights(
    const std::vector<u64>& m_digits
) {
    std::vector<u128> W(m_digits.size() + 1);

    W[0] = 1;

    for (std::size_t i = 0;
         i < m_digits.size();
         ++i) {

        W[i + 1] =
            W[i] * (u128)(m_digits[i] + 1);
    }

    return W;
}

// ============================================================
// MISS test
//
// x is MISS iff x_i <= m_i for every base-p digit.
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

    for (std::size_t i = 0; i < L; ++i) {
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
// Mixed-radix rank
//
// R(x) = sum x_i W_i
// ============================================================

u128 mixed_rank(
    const std::vector<u128>& W,
    const std::vector<u64>& x_digits
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < x_digits.size();
         ++i) {

        rank +=
            (u128)x_digits[i] * W[i];
    }

    return rank;
}

// ============================================================
// Largest MISS <= n
//
// This is independently constructed from base-p digits.
// ============================================================

u64 predecessor_miss(
    u64 p,
    u64 m,
    u64 n
) {
    if (n >= m) {
        return m;
    }

    const auto nd = base_digits(p, n);
    const auto md = base_digits(p, m);

    const std::size_t L =
        std::max(nd.size(), md.size());

    std::vector<u64> x(L, 0);

    // Find most-significant position h with n_h > m_h.
    bool offending = false;
    std::size_t h = 0;

    for (std::size_t pos = L; pos-- > 0;) {
        const u64 ni =
            pos < nd.size() ? nd[pos] : 0;

        const u64 mi =
            pos < md.size() ? md[pos] : 0;

        if (ni > mi) {
            offending = true;
            h = pos;
            break;
        }
    }

    if (!offending) {
        return n;
    }

    // Higher digits remain equal to n.
    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        x[i] =
            i < nd.size() ? nd[i] : 0;
    }

    // h is clamped to m_h.
    x[h] =
        h < md.size() ? md[h] : 0;

    // Lower digits become maximal MISS digits.
    for (std::size_t i = 0;
         i < h;
         ++i) {

        x[i] =
            i < md.size() ? md[i] : 0;
    }

    return (u64)digits_to_value(p, x);
}

// ============================================================
// Direct closed digit formula for HIT prefix count
//
// Counts HIT values x with 0 <= x < n.
//
// Case A: n <=_p m
//
//   C(n) = n - R(n)
//
// Case B: n not <=_p m
//
// Let h be the most-significant digit with n_h > m_h.
// Then
//
//   C(n)
//     = n - [
//           sum_{i>h} n_i W_i
//           + (m_h+1) W_h
//         ].
//
// ============================================================

u128 closed_hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    const auto nd = base_digits(p, n);
    const auto md = base_digits(p, m);

    const std::size_t L =
        std::max(nd.size(), md.size());

    const auto W =
        mixed_weights(md);

    // Find most-significant offending digit.
    bool offending = false;
    std::size_t h = 0;

    for (std::size_t pos = L; pos-- > 0;) {
        const u64 ni =
            pos < nd.size() ? nd[pos] : 0;

        const u64 mi =
            pos < md.size() ? md[pos] : 0;

        if (ni > mi) {
            offending = true;
            h = pos;
            break;
        }
    }

    u128 miss_prefix = 0;

    if (!offending) {
        // n is itself MISS.
        for (std::size_t i = 0;
             i < L;
             ++i) {

            const u64 ni =
                i < nd.size() ? nd[i] : 0;

            miss_prefix +=
                (u128)ni * W[i];
        }

        // Number of MISS values below n is rank(n).
        return (u128)n - miss_prefix;
    }

    // All digits above h stay equal to n.
    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        const u64 ni =
            i < nd.size() ? nd[i] : 0;

        miss_prefix +=
            (u128)ni * W[i];
    }

    // At h, all values 0..m_h are possible.
    miss_prefix +=
        (u128)(
            (h < md.size() ? md[h] : 0) + 1
        ) * W[h];

    return (u128)n - miss_prefix;
}

// ============================================================
// Direct digitwise MISS prefix count
//
// Number of MISS values x with 0 <= x <= y.
//
// This is used as an independent reference.
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
        mixed_weights(md);

    u128 result = 0;

    for (std::size_t pos = L; pos-- > 0;) {
        const u64 yi =
            pos < yd.size() ? yd[pos] : 0;

        const u64 mi =
            pos < md.size() ? md[pos] : 0;

        const u64 choices =
            std::min<u64>(yi, mi + 1);

        result +=
            (u128)choices * W[pos];

        if (yi > mi) {
            return result;
        }
    }

    return result + 1;
}

// ============================================================
// Independent reference HIT prefix
//
// Number of x in [0,n) minus MISS values in [0,n).
// ============================================================

u128 direct_hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    const u128 miss =
        direct_miss_prefix(
            p, m, n - 1
        );

    return (u128)n - miss;
}

// ============================================================
// GAP-BASED PREFIX COUNT
//
// This reconstructs the result through carry-type gaps.
//
// It is intentionally independent of the closed formula.
// ============================================================

u128 gap_length(
    u64 p,
    const std::vector<u64>& md,
    std::size_t r
) {
    u128 power = 1;
    u128 lower = 0;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        lower +=
            (u128)md[i] * power;

        power *= (u128)p;
    }

    return power - 1 - lower;
}

u128 type_count_below_rank(
    const std::vector<u64>& md,
    const std::vector<u128>& W,
    std::size_t r,
    u128 K
) {
    const u128 block =
        W[r + 1];

    const u128 full =
        K / block;

    const u128 rem =
        K % block;

    u128 result =
        full * (u128)md[r];

    const u128 partial =
        rem / W[r];

    result +=
        std::min<u128>(
            (u128)md[r],
            partial
        );

    return result;
}

u128 gap_hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    const auto md =
        base_digits(p, m);

    const auto W =
        mixed_weights(md);

    const u64 pred =
        predecessor_miss(p, m, n);

    const auto pd =
        base_digits(p, pred);

    std::vector<u64> padded =
        pd;

    padded.resize(md.size(), 0);

    const u128 K =
        mixed_rank(W, padded);

    u128 result = 0;

    for (std::size_t r = 0;
         r < md.size();
         ++r) {

        const u128 C =
            type_count_below_rank(
                md, W, r, K
            );

        const u128 G =
            gap_length(p, md, r);

        result += C * G;
    }

    if (pred < n) {
        result +=
            (u128)(n - pred - 1);
    }

    return result;
}

// ============================================================
// Verify one point against THREE formulas:
//
//   1. direct digitwise reference
//   2. carry-gap reconstruction
//   3. closed digit formula
// ============================================================

bool verify_point(
    u64 p,
    u64 m,
    u64 n
) {
    const u128 a =
        direct_hit_prefix(p, m, n);

    const u128 b =
        gap_hit_prefix(p, m, n);

    const u128 c =
        closed_hit_prefix(p, m, n);

    return
        a == b &&
        a == c;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout << "START EXPERIMENT 260\n";
    std::cout << "CLOSED DIGIT FORMULA FOR HIT PREFIX COUNT\n";
    std::cout << "DIRECT VS GAP-DECOMPOSITION VS CLOSED FORM\n\n";

    // --------------------------------------------------------
    // Deterministic cases
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 0},
        {2, 1},
        {2, 2},
        {2, 3},
        {2, 6},
        {2, 31},
        {3, 10},
        {3, 80},
        {5, 124},
        {7, 999},
        {11, 12345},
        {17, 1000000},
        {13, 987654321012345678ULL}
    };

    std::size_t deterministic_fail = 0;

    for (const auto& [p, m] : cases) {
        bool ok = true;

        const std::vector<u64> tests = {
            0,
            1,
            2,
            m / 5,
            m / 3,
            m / 2,
            m > 0 ? m - 1 : 0,
            m,
            m < UINT64_MAX ? m + 1 : m
        };

        for (u64 n : tests) {
            if (n > m + 1) {
                continue;
            }

            if (!verify_point(p, m, n)) {
                ok = false;
            }
        }

        if (!ok) {
            ++deterministic_fail;
        }

        std::cout
            << "p=" << p
            << " m=" << m
            << " pass=" << (ok ? 1 : 0)
            << '\n';
    }

    std::cout
        << "deterministic_cases="
        << cases.size()
        << " fail="
        << deterministic_fail
        << "\n\n";

    // --------------------------------------------------------
    // Exhaustive small verification
    // --------------------------------------------------------

    std::size_t exhaustive_cases = 0;
    std::size_t exhaustive_fail = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL}) {
        for (u64 m = 0; m <= 1500; ++m) {
            for (u64 n = 0; n <= m + 1; ++n) {
                ++exhaustive_cases;

                if (!verify_point(p, m, n)) {
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
    // Random verification
    // --------------------------------------------------------

    std::mt19937_64 rng(0x260260260ULL);

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37, 41, 43, 47
    };

    const std::size_t random_cases = 100000;
    std::size_t random_fail = 0;

    for (std::size_t tc = 0;
         tc < random_cases;
         ++tc) {

        const u64 p =
            primes[rng() % primes.size()];

        const u64 m =
            rng() % 1000000000000000000ULL;

        const u64 n =
            rng() % (m + 2);

        if (!verify_point(p, m, n)) {
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
    // Large case with explicit formula comparison
    // --------------------------------------------------------

    {
        const u64 p = 13;
        const u64 m =
            987654321012345678ULL;
        const u64 n =
            876543210123456789ULL;

        const u128 direct =
            direct_hit_prefix(p, m, n);

        const u128 gap =
            gap_hit_prefix(p, m, n);

        const u128 closed =
            closed_hit_prefix(p, m, n);

        const u64 pred =
            predecessor_miss(p, m, n);

        std::cout << "LARGE_CASE\n";
        std::cout << "p=" << p << '\n';
        std::cout << "m=" << m << '\n';
        std::cout << "n=" << n << '\n';
        std::cout << "predecessor_miss=" << pred << '\n';

        std::cout << "direct=";
        print_u128(direct);

        std::cout << "\ngap=";
        print_u128(gap);

        std::cout << "\nclosed=";
        print_u128(closed);

        std::cout
            << "\nlarge_case_pass="
            << (
                direct == gap &&
                direct == closed
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

    std::cout << "FINISHED EXPERIMENT 260\n";

    return overall ? 0 : 1;
}
