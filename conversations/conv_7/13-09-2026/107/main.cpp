#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

// ============================================================
// OUTPUT
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
// BASE-p DIGITS
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
// DIGITS -> VALUE
// ============================================================

u128 digits_to_value(
    u64 p,
    const std::vector<u64>& d
) {
    u128 value = 0;
    u128 power = 1;

    for (u64 x : d) {
        value += (u128)x * power;
        power *= (u128)p;
    }

    return value;
}

// ============================================================
// MIXED-RADIX WEIGHTS
//
// W_r = product_{i<r}(m_i+1)
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
// MISS RANK
//
// R(x) = sum_i x_i W_i
// ============================================================

u128 mixed_rank(
    const std::vector<u64>& m_digits,
    const std::vector<u128>& W,
    const std::vector<u64>& x_digits
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < m_digits.size();
         ++i) {

        rank +=
            (u128)x_digits[i] * W[i];
    }

    return rank;
}

// ============================================================
// DIGITWISE MISS TEST
//
// x is MISS iff x_i <= m_i for every digit.
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
        const u64 a =
            (i < xd.size()) ? xd[i] : 0;

        const u64 b =
            (i < md.size()) ? md[i] : 0;

        if (a > b) {
            return false;
        }
    }

    return true;
}

// ============================================================
// PREDECESSOR MISS
//
// Greatest MISS value <= n.
//
// If n > m, this is simply m.
//
// If n is already MISS, it is n.
//
// Otherwise find the most-significant digit where
// n_i > m_i and clamp that digit to m_i while setting
// all lower digits to m_i.
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

    for (std::size_t i = 0; i < L; ++i) {
        const u64 a =
            (i < nd.size()) ? nd[i] : 0;

        const u64 b =
            (i < md.size()) ? md[i] : 0;

        x[i] = a;

        if (a > b) {
            // Find the most-significant offending digit.
            std::size_t h = i;

            for (std::size_t j = i + 1;
                 j < L;
                 ++j) {

                const u64 aj =
                    (j < nd.size()) ? nd[j] : 0;

                const u64 bj =
                    (j < md.size()) ? md[j] : 0;

                if (aj > bj) {
                    h = j;
                }
            }

            // Digits above h remain equal to n.
            // Digit h becomes m_h.
            // Digits below h become m_i.
            for (std::size_t j = 0; j < h; ++j) {
                x[j] =
                    (j < md.size()) ? md[j] : 0;
            }

            x[h] =
                (h < md.size()) ? md[h] : 0;

            for (std::size_t j = h + 1;
                 j < L;
                 ++j) {

                x[j] =
                    (j < nd.size()) ? nd[j] : 0;
            }

            return (u64)digits_to_value(p, x);
        }
    }

    return n;
}

// ============================================================
// CARRY-GAP LENGTH
//
// G_r = p^r - 1 - sum_{i<r} m_i p^i
// ============================================================

u128 gap_length(
    u64 p,
    const std::vector<u64>& m_digits,
    std::size_t r
) {
    u128 power = 1;
    u128 lower = 0;

    for (std::size_t i = 0; i < r; ++i) {
        lower +=
            (u128)m_digits[i] * power;

        power *= (u128)p;
    }

    return power - 1 - lower;
}

// ============================================================
// DIRECT DIGITWISE MISS PREFIX COUNT
//
// Number of x with 0 <= x <= y and x_i <= m_i.
//
// This is an independent digit-DP-style count.
// ============================================================

u128 direct_miss_prefix(
    u64 p,
    u64 m,
    u64 y
) {
    if (y < 0) {
        return 0;
    }

    const auto yd = base_digits(p, y);
    const auto md = base_digits(p, m);

    const std::size_t L =
        std::max(yd.size(), md.size());

    // Lower mixed-radix products.
    std::vector<u128> lower(L + 1);
    lower[0] = 1;

    for (std::size_t i = 0; i < L; ++i) {
        const u64 mi =
            (i < md.size()) ? md[i] : 0;

        lower[i + 1] =
            lower[i] * (u128)(mi + 1);
    }

    u128 count = 0;

    // Scan from most-significant digit downward.
    for (std::size_t pos = L; pos-- > 0;) {
        const u64 yi =
            (pos < yd.size()) ? yd[pos] : 0;

        const u64 mi =
            (pos < md.size()) ? md[pos] : 0;

        // Choices strictly below y_i while remaining
        // compatible with m_i.
        const u64 choices =
            std::min<u64>(yi, mi + 1);

        if (choices > 0) {
            count +=
                (u128)choices * lower[pos];
        }

        // If y_i > m_i, the tight prefix can no longer
        // continue into the MISS set.
        if (yi > mi) {
            return count;
        }

        // Otherwise continue with y_i itself.
    }

    // y itself is MISS.
    return count + 1;
}

// ============================================================
// GAP-TYPE COUNT BELOW RANK K
//
// A rank k has carry type r iff:
//
//   k_i = m_i for i < r
//   k_r < m_r
//
// For ranks 0 <= k < K:
//
// N_r(K) = floor(K / W_{r+1}) * m_r
//        + min(m_r, floor((K mod W_{r+1}) / W_r))
// ============================================================

u128 type_count_below_rank(
    const std::vector<u64>& m_digits,
    const std::vector<u128>& W,
    std::size_t r,
    u128 K
) {
    const u128 block =
        W[r + 1];

    const u128 full_blocks =
        K / block;

    const u128 remainder =
        K % block;

    u128 count =
        full_blocks * (u128)m_digits[r];

    const u128 partial =
        remainder / W[r];

    count +=
        std::min<u128>(
            (u128)m_digits[r],
            partial
        );

    return count;
}

// ============================================================
// GAP-BASED PREFIX HIT COUNT
//
// Counts HIT values x with 0 <= x < n.
//
// Find M_K = predecessor MISS(n).
// All complete gaps before M_K are counted algebraically.
// Then count the partial current gap.
// ============================================================

u128 gap_prefix_hit_count(
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
        mixed_rank(md, W, padded);

    u128 result = 0;

    // All complete gaps whose left endpoint has
    // mixed-radix rank < K.
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

    // Partial current gap.
    //
    // pred is the left MISS endpoint.
    // HIT values start at pred+1.
    //
    // Only values strictly below n count.
    if ((u64)pred < n) {
        result +=
            (u128)(n - pred - 1);
    }

    return result;
}

// ============================================================
// DIRECT HIT PREFIX COUNT
//
// HIT values are complement of MISS values.
// ============================================================

u128 direct_hit_prefix(
    u64 p,
    u64 m,
    u64 n
) {
    if (n == 0) {
        return 0;
    }

    // Number of integers x with 0 <= x < n.
    const u128 total =
        (u128)n;

    // Number of MISS values x with 0 <= x < n.
    const u128 miss =
        direct_miss_prefix(
            p, m, n - 1
        );

    return total - miss;
}

// ============================================================
// VERIFY ONE CASE
// ============================================================

bool verify_case(
    u64 p,
    u64 m,
    u64 n
) {
    const u128 direct =
        direct_hit_prefix(p, m, n);

    const u128 reconstructed =
        gap_prefix_hit_count(p, m, n);

    return direct == reconstructed;
}

// ============================================================
// MAIN
// ============================================================

int main() {
    std::cout << "START EXPERIMENT 259\n";
    std::cout << "PREFIX HIT COUNT FROM CARRY-TYPE GAPS\n";
    std::cout << "ARBITRARY p,m,n\n\n";

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
        std::vector<u64> tests = {
            0,
            1,
            2,
            m / 4,
            m / 2,
            m > 0 ? m - 1 : 0,
            m,
            m < UINT64_MAX ? m + 1 : m
        };

        bool ok = true;

        for (u64 n : tests) {
            if (n > m + 1) {
                continue;
            }

            if (!verify_case(p, m, n)) {
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
    // Exhaustive small cases
    // --------------------------------------------------------

    std::size_t exhaustive_cases = 0;
    std::size_t exhaustive_fail = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL}) {
        for (u64 m = 0; m <= 1000; ++m) {
            for (u64 n = 0; n <= m + 1; ++n) {
                ++exhaustive_cases;

                if (!verify_case(p, m, n)) {
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
    // Random cases
    // --------------------------------------------------------

    std::mt19937_64 rng(0x259259259ULL);

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

        if (!verify_case(p, m, n)) {
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
    // Rank-type demonstration
    // --------------------------------------------------------

    {
        const u64 p = 13;
        const u64 m =
            987654321012345678ULL;

        const auto md =
            base_digits(p, m);

        const auto W =
            mixed_weights(md);

        const u64 n =
            876543210123456789ULL;

        const u64 pred =
            predecessor_miss(p, m, n);

        const auto pd =
            base_digits(p, pred);

        std::vector<u64> padded =
            pd;

        padded.resize(md.size(), 0);

        const u128 K =
            mixed_rank(md, W, padded);

        const u128 direct =
            direct_hit_prefix(p, m, n);

        const u128 reconstructed =
            gap_prefix_hit_count(p, m, n);

        std::cout << "LARGE_CASE\n";
        std::cout << "p=" << p << '\n';
        std::cout << "m=" << m << '\n';
        std::cout << "n=" << n << '\n';
        std::cout << "predecessor_miss=" << pred << '\n';

        std::cout << "pred_rank=";
        print_u128(K);

        std::cout << "\ndirect_hit_prefix=";
        print_u128(direct);

        std::cout << "\ngap_reconstructed=";
        print_u128(reconstructed);

        std::cout
            << "\nlarge_case_pass="
            << (direct == reconstructed ? 1 : 0)
            << "\n\n";
    }

    const bool overall =
        deterministic_fail == 0 &&
        exhaustive_fail == 0 &&
        random_fail == 0;

    std::cout
        << "OVERALL PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout << "FINISHED EXPERIMENT 259\n";

    return overall ? 0 : 1;
}
