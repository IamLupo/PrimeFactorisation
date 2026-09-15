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
// Digits -> value
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
    const std::vector<u64>& md
) {
    std::vector<u128> W(
        md.size() + 1
    );

    W[0] = 1;

    for (std::size_t i = 0;
         i < md.size();
         ++i) {

        W[i + 1] =
            W[i] * (u128)(md[i] + 1);
    }

    return W;
}

// ============================================================
// Digitwise MISS test
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
// Direct MISS prefix
//
// Number of MISS values x with 0 <= x <= y.
// ============================================================

u128 direct_miss_prefix(
    u64 p,
    u64 m,
    u64 y
) {
    const auto yd = base_digits(p, y);
    const auto md = base_digits(p, m);

    const std::size_t L =
        std::max(
            yd.size(),
            md.size()
        );

    const auto W =
        mixed_weights(md);

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
// Direct HIT prefix
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

    const u128 miss =
        direct_miss_prefix(
            p,
            m,
            n - 1
        );

    return (u128)n - miss;
}

// ============================================================
// Power
// ============================================================

u64 p_power_u64(
    u64 p,
    std::size_t r
) {
    u64 result = 1;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        result *= p;
    }

    return result;
}

// ============================================================
// Product of lower digit radices
//
// L_r = product_{i<r}(m_i+1)
// ============================================================

u128 lower_product(
    const std::vector<u64>& md,
    std::size_t r
) {
    u128 result = 1;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        result *=
            (u128)(md[i] + 1);
    }

    return result;
}

// ============================================================
// Lower m = m mod p^r
// ============================================================

u64 lower_part(
    u64 p,
    u64 m,
    std::size_t r
) {
    const u64 pr =
        p_power_u64(p, r);

    return m % pr;
}

// ============================================================
// Lower n = n mod p^r
// ============================================================

u64 lower_part_n(
    u64 p,
    u64 n,
    std::size_t r
) {
    const u64 pr =
        p_power_u64(p, r);

    return n % pr;
}

// ============================================================
// CORRECTED NEW MISS SLICE
//
// m' = m + p^r
//
// Assume m_r < p-1, so there is no carry.
//
// New MISS values satisfy:
//
//   x_r = m_r + 1
//
//   x_i <= m_i, i != r
//
// For x < n:
//
//   higher(x) < higher(n):
//       all lower digits are possible.
//
//   higher(x) = higher(n):
//
//       if n_r > m_r+1:
//           all lower digits are possible.
//
//       if n_r = m_r+1:
//           only lower values < n_low.
//
//       if n_r < m_r+1:
//           nothing.
//
// ============================================================

u128 new_miss_slice_count(
    u64 p,
    u64 m,
    std::size_t r,
    u64 n
) {
    const auto md =
        base_digits(p, m);

    const u64 pr =
        p_power_u64(p, r);

    const u64 pr1 =
        pr * p;

    const u64 higher_n =
        n / pr1;

    const u64 digit_n =
        (n / pr) % p;

    const u64 digit_m =
        md[r];

    const u64 new_digit =
        digit_m + 1;

    const u64 higher_m =
        m / pr1;

    const u128 lower_total =
        lower_product(md, r);

    // --------------------------------------------------------
    // Higher prefixes strictly below higher_n.
    // --------------------------------------------------------

    u128 higher_less = 0;

    if (higher_n > 0) {
        higher_less =
            direct_miss_prefix(
                p,
                higher_m,
                higher_n - 1
            );
    }

    u128 result =
        higher_less * lower_total;

    // --------------------------------------------------------
    // Higher prefix exactly equal to higher_n.
    //
    // It must itself be a valid MISS higher prefix.
    // --------------------------------------------------------

    const bool higher_equal_valid =
        is_miss(
            p,
            higher_m,
            higher_n
        );

    if (!higher_equal_valid) {
        return result;
    }

    // --------------------------------------------------------
    // Case 1:
    //
    // n_r > m_r+1
    //
    // Entire lower slice is below n.
    // --------------------------------------------------------

    if (digit_n > new_digit) {
        result += lower_total;
        return result;
    }

    // --------------------------------------------------------
    // Case 2:
    //
    // n_r < m_r+1
    //
    // New slice starts after n.
    // --------------------------------------------------------

    if (digit_n < new_digit) {
        return result;
    }

    // --------------------------------------------------------
    // Case 3:
    //
    // n_r == m_r+1
    //
    // Compare lower digits.
    // --------------------------------------------------------

    const u64 n_low =
        lower_part_n(
            p,
            n,
            r
        );

    if (n_low == 0) {
        return result;
    }

    const u64 m_low =
        lower_part(
            p,
            m,
            r
        );

    // Count lower MISS values < n_low.
    //
    // The lower digit system is bounded by m_low.
    //
    // n_low itself need not be MISS.
    const u128 partial_lower =
        direct_miss_prefix(
            p,
            m_low,
            n_low - 1
        );

    result += partial_lower;

    return result;
}

// ============================================================
// Brute-force new slice
//
// Used only for small cases.
// ============================================================

u64 brute_new_slice(
    u64 p,
    u64 m,
    std::size_t r,
    u64 n
) {
    const u64 step =
        p_power_u64(p, r);

    const u64 mprime =
        m + step;

    u64 count = 0;

    for (u64 x = 0;
         x < n;
         ++x) {

        if (is_miss(p, mprime, x) &&
            !is_miss(p, m, x)) {

            ++count;
        }
    }

    return count;
}

// ============================================================
// Verify recurrence
// ============================================================

bool verify_point(
    u64 p,
    u64 m,
    std::size_t r,
    u64 n
) {
    const auto md =
        base_digits(p, m);

    if (r >= md.size()) {
        return false;
    }

    if (md[r] >= p - 1) {
        return true;
    }

    const u64 step =
        p_power_u64(p, r);

    const u64 mprime =
        m + step;

    const u128 old_hit =
        direct_hit_prefix(
            p,
            m,
            n
        );

    const u128 new_hit =
        direct_hit_prefix(
            p,
            mprime,
            n
        );

    const u128 delta =
        old_hit - new_hit;

    const u128 expected =
        new_miss_slice_count(
            p,
            m,
            r,
            n
        );

    return delta == expected;
}

// ============================================================
// Main
// ============================================================

int main() {
    std::cout << "START EXPERIMENT 262\n";
    std::cout << "CORRECTED LOCAL BASE-p DIGIT UPDATE\n";
    std::cout << "m -> m + p^r WITHOUT CARRY\n";
    std::cout << "ALL THREE r-DIGIT PREFIX CASES\n\n";

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> cases = {
        {2, 2},
        {2, 4},
        {2, 10},
        {3, 10},
        {3, 80},
        {5, 124},
        {7, 999},
        {11, 12345},
        {13, 987654321},
        {17, 1000000},
        {19, 123456789},
        {31, 987654321012ULL}
    };

    std::size_t deterministic_fail = 0;

    for (const auto& [p, m] : cases) {
        const auto md =
            base_digits(p, m);

        bool ok = true;

        for (std::size_t r = 0;
             r < md.size();
             ++r) {

            if (md[r] >= p - 1) {
                continue;
            }

            const u64 step =
                p_power_u64(p, r);

            const u64 mprime =
                m + step;

            const std::vector<u64> tests = {
                0,
                1,
                2,
                m / 4,
                m / 2,
                m,
                mprime - 1,
                mprime,
                mprime + 1
            };

            for (u64 n : tests) {
                if (!verify_point(
                        p,
                        m,
                        r,
                        n)) {

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
    // Exhaustive recurrence test
    // --------------------------------------------------------

    std::size_t exhaustive_cases = 0;
    std::size_t exhaustive_fail = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL}) {
        for (u64 m = 0;
             m <= 500;
             ++m) {

            const auto md =
                base_digits(p, m);

            for (std::size_t r = 0;
                 r < md.size();
                 ++r) {

                if (md[r] >= p - 1) {
                    continue;
                }

                const u64 mprime =
                    m + p_power_u64(p, r);

                for (u64 n = 0;
                     n <= mprime + 1;
                     ++n) {

                    ++exhaustive_cases;

                    if (!verify_point(
                            p,
                            m,
                            r,
                            n)) {

                        ++exhaustive_fail;
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
    // Independent brute-force validation
    // --------------------------------------------------------

    std::size_t brute_cases = 0;
    std::size_t brute_fail = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (u64 m = 0;
             m <= 100;
             ++m) {

            const auto md =
                base_digits(p, m);

            for (std::size_t r = 0;
                 r < md.size();
                 ++r) {

                if (md[r] >= p - 1) {
                    continue;
                }

                const u64 mprime =
                    m + p_power_u64(p, r);

                for (u64 n = 0;
                     n <= mprime + 1;
                     ++n) {

                    ++brute_cases;

                    const u128 formula =
                        new_miss_slice_count(
                            p,
                            m,
                            r,
                            n
                        );

                    const u64 brute =
                        brute_new_slice(
                            p,
                            m,
                            r,
                            n
                        );

                    if (formula !=
                        (u128)brute) {

                        ++brute_fail;
                    }
                }
            }
        }
    }

    std::cout
        << "brute_slice_cases="
        << brute_cases
        << " fail="
        << brute_fail
        << "\n\n";

    // --------------------------------------------------------
    // Random
    // --------------------------------------------------------

    std::mt19937_64 rng(
        0x262262262ULL
    );

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37, 41,
        43, 47
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

        const auto md =
            base_digits(p, m);

        std::vector<std::size_t> valid;

        for (std::size_t r = 0;
             r < md.size();
             ++r) {

            if (md[r] < p - 1) {
                valid.push_back(r);
            }
        }

        if (valid.empty()) {
            continue;
        }

        const std::size_t r =
            valid[rng() % valid.size()];

        const u64 step =
            p_power_u64(p, r);

        const u64 mprime =
            m + step;

        // Include values deliberately targeting
        // all three r-digit boundary cases.
        const u64 pr =
            step;

        const u64 nr =
            (mprime / pr) % p;

        const u64 lower_n =
            rng() % pr;

        const u64 higher_n =
            rng() % (
                mprime / (pr * p) + 3
            );

        const std::vector<u64> tests = {
            0,
            1,
            m / 3,
            m / 2,
            m,
            mprime - 1,
            mprime,
            mprime + 1,

            // r-digit equal case with arbitrary lower part.
            higher_n * pr * p
                + (u64)nr * pr
                + lower_n
        };

        for (u64 n : tests) {
            if (!verify_point(
                    p,
                    m,
                    r,
                    n)) {

                ++random_fail;
            }
        }
    }

    std::cout
        << "random_cases="
        << random_cases
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
            base_digits(p, m);

        std::size_t r = 0;

        while (r < md.size() &&
               md[r] >= p - 1) {

            ++r;
        }

        if (r == md.size()) {
            r = 0;
        }

        const u64 step =
            p_power_u64(p, r);

        const u64 mprime =
            m + step;

        const u64 n =
            876543210123456789ULL;

        const u128 old_hit =
            direct_hit_prefix(
                p,
                m,
                n
            );

        const u128 new_hit =
            direct_hit_prefix(
                p,
                mprime,
                n
            );

        const u128 delta =
            old_hit - new_hit;

        const u128 slice =
            new_miss_slice_count(
                p,
                m,
                r,
                n
            );

        std::cout << "LARGE_CASE\n";
        std::cout << "p=" << p << '\n';
        std::cout << "m=" << m << '\n';
        std::cout << "r=" << r << '\n';
        std::cout << "m_prime=" << mprime << '\n';
        std::cout << "n=" << n << '\n';

        std::cout << "old_hit=";
        print_u128(old_hit);

        std::cout << "\nnew_hit=";
        print_u128(new_hit);

        std::cout << "\ndelta=";
        print_u128(delta);

        std::cout << "\nnew_miss_slice=";
        print_u128(slice);

        std::cout
            << "\nlarge_case_pass="
            << (delta == slice ? 1 : 0)
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

    std::cout << "FINISHED EXPERIMENT 262\n";

    return overall ? 0 : 1;
}
