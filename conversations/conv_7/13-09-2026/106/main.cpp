#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

// ------------------------------------------------------------
// Output
// ------------------------------------------------------------

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

// ------------------------------------------------------------
// Base-p digits, least significant first
// ------------------------------------------------------------

std::vector<u64> base_digits(u64 p, u64 n) {
    std::vector<u64> d;

    while (n > 0) {
        d.push_back(n % p);
        n /= p;
    }

    if (d.empty()) {
        d.push_back(0);
    }

    return d;
}

// ------------------------------------------------------------
// Digits -> integer
// ------------------------------------------------------------

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

// ------------------------------------------------------------
// Product_i (m_i + 1)
// ------------------------------------------------------------

u128 miss_count(
    const std::vector<u64>& m_digits
) {
    u128 result = 1;

    for (u64 d : m_digits) {
        result *= (u128)(d + 1);
    }

    return result;
}

// ------------------------------------------------------------
// Type count
//
// C_r = m_r * product_{i>r}(m_i+1)
// ------------------------------------------------------------

u128 type_count(
    const std::vector<u64>& m_digits,
    std::size_t r
) {
    u128 result = m_digits[r];

    for (std::size_t i = r + 1; i < m_digits.size(); ++i) {
        result *= (u128)(m_digits[i] + 1);
    }

    return result;
}

// ------------------------------------------------------------
// p^r
// ------------------------------------------------------------

u128 p_power(
    u64 p,
    std::size_t r
) {
    u128 result = 1;

    for (std::size_t i = 0; i < r; ++i) {
        result *= (u128)p;
    }

    return result;
}

// ------------------------------------------------------------
// sum_{i<r} m_i p^i
// ------------------------------------------------------------

u128 lower_value(
    u64 p,
    const std::vector<u64>& m_digits,
    std::size_t r
) {
    u128 result = 0;
    u128 power = 1;

    for (std::size_t i = 0; i < r; ++i) {
        result += (u128)m_digits[i] * power;
        power *= (u128)p;
    }

    return result;
}

// ------------------------------------------------------------
// G_r
//
// G_r = p^r - 1 - sum_{i<r} m_i p^i
// ------------------------------------------------------------

u128 gap_length(
    u64 p,
    const std::vector<u64>& m_digits,
    std::size_t r
) {
    return p_power(p, r)
         - 1
         - lower_value(p, m_digits, r);
}

// ------------------------------------------------------------
// Construct the FIRST mixed-radix element of type r.
//
// Digits:
//   i < r : m_i
//   i = r : 0
//   i > r : arbitrary
//
// Here we choose all upper digits = 0.
// ------------------------------------------------------------

u128 first_type_element(
    u64 p,
    const std::vector<u64>& m_digits,
    std::size_t r
) {
    std::vector<u64> d(m_digits.size(), 0);

    for (std::size_t i = 0; i < r; ++i) {
        d[i] = m_digits[i];
    }

    d[r] = 0;

    return digits_to_value(p, d);
}

// ------------------------------------------------------------
// Construct the LAST mixed-radix element of type r.
//
// Digits:
//   i < r : m_i
//   i = r : m_r - 1
//   i > r : m_i
// ------------------------------------------------------------

u128 last_type_element(
    u64 p,
    const std::vector<u64>& m_digits,
    std::size_t r
) {
    std::vector<u64> d = m_digits;

    for (std::size_t i = 0; i < r; ++i) {
        d[i] = m_digits[i];
    }

    d[r] = m_digits[r] - 1;

    return digits_to_value(p, d);
}

// ------------------------------------------------------------
// Exact mixed-radix successor.
//
// The input is represented by bounded digits:
//     0 <= d_i <= m_i
//
// Find the lowest digit that is not maximal, increment it,
// and reset all lower digits to zero.
// ------------------------------------------------------------

u128 successor_digits(
    u64 p,
    const std::vector<u64>& m_digits,
    std::vector<u64> d
) {
    for (std::size_t i = 0; i < d.size(); ++i) {
        if (d[i] < m_digits[i]) {
            ++d[i];

            for (std::size_t j = 0; j < i; ++j) {
                d[j] = 0;
            }

            return digits_to_value(p, d);
        }

        d[i] = 0;
    }

    return digits_to_value(p, d);
}

// ------------------------------------------------------------
// Verify one type.
//
// Only algebraic construction is used.
// ------------------------------------------------------------

bool verify_type(
    u64 p,
    u64 m,
    const std::vector<u64>& digits,
    std::size_t r
) {
    if (digits[r] == 0) {
        return true;
    }

    const u128 expected_gap =
        gap_length(p, digits, r);

    // First element of the type.
    {
        const u128 x =
            first_type_element(p, digits, r);

        std::vector<u64> xd =
            base_digits(p, (u64)x);

        xd.resize(digits.size(), 0);

        const u128 y =
            successor_digits(p, digits, xd);

        const u128 actual_gap =
            y - x - 1;

        if (actual_gap != expected_gap) {
            return false;
        }
    }

    // Last element of the type.
    {
        const u128 x =
            last_type_element(p, digits, r);

        std::vector<u64> xd =
            base_digits(p, (u64)x);

        xd.resize(digits.size(), 0);

        const u128 y =
            successor_digits(p, digits, xd);

        const u128 actual_gap =
            y - x - 1;

        if (actual_gap != expected_gap) {
            return false;
        }
    }

    return true;
}

// ------------------------------------------------------------
// Verify complete algebraic identity
// ------------------------------------------------------------

bool verify_identity(
    u64 p,
    u64 m
) {
    const auto digits = base_digits(p, m);

    const u128 product =
        miss_count(digits);

    u128 count_sum = 0;
    u128 weighted_sum = 0;

    for (std::size_t r = 0; r < digits.size(); ++r) {
        const u128 c =
            type_count(digits, r);

        const u128 g =
            gap_length(p, digits, r);

        count_sum += c;
        weighted_sum += c * g;
    }

    const u128 expected_count =
        product - 1;

    const u128 expected_hits =
        (u128)m + 1 - product;

    return
        count_sum == expected_count &&
        weighted_sum == expected_hits;
}

// ------------------------------------------------------------
// Exhaustive mixed-radix enumeration for small m
//
// This directly walks the MISS set.
// ------------------------------------------------------------

bool exhaustive_small(
    u64 p,
    u64 m
) {
    const auto digits = base_digits(p, m);

    std::vector<u64> d(digits.size(), 0);
    std::vector<u64> observed(digits.size(), 0);

    while (true) {
        bool is_m = true;

        for (std::size_t i = 0; i < digits.size(); ++i) {
            if (d[i] != digits[i]) {
                is_m = false;
                break;
            }
        }

        if (!is_m) {
            std::size_t r = 0;

            while (r < digits.size() &&
                   d[r] == digits[r]) {
                ++r;
            }

            if (r >= digits.size()) {
                return false;
            }

            ++observed[r];
        }

        // Mixed-radix successor.
        std::size_t i = 0;

        while (i < d.size() && d[i] == digits[i]) {
            d[i] = 0;
            ++i;
        }

        if (i == d.size()) {
            break;
        }

        ++d[i];
    }

    for (std::size_t r = 0; r < digits.size(); ++r) {
        if ((u128)observed[r] !=
            type_count(digits, r)) {
            return false;
        }
    }

    return true;
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------

int main() {
    std::cout << "START EXPERIMENT 258\n";
    std::cout << "ARBITRARY-p,m MIXED-RADIX CARRY TYPES\n";
    std::cout << "DIRECT GAP REPRESENTATIVE VERIFICATION\n\n";

    // --------------------------------------------------------
    // Deterministic
    // --------------------------------------------------------

    const std::vector<std::pair<u64, u64>> deterministic = {
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
        {17, 1000000}
    };

    std::size_t deterministic_fail = 0;

    for (const auto& [p, m] : deterministic) {
        const auto digits =
            base_digits(p, m);

        bool ok =
            verify_identity(p, m);

        for (std::size_t r = 0;
             r < digits.size();
             ++r) {

            if (!verify_type(
                    p, m, digits, r)) {
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
        << deterministic.size()
        << " fail="
        << deterministic_fail
        << "\n\n";

    // --------------------------------------------------------
    // Random large cases
    // --------------------------------------------------------

    std::mt19937_64 rng(0x258258258ULL);

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

        bool ok =
            verify_identity(p, m);

        const auto digits =
            base_digits(p, m);

        for (std::size_t r = 0;
             r < digits.size();
             ++r) {

            if (!verify_type(
                    p, m, digits, r)) {
                ok = false;
                break;
            }
        }

        if (!ok) {
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
    // Exhaustive small cases
    // --------------------------------------------------------

    const std::size_t exhaustive_cases = 500;
    std::size_t exhaustive_fail = 0;

    for (std::size_t tc = 0;
         tc < exhaustive_cases;
         ++tc) {

        const u64 p =
            primes[rng() % primes.size()];

        const u64 m =
            rng() % 50000;

        if (!exhaustive_small(p, m)) {
            ++exhaustive_fail;
        }
    }

    std::cout
        << "exhaustive_small_cases="
        << exhaustive_cases
        << " fail="
        << exhaustive_fail
        << "\n\n";

    // --------------------------------------------------------
    // Large symbolic case
    // --------------------------------------------------------

    {
        const u64 p = 13;
        const u64 m =
            987654321012345678ULL;

        const auto digits =
            base_digits(p, m);

        const u128 product =
            miss_count(digits);

        u128 count_sum = 0;
        u128 weighted_sum = 0;

        std::cout << "LARGE_CASE\n";
        std::cout << "p=" << p << '\n';
        std::cout << "m=" << m << '\n';
        std::cout << "digits=" << digits.size() << '\n';

        for (std::size_t r = 0;
             r < digits.size();
             ++r) {

            const u128 c =
                type_count(digits, r);

            const u128 g =
                gap_length(p, digits, r);

            count_sum += c;
            weighted_sum += c * g;

            std::cout
                << "r=" << r
                << " m_r=" << digits[r]
                << " C_r=";

            print_u128(c);

            std::cout
                << " G_r=";

            print_u128(g);

            std::cout << '\n';
        }

        const u128 expected_count =
            product - 1;

        const u128 expected_hits =
            (u128)m + 1 - product;

        std::cout << "type_total=";
        print_u128(count_sum);

        std::cout << "\nexpected_gap_count=";
        print_u128(expected_count);

        std::cout << "\nweighted_total=";
        print_u128(weighted_sum);

        std::cout << "\nexpected_hits=";
        print_u128(expected_hits);

        std::cout
            << "\nlarge_case_pass="
            << (
                count_sum == expected_count &&
                weighted_sum == expected_hits
                    ? 1
                    : 0
            )
            << "\n\n";
    }

    const bool overall =
        deterministic_fail == 0 &&
        random_fail == 0 &&
        exhaustive_fail == 0;

    std::cout
        << "OVERALL PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout << "FINISHED EXPERIMENT 258\n";

    return overall ? 0 : 1;
}
