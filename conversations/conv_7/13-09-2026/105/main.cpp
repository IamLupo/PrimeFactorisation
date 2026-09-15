#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

// ------------------------------------------------------------
// Utility printing
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

std::vector<u64> base_digits(u64 p, u64 m) {
    std::vector<u64> d;

    if (m == 0) {
        d.push_back(0);
        return d;
    }

    while (m > 0) {
        d.push_back(m % p);
        m /= p;
    }

    return d;
}

// ------------------------------------------------------------
// Product of (m_i + 1)
// ------------------------------------------------------------

u128 miss_count(u64 p, u64 m) {
    const auto d = base_digits(p, m);

    u128 product = 1;

    for (u64 x : d) {
        product *= (u128)(x + 1);
    }

    return product;
}

// ------------------------------------------------------------
// Number of carry-type-r gaps
//
// C_r = m_r * product_{i>r}(m_i+1)
// ------------------------------------------------------------

u128 type_count(
    const std::vector<u64>& digits,
    std::size_t r
) {
    u128 result = digits[r];

    for (std::size_t i = r + 1; i < digits.size(); ++i) {
        result *= (u128)(digits[i] + 1);
    }

    return result;
}

// ------------------------------------------------------------
// p^r
// ------------------------------------------------------------

u128 pow_u128(u64 p, std::size_t r) {
    u128 result = 1;

    for (std::size_t i = 0; i < r; ++i) {
        result *= (u128)p;
    }

    return result;
}

// ------------------------------------------------------------
// Lower prefix sum
//
// sum_{i<r} m_i p^i
// ------------------------------------------------------------

u128 lower_value(
    u64 p,
    const std::vector<u64>& digits,
    std::size_t r
) {
    u128 result = 0;
    u128 power = 1;

    for (std::size_t i = 0; i < r; ++i) {
        result += (u128)digits[i] * power;
        power *= (u128)p;
    }

    return result;
}

// ------------------------------------------------------------
// Carry-type gap length
//
// G_r = p^r - 1 - sum_{i<r} m_i p^i
// ------------------------------------------------------------

u128 gap_length(
    u64 p,
    const std::vector<u64>& digits,
    std::size_t r
) {
    return pow_u128(p, r)
         - 1
         - lower_value(p, digits, r);
}

// ------------------------------------------------------------
// Construct x from a digit vector
// ------------------------------------------------------------

u128 digits_to_value(
    u64 p,
    const std::vector<u64>& digits
) {
    u128 result = 0;
    u128 power = 1;

    for (u64 d : digits) {
        result += (u128)d * power;
        power *= (u128)p;
    }

    return result;
}

// ------------------------------------------------------------
// Construct a valid MISS element of carry type r.
//
// lower digits  = m_i
// digit r       = chosen digit < m_r
// upper digits  = chosen values
// ------------------------------------------------------------

u128 construct_type_element(
    u64 p,
    const std::vector<u64>& m_digits,
    std::size_t r,
    u64 digit_r,
    bool upper_max
) {
    std::vector<u64> x = m_digits;

    for (std::size_t i = 0; i < r; ++i) {
        x[i] = m_digits[i];
    }

    x[r] = digit_r;

    for (std::size_t i = r + 1; i < x.size(); ++i) {
        x[i] = upper_max ? m_digits[i] : 0;
    }

    return digits_to_value(p, x);
}

// ------------------------------------------------------------
// Mixed-radix successor of a MISS element.
//
// Radices are m_i + 1.
// ------------------------------------------------------------

u128 mixed_successor(
    u64 p,
    const std::vector<u64>& m_digits,
    u128 x
) {
    std::vector<u64> d = base_digits(p, (u64)x);

    d.resize(m_digits.size(), 0);

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

    return x;
}

// ------------------------------------------------------------
// Determine carry position r from a MISS element x != m.
//
// r is the least significant digit that is not maximal.
// ------------------------------------------------------------

std::size_t carry_position(
    u64 p,
    const std::vector<u64>& m_digits,
    u128 x
) {
    const auto d = base_digits(p, (u64)x);

    for (std::size_t i = 0; i < m_digits.size(); ++i) {
        if (d[i] < m_digits[i]) {
            return i;
        }
    }

    return m_digits.size();
}

// ------------------------------------------------------------
// Verify one arbitrary (p,m) case algebraically
// ------------------------------------------------------------

bool verify_algebraic(
    u64 p,
    u64 m
) {
    const auto d = base_digits(p, m);

    const u128 product = miss_count(p, m);

    u128 type_total = 0;
    u128 weighted_total = 0;

    for (std::size_t r = 0; r < d.size(); ++r) {
        const u128 c = type_count(d, r);
        const u128 g = gap_length(p, d, r);

        type_total += c;
        weighted_total += c * g;
    }

    // There are product-1 gaps between MISS elements.
    const bool count_ok =
        type_total == product - 1;

    // HIT total = m+1 - MISS total.
    const u128 expected_hits =
        (u128)m + 1 - product;

    const bool weighted_ok =
        weighted_total == expected_hits;

    return count_ok && weighted_ok;
}

// ------------------------------------------------------------
// Verify actual representatives for each carry type
// ------------------------------------------------------------

bool verify_type_samples(
    u64 p,
    u64 m
) {
    const auto d = base_digits(p, m);

    for (std::size_t r = 0; r < d.size(); ++r) {
        if (d[r] == 0) {
            continue;
        }

        // First representative of this type.
        {
            const u64 digit_r = 0;

            const u128 x =
                construct_type_element(
                    p, d, r, digit_r, false
                );

            if (x >= (u128)m) {
                return false;
            }

            const u128 y = mixed_successor(p, d, x);

            const u128 actual_gap = y - x - 1;
            const u128 expected_gap = gap_length(p, d, r);

            if (actual_gap != expected_gap) {
                return false;
            }

            if (carry_position(p, d, x) != r) {
                return false;
            }
        }

        // Last representative of this type.
        {
            const u64 digit_r = d[r] - 1;

            const u128 x =
                construct_type_element(
                    p, d, r, digit_r, true
                );

            if (x >= (u128)m) {
                return false;
            }

            const u128 y = mixed_successor(p, d, x);

            const u128 actual_gap = y - x - 1;
            const u128 expected_gap = gap_length(p, d, r);

            if (actual_gap != expected_gap) {
                return false;
            }

            if (carry_position(p, d, x) != r) {
                return false;
            }
        }
    }

    return true;
}

// ------------------------------------------------------------
// Exhaustive small-case verification.
//
// Enumerates every MISS element using the mixed-radix ordering
// and counts its actual carry type.
// ------------------------------------------------------------

bool verify_exhaustive_small(
    u64 p,
    u64 m
) {
    const auto d = base_digits(p, m);

    std::vector<u64> actual_counts(d.size(), 0);

    u128 x = 0;

    while (true) {
        const bool is_last = (x == (u128)m);

        if (!is_last) {
            const std::size_t r =
                carry_position(p, d, x);

            if (r >= d.size()) {
                return false;
            }

            ++actual_counts[r];
        }

        if (is_last) {
            break;
        }

        x = mixed_successor(p, d, x);
    }

    for (std::size_t r = 0; r < d.size(); ++r) {
        const u128 expected = type_count(d, r);

        if ((u128)actual_counts[r] != expected) {
            return false;
        }
    }

    return true;
}

// ------------------------------------------------------------
// One complete case
// ------------------------------------------------------------

bool verify_case(
    u64 p,
    u64 m,
    bool exhaustive
) {
    if (!verify_algebraic(p, m)) {
        return false;
    }

    if (!verify_type_samples(p, m)) {
        return false;
    }

    if (exhaustive) {
        if (!verify_exhaustive_small(p, m)) {
            return false;
        }
    }

    return true;
}

// ------------------------------------------------------------
// Main experiment
// ------------------------------------------------------------

int main() {
    std::cout << "START EXPERIMENT 257\n";
    std::cout << "ARBITRARY-p,m CARRY-TYPE DISTRIBUTION\n";
    std::cout << "AND WEIGHTED HIT-GAP IDENTITY\n\n";

    // --------------------------------------------------------
    // Deterministic cases
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
        const bool exhaustive = (m <= 50000);

        const bool ok =
            verify_case(p, m, exhaustive);

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
    // Random algebraic cases
    // --------------------------------------------------------

    std::mt19937_64 rng(0x257257257ULL);

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37, 41, 43, 47
    };

    const std::size_t random_cases = 100000;
    std::size_t random_fail = 0;

    for (std::size_t tc = 0; tc < random_cases; ++tc) {
        const u64 p =
            primes[rng() % primes.size()];

        const u64 m =
            rng() % 1000000000000000000ULL;

        if (!verify_case(p, m, false)) {
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
    // Exhaustive small random cases
    // --------------------------------------------------------

    const std::size_t exhaustive_cases = 500;
    std::size_t exhaustive_fail = 0;

    for (std::size_t tc = 0; tc < exhaustive_cases; ++tc) {
        const u64 p =
            primes[rng() % primes.size()];

        const u64 m =
            rng() % 50000;

        if (!verify_case(p, m, true)) {
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
    // Direct identity demonstration on one large case
    // --------------------------------------------------------

    {
        const u64 p = 13;
        const u64 m = 987654321012345678ULL;

        const auto d = base_digits(p, m);

        const u128 product = miss_count(p, m);

        u128 type_total = 0;
        u128 weighted_total = 0;

        std::cout << "LARGE_CASE\n";
        std::cout << "p=" << p << "\n";
        std::cout << "m=" << m << "\n";
        std::cout << "digits=" << d.size() << "\n";

        for (std::size_t r = 0; r < d.size(); ++r) {
            const u128 c = type_count(d, r);
            const u128 g = gap_length(p, d, r);

            type_total += c;
            weighted_total += c * g;

            std::cout
                << "r=" << r
                << " m_r=" << d[r]
                << " C_r=";
            print_u128(c);
            std::cout << " G_r=";
            print_u128(g);
            std::cout << '\n';
        }

        const u128 expected_gap_count =
            product - 1;

        const u128 expected_hits =
            (u128)m + 1 - product;

        std::cout << "type_total=";
        print_u128(type_total);
        std::cout << "\nexpected_gap_count=";
        print_u128(expected_gap_count);

        std::cout << "\nweighted_total=";
        print_u128(weighted_total);
        std::cout << "\nexpected_hits=";
        print_u128(expected_hits);

        std::cout
            << "\nlarge_case_pass="
            << ((type_total == expected_gap_count &&
                 weighted_total == expected_hits) ? 1 : 0)
            << "\n\n";
    }

    // --------------------------------------------------------
    // Final result
    // --------------------------------------------------------

    const bool overall =
        deterministic_fail == 0 &&
        random_fail == 0 &&
        exhaustive_fail == 0;

    std::cout
        << "OVERALL PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout << "FINISHED EXPERIMENT 257\n";

    return overall ? 0 : 1;
}
