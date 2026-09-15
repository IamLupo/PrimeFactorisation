#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

static constexpr int EXPERIMENT = 158;

// ============================================================
// Simple primality test
// ============================================================

bool is_prime(u64 n)
{
    if (n < 2)
        return false;

    if (n % 2 == 0)
        return n == 2;

    for (u64 d = 3; d <= n / d; d += 2)
    {
        if (n % d == 0)
            return false;
    }

    return true;
}

// ============================================================
// Random prime in a range
// ============================================================

u64 random_prime(
    u64 low,
    u64 high,
    std::mt19937_64& rng
)
{
    std::uniform_int_distribution<u64> dist(
        low,
        high
    );

    for (;;)
    {
        u64 x = dist(rng);

        if ((x & 1ULL) == 0)
            ++x;

        if (x > high)
            continue;

        if (is_prime(x))
            return x;
    }
}

// ============================================================
// Direct Lucas divisibility test
//
// Returns true iff
//
//     p | C(m,d)
//
// by Lucas' theorem.
// ============================================================

bool lucas_hit(
    u64 p,
    u64 m,
    u64 d
)
{
    while (m > 0 || d > 0)
    {
        u64 md = m % p;
        u64 dd = d % p;

        if (dd > md)
            return true;

        m /= p;
        d /= p;
    }

    return false;
}

// ============================================================
// Kummer carry count
//
// p | C(m,d) iff at least one carry occurs.
// ============================================================

int kummer_carries(
    u64 p,
    u64 m,
    u64 d
)
{
    u64 x = d;
    u64 y = m - d;

    u64 carry = 0;
    int count = 0;

    while (x > 0 || y > 0 || carry > 0)
    {
        u64 xd = x % p;
        u64 yd = y % p;

        u64 sum = xd + yd + carry;

        if (sum >= p)
        {
            ++count;
            carry = 1;
        }
        else
        {
            carry = 0;
        }

        x /= p;
        y /= p;
    }

    return count;
}

// ============================================================
// Proposed general hit predicate
//
// m = a*p + r
//
// HIT intervals:
//
//     [j*p + r + 1, (j+1)*p - 1]
//
// for j = 0 ... a-1.
//
// ============================================================

bool predicted_hit(
    u64 p,
    u64 a,
    u64 r,
    u64 t
)
{
    for (u64 j = 0; j < a; ++j)
    {
        u64 low =
            j * p + r + 1;

        u64 high =
            (j + 1) * p - 1;

        if (low <= high &&
            t >= low &&
            t <= high)
        {
            return true;
        }
    }

    return false;
}

// ============================================================
// Extract contiguous intervals
// ============================================================

std::vector<std::pair<u64, u64>> make_intervals(
    const std::vector<bool>& values,
    u64 offset
)
{
    std::vector<std::pair<u64, u64>> result;

    bool active = false;
    u64 start = 0;

    for (std::size_t i = 0;
         i < values.size();
         ++i)
    {
        u64 t =
            offset +
            static_cast<u64>(i);

        bool h = values[i];

        if (h && !active)
        {
            active = true;
            start = t;
        }

        bool last =
            (i + 1 == values.size());

        if (active && (!h || last))
        {
            u64 end =
                (h && last)
                    ? t
                    : t - 1;

            result.push_back(
                {start, end}
            );

            active = false;
        }
    }

    return result;
}

// ============================================================
// Print intervals
// ============================================================

void print_intervals(
    const std::vector<std::pair<u64, u64>>& intervals
)
{
    if (intervals.empty())
    {
        std::cout << "NONE";
        return;
    }

    for (std::size_t i = 0;
         i < intervals.size();
         ++i)
    {
        if (i != 0)
            std::cout << ",";

        std::cout
            << "["
            << intervals[i].first
            << ","
            << intervals[i].second
            << "]";
    }
}

// ============================================================
// One test
// ============================================================

struct TestResult
{
    bool lucas_match;
    bool kummer_match;
    u64 mismatch_t;
    int hit_intervals;
};

// ============================================================
// Test one p,a,r combination
// ============================================================

TestResult test_configuration(
    u64 p,
    u64 a,
    u64 r
)
{
    const u64 m =
        a * p + r;

    bool lucas_match = true;
    bool kummer_match = true;

    u64 mismatch_t = 0;

    std::vector<bool> actual;
    actual.reserve(
        static_cast<std::size_t>(m)
    );

    // --------------------------------------------------------
    // Test every t = 0,...,m-1.
    // --------------------------------------------------------

    for (u64 t = 0; t < m; ++t)
    {
        u64 d =
            m - t;

        bool lucas =
            lucas_hit(
                p,
                m,
                d
            );

        bool predicted =
            predicted_hit(
                p,
                a,
                r,
                t
            );

        int carries =
            kummer_carries(
                p,
                m,
                d
            );

        bool kummer =
            carries > 0;

        actual.push_back(lucas);

        if (lucas != predicted)
        {
            if (lucas_match)
                mismatch_t = t;

            lucas_match = false;
        }

        if (lucas != kummer)
            kummer_match = false;
    }

    auto intervals =
        make_intervals(
            actual,
            0
        );

    return
    {
        lucas_match,
        kummer_match,
        mismatch_t,
        static_cast<int>(
            intervals.size()
        )
    };
}

// ============================================================
// Find actual hit intervals for display
// ============================================================

std::vector<std::pair<u64, u64>> actual_intervals(
    u64 p,
    u64 m
)
{
    std::vector<bool> values;

    values.reserve(
        static_cast<std::size_t>(m)
    );

    for (u64 t = 0; t < m; ++t)
    {
        u64 d = m - t;

        values.push_back(
            lucas_hit(
                p,
                m,
                d
            )
        );
    }

    return make_intervals(
        values,
        0
    );
}

// ============================================================
// Main
// ============================================================

int main()
{
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n\n";

    std::mt19937_64 rng(
        0x158158158ULL
    );

    int total_tests = 0;
    int lucas_pass = 0;
    int kummer_pass = 0;
    int complete_pass = 0;

    // --------------------------------------------------------
    // Small exhaustive phase.
    //
    // Test every prime p <= 100,
    // every a = 1...8,
    // every r = 0...p-1.
    // --------------------------------------------------------

    std::cout
        << "PHASE 1: EXHAUSTIVE SMALL PRIMES"
        << "\n";

    for (u64 p = 2; p <= 100; ++p)
    {
        if (!is_prime(p))
            continue;

        for (u64 a = 1; a <= 8; ++a)
        {
            for (u64 r = 0; r < p; ++r)
            {
                ++total_tests;

                TestResult result =
                    test_configuration(
                        p,
                        a,
                        r
                    );

                if (result.lucas_match)
                    ++lucas_pass;

                if (result.kummer_match)
                    ++kummer_pass;

                if (result.lucas_match &&
                    result.kummer_match)
                {
                    ++complete_pass;
                }
                else
                {
                    std::cout
                        << "MISMATCH"
                        << "\n";

                    std::cout
                        << "  p="
                        << p
                        << "\n";

                    std::cout
                        << "  a="
                        << a
                        << "\n";

                    std::cout
                        << "  r="
                        << r
                        << "\n";

                    std::cout
                        << "  m="
                        << a * p + r
                        << "\n";

                    std::cout
                        << "  mismatch_t="
                        << result.mismatch_t
                        << "\n";

                    std::cout
                        << "\n";
                }
            }
        }
    }

    std::cout
        << "  total_tests="
        << total_tests
        << "\n";

    std::cout
        << "  lucas_formula_pass="
        << lucas_pass
        << "\n";

    std::cout
        << "  kummers_theorem_pass="
        << kummer_pass
        << "\n";

    std::cout
        << "  complete_pass="
        << complete_pass
        << "\n\n";

    // --------------------------------------------------------
    // Random larger examples.
    // --------------------------------------------------------

    std::cout
        << "PHASE 2: RANDOM LARGER CONFIGURATIONS"
        << "\n";

    const int RANDOM_CASES = 20;

    for (int case_id = 1;
         case_id <= RANDOM_CASES;
         ++case_id)
    {
        u64 p =
            random_prime(
                1000,
                10000,
                rng
            );

        std::uniform_int_distribution<u64>
            a_dist(1, 12);

        u64 a =
            a_dist(rng);

        std::uniform_int_distribution<u64>
            r_dist(0, p - 1);

        u64 r =
            r_dist(rng);

        u64 m =
            a * p + r;

        TestResult result =
            test_configuration(
                p,
                a,
                r
            );

        auto intervals =
            actual_intervals(
                p,
                m
            );

        std::cout
            << "CASE "
            << case_id
            << "\n";

        std::cout
            << "  p="
            << p
            << "\n";

        std::cout
            << "  a="
            << a
            << "\n";

        std::cout
            << "  r="
            << r
            << "\n";

        std::cout
            << "  m="
            << m
            << "\n";

        std::cout
            << "  lucas_formula_match="
            << (result.lucas_match ? 1 : 0)
            << "\n";

        std::cout
            << "  k_lucas_consistent="
            << (result.kummer_match ? 1 : 0)
            << "\n";

        std::cout
            << "  actual_hit_interval_count="
            << intervals.size()
            << "\n";

        std::cout
            << "  actual_hit_intervals=";

        print_intervals(intervals);

        std::cout
            << "\n";

        std::cout
            << "  predicted_hit_interval_count="
            << (
                (r + 1 < p)
                    ? a
                    : 0
            )
            << "\n";

        std::cout
            << "\n";
    }

    // --------------------------------------------------------
    // Explicit example corresponding to Experiment 157 style.
    // --------------------------------------------------------

    std::cout
        << "PHASE 3: EXPLICIT QUOTIENT-THREE EXAMPLE"
        << "\n";

    {
        const u64 p = 25229;
        const u64 m = 90873;

        const u64 a =
            m / p;

        const u64 r =
            m % p;

        auto intervals =
            actual_intervals(
                p,
                m
            );

        std::cout
            << "  p="
            << p
            << "\n";

        std::cout
            << "  m="
            << m
            << "\n";

        std::cout
            << "  a="
            << a
            << "\n";

        std::cout
            << "  r="
            << r
            << "\n";

        std::cout
            << "  m_equals_ap_plus_r="
            << (
                m == a * p + r
                    ? 1
                    : 0
            )
            << "\n";

        std::cout
            << "  actual_hit_intervals=";

        print_intervals(intervals);

        std::cout
            << "\n";

        std::cout
            << "  predicted_hit_intervals=";

        bool first = true;

        for (u64 j = 0; j < a; ++j)
        {
            u64 low =
                j * p + r + 1;

            u64 high =
                (j + 1) * p - 1;

            if (low > high)
                continue;

            if (!first)
                std::cout << ",";

            std::cout
                << "["
                << low
                << ","
                << high
                << "]";

            first = false;
        }

        if (first)
            std::cout << "NONE";

        std::cout
            << "\n";
    }

    std::cout
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
