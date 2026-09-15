#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

static constexpr int EXPERIMENT = 159;

// ============================================================
// Primality
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
// Lucas direct test
//
// true  => p | C(m,d)
// false => p does not divide C(m,d)
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
// Get base-p digits, least significant first.
// ============================================================

std::vector<u64> base_p_digits(
    u64 x,
    u64 p
)
{
    std::vector<u64> digits;

    while (x > 0)
    {
        digits.push_back(x % p);
        x /= p;
    }

    if (digits.empty())
        digits.push_back(0);

    return digits;
}

// ============================================================
// Generalized digitwise prediction
//
// d = m - t.
//
// Lucas says HIT iff some digit of d exceeds the corresponding
// digit of m.
//
// Equivalently:
//
//     MISS iff subtraction m-t causes no borrow.
//
// That is equivalent to:
//
//     every base-p digit of t <= corresponding digit of m.
//
// Therefore:
//
//     HIT iff some t_i > m_i.
// ============================================================

bool digitwise_hit(
    u64 p,
    u64 m,
    u64 t
)
{
    // t > m is outside the intended range.
    if (t > m)
        return false;

    while (m > 0 || t > 0)
    {
        u64 md = m % p;
        u64 td = t % p;

        if (td > md)
            return true;

        m /= p;
        t /= p;
    }

    return false;
}

// ============================================================
// Extract maximal intervals from predicate.
// ============================================================

std::vector<std::pair<u64, u64>> make_intervals(
    u64 begin,
    u64 end,
    const std::vector<bool>& values
)
{
    std::vector<std::pair<u64, u64>> intervals;

    bool active = false;
    u64 start = 0;

    for (u64 t = begin; t <= end; ++t)
    {
        bool hit =
            values[
                static_cast<std::size_t>(
                    t - begin
                )
            ];

        if (hit && !active)
        {
            active = true;
            start = t;
        }

        bool last = (t == end);

        if (active && (!hit || last))
        {
            u64 stop =
                (hit && last)
                    ? t
                    : t - 1;

            intervals.push_back(
                {start, stop}
            );

            active = false;
        }
    }

    return intervals;
}

// ============================================================
// Check predicted interval list against direct evaluation.
// ============================================================

bool interval_contains(
    const std::vector<std::pair<u64, u64>>& intervals,
    u64 x
)
{
    for (const auto& interval : intervals)
    {
        if (x >= interval.first &&
            x <= interval.second)
        {
            return true;
        }

        if (x < interval.first)
            return false;
    }

    return false;
}

// ============================================================
// Build predicted HIT set from digitwise rule.
//
// We deliberately build it by digits, not Lucas.
// For each t, we only test:
//
//     exists digit t_i > m_i
//
// The result is then compared with Lucas.
// ============================================================

std::vector<bool> build_digitwise_pattern(
    u64 p,
    u64 m
)
{
    std::vector<bool> result;

    result.reserve(
        static_cast<std::size_t>(m)
    );

    for (u64 t = 0; t < m; ++t)
    {
        result.push_back(
            digitwise_hit(
                p,
                m,
                t
            )
        );
    }

    return result;
}

// ============================================================
// Convert a digit pattern into intervals.
// ============================================================

std::vector<std::pair<u64, u64>> predicted_intervals(
    u64 p,
    u64 m
)
{
    auto pattern =
        build_digitwise_pattern(
            p,
            m
        );

    return make_intervals(
        0,
        m - 1,
        pattern
    );
}

// ============================================================
// Exhaustive test of one p,m.
//
// Tests every:
//
//     t = 0 ... m-1
//
// against both:
//
//     Lucas
//
// and:
//
//     digitwise condition.
// ============================================================

struct TestResult
{
    bool passed;
    u64 mismatch_t;
    std::size_t interval_count;
};

TestResult test_configuration(
    u64 p,
    u64 m
)
{
    bool passed = true;
    u64 mismatch_t = 0;

    std::vector<bool> pattern;

    pattern.reserve(
        static_cast<std::size_t>(m)
    );

    for (u64 t = 0; t < m; ++t)
    {
        u64 d = m - t;

        bool direct =
            lucas_hit(
                p,
                m,
                d
            );

        bool predicted =
            digitwise_hit(
                p,
                m,
                t
            );

        pattern.push_back(predicted);

        if (direct != predicted)
        {
            if (passed)
                mismatch_t = t;

            passed = false;
        }
    }

    auto intervals =
        make_intervals(
            0,
            m - 1,
            pattern
        );

    return
    {
        passed,
        mismatch_t,
        intervals.size()
    };
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
// Generate selected m values that cross the quotient boundary.
//
// Includes:
//     a < p
//     a = p
//     a > p
//     a much larger than p
// ============================================================

std::vector<u64> interesting_m(
    u64 p
)
{
    std::vector<u64> values;

    const u64 r_values[] =
    {
        0,
        1,
        p / 3,
        p / 2,
        p - 2,
        p - 1
    };

    const u64 a_values[] =
    {
        1,
        2,
        p - 1,
        p,
        p + 1,
        2 * p,
        2 * p + 1,
        3 * p + 2
    };

    for (u64 a : a_values)
    {
        for (u64 r : r_values)
        {
            if (r < p)
                values.push_back(
                    a * p + r
                );
        }
    }

    std::sort(
        values.begin(),
        values.end()
    );

    values.erase(
        std::unique(
            values.begin(),
            values.end()
        ),
        values.end()
    );

    return values;
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

    // --------------------------------------------------------
    // Phase 1:
    // Exhaustive small primes.
    // --------------------------------------------------------

    std::cout
        << "PHASE 1: EXHAUSTIVE SMALL PRIMES"
        << "\n";

    int total = 0;
    int passed = 0;
    int failed = 0;

    for (u64 p = 2; p <= 31; ++p)
    {
        if (!is_prime(p))
            continue;

        // Test m values through several base-p digits.
        for (u64 m = 1;
             m <= 4 * p * p;
             ++m)
        {
            ++total;

            TestResult result =
                test_configuration(
                    p,
                    m
                );

            if (result.passed)
            {
                ++passed;
            }
            else
            {
                ++failed;

                std::cout
                    << "MISMATCH"
                    << "\n";

                std::cout
                    << "  p="
                    << p
                    << "\n";

                std::cout
                    << "  m="
                    << m
                    << "\n";

                std::cout
                    << "  mismatch_t="
                    << result.mismatch_t
                    << "\n";
            }
        }
    }

    std::cout
        << "  total="
        << total
        << "\n";

    std::cout
        << "  passed="
        << passed
        << "\n";

    std::cout
        << "  failed="
        << failed
        << "\n\n";

    // --------------------------------------------------------
    // Phase 2:
    // Explicit quotient-overflow cases.
    //
    // These are exactly the cases where Experiment 158's
    // two-digit formula failed.
    // --------------------------------------------------------

    std::cout
        << "PHASE 2: QUOTIENT OVERFLOW CASES"
        << "\n";

    const u64 primes[] =
    {
        2, 3, 5, 7, 11, 13, 17, 19
    };

    for (u64 p : primes)
    {
        if (!is_prime(p))
            continue;

        auto m_values =
            interesting_m(p);

        for (u64 m : m_values)
        {
            TestResult result =
                test_configuration(
                    p,
                    m
                );

            std::cout
                << "CASE"
                << "\n";

            std::cout
                << "  p="
                << p
                << "\n";

            std::cout
                << "  m="
                << m
                << "\n";

            std::cout
                << "  quotient="
                << m / p
                << "\n";

            std::cout
                << "  remainder="
                << m % p
                << "\n";

            std::cout
                << "  base_p_digits_m=";

            auto digits =
                base_p_digits(
                    m,
                    p
                );

            for (std::size_t i = 0;
                 i < digits.size();
                 ++i)
            {
                if (i != 0)
                    std::cout << ",";

                std::cout
                    << digits[i];
            }

            std::cout
                << "\n";

            std::cout
                << "  formula_pass="
                << (result.passed ? 1 : 0)
                << "\n";

            std::cout
                << "  interval_count="
                << result.interval_count
                << "\n";

            auto intervals =
                predicted_intervals(
                    p,
                    m
                );

            std::cout
                << "  hit_intervals=";

            print_intervals(intervals);

            std::cout
                << "\n\n";
        }
    }

    // --------------------------------------------------------
    // Phase 3:
    // Larger random values.
    // --------------------------------------------------------

    std::cout
        << "PHASE 3: RANDOM LARGER VALUES"
        << "\n";

    std::mt19937_64 rng(
        0x159159159ULL
    );

    const int RANDOM_CASES = 20;

    for (int case_id = 1;
         case_id <= RANDOM_CASES;
         ++case_id)
    {
        std::uniform_int_distribution<u64>
            prime_dist(101, 1000);

        u64 p;

        do
        {
            p = prime_dist(rng);
        }
        while (!is_prime(p));

        std::uniform_int_distribution<u64>
            digit_count_dist(1, 4);

        int digits =
            static_cast<int>(
                digit_count_dist(rng)
            );

        u64 max_m = 1;

        for (int i = 0;
             i < digits;
             ++i)
        {
            if (max_m > 500000 / p)
            {
                max_m = 500000;
                break;
            }

            max_m *= p;
        }

        if (max_m < 10)
            max_m = 10;

        std::uniform_int_distribution<u64>
            m_dist(
                1,
                max_m
            );

        u64 m =
            m_dist(rng);

        TestResult result =
            test_configuration(
                p,
                m
            );

        auto intervals =
            predicted_intervals(
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
            << "  m="
            << m
            << "\n";

        std::cout
            << "  base_p_digit_count="
            << base_p_digits(m, p).size()
            << "\n";

        std::cout
            << "  formula_pass="
            << (result.passed ? 1 : 0)
            << "\n";

        std::cout
            << "  hit_interval_count="
            << intervals.size()
            << "\n";

        std::cout
            << "  hit_intervals=";

        print_intervals(intervals);

        std::cout
            << "\n\n";
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
