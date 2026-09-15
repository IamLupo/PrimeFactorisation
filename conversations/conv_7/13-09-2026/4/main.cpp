#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

static constexpr int EXPERIMENT = 160;

// ============================================================
// Prime test
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
// General Lucas HIT predicate
//
// We test:
//
//     p | C(m,m-t)
//
// via Lucas.
//
// Equivalent:
//
//     HIT iff some base-p digit of t
//     exceeds the corresponding digit of m.
// ============================================================

bool digitwise_hit(
    u64 p,
    u64 m,
    u64 t
)
{
    if (t >= m)
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
// Build maximal HIT intervals
// ============================================================

std::vector<std::pair<u64, u64>>
build_hit_intervals(
    u64 p,
    u64 m
)
{
    std::vector<std::pair<u64, u64>> intervals;

    bool active = false;
    u64 start = 0;

    for (u64 t = 0; t < m; ++t)
    {
        bool hit =
            digitwise_hit(
                p,
                m,
                t
            );

        if (hit && !active)
        {
            active = true;
            start = t;
        }

        bool last =
            (t + 1 == m);

        if (active && (!hit || last))
        {
            u64 end =
                (hit && last)
                    ? t
                    : t - 1;

            intervals.push_back(
                {start, end}
            );

            active = false;
        }
    }

    return intervals;
}

// ============================================================
// GCD of a sequence of differences
//
// For values:
//     x0,x1,x2,...
//
// return:
//     gcd(|xi-x0|).
// ============================================================

u64 gcd_of_differences(
    const std::vector<u64>& values
)
{
    if (values.size() < 2)
        return 0;

    u64 g = 0;
    u64 first = values.front();

    for (std::size_t i = 1;
         i < values.size();
         ++i)
    {
        u64 a = values[i];
        u64 b = first;

        u64 diff =
            (a >= b)
                ? (a - b)
                : (b - a);

        g =
            std::gcd(
                g,
                diff
            );
    }

    return g;
}

// ============================================================
// GCD of arbitrary adjacent differences
// ============================================================

u64 gcd_of_adjacent_differences(
    const std::vector<u64>& values
)
{
    if (values.size() < 2)
        return 0;

    u64 g = 0;

    for (std::size_t i = 1;
         i < values.size();
         ++i)
    {
        u64 diff =
            values[i] -
            values[i - 1];

        g =
            std::gcd(
                g,
                diff
            );
    }

    return g;
}

// ============================================================
// Extract starts
// ============================================================

std::vector<u64> interval_starts(
    const std::vector<std::pair<u64, u64>>& intervals
)
{
    std::vector<u64> result;
    result.reserve(intervals.size());

    for (const auto& interval : intervals)
        result.push_back(interval.first);

    return result;
}

// ============================================================
// Extract ends
// ============================================================

std::vector<u64> interval_ends(
    const std::vector<std::pair<u64, u64>>& intervals
)
{
    std::vector<u64> result;
    result.reserve(intervals.size());

    for (const auto& interval : intervals)
        result.push_back(interval.second);

    return result;
}

// ============================================================
// All interval widths
// ============================================================

std::vector<u64> interval_widths(
    const std::vector<std::pair<u64, u64>>& intervals
)
{
    std::vector<u64> result;
    result.reserve(intervals.size());

    for (const auto& interval : intervals)
    {
        result.push_back(
            interval.second -
            interval.first +
            1
        );
    }

    return result;
}

// ============================================================
// GCD of interval widths
// ============================================================

u64 gcd_vector(
    const std::vector<u64>& values
)
{
    if (values.empty())
        return 0;

    u64 g = values.front();

    for (std::size_t i = 1;
         i < values.size();
         ++i)
    {
        g =
            std::gcd(
                g,
                values[i]
            );
    }

    return g;
}

// ============================================================
// A richer recovery attempt
//
// Candidate sources:
//
//   1. gcd of starts relative to first start
//   2. gcd of adjacent start differences
//   3. gcd of ends relative to first end
//   4. gcd of adjacent end differences
//   5. gcd of start/end differences combined
//
// ============================================================

struct RecoveryResult
{
    u64 start_gcd;
    u64 start_adjacent_gcd;

    u64 end_gcd;
    u64 end_adjacent_gcd;

    u64 combined_gcd;

    bool start_exact;
    bool start_adjacent_exact;

    bool end_exact;
    bool end_adjacent_exact;

    bool combined_exact;
};

// ============================================================
// Recover candidate p
// ============================================================

RecoveryResult recover_from_intervals(
    u64 p,
    const std::vector<std::pair<u64, u64>>& intervals
)
{
    auto starts =
        interval_starts(intervals);

    auto ends =
        interval_ends(intervals);

    u64 start_gcd =
        gcd_of_differences(
            starts
        );

    u64 start_adjacent_gcd =
        gcd_of_adjacent_differences(
            starts
        );

    u64 end_gcd =
        gcd_of_differences(
            ends
        );

    u64 end_adjacent_gcd =
        gcd_of_adjacent_differences(
            ends
        );

    std::vector<u64> all_boundaries;

    for (const auto& interval : intervals)
    {
        all_boundaries.push_back(
            interval.first
        );

        all_boundaries.push_back(
            interval.second
        );
    }

    u64 combined_gcd =
        gcd_of_differences(
            all_boundaries
        );

    return
    {
        start_gcd,
        start_adjacent_gcd,

        end_gcd,
        end_adjacent_gcd,

        combined_gcd,

        start_gcd == p,
        start_adjacent_gcd == p,

        end_gcd == p,
        end_adjacent_gcd == p,

        combined_gcd == p
    };
}

// ============================================================
// Print intervals
// ============================================================

void print_intervals(
    const std::vector<std::pair<u64, u64>>& intervals,
    std::size_t maximum = 20
)
{
    if (intervals.empty())
    {
        std::cout << "NONE";
        return;
    }

    std::size_t count =
        std::min(
            maximum,
            intervals.size()
        );

    for (std::size_t i = 0;
         i < count;
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

    if (count < intervals.size())
    {
        std::cout
            << ",..."
            << "("
            << intervals.size()
            << " total)";
    }
}

// ============================================================
// One configuration
// ============================================================

void run_case(
    u64 p,
    u64 m
)
{
    auto intervals =
        build_hit_intervals(
            p,
            m
        );

    RecoveryResult recovery =
        recover_from_intervals(
            p,
            intervals
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
        << "  interval_count="
        << intervals.size()
        << "\n";

    std::cout
        << "  hit_intervals=";

    print_intervals(
        intervals
    );

    std::cout
        << "\n";

    std::cout
        << "  start_gcd="
        << recovery.start_gcd
        << "\n";

    std::cout
        << "  start_adjacent_gcd="
        << recovery.start_adjacent_gcd
        << "\n";

    std::cout
        << "  end_gcd="
        << recovery.end_gcd
        << "\n";

    std::cout
        << "  end_adjacent_gcd="
        << recovery.end_adjacent_gcd
        << "\n";

    std::cout
        << "  combined_gcd="
        << recovery.combined_gcd
        << "\n";

    std::cout
        << "  start_exact="
        << (recovery.start_exact ? 1 : 0)
        << "\n";

    std::cout
        << "  start_adjacent_exact="
        << (recovery.start_adjacent_exact ? 1 : 0)
        << "\n";

    std::cout
        << "  end_exact="
        << (recovery.end_exact ? 1 : 0)
        << "\n";

    std::cout
        << "  end_adjacent_exact="
        << (recovery.end_adjacent_exact ? 1 : 0)
        << "\n";

    std::cout
        << "  combined_exact="
        << (recovery.combined_exact ? 1 : 0)
        << "\n\n";
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
    // Exhaustive small prime / m sweep.
    // --------------------------------------------------------

    std::cout
        << "PHASE 1: EXHAUSTIVE SMALL CASES"
        << "\n";

    int total_multi = 0;

    int start_exact = 0;
    int start_adjacent_exact = 0;

    int end_exact = 0;
    int end_adjacent_exact = 0;

    int combined_exact = 0;

    int zero_intervals = 0;
    int one_interval = 0;

    for (u64 p = 2; p <= 47; ++p)
    {
        if (!is_prime(p))
            continue;

        for (u64 m = 2;
             m <= 6 * p * p;
             ++m)
        {
            auto intervals =
                build_hit_intervals(
                    p,
                    m
                );

            if (intervals.empty())
            {
                ++zero_intervals;
                continue;
            }

            if (intervals.size() == 1)
            {
                ++one_interval;
                continue;
            }

            ++total_multi;

            RecoveryResult result =
                recover_from_intervals(
                    p,
                    intervals
                );

            if (result.start_exact)
                ++start_exact;

            if (result.start_adjacent_exact)
                ++start_adjacent_exact;

            if (result.end_exact)
                ++end_exact;

            if (result.end_adjacent_exact)
                ++end_adjacent_exact;

            if (result.combined_exact)
                ++combined_exact;
        }
    }

    std::cout
        << "  multi_interval_cases="
        << total_multi
        << "\n";

    std::cout
        << "  zero_interval_cases="
        << zero_intervals
        << "\n";

    std::cout
        << "  one_interval_cases="
        << one_interval
        << "\n";

    std::cout
        << "  start_exact="
        << start_exact
        << "/"
        << total_multi
        << "\n";

    std::cout
        << "  start_adjacent_exact="
        << start_adjacent_exact
        << "/"
        << total_multi
        << "\n";

    std::cout
        << "  end_exact="
        << end_exact
        << "/"
        << total_multi
        << "\n";

    std::cout
        << "  end_adjacent_exact="
        << end_adjacent_exact
        << "/"
        << total_multi
        << "\n";

    std::cout
        << "  combined_exact="
        << combined_exact
        << "/"
        << total_multi
        << "\n\n";

    // --------------------------------------------------------
    // Phase 2:
    // Selected examples from the previous experiments.
    // --------------------------------------------------------

    std::cout
        << "PHASE 2: SELECTED EXAMPLES"
        << "\n";

    run_case(
        25229,
        90873
    );

    run_case(
        43103,
        97721
    );

    run_case(
        29873,
        99073
    );

    run_case(
        5,
        59
    );

    run_case(
        11,
        395
    );

    run_case(
        17,
        917
    );

    // --------------------------------------------------------
    // Phase 3:
    // Larger random cases.
    // --------------------------------------------------------

    std::cout
        << "PHASE 3: RANDOM LARGER CASES"
        << "\n";

    std::mt19937_64 rng(
        0x160160160ULL
    );

    const int RANDOM_CASES = 20;

    for (int case_id = 1;
         case_id <= RANDOM_CASES;
         ++case_id)
    {
        std::uniform_int_distribution<u64>
            prime_dist(
                101,
                5000
            );

        u64 p;

        do
        {
            p = prime_dist(rng);
        }
        while (!is_prime(p));

        std::uniform_int_distribution<u64>
            a_dist(
                1,
                4 * p
            );

        u64 a =
            a_dist(rng);

        std::uniform_int_distribution<u64>
            r_dist(
                0,
                p - 1
            );

        u64 r =
            r_dist(rng);

        u64 m =
            a * p + r;

        std::cout
            << "RANDOM CASE "
            << case_id
            << "\n";

        run_case(
            p,
            m
        );
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
