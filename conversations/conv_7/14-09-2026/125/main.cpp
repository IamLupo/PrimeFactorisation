#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

static constexpr u64 PRIME_LIMIT = 5000;
static constexpr u64 K_LIMIT     = 2000;
static constexpr u64 M_MAX       = 32;

struct Witness
{
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
    bool valid = false;
};

struct Excursion
{
    std::vector<u64> states;
};

struct Stats
{
    u64 prime_count = 0;
    u64 total_points = 0;

    u64 trajectories = 0;

    u64 compressed_trajectory_count = 0;
    u64 total_compressed_states = 0;

    u64 trajectories_start_at_one = 0;
    u64 trajectories_not_start_at_one = 0;

    u64 trajectories_end_at_one = 0;
    u64 trajectories_not_end_at_one = 0;

    u64 excursions = 0;

    u64 excursions_return_to_one = 0;
    u64 excursions_not_return_to_one = 0;

    u64 palindromic_excursions = 0;
    u64 non_palindromic_excursions = 0;

    u64 excursions_with_internal_repeat = 0;
    u64 excursions_without_internal_repeat = 0;

    u64 non_one_cycles = 0;

    u64 max_compressed_trajectory_length = 0;
    u64 max_excursion_length = 0;
    u64 max_excursion_depth = 0;

    u64 total_excursion_states = 0;

    u64 trajectories_with_multiple_excursions = 0;

    bool have_first_start_counterexample = false;
    u64 first_start_counterexample_r = 0;
    std::vector<u64> first_start_counterexample;

    bool have_first_end_counterexample = false;
    u64 first_end_counterexample_r = 0;
    std::vector<u64> first_end_counterexample;

    bool have_first_nonreturn_counterexample = false;
    u64 first_nonreturn_r = 0;
    std::vector<u64> first_nonreturn_excursion;

    bool have_first_nonpal_counterexample = false;
    u64 first_nonpal_r = 0;
    std::vector<u64> first_nonpal_excursion;

    bool have_first_repeat_counterexample = false;
    u64 first_repeat_r = 0;
    std::vector<u64> first_repeat_excursion;

    bool have_first_cycle_counterexample = false;
    u64 first_cycle_r = 0;
    std::vector<u64> first_cycle_trajectory;

    bool have_first_multiple_excursion = false;
    u64 first_multiple_excursion_r = 0;
    std::vector<u64> first_multiple_excursion_trajectory;

    bool have_longest_trajectory = false;
    u64 longest_trajectory_r = 0;
    std::vector<u64> longest_trajectory;

    bool have_longest_excursion = false;
    u64 longest_excursion_r = 0;
    std::vector<u64> longest_excursion;
};

static bool is_prime(u64 n)
{
    if (n < 2)
        return false;

    if (n == 2)
        return true;

    if ((n & 1ULL) == 0)
        return false;

    for (u64 d = 3; d <= n / d; d += 2)
    {
        if (n % d == 0)
            return false;
    }

    return true;
}

static std::vector<u64> generate_primes(u64 limit)
{
    std::vector<u64> primes;

    for (u64 n = 2; n <= limit; ++n)
    {
        if (is_prime(n))
            primes.push_back(n);
    }

    return primes;
}

static bool make_witness(
    u64 r,
    u64 m,
    u64 k,
    int sign,
    Witness& out)
{
    const i128 value =
        static_cast<i128>(m) * r +
        static_cast<i128>(sign);

    if (value <= 0)
        return false;

    if (value % k != 0)
        return false;

    const i128 t = value / k;

    if (t <= 0 ||
        t > std::numeric_limits<u64>::max())
    {
        return false;
    }

    out.valid = true;
    out.m = m;
    out.k = k;
    out.t = static_cast<u64>(t);
    out.sign = sign;

    return true;
}

static void update_best_for_m(
    u64 r,
    u64 m,
    u64 k,
    Witness& best)
{
    Witness candidate;

    if (make_witness(r, m, k, +1, candidate))
    {
        if (!best.valid ||
            candidate.k > best.k ||
            (candidate.k == best.k &&
             candidate.t < best.t))
        {
            best = candidate;
        }
    }

    if (make_witness(r, m, k, -1, candidate))
    {
        if (!best.valid ||
            candidate.k > best.k ||
            (candidate.k == best.k &&
             candidate.t < best.t))
        {
            best = candidate;
        }
    }
}

static bool normalized_better(
    const Witness& a,
    const Witness& b)
{
    return
        static_cast<u128>(a.k) * b.m >
        static_cast<u128>(b.k) * a.m;
}

static bool normalized_equal(
    const Witness& a,
    const Witness& b)
{
    return
        static_cast<u128>(a.k) * b.m ==
        static_cast<u128>(b.k) * a.m;
}

static void consider_winner(
    const Witness& candidate,
    Witness& winner)
{
    if (!candidate.valid)
        return;

    if (!winner.valid)
    {
        winner = candidate;
        return;
    }

    if (normalized_better(candidate, winner))
    {
        winner = candidate;
        return;
    }

    if (normalized_equal(candidate, winner))
    {
        if (candidate.t < winner.t ||
            (candidate.t == winner.t &&
             candidate.m < winner.m))
        {
            winner = candidate;
        }
    }
}

static std::vector<u64> compress_trajectory(
    const std::vector<u64>& raw)
{
    std::vector<u64> result;

    for (u64 m : raw)
    {
        if (result.empty() ||
            result.back() != m)
        {
            result.push_back(m);
        }
    }

    return result;
}

static bool is_palindrome(
    const std::vector<u64>& values)
{
    if (values.size() < 2)
        return true;

    std::size_t left = 0;
    std::size_t right = values.size() - 1;

    while (left < right)
    {
        if (values[left] != values[right])
            return false;

        ++left;
        --right;
    }

    return true;
}

static bool has_internal_repeat(
    const std::vector<u64>& values)
{
    for (std::size_t i = 0;
         i < values.size();
         ++i)
    {
        for (std::size_t j = i + 1;
             j < values.size();
             ++j)
        {
            if (values[i] == values[j])
                return true;
        }
    }

    return false;
}

static bool contains_only_non_one(
    const std::vector<u64>& values)
{
    if (values.empty())
        return false;

    for (u64 m : values)
    {
        if (m == 1)
            return false;
    }

    return true;
}

static u64 max_state(
    const std::vector<u64>& values)
{
    u64 result = 0;

    for (u64 m : values)
    {
        if (m > result)
            result = m;
    }

    return result;
}

static void print_sequence(
    const char* name,
    const std::vector<u64>& values)
{
    std::cout << name << "=";

    if (values.empty())
    {
        std::cout << "EMPTY\n";
        return;
    }

    for (std::size_t i = 0;
         i < values.size();
         ++i)
    {
        if (i != 0)
            std::cout << ",";

        std::cout << values[i];
    }

    std::cout << '\n';
}

static void analyze_trajectory(
    u64 r,
    const std::vector<u64>& trajectory,
    Stats& stats)
{
    ++stats.trajectories;

    if (trajectory.empty())
        return;

    if (trajectory.front() == 1)
    {
        ++stats.trajectories_start_at_one;
    }
    else
    {
        ++stats.trajectories_not_start_at_one;

        if (!stats.have_first_start_counterexample)
        {
            stats.have_first_start_counterexample = true;
            stats.first_start_counterexample_r = r;
            stats.first_start_counterexample =
                trajectory;
        }
    }

    if (trajectory.back() == 1)
    {
        ++stats.trajectories_end_at_one;
    }
    else
    {
        ++stats.trajectories_not_end_at_one;

        if (!stats.have_first_end_counterexample)
        {
            stats.have_first_end_counterexample = true;
            stats.first_end_counterexample_r = r;
            stats.first_end_counterexample =
                trajectory;
        }
    }

    if (trajectory.size() >
        stats.max_compressed_trajectory_length)
    {
        stats.max_compressed_trajectory_length =
            static_cast<u64>(trajectory.size());

        stats.have_longest_trajectory = true;
        stats.longest_trajectory_r = r;
        stats.longest_trajectory = trajectory;
    }

    stats.total_compressed_states +=
        static_cast<u64>(trajectory.size());

    if (trajectory.front() == 1 &&
        trajectory.back() == 1)
    {
        /*
            Split

                1, ..., 1, ..., 1

            into maximal interior excursions.
        */

        std::size_t start = 1;

        while (start + 1 < trajectory.size())
        {
            std::size_t end = start;

            while (end + 1 < trajectory.size() &&
                   trajectory[end + 1] != 1)
            {
                ++end;
            }

            if (end >= start)
            {
                std::vector<u64> excursion(
                    trajectory.begin() +
                        static_cast<std::ptrdiff_t>(start),
                    trajectory.begin() +
                        static_cast<std::ptrdiff_t>(end + 1)
                );

                ++stats.excursions;

                stats.total_excursion_states +=
                    static_cast<u64>(
                        excursion.size()
                    );

                if (excursion.size() >
                    stats.max_excursion_length)
                {
                    stats.max_excursion_length =
                        static_cast<u64>(
                            excursion.size()
                        );

                    stats.have_longest_excursion = true;
                    stats.longest_excursion_r = r;
                    stats.longest_excursion =
                        excursion;
                }

                const u64 depth =
                    max_state(excursion);

                if (depth >
                    stats.max_excursion_depth)
                {
                    stats.max_excursion_depth =
                        depth;
                }

                if (trajectory[end + 1] == 1)
                {
                    ++stats.excursions_return_to_one;
                }
                else
                {
                    ++stats.excursions_not_return_to_one;

                    if (!stats.have_first_nonreturn_counterexample)
                    {
                        stats.have_first_nonreturn_counterexample =
                            true;
                        stats.first_nonreturn_r = r;
                        stats.first_nonreturn_excursion =
                            excursion;
                    }
                }

                if (is_palindrome(excursion))
                {
                    ++stats.palindromic_excursions;
                }
                else
                {
                    ++stats.non_palindromic_excursions;

                    if (!stats.have_first_nonpal_counterexample)
                    {
                        stats.have_first_nonpal_counterexample =
                            true;
                        stats.first_nonpal_r = r;
                        stats.first_nonpal_excursion =
                            excursion;
                    }
                }

                if (has_internal_repeat(excursion))
                {
                    ++stats.excursions_with_internal_repeat;

                    if (!stats.have_first_repeat_counterexample)
                    {
                        stats.have_first_repeat_counterexample =
                            true;
                        stats.first_repeat_r = r;
                        stats.first_repeat_excursion =
                            excursion;
                    }
                }
                else
                {
                    ++stats.excursions_without_internal_repeat;
                }

                if (contains_only_non_one(excursion))
                {
                    ++stats.non_one_cycles;

                    if (!stats.have_first_cycle_counterexample)
                    {
                        stats.have_first_cycle_counterexample =
                            true;
                        stats.first_cycle_r = r;
                        stats.first_cycle_trajectory =
                            trajectory;
                    }
                }
            }

            start = end + 2;
        }
    }

    /*
        Count multiple excursions.

        Each occurrence of a non-1 block corresponds to one
        excursion from and back to state 1.
    */

    u64 excursion_count_for_trajectory = 0;

    for (std::size_t i = 0;
         i < trajectory.size();)
    {
        if (trajectory[i] == 1)
        {
            ++i;
            continue;
        }

        ++excursion_count_for_trajectory;

        while (i < trajectory.size() &&
               trajectory[i] != 1)
        {
            ++i;
        }
    }

    if (excursion_count_for_trajectory > 1)
    {
        ++stats.trajectories_with_multiple_excursions;

        if (!stats.have_first_multiple_excursion)
        {
            stats.have_first_multiple_excursion = true;
            stats.first_multiple_excursion_r = r;
            stats.first_multiple_excursion_trajectory =
                trajectory;
        }
    }
}

int main()
{
    std::cout
        << "START EXPERIMENT 417\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    stats.prime_count =
        static_cast<u64>(primes.size());

    std::vector<Witness> best(
        M_MAX + 1
    );

    for (u64 r : primes)
    {
        for (u64 m = 0; m <= M_MAX; ++m)
            best[m] = Witness{};

        std::vector<u64> raw_trajectory;
        raw_trajectory.reserve(K_LIMIT);

        for (u64 K = 1;
             K <= K_LIMIT;
             ++K)
        {
            ++stats.total_points;

            for (u64 m = 1;
                 m <= M_MAX;
                 ++m)
            {
                update_best_for_m(
                    r,
                    m,
                    K,
                    best[m]
                );
            }

            Witness winner;

            for (u64 m = 1;
                 m <= M_MAX;
                 ++m)
            {
                consider_winner(
                    best[m],
                    winner
                );
            }

            if (winner.valid)
                raw_trajectory.push_back(
                    winner.m
                );
        }

        const std::vector<u64> trajectory =
            compress_trajectory(
                raw_trajectory
            );

        ++stats.compressed_trajectory_count;

        analyze_trajectory(
            r,
            trajectory,
            stats
        );
    }

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT
        << '\n'
        << "K_LIMIT="
        << K_LIMIT
        << '\n'
        << "M_MAX="
        << M_MAX
        << '\n'
        << "PRIME_COUNT="
        << stats.prime_count
        << '\n'
        << '\n';

    std::cout
        << "TOTAL_POINTS="
        << stats.total_points
        << '\n'
        << "TRAJECTORIES="
        << stats.trajectories
        << '\n'
        << "COMPRESSED_TRAJECTORY_COUNT="
        << stats.compressed_trajectory_count
        << '\n'
        << "TOTAL_COMPRESSED_STATES="
        << stats.total_compressed_states
        << '\n'
        << '\n';

    std::cout
        << "TRAJECTORIES_START_AT_ONE="
        << stats.trajectories_start_at_one
        << '\n'
        << "TRAJECTORIES_NOT_START_AT_ONE="
        << stats.trajectories_not_start_at_one
        << '\n'
        << "TRAJECTORIES_END_AT_ONE="
        << stats.trajectories_end_at_one
        << '\n'
        << "TRAJECTORIES_NOT_END_AT_ONE="
        << stats.trajectories_not_end_at_one
        << '\n'
        << '\n';

    std::cout
        << "EXCURSIONS="
        << stats.excursions
        << '\n'
        << "EXCURSIONS_RETURN_TO_ONE="
        << stats.excursions_return_to_one
        << '\n'
        << "EXCURSIONS_NOT_RETURN_TO_ONE="
        << stats.excursions_not_return_to_one
        << '\n'
        << '\n';

    std::cout
        << "PALINDROMIC_EXCURSIONS="
        << stats.palindromic_excursions
        << '\n'
        << "NON_PALINDROMIC_EXCURSIONS="
        << stats.non_palindromic_excursions
        << '\n'
        << '\n';

    std::cout
        << "EXCURSIONS_WITH_INTERNAL_REPEAT="
        << stats.excursions_with_internal_repeat
        << '\n'
        << "EXCURSIONS_WITHOUT_INTERNAL_REPEAT="
        << stats.excursions_without_internal_repeat
        << '\n'
        << '\n';

    std::cout
        << "NON_ONE_CYCLES="
        << stats.non_one_cycles
        << '\n'
        << "TRAJECTORIES_WITH_MULTIPLE_EXCURSIONS="
        << stats.trajectories_with_multiple_excursions
        << '\n'
        << '\n';

    std::cout
        << "MAX_COMPRESSED_TRAJECTORY_LENGTH="
        << stats.max_compressed_trajectory_length
        << '\n'
        << "MAX_EXCURSION_LENGTH="
        << stats.max_excursion_length
        << '\n'
        << "MAX_EXCURSION_DEPTH="
        << stats.max_excursion_depth
        << '\n'
        << '\n';

    if (stats.have_longest_trajectory)
    {
        std::cout
            << "LONGEST_TRAJECTORY\n"
            << "R="
            << stats.longest_trajectory_r
            << '\n';

        print_sequence(
            "TRAJECTORY",
            stats.longest_trajectory
        );

        std::cout << '\n';
    }

    if (stats.have_longest_excursion)
    {
        std::cout
            << "LONGEST_EXCURSION\n"
            << "R="
            << stats.longest_excursion_r
            << '\n';

        print_sequence(
            "EXCURSION",
            stats.longest_excursion
        );

        std::cout << '\n';
    }

    if (stats.have_first_start_counterexample)
    {
        std::cout
            << "FIRST_START_COUNTEREXAMPLE\n"
            << "R="
            << stats.first_start_counterexample_r
            << '\n';

        print_sequence(
            "TRAJECTORY",
            stats.first_start_counterexample
        );

        std::cout << '\n';
    }

    if (stats.have_first_end_counterexample)
    {
        std::cout
            << "FIRST_END_COUNTEREXAMPLE\n"
            << "R="
            << stats.first_end_counterexample_r
            << '\n';

        print_sequence(
            "TRAJECTORY",
            stats.first_end_counterexample
        );

        std::cout << '\n';
    }

    if (stats.have_first_nonreturn_counterexample)
    {
        std::cout
            << "FIRST_NONRETURN_EXCURSION\n"
            << "R="
            << stats.first_nonreturn_r
            << '\n';

        print_sequence(
            "EXCURSION",
            stats.first_nonreturn_excursion
        );

        std::cout << '\n';
    }

    if (stats.have_first_nonpal_counterexample)
    {
        std::cout
            << "FIRST_NONPALINDROMIC_EXCURSION\n"
            << "R="
            << stats.first_nonpal_r
            << '\n';

        print_sequence(
            "EXCURSION",
            stats.first_nonpal_excursion
        );

        std::cout << '\n';
    }

    if (stats.have_first_repeat_counterexample)
    {
        std::cout
            << "FIRST_INTERNAL_REPEAT_EXCURSION\n"
            << "R="
            << stats.first_repeat_r
            << '\n';

        print_sequence(
            "EXCURSION",
            stats.first_repeat_excursion
        );

        std::cout << '\n';
    }

    if (stats.have_first_cycle_counterexample)
    {
        std::cout
            << "FIRST_NON_ONE_CYCLE_TRAJECTORY\n"
            << "R="
            << stats.first_cycle_r
            << '\n';

        print_sequence(
            "TRAJECTORY",
            stats.first_cycle_trajectory
        );

        std::cout << '\n';
    }

    if (stats.have_first_multiple_excursion)
    {
        std::cout
            << "FIRST_MULTIPLE_EXCURSION_TRAJECTORY\n"
            << "R="
            << stats.first_multiple_excursion_r
            << '\n';

        print_sequence(
            "TRAJECTORY",
            stats.first_multiple_excursion_trajectory
        );

        std::cout << '\n';
    }

    std::cout
        << "START_AT_ONE_STATUS="
        << (
            stats.trajectories_not_start_at_one == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "END_AT_ONE_STATUS="
        << (
            stats.trajectories_not_end_at_one == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "ALL_EXCURSIONS_RETURN_STATUS="
        << (
            stats.excursions_not_return_to_one == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "PALINDROMIC_EXCURSION_STATUS="
        << (
            stats.non_palindromic_excursions == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "NO_INTERNAL_REPEAT_STATUS="
        << (
            stats.excursions_with_internal_repeat == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "NO_NON_ONE_CYCLE_STATUS="
        << (
            stats.non_one_cycles == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "SINGLE_EXCURSION_STATUS="
        << (
            stats.trajectories_with_multiple_excursions == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 417\n";

    return 0;
}
