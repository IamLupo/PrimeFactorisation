#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

static constexpr u64 PRIME_LIMIT = 5000;
static constexpr u64 K_LIMIT = 2000;
static constexpr u64 M_MAX = 32;

struct SignedState
{
    u64 m = 0;
    int sign = 0;
};

struct Witness
{
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
    bool valid = false;
};

struct Transition
{
    u64 K = 0;
    Witness old_w;
    Witness new_w;
};

struct Trajectory
{
    std::vector<Witness> winners;
    std::vector<Transition> transitions;
};

struct Excursion
{
    std::vector<Witness> states;
    std::vector<Transition> transitions;
};

struct Stats
{
    u64 prime_count = 0;
    u64 total_points = 0;

    u64 trajectories = 0;

    u64 projected_excursions = 0;
    u64 signed_excursions = 0;

    u64 projected_palindromes = 0;
    u64 projected_nonpalindromes = 0;

    u64 signed_palindromes = 0;
    u64 signed_nonpalindromes = 0;

    u64 signed_state_reverse_matches = 0;
    u64 signed_state_reverse_failures = 0;

    u64 signed_transition_reverse_matches = 0;
    u64 signed_transition_reverse_failures = 0;

    u64 signed_divisor_reverse_matches = 0;
    u64 signed_divisor_reverse_failures = 0;

    u64 signed_K_reverse_matches = 0;
    u64 signed_K_reverse_failures = 0;

    u64 signed_k_reverse_matches = 0;
    u64 signed_k_reverse_failures = 0;

    u64 signed_t_reverse_matches = 0;
    u64 signed_t_reverse_failures = 0;

    u64 signed_delta_reverse_matches = 0;
    u64 signed_delta_reverse_failures = 0;

    u64 signed_correction_reverse_matches = 0;
    u64 signed_correction_reverse_failures = 0;

    u64 projected_only_palindrome = 0;
    u64 signed_only_palindrome = 0;

    u64 signed_internal_repeats = 0;
    u64 signed_no_internal_repeats = 0;

    u64 signed_multiple_excursions = 0;

    u64 max_signed_trajectory_length = 0;
    u64 max_signed_excursion_length = 0;
    u64 max_signed_depth = 0;

    u64 max_signed_state_m = 0;

    bool have_first_signed_nonpalindrome = false;
    u64 first_signed_nonpalindrome_r = 0;
    Excursion first_signed_nonpalindrome;

    bool have_first_signed_reverse_failure = false;
    u64 first_signed_reverse_failure_r = 0;
    Excursion first_signed_reverse_failure;

    bool have_first_signed_divisor_failure = false;
    u64 first_signed_divisor_failure_r = 0;
    Excursion first_signed_divisor_failure;

    bool have_first_signed_cycle = false;
    u64 first_signed_cycle_r = 0;
    Excursion first_signed_cycle;
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
        t > static_cast<i128>(
            std::numeric_limits<u64>::max()))
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

static u64 witness_value(
    u64 r,
    const Witness& w)
{
    const i128 value =
        static_cast<i128>(w.m) * r +
        static_cast<i128>(w.sign);

    return static_cast<u64>(value);
}

static u64 transition_delta(
    const Transition& tr)
{
    const i128 delta =
        static_cast<i128>(tr.old_w.m) *
            tr.new_w.k -
        static_cast<i128>(tr.new_w.m) *
            tr.old_w.k;

    if (delta <= 0)
        return 0;

    return static_cast<u64>(delta);
}

static i128 transition_correction(
    const Transition& tr)
{
    return
        static_cast<i128>(tr.old_w.sign) *
            tr.new_w.k -
        static_cast<i128>(tr.new_w.sign) *
            tr.old_w.k;
}

static std::vector<Witness> get_projected_states(
    const Excursion& excursion)
{
    std::vector<Witness> result =
        excursion.states;

    return result;
}

static bool projected_palindrome(
    const Excursion& excursion)
{
    std::size_t left = 0;
    std::size_t right =
        excursion.states.size();

    if (right == 0)
        return true;

    --right;

    while (left < right)
    {
        if (excursion.states[left].m !=
            excursion.states[right].m)
        {
            return false;
        }

        ++left;
        --right;
    }

    return true;
}

static bool signed_palindrome(
    const Excursion& excursion)
{
    std::size_t left = 0;
    std::size_t right =
        excursion.states.size();

    if (right == 0)
        return true;

    --right;

    while (left < right)
    {
        const Witness& a =
            excursion.states[left];

        const Witness& b =
            excursion.states[right];

        if (a.m != b.m ||
            a.sign != b.sign)
        {
            return false;
        }

        ++left;
        --right;
    }

    return true;
}

static bool has_signed_internal_repeat(
    const Excursion& excursion)
{
    for (std::size_t i = 1;
         i + 1 < excursion.states.size();
         ++i)
    {
        for (std::size_t j = i + 1;
             j + 1 < excursion.states.size();
             ++j)
        {
            if (excursion.states[i].m ==
                    excursion.states[j].m &&
                excursion.states[i].sign ==
                    excursion.states[j].sign)
            {
                return true;
            }
        }
    }

    return false;
}

static u64 max_signed_m(
    const Excursion& excursion)
{
    u64 result = 0;

    for (const Witness& w : excursion.states)
    {
        if (w.m > result)
            result = w.m;
    }

    return result;
}

static bool same_signed_transition_reverse(
    const Transition& outward,
    const Transition& inward)
{
    return
        outward.old_w.m ==
            inward.new_w.m &&
        outward.old_w.sign ==
            inward.new_w.sign &&
        outward.new_w.m ==
            inward.old_w.m &&
        outward.new_w.sign ==
            inward.old_w.sign;
}

static bool same_divisor_pair_reverse(
    u64 r,
    const Transition& outward,
    const Transition& inward)
{
    const u64 outward_old =
        witness_value(
            r,
            outward.old_w
        );

    const u64 outward_new =
        witness_value(
            r,
            outward.new_w
        );

    const u64 inward_old =
        witness_value(
            r,
            inward.old_w
        );

    const u64 inward_new =
        witness_value(
            r,
            inward.new_w
        );

    return
        outward_old == inward_new &&
        outward_new == inward_old;
}

static bool signed_k_reverse(
    const Transition& outward,
    const Transition& inward)
{
    return
        outward.old_w.k ==
            inward.new_w.k &&
        outward.new_w.k ==
            inward.old_w.k;
}

static bool signed_t_reverse(
    const Transition& outward,
    const Transition& inward)
{
    return
        outward.old_w.t ==
            inward.new_w.t &&
        outward.new_w.t ==
            inward.old_w.t;
}

static bool same_delta_reverse(
    const Transition& outward,
    const Transition& inward)
{
    return
        transition_delta(outward) ==
        transition_delta(inward);
}

static bool same_correction_reverse(
    const Transition& outward,
    const Transition& inward)
{
    return
        transition_correction(outward) ==
        transition_correction(inward);
}

static bool same_K_reverse(
    const Transition& outward,
    const Transition& inward)
{
    return outward.K == inward.K;
}

static void print_witness(
    const char* name,
    const Witness& w)
{
    std::cout
        << name << ".m="
        << w.m << '\n'
        << name << ".k="
        << w.k << '\n'
        << name << ".t="
        << w.t << '\n'
        << name << ".sign="
        << w.sign << '\n';
}

static void print_excursion(
    const char* prefix,
    const Excursion& excursion)
{
    std::cout
        << prefix
        << ".STATES=";

    for (std::size_t i = 0;
         i < excursion.states.size();
         ++i)
    {
        if (i != 0)
            std::cout << ",";

        std::cout
            << "("
            << excursion.states[i].m
            << ","
            << excursion.states[i].sign
            << ")";
    }

    std::cout << '\n';

    std::cout
        << prefix
        << ".M_SEQUENCE=";

    for (std::size_t i = 0;
         i < excursion.states.size();
         ++i)
    {
        if (i != 0)
            std::cout << ",";

        std::cout
            << excursion.states[i].m;
    }

    std::cout << '\n';

    for (std::size_t i = 0;
         i < excursion.transitions.size();
         ++i)
    {
        const Transition& tr =
            excursion.transitions[i];

        std::cout
            << prefix
            << ".EVENT_"
            << i
            << ".K="
            << tr.K
            << '\n'
            << prefix
            << ".EVENT_"
            << i
            << ".OLD_VALUE=";

        std::cout
            << witness_value(
                excursion.transitions.empty()
                    ? 0
                    : 0,
                tr.old_w
            );

        std::cout << '\n';

        print_witness(
            "OLD",
            tr.old_w
        );

        print_witness(
            "NEW",
            tr.new_w
        );
    }
}

static void analyze_excursion(
    u64 r,
    const Excursion& excursion,
    Stats& stats)
{
    if (excursion.states.size() < 3)
        return;

    const bool projected =
        projected_palindrome(
            excursion
        );

    const bool signed_pal =
        signed_palindrome(
            excursion
        );

    if (projected)
        ++stats.projected_palindromes;
    else
        ++stats.projected_nonpalindromes;

    if (signed_pal)
        ++stats.signed_palindromes;
    else
    {
        ++stats.signed_nonpalindromes;

        if (!stats.have_first_signed_nonpalindrome)
        {
            stats.have_first_signed_nonpalindrome =
                true;

            stats.first_signed_nonpalindrome_r = r;
            stats.first_signed_nonpalindrome =
                excursion;
        }
    }

    if (projected && !signed_pal)
        ++stats.projected_only_palindrome;

    if (!projected && signed_pal)
        ++stats.signed_only_palindrome;

    if (has_signed_internal_repeat(
            excursion))
    {
        ++stats.signed_internal_repeats;
    }
    else
    {
        ++stats.signed_no_internal_repeats;
    }

    const u64 depth =
        max_signed_m(
            excursion
        );

    if (depth > stats.max_signed_depth)
        stats.max_signed_depth = depth;

    const u64 length =
        static_cast<u64>(
            excursion.states.size()
        );

    if (length >
        stats.max_signed_excursion_length)
    {
        stats.max_signed_excursion_length =
            length;
    }

    const u64 n =
        static_cast<u64>(
            excursion.transitions.size()
        );

    for (u64 i = 0;
         i < n / 2;
         ++i)
    {
        const Transition& outward =
            excursion.transitions[i];

        const Transition& inward =
            excursion.transitions[
                n - 1 - i
            ];

        if (same_signed_transition_reverse(
                outward,
                inward))
        {
            ++stats.signed_state_reverse_matches;
        }
        else
        {
            ++stats.signed_state_reverse_failures;

            if (!stats.have_first_signed_reverse_failure)
            {
                stats.have_first_signed_reverse_failure =
                    true;

                stats.first_signed_reverse_failure_r =
                    r;

                stats.first_signed_reverse_failure =
                    excursion;
            }
        }

        if (same_divisor_pair_reverse(
                r,
                outward,
                inward))
        {
            ++stats.signed_divisor_reverse_matches;
        }
        else
        {
            ++stats.signed_divisor_reverse_failures;

            if (!stats.have_first_signed_divisor_failure)
            {
                stats.have_first_signed_divisor_failure =
                    true;

                stats.first_signed_divisor_failure_r =
                    r;

                stats.first_signed_divisor_failure =
                    excursion;
            }
        }

        if (outward.old_w.k ==
                inward.new_w.k &&
            outward.new_w.k ==
                inward.old_w.k)
        {
            ++stats.signed_k_reverse_matches;
        }
        else
        {
            ++stats.signed_k_reverse_failures;
        }

        if (signed_t_reverse(
                outward,
                inward))
        {
            ++stats.signed_t_reverse_matches;
        }
        else
        {
            ++stats.signed_t_reverse_failures;
        }

        if (same_delta_reverse(
                outward,
                inward))
        {
            ++stats.signed_delta_reverse_matches;
        }
        else
        {
            ++stats.signed_delta_reverse_failures;
        }

        if (same_correction_reverse(
                outward,
                inward))
        {
            ++stats.signed_correction_reverse_matches;
        }
        else
        {
            ++stats.signed_correction_reverse_failures;
        }

        if (same_K_reverse(
                outward,
                inward))
        {
            ++stats.signed_K_reverse_matches;
        }
        else
        {
            ++stats.signed_K_reverse_failures;
        }

        if (same_signed_transition_reverse(
                outward,
                inward))
        {
            ++stats.signed_transition_reverse_matches;
        }
        else
        {
            ++stats.signed_transition_reverse_failures;
        }
    }
}

static void process_prime(
    u64 r,
    Stats& stats)
{
    std::vector<Witness> best(
        M_MAX + 1
    );

    std::vector<Witness> trajectory_states;
    std::vector<Transition> trajectory_events;

    Witness previous_winner;
    bool have_previous = false;

    for (u64 K = 1;
         K <= K_LIMIT;
         ++K)
    {
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

        if (!winner.valid)
            continue;

        if (!have_previous)
        {
            trajectory_states.push_back(
                winner
            );

            previous_winner = winner;
            have_previous = true;
            continue;
        }

        if (winner.m != previous_winner.m)
        {
            Transition tr;

            tr.K = K;
            tr.old_w = previous_winner;
            tr.new_w = winner;

            trajectory_events.push_back(
                tr
            );

            trajectory_states.push_back(
                winner
            );
        }

        previous_winner = winner;
    }

    if (trajectory_states.empty())
        return;

    ++stats.trajectories;

    if (trajectory_states.size() >
        stats.max_signed_trajectory_length)
    {
        stats.max_signed_trajectory_length =
            static_cast<u64>(
                trajectory_states.size()
            );
    }

    /*
        Split into complete excursions from m=1 back to m=1.
    */

    std::size_t start = 0;

    while (start < trajectory_states.size())
    {
        if (trajectory_states[start].m != 1)
        {
            ++start;
            continue;
        }

        std::size_t end = start + 1;

        while (end < trajectory_states.size() &&
               trajectory_states[end].m != 1)
        {
            ++end;
        }

        if (end < trajectory_states.size() &&
            end > start + 1)
        {
            Excursion excursion;

            for (std::size_t i = start;
                 i <= end;
                 ++i)
            {
                excursion.states.push_back(
                    trajectory_states[i]
                );
            }

            /*
                Transition i connects state i to i+1.
            */
            for (std::size_t i = start;
                 i < end;
                 ++i)
            {
                excursion.transitions.push_back(
                    trajectory_events[i]
                );
            }

            ++stats.projected_excursions;
            ++stats.signed_excursions;

            analyze_excursion(
                r,
                excursion,
                stats
            );
        }

        if (end < trajectory_states.size())
            start = end;
        else
            break;
    }

    /*
        Count multiple excursions for this trajectory.
    */
    u64 excursion_count = 0;

    std::size_t p = 0;

    while (p < trajectory_states.size())
    {
        if (trajectory_states[p].m != 1)
        {
            ++p;
            continue;
        }

        std::size_t q = p + 1;

        while (q < trajectory_states.size() &&
               trajectory_states[q].m != 1)
        {
            ++q;
        }

        if (q < trajectory_states.size() &&
            q > p + 1)
        {
            ++excursion_count;
        }

        if (q < trajectory_states.size())
            p = q;
        else
            break;
    }

    if (excursion_count > 1)
        ++stats.signed_multiple_excursions;
}

int main()
{
    std::cout
        << "START EXPERIMENT 420\n";

    const std::vector<u64> primes =
        generate_primes(
            PRIME_LIMIT
        );

    Stats stats;

    stats.prime_count =
        static_cast<u64>(
            primes.size()
        );

    stats.total_points =
        stats.prime_count * K_LIMIT;

    for (u64 r : primes)
    {
        process_prime(
            r,
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
        << "PROJECTED_EXCURSIONS="
        << stats.projected_excursions
        << '\n'
        << "SIGNED_EXCURSIONS="
        << stats.signed_excursions
        << '\n'
        << '\n';

    std::cout
        << "PROJECTED_PALINDROMES="
        << stats.projected_palindromes
        << '\n'
        << "PROJECTED_NONPALINDROMES="
        << stats.projected_nonpalindromes
        << '\n'
        << "SIGNED_PALINDROMES="
        << stats.signed_palindromes
        << '\n'
        << "SIGNED_NONPALINDROMES="
        << stats.signed_nonpalindromes
        << '\n'
        << '\n';

    std::cout
        << "PROJECTED_ONLY_PALINDROME="
        << stats.projected_only_palindrome
        << '\n'
        << "SIGNED_ONLY_PALINDROME="
        << stats.signed_only_palindrome
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_STATE_REVERSE_MATCHES="
        << stats.signed_state_reverse_matches
        << '\n'
        << "SIGNED_STATE_REVERSE_FAILURES="
        << stats.signed_state_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_TRANSITION_REVERSE_MATCHES="
        << stats.signed_transition_reverse_matches
        << '\n'
        << "SIGNED_TRANSITION_REVERSE_FAILURES="
        << stats.signed_transition_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_DIVISOR_REVERSE_MATCHES="
        << stats.signed_divisor_reverse_matches
        << '\n'
        << "SIGNED_DIVISOR_REVERSE_FAILURES="
        << stats.signed_divisor_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_K_REVERSE_MATCHES="
        << stats.signed_k_reverse_matches
        << '\n'
        << "SIGNED_K_REVERSE_FAILURES="
        << stats.signed_k_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_T_REVERSE_MATCHES="
        << stats.signed_t_reverse_matches
        << '\n'
        << "SIGNED_T_REVERSE_FAILURES="
        << stats.signed_t_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_DELTA_REVERSE_MATCHES="
        << stats.signed_delta_reverse_matches
        << '\n'
        << "SIGNED_DELTA_REVERSE_FAILURES="
        << stats.signed_delta_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_CORRECTION_REVERSE_MATCHES="
        << stats.signed_correction_reverse_matches
        << '\n'
        << "SIGNED_CORRECTION_REVERSE_FAILURES="
        << stats.signed_correction_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_K_EVENT_REVERSE_MATCHES="
        << stats.signed_K_reverse_matches
        << '\n'
        << "SIGNED_K_EVENT_REVERSE_FAILURES="
        << stats.signed_K_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "SIGNED_INTERNAL_REPEATS="
        << stats.signed_internal_repeats
        << '\n'
        << "SIGNED_NO_INTERNAL_REPEATS="
        << stats.signed_no_internal_repeats
        << '\n'
        << "SIGNED_MULTIPLE_EXCURSIONS="
        << stats.signed_multiple_excursions
        << '\n'
        << '\n';

    std::cout
        << "MAX_SIGNED_TRAJECTORY_LENGTH="
        << stats.max_signed_trajectory_length
        << '\n'
        << "MAX_SIGNED_EXCURSION_LENGTH="
        << stats.max_signed_excursion_length
        << '\n'
        << "MAX_SIGNED_DEPTH="
        << stats.max_signed_depth
        << '\n';

    if (stats.have_first_signed_nonpalindrome)
    {
        std::cout
            << "FIRST_SIGNED_NONPALINDROME_R="
            << stats.first_signed_nonpalindrome_r
            << '\n';

        print_excursion(
            "FIRST_SIGNED_NONPALINDROME",
            stats.first_signed_nonpalindrome
        );

        std::cout << '\n';
    }

    if (stats.have_first_signed_reverse_failure)
    {
        std::cout
            << "FIRST_SIGNED_REVERSE_FAILURE_R="
            << stats.first_signed_reverse_failure_r
            << '\n';

        print_excursion(
            "FIRST_SIGNED_REVERSE_FAILURE",
            stats.first_signed_reverse_failure
        );

        std::cout << '\n';
    }

    if (stats.have_first_signed_divisor_failure)
    {
        std::cout
            << "FIRST_SIGNED_DIVISOR_FAILURE_R="
            << stats.first_signed_divisor_failure_r
            << '\n';

        print_excursion(
            "FIRST_SIGNED_DIVISOR_FAILURE",
            stats.first_signed_divisor_failure
        );

        std::cout << '\n';
    }

    std::cout
        << "PROJECTED_PALINDROME_STATUS="
        << (
            stats.projected_nonpalindromes == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "SIGNED_PALINDROME_STATUS="
        << (
            stats.signed_nonpalindromes == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "SIGNED_STATE_REVERSAL_STATUS="
        << (
            stats.signed_state_reverse_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "SIGNED_DIVISOR_REVERSAL_STATUS="
        << (
            stats.signed_divisor_reverse_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 420\n";

    return 0;
}
