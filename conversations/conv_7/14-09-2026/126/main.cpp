#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
#include <numeric>

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

struct TransitionEvent
{
    u64 K = 0;

    Witness old_w;
    Witness new_w;

    u64 delta = 0;
    i128 correction = 0;
    u64 threshold = 0;

    u64 exact_drop = 0;
};

struct Excursion
{
    std::vector<u64> states;
    std::vector<TransitionEvent> events;
};

struct Stats
{
    u64 prime_count = 0;
    u64 total_points = 0;

    u64 excursions = 0;
    u64 mirrored_transition_pairs = 0;

    u64 state_reverse_matches = 0;
    u64 state_reverse_failures = 0;

    u64 delta_equal = 0;
    u64 delta_not_equal = 0;

    u64 correction_equal = 0;
    u64 correction_not_equal = 0;

    u64 exact_drop_equal = 0;
    u64 exact_drop_not_equal = 0;

    u64 k_old_equal = 0;
    u64 k_old_not_equal = 0;

    u64 k_new_equal = 0;
    u64 k_new_not_equal = 0;

    u64 sign_exact_reverse_matches = 0;
    u64 sign_exact_reverse_failures = 0;

    u64 threshold_equal = 0;
    u64 threshold_not_equal = 0;

    u64 K_equal = 0;
    u64 K_not_equal = 0;

    u64 divisor_value_reverse_matches = 0;
    u64 divisor_value_reverse_failures = 0;

    u64 product_k_equal = 0;
    u64 product_k_not_equal = 0;

    u64 t_sum_equal = 0;
    u64 t_sum_not_equal = 0;

    u64 all_invariants_match = 0;
    u64 any_invariant_mismatch = 0;

    u64 first_mismatch_r = 0;
    Excursion first_mismatch;

    bool have_first_mismatch = false;

    u64 strongest_relation_matches = 0;
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

    if (make_witness(
            r,
            m,
            k,
            +1,
            candidate))
    {
        if (!best.valid ||
            candidate.k > best.k ||
            (candidate.k == best.k &&
             candidate.t < best.t))
        {
            best = candidate;
        }
    }

    if (make_witness(
            r,
            m,
            k,
            -1,
            candidate))
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

static u64 compute_delta(
    const Witness& old_w,
    const Witness& new_w)
{
    const i128 value =
        static_cast<i128>(old_w.m) * new_w.k -
        static_cast<i128>(new_w.m) * old_w.k;

    if (value <= 0)
        return 0;

    return static_cast<u64>(value);
}

static i128 compute_correction(
    const Witness& old_w,
    const Witness& new_w)
{
    return
        static_cast<i128>(old_w.sign) * new_w.k -
        static_cast<i128>(new_w.sign) * old_w.k;
}

static u64 compute_threshold(
    u64 delta,
    i128 correction)
{
    if (delta == 0)
        return 0;

    if (correction < 0)
        return 1;

    const i128 d =
        static_cast<i128>(delta);

    const i128 threshold =
        correction / d + 1;

    if (threshold <= 1)
        return 1;

    if (threshold >
        static_cast<i128>(
            std::numeric_limits<u64>::max()))
    {
        return std::numeric_limits<u64>::max();
    }

    return static_cast<u64>(threshold);
}

static u64 exact_drop(
    const Witness& old_w,
    const Witness& new_w)
{
    if (old_w.t <= new_w.t)
        return 0;

    return old_w.t - new_w.t;
}

static TransitionEvent build_event(
    u64 K,
    const Witness& old_w,
    const Witness& new_w)
{
    TransitionEvent event;

    event.K = K;
    event.old_w = old_w;
    event.new_w = new_w;

    event.delta =
        compute_delta(
            old_w,
            new_w
        );

    event.correction =
        compute_correction(
            old_w,
            new_w
        );

    event.threshold =
        compute_threshold(
            event.delta,
            event.correction
        );

    event.exact_drop =
        exact_drop(
            old_w,
            new_w
        );

    return event;
}

static u64 plus_value(
    const Witness& w,
    u64 r)
{
    /*
        The represented divisor is

            m*r + sign.
    */

    const i128 value =
        static_cast<i128>(w.m) * r +
        static_cast<i128>(w.sign);

    return static_cast<u64>(value);
}

static bool divisor_relation_matches(
    u64 r,
    const TransitionEvent& outward,
    const TransitionEvent& inward)
{
    /*
        The return transition reverses the m-state edge:

            a -> b

        versus

            b -> a.

        Check whether the divisor attached to the new
        outward witness equals the divisor attached to the
        old return witness, and vice versa.
    */

    const u64 outward_new_value =
        plus_value(
            outward.new_w,
            r
        );

    const u64 inward_old_value =
        plus_value(
            inward.old_w,
            r
        );

    const u64 outward_old_value =
        plus_value(
            outward.old_w,
            r
        );

    const u64 inward_new_value =
        plus_value(
            inward.new_w,
            r
        );

    const bool first =
        outward_new_value ==
        inward_old_value;

    const bool second =
        outward_old_value ==
        inward_new_value;

    return first && second;
}

static bool basic_reverse_state_match(
    const TransitionEvent& outward,
    const TransitionEvent& inward)
{
    return
        outward.old_w.m ==
            inward.new_w.m &&
        outward.new_w.m ==
            inward.old_w.m;
}

static bool sign_reverse_match(
    const TransitionEvent& outward,
    const TransitionEvent& inward)
{
    return
        outward.old_w.sign ==
            inward.new_w.sign &&
        outward.new_w.sign ==
            inward.old_w.sign;
}

static void compare_mirrored_events(
    u64 r,
    const TransitionEvent& outward,
    const TransitionEvent& inward,
    Stats& stats,
    Excursion& excursion)
{
    ++stats.mirrored_transition_pairs;

    const bool state_match =
        basic_reverse_state_match(
            outward,
            inward
        );

    if (state_match)
        ++stats.state_reverse_matches;
    else
        ++stats.state_reverse_failures;

    if (outward.delta == inward.delta)
        ++stats.delta_equal;
    else
        ++stats.delta_not_equal;

    if (outward.correction == inward.correction)
        ++stats.correction_equal;
    else
        ++stats.correction_not_equal;

    if (outward.exact_drop == inward.exact_drop)
        ++stats.exact_drop_equal;
    else
        ++stats.exact_drop_not_equal;

    if (outward.old_w.k == inward.new_w.k)
        ++stats.k_old_equal;
    else
        ++stats.k_old_not_equal;

    if (outward.new_w.k == inward.old_w.k)
        ++stats.k_new_equal;
    else
        ++stats.k_new_not_equal;

    if (sign_reverse_match(outward, inward))
        ++stats.sign_exact_reverse_matches;
    else
        ++stats.sign_exact_reverse_failures;

    if (outward.threshold == inward.threshold)
        ++stats.threshold_equal;
    else
        ++stats.threshold_not_equal;

    if (outward.K == inward.K)
        ++stats.K_equal;
    else
        ++stats.K_not_equal;

    if (divisor_relation_matches(
            r,
            outward,
            inward))
    {
        ++stats.divisor_value_reverse_matches;
    }
    else
    {
        ++stats.divisor_value_reverse_failures;
    }

    const u128 outward_product =
        static_cast<u128>(
            outward.old_w.k) *
        outward.new_w.k;

    const u128 inward_product =
        static_cast<u128>(
            inward.old_w.k) *
        inward.new_w.k;

    if (outward_product == inward_product)
        ++stats.product_k_equal;
    else
        ++stats.product_k_not_equal;

    const u128 outward_t_sum =
        static_cast<u128>(
            outward.old_w.t) +
        outward.new_w.t;

    const u128 inward_t_sum =
        static_cast<u128>(
            inward.old_w.t) +
        inward.new_w.t;

    if (outward_t_sum == inward_t_sum)
        ++stats.t_sum_equal;
    else
        ++stats.t_sum_not_equal;

    const bool all_core =
        state_match &&
        outward.delta == inward.delta &&
        outward.correction == inward.correction &&
        outward.exact_drop == inward.exact_drop &&
        outward.old_w.k == inward.new_w.k &&
        outward.new_w.k == inward.old_w.k &&
        sign_reverse_match(
            outward,
            inward
        ) &&
        divisor_relation_matches(
            r,
            outward,
            inward
        );

    if (all_core)
    {
        ++stats.all_invariants_match;
        ++stats.strongest_relation_matches;
    }
    else
    {
        ++stats.any_invariant_mismatch;

        if (!stats.have_first_mismatch)
        {
            stats.have_first_mismatch = true;
            stats.first_mismatch_r = r;
            stats.first_mismatch = excursion;
        }
    }
}

static void analyze_excursion(
    u64 r,
    const std::vector<u64>& states,
    const std::vector<TransitionEvent>& events,
    Stats& stats)
{
    if (states.size() < 3)
        return;

    if (events.size() + 1 != states.size())
        return;

    /*
        For a palindrome

            1,a,b,c,b,a,1

        the transition sequence is

            1->a
            a->b
            b->c
            c->b
            b->a
            a->1

        so event i should mirror event
        events.size()-1-i.
    */

    const std::size_t n = events.size();

    for (std::size_t i = 0;
         i < n / 2;
         ++i)
    {
        const TransitionEvent& outward =
            events[i];

        const TransitionEvent& inward =
            events[n - 1 - i];

        compare_mirrored_events(
            r,
            outward,
            inward,
            stats,
            const_cast<Excursion&>(
                *new Excursion{
                    states,
                    events
                })
        );
    }
}

static void process_prime(
    u64 r,
    Stats& stats)
{
    std::vector<Witness> best(
        M_MAX + 1
    );

    std::vector<u64> states;
    std::vector<TransitionEvent> events;

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
            states.push_back(winner.m);
            previous_winner = winner;
            have_previous = true;
            continue;
        }

        if (winner.m != previous_winner.m)
        {
            const TransitionEvent event =
                build_event(
                    K,
                    previous_winner,
                    winner
                );

            states.push_back(winner.m);
            events.push_back(event);
        }

        previous_winner = winner;
    }

    /*
        Find maximal excursions between state 1's.
    */

    std::size_t i = 0;

    while (i < states.size())
    {
        if (states[i] != 1)
        {
            ++i;
            continue;
        }

        std::size_t j = i + 1;

        while (j < states.size() &&
               states[j] != 1)
        {
            ++j;
        }

        if (j < states.size() &&
            j > i + 1)
        {
            /*
                State range:

                    states[i] ... states[j]

                Event range:

                    events[i] ... events[j-1]
            */

            std::vector<u64> excursion_states(
                states.begin() +
                    static_cast<std::ptrdiff_t>(i),
                states.begin() +
                    static_cast<std::ptrdiff_t>(j + 1)
            );

            std::vector<TransitionEvent> excursion_events(
                events.begin() +
                    static_cast<std::ptrdiff_t>(i),
                events.begin() +
                    static_cast<std::ptrdiff_t>(j)
            );

            ++stats.excursions;

            analyze_excursion(
                r,
                excursion_states,
                excursion_events,
                stats
            );
        }

        if (j < states.size())
            i = j;
        else
            break;
    }
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

static void print_event(
    const char* name,
    u64 r,
    const TransitionEvent& e)
{
    std::cout
        << name << ".K="
        << e.K << '\n'
        << name << ".delta="
        << e.delta << '\n'
        << name << ".correction="
        << static_cast<long long>(
            e.correction)
        << '\n'
        << name << ".threshold="
        << e.threshold << '\n'
        << name << ".exact_drop="
        << e.exact_drop << '\n';

    print_witness(
        name,
        e.old_w
    );

    print_witness(
        "NEW",
        e.new_w
    );

    std::cout
        << name << ".divisor_old="
        << plus_value(
            e.old_w,
            r
        )
        << '\n'
        << name << ".divisor_new="
        << plus_value(
            e.new_w,
            r
        )
        << '\n';
}

int main()
{
    std::cout
        << "START EXPERIMENT 418\n";

    const std::vector<u64> primes =
        generate_primes(
            PRIME_LIMIT
        );

    Stats stats;

    stats.prime_count =
        static_cast<u64>(primes.size());

    for (u64 r : primes)
    {
        stats.total_points += K_LIMIT;

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
        << "EXCURSIONS="
        << stats.excursions
        << '\n'
        << "MIRRORED_TRANSITION_PAIRS="
        << stats.mirrored_transition_pairs
        << '\n'
        << '\n';

    std::cout
        << "STATE_REVERSE_MATCHES="
        << stats.state_reverse_matches
        << '\n'
        << "STATE_REVERSE_FAILURES="
        << stats.state_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "DELTA_EQUAL="
        << stats.delta_equal
        << '\n'
        << "DELTA_NOT_EQUAL="
        << stats.delta_not_equal
        << '\n'
        << '\n';

    std::cout
        << "CORRECTION_EQUAL="
        << stats.correction_equal
        << '\n'
        << "CORRECTION_NOT_EQUAL="
        << stats.correction_not_equal
        << '\n'
        << '\n';

    std::cout
        << "EXACT_DROP_EQUAL="
        << stats.exact_drop_equal
        << '\n'
        << "EXACT_DROP_NOT_EQUAL="
        << stats.exact_drop_not_equal
        << '\n'
        << '\n';

    std::cout
        << "K_OLD_EQUAL_MIRROR_NEW="
        << stats.k_old_equal
        << '\n'
        << "K_OLD_NOT_EQUAL_MIRROR_NEW="
        << stats.k_old_not_equal
        << '\n'
        << "K_NEW_EQUAL_MIRROR_OLD="
        << stats.k_new_equal
        << '\n'
        << "K_NEW_NOT_EQUAL_MIRROR_OLD="
        << stats.k_new_not_equal
        << '\n'
        << '\n';

    std::cout
        << "SIGN_REVERSE_MATCHES="
        << stats.sign_exact_reverse_matches
        << '\n'
        << "SIGN_REVERSE_FAILURES="
        << stats.sign_exact_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "THRESHOLD_EQUAL="
        << stats.threshold_equal
        << '\n'
        << "THRESHOLD_NOT_EQUAL="
        << stats.threshold_not_equal
        << '\n'
        << '\n';

    std::cout
        << "K_EVENT_EQUAL="
        << stats.K_equal
        << '\n'
        << "K_EVENT_NOT_EQUAL="
        << stats.K_not_equal
        << '\n'
        << '\n';

    std::cout
        << "DIVISOR_VALUE_REVERSE_MATCHES="
        << stats.divisor_value_reverse_matches
        << '\n'
        << "DIVISOR_VALUE_REVERSE_FAILURES="
        << stats.divisor_value_reverse_failures
        << '\n'
        << '\n';

    std::cout
        << "PRODUCT_K_EQUAL="
        << stats.product_k_equal
        << '\n'
        << "PRODUCT_K_NOT_EQUAL="
        << stats.product_k_not_equal
        << '\n'
        << '\n';

    std::cout
        << "T_SUM_EQUAL="
        << stats.t_sum_equal
        << '\n'
        << "T_SUM_NOT_EQUAL="
        << stats.t_sum_not_equal
        << '\n'
        << '\n';

    std::cout
        << "ALL_CORE_INVARIANTS_MATCH="
        << stats.all_invariants_match
        << '\n'
        << "ANY_CORE_INVARIANT_MISMATCH="
        << stats.any_invariant_mismatch
        << '\n'
        << '\n';

    if (stats.have_first_mismatch)
    {
        std::cout
            << "FIRST_CORE_MISMATCH_R="
            << stats.first_mismatch_r
            << '\n';

        /*
            Reconstruct and print the first mismatch excursion.
        */

        const std::vector<u64>& states =
            stats.first_mismatch.states;

        const std::vector<TransitionEvent>& events =
            stats.first_mismatch.events;

        std::cout << "FIRST_MISMATCH_TRAJECTORY=";

        for (std::size_t i = 0;
             i < states.size();
             ++i)
        {
            if (i != 0)
                std::cout << ",";

            std::cout << states[i];
        }

        std::cout << '\n';

        if (!events.empty())
        {
            print_event(
                "OUTWARD",
                stats.first_mismatch_r,
                events.front()
            );

            print_event(
                "MIRROR",
                stats.first_mismatch_r,
                events.back()
            );
        }

        std::cout << '\n';
    }

    std::cout
        << "STATE_REVERSAL_STATUS="
        << (
            stats.state_reverse_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "DIVISOR_PAIR_REVERSAL_STATUS="
        << (
            stats.divisor_value_reverse_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "CORE_INVARIANT_STATUS="
        << (
            stats.any_invariant_mismatch == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 418\n";

    return 0;
}
