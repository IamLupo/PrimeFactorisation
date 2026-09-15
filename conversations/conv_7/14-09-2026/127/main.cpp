#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
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

struct TransitionEvent
{
    u64 K = 0;

    Witness old_w;
    Witness new_w;

    u64 old_value = 0;
    u64 new_value = 0;

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

struct PairTrack
{
    u64 m = 0;
    int sign = 0;
    u64 value = 0;

    /*
        At each K, best divisor <= K for this fixed
        arithmetic value.
    */
    std::vector<u64> best_k;
    std::vector<u64> best_t;
    std::vector<unsigned char> has_record;
};

struct PairAnalysis
{
    bool same_arithmetic_pair = false;

    bool outward_is_record_event = false;
    bool inward_is_record_event = false;

    bool outward_was_actual_two_track_winner = false;
    bool inward_was_actual_two_track_winner = false;

    bool envelope_reconstructed = false;

    u64 outward_record_count = 0;
    u64 inward_record_count = 0;

    u64 first_difference_K = 0;
};

struct Stats
{
    u64 prime_count = 0;
    u64 total_points = 0;

    u64 excursions = 0;
    u64 mirrored_pairs = 0;

    u64 same_arithmetic_pair = 0;
    u64 different_arithmetic_pair = 0;

    u64 outward_record_matches = 0;
    u64 outward_record_failures = 0;

    u64 inward_record_matches = 0;
    u64 inward_record_failures = 0;

    u64 two_track_reconstruction_matches = 0;
    u64 two_track_reconstruction_failures = 0;

    u64 mirrored_sign_matches = 0;
    u64 mirrored_sign_failures = 0;

    u64 mirrored_m_value_matches = 0;
    u64 mirrored_m_value_failures = 0;

    u64 first_difference_K_count = 0;

    u64 exceptional_pair_count = 0;

    u64 max_record_count = 0;
    u64 max_record_r = 0;

    bool have_first_exception = false;
    u64 first_exception_r = 0;
    Excursion first_exception;

    bool have_first_reconstruction_failure = false;
    u64 first_reconstruction_failure_r = 0;
    Excursion first_reconstruction_failure;

    bool have_max_records = false;
    PairTrack max_record_track;
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

static u64 witness_value(
    u64 r,
    const Witness& w)
{
    const i128 value =
        static_cast<i128>(w.m) * r +
        static_cast<i128>(w.sign);

    return static_cast<u64>(value);
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

    const i128 threshold =
        correction /
        static_cast<i128>(delta) + 1;

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

static TransitionEvent build_event(
    u64 r,
    u64 K,
    const Witness& old_w,
    const Witness& new_w)
{
    TransitionEvent event;

    event.K = K;
    event.old_w = old_w;
    event.new_w = new_w;

    event.old_value =
        witness_value(r, old_w);

    event.new_value =
        witness_value(r, new_w);

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

    if (old_w.t > new_w.t)
        event.exact_drop =
            old_w.t - new_w.t;

    return event;
}

static PairTrack build_pair_track(
    u64 r,
    u64 m,
    int sign)
{
    PairTrack track;

    track.m = m;
    track.sign = sign;

    const i128 value =
        static_cast<i128>(m) * r +
        static_cast<i128>(sign);

    track.value =
        static_cast<u64>(value);

    track.best_k.assign(
        K_LIMIT + 1,
        0
    );

    track.best_t.assign(
        K_LIMIT + 1,
        0
    );

    track.has_record.assign(
        K_LIMIT + 1,
        0
    );

    u64 current_k = 0;
    u64 current_t = 0;

    for (u64 K = 1;
         K <= K_LIMIT;
         ++K)
    {
        if (track.value % K == 0)
        {
            const u64 t =
                track.value / K;

            if (K > current_k)
            {
                current_k = K;
                current_t = t;
                track.has_record[K] = 1;
            }
        }

        track.best_k[K] = current_k;
        track.best_t[K] = current_t;
    }

    return track;
}

static bool same_unordered_values(
    const TransitionEvent& a,
    const TransitionEvent& b)
{
    return
        (a.old_value == b.new_value &&
         a.new_value == b.old_value);
}

static bool same_unordered_states(
    const TransitionEvent& a,
    const TransitionEvent& b)
{
    return
        (a.old_w.m == b.new_w.m &&
         a.new_w.m == b.old_w.m);
}

static bool same_unordered_signs(
    const TransitionEvent& a,
    const TransitionEvent& b)
{
    return
        (a.old_w.sign == b.new_w.sign &&
         a.new_w.sign == b.old_w.sign);
}

static bool event_uses_record(
    const PairTrack& old_track,
    const PairTrack& new_track,
    const TransitionEvent& event)
{
    if (event.K > K_LIMIT)
        return false;

    const u64 old_best =
        old_track.best_k[event.K];

    const u64 new_best =
        new_track.best_k[event.K];

    return
        old_best == event.old_w.k &&
        new_best == event.new_w.k;
}

static bool two_track_winner_at_K(
    const PairTrack& old_track,
    const PairTrack& new_track,
    u64 K,
    Witness& winner)
{
    winner = Witness{};

    if (K == 0 ||
        K > K_LIMIT)
        return false;

    const u64 old_k =
        old_track.best_k[K];

    const u64 new_k =
        new_track.best_k[K];

    if (old_k == 0 &&
        new_k == 0)
    {
        return false;
    }

    if (old_k != 0)
    {
        Witness old_w;

        const i128 value =
            static_cast<i128>(
                old_track.m) * 1;

        (void)value;

        old_w.valid = true;
        old_w.m = old_track.m;
        old_w.k = old_k;
        old_w.sign = old_track.sign;
        old_w.t =
            old_track.best_t[K];

        consider_winner(
            old_w,
            winner
        );
    }

    if (new_k != 0)
    {
        Witness new_w;

        new_w.valid = true;
        new_w.m = new_track.m;
        new_w.k = new_k;
        new_w.sign = new_track.sign;
        new_w.t =
            new_track.best_t[K];

        consider_winner(
            new_w,
            winner
        );
    }

    return winner.valid;
}

static PairAnalysis analyze_pair(
    u64 r,
    const TransitionEvent& outward,
    const TransitionEvent& inward)
{
    PairAnalysis result;

    result.same_arithmetic_pair =
        same_unordered_values(
            outward,
            inward
        );

    result.outward_is_record_event = false;
    result.inward_is_record_event = false;

    /*
        Tracks will be built by the caller only for the
        same-arithmetic-pair case.
    */

    (void)r;

    return result;
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
    const TransitionEvent& event)
{
    std::cout
        << name << ".K="
        << event.K << '\n'
        << name << ".old_value="
        << event.old_value << '\n'
        << name << ".new_value="
        << event.new_value << '\n'
        << name << ".delta="
        << event.delta << '\n'
        << name << ".correction="
        << static_cast<long long>(
            event.correction)
        << '\n'
        << name << ".threshold="
        << event.threshold << '\n'
        << name << ".exact_drop="
        << event.exact_drop << '\n';

    print_witness(
        name,
        event.old_w
    );

    print_witness(
        "NEW",
        event.new_w
    );
}

static void compare_pair(
    u64 r,
    const TransitionEvent& outward,
    const TransitionEvent& inward,
    Stats& stats,
    const Excursion& full_excursion)
{
    ++stats.mirrored_pairs;

    const bool same_states =
        same_unordered_states(
            outward,
            inward
        );

    const bool same_values =
        same_unordered_values(
            outward,
            inward
        );

    const bool same_signs =
        same_unordered_signs(
            outward,
            inward
        );

    if (same_values)
        ++stats.same_arithmetic_pair;
    else
        ++stats.different_arithmetic_pair;

    if (same_signs)
        ++stats.mirrored_sign_matches;
    else
        ++stats.mirrored_sign_failures;

    PairAnalysis analysis =
        analyze_pair(
            r,
            outward,
            inward
        );

    if (same_values)
    {
        /*
            The two arithmetic values are fixed:

                A = outward.old_value
                B = outward.new_value

            On the return they appear in reverse order.

            Build the exact divisor-record histories of A and B.
        */

        PairTrack A =
            build_pair_track(
                r,
                outward.old_w.m,
                outward.old_w.sign
            );

        PairTrack B =
            build_pair_track(
                r,
                outward.new_w.m,
                outward.new_w.sign
            );

        const bool outward_record =
            event_uses_record(
                A,
                B,
                outward
            );

        /*
            The mirrored event has A/B reversed.
        */
        const bool inward_record =
            event_uses_record(
                B,
                A,
                inward
            );

        if (outward_record)
            ++stats.outward_record_matches;
        else
            ++stats.outward_record_failures;

        if (inward_record)
            ++stats.inward_record_matches;
        else
            ++stats.inward_record_failures;

        /*
            Check whether the two events are actual changes in
            the two-track normalized winner at their K values.
        */

        Witness outward_two_track;
        Witness inward_two_track;

        const bool outward_available =
            two_track_winner_at_K(
                A,
                B,
                outward.K,
                outward_two_track
            );

        const bool inward_available =
            two_track_winner_at_K(
                B,
                A,
                inward.K,
                inward_two_track
            );

        const bool outward_winner =
            outward_available &&
            outward_two_track.m ==
                outward.new_w.m &&
            outward_two_track.k ==
                outward.new_w.k;

        const bool inward_winner =
            inward_available &&
            inward_two_track.m ==
                inward.new_w.m &&
            inward_two_track.k ==
                inward.new_w.k;

        if (outward_winner &&
            inward_winner)
        {
            ++stats.two_track_reconstruction_matches;
        }
        else
        {
            ++stats.two_track_reconstruction_failures;

            if (!stats.have_first_reconstruction_failure)
            {
                stats.have_first_reconstruction_failure =
                    true;

                stats.first_reconstruction_failure_r = r;
                stats.first_reconstruction_failure =
                    full_excursion;
            }
        }
    }
    else
    {
        ++stats.two_track_reconstruction_failures;
    }

    if (same_values &&
        same_states &&
        same_signs)
    {
        /*
            This is the strongest simple arithmetic-pair
            symmetry we can establish directly.
        */
    }

    /*
        Print/store the first exceptional arithmetic pair.
    */
    if (!same_values)
    {
        ++stats.exceptional_pair_count;

        if (!stats.have_first_exception)
        {
            stats.have_first_exception = true;
            stats.first_exception_r = r;
            stats.first_exception =
                full_excursion;
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
            TransitionEvent event =
                build_event(
                    r,
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
        Identify each complete excursion.
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

            Excursion full_excursion;

            full_excursion.states =
                excursion_states;

            full_excursion.events =
                excursion_events;

            ++stats.excursions;

            const std::size_t n =
                excursion_events.size();

            for (std::size_t p = 0;
                 p < n / 2;
                 ++p)
            {
                const TransitionEvent& outward =
                    excursion_events[p];

                const TransitionEvent& inward =
                    excursion_events[n - 1 - p];

                compare_pair(
                    r,
                    outward,
                    inward,
                    stats,
                    full_excursion
                );
            }
        }

        if (j < states.size())
            i = j;
        else
            break;
    }
}

int main()
{
    std::cout
        << "START EXPERIMENT 419\n";

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
        << "EXCURSIONS="
        << stats.excursions
        << '\n'
        << "MIRRORED_PAIRS="
        << stats.mirrored_pairs
        << '\n'
        << '\n';

    std::cout
        << "SAME_ARITHMETIC_PAIR="
        << stats.same_arithmetic_pair
        << '\n'
        << "DIFFERENT_ARITHMETIC_PAIR="
        << stats.different_arithmetic_pair
        << '\n'
        << '\n';

    std::cout
        << "OUTWARD_RECORD_MATCHES="
        << stats.outward_record_matches
        << '\n'
        << "OUTWARD_RECORD_FAILURES="
        << stats.outward_record_failures
        << '\n'
        << "INWARD_RECORD_MATCHES="
        << stats.inward_record_matches
        << '\n'
        << "INWARD_RECORD_FAILURES="
        << stats.inward_record_failures
        << '\n'
        << '\n';

    std::cout
        << "TWO_TRACK_RECONSTRUCTION_MATCHES="
        << stats.two_track_reconstruction_matches
        << '\n'
        << "TWO_TRACK_RECONSTRUCTION_FAILURES="
        << stats.two_track_reconstruction_failures
        << '\n'
        << '\n';

    std::cout
        << "MIRRORED_SIGN_MATCHES="
        << stats.mirrored_sign_matches
        << '\n'
        << "MIRRORED_SIGN_FAILURES="
        << stats.mirrored_sign_failures
        << '\n'
        << '\n';

    std::cout
        << "EXCEPTIONAL_PAIR_COUNT="
        << stats.exceptional_pair_count
        << '\n'
        << '\n';

    if (stats.have_first_exception)
    {
        std::cout
            << "FIRST_ARITHMETIC_PAIR_EXCEPTION_R="
            << stats.first_exception_r
            << '\n';

        std::cout
            << "TRAJECTORY=";

        for (std::size_t p = 0;
             p < stats.first_exception.states.size();
             ++p)
        {
            if (p != 0)
                std::cout << ",";

            std::cout
                << stats.first_exception.states[p];
        }

        std::cout << '\n';

        for (const TransitionEvent& event :
             stats.first_exception.events)
        {
            print_event(
                "EVENT",
                event
            );
        }

        std::cout << '\n';
    }

    std::cout
        << "SAME_ARITHMETIC_PAIR_STATUS="
        << (
            stats.different_arithmetic_pair == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "OUTWARD_RECORD_STATUS="
        << (
            stats.outward_record_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "INWARD_RECORD_STATUS="
        << (
            stats.inward_record_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "TWO_TRACK_STATUS="
        << (
            stats.two_track_reconstruction_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 419\n";

    return 0;
}
