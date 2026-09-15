#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

static constexpr u64 PRIME_LIMIT = 5000;
static constexpr u64 K_LIMIT     = 500;
static constexpr u64 M_MAX       = 64;

struct Witness
{
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
    bool valid = false;
};

struct Stats
{
    u64 prime_count = 0;

    u64 k_points = 0;

    u64 normalized_transitions = 0;
    u64 normalized_delta_one = 0;
    u64 normalized_delta_gt_one = 0;

    u64 exact_transitions = 0;
    u64 exact_delta_one = 0;
    u64 exact_delta_gt_one = 0;

    u64 normalized_new_record_transitions = 0;
    u64 normalized_transitions_without_new_record = 0;

    u64 exact_new_record_transitions = 0;
    u64 exact_transitions_without_new_record = 0;

    u64 normalized_forward_transitions = 0;
    u64 normalized_backward_transitions = 0;

    u64 exact_forward_transitions = 0;
    u64 exact_backward_transitions = 0;

    u64 normalized_delta_min = std::numeric_limits<u64>::max();
    u64 normalized_delta_max = 0;

    u64 exact_delta_min = std::numeric_limits<u64>::max();
    u64 exact_delta_max = 0;

    u64 normalized_max_k_sum = 0;
    u64 exact_max_k_sum = 0;

    bool have_normalized_delta_counterexample = false;
    u64 counter_r = 0;
    u64 counter_k = 0;
    Witness counter_old;
    Witness counter_new;

    bool have_exact_delta_counterexample = false;
    u64 exact_counter_r = 0;
    u64 exact_counter_k = 0;
    Witness exact_counter_old;
    Witness exact_counter_new;
};

static bool is_prime(u64 n)
{
    if (n < 2)
        return false;

    if (n == 2)
        return true;

    if (n % 2 == 0)
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

static Witness best_witness_for_m(
    u64 r,
    u64 m,
    u64 K)
{
    Witness best;

    for (u64 k = 1; k <= K; ++k)
    {
        Witness plus;
        Witness minus;

        if (make_witness(r, m, k, +1, plus))
        {
            if (!best.valid ||
                plus.k > best.k ||
                (plus.k == best.k &&
                 plus.t < best.t))
            {
                best = plus;
            }
        }

        if (make_witness(r, m, k, -1, minus))
        {
            if (!best.valid ||
                minus.k > best.k ||
                (minus.k == best.k &&
                 minus.t < best.t))
            {
                best = minus;
            }
        }
    }

    return best;
}

static bool normalized_better(
    const Witness& a,
    const Witness& b)
{
    /*
        a wins if

            ka/ma > kb/mb.

        Exact cross multiplication avoids floating point.
    */
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

static bool exact_better(
    const Witness& a,
    const Witness& b)
{
    return a.t < b.t;
}

static bool exact_equal(
    const Witness& a,
    const Witness& b)
{
    return a.t == b.t;
}

static u64 transition_delta(
    const Witness& old_w,
    const Witness& new_w)
{
    /*
        New normalized ratio is assumed to be larger:

            knew/mnew > kold/mold

        Therefore

            Delta =
                mold*knew - mnew*kold
    */
    const i128 delta =
        static_cast<i128>(old_w.m) * new_w.k -
        static_cast<i128>(new_w.m) * old_w.k;

    if (delta <= 0)
        return 0;

    return static_cast<u64>(delta);
}

static void update_normalized_transition(
    u64 r,
    u64 K,
    const Witness& old_w,
    const Witness& new_w,
    bool new_record,
    Stats& stats)
{
    ++stats.normalized_transitions;

    if (new_record)
        ++stats.normalized_new_record_transitions;
    else
        ++stats.normalized_transitions_without_new_record;

    const u64 delta =
        transition_delta(old_w, new_w);

    if (delta == 1)
    {
        ++stats.normalized_delta_one;
    }
    else
    {
        ++stats.normalized_delta_gt_one;

        if (!stats.have_normalized_delta_counterexample)
        {
            stats.have_normalized_delta_counterexample = true;
            stats.counter_r = r;
            stats.counter_k = K;
            stats.counter_old = old_w;
            stats.counter_new = new_w;
        }
    }

    if (delta < stats.normalized_delta_min)
        stats.normalized_delta_min = delta;

    if (delta > stats.normalized_delta_max)
        stats.normalized_delta_max = delta;

    const u64 k_sum =
        old_w.k + new_w.k;

    if (k_sum > stats.normalized_max_k_sum)
        stats.normalized_max_k_sum = k_sum;

    const bool forward =
        normalized_better(new_w, old_w);

    if (forward)
        ++stats.normalized_forward_transitions;
    else
        ++stats.normalized_backward_transitions;
}

static void update_exact_transition(
    u64 r,
    u64 K,
    const Witness& old_w,
    const Witness& new_w,
    bool new_record,
    Stats& stats)
{
    ++stats.exact_transitions;

    if (new_record)
        ++stats.exact_new_record_transitions;
    else
        ++stats.exact_transitions_without_new_record;

    /*
        For the exact winner, the new winner can move in
        either normalized direction. Only compute Delta when
        the new ratio is larger, which is the direction used
        by the normalized transition theorem.
    */
    if (normalized_better(new_w, old_w))
    {
        const u64 delta =
            transition_delta(old_w, new_w);

        if (delta == 1)
            ++stats.exact_delta_one;
        else
        {
            ++stats.exact_delta_gt_one;

            if (!stats.have_exact_delta_counterexample)
            {
                stats.have_exact_delta_counterexample = true;
                stats.exact_counter_r = r;
                stats.exact_counter_k = K;
                stats.exact_counter_old = old_w;
                stats.exact_counter_new = new_w;
            }
        }

        if (delta < stats.exact_delta_min)
            stats.exact_delta_min = delta;

        if (delta > stats.exact_delta_max)
            stats.exact_delta_max = delta;

        const u64 k_sum =
            old_w.k + new_w.k;

        if (k_sum > stats.exact_max_k_sum)
            stats.exact_max_k_sum = k_sum;
    }

    if (exact_better(new_w, old_w))
        ++stats.exact_forward_transitions;
    else
        ++stats.exact_backward_transitions;
}

static Witness choose_normalized_winner(
    const std::vector<Witness>& current)
{
    Witness winner;

    for (u64 m = 1; m < current.size(); ++m)
    {
        const Witness& candidate = current[m];

        if (!candidate.valid)
            continue;

        if (!winner.valid)
        {
            winner = candidate;
            continue;
        }

        if (normalized_better(candidate, winner))
        {
            winner = candidate;
            continue;
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

    return winner;
}

static Witness choose_exact_winner(
    const std::vector<Witness>& current)
{
    Witness winner;

    for (u64 m = 1; m < current.size(); ++m)
    {
        const Witness& candidate = current[m];

        if (!candidate.valid)
            continue;

        if (!winner.valid)
        {
            winner = candidate;
            continue;
        }

        if (exact_better(candidate, winner))
        {
            winner = candidate;
            continue;
        }

        if (exact_equal(candidate, winner))
        {
            if (normalized_better(candidate, winner) ||
                (normalized_equal(candidate, winner) &&
                 candidate.m < winner.m))
            {
                winner = candidate;
            }
        }
    }

    return winner;
}

static void print_witness(
    const char* name,
    const Witness& w)
{
    std::cout
        << name << ".m=" << w.m << '\n'
        << name << ".k=" << w.k << '\n'
        << name << ".t=" << w.t << '\n'
        << name << ".sign=" << w.sign << '\n';
}

int main()
{
    std::cout << "START EXPERIMENT 412\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    stats.prime_count =
        static_cast<u64>(primes.size());

    for (u64 r : primes)
    {
        std::vector<Witness> previous(
            M_MAX + 1
        );

        std::vector<Witness> current(
            M_MAX + 1
        );

        /*
            K grows incrementally.

            For each m, only the newly introduced k=K needs
            to be tested against the existing best witness.
        */

        for (u64 K = 1; K <= K_LIMIT; ++K)
        {
            ++stats.k_points;

            for (u64 m = 1; m <= M_MAX; ++m)
            {
                Witness best =
                    current[m];

                Witness plus;
                Witness minus;

                bool new_record = false;

                if (make_witness(
                        r,
                        m,
                        K,
                        +1,
                        plus))
                {
                    if (!best.valid ||
                        plus.k > best.k ||
                        (plus.k == best.k &&
                         plus.t < best.t))
                    {
                        best = plus;
                        new_record = true;
                    }
                }

                if (make_witness(
                        r,
                        m,
                        K,
                        -1,
                        minus))
                {
                    if (!best.valid ||
                        minus.k > best.k ||
                        (minus.k == best.k &&
                         minus.t < best.t))
                    {
                        best = minus;
                        new_record = true;
                    }
                }

                current[m] = best;

                /*
                    Store whether the current K created a new
                    record for this m.

                    We encode this in the sign of a separate
                    local array below by comparing k values
                    between previous/current winners.
                */
                (void)new_record;
            }

            const Witness normalized_old =
                choose_normalized_winner(previous);

            const Witness normalized_new =
                choose_normalized_winner(current);

            if (normalized_old.valid &&
                normalized_new.valid &&
                normalized_old.m != normalized_new.m)
            {
                bool new_record =
                    normalized_new.k > normalized_old.k;

                update_normalized_transition(
                    r,
                    K,
                    normalized_old,
                    normalized_new,
                    new_record,
                    stats
                );
            }

            const Witness exact_old =
                choose_exact_winner(previous);

            const Witness exact_new =
                choose_exact_winner(current);

            if (exact_old.valid &&
                exact_new.valid &&
                exact_old.m != exact_new.m)
            {
                bool new_record =
                    exact_new.k > exact_old.k;

                update_exact_transition(
                    r,
                    K,
                    exact_old,
                    exact_new,
                    new_record,
                    stats
                );
            }

            previous = current;
        }
    }

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT << '\n'
        << "K_LIMIT="
        << K_LIMIT << '\n'
        << "M_MAX="
        << M_MAX << '\n'
        << "PRIME_COUNT="
        << stats.prime_count << '\n'
        << '\n';

    std::cout
        << "K_POINTS="
        << stats.k_points << '\n'
        << '\n';

    std::cout
        << "NORMALIZED_TRANSITIONS="
        << stats.normalized_transitions << '\n'
        << "NORMALIZED_DELTA_ONE="
        << stats.normalized_delta_one << '\n'
        << "NORMALIZED_DELTA_GT_ONE="
        << stats.normalized_delta_gt_one << '\n'
        << "NORMALIZED_NEW_RECORD_TRANSITIONS="
        << stats.normalized_new_record_transitions << '\n'
        << "NORMALIZED_TRANSITIONS_WITHOUT_NEW_RECORD="
        << stats.normalized_transitions_without_new_record << '\n'
        << "NORMALIZED_FORWARD_TRANSITIONS="
        << stats.normalized_forward_transitions << '\n'
        << "NORMALIZED_BACKWARD_TRANSITIONS="
        << stats.normalized_backward_transitions << '\n'
        << '\n';

    std::cout
        << "EXACT_TRANSITIONS="
        << stats.exact_transitions << '\n'
        << "EXACT_DELTA_ONE="
        << stats.exact_delta_one << '\n'
        << "EXACT_DELTA_GT_ONE="
        << stats.exact_delta_gt_one << '\n'
        << "EXACT_NEW_RECORD_TRANSITIONS="
        << stats.exact_new_record_transitions << '\n'
        << "EXACT_TRANSITIONS_WITHOUT_NEW_RECORD="
        << stats.exact_transitions_without_new_record << '\n'
        << "EXACT_FORWARD_TRANSITIONS="
        << stats.exact_forward_transitions << '\n'
        << "EXACT_BACKWARD_TRANSITIONS="
        << stats.exact_backward_transitions << '\n'
        << '\n';

    std::cout
        << "NORMALIZED_DELTA_MIN="
        << stats.normalized_delta_min << '\n'
        << "NORMALIZED_DELTA_MAX="
        << stats.normalized_delta_max << '\n'
        << "NORMALIZED_MAX_K_SUM="
        << stats.normalized_max_k_sum << '\n'
        << '\n';

    std::cout
        << "EXACT_DELTA_MIN="
        << stats.exact_delta_min << '\n'
        << "EXACT_DELTA_MAX="
        << stats.exact_delta_max << '\n'
        << "EXACT_MAX_K_SUM="
        << stats.exact_max_k_sum << '\n'
        << '\n';

    if (stats.have_normalized_delta_counterexample)
    {
        std::cout
            << "FIRST_NORMALIZED_DELTA_GT_ONE\n"
            << "R="
            << stats.counter_r << '\n'
            << "K="
            << stats.counter_k << '\n';

        print_witness(
            "OLD",
            stats.counter_old
        );

        print_witness(
            "NEW",
            stats.counter_new
        );

        std::cout << '\n';
    }

    if (stats.have_exact_delta_counterexample)
    {
        std::cout
            << "FIRST_EXACT_DELTA_GT_ONE\n"
            << "R="
            << stats.exact_counter_r << '\n'
            << "K="
            << stats.exact_counter_k << '\n';

        print_witness(
            "OLD",
            stats.exact_counter_old
        );

        print_witness(
            "NEW",
            stats.exact_counter_new
        );

        std::cout << '\n';
    }

    std::cout
        << "NORMALIZED_DELTA_ONE_STATUS="
        << (
            stats.normalized_delta_gt_one == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "NORMALIZED_RECORD_STATUS="
        << (
            stats.normalized_transitions_without_new_record == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "EXACT_DELTA_ONE_STATUS="
        << (
            stats.exact_delta_gt_one == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "EXACT_RECORD_STATUS="
        << (
            stats.exact_transitions_without_new_record == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout << "FINISHED EXPERIMENT 412\n";

    return 0;
}
