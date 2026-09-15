#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
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

struct Transition
{
    Witness old_w;
    Witness new_w;

    u64 r = 0;
    u64 K = 0;

    u64 delta = 0;

    u64 gcd_m = 0;
    u64 gcd_k = 0;
    u64 gcd_delta_denominator = 0;

    u64 normalized_jump_num = 0;
    u64 normalized_jump_den = 1;

    u64 exact_drop = 0;

    bool new_divisor_record = false;
};

struct Stats
{
    u64 prime_count = 0;
    u64 k_points = 0;

    u64 normalized_transitions = 0;
    u64 exact_transitions = 0;

    u64 normalized_exact_same_transition = 0;
    u64 normalized_exact_transition_mismatch = 0;

    u64 transitions_with_new_record = 0;
    u64 transitions_without_new_record = 0;

    u64 delta_min = std::numeric_limits<u64>::max();
    u64 delta_max = 0;

    u64 gcd_m_min = std::numeric_limits<u64>::max();
    u64 gcd_m_max = 0;

    u64 gcd_k_min = std::numeric_limits<u64>::max();
    u64 gcd_k_max = 0;

    u64 gcd_delta_denominator_min =
        std::numeric_limits<u64>::max();
    u64 gcd_delta_denominator_max = 0;

    u64 jump_num_min = std::numeric_limits<u64>::max();
    u64 jump_num_max = 0;

    u64 jump_den_min = std::numeric_limits<u64>::max();
    u64 jump_den_max = 0;

    u64 exact_drop_min =
        std::numeric_limits<u64>::max();
    u64 exact_drop_max = 0;

    u64 sign_pp = 0;
    u64 sign_pm = 0;
    u64 sign_mp = 0;
    u64 sign_mm = 0;

    u64 old_m_min = std::numeric_limits<u64>::max();
    u64 old_m_max = 0;

    u64 new_m_min = std::numeric_limits<u64>::max();
    u64 new_m_max = 0;

    u64 old_k_min = std::numeric_limits<u64>::max();
    u64 old_k_max = 0;

    u64 new_k_min = std::numeric_limits<u64>::max();
    u64 new_k_max = 0;

    u64 max_jump_delta = 0;
    bool have_max_jump = false;
    Transition max_jump_transition;

    bool have_first_transition = false;
    Transition first_transition;

    bool have_first_without_record = false;
    Transition first_without_record;
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

static Witness best_for_m(
    u64 r,
    u64 m,
    u64 K)
{
    Witness best;

    for (u64 k = 1; k <= K; ++k)
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

    return best;
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

static u64 transition_delta(
    const Witness& old_w,
    const Witness& new_w)
{
    /*
        New normalized winner is larger:

            k2/m2 > k1/m1

        Therefore

            Delta = m1*k2 - m2*k1 > 0.
    */

    const i128 delta =
        static_cast<i128>(old_w.m) * new_w.k -
        static_cast<i128>(new_w.m) * old_w.k;

    if (delta <= 0)
        return 0;

    return static_cast<u64>(delta);
}

static bool is_new_divisor_record(
    const Witness& previous,
    const Witness& current)
{
    /*
        For the winner transition, a new record means that
        the newly winning witness has a strictly larger k
        than the previous winner.

        Since K increases by one at each step, this means the
        current winning divisor was introduced at this K.
    */

    return current.k > previous.k;
}

static Transition build_transition(
    u64 r,
    u64 K,
    const Witness& old_w,
    const Witness& new_w)
{
    Transition tr;

    tr.old_w = old_w;
    tr.new_w = new_w;
    tr.r = r;
    tr.K = K;

    tr.delta =
        transition_delta(
            old_w,
            new_w
        );

    tr.gcd_m =
        std::gcd(
            old_w.m,
            new_w.m
        );

    tr.gcd_k =
        std::gcd(
            old_w.k,
            new_w.k
        );

    const u128 denominator =
        static_cast<u128>(old_w.m) *
        new_w.m;

    const u64 denominator64 =
        static_cast<u64>(denominator);

    tr.gcd_delta_denominator =
        std::gcd(
            tr.delta,
            denominator64
        );

    /*
        Normalized jump:

            k2/m2 - k1/m1
              = Delta/(m1*m2).
    */

    tr.normalized_jump_num = tr.delta;
    tr.normalized_jump_den = denominator64;

    const i128 exact_drop =
        static_cast<i128>(old_w.t) -
        static_cast<i128>(new_w.t);

    if (exact_drop > 0 &&
        exact_drop <=
            std::numeric_limits<u64>::max())
    {
        tr.exact_drop =
            static_cast<u64>(exact_drop);
    }

    tr.new_divisor_record =
        is_new_divisor_record(
            old_w,
            new_w
        );

    return tr;
}

static void update_sign_stats(
    const Transition& tr,
    Stats& stats)
{
    const int s1 = tr.old_w.sign;
    const int s2 = tr.new_w.sign;

    if (s1 == +1 && s2 == +1)
        ++stats.sign_pp;
    else if (s1 == +1 && s2 == -1)
        ++stats.sign_pm;
    else if (s1 == -1 && s2 == +1)
        ++stats.sign_mp;
    else
        ++stats.sign_mm;
}

static void update_stats(
    const Transition& tr,
    Stats& stats)
{
    ++stats.normalized_transitions;

    if (tr.new_divisor_record)
        ++stats.transitions_with_new_record;
    else
        ++stats.transitions_without_new_record;

    if (tr.delta < stats.delta_min)
        stats.delta_min = tr.delta;

    if (tr.delta > stats.delta_max)
        stats.delta_max = tr.delta;

    if (tr.gcd_m < stats.gcd_m_min)
        stats.gcd_m_min = tr.gcd_m;

    if (tr.gcd_m > stats.gcd_m_max)
        stats.gcd_m_max = tr.gcd_m;

    if (tr.gcd_k < stats.gcd_k_min)
        stats.gcd_k_min = tr.gcd_k;

    if (tr.gcd_k > stats.gcd_k_max)
        stats.gcd_k_max = tr.gcd_k;

    if (tr.gcd_delta_denominator <
        stats.gcd_delta_denominator_min)
    {
        stats.gcd_delta_denominator_min =
            tr.gcd_delta_denominator;
    }

    if (tr.gcd_delta_denominator >
        stats.gcd_delta_denominator_max)
    {
        stats.gcd_delta_denominator_max =
            tr.gcd_delta_denominator;
    }

    if (tr.normalized_jump_num <
        stats.jump_num_min)
    {
        stats.jump_num_min =
            tr.normalized_jump_num;
    }

    if (tr.normalized_jump_num >
        stats.jump_num_max)
    {
        stats.jump_num_max =
            tr.normalized_jump_num;
    }

    if (tr.normalized_jump_den <
        stats.jump_den_min)
    {
        stats.jump_den_min =
            tr.normalized_jump_den;
    }

    if (tr.normalized_jump_den >
        stats.jump_den_max)
    {
        stats.jump_den_max =
            tr.normalized_jump_den;
    }

    if (tr.exact_drop <
        stats.exact_drop_min)
    {
        stats.exact_drop_min =
            tr.exact_drop;
    }

    if (tr.exact_drop >
        stats.exact_drop_max)
    {
        stats.exact_drop_max =
            tr.exact_drop;
    }

    if (tr.old_w.m < stats.old_m_min)
        stats.old_m_min = tr.old_w.m;

    if (tr.old_w.m > stats.old_m_max)
        stats.old_m_max = tr.old_w.m;

    if (tr.new_w.m < stats.new_m_min)
        stats.new_m_min = tr.new_w.m;

    if (tr.new_w.m > stats.new_m_max)
        stats.new_m_max = tr.new_w.m;

    if (tr.old_w.k < stats.old_k_min)
        stats.old_k_min = tr.old_w.k;

    if (tr.old_w.k > stats.old_k_max)
        stats.old_k_max = tr.old_w.k;

    if (tr.new_w.k < stats.new_k_min)
        stats.new_k_min = tr.new_w.k;

    if (tr.new_w.k > stats.new_k_max)
        stats.new_k_max = tr.new_w.k;

    update_sign_stats(tr, stats);

    if (!stats.have_first_transition)
    {
        stats.have_first_transition = true;
        stats.first_transition = tr;
    }

    if (!tr.new_divisor_record &&
        !stats.have_first_without_record)
    {
        stats.have_first_without_record = true;
        stats.first_without_record = tr;
    }

    if (!stats.have_max_jump ||
        tr.delta > stats.max_jump_delta)
    {
        stats.have_max_jump = true;
        stats.max_jump_delta = tr.delta;
        stats.max_jump_transition = tr;
    }
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

static void print_transition(
    const char* prefix,
    const Transition& tr)
{
    std::cout
        << prefix << ".r="
        << tr.r << '\n'
        << prefix << ".K="
        << tr.K << '\n'
        << prefix << ".delta="
        << tr.delta << '\n'
        << prefix << ".gcd_m="
        << tr.gcd_m << '\n'
        << prefix << ".gcd_k="
        << tr.gcd_k << '\n'
        << prefix << ".gcd_delta_denominator="
        << tr.gcd_delta_denominator << '\n'
        << prefix << ".jump_num="
        << tr.normalized_jump_num << '\n'
        << prefix << ".jump_den="
        << tr.normalized_jump_den << '\n'
        << prefix << ".exact_drop="
        << tr.exact_drop << '\n'
        << prefix << ".new_divisor_record="
        << (tr.new_divisor_record ? 1 : 0)
        << '\n';

    print_witness(
        prefix,
        tr.old_w
    );

    print_witness(
        "NEW",
        tr.new_w
    );
}

int main()
{
    std::cout << "START EXPERIMENT 413\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    stats.prime_count =
        static_cast<u64>(primes.size());

    /*
        Reuse fixed-size vectors.

        This avoids the vector assignment pattern that caused
        the GCC 15 optimizer warning in Experiment 412.
    */

    std::vector<Witness> current(
        M_MAX + 1
    );

    std::vector<Witness> previous(
        M_MAX + 1
    );

    for (u64 r : primes)
    {
        for (u64 m = 0; m <= M_MAX; ++m)
        {
            current[m] = Witness{};
            previous[m] = Witness{};
        }

        bool have_previous = false;

        for (u64 K = 1; K <= K_LIMIT; ++K)
        {
            ++stats.k_points;

            for (u64 m = 1; m <= M_MAX; ++m)
            {
                current[m] =
                    best_for_m(
                        r,
                        m,
                        K
                    );
            }

            const Witness normalized_current =
                choose_normalized_winner(
                    current
                );

            const Witness exact_current =
                choose_exact_winner(
                    current
                );

            if (have_previous)
            {
                const Witness normalized_previous =
                    choose_normalized_winner(
                        previous
                    );

                const Witness exact_previous =
                    choose_exact_winner(
                        previous
                    );

                /*
                    Normalized winner changed.
                */

                if (normalized_previous.valid &&
                    normalized_current.valid &&
                    normalized_previous.m !=
                        normalized_current.m)
                {
                    Transition tr =
                        build_transition(
                            r,
                            K,
                            normalized_previous,
                            normalized_current
                        );

                    update_stats(
                        tr,
                        stats
                    );
                }

                /*
                    Exact winner changed.

                    We separately compare the two winner
                    identities to see whether exact and
                    normalized envelopes take the same path.
                */

                if (exact_previous.valid &&
                    exact_current.valid &&
                    exact_previous.m !=
                        exact_current.m)
                {
                    ++stats.exact_transitions;

                    if (normalized_previous.valid &&
                        normalized_current.valid &&
                        exact_previous.m ==
                            normalized_previous.m &&
                        exact_current.m ==
                            normalized_current.m)
                    {
                        ++stats.normalized_exact_same_transition;
                    }
                    else
                    {
                        ++stats.normalized_exact_transition_mismatch;
                    }
                }
            }

            /*
                Swap contents manually to avoid vector
                assignment and the GCC warning observed in 412.
            */

            for (u64 m = 1; m <= M_MAX; ++m)
                previous[m] = current[m];

            have_previous = true;
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
        << "EXACT_TRANSITIONS="
        << stats.exact_transitions << '\n'
        << "NORMALIZED_EXACT_SAME_TRANSITION="
        << stats.normalized_exact_same_transition
        << '\n'
        << "NORMALIZED_EXACT_TRANSITION_MISMATCH="
        << stats.normalized_exact_transition_mismatch
        << '\n'
        << '\n';

    std::cout
        << "TRANSITIONS_WITH_NEW_RECORD="
        << stats.transitions_with_new_record << '\n'
        << "TRANSITIONS_WITHOUT_NEW_RECORD="
        << stats.transitions_without_new_record << '\n'
        << '\n';

    std::cout
        << "DELTA_MIN="
        << stats.delta_min << '\n'
        << "DELTA_MAX="
        << stats.delta_max << '\n'
        << '\n';

    std::cout
        << "GCD_M_MIN="
        << stats.gcd_m_min << '\n'
        << "GCD_M_MAX="
        << stats.gcd_m_max << '\n'
        << "GCD_K_MIN="
        << stats.gcd_k_min << '\n'
        << "GCD_K_MAX="
        << stats.gcd_k_max << '\n'
        << "GCD_DELTA_DENOMINATOR_MIN="
        << stats.gcd_delta_denominator_min << '\n'
        << "GCD_DELTA_DENOMINATOR_MAX="
        << stats.gcd_delta_denominator_max << '\n'
        << '\n';

    std::cout
        << "JUMP_NUM_MIN="
        << stats.jump_num_min << '\n'
        << "JUMP_NUM_MAX="
        << stats.jump_num_max << '\n'
        << "JUMP_DEN_MIN="
        << stats.jump_den_min << '\n'
        << "JUMP_DEN_MAX="
        << stats.jump_den_max << '\n'
        << '\n';

    std::cout
        << "EXACT_DROP_MIN="
        << stats.exact_drop_min << '\n'
        << "EXACT_DROP_MAX="
        << stats.exact_drop_max << '\n'
        << '\n';

    std::cout
        << "OLD_M_MIN="
        << stats.old_m_min << '\n'
        << "OLD_M_MAX="
        << stats.old_m_max << '\n'
        << "NEW_M_MIN="
        << stats.new_m_min << '\n'
        << "NEW_M_MAX="
        << stats.new_m_max << '\n'
        << "OLD_K_MIN="
        << stats.old_k_min << '\n'
        << "OLD_K_MAX="
        << stats.old_k_max << '\n'
        << "NEW_K_MIN="
        << stats.new_k_min << '\n'
        << "NEW_K_MAX="
        << stats.new_k_max << '\n'
        << '\n';

    std::cout
        << "SIGN_PP="
        << stats.sign_pp << '\n'
        << "SIGN_PM="
        << stats.sign_pm << '\n'
        << "SIGN_MP="
        << stats.sign_mp << '\n'
        << "SIGN_MM="
        << stats.sign_mm << '\n'
        << '\n';

    if (stats.have_first_transition)
    {
        std::cout
            << "FIRST_TRANSITION\n";

        print_transition(
            "TR",
            stats.first_transition
        );

        std::cout << '\n';
    }

    if (stats.have_max_jump)
    {
        std::cout
            << "MAX_DELTA_TRANSITION\n";

        print_transition(
            "TR",
            stats.max_jump_transition
        );

        std::cout << '\n';
    }

    if (stats.have_first_without_record)
    {
        std::cout
            << "FIRST_TRANSITION_WITHOUT_RECORD\n";

        print_transition(
            "TR",
            stats.first_without_record
        );

        std::cout << '\n';
    }

    std::cout
        << "ALL_TRANSITIONS_ON_RECORD="
        << (
            stats.transitions_without_new_record == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "NORMALIZED_EXACT_PATH_STATUS="
        << (
            stats.normalized_exact_transition_mismatch == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout << "FINISHED EXPERIMENT 413\n";

    return 0;
}
