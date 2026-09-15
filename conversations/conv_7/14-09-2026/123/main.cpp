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

struct FirstWin
{
    bool found = false;
    u64 K = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0;
    u64 r = 0;
};

struct TransitionCount
{
    u64 count = 0;
};

struct Stats
{
    u64 prime_count = 0;
    u64 total_points = 0;

    u64 winners[33] = {};

    FirstWin first_win[33];

    TransitionCount transitions[33][33];

    u64 max_first_K[33] = {};
    u64 min_first_K[33] = {};

    u64 first_m_ge7_r = 0;
    u64 first_m_ge7_K = 0;
    Witness first_m_ge7;

    bool have_m_ge7 = false;

    u64 max_winner_m = 0;
    u64 max_winner_r = 0;
    u64 max_winner_K = 0;
    u64 max_winner_k = 0;
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

static void consider_normalized(
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

static void record_first_win(
    const Witness& winner,
    u64 r,
    u64 K,
    Stats& stats)
{
    const u64 m = winner.m;

    ++stats.winners[m];

    if (!stats.first_win[m].found)
    {
        stats.first_win[m].found = true;
        stats.first_win[m].K = K;
        stats.first_win[m].k = winner.k;
        stats.first_win[m].t = winner.t;
        stats.first_win[m].sign = winner.sign;
        stats.first_win[m].r = r;
    }

    if (K > stats.max_first_K[m])
        stats.max_first_K[m] = K;

    if (stats.min_first_K[m] == 0 ||
        K < stats.min_first_K[m])
    {
        stats.min_first_K[m] = K;
    }

    if (m > stats.max_winner_m)
    {
        stats.max_winner_m = m;
        stats.max_winner_r = r;
        stats.max_winner_K = K;
        stats.max_winner_k = winner.k;
    }

    if (m >= 7 &&
        !stats.have_m_ge7)
    {
        stats.have_m_ge7 = true;
        stats.first_m_ge7_r = r;
        stats.first_m_ge7_K = K;
        stats.first_m_ge7 = winner;
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

static void process_prime(
    u64 r,
    Stats& stats)
{
    std::vector<Witness> best(
        M_MAX + 1
    );

    Witness previous_winner;
    bool have_previous = false;

    for (u64 K = 1; K <= K_LIMIT; ++K)
    {
        ++stats.total_points;

        for (u64 m = 1; m <= M_MAX; ++m)
        {
            Witness candidate;

            if (make_witness(
                    r,
                    m,
                    K,
                    +1,
                    candidate))
            {
                if (!best[m].valid ||
                    candidate.k > best[m].k ||
                    (candidate.k == best[m].k &&
                     candidate.t < best[m].t))
                {
                    best[m] = candidate;
                }
            }

            if (make_witness(
                    r,
                    m,
                    K,
                    -1,
                    candidate))
            {
                if (!best[m].valid ||
                    candidate.k > best[m].k ||
                    (candidate.k == best[m].k &&
                     candidate.t < best[m].t))
                {
                    best[m] = candidate;
                }
            }
        }

        Witness winner;

        for (u64 m = 1; m <= M_MAX; ++m)
        {
            if (!best[m].valid)
                continue;

            consider_normalized(
                best[m],
                winner
            );
        }

        if (!winner.valid)
            continue;

        record_first_win(
            winner,
            r,
            K,
            stats
        );

        if (have_previous &&
            previous_winner.valid &&
            previous_winner.m != winner.m)
        {
            ++stats.transitions
                [previous_winner.m]
                [winner.m]
                .count;
        }

        previous_winner = winner;
        have_previous = true;
    }
}

int main()
{
    std::cout << "START EXPERIMENT 415\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    stats.prime_count =
        static_cast<u64>(primes.size());

    for (u64 r : primes)
    {
        process_prime(
            r,
            stats
        );
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
        << "TOTAL_POINTS="
        << stats.total_points << '\n'
        << '\n';

    std::cout
        << "WINNER_COUNTS\n";

    for (u64 m = 1; m <= M_MAX; ++m)
    {
        std::cout
            << "M_"
            << m
            << "="
            << stats.winners[m]
            << '\n';
    }

    std::cout << '\n';

    std::cout
        << "FIRST_WIN_BY_M\n";

    for (u64 m = 1; m <= M_MAX; ++m)
    {
        if (!stats.first_win[m].found)
        {
            std::cout
                << "M_"
                << m
                << "=NEVER\n";

            continue;
        }

        const FirstWin& f =
            stats.first_win[m];

        std::cout
            << "M_"
            << m
            << ".FIRST_K="
            << f.K << '\n'
            << "M_"
            << m
            << ".R="
            << f.r << '\n'
            << "M_"
            << m
            << ".K="
            << f.K << '\n'
            << "M_"
            << m
            << ".KDIV="
            << f.k << '\n'
            << "M_"
            << m
            << ".T="
            << f.t << '\n'
            << "M_"
            << m
            << ".SIGN="
            << f.sign << '\n';
    }

    std::cout << '\n';

    std::cout
        << "TRANSITION_MATRIX\n";

    for (u64 m1 = 1; m1 <= M_MAX; ++m1)
    {
        for (u64 m2 = 1; m2 <= M_MAX; ++m2)
        {
            if (stats.transitions[m1][m2].count == 0)
                continue;

            std::cout
                << m1
                << "->"
                << m2
                << "="
                << stats.transitions[m1][m2].count
                << '\n';
        }
    }

    std::cout << '\n';

    std::cout
        << "MAX_FIRST_K_BY_M\n";

    for (u64 m = 1; m <= M_MAX; ++m)
    {
        if (stats.first_win[m].found)
        {
            std::cout
                << "M_"
                << m
                << ".MIN_FIRST_K="
                << stats.min_first_K[m]
                << '\n'
                << "M_"
                << m
                << ".MAX_FIRST_K="
                << stats.max_first_K[m]
                << '\n';
        }
    }

    std::cout << '\n';

    std::cout
        << "MAX_WINNER_M="
        << stats.max_winner_m
        << '\n'
        << "MAX_WINNER_R="
        << stats.max_winner_r
        << '\n'
        << "MAX_WINNER_K="
        << stats.max_winner_K
        << '\n'
        << "MAX_WINNER_DIVISOR="
        << stats.max_winner_k
        << '\n'
        << '\n';

    if (stats.have_m_ge7)
    {
        std::cout
            << "FIRST_M_GE7_WINNER\n"
            << "R="
            << stats.first_m_ge7_r
            << '\n'
            << "K="
            << stats.first_m_ge7_K
            << '\n';

        print_witness(
            "WINNER",
            stats.first_m_ge7
        );

        std::cout << '\n';
    }
    else
    {
        std::cout
            << "FIRST_M_GE7_WINNER=NEVER\n\n";
    }

    std::cout
        << "M_GE7_STATUS="
        << (
            stats.have_m_ge7
                ? "FOUND"
                : "NOT_FOUND"
        )
        << '\n';

    std::cout << "FINISHED EXPERIMENT 415\n";

    return 0;
}
