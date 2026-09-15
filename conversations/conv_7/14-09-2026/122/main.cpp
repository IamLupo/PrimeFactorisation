#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

static constexpr u64 PRIME_LIMIT = 5000;
static constexpr u64 K_LIMIT     = 1000;
static constexpr u64 M_MAX       = 128;
static constexpr u64 SMALL_M     = 6;

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

    u64 global_winner_points = 0;

    u64 winner_m_1 = 0;
    u64 winner_m_2 = 0;
    u64 winner_m_3 = 0;
    u64 winner_m_4 = 0;
    u64 winner_m_5 = 0;
    u64 winner_m_6 = 0;

    u64 winner_m_ge7 = 0;

    u64 restricted_matches_global = 0;
    u64 restricted_mismatches_global = 0;

    u64 exact_restricted_matches_global = 0;
    u64 exact_restricted_mismatches_global = 0;

    u64 max_global_winner_m = 0;
    u64 max_global_winner_k = 0;
    u64 max_global_winner_r = 0;
    u64 max_global_winner_K = 0;

    u64 max_global_winner_t = 0;

    bool have_first_ge7 = false;
    Witness first_ge7;
    u64 first_ge7_r = 0;
    u64 first_ge7_K = 0;

    bool have_first_restricted_mismatch = false;
    Witness mismatch_global;
    Witness mismatch_restricted;
    u64 mismatch_r = 0;
    u64 mismatch_K = 0;

    bool have_first_exact_mismatch = false;
    Witness exact_mismatch_global;
    Witness exact_mismatch_restricted;
    u64 exact_mismatch_r = 0;
    u64 exact_mismatch_K = 0;

    u64 exact_global_points = 0;
    u64 exact_restricted_points = 0;
    u64 exact_same_winner = 0;
    u64 exact_winner_mismatch = 0;

    u64 m_ge7_normalized = 0;
    u64 m_ge7_exact = 0;
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

static void consider_exact(
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

    if (exact_better(candidate, winner))
    {
        winner = candidate;
        return;
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

static void record_winner_m(
    u64 m,
    Stats& stats)
{
    switch (m)
    {
        case 1:
            ++stats.winner_m_1;
            break;

        case 2:
            ++stats.winner_m_2;
            break;

        case 3:
            ++stats.winner_m_3;
            break;

        case 4:
            ++stats.winner_m_4;
            break;

        case 5:
            ++stats.winner_m_5;
            break;

        case 6:
            ++stats.winner_m_6;
            break;

        default:
            ++stats.winner_m_ge7;
            break;
    }
}

static void process_prime(
    u64 r,
    Stats& stats)
{
    /*
        Best witness for each m at the current K.
        Because K increases one at a time, only the newly
        introduced divisor k=K needs to be checked.
    */
    std::vector<Witness> best(
        M_MAX + 1
    );

    Witness global_normalized;
    Witness small_normalized;

    Witness global_exact;
    Witness small_exact;

    for (u64 K = 1; K <= K_LIMIT; ++K)
    {
        ++stats.k_points;
        ++stats.global_winner_points;

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

        global_normalized = Witness{};
        small_normalized = Witness{};

        global_exact = Witness{};
        small_exact = Witness{};

        for (u64 m = 1; m <= M_MAX; ++m)
        {
            if (!best[m].valid)
                continue;

            consider_normalized(
                best[m],
                global_normalized
            );

            consider_exact(
                best[m],
                global_exact
            );

            if (m <= SMALL_M)
            {
                consider_normalized(
                    best[m],
                    small_normalized
                );

                consider_exact(
                    best[m],
                    small_exact
                );
            }
        }

        if (global_normalized.valid)
        {
            record_winner_m(
                global_normalized.m,
                stats
            );

            if (global_normalized.m > stats.max_global_winner_m)
            {
                stats.max_global_winner_m =
                    global_normalized.m;

                stats.max_global_winner_k =
                    global_normalized.k;

                stats.max_global_winner_r = r;
                stats.max_global_winner_K = K;
                stats.max_global_winner_t =
                    global_normalized.t;
            }

            if (global_normalized.m >= 7)
            {
                ++stats.m_ge7_normalized;

                if (!stats.have_first_ge7)
                {
                    stats.have_first_ge7 = true;
                    stats.first_ge7 =
                        global_normalized;
                    stats.first_ge7_r = r;
                    stats.first_ge7_K = K;
                }
            }
        }

        if (global_exact.m >= 7 &&
            global_exact.valid)
        {
            ++stats.m_ge7_exact;
        }

        if (global_normalized.valid &&
            small_normalized.valid)
        {
            if (global_normalized.m ==
                    small_normalized.m &&
                global_normalized.k ==
                    small_normalized.k)
            {
                ++stats.restricted_matches_global;
            }
            else
            {
                ++stats.restricted_mismatches_global;

                if (!stats.have_first_restricted_mismatch)
                {
                    stats.have_first_restricted_mismatch = true;

                    stats.mismatch_global =
                        global_normalized;

                    stats.mismatch_restricted =
                        small_normalized;

                    stats.mismatch_r = r;
                    stats.mismatch_K = K;
                }
            }
        }

        if (global_exact.valid &&
            small_exact.valid)
        {
            ++stats.exact_global_points;
            ++stats.exact_restricted_points;

            if (global_exact.m ==
                    small_exact.m &&
                global_exact.k ==
                    small_exact.k)
            {
                ++stats.exact_restricted_matches_global;
            }
            else
            {
                ++stats.exact_restricted_mismatches_global;

                if (!stats.have_first_exact_mismatch)
                {
                    stats.have_first_exact_mismatch = true;

                    stats.exact_mismatch_global =
                        global_exact;

                    stats.exact_mismatch_restricted =
                        small_exact;

                    stats.exact_mismatch_r = r;
                    stats.exact_mismatch_K = K;
                }
            }
        }

        if (global_normalized.valid &&
            global_exact.valid)
        {
            if (global_normalized.m ==
                    global_exact.m &&
                global_normalized.k ==
                    global_exact.k)
            {
                ++stats.exact_same_winner;
            }
            else
            {
                ++stats.exact_winner_mismatch;
            }
        }
    }
}

int main()
{
    std::cout << "START EXPERIMENT 414\n";

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
        << "SMALL_M="
        << SMALL_M << '\n'
        << "PRIME_COUNT="
        << stats.prime_count << '\n'
        << '\n';

    std::cout
        << "GLOBAL_WINNER_POINTS="
        << stats.global_winner_points << '\n'
        << '\n';

    std::cout
        << "WINNER_M_1="
        << stats.winner_m_1 << '\n'
        << "WINNER_M_2="
        << stats.winner_m_2 << '\n'
        << "WINNER_M_3="
        << stats.winner_m_3 << '\n'
        << "WINNER_M_4="
        << stats.winner_m_4 << '\n'
        << "WINNER_M_5="
        << stats.winner_m_5 << '\n'
        << "WINNER_M_6="
        << stats.winner_m_6 << '\n'
        << "WINNER_M_GE7="
        << stats.winner_m_ge7 << '\n'
        << '\n';

    std::cout
        << "M_GE7_NORMALIZED="
        << stats.m_ge7_normalized << '\n'
        << "M_GE7_EXACT="
        << stats.m_ge7_exact << '\n'
        << '\n';

    std::cout
        << "RESTRICTED_MATCHES_GLOBAL="
        << stats.restricted_matches_global << '\n'
        << "RESTRICTED_MISMATCHES_GLOBAL="
        << stats.restricted_mismatches_global << '\n'
        << '\n';

    std::cout
        << "EXACT_RESTRICTED_MATCHES_GLOBAL="
        << stats.exact_restricted_matches_global << '\n'
        << "EXACT_RESTRICTED_MISMATCHES_GLOBAL="
        << stats.exact_restricted_mismatches_global << '\n'
        << '\n';

    std::cout
        << "GLOBAL_NORMALIZED_EXACT_SAME="
        << stats.exact_same_winner << '\n'
        << "GLOBAL_NORMALIZED_EXACT_MISMATCH="
        << stats.exact_winner_mismatch << '\n'
        << '\n';

    std::cout
        << "MAX_GLOBAL_WINNER_M="
        << stats.max_global_winner_m << '\n'
        << "MAX_GLOBAL_WINNER_K="
        << stats.max_global_winner_k << '\n'
        << "MAX_GLOBAL_WINNER_R="
        << stats.max_global_winner_r << '\n'
        << "MAX_GLOBAL_WINNER_K_POINT="
        << stats.max_global_winner_K << '\n'
        << "MAX_GLOBAL_WINNER_T="
        << stats.max_global_winner_t << '\n'
        << '\n';

    if (stats.have_first_ge7)
    {
        std::cout
            << "FIRST_M_GE7_WINNER\n"
            << "R="
            << stats.first_ge7_r << '\n'
            << "K="
            << stats.first_ge7_K << '\n';

        print_witness(
            "WINNER",
            stats.first_ge7
        );

        std::cout << '\n';
    }

    if (stats.have_first_restricted_mismatch)
    {
        std::cout
            << "FIRST_NORMALIZED_RESTRICTED_MISMATCH\n"
            << "R="
            << stats.mismatch_r << '\n'
            << "K="
            << stats.mismatch_K << '\n';

        print_witness(
            "GLOBAL",
            stats.mismatch_global
        );

        print_witness(
            "SMALL",
            stats.mismatch_restricted
        );

        std::cout << '\n';
    }

    if (stats.have_first_exact_mismatch)
    {
        std::cout
            << "FIRST_EXACT_RESTRICTED_MISMATCH\n"
            << "R="
            << stats.exact_mismatch_r << '\n'
            << "K="
            << stats.exact_mismatch_K << '\n';

        print_witness(
            "GLOBAL",
            stats.exact_mismatch_global
        );

        print_witness(
            "SMALL",
            stats.exact_mismatch_restricted
        );

        std::cout << '\n';
    }

    std::cout
        << "NO_M_GE7_STATUS="
        << (
            stats.winner_m_ge7 == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "M_LE6_GLOBAL_STATUS="
        << (
            stats.restricted_mismatches_global == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "EXACT_M_LE6_GLOBAL_STATUS="
        << (
            stats.exact_restricted_mismatches_global == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "NORMALIZED_EXACT_STATUS="
        << (
            stats.exact_winner_mismatch == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout << "FINISHED EXPERIMENT 414\n";

    return 0;
}
