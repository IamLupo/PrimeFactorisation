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

struct EdgeStats
{
    u64 count = 0;

    u64 min_r = std::numeric_limits<u64>::max();
    u64 max_r = 0;

    u64 min_K = std::numeric_limits<u64>::max();
    u64 max_K = 0;

    u64 min_k_old = std::numeric_limits<u64>::max();
    u64 max_k_old = 0;

    u64 min_k_new = std::numeric_limits<u64>::max();
    u64 max_k_new = 0;

    u64 min_delta = std::numeric_limits<u64>::max();
    u64 max_delta = 0;

    u64 min_abs_correction = std::numeric_limits<u64>::max();
    u64 max_abs_correction = 0;

    u64 min_threshold_gap = std::numeric_limits<u64>::max();
    u64 max_threshold_gap = 0;

    u64 sign_pp = 0;
    u64 sign_pm = 0;
    u64 sign_mp = 0;
    u64 sign_mm = 0;

    u64 gcd_k_min = std::numeric_limits<u64>::max();
    u64 gcd_k_max = 0;

    u64 gcd_m_min = std::numeric_limits<u64>::max();
    u64 gcd_m_max = 0;

    u64 gcd_delta_kprod_min =
        std::numeric_limits<u64>::max();

    u64 gcd_delta_kprod_max = 0;

    bool have_first = false;

    u64 first_r = 0;
    u64 first_K = 0;

    Witness first_old;
    Witness first_new;
};

struct Stats
{
    u64 prime_count = 0;
    u64 total_points = 0;
    u64 total_transitions = 0;

    EdgeStats edges[M_MAX + 1][M_MAX + 1];

    u64 edge_types = 0;
    u64 directed_edges = 0;

    u64 reverse_type_pairs = 0;
    u64 missing_reverse_edges = 0;

    u64 reciprocal_count_match = 0;
    u64 reciprocal_count_mismatch = 0;

    u64 threshold_valid_count = 0;
    u64 threshold_invalid_count = 0;

    u64 actual_r_gt_threshold = 0;
    u64 actual_r_eq_threshold = 0;
    u64 actual_r_lt_threshold = 0;

    u64 max_delta = 0;
    u64 max_delta_edge_old = 0;
    u64 max_delta_edge_new = 0;
    u64 max_delta_r = 0;
    u64 max_delta_K = 0;

    u64 max_threshold_gap = 0;
    u64 max_threshold_gap_r = 0;
    u64 max_threshold_gap_K = 0;
    u64 max_threshold_gap_old_m = 0;
    u64 max_threshold_gap_new_m = 0;

    bool have_first_missing_reverse = false;
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
    std::vector<Witness>& best)
{
    Witness candidate;

    if (make_witness(
            r,
            m,
            k,
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
            k,
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
    /*
        New winner has the larger normalized score:

            k2/m2 > k1/m1

        therefore

            Delta = m1*k2 - m2*k1 > 0.
    */

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
    /*
        t_old < t_new iff

            Delta*r > C

        with

            C = e_old*k_new - e_new*k_old.
    */

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

static void record_sign(
    const Witness& old_w,
    const Witness& new_w,
    EdgeStats& edge)
{
    if (old_w.sign == +1 &&
        new_w.sign == +1)
    {
        ++edge.sign_pp;
    }
    else if (old_w.sign == +1 &&
             new_w.sign == -1)
    {
        ++edge.sign_pm;
    }
    else if (old_w.sign == -1 &&
             new_w.sign == +1)
    {
        ++edge.sign_mp;
    }
    else
    {
        ++edge.sign_mm;
    }
}

static void record_transition(
    u64 r,
    u64 K,
    const Witness& old_w,
    const Witness& new_w,
    Stats& stats)
{
    const u64 m1 = old_w.m;
    const u64 m2 = new_w.m;

    if (m1 == m2)
        return;

    EdgeStats& edge =
        stats.edges[m1][m2];

    ++edge.count;
    ++stats.total_transitions;

    if (!edge.have_first)
    {
        edge.have_first = true;
        edge.first_r = r;
        edge.first_K = K;
        edge.first_old = old_w;
        edge.first_new = new_w;
    }

    if (r < edge.min_r)
        edge.min_r = r;

    if (r > edge.max_r)
        edge.max_r = r;

    if (K < edge.min_K)
        edge.min_K = K;

    if (K > edge.max_K)
        edge.max_K = K;

    if (old_w.k < edge.min_k_old)
        edge.min_k_old = old_w.k;

    if (old_w.k > edge.max_k_old)
        edge.max_k_old = old_w.k;

    if (new_w.k < edge.min_k_new)
        edge.min_k_new = new_w.k;

    if (new_w.k > edge.max_k_new)
        edge.max_k_new = new_w.k;

    const u64 delta =
        compute_delta(
            old_w,
            new_w
        );

    const i128 correction =
        compute_correction(
            old_w,
            new_w
        );

    const u64 threshold =
        compute_threshold(
            delta,
            correction
        );

    if (delta < edge.min_delta)
        edge.min_delta = delta;

    if (delta > edge.max_delta)
        edge.max_delta = delta;

    u64 abs_correction = 0;

    if (correction >= 0)
    {
        abs_correction =
            static_cast<u64>(correction);
    }
    else
    {
        abs_correction =
            static_cast<u64>(-correction);
    }

    if (abs_correction <
        edge.min_abs_correction)
    {
        edge.min_abs_correction =
            abs_correction;
    }

    if (abs_correction >
        edge.max_abs_correction)
    {
        edge.max_abs_correction =
            abs_correction;
    }

    const u64 gcd_k =
        std::gcd(
            old_w.k,
            new_w.k
        );

    const u64 gcd_m =
        std::gcd(
            old_w.m,
            new_w.m
        );

    const u128 k_product_128 =
        static_cast<u128>(old_w.k) *
        new_w.k;

    const u64 k_product =
        static_cast<u64>(k_product_128);

    const u64 gcd_delta_kprod =
        std::gcd(
            delta,
            k_product
        );

    if (gcd_k < edge.gcd_k_min)
        edge.gcd_k_min = gcd_k;

    if (gcd_k > edge.gcd_k_max)
        edge.gcd_k_max = gcd_k;

    if (gcd_m < edge.gcd_m_min)
        edge.gcd_m_min = gcd_m;

    if (gcd_m > edge.gcd_m_max)
        edge.gcd_m_max = gcd_m;

    if (gcd_delta_kprod <
        edge.gcd_delta_kprod_min)
    {
        edge.gcd_delta_kprod_min =
            gcd_delta_kprod;
    }

    if (gcd_delta_kprod >
        edge.gcd_delta_kprod_max)
    {
        edge.gcd_delta_kprod_max =
            gcd_delta_kprod;
    }

    u64 threshold_gap = 0;

    if (r >= threshold)
    {
        threshold_gap = r - threshold;

        ++stats.actual_r_gt_threshold;
    }
    else
    {
        threshold_gap = threshold - r;

        ++stats.actual_r_lt_threshold;
    }

    if (r == threshold)
        ++stats.actual_r_eq_threshold;

    if (threshold_gap <
        edge.min_threshold_gap)
    {
        edge.min_threshold_gap =
            threshold_gap;
    }

    if (threshold_gap >
        edge.max_threshold_gap)
    {
        edge.max_threshold_gap =
            threshold_gap;
    }

    if (threshold > 0)
        ++stats.threshold_valid_count;
    else
        ++stats.threshold_invalid_count;

    if (delta > stats.max_delta)
    {
        stats.max_delta = delta;
        stats.max_delta_edge_old = m1;
        stats.max_delta_edge_new = m2;
        stats.max_delta_r = r;
        stats.max_delta_K = K;
    }

    if (threshold_gap >
        stats.max_threshold_gap)
    {
        stats.max_threshold_gap =
            threshold_gap;

        stats.max_threshold_gap_r = r;
        stats.max_threshold_gap_K = K;
        stats.max_threshold_gap_old_m = m1;
        stats.max_threshold_gap_new_m = m2;
    }

    record_sign(
        old_w,
        new_w,
        edge
    );
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

static void print_edge_summary(
    u64 m1,
    u64 m2,
    const EdgeStats& edge)
{
    std::cout
        << "EDGE "
        << m1 << "->" << m2
        << " COUNT="
        << edge.count
        << " R_RANGE=["
        << edge.min_r
        << ","
        << edge.max_r
        << "]"
        << " K_RANGE=["
        << edge.min_K
        << ","
        << edge.max_K
        << "]"
        << " KOLD_RANGE=["
        << edge.min_k_old
        << ","
        << edge.max_k_old
        << "]"
        << " KNEW_RANGE=["
        << edge.min_k_new
        << ","
        << edge.max_k_new
        << "]"
        << " DELTA_RANGE=["
        << edge.min_delta
        << ","
        << edge.max_delta
        << "]"
        << " C_RANGE=["
        << edge.min_abs_correction
        << ","
        << edge.max_abs_correction
        << "]"
        << " GAP_RANGE=["
        << edge.min_threshold_gap
        << ","
        << edge.max_threshold_gap
        << "]"
        << " SIGN=("
        << edge.sign_pp << ","
        << edge.sign_pm << ","
        << edge.sign_mp << ","
        << edge.sign_mm << ")"
        << " GCDK=["
        << edge.gcd_k_min
        << ","
        << edge.gcd_k_max
        << "]"
        << " GCDM=["
        << edge.gcd_m_min
        << ","
        << edge.gcd_m_max
        << "]"
        << " GCDDC=["
        << edge.gcd_delta_kprod_min
        << ","
        << edge.gcd_delta_kprod_max
        << "]"
        << '\n';

    if (edge.have_first)
    {
        std::cout
            << "  FIRST_R="
            << edge.first_r
            << " FIRST_K="
            << edge.first_K
            << '\n';

        print_witness(
            "  OLD",
            edge.first_old
        );

        print_witness(
            "  NEW",
            edge.first_new
        );
    }
}

int main()
{
    std::cout
        << "START EXPERIMENT 416\n";

    const std::vector<u64> primes =
        generate_primes(
            PRIME_LIMIT
        );

    Stats stats;

    stats.prime_count =
        static_cast<u64>(
            primes.size()
        );

    std::vector<Witness> best(
        M_MAX + 1
    );

    for (u64 r : primes)
    {
        for (u64 m = 0; m <= M_MAX; ++m)
            best[m] = Witness{};

        Witness previous_winner;
        bool have_previous = false;

        for (u64 K = 1;
             K <= K_LIMIT;
             ++K)
        {
            ++stats.total_points;

            /*
                Because K advances by one, only the newly
                introduced divisor k=K can improve a given m.
            */
            for (u64 m = 1;
                 m <= M_MAX;
                 ++m)
            {
                update_best_for_m(
                    r,
                    m,
                    K,
                    best
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

            if (have_previous &&
                previous_winner.valid &&
                previous_winner.m != winner.m)
            {
                record_transition(
                    r,
                    K,
                    previous_winner,
                    winner,
                    stats
                );
            }

            previous_winner = winner;
            have_previous = true;
        }
    }

    /*
        Count directed edge types and test whether every edge
        has a reverse edge.
    */

    for (u64 m1 = 1;
         m1 <= M_MAX;
         ++m1)
    {
        for (u64 m2 = 1;
             m2 <= M_MAX;
             ++m2)
        {
            if (stats.edges[m1][m2].count == 0)
                continue;

            ++stats.directed_edges;
            ++stats.edge_types;
        }
    }

    for (u64 m1 = 1;
         m1 <= M_MAX;
         ++m1)
    {
        for (u64 m2 = m1 + 1;
             m2 <= M_MAX;
             ++m2)
        {
            const u64 forward =
                stats.edges[m1][m2].count;

            const u64 reverse =
                stats.edges[m2][m1].count;

            if (forward == 0 &&
                reverse == 0)
            {
                continue;
            }

            ++stats.reverse_type_pairs;

            if (forward == 0 ||
                reverse == 0)
            {
                ++stats.missing_reverse_edges;

                if (!stats.have_first_missing_reverse)
                {
                    stats.have_first_missing_reverse =
                        true;
                }
            }
            else if (forward == reverse)
            {
                ++stats.reciprocal_count_match;
            }
            else
            {
                ++stats.reciprocal_count_mismatch;
            }
        }
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
        << "TOTAL_TRANSITIONS="
        << stats.total_transitions
        << '\n'
        << '\n';

    std::cout
        << "DIRECTED_EDGE_TYPES="
        << stats.directed_edges
        << '\n'
        << "REVERSE_TYPE_PAIRS="
        << stats.reverse_type_pairs
        << '\n'
        << "MISSING_REVERSE_EDGES="
        << stats.missing_reverse_edges
        << '\n'
        << "RECIPROCAL_COUNT_MATCH="
        << stats.reciprocal_count_match
        << '\n'
        << "RECIPROCAL_COUNT_MISMATCH="
        << stats.reciprocal_count_mismatch
        << '\n'
        << '\n';

    std::cout
        << "THRESHOLD_VALID_COUNT="
        << stats.threshold_valid_count
        << '\n'
        << "THRESHOLD_INVALID_COUNT="
        << stats.threshold_invalid_count
        << '\n'
        << "ACTUAL_R_GT_THRESHOLD="
        << stats.actual_r_gt_threshold
        << '\n'
        << "ACTUAL_R_EQ_THRESHOLD="
        << stats.actual_r_eq_threshold
        << '\n'
        << "ACTUAL_R_LT_THRESHOLD="
        << stats.actual_r_lt_threshold
        << '\n'
        << '\n';

    std::cout
        << "MAX_DELTA="
        << stats.max_delta
        << '\n'
        << "MAX_DELTA_EDGE="
        << stats.max_delta_edge_old
        << "->"
        << stats.max_delta_edge_new
        << '\n'
        << "MAX_DELTA_R="
        << stats.max_delta_r
        << '\n'
        << "MAX_DELTA_K="
        << stats.max_delta_K
        << '\n'
        << '\n';

    std::cout
        << "MAX_THRESHOLD_GAP="
        << stats.max_threshold_gap
        << '\n'
        << "MAX_THRESHOLD_GAP_R="
        << stats.max_threshold_gap_r
        << '\n'
        << "MAX_THRESHOLD_GAP_K="
        << stats.max_threshold_gap_K
        << '\n'
        << "MAX_THRESHOLD_GAP_OLD_M="
        << stats.max_threshold_gap_old_m
        << '\n'
        << "MAX_THRESHOLD_GAP_NEW_M="
        << stats.max_threshold_gap_new_m
        << '\n'
        << '\n';

    std::cout
        << "EDGE_SUMMARY\n";

    for (u64 m1 = 1;
         m1 <= M_MAX;
         ++m1)
    {
        for (u64 m2 = 1;
             m2 <= M_MAX;
             ++m2)
        {
            if (stats.edges[m1][m2].count == 0)
                continue;

            print_edge_summary(
                m1,
                m2,
                stats.edges[m1][m2]
            );
        }
    }

    std::cout
        << '\n'
        << "REVERSE_EDGE_STATUS="
        << (
            stats.missing_reverse_edges == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "RECIPROCAL_COUNT_STATUS="
        << (
            stats.reciprocal_count_mismatch == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 416\n";

    return 0;
}
