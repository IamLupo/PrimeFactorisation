#include <gmpxx.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

using u64 = std::uint64_t;

static constexpr u64 PRIME_LIMIT = 5000;
static constexpr u64 K_LIMIT     = 500;
static constexpr u64 M_MAX       = 64;

struct Witness {
    u64 m = 0;
    u64 k = 0;
    u64 t = 0;
    int sign = 0; // +1 or -1
    bool valid = false;
};

struct Stats {
    u64 pair_tests = 0;

    u64 normalized_first = 0;
    u64 normalized_second = 0;

    u64 exact_agree = 0;
    u64 exact_disagree = 0;

    u64 threshold_tests = 0;
    u64 threshold_bound_failures = 0;

    u64 global_bound_tests = 0;
    u64 global_bound_failures = 0;

    u64 mismatches_r_above_2k = 0;

    u64 max_threshold = 0;
    u64 max_threshold_k_sum = 0;

    u64 max_ratio_num = 0;
    u64 max_ratio_den = 1;

    u64 max_threshold_r = 0;
    u64 max_threshold_r_m = 0;
    u64 max_threshold_r_k1 = 0;
    u64 max_threshold_r_k2 = 0;
    u64 max_threshold_r_m1 = 0;
    u64 max_threshold_r_m2 = 0;
    int max_threshold_r_sign1 = 0;
    int max_threshold_r_sign2 = 0;
};

static mpz_class mpz_from_u64(u64 value)
{
    return mpz_class(std::to_string(value));
}

static bool is_prime(u64 n)
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

static bool valid_witness(
    u64 r,
    u64 m,
    u64 k,
    u64& t,
    int& sign)
{
    const u64 value_plus  = m * r + 1;
    const u64 value_minus = (m * r >= 1) ? (m * r - 1) : 0;

    bool found = false;
    u64 best_t = std::numeric_limits<u64>::max();
    int best_sign = 0;

    if (value_plus % k == 0)
    {
        const u64 candidate_t = value_plus / k;

        if (candidate_t < best_t)
        {
            best_t = candidate_t;
            best_sign = +1;
            found = true;
        }
    }

    if (m * r >= 1 && value_minus % k == 0)
    {
        const u64 candidate_t = value_minus / k;

        if (candidate_t < best_t)
        {
            best_t = candidate_t;
            best_sign = -1;
            found = true;
        }
    }

    if (!found)
        return false;

    t = best_t;
    sign = best_sign;
    return true;
}

static Witness best_for_m(u64 r, u64 m, u64 K_LIMIT_LOCAL)
{
    Witness best;

    for (u64 k = 1; k <= K_LIMIT_LOCAL; ++k)
    {
        u64 t = 0;
        int sign = 0;

        if (!valid_witness(r, m, k, t, sign))
            continue;

        if (!best.valid ||
            k > best.k ||
            (k == best.k && t < best.t))
        {
            best.valid = true;
            best.m = m;
            best.k = k;
            best.t = t;
            best.sign = sign;
        }
    }

    return best;
}

static bool normalized_first(const Witness& a, const Witness& b)
{
    // Compare k1/m1 > k2/m2 without floating point.
    return (__uint128_t)a.k * b.m >
           (__uint128_t)b.k * a.m;
}

static bool normalized_second(const Witness& a, const Witness& b)
{
    return (__uint128_t)a.k * b.m <
           (__uint128_t)b.k * a.m;
}

static bool exact_first(u64 r, const Witness& a, const Witness& b)
{
    // Exact t:
    //   t = (m*r + sign) / k
    //
    // Compare:
    //   (m1*r + e1)/k1 < (m2*r + e2)/k2
    //
    // using 128-bit arithmetic.

    const __uint128_t left =
        (__uint128_t)(a.m * r + static_cast<u64>(a.sign)) * b.k;

    const __uint128_t right =
        (__uint128_t)(b.m * r + static_cast<u64>(b.sign)) * a.k;

    return left < right;
}

static bool exact_second(u64 r, const Witness& a, const Witness& b)
{
    const __uint128_t left =
        (__uint128_t)(a.m * r + static_cast<u64>(a.sign)) * b.k;

    const __uint128_t right =
        (__uint128_t)(b.m * r + static_cast<u64>(b.sign)) * a.k;

    return left > right;
}

static u64 exact_threshold_r(
    const Witness& a,
    const Witness& b)
{
    // Only called when normalized A comes first:
    //
    // Delta = m2*k1 - m1*k2 > 0
    //
    // A < B iff
    //   Delta*r > e1*k2 - e2*k1
    //
    // Therefore:
    //   r >= floor(C/Delta) + 1
    //
    const __int128 delta =
        static_cast<__int128>(b.m) * a.k -
        static_cast<__int128>(a.m) * b.k;

    const __int128 correction =
        static_cast<__int128>(a.sign) * b.k -
        static_cast<__int128>(b.sign) * a.k;

    if (delta <= 0)
        return 0;

    if (correction < 0)
        return 1;

    const __int128 threshold =
        correction / delta + 1;

    if (threshold <= 0)
        return 1;

    if (threshold > std::numeric_limits<u64>::max())
        return std::numeric_limits<u64>::max();

    return static_cast<u64>(threshold);
}

static bool correction_nonnegative(
    const Witness& a,
    const Witness& b)
{
    const __int128 correction =
        static_cast<__int128>(a.sign) * b.k -
        static_cast<__int128>(b.sign) * a.k;

    return correction >= 0;
}

static void analyze_pair(
    u64 r,
    const Witness& a,
    const Witness& b,
    Stats& stats)
{
    if (!a.valid || !b.valid)
        return;

    ++stats.pair_tests;

    const bool norm_a_first = normalized_first(a, b);
    const bool norm_b_first = normalized_second(a, b);

    if (!norm_a_first && !norm_b_first)
        return; // normalized tie

    if (norm_a_first)
        ++stats.normalized_first;
    else
        ++stats.normalized_second;

    bool exact_a_first = exact_first(r, a, b);
    bool exact_b_first = exact_second(r, a, b);

    bool agree = false;

    if (norm_a_first)
        agree = exact_a_first;
    else
        agree = exact_b_first;

    if (agree)
    {
        ++stats.exact_agree;
    }
    else
    {
        ++stats.exact_disagree;

        if (r > a.k + b.k)
            ++stats.mismatches_r_above_2k;
    }

    if (!norm_a_first)
        return;

    const __int128 delta =
        static_cast<__int128>(b.m) * a.k -
        static_cast<__int128>(a.m) * b.k;

    const __int128 correction =
        static_cast<__int128>(a.sign) * b.k -
        static_cast<__int128>(b.sign) * a.k;

    if (delta <= 0)
        return;

    const u64 threshold = exact_threshold_r(a, b);

    ++stats.threshold_tests;

    if (threshold > a.k + b.k)
        ++stats.threshold_bound_failures;

    if (threshold > stats.max_threshold)
    {
        stats.max_threshold = threshold;
        stats.max_threshold_k_sum = a.k + b.k;
    }

    // Check universal bound r > k1+k2.
    if (r > a.k + b.k)
    {
        ++stats.global_bound_tests;

        if (!exact_a_first)
            ++stats.global_bound_failures;
    }

    if (correction > 0)
    {
        const u64 correction_u64 =
            static_cast<u64>(correction);

        const u64 delta_u64 =
            static_cast<u64>(delta);

        const u64 ratio_num =
            correction_u64;

        const u64 ratio_den =
            delta_u64 * (a.k + b.k);

        if (static_cast<__int128>(ratio_num) * stats.max_ratio_den >
            static_cast<__int128>(stats.max_ratio_num) * ratio_den)
        {
            stats.max_ratio_num = ratio_num;
            stats.max_ratio_den = ratio_den;
        }
    }

    if (threshold >= stats.max_threshold_r)
    {
        stats.max_threshold_r = threshold;
        stats.max_threshold_r_k1 = a.k;
        stats.max_threshold_r_k2 = b.k;
        stats.max_threshold_r_m1 = a.m;
        stats.max_threshold_r_m2 = b.m;
        stats.max_threshold_r_sign1 = a.sign;
        stats.max_threshold_r_sign2 = b.sign;
    }

    // Directly verify the threshold formula at threshold-1 and threshold,
    // provided threshold is meaningful and small enough.
    if (threshold > 1 && threshold <= PRIME_LIMIT + 1)
    {
        const u64 below = threshold - 1;

        const bool below_exact = exact_first(below, a, b);
        const bool at_exact    = exact_first(threshold, a, b);

        const bool expected_below =
            (static_cast<__int128>(delta) * below) > correction;

        const bool expected_at =
            (static_cast<__int128>(delta) * threshold) > correction;

        if (below_exact != expected_below ||
            at_exact != expected_at)
        {
            // This should be impossible if the derivation is correct.
            ++stats.threshold_bound_failures;
        }
    }

    (void)correction_nonnegative;
}

int main()
{
    std::cout << "START EXPERIMENT 409\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    for (u64 r : primes)
    {
        std::vector<Witness> witnesses(M_MAX + 1);

        for (u64 m = 1; m <= M_MAX; ++m)
        {
            witnesses[m] =
                best_for_m(r, m, K_LIMIT);
        }

        for (u64 K = 1; K <= K_LIMIT; ++K)
        {
            std::vector<Witness> current(M_MAX + 1);

            for (u64 m = 1; m <= M_MAX; ++m)
            {
                Witness best;

                for (u64 k = 1; k <= K; ++k)
                {
                    u64 t = 0;
                    int sign = 0;

                    if (!valid_witness(r, m, k, t, sign))
                        continue;

                    if (!best.valid ||
                        k > best.k ||
                        (k == best.k && t < best.t))
                    {
                        best.valid = true;
                        best.m = m;
                        best.k = k;
                        best.t = t;
                        best.sign = sign;
                    }
                }

                current[m] = best;
            }

            for (u64 m1 = 1; m1 <= M_MAX; ++m1)
            {
                if (!current[m1].valid)
                    continue;

                for (u64 m2 = m1 + 1; m2 <= M_MAX; ++m2)
                {
                    if (!current[m2].valid)
                        continue;

                    analyze_pair(
                        r,
                        current[m1],
                        current[m2],
                        stats
                    );
                }
            }
        }
    }

    std::cout
        << "PRIME_LIMIT=" << PRIME_LIMIT << '\n'
        << "K_LIMIT=" << K_LIMIT << '\n'
        << "M_MAX=" << M_MAX << '\n'
        << "PRIME_COUNT=" << primes.size() << '\n'
        << '\n';

    std::cout
        << "PAIR_TESTS=" << stats.pair_tests << '\n'
        << "NORMALIZED_FIRST=" << stats.normalized_first << '\n'
        << "NORMALIZED_SECOND=" << stats.normalized_second << '\n'
        << "EXACT_AGREE=" << stats.exact_agree << '\n'
        << "EXACT_DISAGREE=" << stats.exact_disagree << '\n'
        << '\n';

    std::cout
        << "THRESHOLD_TESTS=" << stats.threshold_tests << '\n'
        << "THRESHOLD_BOUND_FAILURES="
        << stats.threshold_bound_failures << '\n'
        << '\n';

    std::cout
        << "GLOBAL_BOUND_TESTS="
        << stats.global_bound_tests << '\n'
        << "GLOBAL_BOUND_FAILURES="
        << stats.global_bound_failures << '\n'
        << '\n';

    std::cout
        << "MISMATCHES_R_ABOVE_2K="
        << stats.mismatches_r_above_2k << '\n'
        << '\n';

    std::cout
        << "MAX_THRESHOLD="
        << stats.max_threshold << '\n'
        << "MAX_THRESHOLD_K_SUM="
        << stats.max_threshold_k_sum << '\n'
        << '\n';

    std::cout
        << "MAX_THRESHOLD_R="
        << stats.max_threshold_r << '\n'
        << "MAX_THRESHOLD_R_K1="
        << stats.max_threshold_r_k1 << '\n'
        << "MAX_THRESHOLD_R_K2="
        << stats.max_threshold_r_k2 << '\n'
        << "MAX_THRESHOLD_R_M1="
        << stats.max_threshold_r_m1 << '\n'
        << "MAX_THRESHOLD_R_M2="
        << stats.max_threshold_r_m2 << '\n'
        << "MAX_THRESHOLD_R_SIGN1="
        << stats.max_threshold_r_sign1 << '\n'
        << "MAX_THRESHOLD_R_SIGN2="
        << stats.max_threshold_r_sign2 << '\n'
        << '\n';

    if (stats.max_ratio_den != 0)
    {
        std::cout
            << "MAX_CORRECTION_OVER_DELTA_KSUM="
            << stats.max_ratio_num
            << "/"
            << stats.max_ratio_den
            << '\n';
    }

    std::cout << '\n';

    if (stats.threshold_bound_failures == 0)
    {
        std::cout
            << "THRESHOLD_BOUND_STATUS=PASS\n";
    }
    else
    {
        std::cout
            << "THRESHOLD_BOUND_STATUS=FAIL\n";
    }

    if (stats.global_bound_failures == 0)
    {
        std::cout
            << "GLOBAL_BOUND_STATUS=PASS\n";
    }
    else
    {
        std::cout
            << "GLOBAL_BOUND_STATUS=FAIL\n";
    }

    if (stats.mismatches_r_above_2k == 0)
    {
        std::cout
            << "R_GT_K1_PLUS_K2_STATUS=PASS\n";
    }
    else
    {
        std::cout
            << "R_GT_K1_PLUS_K2_STATUS=FAIL\n";
    }

    std::cout << "FINISHED EXPERIMENT 409\n";

    return 0;
}