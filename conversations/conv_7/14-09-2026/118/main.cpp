#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using u64 = std::uint64_t;
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
    u64 pair_tests = 0;

    u64 normalized_first = 0;
    u64 normalized_second = 0;

    u64 threshold_tests = 0;

    u64 threshold_bound_plus_one_failures = 0;

    u64 equality_cases = 0;
    u64 equality_delta_not_one = 0;
    u64 equality_c_not_ksum = 0;
    u64 equality_characterization_failures = 0;

    u64 equality_sign_pp = 0;
    u64 equality_sign_pm = 0;
    u64 equality_sign_mp = 0;
    u64 equality_sign_mm = 0;

    u64 max_threshold = 0;
    u64 max_k_sum = 0;
    u64 max_threshold_excess = 0;

    u64 min_delta = std::numeric_limits<u64>::max();
    u64 max_delta = 0;

    u64 min_correction = std::numeric_limits<u64>::max();
    u64 max_correction = 0;

    bool have_first_equality = false;

    Witness first_a;
    Witness first_b;

    u64 first_r = 0;
    u64 first_threshold = 0;
    u64 first_k_sum = 0;
    u64 first_delta = 0;
    i128 first_correction = 0;
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

static bool valid_witness(
    u64 r,
    u64 m,
    u64 k,
    u64& t,
    int& sign)
{
    const u64 mr = m * r;

    bool found = false;
    u64 best_t = std::numeric_limits<u64>::max();
    int best_sign = 0;

    {
        const u64 value = mr + 1;

        if (value % k == 0)
        {
            const u64 candidate_t = value / k;

            if (candidate_t < best_t)
            {
                best_t = candidate_t;
                best_sign = +1;
                found = true;
            }
        }
    }

    if (mr >= 1)
    {
        const u64 value = mr - 1;

        if (value % k == 0)
        {
            const u64 candidate_t = value / k;

            if (candidate_t < best_t)
            {
                best_t = candidate_t;
                best_sign = -1;
                found = true;
            }
        }
    }

    if (!found)
        return false;

    t = best_t;
    sign = best_sign;

    return true;
}

static Witness best_for_m(
    u64 r,
    u64 m,
    u64 k_limit)
{
    Witness best;

    for (u64 k = 1; k <= k_limit; ++k)
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

static bool normalized_first(
    const Witness& a,
    const Witness& b)
{
    return (__uint128_t)a.k * b.m >
           (__uint128_t)b.k * a.m;
}

static u64 threshold_from_pair(
    const Witness& a,
    const Witness& b,
    u64& delta_out,
    i128& correction_out)
{
    const i128 delta =
        static_cast<i128>(b.m) * a.k -
        static_cast<i128>(a.m) * b.k;

    const i128 correction =
        static_cast<i128>(a.sign) * b.k -
        static_cast<i128>(b.sign) * a.k;

    delta_out = static_cast<u64>(delta);
    correction_out = correction;

    if (delta <= 0)
        return 0;

    if (correction < 0)
        return 1;

    const i128 threshold =
        correction / delta + 1;

    return static_cast<u64>(threshold);
}

static void record_equality_sign(
    const Witness& a,
    const Witness& b,
    Stats& stats)
{
    if (a.sign == +1 && b.sign == +1)
        ++stats.equality_sign_pp;
    else if (a.sign == +1 && b.sign == -1)
        ++stats.equality_sign_pm;
    else if (a.sign == -1 && b.sign == +1)
        ++stats.equality_sign_mp;
    else
        ++stats.equality_sign_mm;
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

    if (!normalized_first(a, b))
    {
        ++stats.normalized_second;
        return;
    }

    ++stats.normalized_first;

    u64 delta = 0;
    i128 correction = 0;

    const u64 threshold =
        threshold_from_pair(
            a,
            b,
            delta,
            correction
        );

    ++stats.threshold_tests;

    const u64 k_sum = a.k + b.k;

    if (threshold > k_sum + 1)
        ++stats.threshold_bound_plus_one_failures;

    const u64 excess =
        threshold > k_sum
            ? threshold - k_sum
            : 0;

    if (excess > stats.max_threshold_excess)
        stats.max_threshold_excess = excess;

    if (threshold > stats.max_threshold)
    {
        stats.max_threshold = threshold;
        stats.max_k_sum = k_sum;
    }

    if (delta < stats.min_delta)
        stats.min_delta = delta;

    if (delta > stats.max_delta)
        stats.max_delta = delta;

    if (correction >= 0)
    {
        const u64 c =
            static_cast<u64>(correction);

        if (c < stats.min_correction)
            stats.min_correction = c;

        if (c > stats.max_correction)
            stats.max_correction = c;
    }

    /*
        Sharp boundary:

            threshold = k1 + k2 + 1

        should imply:

            Delta = 1
            C     = k1 + k2
    */
    if (threshold == k_sum + 1)
    {
        ++stats.equality_cases;

        record_equality_sign(a, b, stats);

        const bool delta_ok =
            (delta == 1);

        const bool correction_ok =
            (correction == static_cast<i128>(k_sum));

        if (!delta_ok)
            ++stats.equality_delta_not_one;

        if (!correction_ok)
            ++stats.equality_c_not_ksum;

        if (!delta_ok || !correction_ok)
            ++stats.equality_characterization_failures;

        if (!stats.have_first_equality)
        {
            stats.have_first_equality = true;

            stats.first_a = a;
            stats.first_b = b;

            stats.first_r = r;
            stats.first_threshold = threshold;
            stats.first_k_sum = k_sum;
            stats.first_delta = delta;
            stats.first_correction = correction;
        }
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

int main()
{
    std::cout << "START EXPERIMENT 410\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    for (u64 r : primes)
    {
        for (u64 K = 1; K <= K_LIMIT; ++K)
        {
            std::vector<Witness> current(M_MAX + 1);

            for (u64 m = 1; m <= M_MAX; ++m)
            {
                current[m] =
                    best_for_m(
                        r,
                        m,
                        K
                    );
            }

            for (u64 m1 = 1; m1 <= M_MAX; ++m1)
            {
                if (!current[m1].valid)
                    continue;

                for (u64 m2 = m1 + 1;
                     m2 <= M_MAX;
                     ++m2)
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
        << "PAIR_TESTS="
        << stats.pair_tests << '\n'
        << "NORMALIZED_FIRST="
        << stats.normalized_first << '\n'
        << "NORMALIZED_SECOND="
        << stats.normalized_second << '\n'
        << '\n';

    std::cout
        << "THRESHOLD_TESTS="
        << stats.threshold_tests << '\n'
        << "THRESHOLD_BOUND_PLUS_ONE_FAILURES="
        << stats.threshold_bound_plus_one_failures << '\n'
        << '\n';

    std::cout
        << "EQUALITY_CASES="
        << stats.equality_cases << '\n'
        << "EQUALITY_DELTA_NOT_ONE="
        << stats.equality_delta_not_one << '\n'
        << "EQUALITY_C_NOT_KSUM="
        << stats.equality_c_not_ksum << '\n'
        << "EQUALITY_CHARACTERIZATION_FAILURES="
        << stats.equality_characterization_failures << '\n'
        << '\n';

    std::cout
        << "EQUALITY_SIGN_PP="
        << stats.equality_sign_pp << '\n'
        << "EQUALITY_SIGN_PM="
        << stats.equality_sign_pm << '\n'
        << "EQUALITY_SIGN_MP="
        << stats.equality_sign_mp << '\n'
        << "EQUALITY_SIGN_MM="
        << stats.equality_sign_mm << '\n'
        << '\n';

    std::cout
        << "MAX_THRESHOLD="
        << stats.max_threshold << '\n'
        << "MAX_K_SUM="
        << stats.max_k_sum << '\n'
        << "MAX_THRESHOLD_EXCESS="
        << stats.max_threshold_excess << '\n'
        << '\n';

    std::cout
        << "MIN_DELTA="
        << stats.min_delta << '\n'
        << "MAX_DELTA="
        << stats.max_delta << '\n'
        << "MIN_CORRECTION="
        << stats.min_correction << '\n'
        << "MAX_CORRECTION="
        << stats.max_correction << '\n'
        << '\n';

    if (stats.have_first_equality)
    {
        std::cout << "FIRST_EQUALITY_CASE\n";

        std::cout
            << "R="
            << stats.first_r << '\n'
            << "THRESHOLD="
            << stats.first_threshold << '\n'
            << "K_SUM="
            << stats.first_k_sum << '\n'
            << "DELTA="
            << stats.first_delta << '\n'
            << "CORRECTION="
            << static_cast<long long>(
                stats.first_correction)
            << '\n';

        print_witness("A", stats.first_a);
        print_witness("B", stats.first_b);

        std::cout << '\n';
    }

    std::cout
        << "THRESHOLD_BOUND_PLUS_ONE_STATUS="
        << (
            stats.threshold_bound_plus_one_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "EQUALITY_CHARACTERIZATION_STATUS="
        << (
            stats.equality_characterization_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout << "FINISHED EXPERIMENT 410\n";

    return 0;
}