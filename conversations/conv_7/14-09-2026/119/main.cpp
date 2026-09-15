#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
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
    u64 determinant_one_pairs = 0;

    u64 prime_tests = 0;

    u64 congruence_tests = 0;
    u64 congruence_failures = 0;

    u64 reverse_tests = 0;
    u64 reverse_failures = 0;

    u64 witness_tests = 0;
    u64 witness_failures = 0;

    u64 sharp_threshold_tests = 0;
    u64 sharp_threshold_failures = 0;

    u64 sharp_observed = 0;
    u64 sharp_congruence_failures = 0;

    u64 sharp_sign_pm_failures = 0;

    u64 t_difference_tests = 0;
    u64 t_difference_failures = 0;

    u64 min_lcm = std::numeric_limits<u64>::max();
    u64 max_lcm = 0;

    u64 min_prime_r = std::numeric_limits<u64>::max();
    u64 max_prime_r = 0;

    bool have_first_pair = false;
    u64 first_m1 = 0;
    u64 first_m2 = 0;
    u64 first_k1 = 0;
    u64 first_k2 = 0;
    u64 first_r = 0;
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

static Witness make_witness(
    u64 r,
    u64 m,
    u64 k,
    int sign)
{
    Witness w;

    const i128 value =
        static_cast<i128>(m) * r + sign;

    if (value <= 0)
        return w;

    if (value % k != 0)
        return w;

    const i128 t128 = value / k;

    if (t128 > std::numeric_limits<u64>::max())
        return w;

    w.valid = true;
    w.m = m;
    w.k = k;
    w.t = static_cast<u64>(t128);
    w.sign = sign;

    return w;
}

static bool normalized_first(
    const Witness& a,
    const Witness& b)
{
    return (__uint128_t)a.k * b.m >
           (__uint128_t)b.k * a.m;
}

static bool determinant_one(
    u64 m1,
    u64 m2,
    u64 k1,
    u64 k2)
{
    return static_cast<i128>(m2) * k1 -
           static_cast<i128>(m1) * k2 == 1;
}

static u64 lcm_u64(u64 a, u64 b)
{
    return (a / std::gcd(a, b)) * b;
}

static bool congruence_condition(
    u64 r,
    u64 k1,
    u64 k2)
{
    const u64 L = lcm_u64(k1, k2);
    const u64 residue = (k1 + k2) % L;

    return r % L == residue;
}

static bool congruence_components(
    u64 r,
    u64 k1,
    u64 k2)
{
    return
        (r % k1 == k2 % k1) &&
        (r % k2 == k1 % k2);
}

static bool verify_t_difference(
    u64 r,
    const Witness& a,
    const Witness& b)
{
    if (!a.valid || !b.valid)
        return false;

    const i128 lhs =
        static_cast<i128>(b.t) -
        static_cast<i128>(a.t);

    const i128 delta =
        static_cast<i128>(b.m) * a.k -
        static_cast<i128>(a.m) * b.k;

    const i128 numerator =
        delta * r -
        static_cast<i128>(a.k + b.k);

    const i128 denominator =
        static_cast<i128>(a.k) * b.k;

    if (denominator == 0)
        return false;

    return lhs * denominator == numerator;
}

static void test_determinant_pair(
    u64 m1,
    u64 m2,
    u64 k1,
    u64 k2,
    const std::vector<u64>& primes,
    Stats& stats)
{
    if (!determinant_one(m1, m2, k1, k2))
        return;

    ++stats.determinant_one_pairs;

    const u64 L = lcm_u64(k1, k2);

    if (L < stats.min_lcm)
        stats.min_lcm = L;

    if (L > stats.max_lcm)
        stats.max_lcm = L;

    for (u64 r : primes)
    {
        ++stats.prime_tests;

        const bool congruence =
            congruence_condition(r, k1, k2);

        const bool components =
            congruence_components(r, k1, k2);

        ++stats.congruence_tests;

        if (!congruence || !components)
            ++stats.congruence_failures;

        /*
            Reverse direction:

            determinant = 1
            + congruence

            should imply

                k1 | (m1*r + 1)
                k2 | (m2*r - 1)
        */
        if (!congruence)
            continue;

        ++stats.reverse_tests;

        Witness a =
            make_witness(r, m1, k1, +1);

        Witness b =
            make_witness(r, m2, k2, -1);

        if (!a.valid || !b.valid)
        {
            ++stats.reverse_failures;
            continue;
        }

        ++stats.witness_tests;

        if (!normalized_first(a, b))
        {
            ++stats.witness_failures;
            continue;
        }

        /*
            Since Delta = 1 and C = k1+k2,
            the theoretical sharp threshold is

                r_min = k1+k2+1.
        */

        const u64 threshold =
            k1 + k2 + 1;

        if (r == threshold)
        {
            ++stats.sharp_threshold_tests;

            if (a.t >= b.t)
                ++stats.sharp_threshold_failures;
        }

        /*
            At any valid r, verify

                t2 - t1
                  = [r-(k1+k2)]/(k1*k2)

            because Delta = 1.
        */

        ++stats.t_difference_tests;

        if (!verify_t_difference(r, a, b))
            ++stats.t_difference_failures;

        /*
            Record the first valid prime witness.
        */

        if (!stats.have_first_pair)
        {
            stats.have_first_pair = true;

            stats.first_m1 = m1;
            stats.first_m2 = m2;
            stats.first_k1 = k1;
            stats.first_k2 = k2;
            stats.first_r = r;
        }

        if (r < stats.min_prime_r)
            stats.min_prime_r = r;

        if (r > stats.max_prime_r)
            stats.max_prime_r = r;
    }
}

static void analyze_observed_sharp_case(
    u64 r,
    const Witness& a,
    const Witness& b,
    Stats& stats)
{
    if (!a.valid || !b.valid)
        return;

    const i128 delta =
        static_cast<i128>(b.m) * a.k -
        static_cast<i128>(a.m) * b.k;

    const i128 correction =
        static_cast<i128>(a.sign) * b.k -
        static_cast<i128>(b.sign) * a.k;

    if (delta != 1 ||
        correction != static_cast<i128>(a.k + b.k))
    {
        return;
    }

    ++stats.sharp_observed;

    if (a.sign != +1 ||
        b.sign != -1)
    {
        ++stats.sharp_sign_pm_failures;
        return;
    }

    const u64 L =
        lcm_u64(a.k, b.k);

    if (!congruence_condition(r, a.k, b.k))
        ++stats.sharp_congruence_failures;

    if (r % a.k != b.k % a.k)
        ++stats.sharp_congruence_failures;

    if (r % b.k != a.k % b.k)
        ++stats.sharp_congruence_failures;

    const u64 threshold =
        a.k + b.k + 1;

    if (threshold <= r)
        ++stats.sharp_threshold_failures;

    (void)L;
}

int main()
{
    std::cout << "START EXPERIMENT 411\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    Stats stats;

    /*
        Enumerate all determinant-1 pairs

            m2*k1 - m1*k2 = 1

        with

            m1,m2 <= 64
            k1,k2 <= 500.

        For fixed m1,m2,k1, k2 is determined.
    */

    for (u64 m1 = 1; m1 <= M_MAX; ++m1)
    {
        for (u64 m2 = m1 + 1;
             m2 <= M_MAX;
             ++m2)
        {
            for (u64 k1 = 1; k1 <= K_LIMIT; ++k1)
            {
                const i128 numerator =
                    static_cast<i128>(m2) * k1 - 1;

                if (numerator <= 0)
                    continue;

                if (numerator % m1 != 0)
                    continue;

                const i128 k2_128 =
                    numerator / m1;

                if (k2_128 < 1 ||
                    k2_128 > K_LIMIT)
                {
                    continue;
                }

                const u64 k2 =
                    static_cast<u64>(k2_128);

                if (!determinant_one(
                        m1,
                        m2,
                        k1,
                        k2))
                {
                    continue;
                }

                test_determinant_pair(
                    m1,
                    m2,
                    k1,
                    k2,
                    primes,
                    stats
                );
            }
        }
    }

    /*
        Independently reproduce the observed sharp-bound
        witnesses using the same best-witness construction
        from Experiments 409/410.

        This checks that every observed sharp case also
        satisfies the congruence.
    */

    u64 observed_sharp_cases = 0;

    for (u64 r : primes)
    {
        for (u64 K = 1; K <= K_LIMIT; ++K)
        {
            std::vector<Witness> current(
                M_MAX + 1
            );

            for (u64 m = 1; m <= M_MAX; ++m)
            {
                Witness best;

                for (u64 k = 1; k <= K; ++k)
                {
                    Witness plus =
                        make_witness(
                            r,
                            m,
                            k,
                            +1
                        );

                    Witness minus =
                        make_witness(
                            r,
                            m,
                            k,
                            -1
                        );

                    if (plus.valid &&
                        (!best.valid ||
                         k > best.k ||
                         (k == best.k &&
                          plus.t < best.t)))
                    {
                        best = plus;
                    }

                    if (minus.valid &&
                        (!best.valid ||
                         k > best.k ||
                         (k == best.k &&
                          minus.t < best.t)))
                    {
                        best = minus;
                    }
                }

                current[m] = best;
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

                    const Witness& a = current[m1];
                    const Witness& b = current[m2];

                    if (!normalized_first(a, b))
                        continue;

                    const i128 delta =
                        static_cast<i128>(b.m) * a.k -
                        static_cast<i128>(a.m) * b.k;

                    const i128 correction =
                        static_cast<i128>(a.sign) * b.k -
                        static_cast<i128>(b.sign) * a.k;

                    if (delta != 1)
                        continue;

                    if (correction !=
                        static_cast<i128>(a.k + b.k))
                    {
                        continue;
                    }

                    ++observed_sharp_cases;

                    analyze_observed_sharp_case(
                        r,
                        a,
                        b,
                        stats
                    );
                }
            }
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
        << primes.size() << '\n'
        << '\n';

    std::cout
        << "DETERMINANT_ONE_PAIRS="
        << stats.determinant_one_pairs << '\n'
        << "PRIME_TESTS="
        << stats.prime_tests << '\n'
        << '\n';

    std::cout
        << "CONGRUENCE_TESTS="
        << stats.congruence_tests << '\n'
        << "CONGRUENCE_FAILURES="
        << stats.congruence_failures << '\n'
        << '\n';

    std::cout
        << "REVERSE_TESTS="
        << stats.reverse_tests << '\n'
        << "REVERSE_FAILURES="
        << stats.reverse_failures << '\n'
        << '\n';

    std::cout
        << "WITNESS_TESTS="
        << stats.witness_tests << '\n'
        << "WITNESS_FAILURES="
        << stats.witness_failures << '\n'
        << '\n';

    std::cout
        << "SHARP_THRESHOLD_TESTS="
        << stats.sharp_threshold_tests << '\n'
        << "SHARP_THRESHOLD_FAILURES="
        << stats.sharp_threshold_failures << '\n'
        << '\n';

    std::cout
        << "T_DIFFERENCE_TESTS="
        << stats.t_difference_tests << '\n'
        << "T_DIFFERENCE_FAILURES="
        << stats.t_difference_failures << '\n'
        << '\n';

    std::cout
        << "OBSERVED_SHARP_CASES="
        << observed_sharp_cases << '\n'
        << "OBSERVED_SHARP_CONGRUENCE_FAILURES="
        << stats.sharp_congruence_failures << '\n'
        << "OBSERVED_SHARP_SIGN_PM_FAILURES="
        << stats.sharp_sign_pm_failures << '\n'
        << '\n';

    std::cout
        << "MIN_LCM="
        << stats.min_lcm << '\n'
        << "MAX_LCM="
        << stats.max_lcm << '\n'
        << '\n';

    std::cout
        << "MIN_VALID_PRIME_R="
        << stats.min_prime_r << '\n'
        << "MAX_VALID_PRIME_R="
        << stats.max_prime_r << '\n'
        << '\n';

    if (stats.have_first_pair)
    {
        std::cout
            << "FIRST_VALID_DETERMINANT_ONE_PAIR\n"
            << "M1=" << stats.first_m1 << '\n'
            << "M2=" << stats.first_m2 << '\n'
            << "K1=" << stats.first_k1 << '\n'
            << "K2=" << stats.first_k2 << '\n'
            << "R=" << stats.first_r << '\n'
            << '\n';
    }

    std::cout
        << "CONGRUENCE_STATUS="
        << (
            stats.congruence_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "REVERSE_CONSTRUCTION_STATUS="
        << (
            stats.reverse_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "WITNESS_CONSTRUCTION_STATUS="
        << (
            stats.witness_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "T_DIFFERENCE_STATUS="
        << (
            stats.t_difference_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "OBSERVED_SHARP_CONGRUENCE_STATUS="
        << (
            stats.sharp_congruence_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 411\n";

    return 0;
}
