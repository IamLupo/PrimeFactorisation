#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr int EXPERIMENT = 154;

// ============================================================
// u128 -> string
// ============================================================

std::string u128_to_string(u128 x)
{
    if (x == 0)
        return "0";

    std::string s;

    while (x > 0)
    {
        int digit = static_cast<int>(x % 10);
        s.push_back(static_cast<char>('0' + digit));
        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    return s;
}

// ============================================================
// gcd
// ============================================================

u128 gcd_u128(u128 a, u128 b)
{
    while (b != 0)
    {
        u128 r = a % b;
        a = b;
        b = r;
    }

    return a;
}

// ============================================================
// modular multiplication
// ============================================================

u128 mul_mod(u128 a, u128 b, u128 mod)
{
    return (a * b) % mod;
}

// ============================================================
// integer square root
// ============================================================

u64 isqrt_u64(u64 n)
{
    u64 x = static_cast<u64>(
        std::sqrt(static_cast<long double>(n))
    );

    while ((x + 1) <= n / (x + 1))
        ++x;

    while (x > n / x)
        --x;

    return x;
}

// ============================================================
// primality
// ============================================================

bool is_prime_u64(u64 n)
{
    if (n < 2)
        return false;

    const u64 small_primes[] =
    {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };

    for (u64 p : small_primes)
    {
        if (n == p)
            return true;

        if (n % p == 0)
            return false;
    }

    u64 d = n - 1;
    int s = 0;

    while ((d & 1ULL) == 0)
    {
        d >>= 1;
        ++s;
    }

    auto mod_mul_u64 =
        [](u64 a, u64 b, u64 mod) -> u64
        {
            return static_cast<u64>(
                (static_cast<u128>(a) * b) % mod
            );
        };

    auto mod_pow =
        [&](u64 a, u64 e, u64 mod) -> u64
        {
            u64 result = 1;

            while (e)
            {
                if (e & 1ULL)
                    result = mod_mul_u64(result, a, mod);

                a = mod_mul_u64(a, a, mod);
                e >>= 1;
            }

            return result;
        };

    const u64 bases[] =
    {
        2, 3, 5, 7, 11, 13, 17
    };

    for (u64 a : bases)
    {
        if (a >= n)
            continue;

        u64 x = mod_pow(a, d, n);

        if (x == 1 || x == n - 1)
            continue;

        bool passed = false;

        for (int r = 1; r < s; ++r)
        {
            x = mod_mul_u64(x, x, n);

            if (x == n - 1)
            {
                passed = true;
                break;
            }
        }

        if (!passed)
            return false;
    }

    return true;
}

// ============================================================
// random prime near target
// ============================================================

u64 random_prime_near(
    u64 center,
    u64 radius,
    std::mt19937_64& rng
)
{
    std::uniform_int_distribution<u64> dist(
        0,
        radius
    );

    for (;;)
    {
        u64 offset = dist(rng);

        u64 candidate;

        if (rng() & 1ULL)
        {
            candidate = center + offset;
        }
        else
        {
            if (offset >= center)
                continue;

            candidate = center - offset;
        }

        candidate |= 1ULL;

        if (candidate < 3)
            continue;

        if (is_prime_u64(candidate))
            return candidate;
    }
}

// ============================================================
// generate semiprime
// ============================================================

std::pair<u64, u64> make_semiprime(
    u64 target,
    double ratio,
    std::mt19937_64& rng
)
{
    long double p_real =
        std::sqrt(
            static_cast<long double>(target) /
            ratio
        );

    u64 p_center =
        static_cast<u64>(p_real);

    u64 q_center =
        static_cast<u64>(
            static_cast<long double>(p_center) *
            ratio
        );

    const u64 radius = 100000;

    for (;;)
    {
        u64 p =
            random_prime_near(
                p_center,
                radius,
                rng
            );

        u64 q =
            random_prime_near(
                q_center,
                radius,
                rng
            );

        if (p == q)
            continue;

        if (p > q)
            std::swap(p, q);

        u128 N =
            static_cast<u128>(p) *
            static_cast<u128>(q);

        u128 target128 =
            static_cast<u128>(target);

        if (N >= target128 * 9 / 10 &&
            N <= target128 * 11 / 10)
        {
            return {p, q};
        }
    }
}

// ============================================================
// Probe result
// ============================================================

struct ProbeResult
{
    bool success;
    u64 candidate;
    char reason;
};

// ============================================================
// Binary search inside ONE safe segment.
//
// IMPORTANT:
//
// m = U is fixed for the entire search.
//
// Product:
//
//     P(t) = product(t+1 ... U)
//
// Because every MISS gives a smaller t, we can extend the
// existing product instead of recomputing it.
//
// Exactly one gcd is performed per probe.
// ============================================================

ProbeResult search_segment(
    u128 N,
    u64 L,
    u64 U,
    u64& product_accumulator_steps,
    u64& gcd_calls,
    u64& probes
)
{
    if (L > U)
        return {false, 0, 'X'};

    // Safe condition.
    if (U >= 2 * L)
        return {false, 0, 'U'};

    u64 left = L;
    u64 right = U;

    // This is the fixed m.
    const u64 m = U;

    // We have not yet accumulated any terms.
    //
    // After the first probe at t:
    //
    //     accumulator =
    //         product(t+1 ... U)
    //
    u64 accumulated_from = U + 1;

    u128 accumulator = 1;

    while (left < right)
    {
        ++probes;

        u64 t =
            left +
            (right - left) / 2;

        // ----------------------------------------------------
        // We need:
        //
        //     product(t+1 ... U)
        //
        // If this is the first probe, multiply t+1..U.
        //
        // Later, accumulated_from = previous_t + 1.
        // Since t decreases, only multiply:
        //
        //     t+1 ... previous_t
        // ----------------------------------------------------

        u64 start;

        if (accumulated_from > U)
        {
            start = t + 1;
        }
        else
        {
            start = t + 1;
        }

        if (start <= U)
        {
            // We must only multiply the not-yet-added
            // portion.
            u64 upper =
                U;

            if (accumulated_from <= U)
            {
                // Previously we already have:
                //
                //     accumulated_from ... U
                //
                // Therefore add:
                //
                //     t+1 ... accumulated_from-1
                //
                if (t + 1 >= accumulated_from)
                {
                    start = 0;
                    upper = 0;
                }
                else
                {
                    start = t + 1;
                    upper = accumulated_from - 1;
                }
            }

            if (start != 0 && start <= upper)
            {
                for (u64 x = start; x <= upper; ++x)
                {
                    accumulator =
                        mul_mod(
                            accumulator,
                            static_cast<u128>(x),
                            N
                        );

                    ++product_accumulator_steps;
                }
            }

            accumulated_from = t + 1;
        }

        // ----------------------------------------------------
        // One and only one GCD for this probe.
        // ----------------------------------------------------

        ++gcd_calls;

        u128 g =
            gcd_u128(
                accumulator,
                N
            );

        if (g > 1 && g < N)
        {
            if (g <= static_cast<u128>(UINT64_MAX))
            {
                return
                {
                    true,
                    static_cast<u64>(g),
                    'P'
                };
            }

            return
            {
                true,
                0,
                'G'
            };
        }

        // ----------------------------------------------------
        // MISS
        //
        // Lucas:
        //
        //     p | C(U,U-t)
        //         iff p > t
        //
        // Therefore MISS => p <= t.
        // ----------------------------------------------------

        right = t;
    }

    return {false, 0, 'M'};
}

// ============================================================
// Factorization result
// ============================================================

struct FactorResult
{
    bool success;
    u64 factor;
    char reason;

    u64 probes;
    u64 product_steps;
    u64 gcd_calls;

    int scanned_segments;
    int total_segments;

    double runtime_seconds;
};

// ============================================================
// Factor N
// ============================================================

FactorResult factor_from_n(
    u128 N,
    double assumed_ratio
)
{
    auto begin =
        std::chrono::high_resolution_clock::now();

    u64 sqrt_n =
        isqrt_u64(
            static_cast<u64>(N)
        );

    // --------------------------------------------------------
    // Initial N-only interval:
    //
    //     sqrt(N/C) <= p <= sqrt(N)
    // --------------------------------------------------------

    long double lower_real =
        std::sqrt(
            static_cast<long double>(N) /
            assumed_ratio
        );

    u64 global_L =
        static_cast<u64>(
            std::ceil(lower_real)
        );

    u64 global_R =
        sqrt_n;

    // --------------------------------------------------------
    // Split into safe segments.
    //
    // U <= 1.5L < 2L
    // --------------------------------------------------------

    const long double ALPHA = 1.5L;

    std::vector<std::pair<u64, u64>> segments;

    u64 L = global_L;

    while (L <= global_R)
    {
        u64 U =
            static_cast<u64>(
                std::floor(
                    ALPHA *
                    static_cast<long double>(L)
                )
            );

        if (U < L)
            U = L;

        if (U > global_R)
            U = global_R;

        segments.push_back({L, U});

        if (U == global_R)
            break;

        L = U + 1;
    }

    u64 probes = 0;
    u64 product_steps = 0;
    u64 gcd_calls = 0;

    int scanned_segments = 0;

    // --------------------------------------------------------
    // Search each segment.
    // --------------------------------------------------------

    for (const auto& segment : segments)
    {
        ++scanned_segments;

        const u64 segment_L = segment.first;
        const u64 segment_U = segment.second;

        ProbeResult result =
            search_segment(
                N,
                segment_L,
                segment_U,
                product_steps,
                gcd_calls,
                probes
            );

        if (result.success)
        {
            auto end =
                std::chrono::high_resolution_clock::now();

            double runtime =
                std::chrono::duration<double>(
                    end - begin
                ).count();

            return
            {
                true,
                result.candidate,
                result.reason,
                probes,
                product_steps,
                gcd_calls,
                scanned_segments,
                static_cast<int>(
                    segments.size()
                ),
                runtime
            };
        }
    }

    auto end =
        std::chrono::high_resolution_clock::now();

    double runtime =
        std::chrono::duration<double>(
            end - begin
        ).count();

    return
    {
        false,
        0,
        '?',
        probes,
        product_steps,
        gcd_calls,
        scanned_segments,
        static_cast<int>(
            segments.size()
        ),
        runtime
    };
}

// ============================================================
// Main
// ============================================================

int main()
{
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n\n";

    std::mt19937_64 rng(
        0x154154154ULL
    );

    struct Tier
    {
        std::string name;
        u64 target;
        int cases;
    };

    const std::vector<Tier> tiers =
    {
        {
            "N~1e18",
            1000000000000000000ULL,
            2
        },
        {
            "N~1e19",
            10000000000000000000ULL,
            1
        }
    };

    // Intentionally broad N-only assumption.
    const double ASSUMED_RATIO = 16.0;

    for (const auto& tier : tiers)
    {
        std::cout
            << "TIER "
            << tier.name
            << "\n";

        int successes = 0;
        int failures = 0;

        double total_runtime = 0.0;
        double max_runtime = 0.0;

        double total_probes = 0.0;

        u64 total_product_steps = 0;
        u64 max_product_steps = 0;

        u64 total_gcd_calls = 0;
        u64 max_gcd_calls = 0;

        for (int case_id = 1;
             case_id <= tier.cases;
             ++case_id)
        {
            double requested_ratio =
                1.05 +
                static_cast<double>(
                    rng() % 4500
                ) / 10000.0;

            auto factors =
                make_semiprime(
                    tier.target,
                    requested_ratio,
                    rng
                );

            u64 p = factors.first;
            u64 q = factors.second;

            u128 N =
                static_cast<u128>(p) *
                static_cast<u128>(q);

            double actual_ratio =
                static_cast<double>(q) /
                static_cast<double>(p);

            std::cout
                << "CASE "
                << case_id
                << "\n";

            std::cout
                << "  p="
                << p
                << "\n";

            std::cout
                << "  q="
                << q
                << "\n";

            std::cout
                << "  N="
                << u128_to_string(N)
                << "\n";

            std::cout
                << std::fixed
                << std::setprecision(9)
                << "  ratio="
                << actual_ratio
                << "\n";

            FactorResult result =
                factor_from_n(
                    N,
                    ASSUMED_RATIO
                );

            bool correct =
                result.success &&
                (
                    result.factor == p ||
                    result.factor == q
                );

            if (result.success && correct)
                ++successes;
            else
                ++failures;

            std::cout
                << "  success="
                << (result.success ? 1 : 0)
                << "\n";

            std::cout
                << "  candidate="
                << result.factor
                << "\n";

            std::cout
                << "  correct="
                << (correct ? 1 : 0)
                << "\n";

            std::cout
                << "  reason="
                << result.reason
                << "\n";

            std::cout
                << std::setprecision(9)
                << "  runtime_seconds="
                << result.runtime_seconds
                << "\n";

            std::cout
                << "  probes="
                << result.probes
                << "\n";

            std::cout
                << "  product_steps="
                << result.product_steps
                << "\n";

            std::cout
                << "  gcd_calls="
                << result.gcd_calls
                << "\n";

            std::cout
                << "  segments="
                << result.scanned_segments
                << "/"
                << result.total_segments
                << "\n";

            std::cout << "\n";

            total_runtime +=
                result.runtime_seconds;

            max_runtime =
                std::max(
                    max_runtime,
                    result.runtime_seconds
                );

            total_probes +=
                static_cast<double>(
                    result.probes
                );

            total_product_steps +=
                result.product_steps;

            max_product_steps =
                std::max(
                    max_product_steps,
                    result.product_steps
                );

            total_gcd_calls +=
                result.gcd_calls;

            max_gcd_calls =
                std::max(
                    max_gcd_calls,
                    result.gcd_calls
                );
        }

        std::cout
            << "TIER SUMMARY"
            << "\n";

        std::cout
            << "cases="
            << tier.cases
            << "\n";

        std::cout
            << "successes="
            << successes
            << "\n";

        std::cout
            << "failures="
            << failures
            << "\n";

        std::cout
            << std::fixed
            << std::setprecision(4)
            << "success_rate="
            << static_cast<double>(successes) /
                   static_cast<double>(tier.cases)
            << "\n";

        std::cout
            << "avg_runtime="
            << total_runtime /
                   static_cast<double>(tier.cases)
            << "\n";

        std::cout
            << "max_runtime="
            << max_runtime
            << "\n";

        std::cout
            << "avg_probes="
            << total_probes /
                   static_cast<double>(tier.cases)
            << "\n";

        std::cout
            << "avg_product_steps="
            << static_cast<double>(
                   total_product_steps
               ) /
                   static_cast<double>(tier.cases)
            << "\n";

        std::cout
            << "max_product_steps="
            << max_product_steps
            << "\n";

        std::cout
            << "avg_gcd_calls="
            << static_cast<double>(
                   total_gcd_calls
               ) /
                   static_cast<double>(tier.cases)
            << "\n";

        std::cout
            << "max_gcd_calls="
            << max_gcd_calls
            << "\n\n";
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
