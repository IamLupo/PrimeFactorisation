#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr int EXPERIMENT = 155;

// ============================================================
// u128 string
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
// integer sqrt
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

    const u64 small[] =
    {
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    };

    for (u64 p : small)
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

    auto mod_mul64 =
        [](u64 a, u64 b, u64 mod) -> u64
        {
            return static_cast<u64>(
                (static_cast<u128>(a) * b) % mod
            );
        };

    auto mod_pow =
        [&](u64 a, u64 e, u64 mod) -> u64
        {
            u64 r = 1;

            while (e)
            {
                if (e & 1ULL)
                    r = mod_mul64(r, a, mod);

                a = mod_mul64(a, a, mod);
                e >>= 1;
            }

            return r;
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
            x = mod_mul64(x, x, n);

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
// random prime
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

        if (candidate >= 3 &&
            is_prime_u64(candidate))
        {
            return candidate;
        }
    }
}

// ============================================================
// semiprime generation
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
// Product mode
// ============================================================

enum class ProductMode
{
    SCALAR,
    PAIR
};

// ============================================================
// Multiply a newly added range.
//
// SCALAR:
//     accumulator = accumulator*x mod N
//
// PAIR:
//     accumulator = accumulator*(x*y) mod N
//
// The pair mode is used only when the intermediate u128
// product is provably safe.
// ============================================================

u64 multiply_range(
    u128& accumulator,
    u128 N,
    u64 start,
    u64 end,
    ProductMode mode
)
{
    if (start > end)
        return 0;

    u64 steps = 0;

    if (mode == ProductMode::SCALAR)
    {
        for (u64 x = start; x <= end; ++x)
        {
            accumulator =
                mul_mod(
                    accumulator,
                    static_cast<u128>(x),
                    N
                );

            ++steps;
        }

        return steps;
    }

    u64 x = start;

    while (x <= end)
    {
        // ----------------------------------------------------
        // Safe pair test:
        //
        // We require
        //
        //     accumulator * x * y < 2^128.
        //
        // Instead of relying on approximate bounds, check
        // the division limit explicitly.
        // ----------------------------------------------------

        if (x < end)
        {
            u64 y = x + 1;

            u128 pair =
                static_cast<u128>(x) *
                static_cast<u128>(y);

            const u128 MAX =
                std::numeric_limits<u128>::max();

            if (accumulator <= MAX / pair)
            {
                accumulator =
                    (accumulator * pair) % N;

                steps += 2;
                x += 2;
                continue;
            }
        }

        accumulator =
            mul_mod(
                accumulator,
                static_cast<u128>(x),
                N
            );

        ++steps;
        ++x;
    }

    return steps;
}

// ============================================================
// Search one safe segment
// ============================================================

struct SearchResult
{
    bool success;
    u64 factor;
    char reason;

    u64 probes;
    u64 product_steps;
    u64 gcd_calls;
};

SearchResult search_segment(
    u128 N,
    u64 L,
    u64 U,
    ProductMode mode
)
{
    if (L > U)
        return {false, 0, 'X', 0, 0, 0};

    if (U >= 2 * L)
        return {false, 0, 'U', 0, 0, 0};

    u64 left = L;
    u64 right = U;

    u128 accumulator = 1;

    u64 accumulated_from = U + 1;

    u64 probes = 0;
    u64 product_steps = 0;
    u64 gcd_calls = 0;

    while (left < right)
    {
        ++probes;

        u64 t =
            left +
            (right - left) / 2;

        u64 new_start = t + 1;
        u64 new_end = U;

        if (accumulated_from <= U)
        {
            if (new_start < accumulated_from)
            {
                new_end =
                    accumulated_from - 1;
            }
            else
            {
                new_start = 1;
                new_end = 0;
            }
        }

        if (new_start <= new_end)
        {
            product_steps +=
                multiply_range(
                    accumulator,
                    N,
                    new_start,
                    new_end,
                    mode
                );
        }

        accumulated_from = t + 1;

        ++gcd_calls;

        u128 g =
            gcd_u128(
                accumulator,
                N
            );

        if (g > 1 && g < N)
        {
            if (g <=
                static_cast<u128>(
                    std::numeric_limits<u64>::max()
                ))
            {
                return
                {
                    true,
                    static_cast<u64>(g),
                    'P',
                    probes,
                    product_steps,
                    gcd_calls
                };
            }
        }

        // MISS => p <= t
        right = t;
    }

    return
    {
        false,
        0,
        'M',
        probes,
        product_steps,
        gcd_calls
    };
}

// ============================================================
// Factorization
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

FactorResult factor_from_n(
    u128 N,
    double assumed_ratio,
    ProductMode mode
)
{
    auto begin =
        std::chrono::high_resolution_clock::now();

    u64 sqrt_n =
        isqrt_u64(
            static_cast<u64>(N)
        );

    u64 global_L =
        static_cast<u64>(
            std::ceil(
                std::sqrt(
                    static_cast<long double>(N) /
                    assumed_ratio
                )
            )
        );

    u64 global_R = sqrt_n;

    const long double ALPHA = 1.5L;

    std::vector<std::pair<u64, u64>>
        segments;

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

    u64 total_probes = 0;
    u64 total_steps = 0;
    u64 total_gcds = 0;

    int scanned_segments = 0;

    for (const auto& segment : segments)
    {
        ++scanned_segments;

        SearchResult result =
            search_segment(
                N,
                segment.first,
                segment.second,
                mode
            );

        total_probes += result.probes;
        total_steps += result.product_steps;
        total_gcds += result.gcd_calls;

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
                result.factor,
                result.reason,
                total_probes,
                total_steps,
                total_gcds,
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
        total_probes,
        total_steps,
        total_gcds,
        scanned_segments,
        static_cast<int>(
            segments.size()
        ),
        runtime
    };
}

// ============================================================
// Run one mode
// ============================================================

void run_mode(
    const std::string& name,
    ProductMode mode,
    u128 N,
    u64 p,
    u64 q
)
{
    FactorResult result =
        factor_from_n(
            N,
            16.0,
            mode
        );

    bool correct =
        result.success &&
        (
            result.factor == p ||
            result.factor == q
        );

    std::cout
        << "  MODE "
        << name
        << "\n";

    std::cout
        << "    success="
        << (result.success ? 1 : 0)
        << "\n";

    std::cout
        << "    candidate="
        << result.factor
        << "\n";

    std::cout
        << "    correct="
        << (correct ? 1 : 0)
        << "\n";

    std::cout
        << "    runtime_seconds="
        << std::fixed
        << std::setprecision(9)
        << result.runtime_seconds
        << "\n";

    std::cout
        << "    probes="
        << result.probes
        << "\n";

    std::cout
        << "    product_steps="
        << result.product_steps
        << "\n";

    std::cout
        << "    gcd_calls="
        << result.gcd_calls
        << "\n";

    std::cout
        << "    segments="
        << result.scanned_segments
        << "/"
        << result.total_segments
        << "\n\n";
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
        0x155155155ULL
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

    for (const auto& tier : tiers)
    {
        std::cout
            << "TIER "
            << tier.name
            << "\n";

        for (int case_id = 1;
             case_id <= tier.cases;
             ++case_id)
        {
            double ratio =
                1.05 +
                static_cast<double>(
                    rng() % 4500
                ) / 10000.0;

            auto factors =
                make_semiprime(
                    tier.target,
                    ratio,
                    rng
                );

            u64 p = factors.first;
            u64 q = factors.second;

            u128 N =
                static_cast<u128>(p) *
                static_cast<u128>(q);

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

            run_mode(
                "SCALAR",
                ProductMode::SCALAR,
                N,
                p,
                q
            );

            run_mode(
                "PAIR",
                ProductMode::PAIR,
                N,
                p,
                q
            );
        }

        std::cout << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
