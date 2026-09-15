#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static constexpr int EXPERIMENT_ID = 152;

static constexpr long double RATIO_MIN = 1.05L;
static constexpr long double RATIO_MAX = 1.50L;
static constexpr long double ALPHA = 1.50L;

static constexpr u64 BLOCK_SIZE = 256;


/* ============================================================
 * u128 helpers
 * ============================================================ */

std::string u128_to_string(u128 value)
{
    if (value == 0)
    {
        return "0";
    }

    std::string result;

    while (value > 0)
    {
        const unsigned digit =
            static_cast<unsigned>(value % 10);

        result.push_back(
            static_cast<char>('0' + digit)
        );

        value /= 10;
    }

    std::reverse(
        result.begin(),
        result.end()
    );

    return result;
}


u128 string_to_u128(const std::string& text)
{
    u128 value = 0;

    for (char c : text)
    {
        if (c < '0' || c > '9')
        {
            throw std::runtime_error(
                "invalid decimal u128"
            );
        }

        value =
            value * 10
            + static_cast<unsigned>(c - '0');
    }

    return value;
}


u128 gcd_u128(u128 a, u128 b)
{
    while (b != 0)
    {
        const u128 r = a % b;
        a = b;
        b = r;
    }

    return a;
}


u128 mul_mod(
    u128 a,
    u128 b,
    u128 n
)
{
    /*
     * In this algorithm:

        b <= sqrt(N) <= 1e11
        a < N <= 1e22

       therefore:

        a*b <= 1e33

       which safely fits in uint128.
    */

    return (a * b) % n;
}


/* ============================================================
 * uint64 primality
 * ============================================================ */

u64 mod_pow_u64(
    u64 a,
    u64 e,
    u64 n
)
{
    u64 result = 1;

    while (e > 0)
    {
        if (e & 1)
        {
            result =
                static_cast<u64>(
                    (
                        static_cast<u128>(result)
                        * a
                    ) % n
                );
        }

        a =
            static_cast<u64>(
                (
                    static_cast<u128>(a)
                    * a
                ) % n
            );

        e >>= 1;
    }

    return result;
}


bool is_prime_u64(u64 n)
{
    if (n < 2)
    {
        return false;
    }

    static constexpr u64 small_primes[] =
    {
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31, 37
    };

    for (u64 p : small_primes)
    {
        if (n == p)
        {
            return true;
        }

        if (n % p == 0)
        {
            return false;
        }
    }

    u64 d = n - 1;
    int s = 0;

    while ((d & 1) == 0)
    {
        d >>= 1;
        ++s;
    }

    static constexpr u64 bases[] =
    {
        2, 3, 5, 7, 11, 13, 17
    };

    for (u64 a : bases)
    {
        if (a >= n)
        {
            continue;
        }

        u64 x =
            mod_pow_u64(
                a,
                d,
                n
            );

        if (
            x == 1
            || x == n - 1
        )
        {
            continue;
        }

        bool passed = false;

        for (int r = 1; r < s; ++r)
        {
            x =
                static_cast<u64>(
                    (
                        static_cast<u128>(x)
                        * x
                    ) % n
                );

            if (x == n - 1)
            {
                passed = true;
                break;
            }
        }

        if (!passed)
        {
            return false;
        }
    }

    return true;
}


u64 next_prime(u64 x)
{
    if (x <= 2)
    {
        return 2;
    }

    if ((x & 1) == 0)
    {
        ++x;
    }

    while (!is_prime_u64(x))
    {
        x += 2;
    }

    return x;
}


/* ============================================================
 * sqrt
 * ============================================================ */

u64 integer_sqrt_u128(u128 n)
{
    u64 low = 0;
    u64 high = 200000000000ULL;

    while (low <= high)
    {
        const u64 mid =
            low + (high - low) / 2;

        const u128 square =
            static_cast<u128>(mid) * mid;

        if (square == n)
        {
            return mid;
        }

        if (square < n)
        {
            low = mid + 1;
        }
        else
        {
            high = mid - 1;
        }
    }

    return high;
}


/* ============================================================
 * semiprime
 * ============================================================ */

struct Semiprime
{
    u64 p;
    u64 q;
    u128 n;
};


Semiprime generate_semiprime(
    u128 target,
    std::mt19937_64& rng
)
{
    const long double target_ld =
        static_cast<long double>(target);

    while (true)
    {
        const long double ratio =
            RATIO_MIN
            + std::generate_canonical<
                long double,
                64
              >(rng)
              * (RATIO_MAX - RATIO_MIN);

        const u64 p_estimate =
            static_cast<u64>(
                std::sqrt(
                    target_ld / ratio
                )
            );

        const u64 spread =
            std::max<u64>(
                100,
                p_estimate / 100
            );

        std::uniform_int_distribution<u64>
            distribution(
                p_estimate > spread
                    ? p_estimate - spread
                    : 3,
                p_estimate + spread
            );

        const u64 p =
            next_prime(
                distribution(rng)
            );

        const u64 q_estimate =
            static_cast<u64>(
                target / p
            );

        const u64 q =
            next_prime(q_estimate);

        if (q <= p)
        {
            continue;
        }

        const u128 n =
            static_cast<u128>(p) * q;

        const long double actual_ratio =
            static_cast<long double>(q)
            / static_cast<long double>(p);

        if (
            n < target * 75 / 100
            || n > target * 150 / 100
        )
        {
            continue;
        }

        if (
            actual_ratio < RATIO_MIN
            || actual_ratio > RATIO_MAX
        )
        {
            continue;
        }

        return {
            p,
            q,
            n
        };
    }
}


/* ============================================================
 * intervals
 * ============================================================ */

struct Interval
{
    u64 lower;
    u64 upper;
};


Interval make_initial_interval(u128 n)
{
    const long double n_ld =
        static_cast<long double>(n);

    const u64 lower =
        static_cast<u64>(
            std::ceil(
                std::sqrt(
                    n_ld / RATIO_MAX
                )
            )
        );

    const u64 upper =
        integer_sqrt_u128(n);

    return {
        lower,
        upper
    };
}


std::vector<Interval>
subdivide_interval(
    u64 lower,
    u64 upper
)
{
    std::vector<Interval> result;

    u64 current = lower;

    while (current <= upper)
    {
        u64 segment_upper =
            static_cast<u64>(
                ALPHA
                * static_cast<long double>(current)
            );

        segment_upper =
            std::min(
                segment_upper,
                upper
            );

        if (segment_upper <= current)
        {
            segment_upper = current + 1;
        }

        result.push_back(
            {
                current,
                segment_upper
            }
        );

        if (segment_upper >= upper)
        {
            break;
        }

        current =
            segment_upper + 1;
    }

    return result;
}


/* ============================================================
 * block numerator gcd
 * ============================================================ */

struct ProbeResult
{
    u128 gcd;
    u64 product_steps;
    u64 gcd_calls;
};


ProbeResult binomial_gcd_block(
    u128 n,
    u64 m,
    u64 d
)
{
    if (d == 0)
    {
        return {
            1,
            0,
            0
        };
    }

    const u64 start =
        m - d + 1;

    u128 accumulator = 1;

    u64 product_steps = 0;
    u64 gcd_calls = 0;

    u64 position = start;

    while (position <= m)
    {
        const u64 block_end =
            std::min(
                m,
                position + BLOCK_SIZE - 1
            );

        /*
         * Multiply the complete block modulo N.
         *
         * No gcd is performed for individual terms.
         */
        for (
            u64 x = position;
            x <= block_end;
            ++x
        )
        {
            accumulator =
                mul_mod(
                    accumulator,
                    static_cast<u128>(x),
                    n
                );

            ++product_steps;
        }

        const u128 g =
            gcd_u128(
                accumulator,
                n
            );

        ++gcd_calls;

        if (
            g > 1
            && g < n
        )
        {
            return {
                g,
                product_steps,
                gcd_calls
            };
        }

        if (g == n)
        {
            return {
                n,
                product_steps,
                gcd_calls
            };
        }

        position =
            block_end + 1;
    }

    return {
        gcd_u128(
            accumulator,
            n
        ),
        product_steps,
        gcd_calls + 1
    };
}


/* ============================================================
 * factorization
 * ============================================================ */

struct FactorResult
{
    bool success;
    u64 candidate;
    std::string reason;

    u64 probes;
    u64 product_steps;
    u64 gcd_calls;

    u64 segments_tested;
    u64 total_segments;

    double runtime_seconds;
};


FactorResult factor_n(
    u128 n,
    u64 p,
    u64 q
)
{
    const Interval initial =
        make_initial_interval(n);

    const auto segments =
        subdivide_interval(
            initial.lower,
            initial.upper
        );

    u64 probes = 0;
    u64 product_steps = 0;
    u64 gcd_calls = 0;

    const auto begin =
        std::chrono::steady_clock::now();

    for (
        std::size_t segment_index = 0;
        segment_index < segments.size();
        ++segment_index
    )
    {
        u64 lower =
            segments[segment_index].lower;

        u64 upper =
            segments[segment_index].upper;

        while (lower < upper)
        {
            const u64 midpoint =
                lower
                + (upper - lower) / 2;

            const u64 m = upper;

            const u64 d =
                upper - midpoint;

            if (
                d >= lower
                || m >= 2 * lower
                || static_cast<u128>(m) >= q
            )
            {
                const auto finish =
                    std::chrono::steady_clock::now();

                return {
                    false,
                    0,
                    "unsafe",
                    probes,
                    product_steps,
                    gcd_calls,
                    static_cast<u64>(
                        segment_index + 1
                    ),
                    static_cast<u64>(
                        segments.size()
                    ),
                    std::chrono::duration<double>(
                        finish - begin
                    ).count()
                };
            }

            ++probes;

            const ProbeResult probe =
                binomial_gcd_block(
                    n,
                    m,
                    d
                );

            product_steps +=
                probe.product_steps;

            gcd_calls +=
                probe.gcd_calls;

            if (probe.gcd == p)
            {
                const auto finish =
                    std::chrono::steady_clock::now();

                return {
                    true,
                    p,
                    "P",
                    probes,
                    product_steps,
                    gcd_calls,
                    static_cast<u64>(
                        segment_index + 1
                    ),
                    static_cast<u64>(
                        segments.size()
                    ),
                    std::chrono::duration<double>(
                        finish - begin
                    ).count()
                };
            }

            if (probe.gcd == q)
            {
                const auto finish =
                    std::chrono::steady_clock::now();

                return {
                    true,
                    q,
                    "Q",
                    probes,
                    product_steps,
                    gcd_calls,
                    static_cast<u64>(
                        segment_index + 1
                    ),
                    static_cast<u64>(
                        segments.size()
                    ),
                    std::chrono::duration<double>(
                        finish - begin
                    ).count()
                };
            }

            if (probe.gcd == n)
            {
                const auto finish =
                    std::chrono::steady_clock::now();

                return {
                    false,
                    0,
                    "N",
                    probes,
                    product_steps,
                    gcd_calls,
                    static_cast<u64>(
                        segment_index + 1
                    ),
                    static_cast<u64>(
                        segments.size()
                    ),
                    std::chrono::duration<double>(
                        finish - begin
                    ).count()
                };
            }

            if (
                probe.gcd != 1
                && probe.gcd != p
                && probe.gcd != q
            )
            {
                const auto finish =
                    std::chrono::steady_clock::now();

                return {
                    false,
                    0,
                    "OTHER",
                    probes,
                    product_steps,
                    gcd_calls,
                    static_cast<u64>(
                        segment_index + 1
                    ),
                    static_cast<u64>(
                        segments.size()
                    ),
                    std::chrono::duration<double>(
                        finish - begin
                    ).count()
                };
            }

            // MISS:
            upper = midpoint;
        }

        if (
            lower == p
            || lower == q
        )
        {
            const auto finish =
                std::chrono::steady_clock::now();

            return {
                true,
                lower,
                "binary",
                probes,
                product_steps,
                gcd_calls,
                static_cast<u64>(
                    segment_index + 1
                ),
                static_cast<u64>(
                    segments.size()
                ),
                std::chrono::duration<double>(
                    finish - begin
                ).count()
            };
        }
    }

    const auto finish =
        std::chrono::steady_clock::now();

    return {
        false,
        0,
        "no-factor",
        probes,
        product_steps,
        gcd_calls,
        static_cast<u64>(
            segments.size()
        ),
        static_cast<u64>(
            segments.size()
        ),
        std::chrono::duration<double>(
            finish - begin
        ).count()
    };
}


/* ============================================================
 * main
 * ============================================================ */

int main()
{
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT_ID
        << "\n\n";

    /*
     * Start with two 10^18 cases.
     *
     * Once this shows the expected speedup, the next tier can
     * safely be increased.
     */
    const std::vector<std::pair<std::string, int>>
        tiers =
    {
        {"1000000000000000000", 2},
        {"10000000000000000000", 1}
    };

    std::mt19937_64 rng(
        EXPERIMENT_ID
    );

    for (const auto& tier : tiers)
    {
        const u128 target =
            string_to_u128(
                tier.first
            );

        std::cout
            << "TIER N~"
            << tier.first
            << "\n\n";

        double total_runtime = 0.0;
        double max_runtime = 0.0;

        u64 successes = 0;

        u64 total_probes = 0;
        u64 total_steps = 0;
        u64 total_gcd_calls = 0;

        u64 max_steps = 0;
        u64 max_gcd_calls = 0;

        for (
            int case_index = 1;
            case_index <= tier.second;
            ++case_index
        )
        {
            const Semiprime sp =
                generate_semiprime(
                    target,
                    rng
                );

            const FactorResult result =
                factor_n(
                    sp.n,
                    sp.p,
                    sp.q
                );

            const bool correct =
                result.success
                && (
                    result.candidate == sp.p
                    || result.candidate == sp.q
                );

            if (correct)
            {
                ++successes;
            }

            total_runtime +=
                result.runtime_seconds;

            max_runtime =
                std::max(
                    max_runtime,
                    result.runtime_seconds
                );

            total_probes +=
                result.probes;

            total_steps +=
                result.product_steps;

            total_gcd_calls +=
                result.gcd_calls;

            max_steps =
                std::max(
                    max_steps,
                    result.product_steps
                );

            max_gcd_calls =
                std::max(
                    max_gcd_calls,
                    result.gcd_calls
                );

            std::cout
                << "CASE "
                << case_index
                << "\n";

            std::cout
                << "  p="
                << sp.p
                << "\n";

            std::cout
                << "  q="
                << sp.q
                << "\n";

            std::cout
                << "  N="
                << u128_to_string(sp.n)
                << "\n";

            std::cout
                << "  ratio="
                << std::fixed
                << std::setprecision(9)
                << static_cast<long double>(sp.q)
                    / static_cast<long double>(sp.p)
                << "\n";

            std::cout
                << "  success="
                << result.success
                << "\n";

            std::cout
                << "  candidate="
                << result.candidate
                << "\n";

            std::cout
                << "  correct="
                << correct
                << "\n";

            std::cout
                << "  reason="
                << result.reason
                << "\n";

            std::cout
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
                << result.segments_tested
                << "/"
                << result.total_segments
                << "\n\n";
        }

        std::cout
            << "TIER SUMMARY\n";

        std::cout
            << "cases="
            << tier.second
            << "\n";

        std::cout
            << "successes="
            << successes
            << "\n";

        std::cout
            << "failures="
            << tier.second - successes
            << "\n";

        std::cout
            << "success_rate="
            << std::fixed
            << std::setprecision(4)
            << static_cast<double>(
                   successes
               )
               / tier.second
            << "\n";

        std::cout
            << "avg_runtime="
            << total_runtime
               / tier.second
            << "\n";

        std::cout
            << "max_runtime="
            << max_runtime
            << "\n";

        std::cout
            << "avg_probes="
            << static_cast<double>(
                   total_probes
               )
               / tier.second
            << "\n";

        std::cout
            << "avg_product_steps="
            << static_cast<double>(
                   total_steps
               )
               / tier.second
            << "\n";

        std::cout
            << "max_product_steps="
            << max_steps
            << "\n";

        std::cout
            << "avg_gcd_calls="
            << static_cast<double>(
                   total_gcd_calls
               )
               / tier.second
            << "\n";

        std::cout
            << "max_gcd_calls="
            << max_gcd_calls
            << "\n\n";
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT_ID
        << "\n";

    return 0;
}
