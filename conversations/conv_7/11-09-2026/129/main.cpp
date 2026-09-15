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

static constexpr int EXPERIMENT = 156;

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
// u128 gcd
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
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
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

    auto mod_mul =
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
                    result = mod_mul(result, a, mod);

                a = mod_mul(a, a, mod);
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
            x = mod_mul(x, x, n);

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

        if (candidate >= 3 &&
            is_prime_u64(candidate))
        {
            return candidate;
        }
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
// Montgomery arithmetic
//
// N is odd because N = p*q for odd primes.
//
// R = 2^64.
//
// Montgomery representation:
//     x_bar = x * R mod N
//
// montgomery multiplication:
//     mont_mul(a_bar,b_bar)
//         = a*b*R^(-1) mod N
// ============================================================

struct Montgomery64
{
    u64 n;
    u64 n_prime;
    u64 r_mod_n;
    u64 r2_mod_n;

    // --------------------------------------------------------
    // Compute -N^{-1} mod 2^64.
    //
    // Newton iteration:
    //     x <- x * (2 - N*x)
    // --------------------------------------------------------

    static u64 compute_n_prime(u64 n)
    {
        u64 x = n;

        // Six Newton steps are enough for 64 bits.
        x *= 2ULL - n * x;
        x *= 2ULL - n * x;
        x *= 2ULL - n * x;
        x *= 2ULL - n * x;
        x *= 2ULL - n * x;
        x *= 2ULL - n * x;

        return ~x + 1ULL;
    }

    // --------------------------------------------------------
    // Compute 2^64 mod n and 2^128 mod n.
    // --------------------------------------------------------

    static u64 compute_r_mod_n(u64 n)
    {
        u128 R = (static_cast<u128>(1) << 64);
        return static_cast<u64>(R % n);
    }

    static u64 compute_r2_mod_n(
        u64 n,
        u64 r_mod_n
    )
    {
        return static_cast<u64>(
            (
                static_cast<u128>(r_mod_n) *
                static_cast<u128>(r_mod_n)
            ) % n
        );
    }

    explicit Montgomery64(u64 modulus)
        : n(modulus)
    {
        if ((n & 1ULL) == 0)
        {
            std::cerr
                << "ERROR: Montgomery modulus must be odd.\n";

            std::exit(1);
        }

        n_prime =
            compute_n_prime(n);

        r_mod_n =
            compute_r_mod_n(n);

        r2_mod_n =
            compute_r2_mod_n(
                n,
                r_mod_n
            );
    }

    // --------------------------------------------------------
    // Montgomery reduction.
    //
    // Input:
    //     t < N^2 < 2^128
    //
    // Output:
    //     t * R^(-1) mod N
    // --------------------------------------------------------

    u64 reduce(u128 t) const
    {
        u64 m =
            static_cast<u64>(t) *
            n_prime;

        u128 u =
            t +
            static_cast<u128>(m) *
            static_cast<u128>(n);

        u64 result =
            static_cast<u64>(u >> 64);

        if (result >= n)
            result -= n;

        return result;
    }

    // --------------------------------------------------------
    // a,b are already Montgomery residues.
    // --------------------------------------------------------

    u64 mul(u64 a, u64 b) const
    {
        return reduce(
            static_cast<u128>(a) *
            static_cast<u128>(b)
        );
    }

    // --------------------------------------------------------
    // Convert normal x into Montgomery representation.
    //
    // x_bar = x * R^2 * R^-1 = x*R mod N
    // --------------------------------------------------------

    u64 to_mont(u64 x) const
    {
        return reduce(
            static_cast<u128>(x) *
            static_cast<u128>(r2_mod_n)
        );
    }

    // --------------------------------------------------------
    // Convert Montgomery x back to normal representation.
    // --------------------------------------------------------

    u64 from_mont(u64 x) const
    {
        return reduce(
            static_cast<u128>(x)
        );
    }
};

// ============================================================
// Verify Montgomery implementation.
// ============================================================

bool verify_montgomery(
    u64 n,
    const Montgomery64& mont
)
{
    std::mt19937_64 rng(
        0x156ABCDEFULL
    );

    for (int i = 0; i < 10000; ++i)
    {
        u64 a =
            1 +
            rng() % (n - 1);

        u64 b =
            1 +
            rng() % (n - 1);

        u64 a_m =
            mont.to_mont(a);

        u64 b_m =
            mont.to_mont(b);

        u64 c_m =
            mont.mul(
                a_m,
                b_m
            );

        u64 c =
            mont.from_mont(c_m);

        u64 expected =
            static_cast<u64>(
                (
                    static_cast<u128>(a) *
                    static_cast<u128>(b)
                ) % n
            );

        if (c != expected)
            return false;
    }

    return true;
}

// ============================================================
// Result
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

// ============================================================
// Search using ordinary paired modular multiplication
//
// For x,y <= floor(sqrt(N)):
//
//     x*y < N
//
// because x != y and both are <= floor(sqrt(N)).
//
// Therefore pair itself requires no modulo.
// ============================================================

SearchResult search_pair_u128(
    u128 N,
    u64 L,
    u64 U
)
{
    if (L > U)
        return {false, 0, 'X', 0, 0, 0};

    if (U >= 2 * L)
        return {false, 0, 'U', 0, 0, 0};

    u64 left = L;
    u64 right = U;

    u64 accumulated_from = U + 1;

    u128 accumulator = 1;

    u64 probes = 0;
    u64 product_steps = 0;
    u64 gcd_calls = 0;

    while (left < right)
    {
        ++probes;

        u64 t =
            left +
            (right - left) / 2;

        u64 start = t + 1;
        u64 end = U;

        if (accumulated_from <= U)
        {
            if (start < accumulated_from)
            {
                end =
                    accumulated_from - 1;
            }
            else
            {
                start = 1;
                end = 0;
            }
        }

        // ----------------------------------------------------
        // Pair adjacent terms.
        // ----------------------------------------------------

        u64 x = start;

        while (x <= end)
        {
            if (x < end)
            {
                u64 y = x + 1;

                u64 pair =
                    x * y;

                accumulator =
                    (
                        accumulator *
                        static_cast<u128>(pair)
                    ) % N;

                product_steps += 2;

                x += 2;
            }
            else
            {
                accumulator =
                    (
                        accumulator *
                        static_cast<u128>(x)
                    ) % N;

                ++product_steps;
                ++x;
            }
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
                    UINT64_MAX
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

        // MISS => p <= t.
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
// Search using Montgomery multiplication.
//
// accumulator is stored in Montgomery form.
//
// The paired factor:
//
//     x*y
//
// is < N and therefore directly fits as a normal u64 residue.
//
// We convert the pair to Montgomery representation once,
// then use Montgomery multiplication with no division by N.
// ============================================================

SearchResult search_pair_mont(
    u64 N,
    u64 L,
    u64 U,
    const Montgomery64& mont
)
{
    if (L > U)
        return {false, 0, 'X', 0, 0, 0};

    if (U >= 2 * L)
        return {false, 0, 'U', 0, 0, 0};

    u64 left = L;
    u64 right = U;

    u64 accumulated_from = U + 1;

    // Montgomery representation of 1.
    u64 accumulator =
        mont.to_mont(1);

    u64 probes = 0;
    u64 product_steps = 0;
    u64 gcd_calls = 0;

    while (left < right)
    {
        ++probes;

        u64 t =
            left +
            (right - left) / 2;

        u64 start = t + 1;
        u64 end = U;

        if (accumulated_from <= U)
        {
            if (start < accumulated_from)
            {
                end =
                    accumulated_from - 1;
            }
            else
            {
                start = 1;
                end = 0;
            }
        }

        // ----------------------------------------------------
        // Paired Montgomery product.
        // ----------------------------------------------------

        u64 x = start;

        while (x <= end)
        {
            if (x < end)
            {
                u64 y = x + 1;

                // Since x,y <= floor(sqrt(N))
                // and x != y:
                //
                //     x*y < N
                //
                u64 pair =
                    x * y;

                u64 pair_mont =
                    mont.to_mont(pair);

                accumulator =
                    mont.mul(
                        accumulator,
                        pair_mont
                    );

                product_steps += 2;

                x += 2;
            }
            else
            {
                u64 x_mont =
                    mont.to_mont(x);

                accumulator =
                    mont.mul(
                        accumulator,
                        x_mont
                    );

                ++product_steps;
                ++x;
            }
        }

        accumulated_from = t + 1;

        ++gcd_calls;

        u64 normal =
            mont.from_mont(
                accumulator
            );

        u128 g =
            gcd_u128(
                static_cast<u128>(normal),
                static_cast<u128>(N)
            );

        if (g > 1 && g < N)
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

        // MISS => p <= t.
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
// Factor result
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
// factor using selected mode
// ============================================================

FactorResult factor_from_n(
    u128 N128,
    double assumed_ratio,
    bool use_montgomery
)
{
    auto begin =
        std::chrono::high_resolution_clock::now();

    if (N128 > static_cast<u128>(UINT64_MAX))
    {
        return
        {
            false,
            0,
            'N',
            0,
            0,
            0,
            0,
            0,
            0.0
        };
    }

    u64 N =
        static_cast<u64>(N128);

    u64 sqrt_n =
        isqrt_u64(N);

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

    Montgomery64 mont(N);

    u64 probes = 0;
    u64 product_steps = 0;
    u64 gcd_calls = 0;

    int scanned_segments = 0;

    for (const auto& segment : segments)
    {
        ++scanned_segments;

        SearchResult result;

        if (use_montgomery)
        {
            result =
                search_pair_mont(
                    N,
                    segment.first,
                    segment.second,
                    mont
                );
        }
        else
        {
            result =
                search_pair_u128(
                    N128,
                    segment.first,
                    segment.second
                );
        }

        probes += result.probes;
        product_steps += result.product_steps;
        gcd_calls += result.gcd_calls;

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
// print result
// ============================================================

void print_result(
    const std::string& mode,
    const FactorResult& result,
    u64 p,
    u64 q
)
{
    bool correct =
        result.success &&
        (
            result.factor == p ||
            result.factor == q
        );

    std::cout
        << "  MODE "
        << mode
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
        << "    reason="
        << result.reason
        << "\n";

    std::cout
        << std::fixed
        << std::setprecision(9)
        << "    runtime_seconds="
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
// main
// ============================================================

int main()
{
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n\n";

    std::mt19937_64 rng(
        0x156156156ULL
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

            // ------------------------------------------------
            // Montgomery self-test.
            // ------------------------------------------------

            Montgomery64 mont(
                static_cast<u64>(N)
            );

            bool mont_ok =
                verify_montgomery(
                    static_cast<u64>(N),
                    mont
                );

            std::cout
                << "  montgomery_selftest="
                << (mont_ok ? 1 : 0)
                << "\n";

            if (!mont_ok)
            {
                std::cout
                    << "ERROR: Montgomery self-test failed.\n";

                return 1;
            }

            // ------------------------------------------------
            // Baseline.
            // ------------------------------------------------

            FactorResult normal =
                factor_from_n(
                    N,
                    16.0,
                    false
                );

            print_result(
                "PAIR_U128",
                normal,
                p,
                q
            );

            // ------------------------------------------------
            // Montgomery.
            // ------------------------------------------------

            FactorResult montgomery =
                factor_from_n(
                    N,
                    16.0,
                    true
                );

            print_result(
                "PAIR_MONTGOMERY",
                montgomery,
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