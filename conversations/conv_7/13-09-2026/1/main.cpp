#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using u64 = std::uint64_t;

static constexpr int EXPERIMENT = 157;

// ============================================================
// Prime test
// ============================================================

bool is_prime_u64(u64 n)
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

// ============================================================
// Random prime near target
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
// Integer square root
// ============================================================

u64 isqrt_u64(u64 n)
{
    u64 x =
        static_cast<u64>(
            std::sqrt(
                static_cast<long double>(n)
            )
        );

    while ((x + 1) <= n / (x + 1))
        ++x;

    while (x > n / x)
        --x;

    return x;
}

// ============================================================
// Make semiprime
// ============================================================

std::pair<u64, u64> make_semiprime(
    u64 target,
    double ratio,
    std::mt19937_64& rng
)
{
    u64 p_center =
        static_cast<u64>(
            std::sqrt(
                static_cast<long double>(target) /
                ratio
            )
        );

    u64 q_center =
        static_cast<u64>(
            static_cast<long double>(p_center) *
            ratio
        );

    const u64 radius = 5000;

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

        return {p, q};
    }
}

// ============================================================
// Lucas theorem divisibility test
//
// Returns true iff
//
//     p | C(m,d)
//
// for prime p.
//
// Lucas:
//     C(m,d) != 0 mod p
//
// iff every base-p digit of d is <= corresponding
// base-p digit of m.
//
// Therefore divisibility occurs iff at least one
// digit comparison fails.
// ============================================================

bool lucas_hit(
    u64 p,
    u64 m,
    u64 d
)
{
    while (m > 0 || d > 0)
    {
        u64 m_digit = m % p;
        u64 d_digit = d % p;

        if (d_digit > m_digit)
            return true;

        m /= p;
        d /= p;
    }

    return false;
}

// ============================================================
// Kummer carry count
//
// Number of carries when adding:
//
//     d + (m-d) = m
//
// in base p.
//
// p | C(m,d) iff carry_count > 0.
// ============================================================

int kummers_carry_count(
    u64 p,
    u64 m,
    u64 d
)
{
    u64 x = d;
    u64 y = m - d;

    int carries = 0;
    u64 carry = 0;

    while (x > 0 || y > 0 || carry > 0)
    {
        u64 xd = x % p;
        u64 yd = y % p;

        u64 sum = xd + yd + carry;

        if (sum >= p)
        {
            ++carries;
            carry = 1;
        }
        else
        {
            carry = 0;
        }

        x /= p;
        y /= p;
    }

    return carries;
}

// ============================================================
// Extract contiguous TRUE intervals.
//
// Input:
//     hit[t] for t = L ... R-1
//
// Output:
//     intervals [a,b] where predicate is TRUE.
// ============================================================

std::vector<std::pair<u64, u64>> hit_intervals(
    u64 L,
    u64 R,
    const std::vector<bool>& hit
)
{
    std::vector<std::pair<u64, u64>> result;

    bool active = false;
    u64 start = 0;

    for (u64 t = L; t < R; ++t)
    {
        bool h =
            hit[
                static_cast<std::size_t>(t - L)
            ];

        if (h && !active)
        {
            active = true;
            start = t;
        }

        bool last =
            (t + 1 == R);

        if (active && (!h || last))
        {
            u64 end;

            if (h && last)
                end = t;
            else
                end = t - 1;

            result.push_back(
                {start, end}
            );

            active = false;
        }
    }

    return result;
}

// ============================================================
// Check whether predicate is exactly a suffix:
//
//     MISS MISS ... MISS HIT HIT ... HIT
//
// This is the required structure for the simple binary
// search rule:
//
//     MISS -> left
//     HIT  -> right
// ============================================================

bool is_monotone_suffix(
    const std::vector<bool>& hit
)
{
    bool seen_hit = false;

    for (bool h : hit)
    {
        if (h)
        {
            seen_hit = true;
        }
        else if (seen_hit)
        {
            return false;
        }
    }

    return true;
}

// ============================================================
// Find first HIT in a monotone suffix
// ============================================================

u64 first_hit(
    u64 L,
    u64 R,
    const std::vector<bool>& hit
)
{
    for (u64 t = L; t < R; ++t)
    {
        if (hit[
                static_cast<std::size_t>(t - L)
            ])
        {
            return t;
        }
    }

    return R;
}

// ============================================================
// Simulate the broad-interval binary search.
//
// Assumption:
//
//     HIT  => p > t
//     MISS => p <= t
//
// This assumption is exactly what we are testing.
// ============================================================

u64 simulated_binary_search(
    u64 L,
    u64 R,
    const std::vector<bool>& hit
)
{
    u64 left = L;
    u64 right = R;

    while (left < right)
    {
        u64 t =
            left +
            (right - left) / 2;

        bool h =
            hit[
                static_cast<std::size_t>(t - L)
            ];

        if (h)
        {
            left = t + 1;
        }
        else
        {
            right = t;
        }
    }

    return left;
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
        0x157157157ULL
    );

    const int CASES = 12;

    // Deliberately broad assumption.
    const double ASSUMED_RATIO = 16.0;

    for (int case_id = 1;
         case_id <= CASES;
         ++case_id)
    {
        double requested_ratio =
            1.05 +
            static_cast<double>(
                rng() % 14000
            ) / 1000.0;

        auto factors =
            make_semiprime(
                10000000000ULL,
                requested_ratio,
                rng
            );

        u64 p = factors.first;
        u64 q = factors.second;

        u64 N = p * q;

        u64 L =
            static_cast<u64>(
                std::ceil(
                    std::sqrt(
                        static_cast<long double>(N) /
                        ASSUMED_RATIO
                    )
                )
            );

        u64 R =
            isqrt_u64(N);

        // ----------------------------------------------------
        // We deliberately use ONE m:
        //
        //     m = R
        //
        // across the entire broad interval [L,R].
        // ----------------------------------------------------

        u64 m = R;

        std::vector<bool> hit;

        hit.reserve(
            static_cast<std::size_t>(
                R - L
            )
        );

        int k_validations = 0;

        bool kummer_consistent = true;

        for (u64 t = L; t < R; ++t)
        {
            u64 d =
                m - t;

            bool h =
                lucas_hit(
                    p,
                    m,
                    d
                );

            hit.push_back(h);

            int carries =
                kummers_carry_count(
                    p,
                    m,
                    d
                );

            bool h_kummer =
                carries > 0;

            if (h != h_kummer)
            {
                kummer_consistent = false;
            }

            ++k_validations;
        }

        auto intervals =
            hit_intervals(
                L,
                R,
                hit
            );

        bool monotone =
            is_monotone_suffix(
                hit
            );

        u64 first =
            first_hit(
                L,
                R,
                hit
            );

        u64 simulated =
            simulated_binary_search(
                L,
                R,
                hit
            );

        // ----------------------------------------------------
        // Actual quotient at m=R.
        // ----------------------------------------------------

        u64 quotient =
            m / p;

        u64 remainder =
            m % p;

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
            << N
            << "\n";

        std::cout
            << std::fixed
            << std::setprecision(6)
            << "  actual_ratio="
            << static_cast<double>(q) /
               static_cast<double>(p)
            << "\n";

        std::cout
            << "  L="
            << L
            << "\n";

        std::cout
            << "  R="
            << R
            << "\n";

        std::cout
            << "  m="
            << m
            << "\n";

        std::cout
            << "  floor_m_over_p="
            << quotient
            << "\n";

        std::cout
            << "  m_mod_p="
            << remainder
            << "\n";

        std::cout
            << "  candidate_width="
            << (R - L)
            << "\n";

        std::cout
            << "  lucas_tests="
            << hit.size()
            << "\n";

        std::cout
            << "  k_validations="
            << k_validations
            << "\n";

        std::cout
            << "  k_lucas_consistent="
            << (kummer_consistent ? 1 : 0)
            << "\n";

        std::cout
            << "  monotone_suffix="
            << (monotone ? 1 : 0)
            << "\n";

        std::cout
            << "  first_hit="
            << first
            << "\n";

        std::cout
            << "  simulated_binary_result="
            << simulated
            << "\n";

        std::cout
            << "  binary_result_correct="
            << (simulated == p ? 1 : 0)
            << "\n";

        std::cout
            << "  hit_interval_count="
            << intervals.size()
            << "\n";

        std::cout
            << "  hit_intervals=";

        if (intervals.empty())
        {
            std::cout << "NONE";
        }
        else
        {
            for (std::size_t i = 0;
                 i < intervals.size();
                 ++i)
            {
                if (i != 0)
                    std::cout << ",";

                std::cout
                    << "["
                    << intervals[i].first
                    << ","
                    << intervals[i].second
                    << "]";
            }
        }

        std::cout
            << "\n";

        std::cout << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
