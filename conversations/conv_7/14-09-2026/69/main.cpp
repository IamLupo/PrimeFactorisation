#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        __builtin_sqrtl(
            static_cast<long double>(n)
        )
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(
        static_cast<std::size_t>(limit + 1),
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<long long>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 100000) {
            usable.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        usable.size() - 1
    );

    std::vector<CaseData> result;
    result.reserve(
        static_cast<std::size_t>(count)
    );

    while (
        static_cast<int>(result.size()) < count
    ) {
        const std::size_t i = dist(rng);
        const std::size_t j = dist(rng);

        if (i == j) {
            continue;
        }

        const u64 p =
            static_cast<u64>(
                std::min(
                    usable[i],
                    usable[j]
                )
            );

        const u64 q =
            static_cast<u64>(
                std::max(
                    usable[i],
                    usable[j]
                )
            );

        const u64 N = p * q;

        if (!used_N.insert(N).second) {
            continue;
        }

        result.push_back({
            p,
            q,
            N,
            isqrt_u64(N)
        });
    }

    return result;
}

static u64 abs_i128_to_u64(i128 x) {
    if (x < 0) {
        x = -x;
    }

    return static_cast<u64>(x);
}

static void print_i128(i128 x) {
    if (x == 0) {
        std::cout << "0";
        return;
    }

    if (x < 0) {
        std::cout << '-';
        x = -x;
    }

    std::string out;

    while (x > 0) {
        const unsigned digit =
            static_cast<unsigned>(x % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(
        out.begin(),
        out.end()
    );

    std::cout << out;
}

static u64 gcd_abs(
    i128 x,
    u64 N
) {
    return std::gcd(
        abs_i128_to_u64(x),
        N
    );
}

int main() {
    std::cout << "START EXPERIMENT 358\n";

    constexpr int CASE_COUNT = 2000;
    constexpr int PRIME_LIMIT = 100000;

    /*
     * For each case we also inspect a small deterministic
     * neighborhood around floor(sqrt(Delta)).
     */
    constexpr int NEIGHBOR_RADIUS = 10;

    std::mt19937_64 rng(
        0x358358358ULL
    );

    const std::vector<int> primes =
        sieve_primes(PRIME_LIMIT);

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    u64 delta_positive = 0;
    u64 delta_zero = 0;
    u64 delta_negative = 0;

    u64 delta_square = 0;

    u64 factor_from_A_minus_z = 0;
    u64 factor_from_A_plus_z = 0;

    u64 exact_p_from_minus = 0;
    u64 exact_p_from_plus = 0;

    u64 exact_q_from_minus = 0;
    u64 exact_q_from_plus = 0;

    u64 nontrivial_minus = 0;
    u64 nontrivial_plus = 0;

    u64 any_factor_from_neighborhood = 0;

    u64 p_hits_neighborhood = 0;
    u64 q_hits_neighborhood = 0;

    u64 epsilon_zero = 0;
    u64 epsilon_lt_p = 0;
    u64 epsilon_lt_q = 0;
    u64 epsilon_lt_s = 0;

    u64 min_epsilon = UINT64_MAX;
    u64 max_epsilon = 0;

    u64 min_positive_delta = UINT64_MAX;
    u64 max_positive_delta = 0;

    u64 factor_neighborhood_radius_used =
        static_cast<u64>(NEIGHBOR_RADIUS);

    u64 exact_square_potential_hits = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        /*
         * A = 7s+4.
         *
         * Delta = A^2 - 48N.
         */
        const i128 A =
            static_cast<i128>(7) *
            static_cast<i128>(s)
            +
            static_cast<i128>(4);

        const i128 Delta =
            A * A
            -
            static_cast<i128>(48) *
            static_cast<i128>(N);

        if (Delta < 0) {
            ++delta_negative;
            continue;
        }

        if (Delta == 0) {
            ++delta_zero;
        } else {
            ++delta_positive;
        }

        const u64 delta_u =
            static_cast<u64>(Delta);

        const u64 z =
            isqrt_u64(delta_u);

        const u128 z2 =
            static_cast<u128>(z) *
            static_cast<u128>(z);

        const u128 delta128 =
            static_cast<u128>(delta_u);

        const u64 epsilon =
            static_cast<u64>(
                delta128 - z2
            );

        if (z2 == delta128) {
            ++delta_square;
            ++epsilon_zero;
        }

        min_epsilon =
            std::min(
                min_epsilon,
                epsilon
            );

        max_epsilon =
            std::max(
                max_epsilon,
                epsilon
            );

        min_positive_delta =
            std::min(
                min_positive_delta,
                delta_u
            );

        max_positive_delta =
            std::max(
                max_positive_delta,
                delta_u
            );

        if (epsilon < p) {
            ++epsilon_lt_p;
        }

        if (epsilon < q) {
            ++epsilon_lt_q;
        }

        if (epsilon < s) {
            ++epsilon_lt_s;
        }

        /*
         * ---------------------------------------------------------
         * Difference of squares at z:
         *
         * A^2-z^2 = 48N+epsilon.
         * ---------------------------------------------------------
         */
        const i128 A_minus_z =
            A - static_cast<i128>(z);

        const i128 A_plus_z =
            A + static_cast<i128>(z);

        const u64 g_minus =
            gcd_abs(
                A_minus_z,
                N
            );

        const u64 g_plus =
            gcd_abs(
                A_plus_z,
                N
            );

        if (g_minus != 1) {
            ++nontrivial_minus;
        }

        if (g_plus != 1) {
            ++nontrivial_plus;
        }

        if (g_minus == p) {
            ++exact_p_from_minus;
        }

        if (g_minus == q) {
            ++exact_q_from_minus;
        }

        if (g_plus == p) {
            ++exact_p_from_plus;
        }

        if (g_plus == q) {
            ++exact_q_from_plus;
        }

        if (
            g_minus == p ||
            g_minus == q
        ) {
            ++factor_from_A_minus_z;
        }

        if (
            g_plus == p ||
            g_plus == q
        ) {
            ++factor_from_A_plus_z;
        }

        if (
            g_minus == p ||
            g_minus == q ||
            g_plus == p ||
            g_plus == q
        ) {
            ++exact_square_potential_hits;
        }

        /*
         * ---------------------------------------------------------
         * Search the deterministic neighborhood
         *
         * z+j,  -R <= j <= R.
         *
         * For each y, test
         *
         * gcd(A-y,N)
         * gcd(A+y,N)
         *
         * This tests whether a very small correction from the
         * discriminant square root exposes a factor.
         * ---------------------------------------------------------
         */

        bool case_factor = false;
        bool case_p = false;
        bool case_q = false;

        for (int offset = -NEIGHBOR_RADIUS;
             offset <= NEIGHBOR_RADIUS;
             ++offset) {

            if (offset < 0 &&
                static_cast<u64>(-offset) > z) {
                continue;
            }

            const u64 y =
                offset >= 0
                    ? z + static_cast<u64>(offset)
                    : z - static_cast<u64>(-offset);

            const i128 y128 =
                static_cast<i128>(y);

            const u64 gm =
                gcd_abs(
                    A - y128,
                    N
                );

            const u64 gp =
                gcd_abs(
                    A + y128,
                    N
                );

            if (gm == p || gp == p) {
                case_factor = true;
                case_p = true;
            }

            if (gm == q || gp == q) {
                case_factor = true;
                case_q = true;
            }
        }

        if (case_factor) {
            ++any_factor_from_neighborhood;
        }

        if (case_p) {
            ++p_hits_neighborhood;
        }

        if (case_q) {
            ++q_hits_neighborhood;
        }

        if (case_id < 10) {
            std::cout
                << "\nCASE="
                << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << "\n";

            std::cout
                << "A=";

            print_i128(A);

            std::cout
                << "\nDELTA=";

            print_i128(Delta);

            std::cout
                << "\n";

            if (Delta >= 0) {
                std::cout
                    << "SQRT_FLOOR="
                    << z
                    << "\n";

                std::cout
                    << "EPSILON="
                    << epsilon
                    << "\n";

                std::cout
                    << "GCD_A_MINUS_Z="
                    << g_minus
                    << "\n";

                std::cout
                    << "GCD_A_PLUS_Z="
                    << g_plus
                    << "\n";

                std::cout
                    << "DELTA_IS_SQUARE="
                    << (
                        z2 == delta128
                            ? "YES"
                            : "NO"
                    )
                    << "\n";
            }
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "DELTA_POSITIVE="
        << delta_positive
        << "\n";

    std::cout
        << "DELTA_ZERO="
        << delta_zero
        << "\n";

    std::cout
        << "DELTA_NEGATIVE="
        << delta_negative
        << "\n";

    std::cout
        << "DELTA_SQUARE="
        << delta_square
        << "\n";

    std::cout
        << "EPSILON_ZERO="
        << epsilon_zero
        << "\n";

    std::cout
        << "MIN_EPSILON="
        << min_epsilon
        << "\n";

    std::cout
        << "MAX_EPSILON="
        << max_epsilon
        << "\n";

    std::cout
        << "EPSILON_LT_P="
        << epsilon_lt_p
        << "\n";

    std::cout
        << "EPSILON_LT_Q="
        << epsilon_lt_q
        << "\n";

    std::cout
        << "EPSILON_LT_S="
        << epsilon_lt_s
        << "\n";

    std::cout
        << "NONTRIVIAL_GCD_A_MINUS_Z="
        << nontrivial_minus
        << "\n";

    std::cout
        << "NONTRIVIAL_GCD_A_PLUS_Z="
        << nontrivial_plus
        << "\n";

    std::cout
        << "FACTOR_FROM_A_MINUS_Z="
        << factor_from_A_minus_z
        << "\n";

    std::cout
        << "FACTOR_FROM_A_PLUS_Z="
        << factor_from_A_plus_z
        << "\n";

    std::cout
        << "EXACT_P_FROM_MINUS="
        << exact_p_from_minus
        << "\n";

    std::cout
        << "EXACT_P_FROM_PLUS="
        << exact_p_from_plus
        << "\n";

    std::cout
        << "EXACT_Q_FROM_MINUS="
        << exact_q_from_minus
        << "\n";

    std::cout
        << "EXACT_Q_FROM_PLUS="
        << exact_q_from_plus
        << "\n";

    std::cout
        << "ANY_FACTOR_FROM_NEIGHBORHOOD="
        << any_factor_from_neighborhood
        << "\n";

    std::cout
        << "P_HITS_NEIGHBORHOOD="
        << p_hits_neighborhood
        << "\n";

    std::cout
        << "Q_HITS_NEIGHBORHOOD="
        << q_hits_neighborhood
        << "\n";

    std::cout
        << "NEIGHBOR_RADIUS="
        << factor_neighborhood_radius_used
        << "\n";

    if (delta_square > 0 ||
        factor_from_A_minus_z > 0 ||
        factor_from_A_plus_z > 0) {

        std::cout
            << "DISCRIMINANT_STATUS="
            << "NONTRIVIAL_FACTOR_SIGNAL_FOUND"
            << "\n";
    } else {
        std::cout
            << "DISCRIMINANT_STATUS="
            << "NO_DIRECT_FACTOR_SIGNAL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 358\n";

    return 0;
}
