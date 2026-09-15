#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

static mpz_class mpz_from_u64(u64 value) {
    return mpz_class(std::to_string(value));
}

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        std::sqrt(static_cast<long double>(n))
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static u64 integer_nth_root_4(u64 n) {
    u64 x = static_cast<u64>(
        std::sqrt(
            std::sqrt(
                static_cast<long double>(n)
            )
        )
    );

    auto leq = [n](u64 y) -> bool {
        const u64 a = y * y;
        if (a == 0) {
            return true;
        }

        if (a > n / a) {
            return false;
        }

        return a * a <= n;
    };

    while ((x + 1) > x &&
           leq(x + 1)) {
        ++x;
    }

    while (x > 0 && !leq(x)) {
        --x;
    }

    return x;
}

static bool is_prime_u64(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d <= n / d; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

static u64 random_prime(
    std::mt19937_64& rng,
    u64 lo,
    u64 hi
) {
    std::uniform_int_distribution<u64> dist(lo, hi);

    while (true) {
        u64 x = dist(rng);

        if (x < 2) {
            x = 2;
        }

        if (x > 2 && x % 2 == 0) {
            ++x;
        }

        while (x <= hi && !is_prime_u64(x)) {
            x += 2;
        }

        if (x >= lo && x <= hi) {
            return x;
        }
    }
}

static CaseData make_case(
    std::mt19937_64& rng
) {
    while (true) {
        u64 p = random_prime(rng, 5, 3000);
        u64 q = random_prime(rng, 3001, 10000);

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        const u64 N = p * q;
        const u64 s = isqrt_u64(N);

        if (s * s > N) {
            continue;
        }

        if ((s + 1) * (s + 1) <= N) {
            continue;
        }

        return {p, q, N, s};
    }
}

/*
    Smallest positive root of

        2t^2 == 1 (mod p)

    This is only used as the ground truth.
*/
static u64 smallest_p_root(u64 p) {
    for (u64 t = 1; t < p; ++t) {
        if (((2ULL * t * t) % p) == 1) {
            return t;
        }
    }

    return 0;
}

static bool direct_factor_test(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        2 * T * T - 1;

    const mpz_class g =
        gcd(N, value);

    return g > 1 && g < c.N;
}

/*
    Construct sparse schedule 1,2,4,8,...

    We also include the largest point below each power-of-two
    interval to avoid an excessively coarse schedule.
*/
static std::vector<u64> build_power_schedule(u64 s) {
    std::set<u64> values;

    u64 x = 1;

    while (x < s) {
        values.insert(x);

        if (x > s / 2) {
            break;
        }

        x *= 2;
    }

    return std::vector<u64>(
        values.begin(),
        values.end()
    );
}

/*
    Schedule based on floor(s/k), k=1,...,K.
*/
static std::vector<u64> build_harmonic_schedule(
    u64 s,
    u64 K
) {
    std::set<u64> values;

    for (u64 k = 1; k <= K; ++k) {
        const u64 x = s / k;

        if (x >= 1 && x < s) {
            values.insert(x);
        }
    }

    return std::vector<u64>(
        values.begin(),
        values.end()
    );
}

/*
    Mixed schedule:
      powers of two,
      floor(s/k),
      small integers.

    This is intentionally sparse.
*/
static std::vector<u64> build_mixed_schedule(
    u64 s
) {
    std::set<u64> values;

    /*
        Small deterministic prefix.
    */
    const u64 prefix =
        std::min<u64>(64, s);

    for (u64 t = 1; t <= prefix; ++t) {
        values.insert(t);
    }

    /*
        Powers of two.
    */
    for (u64 t = 1; t < s; ) {
        values.insert(t);

        if (t > s / 2) {
            break;
        }

        t *= 2;
    }

    /*
        Harmonic landmarks.
    */
    const u64 K =
        std::min<u64>(
            256,
            s
        );

    for (u64 k = 1; k <= K; ++k) {
        const u64 t = s / k;

        if (t >= 1 && t < s) {
            values.insert(t);
        }
    }

    return std::vector<u64>(
        values.begin(),
        values.end()
    );
}

/*
    Test a schedule against the known smallest p-root.

    The schedule succeeds if it contains a value >= the root in the
    exact residue progression? No: we require testing the exact root
    itself, because gcd(N,2t²-1) only succeeds at the root positions.
*/
static bool schedule_contains_root(
    const std::vector<u64>& schedule,
    u64 root
) {
    return std::binary_search(
        schedule.begin(),
        schedule.end(),
        root
    );
}

int main() {
    constexpr int EXPERIMENT = 387;
    constexpr int CASES = 2000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38720260914ULL);

    u64 generic_cases = 0;
    u64 no_p_root_cases = 0;

    u64 total_linear_half = 0;
    u64 total_linear_quarter = 0;

    u64 power_schedule_total = 0;
    u64 harmonic_schedule_total = 0;
    u64 mixed_schedule_total = 0;

    u64 power_success = 0;
    u64 harmonic_success = 0;
    u64 mixed_success = 0;

    u64 quarter_success = 0;

    u64 max_power_size = 0;
    u64 max_harmonic_size = 0;
    u64 max_mixed_size = 0;

    u64 min_root = UINT64_MAX;
    u64 max_root = 0;

    u64 total_root = 0;

    /*
        Factor-blind validation:
        when the scheduled t is actually tested, check whether
        gcd(N,2t²-1) recovers a factor.
    */
    u64 mixed_factor_success = 0;
    u64 mixed_factor_attempts = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 s = c.s;

        const u64 root =
            smallest_p_root(p);

        if (root == 0) {
            ++no_p_root_cases;
            continue;
        }

        ++generic_cases;

        min_root =
            std::min(
                min_root,
                root
            );

        max_root =
            std::max(
                max_root,
                root
            );

        total_root += root;

        /*
            Deterministic scan bounds.
        */
        const u64 half =
            s / 2;

        const u64 quarter =
            s / 4;

        total_linear_half += half;
        total_linear_quarter += quarter;

        if (root <= quarter) {
            ++quarter_success;
        }

        /*
            Sparse schedules.
        */
        const std::vector<u64> power =
            build_power_schedule(s);

        const std::vector<u64> harmonic =
            build_harmonic_schedule(
                s,
                256
            );

        const std::vector<u64> mixed =
            build_mixed_schedule(s);

        power_schedule_total +=
            power.size();

        harmonic_schedule_total +=
            harmonic.size();

        mixed_schedule_total +=
            mixed.size();

        max_power_size =
            std::max<u64>(
                max_power_size,
                power.size()
            );

        max_harmonic_size =
            std::max<u64>(
                max_harmonic_size,
                harmonic.size()
            );

        max_mixed_size =
            std::max<u64>(
                max_mixed_size,
                mixed.size()
            );

        if (schedule_contains_root(
                power,
                root
        )) {
            ++power_success;
        }

        if (schedule_contains_root(
                harmonic,
                root
        )) {
            ++harmonic_success;
        }

        if (schedule_contains_root(
                mixed,
                root
        )) {
            ++mixed_success;
        }

        /*
            Actually execute the mixed schedule as a factor-blind
            search, rather than just checking root membership.
        */
        bool found = false;

        for (u64 t : mixed) {
            ++mixed_factor_attempts;

            if (direct_factor_test(c, t)) {
                ++mixed_factor_success;
                found = true;
                break;
            }
        }

        (void)found;
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "GENERIC_P_ROOT_CASES="
        << generic_cases
        << '\n';

    std::cout
        << "NO_P_ROOT_CASES="
        << no_p_root_cases
        << '\n';

    std::cout
        << "TOTAL_SMALLEST_ROOT="
        << total_root
        << '\n';

    std::cout
        << "MIN_ROOT="
        << min_root
        << '\n';

    std::cout
        << "MAX_ROOT="
        << max_root
        << '\n';

    std::cout
        << "AVERAGE_ROOT="
        << (
            generic_cases > 0
                ? static_cast<double>(
                    total_root
                  ) /
                  static_cast<double>(
                    generic_cases
                  )
                : 0.0
        )
        << '\n';

    std::cout
        << "TOTAL_LINEAR_HALF_TESTS="
        << total_linear_half
        << '\n';

    std::cout
        << "TOTAL_LINEAR_QUARTER_TESTS="
        << total_linear_quarter
        << '\n';

    std::cout
        << "QUARTER_SUCCESS="
        << quarter_success
        << '\n';

    std::cout
        << "POWER_SCHEDULE_SUCCESS="
        << power_success
        << '\n';

    std::cout
        << "HARMONIC_SCHEDULE_SUCCESS="
        << harmonic_success
        << '\n';

    std::cout
        << "MIXED_SCHEDULE_SUCCESS="
        << mixed_success
        << '\n';

    std::cout
        << "POWER_SCHEDULE_TOTAL="
        << power_schedule_total
        << '\n';

    std::cout
        << "HARMONIC_SCHEDULE_TOTAL="
        << harmonic_schedule_total
        << '\n';

    std::cout
        << "MIXED_SCHEDULE_TOTAL="
        << mixed_schedule_total
        << '\n';

    std::cout
        << "MAX_POWER_SCHEDULE_SIZE="
        << max_power_size
        << '\n';

    std::cout
        << "MAX_HARMONIC_SCHEDULE_SIZE="
        << max_harmonic_size
        << '\n';

    std::cout
        << "MAX_MIXED_SCHEDULE_SIZE="
        << max_mixed_size
        << '\n';

    std::cout
        << "MIXED_FACTOR_ATTEMPTS="
        << mixed_factor_attempts
        << '\n';

    std::cout
        << "MIXED_FACTOR_SUCCESS="
        << mixed_factor_success
        << '\n';

    const bool quarter_status =
        generic_cases == 0 ||
        quarter_success == generic_cases;

    std::cout
        << "QUARTER_BOUND_STATUS="
        << (
            quarter_status
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
