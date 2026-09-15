#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;

struct PrimePair {
    u64 p;
    u64 q;
    u64 n;
    u64 s;
};

struct RootInfo {
    bool exists = false;
    u64 k = 0;
    u64 t = 0;
};

static u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        u64 r = a % b;
        a = b;
        b = r;
    }
    return a;
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

        if (d > 50000) {
            break;
        }
    }

    return true;
}

static u64 random_prime(std::mt19937_64& rng, u64 lo, u64 hi) {
    std::uniform_int_distribution<u64> dist(lo, hi);

    while (true) {
        u64 x = dist(rng) | 1ULL;

        if (x < lo || x > hi) {
            continue;
        }

        if (is_prime_u64(x)) {
            return x;
        }
    }
}

static u64 integer_sqrt_u64(u64 n) {
    u64 lo = 0;
    u64 hi = 1000000000ULL;

    while (lo <= hi) {
        const u64 mid = lo + (hi - lo) / 2;
        const u128 sq = static_cast<u128>(mid) * mid;

        if (sq <= n) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return hi;
}

static PrimePair generate_case(std::mt19937_64& rng) {
    const u64 p = random_prime(rng, 10000, 1000000);

    u64 q = random_prime(rng, 10000, 1000000);

    while (q == p) {
        q = random_prime(rng, 10000, 1000000);
    }

    const u64 a = std::min(p, q);
    const u64 b = std::max(p, q);

    const u64 n = a * b;
    const u64 s = integer_sqrt_u64(n);

    return {a, b, n, s};
}

static bool inverse_mod(
    u64 a,
    u64 m,
    u64& inverse
) {
    if (m == 0) {
        return false;
    }

    a %= m;

    if (a == 0 || gcd_u64(a, m) != 1) {
        return false;
    }

    i64 old_r = static_cast<i64>(m);
    i64 r = static_cast<i64>(a);

    i64 old_s = 0;
    i64 s = 1;

    while (r != 0) {
        const i64 q = old_r / r;

        const i64 next_r = old_r - q * r;
        old_r = r;
        r = next_r;

        const i64 next_s = old_s - q * s;
        old_s = s;
        s = next_s;
    }

    if (old_r != 1) {
        return false;
    }

    i64 value = old_s % static_cast<i64>(m);

    if (value < 0) {
        value += static_cast<i64>(m);
    }

    inverse = static_cast<u64>(value);

    return true;
}

static RootInfo best_root(
    u64 r,
    u64 k_max
) {
    RootInfo best;

    for (u64 k = 2; k <= k_max; ++k) {
        u64 inv = 0;

        if (!inverse_mod(k, r, inv)) {
            continue;
        }

        const u64 t = std::min(inv, r - inv);

        if (!best.exists ||
            t < best.t ||
            (t == best.t && k < best.k)) {
            best.exists = true;
            best.k = k;
            best.t = t;
        }
    }

    return best;
}

static bool verify_root(
    u64 r,
    const RootInfo& root
) {
    if (!root.exists) {
        return false;
    }

    const u128 kt =
        static_cast<u128>(root.k) * root.t;

    const u64 value_minus =
        static_cast<u64>(kt - 1);

    const u64 value_plus =
        static_cast<u64>(kt + 1);

    return
        (value_minus % r == 0) ||
        (value_plus % r == 0);
}

static void run_kmax(
    const std::vector<PrimePair>& cases,
    u64 k_max
) {
    u64 valid_p = 0;
    u64 valid_q = 0;

    u64 verification_failures = 0;
    u64 lower_bound_failures_p = 0;
    u64 lower_bound_failures_q = 0;

    u64 p_first = 0;
    u64 q_first = 0;
    u64 tied = 0;

    u64 total_p_root = 0;
    u64 total_q_root = 0;

    u64 min_p_root = std::numeric_limits<u64>::max();
    u64 max_p_root = 0;

    u64 min_q_root = std::numeric_limits<u64>::max();
    u64 max_q_root = 0;

    u64 total_p_gap = 0;
    u64 total_q_gap = 0;

    u64 min_p_gap = std::numeric_limits<u64>::max();
    u64 max_p_gap = 0;

    u64 min_q_gap = std::numeric_limits<u64>::max();
    u64 max_q_gap = 0;

    for (const PrimePair& c : cases) {
        const RootInfo p_root = best_root(c.p, k_max);
        const RootInfo q_root = best_root(c.q, k_max);

        if (p_root.exists) {
            ++valid_p;

            if (!verify_root(c.p, p_root)) {
                ++verification_failures;
            }

            const u64 lower =
                (c.p - 1 + k_max - 1) / k_max;

            if (p_root.t < lower) {
                ++lower_bound_failures_p;
            }

            const u64 gap =
                p_root.t * k_max - (c.p - 1);

            ++total_p_root;

            min_p_root = std::min(min_p_root, p_root.t);
            max_p_root = std::max(max_p_root, p_root.t);

            total_p_gap += gap;
            min_p_gap = std::min(min_p_gap, gap);
            max_p_gap = std::max(max_p_gap, gap);
        }

        if (q_root.exists) {
            ++valid_q;

            if (!verify_root(c.q, q_root)) {
                ++verification_failures;
            }

            const u64 lower =
                (c.q - 1 + k_max - 1) / k_max;

            if (q_root.t < lower) {
                ++lower_bound_failures_q;
            }

            const u64 gap =
                q_root.t * k_max - (c.q - 1);

            ++total_q_root;

            min_q_root = std::min(min_q_root, q_root.t);
            max_q_root = std::max(max_q_root, q_root.t);

            total_q_gap += gap;
            min_q_gap = std::min(min_q_gap, gap);
            max_q_gap = std::max(max_q_gap, gap);
        }

        if (p_root.exists && q_root.exists) {
            if (p_root.t < q_root.t) {
                ++p_first;
            } else if (q_root.t < p_root.t) {
                ++q_first;
            } else {
                ++tied;
            }
        }
    }

    std::cout << "K_MAX=" << k_max << '\n';

    std::cout << "  P_VALID=" << valid_p << '\n';
    std::cout << "  Q_VALID=" << valid_q << '\n';

    std::cout << "  P_FIRST=" << p_first << '\n';
    std::cout << "  Q_FIRST=" << q_first << '\n';
    std::cout << "  TIED=" << tied << '\n';

    std::cout << "  P_MIN_ROOT=" << min_p_root << '\n';
    std::cout << "  P_MAX_ROOT=" << max_p_root << '\n';

    std::cout << "  Q_MIN_ROOT=" << min_q_root << '\n';
    std::cout << "  Q_MAX_ROOT=" << max_q_root << '\n';

    std::cout
        << "  P_AVERAGE_ROOT="
        << static_cast<double>(total_p_root == 0 ? 0 : total_p_root) /
           static_cast<double>(valid_p)
        << '\n';

    std::cout
        << "  Q_AVERAGE_ROOT="
        << static_cast<double>(total_q_root == 0 ? 0 : total_q_root) /
           static_cast<double>(valid_q)
        << '\n';

    std::cout
        << "  P_LOWER_BOUND_FAILURES="
        << lower_bound_failures_p
        << '\n';

    std::cout
        << "  Q_LOWER_BOUND_FAILURES="
        << lower_bound_failures_q
        << '\n';

    std::cout
        << "  P_GAP_MIN="
        << min_p_gap
        << '\n';

    std::cout
        << "  P_GAP_MAX="
        << max_p_gap
        << '\n';

    std::cout
        << "  P_GAP_AVERAGE="
        << static_cast<double>(total_p_gap) /
           static_cast<double>(valid_p)
        << '\n';

    std::cout
        << "  Q_GAP_MIN="
        << min_q_gap
        << '\n';

    std::cout
        << "  Q_GAP_MAX="
        << max_q_gap
        << '\n';

    std::cout
        << "  Q_GAP_AVERAGE="
        << static_cast<double>(total_q_gap) /
           static_cast<double>(valid_q)
        << '\n';

    std::cout
        << "  ROOT_VERIFICATION_FAILURES="
        << verification_failures
        << '\n';

    std::cout << '\n';
}

int main() {
    constexpr u64 EXPERIMENT = 395;
    constexpr u64 CASES = 3000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(395123456789ULL);

    std::vector<PrimePair> samples;
    samples.reserve(CASES);

    for (u64 i = 0; i < CASES; ++i) {
        samples.push_back(generate_case(rng));
    }

    std::cout << "CASES=" << CASES << '\n';
    std::cout << "K_MIN=2\n";
    std::cout << "K_MAX=20\n";
    std::cout << "PRIME_MIN=10000\n";
    std::cout << "PRIME_MAX=1000000\n\n";

    for (u64 k_max = 2; k_max <= 20; ++k_max) {
        run_kmax(samples, k_max);
    }

    std::cout
        << "STATUS=PASS_IF_ROOT_VERIFICATION_FAILURES_AND_LOWER_BOUND_FAILURES_ARE_ZERO"
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}