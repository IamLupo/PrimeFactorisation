#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Parameters {
    u64 p;
    u64 m;
    u64 e;
    u64 s0;
    u64 q;
    u64 pe;
};

struct Interval {
    u64 start;
    u64 end;
};

void start_experiment() {
    std::cout << "START EXPERIMENT 183\n";
}

void finish_experiment() {
    std::cout << "FINISHED EXPERIMENT 183\n";
}

u64 pow_u64(u64 p, u64 e) {
    u128 result = 1;

    for (u64 i = 0; i < e; ++i) {
        result *= p;
    }

    return static_cast<u64>(result);
}

bool is_hit(u64 p, u64 m, u64 x) {
    u64 a = x;
    u64 b = m;

    while (a > 0 || b > 0) {
        const u64 ad = a % p;
        const u64 bd = b % p;

        if (ad > bd) {
            return true;
        }

        a /= p;
        b /= p;
    }

    return false;
}

bool digitwise_leq(u64 p, u64 x, u64 q) {
    while (x > 0 || q > 0) {
        const u64 xd = x % p;
        const u64 qd = q % p;

        if (xd > qd) {
            return false;
        }

        x /= p;
        q /= p;
    }

    return true;
}

bool derive_parameters(
    u64 p,
    u64 m,
    Parameters &out
) {
    u64 x = m;
    u64 a = 0;

    while (x > 0 && x % p == p - 1) {
        ++a;
        x /= p;
    }

    if (x == 0) {
        return false;
    }

    const u64 b = x % p;

    x /= p;

    if (b >= p - 1) {
        return false;
    }

    u64 z = 0;

    while (x > 0 && x % p == 0) {
        ++z;
        x /= p;
    }

    const u64 e = a + 1 + z;
    const u64 pa = pow_u64(p, a);
    const u64 pe = pow_u64(p, e);

    const u128 numerator =
        static_cast<u128>(m) +
        1 -
        static_cast<u128>((b + 1) * pa);

    if (numerator == 0 ||
        numerator % pe != 0) {
        return false;
    }

    const u64 q =
        static_cast<u64>(numerator / pe);

    if (q == 0) {
        return false;
    }

    out = {
        p,
        m,
        e,
        (b + 1) * pa,
        q,
        pe
    };

    return true;
}

std::vector<Interval> collect_intervals(
    u64 p,
    u64 m
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 x = 1; x <= m; ++x) {
        const bool hit = is_hit(p, m, x);

        if (hit && !inside) {
            inside = true;
            start = x;
        }

        if (!hit && inside) {
            result.push_back({start, x - 1});
            inside = false;
        }
    }

    if (inside) {
        result.push_back({start, m});
    }

    return result;
}

u64 lowest_r(
    u64 p,
    u64 j,
    u64 q
) {
    u64 r = 0;

    while (j > 0 || q > 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd < qd) {
            return r;
        }

        j /= p;
        q /= p;
        ++r;
    }

    return UINT64_MAX;
}

u64 lower_part(
    u64 p,
    u64 j,
    u64 r
) {
    if (r == 0) {
        return 0;
    }

    return j % pow_u64(p, r);
}

u64 predicted_length(
    u64 p,
    u64 e,
    u64 s0,
    u64 j,
    u64 q
) {
    const u64 r = lowest_r(p, j, q);

    if (r == UINT64_MAX) {
        return UINT64_MAX;
    }

    const u64 pr = pow_u64(p, r);
    const u64 low = j % pr;
    const u64 pe = pow_u64(p, e);

    const u128 length =
        (static_cast<u128>(pr) -
         static_cast<u128>(low)) *
        static_cast<u128>(pe) -
        static_cast<u128>(s0);

    return static_cast<u64>(length);
}

void run_exhaustive(
    u64 max_m,
    u64 &intervals,
    u64 &formula_pass,
    u64 &lower_part_pass,
    u64 &same_lower_same_length,
    u64 &same_lower_pairs
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    /*
     * Key:
     *
     *   (p,e,s0,r,lower_part)
     *
     * should determine the interval length.
     */
    struct Key {
        u64 p;
        u64 e;
        u64 s0;
        u64 r;
        u64 low;

        bool operator<(const Key &other) const {
            if (p != other.p) return p < other.p;
            if (e != other.e) return e < other.e;
            if (s0 != other.s0) return s0 < other.s0;
            if (r != other.r) return r < other.r;
            return low < other.low;
        }
    };

    std::map<Key, u64> length_by_key;

    for (const u64 p : primes) {
        for (u64 m = 2; m <= max_m; ++m) {
            Parameters par{};

            if (!derive_parameters(
                    p,
                    m,
                    par
                )) {
                continue;
            }

            const auto actual =
                collect_intervals(p, m);

            for (const auto &interval : actual) {
                ++intervals;

                const u64 x = interval.start;

                if (x < par.s0) {
                    continue;
                }

                const u64 delta = x - par.s0;

                if (delta % par.pe != 0) {
                    continue;
                }

                const u64 j = delta / par.pe;

                if (!digitwise_leq(
                        p,
                        j,
                        par.q
                    )) {
                    continue;
                }

                const u64 r =
                    lowest_r(
                        p,
                        j,
                        par.q
                    );

                const u64 low =
                    lower_part(
                        p,
                        j,
                        r
                    );

                const u64 predicted =
                    predicted_length(
                        p,
                        par.e,
                        par.s0,
                        j,
                        par.q
                    );

                const u64 actual_length =
                    interval.end - interval.start + 1;

                if (predicted == actual_length) {
                    ++formula_pass;
                }

                /*
                 * Check the algebraic decomposition directly.
                 *
                 * L + s0 = (p^r - low) p^e.
                 */
                const u64 pe = par.pe;

                const u128 lhs =
                    static_cast<u128>(
                        actual_length
                    ) +
                    static_cast<u128>(
                        par.s0
                    );

                const u128 rhs =
                    (static_cast<u128>(
                        pow_u64(p, r)
                    ) -
                     static_cast<u128>(low)) *
                    static_cast<u128>(pe);

                if (lhs == rhs) {
                    ++lower_part_pass;
                }

                const Key key{
                    p,
                    par.e,
                    par.s0,
                    r,
                    low
                };

                auto it = length_by_key.find(key);

                if (it == length_by_key.end()) {
                    length_by_key.emplace(
                        key,
                        actual_length
                    );
                } else {
                    ++same_lower_pairs;

                    if (it->second == actual_length) {
                        ++same_lower_same_length;
                    }
                }
            }
        }
    }
}

void print_examples() {
    struct Example {
        u64 p;
        u64 m;
    };

    const std::vector<Example> examples = {
        {2, 20},
        {2, 40},
        {3, 41},
        {5, 194}
    };

    std::cout << "\nEXAMPLES\n";

    for (const auto &ex : examples) {
        Parameters par{};

        if (!derive_parameters(
                ex.p,
                ex.m,
                par
            )) {
            continue;
        }

        const auto intervals =
            collect_intervals(
                ex.p,
                ex.m
            );

        std::cout
            << "p=" << ex.p
            << " m=" << ex.m
            << "\n";

        for (const auto &interval : intervals) {
            const u64 x = interval.start;

            const u64 j =
                (x - par.s0) / par.pe;

            const u64 r =
                lowest_r(
                    par.p,
                    j,
                    par.q
                );

            const u64 low =
                lower_part(
                    par.p,
                    j,
                    r
                );

            const u64 length =
                interval.end - x + 1;

            std::cout
                << "  j=" << j
                << " r=" << r
                << " low=" << low
                << " length=" << length
                << "\n";
        }
    }
}

int main() {
    print_start();

    u64 intervals = 0;
    u64 formula_pass = 0;
    u64 lower_part_pass = 0;
    u64 same_lower_same_length = 0;
    u64 same_lower_pairs = 0;

    run_exhaustive(
        5000,
        intervals,
        formula_pass,
        lower_part_pass,
        same_lower_same_length,
        same_lower_pairs
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "intervals_tested="
        << intervals
        << "\n";

    std::cout
        << "formula_pass="
        << formula_pass
        << "/" << intervals
        << "\n";

    std::cout
        << "lower_part_identity_pass="
        << lower_part_pass
        << "/" << intervals
        << "\n";

    std::cout
        << "same_lower_same_length="
        << same_lower_same_length
        << "/" << same_lower_pairs
        << "\n";

    print_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "L_j + s0 = "
        << "(p^r - (j mod p^r)) * p^e\n";

    std::cout
        << "where r = lowest digit with j_r < q_r\n";

    std::cout
        << "Therefore the length depends only on "
        << "the lower r digits of j.\n";

    finish_experiment();

    return 0;
}