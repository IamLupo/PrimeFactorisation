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
            result.push_back({
                start,
                x - 1
            });

            inside = false;
        }
    }

    if (inside) {
        result.push_back({
            start,
            m
        });
    }

    return result;
}

/*
 * r = first digit position where j_r < q_r.
 */
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

/*
 * Verify that lower digits of j and q are identical.
 */
bool lower_digits_equal(
    u64 p,
    u64 j,
    u64 q,
    u64 r
) {
    if (r == 0) {
        return true;
    }

    const u64 pr = pow_u64(p, r);

    return (j % pr) == (q % pr);
}

/*
 * New length formula:

 *   L_r =
 *     (p^r - (q mod p^r)) * p^e - s0
 */
u64 predicted_length_from_r(
    u64 p,
    u64 e,
    u64 s0,
    u64 q,
    u64 r
) {
    const u64 pr = pow_u64(p, r);
    const u64 low = q % pr;
    const u64 pe = pow_u64(p, e);

    const u128 result =
        (static_cast<u128>(pr) -
         static_cast<u128>(low)) *
        static_cast<u128>(pe) -
        static_cast<u128>(s0);

    return static_cast<u64>(result);
}

void run_exhaustive(
    u64 max_m,
    u64 &intervals,
    u64 &lower_digits_pass,
    u64 &length_formula_pass,
    u64 &same_r_pairs,
    u64 &same_r_same_length
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    struct RKey {
        u64 p;
        u64 e;
        u64 s0;
        u64 q;
        u64 r;

        bool operator<(const RKey &other) const {
            if (p != other.p) return p < other.p;
            if (e != other.e) return e < other.e;
            if (s0 != other.s0) return s0 < other.s0;
            if (q != other.q) return q < other.q;
            return r < other.r;
        }
    };

    std::map<RKey, u64> first_length;

    for (const u64 p : primes) {
        for (u64 m = 2;
             m <= max_m;
             ++m) {

            Parameters par{};

            if (!derive_parameters(
                    p,
                    m,
                    par
                )) {
                continue;
            }

            const auto actual =
                collect_intervals(
                    p,
                    m
                );

            for (const auto &interval :
                 actual) {

                ++intervals;

                const u64 x =
                    interval.start;

                if (x < par.s0) {
                    continue;
                }

                const u64 delta =
                    x - par.s0;

                if (delta % par.pe != 0) {
                    continue;
                }

                const u64 j =
                    delta / par.pe;

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

                /*
                 * Core new identity.
                 */
                if (lower_digits_equal(
                        p,
                        j,
                        par.q,
                        r
                    )) {
                    ++lower_digits_pass;
                }

                const u64 predicted =
                    predicted_length_from_r(
                        p,
                        par.e,
                        par.s0,
                        par.q,
                        r
                    );

                const u64 actual_length =
                    interval.end -
                    interval.start +
                    1;

                if (predicted ==
                    actual_length) {
                    ++length_formula_pass;
                }

                /*
                 * Check that r alone determines
                 * the length for a fixed system.
                 */
                const RKey key{
                    p,
                    par.e,
                    par.s0,
                    par.q,
                    r
                };

                auto it =
                    first_length.find(key);

                if (it == first_length.end()) {
                    first_length.emplace(
                        key,
                        actual_length
                    );
                } else {
                    ++same_r_pairs;

                    if (it->second ==
                        actual_length) {
                        ++same_r_same_length;
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
            << " q=" << par.q
            << "\n";

        for (const auto &interval :
             intervals) {

            const u64 j =
                (interval.start -
                 par.s0) /
                par.pe;

            const u64 r =
                lowest_r(
                    par.p,
                    j,
                    par.q
                );

            const u64 qlow =
                r == 0
                ? 0
                : par.q % pow_u64(
                    par.p,
                    r
                );

            const u64 length =
                interval.end -
                interval.start +
                1;

            std::cout
                << "  j=" << j
                << " r=" << r
                << " q_low=" << qlow
                << " length=" << length
                << "\n";
        }
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 185\n";

    u64 intervals = 0;
    u64 lower_digits_pass = 0;
    u64 length_formula_pass = 0;
    u64 same_r_pairs = 0;
    u64 same_r_same_length = 0;

    run_exhaustive(
        5000,
        intervals,
        lower_digits_pass,
        length_formula_pass,
        same_r_pairs,
        same_r_same_length
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "intervals_tested="
        << intervals
        << "\n";

    std::cout
        << "lower_digits_equal_pass="
        << lower_digits_pass
        << "/" << intervals
        << "\n";

    std::cout
        << "r_only_length_formula_pass="
        << length_formula_pass
        << "/" << intervals
        << "\n";

    std::cout
        << "same_r_pairs="
        << same_r_pairs
        << "\n";

    std::cout
        << "same_r_same_length="
        << same_r_same_length
        << "/" << same_r_pairs
        << "\n";

    print_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "j mod p^r = q mod p^r\n";

    std::cout
        << "L_r = "
        << "(p^r - (q mod p^r))*p^e - s0\n";

    std::cout
        << "Therefore interval length depends only on r.\n";

    std::cout
        << "FINISHED EXPERIMENT 185\n";

    return 0;
}
