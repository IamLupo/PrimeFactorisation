#include <algorithm>
#include <cstdint>
#include <iostream>
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
 * New endpoint formula:

 *   E =
 *     ( floor(j / p^r) + 1 ) * p^(e+r) - 1
 */
u64 predicted_endpoint(
    u64 p,
    u64 e,
    u64 j,
    u64 q
) {
    const u64 r =
        lowest_r(p, j, q);

    if (r == UINT64_MAX) {
        return UINT64_MAX;
    }

    const u64 pr =
        pow_u64(p, r);

    const u64 per =
        pow_u64(p, e + r);

    const u64 high =
        j / pr;

    const u128 endpoint =
        static_cast<u128>(high + 1) *
        static_cast<u128>(per) -
        1;

    return static_cast<u64>(endpoint);
}

void run_exhaustive(
    u64 max_m,
    u64 &intervals_tested,
    u64 &endpoint_pass,
    u64 &divisibility_pass,
    u64 &length_pass,
    u64 &start_end_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

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

            const auto intervals =
                collect_intervals(p, m);

            for (const auto &interval :
                 intervals) {

                ++intervals_tested;

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

                if (r == UINT64_MAX) {
                    continue;
                }

                const u64 endpoint =
                    predicted_endpoint(
                        p,
                        par.e,
                        j,
                        par.q
                    );

                if (endpoint ==
                    interval.end) {
                    ++endpoint_pass;
                }

                /*
                 * Endpoint + 1 must be exactly divisible
                 * by p^(e+r).
                 */
                const u64 modulus =
                    pow_u64(
                        p,
                        par.e + r
                    );

                if ((interval.end + 1) %
                    modulus == 0) {
                    ++divisibility_pass;
                }

                /*
                 * Check the endpoint formula against
                 * the length formula from Exp. 183.
                 */
                const u64 low =
                    j % pow_u64(p, r);

                const u64 length =
                    static_cast<u64>(
                        (
                            static_cast<u128>(
                                pow_u64(p, r)
                            ) -
                            static_cast<u128>(low)
                        ) *
                        static_cast<u128>(par.pe) -
                        static_cast<u128>(par.s0)
                    );

                const u64 reconstructed_end =
                    x + length - 1;

                if (reconstructed_end ==
                    interval.end) {
                    ++length_pass;
                }

                /*
                 * Complete direct check:
                 *
                 * start = s0 + j*p^e
                 * end   = formula
                 */
                const u64 reconstructed_start =
                    static_cast<u64>(
                        static_cast<u128>(
                            par.s0
                        ) +
                        static_cast<u128>(j) *
                            static_cast<u128>(
                                par.pe
                            )
                    );

                if (reconstructed_start ==
                        interval.start &&
                    endpoint ==
                        interval.end) {
                    ++start_end_pass;
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

        for (const auto &interval :
             intervals) {

            const u64 j =
                (interval.start - par.s0) /
                par.pe;

            const u64 r =
                lowest_r(
                    par.p,
                    j,
                    par.q
                );

            const u64 endpoint =
                predicted_endpoint(
                    par.p,
                    par.e,
                    j,
                    par.q
                );

            const u64 modulus =
                pow_u64(
                    par.p,
                    par.e + r
                );

            std::cout
                << "  j=" << j
                << " r=" << r
                << " start="
                << interval.start
                << " end="
                << interval.end
                << " end+1="
                << interval.end + 1
                << " modulus="
                << modulus
                << "\n";
        }
    }
}

int main() {
    std::cout << "START EXPERIMENT 184\n";

    u64 intervals_tested = 0;
    u64 endpoint_pass = 0;
    u64 divisibility_pass = 0;
    u64 length_pass = 0;
    u64 start_end_pass = 0;

    run_exhaustive(
        5000,
        intervals_tested,
        endpoint_pass,
        divisibility_pass,
        length_pass,
        start_end_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "intervals_tested="
        << intervals_tested
        << "\n";

    std::cout
        << "endpoint_formula_pass="
        << endpoint_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "endpoint_divisibility_pass="
        << divisibility_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "length_equivalence_pass="
        << length_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "complete_start_end_pass="
        << start_end_pass
        << "/" << intervals_tested
        << "\n";

    print_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "r = lowest digit with j_r < q_r\n";

    std::cout
        << "E_j = "
        << "(floor(j / p^r) + 1) * p^(e+r) - 1\n";

    std::cout
        << "equivalently:\n";

    std::cout
        << "E_j + 1 is divisible by p^(e+r)\n";

    std::cout
        << "FINISHED EXPERIMENT 184\n";

    return 0;
}
