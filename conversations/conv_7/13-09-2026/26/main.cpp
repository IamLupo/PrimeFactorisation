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

void print_start() {
    std::cout << "START EXPERIMENT 181\n";
}

void print_finish() {
    std::cout << "FINISHED EXPERIMENT 181\n";
}

u64 pow_u64(u64 p, u64 e) {
    u128 r = 1;

    for (u64 i = 0; i < e; ++i) {
        r *= p;
    }

    return static_cast<u64>(r);
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

bool digitwise_leq(u64 p, u64 x, u64 m) {
    while (x > 0 || m > 0) {
        const u64 xd = x % p;
        const u64 md = m % p;

        if (xd > md) {
            return false;
        }

        x /= p;
        m /= p;
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

/*
 * Collect actual consecutive HIT intervals in ONE scan.
 */
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
 * Determine the lowest base-p digit position k for which:

 *   x_k < m_k

 * while all higher digits satisfy

 *   x_i <= m_i.

 * For an actual HIT start x, such a position exists.
 */
u64 lowest_eligible_position(
    u64 p,
    u64 m,
    u64 x
) {
    /*
     * We only need to know whether all higher digits
     * are valid. Process from the most significant side.
     *
     * First store digits.
     */
    std::vector<u64> xd;
    std::vector<u64> md;

    u64 a = x;
    u64 b = m;

    while (a > 0 || b > 0) {
        xd.push_back(a % p);
        md.push_back(b % p);

        a /= p;
        b /= p;
    }

    const std::size_t n =
        std::max(xd.size(), md.size());

    xd.resize(n, 0);
    md.resize(n, 0);

    bool higher_ok = true;

    /*
     * Scan from MSB to LSB and remember the lowest
     * position where x_k < m_k while higher digits
     * remain valid.
     */
    u64 answer = UINT64_MAX;

    for (std::size_t i = n; i-- > 0;) {
        const u64 xd_i = xd[i];
        const u64 md_i = md[i];

        if (xd_i > md_i) {
            higher_ok = false;
        }

        if (higher_ok && xd_i < md_i) {
            answer = static_cast<u64>(i);
        }
    }

    return answer;
}

/*
 * Proposed direct formula:

 *   nextMiss(x) - x
 *       = p^k - (x mod p^k)

 * where k is the lowest eligible position.
 */
u64 predicted_gap(
    u64 p,
    u64 m,
    u64 x
) {
    const u64 k =
        lowest_eligible_position(
            p,
            m,
            x
        );

    if (k == UINT64_MAX) {
        return UINT64_MAX;
    }

    const u64 pk = pow_u64(p, k);

    return pk - (x % pk);
}

u64 predicted_next_miss(
    u64 p,
    u64 m,
    u64 x
) {
    const u64 gap =
        predicted_gap(
            p,
            m,
            x
        );

    if (gap == UINT64_MAX) {
        return UINT64_MAX;
    }

    return x + gap;
}

/*
 * Exhaustive small validation.

 * Important:
 *   We do NOT brute-force nextMiss for every x.

 * Instead:
 *   - scan each m once
 *   - obtain actual HIT intervals
 *   - test formula only at actual interval starts.
 */
void run_small_exhaustive(
    u64 max_m,
    u64 &parameter_cases,
    u64 &intervals_tested,
    u64 &next_miss_pass,
    u64 &length_pass,
    u64 &boundary_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

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

            ++parameter_cases;

            const auto intervals =
                collect_intervals(
                    p,
                    m
                );

            for (const auto &interval : intervals) {
                ++intervals_tested;

                const u64 x =
                    interval.start;

                const u64 actual_next =
                    interval.end + 1;

                const u64 predicted_next =
                    predicted_next_miss(
                        p,
                        m,
                        x
                    );

                if (predicted_next ==
                    actual_next) {
                    ++next_miss_pass;
                }

                if (predicted_next !=
                    UINT64_MAX) {

                    const u64 predicted_length =
                        predicted_next - x;

                    const u64 actual_length =
                        interval.end - x + 1;

                    if (predicted_length ==
                        actual_length) {
                        ++length_pass;
                    }

                    bool ok = true;

                    if (!is_hit(
                            p,
                            m,
                            x
                        )) {
                        ok = false;
                    }

                    if (predicted_next <= m &&
                        is_hit(
                            p,
                            m,
                            predicted_next
                        )) {
                        ok = false;
                    }

                    if (x > 1 &&
                        is_hit(
                            p,
                            m,
                            x - 1
                        )) {
                        ok = false;
                    }

                    if (ok) {
                        ++boundary_pass;
                    }
                }
            }
        }
    }
}

/*
 * Targeted structural examples.
 */
void run_examples() {
    struct Example {
        u64 p;
        u64 m;
    };

    const std::vector<Example> examples = {
        {2, 5},
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
            const u64 x =
                interval.start;

            const u64 gap =
                predicted_gap(
                    ex.p,
                    ex.m,
                    x
                );

            const u64 predicted_end =
                x + gap - 1;

            const u64 k =
                lowest_eligible_position(
                    ex.p,
                    ex.m,
                    x
                );

            std::cout
                << "  start=" << x
                << " actual_end="
                << interval.end
                << " k=" << k
                << " gap=" << gap
                << " predicted_end="
                << predicted_end
                << " length="
                << (interval.end - x + 1)
                << "\n";
        }
    }
}

int main() {
    print_start();

    u64 parameter_cases = 0;
    u64 intervals_tested = 0;
    u64 next_miss_pass = 0;
    u64 length_pass = 0;
    u64 boundary_pass = 0;

    run_small_exhaustive(
        5000,
        parameter_cases,
        intervals_tested,
        next_miss_pass,
        length_pass,
        boundary_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << parameter_cases
        << "\n";

    std::cout
        << "intervals_tested="
        << intervals_tested
        << "\n";

    std::cout
        << "next_miss_pass="
        << next_miss_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "length_pass="
        << length_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "boundary_pass="
        << boundary_pass
        << "/" << intervals_tested
        << "\n";

    run_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "nextMiss(x)-x = "
        << "p^k - (x mod p^k)\n";

    std::cout
        << "k = lowest eligible base-p digit position\n";

    std::cout
        << "I(x) = "
        << "[x, x + p^k - (x mod p^k) - 1]\n";

    print_finish();

    return 0;
}