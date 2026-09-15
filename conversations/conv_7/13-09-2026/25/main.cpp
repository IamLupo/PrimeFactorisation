#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
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
    std::cout << "START EXPERIMENT 179\n";
}

void finish_experiment() {
    std::cout << "FINISHED EXPERIMENT 179\n";
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

bool is_miss(u64 p, u64 m, u64 x) {
    return !is_hit(p, m, x);
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

    if (numerator == 0 || numerator % pe != 0) {
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

u64 lowest_available_digit(
    u64 p,
    u64 j,
    u64 q
) {
    u64 k = 0;

    while (j > 0 || q > 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd < qd) {
            return k;
        }

        j /= p;
        q /= p;
        ++k;
    }

    return UINT64_MAX;
}

u64 next_valid_index(
    u64 p,
    u64 j,
    u64 q
) {
    const u64 k =
        lowest_available_digit(p, j, q);

    if (k == UINT64_MAX) {
        return q;
    }

    const u64 pk = pow_u64(p, k);

    return j + pk - (j % pk);
}

/*
 * Exact candidate formula for the next MISS.

 * We search for the smallest digit position k such that:
 *
 *   x_k < m_k
 *
 * and all higher digits of x are already <= the
 * corresponding digits of m.
 *
 * Then:
 *
 *   y_i = x_i, i > k
 *   y_k = x_k + 1
 *   y_i = 0, i < k
 *
 * This y is the smallest digitwise-valid number
 * larger than x for that k.
 */
u64 next_miss_formula(
    u64 p,
    u64 m,
    u64 x
) {
    /*
     * Store digits LSF-first.
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

    /*
     * higher_ok[k] means all digits i>k
     * satisfy x_i <= m_i.
     */
    std::vector<bool> higher_ok(n, true);

    bool ok = true;

    for (std::size_t i = n; i-- > 0;) {
        if (i + 1 < n) {
            higher_ok[i] = ok;
        }

        if (xd[i] > md[i]) {
            ok = false;
        }
    }

    bool found = false;
    u64 best = UINT64_MAX;

    u64 place = 1;

    for (std::size_t k = 0; k < n; ++k) {
        if (xd[k] >= md[k]) {
            place *= p;
            continue;
        }

        if (!higher_ok[k]) {
            place *= p;
            continue;
        }

        /*
         * Keep higher digits unchanged.
         */
        u128 higher =
            static_cast<u128>(x) /
            static_cast<u128>(place * p);

        u128 candidate =
            higher *
            static_cast<u128>(place * p);

        /*
         * Increment digit k by one.
         */
        candidate +=
            static_cast<u128>(xd[k] + 1) *
            static_cast<u128>(place);

        /*
         * Lower digits are zero automatically.
         */
        if (candidate > x &&
            candidate <= UINT64_MAX) {
            const u64 y =
                static_cast<u64>(candidate);

            if (!found || y < best) {
                best = y;
                found = true;
            }
        }

        place *= p;
    }

    if (!found) {
        return UINT64_MAX;
    }

    return best;
}

/*
 * Brute-force next MISS. Used only for small values.
 */
u64 next_miss_bruteforce(
    u64 p,
    u64 m,
    u64 x
) {
    for (u64 y = x + 1; y <= m + 1; ++y) {
        if (is_miss(p, m, y)) {
            return y;
        }
    }

    return UINT64_MAX;
}

void run_small_exhaustive(
    u64 &tested,
    u64 &passed,
    u64 &boundary_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (u64 p : primes) {
        for (u64 m = 2; m <= 5000; ++m) {
            Parameters par{};

            if (!derive_parameters(p, m, par)) {
                continue;
            }

            /*
             * Only test actual HIT starts.
             */
            for (u64 j = 0; j < par.q; ++j) {
                /*
                 * The start-index characterization.
                 */
                if (!digitwise_leq(p, j, par.q)) {
                    continue;
                }

                const u64 start =
                    static_cast<u64>(
                        static_cast<u128>(par.s0) +
                        static_cast<u128>(j) * par.pe
                    );

                if (!is_hit(p, m, start)) {
                    continue;
                }

                const u64 actual =
                    next_miss_bruteforce(
                        p,
                        m,
                        start
                    );

                const u64 predicted =
                    next_miss_formula(
                        p,
                        m,
                        start
                    );

                ++tested;

                if (actual == predicted) {
                    ++passed;
                }

                /*
                 * The interval itself is
                 *
                 *   [start, next_miss - 1].
                 */
                if (predicted != UINT64_MAX) {
                    const u64 end =
                        predicted - 1;

                    bool ok =
                        is_hit(p, m, start) &&
                        is_hit(p, m, end);

                    if (start > 1 &&
                        is_hit(
                            p,
                            m,
                            start - 1
                        )) {
                        ok = false;
                    }

                    if (end < m &&
                        is_hit(
                            p,
                            m,
                            end + 1
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

void run_large_random(
    u64 trials,
    u64 &formula_pass,
    u64 &boundary_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    std::mt19937_64 rng(
        0x17920260913ULL
    );

    for (u64 trial = 0;
         trial < trials;
         ++trial) {

        const u64 p =
            primes[trial % primes.size()];

        std::uniform_int_distribution<u64> dist(
            1000000000000ULL,
            1000000000000000000ULL
        );

        const u64 m = dist(rng);

        Parameters par{};

        if (!derive_parameters(p, m, par)) {
            --trial;
            continue;
        }

        /*
         * Generate a random valid start index
         * digitwise <= q.
         */
        std::vector<u64> qdigits;

        u64 q = par.q;

        while (q > 0) {
            qdigits.push_back(q % p);
            q /= p;
        }

        if (qdigits.empty()) {
            qdigits.push_back(0);
        }

        u64 j = 0;

        for (;;) {
            u128 value = 0;
            u128 place = 1;

            for (u64 digit : qdigits) {
                std::uniform_int_distribution<u64> pick(
                    0,
                    digit
                );

                value +=
                    static_cast<u128>(
                        pick(rng)
                    ) *
                    place;

                place *= p;
            }

            j = static_cast<u64>(value);

            if (j < par.q) {
                break;
            }
        }

        const u64 start =
            static_cast<u64>(
                static_cast<u128>(par.s0) +
                static_cast<u128>(j) * par.pe
            );

        const u64 predicted =
            next_miss_formula(
                p,
                m,
                start
            );

        if (predicted == UINT64_MAX) {
            continue;
        }

        ++formula_pass;

        const u64 end =
            predicted - 1;

        bool ok =
            start <= m &&
            end <= m &&
            start <= end &&
            is_hit(p, m, start) &&
            is_miss(p, m, predicted);

        if (start > 1 &&
            is_hit(
                p,
                m,
                start - 1
            )) {
            ok = false;
        }

        /*
         * We cannot scan a huge interval, so test
         * the candidate construction directly:
         *
         * predicted must be digitwise <= m.
         */
        if (!digitwise_leq(
                p,
                predicted,
                m
            )) {
            ok = false;
        }

        if (ok) {
            ++boundary_pass;
        }
    }
}

void show_examples() {
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

    std::cout
        << "\nEXAMPLES\n";

    for (const auto &ex : examples) {
        Parameters par{};

        if (!derive_parameters(
                ex.p,
                ex.m,
                par
            )) {
            continue;
        }

        std::cout
            << "p=" << ex.p
            << " m=" << ex.m
            << "\n";

        for (u64 j = 0; j < par.q; ++j) {
            if (!digitwise_leq(
                    par.p,
                    j,
                    par.q
                )) {
                continue;
            }

            const u64 start =
                static_cast<u64>(
                    static_cast<u128>(par.s0) +
                    static_cast<u128>(j) *
                        par.pe
                );

            const u64 miss =
                next_miss_formula(
                    par.p,
                    par.m,
                    start
                );

            const u64 end =
                miss - 1;

            std::cout
                << "  j=" << j
                << " start=" << start
                << " next_miss=" << miss
                << " end=" << end
                << " length="
                << (end - start + 1)
                << "\n";
        }
    }
}

int main() {
    start_experiment();

    u64 tested = 0;
    u64 passed = 0;
    u64 boundary_pass = 0;

    run_small_exhaustive(
        tested,
        passed,
        boundary_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "next_miss_cases="
        << tested
        << "\n";

    std::cout
        << "next_miss_formula_pass="
        << passed
        << "/" << tested
        << "\n";

    std::cout
        << "interval_boundary_pass="
        << boundary_pass
        << "/" << tested
        << "\n";

    u64 large_formula_pass = 0;
    u64 large_boundary_pass = 0;

    run_large_random(
        200000,
        large_formula_pass,
        large_boundary_pass
    );

    std::cout
        << "\nLARGE RANDOM\n";

    std::cout
        << "next_miss_formula_pass="
        << large_formula_pass
        << "/200000\n";

    std::cout
        << "boundary_pass="
        << large_boundary_pass
        << "/200000\n";

    show_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "next_miss(x) = smallest y>x with y <=_p m\n";

    std::cout
        << "I(x) = [x, next_miss(x)-1]\n";

    finish_experiment();

    return 0;
}
