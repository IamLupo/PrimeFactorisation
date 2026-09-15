#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Interval {
    u64 start;
    u64 end;
};

struct Parameters {
    u64 p;
    u64 m;
    u64 a;
    u64 b;
    u64 z;
    u64 e;
    u64 s0;
    u64 q;
    u64 pe;
};

void print_start() {
    std::cout << "START EXPERIMENT 176\n";
}

void print_finish() {
    std::cout << "FINISHED EXPERIMENT 176\n";
}

u64 pow_u64(u64 base, u64 exp) {
    u128 result = 1;

    for (u64 i = 0; i < exp; ++i) {
        result *= base;
    }

    return static_cast<u64>(result);
}

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d * d <= n; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

std::vector<u64> small_primes() {
    return {2, 3, 5, 7, 11};
}

/*
 * Extract the experimentally observed structure

 *   a = number of initial base-p digits equal to p-1
 *   b = first digit < p-1
 *   z = number of zero digits immediately following b
 *   e = a + 1 + z
 *   s0 = (b+1) * p^a
 *   m = q*p^e + s0 - 1
 */
bool derive_parameters(u64 p, u64 m, Parameters &out) {
    if (p < 2 || m == 0) {
        return false;
    }

    u64 x = m;
    u64 a = 0;

    while (x > 0 && (x % p) == p - 1) {
        ++a;
        x /= p;
    }

    /*
     * If every digit was p-1, there is no first b < p-1.
     * Exclude this degenerate case.
     */
    if (x == 0) {
        return false;
    }

    const u64 b = x % p;
    x /= p;

    if (b >= p - 1) {
        return false;
    }

    u64 z = 0;

    while (x > 0 && (x % p) == 0) {
        ++z;
        x /= p;
    }

    const u64 e = a + 1 + z;
    const u64 pa = pow_u64(p, a);
    const u64 pe = pow_u64(p, e);

    const u128 numerator =
        static_cast<u128>(m) + 1 - static_cast<u128>((b + 1) * pa);

    if (numerator == 0 || numerator % pe != 0) {
        return false;
    }

    const u64 q = static_cast<u64>(numerator / pe);

    if (q == 0) {
        return false;
    }

    out = {p, m, a, b, z, e, (b + 1) * pa, q, pe};
    return true;
}

/*
 * Lucas/Kummer HIT predicate:

 *   MISS iff every base-p digit of t is <=
 *   the corresponding digit of m.

 * Therefore:
 *   HIT iff at least one digit of t exceeds m.
 */
bool is_hit(u64 p, u64 m, u64 t) {
    u64 a = t;
    u64 b = m;

    while (a > 0 || b > 0) {
        const u64 td = a % p;
        const u64 md = b % p;

        if (td > md) {
            return true;
        }

        a /= p;
        b /= p;
    }

    return false;
}

bool digitwise_leq(u64 p, u64 j, u64 q) {
    while (j > 0 || q > 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd > qd) {
            return false;
        }

        j /= p;
        q /= p;
    }

    return true;
}

/*
 * Return the lowest digit position k such that

 *     j_k < q_k

 * assuming j < q and j <=_p q.
 */
u64 lowest_available_digit(u64 p, u64 j, u64 q) {
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

/*
 * Exact next-valid-index formula:

 *   j+ = j + p^k - (j mod p^k)

 * where k is the lowest digit with j_k < q_k.
 */
u64 next_valid_index(u64 p, u64 j, u64 q) {
    const u64 k = lowest_available_digit(p, j, q);

    if (k == UINT64_MAX) {
        return q;
    }

    const u64 pk = pow_u64(p, k);
    const u64 remainder = j % pk;

    return j + (pk - remainder);
}

/*
 * Collect all actual HIT intervals by direct scanning.
 */
std::vector<Interval> collect_actual_intervals(u64 p, u64 m) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit = is_hit(p, m, t);

        if (hit && !inside) {
            inside = true;
            start = t;
        }

        if (!hit && inside) {
            inside = false;
            result.push_back({start, t - 1});
        }
    }

    if (inside) {
        result.push_back({start, m});
    }

    return result;
}

/*
 * Construct all predicted intervals from valid start indices.
 *
 * Valid starts:
 *
 *   j < q
 *   j <=_p q
 *
 * and
 *
 *   start = s0 + j*p^e
 *   end   = s0 + j+*p^e - 1
 */
std::vector<Interval> collect_predicted_intervals(const Parameters &par) {
    std::vector<Interval> result;

    for (u64 j = 0; j < par.q; ++j) {
        if (!digitwise_leq(par.p, j, par.q)) {
            continue;
        }

        const u64 jp = next_valid_index(par.p, j, par.q);

        const u128 start128 =
            static_cast<u128>(par.s0) +
            static_cast<u128>(j) * par.pe;

        const u128 end128 =
            static_cast<u128>(par.s0) +
            static_cast<u128>(jp) * par.pe -
            1;

        result.push_back({
            static_cast<u64>(start128),
            static_cast<u64>(end128)
        });
    }

    return result;
}

/*
 * Direct formula for the gap.
 */
u64 predicted_gap(u64 p, u64 j, u64 q, u64 &k_out) {
    const u64 k = lowest_available_digit(p, j, q);
    k_out = k;

    if (k == UINT64_MAX) {
        return 0;
    }

    const u64 pk = pow_u64(p, k);

    return pk - (j % pk);
}

bool compare_intervals(
    const std::vector<Interval> &actual,
    const std::vector<Interval> &predicted
) {
    if (actual.size() != predicted.size()) {
        return false;
    }

    for (std::size_t i = 0; i < actual.size(); ++i) {
        if (actual[i].start != predicted[i].start ||
            actual[i].end != predicted[i].end) {
            return false;
        }
    }

    return true;
}

/*
 * Exhaustive small test.

 * For every small m:
 *   1. derive (e,s0,q)
 *   2. scan every t directly
 *   3. generate intervals from the digit formula
 *   4. verify gap formula for every interval
 */
void run_small_exhaustive(
    u64 &parameter_cases,
    u64 &intervals_tested,
    u64 &interval_formula_pass,
    u64 &gap_formula_pass
) {
    const auto primes = small_primes();

    for (u64 p : primes) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 2; m <= 2000; ++m) {
            Parameters par{};

            if (!derive_parameters(p, m, par)) {
                continue;
            }

            ++parameter_cases;

            const auto actual = collect_actual_intervals(p, m);
            const auto predicted = collect_predicted_intervals(par);

            if (compare_intervals(actual, predicted)) {
                ++interval_formula_pass;
            }

            for (u64 j = 0; j < par.q; ++j) {
                if (!digitwise_leq(p, j, par.q)) {
                    continue;
                }

                const u64 jp = next_valid_index(p, j, par.q);

                u64 k = 0;
                const u64 gap = predicted_gap(p, j, par.q, k);

                if (jp - j == gap) {
                    ++gap_formula_pass;
                }

                ++intervals_tested;
            }
        }
    }
}

/*
 * Generate a random digitwise-valid j.

 * Every digit of j is independently chosen from [0,q_i].
 * The resulting number automatically satisfies j <=_p q.
 */
u64 random_valid_index(
    u64 p,
    u64 q,
    std::mt19937_64 &rng
) {
    std::vector<u64> qdigits;

    u64 x = q;

    while (x > 0) {
        qdigits.push_back(x % p);
        x /= p;
    }

    if (qdigits.empty()) {
        qdigits.push_back(0);
    }

    for (;;) {
        u128 value = 0;
        u128 place = 1;

        for (u64 digit : qdigits) {
            std::uniform_int_distribution<u64> dist(0, digit);
            const u64 chosen = dist(rng);

            value += static_cast<u128>(chosen) * place;
            place *= p;
        }

        const u64 j = static_cast<u64>(value);

        if (j < q) {
            return j;
        }
    }
}

/*
 * Generate a large random m in a safe 64-bit range.
 */
u64 random_large_m(std::mt19937_64 &rng) {
    std::uniform_int_distribution<u64> dist(
        1000000000000ULL,
        1000000000000000000ULL
    );

    return dist(rng);
}

void run_large_random(
    u64 trials,
    u64 &parameter_pass,
    u64 &gap_pass,
    u64 &boundary_pass
) {
    const auto primes = small_primes();

    std::mt19937_64 rng(0x17620260913ULL);

    for (u64 trial = 0; trial < trials; ++trial) {
        const u64 p = primes[trial % primes.size()];
        const u64 m = random_large_m(rng);

        Parameters par{};

        if (!derive_parameters(p, m, par)) {
            --trial;
            continue;
        }

        ++parameter_pass;

        const u64 j = random_valid_index(p, par.q, rng);

        u64 k = 0;
        const u64 gap = predicted_gap(p, j, par.q, k);
        const u64 jp = next_valid_index(p, j, par.q);

        if (jp - j != gap) {
            continue;
        }

        ++gap_pass;

        const u128 start128 =
            static_cast<u128>(par.s0) +
            static_cast<u128>(j) * par.pe;

        const u128 end128 =
            static_cast<u128>(par.s0) +
            static_cast<u128>(jp) * par.pe -
            1;

        const u64 start = static_cast<u64>(start128);
        const u64 end = static_cast<u64>(end128);

        bool ok = true;

        if (start > m || end > m || start > end) {
            ok = false;
        }

        if (!is_hit(p, m, start)) {
            ok = false;
        }

        if (!is_hit(p, m, end)) {
            ok = false;
        }

        if (end < m) {
            if (is_hit(p, m, end + 1)) {
                ok = false;
            }
        }

        if (ok) {
            ++boundary_pass;
        }
    }
}

void print_gap_examples() {
    struct Example {
        u64 p;
        u64 m;
    };

    const std::vector<Example> examples = {
        {2, 20},
        {2, 40},
        {3, 41},
        {5, 194},
        {7, 917}
    };

    std::cout << "\nGAP EXAMPLES\n";

    for (const auto &ex : examples) {
        Parameters par{};

        if (!derive_parameters(ex.p, ex.m, par)) {
            continue;
        }

        std::cout
            << "p=" << par.p
            << " m=" << par.m
            << " e=" << par.e
            << " s0=" << par.s0
            << " q=" << par.q
            << "\n";

        std::cout << "  ";

        bool first = true;

        for (u64 j = 0; j < par.q; ++j) {
            if (!digitwise_leq(par.p, j, par.q)) {
                continue;
            }

            u64 k = 0;
            const u64 gap = predicted_gap(par.p, j, par.q, k);
            const u64 jp = j + gap;

            const u64 start =
                static_cast<u64>(
                    static_cast<u128>(par.s0) +
                    static_cast<u128>(j) * par.pe
                );

            const u64 length =
                static_cast<u64>(
                    static_cast<u128>(gap) * par.pe
                );

            if (!first) {
                std::cout << " | ";
            }

            first = false;

            std::cout
                << "j=" << j
                << " -> " << jp
                << ", k=" << k
                << ", gap=" << gap
                << ", length=" << length
                << ", start=" << start;
        }

        std::cout << "\n";
    }
}

int main() {
    print_start();

    u64 parameter_cases = 0;
    u64 intervals_tested = 0;
    u64 interval_formula_pass = 0;
    u64 gap_formula_pass = 0;

    run_small_exhaustive(
        parameter_cases,
        intervals_tested,
        interval_formula_pass,
        gap_formula_pass
    );

    std::cout << "\nSMALL EXHAUSTIVE\n";
    std::cout
        << "parameter_cases="
        << parameter_cases
        << "\n";

    std::cout
        << "interval_formula_pass="
        << interval_formula_pass
        << "/" << parameter_cases
        << "\n";

    std::cout
        << "gap_cases="
        << intervals_tested
        << "\n";

    std::cout
        << "gap_formula_pass="
        << gap_formula_pass
        << "/" << intervals_tested
        << "\n";

    u64 parameter_pass = 0;
    u64 gap_pass = 0;
    u64 boundary_pass = 0;

    run_large_random(
        200000,
        parameter_pass,
        gap_pass,
        boundary_pass
    );

    std::cout << "\nLARGE RANDOM\n";

    std::cout
        << "parameter_pass="
        << parameter_pass
        << "/200000\n";

    std::cout
        << "gap_pass="
        << gap_pass
        << "/200000\n";

    std::cout
        << "boundary_pass="
        << boundary_pass
        << "/200000\n";

    print_gap_examples();

    std::cout << "\nIDENTITY UNDER TEST\n";
    std::cout
        << "j_next - j = p^k - (j mod p^k)\n";

    std::cout
        << "|I_j| = (p^k - (j mod p^k)) * p^e\n";

    print_finish();

    return 0;
}
