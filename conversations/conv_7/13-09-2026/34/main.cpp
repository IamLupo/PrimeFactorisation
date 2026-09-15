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

bool is_miss(u64 p, u64 m, u64 x) {
    return !is_hit(p, m, x);
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
 * Produce all predicted HIT starts.

 * j is a valid start index iff:

 *   j < q
 *   j <=_p q
 *
 * and

 *   x_j = s0 + j*p^e.
 */
std::vector<u64> predicted_starts(
    const Parameters &par
) {
    std::vector<u64> result;

    /*
     * q <= 5000 in the exhaustive experiment, so
     * enumerating j is completely bounded here.
     */
    for (u64 j = 0; j < par.q; ++j) {
        if (!digitwise_leq(
                par.p,
                j,
                par.q
            )) {
            continue;
        }

        const u128 x =
            static_cast<u128>(par.s0) +
            static_cast<u128>(j) *
                static_cast<u128>(par.pe);

        if (x <= par.m) {
            result.push_back(
                static_cast<u64>(x)
            );
        }
    }

    return result;
}

bool strictly_increasing(
    const std::vector<u64> &values
) {
    for (std::size_t i = 1;
         i < values.size();
         ++i) {

        if (values[i - 1] >= values[i]) {
            return false;
        }
    }

    return true;
}

bool same_vector(
    const std::vector<u64> &a,
    const std::vector<u64> &b
) {
    return a == b;
}

void run_small_exhaustive(
    u64 max_m,
    u64 &parameter_cases,
    u64 &cases,
    u64 &exact_set_pass,
    u64 &actual_starts_pass,
    u64 &predicted_starts_pass,
    u64 &count_pass,
    u64 &ordering_pass,
    u64 &predecessor_miss_pass
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

            ++parameter_cases;

            const auto actual_intervals =
                collect_intervals(
                    p,
                    m
                );

            std::vector<u64> actual_starts;

            for (const auto &interval :
                 actual_intervals) {

                actual_starts.push_back(
                    interval.start
                );
            }

            const auto predicted =
                predicted_starts(par);

            ++cases;

            /*
             * Entire set equality.
             */
            if (same_vector(
                    actual_starts,
                    predicted
                )) {
                ++exact_set_pass;
            }

            /*
             * Every actual start belongs to the predicted
             * arithmetic/lucas construction.
             */
            bool actual_ok = true;

            for (const u64 x :
                 actual_starts) {

                if (x < par.s0) {
                    actual_ok = false;
                    break;
                }

                const u64 delta =
                    x - par.s0;

                if (delta % par.pe != 0) {
                    actual_ok = false;
                    break;
                }

                const u64 j =
                    delta / par.pe;

                if (j >= par.q ||
                    !digitwise_leq(
                        par.p,
                        j,
                        par.q
                    )) {
                    actual_ok = false;
                    break;
                }
            }

            if (actual_ok) {
                ++actual_starts_pass;
            }

            /*
             * Every predicted start is genuinely a HIT
             * and begins an interval.
             */
            bool predicted_ok = true;

            for (const u64 x :
                 predicted) {

                if (!is_hit(
                        par.p,
                        par.m,
                        x
                    )) {
                    predicted_ok = false;
                    break;
                }

                if (x > 1 &&
                    is_hit(
                        par.p,
                        par.m,
                        x - 1
                    )) {
                    predicted_ok = false;
                    break;
                }
            }

            if (predicted_ok) {
                ++predicted_starts_pass;
            }

            /*
             * Cardinality independently.
             */
            if (actual_starts.size() ==
                predicted.size()) {
                ++count_pass;
            }

            /*
             * The arithmetic progression x_j is strictly
             * increasing because p^e > 0.
             */
            if (strictly_increasing(predicted)) {
                ++ordering_pass;
            }

            /*
             * Every actual start must have a MISS immediately
             * before it, except possibly start=1.
             */
            bool predecessor_ok = true;

            for (const u64 x :
                 actual_starts) {

                if (x > 1 &&
                    !is_miss(
                        par.p,
                        par.m,
                        x - 1
                    )) {
                    predecessor_ok = false;
                    break;
                }
            }

            if (predecessor_ok) {
                ++predecessor_miss_pass;
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

        const auto actual =
            collect_intervals(
                ex.p,
                ex.m
            );

        const auto predicted =
            predicted_starts(par);

        std::cout
            << "p=" << ex.p
            << " m=" << ex.m
            << " e=" << par.e
            << " s0=" << par.s0
            << " q=" << par.q
            << "\n";

        std::cout
            << "  actual_starts:";

        for (const auto &interval : actual) {
            std::cout
                << " "
                << interval.start;
        }

        std::cout
            << "\n";

        std::cout
            << "  predicted_starts:";

        for (const u64 x : predicted) {
            std::cout
                << " "
                << x;
        }

        std::cout
            << "\n";
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 189\n";

    u64 parameter_cases = 0;
    u64 cases = 0;
    u64 exact_set_pass = 0;
    u64 actual_starts_pass = 0;
    u64 predicted_starts_pass = 0;
    u64 count_pass = 0;
    u64 ordering_pass = 0;
    u64 predecessor_miss_pass = 0;

    run_small_exhaustive(
        5000,
        parameter_cases,
        cases,
        exact_set_pass,
        actual_starts_pass,
        predicted_starts_pass,
        count_pass,
        ordering_pass,
        predecessor_miss_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << parameter_cases
        << "\n";

    std::cout
        << "theorem_set_equality_pass="
        << exact_set_pass
        << "/" << cases
        << "\n";

    std::cout
        << "actual_start_characterization_pass="
        << actual_starts_pass
        << "/" << cases
        << "\n";

    std::cout
        << "predicted_start_validity_pass="
        << predicted_starts_pass
        << "/" << cases
        << "\n";

    std::cout
        << "start_count_pass="
        << count_pass
        << "/" << cases
        << "\n";

    std::cout
        << "ordering_pass="
        << ordering_pass
        << "/" << cases
        << "\n";

    std::cout
        << "predecessor_miss_pass="
        << predecessor_miss_pass
        << "/" << cases
        << "\n";

    print_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "HIT_STARTS = "
        << "{ s0 + j*p^e : 0 <= j < q, "
        << "j <=_p q }\n";

    std::cout
        << "FINISHED EXPERIMENT 189\n";

    return 0;
}
