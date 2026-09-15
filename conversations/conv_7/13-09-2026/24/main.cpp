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

void start_experiment() {
    std::cout << "START EXPERIMENT 178\n";
}

void finish_experiment() {
    std::cout << "FINISHED EXPERIMENT 178\n";
}

u64 pow_u64(u64 p, u64 e) {
    u128 result = 1;

    for (u64 i = 0; i < e; ++i) {
        result *= p;
    }

    return static_cast<u64>(result);
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

bool is_hit(u64 p, u64 m, u64 t) {
    u64 a = t;
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

    const u128 remainder =
        static_cast<u128>(m) +
        1 -
        static_cast<u128>((b + 1) * pa);

    if (remainder == 0 || remainder % pe != 0) {
        return false;
    }

    const u64 q =
        static_cast<u64>(remainder / pe);

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

std::vector<Interval> collect_actual_intervals(
    u64 p,
    u64 m
) {
    std::vector<Interval> intervals;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit = is_hit(p, m, t);

        if (hit && !inside) {
            inside = true;
            start = t;
        }

        if (!hit && inside) {
            intervals.push_back({start, t - 1});
            inside = false;
        }
    }

    if (inside) {
        intervals.push_back({start, m});
    }

    return intervals;
}

u64 predicted_start(
    const Parameters &par,
    u64 j
) {
    return static_cast<u64>(
        static_cast<u128>(par.s0) +
        static_cast<u128>(j) * par.pe
    );
}

u64 predicted_gap(
    const Parameters &par,
    u64 j
) {
    const u64 jp =
        next_valid_index(
            par.p,
            j,
            par.q
        );

    return jp - j;
}

void print_local_region(
    const Parameters &par,
    u64 centre
) {
    const u64 begin =
        centre > 10 ? centre - 10 : 1;

    const u64 end =
        std::min(
            par.m,
            centre + 10
        );

    for (u64 t = begin; t <= end; ++t) {
        std::cout
            << (is_hit(par.p, par.m, t) ? "H" : ".")
            << t;

        if (t != end) {
            std::cout << " ";
        }
    }

    std::cout << "\n";
}

void run_diagnostics() {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    u64 cases = 0;
    u64 failures = 0;

    for (u64 p : primes) {
        for (u64 m = 2; m <= 2000; ++m) {
            Parameters par{};

            if (!derive_parameters(p, m, par)) {
                continue;
            }

            ++cases;

            const auto actual =
                collect_actual_intervals(p, m);

            u64 actual_index = 0;

            for (u64 j = 0; j < par.q; ++j) {
                if (!digitwise_leq(p, j, par.q)) {
                    continue;
                }

                if (actual_index >= actual.size()) {
                    ++failures;

                    std::cout
                        << "\nCOUNTEREXAMPLE\n"
                        << "p=" << p
                        << " m=" << m
                        << " e=" << par.e
                        << " s0=" << par.s0
                        << " q=" << par.q
                        << "\n"
                        << "j=" << j
                        << " actual_index_exhausted\n";

                    print_local_region(
                        par,
                        predicted_start(par, j)
                    );

                    return;
                }

                const u64 actual_start =
                    actual[actual_index].start;

                const u64 predicted =
                    predicted_start(par, j);

                if (actual_start != predicted) {
                    ++failures;

                    const u64 jp =
                        next_valid_index(
                            p,
                            j,
                            par.q
                        );

                    const u64 gap =
                        predicted_gap(par, j);

                    std::cout
                        << "\nSTART COUNTEREXAMPLE\n"
                        << "p=" << p
                        << " m=" << m
                        << " e=" << par.e
                        << " s0=" << par.s0
                        << " q=" << par.q
                        << "\n"
                        << "j=" << j
                        << " j_next=" << jp
                        << " gap=" << gap
                        << "\n"
                        << "predicted_start="
                        << predicted
                        << "\n"
                        << "actual_start="
                        << actual_start
                        << "\n";

                    print_local_region(
                        par,
                        predicted
                    );

                    return;
                }

                /*
                 * The start is correct. Now inspect the entire
                 * region until the next predicted start.
                 */
                const u64 jp =
                    next_valid_index(
                        p,
                        j,
                        par.q
                    );

                const u64 next_start =
                    predicted_start(par, jp);

                const u64 predicted_end =
                    next_start - 2;

                const u64 actual_end =
                    actual[actual_index].end;

                if (actual_end != predicted_end) {
                    ++failures;

                    std::cout
                        << "\nEND COUNTEREXAMPLE\n"
                        << "p=" << p
                        << " m=" << m
                        << " e=" << par.e
                        << " s0=" << par.s0
                        << " q=" << par.q
                        << "\n"
                        << "j=" << j
                        << " j_next=" << jp
                        << " gap=" << (jp - j)
                        << "\n"
                        << "start="
                        << actual_start
                        << "\n"
                        << "actual_end="
                        << actual_end
                        << "\n"
                        << "predicted_end="
                        << predicted_end
                        << "\n"
                        << "next_start="
                        << next_start
                        << "\n";

                    std::cout
                        << "REGION:\n";

                    print_local_region(
                        par,
                        actual_start
                    );

                    /*
                     * Print every point between the two starts.
                     */
                    std::cout
                        << "POINTS BETWEEN STARTS:\n";

                    for (
                        u64 t = actual_start;
                        t < next_start;
                        ++t
                    ) {
                        std::cout
                            << t
                            << ": "
                            << (
                                is_hit(
                                    p,
                                    m,
                                    t
                                )
                                ? "HIT"
                                : "MISS"
                            )
                            << "\n";
                    }

                    return;
                }

                ++actual_index;
            }
        }
    }

    std::cout
        << "\nNO COUNTEREXAMPLE FOUND\n"
        << "parameter_cases="
        << cases
        << "\n";

    (void)failures;
}

void verify_gap_identity() {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    u64 tested = 0;
    u64 passed = 0;

    for (u64 p : primes) {
        for (u64 m = 2; m <= 2000; ++m) {
            Parameters par{};

            if (!derive_parameters(p, m, par)) {
                continue;
            }

            for (u64 j = 0; j < par.q; ++j) {
                if (!digitwise_leq(
                    p,
                    j,
                    par.q
                )) {
                    continue;
                }

                const u64 k =
                    lowest_available_digit(
                        p,
                        j,
                        par.q
                    );

                const u64 pk =
                    pow_u64(p, k);

                const u64 formula =
                    pk - (j % pk);

                const u64 actual =
                    next_valid_index(
                        p,
                        j,
                        par.q
                    ) - j;

                ++tested;

                if (formula == actual) {
                    ++passed;
                }
            }
        }
    }

    std::cout
        << "\nGAP IDENTITY\n"
        << "tested=" << tested
        << "\n"
        << "passed="
        << passed
        << "/" << tested
        << "\n";
}

int main() {
    start_experiment();

    verify_gap_identity();

    std::cout
        << "\nINTERVAL DIAGNOSTIC\n";

    run_diagnostics();

    finish_experiment();

    return 0;
}
