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

u64 calc_start(
    const Parameters &par,
    u64 j
) {
    const u128 value =
        static_cast<u128>(par.s0) +
        static_cast<u128>(j) *
        static_cast<u128>(par.pe);

    return static_cast<u64>(value);
}

u64 calc_endpoint(
    const Parameters &par,
    u64 j
) {
    const u64 r =
        lowest_r(
            par.p,
            j,
            par.q
        );

    if (r == UINT64_MAX) {
        return UINT64_MAX;
    }

    const u64 pr =
        pow_u64(par.p, r);

    const u64 per =
        pow_u64(
            par.p,
            par.e + r
        );

    const u64 high =
        j / pr;

    const u128 value =
        static_cast<u128>(high + 1) *
        static_cast<u128>(per) -
        1;

    return static_cast<u64>(value);
}

std::vector<u64> valid_indices(
    const Parameters &par
) {
    std::vector<u64> result;

    for (u64 j = 0; j < par.q; ++j) {
        if (digitwise_leq(
                par.p,
                j,
                par.q
            )) {
            result.push_back(j);
        }
    }

    return result;
}

void run_small_exhaustive(
    u64 max_m,
    u64 &parameter_cases,
    u64 &intervals_tested,
    u64 &endpoint_pass,
    u64 &endpoint_miss_pass,
    u64 &start_hit_pass,
    u64 &complete_interval_pass
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

            const auto intervals =
                collect_intervals(
                    p,
                    m
                );

            const auto indices =
                valid_indices(par);

            if (intervals.size() !=
                indices.size()) {
                continue;
            }

            for (std::size_t i = 0;
                 i < intervals.size();
                 ++i) {

                ++intervals_tested;

                const u64 j =
                    indices[i];

                const u64 predicted_start_value =
                    calc_start(
                        par,
                        j
                    );

                const u64 predicted_end_value =
                    calc_endpoint(
                        par,
                        j
                    );

                const u64 actual_start =
                    intervals[i].start;

                const u64 actual_end =
                    intervals[i].end;

                if (predicted_end_value ==
                    actual_end) {
                    ++endpoint_pass;
                }

                if (predicted_start_value ==
                        actual_start &&
                    is_hit(
                        par.p,
                        par.m,
                        predicted_start_value
                    )) {
                    ++start_hit_pass;
                }

                if (predicted_end_value < par.m) {
                    if (is_miss(
                            par.p,
                            par.m,
                            predicted_end_value + 1
                        )) {
                        ++endpoint_miss_pass;
                    }
                } else {
                    ++endpoint_miss_pass;
                }

                if (predicted_start_value ==
                        actual_start &&
                    predicted_end_value ==
                        actual_end) {
                    ++complete_interval_pass;
                }
            }
        }
    }
}

u64 random_valid_j(
    u64 p,
    u64 q,
    std::mt19937_64 &rng
) {
    std::vector<u64> digits;

    u64 x = q;

    while (x > 0) {
        digits.push_back(x % p);
        x /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    for (;;) {
        u128 value = 0;
        u128 place = 1;

        for (const u64 qd : digits) {
            std::uniform_int_distribution<u64> dist(
                0,
                qd
            );

            value +=
                static_cast<u128>(
                    dist(rng)
                ) *
                place;

            place *= p;
        }

        const u64 j =
            static_cast<u64>(value);

        if (j < q) {
            return j;
        }
    }
}

void run_large_random(
    u64 trials,
    u64 &start_pass,
    u64 &endpoint_hit_pass,
    u64 &endpoint_miss_pass,
    u64 &divisibility_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    std::mt19937_64 rng(
        0x19020260913ULL
    );

    std::uniform_int_distribution<u64> m_dist(
        1000000000000ULL,
        1000000000000000000ULL
    );

    for (u64 trial = 0;
         trial < trials;
         ++trial) {

        const u64 p =
            primes[trial % primes.size()];

        const u64 m =
            m_dist(rng);

        Parameters par{};

        if (!derive_parameters(
                p,
                m,
                par
            )) {
            --trial;
            continue;
        }

        const u64 j =
            random_valid_j(
                p,
                par.q,
                rng
            );

        const u64 start =
            calc_start(
                par,
                j
            );

        const u64 end =
            calc_endpoint(
                par,
                j
            );

        if (is_hit(
                p,
                m,
                start
            )) {
            ++start_pass;
        }

        if (end <= m) {
            if (is_hit(
                    p,
                    m,
                    end
                )) {
                ++endpoint_hit_pass;
            }

            if (end < m) {
                if (is_miss(
                        p,
                        m,
                        end + 1
                    )) {
                    ++endpoint_miss_pass;
                }
            } else {
                ++endpoint_miss_pass;
            }
        }

        const u64 r =
            lowest_r(
                p,
                j,
                par.q
            );

        const u64 modulus =
            pow_u64(
                p,
                par.e + r
            );

        if ((end + 1) %
                modulus ==
            0) {
            ++divisibility_pass;
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

        const auto indices =
            valid_indices(par);

        std::cout
            << "p=" << ex.p
            << " m=" << ex.m
            << " e=" << par.e
            << " s0=" << par.s0
            << " q=" << par.q
            << "\n";

        for (const u64 j : indices) {
            const u64 r =
                lowest_r(
                    par.p,
                    j,
                    par.q
                );

            const u64 start =
                calc_start(
                    par,
                    j
                );

            const u64 end =
                calc_endpoint(
                    par,
                    j
                );

            std::cout
                << "  j=" << j
                << " r=" << r
                << " start=" << start
                << " end=" << end
                << "\n";
        }
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 190\n";

    u64 parameter_cases = 0;
    u64 intervals_tested = 0;
    u64 endpoint_pass = 0;
    u64 endpoint_miss_pass = 0;
    u64 start_hit_pass = 0;
    u64 complete_interval_pass = 0;

    run_small_exhaustive(
        5000,
        parameter_cases,
        intervals_tested,
        endpoint_pass,
        endpoint_miss_pass,
        start_hit_pass,
        complete_interval_pass
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
        << "endpoint_formula_pass="
        << endpoint_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "endpoint_next_miss_pass="
        << endpoint_miss_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "start_hit_pass="
        << start_hit_pass
        << "/" << intervals_tested
        << "\n";

    std::cout
        << "complete_interval_pass="
        << complete_interval_pass
        << "/" << intervals_tested
        << "\n";

    u64 large_start_pass = 0;
    u64 large_endpoint_hit_pass = 0;
    u64 large_endpoint_miss_pass = 0;
    u64 large_divisibility_pass = 0;

    run_large_random(
        200000,
        large_start_pass,
        large_endpoint_hit_pass,
        large_endpoint_miss_pass,
        large_divisibility_pass
    );

    std::cout
        << "\nLARGE RANDOM\n";

    std::cout
        << "start_hit_pass="
        << large_start_pass
        << "/200000\n";

    std::cout
        << "endpoint_hit_pass="
        << large_endpoint_hit_pass
        << "/200000\n";

    std::cout
        << "endpoint_next_miss_pass="
        << large_endpoint_miss_pass
        << "/200000\n";

    std::cout
        << "endpoint_divisibility_pass="
        << large_divisibility_pass
        << "/200000\n";

    print_examples();

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "E_j = "
        << "(floor(j / p^r) + 1) * p^(e+r) - 1\n";

    std::cout
        << "r = lowest digit position with j_r < q_r\n";

    std::cout
        << "E_j + 1 is divisible by p^(e+r)\n";

    std::cout
        << "FINISHED EXPERIMENT 190\n";

    return 0;
}