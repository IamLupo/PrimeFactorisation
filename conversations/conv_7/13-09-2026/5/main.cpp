#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 lo;
    u64 hi;
};

struct PowerCheck {
    u64 value;
    bool divisible_by_p;
    bool pure_power;
    int exponent;
    u64 residual;
};

u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

std::vector<int> base_p_digits(u64 x, u64 p) {
    std::vector<int> digits;

    while (x > 0) {
        digits.push_back(static_cast<int>(x % p));
        x /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

std::string digits_to_string(const std::vector<int>& digits) {
    std::string out = "[";

    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i != 0) {
            out += ",";
        }

        out += std::to_string(digits[i]);
    }

    out += "]";
    return out;
}

/*
    Lucas criterion for

        p | C(m, m-t)

    Since C(m,m-t) = C(m,t), divisibility fails exactly when every
    base-p digit of t is <= the corresponding digit of m.
*/
bool hit_by_lucas(u64 p, u64 m, u64 t) {
    while (m > 0 || t > 0) {
        u64 md = m % p;
        u64 td = t % p;

        if (td > md) {
            return true;
        }

        m /= p;
        t /= p;
    }

    return false;
}

std::vector<Interval> build_hit_intervals(u64 p, u64 m) {
    std::vector<Interval> intervals;

    if (m == 0) {
        return intervals;
    }

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        bool hit = hit_by_lucas(p, m, t);

        if (hit && !inside) {
            start = t;
            inside = true;
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

u64 spacing_gcd_starts(const std::vector<Interval>& intervals) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        u64 diff = intervals[i].lo - intervals[i - 1].lo;
        g = gcd_u64(g, diff);
    }

    return g;
}

u64 spacing_gcd_ends(const std::vector<Interval>& intervals) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        u64 diff = intervals[i].hi - intervals[i - 1].hi;
        g = gcd_u64(g, diff);
    }

    return g;
}

PowerCheck analyze_power(u64 value, u64 p) {
    PowerCheck result{};
    result.value = value;
    result.divisible_by_p = (value != 0 && value % p == 0);
    result.pure_power = false;
    result.exponent = 0;
    result.residual = value;

    if (!result.divisible_by_p) {
        return result;
    }

    while (result.residual % p == 0) {
        result.residual /= p;
        ++result.exponent;
    }

    result.pure_power = (result.residual == 1);
    return result;
}

void print_exception(
    u64 p,
    u64 m,
    const std::vector<Interval>& intervals,
    const PowerCheck& starts,
    const PowerCheck& ends
) {
    std::cout << "\nEXCEPTION\n";
    std::cout << "p=" << p << "\n";
    std::cout << "m=" << m << "\n";
    std::cout << "m_base_p=" << digits_to_string(base_p_digits(m, p)) << "\n";
    std::cout << "interval_count=" << intervals.size() << "\n";

    std::cout << "start_gcd=" << starts.value
              << " pure_power=" << starts.pure_power
              << " exponent=" << starts.exponent
              << " residual=" << starts.residual << "\n";

    std::cout << "end_gcd=" << ends.value
              << " pure_power=" << ends.pure_power
              << " exponent=" << ends.exponent
              << " residual=" << ends.residual << "\n";

    std::cout << "first_intervals=";

    std::size_t limit = std::min<std::size_t>(intervals.size(), 8);

    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << " ";
        }

        std::cout << "[" << intervals[i].lo
                  << "," << intervals[i].hi << "]";
    }

    std::cout << "\n";
}

void run_exhaustive_phase(
    int max_prime,
    u64 max_m
) {
    std::cout << "\nPHASE 1: EXHAUSTIVE SMALL TEST\n";

    u64 total_multi = 0;
    u64 no_start_gcd = 0;
    u64 start_divisible = 0;
    u64 end_divisible = 0;
    u64 start_pure = 0;
    u64 end_pure = 0;
    u64 exponent_equal = 0;
    u64 exponent_difference = 0;

    u64 exceptions = 0;

    for (int p = 2; p <= max_prime; ++p) {
        bool prime = true;

        for (int d = 2; d * d <= p; ++d) {
            if (p % d == 0) {
                prime = false;
                break;
            }
        }

        if (!prime) {
            continue;
        }

        for (u64 m = 1; m <= max_m; ++m) {
            std::vector<Interval> intervals = build_hit_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total_multi;

            u64 gs = spacing_gcd_starts(intervals);
            u64 ge = spacing_gcd_ends(intervals);

            if (gs == 0 || ge == 0) {
                ++no_start_gcd;
                continue;
            }

            PowerCheck starts = analyze_power(gs, p);
            PowerCheck ends = analyze_power(ge, p);

            if (starts.divisible_by_p) {
                ++start_divisible;
            }

            if (ends.divisible_by_p) {
                ++end_divisible;
            }

            if (starts.pure_power) {
                ++start_pure;
            }

            if (ends.pure_power) {
                ++end_pure;
            }

            if (starts.pure_power &&
                ends.pure_power &&
                starts.exponent == ends.exponent) {
                ++exponent_equal;
            } else {
                ++exponent_difference;

                if (exceptions < 20) {
                    print_exception(
                        p,
                        m,
                        intervals,
                        starts,
                        ends
                    );
                }

                ++exceptions;
            }
        }
    }

    std::cout << "\nSUMMARY PHASE 1\n";
    std::cout << "multi_interval_cases=" << total_multi << "\n";
    std::cout << "no_gcd_cases=" << no_start_gcd << "\n";
    std::cout << "start_divisible_by_p="
              << start_divisible << "/" << total_multi << "\n";
    std::cout << "end_divisible_by_p="
              << end_divisible << "/" << total_multi << "\n";
    std::cout << "start_pure_power="
              << start_pure << "/" << total_multi << "\n";
    std::cout << "end_pure_power="
              << end_pure << "/" << total_multi << "\n";
    std::cout << "same_exponent="
              << exponent_equal << "/" << total_multi << "\n";
    std::cout << "exception_cases=" << exponent_difference << "\n";
}

void run_targeted_phase() {
    std::cout << "\nPHASE 2: TARGETED BASE-p STRUCTURES\n";

    struct Case {
        u64 p;
        u64 m;
    };

    std::vector<Case> cases = {
        {5, 59},
        {11, 395},
        {17, 917},
        {13, 181},
        {19, 341},
        {23, 529},
        {29, 899},
        {31, 1023},
        {37, 1369},
        {41, 1681},
        {43, 1848},
        {47, 2208}
    };

    for (const Case& c : cases) {
        std::vector<Interval> intervals =
            build_hit_intervals(c.p, c.m);

        u64 gs = spacing_gcd_starts(intervals);
        u64 ge = spacing_gcd_ends(intervals);

        PowerCheck starts = analyze_power(gs, c.p);
        PowerCheck ends = analyze_power(ge, c.p);

        std::cout << "\np=" << c.p
                  << " m=" << c.m
                  << " digits="
                  << digits_to_string(base_p_digits(c.m, c.p))
                  << " intervals=" << intervals.size()
                  << "\n";

        if (intervals.size() >= 2) {
            std::cout << "start_gcd=" << gs
                      << " = p^" << starts.exponent
                      << " residual=" << starts.residual << "\n";

            std::cout << "end_gcd=" << ge
                      << " = p^" << ends.exponent
                      << " residual=" << ends.residual << "\n";
        } else {
            std::cout << "not enough intervals for spacing GCD\n";
        }
    }
}

void run_random_phase(
    int cases,
    u64 seed
) {
    std::cout << "\nPHASE 3: RANDOM LARGER CASES\n";

    std::mt19937_64 rng(seed);

    std::vector<u64> primes = {
        101, 127, 163, 211, 257,
        331, 367, 431, 509, 743,
        1009, 2003, 3001, 4001, 5003
    };

    u64 valid = 0;
    u64 pure_power_pass = 0;
    u64 exponent_match = 0;
    u64 failures_printed = 0;

    std::uniform_int_distribution<int> prime_dist(
        0,
        static_cast<int>(primes.size()) - 1
    );

    std::uniform_int_distribution<u64> multiplier_dist(
        10,
        200
    );

    for (int c = 0; c < cases; ++c) {
        u64 p = primes[prime_dist(rng)];
        u64 m = p * multiplier_dist(rng) + (rng() % p);

        std::vector<Interval> intervals =
            build_hit_intervals(p, m);

        if (intervals.size() < 2) {
            continue;
        }

        ++valid;

        u64 gs = spacing_gcd_starts(intervals);
        u64 ge = spacing_gcd_ends(intervals);

        PowerCheck starts = analyze_power(gs, p);
        PowerCheck ends = analyze_power(ge, p);

        bool both_pure =
            starts.pure_power &&
            ends.pure_power;

        bool same_exp =
            both_pure &&
            starts.exponent == ends.exponent;

        if (both_pure) {
            ++pure_power_pass;
        }

        if (same_exp) {
            ++exponent_match;
        }

        if ((!both_pure || !same_exp) &&
            failures_printed < 10) {

            print_exception(
                p,
                m,
                intervals,
                starts,
                ends
            );

            ++failures_printed;
        }

        std::cout << "case=" << c + 1
                  << " p=" << p
                  << " m=" << m
                  << " intervals=" << intervals.size()
                  << " Gs=" << gs
                  << " Ge=" << ge
                  << " es=" << starts.exponent
                  << " ee=" << ends.exponent
                  << "\n";
    }

    std::cout << "\nSUMMARY PHASE 3\n";
    std::cout << "valid_multi_interval_cases=" << valid << "\n";
    std::cout << "both_pure_power="
              << pure_power_pass << "/" << valid << "\n";
    std::cout << "same_exponent="
              << exponent_match << "/" << valid << "\n";
}

void run_edge_cases() {
    std::cout << "\nPHASE 4: EDGE CASES\n";

    struct Case {
        u64 p;
        u64 m;
    };

    std::vector<Case> cases = {
        {2, 4},
        {2, 7},
        {3, 8},
        {3, 17},
        {5, 24},
        {5, 25},
        {5, 59},
        {7, 48},
        {7, 342},
        {11, 120},
        {11, 395},
        {17, 917}
    };

    for (const Case& c : cases) {
        std::vector<Interval> intervals =
            build_hit_intervals(c.p, c.m);

        u64 gs = spacing_gcd_starts(intervals);
        u64 ge = spacing_gcd_ends(intervals);

        PowerCheck starts = analyze_power(gs, c.p);
        PowerCheck ends = analyze_power(ge, c.p);

        std::cout << "p=" << c.p
                  << " m=" << c.m
                  << " digits="
                  << digits_to_string(base_p_digits(c.m, c.p))
                  << " intervals=" << intervals.size();

        if (intervals.size() >= 2) {
            std::cout << " Gs=" << gs
                      << " Es=" << starts.exponent
                      << " Rs=" << starts.residual
                      << " Ge=" << ge
                      << " Ee=" << ends.exponent
                      << " Re=" << ends.residual;
        }

        std::cout << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 161;

    std::cout << "START EXPERIMENT " << EXPERIMENT << "\n";

    run_exhaustive_phase(97, 500);
    run_targeted_phase();
    run_random_phase(40, 0x1612026ULL);
    run_edge_cases();

    std::cout << "\nFINISHED EXPERIMENT " << EXPERIMENT << "\n";

    return 0;
}
