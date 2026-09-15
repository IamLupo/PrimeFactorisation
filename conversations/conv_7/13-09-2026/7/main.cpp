#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 lo;
    u64 hi;
};

struct PowerInfo {
    u64 value;
    int exponent;
    u64 residual;
};

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

bool lucas_hit(u64 p, u64 m, u64 t) {
    while (m > 0 || t > 0) {
        const u64 md = m % p;
        const u64 td = t % p;

        if (td > md) {
            return true;
        }

        m /= p;
        t /= p;
    }

    return false;
}

std::vector<Interval> build_intervals(u64 p, u64 m) {
    std::vector<Interval> intervals;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit = lucas_hit(p, m, t);

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
        g = std::gcd(
            g,
            intervals[i].lo - intervals[i - 1].lo
        );
    }

    return g;
}

u64 spacing_gcd_ends(const std::vector<Interval>& intervals) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        g = std::gcd(
            g,
            intervals[i].hi - intervals[i - 1].hi
        );
    }

    return g;
}

PowerInfo analyze_power(u64 value, u64 p) {
    PowerInfo result{};
    result.value = value;
    result.exponent = 0;
    result.residual = value;

    if (value == 0) {
        return result;
    }

    while (result.residual % p == 0) {
        result.residual /= p;
        ++result.exponent;
    }

    return result;
}

std::vector<u64> digits_of(u64 n, u64 p) {
    std::vector<u64> digits;

    while (n > 0) {
        digits.push_back(n % p);
        n /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

std::string digits_to_string(const std::vector<u64>& digits) {
    std::string s = "[";

    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (i != 0) {
            s += ",";
        }

        s += std::to_string(digits[i]);
    }

    s += "]";
    return s;
}

std::string truncations_to_string(u64 m, u64 p) {
    std::string s;

    u64 x = m;
    int index = 0;

    while (x > 0) {
        if (!s.empty()) {
            s += " ; ";
        }

        s += "q" + std::to_string(index) + "=";
        s += std::to_string(x);
        s += "(digit=";
        s += std::to_string(x % p);
        s += ")";

        x /= p;
        ++index;
    }

    return s;
}

bool is_pure_power(u64 value, u64 p) {
    if (value == 0) {
        return false;
    }

    while (value % p == 0) {
        value /= p;
    }

    return value == 1;
}

void inspect_case(u64 p, u64 m) {
    const auto intervals = build_intervals(p, m);

    if (intervals.size() < 2) {
        return;
    }

    const u64 gs = spacing_gcd_starts(intervals);
    const u64 ge = spacing_gcd_ends(intervals);

    const PowerInfo ps = analyze_power(gs, p);
    const PowerInfo pe = analyze_power(ge, p);

    std::cout
        << "p=" << p
        << " m=" << m
        << " intervals=" << intervals.size()
        << " e=" << ps.exponent
        << " digits=" << digits_to_string(digits_of(m, p))
        << "\n";

    std::cout
        << "  start G=" << gs
        << " end G=" << ge
        << " residuals=("
        << ps.residual << ","
        << pe.residual << ")\n";

    std::cout
        << "  "
        << truncations_to_string(m, p)
        << "\n";
}

void run_small_exhaustive() {
    std::cout << "\nPHASE 1: EXHAUSTIVE DIGIT TABLE\n";

    u64 total = 0;
    u64 max_exponent = 0;

    for (u64 p = 2; p <= 31; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 500; ++m) {
            const auto intervals = build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const u64 gs = spacing_gcd_starts(intervals);
            const PowerInfo ps = analyze_power(gs, p);

            max_exponent =
                std::max<u64>(
                    max_exponent,
                    static_cast<u64>(ps.exponent)
                );
        }
    }

    std::cout << "multi_interval_cases=" << total << "\n";
    std::cout << "maximum_exponent=" << max_exponent << "\n";

    /*
        Print representative cases for each exponent.
    */
    for (int target_e = 1; target_e <= 6; ++target_e) {
        bool found = false;

        for (u64 p = 2; p <= 31 && !found; ++p) {
            if (!is_prime(p)) {
                continue;
            }

            for (u64 m = 1; m <= 500 && !found; ++m) {
                const auto intervals = build_intervals(p, m);

                if (intervals.size() < 2) {
                    continue;
                }

                const u64 gs = spacing_gcd_starts(intervals);
                const PowerInfo ps = analyze_power(gs, p);

                if (ps.exponent == target_e) {
                    std::cout
                        << "\nREPRESENTATIVE e="
                        << target_e << "\n";

                    inspect_case(p, m);
                    found = true;
                }
            }
        }
    }
}

void run_digit_pattern_phase() {
    std::cout << "\nPHASE 2: DIGIT PATTERN SEARCH\n";

    /*
        Search larger m while deliberately covering
        different low-order digit patterns.
    */
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31
    };

    for (u64 p : primes) {
        std::cout << "\nPRIME p=" << p << "\n";

        u64 printed = 0;

        for (u64 m = 1; m <= 5000; ++m) {
            const auto intervals = build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            const u64 gs = spacing_gcd_starts(intervals);
            const PowerInfo ps = analyze_power(gs, p);

            /*
                Only show cases where e >= 2.
            */
            if (ps.exponent >= 2) {
                inspect_case(p, m);

                ++printed;

                if (printed >= 20) {
                    break;
                }
            }
        }
    }
}

void run_truncation_comparison() {
    std::cout << "\nPHASE 3: TRUNCATION COMPARISON\n";

    u64 total = 0;
    u64 recursive_matches = 0;

    /*
        For each case compute the exponent of every
        truncation floor(m/p^j).

        This does NOT impose a formula. It records the
        sequence so we can inspect whether e(m) is related
        recursively to e(floor(m/p)).
    */
    for (u64 p = 2; p <= 47; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 3000; ++m) {
            const auto intervals = build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const u64 gs = spacing_gcd_starts(intervals);
            const PowerInfo ps = analyze_power(gs, p);

            const u64 q = m / p;

            const auto q_intervals =
                build_intervals(p, q);

            int e_q = 1;

            if (q_intervals.size() >= 2) {
                const u64 gq =
                    spacing_gcd_starts(q_intervals);

                const PowerInfo pq =
                    analyze_power(gq, p);

                e_q = pq.exponent;
            }

            /*
                Print cases where the relationship between
                e(m) and e(floor(m/p)) is unusual.
            */
            if (ps.exponent != e_q &&
                ps.exponent != e_q + 1) {

                if (recursive_matches < 30) {
                    std::cout
                        << "UNUSUAL "
                        << "p=" << p
                        << " m=" << m
                        << " digits="
                        << digits_to_string(
                            digits_of(m, p)
                        )
                        << " e(m)="
                        << ps.exponent
                        << " e(floor(m/p))="
                        << e_q
                        << "\n";
                }
            } else {
                ++recursive_matches;
            }
        }
    }

    std::cout << "\nSUMMARY PHASE 3\n";
    std::cout << "cases=" << total << "\n";
    std::cout << "simple_recursive_relation_count="
              << recursive_matches << "\n";
}

void run_high_exponent_search() {
    std::cout << "\nPHASE 4: SEARCH FOR LARGE e\n";

    u64 total = 0;
    u64 max_e = 0;

    for (u64 p = 2; p <= 19; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 20000; ++m) {
            const auto intervals = build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const u64 gs = spacing_gcd_starts(intervals);
            const PowerInfo ps = analyze_power(gs, p);

            if (ps.exponent > max_e) {
                max_e = ps.exponent;

                std::cout
                    << "NEW_MAX "
                    << "p=" << p
                    << " m=" << m
                    << " e=" << ps.exponent
                    << " G=" << gs
                    << " digits="
                    << digits_to_string(
                        digits_of(m, p)
                    )
                    << "\n";
            }
        }
    }

    std::cout << "\nSUMMARY PHASE 4\n";
    std::cout << "multi_interval_cases=" << total << "\n";
    std::cout << "maximum_exponent=" << max_e << "\n";
}

int main() {
    constexpr int EXPERIMENT = 163;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_small_exhaustive();
    run_digit_pattern_phase();
    run_truncation_comparison();
    run_high_exponent_search();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
