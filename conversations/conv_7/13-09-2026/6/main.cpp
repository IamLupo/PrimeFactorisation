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

std::vector<Interval> build_intervals(u64 p, u64 m) {
    std::vector<Interval> intervals;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        bool hit = lucas_hit(p, m, t);

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
        u64 d = intervals[i].lo - intervals[i - 1].lo;
        g = std::gcd(g, d);
    }

    return g;
}

u64 spacing_gcd_ends(const std::vector<Interval>& intervals) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        u64 d = intervals[i].hi - intervals[i - 1].hi;
        g = std::gcd(g, d);
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

std::string base_p_digits(u64 n, u64 p) {
    std::vector<u64> digits;

    while (n > 0) {
        digits.push_back(n % p);
        n /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

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

void print_failure(
    const std::string& reason,
    u64 p,
    u64 m,
    const std::vector<Interval>& intervals,
    u64 gs,
    u64 ge,
    const PowerInfo& ps,
    const PowerInfo& pe
) {
    std::cout << "\nFAILURE: " << reason << "\n";
    std::cout << "p=" << p << "\n";
    std::cout << "m=" << m << "\n";
    std::cout << "m_mod_p=" << (m % p) << "\n";
    std::cout << "digits=" << base_p_digits(m, p) << "\n";
    std::cout << "intervals=" << intervals.size() << "\n";

    std::cout << "Gs=" << gs
              << " exponent=" << ps.exponent
              << " residual=" << ps.residual << "\n";

    std::cout << "Ge=" << ge
              << " exponent=" << pe.exponent
              << " residual=" << pe.residual << "\n";

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

void run_exhaustive_phase() {
    std::cout << "\nPHASE 1: EXPONENT CLASSIFICATION\n";

    u64 total = 0;
    u64 predicted = 0;
    u64 exponent_one = 0;
    u64 exponent_two = 0;
    u64 exponent_three_or_more = 0;

    u64 pure_power_pass = 0;
    u64 equal_exponent_pass = 0;
    u64 classification_failures = 0;

    for (u64 p = 2; p <= 149; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 3000; ++m) {
            std::vector<Interval> intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            u64 gs = spacing_gcd_starts(intervals);
            u64 ge = spacing_gcd_ends(intervals);

            PowerInfo ps = analyze_power(gs, p);
            PowerInfo pe = analyze_power(ge, p);

            bool pure =
                ps.residual == 1 &&
                pe.residual == 1;

            bool equal =
                ps.exponent == pe.exponent;

            if (pure) {
                ++pure_power_pass;
            }

            if (equal) {
                ++equal_exponent_pass;
            }

            if (ps.exponent == 1) {
                ++exponent_one;
            }

            if (ps.exponent == 2) {
                ++exponent_two;
            }

            if (ps.exponent >= 3) {
                ++exponent_three_or_more;
            }

            /*
                Proposed rule:

                    m mod p = p-1  -> e=2
                    otherwise      -> e=1
            */
            int predicted_exponent =
                (m % p == p - 1) ? 2 : 1;

            bool ok =
                pure &&
                equal &&
                ps.exponent == predicted_exponent;

            if (ok) {
                ++predicted;
            } else {
                ++classification_failures;

                if (classification_failures <= 20) {
                    print_failure(
                        "exponent classification",
                        p,
                        m,
                        intervals,
                        gs,
                        ge,
                        ps,
                        pe
                    );
                }
            }
        }
    }

    std::cout << "\nSUMMARY PHASE 1\n";
    std::cout << "multi_interval_cases=" << total << "\n";
    std::cout << "pure_power_pass="
              << pure_power_pass << "/" << total << "\n";
    std::cout << "equal_exponent_pass="
              << equal_exponent_pass << "/" << total << "\n";
    std::cout << "e=1=" << exponent_one << "\n";
    std::cout << "e=2=" << exponent_two << "\n";
    std::cout << "e>=3=" << exponent_three_or_more << "\n";
    std::cout << "classification_pass="
              << predicted << "/" << total << "\n";
    std::cout << "classification_failures="
              << classification_failures << "\n";
}

void run_mod_minus_one_phase() {
    std::cout << "\nPHASE 2: m = -1 mod p TARGETS\n";

    u64 tested = 0;
    u64 multi = 0;
    u64 e2 = 0;
    u64 e3 = 0;

    for (u64 p = 2; p <= 200; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        /*
            Test several values with m = kp - 1.
        */
        for (u64 k = 2; k <= 80; ++k) {
            u64 m = k * p - 1;

            ++tested;

            std::vector<Interval> intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++multi;

            u64 gs = spacing_gcd_starts(intervals);
            PowerInfo ps = analyze_power(gs, p);

            if (ps.exponent == 2) {
                ++e2;
            }

            if (ps.exponent >= 3) {
                ++e3;

                std::cout
                    << "HIGH_EXPONENT "
                    << "p=" << p
                    << " m=" << m
                    << " k=" << k
                    << " intervals=" << intervals.size()
                    << " G=" << gs
                    << " e=" << ps.exponent
                    << " residual=" << ps.residual
                    << "\n";
            }
        }
    }

    std::cout << "\nSUMMARY PHASE 2\n";
    std::cout << "tested_m_equals_minus_one_mod_p=" << tested << "\n";
    std::cout << "multi_interval=" << multi << "\n";
    std::cout << "e=2=" << e2 << "\n";
    std::cout << "e>=3=" << e3 << "\n";
}

void run_non_minus_one_phase() {
    std::cout << "\nPHASE 3: m != -1 mod p TARGETS\n";

    u64 total = 0;
    u64 e1 = 0;
    u64 e2_or_more = 0;

    for (u64 p = 2; p <= 200; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 k = 2; k <= 100; ++k) {
            for (u64 r = 0; r + 1 < p; ++r) {
                u64 m = k * p + r;

                std::vector<Interval> intervals =
                    build_intervals(p, m);

                if (intervals.size() < 2) {
                    continue;
                }

                ++total;

                u64 gs = spacing_gcd_starts(intervals);
                PowerInfo ps = analyze_power(gs, p);

                if (ps.exponent == 1) {
                    ++e1;
                } else {
                    ++e2_or_more;

                    if (e2_or_more <= 20) {
                        std::cout
                            << "NON_EXPECTED "
                            << "p=" << p
                            << " m=" << m
                            << " m_mod_p=" << r
                            << " intervals=" << intervals.size()
                            << " G=" << gs
                            << " e=" << ps.exponent
                            << " residual=" << ps.residual
                            << "\n";
                    }
                }
            }
        }
    }

    std::cout << "\nSUMMARY PHASE 3\n";
    std::cout << "multi_interval_cases=" << total << "\n";
    std::cout << "e=1=" << e1 << "\n";
    std::cout << "e>=2=" << e2_or_more << "\n";
}

void run_large_targeted_phase() {
    std::cout << "\nPHASE 4: LARGER TARGETS\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    std::vector<TestCase> cases = {
        {101, 1058},
        {127, 127 * 17 - 1},
        {211, 211 * 37 - 1},
        {431, 431 * 91 - 1},
        {743, 743 * 113 - 1},
        {1009, 1009 * 127 - 1},
        {2003, 2003 * 151 - 1},
        {4001, 4001 * 177 - 1},
        {5003, 5003 * 193 - 1},

        {101, 1057},
        {127, 127 * 17 - 2},
        {211, 211 * 37 - 2},
        {431, 431 * 91 - 2},
        {743, 743 * 113 - 2},
        {1009, 1009 * 127 - 2},
        {2003, 2003 * 151 - 2}
    };

    for (const TestCase& c : cases) {
        std::vector<Interval> intervals =
            build_intervals(c.p, c.m);

        if (intervals.size() < 2) {
            std::cout
                << "p=" << c.p
                << " m=" << c.m
                << " intervals=" << intervals.size()
                << " insufficient_intervals\n";
            continue;
        }

        u64 gs = spacing_gcd_starts(intervals);
        u64 ge = spacing_gcd_ends(intervals);

        PowerInfo ps = analyze_power(gs, c.p);
        PowerInfo pe = analyze_power(ge, c.p);

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " m_mod_p=" << (c.m % c.p)
            << " intervals=" << intervals.size()
            << " Gs=" << gs
            << " es=" << ps.exponent
            << " Ge=" << ge
            << " ee=" << pe.exponent
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 162;

    std::cout << "START EXPERIMENT "
              << EXPERIMENT << "\n";

    run_exhaustive_phase();
    run_mod_minus_one_phase();
    run_non_minus_one_phase();
    run_large_targeted_phase();

    std::cout << "\nFINISHED EXPERIMENT "
              << EXPERIMENT << "\n";

    return 0;
}
