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

u64 start_spacing_gcd(
    const std::vector<Interval>& intervals
) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        const u64 diff =
            intervals[i].lo - intervals[i - 1].lo;

        g = std::gcd(g, diff);
    }

    return g;
}

u64 end_spacing_gcd(
    const std::vector<Interval>& intervals
) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1; i < intervals.size(); ++i) {
        const u64 diff =
            intervals[i].hi - intervals[i - 1].hi;

        g = std::gcd(g, diff);
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

std::vector<int> base_p_digits(
    u64 m,
    u64 p
) {
    std::vector<int> digits;

    while (m > 0) {
        digits.push_back(
            static_cast<int>(m % p)
        );

        m /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

std::string digits_string(
    const std::vector<int>& digits
) {
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
    Return the first i such that

        d_i < d_{i+1}.

    Digits are least-significant first.

    Conjectured exponent:

        e = i + 1.

    Return -1 when no ascent exists.
*/
int first_ascent_position(
    const std::vector<int>& digits
) {
    if (digits.size() < 2) {
        return -1;
    }

    for (std::size_t i = 0;
         i + 1 < digits.size();
         ++i) {

        if (digits[i] < digits[i + 1]) {
            return static_cast<int>(i);
        }
    }

    return -1;
}

void print_case(
    u64 p,
    u64 m,
    const std::vector<int>& digits,
    const std::vector<Interval>& intervals,
    int observed_e,
    int predicted_e,
    int ascent
) {
    std::cout
        << "\nCASE\n";

    std::cout
        << "p=" << p
        << " m=" << m
        << " digits="
        << digits_string(digits)
        << "\n";

    std::cout
        << "intervals="
        << intervals.size()
        << "\n";

    std::cout
        << "first_ascent=";

    if (ascent < 0) {
        std::cout << "NONE";
    } else {
        std::cout << ascent;
    }

    std::cout
        << " observed_e="
        << observed_e
        << " predicted_e="
        << predicted_e
        << "\n";

    std::cout
        << "starts=";

    const std::size_t limit =
        std::min<std::size_t>(
            intervals.size(),
            12
        );

    for (std::size_t i = 0;
         i < limit;
         ++i) {

        if (i != 0) {
            std::cout << ",";
        }

        std::cout << intervals[i].lo;
    }

    if (intervals.size() > limit) {
        std::cout << ",...";
    }

    std::cout << "\n";
}

void run_exhaustive_phase() {
    std::cout
        << "\nPHASE 1: EXHAUSTIVE FIRST-ASCENT TEST\n";

    u64 total = 0;
    u64 predicted = 0;
    u64 failed = 0;

    u64 no_ascent = 0;
    u64 no_ascent_multi = 0;

    u64 pure_power_fail = 0;
    u64 start_end_exponent_mismatch = 0;

    int maximum_e = 0;

    for (u64 p = 2; p <= 101; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 5000; ++m) {
            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const auto digits =
                base_p_digits(m, p);

            const int ascent =
                first_ascent_position(digits);

            if (ascent < 0) {
                ++no_ascent;
                ++no_ascent_multi;
            }

            /*
                Main conjecture.

                For a normal case:
                    e = first_ascent + 1

                For no-ascent cases we temporarily predict
                e=1 so that they are explicitly visible as
                potential exceptions.
            */
            const int predicted_e =
                (ascent >= 0)
                    ? ascent + 1
                    : 1;

            const u64 gs =
                start_spacing_gcd(intervals);

            const u64 ge =
                end_spacing_gcd(intervals);

            const PowerInfo ps =
                analyze_power(gs, p);

            const PowerInfo pe =
                analyze_power(ge, p);

            if (ps.residual != 1 ||
                pe.residual != 1) {

                ++pure_power_fail;
            }

            if (ps.exponent != pe.exponent) {
                ++start_end_exponent_mismatch;
            }

            maximum_e =
                std::max(
                    maximum_e,
                    ps.exponent
                );

            const bool ok =
                ps.residual == 1 &&
                pe.residual == 1 &&
                ps.exponent == pe.exponent &&
                ps.exponent == predicted_e;

            if (ok) {
                ++predicted;
            } else {
                ++failed;

                if (failed <= 40) {
                    print_case(
                        p,
                        m,
                        digits,
                        intervals,
                        ps.exponent,
                        predicted_e,
                        ascent
                    );
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n";

    std::cout
        << "multi_interval_cases="
        << total
        << "\n";

    std::cout
        << "prediction_pass="
        << predicted
        << "/" << total
        << "\n";

    std::cout
        << "prediction_fail="
        << failed
        << "\n";

    std::cout
        << "no_ascent_cases="
        << no_ascent
        << "\n";

    std::cout
        << "pure_power_fail="
        << pure_power_fail
        << "\n";

    std::cout
        << "start_end_exponent_mismatch="
        << start_end_exponent_mismatch
        << "\n";

    std::cout
        << "maximum_e="
        << maximum_e
        << "\n";
}

void run_no_ascent_phase() {
    std::cout
        << "\nPHASE 2: NO-ASCENT CASES\n";

    u64 total = 0;
    u64 multi = 0;
    u64 e1 = 0;
    u64 e2_or_more = 0;

    for (u64 p = 2; p <= 31; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 10000; ++m) {
            const auto digits =
                base_p_digits(m, p);

            const int ascent =
                first_ascent_position(digits);

            if (ascent >= 0) {
                continue;
            }

            ++total;

            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++multi;

            const u64 g =
                start_spacing_gcd(intervals);

            const PowerInfo power =
                analyze_power(g, p);

            if (power.exponent == 1) {
                ++e1;
            } else {
                ++e2_or_more;

                if (e2_or_more <= 30) {
                    print_case(
                        p,
                        m,
                        digits,
                        intervals,
                        power.exponent,
                        1,
                        ascent
                    );
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 2\n";

    std::cout
        << "no_ascent_total="
        << total
        << "\n";

    std::cout
        << "no_ascent_multi_interval="
        << multi
        << "\n";

    std::cout
        << "no_ascent_e1="
        << e1
        << "\n";

    std::cout
        << "no_ascent_e>=2="
        << e2_or_more
        << "\n";
}

void run_base_specific_phase() {
    std::cout
        << "\nPHASE 3: BASE-SPECIFIC CHECKS\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11, 13, 17, 19
    };

    for (u64 p : primes) {
        u64 total = 0;
        u64 passed = 0;

        std::cout
            << "\np=" << p
            << "\n";

        for (u64 m = 1; m <= 20000; ++m) {
            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const auto digits =
                base_p_digits(m, p);

            const int ascent =
                first_ascent_position(digits);

            if (ascent < 0) {
                /*
                    Skip the unresolved no-ascent class
                    in this phase.
                */
                continue;
            }

            const u64 g =
                start_spacing_gcd(intervals);

            const PowerInfo power =
                analyze_power(g, p);

            const int predicted =
                ascent + 1;

            if (power.exponent == predicted &&
                power.residual == 1) {

                ++passed;
            }
        }

        std::cout
            << "multi_interval_cases="
            << total
            << "\n";

        std::cout
            << "first_ascent_pass="
            << passed
            << "\n";
    }
}

void run_targeted_examples() {
    std::cout
        << "\nPHASE 4: TARGETED EXAMPLES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 6},
        {2, 12},
        {2, 13},
        {2, 20},
        {2, 21},
        {2, 24},
        {2, 25},
        {2, 27},
        {2, 28},
        {2, 40},
        {2, 48},
        {2, 96},
        {2, 192},
        {2, 384},
        {3, 18},
        {3, 23},
        {3, 36},
        {3, 54},
        {3, 59},
        {5, 50},
        {5, 59},
        {5, 75},
        {5, 94},
        {7, 98},
        {11, 242},
        {17, 578}
    };

    for (const auto& c : cases) {
        const auto intervals =
            build_intervals(c.p, c.m);

        if (intervals.size() < 2) {
            continue;
        }

        const auto digits =
            base_p_digits(c.m, c.p);

        const int ascent =
            first_ascent_position(digits);

        const u64 g =
            start_spacing_gcd(intervals);

        const PowerInfo power =
            analyze_power(g, c.p);

        const int predicted =
            (ascent >= 0)
                ? ascent + 1
                : 1;

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " digits="
            << digits_string(digits)
            << " ascent=";

        if (ascent < 0) {
            std::cout << "NONE";
        } else {
            std::cout << ascent;
        }

        std::cout
            << " observed_e="
            << power.exponent
            << " predicted_e="
            << predicted
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 165;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_phase();
    run_no_ascent_phase();
    run_base_specific_phase();
    run_targeted_examples();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
