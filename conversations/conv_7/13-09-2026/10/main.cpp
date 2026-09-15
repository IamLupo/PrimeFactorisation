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

bool lucas_hit(
    u64 p,
    u64 m,
    u64 t
) {
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

std::vector<Interval> build_intervals(
    u64 p,
    u64 m
) {
    std::vector<Interval> intervals;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit =
            lucas_hit(p, m, t);

        if (hit && !inside) {
            start = t;
            inside = true;
        }

        if (!hit && inside) {
            intervals.push_back({
                start,
                t - 1
            });

            inside = false;
        }
    }

    if (inside) {
        intervals.push_back({
            start,
            m
        });
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

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        const u64 diff =
            intervals[i].lo -
            intervals[i - 1].lo;

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

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        const u64 diff =
            intervals[i].hi -
            intervals[i - 1].hi;

        g = std::gcd(g, diff);
    }

    return g;
}

PowerInfo analyze_power(
    u64 value,
    u64 p
) {
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
    std::string result = "[";

    for (std::size_t i = 0;
         i < digits.size();
         ++i) {

        if (i != 0) {
            result += ",";
        }

        result +=
            std::to_string(digits[i]);
    }

    result += "]";

    return result;
}

/*
    NEW CONJECTURED RULE

    Digits are least-significant first.

    Repeatedly remove low-order digits equal to p-1.

    At the first remaining digit d < p-1:

        e = 1 + number of consecutive zero digits
            immediately following d.

    Examples:

        p=3, digits [1,1,2]
        first digit is 1 < 2
        zero-run = 0
        e = 1

        p=3, digits [2,1,0,2]
        remove 2
        remaining [1,0,2]
        zero-run = 1
        e = 2

    However, the data showed m=59 has e=3.
    Therefore the recursive interpretation is:

        If d0 = p-1:
            e = 1 + e(high)

        Otherwise:
            e = 1 + zero_run(high)

    This function implements exactly that.
*/
int predicted_exponent(
    u64 m,
    u64 p
) {
    int accumulated = 0;

    while (m > 0) {
        const u64 digit = m % p;
        m /= p;

        if (digit == p - 1) {
            ++accumulated;
            continue;
        }

        ++accumulated;

        while (m > 0 &&
               m % p == 0) {

            ++accumulated;
            m /= p;
        }

        return accumulated;
    }

    /*
        All digits were p-1.

        Such m have no remaining digit below p-1.
        In the observed data these are generally
        zero- or one-interval cases.

        Return accumulated as the formal extension.
    */
    return accumulated;
}

int first_non_pminus1_position(
    u64 m,
    u64 p
) {
    int position = 0;

    while (m > 0) {
        const u64 digit = m % p;

        if (digit != p - 1) {
            return position;
        }

        m /= p;
        ++position;
    }

    return -1;
}

int zero_run_after_position(
    u64 m,
    u64 p,
    int position
) {
    if (position < 0) {
        return 0;
    }

    for (int i = 0; i <= position; ++i) {
        m /= p;
    }

    int zeros = 0;

    while (m > 0 &&
           m % p == 0) {

        ++zeros;
        m /= p;
    }

    return zeros;
}

void print_failure(
    u64 p,
    u64 m,
    const std::vector<int>& digits,
    int observed,
    int predicted,
    const std::vector<Interval>& intervals
) {
    std::cout
        << "\nFAILURE\n";

    std::cout
        << "p=" << p
        << " m=" << m
        << " digits="
        << digits_string(digits)
        << "\n";

    std::cout
        << "observed_e="
        << observed
        << " predicted_e="
        << predicted
        << "\n";

    const int first =
        first_non_pminus1_position(m, p);

    std::cout
        << "first_non_pminus1_position="
        << first
        << "\n";

    std::cout
        << "zero_run_after="
        << zero_run_after_position(
            m,
            p,
            first
        )
        << "\n";

    std::cout
        << "intervals="
        << intervals.size()
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

        std::cout
            << intervals[i].lo;
    }

    if (intervals.size() > limit) {
        std::cout << ",...";
    }

    std::cout << "\n";
}

void run_exhaustive_phase() {
    std::cout
        << "\nPHASE 1: EXHAUSTIVE RECURSIVE RULE\n";

    u64 total = 0;
    u64 passed = 0;
    u64 failed = 0;

    u64 pure_power_fail = 0;
    u64 endpoint_mismatch = 0;

    int maximum_observed = 0;

    for (u64 p = 2; p <= 149; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1;
             m <= 5000;
             ++m) {

            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const u64 gs =
                start_spacing_gcd(intervals);

            const u64 ge =
                end_spacing_gcd(intervals);

            const PowerInfo ps =
                analyze_power(gs, p);

            const PowerInfo pe =
                analyze_power(ge, p);

            maximum_observed =
                std::max(
                    maximum_observed,
                    ps.exponent
                );

            if (ps.residual != 1 ||
                pe.residual != 1) {

                ++pure_power_fail;
            }

            if (ps.exponent !=
                pe.exponent) {

                ++endpoint_mismatch;
            }

            const int predicted =
                predicted_exponent(m, p);

            if (ps.residual == 1 &&
                pe.residual == 1 &&
                ps.exponent ==
                    pe.exponent &&
                ps.exponent ==
                    predicted) {

                ++passed;
            } else {
                ++failed;

                if (failed <= 40) {
                    print_failure(
                        p,
                        m,
                        base_p_digits(m, p),
                        ps.exponent,
                        predicted,
                        intervals
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
        << passed
        << "/" << total
        << "\n";

    std::cout
        << "prediction_fail="
        << failed
        << "\n";

    std::cout
        << "pure_power_fail="
        << pure_power_fail
        << "\n";

    std::cout
        << "start_end_exponent_mismatch="
        << endpoint_mismatch
        << "\n";

    std::cout
        << "maximum_observed_e="
        << maximum_observed
        << "\n";
}

void run_large_phase() {
    std::cout
        << "\nPHASE 2: LARGER CASES\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        17, 31, 101, 211,
        431, 1009, 2003
    };

    u64 total = 0;
    u64 passed = 0;

    for (u64 p : primes) {

        for (u64 m = 1;
             m <= 50000;
             ++m) {

            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const u64 g =
                start_spacing_gcd(intervals);

            const PowerInfo observed =
                analyze_power(g, p);

            const int predicted =
                predicted_exponent(m, p);

            if (observed.residual == 1 &&
                observed.exponent ==
                    predicted) {

                ++passed;
            } else if (passed < 10) {
                print_failure(
                    p,
                    m,
                    base_p_digits(m, p),
                    observed.exponent,
                    predicted,
                    intervals
                );
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 2\n";

    std::cout
        << "large_multi_interval_cases="
        << total
        << "\n";

    std::cout
        << "large_prediction_pass="
        << passed
        << "/" << total
        << "\n";
}

void run_digit_family_phase() {
    std::cout
        << "\nPHASE 3: CONTROLLED DIGIT FAMILIES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    /*
        These families deliberately exercise:

        1. long p-1 prefixes
        2. long zero runs
        3. p-1 followed by zero runs
        4. middle digits
        5. combinations of all three
    */
    std::vector<TestCase> cases = {
        {3, 41},
        {3, 122},
        {3, 125},
        {3, 365},
        {3, 368},
        {3, 377},
        {3, 1094},
        {3, 1097},
        {3, 1106},
        {3, 1133},

        {5, 64},
        {5, 69},
        {5, 94},
        {5, 159},
        {5, 164},
        {5, 169},
        {5, 189},
        {5, 194},
        {5, 219},

        {7, 98},
        {7, 139},
        {11, 242},

        {2, 12},
        {2, 20},
        {2, 24},
        {2, 28},
        {2, 40},
        {2, 48},
        {2, 96},
        {2, 192},
        {2, 384},
        {2, 768},
        {2, 1536}
    };

    for (const auto& c : cases) {
        const auto intervals =
            build_intervals(c.p, c.m);

        if (intervals.size() < 2) {
            std::cout
                << "p=" << c.p
                << " m=" << c.m
                << " insufficient_intervals\n";
            continue;
        }

        const auto digits =
            base_p_digits(c.m, c.p);

        const u64 g =
            start_spacing_gcd(intervals);

        const PowerInfo observed =
            analyze_power(g, c.p);

        const int predicted =
            predicted_exponent(c.m, c.p);

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " digits="
            << digits_string(digits)
            << " observed_e="
            << observed.exponent
            << " predicted_e="
            << predicted
            << " G="
            << g
            << "\n";
    }
}

void run_recursive_trace_phase() {
    std::cout
        << "\nPHASE 4: RECURSIVE TRACES\n";

    const std::vector<std::pair<u64, u64>> cases = {
        {3, 22},
        {3, 41},
        {3, 59},
        {3, 125},
        {5, 59},
        {5, 94},
        {2, 24},
        {2, 48},
        {2, 96}
    };

    for (const auto& [p, original] : cases) {
        u64 m = original;

        std::cout
            << "\np=" << p
            << " m=" << original
            << "\n";

        int accumulated = 0;

        while (m > 0) {
            const u64 digit =
                m % p;

            m /= p;

            std::cout
                << "digit="
                << digit;

            if (digit == p - 1) {
                ++accumulated;

                std::cout
                    << " -> strip p-1"
                    << " accumulated="
                    << accumulated
                    << "\n";

                continue;
            }

            ++accumulated;

            int zeros = 0;

            while (m > 0 &&
                   m % p == 0) {

                ++zeros;
                ++accumulated;
                m /= p;
            }

            std::cout
                << " -> first digit < p-1"
                << " zero_run="
                << zeros
                << " predicted_e="
                << accumulated
                << "\n";

            break;
        }
    }
}

int main() {
    constexpr int EXPERIMENT = 166;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_phase();
    run_large_phase();
    run_digit_family_phase();
    run_recursive_trace_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
