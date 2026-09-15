#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 left;
    u64 right;
};

struct EnumerationResult {
    std::vector<Interval> intervals;
    u64 increment_steps = 0;
    u64 carry_steps = 0;
    u64 reset_steps = 0;
};

u64 ipow_u64(u64 p, std::size_t e) {
    u64 result = 1;

    for (std::size_t i = 0; i < e; ++i) {
        result *= p;
    }

    return result;
}

std::vector<u64> base_p_digits(u64 n, u64 p) {
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

u64 expected_interval_count(u64 q, u64 p) {
    const std::vector<u64> q_digits = base_p_digits(q, p);

    u64 product = 1;

    for (u64 d : q_digits) {
        product *= (d + 1);
    }

    return product - 1;
}

std::vector<u64> make_powers(u64 p, std::size_t count) {
    std::vector<u64> powers(count + 1, 1);

    for (std::size_t i = 1; i <= count; ++i) {
        powers[i] = powers[i - 1] * p;
    }

    return powers;
}

std::vector<u64> make_q_digits(
    u64 q,
    u64 p,
    std::size_t count
) {
    std::vector<u64> digits(count, 0);

    for (std::size_t i = 0; i < count; ++i) {
        digits[i] = q % p;
        q /= p;
    }

    return digits;
}

std::size_t lowest_positive_q_digit(
    const std::vector<u64>& q_digits
) {
    for (std::size_t i = 0; i < q_digits.size(); ++i) {
        if (q_digits[i] > 0) {
            return i;
        }
    }

    return q_digits.size();
}

std::size_t find_r(
    const std::vector<u64>& digits,
    const std::vector<u64>& q_digits
) {
    const std::size_t n = q_digits.size();

    for (std::size_t i = 0; i < n; ++i) {
        if (digits[i] < q_digits[i]) {
            return i;
        }
    }

    return n;
}

u64 digits_to_value(
    const std::vector<u64>& digits,
    const std::vector<u64>& powers
) {
    u64 value = 0;

    for (std::size_t i = 0; i < digits.size(); ++i) {
        value += digits[i] * powers[i];
    }

    return value;
}

bool increment_mixed_radix(
    std::vector<u64>& digits,
    const std::vector<u64>& q_digits,
    EnumerationResult& stats
) {
    const std::size_t n = digits.size();

    ++stats.increment_steps;

    for (std::size_t i = 0; i < n; ++i) {
        if (digits[i] < q_digits[i]) {
            ++digits[i];
            return true;
        }

        digits[i] = 0;
        ++stats.carry_steps;
        ++stats.reset_steps;
    }

    return false;
}

EnumerationResult enumerate_mixed_radix(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    EnumerationResult result;

    const std::vector<u64> q_raw = base_p_digits(q, p);

    const std::size_t digit_count = q_raw.size();

    const std::vector<u64> q_digits =
        make_q_digits(q, p, digit_count);

    const std::vector<u64> powers =
        make_powers(p, e + digit_count + 1);

    const u64 p_e = powers[e];

    const u64 expected =
        expected_interval_count(q, p);

    result.intervals.reserve(
        static_cast<std::size_t>(expected)
    );

    std::vector<u64> digits(digit_count, 0);

    while (true) {
        const u64 j =
            digits_to_value(digits, powers);

        if (j == q) {
            break;
        }

        const std::size_t r =
            find_r(digits, q_digits);

        const u64 left =
            s0 + j * p_e;

        const u64 right =
            (j / powers[r] + 1) *
            powers[e + r] - 1;

        result.intervals.push_back({
            left,
            right
        });

        if (!increment_mixed_radix(
                digits,
                q_digits,
                result
        )) {
            break;
        }
    }

    return result;
}

std::vector<Interval> direct_hit_intervals(
    u64 m,
    u64 p
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        u64 a = t;
        u64 b = m;

        bool miss = true;

        while (a > 0 || b > 0) {
            const u64 da = a % p;
            const u64 db = b % p;

            if (da > db) {
                miss = false;
                break;
            }

            a /= p;
            b /= p;
        }

        const bool hit = !miss;

        if (hit && !inside) {
            start = t;
            inside = true;
        }

        if (!hit && inside) {
            result.push_back({
                start,
                t - 1
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

bool intervals_equal(
    const std::vector<Interval>& a,
    const std::vector<Interval>& b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].left != b[i].left ||
            a[i].right != b[i].right) {
            return false;
        }
    }

    return true;
}

void print_sample(
    const std::vector<Interval>& intervals,
    std::size_t count
) {
    const std::size_t n =
        std::min(count, intervals.size());

    for (std::size_t i = 0; i < n; ++i) {
        std::cout
            << "sample[" << i << "]=["
            << intervals[i].left << ","
            << intervals[i].right << "]\n";
    }
}

void run_small_case(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const u64 m =
        s0 + q * ipow_u64(p, e) - 1;

    const u64 expected =
        expected_interval_count(q, p);

    const auto direct_begin =
        std::chrono::steady_clock::now();

    const std::vector<Interval> direct =
        direct_hit_intervals(m, p);

    const auto direct_end =
        std::chrono::steady_clock::now();

    const auto structural_begin =
        std::chrono::steady_clock::now();

    const EnumerationResult structural =
        enumerate_mixed_radix(p, e, s0, q);

    const auto structural_end =
        std::chrono::steady_clock::now();

    const double direct_seconds =
        std::chrono::duration<double>(
            direct_end - direct_begin
        ).count();

    const double structural_seconds =
        std::chrono::duration<double>(
            structural_end - structural_begin
        ).count();

    std::cout << "\nSMALL CASE\n";

    std::cout
        << "p=" << p
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
        << " m=" << m
        << "\n";

    std::cout
        << "expected_intervals="
        << expected
        << "\n"
        << "direct_intervals="
        << direct.size()
        << "\n"
        << "odometer_intervals="
        << structural.intervals.size()
        << "\n"
        << "set_equality="
        << (intervals_equal(
                direct,
                structural.intervals
            ) ? 1 : 0)
        << "\n"
        << "direct_seconds="
        << std::setprecision(10)
        << direct_seconds
        << "\n"
        << "odometer_seconds="
        << structural_seconds
        << "\n"
        << "increment_steps="
        << structural.increment_steps
        << "\n"
        << "carry_steps="
        << structural.carry_steps
        << "\n"
        << "reset_steps="
        << structural.reset_steps
        << "\n";
}

void run_large_case(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const u64 m =
        s0 + q * ipow_u64(p, e) - 1;

    const u64 expected =
        expected_interval_count(q, p);

    const auto begin =
        std::chrono::steady_clock::now();

    const EnumerationResult result =
        enumerate_mixed_radix(p, e, s0, q);

    const auto end =
        std::chrono::steady_clock::now();

    const double seconds =
        std::chrono::duration<double>(
            end - begin
        ).count();

    const bool pass =
        result.intervals.size() == expected;

    std::cout << "\nLARGE CASE\n";

    std::cout
        << "p=" << p
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
        << "\n"
        << "m=" << m
        << "\n"
        << "expected_intervals="
        << expected
        << "\n"
        << "generated_intervals="
        << result.intervals.size()
        << "\n"
        << "enumeration_pass="
        << (pass ? 1 : 0)
        << "\n"
        << "seconds="
        << std::setprecision(10)
        << seconds
        << "\n";

    if (seconds > 0.0) {
        std::cout
            << "intervals_per_second="
            << static_cast<double>(
                result.intervals.size()
            ) / seconds
            << "\n";
    }

    std::cout
        << "increment_steps="
        << result.increment_steps
        << "\n"
        << "carry_steps="
        << result.carry_steps
        << "\n"
        << "reset_steps="
        << result.reset_steps
        << "\n"
        << "sample_count="
        << std::min<std::size_t>(
            10,
            result.intervals.size()
        )
        << "\n";

    print_sample(result.intervals, 10);

    std::cout
        << "direct_scan_points_would_be="
        << m
        << "\n";
}

int main() {
    std::cout << "START EXPERIMENT 197\n";

    std::size_t small_pass = 0;
    const std::size_t small_cases = 4;

    const std::vector<std::vector<u64>> cases = {
        {2, 3, 1, 63},
        {2, 4, 1, 255},
        {3, 4, 1, 80},
        {5, 3, 1, 100}
    };

    for (const auto& c : cases) {
        const u64 p = c[0];
        const u64 e = c[1];
        const u64 s0 = c[2];
        const u64 q = c[3];

        const u64 m =
            s0 + q * ipow_u64(p, e) - 1;

        const std::vector<Interval> direct =
            direct_hit_intervals(m, p);

        const EnumerationResult structural =
            enumerate_mixed_radix(p, e, s0, q);

        if (intervals_equal(
                direct,
                structural.intervals
            )) {
            ++small_pass;
        }

        run_small_case(
            p,
            e,
            s0,
            q
        );
    }

    std::cout << "\nSMALL SUMMARY\n";

    std::cout
        << "cases="
        << small_cases
        << "\n"
        << "set_equality_pass="
        << small_pass
        << "/"
        << small_cases
        << "\n";

    run_large_case(
        2,
        40,
        1,
        1048575
    );

    run_large_case(
        5,
        20,
        1,
        390624
    );

    run_large_case(
        2,
        20,
        1,
        1073741825ULL
    );

    std::cout << "\nCOMPLEXITY\n";

    std::cout
        << "direct_scan = O(m log_p(m))\n"
        << "structural_recursion = O(I log_p(q))\n"
        << "mixed_radix_enumeration = O(I + carry_work)\n"
        << "I = product_i(q_i+1)-1\n";

    std::cout << "FINISHED EXPERIMENT 197\n";

    return 0;
}