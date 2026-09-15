#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 left;
    u64 right;
};

struct EnumerationResult {
    std::vector<Interval> intervals;
    u64 generated = 0;
    bool invariant_ok = true;
    bool arithmetic_safe = true;
};

struct GapResult {
    u64 checked = 0;
    u64 failures = 0;

    u64 total_gap = 0;
    u64 expected_total_gap = 0;

    u64 min_gap = std::numeric_limits<u64>::max();
    u64 max_gap = 0;

    bool all_equal_s0 = true;
};

u64 ipow_u64(
    u64 p,
    std::size_t e
) {
    u64 result = 1;

    for (std::size_t i = 0;
         i < e;
         ++i) {

        if (result >
            std::numeric_limits<u64>::max() / p) {

            throw std::overflow_error(
                "power exceeds uint64_t"
            );
        }

        result *= p;
    }

    return result;
}

bool multiplication_safe(
    u64 a,
    u64 b
) {
    return a == 0 ||
           b <=
               std::numeric_limits<u64>::max() / a;
}

u64 safe_multiply(
    u64 a,
    u64 b
) {
    if (!multiplication_safe(a, b)) {
        throw std::overflow_error(
            "multiplication exceeds uint64_t"
        );
    }

    return a * b;
}

u64 safe_add(
    u64 a,
    u64 b
) {
    if (a >
        std::numeric_limits<u64>::max() - b) {

        throw std::overflow_error(
            "addition exceeds uint64_t"
        );
    }

    return a + b;
}

std::vector<u64> base_p_digits(
    u64 n,
    u64 p
) {
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

std::vector<u64> make_powers(
    u64 p,
    std::size_t count
) {
    std::vector<u64> powers(
        count + 1,
        1
    );

    for (std::size_t i = 1;
         i <= count;
         ++i) {

        powers[i] =
            ipow_u64(p, i);
    }

    return powers;
}

std::vector<u64> make_q_digits(
    u64 q,
    u64 p,
    std::size_t count
) {
    std::vector<u64> digits(
        count,
        0
    );

    for (std::size_t i = 0;
         i < count;
         ++i) {

        digits[i] = q % p;
        q /= p;
    }

    return digits;
}

u64 expected_interval_count(
    u64 q,
    u64 p
) {
    const std::vector<u64> digits =
        base_p_digits(q, p);

    u64 product = 1;

    for (u64 d : digits) {
        product =
            safe_multiply(
                product,
                d + 1
            );
    }

    return product - 1;
}

u64 expected_type_count(
    const std::vector<u64>& q_digits,
    std::size_t r
) {
    u64 count =
        q_digits[r];

    for (std::size_t i = r + 1;
         i < q_digits.size();
         ++i) {

        count =
            safe_multiply(
                count,
                q_digits[i] + 1
            );
    }

    return count;
}

bool verify_type(
    u64 j,
    std::size_t r,
    u64 p,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers
) {
    const std::size_t n =
        q_digits.size();

    for (std::size_t i = 0;
         i < n;
         ++i) {

        const u64 digit =
            (j / powers[i]) % p;

        if (digit > q_digits[i]) {
            return false;
        }

        if (i < r &&
            digit != q_digits[i]) {
            return false;
        }

        if (i == r &&
            digit >= q_digits[r]) {
            return false;
        }
    }

    return true;
}

void enumerate_high_digits(
    std::size_t r,
    u64 s0,
    u64 e,
    u64 p,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    std::size_t position,
    u64 current_j,
    EnumerationResult& result
) {
    const std::size_t n =
        q_digits.size();

    if (!result.arithmetic_safe) {
        return;
    }

    if (position >= n) {
        const u64 j =
            current_j;

        if (!verify_type(
                j,
                r,
                p,
                q_digits,
                powers
            )) {

            result.invariant_ok = false;
            return;
        }

        if (!multiplication_safe(
                j,
                powers[e]
            )) {

            result.arithmetic_safe = false;
            return;
        }

        const u64 left =
            safe_add(
                s0,
                j * powers[e]
            );

        const u64 quotient =
            j / powers[r];

        const u64 factor =
            quotient + 1;

        if (!multiplication_safe(
                factor,
                powers[e + r]
            )) {

            result.arithmetic_safe = false;
            return;
        }

        const u64 right =
            factor *
            powers[e + r] -
            1;

        result.intervals.push_back({
            left,
            right
        });

        ++result.generated;
        return;
    }

    for (u64 d = 0;
         d <= q_digits[position];
         ++d) {

        if (!multiplication_safe(
                d,
                powers[position]
            )) {

            result.arithmetic_safe = false;
            return;
        }

        const u64 contribution =
            d * powers[position];

        if (current_j >
            std::numeric_limits<u64>::max() -
                contribution) {

            result.arithmetic_safe = false;
            return;
        }

        const u64 next_j =
            current_j +
            contribution;

        enumerate_high_digits(
            r,
            s0,
            e,
            p,
            q_digits,
            powers,
            position + 1,
            next_j,
            result
        );
    }
}

EnumerationResult enumerate_by_type(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    EnumerationResult result;

    const std::vector<u64> raw_q =
        base_p_digits(q, p);

    const std::size_t digit_count =
        raw_q.size();

    const std::vector<u64> q_digits =
        make_q_digits(
            q,
            p,
            digit_count
        );

    const std::vector<u64> powers =
        make_powers(
            p,
            e + digit_count + 1
        );

    const u64 expected =
        expected_interval_count(
            q,
            p
        );

    result.intervals.reserve(
        static_cast<std::size_t>(
            expected
        )
    );

    for (std::size_t r = 0;
         r < digit_count;
         ++r) {

        if (q_digits[r] == 0) {
            continue;
        }

        const u64 expected_type =
            expected_type_count(
                q_digits,
                r
            );

        const u64 before =
            result.generated;

        u64 fixed_low = 0;

        for (std::size_t i = 0;
             i < r;
             ++i) {

            const u64 contribution =
                safe_multiply(
                    q_digits[i],
                    powers[i]
                );

            fixed_low =
                safe_add(
                    fixed_low,
                    contribution
                );
        }

        for (u64 d = 0;
             d < q_digits[r];
             ++d) {

            const u64 contribution =
                safe_multiply(
                    d,
                    powers[r]
                );

            const u64 base =
                safe_add(
                    fixed_low,
                    contribution
                );

            enumerate_high_digits(
                r,
                s0,
                e,
                p,
                q_digits,
                powers,
                r + 1,
                base,
                result
            );

            if (!result.arithmetic_safe) {
                return result;
            }
        }

        const u64 generated_type =
            result.generated -
            before;

        if (generated_type !=
            expected_type) {

            result.invariant_ok = false;
        }
    }

    if (result.generated != expected) {
        result.invariant_ok = false;
    }

    return result;
}

void sort_intervals(
    std::vector<Interval>& intervals
) {
    std::sort(
        intervals.begin(),
        intervals.end(),
        [](const Interval& a,
           const Interval& b) {

            if (a.left != b.left) {
                return a.left < b.left;
            }

            return a.right < b.right;
        }
    );
}

GapResult analyze_gaps(
    const std::vector<Interval>& intervals,
    u64 s0
) {
    GapResult result;

    if (intervals.size() < 2) {
        return result;
    }

    result.expected_total_gap =
        safe_multiply(
            s0,
            intervals.size() - 1
        );

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        const Interval& previous =
            intervals[i - 1];

        const Interval& current =
            intervals[i];

        ++result.checked;

        if (current.left <=
            previous.right) {

            ++result.failures;
            result.all_equal_s0 = false;
            continue;
        }

        const u64 gap =
            current.left -
            previous.right -
            1;

        result.total_gap =
            safe_add(
                result.total_gap,
                gap
            );

        result.min_gap =
            std::min(
                result.min_gap,
                gap
            );

        result.max_gap =
            std::max(
                result.max_gap,
                gap
            );

        if (gap != s0) {
            ++result.failures;
            result.all_equal_s0 = false;
        }
    }

    return result;
}

void print_sample(
    const std::vector<Interval>& intervals,
    std::size_t count
) {
    const std::size_t n =
        std::min(
            count,
            intervals.size()
        );

    for (std::size_t i = 0;
         i < n;
         ++i) {

        std::cout
            << "sample[" << i << "]=["
            << intervals[i].left
            << ","
            << intervals[i].right
            << "]\n";
    }
}

void run_case(
    u64 p,
    u64 e,
    u64 s0,
    u64 q,
    bool show_samples
) {
    std::cout
        << "\nCASE\n"
        << "p=" << p
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
        << "\n";

    try {
        const u64 p_e =
            ipow_u64(
                p,
                e
            );

        if (!multiplication_safe(
                q,
                p_e
            )) {

            std::cout
                << "arithmetic_safe=0\n"
                << "case_skipped_due_to_uint64_overflow=1\n";
            return;
        }

        const u64 q_times_p_e =
            q * p_e;

        if (q_times_p_e >
            std::numeric_limits<u64>::max() -
                s0) {

            std::cout
                << "arithmetic_safe=0\n"
                << "case_skipped_due_to_uint64_overflow=1\n";
            return;
        }

        const u64 m =
            q_times_p_e +
            s0 -
            1;

        const u64 expected =
            expected_interval_count(
                q,
                p
            );

        const auto begin =
            std::chrono::steady_clock::now();

        EnumerationResult result =
            enumerate_by_type(
                p,
                e,
                s0,
                q
            );

        const auto generation_end =
            std::chrono::steady_clock::now();

        if (!result.arithmetic_safe) {
            std::cout
                << "arithmetic_safe=0\n"
                << "enumeration_overflow=1\n";
            return;
        }

        sort_intervals(
            result.intervals
        );

        const auto sort_end =
            std::chrono::steady_clock::now();

        const GapResult gaps =
            analyze_gaps(
                result.intervals,
                s0
            );

        const auto end =
            std::chrono::steady_clock::now();

        const bool count_pass =
            result.generated == expected;

        const bool gap_pass =
            gaps.failures == 0 &&
            gaps.all_equal_s0;

        const bool sum_pass =
            gaps.total_gap ==
            gaps.expected_total_gap;

        const bool final_pass =
            result.invariant_ok &&
            count_pass &&
            gap_pass &&
            sum_pass;

        const double generation_seconds =
            std::chrono::duration<double>(
                generation_end -
                begin
            ).count();

        const double sort_seconds =
            std::chrono::duration<double>(
                sort_end -
                generation_end
            ).count();

        const double total_seconds =
            std::chrono::duration<double>(
                end -
                begin
            ).count();

        std::cout
            << "m="
            << m
            << "\n"
            << "expected_intervals="
            << expected
            << "\n"
            << "generated_intervals="
            << result.generated
            << "\n"
            << "enumeration_pass="
            << (count_pass ? 1 : 0)
            << "\n"
            << "invariant_ok="
            << (result.invariant_ok ? 1 : 0)
            << "\n"
            << "gap_checks="
            << gaps.checked
            << "\n"
            << "gap_failures="
            << gaps.failures
            << "\n"
            << "all_internal_gaps_equal_s0="
            << (gaps.all_equal_s0 ? 1 : 0)
            << "\n"
            << "total_internal_gap="
            << gaps.total_gap
            << "\n"
            << "expected_internal_gap="
            << gaps.expected_total_gap
            << "\n"
            << "gap_sum_pass="
            << (sum_pass ? 1 : 0)
            << "\n"
            << "min_gap="
            << gaps.min_gap
            << "\n"
            << "max_gap="
            << gaps.max_gap
            << "\n"
            << "gap_theorem_pass="
            << (gap_pass ? 1 : 0)
            << "\n"
            << "final_pass="
            << (final_pass ? 1 : 0)
            << "\n"
            << "generation_seconds="
            << std::setprecision(10)
            << generation_seconds
            << "\n"
            << "sort_seconds="
            << sort_seconds
            << "\n"
            << "total_seconds="
            << total_seconds
            << "\n";

        if (show_samples) {
            print_sample(
                result.intervals,
                10
            );
        }
    } catch (const std::exception& ex) {
        std::cout
            << "arithmetic_safe=0\n"
            << "reason="
            << ex.what()
            << "\n";
    }
}

u64 random_q(
    std::mt19937_64& rng,
    u64 p,
    std::size_t digits
) {
    std::uniform_int_distribution<u64> digit_dist(
        0,
        p - 1
    );

    u64 q = 0;
    u64 power = 1;

    for (std::size_t i = 0;
         i < digits;
         ++i) {

        const u64 d =
            digit_dist(rng);

        q =
            safe_add(
                q,
                safe_multiply(
                    d,
                    power
                )
            );

        power =
            safe_multiply(
                power,
                p
            );
    }

    if (q == 0) {
        q = 1;
    }

    return q;
}

void run_random_gap_suite(
    u64 p,
    std::size_t digits,
    u64 e,
    std::size_t cases,
    u64 seed
) {
    std::mt19937_64 rng(seed);

    u64 failures = 0;
    u64 overflowed = 0;
    u64 total_intervals = 0;

    u64 worst_min_gap =
        std::numeric_limits<u64>::max();

    u64 worst_max_gap = 0;

    for (std::size_t case_index = 0;
         case_index < cases;
         ++case_index) {

        const u64 q =
            random_q(
                rng,
                p,
                digits
            );

        /*
            Use s0 values including nontrivial values.
        */
        const u64 s0 =
            1 +
            (case_index % 9);

        try {
            const u64 p_e =
                ipow_u64(p, e);

            if (!multiplication_safe(
                    q,
                    p_e
                )) {
                ++overflowed;
                continue;
            }

            if (q * p_e >
                std::numeric_limits<u64>::max() -
                    s0) {

                ++overflowed;
                continue;
            }

            EnumerationResult result =
                enumerate_by_type(
                    p,
                    e,
                    s0,
                    q
                );

            if (!result.arithmetic_safe) {
                ++overflowed;
                continue;
            }

            sort_intervals(
                result.intervals
            );

            const GapResult gaps =
                analyze_gaps(
                    result.intervals,
                    s0
                );

            const u64 expected =
                expected_interval_count(
                    q,
                    p
                );

            const bool pass =
                result.invariant_ok &&
                result.generated == expected &&
                gaps.failures == 0 &&
                gaps.total_gap ==
                    gaps.expected_total_gap;

            if (!pass) {
                ++failures;
            }

            total_intervals +=
                result.generated;

            if (!result.intervals.empty()) {
                worst_min_gap =
                    std::min(
                        worst_min_gap,
                        gaps.min_gap
                    );

                worst_max_gap =
                    std::max(
                        worst_max_gap,
                        gaps.max_gap
                    );
            }
        } catch (...) {
            ++overflowed;
        }
    }

    if (worst_min_gap ==
        std::numeric_limits<u64>::max()) {

        worst_min_gap = 0;
    }

    std::cout
        << "\nRANDOM GAP SUITE\n"
        << "p=" << p
        << " digits=" << digits
        << " e=" << e
        << " cases=" << cases
        << "\n"
        << "failures="
        << failures
        << "\n"
        << "overflowed_or_skipped="
        << overflowed
        << "\n"
        << "total_intervals="
        << total_intervals
        << "\n"
        << "observed_min_gap="
        << worst_min_gap
        << "\n"
        << "observed_max_gap="
        << worst_max_gap
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 213\n";

    /*
        Exact earlier cases.
    */
    run_case(
        2,
        3,
        1,
        63,
        true
    );

    run_case(
        2,
        4,
        1,
        255,
        true
    );

    run_case(
        3,
        4,
        1,
        80,
        true
    );

    run_case(
        5,
        3,
        1,
        100,
        true
    );

    /*
        Nontrivial s0.
    */
    run_case(
        3,
        5,
        27,
        80,
        true
    );

    /*
        Large safe cases.
    */
    run_case(
        2,
        40,
        1,
        1048575,
        false
    );

    run_case(
        5,
        15,
        1,
        390624,
        false
    );

    run_case(
        2,
        20,
        1,
        1073741825ULL,
        false
    );

    /*
        Random verification with varying s0.
    */
    run_random_gap_suite(
        2,
        16,
        20,
        500,
        0x12345678ULL
    );

    run_random_gap_suite(
        3,
        10,
        10,
        500,
        0x87654321ULL
    );

    run_random_gap_suite(
        5,
        8,
        10,
        500,
        0xCAFEBABEULL
    );

    std::cout
        << "\nIDENTITY\n"
        << "For consecutive maximal HIT intervals I_i,I_{i+1}:\n"
        << "gap_i = left(I_{i+1})-right(I_i)-1\n"
        << "target: gap_i = s0\n"
        << "therefore total internal MISS length = s0*(I-1)\n";

    std::cout
        << "FINISHED EXPERIMENT 213\n";

    return 0;
}
