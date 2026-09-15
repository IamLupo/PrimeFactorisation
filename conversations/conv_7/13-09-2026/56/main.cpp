#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
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

struct CountResult {
    u64 value = 0;
    bool safe = true;
};

u64 ipow_u64(
    u64 p,
    std::size_t e
) {
    u64 result = 1;

    for (std::size_t i = 0; i < e; ++i) {
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

u64 product_digit_plus_one(
    const std::vector<u64>& digits
) {
    u64 product = 1;

    for (u64 d : digits) {
        product =
            safe_multiply(
                product,
                d + 1
            );
    }

    return product;
}

u64 expected_interval_count(
    u64 q,
    u64 p
) {
    const std::vector<u64> digits =
        base_p_digits(q, p);

    return
        product_digit_plus_one(digits)
        - 1;
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

        /*
            left = s0 + j*p^e
        */
        if (!multiplication_safe(
                j,
                powers[e]
            )) {

            result.arithmetic_safe = false;
            return;
        }

        const u64 left_offset =
            j * powers[e];

        if (left_offset >
            std::numeric_limits<u64>::max() -
                s0) {

            result.arithmetic_safe = false;
            return;
        }

        const u64 left =
            s0 + left_offset;

        /*
            right =
                (floor(j/p^r)+1)*p^(e+r)-1
        */
        const u64 quotient =
            j / powers[r];

        if (quotient ==
            std::numeric_limits<u64>::max()) {

            result.arithmetic_safe = false;
            return;
        }

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

        if (!result.arithmetic_safe) {
            return;
        }
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

/*
    Direct digit-count formula:

        MISS in [0,m]
            = product_i(m_i+1)

    HIT in [1,m]:

        H = m+1 - product_i(m_i+1)
*/
u64 direct_digit_hit_count(
    u64 m,
    u64 p
) {
    const std::vector<u64> m_digits =
        base_p_digits(m, p);

    const u64 miss_count =
        product_digit_plus_one(
            m_digits
        );

    return
        m + 1 - miss_count;
}

/*
    First Phase-I closed form:

        H =
            q*p^e
            - s0*(product_i(q_i+1)-1)
*/
u64 closed_form_hit_count(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const std::vector<u64> q_digits =
        base_p_digits(q, p);

    const u64 digit_product =
        product_digit_plus_one(
            q_digits
        );

    const u64 intervals =
        digit_product - 1;

    const u64 p_e =
        ipow_u64(
            p,
            e
        );

    const u64 main_term =
        safe_multiply(
            q,
            p_e
        );

    const u64 correction =
        safe_multiply(
            s0,
            intervals
        );

    if (main_term < correction) {
        throw std::overflow_error(
            "closed-form subtraction underflow"
        );
    }

    return
        main_term -
        correction;
}

/*
    Sum of interval lengths.
*/
u64 summed_interval_length(
    const std::vector<Interval>& intervals
) {
    u64 total = 0;

    for (const Interval& interval :
         intervals) {

        const u64 length =
            interval.right -
            interval.left +
            1;

        total =
            safe_add(
                total,
                length
            );
    }

    return total;
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

struct GapStats {
    u64 overlaps = 0;
    u64 touching = 0;
    u64 proper_gaps = 0;

    u64 min_gap =
        std::numeric_limits<u64>::max();

    u64 max_gap = 0;
};

GapStats analyze_gaps(
    const std::vector<Interval>& intervals
) {
    GapStats stats;

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        const Interval& a =
            intervals[i - 1];

        const Interval& b =
            intervals[i];

        if (b.left <= a.right) {
            ++stats.overlaps;
            continue;
        }

        const u64 gap =
            b.left -
            a.right -
            1;

        if (gap == 0) {
            ++stats.touching;
        } else {
            ++stats.proper_gaps;

            stats.min_gap =
                std::min(
                    stats.min_gap,
                    gap
                );

            stats.max_gap =
                std::max(
                    stats.max_gap,
                    gap
                );
        }
    }

    if (stats.min_gap ==
        std::numeric_limits<u64>::max()) {

        stats.min_gap = 0;
    }

    return stats;
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

    u64 p_e;

    try {
        p_e =
            ipow_u64(
                p,
                e
            );
    } catch (const std::exception& ex) {
        std::cout
            << "arithmetic_safe=0\n"
            << "reason=" << ex.what()
            << "\n";
        return;
    }

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

    const u64 expected_intervals =
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

    const u64 interval_sum =
        summed_interval_length(
            result.intervals
        );

    const u64 digit_formula =
        direct_digit_hit_count(
            m,
            p
        );

    const u64 closed_formula =
        closed_form_hit_count(
            p,
            e,
            s0,
            q
        );

    sort_intervals(
        result.intervals
    );

    const GapStats gaps =
        analyze_gaps(
            result.intervals
        );

    const auto end =
        std::chrono::steady_clock::now();

    const double generation_seconds =
        std::chrono::duration<double>(
            generation_end -
            begin
        ).count();

    const double total_seconds =
        std::chrono::duration<double>(
            end -
            begin
        ).count();

    const bool interval_count_pass =
        result.generated ==
        expected_intervals;

    const bool digit_formula_pass =
        interval_sum ==
        digit_formula;

    const bool closed_formula_pass =
        interval_sum ==
        closed_formula;

    const bool geometry_pass =
        gaps.overlaps == 0 &&
        gaps.touching == 0;

    const bool final_pass =
        result.invariant_ok &&
        interval_count_pass &&
        digit_formula_pass &&
        closed_formula_pass &&
        geometry_pass;

    std::cout
        << "m="
        << m
        << "\n"
        << "expected_intervals="
        << expected_intervals
        << "\n"
        << "generated_intervals="
        << result.generated
        << "\n"
        << "enumeration_pass="
        << (interval_count_pass ? 1 : 0)
        << "\n"
        << "invariant_ok="
        << (result.invariant_ok ? 1 : 0)
        << "\n"
        << "summed_interval_length="
        << interval_sum
        << "\n"
        << "direct_digit_hit_count="
        << digit_formula
        << "\n"
        << "digit_formula_pass="
        << (digit_formula_pass ? 1 : 0)
        << "\n"
        << "closed_form_hit_count="
        << closed_formula
        << "\n"
        << "closed_formula_pass="
        << (closed_formula_pass ? 1 : 0)
        << "\n"
        << "overlaps="
        << gaps.overlaps
        << "\n"
        << "touching="
        << gaps.touching
        << "\n"
        << "proper_gaps="
        << gaps.proper_gaps
        << "\n"
        << "min_gap="
        << gaps.min_gap
        << "\n"
        << "max_gap="
        << gaps.max_gap
        << "\n"
        << "geometry_pass="
        << (geometry_pass ? 1 : 0)
        << "\n"
        << "final_pass="
        << (final_pass ? 1 : 0)
        << "\n"
        << "generation_seconds="
        << std::setprecision(10)
        << generation_seconds
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
}

int main() {
    std::cout
        << "START EXPERIMENT 212\n";

    /*
        Dense binary case.
    */
    run_case(
        2,
        40,
        1,
        1048575,
        true
    );

    /*
        Dense base-5 case that safely fits uint64_t.
    */
    run_case(
        5,
        15,
        1,
        390624,
        true
    );

    /*
        Sparse binary case.
    */
    run_case(
        2,
        20,
        1,
        1073741825ULL,
        true
    );

    /*
        Small exact case with s0 != 1.
        This checks that the correction term really uses s0.
    */
    run_case(
        3,
        5,
        27,
        80,
        true
    );

    /*
        Deliberately overflowing case.
    */
    run_case(
        5,
        20,
        1,
        390624,
        false
    );

    std::cout
        << "\nIDENTITIES\n"
        << "H_interval = sum_j |I_j|\n"
        << "H_digit = m+1-product_i(m_i+1)\n"
        << "H_closed = q*p^e-s0*(product_i(q_i+1)-1)\n"
        << "target: H_interval = H_digit = H_closed\n"
        << "FINISHED EXPERIMENT 212\n";

    return 0;
}
