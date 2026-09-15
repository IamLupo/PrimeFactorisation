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

struct GeometryStats {
    u64 overlaps = 0;
    u64 touching = 0;
    u64 proper_gaps = 0;

    u64 min_gap = std::numeric_limits<u64>::max();
    u64 max_gap = 0;

    u64 total_length = 0;
};

struct EnumerationResult {
    std::vector<Interval> intervals;

    u64 generated = 0;

    bool invariant_ok = true;
    bool arithmetic_safe = true;
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
           b <= std::numeric_limits<u64>::max() / a;
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
        powers[i] = ipow_u64(p, i);
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
        product = safe_multiply(
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
    u64 count = q_digits[r];

    for (std::size_t i = r + 1;
         i < q_digits.size();
         ++i) {
        count = safe_multiply(
            count,
            q_digits[i] + 1
        );
    }

    return count;
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

void analyze_geometry(
    const std::vector<Interval>& intervals,
    GeometryStats& stats
) {
    if (intervals.empty()) {
        stats.min_gap = 0;
        return;
    }

    for (const Interval& interval : intervals) {
        if (interval.right < interval.left) {
            ++stats.overlaps;
            continue;
        }

        stats.total_length +=
            interval.right -
            interval.left +
            1;
    }

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        const Interval& previous =
            intervals[i - 1];

        const Interval& current =
            intervals[i];

        if (current.left <= previous.right) {
            ++stats.overlaps;
            continue;
        }

        const u64 gap =
            current.left -
            previous.right -
            1;

        if (gap == 0) {
            ++stats.touching;
        } else {
            ++stats.proper_gaps;
        }

        stats.min_gap =
            std::min(stats.min_gap, gap);

        stats.max_gap =
            std::max(stats.max_gap, gap);
    }

    if (stats.min_gap ==
        std::numeric_limits<u64>::max()) {
        stats.min_gap = 0;
    }
}

bool verify_type(
    u64 j,
    std::size_t r,
    u64 p,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers
) {
    const std::size_t n = q_digits.size();

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
    const std::size_t n = q_digits.size();

    if (!result.arithmetic_safe) {
        return;
    }

    if (position >= n) {
        const u64 j = current_j;

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
        u64 left_offset;

        if (!multiplication_safe(
                j,
                powers[e]
            )) {
            result.arithmetic_safe = false;
            return;
        }

        left_offset =
            j * powers[e];

        if (left_offset >
            std::numeric_limits<u64>::max() - s0) {
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
            d *
            powers[position];

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

            if (!multiplication_safe(
                    q_digits[i],
                    powers[i]
                )) {
                result.arithmetic_safe = false;
                return result;
            }

            fixed_low +=
                q_digits[i] *
                powers[i];
        }

        for (u64 d = 0;
             d < q_digits[r];
             ++d) {

            if (!multiplication_safe(
                    d,
                    powers[r]
                )) {
                result.arithmetic_safe = false;
                return result;
            }

            const u64 contribution =
                d * powers[r];

            if (fixed_low >
                std::numeric_limits<u64>::max() -
                    contribution) {
                result.arithmetic_safe = false;
                return result;
            }

            const u64 base =
                fixed_low +
                contribution;

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

        if (generated_type != expected_type) {
            result.invariant_ok = false;
        }
    }

    const u64 expected_total =
        expected_interval_count(
            q,
            p
        );

    if (result.generated != expected_total) {
        result.invariant_ok = false;
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

    u64 p_e;

    try {
        p_e = ipow_u64(p, e);
    } catch (const std::exception& ex) {
        std::cout
            << "arithmetic_safe=0\n"
            << "reason=" << ex.what()
            << "\n";
        return;
    }

    const bool m_safe =
        multiplication_safe(q, p_e);

    std::cout
        << "p_power_e="
        << p_e
        << "\n"
        << "q_times_p_power_e_safe="
        << (m_safe ? 1 : 0)
        << "\n";

    if (!m_safe) {
        std::cout
            << "case_skipped_due_to_uint64_overflow=1\n";
        return;
    }

    const u64 q_times_p_e =
        q * p_e;

    if (q_times_p_e >
        std::numeric_limits<u64>::max() -
            s0) {
        std::cout
            << "m_addition_safe=0\n"
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

    sort_intervals(
        result.intervals
    );

    const auto sort_end =
        std::chrono::steady_clock::now();

    GeometryStats geometry;

    analyze_geometry(
        result.intervals,
        geometry
    );

    const auto geometry_end =
        std::chrono::steady_clock::now();

    const bool pass =
        result.arithmetic_safe &&
        result.invariant_ok &&
        result.generated == expected &&
        result.intervals.size() == expected &&
        geometry.overlaps == 0 &&
        geometry.touching == 0;

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

    const double geometry_seconds =
        std::chrono::duration<double>(
            geometry_end -
            sort_end
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
        << (pass ? 1 : 0)
        << "\n"
        << "arithmetic_safe="
        << (result.arithmetic_safe ? 1 : 0)
        << "\n"
        << "invariant_ok="
        << (result.invariant_ok ? 1 : 0)
        << "\n"
        << "overlaps="
        << geometry.overlaps
        << "\n"
        << "touching="
        << geometry.touching
        << "\n"
        << "proper_gaps="
        << geometry.proper_gaps
        << "\n"
        << "min_gap="
        << geometry.min_gap
        << "\n"
        << "max_gap="
        << geometry.max_gap
        << "\n"
        << "total_interval_length="
        << geometry.total_length
        << "\n"
        << "generation_seconds="
        << std::setprecision(10)
        << generation_seconds
        << "\n"
        << "sort_seconds="
        << sort_seconds
        << "\n"
        << "geometry_seconds="
        << geometry_seconds
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
        << "START EXPERIMENT 211\n";

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
        Safe dense base-5 case.
    */
    run_case(
        5,
        15,
        1,
        390624,
        true
    );

    /*
        Sparse case.
    */
    run_case(
        2,
        20,
        1,
        1073741825ULL,
        true
    );

    /*
        Deliberately overflowing case.
        This must be detected and skipped.
    */
    run_case(
        5,
        20,
        1,
        390624,
        false
    );

    std::cout
        << "\nCOMPLEXITY\n"
        << "structural generation = O(I * digit_work)\n"
        << "geometry verification = O(I)\n"
        << "all uint64_t arithmetic is overflow-checked\n"
        << "I = product_i(q_i+1)-1\n";

    std::cout
        << "FINISHED EXPERIMENT 211\n";

    return 0;
}