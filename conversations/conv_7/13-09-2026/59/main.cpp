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

struct BoundaryGapStats {
    bool first_gap_ok = true;
    bool last_gap_ok = true;

    u64 first_gap = 0;
    u64 last_gap = 0;

    u64 internal_checked = 0;
    u64 internal_failures = 0;

    u64 internal_total = 0;
    u64 expected_internal_total = 0;

    u64 min_internal_gap =
        std::numeric_limits<u64>::max();

    u64 max_internal_gap = 0;
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

u64 digit_product_plus_one(
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
        digit_product_plus_one(digits) -
        1;
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

        /*
            Independent per-type count check.
        */
        const u64 generated_type =
            result.generated -
            before;

        u64 expected_type =
            q_digits[r];

        for (std::size_t i = r + 1;
             i < q_digits.size();
             ++i) {

            expected_type =
                safe_multiply(
                    expected_type,
                    q_digits[i] + 1
                );
        }

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

BoundaryGapStats analyze_boundary_gaps(
    const std::vector<Interval>& intervals,
    u64 s0,
    u64 m
) {
    BoundaryGapStats result;

    if (intervals.empty()) {
        result.first_gap_ok = false;
        result.last_gap_ok = false;
        return result;
    }

    /*
        MISS before first HIT:

            [0, first.left-1]

        Its size must be s0.
    */
    result.first_gap =
        intervals.front().left;

    result.first_gap_ok =
        result.first_gap == s0;

    /*
        MISS after last HIT:

            [last.right+1, m]

        Its size must also be s0.
    */
    if (intervals.back().right > m) {
        result.last_gap_ok = false;
        result.last_gap = 0;
    } else {
        result.last_gap =
            m -
            intervals.back().right;

        result.last_gap_ok =
            result.last_gap == s0;
    }

    /*
        Internal gaps.
    */
    if (intervals.size() >= 2) {
        result.expected_internal_total =
            safe_multiply(
                s0,
                intervals.size() - 1
            );
    }

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        const Interval& previous =
            intervals[i - 1];

        const Interval& current =
            intervals[i];

        ++result.internal_checked;

        if (current.left <=
            previous.right) {

            ++result.internal_failures;
            continue;
        }

        const u64 gap =
            current.left -
            previous.right -
            1;

        result.internal_total =
            safe_add(
                result.internal_total,
                gap
            );

        result.min_internal_gap =
            std::min(
                result.min_internal_gap,
                gap
            );

        result.max_internal_gap =
            std::max(
                result.max_internal_gap,
                gap
            );

        if (gap != s0) {
            ++result.internal_failures;
        }
    }

    if (result.min_internal_gap ==
        std::numeric_limits<u64>::max()) {

        result.min_internal_gap = 0;
    }

    return result;
}

u64 geometry_miss_count(
    u64 s0,
    u64 interval_count
) {
    /*
        [0,m] contains:

            s0 MISS points before first HIT
            s0 MISS points in each of I-1 internal gaps
            s0 MISS points after last HIT

        Therefore:

            s0 + s0(I-1) + s0 = s0(I+1)
    */
    const u64 blocks =
        interval_count + 1;

    return
        safe_multiply(
            s0,
            blocks
        );
}

u64 digit_miss_count(
    u64 m,
    u64 p
) {
    const std::vector<u64> digits =
        base_p_digits(
            m,
            p
        );

    return
        digit_product_plus_one(
            digits
        );
}

u64 structural_miss_count(
    u64 s0,
    u64 q,
    u64 p
) {
    const std::vector<u64> q_digits =
        base_p_digits(
            q,
            p
        );

    const u64 product =
        digit_product_plus_one(
            q_digits
        );

    return
        safe_multiply(
            s0,
            product
        );
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

        sort_intervals(
            result.intervals
        );

        const auto sort_end =
            std::chrono::steady_clock::now();

        const BoundaryGapStats gaps =
            analyze_boundary_gaps(
                result.intervals,
                s0,
                m
            );

        const u64 geometry_miss =
            geometry_miss_count(
                s0,
                result.generated
            );

        const u64 digit_miss =
            digit_miss_count(
                m,
                p
            );

        const u64 structural_miss =
            structural_miss_count(
                s0,
                q,
                p
            );

        const auto end =
            std::chrono::steady_clock::now();

        const bool enumeration_pass =
            result.invariant_ok &&
            result.generated ==
                expected_intervals;

        const bool first_gap_pass =
            gaps.first_gap_ok;

        const bool last_gap_pass =
            gaps.last_gap_ok;

        const bool internal_gap_pass =
            gaps.internal_failures == 0;

        const bool internal_sum_pass =
            gaps.internal_total ==
            gaps.expected_internal_total;

        const bool geometry_miss_pass =
            geometry_miss ==
            structural_miss;

        const bool digit_miss_pass =
            digit_miss ==
            structural_miss;

        const bool three_way_pass =
            geometry_miss ==
            digit_miss &&
            digit_miss ==
            structural_miss;

        const bool final_pass =
            enumeration_pass &&
            first_gap_pass &&
            last_gap_pass &&
            internal_gap_pass &&
            internal_sum_pass &&
            geometry_miss_pass &&
            digit_miss_pass &&
            three_way_pass;

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
            << expected_intervals
            << "\n"
            << "generated_intervals="
            << result.generated
            << "\n"
            << "enumeration_pass="
            << (enumeration_pass ? 1 : 0)
            << "\n"
            << "invariant_ok="
            << (result.invariant_ok ? 1 : 0)
            << "\n"
            << "first_gap="
            << gaps.first_gap
            << "\n"
            << "expected_first_gap="
            << s0
            << "\n"
            << "first_gap_pass="
            << (first_gap_pass ? 1 : 0)
            << "\n"
            << "last_gap="
            << gaps.last_gap
            << "\n"
            << "expected_last_gap="
            << s0
            << "\n"
            << "last_gap_pass="
            << (last_gap_pass ? 1 : 0)
            << "\n"
            << "internal_gap_checks="
            << gaps.internal_checked
            << "\n"
            << "internal_gap_failures="
            << gaps.internal_failures
            << "\n"
            << "internal_gap_pass="
            << (internal_gap_pass ? 1 : 0)
            << "\n"
            << "internal_gap_total="
            << gaps.internal_total
            << "\n"
            << "expected_internal_gap_total="
            << gaps.expected_internal_total
            << "\n"
            << "internal_gap_sum_pass="
            << (internal_sum_pass ? 1 : 0)
            << "\n"
            << "min_internal_gap="
            << gaps.min_internal_gap
            << "\n"
            << "max_internal_gap="
            << gaps.max_internal_gap
            << "\n"
            << "geometry_miss_count="
            << geometry_miss
            << "\n"
            << "digit_miss_count="
            << digit_miss
            << "\n"
            << "structural_miss_count="
            << structural_miss
            << "\n"
            << "geometry_miss_pass="
            << (geometry_miss_pass ? 1 : 0)
            << "\n"
            << "digit_miss_pass="
            << (digit_miss_pass ? 1 : 0)
            << "\n"
            << "three_way_miss_pass="
            << (three_way_pass ? 1 : 0)
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
    std::uniform_int_distribution<u64> dist(
        0,
        p - 1
    );

    u64 q = 0;
    u64 power = 1;

    for (std::size_t i = 0;
         i < digits;
         ++i) {

        const u64 d =
            dist(rng);

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

void run_random_suite(
    u64 p,
    std::size_t digits,
    u64 e,
    std::size_t cases,
    u64 seed
) {
    std::mt19937_64 rng(seed);

    u64 failures = 0;
    u64 skipped = 0;

    for (std::size_t i = 0;
         i < cases;
         ++i) {

        const u64 q =
            random_q(
                rng,
                p,
                digits
            );

        const u64 s0 =
            1 + (i % 11);

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
                ++skipped;
                continue;
            }

            const u64 qp =
                q * p_e;

            if (qp >
                std::numeric_limits<u64>::max() -
                    s0) {
                ++skipped;
                continue;
            }

            const u64 m =
                qp + s0 - 1;

            EnumerationResult result =
                enumerate_by_type(
                    p,
                    e,
                    s0,
                    q
                );

            if (!result.arithmetic_safe) {
                ++skipped;
                continue;
            }

            sort_intervals(
                result.intervals
            );

            const BoundaryGapStats gaps =
                analyze_boundary_gaps(
                    result.intervals,
                    s0,
                    m
                );

            const u64 geometry_miss =
                geometry_miss_count(
                    s0,
                    result.generated
                );

            const u64 digit_miss =
                digit_miss_count(
                    m,
                    p
                );

            const u64 structural_miss =
                structural_miss_count(
                    s0,
                    q,
                    p
                );

            const bool pass =
                result.invariant_ok &&
                result.generated ==
                    expected_interval_count(
                        q,
                        p
                    ) &&
                gaps.first_gap_ok &&
                gaps.last_gap_ok &&
                gaps.internal_failures == 0 &&
                gaps.internal_total ==
                    gaps.expected_internal_total &&
                geometry_miss ==
                    digit_miss &&
                digit_miss ==
                    structural_miss;

            if (!pass) {
                ++failures;
            }
        } catch (...) {
            ++skipped;
        }
    }

    std::cout
        << "\nRANDOM SUITE\n"
        << "p=" << p
        << " digits=" << digits
        << " e=" << e
        << " cases=" << cases
        << "\n"
        << "failures="
        << failures
        << "\n"
        << "skipped="
        << skipped
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 215\n";

    run_case(
        2,
        3,
        1,
        63,
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

    run_case(
        3,
        5,
        27,
        80,
        true
    );

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

    run_random_suite(
        2,
        16,
        20,
        300,
        0x12345678ULL
    );

    run_random_suite(
        3,
        10,
        10,
        300,
        0x87654321ULL
    );

    run_random_suite(
        5,
        8,
        10,
        300,
        0xCAFEBABEULL
    );

    std::cout
        << "\nIDENTITY\n"
        << "first MISS block = [0,s0-1], length s0\n"
        << "internal MISS blocks = I-1, each length s0\n"
        << "last MISS block = [last_right+1,m], length s0\n"
        << "total MISS blocks = I+1\n"
        << "MISS_[0,m] = s0*(I+1)\n"
        << "I+1 = product_i(q_i+1)\n"
        << "MISS_[0,m] = s0*product_i(q_i+1)\n"
        << "MISS_[0,m] = product_i(m_i+1)\n"
        << "FINISHED EXPERIMENT 215\n";

    return 0;
}
