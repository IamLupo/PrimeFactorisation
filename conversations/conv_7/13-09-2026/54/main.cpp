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

struct TypeStats {
    std::size_t r = 0;
    u64 expected_count = 0;
    u64 generated_count = 0;
};

struct EnumerationResult {
    std::vector<Interval> intervals;
    std::vector<TypeStats> types;

    u64 generated = 0;
    bool invariant_ok = true;
};

u64 ipow_u64(
    u64 p,
    std::size_t e
) {
    u64 result = 1;

    for (std::size_t i = 0; i < e; ++i) {
        result *= p;
    }

    return result;
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
            powers[i - 1] * p;
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
        product *= d + 1;
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

        count *=
            q_digits[i] + 1;
    }

    return count;
}

bool digitwise_leq(
    u64 a,
    u64 b,
    u64 p
) {
    while (a > 0 || b > 0) {
        const u64 da = a % p;
        const u64 db = b % p;

        if (da > db) {
            return false;
        }

        a /= p;
        b /= p;
    }

    return true;
}

std::vector<Interval> direct_hit_intervals(
    u64 m,
    u64 p
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1;
         t <= m;
         ++t) {

        const bool hit =
            !digitwise_leq(
                t,
                m,
                p
            );

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

    for (std::size_t i = 0;
         i < a.size();
         ++i) {

        if (a[i].left != b[i].left ||
            a[i].right != b[i].right) {

            return false;
        }
    }

    return true;
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

/*
    Verify the exact r-type condition for j.

        j_i = q_i      for i < r
        j_r < q_r
        j_i <= q_i    for i > r

    j=0 is explicitly valid when q_r > 0.
*/
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
            digit >= q_digits[i]) {

            return false;
        }
    }

    return true;
}

void enumerate_high_digits(
    std::size_t r,
    u64 s0,
    u64 p,
    u64 e,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    std::size_t position,
    u64 current_j,
    EnumerationResult& result,
    TypeStats& stats
) {
    const std::size_t n =
        q_digits.size();

    if (position >= n) {
        const u64 j =
            current_j;

        /*
            j=0 is VALID.
            It is precisely the first state of r=0 when q_0>0.
        */
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

        const u64 left =
            s0 +
            j * powers[e];

        const u64 right =
            (j / powers[r] + 1) *
            powers[e + r] -
            1;

        result.intervals.push_back({
            left,
            right
        });

        ++stats.generated_count;
        ++result.generated;

        return;
    }

    for (u64 d = 0;
         d <= q_digits[position];
         ++d) {

        const u64 next_j =
            current_j +
            d * powers[position];

        enumerate_high_digits(
            r,
            s0,
            p,
            e,
            q_digits,
            powers,
            position + 1,
            next_j,
            result,
            stats
        );
    }
}

void enumerate_type(
    std::size_t r,
    u64 s0,
    u64 p,
    u64 e,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    EnumerationResult& result,
    TypeStats& stats
) {
    u64 fixed_low = 0;

    /*
        i < r:
            j_i = q_i
    */
    for (std::size_t i = 0;
         i < r;
         ++i) {

        fixed_low +=
            q_digits[i] *
            powers[i];
    }

    /*
        i = r:
            j_r = 0,...,q_r-1
    */
    for (u64 d = 0;
         d < q_digits[r];
         ++d) {

        const u64 base =
            fixed_low +
            d * powers[r];

        enumerate_high_digits(
            r,
            s0,
            p,
            e,
            q_digits,
            powers,
            r + 1,
            base,
            result,
            stats
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

        TypeStats stats;

        stats.r = r;

        stats.expected_count =
            expected_type_count(
                q_digits,
                r
            );

        enumerate_type(
            r,
            s0,
            p,
            e,
            q_digits,
            powers,
            result,
            stats
        );

        if (stats.generated_count !=
            stats.expected_count) {

            result.invariant_ok = false;
        }

        result.types.push_back(
            stats
        );
    }

    if (result.generated != expected) {
        result.invariant_ok = false;
    }

    return result;
}

void print_type_stats(
    const std::vector<TypeStats>& types
) {
    for (const TypeStats& type : types) {
        std::cout
            << "type_r=" << type.r
            << " expected="
            << type.expected_count
            << " generated="
            << type.generated_count
            << "\n";
    }
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

void run_small_case(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const u64 m =
        s0 +
        q * ipow_u64(p, e) -
        1;

    const u64 expected =
        expected_interval_count(
            q,
            p
        );

    const auto direct_begin =
        std::chrono::steady_clock::now();

    const std::vector<Interval> direct =
        direct_hit_intervals(
            m,
            p
        );

    const auto direct_end =
        std::chrono::steady_clock::now();

    const auto structural_begin =
        std::chrono::steady_clock::now();

    EnumerationResult result =
        enumerate_by_type(
            p,
            e,
            s0,
            q
        );

    const auto structural_generated_end =
        std::chrono::steady_clock::now();

    sort_intervals(
        result.intervals
    );

    const auto structural_end =
        std::chrono::steady_clock::now();

    const double direct_seconds =
        std::chrono::duration<double>(
            direct_end -
            direct_begin
        ).count();

    const double structural_seconds =
        std::chrono::duration<double>(
            structural_end -
            structural_begin
        ).count();

    const bool equality =
        intervals_equal(
            direct,
            result.intervals
        );

    std::cout
        << "\nSMALL CASE\n"
        << "p=" << p
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
        << " m=" << m
        << "\n"
        << "expected_intervals="
        << expected
        << "\n"
        << "direct_intervals="
        << direct.size()
        << "\n"
        << "structural_intervals="
        << result.intervals.size()
        << "\n"
        << "set_equality="
        << (equality ? 1 : 0)
        << "\n"
        << "invariant_ok="
        << (result.invariant_ok ? 1 : 0)
        << "\n"
        << "direct_seconds="
        << std::setprecision(10)
        << direct_seconds
        << "\n"
        << "structural_seconds="
        << structural_seconds
        << "\n"
        << "types="
        << result.types.size()
        << "\n";

    print_type_stats(
        result.types
    );

    std::cout
        << "sample_count="
        << std::min<std::size_t>(
            10,
            result.intervals.size()
        )
        << "\n";

    print_sample(
        result.intervals,
        10
    );

    (void)structural_generated_end;
}

void run_large_case(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
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

    const bool pass =
        result.invariant_ok &&
        result.generated == expected &&
        result.intervals.size() == expected;

    std::cout
        << "\nLARGE CASE\n"
        << "p=" << p
        << " e=" << e
        << " s0=" << s0
        << " q=" << q
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
        << "invariant_ok="
        << (result.invariant_ok ? 1 : 0)
        << "\n"
        << "generation_seconds="
        << std::setprecision(10)
        << generation_seconds
        << "\n"
        << "generation_intervals_per_second="
        << (
            generation_seconds > 0.0
                ? static_cast<double>(
                    result.generated
                ) / generation_seconds
                : 0.0
        )
        << "\n"
        << "total_seconds="
        << total_seconds
        << "\n"
        << "type_count="
        << result.types.size()
        << "\n";

    print_sample(
        result.intervals,
        10
    );
}

int main() {
    std::cout
        << "START EXPERIMENT 209\n";

    run_small_case(
        2,
        3,
        1,
        63
    );

    run_small_case(
        2,
        4,
        1,
        255
    );

    run_small_case(
        3,
        4,
        1,
        80
    );

    run_small_case(
        5,
        3,
        1,
        100
    );

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

    std::cout
        << "\nCOMPLEXITY\n"
        << "C_r = q_r * product_{i>r}(q_i+1)\n"
        << "sum_r C_r = product_i(q_i+1)-1\n"
        << "structural emitted states = I\n"
        << "I = product_i(q_i+1)-1\n"
        << "verification sort = O(I log I)\n";

    std::cout
        << "FINISHED EXPERIMENT 209\n";

    return 0;
}
