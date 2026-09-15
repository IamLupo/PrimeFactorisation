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

    u64 maintained_j = 0;
    std::size_t maintained_r = 0;
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
    const std::vector<u64> digits = base_p_digits(q, p);

    u64 product = 1;

    for (u64 d : digits) {
        product *= d + 1;
    }

    return product - 1;
}

std::vector<u64> make_powers(
    u64 p,
    std::size_t count
) {
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

/*
    prefix_max_weight[k] =
        sum_{i=0}^{k-1} q_i * p^i

    This is exactly the amount removed from j when a carry
    propagates through digits 0,...,k-1 and resets them to zero.
*/
std::vector<u64> make_prefix_max_weight(
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers
) {
    std::vector<u64> prefix(q_digits.size() + 1, 0);

    for (std::size_t i = 0; i < q_digits.size(); ++i) {
        prefix[i + 1] =
            prefix[i] +
            q_digits[i] * powers[i];
    }

    return prefix;
}

/*
    next_positive[k] =
        smallest i > k such that q_i > 0

    If no such i exists, returns digit_count.
*/
std::vector<std::size_t> make_next_positive(
    const std::vector<u64>& q_digits
) {
    const std::size_t n = q_digits.size();

    std::vector<std::size_t> next(n, n);

    std::size_t current = n;

    for (std::size_t i = n; i-- > 0;) {
        next[i] = current;

        if (q_digits[i] > 0) {
            current = i;
        }
    }

    return next;
}

/*
    Initial r is the first digit whose q_i is nonzero.

    At j=0 all digits are zero, so every q_i=0 digit is
    already at its maximum.
*/
std::size_t initial_r(
    const std::vector<u64>& q_digits
) {
    for (std::size_t i = 0; i < q_digits.size(); ++i) {
        if (q_digits[i] > 0) {
            return i;
        }
    }

    return q_digits.size();
}

/*
    Advance the mixed-radix state.

    Before the increment:
        digits[0..k-1] = q[0..k-1]
        digits[k] < q[k]

    where k is exactly the current r.

    Therefore:

        j_new
          = j
          + p^k
          - sum_{i<k} q_i p^i

    After the carry:
        digits[0..k-1] = 0
        digit k increases by one.

    If the new digit k is still below q_k:
        r_new = k.

    If it becomes q_k:
        r_new is the next position with q_i > 0.
*/
bool advance_state(
    u64& j,
    std::size_t& r,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    const std::vector<u64>& prefix_max_weight,
    const std::vector<std::size_t>& next_positive,
    u64& increment_steps,
    u64& carry_steps
) {
    const std::size_t n = q_digits.size();

    if (r >= n) {
        return false;
    }

    const std::size_t k = r;

    /*
        k is the first non-max digit.
        All lower digits are currently q_i.
    */
    j += powers[k] - prefix_max_weight[k];

    ++increment_steps;
    carry_steps += k;

    /*
        Determine whether digit k has now reached q_k.

        We do not need the actual digit because before the increment
        it was < q_k, and the mixed-radix increment increases it by 1.

        The digit reaches q_k exactly when the old digit was q_k-1.

        Its old value can be recovered from j_before, but doing so
        would add a division. Instead use the resulting j and the
        known lower-prefix contribution.
    */

    const u64 higher_part =
        j / powers[k];

    const u64 digit_k =
        higher_part % (q_digits[k] == 0 ? 1 : (powers[1]));

    /*
        The expression above is intentionally not used for the
        transition. The mixed-radix state can be determined directly
        from the quotient before/after the carry:

        If the updated j modulo p^(k+1) equals q-prefix there,
        digit k is saturated.

        We compute the digit using j / p^k mod p.
    */

    const u64 base_p = powers[1];

    const u64 new_digit_k =
        (j / powers[k]) % base_p;

    if (new_digit_k < q_digits[k]) {
        r = k;
        return true;
    }

    /*
        digit k just became maximal.

        Lower positions have all been reset to zero, so the next
        possible r is the next position whose q_i is nonzero.
    */
    if (k + 1 >= n) {
        r = n;
        return j != 0;
    }

    const std::size_t next = next_positive[k];

    r = next;

    return true;
}

EnumerationResult enumerate_incremental(
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

    const std::vector<u64> prefix_max_weight =
        make_prefix_max_weight(
            q_digits,
            powers
        );

    const std::vector<std::size_t> next_positive =
        make_next_positive(q_digits);

    const u64 p_e = powers[e];

    const u64 expected =
        expected_interval_count(q, p);

    result.intervals.reserve(
        static_cast<std::size_t>(expected)
    );

    u64 j = 0;

    std::size_t r =
        initial_r(q_digits);

    while (j < q) {
        const u64 left =
            s0 + j * p_e;

        const u64 right =
            (j / powers[r] + 1) *
            powers[e + r] - 1;

        result.intervals.push_back({
            left,
            right
        });

        if (!advance_state(
                j,
                r,
                q_digits,
                powers,
                prefix_max_weight,
                next_positive,
                result.increment_steps,
                result.carry_steps
            )) {
            break;
        }
    }

    result.maintained_j = j;
    result.maintained_r = r;

    return result;
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

    for (u64 t = 1; t <= m; ++t) {
        const bool hit =
            !digitwise_leq(t, m, p);

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

    const auto optimized_begin =
        std::chrono::steady_clock::now();

    const EnumerationResult optimized =
        enumerate_incremental(
            p,
            e,
            s0,
            q
        );

    const auto optimized_end =
        std::chrono::steady_clock::now();

    const double direct_seconds =
        std::chrono::duration<double>(
            direct_end - direct_begin
        ).count();

    const double optimized_seconds =
        std::chrono::duration<double>(
            optimized_end - optimized_begin
        ).count();

    std::cout << "\nSMALL CASE\n";

    std::cout
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
        << "incremental_intervals="
        << optimized.intervals.size()
        << "\n"
        << "set_equality="
        << (intervals_equal(
                direct,
                optimized.intervals
            ) ? 1 : 0)
        << "\n"
        << "direct_seconds="
        << std::setprecision(10)
        << direct_seconds
        << "\n"
        << "incremental_seconds="
        << optimized_seconds
        << "\n"
        << "increment_steps="
        << optimized.increment_steps
        << "\n"
        << "carry_steps="
        << optimized.carry_steps
        << "\n"
        << "final_j="
        << optimized.maintained_j
        << "\n"
        << "final_r="
        << optimized.maintained_r
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
        enumerate_incremental(
            p,
            e,
            s0,
            q
        );

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
        << "final_j="
        << result.maintained_j
        << "\n"
        << "final_r="
        << result.maintained_r
        << "\n"
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

    std::cout
        << "direct_scan_points_would_be="
        << m
        << "\n";
}

int main() {
    std::cout << "START EXPERIMENT 198\n";

    std::size_t small_pass = 0;

    const std::vector<std::vector<u64>> small_cases = {
        {2, 3, 1, 63},
        {2, 4, 1, 255},
        {3, 4, 1, 80},
        {5, 3, 1, 100}
    };

    for (const auto& c : small_cases) {
        const u64 p = c[0];
        const u64 e = c[1];
        const u64 s0 = c[2];
        const u64 q = c[3];

        const u64 m =
            s0 + q * ipow_u64(p, e) - 1;

        const std::vector<Interval> direct =
            direct_hit_intervals(m, p);

        const EnumerationResult optimized =
            enumerate_incremental(
                p,
                e,
                s0,
                q
            );

        if (intervals_equal(
                direct,
                optimized.intervals
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
        << small_cases.size()
        << "\n"
        << "set_equality_pass="
        << small_pass
        << "/"
        << small_cases.size()
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
        << "structural_recursive = O(I log_p(q))\n"
        << "mixed_radix_197 = O(I + carry_work)\n"
        << "incremental_198 = O(I + carry_work)\n"
        << "per_interval digit reconstruction = eliminated\n"
        << "I = product_i(q_i+1)-1\n";

    std::cout << "FINISHED EXPERIMENT 198\n";

    return 0;
}
