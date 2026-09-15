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
    u64 reset_steps = 0;
    u64 r_scan_steps = 0;
    u64 r_scan_events = 0;

    u64 final_j = 0;
    std::size_t final_r = 0;

    bool invariant_ok = true;
};

u64 ipow_u64(u64 p, std::size_t e) {
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

bool state_invariant_holds(
    u64 j,
    const std::vector<u64>& digits,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    std::size_t r
) {
    u64 reconstructed = 0;

    for (std::size_t i = 0;
         i < digits.size();
         ++i) {

        if (digits[i] > q_digits[i]) {
            return false;
        }

        reconstructed +=
            digits[i] * powers[i];
    }

    if (reconstructed != j) {
        return false;
    }

    if (r >= digits.size()) {
        return j == 0;
    }

    for (std::size_t i = 0;
         i < r;
         ++i) {

        if (digits[i] != q_digits[i]) {
            return false;
        }
    }

    return digits[r] < q_digits[r];
}

/*
    Initial state j=0.

    r is the lowest digit for which q_i > 0,
    because every digit of j is initially zero.
*/
std::size_t initial_r(
    const std::vector<u64>& q_digits
) {
    for (std::size_t i = 0;
         i < q_digits.size();
         ++i) {

        if (q_digits[i] > 0) {
            return i;
        }
    }

    return q_digits.size();
}

/*
    Advance from one valid digitwise state j to the next.

    Invariant:

        r = min { i : digits[i] < q_i }

    Therefore every digit below r is exactly q_i.

    Case r > 0:

        reset digits[0..r-1] to zero,
        increment digit r.

        Since digit 0 becomes zero and q_0 > 0,
        the new r is 0.

    Case r = 0:

        increment digit 0.

        If it remains below q_0, r stays 0.

        If it reaches q_0, scan upward for the next
        digit which is still below its q-bound.
*/
bool advance_state(
    std::vector<u64>& digits,
    u64& j,
    std::size_t& r,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    EnumerationResult& result
) {
    const std::size_t n =
        digits.size();

    if (r >= n) {
        return false;
    }

    ++result.increment_steps;

    /*
        CASE 1:
        r > 0.

        All lower digits are maximal.
    */
    if (r > 0) {
        const std::size_t k = r;

        for (std::size_t i = 0;
             i < k;
             ++i) {

            if (digits[i] != q_digits[i]) {
                result.invariant_ok = false;
                return false;
            }

            j -=
                q_digits[i] *
                powers[i];

            digits[i] = 0;

            ++result.reset_steps;
        }

        /*
            digit k is non-maximal.
        */
        ++digits[k];
        j += powers[k];

        /*
            Because digit 0 was reset to zero,
            it is now the lowest non-maximal digit.
        */
        r = 0;

        return true;
    }

    /*
        CASE 2:
        r == 0.

        Increment digit 0.
    */
    ++digits[0];
    j += powers[0];

    /*
        Still non-maximal.
    */
    if (digits[0] < q_digits[0]) {
        r = 0;
        return true;
    }

    /*
        digit 0 reached q_0.

        The next valid state must find the lowest
        higher digit still below its bound.
    */
    ++result.r_scan_events;

    for (std::size_t i = 1;
         i < n;
         ++i) {

        ++result.r_scan_steps;

        if (digits[i] < q_digits[i]) {
            r = i;
            return true;
        }
    }

    /*
        Every digit has reached its maximum.

        Therefore j=q.
    */
    r = n;

    return true;
}

EnumerationResult enumerate_incremental(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    EnumerationResult result;

    const std::vector<u64> q_raw =
        base_p_digits(q, p);

    const std::size_t digit_count =
        q_raw.size();

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

    const u64 p_e =
        powers[e];

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

    std::vector<u64> digits(
        digit_count,
        0
    );

    u64 j = 0;

    std::size_t r =
        initial_r(q_digits);

    while (j < q) {
        if (!state_invariant_holds(
                j,
                digits,
                q_digits,
                powers,
                r
            )) {

            result.invariant_ok = false;
            break;
        }

        const u64 left =
            s0 +
            j * p_e;

        const u64 right =
            (j / powers[r] + 1) *
            powers[e + r] -
            1;

        result.intervals.push_back({
            left,
            right
        });

        if (!advance_state(
                digits,
                j,
                r,
                q_digits,
                powers,
                result
            )) {

            break;
        }
    }

    result.final_j = j;
    result.final_r = r;

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

    const auto incremental_begin =
        std::chrono::steady_clock::now();

    const EnumerationResult incremental =
        enumerate_incremental(
            p,
            e,
            s0,
            q
        );

    const auto incremental_end =
        std::chrono::steady_clock::now();

    const double direct_seconds =
        std::chrono::duration<double>(
            direct_end -
            direct_begin
        ).count();

    const double incremental_seconds =
        std::chrono::duration<double>(
            incremental_end -
            incremental_begin
        ).count();

    const bool equality =
        intervals_equal(
            direct,
            incremental.intervals
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
        << "incremental_intervals="
        << incremental.intervals.size()
        << "\n"
        << "set_equality="
        << (equality ? 1 : 0)
        << "\n"
        << "invariant_ok="
        << (incremental.invariant_ok ? 1 : 0)
        << "\n"
        << "direct_seconds="
        << std::setprecision(10)
        << direct_seconds
        << "\n"
        << "incremental_seconds="
        << incremental_seconds
        << "\n"
        << "increment_steps="
        << incremental.increment_steps
        << "\n"
        << "reset_steps="
        << incremental.reset_steps
        << "\n"
        << "r_scan_events="
        << incremental.r_scan_events
        << "\n"
        << "r_scan_steps="
        << incremental.r_scan_steps
        << "\n"
        << "final_j="
        << incremental.final_j
        << "\n"
        << "final_r="
        << incremental.final_r
        << "\n";
}

void run_large_case(
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
            end -
            begin
        ).count();

    const bool pass =
        result.intervals.size() == expected &&
        result.invariant_ok &&
        result.final_j == q;

    std::cout
        << "\nLARGE CASE\n"
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
        << "invariant_ok="
        << (result.invariant_ok ? 1 : 0)
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
            ) /
            seconds
            << "\n";
    }

    std::cout
        << "increment_steps="
        << result.increment_steps
        << "\n"
        << "reset_steps="
        << result.reset_steps
        << "\n"
        << "r_scan_events="
        << result.r_scan_events
        << "\n"
        << "r_scan_steps="
        << result.r_scan_steps
        << "\n"
        << "final_j="
        << result.final_j
        << "\n"
        << "final_r="
        << result.final_r
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
    std::cout
        << "START EXPERIMENT 204\n";

    const std::vector<std::vector<u64>> cases = {
        {2, 3, 1, 63},
        {2, 4, 1, 255},
        {3, 4, 1, 80},
        {5, 3, 1, 100}
    };

    std::size_t small_pass = 0;

    for (const auto& c : cases) {
        const u64 p = c[0];
        const u64 e = c[1];
        const u64 s0 = c[2];
        const u64 q = c[3];

        const u64 m =
            s0 +
            q * ipow_u64(p, e) -
            1;

        const std::vector<Interval> direct =
            direct_hit_intervals(
                m,
                p
            );

        const EnumerationResult incremental =
            enumerate_incremental(
                p,
                e,
                s0,
                q
            );

        if (
            intervals_equal(
                direct,
                incremental.intervals
            ) &&
            incremental.invariant_ok &&
            incremental.final_j == q
        ) {
            ++small_pass;
        }

        run_small_case(
            p,
            e,
            s0,
            q
        );
    }

    std::cout
        << "\nSMALL SUMMARY\n"
        << "cases="
        << cases.size()
        << "\n"
        << "set_equality_pass="
        << small_pass
        << "/"
        << cases.size()
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

    std::cout
        << "\nCOMPLEXITY\n"
        << "direct_scan = O(m log_p(m))\n"
        << "incremental = O(I + reset_work + r_scan_work)\n"
        << "r scan occurs only when digit 0 saturates\n"
        << "I = product_i(q_i+1)-1\n";

    std::cout
        << "FINISHED EXPERIMENT 204\n";

    return 0;
}
