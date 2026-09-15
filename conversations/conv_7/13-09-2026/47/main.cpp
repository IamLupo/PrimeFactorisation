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
    u64 r_updates = 0;

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

/*
    Doubly linked list of digit positions satisfying

        digits[i] < q_digits[i].

    The smallest such position is exactly r.
*/
struct NonMaxList {
    std::vector<std::size_t> prev;
    std::vector<std::size_t> next;

    std::size_t head;
    std::size_t tail;

    explicit NonMaxList(
        std::size_t n
    )
        : prev(n, n),
          next(n, n),
          head(n),
          tail(n) {
    }
};

void initialize_nonmax_list(
    NonMaxList& list,
    const std::vector<u64>& q_digits
) {
    const std::size_t n =
        q_digits.size();

    list.head = n;
    list.tail = n;

    for (std::size_t i = 0; i < n; ++i) {
        list.prev[i] = n;
        list.next[i] = n;
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (q_digits[i] == 0) {
            continue;
        }

        if (list.head == n) {
            list.head = i;
            list.tail = i;
        } else {
            list.next[list.tail] = i;
            list.prev[i] = list.tail;
            list.tail = i;
        }
    }
}

void remove_nonmax(
    NonMaxList& list,
    std::size_t i
) {
    const std::size_t n =
        list.prev.size();

    const std::size_t p = list.prev[i];
    const std::size_t q = list.next[i];

    if (p != n) {
        list.next[p] = q;
    } else {
        list.head = q;
    }

    if (q != n) {
        list.prev[q] = p;
    } else {
        list.tail = p;
    }

    list.prev[i] = n;
    list.next[i] = n;
}

void insert_nonmax_sorted(
    NonMaxList& list,
    std::size_t i
) {
    const std::size_t n =
        list.prev.size();

    /*
        Insert into sorted order.

        This is only used for positions reset by a carry.
        Those positions are below the carry position, so they
        are smaller than the current r. Therefore, in the normal
        transition, inserting at the front is sufficient.
    */
    if (list.head == n) {
        list.head = i;
        list.tail = i;
        return;
    }

    list.prev[i] = n;
    list.next[i] = list.head;
    list.prev[list.head] = i;
    list.head = i;
}

std::size_t current_r(
    const NonMaxList& list
) {
    return list.head;
}

bool state_invariant_holds(
    u64 j,
    const std::vector<u64>& digits,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    std::size_t r,
    const NonMaxList& list
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

    /*
        Every position before r must be maximal.
    */
    for (std::size_t i = 0;
         i < r;
         ++i) {

        if (digits[i] != q_digits[i]) {
            return false;
        }
    }

    if (r < digits.size()) {
        if (digits[r] >= q_digits[r]) {
            return false;
        }
    }

    /*
        Verify that the linked-list head really is r.
    */
    if (list.head != r) {
        return false;
    }

    return true;
}

bool advance_state(
    std::vector<u64>& digits,
    u64& j,
    NonMaxList& list,
    std::size_t& r,
    const std::vector<u64>& q_digits,
    const std::vector<u64>& powers,
    const std::vector<u64>& prefix_max_weight,
    EnumerationResult& result
) {
    const std::size_t n =
        digits.size();

    if (r >= n) {
        return false;
    }

    ++result.increment_steps;

    const std::size_t k = r;

    /*
        Case 1:
        r can increase without reaching q_r.
    */
    if (digits[k] + 1 < q_digits[k]) {
        ++digits[k];

        j += powers[k];

        return true;
    }

    /*
        Digit k reaches its maximum.

        All lower digits are currently maximal.
    */
    if (k > 0) {
        /*
            Remove their contribution in one arithmetic operation.
        */
        j -= prefix_max_weight[k];

        /*
            They become zero and therefore non-maximal again.

            They must be reinserted into the list.
        */
        for (std::size_t i = 0;
             i < k;
             ++i) {

            if (q_digits[i] == 0) {
                continue;
            }

            digits[i] = 0;

            insert_nonmax_sorted(
                list,
                i
            );

            ++result.carry_steps;
        }
    }

    /*
        Increment k to q_k.
    */
    ++digits[k];

    j += powers[k];

    /*
        k was the current head, so it is no longer non-maximal.
    */
    remove_nonmax(
        list,
        k
    );

    /*
        After the update, the new r is simply the new head.
    */
    r = current_r(list);

    ++result.r_updates;

    return true;
}

EnumerationResult enumerate_list(
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

    /*
        prefix_max_weight[k]
          =
        sum_{i<k} q_i p^i.
    */
    std::vector<u64> prefix_max_weight(
        digit_count + 1,
        0
    );

    for (std::size_t i = 0;
         i < digit_count;
         ++i) {

        prefix_max_weight[i + 1] =
            prefix_max_weight[i] +
            q_digits[i] * powers[i];
    }

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

    NonMaxList list(
        digit_count
    );

    initialize_nonmax_list(
        list,
        q_digits
    );

    u64 j = 0;

    std::size_t r =
        current_r(list);

    while (j < q) {
        if (!state_invariant_holds(
                j,
                digits,
                q_digits,
                powers,
                r,
                list
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
                list,
                r,
                q_digits,
                powers,
                prefix_max_weight,
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

    const auto optimized_begin =
        std::chrono::steady_clock::now();

    const EnumerationResult optimized =
        enumerate_list(
            p,
            e,
            s0,
            q
        );

    const auto optimized_end =
        std::chrono::steady_clock::now();

    const double direct_seconds =
        std::chrono::duration<double>(
            direct_end -
            direct_begin
        ).count();

    const double optimized_seconds =
        std::chrono::duration<double>(
            optimized_end -
            optimized_begin
        ).count();

    const bool equality =
        intervals_equal(
            direct,
            optimized.intervals
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
        << "list_intervals="
        << optimized.intervals.size()
        << "\n"
        << "set_equality="
        << (equality ? 1 : 0)
        << "\n"
        << "invariant_ok="
        << (optimized.invariant_ok ? 1 : 0)
        << "\n"
        << "direct_seconds="
        << std::setprecision(10)
        << direct_seconds
        << "\n"
        << "list_seconds="
        << optimized_seconds
        << "\n"
        << "increment_steps="
        << optimized.increment_steps
        << "\n"
        << "carry_steps="
        << optimized.carry_steps
        << "\n"
        << "r_updates="
        << optimized.r_updates
        << "\n"
        << "final_j="
        << optimized.final_j
        << "\n"
        << "final_r="
        << optimized.final_r
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
        enumerate_list(
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
        result.intervals.size() ==
            expected &&
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
        << "carry_steps="
        << result.carry_steps
        << "\n"
        << "r_updates="
        << result.r_updates
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
        << "START EXPERIMENT 202\n";

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

        const EnumerationResult optimized =
            enumerate_list(
                p,
                e,
                s0,
                q
            );

        if (
            intervals_equal(
                direct,
                optimized.intervals
            ) &&
            optimized.invariant_ok &&
            optimized.final_j == q
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
        << "reference_odometer = O(I + carry_work + r_scan_work)\n"
        << "list_enumeration = O(I + carry_work)\n"
        << "r lookup = O(1)\n"
        << "I = product_i(q_i+1)-1\n";

    std::cout
        << "FINISHED EXPERIMENT 202\n";

    return 0;
}
