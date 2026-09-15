#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 left;
    u64 right;
};

struct EnumerationResult {
    u64 interval_count = 0;
    u64 increment_steps = 0;
    u64 reset_steps = 0;
    u64 r_scan_steps = 0;
    u64 r_scan_events = 0;
    u64 final_j = 0;
    std::size_t final_r = 0;
    bool invariant_ok = true;
};

u64 ipow_u64(
    u64 p,
    std::size_t e
) {
    u64 result = 1;

    for (std::size_t i = 0;
         i < e;
         ++i) {

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
        return false;
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
        r > 0:

        digits[0..r-1] are all maximal.
        Reset them to zero.
        Increment digit r.

        The first reset position with q_i > 0
        becomes the new r.
    */
    if (r > 0) {
        const std::size_t k = r;

        std::size_t first_reset_positive = n;

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

            if (q_digits[i] > 0 &&
                first_reset_positive == n) {

                first_reset_positive = i;
            }
        }

        ++digits[k];
        j += powers[k];

        if (first_reset_positive < n) {
            r = first_reset_positive;
            return true;
        }

        if (digits[k] < q_digits[k]) {
            r = k;
            return true;
        }

        /*
            k saturated and there is no lower non-max digit.
            Find the next higher non-max digit.
        */
        ++result.r_scan_events;

        for (std::size_t i = k + 1;
             i < n;
             ++i) {

            ++result.r_scan_steps;

            if (digits[i] < q_digits[i]) {
                r = i;
                return true;
            }
        }

        r = n;
        return true;
    }

    /*
        r == 0.
    */
    ++digits[0];
    j += powers[0];

    if (digits[0] < q_digits[0]) {
        return true;
    }

    /*
        Digit zero saturated.
        Search upward.
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

    r = n;

    return true;
}

EnumerationResult enumerate(
    u64 p,
    u64 e,
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

    std::vector<u64> digits(
        digit_count,
        0
    );

    u64 j = 0;

    std::size_t r =
        initial_r(q_digits);

    const u64 expected =
        expected_interval_count(
            q,
            p
        );

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

        ++result.interval_count;

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

    if (result.interval_count != expected) {
        result.invariant_ok = false;
    }

    result.final_j = j;
    result.final_r = r;

    return result;
}

u64 random_q_dense(
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

        const u64 digit =
            dist(rng);

        q +=
            digit *
            power;

        power *= p;
    }

    /*
        Avoid q=0 and q with only one trivial digit.
    */
    if (q == 0) {
        q = 1;
    }

    return q;
}

u64 random_q_sparse(
    std::mt19937_64& rng,
    u64 p,
    std::size_t digits
) {
    std::uniform_int_distribution<u64> digit_dist(
        0,
        p - 1
    );

    std::bernoulli_distribution active(
        0.20
    );

    u64 q = 0;
    u64 power = 1;

    for (std::size_t i = 0;
         i < digits;
         ++i) {

        u64 d = 0;

        if (active(rng)) {
            d = digit_dist(rng);

            if (d == 0) {
                d = 1;
            }
        }

        q += d * power;

        power *= p;
    }

    if (q == 0) {
        q = 1;
    }

    return q;
}

void run_fixed_case(
    u64 p,
    u64 e,
    u64 q
) {
    const auto begin =
        std::chrono::steady_clock::now();

    const EnumerationResult result =
        enumerate(
            p,
            e,
            q
        );

    const auto end =
        std::chrono::steady_clock::now();

    const double seconds =
        std::chrono::duration<double>(
            end - begin
        ).count();

    const u64 I =
        expected_interval_count(
            q,
            p
        );

    const bool scan_bound =
        result.r_scan_steps <= I;

    std::cout
        << "\nFIXED CASE\n"
        << "p=" << p
        << " e=" << e
        << " q=" << q
        << "\n"
        << "I=" << I
        << "\n"
        << "interval_count="
        << result.interval_count
        << "\n"
        << "enumeration_pass="
        << (
            result.invariant_ok &&
            result.interval_count == I &&
            result.final_j == q
                ? 1
                : 0
        )
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
        << "r_scan_leq_I="
        << (scan_bound ? 1 : 0)
        << "\n"
        << "total_work="
        << (
            result.reset_steps +
            result.r_scan_steps
        )
        << "\n"
        << "work_per_interval="
        << std::setprecision(8)
        << static_cast<double>(
            result.reset_steps +
            result.r_scan_steps
        ) /
        static_cast<double>(I)
        << "\n"
        << "seconds="
        << std::setprecision(10)
        << seconds
        << "\n";
}

void run_random_suite(
    u64 p,
    std::size_t digits,
    std::size_t cases,
    bool sparse,
    u64 seed
) {
    std::mt19937_64 rng(seed);

    u64 failed = 0;
    u64 scan_bound_failed = 0;

    u64 total_I = 0;
    u64 total_reset = 0;
    u64 total_scan = 0;

    double max_work_ratio = 0.0;
    double max_scan_ratio = 0.0;

    for (std::size_t case_index = 0;
         case_index < cases;
         ++case_index) {

        const u64 q =
            sparse
                ? random_q_sparse(
                    rng,
                    p,
                    digits
                )
                : random_q_dense(
                    rng,
                    p,
                    digits
                );

        const EnumerationResult result =
            enumerate(
                p,
                digits / 2,
                q
            );

        const u64 I =
            expected_interval_count(
                q,
                p
            );

        const bool pass =
            result.invariant_ok &&
            result.interval_count == I &&
            result.final_j == q;

        const bool scan_ok =
            result.r_scan_steps <= I;

        if (!pass) {
            ++failed;
        }

        if (!scan_ok) {
            ++scan_bound_failed;
        }

        total_I += I;
        total_reset += result.reset_steps;
        total_scan += result.r_scan_steps;

        const double work_ratio =
            static_cast<double>(
                result.reset_steps +
                result.r_scan_steps
            ) /
            static_cast<double>(I);

        const double scan_ratio =
            static_cast<double>(
                result.r_scan_steps
            ) /
            static_cast<double>(I);

        max_work_ratio =
            std::max(
                max_work_ratio,
                work_ratio
            );

        max_scan_ratio =
            std::max(
                max_scan_ratio,
                scan_ratio
            );
    }

    std::cout
        << "\nRANDOM SUITE\n"
        << "p=" << p
        << " digits=" << digits
        << " cases=" << cases
        << " mode="
        << (sparse ? "sparse" : "dense")
        << "\n"
        << "correctness_failures="
        << failed
        << "\n"
        << "scan_bound_failures="
        << scan_bound_failed
        << "\n"
        << "total_I="
        << total_I
        << "\n"
        << "total_reset_steps="
        << total_reset
        << "\n"
        << "total_r_scan_steps="
        << total_scan
        << "\n"
        << "total_work="
        << total_reset + total_scan
        << "\n"
        << "total_work_per_interval="
        << std::setprecision(10)
        << static_cast<double>(
            total_reset + total_scan
        ) /
        static_cast<double>(
            total_I
        )
        << "\n"
        << "max_work_ratio="
        << max_work_ratio
        << "\n"
        << "max_scan_ratio="
        << max_scan_ratio
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 207\n";

    run_fixed_case(
        2,
        40,
        1048575
    );

    run_fixed_case(
        5,
        20,
        390624
    );

    run_fixed_case(
        2,
        20,
        1073741825ULL
    );

    run_random_suite(
        2,
        20,
        1000,
        false,
        0x12345678ULL
    );

    run_random_suite(
        2,
        20,
        1000,
        true,
        0x87654321ULL
    );

    run_random_suite(
        3,
        12,
        1000,
        false,
        0x13579BDFULL
    );

    run_random_suite(
        5,
        10,
        1000,
        false,
        0x2468ACE0ULL
    );

    std::cout
        << "\nCOMPLEXITY\n"
        << "interval enumeration = O(I + reset_work + r_scan_work)\n"
        << "I = product_i(q_i+1)-1\n"
        << "experimental target: r_scan_work <= I\n"
        << "experimental target: total work = O(I)\n";

    std::cout
        << "FINISHED EXPERIMENT 207\n";

    return 0;
}
