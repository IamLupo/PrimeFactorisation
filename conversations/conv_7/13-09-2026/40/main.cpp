#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Parameters {
    u64 p;
    u64 e;
    u64 s0;
    u64 q;
    u64 pe;
    u64 m;
};

struct Interval {
    u64 start;
    u64 end;
};

struct Enumerator {
    Parameters par;
    std::vector<u64> digits_msb;
    std::vector<u64> powers;

    u64 count = 0;

    // For exact small validation.
    std::vector<Interval> intervals;

    // For large runs, only retain a few samples.
    std::vector<Interval> samples;
};

u64 pow_u64(
    u64 p,
    u64 e
) {
    u128 result = 1;

    for (u64 i = 0; i < e; ++i) {
        result *= p;
    }

    return static_cast<u64>(result);
}

Parameters make_parameters(
    u64 p,
    u64 e,
    u64 s0,
    u64 q
) {
    const u64 pe =
        pow_u64(
            p,
            e
        );

    const u128 m =
        static_cast<u128>(s0) +
        static_cast<u128>(q) *
            static_cast<u128>(pe) -
        1;

    return {
        p,
        e,
        s0,
        q,
        pe,
        static_cast<u64>(m)
    };
}

std::vector<u64> digits_msb(
    u64 p,
    u64 q
) {
    std::vector<u64> digits;

    while (q > 0) {
        digits.push_back(q % p);
        q /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    std::reverse(
        digits.begin(),
        digits.end()
    );

    return digits;
}

std::vector<u64> make_powers(
    u64 p,
    u64 max_power
) {
    std::vector<u64> powers(
        max_power + 1,
        1
    );

    for (u64 i = 1;
         i <= max_power;
         ++i) {

        powers[i] =
            static_cast<u64>(
                static_cast<u128>(
                    powers[i - 1]
                ) *
                static_cast<u128>(p)
            );
    }

    return powers;
}

bool is_hit(
    u64 p,
    u64 m,
    u64 x
) {
    while (x > 0 || m > 0) {
        const u64 xd = x % p;
        const u64 md = m % p;

        if (xd > md) {
            return true;
        }

        x /= p;
        m /= p;
    }

    return false;
}

std::vector<Interval> direct_scan(
    const Parameters &par
) {
    std::vector<Interval> result;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1;
         t <= par.m;
         ++t) {

        const bool hit =
            is_hit(
                par.p,
                par.m,
                t
            );

        if (hit && !inside) {
            inside = true;
            start = t;
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
            par.m
        });
    }

    return result;
}

Interval make_interval(
    const Enumerator &en,
    u64 j,
    u64 r
) {
    const u64 pr =
        en.powers[r];

    const u64 per =
        en.powers[
            en.par.e + r
        ];

    const u64 high =
        j / pr;

    const u64 start =
        static_cast<u64>(
            static_cast<u128>(en.par.s0) +
            static_cast<u128>(j) *
                static_cast<u128>(en.par.pe)
        );

    const u64 end =
        static_cast<u64>(
            static_cast<u128>(high + 1) *
                static_cast<u128>(per) -
            1
        );

    return {
        start,
        end
    };
}

/*
 * Exact comparison without relying on
 * Interval::operator==.
 */
bool same_intervals(
    const std::vector<Interval> &a,
    const std::vector<Interval> &b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0;
         i < a.size();
         ++i) {

        if (a[i].start != b[i].start ||
            a[i].end != b[i].end) {
            return false;
        }
    }

    return true;
}

u64 digit_product(
    u64 p,
    u64 q
) {
    u64 result = 1;

    while (q > 0) {
        result *=
            (q % p) + 1;

        q /= p;
    }

    return result;
}

u64 expected_interval_count(
    u64 p,
    u64 q
) {
    return digit_product(
        p,
        q
    ) - 1;
}

/*
 * Store samples without affecting enumeration.
 */
void store_sample(
    Enumerator &en,
    const Interval &interval
) {
    if (en.samples.size() < 10) {
        en.samples.push_back(interval);
        return;
    }

    if (en.count % 100000 == 0) {
        en.samples.push_back(interval);
    }
}

/*
 * Recursive digit-tree enumeration.

 * We enumerate j in increasing numerical order.

 * r is tracked while descending:
 *
 *   base-p position =
 *       number_of_digits - 1 - depth
 *
 * Whenever d < q_digit, this position is a candidate
 * for r. A later recursion level represents a lower digit,
 * so it correctly overwrites the previous candidate.
 */
void enumerate_recursive(
    Enumerator &en,
    std::size_t depth,
    bool equal_prefix,
    u64 current,
    u64 current_r
) {
    const std::size_t digit_count =
        en.digits_msb.size();

    if (depth == digit_count) {
        /*
         * equal_prefix means current == q.
         * q itself is not an interval start.
         */
        if (equal_prefix) {
            return;
        }

        const Interval interval =
            make_interval(
                en,
                current,
                current_r
            );

        ++en.count;

        if (en.intervals.capacity() > 0) {
            en.intervals.push_back(interval);
        }

        store_sample(
            en,
            interval
        );

        return;
    }

    const u64 q_digit =
        en.digits_msb[depth];

    const u64 max_digit =
        equal_prefix
        ? q_digit
        : en.par.p - 1;

    const u64 r_position =
        static_cast<u64>(
            digit_count - 1 - depth
        );

    for (u64 d = 0;
         d <= max_digit;
         ++d) {

        const bool next_equal =
            equal_prefix &&
            d == q_digit;

        u64 next_r =
            current_r;

        if (d < q_digit) {
            next_r = r_position;
        }

        const u64 next_current =
            current * en.par.p + d;

        enumerate_recursive(
            en,
            depth + 1,
            next_equal,
            next_current,
            next_r
        );
    }
}

/*
 * Prepare an enumerator.
 */
Enumerator make_enumerator(
    const Parameters &par,
    bool store_all
) {
    Enumerator en;

    en.par = par;

    en.digits_msb =
        digits_msb(
            par.p,
            par.q
        );

    /*
     * We need powers up to e + number_of_digits(q).
     */
    const u64 max_power =
        par.e +
        static_cast<u64>(
            en.digits_msb.size()
        );

    en.powers =
        make_powers(
            par.p,
            max_power
        );

    if (store_all) {
        en.intervals.reserve(
            static_cast<std::size_t>(
                expected_interval_count(
                    par.p,
                    par.q
                )
            )
        );
    }

    return en;
}

void run_enumeration(
    Enumerator &en
) {
    enumerate_recursive(
        en,
        0,
        true,
        0,
        UINT64_MAX
    );
}

void run_small_exact(
    u64 &cases,
    u64 &pass
) {
    struct TestCase {
        u64 p;
        u64 e;
        u64 s0;
        u64 q;
    };

    /*
     * These produce manageable m and manageable
     * interval counts for exact comparison.
     */
    const std::vector<TestCase> tests = {
        {2, 3, 1, 63},
        {2, 4, 1, 255},
        {3, 4, 1, 80},
        {5, 3, 1, 100}
    };

    for (const auto &test : tests) {
        const Parameters par =
            make_parameters(
                test.p,
                test.e,
                test.s0,
                test.q
            );

        ++cases;

        const auto direct_begin =
            std::chrono::steady_clock::now();

        const auto direct =
            direct_scan(
                par
            );

        const auto direct_end =
            std::chrono::steady_clock::now();

        Enumerator en =
            make_enumerator(
                par,
                true
            );

        const auto structural_begin =
            std::chrono::steady_clock::now();

        run_enumeration(
            en
        );

        const auto structural_end =
            std::chrono::steady_clock::now();

        const bool ok =
            same_intervals(
                direct,
                en.intervals
            );

        if (ok) {
            ++pass;
        }

        const std::chrono::duration<double>
            direct_elapsed =
                direct_end - direct_begin;

        const std::chrono::duration<double>
            structural_elapsed =
                structural_end -
                structural_begin;

        std::cout
            << "\nSMALL CASE\n";

        std::cout
            << "p=" << par.p
            << " e=" << par.e
            << " s0=" << par.s0
            << " q=" << par.q
            << " m=" << par.m
            << "\n";

        std::cout
            << "direct_intervals="
            << direct.size()
            << "\n";

        std::cout
            << "structural_intervals="
            << en.intervals.size()
            << "\n";

        std::cout
            << "set_equality="
            << (ok ? 1 : 0)
            << "\n";

        std::cout
            << "direct_seconds="
            << direct_elapsed.count()
            << "\n";

        std::cout
            << "structural_seconds="
            << structural_elapsed.count()
            << "\n";
    }
}

void run_large_structural() {
    struct TestCase {
        u64 p;
        u64 e;
        u64 s0;
        u64 q;
    };

    const std::vector<TestCase> tests = {
        /*
         * q = 2^20 - 1
         *
         * I = 1,048,575
         */
        {
            2,
            40,
            1,
            (1ULL << 20) - 1
        },

        /*
         * q = 5^8 - 1
         *
         * All eight base-5 digits are 4,
         * giving 5^8 - 1 intervals.
         */
        {
            5,
            15,
            1,
            390624
        },

        /*
         * Sparse q:
         *
         * only a few valid digitwise subnumbers.
         */
        {
            2,
            20,
            1,
            (1ULL << 30) + 1
        }
    };

    for (const auto &test : tests) {
        const Parameters par =
            make_parameters(
                test.p,
                test.e,
                test.s0,
                test.q
            );

        Enumerator en =
            make_enumerator(
                par,
                false
            );

        const u64 expected =
            expected_interval_count(
                par.p,
                par.q
            );

        const auto begin =
            std::chrono::steady_clock::now();

        run_enumeration(
            en
        );

        const auto end =
            std::chrono::steady_clock::now();

        const std::chrono::duration<double>
            elapsed =
                end - begin;

        std::cout
            << "\nLARGE CASE\n";

        std::cout
            << "p=" << par.p
            << " e=" << par.e
            << " s0=" << par.s0
            << " q=" << par.q
            << "\n";

        std::cout
            << "m=" << par.m
            << "\n";

        std::cout
            << "expected_intervals="
            << expected
            << "\n";

        std::cout
            << "generated_intervals="
            << en.count
            << "\n";

        std::cout
            << "enumeration_pass="
            << (
                en.count == expected
                ? 1
                : 0
            )
            << "\n";

        std::cout
            << "seconds="
            << elapsed.count()
            << "\n";

        if (elapsed.count() > 0.0) {
            std::cout
                << "intervals_per_second="
                << static_cast<double>(
                    en.count
                ) /
                elapsed.count()
                << "\n";
        }

        std::cout
            << "samples="
            << en.samples.size()
            << "\n";

        for (std::size_t i = 0;
             i < en.samples.size();
             ++i) {

            std::cout
                << "sample["
                << i
                << "]="
                << "["
                << en.samples[i].start
                << ","
                << en.samples[i].end
                << "]\n";
        }

        std::cout
            << "direct_scan_points_would_be="
            << par.m
            << "\n";
    }
}

int main() {
    std::cout
        << "START EXPERIMENT 195\n";

    u64 small_cases = 0;
    u64 small_pass = 0;

    run_small_exact(
        small_cases,
        small_pass
    );

    std::cout
        << "\nSMALL SUMMARY\n";

    std::cout
        << "cases="
        << small_cases
        << "\n";

    std::cout
        << "set_equality_pass="
        << small_pass
        << "/" << small_cases
        << "\n";

    run_large_structural();

    std::cout
        << "\nCOMPLEXITY\n";

    std::cout
        << "direct_scan = O(m log_p(m))\n";

    std::cout
        << "structural_enumeration = O(I * log_p(q))\n";

    std::cout
        << "I = product_i(q_i+1)-1\n";

    std::cout
        << "FINISHED EXPERIMENT 195\n";

    return 0;
}