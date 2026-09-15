#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

/*
    EXPERIMENT 284

    Goal:

        1. exact global MISS delta
        2. explicit brute-force new-MISS slice
        3. proposed closed local-slice formula

    We want to determine which equality fails:

        (1) == (2)
        (2) == (3)

    The brute-force slice is defined directly from the
    digitwise condition and does NOT call miss_prefix().
*/

std::string to_string_u128(u128 x) {
    if (x == 0) {
        return "0";
    }

    std::string s;

    while (x > 0) {
        const unsigned digit =
            static_cast<unsigned>(x % 10);

        s.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    return s;
}

u128 ipow(u128 base, unsigned exp) {
    u128 result = 1;

    while (exp > 0) {
        if (exp & 1u) {
            result *= base;
        }

        base *= base;
        exp >>= 1u;
    }

    return result;
}

bool digitwise_leq(
    u128 x,
    u128 m,
    u64 p
) {
    while (x > 0 || m > 0) {
        const u64 xd =
            static_cast<u64>(
                x % static_cast<u128>(p)
            );

        const u64 md =
            static_cast<u64>(
                m % static_cast<u128>(p)
            );

        if (xd > md) {
            return false;
        }

        x /= static_cast<u128>(p);
        m /= static_cast<u128>(p);
    }

    return true;
}

/*
    Direct brute-force MISS prefix.

    This is intentionally independent from the closed
    miss_prefix implementation.
*/
u64 brute_miss_prefix(
    u64 m,
    u64 y,
    u64 p
) {
    u64 count = 0;

    for (u64 x = 0; x <= y; ++x) {
        if (digitwise_leq(x, m, p)) {
            ++count;
        }
    }

    return count;
}

/*
    Exact global local delta:

        M_{m+p^r}(n-1) - M_m(n-1)

    calculated entirely by brute force.
*/
u64 brute_global_delta(
    u64 current,
    unsigned r,
    u64 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    const u64 step =
        static_cast<u64>(
            ipow(
                static_cast<u128>(p),
                r
            )
        );

    const u64 next =
        current + step;

    const u64 before =
        brute_miss_prefix(
            current,
            n - 1,
            p
        );

    const u64 after =
        brute_miss_prefix(
            next,
            n - 1,
            p
        );

    return after - before;
}

/*
    Independent direct definition of the NEW MISS slice.

    For a no-carry update m' = m + p^r:

        x_r = m_r + 1

    and every other digit must satisfy the old bound.

    We simply enumerate x < n and check that directly.
*/
u64 brute_new_slice(
    u64 current,
    unsigned r,
    u64 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    const u64 step =
        static_cast<u64>(
            ipow(
                static_cast<u128>(p),
                r
            )
        );

    const u64 next =
        current + step;

    u64 m_temp = current;

    u64 mr =
        static_cast<u64>(
            m_temp /
            step
        ) % p;

    const u64 new_digit = mr + 1;

    u64 count = 0;

    for (u64 x = 0; x < n; ++x) {
        if (!digitwise_leq(x, next, p)) {
            continue;
        }

        if (digitwise_leq(x, current, p)) {
            continue;
        }

        /*
            Explicitly verify the defining changed digit.
        */
        u64 xd =
            static_cast<u64>(
                x / step
            ) % p;

        if (xd != new_digit) {
            std::cerr
                << "INTERNAL SLICE DEFINITION ERROR\n";

            std::abort();
        }

        ++count;
    }

    return count;
}

/*
    Candidate closed formula under investigation.
*/
u128 closed_local_slice(
    u64 current,
    unsigned r,
    u64 n,
    u64 p
) {
    if (n == 0) {
        return 0;
    }

    const u64 y = n - 1;

    const u128 pr =
        ipow(
            static_cast<u128>(p),
            r
        );

    const u128 pr1 =
        pr *
        static_cast<u128>(p);

    const u128 high_m =
        static_cast<u128>(current) /
        pr1;

    const u128 high_y =
        static_cast<u128>(y) /
        pr1;

    const u64 mr =
        static_cast<u64>(
            (static_cast<u128>(current) / pr) %
            static_cast<u128>(p)
        );

    const u64 digit_y =
        static_cast<u64>(
            (static_cast<u128>(y) / pr) %
            static_cast<u128>(p)
        );

    const u128 low_m =
        static_cast<u128>(current) %
        pr;

    const u128 low_y =
        static_cast<u128>(y) %
        pr;

    const u64 new_digit =
        mr + 1;

    /*
        Invalid carry case.
    */
    if (new_digit >= p) {
        return 0;
    }

    /*
        We intentionally use a recursive independent
        brute count for the lower block here only for
        small values.
    */
    u128 lower_total = 0;

    for (u128 x = 0; x < pr; ++x) {
        if (digitwise_leq(x, low_m, p)) {
            ++lower_total;
        }
    }

    u128 result = 0;

    /*
        Count admissible high prefixes strictly below high_y.
    */
    if (high_y > 0) {
        for (u128 h = 0; h < high_y; ++h) {
            if (digitwise_leq(h, high_m, p)) {
                result += lower_total;
            }
        }
    }

    /*
        Equal high prefix.
    */
    if (!digitwise_leq(
            high_y,
            high_m,
            p
        )) {
        return result;
    }

    if (digit_y < new_digit) {
        return result;
    }

    if (digit_y > new_digit) {
        return result + lower_total;
    }

    /*
        Equal changed digit: partial lower block.
    */
    for (u128 x = 0; x <= low_y; ++x) {
        if (digitwise_leq(x, low_m, p)) {
            ++result;
        }
    }

    return result;
}

struct Stats {
    std::size_t cases = 0;

    std::size_t global_vs_slice_fail = 0;
    std::size_t slice_vs_formula_fail = 0;
    std::size_t global_vs_formula_fail = 0;

    std::size_t printed = 0;
};

void print_failure(
    const char* category,
    u64 p,
    u64 current,
    unsigned r,
    u64 n,
    u64 global_delta,
    u64 brute_slice,
    u128 formula
) {
    std::cout
        << "\nFAIL "
        << category
        << '\n';

    std::cout
        << "p="
        << p
        << '\n';

    std::cout
        << "current="
        << current
        << '\n';

    std::cout
        << "r="
        << r
        << '\n';

    std::cout
        << "n="
        << n
        << '\n';

    std::cout
        << "global_delta="
        << global_delta
        << '\n';

    std::cout
        << "brute_slice="
        << brute_slice
        << '\n';

    std::cout
        << "formula="
        << to_string_u128(formula)
        << '\n';
}

void check_case(
    u64 p,
    u64 current,
    unsigned r,
    u64 n,
    Stats& stats
) {
    const u64 step =
        static_cast<u64>(
            ipow(
                static_cast<u128>(p),
                r
            )
        );

    const u64 mr =
        static_cast<u64>(
            (static_cast<u128>(current) / step) %
            static_cast<u128>(p)
        );

    /*
        Only test no-carry updates.
    */
    if (mr + 1 >= p) {
        return;
    }

    ++stats.cases;

    const u64 global_delta =
        brute_global_delta(
            current,
            r,
            n,
            p
        );

    const u64 brute_slice =
        brute_new_slice(
            current,
            r,
            n,
            p
        );

    const u128 formula =
        closed_local_slice(
            current,
            r,
            n,
            p
        );

    const bool global_ok =
        global_delta == brute_slice;

    const bool formula_ok =
        static_cast<u128>(brute_slice) == formula;

    const bool all_ok =
        static_cast<u128>(global_delta) == formula;

    if (!global_ok) {
        ++stats.global_vs_slice_fail;

        if (stats.printed < 5) {
            ++stats.printed;

            print_failure(
                "GLOBAL_VS_BRUTE_SLICE",
                p,
                current,
                r,
                n,
                global_delta,
                brute_slice,
                formula
            );
        }
    }

    if (!formula_ok) {
        ++stats.slice_vs_formula_fail;

        if (stats.printed < 5) {
            ++stats.printed;

            print_failure(
                "BRUTE_SLICE_VS_FORMULA",
                p,
                current,
                r,
                n,
                global_delta,
                brute_slice,
                formula
            );
        }
    }

    if (!all_ok) {
        ++stats.global_vs_formula_fail;
    }
}

void deterministic_tests(
    Stats& stats
) {
    const u64 primes[] = {
        2, 3, 5, 7, 11
    };

    for (u64 p : primes) {
        for (u64 current = 0;
             current <= 80;
             ++current) {

            for (unsigned r = 0;
                 r <= 4;
                 ++r) {

                for (u64 n = 0;
                     n <= 120;
                     ++n) {

                    check_case(
                        p,
                        current,
                        r,
                        n,
                        stats
                    );
                }
            }
        }
    }
}

void boundary_tests(
    Stats& stats
) {
    const u64 primes[] = {
        2, 3, 5, 7, 11, 13,
        17, 19, 23
    };

    for (u64 p : primes) {
        for (unsigned r = 0;
             r <= 5;
             ++r) {

            const u64 step =
                static_cast<u64>(
                    ipow(
                        static_cast<u128>(p),
                        r
                    )
                );

            const u64 currents[] = {
                0,
                1,
                step > 0 ? step - 1 : 0,
                step,
                step + 1,
                2 * step + 3
            };

            for (u64 current : currents) {
                const u64 mr =
                    static_cast<u64>(
                        (static_cast<u128>(current) / step) %
                        static_cast<u128>(p)
                    );

                if (mr + 1 >= p) {
                    continue;
                }

                const u64 next =
                    current + step;

                const u64 ns[] = {
                    0,
                    1,
                    current,
                    current + 1,
                    next,
                    next > 0 ? next - 1 : 0,
                    next + 1,
                    step,
                    step + 1
                };

                for (u64 n : ns) {
                    if (n > 3000) {
                        continue;
                    }

                    check_case(
                        p,
                        current,
                        r,
                        n,
                        stats
                    );
                }
            }
        }
    }
}

void random_tests(
    Stats& stats,
    std::size_t count
) {
    std::mt19937_64 rng(
        0x284284ULL
    );

    /*
        Keep n small enough that the independent
        brute-force definition is cheap.

        This is intentional: we are diagnosing the
        mathematics, not scalability.
    */
    const u64 primes[] = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    for (std::size_t i = 0;
         i < count;
         ++i) {

        const u64 p =
            primes[
                rng() %
                (sizeof(primes) /
                 sizeof(primes[0]))
            ];

        const unsigned r =
            static_cast<unsigned>(
                rng() % 5
            );

        /*
            Keep current sufficiently small so the
            independent closed-form diagnostic itself
            remains practical.
        */
        const u64 current =
            rng() % 5000;

        const u64 n =
            rng() % 3000;

        check_case(
            p,
            current,
            r,
            n,
            stats
        );
    }
}

void print_stats(
    const char* name,
    const Stats& stats
) {
    std::cout
        << '\n'
        << name
        << '\n';

    std::cout
        << "cases="
        << stats.cases
        << '\n';

    std::cout
        << "global_vs_brute_slice_fail="
        << stats.global_vs_slice_fail
        << '\n';

    std::cout
        << "brute_slice_vs_formula_fail="
        << stats.slice_vs_formula_fail
        << '\n';

    std::cout
        << "global_vs_formula_fail="
        << stats.global_vs_formula_fail
        << '\n';
}

int main() {
    std::cout
        << "START EXPERIMENT 284\n";

    std::cout
        << "INDEPENDENT BRUTE-FORCE NEW-MISS SLICE\n";

    std::cout
        << "GLOBAL DELTA VS SLICE VS CLOSED FORM\n\n";

    Stats deterministic;
    Stats boundary;
    Stats random;

    deterministic_tests(
        deterministic
    );

    print_stats(
        "DETERMINISTIC",
        deterministic
    );

    boundary_tests(
        boundary
    );

    print_stats(
        "BOUNDARY",
        boundary
    );

    random_tests(
        random,
        5000
    );

    print_stats(
        "RANDOM",
        random
    );

    const std::size_t total_cases =
        deterministic.cases +
        boundary.cases +
        random.cases;

    const std::size_t global_slice_fail =
        deterministic.global_vs_slice_fail +
        boundary.global_vs_slice_fail +
        random.global_vs_slice_fail;

    const std::size_t slice_formula_fail =
        deterministic.slice_vs_formula_fail +
        boundary.slice_vs_formula_fail +
        random.slice_vs_formula_fail;

    const std::size_t global_formula_fail =
        deterministic.global_vs_formula_fail +
        boundary.global_vs_formula_fail +
        random.global_vs_formula_fail;

    std::cout
        << "\nTOTAL\n";

    std::cout
        << "cases="
        << total_cases
        << '\n';

    std::cout
        << "global_vs_brute_slice_fail="
        << global_slice_fail
        << '\n';

    std::cout
        << "brute_slice_vs_formula_fail="
        << slice_formula_fail
        << '\n';

    std::cout
        << "global_vs_formula_fail="
        << global_formula_fail
        << '\n';

    const bool pass =
        global_slice_fail == 0 &&
        slice_formula_fail == 0 &&
        global_formula_fail == 0;

    std::cout
        << "\nOVERALL PASS="
        << (pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 284\n";

    return pass ? 0 : 1;
}
