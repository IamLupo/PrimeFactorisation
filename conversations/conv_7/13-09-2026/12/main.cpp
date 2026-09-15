#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 lo;
    u64 hi;
};

struct Structure {
    int trailing_pminus1;
    int zero_run;
    int exponent;
    u64 first_start;
    u64 modulus;
};

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d * d <= n; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

bool lucas_hit(
    u64 p,
    u64 m,
    u64 t
) {
    while (m > 0 || t > 0) {
        const u64 md = m % p;
        const u64 td = t % p;

        if (td > md) {
            return true;
        }

        m /= p;
        t /= p;
    }

    return false;
}

/*
    Extract the digit structure predicted by Experiment 166B.

    Digits are least-significant first.

    a = number of initial digits equal to p-1.

    b = first digit below p-1.

    z = number of consecutive zero digits after b.

    e = a + 1 + z

    first_start = (b+1) * p^a

    modulus = p^e
*/
Structure analyze_digit_structure(
    u64 p,
    u64 m
) {
    Structure result{};

    u64 x = m;

    while (x > 0 &&
           x % p == p - 1) {

        ++result.trailing_pminus1;
        x /= p;
    }

    /*
        The all-(p-1) case does not normally produce
        a useful multi-interval structure. Keep the
        formula defined anyway.
    */
    if (x == 0) {
        result.zero_run = 0;
        result.exponent =
            result.trailing_pminus1;

        u64 pa = 1;

        for (int i = 0;
             i < result.trailing_pminus1;
             ++i) {
            pa *= p;
        }

        result.first_start = pa;
        result.modulus = pa;
        return result;
    }

    const u64 b = x % p;

    u64 pa = 1;

    for (int i = 0;
         i < result.trailing_pminus1;
         ++i) {
        pa *= p;
    }

    x /= p;

    while (x > 0 &&
           x % p == 0) {

        ++result.zero_run;
        x /= p;
    }

    result.exponent =
        result.trailing_pminus1 +
        1 +
        result.zero_run;

    result.first_start =
        (b + 1) * pa;

    result.modulus = 1;

    for (int i = 0;
         i < result.exponent;
         ++i) {
        result.modulus *= p;
    }

    return result;
}

std::vector<Interval> build_intervals(
    u64 p,
    u64 m
) {
    std::vector<Interval> intervals;

    bool inside = false;
    u64 start = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit =
            lucas_hit(p, m, t);

        if (hit && !inside) {
            start = t;
            inside = true;
        }

        if (!hit && inside) {
            intervals.push_back({
                start,
                t - 1
            });

            inside = false;
        }
    }

    if (inside) {
        intervals.push_back({
            start,
            m
        });
    }

    return intervals;
}

u64 gcd_start_gaps(
    const std::vector<Interval>& intervals
) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 g = 0;

    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        g = std::gcd(
            g,
            intervals[i].lo -
            intervals[i - 1].lo
        );
    }

    return g;
}

u64 minimum_start_gap(
    const std::vector<Interval>& intervals
) {
    if (intervals.size() < 2) {
        return 0;
    }

    u64 minimum =
        intervals[1].lo -
        intervals[0].lo;

    for (std::size_t i = 2;
         i < intervals.size();
         ++i) {

        minimum =
            std::min(
                minimum,
                intervals[i].lo -
                intervals[i - 1].lo
            );
    }

    return minimum;
}

std::vector<int> base_p_digits(
    u64 m,
    u64 p
) {
    std::vector<int> digits;

    while (m > 0) {
        digits.push_back(
            static_cast<int>(m % p)
        );

        m /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

std::string digits_string(
    const std::vector<int>& digits
) {
    std::string s = "[";

    for (std::size_t i = 0;
         i < digits.size();
         ++i) {

        if (i != 0) {
            s += ",";
        }

        s += std::to_string(digits[i]);
    }

    s += "]";

    return s;
}

void print_failure(
    const std::string& reason,
    u64 p,
    u64 m,
    const Structure& structure,
    const std::vector<Interval>& intervals
) {
    std::cout
        << "\nFAILURE: "
        << reason
        << "\n";

    std::cout
        << "p=" << p
        << " m=" << m
        << " digits="
        << digits_string(
            base_p_digits(m, p)
        )
        << "\n";

    std::cout
        << "a="
        << structure.trailing_pminus1
        << " z="
        << structure.zero_run
        << " e="
        << structure.exponent
        << "\n";

    std::cout
        << "predicted_first_start="
        << structure.first_start
        << "\n";

    std::cout
        << "predicted_modulus="
        << structure.modulus
        << "\n";

    std::cout
        << "actual_interval_count="
        << intervals.size()
        << "\n";

    std::cout
        << "actual_starts=";

    const std::size_t limit =
        std::min<std::size_t>(
            intervals.size(),
            20
        );

    for (std::size_t i = 0;
         i < limit;
         ++i) {

        if (i != 0) {
            std::cout << ",";
        }

        std::cout
            << intervals[i].lo;
    }

    if (intervals.size() > limit) {
        std::cout << ",...";
    }

    std::cout << "\n";
}

void run_exhaustive_phase() {
    std::cout
        << "\nPHASE 1: EXHAUSTIVE START-STRUCTURE TEST\n";

    u64 total = 0;
    u64 first_start_pass = 0;
    u64 residue_pass = 0;
    u64 minimum_gap_pass = 0;
    u64 gcd_pass = 0;

    u64 failure_first_start = 0;
    u64 failure_residue = 0;
    u64 failure_gap = 0;
    u64 failure_gcd = 0;

    for (u64 p = 2; p <= 37; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= 1500; ++m) {
            const auto intervals =
                build_intervals(p, m);

            if (intervals.size() < 2) {
                continue;
            }

            ++total;

            const Structure structure =
                analyze_digit_structure(
                    p,
                    m
                );

            /*
                Test exact first interval start.
            */
            const bool first_ok =
                intervals[0].lo ==
                structure.first_start;

            if (first_ok) {
                ++first_start_pass;
            } else {
                ++failure_first_start;

                if (failure_first_start <= 10) {
                    print_failure(
                        "first start",
                        p,
                        m,
                        structure,
                        intervals
                    );
                }
            }

            /*
                Every interval start should have the
                same residue modulo p^e.
            */
            bool residue_ok = true;

            for (const auto& interval :
                 intervals) {

                if (interval.lo %
                    structure.modulus !=
                    structure.first_start %
                    structure.modulus) {

                    residue_ok = false;
                    break;
                }
            }

            if (residue_ok) {
                ++residue_pass;
            } else {
                ++failure_residue;

                if (failure_residue <= 10) {
                    print_failure(
                        "start residue",
                        p,
                        m,
                        structure,
                        intervals
                    );
                }
            }

            /*
                At least one adjacent start gap should be
                exactly p^e.
            */
            const u64 minimum_gap =
                minimum_start_gap(
                    intervals
                );

            if (minimum_gap ==
                structure.modulus) {

                ++minimum_gap_pass;
            } else {
                ++failure_gap;

                if (failure_gap <= 10) {
                    print_failure(
                        "minimum gap",
                        p,
                        m,
                        structure,
                        intervals
                    );

                    std::cout
                        << "actual_minimum_gap="
                        << minimum_gap
                        << "\n";
                }
            }

            /*
                Existing gcd statement.
            */
            const u64 g =
                gcd_start_gaps(
                    intervals
                );

            if (g == structure.modulus) {
                ++gcd_pass;
            } else {
                ++failure_gcd;

                if (failure_gcd <= 10) {
                    print_failure(
                        "GCD",
                        p,
                        m,
                        structure,
                        intervals
                    );

                    std::cout
                        << "actual_gcd="
                        << g
                        << "\n";
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n";

    std::cout
        << "multi_interval_cases="
        << total
        << "\n";

    std::cout
        << "first_start_pass="
        << first_start_pass
        << "/" << total
        << "\n";

    std::cout
        << "start_residue_pass="
        << residue_pass
        << "/" << total
        << "\n";

    std::cout
        << "minimum_gap_pass="
        << minimum_gap_pass
        << "/" << total
        << "\n";

    std::cout
        << "gcd_pass="
        << gcd_pass
        << "/" << total
        << "\n";

    std::cout
        << "failure_first_start="
        << failure_first_start
        << "\n";

    std::cout
        << "failure_residue="
        << failure_residue
        << "\n";

    std::cout
        << "failure_gap="
        << failure_gap
        << "\n";

    std::cout
        << "failure_gcd="
        << failure_gcd
        << "\n";
}

void run_targeted_phase() {
    std::cout
        << "\nPHASE 2: TARGETED LARGE STRUCTURES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 12},
        {2, 24},
        {2, 48},
        {2, 96},
        {2, 192},
        {2, 384},
        {2, 768},
        {2, 1536},

        {3, 41},
        {3, 59},
        {3, 125},
        {3, 377},
        {3, 1133},
        {3, 3401},

        {5, 59},
        {5, 94},
        {5, 194},
        {5, 219},

        {7, 98},
        {11, 242},
        {17, 578}
    };

    for (const auto& c : cases) {
        const auto intervals =
            build_intervals(
                c.p,
                c.m
            );

        if (intervals.size() < 2) {
            std::cout
                << "p=" << c.p
                << " m=" << c.m
                << " insufficient_intervals\n";

            continue;
        }

        const Structure structure =
            analyze_digit_structure(
                c.p,
                c.m
            );

        const u64 g =
            gcd_start_gaps(intervals);

        const u64 min_gap =
            minimum_start_gap(intervals);

        bool residue_ok = true;

        for (const auto& interval :
             intervals) {

            if (interval.lo %
                structure.modulus !=
                structure.first_start %
                structure.modulus) {

                residue_ok = false;
                break;
            }
        }

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " digits="
            << digits_string(
                base_p_digits(
                    c.m,
                    c.p
                )
            )
            << " a="
            << structure.trailing_pminus1
            << " z="
            << structure.zero_run
            << " e="
            << structure.exponent
            << " first_start="
            << intervals[0].lo
            << "/"
            << structure.first_start
            << " residue="
            << (residue_ok ? "PASS" : "FAIL")
            << " min_gap="
            << min_gap
            << " G="
            << g
            << "\n";
    }
}

void run_random_phase() {
    std::cout
        << "\nPHASE 3: RANDOM VALIDATION\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127
    };

    u64 state =
        0x1672026ULL;

    constexpr int CASES = 100;
    constexpr u64 MAX_M = 100000;

    u64 tested = 0;
    u64 failed = 0;

    for (int i = 0;
         i < CASES;
         ++i) {

        const u64 p =
            primes[
                state %
                primes.size()
            ];

        state =
            state * 6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 m =
            1 + (state % MAX_M);

        state =
            state * 6364136223846793005ULL +
            1442695040888963407ULL;

        const auto intervals =
            build_intervals(p, m);

        if (intervals.size() < 2) {
            continue;
        }

        ++tested;

        const Structure structure =
            analyze_digit_structure(
                p,
                m
            );

        const u64 g =
            gcd_start_gaps(intervals);

        const u64 minimum_gap =
            minimum_start_gap(intervals);

        bool residue_ok = true;

        for (const auto& interval :
             intervals) {

            if (interval.lo %
                structure.modulus !=
                structure.first_start %
                structure.modulus) {

                residue_ok = false;
                break;
            }
        }

        const bool ok =
            intervals[0].lo ==
                structure.first_start &&
            residue_ok &&
            minimum_gap ==
                structure.modulus &&
            g ==
                structure.modulus;

        if (!ok) {
            ++failed;

            if (failed <= 10) {
                print_failure(
                    "random validation",
                    p,
                    m,
                    structure,
                    intervals
                );

                std::cout
                    << "actual_min_gap="
                    << minimum_gap
                    << "\n"
                    << "actual_gcd="
                    << g
                    << "\n";
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n";

    std::cout
        << "multi_interval_tested="
        << tested
        << "\n";

    std::cout
        << "all_properties_pass="
        << (failed == 0 ? "YES" : "NO")
        << "\n";

    std::cout
        << "failures="
        << failed
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 167;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    run_exhaustive_phase();
    run_targeted_phase();
    run_random_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
