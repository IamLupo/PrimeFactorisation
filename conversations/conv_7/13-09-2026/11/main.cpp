#include <algorithm>
#include <chrono>
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

struct PowerInfo {
    u64 value;
    int exponent;
    u64 residual;
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

bool lucas_hit(u64 p, u64 m, u64 t) {
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

u64 observed_start_gcd(
    u64 p,
    u64 m,
    std::size_t* interval_count
) {
    bool inside = false;
    u64 previous_start = 0;
    u64 g = 0;
    std::size_t count = 0;

    for (u64 t = 1; t <= m; ++t) {
        const bool hit = lucas_hit(p, m, t);

        if (hit && !inside) {
            const u64 start = t;

            if (count > 0) {
                g = std::gcd(
                    g,
                    start - previous_start
                );
            }

            previous_start = start;
            ++count;
            inside = true;
        }

        if (!hit && inside) {
            inside = false;
        }
    }

    if (interval_count != nullptr) {
        *interval_count = count;
    }

    return g;
}

PowerInfo analyze_power(u64 value, u64 p) {
    PowerInfo result{};

    result.value = value;
    result.exponent = 0;
    result.residual = value;

    if (value == 0) {
        return result;
    }

    while (result.residual % p == 0) {
        result.residual /= p;
        ++result.exponent;
    }

    return result;
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
    std::string result = "[";

    for (std::size_t i = 0;
         i < digits.size();
         ++i) {

        if (i != 0) {
            result += ",";
        }

        result += std::to_string(digits[i]);
    }

    result += "]";

    return result;
}

/*
    Conjectured exponent rule.

    Read base-p digits least-significant first.

    Each initial digit equal to p-1 contributes one level
    recursively.

    Once the first digit d < p-1 is reached, that digit
    contributes one level, followed by one level for every
    consecutive zero digit immediately above it.
*/
int predicted_exponent(
    u64 m,
    u64 p
) {
    int exponent = 0;

    while (m > 0) {
        const u64 digit = m % p;
        m /= p;

        if (digit == p - 1) {
            ++exponent;
            continue;
        }

        ++exponent;

        while (m > 0 &&
               m % p == 0) {

            ++exponent;
            m /= p;
        }

        return exponent;
    }

    return exponent;
}

void print_failure(
    u64 p,
    u64 m,
    u64 observed_g,
    const PowerInfo& observed,
    int predicted,
    std::size_t intervals
) {
    std::cout
        << "\nFAILURE\n"
        << "p=" << p
        << " m=" << m
        << " digits="
        << digits_string(
            base_p_digits(m, p)
        )
        << "\n"
        << "intervals=" << intervals
        << "\n"
        << "G=" << observed_g
        << "\n"
        << "observed_e="
        << observed.exponent
        << "\n"
        << "predicted_e="
        << predicted
        << "\n"
        << "residual="
        << observed.residual
        << "\n";
}

/*
    Small exhaustive validation.

    This is the expensive part, so the limits are deliberately
    modest and fixed.
*/
void run_small_exhaustive() {
    std::cout
        << "\nPHASE 1: SMALL EXHAUSTIVE VALIDATION\n";

    const u64 max_prime = 37;
    const u64 max_m = 1500;

    u64 total = 0;
    u64 passed = 0;
    u64 failed = 0;

    int max_e = 0;

    for (u64 p = 2; p <= max_prime; ++p) {
        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1; m <= max_m; ++m) {
            std::size_t intervals = 0;

            const u64 g =
                observed_start_gcd(
                    p,
                    m,
                    &intervals
                );

            if (intervals < 2) {
                continue;
            }

            ++total;

            const PowerInfo observed =
                analyze_power(g, p);

            const int predicted =
                predicted_exponent(m, p);

            max_e =
                std::max(
                    max_e,
                    observed.exponent
                );

            const bool ok =
                observed.residual == 1 &&
                observed.exponent ==
                    predicted;

            if (ok) {
                ++passed;
            } else {
                ++failed;

                if (failed <= 20) {
                    print_failure(
                        p,
                        m,
                        g,
                        observed,
                        predicted,
                        intervals
                    );
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "multi_interval_cases="
        << total
        << "\n"
        << "prediction_pass="
        << passed
        << "/" << total
        << "\n"
        << "prediction_fail="
        << failed
        << "\n"
        << "maximum_observed_e="
        << max_e
        << "\n";
}

/*
    Deterministic targeted cases.

    These deliberately exercise the structures that broke
    Experiment 165.
*/
void run_targeted_cases() {
    std::cout
        << "\nPHASE 2: TARGETED VALIDATION\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 6},
        {2, 12},
        {2, 13},
        {2, 20},
        {2, 24},
        {2, 25},
        {2, 27},
        {2, 28},
        {2, 40},
        {2, 48},
        {2, 96},
        {2, 192},
        {2, 384},
        {2, 768},

        {3, 22},
        {3, 31},
        {3, 41},
        {3, 49},
        {3, 59},
        {3, 68},
        {3, 85},
        {3, 95},
        {3, 122},
        {3, 125},
        {3, 206},
        {3, 377},

        {5, 59},
        {5, 64},
        {5, 69},
        {5, 94},
        {5, 159},
        {5, 194},
        {5, 219},

        {7, 98},
        {11, 242},
        {17, 578}
    };

    u64 passed = 0;
    u64 tested = 0;

    for (const auto& c : cases) {
        std::size_t intervals = 0;

        const u64 g =
            observed_start_gcd(
                c.p,
                c.m,
                &intervals
            );

        if (intervals < 2) {
            continue;
        }

        ++tested;

        const PowerInfo observed =
            analyze_power(g, c.p);

        const int predicted =
            predicted_exponent(
                c.m,
                c.p
            );

        const bool ok =
            observed.residual == 1 &&
            observed.exponent ==
                predicted;

        if (ok) {
            ++passed;
        }

        std::cout
            << "p=" << c.p
            << " m=" << c.m
            << " digits="
            << digits_string(
                base_p_digits(c.m, c.p)
            )
            << " intervals="
            << intervals
            << " G="
            << g
            << " observed_e="
            << observed.exponent
            << " predicted_e="
            << predicted
            << " "
            << (ok ? "PASS" : "FAIL")
            << "\n";
    }

    std::cout
        << "\nSUMMARY PHASE 2\n"
        << "tested="
        << tested
        << "\n"
        << "passed="
        << passed
        << "/" << tested
        << "\n";
}

/*
    Deterministic pseudo-random generator.
*/
u64 next_random(u64& state) {
    state =
        state * 6364136223846793005ULL +
        1442695040888963407ULL;

    return state;
}

/*
    Generate a deterministic random m.

    We keep m <= 100000 so the actual interval scan remains
    cheap even in the validation phase.
*/
u64 random_m(
    u64& state,
    u64 limit
) {
    return 1 + (next_random(state) % limit);
}

/*
    Larger validation.

    Only 100 cases are actually evaluated with the expensive
    interval scan.

    This gives an empirical check at much larger m while
    keeping the runtime bounded.
*/
void run_random_large_validation() {
    std::cout
        << "\nPHASE 3: RANDOM LARGE VALIDATION\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67, 101,
        127, 211, 431, 1009
    };

    constexpr int CASES = 100;
    constexpr u64 MAX_M = 100000;

    u64 state =
        0x166B2026ULL;

    u64 tested = 0;
    u64 passed = 0;
    u64 failed = 0;

    int max_e = 0;

    for (int i = 0;
         i < CASES;
         ++i) {

        const u64 p =
            primes[
                next_random(state) %
                primes.size()
            ];

        const u64 m =
            random_m(
                state,
                MAX_M
            );

        std::size_t intervals = 0;

        const u64 g =
            observed_start_gcd(
                p,
                m,
                &intervals
            );

        if (intervals < 2) {
            continue;
        }

        ++tested;

        const PowerInfo observed =
            analyze_power(g, p);

        const int predicted =
            predicted_exponent(m, p);

        max_e =
            std::max(
                max_e,
                observed.exponent
            );

        if (observed.residual == 1 &&
            observed.exponent ==
                predicted) {

            ++passed;
        } else {
            ++failed;

            if (failed <= 10) {
                print_failure(
                    p,
                    m,
                    g,
                    observed,
                    predicted,
                    intervals
                );
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 3\n"
        << "attempted="
        << CASES
        << "\n"
        << "multi_interval_tested="
        << tested
        << "\n"
        << "prediction_pass="
        << passed
        << "/" << tested
        << "\n"
        << "prediction_fail="
        << failed
        << "\n"
        << "maximum_observed_e="
        << max_e
        << "\n";
}

/*
    Pure digit-side search.

    This phase performs no interval construction at all.

    It searches millions of m values for unusually large
    predicted exponents. These results are not proof; they
    tell us which digit structures are worth validating.
*/
void run_digit_search() {
    std::cout
        << "\nPHASE 4: FAST DIGIT SEARCH\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31, 37
    };

    constexpr u64 LIMIT = 5000000;

    int global_max_e = 0;
    u64 global_max_p = 0;
    u64 global_max_m = 0;

    for (u64 p : primes) {
        int max_e = 0;
        u64 max_m = 0;

        for (u64 m = 1;
             m <= LIMIT;
             ++m) {

            const int e =
                predicted_exponent(
                    m,
                    p
                );

            if (e > max_e) {
                max_e = e;
                max_m = m;
            }
        }

        std::cout
            << "p=" << p
            << " predicted_max_e="
            << max_e
            << " at_m="
            << max_m
            << " digits="
            << digits_string(
                base_p_digits(
                    max_m,
                    p
                )
            )
            << "\n";

        if (max_e > global_max_e) {
            global_max_e = max_e;
            global_max_p = p;
            global_max_m = max_m;
        }
    }

    std::cout
        << "\nGLOBAL DIGIT SEARCH MAX\n"
        << "p="
        << global_max_p
        << "\n"
        << "m="
        << global_max_m
        << "\n"
        << "predicted_e="
        << global_max_e
        << "\n"
        << "digits="
        << digits_string(
            base_p_digits(
                global_max_m,
                global_max_p
            )
        )
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 166;

    const auto start =
        std::chrono::steady_clock::now();

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "B\n";

    run_small_exhaustive();
    run_targeted_cases();
    run_random_large_validation();
    run_digit_search();

    const auto finish =
        std::chrono::steady_clock::now();

    const double seconds =
        std::chrono::duration<double>(
            finish - start
        ).count();

    std::cout
        << "\nRUNTIME_SECONDS="
        << seconds
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "B\n";

    return 0;
}
