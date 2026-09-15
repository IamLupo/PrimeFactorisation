#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Structure {
    int a;
    int z;
    int exponent;

    u64 first_start;
    u64 modulus;
    u64 second_start;
};

struct ValidationResult {
    bool first_is_hit;
    bool first_previous_is_miss;
    bool second_is_hit;
    bool second_previous_is_miss;
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
    Recover the experimentally observed structure:

        a = number of initial digits equal to p-1

        b = first digit below p-1

        z = number of zero digits after b

        e = a + 1 + z

        s0 = (b+1) p^a

        modulus = p^e

        s1 = s0 + p^e
*/
Structure analyze_structure(
    u64 p,
    u64 m
) {
    Structure s{};

    u64 x = m;

    while (x > 0 &&
           x % p == p - 1) {

        ++s.a;
        x /= p;
    }

    /*
        The all-(p-1) case normally does not give a useful
        multi-interval configuration. Keep it defined.
    */
    if (x == 0) {
        s.z = 0;
        s.exponent = std::max(1, s.a);

        s.first_start = 1;
        s.modulus = 1;

        for (int i = 0;
             i < s.exponent;
             ++i) {

            s.modulus *= p;
        }

        s.second_start =
            s.first_start +
            s.modulus;

        return s;
    }

    const u64 b = x % p;

    u64 pa = 1;

    for (int i = 0;
         i < s.a;
         ++i) {

        pa *= p;
    }

    x /= p;

    while (x > 0 &&
           x % p == 0) {

        ++s.z;
        x /= p;
    }

    s.exponent =
        s.a + 1 + s.z;

    s.first_start =
        (b + 1) * pa;

    s.modulus = 1;

    for (int i = 0;
         i < s.exponent;
         ++i) {

        s.modulus *= p;
    }

    s.second_start =
        s.first_start +
        s.modulus;

    return s;
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

        result +=
            std::to_string(digits[i]);
    }

    result += "]";

    return result;
}

/*
    Validate the predicted starts without searching.

    We check:

        s0 - 1 : MISS
        s0     : HIT
        s1 - 1 : MISS
        s1     : HIT

    This does not prove that there are no earlier HITs between
    these boundaries. That is tested separately on a small
    exhaustive sample.
*/
ValidationResult validate_predicted_boundaries(
    u64 p,
    u64 m,
    const Structure& s
) {
    ValidationResult result{};

    result.first_previous_is_miss =
        s.first_start > 0 &&
        !lucas_hit(
            p,
            m,
            s.first_start - 1
        );

    result.first_is_hit =
        lucas_hit(
            p,
            m,
            s.first_start
        );

    result.second_previous_is_miss =
        s.second_start > 0 &&
        !lucas_hit(
            p,
            m,
            s.second_start - 1
        );

    result.second_is_hit =
        lucas_hit(
            p,
            m,
            s.second_start
        );

    return result;
}

/*
    Small exhaustive interval scan.

    This is deliberately limited to small m.
    It verifies the stronger claim that the predicted two
    starts really are the first two HIT starts.
*/
bool validate_first_two_exactly(
    u64 p,
    u64 m,
    const Structure& s,
    u64 scan_limit
) {
    if (m > scan_limit) {
        return true;
    }

    int found = 0;

    for (u64 t = 1;
         t <= m;
         ++t) {

        const bool hit =
            lucas_hit(
                p,
                m,
                t
            );

        if (!hit) {
            continue;
        }

        ++found;

        if (found == 1 &&
            t != s.first_start) {

            return false;
        }

        if (found == 2 &&
            t != s.second_start) {

            return false;
        }

        if (found >= 2) {
            return true;
        }
    }

    return false;
}

void print_failure(
    const std::string& reason,
    u64 p,
    u64 m,
    const Structure& s,
    const ValidationResult& v
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
        << "a=" << s.a
        << " z=" << s.z
        << " e=" << s.exponent
        << "\n";

    std::cout
        << "s0=" << s.first_start
        << " s1=" << s.second_start
        << " gap=" << s.modulus
        << "\n";

    std::cout
        << "s0-1_miss="
        << v.first_previous_is_miss
        << " s0_hit="
        << v.first_is_hit
        << " s1-1_miss="
        << v.second_previous_is_miss
        << " s1_hit="
        << v.second_is_hit
        << "\n";
}

void run_small_exhaustive() {
    std::cout
        << "\nPHASE 1: SMALL EXACT FIRST-TWO TEST\n";

    constexpr u64 MAX_PRIME = 37;
    constexpr u64 MAX_M = 1500;

    u64 total = 0;
    u64 pass = 0;
    u64 fail = 0;

    for (u64 p = 2;
         p <= MAX_PRIME;
         ++p) {

        if (!is_prime(p)) {
            continue;
        }

        for (u64 m = 1;
             m <= MAX_M;
             ++m) {

            const Structure s =
                analyze_structure(
                    p,
                    m
                );

            /*
                Skip configurations that do not have enough
                room for the second predicted boundary.
            */
            if (s.second_start > m) {
                continue;
            }

            const ValidationResult v =
                validate_predicted_boundaries(
                    p,
                    m,
                    s
                );

            const bool boundary_ok =
                v.first_previous_is_miss &&
                v.first_is_hit &&
                v.second_previous_is_miss &&
                v.second_is_hit;

            if (!boundary_ok) {
                continue;
            }

            ++total;

            const bool exact =
                validate_first_two_exactly(
                    p,
                    m,
                    s,
                    MAX_M
                );

            if (exact) {
                ++pass;
            } else {
                ++fail;

                if (fail <= 20) {
                    print_failure(
                        "exact first-two mismatch",
                        p,
                        m,
                        s,
                        v
                    );
                }
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 1\n"
        << "boundary_cases="
        << total
        << "\n"
        << "exact_first_two_pass="
        << pass
        << "/" << total
        << "\n"
        << "failures="
        << fail
        << "\n";
}

/*
    Large validation.

    No search over t.

    We only evaluate the four relevant Lucas predicates
    around the predicted boundaries.

    This makes runtime O(number of cases * log_p(m)).
*/
void run_large_boundary_phase() {
    std::cout
        << "\nPHASE 2: LARGE BOUNDARY VALIDATION\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7,
        11, 13, 17, 19,
        23, 29, 31, 37,
        43, 53, 67,
        101, 127, 211,
        431, 1009, 2003,
        4001, 5003
    };

    constexpr int CASES = 10000;
    constexpr u64 MAX_M =
        1000000000000ULL;

    u64 state =
        0x168B2026ULL;

    u64 tested = 0;
    u64 passed = 0;
    u64 failed = 0;

    for (int i = 0;
         i < CASES;
         ++i) {

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 p =
            primes[
                state % primes.size()
            ];

        state =
            state *
            6364136223846793005ULL +
            1442695040888963407ULL;

        const u64 m =
            1 + state % MAX_M;

        const Structure s =
            analyze_structure(
                p,
                m
            );

        if (s.second_start > m) {
            continue;
        }

        ++tested;

        const ValidationResult v =
            validate_predicted_boundaries(
                p,
                m,
                s
            );

        const bool ok =
            v.first_previous_is_miss &&
            v.first_is_hit &&
            v.second_previous_is_miss &&
            v.second_is_hit;

        if (ok) {
            ++passed;
        } else {
            ++failed;

            if (failed <= 20) {
                print_failure(
                    "large boundary",
                    p,
                    m,
                    s,
                    v
                );
            }
        }
    }

    std::cout
        << "\nSUMMARY PHASE 2\n"
        << "attempted="
        << CASES
        << "\n"
        << "usable_cases="
        << tested
        << "\n"
        << "boundary_pass="
        << passed
        << "/" << tested
        << "\n"
        << "boundary_fail="
        << failed
        << "\n";
}

/*
    Specifically attack high-exponent cases.

    These have very large p^e gaps but are still cheap to test
    because we evaluate only four Lucas predicates.
*/
void run_high_exponent_phase() {
    std::cout
        << "\nPHASE 3: HIGH-EXPONENT BOUNDARIES\n";

    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (u64 p : primes) {
        std::cout
            << "\np=" << p
            << "\n";

        u64 power = 1;

        for (int k = 1;
             k <= 30;
             ++k) {

            if (power >
                UINT64_MAX / p) {

                break;
            }

            power *= p;

            if (power <= 1) {
                continue;
            }

            const u64 m =
                power - 1;

            const Structure s =
                analyze_structure(
                    p,
                    m
                );

            if (s.second_start > m) {
                std::cout
                    << "k=" << k
                    << " m=" << m
                    << " insufficient_second_start\n";

                continue;
            }

            const ValidationResult v =
                validate_predicted_boundaries(
                    p,
                    m,
                    s
                );

            const bool ok =
                v.first_previous_is_miss &&
                v.first_is_hit &&
                v.second_previous_is_miss &&
                v.second_is_hit;

            std::cout
                << "k=" << k
                << " m=" << m
                << " e=" << s.exponent
                << " s0=" << s.first_start
                << " s1=" << s.second_start
                << " gap=" << s.modulus
                << " "
                << (ok ? "PASS" : "FAIL")
                << "\n";
        }
    }
}

void run_digit_family_phase() {
    std::cout
        << "\nPHASE 4: DIGIT FAMILIES\n";

    struct TestCase {
        u64 p;
        u64 m;
    };

    const std::vector<TestCase> cases = {
        {2, 4194303},
        {3, 3188645},
        {5, 3906249},
        {7, 1647085},
        {11, 3543121},

        {3, 41},
        {3, 59},
        {3, 125},
        {3, 377},
        {3, 1133},
        {3, 3401},

        {5, 59},
        {5, 94},
        {5, 194},
        {5, 219}
    };

    u64 passed = 0;

    for (const auto& c : cases) {
        const Structure s =
            analyze_structure(
                c.p,
                c.m
            );

        if (s.second_start > c.m) {
            std::cout
                << "p=" << c.p
                << " m=" << c.m
                << " insufficient_second_start\n";

            continue;
        }

        const ValidationResult v =
            validate_predicted_boundaries(
                c.p,
                c.m,
                s
            );

        const bool ok =
            v.first_previous_is_miss &&
            v.first_is_hit &&
            v.second_previous_is_miss &&
            v.second_is_hit;

        if (ok) {
            ++passed;
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
            << " e="
            << s.exponent
            << " s0="
            << s.first_start
            << " s1="
            << s.second_start
            << " "
            << (ok ? "PASS" : "FAIL")
            << "\n";
    }

    std::cout
        << "\nSUMMARY PHASE 4\n"
        << "passed="
        << passed
        << "/" << cases.size()
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 168;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "B\n";

    run_small_exhaustive();
    run_large_boundary_phase();
    run_high_exponent_phase();
    run_digit_family_phase();

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "B\n";

    return 0;
}