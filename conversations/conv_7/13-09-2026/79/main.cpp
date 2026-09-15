#include <cstdint>
#include <iostream>
#include <random>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Digits {
    u64 value;
    u64 p;
    std::vector<u64> d;
};

struct Decomposition {
    u64 q;
    u64 Q;
    u64 a;
    unsigned R;
};

struct Result {
    bool pass;
    u64 checks;
    u64 failures;
};

std::vector<u64> to_digits(u64 n, u64 p)
{
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

u64 from_digits(const std::vector<u64>& digits, u64 p)
{
    u64 value = 0;
    u64 power = 1;

    for (u64 digit : digits) {
        value += digit * power;
        power *= p;
    }

    return value;
}

u64 power_u64(u64 p, unsigned e)
{
    u64 result = 1;

    for (unsigned i = 0; i < e; ++i) {
        result *= p;
    }

    return result;
}

Decomposition decompose(u64 q, u64 p)
{
    unsigned R = 0;
    u64 power = 1;

    while (power <= q / p) {
        power *= p;
        ++R;
    }

    const u64 a = q / power;
    const u64 Q = q - a * power;

    return {q, Q, a, R};
}

/*
 * Direct digit-level definition:
 *
 * C_r = q_r * product_{i>r}(q_i + 1)
 *
 * This function does NOT use the previous type_count implementation.
 */
std::vector<u64> construct_type_counts(const std::vector<u64>& digits)
{
    const unsigned R = static_cast<unsigned>(digits.size()) - 1;

    std::vector<u64> counts(R + 1, 0);

    u64 suffix_product = 1;

    for (int r = static_cast<int>(R); r >= 0; --r) {
        counts[r] = digits[r] * suffix_product;
        suffix_product *= (digits[r] + 1);
    }

    return counts;
}

bool same_vector(const std::vector<u64>& a,
                 const std::vector<u64>& b)
{
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }

    return true;
}

bool test_case(u64 p, u64 q, bool verbose)
{
    if (q == 0) {
        return true;
    }

    const Decomposition dec = decompose(q, p);

    /*
     * Q=0 is represented by the single zero digit.
     */
    const std::vector<u64> Q_digits =
        to_digits(dec.Q, p);

    const std::vector<u64> q_digits =
        to_digits(dec.q, p);

    /*
     * ------------------------------------------------------------
     * 1. Structural digit decomposition
     * ------------------------------------------------------------
     */

    bool pass = true;

    if (q_digits.size() != dec.R + 1) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL DIGIT LENGTH\n"
                << "p=" << p
                << " q=" << q
                << " R=" << dec.R
                << '\n';
        }
    }

    if (q_digits[dec.R] != dec.a) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL TOP DIGIT\n"
                << "p=" << p
                << " q=" << q
                << " a=" << dec.a
                << " actual=" << q_digits[dec.R]
                << '\n';
        }
    }

    /*
     * q = Q + a p^R
     */
    if (from_digits(q_digits, p) != q) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL DIGIT RECONSTRUCTION q\n"
                << "p=" << p
                << " q=" << q
                << '\n';
        }
    }

    if (from_digits(Q_digits, p) != dec.Q) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL DIGIT RECONSTRUCTION Q\n"
                << "p=" << p
                << " Q=" << dec.Q
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 2. Build C(Q) directly from Q's digits
     * ------------------------------------------------------------
     */

    std::vector<u64> C_Q;

    if (dec.Q == 0) {
        /*
         * The zero value has no genuine type positions.
         *
         * Represent it by an empty contribution vector for the
         * lower-digit comparison.
         */
        C_Q.clear();
    } else {
        C_Q = construct_type_counts(Q_digits);
    }

    /*
     * C(q) is always R+1 entries.
     */
    const std::vector<u64> C_q =
        construct_type_counts(q_digits);

    /*
     * ------------------------------------------------------------
     * 3. The core digit-level scaling law
     *
     * For r < R:
     *
     *     C_r(q) = (a+1) C_r(Q)
     *
     * because the higher digits are exactly the digits of Q,
     * multiplied by the newly appended factor (a+1).
     *
     * At r=R:
     *
     *     C_R(q) = a.
     * ------------------------------------------------------------
     */

    for (unsigned r = 0; r < dec.R; ++r) {
        const u64 Cq = C_q[r];

        const u64 CQ =
            (r < C_Q.size()) ? C_Q[r] : 0;

        const u128 expected =
            static_cast<u128>(dec.a + 1) * CQ;

        if (static_cast<u128>(Cq) != expected) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL LOWER COUNT SCALING\n"
                    << "p=" << p
                    << " q=" << q
                    << " Q=" << dec.Q
                    << " a=" << dec.a
                    << " R=" << dec.R
                    << " r=" << r
                    << " Cq=" << Cq
                    << " CQ=" << CQ
                    << " expected="
                    << static_cast<u64>(expected)
                    << '\n';
            }
        }
    }

    const u64 top_count = C_q[dec.R];

    if (top_count != dec.a) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL TOP COUNT\n"
                << "p=" << p
                << " q=" << q
                << " Q=" << dec.Q
                << " a=" << dec.a
                << " R=" << dec.R
                << " C_R=" << top_count
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 4. Verify the full vector is exactly the transformed vector
     *
     *     C(q) =
     *       ((a+1) C(Q)_0,
     *        ...,
     *        (a+1) C(Q)_{R-1},
     *        a)
     * ------------------------------------------------------------
     */

    for (unsigned r = 0; r < dec.R; ++r) {
        const u64 CQ =
            (r < C_Q.size()) ? C_Q[r] : 0;

        const u128 expected =
            static_cast<u128>(dec.a + 1) * CQ;

        if (static_cast<u128>(C_q[r]) != expected) {
            pass = false;
        }
    }

    if (C_q[dec.R] != dec.a) {
        pass = false;
    }

    /*
     * ------------------------------------------------------------
     * 5. Verify the sum of the type counts
     *
     *     sum_r C_r(q)
     *       = product_i(q_i+1) - 1
     * ------------------------------------------------------------
     */

    u128 sum_Cq = 0;

    for (u64 c : C_q) {
        sum_Cq += c;
    }

    u128 digit_product_q = 1;

    for (u64 digit : q_digits) {
        digit_product_q *= digit + 1;
    }

    const u128 expected_sum =
        digit_product_q - 1;

    if (sum_Cq != expected_sum) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL COUNT SUM\n"
                << "p=" << p
                << " q=" << q
                << " sum_C=" << static_cast<u64>(sum_Cq)
                << " expected="
                << static_cast<u64>(expected_sum)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 6. Verify the transformation of the number of HIT intervals
     *
     * I(q) = product(q_i+1)-1
     *
     * I(Q) = product(Q_i+1)-1
     *
     * Therefore:
     *
     * I(q) = (a+1)I(Q) + a
     * ------------------------------------------------------------
     */

    u128 Iq = sum_Cq;

    u128 IQ = 0;

    if (!C_Q.empty()) {
        for (u64 c : C_Q) {
            IQ += c;
        }
    }

    const u128 expected_I =
        static_cast<u128>(dec.a + 1) * IQ + dec.a;

    if (Iq != expected_I) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL INTERVAL COUNT RECURRENCE\n"
                << "p=" << p
                << " q=" << q
                << " Q=" << dec.Q
                << " a=" << dec.a
                << " Iq=" << static_cast<u64>(Iq)
                << " IQ=" << static_cast<u64>(IQ)
                << " expected="
                << static_cast<u64>(expected_I)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 7. Print one small structural example when requested
     * ------------------------------------------------------------
     */

    if (verbose) {
        std::cout
            << "STRUCTURE"
            << " p=" << p
            << " q=" << q
            << " Q=" << dec.Q
            << " a=" << dec.a
            << " R=" << dec.R
            << " Cq=[";

        for (std::size_t i = 0; i < C_q.size(); ++i) {
            if (i != 0) {
                std::cout << ',';
            }

            std::cout << C_q[i];
        }

        std::cout << "] pass=" << (pass ? 1 : 0)
                  << '\n';
    }

    return pass;
}

bool deterministic_tests()
{
    struct Test {
        u64 p;
        u64 q;
    };

    const Test tests[] = {
        {2, 1},
        {2, 2},
        {2, 3},
        {2, 7},
        {2, 15},
        {2, 31},
        {2, 63},

        {3, 1},
        {3, 2},
        {3, 8},
        {3, 26},
        {3, 80},

        {5, 1},
        {5, 4},
        {5, 24},
        {5, 25},
        {5, 100},

        {2, 1048575},
        {3, 987654321},
        {5, 1000007654321ULL},
        {2, 1073741825ULL}
    };

    bool all_pass = true;

    std::cout << "DETERMINISTIC TESTS\n";

    for (const auto& t : tests) {
        const bool pass =
            test_case(t.p, t.q, true);

        if (!pass) {
            all_pass = false;
        }
    }

    return all_pass;
}

Result random_tests(u64 seed, u64 cases)
{
    std::mt19937_64 rng(seed);

    const u64 primes[] = {2, 3, 5};

    u64 failures = 0;

    for (u64 i = 0; i < cases; ++i) {
        const u64 p =
            primes[rng() % 3];

        const u64 q =
            1 + rng() % 5000000000000ULL;

        const bool pass =
            test_case(p, q, false);

        if (!pass) {
            ++failures;

            std::cout
                << "RANDOM FAILURE"
                << " case=" << i
                << " p=" << p
                << " q=" << q
                << '\n';

            test_case(p, q, true);
        }
    }

    return {
        failures == 0,
        cases,
        failures
    };
}

bool recursive_digit_chain(u64 p, u64 q, bool verbose)
{
    u64 current = q;

    while (current > 0) {
        const Decomposition d =
            decompose(current, p);

        const bool pass =
            test_case(p, current, false);

        if (!pass) {
            if (verbose) {
                std::cout
                    << "CHAIN FAILURE"
                    << " p=" << p
                    << " q=" << current
                    << '\n';

                test_case(p, current, true);
            }

            return false;
        }

        if (d.Q == 0) {
            break;
        }

        current = d.Q;
    }

    return true;
}

int main()
{
    std::cout << "START EXPERIMENT 232\n";

    bool overall_pass = true;

    const bool deterministic_pass =
        deterministic_tests();

    std::cout
        << "deterministic_pass="
        << (deterministic_pass ? 1 : 0)
        << '\n';

    if (!deterministic_pass) {
        overall_pass = false;
    }

    /*
     * Explicit Q=0 cases.
     */
    std::cout << "Q=0 TESTS\n";

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (u64 a = 1; a <= p; ++a) {
            const u64 q = a;

            const bool pass =
                test_case(p, q, true);

            if (!pass) {
                overall_pass = false;
            }
        }
    }

    const Result random_result =
        random_tests(232232232ULL, 100000);

    std::cout
        << "random_cases="
        << random_result.checks
        << " random_failures="
        << random_result.failures
        << " random_pass="
        << (random_result.pass ? 1 : 0)
        << '\n';

    if (!random_result.pass) {
        overall_pass = false;
    }

    /*
     * Follow Q recursively and verify the digit-extension law
     * at every level.
     */
    std::cout << "RECURSIVE DIGIT CHAINS\n";

    const struct {
        u64 p;
        u64 q;
    } chains[] = {
        {2, 1048575},
        {3, 987654321},
        {5, 1000007654321ULL},
        {2, 1073741825ULL},
        {3, 4999999999999ULL},
        {5, 4999999999999ULL}
    };

    u64 chain_failures = 0;
    u64 chain_cases = 0;

    for (const auto& c : chains) {
        ++chain_cases;

        const bool pass =
            recursive_digit_chain(c.p, c.q, true);

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            ++chain_failures;
        }
    }

    std::cout
        << "chain_cases=" << chain_cases
        << " chain_failures=" << chain_failures
        << " chain_pass="
        << (chain_failures == 0 ? 1 : 0)
        << '\n';

    if (chain_failures != 0) {
        overall_pass = false;
    }

    /*
     * Final structural summary.
     */
    std::cout
        << "OVERALL_PASS="
        << (overall_pass ? 1 : 0)
        << '\n';

    std::cout << "FINISHED EXPERIMENT 232\n";

    return overall_pass ? 0 : 1;
}
