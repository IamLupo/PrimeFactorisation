#include <cstdint>
#include <iostream>
#include <random>
#include <limits>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct Decomposition {
    u64 q;
    u64 Q;
    u64 a;
    unsigned R;
};

struct TestResult {
    bool pass;
    u64 checks;
    u64 failures;
};

Decomposition decompose_q(u64 q, u64 p)
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

u64 pow_u64(u64 p, unsigned e)
{
    u64 result = 1;

    for (unsigned i = 0; i < e; ++i) {
        result *= p;
    }

    return result;
}

u64 digit_at(u64 n, u64 p, unsigned r)
{
    const u64 pr = pow_u64(p, r);
    return (n / pr) % p;
}

u64 digit_product(u64 n, u64 p)
{
    if (n == 0) {
        return 1;
    }

    u64 product = 1;
    u64 current = n;

    while (current > 0) {
        const u64 digit = current % p;
        product *= (digit + 1);
        current /= p;
    }

    return product;
}

u64 type_count(u64 q, u64 p, unsigned r)
{
    const u64 qr = digit_at(q, p, r);

    u64 higher_product = 1;
    u64 current = q;

    for (unsigned i = 0; i <= r; ++i) {
        current /= p;
    }

    while (current > 0) {
        const u64 digit = current % p;
        higher_product *= (digit + 1);
        current /= p;
    }

    return qr * higher_product;
}

u64 A_term(u64 q, u64 p, unsigned r)
{
    const u64 pr = pow_u64(p, r);
    return pr - (q % pr);
}

u64 S_direct(u64 q, u64 p)
{
    if (q == 0) {
        return 0;
    }

    u64 result = 0;
    const Decomposition d = decompose_q(q, p);

    for (unsigned r = 0; r <= d.R; ++r) {
        const u64 C = type_count(q, p, r);
        const u64 A = A_term(q, p, r);

        result += C * A;
    }

    return result;
}

bool test_one(u64 q, u64 p, bool verbose)
{
    if (q == 0) {
        return true;
    }

    const Decomposition d = decompose_q(q, p);

    bool pass = true;

    /*
     * ------------------------------------------------------------
     * 1. Term-by-term decomposition
     * ------------------------------------------------------------
     */

    for (unsigned r = 0; r < d.R; ++r) {
        const u64 Cq = type_count(d.q, p, r);
        const u64 CQ = type_count(d.Q, p, r);

        const u64 Aq = A_term(d.q, p, r);
        const u64 AQ = A_term(d.Q, p, r);

        const u128 expected_C =
            static_cast<u128>(d.a + 1) * CQ;

        const u128 expected_T =
            static_cast<u128>(d.a + 1) * CQ * AQ;

        const u128 actual_T =
            static_cast<u128>(Cq) * Aq;

        if (static_cast<u128>(Cq) != expected_C) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL LOWER COUNT\n"
                    << "p=" << p
                    << " q=" << d.q
                    << " Q=" << d.Q
                    << " a=" << d.a
                    << " R=" << d.R
                    << " r=" << r
                    << " Cq=" << Cq
                    << " CQ=" << CQ
                    << " expected=" << static_cast<u64>(expected_C)
                    << '\n';
            }
        }

        if (Aq != AQ) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL LOWER A\n"
                    << "p=" << p
                    << " q=" << d.q
                    << " Q=" << d.Q
                    << " a=" << d.a
                    << " R=" << d.R
                    << " r=" << r
                    << " Aq=" << Aq
                    << " AQ=" << AQ
                    << '\n';
            }
        }

        if (actual_T != expected_T) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL LOWER TERM\n"
                    << "p=" << p
                    << " q=" << d.q
                    << " Q=" << d.Q
                    << " a=" << d.a
                    << " R=" << d.R
                    << " r=" << r
                    << " actual=" << static_cast<u64>(actual_T)
                    << " expected=" << static_cast<u64>(expected_T)
                    << '\n';
            }
        }
    }

    /*
     * ------------------------------------------------------------
     * 2. Top term
     * ------------------------------------------------------------
     */

    {
        const u64 CR = type_count(d.q, p, d.R);
        const u64 AR = A_term(d.q, p, d.R);

        const u64 pR = pow_u64(p, d.R);
        const u64 expected_A = pR - d.Q;

        if (CR != d.a) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL TOP COUNT\n"
                    << "p=" << p
                    << " q=" << d.q
                    << " Q=" << d.Q
                    << " a=" << d.a
                    << " R=" << d.R
                    << " CR=" << CR
                    << " expected=" << d.a
                    << '\n';
            }
        }

        if (AR != expected_A) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL TOP A\n"
                    << "p=" << p
                    << " q=" << d.q
                    << " Q=" << d.Q
                    << " a=" << d.a
                    << " R=" << d.R
                    << " AR=" << AR
                    << " expected=" << expected_A
                    << '\n';
            }
        }

        const u128 actual_T =
            static_cast<u128>(CR) * AR;

        const u128 expected_T =
            static_cast<u128>(d.a) * (pR - d.Q);

        if (actual_T != expected_T) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL TOP TERM\n"
                    << "p=" << p
                    << " q=" << d.q
                    << " Q=" << d.Q
                    << " a=" << d.a
                    << " R=" << d.R
                    << " actual=" << static_cast<u64>(actual_T)
                    << " expected=" << static_cast<u64>(expected_T)
                    << '\n';
            }
        }
    }

    /*
     * ------------------------------------------------------------
     * 3. Sum the decomposed terms
     * ------------------------------------------------------------
     */

    const u64 Sq = S_direct(d.q, p);
    const u64 SQ = S_direct(d.Q, p);

    const u64 pR = pow_u64(p, d.R);

    const u128 recurrence =
        static_cast<u128>(d.a + 1) * SQ
        + static_cast<u128>(d.a) * (pR - d.Q);

    if (static_cast<u128>(Sq) != recurrence) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL S RECURRENCE\n"
                << "p=" << p
                << " q=" << d.q
                << " Q=" << d.Q
                << " a=" << d.a
                << " R=" << d.R
                << " S(q)=" << Sq
                << " S(Q)=" << SQ
                << " expected=" << static_cast<u64>(recurrence)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 4. Verify S(q)=q independently
     * ------------------------------------------------------------
     */

    if (Sq != d.q) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL S(q)=q\n"
                << "p=" << p
                << " q=" << d.q
                << " S=" << Sq
                << '\n';
        }
    }

    if (SQ != d.Q) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL S(Q)=Q\n"
                << "p=" << p
                << " Q=" << d.Q
                << " S=" << SQ
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 5. Digit-product recurrence
     *
     * I(q)+1 = prod_i(q_i+1)
     *
     * q = Q + a p^R
     * => prod_i(q_i+1) = (a+1) prod_i(Q_i+1)
     * => I(q) = (a+1)I(Q)+a
     * ------------------------------------------------------------
     */

    const u64 Pq = digit_product(d.q, p);
    const u64 PQ = digit_product(d.Q, p);

    const u128 expected_P =
        static_cast<u128>(d.a + 1) * PQ;

    if (static_cast<u128>(Pq) != expected_P) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL DIGIT PRODUCT\n"
                << "p=" << p
                << " q=" << d.q
                << " Q=" << d.Q
                << " a=" << d.a
                << " product(q)=" << Pq
                << " product(Q)=" << PQ
                << " expected=" << static_cast<u64>(expected_P)
                << '\n';
        }
    }

    return pass;
}

bool run_deterministic()
{
    struct Case {
        u64 p;
        u64 q;
    };

    const Case cases[] = {
        {2, 1},
        {2, 2},
        {2, 3},
        {2, 7},
        {2, 63},
        {2, 255},
        {2, 1023},

        {3, 1},
        {3, 2},
        {3, 8},
        {3, 26},
        {3, 80},
        {3, 242},

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

    for (const auto& c : cases) {
        const bool pass = test_one(c.q, c.p, true);

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            all_pass = false;
        }
    }

    /*
     * Explicit Q=0 cases.
     */
    std::cout << "BASE CASE TESTS\n";

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (u64 a = 1; a <= 4; ++a) {
            const u64 q = a;

            const bool pass = test_one(q, p, true);

            std::cout
                << "p=" << p
                << " q=" << q
                << " Q=0"
                << " pass=" << (pass ? 1 : 0)
                << '\n';

            if (!pass) {
                all_pass = false;
            }
        }
    }

    return all_pass;
}

TestResult run_random(u64 seed, u64 cases)
{
    std::mt19937_64 rng(seed);

    const u64 primes[] = {2, 3, 5};

    u64 checks = 0;
    u64 failures = 0;

    for (u64 i = 0; i < cases; ++i) {
        const u64 p = primes[rng() % 3];

        /*
         * Keep q comfortably inside uint64_t while still spanning
         * multiple base-p digits.
         */
        const u64 q = 1 + (rng() % 5000000000000ULL);

        const bool pass = test_one(q, p, false);

        ++checks;

        if (!pass) {
            ++failures;

            std::cout
                << "RANDOM FAILURE"
                << " case=" << i
                << " p=" << p
                << " q=" << q
                << '\n';

            /*
             * Re-run verbosely so the first structural discrepancy
             * is visible.
             */
            test_one(q, p, true);
        }
    }

    return {failures == 0, checks, failures};
}

int main()
{
    std::cout << "START EXPERIMENT 231\n";

    bool overall_pass = true;

    const bool deterministic_pass = run_deterministic();

    std::cout
        << "deterministic_pass="
        << (deterministic_pass ? 1 : 0)
        << '\n';

    if (!deterministic_pass) {
        overall_pass = false;
    }

    const TestResult random_result =
        run_random(231231231ULL, 100000);

    std::cout
        << "random_cases=" << random_result.checks
        << " random_failures=" << random_result.failures
        << " random_pass=" << (random_result.pass ? 1 : 0)
        << '\n';

    if (!random_result.pass) {
        overall_pass = false;
    }

    /*
     * Direct confirmation of the recurrence itself on a few
     * deliberately chosen large values.
     */
    std::cout << "LARGE RECURSION CHECKS\n";

    const struct {
        u64 p;
        u64 q;
    } large_cases[] = {
        {2, 1000000000000ULL},
        {3, 1000000000000ULL},
        {5, 1000000000000ULL},
        {2, 4999999999999ULL},
        {3, 4999999999999ULL},
        {5, 4999999999999ULL}
    };

    for (const auto& c : large_cases) {
        const Decomposition d = decompose_q(c.q, c.p);

        const u64 Sq = S_direct(d.q, c.p);
        const u64 SQ = S_direct(d.Q, c.p);
        const u64 pR = pow_u64(c.p, d.R);

        const u128 rhs =
            static_cast<u128>(d.a + 1) * SQ
            + static_cast<u128>(d.a) * (pR - d.Q);

        const bool pass =
            static_cast<u128>(Sq) == rhs &&
            Sq == d.q;

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << " Q=" << d.Q
            << " a=" << d.a
            << " R=" << d.R
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            overall_pass = false;
        }
    }

    std::cout
        << "OVERALL_PASS="
        << (overall_pass ? 1 : 0)
        << '\n';

    std::cout << "FINISHED EXPERIMENT 231\n";

    return overall_pass ? 0 : 1;
}
