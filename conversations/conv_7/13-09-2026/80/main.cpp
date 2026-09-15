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

struct StructuralCase {
    u64 p;
    unsigned a0;
    unsigned b;
    unsigned z;
    unsigned e;

    u64 s0;
    u64 q;
};

struct Result {
    bool pass;
    u64 cases;
    u64 failures;
    u64 skipped;
};

bool add_overflow_u64(u64 a, u64 b)
{
    return b > std::numeric_limits<u64>::max() - a;
}

bool mul_overflow_u64(u64 a, u64 b)
{
    if (a == 0 || b == 0) {
        return false;
    }

    return a > std::numeric_limits<u64>::max() / b;
}

u64 pow_u64(u64 p, unsigned e)
{
    u64 result = 1;

    for (unsigned i = 0; i < e; ++i) {
        result *= p;
    }

    return result;
}

bool safe_pow_u64(u64 p, unsigned e, u64& out)
{
    u64 result = 1;

    for (unsigned i = 0; i < e; ++i) {
        if (mul_overflow_u64(result, p)) {
            return false;
        }

        result *= p;
    }

    out = result;
    return true;
}

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

std::vector<u64> digits_of(u64 n, u64 p)
{
    std::vector<u64> digits;

    if (n == 0) {
        digits.push_back(0);
        return digits;
    }

    while (n > 0) {
        digits.push_back(n % p);
        n /= p;
    }

    return digits;
}

/*
 * C_r(q) = q_r * product_{i>r}(q_i + 1)
 *
 * Constructed directly from digits.
 */
std::vector<u64> construct_type_counts(
    const std::vector<u64>& digits)
{
    if (digits.empty()) {
        return {};
    }

    const unsigned R =
        static_cast<unsigned>(digits.size() - 1);

    std::vector<u64> counts(R + 1, 0);

    u64 suffix_product = 1;

    for (int r = static_cast<int>(R); r >= 0; --r) {
        counts[r] =
            digits[r] * suffix_product;

        suffix_product *=
            digits[r] + 1;
    }

    return counts;
}

/*
 * S(q) = sum_r C_r(q) * (p^r - (q mod p^r))
 *
 * The implementation uses u128 internally.
 */
u128 S_direct_u128(u64 q, u64 p)
{
    if (q == 0) {
        return 0;
    }

    const std::vector<u64> digits =
        digits_of(q, p);

    const std::vector<u64> counts =
        construct_type_counts(digits);

    u128 total = 0;
    u64 power = 1;

    for (unsigned r = 0; r < counts.size(); ++r) {
        if (r > 0) {
            power *= p;
        }

        const u64 remainder = q % power;
        const u64 A = power - remainder;

        total +=
            static_cast<u128>(counts[r]) *
            static_cast<u128>(A);
    }

    return total;
}

/*
 * I(q) = product_i(q_i+1)-1
 */
u128 I_direct_u128(u64 q, u64 p)
{
    const std::vector<u64> digits =
        digits_of(q, p);

    u128 product = 1;

    for (u64 digit : digits) {
        product *=
            static_cast<u128>(digit + 1);
    }

    return product - 1;
}

/*
 * H(q) = p^e * S(q) - s0 * I(q)
 */
u128 H_from_definition(
    u64 q,
    u64 p,
    u64 e_power,
    u64 s0)
{
    const u128 S =
        S_direct_u128(q, p);

    const u128 I =
        I_direct_u128(q, p);

    return
        static_cast<u128>(e_power) * S
        - static_cast<u128>(s0) * I;
}

/*
 * Equivalent closed formula:
 *
 * H(q) = q p^e - s0 I(q)
 */
u128 H_closed(
    u64 q,
    u64 e_power,
    u64 s0,
    u64 p)
{
    const u128 I =
        I_direct_u128(q, p);

    return
        static_cast<u128>(q) *
        static_cast<u128>(e_power)
        - static_cast<u128>(s0) * I;
}

/*
 * Structural parameterization:
 *
 * e = a0 + z + 1
 * s0 = (b+1)p^a0
 *
 * with
 *
 * p prime
 * a0,z >= 0
 * 0 <= b <= p-2
 * q >= 1
 *
 * m = s0 + q p^e - 1
 *
 * We only need s0, e and q for H.
 */
bool build_structural_case(
    u64 p,
    unsigned a0,
    unsigned b,
    unsigned z,
    u64 q,
    StructuralCase& out)
{
    if (p != 2 && p != 3 && p != 5) {
        return false;
    }

    if (b >= p - 1) {
        return false;
    }

    const unsigned e =
        a0 + z + 1;

    u64 p_a;

    if (!safe_pow_u64(p, a0, p_a)) {
        return false;
    }

    if (mul_overflow_u64(
            static_cast<u64>(b + 1),
            p_a)) {
        return false;
    }

    const u64 s0 =
        static_cast<u64>(b + 1) * p_a;

    u64 p_e;

    if (!safe_pow_u64(p, e, p_e)) {
        return false;
    }

    /*
     * Ensure q p^e itself is representable.
     */
    if (q != 0 &&
        mul_overflow_u64(q, p_e)) {
        return false;
    }

    out = {
        p,
        a0,
        b,
        z,
        e,
        s0,
        q
    };

    return true;
}

/*
 * Verify H recurrence for one structural case.
 */
bool test_case(
    const StructuralCase& sc,
    bool verbose)
{
    const u64 p = sc.p;

    u64 p_e;

    if (!safe_pow_u64(p, sc.e, p_e)) {
        return false;
    }

    const Decomposition d =
        decompose_q(sc.q, p);

    /*
     * ------------------------------------------------------------
     * 1. Verify S recurrence independently
     *
     * S(q) =
     *   (a+1)S(Q) + a(p^R-Q)
     * ------------------------------------------------------------
     */

    const u128 S_q =
        S_direct_u128(sc.q, p);

    const u128 S_Q =
        S_direct_u128(d.Q, p);

    u64 p_R;

    if (!safe_pow_u64(p, d.R, p_R)) {
        return false;
    }

    const u128 expected_S =
        static_cast<u128>(d.a + 1) * S_Q
        + static_cast<u128>(d.a) *
          (static_cast<u128>(p_R) - d.Q);

    bool pass = true;

    if (S_q != expected_S) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL S RECURRENCE\n"
                << "p=" << p
                << " q=" << sc.q
                << " Q=" << d.Q
                << " a=" << d.a
                << " R=" << d.R
                << " S(q)=" << static_cast<u64>(S_q)
                << " S(Q)=" << static_cast<u64>(S_Q)
                << " expected="
                << static_cast<u64>(expected_S)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 2. Verify S(q)=q
     * ------------------------------------------------------------
     */

    if (S_q != sc.q) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL S(q)=q\n"
                << "p=" << p
                << " q=" << sc.q
                << " S=" << static_cast<u64>(S_q)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 3. Verify I recurrence
     *
     * I(q) = (a+1)I(Q)+a
     * ------------------------------------------------------------
     */

    const u128 I_q =
        I_direct_u128(sc.q, p);

    const u128 I_Q =
        I_direct_u128(d.Q, p);

    const u128 expected_I =
        static_cast<u128>(d.a + 1) * I_Q
        + static_cast<u128>(d.a);

    if (I_q != expected_I) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL I RECURRENCE\n"
                << "p=" << p
                << " q=" << sc.q
                << " Q=" << d.Q
                << " a=" << d.a
                << " I(q)=" << static_cast<u64>(I_q)
                << " I(Q)=" << static_cast<u64>(I_Q)
                << " expected="
                << static_cast<u64>(expected_I)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 4. Verify H definition equals closed form
     *
     * H = p^e S - s0 I
     *
     * and
     *
     * H = q p^e - s0 I
     * ------------------------------------------------------------
     */

    const u128 H_q_definition =
        H_from_definition(
            sc.q,
            p,
            p_e,
            sc.s0);

    const u128 H_q_closed =
        H_closed(
            sc.q,
            p_e,
            sc.s0,
            p);

    if (H_q_definition != H_q_closed) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL H DEFINITION/CLOSED\n"
                << "p=" << p
                << " q=" << sc.q
                << " e=" << sc.e
                << " s0=" << sc.s0
                << " H_definition="
                << static_cast<u64>(H_q_definition)
                << " H_closed="
                << static_cast<u64>(H_q_closed)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 5. Verify the H recurrence directly
     *
     * H(q)
     *
     * =
     * (a+1)H(Q)
     * + a p^e (p^R-Q)
     * - a s0
     * ------------------------------------------------------------
     */

    const u128 H_Q =
        H_from_definition(
            d.Q,
            p,
            p_e,
            sc.s0);

    const u128 correction =
        static_cast<u128>(d.a) *
        static_cast<u128>(p_e) *
        (static_cast<u128>(p_R) - d.Q);

    const u128 expected_H =
        static_cast<u128>(d.a + 1) *
        H_Q
        + correction
        - static_cast<u128>(d.a) *
          static_cast<u128>(sc.s0);

    if (H_q_definition != expected_H) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL H RECURRENCE\n"
                << "p=" << p
                << " q=" << sc.q
                << " Q=" << d.Q
                << " a=" << d.a
                << " R=" << d.R
                << " e=" << sc.e
                << " s0=" << sc.s0
                << " H(q)="
                << static_cast<u64>(H_q_definition)
                << " H(Q)="
                << static_cast<u64>(H_Q)
                << " correction="
                << static_cast<u64>(correction)
                << " expected="
                << static_cast<u64>(expected_H)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 6. Verify the recurrence in expanded form.
     *
     * Starting from
     *
     * H(q)=p^e S(q)-s0 I(q)
     *
     * replace S(q) and I(q):
     *
     * p^e[(a+1)S(Q)+a(p^R-Q)]
     * -
     * s0[(a+1)I(Q)+a]
     *
     * ------------------------------------------------------------
     */

    const u128 expanded =
        static_cast<u128>(p_e) *
        (
            static_cast<u128>(d.a + 1) * S_Q
            + static_cast<u128>(d.a) *
              (static_cast<u128>(p_R) - d.Q)
        )
        -
        static_cast<u128>(sc.s0) *
        (
            static_cast<u128>(d.a + 1) * I_Q
            + static_cast<u128>(d.a)
        );

    if (H_q_definition != expanded) {
        pass = false;

        if (verbose) {
            std::cout
                << "FAIL H EXPANDED IDENTITY\n"
                << "p=" << p
                << " q=" << sc.q
                << " expanded="
                << static_cast<u64>(expanded)
                << " H="
                << static_cast<u64>(H_q_definition)
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * 7. Small structural output.
     * ------------------------------------------------------------
     */

    if (verbose) {
        std::cout
            << "STRUCTURE"
            << " p=" << p
            << " a0=" << sc.a0
            << " b=" << sc.b
            << " z=" << sc.z
            << " e=" << sc.e
            << " s0=" << sc.s0
            << " q=" << sc.q
            << " Q=" << d.Q
            << " a=" << d.a
            << " R=" << d.R
            << " S(q)=" << static_cast<u64>(S_q)
            << " I(q)=" << static_cast<u64>(I_q)
            << " H(q)="
            << static_cast<u64>(H_q_definition)
            << " pass=" << (pass ? 1 : 0)
            << '\n';
    }

    return pass;
}

bool deterministic_tests()
{
    struct Parameters {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    };

    const Parameters tests[] = {
        {2, 0, 0, 0, 1},
        {2, 0, 0, 0, 3},
        {2, 1, 0, 0, 7},
        {2, 2, 1, 1, 63},
        {3, 0, 0, 0, 1},
        {3, 1, 1, 0, 8},
        {3, 2, 1, 1, 80},
        {5, 0, 0, 0, 1},
        {5, 1, 2, 2, 24},
        {5, 3, 1, 4, 100},
        {2, 4, 0, 3, 1048575},
        {3, 4, 1, 5, 987654321ULL},
        {5, 4, 2, 6, 1000007654321ULL},
        {2, 5, 0, 2, 1073741825ULL}
    };

    bool all_pass = true;

    std::cout << "DETERMINISTIC TESTS\n";

    for (const auto& t : tests) {
        StructuralCase sc;

        if (!build_structural_case(
                t.p,
                t.a0,
                t.b,
                t.z,
                t.q,
                sc)) {

            std::cout
                << "SKIP DETERMINISTIC STRUCTURE\n"
                << "p=" << t.p
                << " q=" << t.q
                << '\n';

            continue;
        }

        const bool pass =
            test_case(sc, true);

        if (!pass) {
            all_pass = false;
        }
    }

    return all_pass;
}

Result random_pure_q_tests(
    u64 seed,
    u64 cases)
{
    std::mt19937_64 rng(seed);

    const u64 primes[] = {2, 3, 5};

    u64 failures = 0;

    for (u64 i = 0; i < cases; ++i) {
        const u64 p =
            primes[rng() % 3];

        const unsigned a0 =
            static_cast<unsigned>(rng() % 8);

        const unsigned z =
            static_cast<unsigned>(rng() % 8);

        const unsigned b =
            static_cast<unsigned>(
                rng() % (p - 1));

        const u64 q =
            1 + rng() % 5000000000000ULL;

        StructuralCase sc;

        if (!build_structural_case(
                p, a0, b, z, q, sc)) {
            continue;
        }

        const bool pass =
            test_case(sc, false);

        if (!pass) {
            ++failures;

            std::cout
                << "RANDOM FAILURE"
                << " case=" << i
                << " p=" << p
                << " a0=" << a0
                << " b=" << b
                << " z=" << z
                << " q=" << q
                << '\n';

            test_case(sc, true);
        }
    }

    return {
        failures == 0,
        cases,
        failures,
        0
    };
}

/*
 * Randomly generate q by choosing its highest digit a.
 *
 * This deliberately exercises all possible top-digit values:
 *
 *     1 <= a <= p-1
 *
 * and therefore tests the recurrence for non-binary leading digits.
 */
Result random_digit_extension_tests(
    u64 seed,
    u64 cases)
{
    std::mt19937_64 rng(seed);

    const u64 primes[] = {2, 3, 5};

    u64 failures = 0;

    for (u64 i = 0; i < cases; ++i) {
        const u64 p =
            primes[rng() % 3];

        const unsigned R =
            static_cast<unsigned>(
                1 + rng() % 20);

        const u64 pR =
            pow_u64(p, R);

        const u64 a =
            1 + rng() % (p - 1);

        const u64 Q =
            rng() % pR;

        const u64 q =
            Q + a * pR;

        const unsigned a0 =
            static_cast<unsigned>(rng() % 6);

        const unsigned z =
            static_cast<unsigned>(rng() % 6);

        const unsigned b =
            static_cast<unsigned>(
                rng() % (p - 1));

        StructuralCase sc;

        if (!build_structural_case(
                p, a0, b, z, q, sc)) {
            continue;
        }

        const bool pass =
            test_case(sc, false);

        if (!pass) {
            ++failures;

            std::cout
                << "EXTENSION FAILURE"
                << " case=" << i
                << " p=" << p
                << " Q=" << Q
                << " a=" << a
                << " R=" << R
                << " q=" << q
                << '\n';

            test_case(sc, true);
        }
    }

    return {
        failures == 0,
        cases,
        failures,
        0
    };
}

/*
 * Explicit Q=0 test suite.
 */
bool q_zero_tests()
{
    bool all_pass = true;

    std::cout << "Q=0 TESTS\n";

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (unsigned a = 1; a < p; ++a) {
            for (unsigned a0 = 0; a0 <= 3; ++a0) {
                for (unsigned z = 0; z <= 2; ++z) {
                    const unsigned b =
                        static_cast<unsigned>(
                            (a0 + z) % (p - 1));

                    StructuralCase sc;

                    if (!build_structural_case(
                            p,
                            a0,
                            b,
                            z,
                            a,
                            sc)) {
                        continue;
                    }

                    /*
                     * Since q=a<p, its leading digit is exactly a
                     * and Q=0.
                     */
                    const Decomposition d =
                        decompose_q(sc.q, sc.p);

                    if (d.Q != 0 ||
                        d.a != a) {
                        all_pass = false;

                        std::cout
                            << "FAIL Q=0 DECOMPOSITION\n"
                            << "p=" << p
                            << " q=" << sc.q
                            << " expected_a=" << a
                            << " actual_a=" << d.a
                            << " Q=" << d.Q
                            << '\n';

                        continue;
                    }

                    const bool pass =
                        test_case(sc, true);

                    if (!pass) {
                        all_pass = false;
                    }
                }
            }
        }
    }

    return all_pass;
}

/*
 * Recursive application of the H recurrence:
 *
 * q -> Q -> Q_2 -> ...
 *
 * This tests the recurrence at every digit-removal stage.
 */
bool recursive_chain(
    StructuralCase sc,
    bool verbose)
{
    u64 current = sc.q;

    u64 steps = 0;

    while (true) {
        StructuralCase current_case = sc;
        current_case.q = current;

        const bool pass =
            test_case(
                current_case,
                verbose);

        if (!pass) {
            return false;
        }

        ++steps;

        const Decomposition d =
            decompose_q(current, sc.p);

        if (d.Q == 0) {
            break;
        }

        current = d.Q;
    }

    return steps > 0;
}

int main()
{
    std::cout << "START EXPERIMENT 233\n";

    bool overall_pass = true;

    /*
     * ------------------------------------------------------------
     * Deterministic
     * ------------------------------------------------------------
     */

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
     * ------------------------------------------------------------
     * Q=0 boundary
     * ------------------------------------------------------------
     */

    const bool q_zero_pass =
        q_zero_tests();

    std::cout
        << "q_zero_pass="
        << (q_zero_pass ? 1 : 0)
        << '\n';

    if (!q_zero_pass) {
        overall_pass = false;
    }

    /*
     * ------------------------------------------------------------
     * Random structural parameters
     * ------------------------------------------------------------
     */

    const Result structural_random =
        random_pure_q_tests(
            233233233ULL,
            100000);

    std::cout
        << "structural_random_cases="
        << structural_random.cases
        << " structural_random_failures="
        << structural_random.failures
        << " structural_random_pass="
        << (structural_random.pass ? 1 : 0)
        << '\n';

    if (!structural_random.pass) {
        overall_pass = false;
    }

    /*
     * ------------------------------------------------------------
     * Random digit-extension construction
     * ------------------------------------------------------------
     */

    const Result extension_random =
        random_digit_extension_tests(
            733233233ULL,
            100000);

    std::cout
        << "extension_cases="
        << extension_random.cases
        << " extension_failures="
        << extension_random.failures
        << " extension_pass="
        << (extension_random.pass ? 1 : 0)
        << '\n';

    if (!extension_random.pass) {
        overall_pass = false;
    }

    /*
     * ------------------------------------------------------------
     * Recursive digit chains
     * ------------------------------------------------------------
     */

    std::cout << "RECURSIVE H CHAINS\n";

    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    } chains[] = {
        {2, 0, 0, 0, 1048575},
        {3, 1, 1, 2, 987654321ULL},
        {5, 2, 3, 4, 1000007654321ULL},
        {2, 4, 0, 2, 1073741825ULL},
        {3, 3, 1, 5, 4999999999999ULL},
        {5, 4, 2, 6, 4999999999999ULL}
    };

    u64 chain_cases = 0;
    u64 chain_failures = 0;

    for (const auto& c : chains) {
        StructuralCase sc;

        if (!build_structural_case(
                c.p,
                c.a0,
                c.b,
                c.z,
                c.q,
                sc)) {

            std::cout
                << "CHAIN SKIPPED"
                << " p=" << c.p
                << " q=" << c.q
                << '\n';

            continue;
        }

        ++chain_cases;

        const bool pass =
            recursive_chain(sc, false);

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << " e=" << sc.e
            << " s0=" << sc.s0
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            ++chain_failures;

            recursive_chain(sc, true);
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
     * ------------------------------------------------------------
     * Final
     * ------------------------------------------------------------
     */

    std::cout
        << "OVERALL_PASS="
        << (overall_pass ? 1 : 0)
        << '\n';

    std::cout << "FINISHED EXPERIMENT 233\n";

    return overall_pass ? 0 : 1;
}
