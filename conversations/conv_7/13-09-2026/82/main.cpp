#include <cstdint>
#include <iostream>
#include <random>
#include <limits>
#include <vector>
#include <algorithm>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct StructuralCase {
    u64 p;
    unsigned a0;
    unsigned b;
    unsigned z;
    unsigned e;
    u64 s0;
    u64 q;
};

struct Interval {
    u64 j;
    unsigned r;
    u64 x;
    u64 end;
    u64 length;
};

struct Result {
    bool pass;
    u64 cases;
    u64 failures;
};

bool mul_overflow_u64(u64 a, u64 b)
{
    if (a == 0 || b == 0) {
        return false;
    }

    return a > std::numeric_limits<u64>::max() / b;
}

bool add_overflow_u64(u64 a, u64 b)
{
    return b > std::numeric_limits<u64>::max() - a;
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

u64 pow_u64(u64 p, unsigned e)
{
    u64 result = 1;

    for (unsigned i = 0; i < e; ++i) {
        result *= p;
    }

    return result;
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

u64 digit_at(u64 n, u64 p, unsigned r)
{
    const u64 pr = pow_u64(p, r);
    return (n / pr) % p;
}

unsigned highest_digit_position(u64 n, u64 p)
{
    if (n == 0) {
        return 0;
    }

    unsigned r = 0;
    u64 power = 1;

    while (power <= n / p) {
        power *= p;
        ++r;
    }

    return r;
}

u64 product_digits_plus_one(u64 q, u64 p)
{
    const std::vector<u64> digits =
        digits_of(q, p);

    u64 product = 1;

    for (u64 digit : digits) {
        product *= digit + 1;
    }

    return product;
}

u64 interval_count(u64 q, u64 p)
{
    return product_digits_plus_one(q, p) - 1;
}

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
 * Independently generate every admissible digit tuple.
 *
 * If q has digits q_0,...,q_R, generate all tuples
 *
 *     0 <= j_i <= q_i.
 *
 * This does NOT use mixed-radix unranking.
 */
void generate_admissible_recursive(
    const std::vector<u64>& q_digits,
    std::size_t position,
    u64 p,
    u64 current_j,
    u64 current_power,
    std::vector<u64>& output)
{
    if (position == q_digits.size()) {
        output.push_back(current_j);
        return;
    }

    for (u64 digit = 0;
         digit <= q_digits[position];
         ++digit) {

        if (mul_overflow_u64(
                digit,
                current_power)) {
            continue;
        }

        const u64 contribution =
            digit * current_power;

        if (add_overflow_u64(
                current_j,
                contribution)) {
            continue;
        }

        if (mul_overflow_u64(
                current_power,
                p)) {
            /*
             * No further digit positions can be represented.
             */
            continue;
        }

        generate_admissible_recursive(
            q_digits,
            position + 1,
            p,
            current_j + contribution,
            current_power * p,
            output);
    }
}

/*
 * Wrapper around the independent direct tuple generation.
 */
std::vector<u64> generate_admissible_direct(
    u64 q,
    u64 p)
{
    const std::vector<u64> q_digits =
        digits_of(q, p);

    std::vector<u64> result;

    /*
     * The maximum generated value is q itself, so the number of
     * tuples is product(q_i+1), which is manageable for this
     * exhaustive experiment.
     */
    generate_admissible_recursive(
        q_digits,
        0,
        p,
        0,
        1,
        result);

    std::sort(
        result.begin(),
        result.end());

    return result;
}

/*
 * Mixed-radix unrank from the earlier construction.
 *
 * k = 0,...,I-1 maps to the admissible j values, excluding j=q.
 */
u64 unrank_j(
    u64 k,
    u64 q,
    u64 p)
{
    const std::vector<u64> q_digits =
        digits_of(q, p);

    u64 j = 0;
    u64 remaining = k;
    u64 place = 1;

    for (u64 digit : q_digits) {
        const u64 radix =
            digit + 1;

        const u64 ji =
            remaining % radix;

        if (mul_overflow_u64(
                ji,
                place)) {
            return 0;
        }

        j += ji * place;

        remaining /= radix;

        if (mul_overflow_u64(
                place,
                p)) {
            return 0;
        }

        place *= p;
    }

    return j;
}

/*
 * First digit position with j_r < q_r.
 */
unsigned interval_type(
    u64 j,
    u64 q,
    u64 p)
{
    const unsigned R =
        highest_digit_position(q, p);

    for (unsigned r = 0; r <= R; ++r) {
        const u64 jd =
            digit_at(j, p, r);

        const u64 qd =
            digit_at(q, p, r);

        if (jd < qd) {
            return r;
        }
    }

    return R + 1;
}

/*
 * Construct the actual HIT interval associated with j.
 */
bool build_interval(
    u64 j,
    unsigned r,
    const StructuralCase& sc,
    Interval& out)
{
    u64 p_e;

    if (!safe_pow_u64(
            sc.p,
            sc.e,
            p_e)) {
        return false;
    }

    u64 p_r;

    if (!safe_pow_u64(
            sc.p,
            r,
            p_r)) {
        return false;
    }

    const unsigned e_r =
        sc.e + r;

    u64 p_er;

    if (!safe_pow_u64(
            sc.p,
            e_r,
            p_er)) {
        return false;
    }

    if (mul_overflow_u64(
            j,
            p_e)) {
        return false;
    }

    const u64 j_step =
        j * p_e;

    if (add_overflow_u64(
            sc.s0,
            j_step)) {
        return false;
    }

    const u64 x =
        sc.s0 + j_step;

    const u64 quotient =
        j / p_r;

    if (quotient ==
        std::numeric_limits<u64>::max()) {
        return false;
    }

    const u64 multiplier =
        quotient + 1;

    if (mul_overflow_u64(
            multiplier,
            p_er)) {
        return false;
    }

    const u64 endpoint_plus_one =
        multiplier * p_er;

    if (endpoint_plus_one == 0) {
        return false;
    }

    const u64 end =
        endpoint_plus_one - 1;

    if (end < x) {
        return false;
    }

    const u64 length =
        end - x + 1;

    out = {
        j,
        r,
        x,
        end,
        length
    };

    return true;
}

/*
 * Verify that j is digitwise admissible.
 */
bool is_admissible(
    u64 j,
    u64 q,
    u64 p)
{
    const unsigned R =
        highest_digit_position(q, p);

    const unsigned JR =
        highest_digit_position(j, p);

    const unsigned max_r =
        (JR > R) ? JR : R;

    for (unsigned r = 0; r <= max_r; ++r) {
        if (digit_at(j, p, r) >
            digit_at(q, p, r)) {
            return false;
        }
    }

    return true;
}

/*
 * Main independent bijection test.
 */
bool test_bijection(
    const StructuralCase& sc,
    bool verbose)
{
    const u64 expected_interval_count =
        interval_count(sc.q, sc.p);

    /*
     * ------------------------------------------------------------
     * A. Directly generate every admissible tuple.
     * ------------------------------------------------------------
     */
    std::vector<u64> direct =
        generate_admissible_direct(
            sc.q,
            sc.p);

    /*
     * The direct tuple set includes the terminal tuple j=q.
     *
     * HIT intervals correspond to everything except j=q.
     */
    if (direct.empty()) {
        return false;
    }

    if (direct.back() != sc.q) {
        if (verbose) {
            std::cout
                << "FAIL TERMINAL DIRECT TUPLE"
                << " q=" << sc.q
                << " last=" << direct.back()
                << '\n';
        }

        return false;
    }

    const u64 direct_without_terminal =
        static_cast<u64>(direct.size() - 1);

    if (direct_without_terminal !=
        expected_interval_count) {

        if (verbose) {
            std::cout
                << "FAIL DIRECT COUNT"
                << " q=" << sc.q
                << " direct="
                << direct_without_terminal
                << " expected="
                << expected_interval_count
                << '\n';
        }

        return false;
    }

    /*
     * ------------------------------------------------------------
     * B. Verify every direct tuple is admissible and unique.
     * ------------------------------------------------------------
     */
    for (std::size_t i = 0;
         i < direct.size();
         ++i) {

        if (!is_admissible(
                direct[i],
                sc.q,
                sc.p)) {

            if (verbose) {
                std::cout
                    << "FAIL DIRECT ADMISSIBILITY"
                    << " index=" << i
                    << " j=" << direct[i]
                    << '\n';
            }

            return false;
        }

        if (i > 0 &&
            direct[i - 1] >= direct[i]) {

            if (verbose) {
                std::cout
                    << "FAIL DIRECT ORDER"
                    << " index=" << i
                    << " prev=" << direct[i - 1]
                    << " current=" << direct[i]
                    << '\n';
            }

            return false;
        }
    }

    /*
     * ------------------------------------------------------------
     * C. Compare direct enumeration with mixed-radix unranking.
     *
     * direct[k] must equal unrank_j(k).
     * ------------------------------------------------------------
     */
    for (u64 k = 0;
         k < expected_interval_count;
         ++k) {

        const u64 direct_j =
            direct[k];

        const u64 ranked_j =
            unrank_j(
                k,
                sc.q,
                sc.p);

        if (direct_j != ranked_j) {
            if (verbose) {
                std::cout
                    << "FAIL BIJECTION"
                    << " k=" << k
                    << " direct_j=" << direct_j
                    << " ranked_j=" << ranked_j
                    << '\n';
            }

            return false;
        }
    }

    /*
     * The final direct tuple is exactly q and must be excluded
     * from the HIT interval sequence.
     */
    const u64 terminal =
        direct.back();

    if (terminal != sc.q) {
        return false;
    }

    /*
     * ------------------------------------------------------------
     * D. Verify interval sequence independently from the direct j.
     * ------------------------------------------------------------
     */
    u64 previous_end = 0;
    bool have_previous = false;

    for (u64 k = 0;
         k < expected_interval_count;
         ++k) {

        const u64 j =
            direct[k];

        const unsigned r =
            interval_type(
                j,
                sc.q,
                sc.p);

        if (r > highest_digit_position(
                    sc.q,
                    sc.p)) {

            if (verbose) {
                std::cout
                    << "FAIL INTERVAL TYPE"
                    << " k=" << k
                    << " j=" << j
                    << '\n';
            }

            return false;
        }

        Interval interval;

        if (!build_interval(
                j,
                r,
                sc,
                interval)) {

            if (verbose) {
                std::cout
                    << "FAIL INTERVAL BUILD"
                    << " k=" << k
                    << " j=" << j
                    << " r=" << r
                    << '\n';
            }

            return false;
        }

        /*
         * First interval starts at s0.
         */
        if (!have_previous) {
            if (interval.x != sc.s0) {
                if (verbose) {
                    std::cout
                        << "FAIL FIRST INTERVAL START"
                        << " x=" << interval.x
                        << " expected="
                        << sc.s0
                        << '\n';
                }

                return false;
            }
        } else {
            /*
             * Consecutive HIT intervals are separated by exactly
             * s0 MISS points.
             */
            if (previous_end >
                std::numeric_limits<u64>::max() -
                    sc.s0 - 1) {
                return false;
            }

            const u64 expected_x =
                previous_end + sc.s0 + 1;

            if (interval.x != expected_x) {
                if (verbose) {
                    std::cout
                        << "FAIL INTERVAL ADJACENCY"
                        << " k=" << k
                        << " j=" << j
                        << " previous_end="
                        << previous_end
                        << " x=" << interval.x
                        << " expected="
                        << expected_x
                        << '\n';
                }

                return false;
            }
        }

        previous_end =
            interval.end;

        have_previous = true;

        /*
         * Small sample output.
         */
        if (verbose && k < 5) {
            std::cout
                << "BIJECTION"
                << " k=" << k
                << " j=" << j
                << " r=" << r
                << " x=" << interval.x
                << " end=" << interval.end
                << '\n';
        }
    }

    /*
     * ------------------------------------------------------------
     * E. Verify final interval ends at m-s0.
     * ------------------------------------------------------------
     */
    if (expected_interval_count > 0) {
        u64 p_e;

        if (!safe_pow_u64(
                sc.p,
                sc.e,
                p_e)) {
            return false;
        }

        const u128 expected_final_end =
            static_cast<u128>(sc.q) *
            static_cast<u128>(p_e) - 1;

        if (static_cast<u128>(previous_end) !=
            expected_final_end) {

            if (verbose) {
                std::cout
                    << "FAIL FINAL INTERVAL"
                    << " actual="
                    << previous_end
                    << " expected="
                    << static_cast<u64>(
                           expected_final_end)
                    << '\n';
            }

            return false;
        }
    }

    /*
     * ------------------------------------------------------------
     * F. Verify terminal exclusion:
     *
     * j=q is admissible digitwise but is exactly the unique tuple
     * that has no first digit with j_r<q_r.
     * ------------------------------------------------------------
     */
    if (interval_type(
            terminal,
            sc.q,
            sc.p)
        <= highest_digit_position(
            sc.q,
            sc.p)) {

        if (verbose) {
            std::cout
                << "FAIL TERMINAL TYPE"
                << " q=" << sc.q
                << " type="
                << interval_type(
                       terminal,
                       sc.q,
                       sc.p)
                << '\n';
        }

        return false;
    }

    if (verbose) {
        std::cout
            << "BIJECTION SUMMARY"
            << " p=" << sc.p
            << " q=" << sc.q
            << " admissible_tuples="
            << direct.size()
            << " HIT_intervals="
            << expected_interval_count
            << " terminal=" << terminal
            << " pass=1"
            << '\n';
    }

    return true;
}

Result deterministic_tests()
{
    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    } tests[] = {
        {2, 0, 0, 0, 1},
        {2, 0, 0, 0, 3},
        {2, 1, 0, 0, 7},
        {2, 2, 1, 1, 15},
        {2, 2, 1, 1, 31},
        {2, 3, 0, 2, 63},

        {3, 0, 0, 0, 2},
        {3, 1, 1, 0, 8},
        {3, 2, 1, 1, 26},
        {3, 2, 1, 1, 80},
        {3, 3, 0, 2, 242},

        {5, 0, 0, 0, 4},
        {5, 1, 1, 0, 24},
        {5, 2, 2, 1, 100},
        {5, 2, 1, 2, 124}
    };

    u64 case_count = 0;
    u64 failure_count = 0;

    std::cout
        << "DETERMINISTIC BIJECTION TESTS\n";

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
                << "SKIP"
                << " p=" << t.p
                << " q=" << t.q
                << '\n';

            continue;
        }

        ++case_count;

        const bool pass =
            test_bijection(
                sc,
                true);

        if (!pass) {
            ++failure_count;
        }
    }

    const bool overall =
        failure_count == 0;

    std::cout
        << "deterministic_cases="
        << case_count
        << " deterministic_failures="
        << failure_count
        << " deterministic_pass="
        << (overall ? 1 : 0)
        << '\n';

    return {
        overall,
        case_count,
        failure_count
    };
}

Result random_tests(
    u64 seed,
    u64 cases)
{
    std::mt19937_64 rng(seed);

    const u64 primes[] = {
        2,
        3,
        5
    };

    u64 failure_count = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 p =
            primes[rng() % 3];

        const unsigned a0 =
            static_cast<unsigned>(
                rng() % 4);

        const unsigned z =
            static_cast<unsigned>(
                rng() % 4);

        const unsigned b =
            static_cast<unsigned>(
                rng() % (p - 1));

        /*
         * Small enough for independent exhaustive tuple generation.
         */
        const u64 q =
            1 + rng() % 5000;

        StructuralCase sc;

        if (!build_structural_case(
                p,
                a0,
                b,
                z,
                q,
                sc)) {
            continue;
        }

        const bool pass =
            test_bijection(
                sc,
                false);

        if (!pass) {
            ++failure_count;

            std::cout
                << "RANDOM BIJECTION FAILURE"
                << " case=" << i
                << " p=" << p
                << " a0=" << a0
                << " b=" << b
                << " z=" << z
                << " q=" << q
                << '\n';

            test_bijection(
                sc,
                true);
        }
    }

    return {
        failure_count == 0,
        cases,
        failure_count
    };
}

/*
 * Verify the bijection recursively through the Q chain:
 *
 * q -> Q -> Q_2 -> ...
 */
bool recursive_chain(
    const StructuralCase& base,
    bool verbose,
    u64& steps)
{
    u64 current =
        base.q;

    steps = 0;

    while (true) {
        StructuralCase sc =
            base;

        sc.q =
            current;

        const bool pass =
            test_bijection(
                sc,
                verbose);

        if (!pass) {
            return false;
        }

        ++steps;

        const unsigned R =
            highest_digit_position(
                current,
                base.p);

        if (current == 0) {
            break;
        }

        const u64 p_R =
            pow_u64(
                base.p,
                R);

        const u64 a =
            current / p_R;

        const u64 Q =
            current - a * p_R;

        if (Q == 0) {
            /*
             * The next stage would be q=0, for which there are no
             * HIT intervals.
             */
            break;
        }

        current = Q;
    }

    return true;
}

int main()
{
    std::cout
        << "START EXPERIMENT 235\n";

    bool overall_pass = true;

    /*
     * ------------------------------------------------------------
     * Deterministic independent tuple enumeration
     * ------------------------------------------------------------
     */
    const Result deterministic =
        deterministic_tests();

    if (!deterministic.pass) {
        overall_pass = false;
    }

    /*
     * ------------------------------------------------------------
     * Random independent tuple enumeration
     * ------------------------------------------------------------
     */
    const Result random =
        random_tests(
            235235235ULL,
            10000);

    std::cout
        << "random_cases="
        << random.cases
        << " random_failures="
        << random.failures
        << " random_pass="
        << (random.pass ? 1 : 0)
        << '\n';

    if (!random.pass) {
        overall_pass = false;
    }

    /*
     * ------------------------------------------------------------
     * Recursive chains
     * ------------------------------------------------------------
     */
    std::cout
        << "RECURSIVE BIJECTION CHAINS\n";

    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    } chains[] = {
        {2, 0, 0, 0, 63},
        {3, 1, 1, 2, 242},
        {5, 1, 2, 1, 624},
        {3, 2, 1, 1, 1000},
        {5, 2, 1, 1, 2000}
    };

    u64 chain_cases = 0;
    u64 chain_failures = 0;
    u64 chain_steps = 0;

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

        u64 steps = 0;

        const bool pass =
            recursive_chain(
                sc,
                false,
                steps);

        chain_steps += steps;

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << " steps=" << steps
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            ++chain_failures;

            recursive_chain(
                sc,
                true,
                steps);
        }
    }

    std::cout
        << "chain_cases="
        << chain_cases
        << " chain_steps="
        << chain_steps
        << " chain_failures="
        << chain_failures
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

    std::cout
        << "FINISHED EXPERIMENT 235\n";

    return overall_pass ? 0 : 1;
}
