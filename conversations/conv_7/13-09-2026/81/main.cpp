#include <cstdint>
#include <iostream>
#include <random>
#include <limits>
#include <vector>

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

u64 product_digits_plus_one(u64 q, u64 p)
{
    const std::vector<u64> digits = digits_of(q, p);

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

u64 digit_at(u64 q, u64 p, unsigned r)
{
    const u64 power = pow_u64(p, r);
    return (q / power) % p;
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

/*
 * C_r(q) = q_r * product_{i>r}(q_i + 1)
 */
u64 type_count(u64 q, u64 p, unsigned r)
{
    const unsigned R =
        highest_digit_position(q, p);

    if (r > R) {
        return 0;
    }

    const u64 qr =
        digit_at(q, p, r);

    u64 higher_product = 1;
    u64 current = q;

    for (unsigned i = 0; i <= r; ++i) {
        current /= p;
    }

    while (current > 0) {
        const u64 digit = current % p;
        higher_product *= digit + 1;
        current /= p;
    }

    return qr * higher_product;
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

    const unsigned e = a0 + z + 1;

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
 * Mixed-radix unrank:
 *
 * admissible j satisfy
 *
 *     0 <= j_i <= q_i
 *
 * and rank k is interpreted in mixed radix q_i+1.
 */
u64 unrank_j(u64 k, u64 q, u64 p)
{
    const std::vector<u64> q_digits =
        digits_of(q, p);

    u64 j = 0;
    u64 place = 1;
    u64 remaining = k;

    for (u64 digit : q_digits) {
        const u64 radix = digit + 1;
        const u64 ji = remaining % radix;

        j += ji * place;

        remaining /= radix;
        place *= p;
    }

    return j;
}

/*
 * r = first digit where j_r < q_r.
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
 * Build:
 *
 *     x_j = s0 + j p^e
 *
 *     E_j =
 *       (floor(j/p^r)+1)p^(e+r)-1
 */
bool build_interval(
    u64 j,
    unsigned r,
    const StructuralCase& sc,
    Interval& out)
{
    u64 p_e;

    if (!safe_pow_u64(sc.p, sc.e, p_e)) {
        return false;
    }

    u64 p_r;

    if (!safe_pow_u64(sc.p, r, p_r)) {
        return false;
    }

    const unsigned e_r = sc.e + r;

    u64 p_er;

    if (!safe_pow_u64(sc.p, e_r, p_er)) {
        return false;
    }

    if (mul_overflow_u64(j, p_e)) {
        return false;
    }

    const u64 step = j * p_e;

    if (add_overflow_u64(sc.s0, step)) {
        return false;
    }

    const u64 x =
        sc.s0 + step;

    const u64 quotient =
        j / p_r;

    if (quotient == std::numeric_limits<u64>::max()) {
        return false;
    }

    const u64 multiplier =
        quotient + 1;

    if (mul_overflow_u64(multiplier, p_er)) {
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
 * L_r =
 *
 *     (p^r - (q mod p^r)) p^e - s0
 */
bool expected_type_length(
    unsigned r,
    const StructuralCase& sc,
    u64& out)
{
    u64 p_r;
    u64 p_e;

    if (!safe_pow_u64(sc.p, r, p_r)) {
        return false;
    }

    if (!safe_pow_u64(sc.p, sc.e, p_e)) {
        return false;
    }

    const u64 remainder =
        sc.q % p_r;

    const u64 gap =
        p_r - remainder;

    if (mul_overflow_u64(gap, p_e)) {
        return false;
    }

    const u64 scaled_gap =
        gap * p_e;

    if (scaled_gap < sc.s0) {
        return false;
    }

    out = scaled_gap - sc.s0;

    return true;
}

/*
 * Exhaustively construct all HIT intervals for a bounded q.
 */
bool exhaustive_geometry_check(
    const StructuralCase& sc,
    bool verbose)
{
    const u64 interval_total =
        interval_count(sc.q, sc.p);

    const unsigned R =
        highest_digit_position(sc.q, sc.p);

    std::vector<u64> actual_count(R + 1, 0);
    std::vector<u128> actual_length_sum(R + 1, 0);

    bool pass = true;

    u64 previous_end = 0;
    bool have_previous = false;

    for (u64 k = 0; k < interval_total; ++k) {
        const u64 j =
            unrank_j(k, sc.q, sc.p);

        /*
         * Check j_i <= q_i digitwise.
         */
        bool admissible = true;

        const unsigned jR =
            highest_digit_position(j, sc.p);

        const unsigned max_digit =
            (jR > R) ? jR : R;

        for (unsigned r = 0; r <= max_digit; ++r) {
            if (digit_at(j, sc.p, r) >
                digit_at(sc.q, sc.p, r)) {

                admissible = false;
                break;
            }
        }

        if (!admissible) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL ADMISSIBILITY"
                    << " k=" << k
                    << " j=" << j
                    << '\n';
            }

            continue;
        }

        const unsigned r =
            interval_type(
                j,
                sc.q,
                sc.p);

        if (r > R) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL TERMINAL TYPE"
                    << " k=" << k
                    << " j=" << j
                    << '\n';
            }

            continue;
        }

        Interval interval;

        if (!build_interval(
                j,
                r,
                sc,
                interval)) {

            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL INTERVAL BUILD"
                    << " k=" << k
                    << " j=" << j
                    << " r=" << r
                    << '\n';
            }

            continue;
        }

        /*
         * The first HIT interval begins at s0.
         */
        if (!have_previous) {
            if (interval.x != sc.s0) {
                pass = false;

                if (verbose) {
                    std::cout
                        << "FAIL FIRST START"
                        << " x=" << interval.x
                        << " expected=" << sc.s0
                        << '\n';
                }
            }
        } else {
            /*
             * Every MISS block between HIT intervals has length s0.
             */
            if (previous_end >
                std::numeric_limits<u64>::max() -
                    sc.s0 - 1) {

                pass = false;
                continue;
            }

            const u64 expected_x =
                previous_end + sc.s0 + 1;

            if (interval.x != expected_x) {
                pass = false;

                if (verbose) {
                    std::cout
                        << "FAIL MISS GAP"
                        << " k=" << k
                        << " previous_end=" << previous_end
                        << " x=" << interval.x
                        << " expected=" << expected_x
                        << '\n';
                }
            }
        }

        ++actual_count[r];

        actual_length_sum[r] +=
            static_cast<u128>(interval.length);

        previous_end = interval.end;
        have_previous = true;

        /*
         * Final HIT interval ends at m-s0:
         *
         *     m-s0 = q p^e - 1.
         */
        if (k + 1 == interval_total) {
            u64 p_e;

            if (!safe_pow_u64(
                    sc.p,
                    sc.e,
                    p_e)) {

                pass = false;
                continue;
            }

            const u128 expected_final_end =
                static_cast<u128>(sc.q) *
                static_cast<u128>(p_e) - 1;

            if (static_cast<u128>(interval.end) !=
                expected_final_end) {

                pass = false;

                if (verbose) {
                    std::cout
                        << "FAIL FINAL END"
                        << " actual=" << interval.end
                        << " expected="
                        << static_cast<u64>(
                               expected_final_end)
                        << '\n';
                }
            }
        }

        /*
         * Only show the first few actual intervals.
         */
        if (verbose && k < 5) {
            std::cout
                << "INTERVAL"
                << " k=" << k
                << " j=" << interval.j
                << " r=" << interval.r
                << " x=" << interval.x
                << " end=" << interval.end
                << " length=" << interval.length
                << '\n';
        }
    }

    /*
     * Compare actual geometry against C_r L_r.
     */
    for (unsigned r = 0; r <= R; ++r) {
        u64 expected_length;

        if (!expected_type_length(
                r,
                sc,
                expected_length)) {

            pass = false;
            continue;
        }

        const u64 expected_count =
            type_count(sc.q, sc.p, r);

        const u128 expected_group_sum =
            static_cast<u128>(expected_count) *
            static_cast<u128>(expected_length);

        if (actual_count[r] != expected_count) {
            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL TYPE COUNT"
                    << " r=" << r
                    << " actual=" << actual_count[r]
                    << " expected=" << expected_count
                    << '\n';
            }
        }

        if (actual_length_sum[r] !=
            expected_group_sum) {

            pass = false;

            if (verbose) {
                std::cout
                    << "FAIL TYPE LENGTH SUM"
                    << " r=" << r
                    << " actual="
                    << static_cast<u64>(
                           actual_length_sum[r])
                    << " expected="
                    << static_cast<u64>(
                           expected_group_sum)
                    << " C=" << expected_count
                    << " L=" << expected_length
                    << '\n';
            }
        }
    }

    /*
     * Total HIT length.
     */
    u128 total_hit_length = 0;

    for (unsigned r = 0; r <= R; ++r) {
        total_hit_length +=
            actual_length_sum[r];
    }

    u64 p_e;

    if (!safe_pow_u64(
            sc.p,
            sc.e,
            p_e)) {

        return false;
    }

    const u128 expected_hit_length =
        static_cast<u128>(sc.q) *
        static_cast<u128>(p_e)
        -
        static_cast<u128>(sc.s0) *
        static_cast<u128>(interval_total);

    if (total_hit_length !=
        expected_hit_length) {

        pass = false;

        if (verbose) {
            std::cout
                << "FAIL TOTAL HIT"
                << " actual="
                << static_cast<u64>(
                       total_hit_length)
                << " expected="
                << static_cast<u64>(
                       expected_hit_length)
                << '\n';
        }
    }

    /*
     * Total MISS length.
     *
     * There are I+1 MISS blocks, each length s0.
     */
    const u128 miss_length =
        (static_cast<u128>(interval_total) + 1) *
        static_cast<u128>(sc.s0);

    /*
     * The entire domain is
     *
     *     [0,m]
     *
     * with
     *
     *     m+1 = s0 + q p^e.
     */
    const u128 domain_length =
        static_cast<u128>(sc.s0)
        +
        static_cast<u128>(sc.q) *
        static_cast<u128>(p_e);

    if (total_hit_length + miss_length !=
        domain_length) {

        pass = false;

        if (verbose) {
            std::cout
                << "FAIL COMPLETE PARTITION"
                << " HIT="
                << static_cast<u64>(
                       total_hit_length)
                << " MISS="
                << static_cast<u64>(
                       miss_length)
                << " DOMAIN="
                << static_cast<u64>(
                       domain_length)
                << '\n';
        }
    }

    if (verbose) {
        std::cout
            << "GEOMETRY"
            << " p=" << sc.p
            << " q=" << sc.q
            << " e=" << sc.e
            << " s0=" << sc.s0
            << " intervals=" << interval_total
            << " HIT="
            << static_cast<u64>(
                   total_hit_length)
            << " MISS="
            << static_cast<u64>(
                   miss_length)
            << " pass=" << (pass ? 1 : 0)
            << '\n';
    }

    return pass;
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

    std::cout << "DETERMINISTIC GEOMETRY TESTS\n";

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
            exhaustive_geometry_check(sc, true);

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

    const u64 primes[] = {2, 3, 5};

    u64 failure_count = 0;

    for (u64 i = 0; i < cases; ++i) {
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
         * q is deliberately bounded so that all HIT intervals
         * can be explicitly enumerated.
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
            exhaustive_geometry_check(
                sc,
                false);

        if (!pass) {
            ++failure_count;

            std::cout
                << "RANDOM GEOMETRY FAILURE"
                << " case=" << i
                << " p=" << p
                << " a0=" << a0
                << " b=" << b
                << " z=" << z
                << " q=" << q
                << '\n';

            exhaustive_geometry_check(
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
 * Large cases:
 *
 * Do NOT enumerate individual intervals.
 *
 * Instead verify the grouped geometric formula
 *
 *     H = sum_r C_r L_r
 *
 * against the closed form
 *
 *     H = q p^e - s0 I.
 */
bool large_aggregate_check(
    const StructuralCase& sc,
    bool verbose)
{
    u64 p_e;

    if (!safe_pow_u64(
            sc.p,
            sc.e,
            p_e)) {
        return false;
    }

    const u64 interval_total =
        interval_count(sc.q, sc.p);

    const unsigned R =
        highest_digit_position(sc.q, sc.p);

    u128 grouped_hit_length = 0;

    for (unsigned r = 0; r <= R; ++r) {
        u64 L;

        if (!expected_type_length(
                r,
                sc,
                L)) {
            return false;
        }

        const u64 C =
            type_count(sc.q, sc.p, r);

        grouped_hit_length +=
            static_cast<u128>(C) *
            static_cast<u128>(L);
    }

    const u128 closed_hit_length =
        static_cast<u128>(sc.q) *
        static_cast<u128>(p_e)
        -
        static_cast<u128>(sc.s0) *
        static_cast<u128>(interval_total);

    const bool pass =
        grouped_hit_length ==
        closed_hit_length;

    if (verbose) {
        std::cout
            << "LARGE AGGREGATE"
            << " p=" << sc.p
            << " q=" << sc.q
            << " e=" << sc.e
            << " s0=" << sc.s0
            << " intervals=" << interval_total
            << " grouped="
            << static_cast<u64>(
                   grouped_hit_length)
            << " closed="
            << static_cast<u64>(
                   closed_hit_length)
            << " pass=" << (pass ? 1 : 0)
            << '\n';
    }

    return pass;
}

Result large_tests()
{
    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    } tests[] = {
        {2, 4, 0, 3, 1048575},
        {3, 4, 1, 5, 987654321ULL},
        {5, 2, 1, 4, 1000007654321ULL},
        {2, 5, 0, 2, 1073741825ULL},
        {3, 3, 1, 5, 4999999999999ULL}
    };

    u64 failure_count = 0;
    u64 case_count = 0;

    std::cout << "LARGE AGGREGATE TESTS\n";

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
                << "SKIP LARGE"
                << " p=" << t.p
                << " q=" << t.q
                << '\n';

            continue;
        }

        ++case_count;

        if (!large_aggregate_check(
                sc,
                true)) {
            ++failure_count;
        }
    }

    return {
        failure_count == 0,
        case_count,
        failure_count
    };
}

int main()
{
    std::cout << "START EXPERIMENT 234\n";

    bool overall_pass = true;

    /*
     * Exhaustive small geometry.
     */
    const Result deterministic =
        deterministic_tests();

    if (!deterministic.pass) {
        overall_pass = false;
    }

    /*
     * Random exhaustive geometry.
     */
    const Result random =
        random_tests(
            234234234ULL,
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
     * Large aggregate geometry.
     */
    const Result large =
        large_tests();

    std::cout
        << "large_cases="
        << large.cases
        << " large_failures="
        << large.failures
        << " large_pass="
        << (large.pass ? 1 : 0)
        << '\n';

    if (!large.pass) {
        overall_pass = false;
    }

    std::cout
        << "OVERALL_PASS="
        << (overall_pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 234\n";

    return overall_pass ? 0 : 1;
}