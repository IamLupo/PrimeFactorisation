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
    u64 p_e;
    u64 m;
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

    if (mul_overflow_u64(q, p_e)) {
        return false;
    }

    const u64 q_times_p_e =
        q * p_e;

    /*
     * m = s0 + q p^e - 1
     */
    if (add_overflow_u64(
            s0,
            q_times_p_e)) {
        return false;
    }

    const u64 m_plus_one =
        s0 + q_times_p_e;

    if (m_plus_one == 0) {
        return false;
    }

    const u64 m =
        m_plus_one - 1;

    out = {
        p,
        a0,
        b,
        z,
        e,
        s0,
        q,
        p_e,
        m
    };

    return true;
}

/*
 * Lucas digit condition:
 *
 * n <=_p m iff every base-p digit of n is <= the corresponding
 * base-p digit of m.
 *
 * Therefore:
 *
 *     MISS iff n <=_p m
 *     HIT  iff n not<=_p m
 */
bool digitwise_leq(
    u64 n,
    u64 m,
    u64 p)
{
    const unsigned R =
        highest_digit_position(m, p);

    const unsigned NR =
        highest_digit_position(n, p);

    const unsigned max_r =
        (NR > R) ? NR : R;

    for (unsigned r = 0; r <= max_r; ++r) {
        const u64 nd =
            digit_at(n, p, r);

        const u64 md =
            digit_at(m, p, r);

        if (nd > md) {
            return false;
        }
    }

    return true;
}

bool lucas_hit(
    u64 n,
    u64 m,
    u64 p)
{
    return !digitwise_leq(
        n,
        m,
        p);
}

/*
 * Mixed-radix unrank.
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
 * Type of interval:
 *
 * first r such that j_r < q_r.
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
 * Build the interval:
 *
 * x_j = s0 + j p^e
 *
 * E_j =
 * ((floor(j/p^r)+1)p^(e+r))-1
 */
bool build_interval(
    u64 j,
    unsigned r,
    const StructuralCase& sc,
    Interval& out)
{
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
            sc.p_e)) {
        return false;
    }

    const u64 j_step =
        j * sc.p_e;

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

    out = {
        j,
        r,
        x,
        end,
        end - x + 1
    };

    return true;
}

/*
 * Direct construction of every HIT interval, using the already
 * established mixed-radix indexing.
 */
bool construct_intervals(
    const StructuralCase& sc,
    std::vector<Interval>& intervals)
{
    const u64 count =
        interval_count(
            sc.q,
            sc.p);

    intervals.clear();
    intervals.reserve(
        static_cast<std::size_t>(count));

    for (u64 k = 0;
         k < count;
         ++k) {

        const u64 j =
            unrank_j(
                k,
                sc.q,
                sc.p);

        const unsigned r =
            interval_type(
                j,
                sc.q,
                sc.p);

        const unsigned R =
            highest_digit_position(
                sc.q,
                sc.p);

        if (r > R) {
            return false;
        }

        Interval interval;

        if (!build_interval(
                j,
                r,
                sc,
                interval)) {
            return false;
        }

        intervals.push_back(
            interval);
    }

    return true;
}

/*
 * Check whether n belongs to one of the constructed geometric
 * HIT intervals.
 */
bool geometry_hit(
    u64 n,
    const std::vector<Interval>& intervals)
{
    for (const Interval& interval : intervals) {
        if (n < interval.x) {
            return false;
        }

        if (n <= interval.end) {
            return true;
        }
    }

    return false;
}

/*
 * Faster geometry predicate for the bounded exhaustive loop.
 *
 * Intervals are sorted and non-overlapping, so use a moving
 * interval index rather than scanning from zero for every n.
 */
bool compare_geometry_to_lucas(
    const StructuralCase& sc,
    const std::vector<Interval>& intervals,
    bool verbose)
{
    std::size_t interval_index = 0;

    u64 geometry_hits = 0;
    u64 lucas_hits = 0;

    for (u64 n = 0;
         n <= sc.m;
         ++n) {

        while (interval_index < intervals.size() &&
               n > intervals[interval_index].end) {

            ++interval_index;
        }

        const bool geometry =
            interval_index < intervals.size() &&
            n >= intervals[interval_index].x &&
            n <= intervals[interval_index].end;

        const bool lucas =
            lucas_hit(
                n,
                sc.m,
                sc.p);

        if (geometry) {
            ++geometry_hits;
        }

        if (lucas) {
            ++lucas_hits;
        }

        if (geometry != lucas) {
            if (verbose) {
                std::cout
                    << "MISMATCH"
                    << " n=" << n
                    << " geometry="
                    << (geometry ? 1 : 0)
                    << " lucas="
                    << (lucas ? 1 : 0)
                    << " p=" << sc.p
                    << " m=" << sc.m
                    << '\n';

                const unsigned nr =
                    highest_digit_position(
                        n,
                        sc.p);

                const unsigned mr =
                    highest_digit_position(
                        sc.m,
                        sc.p);

                const unsigned max_r =
                    (nr > mr) ? nr : mr;

                std::cout << "DIGITS n=[";
                for (unsigned r = 0;
                     r <= max_r;
                     ++r) {

                    if (r != 0) {
                        std::cout << ',';
                    }

                    std::cout
                        << digit_at(
                               n,
                               sc.p,
                               r);
                }

                std::cout << "] m=[";

                for (unsigned r = 0;
                     r <= max_r;
                     ++r) {

                    if (r != 0) {
                        std::cout << ',';
                    }

                    std::cout
                        << digit_at(
                               sc.m,
                               sc.p,
                               r);
                }

                std::cout << "]\n";
            }

            return false;
        }
    }

    /*
     * Verify the total counts independently.
     */
    u128 domain_length =
        static_cast<u128>(sc.m) + 1;

    if (geometry_hits != lucas_hits) {
        if (verbose) {
            std::cout
                << "FAIL HIT COUNT"
                << " geometry="
                << geometry_hits
                << " lucas="
                << lucas_hits
                << '\n';
        }

        return false;
    }

    if (geometry_hits >
        static_cast<u128>(
            std::numeric_limits<u64>::max())) {
        return false;
    }

    /*
     * The HIT total should also equal H.
     */
    u128 expected_H = 0;

    const unsigned R =
        highest_digit_position(
            sc.q,
            sc.p);

    for (unsigned r = 0;
         r <= R;
         ++r) {

        u64 p_r;
        u64 L;

        if (!safe_pow_u64(
                sc.p,
                r,
                p_r)) {
            return false;
        }

        if (!safe_pow_u64(
                sc.p,
                sc.e,
                p_r)) {
            return false;
        }

        /*
         * Recompute L_r without reusing an interval result.
         */
        u64 p_r_only;
        u64 p_e_only;

        if (!safe_pow_u64(
                sc.p,
                r,
                p_r_only)) {
            return false;
        }

        if (!safe_pow_u64(
                sc.p,
                sc.e,
                p_e_only)) {
            return false;
        }

        const u64 remainder =
            sc.q % p_r_only;

        const u64 gap =
            p_r_only - remainder;

        if (mul_overflow_u64(
                gap,
                p_e_only)) {
            return false;
        }

        const u64 scaled_gap =
            gap * p_e_only;

        if (scaled_gap < sc.s0) {
            return false;
        }

        L =
            scaled_gap - sc.s0;

        /*
         * Type count.
         */
        const u64 C =
            [&]() -> u64 {
                const u64 qr =
                    digit_at(
                        sc.q,
                        sc.p,
                        r);

                u64 higher =
                    1;

                u64 current =
                    sc.q;

                for (unsigned i = 0;
                     i <= r;
                     ++i) {
                    current /= sc.p;
                }

                while (current > 0) {
                    const u64 d =
                        current % sc.p;

                    higher *= d + 1;
                    current /= sc.p;
                }

                return qr * higher;
            }();

        expected_H +=
            static_cast<u128>(C) *
            static_cast<u128>(L);
    }

    /*
     * Silence unused warning for the domain variable while also
     * retaining the conceptual domain check.
     */
    if (domain_length == 0) {
        return false;
    }

    if (static_cast<u128>(geometry_hits) !=
        expected_H) {

        if (verbose) {
            std::cout
                << "FAIL H COMPARISON"
                << " direct="
                << geometry_hits
                << " grouped="
                << static_cast<u64>(
                       expected_H)
                << '\n';
        }

        return false;
    }

    if (verbose) {
        std::cout
            << "LUCAS GEOMETRY"
            << " p=" << sc.p
            << " q=" << sc.q
            << " m=" << sc.m
            << " intervals="
            << intervals.size()
            << " HIT="
            << geometry_hits
            << " MISS="
            << (static_cast<u64>(
                    domain_length -
                    geometry_hits))
            << " pass=1"
            << '\n';
    }

    return true;
}

bool test_case(
    const StructuralCase& sc,
    bool verbose)
{
    std::vector<Interval> intervals;

    if (!construct_intervals(
            sc,
            intervals)) {

        if (verbose) {
            std::cout
                << "FAIL INTERVAL CONSTRUCTION\n";
        }

        return false;
    }

    /*
     * Verify intervals are ordered and non-overlapping.
     */
    for (std::size_t i = 1;
         i < intervals.size();
         ++i) {

        if (intervals[i - 1].end >=
            intervals[i].x) {

            if (verbose) {
                std::cout
                    << "FAIL INTERVAL ORDER"
                    << " previous_end="
                    << intervals[i - 1].end
                    << " current_x="
                    << intervals[i].x
                    << '\n';
            }

            return false;
        }
    }

    /*
     * First interval must start at s0.
     */
    if (!intervals.empty() &&
        intervals.front().x != sc.s0) {

        if (verbose) {
            std::cout
                << "FAIL FIRST HIT START"
                << " x="
                << intervals.front().x
                << " expected="
                << sc.s0
                << '\n';
        }

        return false;
    }

    /*
     * Final interval must end at m-s0.
     */
    if (!intervals.empty()) {
        const u64 expected_end =
            sc.m - sc.s0;

        if (intervals.back().end !=
            expected_end) {

            if (verbose) {
                std::cout
                    << "FAIL LAST HIT END"
                    << " end="
                    << intervals.back().end
                    << " expected="
                    << expected_end
                    << '\n';
            }

            return false;
        }
    }

    return compare_geometry_to_lucas(
        sc,
        intervals,
        verbose);
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
        << "DETERMINISTIC LUCAS-GEOMETRY TESTS\n";

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
            test_case(
                sc,
                true);

        std::cout
            << "CASE"
            << " p=" << sc.p
            << " q=" << sc.q
            << " m=" << sc.m
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            ++failure_count;
        }
    }

    return {
        failure_count == 0,
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
         * Explicitly bounded because this experiment scans [0,m].
         */
        const u64 q =
            1 + rng() % 1000;

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

        /*
         * Keep the actual scanned domain bounded.
         */
        if (sc.m > 5000000ULL) {
            continue;
        }

        const bool pass =
            test_case(
                sc,
                false);

        if (!pass) {
            ++failure_count;

            std::cout
                << "RANDOM FAILURE"
                << " case=" << i
                << " p=" << p
                << " a0=" << a0
                << " b=" << b
                << " z=" << z
                << " q=" << q
                << " m=" << sc.m
                << '\n';

            test_case(
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

bool recursive_chain(
    const StructuralCase& base,
    bool verbose,
    u64& steps)
{
    u64 current =
        base.q;

    steps = 0;

    while (current > 0) {
        StructuralCase sc =
            base;

        sc.q =
            current;

        /*
         * Rebuild m for the new q.
         */
        StructuralCase rebuilt;

        if (!build_structural_case(
                base.p,
                base.a0,
                base.b,
                base.z,
                current,
                rebuilt)) {
            return false;
        }

        const bool pass =
            test_case(
                rebuilt,
                verbose);

        if (!pass) {
            return false;
        }

        ++steps;

        const unsigned R =
            highest_digit_position(
                current,
                base.p);

        const u64 p_R =
            pow_u64(
                base.p,
                R);

        const u64 a =
            current / p_R;

        const u64 Q =
            current - a * p_R;

        if (Q == 0) {
            break;
        }

        current = Q;
    }

    return true;
}

int main()
{
    std::cout
        << "START EXPERIMENT 236\n";

    bool overall_pass = true;

    /*
     * Deterministic full-domain verification.
     */
    const Result deterministic =
        deterministic_tests();

    std::cout
        << "deterministic_cases="
        << deterministic.cases
        << " deterministic_failures="
        << deterministic.failures
        << " deterministic_pass="
        << (deterministic.pass ? 1 : 0)
        << '\n';

    if (!deterministic.pass) {
        overall_pass = false;
    }

    /*
     * Random full-domain verification.
     */
    const Result random =
        random_tests(
            236236236ULL,
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
     * Recursive q -> Q chains.
     */
    std::cout
        << "RECURSIVE LUCAS-GEOMETRY CHAINS\n";

    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    } chains[] = {
        {2, 0, 0, 0, 63},
        {3, 1, 1, 0, 80},
        {5, 1, 1, 0, 124},
        {3, 2, 1, 1, 242},
        {5, 2, 1, 1, 624}
    };

    u64 chain_cases = 0;
    u64 chain_steps = 0;
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

        if (sc.m > 5000000ULL) {
            std::cout
                << "CHAIN SKIPPED DOMAIN"
                << " p=" << c.p
                << " q=" << c.q
                << " m=" << sc.m
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

    std::cout
        << "OVERALL_PASS="
        << (overall_pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 236\n";

    return overall_pass ? 0 : 1;
}
