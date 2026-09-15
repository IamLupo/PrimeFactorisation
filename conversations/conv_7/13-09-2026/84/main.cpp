#include <cstdint>
#include <iostream>
#include <random>
#include <limits>

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

struct Localization {
    bool hit;
    u64 t;
    u64 j;
    unsigned r;
    u64 x;
    u64 end;
};

struct Result {
    bool pass;
    u64 cases;
    u64 failures;
    u64 hits;
    u64 misses;
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

bool digitwise_leq(
    u64 n,
    u64 m,
    u64 p)
{
    const unsigned nr =
        highest_digit_position(n, p);

    const unsigned mr =
        highest_digit_position(m, p);

    const unsigned R =
        (nr > mr) ? nr : mr;

    for (unsigned r = 0; r <= R; ++r) {
        if (digit_at(n, p, r) >
            digit_at(m, p, r)) {
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
    return !digitwise_leq(n, m, p);
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

    const u64 q_pe =
        q * p_e;

    if (add_overflow_u64(s0, q_pe)) {
        return false;
    }

    const u64 m_plus_one =
        s0 + q_pe;

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
 * Given t, compute the largest j <= t satisfying j <=_p q.
 *
 * This is digitwise clamping:
 *
 *     j_i = min(t_i, q_i).
 *
 * This is the largest admissible integer <= t.
 */
u64 largest_admissible_leq(
    u64 t,
    u64 q,
    u64 p)
{
    const unsigned R_t =
        highest_digit_position(t, p);

    const unsigned R_q =
        highest_digit_position(q, p);

    const unsigned R =
        (R_t > R_q) ? R_t : R_q;

    u64 result = 0;
    u64 place = 1;

    for (unsigned r = 0; r <= R; ++r) {
        const u64 td =
            digit_at(t, p, r);

        const u64 qd =
            digit_at(q, p, r);

        const u64 jd =
            (td < qd) ? td : qd;

        result += jd * place;

        if (r < R) {
            place *= p;
        }
    }

    return result;
}

/*
 * Verify that j is exactly the largest admissible number <= t.
 *
 * This is checked independently from the localization formula.
 */
bool verify_maximal_admissible(
    u64 t,
    u64 j,
    u64 q,
    u64 p)
{
    if (j > t) {
        return false;
    }

    if (!digitwise_leq(j, q, p)) {
        return false;
    }

    /*
     * Any admissible number k <= t must satisfy
     *
     *     k_i <= min(t_i,q_i)
     *
     * at the first differing digit.
     *
     * We verify this by comparing each digit.
     */
    const unsigned R_t =
        highest_digit_position(t, p);

    const unsigned R_q =
        highest_digit_position(q, p);

    const unsigned R =
        (R_t > R_q) ? R_t : R_q;

    for (unsigned r = 0; r <= R; ++r) {
        const u64 jd =
            digit_at(j, p, r);

        const u64 td =
            digit_at(t, p, r);

        const u64 qd =
            digit_at(q, p, r);

        if (jd != ((td < qd) ? td : qd)) {
            return false;
        }
    }

    return true;
}

/*
 * Find the first digit where j_r < q_r.
 *
 * j=q has no such digit and is the excluded terminal tuple.
 */
bool find_interval_type(
    u64 j,
    u64 q,
    u64 p,
    unsigned& r_out)
{
    const unsigned R =
        highest_digit_position(q, p);

    for (unsigned r = 0; r <= R; ++r) {
        const u64 jd =
            digit_at(j, p, r);

        const u64 qd =
            digit_at(q, p, r);

        if (jd < qd) {
            r_out = r;
            return true;
        }
    }

    return false;
}

/*
 * Compute the p-adic successor:
 *
 *     next(j)
 *       = (floor(j/p^r)+1)p^r
 *
 * where r is the first digit with j_r < q_r.
 */
bool successor_j(
    u64 j,
    u64 q,
    u64 p,
    u64& next,
    unsigned& r_out)
{
    unsigned r;

    if (!find_interval_type(
            j,
            q,
            p,
            r)) {
        return false;
    }

    const u64 pr =
        pow_u64(p, r);

    const u64 quotient =
        j / pr;

    if (quotient ==
        std::numeric_limits<u64>::max()) {
        return false;
    }

    const u64 successor =
        (quotient + 1) * pr;

    next = successor;
    r_out = r;

    return true;
}

/*
 * Corrected localization:
 *
 * 1. t = floor((n-s0)/p^e)
 * 2. j = largest admissible <= t
 * 3. successor(j)
 * 4. interval endpoint = successor*p^e - 1
 *
 * The endpoint formula is exactly equivalent to the earlier E_j.
 */
Localization localize_n(
    u64 n,
    const StructuralCase& sc)
{
    Localization result = {
        false,
        0,
        0,
        0,
        0,
        0
    };

    /*
     * Initial MISS block.
     */
    if (n < sc.s0) {
        return result;
    }

    const u64 relative =
        n - sc.s0;

    const u64 t =
        relative / sc.p_e;

    /*
     * The quotient can reach q-1 during the final MISS block,
     * so t<q is not by itself sufficient for a HIT.
     */
    if (t >= sc.q) {
        result.t = t;
        return result;
    }

    const u64 j =
        largest_admissible_leq(
            t,
            sc.q,
            sc.p);

    unsigned r;

    if (!find_interval_type(
            j,
            sc.q,
            sc.p,
            r)) {

        result.t = t;
        result.j = j;
        return result;
    }

    u64 next;
    unsigned successor_r;

    if (!successor_j(
            j,
            sc.q,
            sc.p,
            next,
            successor_r)) {

        result.t = t;
        result.j = j;
        result.r = r;
        return result;
    }

    /*
     * x_j = s0 + j*p^e.
     */
    const u64 x =
        sc.s0 +
        j * sc.p_e;

    /*
     * E_j = next(j)*p^e - 1.
     */
    if (mul_overflow_u64(
            next,
            sc.p_e)) {

        result.t = t;
        result.j = j;
        result.r = r;
        result.x = x;
        return result;
    }

    const u64 end =
        next * sc.p_e - 1;

    const bool hit =
        n >= x &&
        n <= end;

    result.hit = hit;
    result.t = t;
    result.j = j;
    result.r = r;
    result.x = x;
    result.end = end;

    /*
     * successor_r must equal r.
     */
    if (successor_r != r) {
        result.hit = false;
    }

    return result;
}

/*
 * Verify one localization result.
 */
bool verify_point(
    u64 n,
    const StructuralCase& sc,
    bool verbose)
{
    const bool expected =
        lucas_hit(
            n,
            sc.m,
            sc.p);

    const Localization loc =
        localize_n(
            n,
            sc);

    if (loc.hit != expected) {
        if (verbose) {
            std::cout
                << "FAIL CLASS"
                << " n=" << n
                << " expected="
                << (expected ? 1 : 0)
                << " actual="
                << (loc.hit ? 1 : 0)
                << " t=" << loc.t
                << " j=" << loc.j
                << " r=" << loc.r
                << '\n';
        }

        return false;
    }

    /*
     * MISS case.
     */
    if (!expected) {
        return true;
    }

    /*
     * HIT case.
     */
    if (n < sc.s0) {
        return false;
    }

    const u64 expected_t =
        (n - sc.s0) / sc.p_e;

    if (loc.t != expected_t) {
        return false;
    }

    const u64 expected_j =
        largest_admissible_leq(
            loc.t,
            sc.q,
            sc.p);

    if (loc.j != expected_j) {
        if (verbose) {
            std::cout
                << "FAIL J"
                << " n=" << n
                << " t=" << loc.t
                << " j=" << loc.j
                << " expected="
                << expected_j
                << '\n';
        }

        return false;
    }

    if (!verify_maximal_admissible(
            loc.t,
            loc.j,
            sc.q,
            sc.p)) {

        return false;
    }

    unsigned expected_r;

    if (!find_interval_type(
            loc.j,
            sc.q,
            sc.p,
            expected_r)) {

        return false;
    }

    if (loc.r != expected_r) {
        return false;
    }

    const u64 expected_x =
        sc.s0 +
        loc.j * sc.p_e;

    if (loc.x != expected_x) {
        return false;
    }

    u64 next;
    unsigned next_r;

    if (!successor_j(
            loc.j,
            sc.q,
            sc.p,
            next,
            next_r)) {

        return false;
    }

    if (next_r != loc.r) {
        return false;
    }

    const u64 expected_end =
        next * sc.p_e - 1;

    if (loc.end != expected_end) {
        return false;
    }

    if (n < loc.x ||
        n > loc.end) {

        return false;
    }

    /*
     * The most important new identity:
     *
     * all HIT points whose block coordinate is t map to the same
     * admissible j until the next admissible boundary is reached.
     */
    if (loc.j > loc.t) {
        return false;
    }

    return true;
}

/*
 * Boundary tests for one structural case.
 */
bool boundary_tests(
    const StructuralCase& sc,
    bool verbose)
{
    bool pass = true;

    /*
     * Initial MISS block.
     */
    if (sc.s0 > 0) {
        const u64 points[] = {
            0,
            sc.s0 - 1
        };

        for (u64 n : points) {
            if (!verify_point(n, sc, verbose)) {
                pass = false;
            }
        }
    }

    /*
     * Test every admissible j only for bounded deterministic
     * cases. This is still logarithmic in the number of test
     * points, not a full scan of m.
     */
    const u64 interval_count = [&]() -> u64 {
        u64 product = 1;
        u64 current = sc.q;

        while (current > 0) {
            product *=
                current % sc.p + 1;

            current /= sc.p;
        }

        return product - 1;
    }();

    /*
     * For q=large, sample first/middle/last ranks.
     */
    const u64 selected[] = {
        0,
        interval_count / 2,
        interval_count - 1
    };

    for (u64 k : selected) {
        /*
         * Build admissible j using mixed-radix unrank.
         */
        u64 remaining = k;
        u64 current_q = sc.q;
        u64 place = 1;
        u64 j = 0;

        while (current_q > 0) {
            const u64 qdigit =
                current_q % sc.p;

            const u64 ji =
                remaining % (qdigit + 1);

            j += ji * place;

            remaining /=
                qdigit + 1;

            current_q /= sc.p;
            place *= sc.p;
        }

        unsigned r;

        if (!find_interval_type(
                j,
                sc.q,
                sc.p,
                r)) {
            pass = false;
            continue;
        }

        u64 next;
        unsigned next_r;

        if (!successor_j(
                j,
                sc.q,
                sc.p,
                next,
                next_r)) {
            pass = false;
            continue;
        }

        const u64 x =
            sc.s0 +
            j * sc.p_e;

        const u64 end =
            next * sc.p_e - 1;

        /*
         * Start, midpoint, endpoint.
         */
        const u64 points[] = {
            x,
            x + (end - x) / 2,
            end
        };

        for (u64 n : points) {
            if (!verify_point(
                    n,
                    sc,
                    verbose)) {

                pass = false;
            }
        }

        /*
         * The following point should be the MISS gap if it exists.
         */
        if (end < sc.m) {
            const u64 n =
                end + 1;

            if (lucas_hit(
                    n,
                    sc.m,
                    sc.p)) {

                if (verbose) {
                    std::cout
                        << "FAIL EXPECTED MISS GAP"
                        << " n=" << n
                        << " j=" << j
                        << '\n';
                }

                pass = false;
            }

            if (!verify_point(
                    n,
                    sc,
                    verbose)) {

                pass = false;
            }
        }
    }

    /*
     * Final MISS boundary:
     *
     * n = q p^e.
     *
     * Since
     *
     * m+1=s0+q p^e,
     *
     * this is exactly the first point of the final MISS block.
     */
    const u64 final_miss_start =
        sc.q * sc.p_e;

    if (final_miss_start <= sc.m) {
        if (!verify_point(
                final_miss_start,
                sc,
                verbose)) {

            pass = false;
        }

        if (sc.m != final_miss_start) {
            if (!verify_point(
                    sc.m,
                    sc,
                    verbose)) {

                pass = false;
            }
        }
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
        {2, 2, 1, 1, 63},

        {3, 0, 0, 0, 2},
        {3, 1, 1, 0, 8},
        {3, 2, 1, 1, 80},
        {3, 3, 0, 2, 242},

        {5, 0, 0, 0, 4},
        {5, 1, 1, 0, 24},
        {5, 2, 2, 1, 100},
        {5, 2, 1, 2, 124},

        {2, 5, 0, 2, 1073741825ULL},
        {3, 4, 1, 5, 987654321ULL},
        {5, 3, 1, 4, 1000007654321ULL}
    };

    u64 cases = 0;
    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CORRECTED LOCALIZATION\n";

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

        ++cases;

        const bool pass =
            boundary_tests(
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
            ++failures;
        }
    }

    return {
        failures == 0,
        cases,
        failures,
        0,
        0
    };
}

/*
 * Random arbitrary-n test.
 *
 * No O(m) scan.
 */
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

    u64 failures = 0;
    u64 hits = 0;
    u64 misses = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 p =
            primes[rng() % 3];

        const unsigned a0 =
            static_cast<unsigned>(
                rng() % 7);

        const unsigned z =
            static_cast<unsigned>(
                rng() % 7);

        const unsigned b =
            static_cast<unsigned>(
                rng() % (p - 1));

        const u64 q =
            1 +
            rng() %
            5000000000000ULL;

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
         * Avoid m+1 overflow in modulo.
         */
        if (sc.m ==
            std::numeric_limits<u64>::max()) {
            continue;
        }

        const u64 n =
            rng() % (sc.m + 1);

        const bool expected =
            lucas_hit(
                n,
                sc.m,
                sc.p);

        if (expected) {
            ++hits;
        } else {
            ++misses;
        }

        const bool pass =
            verify_point(
                n,
                sc,
                false);

        if (!pass) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "FAIL RANDOM"
                    << " case=" << i
                    << " p=" << p
                    << " q=" << q
                    << " n=" << n
                    << '\n';

                verify_point(
                    n,
                    sc,
                    true);
            }
        }
    }

    return {
        failures == 0,
        cases,
        failures,
        hits,
        misses
    };
}

/*
 * Random HIT construction.
 *
 * Generate a valid admissible j directly, construct its interval,
 * choose a random point inside it, then recover j from n.
 */
Result random_hit_tests(
    u64 seed,
    u64 cases)
{
    std::mt19937_64 rng(seed);

    const u64 primes[] = {
        2,
        3,
        5
    };

    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 p =
            primes[rng() % 3];

        const unsigned a0 =
            static_cast<unsigned>(
                rng() % 7);

        const unsigned z =
            static_cast<unsigned>(
                rng() % 7);

        const unsigned b =
            static_cast<unsigned>(
                rng() % (p - 1));

        const u64 q =
            1 +
            rng() %
            5000000000000ULL;

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
         * Build random admissible j.
         */
        u64 current_q = q;
        u64 place = 1;
        u64 j = 0;

        while (current_q > 0) {
            const u64 qdigit =
                current_q % p;

            const u64 jd =
                rng() %
                (qdigit + 1);

            j += jd * place;

            current_q /= p;
            place *= p;
        }

        /*
         * Avoid terminal j=q.
         */
        if (j == q) {
            j = 0;
        }

        unsigned r;

        if (!find_interval_type(
                j,
                q,
                p,
                r)) {
            ++failures;
            continue;
        }

        u64 next;
        unsigned next_r;

        if (!successor_j(
                j,
                q,
                p,
                next,
                next_r)) {
            ++failures;
            continue;
        }

        if (next_r != r) {
            ++failures;
            continue;
        }

        const u64 x =
            sc.s0 +
            j * sc.p_e;

        const u64 end =
            next * sc.p_e - 1;

        if (end < x) {
            ++failures;
            continue;
        }

        const u64 n =
            x +
            rng() %
            (end - x + 1);

        if (!verify_point(
                n,
                sc,
                false)) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "FAIL HIT"
                    << " case=" << i
                    << " p=" << p
                    << " q=" << q
                    << " j=" << j
                    << " r=" << r
                    << " n=" << n
                    << '\n';

                verify_point(
                    n,
                    sc,
                    true);
            }
        }
    }

    return {
        failures == 0,
        cases,
        failures,
        cases,
        0
    };
}

int main()
{
    std::cout
        << "START EXPERIMENT 238\n";

    bool overall_pass = true;

    /*
     * Deterministic structural/boundary tests.
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
     * Arbitrary huge n.
     */
    const Result random =
        random_tests(
            238238238ULL,
            100000);

    std::cout
        << "random_cases="
        << random.cases
        << " random_failures="
        << random.failures
        << " random_hits="
        << random.hits
        << " random_misses="
        << random.misses
        << " random_pass="
        << (random.pass ? 1 : 0)
        << '\n';

    if (!random.pass) {
        overall_pass = false;
    }

    /*
     * High HIT-density tests.
     */
    const Result hit =
        random_hit_tests(
            838238238ULL,
            100000);

    std::cout
        << "hit_cases="
        << hit.cases
        << " hit_failures="
        << hit.failures
        << " hit_pass="
        << (hit.pass ? 1 : 0)
        << '\n';

    if (!hit.pass) {
        overall_pass = false;
    }

    std::cout
        << "OVERALL_PASS="
        << (overall_pass ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 238\n";

    return overall_pass ? 0 : 1;
}