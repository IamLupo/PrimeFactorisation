#include <cstdint>
#include <iostream>
#include <random>
#include <limits>

using u64 = std::uint64_t;

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
 * Correct largest-admissible <= t.
 *
 * Let h be the most significant digit where
 *
 *     t_h > q_h.
 *
 * Then:
 *
 *   higher digits = t_i
 *   digit h       = q_h
 *   lower digits  = q_i
 *
 * If no such h exists, t itself is admissible.
 *
 * This is the correct lexicographic construction.
 */
u64 largest_admissible_leq(
    u64 t,
    u64 q,
    u64 p)
{
    const unsigned Rt =
        highest_digit_position(t, p);

    const unsigned Rq =
        highest_digit_position(q, p);

    const unsigned R =
        (Rt > Rq) ? Rt : Rq;

    /*
     * Find the MOST significant violation.
     */
    bool found_violation = false;
    unsigned h = 0;

    for (int r = static_cast<int>(R);
         r >= 0;
         --r) {

        const unsigned rr =
            static_cast<unsigned>(r);

        const u64 td =
            digit_at(t, p, rr);

        const u64 qd =
            digit_at(q, p, rr);

        if (td > qd) {
            found_violation = true;
            h = rr;
            break;
        }
    }

    /*
     * No violation: t itself is admissible.
     */
    if (!found_violation) {
        return t;
    }

    u64 result = 0;
    u64 place = 1;

    for (unsigned r = 0; r <= R; ++r) {
        u64 d;

        if (r > h) {
            d = digit_at(t, p, r);
        } else {
            /*
             * At and below h, maximize subject to <=_p q.
             */
            d = digit_at(q, p, r);
        }

        result += d * place;

        if (r < R) {
            place *= p;
        }
    }

    return result;
}

/*
 * Exhaustive small-domain reference:
 *
 * Find the actual largest j satisfying
 *
 *     j <= t
 *     j <=_p q
 *
 * by scanning. This is ONLY used for small q/t cases to validate
 * the fast localization rule itself.
 */
u64 brute_largest_admissible(
    u64 t,
    u64 q,
    u64 p)
{
    for (u64 j = t;; --j) {
        if (digitwise_leq(j, q, p)) {
            return j;
        }

        if (j == 0) {
            break;
        }
    }

    /*
     * 0 is always admissible.
     */
    return 0;
}

bool find_interval_type(
    u64 j,
    u64 q,
    u64 p,
    unsigned& r_out)
{
    const unsigned R =
        highest_digit_position(q, p);

    for (unsigned r = 0; r <= R; ++r) {
        if (digit_at(j, p, r) <
            digit_at(q, p, r)) {

            r_out = r;
            return true;
        }
    }

    return false;
}

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

    next =
        (quotient + 1) * pr;

    r_out = r;

    return true;
}

/*
 * Corrected interval localization.
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

    result.t = t;

    /*
     * Final MISS block:
     *
     * t >= q.
     */
    if (t >= sc.q) {
        return result;
    }

    const u64 j =
        largest_admissible_leq(
            t,
            sc.q,
            sc.p);

    result.j = j;

    unsigned r;

    if (!find_interval_type(
            j,
            sc.q,
            sc.p,
            r)) {
        return result;
    }

    result.r = r;

    const u64 x =
        sc.s0 +
        j * sc.p_e;

    result.x = x;

    u64 next;
    unsigned next_r;

    if (!successor_j(
            j,
            sc.q,
            sc.p,
            next,
            next_r)) {
        return result;
    }

    if (next_r != r) {
        return result;
    }

    const u64 end =
        next * sc.p_e - 1;

    result.end = end;

    result.hit =
        n >= x &&
        n <= end;

    return result;
}

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
        localize_n(n, sc);

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

    if (!expected) {
        return true;
    }

    /*
     * HIT: verify t -> j independently against the brute
     * reference whenever t is small enough.
     */
    if (loc.t <= 100000) {
        const u64 brute =
            brute_largest_admissible(
                loc.t,
                sc.q,
                sc.p);

        if (loc.j != brute) {
            if (verbose) {
                std::cout
                    << "FAIL BRUTE J"
                    << " n=" << n
                    << " t=" << loc.t
                    << " j=" << loc.j
                    << " brute=" << brute
                    << '\n';
            }

            return false;
        }
    }

    /*
     * j must itself be admissible.
     */
    if (!digitwise_leq(
            loc.j,
            sc.q,
            sc.p)) {
        return false;
    }

    if (loc.j > loc.t) {
        return false;
    }

    /*
     * r must be the first digit where j_r < q_r.
     */
    unsigned expected_r;

    if (!find_interval_type(
            loc.j,
            sc.q,
            sc.p,
            expected_r)) {
        return false;
    }

    if (expected_r != loc.r) {
        return false;
    }

    /*
     * Interval reconstruction.
     */
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

    return true;
}

bool sparse_counterexample()
{
    /*
     * This is the exact family that killed Experiment 238:
     *
     *     q = 2^30 + 1
     *
     * and t=2^29.
     *
     * Correct answer is j=1, not j=0.
     */
    const u64 p = 2;
    const u64 q =
        1073741825ULL;

    const u64 t =
        536870912ULL;

    const u64 expected = 1;

    const u64 actual =
        largest_admissible_leq(
            t,
            q,
            p);

    std::cout
        << "SPARSE COUNTEREXAMPLE"
        << " p=" << p
        << " q=" << q
        << " t=" << t
        << " j=" << actual
        << " expected=" << expected
        << " pass="
        << (actual == expected ? 1 : 0)
        << '\n';

    return actual == expected;
}

Result brute_reference_tests()
{
    u64 failures = 0;
    u64 cases = 0;

    /*
     * Small exhaustive q,t reference comparison.
     */
    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (u64 q = 1; q <= 500; ++q) {
            for (u64 t = 0; t <= 1000; ++t) {
                ++cases;

                const u64 fast =
                    largest_admissible_leq(
                        t,
                        q,
                        p);

                const u64 brute =
                    brute_largest_admissible(
                        t,
                        q,
                        p);

                if (fast != brute) {
                    ++failures;

                    if (failures <= 5) {
                        std::cout
                            << "FAIL BRUTE REFERENCE"
                            << " p=" << p
                            << " q=" << q
                            << " t=" << t
                            << " fast=" << fast
                            << " brute=" << brute
                            << '\n';
                    }
                }
            }
        }
    }

    return {
        failures == 0,
        cases,
        failures
    };
}

Result random_localization_tests(
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

        const u64 n =
            rng() % (sc.m + 1);

        if (!verify_point(
                n,
                sc,
                false)) {

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
        failures
    };
}

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
         * Generate a random admissible j.
         */
        u64 j = 0;
        u64 place = 1;
        u64 current = q;

        while (current > 0) {
            const u64 qdigit =
                current % p;

            const u64 jd =
                rng() % (qdigit + 1);

            j += jd * place;

            current /= p;
            place *= p;
        }

        /*
         * Terminal j=q is excluded.
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
                    << "FAIL RANDOM HIT"
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
        failures
    };
}

bool deterministic_localization_tests()
{
    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
    } cases[] = {
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

    bool all_pass = true;

    for (const auto& c : cases) {
        StructuralCase sc;

        if (!build_structural_case(
                c.p,
                c.a0,
                c.b,
                c.z,
                c.q,
                sc)) {
            continue;
        }

        bool pass = true;

        /*
         * Initial and final boundary points.
         */
        if (sc.s0 > 0) {
            pass &=
                verify_point(
                    0,
                    sc,
                    true);

            pass &=
                verify_point(
                    sc.s0 - 1,
                    sc,
                    true);
        }

        /*
         * Deliberately sample points based on large block
         * coordinates.
         */
        const u64 candidates[] = {
            sc.s0,
            sc.s0 + sc.p_e - 1,
            sc.s0 + 2 * sc.p_e,
            sc.m - sc.s0,
            sc.m - sc.s0 + 1,
            sc.m
        };

        for (u64 n : candidates) {
            if (n <= sc.m) {
                pass &=
                    verify_point(
                        n,
                        sc,
                        true);
            }
        }

        std::cout
            << "CASE"
            << " p=" << c.p
            << " q=" << c.q
            << " m=" << sc.m
            << " pass=" << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            all_pass = false;
        }
    }

    return all_pass;
}

int main()
{
    std::cout
        << "START EXPERIMENT 239\n";

    bool overall_pass = true;

    /*
     * Exact counterexample from Experiment 238.
     */
    const bool sparse_pass =
        sparse_counterexample();

    if (!sparse_pass) {
        overall_pass = false;
    }

    /*
     * Exhaustive small reference:
     *
     * fast largest-admissible(t)
     * versus brute-force largest-admissible(t).
     */
    const Result reference =
        brute_reference_tests();

    std::cout
        << "reference_cases="
        << reference.cases
        << " reference_failures="
        << reference.failures
        << " reference_pass="
        << (reference.pass ? 1 : 0)
        << '\n';

    if (!reference.pass) {
        overall_pass = false;
    }

    /*
     * Deterministic large structural cases.
     */
    const bool deterministic_pass =
        deterministic_localization_tests();

    std::cout
        << "deterministic_pass="
        << (deterministic_pass ? 1 : 0)
        << '\n';

    if (!deterministic_pass) {
        overall_pass = false;
    }

    /*
     * Random arbitrary n.
     */
    const Result random =
        random_localization_tests(
            239239239ULL,
            100000);

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
     * Random points deliberately chosen from known HIT intervals.
     */
    const Result hit =
        random_hit_tests(
            939239239ULL,
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
        << "FINISHED EXPERIMENT 239\n";

    return overall_pass ? 0 : 1;
}
