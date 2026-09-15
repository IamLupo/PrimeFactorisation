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
 * Largest admissible j <= t.
 *
 * Find the most significant digit h where t_h > q_h.
 *
 * If none exists:
 *
 *     j=t.
 *
 * Otherwise:
 *
 *     j_i = t_i for i>h
 *     j_i = q_i for i<=h.
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

    bool violation = false;
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
            violation = true;
            h = rr;
            break;
        }
    }

    if (!violation) {
        return t;
    }

    u64 result = 0;
    u64 place = 1;

    for (unsigned r = 0;
         r <= R;
         ++r) {

        const u64 d =
            (r > h)
                ? digit_at(t, p, r)
                : digit_at(q, p, r);

        result += d * place;

        if (r < R) {
            place *= p;
        }
    }

    return result;
}

bool find_interval_type(
    u64 j,
    u64 q,
    u64 p,
    unsigned& r_out)
{
    const unsigned R =
        highest_digit_position(q, p);

    for (unsigned r = 0;
         r <= R;
         ++r) {

        if (digit_at(j, p, r) <
            digit_at(q, p, r)) {

            r_out = r;
            return true;
        }
    }

    return false;
}

/*
 * p-adic / mixed-radix successor:
 *
 * next(j) =
 *     (floor(j/p^r)+1)p^r
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

    next =
        (quotient + 1) * pr;

    r_out = r;

    return true;
}

/*
 * Mixed-radix rank:
 *
 *     rank(j)
 *       = sum_i j_i W_i
 *
 *     W_i = product_{h<i}(q_h+1)
 *
 * This is the rank among admissible tuples.
 */
u64 rank_j(
    u64 j,
    u64 q,
    u64 p)
{
    const unsigned R =
        highest_digit_position(q, p);

    u64 rank = 0;
    u64 weight = 1;

    for (unsigned r = 0;
         r <= R;
         ++r) {

        const u64 jd =
            digit_at(j, p, r);

        rank +=
            jd * weight;

        const u64 qd =
            digit_at(q, p, r);

        weight *=
            qd + 1;
    }

    return rank;
}

/*
 * Mixed-radix unrank:
 *
 * Given k, recover the unique admissible j.
 */
u64 unrank_j(
    u64 k,
    u64 q,
    u64 p)
{
    const unsigned R =
        highest_digit_position(q, p);

    u64 remaining = k;
    u64 j = 0;
    u64 place = 1;

    for (unsigned r = 0;
         r <= R;
         ++r) {

        const u64 qd =
            digit_at(q, p, r);

        const u64 radix =
            qd + 1;

        const u64 jd =
            remaining % radix;

        j +=
            jd * place;

        remaining /=
            radix;

        place *= p;
    }

    return j;
}

/*
 * Complete interval localization.
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

    /*
     * Block coordinate.
     */
    const u64 t =
        relative / sc.p_e;

    result.t = t;

    /*
     * Final MISS block.
     */
    if (t >= sc.q) {
        return result;
    }

    /*
     * Largest admissible interval label <= t.
     */
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

/*
 * Verify one n and return its complete rank information.
 */
bool verify_point(
    u64 n,
    const StructuralCase& sc,
    bool verbose)
{
    const bool expected_hit =
        lucas_hit(
            n,
            sc.m,
            sc.p);

    const Localization loc =
        localize_n(
            n,
            sc);

    if (loc.hit != expected_hit) {
        if (verbose) {
            std::cout
                << "FAIL CLASS"
                << " n=" << n
                << " expected="
                << (expected_hit ? 1 : 0)
                << " actual="
                << (loc.hit ? 1 : 0)
                << " t=" << loc.t
                << " j=" << loc.j
                << '\n';
        }

        return false;
    }

    /*
     * MISS.
     */
    if (!expected_hit) {
        return true;
    }

    /*
     * HIT must have a valid interval.
     */
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
        return false;
    }

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
     * Rank.
     */
    const u64 k =
        rank_j(
            loc.j,
            sc.q,
            sc.p);

    /*
     * Round trip rank -> j.
     */
    const u64 round_trip_j =
        unrank_j(
            k,
            sc.q,
            sc.p);

    if (round_trip_j != loc.j) {
        if (verbose) {
            std::cout
                << "FAIL RANK ROUND TRIP"
                << " n=" << n
                << " j=" << loc.j
                << " k=" << k
                << " round_trip="
                << round_trip_j
                << '\n';
        }

        return false;
    }

    /*
     * Interval type.
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
     * x_j.
     */
    const u64 expected_x =
        sc.s0 +
        loc.j * sc.p_e;

    if (loc.x != expected_x) {
        return false;
    }

    /*
     * Successor.
     */
    u64 next_j;
    unsigned successor_r;

    if (!successor_j(
            loc.j,
            sc.q,
            sc.p,
            next_j,
            successor_r)) {
        return false;
    }

    if (successor_r != loc.r) {
        return false;
    }

    /*
     * Rank successor.
     */
    const u64 successor_rank =
        rank_j(
            next_j,
            sc.q,
            sc.p);

    if (successor_rank != k + 1) {
        if (verbose) {
            std::cout
                << "FAIL SUCCESSOR RANK"
                << " n=" << n
                << " j=" << loc.j
                << " k=" << k
                << " next=" << next_j
                << " next_rank="
                << successor_rank
                << '\n';
        }

        return false;
    }

    /*
     * Endpoint.
     */
    const u64 expected_end =
        next_j * sc.p_e - 1;

    if (loc.end != expected_end) {
        return false;
    }

    if (n < loc.x ||
        n > loc.end) {
        return false;
    }

    /*
     * Every n inside this interval has exactly the same rank.
     */
    const u64 reconstructed_rank =
        rank_j(
            loc.j,
            sc.q,
            sc.p);

    if (reconstructed_rank != k) {
        return false;
    }

    return true;
}

/*
 * Directly test the rank/unrank bijection for bounded q.
 *
 * This is exhaustive over every rank, but never scans the integer
 * domain [0,m].
 */
bool exhaustive_rank_bijection(
    const StructuralCase& sc,
    bool verbose)
{
    u64 interval_count = 1;

    const unsigned R =
        highest_digit_position(
            sc.q,
            sc.p);

    for (unsigned r = 0;
         r <= R;
         ++r) {

        interval_count *=
            digit_at(
                sc.q,
                sc.p,
                r) + 1;
    }

    /*
     * q itself is the terminal tuple, so HIT intervals are
     *
     * 0,...,interval_count-2.
     */
    const u64 hit_count =
        interval_count - 1;

    for (u64 k = 0;
         k < hit_count;
         ++k) {

        const u64 j =
            unrank_j(
                k,
                sc.q,
                sc.p);

        const u64 recovered_rank =
            rank_j(
                j,
                sc.q,
                sc.p);

        if (recovered_rank != k) {
            if (verbose) {
                std::cout
                    << "FAIL UNRANK/RANK"
                    << " k=" << k
                    << " j=" << j
                    << " recovered="
                    << recovered_rank
                    << '\n';
            }

            return false;
        }

        /*
         * Successor rank.
         */
        u64 next_j;
        unsigned r;

        if (!successor_j(
                j,
                sc.q,
                sc.p,
                next_j,
                r)) {
            return false;
        }

        const u64 next_rank =
            rank_j(
                next_j,
                sc.q,
                sc.p);

        if (next_rank != k + 1) {
            if (verbose) {
                std::cout
                    << "FAIL RANK SUCCESSOR"
                    << " k=" << k
                    << " j=" << j
                    << " next=" << next_j
                    << " next_rank="
                    << next_rank
                    << '\n';
            }

            return false;
        }

        /*
         * For the final HIT interval, next_j must equal q.
         */
        if (k + 1 == hit_count) {
            if (next_j != sc.q) {
                return false;
            }

            const u64 terminal_rank =
                rank_j(
                    sc.q,
                    sc.q,
                    sc.p);

            if (terminal_rank != hit_count) {
                if (verbose) {
                    std::cout
                        << "FAIL TERMINAL RANK"
                        << " q=" << sc.q
                        << " rank="
                        << terminal_rank
                        << " expected="
                        << hit_count
                        << '\n';
                }

                return false;
            }
        }
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
        << "DETERMINISTIC RANK LOCALIZATION\n";

    for (const auto& c : tests) {
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

        ++cases;

        bool pass =
            exhaustive_rank_bijection(
                sc,
                true);

        /*
         * A few direct n points.
         */
        const u64 points[] = {
            0,
            sc.s0,
            sc.s0 + sc.p_e - 1,
            sc.m - sc.s0,
            sc.m
        };

        for (u64 n : points) {
            if (n <= sc.m) {
                if (!verify_point(
                        n,
                        sc,
                        true)) {
                    pass = false;
                }
            }
        }

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
        failures
    };
}

/*
 * Random large-n localization plus rank verification.
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

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 p =
            primes[rng() % 3];

        const unsigned a0 =
            static_cast<unsigned>(
                rng() % 8);

        const unsigned z =
            static_cast<unsigned>(
                rng() % 8);

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

/*
 * Random HIT points constructed from admissible j.
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
                rng() % 8);

        const unsigned z =
            static_cast<unsigned>(
                rng() % 8);

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
         * Generate random admissible j.
         */
        u64 j = 0;
        u64 current = q;
        u64 place = 1;

        while (current > 0) {
            const u64 qdigit =
                current % p;

            const u64 jd =
                rng() % (qdigit + 1);

            j +=
                jd * place;

            current /= p;
            place *= p;
        }

        /*
         * Exclude terminal j=q.
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

/*
 * Test the extremely sparse cases that are especially useful for
 * exposing digit-order mistakes.
 */
bool sparse_tests()
{
    bool all_pass = true;

    const struct {
        u64 p;
        u64 q;
    } cases[] = {
        {2, 3},
        {2, 5},
        {2, 9},
        {2, 17},
        {2, 65},
        {2, 257},
        {2, 1025},
        {2, 1048577},
        {2, 1073741825ULL},

        {3, 10},
        {3, 28},
        {3, 82},

        {5, 26},
        {5, 126},
        {5, 626}
    };

    std::cout
        << "SPARSE STRUCTURAL TESTS\n";

    for (const auto& c : cases) {
        /*
         * q = p^R + 1 style values are intentional.
         */
        const unsigned a0 = 0;
        const unsigned b = 0;
        const unsigned z = 0;

        StructuralCase sc;

        if (!build_structural_case(
                c.p,
                a0,
                b,
                z,
                c.q,
                sc)) {
            all_pass = false;
            continue;
        }

        /*
         * Test points exactly at large block boundaries.
         */
        bool pass = true;

        for (unsigned shift = 0;
             shift <= highest_digit_position(
                          c.q,
                          c.p);
             ++shift) {

            const u64 block =
                pow_u64(
                    c.p,
                    shift);

            if (block > c.q) {
                break;
            }

            const u64 base =
                sc.s0 +
                block * sc.p_e;

            if (base <= sc.m) {
                if (!verify_point(
                        base,
                        sc,
                        true)) {
                    pass = false;
                }

                if (base + sc.p_e - 1 <= sc.m) {
                    if (!verify_point(
                            base + sc.p_e - 1,
                            sc,
                            true)) {
                        pass = false;
                    }
                }
            }
        }

        std::cout
            << "p=" << c.p
            << " q=" << c.q
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
        << "START EXPERIMENT 240\n";

    bool overall_pass = true;

    /*
     * Deterministic rank/unrank and localization.
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
     * Sparse digit-pattern stress.
     */
    const bool sparse_pass =
        sparse_tests();

    std::cout
        << "sparse_pass="
        << (sparse_pass ? 1 : 0)
        << '\n';

    if (!sparse_pass) {
        overall_pass = false;
    }

    /*
     * Arbitrary n.
     */
    const Result random =
        random_tests(
            240240240ULL,
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
     * Random HITs.
     */
    const Result hit =
        random_hit_tests(
            940240240ULL,
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
        << "FINISHED EXPERIMENT 240\n";

    return overall_pass ? 0 : 1;
}
