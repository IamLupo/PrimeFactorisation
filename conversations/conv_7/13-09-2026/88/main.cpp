#include <cstdint>
#include <iostream>
#include <random>
#include <limits>

using u64 = std::uint64_t;

struct Digits {
    u64 value;
    u64 p;
    u64 digit[64];
    unsigned size;
};

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
    u64 rank;
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

Digits make_digits(u64 n, u64 p)
{
    Digits d{};
    d.value = n;
    d.p = p;
    d.size = 0;

    if (n == 0) {
        d.digit[0] = 0;
        d.size = 1;
        return d;
    }

    while (n > 0) {
        d.digit[d.size++] = n % p;
        n /= p;
    }

    return d;
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

    if (!safe_pow_u64(
            p,
            a0,
            p_a)) {
        return false;
    }

    if (mul_overflow_u64(
            static_cast<u64>(b + 1),
            p_a)) {
        return false;
    }

    const u64 s0 =
        static_cast<u64>(b + 1) *
        p_a;

    u64 p_e;

    if (!safe_pow_u64(
            p,
            e,
            p_e)) {
        return false;
    }

    if (mul_overflow_u64(
            q,
            p_e)) {
        return false;
    }

    const u64 q_pe =
        q * p_e;

    if (add_overflow_u64(
            s0,
            q_pe)) {
        return false;
    }

    const u64 m_plus_one =
        s0 + q_pe;

    if (m_plus_one == 0) {
        return false;
    }

    out = {
        p,
        a0,
        b,
        z,
        e,
        s0,
        q,
        p_e,
        m_plus_one - 1
    };

    return true;
}

u64 rank_from_digits(
    const Digits& j,
    const Digits& q)
{
    u64 rank = 0;
    u64 weight = 1;

    const unsigned R =
        (j.size > q.size)
            ? j.size
            : q.size;

    for (unsigned r = 0;
         r < R;
         ++r) {

        const u64 jd =
            (r < j.size)
                ? j.digit[r]
                : 0;

        rank +=
            jd * weight;

        const u64 qd =
            (r < q.size)
                ? q.digit[r]
                : 0;

        weight *= qd + 1;
    }

    return rank;
}

u64 rank_j(
    u64 j,
    u64 q,
    u64 p)
{
    const Digits jd =
        make_digits(j, p);

    const Digits qd =
        make_digits(q, p);

    return rank_from_digits(
        jd,
        qd);
}

u64 unrank_k(
    u64 k,
    const Digits& q)
{
    u64 remaining = k;
    u64 place = 1;
    u64 j = 0;

    for (unsigned r = 0;
         r < q.size;
         ++r) {

        const u64 radix =
            q.digit[r] + 1;

        const u64 jd =
            remaining % radix;

        j +=
            jd * place;

        remaining /=
            radix;

        place *= q.p;
    }

    return j;
}

bool find_interval_type(
    u64 j,
    u64 q,
    u64 p,
    unsigned& r_out)
{
    const Digits jd =
        make_digits(j, p);

    const Digits qd =
        make_digits(q, p);

    const unsigned R =
        (jd.size > qd.size)
            ? jd.size
            : qd.size;

    for (unsigned r = 0;
         r < R;
         ++r) {

        const u64 j_digit =
            (r < jd.size)
                ? jd.digit[r]
                : 0;

        const u64 q_digit =
            (r < qd.size)
                ? qd.digit[r]
                : 0;

        if (j_digit < q_digit) {
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

    const u64 p_r =
        pow_u64(p, r);

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
            p_r)) {
        return false;
    }

    next =
        multiplier * p_r;

    r_out = r;

    return true;
}

/*
 * Largest admissible j <= t.
 *
 * Find the most significant digit h with
 *
 *     t_h > q_h.
 *
 * Then:
 *
 *     j_i=t_i  for i>h
 *     j_i=q_i  for i<=h
 *
 * If no violation exists, j=t.
 */
u64 largest_admissible_leq(
    u64 t,
    u64 q,
    u64 p)
{
    const Digits td =
        make_digits(t, p);

    const Digits qd =
        make_digits(q, p);

    const unsigned R =
        (td.size > qd.size)
            ? td.size
            : qd.size;

    bool violation = false;
    unsigned h = 0;

    for (int r =
             static_cast<int>(R) - 1;
         r >= 0;
         --r) {

        const unsigned rr =
            static_cast<unsigned>(r);

        const u64 t_digit =
            (rr < td.size)
                ? td.digit[rr]
                : 0;

        const u64 q_digit =
            (rr < qd.size)
                ? qd.digit[rr]
                : 0;

        if (t_digit > q_digit) {
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
         r < R;
         ++r) {

        const u64 t_digit =
            (r < td.size)
                ? td.digit[r]
                : 0;

        const u64 q_digit =
            (r < qd.size)
                ? qd.digit[r]
                : 0;

        const u64 j_digit =
            (r > h)
                ? t_digit
                : q_digit;

        result +=
            j_digit * place;

        place *= p;
    }

    return result;
}

bool lucas_hit(
    u64 n,
    u64 m,
    u64 p)
{
    const Digits nd =
        make_digits(n, p);

    const Digits md =
        make_digits(m, p);

    const unsigned R =
        (nd.size > md.size)
            ? nd.size
            : md.size;

    for (unsigned r = 0;
         r < R;
         ++r) {

        const u64 n_digit =
            (r < nd.size)
                ? nd.digit[r]
                : 0;

        const u64 m_digit =
            (r < md.size)
                ? md.digit[r]
                : 0;

        if (n_digit > m_digit) {
            return true;
        }
    }

    return false;
}

/*
 * Complete n -> t -> j -> interval -> rank localization.
 *
 * IMPORTANT:
 *
 * Finding j alone does NOT mean n is inside the interval.
 * We additionally require
 *
 *     n < successor(j) * p^e.
 */
bool localize_n(
    u64 n,
    const StructuralCase& sc,
    Localization& out)
{
    out = {
        false,
        0,
        0,
        0,
        0,
        0,
        0
    };

    if (n < sc.s0) {
        return true;
    }

    const u64 relative =
        n - sc.s0;

    const u64 t =
        relative / sc.p_e;

    out.t = t;

    /*
     * Final MISS region.
     */
    if (t >= sc.q) {
        return true;
    }

    const u64 j =
        largest_admissible_leq(
            t,
            sc.q,
            sc.p);

    out.j = j;

    /*
     * Terminal tuple is not a HIT interval.
     */
    if (j >= sc.q) {
        return true;
    }

    unsigned r;

    if (!find_interval_type(
            j,
            sc.q,
            sc.p,
            r)) {
        return true;
    }

    out.r = r;

    const u64 x =
        sc.s0 +
        j * sc.p_e;

    out.x = x;

    u64 next;
    unsigned successor_r;

    if (!successor_j(
            j,
            sc.q,
            sc.p,
            next,
            successor_r)) {
        return true;
    }

    if (successor_r != r) {
        return true;
    }

    if (mul_overflow_u64(
            next,
            sc.p_e)) {
        return true;
    }

    const u64 end_plus_one =
        next * sc.p_e;

    const u64 end =
        end_plus_one - 1;

    out.end = end;

    /*
     * Critical interval-membership test.
     */
    if (n < x ||
        n >= end_plus_one) {
        return true;
    }

    out.rank =
        rank_j(
            j,
            sc.q,
            sc.p);

    out.hit = true;

    return true;
}

bool verify_point(
    u64 n,
    const StructuralCase& sc,
    bool verbose)
{
    Localization loc;

    if (!localize_n(
            n,
            sc,
            loc)) {
        return false;
    }

    const bool expected =
        lucas_hit(
            n,
            sc.m,
            sc.p);

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
                << " rank=" << loc.rank
                << " r=" << loc.r
                << " x=" << loc.x
                << " end=" << loc.end
                << " s0=" << sc.s0
                << " e=" << sc.e
                << " m=" << sc.m
                << '\n';
        }

        return false;
    }

    if (!expected) {
        return true;
    }

    /*
     * Recovered j.
     */
    const u64 expected_j =
        largest_admissible_leq(
            loc.t,
            sc.q,
            sc.p);

    if (loc.j != expected_j) {
        return false;
    }

    if (loc.j >= sc.q) {
        return false;
    }

    /*
     * Rank.
     */
    const u64 expected_rank =
        rank_j(
            loc.j,
            sc.q,
            sc.p);

    if (loc.rank != expected_rank) {
        return false;
    }

    /*
     * Rank round trip.
     */
    const Digits qd =
        make_digits(
            sc.q,
            sc.p);

    if (unrank_k(
            loc.rank,
            qd) != loc.j) {
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

    if (rank_j(
            next_j,
            sc.q,
            sc.p) !=
        loc.rank + 1) {
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

    return true;
}

/*
 * Regression cases using the ORIGINAL structural parameters.
 *
 * This is the important correction to Experiment 242.
 */
bool regression_tests()
{
    const struct {
        u64 p;
        unsigned a0;
        unsigned b;
        unsigned z;
        u64 q;
        u64 n;
    } cases[] = {
        {2, 0, 0, 0, 1, 2},

        {2, 0, 0, 0, 3, 2},

        {2, 1, 0, 0, 7, 5},
        {2, 1, 0, 0, 7, 29},

        {3, 0, 0, 0, 2, 3},
        {3, 0, 0, 0, 2, 6},

        {5, 0, 0, 0, 4, 5},
        {5, 0, 0, 0, 4, 20},

        {5, 3, 1, 4, 100, 62574},

        {2, 5, 0, 2,
         1073741825ULL,
         274877907231ULL},

        {3, 4, 1, 5,
         987654321ULL,
         58320000000890ULL},

        {5, 3, 1, 4,
         1000007654321ULL,
         390627989969140874ULL}
    };

    bool all_pass = true;

    std::cout
        << "REGRESSION TESTS\n";

    for (const auto& c : cases) {
        StructuralCase sc;

        if (!build_structural_case(
                c.p,
                c.a0,
                c.b,
                c.z,
                c.q,
                sc)) {

            all_pass = false;
            continue;
        }

        const bool pass =
            verify_point(
                c.n,
                sc,
                true);

        std::cout
            << "p=" << c.p
            << " a0=" << c.a0
            << " b=" << c.b
            << " z=" << c.z
            << " q=" << c.q
            << " n=" << c.n
            << " e=" << sc.e
            << " s0=" << sc.s0
            << " m=" << sc.m
            << " pass="
            << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            all_pass = false;
        }
    }

    return all_pass;
}

Result exhaustive_reference()
{
    u64 cases = 0;
    u64 failures = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (u64 q = 1;
             q <= 200;
             ++q) {

            StructuralCase sc;

            if (!build_structural_case(
                    p,
                    0,
                    0,
                    0,
                    q,
                    sc)) {
                continue;
            }

            if (sc.m > 200000) {
                continue;
            }

            for (u64 n = 0;
                 n <= sc.m;
                 ++n) {

                ++cases;

                if (!verify_point(
                        n,
                        sc,
                        false)) {

                    ++failures;

                    if (failures <= 5) {
                        std::cout
                            << "FAIL REFERENCE"
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
        }
    }

    return {
        failures == 0,
        cases,
        failures
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
            rng() %
            (sc.m + 1);

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

        const Digits qd =
            make_digits(q, p);

        /*
         * Random admissible j.
         */
        u64 j = 0;
        u64 place = 1;

        for (unsigned r = 0;
             r < qd.size;
             ++r) {

            const u64 jd =
                rng() %
                (qd.digit[r] + 1);

            j +=
                jd * place;

            place *= p;
        }

        /*
         * Terminal tuple q is not a HIT interval.
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

int main()
{
    std::cout
        << "START EXPERIMENT 243\n";

    bool overall_pass = true;

    /*
     * Exact regressions from 241/242 with correct structure.
     */
    const bool regression_pass =
        regression_tests();

    std::cout
        << "regression_pass="
        << (regression_pass ? 1 : 0)
        << '\n';

    if (!regression_pass) {
        overall_pass = false;
    }

    /*
     * Complete small reference.
     */
    const Result reference =
        exhaustive_reference();

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
     * Large arbitrary-n tests.
     */
    const Result random =
        random_tests(
            243243243ULL,
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
     * Random HIT tests.
     */
    const Result hit =
        random_hit_tests(
            943243243ULL,
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
        << "FINISHED EXPERIMENT 243\n";

    return overall_pass ? 0 : 1;
}