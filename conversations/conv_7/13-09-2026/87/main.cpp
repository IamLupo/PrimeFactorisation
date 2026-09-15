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

    for (unsigned r = 0; r < R; ++r) {
        const u64 jd =
            (r < j.size)
                ? j.digit[r]
                : 0;

        rank += jd * weight;

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

    return rank_from_digits(jd, qd);
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

        j += jd * place;

        remaining /= radix;
        place *= q.p;
    }

    return j;
}

/*
 * Find the first digit r where j_r < q_r.
 *
 * This is the interval type.
 */
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

/*
 * p-adic successor of an admissible j.
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
 * Find the MOST significant h where t_h > q_h.
 *
 * If none exists, j=t.
 *
 * Otherwise:
 *
 *   j_i=t_i for i>h
 *   j_i=q_i for i<=h
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

        result += j_digit * place;
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
 * Complete direct rank localization:
 *
 *     n -> t -> j -> rank
 */
bool direct_rank_from_n(
    u64 n,
    const StructuralCase& sc,
    u64& rank_out,
    u64& j_out,
    u64& t_out,
    unsigned& h_out)
{
    if (n < sc.s0) {
        return false;
    }

    const u64 relative =
        n - sc.s0;

    const u64 t =
        relative / sc.p_e;

    if (t >= sc.q) {
        return false;
    }

    const Digits qd =
        make_digits(sc.q, sc.p);

    const Digits td =
        make_digits(t, sc.p);

    const unsigned R =
        (td.size > qd.size)
            ? td.size
            : qd.size;

    bool violation = false;
    unsigned h = 0;

    /*
     * Find most significant violation.
     */
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

    u64 j = 0;
    u64 rank = 0;
    u64 place = 1;
    u64 weight = 1;

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
            (!violation)
                ? t_digit
                : (r > h
                       ? t_digit
                       : q_digit);

        j += j_digit * place;
        rank += j_digit * weight;

        place *= sc.p;
        weight *= q_digit + 1;
    }

    /*
     * j=q is the terminal tuple, not a HIT interval.
     */
    if (j >= sc.q) {
        return false;
    }

    j_out = j;
    rank_out = rank;
    t_out = t;
    h_out = violation ? h : R;

    return true;
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

    u64 direct_rank = 0;
    u64 direct_j = 0;
    u64 t = 0;
    unsigned h = 0;

    const bool localized =
        direct_rank_from_n(
            n,
            sc,
            direct_rank,
            direct_j,
            t,
            h);

    if (localized != expected) {
        if (verbose) {
            std::cout
                << "FAIL LOCALIZATION"
                << " n=" << n
                << " expected="
                << (expected ? 1 : 0)
                << " localized="
                << (localized ? 1 : 0)
                << " t=" << t
                << " j=" << direct_j
                << '\n';
        }

        return false;
    }

    if (!expected) {
        return true;
    }

    /*
     * Independent rank calculation from j.
     */
    const u64 traditional_rank =
        rank_j(
            direct_j,
            sc.q,
            sc.p);

    if (traditional_rank != direct_rank) {
        if (verbose) {
            std::cout
                << "FAIL DIRECT RANK"
                << " n=" << n
                << " t=" << t
                << " j=" << direct_j
                << " direct="
                << direct_rank
                << " traditional="
                << traditional_rank
                << '\n';
        }

        return false;
    }

    /*
     * Rank -> j round trip.
     */
    const Digits qd =
        make_digits(
            sc.q,
            sc.p);

    const u64 round_trip =
        unrank_k(
            direct_rank,
            qd);

    if (round_trip != direct_j) {
        if (verbose) {
            std::cout
                << "FAIL ROUND TRIP"
                << " n=" << n
                << " rank="
                << direct_rank
                << " j="
                << direct_j
                << " round_trip="
                << round_trip
                << '\n';
        }

        return false;
    }

    /*
     * Successor rank property.
     */
    u64 next_j;
    unsigned r;

    if (!successor_j(
            direct_j,
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

    if (next_rank != direct_rank + 1) {
        if (verbose) {
            std::cout
                << "FAIL SUCCESSOR RANK"
                << " n=" << n
                << " j=" << direct_j
                << " rank="
                << direct_rank
                << " next="
                << next_j
                << " next_rank="
                << next_rank
                << '\n';
        }

        return false;
    }

    /*
     * Direct j must be admissible.
     */
    const Digits jd =
        make_digits(
            direct_j,
            sc.p);

    for (unsigned i = 0;
         i < qd.size;
         ++i) {

        const u64 jdigit =
            (i < jd.size)
                ? jd.digit[i]
                : 0;

        if (jdigit > qd.digit[i]) {
            return false;
        }
    }

    /*
     * And j <= t.
     */
    if (direct_j > t) {
        return false;
    }

    return true;
}

u64 brute_largest_admissible(
    u64 t,
    u64 q,
    u64 p)
{
    const Digits qd =
        make_digits(q, p);

    for (u64 j = t;; --j) {
        const Digits jd =
            make_digits(j, p);

        bool admissible = true;

        const unsigned R =
            (jd.size > qd.size)
                ? jd.size
                : qd.size;

        for (unsigned r = 0;
             r < R;
             ++r) {

            const u64 jdigit =
                (r < jd.size)
                    ? jd.digit[r]
                    : 0;

            const u64 qdigit =
                (r < qd.size)
                    ? qd.digit[r]
                    : 0;

            if (jdigit > qdigit) {
                admissible = false;
                break;
            }
        }

        if (admissible) {
            return j;
        }

        if (j == 0) {
            break;
        }
    }

    return 0;
}

Result exhaustive_reference()
{
    u64 cases = 0;
    u64 failures = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (u64 q = 1;
             q <= 500;
             ++q) {

            const Digits qd =
                make_digits(q, p);

            for (u64 t = 0;
                 t <= 1000;
                 ++t) {

                ++cases;

                u64 j;
                u64 rank;
                unsigned h;

                /*
                 * Direct rank routine requires q-domain HIT
                 * semantics, so only compare its j machinery here.
                 */
                const Digits td =
                    make_digits(t, p);

                const unsigned R =
                    (td.size > qd.size)
                        ? td.size
                        : qd.size;

                bool violation = false;
                unsigned highest = 0;

                for (int r =
                         static_cast<int>(R) - 1;
                     r >= 0;
                     --r) {

                    const unsigned rr =
                        static_cast<unsigned>(r);

                    const u64 tdigit =
                        (rr < td.size)
                            ? td.digit[rr]
                            : 0;

                    const u64 qdigit =
                        (rr < qd.size)
                            ? qd.digit[rr]
                            : 0;

                    if (tdigit > qdigit) {
                        violation = true;
                        highest = rr;
                        break;
                    }
                }

                j = 0;
                rank = 0;

                u64 place = 1;
                u64 weight = 1;

                for (unsigned r = 0;
                     r < R;
                     ++r) {

                    const u64 tdigit =
                        (r < td.size)
                            ? td.digit[r]
                            : 0;

                    const u64 qdigit =
                        (r < qd.size)
                            ? qd.digit[r]
                            : 0;

                    const u64 jd =
                        !violation
                            ? tdigit
                            : (r > highest
                                   ? tdigit
                                   : qdigit);

                    j += jd * place;
                    rank += jd * weight;

                    place *= p;
                    weight *= qdigit + 1;
                }

                h = violation
                    ? highest
                    : R;

                (void)h;

                const u64 brute =
                    brute_largest_admissible(
                        t,
                        q,
                        p);

                if (j != brute) {
                    ++failures;

                    if (failures <= 5) {
                        std::cout
                            << "FAIL REFERENCE"
                            << " p=" << p
                            << " q=" << q
                            << " t=" << t
                            << " direct_j=" << j
                            << " brute_j="
                            << brute
                            << '\n';
                    }

                    continue;
                }

                const u64 expected_rank =
                    rank_j(
                        j,
                        q,
                        p);

                if (rank != expected_rank) {
                    ++failures;

                    if (failures <= 5) {
                        std::cout
                            << "FAIL REFERENCE RANK"
                            << " p=" << p
                            << " q=" << q
                            << " t=" << t
                            << " j=" << j
                            << " direct="
                            << rank
                            << " expected="
                            << expected_rank
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

Result deterministic_tests()
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

    u64 failures = 0;

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

        for (u64 seed = 1;
             seed <= 100;
             ++seed) {

            const u64 n =
                (seed * 6364136223846793005ULL +
                 c.q) %
                (sc.m + 1);

            if (!verify_point(
                    n,
                    sc,
                    false)) {
                pass = false;
                break;
            }
        }

        std::cout
            << "CASE"
            << " p=" << sc.p
            << " q=" << sc.q
            << " m=" << sc.m
            << " pass="
            << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            ++failures;
        }
    }

    return {
        failures == 0,
        15,
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
         * Generate random admissible j.
         */
        u64 j = 0;
        u64 place = 1;

        for (unsigned r = 0;
             r < qd.size;
             ++r) {

            const u64 jd =
                rng() %
                (qd.digit[r] + 1);

            j += jd * place;
            place *= p;
        }

        /*
         * Exclude j=q.
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
        << "START EXPERIMENT 241\n";

    bool overall_pass = true;

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

    const Result random =
        random_tests(
            241241241ULL,
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

    const Result hit =
        random_hit_tests(
            941241241ULL,
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
        << "FINISHED EXPERIMENT 241\n";

    return overall_pass ? 0 : 1;
}