#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct Interval {
    u64 left;
    u64 right;
};

struct CaseData {
    u64 p;
    u64 a;
    u64 b;
    u64 z;
    u64 e;
    u64 s0;
    u64 q;
    u64 pe;
    u64 m;

    std::vector<u64> q_digits;
    std::vector<u64> radix;
    std::vector<u64> weights;
    std::vector<u64> powers;

    u64 tuple_count;
    u64 interval_count;
};

bool safe_add(u64 a, u64 b, u64& out) {
    if (b > UINT64_MAX - a) {
        return false;
    }

    out = a + b;
    return true;
}

bool safe_mul(u64 a, u64 b, u64& out) {
    if (a != 0 && b > UINT64_MAX / a) {
        return false;
    }

    out = a * b;
    return true;
}

bool pow_u64(u64 base, u64 exp, u64& out) {
    out = 1;

    for (u64 i = 0; i < exp; ++i) {
        if (!safe_mul(out, base, out)) {
            return false;
        }
    }

    return true;
}

u64 p_digit(u64 x, u64 p, u64 i) {
    while (i > 0) {
        x /= p;
        --i;
    }

    return x % p;
}

bool leq_p(u64 j, u64 q, u64 p) {
    while (j != 0 || q != 0) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd > qd) {
            return false;
        }

        j /= p;
        q /= p;
    }

    return true;
}

bool build_case(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    CaseData& c
) {
    if (p < 2 || q == 0 || b > p - 2) {
        return false;
    }

    u64 b1;
    if (!safe_add(b, 1, b1)) {
        return false;
    }

    u64 pa;
    if (!pow_u64(p, a, pa)) {
        return false;
    }

    u64 s0;
    if (!safe_mul(b1, pa, s0)) {
        return false;
    }

    u64 e;
    if (!safe_add(a, z, e)) {
        return false;
    }

    if (!safe_add(e, 1, e)) {
        return false;
    }

    u64 pe;
    if (!pow_u64(p, e, pe)) {
        return false;
    }

    u64 qpe;
    if (!safe_mul(q, pe, qpe)) {
        return false;
    }

    u64 m_plus_one;
    if (!safe_add(s0, qpe, m_plus_one)) {
        return false;
    }

    if (m_plus_one == 0) {
        return false;
    }

    c = {};
    c.p = p;
    c.a = a;
    c.b = b;
    c.z = z;
    c.e = e;
    c.s0 = s0;
    c.q = q;
    c.pe = pe;
    c.m = m_plus_one - 1;

    u64 temp = q;

    while (true) {
        const u64 d = temp % p;

        c.q_digits.push_back(d);
        c.radix.push_back(d + 1);

        temp /= p;

        if (temp == 0) {
            break;
        }
    }

    c.weights.resize(c.radix.size());
    c.powers.resize(c.radix.size());

    c.weights[0] = 1;
    c.powers[0] = 1;

    for (size_t i = 1; i < c.radix.size(); ++i) {
        if (!safe_mul(
                c.weights[i - 1],
                c.radix[i - 1],
                c.weights[i]
            )) {
            return false;
        }

        if (!safe_mul(
                c.powers[i - 1],
                p,
                c.powers[i]
            )) {
            return false;
        }
    }

    c.tuple_count = 1;

    for (u64 r : c.radix) {
        if (!safe_mul(
                c.tuple_count,
                r,
                c.tuple_count
            )) {
            return false;
        }
    }

    c.interval_count = c.tuple_count - 1;

    return true;
}

/*
    Mixed-radix unrank:

        j_i = floor(k / W_i) mod (q_i+1)

    and

        j = sum_i j_i p^i.
*/
bool unrank_j(
    u64 k,
    const CaseData& c,
    u64& j
) {
    if (k >= c.tuple_count) {
        return false;
    }

    j = 0;

    for (size_t i = 0; i < c.radix.size(); ++i) {
        const u64 d =
            (k / c.weights[i]) % c.radix[i];

        u64 term;

        if (!safe_mul(
                d,
                c.powers[i],
                term
            )) {
            return false;
        }

        if (!safe_add(
                j,
                term,
                j
            )) {
            return false;
        }
    }

    return true;
}

/*
    Rank:

        R(j) = sum_i j_i W_i.
*/
bool rank_j(
    u64 j,
    const CaseData& c,
    u64& rank
) {
    rank = 0;

    for (size_t i = 0; i < c.radix.size(); ++i) {
        const u64 d =
            p_digit(j, c.p, static_cast<u64>(i));

        if (d > c.q_digits[i]) {
            return false;
        }

        u64 term;

        if (!safe_mul(
                d,
                c.weights[i],
                term
            )) {
            return false;
        }

        if (!safe_add(
                rank,
                term,
                rank
            )) {
            return false;
        }
    }

    return true;
}

/*
    r(j) = first digit position with j_r < q_r.
*/
bool first_lower_digit(
    u64 j,
    u64 q,
    u64 p,
    u64& r
) {
    r = 0;

    while (true) {
        const u64 jd = j % p;
        const u64 qd = q % p;

        if (jd < qd) {
            return true;
        }

        if (jd > qd) {
            return false;
        }

        j /= p;
        q /= p;
        ++r;

        if (j == 0 && q == 0) {
            return false;
        }
    }
}

/*
    Reference interval from an already known Lucas index j.
*/
bool interval_from_j(
    u64 j,
    const CaseData& c,
    Interval& out,
    u64& r_out
) {
    if (j >= c.q) {
        return false;
    }

    if (!leq_p(j, c.q, c.p)) {
        return false;
    }

    u64 r;

    if (!first_lower_digit(
            j,
            c.q,
            c.p,
            r
        )) {
        return false;
    }

    r_out = r;

    u64 jpe;

    if (!safe_mul(
            j,
            c.pe,
            jpe
        )) {
        return false;
    }

    u64 left;

    if (!safe_add(
            c.s0,
            jpe,
            left
        )) {
        return false;
    }

    u64 pr;

    if (!pow_u64(
            c.p,
            r,
            pr
        )) {
        return false;
    }

    const u64 prefix = j / pr;

    u64 prefix_plus_one;

    if (!safe_add(
            prefix,
            1,
            prefix_plus_one
        )) {
        return false;
    }

    u64 power_er;

    if (!safe_mul(
            c.pe,
            pr,
            power_er
        )) {
        return false;
    }

    u64 right_plus_one;

    if (!safe_mul(
            prefix_plus_one,
            power_er,
            right_plus_one
        )) {
        return false;
    }

    if (right_plus_one == 0) {
        return false;
    }

    const u64 right =
        right_plus_one - 1;

    if (right < left || right > c.m) {
        return false;
    }

    out = {left, right};

    return true;
}

/*
    Direct rank -> interval construction.

    The rank k produces j(k) directly.
    Then r(k) is found directly from j(k).
*/
bool interval_from_rank(
    u64 k,
    const CaseData& c,
    Interval& out,
    u64& j_out,
    u64& r_out
) {
    if (k >= c.interval_count) {
        return false;
    }

    u64 j;

    if (!unrank_j(
            k,
            c,
            j
        )) {
        return false;
    }

    u64 r;

    if (!first_lower_digit(
            j,
            c.q,
            c.p,
            r
        )) {
        return false;
    }

    u64 jpe;

    if (!safe_mul(
            j,
            c.pe,
            jpe
        )) {
        return false;
    }

    u64 left;

    if (!safe_add(
            c.s0,
            jpe,
            left
        )) {
        return false;
    }

    u64 pr;

    if (!pow_u64(
            c.p,
            r,
            pr
        )) {
        return false;
    }

    const u64 prefix =
        j / pr;

    u64 prefix_plus_one;

    if (!safe_add(
            prefix,
            1,
            prefix_plus_one
        )) {
        return false;
    }

    u64 power_er;

    if (!safe_mul(
            c.pe,
            pr,
            power_er
        )) {
        return false;
    }

    u64 right_plus_one;

    if (!safe_mul(
            prefix_plus_one,
            power_er,
            right_plus_one
        )) {
        return false;
    }

    if (right_plus_one == 0) {
        return false;
    }

    const u64 right =
        right_plus_one - 1;

    out = {left, right};

    j_out = j;
    r_out = r;

    return right >= left &&
           right <= c.m;
}

/*
    Independent explicit construction of the complete reference
    sequence for manageable cases.
*/
std::vector<Interval> reference_intervals(
    const CaseData& c
) {
    std::vector<Interval> result;

    for (u64 j = 0; j < c.q; ++j) {
        if (!leq_p(
                j,
                c.q,
                c.p
            )) {
            continue;
        }

        Interval iv;
        u64 r;

        if (!interval_from_j(
                j,
                c,
                iv,
                r
            )) {
            result.clear();
            return result;
        }

        result.push_back(iv);
    }

    return result;
}

bool interval_equal(
    const Interval& a,
    const Interval& b
) {
    return a.left == b.left &&
           a.right == b.right;
}

/*
    Exhaustive comparison for manageable q.
*/
bool exhaustive_check(
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    if (c.q > 2000000) {
        return true;
    }

    const auto reference =
        reference_intervals(c);

    if (reference.size() != c.interval_count) {
        ++failures;
        return false;
    }

    for (u64 k = 0;
         k < c.interval_count;
         ++k) {

        ++checks;

        Interval direct;
        u64 j;
        u64 r;

        if (!interval_from_rank(
                k,
                c,
                direct,
                j,
                r
            )) {
            ++failures;
            continue;
        }

        if (!interval_equal(
                direct,
                reference[k]
            )) {
            ++failures;
        }

        u64 recovered_rank;

        if (!rank_j(
                j,
                c,
                recovered_rank
            )) {
            ++failures;
        } else if (recovered_rank != k) {
            ++failures;
        }
    }

    return failures == 0;
}

/*
    Random rank checks. This is the main mechanism for huge q.
*/
bool random_rank_check(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 samples,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    if (c.interval_count == 0) {
        ++failures;
        return false;
    }

    for (u64 n = 0;
         n < samples;
         ++n) {

        ++checks;

        const u64 k =
            rng() % c.interval_count;

        Interval direct;
        u64 j;
        u64 r;

        if (!interval_from_rank(
                k,
                c,
                direct,
                j,
                r
            )) {
            ++failures;
            continue;
        }

        /*
            Check Lucas admissibility.
        */
        if (j >= c.q ||
            !leq_p(
                j,
                c.q,
                c.p
            )) {
            ++failures;
        }

        /*
            Check rank recovery.
        */
        u64 recovered_rank;

        if (!rank_j(
                j,
                c,
                recovered_rank
            )) {
            ++failures;
        } else if (recovered_rank != k) {
            ++failures;
        }

        /*
            Recompute the interval through the reference j-formula
            and compare.
        */
        Interval reference;
        u64 reference_r;

        if (!interval_from_j(
                j,
                c,
                reference,
                reference_r
            )) {
            ++failures;
        } else {
            if (!interval_equal(
                    direct,
                    reference
                )) {
                ++failures;
            }

            if (r != reference_r) {
                ++failures;
            }
        }

        /*
            Check interval boundaries.
        */
        if (direct.left < c.s0 ||
            direct.right > c.m ||
            direct.left > direct.right) {
            ++failures;
        }

        /*
            The interval immediately following this one, when it
            exists, must begin after exactly s0 MISS values.
        */
        if (k + 1 < c.interval_count) {
            Interval next;
            u64 next_j;
            u64 next_r;

            if (!interval_from_rank(
                    k + 1,
                    c,
                    next,
                    next_j,
                    next_r
                )) {
                ++failures;
            } else {
                u64 expected_next_left;

                if (!safe_add(
                        direct.right,
                        c.s0,
                        expected_next_left
                    ) ||
                    !safe_add(
                        expected_next_left,
                        1,
                        expected_next_left
                    )) {
                    ++failures;
                } else if (
                    next.left != expected_next_left
                ) {
                    ++failures;
                }
            }
        }
    }

    return failures == 0;
}

bool run_case(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_samples,
    bool print_details
) {
    u64 exhaustive_checks;
    u64 exhaustive_failures;

    const bool exhaustive_pass =
        exhaustive_check(
            c,
            exhaustive_checks,
            exhaustive_failures
        );

    u64 random_checks;
    u64 random_failures;

    const bool random_pass =
        random_rank_check(
            c,
            rng,
            random_samples,
            random_checks,
            random_failures
        );

    /*
        Explicit boundary ranks.
    */
    bool boundary_pass = true;

    const u64 boundary_samples =
        std::min<u64>(
            10,
            c.interval_count
        );

    for (u64 k = 0;
         k < boundary_samples;
         ++k) {

        Interval iv;
        u64 j;
        u64 r;

        if (!interval_from_rank(
                k,
                c,
                iv,
                j,
                r
            )) {
            boundary_pass = false;
        }
    }

    for (u64 offset = 1;
         offset <= boundary_samples;
         ++offset) {

        const u64 k =
            c.interval_count - offset;

        Interval iv;
        u64 j;
        u64 r;

        if (!interval_from_rank(
                k,
                c,
                iv,
                j,
                r
            )) {
            boundary_pass = false;
        }
    }

    /*
        First interval.
    */
    bool first_pass = false;

    {
        Interval first;
        u64 j;
        u64 r;

        if (interval_from_rank(
                0,
                c,
                first,
                j,
                r
            )) {
            first_pass =
                first.left == c.s0;
        }
    }

    /*
        Last interval.
    */
    bool last_pass = false;

    {
        Interval last;
        u64 j;
        u64 r;

        if (interval_from_rank(
                c.interval_count - 1,
                c,
                last,
                j,
                r
            )) {
            last_pass =
                last.right == c.m - c.s0;
        }
    }

    const bool final_pass =
        exhaustive_pass &&
        random_pass &&
        boundary_pass &&
        first_pass &&
        last_pass;

    if (print_details) {
        std::cout
            << "CASE\n";

        std::cout
            << "p=" << c.p
            << " a=" << c.a
            << " b=" << c.b
            << " z=" << c.z
            << " e=" << c.e
            << " s0=" << c.s0
            << " q=" << c.q
            << '\n';

        std::cout
            << "m=" << c.m
            << '\n';

        std::cout
            << "interval_count="
            << c.interval_count
            << '\n';

        std::cout
            << "exhaustive_checks="
            << exhaustive_checks
            << '\n';

        std::cout
            << "exhaustive_failures="
            << exhaustive_failures
            << '\n';

        std::cout
            << "exhaustive_pass="
            << exhaustive_pass
            << '\n';

        std::cout
            << "random_checks="
            << random_checks
            << '\n';

        std::cout
            << "random_failures="
            << random_failures
            << '\n';

        std::cout
            << "random_pass="
            << random_pass
            << '\n';

        std::cout
            << "boundary_pass="
            << boundary_pass
            << '\n';

        std::cout
            << "first_interval_pass="
            << first_pass
            << '\n';

        std::cout
            << "last_interval_pass="
            << last_pass
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << '\n';

        const u64 sample_count =
            std::min<u64>(
                10,
                c.interval_count
            );

        for (u64 k = 0;
             k < sample_count;
             ++k) {

            Interval iv;
            u64 j;
            u64 r;

            interval_from_rank(
                k,
                c,
                iv,
                j,
                r
            );

            std::cout
                << "sample["
                << k
                << "] rank="
                << k
                << " j="
                << j
                << " r="
                << r
                << " interval=["
                << iv.left
                << ","
                << iv.right
                << "]\n";
        }

        std::cout << '\n';
    }

    return final_pass;
}

u64 random_bounded(
    std::mt19937_64& rng,
    u64 lo,
    u64 hi
) {
    std::uniform_int_distribution<u64> dist(lo, hi);
    return dist(rng);
}

bool make_structural_case(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    CaseData& c
) {
    return build_case(
        p,
        a,
        b,
        z,
        q,
        c
    );
}

int main() {
    std::cout
        << "START EXPERIMENT 222\n\n";

    bool all_pass = true;

    std::mt19937_64 rng(
        222222222ULL
    );

    /*
        Small exhaustive cases.
    */
    {
        CaseData c;

        make_structural_case(
            2,
            0,
            0,
            2,
            63,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                1000,
                true
            ) &&
            all_pass;
    }

    {
        CaseData c;

        make_structural_case(
            3,
            0,
            0,
            3,
            80,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                1000,
                true
            ) &&
            all_pass;
    }

    {
        CaseData c;

        make_structural_case(
            5,
            0,
            0,
            2,
            100,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                1000,
                true
            ) &&
            all_pass;
    }

    /*
        Nontrivial s0.
    */
    {
        CaseData c;

        make_structural_case(
            3,
            3,
            2,
            1,
            80,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                1000,
                true
            ) &&
            all_pass;
    }

    /*
        Larger exhaustive binary.
    */
    {
        CaseData c;

        make_structural_case(
            2,
            0,
            0,
            20,
            1048575,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                10000,
                false
            ) &&
            all_pass;
    }

    /*
        Large cases tested through direct rank -> interval only.
    */
    {
        CaseData c;

        make_structural_case(
            2,
            0,
            0,
            40,
            1048575,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                100000,
                true
            ) &&
            all_pass;
    }

    {
        CaseData c;

        make_structural_case(
            3,
            0,
            0,
            10,
            1000001234567ULL,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                100000,
                true
            ) &&
            all_pass;
    }

    {
        CaseData c;

        make_structural_case(
            5,
            0,
            0,
            8,
            1000007654321ULL,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                100000,
                true
            ) &&
            all_pass;
    }

    /*
        Nontrivial s0 + large q.
    */
    {
        CaseData c;

        make_structural_case(
            3,
            4,
            1,
            6,
            987654321ULL,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                100000,
                true
            ) &&
            all_pass;
    }

    /*
        Sparse q.
    */
    {
        CaseData c;

        make_structural_case(
            2,
            0,
            0,
            19,
            1073741825ULL,
            c
        );

        all_pass =
            run_case(
                c,
                rng,
                100000,
                true
            ) &&
            all_pass;
    }

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "k -> j(k) by mixed-radix unrank\n";

    std::cout
        << "j(k) -> r(k) by first digit with j_r < q_r\n";

    std::cout
        << "x_k = s0 + j(k) p^e\n";

    std::cout
        << "E_k = (floor(j(k)/p^r)+1)p^(e+r)-1\n";

    std::cout
        << "direct rank -> interval equals reference interval\n";

    std::cout
        << "rank k equals the ordered HIT interval index\n";

    std::cout
        << "first interval begins at s0\n";

    std::cout
        << "last interval ends at m-s0\n";

    std::cout
        << "x_{k+1} = E_k + s0 + 1\n";

    std::cout
        << "large cases use no O(q) enumeration\n";

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << all_pass
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 222\n";

    return all_pass ? 0 : 1;
}
