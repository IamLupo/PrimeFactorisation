#include <algorithm>
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
};

bool safe_mul(u64 a, u64 b, u64& out) {
    if (a != 0 && b > UINT64_MAX / a) {
        return false;
    }

    out = a * b;
    return true;
}

bool safe_add(u64 a, u64 b, u64& out) {
    if (b > UINT64_MAX - a) {
        return false;
    }

    out = a + b;
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

u64 common_prefix_digits(u64 j, u64 q, u64 p) {
    u64 r = 0;

    while ((j % p) == (q % p)) {
        ++r;

        j /= p;
        q /= p;

        if (j == 0 && q == 0) {
            break;
        }
    }

    return r;
}

bool build_case(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    CaseData& c
) {
    if (p < 2) {
        return false;
    }

    if (b > p - 2) {
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

    u64 mp1;
    if (!safe_add(s0, qpe, mp1)) {
        return false;
    }

    if (mp1 == 0) {
        return false;
    }

    c.p = p;
    c.a = a;
    c.b = b;
    c.z = z;
    c.e = e;
    c.s0 = s0;
    c.q = q;
    c.pe = pe;
    c.m = mp1 - 1;

    return true;
}

bool interval_for_j(
    const CaseData& c,
    u64 j,
    Interval& out
) {
    if (j >= c.q) {
        return false;
    }

    if (!leq_p(j, c.q, c.p)) {
        return false;
    }

    u64 jpe;
    if (!safe_mul(j, c.pe, jpe)) {
        return false;
    }

    u64 left;
    if (!safe_add(c.s0, jpe, left)) {
        return false;
    }

    const u64 r =
        common_prefix_digits(j, c.q, c.p);

    u64 pr;
    if (!pow_u64(c.p, r, pr)) {
        return false;
    }

    u64 shifted_j = j;

    for (u64 i = 0; i < r; ++i) {
        shifted_j /= c.p;
    }

    u64 prefix_plus_one;
    if (!safe_add(shifted_j, 1, prefix_plus_one)) {
        return false;
    }

    u64 power_er;
    if (!safe_mul(c.pe, pr, power_er)) {
        return false;
    }

    u64 right_plus_one;
    if (!safe_mul(prefix_plus_one, power_er, right_plus_one)) {
        return false;
    }

    if (right_plus_one == 0) {
        return false;
    }

    const u64 right = right_plus_one - 1;

    if (right < left || right > c.m) {
        return false;
    }

    out = {left, right};
    return true;
}

std::vector<Interval> reference_intervals(const CaseData& c) {
    std::vector<Interval> result;

    for (u64 j = 0; j < c.q; ++j) {
        if (!leq_p(j, c.q, c.p)) {
            continue;
        }

        Interval interval;

        if (!interval_for_j(c, j, interval)) {
            continue;
        }

        result.push_back(interval);
    }

    std::sort(
        result.begin(),
        result.end(),
        [](const Interval& x, const Interval& y) {
            if (x.left != y.left) {
                return x.left < y.left;
            }

            return x.right < y.right;
        }
    );

    return result;
}

bool geometry_next_interval(
    const CaseData& c,
    const Interval& current,
    Interval& next
) {
    /*
        Geometry recurrence:

            x_next = E_current + s0 + 1

        where E_current = current.right.

        Since x = s0 + j*p^e,

            j_next = (x_next - s0) / p^e.

        The theorem predicts that this reconstructed j is the
        next admissible Lucas index.
    */

    u64 next_left;

    if (!safe_add(current.right, c.s0, next_left)) {
        return false;
    }

    if (!safe_add(next_left, 1, next_left)) {
        return false;
    }

    if (next_left > c.m) {
        return false;
    }

    if (next_left < c.s0) {
        return false;
    }

    const u64 delta = next_left - c.s0;

    if (c.pe == 0 || delta % c.pe != 0) {
        return false;
    }

    const u64 next_j = delta / c.pe;

    if (next_j >= c.q) {
        return false;
    }

    if (!leq_p(next_j, c.q, c.p)) {
        return false;
    }

    return interval_for_j(c, next_j, next);
}

std::vector<Interval> geometry_intervals(
    const CaseData& c,
    bool& construction_ok
) {
    construction_ok = true;

    std::vector<Interval> result;

    if (c.q == 0) {
        construction_ok = false;
        return result;
    }

    /*
        The first interval always starts at s0.
    */
    Interval current;

    if (!interval_for_j(c, 0, current)) {
        construction_ok = false;
        return result;
    }

    result.push_back(current);

    while (true) {
        Interval next;

        if (!geometry_next_interval(c, current, next)) {
            break;
        }

        if (next.left <= current.left) {
            construction_ok = false;
            return result;
        }

        result.push_back(next);
        current = next;
    }

    return result;
}

bool intervals_equal(
    const std::vector<Interval>& a,
    const std::vector<Interval>& b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].left != b[i].left ||
            a[i].right != b[i].right) {
            return false;
        }
    }

    return true;
}

bool check_case(
    const CaseData& c,
    bool print_details
) {
    const auto reference = reference_intervals(c);

    bool geometry_ok = false;
    const auto geometry =
        geometry_intervals(c, geometry_ok);

    const bool count_ok =
        reference.size() == geometry.size();

    const bool sequence_ok =
        intervals_equal(reference, geometry);

    /*
        Verify the recurrence explicitly for every neighboring
        interval in the geometry-generated sequence.
    */
    bool recurrence_ok = true;
    u64 recurrence_checks = 0;
    u64 recurrence_failures = 0;

    for (size_t i = 1; i < geometry.size(); ++i) {
        ++recurrence_checks;

        u64 expected_left;

        if (!safe_add(
                geometry[i - 1].right,
                c.s0,
                expected_left
            ) ||
            !safe_add(
                expected_left,
                1,
                expected_left
            )) {
            recurrence_ok = false;
            ++recurrence_failures;
            continue;
        }

        if (geometry[i].left != expected_left) {
            recurrence_ok = false;
            ++recurrence_failures;
        }
    }

    /*
        Verify the first and last geometric boundaries.
    */
    bool first_ok =
        !geometry.empty() &&
        geometry.front().left == c.s0;

    bool last_ok =
        !geometry.empty() &&
        geometry.back().right == c.m - c.s0;

    /*
        The number of intervals should be

            I = product(q_i + 1) - 1.

        Compute this directly from q's base-p digits.
    */
    u64 digit_product = 1;
    u64 qq = c.q;

    bool digit_product_ok = true;

    while (qq != 0) {
        const u64 digit = qq % c.p;

        u64 factor;
        if (!safe_add(digit, 1, factor) ||
            !safe_mul(digit_product, factor, digit_product)) {
            digit_product_ok = false;
            break;
        }

        qq /= c.p;
    }

    if (c.q == 0) {
        digit_product_ok = false;
    }

    u64 expected_intervals = 0;

    if (digit_product_ok) {
        expected_intervals = digit_product - 1;
    }

    const bool count_formula_ok =
        digit_product_ok &&
        reference.size() == expected_intervals;

    const bool final_pass =
        geometry_ok &&
        count_ok &&
        sequence_ok &&
        recurrence_ok &&
        first_ok &&
        last_ok &&
        count_formula_ok;

    if (print_details) {
        std::cout << "CASE\n";

        std::cout << "p=" << c.p
                  << " a=" << c.a
                  << " b=" << c.b
                  << " z=" << c.z
                  << " e=" << c.e
                  << " s0=" << c.s0
                  << " q=" << c.q << '\n';

        std::cout << "m=" << c.m << '\n';

        std::cout << "reference_intervals="
                  << reference.size() << '\n';

        std::cout << "geometry_intervals="
                  << geometry.size() << '\n';

        std::cout << "geometry_construction_pass="
                  << geometry_ok << '\n';

        std::cout << "count_pass="
                  << count_ok << '\n';

        std::cout << "exact_sequence_pass="
                  << sequence_ok << '\n';

        std::cout << "recurrence_checks="
                  << recurrence_checks << '\n';

        std::cout << "recurrence_failures="
                  << recurrence_failures << '\n';

        std::cout << "recurrence_pass="
                  << recurrence_ok << '\n';

        std::cout << "first_boundary_pass="
                  << first_ok << '\n';

        std::cout << "last_boundary_pass="
                  << last_ok << '\n';

        std::cout << "expected_intervals="
                  << expected_intervals << '\n';

        std::cout << "count_formula_pass="
                  << count_formula_ok << '\n';

        std::cout << "final_pass="
                  << final_pass << '\n';

        const size_t sample_count =
            std::min<size_t>(10, geometry.size());

        for (size_t i = 0; i < sample_count; ++i) {
            std::cout << "sample[" << i << "]="
                      << "[" << geometry[i].left
                      << "," << geometry[i].right
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

bool run_random_suite(
    u64 p,
    int cases,
    u64 max_a,
    u64 max_z,
    u64 max_q,
    std::mt19937_64& rng
) {
    int failures = 0;
    int skipped = 0;

    for (int i = 0; i < cases; ++i) {
        const u64 a =
            random_bounded(rng, 0, max_a);

        const u64 b =
            random_bounded(rng, 0, p - 2);

        const u64 z =
            random_bounded(rng, 0, max_z);

        const u64 q =
            random_bounded(rng, 1, max_q);

        CaseData c;

        if (!build_case(
                p,
                a,
                b,
                z,
                q,
                c
            )) {
            ++skipped;
            continue;
        }

        /*
            Keep the reference construction practical.
        */
        if (q > 2000000) {
            ++skipped;
            continue;
        }

        if (!check_case(c, false)) {
            ++failures;

            if (failures <= 3) {
                std::cout << "RANDOM FAILURE\n";

                std::cout << "p=" << c.p
                          << " a=" << c.a
                          << " b=" << c.b
                          << " z=" << c.z
                          << " e=" << c.e
                          << " s0=" << c.s0
                          << " q=" << c.q
                          << " m=" << c.m
                          << "\n\n";
            }
        }
    }

    std::cout << "RANDOM SUITE\n";

    std::cout << "p=" << p
              << " cases=" << cases
              << " failures=" << failures
              << " skipped=" << skipped
              << '\n';

    std::cout << "suite_pass="
              << (failures == 0)
              << "\n\n";

    return failures == 0;
}

int main() {
    std::cout << "START EXPERIMENT 217\n\n";

    bool all_pass = true;

    /*
        Small dense binary.
    */
    {
        CaseData c;

        build_case(
            2,
            0,
            0,
            2,
            63,
            c
        );

        all_pass =
            check_case(c, true) &&
            all_pass;
    }

    /*
        Small ternary.
    */
    {
        CaseData c;

        build_case(
            3,
            0,
            0,
            3,
            80,
            c
        );

        all_pass =
            check_case(c, true) &&
            all_pass;
    }

    /*
        Small base-5.
    */
    {
        CaseData c;

        build_case(
            5,
            0,
            0,
            2,
            100,
            c
        );

        all_pass =
            check_case(c, true) &&
            all_pass;
    }

    /*
        Nontrivial s0.
    */
    {
        CaseData c;

        build_case(
            3,
            3,
            2,
            1,
            80,
            c
        );

        all_pass =
            check_case(c, true) &&
            all_pass;
    }

    /*
        Large dense binary.
    */
    {
        CaseData c;

        build_case(
            2,
            0,
            0,
            39,
            1048575,
            c
        );

        all_pass =
            check_case(c, false) &&
            all_pass;
    }

    /*
        Large base-5 case that stays inside uint64_t.
    */
    {
        CaseData c;

        build_case(
            5,
            0,
            0,
            14,
            390624,
            c
        );

        all_pass =
            check_case(c, false) &&
            all_pass;
    }

    /*
        Sparse binary case.
    */
    {
        CaseData c;

        build_case(
            2,
            0,
            0,
            19,
            1073741825ULL,
            c
        );

        all_pass =
            check_case(c, true) &&
            all_pass;
    }

    std::mt19937_64 rng(217217217ULL);

    all_pass =
        run_random_suite(
            2,
            300,
            10,
            10,
            1000000,
            rng
        ) &&
        all_pass;

    all_pass =
        run_random_suite(
            3,
            300,
            8,
            8,
            500000,
            rng
        ) &&
        all_pass;

    all_pass =
        run_random_suite(
            5,
            300,
            6,
            6,
            500000,
            rng
        ) &&
        all_pass;

    std::cout << "THEOREM TARGET\n";

    std::cout
        << "x_0 = s0\n";

    std::cout
        << "x_{k+1} = E_k + s0 + 1\n";

    std::cout
        << "j_k = (x_k-s0)/p^e\n";

    std::cout
        << "geometry-generated intervals == Lucas-generated intervals\n";

    std::cout
        << "no sorting is used by the geometry generator\n";

    std::cout << '\n';

    std::cout << "OVERALL_PASS="
              << all_pass
              << '\n';

    std::cout << "FINISHED EXPERIMENT 217\n";

    return all_pass ? 0 : 1;
}
