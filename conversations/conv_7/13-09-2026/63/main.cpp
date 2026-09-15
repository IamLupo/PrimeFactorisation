#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct CaseData {
    u64 p;
    u64 a;
    u64 b;
    u64 z;
    u64 e;
    u64 s0;
    u64 q;
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

bool build_case(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    CaseData& c
) {
    if (p < 2 || b > p - 2 || q == 0) {
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

    c.p = p;
    c.a = a;
    c.b = b;
    c.z = z;
    c.e = e;
    c.s0 = s0;
    c.q = q;

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

/*
    r(j) = first digit position where j_r < q_r.

    For every admissible j<q this exists.
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

        /*
            j=q would mean the original index was q,
            which is outside the admissible interval.
        */
        if (j == 0 && q == 0) {
            return false;
        }
    }
}

/*
    Pure p-adic successor:

        next(j)
        =
        (floor(j / p^r) + 1) p^r.
*/
bool successor_padic(
    u64 j,
    const CaseData& c,
    u64& next,
    u64& r_out
) {
    u64 r;

    if (!first_lower_digit(j, c.q, c.p, r)) {
        return false;
    }

    r_out = r;

    u64 pr;
    if (!pow_u64(c.p, r, pr)) {
        return false;
    }

    const u64 prefix = j / pr;

    u64 prefix_plus_one;
    if (!safe_add(prefix, 1, prefix_plus_one)) {
        return false;
    }

    if (!safe_mul(prefix_plus_one, pr, next)) {
        return false;
    }

    if (next <= j || next >= c.q) {
        return false;
    }

    if (!leq_p(next, c.q, c.p)) {
        return false;
    }

    return true;
}

/*
    Independently enumerate every admissible Lucas index.
*/
std::vector<u64> reference_indices(const CaseData& c) {
    std::vector<u64> result;

    for (u64 j = 0; j < c.q; ++j) {
        if (leq_p(j, c.q, c.p)) {
            result.push_back(j);
        }
    }

    return result;
}

/*
    Generate the same sequence exclusively from the
    p-adic successor map.
*/
std::vector<u64> successor_indices(
    const CaseData& c,
    bool& construction_ok,
    std::vector<u64>& r_values
) {
    construction_ok = true;

    std::vector<u64> result;

    if (c.q == 0) {
        construction_ok = false;
        return result;
    }

    u64 current = 0;

    if (!leq_p(current, c.q, c.p)) {
        construction_ok = false;
        return result;
    }

    result.push_back(current);

    while (true) {
        u64 next;
        u64 r;

        if (!successor_padic(
                current,
                c,
                next,
                r
            )) {
            break;
        }

        result.push_back(next);
        r_values.push_back(r);

        current = next;
    }

    return result;
}

/*
    Direct verification of the successor identity for each
    neighboring admissible pair.
*/
bool verify_successor_formula(
    const std::vector<u64>& sequence,
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    for (size_t i = 0; i + 1 < sequence.size(); ++i) {
        ++checks;

        u64 r;

        if (!first_lower_digit(
                sequence[i],
                c.q,
                c.p,
                r
            )) {
            ++failures;
            continue;
        }

        u64 pr;

        if (!pow_u64(c.p, r, pr)) {
            ++failures;
            continue;
        }

        const u64 prefix = sequence[i] / pr;

        u64 prefix_plus_one;

        if (!safe_add(
                prefix,
                1,
                prefix_plus_one
            )) {
            ++failures;
            continue;
        }

        u64 expected;

        if (!safe_mul(
                prefix_plus_one,
                pr,
                expected
            )) {
            ++failures;
            continue;
        }

        if (sequence[i + 1] != expected) {
            ++failures;
        }
    }

    return failures == 0;
}

/*
    Check the r-value claimed by the successor formula against
    the actual first lower digit.
*/
bool verify_r_values(
    const std::vector<u64>& sequence,
    const std::vector<u64>& r_values,
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    if (sequence.size() != r_values.size() + 1) {
        ++failures;
        return false;
    }

    for (size_t i = 0; i < r_values.size(); ++i) {
        ++checks;

        u64 actual_r;

        if (!first_lower_digit(
                sequence[i],
                c.q,
                c.p,
                actual_r
            )) {
            ++failures;
            continue;
        }

        if (actual_r != r_values[i]) {
            ++failures;
        }
    }

    return failures == 0;
}

bool sequences_equal(
    const std::vector<u64>& a,
    const std::vector<u64>& b
) {
    if (a.size() != b.size()) {
        return false;
    }

    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }

    return true;
}

bool check_case(
    const CaseData& c,
    bool print_details
) {
    const auto reference =
        reference_indices(c);

    std::vector<u64> r_values;

    bool construction_ok;

    const auto generated =
        successor_indices(
            c,
            construction_ok,
            r_values
        );

    bool formula_ok;
    u64 formula_checks;
    u64 formula_failures;

    formula_ok =
        verify_successor_formula(
            generated,
            c,
            formula_checks,
            formula_failures
        );

    bool r_ok;
    u64 r_checks;
    u64 r_failures;

    r_ok =
        verify_r_values(
            generated,
            r_values,
            c,
            r_checks,
            r_failures
        );

    const bool count_ok =
        reference.size() == generated.size();

    const bool sequence_ok =
        sequences_equal(
            reference,
            generated
        );

    const bool first_ok =
        !generated.empty() &&
        generated.front() == 0;

    const bool last_ok =
        !generated.empty() &&
        generated.back() < c.q;

    const bool final_pass =
        construction_ok &&
        count_ok &&
        sequence_ok &&
        formula_ok &&
        r_ok &&
        first_ok &&
        last_ok;

    if (print_details) {
        std::cout << "CASE\n";

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
            << "reference_indices="
            << reference.size()
            << '\n';

        std::cout
            << "successor_indices="
            << generated.size()
            << '\n';

        std::cout
            << "construction_pass="
            << construction_ok
            << '\n';

        std::cout
            << "count_pass="
            << count_ok
            << '\n';

        std::cout
            << "exact_sequence_pass="
            << sequence_ok
            << '\n';

        std::cout
            << "successor_formula_checks="
            << formula_checks
            << '\n';

        std::cout
            << "successor_formula_failures="
            << formula_failures
            << '\n';

        std::cout
            << "successor_formula_pass="
            << formula_ok
            << '\n';

        std::cout
            << "r_checks="
            << r_checks
            << '\n';

        std::cout
            << "r_failures="
            << r_failures
            << '\n';

        std::cout
            << "r_tracking_pass="
            << r_ok
            << '\n';

        std::cout
            << "first_index="
            << (generated.empty() ? 0 : generated.front())
            << '\n';

        std::cout
            << "first_index_pass="
            << first_ok
            << '\n';

        std::cout
            << "last_index="
            << (generated.empty() ? 0 : generated.back())
            << '\n';

        std::cout
            << "last_index_pass="
            << last_ok
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << '\n';

        const size_t sample_count =
            std::min<size_t>(
                15,
                generated.size()
            );

        for (size_t i = 0; i < sample_count; ++i) {
            std::cout
                << "sample_j[" << i << "]="
                << generated[i];

            if (i < r_values.size()) {
                std::cout
                    << " r=" << r_values[i];
            }

            std::cout << '\n';
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
            random_bounded(
                rng,
                0,
                max_a
            );

        const u64 b =
            random_bounded(
                rng,
                0,
                p - 2
            );

        const u64 z =
            random_bounded(
                rng,
                0,
                max_z
            );

        const u64 q =
            random_bounded(
                rng,
                1,
                max_q
            );

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
            The reference construction is O(q), so keep q moderate.
        */
        if (q > 2000000) {
            ++skipped;
            continue;
        }

        if (!check_case(c, false)) {
            ++failures;

            if (failures <= 3) {
                std::cout
                    << "RANDOM FAILURE\n";

                std::cout
                    << "p=" << c.p
                    << " a=" << c.a
                    << " b=" << c.b
                    << " z=" << c.z
                    << " e=" << c.e
                    << " s0=" << c.s0
                    << " q=" << c.q
                    << '\n';

                std::cout << '\n';
            }
        }
    }

    std::cout
        << "RANDOM SUITE\n";

    std::cout
        << "p=" << p
        << " cases=" << cases
        << " failures=" << failures
        << " skipped=" << skipped
        << '\n';

    std::cout
        << "suite_pass="
        << (failures == 0)
        << "\n\n";

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 218\n\n";

    bool all_pass = true;

    /*
        Binary dense.
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
        Ternary dense.
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
        Base-5.
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
        Large base-5 case known to fit uint64_t.
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
        Sparse binary.
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

    std::mt19937_64 rng(
        218218218ULL
    );

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

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "J = {j : 0 <=_p j <= q, j < q}\n";

    std::cout
        << "r(j) = first digit position with j_r < q_r\n";

    std::cout
        << "next(j) = (floor(j/p^r(j)) + 1) p^r(j)\n";

    std::cout
        << "successor map enumerates every admissible j exactly once\n";

    std::cout
        << "successor sequence == independently enumerated Lucas sequence\n";

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << all_pass
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 218\n";

    return all_pass ? 0 : 1;
}
