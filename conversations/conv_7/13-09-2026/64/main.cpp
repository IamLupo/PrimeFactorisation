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

u64 digit_at(u64 x, u64 p, u64 i) {
    while (i > 0) {
        x /= p;
        --i;
    }

    return x % p;
}

u64 set_digit(
    u64 x,
    u64 p,
    u64 pos,
    u64 value
) {
    u64 power;
    pow_u64(p, pos, power);

    const u64 old_digit =
        digit_at(x, p, pos);

    const u64 old_value =
        old_digit * power;

    const u64 new_value =
        value * power;

    return x - old_value + new_value;
}

/*
    Reference successor from the previously established formula.
*/
bool successor_formula(
    u64 j,
    const CaseData& c,
    u64& next,
    u64& r
) {
    if (!first_lower_digit(
            j,
            c.q,
            c.p,
            r
        )) {
        return false;
    }

    u64 pr;

    if (!pow_u64(c.p, r, pr)) {
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

    if (!safe_mul(
            prefix_plus_one,
            pr,
            next
        )) {
        return false;
    }

    return next > j && next < c.q;
}

/*
    Independently construct the successor digit-by-digit:

      digits below r -> 0
      digit r        -> digit r + 1
      digits above r -> unchanged
*/
bool successor_digits(
    u64 j,
    const CaseData& c,
    u64& next,
    u64& r
) {
    if (!first_lower_digit(
            j,
            c.q,
            c.p,
            r
        )) {
        return false;
    }

    next = j;

    const u64 current_r_digit =
        digit_at(j, c.p, r);

    /*
        The digit must be strictly smaller than q_r.
        Therefore incrementing it remains <= q_r.
    */
    const u64 q_r =
        digit_at(c.q, c.p, r);

    if (current_r_digit >= q_r) {
        return false;
    }

    const u64 incremented =
        current_r_digit + 1;

    /*
        Reset all lower digits.
    */
    for (u64 i = 0; i < r; ++i) {
        next = set_digit(
            next,
            c.p,
            i,
            0
        );
    }

    /*
        Increment the first lower digit.
    */
    next = set_digit(
        next,
        c.p,
        r,
        incremented
    );

    return next > j && next < c.q;
}

std::vector<u64> reference_indices(
    const CaseData& c
) {
    std::vector<u64> result;

    for (u64 j = 0; j < c.q; ++j) {
        if (leq_p(j, c.q, c.p)) {
            result.push_back(j);
        }
    }

    return result;
}

std::vector<u64> successor_sequence(
    const CaseData& c,
    bool digit_mode,
    bool& construction_ok,
    std::vector<u64>& r_values
) {
    construction_ok = true;
    r_values.clear();

    std::vector<u64> result;

    u64 current = 0;

    if (!leq_p(
            current,
            c.q,
            c.p
        )) {
        construction_ok = false;
        return result;
    }

    result.push_back(current);

    while (true) {
        u64 next;
        u64 r;

        const bool ok =
            digit_mode
            ? successor_digits(
                  current,
                  c,
                  next,
                  r
              )
            : successor_formula(
                  current,
                  c,
                  next,
                  r
              );

        if (!ok) {
            break;
        }

        if (next <= current) {
            construction_ok = false;
            break;
        }

        result.push_back(next);
        r_values.push_back(r);

        current = next;
    }

    return result;
}

bool equal_sequences(
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

/*
    Verify the digit transition independently for every edge.
*/
bool verify_digit_transition(
    const std::vector<u64>& sequence,
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    for (size_t k = 0; k + 1 < sequence.size(); ++k) {
        ++checks;

        const u64 j = sequence[k];
        const u64 next = sequence[k + 1];

        u64 r;

        if (!first_lower_digit(
                j,
                c.q,
                c.p,
                r
            )) {
            ++failures;
            continue;
        }

        /*
            Check every digit up to the useful range.
        */
        u64 temp_q = c.q;
        u64 digit_count = 0;

        while (temp_q != 0) {
            ++digit_count;
            temp_q /= c.p;
        }

        for (u64 i = 0; i < digit_count; ++i) {
            const u64 jd =
                digit_at(j, c.p, i);

            const u64 nd =
                digit_at(next, c.p, i);

            const u64 expected =
                (i < r)
                    ? 0
                    : (i == r)
                        ? jd + 1
                        : jd;

            if (nd != expected) {
                ++failures;
            }

            const u64 qd =
                digit_at(c.q, c.p, i);

            if (nd > qd) {
                ++failures;
            }
        }
    }

    return failures == 0;
}

bool check_case(
    const CaseData& c,
    bool print_details
) {
    const auto reference =
        reference_indices(c);

    bool formula_construction_ok;
    std::vector<u64> formula_r;

    const auto formula_sequence =
        successor_sequence(
            c,
            false,
            formula_construction_ok,
            formula_r
        );

    bool digit_construction_ok;
    std::vector<u64> digit_r;

    const auto digit_sequence =
        successor_sequence(
            c,
            true,
            digit_construction_ok,
            digit_r
        );

    const bool formula_exact =
        equal_sequences(
            reference,
            formula_sequence
        );

    const bool digit_exact =
        equal_sequences(
            reference,
            digit_sequence
        );

    const bool formula_vs_digit =
        equal_sequences(
            formula_sequence,
            digit_sequence
        );

    u64 digit_checks;
    u64 digit_failures;

    const bool transition_ok =
        verify_digit_transition(
            digit_sequence,
            c,
            digit_checks,
            digit_failures
        );

    bool r_match = true;

    if (formula_r.size() != digit_r.size()) {
        r_match = false;
    } else {
        for (size_t i = 0; i < formula_r.size(); ++i) {
            if (formula_r[i] != digit_r[i]) {
                r_match = false;
                break;
            }
        }
    }

    /*
        Verify the number of possible digit-tuples:

            product_i(q_i+1)

        includes q itself, so admissible j<q count is
            product_i(q_i+1)-1.
    */
    u64 product = 1;
    u64 qtemp = c.q;
    bool product_ok = true;

    while (qtemp != 0) {
        const u64 qdigit =
            qtemp % c.p;

        u64 factor;

        if (!safe_add(
                qdigit,
                1,
                factor
            ) ||
            !safe_mul(
                product,
                factor,
                product
            )) {
            product_ok = false;
            break;
        }

        qtemp /= c.p;
    }

    u64 expected_count = 0;

    if (product_ok) {
        expected_count = product - 1;
    }

    const bool count_ok =
        product_ok &&
        reference.size() == expected_count;

    const bool first_zero_ok =
        !digit_sequence.empty() &&
        digit_sequence.front() == 0;

    const bool final_pass =
        formula_construction_ok &&
        digit_construction_ok &&
        formula_exact &&
        digit_exact &&
        formula_vs_digit &&
        transition_ok &&
        r_match &&
        count_ok &&
        first_zero_ok;

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
            << "formula_indices="
            << formula_sequence.size()
            << '\n';

        std::cout
            << "digit_indices="
            << digit_sequence.size()
            << '\n';

        std::cout
            << "formula_construction_pass="
            << formula_construction_ok
            << '\n';

        std::cout
            << "digit_construction_pass="
            << digit_construction_ok
            << '\n';

        std::cout
            << "reference_vs_formula_pass="
            << formula_exact
            << '\n';

        std::cout
            << "reference_vs_digits_pass="
            << digit_exact
            << '\n';

        std::cout
            << "formula_vs_digits_pass="
            << formula_vs_digit
            << '\n';

        std::cout
            << "digit_transition_checks="
            << digit_checks
            << '\n';

        std::cout
            << "digit_transition_failures="
            << digit_failures
            << '\n';

        std::cout
            << "digit_transition_pass="
            << transition_ok
            << '\n';

        std::cout
            << "r_sequence_pass="
            << r_match
            << '\n';

        std::cout
            << "expected_count="
            << expected_count
            << '\n';

        std::cout
            << "count_pass="
            << count_ok
            << '\n';

        std::cout
            << "first_index_zero_pass="
            << first_zero_ok
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << '\n';

        const size_t sample_count =
            std::min<size_t>(
                15,
                digit_sequence.size()
            );

        for (size_t i = 0; i < sample_count; ++i) {
            std::cout
                << "sample_j["
                << i
                << "]="
                << digit_sequence[i];

            if (i < digit_r.size()) {
                const u64 r =
                    digit_r[i];

                std::cout
                    << " r="
                    << r;

                std::cout
                    << " digits_before=";

                u64 temp =
                    digit_sequence[i];

                for (u64 d = 0; d <= r; ++d) {
                    std::cout
                        << digit_at(
                               temp,
                               c.p,
                               d
                           );

                    if (d != r) {
                        std::cout << ",";
                    }
                }

                std::cout
                    << " digits_after=";

                temp =
                    digit_sequence[i + 1];

                for (u64 d = 0; d <= r; ++d) {
                    std::cout
                        << digit_at(
                               temp,
                               c.p,
                               d
                           );

                    if (d != r) {
                        std::cout << ",";
                    }
                }
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
        << "START EXPERIMENT 219\n\n";

    bool all_pass = true;

    /*
        Binary.
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
        Ternary.
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
        Base 5.
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
        Large base-5.
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
        219219219ULL
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
        << "For j -> j' with r=r(j):\n";

    std::cout
        << "j'_i = 0 for i<r\n";

    std::cout
        << "j'_r = j_r + 1\n";

    std::cout
        << "j'_i = j_i for i>r\n";

    std::cout
        << "j_r < q_r before increment\n";

    std::cout
        << "j'_i <= q_i for every digit\n";

    std::cout
        << "digit successor == floor/power successor\n";

    std::cout
        << "digit successor == independently enumerated Lucas set\n";

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << all_pass
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 219\n";

    return all_pass ? 0 : 1;
}
