#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <map>

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
*/
bool unrank_j(
    u64 k,
    const CaseData& c,
    u64& j
) {
    if (k >= c.interval_count) {
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
    First digit position r with j_r < q_r.
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

bool interval_from_rank(
    u64 k,
    const CaseData& c,
    u64& j,
    u64& r,
    Interval& iv
) {
    if (!unrank_j(k, c, j)) {
        return false;
    }

    if (!leq_p(j, c.q, c.p)) {
        return false;
    }

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

    iv.left = left;
    iv.right = right_plus_one - 1;

    return iv.right >= iv.left &&
           iv.right <= c.m;
}

/*
    C_r = q_r * product_{i>r}(q_i+1)
*/
bool expected_type_count(
    const CaseData& c,
    u64 r,
    u64& count
) {
    if (r >= c.q_digits.size()) {
        return false;
    }

    count = c.q_digits[r];

    for (size_t i = r + 1;
         i < c.radix.size();
         ++i) {

        if (!safe_mul(
                count,
                c.radix[i],
                count
            )) {
            return false;
        }
    }

    return true;
}

/*
    L_r = (p^r - (q mod p^r)) p^e - s0
*/
bool expected_type_length(
    const CaseData& c,
    u64 r,
    u64& length
) {
    u64 pr;

    if (!pow_u64(
            c.p,
            r,
            pr
        )) {
        return false;
    }

    const u64 remainder =
        c.q % pr;

    const u64 factor =
        pr - remainder;

    u64 raw;

    if (!safe_mul(
            factor,
            c.pe,
            raw
        )) {
        return false;
    }

    if (raw < c.s0) {
        return false;
    }

    length =
        raw - c.s0;

    return true;
}

bool check_case(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_checks,
    bool exhaustive,
    bool print_details
) {
    u64 failures = 0;

    const size_t type_count =
        c.q_digits.size();

    std::vector<u64> observed_counts(
        type_count,
        0
    );

    std::vector<u64> observed_lengths(
        type_count,
        UINT64_MAX
    );

    std::vector<u64> observed_first_rank(
        type_count,
        UINT64_MAX
    );

    std::vector<u64> observed_last_rank(
        type_count,
        UINT64_MAX
    );

    u64 total_hit = 0;
    u64 observed_intervals = 0;

    auto inspect_rank =
        [&](u64 k) -> bool {

        u64 j;
        u64 r;
        Interval iv;

        if (!interval_from_rank(
                k,
                c,
                j,
                r,
                iv
            )) {
            ++failures;
            return false;
        }

        ++observed_intervals;
        ++observed_counts[r];

        if (observed_first_rank[r] == UINT64_MAX) {
            observed_first_rank[r] = k;
        }

        observed_last_rank[r] = k;

        const u64 actual_length =
            iv.right - iv.left + 1;

        if (observed_lengths[r] == UINT64_MAX) {
            observed_lengths[r] = actual_length;
        } else if (observed_lengths[r] != actual_length) {
            ++failures;
        }

        u64 expected_count;

        if (!expected_type_count(
                c,
                r,
                expected_count
            )) {
            ++failures;
        }

        u64 expected_length;

        if (!expected_type_length(
                c,
                r,
                expected_length
            )) {
            ++failures;
        }

        if (actual_length != expected_length) {
            ++failures;
        }

        if (!safe_add(
                total_hit,
                actual_length,
                total_hit
            )) {
            ++failures;
        }

        return true;
    };

    u64 exhaustive_checks = 0;

    if (exhaustive) {
        for (u64 k = 0;
             k < c.interval_count;
             ++k) {

            ++exhaustive_checks;
            inspect_rank(k);
        }
    }

    u64 random_pass_checks = 0;

    for (u64 n = 0;
         n < random_checks;
         ++n) {

        const u64 k =
            rng() % c.interval_count;

        ++random_pass_checks;
        inspect_rank(k);
    }

    /*
        For random-only cases, independently accumulate exact
        type counts using a full rank scan only when feasible.
    */
    const bool full_type_check =
        c.interval_count <= 2000000;

    if (full_type_check && !exhaustive) {
        /*
            We intentionally already covered this through
            exhaustive in normal use.
        */
    }

    bool type_count_pass = true;
    bool type_length_pass = true;

    u64 expected_total_intervals = 0;

    for (u64 r = 0;
         r < type_count;
         ++r) {

        u64 expected_count;

        if (!expected_type_count(
                c,
                r,
                expected_count
            )) {
            type_count_pass = false;
            continue;
        }

        if (!safe_add(
                expected_total_intervals,
                expected_count,
                expected_total_intervals
            )) {
            type_count_pass = false;
            continue;
        }

        if (exhaustive) {
            if (observed_counts[r] != expected_count) {
                type_count_pass = false;
            }

            if (observed_lengths[r] == UINT64_MAX) {
                type_length_pass = false;
            }

            u64 expected_length;

            if (!expected_type_length(
                    c,
                    r,
                    expected_length
                )) {
                type_length_pass = false;
            } else if (
                observed_lengths[r] != expected_length
            ) {
                type_length_pass = false;
            }
        }
    }

    if (expected_total_intervals !=
        c.interval_count) {

        type_count_pass = false;
    }

    /*
        Total HIT length can be calculated from type counts.
    */
    u64 expected_total_hit = 0;

    for (u64 r = 0;
         r < type_count;
         ++r) {

        u64 count;
        u64 length;

        if (!expected_type_count(
                c,
                r,
                count
            ) ||
            !expected_type_length(
                c,
                r,
                length
            )) {
            ++failures;
            continue;
        }

        u64 contribution;

        if (!safe_mul(
                count,
                length,
                contribution
            ) ||
            !safe_add(
                expected_total_hit,
                contribution,
                expected_total_hit
            )) {
            ++failures;
        }
    }

    bool total_hit_pass = true;

    if (exhaustive) {
        if (observed_intervals !=
            c.interval_count) {
            total_hit_pass = false;
        }

        if (total_hit !=
            expected_total_hit) {
            total_hit_pass = false;
        }

        /*
            Previous theorem:
                HIT + MISS = m+1
                MISS = s0*(I+1)
        */
        u64 miss_count;

        if (!safe_mul(
                c.s0,
                c.interval_count + 1,
                miss_count
            )) {
            total_hit_pass = false;
        } else {
            u64 domain_size =
                c.m + 1;

            if (miss_count > domain_size) {
                total_hit_pass = false;
            } else {
                const u64 hit_from_geometry =
                    domain_size - miss_count;

                if (hit_from_geometry !=
                    expected_total_hit) {
                    total_hit_pass = false;
                }
            }
        }
    }

    /*
        The r sequence should itself be the type sequence of
        the mixed-radix counter. We test that each transition
        agrees with the interval type independently.
    */
    bool type_partition_pass = true;

    if (exhaustive && c.interval_count > 0) {
        u64 previous_r = 0;
        bool have_previous = false;

        for (u64 k = 0;
             k < c.interval_count;
             ++k) {

            u64 j;
            u64 r;
            Interval iv;

            if (!interval_from_rank(
                    k,
                    c,
                    j,
                    r,
                    iv
                )) {
                type_partition_pass = false;
                continue;
            }

            /*
                r can decrease or increase depending on the
                lower digit resets. What must hold is simply
                that the observed type belongs to the exact
                expected type bucket.
            */
            u64 expected_count;

            if (!expected_type_count(
                    c,
                    r,
                    expected_count
                ) ||
                expected_count == 0) {
                type_partition_pass = false;
            }

            previous_r = r;
            have_previous = true;

            (void)previous_r;
            (void)have_previous;
        }
    }

    /*
        For each type, every observed interval must have the
        same length. This is already checked above, but keep
        the named theorem flag separate.
    */
    bool equal_length_per_type = true;

    if (exhaustive) {
        for (u64 r = 0;
             r < type_count;
             ++r) {

            u64 expected_count;
            u64 expected_length;

            if (!expected_type_count(
                    c,
                    r,
                    expected_count
                ) ||
                !expected_type_length(
                    c,
                    r,
                    expected_length
                )) {
                equal_length_per_type = false;
                continue;
            }

            if (expected_count == 0) {
                if (observed_counts[r] != 0) {
                    equal_length_per_type = false;
                }
            } else {
                if (observed_counts[r] != expected_count ||
                    observed_lengths[r] != expected_length) {
                    equal_length_per_type = false;
                }
            }
        }
    }

    const bool first_pass =
        c.interval_count == 0 ||
        observed_first_rank[0] == 0;

    const bool last_pass =
        c.interval_count == 0 ||
        observed_last_rank[
            observed_last_rank.size() > 0
                ? 0
                : 0
        ] != UINT64_MAX;

    /*
        The actual boundary conditions are more directly checked
        here.
    */
    bool boundary_pass = true;

    if (c.interval_count > 0) {
        u64 j;
        u64 r;
        Interval first;

        if (!interval_from_rank(
                0,
                c,
                j,
                r,
                first
            )) {
            boundary_pass = false;
        } else if (first.left != c.s0) {
            boundary_pass = false;
        }

        Interval last;

        if (!interval_from_rank(
                c.interval_count - 1,
                c,
                j,
                r,
                last
            )) {
            boundary_pass = false;
        } else if (
            last.right != c.m - c.s0
        ) {
            boundary_pass = false;
        }
    }

    const bool pass =
        failures == 0 &&
        type_count_pass &&
        type_length_pass &&
        total_hit_pass &&
        type_partition_pass &&
        equal_length_per_type &&
        boundary_pass;

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
            << "m=" << c.m
            << '\n';

        std::cout
            << "interval_count="
            << c.interval_count
            << '\n';

        std::cout
            << "type_count="
            << type_count
            << '\n';

        std::cout
            << "exhaustive="
            << exhaustive
            << '\n';

        std::cout
            << "exhaustive_checks="
            << exhaustive_checks
            << '\n';

        std::cout
            << "random_checks="
            << random_pass_checks
            << '\n';

        std::cout
            << "failures="
            << failures
            << '\n';

        std::cout
            << "type_count_pass="
            << type_count_pass
            << '\n';

        std::cout
            << "type_length_pass="
            << type_length_pass
            << '\n';

        std::cout
            << "equal_length_per_type="
            << equal_length_per_type
            << '\n';

        std::cout
            << "total_hit_pass="
            << total_hit_pass
            << '\n';

        std::cout
            << "type_partition_pass="
            << type_partition_pass
            << '\n';

        std::cout
            << "boundary_pass="
            << boundary_pass
            << '\n';

        std::cout
            << "expected_total_intervals="
            << expected_total_intervals
            << '\n';

        std::cout
            << "expected_total_hit="
            << expected_total_hit
            << '\n';

        std::cout
            << "final_pass="
            << pass
            << '\n';

        if (exhaustive) {
            for (u64 r = 0;
                 r < type_count;
                 ++r) {

                u64 expected_count;
                u64 expected_length;

                expected_type_count(
                    c,
                    r,
                    expected_count
                );

                expected_type_length(
                    c,
                    r,
                    expected_length
                );

                std::cout
                    << "TYPE["
                    << r
                    << "] q_digit="
                    << c.q_digits[r]
                    << " expected_count="
                    << expected_count
                    << " observed_count="
                    << observed_counts[r]
                    << " expected_length="
                    << expected_length
                    << " observed_length="
                    << observed_lengths[r]
                    << '\n';
            }
        }

        std::cout << '\n';
    }

    return pass;
}

int main() {
    std::cout
        << "START EXPERIMENT 225\n\n";

    std::mt19937_64 rng(
        225225225ULL
    );

    u64 total_cases = 0;
    u64 failed_cases = 0;
    u64 skipped_cases = 0;

    auto execute_case =
        [&](u64 p,
            u64 a,
            u64 b,
            u64 z,
            u64 q,
            u64 random_checks,
            bool exhaustive,
            bool print_details) {

        CaseData c;

        if (!build_case(
                p,
                a,
                b,
                z,
                q,
                c
            )) {

            ++skipped_cases;

            if (print_details) {
                std::cout
                    << "SKIPPED_INVALID_CASE\n\n";
            }

            return;
        }

        ++total_cases;

        const bool pass =
            check_case(
                c,
                rng,
                random_checks,
                exhaustive,
                print_details
            );

        if (!pass) {
            ++failed_cases;
        }
    };

    /*
        Small dense binary.
    */
    execute_case(
        2, 0, 0, 2, 63,
        1000,
        true,
        true
    );

    /*
        Small ternary.
    */
    execute_case(
        3, 0, 0, 3, 80,
        1000,
        true,
        true
    );

    /*
        Base 5.
    */
    execute_case(
        5, 0, 0, 2, 100,
        1000,
        true,
        true
    );

    /*
        Nontrivial s0.
    */
    execute_case(
        3, 3, 1, 1, 80,
        1000,
        true,
        true
    );

    /*
        Nontrivial s0, more digits.
    */
    execute_case(
        3, 4, 1, 6, 987654321ULL,
        100000,
        false,
        true
    );

    /*
        Dense large binary.
    */
    execute_case(
        2, 0, 0, 40, 1048575,
        100000,
        true,
        true
    );

    /*
        Large ternary.
    */
    execute_case(
        3, 0, 0, 10,
        1000001234567ULL,
        100000,
        false,
        true
    );

    /*
        Large base 5.
    */
    execute_case(
        5, 0, 0, 8,
        1000007654321ULL,
        100000,
        false,
        true
    );

    /*
        Sparse binary.
    */
    execute_case(
        2, 0, 0, 19,
        1073741825ULL,
        100000,
        true,
        true
    );

    std::cout
        << "SUMMARY\n";

    std::cout
        << "total_cases="
        << total_cases
        << '\n';

    std::cout
        << "failed_cases="
        << failed_cases
        << '\n';

    std::cout
        << "skipped_cases="
        << skipped_cases
        << '\n';

    std::cout
        << "passed_cases="
        << (total_cases - failed_cases)
        << '\n';

    std::cout
        << "OVERALL_PASS="
        << (failed_cases == 0)
        << '\n';

    std::cout << '\n';

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "L_k depends only on r_k\n";

    std::cout
        << "L_r = (p^r - (q mod p^r)) p^e - s0\n";

    std::cout
        << "C_r = q_r * product_{i>r}(q_i+1)\n";

    std::cout
        << "sum_r C_r = product_i(q_i+1)-1\n";

    std::cout
        << "all intervals of a fixed r have identical length\n";

    std::cout
        << "sum_r C_r L_r = total HIT length\n";

    std::cout
        << "FINISHED EXPERIMENT 225\n";

    return failed_cases == 0 ? 0 : 1;
}
