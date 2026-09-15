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

    if (!first_lower_digit(j, c.q, c.p, r)) {
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

    u64 pr;

    if (!pow_u64(c.p, r, pr)) {
        return false;
    }

    const u64 prefix = j / pr;

    u64 prefix_plus_one;

    if (!safe_add(prefix, 1, prefix_plus_one)) {
        return false;
    }

    u64 power_er;

    if (!safe_mul(c.pe, pr, power_er)) {
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

bool expected_type_length(
    const CaseData& c,
    u64 r,
    u64& length
) {
    u64 pr;

    if (!pow_u64(c.p, r, pr)) {
        return false;
    }

    const u64 remainder = c.q % pr;
    const u64 factor = pr - remainder;

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

    length = raw - c.s0;

    return true;
}

bool check_exhaustive_types(
    const CaseData& c,
    u64& failures,
    std::vector<u64>& observed_counts,
    std::vector<u64>& observed_lengths,
    u64& observed_total_hit
) {
    failures = 0;
    observed_total_hit = 0;

    const size_t T = c.q_digits.size();

    observed_counts.assign(T, 0);
    observed_lengths.assign(T, UINT64_MAX);

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
            ++failures;
            continue;
        }

        const u64 length =
            iv.right - iv.left + 1;

        ++observed_counts[r];

        if (observed_lengths[r] == UINT64_MAX) {
            observed_lengths[r] = length;
        } else if (observed_lengths[r] != length) {
            ++failures;
        }

        if (!safe_add(
                observed_total_hit,
                length,
                observed_total_hit
            )) {
            ++failures;
        }
    }

    for (u64 r = 0;
         r < T;
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
            ++failures;
            continue;
        }

        if (observed_counts[r] != expected_count) {
            ++failures;
        }

        if (expected_count == 0) {
            if (observed_lengths[r] != UINT64_MAX) {
                ++failures;
            }
        } else {
            if (observed_lengths[r] != expected_length) {
                ++failures;
            }
        }
    }

    /*
        Sum_r C_r.
    */
    u64 expected_interval_count = 0;

    for (u64 r = 0;
         r < T;
         ++r) {

        u64 count;

        if (!expected_type_count(
                c,
                r,
                count
            ) ||
            !safe_add(
                expected_interval_count,
                count,
                expected_interval_count
            )) {
            ++failures;
        }
    }

    if (expected_interval_count != c.interval_count) {
        ++failures;
    }

    /*
        Sum_r C_r L_r.
    */
    u64 expected_total_hit = 0;

    for (u64 r = 0;
         r < T;
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

    if (observed_total_hit != expected_total_hit) {
        ++failures;
    }

    /*
        Compare against the complete HIT count from
        the MISS-block theorem.
    */
    u64 miss_count;

    if (!safe_mul(
            c.s0,
            c.interval_count + 1,
            miss_count
        )) {
        ++failures;
    } else {
        const u64 domain_size = c.m + 1;

        if (miss_count > domain_size) {
            ++failures;
        } else {
            const u64 hit_count =
                domain_size - miss_count;

            if (hit_count != expected_total_hit) {
                ++failures;
            }
        }
    }

    return failures == 0;
}

bool check_random_ranks(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_checks,
    u64& failures
) {
    failures = 0;

    /*
        Random checks must NOT modify the exhaustive
        type aggregates.
    */
    for (u64 n = 0;
         n < random_checks;
         ++n) {

        const u64 k =
            rng() % c.interval_count;

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
            continue;
        }

        u64 expected_length;

        if (!expected_type_length(
                c,
                r,
                expected_length
            )) {
            ++failures;
            continue;
        }

        const u64 actual_length =
            iv.right - iv.left + 1;

        if (actual_length != expected_length) {
            ++failures;
        }

        if (!leq_p(
                j,
                c.q,
                c.p
            )) {
            ++failures;
        }

        u64 expected_left;

        u64 jpe;

        if (!safe_mul(
                j,
                c.pe,
                jpe
            ) ||
            !safe_add(
                c.s0,
                jpe,
                expected_left
            )) {
            ++failures;
        } else if (
            iv.left != expected_left
        ) {
            ++failures;
        }
    }

    return failures == 0;
}

bool run_case(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_checks,
    bool do_exhaustive,
    bool print_details
) {
    u64 exhaustive_failures = 0;
    u64 random_failures = 0;

    std::vector<u64> observed_counts;
    std::vector<u64> observed_lengths;

    u64 observed_total_hit = 0;

    bool exhaustive_pass = true;

    if (do_exhaustive) {
        exhaustive_pass =
            check_exhaustive_types(
                c,
                exhaustive_failures,
                observed_counts,
                observed_lengths,
                observed_total_hit
            );
    }

    const bool random_pass =
        check_random_ranks(
            c,
            rng,
            random_checks,
            random_failures
        );

    const bool final_pass =
        exhaustive_pass &&
        random_pass;

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
            << "type_count="
            << c.q_digits.size()
            << '\n';

        std::cout
            << "exhaustive="
            << do_exhaustive
            << '\n';

        std::cout
            << "exhaustive_failures="
            << exhaustive_failures
            << '\n';

        std::cout
            << "random_checks="
            << random_checks
            << '\n';

        std::cout
            << "random_failures="
            << random_failures
            << '\n';

        if (do_exhaustive) {
            std::cout
                << "observed_total_hit="
                << observed_total_hit
                << '\n';

            for (u64 r = 0;
                 r < c.q_digits.size();
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
                    << "]"
                    << " q_digit="
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

        std::cout
            << "exhaustive_pass="
            << exhaustive_pass
            << '\n';

        std::cout
            << "random_pass="
            << random_pass
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << '\n';

        std::cout << '\n';
    }

    return final_pass;
}

int main() {
    std::cout
        << "START EXPERIMENT 226\n\n";

    std::mt19937_64 rng(
        226226226ULL
    );

    u64 total_cases = 0;
    u64 failed_cases = 0;

    auto execute =
        [&](u64 p,
            u64 a,
            u64 b,
            u64 z,
            u64 q,
            u64 random_checks,
            bool exhaustive) {

        CaseData c;

        if (!build_case(
                p,
                a,
                b,
                z,
                q,
                c
            )) {
            std::cout
                << "INVALID_CASE"
                << " p=" << p
                << " a=" << a
                << " b=" << b
                << " z=" << z
                << " q=" << q
                << '\n';

            return;
        }

        ++total_cases;

        const bool pass =
            run_case(
                c,
                rng,
                random_checks,
                exhaustive,
                true
            );

        if (!pass) {
            ++failed_cases;
        }
    };

    /*
        Small exhaustive cases.
    */
    execute(
        2, 0, 0, 2, 63,
        1000,
        true
    );

    execute(
        3, 0, 0, 3, 80,
        1000,
        true
    );

    execute(
        5, 0, 0, 2, 100,
        1000,
        true
    );

    /*
        Valid nontrivial s0.
    */
    execute(
        3, 3, 1, 1, 80,
        1000,
        true
    );

    /*
        Larger nontrivial s0.
    */
    execute(
        3, 4, 1, 6,
        987654321ULL,
        100000,
        false
    );

    /*
        Dense binary.
    */
    execute(
        2, 0, 0, 40,
        1048575,
        100000,
        true
    );

    /*
        Large ternary.
    */
    execute(
        3, 0, 0, 10,
        1000001234567ULL,
        100000,
        false
    );

    /*
        Large base-5.
    */
    execute(
        5, 0, 0, 8,
        1000007654321ULL,
        100000,
        false
    );

    /*
        Sparse binary.
    */
    execute(
        2, 0, 0, 19,
        1073741825ULL,
        100000,
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
        << "C_r = q_r * product_{i>r}(q_i+1)\n";

    std::cout
        << "L_r = (p^r - (q mod p^r)) p^e - s0\n";

    std::cout
        << "sum_r C_r = product_i(q_i+1)-1\n";

    std::cout
        << "sum_r C_r L_r = total HIT count\n";

    std::cout
        << "MISS = s0 * (I+1)\n";

    std::cout
        << "HIT = (m+1) - MISS\n";

    std::cout
        << "random checks do not modify type aggregates\n";

    std::cout
        << "FINISHED EXPERIMENT 226\n";

    return failed_cases == 0 ? 0 : 1;
}
