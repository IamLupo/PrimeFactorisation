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
    u64 pe;

    std::vector<u64> q_digits;
    std::vector<u64> radix;
    std::vector<u64> p_powers;

    u64 digit_product;
    u64 interval_count;
};

/*
    Overflow-safe arithmetic.
*/
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

/*
    Build the structural parameterization:

        s0 = (b+1)p^a
        e  = a+z+1
*/
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

    c = {};

    c.p = p;
    c.a = a;
    c.b = b;
    c.z = z;
    c.e = e;
    c.s0 = s0;
    c.q = q;
    c.pe = pe;

    /*
        Base-p digits of q.
    */
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

    /*
        p^r.
    */
    c.p_powers.resize(c.q_digits.size());

    c.p_powers[0] = 1;

    for (size_t r = 1;
         r < c.p_powers.size();
         ++r) {

        if (!safe_mul(
                c.p_powers[r - 1],
                p,
                c.p_powers[r]
            )) {
            return false;
        }
    }

    /*
        Product_i(q_i+1).
    */
    c.digit_product = 1;

    for (u64 radix : c.radix) {
        if (!safe_mul(
                c.digit_product,
                radix,
                c.digit_product
            )) {
            return false;
        }
    }

    c.interval_count =
        c.digit_product - 1;

    return true;
}

/*
    C_r = q_r * product_{i>r}(q_i+1)
*/
bool type_count(
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
    A_r = p^r - (q mod p^r)
*/
bool type_gap_factor(
    const CaseData& c,
    u64 r,
    u64& factor
) {
    if (r >= c.p_powers.size()) {
        return false;
    }

    const u64 pr =
        c.p_powers[r];

    const u64 remainder =
        c.q % pr;

    factor =
        pr - remainder;

    return true;
}

/*
    L_r =
        A_r * p^e - s0
*/
bool type_length(
    const CaseData& c,
    u64 r,
    u64& length
) {
    u64 factor;

    if (!type_gap_factor(
            c,
            r,
            factor
        )) {
        return false;
    }

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

/*
    The inner identity:

        S(q) =
        sum_r C_r (p^r - q mod p^r)

    Expected result:

        S(q) = q
*/
bool compute_inner_sum(
    const CaseData& c,
    u64& sum
) {
    sum = 0;

    for (u64 r = 0;
         r < c.q_digits.size();
         ++r) {

        u64 count;

        if (!type_count(
                c,
                r,
                count
            )) {
            return false;
        }

        u64 factor;

        if (!type_gap_factor(
                c,
                r,
                factor
            )) {
            return false;
        }

        u64 term;

        if (!safe_mul(
                count,
                factor,
                term
            )) {
            return false;
        }

        if (!safe_add(
                sum,
                term,
                sum
            )) {
            return false;
        }
    }

    return true;
}

/*
    Explicit weighted type sum:

        H_type = sum_r C_r L_r
*/
bool compute_type_hit_sum(
    const CaseData& c,
    u64& sum
) {
    sum = 0;

    for (u64 r = 0;
         r < c.q_digits.size();
         ++r) {

        u64 count;

        if (!type_count(
                c,
                r,
                count
            )) {
            return false;
        }

        u64 length;

        if (!type_length(
                c,
                r,
                length
            )) {
            return false;
        }

        u64 term;

        if (!safe_mul(
                count,
                length,
                term
            )) {
            return false;
        }

        if (!safe_add(
                sum,
                term,
                sum
            )) {
            return false;
        }
    }

    return true;
}

/*
    Closed form:

        H_closed =
        q p^e - s0 (product_i(q_i+1)-1)
*/
bool compute_closed_hit(
    const CaseData& c,
    u64& result
) {
    u64 qpe;

    if (!safe_mul(
            c.q,
            c.pe,
            qpe
        )) {
        return false;
    }

    u64 miss_multiplier =
        c.digit_product - 1;

    u64 miss_term;

    if (!safe_mul(
            c.s0,
            miss_multiplier,
            miss_term
        )) {
        return false;
    }

    if (qpe < miss_term) {
        return false;
    }

    result =
        qpe - miss_term;

    return true;
}

/*
    Alternative derivation:

        H_type
        = p^e * sum C_r A_r
          - s0 * sum C_r

    where

        sum C_r = product(q_i+1)-1

        sum C_r A_r = q.
*/
bool compute_decomposed_hit(
    const CaseData& c,
    u64 inner_sum,
    u64& result
) {
    u64 weighted_q;

    if (!safe_mul(
            inner_sum,
            c.pe,
            weighted_q
        )) {
        return false;
    }

    u64 miss_term;

    if (!safe_mul(
            c.s0,
            c.interval_count,
            miss_term
        )) {
        return false;
    }

    if (weighted_q < miss_term) {
        return false;
    }

    result =
        weighted_q - miss_term;

    return true;
}

/*
    Check one case.
*/
bool check_case(
    const CaseData& c,
    bool print_details
) {
    u64 inner_sum = 0;
    u64 type_hit = 0;
    u64 closed_hit = 0;
    u64 decomposed_hit = 0;

    bool inner_ok =
        compute_inner_sum(
            c,
            inner_sum
        );

    bool type_ok =
        compute_type_hit_sum(
            c,
            type_hit
        );

    bool closed_ok =
        compute_closed_hit(
            c,
            closed_hit
        );

    bool decomposed_ok =
        inner_ok &&
        compute_decomposed_hit(
            c,
            inner_sum,
            decomposed_hit
        );

    const bool inner_identity_pass =
        inner_ok &&
        inner_sum == c.q;

    const bool type_closed_pass =
        type_ok &&
        closed_ok &&
        type_hit == closed_hit;

    const bool decomposition_pass =
        type_ok &&
        decomposed_ok &&
        type_hit == decomposed_hit;

    /*
        Complete HIT count from the previous MISS theorem.
    */
    u64 expected_hit_from_domain = 0;

    bool domain_pass = false;

    u64 miss_count = 0;

    if (safe_mul(
            c.s0,
            c.interval_count + 1,
            miss_count
        )) {

        const u64 domain_size =
            /*
                m+1 = s0 + q*p^e
            */
            c.s0;

        u64 qpe;

        if (safe_mul(
                c.q,
                c.pe,
                qpe
            ) &&
            safe_add(
                domain_size,
                qpe,
                expected_hit_from_domain
            )) {

            if (miss_count <= expected_hit_from_domain) {
                expected_hit_from_domain -= miss_count;
                domain_pass =
                    true;
            }
        }
    }

    const bool domain_hit_pass =
        domain_pass &&
        closed_ok &&
        closed_hit == expected_hit_from_domain;

    const bool final_pass =
        inner_identity_pass &&
        type_closed_pass &&
        decomposition_pass &&
        domain_hit_pass;

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
            << "digit_count="
            << c.q_digits.size()
            << '\n';

        std::cout
            << "digit_product="
            << c.digit_product
            << '\n';

        std::cout
            << "interval_count="
            << c.interval_count
            << '\n';

        std::cout
            << "inner_sum="
            << inner_sum
            << '\n';

        std::cout
            << "inner_expected_q="
            << c.q
            << '\n';

        std::cout
            << "inner_identity_pass="
            << inner_identity_pass
            << '\n';

        std::cout
            << "type_hit="
            << type_hit
            << '\n';

        std::cout
            << "closed_hit="
            << closed_hit
            << '\n';

        std::cout
            << "type_closed_pass="
            << type_closed_pass
            << '\n';

        std::cout
            << "decomposed_hit="
            << decomposed_hit
            << '\n';

        std::cout
            << "decomposition_pass="
            << decomposition_pass
            << '\n';

        std::cout
            << "domain_hit="
            << expected_hit_from_domain
            << '\n';

        std::cout
            << "domain_hit_pass="
            << domain_hit_pass
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << '\n';

        /*
            Print the complete algebraic type table for small q.
        */
        if (c.q_digits.size() <= 12) {
            for (u64 r = 0;
                 r < c.q_digits.size();
                 ++r) {

                u64 count;
                u64 factor;
                u64 length;

                type_count(
                    c,
                    r,
                    count
                );

                type_gap_factor(
                    c,
                    r,
                    factor
                );

                type_length(
                    c,
                    r,
                    length
                );

                std::cout
                    << "TYPE["
                    << r
                    << "]"
                    << " q_digit="
                    << c.q_digits[r]
                    << " C_r="
                    << count
                    << " A_r="
                    << factor
                    << " L_r="
                    << length
                    << '\n';
            }
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
    std::uniform_int_distribution<u64> dist(
        lo,
        hi
    );

    return dist(rng);
}

int main() {
    std::cout
        << "START EXPERIMENT 227\n\n";

    std::mt19937_64 rng(
        227227227ULL
    );

    u64 total_cases = 0;
    u64 failed_cases = 0;
    u64 skipped_cases = 0;

    auto execute =
        [&](u64 p,
            u64 a,
            u64 b,
            u64 z,
            u64 q,
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

            std::cout
                << "SKIPPED_INVALID_CASE"
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
            check_case(
                c,
                print_details
            );

        if (!pass) {
            ++failed_cases;
        }
    };

    /*
        Small exact cases.
    */
    execute(
        2, 0, 0, 2, 63,
        true
    );

    execute(
        3, 0, 0, 3, 80,
        true
    );

    execute(
        5, 0, 0, 2, 100,
        true
    );

    execute(
        3, 3, 1, 1, 80,
        true
    );

    /*
        Nontrivial large s0.
    */
    execute(
        3, 4, 1, 6,
        987654321ULL,
        true
    );

    /*
        Large dense binary.
    */
    execute(
        2, 0, 0, 40,
        1048575,
        false
    );

    /*
        Large ternary.
    */
    execute(
        3, 0, 0, 10,
        1000001234567ULL,
        false
    );

    /*
        Large base 5.
    */
    execute(
        5, 0, 0, 8,
        1000007654321ULL,
        false
    );

    /*
        Sparse binary.
    */
    execute(
        2, 0, 0, 19,
        1073741825ULL,
        true
    );

    /*
        Random structurally valid cases.
    */
    for (int i = 0; i < 1000; ++i) {
        const u64 p_choice =
            (i % 3 == 0)
                ? 2
                : (i % 3 == 1)
                    ? 3
                    : 5;

        const u64 a =
            random_bounded(
                rng,
                0,
                8
            );

        const u64 b =
            random_bounded(
                rng,
                0,
                p_choice - 2
            );

        const u64 z =
            random_bounded(
                rng,
                0,
                10
            );

        /*
            Moderate q so everything remains inside uint64_t.
        */
        const u64 q =
            random_bounded(
                rng,
                1,
                1000000000ULL
            );

        CaseData c;

        if (!build_case(
                p_choice,
                a,
                b,
                z,
                q,
                c
            )) {
            ++skipped_cases;
            continue;
        }

        ++total_cases;

        if (!check_case(
                c,
                false
            )) {
            ++failed_cases;

            if (failed_cases <= 3) {
                std::cout
                    << "RANDOM_FAILURE\n";

                std::cout
                    << "p=" << p_choice
                    << " a=" << a
                    << " b=" << b
                    << " z=" << z
                    << " q=" << q
                    << '\n';
            }
        }
    }

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
        << "S(q) = sum_r C_r (p^r - (q mod p^r))\n";

    std::cout
        << "S(q) = q\n";

    std::cout
        << "C_r = q_r product_{i>r}(q_i+1)\n";

    std::cout
        << "H = sum_r C_r L_r\n";

    std::cout
        << "H = q p^e - s0(product_i(q_i+1)-1)\n";

    std::cout
        << "H = (m+1) - s0(product_i(q_i+1))\n";

    std::cout
        << "no interval enumeration is used\n";

    std::cout
        << "FINISHED EXPERIMENT 227\n";

    return failed_cases == 0 ? 0 : 1;
}
