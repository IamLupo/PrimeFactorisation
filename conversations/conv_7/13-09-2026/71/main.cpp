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

        u64 contribution;

        if (!safe_mul(
                d,
                c.powers[i],
                contribution
            )) {
            return false;
        }

        if (!safe_add(
                j,
                contribution,
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

/*
    Direct length from the interval endpoints:

        L = E-j

    More useful closed form:

        L_k =
        (p^r - (j mod p^r)) p^e - s0
*/
bool length_from_j_formula(
    u64 j,
    u64 r,
    const CaseData& c,
    u64& length
) {
    u64 pr;

    if (!pow_u64(c.p, r, pr)) {
        return false;
    }

    const u64 remainder = j % pr;

    if (remainder > pr) {
        return false;
    }

    const u64 gap_factor =
        pr - remainder;

    u64 raw;

    if (!safe_mul(
            gap_factor,
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

/*
    Alternative closed form using q instead of j:

        L_k =
        (p^r - (q mod p^r)) p^e - s0

    because j_i=q_i for all i<r.
*/
bool length_from_q_formula(
    u64 r,
    const CaseData& c,
    u64& length
) {
    u64 pr;

    if (!pow_u64(c.p, r, pr)) {
        return false;
    }

    const u64 remainder =
        c.q % pr;

    if (remainder > pr) {
        return false;
    }

    const u64 gap_factor =
        pr - remainder;

    u64 raw;

    if (!safe_mul(
            gap_factor,
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

bool check_rank(
    u64 k,
    const CaseData& c,
    u64& failures
) {
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

    const u64 endpoint_length =
        iv.right - iv.left + 1;

    u64 j_length;

    if (!length_from_j_formula(
            j,
            r,
            c,
            j_length
        )) {
        ++failures;
        return false;
    }

    if (endpoint_length != j_length) {
        ++failures;
        return false;
    }

    u64 q_length;

    if (!length_from_q_formula(
            r,
            c,
            q_length
        )) {
        ++failures;
        return false;
    }

    if (j_length != q_length) {
        ++failures;
        return false;
    }

    /*
        Therefore all three descriptions agree:

            endpoint length
            j-based length
            q-based length
    */
    return true;
}

bool run_case(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_checks,
    bool print_details
) {
    u64 failures = 0;

    u64 exhaustive_checks = 0;

    if (c.interval_count <= 2000000) {
        for (u64 k = 0;
             k < c.interval_count;
             ++k) {

            ++exhaustive_checks;

            check_rank(
                k,
                c,
                failures
            );
        }
    }

    u64 random_checks_done = 0;

    for (u64 n = 0;
         n < random_checks;
         ++n) {

        ++random_checks_done;

        const u64 k =
            rng() % c.interval_count;

        check_rank(
            k,
            c,
            failures
        );
    }

    /*
        Check first/last ranks explicitly.
    */
    bool first_pass = true;
    bool last_pass = true;

    {
        u64 j;
        u64 r;
        Interval iv;

        if (!interval_from_rank(
                0,
                c,
                j,
                r,
                iv
            )) {
            first_pass = false;
            ++failures;
        } else {
            if (iv.left != c.s0) {
                first_pass = false;
                ++failures;
            }
        }
    }

    {
        u64 j;
        u64 r;
        Interval iv;

        const u64 k =
            c.interval_count - 1;

        if (!interval_from_rank(
                k,
                c,
                j,
                r,
                iv
            )) {
            last_pass = false;
            ++failures;
        } else {
            if (iv.right != c.m - c.s0) {
                last_pass = false;
                ++failures;
            }
        }
    }

    const bool pass =
        failures == 0;

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
            << "random_checks="
            << random_checks_done
            << '\n';

        std::cout
            << "failures="
            << failures
            << '\n';

        std::cout
            << "first_pass="
            << first_pass
            << '\n';

        std::cout
            << "last_pass="
            << last_pass
            << '\n';

        std::cout
            << "final_pass="
            << pass
            << '\n';

        const u64 sample_count =
            std::min<u64>(
                10,
                c.interval_count
            );

        for (u64 k = 0;
             k < sample_count;
             ++k) {

            u64 j;
            u64 r;
            Interval iv;

            interval_from_rank(
                k,
                c,
                j,
                r,
                iv
            );

            const u64 actual_length =
                iv.right - iv.left + 1;

            u64 j_length = 0;
            u64 q_length = 0;

            length_from_j_formula(
                j,
                r,
                c,
                j_length
            );

            length_from_q_formula(
                r,
                c,
                q_length
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
                << "]"
                << " actual_length="
                << actual_length
                << " j_length="
                << j_length
                << " q_length="
                << q_length
                << '\n';
        }

        std::cout << '\n';
    }

    return pass;
}

bool build_and_run(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    std::mt19937_64& rng,
    u64 random_checks,
    bool print_details
) {
    CaseData c;

    if (!build_case(
            p,
            a,
            b,
            z,
            q,
            c
        )) {
        return false;
    }

    return run_case(
        c,
        rng,
        random_checks,
        print_details
    );
}

int main() {
    std::cout
        << "START EXPERIMENT 224\n\n";

    std::mt19937_64 rng(
        224224224ULL
    );

    u64 total_cases = 0;
    u64 failed_cases = 0;

    auto execute = [&](u64 p,
                       u64 a,
                       u64 b,
                       u64 z,
                       u64 q,
                       u64 random_checks,
                       bool print_details) {
        ++total_cases;

        const bool pass =
            build_and_run(
                p,
                a,
                b,
                z,
                q,
                rng,
                random_checks,
                print_details
            );

        if (!pass) {
            ++failed_cases;
        }
    };

    /*
        Small / exhaustive.
    */
    execute(
        2, 0, 0, 2, 63,
        1000, true
    );

    execute(
        3, 0, 0, 3, 80,
        1000, true
    );

    execute(
        5, 0, 0, 2, 100,
        1000, true
    );

    /*
        Nontrivial s0.
    */
    execute(
        3, 3, 2, 1, 80,
        1000, true
    );

    /*
        Dense large binary.
    */
    execute(
        2, 0, 0, 40, 1048575,
        100000, true
    );

    /*
        Large ternary.
    */
    execute(
        3, 0, 0, 10, 1000001234567ULL,
        100000, true
    );

    /*
        Large base 5.
    */
    execute(
        5, 0, 0, 8, 1000007654321ULL,
        100000, true
    );

    /*
        Large nontrivial s0.
    */
    execute(
        3, 4, 1, 6, 987654321ULL,
        100000, true
    );

    /*
        Sparse binary.
    */
    execute(
        2, 0, 0, 19, 1073741825ULL,
        100000, true
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

    std::cout
        << '\n';

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "L_k = E_k - x_k + 1\n";

    std::cout
        << "L_k = (p^r - (j_k mod p^r)) p^e - s0\n";

    std::cout
        << "L_k = (p^r - (q mod p^r)) p^e - s0\n";

    std::cout
        << "j_k and q agree in digits 0,...,r-1\n";

    std::cout
        << "endpoint length = j-based length = q-based length\n";

    std::cout
        << "FINISHED EXPERIMENT 224\n";

    return failed_cases == 0 ? 0 : 1;
}
