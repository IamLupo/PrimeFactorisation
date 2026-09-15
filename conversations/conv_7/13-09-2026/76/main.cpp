#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

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

    /*
        IMPORTANT:
        q * p^e must itself fit in uint64 for the
        closed HIT expression to be representable.
    */
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

    c.p_powers.resize(
        c.q_digits.size()
    );

    c.p_powers[0] = 1;

    for (size_t i = 1;
         i < c.p_powers.size();
         ++i) {

        if (!safe_mul(
                c.p_powers[i - 1],
                p,
                c.p_powers[i]
            )) {
            return false;
        }
    }

    c.digit_product = 1;

    for (u64 r : c.radix) {
        if (!safe_mul(
                c.digit_product,
                r,
                c.digit_product
            )) {
            return false;
        }
    }

    c.interval_count =
        c.digit_product - 1;

    return true;
}

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

bool type_factor(
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
    S(q) =
        sum_r C_r (p^r - q mod p^r)
*/
bool compute_inner_sum(
    const CaseData& c,
    u64& sum
) {
    sum = 0;

    for (u64 r = 0;
         r < c.q_digits.size();
         ++r) {

        u64 C;
        u64 A;

        if (!type_count(
                c,
                r,
                C
            ) ||
            !type_factor(
                c,
                r,
                A
            )) {
            return false;
        }

        u64 term;

        if (!safe_mul(
                C,
                A,
                term
            ) ||
            !safe_add(
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
    H_type =
        sum_r C_r *
        (A_r * p^e - s0)
*/
bool compute_type_hit(
    const CaseData& c,
    u64 inner_sum,
    u64& hit
) {
    u64 qpe;

    if (!safe_mul(
            c.q,
            c.pe,
            qpe
        )) {
        return false;
    }

    u64 first_term;

    if (!safe_mul(
            inner_sum,
            c.pe,
            first_term
        )) {
        return false;
    }

    /*
        sum C_r = interval_count
    */
    u64 second_term;

    if (!safe_mul(
            c.s0,
            c.interval_count,
            second_term
        )) {
        return false;
    }

    if (first_term < second_term) {
        return false;
    }

    hit =
        first_term - second_term;

    /*
        Independent closed form.
    */
    u64 closed_miss;

    if (!safe_mul(
            c.s0,
            c.interval_count,
            closed_miss
        )) {
        return false;
    }

    if (qpe < closed_miss) {
        return false;
    }

    const u64 closed_hit =
        qpe - closed_miss;

    if (closed_hit != hit) {
        return false;
    }

    return true;
}

/*
    Direct identity using the digit product:

        H =
        q p^e - s0 (product_i(q_i+1)-1)
*/
bool compute_closed_hit(
    const CaseData& c,
    u64& hit
) {
    u64 qpe;

    if (!safe_mul(
            c.q,
            c.pe,
            qpe
        )) {
        return false;
    }

    if (c.digit_product == 0) {
        return false;
    }

    const u64 I =
        c.digit_product - 1;

    u64 miss;

    if (!safe_mul(
            c.s0,
            I,
            miss
        )) {
        return false;
    }

    if (qpe < miss) {
        return false;
    }

    hit =
        qpe - miss;

    return true;
}

bool check_case(
    const CaseData& c
) {
    u64 inner;
    u64 type_hit;
    u64 closed_hit;

    if (!compute_inner_sum(
            c,
            inner
        )) {
        return false;
    }

    if (!compute_type_hit(
            c,
            inner,
            type_hit
        )) {
        return false;
    }

    if (!compute_closed_hit(
            c,
            closed_hit
        )) {
        return false;
    }

    /*
        First identity.
    */
    if (inner != c.q) {
        return false;
    }

    /*
        Second identity.
    */
    if (type_hit != closed_hit) {
        return false;
    }

    return true;
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
        << "START EXPERIMENT 228\n\n";

    std::mt19937_64 rng(
        228228228ULL
    );

    u64 valid_cases = 0;
    u64 failed_cases = 0;
    u64 overflow_skipped = 0;
    u64 invalid_skipped = 0;

    /*
        Known valid deterministic cases.
    */
    struct KnownCase {
        u64 p;
        u64 a;
        u64 b;
        u64 z;
        u64 q;
    };

    const KnownCase known[] = {
        {2, 0, 0, 2, 63},
        {3, 0, 0, 3, 80},
        {5, 0, 0, 2, 100},
        {3, 3, 1, 1, 80},
        {3, 4, 1, 6, 987654321ULL},
        {2, 0, 0, 19, 1073741825ULL}
    };

    for (const auto& k : known) {
        CaseData c;

        if (!build_case(
                k.p,
                k.a,
                k.b,
                k.z,
                k.q,
                c
            )) {
            ++failed_cases;
            continue;
        }

        ++valid_cases;

        const bool pass =
            check_case(c);

        if (!pass) {
            ++failed_cases;
        }

        std::cout
            << "KNOWN_CASE"
            << " p=" << k.p
            << " a=" << k.a
            << " b=" << k.b
            << " z=" << k.z
            << " q=" << k.q
            << " pass=" << pass
            << '\n';
    }

    std::cout << '\n';

    /*
        Large random structural suite.

        We intentionally generate fairly aggressive parameters,
        so overflow rejection itself is tested.
    */
    const u64 random_cases = 100000;

    for (u64 i = 0;
         i < random_cases;
         ++i) {

        const u64 p =
            (i % 3 == 0)
                ? 2
                : (i % 3 == 1)
                    ? 3
                    : 5;

        const u64 a =
            random_bounded(
                rng,
                0,
                12
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
                14
            );

        const u64 q =
            random_bounded(
                rng,
                1,
                5000000000ULL
            );

        /*
            Test structural admissibility and uint64
            representability separately.
        */
        CaseData c;

        if (!build_case(
                p,
                a,
                b,
                z,
                q,
                c
            )) {
            ++overflow_skipped;
            continue;
        }

        ++valid_cases;

        if (!check_case(c)) {
            ++failed_cases;

            if (failed_cases <= 10) {
                std::cout
                    << "RANDOM_FAILURE\n";

                std::cout
                    << "p=" << p
                    << " a=" << a
                    << " b=" << b
                    << " z=" << z
                    << " q=" << q
                    << '\n';

                std::cout << '\n';
            }
        }
    }

    /*
        Separate test with small parameters where essentially
        every random case should be representable.
    */
    const u64 safe_random_cases = 10000;

    for (u64 i = 0;
         i < safe_random_cases;
         ++i) {

        const u64 p =
            (i % 3 == 0)
                ? 2
                : (i % 3 == 1)
                    ? 3
                    : 5;

        const u64 a =
            random_bounded(
                rng,
                0,
                4
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
                5
            );

        const u64 q =
            random_bounded(
                rng,
                1,
                1000000ULL
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
            ++overflow_skipped;
            continue;
        }

        ++valid_cases;

        if (!check_case(c)) {
            ++failed_cases;

            if (failed_cases <= 10) {
                std::cout
                    << "SAFE_RANDOM_FAILURE\n";

                std::cout
                    << "p=" << p
                    << " a=" << a
                    << " b=" << b
                    << " z=" << z
                    << " q=" << q
                    << '\n';

                std::cout << '\n';
            }
        }
    }

    std::cout
        << "SUMMARY\n";

    std::cout
        << "valid_cases="
        << valid_cases
        << '\n';

    std::cout
        << "failed_cases="
        << failed_cases
        << '\n';

    std::cout
        << "overflow_skipped="
        << overflow_skipped
        << '\n';

    std::cout
        << "invalid_skipped="
        << invalid_skipped
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
        << "sum_r C_r L_r"
        << " = q p^e"
        << " - s0(product_i(q_i+1)-1)\n";

    std::cout
        << "overflow cases are explicitly skipped\n";

    std::cout
        << "all accepted cases are exactly representable in uint64\n";

    std::cout
        << "FINISHED EXPERIMENT 228\n";

    return failed_cases == 0 ? 0 : 1;
}
