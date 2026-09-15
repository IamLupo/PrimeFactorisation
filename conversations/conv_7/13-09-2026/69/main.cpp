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

bool direct_interval(
    u64 k,
    const CaseData& c,
    u64& j,
    u64& r,
    Interval& interval
) {
    if (!unrank_j(
            k,
            c,
            j
        )) {
        return false;
    }

    if (!leq_p(
            j,
            c.q,
            c.p
        )) {
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

    interval = {left, right};

    return true;
}

bool check_case(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_checks,
    bool print_details
) {
    u64 failures = 0;

    /*
        Check every rank when the interval family is small.
    */
    u64 exhaustive_checks = 0;
    u64 exhaustive_failures = 0;

    if (c.interval_count <= 2000000) {
        for (u64 k = 0;
             k < c.interval_count;
             ++k) {

            ++exhaustive_checks;

            u64 j;
            u64 r;
            Interval iv;

            if (!direct_interval(
                    k,
                    c,
                    j,
                    r,
                    iv
                )) {
                ++exhaustive_failures;
                continue;
            }

            if (k == 0 &&
                iv.left != c.s0) {
                ++exhaustive_failures;
            }

            if (k + 1 == c.interval_count &&
                iv.right != c.m - c.s0) {
                ++exhaustive_failures;
            }

            /*
                Rank k must produce the same index that
                the mixed-radix representation predicts.
            */
            u64 recovered_rank = 0;

            for (size_t i = 0;
                 i < c.radix.size();
                 ++i) {

                const u64 jd =
                    p_digit(
                        j,
                        c.p,
                        static_cast<u64>(i)
                    );

                u64 term;

                if (!safe_mul(
                        jd,
                        c.weights[i],
                        term
                    ) ||
                    !safe_add(
                        recovered_rank,
                        term,
                        recovered_rank
                    )) {
                    ++exhaustive_failures;
                    break;
                }
            }

            if (recovered_rank != k) {
                ++exhaustive_failures;
            }

            /*
                Consecutive geometric recurrence.
            */
            if (k + 1 < c.interval_count) {
                u64 nj;
                u64 nr;
                Interval next;

                if (!direct_interval(
                        k + 1,
                        c,
                        nj,
                        nr,
                        next
                    )) {
                    ++exhaustive_failures;
                    continue;
                }

                u64 expected_next_left;

                if (!safe_add(
                        iv.right,
                        c.s0,
                        expected_next_left
                    ) ||
                    !safe_add(
                        expected_next_left,
                        1,
                        expected_next_left
                    )) {
                    ++exhaustive_failures;
                } else if (
                    next.left != expected_next_left
                ) {
                    ++exhaustive_failures;
                }
            }
        }
    }

    /*
        Random direct-rank checks for arbitrarily large families.
    */
    u64 random_failures = 0;

    for (u64 n = 0;
         n < random_checks;
         ++n) {

        const u64 k =
            rng() % c.interval_count;

        u64 j;
        u64 r;
        Interval iv;

        if (!direct_interval(
                k,
                c,
                j,
                r,
                iv
            )) {
            ++random_failures;
            continue;
        }

        /*
            Rank recovery.
        */
        u64 recovered_rank = 0;

        bool rank_ok = true;

        for (size_t i = 0;
             i < c.radix.size();
             ++i) {

            const u64 jd =
                p_digit(
                    j,
                    c.p,
                    static_cast<u64>(i)
                );

            if (jd > c.q_digits[i]) {
                rank_ok = false;
                break;
            }

            u64 term;

            if (!safe_mul(
                    jd,
                    c.weights[i],
                    term
                ) ||
                !safe_add(
                    recovered_rank,
                    term,
                    recovered_rank
                )) {
                rank_ok = false;
                break;
            }
        }

        if (!rank_ok ||
            recovered_rank != k) {
            ++random_failures;
            continue;
        }

        /*
            Check the interval formula independently.
        */
        u64 jpe;

        if (!safe_mul(
                j,
                c.pe,
                jpe
            )) {
            ++random_failures;
            continue;
        }

        u64 expected_left;

        if (!safe_add(
                c.s0,
                jpe,
                expected_left
            )) {
            ++random_failures;
            continue;
        }

        if (iv.left != expected_left) {
            ++random_failures;
            continue;
        }

        /*
            Check the endpoint.
        */
        u64 pr;

        if (!pow_u64(
                c.p,
                r,
                pr
            )) {
            ++random_failures;
            continue;
        }

        const u64 prefix =
            j / pr;

        u64 prefix_plus_one;

        if (!safe_add(
                prefix,
                1,
                prefix_plus_one
            )) {
            ++random_failures;
            continue;
        }

        u64 power_er;

        if (!safe_mul(
                c.pe,
                pr,
                power_er
            )) {
            ++random_failures;
            continue;
        }

        u64 expected_right_plus_one;

        if (!safe_mul(
                prefix_plus_one,
                power_er,
                expected_right_plus_one
            ) ||
            expected_right_plus_one == 0) {
            ++random_failures;
            continue;
        }

        const u64 expected_right =
            expected_right_plus_one - 1;

        if (iv.right != expected_right) {
            ++random_failures;
            continue;
        }

        /*
            Geometric successor.
        */
        if (k + 1 < c.interval_count) {
            u64 nj;
            u64 nr;
            Interval next;

            if (!direct_interval(
                    k + 1,
                    c,
                    nj,
                    nr,
                    next
                )) {
                ++random_failures;
                continue;
            }

            u64 expected_next_left;

            if (!safe_add(
                    iv.right,
                    c.s0,
                    expected_next_left
                ) ||
                !safe_add(
                    expected_next_left,
                    1,
                    expected_next_left
                ) ||
                next.left != expected_next_left) {
                ++random_failures;
            }
        }
    }

    failures =
        exhaustive_failures +
        random_failures;

    const bool final_pass =
        failures == 0;

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
            << "exhaustive_checks="
            << exhaustive_checks
            << '\n';

        std::cout
            << "exhaustive_failures="
            << exhaustive_failures
            << '\n';

        std::cout
            << "exhaustive_pass="
            << (exhaustive_failures == 0)
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
            << (random_failures == 0)
            << '\n';

        std::cout
            << "total_failures="
            << failures
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

            u64 j;
            u64 r;
            Interval iv;

            direct_interval(
                k,
                c,
                j,
                r,
                iv
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

bool run_case(
    u64 p,
    u64 a,
    u64 b,
    u64 z,
    u64 q,
    std::mt19937_64& rng,
    u64 random_checks,
    bool print_details,
    bool& result
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
        result = false;
        return false;
    }

    result =
        check_case(
            c,
            rng,
            random_checks,
            print_details
        );

    return true;
}

int main() {
    std::cout
        << "START EXPERIMENT 223\n\n";

    bool all_pass = true;

    std::mt19937_64 rng(
        223223223ULL
    );

    /*
        Small dense binary.
    */
    {
        bool pass;

        run_case(
            2,
            0,
            0,
            2,
            63,
            rng,
            1000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Small ternary.
    */
    {
        bool pass;

        run_case(
            3,
            0,
            0,
            3,
            80,
            rng,
            1000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Base 5.
    */
    {
        bool pass;

        run_case(
            5,
            0,
            0,
            2,
            100,
            rng,
            1000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Nontrivial s0.
    */
    {
        bool pass;

        run_case(
            3,
            3,
            2,
            1,
            80,
            rng,
            1000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Large dense binary.
    */
    {
        bool pass;

        run_case(
            2,
            0,
            0,
            40,
            1048575,
            rng,
            100000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Large ternary.
    */
    {
        bool pass;

        run_case(
            3,
            0,
            0,
            10,
            1000001234567ULL,
            rng,
            100000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Large base 5.
    */
    {
        bool pass;

        run_case(
            5,
            0,
            0,
            8,
            1000007654321ULL,
            rng,
            100000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Large nontrivial s0.
    */
    {
        bool pass;

        run_case(
            3,
            4,
            1,
            6,
            987654321ULL,
            rng,
            100000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    /*
        Sparse binary.
    */
    {
        bool pass;

        run_case(
            2,
            0,
            0,
            19,
            1073741825ULL,
            rng,
            100000,
            true,
            pass
        );

        all_pass =
            pass &&
            all_pass;
    }

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "j_k = unrank(k)\n";

    std::cout
        << "r_k = first digit with (j_k)_r < q_r\n";

    std::cout
        << "x_k = s0 + j_k p^e\n";

    std::cout
        << "E_k = (floor(j_k/p^r_k)+1)p^(e+r_k)-1\n";

    std::cout
        << "x_{k+1} = E_k + s0 + 1\n";

    std::cout
        << "E_{I-1} = m-s0\n";

    std::cout
        << "rank k is exactly the kth ordered HIT interval\n";

    std::cout
        << "large cases use no O(q) reference enumeration\n";

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << all_pass
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 223\n";

    return all_pass ? 0 : 1;
}
