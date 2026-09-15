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

/* ------------------------------------------------------------
   Basic arithmetic helpers
   ------------------------------------------------------------ */

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

/* ------------------------------------------------------------
   Structural case construction
   ------------------------------------------------------------ */

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

/* ------------------------------------------------------------
   Base-p digit operations
   ------------------------------------------------------------ */

u64 digit_at(u64 x, u64 p, u64 i) {
    while (i > 0) {
        x /= p;
        --i;
    }

    return x % p;
}

u64 digit_count(u64 x, u64 p) {
    u64 count = 0;

    do {
        ++count;
        x /= p;
    } while (x != 0);

    return count;
}

/* ------------------------------------------------------------
   Lucas admissibility
   ------------------------------------------------------------ */

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

std::vector<u64> reference_indices(const CaseData& c) {
    std::vector<u64> result;

    for (u64 j = 0; j < c.q; ++j) {
        if (leq_p(j, c.q, c.p)) {
            result.push_back(j);
        }
    }

    return result;
}

/* ------------------------------------------------------------
   Compute mixed-radix digit capacities q_i + 1
   ------------------------------------------------------------ */

bool radix_sizes(
    const CaseData& c,
    std::vector<u64>& radix,
    u64& tuple_count
) {
    radix.clear();

    u64 q = c.q;

    while (true) {
        const u64 q_digit =
            q % c.p;

        u64 base;

        if (!safe_add(q_digit, 1, base)) {
            return false;
        }

        radix.push_back(base);

        q /= c.p;

        if (q == 0) {
            break;
        }
    }

    tuple_count = 1;

    for (u64 base : radix) {
        if (!safe_mul(tuple_count, base, tuple_count)) {
            return false;
        }
    }

    return true;
}

/* ------------------------------------------------------------
   Rank

   R(j) = sum_i j_i * product_{t<i} (q_t + 1)

   Least-significant digit varies fastest.
   ------------------------------------------------------------ */

bool rank_index(
    u64 j,
    const CaseData& c,
    u64& rank
) {
    std::vector<u64> radix;
    u64 tuple_count;

    if (!radix_sizes(
            c,
            radix,
            tuple_count
        )) {
        return false;
    }

    rank = 0;

    u64 weight = 1;

    for (size_t i = 0; i < radix.size(); ++i) {
        const u64 jd =
            digit_at(j, c.p, static_cast<u64>(i));

        if (jd >= radix[i]) {
            return false;
        }

        u64 term;

        if (!safe_mul(
                jd,
                weight,
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

        if (i + 1 < radix.size()) {
            if (!safe_mul(
                    weight,
                    radix[i],
                    weight
                )) {
                return false;
            }
        }
    }

    return true;
}

/* ------------------------------------------------------------
   Unrank

   k -> digits:

       j_i =
       floor(k / weight_i) mod (q_i + 1)

   with weight_i = product_{t<i}(q_t+1).
   ------------------------------------------------------------ */

bool unrank_index(
    u64 rank,
    const CaseData& c,
    u64& j
) {
    std::vector<u64> radix;
    u64 tuple_count;

    if (!radix_sizes(
            c,
            radix,
            tuple_count
        )) {
        return false;
    }

    if (rank >= tuple_count) {
        return false;
    }

    j = 0;

    u64 weight = 1;

    for (size_t i = 0; i < radix.size(); ++i) {
        const u64 digit =
            (rank / weight) % radix[i];

        u64 power;

        if (!pow_u64(
                c.p,
                static_cast<u64>(i),
                power
            )) {
            return false;
        }

        u64 contribution;

        if (!safe_mul(
                digit,
                power,
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

        if (i + 1 < radix.size()) {
            if (!safe_mul(
                    weight,
                    radix[i],
                    weight
                )) {
                return false;
            }
        }
    }

    return true;
}

/* ------------------------------------------------------------
   Generate the j-sequence through the successor map
   ------------------------------------------------------------ */

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

bool successor(
    u64 j,
    const CaseData& c,
    u64& next
) {
    u64 r;

    if (!first_lower_digit(
            j,
            c.q,
            c.p,
            r
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

    if (!safe_mul(
            prefix_plus_one,
            pr,
            next
        )) {
        return false;
    }

    if (next <= j || next >= c.q) {
        return false;
    }

    return leq_p(
        next,
        c.q,
        c.p
    );
}

std::vector<u64> successor_sequence(
    const CaseData& c,
    bool& construction_ok
) {
    construction_ok = true;

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

        if (!successor(
                current,
                c,
                next
            )) {
            break;
        }

        result.push_back(next);
        current = next;
    }

    return result;
}

/* ------------------------------------------------------------
   Sequence comparison
   ------------------------------------------------------------ */

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

/* ------------------------------------------------------------
   Verify:
       rank(j_k) = k
   ------------------------------------------------------------ */

bool verify_rank_sequence(
    const std::vector<u64>& sequence,
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    for (size_t k = 0; k < sequence.size(); ++k) {
        ++checks;

        u64 rank;

        if (!rank_index(
                sequence[k],
                c,
                rank
            )) {
            ++failures;
            continue;
        }

        if (rank != static_cast<u64>(k)) {
            ++failures;
        }
    }

    return failures == 0;
}

/* ------------------------------------------------------------
   Verify:
       rank(next(j)) = rank(j) + 1
   ------------------------------------------------------------ */

bool verify_rank_successor(
    const std::vector<u64>& sequence,
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    for (size_t k = 0;
         k + 1 < sequence.size();
         ++k) {

        ++checks;

        u64 rank_a;
        u64 rank_b;

        if (!rank_index(
                sequence[k],
                c,
                rank_a
            ) ||
            !rank_index(
                sequence[k + 1],
                c,
                rank_b
            )) {
            ++failures;
            continue;
        }

        if (rank_a + 1 != rank_b) {
            ++failures;
        }
    }

    return failures == 0;
}

/* ------------------------------------------------------------
   Verify:
       unrank(rank(j)) = j
   ------------------------------------------------------------ */

bool verify_rank_unrank(
    const std::vector<u64>& sequence,
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    for (u64 j : sequence) {
        ++checks;

        u64 rank;
        u64 reconstructed;

        if (!rank_index(
                j,
                c,
                rank
            ) ||
            !unrank_index(
                rank,
                c,
                reconstructed
            )) {
            ++failures;
            continue;
        }

        if (reconstructed != j) {
            ++failures;
        }
    }

    return failures == 0;
}

/* ------------------------------------------------------------
   Verify:
       unrank(k) belongs to the Lucas-admissible set
       for every k.
   ------------------------------------------------------------ */

bool verify_all_ranks(
    const CaseData& c,
    u64& checks,
    u64& failures
) {
    checks = 0;
    failures = 0;

    std::vector<u64> radix;
    u64 tuple_count;

    if (!radix_sizes(
            c,
            radix,
            tuple_count
        )) {
        ++failures;
        return false;
    }

    /*
        q itself occupies the final tuple rank, so admissible
        ranks are 0 .. tuple_count-2.
    */
    if (tuple_count == 0) {
        ++failures;
        return false;
    }

    const u64 last_rank =
        tuple_count - 1;

    for (u64 rank = 0;
         rank < last_rank;
         ++rank) {

        ++checks;

        u64 j;

        if (!unrank_index(
                rank,
                c,
                j
            )) {
            ++failures;
            continue;
        }

        if (j >= c.q ||
            !leq_p(
                j,
                c.q,
                c.p
            )) {
            ++failures;
        }
    }

    /*
        The excluded final tuple should unrank to q itself.
    */
    u64 q_reconstructed;

    if (!unrank_index(
            last_rank,
            c,
            q_reconstructed
        )) {
        ++failures;
    } else if (q_reconstructed != c.q) {
        ++failures;
    }

    return failures == 0;
}

/* ------------------------------------------------------------
   Main case checker
   ------------------------------------------------------------ */

bool check_case(
    const CaseData& c,
    bool print_details
) {
    const auto reference =
        reference_indices(c);

    bool construction_ok;

    const auto sequence =
        successor_sequence(
            c,
            construction_ok
        );

    u64 rank_checks;
    u64 rank_failures;

    const bool rank_sequence_pass =
        verify_rank_sequence(
            sequence,
            c,
            rank_checks,
            rank_failures
        );

    u64 successor_checks;
    u64 successor_failures;

    const bool rank_successor_pass =
        verify_rank_successor(
            sequence,
            c,
            successor_checks,
            successor_failures
        );

    u64 roundtrip_checks;
    u64 roundtrip_failures;

    const bool roundtrip_pass =
        verify_rank_unrank(
            sequence,
            c,
            roundtrip_checks,
            roundtrip_failures
        );

    u64 all_rank_checks;
    u64 all_rank_failures;

    const bool all_ranks_pass =
        verify_all_ranks(
            c,
            all_rank_checks,
            all_rank_failures
        );

    const bool exact_sequence =
        sequences_equal(
            reference,
            sequence
        );

    std::vector<u64> radix;
    u64 tuple_count;

    const bool radix_pass =
        radix_sizes(
            c,
            radix,
            tuple_count
        );

    u64 expected_count = 0;

    if (radix_pass &&
        tuple_count > 0) {

        expected_count =
            tuple_count - 1;
    }

    const bool count_pass =
        radix_pass &&
        reference.size() == expected_count;

    const bool final_pass =
        construction_ok &&
        exact_sequence &&
        rank_sequence_pass &&
        rank_successor_pass &&
        roundtrip_pass &&
        all_ranks_pass &&
        count_pass;

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
            << sequence.size()
            << '\n';

        std::cout
            << "construction_pass="
            << construction_ok
            << '\n';

        std::cout
            << "exact_sequence_pass="
            << exact_sequence
            << '\n';

        std::cout
            << "tuple_count="
            << tuple_count
            << '\n';

        std::cout
            << "expected_admissible_count="
            << expected_count
            << '\n';

        std::cout
            << "count_pass="
            << count_pass
            << '\n';

        std::cout
            << "rank_sequence_checks="
            << rank_checks
            << '\n';

        std::cout
            << "rank_sequence_failures="
            << rank_failures
            << '\n';

        std::cout
            << "rank_sequence_pass="
            << rank_sequence_pass
            << '\n';

        std::cout
            << "rank_successor_checks="
            << successor_checks
            << '\n';

        std::cout
            << "rank_successor_failures="
            << successor_failures
            << '\n';

        std::cout
            << "rank_successor_pass="
            << rank_successor_pass
            << '\n';

        std::cout
            << "rank_unrank_checks="
            << roundtrip_checks
            << '\n';

        std::cout
            << "rank_unrank_failures="
            << roundtrip_failures
            << '\n';

        std::cout
            << "rank_unrank_pass="
            << roundtrip_pass
            << '\n';

        std::cout
            << "all_rank_checks="
            << all_rank_checks
            << '\n';

        std::cout
            << "all_rank_failures="
            << all_rank_failures
            << '\n';

        std::cout
            << "all_ranks_pass="
            << all_ranks_pass
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << '\n';

        const size_t sample_count =
            std::min<size_t>(
                15,
                sequence.size()
            );

        for (size_t i = 0;
             i < sample_count;
             ++i) {

            u64 rank;
            rank_index(
                sequence[i],
                c,
                rank
            );

            std::cout
                << "sample["
                << i
                << "] j="
                << sequence[i]
                << " rank="
                << rank
                << '\n';
        }

        std::cout << '\n';
    }

    return final_pass;
}

/* ------------------------------------------------------------
   Random test suite
   ------------------------------------------------------------ */

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

    for (int i = 0;
         i < cases;
         ++i) {

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
            The reference enumeration is O(q).
        */
        if (q > 2000000) {
            ++skipped;
            continue;
        }

        if (!check_case(
                c,
                false
            )) {
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

/* ------------------------------------------------------------
   Main
   ------------------------------------------------------------ */

int main() {
    std::cout
        << "START EXPERIMENT 220\n\n";

    bool all_pass = true;

    /*
        Small binary.
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
            check_case(
                c,
                true
            ) &&
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
            check_case(
                c,
                true
            ) &&
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
            check_case(
                c,
                true
            ) &&
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
            check_case(
                c,
                true
            ) &&
            all_pass;
    }

    /*
        Dense binary.
    */
    {
        CaseData c;

        build_case(
            2,
            0,
            0,
            19,
            1048575,
            c
        );

        all_pass =
            check_case(
                c,
                false
            ) &&
            all_pass;
    }

    /*
        Base 5 dense case.
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
            check_case(
                c,
                false
            ) &&
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
            check_case(
                c,
                true
            ) &&
            all_pass;
    }

    std::mt19937_64 rng(
        220220220ULL
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
        << "R(j) = sum_i j_i * product_{t<i}(q_t+1)\n";

    std::cout
        << "rank(j_k) = k\n";

    std::cout
        << "rank(next(j)) = rank(j) + 1\n";

    std::cout
        << "unrank(rank(j)) = j\n";

    std::cout
        << "unrank(k) is Lucas-admissible for 0 <= k < I\n";

    std::cout
        << "unrank(I) = q (excluded terminal tuple)\n";

    std::cout
        << "I = product_i(q_i+1) - 1\n";

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << all_pass
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 220\n";

    return all_pass ? 0 : 1;
}
