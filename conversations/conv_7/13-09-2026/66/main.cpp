#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct CaseData {
    u64 p;
    u64 q;

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

u64 digit_at(u64 x, u64 p, u64 i) {
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
    u64 q,
    CaseData& c
) {
    if (p < 2 || q == 0) {
        return false;
    }

    c = {};
    c.p = p;
    c.q = q;

    u64 temp = q;

    while (true) {
        const u64 d = temp % p;

        c.q_digits.push_back(d);

        u64 radix_value;

        if (!safe_add(d, 1, radix_value)) {
            return false;
        }

        c.radix.push_back(radix_value);

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

    if (c.tuple_count == 0) {
        return false;
    }

    c.interval_count = c.tuple_count - 1;

    return true;
}

/*
    IMPORTANT:

    Base-p digit extraction is

        j_i = floor(j / p^i) mod p

    NOT modulo (q_i + 1).

    The latter is the mixed-radix radix for the rank,
    not the base-p representation of j.
*/
u64 p_digit(
    u64 j,
    const CaseData& c,
    size_t i
) {
    return (j / c.powers[i]) % c.p;
}

/*
    Rank:

        R(j) = sum_i j_i W_i

    where

        W_i = product_{t<i}(q_t+1).
*/
bool rank_index(
    u64 j,
    const CaseData& c,
    u64& rank
) {
    rank = 0;

    for (size_t i = 0; i < c.radix.size(); ++i) {
        const u64 d =
            p_digit(
                j,
                c,
                i
            );

        const u64 qd =
            c.q_digits[i];

        if (d > qd) {
            return false;
        }

        u64 term;

        if (!safe_mul(
                d,
                c.weights[i],
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
    }

    return true;
}

/*
    Unrank:

        j_i =
            floor(k / W_i) mod (q_i+1)

    and

        j = sum_i j_i p^i.
*/
bool unrank_index(
    u64 rank,
    const CaseData& c,
    u64& j
) {
    if (rank >= c.tuple_count) {
        return false;
    }

    j = 0;

    for (size_t i = 0; i < c.radix.size(); ++i) {
        const u64 d =
            (rank / c.weights[i]) %
            c.radix[i];

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

bool check_case(
    const CaseData& c,
    std::mt19937_64& rng,
    u64 random_checks,
    bool print_details
) {
    u64 failures = 0;

    u64 terminal = 0;

    const bool terminal_ok =
        unrank_index(
            c.interval_count,
            c,
            terminal
        ) &&
        terminal == c.q;

    if (!terminal_ok) {
        ++failures;
    }

    u64 q_rank = 0;

    const bool q_rank_ok =
        rank_index(
            c.q,
            c,
            q_rank
        ) &&
        q_rank == c.interval_count;

    if (!q_rank_ok) {
        ++failures;
    }

    for (u64 n = 0;
         n < random_checks;
         ++n) {

        const u64 k =
            (c.interval_count == 1)
                ? 0
                : rng() % c.interval_count;

        u64 j;

        if (!unrank_index(
                k,
                c,
                j
            )) {
            ++failures;
            continue;
        }

        /*
            Lucas digitwise admissibility.
        */
        if (j >= c.q ||
            !leq_p(
                j,
                c.q,
                c.p
            )) {
            ++failures;
        }

        /*
            Every actual base-p digit must satisfy
                j_i <= q_i.
        */
        for (size_t i = 0;
             i < c.q_digits.size();
             ++i) {

            const u64 jd =
                p_digit(
                    j,
                    c,
                    i
                );

            if (jd > c.q_digits[i]) {
                ++failures;
            }
        }

        u64 recovered_rank;

        if (!rank_index(
                j,
                c,
                recovered_rank
            )) {
            ++failures;
            continue;
        }

        if (recovered_rank != k) {
            ++failures;
        }

        u64 roundtrip;

        if (!unrank_index(
                recovered_rank,
                c,
                roundtrip
            )) {
            ++failures;
            continue;
        }

        if (roundtrip != j) {
            ++failures;
        }
    }

    /*
        Check first ranks.
    */
    const u64 low_checks =
        std::min<u64>(
            16,
            c.interval_count
        );

    for (u64 k = 0;
         k < low_checks;
         ++k) {

        u64 j;
        u64 rank;

        if (!unrank_index(
                k,
                c,
                j
            ) ||
            !rank_index(
                j,
                c,
                rank
            )) {
            ++failures;
            continue;
        }

        if (rank != k) {
            ++failures;
        }
    }

    /*
        Check last ranks.
    */
    const u64 high_checks =
        std::min<u64>(
            16,
            c.interval_count
        );

    for (u64 offset = 1;
         offset <= high_checks;
         ++offset) {

        const u64 k =
            c.interval_count - offset;

        u64 j;
        u64 rank;

        if (!unrank_index(
                k,
                c,
                j
            ) ||
            !rank_index(
                j,
                c,
                rank
            )) {
            ++failures;
            continue;
        }

        if (rank != k) {
            ++failures;
        }
    }

    const bool pass =
        failures == 0;

    if (print_details) {
        std::cout
            << "CASE\n";

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << '\n';

        std::cout
            << "digits="
            << c.q_digits.size()
            << '\n';

        std::cout
            << "tuple_count="
            << c.tuple_count
            << '\n';

        std::cout
            << "interval_count="
            << c.interval_count
            << '\n';

        std::cout
            << "terminal_unrank="
            << terminal
            << '\n';

        std::cout
            << "terminal_expected="
            << c.q
            << '\n';

        std::cout
            << "terminal_unrank_pass="
            << terminal_ok
            << '\n';

        std::cout
            << "rank_q="
            << q_rank
            << '\n';

        std::cout
            << "rank_q_expected="
            << c.interval_count
            << '\n';

        std::cout
            << "rank_q_pass="
            << q_rank_ok
            << '\n';

        std::cout
            << "random_checks="
            << random_checks
            << '\n';

        std::cout
            << "failures="
            << failures
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

            unrank_index(
                k,
                c,
                j
            );

            std::cout
                << "sample_low["
                << k
                << "]=rank "
                << k
                << " j "
                << j
                << '\n';
        }

        for (u64 offset = 1;
             offset <= 10 &&
             offset <= c.interval_count;
             ++offset) {

            const u64 k =
                c.interval_count - offset;

            u64 j;

            unrank_index(
                k,
                c,
                j
            );

            std::cout
                << "sample_high["
                << offset
                << "]=rank "
                << k
                << " j "
                << j
                << '\n';
        }

        std::cout << '\n';
    }

    return pass;
}

bool run_suite(
    std::mt19937_64& rng,
    u64 p,
    u64 q,
    u64 random_checks,
    bool print_case
) {
    CaseData c;

    if (!build_case(
            p,
            q,
            c
        )) {
        return false;
    }

    return check_case(
        c,
        rng,
        random_checks,
        print_case
    );
}

int main() {
    std::cout
        << "START EXPERIMENT 221\n\n";

    bool all_pass = true;

    std::mt19937_64 rng(
        221221221ULL
    );

    all_pass =
        run_suite(
            rng,
            2,
            63,
            1000,
            true
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            3,
            80,
            1000,
            true
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            5,
            100,
            1000,
            true
        ) &&
        all_pass;

    /*
        Large cases are tested without O(q) enumeration.
    */
    all_pass =
        run_suite(
            rng,
            2,
            (1ULL << 60) + 12345ULL,
            100000,
            true
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            3,
            1000001234567ULL,
            100000,
            true
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            5,
            1000007654321ULL,
            100000,
            true
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            2,
            UINT64_MAX - 1,
            100000,
            false
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            3,
            UINT64_MAX / 3,
            100000,
            false
        ) &&
        all_pass;

    all_pass =
        run_suite(
            rng,
            5,
            UINT64_MAX / 5,
            100000,
            false
        ) &&
        all_pass;

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "base-p digit j_i = floor(j/p^i) mod p\n";

    std::cout
        << "rank(j) = sum_i j_i W_i\n";

    std::cout
        << "W_i = product_{t<i}(q_t+1)\n";

    std::cout
        << "unrank(k)_i = floor(k/W_i) mod (q_i+1)\n";

    std::cout
        << "rank(unrank(k)) = k\n";

    std::cout
        << "unrank(rank(j)) = j\n";

    std::cout
        << "unrank(k) is Lucas-admissible for 0 <= k < I\n";

    std::cout
        << "unrank(I) = q\n";

    std::cout
        << "no O(q) enumeration is used\n";

    std::cout << '\n';

    std::cout
        << "OVERALL_PASS="
        << all_pass
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 221\n";

    return all_pass ? 0 : 1;
}