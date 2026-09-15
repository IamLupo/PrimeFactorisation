#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct CaseData {
    u64 p;
    u64 q;

    std::vector<u64> digits;
    std::vector<u64> radix;
    std::vector<u64> powers;
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

        c.digits.push_back(d);
        c.radix.push_back(d + 1);

        temp /= p;

        if (temp == 0) {
            break;
        }
    }

    c.powers.resize(c.digits.size());
    c.powers[0] = 1;

    for (size_t i = 1; i < c.powers.size(); ++i) {
        if (!safe_mul(
                c.powers[i - 1],
                p,
                c.powers[i]
            )) {
            return false;
        }
    }

    return true;
}

/*
    S(q) =
        sum_r C_r (p^r - q mod p^r)
*/
bool compute_S(
    const CaseData& c,
    u64& S
) {
    S = 0;

    for (u64 r = 0;
         r < c.digits.size();
         ++r) {

        u64 count = c.digits[r];

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

        const u64 pr =
            c.powers[r];

        const u64 remainder =
            c.q % pr;

        const u64 factor =
            pr - remainder;

        u64 term;

        if (!safe_mul(
                count,
                factor,
                term
            )) {
            return false;
        }

        if (!safe_add(
                S,
                term,
                S
            )) {
            return false;
        }
    }

    return true;
}

/*
    q = Q + a p^R
*/
bool split_high_digit(
    const CaseData& c,
    u64& Q,
    u64& a,
    u64& R,
    u64& pR
) {
    if (c.digits.empty()) {
        return false;
    }

    R =
        static_cast<u64>(
            c.digits.size() - 1
        );

    a = c.digits[R];
    pR = c.powers[R];

    u64 high;

    if (!safe_mul(
            a,
            pR,
            high
        )) {
        return false;
    }

    if (high > c.q) {
        return false;
    }

    Q = c.q - high;

    return true;
}

/*
    Compute S(Q), including the base case Q=0.

        S(0) = 0.
*/
bool compute_S_value(
    u64 p,
    u64 q,
    u64& S
) {
    if (q == 0) {
        S = 0;
        return true;
    }

    CaseData c;

    if (!build_case(
            p,
            q,
            c
        )) {
        return false;
    }

    return compute_S(
        c,
        S
    );
}

/*
    Verify:

        S(q)
        =
        (a+1)S(Q)
        +
        a(p^R-Q)
*/
bool check_case(
    const CaseData& c,
    bool print_details
) {
    u64 Q;
    u64 a;
    u64 R;
    u64 pR;

    if (!split_high_digit(
            c,
            Q,
            a,
            R,
            pR
        )) {
        return false;
    }

    u64 S_q;

    if (!compute_S(
            c,
            S_q
        )) {
        return false;
    }

    u64 S_Q;

    if (!compute_S_value(
            c.p,
            Q,
            S_Q
        )) {
        return false;
    }

    u64 a_plus_one;

    if (!safe_add(
            a,
            1,
            a_plus_one
        )) {
        return false;
    }

    u64 first_term;

    if (!safe_mul(
            a_plus_one,
            S_Q,
            first_term
        )) {
        return false;
    }

    if (pR < Q) {
        return false;
    }

    const u64 second_factor =
        pR - Q;

    u64 second_term;

    if (!safe_mul(
            a,
            second_factor,
            second_term
        )) {
        return false;
    }

    u64 recurrence_value;

    if (!safe_add(
            first_term,
            second_term,
            recurrence_value
        )) {
        return false;
    }

    u64 reconstructed_q;

    if (!safe_mul(
            a,
            pR,
            reconstructed_q
        )) {
        return false;
    }

    if (!safe_add(
            Q,
            reconstructed_q,
            reconstructed_q
        )) {
        return false;
    }

    const bool decomposition_pass =
        reconstructed_q == c.q;

    const bool direct_identity_pass =
        S_q == c.q;

    const bool recurrence_pass =
        S_q == recurrence_value;

    /*
        The one-digit case is now explicitly checked.
    */
    bool base_case_pass = true;

    if (c.digits.size() == 1) {
        base_case_pass =
            Q == 0 &&
            R == 0 &&
            pR == 1 &&
            S_Q == 0 &&
            S_q == a &&
            recurrence_value == a;
    }

    const bool final_pass =
        decomposition_pass &&
        direct_identity_pass &&
        recurrence_pass &&
        base_case_pass;

    if (print_details) {
        std::cout
            << "CASE\n";

        std::cout
            << "p=" << c.p
            << " q=" << c.q
            << '\n';

        std::cout
            << "digit_count="
            << c.digits.size()
            << '\n';

        std::cout
            << "highest_digit_position="
            << R
            << '\n';

        std::cout
            << "highest_digit_a="
            << a
            << '\n';

        std::cout
            << "p^R="
            << pR
            << '\n';

        std::cout
            << "Q="
            << Q
            << '\n';

        std::cout
            << "S(Q)="
            << S_Q
            << '\n';

        std::cout
            << "S(q)="
            << S_q
            << '\n';

        std::cout
            << "recurrence_value="
            << recurrence_value
            << '\n';

        std::cout
            << "expected_S_q="
            << c.q
            << '\n';

        std::cout
            << "q_decomposition_pass="
            << decomposition_pass
            << '\n';

        std::cout
            << "direct_S_q_eq_q="
            << direct_identity_pass
            << '\n';

        std::cout
            << "recurrence_pass="
            << recurrence_pass
            << '\n';

        std::cout
            << "base_case_pass="
            << base_case_pass
            << '\n';

        std::cout
            << "final_pass="
            << final_pass
            << "\n\n";
    }

    return final_pass;
}

/*
    Explicit recursion chain:

        q
        -> Q
        -> next Q
        -> ...
        -> 0
*/
bool check_recursive_chain(
    u64 p,
    u64 initial_q,
    u64& chain_steps
) {
    chain_steps = 0;

    u64 q = initial_q;

    while (true) {
        if (q == 0) {
            /*
                True base case:
                    S(0)=0.
            */
            u64 S0;

            if (!compute_S_value(
                    p,
                    0,
                    S0
                )) {
                return false;
            }

            return S0 == 0;
        }

        CaseData c;

        if (!build_case(
                p,
                q,
                c
            )) {
            return false;
        }

        if (!check_case(
                c,
                false
            )) {
            return false;
        }

        ++chain_steps;

        if (q < p) {
            /*
                The next Q after a single-digit number is 0.
            */
            q = 0;
            continue;
        }

        const u64 R =
            static_cast<u64>(
                c.digits.size() - 1
            );

        const u64 a =
            c.digits[R];

        const u64 pR =
            c.powers[R];

        u64 high;

        if (!safe_mul(
                a,
                pR,
                high
            )) {
            return false;
        }

        if (high > q) {
            return false;
        }

        q -= high;
    }
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
        << "START EXPERIMENT 230\n\n";

    std::mt19937_64 rng(
        230230230ULL
    );

    u64 failed_cases = 0;
    u64 total_cases = 0;

    auto execute =
        [&](u64 p,
            u64 q,
            bool print_details) {

        CaseData c;

        if (!build_case(
                p,
                q,
                c
            )) {
            ++failed_cases;
            ++total_cases;
            return;
        }

        ++total_cases;

        if (!check_case(
                c,
                print_details
            )) {
            ++failed_cases;
        }
    };

    /*
        Deterministic checks.
    */
    execute(2, 63, true);
    execute(3, 80, true);
    execute(5, 100, true);
    execute(2, 1048575, true);
    execute(3, 987654321, true);
    execute(5, 1000007654321ULL, true);
    execute(2, 1073741825ULL, true);

    /*
        Explicit Q=0 base.
    */
    {
        u64 S0;

        const bool base_pass =
            compute_S_value(
                5,
                0,
                S0
            ) &&
            S0 == 0;

        std::cout
            << "BASE_CASE"
            << " p=5"
            << " q=0"
            << " S(0)="
            << S0
            << " pass="
            << base_pass
            << '\n';

        if (!base_pass) {
            ++failed_cases;
        }

        ++total_cases;
    }

    /*
        Random direct cases.
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

        const u64 q =
            random_bounded(
                rng,
                1,
                5000000000000ULL
            );

        CaseData c;

        if (!build_case(
                p,
                q,
                c
            )) {
            ++failed_cases;
            ++total_cases;
            continue;
        }

        ++total_cases;

        if (!check_case(
                c,
                false
            )) {
            ++failed_cases;

            if (failed_cases <= 5) {
                std::cout
                    << "RANDOM_FAILURE"
                    << " p=" << p
                    << " q=" << q
                    << '\n';
            }
        }
    }

    /*
        Recursive chains.
    */
    u64 chain_cases = 0;
    u64 chain_failures = 0;
    u64 chain_steps = 0;

    for (u64 p : {2ULL, 3ULL, 5ULL}) {
        for (int sample = 0;
             sample < 1000;
             ++sample) {

            const u64 q =
                random_bounded(
                    rng,
                    1,
                    1000000000000ULL
                );

            ++chain_cases;

            u64 steps;

            if (!check_recursive_chain(
                    p,
                    q,
                    steps
                )) {

                ++chain_failures;

                if (chain_failures <= 5) {
                    std::cout
                        << "CHAIN_FAILURE"
                        << " p=" << p
                        << " q=" << q
                        << '\n';
                }

            } else {
                chain_steps += steps;
            }
        }
    }

    std::cout
        << "RECURSIVE CHAIN\n";

    std::cout
        << "chain_cases="
        << chain_cases
        << '\n';

    std::cout
        << "chain_failures="
        << chain_failures
        << '\n';

    std::cout
        << "total_chain_steps="
        << chain_steps
        << '\n';

    std::cout
        << "chain_pass="
        << (chain_failures == 0)
        << '\n';

    std::cout << '\n';

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
        << "OVERALL_PASS="
        << (
            failed_cases == 0 &&
            chain_failures == 0
        )
        << '\n';

    std::cout << '\n';

    std::cout
        << "THEOREM TARGET\n";

    std::cout
        << "S(0)=0\n";

    std::cout
        << "q=Q+a p^R\n";

    std::cout
        << "S(q)=(a+1)S(Q)+a(p^R-Q)\n";

    std::cout
        << "S(q)=q\n";

    std::cout
        << "recursive removal reaches S(0)=0\n";

    std::cout
        << "FINISHED EXPERIMENT 230\n";

    return (
        failed_cases == 0 &&
        chain_failures == 0
    ) ? 0 : 1;
}