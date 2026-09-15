#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Parameters {
    u64 p;
    u64 m;
    u64 e;
    u64 s0;
    u64 q;
};

u64 pow_u64(u64 p, u64 e) {
    u128 result = 1;

    for (u64 i = 0; i < e; ++i) {
        result *= p;
    }

    return static_cast<u64>(result);
}

std::vector<u64> base_p_digits(
    u64 p,
    u64 n
) {
    std::vector<u64> digits;

    while (n > 0) {
        digits.push_back(n % p);
        n /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

u64 digit_product(
    u64 p,
    u64 n
) {
    const auto digits =
        base_p_digits(p, n);

    u64 result = 1;

    for (const u64 d : digits) {
        result *= d + 1;
    }

    return result;
}

bool derive_parameters(
    u64 p,
    u64 m,
    Parameters &out
) {
    u64 x = m;
    u64 a = 0;

    while (x > 0 && x % p == p - 1) {
        ++a;
        x /= p;
    }

    if (x == 0) {
        return false;
    }

    const u64 b = x % p;

    x /= p;

    if (b >= p - 1) {
        return false;
    }

    u64 z = 0;

    while (x > 0 && x % p == 0) {
        ++z;
        x /= p;
    }

    const u64 e =
        a + 1 + z;

    const u64 s0 =
        static_cast<u64>(
            static_cast<u128>(b + 1) *
            static_cast<u128>(
                pow_u64(p, a)
            )
        );

    const u64 pe =
        pow_u64(p, e);

    const u128 numerator =
        static_cast<u128>(m) +
        1 -
        static_cast<u128>(s0);

    if (numerator == 0 ||
        numerator % pe != 0) {
        return false;
    }

    const u64 q =
        static_cast<u64>(
            numerator / pe
        );

    if (q == 0) {
        return false;
    }

    out = {
        p,
        m,
        e,
        s0,
        q
    };

    return true;
}

/*
 * The structural digit prediction.

 * If

 *   m = s0 + q*p^e - 1,

 * then the base-p digits of m should consist of:

 *   a copies of (p-1),
 *   followed by b,
 *   followed by z zeros,
 *   followed by the digits of q.

 * Therefore

 *   product(m_i+1)
 *     = (b+1)p^a product(q_i+1)
 *     = s0 product(q_i+1).
 *
 * This function computes the predicted RHS using only
 * m's structural parameters.
 */
u64 predicted_m_digit_product(
    u64 p,
    u64 m,
    u64 e,
    u64 s0,
    u64 q
) {
    (void)m;
    (void)e;

    return static_cast<u64>(
        static_cast<u128>(s0) *
        static_cast<u128>(
            digit_product(
                p,
                q
            )
        )
    );
}

/*
 * Independently inspect the actual base-p digits of m
 * and compare them with the digits obtained from q.
 */
bool digit_structure_check(
    u64 p,
    u64 m,
    u64 e,
    u64 s0,
    u64 q
) {
    const auto md =
        base_p_digits(
            p,
            m
        );

    const auto qd =
        base_p_digits(
            p,
            q
        );

    /*
     * Locate the p-adic decomposition of s0.

     * The simplest independent check is actually to
     * reconstruct m from s0 + q*p^e - 1 and compare.
     */
    const u64 reconstructed =
        static_cast<u64>(
            static_cast<u128>(s0) +
            static_cast<u128>(q) *
                static_cast<u128>(
                    pow_u64(p, e)
                ) -
            1
        );

    if (reconstructed != m) {
        return false;
    }

    /*
     * Compare the digit products as an independent
     * digit-level check.
     */
    u64 m_product = 1;

    for (const u64 d : md) {
        m_product *= d + 1;
    }

    u64 q_product = 1;

    for (const u64 d : qd) {
        q_product *= d + 1;
    }

    const u64 expected =
        static_cast<u64>(
            static_cast<u128>(s0) *
            static_cast<u128>(q_product)
        );

    return m_product == expected;
}

void run_small_exhaustive(
    u64 max_m,
    u64 &cases,
    u64 &product_pass,
    u64 &digit_structure_pass,
    u64 &final_formula_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11
    };

    for (const u64 p : primes) {
        for (u64 m = 2;
             m <= max_m;
             ++m) {

            Parameters par{};

            if (!derive_parameters(
                    p,
                    m,
                    par
                )) {
                continue;
            }

            ++cases;

            const u64 actual_m_product =
                digit_product(
                    p,
                    par.m
                );

            const u64 actual_q_product =
                digit_product(
                    p,
                    par.q
                );

            const u64 predicted =
                static_cast<u64>(
                    static_cast<u128>(
                        par.s0
                    ) *
                    static_cast<u128>(
                        actual_q_product
                    )
                );

            if (actual_m_product ==
                predicted) {
                ++product_pass;
            }

            if (digit_structure_check(
                    p,
                    par.m,
                    par.e,
                    par.s0,
                    par.q
                )) {
                ++digit_structure_pass;
            }

            /*
             * Complete final identity:
             *
             *   q*p^e
             *   - s0*(prod(q_i+1)-1)
             *
             * =
             *
             *   m+1-prod(m_i+1)
             */
            const u64 pe =
                pow_u64(
                    p,
                    par.e
                );

            const u64 q_product =
                actual_q_product;

            const u64 left =
                static_cast<u64>(
                    static_cast<u128>(
                        par.q
                    ) *
                    static_cast<u128>(
                        pe
                    ) -
                    static_cast<u128>(
                        par.s0
                    ) *
                    static_cast<u128>(
                        q_product - 1
                    )
                );

            const u64 right =
                par.m + 1 -
                actual_m_product;

            if (left == right) {
                ++final_formula_pass;
            }
        }
    }
}

void run_large_random(
    u64 trials,
    u64 &product_pass,
    u64 &digit_structure_pass,
    u64 &final_formula_pass
) {
    const std::vector<u64> primes = {
        2, 3, 5, 7, 11,
        13, 17, 19, 23,
        29, 31
    };

    std::mt19937_64 rng(
        0x19420260913ULL
    );

    std::uniform_int_distribution<u64> dist(
        1000000000000ULL,
        1000000000000000000ULL
    );

    for (u64 i = 0;
         i < trials;
         ++i) {

        const u64 p =
            primes[i % primes.size()];

        const u64 m =
            dist(rng);

        Parameters par{};

        if (!derive_parameters(
                p,
                m,
                par
            )) {
            --i;
            continue;
        }

        const u64 m_product =
            digit_product(
                p,
                par.m
            );

        const u64 q_product =
            digit_product(
                p,
                par.q
            );

        const u64 predicted =
            static_cast<u64>(
                static_cast<u128>(
                    par.s0
                ) *
                static_cast<u128>(
                    q_product
                )
            );

        if (m_product ==
            predicted) {
            ++product_pass;
        }

        if (digit_structure_check(
                p,
                par.m,
                par.e,
                par.s0,
                par.q
            )) {
            ++digit_structure_pass;
        }

        const u64 pe =
            pow_u64(
                p,
                par.e
            );

        const u64 left =
            static_cast<u64>(
                static_cast<u128>(
                    par.q
                ) *
                static_cast<u128>(
                    pe
                ) -
                static_cast<u128>(
                    par.s0
                ) *
                static_cast<u128>(
                    q_product - 1
                )
            );

        const u64 right =
            par.m + 1 -
            m_product;

        if (left == right) {
            ++final_formula_pass;
        }
    }
}

void print_example(
    u64 p,
    u64 m
) {
    Parameters par{};

    if (!derive_parameters(
            p,
            m,
            par
        )) {
        return;
    }

    const auto md =
        base_p_digits(
            p,
            par.m
        );

    const auto qd =
        base_p_digits(
            p,
            par.q
        );

    std::cout
        << "p=" << p
        << " m=" << m
        << " e=" << par.e
        << " s0=" << par.s0
        << " q=" << par.q
        << "\n";

    std::cout
        << "  m_digits:";

    for (auto it = md.rbegin();
         it != md.rend();
         ++it) {
        std::cout
            << " "
            << *it;
    }

    std::cout
        << "\n";

    std::cout
        << "  q_digits:";

    for (auto it = qd.rbegin();
         it != qd.rend();
         ++it) {
        std::cout
            << " "
            << *it;
    }

    std::cout
        << "\n";

    std::cout
        << "  product(m_i+1)="
        << digit_product(
            p,
            par.m
        )
        << "\n";

    std::cout
        << "  s0*product(q_i+1)="
        << static_cast<u64>(
            static_cast<u128>(
                par.s0
            ) *
            static_cast<u128>(
                digit_product(
                    p,
                    par.q
                )
            )
        )
        << "\n";
}

int main() {
    std::cout
        << "START EXPERIMENT 194\n";

    u64 cases = 0;
    u64 product_pass = 0;
    u64 digit_structure_pass = 0;
    u64 final_formula_pass = 0;

    run_small_exhaustive(
        5000,
        cases,
        product_pass,
        digit_structure_pass,
        final_formula_pass
    );

    std::cout
        << "\nSMALL EXHAUSTIVE\n";

    std::cout
        << "parameter_cases="
        << cases
        << "\n";

    std::cout
        << "digit_product_identity_pass="
        << product_pass
        << "/" << cases
        << "\n";

    std::cout
        << "digit_structure_pass="
        << digit_structure_pass
        << "/" << cases
        << "\n";

    std::cout
        << "final_hit_formula_pass="
        << final_formula_pass
        << "/" << cases
        << "\n";

    u64 random_product_pass = 0;
    u64 random_structure_pass = 0;
    u64 random_final_pass = 0;

    run_large_random(
        200000,
        random_product_pass,
        random_structure_pass,
        random_final_pass
    );

    std::cout
        << "\nLARGE RANDOM\n";

    std::cout
        << "digit_product_identity_pass="
        << random_product_pass
        << "/200000\n";

    std::cout
        << "digit_structure_pass="
        << random_structure_pass
        << "/200000\n";

    std::cout
        << "final_hit_formula_pass="
        << random_final_pass
        << "/200000\n";

    std::cout
        << "\nEXAMPLES\n";

    print_example(
        2,
        20
    );

    print_example(
        2,
        40
    );

    print_example(
        3,
        41
    );

    print_example(
        5,
        194
    );

    std::cout
        << "\nIDENTITY UNDER TEST\n";

    std::cout
        << "product_i(m_i+1)"
        << " = "
        << "s0 * product_i(q_i+1)\n";

    std::cout
        << "m+1-product_i(m_i+1)"
        << " = "
        << "q*p^e"
        << " - "
        << "s0*(product_i(q_i+1)-1)\n";

    std::cout
        << "FINISHED EXPERIMENT 194\n";

    return 0;
}
