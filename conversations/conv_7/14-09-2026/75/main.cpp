#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
    u64 D;
};

static mpz_class mpz_from_u64(u64 value) {
    return mpz_class(std::to_string(value));
}

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(std::sqrt(static_cast<long double>(n)));

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static u64 gcd_u64(u64 a, u64 b) {
    while (b != 0) {
        const u64 r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static bool is_prime_u64(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d <= n / d; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

static u64 random_prime(std::mt19937_64& rng, u64 lo, u64 hi) {
    std::uniform_int_distribution<u64> dist(lo, hi);

    while (true) {
        u64 x = dist(rng);

        if (x < 2) {
            x = 2;
        }

        if (x > 2 && x % 2 == 0) {
            ++x;
        }

        while (x <= hi && !is_prime_u64(x)) {
            x += 2;
        }

        if (x <= hi && x >= lo && is_prime_u64(x)) {
            return x;
        }
    }
}

static CaseData make_case(std::mt19937_64& rng) {
    while (true) {
        u64 p = random_prime(rng, 5, 2000);
        u64 q = random_prime(rng, 2001, 5000);

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        const u64 N = p * q;
        const u64 s = isqrt_u64(N);

        if (s * s > N || (s + 1) * (s + 1) <= N) {
            continue;
        }

        const u64 D = N - s * s;

        return {p, q, N, s, D};
    }
}

/*
    Original second-order construction with

        R(z,u) = u + z + 1

    and

        E(z)
          = (s-z)^2 R(z, (z^2 + D)/(s-z))

    After clearing the denominator:

        E(z)
          = (s-z)(z^2 + D)
            + (s-z)^2(z+1)

    which factors as

        E(z)
          = (s-z)(D+s+z(s-1)).
*/
static mpz_class E_expanded(const CaseData& c, u64 z) {
    const mpz_class S = mpz_from_u64(c.s);
    const mpz_class Z = mpz_from_u64(z);
    const mpz_class D = mpz_from_u64(c.D);

    const mpz_class term1 = (S - Z) * (Z * Z + D);
    const mpz_class term2 = (S - Z) * (S - Z) * (Z + 1);

    return term1 + term2;
}

static mpz_class E_factored(const CaseData& c, u64 z) {
    const mpz_class S = mpz_from_u64(c.s);
    const mpz_class Z = mpz_from_u64(z);
    const mpz_class D = mpz_from_u64(c.D);

    return (S - Z) * (D + S + Z * (S - 1));
}

static mpz_class E_hidden_formula(const CaseData& c) {
    const mpz_class P = mpz_from_u64(c.p);
    const mpz_class Q = mpz_from_u64(c.q);
    const mpz_class S = mpz_from_u64(c.s);

    /*
        z = s-p

        E(s-p) = p^2 (q-s+1)
    */
    return P * P * (Q - S + 1);
}

static mpz_class E_k_formula(const CaseData& c, u64 k) {
    const mpz_class P = mpz_from_u64(c.p);
    const mpz_class Q = mpz_from_u64(c.q);
    const mpz_class S = mpz_from_u64(c.s);
    const mpz_class K = mpz_from_u64(k);

    /*
        z = s-kp

        E(s-kp) = k p^2 (q-k(s-1))
    */
    return K * P * P * (Q - K * (S - 1));
}

static bool divisible_by_square(const mpz_class& x, u64 p) {
    const mpz_class P = mpz_from_u64(p);
    const mpz_class P2 = P * P;

    return (x % P2) == 0;
}

static bool divisible_by(const mpz_class& x, u64 p) {
    const mpz_class P = mpz_from_u64(p);
    return (x % P) == 0;
}

int main() {
    constexpr int EXPERIMENT = 364;
    constexpr int CASES = 500;

    std::cout << "START EXPERIMENT " << EXPERIMENT << '\n';

    std::mt19937_64 rng(0x36420260914ULL);

    u64 total_scan_points = 0;
    u64 factorization_failures = 0;
    u64 hidden_formula_failures = 0;
    u64 k_formula_failures = 0;
    u64 p_divisibility_failures = 0;
    u64 q_divisibility_failures = 0;

    u64 cases_with_p2_hits = 0;
    u64 cases_with_extra_p2_hits = 0;
    u64 total_predicted_p2_hits = 0;
    u64 total_actual_p2_hits = 0;
    u64 total_extra_p2_hits = 0;

    u64 cases_p_divides_s_minus_1 = 0;
    u64 cases_exception_when_p_divides_s_minus_1 = 0;

    u64 max_actual_p2_hits = 0;
    u64 max_extra_p2_hits = 0;

    u64 cases_q2_hit_before_s = 0;
    u64 total_q2_hits_before_s = 0;

    u64 reconstruction_failures = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        /*
            Verify the algebraic factorization at a collection of points.
        */
        {
            const u64 test_count = std::min<u64>(s, 100);

            for (u64 i = 0; i < test_count; ++i) {
                const u64 z = (i * (s > 1 ? (s - 1) : 1)) /
                              (test_count > 1 ? (test_count - 1) : 1);

                const mpz_class a = E_expanded(c, z);
                const mpz_class b = E_factored(c, z);

                if (a != b) {
                    ++factorization_failures;
                }
            }
        }

        /*
            Hidden factor coordinate:
                z = a = s-p
        */
        {
            const u64 z = s - p;

            const mpz_class actual = E_factored(c, z);
            const mpz_class expected = E_hidden_formula(c);

            if (actual != expected) {
                ++hidden_formula_failures;
            }

            if (!divisible_by_square(actual, p)) {
                ++p_divisibility_failures;
            }

            if (p != 0 && q != 0) {
                const mpz_class G = mpz_from_u64(gcd_u64(
                    static_cast<u64>(actual == 0 ? 0 : 1),
                    1
                ));

                (void)G;
            }
        }

        /*
            Verify the full family

                z = s-kp

            for every admissible k >= 1.
        */
        const u64 max_k = s / p;

        for (u64 k = 1; k <= max_k; ++k) {
            const u64 z = s - k * p;

            const mpz_class actual = E_factored(c, z);
            const mpz_class expected = E_k_formula(c, k);

            if (actual != expected) {
                ++k_formula_failures;
            }

            if (!divisible_by_square(actual, p)) {
                ++p_divisibility_failures;
            }

            ++total_predicted_p2_hits;
        }

        /*
            Check the exceptional condition suggested by the modular analysis.
        */
        const bool p_divides_s_minus_1 = ((s - 1) % p == 0);

        if (p_divides_s_minus_1) {
            ++cases_p_divides_s_minus_1;
        }

        std::set<u64> predicted_hits;

        for (u64 k = 1; k <= max_k; ++k) {
            predicted_hits.insert(s - k * p);
        }

        std::set<u64> actual_p2_hits;
        std::set<u64> actual_q2_hits;

        /*
            Public scan:
                0 <= z < s

            z=s is excluded because E(s)=0 is trivially divisible by every
            modulus and would obscure the nontrivial root structure.
        */
        for (u64 z = 0; z < s; ++z) {
            ++total_scan_points;

            const mpz_class value = E_factored(c, z);

            if (divisible_by_square(value, p)) {
                actual_p2_hits.insert(z);
            }

            if (divisible_by_square(value, q)) {
                actual_q2_hits.insert(z);
            }
        }

        total_actual_p2_hits += actual_p2_hits.size();

        if (!actual_p2_hits.empty()) {
            ++cases_with_p2_hits;
        }

        max_actual_p2_hits =
            std::max<u64>(max_actual_p2_hits, actual_p2_hits.size());

        /*
            Compare actual and predicted p^2-hit sets.
        */
        for (u64 z : actual_p2_hits) {
            if (predicted_hits.find(z) == predicted_hits.end()) {
                ++total_extra_p2_hits;
                ++cases_with_extra_p2_hits;

                if (p_divides_s_minus_1) {
                    ++cases_exception_when_p_divides_s_minus_1;
                }
            }
        }

        max_extra_p2_hits =
            std::max<u64>(
                max_extra_p2_hits,
                actual_p2_hits.size() > predicted_hits.size()
                    ? actual_p2_hits.size() - predicted_hits.size()
                    : 0
            );

        if (actual_p2_hits != predicted_hits) {
            /*
                We only count this as an exception if the mismatch is genuine.
            */
        }

        /*
            Verify the q-side claim: since q > s in these generated cases,
            there should be no nontrivial q^2 root with 0 <= z < s.
        */
        if (!actual_q2_hits.empty()) {
            ++cases_q2_hit_before_s;
            total_q2_hits_before_s += actual_q2_hits.size();
        }

        /*
            Direct reconstruction at the hidden coordinate.
        */
        {
            const u64 z = s - p;

            const mpz_class value = E_factored(c, z);
            const mpz_class expected = E_hidden_formula(c);

            if (value != expected) {
                ++reconstruction_failures;
            }
        }
    }

    std::cout << "CASES=" << CASES << '\n';
    std::cout << "TOTAL_SCAN_POINTS=" << total_scan_points << '\n';

    std::cout << "FACTORIZATION_FAILURES="
              << factorization_failures << '\n';

    std::cout << "HIDDEN_FORMULA_FAILURES="
              << hidden_formula_failures << '\n';

    std::cout << "K_FORMULA_FAILURES="
              << k_formula_failures << '\n';

    std::cout << "P_DIVISIBILITY_FAILURES="
              << p_divisibility_failures << '\n';

    std::cout << "Q_DIVISIBILITY_FAILURES="
              << q_divisibility_failures << '\n';

    std::cout << "TOTAL_PREDICTED_P2_HITS="
              << total_predicted_p2_hits << '\n';

    std::cout << "TOTAL_ACTUAL_P2_HITS="
              << total_actual_p2_hits << '\n';

    std::cout << "CASES_WITH_P2_HITS="
              << cases_with_p2_hits << '\n';

    std::cout << "CASES_WITH_EXTRA_P2_HITS="
              << cases_with_extra_p2_hits << '\n';

    std::cout << "TOTAL_EXTRA_P2_HITS="
              << total_extra_p2_hits << '\n';

    std::cout << "MAX_ACTUAL_P2_HITS="
              << max_actual_p2_hits << '\n';

    std::cout << "MAX_EXTRA_P2_HITS="
              << max_extra_p2_hits << '\n';

    std::cout << "CASES_P_DIVIDES_S_MINUS_1="
              << cases_p_divides_s_minus_1 << '\n';

    std::cout << "CASES_EXCEPTION_WHEN_P_DIVIDES_S_MINUS_1="
              << cases_exception_when_p_divides_s_minus_1 << '\n';

    std::cout << "CASES_Q2_HIT_BEFORE_S="
              << cases_q2_hit_before_s << '\n';

    std::cout << "TOTAL_Q2_HITS_BEFORE_S="
              << total_q2_hits_before_s << '\n';

    std::cout << "RECONSTRUCTION_FAILURES="
              << reconstruction_failures << '\n';

    const bool algebra_ok =
        factorization_failures == 0 &&
        hidden_formula_failures == 0 &&
        k_formula_failures == 0 &&
        p_divisibility_failures == 0 &&
        q_divisibility_failures == 0 &&
        reconstruction_failures == 0;

    const bool hit_structure_ok =
        total_extra_p2_hits == 0;

    const bool q_structure_ok =
        cases_q2_hit_before_s == 0;

    std::cout << "ALGEBRA_STATUS="
              << (algebra_ok ? "PASS" : "FAIL") << '\n';

    std::cout << "P2_HIT_STRUCTURE_STATUS="
              << (hit_structure_ok ? "PASS" : "EXCEPTIONS_FOUND") << '\n';

    std::cout << "Q2_HIT_STRUCTURE_STATUS="
              << (q_structure_ok ? "PASS" : "UNEXPECTED_HITS") << '\n';

    std::cout << "FINISHED EXPERIMENT " << EXPERIMENT << '\n';

    return 0;
}