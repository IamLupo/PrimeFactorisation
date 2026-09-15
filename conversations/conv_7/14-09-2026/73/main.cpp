#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

struct Coeff {
    i64 alpha;
    i64 beta;
    i64 gamma;
    i64 delta;
    i64 epsilon;
    i64 zeta;
};

static mpz_class mpz_from_u64(u64 value) {
    return mpz_class(
        std::to_string(value)
    );
}

static mpz_class mpz_from_i64(i64 value) {
    return mpz_class(
        std::to_string(value)
    );
}

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        __builtin_sqrtl(
            static_cast<long double>(n)
        )
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(
        static_cast<std::size_t>(limit + 1),
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<long long>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 50000) {
            usable.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        usable.size() - 1
    );

    std::vector<CaseData> result;
    result.reserve(
        static_cast<std::size_t>(count)
    );

    while (
        static_cast<int>(result.size()) < count
    ) {
        const std::size_t i = dist(rng);
        const std::size_t j = dist(rng);

        if (i == j) {
            continue;
        }

        const u64 p =
            static_cast<u64>(
                std::min(
                    usable[i],
                    usable[j]
                )
            );

        const u64 q =
            static_cast<u64>(
                std::max(
                    usable[i],
                    usable[j]
                )
            );

        const u64 N = p * q;

        if (!used_N.insert(N).second) {
            continue;
        }

        result.push_back({
            p,
            q,
            N,
            isqrt_u64(N)
        });
    }

    return result;
}

/*
 * R(z,u) =
 *
 * alpha*u^2
 * + beta*z*u
 * + gamma*z^2
 * + delta*u
 * + epsilon*z
 * + zeta
 */
static mpz_class evaluate_R(
    const Coeff& c,
    const mpz_class& z,
    const mpz_class& u
) {
    const mpz_class alpha =
        mpz_from_i64(c.alpha);

    const mpz_class beta =
        mpz_from_i64(c.beta);

    const mpz_class gamma =
        mpz_from_i64(c.gamma);

    const mpz_class delta =
        mpz_from_i64(c.delta);

    const mpz_class epsilon =
        mpz_from_i64(c.epsilon);

    const mpz_class zeta =
        mpz_from_i64(c.zeta);

    return
        alpha * u * u
        +
        beta * z * u
        +
        gamma * z * z
        +
        delta * u
        +
        epsilon * z
        +
        zeta;
}

/*
 * E(z) =
 *
 * (s-z)^2 R(z,(z^2+D)/(s-z)).
 *
 * After clearing the denominator:
 *
 * alpha*(z^2+D)^2
 *
 * + beta*z*(z^2+D)*(s-z)
 *
 * + gamma*z^2*(s-z)^2
 *
 * + delta*(z^2+D)*(s-z)
 *
 * + epsilon*z*(s-z)^2
 *
 * + zeta*(s-z)^2
 */
static mpz_class evaluate_E(
    const Coeff& c,
    const mpz_class& z,
    const mpz_class& s,
    const mpz_class& D
) {
    const mpz_class alpha =
        mpz_from_i64(c.alpha);

    const mpz_class beta =
        mpz_from_i64(c.beta);

    const mpz_class gamma =
        mpz_from_i64(c.gamma);

    const mpz_class delta =
        mpz_from_i64(c.delta);

    const mpz_class epsilon =
        mpz_from_i64(c.epsilon);

    const mpz_class zeta =
        mpz_from_i64(c.zeta);

    const mpz_class z2 =
        z * z;

    const mpz_class s_minus_z =
        s - z;

    const mpz_class z2_plus_D =
        z2 + D;

    return
        alpha *
            z2_plus_D *
            z2_plus_D

        +
        beta *
            z *
            z2_plus_D *
            s_minus_z

        +
        gamma *
            z2 *
            s_minus_z *
            s_minus_z

        +
        delta *
            z2_plus_D *
            s_minus_z

        +
        epsilon *
            z *
            s_minus_z *
            s_minus_z

        +
        zeta *
            s_minus_z *
            s_minus_z;
}

static mpz_class gcd_mpz(
    const mpz_class& value,
    const mpz_class& n
) {
    mpz_class x = value;

    if (x < 0) {
        x = -x;
    }

    mpz_class g;

    mpz_gcd(
        g.get_mpz_t(),
        x.get_mpz_t(),
        n.get_mpz_t()
    );

    return g;
}

static Coeff random_coefficients(
    std::mt19937_64& rng
) {
    std::uniform_int_distribution<int> dist(
        -3,
        3
    );

    Coeff c{
        dist(rng),
        dist(rng),
        dist(rng),
        dist(rng),
        dist(rng),
        dist(rng)
    };

    /*
     * Avoid the zero polynomial.
     */
    if (
        c.alpha == 0 &&
        c.beta == 0 &&
        c.gamma == 0 &&
        c.delta == 0 &&
        c.epsilon == 0 &&
        c.zeta == 0
    ) {
        c.alpha = 1;
    }

    return c;
}

static void print_coeff(
    const Coeff& c
) {
    std::cout
        << "("
        << c.alpha << ","
        << c.beta << ","
        << c.gamma << ","
        << c.delta << ","
        << c.epsilon << ","
        << c.zeta
        << ")";
}

static bool mpz_divisible_by(
    const mpz_class& x,
    const mpz_class& d
) {
    return (x % d) == 0;
}

int main() {
    std::cout << "START EXPERIMENT 362\n";

    constexpr int CASE_COUNT = 500;
    constexpr int POLYNOMIALS_PER_CASE = 20;
    constexpr int PRIME_LIMIT = 50000;

    std::mt19937_64 rng(
        0x362362362ULL
    );

    const std::vector<int> primes =
        sieve_primes(
            PRIME_LIMIT
        );

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    u64 total_tests = 0;

    u64 p_squared_identity_failures = 0;
    u64 q_squared_identity_failures = 0;

    u64 p_squared_divisibility_failures = 0;
    u64 q_squared_divisibility_failures = 0;

    u64 exact_p2_gcd = 0;
    u64 exact_q2_gcd = 0;

    u64 exact_p_gcd = 0;
    u64 exact_q_gcd = 0;

    u64 unexpected_p_gcd = 0;
    u64 unexpected_q_gcd = 0;

    u64 p_valuation_at_least_2 = 0;
    u64 q_valuation_at_least_2 = 0;

    u64 p_valuation_exactly_2 = 0;
    u64 q_valuation_exactly_2 = 0;

    u64 p_valuation_at_least_3 = 0;
    u64 q_valuation_at_least_3 = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c =
            cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s_u64 = c.s;

        const u64 D_u64 =
            N - s_u64 * s_u64;

        const mpz_class pz =
            mpz_from_u64(p);

        const mpz_class qz =
            mpz_from_u64(q);

        const mpz_class Nz =
            mpz_from_u64(N);

        const mpz_class sz =
            mpz_from_u64(s_u64);

        const mpz_class Dz =
            mpz_from_u64(D_u64);

        /*
         * Hidden roots of the auxiliary quadratic:
         *
         * a = s-p
         *
         * -b = s-q
         */
        const u64 a =
            s_u64 - p;

        const i64 b =
            -static_cast<i64>(
                q - s_u64
            );

        const mpz_class az =
            mpz_from_u64(a);

        const mpz_class minus_bz =
            mpz_from_i64(b);

        /*
         * Recover
         *
         * u = (a^2+D)/(s-a).
         */
        const mpz_class numerator_u =
            az * az + Dz;

        const mpz_class denominator_u =
            sz - az;

        if (denominator_u == 0) {
            std::cerr
                << "INTERNAL ERROR: zero denominator\n";
            return 1;
        }

        const mpz_class remainder_u =
            numerator_u % denominator_u;

        if (remainder_u != 0) {
            std::cerr
                << "INTERNAL ERROR: non-integral u\n";
            return 1;
        }

        const mpz_class uz =
            numerator_u /
            denominator_u;

        const mpz_class expected_u =
            pz + qz -
            mpz_class(2) * sz;

        if (uz != expected_u) {
            std::cerr
                << "INTERNAL ERROR: u identity failed\n";
            return 1;
        }

        /*
         * N^2 is small enough for this experiment, but keep it
         * entirely in GMP so there is no fixed-width overflow.
         */
        const mpz_class N2 =
            Nz * Nz;

        const mpz_class p2 =
            pz * pz;

        const mpz_class q2 =
            qz * qz;

        const mpz_class p3 =
            p2 * pz;

        const mpz_class q3 =
            q2 * qz;

        for (int j = 0;
             j < POLYNOMIALS_PER_CASE;
             ++j) {

            ++total_tests;

            const Coeff coeff =
                random_coefficients(
                    rng
                );

            /*
             * ----------------------------------------------------
             * E(a)
             *
             * Since s-a=p:
             *
             * E(a)=p^2 R(a,u)
             * ----------------------------------------------------
             */
            const mpz_class Ea =
                evaluate_E(
                    coeff,
                    az,
                    sz,
                    Dz
                );

            const mpz_class R_a =
                evaluate_R(
                    coeff,
                    az,
                    uz
                );

            const mpz_class expected_Ea =
                p2 * R_a;

            if (Ea != expected_Ea) {
                ++p_squared_identity_failures;
            }

            /*
             * ----------------------------------------------------
             * E(-b)
             *
             * Since s-(-b)=q:
             *
             * E(-b)=q^2 R(-b,u)
             * ----------------------------------------------------
             */
            const mpz_class R_b =
                evaluate_R(
                    coeff,
                    minus_bz,
                    uz
                );

            const mpz_class Eb =
                evaluate_E(
                    coeff,
                    minus_bz,
                    sz,
                    Dz
                );

            const mpz_class expected_Eb =
                q2 * R_b;

            if (Eb != expected_Eb) {
                ++q_squared_identity_failures;
            }

            /*
             * Direct p^2/q^2 divisibility.
             */
            if (!mpz_divisible_by(Ea, p2)) {
                ++p_squared_divisibility_failures;
            }

            if (!mpz_divisible_by(Eb, q2)) {
                ++q_squared_divisibility_failures;
            }

            /*
             * gcd(E(a),N^2)
             */
            const mpz_class gp2 =
                gcd_mpz(Ea, N2);

            const mpz_class gq2 =
                gcd_mpz(Eb, N2);

            if (gp2 == p2) {
                ++exact_p2_gcd;
            }

            if (gq2 == q2) {
                ++exact_q2_gcd;
            }

            if (gp2 == pz) {
                ++exact_p_gcd;
            }

            if (gq2 == qz) {
                ++exact_q_gcd;
            }

            /*
             * Any gcd other than 1,p,p^2 is unexpected on the
             * p side; likewise for q.
             */
            if (
                gp2 != 1 &&
                gp2 != pz &&
                gp2 != p2
            ) {
                ++unexpected_p_gcd;
            }

            if (
                gq2 != 1 &&
                gq2 != qz &&
                gq2 != q2
            ) {
                ++unexpected_q_gcd;
            }

            /*
             * Valuation checks.
             */
            const bool p2_hit =
                mpz_divisible_by(
                    Ea,
                    p2
                );

            const bool q2_hit =
                mpz_divisible_by(
                    Eb,
                    q2
                );

            const bool p3_hit =
                mpz_divisible_by(
                    Ea,
                    p3
                );

            const bool q3_hit =
                mpz_divisible_by(
                    Eb,
                    q3
                );

            if (p2_hit) {
                ++p_valuation_at_least_2;
            }

            if (q2_hit) {
                ++q_valuation_at_least_2;
            }

            if (p2_hit && !p3_hit) {
                ++p_valuation_exactly_2;
            }

            if (q2_hit && !q3_hit) {
                ++q_valuation_exactly_2;
            }

            if (p3_hit) {
                ++p_valuation_at_least_3;
            }

            if (q3_hit) {
                ++q_valuation_at_least_3;
            }

            if (case_id == 0 && j < 5) {
                std::cout
                    << "\nSAMPLE="
                    << (j + 1)
                    << "\n";

                std::cout
                    << "COEFF=";
                print_coeff(coeff);
                std::cout
                    << "\n";

                std::cout
                    << "E(a)="
                    << Ea
                    << "\n";

                std::cout
                    << "E(-b)="
                    << Eb
                    << "\n";

                std::cout
                    << "GCD_EA_N2="
                    << gp2
                    << "\n";

                std::cout
                    << "GCD_EB_N2="
                    << gq2
                    << "\n";

                std::cout
                    << "P2_DIVISIBLE="
                    << (p2_hit ? "YES" : "NO")
                    << "\n";

                std::cout
                    << "Q2_DIVISIBLE="
                    << (q2_hit ? "YES" : "NO")
                    << "\n";
            }
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "POLYNOMIALS_PER_CASE="
        << POLYNOMIALS_PER_CASE
        << "\n";

    std::cout
        << "TOTAL_TESTS="
        << total_tests
        << "\n";

    std::cout
        << "P_SQUARED_IDENTITY_FAILURES="
        << p_squared_identity_failures
        << "\n";

    std::cout
        << "Q_SQUARED_IDENTITY_FAILURES="
        << q_squared_identity_failures
        << "\n";

    std::cout
        << "P_SQUARED_DIVISIBILITY_FAILURES="
        << p_squared_divisibility_failures
        << "\n";

    std::cout
        << "Q_SQUARED_DIVISIBILITY_FAILURES="
        << q_squared_divisibility_failures
        << "\n";

    std::cout
        << "EXACT_P2_GCD="
        << exact_p2_gcd
        << "\n";

    std::cout
        << "EXACT_Q2_GCD="
        << exact_q2_gcd
        << "\n";

    std::cout
        << "EXACT_P_GCD="
        << exact_p_gcd
        << "\n";

    std::cout
        << "EXACT_Q_GCD="
        << exact_q_gcd
        << "\n";

    std::cout
        << "UNEXPECTED_P_GCD="
        << unexpected_p_gcd
        << "\n";

    std::cout
        << "UNEXPECTED_Q_GCD="
        << unexpected_q_gcd
        << "\n";

    std::cout
        << "P_VALUATION_AT_LEAST_2="
        << p_valuation_at_least_2
        << "\n";

    std::cout
        << "Q_VALUATION_AT_LEAST_2="
        << q_valuation_at_least_2
        << "\n";

    std::cout
        << "P_VALUATION_EXACTLY_2="
        << p_valuation_exactly_2
        << "\n";

    std::cout
        << "Q_VALUATION_EXACTLY_2="
        << q_valuation_exactly_2
        << "\n";

    std::cout
        << "P_VALUATION_AT_LEAST_3="
        << p_valuation_at_least_3
        << "\n";

    std::cout
        << "Q_VALUATION_AT_LEAST_3="
        << q_valuation_at_least_3
        << "\n";

    const bool structure_ok =
        p_squared_identity_failures == 0 &&
        q_squared_identity_failures == 0 &&
        p_squared_divisibility_failures == 0 &&
        q_squared_divisibility_failures == 0;

    std::cout
        << "STRUCTURE_STATUS="
        << (
            structure_ok
                ? "SECOND_ORDER_FACTOR_DIVISIBILITY_CONFIRMED"
                : "FAIL"
        )
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 362\n";

    return 0;
}