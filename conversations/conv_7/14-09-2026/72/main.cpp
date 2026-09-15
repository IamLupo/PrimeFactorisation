#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

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
        if (p >= 1009 && p <= 100000) {
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

static void print_i128(i128 x) {
    if (x == 0) {
        std::cout << "0";
        return;
    }

    if (x < 0) {
        std::cout << '-';
        x = -x;
    }

    std::string out;

    while (x > 0) {
        const unsigned digit =
            static_cast<unsigned>(x % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(
        out.begin(),
        out.end()
    );

    std::cout << out;
}

static i128 sq(i128 x) {
    return x * x;
}

static u64 abs_i128_mod(
    i128 x,
    u64 mod
) {
    if (x < 0) {
        x = -x;
    }

    return static_cast<u64>(
        static_cast<u128>(x) %
        static_cast<u128>(mod)
    );
}

int main() {
    std::cout << "START EXPERIMENT 361\n";

    constexpr int CASE_COUNT = 2000;
    constexpr int PRIME_LIMIT = 100000;

    std::mt19937_64 rng(
        0x361361361ULL
    );

    const std::vector<int> primes =
        sieve_primes(PRIME_LIMIT);

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    u64 a_definition_failures = 0;
    u64 b_definition_failures = 0;

    u64 D_identity_failures = 0;

    u64 quadratic_at_a_failures = 0;
    u64 quadratic_at_minus_b_failures = 0;

    u64 root_sum_failures = 0;
    u64 root_product_failures = 0;

    u64 discriminant_failures = 0;

    u64 discriminant_square_failures = 0;
    u64 discriminant_sqrt_failures = 0;

    u64 reconstruction_p_failures = 0;
    u64 reconstruction_q_failures = 0;

    u64 discriminant_equals_q_minus_p = 0;

    u64 total_a = 0;
    u64 total_b = 0;
    u64 total_u = 0;

    u64 min_a = UINT64_MAX;
    u64 max_a = 0;

    u64 min_b = UINT64_MAX;
    u64 max_b = 0;

    u64 min_u = UINT64_MAX;
    u64 max_u = 0;

    /*
     * Compare the quadratic directly with the two previous
     * coordinate systems:
     *
     * a = s-p
     *
     * y = (p+q-2s)/2
     *
     * and the usual Fermat discriminant.
     */
    u64 relation_u_center_failures = 0;
    u64 relation_discriminant_failures = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 D =
            N - s * s;

        /*
         * Hidden deviations:
         *
         * a = s-p
         * b = q-s
         */
        const u64 a =
            s - p;

        const u64 b =
            q - s;

        /*
         * u = b-a = p+q-2s.
         */
        const i128 u_signed =
            static_cast<i128>(b)
            -
            static_cast<i128>(a);

        const u64 u =
            u_signed < 0
                ? 0
                : static_cast<u64>(u_signed);

        /*
         * These are known positive for the generated cases,
         * but verify the definitions explicitly.
         */
        if (
            static_cast<i128>(a) !=
            static_cast<i128>(s) -
            static_cast<i128>(p)
        ) {
            ++a_definition_failures;
        }

        if (
            static_cast<i128>(b) !=
            static_cast<i128>(q) -
            static_cast<i128>(s)
        ) {
            ++b_definition_failures;
        }

        /*
         * ---------------------------------------------------------
         * Fundamental identity:
         *
         * D = su - ab.
         * ---------------------------------------------------------
         */
        const i128 D_from_ab =
            static_cast<i128>(s) *
            static_cast<i128>(u)
            -
            static_cast<i128>(a) *
            static_cast<i128>(b);

        if (
            D_from_ab !=
            static_cast<i128>(D)
        ) {
            ++D_identity_failures;
        }

        /*
         * ---------------------------------------------------------
         * Quadratic:
         *
         * F(z)=z^2 + uz + D-su.
         *
         * Expected factorization:
         *
         * F(z)=(z-a)(z+b).
         * ---------------------------------------------------------
         */
        const i128 U =
            static_cast<i128>(u);

        const i128 A128 =
            static_cast<i128>(a);

        const i128 B128 =
            static_cast<i128>(b);

        const i128 D128 =
            static_cast<i128>(D);

        const i128 constant =
            D128 -
            static_cast<i128>(s) * U;

        const i128 F_at_a =
            sq(A128)
            +
            U * A128
            +
            constant;

        const i128 F_at_minus_b =
            sq(-B128)
            +
            U * (-B128)
            +
            constant;

        if (F_at_a != 0) {
            ++quadratic_at_a_failures;
        }

        if (F_at_minus_b != 0) {
            ++quadratic_at_minus_b_failures;
        }

        /*
         * Root sum/product:
         *
         * roots = a, -b
         *
         * sum = a-b = -u
         * product = -ab = D-su
         */
        const i128 root_sum =
            A128 - B128;

        const i128 root_product =
            -A128 * B128;

        if (
            root_sum != -U
        ) {
            ++root_sum_failures;
        }

        if (
            root_product != constant
        ) {
            ++root_product_failures;
        }

        /*
         * ---------------------------------------------------------
         * Discriminant:
         *
         * Delta = u^2 - 4(D-su)
         *
         *          = (a+b)^2
         *
         *          = (q-p)^2.
         * ---------------------------------------------------------
         */
        const i128 Delta =
            U * U
            -
            static_cast<i128>(4) *
            constant;

        const u64 q_minus_p =
            q - p;

        const i128 expected_delta =
            static_cast<i128>(q_minus_p) *
            static_cast<i128>(q_minus_p);

        if (
            Delta != expected_delta
        ) {
            ++discriminant_failures;
        }

        if (
            Delta !=
            sq(
                static_cast<i128>(a + b)
            )
        ) {
            ++discriminant_square_failures;
        }

        /*
         * Exact square root is q-p.
         */
        if (Delta >= 0) {
            const u64 delta_u =
                static_cast<u64>(Delta);

            const u64 z =
                isqrt_u64(delta_u);

            if (z != q_minus_p) {
                ++discriminant_sqrt_failures;
            }
        }

        /*
         * ---------------------------------------------------------
         * Recover the roots algebraically:
         *
         * z1=(-u + sqrt(Delta))/2 = a
         * z2=(-u - sqrt(Delta))/2 = -b
         * ---------------------------------------------------------
         */
        const i128 delta_sqrt =
            static_cast<i128>(q_minus_p);

        const i128 numerator_a =
            -U + delta_sqrt;

        const i128 numerator_b =
            -U - delta_sqrt;

        if (
            numerator_a !=
            static_cast<i128>(2) *
            static_cast<i128>(a)
        ) {
            ++reconstruction_p_failures;
        }

        if (
            numerator_b !=
            -static_cast<i128>(2) *
            static_cast<i128>(b)
        ) {
            ++reconstruction_q_failures;
        }

        /*
         * ---------------------------------------------------------
         * The quadratic now also gives:
         *
         * p = s-a
         * q = s+b
         *
         * and therefore factors N once the correct root is known.
         * ---------------------------------------------------------
         */

        const u64 reconstructed_p =
            s - a;

        const u64 reconstructed_q =
            s + b;

        if (reconstructed_p != p) {
            ++reconstruction_p_failures;
        }

        if (reconstructed_q != q) {
            ++reconstruction_q_failures;
        }

        /*
         * ---------------------------------------------------------
         * Relation to Fermat:
         *
         * y=(p+q-2s)/2 = u/2.
         *
         * Delta=(q-p)^2.
         *
         * Thus:
         *
         * y^2 + sy + D = a(a+u)
         *
         * and the quadratic is exactly the translated
         * difference-of-squares equation.
         * ---------------------------------------------------------
         */

        if (
            (u & 1ULL) != 0
        ) {
            ++relation_u_center_failures;
        }

        const u64 y =
            u / 2;

        const u128 fermat_left =
            (
                static_cast<u128>(s)
                +
                static_cast<u128>(y)
            )
            *
            (
                static_cast<u128>(s)
                +
                static_cast<u128>(y)
            );

        const u64 B =
            q_minus_p / 2;

        const u128 fermat_right =
            static_cast<u128>(B) *
            static_cast<u128>(B)
            +
            static_cast<u128>(N);

        if (
            fermat_left != fermat_right
        ) {
            ++relation_discriminant_failures;
        }

        /*
         * Statistics.
         */
        total_a += a;
        total_b += b;
        total_u += u;

        min_a = std::min(min_a, a);
        max_a = std::max(max_a, a);

        min_b = std::min(min_b, b);
        max_b = std::max(max_b, b);

        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);

        if (Delta ==
            expected_delta) {
            ++discriminant_equals_q_minus_p;
        }

        if (case_id < 10) {
            std::cout
                << "\nCASE="
                << (case_id + 1)
                << " p=" << p
                << " q=" << q
                << " N=" << N
                << " s=" << s
                << "\n";

            std::cout
                << "a=s-p="
                << a
                << "\n";

            std::cout
                << "b=q-s="
                << b
                << "\n";

            std::cout
                << "u=b-a="
                << u
                << "\n";

            std::cout
                << "D="
                << D
                << "\n";

            std::cout
                << "D_from_su_ab=";

            print_i128(D_from_ab);

            std::cout
                << "\n";

            std::cout
                << "F(a)=";

            print_i128(F_at_a);

            std::cout
                << "\n";

            std::cout
                << "F(-b)=";

            print_i128(F_at_minus_b);

            std::cout
                << "\n";

            std::cout
                << "DISCRIMINANT=";

            print_i128(Delta);

            std::cout
                << "\n";

            std::cout
                << "(q-p)^2="
                << q_minus_p * q_minus_p
                << "\n";

            std::cout
                << "SQRT_DISCRIMINANT="
                << q_minus_p
                << "\n";

            std::cout
                << "ROOT_1=";

            print_i128(
                numerator_a / 2
            );

            std::cout
                << "\n";

            std::cout
                << "ROOT_2=";

            print_i128(
                numerator_b / 2
            );

            std::cout
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "A_DEFINITION_FAILURES="
        << a_definition_failures
        << "\n";

    std::cout
        << "B_DEFINITION_FAILURES="
        << b_definition_failures
        << "\n";

    std::cout
        << "D_IDENTITY_FAILURES="
        << D_identity_failures
        << "\n";

    std::cout
        << "QUADRATIC_AT_A_FAILURES="
        << quadratic_at_a_failures
        << "\n";

    std::cout
        << "QUADRATIC_AT_MINUS_B_FAILURES="
        << quadratic_at_minus_b_failures
        << "\n";

    std::cout
        << "ROOT_SUM_FAILURES="
        << root_sum_failures
        << "\n";

    std::cout
        << "ROOT_PRODUCT_FAILURES="
        << root_product_failures
        << "\n";

    std::cout
        << "DISCRIMINANT_FAILURES="
        << discriminant_failures
        << "\n";

    std::cout
        << "DISCRIMINANT_SQUARE_FAILURES="
        << discriminant_square_failures
        << "\n";

    std::cout
        << "DISCRIMINANT_SQRT_FAILURES="
        << discriminant_sqrt_failures
        << "\n";

    std::cout
        << "RECONSTRUCTION_P_FAILURES="
        << reconstruction_p_failures
        << "\n";

    std::cout
        << "RECONSTRUCTION_Q_FAILURES="
        << reconstruction_q_failures
        << "\n";

    std::cout
        << "U_CENTER_RELATION_FAILURES="
        << relation_u_center_failures
        << "\n";

    std::cout
        << "FERMAT_RELATION_FAILURES="
        << relation_discriminant_failures
        << "\n";

    std::cout
        << "DISCRIMINANT_EQ_Q_MINUS_P_SQUARED="
        << discriminant_equals_q_minus_p
        << "\n";

    std::cout
        << "MIN_A="
        << min_a
        << "\n";

    std::cout
        << "MAX_A="
        << max_a
        << "\n";

    std::cout
        << "MIN_B="
        << min_b
        << "\n";

    std::cout
        << "MAX_B="
        << max_b
        << "\n";

    std::cout
        << "MIN_U="
        << min_u
        << "\n";

    std::cout
        << "MAX_U="
        << max_u
        << "\n";

    std::cout
        << "AVG_A="
        << total_a /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "AVG_B="
        << total_b /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "AVG_U="
        << total_u /
           static_cast<u64>(cases.size())
        << "\n";

    const bool structure_ok =
        a_definition_failures == 0 &&
        b_definition_failures == 0 &&
        D_identity_failures == 0 &&
        quadratic_at_a_failures == 0 &&
        quadratic_at_minus_b_failures == 0 &&
        root_sum_failures == 0 &&
        root_product_failures == 0 &&
        discriminant_failures == 0 &&
        discriminant_square_failures == 0 &&
        discriminant_sqrt_failures == 0 &&
        reconstruction_p_failures == 0 &&
        reconstruction_q_failures == 0 &&
        relation_u_center_failures == 0 &&
        relation_discriminant_failures == 0;

    if (structure_ok) {
        std::cout
            << "STRUCTURE_STATUS="
            << "DEVIATION_QUADRATIC_AND_FERMAT_DISCRIMINANT_EXACT"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 361\n";

    return 0;
}
