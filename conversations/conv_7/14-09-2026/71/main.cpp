#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
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

int main() {
    std::cout << "START EXPERIMENT 360\n";

    constexpr int CASE_COUNT = 2000;
    constexpr int PRIME_LIMIT = 100000;

    std::mt19937_64 rng(
        0x360360360ULL
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

    u64 coordinate_failures = 0;
    u64 quotient_difference_failures = 0;
    u64 quotient_center_failures = 0;

    u64 fermat_identity_failures = 0;
    u64 N_reconstruction_failures = 0;

    u64 polynomial_center_failures = 0;
    u64 polynomial_diff_failures = 0;

    u64 parity_failures = 0;

    u64 total_x = 0;
    u64 total_A_minus_s = 0;
    u64 total_B = 0;

    u64 min_A_minus_s = UINT64_MAX;
    u64 max_A_minus_s = 0;

    u64 min_B = UINT64_MAX;
    u64 max_B = 0;

    u64 cases_A_equals_s_plus_1 = 0;
    u64 cases_B_equals_1 = 0;

    u64 fermat_iterations_sum = 0;
    u64 polynomial_x_search_sum = 0;

    u64 fermat_iterations_lt_10 = 0;
    u64 fermat_iterations_lt_100 = 0;
    u64 fermat_iterations_lt_1000 = 0;

    u64 polynomial_x_lt_10 = 0;
    u64 polynomial_x_lt_100 = 0;
    u64 polynomial_x_lt_1000 = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c =
            cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 D =
            N - s * s;

        /*
         * ---------------------------------------------------------
         * Polynomial coordinate.
         *
         * x_p = s+1-p.
         * ---------------------------------------------------------
         */

        const u64 x =
            s + 1 - p;

        const u64 a =
            s - p;

        const u64 b =
            q - s;

        if (x != a + 1) {
            ++coordinate_failures;
        }

        /*
         * ---------------------------------------------------------
         * Polynomial quotient coordinates.
         * ---------------------------------------------------------
         */

        const u64 difference =
            q - p;

        const i128 center_signed =
            static_cast<i128>(p)
            +
            static_cast<i128>(q)
            -
            static_cast<i128>(2) *
            static_cast<i128>(s);

        if (center_signed < 0) {
            ++coordinate_failures;
            continue;
        }

        const u64 center =
            static_cast<u64>(center_signed);

        /*
         * ---------------------------------------------------------
         * Fermat coordinates.
         *
         * A=(p+q)/2
         * B=(q-p)/2
         * ---------------------------------------------------------
         */

        const u64 sum =
            p + q;

        if ((sum & 1ULL) != 0 ||
            (difference & 1ULL) != 0) {

            ++parity_failures;
            continue;
        }

        const u64 A =
            sum / 2;

        const u64 B =
            difference / 2;

        /*
         * ---------------------------------------------------------
         * Central identities.
         * ---------------------------------------------------------
         */

        if (
            center !=
            2 * (A - s)
        ) {
            ++quotient_center_failures;
        }

        if (
            difference !=
            2 * B
        ) {
            ++quotient_difference_failures;
        }

        /*
         * a,b version:
         *
         * center=b-a
         * difference=a+b
         */
        if (
            center != b - a ||
            difference != a + b
        ) {
            ++coordinate_failures;
        }

        /*
         * ---------------------------------------------------------
         * Fermat identity:
         *
         * A^2-B^2=N.
         * ---------------------------------------------------------
         */

        const u128 A2 =
            static_cast<u128>(A) *
            static_cast<u128>(A);

        const u128 B2 =
            static_cast<u128>(B) *
            static_cast<u128>(B);

        if (A2 - B2 !=
            static_cast<u128>(N)) {

            ++fermat_identity_failures;
        }

        /*
         * ---------------------------------------------------------
         * Reconstruct N from s and the two quotient coordinates:
         *
         * N =
         * s^2 + s*center
         * + (center^2-difference^2)/4.
         * ---------------------------------------------------------
         */

        const i128 reconstruction =
            static_cast<i128>(s) *
                static_cast<i128>(s)
            +
            static_cast<i128>(s) *
                static_cast<i128>(center)
            +
            (
                static_cast<i128>(center) *
                    static_cast<i128>(center)
                -
                static_cast<i128>(difference) *
                    static_cast<i128>(difference)
            ) / 4;

        if (
            reconstruction !=
            static_cast<i128>(N)
        ) {
            ++N_reconstruction_failures;
        }

        /*
         * ---------------------------------------------------------
         * Polynomial identities.
         *
         * P_center(x)=D+(x-1)^2
         *            = p*(p+q-2s)
         *            = 2p(A-s)
         *
         * P_diff(x)=p*(q-p)=2pB.
         * ---------------------------------------------------------
         */

        const i128 center_poly =
            static_cast<i128>(D)
            +
            (
                static_cast<i128>(x) -
                static_cast<i128>(1)
            ) *
            (
                static_cast<i128>(x) -
                static_cast<i128>(1)
            );

        const i128 expected_center_poly =
            static_cast<i128>(p) *
            static_cast<i128>(center);

        if (
            center_poly !=
            expected_center_poly
        ) {
            ++polynomial_center_failures;
        }

        const i128 diff_poly =
            static_cast<i128>(p) *
            static_cast<i128>(difference);

        if (
            abs_i128_mod(
                diff_poly,
                p
            ) != 0
        ) {
            ++polynomial_diff_failures;
        }

        /*
         * ---------------------------------------------------------
         * Search-coordinate comparison.
         *
         * A classical Fermat search starts at s+1 and searches A.
         *
         * A polynomial search in x starts from x=1 upward if
         * looking for the hidden p-coordinate.
         *
         * These are related by
         *
         * A-s = x-1.
         * ---------------------------------------------------------
         */

        const u64 fermat_offset =
            A - s;

        const u64 polynomial_offset =
            x - 1;

        if (
            fermat_offset !=
            polynomial_offset
        ) {
            ++coordinate_failures;
        }

        /*
         * The Fermat square is:
         *
         * B^2 = A^2-N.
         */
        const u64 fermat_iterations =
            fermat_offset;

        const u64 polynomial_iterations =
            polynomial_offset;

        fermat_iterations_sum +=
            fermat_iterations;

        polynomial_x_search_sum +=
            polynomial_iterations;

        if (fermat_iterations < 10) {
            ++fermat_iterations_lt_10;
        }

        if (fermat_iterations < 100) {
            ++fermat_iterations_lt_100;
        }

        if (fermat_iterations < 1000) {
            ++fermat_iterations_lt_1000;
        }

        if (polynomial_iterations < 10) {
            ++polynomial_x_lt_10;
        }

        if (polynomial_iterations < 100) {
            ++polynomial_x_lt_100;
        }

        if (polynomial_iterations < 1000) {
            ++polynomial_x_lt_1000;
        }

        min_A_minus_s =
            std::min(
                min_A_minus_s,
                A - s
            );

        max_A_minus_s =
            std::max(
                max_A_minus_s,
                A - s
            );

        min_B =
            std::min(
                min_B,
                B
            );

        max_B =
            std::max(
                max_B,
                B
            );

        if (A == s + 1) {
            ++cases_A_equals_s_plus_1;
        }

        if (B == 1) {
            ++cases_B_equals_1;
        }

        total_x += x;
        total_A_minus_s += A - s;
        total_B += B;

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
                << "x=s+1-p="
                << x
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
                << "q-p="
                << difference
                << "\n";

            std::cout
                << "p+q-2s="
                << center
                << "\n";

            std::cout
                << "A="
                << A
                << "\n";

            std::cout
                << "B="
                << B
                << "\n";

            std::cout
                << "A-s="
                << (A - s)
                << "\n";

            std::cout
                << "B="
                << B
                << "\n";

            std::cout
                << "A^2-B^2="
                << static_cast<u64>(
                    A2 - B2
                )
                << "\n";

            std::cout
                << "FERMAT_OFFSET="
                << fermat_offset
                << "\n";

            std::cout
                << "POLYNOMIAL_OFFSET="
                << polynomial_offset
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "COORDINATE_FAILURES="
        << coordinate_failures
        << "\n";

    std::cout
        << "QUOTIENT_DIFFERENCE_FAILURES="
        << quotient_difference_failures
        << "\n";

    std::cout
        << "QUOTIENT_CENTER_FAILURES="
        << quotient_center_failures
        << "\n";

    std::cout
        << "FERMAT_IDENTITY_FAILURES="
        << fermat_identity_failures
        << "\n";

    std::cout
        << "N_RECONSTRUCTION_FAILURES="
        << N_reconstruction_failures
        << "\n";

    std::cout
        << "POLYNOMIAL_CENTER_FAILURES="
        << polynomial_center_failures
        << "\n";

    std::cout
        << "POLYNOMIAL_DIFF_FAILURES="
        << polynomial_diff_failures
        << "\n";

    std::cout
        << "PARITY_FAILURES="
        << parity_failures
        << "\n";

    std::cout
        << "MIN_A_MINUS_S="
        << min_A_minus_s
        << "\n";

    std::cout
        << "MAX_A_MINUS_S="
        << max_A_minus_s
        << "\n";

    std::cout
        << "MIN_B="
        << min_B
        << "\n";

    std::cout
        << "MAX_B="
        << max_B
        << "\n";

    std::cout
        << "AVERAGE_FERMAT_OFFSET="
        << fermat_iterations_sum /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "AVERAGE_POLYNOMIAL_OFFSET="
        << polynomial_x_search_sum /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "TOTAL_B="
        << total_B
        << "\n";

    std::cout
        << "FERMAT_OFFSET_LT_10="
        << fermat_iterations_lt_10
        << "\n";

    std::cout
        << "FERMAT_OFFSET_LT_100="
        << fermat_iterations_lt_100
        << "\n";

    std::cout
        << "FERMAT_OFFSET_LT_1000="
        << fermat_iterations_lt_1000
        << "\n";

    std::cout
        << "POLYNOMIAL_OFFSET_LT_10="
        << polynomial_x_lt_10
        << "\n";

    std::cout
        << "POLYNOMIAL_OFFSET_LT_100="
        << polynomial_x_lt_100
        << "\n";

    std::cout
        << "POLYNOMIAL_OFFSET_LT_1000="
        << polynomial_x_lt_1000
        << "\n";

    std::cout
        << "A_EQUALS_S_PLUS_1="
        << cases_A_equals_s_plus_1
        << "\n";

    std::cout
        << "B_EQUALS_1="
        << cases_B_equals_1
        << "\n";

    if (
        coordinate_failures == 0 &&
        quotient_difference_failures == 0 &&
        quotient_center_failures == 0 &&
        fermat_identity_failures == 0 &&
        N_reconstruction_failures == 0 &&
        polynomial_center_failures == 0 &&
        polynomial_diff_failures == 0 &&
        parity_failures == 0
    ) {
        std::cout
            << "STRUCTURE_STATUS="
            << "POLYNOMIAL_FERMAT_EQUIVALENCE_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 360\n";

    return 0;
}
