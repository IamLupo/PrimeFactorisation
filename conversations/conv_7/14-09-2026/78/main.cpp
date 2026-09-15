#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
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
    u64 x = static_cast<u64>(
        std::sqrt(static_cast<long double>(n))
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
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

static u64 random_prime(
    std::mt19937_64& rng,
    u64 lo,
    u64 hi
) {
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

        if (x <= hi && x >= lo) {
            return x;
        }
    }
}

static CaseData make_case(
    std::mt19937_64& rng
) {
    while (true) {
        u64 p = random_prime(rng, 5, 3000);
        u64 q = random_prime(rng, 3001, 10000);

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        const u64 N = p * q;
        const u64 s = isqrt_u64(N);

        if (s * s > N) {
            continue;
        }

        if ((s + 1) * (s + 1) <= N) {
            continue;
        }

        return {
            p,
            q,
            N,
            s,
            N - s * s
        };
    }
}

/*
    E(z) = (s-z)(D+s+z(s-1))
*/
static mpz_class E(
    u64 s,
    u64 D,
    u64 z
) {
    const mpz_class S = mpz_from_u64(s);
    const mpz_class Z = mpz_from_u64(z);
    const mpz_class DD = mpz_from_u64(D);

    return (S - Z) *
           (DD + S + Z * (S - 1));
}

int main() {
    constexpr int EXPERIMENT = 367;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x36720260914ULL);

    u64 reconstruction_s_failures = 0;
    u64 reconstruction_n_failures = 0;

    u64 s_integer_failures = 0;
    u64 s_floor_failures = 0;
    u64 n_positive_failures = 0;

    u64 sequence_failures = 0;
    u64 sampled_identity_failures = 0;

    u64 max_z_tested = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 s = c.s;
        const u64 D = c.D;

        /*
            z = s-1
            => E(s-1) = N-s+1
        */
        const mpz_class E1 =
            E(s, D, s - 1);

        /*
            z = s-2
            => E(s-2) = 2(N-2(s-1))
        */
        const mpz_class E2 =
            E(s, D, s - 2);

        /*
            Recover s from the two public evaluations.

                2E1 - E2 = 2(s-1)

            Therefore

                s = 1 + (2E1-E2)/2
        */
        const mpz_class numerator =
            2 * E1 - E2;

        if ((numerator % 2) != 0) {
            ++s_integer_failures;
            continue;
        }

        const mpz_class recovered_s_mpz =
            numerator / 2 + 1;

        if (recovered_s_mpz <= 0) {
            ++s_integer_failures;
            continue;
        }

        const u64 recovered_s =
            recovered_s_mpz.get_ui();

        if (recovered_s != s) {
            ++reconstruction_s_failures;
        }

        /*
            Recover N from

                N = E(s-1) + s - 1.
        */
        const mpz_class recovered_n_mpz =
            E1 + recovered_s_mpz - 1;

        const u64 recovered_n =
            recovered_n_mpz.get_ui();

        if (recovered_n != c.N) {
            ++reconstruction_n_failures;
        }

        if (recovered_n == 0) {
            ++n_positive_failures;
        }

        /*
            Check that recovered s is actually floor(sqrt(N)).
        */
        if (recovered_s * recovered_s > recovered_n) {
            ++s_floor_failures;
        } else if ((recovered_s + 1) *
                   (recovered_s + 1) <= recovered_n) {
            ++s_floor_failures;
        }

        /*
            Now reconstruct D from recovered N and s.
        */
        const mpz_class reconstructed_D =
            recovered_n_mpz -
            recovered_s_mpz * recovered_s_mpz;

        /*
            The recovered (N,s,D) should reproduce the
            complete E sequence.
        */
        const u64 upper =
            std::min<u64>(s, 100);

        for (u64 z = 0; z <= upper; ++z) {
            max_z_tested =
                std::max(max_z_tested, z);

            const mpz_class original =
                E(s, D, z);

            const mpz_class reconstructed =
                (recovered_s_mpz - mpz_from_u64(z)) *
                (
                    reconstructed_D +
                    recovered_s_mpz +
                    mpz_from_u64(z) *
                    (recovered_s_mpz - 1)
                );

            if (original != reconstructed) {
                ++sequence_failures;
            }
        }

        /*
            Directly verify the two reconstruction identities.
        */
        const mpz_class identity_s =
            1 + (2 * E1 - E2) / 2;

        if (identity_s != mpz_from_u64(s)) {
            ++sampled_identity_failures;
        }

        const mpz_class identity_n =
            E1 + identity_s - 1;

        if (identity_n != mpz_from_u64(c.N)) {
            ++sampled_identity_failures;
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "RECONSTRUCTION_S_FAILURES="
        << reconstruction_s_failures
        << '\n';

    std::cout
        << "RECONSTRUCTION_N_FAILURES="
        << reconstruction_n_failures
        << '\n';

    std::cout
        << "S_INTEGER_FAILURES="
        << s_integer_failures
        << '\n';

    std::cout
        << "S_FLOOR_FAILURES="
        << s_floor_failures
        << '\n';

    std::cout
        << "N_POSITIVE_FAILURES="
        << n_positive_failures
        << '\n';

    std::cout
        << "SEQUENCE_RECONSTRUCTION_FAILURES="
        << sequence_failures
        << '\n';

    std::cout
        << "DIRECT_IDENTITY_FAILURES="
        << sampled_identity_failures
        << '\n';

    std::cout
        << "MAX_Z_TESTED="
        << max_z_tested
        << '\n';

    const bool status =
        reconstruction_s_failures == 0 &&
        reconstruction_n_failures == 0 &&
        s_integer_failures == 0 &&
        s_floor_failures == 0 &&
        n_positive_failures == 0 &&
        sequence_failures == 0 &&
        sampled_identity_failures == 0;

    std::cout
        << "RECONSTRUCTION_STATUS="
        << (status ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
