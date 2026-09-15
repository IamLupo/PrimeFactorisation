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
};

static mpz_class mpz_from_u64(u64 value) {
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
 * E(z) for
 *
 * R(z,u) = u + z + 1
 *
 * and
 *
 * E(z)
 *   = (s-z)^2 R(z,(z^2+D)/(s-z)).
 *
 * Expanded:
 *
 * E(z)
 *   = (s-z)(z^2+D)
 *     + z(s-z)^2
 *     + (s-z)^2.
 */
static mpz_class evaluate_E(
    u64 z,
    u64 s,
    u64 D
) {
    const mpz_class Z =
        mpz_from_u64(z);

    const mpz_class S =
        mpz_from_u64(s);

    const mpz_class Dz =
        mpz_from_u64(D);

    const mpz_class delta =
        S - Z;

    const mpz_class z2 =
        Z * Z;

    const mpz_class z2_plus_D =
        z2 + Dz;

    return
        delta * z2_plus_D
        +
        Z * delta * delta
        +
        delta * delta;
}

static mpz_class gcd_mpz(
    const mpz_class& value,
    const mpz_class& modulus
) {
    mpz_class x = value;

    if (x < 0) {
        x = -x;
    }

    mpz_class g;

    mpz_gcd(
        g.get_mpz_t(),
        x.get_mpz_t(),
        modulus.get_mpz_t()
    );

    return g;
}

int main() {
    std::cout << "START EXPERIMENT 363\n";

    constexpr int CASE_COUNT = 500;
    constexpr int PRIME_LIMIT = 50000;

    std::mt19937_64 rng(
        0x363363363ULL
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

    u64 total_scan_points = 0;

    u64 hidden_coordinate_failures = 0;
    u64 hidden_p2_failures = 0;

    u64 hidden_gcd_N_failures = 0;
    u64 hidden_gcd_N2_failures = 0;

    u64 cases_with_unique_p2_hit = 0;
    u64 cases_with_multiple_p2_hits = 0;
    u64 cases_with_no_p2_hit = 0;

    u64 total_p2_hits = 0;
    u64 total_p_hits = 0;
    u64 total_q_hits = 0;
    u64 total_N_hits = 0;

    u64 total_N2_p2_hits = 0;
    u64 total_N2_q2_hits = 0;
    u64 total_N2_other_hits = 0;

    u64 accidental_p2_hits = 0;
    u64 accidental_q2_hits = 0;

    u64 cases_hidden_is_first_p2 = 0;
    u64 cases_hidden_is_only_p2 = 0;

    u64 max_p2_hits_case = 0;

    u64 gcd_N2_one_count = 0;
    u64 gcd_N2_p_count = 0;
    u64 gcd_N2_q_count = 0;
    u64 gcd_N2_p2_count = 0;
    u64 gcd_N2_q2_count = 0;
    u64 gcd_N2_N_count = 0;
    u64 gcd_N2_N2_count = 0;
    u64 gcd_N2_other_count = 0;

    u64 total_hidden_a = 0;

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

        const u64 a =
            s - p;

        total_hidden_a += a;

        const mpz_class Nz =
            mpz_from_u64(N);

        const mpz_class N2 =
            Nz * Nz;

        const mpz_class pz =
            mpz_from_u64(p);

        const mpz_class qz =
            mpz_from_u64(q);

        const mpz_class p2 =
            pz * pz;

        const mpz_class q2 =
            qz * qz;

        /*
         * Hidden coordinate.
         */
        const mpz_class E_hidden =
            evaluate_E(
                a,
                s,
                D
            );

        const mpz_class gp_hidden =
            gcd_mpz(
                E_hidden,
                Nz
            );

        const mpz_class gp2_hidden =
            gcd_mpz(
                E_hidden,
                N2
            );

        if (gp_hidden != pz) {
            ++hidden_gcd_N_failures;
        }

        if (gp2_hidden != p2) {
            ++hidden_gcd_N2_failures;
        }

        if ((E_hidden % p2) != 0) {
            ++hidden_p2_failures;
        }

        if (a + p != s) {
            ++hidden_coordinate_failures;
        }

        u64 case_p2_hits = 0;
        u64 case_p_hits = 0;
        u64 case_q_hits = 0;

        u64 case_first_p2 =
            UINT64_MAX;

        for (u64 z = 0;
             z <= s;
             ++z) {

            ++total_scan_points;

            const mpz_class E =
                evaluate_E(
                    z,
                    s,
                    D
                );

            const mpz_class gN =
                gcd_mpz(
                    E,
                    Nz
                );

            const mpz_class gN2 =
                gcd_mpz(
                    E,
                    N2
                );

            /*
             * First-level factor hits.
             */
            if (gN == pz) {
                ++total_p_hits;
                ++case_p_hits;
            } else if (gN == qz) {
                ++total_q_hits;
                ++case_q_hits;
            } else if (gN == Nz) {
                ++total_N_hits;
            }

            /*
             * Second-level gcd classification.
             */
            if (gN2 == 1) {
                ++gcd_N2_one_count;
            } else if (gN2 == pz) {
                ++gcd_N2_p_count;
            } else if (gN2 == qz) {
                ++gcd_N2_q_count;
            } else if (gN2 == p2) {
                ++gcd_N2_p2_count;
                ++total_N2_p2_hits;

                ++case_p2_hits;

                if (z != a) {
                    ++accidental_p2_hits;
                }

                if (case_first_p2 == UINT64_MAX) {
                    case_first_p2 = z;
                }
            } else if (gN2 == q2) {
                ++gcd_N2_q2_count;
                ++total_N2_q2_hits;

                if (z != a) {
                    ++accidental_q2_hits;
                }
            } else if (gN2 == Nz) {
                ++gcd_N2_N_count;
            } else if (gN2 == N2) {
                ++gcd_N2_N2_count;
            } else {
                ++gcd_N2_other_count;
                ++total_N2_other_hits;
            }
        }

        total_p2_hits += case_p2_hits;

        if (case_p2_hits == 0) {
            ++cases_with_no_p2_hit;
        } else if (case_p2_hits == 1) {
            ++cases_with_unique_p2_hit;
        } else {
            ++cases_with_multiple_p2_hits;
        }

        max_p2_hits_case =
            std::max(
                max_p2_hits_case,
                case_p2_hits
            );

        if (case_first_p2 == a) {
            ++cases_hidden_is_first_p2;
        }

        if (
            case_p2_hits == 1 &&
            case_first_p2 == a
        ) {
            ++cases_hidden_is_only_p2;
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
                << "E(a)="
                << E_hidden
                << "\n";

            std::cout
                << "GCD_EA_N="
                << gp_hidden
                << "\n";

            std::cout
                << "GCD_EA_N2="
                << gp2_hidden
                << "\n";

            std::cout
                << "P2_DIVISIBLE="
                << (
                    (E_hidden % p2) == 0
                        ? "YES"
                        : "NO"
                )
                << "\n";

            std::cout
                << "P2_HITS_IN_SCAN="
                << case_p2_hits
                << "\n";

            std::cout
                << "FIRST_P2_HIT=";

            if (case_first_p2 == UINT64_MAX) {
                std::cout << -1;
            } else {
                std::cout << case_first_p2;
            }

            std::cout
                << "\n";

            std::cout
                << "HIDDEN_A_IS_FIRST_P2="
                << (
                    case_first_p2 == a
                        ? "YES"
                        : "NO"
                )
                << "\n";

            std::cout
                << "P_HITS="
                << case_p_hits
                << "\n";

            std::cout
                << "Q_HITS="
                << case_q_hits
                << "\n";
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "TOTAL_SCAN_POINTS="
        << total_scan_points
        << "\n";

    std::cout
        << "AVERAGE_HIDDEN_A="
        << total_hidden_a /
           static_cast<u64>(cases.size())
        << "\n";

    std::cout
        << "HIDDEN_COORDINATE_FAILURES="
        << hidden_coordinate_failures
        << "\n";

    std::cout
        << "HIDDEN_P2_FAILURES="
        << hidden_p2_failures
        << "\n";

    std::cout
        << "HIDDEN_GCD_N_FAILURES="
        << hidden_gcd_N_failures
        << "\n";

    std::cout
        << "HIDDEN_GCD_N2_FAILURES="
        << hidden_gcd_N2_failures
        << "\n";

    std::cout
        << "TOTAL_P2_HITS="
        << total_p2_hits
        << "\n";

    std::cout
        << "CASES_NO_P2_HIT="
        << cases_with_no_p2_hit
        << "\n";

    std::cout
        << "CASES_UNIQUE_P2_HIT="
        << cases_with_unique_p2_hit
        << "\n";

    std::cout
        << "CASES_MULTIPLE_P2_HITS="
        << cases_with_multiple_p2_hits
        << "\n";

    std::cout
        << "MAX_P2_HITS_CASE="
        << max_p2_hits_case
        << "\n";

    std::cout
        << "HIDDEN_A_IS_FIRST_P2="
        << cases_hidden_is_first_p2
        << "\n";

    std::cout
        << "HIDDEN_A_IS_ONLY_P2="
        << cases_hidden_is_only_p2
        << "\n";

    std::cout
        << "ACCIDENTAL_P2_HITS="
        << accidental_p2_hits
        << "\n";

    std::cout
        << "ACCIDENTAL_Q2_HITS="
        << accidental_q2_hits
        << "\n";

    std::cout
        << "TOTAL_P_HITS="
        << total_p_hits
        << "\n";

    std::cout
        << "TOTAL_Q_HITS="
        << total_q_hits
        << "\n";

    std::cout
        << "TOTAL_N_HITS="
        << total_N_hits
        << "\n";

    std::cout
        << "N2_GCD_1_COUNT="
        << gcd_N2_one_count
        << "\n";

    std::cout
        << "N2_GCD_P_COUNT="
        << gcd_N2_p_count
        << "\n";

    std::cout
        << "N2_GCD_Q_COUNT="
        << gcd_N2_q_count
        << "\n";

    std::cout
        << "N2_GCD_P2_COUNT="
        << gcd_N2_p2_count
        << "\n";

    std::cout
        << "N2_GCD_Q2_COUNT="
        << gcd_N2_q2_count
        << "\n";

    std::cout
        << "N2_GCD_N_COUNT="
        << gcd_N2_N_count
        << "\n";

    std::cout
        << "N2_GCD_N2_COUNT="
        << gcd_N2_N2_count
        << "\n";

    std::cout
        << "N2_GCD_OTHER_COUNT="
        << gcd_N2_other_count
        << "\n";

    if (
        hidden_coordinate_failures == 0 &&
        hidden_p2_failures == 0 &&
        hidden_gcd_N_failures == 0 &&
        hidden_gcd_N2_failures == 0
    ) {
        std::cout
            << "STRUCTURE_STATUS="
            << "PUBLIC_SCAN_P2_STRUCTURE_CONFIRMED"
            << "\n";
    } else {
        std::cout
            << "STRUCTURE_STATUS=FAIL"
            << "\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 363\n";

    return 0;
}