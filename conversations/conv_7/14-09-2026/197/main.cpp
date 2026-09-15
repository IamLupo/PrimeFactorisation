#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct Representation {
    u64 j;
    u64 k1;
    u64 k2;
    u64 value;
};

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
};

bool is_prime(u64 n) {
    if (n < 2) {
        return false;
    }

    if (n % 2 == 0) {
        return n == 2;
    }

    for (u64 d = 3; d * d <= n; d += 2) {
        if (n % d == 0) {
            return false;
        }
    }

    return true;
}

std::vector<Representation> find_representations(
    u64 p,
    u64 j_limit) {

    std::vector<Representation> result;

    for (u64 j = 1; j <= j_limit; ++j) {

        const u64 n =
            j * p + 1;

        /*
            We only need divisor pairs d1*d2=n with
            d1 == d2 == 1 (mod j).
        */
        for (u64 d1 = 3;
             d1 * d1 <= n;
             ++d1) {

            if (n % d1 != 0) {
                continue;
            }

            const u64 d2 =
                n / d1;

            if (d1 % j != 1 ||
                d2 % j != 1) {
                continue;
            }

            const u64 k1 =
                (d1 - 1) / j;

            const u64 k2 =
                (d2 - 1) / j;

            if (k1 == 0 ||
                k2 == 0) {
                continue;
            }

            const u64 reconstructed =
                k1 +
                k2 +
                j * k1 * k2;

            if (reconstructed != p) {
                continue;
            }

            result.push_back({
                j,
                k1,
                k2,
                p
            });
        }
    }

    std::sort(
        result.begin(),
        result.end(),
        [](const Representation& a,
           const Representation& b) {
            if (a.j != b.j) {
                return a.j < b.j;
            }

            if (a.k1 != b.k1) {
                return a.k1 < b.k1;
            }

            return a.k2 < b.k2;
        }
    );

    return result;
}

CaseData generate_case(
    const std::vector<int>& primes,
    std::mt19937_64& rng) {

    std::uniform_int_distribution<std::size_t> dist(
        0,
        primes.size() - 1
    );

    while (true) {
        u64 p =
            static_cast<u64>(
                primes[dist(rng)]
            );

        u64 q =
            static_cast<u64>(
                primes[dist(rng)]
            );

        if (p == q) {
            continue;
        }

        if (p > q) {
            std::swap(p, q);
        }

        return {
            p,
            q,
            p * q
        };
    }
}

double average(
    u64 total,
    u64 count) {

    if (count == 0) {
        return 0.0;
    }

    return
        static_cast<double>(total) /
        static_cast<double>(count);
}

int main() {
    constexpr int PRIME_LOW = 10000;
    constexpr int PRIME_HIGH = 100000;
    constexpr int CASE_COUNT = 500;

    constexpr u64 J_LIMIT = 50;

    std::cout
        << "START EXPERIMENT 490\n";

    std::vector<int> primes;

    for (int n = PRIME_LOW;
         n <= PRIME_HIGH;
         ++n) {

        if (is_prime(
                static_cast<u64>(n)
            )) {

            primes.push_back(n);
        }
    }

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << "\n";

    std::mt19937_64 rng(
        4901234567ULL
    );

    u64 representation_cases = 0;

    u64 pair_count = 0;

    u64 invariant_failures = 0;

    u64 total_jp = 0;
    u64 total_jq = 0;

    u64 total_j_sum = 0;
    u64 total_j_product = 0;

    u64 first_N = 0;
    u64 first_p = 0;
    u64 first_q = 0;

    const auto start =
        std::chrono::steady_clock::now();

    for (int case_id = 0;
         case_id < CASE_COUNT;
         ++case_id) {

        const CaseData data =
            generate_case(
                primes,
                rng
            );

        const auto reps_p =
            find_representations(
                data.p,
                J_LIMIT
            );

        const auto reps_q =
            find_representations(
                data.q,
                J_LIMIT
            );

        if (case_id == 0) {
            first_N = data.N;
            first_p = data.p;
            first_q = data.q;
        }

        if (reps_p.empty() ||
            reps_q.empty()) {

            continue;
        }

        ++representation_cases;

        for (const auto& rp : reps_p) {

            for (const auto& rq : reps_q) {

                ++pair_count;

                const u64 jp =
                    rp.j;

                const u64 jq =
                    rq.j;

                const u64 X =
                    (jp * data.p + 1) *
                    (jq * data.q + 1);

                const u64 expanded =
                    jp * jq * data.N +
                    jp * data.p +
                    jq * data.q +
                    1;

                if (X != expanded) {
                    ++invariant_failures;
                }

                total_jp += jp;
                total_jq += jq;

                total_j_sum +=
                    jp + jq;

                total_j_product +=
                    jp * jq;
            }
        }

        if ((case_id + 1) % 100 == 0) {
            std::cout
                << "PROGRESS="
                << (case_id + 1)
                << "/"
                << CASE_COUNT
                << "\n";
        }
    }

    const auto end =
        std::chrono::steady_clock::now();

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << "\n";

    std::cout
        << "CASES_WITH_BOTH_REPRESENTATIONS="
        << representation_cases
        << "\n";

    std::cout
        << "TOTAL_J_PAIR_COMBINATIONS="
        << pair_count
        << "\n";

    std::cout
        << "IDENTITY_FAILURES="
        << invariant_failures
        << "\n";

    std::cout
        << "AVG_JP="
        << average(
               total_jp,
               pair_count
           )
        << "\n";

    std::cout
        << "AVG_JQ="
        << average(
               total_jq,
               pair_count
           )
        << "\n";

    std::cout
        << "AVG_J_SUM="
        << average(
               total_j_sum,
               pair_count
           )
        << "\n";

    std::cout
        << "AVG_J_PRODUCT="
        << average(
               total_j_product,
               pair_count
           )
        << "\n";

    std::cout
        << "FIRST_N="
        << first_N
        << "\n";

    std::cout
        << "FIRST_P="
        << first_p
        << "\n";

    std::cout
        << "FIRST_Q="
        << first_q
        << "\n";

    /*
        Detailed first case.
    */
    const auto first_p_reps =
        find_representations(
            first_p,
            J_LIMIT
        );

    const auto first_q_reps =
        find_representations(
            first_q,
            J_LIMIT
        );

    std::cout
        << "FIRST_P_REPRESENTATIONS="
        << first_p_reps.size()
        << "\n";

    std::cout
        << "FIRST_Q_REPRESENTATIONS="
        << first_q_reps.size()
        << "\n";

    for (const auto& rp :
         first_p_reps) {

        for (const auto& rq :
             first_q_reps) {

            const u64 jp =
                rp.j;

            const u64 jq =
                rq.j;

            const u64 X =
                (jp * first_p + 1) *
                (jq * first_q + 1);

            const u64 expanded =
                jp * jq * first_N +
                jp * first_p +
                jq * first_q +
                1;

            std::cout
                << "FIRST_PAIR"
                << " JP="
                << jp
                << " JQ="
                << jq
                << " JP_PLUS_JQ="
                << jp + jq
                << " JP_TIMES_JQ="
                << jp * jq
                << " X="
                << X
                << " EXPANDED="
                << expanded
                << "\n";
        }
    }

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 490\n";

    return 0;
}
