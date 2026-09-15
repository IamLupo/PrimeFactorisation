#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Candidate {
    bool valid = false;
    u64 r = 0;
    u64 k = 0;
    int m = 0;
    int sign = 0;
    u64 d = 0;
};

struct CaseResult {
    bool success = false;
    bool core_success = false;
    bool extended_success = false;

    u64 N = 0;
    u64 p = 0;
    u64 q = 0;
    u64 s = 0;

    Candidate best_core;
    Candidate best_extended;
};

static std::vector<u64> generate_primes(
    int limit
) {
    std::vector<bool> sieve(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    if (limit >= 0) {
        sieve[0] = false;
    }

    if (limit >= 1) {
        sieve[1] = false;
    }

    for (
        int p = 2;
        static_cast<long long>(p) * p <= limit;
        ++p
    ) {
        if (!sieve[p]) {
            continue;
        }

        for (
            int x = p * p;
            x <= limit;
            x += p
        ) {
            sieve[x] = false;
        }
    }

    std::vector<u64> primes;

    for (int x = 2; x <= limit; ++x) {
        if (sieve[x]) {
            primes.push_back(
                static_cast<u64>(x)
            );
        }
    }

    return primes;
}

static u64 integer_sqrt(
    u64 n
) {
    u64 x =
        static_cast<u64>(
            __builtin_sqrt(
                static_cast<long double>(n)
            )
        );

    while (
        static_cast<i128>(x + 1) *
        static_cast<i128>(x + 1) <=
        static_cast<i128>(n)
    ) {
        ++x;
    }

    while (
        static_cast<i128>(x) *
        static_cast<i128>(x) >
        static_cast<i128>(n)
    ) {
        --x;
    }

    return x;
}

static u64 inverse_mod(
    u64 a,
    u64 mod
) {
    i128 t = 0;
    i128 new_t = 1;

    i128 r = static_cast<i128>(mod);
    i128 new_r = static_cast<i128>(a);

    while (new_r != 0) {
        const i128 q = r / new_r;

        const i128 temp_t = t - q * new_t;
        t = new_t;
        new_t = temp_t;

        const i128 temp_r = r - q * new_r;
        r = new_r;
        new_r = temp_r;
    }

    if (r != 1) {
        return 0;
    }

    t %= static_cast<i128>(mod);

    if (t < 0) {
        t += static_cast<i128>(mod);
    }

    return static_cast<u64>(t);
}

static Candidate blind_candidate(
    u64 N,
    u64 s,
    int m,
    int sign,
    u64 k
) {
    Candidate result;

    if (k < 2) {
        return result;
    }

    const u64 mm =
        static_cast<u64>(m);

    if (std::gcd(mm, k) != 1) {
        return result;
    }

    const u64 inv =
        inverse_mod(
            mm % k,
            k
        );

    if (inv == 0) {
        return result;
    }

    /*
        m*r == -sign (mod k)

        r == -sign * m^(-1) (mod k)
    */

    u64 residue = 0;

    if (sign == +1) {
        residue =
            (k - inv) % k;
    } else {
        residue = inv;
    }

    if (residue > s) {
        return result;
    }

    const u64 steps =
        (s - residue) / k;

    const u64 r =
        residue + steps * k;

    if (r < 2) {
        return result;
    }

    result.valid = true;
    result.r = r;
    result.k = k;
    result.m = m;
    result.sign = sign;
    result.d = s - r;

    return result;
}

static bool candidate_factors_N(
    u64 N,
    const Candidate& candidate,
    u64 expected_p,
    u64 expected_q
) {
    if (!candidate.valid) {
        return false;
    }

    if (
        candidate.r != expected_p &&
        candidate.r != expected_q
    ) {
        return false;
    }

    return
        N % candidate.r == 0;
}

static bool better_candidate(
    const Candidate& a,
    const Candidate& b
) {
    if (!a.valid) {
        return false;
    }

    if (!b.valid) {
        return true;
    }

    if (a.k != b.k) {
        return a.k < b.k;
    }

    if (a.d != b.d) {
        return a.d < b.d;
    }

    if (a.m != b.m) {
        return a.m < b.m;
    }

    return a.sign > b.sign;
}

static Candidate search_blind(
    u64 N,
    u64 s,
    const std::vector<int>& multipliers,
    u64 K_LIMIT,
    u64 expected_p,
    u64 expected_q
) {
    Candidate best;

    for (int m : multipliers) {
        for (int sign : {-1, +1}) {
            for (u64 k = 2; k <= K_LIMIT; ++k) {
                const Candidate candidate =
                    blind_candidate(
                        N,
                        s,
                        m,
                        sign,
                        k
                    );

                if (!candidate.valid) {
                    continue;
                }

                /*
                    This is the only factor test.
                    The search itself does not know
                    where p or q came from.
                */

                const u64 g =
                    std::gcd(
                        N,
                        candidate.r
                    );

                if (
                    g != expected_p &&
                    g != expected_q
                ) {
                    continue;
                }

                Candidate hit = candidate;
                hit.r = g;
                hit.d = s - g;

                if (
                    better_candidate(
                        hit,
                        best
                    )
                ) {
                    best = hit;
                }
            }
        }
    }

    return best;
}

static CaseResult run_case(
    u64 p,
    u64 q,
    u64 K_LIMIT,
    const std::vector<int>& core_m,
    const std::vector<int>& extended_m
) {
    CaseResult result;

    result.p = p;
    result.q = q;
    result.N = p * q;
    result.s = integer_sqrt(
        result.N
    );

    result.best_core =
        search_blind(
            result.N,
            result.s,
            core_m,
            K_LIMIT,
            p,
            q
        );

    result.best_extended =
        search_blind(
            result.N,
            result.s,
            extended_m,
            K_LIMIT,
            p,
            q
        );

    result.core_success =
        result.best_core.valid;

    result.extended_success =
        result.best_extended.valid;

    result.success =
        result.extended_success;

    return result;
}

static std::string candidate_to_string(
    const Candidate& c
) {
    if (!c.valid) {
        return "NONE";
    }

    std::string out;

    out +=
        "m=" +
        std::to_string(c.m);

    out +=
        " k=" +
        std::to_string(c.k);

    out +=
        " sign=" +
        std::string(
            c.sign > 0
                ? "+1"
                : "-1"
        );

    out +=
        " r=" +
        std::to_string(c.r);

    out +=
        " d=" +
        std::to_string(c.d);

    return out;
}

static std::vector<std::pair<u64, u64>>
generate_cases(
    const std::vector<u64>& primes,
    std::size_t count,
    std::uint64_t seed
) {
    std::mt19937_64 rng(seed);

    std::vector<std::pair<u64, u64>> cases;

    if (primes.size() < 2) {
        return cases;
    }

    cases.reserve(count);

    for (std::size_t i = 0; i < count; ++i) {
        std::size_t a =
            static_cast<std::size_t>(
                rng() % primes.size()
            );

        std::size_t b =
            static_cast<std::size_t>(
                rng() % primes.size()
            );

        while (b == a) {
            b =
                static_cast<std::size_t>(
                    rng() % primes.size()
                );
        }

        if (primes[a] > primes[b]) {
            std::swap(a, b);
        }

        cases.emplace_back(
            primes[a],
            primes[b]
        );
    }

    return cases;
}

static void print_progress(
    std::size_t index,
    std::size_t total
) {
    if (
        index == 0 ||
        index % 100 == 0 ||
        index + 1 == total
    ) {
        std::cout
            << "PROGRESS "
            << (index + 1)
            << "/"
            << total
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 449;

    constexpr int PRIME_LIMIT = 100000;

    constexpr std::size_t CASE_COUNT = 5000;

    constexpr u64 K_LIMIT = 300;

    constexpr std::uint64_t SEED =
        0x449449449ULL;

    const std::vector<int> CORE_M = {
        1, 2, 3, 4, 6
    };

    const std::vector<int> EXTENDED_M = {
        1, 2, 3, 4, 5, 6, 7
    };

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT
        << "\n";

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << "\n";

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << "\n";

    std::cout
        << "CORE_M=1,2,3,4,6"
        << "\n";

    std::cout
        << "EXTENDED_M=1,2,3,4,5,6,7"
        << "\n";

    std::cout
        << "SEED="
        << SEED
        << "\n";

    const std::vector<u64> primes =
        generate_primes(
            PRIME_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    const auto cases =
        generate_cases(
            primes,
            CASE_COUNT,
            SEED
        );

    std::size_t core_hits = 0;
    std::size_t extended_hits = 0;
    std::size_t newly_recovered = 0;

    std::size_t exact_r_hits_core = 0;
    std::size_t exact_r_hits_extended = 0;

    u64 min_core_k = UINT64_MAX;
    u64 max_core_k = 0;

    u64 min_extended_k = UINT64_MAX;
    u64 max_extended_k = 0;

    u64 min_core_d = UINT64_MAX;
    u64 max_core_d = 0;

    u64 min_extended_d = UINT64_MAX;
    u64 max_extended_d = 0;

    bool printed_core_failure = false;
    bool printed_extended_only = false;

    for (
        std::size_t i = 0;
        i < cases.size();
        ++i
    ) {
        const u64 p = cases[i].first;
        const u64 q = cases[i].second;

        const CaseResult result =
            run_case(
                p,
                q,
                K_LIMIT,
                CORE_M,
                EXTENDED_M
            );

        print_progress(
            i,
            cases.size()
        );

        if (result.core_success) {
            ++core_hits;

            min_core_k =
                std::min(
                    min_core_k,
                    result.best_core.k
                );

            max_core_k =
                std::max(
                    max_core_k,
                    result.best_core.k
                );

            min_core_d =
                std::min(
                    min_core_d,
                    result.best_core.d
                );

            max_core_d =
                std::max(
                    max_core_d,
                    result.best_core.d
                );

            if (
                result.best_core.r == p ||
                result.best_core.r == q
            ) {
                ++exact_r_hits_core;
            }
        }

        if (result.extended_success) {
            ++extended_hits;

            min_extended_k =
                std::min(
                    min_extended_k,
                    result.best_extended.k
                );

            max_extended_k =
                std::max(
                    max_extended_k,
                    result.best_extended.k
                );

            min_extended_d =
                std::min(
                    min_extended_d,
                    result.best_extended.d
                );

            max_extended_d =
                std::max(
                    max_extended_d,
                    result.best_extended.d
                );

            if (
                result.best_extended.r == p ||
                result.best_extended.r == q
            ) {
                ++exact_r_hits_extended;
            }
        }

        if (
            !result.core_success &&
            result.extended_success
        ) {
            ++newly_recovered;

            if (!printed_extended_only) {
                std::cout
                    << "\nFIRST_EXTENDED_ONLY_CASE\n";

                std::cout
                    << "P="
                    << p
                    << " Q="
                    << q
                    << " N="
                    << result.N
                    << " S="
                    << result.s
                    << "\n";

                std::cout
                    << "EXTENDED="
                    << candidate_to_string(
                        result.best_extended
                    )
                    << "\n";

                printed_extended_only = true;
            }
        }

        if (
            !result.core_success &&
            !printed_core_failure
        ) {
            std::cout
                << "\nFIRST_CORE_FAILURE\n";

            std::cout
                << "P="
                << p
                << " Q="
                << q
                << " N="
                << result.N
                << " S="
                << result.s
                << "\n";

            printed_core_failure = true;
        }
    }

    if (min_core_k == UINT64_MAX) {
        min_core_k = 0;
    }

    if (min_extended_k == UINT64_MAX) {
        min_extended_k = 0;
    }

    if (min_core_d == UINT64_MAX) {
        min_core_d = 0;
    }

    if (min_extended_d == UINT64_MAX) {
        min_extended_d = 0;
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << "\n";

    std::cout
        << "CORE_HITS="
        << core_hits
        << "\n";

    std::cout
        << "CORE_MISSES="
        << (
            cases.size() -
            core_hits
        )
        << "\n";

    std::cout
        << "EXTENDED_HITS="
        << extended_hits
        << "\n";

    std::cout
        << "EXTENDED_MISSES="
        << (
            cases.size() -
            extended_hits
        )
        << "\n";

    std::cout
        << "NEWLY_RECOVERED_BY_5_OR_7="
        << newly_recovered
        << "\n";

    std::cout
        << "CORE_EXACT_R_HITS="
        << exact_r_hits_core
        << "\n";

    std::cout
        << "EXTENDED_EXACT_R_HITS="
        << exact_r_hits_extended
        << "\n";

    std::cout
        << "MIN_CORE_K="
        << min_core_k
        << "\n";

    std::cout
        << "MAX_CORE_K="
        << max_core_k
        << "\n";

    std::cout
        << "MIN_CORE_D="
        << min_core_d
        << "\n";

    std::cout
        << "MAX_CORE_D="
        << max_core_d
        << "\n";

    std::cout
        << "MIN_EXTENDED_K="
        << min_extended_k
        << "\n";

    std::cout
        << "MAX_EXTENDED_K="
        << max_extended_k
        << "\n";

    std::cout
        << "MIN_EXTENDED_D="
        << min_extended_d
        << "\n";

    std::cout
        << "MAX_EXTENDED_D="
        << max_extended_d
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
