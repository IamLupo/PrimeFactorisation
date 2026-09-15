#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
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

        if (x >= lo && x <= hi) {
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

        return {p, q, N, s};
    }
}

static mpz_class F(
    const CaseData& c,
    u64 t
) {
    const mpz_class T = mpz_from_u64(t);
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);

    return T * (N - T * A);
}

static mpz_class C(
    const CaseData& c,
    u64 t
) {
    const mpz_class fm = F(c, t - 1);
    const mpz_class f0 = F(c, t);
    const mpz_class fp = F(c, t + 1);

    return fp * fm - f0 * f0;
}

static bool boolean_hit(
    const CaseData& c,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class value = C(c, t);
    const mpz_class g = gcd(N, value);

    return g > 1 && g < N;
}

/*
    Generate pairwise differences from the first K hit positions.

    We deliberately use only the Boolean hit positions.
*/
static std::vector<u64> build_differences(
    const std::vector<u64>& hits,
    std::size_t K
) {
    const std::size_t n =
        std::min(K, hits.size());

    std::vector<u64> differences;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const u64 d =
                hits[j] - hits[i];

            if (d > 0) {
                differences.push_back(d);
            }
        }
    }

    return differences;
}

/*
    For every pair of differences, calculate their gcd.

    A true p-generated pair of same-class differences gives p
    (or a multiple of p). The count of exact gcd==p is therefore
    a factor-blind periodicity statistic.
*/
static std::map<u64, u64> gcd_pair_frequency(
    const std::vector<u64>& differences
) {
    std::map<u64, u64> frequencies;

    for (std::size_t i = 0;
         i < differences.size();
         ++i) {

        for (std::size_t j = i + 1;
             j < differences.size();
             ++j) {

            const u64 g =
                std::gcd(
                    differences[i],
                    differences[j]
                );

            if (g > 1) {
                ++frequencies[g];
            }
        }
    }

    return frequencies;
}

int main() {
    constexpr int EXPERIMENT = 379;
    constexpr int CASES = 500;

    /*
        Limiting the first K hits keeps the experiment bounded.
        K=40 gives at most 780 pairwise differences and about
        303k gcd calculations per case.
    */
    constexpr std::size_t K = 40;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37920260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_with_ge3_hits = 0;
    u64 cases_with_ge6_hits = 0;

    u64 total_difference_count = 0;
    u64 total_gcd_pair_count = 0;

    u64 total_p_gcd_frequency = 0;
    u64 total_q_gcd_frequency = 0;

    u64 max_p_gcd_frequency = 0;
    u64 max_q_gcd_frequency = 0;

    u64 cases_p_gcd_positive = 0;
    u64 cases_q_gcd_positive = 0;

    u64 cases_p_max_unique = 0;
    u64 cases_p_max_tied = 0;
    u64 cases_wrong_unique_max = 0;
    u64 cases_no_candidate = 0;

    u64 cases_p_in_top3 = 0;

    /*
        Structural validation:
        if three positions belong to the same p residue class,
        gcd of two differences from their anchor should be
        divisible by p.
    */
    u64 p_exact_gcd_hits = 0;
    u64 p_multiple_gcd_hits = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        std::vector<u64> hits;

        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;

            if (boolean_hit(c, t)) {
                hits.push_back(t);
                ++total_boolean_hits;
            }
        }

        if (hits.size() >= 3) {
            ++cases_with_ge3_hits;
        }

        if (hits.size() >= 6) {
            ++cases_with_ge6_hits;
        }

        if (hits.size() < 3) {
            ++cases_no_candidate;
            continue;
        }

        const std::vector<u64> differences =
            build_differences(hits, K);

        total_difference_count +=
            differences.size();

        const std::map<u64, u64> frequencies =
            gcd_pair_frequency(differences);

        u64 max_frequency = 0;
        u64 best_candidate = 0;
        u64 best_count = 0;

        std::vector<
            std::pair<u64, u64>
        > ordered;

        ordered.reserve(frequencies.size());

        for (const auto& entry : frequencies) {
            ordered.push_back(entry);
            total_gcd_pair_count += 1;

            const u64 candidate = entry.first;
            const u64 frequency = entry.second;

            if (frequency > max_frequency) {
                max_frequency = frequency;
                best_candidate = candidate;
                best_count = 1;
            } else if (frequency == max_frequency) {
                ++best_count;
            }
        }

        /*
            Exact p/q diagnostics.
        */
        const auto p_it =
            frequencies.find(p);

        const u64 p_frequency =
            p_it == frequencies.end()
                ? 0
                : p_it->second;

        const auto q_it =
            frequencies.find(q);

        const u64 q_frequency =
            q_it == frequencies.end()
                ? 0
                : q_it->second;

        total_p_gcd_frequency +=
            p_frequency;

        total_q_gcd_frequency +=
            q_frequency;

        max_p_gcd_frequency =
            std::max(
                max_p_gcd_frequency,
                p_frequency
            );

        max_q_gcd_frequency =
            std::max(
                max_q_gcd_frequency,
                q_frequency
            );

        if (p_frequency > 0) {
            ++cases_p_gcd_positive;
        }

        if (q_frequency > 0) {
            ++cases_q_gcd_positive;
        }

        /*
            Is p the unique strongest gcd candidate?
        */
        if (best_count == 1) {
            if (best_candidate == p) {
                ++cases_p_max_unique;
            } else {
                ++cases_wrong_unique_max;
            }
        } else {
            if (p_frequency == max_frequency) {
                ++cases_p_max_tied;
            }
        }

        /*
            Top-3 diagnostic.
        */
        std::sort(
            ordered.begin(),
            ordered.end(),
            [](
                const auto& a,
                const auto& b
            ) {
                if (a.second != b.second) {
                    return a.second > b.second;
                }

                return a.first < b.first;
            }
        );

        const std::size_t top =
            std::min<std::size_t>(
                3,
                ordered.size()
            );

        for (std::size_t i = 0; i < top; ++i) {
            if (ordered[i].first == p) {
                ++cases_p_in_top3;
                break;
            }
        }

        /*
            Direct structural validation using the first K hits.

            Look for triples with:

                d1 = x2-x1
                d2 = x3-x1

            and check whether gcd(d1,d2)=p or a multiple of p.

            We do this only for triples whose three positions are
            known to lie in the same residue class modulo p.
        */
        const std::size_t n =
            std::min(K, hits.size());

        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                for (std::size_t k = j + 1; k < n; ++k) {

                    const u64 r1 =
                        hits[i] % p;

                    const u64 r2 =
                        hits[j] % p;

                    const u64 r3 =
                        hits[k] % p;

                    if (r1 != r2 || r1 != r3) {
                        continue;
                    }

                    const u64 d1 =
                        hits[j] - hits[i];

                    const u64 d2 =
                        hits[k] - hits[i];

                    const u64 g =
                        std::gcd(d1, d2);

                    if (g == p) {
                        ++p_exact_gcd_hits;
                    } else if (g % p == 0) {
                        ++p_multiple_gcd_hits;
                    }
                }
            }
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "K="
        << K
        << '\n';

    std::cout
        << "TOTAL_SCAN_POINTS="
        << total_scan_points
        << '\n';

    std::cout
        << "TOTAL_BOOLEAN_HITS="
        << total_boolean_hits
        << '\n';

    std::cout
        << "CASES_GE3_HITS="
        << cases_with_ge3_hits
        << '\n';

    std::cout
        << "CASES_GE6_HITS="
        << cases_with_ge6_hits
        << '\n';

    std::cout
        << "TOTAL_DIFFERENCES="
        << total_difference_count
        << '\n';

    std::cout
        << "TOTAL_GCD_PAIR_ENTRIES="
        << total_gcd_pair_count
        << '\n';

    std::cout
        << "TOTAL_P_GCD_FREQUENCY="
        << total_p_gcd_frequency
        << '\n';

    std::cout
        << "TOTAL_Q_GCD_FREQUENCY="
        << total_q_gcd_frequency
        << '\n';

    std::cout
        << "MAX_P_GCD_FREQUENCY="
        << max_p_gcd_frequency
        << '\n';

    std::cout
        << "MAX_Q_GCD_FREQUENCY="
        << max_q_gcd_frequency
        << '\n';

    std::cout
        << "CASES_P_GCD_POSITIVE="
        << cases_p_gcd_positive
        << '\n';

    std::cout
        << "CASES_Q_GCD_POSITIVE="
        << cases_q_gcd_positive
        << '\n';

    std::cout
        << "CASES_P_MAX_UNIQUE="
        << cases_p_max_unique
        << '\n';

    std::cout
        << "CASES_P_MAX_TIED="
        << cases_p_max_tied
        << '\n';

    std::cout
        << "CASES_WRONG_UNIQUE_MAX="
        << cases_wrong_unique_max
        << '\n';

    std::cout
        << "CASES_NO_CANDIDATE="
        << cases_no_candidate
        << '\n';

    std::cout
        << "CASES_P_IN_TOP3="
        << cases_p_in_top3
        << '\n';

    std::cout
        << "P_EXACT_GCD_TRIPLES="
        << p_exact_gcd_hits
        << '\n';

    std::cout
        << "P_MULTIPLE_GCD_TRIPLES="
        << p_multiple_gcd_hits
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
