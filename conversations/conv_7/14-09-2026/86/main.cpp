#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
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

/*
    Factor-blind Boolean oracle.
*/
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
    Count how often a particular positive difference d occurs
    among all unordered pairs of hit positions.
*/
static u64 difference_frequency(
    const std::vector<u64>& hits,
    u64 d
) {
    u64 count = 0;

    for (std::size_t i = 0; i < hits.size(); ++i) {
        for (std::size_t j = i + 1; j < hits.size(); ++j) {
            const u64 diff = hits[j] - hits[i];

            if (diff == d) {
                ++count;
            }
        }
    }

    return count;
}

int main() {
    constexpr int EXPERIMENT = 375;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37520260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_p_recovered = 0;
    u64 cases_p_unique_max = 0;
    u64 cases_wrong_unique_max = 0;
    u64 cases_tied_max = 0;
    u64 cases_no_candidate = 0;

    u64 total_p_frequency = 0;
    u64 total_2p_frequency = 0;
    u64 total_q_frequency = 0;

    u64 p_frequency_positive = 0;
    u64 p_frequency_above_2p = 0;

    u64 max_p_frequency = 0;
    u64 max_2p_frequency = 0;
    u64 max_q_frequency = 0;

    /*
        How many cases have enough hit positions for a meaningful
        pairwise-difference analysis.
    */
    u64 cases_hits_lt3 = 0;
    u64 cases_hits_ge3 = 0;

    /*
        Verify that p really appears as a same-root difference when
        there are enough p-root positions.
    */
    u64 spacing_structure_failures = 0;
    u64 cases_verified_p_period = 0;

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

        if (hits.size() < 3) {
            ++cases_hits_lt3;
            continue;
        }

        ++cases_hits_ge3;

        /*
            Difference p:
            Same p-root residue classes are separated by p.
        */
        const u64 freq_p =
            difference_frequency(hits, p);

        const u64 freq_2p =
            (2 * p < s)
                ? difference_frequency(hits, 2 * p)
                : 0;

        /*
            q>s, so there cannot be two distinct q-root positions
            in 1 <= t < s. Therefore q cannot occur as a pairwise
            difference inside the scan interval, but it is reported
            explicitly for completeness.
        */
        const u64 freq_q = 0;

        total_p_frequency += freq_p;
        total_2p_frequency += freq_2p;
        total_q_frequency += freq_q;

        max_p_frequency =
            std::max(max_p_frequency, freq_p);

        max_2p_frequency =
            std::max(max_2p_frequency, freq_2p);

        max_q_frequency =
            std::max(max_q_frequency, freq_q);

        if (freq_p > 0) {
            ++p_frequency_positive;
        }

        if (freq_p > freq_2p) {
            ++p_frequency_above_2p;
        }

        /*
            Build frequency table for all observed differences.
        */
        std::map<u64, u64> frequencies;

        for (std::size_t i = 0; i < hits.size(); ++i) {
            for (std::size_t j = i + 1; j < hits.size(); ++j) {
                const u64 diff =
                    hits[j] - hits[i];

                if (diff == 0) {
                    continue;
                }

                ++frequencies[diff];
            }
        }

        /*
            Find the strongest difference candidate.
        */
        u64 maximum_frequency = 0;
        u64 unique_max_candidates = 0;
        u64 best_candidate = 0;

        for (const auto& entry : frequencies) {
            const u64 candidate = entry.first;
            const u64 frequency = entry.second;

            if (frequency > maximum_frequency) {
                maximum_frequency = frequency;
                unique_max_candidates = 1;
                best_candidate = candidate;
            } else if (frequency == maximum_frequency) {
                ++unique_max_candidates;
            }
        }

        if (maximum_frequency == 0) {
            ++cases_no_candidate;
        } else if (unique_max_candidates == 1) {
            if (best_candidate == p) {
                ++cases_p_recovered;
                ++cases_p_unique_max;
            } else {
                ++cases_wrong_unique_max;
            }
        } else {
            ++cases_tied_max;
        }

        /*
            Structural verification:

            Collect only the p-periodic positions from the hit set.
            Any two consecutive positions from the same root class
            differ by exactly p.
        */
        std::vector<u64> p_hits;

        for (u64 t : hits) {
            /*
                Reconstruct the p-specific condition solely for
                validation of the discovered periodicity.
            */
            const u64 x = t % p;

            if (((2ULL * x * x) % p) == 1) {
                p_hits.push_back(t);
            }
        }

        std::sort(p_hits.begin(), p_hits.end());

        if (p_hits.size() >= 2) {
            ++cases_verified_p_period;

            bool found_p_spacing = false;

            for (std::size_t i = 1;
                 i < p_hits.size();
                 ++i) {

                if (p_hits[i] - p_hits[i - 1] == p) {
                    found_p_spacing = true;
                    break;
                }
            }

            /*
                If there are multiple hits from a single residue
                class within the interval, p must occur as a difference.
            */
            if (p_hits.size() >= 3 && !found_p_spacing) {
                ++spacing_structure_failures;
            }
        }
    }

    std::cout
        << "CASES="
        << CASES
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
        << "CASES_HITS_LT3="
        << cases_hits_lt3
        << '\n';

    std::cout
        << "CASES_HITS_GE3="
        << cases_hits_ge3
        << '\n';

    std::cout
        << "CASES_P_RECOVERED="
        << cases_p_recovered
        << '\n';

    std::cout
        << "CASES_P_UNIQUE_MAX="
        << cases_p_unique_max
        << '\n';

    std::cout
        << "CASES_WRONG_UNIQUE_MAX="
        << cases_wrong_unique_max
        << '\n';

    std::cout
        << "CASES_TIED_MAX="
        << cases_tied_max
        << '\n';

    std::cout
        << "CASES_NO_CANDIDATE="
        << cases_no_candidate
        << '\n';

    std::cout
        << "TOTAL_P_FREQUENCY="
        << total_p_frequency
        << '\n';

    std::cout
        << "TOTAL_2P_FREQUENCY="
        << total_2p_frequency
        << '\n';

    std::cout
        << "TOTAL_Q_FREQUENCY="
        << total_q_frequency
        << '\n';

    std::cout
        << "P_FREQUENCY_POSITIVE="
        << p_frequency_positive
        << '\n';

    std::cout
        << "P_FREQUENCY_ABOVE_2P="
        << p_frequency_above_2p
        << '\n';

    std::cout
        << "MAX_P_FREQUENCY="
        << max_p_frequency
        << '\n';

    std::cout
        << "MAX_2P_FREQUENCY="
        << max_2p_frequency
        << '\n';

    std::cout
        << "MAX_Q_FREQUENCY="
        << max_q_frequency
        << '\n';

    std::cout
        << "CASES_VERIFIED_P_PERIOD="
        << cases_verified_p_period
        << '\n';

    std::cout
        << "SPACING_STRUCTURE_FAILURES="
        << spacing_structure_failures
        << '\n';

    std::cout
        << "STATUS="
        << (
            spacing_structure_failures == 0
                ? "PASS"
                : "FAIL"
        )
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
