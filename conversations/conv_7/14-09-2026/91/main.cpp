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

/*
    F(t) = t(N - t(s-1))
*/
static mpz_class F(
    const CaseData& c,
    u64 t
) {
    const mpz_class T = mpz_from_u64(t);
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(c.s - 1);

    return T * (N - T * A);
}

/*
    C(t)=F(t+1)F(t-1)-F(t)^2
*/
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
    Count the frequency of exact lag-2 differences

        h[i+2] - h[i].
*/
static std::map<u64, u64> lag2_frequencies(
    const std::vector<u64>& hits
) {
    std::map<u64, u64> frequencies;

    if (hits.size() < 3) {
        return frequencies;
    }

    for (std::size_t i = 0;
         i + 2 < hits.size();
         ++i) {

        const u64 d =
            hits[i + 2] - hits[i];

        ++frequencies[d];
    }

    return frequencies;
}

int main() {
    constexpr int EXPERIMENT = 380;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38020260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_hits_lt3 = 0;
    u64 cases_hits_ge3 = 0;

    u64 total_lag2_points = 0;

    u64 total_true_p_lag2 = 0;
    u64 max_true_p_lag2 = 0;

    u64 cases_p_frequency_positive = 0;

    /*
        Recovery by mode of lag-2 difference.
    */
    u64 cases_unique_p_mode = 0;
    u64 cases_p_tied_mode = 0;
    u64 cases_wrong_unique_mode = 0;
    u64 cases_no_mode = 0;

    /*
        How often does p equal the global maximum frequency?
    */
    u64 cases_p_is_max_frequency = 0;

    /*
        Exact structural validation using the true p.
    */
    u64 two_step_p_failures = 0;
    u64 two_step_p_checks = 0;

    /*
        q can contribute at most two hits because q>s.
        This gives us a theoretical contamination bound.
    */
    u64 cases_outliers_le4 = 0;
    u64 max_non_p_lag2 = 0;

    /*
        Compare against lag-1 pair sums from the previous experiment.
    */
    u64 total_lag1_pair_sums = 0;
    u64 p_lag1_frequency = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 s = c.s;

        /*
            Build Boolean hit sequence.
        */
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
            ++cases_no_mode;
            continue;
        }

        ++cases_hits_ge3;

        const std::map<u64, u64> frequencies =
            lag2_frequencies(hits);

        total_lag2_points +=
            hits.size() - 2;

        /*
            Frequency at the true p.
        */
        const auto p_it =
            frequencies.find(p);

        const u64 p_frequency =
            p_it == frequencies.end()
                ? 0
                : p_it->second;

        total_true_p_lag2 +=
            p_frequency;

        max_true_p_lag2 =
            std::max(
                max_true_p_lag2,
                p_frequency
            );

        if (p_frequency > 0) {
            ++cases_p_frequency_positive;
        }

        /*
            Find the frequency mode.
        */
        u64 maximum_frequency = 0;
        u64 best_candidate = 0;
        u64 number_of_maxima = 0;

        for (const auto& entry : frequencies) {
            const u64 candidate = entry.first;
            const u64 frequency = entry.second;

            if (frequency > maximum_frequency) {
                maximum_frequency = frequency;
                best_candidate = candidate;
                number_of_maxima = 1;
            } else if (frequency == maximum_frequency) {
                ++number_of_maxima;
            }
        }

        if (maximum_frequency == 0) {
            ++cases_no_mode;
        } else if (number_of_maxima == 1) {
            if (best_candidate == p) {
                ++cases_unique_p_mode;
            } else {
                ++cases_wrong_unique_mode;
            }
        } else {
            if (p_frequency == maximum_frequency) {
                ++cases_p_tied_mode;
            } else {
                ++cases_wrong_unique_mode;
            }
        }

        if (p_frequency == maximum_frequency) {
            ++cases_p_is_max_frequency;
        }

        /*
            ---------------------------------------------------------
            Exact validation of the two-step p structure.
            ---------------------------------------------------------

            Find the actual p-root positions, then verify

                t[i+2] - t[i] = p.
        */
        std::vector<u64> p_hits;

        if ((s - 1) % p != 0) {
            for (u64 t = 1; t < s; ++t) {
                const u64 x = t % p;

                if (((2ULL * x * x) % p) == 1) {
                    p_hits.push_back(t);
                }
            }
        }

        if (p_hits.size() >= 3) {
            for (std::size_t i = 0;
                 i + 2 < p_hits.size();
                 ++i) {

                ++two_step_p_checks;

                if (p_hits[i + 2] -
                    p_hits[i] != p) {

                    ++two_step_p_failures;
                }
            }
        }

        /*
            The Boolean hit sequence contains all p-hits plus at most
            two q-hits. Each additional hit can affect at most two
            lag-2 windows on each side, giving a small contamination
            bound. We record the number of lag-2 values != p.
        */
        u64 non_p_lag2 = 0;

        for (const auto& entry : frequencies) {
            if (entry.first != p) {
                non_p_lag2 += entry.second;
            }
        }

        max_non_p_lag2 =
            std::max(
                max_non_p_lag2,
                non_p_lag2
            );

        if (non_p_lag2 <= 4) {
            ++cases_outliers_le4;
        }

        /*
            Previous-style statistic:
            adjacent gap sums.
        */
        if (hits.size() >= 3) {
            std::vector<u64> gaps;

            gaps.reserve(hits.size() - 1);

            for (std::size_t i = 1;
                 i < hits.size();
                 ++i) {
                gaps.push_back(
                    hits[i] - hits[i - 1]
                );
            }

            for (std::size_t i = 1;
                 i < gaps.size();
                 ++i) {

                ++total_lag1_pair_sums;

                if (gaps[i - 1] + gaps[i] == p) {
                    ++p_lag1_frequency;
                }
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
        << "TOTAL_LAG2_POINTS="
        << total_lag2_points
        << '\n';

    std::cout
        << "TOTAL_TRUE_P_LAG2_FREQUENCY="
        << total_true_p_lag2
        << '\n';

    std::cout
        << "MAX_TRUE_P_LAG2_FREQUENCY="
        << max_true_p_lag2
        << '\n';

    std::cout
        << "CASES_P_FREQUENCY_POSITIVE="
        << cases_p_frequency_positive
        << '\n';

    std::cout
        << "CASES_UNIQUE_P_MODE="
        << cases_unique_p_mode
        << '\n';

    std::cout
        << "CASES_P_TIED_MODE="
        << cases_p_tied_mode
        << '\n';

    std::cout
        << "CASES_P_IS_MAX_FREQUENCY="
        << cases_p_is_max_frequency
        << '\n';

    std::cout
        << "CASES_WRONG_UNIQUE_MODE="
        << cases_wrong_unique_mode
        << '\n';

    std::cout
        << "CASES_NO_MODE="
        << cases_no_mode
        << '\n';

    std::cout
        << "TWO_STEP_P_CHECKS="
        << two_step_p_checks
        << '\n';

    std::cout
        << "TWO_STEP_P_FAILURES="
        << two_step_p_failures
        << '\n';

    std::cout
        << "CASES_OUTLIERS_LE4="
        << cases_outliers_le4
        << '\n';

    std::cout
        << "MAX_NON_P_LAG2="
        << max_non_p_lag2
        << '\n';

    std::cout
        << "TOTAL_LAG1_PAIR_SUMS="
        << total_lag1_pair_sums
        << '\n';

    std::cout
        << "P_LAG1_FREQUENCY="
        << p_lag1_frequency
        << '\n';

    const bool structure_status =
        two_step_p_failures == 0;

    std::cout
        << "P_TWO_STEP_STRUCTURE_STATUS="
        << (
            structure_status
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
