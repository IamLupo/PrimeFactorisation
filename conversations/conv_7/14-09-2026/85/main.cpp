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
    F(t) = t(N-t(s-1))
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
    Return true iff C(t) has a nontrivial factor of N.

    IMPORTANT:
    We deliberately do not inspect which factor was found.
    This is the factor-blind Boolean oracle.
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
    Generic p-root condition.
*/
static bool p_root(
    u64 t,
    u64 p
) {
    const u64 x = t % p;
    return ((2ULL * x * x) % p) == 1;
}

/*
    Given sorted p-hit positions, consecutive gaps should alternate

        a, p-a, a, p-a, ...

    so adjacent gap sums should equal p.
*/
static u64 count_p_gap_sum_matches(
    const std::vector<u64>& hits,
    u64 p
) {
    if (hits.size() < 3) {
        return 0;
    }

    std::vector<u64> gaps;

    gaps.reserve(hits.size() - 1);

    for (std::size_t i = 1; i < hits.size(); ++i) {
        gaps.push_back(hits[i] - hits[i - 1]);
    }

    u64 count = 0;

    for (std::size_t i = 1; i < gaps.size(); ++i) {
        if (gaps[i - 1] + gaps[i] == p) {
            ++count;
        }
    }

    return count;
}

int main() {
    constexpr int EXPERIMENT = 374;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37420260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_p_exception = 0;
    u64 cases_with_no_p_roots = 0;
    u64 cases_with_enough_p_roots = 0;

    u64 gap_identity_failures = 0;

    u64 cases_recovered_p = 0;
    u64 cases_recovered_wrong = 0;
    u64 cases_recovered_none = 0;

    u64 total_candidate_sums = 0;

    u64 max_p_hit_count = 0;
    u64 min_p_hit_count = UINT64_MAX;

    u64 max_q_hit_count = 0;

    u64 total_p_gap_matches = 0;

    /*
        For diagnostics:
        count how often the true p gets the maximum frequency among
        adjacent gap sums in the complete Boolean hit sequence.
    */
    u64 cases_true_p_is_max_frequency = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;

        const bool p_exception =
            ((s - 1) % p == 0);

        if (p_exception) {
            ++cases_p_exception;
        }

        /*
            Factor-blind hit sequence.
        */
        std::vector<u64> all_hits;

        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;

            if (boolean_hit(c, t)) {
                all_hits.push_back(t);
                ++total_boolean_hits;
            }
        }

        /*
            Independently identify the p-root positions only for
            validation of the spacing theorem.
        */
        std::vector<u64> p_hits;

        for (u64 t = 1; t < s; ++t) {
            if (p_root(t, p)) {
                p_hits.push_back(t);
            }
        }

        if (p_exception) {
            /*
                In the exceptional branch p divides every C(t), so
                the two-root arithmetic progression model does not apply.
            */
        } else if (p_hits.empty()) {
            ++cases_with_no_p_roots;
        } else {
            if (p_hits.size() >= 3) {
                ++cases_with_enough_p_roots;
            }

            max_p_hit_count =
                std::max<u64>(
                    max_p_hit_count,
                    p_hits.size()
                );

            min_p_hit_count =
                std::min<u64>(
                    min_p_hit_count,
                    p_hits.size()
                );

            /*
                Verify the spacing identity directly.
            */
            if (p_hits.size() >= 3) {
                for (std::size_t i = 1;
                     i + 1 < p_hits.size();
                     ++i) {

                    const u64 g1 =
                        p_hits[i] - p_hits[i - 1];

                    const u64 g2 =
                        p_hits[i + 1] - p_hits[i];

                    if (g1 + g2 != p) {
                        ++gap_identity_failures;
                    }
                }
            }
        }

        /*
            q has q>s in our construction, so it can contribute at
            most two generic roots in the scanned interval.
        */
        u64 q_hit_count = 0;

        for (u64 t = 1; t < s; ++t) {
            if (p_root(t, q)) {
                ++q_hit_count;
            }
        }

        max_q_hit_count =
            std::max(max_q_hit_count, q_hit_count);

        /*
            ------------------------------------------------------------
            Factor-blind recovery attempt
            ------------------------------------------------------------

            Take the complete Boolean hit sequence and form consecutive
            gaps. For the p-root sequence, adjacent gap sums equal p.

            q can inject at most two anomalous positions.
        */
        if (all_hits.size() >= 3) {
            std::vector<u64> gaps;

            gaps.reserve(all_hits.size() - 1);

            for (std::size_t i = 1;
                 i < all_hits.size();
                 ++i) {
                gaps.push_back(
                    all_hits[i] - all_hits[i - 1]
                );
            }

            std::map<u64, u64> frequency;

            for (std::size_t i = 1;
                 i < gaps.size();
                 ++i) {

                const u64 sum =
                    gaps[i - 1] + gaps[i];

                ++frequency[sum];
                ++total_candidate_sums;
            }

            /*
                Find maximum-frequency candidate.
                Ties are treated as ambiguous.
            */
            u64 max_frequency = 0;
            u64 best_candidate = 0;
            u64 best_count = 0;

            for (const auto& entry : frequency) {
                const u64 candidate = entry.first;
                const u64 count = entry.second;

                if (count > max_frequency) {
                    max_frequency = count;
                    best_candidate = candidate;
                    best_count = 1;
                } else if (count == max_frequency) {
                    ++best_count;
                }
            }

            if (max_frequency > 0) {
                const auto p_it =
                    frequency.find(p);

                const u64 p_frequency =
                    p_it == frequency.end()
                        ? 0
                        : p_it->second;

                if (p_frequency == max_frequency) {
                    ++cases_true_p_is_max_frequency;
                }

                /*
                    Only claim an actual recovery when p is the
                    unique maximum-frequency candidate.
                */
                if (best_count == 1 &&
                    best_candidate == p) {

                    ++cases_recovered_p;
                } else if (best_count == 1) {
                    ++cases_recovered_wrong;
                } else {
                    ++cases_recovered_none;
                }
            } else {
                ++cases_recovered_none;
            }
        } else {
            ++cases_recovered_none;
        }

        /*
            Count the p-gap identity matches within the actual p-only
            sequence.
        */
        if (p_hits.size() >= 3) {
            total_p_gap_matches +=
                count_p_gap_sum_matches(
                    p_hits,
                    p
                );
        }

        /*
            Silence unused warning in some compiler configurations.
        */
        (void)q;
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
        << "CASES_P_EXCEPTION="
        << cases_p_exception
        << '\n';

    std::cout
        << "CASES_WITH_NO_P_ROOTS="
        << cases_with_no_p_roots
        << '\n';

    std::cout
        << "CASES_WITH_ENOUGH_P_ROOTS="
        << cases_with_enough_p_roots
        << '\n';

    std::cout
        << "MIN_P_HIT_COUNT="
        << min_p_hit_count
        << '\n';

    std::cout
        << "MAX_P_HIT_COUNT="
        << max_p_hit_count
        << '\n';

    std::cout
        << "MAX_Q_HIT_COUNT="
        << max_q_hit_count
        << '\n';

    std::cout
        << "P_GAP_IDENTITY_FAILURES="
        << gap_identity_failures
        << '\n';

    std::cout
        << "TOTAL_P_GAP_SUM_MATCHES="
        << total_p_gap_matches
        << '\n';

    std::cout
        << "TOTAL_CANDIDATE_SUMS="
        << total_candidate_sums
        << '\n';

    std::cout
        << "CASES_TRUE_P_IS_MAX_FREQUENCY="
        << cases_true_p_is_max_frequency
        << '\n';

    std::cout
        << "CASES_RECOVERED_P="
        << cases_recovered_p
        << '\n';

    std::cout
        << "CASES_RECOVERED_WRONG="
        << cases_recovered_wrong
        << '\n';

    std::cout
        << "CASES_RECOVERED_NONE="
        << cases_recovered_none
        << '\n';

    const bool spacing_status =
        gap_identity_failures == 0;

    std::cout
        << "P_SPACING_STATUS="
        << (spacing_status ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
