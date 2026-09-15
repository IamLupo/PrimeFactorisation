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
    For one residue class modulo d:

        x_0 < x_1 < ... < x_k

    return the number of points belonging to the longest exact
    arithmetic progression with step d.

    Since all points already have the same residue modulo d, this
    is equivalent to finding the longest run where consecutive
    values differ exactly by d.
*/
static u64 longest_step_run(
    const std::vector<u64>& values,
    u64 d
) {
    if (values.empty()) {
        return 0;
    }

    u64 best = 1;
    u64 current = 1;

    for (std::size_t i = 1; i < values.size(); ++i) {
        if (values[i] - values[i - 1] == d) {
            ++current;
        } else {
            current = 1;
        }

        best = std::max(best, current);
    }

    return best;
}

/*
    Score candidate d.

    For every residue modulo d, collect the hits. For the two
    strongest residue classes, measure the longest exact d-step
    progression in each class.

    Score:
        len1 + len2

    This is deliberately stronger than merely counting how many
    values lie in two residue classes.
*/
static u64 progression_score(
    const std::vector<u64>& hits,
    u64 d
) {
    if (hits.empty()) {
        return 0;
    }

    std::map<u64, std::vector<u64>> classes;

    for (u64 x : hits) {
        classes[x % d].push_back(x);
    }

    u64 best1 = 0;
    u64 best2 = 0;

    for (const auto& entry : classes) {
        const u64 score =
            longest_step_run(entry.second, d);

        if (score > best1) {
            best2 = best1;
            best1 = score;
        } else if (score > best2) {
            best2 = score;
        }
    }

    return best1 + best2;
}

/*
    Stronger exact score:

    For each of the two best residue classes, count ALL points that
    participate in an exact d-step chain.

    This allows isolated q hits to be ignored, but does not count
    arbitrary same-residue points as structured.
*/
static u64 structured_score(
    const std::vector<u64>& hits,
    u64 d
) {
    if (hits.empty()) {
        return 0;
    }

    std::map<u64, std::vector<u64>> classes;

    for (u64 x : hits) {
        classes[x % d].push_back(x);
    }

    std::vector<u64> scores;

    for (const auto& entry : classes) {
        const auto& values = entry.second;

        if (values.empty()) {
            continue;
        }

        u64 score = 1;

        for (std::size_t i = 1; i < values.size(); ++i) {
            if (values[i] - values[i - 1] == d) {
                ++score;
            }
        }

        /*
            The above score can count several disjoint chains as one.
            Recompute using the largest exact chain.
        */
        score = longest_step_run(values, d);

        scores.push_back(score);
    }

    std::sort(
        scores.begin(),
        scores.end(),
        std::greater<u64>()
    );

    u64 result = 0;

    if (!scores.empty()) {
        result += scores[0];
    }

    if (scores.size() >= 2) {
        result += scores[1];
    }

    return result;
}

int main() {
    constexpr int EXPERIMENT = 378;
    constexpr int CASES = 500;
    constexpr u64 MAX_CANDIDATE = 3000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37820260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_p_exception = 0;

    u64 total_true_p_progression_score = 0;
    u64 total_true_p_structured_score = 0;

    u64 cases_true_p_covers_all_but_2 = 0;
    u64 max_true_p_outliers = 0;

    u64 cases_p_unique_best = 0;
    u64 cases_p_tied_best = 0;
    u64 cases_wrong_unique_best = 0;
    u64 cases_no_candidate = 0;

    u64 cases_progression_score_p_unique = 0;
    u64 cases_progression_score_p_tied = 0;
    u64 cases_progression_score_wrong = 0;

    u64 max_progression_score = 0;
    u64 max_structured_score = 0;

    for (int case_id = 0; case_id < CASES; ++case_id) {
        const CaseData c = make_case(rng);

        const u64 p = c.p;
        const u64 s = c.s;

        const bool p_exception =
            ((s - 1) % p == 0);

        if (p_exception) {
            ++cases_p_exception;
        }

        /*
            Factor-blind hit set.
        */
        std::vector<u64> hits;

        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;

            if (boolean_hit(c, t)) {
                hits.push_back(t);
                ++total_boolean_hits;
            }
        }

        /*
            Measure the true p structure.
        */
        const u64 true_progression =
            progression_score(hits, p);

        const u64 true_structured =
            structured_score(hits, p);

        total_true_p_progression_score +=
            true_progression;

        total_true_p_structured_score +=
            true_structured;

        max_progression_score =
            std::max(
                max_progression_score,
                true_progression
            );

        max_structured_score =
            std::max(
                max_structured_score,
                true_structured
            );

        const u64 true_outliers =
            hits.size() > true_structured
                ? hits.size() - true_structured
                : 0;

        max_true_p_outliers =
            std::max(
                max_true_p_outliers,
                true_outliers
            );

        if (true_outliers <= 2) {
            ++cases_true_p_covers_all_but_2;
        }

        /*
            ---------------------------------------------------------
            Search candidate d.
            ---------------------------------------------------------
        */
        u64 best_coverage = 0;
        u64 best_candidate = 0;
        u64 best_count = 0;

        u64 best_progression = 0;
        u64 best_progression_candidate = 0;
        u64 best_progression_count = 0;

        for (u64 d = 2;
             d <= MAX_CANDIDATE;
             ++d) {

            const u64 score =
                structured_score(hits, d);

            if (score > best_coverage) {
                best_coverage = score;
                best_candidate = d;
                best_count = 1;
            } else if (score == best_coverage) {
                ++best_count;
            }

            const u64 progression =
                progression_score(hits, d);

            if (progression > best_progression) {
                best_progression = progression;
                best_progression_candidate = d;
                best_progression_count = 1;
            } else if (progression == best_progression) {
                ++best_progression_count;
            }
        }

        /*
            Exact structured-score recovery.
        */
        if (hits.empty()) {
            ++cases_no_candidate;
        } else if (best_count == 1) {
            if (best_candidate == p) {
                ++cases_p_unique_best;
            } else {
                ++cases_wrong_unique_best;
            }
        } else {
            const u64 p_score =
                structured_score(hits, p);

            if (p_score == best_coverage) {
                ++cases_p_tied_best;
            } else {
                ++cases_wrong_unique_best;
            }
        }

        /*
            Longest-progression-score recovery.
        */
        if (best_progression_count == 1) {
            if (best_progression_candidate == p) {
                ++cases_progression_score_p_unique;
            } else {
                ++cases_progression_score_wrong;
            }
        } else {
            const u64 p_score =
                progression_score(hits, p);

            if (p_score == best_progression) {
                ++cases_progression_score_p_tied;
            } else {
                ++cases_progression_score_wrong;
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
        << "CASES_P_EXCEPTION="
        << cases_p_exception
        << '\n';

    std::cout
        << "TOTAL_TRUE_P_PROGRESSION_SCORE="
        << total_true_p_progression_score
        << '\n';

    std::cout
        << "TOTAL_TRUE_P_STRUCTURED_SCORE="
        << total_true_p_structured_score
        << '\n';

    std::cout
        << "CASES_TRUE_P_COVERS_ALL_BUT_2="
        << cases_true_p_covers_all_but_2
        << '\n';

    std::cout
        << "MAX_TRUE_P_OUTLIERS="
        << max_true_p_outliers
        << '\n';

    std::cout
        << "MAX_PROGRESSION_SCORE="
        << max_progression_score
        << '\n';

    std::cout
        << "MAX_STRUCTURED_SCORE="
        << max_structured_score
        << '\n';

    std::cout
        << "CASES_P_UNIQUE_BEST="
        << cases_p_unique_best
        << '\n';

    std::cout
        << "CASES_P_TIED_BEST="
        << cases_p_tied_best
        << '\n';

    std::cout
        << "CASES_WRONG_UNIQUE_BEST="
        << cases_wrong_unique_best
        << '\n';

    std::cout
        << "CASES_NO_CANDIDATE="
        << cases_no_candidate
        << '\n';

    std::cout
        << "CASES_PROGRESSION_P_UNIQUE="
        << cases_progression_score_p_unique
        << '\n';

    std::cout
        << "CASES_PROGRESSION_P_TIED="
        << cases_progression_score_p_tied
        << '\n';

    std::cout
        << "CASES_PROGRESSION_WRONG="
        << cases_progression_score_wrong
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
