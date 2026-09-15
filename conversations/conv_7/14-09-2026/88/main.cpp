#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
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
    Factor-blind Boolean hit oracle.
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
    For a candidate modulus d, determine the maximum number of hit
    positions that can be covered by TWO residue classes modulo d.

    Since the true p-root set consists of exactly two residue classes
    modulo p, this is the structural statistic we want.

    We deliberately do not require the classes to have any particular
    representatives. We simply count the two most frequent residues.
*/
static u64 best_two_class_coverage(
    const std::vector<u64>& hits,
    u64 d
) {
    if (hits.empty()) {
        return 0;
    }

    std::map<u64, u64> frequency;

    for (u64 x : hits) {
        ++frequency[x % d];
    }

    u64 best1 = 0;
    u64 best2 = 0;

    for (const auto& entry : frequency) {
        const u64 count = entry.second;

        if (count > best1) {
            best2 = best1;
            best1 = count;
        } else if (count > best2) {
            best2 = count;
        }
    }

    return best1 + best2;
}

/*
    Return the number of outliers after selecting the two strongest
    residue classes modulo d.
*/
static u64 outliers_for_modulus(
    const std::vector<u64>& hits,
    u64 d
) {
    const u64 covered =
        best_two_class_coverage(hits, d);

    return hits.size() - covered;
}

int main() {
    constexpr int EXPERIMENT = 377;
    constexpr int CASES = 500;

    /*
        Search candidate moduli only up to max possible p.
        Our generated p is at most 3000.
    */
    constexpr u64 MAX_CANDIDATE = 3000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37720260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_p_exception = 0;

    /*
        Exact structural validation.
    */
    u64 cases_true_p_outliers_le2 = 0;
    u64 cases_true_p_outliers_0 = 0;

    /*
        Recovery statistics.
    */
    u64 cases_p_unique_best = 0;
    u64 cases_p_best_tied = 0;
    u64 cases_wrong_unique_best = 0;
    u64 cases_no_valid_candidate = 0;

    u64 total_best_coverage = 0;
    u64 total_true_p_coverage = 0;

    u64 max_true_p_outliers = 0;
    u64 max_best_outliers = 0;

    /*
        How often does the true p have >= H-2 coverage?
    */
    u64 cases_p_covers_all_but_2 = 0;

    /*
        Diagnostics for hit sets with enough points to contain
        a meaningful two-class structure.
    */
    u64 cases_hits_lt3 = 0;
    u64 cases_hits_ge3 = 0;

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
            Build factor-blind hit set.
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
        } else {
            ++cases_hits_ge3;
        }

        /*
            True-p structural score.

            For a generic case this should leave at most the q-side
            hits as outliers, hence <= 2.
        */
        const u64 true_p_coverage =
            best_two_class_coverage(hits, p);

        const u64 true_p_outliers =
            hits.size() - true_p_coverage;

        total_true_p_coverage +=
            true_p_coverage;

        max_true_p_outliers =
            std::max(
                max_true_p_outliers,
                true_p_outliers
            );

        if (true_p_outliers <= 2) {
            ++cases_true_p_outliers_le2;
        }

        if (true_p_outliers == 0) {
            ++cases_true_p_outliers_0;
        }

        if (hits.size() >= 2 &&
            true_p_coverage + 2 >= hits.size()) {
            ++cases_p_covers_all_but_2;
        }

        /*
            Search all candidate moduli.

            We only care about candidates <= 3000.
            The objective is:

                maximize two-class coverage.

            Ties are retained.
        */
        u64 best_coverage = 0;
        u64 best_candidate = 0;
        u64 number_best = 0;

        if (!hits.empty()) {
            for (u64 d = 2;
                 d <= MAX_CANDIDATE;
                 ++d) {

                const u64 coverage =
                    best_two_class_coverage(
                        hits,
                        d
                    );

                if (coverage > best_coverage) {
                    best_coverage = coverage;
                    best_candidate = d;
                    number_best = 1;
                } else if (coverage == best_coverage) {
                    ++number_best;
                }
            }
        }

        total_best_coverage += best_coverage;

        const u64 best_outliers =
            hits.size() - best_coverage;

        max_best_outliers =
            std::max(
                max_best_outliers,
                best_outliers
            );

        if (hits.empty()) {
            ++cases_no_valid_candidate;
        } else if (number_best == 1) {
            if (best_candidate == p) {
                ++cases_p_unique_best;
            } else {
                ++cases_wrong_unique_best;
            }
        } else {
            /*
                p is counted as a best candidate only when its own
                coverage equals the global optimum.
            */
            if (true_p_coverage == best_coverage) {
                ++cases_p_best_tied;
            } else {
                ++cases_wrong_unique_best;
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
        << "CASES_HITS_LT3="
        << cases_hits_lt3
        << '\n';

    std::cout
        << "CASES_HITS_GE3="
        << cases_hits_ge3
        << '\n';

    std::cout
        << "TRUE_P_TOTAL_COVERAGE="
        << total_true_p_coverage
        << '\n';

    std::cout
        << "CASES_TRUE_P_OUTLIERS_LE2="
        << cases_true_p_outliers_le2
        << '\n';

    std::cout
        << "CASES_TRUE_P_OUTLIERS_0="
        << cases_true_p_outliers_0
        << '\n';

    std::cout
        << "CASES_P_COVERS_ALL_BUT_2="
        << cases_p_covers_all_but_2
        << '\n';

    std::cout
        << "MAX_TRUE_P_OUTLIERS="
        << max_true_p_outliers
        << '\n';

    std::cout
        << "BEST_TOTAL_COVERAGE="
        << total_best_coverage
        << '\n';

    std::cout
        << "MAX_BEST_OUTLIERS="
        << max_best_outliers
        << '\n';

    std::cout
        << "CASES_P_UNIQUE_BEST="
        << cases_p_unique_best
        << '\n';

    std::cout
        << "CASES_P_BEST_TIED="
        << cases_p_best_tied
        << '\n';

    std::cout
        << "CASES_WRONG_UNIQUE_BEST="
        << cases_wrong_unique_best
        << '\n';

    std::cout
        << "CASES_NO_VALID_CANDIDATE="
        << cases_no_valid_candidate
        << '\n';

    const bool structural_status =
        max_true_p_outliers <= 2;

    std::cout
        << "TWO_RESIDUE_CLASS_STATUS="
        << (
            structural_status
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
