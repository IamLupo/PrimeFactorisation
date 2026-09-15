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
    p-root condition:

        2t^2 == 1 (mod p)
*/
static bool p_root(
    u64 t,
    u64 p
) {
    const u64 x = t % p;

    return ((2ULL * x * x) % p) == 1;
}

/*
    Count occurrences of a particular pairwise difference.
*/
static u64 difference_frequency(
    const std::vector<u64>& hits,
    u64 d
) {
    u64 count = 0;

    for (std::size_t i = 0; i < hits.size(); ++i) {
        for (std::size_t j = i + 1; j < hits.size(); ++j) {
            const u64 diff =
                hits[j] - hits[i];

            if (diff == d) {
                ++count;
            }
        }
    }

    return count;
}

int main() {
    constexpr int EXPERIMENT = 376;
    constexpr int CASES = 500;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x37620260914ULL);

    u64 total_scan_points = 0;
    u64 total_boolean_hits = 0;

    u64 cases_p_exception = 0;

    u64 cases_with_p_roots = 0;
    u64 cases_with_ge3_p_roots = 0;

    /*
        Correct structural checks:
        t[i+2]-t[i] must be p.
    */
    u64 p_two_step_failures = 0;
    u64 total_two_step_checks = 0;

    /*
        Pairwise-difference recovery.
    */
    u64 cases_p_recovered = 0;
    u64 cases_wrong_unique_max = 0;
    u64 cases_tied_max = 0;
    u64 cases_no_candidate = 0;

    u64 total_p_frequency = 0;
    u64 total_2p_frequency = 0;

    u64 max_p_frequency = 0;
    u64 max_2p_frequency = 0;

    u64 cases_p_is_max_frequency = 0;

    /*
        Stronger experiment:
        Split p-roots into their two residue classes and determine
        whether each class is an exact arithmetic progression of step p.
    */
    u64 residue_class_failures = 0;
    u64 residue_class_count = 0;

    u64 total_class_points = 0;

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
            Factor-blind hit positions.
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
            Exact p-root positions, only for validating the structure.
        */
        std::vector<u64> p_hits;

        if (!p_exception) {
            for (u64 t = 1; t < s; ++t) {
                if (p_root(t, p)) {
                    p_hits.push_back(t);
                }
            }
        }

        if (!p_hits.empty()) {
            ++cases_with_p_roots;
        }

        if (p_hits.size() >= 3) {
            ++cases_with_ge3_p_roots;

            /*
                Correct relation:
                    p_hits[i+2] - p_hits[i] = p
            */
            for (std::size_t i = 0;
                 i + 2 < p_hits.size();
                 ++i) {

                ++total_two_step_checks;

                if (p_hits[i + 2] - p_hits[i] != p) {
                    ++p_two_step_failures;
                }
            }
        }

        /*
            Separate the two residue classes modulo p.
        */
        if (!p_exception) {
            std::vector<u64> class_a;
            std::vector<u64> class_b;

            for (u64 t : p_hits) {
                const u64 residue = t % p;

                /*
                    Pick the smaller residue as class A.
                    The other root is p-residue.
                */
                if (class_a.empty()) {
                    class_a.push_back(t);
                } else {
                    const u64 first_residue =
                        class_a.front() % p;

                    if (residue == first_residue) {
                        class_a.push_back(t);
                    } else {
                        class_b.push_back(t);
                    }
                }
            }

            /*
                Every class must be an arithmetic progression
                with step p.
            */
            const std::vector<
                std::vector<u64>*
            > classes = {
                &class_a,
                &class_b
            };

            for (const auto* cls : classes) {
                if (cls->empty()) {
                    continue;
                }

                ++residue_class_count;

                total_class_points +=
                    cls->size();

                for (std::size_t i = 1;
                     i < cls->size();
                     ++i) {

                    if ((*cls)[i] - (*cls)[i - 1] != p) {
                        ++residue_class_failures;
                    }
                }
            }
        }

        /*
            Pairwise-difference frequency over Boolean hits.
        */
        if (all_hits.size() >= 3) {
            std::map<u64, u64> frequencies;

            for (std::size_t i = 0;
                 i < all_hits.size();
                 ++i) {

                for (std::size_t j = i + 1;
                     j < all_hits.size();
                     ++j) {

                    const u64 diff =
                        all_hits[j] - all_hits[i];

                    if (diff > 0) {
                        ++frequencies[diff];
                    }
                }
            }

            const u64 freq_p =
                difference_frequency(
                    all_hits,
                    p
                );

            const u64 freq_2p =
                2 * p < s
                    ? difference_frequency(
                        all_hits,
                        2 * p
                    )
                    : 0;

            total_p_frequency += freq_p;
            total_2p_frequency += freq_2p;

            max_p_frequency =
                std::max(
                    max_p_frequency,
                    freq_p
                );

            max_2p_frequency =
                std::max(
                    max_2p_frequency,
                    freq_2p
                );

            /*
                Determine maximum-frequency candidates.
            */
            u64 maximum = 0;
            u64 best = 0;
            u64 number_of_maxima = 0;

            for (const auto& entry : frequencies) {
                const u64 candidate = entry.first;
                const u64 frequency = entry.second;

                if (frequency > maximum) {
                    maximum = frequency;
                    best = candidate;
                    number_of_maxima = 1;
                } else if (frequency == maximum) {
                    ++number_of_maxima;
                }
            }

            if (maximum == 0) {
                ++cases_no_candidate;
            } else if (number_of_maxima == 1) {
                if (best == p) {
                    ++cases_p_recovered;
                    ++cases_p_is_max_frequency;
                } else {
                    ++cases_wrong_unique_max;
                }
            } else {
                ++cases_tied_max;

                /*
                    p is still a maximum if its frequency equals
                    the maximum.
                */
                auto it = frequencies.find(p);

                if (it != frequencies.end() &&
                    it->second == maximum) {
                    ++cases_p_is_max_frequency;
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
        << "CASES_P_EXCEPTION="
        << cases_p_exception
        << '\n';

    std::cout
        << "CASES_WITH_P_ROOTS="
        << cases_with_p_roots
        << '\n';

    std::cout
        << "CASES_WITH_GE3_P_ROOTS="
        << cases_with_ge3_p_roots
        << '\n';

    std::cout
        << "TOTAL_TWO_STEP_CHECKS="
        << total_two_step_checks
        << '\n';

    std::cout
        << "P_TWO_STEP_FAILURES="
        << p_two_step_failures
        << '\n';

    std::cout
        << "RESIDUE_CLASS_COUNT="
        << residue_class_count
        << '\n';

    std::cout
        << "TOTAL_CLASS_POINTS="
        << total_class_points
        << '\n';

    std::cout
        << "RESIDUE_CLASS_FAILURES="
        << residue_class_failures
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
        << "MAX_P_FREQUENCY="
        << max_p_frequency
        << '\n';

    std::cout
        << "MAX_2P_FREQUENCY="
        << max_2p_frequency
        << '\n';

    std::cout
        << "CASES_P_IS_MAX_FREQUENCY="
        << cases_p_is_max_frequency
        << '\n';

    std::cout
        << "CASES_P_RECOVERED="
        << cases_p_recovered
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

    const bool structure_status =
        p_two_step_failures == 0 &&
        residue_class_failures == 0;

    std::cout
        << "P_ROOT_STRUCTURE_STATUS="
        << (structure_status ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
