#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
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

struct CoeffStats {
    u64 cases_with_hit = 0;
    u64 total_first_t = 0;
    u64 min_first_t = UINT64_MAX;
    u64 max_first_t = 0;
    u64 total_tests = 0;
    u64 factor_p = 0;
    u64 factor_q = 0;
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

static u64 gcd_candidate(
    const CaseData& c,
    u64 a,
    u64 t
) {
    const mpz_class N = mpz_from_u64(c.N);
    const mpz_class A = mpz_from_u64(a);
    const mpz_class T = mpz_from_u64(t);

    const mpz_class value =
        A * T * T - 1;

    const mpz_class g =
        gcd(N, value);

    return g.get_ui();
}

int main() {
    constexpr int EXPERIMENT = 391;
    constexpr int CASES = 1000;

    constexpr u64 A_MIN = 2;
    constexpr u64 A_MAX = 31;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x39120260914ULL);

    std::vector<CoeffStats> stats(A_MAX + 1);

    u64 sequential_success = 0;
    u64 sequential_failure = 0;

    u64 sequential_total_a = 0;
    u64 sequential_total_t = 0;

    u64 sequential_min_a = UINT64_MAX;
    u64 sequential_max_a = 0;

    u64 sequential_min_t = UINT64_MAX;
    u64 sequential_max_t = 0;

    u64 sequential_factor_p = 0;
    u64 sequential_factor_q = 0;

    u64 best_coefficient_cases = 0;
    u64 coefficient2_best_cases = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c =
            make_case(rng);

        const u64 max_t =
            c.s / 2;

        /*
            Store the first hit for each coefficient for this case.
        */
        std::vector<u64> case_first_t(
            A_MAX + 1,
            0
        );

        /*
            ---------------------------------------------------------
            Individual coefficient statistics.
            ---------------------------------------------------------
        */
        for (u64 a = A_MIN;
             a <= A_MAX;
             ++a) {

            CoeffStats& st = stats[a];

            for (u64 t = 1;
                 t <= max_t;
                 ++t) {

                ++st.total_tests;

                const u64 g =
                    gcd_candidate(c, a, t);

                if (g <= 1 ||
                    g >= c.N) {
                    continue;
                }

                case_first_t[a] = t;

                ++st.cases_with_hit;
                st.total_first_t += t;

                st.min_first_t =
                    std::min(
                        st.min_first_t,
                        t
                    );

                st.max_first_t =
                    std::max(
                        st.max_first_t,
                        t
                    );

                if (g == c.p) {
                    ++st.factor_p;
                } else if (g == c.q) {
                    ++st.factor_q;
                }

                break;
            }
        }

        /*
            ---------------------------------------------------------
            Sequential coefficient strategy.
            ---------------------------------------------------------
        */
        bool sequential_found = false;

        for (u64 a = A_MIN;
             a <= A_MAX;
             ++a) {

            ++sequential_total_a;

            for (u64 t = 1;
                 t <= max_t;
                 ++t) {

                ++sequential_total_t;

                const u64 g =
                    gcd_candidate(c, a, t);

                if (g <= 1 ||
                    g >= c.N) {
                    continue;
                }

                sequential_found = true;

                ++sequential_success;

                sequential_min_a =
                    std::min(
                        sequential_min_a,
                        a
                    );

                sequential_max_a =
                    std::max(
                        sequential_max_a,
                        a
                    );

                sequential_min_t =
                    std::min(
                        sequential_min_t,
                        t
                    );

                sequential_max_t =
                    std::max(
                        sequential_max_t,
                        t
                    );

                if (g == c.p) {
                    ++sequential_factor_p;
                } else if (g == c.q) {
                    ++sequential_factor_q;
                }

                break;
            }

            if (sequential_found) {
                break;
            }
        }

        if (!sequential_found) {
            ++sequential_failure;
        }

        /*
            ---------------------------------------------------------
            Determine which coefficient has the smallest first hit.
            ---------------------------------------------------------
        */
        u64 best_t = UINT64_MAX;
        std::vector<u64> best_as;

        for (u64 a = A_MIN;
             a <= A_MAX;
             ++a) {

            const u64 t =
                case_first_t[a];

            if (t == 0) {
                continue;
            }

            if (t < best_t) {
                best_t = t;
                best_as.clear();
                best_as.push_back(a);
            } else if (t == best_t) {
                best_as.push_back(a);
            }
        }

        if (best_t != UINT64_MAX) {
            ++best_coefficient_cases;

            if (std::find(
                    best_as.begin(),
                    best_as.end(),
                    2
                ) != best_as.end()) {

                ++coefficient2_best_cases;
            }
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "A_MIN="
        << A_MIN
        << '\n';

    std::cout
        << "A_MAX="
        << A_MAX
        << '\n';

    std::cout
        << "SEQUENTIAL_SUCCESS="
        << sequential_success
        << '\n';

    std::cout
        << "SEQUENTIAL_FAILURE="
        << sequential_failure
        << '\n';

    std::cout
        << "SEQUENTIAL_TOTAL_A_TESTS="
        << sequential_total_a
        << '\n';

    std::cout
        << "SEQUENTIAL_TOTAL_T_TESTS="
        << sequential_total_t
        << '\n';

    std::cout
        << "SEQUENTIAL_MIN_A="
        << sequential_min_a
        << '\n';

    std::cout
        << "SEQUENTIAL_MAX_A="
        << sequential_max_a
        << '\n';

    std::cout
        << "SEQUENTIAL_MIN_T="
        << sequential_min_t
        << '\n';

    std::cout
        << "SEQUENTIAL_MAX_T="
        << sequential_max_t
        << '\n';

    std::cout
        << "SEQUENTIAL_FACTOR_P="
        << sequential_factor_p
        << '\n';

    std::cout
        << "SEQUENTIAL_FACTOR_Q="
        << sequential_factor_q
        << '\n';

    for (u64 a = A_MIN;
         a <= A_MAX;
         ++a) {

        const CoeffStats& st = stats[a];

        std::cout
            << "COEFF="
            << a
            << '\n';

        std::cout
            << "  CASES_WITH_HIT="
            << st.cases_with_hit
            << '\n';

        std::cout
            << "  TOTAL_TESTS="
            << st.total_tests
            << '\n';

        std::cout
            << "  FACTOR_P="
            << st.factor_p
            << '\n';

        std::cout
            << "  FACTOR_Q="
            << st.factor_q
            << '\n';

        std::cout
            << "  MIN_T="
            << st.min_first_t
            << '\n';

        std::cout
            << "  MAX_T="
            << st.max_first_t
            << '\n';

        std::cout
            << "  AVERAGE_FIRST_T="
            << (
                st.cases_with_hit > 0
                    ? static_cast<double>(
                        st.total_first_t
                      ) /
                      static_cast<double>(
                        st.cases_with_hit
                      )
                    : 0.0
            )
            << '\n';
    }

    std::cout
        << "BEST_COEFFICIENT_CASES="
        << best_coefficient_cases
        << '\n';

    std::cout
        << "COEFFICIENT_2_IS_BEST="
        << coefficient2_best_cases
        << '\n';

    const bool status =
        sequential_success > 0;

    std::cout
        << "STATUS="
        << (
            status
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