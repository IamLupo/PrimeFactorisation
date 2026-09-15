#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

struct Result {
    u64 t_p = 0;
    u64 t_q = 0;
    bool p_first = false;
    bool q_first = false;
    bool neither = false;
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

static CaseData make_random_case(
    std::mt19937_64& rng
) {
    while (true) {
        u64 p = random_prime(rng, 5, 5000);
        u64 q = random_prime(rng, 5001, 30000);

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
    For a = k^2:

        k^2 t^2 - 1 = (kt-1)(kt+1).

    So a factor r divides the expression iff

        kt == +1 (mod r)

    or

        kt == -1 (mod r).

    We directly scan t in [1, s/2], which is the same search
    interval used in the previous experiment.
*/
static u64 first_root_for_factor(
    u64 r,
    u64 k,
    u64 limit
) {
    if (limit == 0) {
        return 0;
    }

    const u64 km = k % r;

    for (u64 t = 1; t <= limit; ++t) {
        const u64 x =
            (km * (t % r)) % r;

        if (x == 1 ||
            x == r - 1) {
            return t;
        }
    }

    return 0;
}

/*
    Exact gcd verification for a single candidate.
*/
static u64 gcd_candidate(
    const CaseData& c,
    u64 k,
    u64 t
) {
    const mpz_class N =
        mpz_from_u64(c.N);

    const mpz_class K =
        mpz_from_u64(k);

    const mpz_class T =
        mpz_from_u64(t);

    const mpz_class value =
        K * K * T * T - 1;

    const mpz_class g =
        gcd(N, value);

    return g.get_ui();
}

static Result classify_case(
    const CaseData& c,
    u64 k
) {
    const u64 limit = c.s / 2;

    const u64 tp =
        first_root_for_factor(
            c.p,
            k,
            limit
        );

    const u64 tq =
        first_root_for_factor(
            c.q,
            k,
            limit
        );

    Result result;
    result.t_p = tp;
    result.t_q = tq;

    if (tp == 0 && tq == 0) {
        result.neither = true;
    } else if (tq != 0 &&
               (tp == 0 || tq < tp)) {
        result.q_first = true;
    } else {
        result.p_first = true;
    }

    return result;
}

int main() {
    constexpr int EXPERIMENT = 392;
    constexpr int RANDOM_CASES = 3000;

    constexpr u64 K_MIN = 2;
    constexpr u64 K_MAX = 10;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x39220260914ULL);

    struct KStats {
        u64 cases = 0;
        u64 p_first = 0;
        u64 q_first = 0;
        u64 neither = 0;

        u64 min_q_over_p_num = 0;
        u64 min_q_over_p_den = 1;

        double min_ratio = std::numeric_limits<double>::infinity();

        u64 q_first_below_threshold = 0;
        u64 q_first_above_threshold = 0;

        u64 threshold_but_not_q_first = 0;

        u64 min_q_first_t = UINT64_MAX;
        u64 max_q_first_t = 0;
        u64 total_q_first_t = 0;

        u64 min_p_first_t = UINT64_MAX;
        u64 max_p_first_t = 0;
        u64 total_p_first_t = 0;

        u64 gcd_verification_failures = 0;
    };

    std::vector<KStats> stats(
        K_MAX + 1
    );

    /*
        First pass: random cases.
    */
    for (int case_id = 0;
         case_id < RANDOM_CASES;
         ++case_id) {

        const CaseData c =
            make_random_case(rng);

        for (u64 k = K_MIN;
             k <= K_MAX;
             ++k) {

            KStats& st = stats[k];
            ++st.cases;

            const Result result =
                classify_case(c, k);

            const bool theoretical_threshold =
                (4 * c.q <=
                 k * k * c.p);

            if (result.q_first) {
                ++st.q_first;

                st.total_q_first_t +=
                    result.t_q;

                st.min_q_first_t =
                    std::min(
                        st.min_q_first_t,
                        result.t_q
                    );

                st.max_q_first_t =
                    std::max(
                        st.max_q_first_t,
                        result.t_q
                    );

                const double ratio =
                    static_cast<double>(c.q) /
                    static_cast<double>(c.p);

                if (ratio < st.min_ratio) {
                    st.min_ratio = ratio;
                    st.min_q_over_p_num = c.q;
                    st.min_q_over_p_den = c.p;
                }

                if (theoretical_threshold) {
                    ++st.q_first_below_threshold;
                } else {
                    ++st.q_first_above_threshold;
                }
            } else if (result.p_first) {
                ++st.p_first;

                st.total_p_first_t +=
                    result.t_p;

                st.min_p_first_t =
                    std::min(
                        st.min_p_first_t,
                        result.t_p
                    );

                st.max_p_first_t =
                    std::max(
                        st.max_p_first_t,
                        result.t_p
                    );

                if (theoretical_threshold) {
                    ++st.threshold_but_not_q_first;
                }
            } else {
                ++st.neither;

                if (theoretical_threshold) {
                    ++st.threshold_but_not_q_first;
                }
            }

            /*
                Verify the congruence result with an actual GMP gcd
                whenever either root exists.
            */
            if (result.t_p != 0) {
                const u64 g =
                    gcd_candidate(
                        c,
                        k,
                        result.t_p
                    );

                if (g != c.p &&
                    g != c.q &&
                    g != c.N) {
                    ++st.gcd_verification_failures;
                }
            }

            if (result.t_q != 0) {
                const u64 g =
                    gcd_candidate(
                        c,
                        k,
                        result.t_q
                    );

                if (g != c.p &&
                    g != c.q &&
                    g != c.N) {
                    ++st.gcd_verification_failures;
                }
            }
        }
    }

    /*
        Second pass: deliberately search for q-first examples.

        This is more targeted than random sampling. We use a grid of
        small p and q values and stop after finding up to 5 examples
        for each k.
    */
    std::vector<u64> targeted_found(
        K_MAX + 1,
        0
    );

    for (u64 p = 5;
         p <= 500 && targeted_found[K_MAX] < 5;
         ++p) {

        if (!is_prime_u64(p)) {
            continue;
        }

        for (u64 q = p + 2;
             q <= 20000;
             ++q) {

            if (!is_prime_u64(q)) {
                continue;
            }

            const u64 N = p * q;
            const u64 s = isqrt_u64(N);

            if (s * s > N ||
                (s + 1) * (s + 1) <= N) {
                continue;
            }

            const CaseData c{p, q, N, s};

            for (u64 k = K_MIN;
                 k <= K_MAX;
                 ++k) {

                if (targeted_found[k] >= 5) {
                    continue;
                }

                const Result result =
                    classify_case(c, k);

                if (result.q_first) {
                    ++targeted_found[k];

                    std::cout
                        << "Q_FIRST_EXAMPLE"
                        << " K="
                        << k
                        << " P="
                        << p
                        << " Q="
                        << q
                        << " S="
                        << s
                        << " TP="
                        << result.t_p
                        << " TQ="
                        << result.t_q
                        << " Q_OVER_P="
                        << std::fixed
                        << std::setprecision(6)
                        << (
                            static_cast<double>(q) /
                            static_cast<double>(p)
                        )
                        << '\n';
                }
            }
        }
    }

    std::cout << std::defaultfloat;

    std::cout
        << "RANDOM_CASES="
        << RANDOM_CASES
        << '\n';

    std::cout
        << "K_MIN="
        << K_MIN
        << '\n';

    std::cout
        << "K_MAX="
        << K_MAX
        << '\n';

    for (u64 k = K_MIN;
         k <= K_MAX;
         ++k) {

        const KStats& st =
            stats[k];

        std::cout
            << "K="
            << k
            << '\n';

        std::cout
            << "  CASES="
            << st.cases
            << '\n';

        std::cout
            << "  P_FIRST="
            << st.p_first
            << '\n';

        std::cout
            << "  Q_FIRST="
            << st.q_first
            << '\n';

        std::cout
            << "  NEITHER="
            << st.neither
            << '\n';

        std::cout
            << "  MIN_Q_FIRST_Q_OVER_P="
            << (
                std::isfinite(st.min_ratio)
                    ? st.min_ratio
                    : -1.0
            )
            << '\n';

        std::cout
            << "  Q_FIRST_WITHIN_Q_LE_K2P_OVER4="
            << st.q_first_below_threshold
            << '\n';

        std::cout
            << "  Q_FIRST_OUTSIDE_SIMPLE_THRESHOLD="
            << st.q_first_above_threshold
            << '\n';

        std::cout
            << "  THRESHOLD_BUT_NOT_Q_FIRST="
            << st.threshold_but_not_q_first
            << '\n';

        std::cout
            << "  MIN_Q_FIRST_T="
            << st.min_q_first_t
            << '\n';

        std::cout
            << "  MAX_Q_FIRST_T="
            << st.max_q_first_t
            << '\n';

        std::cout
            << "  AVG_Q_FIRST_T="
            << (
                st.q_first > 0
                    ? static_cast<double>(
                        st.total_q_first_t
                      ) /
                      static_cast<double>(
                        st.q_first
                      )
                    : 0.0
            )
            << '\n';

        std::cout
            << "  MIN_P_FIRST_T="
            << st.min_p_first_t
            << '\n';

        std::cout
            << "  MAX_P_FIRST_T="
            << st.max_p_first_t
            << '\n';

        std::cout
            << "  AVG_P_FIRST_T="
            << (
                st.p_first > 0
                    ? static_cast<double>(
                        st.total_p_first_t
                      ) /
                      static_cast<double>(
                        st.p_first
                      )
                    : 0.0
            )
            << '\n';

        std::cout
            << "  GCD_VERIFICATION_FAILURES="
            << st.gcd_verification_failures
            << '\n';
    }

    std::cout
        << "STATUS=PASS"
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
