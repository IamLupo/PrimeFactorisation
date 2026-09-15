#include <algorithm>
#include <chrono>
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

struct Stats {
    u64 gcd_calls = 0;
    u64 product_terms = 0;
    u64 first_root = 0;
    bool found = false;
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

static mpz_class polynomial_value(u64 t) {
    const mpz_class T = mpz_from_u64(t);
    return 2 * T * T - 1;
}

static u64 direct_first_hit(
    const CaseData& c
) {
    const mpz_class N = mpz_from_u64(c.N);

    for (u64 t = 1; t < c.s; ++t) {
        const mpz_class value = polynomial_value(t);
        const mpz_class g = gcd(N, value);

        const u64 gu = g.get_ui();

        if (gu > 1 && gu < c.N) {
            return t;
        }
    }

    return 0;
}

static mpz_class block_product(
    const mpz_class& N,
    u64 lo,
    u64 hi,
    Stats& stats
) {
    mpz_class product = 1;

    for (u64 t = lo; t <= hi; ++t) {
        product *= polynomial_value(t);
        product %= N;
        ++stats.product_terms;
    }

    return product;
}

/*
    Search [lo,hi] for the first nontrivial factor hit.

    The caller only invokes this when the whole interval is known
    to contain a factor hit.

    We always search the left half first, so the first returned
    position is the smallest t in the interval producing a factor.
*/
static u64 isolate_first_hit(
    const mpz_class& N,
    u64 lo,
    u64 hi,
    const mpz_class& known_product,
    Stats& stats
) {
    if (lo == hi) {
        ++stats.gcd_calls;

        const mpz_class value =
            polynomial_value(lo);

        ++stats.product_terms;

        const mpz_class g =
            gcd(N, value);

        const u64 gu = g.get_ui();

        if (gu > 1 && gu < N) {
            return lo;
        }

        return 0;
    }

    const u64 mid =
        lo + (hi - lo) / 2;

    /*
        We cannot derive a GMP child product from the parent product,
        so compute the left child explicitly.
    */
    const mpz_class left_product =
        block_product(
            N,
            lo,
            mid,
            stats
        );

    ++stats.gcd_calls;

    const mpz_class left_g =
        gcd(N, left_product);

    if (left_g > 1) {
        const u64 result =
            isolate_first_hit(
                N,
                lo,
                mid,
                left_product,
                stats
            );

        if (result != 0) {
            return result;
        }
    }

    /*
        No hit on the left, therefore the known parent hit must
        be in the right half.
    */
    const mpz_class right_product =
        block_product(
            N,
            mid + 1,
            hi,
            stats
        );

    ++stats.gcd_calls;

    const mpz_class right_g =
        gcd(N, right_product);

    if (right_g > 1) {
        return isolate_first_hit(
            N,
            mid + 1,
            hi,
            right_product,
            stats
        );
    }

    return 0;
}

static u64 batch_first_hit(
    const CaseData& c,
    u64 block_size,
    Stats& stats
) {
    const mpz_class N =
        mpz_from_u64(c.N);

    for (u64 lo = 1;
         lo < c.s;
         lo += block_size) {

        const u64 hi =
            std::min(
                c.s - 1,
                lo + block_size - 1
            );

        const mpz_class product =
            block_product(
                N,
                lo,
                hi,
                stats
            );

        ++stats.gcd_calls;

        const mpz_class g =
            gcd(N, product);

        if (g == 1) {
            continue;
        }

        const u64 result =
            isolate_first_hit(
                N,
                lo,
                hi,
                product,
                stats
            );

        if (result != 0) {
            return result;
        }
    }

    return 0;
}

int main() {
    constexpr int EXPERIMENT = 389;
    constexpr int CASES = 500;

    const std::vector<u64> BLOCK_SIZES = {
        8,
        16,
        32,
        64,
        128,
        256
    };

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38920260914ULL);

    u64 cases_with_factor = 0;
    u64 cases_without_factor = 0;

    u64 total_direct_first = 0;

    std::vector<u64> total_gcd_calls(
        BLOCK_SIZES.size(),
        0
    );

    std::vector<u64> total_product_terms(
        BLOCK_SIZES.size(),
        0
    );

    std::vector<u64> recovered_cases(
        BLOCK_SIZES.size(),
        0
    );

    std::vector<u64> mismatch_cases(
        BLOCK_SIZES.size(),
        0
    );

    std::vector<u64> no_result_cases(
        BLOCK_SIZES.size(),
        0
    );

    std::vector<double> total_ms(
        BLOCK_SIZES.size(),
        0.0
    );

    u64 max_direct_scan = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c =
            make_case(rng);

        const u64 direct_root =
            direct_first_hit(c);

        if (direct_root == 0) {
            ++cases_without_factor;
            continue;
        }

        ++cases_with_factor;

        total_direct_first +=
            direct_root;

        max_direct_scan =
            std::max(
                max_direct_scan,
                direct_root
            );

        for (std::size_t i = 0;
             i < BLOCK_SIZES.size();
             ++i) {

            const u64 block_size =
                BLOCK_SIZES[i];

            Stats stats;

            const auto start =
                std::chrono::steady_clock::now();

            const u64 recovered =
                batch_first_hit(
                    c,
                    block_size,
                    stats
                );

            const auto finish =
                std::chrono::steady_clock::now();

            const double elapsed_ms =
                std::chrono::duration<double,
                                      std::milli>(
                    finish - start
                ).count();

            total_ms[i] += elapsed_ms;

            total_gcd_calls[i] +=
                stats.gcd_calls;

            total_product_terms[i] +=
                stats.product_terms;

            if (recovered == direct_root &&
                recovered != 0) {
                ++recovered_cases[i];
            }

            if (recovered == 0) {
                ++no_result_cases[i];
            }

            if (recovered != direct_root) {
                ++mismatch_cases[i];
            }
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "CASES_WITH_FACTOR="
        << cases_with_factor
        << '\n';

    std::cout
        << "CASES_WITHOUT_FACTOR="
        << cases_without_factor
        << '\n';

    std::cout
        << "TOTAL_DIRECT_FIRST_ROOT="
        << total_direct_first
        << '\n';

    std::cout
        << "MAX_DIRECT_FIRST_ROOT="
        << max_direct_scan
        << '\n';

    for (std::size_t i = 0;
         i < BLOCK_SIZES.size();
         ++i) {

        const u64 block_size =
            BLOCK_SIZES[i];

        std::cout
            << "BLOCK_SIZE="
            << block_size
            << '\n';

        std::cout
            << "  RECOVERED_CASES="
            << recovered_cases[i]
            << '\n';

        std::cout
            << "  MISMATCH_CASES="
            << mismatch_cases[i]
            << '\n';

        std::cout
            << "  NO_RESULT_CASES="
            << no_result_cases[i]
            << '\n';

        std::cout
            << "  TOTAL_GCD_CALLS="
            << total_gcd_calls[i]
            << '\n';

        std::cout
            << "  TOTAL_PRODUCT_TERMS="
            << total_product_terms[i]
            << '\n';

        std::cout
            << "  AVERAGE_GCD_CALLS="
            << (
                cases_with_factor > 0
                    ? static_cast<double>(
                        total_gcd_calls[i]
                      ) /
                      static_cast<double>(
                        cases_with_factor
                      )
                    : 0.0
            )
            << '\n';

        std::cout
            << "  AVERAGE_PRODUCT_TERMS="
            << (
                cases_with_factor > 0
                    ? static_cast<double>(
                        total_product_terms[i]
                      ) /
                      static_cast<double>(
                        cases_with_factor
                      )
                    : 0.0
            )
            << '\n';

        std::cout
            << "  AVERAGE_TIME_MS="
            << (
                cases_with_factor > 0
                    ? total_ms[i] /
                      static_cast<double>(
                          cases_with_factor
                      )
                    : 0.0
            )
            << '\n';
    }

    bool all_correct = true;

    for (std::size_t i = 0;
         i < BLOCK_SIZES.size();
         ++i) {

        if (mismatch_cases[i] != 0 ||
            no_result_cases[i] != 0) {
            all_correct = false;
        }
    }

    std::cout
        << "FIRST_HIT_BATCH_STATUS="
        << (
            all_correct
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
