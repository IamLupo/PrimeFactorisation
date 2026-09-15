#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>

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
    const mpz_class N =
        mpz_from_u64(c.N);

    for (u64 t = 1; t < c.s; ++t) {
        const mpz_class value =
            polynomial_value(t);

        const mpz_class g =
            gcd(N, value);

        const u64 gu =
            g.get_ui();

        if (gu > 1 && gu < c.N) {
            return t;
        }
    }

    return 0;
}

static mpz_class product_mod(
    const mpz_class& N,
    u64 lo,
    u64 hi,
    Stats& stats
) {
    mpz_class product = 1;

    if (lo > hi) {
        return product;
    }

    for (u64 t = lo; t <= hi; ++t) {
        product *= polynomial_value(t);
        product %= N;

        ++stats.product_terms;
    }

    return product;
}

/*
    Given a block known to contain a hit, recursively find the
    smallest hit position.
*/
static u64 isolate_first(
    const mpz_class& N,
    u64 lo,
    u64 hi,
    const mpz_class& product,
    Stats& stats
) {
    if (lo == hi) {
        ++stats.gcd_calls;
        ++stats.product_terms;

        const mpz_class value =
            polynomial_value(lo);

        const mpz_class g =
            gcd(N, value);

        return g > 1 ? lo : 0;
    }

    const u64 mid =
        lo + (hi - lo) / 2;

    const mpz_class left_product =
        product_mod(
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
            isolate_first(
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

    const mpz_class right_product =
        product_mod(
            N,
            mid + 1,
            hi,
            stats
        );

    ++stats.gcd_calls;

    const mpz_class right_g =
        gcd(N, right_product);

    if (right_g > 1) {
        return isolate_first(
            N,
            mid + 1,
            hi,
            right_product,
            stats
        );
    }

    return 0;
}

/*
    Adaptive search.

    Start with block size INITIAL_BLOCK. After every miss, double
    the block size. Once a block contains a factor, isolate the first
    position inside that block.

    Example:

        [1,8]
        [9,24]
        [25,56]
        [57,120]
        ...

    The actual covered intervals remain contiguous.
*/
static u64 adaptive_first_hit(
    const CaseData& c,
    u64 initial_block,
    u64 max_block,
    Stats& stats
) {
    const mpz_class N =
        mpz_from_u64(c.N);

    u64 lo = 1;
    u64 block_size = initial_block;

    while (lo < c.s) {
        const u64 hi =
            std::min(
                c.s - 1,
                lo + block_size - 1
            );

        const mpz_class product =
            product_mod(
                N,
                lo,
                hi,
                stats
            );

        ++stats.gcd_calls;

        const mpz_class g =
            gcd(N, product);

        if (g > 1) {
            return isolate_first(
                N,
                lo,
                hi,
                product,
                stats
            );
        }

        lo = hi + 1;

        if (block_size < max_block) {
            const u64 doubled =
                block_size * 2;

            block_size =
                std::min(
                    max_block,
                    doubled
                );
        }
    }

    return 0;
}

int main() {
    constexpr int EXPERIMENT = 390;
    constexpr int CASES = 500;

    const std::vector<u64> INITIAL_BLOCKS = {
        2,
        4,
        8,
        16,
        32
    };

    constexpr u64 MAX_BLOCK = 128;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x39020260914ULL);

    u64 cases_with_factor = 0;
    u64 cases_without_factor = 0;

    std::vector<u64> recovered(
        INITIAL_BLOCKS.size(),
        0
    );

    std::vector<u64> mismatch(
        INITIAL_BLOCKS.size(),
        0
    );

    std::vector<u64> no_result(
        INITIAL_BLOCKS.size(),
        0
    );

    std::vector<u64> total_gcd_calls(
        INITIAL_BLOCKS.size(),
        0
    );

    std::vector<u64> total_product_terms(
        INITIAL_BLOCKS.size(),
        0
    );

    std::vector<double> total_time_ms(
        INITIAL_BLOCKS.size(),
        0.0
    );

    u64 direct_first_total = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c =
            make_case(rng);

        const u64 direct_first =
            direct_first_hit(c);

        if (direct_first == 0) {
            ++cases_without_factor;
            continue;
        }

        ++cases_with_factor;

        direct_first_total +=
            direct_first;

        for (std::size_t i = 0;
             i < INITIAL_BLOCKS.size();
             ++i) {

            const u64 initial =
                INITIAL_BLOCKS[i];

            Stats stats;

            const auto start =
                std::chrono::steady_clock::now();

            const u64 result =
                adaptive_first_hit(
                    c,
                    initial,
                    MAX_BLOCK,
                    stats
                );

            const auto finish =
                std::chrono::steady_clock::now();

            const double elapsed =
                std::chrono::duration<
                    double,
                    std::milli
                >(finish - start).count();

            total_time_ms[i] +=
                elapsed;

            total_gcd_calls[i] +=
                stats.gcd_calls;

            total_product_terms[i] +=
                stats.product_terms;

            if (result == direct_first) {
                ++recovered[i];
            } else {
                ++mismatch[i];
            }

            if (result == 0) {
                ++no_result[i];
            }
        }
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "MAX_BLOCK="
        << MAX_BLOCK
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
        << direct_first_total
        << '\n';

    for (std::size_t i = 0;
         i < INITIAL_BLOCKS.size();
         ++i) {

        std::cout
            << "INITIAL_BLOCK="
            << INITIAL_BLOCKS[i]
            << '\n';

        std::cout
            << "  RECOVERED_CASES="
            << recovered[i]
            << '\n';

        std::cout
            << "  MISMATCH_CASES="
            << mismatch[i]
            << '\n';

        std::cout
            << "  NO_RESULT_CASES="
            << no_result[i]
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
                    ? total_time_ms[i] /
                      static_cast<double>(
                          cases_with_factor
                      )
                    : 0.0
            )
            << '\n';
    }

    bool all_correct = true;

    for (std::size_t i = 0;
         i < INITIAL_BLOCKS.size();
         ++i) {

        if (mismatch[i] != 0 ||
            no_result[i] != 0) {

            all_correct = false;
        }
    }

    std::cout
        << "ADAPTIVE_BATCH_STATUS="
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
