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
    Compute

        2t^2 - 1

    exactly.
*/
static mpz_class polynomial_value(
    u64 t
) {
    const mpz_class T =
        mpz_from_u64(t);

    return 2 * T * T - 1;
}

/*
    Build the product

        product_{t=lo}^{hi} (2t^2-1) mod N

    so the intermediate product remains bounded.
*/
static mpz_class block_product_mod_n(
    const mpz_class& N,
    u64 lo,
    u64 hi
) {
    mpz_class product = 1;

    for (u64 t = lo; t <= hi; ++t) {
        const mpz_class value =
            polynomial_value(t);

        product *= value;
        product %= N;
    }

    return product;
}

/*
    Direct single-value gcd.
*/
static u64 direct_gcd(
    const CaseData& c,
    u64 t
) {
    const mpz_class N =
        mpz_from_u64(c.N);

    const mpz_class value =
        polynomial_value(t);

    const mpz_class g =
        gcd(N, value);

    return g.get_ui();
}

/*
    Recursive isolation of a block containing at least one root.

    The function assumes gcd(N, product(block)) > 1.

    It recursively splits until single positions are reached.

    We count:
      - number of gcd calls
      - number of product terms evaluated
*/
static void isolate_block(
    const mpz_class& N,
    u64 lo,
    u64 hi,
    std::vector<u64>& roots,
    u64& gcd_calls,
    u64& product_terms
) {
    if (lo > hi) {
        return;
    }

    ++gcd_calls;

    if (lo == hi) {
        const mpz_class value =
            polynomial_value(lo);

        ++product_terms;

        const mpz_class g =
            gcd(N, value);

        const u64 gu =
            g.get_ui();

        if (gu > 1) {
            roots.push_back(lo);
        }

        return;
    }

    /*
        Compute the block product modulo N.
    */
    const mpz_class product =
        block_product_mod_n(
            N,
            lo,
            hi
        );

    if (product == 0) {
        /*
            Product is 0 mod N. There is definitely at least one
            factor hit, but we still split normally.
        */
    }

    const mpz_class g =
        gcd(N, product);

    if (g == 1) {
        return;
    }

    const u64 mid =
        lo + (hi - lo) / 2;

    isolate_block(
        N,
        lo,
        mid,
        roots,
        gcd_calls,
        product_terms
    );

    isolate_block(
        N,
        mid + 1,
        hi,
        roots,
        gcd_calls,
        product_terms
    );
}

/*
    More efficient recursive version which receives a block product
    so that it does not recompute the parent product.
*/
static void isolate_block_product(
    const mpz_class& N,
    u64 lo,
    u64 hi,
    const mpz_class& product,
    std::vector<u64>& roots,
    u64& gcd_calls,
    u64& product_terms
) {
    ++gcd_calls;

    if (lo == hi) {
        const mpz_class value =
            polynomial_value(lo);

        ++product_terms;

        const mpz_class g =
            gcd(N, value);

        if (g > 1) {
            roots.push_back(lo);
        }

        return;
    }

    const mpz_class g =
        gcd(N, product);

    if (g == 1) {
        return;
    }

    const u64 mid =
        lo + (hi - lo) / 2;

    const mpz_class left_product =
        block_product_mod_n(
            N,
            lo,
            mid
        );

    const mpz_class right_product =
        block_product_mod_n(
            N,
            mid + 1,
            hi
        );

    isolate_block_product(
        N,
        lo,
        mid,
        left_product,
        roots,
        gcd_calls,
        product_terms
    );

    isolate_block_product(
        N,
        mid + 1,
        hi,
        right_product,
        roots,
        gcd_calls,
        product_terms
    );
}

int main() {
    constexpr int EXPERIMENT = 388;
    constexpr int CASES = 500;

    /*
        Batch size.

        Larger blocks reduce top-level gcd calls but make block
        product construction more expensive. 64 is a reasonable
        baseline for this experiment.
    */
    constexpr u64 BLOCK_SIZE = 64;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x38820260914ULL);

    u64 total_scan_points = 0;
    u64 total_direct_gcd_calls = 0;

    u64 total_batch_gcd_calls = 0;
    u64 total_batch_product_terms = 0;

    u64 total_roots_direct = 0;
    u64 total_roots_batch = 0;

    u64 root_set_mismatches = 0;

    u64 cases_with_root = 0;
    u64 cases_without_root = 0;

    u64 batch_cases_with_root = 0;
    u64 batch_cases_without_root = 0;

    u64 first_root_direct_total = 0;
    u64 first_root_batch_total = 0;

    u64 first_root_mismatches = 0;

    u64 max_direct_gcd_calls = 0;
    u64 max_batch_gcd_calls = 0;

    u64 max_batch_product_terms = 0;

    for (int case_id = 0;
         case_id < CASES;
         ++case_id) {

        const CaseData c =
            make_case(rng);

        const u64 s = c.s;

        /*
            ---------------------------------------------------------
            Baseline: direct gcd for every t.
            ---------------------------------------------------------
        */
        std::vector<u64> direct_roots;

        for (u64 t = 1; t < s; ++t) {
            ++total_scan_points;
            ++total_direct_gcd_calls;

            const u64 g =
                direct_gcd(c, t);

            if (g > 1 &&
                g < c.N) {

                direct_roots.push_back(t);
            }
        }

        if (direct_roots.empty()) {
            ++cases_without_root;
        } else {
            ++cases_with_root;

            first_root_direct_total +=
                direct_roots.front();
        }

        /*
            ---------------------------------------------------------
            Batch GCD.
            ---------------------------------------------------------
        */
        std::vector<u64> batch_roots;

        u64 case_batch_gcd_calls = 0;
        u64 case_batch_product_terms = 0;

        /*
            Process the scan interval in blocks.
        */
        for (u64 lo = 1;
             lo < s;
             lo += BLOCK_SIZE) {

            const u64 hi =
                std::min(
                    s - 1,
                    lo + BLOCK_SIZE - 1
                );

            const mpz_class N =
                mpz_from_u64(c.N);

            const mpz_class product =
                block_product_mod_n(
                    N,
                    lo,
                    hi
                );

            case_batch_product_terms +=
                hi - lo + 1;

            ++case_batch_gcd_calls;

            const mpz_class block_g =
                gcd(N, product);

            if (block_g == 1) {
                continue;
            }

            /*
                Isolate the exact root positions inside this block.
            */
            isolate_block_product(
                N,
                lo,
                hi,
                product,
                batch_roots,
                case_batch_gcd_calls,
                case_batch_product_terms
            );
        }

        std::sort(
            batch_roots.begin(),
            batch_roots.end()
        );

        std::sort(
            direct_roots.begin(),
            direct_roots.end()
        );

        /*
            Remove duplicates defensively.
        */
        direct_roots.erase(
            std::unique(
                direct_roots.begin(),
                direct_roots.end()
            ),
            direct_roots.end()
        );

        batch_roots.erase(
            std::unique(
                batch_roots.begin(),
                batch_roots.end()
            ),
            batch_roots.end()
        );

        total_roots_direct +=
            direct_roots.size();

        total_roots_batch +=
            batch_roots.size();

        if (batch_roots != direct_roots) {
            ++root_set_mismatches;
        }

        if (batch_roots.empty()) {
            ++batch_cases_without_root;
        } else {
            ++batch_cases_with_root;

            first_root_batch_total +=
                batch_roots.front();
        }

        if (!direct_roots.empty() &&
            !batch_roots.empty() &&
            direct_roots.front() !=
                batch_roots.front()) {

            ++first_root_mismatches;
        }

        total_batch_gcd_calls +=
            case_batch_gcd_calls;

        total_batch_product_terms +=
            case_batch_product_terms;

        max_direct_gcd_calls =
            std::max(
                max_direct_gcd_calls,
                s > 0 ? s - 1 : 0
            );

        max_batch_gcd_calls =
            std::max(
                max_batch_gcd_calls,
                case_batch_gcd_calls
            );

        max_batch_product_terms =
            std::max(
                max_batch_product_terms,
                case_batch_product_terms
            );
    }

    std::cout
        << "CASES="
        << CASES
        << '\n';

    std::cout
        << "BLOCK_SIZE="
        << BLOCK_SIZE
        << '\n';

    std::cout
        << "TOTAL_SCAN_POINTS="
        << total_scan_points
        << '\n';

    std::cout
        << "TOTAL_DIRECT_GCD_CALLS="
        << total_direct_gcd_calls
        << '\n';

    std::cout
        << "TOTAL_BATCH_GCD_CALLS="
        << total_batch_gcd_calls
        << '\n';

    std::cout
        << "TOTAL_BATCH_PRODUCT_TERMS="
        << total_batch_product_terms
        << '\n';

    std::cout
        << "TOTAL_DIRECT_ROOTS="
        << total_roots_direct
        << '\n';

    std::cout
        << "TOTAL_BATCH_ROOTS="
        << total_roots_batch
        << '\n';

    std::cout
        << "ROOT_SET_MISMATCHES="
        << root_set_mismatches
        << '\n';

    std::cout
        << "CASES_WITH_ROOT="
        << cases_with_root
        << '\n';

    std::cout
        << "CASES_WITHOUT_ROOT="
        << cases_without_root
        << '\n';

    std::cout
        << "BATCH_CASES_WITH_ROOT="
        << batch_cases_with_root
        << '\n';

    std::cout
        << "BATCH_CASES_WITHOUT_ROOT="
        << batch_cases_without_root
        << '\n';

    std::cout
        << "TOTAL_FIRST_ROOT_DIRECT="
        << first_root_direct_total
        << '\n';

    std::cout
        << "TOTAL_FIRST_ROOT_BATCH="
        << first_root_batch_total
        << '\n';

    std::cout
        << "FIRST_ROOT_MISMATCHES="
        << first_root_mismatches
        << '\n';

    std::cout
        << "MAX_DIRECT_GCD_CALLS="
        << max_direct_gcd_calls
        << '\n';

    std::cout
        << "MAX_BATCH_GCD_CALLS="
        << max_batch_gcd_calls
        << '\n';

    std::cout
        << "MAX_BATCH_PRODUCT_TERMS="
        << max_batch_product_terms
        << '\n';

    const bool status =
        root_set_mismatches == 0 &&
        first_root_mismatches == 0 &&
        cases_with_root == batch_cases_with_root &&
        cases_without_root == batch_cases_without_root;

    std::cout
        << "BATCH_GCD_CORRECTNESS_STATUS="
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
