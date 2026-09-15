#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;
using i64 = std::int64_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 s;
    u64 D;
    u64 h;
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

static u64 mod_u64(u64 a, u64 m) {
    return a % m;
}

static u64 mod_inverse(u64 a, u64 m) {
    i64 t = 0;
    i64 new_t = 1;
    i64 r = static_cast<i64>(m);
    i64 new_r = static_cast<i64>(a % m);

    while (new_r != 0) {
        const i64 q = r / new_r;

        const i64 next_t = t - q * new_t;
        t = new_t;
        new_t = next_t;

        const i64 next_r = r - q * new_r;
        r = new_r;
        new_r = next_r;
    }

    if (r != 1) {
        return 0;
    }

    if (t < 0) {
        t += static_cast<i64>(m);
    }

    return static_cast<u64>(t);
}

static bool is_prime_u64(u64 n) {
    if (n < 2) {
        return false;
    }

    static const u64 small_primes[] = {
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    };

    for (u64 p : small_primes) {
        if (n == p) {
            return true;
        }

        if (n % p == 0) {
            return false;
        }
    }

    const mpz_class N = mpz_from_u64(n);

    return mpz_probab_prime_p(N.get_mpz_t(), 25) > 0;
}

static bool construct_case(
    u64 p,
    u64 h,
    CaseData& out
) {
    const u64 s = 1 + h * p;

    const u64 s2 = s * s;
    const u64 sp1 = s + 1;
    const u64 upper_square_minus_one =
        sp1 * sp1 - 1;

    const u64 q_min =
        (s2 + p - 1) / p;

    const u64 q_max =
        upper_square_minus_one / p;

    for (u64 q = q_min; q <= q_max; ++q) {
        if (q == p) {
            continue;
        }

        if (q <= p) {
            continue;
        }

        if (!is_prime_u64(q)) {
            continue;
        }

        const u64 N = p * q;
        const u64 actual_s = isqrt_u64(N);

        if (actual_s != s) {
            continue;
        }

        const u64 D = N - s * s;

        out = {p, q, s, D, h};
        return true;
    }

    return false;
}

static mpz_class E_value(
    const CaseData& c,
    u64 z
) {
    const mpz_class S = mpz_from_u64(c.s);
    const mpz_class Z = mpz_from_u64(z);
    const mpz_class D = mpz_from_u64(c.D);

    return (S - Z) *
           (D + S + Z * (S - 1));
}

static bool divisible_by_p2(
    const mpz_class& x,
    u64 p
) {
    const mpz_class P = mpz_from_u64(p);
    const mpz_class P2 = P * P;

    return (x % P2) == 0;
}

static std::set<u64> build_predicted_hits(
    const CaseData& c
) {
    std::set<u64> predicted;

    const u64 p = c.p;
    const u64 s = c.s;
    const u64 h = c.h;

    /*
        First root:
            z == 1 (mod p)
    */
    for (u64 z = 1; z < s; z += p) {
        predicted.insert(z);
    }

    /*
        Exceptional second root:
            h z == h-q (mod p)

        i.e.
            z == 1-q*h^{-1} (mod p)
    */
    const u64 h_mod_p = h % p;

    if (h_mod_p != 0) {
        const u64 inv_h = mod_inverse(h_mod_p, p);

        const u64 q_mod_p = c.q % p;

        const u64 term =
            (q_mod_p * inv_h) % p;

        const u64 root2 =
            (1 + p - term) % p;

        for (u64 z = root2; z < s; z += p) {
            predicted.insert(z);
        }
    }

    return predicted;
}

int main() {
    constexpr int EXPERIMENT = 365;
    constexpr int TARGET_CASES = 100;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    /*
        Small primes deliberately chosen so that

            s = 1 + h p

        can be constructed with q > 2000.
    */
    const std::vector<u64> primes = {
        5, 7, 11, 13, 17, 19, 23,
        29, 31, 37, 41, 43, 47
    };

    std::vector<CaseData> cases;

    /*
        First collect cases with h % p == 0.
        These test whether the "second root" disappears.
    */
    for (u64 p : primes) {
        for (u64 multiplier = 5; multiplier <= 100; ++multiplier) {
            const u64 h = p * multiplier;

            CaseData c{};

            if (construct_case(p, h, c)) {
                cases.push_back(c);

                if (cases.size() >= TARGET_CASES / 2) {
                    break;
                }
            }
        }

        if (cases.size() >= TARGET_CASES / 2) {
            break;
        }
    }

    const std::size_t divisible_count = cases.size();

    /*
        Now collect cases with h % p != 0.
        These should exhibit the second root class.
    */
    for (u64 p : primes) {
        for (u64 h = 20; h <= 300; ++h) {
            if (h % p == 0) {
                continue;
            }

            CaseData c{};

            if (construct_case(p, h, c)) {
                cases.push_back(c);

                if (cases.size() >= TARGET_CASES) {
                    break;
                }
            }
        }

        if (cases.size() >= TARGET_CASES) {
            break;
        }
    }

    u64 factorization_failures = 0;
    u64 predicted_formula_failures = 0;

    u64 total_scan_points = 0;
    u64 total_predicted_hits = 0;
    u64 total_actual_hits = 0;

    u64 cases_with_extra_hits = 0;
    u64 total_extra_hits = 0;

    u64 cases_with_two_root_classes = 0;
    u64 cases_with_one_root_class = 0;

    u64 cases_h_divisible_p = 0;
    u64 cases_h_not_divisible_p = 0;

    u64 mismatched_hit_sets = 0;

    u64 max_actual_hits = 0;
    u64 max_extra_hits = 0;

    for (const CaseData& c : cases) {
        const u64 p = c.p;
        const u64 q = c.q;
        const u64 s = c.s;
        const u64 h = c.h;

        /*
            Check:

                s = 1 + hp
        */
        if (s != 1 + h * p) {
            ++predicted_formula_failures;
        }

        /*
            Check the exact exceptional factorization:

                D+s+z(s-1)
                  =
                p(q-h+hz-h^2p)
        */
        for (u64 z = 0;
             z < std::min<u64>(s, 100);
             ++z) {

            const mpz_class S = mpz_from_u64(s);
            const mpz_class Z = mpz_from_u64(z);
            const mpz_class D = mpz_from_u64(c.D);
            const mpz_class P = mpz_from_u64(p);
            const mpz_class Q = mpz_from_u64(q);
            const mpz_class H = mpz_from_u64(h);

            const mpz_class lhs =
                D + S + Z * (S - 1);

            const mpz_class rhs =
                P * (
                    Q - H
                    + H * Z
                    - H * H * P
                );

            if (lhs != rhs) {
                ++factorization_failures;
            }
        }

        const std::set<u64> predicted =
            build_predicted_hits(c);

        std::set<u64> actual;

        for (u64 z = 0; z < s; ++z) {
            ++total_scan_points;

            const mpz_class value =
                E_value(c, z);

            if (divisible_by_p2(value, p)) {
                actual.insert(z);
            }
        }

        total_predicted_hits += predicted.size();
        total_actual_hits += actual.size();

        max_actual_hits =
            std::max<u64>(
                max_actual_hits,
                actual.size()
            );

        /*
            Count the root-class type.
        */
        if (h % p == 0) {
            ++cases_h_divisible_p;

            if (actual.size() == predicted.size()) {
                ++cases_with_one_root_class;
            }
        } else {
            ++cases_h_not_divisible_p;

            if (predicted.size() > 0) {
                ++cases_with_two_root_classes;
            }
        }

        /*
            Exact set comparison.
        */
        if (actual != predicted) {
            ++mismatched_hit_sets;
        }

        /*
            Extra roots = actual roots not explained
            by the two modular classes.
        */
        for (u64 z : actual) {
            if (predicted.find(z) == predicted.end()) {
                ++total_extra_hits;
            }
        }

        if (actual.size() > predicted.size()) {
            ++cases_with_extra_hits;

            total_extra_hits +=
                actual.size() - predicted.size();
        }

        max_extra_hits =
            std::max<u64>(
                max_extra_hits,
                actual.size() > predicted.size()
                    ? actual.size() - predicted.size()
                    : 0
            );
    }

    std::cout
        << "TOTAL_CASES="
        << cases.size()
        << '\n';

    std::cout
        << "CASES_WITH_H_DIVISIBLE_BY_P="
        << cases_h_divisible_p
        << '\n';

    std::cout
        << "CASES_WITH_H_NOT_DIVISIBLE_BY_P="
        << cases_h_not_divisible_p
        << '\n';

    std::cout
        << "TARGET_FIRST_HALF_CASES="
        << divisible_count
        << '\n';

    std::cout
        << "TOTAL_SCAN_POINTS="
        << total_scan_points
        << '\n';

    std::cout
        << "FACTORIZATION_FAILURES="
        << factorization_failures
        << '\n';

    std::cout
        << "PREDICTED_FORMULA_FAILURES="
        << predicted_formula_failures
        << '\n';

    std::cout
        << "TOTAL_PREDICTED_P2_HITS="
        << total_predicted_hits
        << '\n';

    std::cout
        << "TOTAL_ACTUAL_P2_HITS="
        << total_actual_hits
        << '\n';

    std::cout
        << "MISMATCHED_HIT_SETS="
        << mismatched_hit_sets
        << '\n';

    std::cout
        << "CASES_WITH_EXTRA_HITS="
        << cases_with_extra_hits
        << '\n';

    std::cout
        << "TOTAL_EXTRA_HITS="
        << total_extra_hits
        << '\n';

    std::cout
        << "MAX_ACTUAL_HITS="
        << max_actual_hits
        << '\n';

    std::cout
        << "MAX_EXTRA_HITS="
        << max_extra_hits
        << '\n';

    std::cout
        << "CASES_WITH_EXPECTED_TWO_ROOT_CLASSES="
        << cases_with_two_root_classes
        << '\n';

    std::cout
        << "CASES_WITH_EXPECTED_ONE_ROOT_CLASS="
        << cases_with_one_root_class
        << '\n';

    const bool algebra_ok =
        factorization_failures == 0 &&
        predicted_formula_failures == 0;

    const bool root_structure_ok =
        total_actual_hits == total_predicted_hits &&
        mismatched_hit_sets == 0 &&
        total_extra_hits == 0;

    std::cout
        << "ALGEBRA_STATUS="
        << (algebra_ok ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "EXCEPTIONAL_ROOT_STRUCTURE_STATUS="
        << (root_structure_ok ? "PASS" : "FAIL")
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';

    return 0;
}
