#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 465;

constexpr int PRIME_LIMIT = 100000;
constexpr int PRIME_MIN = 10000;
constexpr int CASE_COUNT = 500;

constexpr int K_LIMIT = 1000;
constexpr int M_LIMIT = 7;

struct PrimeCase {
    u64 p;
    u64 q;
    u64 n;
    u64 s;
};

struct Representation {
    bool found = false;

    u64 k1 = 0;
    u64 k2 = 0;
    u64 m = 0;
    u64 j = 0;
    u64 base = 0;

    /*
     * Identity:
     *
     *   (j*k1 + 1)(j*k2 + 1) = j*p + 1
     */
    u64 left_factor = 0;
    u64 right_factor = 0;
};

struct CRTPair {
    int k1;
    int k2;
};

std::vector<int> generate_primes(int limit) {
    std::vector<bool> composite(limit + 1, false);
    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<std::int64_t>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

std::vector<PrimeCase> generate_cases(
    const std::vector<int>& primes,
    std::mt19937_64& rng
) {
    std::vector<int> candidates;

    for (int p : primes) {
        if (p >= PRIME_MIN && p <= PRIME_LIMIT) {
            candidates.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        candidates.size() - 1
    );

    std::vector<PrimeCase> cases;
    cases.reserve(CASE_COUNT);

    for (int i = 0; i < CASE_COUNT; ++i) {
        int p = candidates[dist(rng)];
        int q = candidates[dist(rng)];

        while (q == p) {
            q = candidates[dist(rng)];
        }

        if (p > q) {
            std::swap(p, q);
        }

        PrimeCase c;
        c.p = static_cast<u64>(p);
        c.q = static_cast<u64>(q);
        c.n = c.p * c.q;

        c.s = static_cast<u64>(
            std::sqrt(static_cast<long double>(c.n))
        );

        while ((c.s + 1) * (c.s + 1) <= c.n) {
            ++c.s;
        }

        while (c.s * c.s > c.n) {
            --c.s;
        }

        cases.push_back(c);
    }

    return cases;
}

std::vector<CRTPair> build_determinant_one_pairs() {
    std::vector<CRTPair> pairs;

    /*
     * Determinant-one condition:
     *
     *     m2*k1 - m1*k2 = 1
     *
     * for some 1 <= m1,m2 <= M_LIMIT.
     *
     * Only retain the (k1,k2) pair.
     */
    for (int k1 = 1; k1 <= K_LIMIT; ++k1) {
        for (int k2 = k1; k2 <= K_LIMIT; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            bool valid = false;

            for (int m1 = 1; m1 <= M_LIMIT; ++m1) {
                const int numerator = 1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const int m2 = numerator / k1;

                if (m2 >= 1 && m2 <= M_LIMIT) {
                    valid = true;
                    break;
                }
            }

            if (valid) {
                pairs.push_back({k1, k2});
            }
        }
    }

    return pairs;
}

Representation find_min_j_representation(
    u64 prime,
    const std::vector<CRTPair>& pairs
) {
    Representation best;

    for (const CRTPair& pair : pairs) {
        const u64 k1 = static_cast<u64>(pair.k1);
        const u64 k2 = static_cast<u64>(pair.k2);

        const u64 base = k1 + k2;
        const u64 m = k1 * k2;

        if (base > prime) {
            continue;
        }

        const u64 remainder = prime - base;

        if (remainder % m != 0) {
            continue;
        }

        const u64 j = remainder / m;

        /*
         * Smaller j is primary.
         *
         * Secondary:
         *   smaller M
         *   smaller max(k)
         *   smaller k1
         */
        bool take = false;

        if (!best.found) {
            take = true;
        } else if (j < best.j) {
            take = true;
        } else if (j == best.j && m < best.m) {
            take = true;
        } else if (j == best.j &&
                   m == best.m &&
                   std::max(k1, k2) <
                       std::max(best.k1, best.k2)) {
            take = true;
        } else if (j == best.j &&
                   m == best.m &&
                   std::max(k1, k2) ==
                       std::max(best.k1, best.k2) &&
                   k1 < best.k1) {
            take = true;
        }

        if (!take) {
            continue;
        }

        best.found = true;
        best.k1 = k1;
        best.k2 = k2;
        best.m = m;
        best.j = j;
        best.base = base;

        best.left_factor = j * k1 + 1;
        best.right_factor = j * k2 + 1;
    }

    return best;
}

void print_representation(
    const char* name,
    u64 prime,
    const Representation& r
) {
    std::cout
        << name
        << "_PRIME="
        << prime
        << '\n';

    std::cout
        << name
        << "_FOUND="
        << (r.found ? 1 : 0)
        << '\n';

    if (!r.found) {
        return;
    }

    std::cout
        << name
        << "_K1="
        << r.k1
        << '\n';

    std::cout
        << name
        << "_K2="
        << r.k2
        << '\n';

    std::cout
        << name
        << "_M="
        << r.m
        << '\n';

    std::cout
        << name
        << "_J="
        << r.j
        << '\n';

    std::cout
        << name
        << "_BASE="
        << r.base
        << '\n';

    std::cout
        << name
        << "_LEFT_FACTOR="
        << r.left_factor
        << '\n';

    std::cout
        << name
        << "_RIGHT_FACTOR="
        << r.right_factor
        << '\n';

    std::cout
        << name
        << "_IDENTITY_LHS="
        << r.left_factor * r.right_factor
        << '\n';

    std::cout
        << name
        << "_IDENTITY_RHS="
        << r.j * prime + 1
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x46520260915ULL);

    const std::vector<int> primes =
        generate_primes(PRIME_LIMIT);

    const std::vector<PrimeCase> cases =
        generate_cases(primes, rng);

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << '\n';

    u64 p_found = 0;
    u64 q_found = 0;

    u64 identity_failures = 0;

    u64 j_equal = 0;
    u64 j_p_less_q = 0;
    u64 j_p_greater_q = 0;

    u64 m_equal = 0;
    u64 m_p_less_q = 0;
    u64 m_p_greater_q = 0;

    u64 k1_equal = 0;
    u64 k2_equal = 0;

    u64 mp_divides_q_minus_1 = 0;
    u64 mq_divides_p_minus_1 = 0;

    u64 mp_divides_q = 0;
    u64 mq_divides_p = 0;

    u64 jp_divides_q_minus_1 = 0;
    u64 jq_divides_p_minus_1 = 0;

    u64 cross_product_equal = 0;

    u64 min_j_p = UINT64_MAX;
    u64 min_j_q = UINT64_MAX;

    u64 max_j_p = 0;
    u64 max_j_q = 0;

    long double sum_j_p = 0.0L;
    long double sum_j_q = 0.0L;

    long double sum_m_p = 0.0L;
    long double sum_m_q = 0.0L;

    long double sum_j_product = 0.0L;
    long double sum_m_product = 0.0L;

    u64 first_j_equality = 0;
    u64 first_m_equality = 0;

    bool first_example_printed = false;

    for (std::size_t index = 0; index < cases.size(); ++index) {
        const PrimeCase& c = cases[index];

        /*
         * p and q are independently represented.
         */
        const Representation rp =
            find_min_j_representation(c.p, pairs);

        const Representation rq =
            find_min_j_representation(c.q, pairs);

        if (rp.found) {
            ++p_found;
        }

        if (rq.found) {
            ++q_found;
        }

        if (!(rp.found && rq.found)) {
            continue;
        }

        const u64 lhs_p =
            rp.left_factor * rp.right_factor;

        const u64 rhs_p =
            rp.j * c.p + 1;

        const u64 lhs_q =
            rq.left_factor * rq.right_factor;

        const u64 rhs_q =
            rq.j * c.q + 1;

        if (lhs_p != rhs_p || lhs_q != rhs_q) {
            ++identity_failures;
        }

        if (rp.j == rq.j) {
            ++j_equal;

            if (first_j_equality == 0) {
                first_j_equality = c.p;
            }
        } else if (rp.j < rq.j) {
            ++j_p_less_q;
        } else {
            ++j_p_greater_q;
        }

        if (rp.m == rq.m) {
            ++m_equal;

            if (first_m_equality == 0) {
                first_m_equality = c.p;
            }
        } else if (rp.m < rq.m) {
            ++m_p_less_q;
        } else {
            ++m_p_greater_q;
        }

        if (rp.k1 == rq.k1) {
            ++k1_equal;
        }

        if (rp.k2 == rq.k2) {
            ++k2_equal;
        }

        if (c.q > 1 &&
            (c.q - 1) % rp.m == 0) {
            ++mp_divides_q_minus_1;
        }

        if (c.p > 1 &&
            (c.p - 1) % rq.m == 0) {
            ++mq_divides_p_minus_1;
        }

        if (c.q % rp.m == 0) {
            ++mp_divides_q;
        }

        if (c.p % rq.m == 0) {
            ++mq_divides_p;
        }

        if (c.q > 1 &&
            (c.q - 1) % rp.j == 0) {
            ++jp_divides_q_minus_1;
        }

        if (c.p > 1 &&
            (c.p - 1) % rq.j == 0) {
            ++jq_divides_p_minus_1;
        }

        /*
         * Cross-product checks.
         *
         * These are deliberately simple:
         *
         *   j_p * M_p
         *   j_q * M_q
         */
        const u64 jpmp = rp.j * rp.m;
        const u64 jqmq = rq.j * rq.m;

        if (jpmp == jqmq) {
            ++cross_product_equal;
        }

        min_j_p = std::min(min_j_p, rp.j);
        min_j_q = std::min(min_j_q, rq.j);

        max_j_p = std::max(max_j_p, rp.j);
        max_j_q = std::max(max_j_q, rq.j);

        sum_j_p += static_cast<long double>(rp.j);
        sum_j_q += static_cast<long double>(rq.j);

        sum_m_p += static_cast<long double>(rp.m);
        sum_m_q += static_cast<long double>(rq.m);

        sum_j_product +=
            static_cast<long double>(rp.j) *
            static_cast<long double>(rq.j);

        sum_m_product +=
            static_cast<long double>(rp.m) *
            static_cast<long double>(rq.m);

        if (!first_example_printed) {
            print_representation("FIRST_P", c.p, rp);
            print_representation("FIRST_Q", c.q, rq);

            std::cout
                << "FIRST_N="
                << c.n
                << '\n';

            std::cout
                << "FIRST_JP_TIMES_MP="
                << jpmp
                << '\n';

            std::cout
                << "FIRST_JQ_TIMES_MQ="
                << jqmq
                << '\n';

            first_example_printed = true;
        }

        if ((index + 1) % 100 == 0) {
            std::cout
                << "PROGRESS="
                << (index + 1)
                << "/"
                << cases.size()
                << '\n';
        }
    }

    const u64 paired_cases =
        std::min(p_found, q_found);

    const double avg_j_p =
        paired_cases == 0
            ? 0.0
            : static_cast<double>(
                  sum_j_p /
                  static_cast<long double>(paired_cases)
              );

    const double avg_j_q =
        paired_cases == 0
            ? 0.0
            : static_cast<double>(
                  sum_j_q /
                  static_cast<long double>(paired_cases)
              );

    const double avg_m_p =
        paired_cases == 0
            ? 0.0
            : static_cast<double>(
                  sum_m_p /
                  static_cast<long double>(paired_cases)
              );

    const double avg_m_q =
        paired_cases == 0
            ? 0.0
            : static_cast<double>(
                  sum_m_q /
                  static_cast<long double>(paired_cases)
              );

    const double avg_j_product =
        paired_cases == 0
            ? 0.0
            : static_cast<double>(
                  sum_j_product /
                  static_cast<long double>(paired_cases)
              );

    const double avg_m_product =
        paired_cases == 0
            ? 0.0
            : static_cast<double>(
                  sum_m_product /
                  static_cast<long double>(paired_cases)
              );

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "P_FOUND="
        << p_found
        << '\n';

    std::cout
        << "Q_FOUND="
        << q_found
        << '\n';

    std::cout
        << "PAIRED_CASES="
        << paired_cases
        << '\n';

    std::cout
        << "IDENTITY_FAILURES="
        << identity_failures
        << '\n';

    std::cout
        << "J_EQUAL="
        << j_equal
        << '\n';

    std::cout
        << "J_P_LESS_Q="
        << j_p_less_q
        << '\n';

    std::cout
        << "J_P_GREATER_Q="
        << j_p_greater_q
        << '\n';

    std::cout
        << "M_EQUAL="
        << m_equal
        << '\n';

    std::cout
        << "M_P_LESS_Q="
        << m_p_less_q
        << '\n';

    std::cout
        << "M_P_GREATER_Q="
        << m_p_greater_q
        << '\n';

    std::cout
        << "K1_EQUAL="
        << k1_equal
        << '\n';

    std::cout
        << "K2_EQUAL="
        << k2_equal
        << '\n';

    std::cout
        << "MP_DIVIDES_Q_MINUS_1="
        << mp_divides_q_minus_1
        << '\n';

    std::cout
        << "MQ_DIVIDES_P_MINUS_1="
        << mq_divides_p_minus_1
        << '\n';

    std::cout
        << "MP_DIVIDES_Q="
        << mp_divides_q
        << '\n';

    std::cout
        << "MQ_DIVIDES_P="
        << mq_divides_p
        << '\n';

    std::cout
        << "JP_DIVIDES_Q_MINUS_1="
        << jp_divides_q_minus_1
        << '\n';

    std::cout
        << "JQ_DIVIDES_P_MINUS_1="
        << jq_divides_p_minus_1
        << '\n';

    std::cout
        << "JPMP_EQUALS_JQMQ="
        << cross_product_equal
        << '\n';

    std::cout
        << "MIN_J_P="
        << (min_j_p == UINT64_MAX ? 0 : min_j_p)
        << '\n';

    std::cout
        << "MIN_J_Q="
        << (min_j_q == UINT64_MAX ? 0 : min_j_q)
        << '\n';

    std::cout
        << "MAX_J_P="
        << max_j_p
        << '\n';

    std::cout
        << "MAX_J_Q="
        << max_j_q
        << '\n';

    std::cout
        << "AVG_J_P="
        << avg_j_p
        << '\n';

    std::cout
        << "AVG_J_Q="
        << avg_j_q
        << '\n';

    std::cout
        << "AVG_M_P="
        << avg_m_p
        << '\n';

    std::cout
        << "AVG_M_Q="
        << avg_m_q
        << '\n';

    std::cout
        << "AVG_J_PRODUCT="
        << avg_j_product
        << '\n';

    std::cout
        << "AVG_M_PRODUCT="
        << avg_m_product
        << '\n';

    std::cout
        << "FIRST_J_EQUAL_P="
        << first_j_equality
        << '\n';

    std::cout
        << "FIRST_M_EQUAL_P="
        << first_m_equality
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';
}

int main() {
    main_experiment();
    return 0;
}
