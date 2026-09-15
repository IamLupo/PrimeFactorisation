#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 466;

constexpr int PRIME_LIMIT = 100000;
constexpr int PRIME_MIN = 10000;
constexpr int CASE_COUNT = 500;

constexpr int K_LIMIT = 1000;
constexpr int M_LIMIT = 7;

struct PrimeCase {
    u64 p;
    u64 q;
    u64 r;
    u64 n;
};

struct CRTPair {
    int k1;
    int k2;
};

struct Representation {
    bool found = false;

    u64 prime = 0;

    u64 k1 = 0;
    u64 k2 = 0;

    u64 a = 0;
    u64 m = 0;
    u64 j = 0;

    u64 left_factor = 0;
    u64 right_factor = 0;
};

struct CaseRepresentations {
    Representation p;
    Representation q;
    Representation r;
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

    for (int prime : primes) {
        if (prime >= PRIME_MIN && prime <= PRIME_LIMIT) {
            candidates.push_back(prime);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        candidates.size() - 1
    );

    std::vector<PrimeCase> cases;
    cases.reserve(CASE_COUNT);

    for (int i = 0; i < CASE_COUNT; ++i) {
        u64 p = static_cast<u64>(candidates[dist(rng)]);
        u64 q = static_cast<u64>(candidates[dist(rng)]);
        u64 r = static_cast<u64>(candidates[dist(rng)]);

        while (q == p) {
            q = static_cast<u64>(candidates[dist(rng)]);
        }

        while (r == p || r == q) {
            r = static_cast<u64>(candidates[dist(rng)]);
        }

        std::vector<u64> values;
        values.reserve(3);
        values.push_back(p);
        values.push_back(q);
        values.push_back(r);

        std::sort(values.begin(), values.end());

        PrimeCase c;
        c.p = values[0];
        c.q = values[1];
        c.r = values[2];
        c.n = c.p * c.q * c.r;

        cases.push_back(c);
    }

    return cases;
}

std::vector<CRTPair> build_determinant_one_pairs() {
    std::vector<CRTPair> pairs;

    for (int k1 = 1; k1 <= K_LIMIT; ++k1) {
        for (int k2 = k1; k2 <= K_LIMIT; ++k2) {
            if (std::gcd(k1, k2) != 1) {
                continue;
            }

            bool valid = false;

            /*
             * Need:
             *
             *     m2*k1 - m1*k2 = 1
             *
             * for some
             *
             *     1 <= m1,m2 <= M_LIMIT.
             */
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
    best.prime = prime;

    for (const CRTPair& pair : pairs) {
        const u64 k1 = static_cast<u64>(pair.k1);
        const u64 k2 = static_cast<u64>(pair.k2);

        const u64 a = k1 + k2;
        const u64 m = k1 * k2;

        if (a > prime) {
            continue;
        }

        const u64 remainder = prime - a;

        if (remainder % m != 0) {
            continue;
        }

        const u64 j = remainder / m;

        bool take = false;

        if (!best.found) {
            take = true;
        } else if (j < best.j) {
            take = true;
        } else if (
            j == best.j &&
            m < best.m
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) <
                std::max(best.k1, best.k2)
        ) {
            take = true;
        } else if (
            j == best.j &&
            m == best.m &&
            std::max(k1, k2) ==
                std::max(best.k1, best.k2) &&
            k1 < best.k1
        ) {
            take = true;
        }

        if (!take) {
            continue;
        }

        best.found = true;
        best.k1 = k1;
        best.k2 = k2;
        best.a = a;
        best.m = m;
        best.j = j;

        best.left_factor = j * k1 + 1;
        best.right_factor = j * k2 + 1;
    }

    return best;
}

CaseRepresentations build_case_representations(
    const PrimeCase& c,
    const std::vector<CRTPair>& pairs
) {
    CaseRepresentations result;

    result.p =
        find_min_j_representation(c.p, pairs);

    result.q =
        find_min_j_representation(c.q, pairs);

    result.r =
        find_min_j_representation(c.r, pairs);

    return result;
}

bool divisible_difference(
    u64 divisor,
    u64 value,
    u64 base
) {
    if (divisor == 0) {
        return false;
    }

    /*
     * We only call this with value >= base in this experiment.
     */
    if (value < base) {
        return false;
    }

    return (value - base) % divisor == 0;
}

void print_representation(
    const char* label,
    const Representation& rep
) {
    std::cout
        << label
        << "_PRIME="
        << rep.prime
        << '\n';

    std::cout
        << label
        << "_K1="
        << rep.k1
        << '\n';

    std::cout
        << label
        << "_K2="
        << rep.k2
        << '\n';

    std::cout
        << label
        << "_A="
        << rep.a
        << '\n';

    std::cout
        << label
        << "_M="
        << rep.m
        << '\n';

    std::cout
        << label
        << "_J="
        << rep.j
        << '\n';

    std::cout
        << label
        << "_CHECK="
        << rep.a + rep.j * rep.m
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(0x46620260915ULL);

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

    u64 all_found = 0;
    u64 identity_failures = 0;

    /*
     * Cross residual tests:
     *
     * Mp | q - aq
     * Mp | r - ar
     *
     * and cyclic variants.
     */
    u64 mp_q = 0;
    u64 mp_r = 0;

    u64 mq_p = 0;
    u64 mq_r = 0;

    u64 mr_p = 0;
    u64 mr_q = 0;

    u64 mp_both = 0;
    u64 mq_both = 0;
    u64 mr_both = 0;

    u64 every_other_residual = 0;

    /*
     * M gcd structure.
     */
    u64 gcd_mp_mq_gt1 = 0;
    u64 gcd_mp_mr_gt1 = 0;
    u64 gcd_mq_mr_gt1 = 0;

    u64 all_three_pairwise_coprime = 0;

    /*
     * Pairwise proposed congruences:
     *
     * N == ap*aq (mod lcm(Mp,Mq))
     *
     * etc.
     */
    u64 lcm_pq = 0;
    u64 lcm_pr = 0;
    u64 lcm_qr = 0;

    /*
     * Triple proposed congruence:
     *
     * N == ap*aq*ar
     *     (mod lcm(Mp,Mq,Mr))
     */
    u64 lcm_triple = 0;

    /*
     * Cross-j tests.
     */
    u64 jp_q_minus_aq = 0;
    u64 jp_r_minus_ar = 0;

    u64 jq_p_minus_ap = 0;
    u64 jq_r_minus_ar = 0;

    u64 jr_p_minus_ap = 0;
    u64 jr_q_minus_aq = 0;

    /*
     * Direct product residual tests.
     */
    u64 mp_divides_qr_minus_aqar = 0;
    u64 mq_divides_pr_minus_apar = 0;
    u64 mr_divides_pq_minus_apaq = 0;

    /*
     * Residual statistics.
     */
    u64 max_p_residual = 0;
    u64 max_q_residual = 0;
    u64 max_r_residual = 0;

    long double sum_p_residual = 0.0L;
    long double sum_q_residual = 0.0L;
    long double sum_r_residual = 0.0L;

    bool first_example = false;

    for (int index = 0; index < CASE_COUNT; ++index) {
        const PrimeCase& c = cases[index];

        const CaseRepresentations reps =
            build_case_representations(c, pairs);

        if (!(
            reps.p.found &&
            reps.q.found &&
            reps.r.found
        )) {
            continue;
        }

        ++all_found;

        const Representation& p = reps.p;
        const Representation& q = reps.q;
        const Representation& r = reps.r;

        /*
         * Verify:
         *
         * p = ap + jp*Mp
         * q = aq + jq*Mq
         * r = ar + jr*Mr
         */
        if (
            p.a + p.j * p.m != c.p ||
            q.a + q.j * q.m != c.q ||
            r.a + r.j * r.m != c.r
        ) {
            ++identity_failures;
        }

        const bool mpq =
            divisible_difference(
                p.m,
                c.q,
                q.a
            );

        const bool mpr =
            divisible_difference(
                p.m,
                c.r,
                r.a
            );

        const bool mqp =
            divisible_difference(
                q.m,
                c.p,
                p.a
            );

        const bool mqr =
            divisible_difference(
                q.m,
                c.r,
                r.a
            );

        const bool mrp =
            divisible_difference(
                r.m,
                c.p,
                p.a
            );

        const bool mrq =
            divisible_difference(
                r.m,
                c.q,
                q.a
            );

        if (mpq) ++mp_q;
        if (mpr) ++mp_r;
        if (mqp) ++mq_p;
        if (mqr) ++mq_r;
        if (mrp) ++mr_p;
        if (mrq) ++mr_q;

        if (mpq && mpr) {
            ++mp_both;
        }

        if (mqp && mqr) {
            ++mq_both;
        }

        if (mrp && mrq) {
            ++mr_both;
        }

        if (
            mpq && mpr &&
            mqp && mqr &&
            mrp && mrq
        ) {
            ++every_other_residual;
        }

        /*
         * gcd(Mi, Mj)
         */
        const u64 gcd_pq =
            std::gcd(p.m, q.m);

        const u64 gcd_pr =
            std::gcd(p.m, r.m);

        const u64 gcd_qr =
            std::gcd(q.m, r.m);

        if (gcd_pq > 1) {
            ++gcd_mp_mq_gt1;
        }

        if (gcd_pr > 1) {
            ++gcd_mp_mr_gt1;
        }

        if (gcd_qr > 1) {
            ++gcd_mq_mr_gt1;
        }

        if (
            gcd_pq == 1 &&
            gcd_pr == 1 &&
            gcd_qr == 1
        ) {
            ++all_three_pairwise_coprime;
        }

        /*
         * Pairwise lcm congruences.
         */
        const u64 lcm_pq_value =
            std::lcm(p.m, q.m);

        const u64 lcm_pr_value =
            std::lcm(p.m, r.m);

        const u64 lcm_qr_value =
            std::lcm(q.m, r.m);

        if (
            lcm_pq_value != 0 &&
            c.n >= p.a * q.a &&
            (c.n - p.a * q.a) %
                lcm_pq_value == 0
        ) {
            ++lcm_pq;
        }

        if (
            lcm_pr_value != 0 &&
            c.n >= p.a * r.a &&
            (c.n - p.a * r.a) %
                lcm_pr_value == 0
        ) {
            ++lcm_pr;
        }

        if (
            lcm_qr_value != 0 &&
            c.n >= q.a * r.a &&
            (c.n - q.a * r.a) %
                lcm_qr_value == 0
        ) {
            ++lcm_qr;
        }

        /*
         * Triple lcm.
         */
        const u64 gcd_pq_for_lcm =
            std::gcd(p.m, q.m);

        const u64 lcm_pq_for_triple =
            (p.m / gcd_pq_for_lcm) * q.m;

        const u64 gcd_all =
            std::gcd(lcm_pq_for_triple, r.m);

        const u64 lcm_pqr =
            (lcm_pq_for_triple / gcd_all) * r.m;

        const u64 base_product =
            p.a * q.a * r.a;

        if (
            lcm_pqr != 0 &&
            c.n >= base_product &&
            (c.n - base_product) %
                lcm_pqr == 0
        ) {
            ++lcm_triple;
        }

        /*
         * Cross-j tests.
         */
        if (
            divisible_difference(
                p.j,
                c.q,
                q.a
            )
        ) {
            ++jp_q_minus_aq;
        }

        if (
            divisible_difference(
                p.j,
                c.r,
                r.a
            )
        ) {
            ++jp_r_minus_ar;
        }

        if (
            divisible_difference(
                q.j,
                c.p,
                p.a
            )
        ) {
            ++jq_p_minus_ap;
        }

        if (
            divisible_difference(
                q.j,
                c.r,
                r.a
            )
        ) {
            ++jq_r_minus_ar;
        }

        if (
            divisible_difference(
                r.j,
                c.p,
                p.a
            )
        ) {
            ++jr_p_minus_ap;
        }

        if (
            divisible_difference(
                r.j,
                c.q,
                q.a
            )
        ) {
            ++jr_q_minus_aq;
        }

        /*
         * Product residual tests:
         *
         * Mp | qr - aq*ar
         * Mq | pr - ap*ar
         * Mr | pq - ap*aq
         */
        const u64 qr_base = q.a * r.a;
        const u64 pr_base = p.a * r.a;
        const u64 pq_base = p.a * q.a;

        if (
            c.q * c.r >= qr_base &&
            (c.q * c.r - qr_base) % p.m == 0
        ) {
            ++mp_divides_qr_minus_aqar;
        }

        if (
            c.p * c.r >= pr_base &&
            (c.p * c.r - pr_base) % q.m == 0
        ) {
            ++mq_divides_pr_minus_apar;
        }

        if (
            c.p * c.q >= pq_base &&
            (c.p * c.q - pq_base) % r.m == 0
        ) {
            ++mr_divides_pq_minus_apaq;
        }

        /*
         * Residuals:
         *
         * p-ap = jp*Mp
         * q-aq = jq*Mq
         * r-ar = jr*Mr
         */
        const u64 p_residual = c.p - p.a;
        const u64 q_residual = c.q - q.a;
        const u64 r_residual = c.r - r.a;

        max_p_residual =
            std::max(max_p_residual, p_residual);

        max_q_residual =
            std::max(max_q_residual, q_residual);

        max_r_residual =
            std::max(max_r_residual, r_residual);

        sum_p_residual +=
            static_cast<long double>(p_residual);

        sum_q_residual +=
            static_cast<long double>(q_residual);

        sum_r_residual +=
            static_cast<long double>(r_residual);

        if (!first_example) {
            std::cout
                << "FIRST_N="
                << c.n
                << '\n';

            print_representation("FIRST_P", p);
            print_representation("FIRST_Q", q);
            print_representation("FIRST_R", r);

            std::cout
                << "FIRST_MP_Q="
                << (mpq ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_MP_R="
                << (mpr ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_MQ_P="
                << (mqp ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_MQ_R="
                << (mqr ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_MR_P="
                << (mrp ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_MR_Q="
                << (mrq ? 1 : 0)
                << '\n';

            first_example = true;
        }

        if ((index + 1) % 100 == 0) {
            std::cout
                << "PROGRESS="
                << (index + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const double avg_p_residual =
        all_found == 0
            ? 0.0
            : static_cast<double>(
                sum_p_residual /
                static_cast<long double>(all_found)
            );

    const double avg_q_residual =
        all_found == 0
            ? 0.0
            : static_cast<double>(
                sum_q_residual /
                static_cast<long double>(all_found)
            );

    const double avg_r_residual =
        all_found == 0
            ? 0.0
            : static_cast<double>(
                sum_r_residual /
                static_cast<long double>(all_found)
            );

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "ALL_FOUND="
        << all_found
        << '\n';

    std::cout
        << "IDENTITY_FAILURES="
        << identity_failures
        << '\n';

    std::cout
        << "MP_DIVIDES_Q_MINUS_AQ="
        << mp_q
        << '\n';

    std::cout
        << "MP_DIVIDES_R_MINUS_AR="
        << mp_r
        << '\n';

    std::cout
        << "MQ_DIVIDES_P_MINUS_AP="
        << mq_p
        << '\n';

    std::cout
        << "MQ_DIVIDES_R_MINUS_AR="
        << mq_r
        << '\n';

    std::cout
        << "MR_DIVIDES_P_MINUS_AP="
        << mr_p
        << '\n';

    std::cout
        << "MR_DIVIDES_Q_MINUS_AQ="
        << mr_q
        << '\n';

    std::cout
        << "MP_DIVIDES_BOTH="
        << mp_both
        << '\n';

    std::cout
        << "MQ_DIVIDES_BOTH="
        << mq_both
        << '\n';

    std::cout
        << "MR_DIVIDES_BOTH="
        << mr_both
        << '\n';

    std::cout
        << "EVERY_M_DIVIDES_BOTH_OTHER_RESIDUALS="
        << every_other_residual
        << '\n';

    std::cout
        << "GCD_MP_MQ_GT1="
        << gcd_mp_mq_gt1
        << '\n';

    std::cout
        << "GCD_MP_MR_GT1="
        << gcd_mp_mr_gt1
        << '\n';

    std::cout
        << "GCD_MQ_MR_GT1="
        << gcd_mq_mr_gt1
        << '\n';

    std::cout
        << "ALL_THREE_M_PAIRWISE_COPRIME="
        << all_three_pairwise_coprime
        << '\n';

    std::cout
        << "LCM_PQ_CONGRUENCE="
        << lcm_pq
        << '\n';

    std::cout
        << "LCM_PR_CONGRUENCE="
        << lcm_pr
        << '\n';

    std::cout
        << "LCM_QR_CONGRUENCE="
        << lcm_qr
        << '\n';

    std::cout
        << "LCM_TRIPLE_CONGRUENCE="
        << lcm_triple
        << '\n';

    std::cout
        << "JP_DIVIDES_Q_MINUS_AQ="
        << jp_q_minus_aq
        << '\n';

    std::cout
        << "JP_DIVIDES_R_MINUS_AR="
        << jp_r_minus_ar
        << '\n';

    std::cout
        << "JQ_DIVIDES_P_MINUS_AP="
        << jq_p_minus_ap
        << '\n';

    std::cout
        << "JQ_DIVIDES_R_MINUS_AR="
        << jq_r_minus_ar
        << '\n';

    std::cout
        << "JR_DIVIDES_P_MINUS_AP="
        << jr_p_minus_ap
        << '\n';

    std::cout
        << "JR_DIVIDES_Q_MINUS_AQ="
        << jr_q_minus_aq
        << '\n';

    std::cout
        << "MP_DIVIDES_QR_MINUS_AQAR="
        << mp_divides_qr_minus_aqar
        << '\n';

    std::cout
        << "MQ_DIVIDES_PR_MINUS_APAR="
        << mq_divides_pr_minus_apar
        << '\n';

    std::cout
        << "MR_DIVIDES_PQ_MINUS_APAQ="
        << mr_divides_pq_minus_apaq
        << '\n';

    std::cout
        << "AVG_P_RESIDUAL="
        << avg_p_residual
        << '\n';

    std::cout
        << "AVG_Q_RESIDUAL="
        << avg_q_residual
        << '\n';

    std::cout
        << "AVG_R_RESIDUAL="
        << avg_r_residual
        << '\n';

    std::cout
        << "MAX_P_RESIDUAL="
        << max_p_residual
        << '\n';

    std::cout
        << "MAX_Q_RESIDUAL="
        << max_q_residual
        << '\n';

    std::cout
        << "MAX_R_RESIDUAL="
        << max_r_residual
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