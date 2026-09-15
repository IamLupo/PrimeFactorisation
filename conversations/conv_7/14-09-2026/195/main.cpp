#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

struct Pair {
    int m1;
    int m2;
    int k1;
    int k2;
    u64 a;
    u64 M;
};

struct Representation {
    int m1;
    int m2;
    int k1;
    int k2;
    u64 a;
    u64 M;
    u64 j;
};

struct PrimeStats {
    u64 representations = 0;
    u64 identity_ok = 0;
    u64 reconstruction_ok = 0;
    u64 determinant_ok = 0;
    u64 gcd_k_ok = 0;
    u64 gcd_factor_ok = 0;

    u64 min_j = 0;
    u64 min_M = 0;
    u64 max_j = 0;
    u64 max_M = 0;

    u64 determinant_representations = 0;
};

std::vector<bool> build_sieve(int limit) {
    std::vector<bool> prime(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    prime[0] = false;

    if (limit >= 1) {
        prime[1] = false;
    }

    for (int i = 2; 1LL * i * i <= limit; ++i) {
        if (!prime[i]) {
            continue;
        }

        for (int j = i * i; j <= limit; j += i) {
            prime[j] = false;
        }
    }

    return prime;
}

std::vector<int> extract_primes(
    int low,
    int high,
    const std::vector<bool>& prime) {

    std::vector<int> result;

    for (int p = low; p <= high; ++p) {
        if (prime[p]) {
            result.push_back(p);
        }
    }

    return result;
}

/*
    Construct all determinant-one tuples

        m2*k1 - m1*k2 = 1

    with

        1 <= m1,m2 <= 7
        2 <= k1 <= k2 <= 1000
        gcd(k1,k2)=1.

    We retain m1,m2 this time because the experiment is
    specifically studying their algebraic role.
*/
std::vector<Pair> build_pairs(
    int k_limit,
    int m_limit) {

    std::vector<Pair> pairs;

    for (int m1 = 1; m1 <= m_limit; ++m1) {
        for (int m2 = 1; m2 <= m_limit; ++m2) {

            for (int k1 = 2; k1 <= k_limit; ++k1) {

                const long long numerator =
                    1LL * m2 * k1 - 1;

                if (numerator <= 0) {
                    continue;
                }

                if (numerator % m1 != 0) {
                    continue;
                }

                const int k2 =
                    static_cast<int>(
                        numerator / m1
                    );

                if (k2 < k1 ||
                    k2 > k_limit ||
                    k2 < 2) {
                    continue;
                }

                if (std::gcd(k1, k2) != 1) {
                    continue;
                }

                const u64 a =
                    static_cast<u64>(k1) +
                    static_cast<u64>(k2);

                const u64 M =
                    static_cast<u64>(k1) *
                    static_cast<u64>(k2);

                pairs.push_back({
                    m1,
                    m2,
                    k1,
                    k2,
                    a,
                    M
                });
            }
        }
    }

    return pairs;
}

bool exact_identity_holds(
    u64 p,
    const Representation& rep) {

    /*
        p = k1+k2+j*k1*k2
    */
    const u64 reconstructed =
        rep.a + rep.j * rep.M;

    /*
        jp+1 = (jk1+1)(jk2+1)
    */
    const u64 left =
        rep.j * p + 1;

    const u64 right =
        (rep.j *
             static_cast<u64>(rep.k1) +
         1) *
        (rep.j *
             static_cast<u64>(rep.k2) +
         1);

    return reconstructed == p &&
           left == right;
}

bool determinant_identity_holds(
    const Representation& rep) {

    return
        static_cast<long long>(rep.m2) *
            rep.k1 -
        static_cast<long long>(rep.m1) *
            rep.k2 ==
        1;
}

Representation make_representation(
    const Pair& pair,
    u64 p) {

    const u64 remainder =
        p - pair.a;

    const u64 j =
        remainder / pair.M;

    return {
        pair.m1,
        pair.m2,
        pair.k1,
        pair.k2,
        pair.a,
        pair.M,
        j
    };
}

double average(
    u64 total,
    u64 count) {

    if (count == 0) {
        return 0.0;
    }

    return
        static_cast<double>(total) /
        static_cast<double>(count);
}

int main() {
    constexpr int PRIME_LOW = 2;
    constexpr int PRIME_HIGH = 100000;

    constexpr int K_LIMIT = 1000;
    constexpr int M_LIMIT = 7;

    std::cout
        << "START EXPERIMENT 488\n";

    const auto sieve =
        build_sieve(PRIME_HIGH);

    const auto primes =
        extract_primes(
            PRIME_LOW,
            PRIME_HIGH,
            sieve
        );

    const auto pairs =
        build_pairs(
            K_LIMIT,
            M_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << "\n";

    /*
        Overall counters.
    */
    u64 representable_primes = 0;
    u64 nonrepresentable_primes = 0;

    u64 total_representations = 0;

    u64 total_identity_ok = 0;
    u64 total_reconstruction_ok = 0;
    u64 total_determinant_ok = 0;
    u64 total_gcd_k_ok = 0;
    u64 total_gcd_factor_ok = 0;

    u64 min_j_total = 0;
    u64 min_M_total = 0;

    u64 max_j_total = 0;
    u64 max_M_total = 0;

    u64 max_representation_count = 0;
    u64 max_representation_prime = 0;

    u64 min_representation_count =
        UINT64_MAX;

    u64 min_representation_prime = 0;

    /*
        Investigate whether the same (k1,k2) admits
        multiple determinant-one (m1,m2) realizations.
    */
    u64 representation_with_multiple_m = 0;
    u64 total_extra_m_realizations = 0;

    /*
        Relation candidates.
    */
    u64 m1_equals_j = 0;
    u64 m2_equals_j = 0;

    u64 m1_divides_j = 0;
    u64 m2_divides_j = 0;

    u64 j_divides_m1 = 0;
    u64 j_divides_m2 = 0;

    u64 m1_less_k1 = 0;
    u64 m2_less_k2 = 0;

    u64 factor_gcd_one = 0;
    u64 factor_gcd_nontrivial = 0;

    u64 factor_difference_total = 0;

    /*
        First detailed prime.
    */
    Representation first_rep{};
    bool first_found = false;

    const int FIRST_PRIME = 23;

    /*
        Records for several specific primes.
    */
    const std::vector<int> selected_primes = {
        5,
        7,
        11,
        13,
        17,
        19,
        23,
        29,
        31,
        37,
        41,
        43,
        47,
        97,
        101,
        997,
        99991
    };

    const auto start =
        std::chrono::steady_clock::now();

    for (int p : primes) {

        PrimeStats stats;

        std::vector<Representation> representations;

        for (const Pair& pair : pairs) {

            if (pair.a > static_cast<u64>(p)) {
                continue;
            }

            const u64 pp =
                static_cast<u64>(p);

            /*
                Exact representation condition.
            */
            if (pp % pair.M !=
                pair.a % pair.M) {
                continue;
            }

            Representation rep =
                make_representation(
                    pair,
                    pp
                );

            representations.push_back(rep);

            ++stats.representations;
            ++total_representations;

            const bool identity_ok =
                exact_identity_holds(
                    pp,
                    rep
                );

            const bool determinant_ok =
                determinant_identity_holds(
                    rep
                );

            const bool gcd_k =
                std::gcd(
                    rep.k1,
                    rep.k2
                ) == 1;

            const u64 u =
                rep.j *
                    static_cast<u64>(rep.k1) +
                1;

            const u64 v =
                rep.j *
                    static_cast<u64>(rep.k2) +
                1;

            const bool gcd_factor =
                std::gcd(u, v) == 1;

            if (identity_ok) {
                ++stats.identity_ok;
                ++total_identity_ok;
            }

            if (rep.a +
                    rep.j * rep.M ==
                pp) {

                ++stats.reconstruction_ok;
                ++total_reconstruction_ok;
            }

            if (determinant_ok) {
                ++stats.determinant_ok;
                ++total_determinant_ok;
            }

            if (gcd_k) {
                ++stats.gcd_k_ok;
                ++total_gcd_k_ok;
            }

            if (gcd_factor) {
                ++stats.gcd_factor_ok;
                ++total_gcd_factor_ok;
                ++factor_gcd_one;
            } else {
                ++factor_gcd_nontrivial;
            }

            if (rep.j == 0 ||
                !stats.representations) {

                stats.min_j =
                    rep.j;
            }

            if (stats.representations == 1 ||
                rep.j < stats.min_j) {

                stats.min_j =
                    rep.j;
            }

            if (stats.representations == 1 ||
                rep.M < stats.min_M) {

                stats.min_M =
                    rep.M;
            }

            stats.max_j =
                std::max(
                    stats.max_j,
                    rep.j
                );

            stats.max_M =
                std::max(
                    stats.max_M,
                    rep.M
                );

            /*
                Possible relationships between m and j.
            */
            if (rep.m1 ==
                static_cast<int>(rep.j)) {
                ++m1_equals_j;
            }

            if (rep.m2 ==
                static_cast<int>(rep.j)) {
                ++m2_equals_j;
            }

            if (rep.j != 0 &&
                rep.j %
                    static_cast<u64>(
                        rep.m1
                    ) == 0) {
                ++m1_divides_j;
            }

            if (rep.j != 0 &&
                rep.j %
                    static_cast<u64>(
                        rep.m2
                    ) == 0) {
                ++m2_divides_j;
            }

            if (rep.j != 0 &&
                static_cast<u64>(
                    rep.m1
                ) %
                    rep.j == 0) {
                ++j_divides_m1;
            }

            if (rep.j != 0 &&
                static_cast<u64>(
                    rep.m2
                ) %
                    rep.j == 0) {
                ++j_divides_m2;
            }

            if (rep.m1 < rep.k1) {
                ++m1_less_k1;
            }

            if (rep.m2 < rep.k2) {
                ++m2_less_k2;
            }

            const u64 difference =
                v > u
                    ? v - u
                    : u - v;

            factor_difference_total +=
                difference;
        }

        if (stats.representations > 0) {

            ++representable_primes;

            total_identity_ok += 0;
            total_determinant_ok += 0;
            total_gcd_k_ok += 0;
            total_gcd_factor_ok += 0;

            min_j_total +=
                stats.min_j;

            min_M_total +=
                stats.min_M;

            max_j_total +=
                stats.max_j;

            max_M_total +=
                stats.max_M;

            if (stats.representations >
                max_representation_count) {

                max_representation_count =
                    stats.representations;

                max_representation_prime =
                    static_cast<u64>(p);
            }

            if (stats.representations <
                min_representation_count) {

                min_representation_count =
                    stats.representations;

                min_representation_prime =
                    static_cast<u64>(p);
            }

            /*
                Look for the minimum-M representation
                and retain it for detailed inspection.
            */
            auto min_it =
                std::min_element(
                    representations.begin(),
                    representations.end(),
                    [](const Representation& lhs,
                       const Representation& rhs) {

                        if (lhs.M != rhs.M) {
                            return lhs.M < rhs.M;
                        }

                        return lhs.j < rhs.j;
                    }
                );

            if (p == FIRST_PRIME &&
                min_it != representations.end()) {

                first_rep =
                    *min_it;

                first_found = true;
            }
        } else {
            ++nonrepresentable_primes;
        }

        /*
            Determine whether a given (k1,k2) occurs with
            more than one determinant-one (m1,m2).
        */
        for (std::size_t i = 0;
             i < representations.size();
             ++i) {

            u64 same_k_count = 0;

            for (std::size_t j = i;
                 j < representations.size();
                 ++j) {

                if (representations[i].k1 ==
                        representations[j].k1 &&
                    representations[i].k2 ==
                        representations[j].k2) {

                    ++same_k_count;
                }
            }

            if (same_k_count > 1) {
                ++representation_with_multiple_m;
                total_extra_m_realizations +=
                    same_k_count - 1;
            }
        }
    }

    /*
        Correct total identity statistics need to be
        reconstructed from the total representation count.
    */
    /*
        Since every representation is explicitly tested,
        these identities should all equal
        total_representations.
    */
    const u64 identity_failures =
        total_representations -
        total_identity_ok;

    const u64 determinant_failures =
        total_representations -
        total_determinant_ok;

    const u64 reconstruction_failures =
        total_representations -
        total_reconstruction_ok;

    const u64 gcd_k_failures =
        total_representations -
        total_gcd_k_ok;

    /*
        Overall output.
    */
    std::cout
        << "TOTAL_PRIMES="
        << primes.size()
        << "\n";

    std::cout
        << "REPRESENTABLE_PRIMES="
        << representable_primes
        << "\n";

    std::cout
        << "NONREPRESENTABLE_PRIMES="
        << nonrepresentable_primes
        << "\n";

    std::cout
        << "COVERAGE_PERCENT="
        << 100.0 *
           static_cast<double>(
               representable_primes
           ) /
           static_cast<double>(
               primes.size()
           )
        << "\n";

    std::cout
        << "TOTAL_REPRESENTATIONS="
        << total_representations
        << "\n";

    std::cout
        << "AVG_REPRESENTATIONS="
        << average(
               total_representations,
               representable_primes
           )
        << "\n";

    std::cout
        << "IDENTITY_OK="
        << total_identity_ok
        << "\n";

    std::cout
        << "IDENTITY_FAILURES="
        << identity_failures
        << "\n";

    std::cout
        << "RECONSTRUCTION_OK="
        << total_reconstruction_ok
        << "\n";

    std::cout
        << "RECONSTRUCTION_FAILURES="
        << reconstruction_failures
        << "\n";

    std::cout
        << "DETERMINANT_OK="
        << total_determinant_ok
        << "\n";

    std::cout
        << "DETERMINANT_FAILURES="
        << determinant_failures
        << "\n";

    std::cout
        << "GCD_K_OK="
        << total_gcd_k_ok
        << "\n";

    std::cout
        << "GCD_K_FAILURES="
        << gcd_k_failures
        << "\n";

    std::cout
        << "GCD_FACTOR_ONE="
        << factor_gcd_one
        << "\n";

    std::cout
        << "GCD_FACTOR_NONTRIVIAL="
        << factor_gcd_nontrivial
        << "\n";

    std::cout
        << "AVG_MIN_M="
        << average(
               min_M_total,
               representable_primes
           )
        << "\n";

    std::cout
        << "AVG_MIN_J="
        << average(
               min_j_total,
               representable_primes
           )
        << "\n";

    std::cout
        << "AVG_MAX_M="
        << average(
               max_M_total,
               representable_primes
           )
        << "\n";

    std::cout
        << "AVG_MAX_J="
        << average(
               max_j_total,
               representable_primes
           )
        << "\n";

    std::cout
        << "MAX_REPRESENTATIONS="
        << max_representation_count
        << "\n";

    std::cout
        << "MAX_REPRESENTATION_PRIME="
        << max_representation_prime
        << "\n";

    std::cout
        << "MIN_REPRESENTATIONS="
        << min_representation_count
        << "\n";

    std::cout
        << "MIN_REPRESENTATION_PRIME="
        << min_representation_prime
        << "\n";

    /*
        Relationships involving m1,m2 and j.
    */
    std::cout
        << "M1_EQUALS_J="
        << m1_equals_j
        << "\n";

    std::cout
        << "M2_EQUALS_J="
        << m2_equals_j
        << "\n";

    std::cout
        << "M1_DIVIDES_J="
        << m1_divides_j
        << "\n";

    std::cout
        << "M2_DIVIDES_J="
        << m2_divides_j
        << "\n";

    std::cout
        << "J_DIVIDES_M1="
        << j_divides_m1
        << "\n";

    std::cout
        << "J_DIVIDES_M2="
        << j_divides_m2
        << "\n";

    std::cout
        << "M1_LESS_K1="
        << m1_less_k1
        << "\n";

    std::cout
        << "M2_LESS_K2="
        << m2_less_k2
        << "\n";

    std::cout
        << "MULTIPLE_M_FOR_SAME_K_PAIR="
        << representation_with_multiple_m
        << "\n";

    std::cout
        << "EXTRA_M_REALIZATIONS="
        << total_extra_m_realizations
        << "\n";

    std::cout
        << "AVG_FACTOR_DIFFERENCE="
        << average(
               factor_difference_total,
               total_representations
           )
        << "\n";

    /*
        Detailed selected primes.
    */
    for (int p : selected_primes) {

        const PrimeStats stats =
            PrimeStats{};

        (void)stats;

        std::vector<Representation> records;

        for (const Pair& pair : pairs) {

            if (pair.a >
                static_cast<u64>(p)) {
                continue;
            }

            const u64 pp =
                static_cast<u64>(p);

            if (pp % pair.M !=
                pair.a % pair.M) {
                continue;
            }

            records.push_back(
                make_representation(
                    pair,
                    pp
                )
            );
        }

        if (records.empty()) {
            std::cout
                << "PRIME="
                << p
                << " REPRESENTATIONS=0\n";

            continue;
        }

        std::sort(
            records.begin(),
            records.end(),
            [](const Representation& lhs,
               const Representation& rhs) {

                if (lhs.M != rhs.M) {
                    return lhs.M < rhs.M;
                }

                if (lhs.j != rhs.j) {
                    return lhs.j < rhs.j;
                }

                if (lhs.k1 != rhs.k1) {
                    return lhs.k1 < rhs.k1;
                }

                return lhs.k2 < rhs.k2;
            }
        );

        const Representation& first =
            records.front();

        const u64 u =
            first.j *
                static_cast<u64>(first.k1) +
            1;

        const u64 v =
            first.j *
                static_cast<u64>(first.k2) +
            1;

        std::cout
            << "PRIME="
            << p
            << " REPRESENTATIONS="
            << records.size()
            << " MIN_M="
            << first.M
            << " MIN_J="
            << first.j
            << " K1="
            << first.k1
            << " K2="
            << first.k2
            << " M1="
            << first.m1
            << " M2="
            << first.m2
            << " U="
            << u
            << " V="
            << v
            << " JP_PLUS_1="
            << first.j *
                   static_cast<u64>(p) +
               1
            << "\n";
    }

    if (first_found) {

        const u64 u =
            first_rep.j *
                static_cast<u64>(
                    first_rep.k1
                ) +
            1;

        const u64 v =
            first_rep.j *
                static_cast<u64>(
                    first_rep.k2
                ) +
            1;

        std::cout
            << "FIRST_SPECIAL_PRIME="
            << FIRST_PRIME
            << "\n";

        std::cout
            << "FIRST_SPECIAL_K1="
            << first_rep.k1
            << "\n";

        std::cout
            << "FIRST_SPECIAL_K2="
            << first_rep.k2
            << "\n";

        std::cout
            << "FIRST_SPECIAL_M1="
            << first_rep.m1
            << "\n";

        std::cout
            << "FIRST_SPECIAL_M2="
            << first_rep.m2
            << "\n";

        std::cout
            << "FIRST_SPECIAL_J="
            << first_rep.j
            << "\n";

        std::cout
            << "FIRST_SPECIAL_U="
            << u
            << "\n";

        std::cout
            << "FIRST_SPECIAL_V="
            << v
            << "\n";

        std::cout
            << "FIRST_SPECIAL_PRODUCT="
            << u * v
            << "\n";

        std::cout
            << "FIRST_SPECIAL_JP_PLUS_1="
            << first_rep.j *
                   static_cast<u64>(
                       FIRST_PRIME
                   ) +
               1
            << "\n";
    }

    const auto end =
        std::chrono::steady_clock::now();

    const double elapsed_ms =
        std::chrono::duration<double, std::milli>(
            end - start
        ).count();

    std::cout
        << "ELAPSED_TIME_MS="
        << elapsed_ms
        << "\n";

    std::cout
        << "FINISHED EXPERIMENT 488\n";

    return 0;
}
