#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 474;

constexpr int PRIME_LIMIT = 100000;
constexpr int PRIME_MIN = 10000;

constexpr int CASE_COUNT = 500;

constexpr int K_LIMIT = 1000;
constexpr int M_LIMIT = 7;

struct PrimeCase {
    u64 p;
    u64 q;
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
};

struct CandidateInfo {
    u64 lcm_m = 0;
    u64 difference = 0;
    u64 quotient = 0;
    u64 remainder = 0;

    bool lcm_divides_difference = false;
    bool difference_below_lcm = false;
    bool forced_zero = false;
};

std::vector<int> generate_primes(int limit) {
    std::vector<bool> composite(
        limit + 1,
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (
            static_cast<std::int64_t>(i) * i <=
            limit
        ) {
            for (
                int j = i * i;
                j <= limit;
                j += i
            ) {
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
        if (
            prime >= PRIME_MIN &&
            prime <= PRIME_LIMIT
        ) {
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
        u64 p =
            static_cast<u64>(
                candidates[dist(rng)]
            );

        u64 q =
            static_cast<u64>(
                candidates[dist(rng)]
            );

        while (q == p) {
            q =
                static_cast<u64>(
                    candidates[dist(rng)]
                );
        }

        if (p > q) {
            std::swap(p, q);
        }

        PrimeCase c;
        c.p = p;
        c.q = q;
        c.n = p * q;

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
             * Determinant-one condition:
             *
             *     m2*k1 - m1*k2 = 1
             *
             * with
             *
             *     1 <= m1,m2 <= M_LIMIT.
             */
            for (int m1 = 1; m1 <= M_LIMIT; ++m1) {
                const int numerator =
                    1 + m1 * k2;

                if (numerator % k1 != 0) {
                    continue;
                }

                const int m2 =
                    numerator / k1;

                if (
                    m2 >= 1 &&
                    m2 <= M_LIMIT
                ) {
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
        const u64 k1 =
            static_cast<u64>(pair.k1);

        const u64 k2 =
            static_cast<u64>(pair.k2);

        const u64 a =
            k1 + k2;

        const u64 m =
            k1 * k2;

        if (a > prime) {
            continue;
        }

        const u64 remainder =
            prime - a;

        if (remainder % m != 0) {
            continue;
        }

        const u64 j =
            remainder / m;

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
    }

    return best;
}

bool one_sided_match(
    u64 n,
    const Representation& p,
    u64 candidate
) {
    if (p.m == 0) {
        return false;
    }

    const u64 a_mod =
        p.a % p.m;

    const u64 candidate_mod =
        candidate % p.m;

    const u64 n_mod =
        n % p.m;

    const u64 product_mod =
        (
            a_mod *
            candidate_mod
        ) %
        p.m;

    return product_mod == n_mod;
}

bool two_sided_match(
    u64 n,
    const Representation& p,
    const Representation& r
) {
    if (
        p.m == 0 ||
        r.m == 0
    ) {
        return false;
    }

    if (
        !one_sided_match(
            n,
            p,
            r.prime
        )
    ) {
        return false;
    }

    const u64 a_mod =
        r.a % r.m;

    const u64 p_mod =
        p.prime % r.m;

    const u64 n_mod =
        n % r.m;

    const u64 product_mod =
        (
            a_mod *
            p_mod
        ) %
        r.m;

    return product_mod == n_mod;
}

CandidateInfo analyze_candidate(
    const PrimeCase& c,
    const Representation& p,
    const Representation& r
) {
    CandidateInfo result;

    const u64 g =
        std::gcd(p.m, r.m);

    result.lcm_m =
        (p.m / g) * r.m;

    if (c.n >= p.prime * r.prime) {
        result.difference =
            c.n - p.prime * r.prime;
    } else {
        result.difference =
            p.prime * r.prime - c.n;
    }

    if (result.lcm_m != 0) {
        result.quotient =
            result.difference /
            result.lcm_m;

        result.remainder =
            result.difference %
            result.lcm_m;

        result.lcm_divides_difference =
            result.remainder == 0;

        result.difference_below_lcm =
            result.difference <
            result.lcm_m;

        result.forced_zero =
            result.lcm_divides_difference &&
            result.difference_below_lcm;
    }

    return result;
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
        << rep.a +
           rep.j * rep.m
        << '\n';
}

void print_candidate(
    const char* label,
    const CandidateInfo& info
) {
    std::cout
        << label
        << "_LCM_M="
        << info.lcm_m
        << '\n';

    std::cout
        << label
        << "_DIFFERENCE="
        << info.difference
        << '\n';

    std::cout
        << label
        << "_QUOTIENT="
        << info.quotient
        << '\n';

    std::cout
        << label
        << "_REMAINDER="
        << info.remainder
        << '\n';

    std::cout
        << label
        << "_LCM_DIVIDES="
        << (info.lcm_divides_difference ? 1 : 0)
        << '\n';

    std::cout
        << label
        << "_DIFFERENCE_BELOW_LCM="
        << (info.difference_below_lcm ? 1 : 0)
        << '\n';

    std::cout
        << label
        << "_FORCED_ZERO="
        << (info.forced_zero ? 1 : 0)
        << '\n';
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47420260915ULL
    );

    const auto total_start =
        std::chrono::steady_clock::now();

    const std::vector<int> primes =
        generate_primes(
            PRIME_LIMIT
        );

    const std::vector<PrimeCase> cases =
        generate_cases(
            primes,
            rng
        );

    const std::vector<CRTPair> pairs =
        build_determinant_one_pairs();

    std::vector<Representation> cache(
        PRIME_LIMIT + 1
    );

    for (int prime : primes) {
        if (prime < PRIME_MIN) {
            continue;
        }

        cache[prime] =
            find_min_j_representation(
                static_cast<u64>(prime),
                pairs
            );
    }

    std::cout
        << "PAIR_COUNT="
        << pairs.size()
        << '\n';

    std::cout
        << "PRIME_POOL_SIZE="
        << primes.size()
        << '\n';

    /*
     * Candidate totals.
     */
    u64 valid_cases = 0;

    u64 one_sided_total = 0;
    u64 two_sided_total = 0;

    u64 two_sided_lcm_divides = 0;
    u64 two_sided_forced_zero = 0;

    /*
     * For false positives after the two-sided filter.
     */
    u64 false_positive_two_sided = 0;
    u64 false_positive_lcm_divides = 0;
    u64 false_positive_forced_zero = 0;

    /*
     * Unique counts.
     */
    u64 two_sided_unique = 0;
    u64 forced_zero_unique = 0;

    /*
     * Candidate distribution.
     */
    u64 two_sided_1 = 0;
    u64 two_sided_2_5 = 0;
    u64 two_sided_6_10 = 0;
    u64 two_sided_11_100 = 0;
    u64 two_sided_gt100 = 0;

    /*
     * LCM and difference statistics.
     */
    long double sum_lcm =
        0.0L;

    long double sum_difference =
        0.0L;

    long double sum_false_difference =
        0.0L;

    long double sum_false_lcm =
        0.0L;

    long double sum_quotient =
        0.0L;

    u64 false_quotient_count = 0;

    u64 min_lcm =
        UINT64_MAX;

    u64 max_lcm = 0;

    /*
     * Timing.
     */
    double total_candidate_time_ms = 0.0;

    bool first_example = false;

    for (
        int case_index = 0;
        case_index < CASE_COUNT;
        ++case_index
    ) {
        const PrimeCase& c =
            cases[case_index];

        const Representation& rp =
            cache[c.p];

        const Representation& rq =
            cache[c.q];

        if (!rp.found || !rq.found) {
            continue;
        }

        ++valid_cases;

        /*
         * True q must survive both conditions.
         */
        const bool true_two_sided =
            two_sided_match(
                c.n,
                rp,
                rq
            );

        if (!true_two_sided) {
            std::cerr
                << "ERROR: true q failed two-sided filter"
                << '\n';

            continue;
        }

        const auto scan_start =
            std::chrono::steady_clock::now();

        u64 one_sided_count = 0;
        u64 two_sided_count = 0;

        u64 case_lcm_divides = 0;
        u64 case_forced_zero = 0;

        u64 case_false_positive = 0;
        u64 case_false_lcm = 0;
        u64 case_false_forced = 0;

        /*
         * Scan only the prime pool.
         */
        for (int prime : primes) {
            if (
                prime < PRIME_MIN ||
                prime > PRIME_LIMIT
            ) {
                continue;
            }

            const u64 r_value =
                static_cast<u64>(prime);

            if (r_value == c.p) {
                continue;
            }

            if (
                !one_sided_match(
                    c.n,
                    rp,
                    r_value
                )
            ) {
                continue;
            }

            ++one_sided_count;

            const Representation& rr =
                cache[prime];

            if (!rr.found) {
                continue;
            }

            if (
                !two_sided_match(
                    c.n,
                    rp,
                    rr
                )
            ) {
                continue;
            }

            ++two_sided_count;

            const CandidateInfo info =
                analyze_candidate(
                    c,
                    rp,
                    rr
                );

            if (info.lcm_divides_difference) {
                ++case_lcm_divides;
            }

            if (info.forced_zero) {
                ++case_forced_zero;
            }

            sum_lcm +=
                static_cast<long double>(
                    info.lcm_m
                );

            sum_difference +=
                static_cast<long double>(
                    info.difference
                );

            min_lcm =
                std::min(
                    min_lcm,
                    info.lcm_m
                );

            max_lcm =
                std::max(
                    max_lcm,
                    info.lcm_m
                );

            /*
             * q is the unique candidate with
             * difference == 0.
             */
            if (
                r_value != c.q
            ) {
                ++case_false_positive;

                sum_false_difference +=
                    static_cast<long double>(
                        info.difference
                    );

                sum_false_lcm +=
                    static_cast<long double>(
                        info.lcm_m
                    );

                if (info.lcm_divides_difference) {
                    ++case_false_lcm;
                }

                if (info.forced_zero) {
                    ++case_false_forced;
                }

                if (info.lcm_divides_difference) {
                    sum_quotient +=
                        static_cast<long double>(
                            info.quotient
                        );

                    ++false_quotient_count;
                }
            }

            /*
             * First candidate example.
             */
            if (!first_example) {
                std::cout
                    << "FIRST_N="
                    << c.n
                    << '\n';

                print_representation(
                    "FIRST_P",
                    rp
                );

                print_representation(
                    "FIRST_Q",
                    rq
                );

                std::cout
                    << "FIRST_CANDIDATE="
                    << r_value
                    << '\n';

                print_candidate(
                    "FIRST_CANDIDATE_INFO",
                    info
                );

                first_example = true;
            }
        }

        const auto scan_end =
            std::chrono::steady_clock::now();

        total_candidate_time_ms +=
            std::chrono::duration<
                double,
                std::milli
            >(
                scan_end - scan_start
            ).count();

        one_sided_total +=
            one_sided_count;

        two_sided_total +=
            two_sided_count;

        two_sided_lcm_divides +=
            case_lcm_divides;

        two_sided_forced_zero +=
            case_forced_zero;

        false_positive_two_sided +=
            case_false_positive;

        false_positive_lcm_divides +=
            case_false_lcm;

        false_positive_forced_zero +=
            case_false_forced;

        if (two_sided_count == 1) {
            ++two_sided_unique;
        }

        if (case_forced_zero == 1) {
            ++forced_zero_unique;
        }

        if (two_sided_count == 1) {
            ++two_sided_1;
        } else if (two_sided_count <= 5) {
            ++two_sided_2_5;
        } else if (two_sided_count <= 10) {
            ++two_sided_6_10;
        } else if (two_sided_count <= 100) {
            ++two_sided_11_100;
        } else {
            ++two_sided_gt100;
        }

        if (
            (case_index + 1) % 100 ==
            0
        ) {
            std::cout
                << "PROGRESS="
                << (case_index + 1)
                << "/"
                << CASE_COUNT
                << '\n';
        }
    }

    const auto total_end =
        std::chrono::steady_clock::now();

    const double total_time_ms =
        std::chrono::duration<
            double,
            std::milli
        >(
            total_end - total_start
        ).count();

    const double avg_one_sided =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    one_sided_total
                ) /
                static_cast<long double>(
                    valid_cases
                )
            );

    const double avg_two_sided =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    two_sided_total
                ) /
                static_cast<long double>(
                    valid_cases
                )
            );

    const double avg_case_survival =
        one_sided_total == 0
            ? 0.0
            : static_cast<double>(
                static_cast<long double>(
                    two_sided_total
                ) /
                static_cast<long double>(
                    one_sided_total
                )
            );

    const double avg_lcm =
        two_sided_total == 0
            ? 0.0
            : static_cast<double>(
                sum_lcm /
                static_cast<long double>(
                    two_sided_total
                )
            );

    const double avg_difference =
        two_sided_total == 0
            ? 0.0
            : static_cast<double>(
                sum_difference /
                static_cast<long double>(
                    two_sided_total
                )
            );

    const double avg_false_difference =
        false_positive_two_sided == 0
            ? 0.0
            : static_cast<double>(
                sum_false_difference /
                static_cast<long double>(
                    false_positive_two_sided
                )
            );

    const double avg_false_lcm =
        false_positive_two_sided == 0
            ? 0.0
            : static_cast<double>(
                sum_false_lcm /
                static_cast<long double>(
                    false_positive_two_sided
                )
            );

    const double avg_false_quotient =
        false_quotient_count == 0
            ? 0.0
            : static_cast<double>(
                sum_quotient /
                static_cast<long double>(
                    false_quotient_count
                )
            );

    const double forced_zero_rate =
        two_sided_total == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  two_sided_forced_zero
              ) /
              static_cast<double>(
                  two_sided_total
              );

    const double lcm_divides_rate =
        two_sided_total == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  two_sided_lcm_divides
              ) /
              static_cast<double>(
                  two_sided_total
              );

    const double false_forced_rate =
        false_positive_two_sided == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  false_positive_forced_zero
              ) /
              static_cast<double>(
                  false_positive_two_sided
              );

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    std::cout
        << "VALID_CASES="
        << valid_cases
        << '\n';

    std::cout
        << "TRUE_TWO_SIDED="
        << valid_cases
        << '\n';

    std::cout
        << "TOTAL_ONE_SIDED_CANDIDATES="
        << one_sided_total
        << '\n';

    std::cout
        << "TOTAL_TWO_SIDED_CANDIDATES="
        << two_sided_total
        << '\n';

    std::cout
        << "AVG_ONE_SIDED_CANDIDATES="
        << avg_one_sided
        << '\n';

    std::cout
        << "AVG_TWO_SIDED_CANDIDATES="
        << avg_two_sided
        << '\n';

    std::cout
        << "AVG_TWO_SIDED_SURVIVAL_RATIO="
        << avg_case_survival
        << '\n';

    std::cout
        << "TWO_SIDED_LCM_DIVIDES="
        << two_sided_lcm_divides
        << '\n';

    std::cout
        << "TWO_SIDED_LCM_DIVIDES_PERCENT="
        << lcm_divides_rate
        << '\n';

    std::cout
        << "TWO_SIDED_FORCED_ZERO="
        << two_sided_forced_zero
        << '\n';

    std::cout
        << "FORCED_ZERO_PERCENT="
        << forced_zero_rate
        << '\n';

    std::cout
        << "TWO_SIDED_UNIQUE="
        << two_sided_unique
        << '\n';

    std::cout
        << "FORCED_ZERO_UNIQUE="
        << forced_zero_unique
        << '\n';

    std::cout
        << "TWO_SIDED_1="
        << two_sided_1
        << '\n';

    std::cout
        << "TWO_SIDED_2_5="
        << two_sided_2_5
        << '\n';

    std::cout
        << "TWO_SIDED_6_10="
        << two_sided_6_10
        << '\n';

    std::cout
        << "TWO_SIDED_11_100="
        << two_sided_11_100
        << '\n';

    std::cout
        << "TWO_SIDED_GT100="
        << two_sided_gt100
        << '\n';

    std::cout
        << "FALSE_POSITIVE_TWO_SIDED="
        << false_positive_two_sided
        << '\n';

    std::cout
        << "FALSE_POSITIVE_LCM_DIVIDES="
        << false_positive_lcm_divides
        << '\n';

    std::cout
        << "FALSE_POSITIVE_FORCED_ZERO="
        << false_positive_forced_zero
        << '\n';

    std::cout
        << "FALSE_FORCED_ZERO_PERCENT="
        << false_forced_rate
        << '\n';

    std::cout
        << "AVG_LCM_M="
        << avg_lcm
        << '\n';

    std::cout
        << "MIN_LCM_M="
        << (
            min_lcm == UINT64_MAX
                ? 0
                : min_lcm
        )
        << '\n';

    std::cout
        << "MAX_LCM_M="
        << max_lcm
        << '\n';

    std::cout
        << "AVG_DIFFERENCE="
        << avg_difference
        << '\n';

    std::cout
        << "AVG_FALSE_DIFFERENCE="
        << avg_false_difference
        << '\n';

    std::cout
        << "AVG_FALSE_LCM="
        << avg_false_lcm
        << '\n';

    std::cout
        << "FALSE_QUOTIENT_COUNT="
        << false_quotient_count
        << '\n';

    std::cout
        << "AVG_FALSE_QUOTIENT="
        << avg_false_quotient
        << '\n';

    std::cout
        << "CANDIDATE_SCAN_TIME_MS="
        << total_candidate_time_ms
        << '\n';

    std::cout
        << "TOTAL_TIME_MS="
        << total_time_ms
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
