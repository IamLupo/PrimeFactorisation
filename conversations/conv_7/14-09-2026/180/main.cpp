#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 473;

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

struct ResidueFilter {
    bool valid = false;

    u64 inverse_a = 0;
    u64 target_residue = 0;
    u64 modulus = 0;
};

struct CandidateStats {
    u64 one_sided = 0;
    u64 two_sided = 0;
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
             * m2*k1 - m1*k2 = 1
             *
             * with
             *
             * 1 <= m1,m2 <= M_LIMIT.
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

u64 extended_gcd(
    u64 a,
    u64 b,
    u64& x,
    u64& y
) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }

    u64 x1 = 0;
    u64 y1 = 0;

    const u64 g =
        extended_gcd(
            b,
            a % b,
            x1,
            y1
        );

    x = y1;

    /*
     * This helper is only used for the small
     * moduli encountered here. To avoid signed
     * arithmetic complications, use a separate
     * signed implementation below instead.
     */
    (void)x;
    (void)y;

    return g;
}

std::int64_t extended_gcd_signed(
    std::int64_t a,
    std::int64_t b,
    std::int64_t& x,
    std::int64_t& y
) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }

    std::int64_t x1 = 0;
    std::int64_t y1 = 0;

    const std::int64_t g =
        extended_gcd_signed(
            b,
            a % b,
            x1,
            y1
        );

    x = y1;
    y = x1 - (a / b) * y1;

    return g;
}

u64 modular_inverse(
    u64 a,
    u64 modulus,
    bool& valid
) {
    valid = false;

    if (modulus == 0) {
        return 0;
    }

    if (modulus == 1) {
        valid = true;
        return 0;
    }

    std::int64_t x = 0;
    std::int64_t y = 0;

    const std::int64_t g =
        extended_gcd_signed(
            static_cast<std::int64_t>(
                a % modulus
            ),
            static_cast<std::int64_t>(
                modulus
            ),
            x,
            y
        );

    if (g != 1) {
        return 0;
    }

    std::int64_t result =
        x %
        static_cast<std::int64_t>(
            modulus
        );

    if (result < 0) {
        result +=
            static_cast<std::int64_t>(
                modulus
            );
    }

    valid = true;

    return static_cast<u64>(result);
}

ResidueFilter build_residue_filter(
    const PrimeCase& c,
    const Representation& rp
) {
    ResidueFilter filter;

    filter.modulus = rp.m;

    if (rp.m == 0) {
        return filter;
    }

    if (std::gcd(rp.a, rp.m) != 1) {
        return filter;
    }

    bool valid = false;

    filter.inverse_a =
        modular_inverse(
            rp.a,
            rp.m,
            valid
        );

    if (!valid) {
        return filter;
    }

    filter.target_residue =
        (
            (c.n % rp.m) *
            filter.inverse_a
        ) %
        rp.m;

    filter.valid = true;

    return filter;
}

bool one_sided_match(
    u64 n,
    const Representation& p,
    u64 candidate
) {
    if (p.m == 0) {
        return false;
    }

    /*
     * M_p | N - a_p*r
     */
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

    /*
     * First side:
     *
     * M_p | N - a_p*r
     */
    if (!one_sided_match(
            n,
            p,
            r.prime
        )) {
        return false;
    }

    /*
     * Second side:
     *
     * M_r | N - a_r*p
     */
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

u64 count_one_sided_candidates(
    const std::vector<int>& primes,
    const PrimeCase& c,
    const Representation& p
) {
    u64 count = 0;

    for (int prime : primes) {
        if (
            prime < PRIME_MIN ||
            prime > PRIME_LIMIT
        ) {
            continue;
        }

        if (
            static_cast<u64>(prime) == c.p
        ) {
            continue;
        }

        if (
            one_sided_match(
                c.n,
                p,
                static_cast<u64>(prime)
            )
        ) {
            ++count;
        }
    }

    return count;
}

u64 count_two_sided_candidates(
    const std::vector<int>& primes,
    const PrimeCase& c,
    const Representation& p,
    const std::vector<Representation>& cache
) {
    u64 count = 0;

    for (int prime : primes) {
        if (
            prime < PRIME_MIN ||
            prime > PRIME_LIMIT
        ) {
            continue;
        }

        if (
            static_cast<u64>(prime) == c.p
        ) {
            continue;
        }

        const Representation& r =
            cache[prime];

        if (!r.found) {
            continue;
        }

        if (
            two_sided_match(
                c.n,
                p,
                r
            )
        ) {
            ++count;
        }
    }

    return count;
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

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47320260915ULL
    );

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

    u64 valid_cases = 0;

    u64 true_one_sided = 0;
    u64 true_two_sided = 0;

    /*
     * Candidate totals over the complete prime pool.
     */
    u64 one_sided_total = 0;
    u64 two_sided_total = 0;

    u64 min_one_sided = UINT64_MAX;
    u64 max_one_sided = 0;

    u64 min_two_sided = UINT64_MAX;
    u64 max_two_sided = 0;

    long double sum_one_sided = 0.0L;
    long double sum_two_sided = 0.0L;

    /*
     * Uniqueness.
     */
    u64 one_sided_unique = 0;
    u64 two_sided_unique = 0;

    /*
     * Reduction:
     *
     * two_sided / one_sided
     */
    long double sum_survival_ratio = 0.0L;

    /*
     * Does the second filter eliminate all
     * false positives while keeping q?
     */
    u64 cases_with_false_positive_after_two_sided = 0;

    /*
     * Nearby controls.
     */
    u64 nearby_controls = 0;
    u64 nearby_one_sided_hits = 0;
    u64 nearby_two_sided_hits = 0;

    /*
     * Candidate multiplicity distribution.
     */
    u64 two_sided_count_1 = 0;
    u64 two_sided_count_2_5 = 0;
    u64 two_sided_count_6_10 = 0;
    u64 two_sided_count_11_100 = 0;
    u64 two_sided_count_gt100 = 0;

    /*
     * First example.
     */
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
         * Verify that the true q passes both.
         */
        const bool true_first =
            one_sided_match(
                c.n,
                rp,
                c.q
            );

        const bool true_second =
            two_sided_match(
                c.n,
                rp,
                rq
            );

        if (true_first) {
            ++true_one_sided;
        }

        if (true_second) {
            ++true_two_sided;
        }

        /*
         * Count complete candidate sets.
         */
        const u64 one_sided_count =
            count_one_sided_candidates(
                primes,
                c,
                rp
            );

        const u64 two_sided_count =
            count_two_sided_candidates(
                primes,
                c,
                rp,
                cache
            );

        one_sided_total +=
            one_sided_count;

        two_sided_total +=
            two_sided_count;

        min_one_sided =
            std::min(
                min_one_sided,
                one_sided_count
            );

        max_one_sided =
            std::max(
                max_one_sided,
                one_sided_count
            );

        min_two_sided =
            std::min(
                min_two_sided,
                two_sided_count
            );

        max_two_sided =
            std::max(
                max_two_sided,
                two_sided_count
            );

        sum_one_sided +=
            static_cast<long double>(
                one_sided_count
            );

        sum_two_sided +=
            static_cast<long double>(
                two_sided_count
            );

        if (one_sided_count == 1) {
            ++one_sided_unique;
        }

        if (two_sided_count == 1) {
            ++two_sided_unique;
        }

        if (
            two_sided_count >
            1
        ) {
            /*
             * q is one of the survivors.
             * Therefore all remaining candidates
             * after the second filter are false
             * positives except q.
             */
            ++cases_with_false_positive_after_two_sided;
        }

        if (one_sided_count != 0) {
            sum_survival_ratio +=
                static_cast<long double>(
                    two_sided_count
                ) /
                static_cast<long double>(
                    one_sided_count
                );
        }

        if (two_sided_count == 1) {
            ++two_sided_count_1;
        } else if (two_sided_count <= 5) {
            ++two_sided_count_2_5;
        } else if (two_sided_count <= 10) {
            ++two_sided_count_6_10;
        } else if (two_sided_count <= 100) {
            ++two_sided_count_11_100;
        } else {
            ++two_sided_count_gt100;
        }

        /*
         * Nearby controls.
         *
         * Use the same nearest-prime neighborhood
         * as Experiment 468-472.
         */
        auto lower =
            std::lower_bound(
                primes.begin(),
                primes.end(),
                static_cast<int>(c.q)
            );

        int right =
            static_cast<int>(
                lower - primes.begin()
            );

        int left = right - 1;

        int controls_seen = 0;

        while (
            controls_seen < 100 &&
            (
                left >= 0 ||
                right <
                    static_cast<int>(
                        primes.size()
                    )
            )
        ) {
            bool take_left;

            if (left < 0) {
                take_left = false;
            } else if (
                right >=
                static_cast<int>(
                    primes.size()
                )
            ) {
                take_left = true;
            } else {
                const u64 left_value =
                    static_cast<u64>(
                        primes[left]
                    );

                const u64 right_value =
                    static_cast<u64>(
                        primes[right]
                    );

                const u64 left_distance =
                    c.q > left_value
                        ? c.q - left_value
                        : left_value - c.q;

                const u64 right_distance =
                    c.q > right_value
                        ? c.q - right_value
                        : right_value - c.q;

                take_left =
                    left_distance <=
                    right_distance;
            }

            int candidate;

            if (take_left) {
                candidate = primes[left];
                --left;
            } else {
                candidate = primes[right];
                ++right;
            }

            if (
                static_cast<u64>(candidate) ==
                    c.p ||
                static_cast<u64>(candidate) ==
                    c.q
            ) {
                continue;
            }

            if (candidate < PRIME_MIN) {
                continue;
            }

            ++controls_seen;
            ++nearby_controls;

            const Representation& rr =
                cache[candidate];

            if (!rr.found) {
                continue;
            }

            if (
                one_sided_match(
                    c.n,
                    rp,
                    static_cast<u64>(candidate)
                )
            ) {
                ++nearby_one_sided_hits;
            }

            if (
                two_sided_match(
                    c.n,
                    rp,
                    rr
                )
            ) {
                ++nearby_two_sided_hits;
            }
        }

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
                << "FIRST_TRUE_ONE_SIDED="
                << (true_first ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_TRUE_TWO_SIDED="
                << (true_second ? 1 : 0)
                << '\n';

            std::cout
                << "FIRST_ONE_SIDED_COUNT="
                << one_sided_count
                << '\n';

            std::cout
                << "FIRST_TWO_SIDED_COUNT="
                << two_sided_count
                << '\n';

            std::cout
                << "FIRST_TWO_SIDED_FALSE_POSITIVES="
                << (
                    two_sided_count > 0
                        ? two_sided_count - 1
                        : 0
                )
                << '\n';

            first_example = true;
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

    const double avg_one_sided =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_one_sided /
                static_cast<long double>(
                    valid_cases
                )
            );

    const double avg_two_sided =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_two_sided /
                static_cast<long double>(
                    valid_cases
                )
            );

    const double avg_survival_ratio =
        valid_cases == 0
            ? 0.0
            : static_cast<double>(
                sum_survival_ratio /
                static_cast<long double>(
                    valid_cases
                )
            );

    const double true_one_sided_rate =
        valid_cases == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  true_one_sided
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double true_two_sided_rate =
        valid_cases == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  true_two_sided
              ) /
              static_cast<double>(
                  valid_cases
              );

    const double nearby_one_sided_rate =
        nearby_controls == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  nearby_one_sided_hits
              ) /
              static_cast<double>(
                  nearby_controls
              );

    const double nearby_two_sided_rate =
        nearby_controls == 0
            ? 0.0
            : 100.0 *
              static_cast<double>(
                  nearby_two_sided_hits
              ) /
              static_cast<double>(
                  nearby_controls
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
        << "TRUE_ONE_SIDED="
        << true_one_sided
        << '\n';

    std::cout
        << "TRUE_TWO_SIDED="
        << true_two_sided
        << '\n';

    std::cout
        << "TRUE_ONE_SIDED_RATE="
        << true_one_sided_rate
        << '\n';

    std::cout
        << "TRUE_TWO_SIDED_RATE="
        << true_two_sided_rate
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
        << "MIN_ONE_SIDED_CANDIDATES="
        << (
            min_one_sided == UINT64_MAX
                ? 0
                : min_one_sided
        )
        << '\n';

    std::cout
        << "MAX_ONE_SIDED_CANDIDATES="
        << max_one_sided
        << '\n';

    std::cout
        << "MIN_TWO_SIDED_CANDIDATES="
        << (
            min_two_sided == UINT64_MAX
                ? 0
                : min_two_sided
        )
        << '\n';

    std::cout
        << "MAX_TWO_SIDED_CANDIDATES="
        << max_two_sided
        << '\n';

    std::cout
        << "ONE_SIDED_UNIQUE="
        << one_sided_unique
        << '\n';

    std::cout
        << "TWO_SIDED_UNIQUE="
        << two_sided_unique
        << '\n';

    std::cout
        << "CASES_WITH_FALSE_POSITIVE_AFTER_TWO_SIDED="
        << cases_with_false_positive_after_two_sided
        << '\n';

    std::cout
        << "AVG_TWO_SIDED_SURVIVAL_RATIO="
        << avg_survival_ratio
        << '\n';

    std::cout
        << "TWO_SIDED_COUNT_1="
        << two_sided_count_1
        << '\n';

    std::cout
        << "TWO_SIDED_COUNT_2_5="
        << two_sided_count_2_5
        << '\n';

    std::cout
        << "TWO_SIDED_COUNT_6_10="
        << two_sided_count_6_10
        << '\n';

    std::cout
        << "TWO_SIDED_COUNT_11_100="
        << two_sided_count_11_100
        << '\n';

    std::cout
        << "TWO_SIDED_COUNT_GT100="
        << two_sided_count_gt100
        << '\n';

    std::cout
        << "NEARBY_CONTROL_PAIRS="
        << nearby_controls
        << '\n';

    std::cout
        << "NEARBY_ONE_SIDED_HITS="
        << nearby_one_sided_hits
        << '\n';

    std::cout
        << "NEARBY_TWO_SIDED_HITS="
        << nearby_two_sided_hits
        << '\n';

    std::cout
        << "NEARBY_ONE_SIDED_HIT_RATE="
        << nearby_one_sided_rate
        << '\n';

    std::cout
        << "NEARBY_TWO_SIDED_HIT_RATE="
        << nearby_two_sided_rate
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
