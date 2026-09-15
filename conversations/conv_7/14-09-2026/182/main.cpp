#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

using u64 = std::uint64_t;

constexpr int EXPERIMENT = 475;

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

struct Candidate {
    u64 r = 0;
    Representation rep;

    u64 lcm_m = 0;
    u64 difference = 0;

    bool forced_zero = false;
};

struct SearchResult {
    bool found = false;

    u64 examined = 0;
    u64 candidate = 0;
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
    u64 r
) {
    if (p.m == 0) {
        return false;
    }

    const u64 value =
        (
            (p.a % p.m) *
            (r % p.m)
        ) %
        p.m;

    return value ==
           (n % p.m);
}

bool two_sided_match(
    u64 n,
    const Representation& p,
    const Representation& r
) {
    if (!one_sided_match(
            n,
            p,
            r.prime
        )) {
        return false;
    }

    if (r.m == 0) {
        return false;
    }

    const u64 value =
        (
            (r.a % r.m) *
            (p.prime % r.m)
        ) %
        r.m;

    return value ==
           (n % r.m);
}

Candidate build_candidate(
    const PrimeCase& c,
    const Representation& p,
    const Representation& r
) {
    Candidate candidate;

    candidate.r = r.prime;
    candidate.rep = r;

    const u64 g =
        std::gcd(
            p.m,
            r.m
        );

    candidate.lcm_m =
        (p.m / g) * r.m;

    const u64 product =
        p.prime * r.prime;

    candidate.difference =
        c.n >= product
            ? c.n - product
            : product - c.n;

    candidate.forced_zero =
        candidate.difference <
            candidate.lcm_m &&
        candidate.difference %
            candidate.lcm_m == 0;

    return candidate;
}

SearchResult search_candidates(
    std::vector<Candidate> candidates,
    int mode
) {
    /*
     * mode:
     *
     * 0 = descending r
     * 1 = ascending r
     * 2 = descending LCM
     * 3 = descending LCM/(difference+1)
     */
    if (mode == 0) {
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](
                const Candidate& a,
                const Candidate& b
            ) {
                return a.r > b.r;
            }
        );
    }

    if (mode == 1) {
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](
                const Candidate& a,
                const Candidate& b
            ) {
                return a.r < b.r;
            }
        );
    }

    if (mode == 2) {
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](
                const Candidate& a,
                const Candidate& b
            ) {
                if (a.lcm_m != b.lcm_m) {
                    return a.lcm_m > b.lcm_m;
                }

                return a.r > b.r;
            }
        );
    }

    if (mode == 3) {
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](
                const Candidate& a,
                const Candidate& b
            ) {
                /*
                 * Compare:
                 *
                 * LCM_a / (D_a+1)
                 *
                 * without floating point:
                 *
                 * LCM_a*(D_b+1)
                 * >
                 * LCM_b*(D_a+1)
                 */
                const u64 left =
                    a.lcm_m *
                    (b.difference + 1);

                const u64 right =
                    b.lcm_m *
                    (a.difference + 1);

                if (left != right) {
                    return left > right;
                }

                return a.r > b.r;
            }
        );
    }

    SearchResult result;

    for (
        std::size_t i = 0;
        i < candidates.size();
        ++i
    ) {
        ++result.examined;

        if (candidates[i].forced_zero) {
            result.found = true;
            result.candidate =
                candidates[i].r;
            return result;
        }
    }

    return result;
}

void main_experiment() {
    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << '\n';

    std::mt19937_64 rng(
        0x47520260915ULL
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

    const char* mode_names[4] = {
        "DESCENDING_R",
        "ASCENDING_R",
        "DESCENDING_LCM",
        "DESCENDING_LCM_OVER_DIFFERENCE"
    };

    u64 hits[4] = {0, 0, 0, 0};
    u64 total_examined[4] = {0, 0, 0, 0};
    u64 max_examined[4] = {0, 0, 0, 0};
    u64 min_examined[4] = {
        UINT64_MAX,
        UINT64_MAX,
        UINT64_MAX,
        UINT64_MAX
    };

    u64 exact_rank_1[4] = {0, 0, 0, 0};
    u64 rank_le_2[4] = {0, 0, 0, 0};
    u64 rank_le_5[4] = {0, 0, 0, 0};
    u64 rank_le_10[4] = {0, 0, 0, 0};

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

        if (!rp.found) {
            continue;
        }

        std::vector<Candidate> candidates;

        for (int prime : primes) {
            if (
                prime < PRIME_MIN ||
                prime > PRIME_LIMIT
            ) {
                continue;
            }

            const u64 r =
                static_cast<u64>(prime);

            if (r == c.p) {
                continue;
            }

            if (
                !one_sided_match(
                    c.n,
                    rp,
                    r
                )
            ) {
                continue;
            }

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

            candidates.push_back(
                build_candidate(
                    c,
                    rp,
                    rr
                )
            );
        }

        bool true_survived = false;

        for (
            const Candidate& candidate :
            candidates
        ) {
            if (candidate.r == c.q) {
                true_survived = true;
                break;
            }
        }

        if (!true_survived) {
            std::cerr
                << "ERROR: true q missing from survivors"
                << '\n';

            continue;
        }

        if (!first_example) {
            std::cout
                << "FIRST_N="
                << c.n
                << '\n';

            std::cout
                << "FIRST_P="
                << c.p
                << '\n';

            std::cout
                << "FIRST_Q="
                << c.q
                << '\n';

            std::cout
                << "FIRST_SURVIVOR_COUNT="
                << candidates.size()
                << '\n';

            first_example = true;
        }

        for (int mode = 0; mode < 4; ++mode) {
            const SearchResult result =
                search_candidates(
                    candidates,
                    mode
                );

            if (!result.found) {
                std::cerr
                    << "ERROR: certificate not found"
                    << '\n';

                continue;
            }

            ++hits[mode];

            total_examined[mode] +=
                result.examined;

            max_examined[mode] =
                std::max(
                    max_examined[mode],
                    result.examined
                );

            min_examined[mode] =
                std::min(
                    min_examined[mode],
                    result.examined
                );

            if (result.examined == 1) {
                ++exact_rank_1[mode];
            }

            if (result.examined <= 2) {
                ++rank_le_2[mode];
            }

            if (result.examined <= 5) {
                ++rank_le_5[mode];
            }

            if (result.examined <= 10) {
                ++rank_le_10[mode];
            }
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

    std::cout
        << "CASE_COUNT="
        << CASE_COUNT
        << '\n';

    for (int mode = 0; mode < 4; ++mode) {
        const double avg_examined =
            hits[mode] == 0
                ? 0.0
                : static_cast<double>(
                    static_cast<long double>(
                        total_examined[mode]
                    ) /
                    static_cast<long double>(
                        hits[mode]
                    )
                );

        std::cout
            << mode_names[mode]
            << "_HITS="
            << hits[mode]
            << '\n';

        std::cout
            << mode_names[mode]
            << "_AVG_EXAMINED="
            << avg_examined
            << '\n';

        std::cout
            << mode_names[mode]
            << "_MIN_EXAMINED="
            << (
                min_examined[mode] ==
                    UINT64_MAX
                    ? 0
                    : min_examined[mode]
            )
            << '\n';

        std::cout
            << mode_names[mode]
            << "_MAX_EXAMINED="
            << max_examined[mode]
            << '\n';

        std::cout
            << mode_names[mode]
            << "_RANK_1="
            << exact_rank_1[mode]
            << '\n';

        std::cout
            << mode_names[mode]
            << "_RANK_LE_2="
            << rank_le_2[mode]
            << '\n';

        std::cout
            << mode_names[mode]
            << "_RANK_LE_5="
            << rank_le_5[mode]
            << '\n';

        std::cout
            << mode_names[mode]
            << "_RANK_LE_10="
            << rank_le_10[mode]
            << '\n';
    }

    std::cout
        << "FINISHED EXPERIMENT "
        << EXPERIMENT
        << '\n';
}

int main() {
    main_experiment();
    return 0;
}
