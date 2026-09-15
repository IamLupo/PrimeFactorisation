#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using u64 = std::uint64_t;

struct Witness {
    int m = 1;
    u64 k = 1;
    u64 t = 0;
    int sign = +1;
};

struct Segment {
    u64 begin = 0;
    u64 end = 0;
    Witness witness;
};

struct PairInfo {
    u64 r = 0;
    int source_sign = +1;
    u64 A = 0;
    u64 k = 0;
    u64 t = 0;
};

struct Permutation {
    int map[5] = {1, 2, 3, 4, 6};
};

static std::vector<u64> generate_primes(
    int limit
) {
    std::vector<bool> sieve(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    sieve[0] = false;
    sieve[1] = false;

    for (
        int p = 2;
        static_cast<long long>(p) * p <= limit;
        ++p
    ) {
        if (!sieve[p]) {
            continue;
        }

        for (
            int x = p * p;
            x <= limit;
            x += p
        ) {
            sieve[x] = false;
        }
    }

    std::vector<u64> primes;

    for (int x = 2; x <= limit; ++x) {
        if (sieve[x]) {
            primes.push_back(
                static_cast<u64>(x)
            );
        }
    }

    return primes;
}

static bool divides_stream(
    u64 r,
    int m,
    u64 k,
    int sign
) {
    const __int128 value =
        static_cast<__int128>(m) *
        static_cast<__int128>(r) +
        static_cast<__int128>(sign);

    return
        value > 0 &&
        value %
            static_cast<__int128>(k) == 0;
}

static Witness make_witness(
    u64 r,
    int m,
    u64 k,
    int sign
) {
    Witness w;

    w.m = m;
    w.k = k;
    w.sign = sign;

    const __int128 value =
        static_cast<__int128>(m) *
        static_cast<__int128>(r) +
        static_cast<__int128>(sign);

    w.t =
        static_cast<u64>(
            value /
            static_cast<__int128>(k)
        );

    return w;
}

static bool better_normalized(
    const Witness& a,
    const Witness& b
) {
    const __int128 lhs =
        static_cast<__int128>(a.k) *
        static_cast<__int128>(b.m);

    const __int128 rhs =
        static_cast<__int128>(b.k) *
        static_cast<__int128>(a.m);

    if (lhs != rhs) {
        return lhs > rhs;
    }

    if (a.m != b.m) {
        return a.m < b.m;
    }

    if (a.k != b.k) {
        return a.k > b.k;
    }

    return a.sign > b.sign;
}

static Witness strongest_allowed(
    const std::vector<Witness>& best
) {
    Witness result = best[1];

    for (int m : {2, 3, 4, 6}) {
        const Witness& candidate =
            best[
                static_cast<std::size_t>(m)
            ];

        if (
            better_normalized(
                candidate,
                result
            )
        ) {
            result = candidate;
        }
    }

    return result;
}

static u64 margin_for(
    u64 k,
    const Witness& allowed
) {
    const __int128 value =
        static_cast<__int128>(k) *
        static_cast<__int128>(allowed.m) -
        static_cast<__int128>(5) *
        static_cast<__int128>(allowed.k);

    if (value <= 0) {
        return 0;
    }

    return static_cast<u64>(value);
}

static std::vector<Witness> build_envelope(
    u64 r,
    int K_LIMIT
) {
    std::vector<Witness> best_by_m(7);

    for (int m = 1; m <= 6; ++m) {
        best_by_m[
            static_cast<std::size_t>(m)
        ] =
            make_witness(
                r,
                m,
                1,
                +1
            );
    }

    std::vector<Witness> envelope(
        static_cast<std::size_t>(
            K_LIMIT + 1
        )
    );

    for (
        u64 K = 1;
        K <= static_cast<u64>(K_LIMIT);
        ++K
    ) {
        for (int m : {1, 2, 3, 4, 6}) {
            const std::size_t mi =
                static_cast<std::size_t>(m);

            if (
                divides_stream(
                    r,
                    m,
                    K,
                    +1
                )
            ) {
                const Witness candidate =
                    make_witness(
                        r,
                        m,
                        K,
                        +1
                    );

                if (
                    candidate.k >
                    best_by_m[mi].k
                ) {
                    best_by_m[mi] =
                        candidate;
                }
            }

            if (
                divides_stream(
                    r,
                    m,
                    K,
                    -1
                )
            ) {
                const Witness candidate =
                    make_witness(
                        r,
                        m,
                        K,
                        -1
                    );

                if (
                    candidate.k >
                    best_by_m[mi].k
                ) {
                    best_by_m[mi] =
                        candidate;
                }
            }
        }

        envelope[
            static_cast<std::size_t>(K)
        ] =
            strongest_allowed(
                best_by_m
            );
    }

    return envelope;
}

static std::vector<PairInfo>
find_exceptional_pairs(
    u64 r,
    int source_sign,
    int K_LIMIT
) {
    const __int128 value =
        static_cast<__int128>(5) *
        static_cast<__int128>(r) +
        static_cast<__int128>(source_sign);

    const u64 A =
        static_cast<u64>(value);

    const std::vector<Witness> envelope =
        build_envelope(
            r,
            K_LIMIT
        );

    std::vector<u64> winners;

    for (
        u64 k = 1;
        k <= static_cast<u64>(K_LIMIT);
        ++k
    ) {
        if (A % k != 0) {
            continue;
        }

        const Witness& allowed =
            envelope[
                static_cast<std::size_t>(k)
            ];

        if (
            margin_for(
                k,
                allowed
            ) > 0
        ) {
            winners.push_back(k);
        }
    }

    std::vector<PairInfo> pairs;

    for (u64 k : winners) {
        const u64 t = A / k;

        if (k >= t) {
            continue;
        }

        if (
            t > static_cast<u64>(K_LIMIT)
        ) {
            continue;
        }

        PairInfo pair;

        pair.r = r;
        pair.source_sign = source_sign;
        pair.A = A;
        pair.k = k;
        pair.t = t;

        pairs.push_back(pair);
    }

    return pairs;
}

static std::vector<int> compress_path(
    const std::vector<int>& raw
) {
    std::vector<int> path;

    for (int m : raw) {
        if (
            path.empty() ||
            path.back() != m
        ) {
            path.push_back(m);
        }
    }

    return path;
}

static int multiplier_index(int m) {
    switch (m) {
        case 1:
            return 0;
        case 2:
            return 1;
        case 3:
            return 2;
        case 4:
            return 3;
        case 6:
            return 4;
        default:
            return -1;
    }
}

static int permutation_apply(
    const Permutation& permutation,
    int m
) {
    const int index =
        multiplier_index(m);

    if (index < 0) {
        return -1;
    }

    return permutation.map[index];
}

static bool permutation_is_involution(
    const Permutation& permutation
) {
    const int values[5] =
        {1, 2, 3, 4, 6};

    for (int m : values) {
        const int first =
            permutation_apply(
                permutation,
                m
            );

        const int second =
            permutation_apply(
                permutation,
                first
            );

        if (second != m) {
            return false;
        }
    }

    return true;
}

static bool path_matches_permutation(
    const std::vector<int>& path,
    const Permutation& permutation
) {
    if (path.empty()) {
        return true;
    }

    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        const int reverse_value =
            path[
                path.size() - 1 - i
            ];

        const int mapped =
            permutation_apply(
                permutation,
                reverse_value
            );

        if (mapped != path[i]) {
            return false;
        }
    }

    return true;
}

static bool path_matches_identity(
    const std::vector<int>& path
) {
    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        if (
            path[i] !=
            path[path.size() - 1 - i]
        ) {
            return false;
        }
    }

    return true;
}

static std::string path_string(
    const std::vector<int>& path
) {
    if (path.empty()) {
        return "";
    }

    std::string result =
        std::to_string(path[0]);

    for (
        std::size_t i = 1;
        i < path.size();
        ++i
    ) {
        result += "->";
        result +=
            std::to_string(path[i]);
    }

    return result;
}

static std::string permutation_string(
    const Permutation& permutation
) {
    const int values[5] =
        {1, 2, 3, 4, 6};

    std::string result;

    for (int i = 0; i < 5; ++i) {
        if (i != 0) {
            result += ",";
        }

        result +=
            std::to_string(values[i]);

        result += "->";

        result +=
            std::to_string(
                permutation.map[i]
            );
    }

    return result;
}

static bool next_permutation_manual(
    Permutation& permutation
) {
    int i = 3;

    while (
        i >= 0 &&
        permutation.map[i] >=
        permutation.map[i + 1]
    ) {
        --i;
    }

    if (i < 0) {
        return false;
    }

    int j = 4;

    while (
        permutation.map[j] <=
        permutation.map[i]
    ) {
        --j;
    }

    const int tmp =
        permutation.map[i];

    permutation.map[i] =
        permutation.map[j];

    permutation.map[j] =
        tmp;

    int left = i + 1;
    int right = 4;

    while (left < right) {
        const int x =
            permutation.map[left];

        permutation.map[left] =
            permutation.map[right];

        permutation.map[right] = x;

        ++left;
        --right;
    }

    return true;
}

int main() {
    constexpr int EXPERIMENT = 447;
    constexpr int PRIME_LIMIT = 10000;
    constexpr int K_LIMIT = 3000;

    std::cout
        << "START EXPERIMENT "
        << EXPERIMENT
        << "\n";

    std::cout
        << "PRIME_LIMIT="
        << PRIME_LIMIT
        << "\n";

    std::cout
        << "K_LIMIT="
        << K_LIMIT
        << "\n";

    const std::vector<u64> primes =
        generate_primes(
            PRIME_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    const int base[5] =
        {1, 2, 3, 4, 6};

    std::map<std::string, u64>
        exact_path_histogram;

    std::map<std::string, u64>
        path_length_histogram;

    u64 total_pairs = 0;

    u64 ordinary_palindromes = 0;

    u64 involution_matches_total = 0;

    u64 best_involution_frequency[120] = {};

    u64 involution_matching_paths = 0;

    u64 best_match_sum = 0;

    u64 max_matches_for_any_permutation = 0;

    std::string best_permutation;

    std::map<std::string, u64>
        involution_path_histogram;

    std::vector<std::pair<
        u64,
        std::string
    >> path_examples;

    /*
     * Evaluate every exceptional reciprocal pair.
     */
    for (u64 r : primes) {
        for (int source_sign : {-1, +1}) {
            const std::vector<PairInfo> pairs =
                find_exceptional_pairs(
                    r,
                    source_sign,
                    K_LIMIT
                );

            if (pairs.empty()) {
                continue;
            }

            const std::vector<Witness> envelope =
                build_envelope(
                    r,
                    K_LIMIT
                );

            for (
                const PairInfo& pair :
                pairs
            ) {
                ++total_pairs;

                std::vector<int> raw_path;

                for (
                    u64 K = pair.k;
                    K <= pair.t;
                    ++K
                ) {
                    raw_path.push_back(
                        envelope[
                            static_cast<std::size_t>(K)
                        ].m
                    );
                }

                const std::vector<int> path =
                    compress_path(
                        raw_path
                    );

                const std::string path_text =
                    path_string(path);

                ++exact_path_histogram[
                    path_text
                ];

                ++path_length_histogram[
                    std::to_string(path.size())
                ];

                if (
                    path_matches_identity(
                        path
                    )
                ) {
                    ++ordinary_palindromes;
                }

                /*
                 * Search all 120 permutations.
                 */
                Permutation permutation;

                for (int i = 0; i < 5; ++i) {
                    permutation.map[i] =
                        base[i];
                }

                u64 local_best_matches = 0;
                std::vector<std::string>
                    local_matching_permutations;

                int permutation_index = 0;

                do {
                    u64 matches = 0;

                    bool matches_path =
                        path_matches_permutation(
                            path,
                            permutation
                        );

                    if (matches_path) {
                        matches = 1;
                    }

                    if (
                        matches >
                        local_best_matches
                    ) {
                        local_best_matches =
                            matches;

                        local_matching_permutations.clear();

                        local_matching_permutations.push_back(
                            permutation_string(
                                permutation
                            )
                        );
                    } else if (
                        matches ==
                        local_best_matches &&
                        matches == 1
                    ) {
                        local_matching_permutations.push_back(
                            permutation_string(
                                permutation
                            )
                        );
                    }

                    if (
                        matches == 1 &&
                        permutation_is_involution(
                            permutation
                        )
                    ) {
                        ++involution_matches_total;

                        ++involution_path_histogram[
                            path_text
                        ];

                        involution_matching_paths++;

                        std::cout
                            << "\nINVOLUTION_MATCH\n";

                        std::cout
                            << "R="
                            << pair.r
                            << " SIGN="
                            << (
                                pair.source_sign > 0
                                    ? "+1"
                                    : "-1"
                            )
                            << " A="
                            << pair.A
                            << " K="
                            << pair.k
                            << " T="
                            << pair.t
                            << "\n";

                        std::cout
                            << "PATH="
                            << path_text
                            << "\n";

                        std::cout
                            << "PERMUTATION="
                            << permutation_string(
                                permutation
                            )
                            << "\n";
                    }

                    ++permutation_index;

                } while (
                    next_permutation_manual(
                        permutation
                    )
                );

                if (
                    local_best_matches >
                    max_matches_for_any_permutation
                ) {
                    max_matches_for_any_permutation =
                        local_best_matches;

                    if (
                        !local_matching_permutations.empty()
                    ) {
                        best_permutation =
                            local_matching_permutations[0];
                    }
                }

                best_match_sum +=
                    local_best_matches;

                path_examples.push_back({
                    pair.r,
                    path_text
                });
            }
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_PAIRS="
        << total_pairs
        << "\n";

    std::cout
        << "ORDINARY_PALINDROMES="
        << ordinary_palindromes
        << "\n";

    std::cout
        << "INVOLUTION_MATCHES_TOTAL="
        << involution_matches_total
        << "\n";

    std::cout
        << "PATHS_WITH_INVOLUTION_MATCH="
        << involution_matching_paths
        << "\n";

    std::cout
        << "MAX_MATCHES_FOR_ANY_PERMUTATION="
        << max_matches_for_any_permutation
        << "\n";

    std::cout
        << "BEST_PERMUTATION="
        << (
            best_permutation.empty()
                ? "NONE"
                : best_permutation
        )
        << "\n";

    std::cout
        << "BEST_MATCH_SUM="
        << best_match_sum
        << "\n";

    std::cout
        << "\nEXACT_PATH_HISTOGRAM\n";

    for (
        const auto& entry :
        exact_path_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nPATH_LENGTH_HISTOGRAM\n";

    for (
        const auto& entry :
        path_length_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nINVOLUTION_PATH_HISTOGRAM\n";

    for (
        const auto& entry :
        involution_path_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
