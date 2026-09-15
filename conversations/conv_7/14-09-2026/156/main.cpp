#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Witness {
    int m = 1;
    u64 k = 1;
    u64 t = 0;
    int sign = +1;
};

struct PairInfo {
    u64 r = 0;
    int source_sign = +1;
    u64 A = 0;
    u64 k = 0;
    u64 t = 0;
};

static std::vector<u64> generate_primes(int limit) {
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
    const i128 value =
        static_cast<i128>(m) *
        static_cast<i128>(r) +
        static_cast<i128>(sign);

    return
        value > 0 &&
        value % static_cast<i128>(k) == 0;
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

    const i128 value =
        static_cast<i128>(m) *
        static_cast<i128>(r) +
        static_cast<i128>(sign);

    w.t = static_cast<u64>(
        value /
        static_cast<i128>(k)
    );

    return w;
}

static bool better_normalized(
    const Witness& a,
    const Witness& b
) {
    const i128 lhs =
        static_cast<i128>(a.k) *
        static_cast<i128>(b.m);

    const i128 rhs =
        static_cast<i128>(b.k) *
        static_cast<i128>(a.m);

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
    const i128 value =
        static_cast<i128>(k) *
        static_cast<i128>(allowed.m) -
        static_cast<i128>(5) *
        static_cast<i128>(allowed.k);

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

static std::vector<PairInfo> find_exceptional_pairs(
    u64 r,
    int source_sign,
    int K_LIMIT
) {
    const i128 value =
        static_cast<i128>(5) *
        static_cast<i128>(r) +
        static_cast<i128>(source_sign);

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

        if (
            margin_for(
                k,
                envelope[
                    static_cast<std::size_t>(k)
                ]
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

static std::vector<int> compressed_path(
    const std::vector<Witness>& envelope,
    u64 begin,
    u64 end
) {
    std::vector<int> path;

    for (
        u64 K = begin;
        K <= end;
        ++K
    ) {
        const int m =
            envelope[
                static_cast<std::size_t>(K)
            ].m;

        if (
            path.empty() ||
            path.back() != m
        ) {
            path.push_back(m);
        }
    }

    return path;
}

static int endpoint_map(
    int left,
    int right,
    int x
) {
    if (x == left) {
        return right;
    }

    if (x == right) {
        return left;
    }

    return x;
}

static bool endpoint_path_symmetry(
    const std::vector<int>& path
) {
    if (path.empty()) {
        return true;
    }

    const int left =
        path.front();

    const int right =
        path.back();

    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        const int mirrored =
            path[
                path.size() - 1 - i
            ];

        const int expected =
            endpoint_map(
                left,
                right,
                mirrored
            );

        if (path[i] != expected) {
            return false;
        }
    }

    return true;
}

static std::size_t first_symmetry_failure(
    const std::vector<int>& path
) {
    if (path.empty()) {
        return path.size();
    }

    const int left =
        path.front();

    const int right =
        path.back();

    for (
        std::size_t i = 0;
        i < path.size();
        ++i
    ) {
        const int mirrored =
            path[
                path.size() - 1 - i
            ];

        const int expected =
            endpoint_map(
                left,
                right,
                mirrored
            );

        if (path[i] != expected) {
            return i;
        }
    }

    return path.size();
}

static std::string path_to_string(
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

static std::string mapping_to_string(
    int left,
    int right
) {
    const int values[5] =
        {1, 2, 3, 4, 6};

    std::string result;

    for (int i = 0; i < 5; ++i) {
        if (i != 0) {
            result += ",";
        }

        const int m = values[i];

        result +=
            std::to_string(m);

        result += "->";

        result +=
            std::to_string(
                endpoint_map(
                    left,
                    right,
                    m
                )
            );
    }

    return result;
}

int main() {
    constexpr int EXPERIMENT = 448;
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

    u64 total_pairs = 0;
    u64 symmetry_success = 0;
    u64 symmetry_failure = 0;

    u64 length1_success = 0;
    u64 length2_success = 0;
    u64 length3_success = 0;
    u64 length4plus_success = 0;

    std::map<std::string, u64>
        endpoint_mapping_histogram;

    std::map<std::string, u64>
        successful_path_histogram;

    std::map<std::string, u64>
        failed_path_histogram;

    std::map<int, u64>
        first_failure_position_histogram;

    std::vector<std::string>
        failure_examples;

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

                const std::vector<int> path =
                    compressed_path(
                        envelope,
                        pair.k,
                        pair.t
                    );

                const std::string path_text =
                    path_to_string(path);

                const std::string mapping =
                    mapping_to_string(
                        path.front(),
                        path.back()
                    );

                ++endpoint_mapping_histogram[
                    mapping
                ];

                const bool success =
                    endpoint_path_symmetry(
                        path
                    );

                std::cout
                    << "\nPAIR\n";

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
                    << "ENDPOINT_MAPPING="
                    << mapping
                    << "\n";

                std::cout
                    << "ENDPOINT_SYMMETRY="
                    << (
                        success
                            ? "PASS"
                            : "FAIL"
                    )
                    << "\n";

                if (success) {
                    ++symmetry_success;

                    ++successful_path_histogram[
                        path_text
                    ];

                    if (path.size() == 1) {
                        ++length1_success;
                    } else if (path.size() == 2) {
                        ++length2_success;
                    } else if (path.size() == 3) {
                        ++length3_success;
                    } else {
                        ++length4plus_success;
                    }
                } else {
                    ++symmetry_failure;

                    ++failed_path_histogram[
                        path_text
                    ];

                    const std::size_t failure =
                        first_symmetry_failure(
                            path
                        );

                    ++first_failure_position_histogram[
                        static_cast<int>(failure)
                    ];

                    const int actual =
                        path[failure];

                    const int mirrored =
                        path[
                            path.size() -
                            1 -
                            failure
                        ];

                    const int expected =
                        endpoint_map(
                            path.front(),
                            path.back(),
                            mirrored
                        );

                    std::cout
                        << "FIRST_FAILURE_INDEX="
                        << failure
                        << "\n";

                    std::cout
                        << "FIRST_FAILURE_ACTUAL="
                        << actual
                        << "\n";

                    std::cout
                        << "FIRST_FAILURE_MIRROR="
                        << mirrored
                        << "\n";

                    std::cout
                        << "FIRST_FAILURE_EXPECTED="
                        << expected
                        << "\n";

                    if (
                        failure_examples.size() < 20
                    ) {
                        std::string example;

                        example +=
                            "R=" +
                            std::to_string(
                                pair.r
                            );

                        example +=
                            " SIGN=" +
                            std::string(
                                source_sign > 0
                                    ? "+1"
                                    : "-1"
                            );

                        example +=
                            " K=" +
                            std::to_string(
                                pair.k
                            );

                        example +=
                            " T=" +
                            std::to_string(
                                pair.t
                            );

                        example +=
                            " PATH=" +
                            path_text;

                        example +=
                            " MAPPING=" +
                            mapping;

                        example +=
                            " FAILURE_INDEX=" +
                            std::to_string(
                                failure
                            );

                        example +=
                            " ACTUAL=" +
                            std::to_string(
                                actual
                            );

                        example +=
                            " MIRROR=" +
                            std::to_string(
                                mirrored
                            );

                        example +=
                            " EXPECTED=" +
                            std::to_string(
                                expected
                            );

                        failure_examples.push_back(
                            example
                        );
                    }
                }
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
        << "ENDPOINT_SYMMETRY_SUCCESS="
        << symmetry_success
        << "\n";

    std::cout
        << "ENDPOINT_SYMMETRY_FAILURE="
        << symmetry_failure
        << "\n";

    std::cout
        << "LENGTH1_SUCCESS="
        << length1_success
        << "\n";

    std::cout
        << "LENGTH2_SUCCESS="
        << length2_success
        << "\n";

    std::cout
        << "LENGTH3_SUCCESS="
        << length3_success
        << "\n";

    std::cout
        << "LENGTH4PLUS_SUCCESS="
        << length4plus_success
        << "\n";

    std::cout
        << "\nENDPOINT_MAPPING_HISTOGRAM\n";

    for (
        const auto& entry :
        endpoint_mapping_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nSUCCESSFUL_PATH_HISTOGRAM\n";

    for (
        const auto& entry :
        successful_path_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFAILED_PATH_HISTOGRAM\n";

    for (
        const auto& entry :
        failed_path_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFIRST_FAILURE_POSITION_HISTOGRAM\n";

    for (
        const auto& entry :
        first_failure_position_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFAILURE_EXAMPLES\n";

    for (
        const std::string& example :
        failure_examples
    ) {
        std::cout
            << example
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}