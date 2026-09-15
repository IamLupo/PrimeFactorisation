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

    Witness left_allowed;
    Witness right_allowed;

    u64 left_margin = 0;
    u64 right_margin = 0;
};

static std::vector<u64> generate_primes(int limit) {
    std::vector<bool> sieve(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    if (limit >= 0) {
        sieve[0] = false;
    }

    if (limit >= 1) {
        sieve[1] = false;
    }

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

static bool same_envelope_state(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.k == b.k &&
        a.sign == b.sign;
}

static void print_witness(
    const std::string& label,
    const Witness& w
) {
    std::cout
        << label
        << "(m="
        << w.m
        << ",k="
        << w.k
        << ",t="
        << w.t
        << ",sign="
        << (
            w.sign > 0
                ? "+1"
                : "-1"
        )
        << ")";
}

static void print_segment(
    const Segment& segment
) {
    std::cout
        << "["
        << segment.begin
        << ","
        << segment.end
        << "] ";

    print_witness(
        "ENVELOPE",
        segment.witness
    );

    std::cout
        << "\n";
}

static std::vector<Segment> compress_segments(
    const std::vector<Witness>& envelope,
    u64 begin,
    u64 end
) {
    std::vector<Segment> result;

    if (begin > end) {
        return result;
    }

    Segment current;

    current.begin = begin;
    current.end = begin;
    current.witness =
        envelope[
            static_cast<std::size_t>(begin)
        ];

    for (
        u64 K = begin + 1;
        K <= end;
        ++K
    ) {
        const Witness& w =
            envelope[
                static_cast<std::size_t>(K)
            ];

        if (
            same_envelope_state(
                w,
                current.witness
            )
        ) {
            current.end = K;
        } else {
            result.push_back(current);

            current.begin = K;
            current.end = K;
            current.witness = w;
        }
    }

    result.push_back(current);

    return result;
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

        pair.left_allowed =
            envelope[
                static_cast<std::size_t>(k)
            ];

        pair.right_allowed =
            envelope[
                static_cast<std::size_t>(t)
            ];

        pair.left_margin =
            margin_for(
                k,
                pair.left_allowed
            );

        pair.right_margin =
            margin_for(
                t,
                pair.right_allowed
            );

        pairs.push_back(pair);
    }

    return pairs;
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

int main() {
    constexpr int EXPERIMENT = 446;
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

    u64 one_segment = 0;
    u64 two_segments = 0;
    u64 three_segments = 0;
    u64 more_segments = 0;

    u64 zero_multiplier_changes = 0;
    u64 exactly_one_multiplier_change = 0;
    u64 multiple_multiplier_changes = 0;

    u64 same_endpoint_multiplier = 0;
    u64 different_endpoint_multiplier = 0;

    u64 endpoint_same_segment = 0;
    u64 endpoint_different_segment = 0;

    std::map<std::string, u64>
        path_histogram;

    std::map<int, u64>
        first_transition_histogram;

    std::map<int, u64>
        last_transition_histogram;

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

                const std::vector<Segment> segments =
                    compress_segments(
                        envelope,
                        pair.k,
                        pair.t
                    );

                const std::size_t count =
                    segments.size();

                if (count == 1) {
                    ++one_segment;
                } else if (count == 2) {
                    ++two_segments;
                } else if (count == 3) {
                    ++three_segments;
                } else {
                    ++more_segments;
                }

                if (
                    pair.left_allowed.m ==
                    pair.right_allowed.m
                ) {
                    ++same_endpoint_multiplier;
                } else {
                    ++different_endpoint_multiplier;
                }

                if (
                    segments.front().witness.m ==
                    segments.back().witness.m
                ) {
                    ++endpoint_same_segment;
                } else {
                    ++endpoint_different_segment;
                }

                /*
                 * Compress the segment sequence to multiplier
                 * transitions only.
                 */
                std::vector<int> path;

                for (
                    const Segment& segment :
                    segments
                ) {
                    const int m =
                        segment.witness.m;

                    if (
                        path.empty() ||
                        path.back() != m
                    ) {
                        path.push_back(m);
                    }
                }

                const std::size_t
                    multiplier_changes =
                        path.empty()
                            ? 0
                            : path.size() - 1;

                if (
                    multiplier_changes == 0
                ) {
                    ++zero_multiplier_changes;
                } else if (
                    multiplier_changes == 1
                ) {
                    ++exactly_one_multiplier_change;
                } else {
                    ++multiple_multiplier_changes;
                }

                const std::string path_string =
                    path_to_string(path);

                ++path_histogram[
                    path_string
                ];

                if (path.size() >= 2) {
                    ++first_transition_histogram[
                        path[1]
                    ];

                    ++last_transition_histogram[
                        path[
                            path.size() - 2
                        ]
                    ];
                } else if (!path.empty()) {
                    ++first_transition_histogram[
                        path[0]
                    ];

                    ++last_transition_histogram[
                        path[0]
                    ];
                }

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

                print_witness(
                    "LEFT_ALLOWED",
                    pair.left_allowed
                );

                std::cout
                    << " MARGIN="
                    << pair.left_margin
                    << "\n";

                print_witness(
                    "RIGHT_ALLOWED",
                    pair.right_allowed
                );

                std::cout
                    << " MARGIN="
                    << pair.right_margin
                    << "\n";

                std::cout
                    << "SEGMENT_COUNT="
                    << count
                    << "\n";

                std::cout
                    << "MULTIPLIER_CHANGE_COUNT="
                    << multiplier_changes
                    << "\n";

                std::cout
                    << "PATH="
                    << path_string
                    << "\n";

                std::cout
                    << "ENVELOPE_SEGMENTS\n";

                for (
                    const Segment& segment :
                    segments
                ) {
                    print_segment(
                        segment
                    );
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
        << "ONE_SEGMENT="
        << one_segment
        << "\n";

    std::cout
        << "TWO_SEGMENTS="
        << two_segments
        << "\n";

    std::cout
        << "THREE_SEGMENTS="
        << three_segments
        << "\n";

    std::cout
        << "MORE_SEGMENTS="
        << more_segments
        << "\n";

    std::cout
        << "ZERO_MULTIPLIER_CHANGES="
        << zero_multiplier_changes
        << "\n";

    std::cout
        << "EXACTLY_ONE_MULTIPLIER_CHANGE="
        << exactly_one_multiplier_change
        << "\n";

    std::cout
        << "MULTIPLE_MULTIPLIER_CHANGES="
        << multiple_multiplier_changes
        << "\n";

    std::cout
        << "SAME_ENDPOINT_MULTIPLIER="
        << same_endpoint_multiplier
        << "\n";

    std::cout
        << "DIFFERENT_ENDPOINT_MULTIPLIER="
        << different_endpoint_multiplier
        << "\n";

    std::cout
        << "ENDPOINT_SAME_SEGMENT="
        << endpoint_same_segment
        << "\n";

    std::cout
        << "ENDPOINT_DIFFERENT_SEGMENT="
        << endpoint_different_segment
        << "\n";

    std::cout
        << "\nPATH_HISTOGRAM\n";

    for (
        const auto& entry :
        path_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFIRST_TRANSITION_HISTOGRAM\n";

    for (
        const auto& entry :
        first_transition_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nLAST_TRANSITION_HISTOGRAM\n";

    for (
        const auto& entry :
        last_transition_histogram
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