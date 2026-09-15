#include <cstdint>
#include <iostream>
#include <limits>
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

static std::vector<u64> generate_primes(int limit) {
    std::vector<bool> is_prime(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    is_prime[0] = false;
    is_prime[1] = false;

    for (
        int p = 2;
        static_cast<long long>(p) * p <= limit;
        ++p
    ) {
        if (!is_prime[p]) {
            continue;
        }

        for (
            int x = p * p;
            x <= limit;
            x += p
        ) {
            is_prime[x] = false;
        }
    }

    std::vector<u64> primes;

    for (int x = 2; x <= limit; ++x) {
        if (is_prime[x]) {
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

    if (value <= 0) {
        return false;
    }

    return
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

static bool allowed_m(int m) {
    return
        m == 1 ||
        m == 2 ||
        m == 3 ||
        m == 4 ||
        m == 6;
}

static Witness global_winner(
    const std::vector<Witness>& best_by_m
) {
    Witness best = best_by_m[1];

    for (
        std::size_t m = 2;
        m < best_by_m.size();
        ++m
    ) {
        if (
            better_normalized(
                best_by_m[m],
                best
            )
        ) {
            best = best_by_m[m];
        }
    }

    return best;
}

static Witness allowed_winner(
    const std::vector<Witness>& best_by_m
) {
    Witness best = best_by_m[1];

    for (int m : {2, 3, 4, 6}) {
        const Witness& candidate =
            best_by_m[
                static_cast<std::size_t>(m)
            ];

        if (
            better_normalized(
                candidate,
                best
            )
        ) {
            best = candidate;
        }
    }

    return best;
}

static void print_witness(
    u64 r,
    const std::string& label,
    const Witness& w
) {
    const i128 value =
        static_cast<i128>(w.m) *
        static_cast<i128>(r) +
        static_cast<i128>(w.sign);

    std::cout
        << label
        << "=(m=" << w.m
        << ",k=" << w.k
        << ",t=" << w.t
        << ",sign="
        << (w.sign > 0 ? "+1" : "-1")
        << ",value="
        << static_cast<u64>(value)
        << ")";
}

int main() {
    constexpr int EXPERIMENT = 437;

    constexpr int PRIME_LIMIT = 10000;
    constexpr int K_LIMIT = 3000;
    constexpr int M_MAX = 64;

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

    std::cout
        << "M_MAX="
        << M_MAX
        << "\n";

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    u64 total_points = 0;

    u64 allowed_winners = 0;
    u64 forbidden_winners = 0;

    u64 forbidden_track_records = 0;
    u64 forbidden_track_candidates = 0;

    u64 forbidden_beats_m1 = 0;
    u64 forbidden_beats_allowed = 0;

    u64 allowed_matches_global = 0;
    u64 allowed_mismatch = 0;

    int max_global_m = 0;
    int max_forbidden_track_m = 0;

    std::map<int, u64>
        winner_m_histogram;

    std::map<int, u64>
        forbidden_winner_histogram;

    std::map<int, u64>
        forbidden_record_histogram;

    std::vector<Witness>
        first_forbidden_winner;

    std::vector<Witness>
        first_forbidden_beats_allowed;

    u64 first_forbidden_prime = 0;
    u64 first_forbidden_beats_allowed_prime = 0;

    u64 first_forbidden_K = 0;
    u64 first_forbidden_beats_allowed_K = 0;

    bool have_first_forbidden = false;
    bool have_first_forbidden_beats_allowed = false;

    for (u64 r : primes) {
        std::vector<Witness> best_by_m(
            static_cast<std::size_t>(M_MAX) + 1
        );

        for (int m = 1; m <= M_MAX; ++m) {
            best_by_m[
                static_cast<std::size_t>(m)
            ] = make_witness(
                r,
                m,
                1,
                +1
            );
        }

        /*
         * K = 1.
         */
        {
            ++total_points;

            const Witness global =
                global_winner(best_by_m);

            winner_m_histogram[
                global.m
            ]++;

            if (allowed_m(global.m)) {
                ++allowed_winners;
                ++allowed_matches_global;
            } else {
                ++forbidden_winners;

                forbidden_winner_histogram[
                    global.m
                ]++;

                if (!have_first_forbidden) {
                    first_forbidden_winner.clear();
                    first_forbidden_winner.push_back(
                        global
                    );

                    first_forbidden_prime = r;
                    first_forbidden_K = 1;

                    have_first_forbidden = true;
                }
            }
        }

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            /*
             * Update every m-track.
             */
            for (int m = 1; m <= M_MAX; ++m) {
                const std::size_t mi =
                    static_cast<std::size_t>(m);

                bool found_candidate = false;
                Witness best_candidate;

                if (
                    divides_stream(
                        r,
                        m,
                        K,
                        +1
                    )
                ) {
                    best_candidate =
                        make_witness(
                            r,
                            m,
                            K,
                            +1
                        );

                    found_candidate = true;

                    ++forbidden_track_candidates;
                }

                if (
                    divides_stream(
                        r,
                        m,
                        K,
                        -1
                    )
                ) {
                    Witness minus_candidate =
                        make_witness(
                            r,
                            m,
                            K,
                            -1
                        );

                    if (
                        !found_candidate ||
                        better_normalized(
                            minus_candidate,
                            best_candidate
                        )
                    ) {
                        best_candidate =
                            minus_candidate;
                    }

                    found_candidate = true;

                    ++forbidden_track_candidates;
                }

                if (!found_candidate) {
                    continue;
                }

                if (
                    best_candidate.k >
                    best_by_m[mi].k
                ) {
                    best_by_m[mi] =
                        best_candidate;

                    /*
                     * Only forbidden m are counted here.
                     */
                    if (!allowed_m(m)) {
                        ++forbidden_track_records;

                        forbidden_record_histogram[
                            m
                        ]++;

                        if (
                            m >
                            max_forbidden_track_m
                        ) {
                            max_forbidden_track_m =
                                m;
                        }
                    }
                }
            }

            ++total_points;

            const Witness global =
                global_winner(best_by_m);

            const Witness allowed =
                allowed_winner(best_by_m);

            winner_m_histogram[
                global.m
            ]++;

            if (global.m > max_global_m) {
                max_global_m = global.m;
            }

            if (allowed.m == global.m) {
                ++allowed_matches_global;
            } else {
                ++allowed_mismatch;
            }

            if (allowed_m(global.m)) {
                ++allowed_winners;
            } else {
                ++forbidden_winners;

                forbidden_winner_histogram[
                    global.m
                ]++;

                if (!have_first_forbidden) {
                    first_forbidden_winner.clear();
                    first_forbidden_winner.push_back(
                        global
                    );

                    first_forbidden_prime = r;
                    first_forbidden_K = K;

                    have_first_forbidden = true;
                }
            }

            /*
             * Compare the strongest forbidden track
             * against the strongest allowed track.
             */
            Witness strongest_forbidden =
                best_by_m[1];

            bool have_forbidden = false;

            for (
                int m = 5;
                m <= M_MAX;
                ++m
            ) {
                if (allowed_m(m)) {
                    continue;
                }

                const Witness& candidate =
                    best_by_m[
                        static_cast<std::size_t>(m)
                    ];

                if (
                    !have_forbidden ||
                    better_normalized(
                        candidate,
                        strongest_forbidden
                    )
                ) {
                    strongest_forbidden =
                        candidate;

                    have_forbidden = true;
                }
            }

            if (have_forbidden) {
                if (
                    better_normalized(
                        strongest_forbidden,
                        best_by_m[1]
                    )
                ) {
                    ++forbidden_beats_m1;
                }

                if (
                    better_normalized(
                        strongest_forbidden,
                        allowed
                    )
                ) {
                    ++forbidden_beats_allowed;

                    if (
                        !have_first_forbidden_beats_allowed
                    ) {
                        first_forbidden_beats_allowed.clear();
                        first_forbidden_beats_allowed.push_back(
                            allowed
                        );
                        first_forbidden_beats_allowed.push_back(
                            strongest_forbidden
                        );

                        first_forbidden_beats_allowed_prime =
                            r;

                        first_forbidden_beats_allowed_K =
                            K;

                        have_first_forbidden_beats_allowed =
                            true;
                    }
                }
            }
        }
    }

    std::cout
        << "\nTOTAL_POINTS="
        << total_points
        << "\n";

    std::cout
        << "ALLOWED_WINNERS="
        << allowed_winners
        << "\n";

    std::cout
        << "FORBIDDEN_WINNERS="
        << forbidden_winners
        << "\n";

    std::cout
        << "ALLOWED_MATCHES_GLOBAL="
        << allowed_matches_global
        << "\n";

    std::cout
        << "ALLOWED_MISMATCH="
        << allowed_mismatch
        << "\n";

    std::cout
        << "\nGLOBAL_M\n";

    std::cout
        << "MAX_GLOBAL_M="
        << max_global_m
        << "\n";

    std::cout
        << "MAX_FORBIDDEN_TRACK_M="
        << max_forbidden_track_m
        << "\n";

    std::cout
        << "FORBIDDEN_TRACK_CANDIDATES="
        << forbidden_track_candidates
        << "\n";

    std::cout
        << "FORBIDDEN_TRACK_RECORDS="
        << forbidden_track_records
        << "\n";

    std::cout
        << "FORBIDDEN_BEATS_M1="
        << forbidden_beats_m1
        << "\n";

    std::cout
        << "FORBIDDEN_BEATS_ALLOWED="
        << forbidden_beats_allowed
        << "\n";

    std::cout
        << "\nWINNER_M_HISTOGRAM\n";

    for (
        const auto& entry :
        winner_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFORBIDDEN_WINNER_HISTOGRAM\n";

    for (
        const auto& entry :
        forbidden_winner_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFORBIDDEN_RECORD_M_HISTOGRAM\n";

    for (
        const auto& entry :
        forbidden_record_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_forbidden) {
        std::cout
            << "\nFIRST_FORBIDDEN_GLOBAL_WINNER\n";

        std::cout
            << "PRIME="
            << first_forbidden_prime
            << "\n";

        std::cout
            << "K="
            << first_forbidden_K
            << "\n";

        print_witness(
            first_forbidden_prime,
            "WINNER",
            first_forbidden_winner.front()
        );

        std::cout
            << "\n";
    } else {
        std::cout
            << "\nFIRST_FORBIDDEN_GLOBAL_WINNER NONE\n";
    }

    if (have_first_forbidden_beats_allowed) {
        std::cout
            << "\nFIRST_FORBIDDEN_BEATS_ALLOWED\n";

        std::cout
            << "PRIME="
            << first_forbidden_beats_allowed_prime
            << "\n";

        std::cout
            << "K="
            << first_forbidden_beats_allowed_K
            << "\n";

        print_witness(
            first_forbidden_beats_allowed_prime,
            "ALLOWED",
            first_forbidden_beats_allowed.front()
        );

        std::cout << "\n";

        print_witness(
            first_forbidden_beats_allowed_prime,
            "FORBIDDEN",
            first_forbidden_beats_allowed.back()
        );

        std::cout << "\n";
    } else {
        std::cout
            << "\nFIRST_FORBIDDEN_BEATS_ALLOWED NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
