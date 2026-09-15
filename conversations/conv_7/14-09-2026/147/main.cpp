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

struct ExceptionalPoint {
    u64 r = 0;
    u64 K = 0;
    Witness winner;
    Witness allowed;
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

static i128 arithmetic_value(
    u64 r,
    const Witness& w
) {
    return
        static_cast<i128>(w.m) *
        static_cast<i128>(r) +
        static_cast<i128>(w.sign);
}

static bool better_than(
    const Witness& a,
    const Witness& b
) {
    return better_normalized(a, b);
}

static u64 exact_threshold(
    u64 r,
    const Witness& a,
    const Witness& b
) {
    /*
     * Compare:
     *
     * t_a < t_b
     *
     * where
     *
     * t_a=(m_a r+e_a)/k_a
     * t_b=(m_b r+e_b)/k_b.
     *
     * Assuming a is normalized-first:
     *
     * Delta = m_b*k_a - m_a*k_b > 0
     * C     = e_a*k_b - e_b*k_a
     *
     * r_min=floor(C/Delta)+1.
     */

    const i128 delta =
        static_cast<i128>(b.m) *
        static_cast<i128>(a.k) -
        static_cast<i128>(a.m) *
        static_cast<i128>(b.k);

    const i128 c =
        static_cast<i128>(a.sign) *
        static_cast<i128>(b.k) -
        static_cast<i128>(b.sign) *
        static_cast<i128>(a.k);

    if (delta <= 0) {
        return 0;
    }

    if (c < 0) {
        return 0;
    }

    const i128 q =
        c / delta;

    return static_cast<u64>(q + 1);
}

static u64 score_cross_margin(
    const Witness& a,
    const Witness& b
) {
    /*
     * Positive means a has the larger normalized score.
     *
     * a.k/a.m - b.k/b.m
     *
     * represented without division.
     */
    const i128 value =
        static_cast<i128>(a.k) *
        static_cast<i128>(b.m) -
        static_cast<i128>(b.k) *
        static_cast<i128>(a.m);

    if (value <= 0) {
        return 0;
    }

    return static_cast<u64>(value);
}

static void print_witness(
    u64 r,
    const std::string& label,
    const Witness& w
) {
    std::cout
        << label
        << "=(m=" << w.m
        << ",k=" << w.k
        << ",t=" << w.t
        << ",sign="
        << (w.sign > 0 ? "+1" : "-1")
        << ",value="
        << static_cast<u64>(
            arithmetic_value(r, w)
        )
        << ")";
}

int main() {
    constexpr int EXPERIMENT = 439;
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

    u64 exceptional_points = 0;

    u64 m5_points = 0;
    u64 m7_points = 0;

    u64 m5_beat_m1 = 0;
    u64 m5_beat_m2 = 0;
    u64 m5_beat_m3 = 0;
    u64 m5_beat_m4 = 0;
    u64 m5_beat_m6 = 0;

    u64 m7_beat_m1 = 0;
    u64 m7_beat_m2 = 0;
    u64 m7_beat_m3 = 0;
    u64 m7_beat_m4 = 0;
    u64 m7_beat_m6 = 0;

    u64 m5_all_allowed_beaten = 0;
    u64 m7_all_allowed_beaten = 0;

    u64 m5_new_record_points = 0;
    u64 m7_new_record_points = 0;

    u64 m5_same_record_as_previous = 0;
    u64 m7_same_record_as_previous = 0;

    u64 m5_threshold_consistent = 0;
    u64 m7_threshold_consistent = 0;

    u64 m5_sharp_boundary = 0;
    u64 m7_sharp_boundary = 0;

    u64 m5_threshold_max = 0;
    u64 m7_threshold_max = 0;

    u64 m5_margin_min = 0;
    u64 m7_margin_min = 0;

    u64 m5_margin_max = 0;
    u64 m7_margin_max = 0;

    std::map<u64, u64>
        m5_threshold_histogram;

    std::map<u64, u64>
        m7_threshold_histogram;

    std::map<u64, u64>
        m5_score_margin_histogram;

    std::map<u64, u64>
        m7_score_margin_histogram;

    std::map<int, u64>
        m5_competitor_histogram;

    std::map<int, u64>
        m7_competitor_histogram;

    ExceptionalPoint first_m5;
    ExceptionalPoint first_m7;

    bool have_first_m5 = false;
    bool have_first_m7 = false;

    u64 previous_m5_k = 0;
    u64 previous_m7_k = 0;

    /*
     * These track the per-prime plateau state.
     */
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

        u64 last_m5_k = 0;
        u64 last_m7_k = 0;

        for (
            u64 K = 1;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            if (K > 1) {
                for (int m = 1; m <= M_MAX; ++m) {
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
            }

            const Witness global =
                global_winner(best_by_m);

            if (
                global.m != 5 &&
                global.m != 7
            ) {
                continue;
            }

            ++exceptional_points;

            const Witness allowed =
                allowed_winner(best_by_m);

            if (global.m == 5) {
                ++m5_points;

                if (!have_first_m5) {
                    first_m5.r = r;
                    first_m5.K = K;
                    first_m5.winner = global;
                    first_m5.allowed = allowed;
                    have_first_m5 = true;
                }

                /*
                 * Check whether k changed at this point.
                 */
                if (global.k != last_m5_k) {
                    ++m5_new_record_points;
                    last_m5_k = global.k;
                } else {
                    ++m5_same_record_as_previous;
                }

                const Witness& w1 =
                    best_by_m[1];

                const Witness& w2 =
                    best_by_m[2];

                const Witness& w3 =
                    best_by_m[3];

                const Witness& w4 =
                    best_by_m[4];

                const Witness& w6 =
                    best_by_m[6];

                if (better_than(global, w1)) {
                    ++m5_beat_m1;
                }

                if (better_than(global, w2)) {
                    ++m5_beat_m2;
                }

                if (better_than(global, w3)) {
                    ++m5_beat_m3;
                }

                if (better_than(global, w4)) {
                    ++m5_beat_m4;
                }

                if (better_than(global, w6)) {
                    ++m5_beat_m6;
                }

                /*
                 * Count strongest allowed competitor.
                 */
                Witness strongest = w1;
                int strongest_m = 1;

                for (int m : {2, 3, 4, 6}) {
                    const Witness& candidate =
                        best_by_m[
                            static_cast<std::size_t>(m)
                        ];

                    if (
                        better_normalized(
                            candidate,
                            strongest
                        )
                    ) {
                        strongest = candidate;
                        strongest_m = m;
                    }
                }

                m5_competitor_histogram[
                    strongest_m
                ]++;

                if (
                    better_than(
                        global,
                        strongest
                    )
                ) {
                    ++m5_all_allowed_beaten;
                }

                const u64 margin =
                    score_cross_margin(
                        global,
                        strongest
                    );

                m5_score_margin_histogram[
                    margin
                ]++;

                if (
                    m5_margin_min == 0 ||
                    margin < m5_margin_min
                ) {
                    m5_margin_min = margin;
                }

                if (margin > m5_margin_max) {
                    m5_margin_max = margin;
                }

                /*
                 * Compare threshold with the actual r.
                 */
                if (
                    better_normalized(
                        global,
                        strongest
                    )
                ) {
                    const u64 threshold =
                        exact_threshold(
                            r,
                            global,
                            strongest
                        );

                    m5_threshold_histogram[
                        threshold
                    ]++;

                    if (
                        threshold <= r
                    ) {
                        ++m5_threshold_consistent;
                    }

                    if (
                        threshold > m5_threshold_max
                    ) {
                        m5_threshold_max =
                            threshold;
                    }

                    if (
                        threshold > r
                    ) {
                        ++m5_sharp_boundary;
                    }
                }
            } else {
                ++m7_points;

                if (!have_first_m7) {
                    first_m7.r = r;
                    first_m7.K = K;
                    first_m7.winner = global;
                    first_m7.allowed = allowed;
                    have_first_m7 = true;
                }

                if (global.k != last_m7_k) {
                    ++m7_new_record_points;
                    last_m7_k = global.k;
                } else {
                    ++m7_same_record_as_previous;
                }

                const Witness& w1 =
                    best_by_m[1];

                const Witness& w2 =
                    best_by_m[2];

                const Witness& w3 =
                    best_by_m[3];

                const Witness& w4 =
                    best_by_m[4];

                const Witness& w6 =
                    best_by_m[6];

                if (better_than(global, w1)) {
                    ++m7_beat_m1;
                }

                if (better_than(global, w2)) {
                    ++m7_beat_m2;
                }

                if (better_than(global, w3)) {
                    ++m7_beat_m3;
                }

                if (better_than(global, w4)) {
                    ++m7_beat_m4;
                }

                if (better_than(global, w6)) {
                    ++m7_beat_m6;
                }

                Witness strongest = w1;
                int strongest_m = 1;

                for (int m : {2, 3, 4, 6}) {
                    const Witness& candidate =
                        best_by_m[
                            static_cast<std::size_t>(m)
                        ];

                    if (
                        better_normalized(
                            candidate,
                            strongest
                        )
                    ) {
                        strongest = candidate;
                        strongest_m = m;
                    }
                }

                m7_competitor_histogram[
                    strongest_m
                ]++;

                if (
                    better_than(
                        global,
                        strongest
                    )
                ) {
                    ++m7_all_allowed_beaten;
                }

                const u64 margin =
                    score_cross_margin(
                        global,
                        strongest
                    );

                m7_score_margin_histogram[
                    margin
                ]++;

                if (
                    m7_margin_min == 0 ||
                    margin < m7_margin_min
                ) {
                    m7_margin_min = margin;
                }

                if (margin > m7_margin_max) {
                    m7_margin_max = margin;
                }

                if (
                    better_normalized(
                        global,
                        strongest
                    )
                ) {
                    const u64 threshold =
                        exact_threshold(
                            r,
                            global,
                            strongest
                        );

                    m7_threshold_histogram[
                        threshold
                    ]++;

                    if (
                        threshold <= r
                    ) {
                        ++m7_threshold_consistent;
                    }

                    if (
                        threshold > m7_threshold_max
                    ) {
                        m7_threshold_max =
                            threshold;
                    }

                    if (
                        threshold > r
                    ) {
                        ++m7_sharp_boundary;
                    }
                }
            }
        }
    }

    std::cout
        << "\nTOTAL_EXCEPTIONAL_POINTS="
        << exceptional_points
        << "\n";

    std::cout
        << "M5_POINTS="
        << m5_points
        << "\n";

    std::cout
        << "M7_POINTS="
        << m7_points
        << "\n";

    std::cout
        << "\nM5_COMPETITION\n";

    std::cout
        << "BEATS_M1="
        << m5_beat_m1
        << "\n";

    std::cout
        << "BEATS_M2="
        << m5_beat_m2
        << "\n";

    std::cout
        << "BEATS_M3="
        << m5_beat_m3
        << "\n";

    std::cout
        << "BEATS_M4="
        << m5_beat_m4
        << "\n";

    std::cout
        << "BEATS_M6="
        << m5_beat_m6
        << "\n";

    std::cout
        << "BEATS_ALL_ALLOWED="
        << m5_all_allowed_beaten
        << "\n";

    std::cout
        << "NEW_RECORD_POINTS="
        << m5_new_record_points
        << "\n";

    std::cout
        << "SAME_RECORD_AS_PREVIOUS="
        << m5_same_record_as_previous
        << "\n";

    std::cout
        << "THRESHOLD_CONSISTENT="
        << m5_threshold_consistent
        << "\n";

    std::cout
        << "SHARP_BOUNDARY="
        << m5_sharp_boundary
        << "\n";

    std::cout
        << "THRESHOLD_MAX="
        << m5_threshold_max
        << "\n";

    std::cout
        << "MARGIN_MIN="
        << m5_margin_min
        << "\n";

    std::cout
        << "MARGIN_MAX="
        << m5_margin_max
        << "\n";

    std::cout
        << "\nM7_COMPETITION\n";

    std::cout
        << "BEATS_M1="
        << m7_beat_m1
        << "\n";

    std::cout
        << "BEATS_M2="
        << m7_beat_m2
        << "\n";

    std::cout
        << "BEATS_M3="
        << m7_beat_m3
        << "\n";

    std::cout
        << "BEATS_M4="
        << m7_beat_m4
        << "\n";

    std::cout
        << "BEATS_M6="
        << m7_beat_m6
        << "\n";

    std::cout
        << "BEATS_ALL_ALLOWED="
        << m7_all_allowed_beaten
        << "\n";

    std::cout
        << "NEW_RECORD_POINTS="
        << m7_new_record_points
        << "\n";

    std::cout
        << "SAME_RECORD_AS_PREVIOUS="
        << m7_same_record_as_previous
        << "\n";

    std::cout
        << "THRESHOLD_CONSISTENT="
        << m7_threshold_consistent
        << "\n";

    std::cout
        << "SHARP_BOUNDARY="
        << m7_sharp_boundary
        << "\n";

    std::cout
        << "THRESHOLD_MAX="
        << m7_threshold_max
        << "\n";

    std::cout
        << "MARGIN_MIN="
        << m7_margin_min
        << "\n";

    std::cout
        << "MARGIN_MAX="
        << m7_margin_max
        << "\n";

    std::cout
        << "\nSTRONGEST_ALLOWED_COMPETITOR_M5\n";

    for (
        const auto& entry :
        m5_competitor_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nSTRONGEST_ALLOWED_COMPETITOR_M7\n";

    for (
        const auto& entry :
        m7_competitor_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM5_THRESHOLD_HISTOGRAM\n";

    for (
        const auto& entry :
        m5_threshold_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM7_THRESHOLD_HISTOGRAM\n";

    for (
        const auto& entry :
        m7_threshold_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM5_SCORE_MARGIN_HISTOGRAM\n";

    for (
        const auto& entry :
        m5_score_margin_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM7_SCORE_MARGIN_HISTOGRAM\n";

    for (
        const auto& entry :
        m7_score_margin_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_m5) {
        std::cout
            << "\nFIRST_M5_EXCEPTION\n";

        std::cout
            << "PRIME="
            << first_m5.r
            << "\n";

        std::cout
            << "K="
            << first_m5.K
            << "\n";

        print_witness(
            first_m5.r,
            "WINNER",
            first_m5.winner
        );

        std::cout << "\n";

        print_witness(
            first_m5.r,
            "ALLOWED",
            first_m5.allowed
        );

        std::cout << "\n";
    } else {
        std::cout
            << "\nFIRST_M5_EXCEPTION NONE\n";
    }

    if (have_first_m7) {
        std::cout
            << "\nFIRST_M7_EXCEPTION\n";

        std::cout
            << "PRIME="
            << first_m7.r
            << "\n";

        std::cout
            << "K="
            << first_m7.K
            << "\n";

        print_witness(
            first_m7.r,
            "WINNER",
            first_m7.winner
        );

        std::cout << "\n";

        print_witness(
            first_m7.r,
            "ALLOWED",
            first_m7.allowed
        );

        std::cout << "\n";
    } else {
        std::cout
            << "\nFIRST_M7_EXCEPTION NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
