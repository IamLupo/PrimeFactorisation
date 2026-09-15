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

struct Transition {
    u64 k = 0;
    Witness old_state;
    Witness new_state;
};

struct Excursion {
    u64 prime = 0;
    std::vector<Transition> transitions;
    std::vector<Witness> states;
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

static i128 arithmetic_value(
    u64 r,
    const Witness& w
) {
    return
        static_cast<i128>(w.m) *
        static_cast<i128>(r) +
        static_cast<i128>(w.sign);
}

static bool same_stream(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.sign == b.sign;
}

static bool same_m(
    const Witness& a,
    const Witness& b
) {
    return a.m == b.m;
}

static u64 abs_diff(
    u64 a,
    u64 b
) {
    return a >= b ? a - b : b - a;
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

static void print_transition(
    u64 r,
    const Transition& tr
) {
    std::cout
        << "K="
        << tr.k
        << " ";

    print_witness(
        r,
        "OLD",
        tr.old_state
    );

    std::cout
        << " -> ";

    print_witness(
        r,
        "NEW",
        tr.new_state
    );

    std::cout
        << "\n";
}

static void print_excursion(
    const Excursion& ex
) {
    std::cout
        << "PRIME="
        << ex.prime
        << "\n";

    std::cout
        << "STATES\n";

    for (
        std::size_t i = 0;
        i < ex.states.size();
        ++i
    ) {
        std::cout
            << "  STATE["
            << i
            << "] ";

        print_witness(
            ex.prime,
            "W",
            ex.states[i]
        );

        std::cout
            << "\n";
    }

    std::cout
        << "TRANSITIONS\n";

    for (
        const Transition& tr :
        ex.transitions
    ) {
        std::cout
            << "  ";

        print_transition(
            ex.prime,
            tr
        );
    }
}

int main() {
    constexpr int EXPERIMENT = 430;
    constexpr int PRIME_LIMIT = 5000;
    constexpr int K_LIMIT = 2000;
    constexpr int M_MAX = 32;

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

    u64 total_excursions = 0;
    u64 total_mirrored_transition_pairs = 0;

    u64 same_old_m = 0;
    u64 same_new_m = 0;

    u64 same_old_stream = 0;
    u64 same_new_stream = 0;

    u64 same_old_value = 0;
    u64 same_new_value = 0;

    u64 equal_old_values = 0;
    u64 equal_new_values = 0;

    u64 transition_k_gap_equal = 0;
    u64 transition_t_gap_equal = 0;

    u64 product_identity_exact = 0;

    u64 complement_cross_identity = 0;

    u64 same_arithmetic_pair = 0;
    u64 different_arithmetic_pair = 0;

    u64 transition_interval_match = 0;

    std::map<u64, u64>
        k_gap_histogram;

    std::map<u64, u64>
        t_gap_histogram;

    std::map<u64, u64>
        old_value_difference_histogram;

    std::map<u64, u64>
        new_value_difference_histogram;

    Excursion first_transition_failure;
    bool have_first_transition_failure = false;

    Excursion first_complement_failure;
    bool have_first_complement_failure = false;

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

        std::vector<Witness> winner_by_k(
            static_cast<std::size_t>(K_LIMIT) + 1
        );

        winner_by_k[1] =
            global_winner(best_by_m);

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
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

            winner_by_k[
                static_cast<std::size_t>(K)
            ] =
                global_winner(best_by_m);
        }

        /*
         * Extract only m-transitions.
         */
        std::vector<Transition> all_transitions;
        all_transitions.reserve(256);

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            const Witness& old_state =
                winner_by_k[
                    static_cast<std::size_t>(K - 1)
                ];

            const Witness& new_state =
                winner_by_k[
                    static_cast<std::size_t>(K)
                ];

            if (
                old_state.m !=
                new_state.m
            ) {
                Transition tr;
                tr.k = K;
                tr.old_state = old_state;
                tr.new_state = new_state;

                all_transitions.push_back(tr);
            }
        }

        /*
         * Split transitions into 1 -> ... -> 1 excursions.
         */
        std::size_t transition_start = 0;

        while (
            transition_start <
            all_transitions.size()
        ) {
            /*
             * Find a transition whose NEW state
             * enters non-1 territory.
             */
            while (
                transition_start <
                all_transitions.size() &&
                all_transitions[
                    transition_start
                ].old_state.m != 1
            ) {
                ++transition_start;
            }

            if (
                transition_start >=
                all_transitions.size()
            ) {
                break;
            }

            std::size_t transition_end =
                transition_start;

            while (
                transition_end <
                all_transitions.size()
            ) {
                const Transition& tr =
                    all_transitions[
                        transition_end
                    ];

                if (
                    transition_end >
                    transition_start &&
                    tr.new_state.m == 1
                ) {
                    break;
                }

                ++transition_end;
            }

            if (
                transition_end <=
                transition_start ||
                transition_end >
                    all_transitions.size()
            ) {
                break;
            }

            /*
             * Recover the state sequence corresponding
             * to the transition sequence.
             */
            Excursion ex;
            ex.prime = r;

            const Transition& first_tr =
                all_transitions[
                    transition_start
                ];

            ex.states.push_back(
                first_tr.old_state
            );

            for (
                std::size_t p =
                    transition_start;
                p < transition_end;
                ++p
            ) {
                ex.transitions.push_back(
                    all_transitions[p]
                );

                ex.states.push_back(
                    all_transitions[p].new_state
                );
            }

            if (
                ex.states.size() < 3
            ) {
                ++transition_end;
                transition_start =
                    transition_end;
                continue;
            }

            ++total_excursions;

            /*
             * The state sequence itself is the familiar
             * projected palindrome.
             */
            const std::size_t n =
                ex.transitions.size();

            if (n >= 2) {
                for (
                    std::size_t i = 0;
                    i < n / 2;
                    ++i
                ) {
                    const Transition& outward =
                        ex.transitions[i];

                    const Transition& inward =
                        ex.transitions[n - 1 - i];

                    ++total_mirrored_transition_pairs;

                    /*
                     * Compare OLD states.
                     */
                    if (
                        outward.old_state.m ==
                        inward.new_state.m
                    ) {
                        ++same_old_m;
                    } else if (
                        !have_first_transition_failure
                    ) {
                        first_transition_failure =
                            ex;
                        have_first_transition_failure =
                            true;
                    }

                    /*
                     * Compare NEW states.
                     */
                    if (
                        outward.new_state.m ==
                        inward.old_state.m
                    ) {
                        ++same_new_m;
                    }

                    if (
                        same_stream(
                            outward.old_state,
                            inward.new_state
                        )
                    ) {
                        ++same_old_stream;
                    }

                    if (
                        same_stream(
                            outward.new_state,
                            inward.old_state
                        )
                    ) {
                        ++same_new_stream;
                    }

                    const i128 outward_old_value =
                        arithmetic_value(
                            r,
                            outward.old_state
                        );

                    const i128 inward_new_value =
                        arithmetic_value(
                            r,
                            inward.new_state
                        );

                    const i128 outward_new_value =
                        arithmetic_value(
                            r,
                            outward.new_state
                        );

                    const i128 inward_old_value =
                        arithmetic_value(
                            r,
                            inward.old_state
                        );

                    if (
                        outward_old_value ==
                        inward_new_value
                    ) {
                        ++equal_old_values;
                    } else {
                        const u64 d =
                            outward_old_value >=
                            inward_new_value
                                ? static_cast<u64>(
                                    outward_old_value -
                                    inward_new_value
                                )
                                : static_cast<u64>(
                                    inward_new_value -
                                    outward_old_value
                                );

                        old_value_difference_histogram[
                            d
                        ]++;
                    }

                    if (
                        outward_new_value ==
                        inward_old_value
                    ) {
                        ++equal_new_values;
                    } else {
                        const u64 d =
                            outward_new_value >=
                            inward_old_value
                                ? static_cast<u64>(
                                    outward_new_value -
                                    inward_old_value
                                )
                                : static_cast<u64>(
                                    inward_old_value -
                                    outward_new_value
                                );

                        new_value_difference_histogram[
                            d
                        ]++;
                    }

                    if (
                        outward.old_state.k <
                        inward.new_state.k
                    ) {
                        const u64 dg =
                            inward.new_state.k -
                            outward.old_state.k;

                        k_gap_histogram[dg]++;
                    }

                    if (
                        outward.old_state.t >
                        inward.new_state.t
                    ) {
                        const u64 dg =
                            outward.old_state.t -
                            inward.new_state.t;

                        t_gap_histogram[dg]++;
                    }

                    /*
                     * Arithmetic pair preservation:
                     *
                     * outward:
                     *     A -> B
                     *
                     * inward:
                     *     B' -> A'
                     *
                     * We test whether A=A' and B=B'.
                     */
                    if (
                        outward_old_value ==
                        inward_new_value &&
                        outward_new_value ==
                        inward_old_value
                    ) {
                        ++same_arithmetic_pair;
                    } else {
                        ++different_arithmetic_pair;
                    }

                    /*
                     * Exact divisor identity.
                     */
                    const i128 p1 =
                        static_cast<i128>(
                            outward.old_state.k
                        ) *
                        static_cast<i128>(
                            outward.old_state.t
                        );

                    const i128 p2 =
                        static_cast<i128>(
                            outward.new_state.k
                        ) *
                        static_cast<i128>(
                            outward.new_state.t
                        );

                    const i128 p3 =
                        static_cast<i128>(
                            inward.old_state.k
                        ) *
                        static_cast<i128>(
                            inward.old_state.t
                        );

                    const i128 p4 =
                        static_cast<i128>(
                            inward.new_state.k
                        ) *
                        static_cast<i128>(
                            inward.new_state.t
                        );

                    if (
                        p1 ==
                        outward_old_value &&
                        p2 ==
                        outward_new_value &&
                        p3 ==
                        inward_old_value &&
                        p4 ==
                        inward_new_value
                    ) {
                        ++product_identity_exact;
                    }

                    /*
                     * Complementary cross test:
                     *
                     * outward OLD value / outward OLD k
                     * equals outward OLD t.
                     *
                     * Compare that t with the mirrored
                     * complementary divisor.
                     */
                    if (
                        outward.old_state.t ==
                            inward.new_state.t &&
                        outward.new_state.t ==
                            inward.old_state.t
                    ) {
                        ++complement_cross_identity;
                    }

                    /*
                     * Transition interval comparison.
                     */
                    const u64 outward_interval =
                        outward.new_state.k -
                        outward.old_state.k;

                    const u64 inward_interval =
                        inward.new_state.k -
                        inward.old_state.k;

                    if (
                        outward_interval ==
                        inward_interval
                    ) {
                        ++transition_interval_match;
                    }
                }
            }

            transition_start =
                transition_end;
        }
    }

    std::cout
        << "\nTOTAL_EXCURSIONS="
        << total_excursions
        << "\n";

    std::cout
        << "TOTAL_MIRRORED_TRANSITION_PAIRS="
        << total_mirrored_transition_pairs
        << "\n";

    std::cout
        << "\nM_STATE_MIRROR\n";

    std::cout
        << "OLD_M_TO_MIRROR_NEW="
        << same_old_m
        << "\n";

    std::cout
        << "NEW_M_TO_MIRROR_OLD="
        << same_new_m
        << "\n";

    std::cout
        << "\nSIGNED_STREAM_MIRROR\n";

    std::cout
        << "OLD_STREAM_MATCH="
        << same_old_stream
        << "\n";

    std::cout
        << "NEW_STREAM_MATCH="
        << same_new_stream
        << "\n";

    std::cout
        << "\nARITHMETIC_VALUE_MIRROR\n";

    std::cout
        << "OLD_VALUES_EQUAL="
        << equal_old_values
        << "\n";

    std::cout
        << "NEW_VALUES_EQUAL="
        << equal_new_values
        << "\n";

    std::cout
        << "SAME_ARITHMETIC_PAIR="
        << same_arithmetic_pair
        << "\n";

    std::cout
        << "DIFFERENT_ARITHMETIC_PAIR="
        << different_arithmetic_pair
        << "\n";

    std::cout
        << "\nALGEBRAIC_CHECKS\n";

    std::cout
        << "PRODUCT_IDENTITY_EXACT="
        << product_identity_exact
        << "\n";

    std::cout
        << "COMPLEMENT_CROSS_IDENTITY="
        << complement_cross_identity
        << "\n";

    std::cout
        << "TRANSITION_INTERVAL_MATCH="
        << transition_interval_match
        << "\n";

    std::cout
        << "\nK_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        k_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nT_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        t_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nOLD_VALUE_DIFFERENCE_HISTOGRAM\n";

    for (
        const auto& entry :
        old_value_difference_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nNEW_VALUE_DIFFERENCE_HISTOGRAM\n";

    for (
        const auto& entry :
        new_value_difference_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_transition_failure) {
        std::cout
            << "\nFIRST_TRANSITION_MIRROR_FAILURE\n";

        print_excursion(
            first_transition_failure
        );
    } else {
        std::cout
            << "\nFIRST_TRANSITION_MIRROR_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
