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
    u64 K = 0;
    Witness old_state;
    Witness new_state;
};

struct Excursion {
    u64 prime = 0;
    std::vector<Witness> states;
    std::vector<Transition> transitions;
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
    primes.reserve(static_cast<std::size_t>(limit / 2));

    for (int x = 2; x <= limit; ++x) {
        if (is_prime[x]) {
            primes.push_back(static_cast<u64>(x));
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

static bool projected_palindrome(
    const std::vector<Witness>& states
) {
    if (states.empty()) {
        return true;
    }

    for (
        std::size_t i = 0,
        j = states.size() - 1;
        i < j;
        ++i, --j
    ) {
        if (states[i].m != states[j].m) {
            return false;
        }
    }

    return true;
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

static bool factor_reversal(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.sign == b.sign &&
        a.k == b.t &&
        a.t == b.k;
}

static bool same_stream(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.sign == b.sign;
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
    const i128 value =
        arithmetic_value(r, w);

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

static void print_transition(
    u64 r,
    const Transition& tr
) {
    std::cout
        << "  K="
        << tr.K
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

    std::cout << "\n";
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

        std::cout << "\n";
    }

    std::cout
        << "TRANSITIONS\n";

    for (
        const Transition& tr :
        ex.transitions
    ) {
        print_transition(
            ex.prime,
            tr
        );
    }
}

static bool transition_pair_factor_reversal(
    const Transition& outward,
    const Transition& inward
) {
    return
        factor_reversal(
            outward.new_state,
            inward.old_state
        ) &&
        factor_reversal(
            outward.old_state,
            inward.new_state
        );
}

static bool transition_pair_stream_reversal(
    const Transition& outward,
    const Transition& inward
) {
    return
        same_stream(
            outward.new_state,
            inward.old_state
        ) &&
        same_stream(
            outward.old_state,
            inward.new_state
        );
}

static bool transition_pair_m_reversal(
    const Transition& outward,
    const Transition& inward
) {
    return
        outward.new_state.m ==
            inward.old_state.m &&
        outward.old_state.m ==
            inward.new_state.m;
}

static bool transition_pair_sign_reversal(
    const Transition& outward,
    const Transition& inward
) {
    return
        outward.new_state.sign ==
            inward.old_state.sign &&
        outward.old_state.sign ==
            inward.new_state.sign;
}

int main() {
    constexpr int EXPERIMENT = 431;
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

    u64 factor_reversal_matches = 0;
    u64 factor_reversal_failures = 0;

    u64 stream_matches = 0;
    u64 stream_failures = 0;

    u64 m_reversal_matches = 0;
    u64 m_reversal_failures = 0;

    u64 sign_matches = 0;
    u64 sign_failures = 0;

    u64 k_reversal_matches = 0;
    u64 k_reversal_failures = 0;

    u64 t_reversal_matches = 0;
    u64 t_reversal_failures = 0;

    u64 exact_arithmetic_pair_matches = 0;
    u64 arithmetic_pair_failures = 0;

    u64 value_match_old = 0;
    u64 value_match_new = 0;

    u64 value_difference_two_old = 0;
    u64 value_difference_two_new = 0;

    u64 value_other_difference = 0;

    u64 interval_gap_equal = 0;
    u64 interval_gap_different = 0;

    std::map<u64, u64>
        k_difference_histogram;

    std::map<u64, u64>
        t_difference_histogram;

    std::map<u64, u64>
        K_gap_histogram;

    std::map<u64, u64>
        old_value_difference_histogram;

    std::map<u64, u64>
        new_value_difference_histogram;

    Excursion first_factor_failure;
    bool have_first_factor_failure = false;

    Excursion first_stream_failure;
    bool have_first_stream_failure = false;

    Excursion first_k_failure;
    bool have_first_k_failure = false;

    Excursion first_value_failure;
    bool have_first_value_failure = false;

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
         * Complete winner state for every K.
         */
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

                /*
                 * Both signs must be considered independently.
                 */
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
         * Extract all m-transitions.
         */
        std::vector<Transition> transitions;
        transitions.reserve(256);

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

                tr.K = K;
                tr.old_state = old_state;
                tr.new_state = new_state;

                transitions.push_back(tr);
            }
        }

        /*
         * Find every complete 1 -> ... -> 1 excursion.
         */
        std::size_t p = 0;

        while (p < transitions.size()) {
            while (
                p < transitions.size() &&
                transitions[p].old_state.m != 1
            ) {
                ++p;
            }

            if (p >= transitions.size()) {
                break;
            }

            const std::size_t begin = p;

            while (
                p < transitions.size() &&
                transitions[p].new_state.m != 1
            ) {
                ++p;
            }

            if (p >= transitions.size()) {
                break;
            }

            const std::size_t end = p;

            Excursion ex;
            ex.prime = r;

            for (
                std::size_t q = begin;
                q <= end;
                ++q
            ) {
                ex.transitions.push_back(
                    transitions[q]
                );
            }

            if (ex.transitions.empty()) {
                p = end + 1;
                continue;
            }

            ex.states.push_back(
                ex.transitions.front().old_state
            );

            for (
                const Transition& tr :
                ex.transitions
            ) {
                ex.states.push_back(
                    tr.new_state
                );
            }

            if (ex.states.size() < 3) {
                p = end + 1;
                continue;
            }

            /*
             * All excursions in the previous experiments
             * satisfy this projected palindrome.
             */
            if (
                !projected_palindrome(
                    ex.states
                )
            ) {
                p = end + 1;
                continue;
            }

            ++total_excursions;

            const std::size_t n =
                ex.transitions.size();

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

                const bool m_ok =
                    transition_pair_m_reversal(
                        outward,
                        inward
                    );

                if (m_ok) {
                    ++m_reversal_matches;
                } else {
                    ++m_reversal_failures;
                }

                const bool stream_ok =
                    transition_pair_stream_reversal(
                        outward,
                        inward
                    );

                if (stream_ok) {
                    ++stream_matches;
                } else {
                    ++stream_failures;

                    if (!have_first_stream_failure) {
                        first_stream_failure =
                            ex;

                        have_first_stream_failure =
                            true;
                    }
                }

                const bool sign_ok =
                    transition_pair_sign_reversal(
                        outward,
                        inward
                    );

                if (sign_ok) {
                    ++sign_matches;
                } else {
                    ++sign_failures;
                }

                const bool factor_ok =
                    transition_pair_factor_reversal(
                        outward,
                        inward
                    );

                if (factor_ok) {
                    ++factor_reversal_matches;
                } else {
                    ++factor_reversal_failures;

                    if (!have_first_factor_failure) {
                        first_factor_failure =
                            ex;

                        have_first_factor_failure =
                            true;
                    }
                }

                /*
                 * k,t factor reversal separately.
                 */
                const bool k_ok =
                    outward.new_state.k ==
                        inward.old_state.t &&
                    outward.old_state.k ==
                        inward.new_state.t;

                if (k_ok) {
                    ++k_reversal_matches;
                } else {
                    ++k_reversal_failures;

                    if (!have_first_k_failure) {
                        first_k_failure =
                            ex;

                        have_first_k_failure =
                            true;
                    }
                }

                const bool t_ok =
                    outward.new_state.t ==
                        inward.old_state.k &&
                    outward.old_state.t ==
                        inward.new_state.k;

                if (t_ok) {
                    ++t_reversal_matches;
                } else {
                    ++t_reversal_failures;
                }

                /*
                 * Arithmetic values.
                 *
                 * The correct mirrored comparison is:
                 *
                 * outward OLD <-> inward NEW
                 * outward NEW <-> inward OLD
                 */
                const i128 outward_old_value =
                    arithmetic_value(
                        r,
                        outward.old_state
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

                const i128 inward_new_value =
                    arithmetic_value(
                        r,
                        inward.new_state
                    );

                const u64 old_difference =
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

                const u64 new_difference =
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

                if (old_difference == 0) {
                    ++value_match_old;
                } else if (old_difference == 2) {
                    ++value_difference_two_old;
                } else {
                    old_value_difference_histogram[
                        old_difference
                    ]++;

                    ++value_other_difference;
                }

                if (new_difference == 0) {
                    ++value_match_new;
                } else if (new_difference == 2) {
                    ++value_difference_two_new;
                } else {
                    new_value_difference_histogram[
                        new_difference
                    ]++;

                    ++value_other_difference;
                }

                if (
                    old_difference == 0 &&
                    new_difference == 0
                ) {
                    ++exact_arithmetic_pair_matches;
                } else {
                    ++arithmetic_pair_failures;

                    if (!have_first_value_failure) {
                        first_value_failure = ex;
                        have_first_value_failure = true;
                    }
                }

                /*
                 * k/t difference geometry.
                 */
                const u64 k_difference =
                    abs_diff(
                        outward.new_state.k,
                        inward.old_state.k
                    );

                const u64 t_difference =
                    abs_diff(
                        outward.new_state.t,
                        inward.old_state.t
                    );

                k_difference_histogram[
                    k_difference
                ]++;

                t_difference_histogram[
                    t_difference
                ]++;

                /*
                 * Distance between the transition K
                 * coordinates.
                 */
                const u64 K_gap =
                    inward.K >= outward.K
                        ? inward.K - outward.K
                        : outward.K - inward.K;

                K_gap_histogram[K_gap]++;

                if (
                    outward.K >= outward.old_state.k &&
                    inward.K >= inward.old_state.k
                ) {
                    const u64 outward_interval =
                        outward.K -
                        outward.old_state.k;

                    const u64 inward_interval =
                        inward.K -
                        inward.old_state.k;

                    if (
                        outward_interval ==
                        inward_interval
                    ) {
                        ++interval_gap_equal;
                    } else {
                        ++interval_gap_different;
                    }
                }
            }

            p = end + 1;
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
        << "\nM_REVERSAL\n";

    std::cout
        << "M_REVERSAL_MATCHES="
        << m_reversal_matches
        << "\n";

    std::cout
        << "M_REVERSAL_FAILURES="
        << m_reversal_failures
        << "\n";

    std::cout
        << "\nSTREAM_REVERSAL\n";

    std::cout
        << "STREAM_MATCHES="
        << stream_matches
        << "\n";

    std::cout
        << "STREAM_FAILURES="
        << stream_failures
        << "\n";

    std::cout
        << "\nSIGN_REVERSAL\n";

    std::cout
        << "SIGN_MATCHES="
        << sign_matches
        << "\n";

    std::cout
        << "SIGN_FAILURES="
        << sign_failures
        << "\n";

    std::cout
        << "\nFACTOR_REVERSAL\n";

    std::cout
        << "FACTOR_REVERSAL_MATCHES="
        << factor_reversal_matches
        << "\n";

    std::cout
        << "FACTOR_REVERSAL_FAILURES="
        << factor_reversal_failures
        << "\n";

    std::cout
        << "\nK_REVERSAL\n";

    std::cout
        << "K_REVERSAL_MATCHES="
        << k_reversal_matches
        << "\n";

    std::cout
        << "K_REVERSAL_FAILURES="
        << k_reversal_failures
        << "\n";

    std::cout
        << "\nT_REVERSAL\n";

    std::cout
        << "T_REVERSAL_MATCHES="
        << t_reversal_matches
        << "\n";

    std::cout
        << "T_REVERSAL_FAILURES="
        << t_reversal_failures
        << "\n";

    std::cout
        << "\nARITHMETIC_PAIR\n";

    std::cout
        << "EXACT_ARITHMETIC_PAIR_MATCHES="
        << exact_arithmetic_pair_matches
        << "\n";

    std::cout
        << "ARITHMETIC_PAIR_FAILURES="
        << arithmetic_pair_failures
        << "\n";

    std::cout
        << "VALUE_MATCH_OLD="
        << value_match_old
        << "\n";

    std::cout
        << "VALUE_MATCH_NEW="
        << value_match_new
        << "\n";

    std::cout
        << "VALUE_DIFFERENCE_TWO_OLD="
        << value_difference_two_old
        << "\n";

    std::cout
        << "VALUE_DIFFERENCE_TWO_NEW="
        << value_difference_two_new
        << "\n";

    std::cout
        << "VALUE_OTHER_DIFFERENCE="
        << value_other_difference
        << "\n";

    std::cout
        << "\nTRANSITION_INTERVALS\n";

    std::cout
        << "INTERVAL_GAP_EQUAL="
        << interval_gap_equal
        << "\n";

    std::cout
        << "INTERVAL_GAP_DIFFERENT="
        << interval_gap_different
        << "\n";

    std::cout
        << "\nK_DIFFERENCE_HISTOGRAM\n";

    for (
        const auto& entry :
        k_difference_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nT_DIFFERENCE_HISTOGRAM\n";

    for (
        const auto& entry :
        t_difference_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nK_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        K_gap_histogram
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

    if (have_first_factor_failure) {
        std::cout
            << "\nFIRST_FACTOR_REVERSAL_FAILURE\n";

        print_excursion(
            first_factor_failure
        );
    } else {
        std::cout
            << "\nFIRST_FACTOR_REVERSAL_FAILURE NONE\n";
    }

    if (have_first_stream_failure) {
        std::cout
            << "\nFIRST_STREAM_REVERSAL_FAILURE\n";

        print_excursion(
            first_stream_failure
        );
    } else {
        std::cout
            << "\nFIRST_STREAM_REVERSAL_FAILURE NONE\n";
    }

    if (have_first_k_failure) {
        std::cout
            << "\nFIRST_K_REVERSAL_FAILURE\n";

        print_excursion(
            first_k_failure
        );
    } else {
        std::cout
            << "\nFIRST_K_REVERSAL_FAILURE NONE\n";
    }

    if (have_first_value_failure) {
        std::cout
            << "\nFIRST_ARITHMETIC_PAIR_FAILURE\n";

        print_excursion(
            first_value_failure
        );
    } else {
        std::cout
            << "\nFIRST_ARITHMETIC_PAIR_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}