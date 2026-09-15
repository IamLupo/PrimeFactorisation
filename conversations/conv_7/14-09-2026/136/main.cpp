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

struct Excursion {
    u64 prime = 0;
    std::vector<Witness> states;
};

struct KRecordInfo {
    u64 candidate_count = 0;
    u64 record_count = 0;
    u64 global_winner_changes = 0;
    u64 global_winner_m = 1;
    int global_winner_sign = +1;
    std::uint64_t record_mask = 0;
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

static void print_excursion(
    const Excursion& ex
) {
    std::cout
        << "PRIME="
        << ex.prime
        << " LENGTH="
        << ex.states.size()
        << "\n";

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
}

static std::string mask_to_string(
    std::uint64_t mask,
    int m_max
) {
    std::string out;

    for (int m = 1; m <= m_max; ++m) {
        if (
            (mask &
             (std::uint64_t(1) << (m - 1))) != 0
        ) {
            if (!out.empty()) {
                out += ",";
            }

            out += std::to_string(m);
        }
    }

    return out;
}

int main() {
    constexpr int EXPERIMENT = 428;
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
    u64 total_mirrored_pairs = 0;

    u64 hidden_candidate_points = 0;
    u64 hidden_per_m_records = 0;
    u64 hidden_global_winner_changes = 0;

    u64 hidden_candidate_not_record = 0;
    u64 hidden_record_not_global = 0;

    u64 same_value_pairs = 0;
    u64 difference_two_pairs = 0;

    u64 same_value_record_points = 0;
    u64 difference_two_record_points = 0;

    u64 hidden_m_same_as_endpoint = 0;
    u64 hidden_m_other = 0;

    u64 hidden_global_m_matches_endpoint = 0;
    u64 hidden_global_m_other = 0;

    std::map<u64, u64>
        candidate_count_histogram;

    std::map<u64, u64>
        record_count_histogram;

    std::map<u64, u64>
        global_change_count_histogram;

    std::map<int, u64>
        hidden_record_m_histogram;

    std::map<int, u64>
        hidden_candidate_m_histogram;

    Excursion first_hidden_candidate;
    bool have_first_hidden_candidate = false;

    Excursion first_hidden_record;
    bool have_first_hidden_record = false;

    Excursion first_hidden_global_change;
    bool have_first_hidden_global_change = false;

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
         * We keep the entire winner history for this prime.
         */
        std::vector<Witness> winner_by_k(
            static_cast<std::size_t>(K_LIMIT) + 1
        );

        std::vector<KRecordInfo> info_by_k(
            static_cast<std::size_t>(K_LIMIT) + 1
        );

        Witness previous =
            global_winner(best_by_m);

        winner_by_k[1] = previous;

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            KRecordInfo info;

            Witness previous_winner =
                previous;

            for (int m = 1; m <= M_MAX; ++m) {
                bool candidate_found = false;
                Witness candidate;

                if (
                    divides_stream(
                        r,
                        m,
                        K,
                        +1
                    )
                ) {
                    candidate =
                        make_witness(
                            r,
                            m,
                            K,
                            +1
                        );

                    candidate_found = true;
                    ++info.candidate_count;
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

                    ++info.candidate_count;

                    if (
                        !candidate_found ||
                        better_normalized(
                            minus_candidate,
                            candidate
                        )
                    ) {
                        candidate =
                            minus_candidate;

                        candidate_found = true;
                    }
                }

                if (!candidate_found) {
                    continue;
                }

                const std::size_t mi =
                    static_cast<std::size_t>(m);

                if (
                    candidate.k >
                    best_by_m[mi].k
                ) {
                    best_by_m[mi] =
                        candidate;

                    ++info.record_count;

                    info.record_mask |=
                        (std::uint64_t(1)
                         << (m - 1));
                }
            }

            const Witness current =
                global_winner(best_by_m);

            if (
                current.m !=
                    previous_winner.m ||
                current.k !=
                    previous_winner.k ||
                current.sign !=
                    previous_winner.sign
            ) {
                ++info.global_winner_changes;
            }

            info.global_winner_m =
                static_cast<u64>(current.m);

            info.global_winner_sign =
                current.sign;

            winner_by_k[
                static_cast<std::size_t>(K)
            ] = current;

            info_by_k[
                static_cast<std::size_t>(K)
            ] = info;

            previous = current;
        }

        /*
         * Extract projected winner transitions.
         */
        std::vector<Witness> transition_states;
        transition_states.reserve(128);

        transition_states.push_back(
            winner_by_k[1]
        );

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            const Witness& a =
                winner_by_k[
                    static_cast<std::size_t>(K - 1)
                ];

            const Witness& b =
                winner_by_k[
                    static_cast<std::size_t>(K)
                ];

            if (a.m != b.m) {
                transition_states.push_back(b);
            }
        }

        std::size_t start = 0;

        while (
            start < transition_states.size()
        ) {
            while (
                start < transition_states.size() &&
                transition_states[start].m != 1
            ) {
                ++start;
            }

            if (
                start >=
                transition_states.size()
            ) {
                break;
            }

            std::size_t end = start + 1;

            while (
                end < transition_states.size() &&
                transition_states[end].m != 1
            ) {
                ++end;
            }

            if (
                end >=
                transition_states.size()
            ) {
                break;
            }

            if (end > start + 1) {
                Excursion ex;
                ex.prime = r;

                ex.states.assign(
                    transition_states.begin() +
                        static_cast<std::ptrdiff_t>(start),
                    transition_states.begin() +
                        static_cast<std::ptrdiff_t>(end) +
                        1
                );

                ++total_excursions;

                if (
                    !projected_palindrome(
                        ex.states
                    )
                ) {
                    start = end;
                    continue;
                }

                for (
                    std::size_t i = 0;
                    i < ex.states.size() / 2;
                    ++i
                ) {
                    const std::size_t j =
                        ex.states.size() - 1 - i;

                    const Witness& left =
                        ex.states[i];

                    const Witness& right =
                        ex.states[j];

                    ++total_mirrored_pairs;

                    const u64 low_k =
                        left.k;

                    const u64 high_k =
                        right.k;

                    const i128 left_value =
                        arithmetic_value(
                            r,
                            left
                        );

                    const i128 right_value =
                        arithmetic_value(
                            r,
                            right
                        );

                    const u64 difference =
                        left_value >= right_value
                            ? static_cast<u64>(
                                left_value -
                                right_value
                            )
                            : static_cast<u64>(
                                right_value -
                                left_value
                            );

                    if (difference == 0) {
                        ++same_value_pairs;
                    } else if (difference == 2) {
                        ++difference_two_pairs;
                    }

                    /*
                     * Inspect EVERY K between the mirrored
                     * visible winner states.
                     */
                    for (
                        u64 K = low_k + 1;
                        K < high_k;
                        ++K
                    ) {
                        const KRecordInfo& info =
                            info_by_k[
                                static_cast<std::size_t>(K)
                            ];

                        if (
                            info.candidate_count == 0
                        ) {
                            continue;
                        }

                        ++hidden_candidate_points;

                        if (
                            info.record_count > 0
                        ) {
                            hidden_per_m_records +=
                                info.record_count;

                            if (
                                !have_first_hidden_record
                            ) {
                                first_hidden_record = ex;
                                have_first_hidden_record = true;
                            }
                        }

                        /*
                         * A candidate that is not a new
                         * per-m record.
                         */
                        if (
                            info.candidate_count >
                            info.record_count
                        ) {
                            hidden_candidate_not_record +=
                                info.candidate_count -
                                info.record_count;
                        }

                        /*
                         * Actual global winner change.
                         */
                        if (
                            info.global_winner_changes > 0
                        ) {
                            ++hidden_global_winner_changes;

                            if (
                                !have_first_hidden_global_change
                            ) {
                                first_hidden_global_change =
                                    ex;
                                have_first_hidden_global_change =
                                    true;
                            }
                        }

                        /*
                         * Every per-m record which does not
                         * produce a global winner change is
                         * a local record only.
                         */
                        if (
                            info.record_count > 0 &&
                            info.global_winner_changes == 0
                        ) {
                            hidden_record_not_global +=
                                info.record_count;
                        }

                        candidate_count_histogram[
                            info.candidate_count
                        ]++;

                        record_count_histogram[
                            info.record_count
                        ]++;

                        global_change_count_histogram[
                            info.global_winner_changes
                        ]++;

                        for (int m = 1;
                             m <= M_MAX;
                             ++m) {
                            if (
                                (
                                    info.record_mask &
                                    (std::uint64_t(1)
                                     << (m - 1))
                                ) == 0
                            ) {
                                continue;
                            }

                            ++hidden_record_m_histogram[m];

                            if (
                                m == left.m ||
                                m == right.m
                            ) {
                                ++hidden_m_same_as_endpoint;
                            } else {
                                ++hidden_m_other;
                            }

                            if (
                                difference == 0
                            ) {
                                ++same_value_record_points;
                            } else if (
                                difference == 2
                            ) {
                                ++difference_two_record_points;
                            }
                        }

                        /*
                         * The actual winner at this K.
                         */
                        const int winner_m =
                            static_cast<int>(
                                info.global_winner_m
                            );

                        if (
                            winner_m == left.m ||
                            winner_m == right.m
                        ) {
                            ++hidden_global_m_matches_endpoint;
                        } else {
                            ++hidden_global_m_other;
                        }

                        if (
                            !have_first_hidden_candidate
                        ) {
                            first_hidden_candidate = ex;
                            have_first_hidden_candidate =
                                true;
                        }
                    }
                }
            }

            start = end;
        }
    }

    std::cout
        << "\nTOTAL_EXCURSIONS="
        << total_excursions
        << "\n";

    std::cout
        << "TOTAL_MIRRORED_PAIRS="
        << total_mirrored_pairs
        << "\n";

    std::cout
        << "SAME_VALUE_PAIRS="
        << same_value_pairs
        << "\n";

    std::cout
        << "DIFFERENCE_TWO_PAIRS="
        << difference_two_pairs
        << "\n";

    std::cout
        << "\nCORRECTED_HIDDEN_K_ANALYSIS\n";

    std::cout
        << "HIDDEN_K_CANDIDATE_POINTS="
        << hidden_candidate_points
        << "\n";

    std::cout
        << "HIDDEN_PER_M_RECORDS="
        << hidden_per_m_records
        << "\n";

    std::cout
        << "HIDDEN_CANDIDATE_NOT_RECORD="
        << hidden_candidate_not_record
        << "\n";

    std::cout
        << "HIDDEN_RECORD_NOT_GLOBAL="
        << hidden_record_not_global
        << "\n";

    std::cout
        << "HIDDEN_GLOBAL_WINNER_CHANGES="
        << hidden_global_winner_changes
        << "\n";

    std::cout
        << "HIDDEN_M_SAME_AS_ENDPOINT="
        << hidden_m_same_as_endpoint
        << "\n";

    std::cout
        << "HIDDEN_M_OTHER="
        << hidden_m_other
        << "\n";

    std::cout
        << "HIDDEN_GLOBAL_M_MATCHES_ENDPOINT="
        << hidden_global_m_matches_endpoint
        << "\n";

    std::cout
        << "HIDDEN_GLOBAL_M_OTHER="
        << hidden_global_m_other
        << "\n";

    std::cout
        << "SAME_VALUE_RECORD_POINTS="
        << same_value_record_points
        << "\n";

    std::cout
        << "DIFFERENCE_TWO_RECORD_POINTS="
        << difference_two_record_points
        << "\n";

    std::cout
        << "\nCANDIDATE_COUNT_HISTOGRAM\n";

    for (
        const auto& entry :
        candidate_count_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nRECORD_COUNT_HISTOGRAM\n";

    for (
        const auto& entry :
        record_count_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nGLOBAL_CHANGE_COUNT_HISTOGRAM\n";

    for (
        const auto& entry :
        global_change_count_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nHIDDEN_RECORD_M_HISTOGRAM\n";

    for (
        const auto& entry :
        hidden_record_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_hidden_candidate) {
        std::cout
            << "\nFIRST_HIDDEN_CANDIDATE\n";

        print_excursion(
            first_hidden_candidate
        );
    } else {
        std::cout
            << "\nFIRST_HIDDEN_CANDIDATE NONE\n";
    }

    if (have_first_hidden_record) {
        std::cout
            << "\nFIRST_HIDDEN_PER_M_RECORD\n";

        print_excursion(
            first_hidden_record
        );
    } else {
        std::cout
            << "\nFIRST_HIDDEN_PER_M_RECORD NONE\n";
    }

    if (have_first_hidden_global_change) {
        std::cout
            << "\nFIRST_HIDDEN_GLOBAL_WINNER_CHANGE\n";

        print_excursion(
            first_hidden_global_change
        );
    } else {
        std::cout
            << "\nFIRST_HIDDEN_GLOBAL_WINNER_CHANGE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
