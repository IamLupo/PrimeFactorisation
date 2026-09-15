#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
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

struct Record {
    u64 k = 0;
    u64 t = 0;
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

static bool same_stream(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.sign == b.sign;
}

static std::vector<Record> divisor_records(
    u64 value
) {
    std::vector<Record> records;

    u64 best_k = 0;

    /*
     * Every divisor k <= value is checked.
     * Values here are <= 160001, so this is tiny.
     */
    for (u64 k = 1; k <= value; ++k) {
        if (value % k != 0) {
            continue;
        }

        if (k > best_k) {
            best_k = k;

            Record record;
            record.k = k;
            record.t = value / k;

            records.push_back(record);
        }
    }

    return records;
}

static int find_record_index(
    const std::vector<Record>& records,
    u64 k
) {
    for (
        std::size_t i = 0;
        i < records.size();
        ++i
    ) {
        if (records[i].k == k) {
            return static_cast<int>(i);
        }
    }

    return -1;
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

static void print_pair(
    u64 r,
    const Witness& left,
    const Witness& right
) {
    print_witness(
        r,
        "LEFT",
        left
    );

    std::cout << "\n";

    print_witness(
        r,
        "RIGHT",
        right
    );

    std::cout << "\n";
}

int main() {
    constexpr int EXPERIMENT = 434;
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

    u64 total_pairs = 0;
    u64 same_stream_pairs = 0;
    u64 sign_switch_pairs = 0;

    u64 all_four_are_records = 0;

    u64 left_k_record_matches = 0;
    u64 right_k_record_matches = 0;
    u64 left_t_record_matches = 0;
    u64 right_t_record_matches = 0;

    u64 reciprocal_index_matches = 0;

    u64 left_right_consecutive = 0;
    u64 t_left_right_consecutive = 0;
    u64 full_four_consecutive = 0;

    u64 record_gap_one = 0;
    u64 record_gap_two = 0;
    u64 record_gap_more = 0;

    u64 exact_complement_all = 0;

    std::map<u64, u64>
        record_count_histogram;

    std::map<u64, u64>
        reciprocal_index_gap_histogram;

    std::map<u64, u64>
        record_gap_histogram;

    std::vector<Witness> first_record_failure_pair;
    u64 first_record_failure_prime = 0;
    bool have_first_record_failure = false;

    std::vector<Witness> first_reciprocal_failure_pair;
    u64 first_reciprocal_failure_prime = 0;
    bool have_first_reciprocal_failure = false;

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

        std::vector<Witness> transition_states;
        transition_states.reserve(256);

        transition_states.push_back(
            winner_by_k[1]
        );

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            const Witness& previous =
                winner_by_k[
                    static_cast<std::size_t>(K - 1)
                ];

            const Witness& current =
                winner_by_k[
                    static_cast<std::size_t>(K)
                ];

            if (
                previous.m !=
                current.m
            ) {
                transition_states.push_back(
                    current
                );
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
                end >= transition_states.size()
            ) {
                break;
            }

            if (end > start + 1) {
                std::vector<Witness> excursion;

                excursion.assign(
                    transition_states.begin() +
                        static_cast<std::ptrdiff_t>(start),
                    transition_states.begin() +
                        static_cast<std::ptrdiff_t>(end) +
                        1
                );

                if (
                    !projected_palindrome(
                        excursion
                    )
                ) {
                    start = end;
                    continue;
                }

                const std::size_t n =
                    excursion.size();

                for (
                    std::size_t i = 0;
                    i < n / 2;
                    ++i
                ) {
                    const std::size_t j =
                        n - 1 - i;

                    const Witness& left =
                        excursion[i];

                    const Witness& right =
                        excursion[j];

                    ++total_pairs;

                    if (
                        !same_stream(
                            left,
                            right
                        )
                    ) {
                        ++sign_switch_pairs;
                        continue;
                    }

                    ++same_stream_pairs;

                    const u64 A =
                        static_cast<u64>(
                            arithmetic_value(
                                r,
                                left
                            )
                        );

                    const std::vector<Record> records =
                        divisor_records(A);

                    record_count_histogram[
                        static_cast<u64>(
                            records.size()
                        )
                    ]++;

                    const int left_k_index =
                        find_record_index(
                            records,
                            left.k
                        );

                    const int right_k_index =
                        find_record_index(
                            records,
                            right.k
                        );

                    const int left_t_index =
                        find_record_index(
                            records,
                            left.t
                        );

                    const int right_t_index =
                        find_record_index(
                            records,
                            right.t
                        );

                    const bool all_records =
                        left_k_index >= 0 &&
                        right_k_index >= 0 &&
                        left_t_index >= 0 &&
                        right_t_index >= 0;

                    if (all_records) {
                        ++all_four_are_records;
                    } else {
                        if (!have_first_record_failure) {
                            first_record_failure_pair =
                                excursion;

                            first_record_failure_prime =
                                r;

                            have_first_record_failure =
                                true;
                        }
                    }

                    if (left_k_index >= 0) {
                        ++left_k_record_matches;
                    }

                    if (right_k_index >= 0) {
                        ++right_k_record_matches;
                    }

                    if (left_t_index >= 0) {
                        ++left_t_record_matches;
                    }

                    if (right_t_index >= 0) {
                        ++right_t_record_matches;
                    }

                    if (!all_records) {
                        continue;
                    }

                    /*
                     * Reciprocal divisor involution:
                     *
                     * A / left.k = left.t
                     * A / right.k = right.t
                     *
                     * Therefore the record containing left.t
                     * should be the reciprocal record of left.k,
                     * and similarly for right.
                     */
                    const bool reciprocal_ok =
                        left_t_index == right_k_index &&
                        right_t_index == left_k_index;

                    if (reciprocal_ok) {
                        ++reciprocal_index_matches;
                    } else if (
                        !have_first_reciprocal_failure
                    ) {
                        first_reciprocal_failure_pair =
                            excursion;

                        first_reciprocal_failure_prime =
                            r;

                        have_first_reciprocal_failure =
                            true;
                    }

                    const u64 reciprocal_gap =
                        abs_diff(
                            static_cast<u64>(
                                left_k_index
                            ),
                            static_cast<u64>(
                                right_k_index
                            )
                        );

                    reciprocal_index_gap_histogram[
                        reciprocal_gap
                    ]++;

                    /*
                     * Are the two visible k values
                     * consecutive divisor records?
                     */
                    if (
                        right_k_index ==
                        left_k_index + 1
                    ) {
                        ++left_right_consecutive;
                        ++record_gap_one;
                    } else if (
                        right_k_index ==
                        left_k_index + 2
                    ) {
                        ++record_gap_two;
                    } else if (
                        right_k_index >
                        left_k_index + 2
                    ) {
                        ++record_gap_more;
                    }

                    /*
                     * The reciprocal t-pair should have
                     * the same record separation.
                     */
                    if (
                        left_t_index ==
                        right_k_index &&
                        right_t_index ==
                        left_k_index
                    ) {
                        ++t_left_right_consecutive;
                    }

                    if (
                        right_k_index ==
                        left_k_index + 1 &&
                        left_t_index ==
                        right_k_index &&
                        right_t_index ==
                        left_k_index
                    ) {
                        ++full_four_consecutive;
                    }

                    const u64 gap =
                        abs_diff(
                            static_cast<u64>(
                                left_k_index
                            ),
                            static_cast<u64>(
                                right_k_index
                            )
                        );

                    record_gap_histogram[
                        gap
                    ]++;

                    /*
                     * Both visible states are exact
                     * complementary divisor pairs.
                     */
                    const bool complement_ok =
                        A % left.k == 0 &&
                        A / left.k == left.t &&
                        A % right.k == 0 &&
                        A / right.k == right.t;

                    if (complement_ok) {
                        ++exact_complement_all;
                    }
                }
            }

            start = end;
        }
    }

    std::cout
        << "\nTOTAL_MIRRORED_PAIRS="
        << total_pairs
        << "\n";

    std::cout
        << "SAME_STREAM_PAIRS="
        << same_stream_pairs
        << "\n";

    std::cout
        << "SIGN_SWITCH_PAIRS="
        << sign_switch_pairs
        << "\n";

    std::cout
        << "\nRECORD_MEMBERSHIP\n";

    std::cout
        << "ALL_FOUR_ARE_RECORDS="
        << all_four_are_records
        << "\n";

    std::cout
        << "LEFT_K_RECORDS="
        << left_k_record_matches
        << "\n";

    std::cout
        << "RIGHT_K_RECORDS="
        << right_k_record_matches
        << "\n";

    std::cout
        << "LEFT_T_RECORDS="
        << left_t_record_matches
        << "\n";

    std::cout
        << "RIGHT_T_RECORDS="
        << right_t_record_matches
        << "\n";

    std::cout
        << "\nRECIPROCAL_INVOLUTION\n";

    std::cout
        << "RECIPROCAL_INDEX_MATCHES="
        << reciprocal_index_matches
        << "\n";

    std::cout
        << "\nRECORD_ADJACENCY\n";

    std::cout
        << "LEFT_RIGHT_CONSECUTIVE="
        << left_right_consecutive
        << "\n";

    std::cout
        << "T_LEFT_RIGHT_CONSECUTIVE="
        << t_left_right_consecutive
        << "\n";

    std::cout
        << "FULL_FOUR_CONSECUTIVE="
        << full_four_consecutive
        << "\n";

    std::cout
        << "RECORD_GAP_ONE="
        << record_gap_one
        << "\n";

    std::cout
        << "RECORD_GAP_TWO="
        << record_gap_two
        << "\n";

    std::cout
        << "RECORD_GAP_MORE="
        << record_gap_more
        << "\n";

    std::cout
        << "\nEXACT_COMPLEMENT\n";

    std::cout
        << "EXACT_COMPLEMENT_ALL="
        << exact_complement_all
        << "\n";

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
        << "\nRECIPROCAL_INDEX_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        reciprocal_index_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nRECORD_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        record_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_record_failure) {
        std::cout
            << "\nFIRST_RECORD_MEMBERSHIP_FAILURE\n";

        print_pair(
            first_record_failure_prime,
            first_record_failure_pair.front(),
            first_record_failure_pair.back()
        );
    } else {
        std::cout
            << "\nFIRST_RECORD_MEMBERSHIP_FAILURE NONE\n";
    }

    if (have_first_reciprocal_failure) {
        std::cout
            << "\nFIRST_RECIPROCAL_FAILURE\n";

        print_pair(
            first_reciprocal_failure_prime,
            first_reciprocal_failure_pair.front(),
            first_reciprocal_failure_pair.back()
        );
    } else {
        std::cout
            << "\nFIRST_RECIPROCAL_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}