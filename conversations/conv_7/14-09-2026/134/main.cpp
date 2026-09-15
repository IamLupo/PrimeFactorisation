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

static u64 count_divisors_between(
    u64 value,
    u64 low,
    u64 high
) {
    if (high <= low + 1) {
        return 0;
    }

    u64 count = 0;

    for (
        u64 d = low + 1;
        d < high;
        ++d
    ) {
        if (value % d == 0) {
            ++count;
        }
    }

    return count;
}

static u64 count_divisors_up_to(
    u64 value,
    u64 high
) {
    u64 count = 0;

    for (u64 d = 1; d <= high; ++d) {
        if (value % d == 0) {
            ++count;
        }
    }

    return count;
}

static u64 divisor_rank(
    u64 value,
    u64 k
) {
    if (
        k == 0 ||
        value % k != 0
    ) {
        return 0;
    }

    u64 rank = 0;

    for (u64 d = 1; d <= k; ++d) {
        if (value % d == 0) {
            ++rank;
        }
    }

    return rank;
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

int main() {
    constexpr int EXPERIMENT = 426;
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

    u64 same_value_pairs = 0;
    u64 difference_two_pairs = 0;

    u64 same_hidden_divisor_zero = 0;
    u64 same_hidden_divisor_nonzero = 0;

    u64 diff_hidden_left_zero = 0;
    u64 diff_hidden_left_nonzero = 0;

    u64 diff_hidden_right_zero = 0;
    u64 diff_hidden_right_nonzero = 0;

    u64 diff_hidden_both_zero = 0;

    u64 total_hidden_divisors_same = 0;
    u64 total_hidden_divisors_diff_left = 0;
    u64 total_hidden_divisors_diff_right = 0;

    u64 same_rank_gap_sum = 0;
    u64 same_visible_transition_gap_sum = 0;

    u64 same_rank_gap_equals_hidden_plus_one = 0;

    u64 same_hidden_equals_visible_states = 0;
    u64 same_hidden_less_than_visible_states = 0;
    u64 same_hidden_greater_than_visible_states = 0;

    u64 diff_left_rank_gap_sum = 0;
    u64 diff_right_rank_gap_sum = 0;

    std::map<u64, u64>
        same_hidden_divisor_histogram;

    std::map<u64, u64>
        diff_left_hidden_divisor_histogram;

    std::map<u64, u64>
        diff_right_hidden_divisor_histogram;

    std::map<u64, u64>
        same_rank_gap_histogram;

    Excursion first_same_hidden;
    bool have_first_same_hidden = false;

    Excursion first_diff_hidden;
    bool have_first_diff_hidden = false;

    Excursion first_rank_mismatch;
    bool have_first_rank_mismatch = false;

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

        std::vector<Witness> transition_states;
        transition_states.reserve(128);

        Witness previous =
            global_winner(best_by_m);

        transition_states.push_back(
            previous
        );

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
                    best_by_m[mi] =
                        make_witness(
                            r,
                            m,
                            K,
                            +1
                        );
                } else if (
                    divides_stream(
                        r,
                        m,
                        K,
                        -1
                    )
                ) {
                    best_by_m[mi] =
                        make_witness(
                            r,
                            m,
                            K,
                            -1
                        );
                }
            }

            const Witness current =
                global_winner(best_by_m);

            if (current.m != previous.m) {
                transition_states.push_back(
                    current
                );
            }

            previous = current;
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

                const std::size_t n =
                    ex.states.size();

                for (
                    std::size_t i = 0;
                    i < n / 2;
                    ++i
                ) {
                    const std::size_t j =
                        n - 1 - i;

                    const Witness& left =
                        ex.states[i];

                    const Witness& right =
                        ex.states[j];

                    ++total_mirrored_pairs;

                    const u64 left_value =
                        static_cast<u64>(
                            arithmetic_value(
                                r,
                                left
                            )
                        );

                    const u64 right_value =
                        static_cast<u64>(
                            arithmetic_value(
                                r,
                                right
                            )
                        );

                    const u64 delta_value =
                        left_value >= right_value
                            ? left_value - right_value
                            : right_value - left_value;

                    /*
                     * Number of visible global transitions
                     * strictly between the mirrored states.
                     */
                    const u64 visible_gap =
                        static_cast<u64>(
                            j - i - 1
                        );

                    if (delta_value == 0) {
                        ++same_value_pairs;

                        const u64 hidden =
                            count_divisors_between(
                                left_value,
                                left.k,
                                right.k
                            );

                        same_hidden_divisor_histogram[
                            hidden
                        ]++;

                        total_hidden_divisors_same +=
                            hidden;

                        if (hidden == 0) {
                            ++same_hidden_divisor_zero;
                        } else {
                            ++same_hidden_divisor_nonzero;

                            if (!have_first_same_hidden) {
                                first_same_hidden = ex;
                                have_first_same_hidden = true;
                            }
                        }

                        const u64 rank_left =
                            divisor_rank(
                                left_value,
                                left.k
                            );

                        const u64 rank_right =
                            divisor_rank(
                                right_value,
                                right.k
                            );

                        if (rank_right >= rank_left) {
                            const u64 rank_gap =
                                rank_right -
                                rank_left;

                            same_rank_gap_sum +=
                                rank_gap;

                            same_rank_gap_histogram[
                                rank_gap
                            ]++;

                            if (
                                rank_gap ==
                                hidden + 1
                            ) {
                                ++same_rank_gap_equals_hidden_plus_one;
                            }

                            if (
                                rank_gap <
                                visible_gap
                            ) {
                                ++same_hidden_less_than_visible_states;
                            } else if (
                                rank_gap ==
                                visible_gap
                            ) {
                                ++same_hidden_equals_visible_states;
                            } else {
                                ++same_hidden_greater_than_visible_states;
                            }

                            same_visible_transition_gap_sum +=
                                visible_gap;
                        }

                        if (
                            hidden + 1 !=
                            visible_gap + 1
                        ) {
                            if (!have_first_rank_mismatch) {
                                first_rank_mismatch = ex;
                                have_first_rank_mismatch = true;
                            }
                        }
                    } else if (delta_value == 2) {
                        ++difference_two_pairs;

                        const u64 hidden_left =
                            count_divisors_between(
                                left_value,
                                left.k,
                                right.k
                            );

                        const u64 hidden_right =
                            count_divisors_between(
                                right_value,
                                left.k,
                                right.k
                            );

                        diff_left_hidden_divisor_histogram[
                            hidden_left
                        ]++;

                        diff_right_hidden_divisor_histogram[
                            hidden_right
                        ]++;

                        total_hidden_divisors_diff_left +=
                            hidden_left;

                        total_hidden_divisors_diff_right +=
                            hidden_right;

                        if (hidden_left == 0) {
                            ++diff_hidden_left_zero;
                        } else {
                            ++diff_hidden_left_nonzero;
                        }

                        if (hidden_right == 0) {
                            ++diff_hidden_right_zero;
                        } else {
                            ++diff_hidden_right_nonzero;
                        }

                        if (
                            hidden_left == 0 &&
                            hidden_right == 0
                        ) {
                            ++diff_hidden_both_zero;
                        }

                        if (!have_first_diff_hidden) {
                            if (
                                hidden_left > 0 ||
                                hidden_right > 0
                            ) {
                                first_diff_hidden = ex;
                                have_first_diff_hidden = true;
                            }
                        }

                        const u64 rank_left =
                            divisor_rank(
                                left_value,
                                left.k
                            );

                        const u64 rank_right =
                            divisor_rank(
                                right_value,
                                right.k
                            );

                        if (rank_left > 0) {
                            diff_left_rank_gap_sum +=
                                rank_left;
                        }

                        if (rank_right > 0) {
                            diff_right_rank_gap_sum +=
                                rank_right;
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
        << "\nSAME_VALUE_HIDDEN_DIVISORS\n";

    std::cout
        << "HIDDEN_ZERO="
        << same_hidden_divisor_zero
        << "\n";

    std::cout
        << "HIDDEN_NONZERO="
        << same_hidden_divisor_nonzero
        << "\n";

    std::cout
        << "TOTAL_HIDDEN="
        << total_hidden_divisors_same
        << "\n";

    std::cout
        << "SAME_RANK_GAP_SUM="
        << same_rank_gap_sum
        << "\n";

    std::cout
        << "SAME_VISIBLE_TRANSITION_GAP_SUM="
        << same_visible_transition_gap_sum
        << "\n";

    std::cout
        << "RANK_GAP_EQUALS_HIDDEN_PLUS_ONE="
        << same_rank_gap_equals_hidden_plus_one
        << "\n";

    std::cout
        << "HIDDEN_LESS_THAN_VISIBLE="
        << same_hidden_less_than_visible_states
        << "\n";

    std::cout
        << "HIDDEN_EQUALS_VISIBLE="
        << same_hidden_equals_visible_states
        << "\n";

    std::cout
        << "HIDDEN_GREATER_THAN_VISIBLE="
        << same_hidden_greater_than_visible_states
        << "\n";

    std::cout
        << "\nDIFFERENCE_TWO_HIDDEN_DIVISORS\n";

    std::cout
        << "LEFT_ZERO="
        << diff_hidden_left_zero
        << "\n";

    std::cout
        << "LEFT_NONZERO="
        << diff_hidden_left_nonzero
        << "\n";

    std::cout
        << "RIGHT_ZERO="
        << diff_hidden_right_zero
        << "\n";

    std::cout
        << "RIGHT_NONZERO="
        << diff_hidden_right_nonzero
        << "\n";

    std::cout
        << "BOTH_ZERO="
        << diff_hidden_both_zero
        << "\n";

    std::cout
        << "TOTAL_LEFT_HIDDEN="
        << total_hidden_divisors_diff_left
        << "\n";

    std::cout
        << "TOTAL_RIGHT_HIDDEN="
        << total_hidden_divisors_diff_right
        << "\n";

    std::cout
        << "LEFT_RANK_GAP_SUM="
        << diff_left_rank_gap_sum
        << "\n";

    std::cout
        << "RIGHT_RANK_GAP_SUM="
        << diff_right_rank_gap_sum
        << "\n";

    std::cout
        << "\nSAME_HIDDEN_DIVISOR_HISTOGRAM\n";

    for (
        const auto& entry :
        same_hidden_divisor_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nDIFF_LEFT_HIDDEN_DIVISOR_HISTOGRAM\n";

    for (
        const auto& entry :
        diff_left_hidden_divisor_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nDIFF_RIGHT_HIDDEN_DIVISOR_HISTOGRAM\n";

    for (
        const auto& entry :
        diff_right_hidden_divisor_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nSAME_RANK_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        same_rank_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_same_hidden) {
        std::cout
            << "\nFIRST_SAME_VALUE_WITH_HIDDEN_DIVISOR\n";

        print_excursion(
            first_same_hidden
        );
    } else {
        std::cout
            << "\nFIRST_SAME_VALUE_WITH_HIDDEN_DIVISOR NONE\n";
    }

    if (have_first_diff_hidden) {
        std::cout
            << "\nFIRST_DIFFERENCE_TWO_WITH_HIDDEN_DIVISOR\n";

        print_excursion(
            first_diff_hidden
        );
    } else {
        std::cout
            << "\nFIRST_DIFFERENCE_TWO_WITH_HIDDEN_DIVISOR NONE\n";
    }

    if (have_first_rank_mismatch) {
        std::cout
            << "\nFIRST_RANK_MISMATCH\n";

        print_excursion(
            first_rank_mismatch
        );
    } else {
        std::cout
            << "\nFIRST_RANK_MISMATCH NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
