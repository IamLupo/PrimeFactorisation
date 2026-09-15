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

static u64 abs_diff(
    u64 a,
    u64 b
) {
    return a >= b ? a - b : b - a;
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

static u64 first_divisor_between(
    u64 value,
    u64 low,
    u64 high
) {
    if (high <= low + 1) {
        return 0;
    }

    for (
        u64 d = low + 1;
        d < high;
        ++d
    ) {
        if (value % d == 0) {
            return d;
        }
    }

    return 0;
}

static u64 last_divisor_between(
    u64 value,
    u64 low,
    u64 high
) {
    if (high <= low + 1) {
        return 0;
    }

    for (
        u64 d = high - 1;
        d > low;
        --d
    ) {
        if (value % d == 0) {
            return d;
        }
    }

    return 0;
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
    constexpr int EXPERIMENT = 424;
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

    u64 adjacent_on_left_value = 0;
    u64 adjacent_on_right_value = 0;

    u64 adjacent_on_both_values = 0;

    u64 intermediate_left_divisor = 0;
    u64 intermediate_right_divisor = 0;

    u64 intermediate_on_either = 0;
    u64 intermediate_on_both = 0;

    u64 same_value_adjacent = 0;
    u64 same_value_non_adjacent = 0;

    u64 difference_two_adjacent_left = 0;
    u64 difference_two_adjacent_right = 0;
    u64 difference_two_adjacent_both = 0;

    u64 sign_pp = 0;
    u64 sign_pm = 0;
    u64 sign_mp = 0;
    u64 sign_mm = 0;

    std::map<u64, u64>
        left_divisor_count_histogram;

    std::map<u64, u64>
        right_divisor_count_histogram;

    std::map<u64, u64>
        k_gap_histogram;

    Excursion first_non_adjacent;
    bool have_first_non_adjacent = false;

    Excursion first_same_value_non_adjacent;
    bool have_first_same_value_non_adjacent = false;

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
                end >= transition_states.size()
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

                    const u64 low_k = left.k;
                    const u64 high_k = right.k;

                    const u64 left_count =
                        count_divisors_between(
                            left_value,
                            low_k,
                            high_k
                        );

                    const u64 right_count =
                        count_divisors_between(
                            right_value,
                            low_k,
                            high_k
                        );

                    left_divisor_count_histogram[
                        left_count
                    ]++;

                    right_divisor_count_histogram[
                        right_count
                    ]++;

                    if (left_count == 0) {
                        ++adjacent_on_left_value;
                    } else {
                        ++intermediate_left_divisor;
                    }

                    if (right_count == 0) {
                        ++adjacent_on_right_value;
                    } else {
                        ++intermediate_right_divisor;
                    }

                    if (
                        left_count == 0 &&
                        right_count == 0
                    ) {
                        ++adjacent_on_both_values;
                    }

                    if (
                        left_count > 0 ||
                        right_count > 0
                    ) {
                        ++intermediate_on_either;

                        if (
                            left_count > 0 &&
                            right_count > 0
                        ) {
                            ++intermediate_on_both;
                        }

                        if (!have_first_non_adjacent) {
                            first_non_adjacent = ex;
                            have_first_non_adjacent = true;
                        }
                    }

                    const u64 delta =
                        abs_diff(
                            left_value,
                            right_value
                        );

                    if (delta == 0) {
                        ++same_value_pairs;

                        if (
                            left_count == 0 &&
                            right_count == 0
                        ) {
                            ++same_value_adjacent;
                        } else {
                            ++same_value_non_adjacent;

                            if (
                                !have_first_same_value_non_adjacent
                            ) {
                                first_same_value_non_adjacent =
                                    ex;
                                have_first_same_value_non_adjacent =
                                    true;
                            }
                        }
                    } else if (delta == 2) {
                        ++difference_two_pairs;

                        if (left_count == 0) {
                            ++difference_two_adjacent_left;
                        }

                        if (right_count == 0) {
                            ++difference_two_adjacent_right;
                        }

                        if (
                            left_count == 0 &&
                            right_count == 0
                        ) {
                            ++difference_two_adjacent_both;
                        }
                    }

                    if (
                        left.sign > 0 &&
                        right.sign > 0
                    ) {
                        ++sign_pp;
                    } else if (
                        left.sign > 0 &&
                        right.sign < 0
                    ) {
                        ++sign_pm;
                    } else if (
                        left.sign < 0 &&
                        right.sign > 0
                    ) {
                        ++sign_mp;
                    } else {
                        ++sign_mm;
                    }

                    k_gap_histogram[
                        right.k - left.k
                    ]++;

                    /*
                     * Keep this call because it verifies
                     * that the interval really contains no
                     * hidden divisor.
                     */
                    if (
                        left_count > 0 ||
                        right_count > 0
                    ) {
                        const u64 dl =
                            first_divisor_between(
                                left_value,
                                low_k,
                                high_k
                            );

                        const u64 dr =
                            last_divisor_between(
                                right_value,
                                low_k,
                                high_k
                            );

                        (void)dl;
                        (void)dr;
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
        << "\nADJACENCY\n";

    std::cout
        << "ADJACENT_ON_LEFT_VALUE="
        << adjacent_on_left_value
        << "\n";

    std::cout
        << "ADJACENT_ON_RIGHT_VALUE="
        << adjacent_on_right_value
        << "\n";

    std::cout
        << "ADJACENT_ON_BOTH_VALUES="
        << adjacent_on_both_values
        << "\n";

    std::cout
        << "INTERMEDIATE_LEFT_DIVISOR="
        << intermediate_left_divisor
        << "\n";

    std::cout
        << "INTERMEDIATE_RIGHT_DIVISOR="
        << intermediate_right_divisor
        << "\n";

    std::cout
        << "INTERMEDIATE_ON_EITHER="
        << intermediate_on_either
        << "\n";

    std::cout
        << "INTERMEDIATE_ON_BOTH="
        << intermediate_on_both
        << "\n";

    std::cout
        << "\nSAME_VALUE_ADJACENCY\n";

    std::cout
        << "SAME_VALUE_ADJACENT="
        << same_value_adjacent
        << "\n";

    std::cout
        << "SAME_VALUE_NON_ADJACENT="
        << same_value_non_adjacent
        << "\n";

    std::cout
        << "\nDIFFERENCE_TWO_ADJACENCY\n";

    std::cout
        << "DIFFERENCE_TWO_ADJACENT_LEFT="
        << difference_two_adjacent_left
        << "\n";

    std::cout
        << "DIFFERENCE_TWO_ADJACENT_RIGHT="
        << difference_two_adjacent_right
        << "\n";

    std::cout
        << "DIFFERENCE_TWO_ADJACENT_BOTH="
        << difference_two_adjacent_both
        << "\n";

    std::cout
        << "\nSIGN_PAIR_COUNTS\n";

    std::cout
        << "PP="
        << sign_pp
        << "\n";

    std::cout
        << "PM="
        << sign_pm
        << "\n";

    std::cout
        << "MP="
        << sign_mp
        << "\n";

    std::cout
        << "MM="
        << sign_mm
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
        << "\nLEFT_DIVISOR_COUNT_HISTOGRAM\n";

    for (
        const auto& entry :
        left_divisor_count_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nRIGHT_DIVISOR_COUNT_HISTOGRAM\n";

    for (
        const auto& entry :
        right_divisor_count_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_non_adjacent) {
        std::cout
            << "\nFIRST_NON_ADJACENT\n";

        print_excursion(
            first_non_adjacent
        );
    } else {
        std::cout
            << "\nFIRST_NON_ADJACENT NONE\n";
    }

    if (have_first_same_value_non_adjacent) {
        std::cout
            << "\nFIRST_SAME_VALUE_NON_ADJACENT\n";

        print_excursion(
            first_same_value_non_adjacent
        );
    } else {
        std::cout
            << "\nFIRST_SAME_VALUE_NON_ADJACENT NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
