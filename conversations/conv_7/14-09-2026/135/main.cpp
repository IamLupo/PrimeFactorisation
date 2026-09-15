#include <cstdint>
#include <iostream>
#include <map>
#include <set>
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

struct HiddenOwner {
    u64 divisor = 0;
    int m = 0;
    int sign = 0;
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

static u64 count_visible_transitions(
    const std::vector<Witness>& states,
    std::size_t left,
    std::size_t right
) {
    if (right <= left + 1) {
        return 0;
    }

    return static_cast<u64>(
        right - left - 1
    );
}

static bool hidden_divisor_owner(
    u64 r,
    u64 divisor,
    int m_max,
    HiddenOwner& owner
) {
    for (int m = 1; m <= m_max; ++m) {
        if (
            divides_stream(
                r,
                m,
                divisor,
                +1
            )
        ) {
            owner.divisor = divisor;
            owner.m = m;
            owner.sign = +1;
            return true;
        }

        if (
            divides_stream(
                r,
                m,
                divisor,
                -1
            )
        ) {
            owner.divisor = divisor;
            owner.m = m;
            owner.sign = -1;
            return true;
        }
    }

    return false;
}

static bool is_global_winner_at_k(
    u64 r,
    u64 k,
    int m_max,
    int& winner_m,
    int& winner_sign
) {
    bool have = false;
    Witness best;

    for (int m = 1; m <= m_max; ++m) {
        if (
            divides_stream(
                r,
                m,
                k,
                +1
            )
        ) {
            const Witness w =
                make_witness(
                    r,
                    m,
                    k,
                    +1
                );

            if (
                !have ||
                better_normalized(
                    w,
                    best
                )
            ) {
                best = w;
                have = true;
            }
        }

        if (
            divides_stream(
                r,
                m,
                k,
                -1
            )
        ) {
            const Witness w =
                make_witness(
                    r,
                    m,
                    k,
                    -1
                );

            if (
                !have ||
                better_normalized(
                    w,
                    best
                )
            ) {
                best = w;
                have = true;
            }
        }
    }

    if (!have) {
        return false;
    }

    winner_m = best.m;
    winner_sign = best.sign;

    return true;
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
    constexpr int EXPERIMENT = 427;
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

    u64 total_hidden_divisors = 0;
    u64 hidden_divisors_owned_by_winner = 0;
    u64 hidden_divisors_owned_by_nonwinner = 0;

    u64 hidden_divisors_with_owner = 0;
    u64 hidden_divisors_without_owner = 0;

    u64 hidden_divisors_that_are_global_records = 0;
    u64 hidden_divisors_that_are_not_global_records = 0;

    u64 hidden_divisors_same_m_as_mirror = 0;
    u64 hidden_divisors_other_m = 0;

    std::map<int, u64>
        hidden_m_histogram;

    std::map<int, u64>
        hidden_sign_histogram;

    std::map<int, u64>
        hidden_m_same_as_endpoint_histogram;

    std::map<u64, u64>
        hidden_count_per_pair_histogram;

    std::map<u64, u64>
        hidden_record_count_per_pair_histogram;

    std::map<u64, u64>
        visible_transition_count_histogram;

    Excursion first_hidden_winner;
    bool have_first_hidden_winner = false;

    Excursion first_hidden_nonwinner;
    bool have_first_hidden_nonwinner = false;

    Excursion first_hidden_record;
    bool have_first_hidden_record = false;

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

                    const u64 difference =
                        left_value >= right_value
                            ? left_value - right_value
                            : right_value - left_value;

                    if (difference == 0) {
                        ++same_value_pairs;
                    } else if (difference == 2) {
                        ++difference_two_pairs;
                    }

                    const u64 hidden_low =
                        left.k + 1;

                    const u64 hidden_high =
                        right.k;

                    u64 pair_hidden_count = 0;
                    u64 pair_record_count = 0;

                    for (
                        u64 d = hidden_low;
                        d < hidden_high;
                        ++d
                    ) {
                        HiddenOwner owner;

                        if (
                            !hidden_divisor_owner(
                                r,
                                d,
                                M_MAX,
                                owner
                            )
                        ) {
                            ++hidden_divisors_without_owner;
                            continue;
                        }

                        ++hidden_divisors_with_owner;
                        ++total_hidden_divisors;
                        ++pair_hidden_count;

                        hidden_m_histogram[
                            owner.m
                        ]++;

                        hidden_sign_histogram[
                            owner.sign
                        ]++;

                        if (
                            owner.m == left.m ||
                            owner.m == right.m
                        ) {
                            ++hidden_divisors_same_m_as_mirror;

                            hidden_m_same_as_endpoint_histogram[
                                owner.m
                            ]++;
                        } else {
                            ++hidden_divisors_other_m;
                        }

                        /*
                         * Is this divisor itself the winner
                         * at its K-coordinate?
                         */
                        int winner_m = 0;
                        int winner_sign = 0;

                        const bool has_winner =
                            is_global_winner_at_k(
                                r,
                                d,
                                M_MAX,
                                winner_m,
                                winner_sign
                            );

                        if (
                            has_winner &&
                            winner_m == owner.m &&
                            winner_sign == owner.sign
                        ) {
                            ++hidden_divisors_that_are_global_records;
                            ++pair_record_count;

                            if (
                                !have_first_hidden_record
                            ) {
                                first_hidden_record = ex;
                                have_first_hidden_record = true;
                            }
                        } else {
                            ++hidden_divisors_that_are_not_global_records;
                        }

                        /*
                         * Was the hidden divisor itself
                         * owned by the same multiplier as
                         * the endpoint?
                         */
                        if (
                            owner.m == left.m ||
                            owner.m == right.m
                        ) {
                            /*
                             * No separate counter needed.
                             */
                        }

                        /*
                         * Is this hidden candidate actually
                         * better than the endpoint winner
                         * in normalized score?
                         */
                        const i128 left_score =
                            static_cast<i128>(left.k) *
                            static_cast<i128>(
                                owner.m
                            );

                        const i128 candidate_score =
                            static_cast<i128>(d) *
                            static_cast<i128>(left.m);

                        if (
                            candidate_score >
                            left_score
                        ) {
                            ++hidden_divisors_owned_by_winner;

                            if (
                                !have_first_hidden_winner
                            ) {
                                first_hidden_winner = ex;
                                have_first_hidden_winner = true;
                            }
                        } else {
                            ++hidden_divisors_owned_by_nonwinner;

                            if (
                                !have_first_hidden_nonwinner
                            ) {
                                first_hidden_nonwinner = ex;
                                have_first_hidden_nonwinner = true;
                            }
                        }
                    }

                    hidden_count_per_pair_histogram[
                        pair_hidden_count
                    ]++;

                    hidden_record_count_per_pair_histogram[
                        pair_record_count
                    ]++;

                    visible_transition_count_histogram[
                        count_visible_transitions(
                            ex.states,
                            i,
                            j
                        )
                    ];
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
        << "\nHIDDEN_DIVISOR_TOTALS\n";

    std::cout
        << "TOTAL_HIDDEN_DIVISORS="
        << total_hidden_divisors
        << "\n";

    std::cout
        << "HIDDEN_WITH_OWNER="
        << hidden_divisors_with_owner
        << "\n";

    std::cout
        << "HIDDEN_WITHOUT_OWNER="
        << hidden_divisors_without_owner
        << "\n";

    std::cout
        << "HIDDEN_SAME_M_AS_ENDPOINT="
        << hidden_divisors_same_m_as_mirror
        << "\n";

    std::cout
        << "HIDDEN_OTHER_M="
        << hidden_divisors_other_m
        << "\n";

    std::cout
        << "HIDDEN_GLOBAL_RECORDS="
        << hidden_divisors_that_are_global_records
        << "\n";

    std::cout
        << "HIDDEN_NOT_GLOBAL_RECORDS="
        << hidden_divisors_that_are_not_global_records
        << "\n";

    std::cout
        << "\nHIDDEN_M_HISTOGRAM\n";

    for (
        const auto& entry :
        hidden_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nHIDDEN_SIGN_HISTOGRAM\n";

    for (
        const auto& entry :
        hidden_sign_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nHIDDEN_COUNT_PER_PAIR\n";

    for (
        const auto& entry :
        hidden_count_per_pair_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nHIDDEN_RECORD_COUNT_PER_PAIR\n";

    for (
        const auto& entry :
        hidden_record_count_per_pair_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nVISIBLE_TRANSITION_COUNT_PER_PAIR\n";

    for (
        const auto& entry :
        visible_transition_count_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_hidden_winner) {
        std::cout
            << "\nFIRST_HIDDEN_DIVISOR_RECORD\n";

        print_excursion(
            first_hidden_winner
        );
    } else {
        std::cout
            << "\nFIRST_HIDDEN_DIVISOR_RECORD NONE\n";
    }

    if (have_first_hidden_nonwinner) {
        std::cout
            << "\nFIRST_HIDDEN_NONRECORD\n";

        print_excursion(
            first_hidden_nonwinner
        );
    } else {
        std::cout
            << "\nFIRST_HIDDEN_NONRECORD NONE\n";
    }

    if (have_first_hidden_record) {
        std::cout
            << "\nFIRST_HIDDEN_GLOBAL_RECORD\n";

        print_excursion(
            first_hidden_record
        );
    } else {
        std::cout
            << "\nFIRST_HIDDEN_GLOBAL_RECORD NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
