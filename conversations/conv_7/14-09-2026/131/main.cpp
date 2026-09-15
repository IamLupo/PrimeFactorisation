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

static bool lower_side(
    const Witness& w,
    i128 value
) {
    const i128 k =
        static_cast<i128>(w.k);

    return k * k <= value;
}

static bool upper_side(
    const Witness& w,
    i128 value
) {
    const i128 k =
        static_cast<i128>(w.k);

    return k * k >= value;
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
    constexpr int EXPERIMENT = 423;
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

    u64 projected_palindromes = 0;

    u64 k_strictly_increases = 0;
    u64 k_not_strictly_increases = 0;

    u64 t_strictly_decreases = 0;
    u64 t_not_strictly_decreases = 0;

    u64 opposite_k_order = 0;

    u64 left_lower_right_upper = 0;
    u64 left_lower_right_not_upper = 0;
    u64 left_not_lower_right_upper = 0;
    u64 neither_sqrt_side = 0;

    u64 left_strict_lower_right_strict_upper = 0;

    u64 left_equal_sqrt = 0;
    u64 right_equal_sqrt = 0;

    u64 same_value_pairs = 0;
    u64 difference_two_pairs = 0;
    u64 other_value_difference = 0;

    u64 same_value_k_order_ok = 0;
    u64 difference_two_k_order_ok = 0;

    u64 divisor_product_consistent = 0;
    u64 divisor_product_inconsistent = 0;

    std::map<u64, u64>
        k_gap_histogram;

    std::map<u64, u64>
        t_gap_histogram;

    Excursion first_k_failure;
    bool have_first_k_failure = false;

    Excursion first_t_failure;
    bool have_first_t_failure = false;

    Excursion first_sqrt_failure;
    bool have_first_sqrt_failure = false;

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

                ++projected_palindromes;

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

                    /*
                     * Arithmetic stream relation.
                     */
                    if (
                        left_value ==
                        right_value
                    ) {
                        ++same_value_pairs;
                    } else if (
                        abs_diff(
                            static_cast<u64>(
                                left_value
                            ),
                            static_cast<u64>(
                                right_value
                            )
                        ) == 2
                    ) {
                        ++difference_two_pairs;
                    } else {
                        ++other_value_difference;
                    }

                    /*
                     * k monotonicity.
                     */
                    if (left.k < right.k) {
                        ++k_strictly_increases;
                    } else {
                        ++k_not_strictly_increases;

                        if (!have_first_k_failure) {
                            first_k_failure = ex;
                            have_first_k_failure = true;
                        }
                    }

                    /*
                     * t monotonicity.
                     */
                    if (left.t > right.t) {
                        ++t_strictly_decreases;
                    } else {
                        ++t_not_strictly_decreases;

                        if (!have_first_t_failure) {
                            first_t_failure = ex;
                            have_first_t_failure = true;
                        }
                    }

                    if (
                        left.k < right.k &&
                        left.t > right.t
                    ) {
                        ++opposite_k_order;
                    }

                    /*
                     * Square-root location.
                     */
                    const bool left_lower =
                        lower_side(
                            left,
                            left_value
                        );

                    const bool right_upper =
                        upper_side(
                            right,
                            right_value
                        );

                    if (
                        left_lower &&
                        right_upper
                    ) {
                        ++left_lower_right_upper;

                        if (
                            left.k * left.k <
                            static_cast<u64>(
                                left_value
                            )
                        ) {
                            // Strict lower side.
                            if (
                                static_cast<i128>(
                                    right.k
                                ) *
                                static_cast<i128>(
                                    right.k
                                ) >
                                right_value
                            ) {
                                ++left_strict_lower_right_strict_upper;
                            }
                        }
                    } else if (
                        left_lower &&
                        !right_upper
                    ) {
                        ++left_lower_right_not_upper;
                    } else if (
                        !left_lower &&
                        right_upper
                    ) {
                        ++left_not_lower_right_upper;
                    } else {
                        ++neither_sqrt_side;

                        if (!have_first_sqrt_failure) {
                            first_sqrt_failure = ex;
                            have_first_sqrt_failure = true;
                        }
                    }

                    if (
                        static_cast<i128>(
                            left.k
                        ) *
                        static_cast<i128>(
                            left.k
                        ) ==
                        left_value
                    ) {
                        ++left_equal_sqrt;
                    }

                    if (
                        static_cast<i128>(
                            right.k
                        ) *
                        static_cast<i128>(
                            right.k
                        ) ==
                        right_value
                    ) {
                        ++right_equal_sqrt;
                    }

                    /*
                     * Product consistency.
                     *
                     * These are exact integer checks:
                     * k*t = mr + sign.
                     */
                    const i128 left_product =
                        static_cast<i128>(
                            left.k
                        ) *
                        static_cast<i128>(
                            left.t
                        );

                    const i128 right_product =
                        static_cast<i128>(
                            right.k
                        ) *
                        static_cast<i128>(
                            right.t
                        );

                    if (
                        left_product ==
                        left_value &&
                        right_product ==
                        right_value
                    ) {
                        ++divisor_product_consistent;
                    } else {
                        ++divisor_product_inconsistent;
                    }

                    /*
                     * Gap statistics.
                     */
                    k_gap_histogram[
                        right.k - left.k
                    ]++;

                    t_gap_histogram[
                        left.t - right.t
                    ]++;

                    /*
                     * Stream-specific k-order checks.
                     */
                    if (
                        left_value ==
                        right_value &&
                        left.k < right.k
                    ) {
                        ++same_value_k_order_ok;
                    }

                    if (
                        abs_diff(
                            static_cast<u64>(
                                left_value
                            ),
                            static_cast<u64>(
                                right_value
                            )
                        ) == 2 &&
                        left.k < right.k
                    ) {
                        ++difference_two_k_order_ok;
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
        << "PROJECTED_PALINDROMES="
        << projected_palindromes
        << "\n";

    std::cout
        << "TOTAL_MIRRORED_PAIRS="
        << total_mirrored_pairs
        << "\n";

    std::cout
        << "K_STRICTLY_INCREASES="
        << k_strictly_increases
        << "\n";

    std::cout
        << "K_NOT_STRICTLY_INCREASES="
        << k_not_strictly_increases
        << "\n";

    std::cout
        << "T_STRICTLY_DECREASES="
        << t_strictly_decreases
        << "\n";

    std::cout
        << "T_NOT_STRICTLY_DECREASES="
        << t_not_strictly_decreases
        << "\n";

    std::cout
        << "K_INCREASE_T_DECREASE="
        << opposite_k_order
        << "\n";

    std::cout
        << "\nSQRT_SIDE_CLASSIFICATION\n";

    std::cout
        << "LEFT_LOWER_RIGHT_UPPER="
        << left_lower_right_upper
        << "\n";

    std::cout
        << "LEFT_LOWER_RIGHT_NOT_UPPER="
        << left_lower_right_not_upper
        << "\n";

    std::cout
        << "LEFT_NOT_LOWER_RIGHT_UPPER="
        << left_not_lower_right_upper
        << "\n";

    std::cout
        << "NEITHER_SQRT_SIDE="
        << neither_sqrt_side
        << "\n";

    std::cout
        << "LEFT_STRICT_LOWER_RIGHT_STRICT_UPPER="
        << left_strict_lower_right_strict_upper
        << "\n";

    std::cout
        << "LEFT_EQUAL_SQRT="
        << left_equal_sqrt
        << "\n";

    std::cout
        << "RIGHT_EQUAL_SQRT="
        << right_equal_sqrt
        << "\n";

    std::cout
        << "\nARITHMETIC_RELATION\n";

    std::cout
        << "SAME_VALUE_PAIRS="
        << same_value_pairs
        << "\n";

    std::cout
        << "DIFFERENCE_TWO_PAIRS="
        << difference_two_pairs
        << "\n";

    std::cout
        << "OTHER_VALUE_DIFFERENCE="
        << other_value_difference
        << "\n";

    std::cout
        << "SAME_VALUE_K_ORDER_OK="
        << same_value_k_order_ok
        << "\n";

    std::cout
        << "DIFFERENCE_TWO_K_ORDER_OK="
        << difference_two_k_order_ok
        << "\n";

    std::cout
        << "\nDIVISOR_PRODUCT_CHECK\n";

    std::cout
        << "DIVISOR_PRODUCT_CONSISTENT="
        << divisor_product_consistent
        << "\n";

    std::cout
        << "DIVISOR_PRODUCT_INCONSISTENT="
        << divisor_product_inconsistent
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

    if (have_first_k_failure) {
        std::cout
            << "\nFIRST_K_ORDER_FAILURE\n";

        print_excursion(
            first_k_failure
        );
    } else {
        std::cout
            << "\nFIRST_K_ORDER_FAILURE NONE\n";
    }

    if (have_first_t_failure) {
        std::cout
            << "\nFIRST_T_ORDER_FAILURE\n";

        print_excursion(
            first_t_failure
        );
    } else {
        std::cout
            << "\nFIRST_T_ORDER_FAILURE NONE\n";
    }

    if (have_first_sqrt_failure) {
        std::cout
            << "\nFIRST_SQRT_SIDE_FAILURE\n";

        print_excursion(
            first_sqrt_failure
        );
    } else {
        std::cout
            << "\nFIRST_SQRT_SIDE_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
