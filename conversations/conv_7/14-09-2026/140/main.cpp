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

struct Transition {
    u64 K = 0;
    Witness old_state;
    Witness new_state;
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

static u64 safe_gcd(
    u64 a,
    u64 b
) {
    return std::gcd(a, b);
}

static bool divisible_u64(
    u64 a,
    u64 b
) {
    return b != 0 && a % b == 0;
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
    u64 r,
    const std::vector<Witness>& states
) {
    std::cout
        << "PRIME="
        << r
        << " LENGTH="
        << states.size()
        << "\n";

    for (
        std::size_t i = 0;
        i < states.size();
        ++i
    ) {
        std::cout
            << "  STATE["
            << i
            << "] ";

        print_witness(
            r,
            "W",
            states[i]
        );

        std::cout
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 432;
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

    u64 gcd_one = 0;
    u64 gcd_gt_one = 0;

    u64 product_less_A = 0;
    u64 product_equal_A = 0;
    u64 product_greater_A = 0;

    u64 product_divides_A = 0;
    u64 A_divides_product = 0;

    u64 k_left_divides_A = 0;
    u64 k_right_divides_A = 0;

    u64 k_left_divides_t_right = 0;
    u64 k_right_divides_t_left = 0;

    u64 t_left_divides_k_right = 0;
    u64 t_right_divides_k_left = 0;

    u64 coprime_reduced_pair = 0;

    u64 same_reduced_pair = 0;

    u64 equal_cross_products = 0;

    u64 product_gap_zero = 0;
    u64 product_gap_two = 0;
    u64 product_gap_other = 0;

    u64 symmetric_divisor_ratio = 0;

    std::map<u64, u64>
        gcd_histogram;

    std::map<u64, u64>
        normalized_left_histogram;

    std::map<u64, u64>
        normalized_right_histogram;

    std::map<u64, u64>
        product_ratio_numerator_histogram;

    std::map<u64, u64>
        product_ratio_denominator_histogram;

    std::vector<Witness> first_same_stream_pair;
    std::vector<Witness> first_sign_switch_pair;

    u64 first_same_stream_prime = 0;
    u64 first_sign_switch_prime = 0;

    bool have_first_same_stream = false;
    bool have_first_sign_switch = false;

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
                transition_states.push_back(
                    new_state
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
                start >= transition_states.size()
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

                    const bool same_stream =
                        left.m == right.m &&
                        left.sign == right.sign;

                    if (same_stream) {
                        ++same_stream_pairs;

                        if (
                            !have_first_same_stream
                        ) {
                            first_same_stream_pair =
                                excursion;
                            first_same_stream_prime =
                                r;
                            have_first_same_stream =
                                true;
                        }
                    } else {
                        ++sign_switch_pairs;

                        if (
                            !have_first_sign_switch
                        ) {
                            first_sign_switch_pair =
                                excursion;
                            first_sign_switch_prime =
                                r;
                            have_first_sign_switch =
                                true;
                        }
                    }

                    /*
                     * For same-stream pairs:
                     *
                     * A = mr + sign
                     * kL*tL = A = kR*tR
                     */
                    if (!same_stream) {
                        continue;
                    }

                    const u64 A =
                        static_cast<u64>(
                            arithmetic_value(
                                r,
                                left
                            )
                        );

                    const u64 g =
                        safe_gcd(
                            left.k,
                            right.k
                        );

                    gcd_histogram[g]++;

                    if (g == 1) {
                        ++gcd_one;
                    } else {
                        ++gcd_gt_one;
                    }

                    const i128 product =
                        static_cast<i128>(
                            left.k
                        ) *
                        static_cast<i128>(
                            right.k
                        );

                    if (
                        product <
                        static_cast<i128>(A)
                    ) {
                        ++product_less_A;
                    } else if (
                        product ==
                        static_cast<i128>(A)
                    ) {
                        ++product_equal_A;
                    } else {
                        ++product_greater_A;
                    }

                    if (
                        product != 0 &&
                        static_cast<i128>(A) %
                            product ==
                        0
                    ) {
                        ++product_divides_A;
                    }

                    if (
                        product != 0 &&
                        product %
                            static_cast<i128>(A) ==
                        0
                    ) {
                        ++A_divides_product;
                    }

                    if (
                        A % left.k == 0
                    ) {
                        ++k_left_divides_A;
                    }

                    if (
                        A % right.k == 0
                    ) {
                        ++k_right_divides_A;
                    }

                    if (
                        left.t % right.k == 0
                    ) {
                        ++k_right_divides_t_left;
                    }

                    if (
                        right.t % left.k == 0
                    ) {
                        ++k_left_divides_t_right;
                    }

                    if (
                        right.k % left.t == 0
                    ) {
                        ++t_left_divides_k_right;
                    }

                    if (
                        left.k % right.t == 0
                    ) {
                        ++t_right_divides_k_left;
                    }

                    /*
                     * Reduced factor pair.
                     */
                    const u64 left_reduced =
                        left.k / g;

                    const u64 right_reduced =
                        right.k / g;

                    normalized_left_histogram[
                        left_reduced
                    ]++;

                    normalized_right_histogram[
                        right_reduced
                    ]++;

                    if (
                        safe_gcd(
                            left_reduced,
                            right_reduced
                        ) == 1
                    ) {
                        ++coprime_reduced_pair;
                    }

                    if (
                        left_reduced ==
                        right_reduced
                    ) {
                        ++same_reduced_pair;
                    }

                    /*
                     * Cross-product relation:
                     *
                     * kL*tR versus kR*tL.
                     */
                    const i128 cross_left =
                        static_cast<i128>(
                            left.k
                        ) *
                        static_cast<i128>(
                            right.t
                        );

                    const i128 cross_right =
                        static_cast<i128>(
                            right.k
                        ) *
                        static_cast<i128>(
                            left.t
                        );

                    if (
                        cross_left ==
                        cross_right
                    ) {
                        ++equal_cross_products;
                    }

                    /*
                     * Product ratio:
                     *
                     * (kL*kR) / A
                     */
                    if (
                        product != 0
                    ) {
                        const u64 numerator =
                            static_cast<u64>(
                                product
                            );

                        const u64 denominator =
                            A;

                        const u64 rg =
                            safe_gcd(
                                numerator,
                                denominator
                            );

                        product_ratio_numerator_histogram[
                            numerator / rg
                        ]++;

                        product_ratio_denominator_histogram[
                            denominator / rg
                        ]++;
                    }

                    /*
                     * Since both divide A, the natural
                     * complementary identities can be tested.
                     */
                    if (
                        left.k *
                            left.t ==
                        right.k *
                            right.t
                    ) {
                        ++symmetric_divisor_ratio;
                    }

                    /*
                     * Difference in arithmetic products.
                     * Same-stream should be exactly zero.
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

                    const i128 gap =
                        left_product >=
                        right_product
                            ? left_product -
                                right_product
                            : right_product -
                                left_product;

                    if (gap == 0) {
                        ++product_gap_zero;
                    } else if (gap == 2) {
                        ++product_gap_two;
                    } else {
                        ++product_gap_other;
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
        << "\nGCD\n";

    std::cout
        << "GCD_ONE="
        << gcd_one
        << "\n";

    std::cout
        << "GCD_GT_ONE="
        << gcd_gt_one
        << "\n";

    std::cout
        << "\nPRODUCT_VS_A\n";

    std::cout
        << "PRODUCT_LESS_A="
        << product_less_A
        << "\n";

    std::cout
        << "PRODUCT_EQUAL_A="
        << product_equal_A
        << "\n";

    std::cout
        << "PRODUCT_GREATER_A="
        << product_greater_A
        << "\n";

    std::cout
        << "PRODUCT_DIVIDES_A="
        << product_divides_A
        << "\n";

    std::cout
        << "A_DIVIDES_PRODUCT="
        << A_divides_product
        << "\n";

    std::cout
        << "\nDIVISIBILITY\n";

    std::cout
        << "K_LEFT_DIVIDES_A="
        << k_left_divides_A
        << "\n";

    std::cout
        << "K_RIGHT_DIVIDES_A="
        << k_right_divides_A
        << "\n";

    std::cout
        << "K_RIGHT_DIVIDES_T_LEFT="
        << k_right_divides_t_left
        << "\n";

    std::cout
        << "K_LEFT_DIVIDES_T_RIGHT="
        << k_left_divides_t_right
        << "\n";

    std::cout
        << "T_LEFT_DIVIDES_K_RIGHT="
        << t_left_divides_k_right
        << "\n";

    std::cout
        << "T_RIGHT_DIVIDES_K_LEFT="
        << t_right_divides_k_left
        << "\n";

    std::cout
        << "\nREDUCED_PAIR\n";

    std::cout
        << "COPRIME_REDUCED_PAIR="
        << coprime_reduced_pair
        << "\n";

    std::cout
        << "SAME_REDUCED_PAIR="
        << same_reduced_pair
        << "\n";

    std::cout
        << "\nCROSS_PRODUCTS\n";

    std::cout
        << "EQUAL_CROSS_PRODUCTS="
        << equal_cross_products
        << "\n";

    std::cout
        << "SYMMETRIC_DIVISOR_RATIO="
        << symmetric_divisor_ratio
        << "\n";

    std::cout
        << "\nPRODUCT_GAP\n";

    std::cout
        << "PRODUCT_GAP_ZERO="
        << product_gap_zero
        << "\n";

    std::cout
        << "PRODUCT_GAP_TWO="
        << product_gap_two
        << "\n";

    std::cout
        << "PRODUCT_GAP_OTHER="
        << product_gap_other
        << "\n";

    std::cout
        << "\nGCD_HISTOGRAM\n";

    for (
        const auto& entry :
        gcd_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nREDUCED_LEFT_HISTOGRAM\n";

    for (
        const auto& entry :
        normalized_left_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nREDUCED_RIGHT_HISTOGRAM\n";

    for (
        const auto& entry :
        normalized_right_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_same_stream) {
        std::cout
            << "\nFIRST_SAME_STREAM_PAIR\n";

        print_excursion(
            first_same_stream_prime,
            first_same_stream_pair
        );
    } else {
        std::cout
            << "\nFIRST_SAME_STREAM_PAIR NONE\n";
    }

    if (have_first_sign_switch) {
        std::cout
            << "\nFIRST_SIGN_SWITCH_PAIR\n";

        print_excursion(
            first_sign_switch_prime,
            first_sign_switch_pair
        );
    } else {
        std::cout
            << "\nFIRST_SIGN_SWITCH_PAIR NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
