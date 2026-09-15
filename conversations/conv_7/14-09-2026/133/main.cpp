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

struct PairStats {
    u64 total = 0;

    u64 same_value = 0;
    u64 diff_two = 0;
    u64 other = 0;

    u64 same_below_above = 0;
    u64 same_below_below = 0;
    u64 same_above_above = 0;

    u64 diff_below_above = 0;
    u64 diff_below_below = 0;
    u64 diff_above_above = 0;

    u64 same_strict_straddle = 0;
    u64 diff_strict_straddle = 0;

    u64 same_k_increase = 0;
    u64 diff_k_increase = 0;

    u64 same_t_decrease = 0;
    u64 diff_t_decrease = 0;

    u64 same_exact_product = 0;
    u64 diff_product_plus_or_minus_two = 0;

    u64 same_ratio_identity = 0;

    u64 same_gap_identity = 0;
    u64 diff_gap_identity = 0;
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

static void classify_pair(
    u64 r,
    const Witness& left,
    const Witness& right,
    PairStats& stats
) {
    ++stats.total;

    const i128 lv128 =
        arithmetic_value(r, left);

    const i128 rv128 =
        arithmetic_value(r, right);

    const u64 lv =
        static_cast<u64>(lv128);

    const u64 rv =
        static_cast<u64>(rv128);

    const u64 low_value =
        lv < rv ? lv : rv;

    const u64 high_value =
        lv < rv ? rv : lv;

    const u64 difference =
        lv >= rv ? lv - rv : rv - lv;

    const bool left_below =
        static_cast<i128>(left.k) *
        static_cast<i128>(left.k) <
        lv128;

    const bool left_above =
        static_cast<i128>(left.k) *
        static_cast<i128>(left.k) >
        lv128;

    const bool right_below =
        static_cast<i128>(right.k) *
        static_cast<i128>(right.k) <
        rv128;

    const bool right_above =
        static_cast<i128>(right.k) *
        static_cast<i128>(right.k) >
        rv128;

    const bool same =
        difference == 0;

    const bool diff2 =
        difference == 2;

    if (same) {
        ++stats.same_value;

        if (left_below && right_above) {
            ++stats.same_below_above;
            ++stats.same_strict_straddle;
        } else if (left_below && right_below) {
            ++stats.same_below_below;
        } else if (left_above && right_above) {
            ++stats.same_above_above;
        }

        const i128 left_product =
            static_cast<i128>(left.k) *
            static_cast<i128>(left.t);

        const i128 right_product =
            static_cast<i128>(right.k) *
            static_cast<i128>(right.t);

        if (left_product == right_product) {
            ++stats.same_exact_product;
        }

        /*
         * For equal values:
         *
         * k_L t_L = k_R t_R
         *
         * hence
         *
         * k_R / k_L = t_L / t_R.
         */
        if (
            static_cast<i128>(right.k) *
            static_cast<i128>(right.t) ==
            static_cast<i128>(left.k) *
            static_cast<i128>(left.t)
        ) {
            ++stats.same_ratio_identity;
        }

        /*
         * Rearranged gap identity:
         *
         * Δk * t_L - Δt * k_L
         * - Δk * Δt = 0.
         */
        const i128 dk =
            static_cast<i128>(right.k) -
            static_cast<i128>(left.k);

        const i128 dt =
            static_cast<i128>(left.t) -
            static_cast<i128>(right.t);

        const i128 identity =
            dk * static_cast<i128>(left.t) -
            dt * static_cast<i128>(left.k) -
            dk * dt;

        if (identity == 0) {
            ++stats.same_gap_identity;
        }

        if (right.k > left.k) {
            ++stats.same_k_increase;
        }

        if (left.t > right.t) {
            ++stats.same_t_decrease;
        }
    } else if (diff2) {
        ++stats.diff_two;

        if (left_below && right_above) {
            ++stats.diff_below_above;
            ++stats.diff_strict_straddle;
        } else if (left_below && right_below) {
            ++stats.diff_below_below;
        } else if (left_above && right_above) {
            ++stats.diff_above_above;
        }

        const i128 left_product =
            static_cast<i128>(left.k) *
            static_cast<i128>(left.t);

        const i128 right_product =
            static_cast<i128>(right.k) *
            static_cast<i128>(right.t);

        const i128 product_difference =
            right_product - left_product;

        if (
            product_difference == 2 ||
            product_difference == -2
        ) {
            ++stats.diff_product_plus_or_minus_two;
        }

        /*
         * For difference two:
         *
         * Δk * t_L - Δt * k_L
         * - Δk * Δt = ±2.
         */
        const i128 dk =
            static_cast<i128>(right.k) -
            static_cast<i128>(left.k);

        const i128 dt =
            static_cast<i128>(left.t) -
            static_cast<i128>(right.t);

        const i128 identity =
            dk * static_cast<i128>(left.t) -
            dt * static_cast<i128>(left.k) -
            dk * dt;

        if (
            identity == 2 ||
            identity == -2
        ) {
            ++stats.diff_gap_identity;
        }

        if (right.k > left.k) {
            ++stats.diff_k_increase;
        }

        if (left.t > right.t) {
            ++stats.diff_t_decrease;
        }
    } else {
        ++stats.other;
    }

    (void)low_value;
    (void)high_value;
}

int main() {
    constexpr int EXPERIMENT = 425;
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

    PairStats stats;

    std::map<u64, u64>
        same_k_gap_histogram;

    std::map<u64, u64>
        diff_k_gap_histogram;

    Excursion first_same_non_straddle;
    bool have_first_same_non_straddle = false;

    Excursion first_diff_non_straddle;
    bool have_first_diff_non_straddle = false;

    Excursion first_same_identity_failure;
    bool have_first_same_identity_failure = false;

    Excursion first_diff_identity_failure;
    bool have_first_diff_identity_failure = false;

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

                    PairStats before = stats;

                    classify_pair(
                        r,
                        left,
                        right,
                        stats
                    );

                    const i128 lv =
                        arithmetic_value(
                            r,
                            left
                        );

                    const i128 rv =
                        arithmetic_value(
                            r,
                            right
                        );

                    const u64 difference =
                        lv >= rv
                            ? static_cast<u64>(
                                lv - rv
                              )
                            : static_cast<u64>(
                                rv - lv
                              );

                    const bool left_below =
                        static_cast<i128>(left.k) *
                        static_cast<i128>(left.k) <
                        lv;

                    const bool right_above =
                        static_cast<i128>(right.k) *
                        static_cast<i128>(right.k) >
                        rv;

                    if (
                        difference == 0
                    ) {
                        same_k_gap_histogram[
                            right.k - left.k
                        ]++;

                        if (
                            !(left_below &&
                              right_above)
                        ) {
                            if (
                                !have_first_same_non_straddle
                            ) {
                                first_same_non_straddle =
                                    ex;
                                have_first_same_non_straddle =
                                    true;
                            }
                        }

                        if (
                            stats.same_ratio_identity ==
                            before.same_ratio_identity
                        ) {
                            if (
                                !have_first_same_identity_failure
                            ) {
                                first_same_identity_failure =
                                    ex;
                                have_first_same_identity_failure =
                                    true;
                            }
                        }
                    } else if (
                        difference == 2
                    ) {
                        diff_k_gap_histogram[
                            right.k - left.k
                        ]++;

                        if (
                            !(left_below &&
                              right_above)
                        ) {
                            if (
                                !have_first_diff_non_straddle
                            ) {
                                first_diff_non_straddle =
                                    ex;
                                have_first_diff_non_straddle =
                                    true;
                            }
                        }

                        if (
                            stats.diff_gap_identity ==
                            before.diff_gap_identity
                        ) {
                            if (
                                !have_first_diff_identity_failure
                            ) {
                                first_diff_identity_failure =
                                    ex;
                                have_first_diff_identity_failure =
                                    true;
                            }
                        }
                    }
                }
            }

            start = end;
        }
    }

    std::cout
        << "\nTOTAL_MIRRORED_PAIRS="
        << stats.total
        << "\n";

    std::cout
        << "SAME_VALUE="
        << stats.same_value
        << "\n";

    std::cout
        << "DIFFERENCE_TWO="
        << stats.diff_two
        << "\n";

    std::cout
        << "OTHER="
        << stats.other
        << "\n";

    std::cout
        << "\nSAME_VALUE_SQRT_CLASSIFICATION\n";

    std::cout
        << "SAME_BELOW_ABOVE="
        << stats.same_below_above
        << "\n";

    std::cout
        << "SAME_BELOW_BELOW="
        << stats.same_below_below
        << "\n";

    std::cout
        << "SAME_ABOVE_ABOVE="
        << stats.same_above_above
        << "\n";

    std::cout
        << "SAME_STRICT_STRADDLE="
        << stats.same_strict_straddle
        << "\n";

    std::cout
        << "\nDIFFERENCE_TWO_SQRT_CLASSIFICATION\n";

    std::cout
        << "DIFF_BELOW_ABOVE="
        << stats.diff_below_above
        << "\n";

    std::cout
        << "DIFF_BELOW_BELOW="
        << stats.diff_below_below
        << "\n";

    std::cout
        << "DIFF_ABOVE_ABOVE="
        << stats.diff_above_above
        << "\n";

    std::cout
        << "DIFF_STRICT_STRADDLE="
        << stats.diff_strict_straddle
        << "\n";

    std::cout
        << "\nORDER_CHECKS\n";

    std::cout
        << "SAME_K_INCREASE="
        << stats.same_k_increase
        << "\n";

    std::cout
        << "DIFF_K_INCREASE="
        << stats.diff_k_increase
        << "\n";

    std::cout
        << "SAME_T_DECREASE="
        << stats.same_t_decrease
        << "\n";

    std::cout
        << "DIFF_T_DECREASE="
        << stats.diff_t_decrease
        << "\n";

    std::cout
        << "\nALGEBRAIC_IDENTITIES\n";

    std::cout
        << "SAME_EXACT_PRODUCT="
        << stats.same_exact_product
        << "\n";

    std::cout
        << "SAME_RATIO_IDENTITY="
        << stats.same_ratio_identity
        << "\n";

    std::cout
        << "SAME_GAP_IDENTITY="
        << stats.same_gap_identity
        << "\n";

    std::cout
        << "DIFF_PRODUCT_PLUS_OR_MINUS_TWO="
        << stats.diff_product_plus_or_minus_two
        << "\n";

    std::cout
        << "DIFF_GAP_IDENTITY="
        << stats.diff_gap_identity
        << "\n";

    std::cout
        << "\nSAME_K_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        same_k_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nDIFF_TWO_K_GAP_HISTOGRAM\n";

    for (
        const auto& entry :
        diff_k_gap_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_same_non_straddle) {
        std::cout
            << "\nFIRST_SAME_VALUE_NON_STRADDLE\n";

        print_excursion(
            first_same_non_straddle
        );
    } else {
        std::cout
            << "\nFIRST_SAME_VALUE_NON_STRADDLE NONE\n";
    }

    if (have_first_diff_non_straddle) {
        std::cout
            << "\nFIRST_DIFFERENCE_TWO_NON_STRADDLE\n";

        print_excursion(
            first_diff_non_straddle
        );
    } else {
        std::cout
            << "\nFIRST_DIFFERENCE_TWO_NON_STRADDLE NONE\n";
    }

    if (have_first_same_identity_failure) {
        std::cout
            << "\nFIRST_SAME_IDENTITY_FAILURE\n";

        print_excursion(
            first_same_identity_failure
        );
    } else {
        std::cout
            << "\nFIRST_SAME_IDENTITY_FAILURE NONE\n";
    }

    if (have_first_diff_identity_failure) {
        std::cout
            << "\nFIRST_DIFF_IDENTITY_FAILURE\n";

        print_excursion(
            first_diff_identity_failure
        );
    } else {
        std::cout
            << "\nFIRST_DIFF_IDENTITY_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
