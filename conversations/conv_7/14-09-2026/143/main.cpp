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

static u64 gcd3(
    u64 a,
    u64 b,
    u64 c
) {
    return std::gcd(
        std::gcd(a, b),
        c
    );
}

int main() {
    constexpr int EXPERIMENT = 435;
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

    u64 canonical_pairs = 0;

    u64 q_equals_g = 0;
    u64 q_less_g = 0;
    u64 q_greater_g = 0;

    u64 ratio_integer = 0;
    u64 ratio_noninteger = 0;

    u64 ratio_one = 0;
    u64 ratio_gt_one = 0;

    u64 interior_state_pairs = 0;

    u64 interior_m_divides_a = 0;
    u64 interior_m_divides_b = 0;
    u64 interior_m_divides_q = 0;
    u64 interior_m_divides_g = 0;

    u64 interior_m_equals_a = 0;
    u64 interior_m_equals_b = 0;
    u64 interior_m_equals_q = 0;
    u64 interior_m_equals_g = 0;

    u64 interior_m_coprime_a = 0;
    u64 interior_m_coprime_b = 0;

    u64 transition_m_in_a_interval = 0;
    u64 transition_m_in_b_interval = 0;

    u64 all_interior_m_divide_abq = 0;

    u64 endpoint_a_b_order_ok = 0;
    u64 endpoint_q_range_ok = 0;

    std::map<u64, u64>
        q_over_g_histogram;

    std::map<u64, u64>
        max_interior_m_histogram;

    std::map<int, u64>
        interior_m_histogram;

    std::map<u64, u64>
        endpoint_product_ratio_histogram;

    std::map<u64, u64>
        q_histogram;

    std::map<u64, u64>
        g_histogram;

    std::vector<Witness>
        first_interior_counterexample;

    u64 first_interior_counterexample_prime = 0;
    bool have_first_interior_counterexample = false;

    std::vector<Witness>
        first_abq_failure;

    u64 first_abq_failure_prime = 0;
    bool have_first_abq_failure = false;

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
                end >=
                transition_states.size()
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
                    ++canonical_pairs;

                    const u64 A =
                        static_cast<u64>(
                            arithmetic_value(
                                r,
                                left
                            )
                        );

                    const u64 g =
                        std::gcd(
                            left.k,
                            right.k
                        );

                    if (g == 0) {
                        continue;
                    }

                    const u64 a =
                        left.k / g;

                    const u64 b =
                        right.k / g;

                    const u64 gab =
                        g * a * b;

                    if (
                        gab == 0 ||
                        A % gab != 0
                    ) {
                        continue;
                    }

                    const u64 q =
                        A / gab;

                    g_histogram[g]++;
                    q_histogram[q]++;

                    const u64 ratio_integer_part =
                        q / g;

                    if (q % g == 0) {
                        ++ratio_integer;
                        ++q_over_g_histogram[
                            ratio_integer_part
                        ];

                        if (q == g) {
                            ++q_equals_g;
                            ++ratio_one;
                        } else {
                            ++q_greater_g;
                            ++ratio_gt_one;
                        }
                    } else {
                        ++ratio_noninteger;

                        if (q < g) {
                            ++q_less_g;
                        } else {
                            ++q_greater_g;
                        }
                    }

                    /*
                     * Inspect every interior projected m-state.
                     */
                    u64 max_interior_m = 0;

                    for (
                        std::size_t p = i + 1;
                        p < j;
                        ++p
                    ) {
                        const int m =
                            excursion[p].m;

                        ++interior_state_pairs;
                        ++interior_m_histogram[m];

                        if (
                            static_cast<u64>(m) >
                            max_interior_m
                        ) {
                            max_interior_m =
                                static_cast<u64>(m);
                        }

                        if (
                            a % static_cast<u64>(m) ==
                            0
                        ) {
                            ++interior_m_divides_a;
                        }

                        if (
                            b % static_cast<u64>(m) ==
                            0
                        ) {
                            ++interior_m_divides_b;
                        }

                        if (
                            q % static_cast<u64>(m) ==
                            0
                        ) {
                            ++interior_m_divides_q;
                        }

                        if (
                            g % static_cast<u64>(m) ==
                            0
                        ) {
                            ++interior_m_divides_g;
                        }

                        if (
                            static_cast<u64>(m) ==
                            a
                        ) {
                            ++interior_m_equals_a;
                        }

                        if (
                            static_cast<u64>(m) ==
                            b
                        ) {
                            ++interior_m_equals_b;
                        }

                        if (
                            static_cast<u64>(m) ==
                            q
                        ) {
                            ++interior_m_equals_q;
                        }

                        if (
                            static_cast<u64>(m) ==
                            g
                        ) {
                            ++interior_m_equals_g;
                        }

                        if (
                            std::gcd(
                                static_cast<u64>(m),
                                a
                            ) == 1
                        ) {
                            ++interior_m_coprime_a;
                        }

                        if (
                            std::gcd(
                                static_cast<u64>(m),
                                b
                            ) == 1
                        ) {
                            ++interior_m_coprime_b;
                        }

                        /*
                         * Does this interior multiplier sit
                         * numerically between a and b?
                         */
                        const u64 low_ab =
                            a < b ? a : b;

                        const u64 high_ab =
                            a < b ? b : a;

                        if (
                            static_cast<u64>(m) >
                            low_ab &&
                            static_cast<u64>(m) <
                            high_ab
                        ) {
                            ++transition_m_in_a_interval;
                        }

                        /*
                         * Does m divide abq?
                         */
                        const u64 abq =
                            a * b * q;

                        if (
                            abq %
                            static_cast<u64>(m) ==
                            0
                        ) {
                            ++all_interior_m_divide_abq;
                        } else if (
                            !have_first_abq_failure
                        ) {
                            first_abq_failure =
                                excursion;

                            first_abq_failure_prime =
                                r;

                            have_first_abq_failure =
                                true;
                        }

                        /*
                         * A strict numerical range test against q.
                         */
                        if (
                            static_cast<u64>(m) <= q
                        ) {
                            ++endpoint_q_range_ok;
                        }
                    }

                    max_interior_m_histogram[
                        max_interior_m
                    ]++;

                    /*
                     * The reduced endpoint pair has b>a
                     * because k_L<k_R.
                     */
                    if (a < b) {
                        ++endpoint_a_b_order_ok;
                    }

                    /*
                     * Measure the normalized endpoint product:
                     *
                     * kL*kR / A = g/q.
                     */
                    if (q != 0) {
                        endpoint_product_ratio_histogram[
                            g / std::gcd(g, q)
                        ]++;
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
        << "\nCANONICAL_PARAMETERIZATION\n";

    std::cout
        << "CANONICAL_PAIRS="
        << canonical_pairs
        << "\n";

    std::cout
        << "Q_EQUALS_G="
        << q_equals_g
        << "\n";

    std::cout
        << "Q_LESS_G="
        << q_less_g
        << "\n";

    std::cout
        << "Q_GREATER_G="
        << q_greater_g
        << "\n";

    std::cout
        << "RATIO_INTEGER="
        << ratio_integer
        << "\n";

    std::cout
        << "RATIO_NONINTEGER="
        << ratio_noninteger
        << "\n";

    std::cout
        << "RATIO_ONE="
        << ratio_one
        << "\n";

    std::cout
        << "RATIO_GREATER_THAN_ONE="
        << ratio_gt_one
        << "\n";

    std::cout
        << "\nINTERIOR_M_ANALYSIS\n";

    std::cout
        << "INTERIOR_STATE_COUNT="
        << interior_state_pairs
        << "\n";

    std::cout
        << "INTERIOR_M_DIVIDES_A="
        << interior_m_divides_a
        << "\n";

    std::cout
        << "INTERIOR_M_DIVIDES_B="
        << interior_m_divides_b
        << "\n";

    std::cout
        << "INTERIOR_M_DIVIDES_Q="
        << interior_m_divides_q
        << "\n";

    std::cout
        << "INTERIOR_M_DIVIDES_G="
        << interior_m_divides_g
        << "\n";

    std::cout
        << "INTERIOR_M_EQUALS_A="
        << interior_m_equals_a
        << "\n";

    std::cout
        << "INTERIOR_M_EQUALS_B="
        << interior_m_equals_b
        << "\n";

    std::cout
        << "INTERIOR_M_EQUALS_Q="
        << interior_m_equals_q
        << "\n";

    std::cout
        << "INTERIOR_M_EQUALS_G="
        << interior_m_equals_g
        << "\n";

    std::cout
        << "INTERIOR_M_COPRIME_A="
        << interior_m_coprime_a
        << "\n";

    std::cout
        << "INTERIOR_M_COPRIME_B="
        << interior_m_coprime_b
        << "\n";

    std::cout
        << "INTERIOR_M_IN_A_B_INTERVAL="
        << transition_m_in_a_interval
        << "\n";

    std::cout
        << "INTERIOR_M_DIVIDES_ABQ="
        << all_interior_m_divide_abq
        << "\n";

    std::cout
        << "\nMAX_INTERIOR_M_HISTOGRAM\n";

    for (
        const auto& entry :
        max_interior_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nINTERIOR_M_HISTOGRAM\n";

    for (
        const auto& entry :
        interior_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nQ_OVER_G_HISTOGRAM\n";

    for (
        const auto& entry :
        q_over_g_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nG_HISTOGRAM\n";

    for (
        const auto& entry :
        g_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nQ_HISTOGRAM\n";

    for (
        const auto& entry :
        q_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nENDPOINT_ORDER\n";

    std::cout
        << "A_LESS_THAN_B="
        << endpoint_a_b_order_ok
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
