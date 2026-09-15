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

static bool divides_u64(
    u64 a,
    u64 b
) {
    return b != 0 && a % b == 0;
}

int main() {
    constexpr int EXPERIMENT = 436;
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
    u64 same_stream_pairs = 0;

    u64 interior_states = 0;

    std::map<int, u64>
        interior_m_histogram;

    std::map<int, u64>
        gcd_m_g_histogram;

    std::map<int, u64>
        gcd_m_q_histogram;

    std::map<int, u64>
        gcd_m_a_histogram;

    std::map<int, u64>
        gcd_m_b_histogram;

    std::map<int, u64>
        m_vs_q_sign_histogram;

    std::map<int, u64>
        m_vs_g_sign_histogram;

    u64 m_divides_q = 0;
    u64 m_not_divides_q = 0;

    u64 m_divides_g = 0;
    u64 m_not_divides_g = 0;

    u64 m_divides_ab = 0;
    u64 m_not_divides_ab = 0;

    u64 m_divides_abq = 0;
    u64 m_not_divides_abq = 0;

    u64 m_less_than_a = 0;
    u64 m_equal_a = 0;
    u64 m_between_a_b = 0;
    u64 m_equal_b = 0;
    u64 m_greater_than_b = 0;

    u64 rule_m_le_6_failures = 0;
    u64 rule_m_divides_12q_failures = 0;
    u64 rule_gcd_m_q_gt1_failures = 0;
    u64 rule_gcd_m_g_gt1_failures = 0;

    u64 m_is_even = 0;
    u64 m_is_odd = 0;

    u64 m_divides_2q = 0;
    u64 m_divides_3q = 0;
    u64 m_divides_4q = 0;
    u64 m_divides_6q = 0;
    u64 m_divides_12q = 0;

    u64 max_interior_m = 0;

    std::map<int, u64>
        max_interior_m_histogram;

    std::map<u64, u64>
        q_over_g_histogram;

    std::map<u64, u64>
        q_mod_m_histogram;

    std::map<u64, u64>
        g_mod_m_histogram;

    std::vector<Witness>
        first_m_gt_6;

    std::vector<Witness>
        first_12q_failure;

    u64 first_m_gt_6_prime = 0;
    u64 first_12q_failure_prime = 0;

    bool have_first_m_gt_6 = false;
    bool have_first_12q_failure = false;

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

                ++total_excursions;

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

                    ++total_mirrored_pairs;

                    if (
                        !same_stream(
                            left,
                            right
                        )
                    ) {
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

                    const u64 g =
                        std::gcd(
                            left.k,
                            right.k
                        );

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

                    q_over_g_histogram[
                        g == 0 ? 0 : q / g
                    ]++;

                    u64 pair_max_m = 0;

                    for (
                        std::size_t p = i + 1;
                        p < j;
                        ++p
                    ) {
                        const u64 m =
                            static_cast<u64>(
                                excursion[p].m
                            );

                        ++interior_states;

                        interior_m_histogram[
                            static_cast<int>(m)
                        ]++;

                        if (m > pair_max_m) {
                            pair_max_m = m;
                        }

                        if (m > max_interior_m) {
                            max_interior_m = m;
                        }

                        if (
                            m > 6
                        ) {
                            ++rule_m_le_6_failures;

                            if (
                                !have_first_m_gt_6
                            ) {
                                first_m_gt_6 =
                                    excursion;

                                first_m_gt_6_prime =
                                    r;

                                have_first_m_gt_6 =
                                    true;
                            }
                        }

                        const u64 gm =
                            std::gcd(
                                m,
                                g
                            );

                        const u64 qm =
                            std::gcd(
                                m,
                                q
                            );

                        const u64 am =
                            std::gcd(
                                m,
                                a
                            );

                        const u64 bm =
                            std::gcd(
                                m,
                                b
                            );

                        gcd_m_g_histogram[
                            static_cast<int>(gm)
                        ]++;

                        gcd_m_q_histogram[
                            static_cast<int>(qm)
                        ]++;

                        gcd_m_a_histogram[
                            static_cast<int>(am)
                        ]++;

                        gcd_m_b_histogram[
                            static_cast<int>(bm)
                        ]++;

                        if (
                            divides_u64(q, m)
                        ) {
                            ++m_divides_q;
                        } else {
                            ++m_not_divides_q;
                        }

                        if (
                            divides_u64(g, m)
                        ) {
                            ++m_divides_g;
                        } else {
                            ++m_not_divides_g;
                        }

                        const u64 ab = a * b;

                        if (
                            divides_u64(ab, m)
                        ) {
                            ++m_divides_ab;
                        } else {
                            ++m_not_divides_ab;
                        }

                        const u64 abq =
                            ab * q;

                        if (
                            divides_u64(abq, m)
                        ) {
                            ++m_divides_abq;
                        } else {
                            ++m_not_divides_abq;
                        }

                        /*
                         * Relative position of m
                         * with respect to a,b.
                         */
                        const u64 low =
                            a < b ? a : b;

                        const u64 high =
                            a < b ? b : a;

                        if (m < low) {
                            ++m_less_than_a;
                        } else if (m == a) {
                            ++m_equal_a;
                        } else if (
                            m > low &&
                            m < high
                        ) {
                            ++m_between_a_b;
                        } else if (m == b) {
                            ++m_equal_b;
                        } else {
                            ++m_greater_than_b;
                        }

                        if ((m % 2) == 0) {
                            ++m_is_even;
                        } else {
                            ++m_is_odd;
                        }

                        if (
                            divides_u64(2 * q, m)
                        ) {
                            ++m_divides_2q;
                        }

                        if (
                            divides_u64(3 * q, m)
                        ) {
                            ++m_divides_3q;
                        }

                        if (
                            divides_u64(4 * q, m)
                        ) {
                            ++m_divides_4q;
                        }

                        if (
                            divides_u64(6 * q, m)
                        ) {
                            ++m_divides_6q;
                        }

                        if (
                            divides_u64(12 * q, m)
                        ) {
                            ++m_divides_12q;
                        } else {
                            ++rule_m_divides_12q_failures;

                            if (
                                !have_first_12q_failure
                            ) {
                                first_12q_failure =
                                    excursion;

                                first_12q_failure_prime =
                                    r;

                                have_first_12q_failure =
                                    true;
                            }
                        }

                        if (
                            qm > 1
                        ) {
                            ++rule_gcd_m_q_gt1_failures;
                        }

                        if (
                            gm > 1
                        ) {
                            ++rule_gcd_m_g_gt1_failures;
                        }

                        q_mod_m_histogram[
                            m == 0 ? 0 : q % m
                        ]++;

                        g_mod_m_histogram[
                            m == 0 ? 0 : g % m
                        ]++;
                    }

                    max_interior_m_histogram[
                        pair_max_m
                    ]++;
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
        << "SAME_STREAM_PAIRS="
        << same_stream_pairs
        << "\n";

    std::cout
        << "\nINTERIOR_MULTIPLIERS\n";

    std::cout
        << "INTERIOR_STATES="
        << interior_states
        << "\n";

    std::cout
        << "MAX_INTERIOR_M="
        << max_interior_m
        << "\n";

    std::cout
        << "M_LE_6_FAILURES="
        << rule_m_le_6_failures
        << "\n";

    std::cout
        << "M_DIVIDES_12Q_FAILURES="
        << rule_m_divides_12q_failures
        << "\n";

    std::cout
        << "\nDIVISIBILITY_RELATIONS\n";

    std::cout
        << "M_DIVIDES_Q="
        << m_divides_q
        << "\n";

    std::cout
        << "M_NOT_DIVIDES_Q="
        << m_not_divides_q
        << "\n";

    std::cout
        << "M_DIVIDES_G="
        << m_divides_g
        << "\n";

    std::cout
        << "M_NOT_DIVIDES_G="
        << m_not_divides_g
        << "\n";

    std::cout
        << "M_DIVIDES_AB="
        << m_divides_ab
        << "\n";

    std::cout
        << "M_NOT_DIVIDES_AB="
        << m_not_divides_ab
        << "\n";

    std::cout
        << "M_DIVIDES_ABQ="
        << m_divides_abq
        << "\n";

    std::cout
        << "M_NOT_DIVIDES_ABQ="
        << m_not_divides_abq
        << "\n";

    std::cout
        << "\nQ_GCD_STRUCTURE\n";

    std::cout
        << "GCD_M_Q_GT1="
        << rule_gcd_m_q_gt1_failures
        << "\n";

    std::cout
        << "GCD_M_G_GT1="
        << rule_gcd_m_g_gt1_failures
        << "\n";

    std::cout
        << "\nM_VS_AB\n";

    std::cout
        << "M_LESS_THAN_A="
        << m_less_than_a
        << "\n";

    std::cout
        << "M_EQUAL_A="
        << m_equal_a
        << "\n";

    std::cout
        << "M_BETWEEN_A_B="
        << m_between_a_b
        << "\n";

    std::cout
        << "M_EQUAL_B="
        << m_equal_b
        << "\n";

    std::cout
        << "M_GREATER_THAN_B="
        << m_greater_than_b
        << "\n";

    std::cout
        << "\nM_DIVIDES_MULTIPLES_OF_Q\n";

    std::cout
        << "M_DIVIDES_2Q="
        << m_divides_2q
        << "\n";

    std::cout
        << "M_DIVIDES_3Q="
        << m_divides_3q
        << "\n";

    std::cout
        << "M_DIVIDES_4Q="
        << m_divides_4q
        << "\n";

    std::cout
        << "M_DIVIDES_6Q="
        << m_divides_6q
        << "\n";

    std::cout
        << "M_DIVIDES_12Q="
        << m_divides_12q
        << "\n";

    std::cout
        << "\nPARITY\n";

    std::cout
        << "M_EVEN="
        << m_is_even
        << "\n";

    std::cout
        << "M_ODD="
        << m_is_odd
        << "\n";

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
        << "\nGCD_M_Q_HISTOGRAM\n";

    for (
        const auto& entry :
        gcd_m_q_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nGCD_M_G_HISTOGRAM\n";

    for (
        const auto& entry :
        gcd_m_g_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nGCD_M_A_HISTOGRAM\n";

    for (
        const auto& entry :
        gcd_m_a_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nGCD_M_B_HISTOGRAM\n";

    for (
        const auto& entry :
        gcd_m_b_histogram
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
        << "\nQ_MOD_M_HISTOGRAM\n";

    for (
        const auto& entry :
        q_mod_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nG_MOD_M_HISTOGRAM\n";

    for (
        const auto& entry :
        g_mod_m_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    if (have_first_m_gt_6) {
        std::cout
            << "\nFIRST_M_GREATER_THAN_6\n";

        print_pair(
            first_m_gt_6_prime,
            first_m_gt_6.front(),
            first_m_gt_6.back()
        );
    } else {
        std::cout
            << "\nFIRST_M_GREATER_THAN_6 NONE\n";
    }

    if (have_first_12q_failure) {
        std::cout
            << "\nFIRST_M_NOT_DIVIDING_12Q\n";

        print_pair(
            first_12q_failure_prime,
            first_12q_failure.front(),
            first_12q_failure.back()
        );
    } else {
        std::cout
            << "\nFIRST_M_NOT_DIVIDING_12Q NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
