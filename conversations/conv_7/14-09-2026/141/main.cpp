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
    constexpr int EXPERIMENT = 433;
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

    u64 canonical_success = 0;
    u64 canonical_failure = 0;

    u64 q_equals_g = 0;
    u64 q_less_g = 0;
    u64 q_greater_g = 0;

    u64 g_divides_q = 0;
    u64 q_divides_g = 0;

    u64 ratio_integer = 0;
    u64 ratio_noninteger = 0;

    u64 ratio_one = 0;
    u64 ratio_less_than_one = 0;
    u64 ratio_greater_than_one = 0;

    u64 a_equals_one = 0;
    u64 b_equals_one = 0;

    u64 g_equals_one = 0;
    u64 g_gt_one = 0;

    u64 q_equals_one = 0;
    u64 q_gt_one = 0;

    u64 normalized_product_identity = 0;

    u64 equation_a_b_q = 0;

    std::map<u64, u64>
        g_histogram;

    std::map<u64, u64>
        a_histogram;

    std::map<u64, u64>
        b_histogram;

    std::map<u64, u64>
        q_histogram;

    std::map<u64, u64>
        q_over_g_histogram;

    std::map<u64, u64>
        g_over_q_histogram;

    std::map<u64, u64>
        gcd_ab_histogram;

    std::vector<Witness> first_same_stream_pair;
    u64 first_same_stream_prime = 0;
    bool have_first_same_stream = false;

    std::vector<Witness> first_failure_pair;
    u64 first_failure_prime = 0;
    bool have_first_failure = false;

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

                    const bool same_stream =
                        left.m == right.m &&
                        left.sign == right.sign;

                    if (!same_stream) {
                        ++sign_switch_pairs;
                        continue;
                    }

                    ++same_stream_pairs;

                    if (!have_first_same_stream) {
                        first_same_stream_pair =
                            excursion;
                        first_same_stream_prime =
                            r;
                        have_first_same_stream = true;
                    }

                    const u64 A =
                        static_cast<u64>(
                            arithmetic_value(
                                r,
                                left
                            )
                        );

                    /*
                     * Canonical decomposition:
                     *
                     * g = gcd(kL,kR)
                     * kL = g*a
                     * kR = g*b
                     *
                     * gcd(a,b)=1
                     *
                     * A = g*a*b*q
                     *
                     * tL = b*q
                     * tR = a*q
                     */
                    const u64 g =
                        std::gcd(
                            left.k,
                            right.k
                        );

                    if (
                        g == 0 ||
                        left.k % g != 0 ||
                        right.k % g != 0
                    ) {
                        ++canonical_failure;

                        if (!have_first_failure) {
                            first_failure_pair =
                                excursion;
                            first_failure_prime =
                                r;
                            have_first_failure = true;
                        }

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
                        ++canonical_failure;

                        if (!have_first_failure) {
                            first_failure_pair =
                                excursion;
                            first_failure_prime =
                                r;
                            have_first_failure = true;
                        }

                        continue;
                    }

                    const u64 q =
                        A / gab;

                    const bool equations_ok =
                        std::gcd(a, b) == 1 &&
                        left.t == b * q &&
                        right.t == a * q &&
                        A == g * a * b * q;

                    if (equations_ok) {
                        ++canonical_success;
                    } else {
                        ++canonical_failure;

                        if (!have_first_failure) {
                            first_failure_pair =
                                excursion;
                            first_failure_prime =
                                r;
                            have_first_failure = true;
                        }

                        continue;
                    }

                    ++equation_a_b_q;

                    g_histogram[g]++;
                    a_histogram[a]++;
                    b_histogram[b]++;
                    q_histogram[q]++;

                    gcd_ab_histogram[
                        std::gcd(a, b)
                    ]++;

                    if (g == 1) {
                        ++g_equals_one;
                    } else {
                        ++g_gt_one;
                    }

                    if (q == 1) {
                        ++q_equals_one;
                    } else {
                        ++q_gt_one;
                    }

                    if (a == 1) {
                        ++a_equals_one;
                    }

                    if (b == 1) {
                        ++b_equals_one;
                    }

                    if (q == g) {
                        ++q_equals_g;
                    } else if (q < g) {
                        ++q_less_g;
                    } else {
                        ++q_greater_g;
                    }

                    if (q % g == 0) {
                        ++g_divides_q;

                        const u64 ratio =
                            q / g;

                        q_over_g_histogram[
                            ratio
                        ]++;

                        ++ratio_integer;

                        if (ratio == 1) {
                            ++ratio_one;
                        } else if (ratio < 1) {
                            ++ratio_less_than_one;
                        } else {
                            ++ratio_greater_than_one;
                        }
                    } else {
                        ++ratio_noninteger;
                    }

                    if (g % q == 0) {
                        ++q_divides_g;

                        const u64 ratio =
                            g / q;

                        g_over_q_histogram[
                            ratio
                        ]++;
                    }

                    const i128 product =
                        static_cast<i128>(
                            left.k
                        ) *
                        static_cast<i128>(
                            right.k
                        );

                    if (
                        product ==
                        static_cast<i128>(A)
                    ) {
                        ++normalized_product_identity;
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
        << "CANONICAL_SUCCESS="
        << canonical_success
        << "\n";

    std::cout
        << "CANONICAL_FAILURE="
        << canonical_failure
        << "\n";

    std::cout
        << "EQUATION_A_B_Q="
        << equation_a_b_q
        << "\n";

    std::cout
        << "\nG_RELATION\n";

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
        << "G_DIVIDES_Q="
        << g_divides_q
        << "\n";

    std::cout
        << "Q_DIVIDES_G="
        << q_divides_g
        << "\n";

    std::cout
        << "\nRATIO_Q_OVER_G\n";

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
        << "RATIO_LESS_THAN_ONE="
        << ratio_less_than_one
        << "\n";

    std::cout
        << "RATIO_GREATER_THAN_ONE="
        << ratio_greater_than_one
        << "\n";

    std::cout
        << "\nBASIC_PARAMETERS\n";

    std::cout
        << "G_EQUALS_ONE="
        << g_equals_one
        << "\n";

    std::cout
        << "G_GT_ONE="
        << g_gt_one
        << "\n";

    std::cout
        << "Q_EQUALS_ONE="
        << q_equals_one
        << "\n";

    std::cout
        << "Q_GT_ONE="
        << q_gt_one
        << "\n";

    std::cout
        << "A_EQUALS_ONE="
        << a_equals_one
        << "\n";

    std::cout
        << "B_EQUALS_ONE="
        << b_equals_one
        << "\n";

    std::cout
        << "\nPRODUCT_IDENTITY\n";

    std::cout
        << "K_LEFT_K_RIGHT_EQUALS_A="
        << normalized_product_identity
        << "\n";

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
        << "\nA_HISTOGRAM\n";

    for (
        const auto& entry :
        a_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nB_HISTOGRAM\n";

    for (
        const auto& entry :
        b_histogram
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
        << "\nG_OVER_Q_HISTOGRAM\n";

    for (
        const auto& entry :
        g_over_q_histogram
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

        print_pair(
            first_same_stream_prime,
            first_same_stream_pair.front(),
            first_same_stream_pair.back()
        );
    } else {
        std::cout
            << "\nFIRST_SAME_STREAM_PAIR NONE\n";
    }

    if (have_first_failure) {
        std::cout
            << "\nFIRST_CANONICAL_FAILURE\n";

        print_pair(
            first_failure_prime,
            first_failure_pair.front(),
            first_failure_pair.back()
        );
    } else {
        std::cout
            << "\nFIRST_CANONICAL_FAILURE NONE\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
