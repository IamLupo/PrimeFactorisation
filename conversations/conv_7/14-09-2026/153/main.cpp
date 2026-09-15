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

struct PairData {
    u64 r = 0;
    int source_sign = +1;
    u64 A = 0;

    u64 k = 0;
    u64 t = 0;

    Witness left;
    Witness right;

    u64 left_margin = 0;
    u64 right_margin = 0;
};

static std::string i128_to_string(i128 value) {
    if (value == 0) {
        return "0";
    }

    bool negative = false;

    if (value < 0) {
        negative = true;
        value = -value;
    }

    char buffer[64];
    int pos = 0;

    while (value > 0) {
        const int digit =
            static_cast<int>(value % 10);

        buffer[pos++] =
            static_cast<char>('0' + digit);

        value /= 10;
    }

    std::string result;

    if (negative) {
        result.push_back('-');
    }

    while (pos > 0) {
        --pos;
        result.push_back(buffer[pos]);
    }

    return result;
}

static std::vector<u64> generate_primes(int limit) {
    std::vector<bool> sieve(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    if (limit >= 0) {
        sieve[0] = false;
    }

    if (limit >= 1) {
        sieve[1] = false;
    }

    for (
        int p = 2;
        static_cast<long long>(p) * p <= limit;
        ++p
    ) {
        if (!sieve[p]) {
            continue;
        }

        for (
            int x = p * p;
            x <= limit;
            x += p
        ) {
            sieve[x] = false;
        }
    }

    std::vector<u64> primes;

    for (int x = 2; x <= limit; ++x) {
        if (sieve[x]) {
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

    return
        value > 0 &&
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

static Witness strongest_allowed(
    const std::vector<Witness>& best
) {
    Witness result = best[1];

    for (int m : {2, 3, 4, 6}) {
        const Witness& candidate =
            best[
                static_cast<std::size_t>(m)
            ];

        if (
            better_normalized(
                candidate,
                result
            )
        ) {
            result = candidate;
        }
    }

    return result;
}

static u64 margin_for(
    u64 k,
    const Witness& allowed
) {
    const i128 value =
        static_cast<i128>(k) *
        static_cast<i128>(allowed.m) -
        static_cast<i128>(5) *
        static_cast<i128>(allowed.k);

    if (value <= 0) {
        return 0;
    }

    return static_cast<u64>(value);
}

static void print_witness(
    const std::string& label,
    const Witness& w
) {
    std::cout
        << label
        << "(m="
        << w.m
        << ",k="
        << w.k
        << ",t="
        << w.t
        << ",sign="
        << (
            w.sign > 0
                ? "+1"
                : "-1"
        )
        << ")";
}

static std::string pair_key(
    int a,
    int b
) {
    return
        std::to_string(a) +
        "->" +
        std::to_string(b);
}

static std::string sign_key(
    int a,
    int b
) {
    return
        std::string(
            a > 0 ? "+" : "-"
        ) +
        "->" +
        std::string(
            b > 0 ? "+" : "-"
        );
}

static i128 product_difference(
    u64 a,
    u64 b,
    u64 c,
    u64 d
) {
    return
        static_cast<i128>(a) *
        static_cast<i128>(b) -
        static_cast<i128>(c) *
        static_cast<i128>(d);
}

int main() {
    constexpr int EXPERIMENT = 445;
    constexpr int PRIME_LIMIT = 10000;
    constexpr int K_LIMIT = 3000;

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

    const std::vector<u64> primes =
        generate_primes(PRIME_LIMIT);

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    u64 total_pairs = 0;

    u64 same_m = 0;
    u64 different_m = 0;

    u64 same_sign = 0;
    u64 different_sign = 0;

    u64 m_sum_5 = 0;
    u64 m_sum_4 = 0;
    u64 m_sum_7 = 0;
    u64 m_sum_8 = 0;
    u64 m_sum_other = 0;

    u64 d_left_divides_right_k = 0;
    u64 d_right_divides_left_k = 0;

    u64 d_left_divides_right_t = 0;
    u64 d_right_divides_left_t = 0;

    u64 cross_product_equal = 0;

    u64 left_margin_less = 0;
    u64 left_margin_equal = 0;
    u64 left_margin_greater = 0;

    u64 margin_sum_equal = 0;

    u64 m_product_4 = 0;
    u64 m_product_6 = 0;
    u64 m_product_7 = 0;
    u64 m_product_12 = 0;
    u64 m_product_other = 0;

    std::map<std::string, u64>
        multiplier_pairs;

    std::map<std::string, u64>
        sign_pairs;

    std::map<int, u64>
        multiplier_sum_histogram;

    u64 max_margin_difference = 0;

    std::vector<PairData> pairs;

    for (u64 r : primes) {
        for (int source_sign : {-1, +1}) {
            const i128 value =
                static_cast<i128>(5) *
                static_cast<i128>(r) +
                static_cast<i128>(source_sign);

            const u64 A =
                static_cast<u64>(value);

            std::vector<
                std::vector<Witness>
            > envelope(
                static_cast<std::size_t>(
                    K_LIMIT + 1
                ),
                std::vector<Witness>(7)
            );

            std::vector<Witness> current(7);

            for (int m = 1; m <= 6; ++m) {
                current[
                    static_cast<std::size_t>(m)
                ] =
                    make_witness(
                        r,
                        m,
                        1,
                        +1
                    );
            }

            /*
             * Build the allowed envelope.
             */
            for (
                u64 K = 1;
                K <= static_cast<u64>(K_LIMIT);
                ++K
            ) {
                for (int m : {1, 2, 3, 4, 6}) {
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
                            current[mi].k
                        ) {
                            current[mi] =
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
                            current[mi].k
                        ) {
                            current[mi] =
                                candidate;
                        }
                    }
                }

                envelope[
                    static_cast<std::size_t>(K)
                ] = current;
            }

            /*
             * Find every winning m=5 divisor.
             */
            std::vector<u64> winners;

            for (
                u64 k = 1;
                k <= static_cast<u64>(K_LIMIT);
                ++k
            ) {
                if (A % k != 0) {
                    continue;
                }

                const Witness allowed =
                    strongest_allowed(
                        envelope[
                            static_cast<std::size_t>(k)
                        ]
                    );

                const u64 margin =
                    margin_for(
                        k,
                        allowed
                    );

                if (margin > 0) {
                    winners.push_back(k);
                }
            }

            /*
             * Pair k with A/k exactly once.
             */
            for (u64 k : winners) {
                const u64 t = A / k;

                if (k >= t) {
                    continue;
                }

                if (t > static_cast<u64>(K_LIMIT)) {
                    continue;
                }

                const Witness left_allowed =
                    strongest_allowed(
                        envelope[
                            static_cast<std::size_t>(k)
                        ]
                    );

                const Witness right_allowed =
                    strongest_allowed(
                        envelope[
                            static_cast<std::size_t>(t)
                        ]
                    );

                const u64 left_margin =
                    margin_for(
                        k,
                        left_allowed
                    );

                const u64 right_margin =
                    margin_for(
                        t,
                        right_allowed
                    );

                ++total_pairs;

                if (
                    left_allowed.m ==
                    right_allowed.m
                ) {
                    ++same_m;
                } else {
                    ++different_m;
                }

                if (
                    left_allowed.sign ==
                    right_allowed.sign
                ) {
                    ++same_sign;
                } else {
                    ++different_sign;
                }

                const int sum_m =
                    left_allowed.m +
                    right_allowed.m;

                ++multiplier_sum_histogram[
                    sum_m
                ];

                if (sum_m == 4) {
                    ++m_sum_4;
                } else if (sum_m == 5) {
                    ++m_sum_5;
                } else if (sum_m == 7) {
                    ++m_sum_7;
                } else if (sum_m == 8) {
                    ++m_sum_8;
                } else {
                    ++m_sum_other;
                }

                const int product_m =
                    left_allowed.m *
                    right_allowed.m;

                if (product_m == 4) {
                    ++m_product_4;
                } else if (product_m == 6) {
                    ++m_product_6;
                } else if (product_m == 7) {
                    ++m_product_7;
                } else if (product_m == 12) {
                    ++m_product_12;
                } else {
                    ++m_product_other;
                }

                ++multiplier_pairs[
                    pair_key(
                        left_allowed.m,
                        right_allowed.m
                    )
                ];

                ++sign_pairs[
                    sign_key(
                        left_allowed.sign,
                        right_allowed.sign
                    )
                ];

                if (
                    k % right_allowed.k == 0
                ) {
                    ++d_right_divides_left_k;
                }

                if (
                    t % left_allowed.k == 0
                ) {
                    ++d_left_divides_right_t;
                }

                if (
                    k % right_allowed.k == 0
                ) {
                    ++d_right_divides_left_k;
                }

                if (
                    t % left_allowed.k == 0
                ) {
                    ++d_left_divides_right_t;
                }

                if (
                    k % right_allowed.t == 0
                ) {
                    ++d_right_divides_left_k;
                }

                if (
                    t % left_allowed.t == 0
                ) {
                    ++d_left_divides_right_k;
                }

                const i128 product_left =
                    static_cast<i128>(k) *
                    static_cast<i128>(left_allowed.m);

                const i128 product_right =
                    static_cast<i128>(left_allowed.k) *
                    static_cast<i128>(5);

                const i128 product_left_2 =
                    static_cast<i128>(t) *
                    static_cast<i128>(right_allowed.m);

                const i128 product_right_2 =
                    static_cast<i128>(right_allowed.k) *
                    static_cast<i128>(5);

                if (
                    product_left ==
                    product_left_2 &&
                    product_right ==
                    product_right_2
                ) {
                    ++cross_product_equal;
                }

                if (
                    left_margin < right_margin
                ) {
                    ++left_margin_less;
                } else if (
                    left_margin == right_margin
                ) {
                    ++left_margin_equal;
                    margin_sum_equal +=
                        left_margin;
                } else {
                    ++left_margin_greater;
                }

                const u64 margin_difference =
                    left_margin >= right_margin
                        ? left_margin -
                          right_margin
                        : right_margin -
                          left_margin;

                if (
                    margin_difference >
                    max_margin_difference
                ) {
                    max_margin_difference =
                        margin_difference;
                }

                PairData pair;

                pair.r = r;
                pair.source_sign = source_sign;
                pair.A = A;
                pair.k = k;
                pair.t = t;
                pair.left = left_allowed;
                pair.right = right_allowed;
                pair.left_margin =
                    left_margin;
                pair.right_margin =
                    right_margin;

                pairs.push_back(pair);
            }
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_PAIRS="
        << total_pairs
        << "\n";

    std::cout
        << "SAME_M="
        << same_m
        << "\n";

    std::cout
        << "DIFFERENT_M="
        << different_m
        << "\n";

    std::cout
        << "SAME_SIGN="
        << same_sign
        << "\n";

    std::cout
        << "DIFFERENT_SIGN="
        << different_sign
        << "\n";

    std::cout
        << "\nMULTIPLIER_SUMS\n";

    std::cout
        << "SUM_4="
        << m_sum_4
        << "\n";

    std::cout
        << "SUM_5="
        << m_sum_5
        << "\n";

    std::cout
        << "SUM_7="
        << m_sum_7
        << "\n";

    std::cout
        << "SUM_8="
        << m_sum_8
        << "\n";

    std::cout
        << "SUM_OTHER="
        << m_sum_other
        << "\n";

    std::cout
        << "\nMULTIPLIER_PRODUCTS\n";

    std::cout
        << "PRODUCT_4="
        << m_product_4
        << "\n";

    std::cout
        << "PRODUCT_6="
        << m_product_6
        << "\n";

    std::cout
        << "PRODUCT_7="
        << m_product_7
        << "\n";

    std::cout
        << "PRODUCT_12="
        << m_product_12
        << "\n";

    std::cout
        << "PRODUCT_OTHER="
        << m_product_other
        << "\n";

    std::cout
        << "\nDIVISIBILITY_RELATIONS\n";

    std::cout
        << "RIGHT_K_DIVIDES_LEFT_K="
        << d_right_divides_left_k
        << "\n";

    std::cout
        << "LEFT_K_DIVIDES_RIGHT_T="
        << d_left_divides_right_t
        << "\n";

    std::cout
        << "RIGHT_T_DIVIDES_LEFT_K="
        << d_right_divides_left_k
        << "\n";

    std::cout
        << "LEFT_T_DIVIDES_RIGHT_K="
        << d_left_divides_right_k
        << "\n";

    std::cout
        << "CROSS_PRODUCT_EQUAL="
        << cross_product_equal
        << "\n";

    std::cout
        << "\nMARGIN_ORDER\n";

    std::cout
        << "LEFT_MARGIN_LESS="
        << left_margin_less
        << "\n";

    std::cout
        << "MARGIN_EQUAL="
        << left_margin_equal
        << "\n";

    std::cout
        << "LEFT_MARGIN_GREATER="
        << left_margin_greater
        << "\n";

    std::cout
        << "MAX_MARGIN_DIFFERENCE="
        << max_margin_difference
        << "\n";

    std::cout
        << "\nMULTIPLIER_PAIR_HISTOGRAM\n";

    for (
        const auto& entry :
        multiplier_pairs
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nSIGN_PAIR_HISTOGRAM\n";

    for (
        const auto& entry :
        sign_pairs
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nMULTIPLIER_SUM_HISTOGRAM\n";

    for (
        const auto& entry :
        multiplier_sum_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nALL_RECIPROCAL_PAIRS\n";

    for (const PairData& pair : pairs) {
        std::cout
            << "R="
            << pair.r
            << " SIGN="
            << (
                pair.source_sign > 0
                    ? "+1"
                    : "-1"
            )
            << " A="
            << pair.A
            << " K="
            << pair.k
            << " T="
            << pair.t
            << " ";

        print_witness(
            "LEFT_ALLOWED",
            pair.left
        );

        std::cout
            << " ";

        print_witness(
            "RIGHT_ALLOWED",
            pair.right
        );

        std::cout
            << " LEFT_MARGIN="
            << pair.left_margin
            << " RIGHT_MARGIN="
            << pair.right_margin
            << " M_SUM="
            << (
                pair.left.m +
                pair.right.m
            )
            << " M_PRODUCT="
            << (
                pair.left.m *
                pair.right.m
            )
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
