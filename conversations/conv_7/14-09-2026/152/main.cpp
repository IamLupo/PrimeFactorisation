#include <cstdint>
#include <iostream>
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

struct WinningPoint {
    u64 k = 0;
    u64 t = 0;
    u64 margin = 0;
    int competitor = 0;
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

static std::vector<u64> generate_primes(
    int limit
) {
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

static int strongest_allowed_multiplier(
    const std::vector<Witness>& best
) {
    Witness strongest = best[1];
    int strongest_m = 1;

    for (int m : {2, 3, 4, 6}) {
        const Witness& candidate =
            best[
                static_cast<std::size_t>(m)
            ];

        if (
            better_normalized(
                candidate,
                strongest
            )
        ) {
            strongest = candidate;
            strongest_m = m;
        }
    }

    return strongest_m;
}

static Witness strongest_allowed_witness(
    const std::vector<Witness>& best
) {
    Witness strongest = best[1];

    for (int m : {2, 3, 4, 6}) {
        const Witness& candidate =
            best[
                static_cast<std::size_t>(m)
            ];

        if (
            better_normalized(
                candidate,
                strongest
            )
        ) {
            strongest = candidate;
        }
    }

    return strongest;
}

static i128 normalized_cross(
    u64 k,
    int exceptional_m,
    const Witness& allowed
) {
    return
        static_cast<i128>(k) *
        static_cast<i128>(allowed.m) -
        static_cast<i128>(exceptional_m) *
        static_cast<i128>(allowed.k);
}

static bool is_winner_at(
    u64 r,
    int sign,
    u64 A,
    u64 k,
    const std::vector<Witness>& envelope,
    WinningPoint& output
) {
    if (k == 0 || k > 3000) {
        return false;
    }

    if (A % k != 0) {
        return false;
    }

    const Witness allowed =
        strongest_allowed_witness(
            envelope
        );

    const i128 cross =
        normalized_cross(
            k,
            5,
            allowed
        );

    if (cross <= 0) {
        return false;
    }

    output.k = k;
    output.t = A / k;
    output.margin =
        static_cast<u64>(cross);
    output.competitor =
        allowed.m;

    (void)r;
    (void)sign;

    return true;
}

static void print_winning_point(
    const WinningPoint& point,
    const std::string& label
) {
    std::cout
        << label
        << "(k="
        << point.k
        << ",t="
        << point.t
        << ",margin="
        << point.margin
        << ",competitor=m"
        << point.competitor
        << ")";
}

int main() {
    constexpr int EXPERIMENT = 444;
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
        generate_primes(
            PRIME_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    u64 total_A = 0;

    u64 total_winning_points = 0;
    u64 reciprocal_pairs = 0;
    u64 reciprocal_failures = 0;

    u64 margin_equal = 0;
    u64 margin_not_equal = 0;

    u64 competitor_equal = 0;
    u64 competitor_not_equal = 0;

    u64 first_is_below = 0;
    u64 first_is_above = 0;

    u64 pair_sum_margin_equal = 0;
    u64 pair_sum_margin_not_equal = 0;

    u64 strict_margin_reverse = 0;

    u64 prime_with_exception = 0;

    u64 max_margin_difference = 0;

    std::vector<std::string>
        failure_examples;

    for (u64 r : primes) {
        bool prime_had_exception = false;

        for (int sign : {-1, +1}) {
            ++total_A;

            const i128 value =
                static_cast<i128>(5) *
                static_cast<i128>(r) +
                static_cast<i128>(sign);

            const u64 A =
                static_cast<u64>(value);

            /*
             * The allowed envelope is a function of K.
             *
             * To test both k and t we need all allowed records
             * up to K_LIMIT.
             */
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
             * Collect all m=5 winning divisors.
             */
            std::vector<WinningPoint> winners;

            for (
                u64 k = 1;
                k <= static_cast<u64>(K_LIMIT);
                ++k
            ) {
                if (A % k != 0) {
                    continue;
                }

                const Witness& allowed =
                    strongest_allowed_witness(
                        envelope[
                            static_cast<std::size_t>(k)
                        ]
                    );

                const i128 cross =
                    normalized_cross(
                        k,
                        5,
                        allowed
                    );

                if (cross <= 0) {
                    continue;
                }

                WinningPoint point;

                point.k = k;
                point.t = A / k;
                point.margin =
                    static_cast<u64>(cross);
                point.competitor =
                    allowed.m;

                winners.push_back(point);
            }

            if (winners.empty()) {
                continue;
            }

            prime_had_exception = true;

            /*
             * Test every winning divisor against its complement.
             */
            for (
                const WinningPoint& left :
                winners
            ) {
                ++total_winning_points;

                const u64 complement =
                    left.t;

                WinningPoint right;

                const bool reverse_exists =
                    is_winner_at(
                        r,
                        sign,
                        A,
                        complement,
                        envelope[
                            static_cast<std::size_t>(
                                complement <= K_LIMIT
                                    ? complement
                                    : 0
                            )
                        ],
                        right
                    );

                if (
                    complement == 0 ||
                    complement > K_LIMIT
                ) {
                    ++reciprocal_failures;

                    if (
                        failure_examples.size() < 10
                    ) {
                        failure_examples.push_back(
                            "R=" +
                            std::to_string(r) +
                            " SIGN=" +
                            (
                                sign > 0
                                    ? "+1"
                                    : "-1"
                            ) +
                            " K=" +
                            std::to_string(left.k) +
                            " T=" +
                            std::to_string(complement) +
                            " FAILURE=COMPLEMENT_OUTSIDE_LIMIT"
                        );
                    }

                    continue;
                }

                if (!reverse_exists) {
                    ++reciprocal_failures;

                    if (
                        failure_examples.size() < 10
                    ) {
                        failure_examples.push_back(
                            "R=" +
                            std::to_string(r) +
                            " SIGN=" +
                            (
                                sign > 0
                                    ? "+1"
                                    : "-1"
                            ) +
                            " K=" +
                            std::to_string(left.k) +
                            " T=" +
                            std::to_string(complement) +
                            " FAILURE=REVERSE_NOT_WINNING"
                        );
                    }

                    continue;
                }

                ++reciprocal_pairs;

                if (left.k < left.t) {
                    ++first_is_below;
                } else if (left.k > left.t) {
                    ++first_is_above;
                }

                const u64 margin_difference =
                    left.margin >= right.margin
                        ? left.margin -
                          right.margin
                        : right.margin -
                          left.margin;

                if (
                    margin_difference >
                    max_margin_difference
                ) {
                    max_margin_difference =
                        margin_difference;
                }

                if (
                    left.margin ==
                    right.margin
                ) {
                    ++margin_equal;
                    pair_sum_margin_equal +=
                        left.margin;
                } else {
                    ++margin_not_equal;
                    pair_sum_margin_not_equal +=
                        left.margin +
                        right.margin;
                }

                if (
                    left.competitor ==
                    right.competitor
                ) {
                    ++competitor_equal;
                } else {
                    ++competitor_not_equal;
                }

                /*
                 * The pair should be counted only once when k<t.
                 */
                if (left.k < left.t) {
                    /*
                     * Check whether the reverse point actually
                     * contains exactly the complementary divisor.
                     */
                    if (right.t == left.k) {
                        ++strict_margin_reverse;
                    }
                }
            }

            /*
             * Print one compact block per exceptional A.
             */
            std::cout
                << "\nEXCEPTIONAL_A\n";

            std::cout
                << "R="
                << r
                << " SIGN="
                << (
                    sign > 0
                        ? "+1"
                        : "-1"
                )
                << " A="
                << A
                << " WINNERS="
                << winners.size()
                << "\n";

            for (
                const WinningPoint& point :
                winners
            ) {
                const u64 complement =
                    point.t;

                const bool reverse =
                    complement <=
                        static_cast<u64>(
                            K_LIMIT
                        );

                if (!reverse) {
                    std::cout
                        << "POINT ";
                    print_winning_point(
                        point,
                        "LEFT"
                    );
                    std::cout
                        << " REVERSE=OUT_OF_RANGE\n";
                    continue;
                }

                WinningPoint reverse_point;

                const bool reverse_wins =
                    is_winner_at(
                        r,
                        sign,
                        A,
                        complement,
                        envelope[
                            static_cast<std::size_t>(
                                complement
                            )
                        ],
                        reverse_point
                    );

                std::cout
                    << "POINT ";

                print_winning_point(
                    point,
                    "LEFT"
                );

                if (reverse_wins) {
                    std::cout
                        << " RIGHT=";

                    print_winning_point(
                        reverse_point,
                        ""
                    );

                    std::cout
                        << " MARGIN_DIFF="
                        << (
                            point.margin >=
                            reverse_point.margin
                                ? point.margin -
                                  reverse_point.margin
                                : reverse_point.margin -
                                  point.margin
                        )
                        << "\n";
                } else {
                    std::cout
                        << " RIGHT=NOT_WINNING\n";
                }
            }
        }

        if (prime_had_exception) {
            ++prime_with_exception;
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_A="
        << total_A
        << "\n";

    std::cout
        << "PRIMES_WITH_EXCEPTION="
        << prime_with_exception
        << "\n";

    std::cout
        << "TOTAL_WINNING_POINTS="
        << total_winning_points
        << "\n";

    std::cout
        << "RECIPROCAL_PAIRS="
        << reciprocal_pairs
        << "\n";

    std::cout
        << "RECIPROCAL_FAILURES="
        << reciprocal_failures
        << "\n";

    std::cout
        << "MARGIN_EQUAL="
        << margin_equal
        << "\n";

    std::cout
        << "MARGIN_NOT_EQUAL="
        << margin_not_equal
        << "\n";

    std::cout
        << "COMPETITOR_EQUAL="
        << competitor_equal
        << "\n";

    std::cout
        << "COMPETITOR_NOT_EQUAL="
        << competitor_not_equal
        << "\n";

    std::cout
        << "FIRST_POINT_BELOW_SQRT="
        << first_is_below
        << "\n";

    std::cout
        << "FIRST_POINT_ABOVE_SQRT="
        << first_is_above
        << "\n";

    std::cout
        << "MAX_MARGIN_DIFFERENCE="
        << max_margin_difference
        << "\n";

    std::cout
        << "REVERSE_COMPLEMENT_MATCH="
        << strict_margin_reverse
        << "\n";

    std::cout
        << "\nFAILURE_EXAMPLES\n";

    for (
        const std::string& example :
        failure_examples
    ) {
        std::cout
            << example
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
