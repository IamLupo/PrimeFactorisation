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

struct AAnalysis {
    u64 r = 0;
    int sign = 0;
    u64 A = 0;

    u64 divisor_count = 0;

    u64 winning_divisor_count = 0;
    u64 first_winning_k = 0;
    u64 last_winning_k = 0;

    u64 max_excess_numerator = 0;
    u64 min_winning_excess_numerator = 0;

    u64 winning_below_sqrt = 0;
    u64 winning_above_sqrt = 0;

    u64 winning_records = 0;

    int first_competitor = 0;
    int last_competitor = 0;

    std::vector<u64> winning_k;
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

static u64 normalized_excess_numerator(
    u64 k,
    const Witness& allowed
) {
    /*
     * Compare k/5 against allowed.k/allowed.m.
     *
     * k/5 - allowed.k/allowed.m
     *
     * numerator:
     *
     * k*allowed.m - 5*allowed.k
     */
    const i128 numerator =
        static_cast<i128>(k) *
        static_cast<i128>(allowed.m) -
        static_cast<i128>(5) *
        static_cast<i128>(allowed.k);

    if (numerator <= 0) {
        return 0;
    }

    return static_cast<u64>(numerator);
}

static u64 integer_sqrt(
    u64 n
) {
    u64 x = 1;

    while (
        x <= n / x &&
        (x + 1) <= n / (x + 1)
    ) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static AAnalysis analyze_A(
    u64 r,
    int sign,
    int K_LIMIT
) {
    AAnalysis result;

    result.r = r;
    result.sign = sign;

    const i128 value =
        static_cast<i128>(5) *
        static_cast<i128>(r) +
        static_cast<i128>(sign);

    result.A =
        static_cast<u64>(value);

    /*
     * Track the divisor-record envelope for the
     * allowed multipliers.
     */
    std::vector<Witness> best(
        7
    );

    for (int m = 1; m <= 6; ++m) {
        best[
            static_cast<std::size_t>(m)
        ] =
            make_witness(
                r,
                m,
                1,
                +1
            );
    }

    const u64 root =
        integer_sqrt(result.A);

    /*
     * We need to distinguish:
     *
     * 1. any divisor of A
     * 2. a divisor which actually creates a new
     *    m=5 record
     * 3. a divisor which beats the allowed envelope.
     */
    u64 last_m5_record = 1;

    for (
        u64 k = 1;
        k <= static_cast<u64>(K_LIMIT);
        ++k
    ) {
        /*
         * First update the allowed streams at this K.
         */
        for (int m : {1, 2, 3, 4, 6}) {
            const std::size_t mi =
                static_cast<std::size_t>(m);

            if (
                divides_stream(
                    r,
                    m,
                    k,
                    +1
                )
            ) {
                const Witness candidate =
                    make_witness(
                        r,
                        m,
                        k,
                        +1
                    );

                if (
                    candidate.k >
                    best[mi].k
                ) {
                    best[mi] =
                        candidate;
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
                const Witness candidate =
                    make_witness(
                        r,
                        m,
                        k,
                        -1
                    );

                if (
                    candidate.k >
                    best[mi].k
                ) {
                    best[mi] =
                        candidate;
                }
            }
        }

        if (result.A % k != 0) {
            continue;
        }

        ++result.divisor_count;

        const u64 t =
            result.A / k;

        if (k > last_m5_record) {
            last_m5_record = k;
        }

        const Witness m5 =
            make_witness(
                r,
                5,
                k,
                sign
            );

        Witness strongest =
            best[1];

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

        const u64 excess =
            normalized_excess_numerator(
                k,
                strongest
            );

        if (excess == 0) {
            continue;
        }

        ++result.winning_divisor_count;

        result.winning_k.push_back(k);

        if (result.first_winning_k == 0) {
            result.first_winning_k = k;

            result.first_competitor =
                strongest.m;
        }

        result.last_winning_k = k;
        result.last_competitor =
            strongest.m;

        if (
            result.max_excess_numerator == 0 ||
            excess >
            result.max_excess_numerator
        ) {
            result.max_excess_numerator =
                excess;
        }

        if (
            result.min_winning_excess_numerator ==
            0 ||
            excess <
            result.min_winning_excess_numerator
        ) {
            result.min_winning_excess_numerator =
                excess;
        }

        if (
            k <= root
        ) {
            ++result.winning_below_sqrt;
        } else {
            ++result.winning_above_sqrt;
        }

        if (
            m5.k == k
        ) {
            ++result.winning_records;
        }

        (void)t;
    }

    return result;
}

static void print_vector(
    const std::vector<u64>& values
) {
    for (std::size_t i = 0;
         i < values.size();
         ++i) {

        if (i != 0) {
            std::cout << ",";
        }

        std::cout << values[i];
    }
}

int main() {
    constexpr int EXPERIMENT = 443;
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
    u64 exceptional_A = 0;
    u64 nonexceptional_A = 0;

    u64 total_winning_divisors = 0;

    u64 first_k_sum_exceptional = 0;
    u64 first_k_sum_normal = 0;

    u64 winning_count_sum_exceptional = 0;
    u64 winning_count_sum_normal = 0;

    u64 max_winning_count_exceptional = 0;
    u64 max_winning_count_normal = 0;

    u64 max_excess_exceptional = 0;
    u64 max_excess_normal = 0;

    u64 min_excess_exceptional = 0;
    u64 min_excess_normal = 0;

    u64 winning_below_exceptional = 0;
    u64 winning_above_exceptional = 0;

    u64 exceptional_single_winner = 0;
    u64 exceptional_multiple_winners = 0;

    u64 exceptional_first_competitor_1 = 0;
    u64 exceptional_first_competitor_2 = 0;
    u64 exceptional_first_competitor_3 = 0;
    u64 exceptional_first_competitor_4 = 0;
    u64 exceptional_first_competitor_6 = 0;

    std::vector<AAnalysis>
        exceptional_cases;

    for (u64 r : primes) {
        for (int sign : {-1, +1}) {
            ++total_A;

            const AAnalysis analysis =
                analyze_A(
                    r,
                    sign,
                    K_LIMIT
                );

            if (
                analysis.winning_divisor_count == 0
            ) {
                ++nonexceptional_A;

                if (analysis.first_winning_k != 0) {
                    first_k_sum_normal +=
                        analysis.first_winning_k;
                }

                winning_count_sum_normal +=
                    analysis.winning_divisor_count;

                if (
                    analysis.winning_divisor_count >
                    max_winning_count_normal
                ) {
                    max_winning_count_normal =
                        analysis.winning_divisor_count;
                }

                if (
                    analysis.max_excess_numerator >
                    max_excess_normal
                ) {
                    max_excess_normal =
                        analysis.max_excess_numerator;
                }

                if (
                    analysis.min_winning_excess_numerator >
                    0 &&
                    (
                        min_excess_normal == 0 ||
                        analysis.min_winning_excess_numerator <
                        min_excess_normal
                    )
                ) {
                    min_excess_normal =
                        analysis.min_winning_excess_numerator;
                }

                continue;
            }

            ++exceptional_A;

            exceptional_cases.push_back(
                analysis
            );

            total_winning_divisors +=
                analysis.winning_divisor_count;

            first_k_sum_exceptional +=
                analysis.first_winning_k;

            winning_count_sum_exceptional +=
                analysis.winning_divisor_count;

            if (
                analysis.winning_divisor_count >
                max_winning_count_exceptional
            ) {
                max_winning_count_exceptional =
                    analysis.winning_divisor_count;
            }

            if (
                analysis.max_excess_numerator >
                max_excess_exceptional
            ) {
                max_excess_exceptional =
                    analysis.max_excess_numerator;
            }

            if (
                min_excess_exceptional == 0 ||
                analysis.min_winning_excess_numerator <
                min_excess_exceptional
            ) {
                min_excess_exceptional =
                    analysis.min_winning_excess_numerator;
            }

            winning_below_exceptional +=
                analysis.winning_below_sqrt;

            winning_above_exceptional +=
                analysis.winning_above_sqrt;

            if (
                analysis.winning_divisor_count == 1
            ) {
                ++exceptional_single_winner;
            } else {
                ++exceptional_multiple_winners;
            }

            switch (
                analysis.first_competitor
            ) {
                case 1:
                    ++exceptional_first_competitor_1;
                    break;

                case 2:
                    ++exceptional_first_competitor_2;
                    break;

                case 3:
                    ++exceptional_first_competitor_3;
                    break;

                case 4:
                    ++exceptional_first_competitor_4;
                    break;

                case 6:
                    ++exceptional_first_competitor_6;
                    break;

                default:
                    break;
            }
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "TOTAL_A="
        << total_A
        << "\n";

    std::cout
        << "EXCEPTIONAL_A="
        << exceptional_A
        << "\n";

    std::cout
        << "NONEXCEPTIONAL_A="
        << nonexceptional_A
        << "\n";

    std::cout
        << "TOTAL_WINNING_DIVISORS="
        << total_winning_divisors
        << "\n";

    std::cout
        << "\nEXCEPTIONAL_WINNING_SHAPE\n";

    std::cout
        << "SINGLE_WINNER_A="
        << exceptional_single_winner
        << "\n";

    std::cout
        << "MULTIPLE_WINNER_A="
        << exceptional_multiple_winners
        << "\n";

    std::cout
        << "WINNING_BELOW_SQRT="
        << winning_below_exceptional
        << "\n";

    std::cout
        << "WINNING_ABOVE_SQRT="
        << winning_above_exceptional
        << "\n";

    std::cout
        << "FIRST_COMPETITOR_M1="
        << exceptional_first_competitor_1
        << "\n";

    std::cout
        << "FIRST_COMPETITOR_M2="
        << exceptional_first_competitor_2
        << "\n";

    std::cout
        << "FIRST_COMPETITOR_M3="
        << exceptional_first_competitor_3
        << "\n";

    std::cout
        << "FIRST_COMPETITOR_M4="
        << exceptional_first_competitor_4
        << "\n";

    std::cout
        << "FIRST_COMPETITOR_M6="
        << exceptional_first_competitor_6
        << "\n";

    std::cout
        << "\nAVERAGES\n";

    std::cout
        << "EXCEPTIONAL_FIRST_K_AVG="
        << (
            exceptional_A == 0
                ? 0.0
                : static_cast<double>(
                    first_k_sum_exceptional
                ) /
                  static_cast<double>(
                    exceptional_A
                )
        )
        << "\n";

    std::cout
        << "EXCEPTIONAL_WINNING_COUNT_AVG="
        << (
            exceptional_A == 0
                ? 0.0
                : static_cast<double>(
                    winning_count_sum_exceptional
                ) /
                  static_cast<double>(
                    exceptional_A
                )
        )
        << "\n";

    std::cout
        << "NORMAL_WINNING_COUNT_AVG="
        << (
            nonexceptional_A == 0
                ? 0.0
                : static_cast<double>(
                    winning_count_sum_normal
                ) /
                  static_cast<double>(
                    nonexceptional_A
                )
        )
        << "\n";

    std::cout
        << "\nEXTREMES\n";

    std::cout
        << "MAX_EXCEPTIONAL_WINNING_COUNT="
        << max_winning_count_exceptional
        << "\n";

    std::cout
        << "MAX_NORMAL_WINNING_COUNT="
        << max_winning_count_normal
        << "\n";

    std::cout
        << "MAX_EXCEPTIONAL_EXCESS_NUMERATOR="
        << max_excess_exceptional
        << "\n";

    std::cout
        << "MAX_NORMAL_EXCESS_NUMERATOR="
        << max_excess_normal
        << "\n";

    std::cout
        << "MIN_EXCEPTIONAL_EXCESS_NUMERATOR="
        << min_excess_exceptional
        << "\n";

    std::cout
        << "\nEXCEPTIONAL_CASES\n";

    for (
        const AAnalysis& analysis :
        exceptional_cases
    ) {
        std::cout
            << "R="
            << analysis.r
            << " SIGN="
            << (
                analysis.sign > 0
                    ? "+1"
                    : "-1"
            )
            << " A="
            << analysis.A
            << " WINNING_COUNT="
            << analysis.winning_divisor_count
            << " FIRST_K="
            << analysis.first_winning_k
            << " LAST_K="
            << analysis.last_winning_k
            << " MAX_EXCESS="
            << analysis.max_excess_numerator
            << " MIN_EXCESS="
            << analysis.min_winning_excess_numerator
            << " BELOW="
            << analysis.winning_below_sqrt
            << " ABOVE="
            << analysis.winning_above_sqrt
            << " FIRST_COMPETITOR="
            << analysis.first_competitor
            << " LAST_COMPETITOR="
            << analysis.last_competitor
            << " WINNING_K=";

        print_vector(
            analysis.winning_k
        );

        std::cout
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
