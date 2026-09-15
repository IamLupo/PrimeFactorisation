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

struct ExceptionalPoint {
    u64 prime = 0;
    u64 K = 0;
    Witness winner;
    Witness allowed_winner;
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

static bool allowed_m(int m) {
    return
        m == 1 ||
        m == 2 ||
        m == 3 ||
        m == 4 ||
        m == 6;
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

static Witness allowed_winner(
    const std::vector<Witness>& best_by_m
) {
    Witness best = best_by_m[1];

    for (int m : {2, 3, 4, 6}) {
        const Witness& candidate =
            best_by_m[
                static_cast<std::size_t>(m)
            ];

        if (
            better_normalized(
                candidate,
                best
            )
        ) {
            best = candidate;
        }
    }

    return best;
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

static u64 gcd_u64(
    u64 a,
    u64 b
) {
    while (b != 0) {
        const u64 t = a % b;
        a = b;
        b = t;
    }

    return a;
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

int main() {
    constexpr int EXPERIMENT = 438;
    constexpr int PRIME_LIMIT = 10000;
    constexpr int K_LIMIT = 3000;
    constexpr int M_MAX = 64;

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

    u64 exceptional_points = 0;
    u64 m5_points = 0;
    u64 m7_points = 0;

    u64 m5_prime_count = 0;
    u64 m7_prime_count = 0;

    u64 m5_run_count = 0;
    u64 m7_run_count = 0;

    u64 m5_multi_point_runs = 0;
    u64 m7_multi_point_runs = 0;

    u64 max_m5_run = 0;
    u64 max_m7_run = 0;

    u64 m5_first_K_sum = 0;
    u64 m7_first_K_sum = 0;

    u64 m5_k_equals_t = 0;
    u64 m7_k_equals_t = 0;

    u64 m5_k_less_t = 0;
    u64 m5_k_greater_t = 0;

    u64 m7_k_less_t = 0;
    u64 m7_k_greater_t = 0;

    u64 m5_sign_plus = 0;
    u64 m5_sign_minus = 0;

    u64 m7_sign_plus = 0;
    u64 m7_sign_minus = 0;

    u64 m5_allowed_margin_positive = 0;
    u64 m7_allowed_margin_positive = 0;

    u64 m5_allowed_margin_sum = 0;
    u64 m7_allowed_margin_sum = 0;

    u64 m5_gcd_kr = 0;
    u64 m7_gcd_kr = 0;

    std::map<u64, u64>
        m5_K_histogram;

    std::map<u64, u64>
        m7_K_histogram;

    std::map<u64, u64>
        m5_k_histogram;

    std::map<u64, u64>
        m7_k_histogram;

    std::map<u64, u64>
        m5_t_histogram;

    std::map<u64, u64>
        m7_t_histogram;

    std::map<u64, u64>
        m5_margin_histogram;

    std::map<u64, u64>
        m7_margin_histogram;

    std::vector<ExceptionalPoint>
        first_exceptions;

    first_exceptions.reserve(10);

    /*
     * A run is consecutive K values for the same prime
     * where the global winner remains m=5 or m=7.
     */
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

        int previous_exceptional_m = 0;
        u64 current_run_length = 0;

        for (
            u64 K = 1;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            if (K > 1) {
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
            }

            const Witness global =
                global_winner(best_by_m);

            if (!allowed_m(global.m)) {
                ++exceptional_points;

                const Witness allowed =
                    allowed_winner(best_by_m);

                const i128 global_value =
                    arithmetic_value(
                        r,
                        global
                    );

                const i128 allowed_value =
                    arithmetic_value(
                        r,
                        allowed
                    );

                const i128 global_score_num =
                    static_cast<i128>(
                        global.k
                    ) *
                    static_cast<i128>(
                        allowed.m
                    );

                const i128 allowed_score_num =
                    static_cast<i128>(
                        allowed.k
                    ) *
                    static_cast<i128>(
                        global.m
                    );

                const i128 margin =
                    global_score_num -
                    allowed_score_num;

                if (margin > 0) {
                    if (global.m == 5) {
                        ++m5_allowed_margin_positive;
                        m5_allowed_margin_sum +=
                            static_cast<u64>(margin);
                        m5_margin_histogram[
                            static_cast<u64>(margin)
                        ]++;
                    } else if (global.m == 7) {
                        ++m7_allowed_margin_positive;
                        m7_allowed_margin_sum +=
                            static_cast<u64>(margin);
                        m7_margin_histogram[
                            static_cast<u64>(margin)
                        ]++;
                    }
                }

                if (
                    global.m == 5
                ) {
                    ++m5_points;

                    if (
                        previous_exceptional_m != 5
                    ) {
                        ++m5_run_count;
                        current_run_length = 1;
                    } else {
                        ++current_run_length;
                    }

                    previous_exceptional_m = 5;

                    m5_K_histogram[K]++;
                    m5_k_histogram[global.k]++;
                    m5_t_histogram[global.t]++;

                    m5_first_K_sum += K;

                    if (global.k == global.t) {
                        ++m5_k_equals_t;
                    } else if (global.k < global.t) {
                        ++m5_k_less_t;
                    } else {
                        ++m5_k_greater_t;
                    }

                    if (global.sign > 0) {
                        ++m5_sign_plus;
                    } else {
                        ++m5_sign_minus;
                    }

                    if (
                        gcd_u64(
                            global.k,
                            static_cast<u64>(global.m) *
                            r +
                            static_cast<u64>(
                                global.sign > 0 ? 1 : 0
                            )
                        ) > 1
                    ) {
                        ++m5_gcd_kr;
                    }

                    if (
                        first_exceptions.size() < 10
                    ) {
                        ExceptionalPoint point;
                        point.prime = r;
                        point.K = K;
                        point.winner = global;
                        point.allowed_winner = allowed;

                        first_exceptions.push_back(point);
                    }
                } else if (
                    global.m == 7
                ) {
                    ++m7_points;

                    if (
                        previous_exceptional_m != 7
                    ) {
                        ++m7_run_count;
                        current_run_length = 1;
                    } else {
                        ++current_run_length;
                    }

                    previous_exceptional_m = 7;

                    m7_K_histogram[K]++;
                    m7_k_histogram[global.k]++;
                    m7_t_histogram[global.t]++;

                    m7_first_K_sum += K;

                    if (global.k == global.t) {
                        ++m7_k_equals_t;
                    } else if (global.k < global.t) {
                        ++m7_k_less_t;
                    } else {
                        ++m7_k_greater_t;
                    }

                    if (global.sign > 0) {
                        ++m7_sign_plus;
                    } else {
                        ++m7_sign_minus;
                    }

                    if (
                        gcd_u64(
                            global.k,
                            static_cast<u64>(global.m) *
                            r +
                            static_cast<u64>(
                                global.sign > 0 ? 1 : 0
                            )
                        ) > 1
                    ) {
                        ++m7_gcd_kr;
                    }

                    if (
                        first_exceptions.size() < 10
                    ) {
                        ExceptionalPoint point;
                        point.prime = r;
                        point.K = K;
                        point.winner = global;
                        point.allowed_winner = allowed;

                        first_exceptions.push_back(point);
                    }
                }
            } else {
                if (
                    previous_exceptional_m == 5 &&
                    current_run_length > max_m5_run
                ) {
                    max_m5_run =
                        current_run_length;
                }

                if (
                    previous_exceptional_m == 7 &&
                    current_run_length > max_m7_run
                ) {
                    max_m7_run =
                        current_run_length;
                }

                previous_exceptional_m = 0;
                current_run_length = 0;
            }
        }

        if (
            previous_exceptional_m == 5 &&
            current_run_length > max_m5_run
        ) {
            max_m5_run =
                current_run_length;
        }

        if (
            previous_exceptional_m == 7 &&
            current_run_length > max_m7_run
        ) {
            max_m7_run =
                current_run_length;
        }

        if (m5_K_histogram.empty()) {
            // Nothing.
        }

        /*
         * Count whether this prime ever contains
         * an m=5 or m=7 winner.
         */
        bool has_m5 = false;
        bool has_m7 = false;

        for (const auto& entry : m5_K_histogram) {
            (void)entry;
            has_m5 = true;
            break;
        }

        for (const auto& entry : m7_K_histogram) {
            (void)entry;
            has_m7 = true;
            break;
        }

        if (has_m5) {
            ++m5_prime_count;
        }

        if (has_m7) {
            ++m7_prime_count;
        }
    }

    /*
     * The following two counters are run counts by multiplier.
     * Keep the values above, but the precise multi-point
     * classification is reconstructed from the observed K
     * histogram below.
     */
    std::map<u64, u64> dummy;
    (void)dummy;

    std::cout
        << "\nTOTAL_EXCEPTIONAL_POINTS="
        << exceptional_points
        << "\n";

    std::cout
        << "M5_POINTS="
        << m5_points
        << "\n";

    std::cout
        << "M7_POINTS="
        << m7_points
        << "\n";

    std::cout
        << "M5_PRIME_COUNT="
        << m5_prime_count
        << "\n";

    std::cout
        << "M7_PRIME_COUNT="
        << m7_prime_count
        << "\n";

    std::cout
        << "\nRUNS\n";

    std::cout
        << "M5_RUN_COUNT="
        << m5_run_count
        << "\n";

    std::cout
        << "M7_RUN_COUNT="
        << m7_run_count
        << "\n";

    std::cout
        << "M5_MAX_RUN="
        << max_m5_run
        << "\n";

    std::cout
        << "M7_MAX_RUN="
        << max_m7_run
        << "\n";

    std::cout
        << "\nM5_STRUCTURE\n";

    std::cout
        << "M5_K_LESS_T="
        << m5_k_less_t
        << "\n";

    std::cout
        << "M5_K_EQUAL_T="
        << m5_k_equals_t
        << "\n";

    std::cout
        << "M5_K_GREATER_T="
        << m5_k_greater_t
        << "\n";

    std::cout
        << "M5_SIGN_PLUS="
        << m5_sign_plus
        << "\n";

    std::cout
        << "M5_SIGN_MINUS="
        << m5_sign_minus
        << "\n";

    std::cout
        << "M5_GCD_K_VALUE_GT1="
        << m5_gcd_kr
        << "\n";

    std::cout
        << "\nM7_STRUCTURE\n";

    std::cout
        << "M7_K_LESS_T="
        << m7_k_less_t
        << "\n";

    std::cout
        << "M7_K_EQUAL_T="
        << m7_k_equals_t
        << "\n";

    std::cout
        << "M7_K_GREATER_T="
        << m7_k_greater_t
        << "\n";

    std::cout
        << "M7_SIGN_PLUS="
        << m7_sign_plus
        << "\n";

    std::cout
        << "M7_SIGN_MINUS="
        << m7_sign_minus
        << "\n";

    std::cout
        << "M7_GCD_K_VALUE_GT1="
        << m7_gcd_kr
        << "\n";

    std::cout
        << "\nMARGIN\n";

    std::cout
        << "M5_MARGIN_POSITIVE="
        << m5_allowed_margin_positive
        << "\n";

    std::cout
        << "M7_MARGIN_POSITIVE="
        << m7_allowed_margin_positive
        << "\n";

    std::cout
        << "M5_MARGIN_SUM="
        << m5_allowed_margin_sum
        << "\n";

    std::cout
        << "M7_MARGIN_SUM="
        << m7_allowed_margin_sum
        << "\n";

    std::cout
        << "\nM5_K_HISTOGRAM\n";

    for (
        const auto& entry :
        m5_K_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM7_K_HISTOGRAM\n";

    for (
        const auto& entry :
        m7_K_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM5_K_VALUE_HISTOGRAM\n";

    for (
        const auto& entry :
        m5_k_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM7_K_VALUE_HISTOGRAM\n";

    for (
        const auto& entry :
        m7_k_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM5_T_HISTOGRAM\n";

    for (
        const auto& entry :
        m5_t_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM7_T_HISTOGRAM\n";

    for (
        const auto& entry :
        m7_t_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM5_MARGIN_HISTOGRAM\n";

    for (
        const auto& entry :
        m5_margin_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nM7_MARGIN_HISTOGRAM\n";

    for (
        const auto& entry :
        m7_margin_histogram
    ) {
        std::cout
            << entry.first
            << " "
            << entry.second
            << "\n";
    }

    std::cout
        << "\nFIRST_EXCEPTIONAL_POINTS\n";

    for (
        std::size_t i = 0;
        i < first_exceptions.size();
        ++i
    ) {
        const ExceptionalPoint& point =
            first_exceptions[i];

        std::cout
            << "POINT="
            << (i + 1)
            << " PRIME="
            << point.prime
            << " K="
            << point.K
            << "\n";

        print_witness(
            point.prime,
            "WINNER",
            point.winner
        );

        std::cout << "\n";

        print_witness(
            point.prime,
            "ALLOWED",
            point.allowed_winner
        );

        std::cout << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
