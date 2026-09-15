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

struct Factor {
    u64 p = 0;
    int exponent = 0;
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
    std::vector<bool> is_prime(
        static_cast<std::size_t>(limit) + 1,
        true
    );

    if (limit >= 0) {
        is_prime[0] = false;
    }

    if (limit >= 1) {
        is_prime[1] = false;
    }

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

static bool same_witness(
    const Witness& a,
    const Witness& b
) {
    return
        a.m == b.m &&
        a.k == b.k &&
        a.t == b.t &&
        a.sign == b.sign;
}

static u64 arithmetic_value(
    u64 r,
    const Witness& w
) {
    const i128 value =
        static_cast<i128>(w.m) *
        static_cast<i128>(r) +
        static_cast<i128>(w.sign);

    return static_cast<u64>(value);
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

static std::vector<Factor> factorize(
    u64 n
) {
    std::vector<Factor> result;

    if (n < 2) {
        return result;
    }

    for (
        u64 p = 2;
        p <= n / p;
        ++p
    ) {
        if (n % p != 0) {
            continue;
        }

        int exponent = 0;

        while (n % p == 0) {
            n /= p;
            ++exponent;
        }

        result.push_back({
            p,
            exponent
        });
    }

    if (n > 1) {
        result.push_back({
            n,
            1
        });
    }

    return result;
}

static void print_factorization(
    const std::vector<Factor>& factors
) {
    if (factors.empty()) {
        std::cout << "1";
        return;
    }

    bool first = true;

    for (const Factor& factor : factors) {
        if (!first) {
            std::cout << "*";
        }

        first = false;

        std::cout << factor.p;

        if (factor.exponent > 1) {
            std::cout
                << "^"
                << factor.exponent;
        }
    }
}

static void print_witness(
    u64 r,
    const std::string& label,
    const Witness& w
) {
    std::cout
        << label
        << "=(m="
        << w.m
        << ",k="
        << w.k
        << ",t="
        << w.t
        << ",sign="
        << (w.sign > 0 ? "+1" : "-1")
        << ",value="
        << arithmetic_value(r, w)
        << ",norm="
        << w.k
        << "/"
        << w.m
        << ")";
}

static i128 normalized_cross_gap(
    const Witness& a,
    const Witness& b
) {
    return
        static_cast<i128>(a.k) *
        static_cast<i128>(b.m) -
        static_cast<i128>(b.k) *
        static_cast<i128>(a.m);
}

static void print_allowed_comparisons(
    const Witness& global,
    const std::vector<Witness>& best_by_m
) {
    std::cout
        << "ALLOWED_COMPARISONS\n";

    for (int m : {1, 2, 3, 4, 6}) {
        const Witness& competitor =
            best_by_m[
                static_cast<std::size_t>(m)
            ];

        const i128 cross =
            normalized_cross_gap(
                global,
                competitor
            );

        std::cout
            << "M="
            << m
            << " CROSS="
            << i128_to_string(cross)
            << " K="
            << competitor.k
            << " T="
            << competitor.t
            << " SIGN="
            << (
                competitor.sign > 0
                    ? "+1"
                    : "-1"
            )
            << "\n";
    }
}

int main() {
    constexpr int EXPERIMENT = 440;
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
        generate_primes(
            PRIME_LIMIT
        );

    std::cout
        << "PRIME_COUNT="
        << primes.size()
        << "\n";

    u64 m5_events = 0;
    u64 m7_events = 0;

    u64 m5_prime_count = 0;
    u64 m7_prime_count = 0;

    for (u64 r : primes) {
        std::vector<Witness> best_by_m(
            static_cast<std::size_t>(M_MAX) + 1
        );

        /*
         * Initial divisor record k=1.
         */
        for (int m = 1; m <= M_MAX; ++m) {
            best_by_m[
                static_cast<std::size_t>(m)
            ] =
                make_witness(
                    r,
                    m,
                    1,
                    +1
                );
        }

        bool prime_has_m5_event = false;
        bool prime_has_m7_event = false;

        Witness previous_m5;
        Witness previous_m7;

        bool have_previous_m5 = false;
        bool have_previous_m7 = false;

        for (
            u64 K = 2;
            K <= static_cast<u64>(K_LIMIT);
            ++K
        ) {
            /*
             * Check whether K creates a new divisor record
             * for each multiplier.
             */
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

            const Witness global =
                global_winner(
                    best_by_m
                );

            /*
             * We only care about m=5 and m=7 global winners.
             */
            if (
                global.m != 5 &&
                global.m != 7
            ) {
                continue;
            }

            bool new_exceptional_record = false;

            if (global.m == 5) {
                if (
                    !have_previous_m5 ||
                    !same_witness(
                        global,
                        previous_m5
                    )
                ) {
                    previous_m5 = global;
                    have_previous_m5 = true;
                    new_exceptional_record = true;
                }
            } else {
                if (
                    !have_previous_m7 ||
                    !same_witness(
                        global,
                        previous_m7
                    )
                ) {
                    previous_m7 = global;
                    have_previous_m7 = true;
                    new_exceptional_record = true;
                }
            }

            if (!new_exceptional_record) {
                continue;
            }

            if (global.m == 5) {
                ++m5_events;

                if (!prime_has_m5_event) {
                    ++m5_prime_count;
                    prime_has_m5_event = true;
                }
            } else {
                ++m7_events;

                if (!prime_has_m7_event) {
                    ++m7_prime_count;
                    prime_has_m7_event = true;
                }
            }

            std::cout
                << "\n"
                << (
                    global.m == 5
                        ? "M5_EVENT"
                        : "M7_EVENT"
                )
                << "\n";

            std::cout
                << "PRIME="
                << r
                << "\n";

            std::cout
                << "K="
                << K
                << "\n";

            print_witness(
                r,
                "WINNER",
                global
            );

            std::cout
                << "\n";

            const u64 A =
                arithmetic_value(
                    r,
                    global
                );

            std::cout
                << "A="
                << A
                << "\n";

            std::cout
                << "A_CHECK="
                << global.k
                << "*"
                << global.t
                << "="
                << global.k * global.t
                << "\n";

            std::cout
                << "FACTOR_A=";

            print_factorization(
                factorize(A)
            );

            std::cout
                << "\n";

            std::cout
                << "FACTOR_K=";

            print_factorization(
                factorize(global.k)
            );

            std::cout
                << "\n";

            std::cout
                << "FACTOR_T=";

            print_factorization(
                factorize(global.t)
            );

            std::cout
                << "\n";

            std::cout
                << "GCD_K_T="
                << gcd_u64(
                    global.k,
                    global.t
                )
                << "\n";

            std::cout
                << "GCD_K_R="
                << gcd_u64(
                    global.k,
                    r
                )
                << "\n";

            std::cout
                << "GCD_T_R="
                << gcd_u64(
                    global.t,
                    r
                )
                << "\n";

            const Witness allowed =
                allowed_winner(
                    best_by_m
                );

            std::cout
                << "STRONGEST_ALLOWED_M="
                << allowed.m
                << "\n";

            print_witness(
                r,
                "STRONGEST_ALLOWED",
                allowed
            );

            std::cout
                << "\n";

            const i128 gap =
                normalized_cross_gap(
                    global,
                    allowed
                );

            std::cout
                << "NORMALIZED_CROSS_GAP="
                << i128_to_string(gap)
                << "\n";

            print_allowed_comparisons(
                global,
                best_by_m
            );

            std::cout
                << "M5_STREAM\n";

            print_witness(
                r,
                "M5",
                best_by_m[5]
            );

            std::cout
                << "\n";

            std::cout
                << "M7_STREAM\n";

            print_witness(
                r,
                "M7",
                best_by_m[7]
            );

            std::cout
                << "\n";

            std::cout
                << "R_MOD_K="
                << (r % global.k)
                << "\n";

            std::cout
                << "R_MOD_T="
                << (r % global.t)
                << "\n";

            std::cout
                << "K_MOD_T="
                << (global.k % global.t)
                << "\n";

            std::cout
                << "T_MOD_K="
                << (global.t % global.k)
                << "\n";
        }
    }

    std::cout
        << "\nSUMMARY\n";

    std::cout
        << "M5_EVENTS="
        << m5_events
        << "\n";

    std::cout
        << "M7_EVENTS="
        << m7_events
        << "\n";

    std::cout
        << "M5_EVENT_PRIME_COUNT="
        << m5_prime_count
        << "\n";

    std::cout
        << "M7_EVENT_PRIME_COUNT="
        << m7_prime_count
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}