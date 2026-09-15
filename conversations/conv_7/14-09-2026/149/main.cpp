#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Factor {
    u64 p = 0;
    int exponent = 0;
};

struct Witness {
    int m = 1;
    u64 k = 1;
    u64 t = 0;
    int sign = +1;
};

struct ARecord {
    u64 r = 0;
    int sign = +1;
    u64 A = 0;

    u64 divisor_count = 0;
    u64 distinct_prime_count = 0;
    u64 largest_prime_factor = 0;
    u64 squarefree_kernel = 1;

    u64 divisor_record_count = 0;
    u64 exceptional_record_count = 0;

    u64 max_record_k = 0;
    u64 max_record_t = 0;

    u64 max_margin = 0;
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

    sieve[0] = false;
    sieve[1] = false;

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

    for (int n = 2; n <= limit; ++n) {
        if (sieve[n]) {
            primes.push_back(
                static_cast<u64>(n)
            );
        }
    }

    return primes;
}

static std::vector<Factor> factorize(u64 n) {
    std::vector<Factor> factors;

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

        factors.push_back({
            p,
            exponent
        });
    }

    if (n > 1) {
        factors.push_back({
            n,
            1
        });
    }

    return factors;
}

static u64 divisor_count(
    const std::vector<Factor>& factors
) {
    u64 result = 1;

    for (const Factor& f : factors) {
        result *=
            static_cast<u64>(f.exponent + 1);
    }

    return result;
}

static u64 squarefree_kernel(
    const std::vector<Factor>& factors
) {
    u64 result = 1;

    for (const Factor& f : factors) {
        result *= f.p;
    }

    return result;
}

static u64 largest_prime_factor(
    const std::vector<Factor>& factors
) {
    if (factors.empty()) {
        return 1;
    }

    return factors.back().p;
}

static bool divides(
    u64 A,
    u64 k
) {
    return k != 0 && A % k == 0;
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

static Witness allowed_winner(
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

static i128 normalized_margin(
    const Witness& exceptional,
    const Witness& allowed
) {
    return
        static_cast<i128>(exceptional.k) *
        static_cast<i128>(allowed.m) -
        static_cast<i128>(allowed.k) *
        static_cast<i128>(exceptional.m);
}

static void print_factorization(
    const std::vector<Factor>& factors
) {
    bool first = true;

    for (const Factor& f : factors) {
        if (!first) {
            std::cout << "*";
        }

        first = false;

        std::cout << f.p;

        if (f.exponent > 1) {
            std::cout
                << "^"
                << f.exponent;
        }
    }

    if (first) {
        std::cout << "1";
    }
}

static bool factor_divisor_test(
    u64 A,
    int m,
    u64 k
) {
    return
        divides(
            A,
            k
        ) &&
        A / k >= 1 &&
        m >= 1;
}

static ARecord analyze_A(
    u64 r,
    int sign,
    u64 A,
    int K_LIMIT
) {
    ARecord result;

    result.r = r;
    result.sign = sign;
    result.A = A;

    const std::vector<Factor> factors =
        factorize(A);

    result.divisor_count =
        divisor_count(factors);

    result.distinct_prime_count =
        static_cast<u64>(factors.size());

    result.largest_prime_factor =
        largest_prime_factor(factors);

    result.squarefree_kernel =
        squarefree_kernel(factors);

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

    u64 last_k = 1;

    for (
        u64 k = 2;
        k <= static_cast<u64>(K_LIMIT);
        ++k
    ) {
        if (!factor_divisor_test(A, 5, k)) {
            continue;
        }

        ++result.divisor_record_count;

        const u64 t = A / k;

        const Witness exceptional =
            make_witness(
                r,
                5,
                k,
                sign
            );

        /*
         * Reconstruct the allowed divisor records
         * at the same r.
         */
        for (int m : {1, 2, 3, 4, 6}) {
            Witness current =
                best[
                    static_cast<std::size_t>(m)
                ];

            /*
             * Search only up to k. This is enough to determine
             * whether this divisor produces a new record.
             */
            for (u64 d = current.k + 1; d <= k; ++d) {
                if (d > K_LIMIT) {
                    break;
                }

                if (
                    (static_cast<i128>(m) *
                     static_cast<i128>(r) + 1) %
                    static_cast<i128>(d) == 0
                ) {
                    current =
                        make_witness(
                            r,
                            m,
                            d,
                            +1
                        );
                }

                if (
                    (static_cast<i128>(m) *
                     static_cast<i128>(r) - 1) %
                    static_cast<i128>(d) == 0
                ) {
                    const Witness candidate =
                        make_witness(
                            r,
                            m,
                            d,
                            -1
                        );

                    if (
                        candidate.k >
                        current.k
                    ) {
                        current = candidate;
                    }
                }
            }

            best[
                static_cast<std::size_t>(m)
            ] = current;
        }

        if (k < last_k) {
            continue;
        }

        last_k = k;

        const Witness allowed =
            allowed_winner(best);

        if (
            better_normalized(
                exceptional,
                allowed
            )
        ) {
            ++result.exceptional_record_count;

            const i128 margin =
                normalized_margin(
                    exceptional,
                    allowed
                );

            if (margin > 0) {
                const u64 margin_u64 =
                    static_cast<u64>(margin);

                if (
                    margin_u64 >
                    result.max_margin
                ) {
                    result.max_margin =
                        margin_u64;
                }
            }
        }

        if (k >= result.max_record_k) {
            result.max_record_k = k;
            result.max_record_t = t;
        }
    }

    return result;
}

static void print_record(
    const ARecord& record
) {
    std::cout
        << "R="
        << record.r
        << " SIGN="
        << (
            record.sign > 0
                ? "+1"
                : "-1"
        )
        << " A="
        << record.A
        << " FACTOR=";

    print_factorization(
        factorize(record.A)
    );

    std::cout
        << " TAU="
        << record.divisor_count
        << " OMEGA="
        << record.distinct_prime_count
        << " LARGEST_P="
        << record.largest_prime_factor
        << " RAD="
        << record.squarefree_kernel
        << " RECORDS="
        << record.divisor_record_count
        << " EXCEPTIONAL_RECORDS="
        << record.exceptional_record_count
        << " MAX_K="
        << record.max_record_k
        << " MAX_T="
        << record.max_record_t
        << " MAX_MARGIN="
        << record.max_margin
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 441;
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

    u64 m5_plus_exceptional = 0;
    u64 m5_minus_exceptional = 0;

    u64 total_exceptional_records = 0;

    u64 tau_sum_exceptional = 0;
    u64 tau_sum_normal = 0;

    u64 max_exceptional_tau = 0;
    u64 max_normal_tau = 0;

    u64 max_exceptional_omega = 0;
    u64 max_normal_omega = 0;

    std::vector<ARecord> exceptional_records;

    for (u64 r : primes) {
        for (int sign : {-1, +1}) {
            const i128 value =
                static_cast<i128>(5) *
                static_cast<i128>(r) +
                static_cast<i128>(sign);

            const u64 A =
                static_cast<u64>(value);

            ++total_A;

            const ARecord record =
                analyze_A(
                    r,
                    sign,
                    A,
                    K_LIMIT
                );

            if (
                record.exceptional_record_count >
                0
            ) {
                ++exceptional_A;

                total_exceptional_records +=
                    record.exceptional_record_count;

                tau_sum_exceptional +=
                    record.divisor_count;

                if (
                    record.divisor_count >
                    max_exceptional_tau
                ) {
                    max_exceptional_tau =
                        record.divisor_count;
                }

                if (
                    record.distinct_prime_count >
                    max_exceptional_omega
                ) {
                    max_exceptional_omega =
                        record.distinct_prime_count;
                }

                if (sign > 0) {
                    ++m5_plus_exceptional;
                } else {
                    ++m5_minus_exceptional;
                }

                exceptional_records.push_back(
                    record
                );
            } else {
                tau_sum_normal +=
                    record.divisor_count;

                if (
                    record.divisor_count >
                    max_normal_tau
                ) {
                    max_normal_tau =
                        record.divisor_count;
                }

                if (
                    record.distinct_prime_count >
                    max_normal_omega
                ) {
                    max_normal_omega =
                        record.distinct_prime_count;
                }
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
        << "NON_EXCEPTIONAL_A="
        << (total_A - exceptional_A)
        << "\n";

    std::cout
        << "M5_PLUS_EXCEPTIONAL="
        << m5_plus_exceptional
        << "\n";

    std::cout
        << "M5_MINUS_EXCEPTIONAL="
        << m5_minus_exceptional
        << "\n";

    std::cout
        << "TOTAL_EXCEPTIONAL_RECORDS="
        << total_exceptional_records
        << "\n";

    std::cout
        << "MAX_EXCEPTIONAL_TAU="
        << max_exceptional_tau
        << "\n";

    std::cout
        << "MAX_NORMAL_TAU="
        << max_normal_tau
        << "\n";

    std::cout
        << "MAX_EXCEPTIONAL_OMEGA="
        << max_exceptional_omega
        << "\n";

    std::cout
        << "MAX_NORMAL_OMEGA="
        << max_normal_omega
        << "\n";

    std::cout
        << "\nEXCEPTIONAL_A_RECORDS\n";

    for (
        const ARecord& record :
        exceptional_records
    ) {
        print_record(record);
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
