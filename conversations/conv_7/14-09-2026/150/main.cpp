#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Factor {
    u64 p;
    int exponent;
};

struct Witness {
    int m = 1;
    u64 k = 1;
    u64 t = 0;
    int sign = +1;
};

struct Geometry {
    u64 A = 0;
    u64 sqrt_floor = 0;

    u64 d_below = 0;
    u64 d_above = 0;

    u64 gap = 0;
    u64 product_error = 0;

    u64 divisors_le_3000 = 0;

    u64 window_90_count = 0;
    u64 window_95_count = 0;
    u64 window_99_count = 0;

    u64 nearest_divisor_distance = 0;

    double relative_gap = 0.0;
    double ratio = 0.0;
};

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

static void generate_divisors_recursive(
    const std::vector<Factor>& factors,
    std::size_t index,
    u64 current,
    std::vector<u64>& divisors
) {
    if (index == factors.size()) {
        divisors.push_back(current);
        return;
    }

    const Factor& factor = factors[index];

    u64 value = current;

    for (
        int e = 0;
        e <= factor.exponent;
        ++e
    ) {
        generate_divisors_recursive(
            factors,
            index + 1,
            value,
            divisors
        );

        value *= factor.p;
    }
}

static std::vector<u64> generate_divisors(
    u64 n
) {
    const std::vector<Factor> factors =
        factorize(n);

    std::vector<u64> divisors;

    generate_divisors_recursive(
        factors,
        0,
        1,
        divisors
    );

    return divisors;
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

static u64 isqrt_u64(u64 n) {
    u64 x =
        static_cast<u64>(
            std::sqrt(
                static_cast<long double>(n)
            )
        );

    while (
        (x + 1) <= n / (x + 1)
    ) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static Geometry analyze_geometry(
    u64 A,
    const std::vector<u64>& divisors
) {
    Geometry g;

    g.A = A;

    const long double root =
        std::sqrt(
            static_cast<long double>(A)
        );

    g.sqrt_floor =
        isqrt_u64(A);

    u64 below = 1;
    u64 above = 0;

    for (u64 d : divisors) {
        if (d > 3000) {
            continue;
        }

        ++g.divisors_le_3000;

        const long double ratio =
            static_cast<long double>(d) /
            root;

        if (d <= g.sqrt_floor) {
            if (d > below) {
                below = d;
            }
        }

        if (
            d * d >= A &&
            (above == 0 || d < above)
        ) {
            above = d;
        }

        if (ratio >= 0.90L &&
            ratio <= 1.10L) {
            ++g.window_90_count;
        }

        if (ratio >= 0.95L &&
            ratio <= 1.05L) {
            ++g.window_95_count;
        }

        if (ratio >= 0.99L &&
            ratio <= 1.01L) {
            ++g.window_99_count;
        }
    }

    g.d_below = below;
    g.d_above = above;

    if (above != 0) {
        g.gap = above - below;

        const i128 product =
            static_cast<i128>(below) *
            static_cast<i128>(above);

        const i128 error =
            product -
            static_cast<i128>(A);

        if (error < 0) {
            g.product_error =
                static_cast<u64>(-error);
        } else {
            g.product_error =
                static_cast<u64>(error);
        }

        g.relative_gap =
            static_cast<double>(g.gap) /
            static_cast<double>(root);

        g.ratio =
            static_cast<double>(above) /
            static_cast<double>(below);

        const u64 below_distance =
            g.sqrt_floor >= below
                ? g.sqrt_floor - below
                : below - g.sqrt_floor;

        const u64 above_distance =
            above >= g.sqrt_floor
                ? above - g.sqrt_floor
                : g.sqrt_floor - above;

        g.nearest_divisor_distance =
            below_distance < above_distance
                ? below_distance
                : above_distance;
    }

    return g;
}

static i128 normalized_gap(
    const Witness& exceptional,
    const Witness& allowed
) {
    return
        static_cast<i128>(exceptional.k) *
        static_cast<i128>(allowed.m) -
        static_cast<i128>(allowed.k) *
        static_cast<i128>(exceptional.m);
}

static std::string factorization_string(
    const std::vector<Factor>& factors
) {
    std::string result;

    bool first = true;

    for (const Factor& f : factors) {
        if (!first) {
            result += "*";
        }

        first = false;

        result +=
            std::to_string(f.p);

        if (f.exponent > 1) {
            result += "^";
            result +=
                std::to_string(f.exponent);
        }
    }

    if (result.empty()) {
        result = "1";
    }

    return result;
}

static void print_geometry(
    const Geometry& g
) {
    std::cout
        << "A="
        << g.A
        << " SQRT_FLOOR="
        << g.sqrt_floor
        << " BELOW="
        << g.d_below
        << " ABOVE="
        << g.d_above
        << " GAP="
        << g.gap
        << " PRODUCT_ERROR="
        << g.product_error
        << " REL_GAP="
        << g.relative_gap
        << " RATIO="
        << g.ratio
        << " DIVISORS="
        << g.divisors_le_3000
        << " W90="
        << g.window_90_count
        << " W95="
        << g.window_95_count
        << " W99="
        << g.window_99_count
        << " NEAREST_DIST="
        << g.nearest_divisor_distance
        << "\n";
}

int main() {
    constexpr int EXPERIMENT = 442;
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

    u64 central_gap_sum_exceptional = 0;
    u64 central_gap_sum_normal = 0;

    u64 window90_exceptional = 0;
    u64 window90_normal = 0;

    u64 window95_exceptional = 0;
    u64 window95_normal = 0;

    u64 window99_exceptional = 0;
    u64 window99_normal = 0;

    u64 nearest_distance_sum_exceptional = 0;
    u64 nearest_distance_sum_normal = 0;

    u64 exceptional_first_k_below_sqrt = 0;
    u64 exceptional_first_k_above_sqrt = 0;

    u64 exceptional_first_k_distance_sum = 0;

    u64 exceptional_record_count = 0;

    u64 exceptional_close_pair_10 = 0;
    u64 exceptional_close_pair_5 = 0;
    u64 exceptional_close_pair_1 = 0;

    u64 normal_close_pair_10 = 0;
    u64 normal_close_pair_5 = 0;
    u64 normal_close_pair_1 = 0;

    u64 largest_exceptional_gap = 0;
    u64 largest_normal_gap = 0;

    std::vector<std::string>
        exceptional_output;

    for (u64 r : primes) {
        for (int sign : {-1, +1}) {
            const i128 value =
                static_cast<i128>(5) *
                static_cast<i128>(r) +
                static_cast<i128>(sign);

            const u64 A =
                static_cast<u64>(value);

            ++total_A;

            const std::vector<u64> divisors =
                generate_divisors(A);

            Geometry geometry =
                analyze_geometry(
                    A,
                    divisors
                );

            /*
             * Build the allowed envelope up to K_LIMIT.
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

            std::vector<Witness>
                allowed_at_K(
                    static_cast<std::size_t>(
                        K_LIMIT + 1
                    )
                );

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
                            best[mi].k ||
                            (
                                candidate.k ==
                                best[mi].k &&
                                candidate.sign >
                                best[mi].sign
                            )
                        ) {
                            best[mi] =
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
                            best[mi].k ||
                            (
                                candidate.k ==
                                best[mi].k &&
                                candidate.sign >
                                best[mi].sign
                            )
                        ) {
                            best[mi] =
                                candidate;
                        }
                    }
                }

                allowed_at_K[
                    static_cast<std::size_t>(K)
                ] =
                    allowed_winner(best);
            }

            u64 first_exception_k = 0;
            u64 first_exception_t = 0;
            u64 first_exception_margin = 0;

            bool A_is_exceptional = false;

            /*
             * Inspect every divisor k of A up to K_LIMIT.
             */
            for (u64 k : divisors) {
                if (k > K_LIMIT) {
                    continue;
                }

                const u64 t =
                    A / k;

                const Witness exceptional =
                    make_witness(
                        r,
                        5,
                        k,
                        sign
                    );

                const Witness& allowed =
                    allowed_at_K[
                        static_cast<std::size_t>(k)
                    ];

                if (
                    !better_normalized(
                        exceptional,
                        allowed
                    )
                ) {
                    continue;
                }

                ++exceptional_record_count;

                const i128 margin =
                    normalized_gap(
                        exceptional,
                        allowed
                    );

                if (!A_is_exceptional) {
                    A_is_exceptional = true;

                    first_exception_k = k;
                    first_exception_t = t;

                    if (margin > 0) {
                        first_exception_margin =
                            static_cast<u64>(
                                margin
                            );
                    }

                    const u64 sqrt_floor =
                        geometry.sqrt_floor;

                    if (
                        first_exception_k <=
                        sqrt_floor
                    ) {
                        ++exceptional_first_k_below_sqrt;
                    } else {
                        ++exceptional_first_k_above_sqrt;
                    }

                    const u64 distance =
                        first_exception_k >=
                        sqrt_floor
                            ? first_exception_k -
                              sqrt_floor
                            : sqrt_floor -
                              first_exception_k;

                    exceptional_first_k_distance_sum +=
                        distance;
                }
            }

            if (A_is_exceptional) {
                ++exceptional_A;

                central_gap_sum_exceptional +=
                    geometry.gap;

                nearest_distance_sum_exceptional +=
                    geometry.nearest_divisor_distance;

                window90_exceptional +=
                    geometry.window_90_count;

                window95_exceptional +=
                    geometry.window_95_count;

                window99_exceptional +=
                    geometry.window_99_count;

                if (geometry.gap >
                    largest_exceptional_gap) {
                    largest_exceptional_gap =
                        geometry.gap;
                }

                if (
                    geometry.relative_gap <= 0.10
                ) {
                    ++exceptional_close_pair_10;
                }

                if (
                    geometry.relative_gap <= 0.05
                ) {
                    ++exceptional_close_pair_5;
                }

                if (
                    geometry.relative_gap <= 0.01
                ) {
                    ++exceptional_close_pair_1;
                }

                std::string line;

                line +=
                    "EXCEPTION R=" +
                    std::to_string(r);

                line +=
                    " SIGN=" +
                    std::string(
                        sign > 0 ? "+1" : "-1"
                    );

                line +=
                    " A=" +
                    std::to_string(A);

                line +=
                    " FACTOR=" +
                    factorization_string(
                        factorize(A)
                    );

                line +=
                    " FIRST_K=" +
                    std::to_string(
                        first_exception_k
                    );

                line +=
                    " FIRST_T=" +
                    std::to_string(
                        first_exception_t
                    );

                line +=
                    " FIRST_MARGIN=" +
                    std::to_string(
                        first_exception_margin
                    );

                line +=
                    " BELOW=" +
                    std::to_string(
                        geometry.d_below
                    );

                line +=
                    " ABOVE=" +
                    std::to_string(
                        geometry.d_above
                    );

                line +=
                    " GAP=" +
                    std::to_string(
                        geometry.gap
                    );

                line +=
                    " REL_GAP=" +
                    std::to_string(
                        geometry.relative_gap
                    );

                line +=
                    " RATIO=" +
                    std::to_string(
                        geometry.ratio
                    );

                line +=
                    " W90=" +
                    std::to_string(
                        geometry.window_90_count
                    );

                line +=
                    " W95=" +
                    std::to_string(
                        geometry.window_95_count
                    );

                line +=
                    " W99=" +
                    std::to_string(
                        geometry.window_99_count
                    );

                exceptional_output.push_back(line);
            } else {
                central_gap_sum_normal +=
                    geometry.gap;

                nearest_distance_sum_normal +=
                    geometry.nearest_divisor_distance;

                window90_normal +=
                    geometry.window_90_count;

                window95_normal +=
                    geometry.window_95_count;

                window99_normal +=
                    geometry.window_99_count;

                if (geometry.gap >
                    largest_normal_gap) {
                    largest_normal_gap =
                        geometry.gap;
                }

                if (
                    geometry.relative_gap <= 0.10
                ) {
                    ++normal_close_pair_10;
                }

                if (
                    geometry.relative_gap <= 0.05
                ) {
                    ++normal_close_pair_5;
                }

                if (
                    geometry.relative_gap <= 0.01
                ) {
                    ++normal_close_pair_1;
                }
            }
        }
    }

    const u64 normal_A =
        total_A - exceptional_A;

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
        << "NORMAL_A="
        << normal_A
        << "\n";

    std::cout
        << "EXCEPTIONAL_RECORD_COUNT="
        << exceptional_record_count
        << "\n";

    std::cout
        << "\nCENTRAL_GAP\n";

    std::cout
        << "EXCEPTIONAL_SUM="
        << central_gap_sum_exceptional
        << "\n";

    std::cout
        << "NORMAL_SUM="
        << central_gap_sum_normal
        << "\n";

    std::cout
        << "EXCEPTIONAL_AVG="
        << (
            exceptional_A == 0
                ? 0.0
                : static_cast<double>(
                    central_gap_sum_exceptional
                ) /
                  static_cast<double>(
                    exceptional_A
                )
        )
        << "\n";

    std::cout
        << "NORMAL_AVG="
        << (
            normal_A == 0
                ? 0.0
                : static_cast<double>(
                    central_gap_sum_normal
                ) /
                  static_cast<double>(
                    normal_A
                )
        )
        << "\n";

    std::cout
        << "EXCEPTIONAL_MAX="
        << largest_exceptional_gap
        << "\n";

    std::cout
        << "NORMAL_MAX="
        << largest_normal_gap
        << "\n";

    std::cout
        << "\nWINDOW_COUNTS\n";

    std::cout
        << "EXCEPTIONAL_W90_SUM="
        << window90_exceptional
        << "\n";

    std::cout
        << "NORMAL_W90_SUM="
        << window90_normal
        << "\n";

    std::cout
        << "EXCEPTIONAL_W95_SUM="
        << window95_exceptional
        << "\n";

    std::cout
        << "NORMAL_W95_SUM="
        << window95_normal
        << "\n";

    std::cout
        << "EXCEPTIONAL_W99_SUM="
        << window99_exceptional
        << "\n";

    std::cout
        << "NORMAL_W99_SUM="
        << window99_normal
        << "\n";

    std::cout
        << "\nCLOSE_PAIR_FRACTIONS\n";

    std::cout
        << "EXCEPTIONAL_REL_GAP_LE_10PCT="
        << exceptional_close_pair_10
        << "\n";

    std::cout
        << "EXCEPTIONAL_REL_GAP_LE_5PCT="
        << exceptional_close_pair_5
        << "\n";

    std::cout
        << "EXCEPTIONAL_REL_GAP_LE_1PCT="
        << exceptional_close_pair_1
        << "\n";

    std::cout
        << "NORMAL_REL_GAP_LE_10PCT="
        << normal_close_pair_10
        << "\n";

    std::cout
        << "NORMAL_REL_GAP_LE_5PCT="
        << normal_close_pair_5
        << "\n";

    std::cout
        << "NORMAL_REL_GAP_LE_1PCT="
        << normal_close_pair_1
        << "\n";

    std::cout
        << "\nFIRST_EXCEPTION_POSITION\n";

    std::cout
        << "BELOW_SQRT="
        << exceptional_first_k_below_sqrt
        << "\n";

    std::cout
        << "ABOVE_SQRT="
        << exceptional_first_k_above_sqrt
        << "\n";

    std::cout
        << "DISTANCE_SUM="
        << exceptional_first_k_distance_sum
        << "\n";

    std::cout
        << "DISTANCE_AVG="
        << (
            exceptional_A == 0
                ? 0.0
                : static_cast<double>(
                    exceptional_first_k_distance_sum
                ) /
                  static_cast<double>(
                    exceptional_A
                )
        )
        << "\n";

    std::cout
        << "\nEXCEPTIONAL_GEOMETRY\n";

    for (
        const std::string& line :
        exceptional_output
    ) {
        std::cout
            << line
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT "
        << EXPERIMENT
        << "\n";

    return 0;
}
