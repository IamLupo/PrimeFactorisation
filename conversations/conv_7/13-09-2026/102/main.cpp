#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Data {
    u64 p;
    u128 m;
    std::vector<u64> m_digits;
    std::vector<u128> p_powers;
    std::vector<u128> weights;
    u128 miss_count;
};

struct CarryInfo {
    std::size_t r;
    bool overflow;
};

void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    std::string s;

    while (x > 0) {
        unsigned digit =
            static_cast<unsigned>(x % 10);

        s.push_back(
            static_cast<char>('0' + digit)
        );

        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    std::cout << s;
}

u128 pow_u128(
    u64 base,
    u64 exp
) {
    u128 result = 1;
    u128 b = base;

    while (exp > 0) {
        if (exp & 1ULL) {
            result *= b;
        }

        exp >>= 1ULL;

        if (exp != 0) {
            b *= b;
        }
    }

    return result;
}

std::vector<u64> digits_base(
    u128 value,
    u64 p
) {
    std::vector<u64> digits;

    while (value > 0) {
        digits.push_back(
            static_cast<u64>(
                value % p
            )
        );

        value /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

Data build_data(
    u64 p,
    u128 m
) {
    Data data;

    data.p = p;
    data.m = m;

    data.m_digits =
        digits_base(
            m,
            p
        );

    const std::size_t L =
        data.m_digits.size();

    data.p_powers.resize(L);
    data.weights.resize(L);

    u128 weight = 1;
    u128 count = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        data.p_powers[i] =
            pow_u128(
                p,
                static_cast<u64>(i)
            );

        data.weights[i] =
            weight;

        weight *=
            static_cast<u128>(
                data.m_digits[i] + 1
            );

        count *=
            static_cast<u128>(
                data.m_digits[i] + 1
            );
    }

    data.miss_count = count;

    return data;
}

/*
 * Mixed-radix unrank of the MISS set.
 *
 * Every digit satisfies
 *
 *   0 <= x_i <= m_i.
 */
u128 unrank_miss(
    u128 k,
    const Data &data
) {
    u128 value = 0;
    u128 place = 1;

    for (std::size_t i = 0;
         i < data.m_digits.size();
         ++i) {

        u128 radix =
            static_cast<u128>(
                data.m_digits[i] + 1
            );

        u128 digit =
            k % radix;

        k /= radix;

        value +=
            digit * place;

        place *=
            static_cast<u128>(
                data.p
            );
    }

    return value;
}

u128 rank_miss(
    u128 value,
    const Data &data
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < data.m_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(
                value % data.p
            );

        value /= data.p;

        rank +=
            static_cast<u128>(
                digit
            ) *
            data.weights[i];
    }

    return rank;
}

bool is_miss(
    u128 value,
    const Data &data
) {
    std::vector<u64> digits =
        digits_base(
            value,
            data.p
        );

    const std::size_t L =
        std::max(
            digits.size(),
            data.m_digits.size()
        );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 a =
            i < digits.size()
                ? digits[i]
                : 0;

        u64 b =
            i < data.m_digits.size()
                ? data.m_digits[i]
                : 0;

        if (a > b) {
            return false;
        }
    }

    return true;
}

/*
 * Carry position when k -> k+1.
 *
 * This is entirely mixed-radix and makes no
 * reference to the arithmetic-sequence structure.
 */
CarryInfo carry_position(
    u128 k,
    const Data &data
) {
    for (std::size_t i = 0;
         i < data.m_digits.size();
         ++i) {

        u128 digit =
            k %
            static_cast<u128>(
                data.m_digits[i] + 1
            );

        k /=
            static_cast<u128>(
                data.m_digits[i] + 1
            );

        if (digit <
            static_cast<u128>(
                data.m_digits[i]
            )) {

            return {
                i,
                false
            };
        }
    }

    return {
        data.m_digits.size(),
        true
    };
}

/*
 * General carry-gap formula.
 *
 * When digit r is incremented:
 *
 *   lower digits:
 *       m_0,...,m_{r-1}
 *
 * become
 *
 *       0,...,0.
 *
 * Therefore
 *
 *   M_{k+1}-M_k
 *     = p^r - sum_{i<r} m_i p^i.
 *
 * The open gap has one fewer element.
 */
u128 predicted_gap(
    const Data &data,
    std::size_t r
) {
    u128 lower = 0;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        lower +=
            static_cast<u128>(
                data.m_digits[i]
            ) *
            data.p_powers[i];
    }

    return
        data.p_powers[r] -
        lower -
        1;
}

u128 actual_gap(
    u128 left,
    u128 right
) {
    if (right <= left + 1) {
        return 0;
    }

    return right - left - 1;
}

/*
 * Independently calculate the lower-digit
 * expression from the actual m digits.
 */
u128 direct_difference_formula(
    const Data &data,
    std::size_t r
) {
    u128 lower = 0;

    u128 p_power = 1;

    for (std::size_t i = 0;
         i < r;
         ++i) {

        lower +=
            static_cast<u128>(
                data.m_digits[i]
            ) *
            p_power;

        p_power *=
            static_cast<u128>(
                data.p
            );
    }

    return
        p_power -
        lower -
        1;
}

/*
 * For a nonempty gap, every interior point
 * must be a HIT, because the two endpoints are
 * consecutive MISS values.
 *
 * We do not scan the interval.
 *
 * Instead verify a few boundary points:
 *
 *   left      = MISS
 *   left + 1  = HIT
 *   right - 1 = HIT
 *   right     = MISS
 */
bool verify_gap_boundaries(
    const Data &data,
    u128 left,
    u128 right
) {
    if (!is_miss(left, data)) {
        return false;
    }

    if (!is_miss(right, data)) {
        return false;
    }

    if (right <= left + 1) {
        return true;
    }

    if (is_miss(
            left + 1,
            data
        )) {
        return false;
    }

    if (is_miss(
            right - 1,
            data
        )) {
        return false;
    }

    return true;
}

void print_case(
    const Data &data
) {
    std::cout
        << "p="
        << data.p
        << " m=";

    print_u128(data.m);

    std::cout
        << " miss_count=";

    print_u128(
        data.miss_count
    );

    std::cout
        << '\n';
}

/*
 * --------------------------------------------------------------------------
 * Deterministic arbitrary-m tests
 * --------------------------------------------------------------------------
 */
bool deterministic_tests() {
    struct Test {
        u64 p;
        u128 m;
    };

    const std::vector<Test> tests = {
        {2, 0},
        {2, 1},
        {2, 2},
        {2, 7},
        {2, 29},
        {3, 6},
        {3, 17},
        {5, 20},
        {5, 100},
        {7, 1234},
        {11, 987654321ULL},
        {13, 1000000000000ULL}
    };

    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CASES\n";

    for (const Test &test : tests) {
        Data data =
            build_data(
                test.p,
                test.m
            );

        print_case(data);

        bool pass = true;

        if (data.miss_count <=
            static_cast<u128>(100000)) {

            for (u64 k = 0;
                 k + 1 <
                 static_cast<u64>(
                     data.miss_count
                 );
                 ++k) {

                CarryInfo carry =
                    carry_position(
                        static_cast<u128>(k),
                        data
                    );

                if (carry.overflow) {
                    pass = false;
                    break;
                }

                u128 left =
                    unrank_miss(
                        static_cast<u128>(k),
                        data
                    );

                u128 right =
                    unrank_miss(
                        static_cast<u128>(k + 1),
                        data
                    );

                u128 actual =
                    actual_gap(
                        left,
                        right
                    );

                u128 predicted =
                    predicted_gap(
                        data,
                        carry.r
                    );

                u128 direct =
                    direct_difference_formula(
                        data,
                        carry.r
                    );

                if (actual != predicted ||
                    predicted != direct) {

                    pass = false;
                    break;
                }

                if (!verify_gap_boundaries(
                        data,
                        left,
                        right
                    )) {

                    pass = false;
                    break;
                }
            }

        } else {
            std::vector<u128> ranks;

            ranks.push_back(0);
            ranks.push_back(1);
            ranks.push_back(
                data.miss_count / 4
            );
            ranks.push_back(
                data.miss_count / 2
            );
            ranks.push_back(
                data.miss_count - 2
            );

            for (u128 k : ranks) {
                if (k + 1 >= data.miss_count) {
                    continue;
                }

                CarryInfo carry =
                    carry_position(
                        k,
                        data
                    );

                u128 left =
                    unrank_miss(
                        k,
                        data
                    );

                u128 right =
                    unrank_miss(
                        k + 1,
                        data
                    );

                u128 actual =
                    actual_gap(
                        left,
                        right
                    );

                u128 predicted =
                    predicted_gap(
                        data,
                        carry.r
                    );

                if (actual != predicted) {
                    pass = false;
                    break;
                }

                if (!verify_gap_boundaries(
                        data,
                        left,
                        right
                    )) {

                    pass = false;
                    break;
                }
            }
        }

        std::cout
            << "pass="
            << (pass ? 1 : 0)
            << '\n';

        if (!pass) {
            ++failures;
        }
    }

    std::cout
        << "deterministic_cases="
        << tests.size()
        << " deterministic_failures="
        << failures
        << " deterministic_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * --------------------------------------------------------------------------
 * Random arbitrary-m tests
 * --------------------------------------------------------------------------
 */
bool random_tests() {
    std::mt19937_64 rng(
        0x254123456789ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 bases[] = {
            2, 3, 5, 7, 11, 13
        };

        u64 p =
            bases[
                rng() % 6
            ];

        /*
         * Generate genuinely arbitrary m.
         *
         * Mix small, medium and huge values.
         */
        u128 m;

        switch (i % 5) {
            case 0:
                m =
                    static_cast<u128>(
                        rng() % 1000
                    );
                break;

            case 1:
                m =
                    static_cast<u128>(
                        rng() % 1000000000ULL
                    );
                break;

            case 2:
                m =
                    static_cast<u128>(rng());
                break;

            case 3:
                m =
                    (static_cast<u128>(rng()) << 32) |
                    static_cast<u128>(
                        rng()
                    );
                break;

            default:
                m =
                    (static_cast<u128>(rng()) << 64) |
                    static_cast<u128>(rng());
                break;
        }

        Data data =
            build_data(
                p,
                m
            );

        if (data.miss_count <= 1) {
            continue;
        }

        u128 k;

        if (data.miss_count <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            k =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(
                        data.miss_count - 1
                    )
                );

        } else {
            u128 random128 =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            k =
                random128 %
                (data.miss_count - 1);
        }

        CarryInfo carry =
            carry_position(
                k,
                data
            );

        if (carry.overflow) {
            ++failures;
            continue;
        }

        u128 left =
            unrank_miss(
                k,
                data
            );

        u128 right =
            unrank_miss(
                k + 1,
                data
            );

        u128 actual =
            actual_gap(
                left,
                right
            );

        u128 predicted =
            predicted_gap(
                data,
                carry.r
            );

        u128 direct =
            direct_difference_formula(
                data,
                carry.r
            );

        if (actual != predicted ||
            predicted != direct) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(data);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout
                    << " r="
                    << carry.r
                    << " actual=";

                print_u128(actual);

                std::cout
                    << " predicted=";

                print_u128(predicted);

                std::cout
                    << '\n';
            }

            continue;
        }

        if (!verify_gap_boundaries(
                data,
                left,
                right
            )) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "BOUNDARY FAILURE\n";

                print_case(data);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "random_cases="
        << cases
        << " random_failures="
        << failures
        << " random_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * --------------------------------------------------------------------------
 * Successor identity
 * --------------------------------------------------------------------------
 *
 * Directly verify:
 *
 *   M_{k+1}-M_k
 *     = p^r - lower_m
 *
 * rather than subtracting 1 afterward.
 */
bool successor_tests() {
    std::mt19937_64 rng(
        0x254777888ULL
    );

    const u64 cases = 50000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 bases[] = {
            2, 3, 5, 7, 11
        };

        u64 p =
            bases[
                rng() % 5
            ];

        u128 m =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        Data data =
            build_data(
                p,
                m
            );

        if (data.miss_count <= 1) {
            continue;
        }

        u128 k;

        if (data.miss_count <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            k =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(
                        data.miss_count - 1
                    )
                );

        } else {
            u128 random128 =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            k =
                random128 %
                (data.miss_count - 1);
        }

        CarryInfo carry =
            carry_position(
                k,
                data
            );

        u128 left =
            unrank_miss(
                k,
                data
            );

        u128 right =
            unrank_miss(
                k + 1,
                data
            );

        u128 lower = 0;

        for (std::size_t j = 0;
             j < carry.r;
             ++j) {

            lower +=
                static_cast<u128>(
                    data.m_digits[j]
                ) *
                data.p_powers[j];
        }

        u128 expected_difference =
            data.p_powers[
                carry.r
            ] - lower;

        u128 actual_difference =
            right - left;

        if (actual_difference !=
            expected_difference) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "SUCCESSOR FAILURE\n";

                print_case(data);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout
                    << " r=";

                std::cout
                    << carry.r;

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "successor_cases="
        << cases
        << " successor_failures="
        << failures
        << " successor_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 254\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool successor_ok =
        successor_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        successor_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 254\n";

    return overall ? 0 : 1;
}
