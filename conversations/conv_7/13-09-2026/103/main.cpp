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

void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    std::string s;

    while (x > 0) {
        unsigned d =
            static_cast<unsigned>(x % 10);

        s.push_back(
            static_cast<char>('0' + d)
        );

        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    std::cout << s;
}

std::vector<u64> digits_base(
    u128 value,
    u64 p
) {
    std::vector<u64> digits;

    while (value > 0) {
        digits.push_back(
            static_cast<u64>(value % p)
        );

        value /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
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

Data build_data(
    u64 p,
    u128 m
) {
    Data data;

    data.p = p;
    data.m = m;

    data.m_digits =
        digits_base(m, p);

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
 * Digitwise MISS predicate:
 *
 *     n <=_p m
 */
bool is_miss(
    u128 n,
    const Data &data
) {
    std::vector<u64> nd =
        digits_base(n, data.p);

    const std::size_t L =
        std::max(
            nd.size(),
            data.m_digits.size()
        );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 a =
            i < nd.size()
                ? nd[i]
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
 * Mixed-radix rank of a MISS number.
 *
 *     R(x) = sum_i x_i W_i
 */
u128 rank_miss(
    u128 x,
    const Data &data
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < data.m_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(
                x % data.p
            );

        x /= data.p;

        rank +=
            static_cast<u128>(digit) *
            data.weights[i];
    }

    return rank;
}

/*
 * Mixed-radix unranking.
 */
u128 unrank_miss(
    u128 rank,
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
            rank % radix;

        rank /= radix;

        value +=
            digit * place;

        place *=
            static_cast<u128>(
                data.p
            );
    }

    return value;
}

/*
 * Largest MISS <= n.
 *
 * If n itself is admissible, return n.
 *
 * Otherwise let h be the most significant digit
 * with n_h > m_h.
 *
 * Then:
 *
 *   x_i = n_i  for i > h
 *   x_i = m_i  for i <= h
 *
 * gives the largest admissible predecessor.
 */
u128 miss_predecessor(
    u128 n,
    const Data &data
) {
    if (n > data.m) {
        n = data.m;
    }

    std::vector<u64> nd =
        digits_base(n, data.p);

    const std::size_t L =
        std::max(
            nd.size(),
            data.m_digits.size()
        );

    nd.resize(L, 0);

    std::vector<u64> md =
        data.m_digits;

    md.resize(L, 0);

    /*
     * Find most significant violation.
     */
    std::size_t h = L;

    for (std::size_t i = L;
         i-- > 0;) {

        if (nd[i] > md[i]) {
            h = i;
            break;
        }
    }

    /*
     * n itself is a MISS.
     */
    if (h == L) {
        return n;
    }

    /*
     * Clamp at h and below.
     */
    u128 value = 0;
    u128 place = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 digit;

        if (i > h) {
            digit = nd[i];
        } else {
            digit = md[i];
        }

        value +=
            static_cast<u128>(digit) *
            place;

        place *=
            static_cast<u128>(
                data.p
            );
    }

    return value;
}

/*
 * Direct digit-prefix MISS count:
 *
 * M(n) = #{x <= n : x <=_p m}
 */
u128 miss_prefix_count(
    u128 n,
    const Data &data
) {
    if (n > data.m) {
        n = data.m;
    }

    std::vector<u64> nd =
        digits_base(n, data.p);

    const std::size_t L =
        std::max(
            nd.size(),
            data.m_digits.size()
        );

    nd.resize(L, 0);

    std::vector<u64> md =
        data.m_digits;

    md.resize(L, 0);

    /*
     * Lower-digit products:
     *
     * W_i = prod_{j<i}(m_j+1)
     */
    std::vector<u128> W(
        L + 1,
        1
    );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        W[i + 1] =
            W[i] *
            static_cast<u128>(
                md[i] + 1
            );
    }

    std::size_t h = L;

    for (std::size_t i = L;
         i-- > 0;) {

        if (nd[i] > md[i]) {
            h = i;
            break;
        }
    }

    u128 misses = 0;

    if (h == L) {
        /*
         * n itself is admissible.
         */
        misses = 1;

        for (std::size_t i = 0;
             i < L;
             ++i) {

            misses +=
                static_cast<u128>(
                    nd[i]
                ) *
                W[i];
        }

        return misses;
    }

    /*
     * Higher prefix is fixed.
     */
    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        misses +=
            static_cast<u128>(
                nd[i]
            ) *
            W[i];
    }

    /*
     * At the first violating digit,
     * x_h may take 0,...,m_h.
     */
    misses +=
        static_cast<u128>(
            md[h] + 1
        ) *
        W[h];

    return misses;
}

/*
 * Closed HIT prefix count from the digit formula.
 */
u128 hit_count_digits(
    u128 n,
    const Data &data
) {
    if (n > data.m) {
        n = data.m;
    }

    u128 misses =
        miss_prefix_count(
            n,
            data
        );

    return n + 1 - misses;
}

/*
 * New formula:
 *
 *     C(n) = n - R(pred_M(n))
 */
u128 hit_count_predecessor_rank(
    u128 n,
    const Data &data
) {
    if (n > data.m) {
        n = data.m;
    }

    u128 predecessor =
        miss_predecessor(
            n,
            data
        );

    u128 rank =
        rank_miss(
            predecessor,
            data
        );

    return n - rank;
}

/*
 * Completely independent brute-force count
 * for small domains.
 */
u128 brute_hit_count(
    u128 n,
    const Data &data
) {
    u128 count = 0;

    for (u128 x = 0;
         x <= n;
         ++x) {

        if (!is_miss(
                x,
                data
            )) {

            ++count;
        }
    }

    return count;
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
 * Deterministic tests
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

        /*
         * Exhaust tiny domains.
         */
        if (data.m <=
            static_cast<u128>(10000)) {

            for (u64 n = 0;
                 n <=
                 static_cast<u64>(
                     data.m
                 );
                 ++n) {

                u128 a =
                    hit_count_digits(
                        static_cast<u128>(n),
                        data
                    );

                u128 b =
                    hit_count_predecessor_rank(
                        static_cast<u128>(n),
                        data
                    );

                u128 brute =
                    brute_hit_count(
                        static_cast<u128>(n),
                        data
                    );

                if (a != b ||
                    a != brute) {

                    pass = false;
                    break;
                }
            }

        } else {
            /*
             * Large domains:
             * test meaningful points.
             */
            std::vector<u128> points = {
                0,
                1,
                data.m / 4,
                data.m / 2,
                data.m - 1,
                data.m
            };

            for (u128 n : points) {
                u128 a =
                    hit_count_digits(
                        n,
                        data
                    );

                u128 b =
                    hit_count_predecessor_rank(
                        n,
                        data
                    );

                if (a != b) {
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
        0x255123456789ULL
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

        u128 m;

        switch (i % 5) {
            case 0:
                m =
                    static_cast<u128>(
                        rng() % 100000
                    );
                break;

            case 1:
                m =
                    static_cast<u128>(
                        rng() %
                        1000000000000ULL
                    );
                break;

            case 2:
                m =
                    static_cast<u128>(
                        rng()
                    );
                break;

            case 3:
                m =
                    (static_cast<u128>(rng()) << 32) |
                    static_cast<u128>(rng());
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

        u128 n;

        if (m <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            n =
                static_cast<u128>(
                    rng() %
                    (static_cast<u64>(m) + 1)
                );

        } else {
            u128 r =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            n =
                r % (m + 1);
        }

        u128 digit_formula =
            hit_count_digits(
                n,
                data
            );

        u128 predecessor_formula =
            hit_count_predecessor_rank(
                n,
                data
            );

        if (digit_formula !=
            predecessor_formula) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(data);

                std::cout
                    << "n=";

                print_u128(n);

                std::cout
                    << " digit=";

                print_u128(
                    digit_formula
                );

                std::cout
                    << " predecessor=";

                print_u128(
                    predecessor_formula
                );

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
 * Random predecessor structure tests
 * --------------------------------------------------------------------------
 *
 * Verify directly:
 *
 *   pred(n) <= n
 *
 *   pred(n) is a MISS
 *
 *   every x with pred(n) < x <= n is HIT.
 *
 * We only test the endpoints and n,
 * avoiding any gap scan.
 */
bool predecessor_tests() {
    std::mt19937_64 rng(
        0x255777888ULL
    );

    const u64 cases = 100000;
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

        u128 n =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        n %= (m + 1);

        u128 pred =
            miss_predecessor(
                n,
                data
            );

        if (pred > n) {
            ++failures;
            continue;
        }

        if (!is_miss(
                pred,
                data
            )) {

            ++failures;
            continue;
        }

        /*
         * If n itself is MISS, predecessor must be n.
         */
        bool n_is_miss =
            is_miss(
                n,
                data
            );

        if (n_is_miss &&
            pred != n) {

            ++failures;
            continue;
        }

        /*
         * If there is a gap, the immediate point
         * after predecessor must be HIT.
         */
        if (pred < n) {
            if (is_miss(
                    pred + 1,
                    data
                )) {

                ++failures;
                continue;
            }

            /*
             * The predecessor rank gives the exact
             * MISS prefix count:
             *
             * M(n)=rank(pred)+1.
             */
            u128 rank =
                rank_miss(
                    pred,
                    data
                );

            u128 prefix =
                miss_prefix_count(
                    n,
                    data
                );

            if (prefix != rank + 1) {
                ++failures;
                continue;
            }
        }
    }

    std::cout
        << "predecessor_cases="
        << cases
        << " predecessor_failures="
        << failures
        << " predecessor_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * --------------------------------------------------------------------------
 * Endpoint/rank identity
 * --------------------------------------------------------------------------
 *
 * Test:
 *
 *   C(n) = n - R(pred(n))
 *
 * together with:
 *
 *   C(m)=m+1-prod_i(m_i+1)
 *
 * and
 *
 *   R(m)=prod_i(m_i+1)-1.
 */
bool endpoint_tests() {
    std::mt19937_64 rng(
        0x255424242ULL
    );

    const u64 cases = 50000;
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

        u128 m =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        Data data =
            build_data(
                p,
                m
            );

        u128 pred =
            miss_predecessor(
                m,
                data
            );

        if (pred != m) {
            ++failures;
            continue;
        }

        u128 rank =
            rank_miss(
                m,
                data
            );

        if (rank + 1 !=
            data.miss_count) {

            ++failures;
            continue;
        }

        u128 hit_count =
            hit_count_predecessor_rank(
                m,
                data
            );

        u128 expected =
            m + 1 -
            data.miss_count;

        if (hit_count != expected) {
            ++failures;
        }
    }

    std::cout
        << "endpoint_cases="
        << cases
        << " endpoint_failures="
        << failures
        << " endpoint_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 255\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool predecessor_ok =
        predecessor_tests();

    bool endpoint_ok =
        endpoint_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        predecessor_ok &&
        endpoint_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 255\n";

    return overall ? 0 : 1;
}
