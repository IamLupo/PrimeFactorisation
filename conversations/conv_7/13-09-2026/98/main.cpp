#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Case {
    u64 p;
    u64 a0;
    u64 b;
    u64 z;
    u64 q;
};

struct BuiltCase {
    u128 e;
    u128 s0;
    u128 p_pow_e;
    u128 m;
};

struct MixedRadix {
    std::vector<u64> m_digits;
    std::vector<u128> p_powers;
    std::vector<u128> weights;
    u128 count;
};

void print_u128(u128 value) {
    if (value == 0) {
        std::cout << '0';
        return;
    }

    std::string s;

    while (value > 0) {
        unsigned digit =
            static_cast<unsigned>(value % 10);

        s.push_back(
            static_cast<char>('0' + digit)
        );

        value /= 10;
    }

    std::reverse(s.begin(), s.end());
    std::cout << s;
}

u128 pow_u128(u64 base, u64 exp) {
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

BuiltCase build_case(const Case &c) {
    BuiltCase bc;

    bc.e =
        static_cast<u128>(c.a0) +
        static_cast<u128>(c.z) +
        1;

    bc.s0 =
        static_cast<u128>(c.b + 1) *
        pow_u128(c.p, c.a0);

    bc.p_pow_e =
        pow_u128(
            c.p,
            static_cast<u64>(bc.e)
        );

    bc.m =
        bc.s0 +
        static_cast<u128>(c.q) *
        bc.p_pow_e -
        1;

    return bc;
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

MixedRadix build_mixed_radix(
    const Case &c,
    const BuiltCase &bc
) {
    MixedRadix mr;

    mr.m_digits =
        digits_base(
            bc.m,
            c.p
        );

    std::size_t L =
        mr.m_digits.size();

    mr.p_powers.resize(L);
    mr.weights.resize(L);

    u128 weight = 1;
    u128 count = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        mr.p_powers[i] =
            pow_u128(
                c.p,
                static_cast<u64>(i)
            );

        mr.weights[i] =
            weight;

        weight *=
            static_cast<u128>(
                mr.m_digits[i] + 1
            );

        count *=
            static_cast<u128>(
                mr.m_digits[i] + 1
            );
    }

    mr.count = count;

    return mr;
}

bool is_miss(
    u128 n,
    u128 m,
    u64 p
) {
    std::vector<u64> nd =
        digits_base(n, p);

    std::vector<u64> md =
        digits_base(m, p);

    std::size_t L =
        std::max(
            nd.size(),
            md.size()
        );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 a =
            i < nd.size()
                ? nd[i]
                : 0;

        u64 b =
            i < md.size()
                ? md[i]
                : 0;

        if (a > b) {
            return false;
        }
    }

    return true;
}

/*
 * Rank of an admissible/MISS number:
 *
 * R(n)=sum_i n_i W_i
 *
 * This is mixed-radix rank.
 */
u128 rank_miss(
    u128 n,
    u64 p,
    const MixedRadix &mr
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < mr.m_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(
                n % p
            );

        n /= p;

        rank +=
            static_cast<u128>(digit) *
            mr.weights[i];
    }

    return rank;
}

/*
 * Mixed-radix unranking.
 *
 * Each output digit satisfies
 *
 *   0 <= digit <= m_i.
 *
 * Therefore every result is a MISS number.
 */
u128 unrank_miss(
    u128 rank,
    u64 p,
    const MixedRadix &mr
) {
    u128 value = 0;
    u128 place = 1;

    for (std::size_t i = 0;
         i < mr.m_digits.size();
         ++i) {

        u128 radix =
            static_cast<u128>(
                mr.m_digits[i] + 1
            );

        u128 digit =
            rank % radix;

        rank /= radix;

        value +=
            digit * place;

        place *=
            static_cast<u128>(p);
    }

    return value;
}

/*
 * Direct prefix count M(n):
 *
 * number of MISS values x <= n.
 *
 * This is the digit formula from Experiment 250.
 */
u128 miss_prefix_count(
    u128 n,
    u128 m,
    u64 p
) {
    if (n > m) {
        n = m;
    }

    std::vector<u64> nd =
        digits_base(n, p);

    std::vector<u64> md =
        digits_base(m, p);

    std::size_t L =
        std::max(
            nd.size(),
            md.size()
        );

    std::vector<u64> n_digits(
        L,
        0
    );

    std::vector<u64> m_digits(
        L,
        0
    );

    for (std::size_t i = 0;
         i < nd.size();
         ++i) {

        n_digits[i] =
            nd[i];
    }

    for (std::size_t i = 0;
         i < md.size();
         ++i) {

        m_digits[i] =
            md[i];
    }

    std::vector<u128> weights(
        L + 1,
        1
    );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        weights[i + 1] =
            weights[i] *
            static_cast<u128>(
                m_digits[i] + 1
            );
    }

    std::size_t h = L;

    for (std::size_t i = L;
         i-- > 0;) {

        if (n_digits[i] >
            m_digits[i]) {

            h = i;
            break;
        }
    }

    u128 misses = 0;

    if (h == L) {
        misses = 1;

        for (std::size_t i = 0;
             i < L;
             ++i) {

            misses +=
                static_cast<u128>(
                    n_digits[i]
                ) *
                weights[i];
        }

        return misses;
    }

    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        misses +=
            static_cast<u128>(
                n_digits[i]
            ) *
            weights[i];
    }

    misses +=
        static_cast<u128>(
            m_digits[h] + 1
        ) *
        weights[h];

    return misses;
}

/*
 * For any n:
 *
 * M(n)-1 is the mixed-radix rank of
 * the largest MISS number <= n.
 *
 * This function constructs that predecessor
 * explicitly for verification.
 */
u128 largest_miss_leq(
    u128 n,
    u128 m,
    u64 p
) {
    if (n > m) {
        n = m;
    }

    std::vector<u64> nd =
        digits_base(n, p);

    std::vector<u64> md =
        digits_base(m, p);

    std::size_t L =
        std::max(
            nd.size(),
            md.size()
        );

    nd.resize(L, 0);
    md.resize(L, 0);

    std::size_t h = L;

    for (std::size_t i = L;
         i-- > 0;) {

        if (nd[i] >
            md[i]) {

            h = i;
            break;
        }
    }

    /*
     * n itself is a MISS.
     */
    if (h == L) {
        u128 value = 0;
        u128 place = 1;

        for (std::size_t i = 0;
             i < L;
             ++i) {

            value +=
                static_cast<u128>(
                    nd[i]
                ) *
                place;

            place *=
                static_cast<u128>(p);
        }

        return value;
    }

    /*
     * Most significant violating digit:
     *
     * lower digits are set to p-1 conceptually,
     * but the actual largest MISS must instead
     * keep all higher prefix digits equal and
     * choose m_h at the violating position.
     *
     * Lower digits then become m_i to maximize
     * the admissible value.
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
            static_cast<u128>(p);
    }

    return value;
}

void print_case(const Case &c) {
    BuiltCase bc =
        build_case(c);

    MixedRadix mr =
        build_mixed_radix(
            c,
            bc
        );

    std::cout
        << "p=" << c.p
        << " a0=" << c.a0
        << " b=" << c.b
        << " z=" << c.z
        << " q=" << c.q
        << " e=";

    print_u128(bc.e);

    std::cout
        << " s0=";

    print_u128(bc.s0);

    std::cout
        << " m=";

    print_u128(bc.m);

    std::cout
        << " miss_count=";

    print_u128(mr.count);

    std::cout
        << '\n';
}

/*
 * --------------------------------------------------------------------------
 * Deterministic exhaustive rank/unrank tests
 * --------------------------------------------------------------------------
 */
bool deterministic_tests() {
    const std::vector<Case> cases = {
        {2, 0, 0, 0, 1},
        {2, 0, 0, 0, 3},
        {2, 1, 0, 0, 7},
        {3, 0, 0, 0, 2},
        {5, 0, 0, 0, 4},
        {5, 3, 1, 4, 100}
    };

    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        BuiltCase bc =
            build_case(c);

        MixedRadix mr =
            build_mixed_radix(
                c,
                bc
            );

        bool pass = true;

        /*
         * Exhaust all MISS ranks when feasible.
         */
        if (mr.count <=
            static_cast<u128>(1000000)) {

            u64 count =
                static_cast<u64>(
                    mr.count
                );

            for (u64 k = 0;
                 k < count;
                 ++k) {

                u128 value =
                    unrank_miss(
                        static_cast<u128>(k),
                        c.p,
                        mr
                    );

                if (!is_miss(
                        value,
                        bc.m,
                        c.p
                    )) {

                    pass = false;
                    break;
                }

                u128 rank =
                    rank_miss(
                        value,
                        c.p,
                        mr
                    );

                if (rank !=
                    static_cast<u128>(k)) {

                    pass = false;
                    break;
                }

                u128 prefix =
                    miss_prefix_count(
                        value,
                        bc.m,
                        c.p
                    );

                if (prefix !=
                    static_cast<u128>(k) + 1) {

                    pass = false;
                    break;
                }
            }
        } else {
            /*
             * Large cases:
             * first, middle, last ranks.
             */
            std::vector<u128> ranks;

            ranks.push_back(0);
            ranks.push_back(
                mr.count / 2
            );
            ranks.push_back(
                mr.count - 1
            );

            for (u128 k : ranks) {
                u128 value =
                    unrank_miss(
                        k,
                        c.p,
                        mr
                    );

                if (!is_miss(
                        value,
                        bc.m,
                        c.p
                    )) {

                    pass = false;
                    break;
                }

                if (rank_miss(
                        value,
                        c.p,
                        mr
                    ) != k) {

                    pass = false;
                    break;
                }

                if (miss_prefix_count(
                        value,
                        bc.m,
                        c.p
                    ) != k + 1) {

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
        << cases.size()
        << " deterministic_failures="
        << failures
        << " deterministic_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * --------------------------------------------------------------------------
 * Random rank/unrank tests
 * --------------------------------------------------------------------------
 */
bool random_rank_tests() {
    std::mt19937_64 rng(
        0x251123456789ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 choices[] = {
            2, 3, 5, 7, 11
        };

        u64 p =
            choices[rng() % 5];

        u64 a0 =
            rng() % 8;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 8;

        u64 q =
            1 +
            rng() %
            1000000000000000ULL;

        Case c{
            p,
            a0,
            b,
            z,
            q
        };

        BuiltCase bc =
            build_case(c);

        MixedRadix mr =
            build_mixed_radix(
                c,
                bc
            );

        if (mr.count == 0) {
            continue;
        }

        u128 k;

        if (mr.count <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            k =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(
                        mr.count
                    )
                );

        } else {
            u128 r =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            k =
                r % mr.count;
        }

        u128 value =
            unrank_miss(
                k,
                c.p,
                mr
            );

        bool ok =
            is_miss(
                value,
                bc.m,
                c.p
            );

        if (ok) {
            ok =
                rank_miss(
                    value,
                    c.p,
                    mr
                ) == k;
        }

        if (ok) {
            ok =
                miss_prefix_count(
                    value,
                    bc.m,
                    c.p
                ) == k + 1;
        }

        if (!ok) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANK FAILURE\n";

                print_case(c);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout
                    << " value=";

                print_u128(value);

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "rank_cases="
        << cases
        << " rank_failures="
        << failures
        << " rank_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * --------------------------------------------------------------------------
 * Test arbitrary n:
 *
 * M(n)-1 must equal the rank of the largest MISS <= n.
 * --------------------------------------------------------------------------
 */
bool predecessor_rank_tests() {
    std::mt19937_64 rng(
        0x251777888ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 choices[] = {
            2, 3, 5, 7, 11
        };

        u64 p =
            choices[rng() % 5];

        u64 a0 =
            rng() % 8;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 8;

        u64 q =
            1 +
            rng() %
            1000000000000000ULL;

        Case c{
            p,
            a0,
            b,
            z,
            q
        };

        BuiltCase bc =
            build_case(c);

        MixedRadix mr =
            build_mixed_radix(
                c,
                bc
            );

        if (bc.m >
            static_cast<u128>(
                UINT64_MAX
            )) {

            continue;
        }

        u64 m =
            static_cast<u64>(
                bc.m
            );

        u64 n =
            rng() %
            (m + 1);

        u128 prefix =
            miss_prefix_count(
                static_cast<u128>(n),
                bc.m,
                c.p
            );

        if (prefix == 0) {
            ++failures;
            continue;
        }

        u128 predecessor =
            largest_miss_leq(
                static_cast<u128>(n),
                bc.m,
                c.p
            );

        if (!is_miss(
                predecessor,
                bc.m,
                c.p
            )) {

            ++failures;
            continue;
        }

        u128 predecessor_rank =
            rank_miss(
                predecessor,
                c.p,
                mr
            );

        if (predecessor_rank + 1 != prefix) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "PREDECESSOR FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
                    << " predecessor=";

                print_u128(predecessor);

                std::cout
                    << " prefix=";

                print_u128(prefix);

                std::cout
                    << " rank=";

                print_u128(predecessor_rank);

                std::cout
                    << '\n';
            }
        }

        /*
         * If n itself is a MISS, predecessor must equal n.
         */
        if (is_miss(
                static_cast<u128>(n),
                bc.m,
                c.p
            )) {

            if (predecessor !=
                static_cast<u128>(n)) {

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
 * Total bridge:
 *
 * number of MISS values = product_i(m_i+1)
 *
 * and therefore their ranks are exactly
 *
 * 0,...,product_i(m_i+1)-1.
 * --------------------------------------------------------------------------
 */
bool total_bijection_tests() {
    std::mt19937_64 rng(
        0x251424242ULL
    );

    const u64 cases = 50000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 choices[] = {
            2, 3, 5, 7
        };

        u64 p =
            choices[rng() % 4];

        u128 m =
            static_cast<u128>(
                rng()
            );

        MixedRadix mr;

        mr.m_digits =
            digits_base(
                m,
                p
            );

        std::size_t L =
            mr.m_digits.size();

        mr.weights.resize(L);

        u128 weight = 1;
        u128 count = 1;

        for (std::size_t j = 0;
             j < L;
             ++j) {

            mr.weights[j] =
                weight;

            weight *=
                static_cast<u128>(
                    mr.m_digits[j] + 1
                );

            count *=
                static_cast<u128>(
                    mr.m_digits[j] + 1
                );
        }

        mr.count = count;

        if (count == 0) {
            ++failures;
            continue;
        }

        /*
         * Check the endpoints.
         */
        u128 first =
            unrank_miss(
                0,
                p,
                mr
            );

        u128 last =
            unrank_miss(
                count - 1,
                p,
                mr
            );

        if (first != 0) {
            ++failures;
            continue;
        }

        if (last != m) {
            ++failures;
            continue;
        }

        if (rank_miss(
                first,
                p,
                mr
            ) != 0) {

            ++failures;
            continue;
        }

        if (rank_miss(
                last,
                p,
                mr
            ) != count - 1) {

            ++failures;
            continue;
        }
    }

    std::cout
        << "total_cases="
        << cases
        << " total_failures="
        << failures
        << " total_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 251\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_rank_ok =
        random_rank_tests();

    bool predecessor_ok =
        predecessor_rank_tests();

    bool total_ok =
        total_bijection_tests();

    bool overall =
        deterministic_ok &&
        random_rank_ok &&
        predecessor_ok &&
        total_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 251\n";

    return overall ? 0 : 1;
}
