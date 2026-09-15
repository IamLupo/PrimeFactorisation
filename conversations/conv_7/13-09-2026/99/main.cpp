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

struct Data {
    std::vector<u64> m_digits;
    std::vector<u128> p_powers;
    std::vector<u128> weights;
};

struct Interval {
    u128 start;
    u128 end;
};

void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    std::string s;

    while (x > 0) {
        s.push_back(
            static_cast<char>(
                '0' + static_cast<unsigned>(x % 10)
            )
        );

        x /= 10;
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

Data build_data(
    const Case &c,
    const BuiltCase &bc
) {
    Data data;

    data.m_digits =
        digits_base(bc.m, c.p);

    const std::size_t L =
        data.m_digits.size();

    data.p_powers.resize(L);
    data.weights.resize(L);

    u128 weight = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        data.p_powers[i] =
            pow_u128(
                c.p,
                static_cast<u64>(i)
            );

        data.weights[i] = weight;

        weight *=
            static_cast<u128>(
                data.m_digits[i] + 1
            );
    }

    return data;
}

u128 miss_count(
    const Data &data
) {
    u128 result = 1;

    for (u64 digit : data.m_digits) {
        result *=
            static_cast<u128>(digit + 1);
    }

    return result;
}

u128 unrank_miss(
    u128 rank,
    u64 p,
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

        value += digit * place;

        place *=
            static_cast<u128>(p);
    }

    return value;
}

u128 rank_miss(
    u128 value,
    u64 p,
    const Data &data
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < data.m_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(value % p);

        value /= p;

        rank +=
            static_cast<u128>(digit) *
            data.weights[i];
    }

    return rank;
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

    const std::size_t L =
        std::max(
            nd.size(),
            md.size()
        );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 a =
            i < nd.size() ? nd[i] : 0;

        u64 b =
            i < md.size() ? md[i] : 0;

        if (a > b) {
            return false;
        }
    }

    return true;
}

bool is_hit(
    u128 n,
    u128 m,
    u64 p
) {
    return !is_miss(n, m, p);
}

/*
 * Recover the HIT interval associated with
 * HIT rank h.
 *
 * This uses the previously established q-based
 * interval construction.
 */
Interval hit_interval_from_rank(
    const Case &c,
    const BuiltCase &bc,
    u128 h
) {
    u128 q =
        (bc.m + 1 - bc.s0) /
        bc.p_pow_e;

    std::vector<u64> q_digits =
        digits_base(q, c.p);

    u128 j = 0;
    u128 place = 1;
    u128 remaining = h;

    for (std::size_t i = 0;
         i < q_digits.size();
         ++i) {

        u128 radix =
            static_cast<u128>(
                q_digits[i] + 1
            );

        u128 digit =
            remaining % radix;

        remaining /= radix;

        j += digit * place;

        place *=
            static_cast<u128>(c.p);
    }

    std::size_t r =
        q_digits.size();

    u128 tmp = j;

    for (std::size_t i = 0;
         i < q_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(
                tmp % c.p
            );

        tmp /= c.p;

        if (digit < q_digits[i]) {
            r = i;
            break;
        }
    }

    u128 p_r =
        pow_u128(
            c.p,
            static_cast<u64>(r)
        );

    u128 successor =
        (j / p_r + 1) * p_r;

    u128 start =
        bc.s0 +
        j * bc.p_pow_e;

    u128 end =
        successor *
        bc.p_pow_e - 1;

    return {start, end};
}

/*
 * Core unified test:
 *
 * MISS_k, MISS_{k+1}
 *
 * gap = [MISS_k+1, MISS_{k+1}-1]
 *
 * must equal HIT interval k/sampled rank.
 */
bool verify_rank_gap(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 k
) {
    u128 total_misses =
        miss_count(data);

    if (k + 1 >= total_misses) {
        return true;
    }

    u128 left =
        unrank_miss(
            k,
            c.p,
            data
        );

    u128 right =
        unrank_miss(
            k + 1,
            c.p,
            data
        );

    if (!is_miss(left, bc.m, c.p)) {
        return false;
    }

    if (!is_miss(right, bc.m, c.p)) {
        return false;
    }

    if (rank_miss(left, c.p, data) != k) {
        return false;
    }

    if (rank_miss(right, c.p, data) != k + 1) {
        return false;
    }

    if (right <= left) {
        return false;
    }

    u128 gap_start = left + 1;
    u128 gap_end = right - 1;

    /*
     * A zero-sized gap is simply consecutive MISS values.
     */
    if (gap_start > gap_end) {
        return true;
    }

    /*
     * The gap is a HIT interval.
     *
     * Its start and end must agree with
     * the independent q/r interval construction.
     */
    u128 q =
        (bc.m + 1 - bc.s0) /
        bc.p_pow_e;

    std::vector<u64> q_digits =
        digits_base(q, c.p);

    u128 hit_count = 1;

    for (u64 digit : q_digits) {
        hit_count *=
            static_cast<u128>(digit + 1);
    }

    hit_count -= 1;

    /*
     * The number of nonempty gaps equals hit_count,
     * but several consecutive MISS ranks may have
     * zero gaps inside MISS blocks.
     *
     * Therefore only test the correspondence when
     * the gap is nonempty.
     *
     * Identify its HIT rank from the number of
     * complete MISS blocks before it.
     */
    u128 preceding_misses = k + 1;

    if (bc.s0 == 0) {
        return false;
    }

    /*
     * A nonempty gap after MISS_k means that k is
     * the last MISS of one MISS block.
     *
     * Since every MISS block has length s0:
     *
     *   k+1 = (h+1)*s0
     *
     * hence
     *
     *   h = (k+1)/s0 - 1.
     */
    if (preceding_misses % bc.s0 != 0) {
        return false;
    }

    u128 h =
        preceding_misses /
        bc.s0 - 1;

    if (h >= hit_count) {
        return false;
    }

    Interval expected =
        hit_interval_from_rank(
            c,
            bc,
            h
        );

    if (gap_start != expected.start) {
        return false;
    }

    if (gap_end != expected.end) {
        return false;
    }

    /*
     * Boundary predicate:
     */
    if (!is_hit(
            gap_start,
            bc.m,
            c.p
        )) {
        return false;
    }

    if (!is_hit(
            gap_end,
            bc.m,
            c.p
        )) {
        return false;
    }

    if (gap_start > 0 &&
        is_hit(
            gap_start - 1,
            bc.m,
            c.p
        )) {
        return false;
    }

    if (gap_end < bc.m &&
        is_hit(
            gap_end + 1,
            bc.m,
            c.p
        )) {
        return false;
    }

    return true;
}

void print_case(const Case &c) {
    BuiltCase bc =
        build_case(c);

    Data data =
        build_data(c, bc);

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

    print_u128(miss_count(data));

    std::cout
        << '\n';
}

bool deterministic_tests() {
    const std::vector<Case> cases = {
        {2, 0, 0, 0, 1},
        {2, 0, 0, 0, 3},
        {2, 1, 0, 0, 7},
        {3, 0, 0, 0, 2},
        {5, 0, 0, 0, 4},
        {5, 3, 1, 4, 100},
        {3, 4, 1, 5, 987654321ULL},
        {5, 3, 1, 4, 1000007654321ULL}
    };

    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        BuiltCase bc =
            build_case(c);

        Data data =
            build_data(c, bc);

        u128 count =
            miss_count(data);

        bool pass = true;

        /*
         * Sample ranks only.
         */
        std::vector<u128> ranks;

        ranks.push_back(0);

        if (count > 2) {
            ranks.push_back(1);
            ranks.push_back(count / 4);
            ranks.push_back(count / 2);
            ranks.push_back(
                count - 2
            );
        }

        if (count > 1) {
            ranks.push_back(
                count - 1
            );
        }

        for (u128 k : ranks) {
            if (!verify_rank_gap(
                    c,
                    bc,
                    data,
                    k
                )) {

                pass = false;
                break;
            }
        }

        /*
         * For small MISS sets, exhaust the gaps.
         */
        if (pass &&
            count <=
                static_cast<u128>(50000)) {

            for (u64 k = 0;
                 k + 1 <
                 static_cast<u64>(count);
                 ++k) {

                if (!verify_rank_gap(
                        c,
                        bc,
                        data,
                        static_cast<u128>(k)
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
        << cases.size()
        << " deterministic_failures="
        << failures
        << " deterministic_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

bool random_tests() {
    std::mt19937_64 rng(
        0x252123456789ULL
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

        Data data =
            build_data(c, bc);

        u128 count =
            miss_count(data);

        if (count <= 1) {
            continue;
        }

        u128 k;

        if (count <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            k =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(
                        count - 1
                    )
                );

        } else {
            u128 random128 =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            k =
                random128 % (count - 1);
        }

        if (!verify_rank_gap(
                c,
                bc,
                data,
                k
            )) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(c);

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
 * Specifically test that every sampled nonempty gap
 * has length exactly the previously derived HIT interval
 * length.
 */
bool length_tests() {
    std::mt19937_64 rng(
        0x252777888ULL
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

        u64 a0 =
            rng() % 7;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 7;

        u64 q =
            1 +
            rng() %
            1000000000000ULL;

        Case c{
            p,
            a0,
            b,
            z,
            q
        };

        BuiltCase bc =
            build_case(c);

        Data data =
            build_data(c, bc);

        u128 count =
            miss_count(data);

        if (count <= 1) {
            continue;
        }

        /*
         * Pick a rank likely to be the end of a MISS block.
         *
         * This directly targets nonempty gaps.
         */
        u128 miss_block_index =
            static_cast<u128>(
                rng() %
                1000000000ULL
            );

        u128 k =
            (miss_block_index + 1) *
            bc.s0;

        if (k == 0 ||
            k >= count) {
            continue;
        }

        --k;

        u128 left =
            unrank_miss(
                k,
                p,
                data
            );

        u128 right =
            unrank_miss(
                k + 1,
                p,
                data
            );

        if (right != left + 1 &&
            !verify_rank_gap(
                c,
                bc,
                data,
                k
            )) {

            ++failures;
        }
    }

    std::cout
        << "length_cases="
        << cases
        << " length_failures="
        << failures
        << " length_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 252R\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool length_ok =
        length_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        length_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 252R\n";

    return overall ? 0 : 1;
}