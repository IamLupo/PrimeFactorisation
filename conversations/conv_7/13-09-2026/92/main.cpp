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
    std::vector<u64> q_digits;
    std::vector<u128> p_powers;
    std::vector<u128> weights;
    std::vector<u128> q_prefix_rank;
};

struct Interval {
    u128 k;
    u128 j;
    u128 r;
    u128 start;
    u128 end;
};

struct Result {
    bool hit;
    u128 rank;
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

    u128 p_to_a =
        pow_u128(c.p, c.a0);

    bc.s0 =
        static_cast<u128>(c.b + 1) *
        p_to_a;

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
    u64 value,
    u64 p
) {
    std::vector<u64> digits;

    while (value > 0) {
        digits.push_back(value % p);
        value /= p;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    return digits;
}

Data build_data(const Case &c) {
    Data data;

    data.q_digits =
        digits_base(c.q, c.p);

    const std::size_t L =
        data.q_digits.size();

    data.p_powers.resize(L);
    data.weights.resize(L);
    data.q_prefix_rank.resize(L);

    u128 weight = 1;
    u128 prefix = 0;

    for (std::size_t i = 0; i < L; ++i) {
        data.p_powers[i] =
            pow_u128(
                c.p,
                static_cast<u64>(i)
            );

        data.weights[i] = weight;

        prefix +=
            static_cast<u128>(
                data.q_digits[i]
            ) * weight;

        data.q_prefix_rank[i] =
            prefix;

        weight *=
            static_cast<u128>(
                data.q_digits[i] + 1
            );
    }

    return data;
}

u128 interval_count(
    const Data &data
) {
    u128 product = 1;

    for (u64 digit : data.q_digits) {
        product *=
            static_cast<u128>(
                digit + 1
            );
    }

    return product - 1;
}

/*
 * k -> j using mixed-radix unranking.
 */
u128 unrank_j(
    u128 k,
    u64 p,
    const Data &data
) {
    u128 j = 0;
    u128 place = 1;

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        u128 radix =
            static_cast<u128>(
                data.q_digits[i] + 1
            );

        u128 digit =
            k % radix;

        k /= radix;

        j +=
            digit * place;

        place *=
            static_cast<u128>(p);
    }

    return j;
}

/*
 * Rank(j) using the mixed-radix weights.
 */
u128 rank_j(
    u128 j,
    u64 p,
    const Data &data
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(j % p);

        j /= p;

        rank +=
            static_cast<u128>(digit) *
            data.weights[i];
    }

    return rank;
}

/*
 * Exact interval generated from k.
 */
bool interval_from_rank(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 k,
    Interval &out
) {
    u128 I =
        interval_count(data);

    if (k >= I) {
        return false;
    }

    u128 j =
        unrank_j(
            k,
            c.p,
            data
        );

    std::size_t r =
        data.q_digits.size();

    u128 temp_j = j;

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(
                temp_j % c.p
            );

        temp_j /= c.p;

        if (digit <
            data.q_digits[i]) {

            r = i;
            break;
        }
    }

    if (r == data.q_digits.size()) {
        return false;
    }

    u128 successor =
        (j / data.p_powers[r] + 1) *
        data.p_powers[r];

    u128 start =
        bc.s0 +
        j * bc.p_pow_e;

    u128 end =
        successor * bc.p_pow_e - 1;

    out.k = k;
    out.j = j;
    out.r = static_cast<u128>(r);
    out.start = start;
    out.end = end;

    return true;
}

/*
 * Direct n -> k localization.
 *
 * This is the compressed formula from Experiment 246,
 * with the final interval-membership test in n-space.
 */
Result localize_n(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 n
) {
    if (n < bc.s0 ||
        n > bc.m) {

        return {false, 0};
    }

    u128 t =
        (n - bc.s0) /
        bc.p_pow_e;

    const std::size_t L =
        data.q_digits.size();

    std::vector<u64> t_digits(L, 0);

    u128 temp = t;

    for (std::size_t i = 0; i < L; ++i) {
        t_digits[i] =
            static_cast<u64>(
                temp % c.p
            );

        temp /= c.p;
    }

    bool violation = false;
    std::size_t h = 0;

    for (std::size_t i = L; i-- > 0;) {
        if (t_digits[i] >
            data.q_digits[i]) {

            violation = true;
            h = i;
            break;
        }
    }

    u128 rank = 0;

    if (!violation) {
        for (std::size_t i = 0; i < L; ++i) {
            rank +=
                static_cast<u128>(
                    t_digits[i]
                ) *
                data.weights[i];
        }
    } else {
        rank =
            data.q_prefix_rank[h];

        for (std::size_t i = h + 1;
             i < L;
             ++i) {

            rank +=
                static_cast<u128>(
                    t_digits[i]
                ) *
                data.weights[i];
        }
    }

    /*
     * Find the first j-digit below q.
     */
    std::size_t r = L;

    if (!violation) {
        for (std::size_t i = 0; i < L; ++i) {
            if (t_digits[i] <
                data.q_digits[i]) {

                r = i;
                break;
            }
        }
    } else {
        for (std::size_t i = h + 1;
             i < L;
             ++i) {

            if (t_digits[i] <
                data.q_digits[i]) {

                r = i;
                break;
            }
        }
    }

    if (r == L) {
        return {false, 0};
    }

    /*
     * Compute floor(j/p^r) without constructing j.
     */
    u128 high = 0;
    u128 place = 1;

    for (std::size_t i = r;
         i < L;
         ++i) {

        u64 digit;

        if (violation && i <= h) {
            digit =
                data.q_digits[i];
        } else {
            digit =
                t_digits[i];
        }

        high +=
            static_cast<u128>(digit) *
            place;

        place *=
            static_cast<u128>(c.p);
    }

    u128 successor =
        (high + 1) *
        data.p_powers[r];

    u128 endpoint_exclusive =
        successor *
        bc.p_pow_e;

    if (n >= endpoint_exclusive) {
        return {false, 0};
    }

    return {true, rank};
}

bool same_result(
    const Result &a,
    const Result &b
) {
    if (a.hit != b.hit) {
        return false;
    }

    if (a.hit &&
        a.rank != b.rank) {

        return false;
    }

    return true;
}

void print_case(const Case &c) {
    BuiltCase bc =
        build_case(c);

    Data data =
        build_data(c);

    std::cout
        << "p=" << c.p
        << " a0=" << c.a0
        << " b=" << c.b
        << " z=" << c.z
        << " q=" << c.q
        << " e=";

    print_u128(bc.e);

    std::cout << " s0=";
    print_u128(bc.s0);

    std::cout << " m=";
    print_u128(bc.m);

    std::cout << " I=";
    print_u128(interval_count(data));

    std::cout << '\n';
}

/*
 * For a selected rank k:
 *
 *   k -> interval -> sampled n -> localize(n)
 *
 * must recover exactly k.
 *
 * We also test:
 *
 *   start - 1 = MISS
 *   start     = HIT(k)
 *   end       = HIT(k)
 *   end + 1   = MISS
 *
 * whenever those points exist.
 */
bool verify_rank_inverse(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 k,
    std::mt19937_64 &rng
) {
    Interval interval;

    if (!interval_from_rank(
            c,
            bc,
            data,
            k,
            interval)) {

        return false;
    }

    Result at_start =
        localize_n(
            c,
            bc,
            data,
            interval.start
        );

    if (!at_start.hit ||
        at_start.rank != k) {

        return false;
    }

    Result at_end =
        localize_n(
            c,
            bc,
            data,
            interval.end
        );

    if (!at_end.hit ||
        at_end.rank != k) {

        return false;
    }

    if (interval.start > bc.s0) {
        Result before =
            localize_n(
                c,
                bc,
                data,
                interval.start - 1
            );

        if (before.hit) {
            return false;
        }
    }

    if (interval.end < bc.m) {
        Result after =
            localize_n(
                c,
                bc,
                data,
                interval.end + 1
            );

        if (after.hit) {
            return false;
        }
    }

    /*
     * Random interior point.
     */
    u128 length =
        interval.end -
        interval.start +
        1;

    u128 offset =
        static_cast<u128>(rng()) %
        length;

    u128 n =
        interval.start +
        offset;

    Result interior =
        localize_n(
            c,
            bc,
            data,
            n
        );

    if (!interior.hit ||
        interior.rank != k) {

        return false;
    }

    return true;
}

bool deterministic_tests() {
    const std::vector<Case> cases = {
        {2, 0, 0, 0, 1},
        {2, 0, 0, 0, 3},
        {2, 1, 0, 0, 7},
        {2, 5, 0, 2, 1073741825ULL},

        {3, 0, 0, 0, 2},
        {3, 4, 1, 5, 987654321ULL},

        {5, 0, 0, 0, 4},
        {5, 3, 1, 4, 100},
        {5, 3, 1, 4, 1000007654321ULL},

        {7, 2, 3, 1, 123456789ULL},
        {11, 3, 4, 2, 987654321ULL}
    };

    std::mt19937_64 rng(
        0x247123456789ULL
    );

    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        BuiltCase bc =
            build_case(c);

        Data data =
            build_data(c);

        u128 I =
            interval_count(data);

        bool pass = true;

        /*
         * Exhaust every rank for small cases.
         */
        if (I <= static_cast<u128>(1000000)) {
            u64 count =
                static_cast<u64>(I);

            for (u64 k = 0;
                 k < count;
                 ++k) {

                if (!verify_rank_inverse(
                        c,
                        bc,
                        data,
                        static_cast<u128>(k),
                        rng)) {

                    pass = false;
                    break;
                }
            }

        } else {
            /*
             * Sample the first, last, and interior ranks.
             */
            std::vector<u128> ranks;

            ranks.push_back(0);

            if (I > 1) {
                ranks.push_back(1);
            }

            ranks.push_back(I / 2);

            if (I > 1) {
                ranks.push_back(I - 2);
            }

            ranks.push_back(I - 1);

            for (u128 k : ranks) {
                if (!verify_rank_inverse(
                        c,
                        bc,
                        data,
                        k,
                        rng)) {

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

bool random_rank_tests() {
    std::mt19937_64 rng(
        0x247AAA123ULL
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

        u64 q;

        switch (i % 5) {
            case 0:
                q =
                    1 +
                    rng() % 1000;
                break;

            case 1:
                q =
                    1000000ULL +
                    rng() % 1000000ULL;
                break;

            case 2:
                q =
                    1000000000000ULL +
                    rng() % 1000000000000ULL;
                break;

            case 3:
                q =
                    1000000000000000ULL +
                    rng() % 1000000000000000ULL;
                break;

            default:
                q =
                    1 + rng();
                break;
        }

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
            build_data(c);

        u128 I =
            interval_count(data);

        if (I == 0) {
            continue;
        }

        u128 k;

        if (I <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            k =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(I)
                );

        } else {
            u128 random128 =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            k =
                random128 % I;
        }

        if (!verify_rank_inverse(
                c,
                bc,
                data,
                k,
                rng)) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANK FAILURE\n";

                print_case(c);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout << '\n';
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
 * Test that randomly selected intervals are mutually ordered:
 *
 * x_k < E_k < x_{k+1}
 *
 * and
 *
 * x_{k+1} = E_k + s0 + 1.
 */
bool adjacency_tests() {
    std::mt19937_64 rng(
        0x247555888ULL
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
            build_data(c);

        u128 I =
            interval_count(data);

        if (I <= 1) {
            continue;
        }

        u128 k;

        if (I <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            k =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(
                        I - 1
                    )
                );

        } else {
            u128 random128 =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            k =
                random128 % (I - 1);
        }

        Interval a;
        Interval b_interval;

        if (!interval_from_rank(
                c,
                bc,
                data,
                k,
                a)) {

            ++failures;
            continue;
        }

        if (!interval_from_rank(
                c,
                bc,
                data,
                k + 1,
                b_interval)) {

            ++failures;
            continue;
        }

        if (!(a.start < a.end + 1)) {
            ++failures;
            continue;
        }

        if (b_interval.start !=
            a.end + bc.s0 + 1) {

            ++failures;
            continue;
        }

        if (!(a.end < b_interval.start)) {
            ++failures;
            continue;
        }
    }

    std::cout
        << "adjacency_cases="
        << cases
        << " adjacency_failures="
        << failures
        << " adjacency_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 247\n";

    bool deterministic_ok =
        deterministic_tests();

    bool rank_ok =
        random_rank_tests();

    bool adjacency_ok =
        adjacency_tests();

    bool overall =
        deterministic_ok &&
        rank_ok &&
        adjacency_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 247\n";

    return overall ? 0 : 1;
}
