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

    std::size_t L =
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
 * Full reference implementation.
 *
 * Explicitly constructs the admissible j.
 */
Result reference_localize(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 n
) {
    if (n < bc.s0 || n > bc.m) {
        return {false, 0};
    }

    u128 t =
        (n - bc.s0) /
        bc.p_pow_e;

    std::size_t L =
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

    std::vector<u64> j_digits(L, 0);

    for (std::size_t i = 0; i < L; ++i) {
        if (violation && i <= h) {
            j_digits[i] =
                data.q_digits[i];
        } else {
            j_digits[i] =
                t_digits[i];
        }
    }

    u128 j = 0;

    for (std::size_t i = 0; i < L; ++i) {
        j +=
            static_cast<u128>(
                j_digits[i]
            ) *
            data.p_powers[i];
    }

    std::size_t r = L;

    for (std::size_t i = 0; i < L; ++i) {
        if (j_digits[i] <
            data.q_digits[i]) {

            r = i;
            break;
        }
    }

    if (r == L) {
        return {false, 0};
    }

    u128 successor =
        (j / data.p_powers[r] + 1) *
        data.p_powers[r];

    u128 endpoint_exclusive =
        successor * bc.p_pow_e;

    if (n >= endpoint_exclusive) {
        return {false, 0};
    }

    u128 rank = 0;

    for (std::size_t i = 0; i < L; ++i) {
        rank +=
            static_cast<u128>(
                j_digits[i]
            ) *
            data.weights[i];
    }

    return {true, rank};
}

/*
 * Compressed prefix-rank implementation.
 *
 * No j is constructed.
 *
 * If t <=_p q:
 *
 *     k = sum_i t_i W_i
 *
 * Otherwise let h be the most significant
 * digit satisfying t_h > q_h.
 *
 * Then
 *
 *     k =
 *       sum_{i>h} t_i W_i
 *       +
 *       sum_{i<=h} q_i W_i.
 *
 * The second sum is precomputed as q_prefix_rank[h].
 */
Result compressed_localize(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 n
) {
    if (n < bc.s0 || n > bc.m) {
        return {false, 0};
    }

    u128 t =
        (n - bc.s0) /
        bc.p_pow_e;

    std::size_t L =
        data.q_digits.size();

    /*
     * Extract t digits and locate the
     * most significant violation.
     */
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
        /*
         * t itself is admissible.
         */
        for (std::size_t i = 0; i < L; ++i) {
            rank +=
                static_cast<u128>(
                    t_digits[i]
                ) *
                data.weights[i];
        }

        /*
         * Since j=t, find the first
         * digit t_r < q_r.
         */
        std::size_t r = L;

        for (std::size_t i = 0; i < L; ++i) {
            if (t_digits[i] <
                data.q_digits[i]) {

                r = i;
                break;
            }
        }

        if (r == L) {
            return {false, 0};
        }

        /*
         * Compute successor(j)=successor(t).
         */
        u128 high = 0;
        u128 place = 1;

        for (std::size_t i = r;
             i < L;
             ++i) {

            high +=
                static_cast<u128>(
                    t_digits[i]
                ) *
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

    /*
     * Violating case:
     *
     * k =
     *
     *   q-prefix through h
     *   +
     *   t-suffix above h.
     */
    rank =
        data.q_prefix_rank[h];

    u128 suffix =
        0;

    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        suffix +=
            static_cast<u128>(
                t_digits[i]
            ) *
            data.weights[i];
    }

    rank += suffix;

    /*
     * For i <= h, j_i=q_i.
     *
     * Therefore the first digit r with
     * j_r < q_r must be the first digit
     * below h where q_i has not been
     * copied as a prefix.
     *
     * Since j_i=q_i for every i<=h,
     * r must be > h unless a higher
     * t-digit is already below q.
     *
     * Above h, j_i=t_i.
     */
    std::size_t r = L;

    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        if (t_digits[i] <
            data.q_digits[i]) {

            r = i;
            break;
        }
    }

    /*
     * It is possible that all higher digits
     * equal q. Then r does not occur above h,
     * which means the clamped j is terminal
     * and therefore does not label a HIT interval.
     */
    if (r == L) {
        return {false, 0};
    }

    /*
     * Compute floor(j/p^r) directly.
     */
    u128 high = 0;
    u128 place = 1;

    for (std::size_t i = r;
         i < L;
         ++i) {

        u64 j_digit;

        if (i <= h) {
            j_digit =
                data.q_digits[i];
        } else {
            j_digit =
                t_digits[i];
        }

        high +=
            static_cast<u128>(j_digit) *
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

    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        BuiltCase bc =
            build_case(c);

        Data data =
            build_data(c);

        bool pass = true;

        if (bc.m <=
            static_cast<u128>(50000)) {

            u64 m =
                static_cast<u64>(bc.m);

            for (u64 n = 0;
                 n <= m;
                 ++n) {

                Result ref =
                    reference_localize(
                        c,
                        bc,
                        data,
                        static_cast<u128>(n)
                    );

                Result compressed =
                    compressed_localize(
                        c,
                        bc,
                        data,
                        static_cast<u128>(n)
                    );

                if (!same_result(
                        ref,
                        compressed)) {

                    pass = false;
                    break;
                }
            }

        } else {
            const u64 samples[] = {
                0,
                1,
                2,
                3,
                5,
                7,
                16,
                32,
                100,
                1000,
                10000
            };

            for (u64 n : samples) {
                if (static_cast<u128>(n) > bc.m) {
                    continue;
                }

                Result ref =
                    reference_localize(
                        c,
                        bc,
                        data,
                        static_cast<u128>(n)
                    );

                Result compressed =
                    compressed_localize(
                        c,
                        bc,
                        data,
                        static_cast<u128>(n)
                    );

                if (!same_result(
                        ref,
                        compressed)) {

                    pass = false;
                    break;
                }
            }

            /*
             * Explicitly test the final points.
             */
            const u128 special[] = {
                bc.s0,
                bc.s0 + 1,
                bc.m - 1,
                bc.m
            };

            for (u128 n : special) {
                if (n > bc.m) {
                    continue;
                }

                Result ref =
                    reference_localize(
                        c,
                        bc,
                        data,
                        n
                    );

                Result compressed =
                    compressed_localize(
                        c,
                        bc,
                        data,
                        n
                    );

                if (!same_result(
                        ref,
                        compressed)) {

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
        0x246123456789ULL
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

        if (bc.m >
            static_cast<u128>(
                UINT64_MAX
            )) {

            continue;
        }

        Data data =
            build_data(c);

        u64 m =
            static_cast<u64>(bc.m);

        u64 n =
            rng() % (m + 1);

        Result ref =
            reference_localize(
                c,
                bc,
                data,
                static_cast<u128>(n)
            );

        Result compressed =
            compressed_localize(
                c,
                bc,
                data,
                static_cast<u128>(n)
            );

        if (!same_result(
                ref,
                compressed)) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
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

bool direct_formula_tests() {
    std::mt19937_64 rng(
        0x246555AAAULL
    );

    const u64 cases = 100000;
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
            rng() % 1000000000000ULL;

        Case c{
            p,
            a0,
            b,
            z,
            q
        };

        BuiltCase bc =
            build_case(c);

        if (bc.m >
            static_cast<u128>(
                UINT64_MAX
            )) {

            continue;
        }

        Data data =
            build_data(c);

        u64 m =
            static_cast<u64>(bc.m);

        u64 n =
            rng() % (m + 1);

        Result ref =
            reference_localize(
                c,
                bc,
                data,
                static_cast<u128>(n)
            );

        Result compressed =
            compressed_localize(
                c,
                bc,
                data,
                static_cast<u128>(n)
            );

        if (!same_result(
                ref,
                compressed)) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "FORMULA FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
                    << '\n';
            }
        }
    }

    std::cout
        << "formula_cases="
        << cases
        << " formula_failures="
        << failures
        << " formula_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 246\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool formula_ok =
        direct_formula_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        formula_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 246\n";

    return overall ? 0 : 1;
}
