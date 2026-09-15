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

struct Data {
    std::vector<u64> q_digits;
    std::vector<u128> p_powers;
    std::vector<u128> weights;
};

struct BuiltCase {
    u128 e;
    u128 s0;
    u128 p_pow_e;
    u128 m;
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
        unsigned digit = static_cast<unsigned>(value % 10);
        s.push_back(static_cast<char>('0' + digit));
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

    u128 p_to_a = pow_u128(c.p, c.a0);

    bc.s0 =
        static_cast<u128>(c.b + 1) *
        p_to_a;

    bc.p_pow_e =
        pow_u128(c.p, static_cast<u64>(bc.e));

    bc.m =
        bc.s0 +
        static_cast<u128>(c.q) * bc.p_pow_e -
        1;

    return bc;
}

std::vector<u64> digits_base(u64 value, u64 p) {
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

    data.q_digits = digits_base(c.q, c.p);

    data.p_powers.resize(data.q_digits.size());

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        data.p_powers[i] =
            pow_u128(c.p, static_cast<u64>(i));
    }

    data.weights.resize(data.q_digits.size());

    u128 weight = 1;

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        data.weights[i] = weight;

        weight *=
            static_cast<u128>(
                data.q_digits[i] + 1
            );
    }

    return data;
}

u128 interval_count(const Data &data) {
    u128 product = 1;

    for (u64 digit : data.q_digits) {
        product *=
            static_cast<u128>(digit + 1);
    }

    return product - 1;
}

/*
 * Reference:
 *
 * n -> t -> largest admissible j <= t
 *   -> interval membership
 *   -> rank(j)
 *
 * This explicitly constructs j.
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

    std::vector<u64> t_digits(
        data.q_digits.size(),
        0
    );

    u128 temp = t;

    for (std::size_t i = 0;
         i < t_digits.size();
         ++i) {

        t_digits[i] =
            static_cast<u64>(temp % c.p);

        temp /= c.p;
    }

    bool violation = false;
    std::size_t h = 0;

    for (std::size_t i = data.q_digits.size();
         i-- > 0;) {

        if (t_digits[i] > data.q_digits[i]) {
            violation = true;
            h = i;
            break;
        }
    }

    std::vector<u64> j_digits(
        data.q_digits.size(),
        0
    );

    for (std::size_t i = 0;
         i < j_digits.size();
         ++i) {

        if (violation && i <= h) {
            j_digits[i] =
                data.q_digits[i];
        } else {
            j_digits[i] =
                t_digits[i];
        }
    }

    u128 j = 0;

    for (std::size_t i = 0;
         i < j_digits.size();
         ++i) {

        j +=
            static_cast<u128>(j_digits[i]) *
            data.p_powers[i];
    }

    std::size_t r =
        data.q_digits.size();

    for (std::size_t i = 0;
         i < j_digits.size();
         ++i) {

        if (j_digits[i] <
            data.q_digits[i]) {

            r = i;
            break;
        }
    }

    if (r == data.q_digits.size()) {
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

    for (std::size_t i = 0;
         i < j_digits.size();
         ++i) {

        rank +=
            static_cast<u128>(j_digits[i]) *
            data.weights[i];
    }

    return {true, rank};
}

/*
 * Direct:
 *
 * n -> t
 *
 * Determine the most significant digit violation.
 * Accumulate rank directly from the resulting j-digits.
 *
 * IMPORTANT:
 *
 * membership is checked in the original n-coordinate:
 *
 *     n < successor(j) * p^e
 *
 * NOT against t.
 */
Result direct_localize(
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

    std::vector<u64> t_digits(
        data.q_digits.size(),
        0
    );

    u128 temp = t;

    for (std::size_t i = 0;
         i < t_digits.size();
         ++i) {

        t_digits[i] =
            static_cast<u64>(temp % c.p);

        temp /= c.p;
    }

    bool violation = false;
    std::size_t h = 0;

    for (std::size_t i = data.q_digits.size();
         i-- > 0;) {

        if (t_digits[i] > data.q_digits[i]) {
            violation = true;
            h = i;
            break;
        }
    }

    /*
     * The rank can be computed without constructing j.
     */
    u128 rank = 0;

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        u64 j_digit;

        if (violation && i <= h) {
            j_digit = data.q_digits[i];
        } else {
            j_digit = t_digits[i];
        }

        rank +=
            static_cast<u128>(j_digit) *
            data.weights[i];
    }

    /*
     * Find r, the first digit where j_r < q_r.
     */
    std::size_t r =
        data.q_digits.size();

    for (std::size_t i = 0;
         i < data.q_digits.size();
         ++i) {

        u64 j_digit;

        if (violation && i <= h) {
            j_digit = data.q_digits[i];
        } else {
            j_digit = t_digits[i];
        }

        if (j_digit < data.q_digits[i]) {
            r = i;
            break;
        }
    }

    if (r == data.q_digits.size()) {
        return {false, 0};
    }

    /*
     * Compute floor(j / p^r) directly
     * from the digits above r.
     *
     * This is the quotient whose successor gives
     * the exclusive endpoint.
     */
    u128 high = 0;
    u128 place = 1;

    for (std::size_t i = r;
         i < data.q_digits.size();
         ++i) {

        u64 j_digit;

        if (violation && i <= h) {
            j_digit = data.q_digits[i];
        } else {
            j_digit = t_digits[i];
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

    /*
     * Correct coordinate system:
     *
     * interval:
     *
     * [s0 + j*p^e,
     *  successor(j)*p^e - 1]
     *
     * Therefore:
     *
     * n < successor(j)*p^e
     */
    u128 endpoint_exclusive =
        successor * bc.p_pow_e;

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

    if (a.hit && a.rank != b.rank) {
        return false;
    }

    return true;
}

void print_case(const Case &c) {
    BuiltCase bc = build_case(c);
    Data data = build_data(c);

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
    std::vector<Case> cases = {
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

    std::cout << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        BuiltCase bc = build_case(c);
        Data data = build_data(c);

        std::vector<u128> points;

        points.push_back(0);

        if (bc.s0 > 0) {
            points.push_back(bc.s0 - 1);
        }

        points.push_back(bc.s0);
        points.push_back(bc.m);

        if (bc.m > 0) {
            points.push_back(bc.m - 1);
        }

        if (bc.m <= 50000) {
            for (u64 n = 0;
                 n <= static_cast<u64>(bc.m);
                 ++n) {

                points.push_back(
                    static_cast<u128>(n)
                );
            }
        } else {
            const u64 samples[] = {
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
                if (static_cast<u128>(n) <= bc.m) {
                    points.push_back(
                        static_cast<u128>(n)
                    );
                }
            }
        }

        bool pass = true;

        for (u128 n : points) {
            Result ref =
                reference_localize(
                    c, bc, data, n
                );

            Result direct =
                direct_localize(
                    c, bc, data, n
                );

            if (!same_result(ref, direct)) {
                pass = false;
                break;
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

bool random_arbitrary_tests() {
    std::mt19937_64 rng(
        0x245123456789ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0; i < cases; ++i) {
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
                q = 1 + rng() % 1000;
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
                q = 1 + rng();
                break;
        }

        Case c{p, a0, b, z, q};

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

        u64 n64 =
            rng() % (m + 1);

        u128 n =
            static_cast<u128>(n64);

        Result ref =
            reference_localize(
                c, bc, data, n
            );

        Result direct =
            direct_localize(
                c, bc, data, n
            );

        if (!same_result(ref, direct)) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n64
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

bool random_hit_tests() {
    std::mt19937_64 rng(
        0x245AAA555ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    u64 generated = 0;

    while (generated < cases) {
        const u64 choices[] = {
            2, 3, 5, 7, 11
        };

        u64 p =
            choices[rng() % 5];

        u64 a0 =
            rng() % 7;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 7;

        u64 q;

        switch (generated % 4) {
            case 0:
                q =
                    1 +
                    rng() % 10000;
                break;

            case 1:
                q =
                    1000000ULL +
                    rng() % 1000000ULL;
                break;

            case 2:
                q =
                    1000000000ULL +
                    rng() % 1000000000ULL;
                break;

            default:
                q =
                    1 +
                    rng() %
                    1000000000000ULL;
                break;
        }

        Case c{p, a0, b, z, q};

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

        u128 I =
            interval_count(data);

        if (I == 0 ||
            I >
                static_cast<u128>(
                    UINT64_MAX
                )) {
            continue;
        }

        u64 i_count =
            static_cast<u64>(I);

        u64 k =
            rng() % i_count;

        /*
         * Unrank k.
         */
        u64 remaining = k;

        u128 j = 0;
        u128 place = 1;

        for (std::size_t i = 0;
             i < data.q_digits.size();
             ++i) {

            u64 radix =
                data.q_digits[i] + 1;

            u64 digit =
                remaining % radix;

            remaining /= radix;

            j +=
                static_cast<u128>(digit) *
                place;

            place *=
                static_cast<u128>(c.p);
        }

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
            continue;
        }

        u128 next_j =
            (j / data.p_powers[r] + 1) *
            data.p_powers[r];

        u128 start =
            bc.s0 +
            j * bc.p_pow_e;

        u128 end =
            next_j * bc.p_pow_e - 1;

        if (end < start ||
            end > bc.m) {
            continue;
        }

        u128 n =
            start +
            static_cast<u128>(
                rng()
            ) %
                (end - start + 1);

        Result ref =
            reference_localize(
                c, bc, data, n
            );

        Result direct =
            direct_localize(
                c, bc, data, n
            );

        if (!ref.hit ||
            ref.rank !=
                static_cast<u128>(k) ||
            !same_result(ref, direct)) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "HIT FAILURE\n";

                print_case(c);

                std::cout
                    << "k="
                    << k
                    << '\n';

                std::cout
                    << "n=";

                print_u128(n);

                std::cout << '\n';
            }
        }

        ++generated;
    }

    std::cout
        << "hit_cases="
        << cases
        << " hit_failures="
        << failures
        << " hit_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

bool rank_consistency_tests() {
    std::mt19937_64 rng(
        0x245777888ULL
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
            rng() % 6;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 6;

        u64 q =
            1 +
            rng() % 1000000000ULL;

        Case c{
            p, a0, b, z, q
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

        u128 I =
            interval_count(data);

        if (I == 0 ||
            I >
                static_cast<u128>(
                    UINT64_MAX
                )) {
            continue;
        }

        u64 i_count =
            static_cast<u64>(I);

        u64 k =
            rng() % i_count;

        /*
         * Unrank k.
         */
        u64 remaining = k;

        u128 j = 0;
        u128 place = 1;

        for (std::size_t d = 0;
             d < data.q_digits.size();
             ++d) {

            u64 radix =
                data.q_digits[d] + 1;

            u64 digit =
                remaining % radix;

            remaining /= radix;

            j +=
                static_cast<u128>(digit) *
                place;

            place *=
                static_cast<u128>(c.p);
        }

        std::size_t r =
            data.q_digits.size();

        u128 temp_j = j;

        for (std::size_t d = 0;
             d < data.q_digits.size();
             ++d) {

            u64 digit =
                static_cast<u64>(
                    temp_j % c.p
                );

            temp_j /= c.p;

            if (digit <
                data.q_digits[d]) {

                r = d;
                break;
            }
        }

        if (r == data.q_digits.size()) {
            ++failures;
            continue;
        }

        u128 next_j =
            (j / data.p_powers[r] + 1) *
            data.p_powers[r];

        u128 start =
            bc.s0 +
            j * bc.p_pow_e;

        u128 end =
            next_j *
            bc.p_pow_e -
            1;

        if (end < start ||
            end > bc.m) {

            ++failures;
            continue;
        }

        u128 n =
            start +
            (end - start) / 2;

        Result direct =
            direct_localize(
                c, bc, data, n
            );

        if (!direct.hit ||
            direct.rank !=
                static_cast<u128>(k)) {

            ++failures;
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

int main() {
    std::cout
        << "START EXPERIMENT 245R\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_arbitrary_tests();

    bool hit_ok =
        random_hit_tests();

    bool rank_ok =
        rank_consistency_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        hit_ok &&
        rank_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 245R\n";

    return overall ? 0 : 1;
}