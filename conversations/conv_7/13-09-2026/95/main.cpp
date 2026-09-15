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
};

struct Localization {
    bool hit;
    u128 rank;
    u128 j;
    u128 interval_start;
    u128 interval_end;
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

Data build_data(const Case &c) {
    Data data;

    data.q_digits =
        digits_base(
            static_cast<u128>(c.q),
            c.p
        );

    const std::size_t L =
        data.q_digits.size();

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
 * Existing localization reference.
 */
Localization localize(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 n
) {
    if (n < bc.s0 ||
        n > bc.m) {

        return {
            false, 0, 0, 0, 0
        };
    }

    u128 t =
        (n - bc.s0) /
        bc.p_pow_e;

    const std::size_t L =
        data.q_digits.size();

    std::vector<u64> t_digits(
        L,
        0
    );

    u128 temp = t;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        t_digits[i] =
            static_cast<u64>(
                temp % c.p
            );

        temp /= c.p;
    }

    bool violation = false;
    std::size_t h = 0;

    for (std::size_t i = L;
         i-- > 0;) {

        if (t_digits[i] >
            data.q_digits[i]) {

            violation = true;
            h = i;
            break;
        }
    }

    std::vector<u64> j_digits(
        L,
        0
    );

    for (std::size_t i = 0;
         i < L;
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
    u128 rank = 0;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        j +=
            static_cast<u128>(
                j_digits[i]
            ) *
            data.p_powers[i];

        rank +=
            static_cast<u128>(
                j_digits[i]
            ) *
            data.weights[i];
    }

    bool terminal = true;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        if (j_digits[i] !=
            data.q_digits[i]) {

            terminal = false;
            break;
        }
    }

    if (terminal) {
        u128 I =
            interval_count(data);

        if (I == 0) {
            return {
                false, 0, j, 0, 0
            };
        }

        u128 last_j =
            unrank_j(
                I - 1,
                c.p,
                data
            );

        std::size_t r =
            L;

        u128 temp_j = last_j;

        for (std::size_t i = 0;
             i < L;
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

        u128 successor =
            (last_j /
             data.p_powers[r] + 1) *
            data.p_powers[r];

        u128 start =
            bc.s0 +
            last_j * bc.p_pow_e;

        u128 end =
            successor *
            bc.p_pow_e -
            1;

        return {
            false,
            I - 1,
            last_j,
            start,
            end
        };
    }

    std::size_t r = L;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        if (j_digits[i] <
            data.q_digits[i]) {

            r = i;
            break;
        }
    }

    if (r == L) {
        return {
            false,
            rank,
            j,
            0,
            0
        };
    }

    u128 successor =
        (j /
         data.p_powers[r] + 1) *
        data.p_powers[r];

    u128 start =
        bc.s0 +
        j * bc.p_pow_e;

    u128 end =
        successor *
        bc.p_pow_e -
        1;

    bool hit =
        n >= start &&
        n <= end;

    return {
        hit,
        rank,
        j,
        start,
        end
    };
}

/*
 * Existing interval/rank counting formula.
 */
u128 count_hits_localized(
    const Case &c,
    const BuiltCase &bc,
    const Data &data,
    u128 n
) {
    if (n < bc.s0) {
        return 0;
    }

    if (n > bc.m) {
        n = bc.m;
    }

    Localization loc =
        localize(
            c,
            bc,
            data,
            n
        );

    if (loc.hit) {
        return
            n + 1 -
            bc.s0 *
            (loc.rank + 1);
    }

    u128 I =
        interval_count(data);

    if (I == 0) {
        return 0;
    }

    return
        loc.interval_end + 1 -
        bc.s0 *
        (loc.rank + 1);
}

/*
 * NEW:
 *
 * Direct digit formula for
 *
 * M(n) = #{x <= n : x <=_p m}.
 *
 * Then
 *
 * C(n) = n+1-M(n).
 *
 * This does not use localization, ranks,
 * intervals, or s0/q decomposition.
 */
u128 count_hits_digits(
    const Case &c,
    const BuiltCase &bc,
    u128 n
) {
    if (n == static_cast<u128>(-1)) {
        return 0;
    }

    if (n > bc.m) {
        n = bc.m;
    }

    std::vector<u64> n_digits =
        digits_base(n, c.p);

    std::vector<u64> m_digits =
        digits_base(bc.m, c.p);

    const std::size_t L =
        std::max(
            n_digits.size(),
            m_digits.size()
        );

    /*
     * Number of assignments for the
     * lower suffix satisfying x_i <= m_i.
     */
    u128 suffix_choices = 1;

    /*
     * We need the product
     *
     *   prod_{j<i} (m_j+1)
     *
     * while scanning from MSB to LSB.
     *
     * Build all suffix products first.
     */
    std::vector<u128> suffix_product(
        L + 1,
        1
    );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 md =
            i < m_digits.size()
                ? m_digits[i]
                : 0;

        suffix_product[i + 1] =
            suffix_product[i] *
            static_cast<u128>(md + 1);
    }

    u128 misses = 0;

    /*
     * Scan from MSB to LSB.
     *
     * At position i:
     *
     *   n_i choices are available for x_i<n_i
     *   if n_i <= m_i.
     *
     *   If n_i > m_i, there are m_i+1
     *   admissible values, and the tight path dies.
     */
    for (std::size_t pos = L;
         pos-- > 0;) {

        u64 nd =
            pos < n_digits.size()
                ? n_digits[pos]
                : 0;

        u64 md =
            pos < m_digits.size()
                ? m_digits[pos]
                : 0;

        u128 lower_suffix =
            suffix_product[pos];

        if (nd <= md) {
            misses +=
                static_cast<u128>(nd) *
                lower_suffix;
        } else {
            misses +=
                static_cast<u128>(md + 1) *
                lower_suffix;

            /*
             * No x with the same prefix can
             * remain equal to n.
             */
            return
                n + 1 - misses;
        }
    }

    /*
     * n itself is valid because every digit
     * satisfies n_i <= m_i.
     */
    misses += 1;

    return
        n + 1 - misses;
}

/*
 * Independent Lucas predicate.
 */
bool lucas_hit(
    const Case &c,
    const BuiltCase &bc,
    u128 n
) {
    if (n > bc.m) {
        return false;
    }

    std::vector<u64> nd =
        digits_base(n, c.p);

    std::vector<u64> md =
        digits_base(bc.m, c.p);

    const std::size_t L =
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
            return true;
        }
    }

    return false;
}

u128 brute_lucas_count(
    const Case &c,
    const BuiltCase &bc,
    u64 n
) {
    u128 count = 0;

    for (u64 x = 0;
         x <= n;
         ++x) {

        if (lucas_hit(
                c,
                bc,
                static_cast<u128>(x)
            )) {

            ++count;
        }
    }

    return count;
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
 * Exhaustive cases small enough for
 * direct Lucas counting.
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

        Data data =
            build_data(c);

        bool pass = true;

        if (bc.m <=
            static_cast<u128>(10000)) {

            u64 m =
                static_cast<u64>(bc.m);

            for (u64 n = 0;
                 n <= m;
                 ++n) {

                u128 direct =
                    count_hits_digits(
                        c,
                        bc,
                        static_cast<u128>(n)
                    );

                u128 reference =
                    brute_lucas_count(
                        c,
                        bc,
                        n
                    );

                if (direct != reference) {
                    pass = false;

                    std::cout
                        << "failure_n="
                        << n
                        << '\n';

                    break;
                }
            }

        } else {
            /*
             * Large case: compare the new
             * digit formula with the
             * established localized formula.
             */
            const u128 points[] = {
                0,
                1,
                bc.s0 - 1,
                bc.s0,
                bc.s0 + 1,
                bc.m - 1,
                bc.m
            };

            for (u128 n : points) {
                if (n > bc.m) {
                    continue;
                }

                u128 direct =
                    count_hits_digits(
                        c,
                        bc,
                        n
                    );

                u128 localized =
                    count_hits_localized(
                        c,
                        bc,
                        data,
                        n
                    );

                if (direct != localized) {
                    pass = false;

                    std::cout
                        << "failure_n=";

                    print_u128(n);

                    std::cout
                        << '\n';

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
 * Random formula-vs-localization test.
 *
 * 100000 cases, O(log_p n) each.
 */
bool random_tests() {
    std::mt19937_64 rng(
        0x249123456789ULL
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

        u128 direct =
            count_hits_digits(
                c,
                bc,
                static_cast<u128>(n)
            );

        u128 localized =
            count_hits_localized(
                c,
                bc,
                data,
                static_cast<u128>(n)
            );

        if (direct != localized) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
                    << '\n';

                std::cout
                    << "direct=";

                print_u128(direct);

                std::cout
                    << " localized=";

                print_u128(localized);

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
 * Random small independent Lucas tests.
 *
 * This is the strongest independent check:
 *
 *   digit_count_formula
 *       vs
 *   direct Lucas enumeration
 */
bool independent_small_tests() {
    std::mt19937_64 rng(
        0x249777888ULL
    );

    const u64 cases = 10000;
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
            rng() % 4;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 4;

        u64 q =
            1 +
            rng() % 100;

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
            static_cast<u128>(50000)) {

            continue;
        }

        u64 m =
            static_cast<u64>(bc.m);

        u64 n =
            rng() % (m + 1);

        u128 direct =
            count_hits_digits(
                c,
                bc,
                static_cast<u128>(n)
            );

        u128 reference =
            brute_lucas_count(
                c,
                bc,
                n
            );

        if (direct != reference) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "INDEPENDENT FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
                    << '\n';
            }
        }
    }

    std::cout
        << "independent_cases="
        << cases
        << " independent_failures="
        << failures
        << " independent_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * Monotonicity and one-step derivative:
 *
 * C(n)-C(n-1) = HIT(n).
 */
bool derivative_tests() {
    std::mt19937_64 rng(
        0x249424242ULL
    );

    const u64 cases = 50000;
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

        if (bc.m == 0 ||
            bc.m >
                static_cast<u128>(
                    UINT64_MAX
                )) {

            continue;
        }

        u64 m =
            static_cast<u64>(bc.m);

        u64 n =
            1 +
            rng() % m;

        u128 now =
            count_hits_digits(
                c,
                bc,
                static_cast<u128>(n)
            );

        u128 previous =
            count_hits_digits(
                c,
                bc,
                static_cast<u128>(n - 1)
            );

        u128 delta =
            now - previous;

        u128 expected =
            lucas_hit(
                c,
                bc,
                static_cast<u128>(n)
            )
                ? 1
                : 0;

        if (delta != expected) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "DERIVATIVE FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
                    << '\n';
            }
        }
    }

    std::cout
        << "derivative_cases="
        << cases
        << " derivative_failures="
        << failures
        << " derivative_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 249\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool independent_ok =
        independent_small_tests();

    bool derivative_ok =
        derivative_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        independent_ok &&
        derivative_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 249\n";

    return overall ? 0 : 1;
}
