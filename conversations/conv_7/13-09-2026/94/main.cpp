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

    for (std::size_t i = 0; i < L; ++i) {
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

        j += digit * place;

        place *=
            static_cast<u128>(p);
    }

    return j;
}

/*
 * Direct n -> interval localization.
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
            false,
            0,
            0,
            0,
            0
        };
    }

    u128 t =
        (n - bc.s0) /
        bc.p_pow_e;

    const std::size_t L =
        data.q_digits.size();

    std::vector<u64> t_digits(L, 0);

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

    std::vector<u64> j_digits(L, 0);

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
                false,
                0,
                j,
                0,
                0
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
        (j / data.p_powers[r] + 1) *
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
 * Closed-form count C(n).
 */
u128 count_hits_formula(
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
 * Completely independent Lucas predicate:
 *
 * HIT iff n is not digitwise <= m in base p.
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
        std::max(nd.size(), md.size());

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 a =
            i < nd.size() ? nd[i] : 0;

        u64 b =
            i < md.size() ? md[i] : 0;

        if (a > b) {
            return true;
        }
    }

    return false;
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
 * Small exhaustive test.
 *
 * Only cases with m <= 10000 are exhaustive.
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

            u128 expected = 0;

            for (u64 n = 0;
                 n <= m;
                 ++n) {

                if (lucas_hit(
                        c,
                        bc,
                        static_cast<u128>(n)
                    )) {

                    ++expected;
                }

                u128 actual =
                    count_hits_formula(
                        c,
                        bc,
                        data,
                        static_cast<u128>(n)
                    );

                if (actual != expected) {
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
             * Large deterministic case:
             * test boundaries and random-like fixed points.
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

                u128 c0 =
                    count_hits_formula(
                        c,
                        bc,
                        data,
                        n
                    );

                u128 c1 =
                    n == 0
                        ? 0
                        : count_hits_formula(
                              c,
                              bc,
                              data,
                              n - 1
                          );

                u128 delta =
                    c0 - c1;

                u128 expected_delta =
                    lucas_hit(
                        c,
                        bc,
                        n
                    )
                        ? 1
                        : 0;

                if (delta != expected_delta) {
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
 * 10000 random points.
 *
 * This is O(log_p m) per test.
 */
bool random_local_tests() {
    std::mt19937_64 rng(
        0x248123456789ULL
    );

    const u64 cases = 10000;
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
                UINT64_MAX - 2
            )) {

            continue;
        }

        Data data =
            build_data(c);

        u64 m =
            static_cast<u64>(bc.m);

        u64 n =
            rng() % (m + 1);

        u128 current =
            count_hits_formula(
                c,
                bc,
                data,
                static_cast<u128>(n)
            );

        u128 previous =
            n == 0
                ? 0
                : count_hits_formula(
                      c,
                      bc,
                      data,
                      static_cast<u128>(n - 1)
                  );

        u128 delta =
            current - previous;

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
                    << "RANDOM FAILURE\n";

                print_case(c);

                std::cout
                    << "n="
                    << n
                    << '\n';

                std::cout
                    << "delta=";

                print_u128(delta);

                std::cout
                    << " expected=";

                print_u128(expected);

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
 * Boundary-specific testing.
 *
 * For random k:
 *
 *   x_k      -> HIT(k)
 *   E_k      -> HIT(k)
 *   x_k - 1  -> MISS
 *   E_k + 1  -> MISS
 */
bool interval_boundary_tests() {
    std::mt19937_64 rng(
        0x248777888ULL
    );

    const u64 cases = 10000;
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

        u128 j =
            unrank_j(
                k,
                c.p,
                data
            );

        std::size_t r =
            data.q_digits.size();

        u128 temp = j;

        for (std::size_t d = 0;
             d < data.q_digits.size();
             ++d) {

            u64 digit =
                static_cast<u64>(
                    temp % c.p
                );

            temp /= c.p;

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

        u128 successor =
            (j / data.p_powers[r] + 1) *
            data.p_powers[r];

        u128 start =
            bc.s0 +
            j * bc.p_pow_e;

        u128 end =
            successor *
            bc.p_pow_e -
            1;

        /*
         * Start.
         */
        if (!localize(
                c,
                bc,
                data,
                start
            ).hit) {

            ++failures;
            continue;
        }

        /*
         * End.
         */
        if (!localize(
                c,
                bc,
                data,
                end
            ).hit) {

            ++failures;
            continue;
        }

        /*
         * Before.
         */
        if (start > 0) {
            if (lucas_hit(
                    c,
                    bc,
                    start - 1
                )) {

                ++failures;
                continue;
            }
        }

        /*
         * After.
         */
        if (end < bc.m) {
            if (lucas_hit(
                    c,
                    bc,
                    end + 1
                )) {

                ++failures;
                continue;
            }
        }
    }

    std::cout
        << "boundary_cases="
        << cases
        << " boundary_failures="
        << failures
        << " boundary_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 248T\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_local_tests();

    bool boundary_ok =
        interval_boundary_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        boundary_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 248T\n";

    return overall ? 0 : 1;
}
