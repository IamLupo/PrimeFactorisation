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

        data.weights[i] =
            weight;

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

/*
 * First mixed-radix digit that is not at its maximum.
 *
 * This is the carry position when k -> k+1.
 */
CarryInfo carry_position(
    u128 k,
    const Data &data
) {
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
 * General gap formula from the full base-p
 * digits of m.
 *
 * If carry occurs at r then every lower digit
 * is at its maximum and resets to zero.
 *
 * Therefore
 *
 *   M_{k+1} - M_k
 *     = p^r - sum_{i<r} m_i p^i.
 *
 * HIT length is one less.
 */
u128 predicted_gap_from_m(
    const Case &c,
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

    u128 p_r =
        data.p_powers[r];

    return
        p_r - lower - 1;
}

/*
 * Original q-based formula.
 *
 * Only valid after translating
 *
 *   r_full = e + r_q.
 */
u128 predicted_gap_from_q(
    const Case &c,
    const BuiltCase &bc,
    std::size_t r_full
) {
    if (
        r_full <
        static_cast<std::size_t>(bc.e)
    ) {
        return 0;
    }

    std::size_t rq =
        r_full -
        static_cast<std::size_t>(bc.e);

    u128 q =
        (bc.m + 1 - bc.s0) /
        bc.p_pow_e;

    u128 p_rq =
        pow_u128(
            c.p,
            static_cast<u64>(rq)
        );

    u128 q_mod =
        q % p_rq;

    return
        (p_rq - q_mod) *
        bc.p_pow_e -
        bc.s0;
}

/*
 * Direct actual gap.
 */
u128 actual_gap(
    u128 left,
    u128 right
) {
    if (right <= left + 1) {
        return 0;
    }

    return right - left - 1;
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

/*
 * --------------------------------------------------------------------------
 * Deterministic tests
 * --------------------------------------------------------------------------
 */
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

        std::vector<u128> ranks;

        ranks.push_back(0);

        if (count > 1) {
            ranks.push_back(1);
        }

        if (count > 4) {
            ranks.push_back(count / 4);
            ranks.push_back(count / 2);
            ranks.push_back(count - 2);
        }

        for (u128 k : ranks) {
            if (k + 1 >= count) {
                continue;
            }

            CarryInfo carry =
                carry_position(
                    k,
                    data
                );

            if (carry.overflow) {
                pass = false;
                break;
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

            u128 actual =
                actual_gap(
                    left,
                    right
                );

            u128 general =
                predicted_gap_from_m(
                    c,
                    data,
                    carry.r
                );

            if (actual != general) {
                pass = false;
                break;
            }

            /*
             * The q formula is compared only when
             * the gap is nonempty.
             */
            if (actual != 0) {
                u128 q_formula =
                    predicted_gap_from_q(
                        c,
                        bc,
                        carry.r
                    );

                if (actual != q_formula) {
                    pass = false;
                    break;
                }

                /*
                 * HIT-type index is:
                 *
                 *   r_q = r_full - e
                 */
                if (
                    carry.r <
                    static_cast<std::size_t>(
                        bc.e
                    )
                ) {
                    pass = false;
                    break;
                }
            }
        }

        /*
         * Exhaust small rank spaces.
         */
        if (pass &&
            count <=
                static_cast<u128>(50000)) {

            for (u64 k = 0;
                 k + 1 <
                 static_cast<u64>(count);
                 ++k) {

                CarryInfo carry =
                    carry_position(
                        static_cast<u128>(k),
                        data
                    );

                u128 left =
                    unrank_miss(
                        static_cast<u128>(k),
                        c.p,
                        data
                    );

                u128 right =
                    unrank_miss(
                        static_cast<u128>(k + 1),
                        c.p,
                        data
                    );

                u128 actual =
                    actual_gap(
                        left,
                        right
                    );

                u128 general =
                    predicted_gap_from_m(
                        c,
                        data,
                        carry.r
                    );

                if (actual != general) {
                    pass = false;
                    break;
                }

                if (actual != 0) {
                    u128 q_formula =
                        predicted_gap_from_q(
                            c,
                            bc,
                            carry.r
                        );

                    if (actual != q_formula) {
                        pass = false;
                        break;
                    }
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
 * Random tests
 * --------------------------------------------------------------------------
 */
bool random_tests() {
    std::mt19937_64 rng(
        0x253123456789ULL
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

        /*
         * Strongly sample MISS-block boundaries.
         */
        if (
            bc.s0 > 0 &&
            i % 2 == 0
        ) {
            u128 blocks =
                count / bc.s0;

            if (blocks > 0) {
                u128 block =
                    static_cast<u128>(
                        rng()
                    ) % blocks;

                k =
                    (block + 1) *
                    bc.s0 - 1;

                if (k + 1 >= count) {
                    k =
                        static_cast<u128>(
                            rng()
                        ) % (count - 1);
                }
            } else {
                k =
                    static_cast<u128>(
                        rng()
                    ) % (count - 1);
            }
        } else {
            k =
                static_cast<u128>(
                    rng()
                ) % (count - 1);
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
                c.p,
                data
            );

        u128 right =
            unrank_miss(
                k + 1,
                c.p,
                data
            );

        u128 actual =
            actual_gap(
                left,
                right
            );

        u128 general =
            predicted_gap_from_m(
                c,
                data,
                carry.r
            );

        if (actual != general) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "GENERAL FAILURE\n";

                print_case(c);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout
                    << " r=";

                std::cout
                    << carry.r;

                std::cout
                    << " actual=";

                print_u128(actual);

                std::cout
                    << " general=";

                print_u128(general);

                std::cout
                    << '\n';
            }

            continue;
        }

        /*
         * A nonempty gap is a HIT interval,
         * so the q-based formula should also agree.
         */
        if (actual != 0) {
            u128 q_formula =
                predicted_gap_from_q(
                    c,
                    bc,
                    carry.r
                );

            if (actual != q_formula) {
                ++failures;

                if (failures <= 5) {
                    std::cout
                        << "Q FORMULA FAILURE\n";

                    print_case(c);

                    std::cout
                        << "k=";

                    print_u128(k);

                    std::cout
                        << " full_r=";

                    std::cout
                        << carry.r;

                    std::cout
                        << " q_r=";

                    std::cout
                        << (
                            carry.r -
                            static_cast<std::size_t>(
                                bc.e
                            )
                        );

                    std::cout
                        << " actual=";

                    print_u128(actual);

                    std::cout
                        << " q_formula=";

                    print_u128(q_formula);

                    std::cout
                        << '\n';
                }
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
 * Explicit structural digit test
 * --------------------------------------------------------------------------
 *
 * Verify that m has:
 *
 *   p-1 ... p-1 | b | 0 ... 0 | q_digits
 *    a digits      z zeros
 *
 * and therefore the full carry index is
 *
 *   r_full = e + r_q
 *
 * for nonempty HIT gaps.
 */
bool structural_tests() {
    std::mt19937_64 rng(
        0x253777888ULL
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

        Data data =
            build_data(c, bc);

        u128 count =
            miss_count(data);

        if (count <= 1) {
            continue;
        }

        /*
         * Pick a random rank immediately before a
         * potentially nonempty gap by selecting a
         * MISS block boundary.
         */
        u128 blocks =
            count / bc.s0;

        if (blocks == 0) {
            continue;
        }

        u128 block =
            static_cast<u128>(rng()) % blocks;

        u128 k =
            (block + 1) *
            bc.s0 - 1;

        if (k + 1 >= count) {
            continue;
        }

        CarryInfo carry =
            carry_position(
                k,
                data
            );

        if (carry.overflow ||
            carry.r <
            static_cast<std::size_t>(bc.e)) {

            ++failures;
            continue;
        }

        std::size_t rq =
            carry.r -
            static_cast<std::size_t>(bc.e);

        u128 q_formula =
            predicted_gap_from_q(
                c,
                bc,
                carry.r
            );

        u128 m_formula =
            predicted_gap_from_m(
                c,
                data,
                carry.r
            );

        if (q_formula != m_formula) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "STRUCTURAL FAILURE\n";

                print_case(c);

                std::cout
                    << "full_r=";

                std::cout << carry.r;

                std::cout
                    << " q_r=";

                std::cout << rq;

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "structural_cases="
        << cases
        << " structural_failures="
        << failures
        << " structural_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 253R\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool structural_ok =
        structural_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        structural_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 253R\n";

    return overall ? 0 : 1;
}
