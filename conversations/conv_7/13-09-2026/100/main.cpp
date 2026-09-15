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
    u128 count = 1;

    for (u64 digit : data.m_digits) {
        count *=
            static_cast<u128>(
                digit + 1
            );
    }

    return count;
}

u128 unrank_miss(
    u128 k,
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
            k % radix;

        k /= radix;

        value +=
            digit * place;

        place *=
            static_cast<u128>(p);
    }

    return value;
}

/*
 * Mixed-radix carry position when incrementing k.
 *
 * r is the first digit for which k_r < m_r.
 *
 * Digits below r are at their maxima and reset to 0.
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
 * r from the q-based p-adic interval construction:
 *
 * j_r < q_r,
 * with lower digits equal to q's lower digits
 * after the mixed-radix carry interpretation.
 *
 * For a nonempty HIT gap, this is exactly the
 * carry position.
 */
std::size_t q_gap_type(
    const Case &c,
    const BuiltCase &bc,
    u128 k
) {
    u128 q =
        (bc.m + 1 - bc.s0) /
        bc.p_pow_e;

    std::vector<u64> q_digits =
        digits_base(q, c.p);

    u128 j = 0;
    u128 place = 1;
    u128 remaining = k;

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

        j +=
            digit * place;

        place *=
            static_cast<u128>(c.p);
    }

    for (std::size_t i = 0;
         i < q_digits.size();
         ++i) {

        u64 qd =
            q_digits[i];

        u64 jd =
            static_cast<u64>(
                j /
                pow_u128(
                    c.p,
                    static_cast<u64>(i)
                ) %
                c.p
            );

        if (jd < qd) {
            return i;
        }
    }

    return q_digits.size();
}

u128 predicted_gap_length(
    const Case &c,
    const BuiltCase &bc,
    std::size_t r
) {
    u128 q =
        (bc.m + 1 - bc.s0) /
        bc.p_pow_e;

    u128 p_r =
        pow_u128(
            c.p,
            static_cast<u64>(r)
        );

    u128 q_mod =
        q % p_r;

    return
        (p_r - q_mod) *
        bc.p_pow_e -
        bc.s0;
}

/*
 * Direct gap from consecutive MISS numbers.
 */
u128 actual_gap_length(
    u128 left,
    u128 right
) {
    if (right <= left + 1) {
        return 0;
    }

    return right - left - 1;
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

        /*
         * Exhaust small rank spaces.
         */
        if (count <=
            static_cast<u128>(100000)) {

            for (u64 k = 0;
                 k + 1 <
                 static_cast<u64>(count);
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
                    actual_gap_length(
                        left,
                        right
                    );

                u128 predicted =
                    predicted_gap_length(
                        c,
                        bc,
                        carry.r
                    );

                /*
                 * Zero gap means the two MISS values
                 * are consecutive and there is no HIT.
                 */
                if (actual != 0 &&
                    actual != predicted) {

                    pass = false;
                    break;
                }

                if (actual == 0 &&
                    predicted != 0) {

                    /*
                     * Only compare the p-adic HIT length
                     * when the carry actually crosses a
                     * MISS-block boundary.
                     */
                    u128 block_position =
                        static_cast<u128>(k + 1) %
                        bc.s0;

                    if (block_position == 0) {
                        pass = false;
                        break;
                    }
                }

                if (actual != 0) {
                    std::size_t q_r =
                        q_gap_type(
                            c,
                            bc,
                            static_cast<u128>(k)
                        );

                    if (q_r != carry.r) {
                        pass = false;
                        break;
                    }
                }
            }

        } else {
            /*
             * Large cases: sample ranks.
             *
             * Explicitly target block boundaries because
             * those are the nonzero HIT gaps.
             */
            std::vector<u128> ranks;

            ranks.push_back(0);
            ranks.push_back(1);

            if (bc.s0 > 0) {
                ranks.push_back(
                    bc.s0 - 1
                );

                ranks.push_back(
                    bc.s0 * 2 - 1
                );
            }

            ranks.push_back(
                count / 2
            );

            if (count > 2) {
                ranks.push_back(
                    count - 2
                );
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
                    actual_gap_length(
                        left,
                        right
                    );

                if (actual != 0) {
                    u128 predicted =
                        predicted_gap_length(
                            c,
                            bc,
                            carry.r
                        );

                    if (actual != predicted) {
                        pass = false;
                        break;
                    }

                    if (q_gap_type(
                            c,
                            bc,
                            k
                        ) != carry.r) {

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
 * Random carry-vs-gap tests
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
         * Deliberately bias samples toward
         * possible MISS-block endings.
         */
        if ((i % 3) == 0 &&
            bc.s0 > 0) {

            u128 block =
                1 +
                static_cast<u128>(
                    rng() % 1000000ULL
                );

            k =
                block * bc.s0 - 1;

            if (k + 1 >= count) {
                k =
                    static_cast<u128>(
                        rng()
                    ) % (count - 1);
            }

        } else {
            if (count - 1 <=
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
                u128 r =
                    (static_cast<u128>(rng()) << 64) |
                    static_cast<u128>(rng());

                k =
                    r % (count - 1);
            }
        }

        CarryInfo carry =
            carry_position(k, data);

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
            actual_gap_length(
                left,
                right
            );

        if (actual == 0) {
            continue;
        }

        u128 predicted =
            predicted_gap_length(
                c,
                bc,
                carry.r
            );

        std::size_t q_r =
            q_gap_type(
                c,
                bc,
                k
            );

        if (actual != predicted ||
            q_r != carry.r) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "RANDOM FAILURE\n";

                print_case(c);

                std::cout
                    << "k=";

                print_u128(k);

                std::cout
                    << " r=";

                print_u128(
                    static_cast<u128>(
                        carry.r
                    )
                );

                std::cout
                    << " actual=";

                print_u128(actual);

                std::cout
                    << " predicted=";

                print_u128(predicted);

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
 * Carry transformation test
 * --------------------------------------------------------------------------
 *
 * Explicitly verify that incrementing mixed-radix rank k:
 *
 *   - resets digits below r to zero;
 *   - increments digit r;
 *   - leaves higher digits unchanged.
 *
 * Then compare the induced MISS-number difference
 * against the p-adic gap formula.
 */
bool carry_structure_tests() {
    std::mt19937_64 rng(
        0x253777888ULL
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

        u128 k =
            static_cast<u128>(
                rng()
            ) % (count - 1);

        CarryInfo carry =
            carry_position(k, data);

        if (carry.overflow) {
            ++failures;
            continue;
        }

        /*
         * Extract mixed-radix digits of k and k+1.
         */
        std::size_t L =
            data.m_digits.size();

        std::vector<u64> before(L, 0);
        std::vector<u64> after(L, 0);

        u128 kb = k;
        u128 ka = k + 1;

        for (std::size_t d = 0;
             d < L;
             ++d) {

            u128 radix =
                static_cast<u128>(
                    data.m_digits[d] + 1
                );

            before[d] =
                static_cast<u64>(
                    kb % radix
                );

            after[d] =
                static_cast<u64>(
                    ka % radix
                );

            kb /= radix;
            ka /= radix;
        }

        bool structure_ok = true;

        if (after[carry.r] !=
            before[carry.r] + 1) {

            structure_ok = false;
        }

        for (std::size_t d = 0;
             d < carry.r;
             ++d) {

            if (before[d] !=
                    data.m_digits[d] ||
                after[d] != 0) {

                structure_ok = false;
                break;
            }
        }

        for (std::size_t d = carry.r + 1;
             d < L;
             ++d) {

            if (after[d] != before[d]) {
                structure_ok = false;
                break;
            }
        }

        if (!structure_ok) {
            ++failures;
            continue;
        }

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

        u128 actual =
            actual_gap_length(
                left,
                right
            );

        /*
         * If the carry occurs at the end of a MISS block,
         * this gap must be the HIT interval of type r.
         */
        if (actual != 0) {
            u128 predicted =
                predicted_gap_length(
                    c,
                    bc,
                    carry.r
                );

            if (actual != predicted) {
                ++failures;
            }
        }
    }

    std::cout
        << "carry_cases="
        << cases
        << " carry_failures="
        << failures
        << " carry_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 253\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool carry_ok =
        carry_structure_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        carry_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 253\n";

    return overall ? 0 : 1;
}
