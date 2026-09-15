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

/*
 * --------------------------------------------------------------------------
 * DIRECT CLOSED FORM
 * --------------------------------------------------------------------------
 *
 * M(n) = #{0 <= x <= n : x <=_p m}
 *
 * Let
 *
 *   P_i = prod_{j<i} (m_j + 1)
 *
 * and let h be the most significant index for which
 *
 *   n_h > m_h.
 *
 * Then:
 *
 *   M(n)
 *     = sum_{i>h} n_i P_i
 *       + (m_h+1) P_h.
 *
 * If no violating digit exists:
 *
 *   M(n)
 *     = 1 + sum_i n_i P_i.
 *
 * Therefore:
 *
 *   C(n) = n+1-M(n).
 */
u128 count_misses_closed_form(
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

    std::vector<u64> n_digits(
        L,
        0
    );

    std::vector<u64> m_digits(
        L,
        0
    );

    for (std::size_t i = 0; i < nd.size(); ++i) {
        n_digits[i] = nd[i];
    }

    for (std::size_t i = 0; i < md.size(); ++i) {
        m_digits[i] = md[i];
    }

    /*
     * P_i = product of lower m-digit choices.
     */
    std::vector<u128> P(
        L + 1,
        1
    );

    for (std::size_t i = 0; i < L; ++i) {
        P[i + 1] =
            P[i] *
            static_cast<u128>(
                m_digits[i] + 1
            );
    }

    /*
     * Find most significant violation.
     */
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
        /*
         * n <=_p m.
         * n itself is valid.
         */
        for (std::size_t i = 0;
             i < L;
             ++i) {

            misses +=
                static_cast<u128>(
                    n_digits[i]
                ) *
                P[i];
        }

        misses += 1;

        return misses;
    }

    /*
     * Higher digits remain tight.
     */
    for (std::size_t i = h + 1;
         i < L;
         ++i) {

        misses +=
            static_cast<u128>(
                n_digits[i]
            ) *
            P[i];
    }

    /*
     * At the first violating digit,
     * x_h can take 0,...,m_h.
     * The tight path terminates here.
     */
    misses +=
        static_cast<u128>(
            m_digits[h] + 1
        ) *
        P[h];

    return misses;
}

u128 count_hits_closed_form(
    u128 n,
    u128 m,
    u64 p
) {
    if (n > m) {
        n = m;
    }

    u128 misses =
        count_misses_closed_form(
            n,
            m,
            p
        );

    return n + 1 - misses;
}

/*
 * --------------------------------------------------------------------------
 * INDEPENDENT RECURSIVE REFERENCE
 * --------------------------------------------------------------------------
 *
 * F(n,m) recursively counts digitwise-dominated x <= n.
 *
 * For the most significant digit:
 *
 *   n_h < m_h:
 *       n_h * lower_product + F(n_lower,m_lower)
 *
 *   n_h = m_h:
 *       n_h * lower_product + F(n_lower,m_lower)
 *
 *   n_h > m_h:
 *       (m_h+1) * lower_product
 *
 * Base case:
 *
 *   F(0,0)=1.
 */
u128 count_misses_recursive_digits(
    const std::vector<u64> &n_digits,
    const std::vector<u64> &m_digits,
    std::size_t pos
) {
    if (pos == 0) {
        return 1;
    }

    std::size_t i =
        pos - 1;

    u64 n_i =
        i < n_digits.size()
            ? n_digits[i]
            : 0;

    u64 m_i =
        i < m_digits.size()
            ? m_digits[i]
            : 0;

    /*
     * Product for all lower positions.
     */
    u128 lower_product = 1;

    for (std::size_t j = 0;
         j < i;
         ++j) {

        u64 md =
            j < m_digits.size()
                ? m_digits[j]
                : 0;

        lower_product *=
            static_cast<u128>(md + 1);
    }

    if (n_i > m_i) {
        return
            static_cast<u128>(m_i + 1) *
            lower_product;
    }

    u128 prefix =
        static_cast<u128>(n_i) *
        lower_product;

    return
        prefix +
        count_misses_recursive_digits(
            n_digits,
            m_digits,
            i
        );
}

u128 count_misses_recursive(
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

    return
        count_misses_recursive_digits(
            nd,
            md,
            L
        );
}

/*
 * Independent exhaustive reference for small domains.
 */
u128 brute_miss_count(
    u64 p,
    u128 m,
    u64 n
) {
    u128 count = 0;

    for (u64 x = 0;
         x <= n;
         ++x) {

        std::vector<u64> xd =
            digits_base(
                static_cast<u128>(x),
                p
            );

        std::vector<u64> md =
            digits_base(
                m,
                p
            );

        std::size_t L =
            std::max(
                xd.size(),
                md.size()
            );

        bool admissible = true;

        for (std::size_t i = 0;
             i < L;
             ++i) {

            u64 a =
                i < xd.size()
                    ? xd[i]
                    : 0;

            u64 b =
                i < md.size()
                    ? md[i]
                    : 0;

            if (a > b) {
                admissible = false;
                break;
            }
        }

        if (admissible) {
            ++count;
        }
    }

    return count;
}

void print_case(const Case &c) {
    BuiltCase bc =
        build_case(c);

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

    std::cout << '\n';
}

/*
 * --------------------------------------------------------------------------
 * 1. Deterministic tests
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
        {5, 3, 1, 4, 1000007654321ULL},
        {2, 5, 0, 2, 1073741825ULL},
        {7, 2, 3, 1, 123456789ULL}
    };

    u64 failures = 0;

    std::cout
        << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        BuiltCase bc =
            build_case(c);

        bool pass = true;

        /*
         * Always test n=0, n=m, and the structural
         * points around s0.
         */
        std::vector<u128> points;

        points.push_back(0);
        points.push_back(bc.m);

        if (bc.s0 > 0) {
            points.push_back(bc.s0 - 1);
        }

        points.push_back(bc.s0);

        if (bc.s0 < bc.m) {
            points.push_back(bc.s0 + 1);
        }

        if (bc.m > 0) {
            points.push_back(bc.m - 1);
        }

        /*
         * Small domains are exhaustively checked.
         */
        if (bc.m <=
            static_cast<u128>(10000)) {

            for (u64 n = 0;
                 n <= static_cast<u64>(bc.m);
                 ++n) {

                points.push_back(
                    static_cast<u128>(n)
                );
            }
        }

        for (u128 n : points) {
            if (n > bc.m) {
                continue;
            }

            u128 closed =
                count_misses_closed_form(
                    n,
                    bc.m,
                    c.p
                );

            u128 recursive =
                count_misses_recursive(
                    n,
                    bc.m,
                    c.p
                );

            if (closed != recursive) {
                pass = false;
                break;
            }

            if (bc.m <=
                static_cast<u128>(10000)) {

                u128 brute =
                    brute_miss_count(
                        c.p,
                        bc.m,
                        static_cast<u64>(n)
                    );

                if (closed != brute) {
                    pass = false;
                    break;
                }
            }
        }

        /*
         * Total product identity at n=m:
         *
         * M(m)=prod_i(m_i+1).
         */
        std::vector<u64> md =
            digits_base(
                bc.m,
                c.p
            );

        u128 product = 1;

        for (u64 digit : md) {
            product *=
                static_cast<u128>(
                    digit + 1
                );
        }

        u128 total_at_m =
            count_misses_closed_form(
                bc.m,
                bc.m,
                c.p
            );

        if (total_at_m != product) {
            pass = false;
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
 * 2. Large random closed-vs-recursive test
 * --------------------------------------------------------------------------
 */
bool random_tests() {
    std::mt19937_64 rng(
        0x250123456789ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 choices[] = {
            2, 3, 5, 7, 11, 13
        };

        u64 p =
            choices[rng() % 6];

        u64 a0 =
            rng() % 10;

        u64 b =
            rng() % (p - 1);

        u64 z =
            rng() % 10;

        u64 q;

        switch (i % 5) {
            case 0:
                q =
                    1 +
                    rng() % 1000;
                break;

            case 1:
                q =
                    1 +
                    rng() % 1000000ULL;
                break;

            case 2:
                q =
                    1 +
                    rng() % 1000000000ULL;
                break;

            case 3:
                q =
                    1 +
                    rng() % 1000000000000ULL;
                break;

            default:
                q =
                    1 +
                    rng() %
                    1000000000000000000ULL;
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

        /*
         * Keep n within u64 for convenient random generation.
         */
        if (bc.m >
            static_cast<u128>(
                UINT64_MAX
            )) {
            continue;
        }

        u64 m =
            static_cast<u64>(bc.m);

        u64 n =
            rng() % (m + 1);

        u128 closed =
            count_misses_closed_form(
                static_cast<u128>(n),
                bc.m,
                c.p
            );

        u128 recursive =
            count_misses_recursive(
                static_cast<u128>(n),
                bc.m,
                c.p
            );

        if (closed != recursive) {
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
                    << "closed=";

                print_u128(closed);

                std::cout
                    << " recursive=";

                print_u128(recursive);

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
 * 3. Digit-prefix identity
 * --------------------------------------------------------------------------
 *
 * Explicitly verify:
 *
 * if h is the most significant digit with n_h > m_h,
 *
 * M(n)
 * =
 * sum_{i>h} n_i P_i
 * +
 * (m_h+1)P_h.
 *
 * Otherwise:
 *
 * M(n)
 * =
 * 1 + sum_i n_i P_i.
 */
bool prefix_identity_tests() {
    std::mt19937_64 rng(
        0x250777888ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 cidx = 0;
         cidx < cases;
         ++cidx) {

        const u64 choices[] = {
            2, 3, 5, 7, 11
        };

        u64 p =
            choices[rng() % 5];

        u128 m =
            static_cast<u128>(
                rng()
            );

        u64 extra =
            rng() % 8;

        m +=
            pow_u128(
                p,
                18 + extra
            );

        u128 n =
            (static_cast<u128>(rng()) << 32) ^
            static_cast<u128>(rng());

        if (n > m) {
            n = m;
        }

        std::vector<u64> nd =
            digits_base(n, p);

        std::vector<u64> md =
            digits_base(m, p);

        const std::size_t L =
            std::max(
                nd.size(),
                md.size()
            );

        std::vector<u64> n_digits(L, 0);
        std::vector<u64> m_digits(L, 0);

        for (std::size_t i = 0;
             i < nd.size();
             ++i) {

            n_digits[i] = nd[i];
        }

        for (std::size_t i = 0;
             i < md.size();
             ++i) {

            m_digits[i] = md[i];
        }

        std::vector<u128> P(
            L + 1,
            1
        );

        for (std::size_t i = 0;
             i < L;
             ++i) {

            P[i + 1] =
                P[i] *
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

        u128 explicit_formula = 0;

        if (h == L) {
            explicit_formula = 1;

            for (std::size_t i = 0;
                 i < L;
                 ++i) {

                explicit_formula +=
                    static_cast<u128>(
                        n_digits[i]
                    ) *
                    P[i];
            }

        } else {
            for (std::size_t i = h + 1;
                 i < L;
                 ++i) {

                explicit_formula +=
                    static_cast<u128>(
                        n_digits[i]
                    ) *
                    P[i];
            }

            explicit_formula +=
                static_cast<u128>(
                    m_digits[h] + 1
                ) *
                P[h];
        }

        u128 closed =
            count_misses_closed_form(
                n,
                m,
                p
            );

        if (explicit_formula != closed) {
            ++failures;

            if (failures <= 5) {
                std::cout
                    << "PREFIX FAILURE\n";

                std::cout
                    << "p="
                    << p
                    << " n=";

                print_u128(n);

                std::cout
                    << " m=";

                print_u128(m);

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "prefix_cases="
        << cases
        << " prefix_failures="
        << failures
        << " prefix_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * --------------------------------------------------------------------------
 * 4. Total-product bridge
 * --------------------------------------------------------------------------
 *
 * For n=m:
 *
 * M(m)=prod_i(m_i+1)
 *
 * and therefore
 *
 * C(m)=m+1-prod_i(m_i+1).
 */
bool total_product_tests() {
    std::mt19937_64 rng(
        0x250424242ULL
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

        u128 m =
            static_cast<u128>(rng());

        std::vector<u64> md =
            digits_base(m, p);

        u128 product = 1;

        for (u64 digit : md) {
            product *=
                static_cast<u128>(
                    digit + 1
                );
        }

        u128 misses =
            count_misses_closed_form(
                m,
                m,
                p
            );

        u128 hits =
            count_hits_closed_form(
                m,
                m,
                p
            );

        if (misses != product ||
            hits !=
                m + 1 - product) {

            ++failures;

            if (failures <= 5) {
                std::cout
                    << "PRODUCT FAILURE\n";

                std::cout
                    << "p="
                    << p
                    << " m=";

                print_u128(m);

                std::cout
                    << '\n';
            }
        }
    }

    std::cout
        << "product_cases="
        << cases
        << " product_failures="
        << failures
        << " product_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 250\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_ok =
        random_tests();

    bool prefix_ok =
        prefix_identity_tests();

    bool product_ok =
        total_product_tests();

    bool overall =
        deterministic_ok &&
        random_ok &&
        prefix_ok &&
        product_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 250\n";

    return overall ? 0 : 1;
}
