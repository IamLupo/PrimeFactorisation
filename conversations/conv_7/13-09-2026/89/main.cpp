#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
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

struct RankData {
    std::vector<u64> digits;
    std::vector<u128> weights;
};

struct Interval {
    u128 j;
    u128 r;
    u128 x;
    u128 end;
    u128 next_j;
};

bool u128_to_u64(u128 x, u64 &out) {
    if (x > static_cast<u128>(UINT64_MAX)) {
        return false;
    }
    out = static_cast<u64>(x);
    return true;
}

u128 pow_u128(u64 base, u64 exp) {
    u128 result = 1;
    u128 b = base;

    while (exp > 0) {
        if (exp & 1ULL) {
            result *= b;
        }
        exp >>= 1ULL;

        if (exp > 0) {
            b *= b;
        }
    }

    return result;
}

std::vector<u64> digits_base(u64 value, u64 p) {
    std::vector<u64> digits;

    if (value == 0) {
        digits.push_back(0);
        return digits;
    }

    while (value > 0) {
        digits.push_back(value % p);
        value /= p;
    }

    return digits;
}

RankData build_rank_data(u64 p, u64 q) {
    RankData data;
    data.digits = digits_base(q, p);

    data.weights.resize(data.digits.size());
    u128 weight = 1;

    for (std::size_t i = 0; i < data.digits.size(); ++i) {
        data.weights[i] = weight;
        weight *= static_cast<u128>(data.digits[i] + 1);
    }

    return data;
}

u128 interval_count(const RankData &data) {
    u128 product = 1;

    for (u64 digit : data.digits) {
        product *= static_cast<u128>(digit + 1);
    }

    return product - 1;
}

u128 rank_j(u128 j, u64 p, const RankData &data) {
    u128 rank = 0;

    for (std::size_t i = 0; i < data.digits.size(); ++i) {
        u64 digit = static_cast<u64>(j % p);
        rank += static_cast<u128>(digit) * data.weights[i];
        j /= p;
    }

    return rank;
}

u128 unrank_j(u128 k, u64 p, const RankData &data) {
    u128 j = 0;
    u128 place = 1;

    for (std::size_t i = 0; i < data.digits.size(); ++i) {
        u128 radix = static_cast<u128>(data.digits[i] + 1);
        u128 digit = k % radix;
        k /= radix;

        j += digit * place;
        place *= static_cast<u128>(p);
    }

    return j;
}

u64 first_deficient_digit(u128 j, u64 p, const RankData &data) {
    u64 q_digit;

    for (std::size_t i = 0; i < data.digits.size(); ++i) {
        u64 j_digit = static_cast<u64>(j % p);
        j /= p;

        q_digit = data.digits[i];

        if (j_digit < q_digit) {
            return static_cast<u64>(i);
        }
    }

    return static_cast<u64>(data.digits.size());
}

u128 successor_j(u128 j, u64 p, const RankData &data) {
    u64 r = first_deficient_digit(j, p, data);

    u128 divisor = pow_u128(p, r);
    u128 higher = j / divisor;

    return (higher + 1) * divisor;
}

bool build_case(
    const Case &c,
    u128 &e,
    u128 &s0,
    u128 &m,
    u128 &p_pow_e
) {
    e = static_cast<u128>(c.a0) +
        static_cast<u128>(c.z) +
        1;

    u128 p_to_a = pow_u128(c.p, c.a0);

    s0 = static_cast<u128>(c.b + 1) * p_to_a;
    p_pow_e = pow_u128(c.p, static_cast<u64>(e));

    m = s0 + static_cast<u128>(c.q) * p_pow_e - 1;

    return true;
}

bool verify_interval(
    const Case &c,
    const RankData &data,
    u128 k,
    u128 I,
    u128 e,
    u128 s0,
    u128 m,
    u128 p_pow_e,
    bool verbose = false
) {
    u128 j = unrank_j(k, c.p, data);

    if (rank_j(j, c.p, data) != k) {
        return false;
    }

    u64 r = first_deficient_digit(j, c.p, data);

    if (r >= data.digits.size()) {
        return false;
    }

    u128 next_j = successor_j(j, c.p, data);

    if (next_j <= j) {
        return false;
    }

    u128 x = s0 + j * p_pow_e;

    u128 end = next_j * p_pow_e - 1;

    if (end < x) {
        return false;
    }

    if (k + 1 < I) {
        u128 next_x = s0 + unrank_j(k + 1, c.p, data) * p_pow_e;

        if (next_x != end + s0 + 1) {
            return false;
        }
    }

    if (k == I - 1) {
        if (end != m - s0) {
            return false;
        }
    }

    if (verbose) {
        u64 j64 = 0;
        u64 next64 = 0;

        u128_to_u64(j, j64);
        u128_to_u64(next_j, next64);

        std::cout
            << "  k=" << static_cast<u64>(k)
            << " j=" << j64
            << " r=" << r
            << " next_j=" << next64
            << '\n';
    }

    return true;
}

bool run_case(
    const Case &c,
    std::mt19937_64 &rng,
    u64 exhaustive_limit,
    u64 random_checks
) {
    u128 e;
    u128 s0;
    u128 m;
    u128 p_pow_e;

    build_case(c, e, s0, m, p_pow_e);

    RankData data = build_rank_data(c.p, c.q);
    u128 I = interval_count(data);

    if (I == 0) {
        return true;
    }

    u64 I64 = 0;

    if (u128_to_u64(I, I64) && I64 <= exhaustive_limit) {
        for (u64 k = 0; k < I64; ++k) {
            if (!verify_interval(
                    c, data, k, I, e, s0, m, p_pow_e)) {
                return false;
            }
        }
        return true;
    }

    for (u64 sample = 0; sample < random_checks; ++sample) {
        u128 k;

        if (I <= static_cast<u128>(UINT64_MAX)) {
            u64 limit = static_cast<u64>(I);
            k = static_cast<u128>(rng() % limit);
        } else {
            u128 hi = static_cast<u128>(rng());
            u128 lo = static_cast<u128>(rng());
            u128 value = (hi << 64) | lo;
            k = value % I;
        }

        if (!verify_interval(
                c, data, k, I, e, s0, m, p_pow_e)) {
            return false;
        }
    }

    return true;
}

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

void print_case(const Case &c) {
    u128 e;
    u128 s0;
    u128 m;
    u128 p_pow_e;

    build_case(c, e, s0, m, p_pow_e);

    RankData data = build_rank_data(c.p, c.q);
    u128 I = interval_count(data);

    std::cout
        << "p=" << c.p
        << " a0=" << c.a0
        << " b=" << c.b
        << " z=" << c.z
        << " q=" << c.q
        << " e=";
    print_u128(e);

    std::cout << " s0=";
    print_u128(s0);

    std::cout << " m=";
    print_u128(m);

    std::cout << " I=";
    print_u128(I);

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

    std::mt19937_64 rng(0x244244244ULL);

    u64 failures = 0;

    std::cout << "DETERMINISTIC CASES\n";

    for (const Case &c : cases) {
        print_case(c);

        bool pass = run_case(c, rng, 1000000, 100000);

        std::cout << "pass=" << (pass ? 1 : 0) << '\n';

        if (!pass) {
            ++failures;
        }
    }

    std::cout << "deterministic_cases=" << cases.size()
              << " deterministic_failures=" << failures
              << " deterministic_pass=" << (failures == 0 ? 1 : 0)
              << '\n';

    return failures == 0;
}

bool random_tests() {
    std::mt19937_64 rng(0x244123456789ULL);

    u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0; i < cases; ++i) {
        u64 p_choices[] = {2, 3, 5, 7, 11};

        u64 p = p_choices[rng() % 5];

        u64 a0 = rng() % 8;
        u64 b = rng() % (p - 1);
        u64 z = rng() % 8;

        u64 q;

        switch (i % 5) {
            case 0:
                q = rng() % 1000;
                break;

            case 1:
                q = 1000000ULL + rng() % 1000000ULL;
                break;

            case 2:
                q = 1000000000000ULL + rng() % 1000000000000ULL;
                break;

            case 3:
                q = 1000000000000000ULL + rng() % 1000000000000000ULL;
                break;

            default:
                q = rng();
                break;
        }

        Case c{p, a0, b, z, q};

        u128 e;
        u128 s0;
        u128 m;
        u128 p_pow_e;

        build_case(c, e, s0, m, p_pow_e);

        if (m > static_cast<u128>(UINT64_MAX)) {
            continue;
        }

        if (!run_case(c, rng, 50000, 1000)) {
            ++failures;

            if (failures <= 5) {
                std::cout << "RANDOM FAILURE\n";
                print_case(c);
            }
        }
    }

    std::cout
        << "random_cases=" << cases
        << " random_failures=" << failures
        << " random_pass=" << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

bool rank_successor_tests() {
    std::mt19937_64 rng(0x244555AAAULL);

    u64 cases = 10000;
    u64 failures = 0;

    for (u64 cidx = 0; cidx < cases; ++cidx) {
        u64 p_choices[] = {2, 3, 5, 7};
        u64 p = p_choices[rng() % 4];

        u64 q = 1 + rng() % 1000000000ULL;

        RankData data = build_rank_data(p, q);
        u128 I = interval_count(data);

        if (I == 0) {
            continue;
        }

        for (unsigned sample = 0; sample < 10; ++sample) {
            u128 k;

            if (I <= static_cast<u128>(UINT64_MAX)) {
                k = static_cast<u128>(rng() % static_cast<u64>(I));
            } else {
                u128 value =
                    (static_cast<u128>(rng()) << 64) |
                    static_cast<u128>(rng());

                k = value % I;
            }

            u128 j = unrank_j(k, p, data);
            u128 j_rank = rank_j(j, p, data);

            if (j_rank != k) {
                ++failures;
                continue;
            }

            u128 next_j = successor_j(j, p, data);

            u128 next_rank = rank_j(next_j, p, data);

            if (next_rank != k + 1) {
                ++failures;
            }
        }
    }

    std::cout
        << "successor_cases=" << cases
        << " successor_failures=" << failures
        << " successor_pass=" << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout << "START EXPERIMENT 244\n";

    bool deterministic_ok = deterministic_tests();
    bool random_ok = random_tests();
    bool successor_ok = rank_successor_tests();

    bool overall = deterministic_ok &&
                   random_ok &&
                   successor_ok;

    std::cout << "OVERALL_PASS=" << (overall ? 1 : 0) << '\n';
    std::cout << "FINISHED EXPERIMENT 244\n";

    return overall ? 0 : 1;
}
