#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

struct Data {
    u64 p;
    u128 m;

    std::vector<u64> m_digits;
    std::vector<u128> p_powers;
    std::vector<u128> weights;

    u128 miss_count;
};

struct StepResult {
    bool valid;
    u128 value;
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
    u64 p,
    u128 m
) {
    Data data;

    data.p = p;
    data.m = m;

    data.m_digits =
        digits_base(m, p);

    const std::size_t L =
        data.m_digits.size();

    data.p_powers.resize(L);
    data.weights.resize(L);

    u128 weight = 1;
    u128 count = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        data.p_powers[i] =
            pow_u128(
                p,
                static_cast<u64>(i)
            );

        data.weights[i] =
            weight;

        weight *=
            static_cast<u128>(
                data.m_digits[i] + 1
            );

        count *=
            static_cast<u128>(
                data.m_digits[i] + 1
            );
    }

    data.miss_count = count;

    return data;
}

bool is_miss(
    u128 n,
    const Data &data
) {
    if (n > data.m) {
        return false;
    }

    std::vector<u64> nd =
        digits_base(n, data.p);

    const std::size_t L =
        std::max(
            nd.size(),
            data.m_digits.size()
        );

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 a =
            i < nd.size() ? nd[i] : 0;

        u64 b =
            i < data.m_digits.size()
                ? data.m_digits[i]
                : 0;

        if (a > b) {
            return false;
        }
    }

    return true;
}

u128 rank_miss(
    u128 value,
    const Data &data
) {
    u128 rank = 0;

    for (std::size_t i = 0;
         i < data.m_digits.size();
         ++i) {

        u64 digit =
            static_cast<u64>(
                value % data.p
            );

        value /= data.p;

        rank +=
            static_cast<u128>(digit) *
            data.weights[i];
    }

    return rank;
}

u128 unrank_miss(
    u128 rank,
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

        value +=
            digit * place;

        place *=
            static_cast<u128>(
                data.p
            );
    }

    return value;
}

/*
 * Direct successor of a MISS value.
 *
 * Increment its mixed-radix rank by one,
 * implemented directly on the digits.
 */
StepResult miss_successor(
    u128 value,
    const Data &data
) {
    if (!is_miss(value, data)) {
        return {false, 0};
    }

    u128 rank =
        rank_miss(
            value,
            data
        );

    if (rank + 1 >= data.miss_count) {
        return {false, 0};
    }

    return {
        true,
        unrank_miss(
            rank + 1,
            data
        )
    };
}

/*
 * Direct predecessor of a MISS value.
 */
StepResult miss_predecessor(
    u128 value,
    const Data &data
) {
    if (!is_miss(value, data)) {
        return {false, 0};
    }

    u128 rank =
        rank_miss(
            value,
            data
        );

    if (rank == 0) {
        return {false, 0};
    }

    return {
        true,
        unrank_miss(
            rank - 1,
            data
        )
    };
}

/*
 * Largest MISS <= n.
 *
 * Uses the digit clamp directly.
 */
u128 predecessor_at(
    u128 n,
    const Data &data
) {
    if (n > data.m) {
        n = data.m;
    }

    std::vector<u64> nd =
        digits_base(n, data.p);

    const std::size_t L =
        std::max(
            nd.size(),
            data.m_digits.size()
        );

    nd.resize(L, 0);

    std::vector<u64> md =
        data.m_digits;

    md.resize(L, 0);

    std::size_t h = L;

    for (std::size_t i = L;
         i-- > 0;) {

        if (nd[i] > md[i]) {
            h = i;
            break;
        }
    }

    if (h == L) {
        return n;
    }

    u128 value = 0;
    u128 place = 1;

    for (std::size_t i = 0;
         i < L;
         ++i) {

        u64 digit =
            i > h ? nd[i] : md[i];

        value +=
            static_cast<u128>(digit) *
            place;

        place *=
            static_cast<u128>(
                data.p
            );
    }

    return value;
}

/*
 * Small exhaustive order algebra.
 */
bool deterministic_tests() {
    struct Test {
        u64 p;
        u128 m;
    };

    const std::vector<Test> tests = {
        {2, 0},
        {2, 1},
        {2, 2},
        {2, 7},
        {2, 29},
        {3, 6},
        {3, 17},
        {5, 20},
        {5, 100},
        {7, 1234},
        {11, 987654321ULL}
    };

    u64 failures = 0;

    std::cout << "DETERMINISTIC CASES\n";

    for (const Test &t : tests) {
        Data data =
            build_data(
                t.p,
                t.m
            );

        std::cout
            << "p=" << t.p
            << " m=";

        print_u128(t.m);

        std::cout
            << " miss_count=";

        print_u128(
            data.miss_count
        );

        std::cout << '\n';

        bool pass = true;

        if (data.miss_count <=
            static_cast<u128>(100000)) {

            u64 count =
                static_cast<u64>(
                    data.miss_count
                );

            for (u64 k = 0;
                 k < count;
                 ++k) {

                u128 x =
                    unrank_miss(
                        static_cast<u128>(k),
                        data
                    );

                if (!is_miss(x, data)) {
                    pass = false;
                    break;
                }

                if (rank_miss(x, data) !=
                    static_cast<u128>(k)) {
                    pass = false;
                    break;
                }

                StepResult s =
                    miss_successor(
                        x,
                        data
                    );

                if (k + 1 < count) {
                    if (!s.valid ||
                        s.value !=
                            unrank_miss(
                                static_cast<u128>(k + 1),
                                data
                            )) {

                        pass = false;
                        break;
                    }
                } else {
                    if (s.valid) {
                        pass = false;
                        break;
                    }
                }

                StepResult p =
                    miss_predecessor(
                        x,
                        data
                    );

                if (k > 0) {
                    if (!p.valid ||
                        p.value !=
                            unrank_miss(
                                static_cast<u128>(k - 1),
                                data
                            )) {

                        pass = false;
                        break;
                    }
                } else {
                    if (p.valid) {
                        pass = false;
                        break;
                    }
                }

                /*
                 * Inverse laws.
                 */
                if (s.valid) {
                    StepResult back =
                        miss_predecessor(
                            s.value,
                            data
                        );

                    if (!back.valid ||
                        back.value != x) {

                        pass = false;
                        break;
                    }
                }

                if (p.valid) {
                    StepResult forward =
                        miss_successor(
                            p.value,
                            data
                        );

                    if (!forward.valid ||
                        forward.value != x) {

                        pass = false;
                        break;
                    }
                }
            }
        } else {
            /*
             * Large case:
             * test selected ranks.
             */
            std::vector<u128> ranks = {
                0,
                1,
                data.miss_count / 2,
                data.miss_count - 2,
                data.miss_count - 1
            };

            for (u128 k : ranks) {
                if (k >= data.miss_count) {
                    continue;
                }

                u128 x =
                    unrank_miss(
                        k,
                        data
                    );

                if (!is_miss(x, data) ||
                    rank_miss(x, data) != k) {

                    pass = false;
                    break;
                }

                StepResult s =
                    miss_successor(
                        x,
                        data
                    );

                if (k + 1 <
                    data.miss_count) {

                    if (!s.valid ||
                        rank_miss(
                            s.value,
                            data
                        ) != k + 1) {

                        pass = false;
                        break;
                    }

                } else if (s.valid) {
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
        << tests.size()
        << " deterministic_failures="
        << failures
        << " deterministic_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * Random successor/predecessor inverse tests.
 */
bool random_order_tests() {
    std::mt19937_64 rng(
        0x256123456789ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 bases[] = {
            2, 3, 5, 7, 11, 13
        };

        u64 p =
            bases[
                rng() % 6
            ];

        u128 m =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        Data data =
            build_data(
                p,
                m
            );

        if (data.miss_count <= 1) {
            continue;
        }

        u128 rank;

        if (data.miss_count <=
            static_cast<u128>(
                UINT64_MAX
            )) {

            rank =
                static_cast<u128>(
                    rng() %
                    static_cast<u64>(
                        data.miss_count
                    )
                );

        } else {
            u128 r =
                (static_cast<u128>(rng()) << 64) |
                static_cast<u128>(rng());

            rank =
                r %
                data.miss_count;
        }

        u128 x =
            unrank_miss(
                rank,
                data
            );

        if (!is_miss(x, data) ||
            rank_miss(x, data) != rank) {

            ++failures;
            continue;
        }

        StepResult s =
            miss_successor(
                x,
                data
            );

        StepResult p_result =
            miss_predecessor(
                x,
                data
            );

        if (rank + 1 <
            data.miss_count) {

            if (!s.valid ||
                rank_miss(
                    s.value,
                    data
                ) != rank + 1) {

                ++failures;
                continue;
            }

            StepResult back =
                miss_predecessor(
                    s.value,
                    data
                );

            if (!back.valid ||
                back.value != x) {

                ++failures;
                continue;
            }
        } else {
            if (s.valid) {
                ++failures;
                continue;
            }
        }

        if (rank > 0) {
            if (!p_result.valid ||
                rank_miss(
                    p_result.value,
                    data
                ) != rank - 1) {

                ++failures;
                continue;
            }

            StepResult forward =
                miss_successor(
                    p_result.value,
                    data
                );

            if (!forward.valid ||
                forward.value != x) {

                ++failures;
                continue;
            }
        } else {
            if (p_result.valid) {
                ++failures;
                continue;
            }
        }
    }

    std::cout
        << "random_order_cases="
        << cases
        << " random_order_failures="
        << failures
        << " random_order_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * For arbitrary n:
 *
 * x = pred_M(n)
 *
 * If n is MISS:
 *     x=n
 *
 * If n is HIT:
 *     S(x)>n
 *
 * and the HIT component containing n is
 *
 *     [x+1,S(x)-1].
 */
bool arbitrary_n_tests() {
    std::mt19937_64 rng(
        0x256777888ULL
    );

    const u64 cases = 100000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 bases[] = {
            2, 3, 5, 7, 11
        };

        u64 p =
            bases[
                rng() % 5
            ];

        u128 m =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        Data data =
            build_data(
                p,
                m
            );

        u128 n =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        n %=
            data.m + 1;

        u128 pred =
            predecessor_at(
                n,
                data
            );

        if (pred > n ||
            !is_miss(pred, data)) {

            ++failures;
            continue;
        }

        bool n_miss =
            is_miss(
                n,
                data
            );

        if (n_miss) {
            if (pred != n) {
                ++failures;
                continue;
            }
        } else {
            if (pred >= n) {
                ++failures;
                continue;
            }

            StepResult successor =
                miss_successor(
                    pred,
                    data
                );

            if (!successor.valid) {
                ++failures;
                continue;
            }

            if (successor.value <= n) {
                ++failures;
                continue;
            }

            if (!is_miss(
                    successor.value,
                    data
                )) {

                ++failures;
                continue;
            }

            /*
             * The entire HIT component around n
             * is bounded by these two MISS values.
             */
            if (is_miss(
                    pred + 1,
                    data
                )) {

                ++failures;
                continue;
            }

            if (successor.value > pred + 1) {
                if (is_miss(
                        successor.value - 1,
                        data
                    )) {

                    ++failures;
                    continue;
                }
            }
        }

        /*
         * Rank/prefix identity.
         */
        u128 prefix_rank =
            rank_miss(
                pred,
                data
            );

        if (prefix_rank >=
            data.miss_count) {

            ++failures;
            continue;
        }
    }

    std::cout
        << "arbitrary_n_cases="
        << cases
        << " arbitrary_n_failures="
        << failures
        << " arbitrary_n_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

/*
 * Verify the rank-step equals the actual numerical
 * successor difference.
 *
 * This is a direct ordered-set statement:
 *
 *     successor(rank+1) > successor(rank)
 *
 * and no MISS lies between them.
 */
bool successor_gap_tests() {
    std::mt19937_64 rng(
        0x256424242ULL
    );

    const u64 cases = 50000;
    u64 failures = 0;

    for (u64 i = 0;
         i < cases;
         ++i) {

        const u64 bases[] = {
            2, 3, 5, 7
        };

        u64 p =
            bases[
                rng() % 4
            ];

        u128 m =
            (static_cast<u128>(rng()) << 64) |
            static_cast<u128>(rng());

        Data data =
            build_data(
                p,
                m
            );

        if (data.miss_count <= 1) {
            continue;
        }

        u128 k =
            static_cast<u128>(
                rng()
            ) %
            (data.miss_count - 1);

        u128 left =
            unrank_miss(
                k,
                data
            );

        u128 right =
            unrank_miss(
                k + 1,
                data
            );

        if (right <= left) {
            ++failures;
            continue;
        }

        /*
         * There can be no MISS between them because
         * they are consecutive mixed-radix ranks.
         *
         * Verify the endpoints and immediate interior.
         */
        if (!is_miss(left, data) ||
            !is_miss(right, data)) {

            ++failures;
            continue;
        }

        if (right > left + 1) {
            if (is_miss(
                    left + 1,
                    data
                ) ||
                is_miss(
                    right - 1,
                    data
                )) {

                ++failures;
                continue;
            }
        }

        /*
         * Reconstruct the same pair through
         * successor/predecessor.
         */
        StepResult s =
            miss_successor(
                left,
                data
            );

        StepResult p_result =
            miss_predecessor(
                right,
                data
            );

        if (!s.valid ||
            !p_result.valid ||
            s.value != right ||
            p_result.value != left) {

            ++failures;
            continue;
        }
    }

    std::cout
        << "successor_gap_cases="
        << cases
        << " successor_gap_failures="
        << failures
        << " successor_gap_pass="
        << (failures == 0 ? 1 : 0)
        << '\n';

    return failures == 0;
}

int main() {
    std::cout
        << "START EXPERIMENT 256\n";

    bool deterministic_ok =
        deterministic_tests();

    bool random_order_ok =
        random_order_tests();

    bool arbitrary_n_ok =
        arbitrary_n_tests();

    bool successor_gap_ok =
        successor_gap_tests();

    bool overall =
        deterministic_ok &&
        random_order_ok &&
        arbitrary_n_ok &&
        successor_gap_ok;

    std::cout
        << "OVERALL_PASS="
        << (overall ? 1 : 0)
        << '\n';

    std::cout
        << "FINISHED EXPERIMENT 256\n";

    return overall ? 0 : 1;
}
