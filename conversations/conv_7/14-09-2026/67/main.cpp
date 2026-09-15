#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

struct CaseData {
    u64 p;
    u64 q;
    u64 N;
    u64 s;
};

struct Polynomial {
    std::vector<i64> coeff;
};

static u64 isqrt_u64(u64 n) {
    u64 x = static_cast<u64>(
        __builtin_sqrtl(
            static_cast<long double>(n)
        )
    );

    while ((x + 1) <= n / (x + 1)) {
        ++x;
    }

    while (x > n / x) {
        --x;
    }

    return x;
}

static std::vector<int> sieve_primes(int limit) {
    std::vector<bool> composite(
        static_cast<std::size_t>(limit + 1),
        false
    );

    std::vector<int> primes;

    for (int i = 2; i <= limit; ++i) {
        if (composite[i]) {
            continue;
        }

        primes.push_back(i);

        if (static_cast<long long>(i) * i <= limit) {
            for (int j = i * i; j <= limit; j += i) {
                composite[j] = true;
            }
        }
    }

    return primes;
}

static std::vector<CaseData> generate_cases(
    const std::vector<int>& primes,
    int count,
    std::mt19937_64& rng,
    std::set<u64>& used_N
) {
    std::vector<int> usable;

    for (int p : primes) {
        if (p >= 1009 && p <= 50000) {
            usable.push_back(p);
        }
    }

    std::uniform_int_distribution<std::size_t> dist(
        0,
        usable.size() - 1
    );

    std::vector<CaseData> result;
    result.reserve(
        static_cast<std::size_t>(count)
    );

    while (static_cast<int>(result.size()) < count) {
        const std::size_t i = dist(rng);
        const std::size_t j = dist(rng);

        if (i == j) {
            continue;
        }

        const u64 p =
            static_cast<u64>(
                std::min(usable[i], usable[j])
            );

        const u64 q =
            static_cast<u64>(
                std::max(usable[i], usable[j])
            );

        const u64 N = p * q;

        if (!used_N.insert(N).second) {
            continue;
        }

        result.push_back({
            p,
            q,
            N,
            isqrt_u64(N)
        });
    }

    return result;
}

static i128 eval_polynomial(
    const Polynomial& P,
    i128 x
) {
    i128 result = 0;

    for (auto it = P.coeff.rbegin();
         it != P.coeff.rend();
         ++it) {

        result =
            result * x +
            static_cast<i128>(*it);
    }

    return result;
}

static u64 abs_i128_mod(
    i128 value,
    u64 mod
) {
    if (value < 0) {
        value = -value;
    }

    return static_cast<u64>(
        static_cast<u128>(value) %
        static_cast<u128>(mod)
    );
}

static void print_i128(i128 value) {
    if (value == 0) {
        std::cout << "0";
        return;
    }

    if (value < 0) {
        std::cout << '-';
        value = -value;
    }

    std::string out;

    while (value > 0) {
        const unsigned digit =
            static_cast<unsigned>(value % 10);

        out.push_back(
            static_cast<char>('0' + digit)
        );

        value /= 10;
    }

    std::reverse(
        out.begin(),
        out.end()
    );

    std::cout << out;
}

static Polynomial make_random_R(
    int degree,
    std::mt19937_64& rng
) {
    std::uniform_int_distribution<int> dist(
        -5,
        5
    );

    Polynomial R;

    R.coeff.resize(
        static_cast<std::size_t>(degree + 1)
    );

    for (int i = 0; i <= degree; ++i) {
        R.coeff[
            static_cast<std::size_t>(i)
        ] = static_cast<i64>(
            dist(rng)
        );
    }

    /*
     * Avoid the completely zero polynomial.
     */
    bool all_zero = true;

    for (i64 c : R.coeff) {
        if (c != 0) {
            all_zero = false;
            break;
        }
    }

    if (all_zero) {
        R.coeff[0] = 1;
    }

    return R;
}

static Polynomial make_constructed_P(
    const Polynomial& R,
    i64 c,
    u64 N,
    u64 s
) {
    /*
     * Construct
     *
     * P(x) = cN + (x-(s+1))R(x).
     *
     * Write
     *
     * (x-(s+1))R(x)
     *
     * coefficient-wise.
     */

    const int degree =
        static_cast<int>(R.coeff.size()) - 1;

    Polynomial P;

    P.coeff.assign(
        static_cast<std::size_t>(degree + 2),
        0
    );

    for (int j = 0; j <= degree + 1; ++j) {
        i128 value = 0;

        /*
         * x * R(x)
         */
        if (j >= 1) {
            value += static_cast<i128>(
                R.coeff[
                    static_cast<std::size_t>(j - 1)
                ]
            );
        }

        /*
         * -(s+1) * R(x)
         */
        if (j <= degree) {
            value -=
                static_cast<i128>(s + 1) *
                static_cast<i128>(
                    R.coeff[
                        static_cast<std::size_t>(j)
                    ]
                );
        }

        P.coeff[
            static_cast<std::size_t>(j)
        ] = static_cast<i64>(value);
    }

    /*
     * + cN
     */
    const i128 constant =
        static_cast<i128>(c) *
        static_cast<i128>(N);

    P.coeff[0] +=
        static_cast<i64>(constant);

    return P;
}

static Polynomial make_original_H(
    u64 s,
    u64 D
) {
    /*
     * H(x) = 4x^2 - (s+4)x + 3D - 3s.
     */

    Polynomial H;

    H.coeff.resize(3);

    H.coeff[0] =
        static_cast<i64>(
            static_cast<i128>(3) *
            static_cast<i128>(D)
            -
            static_cast<i128>(3) *
            static_cast<i128>(s)
        );

    H.coeff[1] =
        -static_cast<i64>(s + 4);

    H.coeff[2] = 4;

    return H;
}

int main() {
    std::cout << "START EXPERIMENT 356\n";

    constexpr int CASE_COUNT = 500;
    constexpr int POLYNOMIALS_PER_CASE = 20;
    constexpr int PRIME_LIMIT = 50000;
    constexpr int MAX_R_DEGREE = 5;

    std::mt19937_64 rng(
        0x356356356ULL
    );

    const std::vector<int> primes =
        sieve_primes(PRIME_LIMIT);

    std::set<u64> used_N;

    const std::vector<CaseData> cases =
        generate_cases(
            primes,
            CASE_COUNT,
            rng,
            used_N
        );

    u64 total_constructed = 0;

    u64 constructed_p_failures = 0;
    u64 constructed_q_failures = 0;

    u64 constructed_p_nonzero = 0;
    u64 constructed_q_nonzero = 0;

    u64 gcd_p_exact = 0;
    u64 gcd_q_exact = 0;

    u64 gcd_one = 0;
    u64 gcd_other = 0;

    u64 original_H_p_failures = 0;
    u64 original_H_q_failures = 0;

    u64 degree_0 = 0;
    u64 degree_1 = 0;
    u64 degree_2 = 0;
    u64 degree_3 = 0;
    u64 degree_4 = 0;
    u64 degree_5 = 0;

    for (std::size_t case_id = 0;
         case_id < cases.size();
         ++case_id) {

        const CaseData& c = cases[case_id];

        const u64 p = c.p;
        const u64 q = c.q;
        const u64 N = c.N;
        const u64 s = c.s;

        const u64 D =
            N - s * s;

        /*
         * ---------------------------------------------------------
         * Original H
         * ---------------------------------------------------------
         */

        const Polynomial H =
            make_original_H(
                s,
                D
            );

        const u64 xp =
            s + 1 - p;

        /*
         * q can be larger than s+1, so this is signed.
         */
        const i128 xq =
            static_cast<i128>(s)
            +
            static_cast<i128>(1)
            -
            static_cast<i128>(q);

        const i128 Hp =
            eval_polynomial(
                H,
                static_cast<i128>(xp)
            );

        const i128 Hq =
            eval_polynomial(
                H,
                xq
            );

        if (abs_i128_mod(Hp, p) != 0) {
            ++original_H_p_failures;
        }

        if (abs_i128_mod(Hq, q) != 0) {
            ++original_H_q_failures;
        }

        /*
         * ---------------------------------------------------------
         * Construct arbitrary polynomial families.
         * ---------------------------------------------------------
         */

        std::uniform_int_distribution<int> degree_dist(
            0,
            MAX_R_DEGREE
        );

        std::uniform_int_distribution<int> c_dist(
            -5,
            5
        );

        for (int j = 0;
             j < POLYNOMIALS_PER_CASE;
             ++j) {

            const int degree =
                degree_dist(rng);

            switch (degree) {
                case 0: ++degree_0; break;
                case 1: ++degree_1; break;
                case 2: ++degree_2; break;
                case 3: ++degree_3; break;
                case 4: ++degree_4; break;
                case 5: ++degree_5; break;
            }

            const Polynomial R =
                make_random_R(
                    degree,
                    rng
                );

            const i64 c_constant =
                static_cast<i64>(
                    c_dist(rng)
                );

            const Polynomial P =
                make_constructed_P(
                    R,
                    c_constant,
                    N,
                    s
                );

            ++total_constructed;

            /*
             * The fundamental points:
             *
             * x_p = s+1-p
             * x_q = s+1-q
             *
             * and
             *
             * x_p-(s+1) = -p
             * x_q-(s+1) = -q.
             */

            const i128 Px_p =
                eval_polynomial(
                    P,
                    static_cast<i128>(xp)
                );

            const i128 Px_q =
                eval_polynomial(
                    P,
                    xq
                );

            if (abs_i128_mod(Px_p, p) != 0) {
                ++constructed_p_failures;
            }

            if (abs_i128_mod(Px_q, q) != 0) {
                ++constructed_q_failures;
            }

            /*
             * Is the constructed value actually nonzero?
             */
            if (Px_p != 0) {
                ++constructed_p_nonzero;
            }

            if (Px_q != 0) {
                ++constructed_q_nonzero;
            }

            /*
             * GCD with N.
             *
             * This tells us whether the generic construction
             * automatically isolates a factor or only guarantees
             * divisibility.
             */
            const u64 gp =
                std::gcd(
                    abs_i128_mod(Px_p, N),
                    N
                );

            const u64 gq =
                std::gcd(
                    abs_i128_mod(Px_q, N),
                    N
                );

            if (gp == p) {
                ++gcd_p_exact;
            } else if (gp == 1) {
                ++gcd_one;
            } else {
                ++gcd_other;
            }

            if (gq == q) {
                ++gcd_q_exact;
            } else if (gq == 1) {
                ++gcd_one;
            } else {
                ++gcd_other;
            }

            /*
             * Print a few examples.
             */
            if (case_id == 0 && j < 5) {
                std::cout
                    << "\nPOLYNOMIAL_SAMPLE="
                    << (j + 1)
                    << "\n";

                std::cout
                    << "R_DEGREE="
                    << degree
                    << "\n";

                std::cout
                    << "C="
                    << c_constant
                    << "\n";

                const i128 P_s1 =
                    eval_polynomial(
                        P,
                        static_cast<i128>(
                            s + 1
                        )
                    );

                std::cout
                    << "P(s+1)=";

                print_i128(P_s1);

                std::cout
                    << "\n";

                std::cout
                    << "P(s+1-p)=";

                print_i128(Px_p);

                std::cout
                    << "\n";

                std::cout
                    << "P(s+1-q)=";

                print_i128(Px_q);

                std::cout
                    << "\n";

                std::cout
                    << "GCD_P="
                    << gp
                    << "\n";

                std::cout
                    << "GCD_Q="
                    << gq
                    << "\n";
            }
        }
    }

    std::cout << "\n";

    std::cout
        << "CASES="
        << cases.size()
        << "\n";

    std::cout
        << "POLYNOMIALS_PER_CASE="
        << POLYNOMIALS_PER_CASE
        << "\n";

    std::cout
        << "TOTAL_CONSTRUCTED="
        << total_constructed
        << "\n";

    std::cout
        << "CONSTRUCTED_P_FAILURES="
        << constructed_p_failures
        << "\n";

    std::cout
        << "CONSTRUCTED_Q_FAILURES="
        << constructed_q_failures
        << "\n";

    std::cout
        << "CONSTRUCTED_P_NONZERO="
        << constructed_p_nonzero
        << "\n";

    std::cout
        << "CONSTRUCTED_Q_NONZERO="
        << constructed_q_nonzero
        << "\n";

    std::cout
        << "GCD_P_EXACT="
        << gcd_p_exact
        << "\n";

    std::cout
        << "GCD_Q_EXACT="
        << gcd_q_exact
        << "\n";

    std::cout
        << "GCD_ONE="
        << gcd_one
        << "\n";

    std::cout
        << "GCD_OTHER="
        << gcd_other
        << "\n";

    std::cout
        << "ORIGINAL_H_P_FAILURES="
        << original_H_p_failures
        << "\n";

    std::cout
        << "ORIGINAL_H_Q_FAILURES="
        << original_H_q_failures
        << "\n";

    std::cout
        << "R_DEGREE_0="
        << degree_0
        << "\n";

    std::cout
        << "R_DEGREE_1="
        << degree_1
        << "\n";

    std::cout
        << "R_DEGREE_2="
        << degree_2
        << "\n";

    std::cout
        << "R_DEGREE_3="
        << degree_3
        << "\n";

    std::cout
        << "R_DEGREE_4="
        << degree_4
        << "\n";

    std::cout
        << "R_DEGREE_5="
        << degree_5
        << "\n";

    if (
        constructed_p_failures == 0 &&
        constructed_q_failures == 0 &&
        original_H_p_failures == 0 &&
        original_H_q_failures == 0
    ) {
        std::cout
            << "GENERIC_DIVISIBILITY_STATUS=PASS\n";
    } else {
        std::cout
            << "GENERIC_DIVISIBILITY_STATUS=FAIL\n";
    }

    std::cout
        << "FINISHED EXPERIMENT 356\n";

    return 0;
}