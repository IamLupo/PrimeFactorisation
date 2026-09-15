#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

struct Case {
    u64 p;
    u64 q;
};

static const std::array<int, 6> BASES = {2, 3, 5, 7, 11, 13};
static const std::array<int, 5> T_VALUES = {2, 3, 5, 7, 11};

static u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

static u64 isqrt_u64(u64 n) {
    u64 lo = 0;
    u64 hi = std::min<u64>(n, 1ULL << 32);

    while (lo <= hi) {
        const u64 mid = lo + (hi - lo) / 2;

        if (mid == 0) {
            lo = 1;
            continue;
        }

        if (mid <= n / mid) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return hi;
}

static u64 mul_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((u128(a) * u128(b)) % mod);
}

static u64 pow_mod(u64 a, u64 e, u64 mod) {
    u64 result = 1 % mod;

    while (e > 0) {
        if (e & 1ULL) {
            result = mul_mod(result, a, mod);
        }

        a = mul_mod(a, a, mod);
        e >>= 1ULL;
    }

    return result;
}

static u64 prefix_miss(
    u64 m,
    u64 y,
    int base
) {
    if (y >= m) {
        u64 result = 1;
        u64 x = m;

        if (x == 0) {
            return 1;
        }

        while (x > 0) {
            result *= (x % base) + 1;
            x /= base;
        }

        return result;
    }

    std::vector<u64> mdigits;
    std::vector<u64> ydigits;

    u64 mx = m;
    u64 yx = y;

    while (mx > 0 || yx > 0) {
        mdigits.push_back(mx % base);
        ydigits.push_back(yx % base);

        mx /= base;
        yx /= base;
    }

    const std::size_t L =
        std::max(mdigits.size(), ydigits.size());

    while (mdigits.size() < L) {
        mdigits.push_back(0);
    }

    while (ydigits.size() < L) {
        ydigits.push_back(0);
    }

    std::vector<u64> W(L + 1, 1);

    for (std::size_t i = 0; i < L; ++i) {
        W[i + 1] = W[i] * (mdigits[i] + 1);
    }

    int h = static_cast<int>(L) - 1;

    while (h >= 0 && ydigits[h] == mdigits[h]) {
        --h;
    }

    if (h < 0) {
        return W[L];
    }

    u64 result = 0;

    for (int i = static_cast<int>(L) - 1; i > h; --i) {
        result += ydigits[i] * W[i];
    }

    result += ydigits[h] * W[h];

    if (h == 0) {
        return result;
    }

    u64 low_m = 0;
    u64 low_y = 0;
    u64 power = 1;

    for (int i = 0; i < h; ++i) {
        low_m += mdigits[i] * power;
        low_y += ydigits[i] * power;
        power *= static_cast<u64>(base);
    }

    result += prefix_miss(low_m, low_y, base);

    return result;
}

static std::vector<i128> coefficient_formula(
    u64 m,
    u64 y,
    int base
) {
    std::vector<u64> digits;

    u64 x = m;

    while (x > 0) {
        digits.push_back(x % base);
        x /= base;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    const int L = static_cast<int>(digits.size());

    std::vector<u64> W(L + 1, 1);

    for (int i = 0; i < L; ++i) {
        W[i + 1] = W[i] * (digits[i] + 1);
    }

    std::vector<i128> coeff(L, 0);

    u64 power_r = 1;

    for (int r = 0; r < L; ++r) {
        const u64 Q = y / power_r;
        const u64 R = y % power_r;

        if (Q == 0) {
            coeff[r] = 0;
        } else if (Q > digits[r]) {
            coeff[r] =
                static_cast<i128>(digits[r]) *
                static_cast<i128>(W[r]);
        } else {
            const u64 low_m = m % power_r;
            const u64 A = prefix_miss(low_m, R, base);

            coeff[r] =
                static_cast<i128>(Q - 1) *
                static_cast<i128>(W[r]) +
                static_cast<i128>(A);
        }

        if (r + 1 < L) {
            power_r *= static_cast<u64>(base);
        }
    }

    return coeff;
}

static i128 evaluate_polynomial(
    const std::vector<i128>& coeff,
    u64 t
) {
    i128 result = 0;
    i128 power = 1;

    for (i128 c : coeff) {
        result += c * power;
        power *= static_cast<i128>(t);
    }

    return result;
}

static u64 signed_mod(i128 x, u64 mod) {
    const i128 M = static_cast<i128>(mod);

    i128 r = x % M;

    if (r < 0) {
        r += M;
    }

    return static_cast<u64>(r);
}

static u64 i128_abs_to_u64(i128 x) {
    if (x >= 0) {
        return static_cast<u64>(x);
    }

    return static_cast<u64>(-x);
}

static u64 lcm_safe(u64 a, u64 b) {
    const u64 g = gcd_u64(a, b);
    const u64 q = a / g;

    if (q > UINT64_MAX / b) {
        return UINT64_MAX;
    }

    return q * b;
}

static u64 lcm_1_to(u64 B) {
    u64 result = 1;

    for (u64 x = 2; x <= B; ++x) {
        result = lcm_safe(result, x);

        if (result == UINT64_MAX) {
            return 0;
        }
    }

    return result;
}

static int digit_count(
    u64 x,
    int base
) {
    int count = 1;

    while (x >= static_cast<u64>(base)) {
        x /= static_cast<u64>(base);
        ++count;
    }

    return count;
}

static std::vector<u64> build_factor_blind_exponents(
    u64 m,
    u64 D,
    int base,
    u64 t,
    const std::vector<i128>& coeff
) {
    std::set<u64> unique;

    const u64 len =
        static_cast<u64>(coeff.size());

    u64 support = 0;
    u64 sum_mod = 0;
    u64 weighted_mod = 0;

    for (std::size_t r = 0; r < coeff.size(); ++r) {
        if (coeff[r] != 0) {
            ++support;
        }

        const u64 abs_c =
            i128_abs_to_u64(coeff[r]);

        sum_mod =
            (sum_mod + abs_c % 1000003ULL)
            % 1000003ULL;

        const u64 weighted =
            static_cast<u64>(
                (u128(r) * u128(abs_c % 1000003ULL))
                % 1000003ULL
            );

        weighted_mod =
            (weighted_mod + weighted)
            % 1000003ULL;
    }

    const int digits =
        digit_count(m, base);

    /*
     * Direct small structural exponents.
     */
    unique.insert(2);
    unique.insert(3);
    unique.insert(4);
    unique.insert(5);
    unique.insert(6);
    unique.insert(7);
    unique.insert(8);

    /*
     * Quantities derived from the actual radix polynomial.
     */
    unique.insert(1 + support);
    unique.insert(1 + len);
    unique.insert(2 + support * 2);
    unique.insert(3 + digits * 3);
    unique.insert(5 + static_cast<u64>(base) + t);

    unique.insert(
        2 + (sum_mod % 31ULL)
    );

    unique.insert(
        2 + (weighted_mod % 31ULL)
    );

    unique.insert(
        2 + (D % 31ULL)
    );

    unique.insert(
        2 + (m % 31ULL)
    );

    /*
     * Factor-blind smoothness exponents.
     *
     * These are analogous to a p-1 type exponent,
     * but B is chosen solely from the representation.
     */
    const std::array<u64, 10> fixed_B = {
        8, 12, 16, 20, 24,
        28, 32, 36, 40, 42
    };

    for (u64 B : fixed_B) {
        const u64 k = lcm_1_to(B);

        if (k != 0 && k > 1) {
            unique.insert(k);
        }
    }

    /*
     * Representation-derived B values.
     */
    const std::array<u64, 6> derived_B = {
        8 + support,
        8 + len,
        8 + static_cast<u64>(digits),
        8 + (sum_mod % 20ULL),
        8 + (weighted_mod % 20ULL),
        8 + ((m + D) % 20ULL)
    };

    for (u64 B : derived_B) {
        if (B > 42) {
            B = 42;
        }

        const u64 k = lcm_1_to(B);

        if (k != 0 && k > 1) {
            unique.insert(k);
        }
    }

    return std::vector<u64>(
        unique.begin(),
        unique.end()
    );
}

static void classify_gcd(
    u64 g,
    u64 p,
    u64 q,
    u64& gcd1,
    u64& gcdN,
    u64& nontrivial,
    u64& p_only,
    u64& q_only,
    u64& other
) {
    if (g == 1) {
        ++gcd1;
        return;
    }

    const u64 N = p * q;

    if (g == N) {
        ++gcdN;
        return;
    }

    ++nontrivial;

    if (g == p) {
        ++p_only;
    } else if (g == q) {
        ++q_only;
    } else {
        ++other;
    }
}

int main() {
    std::cout << "START EXPERIMENT 315\n";
    std::cout << "FACTOR-BLIND POWER COLLISIONS\n";
    std::cout << "CAN RADIX-DERIVED EXPONENTS EXPOSE A HIDDEN FACTOR?\n\n";

    const std::array<Case, 30> cases = {{
        {81077, 162749},
        {125017, 174259},
        {57107, 88261},
        {52757, 148457},
        {19483, 34123},
        {123001, 188291},
        {97987, 112583},
        {76129, 192113},
        {41257, 73643},
        {137239, 140419},
        {64187, 115319},
        {162527, 184087},
        {66179, 124753},
        {87943, 158047},
        {42683, 52529},
        {107473, 149711},
        {75853, 82759},
        {150131, 175267},
        {99859, 124769},
        {75337, 162229},
        {124471, 168043},
        {51563, 87683},
        {52237, 142123},
        {19069, 28751},
        {117043, 187637},
        {92203, 112067},
        {75617, 185869},
        {40823, 67843},
        {134369, 136709},
        {63667, 109211}
    }};

    u64 total_values = 0;
    u64 total_exponents = 0;
    u64 total_tests = 0;

    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;
    u64 p_only = 0;
    u64 q_only = 0;
    u64 other = 0;

    u64 zero_x_mod_N = 0;

    u64 printed_hit = 0;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_values = 0;
        u64 case_exponents = 0;
        u64 case_tests = 0;
        u64 case_p_hits = 0;
        u64 case_q_hits = 0;
        u64 case_other_hits = 0;

        /*
         * Same D-boundary region as previous experiments.
         */
        for (int offset = 1; offset <= 5; ++offset) {
            if (s < static_cast<u64>(offset)) {
                continue;
            }

            const u64 m =
                s - static_cast<u64>(offset) + 1;

            const u64 D =
                N - m * m;

            if (D == 0) {
                continue;
            }

            for (int bi = 0; bi < 6; ++bi) {
                const int base = BASES[bi];

                const auto coeff =
                    coefficient_formula(m, N, base);

                for (int ti = 0; ti < 5; ++ti) {
                    const u64 t =
                        static_cast<u64>(T_VALUES[ti]);

                    const i128 raw_value =
                        evaluate_polynomial(coeff, t);

                    /*
                     * Test both:
                     *
                     *   X = P_b(t)
                     *   X = D + P_b(t)
                     *
                     * The exponent never uses p or q.
                     */
                    const std::array<i128, 2> values = {
                        raw_value,
                        raw_value + static_cast<i128>(D)
                    };

                    const auto exponents =
                        build_factor_blind_exponents(
                            m,
                            D,
                            base,
                            t,
                            coeff
                        );

                    for (i128 value : values) {
                        ++case_values;
                        ++total_values;

                        const u64 X =
                            signed_mod(value, N);

                        if (X == 0) {
                            ++zero_x_mod_N;
                            continue;
                        }

                        for (u64 k : exponents) {
                            ++case_exponents;
                            ++total_exponents;

                            /*
                             * X^k - 1 modulo N.
                             *
                             * Crucially k was constructed without
                             * knowing p or q.
                             */
                            const u64 residue =
                                pow_mod(X, k, N);

                            u64 delta;

                            if (residue >= 1) {
                                delta = residue - 1;
                            } else {
                                delta = N - (1 - residue);
                            }

                            const u64 g =
                                gcd_u64(delta, N);

                            ++case_tests;
                            ++total_tests;

                            const u64 old_p =
                                p_only;

                            const u64 old_q =
                                q_only;

                            const u64 old_other =
                                other;

                            classify_gcd(
                                g,
                                p,
                                q,
                                gcd1,
                                gcdN,
                                nontrivial,
                                p_only,
                                q_only,
                                other
                            );

                            if (p_only != old_p) {
                                ++case_p_hits;

                                if (printed_hit < 20) {
                                    std::cout
                                        << "P_HIT"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " k=" << k
                                        << " gcd=" << g
                                        << "\n";

                                    ++printed_hit;
                                }
                            }

                            if (q_only != old_q) {
                                ++case_q_hits;

                                if (printed_hit < 20) {
                                    std::cout
                                        << "Q_HIT"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " k=" << k
                                        << " gcd=" << g
                                        << "\n";

                                    ++printed_hit;
                                }
                            }

                            if (other != old_other) {
                                ++case_other_hits;

                                if (printed_hit < 20) {
                                    std::cout
                                        << "OTHER_HIT"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " k=" << k
                                        << " gcd=" << g
                                        << "\n";

                                    ++printed_hit;
                                }
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "CASE " << case_idx
            << " p=" << p
            << " q=" << q
            << "\n";

        std::cout
            << "values=" << case_values << "\n"
            << "exponents=" << case_exponents << "\n"
            << "tests=" << case_tests << "\n"
            << "p_only=" << case_p_hits << "\n"
            << "q_only=" << case_q_hits << "\n"
            << "other=" << case_other_hits << "\n"
            << "\n";
    }

    std::cout << "============================\n";
    std::cout << "TOTAL\n";

    std::cout
        << "values=" << total_values << "\n"
        << "exponents=" << total_exponents << "\n"
        << "tests=" << total_tests << "\n"
        << "zero_x_mod_N=" << zero_x_mod_N << "\n"
        << "\n";

    std::cout
        << "gcd1=" << gcd1 << "\n"
        << "gcdN=" << gcdN << "\n"
        << "nontrivial=" << nontrivial << "\n"
        << "p_only=" << p_only << "\n"
        << "q_only=" << q_only << "\n"
        << "other=" << other << "\n";

    std::cout << "\nFINISHED EXPERIMENT 315\n";

    return 0;
}
