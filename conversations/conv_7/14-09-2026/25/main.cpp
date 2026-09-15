#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;
using u128 = __uint128_t;

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

static std::vector<u64> distinct_prime_factors(u64 n) {
    std::vector<u64> factors;

    while ((n & 1ULL) == 0) {
        if (factors.empty() || factors.back() != 2) {
            factors.push_back(2);
        }
        n /= 2;
    }

    for (u64 d = 3; d <= n / d; d += 2) {
        if (n % d != 0) {
            continue;
        }

        factors.push_back(d);

        while (n % d == 0) {
            n /= d;
        }
    }

    if (n > 1) {
        factors.push_back(n);
    }

    return factors;
}

static u64 multiplicative_order_prime_modulus(
    u64 x,
    u64 prime
) {
    x %= prime;

    if (x == 0) {
        return 0;
    }

    const u64 phi = prime - 1;
    u64 order = phi;

    const auto factors = distinct_prime_factors(phi);

    for (u64 r : factors) {
        while (order % r == 0 &&
               pow_mod(x, order / r, prime) == 1) {
            order /= r;
        }
    }

    return order;
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

static u64 gcd_power_minus_one(
    u64 x,
    u64 exponent,
    u64 N
) {
    if (x == 0) {
        return 1;
    }

    const u64 residue = pow_mod(x, exponent, N);

    if (residue >= 1) {
        const u64 delta = residue - 1;
        return gcd_u64(delta, N);
    }

    return N;
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
    std::cout << "START EXPERIMENT 314\n";
    std::cout << "MULTIPLICATIVE ORDER FINGERPRINTS\n";
    std::cout << "CAN RADIX FINGERPRINTS HAVE DIFFERENT MULTIPLICATIVE DYNAMICS MODULO p AND q?\n\n";

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

    u64 total_expressions = 0;
    u64 total_values = 0;

    u64 zero_p = 0;
    u64 zero_q = 0;
    u64 nonzero_both = 0;

    u64 order_equal = 0;
    u64 order_diff = 0;

    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;
    u64 p_only = 0;
    u64 q_only = 0;
    u64 other = 0;

    u64 fermat_p_hits = 0;
    u64 fermat_q_hits = 0;

    u64 order_p_hits = 0;
    u64 order_q_hits = 0;

    u64 cross_order_p_hits = 0;
    u64 cross_order_q_hits = 0;

    u64 printed_order_difference = 0;
    u64 printed_factor_hit = 0;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_expressions = 0;
        u64 case_values = 0;
        u64 case_order_equal = 0;
        u64 case_order_diff = 0;
        u64 case_factor_hits = 0;

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

            ++case_expressions;
            ++total_expressions;

            for (int bi = 0; bi < 6; ++bi) {
                const int base = BASES[bi];

                const auto coeff =
                    coefficient_formula(m, N, base);

                for (int ti = 0; ti < 5; ++ti) {
                    const u64 t =
                        static_cast<u64>(T_VALUES[ti]);

                    const i128 raw_value =
                        evaluate_polynomial(coeff, t);

                    const std::array<i128, 2> values = {
                        raw_value,
                        raw_value + static_cast<i128>(D)
                    };

                    for (i128 value : values) {
                        ++case_values;
                        ++total_values;

                        const u64 xp =
                            signed_mod(value, p);

                        const u64 xq =
                            signed_mod(value, q);

                        if (xp == 0) {
                            ++zero_p;
                        }

                        if (xq == 0) {
                            ++zero_q;
                        }

                        if (xp == 0 || xq == 0) {
                            continue;
                        }

                        ++nonzero_both;

                        const u64 ord_p =
                            multiplicative_order_prime_modulus(
                                xp, p);

                        const u64 ord_q =
                            multiplicative_order_prime_modulus(
                                xq, q);

                        if (ord_p == ord_q) {
                            ++order_equal;
                            ++case_order_equal;
                        } else {
                            ++order_diff;
                            ++case_order_diff;

                            if (printed_order_difference < 12) {
                                std::cout
                                    << "ORDER_DIFF"
                                    << " case=" << case_idx
                                    << " base=" << base
                                    << " t=" << t
                                    << " ord_p=" << ord_p
                                    << " ord_q=" << ord_q
                                    << "\n";

                                ++printed_order_difference;
                            }
                        }

                        /*
                         * Direct order exponents.
                         */
                        {
                            const u64 g1 =
                                gcd_power_minus_one(
                                    signed_mod(value, N),
                                    ord_p,
                                    N);

                            if (g1 == p) {
                                ++order_p_hits;
                            }

                            if (g1 == q) {
                                ++order_q_hits;
                            }

                            classify_gcd(
                                g1, p, q,
                                gcd1, gcdN,
                                nontrivial,
                                p_only, q_only, other);

                            if (g1 != 1 && g1 != N) {
                                ++case_factor_hits;

                                if (printed_factor_hit < 12) {
                                    std::cout
                                        << "FACTOR_HIT"
                                        << " exponent=ord_p"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " gcd=" << g1
                                        << "\n";

                                    ++printed_factor_hit;
                                }
                            }
                        }

                        {
                            const u64 g2 =
                                gcd_power_minus_one(
                                    signed_mod(value, N),
                                    ord_q,
                                    N);

                            if (g2 == p) {
                                ++order_p_hits;
                            }

                            if (g2 == q) {
                                ++order_q_hits;
                            }

                            classify_gcd(
                                g2, p, q,
                                gcd1, gcdN,
                                nontrivial,
                                p_only, q_only, other);

                            if (g2 != 1 && g2 != N) {
                                ++case_factor_hits;

                                if (printed_factor_hit < 12) {
                                    std::cout
                                        << "FACTOR_HIT"
                                        << " exponent=ord_q"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " gcd=" << g2
                                        << "\n";

                                    ++printed_factor_hit;
                                }
                            }
                        }

                        /*
                         * Fermat exponents.
                         */
                        {
                            const u64 g3 =
                                gcd_power_minus_one(
                                    signed_mod(value, N),
                                    p - 1,
                                    N);

                            if (g3 == p) {
                                ++fermat_p_hits;
                            }

                            if (g3 == q) {
                                ++fermat_q_hits;
                            }

                            classify_gcd(
                                g3, p, q,
                                gcd1, gcdN,
                                nontrivial,
                                p_only, q_only, other);

                            if (g3 != 1 && g3 != N) {
                                ++case_factor_hits;

                                if (printed_factor_hit < 12) {
                                    std::cout
                                        << "FACTOR_HIT"
                                        << " exponent=p-1"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " gcd=" << g3
                                        << "\n";

                                    ++printed_factor_hit;
                                }
                            }
                        }

                        {
                            const u64 g4 =
                                gcd_power_minus_one(
                                    signed_mod(value, N),
                                    q - 1,
                                    N);

                            if (g4 == p) {
                                ++fermat_p_hits;
                            }

                            if (g4 == q) {
                                ++fermat_q_hits;
                            }

                            classify_gcd(
                                g4, p, q,
                                gcd1, gcdN,
                                nontrivial,
                                p_only, q_only, other);

                            if (g4 != 1 && g4 != N) {
                                ++case_factor_hits;

                                if (printed_factor_hit < 12) {
                                    std::cout
                                        << "FACTOR_HIT"
                                        << " exponent=q-1"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " gcd=" << g4
                                        << "\n";

                                    ++printed_factor_hit;
                                }
                            }
                        }

                        /*
                         * Cross-order exponents:
                         *
                         * gcd(X^ord_p - 1, N)
                         * gcd(X^ord_q - 1, N)
                         *
                         * are already above, but we explicitly record
                         * which hidden field the corresponding order
                         * annihilates.
                         */
                        const u64 xp_check =
                            pow_mod(
                                signed_mod(value, p),
                                ord_p,
                                p);

                        const u64 xq_check =
                            pow_mod(
                                signed_mod(value, q),
                                ord_q,
                                q);

                        if (xp_check == 1 &&
                            xq_check != 1) {
                            ++cross_order_p_hits;
                        }

                        if (xq_check == 1 &&
                            xp_check != 1) {
                            ++cross_order_q_hits;
                        }
                    }
                }
            }
        }

        std::cout << "CASE " << case_idx
                  << " p=" << p
                  << " q=" << q
                  << "\n";

        std::cout
            << "expressions=" << case_expressions << "\n"
            << "values=" << case_values << "\n"
            << "order_equal=" << case_order_equal << "\n"
            << "order_diff=" << case_order_diff << "\n"
            << "factor_hits=" << case_factor_hits << "\n"
            << "\n";
    }

    std::cout << "============================\n";
    std::cout << "TOTAL\n";

    std::cout
        << "expressions=" << total_expressions << "\n"
        << "values=" << total_values << "\n"
        << "zero_p=" << zero_p << "\n"
        << "zero_q=" << zero_q << "\n"
        << "nonzero_both=" << nonzero_both << "\n"
        << "\n";

    std::cout
        << "MULTIPLICATIVE ORDERS\n"
        << "order_equal=" << order_equal << "\n"
        << "order_diff=" << order_diff << "\n"
        << "\n";

    std::cout
        << "GCD RESULTS\n"
        << "gcd1=" << gcd1 << "\n"
        << "gcdN=" << gcdN << "\n"
        << "nontrivial=" << nontrivial << "\n"
        << "p_only=" << p_only << "\n"
        << "q_only=" << q_only << "\n"
        << "other=" << other << "\n"
        << "\n";

    std::cout
        << "DIRECT ORDER HITS\n"
        << "order_p_hits=" << order_p_hits << "\n"
        << "order_q_hits=" << order_q_hits << "\n"
        << "\n";

    std::cout
        << "FERMAT HITS\n"
        << "fermat_p_hits=" << fermat_p_hits << "\n"
        << "fermat_q_hits=" << fermat_q_hits << "\n"
        << "\n";

    std::cout
        << "CROSS-FIELD ORDER ASYMMETRY\n"
        << "cross_order_p_hits=" << cross_order_p_hits << "\n"
        << "cross_order_q_hits=" << cross_order_q_hits << "\n";

    std::cout << "\nFINISHED EXPERIMENT 314\n";

    return 0;
}
