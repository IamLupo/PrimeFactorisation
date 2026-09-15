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

static u64 signed_mod(i128 x, u64 mod) {
    const i128 M = static_cast<i128>(mod);

    i128 r = x % M;

    if (r < 0) {
        r += M;
    }

    return static_cast<u64>(r);
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

static u64 support_size(
    const std::vector<i128>& coeff
) {
    u64 count = 0;

    for (i128 c : coeff) {
        if (c != 0) {
            ++count;
        }
    }

    return count;
}

static u64 coefficient_sum_mod(
    const std::vector<i128>& coeff,
    u64 mod
) {
    u64 result = 0;

    for (i128 c : coeff) {
        if (c < 0) {
            c = -c;
        }

        result += signed_mod(c, mod);
        result %= mod;
    }

    return result;
}

static u64 weighted_sum_mod(
    const std::vector<i128>& coeff,
    u64 mod
) {
    u64 result = 0;

    for (std::size_t r = 0; r < coeff.size(); ++r) {
        i128 c = coeff[r];

        if (c < 0) {
            c = -c;
        }

        const u64 cm = signed_mod(c, mod);

        const u64 term =
            static_cast<u64>(
                (u128(r) * u128(cm)) % mod
            );

        result += term;
        result %= mod;
    }

    return result;
}

static std::vector<u64> intrinsic_exponents(
    const std::vector<i128>& coeff,
    u64 m,
    u64 D,
    int base,
    u64 t
) {
    std::set<u64> out;

    const u64 deg =
        coeff.empty()
            ? 0
            : static_cast<u64>(coeff.size() - 1);

    const u64 support =
        support_size(coeff);

    const u64 s97 =
        coefficient_sum_mod(coeff, 97);

    const u64 s193 =
        coefficient_sum_mod(coeff, 193);

    const u64 w97 =
        weighted_sum_mod(coeff, 97);

    const u64 w193 =
        weighted_sum_mod(coeff, 193);

    out.insert(2 + support);
    out.insert(2 + deg);
    out.insert(2 + static_cast<u64>(base));
    out.insert(2 + t);

    out.insert(2 + s97 % 31);
    out.insert(2 + s193 % 37);
    out.insert(2 + w97 % 31);
    out.insert(2 + w193 % 37);

    out.insert(2 + m % 41);
    out.insert(2 + D % 41);

    out.insert(
        2 + (support * static_cast<u64>(base)) % 43
    );

    out.insert(
        2 + (deg * t) % 43
    );

    out.insert(
        2 + (s97 + w97) % 43
    );

    out.insert(
        2 + (s193 + w193) % 43
    );

    return std::vector<u64>(out.begin(), out.end());
}

static u64 exact_order_from_factor(
    u64 x,
    u64 prime
) {
    if (x % prime == 0) {
        return 0;
    }

    u64 order = prime - 1;

    /*
     * This is only an audit step.
     * It is NOT used to construct the factoring candidate.
     */
    std::vector<u64> factors;

    u64 n = order;

    if ((n & 1ULL) == 0) {
        factors.push_back(2);

        while ((n & 1ULL) == 0) {
            n >>= 1ULL;
        }
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

    for (u64 r : factors) {
        while (order % r == 0 &&
               pow_mod(
                   x % prime,
                   order / r,
                   prime
               ) == 1) {
            order /= r;
        }
    }

    return order;
}

struct HitInfo {
    u64 case_idx;
    int base;
    u64 t;
    u64 k;
    u64 gcd_value;
    u64 factor;
    u64 order;
};

int main() {
    std::cout << "START EXPERIMENT 317\n";
    std::cout << "INTRINSIC HIT ORDER AUDIT\n";
    std::cout << "ARE THE CLEAN HITS JUST SMALL MULTIPLICATIVE-ORDER COLLISIONS?\n\n";

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

    u64 total_tests = 0;
    u64 total_hits = 0;

    u64 p_hits = 0;
    u64 q_hits = 0;

    u64 explained_by_order = 0;
    u64 unexplained_hits = 0;

    u64 p_small_order = 0;
    u64 q_small_order = 0;

    std::vector<HitInfo> hits;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_tests = 0;
        u64 case_hits = 0;
        u64 case_explained = 0;
        u64 case_unexplained = 0;

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

                    const i128 raw =
                        evaluate_polynomial(coeff, t);

                    const auto exponents =
                        intrinsic_exponents(
                            coeff,
                            m,
                            D,
                            base,
                            t
                        );

                    /*
                     * Audit both P_b(t) and D+P_b(t).
                     */
                    const std::array<i128, 2> values = {
                        raw,
                        raw + static_cast<i128>(D)
                    };

                    for (i128 value : values) {
                        const u64 X =
                            signed_mod(value, N);

                        for (u64 k : exponents) {
                            ++total_tests;
                            ++case_tests;

                            const u64 power =
                                pow_mod(X, k, N);

                            const u64 delta =
                                (power >= 1)
                                    ? power - 1
                                    : N - (1 - power);

                            const u64 g =
                                gcd_u64(delta, N);

                            if (g == p || g == q) {
                                ++total_hits;
                                ++case_hits;

                                u64 factor;
                                u64 other;

                                if (g == p) {
                                    factor = p;
                                    other = q;
                                    ++p_hits;
                                } else {
                                    factor = q;
                                    other = p;
                                    ++q_hits;
                                }

                                /*
                                 * Determine whether k is exactly
                                 * explained by the multiplicative
                                 * order modulo the recovered factor.
                                 */
                                const u64 xm =
                                    X % factor;

                                const u64 order =
                                    exact_order_from_factor(
                                        xm,
                                        factor
                                    );

                                const bool explained =
                                    order != 0 &&
                                    k % order == 0;

                                if (explained) {
                                    ++explained_by_order;
                                    ++case_explained;
                                } else {
                                    ++unexplained_hits;
                                    ++case_unexplained;
                                }

                                if (order <= 60) {
                                    if (factor == p) {
                                        ++p_small_order;
                                    } else {
                                        ++q_small_order;
                                    }
                                }

                                HitInfo info{};
                                info.case_idx =
                                    static_cast<u64>(case_idx);
                                info.base = base;
                                info.t = t;
                                info.k = k;
                                info.gcd_value = g;
                                info.factor = factor;
                                info.order = order;

                                hits.push_back(info);

                                std::cout
                                    << "HIT"
                                    << " case=" << case_idx
                                    << " base=" << base
                                    << " t=" << t
                                    << " k=" << k
                                    << " gcd=" << g
                                    << " order_mod_factor="
                                    << order
                                    << " explained="
                                    << (explained ? 1 : 0)
                                    << "\n";
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
            << "\n"
            << "tests=" << case_tests << "\n"
            << "hits=" << case_hits << "\n"
            << "explained=" << case_explained << "\n"
            << "unexplained=" << case_unexplained << "\n\n";
    }

    std::cout << "============================\n";
    std::cout << "TOTAL\n";

    std::cout
        << "tests=" << total_tests << "\n"
        << "hits=" << total_hits << "\n"
        << "p_hits=" << p_hits << "\n"
        << "q_hits=" << q_hits << "\n"
        << "explained_by_order=" << explained_by_order << "\n"
        << "unexplained_hits=" << unexplained_hits << "\n"
        << "p_small_order=" << p_small_order << "\n"
        << "q_small_order=" << q_small_order << "\n";

    std::cout << "\nFINISHED EXPERIMENT 317\n";

    return 0;
}
