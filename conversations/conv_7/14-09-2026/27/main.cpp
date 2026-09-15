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

static u64 abs_mod(
    i128 x,
    u64 mod
) {
    if (x < 0) {
        x = -x;
    }

    return signed_mod(x, mod);
}

static u64 support_size(
    const std::vector<i128>& coeff
) {
    u64 result = 0;

    for (i128 c : coeff) {
        if (c != 0) {
            ++result;
        }
    }

    return result;
}

static u64 coefficient_sum_mod(
    const std::vector<i128>& coeff,
    u64 mod
) {
    u64 result = 0;

    for (i128 c : coeff) {
        result += abs_mod(c, mod);

        if (result >= mod) {
            result %= mod;
        }
    }

    return result % mod;
}

static u64 weighted_sum_mod(
    const std::vector<i128>& coeff,
    u64 mod
) {
    u64 result = 0;

    for (std::size_t r = 0; r < coeff.size(); ++r) {
        const u64 c = abs_mod(coeff[r], mod);

        const u64 term =
            static_cast<u64>(
                (u128(r) * u128(c)) % mod
            );

        result += term;

        if (result >= mod) {
            result %= mod;
        }
    }

    return result % mod;
}

static std::vector<u64> build_intrinsic_exponents(
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

    const u64 sum97 =
        coefficient_sum_mod(coeff, 97);

    const u64 sum193 =
        coefficient_sum_mod(coeff, 193);

    const u64 weighted97 =
        weighted_sum_mod(coeff, 97);

    const u64 weighted193 =
        weighted_sum_mod(coeff, 193);

    /*
     * Every exponent here is derived only from:
     *
     *   coefficient vector
     *   m, D
     *   base, t
     *
     * No p, q, p-1 or q-1 information is used.
     */
    out.insert(2 + support);
    out.insert(2 + deg);
    out.insert(2 + static_cast<u64>(base));
    out.insert(2 + t);

    out.insert(2 + (sum97 % 31));
    out.insert(2 + (sum193 % 37));

    out.insert(2 + (weighted97 % 31));
    out.insert(2 + (weighted193 % 37));

    out.insert(2 + (m % 41));
    out.insert(2 + (D % 41));

    /*
     * Mixed polynomial/radix quantities.
     */
    out.insert(
        2 + ((support * static_cast<u64>(base)) % 43)
    );

    out.insert(
        2 + ((deg * t) % 43)
    );

    out.insert(
        2 + ((sum97 + weighted97) % 43)
    );

    out.insert(
        2 + ((sum193 + weighted193) % 43)
    );

    return std::vector<u64>(out.begin(), out.end());
}

static std::vector<u64> build_control_exponents(
    const std::vector<u64>& intrinsic
) {
    /*
     * Matched controls: for each intrinsic exponent k,
     * use nearby integers k-2,...,k+2.
     *
     * These controls know nothing about the polynomial.
     */
    std::set<u64> out;

    for (u64 k : intrinsic) {
        for (int delta = -2; delta <= 2; ++delta) {
            const i128 candidate =
                static_cast<i128>(k) + delta;

            if (candidate >= 2 &&
                candidate <= 200) {
                out.insert(static_cast<u64>(candidate));
            }
        }
    }

    return std::vector<u64>(out.begin(), out.end());
}

static void classify(
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

static bool divides(u64 a, u64 b) {
    return a != 0 && b % a == 0;
}

int main() {
    std::cout << "START EXPERIMENT 316\n";
    std::cout << "INTRINSIC RADIX EXPONENTS\n";
    std::cout << "CAN POLYNOMIAL-DERIVED SMALL EXPONENTS EXPOSE A HIDDEN FACTOR?\n\n";

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

    u64 intrinsic_tests = 0;
    u64 control_tests = 0;

    u64 intrinsic_gcd1 = 0;
    u64 intrinsic_gcdN = 0;
    u64 intrinsic_nontrivial = 0;
    u64 intrinsic_p_only = 0;
    u64 intrinsic_q_only = 0;
    u64 intrinsic_other = 0;

    u64 control_gcd1 = 0;
    u64 control_gcdN = 0;
    u64 control_nontrivial = 0;
    u64 control_p_only = 0;
    u64 control_q_only = 0;
    u64 control_other = 0;

    u64 intrinsic_p1_contaminated = 0;
    u64 intrinsic_q1_contaminated = 0;

    u64 intrinsic_clean_p_hits = 0;
    u64 intrinsic_clean_q_hits = 0;

    u64 printed_hit = 0;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_intrinsic_tests = 0;
        u64 case_control_tests = 0;
        u64 case_clean_p_hits = 0;
        u64 case_clean_q_hits = 0;

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

                const auto intrinsic_cache =
                    std::array<std::vector<u64>, 5>{};

                (void)intrinsic_cache;

                for (int ti = 0; ti < 5; ++ti) {
                    const u64 t =
                        static_cast<u64>(T_VALUES[ti]);

                    const i128 raw =
                        evaluate_polynomial(coeff, t);

                    const std::array<i128, 2> values = {
                        raw,
                        raw + static_cast<i128>(D)
                    };

                    const auto intrinsic =
                        build_intrinsic_exponents(
                            coeff,
                            m,
                            D,
                            base,
                            t
                        );

                    const auto controls =
                        build_control_exponents(intrinsic);

                    for (i128 value : values) {
                        const u64 X =
                            signed_mod(value, N);

                        /*
                         * Intrinsic exponents.
                         */
                        for (u64 k : intrinsic) {
                            ++intrinsic_tests;
                            ++intrinsic_tests;

                            ++case_intrinsic_tests;

                            const bool p1_divides =
                                divides(p - 1, k);

                            const bool q1_divides =
                                divides(q - 1, k);

                            if (p1_divides) {
                                ++intrinsic_p1_contaminated;
                            }

                            if (q1_divides) {
                                ++intrinsic_q1_contaminated;
                            }

                            const u64 power =
                                pow_mod(X, k, N);

                            const u64 delta =
                                (power >= 1)
                                    ? power - 1
                                    : N - (1 - power);

                            const u64 g =
                                gcd_u64(delta, N);

                            const u64 old_p =
                                intrinsic_p_only;

                            const u64 old_q =
                                intrinsic_q_only;

                            classify(
                                g,
                                p,
                                q,
                                intrinsic_gcd1,
                                intrinsic_gcdN,
                                intrinsic_nontrivial,
                                intrinsic_p_only,
                                intrinsic_q_only,
                                intrinsic_other
                            );

                            /*
                             * A clean hit is one where the exponent
                             * does NOT contain the full group order
                             * of either prime.
                             */
                            if (g == p &&
                                !p1_divides &&
                                !q1_divides) {
                                ++intrinsic_clean_p_hits;
                                ++case_clean_p_hits;

                                if (printed_hit < 20) {
                                    std::cout
                                        << "CLEAN_P_HIT"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " k=" << k
                                        << " gcd=" << g
                                        << "\n";

                                    ++printed_hit;
                                }
                            }

                            if (g == q &&
                                !p1_divides &&
                                !q1_divides) {
                                ++intrinsic_clean_q_hits;
                                ++case_clean_q_hits;

                                if (printed_hit < 20) {
                                    std::cout
                                        << "CLEAN_Q_HIT"
                                        << " case=" << case_idx
                                        << " base=" << base
                                        << " t=" << t
                                        << " k=" << k
                                        << " gcd=" << g
                                        << "\n";

                                    ++printed_hit;
                                }
                            }

                            (void)old_p;
                            (void)old_q;
                        }

                        /*
                         * Matched generic controls.
                         */
                        for (u64 k : controls) {
                            ++control_tests;
                            ++case_control_tests;

                            const u64 power =
                                pow_mod(X, k, N);

                            const u64 delta =
                                (power >= 1)
                                    ? power - 1
                                    : N - (1 - power);

                            const u64 g =
                                gcd_u64(delta, N);

                            classify(
                                g,
                                p,
                                q,
                                control_gcd1,
                                control_gcdN,
                                control_nontrivial,
                                control_p_only,
                                control_q_only,
                                control_other
                            );
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
            << "intrinsic_tests=" << case_intrinsic_tests
            << "\n"
            << "control_tests=" << case_control_tests
            << "\n"
            << "clean_p_hits=" << case_clean_p_hits
            << "\n"
            << "clean_q_hits=" << case_clean_q_hits
            << "\n\n";
    }

    std::cout << "============================\n";
    std::cout << "TOTAL\n";

    std::cout
        << "intrinsic_tests=" << intrinsic_tests
        << "\n"
        << "control_tests=" << control_tests
        << "\n\n";

    std::cout
        << "INTRINSIC\n"
        << "gcd1=" << intrinsic_gcd1 << "\n"
        << "gcdN=" << intrinsic_gcdN << "\n"
        << "nontrivial=" << intrinsic_nontrivial << "\n"
        << "p_only=" << intrinsic_p_only << "\n"
        << "q_only=" << intrinsic_q_only << "\n"
        << "other=" << intrinsic_other << "\n"
        << "p1_contaminated=" << intrinsic_p1_contaminated << "\n"
        << "q1_contaminated=" << intrinsic_q1_contaminated << "\n"
        << "clean_p_hits=" << intrinsic_clean_p_hits << "\n"
        << "clean_q_hits=" << intrinsic_clean_q_hits << "\n\n";

    std::cout
        << "CONTROL\n"
        << "gcd1=" << control_gcd1 << "\n"
        << "gcdN=" << control_gcdN << "\n"
        << "nontrivial=" << control_nontrivial << "\n"
        << "p_only=" << control_p_only << "\n"
        << "q_only=" << control_q_only << "\n"
        << "other=" << control_other << "\n";

    std::cout << "\nFINISHED EXPERIMENT 316\n";

    return 0;
}
