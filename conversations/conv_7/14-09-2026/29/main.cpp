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

struct Cyclotomic {
    int k;

    // Ascending coefficients:
    // phi(x) = coeff[0] + coeff[1] x + ... + coeff[d] x^d
    std::vector<i128> coeff;

    int degree() const {
        return static_cast<int>(coeff.size()) - 1;
    }
};

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

static i128 eval_polynomial(
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

static Cyclotomic make_cyclotomic(int k) {
    switch (k) {
        case 2:
            // x + 1
            return {2, {1, 1}};

        case 3:
            // x^2 + x + 1
            return {3, {1, 1, 1}};

        case 4:
            // x^2 + 1
            return {4, {1, 0, 1}};

        case 5:
            // x^4 + x^3 + x^2 + x + 1
            return {5, {1, 1, 1, 1, 1}};

        case 6:
            // x^2 - x + 1
            return {6, {1, -1, 1}};

        case 8:
            // x^4 + 1
            return {8, {1, 0, 0, 0, 1}};

        case 10:
            // x^4 - x^3 + x^2 - x + 1
            return {10, {1, -1, 1, -1, 1}};

        case 12:
            // x^4 - x^2 + 1
            return {12, {1, 0, -1, 0, 1}};

        default:
            return {0, {}};
    }
}

static std::vector<i128> poly_trim(
    std::vector<i128> a
) {
    while (a.size() > 1 && a.back() == 0) {
        a.pop_back();
    }

    return a;
}

static std::vector<i128> poly_mul_x_mod(
    const std::vector<i128>& a,
    const Cyclotomic& phi
) {
    const int d = phi.degree();

    std::vector<i128> b(d + 1, 0);

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (i + 1 <= static_cast<std::size_t>(d)) {
            b[i + 1] += a[i];
        } else {
            // x^d = -(phi[0] + ... + phi[d-1] x^(d-1))
            for (int j = 0; j < d; ++j) {
                b[j] -= a[i] * phi.coeff[j];
            }
        }
    }

    b.resize(d);
    return poly_trim(b);
}

static std::vector<i128> poly_reduce_mod_cyclotomic(
    const std::vector<i128>& f,
    const Cyclotomic& phi
) {
    const int d = phi.degree();

    if (static_cast<int>(f.size()) <= d) {
        std::vector<i128> result = f;
        result.resize(d, 0);
        return result;
    }

    std::vector<std::vector<i128>> powers(
        f.size(),
        std::vector<i128>(d, 0)
    );

    powers[0][0] = 1;

    for (std::size_t i = 1; i < f.size(); ++i) {
        powers[i] =
            poly_mul_x_mod(powers[i - 1], phi);

        powers[i].resize(d, 0);
    }

    std::vector<i128> result(d, 0);

    for (std::size_t i = 0; i < f.size(); ++i) {
        for (int j = 0; j < d; ++j) {
            result[j] += f[i] * powers[i][j];
        }
    }

    return result;
}

static std::vector<i128> multiply_residue_mod_phi(
    const std::vector<i128>& a,
    const std::vector<i128>& b,
    const Cyclotomic& phi
) {
    const int d = phi.degree();

    std::vector<i128> product(
        a.size() + b.size() - 1,
        0
    );

    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            product[i + j] += a[i] * b[j];
        }
    }

    return poly_reduce_mod_cyclotomic(product, phi);
}

static i128 determinant_small(
    const std::vector<std::vector<i128>>& matrix
) {
    const int n = static_cast<int>(matrix.size());

    if (n == 0) {
        return 1;
    }

    std::vector<int> p(n);

    for (int i = 0; i < n; ++i) {
        p[i] = i;
    }

    i128 det = 0;

    do {
        i128 term = 1;
        int inversions = 0;

        for (int i = 0; i < n; ++i) {
            term *= matrix[i][p[i]];

            for (int j = i + 1; j < n; ++j) {
                if (p[i] > p[j]) {
                    ++inversions;
                }
            }
        }

        if ((inversions & 1) == 0) {
            det += term;
        } else {
            det -= term;
        }
    } while (std::next_permutation(p.begin(), p.end()));

    return det;
}

static i128 cyclotomic_resultant(
    const std::vector<i128>& f,
    const Cyclotomic& phi
) {
    const int d = phi.degree();

    /*
     * In Z[x]/(Phi):
     *
     *   r(x) = f(x) mod Phi(x)
     *
     * The resultant Res(f,Phi) is the determinant
     * of multiplication-by-r on this d-dimensional
     * quotient algebra.
     */
    const auto r =
        poly_reduce_mod_cyclotomic(f, phi);

    std::vector<std::vector<i128>> M(
        d,
        std::vector<i128>(d, 0)
    );

    std::vector<i128> x(d, 0);
    x[0] = 1;

    for (int col = 0; col < d; ++col) {
        const auto product =
            multiply_residue_mod_phi(r, x, phi);

        for (int row = 0; row < d; ++row) {
            M[row][col] =
                row < static_cast<int>(product.size())
                    ? product[row]
                    : 0;
        }

        x = poly_mul_x_mod(x, phi);
        x.resize(d, 0);
    }

    return determinant_small(M);
}

static u64 abs_i128_mod(
    i128 x,
    u64 mod
) {
    if (x < 0) {
        x = -x;
    }

    return static_cast<u64>(
        static_cast<u128>(x) % mod
    );
}

static std::string to_string_i128(i128 x) {
    if (x == 0) {
        return "0";
    }

    bool neg = x < 0;

    if (neg) {
        x = -x;
    }

    std::string s;

    while (x > 0) {
        s.push_back(
            static_cast<char>(
                '0' + static_cast<int>(x % 10)
            )
        );

        x /= 10;
    }

    if (neg) {
        s.push_back('-');
    }

    std::reverse(s.begin(), s.end());

    return s;
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

int main() {
    std::cout << "START EXPERIMENT 318\n";
    std::cout << "CYCLOTOMIC RESULTANT FINGERPRINTS\n";
    std::cout << "CAN THE PATH POLYNOMIAL HAVE A PRIME-SPECIFIC ROOT OF UNITY?\n\n";

    const std::array<Cyclotomic, 8> cyclotomics = {
        make_cyclotomic(2),
        make_cyclotomic(3),
        make_cyclotomic(4),
        make_cyclotomic(5),
        make_cyclotomic(6),
        make_cyclotomic(8),
        make_cyclotomic(10),
        make_cyclotomic(12)
    };

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

    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;
    u64 p_only = 0;
    u64 q_only = 0;
    u64 other = 0;

    u64 zero_resultants = 0;

    u64 raw_tests = 0;
    u64 affine_tests = 0;

    std::vector<int> zero_k;

    u64 printed_hit = 0;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_tests = 0;
        u64 case_p_hits = 0;
        u64 case_q_hits = 0;
        u64 case_zero = 0;

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

                const auto P =
                    coefficient_formula(m, N, base);

                for (int ti = 0; ti < 5; ++ti) {
                    const u64 t =
                        static_cast<u64>(T_VALUES[ti]);

                    /*
                     * We keep the two objects that appeared
                     * in previous experiments:
                     *
                     *   P_b(t)
                     *   D + P_b(t)
                     *
                     * but the resultant is computed from the
                     * entire polynomial P_b(T), not its value
                     * at t.
                     */
                    const std::array<
                        std::vector<i128>, 2
                    > polys = {{
                        P,
                        P
                    }};

                    (void)t;
                    (void)polys;

                    /*
                     * Construct D + P(T) by modifying only
                     * the constant coefficient.
                     */
                    std::vector<i128> affine = P;

                    if (affine.empty()) {
                        affine.push_back(
                            static_cast<i128>(D)
                        );
                    } else {
                        affine[0] +=
                            static_cast<i128>(D);
                    }

                    for (const Cyclotomic& phi :
                         cyclotomics) {

                        const i128 raw_resultant =
                            cyclotomic_resultant(
                                P,
                                phi
                            );

                        const i128 affine_resultant =
                            cyclotomic_resultant(
                                affine,
                                phi
                            );

                        ++raw_tests;
                        ++affine_tests;

                        ++total_tests;
                        ++total_tests;

                        ++case_tests;

                        const u64 raw_abs =
                            abs_i128_mod(
                                raw_resultant,
                                N
                            );

                        const u64 affine_abs =
                            abs_i128_mod(
                                affine_resultant,
                                N
                            );

                        if (raw_resultant == 0) {
                            ++zero_resultants;
                            ++case_zero;
                        }

                        if (affine_resultant == 0) {
                            ++zero_resultants;
                            ++case_zero;
                        }

                        const u64 g_raw =
                            gcd_u64(raw_abs, N);

                        const u64 old_p =
                            p_only;

                        const u64 old_q =
                            q_only;

                        classify(
                            g_raw,
                            p,
                            q,
                            gcd1,
                            gcdN,
                            nontrivial,
                            p_only,
                            q_only,
                            other
                        );

                        if (g_raw == p ||
                            g_raw == q) {

                            if (g_raw == p) {
                                ++case_p_hits;
                            } else {
                                ++case_q_hits;
                            }

                            if (printed_hit < 20) {
                                std::cout
                                    << "RAW_HIT"
                                    << " case=" << case_idx
                                    << " base=" << base
                                    << " phi=" << phi.k
                                    << " gcd=" << g_raw
                                    << " resultant="
                                    << to_string_i128(
                                           raw_resultant)
                                    << "\n";

                                ++printed_hit;
                            }
                        }

                        const u64 g_affine =
                            gcd_u64(affine_abs, N);

                        classify(
                            g_affine,
                            p,
                            q,
                            gcd1,
                            gcdN,
                            nontrivial,
                            p_only,
                            q_only,
                            other
                        );

                        if (g_affine == p ||
                            g_affine == q) {

                            if (g_affine == p) {
                                ++case_p_hits;
                            } else {
                                ++case_q_hits;
                            }

                            if (printed_hit < 20) {
                                std::cout
                                    << "AFFINE_HIT"
                                    << " case=" << case_idx
                                    << " base=" << base
                                    << " phi=" << phi.k
                                    << " gcd=" << g_affine
                                    << " resultant="
                                    << to_string_i128(
                                           affine_resultant)
                                    << "\n";

                                ++printed_hit;
                            }
                        }

                        (void)old_p;
                        (void)old_q;
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
            << "p_hits=" << case_p_hits << "\n"
            << "q_hits=" << case_q_hits << "\n"
            << "zero_resultants=" << case_zero << "\n\n";
    }

    std::cout << "============================\n";
    std::cout << "TOTAL\n";

    std::cout
        << "tests=" << total_tests << "\n"
        << "raw_tests=" << raw_tests << "\n"
        << "affine_tests=" << affine_tests << "\n"
        << "gcd1=" << gcd1 << "\n"
        << "gcdN=" << gcdN << "\n"
        << "nontrivial=" << nontrivial << "\n"
        << "p_only=" << p_only << "\n"
        << "q_only=" << q_only << "\n"
        << "other=" << other << "\n"
        << "zero_resultants=" << zero_resultants << "\n";

    std::cout << "\nFINISHED EXPERIMENT 318\n";

    return 0;
}
