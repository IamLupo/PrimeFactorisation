#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using i128 = __int128_t;

struct Case {
    u64 p;
    u64 q;
};

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

struct Cyclotomic {
    int k;
    std::vector<i128> coeff;
};

static Cyclotomic make_cyclotomic(int k) {
    switch (k) {
        case 2:
            return {2, {1, 1}};

        case 3:
            return {3, {1, 1, 1}};

        case 4:
            return {4, {1, 0, 1}};

        case 5:
            return {5, {1, 1, 1, 1, 1}};

        case 6:
            return {6, {1, -1, 1}};

        case 8:
            return {8, {1, 0, 0, 0, 1}};

        case 10:
            return {10, {1, -1, 1, -1, 1}};

        case 12:
            return {12, {1, 0, -1, 0, 1}};

        default:
            return {0, {}};
    }
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

    const int L =
        static_cast<int>(digits.size());

    std::vector<u64> W(L + 1, 1);

    for (int i = 0; i < L; ++i) {
        W[i + 1] =
            W[i] * (digits[i] + 1);
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
            const u64 low_m =
                m % power_r;

            const u64 A =
                prefix_miss(low_m, R, base);

            coeff[r] =
                static_cast<i128>(Q - 1) *
                    static_cast<i128>(W[r])
                + static_cast<i128>(A);
        }

        if (r + 1 < L) {
            power_r *= static_cast<u64>(base);
        }
    }

    return coeff;
}

static void trim(
    std::vector<i128>& a
) {
    while (a.size() > 1 && a.back() == 0) {
        a.pop_back();
    }
}

static bool divisible_by_cyclotomic(
    const std::vector<i128>& poly,
    const Cyclotomic& phi
) {
    if (poly.size() < phi.coeff.size()) {
        return false;
    }

    std::vector<i128> work = poly;

    trim(work);

    const int d =
        static_cast<int>(phi.coeff.size()) - 1;

    /*
     * Phi_d is monic, so ordinary integer
     * long division is exact when divisible.
     */
    for (int i =
             static_cast<int>(work.size()) - 1;
         i >= d;
         --i) {

        const i128 factor = work[i];

        if (factor == 0) {
            continue;
        }

        for (int j = 0; j <= d; ++j) {
            work[i - d + j] -=
                factor * phi.coeff[j];
        }
    }

    trim(work);

    return work.size() == 1 &&
           work[0] == 0;
}

static std::vector<i128> quotient_by_cyclotomic(
    const std::vector<i128>& poly,
    const Cyclotomic& phi
) {
    std::vector<i128> work = poly;

    trim(work);

    const int d =
        static_cast<int>(phi.coeff.size()) - 1;

    if (static_cast<int>(work.size()) <= d) {
        return {};
    }

    std::vector<i128> quotient(
        work.size() - d,
        0
    );

    for (int i =
             static_cast<int>(work.size()) - 1;
         i >= d;
         --i) {

        const i128 factor = work[i];

        quotient[i - d] = factor;

        if (factor == 0) {
            continue;
        }

        for (int j = 0; j <= d; ++j) {
            work[i - d + j] -=
                factor * phi.coeff[j];
        }
    }

    trim(work);

    if (!(work.size() == 1 &&
          work[0] == 0)) {
        return {};
    }

    trim(quotient);

    return quotient;
}

static bool exact_equal(
    const std::vector<i128>& a,
    const std::vector<i128>& b
) {
    std::size_t n =
        std::max(a.size(), b.size());

    for (std::size_t i = 0; i < n; ++i) {
        const i128 x =
            i < a.size() ? a[i] : 0;

        const i128 y =
            i < b.size() ? b[i] : 0;

        if (x != y) {
            return false;
        }
    }

    return true;
}

static void print_poly_head(
    const std::vector<i128>& p
) {
    const std::size_t limit =
        std::min<std::size_t>(p.size(), 8);

    std::cout << "[";

    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << ",";
        }

        /*
         * Only small coefficients are expected here.
         * Convert manually.
         */
        i128 x = p[i];

        if (x == 0) {
            std::cout << "0";
            continue;
        }

        if (x < 0) {
            std::cout << "-";
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

        std::reverse(s.begin(), s.end());

        std::cout << s;
    }

    if (p.size() > limit) {
        std::cout << ",...";
    }

    std::cout << "]";
}

int main() {
    std::cout << "START EXPERIMENT 319\n";
    std::cout << "EXACT CYCLOTOMIC DIVISIBILITY\n";
    std::cout << "ARE THE ZERO RESULTANTS CAUSED BY EXACT POLYNOMIAL FACTORS?\n\n";

    const std::array<Cyclotomic, 8> phis = {
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

    std::array<u64, 8> raw_divisible{};
    std::array<u64, 8> affine_divisible{};

    std::array<u64, 8> raw_count{};
    std::array<u64, 8> affine_count{};

    u64 total_polynomials = 0;
    u64 raw_exact_factorizations = 0;
    u64 affine_exact_factorizations = 0;

    u64 printed = 0;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_raw = 0;
        u64 case_affine = 0;

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
                const int base =
                    BASES[bi];

                const auto P =
                    coefficient_formula(m, N, base);

                std::vector<i128> A = P;

                A[0] += static_cast<i128>(D);

                ++total_polynomials;

                for (std::size_t pi = 0;
                     pi < phis.size();
                     ++pi) {

                    const auto& phi = phis[pi];

                    ++raw_count[pi];
                    ++affine_count[pi];

                    const bool raw =
                        divisible_by_cyclotomic(
                            P,
                            phi
                        );

                    const bool affine =
                        divisible_by_cyclotomic(
                            A,
                            phi
                        );

                    if (raw) {
                        ++raw_divisible[pi];
                        ++case_raw;
                        ++raw_exact_factorizations;

                        if (printed < 20) {
                            const auto quotient =
                                quotient_by_cyclotomic(
                                    P,
                                    phi
                                );

                            std::cout
                                << "RAW_DIVISIBLE"
                                << " case=" << case_idx
                                << " base=" << base
                                << " phi=" << phi.k
                                << " degree="
                                << (P.size() - 1)
                                << " quotient_degree=";

                            if (quotient.empty()) {
                                std::cout << -1;
                            } else {
                                std::cout
                                    << (quotient.size() - 1);
                            }

                            std::cout
                                << " P_head=";

                            print_poly_head(P);

                            std::cout << "\n";

                            ++printed;
                        }
                    }

                    if (affine) {
                        ++affine_divisible[pi];
                        ++case_affine;
                        ++affine_exact_factorizations;

                        if (printed < 20) {
                            std::cout
                                << "AFFINE_DIVISIBLE"
                                << " case=" << case_idx
                                << " base=" << base
                                << " phi=" << phi.k
                                << " degree="
                                << (A.size() - 1)
                                << "\n";

                            ++printed;
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
            << "raw_divisibility_hits="
            << case_raw
            << "\n"
            << "affine_divisibility_hits="
            << case_affine
            << "\n\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n";

    std::cout
        << "polynomials="
        << total_polynomials
        << "\n"
        << "raw_exact_factorizations="
        << raw_exact_factorizations
        << "\n"
        << "affine_exact_factorizations="
        << affine_exact_factorizations
        << "\n\n";

    std::cout << "PER CYCLOTOMIC\n";

    for (std::size_t i = 0;
         i < phis.size();
         ++i) {

        std::cout
            << "phi_" << phis[i].k
            << " raw="
            << raw_divisible[i]
            << "/"
            << raw_count[i]
            << " affine="
            << affine_divisible[i]
            << "/"
            << affine_count[i]
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT 319\n";

    return 0;
}
