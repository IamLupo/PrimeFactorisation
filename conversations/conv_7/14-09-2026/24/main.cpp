#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

using u64 = std::uint64_t;
using i64 = std::int64_t;
using u128 = __uint128_t;
using i128 = __int128_t;

struct Case {
    u64 p;
    u64 q;
};

static const std::array<int, 6> BASES = {2, 3, 5, 7, 11, 13};
static const std::array<int, 5> T_VALUES = {2, 3, 5, 7, 11};

static void print_u128(u128 x) {
    if (x == 0) {
        std::cout << '0';
        return;
    }

    std::string s;
    while (x > 0) {
        s.push_back(char('0' + (x % 10)));
        x /= 10;
    }

    std::reverse(s.begin(), s.end());
    std::cout << s;
}

static u64 gcd_u64(u64 a, u64 b) {
    return std::gcd(a, b);
}

static u64 isqrt_u64(u64 n) {
    u64 lo = 0;
    u64 hi = std::min<u64>(n, 1ULL << 32);

    while (lo <= hi) {
        u64 mid = lo + (hi - lo) / 2;

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

static u64 add_mod(u64 a, u64 b, u64 mod) {
    return static_cast<u64>((u128(a) + u128(b)) % mod);
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

static u64 mod_inverse_prime(u64 a, u64 p) {
    return pow_mod(a, p - 2, p);
}

static u64 lower_weight(
    u64 m,
    int base,
    int digit_pos
) {
    u64 weight = 1;
    u64 x = m;

    for (int i = 0; i < digit_pos; ++i) {
        const u64 digit = x % base;
        weight *= (digit + 1);
        x /= base;
    }

    return weight;
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

        if (m == 0) {
            return 1;
        }

        return result;
    }

    u64 result = 0;
    u64 m_hi = m;
    u64 y_hi = y;

    std::vector<u64> mdigits;
    std::vector<u64> ydigits;

    while (m_hi > 0 || y_hi > 0) {
        mdigits.push_back(m_hi % base);
        ydigits.push_back(y_hi % base);

        m_hi /= base;
        y_hi /= base;
    }

    const std::size_t L = std::max(mdigits.size(), ydigits.size());

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
        u64 total = 1;
        for (std::size_t i = 0; i < L; ++i) {
            total *= mdigits[i] + 1;
        }
        return total;
    }

    for (int i = static_cast<int>(L) - 1; i > h; --i) {
        result += ydigits[i] * W[i];
    }

    result += ydigits[h] * W[h];

    u64 low_m = 0;
    u64 low_y = 0;
    u64 power = 1;

    for (int i = 0; i < h; ++i) {
        low_m += mdigits[i] * power;
        low_y += ydigits[i] * power;
        power *= base;
    }

    if (h == 0) {
        return result;
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

    for (int r = 0; r < L; ++r) {
        u64 power_r = 1;
        for (int i = 0; i < r; ++i) {
            power_r *= static_cast<u64>(base);
        }

        const u64 Q = y / power_r;
        const u64 R = y % power_r;

        if (Q == 0) {
            coeff[r] = 0;
            continue;
        }

        if (Q > digits[r]) {
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
    }

    return coeff;
}

static i128 evaluate_polynomial(
    const std::vector<i128>& coeff,
    i64 t
) {
    i128 result = 0;
    i128 power = 1;

    for (const i128 c : coeff) {
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

static u64 determinant_mod_permutation(
    const std::array<std::array<u64, 5>, 5>& a,
    u64 mod
) {
    static const std::array<std::array<int, 5>, 120> perms = [] {
        std::array<std::array<int, 5>, 120> out{};
        std::array<int, 5> p = {0, 1, 2, 3, 4};

        int idx = 0;

        do {
            out[idx++] = p;
        } while (std::next_permutation(p.begin(), p.end()));

        return out;
    }();

    u64 determinant = 0;

    for (const auto& p : perms) {
        u64 term = 1;

        for (int i = 0; i < 5; ++i) {
            term = mul_mod(term, a[i][p[i]], mod);
        }

        int inversions = 0;

        for (int i = 0; i < 5; ++i) {
            for (int j = i + 1; j < 5; ++j) {
                if (p[i] > p[j]) {
                    ++inversions;
                }
            }
        }

        if ((inversions & 1) == 0) {
            determinant = add_mod(determinant, term, mod);
        } else {
            if (determinant >= term) {
                determinant -= term;
            } else {
                determinant = mod - (term - determinant);
            }
        }
    }

    return determinant;
}

static int rank_mod(
    const std::array<std::array<u64, 5>, 6>& input,
    u64 prime
) {
    std::array<std::array<u64, 5>, 6> a = input;

    int rank = 0;

    for (int col = 0; col < 5 && rank < 6; ++col) {
        int pivot = -1;

        for (int row = rank; row < 6; ++row) {
            if (a[row][col] % prime != 0) {
                pivot = row;
                break;
            }
        }

        if (pivot == -1) {
            continue;
        }

        std::swap(a[rank], a[pivot]);

        const u64 inv = mod_inverse_prime(a[rank][col], prime);

        for (int j = col; j < 5; ++j) {
            a[rank][j] = mul_mod(a[rank][j], inv, prime);
        }

        for (int row = 0; row < 6; ++row) {
            if (row == rank) {
                continue;
            }

            const u64 factor = a[row][col] % prime;

            if (factor == 0) {
                continue;
            }

            for (int j = col; j < 5; ++j) {
                const u64 sub = mul_mod(factor, a[rank][j], prime);

                if (a[row][j] >= sub) {
                    a[row][j] -= sub;
                } else {
                    a[row][j] = prime - (sub - a[row][j]);
                }
            }
        }

        ++rank;
    }

    return rank;
}

static std::array<std::array<u64, 5>, 5> select_rows(
    const std::array<std::array<u64, 5>, 6>& matrix,
    int omitted
) {
    std::array<std::array<u64, 5>, 5> result{};

    int dst = 0;

    for (int src = 0; src < 6; ++src) {
        if (src == omitted) {
            continue;
        }

        result[dst++] = matrix[src];
    }

    return result;
}

static void analyze_matrix(
    const std::array<std::array<i128, 5>, 6>& matrix,
    u64 p,
    u64 q,
    bool add_D,
    u64 D,
    u64 N,
    u64& minor_gcd1,
    u64& minor_gcdN,
    u64& minor_nontrivial,
    u64& minor_p_only,
    u64& minor_q_only,
    u64& rank_diff,
    u64& rank_equal
) {
    std::array<std::array<u64, 5>, 6> M{};

    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 5; ++c) {
            i128 value = matrix[r][c];

            if (add_D) {
                value += static_cast<i128>(D);
            }

            M[r][c] = signed_mod(value, N);
        }
    }

    const int rp = rank_mod(M, p);
    const int rq = rank_mod(M, q);

    if (rp != rq) {
        ++rank_diff;
    } else {
        ++rank_equal;
    }

    for (int omitted = 0; omitted < 6; ++omitted) {
        const auto minor = select_rows(M, omitted);

        const u64 detN = determinant_mod_permutation(minor, N);
        const u64 g = gcd_u64(N, detN);

        if (g == 1) {
            ++minor_gcd1;
        } else if (g == N) {
            ++minor_gcdN;
        } else {
            ++minor_nontrivial;

            if (g == p) {
                ++minor_p_only;
            } else if (g == q) {
                ++minor_q_only;
            }
        }

        (void)p;
        (void)q;
    }
}

int main() {
    std::cout << "START EXPERIMENT 313\n";
    std::cout << "FINITE-FIELD RANK FINGERPRINTS\n";
    std::cout << "CAN RADIX FINGERPRINT MATRICES HAVE DIFFERENT RANK MODULO p AND q?\n\n";

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

    u64 raw_rank_diff = 0;
    u64 raw_rank_equal = 0;

    u64 raw_gcd1 = 0;
    u64 raw_gcdN = 0;
    u64 raw_nontrivial = 0;
    u64 raw_p_only = 0;
    u64 raw_q_only = 0;

    u64 affine_rank_diff = 0;
    u64 affine_rank_equal = 0;

    u64 affine_gcd1 = 0;
    u64 affine_gcdN = 0;
    u64 affine_nontrivial = 0;
    u64 affine_p_only = 0;
    u64 affine_q_only = 0;

    for (std::size_t case_idx = 0; case_idx < cases.size(); ++case_idx) {
        const u64 p = cases[case_idx].p;
        const u64 q = cases[case_idx].q;
        const u64 N = p * q;

        const u64 s = isqrt_u64(N);

        u64 case_raw_rank_diff = 0;
        u64 case_affine_rank_diff = 0;

        u64 case_raw_nontrivial = 0;
        u64 case_affine_nontrivial = 0;

        u64 expressions = 0;

        for (int offset = 1; offset <= 5; ++offset) {
            if (s < static_cast<u64>(offset)) {
                continue;
            }

            const u64 m = s - static_cast<u64>(offset) + 1;
            const u64 D = N - m * m;

            if (D == 0) {
                continue;
            }

            std::array<std::array<i128, 5>, 6> matrix{};

            for (int bi = 0; bi < 6; ++bi) {
                const int base = BASES[bi];

                const auto coeff =
                    coefficient_formula(m, N, base);

                for (int ti = 0; ti < 5; ++ti) {
                    matrix[bi][ti] =
                        evaluate_polynomial(coeff, T_VALUES[ti]);
                }
            }

            const u64 old_raw_rank_diff = raw_rank_diff;
            const u64 old_affine_rank_diff = affine_rank_diff;

            const u64 old_raw_nontrivial = raw_nontrivial;
            const u64 old_affine_nontrivial = affine_nontrivial;

            analyze_matrix(
                matrix,
                p,
                q,
                false,
                D,
                N,
                raw_gcd1,
                raw_gcdN,
                raw_nontrivial,
                raw_p_only,
                raw_q_only,
                raw_rank_diff,
                raw_rank_equal
            );

            analyze_matrix(
                matrix,
                p,
                q,
                true,
                D,
                N,
                affine_gcd1,
                affine_gcdN,
                affine_nontrivial,
                affine_p_only,
                affine_q_only,
                affine_rank_diff,
                affine_rank_equal
            );

            if (raw_rank_diff != old_raw_rank_diff) {
                ++case_raw_rank_diff;
            }

            if (affine_rank_diff != old_affine_rank_diff) {
                ++case_affine_rank_diff;
            }

            if (raw_nontrivial != old_raw_nontrivial) {
                ++case_raw_nontrivial;
            }

            if (affine_nontrivial != old_affine_nontrivial) {
                ++case_affine_nontrivial;
            }

            ++expressions;
            ++total_expressions;
        }

        std::cout << "CASE " << case_idx
                  << " p=" << p
                  << " q=" << q << "\n";

        std::cout << "expressions=" << expressions << "\n";
        std::cout << "raw_rank_diff=" << case_raw_rank_diff << "\n";
        std::cout << "affine_rank_diff=" << case_affine_rank_diff << "\n";
        std::cout << "raw_nontrivial_minor=" << case_raw_nontrivial << "\n";
        std::cout << "affine_nontrivial_minor=" << case_affine_nontrivial << "\n";
        std::cout << "\n";
    }

    std::cout << "============================\n";
    std::cout << "TOTAL\n";

    std::cout << "expressions="
              << total_expressions << "\n\n";

    std::cout << "RAW P_b(t) MATRIX\n";
    std::cout << "rank_diff=" << raw_rank_diff << "\n";
    std::cout << "rank_equal=" << raw_rank_equal << "\n";
    std::cout << "minor_gcd1=" << raw_gcd1 << "\n";
    std::cout << "minor_gcdN=" << raw_gcdN << "\n";
    std::cout << "minor_nontrivial=" << raw_nontrivial << "\n";
    std::cout << "minor_p_only=" << raw_p_only << "\n";
    std::cout << "minor_q_only=" << raw_q_only << "\n\n";

    std::cout << "AFFINE D + P_b(t) MATRIX\n";
    std::cout << "rank_diff=" << affine_rank_diff << "\n";
    std::cout << "rank_equal=" << affine_rank_equal << "\n";
    std::cout << "minor_gcd1=" << affine_gcd1 << "\n";
    std::cout << "minor_gcdN=" << affine_gcdN << "\n";
    std::cout << "minor_nontrivial=" << affine_nontrivial << "\n";
    std::cout << "minor_p_only=" << affine_p_only << "\n";
    std::cout << "minor_q_only=" << affine_q_only << "\n";

    std::cout << "\nFINISHED EXPERIMENT 313\n";
    return 0;
}
