#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;

struct Case {
    u64 p;
    u64 q;
};

struct Poly {
    std::vector<mpz_class> c; // ascending
};

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

static u64 isqrt_u64(u64 n) {
    u64 lo = 0;
    u64 hi = std::min<u64>(
        n,
        1ULL << 32
    );

    while (lo <= hi) {
        const u64 mid =
            lo + (hi - lo) / 2;

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

static u64 gcd_u64(
    u64 a,
    u64 b
) {
    return std::gcd(a, b);
}

static Poly trim(Poly a) {
    while (
        a.c.size() > 1 &&
        a.c.back() == 0
    ) {
        a.c.pop_back();
    }

    if (a.c.empty()) {
        a.c.push_back(0);
    }

    return a;
}

static int degree(
    const Poly& p
) {
    if (
        p.c.size() == 1 &&
        p.c[0] == 0
    ) {
        return -1;
    }

    return static_cast<int>(
        p.c.size()
    ) - 1;
}

static Poly multiply(
    const Poly& a,
    const Poly& b
) {
    Poly result;

    result.c.assign(
        a.c.size() +
        b.c.size() -
        1,
        0
    );

    for (
        std::size_t i = 0;
        i < a.c.size();
        ++i
    ) {
        for (
            std::size_t j = 0;
            j < b.c.size();
            ++j
        ) {
            result.c[i + j] +=
                a.c[i] * b.c[j];
        }
    }

    return trim(
        std::move(result)
    );
}

static Poly derivative(
    const Poly& p
) {
    const int d =
        degree(p);

    if (d <= 0) {
        return Poly{{0}};
    }

    Poly result;

    result.c.resize(
        static_cast<std::size_t>(d)
    );

    for (int i = 1; i <= d; ++i) {
        result.c[i - 1] =
            p.c[i] * i;
    }

    return trim(
        std::move(result)
    );
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
            result *=
                (x % base) + 1;
            x /= base;
        }

        return result;
    }

    std::vector<u64> mdigits;
    std::vector<u64> ydigits;

    u64 mx = m;
    u64 yx = y;

    while (
        mx > 0 ||
        yx > 0
    ) {
        mdigits.push_back(
            mx % base
        );

        ydigits.push_back(
            yx % base
        );

        mx /= base;
        yx /= base;
    }

    const std::size_t L =
        std::max(
            mdigits.size(),
            ydigits.size()
        );

    while (
        mdigits.size() < L
    ) {
        mdigits.push_back(0);
    }

    while (
        ydigits.size() < L
    ) {
        ydigits.push_back(0);
    }

    std::vector<u64> W(
        L + 1,
        1
    );

    for (
        std::size_t i = 0;
        i < L;
        ++i
    ) {
        W[i + 1] =
            W[i] *
            (mdigits[i] + 1);
    }

    int h =
        static_cast<int>(L) - 1;

    while (
        h >= 0 &&
        ydigits[h] ==
            mdigits[h]
    ) {
        --h;
    }

    if (h < 0) {
        return W[L];
    }

    u64 result = 0;

    for (
        int i =
            static_cast<int>(L) - 1;
        i > h;
        --i
    ) {
        result +=
            ydigits[i] * W[i];
    }

    result +=
        ydigits[h] * W[h];

    if (h == 0) {
        return result;
    }

    u64 low_m = 0;
    u64 low_y = 0;
    u64 power = 1;

    for (int i = 0; i < h; ++i) {
        low_m +=
            mdigits[i] * power;

        low_y +=
            ydigits[i] * power;

        power *=
            static_cast<u64>(base);
    }

    result +=
        prefix_miss(
            low_m,
            low_y,
            base
        );

    return result;
}

static Poly coefficient_polynomial(
    u64 m,
    u64 y,
    int base
) {
    std::vector<u64> digits;

    u64 x = m;

    while (x > 0) {
        digits.push_back(
            x % base
        );

        x /= base;
    }

    if (digits.empty()) {
        digits.push_back(0);
    }

    const int L =
        static_cast<int>(
            digits.size()
        );

    std::vector<u64> W(
        L + 1,
        1
    );

    for (int i = 0; i < L; ++i) {
        W[i + 1] =
            W[i] *
            (digits[i] + 1);
    }

    Poly result;

    result.c.resize(
        static_cast<std::size_t>(L)
    );

    u64 power_r = 1;

    for (int r = 0; r < L; ++r) {
        const u64 Q =
            y / power_r;

        const u64 R =
            y % power_r;

        if (Q == 0) {
            result.c[r] = 0;
        } else if (
            Q > digits[r]
        ) {
            result.c[r] =
                mpz_class(
                    digits[r]
                ) *
                mpz_class(
                    W[r]
                );
        } else {
            const u64 low_m =
                m % power_r;

            const u64 A =
                prefix_miss(
                    low_m,
                    R,
                    base
                );

            result.c[r] =
                mpz_class(
                    Q - 1
                ) *
                mpz_class(
                    W[r]
                )
                +
                mpz_class(A);
        }

        if (r + 1 < L) {
            power_r *=
                static_cast<u64>(base);
        }
    }

    return trim(
        std::move(result)
    );
}

/*
 * Exact determinant via Bareiss over GMP integers.
 *
 * This is O(n^3), not permutation enumeration.
 */
static mpz_class determinant_bareiss(
    std::vector<
        std::vector<mpz_class>
    > a
) {
    const int n =
        static_cast<int>(a.size());

    if (n == 0) {
        return 1;
    }

    if (n == 1) {
        return a[0][0];
    }

    mpz_class previous = 1;
    int sign = 1;

    for (int k = 0; k < n - 1; ++k) {
        int pivot = k;

        while (
            pivot < n &&
            a[pivot][k] == 0
        ) {
            ++pivot;
        }

        if (pivot == n) {
            return 0;
        }

        if (pivot != k) {
            std::swap(
                a[pivot],
                a[k]
            );

            sign = -sign;
        }

        const mpz_class pivot_value =
            a[k][k];

        for (
            int i = k + 1;
            i < n;
            ++i
        ) {
            for (
                int j = k + 1;
                j < n;
                ++j
            ) {
                mpz_class numerator =
                    a[i][j] *
                    pivot_value;

                numerator -=
                    a[i][k] *
                    a[k][j];

                a[i][j] =
                    numerator /
                    previous;
            }
        }

        for (
            int i = k + 1;
            i < n;
            ++i
        ) {
            a[i][k] = 0;
        }

        previous =
            pivot_value;
    }

    return sign < 0
        ? -a[n - 1][n - 1]
        : a[n - 1][n - 1];
}

/*
 * Exact Sylvester resultant.
 */
static mpz_class resultant(
    const Poly& A,
    const Poly& B
) {
    const int m =
        degree(A);

    const int n =
        degree(B);

    if (m < 0 || n < 0) {
        return 0;
    }

    if (m == 0) {
        mpz_class r = 1;

        for (int i = 0; i < n; ++i) {
            r *= A.c[0];
        }

        return r;
    }

    if (n == 0) {
        mpz_class r = 1;

        for (int i = 0; i < m; ++i) {
            r *= B.c[0];
        }

        return r;
    }

    const int size = m + n;

    std::vector<
        std::vector<mpz_class>
    > S(
        size,
        std::vector<mpz_class>(
            size,
            0
        )
    );

    for (int row = 0;
         row < n;
         ++row) {

        for (int j = 0;
             j <= m;
             ++j) {

            S[row][row + j] =
                A.c[j];
        }
    }

    for (int row = 0;
         row < m;
         ++row) {

        for (int j = 0;
             j <= n;
             ++j) {

            S[n + row][row + j] =
                B.c[j];
        }
    }

    return determinant_bareiss(
        std::move(S)
    );
}

static u64 gcd_mpz_u64(
    const mpz_class& x,
    u64 N
) {
    mpz_class ax = x;

    if (ax < 0) {
        ax = -ax;
    }

    mpz_class modN = ax % N;

    return modN.get_ui()
        ? gcd_u64(
              modN.get_ui(),
              N
          )
        : N;
}

static std::string
mpz_to_string(
    const mpz_class& x
) {
    return x.get_str();
}

int main() {
    std::cout
        << "START EXPERIMENT 323\n"
        << "DISCRIMINANT / REPEATED-ROOT TEST\n"
        << "CAN THE RADIX PATH POLYNOMIAL BECOME SINGULAR MODULO A HIDDEN PRIME?\n\n";

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

    u64 tests = 0;

    u64 gcd1 = 0;
    u64 gcdN = 0;
    u64 nontrivial = 0;
    u64 p_only = 0;
    u64 q_only = 0;
    u64 other = 0;

    u64 zero_resultant = 0;

    u64 printed = 0;

    for (std::size_t case_idx = 0;
         case_idx < cases.size();
         ++case_idx) {

        const u64 p =
            cases[case_idx].p;

        const u64 q =
            cases[case_idx].q;

        const u64 N =
            p * q;

        const u64 s =
            isqrt_u64(N);

        u64 case_tests = 0;
        u64 case_p = 0;
        u64 case_q = 0;

        for (int offset = 1;
             offset <= 5;
             ++offset) {

            const u64 off =
                static_cast<u64>(offset);

            if (s < off) {
                continue;
            }

            const u64 m =
                s - off + 1;

            const u64 D =
                N - m * m;

            if (D == 0) {
                continue;
            }

            for (int bi = 0;
                 bi < static_cast<int>(
                     BASES.size());
                 ++bi) {

                const int base =
                    BASES[bi];

                const Poly P =
                    coefficient_polynomial(
                        m,
                        N,
                        base
                    );

                Poly A = P;

                /*
                 * D + P(T)
                 */
                A.c[0] += D;

                const std::array<
                    const Poly*,
                    2
                > polys = {
                    &P,
                    &A
                };

                for (int affine = 0;
                     affine < 2;
                     ++affine) {

                    const Poly& F =
                        *polys[affine];

                    const Poly dF =
                        derivative(F);

                    const mpz_class R =
                        resultant(
                            F,
                            dF
                        );

                    ++tests;
                    ++case_tests;

                    if (R == 0) {
                        ++zero_resultant;
                    }

                    const u64 g =
                        gcd_mpz_u64(
                            R,
                            N
                        );

                    if (g == 1) {
                        ++gcd1;
                    } else if (g == N) {
                        ++gcdN;
                    } else {
                        ++nontrivial;

                        if (g == p) {
                            ++p_only;
                            ++case_p;

                            if (printed < 20) {
                                std::cout
                                    << "P_HIT"
                                    << " case="
                                    << case_idx
                                    << " base="
                                    << base
                                    << " affine="
                                    << affine
                                    << " degree="
                                    << degree(F)
                                    << "\n";

                                ++printed;
                            }
                        } else if (g == q) {
                            ++q_only;
                            ++case_q;

                            if (printed < 20) {
                                std::cout
                                    << "Q_HIT"
                                    << " case="
                                    << case_idx
                                    << " base="
                                    << base
                                    << " affine="
                                    << affine
                                    << " degree="
                                    << degree(F)
                                    << "\n";

                                ++printed;
                            }
                        } else {
                            ++other;
                        }
                    }
                }
            }
        }

        std::cout
            << "CASE "
            << case_idx
            << " p="
            << p
            << " q="
            << q
            << "\n"
            << "tests="
            << case_tests
            << "\n"
            << "p_only="
            << case_p
            << "\n"
            << "q_only="
            << case_q
            << "\n\n";
    }

    std::cout
        << "============================\n"
        << "TOTAL\n"
        << "tests="
        << tests
        << "\n"
        << "gcd1="
        << gcd1
        << "\n"
        << "gcdN="
        << gcdN
        << "\n"
        << "nontrivial="
        << nontrivial
        << "\n"
        << "p_only="
        << p_only
        << "\n"
        << "q_only="
        << q_only
        << "\n"
        << "other="
        << other
        << "\n"
        << "zero_resultant="
        << zero_resultant
        << "\n";

    std::cout
        << "\nFINISHED EXPERIMENT 323\n";

    return 0;
}