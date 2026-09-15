#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <gmpxx.h>

using u64 = std::uint64_t;

struct Case {
    u64 p;
    u64 q;
};

struct PolyZ {
    std::vector<mpz_class> c;
};

struct PolyQ {
    std::vector<mpq_class> c;
};

static const std::array<int, 6> BASES = {
    2, 3, 5, 7, 11, 13
};

/* =========================================================
 * GMP helpers
 * ========================================================= */

static mpz_class abs_mpz(mpz_class x) {
    return x < 0 ? -x : x;
}

static mpz_class gcd_mpz(
    mpz_class a,
    mpz_class b
) {
    a = abs_mpz(a);
    b = abs_mpz(b);

    mpz_class g;

    mpz_gcd(
        g.get_mpz_t(),
        a.get_mpz_t(),
        b.get_mpz_t()
    );

    return g;
}

static mpz_class lcm_mpz(
    const mpz_class& a,
    const mpz_class& b
) {
    if (a == 0 || b == 0) {
        return 0;
    }

    return abs_mpz(
        (a / gcd_mpz(a, b)) * b
    );
}

/* =========================================================
 * Polynomial helpers
 * ========================================================= */

static PolyZ trim_z(PolyZ a) {
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

static PolyQ trim_q(PolyQ a) {
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

static bool zero_z(const PolyZ& a) {
    return a.c.size() == 1 &&
           a.c[0] == 0;
}

static bool zero_q(const PolyQ& a) {
    return a.c.size() == 1 &&
           a.c[0] == 0;
}

static int degree_z(const PolyZ& a) {
    return zero_z(a)
        ? -1
        : static_cast<int>(a.c.size()) - 1;
}

static int degree_q(const PolyQ& a) {
    return zero_q(a)
        ? -1
        : static_cast<int>(a.c.size()) - 1;
}

static PolyQ to_q(const PolyZ& p) {
    std::vector<mpq_class> c(
        p.c.size()
    );

    for (
        std::size_t i = 0;
        i < p.c.size();
        ++i
    ) {
        c[i] = p.c[i];
    }

    return trim_q(
        PolyQ{std::move(c)}
    );
}

static PolyQ derivative_q(
    const PolyQ& p
) {
    const int d =
        degree_q(p);

    if (d <= 0) {
        return PolyQ{
            {mpq_class(0)}
        };
    }

    std::vector<mpq_class> c(d);

    for (int i = 1; i <= d; ++i) {
        c[i - 1] =
            p.c[i] * i;
    }

    return trim_q(
        PolyQ{std::move(c)}
    );
}

/* =========================================================
 * Polynomial division / GCD over Q
 * ========================================================= */

static PolyQ poly_divrem_q(
    const PolyQ& A,
    const PolyQ& B,
    PolyQ& remainder
) {
    PolyQ R =
        trim_q(A);

    if (zero_q(B)) {
        remainder = R;
        return PolyQ{
            {mpq_class(0)}
        };
    }

    const int db =
        degree_q(B);

    if (degree_q(R) < db) {
        remainder = R;
        return PolyQ{
            {mpq_class(0)}
        };
    }

    std::vector<mpq_class> quotient(
        static_cast<std::size_t>(
            degree_q(R) - db + 1
        ),
        0
    );

    while (
        !zero_q(R) &&
        degree_q(R) >= db
    ) {
        const int dr =
            degree_q(R);

        const int shift =
            dr - db;

        const mpq_class factor =
            R.c[dr] / B.c[db];

        quotient[shift] += factor;

        for (int j = 0; j <= db; ++j) {
            R.c[j + shift] -=
                factor * B.c[j];
        }

        R =
            trim_q(std::move(R));
    }

    remainder = R;

    return trim_q(
        PolyQ{std::move(quotient)}
    );
}

static PolyQ poly_mod_q(
    const PolyQ& A,
    const PolyQ& B
) {
    PolyQ remainder;

    poly_divrem_q(
        A,
        B,
        remainder
    );

    return remainder;
}

static PolyQ monic_q(
    PolyQ p
) {
    p = trim_q(
        std::move(p)
    );

    if (zero_q(p)) {
        return p;
    }

    const mpq_class lead =
        p.c.back();

    for (auto& x : p.c) {
        x /= lead;
    }

    return trim_q(
        std::move(p)
    );
}

static PolyQ gcd_q(
    PolyQ A,
    PolyQ B
) {
    A =
        monic_q(std::move(A));

    B =
        monic_q(std::move(B));

    while (!zero_q(B)) {
        PolyQ R =
            poly_mod_q(A, B);

        A =
            std::move(B);

        B =
            monic_q(std::move(R));
    }

    return monic_q(
        std::move(A)
    );
}

/* =========================================================
 * Rational polynomial -> primitive integer polynomial
 * ========================================================= */

static std::vector<mpz_class>
primitive_integer_coefficients(
    const PolyQ& p
) {
    mpz_class denominator_lcm = 1;

    for (const auto& x : p.c) {
        denominator_lcm =
            lcm_mpz(
                denominator_lcm,
                x.get_den()
            );
    }

    std::vector<mpz_class> c(
        p.c.size()
    );

    for (
        std::size_t i = 0;
        i < p.c.size();
        ++i
    ) {
        /*
         * Correct GMP conversion:
         * no xmpz() wrapper.
         */
        c[i] =
            p.c[i].get_num() *
            (
                denominator_lcm /
                p.c[i].get_den()
            );
    }

    mpz_class content = 0;

    for (const auto& x : c) {
        content =
            gcd_mpz(
                content,
                abs_mpz(x)
            );
    }

    if (content != 0) {
        for (auto& x : c) {
            x /= content;
        }
    }

    if (
        !c.empty() &&
        c.back() < 0
    ) {
        for (auto& x : c) {
            x = -x;
        }
    }

    return c;
}

/* =========================================================
 * Polynomial pretty printer
 * ========================================================= */

static std::string poly_string(
    const std::vector<mpz_class>& c
) {
    bool first = true;
    std::string out;

    for (
        int i =
            static_cast<int>(c.size()) - 1;
        i >= 0;
        --i
    ) {
        if (c[i] == 0) {
            continue;
        }

        const bool negative =
            c[i] < 0;

        const mpz_class magnitude =
            abs_mpz(c[i]);

        if (!first) {
            out +=
                negative
                    ? " - "
                    : " + ";
        } else if (negative) {
            out += "-";
        }

        const bool show_coeff =
            magnitude != 1 ||
            i == 0;

        if (show_coeff) {
            out +=
                magnitude.get_str();
        }

        if (i >= 1) {
            out += "T";

            if (i >= 2) {
                out += "^";
                out +=
                    std::to_string(i);
            }
        }

        first = false;
    }

    return first ? "0" : out;
}

/* =========================================================
 * Integer square root
 * ========================================================= */

static u64 isqrt_u64(
    u64 n
) {
    u64 lo = 0;

    u64 hi =
        std::min<u64>(
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

/* =========================================================
 * Correct prefix MISS oracle
 * ========================================================= */

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

    while (mdigits.size() < L) {
        mdigits.push_back(0);
    }

    while (ydigits.size() < L) {
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
            ydigits[i] *
            W[i];
    }

    result +=
        ydigits[h] *
        W[h];

    if (h == 0) {
        return result;
    }

    u64 low_m = 0;
    u64 low_y = 0;
    u64 power = 1;

    for (int i = 0; i < h; ++i) {
        low_m +=
            mdigits[i] *
            power;

        low_y +=
            ydigits[i] *
            power;

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

/* =========================================================
 * Exact coefficient polynomial
 * ========================================================= */

static PolyZ coefficient_polynomial(
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

    std::vector<mpz_class> coeff(
        L,
        0
    );

    u64 power_r = 1;

    for (int r = 0; r < L; ++r) {
        const u64 Q =
            y / power_r;

        const u64 R =
            y % power_r;

        if (Q == 0) {
            coeff[r] = 0;
        } else if (
            Q > digits[r]
        ) {
            coeff[r] =
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

            coeff[r] =
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

    return trim_z(
        PolyZ{
            std::move(coeff)
        }
    );
}

/* =========================================================
 * Main experiment
 * ========================================================= */

int main() {
    std::cout
        << "START EXPERIMENT 324\n"
        << "EXACT REPEATED-FACTOR EXTRACTION\n"
        << "WHAT IS THE COMMON FACTOR OF P(T) AND P'(T)?\n\n";

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
    u64 repeated_cases = 0;

    u64 repeated_degree_1 = 0;
    u64 repeated_degree_2 = 0;
    u64 repeated_degree_3plus = 0;

    std::map<int, u64> base_zero_count;
    std::map<int, u64> offset_zero_count;

    std::map<std::string, u64> factor_shape_count;

    u64 printed = 0;

    for (
        std::size_t case_idx = 0;
        case_idx < cases.size();
        ++case_idx
    ) {
        const u64 p =
            cases[case_idx].p;

        const u64 q =
            cases[case_idx].q;

        const u64 N =
            p * q;

        const u64 s =
            isqrt_u64(N);

        for (
            int offset = 1;
            offset <= 5;
            ++offset
        ) {
            const u64 off =
                static_cast<u64>(
                    offset
                );

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

            for (
                int bi = 0;
                bi < static_cast<int>(
                    BASES.size()
                );
                ++bi
            ) {
                const int base =
                    BASES[bi];

                const PolyZ P =
                    coefficient_polynomial(
                        m,
                        N,
                        base
                    );

                for (
                    int affine = 0;
                    affine < 2;
                    ++affine
                ) {
                    PolyZ F = P;

                    if (affine) {
                        F.c[0] += D;

                        F =
                            trim_z(
                                std::move(F)
                            );
                    }

                    ++total_tests;

                    const PolyQ Fq =
                        to_q(F);

                    const PolyQ Fprime =
                        derivative_q(Fq);

                    const PolyQ G =
                        gcd_q(
                            Fq,
                            Fprime
                        );

                    const int gd =
                        degree_q(G);

                    if (gd <= 0) {
                        continue;
                    }

                    ++repeated_cases;

                    if (gd == 1) {
                        ++repeated_degree_1;
                    } else if (gd == 2) {
                        ++repeated_degree_2;
                    } else {
                        ++repeated_degree_3plus;
                    }

                    ++base_zero_count[base];
                    ++offset_zero_count[offset];

                    const auto primitive =
                        primitive_integer_coefficients(
                            G
                        );

                    const std::string shape =
                        poly_string(
                            primitive
                        );

                    ++factor_shape_count[
                        shape
                    ];

                    if (printed < 40) {
                        std::cout
                            << "REPEATED_FACTOR"
                            << " case="
                            << case_idx
                            << " base="
                            << base
                            << " offset="
                            << offset
                            << " affine="
                            << affine
                            << " degree_F="
                            << degree_q(Fq)
                            << " gcd_degree="
                            << gd
                            << "\n"
                            << "factor="
                            << shape
                            << "\n\n";

                        ++printed;
                    }
                }
            }
        }
    }

    std::cout
        << "============================\n"
        << "TOTAL\n"
        << "tests="
        << total_tests
        << "\n"
        << "repeated_cases="
        << repeated_cases
        << "\n"
        << "repeated_degree_1="
        << repeated_degree_1
        << "\n"
        << "repeated_degree_2="
        << repeated_degree_2
        << "\n"
        << "repeated_degree_3plus="
        << repeated_degree_3plus
        << "\n\n";

    std::cout
        << "BY BASE\n";

    for (int base : BASES) {
        std::cout
            << "base="
            << base
            << " repeated="
            << base_zero_count[base]
            << "\n";
    }

    std::cout
        << "\nBY OFFSET\n";

    for (int offset = 1;
         offset <= 5;
         ++offset) {

        std::cout
            << "offset="
            << offset
            << " repeated="
            << offset_zero_count[offset]
            << "\n";
    }

    std::cout
        << "\nREPEATED-FACTOR SHAPES\n";

    for (
        const auto& entry :
        factor_shape_count
    ) {
        std::cout
            << "count="
            << entry.second
            << " factor="
            << entry.first
            << "\n";
    }

    std::cout
        << "\nFINISHED EXPERIMENT 324\n";

    return 0;
}