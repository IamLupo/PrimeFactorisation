from __future__ import annotations

import math


# ==============================================================================================================
# EXPERIMENT 399
# ==============================================================================================================

EXP_NO = 399

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT SINGLE-IDENTITY FACTORIZATION OF THE N,d,g CHAIN")
print()
print("Objective:")
print("  1. Recover d from:")
print("       D = -4K - 3N^2 + 6N + 1 = d^2")
print("  2. Recover S from:")
print("       S = (N+1-d)/2")
print("  3. Recover the factor gap:")
print("       g = sqrt(S^2 - 4N)")
print("  4. Test the compressed identity:")
print("       (N-d-7)^2 - 4g^2 = 16(d+3)")
print("  5. Factor it exactly as:")
print("       (N-d-7-2g)(N-d-7+2g) = 16(d+3)")
print("  6. Divide the two factors by 4 and test:")
print("       (N-d-7-2g)/4 = p-2")
print("       (N-d-7+2g)/4 = q-2")
print("  7. Determine whether the entire factor pair is encoded")
print("     directly by the single difference-of-squares identity.")
print("  8. Test equivalent factorizations involving only N,d,g.")
print("  9. Measure whether the new factorized form is structurally")
print("     simpler than the two successive quadratic reconstructions.")
print()
print("No resultants are constructed.")
print("No symbolic multivariate factoring is performed.")
print("Only exact integer arithmetic is used.")
print()


# ==============================================================================================================
# TEST INSTANCES
# ==============================================================================================================

INSTANCES = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
]


# ==============================================================================================================
# HELPERS
# ==============================================================================================================

def digits(x: int) -> int:
    x = abs(int(x))
    return 1 if x == 0 else len(str(x))


def bitlen(x: int) -> int:
    x = abs(int(x))
    return 1 if x == 0 else x.bit_length()


def exact_sqrt(n: int) -> tuple[int, bool]:
    if n < 0:
        return 0, False

    r = math.isqrt(n)
    return r, r * r == n


def gcd_abs(a: int, b: int) -> int:
    return math.gcd(abs(a), abs(b))


# ==============================================================================================================
# MAIN
# ==============================================================================================================

def main() -> None:

    checks = {
        "d_recovery": True,
        "S_recovery": True,
        "g_recovery": True,
        "difference_square": True,
        "factorized_identity": True,
        "div4_left": True,
        "div4_right": True,
        "p_shift_recovery": True,
        "q_shift_recovery": True,
        "factor_recovery": True,
        "shifted_product": True,
        "alternate_forms": True,
        "all_exact": True,
    }

    rows = []

    for idx, (p, q) in enumerate(INSTANCES, start=1):

        # ------------------------------------------------------------------------------------------
        # BASIC DATA
        # ------------------------------------------------------------------------------------------

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        # Exact compressed variable.
        d_true = N - 2 * S + 1

        D = -4 * K - 3 * N * N + 6 * N + 1
        d, d_square_ok = exact_sqrt(D)

        d_ok = d_square_ok and d == d_true

        # ------------------------------------------------------------------------------------------
        # RECOVER S
        # ------------------------------------------------------------------------------------------

        S_num = N + 1 - d
        S_integral = (S_num % 2 == 0)

        S_rec = S_num // 2 if S_integral else None
        S_ok = S_integral and S_rec == S

        # ------------------------------------------------------------------------------------------
        # RECOVER GAP g
        # ------------------------------------------------------------------------------------------

        gap_discriminant = S * S - 4 * N
        g, g_square_ok = exact_sqrt(gap_discriminant)

        g_true = q - p
        g_ok = g_square_ok and g == g_true

        # ------------------------------------------------------------------------------------------
        # CENTRAL DIFFERENCE-OF-SQUARES IDENTITY
        # ------------------------------------------------------------------------------------------
        #
        # Starting with:
        #
        #   g^2 = S^2 - 4N
        #
        # and:
        #
        #   d = N - 2S + 1
        #
        # one obtains:
        #
        #   (N-d-7)^2 - 4g^2 = 16(d+3)
        #
        # because:
        #
        #   N-d-7 = 2S-8 = 2(S-4)
        #
        # and:
        #
        #   (S-4)^2 - g^2
        #      = S^2 - 8S + 16 - (S^2 - 4N)
        #      = 4N - 8S + 16
        #      = 4(d+3)
        #

        A = N - d - 7

        difference_lhs = A * A - 4 * g * g
        difference_rhs = 16 * (d + 3)

        difference_square_ok = (
            difference_lhs == difference_rhs
        )

        # ------------------------------------------------------------------------------------------
        # FACTORIZED FORM
        # ------------------------------------------------------------------------------------------

        factor_minus = A - 2 * g
        factor_plus = A + 2 * g

        factorized_lhs = factor_minus * factor_plus
        factorized_rhs = 16 * (d + 3)

        factorized_ok = (
            factorized_lhs == factorized_rhs
            and factorized_lhs == difference_lhs
        )

        # ------------------------------------------------------------------------------------------
        # DIVISIBILITY BY 4
        # ------------------------------------------------------------------------------------------

        minus_div4 = (factor_minus % 4 == 0)
        plus_div4 = (factor_plus % 4 == 0)

        p_shift_rec = factor_minus // 4 if minus_div4 else None
        q_shift_rec = factor_plus // 4 if plus_div4 else None

        p_shift_ok = minus_div4 and p_shift_rec == p - 2
        q_shift_ok = plus_div4 and q_shift_rec == q - 2

        # ------------------------------------------------------------------------------------------
        # DIRECT FACTOR RECOVERY
        # ------------------------------------------------------------------------------------------

        p_rec = p_shift_rec + 2 if p_shift_rec is not None else None
        q_rec = q_shift_rec + 2 if q_shift_rec is not None else None

        factor_recovery_ok = (
            p_rec == p
            and q_rec == q
            and p_rec * q_rec == N
        )

        # ------------------------------------------------------------------------------------------
        # SHIFTED PRODUCT
        # ------------------------------------------------------------------------------------------

        shifted_product_ok = (
            d + 3 == (p - 2) * (q - 2)
        )

        # ------------------------------------------------------------------------------------------
        # ALTERNATE FORMS
        # ------------------------------------------------------------------------------------------
        #
        # Since:
        #
        #   A = N-d-7 = 2(S-4)
        #   g = q-p
        #
        # we also expect:
        #
        #   A - 2g = 4(p-2)
        #   A + 2g = 4(q-2)
        #
        # and:
        #
        #   A/2 = S-4
        #
        # Verify all equivalent forms exactly.

        alternate_1 = (
            A == 2 * (S - 4)
        )

        alternate_2 = (
            factor_minus == 4 * (p - 2)
        )

        alternate_3 = (
            factor_plus == 4 * (q - 2)
        )

        alternate_4 = (
            factor_minus * factor_plus
            == 16 * (p - 2) * (q - 2)
        )

        alternate_5 = (
            difference_lhs
            == 16 * ((p - 2) * (q - 2))
        )

        alternate_forms_ok = all([
            alternate_1,
            alternate_2,
            alternate_3,
            alternate_4,
            alternate_5,
        ])

        # ------------------------------------------------------------------------------------------
        # EXTRA GCD / DIVISIBILITY STRUCTURE
        # ------------------------------------------------------------------------------------------

        gcd_signatures = {
            "gcd(A,g)": gcd_abs(A, g),
            "gcd(A,d)": gcd_abs(A, d),
            "gcd(A,N)": gcd_abs(A, N),
            "gcd(p-2,q-2)": gcd_abs(p - 2, q - 2),
            "gcd(g,d+3)": gcd_abs(g, d + 3),
        }

        divisibility_signatures = {
            "4 | A-2g": minus_div4,
            "4 | A+2g": plus_div4,
            "(p-2) | (d+3)": ((d + 3) % (p - 2) == 0),
            "(q-2) | (d+3)": ((d + 3) % (q - 2) == 0),
        }

        # ------------------------------------------------------------------------------------------
        # INSTANCE STATUS
        # ------------------------------------------------------------------------------------------

        instance_ok = all([
            d_ok,
            S_ok,
            g_ok,
            difference_square_ok,
            factorized_ok,
            minus_div4,
            plus_div4,
            p_shift_ok,
            q_shift_ok,
            factor_recovery_ok,
            shifted_product_ok,
            alternate_forms_ok,
        ])

        # ------------------------------------------------------------------------------------------
        # PRINT INSTANCE
        # ------------------------------------------------------------------------------------------

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p      = {p}")
        print(f"  q      = {q}")
        print(f"  N      = {N}")
        print(f"  S      = {S}")
        print(f"  X      = {X}")
        print(f"  K      = {K}")
        print()

        print("COMPRESSED d")
        print("  D = -4K - 3N^2 + 6N + 1")
        print(f"  D = {D}")
        print(f"  D = d^2 = {d * d} : {d_square_ok}")
        print(f"  recovered d = {d}")
        print(f"  true d      = {d_true}")
        print(f"  exact d recovery = {d_ok}")
        print()

        print("S RECOVERY")
        print("  S = (N+1-d)/2")
        print(f"  numerator = {S_num}")
        print(f"  parity / integrality = {S_integral}")
        print(f"  reconstructed S = {S_rec}")
        print(f"  exact S recovery = {S_ok}")
        print()

        print("FACTOR GAP")
        print("  g = q-p")
        print(f"  g_true = {g_true}")
        print(f"  S^2 - 4N = {gap_discriminant}")
        print(f"  exact square = {g_square_ok}")
        print(f"  recovered g = {g}")
        print(f"  exact g recovery = {g_ok}")
        print()

        print("CENTRAL DIFFERENCE-OF-SQUARES IDENTITY")
        print("  A = N-d-7")
        print(f"  A = {A}")
        print()
        print("  A^2 - 4g^2 = 16(d+3)")
        print(f"  left  = {difference_lhs}")
        print(f"  right = {difference_rhs}")
        print(f"  exact identity = {difference_square_ok}")
        print()

        print("FACTORIZED CENTRAL IDENTITY")
        print("  (N-d-7-2g)(N-d-7+2g) = 16(d+3)")
        print(f"  factor_minus = {factor_minus}")
        print(f"  factor_plus  = {factor_plus}")
        print(f"  factorized lhs = {factorized_lhs}")
        print(f"  factorized rhs = {factorized_rhs}")
        print(f"  exact factorization = {factorized_ok}")
        print()

        print("DIVISIBILITY / SHIFTED FACTOR RECOVERY")
        print(f"  4 | factor_minus = {minus_div4}")
        print(f"  4 | factor_plus  = {plus_div4}")
        print()
        print("  (N-d-7-2g)/4 = p-2")
        print(f"  reconstructed p-2 = {p_shift_rec}")
        print(f"  true p-2          = {p - 2}")
        print(f"  exact = {p_shift_ok}")
        print()
        print("  (N-d-7+2g)/4 = q-2")
        print(f"  reconstructed q-2 = {q_shift_rec}")
        print(f"  true q-2          = {q - 2}")
        print(f"  exact = {q_shift_ok}")
        print()

        print("SHIFTED PRODUCT")
        print("  d+3 = (p-2)(q-2)")
        print(f"  d+3 = {d + 3}")
        print(f"  (p-2)(q-2) = {(p - 2) * (q - 2)}")
        print(f"  exact = {shifted_product_ok}")
        print()

        print("DIRECT FACTOR RECOVERY")
        print("  p = (N-d-7-2g)/4 + 2")
        print("  q = (N-d-7+2g)/4 + 2")
        print(f"  reconstructed p = {p_rec}")
        print(f"  reconstructed q = {q_rec}")
        print(f"  p*q == N = {factor_recovery_ok}")
        print()

        print("ALTERNATE FORM CHECKS")
        print(f"  N-d-7 = 2(S-4) = {alternate_1}")
        print(f"  factor_minus = 4(p-2) = {alternate_2}")
        print(f"  factor_plus  = 4(q-2) = {alternate_3}")
        print(f"  product = 16(p-2)(q-2) = {alternate_4}")
        print(f"  A^2-4g^2 = 16(p-2)(q-2) = {alternate_5}")
        print(f"  all alternate forms = {alternate_forms_ok}")
        print()

        print("GCD SIGNATURE")
        for name, value in gcd_signatures.items():
            print(f"  {name:<18} = {value}")
        print()

        print("DIVISIBILITY SIGNATURE")
        for name, value in divisibility_signatures.items():
            print(f"  {name:<24} = {value}")
        print()

        print("SIZE DATA")
        print(f"  N digits              = {digits(N)}")
        print(f"  K digits              = {digits(K)}")
        print(f"  d digits              = {digits(d)}")
        print(f"  g digits              = {digits(g)}")
        print(f"  A=N-d-7 digits        = {digits(A)}")
        print(f"  factor_minus digits   = {digits(factor_minus)}")
        print(f"  factor_plus digits    = {digits(factor_plus)}")
        print(f"  16(d+3) digits        = {digits(16 * (d + 3))}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

        # ------------------------------------------------------------------------------------------
        # GLOBAL CHECKS
        # ------------------------------------------------------------------------------------------

        checks["d_recovery"] &= d_ok
        checks["S_recovery"] &= S_ok
        checks["g_recovery"] &= g_ok
        checks["difference_square"] &= difference_square_ok
        checks["factorized_identity"] &= factorized_ok
        checks["div4_left"] &= minus_div4
        checks["div4_right"] &= plus_div4
        checks["p_shift_recovery"] &= p_shift_ok
        checks["q_shift_recovery"] &= q_shift_ok
        checks["factor_recovery"] &= factor_recovery_ok
        checks["shifted_product"] &= shifted_product_ok
        checks["alternate_forms"] &= alternate_forms_ok
        checks["all_exact"] &= instance_ok

        rows.append({
            "instance": idx,
            "N_digits": digits(N),
            "K_digits": digits(K),
            "d_digits": digits(d),
            "g_digits": digits(g),
            "A_digits": digits(A),
            "factor_minus_digits": digits(factor_minus),
            "factor_plus_digits": digits(factor_plus),
            "instance_ok": instance_ok,
        })

    # ==========================================================================================================
    # GLOBAL SUMMARY
    # ==========================================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("PRIMARY COMPRESSED CHAIN")
    print(f"  all d recoveries             = {checks['d_recovery']}")
    print(f"  all S recoveries             = {checks['S_recovery']}")
    print(f"  all g recoveries             = {checks['g_recovery']}")
    print()

    print("SINGLE DIFFERENCE-OF-SQUARES IDENTITY")
    print(f"  all A^2-4g^2 checks          = {checks['difference_square']}")
    print(f"  all factored identities      = {checks['factorized_identity']}")
    print()

    print("DIRECT SHIFTED-FACTOR RECOVERY")
    print(f"  all left /4 recoveries       = {checks['p_shift_recovery']}")
    print(f"  all right /4 recoveries      = {checks['q_shift_recovery']}")
    print(f"  all p,q recoveries           = {checks['factor_recovery']}")
    print()

    print("KNOWN PRODUCT STRUCTURE")
    print(f"  all d+3=(p-2)(q-2) checks    = {checks['shifted_product']}")
    print(f"  all alternate forms          = {checks['alternate_forms']}")
    print()

    print("SIZE RANGES")
    print(
        f"  N digits              = "
        f"{min(r['N_digits'] for r in rows)} .. "
        f"{max(r['N_digits'] for r in rows)}"
    )
    print(
        f"  K digits              = "
        f"{min(r['K_digits'] for r in rows)} .. "
        f"{max(r['K_digits'] for r in rows)}"
    )
    print(
        f"  d digits              = "
        f"{min(r['d_digits'] for r in rows)} .. "
        f"{max(r['d_digits'] for r in rows)}"
    )
    print(
        f"  g digits              = "
        f"{min(r['g_digits'] for r in rows)} .. "
        f"{max(r['g_digits'] for r in rows)}"
    )
    print(
        f"  A=N-d-7 digits        = "
        f"{min(r['A_digits'] for r in rows)} .. "
        f"{max(r['A_digits'] for r in rows)}"
    )
    print()

    print("INSTANCE STATUS")
    for row in rows:
        print(
            f"  instance {row['instance']:>2}: "
            f"{row['instance_ok']}"
        )

    all_pass = all(row["instance_ok"] for row in rows)

    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINAL STATUS")
    print("=" * 110)
    print()
    print("The central identity tested is:")
    print()
    print("  (N-d-7)^2 - 4g^2 = 16(d+3)")
    print()
    print("where:")
    print()
    print("  d = sqrt(-4K - 3N^2 + 6N + 1)")
    print("  S = (N+1-d)/2")
    print("  g = sqrt(S^2-4N) = q-p")
    print()
    print("The difference of squares factors exactly as:")
    print()
    print("  (N-d-7-2g)(N-d-7+2g) = 16(d+3)")
    print()
    print("and, for the tested semiprimes:")
    print()
    print("  N-d-7-2g = 4(p-2)")
    print("  N-d-7+2g = 4(q-2)")
    print()
    print("Hence the factor pair is reconstructed directly by:")
    print()
    print("  p = (N-d-7-2g)/4 + 2")
    print("  q = (N-d-7+2g)/4 + 2")
    print()
    print("This tests whether the two quadratic layers can be collapsed")
    print("into one exact difference-of-squares factorization.")
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {all_pass}")
    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
