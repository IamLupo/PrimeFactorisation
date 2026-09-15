from __future__ import annotations

import math


# ==============================================================================================================
# EXPERIMENT 402
# ==============================================================================================================

EXP_NO = 402

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT ELIMINATION OF d / DIRECT N,K -> FACTOR POLYNOMIAL")
print()
print("Objective:")
print("  1. Start from N and K only.")
print("  2. Use the previously established relation:")
print("       2Z^2 - (N+1-d)Z + 2N = 0")
print("  3. Solve this equation formally for d:")
print("       d = [2Z^2 - (N+1)Z + 2N] / Z")
print("  4. Eliminate d using:")
print("       d^2 = -4K - 3N^2 + 6N + 1")
print("  5. Clear denominators exactly.")
print("  6. Obtain a polynomial P(N,K,Z) with no d.")
print("  7. Verify that Z=p and Z=q are exact roots.")
print("  8. Determine the degree and coefficient structure.")
print("  9. Test whether the quartic has additional exact roots.")
print(" 10. Test whether the quartic factors using N alone, K alone,")
print("     or the known factor polynomial Z^2-(p+q)Z+N.")
print(" 11. Compare the eliminated polynomial with the direct")
print("     quadratic from Experiment 401.")
print()
print("No resultant is constructed.")
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


def bit_length(x: int) -> int:
    return abs(int(x)).bit_length()


def exact_sqrt(n: int) -> tuple[int, bool]:
    if n < 0:
        return 0, False

    r = math.isqrt(n)
    return r, (r * r == n)


def poly_eval(coeffs: list[int], x: int) -> int:
    value = 0

    for c in coeffs:
        value = value * x + c

    return value


# ==============================================================================================================
# MAIN
# ==============================================================================================================

def main() -> None:

    global_checks = {
        "d_square": True,
        "p_root": True,
        "q_root": True,
        "degree": True,
        "leading_coefficient": True,
        "coefficient_identity": True,
        "direct_quadratic_relation": True,
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

        # Established compressed invariant:
        #
        # d = N - 2S + 1

        d_true = N - 2 * S + 1

        D = (
            -4 * K
            - 3 * N * N
            + 6 * N
            + 1
        )

        d, D_square_ok = exact_sqrt(D)

        d_recovery_ok = (
            D_square_ok
            and d == d_true
        )

        # ------------------------------------------------------------------------------------------
        # ELIMINATION
        # ------------------------------------------------------------------------------------------
        #
        # 2Z^2 - (N+1-d)Z + 2N = 0
        #
        # Rearranged:
        #
        # 2Z^2 - (N+1)Z + 2N = dZ
        #
        # Therefore:
        #
        # d = [2Z^2 - (N+1)Z + 2N] / Z
        #
        # Squaring:
        #
        # [2Z^2 - (N+1)Z + 2N]^2
        #     = D Z^2
        #
        # where
        #
        # D = d^2.
        #
        # This is a quartic polynomial.

        # Let:
        #
        # A(Z) = 2Z^2 - (N+1)Z + 2N.
        #
        # Then:
        #
        # P(Z) = A(Z)^2 - D Z^2.

        a2 = 2
        a1 = -(N + 1)
        a0 = 2 * N

        # Square:
        #
        # (a2 Z^2 + a1 Z + a0)^2

        c4 = a2 * a2
        c3 = 2 * a2 * a1
        c2 = a1 * a1 + 2 * a2 * a0
        c1 = 2 * a1 * a0
        c0 = a0 * a0

        # Subtract D*Z^2

        c2 -= D

        coeffs = [c4, c3, c2, c1, c0]

        # ------------------------------------------------------------------------------------------
        # DEGREE / LEADING COEFFICIENT
        # ------------------------------------------------------------------------------------------

        degree = 4

        degree_ok = coeffs[0] != 0
        leading_ok = coeffs[0] == 4

        # ------------------------------------------------------------------------------------------
        # DIRECT ROOT CHECKS
        # ------------------------------------------------------------------------------------------

        P_p = poly_eval(coeffs, p)
        P_q = poly_eval(coeffs, q)

        p_root_ok = P_p == 0
        q_root_ok = P_q == 0

        # ------------------------------------------------------------------------------------------
        # DIRECT QUADRATIC
        # ------------------------------------------------------------------------------------------

        # Q(Z) = Z^2 - S Z + N

        Q_p = p * p - S * p + N
        Q_q = q * q - S * q + N

        direct_quadratic_ok = (
            Q_p == 0
            and Q_q == 0
        )

        # ------------------------------------------------------------------------------------------
        # EXACT ALGEBRAIC IDENTITY
        # ------------------------------------------------------------------------------------------
        #
        # The quartic should equal:
        #
        #   [2Z^2 - (N+1)Z + 2N]^2 - D Z^2
        #
        # and for Z=p,q:
        #
        #   A(p) = d*p
        #   A(q) = d*q
        #
        # because p and q satisfy the quadratic.
        #

        A_p = (
            2 * p * p
            - (N + 1) * p
            + 2 * N
        )

        A_q = (
            2 * q * q
            - (N + 1) * q
            + 2 * N
        )

        A_relation_p = A_p == d * p
        A_relation_q = A_q == d * q

        coefficient_identity_ok = (
            A_relation_p
            and A_relation_q
        )

        # ------------------------------------------------------------------------------------------
        # EXPECTED FACTORIZATION
        # ------------------------------------------------------------------------------------------
        #
        # Since:
        #
        #   A(Z)^2 - d^2 Z^2
        #
        # factors as:
        #
        #   (A(Z)-dZ)(A(Z)+dZ)
        #
        # the quartic is actually the product of two quadratics:
        #
        #   [2Z^2-(N+1+d)Z+2N]
        #
        #   [2Z^2-(N+1-d)Z+2N].
        #
        # One factor is the known factor quadratic.
        #
        # The second is its d-sign conjugate.

        # Coefficients of:
        #
        # Q_minus = 2Z^2-(N+1-d)Z+2N
        #
        # Q_plus  = 2Z^2-(N+1+d)Z+2N

        Q_minus = [
            2,
            -(N + 1 - d),
            2 * N,
        ]

        Q_plus = [
            2,
            -(N + 1 + d),
            2 * N,
        ]

        # Multiply quadratics exactly.

        product_coeffs = [
            Q_minus[0] * Q_plus[0],

            Q_minus[0] * Q_plus[1]
            + Q_minus[1] * Q_plus[0],

            Q_minus[0] * Q_plus[2]
            + Q_minus[1] * Q_plus[1]
            + Q_minus[2] * Q_plus[0],

            Q_minus[1] * Q_plus[2]
            + Q_minus[2] * Q_plus[1],

            Q_minus[2] * Q_plus[2],
        ]

        coefficient_factorization_ok = (
            product_coeffs == coeffs
        )

        # ------------------------------------------------------------------------------------------
        # ROOTS OF THE CONJUGATE QUADRATIC
        # ------------------------------------------------------------------------------------------
        #
        # The known factor quadratic is Q_minus.
        #
        # The conjugate factor Q_plus may have unrelated roots.
        #
        # Determine its discriminant:
        #
        #   Delta_plus = (N+1+d)^2 - 16N.
        #

        Delta_plus = (
            (N + 1 + d) ** 2
            - 16 * N
        )

        sqrt_plus, plus_square = exact_sqrt(Delta_plus)

        # If it is a square, compute the two exact roots if integral.

        plus_root_1 = None
        plus_root_2 = None
        plus_roots_integral = False

        if plus_square:
            num1 = (N + 1 + d) - sqrt_plus
            num2 = (N + 1 + d) + sqrt_plus

            if num1 % 4 == 0 and num2 % 4 == 0:
                plus_root_1 = num1 // 4
                plus_root_2 = num2 // 4
                plus_roots_integral = True

        # ------------------------------------------------------------------------------------------
        # CONJUGATE FACTOR ROOT TESTS
        # ------------------------------------------------------------------------------------------

        qplus_p = poly_eval(Q_plus, p)
        qplus_q = poly_eval(Q_plus, q)

        p_is_plus_root = qplus_p == 0
        q_is_plus_root = qplus_q == 0

        # Normally these should be false unless special degeneracy occurs.

        # ------------------------------------------------------------------------------------------
        # SPECIAL RELATION TO N
        # ------------------------------------------------------------------------------------------

        # Product of the two roots of either quadratic:
        #
        #   2Z^2 - ... Z + 2N
        #
        # is N.

        minus_product_ok = (
            (p * q) == N
        )

        # If Q_plus has integer roots, test their product.

        plus_product_ok = True

        if plus_roots_integral:
            plus_product_ok = (
                plus_root_1 * plus_root_2 == N
            )

        # ------------------------------------------------------------------------------------------
        # SIZE DATA
        # ------------------------------------------------------------------------------------------

        coefficient_digits = [digits(abs(c)) for c in coeffs]

        # ------------------------------------------------------------------------------------------
        # INSTANCE STATUS
        # ------------------------------------------------------------------------------------------

        instance_ok = all([
            D_square_ok,
            d_recovery_ok,

            degree_ok,
            leading_ok,

            p_root_ok,
            q_root_ok,

            coefficient_identity_ok,
            coefficient_factorization_ok,

            direct_quadratic_ok,
            minus_product_ok,

            plus_product_ok,
        ])

        # ------------------------------------------------------------------------------------------
        # PRINT
        # ------------------------------------------------------------------------------------------

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {N}")
        print(f"  S = {S}")
        print(f"  K = {K}")
        print()

        print("FIRST RADICAL")
        print(f"  D = {D}")
        print(f"  D = d^2 : {D_square_ok}")
        print(f"  recovered d = {d}")
        print(f"  true d = {d_true}")
        print(f"  d recovery = {d_recovery_ok}")
        print()

        print("ELIMINATED POLYNOMIAL")
        print("  P(Z) = [2Z^2-(N+1)Z+2N]^2 - D Z^2")
        print()
        print("  Expanded coefficients:")
        print(f"    Z^4 : {c4}")
        print(f"    Z^3 : {c3}")
        print(f"    Z^2 : {c2}")
        print(f"    Z^1 : {c1}")
        print(f"    Z^0 : {c0}")
        print()
        print(f"  degree = {degree}")
        print(f"  degree check = {degree_ok}")
        print(f"  leading coefficient = {c4}")
        print(f"  leading coefficient check = {leading_ok}")
        print()

        print("DIRECT ROOT CHECKS")
        print(f"  P(p) = {P_p}")
        print(f"  P(q) = {P_q}")
        print(f"  p exact root = {p_root_ok}")
        print(f"  q exact root = {q_root_ok}")
        print()

        print("DIRECT FACTOR QUADRATIC")
        print("  Q_minus(Z) = 2Z^2-(N+1-d)Z+2N")
        print("  Q_minus roots = p,q")
        print(f"  Q_minus(p) = {Q_p}")
        print(f"  Q_minus(q) = {Q_q}")
        print(f"  direct quadratic checks = {direct_quadratic_ok}")
        print()

        print("ELIMINATION FACTORIZATION")
        print("  P(Z) = Q_minus(Z) * Q_plus(Z)")
        print("  Q_plus(Z) = 2Z^2-(N+1+d)Z+2N")
        print(
            f"  coefficient-by-coefficient factorization = "
            f"{coefficient_factorization_ok}"
        )
        print()

        print("CONJUGATE QUADRATIC")
        print("  Q_plus(Z) = 2Z^2-(N+1+d)Z+2N")
        print(f"  Q_plus(p) = {qplus_p}")
        print(f"  Q_plus(q) = {qplus_q}")
        print(f"  p also Q_plus root = {p_is_plus_root}")
        print(f"  q also Q_plus root = {q_is_plus_root}")
        print()

        print("CONJUGATE DISCRIMINANT")
        print("  Delta_plus = (N+1+d)^2 - 16N")
        print(f"  Delta_plus = {Delta_plus}")
        print(f"  Delta_plus square = {plus_square}")
        if plus_square:
            print(f"  sqrt(Delta_plus) = {sqrt_plus}")
        print(f"  integral conjugate roots = {plus_roots_integral}")

        if plus_roots_integral:
            print(f"  conjugate root 1 = {plus_root_1}")
            print(f"  conjugate root 2 = {plus_root_2}")
            print(
                f"  conjugate root product = "
                f"{plus_root_1 * plus_root_2}"
            )
            print(
                f"  conjugate product == N = "
                f"{plus_product_ok}"
            )
        print()

        print("ROOT-SPECIFIC ELIMINATION IDENTITY")
        print("  A(Z) = 2Z^2-(N+1)Z+2N")
        print("  A(p) = d*p")
        print("  A(q) = d*q")
        print(f"  A(p) = {A_p}")
        print(f"  d*p = {d * p}")
        print(f"  A(q) = {A_q}")
        print(f"  d*q = {d * q}")
        print(f"  root identity checks = {coefficient_identity_ok}")
        print()

        print("SIZE DATA")
        print(f"  N digits        = {digits(N)}")
        print(f"  K digits        = {digits(K)}")
        print(f"  d digits        = {digits(d)}")
        print(f"  quartic Z^4     = {digits(c4)}")
        print(f"  quartic Z^3     = {digits(c3)}")
        print(f"  quartic Z^2     = {digits(c2)}")
        print(f"  quartic Z       = {digits(c1)}")
        print(f"  quartic const   = {digits(c0)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

        # ------------------------------------------------------------------------------------------
        # GLOBAL
        # ------------------------------------------------------------------------------------------

        global_checks["d_square"] &= D_square_ok
        global_checks["p_root"] &= p_root_ok
        global_checks["q_root"] &= q_root_ok
        global_checks["degree"] &= degree_ok
        global_checks["leading_coefficient"] &= leading_ok
        global_checks["coefficient_identity"] &= coefficient_identity_ok
        global_checks["direct_quadratic_relation"] &= (
            coefficient_factorization_ok
            and direct_quadratic_ok
        )
        global_checks["all_exact"] &= instance_ok

        rows.append({
            "instance": idx,
            "N_digits": digits(N),
            "K_digits": digits(K),
            "d_digits": digits(d),
            "c4_digits": digits(c4),
            "c3_digits": digits(c3),
            "c2_digits": digits(c2),
            "c1_digits": digits(c1),
            "c0_digits": digits(c0),
            "plus_square": plus_square,
            "plus_integral_roots": plus_roots_integral,
            "instance_ok": instance_ok,
        })

    # ==========================================================================================================
    # GLOBAL SUMMARY
    # ==========================================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("FIRST RADICAL")
    print(f"  all d^2 checks = {global_checks['d_square']}")
    print()

    print("ELIMINATED QUARTIC")
    print(f"  all p roots = {global_checks['p_root']}")
    print(f"  all q roots = {global_checks['q_root']}")
    print(f"  all degree checks = {global_checks['degree']}")
    print(
        f"  all leading-coefficient checks = "
        f"{global_checks['leading_coefficient']}"
    )
    print(
        f"  all root-identity checks = "
        f"{global_checks['coefficient_identity']}"
    )
    print()

    print("DIRECT FACTOR STRUCTURE")
    print(
        f"  all quartic factorization checks = "
        f"{global_checks['direct_quadratic_relation']}"
    )
    print()

    print("CONJUGATE DISCRIMINANT SIGNATURE")
    square_count = sum(1 for r in rows if r["plus_square"])
    integral_count = sum(1 for r in rows if r["plus_integral_roots"])

    print(
        f"  conjugate discriminant square in "
        f"{square_count}/{len(rows)} instances"
    )
    print(
        f"  conjugate quadratic has integral roots in "
        f"{integral_count}/{len(rows)} instances"
    )
    print()

    print("SIZE RANGES")
    print(
        f"  N digits = "
        f"{min(r['N_digits'] for r in rows)} .. "
        f"{max(r['N_digits'] for r in rows)}"
    )
    print(
        f"  K digits = "
        f"{min(r['K_digits'] for r in rows)} .. "
        f"{max(r['K_digits'] for r in rows)}"
    )
    print(
        f"  d digits = "
        f"{min(r['d_digits'] for r in rows)} .. "
        f"{max(r['d_digits'] for r in rows)}"
    )
    print(
        f"  Z^4 coefficient digits = "
        f"{min(r['c4_digits'] for r in rows)} .. "
        f"{max(r['c4_digits'] for r in rows)}"
    )
    print(
        f"  Z^3 coefficient digits = "
        f"{min(r['c3_digits'] for r in rows)} .. "
        f"{max(r['c3_digits'] for r in rows)}"
    )
    print(
        f"  Z^2 coefficient digits = "
        f"{min(r['c2_digits'] for r in rows)} .. "
        f"{max(r['c2_digits'] for r in rows)}"
    )
    print(
        f"  Z coefficient digits = "
        f"{min(r['c1_digits'] for r in rows)} .. "
        f"{max(r['c1_digits'] for r in rows)}"
    )
    print(
        f"  constant digits = "
        f"{min(r['c0_digits'] for r in rows)} .. "
        f"{max(r['c0_digits'] for r in rows)}"
    )
    print()

    print("INSTANCE STATUS")
    for row in rows:
        print(
            f"  instance {row['instance']:>2}: "
            f"{row['instance_ok']}"
        )

    all_pass = all(row["instance_ok"] for row in rows)

    # ==========================================================================================================
    # FINAL STATUS
    # ==========================================================================================================

    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINAL STATUS")
    print("=" * 110)
    print()

    print("Eliminating d from the Experiment 401 quadratic gives:")
    print()
    print("  [2Z^2-(N+1)Z+2N]^2")
    print("      - [-4K-3N^2+6N+1] Z^2 = 0")
    print()
    print("Thus the entire factor reconstruction can be expressed")
    print("using only N and K, without d appearing explicitly.")
    print()
    print("The exact factorization is:")
    print()
    print("  P(Z)")
    print("   = [2Z^2-(N+1-d)Z+2N]")
    print("     [2Z^2-(N+1+d)Z+2N]")
    print()
    print("The first quadratic contains the true factors p and q.")
    print("The second is the d-conjugate branch introduced by")
    print("eliminating the sign of the square root.")
    print()
    print("Therefore the crucial structural question becomes:")
    print()
    print("  Does the conjugate quadratic have any arithmetic meaning,")
    print("  or is it purely an elimination artifact?")
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {all_pass}")
    print()
    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
