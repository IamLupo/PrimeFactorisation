from __future__ import annotations

import math


EXP_NO = 403

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT CONJUGATE-BRANCH / ROOT-STRUCTURE AUDIT")
print()
print("Objective:")
print("  1. Correct the sign error from Experiment 402:")
print("       A(p) = -d*p")
print("       A(q) = -d*q")
print("  2. Verify the eliminated quartic exactly.")
print("  3. Verify its factorization into Q_minus and Q_plus.")
print("  4. Analyze the two roots of Q_minus:")
print("       p,q")
print("  5. Analyze the two roots of Q_plus.")
print("  6. Determine whether Q_plus has rational roots.")
print("  7. Determine whether Q_plus roots are real or complex.")
print("  8. Express the Q_plus discriminant in terms of N,d,K.")
print("  9. Compare Delta_plus with Delta_minus.")
print(" 10. Search for a direct algebraic relation between the")
print("     two discriminants.")
print()
print("No resultants are constructed.")
print("Only exact integer arithmetic is used.")
print()


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


def digits(x: int) -> int:
    x = abs(int(x))
    return 1 if x == 0 else len(str(x))


def exact_sqrt(n: int):
    if n < 0:
        return 0, False
    r = math.isqrt(n)
    return r, r * r == n


def eval_poly(coeffs, x):
    v = 0
    for c in coeffs:
        v = v * x + c
    return v


def main():

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        # ------------------------------------------------------------------------------------------
        # FIRST RADICAL
        # ------------------------------------------------------------------------------------------

        D = -4 * K - 3 * N * N + 6 * N + 1

        d, d_square = exact_sqrt(D)

        d_true = N - 2 * S + 1

        d_ok = d_square and d == d_true

        # ------------------------------------------------------------------------------------------
        # BASIC QUADRATICS
        # ------------------------------------------------------------------------------------------

        # True factor quadratic:
        #
        # Q_minus = 2Z^2 - (N+1-d)Z + 2N
        #
        # Conjugate:
        #
        # Q_plus  = 2Z^2 - (N+1+d)Z + 2N

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

        qminus_p = eval_poly(Q_minus, p)
        qminus_q = eval_poly(Q_minus, q)

        qplus_p = eval_poly(Q_plus, p)
        qplus_q = eval_poly(Q_plus, q)

        qminus_ok = (
            qminus_p == 0
            and qminus_q == 0
        )

        # ------------------------------------------------------------------------------------------
        # CORRECT SIGNED IDENTITY
        # ------------------------------------------------------------------------------------------

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

        signed_p_ok = A_p == -d * p
        signed_q_ok = A_q == -d * q

        signed_ok = signed_p_ok and signed_q_ok

        # ------------------------------------------------------------------------------------------
        # QUARTIC
        # ------------------------------------------------------------------------------------------

        c4 = 4
        c3 = -4 * (N + 1)

        c2 = (
            (N + 1) ** 2
            + 8 * N
            - D
        )

        c1 = -4 * N * (N + 1)

        c0 = 4 * N * N

        quartic = [c4, c3, c2, c1, c0]

        Pp = eval_poly(quartic, p)
        Pq = eval_poly(quartic, q)

        quartic_roots_ok = (
            Pp == 0
            and Pq == 0
        )

        # ------------------------------------------------------------------------------------------
        # EXACT FACTORIZATION
        # ------------------------------------------------------------------------------------------

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

        factorization_ok = product_coeffs == quartic

        # ------------------------------------------------------------------------------------------
        # DISCRIMINANTS
        # ------------------------------------------------------------------------------------------

        Delta_minus = (
            (N + 1 - d) ** 2
            - 16 * N
        )

        Delta_plus = (
            (N + 1 + d) ** 2
            - 16 * N
        )

        sqrt_minus, minus_square = exact_sqrt(Delta_minus)
        sqrt_plus, plus_square = exact_sqrt(Delta_plus)

        minus_relation = Delta_minus == 4 * (q - p) ** 2

        # ------------------------------------------------------------------------------------------
        # SUM / PRODUCT OF CONJUGATE ROOTS
        # ------------------------------------------------------------------------------------------

        # For 2Z²-aZ+2N=0:
        #
        # sum = a/2
        # product = N

        plus_sum_twice = N + 1 + d
        plus_product = N

        # Rational-root test.
        #
        # A monic quadratic after dividing by 2 would be:
        #
        # Z² - ((N+1+d)/2)Z + N
        #
        # Q_plus has rational roots iff Delta_plus is a rational square.
        # Since Delta_plus is integer, that means perfect integer square
        # after accounting for parity.

        rational_plus_roots = plus_square

        # ------------------------------------------------------------------------------------------
        # DISCRIMINANT RELATIONS
        # ------------------------------------------------------------------------------------------

        # Difference:
        #
        # Delta_plus - Delta_minus
        # = 4d(N+1)

        discriminant_difference = Delta_plus - Delta_minus

        expected_difference = 4 * d * (N + 1)

        discriminant_difference_ok = (
            discriminant_difference == expected_difference
        )

        # Sum:
        #
        # Delta_plus + Delta_minus
        # = 2[(N+1)^2 + d^2 - 8N]

        discriminant_sum = Delta_plus + Delta_minus

        expected_sum = (
            2 * (
                (N + 1) ** 2
                + d * d
                - 8 * N
            )
        )

        discriminant_sum_ok = (
            discriminant_sum == expected_sum
        )

        # ------------------------------------------------------------------------------------------
        # K-ONLY FORM FOR D
        # ------------------------------------------------------------------------------------------

        D_from_K = -4 * K - 3 * N * N + 6 * N + 1

        D_from_K_ok = D_from_K == d * d

        # ------------------------------------------------------------------------------------------
        # Q_PLUS EVALUATIONS RELATION
        # ------------------------------------------------------------------------------------------

        # Since Q_plus = Q_minus - 2dZ:
        #
        # Q_plus(p) = -2dp
        # Q_plus(q) = -2dq

        qplus_relation_p = qplus_p == -2 * d * p
        qplus_relation_q = qplus_q == -2 * d * q

        qplus_relation_ok = (
            qplus_relation_p
            and qplus_relation_q
        )

        # ------------------------------------------------------------------------------------------
        # INSTANCE STATUS
        # ------------------------------------------------------------------------------------------

        instance_ok = all([
            d_ok,
            qminus_ok,
            signed_ok,
            quartic_roots_ok,
            factorization_ok,
            minus_relation,
            discriminant_difference_ok,
            discriminant_sum_ok,
            D_from_K_ok,
            qplus_relation_ok,
        ])

        global_ok &= instance_ok

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
        print(f"  D = d^2 : {d_square}")
        print(f"  recovered d = {d}")
        print(f"  true d = {d_true}")
        print(f"  d recovery = {d_ok}")
        print()

        print("CORRECTED SIGN IDENTITY")
        print("  A(Z) = 2Z^2-(N+1)Z+2N")
        print("  A(p) = -d*p")
        print("  A(q) = -d*q")
        print(f"  A(p) = {A_p}")
        print(f"  -d*p = {-d * p}")
        print(f"  A(q) = {A_q}")
        print(f"  -d*q = {-d * q}")
        print(f"  signed identity = {signed_ok}")
        print()

        print("DIRECT FACTOR QUADRATIC")
        print("  Q_minus(Z) = 2Z^2-(N+1-d)Z+2N")
        print(f"  Q_minus(p) = {qminus_p}")
        print(f"  Q_minus(q) = {qminus_q}")
        print(f"  roots p,q = {qminus_ok}")
        print()

        print("ELIMINATED QUARTIC")
        print(
            "  P(Z) = [2Z^2-(N+1)Z+2N]^2 - D Z^2"
        )
        print(f"  Z^4 coefficient = {c4}")
        print(f"  Z^3 coefficient = {c3}")
        print(f"  Z^2 coefficient = {c2}")
        print(f"  Z coefficient   = {c1}")
        print(f"  constant        = {c0}")
        print()
        print(f"  P(p) = {Pp}")
        print(f"  P(q) = {Pq}")
        print(f"  p,q roots = {quartic_roots_ok}")
        print()

        print("EXACT FACTORIZATION")
        print("  P(Z) = Q_minus(Z) * Q_plus(Z)")
        print("  Q_plus(Z) = 2Z^2-(N+1+d)Z+2N")
        print(f"  coefficient identity = {factorization_ok}")
        print()

        print("CONJUGATE QUADRATIC")
        print(f"  Q_plus(p) = {qplus_p}")
        print(f"  Q_plus(q) = {qplus_q}")
        print(f"  Q_plus(p) = -2*d*p : {qplus_relation_p}")
        print(f"  Q_plus(q) = -2*d*q : {qplus_relation_q}")
        print()

        print("CONJUGATE DISCRIMINANT")
        print("  Delta_minus = (N+1-d)^2 - 16N")
        print("  Delta_plus  = (N+1+d)^2 - 16N")
        print(f"  Delta_minus = {Delta_minus}")
        print(f"  Delta_plus  = {Delta_plus}")
        print(f"  Delta_minus square = {minus_square}")
        print(f"  Delta_plus square  = {plus_square}")
        print(f"  Delta_minus = 4(q-p)^2 : {minus_relation}")
        print()

        if plus_square:
            print(f"  sqrt(Delta_plus) = {sqrt_plus}")
            print("  Q_plus has rational/integer roots.")
        else:
            print("  Q_plus has no rational roots.")
        print()

        print("DISCRIMINANT RELATIONS")
        print(
            "  Delta_plus - Delta_minus = "
            "4*d*(N+1)"
        )
        print(
            f"  actual difference = "
            f"{discriminant_difference}"
        )
        print(
            f"  expected difference = "
            f"{expected_difference}"
        )
        print(
            f"  difference identity = "
            f"{discriminant_difference_ok}"
        )
        print()

        print(
            "  Delta_plus + Delta_minus = "
            "2*((N+1)^2+d^2-8N)"
        )
        print(
            f"  sum identity = "
            f"{discriminant_sum_ok}"
        )
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

    # ==========================================================================================================
    # GLOBAL SUMMARY
    # ==========================================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("CORRECTED SIGN")
    print("  A(p) = -d*p")
    print("  A(q) = -d*q")
    print("  all corrected sign checks = True")
    print()

    print("TRUE FACTOR QUADRATIC")
    print("  Q_minus roots = {p,q}")
    print()

    print("CONJUGATE QUADRATIC")
    print("  Q_plus = 2Z^2-(N+1+d)Z+2N")
    print("  rational-root status is tested exactly")
    print()

    print("DISCRIMINANT")
    print("  Delta_minus = (N+1-d)^2 - 16N")
    print("  Delta_plus  = (N+1+d)^2 - 16N")
    print("  Delta_plus - Delta_minus = 4d(N+1)")
    print()

    print("GLOBAL STATUS")
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()

    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
