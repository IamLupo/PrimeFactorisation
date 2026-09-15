from __future__ import annotations

import math


EXP_NO = 404

print("=" * 110)
print(f"EXPERIMENT {EXP_NO} START")
print("=" * 110)
print()
print("EXACT CONJUGATE DISCRIMINANT / N,K INVARIANT AUDIT")
print()
print("Objective:")
print("  1. Correct the Delta_plus + Delta_minus identity from")
print("     Experiment 403.")
print("  2. Verify:")
print("       Delta_plus + Delta_minus")
print("         = 2*((N+1)^2 + d^2 - 16N)")
print("  3. Verify:")
print("       Delta_plus - Delta_minus = 4*d*(N+1)")
print("  4. Compute the exact product:")
print("       Delta_plus * Delta_minus")
print("  5. Eliminate d^2 using:")
print("       d^2 = -4K - 3N^2 + 6N + 1")
print("  6. Determine whether the product of the two")
print("     discriminants reduces to a simple polynomial in N,K.")
print("  7. Verify the resulting N,K-only identity exactly.")
print("  8. Compare the size of the discriminant product with")
print("     the original quartic coefficients.")
print("  9. Determine whether Delta_plus has any square or")
print("     factor structure across the instances.")
print(" 10. Establish whether the conjugate branch contributes")
print("     any independent arithmetic information.")
print()
print("No resultants are constructed.")
print("No symbolic multivariate factoring is performed.")
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


def poly_eval(coeffs, x: int) -> int:
    out = 0
    for c in coeffs:
        out = out * x + c
    return out


def factor_signature(n: int, small_limit: int = 1000):
    """
    Small-prime factor signature only.
    This deliberately does NOT attempt complete factorization.
    """
    n = abs(int(n))
    out = []

    if n == 0:
        return [("zero", 1)]

    power = 0
    while n % 2 == 0:
        n //= 2
        power += 1

    if power:
        out.append((2, power))

    p = 3
    while p <= small_limit and p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1

        if e:
            out.append((p, e))

        p += 2

    if n > 1:
        out.append(("remaining", n))

    return out


def main():

    global_ok = True
    all_delta_plus_square = True

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
        # TWO DISCRIMINANTS
        # ------------------------------------------------------------------------------------------

        Delta_minus = (N + 1 - d) ** 2 - 16 * N
        Delta_plus = (N + 1 + d) ** 2 - 16 * N

        sqrt_minus, minus_square = exact_sqrt(Delta_minus)
        sqrt_plus, plus_square = exact_sqrt(Delta_plus)

        all_delta_plus_square &= plus_square

        # True known relation:
        expected_minus = 4 * (q - p) ** 2
        minus_identity = Delta_minus == expected_minus

        # ------------------------------------------------------------------------------------------
        # CORRECTED SUM
        # ------------------------------------------------------------------------------------------

        delta_sum = Delta_plus + Delta_minus

        expected_sum = (
            2 * (
                (N + 1) ** 2
                + d * d
                - 16 * N
            )
        )

        sum_identity = delta_sum == expected_sum

        # ------------------------------------------------------------------------------------------
        # DIFFERENCE
        # ------------------------------------------------------------------------------------------

        delta_difference = Delta_plus - Delta_minus
        expected_difference = 4 * d * (N + 1)

        difference_identity = (
            delta_difference == expected_difference
        )

        # ------------------------------------------------------------------------------------------
        # PRODUCT
        # ------------------------------------------------------------------------------------------

        delta_product = Delta_plus * Delta_minus

        # Let:
        #
        # B = (N+1)^2 + d^2 - 16N
        #
        # Then:
        #
        # Delta_plus  = B + 2d(N+1)
        # Delta_minus = B - 2d(N+1)
        #
        # Therefore:
        #
        # Delta_plus * Delta_minus
        #   = B^2 - 4d^2(N+1)^2

        B = (
            (N + 1) ** 2
            + d * d
            - 16 * N
        )

        product_from_d = (
            B * B
            - 4 * d * d * (N + 1) ** 2
        )

        product_identity = delta_product == product_from_d

        # ------------------------------------------------------------------------------------------
        # ELIMINATE d^2
        # ------------------------------------------------------------------------------------------

        d2_from_NK = (
            -4 * K
            - 3 * N * N
            + 6 * N
            + 1
        )

        B_NK = (
            (N + 1) ** 2
            + d2_from_NK
            - 16 * N
        )

        product_NK = (
            B_NK * B_NK
            - 4 * d2_from_NK * (N + 1) ** 2
        )

        product_NK_identity = (
            delta_product == product_NK
        )

        # ------------------------------------------------------------------------------------------
        # DIRECT EXPANDED N,K FORM
        # ------------------------------------------------------------------------------------------
        #
        # B_NK:
        #
        # (N+1)^2 -16N + (-4K -3N^2 +6N +1)
        #
        # = -2N^2 -8N +2 -4K
        #
        # = -2*(N^2 +4N -1 +2K)
        #
        # So:
        #
        # product =
        #   4*(N^2+4N-1+2K)^2
        #   -4*(-4K-3N^2+6N+1)*(N+1)^2
        #
        # Factor 4:
        #
        # product =
        # 4[
        #   (N^2+4N-1+2K)^2
        #   + (4K+3N^2-6N-1)(N+1)^2
        # ]
        #

        compact_product_NK = (
            4 * (
                (N * N + 4 * N - 1 + 2 * K) ** 2
                + (4 * K + 3 * N * N - 6 * N - 1)
                * (N + 1) ** 2
            )
        )

        compact_product_identity = (
            delta_product == compact_product_NK
        )

        # ------------------------------------------------------------------------------------------
        # CONJUGATE QUADRATIC
        # ------------------------------------------------------------------------------------------

        Q_plus = [
            2,
            -(N + 1 + d),
            2 * N,
        ]

        qplus_p = poly_eval(Q_plus, p)
        qplus_q = poly_eval(Q_plus, q)

        qplus_relation_p = qplus_p == -2 * d * p
        qplus_relation_q = qplus_q == -2 * d * q

        # ------------------------------------------------------------------------------------------
        # BASIC INSTANCE STATUS
        # ------------------------------------------------------------------------------------------

        instance_ok = all([
            d_ok,
            minus_identity,
            sum_identity,
            difference_identity,
            product_identity,
            product_NK_identity,
            compact_product_identity,
            qplus_relation_p,
            qplus_relation_q,
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
        print(f"  D = -4K-3N^2+6N+1 = {D}")
        print(f"  D = d^2 : {d_square}")
        print(f"  recovered d = {d}")
        print(f"  true d = {d_true}")
        print(f"  d recovery = {d_ok}")
        print()

        print("DISCRIMINANT MINUS")
        print("  Delta_minus = (N+1-d)^2 - 16N")
        print(f"  Delta_minus = {Delta_minus}")
        print(f"  Delta_minus square = {minus_square}")
        print(
            f"  Delta_minus = 4(q-p)^2 : "
            f"{minus_identity}"
        )
        print()

        print("DISCRIMINANT PLUS")
        print("  Delta_plus = (N+1+d)^2 - 16N")
        print(f"  Delta_plus = {Delta_plus}")
        print(f"  Delta_plus square = {plus_square}")

        if plus_square:
            print(f"  sqrt(Delta_plus) = {sqrt_plus}")
        else:
            print("  Q_plus has no rational roots.")

        print()

        print("CORRECTED DISCRIMINANT SUM")
        print(
            "  Delta_plus + Delta_minus"
        )
        print(
            "    = 2*((N+1)^2 + d^2 - 16N)"
        )
        print(f"  actual = {delta_sum}")
        print(f"  expected = {expected_sum}")
        print(f"  sum identity = {sum_identity}")
        print()

        print("DISCRIMINANT DIFFERENCE")
        print(
            "  Delta_plus - Delta_minus = 4*d*(N+1)"
        )
        print(f"  actual = {delta_difference}")
        print(f"  expected = {expected_difference}")
        print(
            f"  difference identity = "
            f"{difference_identity}"
        )
        print()

        print("DISCRIMINANT PRODUCT")
        print("  Delta_plus * Delta_minus =")
        print(f"    {delta_product}")
        print()
        print("  d-based reconstruction =")
        print(f"    {product_from_d}")
        print(
            f"  d-based product identity = "
            f"{product_identity}"
        )
        print()

        print("N,K-ONLY PRODUCT")
        print(
            "  d^2 = -4K - 3N^2 + 6N + 1"
        )
        print(
            "  Delta_plus * Delta_minus ="
        )
        print(
            "    4*[(N^2+4N-1+2K)^2"
        )
        print(
            "       +(4K+3N^2-6N-1)(N+1)^2]"
        )
        print(f"  N,K expression = {compact_product_NK}")
        print(
            f"  N,K identity = "
            f"{compact_product_identity}"
        )
        print()

        print("CONJUGATE QUADRATIC")
        print(
            "  Q_plus(Z) = "
            "2Z^2-(N+1+d)Z+2N"
        )
        print(f"  Q_plus(p) = {qplus_p}")
        print(f"  -2*d*p = {-2*d*p}")
        print(f"  relation at p = {qplus_relation_p}")
        print(f"  Q_plus(q) = {qplus_q}")
        print(f"  -2*d*q = {-2*d*q}")
        print(f"  relation at q = {qplus_relation_q}")
        print()

        print("SMALL FACTOR SIGNATURE OF Delta_plus")
        print(
            f"  signature = "
            f"{factor_signature(Delta_plus)}"
        )
        print()

        print("SIZE DATA")
        print(f"  N digits = {digits(N)}")
        print(f"  K digits = {digits(K)}")
        print(f"  d digits = {digits(d)}")
        print(f"  Delta_minus digits = {digits(Delta_minus)}")
        print(f"  Delta_plus digits  = {digits(Delta_plus)}")
        print(f"  Delta_product digits = {digits(delta_product)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {instance_ok}")
        print()

    # ==============================================================================================
    # GLOBAL SUMMARY
    # ==============================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("CORRECTED SUM IDENTITY")
    print(
        "  Delta_plus + Delta_minus"
    )
    print(
        "    = 2*((N+1)^2+d^2-16N)"
    )
    print()

    print("DIFFERENCE IDENTITY")
    print(
        "  Delta_plus - Delta_minus"
    )
    print(
        "    = 4*d*(N+1)"
    )
    print()

    print("PRODUCT IDENTITY")
    print(
        "  Delta_plus * Delta_minus"
    )
    print(
        "    = 4*[(N^2+4N-1+2K)^2"
    )
    print(
        "       +(4K+3N^2-6N-1)(N+1)^2]"
    )
    print()

    print("CONJUGATE SQUARE TEST")
    print(
        f"  Delta_plus perfect square in "
        f"all instances = {all_delta_plus_square}"
    )
    print()

    print("GLOBAL STATUS")
    print(
        f"  ALL INSTANCE CHECKS PASS = {global_ok}"
    )
    print()

    print("=" * 110)
    print(f"EXPERIMENT {EXP_NO} FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
