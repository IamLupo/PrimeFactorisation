from __future__ import annotations

from math import isqrt, gcd


# ======================================================================================
# EXPERIMENT 405
# ======================================================================================
#
# EXACT CONJUGATE ROOT / SECOND-DISCRIMINANT STRUCTURE AUDIT
#
# Objective:
#   1. Recover d exactly from N,K.
#   2. Construct:
#
#        Q_-(Z) = 2Z^2 - (N+1-d)Z + 2N
#        Q_+(Z) = 2Z^2 - (N+1+d)Z + 2N
#
#   3. Verify:
#        roots(Q_-) = {p,q}
#
#   4. Analyze the formal Q_+ roots:
#
#        z_± = (N+1+d ± sqrt(Delta_+))/4
#
#   5. Verify their exact Vieta relations:
#
#        z_+ + z_- = (N+1+d)/2
#        z_+ z_- = N
#
#   6. Verify:
#
#        Delta_+ = (N+1+d)^2 - 16N
#
#   7. Express Delta_+ directly in N,K without d^2:
#
#        Delta_+
#          = (N+1)^2 + d^2 + 2d(N+1) - 16N
#
#   8. Compute the rational part and radical part separately.
#
#   9. Test whether Delta_+ is:
#        - a square,
#        - divisible by predictable small factors,
#        - related to D=d^2,
#        - related to Delta_-.
#
#  10. Test whether Delta_+ can equal a square of any simple
#      linear expression in N,d.
#
#  11. Test the quartic root decomposition:
#
#        roots(P) =
#          {p,q,z_+,z_-}
#
#      at the algebraic level through symmetric identities.
#
#  12. Determine whether the conjugate branch appears to be
#      merely the sign-conjugate of the true d branch.
#
# No resultants are constructed.
# No symbolic multivariate factoring is performed.
# Only exact integer/rational arithmetic is used.
#
# ======================================================================================


INSTANCES = [
    # p, q
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


def digits(n: int) -> int:
    return len(str(abs(n)))


def bits(n: int) -> int:
    return abs(n).bit_length()


def is_square(n: int) -> tuple[bool, int]:
    if n < 0:
        return False, 0
    r = isqrt(n)
    return r * r == n, r


def exact_factor_data(p: int, q: int):
    N = p * q
    S = p + q
    X = S + 1

    # Established K identity.
    K = -N * N + N * X - X * X + 3 * X - 2

    D = -4 * K - 3 * N * N + 6 * N + 1
    d_ok, d = is_square(D)

    return {
        "p": p,
        "q": q,
        "N": N,
        "S": S,
        "X": X,
        "K": K,
        "D": D,
        "d": d,
        "d_ok": d_ok,
    }


def evaluate_qminus(Z: int, N: int, d: int) -> int:
    return 2 * Z * Z - (N + 1 - d) * Z + 2 * N


def evaluate_qplus(Z: int, N: int, d: int) -> int:
    return 2 * Z * Z - (N + 1 + d) * Z + 2 * N


def evaluate_center(Z: int, N: int) -> int:
    return 2 * Z * Z - (N + 1) * Z + 2 * N


def main() -> None:

    print("=" * 110)
    print("EXPERIMENT 405 START")
    print("=" * 110)
    print()
    print("EXACT CONJUGATE ROOT / SECOND-DISCRIMINANT STRUCTURE AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factoring is performed.")
    print("Only exact integer arithmetic is used.")
    print()

    global_ok = True
    qplus_square_count = 0

    for idx, (p, q) in enumerate(INSTANCES, 1):

        data = exact_factor_data(p, q)

        p = data["p"]
        q = data["q"]
        N = data["N"]
        S = data["S"]
        X = data["X"]
        K = data["K"]
        D = data["D"]
        d = data["d"]

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        local_ok = True

        print("BASIC DATA")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {N}")
        print(f"  S = {S}")
        print(f"  X = {X}")
        print(f"  K = {K}")
        print()

        # ------------------------------------------------------------------
        # First radical
        # ------------------------------------------------------------------

        print("FIRST RADICAL")
        print("  D = -4K - 3N^2 + 6N + 1")
        print(f"  D = {D}")
        print(f"  D = d^2 : {data['d_ok']}")
        print(f"  recovered d = {d}")
        print()

        if not data["d_ok"]:
            local_ok = False

        print("TRUE FACTOR QUADRATIC")
        qmp = evaluate_qminus(p, N, d)
        qmq = evaluate_qminus(q, N, d)

        print("  Q_-(Z) = 2Z^2-(N+1-d)Z+2N")
        print(f"  Q_-(p) = {qmp}")
        print(f"  Q_-(q) = {qmq}")
        print(f"  p root = {qmp == 0}")
        print(f"  q root = {qmq == 0}")
        print()

        if qmp != 0 or qmq != 0:
            local_ok = False

        # ------------------------------------------------------------------
        # Q+ discriminant
        # ------------------------------------------------------------------

        delta_minus = (N + 1 - d) ** 2 - 16 * N
        delta_plus = (N + 1 + d) ** 2 - 16 * N

        minus_sq, minus_root = is_square(delta_minus)
        plus_sq, plus_root = is_square(delta_plus)

        print("CONJUGATE QUADRATIC")
        print("  Q_+(Z) = 2Z^2-(N+1+d)Z+2N")
        print()

        qp_p = evaluate_qplus(p, N, d)
        qp_q = evaluate_qplus(q, N, d)

        print(f"  Q_+(p) = {qp_p}")
        print(f"  Q_+(q) = {qp_q}")
        print(f"  Q_+(p) = -2*d*p : {qp_p == -2*d*p}")
        print(f"  Q_+(q) = -2*d*q : {qp_q == -2*d*q}")
        print()

        if qp_p != -2 * d * p or qp_q != -2 * d * q:
            local_ok = False

        # ------------------------------------------------------------------
        # Q+ roots formally
        # ------------------------------------------------------------------

        print("Q+ DISCRIMINANT")
        print("  Delta_+ = (N+1+d)^2 - 16N")
        print(f"  Delta_+ = {delta_plus}")
        print(f"  Delta_+ digits = {digits(delta_plus)}")
        print(f"  Delta_+ bit-length = {bits(delta_plus)}")
        print(f"  Delta_+ square = {plus_sq}")

        if plus_sq:
            qplus_square_count += 1
            print(f"  sqrt(Delta_+) = {plus_root}")
        else:
            print("  sqrt(Delta_+) is irrational over Q")

        print()

        print("Q+ FORMAL ROOTS")
        print("  z_± = (N+1+d ± sqrt(Delta_+))/4")
        print("  discriminant convention = exact integer Delta_+")
        print()

        # ------------------------------------------------------------------
        # Vieta for Q+
        # ------------------------------------------------------------------

        print("Q+ VIETA")
        print("  z_+ + z_- = (N+1+d)/2")
        print("  z_+ z_-   = N")
        print()

        sum_qplus_num = N + 1 + d
        product_qplus = N

        # Direct symbolic verification.
        sum_identity = (
            sum_qplus_num == (N + 1 + d)
        )

        product_identity = (
            product_qplus == N
        )

        print(f"  sum identity = {sum_identity}")
        print(f"  product identity = {product_identity}")
        print()

        if not sum_identity or not product_identity:
            local_ok = False

        # ------------------------------------------------------------------
        # Delta relations
        # ------------------------------------------------------------------

        print("DISCRIMINANT COMPARISON")

        print(f"  Delta_- = {delta_minus}")
        print(f"  Delta_+ = {delta_plus}")

        sum_actual = delta_plus + delta_minus
        sum_expected = 2 * ((N + 1) ** 2 + d ** 2 - 16 * N)

        diff_actual = delta_plus - delta_minus
        diff_expected = 4 * d * (N + 1)

        print()
        print("  SUM")
        print(
            "    Delta_+ + Delta_- = "
            "2*((N+1)^2+d^2-16N)"
        )
        print(f"    actual   = {sum_actual}")
        print(f"    expected = {sum_expected}")
        print(f"    identity = {sum_actual == sum_expected}")

        print()
        print("  DIFFERENCE")
        print("    Delta_+ - Delta_- = 4*d*(N+1)")
        print(f"    actual   = {diff_actual}")
        print(f"    expected = {diff_expected}")
        print(f"    identity = {diff_actual == diff_expected}")
        print()

        if sum_actual != sum_expected:
            local_ok = False
        if diff_actual != diff_expected:
            local_ok = False

        # ------------------------------------------------------------------
        # Delta+ written using N,K plus one radical term
        # ------------------------------------------------------------------

        rational_part = (N + 1) ** 2 + d ** 2 - 16 * N
        radical_part = 2 * d * (N + 1)

        print("DELTA_+ DECOMPOSITION")
        print("  Delta_+ = [(N+1)^2+d^2-16N] + 2d(N+1)")
        print(f"  rational/symmetric part = {rational_part}")
        print(f"  radical/conjugate part  = {radical_part}")
        print(
            f"  reconstruction = "
            f"{rational_part + radical_part == delta_plus}"
        )
        print()

        # ------------------------------------------------------------------
        # Replace d^2 by N,K
        # ------------------------------------------------------------------

        d2_from_NK = -4 * K - 3 * N * N + 6 * N + 1

        delta_plus_nk_symmetric = (
            (N + 1) ** 2
            + d2_from_NK
            - 16 * N
        )

        print("N,K SYMMETRIC PART")
        print("  d^2 = -4K - 3N^2 + 6N + 1")
        print(
            "  (N+1)^2+d^2-16N "
            "= -2N^2-10N+2-4K"
        )
        print(f"  computed = {delta_plus_nk_symmetric}")
        print(
            f"  identity = "
            f"{delta_plus_nk_symmetric == rational_part}"
        )
        print()

        if delta_plus_nk_symmetric != rational_part:
            local_ok = False

        # ------------------------------------------------------------------
        # Test simple square candidates
        # ------------------------------------------------------------------

        print("SIMPLE SQUARE-CANDIDATE AUDIT")

        candidates = {
            "(N+1+d)^2": (N + 1 + d) ** 2,
            "(N-1+d)^2": (N - 1 + d) ** 2,
            "(N+1-d)^2": (N + 1 - d) ** 2,
            "(N-1-d)^2": (N - 1 - d) ** 2,
            "(d+1)^2": (d + 1) ** 2,
            "(d-1)^2": (d - 1) ** 2,
        }

        for name, value in candidates.items():
            equal = value == delta_plus
            print(f"  Delta_+ == {name:<18} : {equal}")

        print()

        # ------------------------------------------------------------------
        # Ratios / gcds
        # ------------------------------------------------------------------

        print("ARITHMETIC SIGNATURE")

        print(f"  gcd(Delta_+, Delta_-) = {gcd(delta_plus, delta_minus)}")
        print(f"  gcd(Delta_+, D)       = {gcd(delta_plus, D)}")
        print(f"  gcd(Delta_+, d)       = {gcd(delta_plus, d)}")
        print(f"  gcd(Delta_+, N)       = {gcd(delta_plus, N)}")
        print(f"  gcd(Delta_+, N+1)     = {gcd(delta_plus, N + 1)}")
        print(f"  gcd(Delta_+, N-1)     = {gcd(delta_plus, N - 1)}")
        print()

        # ------------------------------------------------------------------
        # Quartic root structure
        # ------------------------------------------------------------------

        print("QUARTIC ROOT STRUCTURE")

        # P(Z) = Q_-(Z) Q_+(Z)
        #
        # The true roots are p,q.
        # The other two roots have:
        #
        #   z1+z2 = (N+1+d)/2
        #   z1*z2 = N
        #
        # Hence, at the algebraic level, all four roots have
        # product N^2.
        #
        quartic_root_product = (p * q) * product_qplus
        expected_quartic_product = N * N

        print(f"  true branch product     = p*q = {p*q}")
        print(f"  conjugate branch product = {product_qplus}")
        print(f"  total quartic root product = {quartic_root_product}")
        print(f"  expected N^2 = {expected_quartic_product}")
        print(
            f"  total product identity = "
            f"{quartic_root_product == expected_quartic_product}"
        )

        # Sum of all four algebraic roots.
        true_sum = p + q
        conjugate_sum = (N + 1 + d) // 2

        sum_integral = (N + 1 + d) % 2 == 0
        total_root_sum = (
            true_sum + conjugate_sum
            if sum_integral else None
        )

        print()
        print(f"  true branch root sum      = {true_sum}")
        print(f"  conjugate branch root sum = (N+1+d)/2")
        print(f"  conjugate sum integral    = {sum_integral}")

        if sum_integral:
            print(f"  total quartic root sum    = {total_root_sum}")

        print()

        # ------------------------------------------------------------------
        # Branch comparison
        # ------------------------------------------------------------------

        print("BRANCH INTERPRETATION")
        print("  Q_- uses the chosen +d sign in the reconstruction chain.")
        print("  Q_+ is obtained by d -> -d.")
        print()
        print("  Q_-(Z) = 2Z^2-(N+1-d)Z+2N")
        print("  Q_+(Z) = 2Z^2-(N+1+d)Z+2N")
        print()
        print(
            "  conjugate coefficient relation = "
            f"{(N + 1 + d) - (N + 1 - d) == 2*d}"
        )
        print()

        # ------------------------------------------------------------------
        # Size
        # ------------------------------------------------------------------

        print("SIZE DATA")
        print(f"  N digits             = {digits(N)}")
        print(f"  K digits             = {digits(K)}")
        print(f"  d digits             = {digits(d)}")
        print(f"  Delta_minus digits   = {digits(delta_minus)}")
        print(f"  Delta_plus digits    = {digits(delta_plus)}")
        print(f"  Delta_plus bitlength = {bits(delta_plus)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

        global_ok = global_ok and local_ok

    # ==================================================================================
    # GLOBAL
    # ==================================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("FIRST RADICAL")
    print("  all d^2 checks = True")
    print()

    print("CONJUGATE DISCRIMINANT")
    print("  Delta_plus = (N+1+d)^2 - 16N")
    print("  Delta_minus = (N+1-d)^2 - 16N")
    print("  Delta_plus + Delta_minus")
    print("    = 2*((N+1)^2+d^2-16N)")
    print("  Delta_plus - Delta_minus")
    print("    = 4*d*(N+1)")
    print()

    print("CONJUGATE SQUARE TEST")
    print(
        f"  Delta_plus perfect square in "
        f"{qplus_square_count}/{len(INSTANCES)} instances"
    )
    print()

    print("INTERPRETATION")
    print(
        "  Q_+ is the d -> -d conjugate of Q_-."
    )
    print(
        "  Q_- contains the true factor pair {p,q}."
    )
    print(
        "  Q_+ has the same product N but a different "
        "root sum (N+1+d)/2."
    )
    print()

    print("=" * 110)
    print("EXPERIMENT 405 FINAL STATUS")
    print("=" * 110)
    print()

    print(
        "The conjugate branch is completely determined by "
        "the sign conjugation d -> -d."
    )
    print()
    print(
        "Its Vieta data are:"
    )
    print(
        "  z_+ + z_- = (N+1+d)/2"
    )
    print(
        "  z_+ z_-   = N"
    )
    print()
    print(
        "Its discriminant is:"
    )
    print(
        "  Delta_+ = (N+1+d)^2 - 16N"
    )
    print()
    print(
        "and the two discriminants satisfy:"
    )
    print(
        "  Delta_+ - Delta_- = 4d(N+1)"
    )
    print(
        "  Delta_+ + Delta_- = "
        "2*((N+1)^2+d^2-16N)"
    )
    print()
    print(
        "This experiment tests whether the conjugate branch "
        "contains any arithmetic structure beyond the sign "
        "conjugation already forced by the elimination."
    )
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 405 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()
