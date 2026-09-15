# ==============================================================================
# EXPERIMENT 410
# ==============================================================================
#
# EXACT CORRECTED NORMALIZATION / SECONDARY U-DISCRIMINANT STRUCTURE AUDIT
#
# Purpose:
#
#   1. Correct the normalization error from EXPERIMENT 409.
#   2. Verify R(Nx)/N^4 exactly.
#   3. Verify the reciprocal-palindromic quartic directly.
#   4. Verify the U-quadratic:
#
#        N U^2 -(N^2-N+K)U -(N^2-1+2K) = 0
#
#   5. Verify its discriminant:
#
#        Delta_U
#          = (N^2-N+K)^2
#            + 4N(N^2-1+2K)
#
#   6. Verify the alternate simplification:
#
#        Delta_U
#          = (N^2+N+K)^2 + 4N(K-1)
#
#   7. Verify:
#
#        4 Delta_U = (q-p)^2 Delta_+
#
#        16 Delta_U = Delta_- Delta_+
#
#   8. Search for additional low-complexity factor forms of Delta_U.
#   9. Test divisibility / gcd signatures.
#  10. Test whether Delta_U has an independent square structure.
#  11. Reconstruct the original quartic by composition:
#
#        N U^2 - A U + C = 0
#
#        x^2 - Ux + 1 = 0
#
#        W = Nx
#
#   12. Verify the entire chain is exact with no d, S or T appearing
#       in the final quartic or U-quadratic.
#
# No resultants.
# No symbolic multivariate factoring.
# Exact integer arithmetic only.
# ==============================================================================

from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt


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


def exact_sqrt(n: int):
    if n < 0:
        return False, 0
    r = isqrt(n)
    return r * r == n, r


def digits(n: int) -> int:
    return len(str(abs(n)))


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def poly_eval(poly, x):
    total = 0
    for c in reversed(poly):
        total = total * x + c
    return total


def factor_small(n: int, limit: int = 100):
    """
    Extract small prime factors <= limit.
    Returns [(prime, exponent), ...], remainder.
    """
    if n == 0:
        return [], 0

    n = abs(n)
    factors = []

    p = 2
    while p <= limit and p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            factors.append((p, e))

        p = 3 if p == 2 else p + 2

    return factors, n


def main():

    print("=" * 110)
    print("EXPERIMENT 410 START")
    print("=" * 110)
    print()
    print("EXACT CORRECTED NORMALIZATION / SECONDARY U-DISCRIMINANT STRUCTURE AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact integer arithmetic is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        # ------------------------------------------------------------------
        # Base quantities
        # ------------------------------------------------------------------

        N = p * q
        S = p + q

        # The experiment family uses:
        #
        #   X = S+1
        #
        X = S + 1

        # Original K relation from the experiments.
        K = -N * N + N * X - X * X + 3 * X - 2

        D = -4 * K - 3 * N * N + 6 * N + 1
        d_ok, d = exact_sqrt(D)

        T = (N + 1 + d) // 2

        A = N * N - N + K
        C = -N * N + 1 - 2 * K

        # ------------------------------------------------------------------
        # Construct cross-branch quartic directly.
        #
        # z1+z2 = T
        # z1z2 = N
        #
        # (W-pz1)(W-pz2) = W²-pT W+p²N
        # (W-qz1)(W-qz2) = W²-qT W+q²N
        # ------------------------------------------------------------------

        R = poly_mul(
            [p * p * N, -p * T, 1],
            [q * q * N, -q * T, 1],
        )

        # Expected quartic in low -> high order.
        expected_R = [
            N**4,
            -N * N * A,
            N * (-N * N + 2 * N + 1 - 2 * K),
            -A,
            1,
        ]

        quartic_ok = (R == expected_R)

        # ------------------------------------------------------------------
        # Correct normalization:
        #
        # R(Nx)/N^4
        #
        # If R(W)=sum a_i W^i, then:
        #
        #   R(Nx)/N^4
        #      = sum a_i N^i / N^4 * x^i
        #      = sum a_i/N^(4-i) * x^i
        # ------------------------------------------------------------------

        normalized = [
            Fraction(R[i], N ** (4 - i))
            for i in range(5)
        ]

        expected_normalized = [
            Fraction(1, 1),
            Fraction(-A, N),
            Fraction(
                N * (-N * N + 2 * N + 1 - 2 * K),
                N * N,
            ),
            Fraction(-A, N),
            Fraction(1, 1),
        ]

        normalization_ok = normalized == expected_normalized

        # ------------------------------------------------------------------
        # Reciprocal symmetry
        # ------------------------------------------------------------------

        reciprocal_0 = R[0] == N**4 * R[4]
        reciprocal_1 = R[1] == N**2 * R[3]

        # ------------------------------------------------------------------
        # U quadratic
        #
        #   U = x + 1/x
        #
        # gives:
        #
        #   N U² - A U + C = 0
        # ------------------------------------------------------------------

        u_poly = [C, -A, N]

        # ------------------------------------------------------------------
        # U discriminant
        # ------------------------------------------------------------------

        Delta_U = A * A - 4 * N * C

        # Direct expanded form.
        Delta_U_expanded = (
            (N * N - N + K) ** 2
            + 4 * N * (N * N - 1 + 2 * K)
        )

        # New compact form discovered algebraically.
        Delta_U_compact = (
            (N * N + N + K) ** 2
            + 4 * N * (K - 1)
        )

        disc_form_1_ok = Delta_U == Delta_U_expanded
        disc_form_2_ok = Delta_U == Delta_U_compact

        # ------------------------------------------------------------------
        # Relation with the two branch discriminants.
        # ------------------------------------------------------------------

        Delta_minus = (N + 1 - d) ** 2 - 16 * N
        Delta_plus = (N + 1 + d) ** 2 - 16 * N

        gap = q - p

        relation_1 = 4 * Delta_U == gap * gap * Delta_plus
        relation_2 = 16 * Delta_U == Delta_minus * Delta_plus

        # ------------------------------------------------------------------
        # Independent direct expression of Delta_-.
        # Since:
        #
        #   Delta_- = 4(q-p)^2
        #
        # ------------------------------------------------------------------

        delta_minus_ok = Delta_minus == 4 * gap * gap

        # ------------------------------------------------------------------
        # Reconstruct Delta_U from Delta_+ and gap.
        # ------------------------------------------------------------------

        reconstructed_Delta_U = (gap * gap * Delta_plus) // 4
        reconstruction_ok = (
            gap * gap * Delta_plus % 4 == 0
            and reconstructed_Delta_U == Delta_U
        )

        # ------------------------------------------------------------------
        # Perfect-square test.
        # ------------------------------------------------------------------

        du_square, du_root = exact_sqrt(Delta_U)

        # ------------------------------------------------------------------
        # Small factor signature.
        # ------------------------------------------------------------------

        small_factors, remaining = factor_small(Delta_U, limit=100)

        # ------------------------------------------------------------------
        # GCD signatures.
        # ------------------------------------------------------------------

        gcd_N = gcd(Delta_U, N)
        gcd_A = gcd(Delta_U, A)
        gcd_K = gcd(Delta_U, K)
        gcd_d2 = gcd(Delta_U, d * d)
        gcd_C = gcd(Delta_U, abs(C))

        # ------------------------------------------------------------------
        # Special square candidates.
        # ------------------------------------------------------------------

        candidate_squares = {
            "(N^2+N+K)^2": (N * N + N + K) ** 2,
            "(N^2-N+K)^2": (N * N - N + K) ** 2,
            "(N^2+1+K)^2": (N * N + 1 + K) ** 2,
            "(N^2+N-1+K)^2": (N * N + N - 1 + K) ** 2,
            "N^2*d^2": N * N * d * d,
            "N*d^2": N * d * d,
        }

        candidate_hits = {
            label: (Delta_U == value)
            for label, value in candidate_squares.items()
        }

        # ------------------------------------------------------------------
        # Composition reconstruction.
        #
        # From:
        #
        #   N U² - A U + C = 0
        #
        # and U=x+1/x:
        #
        # multiply by x²:
        #
        #   N(x²+1)² - A x(x²+1) + Cx² = 0
        #
        # giving:
        #
        #   N x^4 - A x^3 + (2N+C)x² - Ax + N
        #
        # Since:
        #
        #   2N+C = -N²+2N+1-2K
        #
        # this must match the normalized quartic after multiplying by N.
        # ------------------------------------------------------------------

        composed_x = [
            N,
            -A,
            2 * N + C,
            -A,
            N,
        ]

        target_x = [
            N,
            -A,
            -N * N + 2 * N + 1 - 2 * K,
            -A,
            N,
        ]

        composition_ok = composed_x == target_x

        # ------------------------------------------------------------------
        # Complete local status.
        #
        # NOTE:
        # We deliberately do NOT require d to be present in the final
        # polynomial identities except where branch discriminants are tested.
        # ------------------------------------------------------------------

        local_ok = all([
            d_ok,
            quartic_ok,
            normalization_ok,
            reciprocal_0,
            reciprocal_1,
            disc_form_1_ok,
            disc_form_2_ok,
            relation_1,
            relation_2,
            delta_minus_ok,
            reconstruction_ok,
            composition_ok,
        ])

        global_ok = global_ok and local_ok

        # ------------------------------------------------------------------
        # OUTPUT
        # ------------------------------------------------------------------

        print("=" * 110)
        print(f"INSTANCE {idx}")
        print("=" * 110)
        print()

        print("BASIC DATA")
        print(f"  p = {p}")
        print(f"  q = {q}")
        print(f"  N = {N}")
        print(f"  S = {S}")
        print(f"  T = {T}")
        print(f"  K = {K}")
        print()

        print("FIRST RADICAL")
        print(f"  D = {D}")
        print(f"  D = d^2 : {d_ok}")
        print(f"  d = {d}")
        print()

        print("CROSS-BRANCH QUARTIC")
        print(f"  a4 = {R[4]}")
        print(f"  a3 = {R[3]}")
        print(f"  a2 = {R[2]}")
        print(f"  a1 = {R[1]}")
        print(f"  a0 = {R[0]}")
        print()

        print("RECIPROCAL SYMMETRY")
        print(f"  a0 = N^4*a4 : {reciprocal_0}")
        print(f"  a1 = N^2*a3 : {reciprocal_1}")
        print()

        print("CORRECTED NORMALIZATION")
        print("  R(Nx)/N^4 =")
        print("    x^4")
        print("    -(N^2-N+K)/N * x^3")
        print("    +(-N^2+2N+1-2K)/N^2 * x^2")
        print("    -(N^2-N+K)/N * x")
        print("    +1")
        print()
        print(f"  exact normalization = {normalization_ok}")
        print()

        print("U-SUBSTITUTION")
        print("  U = x + 1/x")
        print()

        print("N,K-ONLY U-QUADRATIC")
        print("  N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
        print()
        print(f"  coefficients = {u_poly}")
        print()

        print("U-DISCRIMINANT")
        print(f"  Delta_U = {Delta_U}")
        print(f"  Delta_U digits = {digits(Delta_U)}")
        print(f"  Delta_U square = {du_square}")
        print()

        print("DISCRIMINANT EXPANSION 1")
        print("  Delta_U = (N^2-N+K)^2 + 4N(N^2-1+2K)")
        print(f"  computed = {Delta_U_expanded}")
        print(f"  identity = {disc_form_1_ok}")
        print()

        print("DISCRIMINANT COMPACT FORM")
        print("  Delta_U = (N^2+N+K)^2 + 4N(K-1)")
        print(f"  computed = {Delta_U_compact}")
        print(f"  identity = {disc_form_2_ok}")
        print()

        print("DOUBLE-DISCRIMINANT RELATIONS")
        print("  Delta_- = (N+1-d)^2 - 16N")
        print("  Delta_+ = (N+1+d)^2 - 16N")
        print(f"  Delta_- = {Delta_minus}")
        print(f"  Delta_+ = {Delta_plus}")
        print()

        print("  4*Delta_U = (q-p)^2*Delta_+")
        print(f"    lhs = {4 * Delta_U}")
        print(f"    rhs = {gap * gap * Delta_plus}")
        print(f"    identity = {relation_1}")
        print()

        print("  16*Delta_U = Delta_-*Delta_+")
        print(f"    lhs = {16 * Delta_U}")
        print(f"    rhs = {Delta_minus * Delta_plus}")
        print(f"    identity = {relation_2}")
        print()

        print("DIRECT DELTA_U RECONSTRUCTION")
        print("  Delta_U = (q-p)^2*Delta_+/4")
        print(f"  reconstructed = {reconstructed_Delta_U}")
        print(f"  exact reconstruction = {reconstruction_ok}")
        print()

        print("SMALL FACTOR SIGNATURE")
        print(f"  small factors = {small_factors}")
        print(f"  remaining = {remaining}")
        print()

        print("GCD SIGNATURE")
        print(f"  gcd(Delta_U,N) = {gcd_N}")
        print(f"  gcd(Delta_U,N^2-N+K) = {gcd_A}")
        print(f"  gcd(Delta_U,K) = {gcd_K}")
        print(f"  gcd(Delta_U,d^2) = {gcd_d2}")
        print(f"  gcd(Delta_U,N^2-1+2K) = {gcd_C}")
        print()

        print("SQUARE-CANDIDATE AUDIT")
        for label, hit in candidate_hits.items():
            print(f"  Delta_U == {label:<25}: {hit}")
        print()

        print("COMPOSITION")
        print("  outer quadratic:")
        print("    N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
        print()
        print("  inner quadratic:")
        print("    x^2-Ux+1=0")
        print()
        print("  W = N*x")
        print()

        print("COMPOSITION COEFFICIENT AUDIT")
        print(f"  composed = {composed_x}")
        print(f"  target   = {target_x}")
        print(f"  exact composition = {composition_ok}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

    # ==========================================================================
    # GLOBAL SUMMARY
    # ==========================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("CORRECTED NORMALIZATION")
    print("  R(Nx)/N^4 is now audited with the required N^i scaling.")
    print()

    print("U-QUADRATIC")
    print("  N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
    print()

    print("U-DISCRIMINANT")
    print("  Delta_U = (N^2-N+K)^2 + 4N(N^2-1+2K)")
    print()
    print("  equivalently:")
    print("  Delta_U = (N^2+N+K)^2 + 4N(K-1)")
    print()

    print("DOUBLE-DISCRIMINANT")
    print("  4 Delta_U = (q-p)^2 Delta_+")
    print("  16 Delta_U = Delta_- Delta_+")
    print()

    print("COMPOSITION")
    print("  N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
    print("  x^2-Ux+1=0")
    print("  W=Nx")
    print()

    print("INTERPRETATION")
    print("  The apparent quartic degree comes from composing two quadratic")
    print("  layers under the reciprocal substitution U=x+x^{-1}.")
    print("  Both quadratic layers are determined entirely by N and K.")
    print()

    print("=" * 110)
    print("EXPERIMENT 410 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 410 FINISHED")
    print()


if __name__ == "__main__":
    main()
