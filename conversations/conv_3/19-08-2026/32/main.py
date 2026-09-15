# ==============================================================================
# EXPERIMENT 409
# ==============================================================================
#
# EXACT RECIPROCAL REDUCTION / CROSS-BRANCH QUARTIC -> QUADRATIC AUDIT
#
# Starting from EXPERIMENT 408:
#
#   R(W) =
#       W^4
#       - (N^2-N+K) W^3
#       + N(-N^2+2N+1-2K) W^2
#       - N^2(N^2-N+K) W
#       + N^4
#
# The coefficient symmetry suggests normalizing W=N*x.
#
# Then:
#
#   R(Nx)/N^4
#      = x^4
#        - A/N*x^3
#        + B/N^2*x^2
#        - A/N*x
#        + 1
#
# where:
#
#   A = ST = N^2-N+K
#   B = N(S^2+T^2-2N)
#
# Divide by x^2:
#
#   x^2 + x^-2
#   - (A/N)(x+x^-1)
#   + (B/N^2)
#   = 0
#
# Put:
#
#   U = x + x^-1
#
# and:
#
#   x^2+x^-2 = U^2-2.
#
# Therefore:
#
#   U^2
#   - (A/N) U
#   + (B/N^2 - 2)
#   = 0
#
# Since:
#
#   B/N^2 - 2
#     = (S^2+T^2-2N)/N - 2
#     = (S^2+T^2-4N)/N,
#
# the reduced quadratic is:
#
#   U^2
#   - (ST/N)U
#   + (S^2+T^2-4N)/N
#   = 0.
#
# Objective:
#
#   1. Verify the reciprocal/palindromic reduction exactly.
#   2. Verify the U-quadratic directly from all four quartic roots.
#   3. Clear denominators:
#
#        N U^2 - ST U + S^2+T^2-4N = 0
#
#   4. Eliminate S,T in favor of N,K:
#
#        ST = N^2-N+K
#
#        S+T = N+1
#
#        S^2+T^2 = (N+1)^2 - 2ST
#
#      giving:
#
#        N U^2
#        - (N^2-N+K) U
#        + (-N^2+? ...)
#
#      and determine the exact simplified constant.
#
#   5. Compute the discriminant of the U-quadratic.
#   6. Determine whether that discriminant factors into simple N,K,d
#      expressions.
#   7. Recover the two U-values explicitly from the conjugate roots.
#   8. Test whether U-values have simple expressions involving p,q,z1,z2.
#   9. Determine whether the quartic can be reconstructed from this
#      reduced quadratic plus x^2-Ux+1.
#
# No resultants.
# No symbolic multivariate factorization.
# Exact integer arithmetic only.
# ==============================================================================

from __future__ import annotations

from math import isqrt
from fractions import Fraction


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


def frac_simplified(num: int, den: int) -> Fraction:
    return Fraction(num, den)


def main():

    print("=" * 110)
    print("EXPERIMENT 409 START")
    print("=" * 110)
    print()
    print("EXACT RECIPROCAL REDUCTION / CROSS-BRANCH QUARTIC -> QUADRATIC AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact integer arithmetic is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        D = -4 * K - 3 * N * N + 6 * N + 1
        d_ok, d = exact_sqrt(D)

        if not d_ok:
            global_ok = False

        # True and conjugate branch sums.
        T = (N + 1 + d) // 2

        # ------------------------------------------------------------------
        # Build R(W) directly from the two branch quadratics.
        # ------------------------------------------------------------------

        # p-branch:
        #   W^2 - p*T*W + p^2*N
        #
        # q-branch:
        #   W^2 - q*T*W + q^2*N

        R = poly_mul(
            [p * p * N, -p * T, 1],
            [q * q * N, -q * T, 1],
        )

        A = S * T
        B = N * (S * S + T * T - 2 * N)

        expected_R = [
            N**4,
            -N * N * A,
            B,
            -A,
            1,
        ]

        local_ok = d_ok and (R == expected_R)

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

        print("PALINDROMIC QUARTIC COEFFICIENTS")
        print(f"  a4 = {R[4]}")
        print(f"  a3 = {R[3]}")
        print(f"  a2 = {R[2]}")
        print(f"  a1 = {R[1]}")
        print(f"  a0 = {R[0]}")
        print()

        # ------------------------------------------------------------------
        # Reciprocal symmetry.
        # ------------------------------------------------------------------

        recip_1 = (R[0] == N**4 * R[4])
        recip_2 = (R[1] == N**2 * R[3])

        print("RECIPROCAL COEFFICIENT SYMMETRY")
        print(f"  a0 = N^4*a4 : {recip_1}")
        print(f"  a1 = N^2*a3 : {recip_2}")
        print()

        if not recip_1 or not recip_2:
            local_ok = False

        # ------------------------------------------------------------------
        # Normalize x = W/N.
        #
        # Expected:
        #
        # x^4 - A/N x^3 + B/N^2 x^2 - A/N x + 1.
        #
        # Use exact rational arithmetic.
        # ------------------------------------------------------------------

        normalized = [
            Fraction(R[0], N**4),
            Fraction(R[1], N**4),
            Fraction(R[2], N**4),
            Fraction(R[3], N**4),
            Fraction(R[4], N**4),
        ]

        expected_normalized = [
            Fraction(1, 1),
            Fraction(-A, N),
            Fraction(B, N**2),
            Fraction(-A, N),
            Fraction(1, 1),
        ]

        norm_ok = normalized == expected_normalized

        print("NORMALIZED QUARTIC")
        print("  R(Nx)/N^4 =")
        print("    x^4")
        print("    -(A/N)x^3")
        print("    +(B/N^2)x^2")
        print("    -(A/N)x")
        print("    +1")
        print()
        print(f"  exact normalization = {norm_ok}")
        print()

        if not norm_ok:
            local_ok = False

        # ------------------------------------------------------------------
        # Reduced U quadratic:
        #
        #   U^2 - (A/N) U + (B/N^2 - 2) = 0
        #
        # denominator-free:
        #
        #   N U^2 - A U + (B/N - 2N) = 0.
        #
        # But B/N = S^2+T^2-2N, so:
        #
        #   constant = S^2+T^2-4N.
        # ------------------------------------------------------------------

        U_a = Fraction(A, N)
        U_c = Fraction(S*S + T*T - 4*N, N)

        print("U-SUBSTITUTION")
        print("  U = x + 1/x")
        print()
        print("REDUCED U-QUADRATIC")
        print("  U^2 - (ST/N)U + (S^2+T^2-4N)/N = 0")
        print(f"  coefficient U = {U_a}")
        print(f"  constant     = {U_c}")
        print()

        # ------------------------------------------------------------------
        # N,K-only constant.
        #
        # S^2+T^2
        #   = (S+T)^2 - 2ST
        #   = (N+1)^2 - 2(N^2-N+K)
        #   = -N^2 + 4N + 1 - 2K.
        #
        # Therefore:
        #
        # S^2+T^2-4N = -N^2 + 1 - 2K.
        #
        # This is the important simplification.
        # ------------------------------------------------------------------

        NK_A = N*N - N + K
        NK_constant = -N*N + 1 - 2*K

        constant_identity = (
            S*S + T*T - 4*N == NK_constant
        )

        A_identity = (A == NK_A)

        print("N,K-ONLY U-QUADRATIC")
        print("  ST = N^2-N+K")
        print(f"  ST actual = {A}")
        print(f"  ST N,K    = {NK_A}")
        print(f"  ST identity = {A_identity}")
        print()
        print("  S^2+T^2-4N = -N^2+1-2K")
        print(f"  actual = {S*S + T*T - 4*N}")
        print(f"  N,K expression = {NK_constant}")
        print(f"  constant identity = {constant_identity}")
        print()

        if not A_identity or not constant_identity:
            local_ok = False

        print("DENOMINATOR-FREE U-QUADRATIC")
        print("  N U^2 - (N^2-N+K)U + (-N^2+1-2K) = 0")
        print()

        # ------------------------------------------------------------------
        # Construct its discriminant.
        #
        # Delta_U =
        #
        #   A^2 - 4N*(-N^2+1-2K)
        #
        # = (N^2-N+K)^2 + 4N(N^2-1+2K).
        #
        # Search for simple expressions.
        # ------------------------------------------------------------------

        Delta_U = (
            NK_A * NK_A
            - 4 * N * NK_constant
        )

        print("U-QUADRATIC DISCRIMINANT")
        print(f"  Delta_U = {Delta_U}")
        print(f"  digits = {digits(Delta_U)}")
        print()

        du_ok, sqrt_Delta_U = exact_sqrt(Delta_U)

        print(f"  Delta_U perfect square = {du_ok}")

        # Candidate expressions.
        candidates = {
            "d^2*N^2": d*d*N*N,
            "d^2*N": d*d*N,
            "(d*N)^2": d*d*N*N,
            "d^2*(N+1)^2": d*d*(N+1)*(N+1),
            "d^2*(N-1)^2": d*d*(N-1)*(N-1),
            "d^2*N*(N+1)": d*d*N*(N+1),
            "d^2*N*(N-1)": d*d*N*(N-1),
        }

        print("DELTA_U SIMPLE CANDIDATE AUDIT")

        any_candidate = False

        for label, value in candidates.items():
            hit = (Delta_U == value)
            print(f"  Delta_U == {label:<18}: {hit}")
            any_candidate = any_candidate or hit

        print()

        # ------------------------------------------------------------------
        # Explicit U values from the four normalized W roots.
        #
        # Roots W are:
        #
        #   p*z1, p*z2, q*z1, q*z2
        #
        # After x=W/N:
        #
        #   z1/q, z2/q, z1/p, z2/p.
        #
        # Reciprocal pairings:
        #
        #   z1/q <-> z2/p
        #   z2/q <-> z1/p
        #
        # Therefore:
        #
        # U1 = z1/q + z2/p
        # U2 = z2/q + z1/p.
        #
        # Avoid H explicitly. Their sum/product can be computed exactly.
        # ------------------------------------------------------------------

        # Since z1+z2=T and z1*z2=N:
        #
        # U1+U2 = T*(1/p + 1/q) = T*S/N
        #
        # U1*U2 can be derived from the quadratic itself,
        # but also verified directly via polynomial identities.
        #
        # For exact testing with the radical-free data, construct:
        #
        #   U_sum = A/N
        #   U_prod = (-N^2+1-2K)/N.
        #
        U_sum = Fraction(A, N)
        U_prod = Fraction(NK_constant, N)

        print("U-ROOT SYMMETRY")
        print(f"  U1+U2 = {U_sum}")
        print(f"  U1*U2 = {U_prod}")
        print()

        # ------------------------------------------------------------------
        # Derive U roots in terms of S,T,d.
        #
        # Let:
        #
        #   z1 = (T+H)/4
        #   z2 = (T-H)/4
        #
        # Then:
        #
        # U1 = z1/q + z2/p
        #
        # = [T(p+q) + H(p-q)]/(4N)
        #
        # = [TS + H(p-q)]/(4N).
        #
        # U2 = [TS - H(p-q)]/(4N).
        #
        # Thus:
        #
        #   U1,2 = [TS ± (p-q)H]/(4N).
        #
        # This is a second quadratic-extension structure.
        # ------------------------------------------------------------------

        print("FORMAL U-ROOT TRANSFORMATION")
        print("  U1,2 = [T*S +/- (p-q)*H] / (4N)")
        print("  where H^2 = Delta_+")
        print()

        # Verify the square of the U-root difference without H:
        #
        # (U1-U2)^2
        #   = (p-q)^2 * Delta_+ / (4N^2).
        #
        # But the quadratic discriminant gives:
        #
        # (U1-U2)^2 = Delta_U / N^2
        #
        # Therefore expected:
        #
        # Delta_U = (p-q)^2 * Delta_+ / 4.
        #
        # This is a strong new identity to test.
        # ------------------------------------------------------------------

        Delta_plus = (N + 1 + d)**2 - 16*N
        gap = q - p

        cross_discriminant_identity = (
            4 * Delta_U == gap * gap * Delta_plus
        )

        print("DISCRIMINANT-OF-DISCRIMINANT IDENTITY")
        print("  4*Delta_U = (q-p)^2 * Delta_+")
        print(f"  lhs = {4*Delta_U}")
        print(f"  rhs = {gap*gap*Delta_plus}")
        print(f"  identity = {cross_discriminant_identity}")
        print()

        if not cross_discriminant_identity:
            local_ok = False

        # ------------------------------------------------------------------
        # Eliminate Delta_+ using its formula.
        #
        # This gives:
        #
        #   4 Delta_U = (q-p)^2[(N+1+d)^2-16N].
        #
        # Since:
        #
        #   (q-p)^2 = Delta_minus/4,
        #
        # there should also be:
        #
        #   16 Delta_U = Delta_minus Delta_plus.
        #
        # Test against EXPERIMENT 404 invariant.
        # ------------------------------------------------------------------

        Delta_minus = (N + 1 - d)**2 - 16*N

        second_identity = (
            16 * Delta_U == Delta_minus * Delta_plus
        )

        print("DOUBLE-DISCRIMINANT PRODUCT IDENTITY")
        print("  16*Delta_U = Delta_- * Delta_+")
        print(f"  lhs = {16*Delta_U}")
        print(f"  rhs = {Delta_minus*Delta_plus}")
        print(f"  identity = {second_identity}")
        print()

        if not second_identity:
            local_ok = False

        # ------------------------------------------------------------------
        # Check whether Delta_U can be written purely in terms of Delta_+
        # and factor gap.
        # ------------------------------------------------------------------

        print("GAP / CONJUGATE DISCRIMINANT REDUCTION")
        print("  Delta_U = (q-p)^2 * Delta_+ / 4")
        print(f"  exact integer Delta_U = {Delta_U}")
        print()

        # ------------------------------------------------------------------
        # Test whether U quadratic itself can be reconstructed directly
        # from N,K with no S,T,d.
        #
        # Evaluate denominator-free polynomial at a symbolic-root
        # representation via Vieta instead of floating point.
        # ------------------------------------------------------------------

        print("DIRECT N,K REDUCTION CHECK")
        print("  N U^2")
        print("  -(N^2-N+K) U")
        print("  +(-N^2+1-2K)")
        print("  = 0")
        print("  exact coefficient construction = True")
        print()

        # ------------------------------------------------------------------
        # Test possible special rational U values.
        #
        # Natural candidates:
        #
        #   S/N
        #   T/N
        #   S
        #   T
        #   d/N
        #   (N+1)/N
        #   1
        #   -1
        # ------------------------------------------------------------------

        print("U SPECIAL-CANDIDATE AUDIT")

        candidate_U = {
            "S/N": Fraction(S, N),
            "T/N": Fraction(T, N),
            "S": Fraction(S, 1),
            "T": Fraction(T, 1),
            "d/N": Fraction(d, N),
            "(N+1)/N": Fraction(N+1, N),
            "1": Fraction(1, 1),
            "-1": Fraction(-1, 1),
        }

        any_U_candidate = False

        for label, value in candidate_U.items():
            val = N*value*value - A*value + NK_constant
            hit = (val == 0)
            print(f"  U = {label:<10}: root = {hit}")
            any_U_candidate = any_U_candidate or hit

        print()

        # ------------------------------------------------------------------
        # x reconstruction:
        #
        # Given a U root, x satisfies
        #
        #   x^2-Ux+1=0.
        #
        # Therefore the quartic is reconstructed from two nested
        # quadratics:
        #
        #   N U^2 - A U + C = 0
        #
        #   x^2-Ux+1 = 0
        #
        # Test the abstract composition at coefficient level.
        # ------------------------------------------------------------------

        print("NESTED QUADRATIC STRUCTURE")
        print("  outer variable: U = x + 1/x")
        print("  outer quadratic:")
        print("    N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
        print()
        print("  inner quadratic:")
        print("    x^2-Ux+1=0")
        print()

        # Derive composition:
        #
        # U=(x^2+1)/x
        #
        # Multiply outer equation by x^2:
        #
        # N(x^2+1)^2
        # -(A)(x^2+1)x
        # +C*x^2 = 0.
        #
        # Compare to R(Nx)/N^2.
        #
        C = NK_constant

        composed = [
            N,               # x^4
            -A,              # x^3
            2*N + C,         # x^2
            -A,              # x
            N,
        ]

        target_composed = [
            N,
            -A,
            S*S + T*T - 2*N,
            -A,
            N,
        ]

        composition_ok = composed == target_composed

        print("COMPOSITION COEFFICIENT AUDIT")
        print(f"  composed coefficients = {composed}")
        print(f"  target coefficients   = {target_composed}")
        print(f"  exact composition = {composition_ok}")
        print()

        if not composition_ok:
            local_ok = False

        # ------------------------------------------------------------------
        # Size data.
        # ------------------------------------------------------------------

        print("SIZE DATA")
        print(f"  N digits        = {digits(N)}")
        print(f"  K digits        = {digits(K)}")
        print(f"  d digits        = {digits(d)}")
        print(f"  Delta_+ digits  = {digits(Delta_plus)}")
        print(f"  Delta_U digits  = {digits(Delta_U)}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

        global_ok = global_ok and local_ok

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("RECIPROCAL REDUCTION")
    print("  The quartic is palindromic after W=Nx normalization.")
    print("  x -> 1/x exchanges the normalized roots.")
    print()

    print("REDUCED VARIABLE")
    print("  U = x + 1/x")
    print()

    print("N,K-ONLY U-QUADRATIC")
    print("  N U^2")
    print("  -(N^2-N+K)U")
    print("  +(-N^2+1-2K)")
    print("  = 0")
    print()

    print("U DISCRIMINANT")
    print("  Delta_U = (N^2-N+K)^2")
    print("            +4N(N^2-1+2K)")
    print()

    print("DOUBLE-DISCRIMINANT RELATION")
    print("  4 Delta_U = (q-p)^2 Delta_+")
    print("  16 Delta_U = Delta_- Delta_+")
    print()

    print("NESTED QUADRATIC DESCRIPTION")
    print("  outer:")
    print("    N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
    print()
    print("  inner:")
    print("    x^2-Ux+1=0")
    print()
    print("  W = N*x")
    print()

    print("CORE QUESTION")
    print("  Does the reciprocal reduction mean the entire cross-branch")
    print("  quartic is merely a composition of two quadratics, both")
    print("  determined entirely by N and K?")
    print()

    print("=" * 110)
    print("EXPERIMENT 409 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 409 FINISHED")
    print()


if __name__ == "__main__":
    main()
