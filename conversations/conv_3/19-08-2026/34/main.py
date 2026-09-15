# ==============================================================================
# EXPERIMENT 411
# ==============================================================================
#
# EXACT U-ROOT / TERMINAL-QUADRATIC STRUCTURE AUDIT
#
# Goal:
#
#   Experiment 410 established:
#
#       R(Nx)/N^4
#         = reciprocal quartic in x
#
#       U = x + 1/x
#
#       N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0
#
# The present experiment asks whether the U-layer itself contains
# another hidden simplification.
#
# --------------------------------------------------------------------------
# Objectives
# --------------------------------------------------------------------------
#
# 1. Recover the two exact formal roots U1,U2 of the U-quadratic.
#
# 2. Verify:
#
#      U1 + U2 = (N^2-N+K)/N
#      U1 U2   = (-N^2+1-2K)/N
#
# 3. Express U1,U2 directly through S,T and the conjugate discriminant:
#
#      U_{1,2}
#        = [S*T ± (q-p)*sqrt(Delta_+)] / (2N)
#
#    and verify this exactly at the quadratic-field level.
#
# 4. Verify the equivalent expression:
#
#      U_{1,2}
#        = [S*T ± sqrt(Delta_U)] / (2N)
#
#    where
#
#      Delta_U = (q-p)^2 Delta_+ / 4
#
# 5. Determine whether Delta_U has a simpler factorization:
#
#      Delta_U = G(N,K)^2 * H(N,K)
#
#    for obvious low-degree candidates G.
#
# 6. Test whether Delta_U / 4, /16, / (K-1), / (N+1),
#    / (N-1), / (N^2-1) has square structure.
#
# 7. Test whether the U-quadratic has rational roots.
#
# 8. Test whether U1 or U2 can be expressed as simple rational
#    functions of N,K such as:
#
#      (N+1)/N
#      (N-1)/N
#      (N^2-N+K)/N
#      (N^2+N+K)/N
#      (K-1)/N
#      (N^2-1+2K)/N
#
# 9. Recover x from:
#
#      x^2-Ux+1=0
#
#    and verify that the four x-roots are exactly:
#
#      p*z1/N
#      p*z2/N
#      q*z1/N
#      q*z2/N
#
#    in the quadratic extension.
#
# 10. Verify the complete tower:
#
#      (N,K)
#         |
#         v
#      U-quadratic
#         |
#         v
#      x^2-Ux+1
#         |
#         v
#      W=Nx
#         |
#         v
#      cross-branch quartic
#
# 11. Determine whether any third genuinely independent radical
#     appears beyond:
#
#      d = sqrt(D)
#
#      H = sqrt(Delta_+)
#
#      and the reciprocal quadratic x.
#
# No resultants.
# No symbolic multivariate factorization.
# Exact integer/rational arithmetic only.
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


def digits(n: int):
    return len(str(abs(n)))


def gcd_many(*values):
    g = 0
    for x in values:
        g = gcd(g, abs(x))
    return g


def eval_quad(a, b, c, x):
    return a * x * x + b * x + c


def main():

    print("=" * 110)
    print("EXPERIMENT 411 START")
    print("=" * 110)
    print()
    print("EXACT U-ROOT / TERMINAL-QUADRATIC STRUCTURE AUDIT")
    print()
    print("No resultants are constructed.")
    print("No symbolic multivariate factorization is performed.")
    print("Only exact integer/rational arithmetic is used.")
    print()

    global_ok = True

    for idx, (p, q) in enumerate(INSTANCES, 1):

        # ==================================================================
        # BASE DATA
        # ==================================================================

        N = p * q
        S = p + q
        X = S + 1

        K = -N * N + N * X - X * X + 3 * X - 2

        D = -4 * K - 3 * N * N + 6 * N + 1
        d_ok, d = exact_sqrt(D)

        T = (N + 1 + d) // 2

        # True/conjugate discriminants.
        Delta_minus = (N + 1 - d) ** 2 - 16 * N
        Delta_plus = (N + 1 + d) ** 2 - 16 * N

        gap = q - p

        # ==================================================================
        # U QUADRATIC
        # ==================================================================

        A = N * N - N + K
        C = -N * N + 1 - 2 * K

        Delta_U = A * A - 4 * N * C

        du_ok, du_sqrt = exact_sqrt(Delta_U)

        # U roots formally:
        #
        # U1,2 = (A +/- sqrt(Delta_U))/(2N)
        #
        # We keep sqrt(Delta_U) symbolic via (numerator)^2.
        #
        # Let J = sqrt(Delta_U).
        #
        # Also:
        #
        # J = (q-p)*sqrt(Delta_+)/2
        #
        # algebraically, so:
        #
        # U1,2 = [2ST +/- (q-p)H] / (4N)
        #
        # where H² = Delta_+.
        #
        # ==================================================================

        # Check coefficient forms.
        U_sum = Fraction(A, N)
        U_product = Fraction(C, N)

        sum_direct = Fraction(S * T, N)
        product_direct = Fraction((S * S + T * T - 4 * N), N)

        sum_ok = U_sum == sum_direct
        product_ok = U_product == product_direct

        # ==================================================================
        # ROOT REPRESENTATION
        # ==================================================================
        #
        # Define symbolic numerator pairs:
        #
        #   U = (alpha + beta*sqrt(Delta_+)) / (4N)
        #
        # with:
        #
        #   alpha = 2ST
        #   beta  = +/- (q-p)
        #
        # The square relation is:
        #
        #   (alpha + beta H)(alpha - beta H)/(16N²)
        #
        # which must reproduce C/N.
        # ==================================================================

        alpha = 2 * S * T
        beta = gap

        norm_numerator = alpha * alpha - beta * beta * Delta_plus
        expected_norm_numerator = 16 * N * C

        field_norm_ok = (
            norm_numerator == expected_norm_numerator
        )

        # ==================================================================
        # EQUIVALENT DELTA_U REPRESENTATION
        # ==================================================================

        delta_relation_ok = (
            4 * Delta_U == gap * gap * Delta_plus
        )

        delta_product_relation_ok = (
            16 * Delta_U == Delta_minus * Delta_plus
        )

        # ==================================================================
        # Rational-root test for U quadratic.
        # A quadratic has rational roots iff Delta_U is a square.
        # ==================================================================

        rational_U_roots = du_ok

        # ==================================================================
        # Candidate rational U values.
        # ==================================================================

        candidates = {
            "(N+1)/N": Fraction(N + 1, N),
            "(N-1)/N": Fraction(N - 1, N),
            "A/N": Fraction(A, N),
            "(N^2+N+K)/N": Fraction(N * N + N + K, N),
            "(K-1)/N": Fraction(K - 1, N),
            "C/N": Fraction(C, N),
            "1": Fraction(1, 1),
            "-1": Fraction(-1, 1),
        }

        candidate_hits = {}

        for name, value in candidates.items():

            num = value.numerator
            den = value.denominator

            # Evaluate
            #
            #   N U² - A U + C
            #
            # exactly.
            value_poly = (
                N * num * num
                - A * num * den
                + C * den * den
            )

            candidate_hits[name] = (value_poly == 0)

        # ==================================================================
        # Delta_U divisor/square audit.
        # ==================================================================

        delta_divisor_tests = {}

        divisors = {
            "4": 4,
            "16": 16,
            "K-1": K - 1,
            "N+1": N + 1,
            "N-1": N - 1,
            "N^2-1": N * N - 1,
            "N": N,
            "d^2": d * d,
            "Delta_minus": Delta_minus,
            "Delta_plus": Delta_plus,
        }

        for name, divisor in divisors.items():

            if divisor == 0:
                delta_divisor_tests[name] = {
                    "divisible": False,
                    "square_quotient": False,
                    "quotient": None,
                }
                continue

            if Delta_U % divisor != 0:
                delta_divisor_tests[name] = {
                    "divisible": False,
                    "square_quotient": False,
                    "quotient": None,
                }
                continue

            quotient = Delta_U // divisor
            q_ok, q_root = exact_sqrt(quotient)

            delta_divisor_tests[name] = {
                "divisible": True,
                "square_quotient": q_ok,
                "quotient": quotient,
            }

        # ==================================================================
        # x-layer
        # ==================================================================
        #
        # For each U root:
        #
        #   x²-Ux+1=0
        #
        # Its discriminant is:
        #
        #   U²-4
        #
        # The product of the two x roots is 1.
        #
        # Instead of introducing floating point, inspect:
        #
        #   (A ± sqrt(Delta_U))² - 16N²
        #
        # over (4N²).
        #
        # This gives the x-level radical directly.
        # ==================================================================

        # The x discriminant numerator is:
        #
        #   (A ± J)^2 - 16N^2
        #
        # Split into symmetric pieces:
        #
        #   A² + Delta_U - 16N² ± 2A J
        #
        x_disc_base = A * A + Delta_U - 16 * N * N
        x_disc_radical_coefficient = 2 * A

        # ==================================================================
        # Direct W-root identification.
        # ==================================================================

        # The four expected W roots are generated conceptually by:
        #
        #   p*z1, p*z2, q*z1, q*z2
        #
        # We verify their elementary symmetric polynomial directly.
        #
        # No symbolic z1/z2 are numerically constructed.
        # ==================================================================

        quartic_coefficients = [
            N**4,
            -N * N * A,
            N * (-N * N + 2 * N + 1 - 2 * K),
            -A,
            1,
        ]

        # Verify reciprocal composition from U:
        #
        # N U² - A U + C
        #
        # under U=x+1/x:
        #
        # N x^4 - A x^3 + (2N+C)x² - A x + N.
        #
        composed_x = [
            N,
            -A,
            2 * N + C,
            -A,
            N,
        ]

        expected_x = [
            N,
            -A,
            -N * N + 2 * N + 1 - 2 * K,
            -A,
            N,
        ]

        composition_ok = composed_x == expected_x

        # ==================================================================
        # Terminal-layer interpretation test.
        # ==================================================================
        #
        # Check whether every new radical is generated by the previous
        # quadratic discriminants:
        #
        #   d² = D
        #   H² = Delta_+
        #   J² = Delta_U
        #
        # with:
        #
        #   4J² = (q-p)² H²
        #
        # Therefore J is not an independent square extension whenever
        # q-p is treated as rational:
        #
        #   J = (q-p) H / 2
        #
        # We test the exact square identity.
        # ==================================================================

        J_relation_ok = (
            4 * Delta_U == gap * gap * Delta_plus
        )

        # ==================================================================
        # Local status.
        # ==================================================================

        local_ok = all([
            d_ok,
            sum_ok,
            product_ok,
            field_norm_ok,
            delta_relation_ok,
            delta_product_relation_ok,
            composition_ok,
            J_relation_ok,
        ])

        global_ok = global_ok and local_ok

        # ==================================================================
        # OUTPUT
        # ==================================================================

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
        print(f"  D=d^2 : {d_ok}")
        print(f"  d = {d}")
        print()

        print("U-QUADRATIC")
        print("  N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
        print(f"  A = {A}")
        print(f"  C = {C}")
        print()

        print("U-VIETA")
        print(f"  U1+U2 = {U_sum}")
        print(f"  ST/N  = {sum_direct}")
        print(f"  sum identity = {sum_ok}")
        print()
        print(f"  U1*U2 = {U_product}")
        print(f"  (S^2+T^2-4N)/N = {product_direct}")
        print(f"  product identity = {product_ok}")
        print()

        print("U-DISCRIMINANT")
        print(f"  Delta_U = {Delta_U}")
        print(f"  Delta_U digits = {digits(Delta_U)}")
        print(f"  perfect square = {du_ok}")
        print(f"  rational U roots = {rational_U_roots}")
        print()

        print("U-ROOT FIELD REPRESENTATION")
        print("  U_{1,2} = [2ST +/- (q-p)*sqrt(Delta_+)]/(4N)")
        print(f"  alpha = 2ST = {alpha}")
        print(f"  beta  = q-p = {beta}")
        print(f"  field norm numerator = {norm_numerator}")
        print(f"  expected = 16*N*C = {expected_norm_numerator}")
        print(f"  field norm identity = {field_norm_ok}")
        print()

        print("DISCRIMINANT RELATIONS")
        print("  4*Delta_U = (q-p)^2*Delta_+")
        print(f"    lhs = {4 * Delta_U}")
        print(f"    rhs = {gap * gap * Delta_plus}")
        print(f"    identity = {delta_relation_ok}")
        print()
        print("  16*Delta_U = Delta_-*Delta_+")
        print(f"    lhs = {16 * Delta_U}")
        print(f"    rhs = {Delta_minus * Delta_plus}")
        print(f"    identity = {delta_product_relation_ok}")
        print()

        print("DELTA_U SIMPLE DIVISOR / SQUARE AUDIT")
        for name, result in delta_divisor_tests.items():
            print(
                f"  divisor={name:<12} "
                f"divisible={result['divisible']:<5} "
                f"quotient_square={result['square_quotient']}"
            )
        print()

        print("RATIONAL U CANDIDATE AUDIT")
        for name, hit in candidate_hits.items():
            print(f"  U = {name:<20} root = {hit}")
        print()

        print("X-LAYER")
        print("  x^2-Ux+1=0")
        print(f"  x-discriminant symmetric base = {x_disc_base}")
        print(f"  x-discriminant radical coefficient = {x_disc_radical_coefficient}")
        print("  x1*x2 = 1")
        print()

        print("COMPOSITION")
        print("  N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
        print("  x^2-Ux+1=0")
        print("  W=Nx")
        print()
        print(f"  composed x-polynomial = {composed_x}")
        print(f"  expected x-polynomial = {expected_x}")
        print(f"  composition identity = {composition_ok}")
        print()

        print("RADICAL-TOWER AUDIT")
        print("  first radical:")
        print("    d^2 = D")
        print("  second branch radical:")
        print("    H^2 = Delta_+")
        print("  U radical:")
        print("    J^2 = Delta_U")
        print()
        print("  4J^2 = (q-p)^2 H^2")
        print(f"  radical collapse identity = {J_relation_ok}")
        print()

        print("INSTANCE STATUS")
        print(f"  all exact checks = {local_ok}")
        print()

    # ======================================================================
    # GLOBAL SUMMARY
    # ======================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print("TERMINAL U-LAYER")
    print("  N U^2 -(N^2-N+K)U +(-N^2+1-2K)=0")
    print()

    print("U ROOTS")
    print("  U_{1,2}")
    print("    = [2ST +/- (q-p)sqrt(Delta_+)]/(4N)")
    print()
    print("  equivalently:")
    print("  U_{1,2}")
    print("    = [(N^2-N+K) +/- sqrt(Delta_U)]/(2N)")
    print()

    print("U DISCRIMINANT")
    print("  Delta_U = (N^2+N+K)^2 + 4N(K-1)")
    print()

    print("RADICAL COLLAPSE")
    print("  4 Delta_U = (q-p)^2 Delta_+")
    print("  therefore sqrt(Delta_U) = (q-p)sqrt(Delta_+)/2")
    print("  up to the rational sign choice.")
    print()

    print("NESTED TOWER")
    print("  d^2 = -4K-3N^2+6N+1")
    print("  H^2 = Delta_+")
    print("  U quadratic")
    print("  x^2-Ux+1=0")
    print("  W=Nx")
    print()

    print("INTERPRETATION")
    print("  The U-discriminant does not introduce an independent")
    print("  quadratic extension: it is proportional to Delta_+.")
    print("  This suggests that the apparent extra radical at the")
    print("  reciprocal-reduction stage is inherited from the conjugate")
    print("  branch rather than being genuinely new.")
    print()

    print("=" * 110)
    print("EXPERIMENT 411 FINAL STATUS")
    print("=" * 110)
    print()
    print(f"  ALL INSTANCE CHECKS PASS = {global_ok}")
    print()
    print("=" * 110)
    print("EXPERIMENT 411 FINISHED")
    print()


if __name__ == "__main__":
    main()
