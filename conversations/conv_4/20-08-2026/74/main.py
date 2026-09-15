#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 551
# 2020 x,y COORDINATES -> MODERN KAPPA CROSSWALK
#
# Main point:
#
#   U = 2y
#   V = 2x - 3
#
#   A0 = U - V = 2y - 2x + 3
#   B0 = U + V = 2y + 2x - 3
#
#   N     = (U^2 - V^2)/r^2
#   S     = 2U/r
#   Delta = 4V^2/r^2
#
# Hence:
#
#   S = 4y/r
#
# The purpose is to rewrite the already-known KAPPA quantities directly
# in the historical x,y coordinate system.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 551 START")
print("=" * 78)
print("2020 x,y COORDINATES -> MODERN KAPPA CROSSWALK")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

x, y, r, z = sp.symbols(
    "x y r z",
    nonzero=True,
)

failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

def simplify_exact(expr):
    return sp.factor(
        sp.cancel(
            sp.together(
                sp.expand(
                    sp.simplify(expr)
                )
            )
        )
    )


def certificate(label, actual, expected=0):
    global failures

    difference = simplify_exact(actual - expected)
    passed = difference == 0

    print(f"  {label}")
    print("    difference =")
    sp.pprint(difference)
    print(f"    PASS = {passed}")
    print()

    if not passed:
        failures += 1

    return passed


def show(label, expr):
    print(f"  {label} =")
    sp.pprint(simplify_exact(expr))
    print()


# ==============================================================================
# [1] HISTORICAL COORDINATES
# ==============================================================================

print("[1] HISTORICAL COORDINATES")
print("-" * 78)

U = 2 * y
V = 2 * x - 3

A0 = simplify_exact(U - V)
B0 = simplify_exact(U + V)

show("U", U)
show("V", V)
show("A0", A0)
show("B0", B0)

certificate(
    "A0 = 2*y - 2*x + 3",
    A0,
    2 * y - 2 * x + 3,
)

certificate(
    "B0 = 2*y + 2*x - 3",
    B0,
    2 * y + 2 * x - 3,
)


# ==============================================================================
# [2] HISTORICAL CONIC
# ==============================================================================

print("[2] HISTORICAL CONIC")
print("-" * 78)

v_xy = y**2 - x**2 + 3*x - 2

K_xy = simplify_exact(
    4 * v_xy - 1
)

show("v(x,y)", v_xy)
show("4v-1", K_xy)

certificate(
    "4v-1 = U^2-V^2",
    K_xy,
    U**2 - V**2,
)

certificate(
    "4v-1 = A0*B0",
    K_xy,
    A0 * B0,
)


# ==============================================================================
# [3] MODERN KAPPA COORDINATES
# ==============================================================================

print("[3] MODERN KAPPA COORDINATES IN x,y")
print("-" * 78)

N_xy = simplify_exact(
    (U**2 - V**2) / r**2
)

S_xy = simplify_exact(
    2 * U / r
)

Delta_xy = simplify_exact(
    4 * V**2 / r**2
)

M1_xy = simplify_exact(
    (U + r)**2 / r**2 - V**2 / r**2
)

show("N(x,y)", N_xy)
show("S(x,y)", S_xy)
show("Delta(x,y)", Delta_xy)
show("M1(x,y)", M1_xy)

certificate(
    "N",
    N_xy,
    (U**2 - V**2) / r**2,
)

certificate(
    "S",
    S_xy,
    4 * y / r,
)

certificate(
    "Delta",
    Delta_xy,
    4 * (2*x - 3)**2 / r**2,
)

certificate(
    "M1",
    M1_xy,
    ((U + r)**2 - V**2) / r**2,
)


# ==============================================================================
# [4] HISTORICAL r=1 SPECIALIZATION
# ==============================================================================

print("[4] HISTORICAL r=1 SPECIALIZATION")
print("-" * 78)

N1 = simplify_exact(N_xy.subs(r, 1))
S1 = simplify_exact(S_xy.subs(r, 1))
Delta1 = simplify_exact(Delta_xy.subs(r, 1))
M11 = simplify_exact(M1_xy.subs(r, 1))

show("N", N1)
show("S", S1)
show("Delta", Delta1)
show("M1", M11)

certificate(
    "N = 4y^2-(2x-3)^2",
    N1,
    4*y**2 - (2*x - 3)**2,
)

certificate(
    "S = 4y",
    S1,
    4*y,
)

certificate(
    "Delta = 4(2x-3)^2",
    Delta1,
    4 * (2*x - 3)**2,
)

certificate(
    "N = 4v-1",
    N1,
    4 * v_xy - 1,
)


# ==============================================================================
# [5] BASIC KAPPA LAYERS
# ==============================================================================

print("[5] KAPPA HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = simplify_exact(
    6 * N_xy - S_xy**2 + S_xy
)

F3 = simplify_exact(
    (S_xy + 1) * F2
)

show("F2(x,y)", F2)
show("F3(x,y)", F3)

certificate(
    "F3=(S+1)F2",
    F3,
    (S_xy + 1) * F2,
)


# ==============================================================================
# [6] F2 IN HISTORICAL COORDINATES
# ==============================================================================

print("[6] F2 DIRECT x,y FORM")
print("-" * 78)

F2_1 = simplify_exact(
    F2.subs(r, 1)
)

show("F2(r=1)", F2_1)

expected_F2_xy = simplify_exact(
    -24*x**2 + 72*x + 4*y - 54
)

certificate(
    "F2 = -24*x^2 + 72*x + 4*y - 54",
    F2_1,
    expected_F2_xy,
)


# ==============================================================================
# [7] F3 DIRECT HISTORICAL FORM
# ==============================================================================

print("[7] F3 DIRECT x,y FORM")
print("-" * 78)

F3_1 = simplify_exact(
    F3.subs(r, 1)
)

show("F3(r=1)", F3_1)

certificate(
    "F3=(4y+1)F2",
    F3_1,
    (4*y + 1) * F2_1,
)


# ==============================================================================
# [8] F2 = 2N + S - Delta
# ==============================================================================

print("[8] F2 SYMMETRIC DECOMPOSITION")
print("-" * 78)

certificate(
    "F2 = 2N + S - Delta",
    F2,
    2*N_xy + S_xy - Delta_xy,
)

print("  Therefore:")
print("      F2 = 2N + S - Delta")
print()
print("  In historical coordinates this is:")
print("      N     -> conic/product channel")
print("      S     -> 4y/r")
print("      Delta -> 4(2x-3)^2/r^2")
print()


# ==============================================================================
# [9] HISTORICAL FACTOR SUM / GAP / PRODUCT
# ==============================================================================

print("[9] A0 / B0 TRIAD")
print("-" * 78)

factor_sum = simplify_exact(A0 + B0)
factor_gap = simplify_exact(B0 - A0)
factor_product = simplify_exact(A0 * B0)

show("A0+B0", factor_sum)
show("B0-A0", factor_gap)
show("A0*B0", factor_product)

certificate(
    "A0+B0 = 4y",
    factor_sum,
    4*y,
)

certificate(
    "B0-A0 = 4x-6",
    factor_gap,
    4*x - 6,
)

certificate(
    "A0*B0 = 4v-1",
    factor_product,
    4*v_xy - 1,
)


# ==============================================================================
# [10] MODERN COORDINATES FROM A0/B0
# ==============================================================================

print("[10] MODERN COORDINATES FROM HISTORICAL FACTOR COORDINATES")
print("-" * 78)

S_from_A0B0 = simplify_exact(
    (A0 + B0) / r
)

N_from_A0B0 = simplify_exact(
    A0 * B0 / r**2
)

Delta_from_A0B0 = simplify_exact(
    (A0 - B0)**2 / r**2
)

show("S from A0+B0", S_from_A0B0)
show("N from A0*B0", N_from_A0B0)
show("Delta from (A0-B0)^2", Delta_from_A0B0)

certificate(
    "S from A0+B0",
    S_from_A0B0,
    S_xy,
)

certificate(
    "N from A0*B0",
    N_from_A0B0,
    N_xy,
)

certificate(
    "Delta from A0-B0",
    Delta_from_A0B0,
    Delta_xy,
)


# ==============================================================================
# [11] KAPPA QUADRATIC DIRECTLY IN x,y
# ==============================================================================

print("[11] KAPPA QUADRATIC")
print("-" * 78)

Q_xy = simplify_exact(
    z**2 - S_xy*z + N_xy
)

Q_factor = simplify_exact(
    (z - A0/r) * (z - B0/r)
)

show("Q(z)", Q_xy)
show("factorized Q(z)", Q_factor)

certificate(
    "Q=(z-A0/r)(z-B0/r)",
    Q_xy,
    Q_factor,
)


# ==============================================================================
# [12] KAPPA +1 TRANSLATION
# ==============================================================================

print("[12] Q(z-1)")
print("-" * 78)

Q_shift = simplify_exact(
    Q_xy.subs(z, z - 1)
)

Q_shift_expected = simplify_exact(
    (z - A0/r - 1)
    * (z - B0/r - 1)
)

show("Q(z-1)", Q_shift)

certificate(
    "Q(z-1)",
    Q_shift,
    Q_shift_expected,
)


# ==============================================================================
# [13] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ==============================================================================

print("[13] FOUR-BASE CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi = simplify_exact(
    Q_xy * Q_shift
)

show("chi(z)", chi)


# ==============================================================================
# [14] DISCRIMINANT DIRECTLY IN x,y
# ==============================================================================

print("[14] DISCRIMINANT")
print("-" * 78)

disc_Q = simplify_exact(
    sp.discriminant(
        sp.Poly(Q_xy, z),
        z,
    )
)

show("disc(Q)", disc_Q)

certificate(
    "disc(Q)=Delta",
    disc_Q,
    Delta_xy,
)


# ==============================================================================
# [15] S^2 - 4N DIRECTLY IN HISTORICAL COORDINATES
# ==============================================================================

print("[15] DISCRIMINANT IDENTITY")
print("-" * 78)

certificate(
    "S^2-4N=Delta",
    S_xy**2 - 4*N_xy,
    Delta_xy,
)


# ==============================================================================
# [16] M1 HISTORICAL FORM
# ==============================================================================

print("[16] M1")
print("-" * 78)

M1_factor = simplify_exact(
    (A0 + r) * (B0 + r) / r**2
)

show("M1 factor form", M1_factor)

certificate(
    "M1=(A0+r)(B0+r)/r^2",
    M1_xy,
    M1_factor,
)


# ==============================================================================
# [17] HISTORICAL TRANSLATION y -> y + h*r/2
# ==============================================================================

print("[17] HISTORICAL TRANSLATION")
print("-" * 78)

h = sp.symbols("h")

x_prime = x
y_prime = y + h*r/2

N_prime = simplify_exact(
    (
        (2*y_prime)**2
        - (2*x_prime - 3)**2
    )
    / r**2
)

S_prime = simplify_exact(
    4*y_prime / r
)

F2_prime = simplify_exact(
    6*N_prime
    - S_prime**2
    + S_prime
)

show("N'", N_prime)
show("S'", S_prime)
show("F2'", F2_prime)

certificate(
    "N'=N+hS+h^2",
    N_prime,
    N_xy + h*S_xy + h**2,
)

certificate(
    "S'=S+2h",
    S_prime,
    S_xy + 2*h,
)


# ==============================================================================
# [18] F2 TRANSLATION RESPONSE
# ==============================================================================

print("[18] F2 TRANSLATION RESPONSE")
print("-" * 78)

F2_difference = simplify_exact(
    F2_prime - F2
)

show("F2'-F2", F2_difference)

certificate(
    "F2'-F2 = 2h(S+h+1)",
    F2_difference,
    2*h*(S_xy + h + 1),
)

S_recovered_from_translation = simplify_exact(
    F2_difference / (2*h) - h - 1
)

show(
    "S recovered from F2 translation",
    S_recovered_from_translation,
)

certificate(
    "S translation recovery",
    S_recovered_from_translation,
    S_xy,
)


# ==============================================================================
# [19] TRANSLATION IN A0/B0
# ==============================================================================

print("[19] HISTORICAL FACTOR TRANSLATION")
print("-" * 78)

A0_shift = simplify_exact(
    A0 + h*r
)

B0_shift = simplify_exact(
    B0 + h*r
)

show("A0'", A0_shift)
show("B0'", B0_shift)

certificate(
    "A0' = A0 + hr",
    A0_shift,
    A0 + h*r,
)

certificate(
    "B0' = B0 + hr",
    B0_shift,
    B0 + h*r,
)


# ==============================================================================
# [20] TRANSLATED PRODUCT
# ==============================================================================

print("[20] TRANSLATED HISTORICAL PRODUCT")
print("-" * 78)

translated_product = simplify_exact(
    A0_shift * B0_shift / r**2
)

show("translated L(h)", translated_product)

certificate(
    "L(h)=N+hS+h^2",
    translated_product,
    N_xy + h*S_xy + h**2,
)


# ==============================================================================
# [21] TRANSLATION QUADRATIC
# ==============================================================================

print("[21] TRANSLATION QUADRATIC")
print("-" * 78)

translation_quadratic = simplify_exact(
    h**2 + S_xy*h + N_xy
)

show("L(h)", translation_quadratic)

certificate(
    "translated product = L(h)",
    translated_product,
    translation_quadratic,
)

Q_from_orbit = simplify_exact(
    translation_quadratic.subs(h, -z)
)

certificate(
    "L(-z)=Q(z)",
    Q_from_orbit,
    Q_xy,
)


# ==============================================================================
# [22] DIRECT HISTORICAL COORDINATE PACKAGE
# ==============================================================================

print("[22] COMPLETE x,y -> KAPPA PACKAGE")
print("-" * 78)

print(
    """
  Historical coordinates:

      A0 = 2y - 2x + 3
      B0 = 2y + 2x - 3

  Therefore:

      A0+B0 = 4y
      B0-A0 = 4x-6
      A0*B0 = 4v-1

  Modern KAPPA coordinates:

      N     = A0*B0/r^2
      S     = (A0+B0)/r
      Delta = (A0-B0)^2/r^2

  Hence:

      S = 4y/r

      Delta = 4(2x-3)^2/r^2

      N = [4y^2-(2x-3)^2]/r^2

  and:

      Q(z) = z^2-Sz+N

            = (z-A0/r)(z-B0/r).

  The translation y -> y+hr/2 gives:

      A0 -> A0+hr
      B0 -> B0+hr

  and therefore:

      L(h) = N+hS+h^2.

  So the historical coordinates already contain the complete
  modern symmetric KAPPA coordinate system.
"""
)


# ==============================================================================
# [23] IMPORTANT STRUCTURAL COMPARISON
# ==============================================================================

print("[23] STRUCTURAL COMPARISON")
print("-" * 78)

print(
    """
  The relevant channels are now explicitly separated:

      symmetric:
          A0+B0 = 4y

      antisymmetric:
          B0-A0 = 4x-6

      product:
          A0*B0 = 4v-1

  Consequently:

      S     comes from y,
      Delta comes from x,
      N     comes from the conic/product.

  This is the direct 2020 -> 2026 coordinate crosswalk.

  No N-only reconstruction of S is involved.
"""
)


# ==============================================================================
# [24] FINAL AUDIT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 551 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print(
    """
The next experiment should use the historical coordinates as the
primary variables rather than treating S as an unknown.

The strongest established identities are:

    S = 4y/r

    Delta = 4(2x-3)^2/r^2

    N = [4y^2-(2x-3)^2]/r^2

    Q(z) = (z-A0/r)(z-B0/r)

The next useful search is therefore:

    rewrite the remaining homogeneous-layer quantities directly
    in terms of:

        y
        2x-3
        A0
        B0
        4v-1

and look for factorization or cancellation that is invisible in
the N,S,Delta representation.

In particular, test whether additional KAPPA layers contain simple
factors of:

    A0
    B0
    A0+B0
    A0-B0
    A0*B0

rather than searching for S indirectly.
"""
)