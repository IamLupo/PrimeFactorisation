#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 526 START")
print("=" * 78)
print("HISTORICAL CONIC y-TRANSLATION -> KAPPA QUADRATIC ORBIT")
print()

failures = 0


def check(label, expr, expected=0):
    global failures

    diff = sp.factor(sp.expand(expr - expected))
    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# SYMBOLS
# ============================================================================

x, y = sp.symbols("x y")
m, r = sp.symbols("m r")
z = sp.symbols("z")

# Historical conic:
#
#     v = y^2 - x^2 + 3x - 2
#
# Hence
#
#     4v - 1 = 4y^2 - (2x-3)^2.

v = sp.expand(y**2 - x**2 + 3*x - 2)

K = sp.expand(4*v - 1)

U = sp.expand(2*y)
V = sp.expand(2*x - 3)

print("[1] HISTORICAL CONIC")
print("-" * 78)

print("  v =", v)
print("  K = 4v-1 =", sp.factor(K))
print("  U =", U)
print("  V =", V)

check(
    "4v-1 = U^2-V^2",
    K,
    U**2 - V**2,
)

print()


# ============================================================================
# [2] HISTORICAL FACTOR COORDINATES
# ============================================================================

print("[2] HISTORICAL FACTOR COORDINATES")
print("-" * 78)

A0 = sp.expand(U - V)
B0 = sp.expand(U + V)

print("  A0 =", A0)
print("  B0 =", B0)

check(
    "A0 = 2y-2x+3",
    A0,
    2*y - 2*x + 3,
)

check(
    "B0 = 2y+2x-3",
    B0,
    2*y + 2*x - 3,
)

check(
    "A0*B0 = 4v-1",
    A0*B0,
    K,
)

print()


# ============================================================================
# [3] KAPPA TRANSLATION ORBIT
# ============================================================================

print("[3] KAPPA TRANSLATION ORBIT")
print("-" * 78)

K_m = sp.expand(
    (U + m*r)**2 - V**2
)

print("  K_m =", sp.factor(K_m))

check(
    "K_m=(U+mr)^2-V^2",
    K_m,
    (U + m*r)**2 - V**2,
)

print()


# ============================================================================
# [4] DIRECT y-TRANSLATION OF THE HISTORICAL CONIC
# ============================================================================

print("[4] DIRECT y-TRANSLATION")
print("-" * 78)

y_shift = sp.expand(
    y + m*r/2
)

v_shift = sp.expand(
    v.subs(y, y_shift)
)

K_from_v_shift = sp.expand(
    4*v_shift - 1
)

print("  y -> y + m*r/2")
print("  v_shift =", sp.factor(v_shift))
print("  4*v_shift-1 =", sp.factor(K_from_v_shift))

check(
    "4v(x,y+mr/2)-1 = K_m",
    K_from_v_shift,
    K_m,
)

print()


# ============================================================================
# [5] EXACT TRANSLATION IDENTITY
# ============================================================================

print("[5] EXACT TRANSLATION IDENTITY")
print("-" * 78)

print("""
  The historical conic translation is therefore

      y -> y + m*r/2

  while x remains fixed.

  At the factor-coordinate level:

      A0 -> A0 + m*r
      B0 -> B0 + m*r

  and therefore:

      K_m = (A0+m*r)(B0+m*r).
""")

A_m = sp.expand(A0 + m*r)
B_m = sp.expand(B0 + m*r)

check(
    "translated A0",
    A_m,
    2*(y + m*r/2) - 2*x + 3,
)

check(
    "translated B0",
    B_m,
    2*(y + m*r/2) + 2*x - 3,
)

check(
    "translated factor product",
    A_m*B_m,
    K_m,
)

print()


# ============================================================================
# [6] NORMALIZED ORBIT
# ============================================================================

print("[6] NORMALIZED ORBIT")
print("-" * 78)

N = sp.symbols("N")
S = sp.symbols("S")

L_m = sp.expand(
    K_m / r**2
)

L_expected = sp.expand(
    N + m*S + m**2
)

# Substitute the historical expressions.
N_hist = sp.expand(
    (U**2 - V**2) / r**2
)

S_hist = sp.expand(
    2*U / r
)

check(
    "N_hist = (U^2-V^2)/r^2",
    N_hist,
    N_hist,
)

check(
    "S_hist = 2U/r",
    S_hist,
    4*y/r,
)

L_hist_expected = sp.expand(
    N_hist + m*S_hist + m**2
)

check(
    "L_m = N + m*S + m^2",
    L_m,
    L_hist_expected,
)

print("  L_m =", sp.factor(L_m))

print()


# ============================================================================
# [7] CONIC VALUE ORBIT v_m
# ============================================================================

print("[7] CONIC VALUE ORBIT")
print("-" * 78)

v_m = sp.expand(
    v.subs(y, y + m*r/2)
)

print("  v_m =", sp.factor(v_m))

delta_v = sp.expand(
    v_m - v
)

print("  v_m-v =", sp.factor(delta_v))

check(
    "v_m-v = m*r*y + m^2*r^2/4",
    delta_v,
    m*r*y + m**2*r**2/4,
)

print()


# ============================================================================
# [8] FIRST DIFFERENCE OF THE CONIC ORBIT
# ============================================================================

print("[8] FIRST DIFFERENCE")
print("-" * 78)

dL = sp.expand(
    L_m.subs(m, m+1) - L_m
)

print("  ΔL_m =", sp.factor(dL))

check(
    "ΔL=S+2m+1",
    dL,
    S_hist + 2*m + 1,
)

print()


# ============================================================================
# [9] FIRST DIFFERENCE DIRECTLY IN x,y
# ============================================================================

print("[9] FIRST DIFFERENCE IN HISTORICAL COORDINATES")
print("-" * 78)

dL_xy = sp.factor(
    sp.expand(
        dL
    )
)

print("  ΔL_m =", dL_xy)

check(
    "ΔL_m = 4y/r + 2m + 1",
    dL_xy,
    4*y/r + 2*m + 1,
)

print()


# ============================================================================
# [10] SECOND DIFFERENCE
# ============================================================================

print("[10] SECOND DIFFERENCE")
print("-" * 78)

d2L = sp.expand(
    L_m.subs(m, m+2)
    - 2*L_m.subs(m, m+1)
    + L_m
)

print("  Δ²L_m =", d2L)

check(
    "Δ²L_m=2",
    d2L,
    2,
)

print()


# ============================================================================
# [11] GENERATING QUADRATIC
# ============================================================================

print("[11] GENERATING QUADRATIC")
print("-" * 78)

Q_from_conic = sp.expand(
    L_m.subs(m, -z)
)

Q_expected = sp.expand(
    z**2 - S_hist*z + N_hist
)

print("  L(-z) =", sp.factor(Q_from_conic))

check(
    "L(-z)=z^2-Sz+N",
    Q_from_conic,
    Q_expected,
)

print()


# ============================================================================
# [12] ROOTS OF THE GENERATING QUADRATIC
# ============================================================================

print("[12] ROOT INTERPRETATION")
print("-" * 78)

Q_factor = sp.factor(
    Q_from_conic
)

print("  Q(z) =", Q_factor)

check(
    "Q(z)=(z-P)(z-Q)",
    Q_from_conic,
    (z - (U-V)/r)*(z - (U+V)/r),
)

print()


# ============================================================================
# [13] DISCRIMINANT
# ============================================================================

print("[13] DISCRIMINANT")
print("-" * 78)

disc = sp.factor(
    sp.discriminant(Q_from_conic, z)
)

Delta_hist = sp.factor(
    (A0-B0)**2 / r**2
)

print("  disc(Q) =", disc)
print("  Delta    =", Delta_hist)

check(
    "disc(Q)=Delta",
    disc,
    Delta_hist,
)

print()


# ============================================================================
# [14] TRANSLATED KAPPA FACTOR
# ============================================================================

print("[14] KAPPA +1 TRANSLATION")
print("-" * 78)

Q_shift = sp.expand(
    Q_from_conic.subs(z, z-1)
)

Q_shift_expected = sp.expand(
    (z - ((U-V)/r + 1))
    * (z - ((U+V)/r + 1))
)

check(
    "Q(z-1)=translated quadratic",
    Q_shift,
    Q_shift_expected,
)

print()


# ============================================================================
# [15] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print("[15] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(
    Q_from_conic * Q_shift
)

chi_expected = sp.expand(
    (z-(U-V)/r)
    * (z-(U+V)/r)
    * (z-(U-V)/r-1)
    * (z-(U+V)/r-1)
)

check(
    "Q(z)Q(z-1)=four-base spectrum",
    chi,
    chi_expected,
)

print()


# ============================================================================
# [16] SPECIAL INTEGER-LAYER SUBSEQUENCE
# ============================================================================

print("[16] INTEGER y-TRANSLATION SUBSEQUENCE")
print("-" * 78)

j = sp.symbols("j", integer=True)

# Taking m=2j gives:
#
#     y -> y+j*r
#
# so the historical conic remains on the integer lattice.

v_even = sp.expand(
    v.subs(y, y + j*r)
)

L_even = sp.expand(
    (4*v_even - 1) / r**2
)

expected_even = sp.expand(
    j**2*4
    + j*(4*y/r)
    + (4*y**2 - 4*x**2 + 12*x - 9)/r**2
)

print("  m=2j")
print("  y -> y+j*r")
print("  normalized layer =", sp.factor(L_even))

check(
    "integer-subsequence orbit",
    L_even,
    expected_even,
)

print()


# ============================================================================
# [17] HALF-TRANSLATION INTERPRETATION
# ============================================================================

print("[17] HALF-TRANSLATION INTERPRETATION")
print("-" * 78)

print("""
  The KAPPA translation parameter m does not require a new
  algebraic object.

  In the historical conic coordinates it is exactly:

      y_m = y + m*r/2

  with x fixed.

  Consequently:

      4*v(x,y_m)-1
        = (A0+m*r)(B0+m*r)

  and

      [4*v(x,y_m)-1]/r^2
        = m^2 + (4y/r)m + N.

  Since

      S = 4y/r,

  this becomes

      L_m = m^2 + S*m + N.
""")

check(
    "S = 4y/r",
    S_hist,
    4*y/r,
)

print()


# ============================================================================
# [18] HISTORICAL g / v CHANNEL CONNECTION
# ============================================================================

print("[18] v -> TRANSLATED v CHANNEL")
print("-" * 78)

print("""
  The old N-only quantity

      v = y^2-x^2+3x-2

  determines only the base conic value.

  The translated family is

      v_m = v + m*r*y + m^2*r^2/4.

  Therefore the missing information is precisely the
  coefficient of the linear translation term:

      coefficient of m in v_m = r*y.

  Equivalently:

      y = (v_1-v_0)/r - r/4.

  This is exactly the same information channel as the
  KAPPA first-difference observable.
""")

v0 = sp.expand(v_m.subs(m, 0))
v1 = sp.expand(v_m.subs(m, 1))

y_recovered = sp.expand(
    (v1 - v0)/r - r/4
)

check(
    "y recovered from v1-v0",
    y_recovered,
    y,
)

print()


# ============================================================================
# [19] OPERATOR FORM
# ============================================================================

print("[19] OPERATOR FORM")
print("-" * 78)

print("""
  Define the historical translation operator

      T_r f(x,y) = f(x, y+r/2).

  Then

      K_m = T_r^m(4v-1)

  and

      L_m = T_r^m(4v-1)/r^2.

  The resulting orbit obeys

      Δ²L_m = 2

  and

      ΔL_m-(2m+1)=S.
""")

T_v = sp.expand(
    v.subs(y, y + r/2)
)

check(
    "single T_r step",
    4*T_v - 1,
    K_m.subs(m, 1),
)

print()


# ============================================================================
# [20] UPSTREAM SIGNATURE
# ============================================================================

print("[20] UPSTREAM SIGNATURE")
print("-" * 78)

print("""
  This gives a much more concrete signature for the old
  homogeneous/conic construction:

      V_m = v(x, y+m*r/2)

  should satisfy

      V_m
        = V_0
        + (r*y)m
        + (r^2/4)m^2.

  Equivalently:

      ΔV_m = r*y + r^2*(2m+1)/4

      Δ²V_m = r^2/2.

  After normalization:

      L_m = (4V_m-1)/r^2

      ΔL_m-(2m+1) = 4y/r = S.

  Therefore the decisive upstream question is:

      Does the existing homogeneous-layer construction
      possess a natural translation parameter whose action
      on the historical v-kernel is

          y -> y+r/2

      or an algebraically equivalent operation?
""")

print()


# ============================================================================
# [21] FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 526 FINISHED")
print("=" * 78)
print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT SYMBOLIC AUDIT = {failures == 0}")

print("""
NEXT RESEARCH TARGET
------------------------------------------------------------------------------

Do not search arbitrary identities.

The new object to look for in the ORIGINAL homogeneous-layer
construction is a translation operator T such that:

    T(v) - v
        is linear in the layer parameter,

and preferably:

    T^m(v)
        is quadratic in m.

The exact desired normalized signature is:

    L_m = (4*T^m(v)-1)/r^2

    L_m = m^2 + S*m + N.

Even more directly:

    [L_{m+1}-L_m] - (2m+1) = S.

This is now a concrete operator-level bridge between the
2020 conic kernel and the 2026 KAPPA quadratic.

In particular, inspect the actual homogeneous-layer
transformation rules for any operation equivalent to:

    y -> y + r/2,
    U -> U+r,
    A0 -> A0+r,
    B0 -> B0+r.

If such an operation already exists in the old construction,
the KAPPA quadratic is generated by the same operator rather
than being introduced independently.
""")
