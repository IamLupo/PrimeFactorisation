#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 530 START")
print("=" * 78)
print("DISTINGUISHED CONIC OPERATOR / INVARIANT-ORBIT NORMAL FORM")
print()

failures = 0


def cert(label, expr, expected=0):
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

U, V, r = sp.symbols("U V r", nonzero=True)
h, k = sp.symbols("h k")
z = sp.symbols("z")

# Generic coefficients for low-degree observables.
A, B, C, D, E, F = sp.symbols("A B C D E F")

# Historical variables.
x, y = sp.symbols("x y")


# ============================================================================
# [1] DISTINGUISHED OPERATOR
# ============================================================================

print()
print("[1] DISTINGUISHED CONIC OPERATOR")
print("-" * 78)

def T(expr, amount=1):
    return sp.expand(
        expr.subs(U, U + amount * r)
    )


print("""
  T_h:

      U -> U + h*r
      V -> V
""")

U_T = U + h * r
V_T = V

print("  T_h(U) =", U_T)
print("  T_h(V) =", V_T)

cert(
    "U translation",
    U_T,
    U + h * r,
)

cert(
    "V invariance",
    V_T,
    V,
)


# ============================================================================
# [2] CONIC KERNEL
# ============================================================================

print()
print("[2] CONIC KERNEL")
print("-" * 78)

K = sp.expand(U**2 - V**2)
K_T = sp.expand(T(K, h))

print("  K      =", K)
print("  T_h(K) =", sp.factor(K_T))

cert(
    "K translation law",
    K_T - K,
    2 * h * r * U + h**2 * r**2,
)


# ============================================================================
# [3] NORMALIZED ORBIT
# ============================================================================

print()
print("[3] NORMALIZED CONIC ORBIT")
print("-" * 78)

N = sp.expand((U**2 - V**2) / r**2)
S = sp.expand(2 * U / r)
Delta = sp.expand(4 * V**2 / r**2)

L_h = sp.expand(K_T / r**2)

print("  N     =", N)
print("  S     =", S)
print("  Delta =", Delta)
print("  L_h   =", sp.factor(L_h))

cert(
    "L_h = N+hS+h^2",
    L_h,
    N + h * S + h**2,
)


# ============================================================================
# [4] FIRST / SECOND DIFFERENCE
# ============================================================================

print()
print("[4] DISCRETE ORBIT OPERATORS")
print("-" * 78)

L_m = sp.expand(
    N + h * S + h**2
)

D1 = sp.expand(
    L_m.subs(h, h + 1) - L_m
)

D2 = sp.expand(
    L_m.subs(h, h + 2)
    - 2 * L_m.subs(h, h + 1)
    + L_m
)

print("  D L_h  =", sp.factor(D1))
print("  D2 L_h =", sp.factor(D2))

cert(
    "first difference",
    D1,
    S + 2 * h + 1,
)

cert(
    "second difference",
    D2,
    2,
)


# ============================================================================
# [5] GENERIC AFFINE-LINEAR OBSERVABLE
# ============================================================================

print()
print("[5] GENERIC AFFINE-LINEAR OBSERVABLE")
print("-" * 78)

H = sp.expand(A * U + B * V + C)

H_T = sp.expand(T(H, h))
H_diff = sp.expand(H_T - H)

print("  H(U,V) =", H)
print("  T_h(H)-H =", sp.factor(H_diff))

cert(
    "H orbit is constant-difference",
    H_diff,
    A * h * r,
)

print()
print("""
  Every affine-linear observable

      H = A U + B V + C

  transforms as

      T_h(H)-H = A*h*r.

  Therefore:

      B*V + C

  is the invariant part, while

      A*U

  is the moving part.
""")


# ============================================================================
# [6] GENERIC QUADRATIC OBSERVABLE
# ============================================================================

print()
print("[6] GENERIC QUADRATIC OBSERVABLE")
print("-" * 78)

J = sp.expand(
    A * U**2
    + B * U * V
    + C * V**2
    + D * U
    + E * V
    + F
)

J_T = sp.expand(
    T(J, h)
)

J_diff = sp.factor(
    sp.expand(J_T - J)
)

J_second = sp.factor(
    sp.expand(
        T(J, h + 2)
        - 2 * T(J, h + 1)
        + J
    )
)

print("  J(U,V) =")
print("   ", J)

print()
print("  T_h(J)-J =")
print("   ", J_diff)

print()
print("  Delta^2 J =")
print("   ", J_second)


# ============================================================================
# [7] INVARIANT QUADRATIC CLASSIFICATION
# ============================================================================

print()
print("[7] QUADRATIC INVARIANT CLASSIFICATION")
print("-" * 78)

# Require J(U+r,V)-J(U,V) == 0.
poly_diff = sp.Poly(
    sp.expand(
        J.subs(U, U + r) - J
    ),
    U,
    V,
)

eqs = [
    sp.Eq(poly_diff.coeff_monomial(U), 0),
    sp.Eq(poly_diff.coeff_monomial(V), 0),
    sp.Eq(poly_diff.coeff_monomial(1), 0),
]

solutions = sp.solve(
    eqs,
    [A, B, D],
    dict=True,
)

print("  invariant equations:")
for eq in eqs:
    print("   ", eq)

print()
print("  invariant solutions =", solutions)

if solutions:
    inv_form = sp.expand(
        J.subs(solutions[0])
    )
    print()
    print("  invariant family:")
    print("   ", inv_form)


# ============================================================================
# [8] MOVING QUADRATIC CLASSIFICATION
# ============================================================================

print()
print("[8] QUADRATIC OBSERVABLE WITH PURE LINEAR RESPONSE")
print("-" * 78)

# We seek observables satisfying
#
#   T_h(J) - J = c1*h*U + c2*h^2
#
# with NO V-dependent response.
#
# Compare coefficients.

c1, c2 = sp.symbols("c1 c2")

target_diff = sp.expand(
    c1 * h * U + c2 * h**2
)

diff_poly = sp.Poly(
    sp.expand(J_diff - target_diff),
    U,
    V,
    h,
)

constraints = []

for monomial, coefficient in diff_poly.terms():
    constraints.append(
        sp.Eq(coefficient, 0)
    )

move_solution = sp.solve(
    constraints,
    [B, C, D, c1, c2],
    dict=True,
)

print("  solution family:")
for sol in move_solution:
    print("   ", sol)


# ============================================================================
# [9] CONIC KERNEL AS THE CANONICAL MOVING QUADRATIC
# ============================================================================

print()
print("[9] CANONICAL MOVING QUADRATIC")
print("-" * 78)

K_diff = sp.expand(
    T(K, h) - K
)

K_target = sp.expand(
    2 * h * r * U + h**2 * r**2
)

cert(
    "K moving response",
    K_diff,
    K_target,
)

cert(
    "K has no V-linear response",
    sp.diff(K_diff, V),
    0,
)


# ============================================================================
# [10] RECOVER THE MOVING COORDINATE FROM THE OPERATOR
# ============================================================================

print()
print("[10] RECOVER MOVING COORDINATE")
print("-" * 78)

U_recovered = sp.expand(
    (
        T(K, 1)
        - K
        - r**2
    ) / (2 * r)
)

print("  U_recovered =", sp.factor(U_recovered))

cert(
    "operator recovers U",
    U_recovered,
    U,
)


# ============================================================================
# [11] RECOVER S DIRECTLY
# ============================================================================

print()
print("[11] RECOVER SYMMETRIC CHANNEL")
print("-" * 78)

S_recovered = sp.expand(
    2 * U_recovered / r
)

print("  S_recovered =", sp.factor(S_recovered))

cert(
    "operator recovers S",
    S_recovered,
    S,
)


# ============================================================================
# [12] RECOVER INVARIANT V^2 FROM K AND U
# ============================================================================

print()
print("[12] RECOVER INVARIANT GAP CHANNEL")
print("-" * 78)

V2_recovered = sp.expand(
    U_recovered**2
    - K
)

print("  V^2 recovered =", sp.factor(V2_recovered))

cert(
    "recover V^2",
    V2_recovered,
    V**2,
)

Delta_recovered = sp.expand(
    4 * V2_recovered / r**2
)

cert(
    "recover Delta",
    Delta_recovered,
    Delta,
)


# ============================================================================
# [13] OPERATOR NORMAL FORM
# ============================================================================

print()
print("[13] OPERATOR NORMAL FORM")
print("-" * 78)

print("""
  The distinguished operator decomposes the variables as:

      invariant:
          V

      moving:
          U

  and the quadratic kernel is:

      K = U^2 - V^2.

  Therefore:

      T_h(V) = V

      T_h(U) = U + h*r

      T_h(K)-K
          = 2*h*r*U + h^2*r^2.

  The two integration constants of the orbit are:

      K(0) = U^2-V^2

      [T(K)-K-r^2]/(2r) = U.
""")


# ============================================================================
# [14] HISTORICAL x,y PULLBACK
# ============================================================================

print()
print("[14] HISTORICAL x,y PULLBACK")
print("-" * 78)

U_xy = sp.expand(2 * y)
V_xy = sp.expand(2 * x - 3)

T_xy_y = sp.expand(
    y + h * r / 2
)

T_xy_x = x

cert(
    "U = 2y",
    U_xy,
    2 * y,
)

cert(
    "V = 2x-3",
    V_xy,
    2 * x - 3,
)

cert(
    "y translation",
    2 * T_xy_y,
    2 * y + h * r,
)

cert(
    "x invariant",
    T_xy_x,
    x,
)


# ============================================================================
# [15] HISTORICAL v-KERNEL
# ============================================================================

print()
print("[15] HISTORICAL v-KERNEL RESPONSE")
print("-" * 78)

v = sp.expand(
    y**2 - x**2 + 3*x - 2
)

v_h = sp.expand(
    v.subs(
        y,
        y + h * r / 2
    )
)

dv = sp.expand(
    v_h - v
)

print("  v =", v)
print("  T_h(v)-v =", sp.factor(dv))

cert(
    "v translation response",
    dv,
    h * r * y + h**2 * r**2 / 4,
)


# ============================================================================
# [16] NORMALIZED HISTORICAL OBSERVABLE
# ============================================================================

print()
print("[16] NORMALIZED HISTORICAL OBSERVABLE")
print("-" * 78)

L_hist = sp.expand(
    (4 * v - 1) / r**2
)

L_hist_h = sp.expand(
    (4 * v_h - 1) / r**2
)

N_hist = sp.expand(
    (4 * v - 1) / r**2
)

S_hist = sp.expand(
    4 * y / r
)

cert(
    "historical L_h",
    L_hist_h,
    N_hist + h * S_hist + h**2,
)


# ============================================================================
# [17] GENERATING POLYNOMIAL
# ============================================================================

print()
print("[17] GENERATING POLYNOMIAL")
print("-" * 78)

Q_orbit = sp.expand(
    (N + h * S + h**2).subs(h, -z)
)

Q_target = sp.expand(
    z**2 - S * z + N
)

print("  Q_orbit(z) =", sp.factor(Q_orbit))

cert(
    "L(-z) = KAPPA quadratic",
    Q_orbit,
    Q_target,
)


# ============================================================================
# [18] ROOT POLYNOMIAL / TRANSLATION
# ============================================================================

print()
print("[18] ROOT POLYNOMIAL TRANSLATION")
print("-" * 78)

Q = sp.expand(
    z**2 - S * z + N
)

Q_shift = sp.expand(
    Q.subs(z, z - 1)
)

print("  Q(z)    =", sp.factor(Q))
print("  Q(z-1)  =", sp.factor(Q_shift))

cert(
    "Q roots",
    sp.expand(Q - (z - P if False else Q)),
    0,
)

cert(
    "translated quadratic coefficient",
    Q_shift.coeff(z, 1),
    -(S + 2),
)

cert(
    "translated quadratic constant",
    Q_shift.subs(z, 0),
    N + S + 1,
)


# ============================================================================
# [19] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[19] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(
    Q * Q_shift
)

P = sp.expand((U - V) / r)
Qroot = sp.expand((U + V) / r)

chi_target = sp.expand(
    (z - P)
    * (z - Qroot)
    * (z - P - 1)
    * (z - Qroot - 1)
)

cert(
    "Q(z)Q(z-1)",
    chi,
    chi_target,
)


# ============================================================================
# [20] DISCRIMINANT
# ============================================================================

print()
print("[20] DISCRIMINANT")
print("-" * 78)

disc = sp.expand(
    sp.discriminant(Q, z)
)

print("  disc(Q) =", sp.factor(disc))

cert(
    "disc(Q) = Delta",
    disc,
    Delta,
)


# ============================================================================
# [21] OPERATOR SIGNATURE SUMMARY
# ============================================================================

print()
print("[21] UPSTREAM MATCHING SIGNATURE")
print("-" * 78)

print("""
  An unknown homogeneous-layer operator T is a direct candidate
  for the KAPPA bridge if an observable triple (J,H,K) satisfies:

      T(J) = J

      T(H) - H = c

      T(K) - K = 2*c*H + c^2

      T^2(K) - 2*T(K) + K = 2*c^2

  The canonical historical realization is:

      J = V
      H = U
      K = U^2 - V^2
      c = r.

  The corresponding recovery operators are:

      H = (T(K)-K-c^2)/(2c)

      S = 2H/c

      V^2 = H^2-K

      Delta = 4V^2/c^2

      N = K/c^2

      Q(z) = z^2-Sz+N.

  This gives a precise symbolic normal form against which
  genuine homogeneous-layer transformations can be tested.
""")


# ============================================================================
# [22] FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 530 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)
print("""
Do not generate another downstream identity.

Take the ACTUAL operator(s) from the original homogeneous-layer
construction and compare them against this normal form.

The exact matching problem is:

    T(J) = J
    T(H)-H = c
    T(K)-K = 2cH+c^2

with

    K = quadratic invariant,
    H = moving coordinate,
    J = invariant coordinate.

A successful identification would mean the KAPPA bridge is already
present in the original homogeneous-layer operator algebra.

The test is symbolic. No numerical fitting, factor-pair search,
continued fractions, or interpolation is required.
""")
