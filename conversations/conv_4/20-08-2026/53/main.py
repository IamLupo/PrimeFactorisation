#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 529 START")
print("=" * 78)
print("AFFINE CONIC OPERATOR NORMAL FORM / KAPPA TRANSLATION CONJUGACY")
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

U, V = sp.symbols("U V")
r = sp.symbols("r", nonzero=True)

a, b, c, d = sp.symbols("a b c d")
alpha, beta = sp.symbols("alpha beta")
h, k = sp.symbols("h k")

# IMPORTANT: spectral polynomial variable
z = sp.symbols("z")

x, y = sp.symbols("x y")

K = sp.expand(U**2 - V**2)

S = sp.expand(2 * U / r)
Delta = sp.expand(4 * V**2 / r**2)
N = sp.expand(K / r**2)


# ============================================================================
# [1] GENERAL AFFINE CONIC OPERATOR
# ============================================================================

print()
print("[1] GENERAL AFFINE CONIC OPERATOR")
print("-" * 78)

U1 = a * U + b * V + alpha
V1 = c * U + d * V + beta

K1 = sp.expand(U1**2 - V1**2)

print("  U' = a U + b V + alpha")
print("  V' = c U + d V + beta")
print()
print("  K' =", sp.factor(K1))
print()


# ============================================================================
# [2] HOMOGENEOUS CONIC-PRESERVING CONDITIONS
# ============================================================================

print("[2] HOMOGENEOUS CONIC-PRESERVING CONDITIONS")
print("-" * 78)

linear_map_K = sp.expand(
    (a * U + b * V)**2 -
    (c * U + d * V)**2
)

poly_K = sp.Poly(
    sp.expand(linear_map_K - K),
    U,
    V,
)

coeff_U2 = sp.expand(poly_K.coeff_monomial(U**2))
coeff_UV = sp.expand(poly_K.coeff_monomial(U * V))
coeff_V2 = sp.expand(poly_K.coeff_monomial(V**2))

print("  coefficient(U^2) =", coeff_U2)
print("  coefficient(UV)  =", coeff_UV)
print("  coefficient(V^2) =", coeff_V2)

cert("U^2 preservation condition", coeff_U2, 0)
cert("UV preservation condition", coeff_UV, 0)
cert("V^2 preservation condition", coeff_V2, 0)

print()
print("""
  The homogeneous part therefore satisfies:

      a^2-c^2 = 1
      ab-cd   = 0
      b^2-d^2 = -1
""")
print()


# ============================================================================
# [3] GENERAL PURE TRANSLATION
# ============================================================================

print("[3] GENERAL PURE TRANSLATION")
print("-" * 78)

Utr = U + alpha * r
Vtr = V + beta * r

Ktr = sp.expand(Utr**2 - Vtr**2)
increment = sp.expand(Ktr - K)

expected_increment = sp.expand(
    2 * alpha * r * U
    - 2 * beta * r * V
    + alpha**2 * r**2
    - beta**2 * r**2
)

print("  T(alpha,beta):")
print("    U -> U + alpha*r")
print("    V -> V + beta*r")
print()
print("  Delta K =", sp.factor(increment))

cert(
    "translation increment",
    increment,
    expected_increment,
)

print()


# ============================================================================
# [4] V-PRESERVING TRANSLATION
# ============================================================================

print("[4] V-PRESERVING TRANSLATION")
print("-" * 78)

V_shift = sp.expand(Vtr - V)

print("  V' - V =", V_shift)

cert(
    "beta=0 preserves V",
    V_shift.subs(beta, 0),
    0,
)

print()


# ============================================================================
# [5] PURE S CHANNEL
# ============================================================================

print("[5] PURE S CHANNEL")
print("-" * 78)

linear_channel = sp.expand(
    (
        increment
        - alpha**2 * r**2
        + beta**2 * r**2
    ) / r**2
)

expected_linear_channel = sp.expand(
    alpha * S - beta * 2 * V / r
)

print("  normalized linear channel =")
print("   ", sp.factor(linear_channel))

cert(
    "general linear channel",
    linear_channel,
    expected_linear_channel,
)

cert(
    "beta=0 gives alpha*S",
    linear_channel.subs(beta, 0),
    alpha * S,
)

print()


# ============================================================================
# [6] NORMALIZED U TRANSLATION ORBIT
# ============================================================================

print("[6] NORMALIZED U-TRANSLATION ORBIT")
print("-" * 78)

Uh = U + h * r

Kh = sp.expand(
    Uh**2 - V**2
)

Lh = sp.expand(
    Kh / r**2
)

print("  K_h =", sp.factor(Kh))
print("  L_h =", sp.factor(Lh))

cert(
    "L_h = N+h*S+h^2",
    Lh,
    N + h * S + h**2,
)

print()


# ============================================================================
# [7] FIRST AND SECOND FINITE DIFFERENCE
# ============================================================================

print("[7] FINITE DIFFERENCE STRUCTURE")
print("-" * 78)

D1 = sp.expand(
    Lh.subs(h, h + 1) - Lh
)

D2 = sp.expand(
    Lh.subs(h, h + 2)
    - 2 * Lh.subs(h, h + 1)
    + Lh
)

print("  Delta L =", sp.factor(D1))
print("  Delta^2 L =", sp.factor(D2))

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

print()


# ============================================================================
# [8] V-FIXED AFFINE CLASSIFICATION
# ============================================================================

print("[8] V-FIXED AFFINE CLASSIFICATION")
print("-" * 78)

A_U, B_U, shift = sp.symbols("A_U B_U shift")

U_fixed = A_U * U + B_U * V + shift

K_fixed = sp.expand(
    U_fixed**2 - V**2
)

poly_fixed = sp.Poly(
    sp.expand(K_fixed - K),
    U,
    V,
)

f_U2 = sp.expand(poly_fixed.coeff_monomial(U**2))
f_UV = sp.expand(poly_fixed.coeff_monomial(U * V))
f_V2 = sp.expand(poly_fixed.coeff_monomial(V**2))

print("  residual U^2 =", f_U2)
print("  residual UV  =", f_UV)
print("  residual V^2 =", f_V2)

# Positive orientation branch
cert(
    "V-fixed + orientation",
    f_U2.subs({A_U: 1, B_U: 0}),
    0,
)

cert(
    "V-fixed + orientation UV",
    f_UV.subs({A_U: 1, B_U: 0}),
    0,
)

# Negative orientation branch
cert(
    "V-fixed - orientation",
    f_U2.subs({A_U: -1, B_U: 0}),
    0,
)

cert(
    "V-fixed - orientation UV",
    f_UV.subs({A_U: -1, B_U: 0}),
    0,
)

print()


# ============================================================================
# [9] ORIENTATION-PRESERVING / REVERSING ORBITS
# ============================================================================

print("[9] TWO V-FIXED OPERATOR BRANCHES")
print("-" * 78)

K_plus = sp.expand(
    (U + alpha * r)**2 - V**2
)

K_minus = sp.expand(
    (-U + alpha * r)**2 - V**2
)

L_plus = sp.expand(K_plus / r**2)
L_minus = sp.expand(K_minus / r**2)

print("  + branch:")
print("    U' = U + alpha*r")
print("    L =", sp.factor(L_plus))

print()
print("  - branch:")
print("    U' = -U + alpha*r")
print("    L =", sp.factor(L_minus))

cert(
    "positive branch",
    L_plus,
    N + alpha * S + alpha**2,
)

cert(
    "negative branch",
    L_minus,
    N - alpha * S + alpha**2,
)

print()


# ============================================================================
# [10] TRANSLATION GROUP COMPOSITION
# ============================================================================

print("[10] TRANSLATION GROUP COMPOSITION")
print("-" * 78)

T_alpha_U = U + alpha * r
T_beta_U = U + beta * r

composition = sp.expand(
    T_alpha_U.subs(U, T_beta_U)
)

direct = sp.expand(
    U + (alpha + beta) * r
)

print("  T_alpha(T_beta(U)) =", composition)
print("  T_(alpha+beta)(U) =", direct)

cert(
    "translation composition",
    composition,
    direct,
)

print()


# ============================================================================
# [11] REFLECTION INVOLUTION
# ============================================================================

print("[11] REFLECTION / ORIENTATION REVERSAL")
print("-" * 78)

R_U = -U
RR_U = sp.expand(-R_U)

print("  R(U)  =", R_U)
print("  R^2(U) =", RR_U)

cert(
    "reflection involution",
    RR_U,
    U,
)

print()


# ============================================================================
# [12] REFLECTION CONJUGACY
# ============================================================================

print("[12] REFLECTION CONJUGACY")
print("-" * 78)

# R T_h R(U) = U-h*r
R_T_R = sp.expand(
    -((-U) + h * r)
)

T_minus_h = sp.expand(
    U - h * r
)

print("  R T_h R(U) =", R_T_R)
print("  T_-h(U)   =", T_minus_h)

cert(
    "R T_h R = T_-h",
    R_T_R,
    T_minus_h,
)

print()


# ============================================================================
# [13] ACTION ON FACTOR COORDINATES
# ============================================================================

print("[13] FACTOR-COORDINATE ACTION")
print("-" * 78)

P = sp.expand((U - V) / r)
Q = sp.expand((U + V) / r)

P_h = sp.expand((U + h * r - V) / r)
Q_h = sp.expand((U + h * r + V) / r)

print("  P =", P)
print("  Q =", Q)
print("  T_h(P) =", P_h)
print("  T_h(Q) =", Q_h)

cert(
    "P translation",
    P_h,
    P + h,
)

cert(
    "Q translation",
    Q_h,
    Q + h,
)

print()


# ============================================================================
# [14] KAPPA QUADRATIC
# ============================================================================

print("[14] KAPPA QUADRATIC")
print("-" * 78)

Qpoly = sp.expand(
    z**2 - S*z + N
)

translated_poly = sp.expand(
    Qpoly.subs(z, z - 1)
)

print("  Q(z) =", sp.factor(Qpoly))
print("  Q(z-1) =", sp.factor(translated_poly))

cert(
    "Q linear coefficient",
    Qpoly.coeff(z, 1),
    -S,
)

cert(
    "Q constant coefficient",
    Qpoly.subs(z, 0),
    N,
)

cert(
    "Q(z-1) linear coefficient",
    translated_poly.coeff(z, 1),
    -(S + 2),
)

cert(
    "Q(z-1) constant coefficient",
    translated_poly.subs(z, 0),
    N + S + 1,
)

print()


# ============================================================================
# [15] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print("[15] FOUR-BASE CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi = sp.expand(
    Qpoly * translated_poly
)

target_chi = sp.expand(
    (z - P)
    * (z - Q)
    * (z - P - 1)
    * (z - Q - 1)
)

cert(
    "Q(z)Q(z-1) = four-base spectrum",
    chi,
    target_chi,
)

print("  chi(z) =", sp.factor(chi))
print()


# ============================================================================
# [16] DISCRIMINANT INVARIANT
# ============================================================================

print("[16] DISCRIMINANT INVARIANT")
print("-" * 78)

disc_Q = sp.expand(
    sp.discriminant(Qpoly, z)
)

print("  discriminant =", sp.factor(disc_Q))

cert(
    "Q discriminant = Delta",
    disc_Q,
    Delta,
)

print()


# ============================================================================
# [17] HISTORICAL x,y COORDINATES
# ============================================================================

print("[17] HISTORICAL x,y CONJUGACY")
print("-" * 78)

Uxy = 2 * y
Vxy = 2 * x - 3

Uxy_shift = sp.expand(
    2 * (y + h * r / 2)
)

Vxy_shift = sp.expand(
    2 * x - 3
)

print("  U = 2y")
print("  V = 2x-3")

cert(
    "U shift",
    Uxy_shift,
    Uxy + h * r,
)

cert(
    "V invariant",
    Vxy_shift,
    Vxy,
)

print()


# ============================================================================
# [18] HISTORICAL v KERNEL
# ============================================================================

print("[18] HISTORICAL v KERNEL")
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
    "v translation law",
    dv,
    h * r * y + h**2 * r**2 / 4,
)

print()


# ============================================================================
# [19] NORMALIZED HISTORICAL ORBIT
# ============================================================================

print("[19] NORMALIZED HISTORICAL ORBIT")
print("-" * 78)

Lhist = sp.expand(
    (4 * v - 1) / r**2
)

Lhist_h = sp.expand(
    (4 * v_h - 1) / r**2
)

N_xy = sp.expand(
    N.subs({
        U: Uxy,
        V: Vxy,
    })
)

S_xy = sp.expand(
    S.subs({
        U: Uxy,
    })
)

expected_Lhist_h = sp.expand(
    N_xy + h * S_xy + h**2
)

print("  L_hist =", sp.factor(Lhist))
print("  T_h(L_hist) =", sp.factor(Lhist_h))

cert(
    "historical orbit = N+hS+h^2",
    Lhist_h,
    expected_Lhist_h,
)

print()


# ============================================================================
# [20] HISTORICAL SECOND DIFFERENCE
# ============================================================================

print("[20] HISTORICAL DISCRETE CURVATURE")
print("-" * 78)

L0 = sp.expand(
    (4 * v - 1) / r**2
)

L1 = sp.expand(
    (4 * v.subs(y, y + r / 2) - 1) / r**2
)

L2 = sp.expand(
    (4 * v.subs(y, y + r) - 1) / r**2
)

hist_curvature = sp.expand(
    L2 - 2 * L1 + L0
)

print("  L0 =", sp.factor(L0))
print("  L1 =", sp.factor(L1))
print("  L2 =", sp.factor(L2))
print("  curvature =", sp.factor(hist_curvature))

cert(
    "historical second difference = 2",
    hist_curvature,
    2,
)

print()


# ============================================================================
# [21] OPERATOR RECONSTRUCTION OF U AND S
# ============================================================================

print("[21] OPERATOR RECONSTRUCTION")
print("-" * 78)

K0 = sp.expand(U**2 - V**2)
K1 = sp.expand((U + r)**2 - V**2)

U_from_operator = sp.expand(
    (K1 - K0 - r**2) / (2 * r)
)

S_from_operator = sp.expand(
    2 * U_from_operator / r
)

cert(
    "U reconstructed from K0,K1",
    U_from_operator,
    U,
)

cert(
    "S reconstructed from K0,K1",
    S_from_operator,
    S,
)

print()


# ============================================================================
# [22] OPERATOR INVARIANT / MOVING COORDINATE
# ============================================================================

print("[22] INVARIANT / MOVING DECOMPOSITION")
print("-" * 78)

V_after = V
U_after = U + h * r

Delta_after = sp.expand(
    4 * V_after**2 / r**2
)

S_after = sp.expand(
    2 * U_after / r
)

print("  invariant coordinate V' =", V_after)
print("  moving coordinate U'    =", U_after)

cert(
    "Delta invariant",
    Delta_after,
    Delta,
)

cert(
    "S affine response",
    S_after,
    S + 2 * h,
)

print()


# ============================================================================
# [23] SEARCHABLE UPSTREAM OPERATOR SIGNATURE
# ============================================================================

print("[23] SEARCHABLE UPSTREAM OPERATOR SIGNATURE")
print("-" * 78)

print("""
  For an unknown operator T, the desired normal form is:

      T(J) = J

      T(H) - H = c

      T(K) - K = 2*c*H + c^2

      T^2(K) - 2*T(K) + K = 2*c^2

  with the historical realization:

      J = V
      H = U
      K = U^2 - V^2
      c = r.

  Equivalently:

      U = (T(K)-K-r^2)/(2r)

      S = 2U/r

      Delta = 4V^2/r^2

      Q(z) = z^2-Sz+N.
""")

print()


# ============================================================================
# [24] FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 529 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)
print("""
The affine normal form is now fully explicit.

The next experiment should take an ACTUAL transformation from
the original homogeneous-layer construction and compare it with:

    T(U,V) = (U+r,V)

using the operator signature:

    invariant J:
        T(J)=J

    translated coordinate H:
        T(H)-H=c

    quadratic observable K:
        T(K)-K=2*c*H+c^2

    curvature:
        T^2(K)-2T(K)+K=2*c^2.

The strongest possible identification is:

    J  -> V
    H  -> U
    K  -> U^2-V^2.

Do not insert p or q during this test.

The goal is to identify the transformation itself before
factor recovery.
""")