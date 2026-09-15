#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 528 START")
print("=" * 78)
print("AFFINE CONIC OPERATOR CLASSIFICATION / KAPPA S-CHANNEL")
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

U, V = sp.symbols("U V")
r = sp.symbols("r", nonzero=True)

a, b, c, d = sp.symbols("a b c d")
alpha, beta = sp.symbols("alpha beta")

t = sp.symbols("t")

K = sp.expand(U**2 - V**2)

S = sp.expand(2 * U / r)
N = sp.expand(K / r**2)
Delta = sp.expand(4 * V**2 / r**2)

print()


# ============================================================================
# [1] BASE CONIC
# ============================================================================

print("[1] BASE QUADRATIC FORM")
print("-" * 78)

print("  K(U,V) = U^2 - V^2")
print("  N =", N)
print("  S =", S)
print("  Delta =", Delta)

check("K definition", K, U**2 - V**2)

print()


# ============================================================================
# [2] GENERAL AFFINE OPERATOR
# ============================================================================

print("[2] GENERAL AFFINE OPERATOR")
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
# [3] HOMOGENEOUS QUADRATIC PART
# ============================================================================

print("[3] HOMOGENEOUS QUADRATIC PART")
print("-" * 78)

K1_hom = sp.expand(
    K1.subs({alpha: 0, beta: 0})
)

print("  homogeneous part =", sp.factor(K1_hom))

target_matrix_conditions = [
    sp.expand(a**2 - c**2 - 1),
    sp.expand(b**2 - d**2 + 1),
    sp.expand(a*b - c*d),
]

print()
print("  K-preserving linear conditions:")
print("    a^2-c^2 = 1")
print("    b^2-d^2 = -1")
print("    ab-cd = 0")

print()


# ============================================================================
# [4] SPECIAL KAPPA TRANSLATION
# ============================================================================

print("[4] PURE U TRANSLATION")
print("-" * 78)

U_T = U + t * r
V_T = V

K_T = sp.expand(U_T**2 - V_T**2)

print("  T_t:")
print("    U -> U+t*r")
print("    V -> V")

print("  K_T-K =", sp.factor(K_T - K))

check(
    "pure U translation increment",
    K_T - K,
    2 * t * r * U + t**2 * r**2,
)

print()


# ============================================================================
# [5] NORMALIZED ORBIT
# ============================================================================

print("[5] NORMALIZED U-TRANSLATION ORBIT")
print("-" * 78)

L_t = sp.expand(K_T / r**2)

print("  L_t =", sp.factor(L_t))

check(
    "L_t = N+tS+t^2",
    L_t,
    N + t * S + t**2,
)

dL = sp.expand(
    L_t.subs(t, t + 1) - L_t
)

d2L = sp.expand(
    L_t.subs(t, t + 2)
    - 2 * L_t.subs(t, t + 1)
    + L_t
)

print("  ΔL_t =", sp.factor(dL))
print("  Δ²L_t =", sp.factor(d2L))

check(
    "first difference",
    dL,
    S + 2 * t + 1,
)

check(
    "second difference",
    d2L,
    2,
)

print()


# ============================================================================
# [6] OPERATOR CHARACTERIZATION
# ============================================================================

print("[6] OPERATOR CHARACTERIZATION OF THE S-CHANNEL")
print("-" * 78)

linear_information = sp.expand(
    (K_T - K) / r**2 - t**2
)

print("  normalized increment minus t^2 =")
print("   ", sp.factor(linear_information))

check(
    "linear channel = t*S",
    linear_information,
    t * S,
)

print("""
  Hence the directional derivative of the conic along the
  U-translation is the S-channel:

      (T_t K - K)/r^2 - t^2 = t*S.

  No V-dependent term survives.
""")

print()


# ============================================================================
# [7] GENERAL TRANSLATION DIRECTION
# ============================================================================

print("[7] GENERAL TRANSLATION DIRECTION")
print("-" * 78)

U_ab = U + alpha * r
V_ab = V + beta * r

K_ab = sp.expand(U_ab**2 - V_ab**2)

normalized_increment = sp.expand(
    (K_ab - K) / r**2
)

print("  normalized increment =")
print("   ", sp.factor(normalized_increment))

expected_general = sp.expand(
    alpha * S
    - beta * (2 * V / r)
    + alpha**2
    - beta**2
)

check(
    "general normalized translation law",
    normalized_increment,
    expected_general,
)

print()


# ============================================================================
# [8] PURE-SUM CONDITION
# ============================================================================

print("[8] PURE-SUM CONDITION")
print("-" * 78)

linear_general = sp.expand(
    normalized_increment - alpha**2 + beta**2
)

print("  linear channel =", sp.factor(linear_general))

print("""
  General linear channel:

      alpha*S - beta*(2V/r)

  For the channel to contain S but no independent V-term,
  the necessary condition is:

      beta = 0.

  Therefore the distinguished KAPPA direction is unique
  among translations up to scale.
""")

check(
    "V-free condition beta=0",
    linear_general.subs(beta, 0),
    alpha * S,
)

print()


# ============================================================================
# [9] GAP CHANNEL
# ============================================================================

print("[9] GAP CHANNEL")
print("-" * 78)

gap = sp.expand(2 * V)

gap_translated = sp.expand(
    2 * (V + beta * r)
)

print("  original gap =", gap)
print("  translated gap =", gap_translated)

check(
    "gap change",
    gap_translated - gap,
    2 * beta * r,
)

print("""
  Therefore:

      beta = 0

  simultaneously means:

      V is invariant,
      Delta is invariant,
      |P-Q| is invariant.

  The KAPPA translation is consequently the unique
  translation direction that changes S while preserving
  the gap coordinate.
""")

print()


# ============================================================================
# [10] DISCRIMINANT UNDER THE DISTINGUISHED OPERATOR
# ============================================================================

print("[10] DISCRIMINANT RESPONSE")
print("-" * 78)

Delta_t = sp.expand(
    4 * V_T**2 / r**2
)

check(
    "Delta invariant under U translation",
    Delta_t,
    Delta,
)

print()


# ============================================================================
# [11] CONIC / KAPPA QUADRATIC CONNECTION
# ============================================================================

print("[11] GENERATING POLYNOMIAL CONNECTION")
print("-" * 78)

z = sp.symbols("z")

Q = sp.expand(
    z**2 - S*z + N
)

L_minus_z = sp.expand(
    L_t.subs(t, -z)
)

print("  Q(z) =", sp.factor(Q))
print("  L(-z) =", sp.factor(L_minus_z))

check(
    "L(-z) = Q(z)",
    L_minus_z,
    Q,
)

print()


# ============================================================================
# [12] ROOT STRUCTURE
# ============================================================================

print("[12] ROOT STRUCTURE")
print("-" * 78)

Q_roots = sp.factor(
    Q.subs(z, z)
)

print("  Q(z) =", Q_roots)

check(
    "Q factorization",
    Q,
    (z - (U - V) / r) *
    (z - (U + V) / r),
)

print()


# ============================================================================
# [13] DISCRIMINANT FROM OPERATOR ORBIT
# ============================================================================

print("[13] DISCRIMINANT FROM THE ORBIT")
print("-" * 78)

Q_disc = sp.expand(
    S**2 - 4 * N
)

check(
    "Q discriminant = Delta",
    Q_disc,
    Delta,
)

print("  discriminant =", sp.factor(Q_disc))

print()


# ============================================================================
# [14] OPERATOR COMPOSITION
# ============================================================================

print("[14] TRANSLATION COMPOSITION")
print("-" * 78)

def T_U(expr, h):
    return sp.expand(
        expr.subs(U, U + h * r)
    )


T1T2 = T_U(
    T_U(K, alpha),
    beta,
)

Tsum = T_U(
    K,
    alpha + beta,
)

print("  T_alpha T_beta(K) =", sp.factor(T1T2))
print("  T_(alpha+beta)(K) =", sp.factor(Tsum))

check(
    "translation composition law",
    T1T2,
    Tsum,
)

print()


# ============================================================================
# [15] DISCRETE GENERATOR
# ============================================================================

print("[15] DISCRETE GENERATOR")
print("-" * 78)

K0 = K
K1 = T_U(K, 1)
K2 = T_U(K, 2)

first = sp.expand(
    K1 - K0
)

second = sp.expand(
    K2 - 2*K1 + K0
)

print("  K0 =", sp.factor(K0))
print("  K1 =", sp.factor(K1))
print("  K2 =", sp.factor(K2))

check(
    "first U difference",
    first,
    2 * U * r + r**2,
)

check(
    "second U difference",
    second,
    2 * r**2,
)

print()


# ============================================================================
# [16] NORMALIZED FIRST DIFFERENCE AS OPERATOR EIGENCHANNEL
# ============================================================================

print("[16] FIRST-DIFFERENCE OPERATOR CHANNEL")
print("-" * 78)

L0 = sp.expand(K0 / r**2)
L1 = sp.expand(K1 / r**2)

channel = sp.expand(
    L1 - L0 - 1
)

print("  L1-L0-1 =", sp.factor(channel))

check(
    "one-step S extraction",
    channel,
    S,
)

print()


# ============================================================================
# [17] SECOND DIFFERENCE AS UNIVERSAL CURVATURE")
# ============================================================================

print("[17] UNIVERSAL CURVATURE")
print("-" * 78)

curvature = sp.expand(
    K2 / r**2
    - 2 * K1 / r**2
    + K0 / r**2
)

print("  normalized curvature =", curvature)

check(
    "curvature independent of U,V",
    curvature,
    2,
)

print()


# ============================================================================
# [18] AFFINE OPERATOR WITH PURE S CHANNEL
# ============================================================================

print("[18] AFFINE PURE-S CLASSIFICATION")
print("-" * 78)

# General affine translation only:
#
#   U' = U + alpha r
#   V' = V + beta r
#
# A pure S channel requires beta=0.

pure_s_solution = sp.solve(
    sp.Eq(beta, 0),
    beta,
    dict=True,
)

print("  condition:", pure_s_solution)

print("""
  Translation class:

      (U,V) -> (U+alpha*r,V)

  produces

      ΔN = alpha*S + alpha^2.

  The linear term is exactly the symmetric coordinate.
""")

check(
    "pure S classification",
    linear_general.subs(beta, 0),
    alpha * S,
)

print()


# ============================================================================
# [19] HISTORICAL x,y FORM
# ============================================================================

print("[19] HISTORICAL x,y FORM")
print("-" * 78)

x, y = sp.symbols("x y")

Uxy = 2 * y
Vxy = 2 * x - 3

Kxy = sp.expand(
    Uxy**2 - Vxy**2
)

Kxy_shifted = sp.expand(
    Kxy.subs(
        y,
        y + alpha * r / 2
    )
)

KU_expected_xy = sp.expand(
    (Uxy + alpha * r)**2 - Vxy**2
)

check(
    "U translation = y translation",
    Kxy_shifted,
    KU_expected_xy,
)

print("  U -> U+alpha*r corresponds to:")
print("    y -> y + alpha*r/2")
print("    x -> x")

print()


# ============================================================================
# [20] HISTORICAL CONIC VALUE
# ============================================================================

print("[20] HISTORICAL v-KERNEL RESPONSE")
print("-" * 78)

v = sp.expand(
    y**2 - x**2 + 3*x - 2
)

v_shift = sp.expand(
    v.subs(
        y,
        y + alpha * r / 2
    )
)

v_delta = sp.expand(
    v_shift - v
)

print("  v =", v)
print("  T(v)-v =", sp.factor(v_delta))

check(
    "v translation law",
    v_delta,
    alpha * r * y + alpha**2 * r**2 / 4,
)

print("""
  Thus the historical v-kernel itself carries a quadratic
  translation orbit:

      T_alpha(v)-v
        = alpha*r*y + alpha^2*r^2/4.

  The coefficient of the linear term is y, which is exactly
  the historical representation of the missing symmetric
  coordinate.
""")

print()


# ============================================================================
# [21] OPERATOR-LEVEL BRIDGE SIGNATURE
# ============================================================================

print("[21] OPERATOR-LEVEL BRIDGE SIGNATURE")
print("-" * 78)

# For normalized L=(4v-1)/r^2:
L_hist = sp.expand(
    (4 * v - 1) / r**2
)

L_hist_shift = sp.expand(
    L_hist.subs(
        y,
        y + alpha * r / 2
    )
)

hist_increment = sp.expand(
    L_hist_shift - L_hist
)

print("  L_hist =", sp.factor(L_hist))
print("  T(L_hist)-L_hist =")
print("   ", sp.factor(hist_increment))

check(
    "historical normalized translation",
    hist_increment,
    alpha * (4*y/r) + alpha**2,
)

check(
    "historical coefficient = S",
    (hist_increment - alpha**2) / alpha,
    4*y/r,
)

print()


# ============================================================================
# [22] OPERATOR CONJUGACY TEST
# ============================================================================

print("[22] OPERATOR CONJUGACY")
print("-" * 78)

# Historical:
#
#    U = 2y
#
# Therefore:
#
#    y -> y + r/2
#
# is exactly:
#
#    U -> U+r.

U_shift_from_y = sp.expand(
    2 * (y + r / 2)
)

check(
    "y -> y+r/2 equals U -> U+r",
    U_shift_from_y,
    Uxy + r,
)

print()


# ============================================================================
# [23] OPERATOR INVARIANT / MOVING PART DECOMPOSITION
# ============================================================================

print("[23] INVARIANT / MOVING PART DECOMPOSITION")
print("-" * 78)

K_U_orbit = sp.expand(
    (U + t*r)**2 - V**2
)

moving_part = sp.expand(
    K_U_orbit - (U**2 - V**2)
)

print("  invariant coordinate: V")
print("  moving coordinate: U")
print("  moving part =", sp.factor(moving_part))

check(
    "moving part",
    moving_part,
    2*t*r*U + t**2*r**2,
)

print("""
  Under the distinguished operator:

      V is invariant,
      U moves affinely.

  Therefore:

      Delta = 4V^2/r^2

  is an operator invariant, while

      S = 2U/r

  is the affine coordinate of the orbit.
""")

print()


# ============================================================================
# [24] FINAL STRUCTURAL CERTIFICATE
# ============================================================================

print("[24] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print("""
  Distinguished operator:

      T_h(U,V) = (U+h*r,V).

  Its normalized conic orbit is:

      L_h
        = (U+h*r)^2/r^2 - V^2/r^2
        = N + h*S + h^2.

  Therefore:

      L_{h+1}-L_h-(2h+1) = S

      L_{h+2}-2L_{h+1}+L_h = 2.

  The operator simultaneously satisfies:

      V invariant
      Delta invariant
      S affine
      N = L_0
      Q(z) = L_{-z}

  and hence:

      Q(z)=z^2-Sz+N.

  This makes the KAPPA quadratic an operator orbit,
  rather than an independently introduced polynomial.
""")

print()


# ============================================================================
# FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 528 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT SYMBOLIC AUDIT = {failures == 0}")

print("""
NEXT RESEARCH TARGET
------------------------------------------------------------------------------

The affine classification says that the relevant upstream object
must realize the operator

    T(U,V) = (U+r,V)

or an algebraically conjugate operator.

Therefore the next experiment should search the ORIGINAL
homogeneous construction for an operator T satisfying:

    T^2(K) - 2T(K) + K = 2*r^2

and

    T(K)-K-r^2

is linear in the hidden symmetric coordinate.

More importantly, search for an operator where one coordinate
is invariant and the other is translated:

    invariant coordinate  -> unchanged
    moving coordinate     -> + constant.

The historical coordinates give the exact template:

    U = 2y
    V = 2x-3

so the target transformation is:

    y -> y + r/2
    x -> x.

If an existing homogeneous-layer transformation can be
conjugated to this operation, that would provide the first
operator-level upstream bridge without reconstructing p,q,
S, Delta, or the KAPPA sequence first.
""")
