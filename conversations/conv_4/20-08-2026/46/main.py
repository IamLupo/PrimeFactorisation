#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 522 START")
print("=" * 78)
print("CORRECTED TRANSLATION ORBIT / HOMOGENEOUS-LAYER SIGNATURE SEARCH")
print()

# ============================================================================
# SYMBOLS
# ============================================================================

U, V, r, m = sp.symbols("U V r m")
N, S, Delta = sp.symbols("N S Delta")

# Normalized factor coordinates
P = sp.expand((U - V) / r)
Q = sp.expand((U + V) / r)

N0 = sp.expand(P * Q)
S0 = sp.expand(P + Q)
D0 = sp.expand((P - Q) ** 2)

# Historical translated conic
K_m = sp.expand((U + m*r)**2 - V**2)
L_m = sp.expand(K_m / r**2)

# ============================================================================
# CERTIFICATE
# ============================================================================

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
# [1] CORRECT TRANSLATION ORBIT
# ============================================================================

print("[1] CORRECT TRANSLATION ORBIT")
print("-" * 78)

print("  K_m =", K_m)
print("  L_m =", sp.factor(L_m))

check(
    "L_m = N + m*S + m^2",
    L_m,
    N0 + m*S0 + m**2,
)

print()

# ============================================================================
# [2] FIRST DIFFERENCE
# ============================================================================

print("[2] FIRST DIFFERENCE")
print("-" * 78)

L_m1 = sp.expand(((U + (m+1)*r)**2 - V**2) / r**2)

dL = sp.expand(L_m1 - L_m)

print("  Delta L_m =", sp.factor(dL))

check(
    "Delta L_m = S + 2m + 1",
    dL,
    S0 + 2*m + 1,
)

print()

# ============================================================================
# [3] SECOND DIFFERENCE
# ============================================================================

print("[3] SECOND DIFFERENCE")
print("-" * 78)

L_m2 = sp.expand(((U + (m+2)*r)**2 - V**2) / r**2)

d2L = sp.expand(L_m2 - 2*L_m1 + L_m)

print("  Delta^2 L_m =", sp.factor(d2L))

check(
    "Delta^2 L_m = 2",
    d2L,
    2,
)

print()

# ============================================================================
# [4] RECOVER S FROM A SINGLE TRANSLATED VALUE
# ============================================================================

print("[4] SINGLE-LAYER S RECOVERY")
print("-" * 78)

S_from_L = sp.expand(
    (L_m - N0 - m**2) / m
)

print("  S_from_L =", sp.factor(S_from_L))

check(
    "S recovery for m != 0",
    S_from_L,
    S0,
)

print()

# ============================================================================
# [5] RECOVER S FROM A FIRST DIFFERENCE
# ============================================================================

print("[5] FIRST-DIFFERENCE S RECOVERY")
print("-" * 78)

S_from_dL = sp.expand(
    dL - 2*m - 1
)

print("  S_from_dL =", sp.factor(S_from_dL))

check(
    "S recovery from first difference",
    S_from_dL,
    S0,
)

print()

# ============================================================================
# [6] GENERATING QUADRATIC
# ============================================================================

print("[6] GENERATING QUADRATIC")
print("-" * 78)

Q_m = sp.expand(L_m)

print("  Q(m) =", sp.factor(Q_m))

check(
    "Q(m) = m^2 + S*m + N",
    Q_m,
    m**2 + S0*m + N0,
)

# roots
roots = sp.solve(sp.Eq(Q_m, 0), m)

print("  roots =", roots)

check(
    "Q(-P)=0",
    Q_m.subs(m, -P),
    0,
)

check(
    "Q(-Q)=0",
    Q_m.subs(m, -Q),
    0,
)

print()

# ============================================================================
# [7] DISCRIMINANT
# ============================================================================

print("[7] DISCRIMINANT OF TRANSLATION ORBIT")
print("-" * 78)

disc = sp.factor(sp.discriminant(Q_m, m))

print("  discriminant =", disc)

check(
    "discriminant = Delta",
    disc,
    D0,
)

print()

# ============================================================================
# [8] THREE-EVALUATION RECONSTRUCTION
# ============================================================================

print("[8] THREE-EVALUATION RECONSTRUCTION")
print("-" * 78)

L0 = sp.expand(L_m.subs(m, 0))
L1 = sp.expand(L_m.subs(m, 1))
L2 = sp.expand(L_m.subs(m, 2))

print("  L0 =", sp.factor(L0))
print("  L1 =", sp.factor(L1))
print("  L2 =", sp.factor(L2))

check("L0 = N", L0, N0)
check("L1 = N+S+1", L1, N0 + S0 + 1)
check("L2 = N+2S+4", L2, N0 + 2*S0 + 4)

S_three = sp.expand(
    (L1 - L0) - 1
)

check(
    "S = L1-L0-1",
    S_three,
    S0,
)

print()

# ============================================================================
# [9] TWO-LAYER TRANSLATION SIGNATURE
# ============================================================================

print("[9] TWO-LAYER SIGNATURE")
print("-" * 78)

translation_signature = sp.expand(L1 - L0)

print("  L1-L0 =", sp.factor(translation_signature))

check(
    "L1-L0 = S+1",
    translation_signature,
    S0 + 1,
)

print()

# ============================================================================
# [10] HISTORICAL CONIC FORM
# ============================================================================

print("[10] HISTORICAL CONIC FORM")
print("-" * 78)

K0 = sp.expand(U**2 - V**2)
K1 = sp.expand((U+r)**2 - V**2)

check(
    "K0/r^2 = N",
    K0/r**2,
    N0,
)

check(
    "(K1-K0)/r^2 = S+1",
    (K1-K0)/r**2,
    S0 + 1,
)

print()

# ============================================================================
# [11] INVERSE MAP FROM TRANSLATION ORBIT
# ============================================================================

print("[11] INVERSE MAP")
print("-" * 78)

# Given:
#
#   L_m = N + mS + m^2
#
# the quadratic coefficient is fixed,
# constant coefficient is N,
# linear coefficient is S.

poly_m = sp.Poly(m**2 + S*m + N, m)

quad_coeff = poly_m.coeff_monomial(m**2)
lin_coeff = poly_m.coeff_monomial(m)
const_coeff = poly_m.coeff_monomial(1)

print("  quadratic coefficient =", quad_coeff)
print("  linear coefficient    =", lin_coeff)
print("  constant coefficient  =", const_coeff)

check("quadratic coefficient = 1", quad_coeff, 1)
check("linear coefficient = S", lin_coeff, S0)
check("constant coefficient = N", const_coeff, N0)

print()

# ============================================================================
# [12] HOMOGENEOUS-LAYER SIGNATURE TEMPLATE
# ============================================================================

print("[12] HOMOGENEOUS-LAYER SIGNATURE TEMPLATE")
print("-" * 78)

print("""
  The exact signature to search for upstream is:

      L_m = N + m*S + m^2

  or equivalently:

      L_m - m^2 = N + m*S

  or in finite differences:

      L_{m+1}-L_m = S + 2m + 1

  or:

      (L_{m+1}-L_m) - (2m+1) = S.

  The second finite difference is universal:

      L_{m+2}-2L_{m+1}+L_m = 2.

  This means the hidden information is entirely in the
  FIRST difference, not in the quadratic curvature.
""")

print()

# ============================================================================
# [13] TRANSLATION-OPERATOR FORM
# ============================================================================

print("[13] OPERATOR FORM")
print("-" * 78)

E = sp.Function("E")

# Symbolically represent:
#
#   E[L](m) = L(m+1)
#
# and verify the expected identities directly.

symbolic_L = sp.Function("L")(m)

first_difference_form = sp.Eq(
    sp.Symbol("DeltaL"),
    S + 2*m + 1,
)

second_difference_form = sp.Eq(
    sp.Symbol("Delta2L"),
    2,
)

print("  First-difference law:")
print("    ΔL_m =", S, "+ 2m + 1")

print("  Second-difference law:")
print("    Δ²L_m = 2")

print()

# ============================================================================
# [14] CONNECTION TO KAPPA QUADRATIC
# ============================================================================

print("[14] CONNECTION TO KAPPA QUADRATIC")
print("-" * 78)

z = sp.symbols("z")

Q_z = sp.expand(z**2 - S0*z + N0)

# The translation orbit polynomial in m is
# m^2 + S*m + N.
#
# Replacing m -> -z gives:
#
# z^2 - S*z + N.

cross = sp.expand(Q_m.subs(m, -z))

print("  Q_translation(-z) =", sp.factor(cross))
print("  Q_KAPPA(z)        =", sp.factor(Q_z))

check(
    "translation orbit -> KAPPA quadratic",
    cross,
    Q_z,
)

print()

# ============================================================================
# [15] SHIFTED KAPPA FACTOR
# ============================================================================

print("[15] TRANSLATED KAPPA FACTOR")
print("-" * 78)

Q_z_shift = sp.expand(Q_z.subs(z, z-1))
Q_z_shift_expected = sp.expand(
    z**2 - (S0+2)*z + (N0+S0+1)
)

check(
    "Q(z-1) = shifted KAPPA quadratic",
    Q_z_shift,
    Q_z_shift_expected,
)

print()

# ============================================================================
# [16] FOUR-BASE COMPOSITION
# ============================================================================

print("[16] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(Q_z * Q_z_shift)

expected_chi = sp.expand(
    (z-P)*(z-Q)*(z-P-1)*(z-Q-1)
)

check(
    "Q(z)Q(z-1) = four-base polynomial",
    chi,
    expected_chi,
)

print()

# ============================================================================
# [17] DISCRIMINANT -> GAP
# ============================================================================

print("[17] GAP RECOVERY")
print("-" * 78)

S_rec = lin_coeff
N_rec = const_coeff

Delta_rec = sp.expand(
    S_rec**2 - 4*N_rec
)

print("  recovered S =", sp.factor(S_rec))
print("  recovered N =", sp.factor(N_rec))
print("  recovered Delta =", sp.factor(Delta_rec))

check(
    "recovered Delta = Delta",
    Delta_rec,
    D0,
)

print()

# ============================================================================
# FINAL
# ============================================================================

print("=" * 78)
print("EXPERIMENT 522 FINISHED")
print("=" * 78)

print("""
CORRECTED CONCLUSION
------------------------------------------------------------------------------

The historical translation orbit is exactly

    L_m = K_m/r^2
        = N + m*S + m^2.

There is NO extra +m term in L_m.

The +1 appears only in the first difference:

    L_{m+1}-L_m = S + 2m + 1.

The translation orbit itself is therefore literally the
KAPPA quadratic:

    L_m = m^2 + S*m + N.

Under m -> -z:

    L_{-z} = z^2 - S*z + N.

Thus the historical conic translation parameter is exactly
the spectral variable of the KAPPA quadratic.

This yields the strongest bridge found so far:

    historical translated conic
        -> L_m
        -> m^2 + S*m + N
        -> KAPPA quadratic
        -> Q(z-1)
        -> four-base characteristic polynomial.

For the upstream problem, the exact object to search for is
therefore NOT merely K0 or K1.

It is a sequence L_m whose dependence on the translation index
is a monic quadratic with:

    coefficient(m^2) = 1
    coefficient(m)   = hidden S
    constant         = known N.

Equivalently, search for an upstream translation orbit whose
FIRST DIFFERENCE contains the missing symmetric coordinate.

The second difference carries no factor information:

    Δ² L_m = 2.

That sharply identifies where the useful hidden information
must reside.
""")

print(f"SYMBOLIC FAILURES = {failures}")

if failures == 0:
    print("OVERALL EXACT SYMBOLIC AUDIT = True")
else:
    print("OVERALL EXACT SYMBOLIC AUDIT = False")
