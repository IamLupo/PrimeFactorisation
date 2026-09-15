#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 523 START")
print("=" * 78)
print("HOMOGENEOUS-LAYER TRANSLATION SIGNATURE / QUADRATIC ORBIT SEARCH")
print()

# ============================================================================
# SYMBOLS
# ============================================================================

U, V, r, m, j = sp.symbols("U V r m j")
z = sp.symbols("z")

# Correct dependent definitions.
N = sp.expand((U**2 - V**2) / r**2)
S = sp.expand(2*U / r)
Delta = sp.expand(4*V**2 / r**2)

P = sp.expand((U - V) / r)
Q = sp.expand((U + V) / r)

failures = 0


def check(label, expr, expected=0):
    global failures

    diff = sp.factor(sp.together(sp.expand(expr - expected)))

    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# [1] EXACT HISTORICAL ORBIT
# ============================================================================

print("[1] EXACT HISTORICAL TRANSLATION ORBIT")
print("-" * 78)

L = sp.expand(
    ((U - V + m*r) * (U + V + m*r)) / r**2
)

print("  L(m) =", sp.factor(L))

check(
    "L(m) = m^2 + S*m + N",
    L,
    m**2 + S*m + N,
)

print()


# ============================================================================
# [2] COEFFICIENT EXTRACTION
# ============================================================================

print("[2] ORBIT COEFFICIENT EXTRACTION")
print("-" * 78)

poly_L = sp.Poly(L, m)

c2 = sp.expand(poly_L.coeff_monomial(m**2))
c1 = sp.expand(poly_L.coeff_monomial(m))
c0 = sp.expand(poly_L.coeff_monomial(1))

print("  coefficient m^2 =", c2)
print("  coefficient m   =", c1)
print("  coefficient 1   =", c0)

check("quadratic coefficient = 1", c2, 1)
check("linear coefficient = S", c1, S)
check("constant coefficient = N", c0, N)

print()


# ============================================================================
# [3] FIRST DIFFERENCE
# ============================================================================

print("[3] FIRST-DIFFERENCE SIGNATURE")
print("-" * 78)

L1 = sp.expand(L.subs(m, m + 1))

dL = sp.expand(L1 - L)

print("  Delta L(m) =", sp.factor(dL))

check(
    "Delta L(m) = S + 2m + 1",
    dL,
    S + 2*m + 1,
)

print()


# ============================================================================
# [4] SECOND DIFFERENCE
# ============================================================================

print("[4] SECOND-DIFFERENCE SIGNATURE")
print("-" * 78)

L2 = sp.expand(L.subs(m, m + 2))

d2L = sp.expand(L2 - 2*L1 + L)

print("  Delta^2 L(m) =", d2L)

check(
    "Delta^2 L(m) = 2",
    d2L,
    2,
)

print()


# ============================================================================
# [5] RECOVER S FROM TWO ADJACENT LAYERS
# ============================================================================

print("[5] TWO-LAYER S RECOVERY")
print("-" * 78)

S_from_layers = sp.expand(
    (L1 - L) - (2*m + 1)
)

print("  S_from_layers =", sp.factor(S_from_layers))

check(
    "S = (L(m+1)-L(m))-(2m+1)",
    S_from_layers,
    S,
)

print()


# ============================================================================
# [6] RECOVER N FROM ONE LAYER
# ============================================================================

print("[6] N RECOVERY")
print("-" * 78)

N_from_layer = sp.expand(
    L - m*S - m**2
)

print("  N_from_layer =", sp.factor(N_from_layer))

check(
    "N = L(m)-mS-m^2",
    N_from_layer,
    N,
)

print()


# ============================================================================
# [7] TRANSLATION ORBIT AS KAPPA QUADRATIC
# ============================================================================

print("[7] TRANSLATION ORBIT -> KAPPA QUADRATIC")
print("-" * 78)

Q_translation = sp.expand(L.subs(m, -z))

Q_kappa = sp.expand(
    (z - P) * (z - Q)
)

print("  Q_translation(z) =", sp.factor(Q_translation))
print("  Q_kappa(z)       =", sp.factor(Q_kappa))

check(
    "L(-z) = (z-P)(z-Q)",
    Q_translation,
    Q_kappa,
)

print()


# ============================================================================
# [8] DISCRIMINANT
# ============================================================================

print("[8] DISCRIMINANT")
print("-" * 78)

disc_L = sp.factor(
    sp.discriminant(sp.Poly(L, m))
)

print("  discriminant =", disc_L)
print("  target Delta =", Delta)

check(
    "disc(L) = Delta",
    disc_L,
    Delta,
)

print()


# ============================================================================
# [9] KAPPA +1 TRANSLATION
# ============================================================================

print("[9] KAPPA +1 TRANSLATION")
print("-" * 78)

Q1 = sp.expand(Q_kappa.subs(z, z - 1))

expected_Q1 = sp.expand(
    (z - (P + 1)) * (z - (Q + 1))
)

check(
    "Q(z-1) = (z-(P+1))(z-(Q+1))",
    Q1,
    expected_Q1,
)

print()


# ============================================================================
# [10] FOUR-BASE COMPOSITION
# ============================================================================

print("[10] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(Q_kappa * Q1)

expected_chi = sp.expand(
    (z-P) * (z-Q) * (z-P-1) * (z-Q-1)
)

check(
    "Q(z)Q(z-1) = four-base characteristic",
    chi,
    expected_chi,
)

print()


# ============================================================================
# [11] CONIC LAYERS K_m
# ============================================================================

print("[11] CONIC LAYER REPRESENTATION")
print("-" * 78)

K = sp.expand(
    (U + m*r)**2 - V**2
)

check(
    "K(m)/r^2 = L(m)",
    K / r**2,
    L,
)

print()


# ============================================================================
# [12] TRANSLATION OF THE CONIC INVARIANT
# ============================================================================

print("[12] CONIC TRANSLATION LAW")
print("-" * 78)

K0 = sp.expand(K.subs(m, 0))
K1 = sp.expand(K.subs(m, 1))

print("  K0 =", sp.factor(K0))
print("  K1 =", sp.factor(K1))

check(
    "K0/r^2 = N",
    K0 / r**2,
    N,
)

check(
    "(K1-K0)/r^2 = S+1",
    (K1 - K0) / r**2,
    S + 1,
)

print()


# ============================================================================
# [13] SPECTRAL COEFFICIENT FROM THE TWO-LAYER SYSTEM
# ============================================================================

print("[13] SPECTRAL LINEAR COEFFICIENT")
print("-" * 78)

a = sp.expand(-S)

a_from_layers = sp.expand(
    1 - (K1 - K0) / r**2
)

print("  a_from_layers =", sp.factor(a_from_layers))
print("  expected a    =", sp.factor(a))

check(
    "a = 1-(K1-K0)/r^2",
    a_from_layers,
    a,
)

print()


# ============================================================================
# [14] ABSTRACT LAYER SIGNATURE
# ============================================================================

print("[14] ABSTRACT UPSTREAM SIGNATURE")
print("-" * 78)

print("""
  A candidate upstream family L_m is sufficient if it satisfies

      L_m = m^2 + A*m + B

  with:

      A = hidden symmetric coordinate,
      B = N.

  Operationally this is equivalent to

      Delta L_m = A + 2m + 1

  and therefore

      A = Delta L_m - 2m - 1.

  The universal curvature condition is

      Delta^2 L_m = 2.

  Thus the homogeneous-layer search should look for a parameter
  m for which the layer observable has:

      second difference = constant 2,

  while the first difference retains the unknown linear coefficient.
""")

print()


# ============================================================================
# [15] OPERATOR TEST
# ============================================================================

print("[15] FINITE-DIFFERENCE OPERATOR")
print("-" * 78)

# Forward difference operator
def forward_difference(expr, var):
    return sp.expand(expr.subs(var, var + 1) - expr)


D1 = forward_difference(L, m)
D2 = forward_difference(D1, m)

check(
    "D(L) = S + 2m + 1",
    D1,
    S + 2*m + 1,
)

check(
    "D^2(L) = 2",
    D2,
    2,
)

print()


# ============================================================================
# [16] GENERAL QUADRATIC RECONSTRUCTION
# ============================================================================

print("[16] GENERAL THREE-LAYER RECONSTRUCTION")
print("-" * 78)

# Treat L0,L1,L2 as abstract observations.
ell0, ell1, ell2 = sp.symbols("L0 L1 L2")

A_from_three = sp.expand(
    ell1 - ell0 - 1
)

B_from_three = sp.expand(
    ell0
)

curvature = sp.expand(
    ell2 - 2*ell1 + ell0
)

print("  A =", A_from_three)
print("  B =", B_from_three)
print("  curvature =", curvature)

print()


# ============================================================================
# [17] CRITICAL UPSTREAM TEST
# ============================================================================

print("[17] CRITICAL UPSTREAM TEST")
print("-" * 78)

print("""
  The actual remaining problem is now:

      HOMOGENEOUS_LAYER(m)
          |
          v
      candidate L_m
          |
          +--> L_{m+1}-L_m
          |
          +--> subtract (2m+1)
          |
          v
          hidden S

  A successful upstream object therefore does NOT need to
  reproduce P, Q, Delta, or the four KAPPA bases directly.

  It only needs to expose a translated quadratic orbit.

  The three strongest signatures are:

      Signature A:
          L_m = m^2 + A*m + B

      Signature B:
          ΔL_m - (2m+1) = constant

      Signature C:
          Δ²L_m = 2

  Among these, Signature B is the most informative because
  it directly exposes the missing linear coefficient.
""")

print()


# ============================================================================
# [18] FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 523 FINISHED")
print("=" * 78)

print()
print(f"SYMBOLIC FAILURES = {failures}")

if failures == 0:
    print("OVERALL EXACT SYMBOLIC AUDIT = True")
else:
    print("OVERALL EXACT SYMBOLIC AUDIT = False")

print("""
NEXT RESEARCH TARGET
------------------------------------------------------------------------------

Do not search more identities among N, S, Delta, F_n, or the
historical factor coordinates.

The next upstream task is:

    enumerate the ACTUAL homogeneous-layer formulas already
    developed in the earlier construction and test whether
    one of their layer indices acts as m in

        L_m = m^2 + A*m + B.

The decisive test is not fitting values.

It is exact symbolic operator structure:

    L_{m+1}-L_m-(2m+1) = A
    L_{m+2}-2L_{m+1}+L_m = 2.

If such an object exists upstream and B=N, then the KAPPA
quadratic is generated directly by the homogeneous layers:

    L_m -> m^2 + S*m + N
        -> z^2-Sz+N
        -> factor pair.

That is the bridge we should attack next.
""")
