#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
EXPERIMENT 519
==============================================================================
HOMOGENEOUS-LAYER TRANSLATION / CONIC OPERATOR BRIDGE

Goal
----
The modern KAPPA recurrence has characteristic polynomial

    chi_F(z)
      = (z-P)(z-Q)(z-(P+1))(z-(Q+1))

with

    P = (U-V)/r
    Q = (U+V)/r.

The historical conic is

    U^2 - V^2 = K

and the +1 shift corresponds to

    U -> U+r,   V -> V.

This experiment studies that translation structurally.

No numerical datasets.
No factor-pair enumeration.
No continued fractions.
No interpolation / data fitting.

The experiment asks:

    1. Can the four-base characteristic polynomial be expressed
       as two translated quadratic factors?

    2. Is the second quadratic obtained from the first by a
       pure U-translation?

    3. What coefficient transformation does U -> U+r induce?

    4. Can the recurrence coefficients be written entirely as
       translation invariants?

    5. What exact operator signature would an upstream
       homogeneous-layer construction need to reproduce?
==============================================================================
"""

import sympy as sp


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def cert(name, expr):
    expr = sp.factor(sp.expand(expr))
    ok = (expr == 0)
    print(f"  {name}")
    print(f"    difference = {expr}")
    print(f"    PASS = {ok}")
    print()
    return ok


def section(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# -----------------------------------------------------------------------------
# Symbols
# -----------------------------------------------------------------------------

z, xi = sp.symbols("z xi")
U, V, r = sp.symbols("U V r", nonzero=True)

P = (U - V) / r
Q = (U + V) / r

N = sp.expand(P * Q)
S = sp.expand(P + Q)
Delta = sp.expand((P - Q) ** 2)

P1 = P + 1
Q1 = Q + 1


# -----------------------------------------------------------------------------
# [1] HISTORICAL CONIC COORDINATES
# -----------------------------------------------------------------------------

section("[1] HISTORICAL / NORMALIZED COORDINATES")

print("  P =", P)
print("  Q =", Q)

print("  N =", sp.factor(N))
print("  S =", sp.factor(S))
print("  Delta =", sp.factor(Delta))

cert(
    "PQ = (U^2-V^2)/r^2",
    P * Q - (U**2 - V**2) / r**2,
)

cert(
    "P+Q = 2U/r",
    P + Q - 2 * U / r,
)

cert(
    "P-Q = -2V/r",
    P - Q + 2 * V / r,
)


# -----------------------------------------------------------------------------
# [2] TARGET KAPPA QUADRATIC FACTORS
# -----------------------------------------------------------------------------

section("[2] TARGET TWO-QUADRATIC FACTORIZATION")

Q0 = sp.expand((z - P) * (z - Q))
Q1 = sp.expand((z - P1) * (z - Q1))

print("  Q0(z) =", sp.factor(Q0))
print("  Q1(z) =", sp.factor(Q1))
print()

cert(
    "Q0 = z^2 - S*z + N",
    Q0 - (z**2 - S * z + N),
)

cert(
    "Q1 = z^2 - (S+2)z + (N+S+1)",
    Q1 - (z**2 - (S + 2) * z + (N + S + 1)),
)


# -----------------------------------------------------------------------------
# [3] TRANSLATION U -> U+r
# -----------------------------------------------------------------------------

section("[3] PURE U-TRANSLATION")

def translate_U(expr):
    return sp.expand(expr.subs(U, U + r))


Q0_translated = sp.factor(translate_U(Q0))

print("  T_U[r](Q0) =", Q0_translated)
print()

cert(
    "T_U[r](Q0) = Q1",
    Q0_translated - Q1,
)


# -----------------------------------------------------------------------------
# [4] COEFFICIENT TRANSLATION LAW
# -----------------------------------------------------------------------------

section("[4] COEFFICIENT TRANSLATION LAW")

a0 = -S
b0 = N

a1 = -(S + 2)
b1 = N + S + 1

print("  Q0(z) = z^2 + a0*z + b0")
print("  Q1(z) = z^2 + a1*z + b1")
print()

print("  a0 =", a0)
print("  b0 =", b0)
print("  a1 =", a1)
print("  b1 =", b1)
print()

cert(
    "linear coefficient translation",
    a1 - a0 + 2,
)

cert(
    "constant coefficient translation",
    b1 - b0 - (S + 1),
)


# -----------------------------------------------------------------------------
# [5] TRANSLATION AS A COEFFICIENT OPERATOR
# -----------------------------------------------------------------------------

section("[5] COEFFICIENT OPERATOR")

# For a generic monic quadratic z^2 + a z + b:
#
# Translating the roots by +1:
#
#   (z-(P+1))(z-(Q+1))
#
# gives
#
#   z^2 + (a-2)z + (b-a+1)
#
# because a = -(P+Q), b = PQ.

a, b = sp.symbols("a b")

generic_q = sp.expand(z**2 + a * z + b)

generic_shift = sp.expand(
    generic_q.subs(z, z - 1)
)

expected_shift = sp.expand(
    z**2 + (a - 2) * z + (b - a + 1)
)

print("  generic quadratic      =", generic_q)
print("  translated quadratic   =", sp.factor(generic_shift))
print("  expected coefficient law:")
print("    a' = a - 2")
print("    b' = b - a + 1")
print()

cert(
    "generic quadratic translation law",
    generic_shift - expected_shift,
)


# -----------------------------------------------------------------------------
# [6] APPLY THE OPERATOR TO Q0
# -----------------------------------------------------------------------------

section("[6] APPLY TRANSLATION OPERATOR TO Q0")

Q0_coeff_a = -S
Q0_coeff_b = N

translated_a = sp.expand(Q0_coeff_a - 2)
translated_b = sp.expand(Q0_coeff_b - Q0_coeff_a + 1)

print("  original coefficients:")
print("    a =", Q0_coeff_a)
print("    b =", Q0_coeff_b)
print()

print("  translated coefficients:")
print("    a' =", translated_a)
print("    b' =", translated_b)
print()

cert(
    "a' = -(S+2)",
    translated_a + (S + 2),
)

cert(
    "b' = N+S+1",
    translated_b - (N + S + 1),
)


# -----------------------------------------------------------------------------
# [7] FOURTH-ORDER CHARACTERISTIC POLYNOMIAL
# -----------------------------------------------------------------------------

section("[7] FOUR-BASE CHARACTERISTIC POLYNOMIAL")

chi = sp.expand(Q0 * Q1)

print("  chi(z) =", sp.factor(chi))
print()

chi_expected = sp.expand(
    (z**2 - S * z + N)
    *
    (z**2 - (S + 2) * z + (N + S + 1))
)

cert(
    "chi = Q0*Q1",
    chi - chi_expected,
)


# -----------------------------------------------------------------------------
# [8] RECURRENCE COEFFICIENTS
# -----------------------------------------------------------------------------

section("[8] RECURRENCE COEFFICIENTS")

poly = sp.Poly(sp.expand(chi), z)

coeffs = poly.all_coeffs()

a1_target = sp.expand(-coeffs[1])
a2_target = sp.expand(coeffs[2])
a3_target = sp.expand(-coeffs[3])
a4_target = sp.expand(coeffs[4])

print("  chi(z) = z^4 - a1*z^3 + a2*z^2 - a3*z + a4")
print()
print("  a1 =", sp.factor(a1_target))
print("  a2 =", sp.factor(a2_target))
print("  a3 =", sp.factor(a3_target))
print("  a4 =", sp.factor(a4_target))
print()

cert(
    "a1 target",
    a1_target - 2 * (S + 1),
)

cert(
    "a2 target",
    a2_target - (2 * N + S**2 + 3 * S + 1),
)

cert(
    "a3 target",
    a3_target - (S + 1) * (2 * N + S),
)

cert(
    "a4 target",
    a4_target - N * (N + S + 1),
)


# -----------------------------------------------------------------------------
# [9] RECURRENCE COEFFICIENTS AS TRANSLATION INVARIANTS
# -----------------------------------------------------------------------------

section("[9] TRANSLATION INVARIANT STRUCTURE")

# The key quantities are:
#
#   a1 = 2(S+1)
#   a2 = N + (N+S+1) + S(S+2)
#
# Observe:
#
#   Q0 constant term = N
#   Q1 constant term = M1 = N+S+1
#
# and
#
#   Q0 linear coefficient = -S
#   Q1 linear coefficient = -(S+2)

M1 = sp.expand(N + S + 1)

print("  M1 =", sp.factor(M1))
print()

cert(
    "M1 = (P+1)(Q+1)",
    M1 - P1 * Q1,
)

cert(
    "M1-N = S+1",
    M1 - N - S - 1,
)

cert(
    "a1 = 2*(S+1)",
    a1_target - 2 * (M1 - N),
)


# -----------------------------------------------------------------------------
# [10] RECURRENCE COEFFICIENTS FROM THE TWO QUADRATICS
# -----------------------------------------------------------------------------

section("[10] COMPOSITION LAW")

# Let
#
#   Q0 = z^2 + a*z + b
#   Q1 = z^2 + c*z + d.
#
# Their product has:
#
#   z^3 coefficient = a+c
#   z^2 coefficient = b+d+ac
#
# For the shifted pair:
#
#   a = -S
#   c = -(S+2)
#   b = N
#   d = N+S+1.

qa = -S
qb = N
qc = -(S + 2)
qd = N + S + 1

composition_a1 = sp.expand(-(qa + qc))
composition_a2 = sp.expand(qb + qd + qa * qc)

print("  From quadratic composition:")
print("    a1 =", sp.factor(composition_a1))
print("    a2 =", sp.factor(composition_a2))
print()

cert(
    "composition a1",
    composition_a1 - a1_target,
)

cert(
    "composition a2",
    composition_a2 - a2_target,
)


# -----------------------------------------------------------------------------
# [11] CENTERING REVISITED
# -----------------------------------------------------------------------------

section("[11] CENTERING OF THE TWO-QUADRATIC SYSTEM")

center = (S + 1) / 2

centered = sp.expand(
    chi.subs(z, xi + center)
)

print("  center =", center)
print()
print("  centered chi =", sp.factor(centered))
print()

centered_expected = sp.expand(
    ((xi + sp.Rational(1, 2))**2 - Delta / 4)
    *
    ((xi - sp.Rational(1, 2))**2 - Delta / 4)
)

cert(
    "centered quartic",
    centered - centered_expected,
)


# -----------------------------------------------------------------------------
# [12] CENTERED COEFFICIENT -> DELTA
# -----------------------------------------------------------------------------

section("[12] CENTERED COEFFICIENT -> GAP")

centered_poly = sp.Poly(centered, xi)
centered_coeffs = centered_poly.all_coeffs()

c2 = sp.expand(centered_coeffs[2])

print("  xi^4 coefficient =", centered_coeffs[0])
print("  xi^3 coefficient =", centered_coeffs[1])
print("  xi^2 coefficient =", c2)
print()

cert(
    "centered xi^3 vanishes",
    centered_coeffs[1],
)

cert(
    "centered xi^2 coefficient",
    c2 + Delta / 2 + sp.Rational(1, 2),
)


# -----------------------------------------------------------------------------
# [13] TRANSLATION-DECOMPOSED OBSERVABLES
# -----------------------------------------------------------------------------

section("[13] TRANSLATION-DECOMPOSED OBSERVABLES")

# Define:
#
#   L0 = PQ
#   L1 = (P+1)(Q+1)
#
#   C0 = P+Q
#   C1 = (P+1)+(Q+1)
#
# Then:
#
#   L1-L0 = C0+1
#   C1-C0 = 2
#
# This is the minimal algebraic translation signature.

L0 = sp.expand(P * Q)
L1 = sp.expand(P1 * Q1)

C0 = sp.expand(P + Q)
C1 = sp.expand(P1 + Q1)

print("  L0 = PQ")
print("  L1 = (P+1)(Q+1)")
print("  C0 = P+Q")
print("  C1 = (P+1)+(Q+1)")
print()

print("  L1-L0 =", sp.factor(L1 - L0))
print("  C1-C0 =", sp.factor(C1 - C0))
print()

cert(
    "L1-L0 = C0+1",
    L1 - L0 - C0 - 1,
)

cert(
    "C1-C0 = 2",
    C1 - C0 - 2,
)


# -----------------------------------------------------------------------------
# [14] CONIC TRANSLATION
# -----------------------------------------------------------------------------

section("[14] HISTORICAL CONIC TRANSLATION")

K0 = sp.expand(U**2 - V**2)
K1 = sp.expand((U + r)**2 - V**2)

print("  K0 =", K0)
print("  K1 =", K1)
print()

print("  K1-K0 =", sp.factor(K1 - K0))
print()

cert(
    "K1-K0 = r*(2U+r)",
    K1 - K0 - r * (2 * U + r),
)

cert(
    "K0/r^2 = N",
    K0 / r**2 - N,
)

cert(
    "K1/r^2 = M1",
    K1 / r**2 - M1,
)


# -----------------------------------------------------------------------------
# [15] BRIDGE SIGNATURE
# -----------------------------------------------------------------------------

section("[15] TRANSLATION BRIDGE SIGNATURE")

bridge_1 = sp.expand((K1 - K0) / r**2)
bridge_2 = sp.expand(K1 / r**2 - K0 / r**2)
bridge_3 = sp.expand(bridge_1 - (S + 1))

print("  (K1-K0)/r^2 =", sp.factor(bridge_1))
print("  K1/r^2-K0/r^2 =", sp.factor(bridge_2))
print()

cert(
    "(K1-K0)/r^2 = S+1",
    bridge_3,
)


# -----------------------------------------------------------------------------
# [16] WHAT AN UPSTREAM HOMOGENEOUS LAYER MUST REPRODUCE
# -----------------------------------------------------------------------------

section("[16] MINIMAL UPSTREAM SIGNATURE")

print(
"""
  The downstream KAPPA spectrum is completely determined by two
  quadratic factors:

      Q0(z) = z^2 - S*z + N
      Q1(z) = z^2 - (S+2)*z + (N+S+1)

  The second factor is not independent.

  It is obtained from the first by the universal coefficient operator

      (a,b) -> (a-2, b-a+1)

  for a monic quadratic z^2+a*z+b.

  Equivalently, at the root level:

      {P,Q} -> {P+1,Q+1}.

  In historical coordinates this is exactly

      U -> U+r
      V -> V.

  Therefore an upstream homogeneous-layer construction does NOT
  need to reproduce all four KAPPA bases independently.

  It is sufficient for it to produce a quadratic object

      Q(z) = z^2 + a(N)*z + b(N)

  together with an intrinsic translation operator

      T(Q)(z) = Q(z-1).

  The desired bridge would then be:

      homogeneous layer
          -> Q(z)
          -> T(Q)(z)
          -> Q*T(Q)
          -> a1,a2
          -> S, Delta.

  The critical obstruction is immediate:

      a(N) = -S
      b(N) = N.

  Thus b is already N-only, but a contains the missing symmetric
  coordinate.

  This isolates the upstream problem to one coefficient:

      N-only homogeneous layer
          -> coefficient carrying S

  rather than eight independent F-values.
"""
)


# -----------------------------------------------------------------------------
# [17] FINAL SYMBOLIC AUDIT
# -----------------------------------------------------------------------------

section("[17] FINAL EXACT AUDIT")

checks = []

checks.append(cert(
    "quadratic factor Q0",
    Q0 - (z**2 - S*z + N),
))

checks.append(cert(
    "quadratic factor Q1",
    Q1 - (z**2 - (S + 2)*z + (N + S + 1)),
))

checks.append(cert(
    "translation operator",
    translate_U(Q0) - Q1,
))

checks.append(cert(
    "generic coefficient translation",
    generic_shift - expected_shift,
))

checks.append(cert(
    "four-base characteristic",
    chi - Q0 * Q1,
))

checks.append(cert(
    "a1 recurrence coefficient",
    a1_target - 2*(S + 1),
))

checks.append(cert(
    "a2 recurrence coefficient",
    a2_target - (2*N + S**2 + 3*S + 1),
))

checks.append(cert(
    "centered quartic",
    centered - centered_expected,
))

checks.append(cert(
    "conic translation",
    K1 - K0 - r*(2*U + r),
))

checks.append(cert(
    "K0/r^2 = N",
    K0/r**2 - N,
))

checks.append(cert(
    "K1/r^2 = M1",
    K1/r**2 - M1,
))

checks.append(cert(
    "translation difference = S+1",
    (K1 - K0)/r**2 - (S + 1),
))

overall = all(checks)

print()
print("  OVERALL EXACT AUDIT =", overall)

print(
"""
==============================================================================
EXPERIMENT 519 CONCLUSION
==============================================================================

The four-base KAPPA spectrum is generated by a single quadratic
plus a universal translation operator.

At the root level:

    {P,Q}
       --(+1)-->
    {P+1,Q+1}

At the historical conic level:

    (U,V)
       --U->U+r-->
    (U+r,V)

At the polynomial level:

    Q(z)=z^2+a*z+b

    T(Q)(z)=Q(z-1)

    (a,b) -> (a-2, b-a+1).

For the KAPPA system:

    a = -S
    b = N

so the entire order-4 characteristic polynomial is generated
from the single missing coefficient S plus the already-known N.

This gives the new upstream target:

    homogeneous layer
        -> quadratic coefficient a
        -> S=-a
        -> translation Q(z)->Q(z-1)
        -> KAPPA characteristic polynomial.

The important unresolved question is therefore:

    Does the original N-only homogeneous-layer construction contain
    a quadratic/second-order object whose coefficient is mathematically
    equivalent to -S, even though S is not explicitly present?

That is the next bridge to attack.
==============================================================================
"""
)
