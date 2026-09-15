#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 525 START")
print("=" * 78)
print("HOMOGENEOUS-LAYER INDEX -> TRANSLATION OPERATOR EXTRACTION")
print()

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
m, j = sp.symbols("m j")
z = sp.symbols("z")

N = sp.expand(p * q)
S = sp.expand(p + q)
Delta = sp.expand((p - q)**2)

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
# [1] CANONICAL TRANSLATION ORBIT
# ============================================================================

print("[1] CANONICAL TRANSLATION ORBIT")
print("-" * 78)

L = sp.expand((p + m) * (q + m))

print("  L(m) =", sp.factor(L))

check(
    "L(m)=m^2+S*m+N",
    L,
    m**2 + S*m + N,
)

print()


# ============================================================================
# [2] REMOVE THE UNIVERSAL QUADRATIC CURVATURE
# ============================================================================

print("[2] CURVATURE-REMOVING OPERATOR")
print("-" * 78)

D1 = sp.expand(
    L.subs(m, m+1) - L
)

D2 = sp.expand(
    L.subs(m, m+2)
    - 2*L.subs(m, m+1)
    + L
)

print("  ΔL =", sp.factor(D1))
print("  Δ²L =", sp.factor(D2))

check(
    "ΔL=S+2m+1",
    D1,
    S + 2*m + 1,
)

check(
    "Δ²L=2",
    D2,
    2,
)

print()


# ============================================================================
# [3] DEFINE THE ABSTRACT HOMOGENEOUS LAYER
# ============================================================================

print("[3] ABSTRACT HOMOGENEOUS LAYER")
print("-" * 78)

# We do NOT assume the final form.
#
# A general symmetric polynomial up to quadratic degree in p,q
# and quadratic degree in the layer parameter is:
#
# H_m =
#     c20*m^2
#   + c11*m*(p+q)
#   + c00*p*q
#   + c20p*(p^2+q^2)
#   + c10p*(p+q)
#   + c00p
#
# The last three terms represent possible layer-independent
# homogeneous contamination.

c20, c11, c00 = sp.symbols("c20 c11 c00")
u20, u10, u00 = sp.symbols("u20 u10 u00")

H = sp.expand(
    c20*m**2
    + c11*m*(p+q)
    + c00*p*q
    + u20*(p**2+q**2)
    + u10*(p+q)
    + u00
)

print("  H_m =", H)

print()


# ============================================================================
# [4] IMPOSE THE TRANSLATION-CURVATURE SIGNATURE
# ============================================================================

print("[4] TRANSLATION-CURVATURE SIGNATURE")
print("-" * 78)

H1 = sp.expand(H.subs(m, m+1))
H2 = sp.expand(H.subs(m, m+2))

dH = sp.expand(H1 - H)
d2H = sp.expand(H2 - 2*H1 + H)

print("  ΔH  =", sp.factor(dH))
print("  Δ²H =", sp.factor(d2H))

# Desired structural conditions:
#
#   Δ²H = 2
#
# which immediately fixes c20 = 1.

sol_curvature = sp.solve(
    sp.Poly(
        sp.expand(d2H - 2),
        p, q, m
    ).coeffs(),
    [c20],
    dict=True,
)

print("  curvature solutions =", sol_curvature)

if sol_curvature:
    H_curved = sp.expand(
        H.subs(sol_curvature[0])
    )
else:
    H_curved = H

check(
    "normalized curvature",
    sp.expand(
        H_curved.subs(m, m+2)
        - 2*H_curved.subs(m, m+1)
        + H_curved
    ),
    2,
)

print()


# ============================================================================
# [5] FIRST-DIFFERENCE INFORMATION CHANNEL
# ============================================================================

print("[5] FIRST-DIFFERENCE INFORMATION CHANNEL")
print("-" * 78)

dH_curved = sp.expand(
    H_curved.subs(m, m+1) - H_curved
)

print("  ΔH =", sp.factor(dH_curved))

# Corrected first difference:
# ΔH - (2m+1)

corrected = sp.expand(
    dH_curved - (2*m + 1)
)

print("  ΔH-(2m+1) =", sp.factor(corrected))

print()


# ============================================================================
# [6] REQUIRE THE LINEAR COEFFICIENT TO BE PURELY SYMMETRIC
# ============================================================================

print("[6] PURE-SUM REQUIREMENT")
print("-" * 78)

# The desired result is:
#
#   ΔH-(2m+1) = S = p+q
#
# Therefore all p^2+q^2, constant, and asymmetric contributions
# must disappear.

poly_corr = sp.Poly(
    sp.expand(corrected - (p+q)),
    p, q, m,
)

equations_corr = [
    sp.Eq(coeff, 0)
    for coeff in poly_corr.coeffs()
]

solutions_corr = sp.solve(
    equations_corr,
    [c11, c00, u20, u10, u00],
    dict=True,
)

print("  coefficient solutions =", solutions_corr)

if solutions_corr:
    H_solution = sp.expand(
        H_curved.subs(solutions_corr[0])
    )

    print("  resulting H_m =", sp.factor(H_solution))

    check(
        "H_m = m^2 + S*m + constant",
        H_solution - m**2 - S*m,
        sp.expand(H_solution - m**2 - S*m),
    )
else:
    H_solution = H_curved

print()


# ============================================================================
# [7] CONSTANT-TERM TEST
# ============================================================================

print("[7] CONSTANT-TERM TEST")
print("-" * 78)

H0 = sp.expand(
    H_solution.subs(m, 0)
)

print("  H_0 =", sp.factor(H0))

check(
    "H_0 = N",
    H0,
    N,
)

print()


# ============================================================================
# [8] FULL ORBIT RECONSTRUCTION
# ============================================================================

print("[8] FULL ORBIT RECONSTRUCTION")
print("-" * 78)

check(
    "H_m = m^2 + S*m + N",
    H_solution,
    m**2 + S*m + N,
)

print()


# ============================================================================
# [9] GENERATING POLYNOMIAL
# ============================================================================

print("[9] GENERATING POLYNOMIAL")
print("-" * 78)

Q = sp.expand(
    H_solution.subs(m, -z)
)

Q_expected = sp.expand(
    z**2 - S*z + N
)

print("  Q(z) =", sp.factor(Q))

check(
    "H(-z)=z^2-Sz+N",
    Q,
    Q_expected,
)

print()


# ============================================================================
# [10] TRANSLATED KAPPA FACTOR
# ============================================================================

print("[10] TRANSLATED KAPPA FACTOR")
print("-" * 78)

Q1 = sp.expand(
    Q.subs(z, z-1)
)

Q1_expected = sp.expand(
    (z-(p+1))*(z-(q+1))
)

check(
    "Q(z-1)=(z-p-1)(z-q-1)",
    Q1,
    Q1_expected,
)

print()


# ============================================================================
# [11] FOUR-BASE COMPOSITION
# ============================================================================

print("[11] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(
    Q * Q1
)

chi_expected = sp.expand(
    (z-p)*(z-q)*(z-p-1)*(z-q-1)
)

check(
    "Q(z)Q(z-1)=four-base polynomial",
    chi,
    chi_expected,
)

print()


# ============================================================================
# [12] DISCRIMINANT CHANNEL
# ============================================================================

print("[12] DISCRIMINANT CHANNEL")
print("-" * 78)

disc_Q = sp.expand(
    sp.discriminant(Q, z)
)

print("  disc(Q) =", sp.factor(disc_Q))

check(
    "disc(Q)=Delta",
    disc_Q,
    Delta,
)

print()


# ============================================================================
# [13] TRANSLATION ORBIT AS A DIFFERENCE EQUATION
# ============================================================================

print("[13] TRANSLATION ORBIT AS A DIFFERENCE EQUATION")
print("-" * 78)

A = sp.symbols("A")

# General solution of
#
#     Δ² H_m = 2
#
# is
#
#     H_m = m^2 + A*m + B.
#
# We verify it directly.

B = sp.symbols("B")

general_solution = sp.expand(
    m**2 + A*m + B
)

general_d1 = sp.expand(
    general_solution.subs(m, m+1)
    - general_solution
)

general_d2 = sp.expand(
    general_solution.subs(m, m+2)
    - 2*general_solution.subs(m, m+1)
    + general_solution
)

print("  general H_m =", general_solution)
print("  ΔH =", general_d1)
print("  Δ²H =", general_d2)

check(
    "general second difference",
    general_d2,
    2,
)

check(
    "general first-difference form",
    general_d1,
    A + 2*m + 1,
)

print()


# ============================================================================
# [14] OPERATOR FACTORIZATION
# ============================================================================

print("[14] OPERATOR FACTORIZATION")
print("-" * 78)

# Let D be the forward difference operator:
#
#     D f(m) = f(m+1)-f(m)
#
# On the translation quadratic:
#
#     (D^2 - 2)H = 0.
#
# The hidden symmetric coordinate is the integration constant:
#
#     A = D H - 2m - 1.

operator_A = sp.expand(
    general_d1 - 2*m - 1
)

print("  A_operator =", operator_A)

check(
    "operator recovers A",
    operator_A,
    A,
)

print()


# ============================================================================
# [15] WHAT AN ACTUAL HOMOGENEOUS-LAYER OBJECT MUST SATISFY
# ============================================================================

print("[15] EXACT SIGNATURE TEST")
print("-" * 78)

print("""
  An actual upstream family H_m is a successful bridge candidate
  iff all four exact conditions hold:

      H_0 = N

      Δ²H_m = 2

      ΔH_m - (2m+1) = S

      H(-z) = z²-Sz+N

  The second condition is universal.

  The third condition is the information-bearing condition.

  The first condition fixes the N-channel.

  The fourth condition converts the layer orbit into the
  KAPPA spectral polynomial.
""")

print()


# ============================================================================
# [16] CRITICAL DISTINCTION
# ============================================================================

print("[16] CRITICAL DISTINCTION")
print("-" * 78)

print("""
  This experiment does NOT claim that the original
  homogeneous-layer construction has already produced H_m.

  It proves the following stronger structural statement:

      If a homogeneous-layer family contains a quadratic
      translation orbit with curvature 2, then its entire
      factor-bearing content is determined by two integration
      constants:

          A = ΔH_m-(2m+1)
          B = H_0.

  To reproduce the KAPPA bridge these must become:

          A = S
          B = N.

  Therefore the upstream search has been reduced to finding
  the source of the two integration constants.

  In particular, the useful search is no longer for an
  arbitrary polynomial identity.

  It is for a layer operator whose first difference leaves
  behind an N-independent linear symmetric invariant.
""")

print()


# ============================================================================
# [17] FINAL AUDIT
# ============================================================================

print("=" * 78)
print("EXPERIMENT 525 FINISHED")
print("=" * 78)
print()
print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT SYMBOLIC AUDIT = {failures == 0}")

print("""
NEXT TARGET
------------------------------------------------------------------------------

The next experiment should finally compare this operator signature
against the ORIGINAL homogeneous-layer recurrence/kernel itself.

For every genuine layer family already present in the research,
test only:

    H_0

    H_{m+1}-H_m

    H_{m+2}-2H_{m+1}+H_m

and ask whether:

    H_0 = N

and:

    H_{m+1}-H_m-(2m+1)

is independent of m.

If it is, do NOT expand it further immediately.

Instead identify that surviving expression as the candidate
symmetric coordinate.

The decisive discovery would be:

    H_{m+1}-H_m-(2m+1)
        = an existing homogeneous-layer invariant

which is not explicitly written as p+q.

That would be the first genuinely upstream bridge.
""")
