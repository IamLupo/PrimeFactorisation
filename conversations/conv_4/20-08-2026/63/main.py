#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 539 START")
print("=" * 78)
print("DISCRIMINANT-CHANNEL CLASSIFICATION / MINKOWSKI KAPPA RECOGNITION")
print()


# ============================================================================
# SYMBOLS
# ============================================================================

x1, x2 = sp.symbols("x1 x2")
t1, t2 = sp.symbols("t1 t2", nonzero=True)

A, B, C = sp.symbols("A B C")
D, E, F = sp.symbols("D E F")

h, z = sp.symbols("h z")

# Canonical moving/invariant coordinates will later be represented by
# linear forms H_raw and J_raw.


# ============================================================================
# HELPERS
# ============================================================================

failures = 0


def simp(expr):
    return sp.factor(sp.expand(sp.simplify(expr)))


def cert(label, expr, expected=0):
    global failures

    diff = simp(expr - expected)

    if isinstance(diff, sp.MatrixBase):
        ok = all(simp(v) == 0 for v in diff)
    else:
        ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# ============================================================================
# 1. GENERAL QUADRATIC FORM
# ============================================================================

print("[1] GENERAL QUADRATIC FORM")
print("-" * 78)

Q = sp.Matrix([
    [A, B / 2],
    [B / 2, C],
])

x = sp.Matrix([x1, x2])
t = sp.Matrix([t1, t2])

lvec = sp.Matrix([D, E])

K = sp.expand(
    (x.T * Q * x)[0]
    + (lvec.T * x)[0]
    + F
)

print("  quadratic matrix Q =")
sp.pprint(Q)

print()
print("  K(x) =")
sp.pprint(K)


# ============================================================================
# 2. TRANSLATION DATA
# ============================================================================

print()
print("[2] TRANSLATION DATA")
print("-" * 78)

qtt = sp.expand(
    (t.T * Q * t)[0]
)

moving_raw = sp.expand(
    2 * (Q * t).dot(x) + (lvec.T * t)[0]
)

print("  q(t,t) =")
sp.pprint(qtt)

print()
print("  first-difference coefficient =")
sp.pprint(moving_raw)

cert(
    "T(K)-K = 2 q(t,x) + q(t,t) + linear_translation",
    K.subs(
        {
            x1: x1 + t1,
            x2: x2 + t2,
        },
        simultaneous=True,
    ) - K,
    moving_raw,
)


# ============================================================================
# 3. NORMALIZED ORBIT
# ============================================================================

print()
print("[3] NORMALIZED QUADRATIC ORBIT")
print("-" * 78)

Kh = sp.expand(
    K.subs(
        {
            x1: x1 + h * t1,
            x2: x2 + h * t2,
        },
        simultaneous=True,
    )
)

Lh = sp.cancel(Kh / qtt)

Lpoly = sp.Poly(sp.expand(Lh), h)

c2 = simp(Lpoly.coeff_monomial(h**2))
c1 = simp(Lpoly.coeff_monomial(h))
c0 = simp(Lpoly.coeff_monomial(1))

print("  L(h) = c2*h² + c1*h + c0")

print()
print("  c2 =")
sp.pprint(c2)

print()
print("  c1 =")
sp.pprint(c1)

print()
print("  c0 =")
sp.pprint(c0)

cert("orbit is monic", c2, 1)


# ============================================================================
# 4. GENERATING SPECTRAL POLYNOMIAL
# ============================================================================

print()
print("[4] GENERATING POLYNOMIAL")
print("-" * 78)

Qorbit = sp.expand(
    Lh.subs(h, -z)
)

disc = sp.factor(
    sp.discriminant(Qorbit, z)
)

print("  Q(z)=L(-z)")
sp.pprint(Qorbit)

print()
print("  discriminant(Q) =")
sp.pprint(disc)


# ============================================================================
# 5. NATURAL MOVING / INVARIANT COORDINATES
# ============================================================================

print()
print("[5] NATURAL MOVING / INVARIANT COORDINATES")
print("-" * 78)

# The moving linear form is q(t,x).
#
# An invariant linear form annihilates t.
#
# In two dimensions:
#
#   J_raw = -t2*x1 + t1*x2

Hraw = sp.expand(
    (2 * Q * t).dot(x)
)

Jraw = sp.expand(
    -t2 * x1 + t1 * x2
)

print("  H_raw =")
sp.pprint(Hraw)

print()
print("  J_raw =")
sp.pprint(Jraw)

cert(
    "J_raw invariant",
    Jraw.subs(
        {
            x1: x1 + t1,
            x2: x2 + t2,
        },
        simultaneous=True,
    ) - Jraw,
    0,
)


# ============================================================================
# 6. CHANGE OF VARIABLES
# ============================================================================

print()
print("[6] CHANGE TO MOVING / TRANSVERSE BASIS")
print("-" * 78)

# Write:
#
#   Hraw = h1*x1 + h2*x2
#   Jraw = j1*x1 + j2*x2
#
# and solve for x1,x2 in terms of Hraw,Jraw.

h1 = sp.expand(Hraw).coeff(x1, 1)
h2 = sp.expand(Hraw).coeff(x2, 1)

j1 = -t2
j2 = t1

basis_matrix = sp.Matrix([
    [h1, h2],
    [j1, j2],
])

det_basis = sp.factor(
    basis_matrix.det()
)

print("  basis matrix =")
sp.pprint(basis_matrix)

print()
print("  basis determinant =")
sp.pprint(det_basis)

cert(
    "basis determinant = 2*q(t,t)",
    det_basis,
    2 * qtt,
)


# ============================================================================
# 7. INDEPENDENCE OF THE DISCRIMINANT FROM THE MOVING COORDINATE
# ============================================================================

print()
print("[7] DISCRIMINANT MOVING-CHANNEL TEST")
print("-" * 78)

# The discriminant is a polynomial in x1,x2.
#
# Replace x1,x2 by symbolic H,J coordinates by solving the linear basis.
Hsym, Jsym = sp.symbols("H J")

HJ = sp.Matrix([Hsym, Jsym])

sol_x = sp.simplify(
    basis_matrix.inv() * HJ
)

x1_HJ = sp.factor(sol_x[0])
x2_HJ = sp.factor(sol_x[1])

print("  x1(H,J) =")
sp.pprint(x1_HJ)

print()
print("  x2(H,J) =")
sp.pprint(x2_HJ)

disc_HJ = sp.factor(
    sp.together(
        disc.subs(
            {
                x1: x1_HJ,
                x2: x2_HJ,
            },
            simultaneous=True,
        )
    )
)

print()
print("  discriminant(H,J) =")
sp.pprint(disc_HJ)


# Polynomial in H.
disc_H_poly = sp.Poly(
    sp.factor(sp.together(disc_HJ)),
    Hsym,
)

print()
print("  coefficient of H² =")
sp.pprint(
    sp.factor(disc_H_poly.coeff_monomial(Hsym**2))
)

print()
print("  coefficient of H =")
sp.pprint(
    sp.factor(disc_H_poly.coeff_monomial(Hsym))
)

print()
print("  H-independent term =")
sp.pprint(
    sp.factor(disc_H_poly.coeff_monomial(1))
)


# ============================================================================
# 8. GENERAL CONDITION FOR H-INDEPENDENCE
# ============================================================================

print()
print("[8] H-INDEPENDENT DISCRIMINANT CONDITIONS")
print("-" * 78)

disc_H2 = sp.factor(
    disc_H_poly.coeff_monomial(Hsym**2)
)

disc_H1 = sp.factor(
    disc_H_poly.coeff_monomial(Hsym)
)

print("""
  A KAPPA-type discriminant must have:

      coefficient(H²) = 0
      coefficient(H)  = 0

  so that

      disc(Q) = function(J)

  rather than depending on the moving coordinate.
""")

print("  condition H² = 0:")
sp.pprint(disc_H2)

print()
print("  condition H¹ = 0:")
sp.pprint(disc_H1)


# ============================================================================
# 9. PURE QUADRATIC-FORM SPECIALIZATION
# ============================================================================

print()
print("[9] HOMOGENEOUS QUADRATIC FORM")
print("-" * 78)

# Remove the linear and constant terms.
Khom = sp.expand(
    (x.T * Q * x)[0]
)

Khom_h = sp.expand(
    Khom.subs(
        {
            x1: x1 + h * t1,
            x2: x2 + h * t2,
        },
        simultaneous=True,
    )
)

Lhom = sp.cancel(
    Khom_h / qtt
)

Qhom = sp.expand(
    Lhom.subs(h, -z)
)

dischom = sp.factor(
    sp.discriminant(Qhom, z)
)

print("  homogeneous orbit discriminant =")
sp.pprint(dischom)


# ============================================================================
# 10. QUADRATIC-FORM INVARIANT
# ============================================================================

print()
print("[10] TRANSVERSE QUADRATIC INVARIANT")
print("-" * 78)

# In two dimensions, the transverse invariant should be J_raw².
#
# Compute the scalar proportionality factor by comparing
# dischom with J_raw².

J2 = sp.expand(
    Jraw**2
)

ratio = sp.factor(
    sp.cancel(dischom / J2)
)

print("  J² =")
sp.pprint(J2)

print()
print("  discriminant / J² =")
sp.pprint(ratio)


# ============================================================================
# 11. MINKOWSKI SPECIALIZATION
# ============================================================================

print()
print("[11] MINKOWSKI SPECIALIZATION")
print("-" * 78)

U, V, c = sp.symbols("U V c", nonzero=True)

mink_subs = {
    A: 1,
    B: 0,
    C: -1,
    D: 0,
    E: 0,
    F: 0,
    x1: U,
    x2: V,
    t1: c,
    t2: 0,
}

Km = sp.expand(
    K.subs(mink_subs)
)

Lhm = sp.expand(
    (
        K.subs(
            mink_subs
            | {
                x1: U + h * c,
                x2: V,
            }
        )
    ) / c**2
)

Qm = sp.expand(
    Lhm.subs(h, -z)
)

discm = sp.factor(
    sp.discriminant(Qm, z)
)

print("  K =")
sp.pprint(Km)

print()
print("  L(h) =")
sp.pprint(Lhm)

print()
print("  Q(z) =")
sp.pprint(Qm)

print()
print("  discriminant =")
sp.pprint(discm)

cert(
    "Minkowski discriminant = 4V²/c²",
    discm,
    4 * V**2 / c**2,
)


# ============================================================================
# 12. HISTORICAL PULLBACK
# ============================================================================

print()
print("[12] HISTORICAL CONIC PULLBACK")
print("-" * 78)

xh, yh = sp.symbols("x_hist y_hist")

Uh = 2 * yh
Vh = 2 * xh - 3

vh = sp.expand(
    yh**2 - xh**2 + 3 * xh - 2
)

Khist = sp.expand(
    Uh**2 - Vh**2
)

cert(
    "historical K = 4v-1",
    Khist,
    4 * vh - 1,
)


# ============================================================================
# 13. HISTORICAL DISCRIMINANT CHANNEL
# ============================================================================

print()
print("[13] HISTORICAL DISCRIMINANT CHANNEL")
print("-" * 78)

historical_translation = sp.expand(
    (Uh + h * c)**2 - Vh**2
)

historical_L = sp.expand(
    historical_translation / c**2
)

historical_Q = sp.expand(
    historical_L.subs(h, -z)
)

historical_disc = sp.factor(
    sp.discriminant(historical_Q, z)
)

print("  historical Q(z) =")
sp.pprint(historical_Q)

print()
print("  historical discriminant =")
sp.pprint(historical_disc)

cert(
    "historical discriminant = 4V²/c²",
    historical_disc,
    4 * Vh**2 / c**2,
)


# ============================================================================
# 14. PERFECT-SQUARE TEST
# ============================================================================

print()
print("[14] PERFECT-SQUARE DISCRIMINANT TEST")
print("-" * 78)

# For the canonical system:
#
#   disc(Q) = 4 J² / scale².
#
# We therefore test whether the transverse-only discriminant
# becomes a square of an invariant linear form.

G0, G1, G2 = sp.symbols("G0 G1 G2")

candidate_square = (G0 * Jsym + G1)**2 + G2

print("""
  Canonical KAPPA requirement:

      disc(Q) = rho * J²

  after an invariant affine gauge.

  A nonzero linear-in-J term can be removed by shifting J
  when appropriate; a residual constant must then also vanish
  for the exact historical Minkowski form.

  Therefore the strongest signature is:

      disc(Q) = rho * J².
""")


# ============================================================================
# 15. OPERATOR INTERPRETATION
# ============================================================================

print()
print("[15] OPERATOR INTERPRETATION")
print("-" * 78)

print("""
  The generic quadratic orbit gives:

      L(h) = h² + A*h + B.

  The spectral discriminant is:

      A² - 4B.

  The KAPPA-specific situation is stronger:

      A  = moving coordinate,
      B  = quadratic invariant,

  and the combination

      A² - 4B

  cancels every moving-coordinate contribution.

  What survives is entirely the transverse invariant.

  Thus the decisive upstream signature is:

      moving channel  ---> cancels out of discriminant
      transverse channel ---> survives as a perfect square.

  This distinguishes the historical/KAPPA geometry from a
  generic quadratic-translation orbit.
""")


# ============================================================================
# 16. RECOGNITION CONDITIONS
# ============================================================================

print()
print("[16] FINAL RECOGNITION CONDITIONS")
print("-" * 78)

print("""
  For an actual homogeneous-layer operator T, candidate K:

      1. T²(K)-2T(K)+K = constant != 0

      2. T(K)-K is affine-linear

      3. the affine-linear response defines H with

             T(H)-H = constant != 0

      4. an independent invariant J exists:

             T(J)=J

      5. the orbit polynomial Q(z)=L(-z) is monic quadratic

      6. its discriminant is independent of H

      7. ideally:

             disc(Q) = rho * J²

         for an invariant scalar rho != 0.

  Conditions 6-7 are the strongest discriminator discovered
  so far.

  A generic quadratic orbit satisfies 1-5.

  The historical/KAPPA conjugacy is distinguished by 6-7.
""")


# ============================================================================
# 17. FINAL AUDIT
# ============================================================================

print()
print("=" * 78)
print("EXPERIMENT 539 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The generic quadratic-plus-translation mechanism is now separated
from the genuinely special KAPPA property.

The next upstream search should therefore inspect the ORIGINAL
homogeneous-layer operator for a quadratic observable K satisfying:

    T²(K)-2T(K)+K = constant

and then compute its orbit polynomial.

The decisive test is:

    Is disc(Q) independent of the moving coordinate?

and, even more strongly:

    Is disc(Q) a scalar multiple of the square of an invariant
    observable?

In the historical realization:

    moving coordinate  = U
    invariant          = V
    quadratic          = U²-V²

and:

    disc(Q) = 4V²/c².

If an actual homogeneous-layer operator produces the same
discriminant-collapse pattern, that is an operator-level
identification of the historical/KAPPA structure.

Do not introduce p or q.
Do not enumerate factors.
Do not use continued fractions.
Do not perform numerical fitting.
""")
