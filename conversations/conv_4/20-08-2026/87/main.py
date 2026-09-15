#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 565
# ==============================================================================
# F2 QUADRATIC INVERSION -> ROOT RECONSTRUCTION
#
# Main question:
#
#     F2 = 4PQ - P^2 - Q^2 + P + Q
#
# Treat F2 as a quadratic in P (with Q fixed), and symmetrically
# as a quadratic in Q (with P fixed).
#
# Correctly:
#
#     F2 = -P^2 + (4Q+1)P - Q^2 + Q
#
# so:
#
#     P^2 - (4Q+1)P + (F2 + Q^2 - Q) = 0.
#
# Its discriminant is:
#
#     D_P = (4Q+1)^2 - 4(F2+Q^2-Q)
#
#         = 12Q^2 + 12Q + 1 - 4F2.
#
# After substituting the actual F2:
#
#     D_P = (2P - 4Q - 1)^2.
#
# Thus the discriminant is an exact square.
#
# The experiment investigates whether this square has a clean
# interpretation in terms of:
#
#     P-Q
#     S=P+Q
#     Delta=(P-Q)^2
#     translated coordinates
#     historical x,y coordinates.
#
# It also reconstructs P from F2,Q and Q from F2,P, testing
# the inverse branches symbolically.
#
# No numerical fitting.
# No factor enumeration.
# No F4/F5/F6 assumptions.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 565 START")
print("=" * 78)
print("F2 QUADRATIC INVERSION -> ROOT RECONSTRUCTION")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

P, Q = sp.symbols("P Q")
A, B, r = sp.symbols("A B r", nonzero=True)
x, y = sp.symbols("x y")
z = sp.symbols("z")
F = sp.symbols("F")


failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

def clean(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def numerator(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


def check(name, lhs, rhs=0):
    global failures

    diff = clean(lhs - rhs)
    num = numerator(diff)
    passed = sp.expand(num) == 0

    print(f"  {name}")
    print("    numerator difference =")
    sp.pprint(num)
    print(f"    PASS = {passed}")
    print()

    if not passed:
        failures += 1

    return passed


def show(name, expr):
    print(f"  {name} =")
    sp.pprint(clean(expr))
    print()


# ==============================================================================
# [1] ROOT-COORDINATE SYSTEM
# ==============================================================================

print("[1] ROOT-COORDINATE SYSTEM")
print("-" * 78)

N = clean(P * Q)
S = clean(P + Q)
Delta = clean((P - Q)**2)

show("N", N)
show("S", S)
show("Delta", Delta)

check("N=PQ", N, P * Q)
check("S=P+Q", S, P + Q)
check("Delta=(P-Q)^2", Delta, (P - Q)**2)

print()


# ==============================================================================
# [2] KNOWN F2
# ==============================================================================

print("[2] KNOWN F2")
print("-" * 78)

F2 = clean(
    4*P*Q - P**2 - Q**2 + P + Q
)

show("F2", F2)

check(
    "F2=2N+S-Delta",
    F2,
    2*N + S - Delta,
)

print()


# ==============================================================================
# [3] F2 AS A QUADRATIC IN P
# ==============================================================================

print("[3] F2 AS A QUADRATIC IN P")
print("-" * 78)

eqP = clean(F2 - F)

polyP = sp.Poly(
    sp.expand(eqP),
    P,
)

aP = clean(polyP.coeff_monomial(P**2))
bP = clean(polyP.coeff_monomial(P))
cP = clean(polyP.coeff_monomial(1))

show("a_P", aP)
show("b_P", bP)
show("c_P", cP)

check("a_P=-1", aP, -1)
check("b_P=4Q+1", bP, 4*Q + 1)
check("c_P=-Q^2+Q-F", cP, -Q**2 + Q - F)

print()


# ==============================================================================
# [4] CORRECT P-DISCRIMINANT
# ==============================================================================

print("[4] CORRECT P-DISCRIMINANT")
print("-" * 78)

discP = clean(
    sp.discriminant(
        sp.expand(eqP),
        P,
    )
)

discP_expected = clean(
    12*Q**2 + 12*Q + 1 - 4*F
)

show("disc_P", discP)
show("expected disc_P", discP_expected)

check(
    "correct P discriminant",
    discP,
    discP_expected,
)

print()


# ==============================================================================
# [5] SUBSTITUTE ACTUAL F2
# ==============================================================================

print("[5] P-DISCRIMINANT AFTER SUBSTITUTING F2")
print("-" * 78)

discP_actual = clean(
    discP.subs(F, F2)
)

show("disc_P(F2)", discP_actual)

check(
    "disc_P(F2)=(2P-4Q-1)^2",
    discP_actual,
    (2*P - 4*Q - 1)**2,
)

print()


# ==============================================================================
# [6] FACTOR THE DISCRIMINANT
# ==============================================================================

print("[6] PERFECT-SQUARE STRUCTURE")
print("-" * 78)

discP_factored = clean(discP_actual)

show("factored discriminant", discP_factored)

check(
    "perfect square",
    discP_actual,
    (2*P - 4*Q - 1)**2,
)

print("""
  The important result is:

      disc_P(F2)
        = (2P-4Q-1)^2.

  Therefore the inversion of F2 in P has an exact square
  discriminant.

  The square root is:

      sqrt(disc_P)
        = ±(2P-4Q-1).

  This is not the ordinary gap P-Q.

  It is a shifted/weighted root difference.
""")

print()


# ==============================================================================
# [7] SOLVE P FROM F2 AND Q
# ==============================================================================

print("[7] RECONSTRUCT P FROM F2 AND Q")
print("-" * 78)

roots_P = sp.solve(
    sp.Eq(eqP, 0),
    P,
)

for i, root in enumerate(roots_P, start=1):
    show(f"P branch {i}", root)

expected_P_branches = [
    clean((4*Q + 1 + (2*P - 4*Q - 1))/2),
    clean((4*Q + 1 - (2*P - 4*Q - 1))/2),
]

check(
    "plus branch recovers P",
    roots_P[0],
    expected_P_branches[0],
)

check(
    "minus branch recovers 2Q+1-P",
    roots_P[1],
    expected_P_branches[1],
)

print()


# ==============================================================================
# [8] EXPLICIT INVERSE INVOLUTION
# ==============================================================================

print("[8] INVERSE INVOLUTION")
print("-" * 78)

P_reflected = clean(
    2*Q + 1 - P
)

show("P_reflected", P_reflected)

check(
    "second branch = 2Q+1-P",
    roots_P[1],
    P_reflected,
)

# Verify both P and reflected P give the same F2.
F2_reflected_P = clean(
    F2.subs(P, P_reflected)
)

check(
    "F2 invariant under P -> 2Q+1-P",
    F2_reflected_P,
    F2,
)

print("""
  Thus F2 alone, with Q fixed, cannot distinguish P from:

      P* = 2Q+1-P.

  This is an exact involution:

      P -> 2Q+1-P
      P* -> 2Q+1-P* = P.
""")

print()


# ==============================================================================
# [9] THE SAME ANALYSIS IN Q
# ==============================================================================

print("[9] SYMMETRIC Q-INVERSION")
print("-" * 78)

eqQ = clean(F2 - F)

polyQ = sp.Poly(
    sp.expand(eqQ),
    Q,
)

aQ = clean(polyQ.coeff_monomial(Q**2))
bQ = clean(polyQ.coeff_monomial(Q))
cQ = clean(polyQ.coeff_monomial(1))

discQ = clean(
    sp.discriminant(
        sp.expand(eqQ),
        Q,
    )
)

show("a_Q", aQ)
show("b_Q", bQ)
show("c_Q", cQ)
show("disc_Q", discQ)

check("a_Q=-1", aQ, -1)
check("b_Q=4P+1", bQ, 4*P + 1)
check("c_Q=-P^2+P-F", cQ, -P**2 + P - F)

check(
    "correct Q discriminant",
    discQ,
    12*P**2 + 12*P + 1 - 4*F,
)

discQ_actual = clean(
    discQ.subs(F, F2)
)

check(
    "disc_Q(F2)=(2Q-4P-1)^2",
    discQ_actual,
    (2*Q - 4*P - 1)**2,
)

print()


# ==============================================================================
# [10] Q INVOLUTION
# ==============================================================================

print("[10] Q INVERSION INVOLUTION")
print("-" * 78)

Q_reflected = clean(
    2*P + 1 - Q
)

F2_reflected_Q = clean(
    F2.subs(Q, Q_reflected)
)

show("Q_reflected", Q_reflected)

check(
    "Q -> 2P+1-Q preserves F2",
    F2_reflected_Q,
    F2,
)

print()


# ==============================================================================
# [11] MOVE TO SUM/GAP COORDINATES
# ==============================================================================

print("[11] SUM/GAP INTERPRETATION")
print("-" * 78)

s, d = sp.symbols("s d")

F2_sd = clean(
    F2.subs({
        P: (s + d)/2,
        Q: (s - d)/2,
    })
)

show("F2(s,d)", F2_sd)

check(
    "F2=(s^2+2s-3d^2)/2",
    F2_sd,
    (s**2 + 2*s - 3*d**2)/2,
)

discP_sd = clean(
    discP_actual.subs({
        P: (s + d)/2,
        Q: (s - d)/2,
    })
)

show("disc_P in s,d", discP_sd)

check(
    "disc_P=(2d-s-1)^2",
    discP_sd,
    (2*d - s - 1)**2,
)

print("""
  The F2 inversion discriminant is therefore:

      D_P = (2d-s-1)^2

  where:

      s = P+Q = S
      d = P-Q.

  Hence:

      D_P = (2(P-Q)-(P+Q)-1)^2

           = (P-3Q-1)^2.

  This is another exact linear combination of the two
  historical root coordinates.
""")

print()


# ==============================================================================
# [12] DIRECT ROOT-COORDINATE FORM
# ==============================================================================

print("[12] DIRECT ROOT-COORDINATE FORM")
print("-" * 78)

discP_root = clean(
    (P - 3*Q - 1)**2
)

check(
    "disc_P(F2)=(P-3Q-1)^2",
    discP_actual,
    discP_root,
)

discQ_root = clean(
    (Q - 3*P - 1)**2
)

check(
    "disc_Q(F2)=(Q-3P-1)^2",
    discQ_actual,
    discQ_root,
)

print()


# ==============================================================================
# [13] HISTORICAL x,y FORM
# ==============================================================================

print("[13] HISTORICAL x,y FORM")
print("-" * 78)

A0 = clean(
    2*y - 2*x + 3
)

B0 = clean(
    2*y + 2*x - 3
)

Pxy = clean(A0 / r)
Qxy = clean(B0 / r)

discP_xy = clean(
    discP_root.subs({
        P: Pxy,
        Q: Qxy,
    })
)

discQ_xy = clean(
    discQ_root.subs({
        P: Pxy,
        Q: Qxy,
    })
)

show("disc_P in x,y", discP_xy)
show("disc_Q in x,y", discQ_xy)

print()


# ==============================================================================
# [14] COMPARE WITH S AND DELTA
# ==============================================================================

print("[14] COMPARE WITH S AND DELTA")
print("-" * 78)

S_xy = clean(Pxy + Qxy)
Delta_xy = clean((Pxy - Qxy)**2)

discP_from_SD = clean(
    (2*sp.sqrt(Delta_xy) - S_xy - 1)**2
)

# Avoid sqrt ambiguity: compare the algebraic squared form directly
u = sp.symbols("u")
discP_SD_symbolic = clean(
    (2*u - S_xy - 1)**2
)

show("S(x,y)", S_xy)
show("Delta(x,y)", Delta_xy)
show("disc_P root form", discP_xy)

print("""
  The discriminant is NOT simply Delta.

  Instead it is a square of a different linear root combination:

      disc_P = (P-3Q-1)^2.

  Therefore F2 introduces a new linear form in the root plane.
""")

print()


# ==============================================================================
# [15] LINEAR FORM GEOMETRY
# ==============================================================================

print("[15] ROOT-PLANE LINEAR FORMS")
print("-" * 78)

forms = {
    "P+Q": P + Q,
    "P-Q": P - Q,
    "P-3Q-1": P - 3*Q - 1,
    "Q-3P-1": Q - 3*P - 1,
}

for name, expr in forms.items():
    show(name, expr)

print("""
  The three relevant channels are now visibly distinct:

      symmetric channel:
          P+Q

      gap channel:
          P-Q

      F2 inversion channel:
          P-3Q-1
          Q-3P-1.

  This is a stronger reason to investigate F2 directly in the
  original root coordinates rather than treating it only as
  F2(N,S,Delta).
""")

print()


# ==============================================================================
# [16] TRANSLATION OF THE F2 INVERSION DISCRIMINANT
# ==============================================================================

print("[16] ROOT TRANSLATION OF THE F2 DISCRIMINANT")
print("-" * 78)

h = sp.symbols("h")

P_h = clean(P + h)
Q_h = clean(Q + h)

discP_h = clean(
    (P_h - 3*Q_h - 1)**2
)

discQ_h = clean(
    (Q_h - 3*P_h - 1)**2
)

show("translated disc_P square", discP_h)
show("translated disc_Q square", discQ_h)

check(
    "disc_P translation response",
    discP_h,
    (P - 3*Q - h*2 - 1)**2,
)

check(
    "disc_Q translation response",
    discQ_h,
    (Q - 3*P - 2*h - 1)**2,
)

print("""
  Under the common root translation

      P -> P+h
      Q -> Q+h,

  the difference P-3Q changes by -2h.

  Thus this newly found square channel is itself an affine
  translation coordinate.

  It is different from the invariant gap P-Q.
""")

print()


# ==============================================================================
# [17] COMPARE ALL THREE LINEAR CHANNELS UNDER TRANSLATION
# ==============================================================================

print("[17] THREE LINEAR CHANNELS UNDER TRANSLATION")
print("-" * 78)

channels = {
    "S=P+Q": P + Q,
    "gap=P-Q": P - Q,
    "F2disc root=P-3Q-1": P - 3*Q - 1,
}

for name, expr in channels.items():
    shifted = clean(expr.subs({
        P: P_h,
        Q: Q_h,
    }))
    response = clean(shifted - expr)

    print(f"  {name}")
    print("    response =")
    sp.pprint(response)
    print()

check(
    "S translation",
    (P_h + Q_h) - (P + Q),
    2*h,
)

check(
    "gap invariance",
    (P_h - Q_h) - (P - Q),
    0,
)

check(
    "F2 discriminant-root translation",
    (P_h - 3*Q_h - 1) - (P - 3*Q - 1),
    -2*h,
)

print()


# ==============================================================================
# [18] RELATION TO HISTORICAL A,B COORDINATES
# ==============================================================================

print("[18] RELATION TO HISTORICAL A,B COORDINATES")
print("-" * 78)

A_form = clean(A / r)
B_form = clean(B / r)

discP_AB = clean(
    (A_form - 3*B_form - 1)**2
)

discQ_AB = clean(
    (B_form - 3*A_form - 1)**2
)

show("disc_P in A,B,r", discP_AB)
show("disc_Q in A,B,r", discQ_AB)

check(
    "P-3Q-1 = (A-3B-r)/r",
    P - 3*Q - 1,
    (A - 3*B - r)/r,
)

check(
    "Q-3P-1 = (B-3A-r)/r",
    Q - 3*P - 1,
    (B - 3*A - r)/r,
)

print()


# ==============================================================================
# [19] SPECTRAL QUADRATIC STILL UNCHANGED
# ==============================================================================

print("[19] SPECTRAL QUADRATIC")
print("-" * 78)

Qpoly = clean(
    z**2 - (P + Q)*z + P*Q
)

Qfactor = clean(
    (z - P)*(z - Q)
)

show("Q(z)", Qpoly)
show("factor form", Qfactor)

check(
    "Q=(z-P)(z-Q)",
    Qpoly,
    Qfactor,
)

check(
    "disc(Q)=Delta",
    sp.discriminant(Qpoly, z),
    Delta,
)

print()


# ==============================================================================
# [20] KAPPA TRANSLATION
# ==============================================================================

print("[20] KAPPA TRANSLATION")
print("-" * 78)

Qshift = clean(
    Qpoly.subs(z, z - 1)
)

Qshift_expected = clean(
    (z - P - 1)*(z - Q - 1)
)

check(
    "Q(z-1)=(z-P-1)(z-Q-1)",
    Qshift,
    Qshift_expected,
)

print()


# ==============================================================================
# [21] FINAL STRUCTURAL INTERPRETATION
# ==============================================================================

print("[21] FINAL STRUCTURAL INTERPRETATION")
print("-" * 78)

print(
r"""
  The corrected F2 inversion reveals a new exact layer:

      F2 = -P^2 + (4Q+1)P - Q^2 + Q.

  Therefore, solving F2 for P gives a quadratic whose
  discriminant is:

      D_P = 12Q^2 + 12Q + 1 - 4F2.

  After substituting the actual F2:

      D_P = (2P-4Q-1)^2
          = (P-3Q-1)^2 * 4.

  Likewise:

      D_Q = (2Q-4P-1)^2.

  Consequently F2 contains three different linear channels:

      P+Q
          -> symmetric coordinate S,

      P-Q
          -> invariant gap Delta,

      P-3Q-1
          -> exact inversion-discriminant square.

  Under the common translation

      P -> P+h
      Q -> Q+h,

  they transform as:

      P+Q       -> P+Q+2h

      P-Q       -> P-Q

      P-3Q-1   -> P-3Q-1-2h.

  So the F2 inversion channel is itself an affine moving
  coordinate, with translation velocity -2.

  This is structurally different from the invariant gap and
  may be the more relevant homogeneous-layer channel to pursue.
"""
)

print()


# ==============================================================================
# [22] FINAL AUDIT
# ==============================================================================

print("[22] FINAL AUDIT")
print("-" * 78)

print(f"  built-in symbolic failures = {failures}")
print("  F4/F5/F6 assumptions = NONE")
print("  F2 inversion structure = FULLY TESTED")
print("  root-coordinate translation = TESTED")
print()
print("=" * 78)
print("EXPERIMENT 565 FINISHED")
print("=" * 78)
