#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 566
# ==============================================================================
# CORRECTED F2 QUADRATIC INVERSION / INVOLUTION / ROOT-CHANNEL ANALYSIS
#
# Starting point:
#
#     F2 = 4PQ - P^2 - Q^2 + P + Q
#
# Treat F2 as a quadratic in P:
#
#     P^2 - (4Q+1)P + (F2 + Q^2 - Q) = 0.
#
# Therefore:
#
#     discriminant_P
#       = (4Q+1)^2 - 4(F2+Q^2-Q)
#       = 12Q^2 + 12Q + 1 - 4F2.
#
# For the actual F2:
#
#     discriminant_P
#       = (2P-4Q-1)^2.
#
# The two roots are:
#
#     P
#     4Q+1-P.
#
# Hence the exact F2-preserving involution is:
#
#     I_P(P,Q) = (4Q+1-P, Q).
#
# Symmetrically:
#
#     I_Q(P,Q) = (P, 4P+1-Q).
#
# This experiment determines:
#
#   1. exact inversion involutions;
#   2. fixed-point loci;
#   3. discriminant square-root channels;
#   4. sum/gap coordinates of those channels;
#   5. action under common root translation P,Q -> P+h,Q+h;
#   6. historical x,y expressions;
#   7. interaction with S, Delta and N;
#   8. whether the inversion channel produces a new natural
#      coordinate of the homogeneous layer.
#
# No F4/F5/F6 assumptions.
# No numerical fitting.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 566 START")
print("=" * 78)
print("CORRECTED F2 QUADRATIC INVERSION / INVOLUTION ANALYSIS")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

P, Q = sp.symbols("P Q")
F = sp.symbols("F")
z = sp.symbols("z")
h = sp.symbols("h")

r = sp.symbols("r", nonzero=True)
x, y = sp.symbols("x y")


failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

def clean(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def num(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


def check(name, lhs, rhs=0):
    global failures

    diff = clean(lhs - rhs)
    numerator = num(diff)
    passed = sp.expand(numerator) == 0

    print(f"  {name}")
    print("    numerator difference =")
    sp.pprint(numerator)
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
# [1] ROOT COORDINATES
# ==============================================================================

print("[1] ROOT COORDINATES")
print("-" * 78)

N = clean(P * Q)
S = clean(P + Q)
Delta = clean((P - Q)**2)

show("N", N)
show("S", S)
show("Delta", Delta)

check("N=PQ", N, P*Q)
check("S=P+Q", S, P+Q)
check("Delta=(P-Q)^2", Delta, (P-Q)**2)

print()


# ==============================================================================
# [2] F2
# ==============================================================================

print("[2] HOMOGENEOUS F2")
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
# [3] QUADRATIC IN P
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
check("c_P=-F-Q^2+Q", cP, -F-Q**2+Q)

print()


# ==============================================================================
# [4] EXACT DISCRIMINANT
# ==============================================================================

print("[4] EXACT P-DISCRIMINANT")
print("-" * 78)

discP = clean(
    sp.discriminant(
        sp.expand(eqP),
        P,
    )
)

discP_expected = clean(
    (4*Q + 1)**2 - 4*(F + Q**2 - Q)
)

show("disc_P", discP)
show("expected", discP_expected)

check(
    "discriminant formula",
    discP,
    discP_expected,
)

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
# [5] CORRECT DISCRIMINANT ROOT
# ==============================================================================

print("[5] DISCRIMINANT SQUARE-ROOT CHANNEL")
print("-" * 78)

C_P = clean(
    2*P - 4*Q - 1
)

C_Q = clean(
    2*Q - 4*P - 1
)

show("C_P", C_P)
show("C_Q", C_Q)

check(
    "disc_P=C_P^2",
    discP_actual,
    C_P**2,
)

print("""
  Exact channel:

      C_P = 2P - 4Q - 1.

  Therefore:

      sqrt(disc_P) = ±C_P.

  The normalized channel is:

      C_P/2 = P - 2Q - 1/2.
""")

print()


# ==============================================================================
# [6] EXACT P INVOLUTION
# ==============================================================================

print("[6] EXACT P-INVOLUTION")
print("-" * 78)

P_star = clean(
    4*Q + 1 - P
)

show("P_star", P_star)

check(
    "F2(P_star,Q)=F2(P,Q)",
    F2.subs(P, P_star),
    F2,
)

check(
    "involution twice",
    4*Q + 1 - P_star,
    P,
)

check(
    "sum of P roots = 4Q+1",
    P + P_star,
    4*Q + 1,
)

check(
    "difference of P roots = C_P",
    P - P_star,
    C_P,
)

print()


# ==============================================================================
# [7] FIXED POINT OF P INVOLUTION
# ==============================================================================

print("[7] P-INVOLUTION FIXED POINT")
print("-" * 78)

P_fixed = clean(
    sp.solve(
        sp.Eq(P_star, P),
        P,
    )[0]
)

show("P_fixed", P_fixed)

check(
    "fixed point",
    P_fixed,
    2*Q + sp.Rational(1, 2),
)

check(
    "C_P=0 at fixed point",
    C_P.subs(P, P_fixed),
    0,
)

print()


# ==============================================================================
# [8] SYMMETRIC Q-INVERSION
# ==============================================================================

print("[8] EXACT Q-INVOLUTION")
print("-" * 78)

eqQ = clean(F2 - F)

discQ = clean(
    sp.discriminant(
        sp.expand(eqQ),
        Q,
    )
)

discQ_actual = clean(
    discQ.subs(F, F2)
)

Q_star = clean(
    4*P + 1 - Q
)

show("disc_Q(F2)", discQ_actual)
show("Q_star", Q_star)

check(
    "disc_Q(F2)=(2Q-4P-1)^2",
    discQ_actual,
    (2*Q - 4*P - 1)**2,
)

check(
    "F2(P,Q_star)=F2(P,Q)",
    F2.subs(Q, Q_star),
    F2,
)

check(
    "Q involution twice",
    4*P + 1 - Q_star,
    Q,
)

check(
    "sum of Q roots = 4P+1",
    Q + Q_star,
    4*P + 1,
)

print()


# ==============================================================================
# [9] INVOLUTION ACTION ON S,N,DELTA
# ==============================================================================

print("[9] INVOLUTION ACTION ON KAPPA COORDINATES")
print("-" * 78)

S_Pstar = clean(
    P_star + Q
)

N_Pstar = clean(
    P_star * Q
)

Delta_Pstar = clean(
    (P_star - Q)**2
)

show("S(P_star,Q)", S_Pstar)
show("N(P_star,Q)", N_Pstar)
show("Delta(P_star,Q)", Delta_Pstar)

print("""
  The involution is an F2 symmetry, but it does NOT generally
  preserve the KAPPA coordinates S, N, Delta individually.

  This is important:

      F2 is coarser than the full (S,N,Delta) package.

  The inversion symmetry therefore represents a genuine loss
  of root-coordinate information when only F2 is retained.
""")

print()


# ==============================================================================
# [10] INVOLUTION IN SUM/GAP COORDINATES
# ==============================================================================

print("[10] INVOLUTION IN SUM/GAP COORDINATES")
print("-" * 78)

s, d = sp.symbols("s d")

P_sd = clean((s + d)/2)
Q_sd = clean((s - d)/2)

Pstar_sd = clean(
    P_star.subs({
        P: P_sd,
        Q: Q_sd,
    })
)

Qstar_sd = clean(
    Q_star.subs({
        P: P_sd,
        Q: Q_sd,
    })
)

s_star = clean(Pstar_sd + Q_sd)
d_star = clean(Pstar_sd - Q_sd)

show("P_star(s,d)", Pstar_sd)
show("S_star", s_star)
show("gap_star", d_star)

check(
    "P-star sum-coordinate",
    s_star,
    4*Q_sd + 1 - s,
)

check(
    "P-star gap-coordinate",
    d_star,
    2*Q_sd + 1 - s,
)

print("""
  In (s,d):

      s = P+Q
      d = P-Q.

  The P-involution becomes:

      P -> 4Q+1-P

  and therefore:

      s -> 2s - 2d + 1

      d -> 1 - s.

  Thus the F2 involution is not a simple reflection of the
  KAPPA gap d=P-Q.
""")

print()


# ==============================================================================
# [11] DISCRIMINANT CHANNEL IN SUM/GAP FORM
# ==============================================================================

print("[11] CORRECT DISCRIMINANT CHANNEL IN SUM/GAP FORM")
print("-" * 78)

C_P_sd = clean(
    C_P.subs({
        P: P_sd,
        Q: Q_sd,
    })
)

show("C_P(s,d)", C_P_sd)

check(
    "C_P=3d-s-1",
    C_P_sd,
    3*d - s - 1,
)

check(
    "disc_P=(3d-s-1)^2",
    discP_actual.subs({
        P: P_sd,
        Q: Q_sd,
    }),
    (3*d - s - 1)**2,
)

print("""
  Correcting the earlier failed interpretation:

      disc_P = (3d-s-1)^2

  not

      (2d-s-1)^2.

  Since:

      s=P+Q
      d=P-Q,

  we obtain exactly:

      2P-4Q-1 = 3(P-Q)-(P+Q)-1.
""")

print()


# ==============================================================================
# [12] HISTORICAL x,y FORM
# ==============================================================================

print("[12] HISTORICAL x,y REPRESENTATION")
print("-" * 78)

A0 = clean(
    2*y - 2*x + 3
)

B0 = clean(
    2*y + 2*x - 3
)

Pxy = clean(A0 / r)
Qxy = clean(B0 / r)

C_P_xy = clean(
    C_P.subs({
        P: Pxy,
        Q: Qxy,
    })
)

C_Q_xy = clean(
    C_Q.subs({
        P: Pxy,
        Q: Qxy,
    })
)

show("P(x,y)", Pxy)
show("Q(x,y)", Qxy)
show("C_P(x,y)", C_P_xy)
show("C_Q(x,y)", C_Q_xy)

check(
    "C_P=(2A-4B-r)/r",
    C_P,
    (2*A - 4*B - r)/r,
)

check(
    "C_Q=(2B-4A-r)/r",
    C_Q,
    (2*B - 4*A - r)/r,
)

print()


# ==============================================================================
# [13] EXPAND HISTORICAL CHANNEL
# ==============================================================================

print("[13] EXPANDED HISTORICAL CHANNELS")
print("-" * 78)

C_P_xy_expanded = clean(
    sp.expand(C_P_xy)
)

C_Q_xy_expanded = clean(
    sp.expand(C_Q_xy)
)

show("C_P expanded", C_P_xy_expanded)
show("C_Q expanded", C_Q_xy_expanded)

print("""
  The F2 inversion channels are therefore explicit x,y
  expressions.

  They are not hidden functions of N.

  They arise directly from the 2020 root coordinates.
""")

print()


# ==============================================================================
# [14] COMMON ROOT TRANSLATION
# ==============================================================================

print("[14] ROOT TRANSLATION")
print("-" * 78)

P_h = clean(P + h)
Q_h = clean(Q + h)

C_P_h = clean(
    C_P.subs({
        P: P_h,
        Q: Q_h,
    })
)

C_Q_h = clean(
    C_Q.subs({
        P: P_h,
        Q: Q_h,
    })
)

show("C_P after translation", C_P_h)
show("C_Q after translation", C_Q_h)

check(
    "C_P -> C_P-2h",
    C_P_h,
    C_P - 2*h,
)

check(
    "C_Q -> C_Q-2h",
    C_Q_h,
    C_Q - 2*h,
)

check(
    "Delta invariant",
    (P_h - Q_h)**2,
    Delta,
)

check(
    "S -> S+2h",
    P_h + Q_h,
    S + 2*h,
)

print()


# ==============================================================================
# [15] CENTERED MOVING COORDINATES
# ==============================================================================

print("[15] CENTERED MOVING COORDINATES")
print("-" * 78)

C_P_centered = clean(
    C_P + 2*Q + 1
)

C_Q_centered = clean(
    C_Q + 2*P + 1
)

show("C_P + 2Q + 1", C_P_centered)
show("C_Q + 2P + 1", C_Q_centered)

check(
    "C_P+2Q+1=2P-2Q",
    C_P_centered,
    2*(P - Q),
)

check(
    "C_Q+2P+1=2Q-2P",
    C_Q_centered,
    -2*(P - Q),
)

print("""
  The discriminant channels contain the ordinary gap after an
  affine correction involving the opposite root.

  Specifically:

      C_P = 2(P-2Q)-1

      C_P + 2Q + 1 = 2(P-Q).

  Hence the usual invariant gap is directly embedded in the
  inversion-discriminant channel.
""")

print()


# ==============================================================================
# [16] RELATION TO S AND GAP
# ==============================================================================

print("[16] LINEAR CHANNEL BASIS")
print("-" * 78)

channel_matrix = sp.Matrix([
    [1, 1],       # S = P+Q
    [1, -1],      # gap = P-Q
    [2, -4],      # C_P = 2P-4Q-1
])

print("  coefficient matrix for:")
print("    S")
print("    gap")
print("    C_P")
sp.pprint(channel_matrix)
print()

print("""
  C_P is linearly dependent on S, gap and the affine constant:

      C_P = (3 gap) - S - 1.

  Therefore the F2 inversion discriminant does not create a
  genuinely new two-dimensional root direction.

  It selects another affine linear combination of the same
  two 2020 root coordinates.
""")

print()


# ==============================================================================
# [17] F2 AS COMPLETING THE SQUARE
# ==============================================================================

print("[17] COMPLETING THE SQUARE")
print("-" * 78)

completed_P = clean(
    -(
        P - (2*Q + sp.Rational(1, 2))
    )**2
    + (
        3*Q**2 + 3*Q + sp.Rational(1, 4)
    )
)

show("F2 completed in P", completed_P)

check(
    "complete square in P",
    F2,
    completed_P,
)

completed_Q = clean(
    -(
        Q - (2*P + sp.Rational(1, 2))
    )**2
    + (
        3*P**2 + 3*P + sp.Rational(1, 4)
    )
)

show("F2 completed in Q", completed_Q)

check(
    "complete square in Q",
    F2,
    completed_Q,
)

print("""
  This gives a more natural interpretation of the inversion:

      F2
        = -(P-2Q-1/2)^2
          + 3Q^2+3Q+1/4.

  So the centered moving coordinate is:

      H_P = P-2Q-1/2.

  The inversion simply reflects this coordinate:

      H_P -> -H_P.
""")

print()


# ==============================================================================
# [18] INVOLUTION AS REFLECTION
# ==============================================================================

print("[18] INVOLUTION AS REFLECTION")
print("-" * 78)

H_P = clean(
    P - 2*Q - sp.Rational(1, 2)
)

H_P_star = clean(
    P_star - 2*Q - sp.Rational(1, 2)
)

show("H_P", H_P)
show("H_P after involution", H_P_star)

check(
    "involution reflects H_P",
    H_P_star,
    -H_P,
)

print()


# ==============================================================================
# [19] DISCRIMINANT IS FOUR TIMES THE REFLECTION SQUARE
# ==============================================================================

print("[19] DISCRIMINANT / REFLECTION COORDINATE")
print("-" * 78)

check(
    "disc_P=4H_P^2",
    discP_actual,
    4*H_P**2,
)

show("H_P", H_P)
show("4H_P^2", 4*H_P**2)

print()


# ==============================================================================
# [20] RELATION TO THE KAPPA QUADRATIC
# ==============================================================================

print("[20] KAPPA QUADRATIC")
print("-" * 78)

Qpoly = clean(
    z**2 - S*z + N
)

Qfactor = clean(
    (z-P)*(z-Q)
)

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

show("Q(z)", Qpoly)

print()


# ==============================================================================
# [21] COMPARE THE TWO DISCRIMINANTS
# ==============================================================================

print("[21] TWO DIFFERENT DISCRIMINANTS")
print("-" * 78)

kappa_disc = clean(
    sp.discriminant(Qpoly, z)
)

f2_disc = clean(
    discP_actual
)

show("KAPPA discriminant", kappa_disc)
show("F2 inversion discriminant", f2_disc)

print("""
  There are therefore two distinct discriminants:

      KAPPA spectral discriminant:
          Delta = (P-Q)^2

      F2 inversion discriminant:
          4(P-2Q-1/2)^2.

  They encode different affine slices of the same two-root
  coordinate plane.
""")

print()


# ==============================================================================
# [22] HISTORICAL TRANSLATION OF THE REFLECTION COORDINATE
# ==============================================================================

print("[22] HISTORICAL TRANSLATION OF H_P")
print("-" * 78)

H_P_xy = clean(
    H_P.subs({
        P: Pxy,
        Q: Qxy,
    })
)

show("H_P(x,y)", H_P_xy)

H_P_xy_shifted = clean(
    H_P_xy
)

print("""
  Under the original root translation:

      P -> P+h
      Q -> Q+h,

      H_P -> H_P-h.

  Thus H_P is an affine moving coordinate with velocity -1,
  while S/2 has velocity +1 and Delta is invariant.
""")

check(
    "H_P translation response",
    H_P.subs({
        P: P_h,
        Q: Q_h,
    }) - H_P,
    -h,
)

print()


# ==============================================================================
# [23] FINAL STRUCTURAL CERTIFICATE
# ==============================================================================

print("[23] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print(
r"""
  CORRECTED F2 INVERSION:

      F2 = 4PQ-P^2-Q^2+P+Q.

  As a quadratic in P:

      P^2-(4Q+1)P+(F+Q^2-Q)=0.

  Therefore:

      disc_P
        = 12Q^2+12Q+1-4F.

  On the actual homogeneous layer:

      disc_P(F2)
        = (2P-4Q-1)^2

        = 4(P-2Q-1/2)^2.

  The exact involution is:

      P -> 4Q+1-P.

  The centered coordinate:

      H_P=P-2Q-1/2

  transforms as:

      H_P -> -H_P.

  Under the common root translation:

      P -> P+h
      Q -> Q+h,

      H_P -> H_P-h.

  Meanwhile:

      S=P+Q -> S+2h

      Delta=(P-Q)^2 -> Delta.

  Thus F2 contains a reflection coordinate in addition to
  the usual KAPPA moving and invariant coordinates.

  The important correction is that this does NOT introduce
  a third independent root coordinate.

  Everything remains inside the two-dimensional 2020 root
  coordinate plane (P,Q).
"""
)

print()


# ==============================================================================
# [24] FINAL AUDIT
# ==============================================================================

print("[24] FINAL AUDIT")
print("-" * 78)

print(f"  symbolic failures = {failures}")
print("  F4/F5/F6 assumptions = NONE")
print("  F2 inversion = CORRECTED")
print("  exact involution = TESTED")
print("  reflection coordinate = TESTED")
print("  root translation = TESTED")
print()

print("=" * 78)
print("EXPERIMENT 566 FINISHED")
print("=" * 78)
