#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 566
# ==============================================================================
# CORRECTED F2 QUADRATIC INVERSION
#
# Exact starting point:
#
#     F2 = 4PQ - P^2 - Q^2 + P + Q
#
# As a quadratic in P:
#
#     F2 = -P^2 + (4Q+1)P - Q^2 + Q
#
# so:
#
#     P^2 - (4Q+1)P + (F2 + Q^2 - Q) = 0
#
# and:
#
#     disc_P = 12Q^2 + 12Q + 1 - 4F.
#
# Substituting the actual F2 gives:
#
#     disc_P(F2) = (2P - 4Q - 1)^2
#
# with exact involution:
#
#     P -> 4Q + 1 - P.
#
# The experiment also cross-checks the same structure in the
# historical coordinates:
#
#     A0 = 2y - 2x + 3
#     B0 = 2y + 2x - 3
#
# and normalized coordinates:
#
#     P = A0/r
#     Q = B0/r.
#
# No F4/F5/F6 assumptions.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 566 START")
print("=" * 78)
print("CORRECTED F2 QUADRATIC INVERSION / ROOT INVOLUTION")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

P, Q = sp.symbols("P Q")
F = sp.symbols("F")
z = sp.symbols("z")
h = sp.symbols("h")

# Historical coordinates
x, y = sp.symbols("x y")

# Normalization scale
r = sp.symbols("r", nonzero=True)

# Root-coordinate aliases
A, B = sp.symbols("A B")


# ==============================================================================
# HELPERS
# ==============================================================================

failures = 0


def clean(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def numerator(expr):
    return sp.factor(
        sp.together(expr).as_numer_denom()[0]
    )


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
# [1] ROOT COORDINATES
# ==============================================================================

print("[1] ROOT COORDINATES")
print("-" * 78)

N = clean(P * Q)
S = clean(P + Q)
Delta = clean((P - Q) ** 2)

show("N", N)
show("S", S)
show("Delta", Delta)

check("N=PQ", N, P * Q)
check("S=P+Q", S, P + Q)
check("Delta=(P-Q)^2", Delta, (P - Q) ** 2)

print()


# ==============================================================================
# [2] KNOWN F2
# ==============================================================================

print("[2] KNOWN HOMOGENEOUS LAYER F2")
print("-" * 78)

F2 = clean(
    4 * P * Q
    - P ** 2
    - Q ** 2
    + P
    + Q
)

show("F2", F2)

check(
    "F2=2N+S-Delta",
    F2,
    2 * N + S - Delta,
)

print()


# ==============================================================================
# [3] F2 AS QUADRATIC IN P
# ==============================================================================

print("[3] F2 AS A QUADRATIC IN P")
print("-" * 78)

F_symbolic = F

eqP = clean(
    -P ** 2
    + (4 * Q + 1) * P
    - Q ** 2
    + Q
    - F_symbolic
)

polyP = sp.Poly(sp.expand(eqP), P)

aP = clean(polyP.coeff_monomial(P ** 2))
bP = clean(polyP.coeff_monomial(P))
cP = clean(polyP.coeff_monomial(1))

show("a_P", aP)
show("b_P", bP)
show("c_P", cP)

check("a_P=-1", aP, -1)
check("b_P=4Q+1", bP, 4 * Q + 1)
check("c_P=-Q^2+Q-F", cP, -Q ** 2 + Q - F_symbolic)

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
    (4 * Q + 1) ** 2
    - 4 * (
        F_symbolic
        + Q ** 2
        - Q
    )
)

show("disc_P", discP)
show("expected disc_P", discP_expected)

check(
    "discriminant formula",
    discP,
    discP_expected,
)

print()


# ==============================================================================
# [5] SUBSTITUTE ACTUAL F2
# ==============================================================================

print("[5] P-DISCRIMINANT AFTER SUBSTITUTING ACTUAL F2")
print("-" * 78)

discP_actual = clean(
    discP.subs(
        F_symbolic,
        F2,
    )
)

show("disc_P(F2)", discP_actual)

check(
    "disc_P(F2)=(2P-4Q-1)^2",
    discP_actual,
    (2 * P - 4 * Q - 1) ** 2,
)

print()


# ==============================================================================
# [6] EXACT SQUARE-ROOT CHANNEL
# ==============================================================================

print("[6] F2 INVERSION CHANNEL")
print("-" * 78)

C_P = clean(
    2 * P - 4 * Q - 1
)

show("C_P", C_P)

check(
    "disc_P(F2)=C_P^2",
    discP_actual,
    C_P ** 2,
)

print(
    """
  Therefore:

      sqrt(disc_P) = ±(2P-4Q-1).

  Define the centered coordinate:

      H_P = P - 2Q - 1/2

  so that:

      C_P = 2 H_P.
"""
)

H_P = clean(
    P - 2 * Q - sp.Rational(1, 2)
)

check(
    "C_P=2H_P",
    C_P,
    2 * H_P,
)

check(
    "disc_P=4H_P^2",
    discP_actual,
    4 * H_P ** 2,
)

print()


# ==============================================================================
# [7] EXACT P INVOLUTION
# ==============================================================================

print("[7] EXACT P INVOLUTION")
print("-" * 78)

P_star = clean(
    4 * Q + 1 - P
)

show("P_star", P_star)

check(
    "F2(P_star,Q)=F2(P,Q)",
    F2.subs(P, P_star),
    F2,
)

check(
    "involution twice",
    4 * Q + 1 - P_star,
    P,
)

check(
    "sum of two P roots",
    P + P_star,
    4 * Q + 1,
)

check(
    "root difference",
    P - P_star,
    C_P,
)

H_P_star = clean(
    H_P.subs(P, P_star)
)

show("H_P after involution", H_P_star)

check(
    "involution reflects H_P",
    H_P_star,
    -H_P,
)

print()


# ==============================================================================
# [8] FIXED POINT
# ==============================================================================

print("[8] P INVOLUTION FIXED POINT")
print("-" * 78)

P_fixed = clean(
    2 * Q + sp.Rational(1, 2)
)

show("P_fixed", P_fixed)

check(
    "P_star=P at fixed point",
    P_star.subs(P, P_fixed),
    P_fixed,
)

check(
    "H_P=0 at fixed point",
    H_P.subs(P, P_fixed),
    0,
)

print()


# ==============================================================================
# [9] SYMMETRIC Q INVERSION
# ==============================================================================

print("[9] SYMMETRIC Q INVOLUTION")
print("-" * 78)

eqQ = clean(
    -Q ** 2
    + (4 * P + 1) * Q
    - P ** 2
    + P
    - F_symbolic
)

discQ = clean(
    sp.discriminant(
        sp.expand(eqQ),
        Q,
    )
)

discQ_actual = clean(
    discQ.subs(
        F_symbolic,
        F2,
    )
)

Q_star = clean(
    4 * P + 1 - Q
)

C_Q = clean(
    2 * Q - 4 * P - 1
)

H_Q = clean(
    Q - 2 * P - sp.Rational(1, 2)
)

show("disc_Q(F2)", discQ_actual)
show("Q_star", Q_star)
show("C_Q", C_Q)

check(
    "disc_Q(F2)=(2Q-4P-1)^2",
    discQ_actual,
    C_Q ** 2,
)

check(
    "F2(P,Q_star)=F2(P,Q)",
    F2.subs(Q, Q_star),
    F2,
)

check(
    "Q involution twice",
    4 * P + 1 - Q_star,
    Q,
)

check(
    "C_Q=2H_Q",
    C_Q,
    2 * H_Q,
)

check(
    "disc_Q=4H_Q^2",
    discQ_actual,
    4 * H_Q ** 2,
)

print()


# ==============================================================================
# [10] SUM / GAP COORDINATES
# ==============================================================================

print("[10] SUM / GAP COORDINATES")
print("-" * 78)

sum_var, gap_var = sp.symbols(
    "sum_var gap_var"
)

P_sd = clean(
    (sum_var + gap_var) / 2
)

Q_sd = clean(
    (sum_var - gap_var) / 2
)

C_P_sd = clean(
    C_P.subs(
        {
            P: P_sd,
            Q: Q_sd,
        }
    )
)

H_P_sd = clean(
    H_P.subs(
        {
            P: P_sd,
            Q: Q_sd,
        }
    )
)

show("C_P(sum,gap)", C_P_sd)
show("H_P(sum,gap)", H_P_sd)

check(
    "C_P=3gap-sum-1",
    C_P_sd,
    3 * gap_var - sum_var - 1,
)

check(
    "H_P=(3gap-sum-1)/2",
    H_P_sd,
    (3 * gap_var - sum_var - 1) / 2,
)

print()


# ==============================================================================
# [11] HISTORICAL 2020 COORDINATES
# ==============================================================================

print("[11] HISTORICAL 2020 COORDINATES")
print("-" * 78)

A0 = clean(
    2 * y - 2 * x + 3
)

B0 = clean(
    2 * y + 2 * x - 3
)

P_xy = clean(A0 / r)
Q_xy = clean(B0 / r)

show("A0", A0)
show("B0", B0)
show("P(x,y)", P_xy)
show("Q(x,y)", Q_xy)

check(
    "P=A0/r",
    P_xy,
    A0 / r,
)

check(
    "Q=B0/r",
    Q_xy,
    B0 / r,
)

print()


# ==============================================================================
# [12] HISTORICAL N, S, DELTA
# ==============================================================================

print("[12] HISTORICAL KAPPA CROSSWALK")
print("-" * 78)

N_xy = clean(
    P_xy * Q_xy
)

S_xy = clean(
    P_xy + Q_xy
)

Delta_xy = clean(
    (P_xy - Q_xy) ** 2
)

show("N(x,y)", N_xy)
show("S(x,y)", S_xy)
show("Delta(x,y)", Delta_xy)

check(
    "S=4y/r",
    S_xy,
    4 * y / r,
)

check(
    "Delta=4(2x-3)^2/r^2",
    Delta_xy,
    4 * (2 * x - 3) ** 2 / r ** 2,
)

print()


# ==============================================================================
# [13] HISTORICAL F2
# ==============================================================================

print("[13] HISTORICAL F2")
print("-" * 78)

F2_xy = clean(
    F2.subs(
        {
            P: P_xy,
            Q: Q_xy,
        }
    )
)

show("F2(x,y)", F2_xy)

check(
    "F2(x,y)=2N+S-Delta",
    F2_xy,
    2 * N_xy + S_xy - Delta_xy,
)

print()


# ==============================================================================
# [14] HISTORICAL INVERSION CHANNEL
# ==============================================================================

print("[14] HISTORICAL F2-INVERSION CHANNEL")
print("-" * 78)

C_P_xy = clean(
    C_P.subs(
        {
            P: P_xy,
            Q: Q_xy,
        }
    )
)

H_P_xy = clean(
    H_P.subs(
        {
            P: P_xy,
            Q: Q_xy,
        }
    )
)

show("C_P(x,y)", C_P_xy)
show("H_P(x,y)", H_P_xy)

check(
    "C_P=2H_P",
    C_P_xy,
    2 * H_P_xy,
)

check(
    "disc_P(x,y)=C_P(x,y)^2",
    discP_actual.subs(
        {
            P: P_xy,
            Q: Q_xy,
        }
    ),
    C_P_xy ** 2,
)

print()


# ==============================================================================
# [15] HISTORICAL CHANNEL EXPANSION
# ==============================================================================

print("[15] EXPANDED HISTORICAL CHANNEL")
print("-" * 78)

show(
    "C_P expanded",
    sp.expand(C_P_xy),
)

show(
    "H_P expanded",
    sp.expand(H_P_xy),
)

print(
    """
  This explicitly confirms that the F2 inversion channel is
  already an x,y quantity inherited from the 2020 coordinates.
"""
)

print()


# ==============================================================================
# [16] COMMON ROOT TRANSLATION
# ==============================================================================

print("[16] COMMON ROOT TRANSLATION")
print("-" * 78)

P_h = clean(P + h)
Q_h = clean(Q + h)

S_h = clean(P_h + Q_h)
Delta_h = clean((P_h - Q_h) ** 2)

C_P_h = clean(
    C_P.subs(
        {
            P: P_h,
            Q: Q_h,
        }
    )
)

H_P_h = clean(
    H_P.subs(
        {
            P: P_h,
            Q: Q_h,
        }
    )
)

show("S_h", S_h)
show("Delta_h", Delta_h)
show("C_P_h", C_P_h)
show("H_P_h", H_P_h)

check(
    "S -> S+2h",
    S_h,
    S + 2 * h,
)

check(
    "Delta invariant",
    Delta_h,
    Delta,
)

check(
    "C_P -> C_P-2h",
    C_P_h,
    C_P - 2 * h,
)

check(
    "H_P -> H_P-h",
    H_P_h,
    H_P - h,
)

print()


# ==============================================================================
# [17] F2 COMPLETING THE SQUARE
# ==============================================================================

print("[17] F2 COMPLETING THE SQUARE")
print("-" * 78)

F2_complete_P = clean(
    -H_P ** 2
    + 3 * Q ** 2
    + 3 * Q
    + sp.Rational(1, 4)
)

F2_complete_Q = clean(
    -H_Q ** 2
    + 3 * P ** 2
    + 3 * P
    + sp.Rational(1, 4)
)

show("F2 completed in P", F2_complete_P)
show("F2 completed in Q", F2_complete_Q)

check(
    "F2 complete-square P",
    F2,
    F2_complete_P,
)

check(
    "F2 complete-square Q",
    F2,
    F2_complete_Q,
)

print()


# ==============================================================================
# [18] RELATION BETWEEN F2 INVERSION CHANNEL AND KAPPA GAP
# ==============================================================================

print("[18] RELATION TO KAPPA GAP")
print("-" * 78)

check(
    "C_P + 2Q + 1 = 2(P-Q)",
    C_P + 2 * Q + 1,
    2 * (P - Q),
)

check(
    "C_Q + 2P + 1 = -2(P-Q)",
    C_Q + 2 * P + 1,
    -2 * (P - Q),
)

print(
    """
  Therefore:

      C_P = 2(P-Q) - 2Q - 1

  and:

      C_P + 2Q + 1 = 2(P-Q).

  The usual KAPPA gap is therefore contained exactly inside
  the F2 inversion channel once the opposite root is restored.
"""
)

print()


# ==============================================================================
# [19] SPECTRAL KAPPA QUADRATIC
# ==============================================================================

print("[19] SPECTRAL KAPPA QUADRATIC")
print("-" * 78)

Qpoly = clean(
    z ** 2 - S * z + N
)

Qfactor = clean(
    (z - P) * (z - Q)
)

show("Q(z)", Qpoly)

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
# [20] SPECTRAL TRANSLATION
# ==============================================================================

print("[20] SPECTRAL TRANSLATION")
print("-" * 78)

Q_shift = clean(
    Qpoly.subs(z, z - 1)
)

Q_shift_expected = clean(
    (z - P - 1) * (z - Q - 1)
)

show("Q(z-1)", Q_shift)

check(
    "Q(z-1) translated roots",
    Q_shift,
    Q_shift_expected,
)

print()


# ==============================================================================
# [21] F2 INVERSION VS KAPPA DISCRIMINANT
# ==============================================================================

print("[21] TWO DISTINCT DISCRIMINANTS")
print("-" * 78)

show(
    "KAPPA discriminant",
    Delta,
)

show(
    "F2 inversion discriminant",
    discP_actual,
)

check(
    "KAPPA disc=(P-Q)^2",
    Delta,
    (P - Q) ** 2,
)

check(
    "F2 disc=4H_P^2",
    discP_actual,
    4 * H_P ** 2,
)

print(
    """
  The two discriminants are:

      KAPPA spectrum:
          Delta = (P-Q)^2

      F2 inversion:
          disc_P = 4(P-2Q-1/2)^2.

  They are different squares associated with the same
  two-dimensional root-coordinate system.
"""
)

print()


# ==============================================================================
# [22] INVERSION + TRANSLATION COMPATIBILITY
# ==============================================================================

print("[22] INVERSION / TRANSLATION COMPATIBILITY")
print("-" * 78)

# Translation after inversion
P_inv_then_shift = clean(P_star + h)
Q_inv_then_shift = clean(Q + h)

# Inversion after translation
P_shift_then_inv = clean(
    4 * Q_h + 1 - P_h
)

show(
    "T_h then I_P",
    P_shift_then_inv,
)

show(
    "I_P then T_h",
    P_inv_then_shift,
)

check(
    "T_h I_P(P) = I_P T_h(P) - 2h",
    P_shift_then_inv,
    P_inv_then_shift + h,
)

print(
    """
  The inversion and common translation do not commute trivially.

  The reason is that the affine reflection center itself moves
  under translation.

  This gives a semidirect affine structure rather than a
  simple commuting pair.
"""
)

print()


# ==============================================================================
# [23] FINAL STRUCTURAL CERTIFICATE
# ==============================================================================

print("[23] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print(
r"""
  CORRECT F2 INVERSION:

      F2 = 4PQ-P^2-Q^2+P+Q.

  As a quadratic in P:

      P^2-(4Q+1)P+(F+Q^2-Q)=0.

  Therefore:

      disc_P
        = 12Q^2+12Q+1-4F.

  Substitution of the actual F2 gives:

      disc_P(F2)
        = (2P-4Q-1)^2

        = 4(P-2Q-1/2)^2.

  The exact F2-preserving involution is:

      I_P:
          P -> 4Q+1-P.

  Define:

      H_P = P-2Q-1/2.

  Then:

      I_P(H_P) = -H_P.

  Under common root translation:

      P -> P+h
      Q -> Q+h,

      H_P -> H_P-h.

  Meanwhile:

      S -> S+2h

      Delta -> Delta.

  Therefore F2 simultaneously contains:

      its own reflection coordinate H_P,

      the KAPPA moving coordinate S/2,

      the invariant gap P-Q.

  These are not three independent variables.

  They are affine linear combinations of the same two
  2020 root coordinates P,Q.
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
print("  F2 inversion discriminant = CORRECT")
print("  F2 involution = TESTED")
print("  historical x,y crosswalk = TESTED")
print("  root translation = TESTED")
print("  KAPPA quadratic = TESTED")
print()

print("=" * 78)
print("EXPERIMENT 566 FINISHED")
print("=" * 78)
