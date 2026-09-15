#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 552
# ==============================================================================
# 2020 FACTOR COORDINATES -> p,q -> S,N,DELTA -> KAPPA LAYERS
#
# Main purpose:
#
#   Stop treating S as an independently hidden quantity.
#
#   The historical coordinates give:
#
#       A0 = 2y - 2x + 3
#       B0 = 2y + 2x - 3
#
#   and, after branch normalization:
#
#       P = A0/r
#       Q = B0/r
#
#   Therefore:
#
#       N     = P Q
#       S     = P + Q
#       Delta = (P-Q)^2
#
#   The experiment rewrites the KAPPA layers DIRECTLY in A0,B0 and
#   x,y, then searches for factorization patterns.
#
# No numerical fitting.
# No divisor enumeration.
# No p/q reconstruction from N.
#
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 552 START")
print("=" * 78)
print("2020 FACTOR COORDINATES -> p,q -> KAPPA LAYERS")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

x, y, r, z = sp.symbols(
    "x y r z",
    nonzero=True,
)

failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

def exact(expr):
    return sp.factor(
        sp.cancel(
            sp.together(
                sp.expand(
                    sp.simplify(expr)
                )
            )
        )
    )


def certificate(label, actual, expected=0):
    global failures

    difference = exact(actual - expected)
    passed = difference == 0

    print(f"  {label}")
    print("    difference =")
    sp.pprint(difference)
    print(f"    PASS = {passed}")
    print()

    if not passed:
        failures += 1

    return passed


def show(label, expr):
    print(f"  {label} =")
    sp.pprint(exact(expr))
    print()


# ==============================================================================
# [1] HISTORICAL FACTOR COORDINATES
# ==============================================================================

print("[1] 2020 FACTOR COORDINATES")
print("-" * 78)

A0 = 2*y - 2*x + 3
B0 = 2*y + 2*x - 3

print("  A0 =")
sp.pprint(A0)

print()
print("  B0 =")
sp.pprint(B0)

print()

certificate(
    "A0 = 2y-2x+3",
    A0,
    2*y - 2*x + 3,
)

certificate(
    "B0 = 2y+2x-3",
    B0,
    2*y + 2*x - 3,
)


# ==============================================================================
# [2] NORMALIZED FACTORS
# ==============================================================================

print("[2] NORMALIZED FACTORS")
print("-" * 78)

P = exact(A0 / r)
Q = exact(B0 / r)

print("  P = A0/r =")
sp.pprint(P)

print()
print("  Q = B0/r =")
sp.pprint(Q)

print()

certificate(
    "A0 = rP",
    A0,
    r*P,
)

certificate(
    "B0 = rQ",
    B0,
    r*Q,
)


# ==============================================================================
# [3] MODERN SYMMETRIC VARIABLES DIRECTLY FROM x,y
# ==============================================================================

print("[3] N,S,DELTA DIRECTLY FROM x,y")
print("-" * 78)

N = exact(P*Q)
S = exact(P + Q)
Delta = exact((P-Q)**2)

show("N", N)
show("S", S)
show("Delta", Delta)

certificate(
    "N = P Q",
    N,
    P*Q,
)

certificate(
    "S = P+Q",
    S,
    P+Q,
)

certificate(
    "Delta = (P-Q)^2",
    Delta,
    (P-Q)**2,
)


# ==============================================================================
# [4] DIRECT HISTORICAL FORMS
# ==============================================================================

print("[4] DIRECT HISTORICAL FORMS")
print("-" * 78)

N_xy = exact(
    (
        4*y**2
        - (2*x-3)**2
    ) / r**2
)

S_xy = exact(
    4*y / r
)

Delta_xy = exact(
    4*(2*x-3)**2 / r**2
)

show("N(x,y)", N_xy)
show("S(x,y)", S_xy)
show("Delta(x,y)", Delta_xy)

certificate(
    "N(x,y)",
    N,
    N_xy,
)

certificate(
    "S(x,y)",
    S,
    S_xy,
)

certificate(
    "Delta(x,y)",
    Delta,
    Delta_xy,
)


# ==============================================================================
# [5] SUM / GAP / PRODUCT TRIAD
# ==============================================================================

print("[5] HISTORICAL SUM / GAP / PRODUCT")
print("-" * 78)

historical_sum = exact(A0 + B0)
historical_gap = exact(B0 - A0)
historical_product = exact(A0 * B0)

show("A0+B0", historical_sum)
show("B0-A0", historical_gap)
show("A0*B0", historical_product)

certificate(
    "A0+B0 = 4y",
    historical_sum,
    4*y,
)

certificate(
    "B0-A0 = 4x-6",
    historical_gap,
    4*x - 6,
)

certificate(
    "A0*B0 = 4v-1",
    historical_product,
    4*(y**2 - x**2 + 3*x - 2) - 1,
)


# ==============================================================================
# [6] Q(z) DIRECTLY FROM 2020 COORDINATES
# ==============================================================================

print("[6] KAPPA QUADRATIC FROM A0,B0")
print("-" * 78)

Qz = exact(
    z**2 - S*z + N
)

Qz_factor = exact(
    (z-P)*(z-Q)
)

show("Q(z)", Qz)
show("(z-P)(z-Q)", Qz_factor)

certificate(
    "Q(z)=(z-P)(z-Q)",
    Qz,
    Qz_factor,
)


# ==============================================================================
# [7] COEFFICIENTS OF Q(z)
# ==============================================================================

print("[7] Q(z) COEFFICIENTS")
print("-" * 78)

Q_poly = sp.Poly(
    sp.expand(Qz),
    z,
)

q2 = Q_poly.coeff_monomial(z**2)
q1 = Q_poly.coeff_monomial(z)
q0 = Q_poly.coeff_monomial(1)

print("  coefficient z^2 =")
sp.pprint(q2)

print()
print("  coefficient z =")
sp.pprint(q1)

print()
print("  constant =")
sp.pprint(q0)

print()

certificate(
    "Q is monic",
    q2,
    1,
)

certificate(
    "linear coefficient = -S",
    q1,
    -S,
)

certificate(
    "constant = N",
    q0,
    N,
)


# ==============================================================================
# [8] KAPPA F2 DIRECTLY IN P,Q
# ==============================================================================

print("[8] F2 IN P,Q")
print("-" * 78)

F2_PQ = exact(
    6*P*Q
    - (P+Q)**2
    + (P+Q)
)

F2_expected_PQ = exact(
    4*P*Q
    - P**2
    - Q**2
    + P
    + Q
)

show("F2", F2_PQ)
show("expanded P,Q form", F2_expected_PQ)

certificate(
    "F2 = 4PQ-P^2-Q^2+P+Q",
    F2_PQ,
    F2_expected_PQ,
)


# ==============================================================================
# [9] F2 FACTORIZATION STRUCTURE
# ==============================================================================

print("[9] F2 FACTORIZATION SEARCH")
print("-" * 78)

#
# Rewrite:
#
#   4PQ-P^2-Q^2
#     = 2PQ-(P-Q)^2
#
# hence:
#
#   F2 = 2PQ-(P-Q)^2+P+Q
#

F2_structural = exact(
    2*P*Q
    - (P-Q)**2
    + P
    + Q
)

show("2PQ-(P-Q)^2+P+Q", F2_structural)

certificate(
    "F2 structural decomposition",
    F2_PQ,
    F2_structural,
)


# ==============================================================================
# [10] F2 IN A0,B0
# ==============================================================================

print("[10] F2 DIRECTLY IN A0,B0")
print("-" * 78)

F2_AB = exact(
    2*A0*B0/r**2
    - (A0-B0)**2/r**2
    + (A0+B0)/r
)

show("F2(A0,B0)", F2_AB)

certificate(
    "F2 = 2AB/r^2 -(A-B)^2/r^2 +(A+B)/r",
    F2_AB,
    F2_PQ,
)


# ==============================================================================
# [11] F2 DIRECTLY IN x,y
# ==============================================================================

print("[11] F2 DIRECT x,y FORM")
print("-" * 78)

F2_xy = exact(
    F2_AB
)

show("F2(x,y)", F2_xy)

F2_expected_xy = exact(
    (
        8*y**2
        - 24*x**2
        + 72*x
        + 4*r*y
        - 54
    ) / r**2
    # Note:
    # the +4*r*y term produces 4y/r after division by r^2.
)

show("expanded expected F2(x,y)", F2_expected_xy)

certificate(
    "F2 direct x,y form",
    F2_xy,
    F2_expected_xy,
)


# ==============================================================================
# [12] IMPORTANT F2 CORRECTION
# ==============================================================================

print("[12] CORRECTED r=1 F2 CHECK")
print("-" * 78)

F2_r1 = exact(F2_xy.subs(r, 1))

show("F2(r=1)", F2_r1)

expected_r1 = (
    8*y**2
    - 24*x**2
    + 72*x
    + 4*y
    - 54
)

certificate(
    "F2(r=1)",
    F2_r1,
    expected_r1,
)

print(
    """
  The complete r=1 polynomial is:

      F2 =
          8y^2
        - 24x^2
        + 72x
        + 4y
        - 54.

  The y^2 term must NOT be omitted.
"""
)


# ==============================================================================
# [13] F2 = 2N + S - DELTA
# ==============================================================================

print("[13] F2 IN N,S,DELTA")
print("-" * 78)

F2_NSD = exact(
    2*N + S - Delta
)

certificate(
    "F2 = 2N+S-Delta",
    F2_PQ,
    F2_NSD,
)


# ==============================================================================
# [14] F3 IN P,Q
# ==============================================================================

print("[14] F3 IN P,Q")
print("-" * 78)

F3_PQ = exact(
    (S + 1) * F2_PQ
)

show("F3", F3_PQ)

certificate(
    "F3=(P+Q+1)F2",
    F3_PQ,
    (P+Q+1)*F2_PQ,
)


# ==============================================================================
# [15] F3 FACTOR STRUCTURE
# ==============================================================================

print("[15] F3 FACTOR STRUCTURE")
print("-" * 78)

F3_AB = exact(
    (A0+B0+r)/r * F2_AB
)

show("F3(A0,B0)", F3_AB)

certificate(
    "F3=(A0+B0+r)/r * F2",
    F3_AB,
    F3_PQ,
)


# ==============================================================================
# [16] F3/F2 WITHOUT RECONSTRUCTING S
# ==============================================================================

print("[16] F3/F2 DIRECT HISTORICAL SUM")
print("-" * 78)

#
# Formally cancel F2 only as a symbolic identity.
# This is not a numerical division procedure.
#

ratio_F3_F2 = exact(
    F3_PQ / F2_PQ
)

show("F3/F2", ratio_F3_F2)

certificate(
    "F3/F2 = S+1",
    ratio_F3_F2,
    S+1,
)

certificate(
    "F3/F2 = (A0+B0+r)/r",
    ratio_F3_F2,
    (A0+B0+r)/r,
)


# ==============================================================================
# [17] F2/F3 AS HISTORICAL LINEAR COORDINATE
# ==============================================================================

print("[17] DIRECT S FROM x,y")
print("-" * 78)

S_from_xy = exact(
    4*y/r
)

S_from_A0B0 = exact(
    (A0+B0)/r
)

show("S = 4y/r", S_from_xy)
show("S = (A0+B0)/r", S_from_A0B0)

certificate(
    "S from y",
    S_from_xy,
    S,
)

certificate(
    "S from A0+B0",
    S_from_A0B0,
    S,
)


# ==============================================================================
# [18] DELTA DIRECTLY FROM x
# ==============================================================================

print("[18] DIRECT DELTA FROM HISTORICAL x")
print("-" * 78)

Delta_from_x = exact(
    4*(2*x-3)**2/r**2
)

show("Delta", Delta_from_x)

certificate(
    "Delta = 4(2x-3)^2/r^2",
    Delta,
    Delta_from_x,
)


# ==============================================================================
# [19] N DIRECTLY FROM HISTORICAL CONIC
# ==============================================================================

print("[19] DIRECT N FROM HISTORICAL CONIC")
print("-" * 78)

v_xy = y**2 - x**2 + 3*x - 2

N_from_v = exact(
    (4*v_xy - 1)/r**2
)

show("N", N_from_v)

certificate(
    "N=(4v-1)/r^2",
    N,
    N_from_v,
)


# ==============================================================================
# [20] COMPLETE DIRECT CROSSWALK
# ==============================================================================

print("[20] COMPLETE 2020 -> 2026 CROSSWALK")
print("-" * 78)

print(
    """
  2020:

      A0 = 2y-2x+3
      B0 = 2y+2x-3

  Therefore:

      P = A0/r
      Q = B0/r

  and directly:

      N = A0 B0/r^2

      S = (A0+B0)/r

      Delta = (A0-B0)^2/r^2.

  Thus:

      S = 4y/r

  is already an explicit x,y quantity.

  No N-only inference is involved.
"""
)


# ==============================================================================
# [21] M1
# ==============================================================================

print("[21] M1 IN A0,B0")
print("-" * 78)

M1_PQ = exact(
    (P+1)*(Q+1)
)

M1_AB = exact(
    (A0+r)*(B0+r)/r**2
)

show("M1(P,Q)", M1_PQ)
show("M1(A0,B0)", M1_AB)

certificate(
    "M1=(A0+r)(B0+r)/r^2",
    M1_PQ,
    M1_AB,
)


# ==============================================================================
# [22] F2 UNDER HISTORICAL TRANSLATION
# ==============================================================================

print("[22] HISTORICAL TRANSLATION OF F2")
print("-" * 78)

h = sp.symbols("h")

A_shift = exact(A0 + h*r)
B_shift = exact(B0 + h*r)

P_shift = exact(A_shift/r)
Q_shift = exact(B_shift/r)

N_shift = exact(P_shift*Q_shift)
S_shift = exact(P_shift + Q_shift)

F2_shift = exact(
    6*N_shift
    - S_shift**2
    + S_shift
)

show("A0'", A_shift)
show("B0'", B_shift)
show("F2'", F2_shift)

certificate(
    "S'=S+2h",
    S_shift,
    S+2*h,
)

certificate(
    "N'=N+hS+h^2",
    N_shift,
    N+h*S+h**2,
)


# ==============================================================================
# [23] F2 TRANSLATION DIFFERENCE
# ==============================================================================

print("[23] F2 FIRST DIFFERENCE")
print("-" * 78)

F2_diff = exact(
    F2_shift - F2_PQ
)

show("F2'-F2", F2_diff)

certificate(
    "F2'-F2 = 2h(S+h+1)",
    F2_diff,
    2*h*(S+h+1),
)

S_from_F2_shift = exact(
    F2_diff/(2*h) - h - 1
)

certificate(
    "S recovered from F2 translation",
    S_from_F2_shift,
    S,
)


# ==============================================================================
# [24] HISTORICAL FACTOR TRANSLATION
# ==============================================================================

print("[24] FACTOR TRANSLATION")
print("-" * 78)

L_h = exact(
    A_shift*B_shift/r**2
)

show("L(h)", L_h)

certificate(
    "L(h)=N+hS+h^2",
    L_h,
    N+h*S+h**2,
)

certificate(
    "L(-z)=Q(z)",
    exact(L_h.subs(h, -z)),
    Qz,
)


# ==============================================================================
# [25] SEARCH FOR LOW-DEGREE FACTOR STRUCTURES
# ==============================================================================

print("[25] LOW-DEGREE FACTOR SEARCH")
print("-" * 78)

#
# Search several natural combinations.
#

candidates = {
    "F2": F2_PQ,
    "F3": F3_PQ,
    "F2-2N": exact(F2_PQ - 2*N),
    "F2-S": exact(F2_PQ - S),
    "F2+Delta": exact(F2_PQ + Delta),
    "F3/(S+1)": exact(F3_PQ/(S+1)),
    "F2/(P+Q+1)": exact(F2_PQ/(P+Q+1)),
}

for name, expr in candidates.items():
    print(f"  {name}")
    sp.pprint(exact(expr))
    print()


# ==============================================================================
# [26] SIMPLE LINEAR FACTORIZATIONS
# ==============================================================================

print("[26] FACTORIZATION OVER A0,B0")
print("-" * 78)

for name, expr in candidates.items():
    expr_ab = exact(
        expr.subs({
            # Replace P,Q structurally by A0/r,B0/r
            P: A0/r,
            Q: B0/r,
        })
    )

    print(f"  {name} ->")
    sp.pprint(sp.factor(expr_ab))
    print()


# ==============================================================================
# [27] DIRECT KAPPA QUADRATIC FACTORIZATION
# ==============================================================================

print("[27] ROOTS DIRECTLY FROM x,y")
print("-" * 78)

print(
    """
  Q(z) factors as:

      Q(z)
        = (z-P)(z-Q)

        = (z-A0/r)(z-B0/r).

  Therefore the two historical linear forms themselves are
  the two KAPPA roots after normalization.
"""
)

certificate(
    "Q factorization",
    Qz,
    (z-A0/r)*(z-B0/r),
)


# ==============================================================================
# [28] DISCRIMINANT
# ==============================================================================

print("[28] DISCRIMINANT FROM x,y")
print("-" * 78)

disc = exact(
    sp.discriminant(
        sp.Poly(Qz, z),
        z,
    )
)

show("disc(Q)", disc)

certificate(
    "disc(Q)=Delta",
    disc,
    Delta,
)


# ==============================================================================
# [29] FINAL STRUCTURAL RESULT
# ==============================================================================

print("[29] FINAL STRUCTURAL RESULT")
print("-" * 78)

print(
    """
  The direct coordinate map is:

      x,y
       |
       +--> A0 = 2y-2x+3
       +--> B0 = 2y+2x-3
                |
                v
             P=A0/r
             Q=B0/r
                |
        +-------+-------+
        |       |       |
        v       v       v
        N       S     Delta
        |       |       |
        |       |       |
        |       |       +--> (P-Q)^2
        |       |
        |       +----------> P+Q
        |
        +------------------> P Q

  Therefore:

      S = 4y/r

  and:

      Q(z)
        = z^2-Sz+N
        = (z-A0/r)(z-B0/r).

  This means the 2020 coordinate system already contains the
  complete symmetric KAPPA factor data.

  The remaining interesting problem is not:

      "Can S be expressed using x,y?"

  It can.

  The interesting problem is:

      "Which ORIGINAL homogeneous-layer expression, when written
       in A0 and B0, becomes simplest?"

  That is what the next experiment should investigate.
"""
)


# ==============================================================================
# FINAL AUDIT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 552 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print(
    """
The historical coordinates explicitly provide:

    P = A0/r
    Q = B0/r

    S = P+Q
    N = P Q
    Delta = (P-Q)^2.

Therefore the next experiment should NOT try to recover S from N.

Instead, take the full known homogeneous-layer algebra and rewrite
each layer in:

    P,Q

then in:

    A0,B0

and finally in:

    x,y.

The specific target is to discover whether the layer hierarchy has
simple factors or recurrences in the two historical root coordinates:

    A0
    B0

For example, search for structures involving:

    A0*B0
    A0+B0
    A0-B0
    (A0+1)(B0+1)
    A0^2
    B0^2
    A0+B0+r
    A0-B0.

This directly tests whether the 2020 construction already contains
the modern KAPPA coordinates at the algebraic root level.
"""
)
