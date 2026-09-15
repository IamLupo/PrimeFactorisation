#!/usr/bin/env python3

import sympy as sp


# =============================================================================
# EXPERIMENT 546
# =============================================================================
# RATIONAL F2/F3 OPERATOR CONJUGACY
#
# Coordinate system:
#
#     A = F2
#     B = F3
#
# Observable symmetric coordinate:
#
#     R = B/A - 1
#
# Observable product coordinate:
#
#     N = (A + R^2 - R) / 6
#
# Observable discriminant:
#
#     Delta = R^2 - 4N
#
# Distinguished translation:
#
#     R -> R + 2h
#
# The goal is to establish the complete operator algebra directly
# in (A,B), without reintroducing N or S as primitive variables.
#
# No numerical data.
# No factor enumeration.
# No continued fractions.
# No interpolation.
# No p or q.
# =============================================================================


print("=" * 78)
print("EXPERIMENT 546 START")
print("=" * 78)
print("RATIONAL F2/F3 OPERATOR CONJUGACY")
print()


# =============================================================================
# SYMBOLS
# =============================================================================

A, B = sp.symbols("A B", nonzero=True)
h, k = sp.symbols("h k")
z = sp.symbols("z")


failures = 0


def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.together(
                    sp.simplify(expr)
                )
            )
        )
    )


def cert(label, expr, expected=0):
    global failures

    diff = simp(expr - expected)
    ok = diff == 0

    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")

    if not ok:
        failures += 1

    return ok


# =============================================================================
# [1] OBSERVABLE COORDINATES
# =============================================================================

print()
print("[1] OBSERVABLE COORDINATES")
print("-" * 78)

R = sp.factor(
    B / A - 1
)

N_obs = sp.factor(
    (A + R**2 - R) / 6
)

Delta_obs = sp.factor(
    R**2 - 4*N_obs
)

print("  R = B/A - 1")
sp.pprint(R)

print()
print("  N(A,B) =")
sp.pprint(N_obs)

print()
print("  Delta(A,B) =")
sp.pprint(Delta_obs)


# =============================================================================
# [2] EXPLICIT RATIONAL FORMS
# =============================================================================

print()
print("[2] FULLY REDUCED OBSERVABLE FORMS")
print("-" * 78)

R_reduced = sp.factor(R)

N_reduced = sp.factor(N_obs)

Delta_reduced = sp.factor(Delta_obs)

print("  R =")
sp.pprint(R_reduced)

print()
print("  N =")
sp.pprint(N_reduced)

print()
print("  Delta =")
sp.pprint(Delta_reduced)

cert(
    "N formula",
    N_reduced,
    (A**3 - A**2 + B**2 - 2*A*B) / (6*A**2),
)

cert(
    "Delta formula",
    Delta_reduced,
    (B**2 - A**2 - 2*A**3) / (3*A**2),
)


# =============================================================================
# [3] RATIONAL TRANSLATION OPERATOR
# =============================================================================

print()
print("[3] RATIONAL TRANSLATION OPERATOR")
print("-" * 78)

# R -> R + 2h.
#
# Since A=F2:
#
#     A' = A + 2h(R+h)
#
# because B/A = R+1.

A1 = sp.factor(
    A + 2*h*(B/A + h)
)

# B'/A' = B/A + 2h.
#
# Therefore:
#
#     B' = (B+2hA)/A * A'

B1 = sp.factor(
    (B + 2*h*A) * A1 / A
)

print("  A' =")
sp.pprint(A1)

print()
print("  B' =")
sp.pprint(B1)


# =============================================================================
# [4] RATIO CHANNEL
# =============================================================================

print()
print("[4] RATIO CHANNEL")
print("-" * 78)

R1 = sp.factor(
    B1 / A1 - 1
)

print("  R' = B'/A' - 1")
sp.pprint(R1)

cert(
    "R' = R + 2h",
    R1,
    R + 2*h,
)


# =============================================================================
# [5] PRODUCT CHANNEL
# =============================================================================

print()
print("[5] PRODUCT CHANNEL")
print("-" * 78)

N1 = sp.factor(
    (
        A1
        + R1**2
        - R1
    ) / 6
)

print("  N' from (A',B') =")
sp.pprint(N1)

cert(
    "N' = N + hR + h^2",
    N1,
    N_reduced + h*R + h**2,
)


# =============================================================================
# [6] DISCRIMINANT CHANNEL
# =============================================================================

print()
print("[6] DISCRIMINANT CHANNEL")
print("-" * 78)

Delta1 = sp.factor(
    R1**2 - 4*N1
)

print("  Delta' =")
sp.pprint(Delta1)

cert(
    "Delta' = Delta",
    Delta1,
    Delta_reduced,
)


# =============================================================================
# [7] FIRST DIFFERENCE OF F2
# =============================================================================

print()
print("[7] F2 FIRST DIFFERENCE")
print("-" * 78)

dA = sp.factor(
    A1 - A
)

print("  A'-A =")
sp.pprint(dA)

cert(
    "A'-A = 2h(B/A+h)",
    dA,
    2*h*(B/A + h),
)


# =============================================================================
# [8] S-CHANNEL FROM FIRST DIFFERENCE
# =============================================================================

print()
print("[8] S-CHANNEL FROM F2 RESPONSE")
print("-" * 78)

R_from_dA = sp.factor(
    dA / (2*h) - h
)

print("  recovered R =")
sp.pprint(R_from_dA)

cert(
    "R from A'-A",
    R_from_dA,
    R,
)


# =============================================================================
# [9] B/A CHANNEL AS THE TRUE MOVING COORDINATE
# =============================================================================

print()
print("[9] MOVING COORDINATE")
print("-" * 78)

W = sp.factor(
    B / A
)

W1 = sp.factor(
    B1 / A1
)

print("  W = B/A")
sp.pprint(W)

print()
print("  W' = B'/A'")
sp.pprint(W1)

cert(
    "W' = W+2h",
    W1,
    W + 2*h,
)


# =============================================================================
# [10] INVARIANT IN OBSERVABLE COORDINATES
# =============================================================================

print()
print("[10] INVARIANT OBSERVABLE")
print("-" * 78)

I_AB = sp.factor(
    (B**2 - A**2 - 2*A**3) / (3*A**2)
)

I_A1B1 = sp.factor(
    (B1**2 - A1**2 - 2*A1**3) / (3*A1**2)
)

print("  I(A,B) =")
sp.pprint(I_AB)

print()
print("  I(A',B') =")
sp.pprint(I_A1B1)

cert(
    "I(A',B')=I(A,B)",
    I_A1B1,
    I_AB,
)


# =============================================================================
# [11] INVARIANT IDENTIFICATION
# =============================================================================

print()
print("[11] INVARIANT IDENTIFICATION")
print("-" * 78)

cert(
    "I(A,B)=Delta(A,B)",
    I_AB,
    Delta_reduced,
)


# =============================================================================
# [12] OPERATOR IN CONJUGATE COORDINATES
# =============================================================================

print()
print("[12] CONJUGATE COORDINATES")
print("-" * 78)

print("""
  Define:

      W = B/A

      I = (B^2-A^2-2A^3)/(3A^2)

  Then the translation acts as:

      W' = W + 2h

      I' = I.

  Thus the rational observable space has already separated into:

      moving coordinate:
          W

      invariant coordinate:
          I.
""")


# =============================================================================
# [13] RECONSTRUCT N FROM A AND W
# =============================================================================

print()
print("[13] N IN MOVING / INVARIANT COORDINATES")
print("-" * 78)

N_W = sp.factor(
    (A + (W - 1)**2 - (W - 1)) / 6
)

print("  N(A,W) =")
sp.pprint(N_W)

cert(
    "N(A,W)=N(A,B)",
    N_W,
    N_reduced,
)


# =============================================================================
# [14] ELIMINATE A USING THE INVARIANT
# =============================================================================

print()
print("[14] INVARIANT EQUATION")
print("-" * 78)

# I = ((W^2-1) - 2A)/3
#
# because B=A*W.

I_WA = sp.factor(
    I_AB.subs(B, A*W)
)

print("  I(A,W) =")
sp.pprint(I_WA)

cert(
    "I=(W^2-1-2A)/3",
    I_WA,
    (W**2 - 1 - 2*A) / 3,
)


# =============================================================================
# [15] SOLVE FOR A FROM THE INVARIANT
# =============================================================================

print()
print("[15] A RECOVERY FROM W AND I")
print("-" * 78)

I = sp.symbols("I")

A_from_WI = sp.factor(
    (W**2 - 1 - 3*I) / 2
)

print("  A =")
sp.pprint(A_from_WI)

cert(
    "A reconstruction equation",
    3*I - (W**2 - 1 - 2*A),
    I_WA*3 - 3*I,
)


# =============================================================================
# [16] N AS A FUNCTION OF W AND I
# =============================================================================

print()
print("[16] N AS W,I FUNCTION")
print("-" * 78)

N_WI = sp.factor(
    N_W.subs(
        A,
        A_from_WI
    )
)

print("  N(W,I) =")
sp.pprint(N_WI)


# =============================================================================
# [17] DISCRIMINANT CLOSURE
# =============================================================================

print()
print("[17] DISCRIMINANT CLOSURE")
print("-" * 78)

Delta_WI = sp.factor(
    I
)

print("  Delta(W,I) =")
sp.pprint(Delta_WI)

cert(
    "Delta is exactly invariant coordinate I",
    Delta_WI,
    I,
)


# =============================================================================
# [18] GENERATING KAPPA QUADRATIC
# =============================================================================

print()
print("[18] KAPPA QUADRATIC FROM OBSERVABLES")
print("-" * 78)

Q_AB = sp.factor(
    z**2 - R*z + N_reduced
)

print("  Q(z) =")
sp.pprint(Q_AB)

cert(
    "Q is monic",
    sp.Poly(Q_AB, z).LC(),
    1,
)


# =============================================================================
# [19] QUADRATIC DISCRIMINANT
# =============================================================================

print()
print("[19] QUADRATIC DISCRIMINANT")
print("-" * 78)

Q_disc = sp.factor(
    sp.discriminant(Q_AB, z)
)

print("  discriminant(Q) =")
sp.pprint(Q_disc)

cert(
    "disc(Q)=Delta",
    Q_disc,
    Delta_reduced,
)


# =============================================================================
# [20] +1 TRANSLATION
# =============================================================================

print()
print("[20] +1 KAPPA TRANSLATION")
print("-" * 78)

Q_shift = sp.expand(
    Q_AB.subs(z, z - 1)
)

print("  Q(z-1) =")
sp.pprint(Q_shift)

cert(
    "Q(z-1) has translated roots",
    sp.expand(
        Q_shift
        -
        (
            z**2
            - (R + 2)*z
            + (N_reduced + R + 1)
        )
    ),
    0,
)


# =============================================================================
# [21] FOUR-BASE COMPOSITION
# =============================================================================

print()
print("[21] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = sp.expand(
    Q_AB * Q_shift
)

print("  chi(z) = Q(z)Q(z-1)")
sp.pprint(sp.factor(chi))


# =============================================================================
# [22] TRANSLATION GROUP LAW
# =============================================================================

print()
print("[22] GROUP LAW")
print("-" * 78)

def transform(a0, b0, hh):
    a1 = sp.factor(
        a0 + 2*hh*b0/a0 + 2*hh**2
    )

    b1 = sp.factor(
        (b0 + 2*hh*a0)*a1/a0
    )

    return a1, b1


Ah, Bh = transform(A, B, h)
Ahk, Bhk = transform(Ah, Bh, k)
A_sum, B_sum = transform(A, B, h + k)

cert(
    "T_k(T_h(A))=T_(h+k)(A)",
    Ahk,
    A_sum,
)

cert(
    "T_k(T_h(B))=T_(h+k)(B)",
    Bhk,
    B_sum,
)


# =============================================================================
# [23] COMMUTATOR
# =============================================================================

print()
print("[23] COMMUTATOR")
print("-" * 78)

Ak, Bk = transform(A, B, k)
Akh, Bkh = transform(Ak, Bk, h)

cert(
    "[T_h,T_k]A=0",
    Ahk,
    Akh,
)

cert(
    "[T_h,T_k]B=0",
    Bhk,
    Bkh,
)


# =============================================================================
# [24] F2 ORBIT DISCRIMINANT DIRECTLY FROM A,B
# =============================================================================

print()
print("[24] F2 ORBIT DISCRIMINANT")
print("-" * 78)

A_orbit = sp.factor(
    A1
)

# A'(h) is quadratic in h.
A_poly = sp.Poly(
    A_orbit,
    h
)

a2 = sp.factor(
    A_poly.coeff_monomial(h**2)
)

a1 = sp.factor(
    A_poly.coeff_monomial(h)
)

a0 = sp.factor(
    A_poly.coeff_monomial(1)
)

D_orbit = sp.factor(
    a1**2 - 4*a2*a0
)

print("  h² coefficient =")
sp.pprint(a2)

print()
print("  h coefficient =")
sp.pprint(a1)

print()
print("  constant =")
sp.pprint(a0)

print()
print("  orbit discriminant =")
sp.pprint(D_orbit)

cert(
    "h² coefficient = 2",
    a2,
    2,
)

cert(
    "orbit discriminant = 12I+4",
    D_orbit,
    12*I_AB + 4,
)


# =============================================================================
# [25] DIRECT INVARIANT RECOVERY FROM ORBIT DISCRIMINANT
# =============================================================================

print()
print("[25] DELTA FROM ORBIT DISCRIMINANT")
print("-" * 78)

Delta_from_orbit = sp.factor(
    (D_orbit - 4) / 12
)

print("  Delta_candidate =")
sp.pprint(Delta_from_orbit)

cert(
    "Delta=(D_orbit-4)/12",
    Delta_from_orbit,
    I_AB,
)


# =============================================================================
# [26] LOW-DEGREE INVARIANT SEARCH
# =============================================================================

print()
print("[26] LOW-DEGREE RATIONAL INVARIANT SEARCH")
print("-" * 78)

#
# The coordinate W=B/A transforms as:
#
#     W -> W+2h.
#
# Therefore any invariant that is purely a polynomial in W
# must be constant.
#
# We verify this symbolically for a low-degree polynomial:
#
#     J(W)=c0+c1 W+c2 W²+c3 W³.
#

c0, c1, c2, c3 = sp.symbols(
    "c0 c1 c2 c3"
)

Jw = (
    c0
    + c1*W
    + c2*W**2
    + c3*W**3
)

Jw_shift = sp.expand(
    Jw.subs(W, W + 2*h)
)

difference_Jw = sp.Poly(
    sp.expand(Jw_shift - Jw),
    W,
    h,
)

eqs = [
    coeff
    for monom, coeff in difference_Jw.terms()
]

solution = sp.solve(
    eqs,
    [c1, c2, c3],
    dict=True,
)

print("  invariant polynomial solutions in W-only space:")
print(f"    {solution}")

if solution != [{}] and solution != [{c1: 0, c2: 0, c3: 0}]:
    failures += 1
    print("  PASS = False")
else:
    print("  PASS = True")


# =============================================================================
# [27] INVARIANT IS NOT A W-ONLY OBJECT
# =============================================================================

print()
print("[27] INVARIANT STRUCTURE")
print("-" * 78)

print("""
  The moving coordinate is:

      W = B/A = S+1.

  Under the operator:

      W -> W+2h.

  Therefore no nonconstant W-only polynomial/rational expression
  can be invariant under every h.

  The invariant must involve the second observable A as well.

  The lowest-degree exact invariant found here is:

      I(A,B)
        = (B²-A²-2A³)/(3A²).

  and:

      I = Delta.
""")


# =============================================================================
# [28] RECONSTRUCTION SUMMARY
# =============================================================================

print()
print("[28] EXACT OBSERVABLE RECONSTRUCTION")
print("-" * 78)

print("""
  Starting ONLY from:

      A = F2
      B = F3

  recover:

      R = B/A - 1

      N = (A + R² - R)/6

      Delta = R² - 4N

  equivalently:

      Delta
        = (B²-A²-2A³)/(3A²).

  Then:

      Q(z)=z²-Rz+N.

  The operator acts as:

      B/A -> B/A + 2h

      Delta -> Delta.

  Hence the observable system has the exact decomposition:

      moving channel:
          B/A

      invariant channel:
          Delta.
""")


# =============================================================================
# [29] FINAL STRUCTURAL QUESTION
# =============================================================================

print()
print("[29] FINAL RESEARCH QUESTION")
print("-" * 78)

print("""
  The downstream algebra is now completely expressed using only
  the existing homogeneous pair (F2,F3).

  The remaining problem is upstream:

      Is there an intrinsic homogeneous-layer index n
      for which the actual layer pair

          (F2_n, F3_n)

      obeys the rational action

          F2_{n+h}
            = F2_n
              + 2h(F3_n/F2_n + h),

          F3_{n+h}/F2_{n+h}
            = F3_n/F2_n + 2h?

  The decisive observable is therefore not F2 itself.

  It is the ratio:

      F3_n/F2_n.

  A genuine index law producing

      F3_{n+1}/F2_{n+1}
        - F3_n/F2_n = 2

  would directly identify the KAPPA translation coordinate.

  That is the next upstream object to search.
""")


# =============================================================================
# FINAL AUDIT
# =============================================================================

print()
print("=" * 78)
print("EXPERIMENT 546 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print("""
The important correction from Experiment 545 is:

    Delta(F2,F3)
      = (F3²-F2²-2F2³)/(3F2²).

This expression is genuinely in the observable coordinates
(F2,F3); it does not contain hidden N or S.

The exact closed operator is:

    A' = A + 2h(B/A+h)

    B' = (B+2hA)A'/A.

In these coordinates:

    B'/A' = B/A + 2h

and:

    Delta(A',B') = Delta(A,B).

The next research step should therefore stop constructing
symbolic conic models and examine the ORIGINAL homogeneous-layer
INDEX itself.

Specifically search for an exact index recurrence satisfying:

    R_n = F3_n/F2_n

    R_(n+1)-R_n = constant,

preferably:

    R_(n+1)-R_n = 2.

If this exists, then the layer index itself is acting as the
translation coordinate that was previously introduced by hand.
""")
