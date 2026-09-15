#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 555
# ==============================================================================
# 2020 ROOT-ORBIT OPERATOR -> HOMOGENEOUS-LAYER RECOGNITION
#
# IMPORTANT:
#   - symbolic only
#   - no numerical fitting
#   - no factor enumeration
#   - no p,q reconstruction
#   - no invented F4/F5/F6 formulas
#
# Primitive historical coordinates:
#
#     A = A0
#     B = B0
#
# with
#
#     A0 = 2y - 2x + 3
#     B0 = 2y + 2x - 3
#
# Normalize by r:
#
#     P = A/r
#     Q = B/r
#
# Exact orbit:
#
#     L_m = (A+mr)(B+mr)/r^2
#         = N + mS + m^2
#
# Known homogeneous layers:
#
#     F2 = 2N + S - Delta
#        = L0 + L1 - 1 - Delta
#
#     F3 = (S+1)F2
#        = (L1-L0)F2
#
# The unresolved object is the ACTUAL F4/F5/F6 from the original
# homogeneous-layer construction. They are deliberately not guessed.
#
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 555 START")
print("=" * 78)
print("2020 ROOT-ORBIT OPERATOR -> HOMOGENEOUS-LAYER RECOGNITION")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
m, h, z = sp.symbols("m h z")

failures = 0


# ==============================================================================
# EXACT HELPERS
# ==============================================================================

def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.together(
                sp.expand(
                    sp.simplify(expr)
                )
            )
        )
    )


def num(expr):
    """
    Clear rational denominators before polynomial operations.
    """
    return sp.expand(
        sp.together(expr).as_numer_denom()[0]
    )


def cert(name, actual, expected=0):
    global failures

    diff = simp(actual - expected)
    passed = diff == 0

    print(f"  {name}")
    print("    difference =")
    sp.pprint(diff)
    print(f"    PASS = {passed}")
    print()

    if not passed:
        failures += 1

    return passed


def show(name, expr):
    print(f"  {name} =")
    sp.pprint(simp(expr))
    print()


# ==============================================================================
# [1] HISTORICAL ROOT COORDINATES
# ==============================================================================

print("[1] HISTORICAL ROOT COORDINATES")
print("-" * 78)

P = simp(A / r)
Q = simp(B / r)

N = simp(P * Q)
S = simp(P + Q)
Delta = simp((P - Q)**2)
M1 = simp((P + 1) * (Q + 1))

show("P", P)
show("Q", Q)
show("N", N)
show("S", S)
show("Delta", Delta)
show("M1", M1)

cert("N = AB/r^2", N, A * B / r**2)
cert("S = (A+B)/r", S, (A + B) / r)
cert("Delta = (A-B)^2/r^2", Delta, (A - B)**2 / r**2)
cert("M1 = (A+r)(B+r)/r^2", M1, (A + r)*(B + r)/r**2)


# ==============================================================================
# [2] TRANSLATION ORBIT
# ==============================================================================

print("[2] 2020 TRANSLATION ORBIT")
print("-" * 78)


def L(j):
    return simp(
        (A + j*r) * (B + j*r) / r**2
    )


L0 = L(0)
L1 = L(1)
L2 = L(2)
L3 = L(3)
L4 = L(4)
L5 = L(5)

show("L0", L0)
show("L1", L1)
show("L2", L2)

cert("L0=N", L0, N)
cert("L1=N+S+1", L1, N + S + 1)
cert("L2=N+2S+4", L2, N + 2*S + 4)


# ==============================================================================
# [3] ORBIT POLYNOMIAL
# ==============================================================================

print("[3] ORBIT POLYNOMIAL")
print("-" * 78)

L_m = simp(
    (A + m*r) * (B + m*r) / r**2
)

show("L(m)", L_m)

cert(
    "L(m)=m^2+S*m+N",
    L_m,
    m**2 + S*m + N,
)


# ==============================================================================
# [4] ORBIT FIRST DIFFERENCE
# ==============================================================================

print("[4] ORBIT FIRST DIFFERENCE")
print("-" * 78)

D_L = simp(
    L(m + 1) - L(m)
)

show("L(m+1)-L(m)", D_L)

cert(
    "first difference",
    D_L,
    S + 2*m + 1,
)


# ==============================================================================
# [5] ORBIT SECOND DIFFERENCE
# ==============================================================================

print("[5] ORBIT SECOND DIFFERENCE")
print("-" * 78)

D2_L = simp(
    L(m + 2)
    - 2*L(m + 1)
    + L(m)
)

show("second difference", D2_L)

cert(
    "second difference = 2",
    D2_L,
    2,
)


# ==============================================================================
# [6] CENTRAL DIFFERENCE
# ==============================================================================

print("[6] CENTRAL DIFFERENCE")
print("-" * 78)

central = simp(
    L(m + h) - L(m - h)
)

cert(
    "central difference",
    central,
    2*h*(S + 2*m),
)


# ==============================================================================
# [7] S RECOVERY FROM ORBIT
# ==============================================================================

print("[7] S RECOVERY FROM ROOT ORBIT")
print("-" * 78)

S_orbit = simp(
    L1 - L0 - 1
)

show("S_orbit", S_orbit)

cert(
    "S=L1-L0-1",
    S_orbit,
    S,
)


# ==============================================================================
# [8] KAPPA QUADRATIC
# ==============================================================================

print("[8] KAPPA QUADRATIC")
print("-" * 78)

Q_orbit = simp(
    L_m.subs(m, -z)
)

Q_target = simp(
    z**2 - S*z + N
)

show("Q_orbit(z)", Q_orbit)

cert(
    "L(-z)=Q(z)",
    Q_orbit,
    Q_target,
)

Q_factor = simp(
    (z - A/r) * (z - B/r)
)

cert(
    "Q=(z-A/r)(z-B/r)",
    Q_orbit,
    Q_factor,
)


# ==============================================================================
# [9] KAPPA TRANSLATION
# ==============================================================================

print("[9] KAPPA +1 TRANSLATION")
print("-" * 78)

Q_shift = simp(
    Q_orbit.subs(z, z - 1)
)

Q_shift_expected = simp(
    (z - (A+r)/r)
    * (z - (B+r)/r)
)

cert(
    "Q(z-1)",
    Q_shift,
    Q_shift_expected,
)


# ==============================================================================
# [10] F2
# ==============================================================================

print("[10] F2 ROOT-ORBIT REPRESENTATION")
print("-" * 78)

F2 = simp(
    6*N - S**2 + S
)

F2_alt = simp(
    L0 + L1 - 1 - Delta
)

show("F2", F2)
show("F2_alt", F2_alt)

cert(
    "F2 = L0+L1-1-Delta",
    F2_alt,
    F2,
)

cert(
    "F2 = 2N+S-Delta",
    F2,
    2*N + S - Delta,
)


# ==============================================================================
# [11] F3
# ==============================================================================

print("[11] F3 ROOT-ORBIT REPRESENTATION")
print("-" * 78)

F3 = simp(
    (S + 1) * F2
)

F3_alt = simp(
    (L1 - L0) * F2
)

show("F3", F3)
show("F3_alt", F3_alt)

cert(
    "F3=(L1-L0)F2",
    F3_alt,
    F3,
)


# ==============================================================================
# [12] F3/F2
# ==============================================================================

print("[12] F3/F2")
print("-" * 78)

ratio23 = simp(
    F3/F2
)

cert(
    "F3/F2 = L1-L0",
    ratio23,
    L1-L0,
)

cert(
    "F3/F2 = S+1",
    ratio23,
    S+1,
)


# ==============================================================================
# [13] F2 TRANSLATION
# ==============================================================================

print("[13] F2 TRANSLATION UNDER ROOT ORBIT")
print("-" * 78)

A_h = A + h*r
B_h = B + h*r

N_h = simp(A_h*B_h/r**2)
S_h = simp((A_h+B_h)/r)
Delta_h = simp((A_h-B_h)**2/r**2)

F2_h = simp(
    6*N_h - S_h**2 + S_h
)

cert(
    "N_h=N+hS+h^2",
    N_h,
    N + h*S + h**2,
)

cert(
    "S_h=S+2h",
    S_h,
    S + 2*h,
)

cert(
    "Delta_h=Delta",
    Delta_h,
    Delta,
)

cert(
    "F2_h-F2",
    F2_h-F2,
    2*h*(S+h+1),
)


# ==============================================================================
# [14] F2 ORBIT DISCRIMINANT
# ==============================================================================

print("[14] F2 ORBIT DISCRIMINANT")
print("-" * 78)

F2_num, F2_den = sp.together(F2_h).as_numer_denom()

F2_poly = sp.Poly(
    sp.expand(F2_num),
    h,
)

c2 = simp(
    F2_poly.coeff_monomial(h**2)
    / F2_den
)

c1 = simp(
    F2_poly.coeff_monomial(h)
    / F2_den
)

c0 = simp(
    F2_poly.coeff_monomial(1)
    / F2_den
)

D_F2 = simp(
    c1**2 - 4*c2*c0
)

show("c2", c2)
show("c1", c1)
show("c0", c0)
show("disc_h(F2)", D_F2)

cert("F2 h^2 coefficient", c2, 2)
cert("F2 h coefficient", c1, 2*(S+1))
cert("F2 constant", c0, F2)
cert(
    "disc(F2 orbit)=12Delta+4",
    D_F2,
    12*Delta + 4,
)


# ==============================================================================
# [15] F2/F3 AS ROOT-ORBIT DATA
# ==============================================================================

print("[15] F2/F3 AS ROOT-ORBIT DATA")
print("-" * 78)

print(
    r"""
  The exact known layer structure is:

      F2 = L0 + L1 - 1 - Delta

      F3 = (L1-L0)F2.

  Therefore:

      F3/F2 = L1-L0 = S+1.

  This is an important distinction:

      L1-L0
  is additive orbit information,

      F3/F2
  is multiplicative layer information.

  They happen to coincide exactly at n=2.
"""
)


# ==============================================================================
# [16] ABSTRACT HIGHER-LAYER TEST
# ==============================================================================

print("[16] ABSTRACT HIGHER-LAYER TEST")
print("-" * 78)

#
# Do NOT define F4/F5/F6.
#
# The purpose of this experiment is to derive the exact identities
# that genuine formulas must satisfy.
#

F4_actual = sp.Symbol("F4_actual")
F5_actual = sp.Symbol("F5_actual")
F6_actual = sp.Symbol("F6_actual")

E2 = simp(
    F4_actual*F2
    - F3**2
    - 2*F3*F2
)

E3 = simp(
    F5_actual*F3
    - F4_actual**2
    - 2*F4_actual*F3
)

E4 = simp(
    F6_actual*F4_actual
    - F5_actual**2
    - 2*F5_actual*F4_actual
)

print("  division-free defects:")
print()
print("    E2 = F4*F2 - F3^2 - 2F3F2")
print("    E3 = F5*F3 - F4^2 - 2F4F3")
print("    E4 = F6*F4 - F5^2 - 2F5F4")
print()


# ==============================================================================
# [17] PREDICTED CONSEQUENCES — NOT ASSERTED AS ACTUAL LAYERS
# ==============================================================================

print("[17] PREDICTED CONSEQUENCES")
print("-" * 78)

F4_pred = simp(
    F3 * (S + 3)
)

F5_pred = simp(
    F4_pred * (S + 5)
)

F6_pred = simp(
    F5_pred * (S + 7)
)

show("F4_pred", F4_pred)
show("F5_pred", F5_pred)
show("F6_pred", F6_pred)

print(
    r"""
  These are hypotheses generated by the root-orbit/index model.

  They are NOT treated as genuine homogeneous layers.

  Genuine formulas must be substituted separately.
"""
)


# ==============================================================================
# [18] PREDICTED VALUES SATISFY THE DIVISION-FREE DEFECTS
# ==============================================================================

print("[18] HYPOTHETICAL LADDER CONSISTENCY")
print("-" * 78)

cert(
    "predicted E2=0",
    F4_pred*F2 - F3**2 - 2*F3*F2,
    0,
)

cert(
    "predicted E3=0",
    F5_pred*F3 - F4_pred**2 - 2*F4_pred*F3,
    0,
)

cert(
    "predicted E4=0",
    F6_pred*F4_pred - F5_pred**2 - 2*F5_pred*F4_pred,
    0,
)


# ==============================================================================
# [19] ROOT-ORBIT ADDITIVE LADDER
# ==============================================================================

print("[19] ROOT-ORBIT ADDITIVE LADDER")
print("-" * 78)

for j in range(0, 6):
    lhs = simp(
        L(j+1) - L(j)
    )
    rhs = simp(
        S + 2*j + 1
    )

    cert(
        f"L{j+1}-L{j}=S+{2*j+1}",
        lhs,
        rhs,
    )


# ==============================================================================
# [20] ROOT-ORBIT / LAYER INDEX COMPARISON
# ==============================================================================

print("[20] ROOT-ORBIT / LAYER INDEX COMPARISON")
print("-" * 78)

print(
    r"""
  Root-orbit law:

      L_{m+1}-L_m = S+2m+1.

  Proposed homogeneous-layer law:

      F_{n+1}/F_n = S+2n-3.

  At n=2:

      F3/F2 = S+1
            = L1-L0.

  Therefore the known n=2 layer ratio coincides with the
  first root-orbit difference.

  The next question is whether this coincidence continues:

      F4/F3 ?= L2-L1 = S+3

      F5/F4 ?= L3-L2 = S+5

      F6/F5 ?= L4-L3 = S+7.

  This is the exact structural comparison to make once the
  genuine higher homogeneous-layer formulas are available.
"""
)


# ==============================================================================
# [21] RECONSTRUCT P,Q FROM ROOT POLYNOMIAL — SYMBOLIC ONLY
# ==============================================================================

print("[21] ROOT POLYNOMIAL STRUCTURE")
print("-" * 78)

Q_poly = sp.expand(
    z**2 - S*z + N
)

Q_factored = sp.expand(
    (z-P)*(z-Q)
)

cert(
    "Q(z)=(z-P)(z-Q)",
    Q_poly,
    Q_factored,
)

disc_Q = simp(
    sp.discriminant(Q_poly, z)
)

cert(
    "disc(Q)=Delta",
    disc_Q,
    Delta,
)


# ==============================================================================
# [22] OPERATOR INTERPRETATION
# ==============================================================================

print("[22] OPERATOR INTERPRETATION")
print("-" * 78)

print(
    r"""
  The exact 2020 root-coordinate operator is:

      T_h:
          A -> A + h r
          B -> B + h r.

  Therefore:

      L_m
        = T_m(AB/r^2)

        = N + mS + m^2.

  The operator supplies:

      moving channel:
          A+B

      invariant channel:
          A-B

      product channel:
          AB.

  The modern KAPPA quadratic is simply:

      Q(z)=L_{-z}.

  The known F2 layer is:

      F2=L0+L1-1-Delta.

  The known F3 layer is:

      F3=(L1-L0)F2.

  Thus the precise unresolved question is:

      Do higher original homogeneous layers continue to be
      generated by this same A,B translation orbit?
"""
)


# ==============================================================================
# [23] FINAL AUDIT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 555 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print(
    r"""
The exact root-coordinate structure is now:

    A0 = 2y-2x+3
    B0 = 2y+2x-3

    L_m = (A0+mr)(B0+mr)/r^2

    L_m = N+mS+m^2

    L_{m+1}-L_m = S+2m+1

    L_{m+2}-2L_{m+1}+L_m = 2

and:

    Q(z)=L_{-z}.

The known homogeneous layers satisfy:

    F2 = L0+L1-1-Delta

    F3 = (L1-L0)F2.

This establishes an exact correspondence at the first
nontrivial layer.

The next decisive test requires the GENUINE original formulas
for F4, F5 and F6.

For those formulas, test only:

    F4/F3 - (L2-L1)

    F5/F4 - (L3-L2)

    F6/F5 - (L4-L3)

or, without division:

    F4*F2 - F3^2 - 2F3F2

    F5*F3 - F4^2 - 2F4F3

    F6*F4 - F5^2 - 2F5F4.

A zero result would identify the homogeneous-layer index with
the same translation coordinate already present in the 2020
factor-coordinate orbit.

This experiment intentionally does not manufacture those
higher layers.
"""
)
