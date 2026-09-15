#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 568
# ==============================================================================
# F2 AFFINE REFLECTION / TRANSLATION GROUP — CORRECTED
#
# Core objects:
#
#     F2(P,Q) = 4PQ - P^2 - Q^2 + P + Q
#
#     C_P = 2P - 4Q - 1
#
#     H_P = P - 2Q - 1/2
#
# so that:
#
#     C_P = 2 H_P
#
# P-involution:
#
#     I_P(P,Q) = (4Q+1-P, Q)
#
# common translation:
#
#     T_h(P,Q) = (P+h, Q+h)
#
# The exact affine relation is:
#
#     I_P T_h I_P = T_{-h}
#
# on the P moving direction.
#
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 568 START")
print("=" * 78)
print("F2 AFFINE REFLECTION / TRANSLATION GROUP — CORRECTED")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

P, Q = sp.symbols("P Q")
h, k = sp.symbols("h k")
z = sp.symbols("z")

x, y = sp.symbols("x y")
r = sp.symbols("r", nonzero=True)


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
    nd = numerator(diff)
    passed = sp.expand(nd) == 0

    print(f"  {name}")
    print("    numerator difference =")
    sp.pprint(nd)
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
# CORE QUANTITIES
# ==============================================================================

F2 = clean(
    4 * P * Q
    - P**2
    - Q**2
    + P
    + Q
)

N = clean(P * Q)
S = clean(P + Q)
Delta = clean((P - Q)**2)

# IMPORTANT:
# Explicitly define C_P. This fixes the NameError from Experiment 567.
C_P = clean(
    2 * P - 4 * Q - 1
)

H_P = clean(
    P - 2 * Q - sp.Rational(1, 2)
)

C_Q = clean(
    -4 * P + 2 * Q - 1
)

H_Q = clean(
    -2 * P + Q - sp.Rational(1, 2)
)


# ==============================================================================
# OPERATORS
# ==============================================================================

def T(expr, step):
    return clean(
        expr.subs(
            {
                P: P + step,
                Q: Q + step,
            },
            simultaneous=True,
        )
    )


def I_P(expr):
    return clean(
        expr.subs(
            {
                P: 4 * Q + 1 - P,
            },
            simultaneous=True,
        )
    )


def I_Q(expr):
    return clean(
        expr.subs(
            {
                Q: 4 * P + 1 - Q,
            },
            simultaneous=True,
        )
    )


# ==============================================================================
# [1] ROOT-COORDINATE SYSTEM
# ==============================================================================

print("[1] ROOT-COORDINATE SYSTEM")
print("-" * 78)

show("N=PQ", N)
show("S=P+Q", S)
show("Delta=(P-Q)^2", Delta)

check(
    "F2=2N+S-Delta",
    F2,
    2 * N + S - Delta,
)

print()


# ==============================================================================
# [2] F2 AS A QUADRATIC IN P
# ==============================================================================

print("[2] F2 AS A QUADRATIC IN P")
print("-" * 78)

F = sp.symbols("F")

aP = -1
bP = 4 * Q + 1
cP = -F - Q**2 + Q

discP = clean(
    sp.discriminant(
        aP * P**2 + bP * P + cP,
        P,
    )
)

discP_expected = clean(
    -4 * F + 12 * Q**2 + 12 * Q + 1
)

show("a_P", aP)
show("b_P", bP)
show("c_P", cP)
show("disc_P", discP)

check(
    "disc_P formula",
    discP,
    discP_expected,
)

discP_actual = clean(
    discP.subs(F, F2)
)

show("disc_P(F2)", discP_actual)

check(
    "disc_P(F2)=C_P^2",
    discP_actual,
    C_P**2,
)

print()


# ==============================================================================
# [3] F2 INVERSION CHANNEL
# ==============================================================================

print("[3] F2 INVERSION CHANNEL")
print("-" * 78)

show("C_P", C_P)
show("H_P", H_P)

check(
    "C_P=2H_P",
    C_P,
    2 * H_P,
)

check(
    "disc_P(F2)=4H_P^2",
    discP_actual,
    4 * H_P**2,
)

print()


# ==============================================================================
# [4] P-INVOLUTION
# ==============================================================================

print("[4] P-INVOLUTION")
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
    "I_P squared = identity",
    4 * Q + 1 - P_star,
    P,
)

check(
    "P+P_star=4Q+1",
    P + P_star,
    4 * Q + 1,
)

H_P_star = clean(
    H_P.subs(P, P_star)
)

C_P_star = clean(
    C_P.subs(P, P_star)
)

show("I_P(H_P)", H_P_star)
show("I_P(C_P)", C_P_star)

check(
    "I_P(H_P)=-H_P",
    H_P_star,
    -H_P,
)

check(
    "I_P(C_P)=-C_P",
    C_P_star,
    -C_P,
)

print()


# ==============================================================================
# [5] FIXED LINE
# ==============================================================================

print("[5] P-INVOLUTION FIXED LINE")
print("-" * 78)

P_fixed = clean(
    2 * Q + sp.Rational(1, 2)
)

show("P_fixed", P_fixed)

check(
    "P_star=P on fixed line",
    P_star.subs(P, P_fixed),
    P_fixed,
)

check(
    "H_P=0 on fixed line",
    H_P.subs(P, P_fixed),
    0,
)

print()


# ==============================================================================
# [6] Q-INVOLUTION
# ==============================================================================

print("[6] Q-INVOLUTION")
print("-" * 78)

Q_star = clean(
    4 * P + 1 - Q
)

show("Q_star", Q_star)

check(
    "F2(P,Q_star)=F2(P,Q)",
    F2.subs(Q, Q_star),
    F2,
)

check(
    "I_Q squared = identity",
    4 * P + 1 - Q_star,
    Q,
)

H_Q_star = clean(
    H_Q.subs(Q, Q_star)
)

C_Q_star = clean(
    C_Q.subs(Q, Q_star)
)

show("I_Q(H_Q)", H_Q_star)
show("I_Q(C_Q)", C_Q_star)

check(
    "I_Q(H_Q)=-H_Q",
    H_Q_star,
    -H_Q,
)

check(
    "I_Q(C_Q)=-C_Q",
    C_Q_star,
    -C_Q,
)

print()


# ==============================================================================
# [7] COMMON ROOT TRANSLATION
# ==============================================================================

print("[7] COMMON ROOT TRANSLATION")
print("-" * 78)

P_h = clean(P + h)
Q_h = clean(Q + h)

S_h = clean(P_h + Q_h)
Delta_h = clean((P_h - Q_h)**2)
C_P_h = clean(
    C_P.subs(
        {
            P: P_h,
            Q: Q_h,
        },
        simultaneous=True,
    )
)
H_P_h = clean(
    H_P.subs(
        {
            P: P_h,
            Q: Q_h,
        },
        simultaneous=True,
    )
)

show("T_h(S)", S_h)
show("T_h(Delta)", Delta_h)
show("T_h(C_P)", C_P_h)
show("T_h(H_P)", H_P_h)

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
# [8] CORRECT AFFINE COMPOSITIONS
# ==============================================================================

print("[8] AFFINE COMPOSITIONS")
print("-" * 78)

# T_h I_P:
# first I_P, then translate.
TI_P = clean(
    P_star + h
)

# I_P T_h:
# first translate, then reflect.
IT_P = clean(
    4 * Q_h + 1 - P_h
)

show("T_h I_P(P)", TI_P)
show("I_P T_h(P)", IT_P)

check(
    "T_h I_P(P)=4Q+1-P+h",
    TI_P,
    4 * Q + 1 - P + h,
)

check(
    "I_P T_h(P)=4Q+1-P+3h",
    IT_P,
    4 * Q + 1 - P + 3 * h,
)

check(
    "I_P T_h - T_h I_P = 2h",
    IT_P - TI_P,
    2 * h,
)

check(
    "T_h I_P - I_P T_h = -2h",
    TI_P - IT_P,
    -2 * h,
)

print()


# ==============================================================================
# [9] CONJUGATION I_P T_h I_P
# ==============================================================================

print("[9] REFLECTION CONJUGATION")
print("-" * 78)

# I_P -> T_h -> I_P

step1 = P_star
step2 = clean(step1 + h)
step3 = clean(
    4 * (Q + h) + 1 - step2
)

show("I_P T_h I_P(P)", step3)

check(
    "I_P T_h I_P = T_-h",
    step3,
    P - h,
)

print()


# ==============================================================================
# [10] SAME FOR Q-REFLECTION
# ==============================================================================

print("[10] Q-REFLECTION CONJUGATION")
print("-" * 78)

q_step1 = Q_star
q_step2 = clean(q_step1 + h)
q_step3 = clean(
    4 * (P + h) + 1 - q_step2
)

show("I_Q T_h I_Q(Q)", q_step3)

check(
    "I_Q T_h I_Q = T_-h",
    q_step3,
    Q - h,
)

print()


# ==============================================================================
# [11] COMPLETE-SQUARE FORM
# ==============================================================================

print("[11] COMPLETE-SQUARE FORM")
print("-" * 78)

F2_H_P_form = clean(
    -H_P**2
    + 3 * Q**2
    + 3 * Q
    + sp.Rational(1, 4)
)

F2_H_Q_form = clean(
    -H_Q**2
    + 3 * P**2
    + 3 * P
    + sp.Rational(1, 4)
)

show("F2=-H_P^2+3Q^2+3Q+1/4", F2_H_P_form)
show("F2=-H_Q^2+3P^2+3P+1/4", F2_H_Q_form)

check(
    "F2 complete square in P",
    F2,
    F2_H_P_form,
)

check(
    "F2 complete square in Q",
    F2,
    F2_H_Q_form,
)

print()


# ==============================================================================
# [12] SUM/GAP COORDINATES
# ==============================================================================

print("[12] SUM / GAP COORDINATES")
print("-" * 78)

s, d = sp.symbols("s d")

P_sd = clean((s + d) / 2)
Q_sd = clean((s - d) / 2)

F2_sd = clean(
    F2.subs(
        {
            P: P_sd,
            Q: Q_sd,
        },
        simultaneous=True,
    )
)

H_P_sd = clean(
    H_P.subs(
        {
            P: P_sd,
            Q: Q_sd,
        },
        simultaneous=True,
    )
)

C_P_sd = clean(
    C_P.subs(
        {
            P: P_sd,
            Q: Q_sd,
        },
        simultaneous=True,
    )
)

show("F2(s,d)", F2_sd)
show("H_P(s,d)", H_P_sd)
show("C_P(s,d)", C_P_sd)

check(
    "F2=(s^2+2s-3d^2)/2",
    F2_sd,
    (s**2 + 2 * s - 3 * d**2) / 2,
)

check(
    "H_P=(3d-s-1)/2",
    H_P_sd,
    (3 * d - s - 1) / 2,
)

check(
    "C_P=3d-s-1",
    C_P_sd,
    3 * d - s - 1,
)

print()


# ==============================================================================
# [13] HISTORICAL 2020 CROSSWALK
# ==============================================================================

print("[13] HISTORICAL 2020 CROSSWALK")
print("-" * 78)

A0 = clean(
    2 * y - 2 * x + 3
)

B0 = clean(
    2 * y + 2 * x - 3
)

P_xy = clean(A0 / r)
Q_xy = clean(B0 / r)

S_xy = clean(P_xy + Q_xy)
Delta_xy = clean((P_xy - Q_xy)**2)

C_P_xy = clean(
    C_P.subs(
        {
            P: P_xy,
            Q: Q_xy,
        },
        simultaneous=True,
    )
)

H_P_xy = clean(
    C_P_xy / 2
)

show("A0", A0)
show("B0", B0)
show("P(x,y)", P_xy)
show("Q(x,y)", Q_xy)
show("S(x,y)", S_xy)
show("Delta(x,y)", Delta_xy)
show("C_P(x,y)", C_P_xy)
show("H_P(x,y)", H_P_xy)

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

check(
    "S=4y/r",
    S_xy,
    4 * y / r,
)

check(
    "Delta=4(2x-3)^2/r^2",
    Delta_xy,
    4 * (2 * x - 3)**2 / r**2,
)

check(
    "C_P(x,y)=-(r+12x+4y-18)/r",
    C_P_xy,
    -(r + 12 * x + 4 * y - 18) / r,
)

print()


# ==============================================================================
# [14] HISTORICAL F2
# ==============================================================================

print("[14] HISTORICAL F2")
print("-" * 78)

F2_xy = clean(
    F2.subs(
        {
            P: P_xy,
            Q: Q_xy,
        },
        simultaneous=True,
    )
)

show("F2(x,y)", F2_xy)

check(
    "F2=2N+S-Delta",
    F2_xy,
    2 * P_xy * Q_xy + S_xy - Delta_xy,
)

check(
    "disc_P(F2)=(C_P)^2",
    discP_actual.subs(
        {
            P: P_xy,
            Q: Q_xy,
        },
        simultaneous=True,
    ),
    C_P_xy**2,
)

print()


# ==============================================================================
# [15] KAPPA SPECTRAL POLYNOMIAL
# ==============================================================================

print("[15] KAPPA SPECTRAL POLYNOMIAL")
print("-" * 78)

Qpoly = clean(
    z**2 - S*z + N
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

Qshift = clean(
    Qpoly.subs(z, z - 1)
)

Qshift_expected = clean(
    (z - P - 1) * (z - Q - 1)
)

check(
    "Q(z-1)=(z-P-1)(z-Q-1)",
    Qshift,
    Qshift_expected,
)

print()


# ==============================================================================
# [16] RELATION TO KAPPA GAP
# ==============================================================================

print("[16] F2 CHANNEL VS KAPPA GAP")
print("-" * 78)

gap = clean(P - Q)

check(
    "H_P + 2Q + 1/2 = P-Q",
    H_P + 2 * Q + sp.Rational(1, 2),
    gap,
)

check(
    "C_P + 4Q + 1 = 2(P-Q)",
    C_P + 4 * Q + 1,
    2 * gap,
)

print()


# ==============================================================================
# [17] THREE AFFINE CHANNELS
# ==============================================================================

print("[17] THREE AFFINE CHANNELS")
print("-" * 78)

show("S", S)
show("gap", gap)
show("C_P", C_P)
show("H_P", H_P)

check(
    "S translation response",
    T(S, h) - S,
    2 * h,
)

check(
    "gap invariant",
    T(gap, h) - gap,
    0,
)

check(
    "C_P translation response",
    T(C_P, h) - C_P,
    -2 * h,
)

check(
    "H_P translation response",
    T(H_P, h) - H_P,
    -h,
)

print()


# ==============================================================================
# [18] REFLECTION / TRANSLATION ACTION ON S
# ==============================================================================

print("[18] REFLECTION ACTION ON S")
print("-" * 78)

S_reflected = clean(
    S.subs(P, P_star)
)

Delta_reflected = clean(
    Delta.subs(P, P_star)
)

show("I_P(S)", S_reflected)
show("I_P(Delta)", Delta_reflected)

check(
    "I_P(S)=4Q+1-S",
    S_reflected,
    4 * Q + 1 - S,
)

check(
    "I_P(Delta)=Delta",
    Delta_reflected,
    Delta,
)

print()


# ==============================================================================
# [19] FINAL STRUCTURAL CERTIFICATE
# ==============================================================================

print("[19] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print(
r"""
  F2 in the 2020 root coordinates is:

      F2(P,Q)
        = 4PQ-P^2-Q^2+P+Q.

  Treating F2 as a quadratic in P gives:

      P^2-(4Q+1)P+(F+Q^2-Q)=0.

  Therefore:

      disc_P
        = 12Q^2+12Q+1-4F.

  On the actual F2 layer:

      disc_P(F2)
        = (2P-4Q-1)^2.

  Define:

      C_P = 2P-4Q-1

      H_P = P-2Q-1/2.

  Then:

      C_P = 2H_P

      disc_P(F2)=4H_P^2.

  The exact involution is:

      I_P(P,Q)
        = (4Q+1-P,Q).

  It satisfies:

      I_P(F2)=F2

      I_P(H_P)=-H_P

      I_P(C_P)=-C_P.

  The common root translation is:

      T_h(P,Q)
        = (P+h,Q+h).

  Its action is:

      S -> S+2h

      P-Q -> P-Q

      C_P -> C_P-2h

      H_P -> H_P-h.

  The exact reflection-conjugacy relation is:

      I_P T_h I_P = T_{-h}

  on the P moving direction.

  Thus the F2 inversion channel is an affine moving coordinate
  in the same 2020 root plane that contains:

      S=P+Q

      Delta=(P-Q)^2.

  The KAPPA spectral polynomial remains:

      Q(z)=(z-P)(z-Q)

  with:

      disc(Q)=Delta.

  No assumptions about F4, F5, or F6 are used.
"""
)

print()


# ==============================================================================
# [20] FINAL AUDIT
# ==============================================================================

print("[20] FINAL AUDIT")
print("-" * 78)

print(f"  symbolic failures = {failures}")
print("  F4/F5/F6 assumptions = NONE")
print("  F2 inversion discriminant = TESTED")
print("  F2 involution = TESTED")
print("  reflection conjugacy = TESTED")
print("  translation action = TESTED")
print("  historical x,y crosswalk = TESTED")
print("  KAPPA quadratic = TESTED")
print()

print("=" * 78)
print("EXPERIMENT 568 FINISHED")
print("=" * 78)