#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 562
# ==============================================================================
# 2020 ROOT COORDINATES -> ACTUAL HOMOGENEOUS-LAYER STEP TEST
#
# Goal:
#   Work directly with the historical factor coordinates A,B and test the
#   unique continuation implied by the already-proven F2/F3 correspondence.
#
# IMPORTANT:
#   F4_ACTUAL, F5_ACTUAL, F6_ACTUAL are placeholders.
#   Replace ONLY those definitions with genuine formulas from the original
#   homogeneous-layer construction.
#
#   Nothing below invents or substitutes them.
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 562 START")
print("=" * 78)
print("2020 ROOT COORDINATES -> ACTUAL HOMOGENEOUS-LAYER STEP TEST")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
F4_ACTUAL, F5_ACTUAL, F6_ACTUAL = sp.symbols(
    "F4_actual F5_actual F6_actual"
)


# ==============================================================================
# HELPERS
# ==============================================================================

failures = 0


def clean(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def numerator(expr):
    return sp.factor(sp.together(expr).as_numer_denom()[0])


def check(name, lhs, rhs=0):
    global failures

    diff = clean(lhs - rhs)
    num = numerator(diff)

    passed = num == 0

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

print("[1] 2020 ROOT COORDINATES")
print("-" * 78)

P = clean(A / r)
Q = clean(B / r)

N = clean(A * B / r**2)
S = clean((A + B) / r)
Delta = clean((A - B)**2 / r**2)

show("P", P)
show("Q", Q)
show("N", N)
show("S", S)
show("Delta", Delta)

check("N=PQ", N, P * Q)
check("S=P+Q", S, P + Q)
check("Delta=(P-Q)^2", Delta, (P - Q)**2)

print()


# ==============================================================================
# [2] ROOT-TRANSLATION ORBIT
# ==============================================================================

print("[2] ROOT-TRANSLATION ORBIT")
print("-" * 78)


def L(m):
    return clean(
        (A + m*r) * (B + m*r) / r**2
    )


def D(m):
    return clean(
        L(m + 1) - L(m)
    )


for m in range(6):
    check(
        f"L{m}=N+{m}S+{m*m}",
        L(m),
        N + m*S + m**2,
    )

for m in range(5):
    check(
        f"D{m}=S+{2*m+1}",
        D(m),
        S + 2*m + 1,
    )

for m in range(4):
    check(
        f"D{m+1}-D{m}=2",
        D(m + 1) - D(m),
        2,
    )

print()


# ==============================================================================
# [3] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[3] KNOWN HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = clean(
    2*N + S - Delta
)

F3 = clean(
    (S + 1) * F2
)

check("F2=2N+S-Delta", F2, 2*N + S - Delta)
check("F3=(S+1)F2", F3, (S + 1) * F2)
check("F3=D0*F2", F3, D(0) * F2)

print()


# ==============================================================================
# [4] ROOT-ORBIT REPRESENTATION
# ==============================================================================

print("[4] ROOT-ORBIT REPRESENTATION")
print("-" * 78)

check(
    "F2=L0+L1-1-Delta",
    F2,
    L(0) + L(1) - 1 - Delta,
)

check(
    "F3=(L1-L0)F2",
    F3,
    (L(1) - L(0)) * F2,
)

print()


# ==============================================================================
# [5] UNIQUE STEP OPERATOR
# ==============================================================================

print("[5] UNIQUE ROOT-ORBIT STEP OPERATOR")
print("-" * 78)

G2 = clean(D(0))
G3 = clean(D(1))
G4 = clean(D(2))
G5 = clean(D(3))

show("G2", G2)
show("G3", G3)
show("G4", G4)
show("G5", G5)

check("G2=S+1", G2, S + 1)
check("G3=S+3", G3, S + 3)
check("G4=S+5", G4, S + 5)
check("G5=S+7", G5, S + 7)

print()


# ==============================================================================
# [6] UNIQUE HIGHER-LAYER PREDICTIONS
# ==============================================================================

print("[6] UNIQUE HIGHER-LAYER PREDICTIONS")
print("-" * 78)

F4_EXPECTED = clean(G3 * F3)
F5_EXPECTED = clean(G4 * F4_EXPECTED)
F6_EXPECTED = clean(G5 * F5_EXPECTED)

show("F4_EXPECTED", F4_EXPECTED)
show("F5_EXPECTED", F5_EXPECTED)
show("F6_EXPECTED", F6_EXPECTED)

check(
    "F4_EXPECTED=(S+3)F3",
    F4_EXPECTED,
    (S + 3) * F3,
)

check(
    "F5_EXPECTED=(S+5)F4_EXPECTED",
    F5_EXPECTED,
    (S + 5) * F4_EXPECTED,
)

check(
    "F6_EXPECTED=(S+7)F5_EXPECTED",
    F6_EXPECTED,
    (S + 7) * F5_EXPECTED,
)

print()


# ==============================================================================
# [7] PREDICTIONS DIRECTLY IN A,B,r
# ==============================================================================

print("[7] PREDICTIONS DIRECTLY IN A,B,r")
print("-" * 78)

show("F4_EXPECTED(A,B,r)", F4_EXPECTED)
show("F5_EXPECTED(A,B,r)", F5_EXPECTED)
show("F6_EXPECTED(A,B,r)", F6_EXPECTED)

print()


# ==============================================================================
# [8] ACTUAL LAYER PLACEHOLDERS
# ==============================================================================

print("[8] ACTUAL HIGHER-LAYER INPUT")
print("-" * 78)

print(
"""
  Replace these three symbols with the REAL formulas:

      F4_ACTUAL = ...
      F5_ACTUAL = ...
      F6_ACTUAL = ...

  They must come from the original homogeneous-layer construction.

  They must NOT be copied from F4_EXPECTED/F5_EXPECTED/F6_EXPECTED.
"""
)

show("F4_ACTUAL", F4_ACTUAL)
show("F5_ACTUAL", F5_ACTUAL)
show("F6_ACTUAL", F6_ACTUAL)

print()


# ==============================================================================
# [9] DIRECT ROOT-STEP TESTS
# ==============================================================================

print("[9] DIRECT ROOT-STEP TESTS")
print("-" * 78)

R4 = clean(
    F4_ACTUAL - G3 * F3
)

R5 = clean(
    F5_ACTUAL - G4 * F4_ACTUAL
)

R6 = clean(
    F6_ACTUAL - G5 * F5_ACTUAL
)

show("F4_ACTUAL-D1*F3", R4)
show("F5_ACTUAL-D2*F4_ACTUAL", R5)
show("F6_ACTUAL-D3*F5_ACTUAL", R6)

print()


# ==============================================================================
# [10] RATIO-FREE STEP TESTS
# ==============================================================================

print("[10] RATIO-FREE STEP TESTS")
print("-" * 78)

check(
    "F4/F3 = S+3",
    F4_ACTUAL,
    (S + 3) * F3,
)

check(
    "F5/F4 = S+5",
    F5_ACTUAL,
    (S + 5) * F4_ACTUAL,
)

check(
    "F6/F5 = S+7",
    F6_ACTUAL,
    (S + 7) * F5_ACTUAL,
)

print()


# ==============================================================================
# [11] DIVISION-FREE SECOND-ORDER TESTS
# ==============================================================================

print("[11] DIVISION-FREE SECOND-ORDER TESTS")
print("-" * 78)

E2 = clean(
    F4_ACTUAL * F2
    - F3**2
    - 2 * F3 * F2
)

E3 = clean(
    F5_ACTUAL * F3
    - F4_ACTUAL**2
    - 2 * F4_ACTUAL * F3
)

E4 = clean(
    F6_ACTUAL * F4_ACTUAL
    - F5_ACTUAL**2
    - 2 * F5_ACTUAL * F4_ACTUAL
)

show("E2", E2)
show("E3", E3)
show("E4", E4)

print(
"""
  Expected if the original homogeneous hierarchy continues the
  root-orbit step:

      E2 = 0
      E3 = 0
      E4 = 0.
"""
)

print()


# ==============================================================================
# [12] RATIO-LADDER FORM WITHOUT DIVISION
# ==============================================================================

print("[12] RATIO-LADDER FORM")
print("-" * 78)

# The quantity
#
#     F_(n+2) F_n - F_(n+1)^2
#
# is the numerator corresponding to R_(n+1)-R_n.
#
# The KAPPA translation hypothesis says this equals:
#
#     2 F_(n+1) F_n.

Q2_defect = clean(
    F4_ACTUAL * F2 - F3**2
)

Q3_defect = clean(
    F5_ACTUAL * F3 - F4_ACTUAL**2
)

Q4_defect = clean(
    F6_ACTUAL * F4_ACTUAL - F5_ACTUAL**2
)

show(
    "F4F2-F3^2",
    Q2_defect,
)

show(
    "F5F3-F4^2",
    Q3_defect,
)

show(
    "F6F4-F5^2",
    Q4_defect,
)

print()


# ==============================================================================
# [13] ROOT-ORBIT CENTRAL DIFFERENCE
# ==============================================================================

print("[13] ROOT-ORBIT CENTRAL DIFFERENCE")
print("-" * 78)

m = sp.symbols("m")

Lm = L(m)

check(
    "L(h)-L(-h)=2hS",
    Lm.subs(m, sp.Symbol("h")) - Lm.subs(m, -sp.Symbol("h")),
    2 * sp.Symbol("h") * S,
)

check(
    "second difference = 2",
    L(2) - 2*L(1) + L(0),
    2,
)

print()


# ==============================================================================
# [14] KAPPA SPECTRAL POLYNOMIAL
# ==============================================================================

print("[14] KAPPA SPECTRAL POLYNOMIAL")
print("-" * 78)

z = sp.symbols("z")

Q = clean(
    z**2 - S*z + N
)

Q_factored = clean(
    (z - P) * (z - Q)
)

show("Q(z)", Q)
show("Q factorization", Q_factored)

check(
    "Q=(z-P)(z-Q)",
    Q,
    Q_factored,
)

check(
    "disc(Q)=Delta",
    sp.discriminant(Q, z),
    Delta,
)

print()


# ==============================================================================
# [15] FOUR-BASE COMPOSITION
# ==============================================================================

print("[15] FOUR-BASE COMPOSITION")
print("-" * 78)

chi = clean(
    Q * Q.subs(z, z - 1)
)

check(
    "Q(z-1) roots translated by +1",
    Q.subs(z, z - 1),
    (z - P - 1) * (z - Q - 1),
)

show("chi(z)", chi)

print()


# ==============================================================================
# [16] EXACT OBSERVABLE INTERPRETATION
# ==============================================================================

print("[16] EXACT OBSERVABLE INTERPRETATION")
print("-" * 78)

print(
r"""
  The proven part is now:

      A0,B0
        |
        v
      L_m=(A0+mr)(B0+mr)/r^2
        |
        +--> N
        +--> S
        +--> Delta
        |
        v
      D_m=L_(m+1)-L_m
        = S+2m+1

  and the known homogeneous pair satisfies:

      F2 = L0+L1-1-Delta

      F3 = D0 F2.

  Therefore:

      F3/F2 = D0 = S+1.

  The only open statement is whether the ORIGINAL higher layers
  continue:

      F4 = D1 F3

      F5 = D2 F4

      F6 = D3 F5.

  This experiment tests exactly that statement once genuine
  F4/F5/F6 formulas are supplied.
"""
)

print()


# ==============================================================================
# [17] FINAL AUDIT
# ==============================================================================

print("[17] FINAL AUDIT")
print("-" * 78)

print(f"  built-in symbolic failures = {failures}")
print("  F4/F5/F6 genuine formulas supplied = NO")
print("  higher-layer conclusion = OPEN")
print()

print("=" * 78)
print("EXPERIMENT 562 FINISHED")
print("=" * 78)
