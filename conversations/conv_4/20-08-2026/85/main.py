#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 563
# ==============================================================================
# 2020 ROOT COORDINATES -> HOMOGENEOUS LAYER CONTINUATION
#
# Important:
#   P_root / Q_root are the historical normalized factor coordinates.
#   Qpoly is the KAPPA spectral polynomial.
#
#   F4_actual / F5_actual / F6_actual remain placeholders until genuine
#   higher homogeneous-layer formulas are supplied.
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 563 START")
print("=" * 78)
print("2020 ROOT COORDINATES -> HOMOGENEOUS LAYER CONTINUATION")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r, z, h = sp.symbols(
    "A B r z h",
    nonzero=True,
)

F4_actual, F5_actual, F6_actual = sp.symbols(
    "F4_actual F5_actual F6_actual"
)

failures = 0


# ==============================================================================
# EXACT RATIONAL HELPERS
# ==============================================================================

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
# [1] HISTORICAL ROOT COORDINATES
# ==============================================================================

print("[1] HISTORICAL ROOT COORDINATES")
print("-" * 78)

P_root = clean(A / r)
Q_root = clean(B / r)

N = clean(A * B / r**2)
S = clean((A + B) / r)
Delta = clean((A - B)**2 / r**2)

show("P_root", P_root)
show("Q_root", Q_root)
show("N", N)
show("S", S)
show("Delta", Delta)

check("N=P_root*Q_root", N, P_root * Q_root)
check("S=P_root+Q_root", S, P_root + Q_root)
check("Delta=(P_root-Q_root)^2", Delta, (P_root - Q_root)**2)

print()


# ==============================================================================
# [2] ROOT TRANSLATION ORBIT
# ==============================================================================

print("[2] ROOT TRANSLATION ORBIT")
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
        f"L({m})=N+{m}S+{m*m}",
        L(m),
        N + m*S + m**2,
    )

for m in range(5):
    check(
        f"D({m})=S+{2*m+1}",
        D(m),
        S + 2*m + 1,
    )

for m in range(4):
    check(
        f"D({m+1})-D({m})=2",
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
check("F3=D(0)*F2", F3, D(0) * F2)

print()


# ==============================================================================
# [4] ROOT-ORBIT REPRESENTATION
# ==============================================================================

print("[4] ROOT-ORBIT REPRESENTATION")
print("-" * 78)

check(
    "F2=L(0)+L(1)-1-Delta",
    F2,
    L(0) + L(1) - 1 - Delta,
)

check(
    "F3=(L(1)-L(0))*F2",
    F3,
    (L(1) - L(0)) * F2,
)

check(
    "F3/F2=D(0)",
    F3,
    D(0) * F2,
)

print()


# ==============================================================================
# [5] ROOT-ORBIT STEP OPERATOR
# ==============================================================================

print("[5] ROOT-ORBIT STEP OPERATOR")
print("-" * 78)

G2 = clean(D(0))
G3 = clean(D(1))
G4 = clean(D(2))
G5 = clean(D(3))
G6 = clean(D(4))

show("G2", G2)
show("G3", G3)
show("G4", G4)
show("G5", G5)
show("G6", G6)

check("G2=S+1", G2, S + 1)
check("G3=S+3", G3, S + 3)
check("G4=S+5", G4, S + 5)
check("G5=S+7", G5, S + 7)
check("G6=S+9", G6, S + 9)

print()


# ==============================================================================
# [6] FORMAL HIGHER-LAYER PREDICTIONS
# ==============================================================================

print("[6] FORMAL HIGHER-LAYER PREDICTIONS")
print("-" * 78)

F4_expected = clean(G3 * F3)
F5_expected = clean(G4 * F4_expected)
F6_expected = clean(G5 * F5_expected)

show("F4_expected", F4_expected)
show("F5_expected", F5_expected)
show("F6_expected", F6_expected)

check(
    "F4_expected=(S+3)F3",
    F4_expected,
    (S + 3) * F3,
)

check(
    "F5_expected=(S+5)F4_expected",
    F5_expected,
    (S + 5) * F4_expected,
)

check(
    "F6_expected=(S+7)F5_expected",
    F6_expected,
    (S + 7) * F5_expected,
)

print()


# ==============================================================================
# [7] ACTUAL HIGHER-LAYER INPUT
# ==============================================================================

print("[7] ACTUAL HIGHER-LAYER INPUT")
print("-" * 78)

print("""
  F4_actual, F5_actual and F6_actual are placeholders.

  Replace ONLY these assignments with the genuine formulas
  from the original homogeneous-layer construction.

  Do NOT copy the formal predictions.
""")

show("F4_actual", F4_actual)
show("F5_actual", F5_actual)
show("F6_actual", F6_actual)

print()


# ==============================================================================
# [8] DIRECT HIGHER-LAYER STEP RESIDUALS
# ==============================================================================

print("[8] DIRECT HIGHER-LAYER STEP RESIDUALS")
print("-" * 78)

R4 = clean(
    F4_actual - G3 * F3
)

R5 = clean(
    F5_actual - G4 * F4_actual
)

R6 = clean(
    F6_actual - G5 * F5_actual
)

show("R4=F4_actual-D1*F3", R4)
show("R5=F5_actual-D2*F4_actual", R5)
show("R6=F6_actual-D3*F5_actual", R6)

print()


# ==============================================================================
# [9] DIVISION-FREE HIGHER-LAYER DEFECTS
# ==============================================================================

print("[9] DIVISION-FREE HIGHER-LAYER DEFECTS")
print("-" * 78)

E2 = clean(
    F4_actual * F2
    - F3**2
    - 2 * F3 * F2
)

E3 = clean(
    F5_actual * F3
    - F4_actual**2
    - 2 * F4_actual * F3
)

E4 = clean(
    F6_actual * F4_actual
    - F5_actual**2
    - 2 * F5_actual * F4_actual
)

show("E2", E2)
show("E3", E3)
show("E4", E4)

print("""
  Required for continuation of the root-orbit step:

      E2 = 0
      E3 = 0
      E4 = 0.
""")

print()


# ==============================================================================
# [10] RATIO-LADDER EQUIVALENCE
# ==============================================================================

print("[10] RATIO-LADDER EQUIVALENCE")
print("-" * 78)

# E2 = 0 is equivalent to
#
#     F4/F3 - F3/F2 = 2
#
# after division by F2*F3.
#
# We verify the numerator identity symbolically.

ratio_defect_2 = clean(
    F4_actual * F2
    - F3 * (F3 + 2*F2)
)

ratio_defect_3 = clean(
    F5_actual * F3
    - F4_actual * (F4_actual + 2*F3)
)

ratio_defect_4 = clean(
    F6_actual * F4_actual
    - F5_actual * (F5_actual + 2*F4_actual)
)

check(
    "E2 equals ratio-step numerator",
    E2,
    ratio_defect_2,
)

check(
    "E3 equals ratio-step numerator",
    E3,
    ratio_defect_3,
)

check(
    "E4 equals ratio-step numerator",
    E4,
    ratio_defect_4,
)

print()


# ==============================================================================
# [11] ROOT-ORBIT STEP COMPARISON
# ==============================================================================

print("[11] ROOT-ORBIT STEP COMPARISON")
print("-" * 78)

check(
    "D(0)=F3/F2",
    F3,
    D(0) * F2,
)

print("""
  The genuine higher-layer continuation would require:

      F4/F3 = D(1)
      F5/F4 = D(2)
      F6/F5 = D(3)

  namely:

      F4/F3 = S+3
      F5/F4 = S+5
      F6/F5 = S+7.
""")

print()


# ==============================================================================
# [12] KAPPA SPECTRAL POLYNOMIAL
# ==============================================================================

print("[12] KAPPA SPECTRAL POLYNOMIAL")
print("-" * 78)

Qpoly = clean(
    z**2 - S*z + N
)

Q_factored = clean(
    (z - P_root) * (z - Q_root)
)

show("Qpoly(z)", Qpoly)
show("factor form", Q_factored)

check(
    "Qpoly=(z-P_root)(z-Q_root)",
    Qpoly,
    Q_factored,
)

check(
    "discriminant(Qpoly)=Delta",
    sp.discriminant(Qpoly, z),
    Delta,
)

print()


# ==============================================================================
# [13] TRANSLATED KAPPA FACTOR
# ==============================================================================

print("[13] TRANSLATED KAPPA FACTOR")
print("-" * 78)

Q_shift = clean(
    Qpoly.subs(z, z - 1)
)

Q_shift_expected = clean(
    (z - P_root - 1) * (z - Q_root - 1)
)

check(
    "Q(z-1)=(z-P-1)(z-Q-1)",
    Q_shift,
    Q_shift_expected,
)

show("Qpoly(z-1)", Q_shift)

print()


# ==============================================================================
# [14] FOUR-BASE CHARACTERISTIC
# ==============================================================================

print("[14] FOUR-BASE CHARACTERISTIC")
print("-" * 78)

chi = clean(
    Qpoly * Q_shift
)

show("chi(z)", chi)

check(
    "chi=Q(z)Q(z-1)",
    chi,
    Qpoly * Q_shift,
)

print()


# ==============================================================================
# [15] ROOT-ORBIT GENERATING POLYNOMIAL
# ==============================================================================

print("[15] ROOT-ORBIT GENERATING POLYNOMIAL")
print("-" * 78)

check(
    "L(-z)=Qpoly(z)",
    L(-z),
    Qpoly,
)

check(
    "L(h)-L(-h)=2hS",
    L(h) - L(-h),
    2*h*S,
)

check(
    "L(h+2)-2L(h+1)+L(h)=2",
    L(h + 2) - 2*L(h + 1) + L(h),
    2,
)

print()


# ==============================================================================
# [16] KNOWN PAIR AS ROOT-ORBIT DATA
# ==============================================================================

print("[16] KNOWN PAIR AS ROOT-ORBIT DATA")
print("-" * 78)

W = clean(F3 / F2)

check(
    "W=S+1",
    W,
    S + 1,
)

check(
    "W=D(0)",
    W,
    D(0),
)

print("""
  The known homogeneous pair therefore gives:

      W = F3/F2
        = D(0)
        = S+1.

  The root-orbit step sequence is:

      D(0), D(1), D(2), ...

      S+1, S+3, S+5, ...

  The open question is whether the genuine higher layers
  multiply successively by these same steps.
""")

print()


# ==============================================================================
# [17] 2020 COORDINATE CROSSWALK
# ==============================================================================

print("[17] 2020 COORDINATE CROSSWALK")
print("-" * 78)

A0 = clean(A)
B0 = clean(B)

show("A0", A0)
show("B0", B0)
show("A0+B0", A0 + B0)
show("B0-A0", B0 - A0)
show("A0*B0", A0 * B0)

check(
    "S=(A0+B0)/r",
    S,
    (A0 + B0) / r,
)

check(
    "N=A0*B0/r^2",
    N,
    A0 * B0 / r**2,
)

check(
    "Delta=(A0-B0)^2/r^2",
    Delta,
    (A0 - B0)**2 / r**2,
)

print()


# ==============================================================================
# [18] STRUCTURAL SUMMARY
# ==============================================================================

print("[18] STRUCTURAL SUMMARY")
print("-" * 78)

print(
r"""
  EXACTLY ESTABLISHED:

      P=A0/r
      Q=B0/r

      N=A0B0/r^2

      S=(A0+B0)/r

      Delta=(A0-B0)^2/r^2

      L_m=(A0+mr)(B0+mr)/r^2

      L_m=m^2+Sm+N

      D_m=L_(m+1)-L_m=S+2m+1

      F2=L0+L1-1-Delta

      F3=D0 F2

      F3/F2=D0=S+1.

  NOT YET ESTABLISHED:

      F4=D1 F3

      F5=D2 F4

      F6=D3 F5.

  Therefore the next meaningful input is still the genuine
  original F4/F5/F6 construction.
"""
)

print()


# ==============================================================================
# [19] FINAL AUDIT
# ==============================================================================

print("[19] FINAL AUDIT")
print("-" * 78)

print(f"  built-in symbolic failures = {failures}")
print("  genuine F4/F5/F6 supplied = NO")
print("  higher-layer continuation = OPEN")

print()
print("=" * 78)
print("EXPERIMENT 563 FINISHED")
print("=" * 78)
