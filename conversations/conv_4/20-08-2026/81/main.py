#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 559
# ==============================================================================
# ACTUAL HIGHER HOMOGENEOUS LAYERS -> 2020 ROOT-ORBIT COMPARISON
#
# Purpose:
#
#   Experiment 558 established:
#
#       L_m = (A + m*r)(B + m*r) / r^2
#
#       D_m = L_(m+1)-L_m = S + 2m + 1
#
#       F2 = (L0 + L1 - 1 - Delta)
#
#       F3 = D0 * F2
#
# What remains unproved:
#
#       F4 ?= D1 F3
#       F5 ?= D2 F4
#       F6 ?= D3 F5
#
# This experiment deliberately separates:
#
#   1. exact root-orbit identities
#   2. formal predictions
#   3. genuine higher-layer formulas
#
# Replace the placeholder F4_actual/F5_actual/F6_actual expressions
# below with the REAL formulas from the original homogeneous-layer
# construction.
#
# No numerical fitting.
# No factor enumeration.
# No p/q insertion.
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 559 START")
print("=" * 78)
print("ACTUAL HIGHER HOMOGENEOUS LAYERS -> 2020 ROOT-ORBIT COMPARISON")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
z, m, n, h = sp.symbols("z m n h")

# Placeholder symbols for genuine higher layers.
# These must be replaced by actual formulas before a genuine
# higher-layer audit can be claimed.
F4_actual = sp.Symbol("F4_actual")
F5_actual = sp.Symbol("F5_actual")
F6_actual = sp.Symbol("F6_actual")


# ==============================================================================
# AUDIT STATE
# ==============================================================================

failures = 0
validation_pending = False


# ==============================================================================
# EXACT SIMPLIFICATION
# ==============================================================================

def simp(expr):
    return sp.factor(
        sp.cancel(
            sp.together(
                sp.simplify(expr)
            )
        )
    )


def numerator_difference(actual, expected):
    """
    Convert a rational symbolic equality into a polynomial numerator.

    This prevents errors such as:

        PolynomialError: r**(-2) contains an element ...

    """
    diff = simp(actual - expected)
    num, den = sp.fraction(sp.together(diff))
    return sp.factor(num), sp.factor(den)


def cert(name, actual, expected=0, count_failure=True):
    global failures

    num, den = numerator_difference(actual, expected)
    passed = (num == 0)

    print(f"  {name}")
    print("    numerator difference =")
    sp.pprint(num)

    if den != 1 and num != 0:
        print("    denominator =")
        sp.pprint(den)

    print(f"    PASS = {passed}")
    print()

    if not passed and count_failure:
        failures += 1

    return passed


def show(name, expr):
    print(f"  {name} =")
    sp.pprint(simp(expr))
    print()


# ==============================================================================
# [1] ROOT COORDINATES
# ==============================================================================

print("[1] HISTORICAL ROOT COORDINATES")
print("-" * 78)

P = simp(A / r)
Q = simp(B / r)

N = simp(A * B / r**2)
S = simp((A + B) / r)
Delta = simp((A - B)**2 / r**2)

show("P", P)
show("Q", Q)
show("N", N)
show("S", S)
show("Delta", Delta)

cert("N=PQ", N, P * Q)
cert("S=P+Q", S, P + Q)
cert("Delta=(P-Q)^2", Delta, (P - Q)**2)

print()


# ==============================================================================
# [2] ROOT TRANSLATION ORBIT
# ==============================================================================

print("[2] ROOT TRANSLATION ORBIT")
print("-" * 78)


def L(k):
    return simp(
        (A + k*r) * (B + k*r) / r**2
    )


def D(k):
    return simp(
        L(k + 1) - L(k)
    )


for k0 in range(6):
    show(f"L{k0}", L(k0))

for k0 in range(5):
    cert(
        f"D{k0}=S+{2*k0+1}",
        D(k0),
        S + 2*k0 + 1,
    )

cert(
    "L(m)=m^2+S*m+N",
    L(m),
    m**2 + S*m + N,
)

cert(
    "D(m)=S+2m+1",
    D(m),
    S + 2*m + 1,
)

cert(
    "D(m+1)-D(m)=2",
    D(m + 1) - D(m),
    2,
)

print()


# ==============================================================================
# [3] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[3] KNOWN HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = simp(
    6*N - S**2 + S
)

F3 = simp(
    (S + 1) * F2
)

cert(
    "F2=2N+S-Delta",
    F2,
    2*N + S - Delta,
)

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1)*F2,
)

cert(
    "F3=D0*F2",
    F3,
    D(0)*F2,
)

print()


# ==============================================================================
# [4] KNOWN ROOT-ORBIT REPRESENTATION
# ==============================================================================

print("[4] KNOWN ROOT-ORBIT REPRESENTATION")
print("-" * 78)

F2_orbit = simp(
    L(0) + L(1) - 1 - Delta
)

F3_orbit = simp(
    D(0) * F2
)

cert(
    "F2=L0+L1-1-Delta",
    F2_orbit,
    F2,
)

cert(
    "F3=D0*F2",
    F3_orbit,
    F3,
)

print()


# ==============================================================================
# [5] FORMAL HIGHER-LAYER PREDICTIONS
# ==============================================================================

print("[5] FORMAL HIGHER-LAYER PREDICTIONS")
print("-" * 78)

F4_pred = simp(
    D(1) * F3
)

F5_pred = simp(
    D(2) * F4_pred
)

F6_pred = simp(
    D(3) * F5_pred
)

show("F4_pred", F4_pred)
show("F5_pred", F5_pred)
show("F6_pred", F6_pred)

cert(
    "F4_pred=(S+3)F3",
    F4_pred,
    (S + 3)*F3,
)

cert(
    "F5_pred=(S+5)F4_pred",
    F5_pred,
    (S + 5)*F4_pred,
)

cert(
    "F6_pred=(S+7)F5_pred",
    F6_pred,
    (S + 7)*F5_pred,
)

print()


# ==============================================================================
# [6] GENUINE HIGHER-LAYER INPUT STATUS
# ==============================================================================

print("[6] GENUINE HIGHER-LAYER INPUT STATUS")
print("-" * 78)

print(
    r"""
  The following symbols are placeholders:

      F4_actual
      F5_actual
      F6_actual

  They are NOT predictions.

  A genuine validation requires replacing them with the
  actual expressions generated by the original homogeneous
  construction.
"""
)

print("  F4_actual placeholder =", F4_actual)
print("  F5_actual placeholder =", F5_actual)
print("  F6_actual placeholder =", F6_actual)
print()

validation_pending = True


# ==============================================================================
# [7] FIRST ACTUAL-LAYER TEST
# ==============================================================================

print("[7] FIRST ACTUAL-LAYER TEST")
print("-" * 78)

E4 = simp(
    F4_actual - D(1)*F3
)

print("  Required identity:")
print()
print("      F4_actual - D1*F3 = 0")
print()

show("E4_actual", E4)

print("  STATUS = PENDING ACTUAL F4 FORMULA")
print()


# ==============================================================================
# [8] SECOND ACTUAL-LAYER TEST
# ==============================================================================

print("[8] SECOND ACTUAL-LAYER TEST")
print("-" * 78)

E5 = simp(
    F5_actual - D(2)*F4_actual
)

print("  Required identity:")
print()
print("      F5_actual - D2*F4_actual = 0")
print()

show("E5_actual", E5)

print("  STATUS = PENDING ACTUAL F5 FORMULA")
print()


# ==============================================================================
# [9] THIRD ACTUAL-LAYER TEST
# ==============================================================================

print("[9] THIRD ACTUAL-LAYER TEST")
print("-" * 78)

E6 = simp(
    F6_actual - D(3)*F5_actual
)

print("  Required identity:")
print()
print("      F6_actual - D3*F5_actual = 0")
print()

show("E6_actual", E6)

print("  STATUS = PENDING ACTUAL F6 FORMULA")
print()


# ==============================================================================
# [10] DIVISION-FREE DEFECTS
# ==============================================================================

print("[10] DIVISION-FREE HIGHER-LAYER DEFECTS")
print("-" * 78)

print(
    r"""
  The ratio tests can be converted into division-free forms.

      E2 = F4*F2 - F3^2 - 2F3F2

      E3 = F5*F3 - F4^2 - 2F4F3

      E4 = F6*F4 - F5^2 - 2F5F4
"""
)

E2_actual = simp(
    F4_actual * F2
    - F3**2
    - 2*F3*F2
)

E3_actual = simp(
    F5_actual * F3
    - F4_actual**2
    - 2*F4_actual*F3
)

E4_actual_divfree = simp(
    F6_actual * F4_actual
    - F5_actual**2
    - 2*F5_actual*F4_actual
)

show("E2_divfree", E2_actual)
show("E3_divfree", E3_actual)
show("E4_divfree", E4_actual_divfree)

print()


# ==============================================================================
# [11] FORMAL DEFECT CONTROL
# ==============================================================================

print("[11] FORMAL DEFECT CONTROL")
print("-" * 78)

formal_E2 = simp(
    F4_pred * F2
    - F3**2
    - 2*F3*F2
)

formal_E3 = simp(
    F5_pred * F3
    - F4_pred**2
    - 2*F4_pred*F3
)

formal_E4 = simp(
    F6_pred * F4_pred
    - F5_pred**2
    - 2*F5_pred*F4_pred
)

cert(
    "formal E2=0",
    formal_E2,
    0,
)

cert(
    "formal E3=0",
    formal_E3,
    0,
)

cert(
    "formal E4=0",
    formal_E4,
    0,
)

print()


# ==============================================================================
# [12] ROOT-ORBIT STEP IDENTIFICATION
# ==============================================================================

print("[12] ROOT-ORBIT STEP IDENTIFICATION")
print("-" * 78)

G2 = D(0)
G3 = D(1)
G4 = D(2)
G5 = D(3)

cert("G2=S+1", G2, S + 1)
cert("G3=S+3", G3, S + 3)
cert("G4=S+5", G4, S + 5)
cert("G5=S+7", G5, S + 7)

cert(
    "G3-G2=2",
    G3 - G2,
    2,
)

cert(
    "G4-G3=2",
    G4 - G3,
    2,
)

cert(
    "G5-G4=2",
    G5 - G4,
    2,
)

print()


# ==============================================================================
# [13] INDEX-LAW FORM
# ==============================================================================

print("[13] INDEX-LAW FORM")
print("-" * 78)

R = simp(
    S + 2*n - 3
)

R_shift = simp(
    R.subs(n, n + h)
)

cert(
    "R_(n+h)=R_n+2h",
    R_shift,
    R + 2*h,
)

cert(
    "R_n-2n=S-3",
    R - 2*n,
    S - 3,
)

print()


# ==============================================================================
# [14] SPECTRAL CONNECTION
# ==============================================================================

print("[14] SPECTRAL CONNECTION")
print("-" * 78)

Q_poly = simp(
    z**2 - S*z + N
)

Q_factor = simp(
    (z - A/r) * (z - B/r)
)

Q_shift = simp(
    Q_poly.subs(z, z - 1)
)

chi = simp(
    Q_poly * Q_shift
)

cert(
    "Q=(z-A/r)(z-B/r)",
    Q_poly,
    Q_factor,
)

cert(
    "disc(Q)=Delta",
    sp.discriminant(Q_poly, z),
    Delta,
)

cert(
    "Q(z-1)",
    Q_shift,
    (z - (A+r)/r) * (z - (B+r)/r),
)

cert(
    "chi=Q(z)Q(z-1)",
    chi,
    Q_factor *
    ((z - (A+r)/r) * (z - (B+r)/r)),
)

print()


# ==============================================================================
# [15] CENTRAL ROOT-ORBIT CHANNEL
# ==============================================================================

print("[15] CENTRAL ROOT-ORBIT CHANNEL")
print("-" * 78)

cert(
    "L(h)-L(-h)=2hS",
    L(h) - L(-h),
    2*h*S,
)

cert(
    "Delta^2 L=2",
    L(m + 2) - 2*L(m + 1) + L(m),
    2,
)

print()


# ==============================================================================
# [16] ACTUAL FORMULA INSERTION GUIDE
# ==============================================================================

print("[16] ACTUAL FORMULA INSERTION GUIDE")
print("-" * 78)

print(
    r"""
  Replace ONLY these three assignments:

      F4_actual = ...
      F5_actual = ...
      F6_actual = ...

  with the genuine formulas from the original homogeneous-layer
  construction.

  Then the script automatically checks:

      F4_actual = D1 F3

      F5_actual = D2 F4_actual

      F6_actual = D3 F5_actual

  and the equivalent division-free identities:

      F4 F2 - F3^2 - 2F3F2 = 0

      F5 F3 - F4^2 - 2F4F3 = 0

      F6 F4 - F5^2 - 2F5F4 = 0.
"""
)

print()


# ==============================================================================
# [17] FINAL STRUCTURAL CERTIFICATE
# ==============================================================================

print("[17] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print(
    r"""
  EXACT:

      P=A/r
      Q=B/r

      N=AB/r^2
      S=(A+B)/r
      Delta=(A-B)^2/r^2

      L_m=(A+mr)(B+mr)/r^2

      D_m=L_(m+1)-L_m
         = S+2m+1

      F2=L0+L1-1-Delta

      F3=D0 F2.

  FORMAL HYPOTHESIS:

      F_(n+1)=D_(n-2)F_n.

  FIRST UNKNOWN CONSEQUENCES:

      F4=D1F3
      F5=D2F4
      F6=D3F5.

  THIS SCRIPT DOES NOT CLAIM THESE ARE TRUE.

  They become established only when the genuine original
  homogeneous-layer formulas are inserted and the exact
  residuals vanish.
"""
)

print()


# ==============================================================================
# FINAL STATUS
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 559 FINISHED")
print("=" * 78)
print()

# Only exact built-in identities are counted in the audit.
# The actual higher-layer validation is explicitly pending.
print(f"BUILT-IN SYMBOLIC FAILURES = {failures}")
print(f"ACTUAL F4/F5/F6 VALIDATION = PENDING")
print(f"OVERALL EXPERIMENT STATUS = {'PENDING ACTUAL HIGHER LAYERS'}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print(
    r"""
Use the genuine F4, F5 and F6 formulas from the original
homogeneous-layer construction.

Do not replace them with the formal predictions.

The decisive identities are:

    F4 - D1*F3 = 0

    F5 - D2*F4 = 0

    F6 - D3*F5 = 0

or equivalently:

    F4F2 - F3^2 - 2F3F2 = 0

    F5F3 - F4^2 - 2F4F3 = 0

    F6F4 - F5^2 - 2F5F4 = 0.

A zero result would establish that the original homogeneous
layer index continues the 2020 root-coordinate translation
operator.
"""
)
