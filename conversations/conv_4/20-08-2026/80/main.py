#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 558
# ==============================================================================
# 2020 ROOT-ORBIT -> HOMOGENEOUS LAYER-STEP OPERATOR
#
# Corrected version of Experiment 557.
#
# Important:
#   The formal hierarchy begins at F2.
#   Therefore recurrence checks involving F_(n-1) must begin at n=3.
#
# Proven:
#   L_m = N + m*S + m^2
#   D_m = L_(m+1)-L_m = S + 2m + 1
#   F2 = L0 + L1 - 1 - Delta
#   F3 = D0 * F2
#
# Hypothesis:
#   F_(n+1) = D_(n-2) F_n
#
# This script does NOT claim genuine F4/F5/F6 formulas.
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 558 START")
print("=" * 78)
print("2020 ROOT-ORBIT -> HOMOGENEOUS LAYER-STEP OPERATOR")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
m, n, h, k, z = sp.symbols("m n h k z")
alpha, beta = sp.symbols("alpha beta")

failures = 0


# ==============================================================================
# SAFE EXACT SIMPLIFICATION
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
# [1] ROOT COORDINATES
# ==============================================================================

print("[1] ROOT COORDINATES")
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


# ==============================================================================
# [2] ROOT TRANSLATION ORBIT
# ==============================================================================

print("[2] ROOT TRANSLATION ORBIT")
print("-" * 78)


def L(j):
    return simp(
        (A + j*r) * (B + j*r) / r**2
    )


def D(j):
    return simp(
        L(j + 1) - L(j)
    )


for j in range(4):
    show(f"L{j}", L(j))

cert("L0=N", L(0), N)
cert("L1=N+S+1", L(1), N + S + 1)
cert("L2=N+2S+4", L(2), N + 2*S + 4)
cert("L3=N+3S+9", L(3), N + 3*S + 9)

cert("D0=S+1", D(0), S + 1)
cert("D1=S+3", D(1), S + 3)
cert("D2=S+5", D(2), S + 5)
cert("D3=S+7", D(3), S + 7)


# ==============================================================================
# [3] GENERAL ORBIT LAW
# ==============================================================================

print("[3] GENERAL ORBIT LAW")
print("-" * 78)

L_general = simp(L(m))
D_general = simp(D(m))

show("L(m)", L_general)
show("D(m)", D_general)

cert(
    "L(m)=m^2+S*m+N",
    L_general,
    m**2 + S*m + N,
)

cert(
    "D(m)=S+2m+1",
    D_general,
    S + 2*m + 1,
)

cert(
    "D(m+1)-D(m)=2",
    D(m + 1) - D(m),
    2,
)


# ==============================================================================
# [4] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[4] KNOWN HOMOGENEOUS LAYERS")
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


# ==============================================================================
# [5] ROOT-ORBIT REPRESENTATION OF F2
# ==============================================================================

print("[5] ROOT-ORBIT REPRESENTATION OF F2")
print("-" * 78)

F2_orbit = simp(
    L(0) + L(1) - 1 - Delta
)

cert(
    "F2=L0+L1-1-Delta",
    F2_orbit,
    F2,
)


# ==============================================================================
# [6] ROOT-ORBIT REPRESENTATION OF F3
# ==============================================================================

print("[6] ROOT-ORBIT REPRESENTATION OF F3")
print("-" * 78)

F3_orbit = simp(
    D(0) * F2
)

cert(
    "F3=D0*F2",
    F3_orbit,
    F3,
)

cert(
    "F3/F2=D0",
    sp.cancel(F3 / F2),
    D(0),
)


# ==============================================================================
# [7] AFFINE STEP CALIBRATION
# ==============================================================================

print("[7] AFFINE STEP CALIBRATION")
print("-" * 78)

G = alpha*n + beta

eq_2 = simp(
    G.subs(n, 2) - (S + 1)
)

eq_3 = simp(
    G.subs(n, 3) - (S + 3)
)

solution = sp.solve(
    [eq_2, eq_3],
    [alpha, beta],
    dict=True,
)

print("  affine solutions =")
print(solution)
print()

if solution:
    G_calibrated = simp(
        G.subs(solution[0])
    )

    show("G_n", G_calibrated)

    cert(
        "G_n=S+2n-3",
        G_calibrated,
        S + 2*n - 3,
    )

    cert(
        "G_n=D_(n-2)",
        G_calibrated,
        D(n - 2),
    )
else:
    failures += 1
    print("  FAILED TO SOLVE AFFINE CALIBRATION")
    print()


# ==============================================================================
# [8] ROOT-ORBIT STEP INVARIANT
# ==============================================================================

print("[8] ROOT-ORBIT STEP INVARIANT")
print("-" * 78)

G_n = simp(
    S + 2*n - 3
)

cert(
    "G_n-2n=S-3",
    G_n - 2*n,
    S - 3,
)


# ==============================================================================
# [9] ROOT TRANSLATION OF THE STEP
# ==============================================================================

print("[9] ROOT TRANSLATION OF THE STEP")
print("-" * 78)

S_h = simp(
    S + 2*h
)

G_translated = simp(
    S_h + 2*n - 3
)

cert(
    "S_h=S+2h",
    S_h,
    S + 2*h,
)

cert(
    "G_translated=G_n+2h",
    G_translated,
    G_n + 2*h,
)


# ==============================================================================
# [10] FORMAL HOMOGENEOUS HIERARCHY
# ==============================================================================

print("[10] FORMAL HOMOGENEOUS HIERARCHY")
print("-" * 78)

F2_seed = sp.Symbol("F2_seed", nonzero=True)

formal_layers = {
    2: F2_seed,
}

# F_(n+1) = D_(n-2) F_n
for layer in range(2, 8):
    step = simp(
        D(layer - 2)
    )

    formal_layers[layer + 1] = simp(
        step * formal_layers[layer]
    )

for layer in range(2, 8):
    show(
        f"formal F{layer}",
        formal_layers[layer],
    )


# ==============================================================================
# [11] FORMAL RATIO LADDER
# ==============================================================================

print("[11] FORMAL RATIO LADDER")
print("-" * 78)

for layer in range(2, 7):

    ratio = simp(
        sp.cancel(
            formal_layers[layer + 1]
            / formal_layers[layer]
        )
    )

    expected = simp(
        S + 2*layer - 3
    )

    cert(
        f"F{layer+1}/F{layer}=S+2n-3",
        ratio,
        expected,
    )


# ==============================================================================
# [12] FORMAL DIVISION-FREE RECURRENCE
# ==============================================================================

print("[12] FORMAL DIVISION-FREE RECURRENCE")
print("-" * 78)

def defect(F_prev, F_curr, F_next):
    return simp(
        F_next * F_prev
        - F_curr**2
        - 2*F_curr*F_prev
    )


# IMPORTANT:
# The formal sequence starts at F2.
#
# Therefore the first valid recurrence comparison is:
#
#     F4*F2 - F3^2 - 2F3F2
#
# which corresponds to the transition centered on F2/F3/F4.
#
# We deliberately DO NOT attempt to access F1.

E2 = defect(
    formal_layers[2],
    formal_layers[3],
    formal_layers[4],
)

E3 = defect(
    formal_layers[3],
    formal_layers[4],
    formal_layers[5],
)

E4 = defect(
    formal_layers[4],
    formal_layers[5],
    formal_layers[6],
)

E5 = defect(
    formal_layers[5],
    formal_layers[6],
    formal_layers[7],
)

cert("E2=0", E2, 0)
cert("E3=0", E3, 0)
cert("E4=0", E4, 0)
cert("E5=0", E5, 0)


# ==============================================================================
# [13] SAME RECURRENCE IN STEP-FACTOR FORM
# ==============================================================================

print("[13] STEP-FACTOR RECURRENCE")
print("-" * 78)

for layer in range(3, 7):

    lhs = formal_layers[layer + 1]

    rhs = simp(
        D(layer - 2) * formal_layers[layer]
    )

    cert(
        f"F{layer+1}=D{layer-2}*F{layer}",
        lhs,
        rhs,
    )


# ==============================================================================
# [14] OBSERVABLE FORM OF THE FIRST UNKNOWN LAYER
# ==============================================================================

print("[14] FIRST UNKNOWN LAYER")
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
    "F5_pred=(S+5)F4",
    F5_pred,
    (S + 5)*F4_pred,
)

cert(
    "F6_pred=(S+7)F5",
    F6_pred,
    (S + 7)*F5_pred,
)


# ==============================================================================
# [15] ROOT-ORBIT COMPARISON
# ==============================================================================

print("[15] ROOT-ORBIT / HOMOGENEOUS STEP COMPARISON")
print("-" * 78)

for layer in range(2, 6):

    step = simp(
        D(layer - 2)
    )

    expected = simp(
        S + 2*layer - 3
    )

    cert(
        f"D{layer-2}=S+{2*layer-3}",
        step,
        expected,
    )


# ==============================================================================
# [16] FORMAL INDEX TRANSLATION
# ==============================================================================

print("[16] FORMAL INDEX TRANSLATION")
print("-" * 78)

R_n = simp(
    S + 2*n - 3
)

R_nh = simp(
    R_n.subs(n, n + h)
)

cert(
    "R_(n+h)=R_n+2h",
    R_nh,
    R_n + 2*h,
)

cert(
    "R_n-2n is invariant",
    R_n - 2*n,
    S - 3,
)


# ==============================================================================
# [17] KAPPA GENERATING QUADRATIC
# ==============================================================================

print("[17] KAPPA GENERATING QUADRATIC")
print("-" * 78)

Q_poly = simp(
    z**2 - S*z + N
)

Q_factor = simp(
    (z - A/r) * (z - B/r)
)

Q_orbit = simp(
    L(-z)
)

show("Q(z)", Q_poly)

cert(
    "Q=(z-A/r)(z-B/r)",
    Q_poly,
    Q_factor,
)

cert(
    "L(-z)=Q(z)",
    Q_orbit,
    Q_poly,
)

cert(
    "disc(Q)=Delta",
    sp.discriminant(Q_poly, z),
    Delta,
)


# ==============================================================================
# [18] FOUR-BASE COMPOSITION
# ==============================================================================

print("[18] FOUR-BASE COMPOSITION")
print("-" * 78)

Q_shift = simp(
    Q_poly.subs(z, z - 1)
)

chi = simp(
    Q_poly * Q_shift
)

chi_factor = simp(
    ((z - A/r) * (z - B/r))
    *
    ((z - (A+r)/r) * (z - (B+r)/r))
)

cert(
    "Q(z-1) is translated-root factor",
    Q_shift,
    ((z - (A+r)/r) * (z - (B+r)/r)),
)

cert(
    "four-base composition",
    chi,
    chi_factor,
)


# ==============================================================================
# [19] CENTRAL ROOT-ORBIT DIFFERENCE
# ==============================================================================

print("[19] CENTRAL ROOT-ORBIT DIFFERENCE")
print("-" * 78)

central = simp(
    L(h) - L(-h)
)

cert(
    "L(h)-L(-h)=2hS",
    central,
    2*h*S,
)


# ==============================================================================
# [20] CURVATURE OF THE ROOT ORBIT
# ==============================================================================

print("[20] ROOT-ORBIT CURVATURE")
print("-" * 78)

curvature = simp(
    L(m + 2)
    - 2*L(m + 1)
    + L(m)
)

cert(
    "Delta^2 L=2",
    curvature,
    2,
)


# ==============================================================================
# [21] PROOF / HYPOTHESIS BOUNDARY
# ==============================================================================

print("[21] PROOF / HYPOTHESIS BOUNDARY")
print("-" * 78)

print(
    r"""
  EXACTLY ESTABLISHED:

      A0 = historical first root coordinate
      B0 = historical second root coordinate

      L_m = (A0+mr)(B0+mr)/r^2

      L_m = N + mS + m^2

      D_m = L_(m+1)-L_m
          = S+2m+1

      F2 = L0+L1-1-Delta

      F3 = D0*F2

      F3/F2 = D0 = S+1.

  FORMALLY DERIVED FROM THE AFFINE STEP HYPOTHESIS:

      F_(n+1) = D_(n-2) F_n

      F_(n+1)/F_n = S+2n-3

      F4=(S+3)F3

      F5=(S+5)F4

      F6=(S+7)F5

      E_n=0.

  NOT ESTABLISHED:

      That the genuine original 2020 homogeneous construction
      actually produces those F4/F5/F6 expressions.

  Therefore the next experiment must use genuine higher-layer
  formulas, rather than these formal predictions.
"""
)


# ==============================================================================
# [22] FINAL STRUCTURAL CERTIFICATE
# ==============================================================================

print("[22] FINAL STRUCTURAL CERTIFICATE")
print("-" * 78)

print(
    r"""
  The strongest proven connection currently available is:

      2020 ROOT COORDINATES
          |
          v
      A0, B0
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
          |
          v
      D_m=S+2m+1

  The known homogeneous layers satisfy:

      F2=L0+L1-1-Delta

      F3=D0 F2.

  Thus the first homogeneous layer transition is exactly
  the first root-orbit step.

  The unresolved question remains:

      Does the ORIGINAL hierarchy continue as

          F4=D1 F3
          F5=D2 F4
          F6=D3 F5 ?

  This experiment does not assume that it does.
"""
)


# ==============================================================================
# FINAL AUDIT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 558 FINISHED")
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
The next experiment should insert the ACTUAL F4/F5/F6 formulas
from the original homogeneous-layer construction.

Then test:

    F4 - D1*F3

    F5 - D2*F4

    F6 - D3*F5.

Equivalently:

    F4/F3 = S+3
    F5/F4 = S+5
    F6/F5 = S+7.

Do not replace genuine layers by the formal predictions.

The purpose of the next step is to determine whether the
root-coordinate translation operator actually generates the
higher homogeneous-layer hierarchy.
"""
)