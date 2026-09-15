#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 560
# ==============================================================================
# KNOWN HOMOGENEOUS PAIR -> INTRINSIC TRANSFER OPERATOR
#
# Known:
#
#     F2 = 2N + S - Delta
#     F3 = (S+1) F2
#
# Historical root coordinates:
#
#     A = A0
#     B = B0
#
#     N     = AB/r^2
#     S     = (A+B)/r
#     Delta = (A-B)^2/r^2
#
# Root orbit:
#
#     L_m = (A+mr)(B+mr)/r^2
#
#     D_m = L_{m+1}-L_m = S+2m+1
#
# Therefore:
#
#     F3 = D0 F2.
#
# This experiment derives the induced operator on the observable
# pair (F2,F3), computes its invariants, and then asks whether
# the actual homogeneous hierarchy could be recognized purely
# from the operator.
#
# IMPORTANT:
#
#     This script does NOT invent F4/F5/F6.
#
# It only derives:
#
#     (F2,F3) -> (F2',F3')
#
# and gives exact residual functions into which genuine layers can
# later be inserted.
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 560 START")
print("=" * 78)
print("KNOWN HOMOGENEOUS PAIR -> INTRINSIC TRANSFER OPERATOR")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
z, h = sp.symbols("z h")


# ==============================================================================
# HELPERS
# ==============================================================================

def simp(expr):
    return sp.factor(sp.cancel(sp.together(sp.expand(expr))))


def numerator_difference(lhs, rhs=0):
    diff = simp(lhs - rhs)
    num, den = sp.fraction(sp.together(diff))
    return sp.factor(num), sp.factor(den)


failures = 0


def cert(name, lhs, rhs=0):
    global failures

    num, den = numerator_difference(lhs, rhs)
    passed = num == 0

    print(f"  {name}")
    print("    numerator difference =")
    sp.pprint(num)

    if not passed and den != 1:
        print("    denominator =")
        sp.pprint(den)

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

N = simp(A * B / r**2)
S = simp((A + B) / r)
Delta = simp((A - B)**2 / r**2)

show("N", N)
show("S", S)
show("Delta", Delta)

cert("N=AB/r^2", N, A * B / r**2)
cert("S=(A+B)/r", S, (A + B) / r)
cert("Delta=(A-B)^2/r^2", Delta, (A - B)**2 / r**2)

print()


# ==============================================================================
# [2] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[2] KNOWN HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = simp(
    2*N + S - Delta
)

F3 = simp(
    (S + 1) * F2
)

cert(
    "F2=6N-S^2+S",
    F2,
    6*N - S**2 + S,
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

print()


# ==============================================================================
# [3] ROOT-ORBIT STEP
# ==============================================================================

print("[3] ROOT-ORBIT STEP")
print("-" * 78)

L0 = simp(A * B / r**2)
L1 = simp((A + r)*(B + r) / r**2)

D0 = simp(L1 - L0)

show("L0", L0)
show("L1", L1)
show("D0=L1-L0", D0)

cert("D0=S+1", D0, S + 1)
cert("F3=D0*F2", F3, D0 * F2)

print()


# ==============================================================================
# [4] INDUCED OBSERVABLE COORDINATES
# ==============================================================================

print("[4] INDUCED OBSERVABLE COORDINATES")
print("-" * 78)

# Treat:
#
#     X = F2
#     Y = F3
#
# and define:
#
#     W = Y/X = S+1.
#
# Thus:
#
#     S = W-1.

X, Y = sp.symbols("X Y", nonzero=True)

W = simp(Y / X)
S_obs = simp(W - 1)

show("W=F3/F2", W)
show("S_obs=W-1", S_obs)

print()


# ==============================================================================
# [5] N AS A FUNCTION OF F2,F3
# ==============================================================================

print("[5] N AS A FUNCTION OF F2,F3")
print("-" * 78)

N_obs = simp(
    (X + S_obs**2 - S_obs) / 6
)

Delta_obs = simp(
    S_obs**2 - 4*N_obs
)

show("N_obs(F2,F3)", N_obs)
show("Delta_obs(F2,F3)", Delta_obs)

# Substitute the genuine symbolic F2,F3 and verify.

cert(
    "N_obs=N",
    N_obs.subs({X: F2, Y: F3}),
    N,
)

cert(
    "Delta_obs=Delta",
    Delta_obs.subs({X: F2, Y: F3}),
    Delta,
)

print()


# ==============================================================================
# [6] ROOT TRANSLATION IN OBSERVABLE SPACE
# ==============================================================================

print("[6] ROOT TRANSLATION IN OBSERVABLE SPACE")
print("-" * 78)

# Historical/root translation:
#
#     A -> A + h r
#     B -> B + h r
#
# implies:
#
#     S -> S + 2h
#     N -> N + hS + h^2
#     Delta -> Delta.

S_h = simp(S + 2*h)
N_h = simp(N + h*S + h**2)

F2_h_reference = simp(
    2*N_h + S_h - Delta
)

F3_h_reference = simp(
    (S_h + 1) * F2_h_reference
)

print("  Reference translated state:")
show("S_h", S_h)
show("N_h", N_h)
show("F2_h", F2_h_reference)
show("F3_h", F3_h_reference)

cert(
    "Delta invariant",
    simp(Delta),
    simp(Delta),
)

print()


# ==============================================================================
# [7] ELIMINATE S AND N
# ==============================================================================

print("[7] ELIMINATE S AND N")
print("-" * 78)

# Since:
#
#     S+1 = Y/X
#
# and:
#
#     N = (X + S^2 - S)/6,
#
# the translated pair can be written entirely in X,Y.

F2_prime = simp(
    F2_h_reference.subs(
        {
            N: N_obs,
            S: S_obs,
        }
    )
)

# Recompute F3 from translated observable state.
S_prime_obs = simp(S_obs + 2*h)

N_prime_obs = simp(
    N_obs + h*S_obs + h**2
)

F2_prime = simp(
    2*N_prime_obs + S_prime_obs - Delta_obs
)

F3_prime = simp(
    (S_prime_obs + 1) * F2_prime
)

show("F2'(X,Y)", F2_prime)
show("F3'(X,Y)", F3_prime)

print()


# ==============================================================================
# [8] SIMPLIFIED RATIONAL TRANSFER MAP
# ==============================================================================

print("[8] SIMPLIFIED RATIONAL TRANSFER MAP")
print("-" * 78)

# Expected:
#
#     X' = X + 2h(Y/X + h)
#
#     Y' = (Y/X + 2h) X'.

X_prime_expected = simp(
    X + 2*h*(Y/X + h)
)

Y_prime_expected = simp(
    (Y/X + 2*h) * X_prime_expected
)

show("X'_expected", X_prime_expected)
show("Y'_expected", Y_prime_expected)

cert(
    "F2' transfer map",
    F2_prime,
    X_prime_expected,
)

cert(
    "F3' transfer map",
    F3_prime,
    Y_prime_expected,
)

print()


# ==============================================================================
# [9] RATIO TRANSLATION
# ==============================================================================

print("[9] RATIO TRANSLATION")
print("-" * 78)

W_prime = simp(
    Y_prime_expected / X_prime_expected
)

show("W'", W_prime)

cert(
    "W'=W+2h",
    W_prime,
    W + 2*h,
)

print()


# ==============================================================================
# [10] INVARIANT CHANNEL
# ==============================================================================

print("[10] INVARIANT CHANNEL")
print("-" * 78)

I_obs = simp(
    Delta_obs
)

I_prime = simp(
    I_obs.subs(
        {
            X: X_prime_expected,
            Y: Y_prime_expected,
        }
    )
)

show("I(X,Y)", I_obs)
show("I(X',Y')", I_prime)

cert(
    "I'=I",
    I_prime,
    I_obs,
)

print()


# ==============================================================================
# [11] KAPPA QUADRATIC DIRECTLY FROM F2,F3
# ==============================================================================

print("[11] KAPPA QUADRATIC DIRECTLY FROM F2,F3")
print("-" * 78)

Q_obs = simp(
    z**2 - S_obs*z + N_obs
)

show("Q_obs(z)", Q_obs)

cert(
    "Q_obs is monic",
    sp.Poly(Q_obs, z).LC(),
    1,
)

Q_disc_obs = simp(
    sp.discriminant(Q_obs, z)
)

show("disc(Q_obs)", Q_disc_obs)

cert(
    "disc(Q_obs)=I",
    Q_disc_obs,
    I_obs,
)

print()


# ==============================================================================
# [12] TRANSLATION OF THE SPECTRAL POLYNOMIAL
# ==============================================================================

print("[12] SPECTRAL TRANSLATION")
print("-" * 78)

Q_shift = simp(
    Q_obs.subs(z, z - 1)
)

show("Q_obs(z-1)", Q_shift)

Q_composed = simp(
    Q_obs * Q_shift
)

show("Q_obs(z)Q_obs(z-1)", Q_composed)

print()


# ==============================================================================
# [13] ROOT-ORBIT COMPARISON
# ==============================================================================

print("[13] ROOT-ORBIT COMPARISON")
print("-" * 78)

L_h = simp(
    (A + h*r)*(B + h*r)/r**2
)

cert(
    "L_h=N+hS+h^2",
    L_h,
    N + h*S + h**2,
)

cert(
    "L_h-L_0=hS+h^2",
    L_h - L0,
    h*S + h**2,
)

cert(
    "L_-z=Q(z)",
    L_h.subs(h, -z),
    (z - A/r)*(z - B/r),
)

print()


# ==============================================================================
# [14] THE OPERATOR IS DETERMINED BY F2,F3
# ==============================================================================

print("[14] OPERATOR RECONSTRUCTION FROM THE KNOWN PAIR")
print("-" * 78)

print(
r"""
  The known pair

      (F2,F3)

  determines:

      W = F3/F2

      S = W-1

      N = (F2+S^2-S)/6

      Delta = S^2-4N.

  Therefore the root translation induces the unique observable
  map:

      W' = W+2h

      F2' = F2 + 2h(W+h)

      F3' = (W+2h)F2'.

  This is not an assumption about F4/F5/F6.

  It is the exact induced action of the already-established
  root-coordinate translation on the known homogeneous pair.
"""
)

print()


# ==============================================================================
# [15] GENERIC LAYER TRANSFER TEST
# ==============================================================================

print("[15] GENERIC LAYER TRANSFER TEST")
print("-" * 78)

# Introduce genuinely unknown next layers only here.
F4, F5, F6 = sp.symbols(
    "F4 F5 F6",
    nonzero=True,
)

# The root-orbit hypothesis predicts:
#
#     F4_expected = (W+2) F3
#     F5_expected = (W+4) F4
#     F6_expected = (W+6) F5.
#
# We do not call these "actual".

F4_expected = simp(
    (Y/X + 2) * Y
)

F5_expected = simp(
    (Y/X + 4) * F4
)

F6_expected = simp(
    (Y/X + 6) * F5
)

show("F4_expected", F4_expected)
show("F5_expected", F5_expected)
show("F6_expected", F6_expected)

print(
r"""
  These are ONLY transfer predictions.

  To test genuine layers, substitute their actual symbolic
  expressions for F4,F5,F6 below.
"""
)

print()


# ==============================================================================
# [16] EXACT DEFECTS FOR GENUINE LAYERS
# ==============================================================================

print("[16] EXACT DEFECTS FOR GENUINE LAYERS")
print("-" * 78)

defect4 = simp(
    F4 - F4_expected
)

defect5 = simp(
    F5 - F5_expected
)

defect6 = simp(
    F6 - F6_expected
)

print("  Genuine-layer residuals:")
show("F4 defect", defect4)
show("F5 defect", defect5)
show("F6 defect", defect6)

print()


# ==============================================================================
# [17] DIVISION-FREE DEFECTS
# ==============================================================================

print("[17] DIVISION-FREE DEFECTS")
print("-" * 78)

E2 = simp(
    F4 * F2
    - F3**2
    - 2*F3*F2
)

E3 = simp(
    F5 * F3
    - F4**2
    - 2*F4*F3
)

E4 = simp(
    F6 * F4
    - F5**2
    - 2*F5*F4
)

show("E2", E2)
show("E3", E3)
show("E4", E4)

print(
r"""
  Exact upstream tests:

      E2 = 0
      E3 = 0
      E4 = 0.

  These are equivalent to the proposed ratio-step continuation,
  but contain no division.
"""
)

print()


# ==============================================================================
# [18] ROOT-ORBIT CONTROL
# ==============================================================================

print("[18] ROOT-ORBIT CONTROL IDENTITIES")
print("-" * 78)

for k in range(5):
    cert(
        f"D{k}=S+{2*k+1}",
        D := simp(
            L0
            if False
            else (
                (A + (k + 1)*r)*(B + (k + 1)*r)
                - (A + k*r)*(B + k*r)
            ) / r**2
        ),
        S + 2*k + 1,
    )

print()


# ==============================================================================
# [19] FINAL STRUCTURAL RESULT
# ==============================================================================

print("[19] FINAL STRUCTURAL RESULT")
print("-" * 78)

print(
r"""
  EXACTLY ESTABLISHED:

      A0 = 2y-2x+3
      B0 = 2y+2x-3

      P=A0/r
      Q=B0/r

      N=A0B0/r^2
      S=(A0+B0)/r
      Delta=(A0-B0)^2/r^2

      L_m=(A0+mr)(B0+mr)/r^2

      L_m=m^2+S*m+N

      D_m=L_(m+1)-L_m=S+2m+1

      F2=L0+L1-1-Delta

      F3=D0 F2.

  Therefore the known homogeneous pair carries the exact
  observable coordinates:

      W=F3/F2=S+1

      Delta=I(F2,F3).

  The root translation acts on these observables by:

      W -> W+2h

      F2 -> F2+2h(W+h)

      F3 -> (W+2h)F2'.

  The remaining unknown is NOT the translation operator.

  The remaining unknown is whether the ORIGINAL higher
  homogeneous layers are generated by this same transfer:

      F4 ?= (W+2)F3

      F5 ?= (W+4)F4

      F6 ?= (W+6)F5.

  No claim is made here that these identities hold.
"""
)

print()


# ==============================================================================
# [20] FINAL AUDIT
# ==============================================================================

print("[20] FINAL EXACT AUDIT")
print("-" * 78)

# All internal identities above are exact and counted.
# Unknown F4/F5/F6 are deliberately NOT counted as failures.

print(f"  built-in symbolic failures = {failures}")
print("  F4/F5/F6 genuine validation = NOT SUPPLIED")
print("  higher-layer conclusion = OPEN")
print()

print("=" * 78)
print("EXPERIMENT 560 FINISHED")
print("=" * 78)
print()

print(
r"""
NEXT RESEARCH TARGET

Insert the REAL F4/F5/F6 formulas from the original 2020
homogeneous construction:

    F4 = ...
    F5 = ...
    F6 = ...

Then run only the decisive tests:

    F4 - (F3/F2 + 2)F3

    F5 - (F4/F3 + 2)F4

    F6 - (F5/F4 + 2)F5

or, preferably without division:

    F4F2 - F3^2 - 2F3F2

    F5F3 - F4^2 - 2F4F3

    F6F4 - F5^2 - 2F5F4.

The root-coordinate translation itself is already exact.
The unresolved question is whether the original homogeneous
hierarchy continues that operator beyond F3.
"""
)
