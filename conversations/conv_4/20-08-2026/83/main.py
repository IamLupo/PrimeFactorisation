#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 561
# ==============================================================================
# CORRECTED F2/F3 OBSERVABLE OPERATOR
# -> INVARIANT CHANNEL
# -> ROOT-ORBIT TRANSFER
# -> HIGHER-LAYER DEFECT TEST
#
# This experiment fixes the substitution bug from Experiment 560.
#
# IMPORTANT:
#     substitutions involving X' and Y' are performed simultaneously.
#
# No F4/F5/F6 formulas are invented.
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 561 START")
print("=" * 78)
print("CORRECTED F2/F3 OBSERVABLE OPERATOR")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r = sp.symbols("A B r", nonzero=True)
X, Y = sp.symbols("X Y", nonzero=True)
h, z = sp.symbols("h z")


# ==============================================================================
# HELPERS
# ==============================================================================

failures = 0


def clean(expr):
    return sp.factor(sp.cancel(sp.together(expr)))


def numdiff(lhs, rhs=0):
    diff = clean(lhs - rhs)
    num, den = sp.fraction(sp.together(diff))
    return sp.factor(num), sp.factor(den)


def check(name, lhs, rhs=0):
    global failures

    num, den = numdiff(lhs, rhs)
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
    sp.pprint(clean(expr))
    print()


# ==============================================================================
# [1] HISTORICAL ROOT COORDINATES
# ==============================================================================

print("[1] HISTORICAL ROOT COORDINATES")
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
# [2] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[2] KNOWN HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = clean(6*N - S**2 + S)
F3 = clean((S + 1) * F2)

check("F2=2N+S-Delta", F2, 2*N + S - Delta)
check("F3=(S+1)F2", F3, (S + 1) * F2)

print()


# ==============================================================================
# [3] ROOT-ORBIT CONTROL
# ==============================================================================

print("[3] ROOT-ORBIT CONTROL")
print("-" * 78)

def L(m):
    return clean(
        (A + m*r) * (B + m*r) / r**2
    )


def D(m):
    return clean(
        L(m + 1) - L(m)
    )


for k in range(5):
    check(
        f"D{k}=S+{2*k+1}",
        D(k),
        S + 2*k + 1,
    )

check(
    "L(m)=m^2+S*m+N",
    L(sp.Symbol("m")),
    sp.Symbol("m")**2 + S*sp.Symbol("m") + N,
)

print()


# ==============================================================================
# [4] OBSERVABLE COORDINATES
# ==============================================================================

print("[4] OBSERVABLE COORDINATES")
print("-" * 78)

W = clean(Y / X)
S_obs = clean(W - 1)

N_obs = clean(
    (X + S_obs**2 - S_obs) / 6
)

Delta_obs = clean(
    S_obs**2 - 4*N_obs
)

show("W=F3/F2", W)
show("S_obs=W-1", S_obs)
show("N_obs", N_obs)
show("Delta_obs", Delta_obs)

check(
    "N_obs(F2,F3)=N",
    N_obs.subs(
        {X: F2, Y: F3},
        simultaneous=True,
    ),
    N,
)

check(
    "Delta_obs(F2,F3)=Delta",
    Delta_obs.subs(
        {X: F2, Y: F3},
        simultaneous=True,
    ),
    Delta,
)

print()


# ==============================================================================
# [5] OBSERVABLE TRANSLATION OPERATOR
# ==============================================================================

print("[5] OBSERVABLE TRANSLATION OPERATOR")
print("-" * 78)

Xprime = clean(
    X + 2*h*(Y/X + h)
)

Yprime = clean(
    (Y/X + 2*h) * Xprime
)

show("X'", Xprime)
show("Y'", Yprime)

Wprime = clean(
    Yprime / Xprime
)

show("W'", Wprime)

check(
    "W'=W+2h",
    Wprime,
    W + 2*h,
)

print()


# ==============================================================================
# [6] CORRECT INVARIANT SUBSTITUTION
# ==============================================================================

print("[6] CORRECT INVARIANT SUBSTITUTION")
print("-" * 78)

# IMPORTANT:
# simultaneous=True prevents SymPy from replacing X inside Yprime
# after X has already been replaced by Xprime.

Delta_prime = clean(
    Delta_obs.subs(
        {
            X: Xprime,
            Y: Yprime,
        },
        simultaneous=True,
    )
)

show("Delta(X,Y)", Delta_obs)
show("Delta(X',Y')", Delta_prime)

check(
    "Delta'=Delta",
    Delta_prime,
    Delta_obs,
)

print()


# ==============================================================================
# [7] N TRANSLATION
# ==============================================================================

print("[7] N TRANSLATION")
print("-" * 78)

N_prime = clean(
    N_obs.subs(
        {
            X: Xprime,
            Y: Yprime,
        },
        simultaneous=True,
    )
)

show("N(X',Y')", N_prime)

check(
    "N'=N+hS+h^2",
    N_prime,
    N_obs + h*S_obs + h**2,
)

print()


# ==============================================================================
# [8] S TRANSLATION
# ==============================================================================

print("[8] S TRANSLATION")
print("-" * 78)

S_prime = clean(
    S_obs.subs(
        {
            X: Xprime,
            Y: Yprime,
        },
        simultaneous=True,
    )
)

show("S(X',Y')", S_prime)

check(
    "S'=S+2h",
    S_prime,
    S_obs + 2*h,
)

print()


# ==============================================================================
# [9] KAPPA QUADRATIC IN OBSERVABLE COORDINATES
# ==============================================================================

print("[9] KAPPA QUADRATIC")
print("-" * 78)

Qobs = clean(
    z**2 - S_obs*z + N_obs
)

disc_Qobs = clean(
    sp.discriminant(Qobs, z)
)

show("Q_obs(z)", Qobs)
show("disc(Q_obs)", disc_Qobs)

check(
    "Q_obs monic",
    sp.Poly(Qobs, z).LC(),
    1,
)

check(
    "disc(Q_obs)=Delta",
    disc_Qobs,
    Delta_obs,
)

print()


# ==============================================================================
# [10] SPECTRAL TRANSLATION
# ==============================================================================

print("[10] SPECTRAL TRANSLATION")
print("-" * 78)

Qshift = clean(
    Qobs.subs(z, z - 1)
)

chi = clean(
    Qobs * Qshift
)

show("Q(z-1)", Qshift)
show("chi(z)=Q(z)Q(z-1)", chi)

print()


# ==============================================================================
# [11] ROOT-ORBIT REPRESENTATION
# ==============================================================================

print("[11] ROOT-ORBIT REPRESENTATION")
print("-" * 78)

Lh = clean(
    (A + h*r)*(B + h*r)/r**2
)

check(
    "L(h)=N+hS+h^2",
    Lh,
    N + h*S + h**2,
)

check(
    "L(-z)=Q(z)",
    Lh.subs(h, -z),
    Qobs.subs(
        {
            X: F2,
            Y: F3,
        },
        simultaneous=True,
    ),
)

print()


# ==============================================================================
# [12] F2 TRANSLATION
# ==============================================================================

print("[12] F2 TRANSLATION")
print("-" * 78)

F2prime_from_NS = clean(
    2*N_prime + S_prime - Delta_prime
)

F3prime_from_NS = clean(
    (S_prime + 1) * F2prime_from_NS
)

Xprime_expected = clean(Xprime)
Yprime_expected = clean(Yprime)

check(
    "F2'=observable X'",
    F2prime_from_NS.subs(
        {
            X: F2,
            Y: F3,
        },
        simultaneous=True,
    ),
    Xprime_expected.subs(
        {
            X: F2,
            Y: F3,
        },
        simultaneous=True,
    ),
)

print()


# ==============================================================================
# [13] OPERATOR INVARIANT IN OBSERVABLE SPACE
# ==============================================================================

print("[13] OPERATOR INVARIANT")
print("-" * 78)

I = clean(Delta_obs)

Iprime = clean(
    I.subs(
        {
            X: Xprime,
            Y: Yprime,
        },
        simultaneous=True,
    )
)

check(
    "I'=I",
    Iprime,
    I,
)

print()


# ==============================================================================
# [14] MOVING / INVARIANT COORDINATE SPLIT
# ==============================================================================

print("[14] MOVING / INVARIANT COORDINATE SPLIT")
print("-" * 78)

show("moving W", W)
show("invariant Delta", I)

check(
    "moving coordinate W'=W+2h",
    Wprime,
    W + 2*h,
)

check(
    "invariant coordinate Delta'=Delta",
    Iprime,
    I,
)

print()


# ==============================================================================
# [15] HIGHER-LAYER TRANSFER PREDICTIONS
# ==============================================================================

print("[15] HIGHER-LAYER TRANSFER PREDICTIONS")
print("-" * 78)

F4, F5, F6 = sp.symbols(
    "F4 F5 F6",
    nonzero=True,
)

# From:
#
#     R_n = S + 2n - 3
#
# the first unknown ratios would be:
#
#     F4/F3 = S+3
#     F5/F4 = S+5
#     F6/F5 = S+7.
#
# These remain TEST TARGETS.

F4_pred = clean(
    (W + 2) * Y
)

F5_pred = clean(
    (W + 4) * F4
)

F6_pred = clean(
    (W + 6) * F5
)

show("F4_pred", F4_pred)
show("F5_pred", F5_pred)
show("F6_pred", F6_pred)

print()


# ==============================================================================
# [16] GENUINE HIGHER-LAYER RESIDUALS
# ==============================================================================

print("[16] GENUINE HIGHER-LAYER RESIDUALS")
print("-" * 78)

R4 = clean(F4 - F4_pred)
R5 = clean(F5 - F5_pred)
R6 = clean(F6 - F6_pred)

show("R4=F4-F4_pred", R4)
show("R5=F5-F5_pred", R5)
show("R6=F6-F6_pred", R6)

print(
"""
  These are intentionally symbolic.

  Do NOT substitute the formal predictions for F4/F5/F6.

  Insert only the genuine homogeneous-layer formulas here.
"""
)
print()


# ==============================================================================
# [17] DIVISION-FREE TESTS
# ==============================================================================

print("[17] DIVISION-FREE HIGHER-LAYER TESTS")
print("-" * 78)

E2 = clean(
    F4*F2 - F3**2 - 2*F3*F2
)

E3 = clean(
    F5*F3 - F4**2 - 2*F4*F3
)

E4 = clean(
    F6*F4 - F5**2 - 2*F5*F4
)

show("E2", E2)
show("E3", E3)
show("E4", E4)

print()


# ==============================================================================
# [18] OPERATOR INTERPRETATION
# ==============================================================================

print("[18] OPERATOR INTERPRETATION")
print("-" * 78)

print(
r"""
  The exact known pair is:

      F2
      F3

  with:

      F3/F2 = S+1.

  The historical root translation induces:

      W = F3/F2

      W -> W+2h

      Delta -> Delta

      N -> N+hS+h^2.

  Therefore the observable pair has the same structure as the
  2020 root coordinates:

      moving coordinate:
          W

      invariant coordinate:
          Delta.

  The higher-layer question is now isolated completely:

      Does the genuine homogeneous sequence continue with

          F4/F3 = W+2

          F5/F4 = W+4

          F6/F5 = W+6

      at the original state?

  Equivalently:

      E2=0
      E3=0
      E4=0.
"""
)

print()


# ==============================================================================
# [19] FINAL AUDIT
# ==============================================================================

print("[19] FINAL EXACT AUDIT")
print("-" * 78)

print(f"  built-in symbolic failures = {failures}")
print("  genuine F4/F5/F6 supplied = NO")
print("  higher-layer continuation = OPEN")
print()

print("=" * 78)
print("EXPERIMENT 561 FINISHED")
print("=" * 78)
