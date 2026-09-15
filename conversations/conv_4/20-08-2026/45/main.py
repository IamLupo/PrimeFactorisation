#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 521 START")
print("=" * 78)
print("TWO-LAYER TRANSLATION OPERATOR / UPSTREAM QUADRATIC SEARCH")
print()

# ============================================================================
# SYMBOLS
# ============================================================================

U, V, r, z, t = sp.symbols("U V r z t", nonzero=True)
N, S, Delta = sp.symbols("N S Delta")

# Historical normalized coordinates
P = (U - V) / r
Q = (U + V) / r

N_hist = sp.expand(P * Q)
S_hist = sp.expand(P + Q)
D_hist = sp.expand((P - Q) ** 2)

# Historical conic layers
K0 = sp.expand(U**2 - V**2)
K1 = sp.expand((U + r)**2 - V**2)

# ============================================================================
# CERTIFICATE
# ============================================================================

def cert(label, expr, expected=0):
    diff = sp.factor(sp.expand(expr - expected))
    ok = diff == 0
    print(f"  {label}")
    print(f"    difference = {diff}")
    print(f"    PASS = {ok}")
    return ok

failures = 0

def check(label, expr, expected=0):
    global failures
    ok = cert(label, expr, expected)
    if not ok:
        failures += 1
    return ok


# ============================================================================
# [1] HISTORICAL TWO-LAYER SYSTEM
# ============================================================================

print("[1] HISTORICAL TWO-LAYER SYSTEM")
print("-" * 78)

print("  K0 =", K0)
print("  K1 =", K1)

check("K0/r^2 = N", K0 / r**2, N_hist)
check(
    "(K1-K0)/r^2 = S+1",
    (K1 - K0) / r**2,
    S_hist + 1,
)

print()


# ============================================================================
# [2] FIRST DIFFERENCE AS THE MISSING LINEAR COEFFICIENT
# ============================================================================

print("[2] FIRST DIFFERENCE AS THE MISSING LINEAR COEFFICIENT")
print("-" * 78)

dK = sp.expand((K1 - K0) / r**2)

candidate_S = sp.expand(dK - 1)
candidate_a = sp.expand(-candidate_S)

print("  dK =", sp.factor(dK))
print("  S_candidate =", sp.factor(candidate_S))
print("  a_candidate =", sp.factor(candidate_a))

check("S_candidate = S", candidate_S, S_hist)
check("a_candidate = -S", candidate_a, -S_hist)

print()


# ============================================================================
# [3] RECONSTRUCT THE MISSING QUADRATIC
# ============================================================================

print("[3] RECONSTRUCT THE MISSING QUADRATIC")
print("-" * 78)

a = candidate_a
b = sp.expand(K0 / r**2)

Q0 = sp.expand(z**2 + a*z + b)
Q0_expected = sp.expand((z - P) * (z - Q))

print("  Q0(z) =", sp.factor(Q0))

check(
    "Q0 = (z-P)(z-Q)",
    Q0,
    Q0_expected,
)

print()


# ============================================================================
# [4] UNIVERSAL TRANSLATION
# ============================================================================

print("[4] UNIVERSAL TRANSLATION")
print("-" * 78)

Q1 = sp.expand(Q0.subs(z, z - 1))
Q1_expected = sp.expand((z - P - 1) * (z - Q - 1))

check(
    "Q1(z)=Q0(z-1)",
    Q1,
    Q1_expected,
)

chi = sp.expand(Q0 * Q1)
chi_expected = sp.expand(
    (z-P)*(z-Q)*(z-P-1)*(z-Q-1)
)

check(
    "Q0*Q1 = four-base characteristic polynomial",
    chi,
    chi_expected,
)

print()


# ============================================================================
# [5] COEFFICIENT-LEVEL TRANSLATION
# ============================================================================

print("[5] COEFFICIENT-LEVEL TRANSLATION")
print("-" * 78)

a0 = sp.expand(a)
b0 = sp.expand(b)

a1 = sp.expand(a0 - 2)
b1 = sp.expand(b0 - a0 + 1)

print("  a0 =", sp.factor(a0))
print("  b0 =", sp.factor(b0))
print("  a1 =", sp.factor(a1))
print("  b1 =", sp.factor(b1))

check("a1 = -(S+2)", a1, -(S_hist + 2))
check("b1 = N+S+1", b1, N_hist + S_hist + 1)

print()


# ============================================================================
# [6] SECOND DIFFERENCE / TRANSLATION CLOSURE
# ============================================================================

print("[6] SECOND DIFFERENCE / TRANSLATION CLOSURE")
print("-" * 78)

# Translation of the conic invariant:
#
# K_m = (U+m r)^2 - V^2
#
# The first difference is affine in m.
# The second difference should be exactly 2 r^2.

m = sp.symbols("m", integer=True)

Km = sp.expand((U + m*r)**2 - V**2)
Km1 = sp.expand((U + (m+1)*r)**2 - V**2)
Km2 = sp.expand((U + (m+2)*r)**2 - V**2)

d1 = sp.expand(Km1 - Km)
d2 = sp.expand(Km2 - 2*Km1 + Km)

print("  K_m =", Km)
print("  ΔK_m =", sp.factor(d1))
print("  Δ²K_m =", sp.factor(d2))

check(
    "first difference",
    d1,
    r*(2*U + (2*m+1)*r),
)

check(
    "second difference = 2r^2",
    d2,
    2*r**2,
)

print()


# ============================================================================
# [7] GENERAL TRANSLATED QUADRATIC FAMILY
# ============================================================================

print("[7] GENERAL TRANSLATED QUADRATIC FAMILY")
print("-" * 78)

# For
#
#   K_m/r^2 = N + m(S+1) + m^2
#
# the successive values form a quadratic sequence in m.

Km_norm = sp.expand(Km / r**2)
Km_norm_expected = sp.expand(
    N_hist + m*(S_hist + 1) + m**2
)

print("  K_m/r^2 =", sp.factor(Km_norm))

check(
    "normalized conic translation law",
    Km_norm,
    Km_norm_expected,
)

print()


# ============================================================================
# [8] INTERPRETATION AS A GENERATING QUADRATIC
# ============================================================================

print("[8] INTERPRETATION AS A GENERATING QUADRATIC")
print("-" * 78)

# Observe:
#
#   K_m/r^2 = Q(m) + m
#
# where Q(m) = m^2 + S*m + N.
#
# Therefore
#
#   Q(m) = K_m/r^2 - m.
#
# This is particularly important because it says the entire
# missing quadratic can be reconstructed from the translation
# orbit of the historical conic.

Qm_from_K = sp.expand(Km_norm - m)
Qm_target = sp.expand(m**2 + S_hist*m + N_hist)

print("  Q(m) from conic =", sp.factor(Qm_from_K))
print("  target Q(m)     =", sp.factor(Qm_target))

check(
    "Q(m)=K_m/r^2-m",
    Qm_from_K,
    Qm_target,
)

print()


# ============================================================================
# [9] ROOT INTERPRETATION
# ============================================================================

print("[9] ROOT INTERPRETATION")
print("-" * 78)

Qm_roots = sp.solve(sp.Eq(Qm_target, 0), m)

print("  Q(m) =", sp.factor(Qm_target))
print("  roots =", Qm_roots)

check(
    "Q(-P)=0",
    Qm_target.subs(m, -P),
    0,
)

check(
    "Q(-Q)=0",
    Qm_target.subs(m, -Q),
    0,
)

print()


# ============================================================================
# [10] DISCRIMINANT OF THE TRANSLATION ORBIT
# ============================================================================

print("[10] DISCRIMINANT OF THE TRANSLATION ORBIT")
print("-" * 78)

disc = sp.expand(
    sp.discriminant(Qm_target, m)
)

print("  discriminant =", sp.factor(disc))

check(
    "translation-orbit discriminant = Delta",
    disc,
    D_hist,
)

print()


# ============================================================================
# [11] CAN A SINGLE TRANSLATION VALUE CARRY S?
# ============================================================================

print("[11] SINGLE TRANSLATION VALUE")
print("-" * 78)

# From
#
#   K_m/r^2 = N + m(S+1) + m^2
#
# we have, for any nonzero integer m,
#
#   S = (K_m/r^2 - N - m^2)/m - 1.
#
# This shows that ONE nonzero translation layer is enough
# once N is already known.

S_from_Km = sp.expand(
    (Km_norm - N_hist - m**2) / m - 1
)

print("  S_from_Km =", sp.factor(S_from_Km))

check(
    "S recovered from one nonzero translated conic layer",
    S_from_Km,
    S_hist,
)

print()


# ============================================================================
# [12] TWO-LAYER MINIMALITY
# ============================================================================

print("[12] TWO-LAYER MINIMALITY")
print("-" * 78)

# With only K0:
#
#   K0/r^2 = N
#
# no linear-in-U information exists.
#
# K1 adds exactly the odd/U-sensitive component.

K0_even = sp.expand((K0 + K0.subs(U, -U)) / 2)
K1K0 = sp.expand(K1 - K0)

print("  even part K0 =", sp.factor(K0_even))
print("  K1-K0        =", sp.factor(K1K0))

check(
    "K0 is U-even",
    K0.subs(U, -U),
    K0,
)

check(
    "K1-K0 contains U",
    sp.expand(sp.diff(K1K0, U)),
    2*r,
)

print()


# ============================================================================
# [13] ABSTRACT UPSTREAM SIGNATURE
# ============================================================================

print("[13] ABSTRACT UPSTREAM SIGNATURE")
print("-" * 78)

print("""
  The historical system has the exact orbit

      L_m = K_m/r^2
          = N + m(S+1) + m^2.

  Therefore:

      ΔL_m = S+1 + 2m+1
      Δ²L_m = 2.

  Equivalently,

      Q(m) = L_m - m
           = m^2 + S m + N.

  Hence:

      coefficient(m^2) = 1
      coefficient(m)   = S
      coefficient(1)   = N.

  The missing KAPPA quadratic is literally the
  generating polynomial of the translated conic orbit.
""")

print()


# ============================================================================
# [14] SYMBOLIC SEARCH FOR EQUIVALENT ORBITS
# ============================================================================

print("[14] EQUIVALENT ORBIT FORMS")
print("-" * 78)

# Three equivalent representations:
#
#   (A) K_m
#   (B) normalized K_m
#   (C) Q(m) = K_m/r^2 - m

orbit_A = sp.expand(Km)
orbit_B = sp.expand(Km / r**2)
orbit_C = sp.expand(Km / r**2 - m)

check("A -> B normalization", orbit_A / r**2, orbit_B)
check(
    "B -> C subtraction",
    orbit_B - m,
    orbit_C,
)

check(
    "C -> generating quadratic",
    orbit_C,
    m**2 + S_hist*m + N_hist,
)

print()


# ============================================================================
# [15] FINAL BRIDGE AUDIT
# ============================================================================

print("[15] FINAL BRIDGE AUDIT")
print("-" * 78)

bridge_checks = [
    ("K0 -> N", K0/r**2, N_hist),
    ("K1-K0 -> S+1", (K1-K0)/r**2, S_hist+1),
    ("Q0 linear coefficient -> -S", a0, -S_hist),
    ("Q0 constant coefficient -> N", b0, N_hist),
    ("Q1 translation", Q1, Q0.subs(z, z-1)),
    ("K_m orbit", orbit_B, N_hist + m*(S_hist+1) + m**2),
    ("Q(m) orbit", orbit_C, m**2 + S_hist*m + N_hist),
    ("Q discriminant", disc, D_hist),
]

for label, lhs, rhs in bridge_checks:
    check(label, lhs, rhs)

print()


# ============================================================================
# CONCLUSION
# ============================================================================

print("=" * 78)
print("EXPERIMENT 521 FINISHED")
print("=" * 78)

print("""
CONCLUSION
------------------------------------------------------------------------------

The historical conic does more than provide two isolated objects.

It possesses an exact translation orbit

    K_m = (U+m r)^2 - V^2.

After normalization,

    L_m = K_m/r^2
        = N + m(S+1) + m^2.

Subtracting the known translation term m gives

    Q(m) = L_m - m
         = m^2 + S m + N.

This is precisely the KAPPA quadratic

    Q(m) = m^2 + S m + N.

Its roots are

    m = -P, -Q,

and its discriminant is

    S^2 - 4N = Delta.

This produces a much tighter structural bridge than merely
saying that K0 and K1 contain N and S.

The historical conic translation orbit itself is a generating
mechanism for the missing quadratic.

The new upstream target is therefore:

    N-only homogeneous construction
        -> a translated quadratic orbit L_m
        -> Q(m)=L_m-m
        -> coefficient of m
        -> S
        -> Delta
        -> p,q.

The critical question for the next experiment is now:

    Does the existing homogeneous-layer construction contain
    an intrinsic translation parameter whose successive layers
    obey a quadratic law

        L_m = N + m(S+1) + m^2

    or an equivalent affine-translated form?

That is the next place to attack. It is closer to the original
homogeneous construction and does not require reconstructing
F0,...,F7.
""")

if failures:
    print(f"SYMBOLIC FAILURES = {failures}")
    raise SystemExit(1)

print("OVERALL EXACT SYMBOLIC AUDIT = True")
