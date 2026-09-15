#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 554
# ==============================================================================
# 2020 ROOT-COORDINATE HOMOGENEOUS ALGEBRA
#
# Denominator-safe version.
#
# Primitive historical coordinates:
#
#     A = A0 = 2y - 2x + 3
#     B = B0 = 2y + 2x - 3
#
# Normalization:
#
#     P = A/r
#     Q = B/r
#
# Therefore:
#
#     N     = AB/r^2
#     S     = (A+B)/r
#     Delta = (A-B)^2/r^2
#     M1    = (A+r)(B+r)/r^2
#
# Known layers:
#
#     F2 = 6N - S^2 + S
#        = 2N + S - Delta
#
#     F3 = (S+1)F2
#
# Main target:
#
#     L_m = (A+mr)(B+mr)/r^2
#
# and determine whether the homogeneous layers arise naturally
# from the translated A,B orbit.
#
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 554 START")
print("=" * 78)
print("2020 ROOT-COORDINATE HOMOGENEOUS ALGEBRA")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r, h, k, z, m = sp.symbols(
    "A B r h k z m",
    nonzero=True,
)

failures = 0


# ==============================================================================
# HELPERS
# ==============================================================================

def simp(expr):
    """Exact simplification while keeping rational expressions exact."""
    return sp.factor(
        sp.cancel(
            sp.together(
                sp.expand(
                    sp.simplify(expr)
                )
            )
        )
    )


def numerator(expr):
    """Return numerator after exact rational combination."""
    return sp.expand(sp.together(expr).as_numer_denom()[0])


def denominator(expr):
    """Return denominator after exact rational combination."""
    return sp.expand(sp.together(expr).as_numer_denom()[1])


def clear_denominator(expr):
    """
    Convert a rational symbolic identity into a polynomial numerator.

    expr == 0  <=>  numerator(together(expr)) == 0
    provided the denominator is nonzero.
    """
    return sp.factor(numerator(expr))


def poly_coeff_equations(expr, *gens):
    """
    Safely obtain coefficient equations from a rational expression.

    The original failure happened because Poly() was given expressions
    containing r**(-1), r**(-2), etc.

    We first combine fractions and extract the numerator.
    """
    num = numerator(expr)
    poly = sp.Poly(num, *gens)
    return poly, poly.coeffs()


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


# ==============================================================================
# [2] BASIC ROOT ALGEBRA
# ==============================================================================

print("[2] BASIC ROOT ALGEBRA")
print("-" * 78)

cert("N = AB/r^2", N, A * B / r**2)
cert("S = (A+B)/r", S, (A + B) / r)
cert("Delta = (A-B)^2/r^2", Delta, (A - B)**2 / r**2)
cert("M1 = (A+r)(B+r)/r^2", M1, (A + r) * (B + r) / r**2)


# ==============================================================================
# [3] KNOWN HOMOGENEOUS LAYERS
# ==============================================================================

print("[3] KNOWN HOMOGENEOUS LAYERS")
print("-" * 78)

F2 = simp(6 * N - S**2 + S)
F3 = simp((S + 1) * F2)

show("F2", F2)
show("F3", F3)

cert("F2 = 2N+S-Delta", F2, 2 * N + S - Delta)
cert("F3 = (S+1)F2", F3, (S + 1) * F2)


# ==============================================================================
# [4] TRANSLATED ROOTS
# ==============================================================================

print("[4] TRANSLATED ROOTS")
print("-" * 78)

A_h = A + h * r
B_h = B + h * r

show("A_h", A_h)
show("B_h", B_h)

cert("A_h = A+hr", A_h, A + h * r)
cert("B_h = B+hr", B_h, B + h * r)


# ==============================================================================
# [5] TRANSLATED PRODUCT ORBIT
# ==============================================================================

print("[5] TRANSLATED PRODUCT ORBIT")
print("-" * 78)

K_h = simp(A_h * B_h)

K_expected = simp(
    A * B
    + h * r * (A + B)
    + h**2 * r**2
)

show("K_h", K_h)
show("K_expected", K_expected)

cert("K_h expansion", K_h, K_expected)

L_h = simp(K_h / r**2)

cert(
    "L_h = N+hS+h^2",
    L_h,
    N + h * S + h**2,
)


# ==============================================================================
# [6] TRANSLATED SUM
# ==============================================================================

print("[6] TRANSLATED SUM")
print("-" * 78)

C_h = simp(A_h + B_h)
C_norm = simp(C_h / r)

cert("C_h = A+B+2hr", C_h, A + B + 2 * h * r)
cert("normalized sum = S+2h", C_norm, S + 2 * h)


# ==============================================================================
# [7] GAP INVARIANCE
# ==============================================================================

print("[7] GAP INVARIANCE")
print("-" * 78)

G = simp(A - B)
G_h = simp(A_h - B_h)

cert("gap invariant", G_h, G)
cert("Delta invariant", G_h**2 / r**2, Delta)


# ==============================================================================
# [8] FIRST TRANSLATED PRODUCTS
# ==============================================================================

print("[8] FIRST TRANSLATED PRODUCTS")
print("-" * 78)

L0 = simp(A * B / r**2)
L1 = simp((A + r) * (B + r) / r**2)
L2 = simp((A + 2 * r) * (B + 2 * r) / r**2)

show("L0", L0)
show("L1", L1)
show("L2", L2)

cert("L0=N", L0, N)
cert("L1=M1", L1, M1)
cert("L1=N+S+1", L1, N + S + 1)
cert("L2=N+2S+4", L2, N + 2 * S + 4)


# ==============================================================================
# [9] TRANSLATED PRODUCT CURVATURE
# ==============================================================================

print("[9] TRANSLATED PRODUCT CURVATURE")
print("-" * 78)

curvature = simp(L2 - 2 * L1 + L0)

show("Delta^2 L", curvature)

cert("second difference = 2", curvature, 2)


# ==============================================================================
# [10] F2 FROM L0, L1 AND GAP
# ==============================================================================

print("[10] F2 FROM L0, L1 AND GAP")
print("-" * 78)

candidate_F2 = simp(L0 + L1 - 1 - Delta)

show("L0+L1-1-Delta", candidate_F2)

cert(
    "F2 = L0+L1-1-Delta",
    F2,
    candidate_F2,
)


# ==============================================================================
# [11] F2 PURE ROOT-ORBIT FORM
# ==============================================================================

print("[11] F2 ROOT-ORBIT REPRESENTATION")
print("-" * 78)

candidate_F2_2 = simp(
    (
        A * B
        + (A + r) * (B + r)
        - r**2
        - (A - B)**2
    ) / r**2
)

show("candidate", candidate_F2_2)

cert(
    "F2 root-orbit representation",
    candidate_F2_2,
    F2,
)


# ==============================================================================
# [12] TRANSLATED F2
# ==============================================================================

print("[12] TRANSLATED F2")
print("-" * 78)

P_h = simp(A_h / r)
Q_h = simp(B_h / r)

N_h = simp(P_h * Q_h)
S_h = simp(P_h + Q_h)
Delta_h = simp((P_h - Q_h)**2)

F2_h = simp(
    6 * N_h
    - S_h**2
    + S_h
)

show("F2_h", F2_h)
show("F2_h-F2", F2_h - F2)

cert(
    "F2_h-F2 = 2h(S+h+1)",
    F2_h - F2,
    2 * h * (S + h + 1),
)


# ==============================================================================
# [13] F2 ORBIT COEFFICIENTS
# ==============================================================================

print("[13] F2 ORBIT COEFFICIENTS")
print("-" * 78)

F2_poly = sp.Poly(
    sp.expand(
        numerator(F2_h)
    ),
    h,
)

# Since the expression has an r^2 denominator, it is safer to obtain
# coefficients from the original rational expression using .coeff()
# after together(). The polynomial itself is recovered from numerator,
# then we normalize by the common denominator.

F2_h_num, F2_h_den = sp.together(F2_h).as_numer_denom()

F2_h_poly = sp.Poly(
    sp.expand(F2_h_num),
    h,
)

c2_num = simp(F2_h_poly.coeff_monomial(h**2))
c1_num = simp(F2_h_poly.coeff_monomial(h))
c0_num = simp(F2_h_poly.coeff_monomial(h**0))

c2 = simp(c2_num / F2_h_den)
c1 = simp(c1_num / F2_h_den)
c0 = simp(c0_num / F2_h_den)

show("coefficient h^2", c2)
show("coefficient h", c1)
show("constant", c0)

cert("h^2 coefficient", c2, 2)
cert("h coefficient", c1, 2 * (S + 1))
cert("constant", c0, F2)


# ==============================================================================
# [14] F2 ORBIT DISCRIMINANT
# ==============================================================================

print("[14] F2 ORBIT DISCRIMINANT")
print("-" * 78)

D_F2 = simp(
    c1**2
    - 4 * c2 * c0
)

show("disc_h(F2_h)", D_F2)

cert(
    "disc(F2_h)=12Delta+4",
    D_F2,
    12 * Delta + 4,
)


# ==============================================================================
# [15] F3/F2 ROOT RATIO
# ==============================================================================

print("[15] F3/F2 ROOT RATIO")
print("-" * 78)

ratio = simp(F3 / F2)

show("F3/F2", ratio)

cert(
    "F3/F2=(A+B+r)/r",
    ratio,
    (A + B + r) / r,
)

cert(
    "F3/F2-1=(A+B)/r",
    ratio - 1,
    S,
)


# ==============================================================================
# [16] F3 FROM TRANSLATED SUM
# ==============================================================================

print("[16] F3 TRANSLATED-SUM STRUCTURE")
print("-" * 78)

half_shift = sp.Rational(1, 2)

C_half_norm = simp(
    (A + half_shift * r + B + half_shift * r) / r
)

candidate_F3 = simp(
    C_half_norm * F2
)

show("normalized half-shift sum", C_half_norm)
show("F3 candidate", candidate_F3)

cert(
    "F3=(S+1)F2",
    candidate_F3,
    F3,
)


# ==============================================================================
# [17] PRODUCT FIRST DIFFERENCE
# ==============================================================================

print("[17] PRODUCT FIRST DIFFERENCE")
print("-" * 78)

dL = simp(L1 - L0)

show("L1-L0", dL)

cert(
    "L1-L0=S+1",
    dL,
    S + 1,
)


# ==============================================================================
# [18] ROOT-ORBIT GENERATING POLYNOMIAL
# ==============================================================================

print("[18] ROOT-ORBIT GENERATING POLYNOMIAL")
print("-" * 78)

Qz = simp(L_h.subs(h, -z))
Q_expected = simp(z**2 - S * z + N)

show("Q(z)=L(-z)", Qz)
show("z^2-Sz+N", Q_expected)

cert(
    "L(-z)=Q(z)",
    Qz,
    Q_expected,
)


# ==============================================================================
# [19] ROOT FACTORIZATION
# ==============================================================================

print("[19] ROOT FACTORIZATION")
print("-" * 78)

Q_factored = simp(
    (z - A / r) * (z - B / r)
)

show("(z-A/r)(z-B/r)", Q_factored)

cert(
    "Q=(z-A/r)(z-B/r)",
    Qz,
    Q_factored,
)


# ==============================================================================
# [20] TRANSLATED KAPPA FACTOR
# ==============================================================================

print("[20] TRANSLATED KAPPA FACTOR")
print("-" * 78)

Q_shift = simp(Qz.subs(z, z - 1))

Q_shift_expected = simp(
    (z - (A + r) / r)
    * (z - (B + r) / r)
)

cert(
    "Q(z-1) translated roots",
    Q_shift,
    Q_shift_expected,
)


# ==============================================================================
# [21] FOUR-BASE CHARACTERISTIC
# ==============================================================================

print("[21] FOUR-BASE CHARACTERISTIC")
print("-" * 78)

chi = simp(Qz * Q_shift)

show("chi(z)", chi)

cert(
    "four-base composition",
    chi,
    Qz * Q_shift,
)


# ==============================================================================
# [22] TRANSLATED F3
# ==============================================================================

print("[22] TRANSLATED F3")
print("-" * 78)

F3_h = simp((S_h + 1) * F2_h)

show("F3_h", F3_h)

cert(
    "F3_h=(S_h+1)F2_h",
    F3_h,
    (S + 2 * h + 1) * F2_h,
)


# ==============================================================================
# [23] RATIO TRANSLATION
# ==============================================================================

print("[23] RATIO TRANSLATION")
print("-" * 78)

ratio_h = simp(F3_h / F2_h)

show("F3_h/F2_h", ratio_h)

cert(
    "ratio shift",
    ratio_h,
    S + 2 * h + 1,
)


# ==============================================================================
# [24] DIRECT LAYER RECONSTRUCTION
# ==============================================================================

print("[24] LAYER RECONSTRUCTION")
print("-" * 78)

F2_from_orbit = simp(
    L0 + L1 - 1 - Delta
)

F3_from_orbit = simp(
    (S + 1) * F2_from_orbit
)

cert(
    "F2 from root orbit",
    F2_from_orbit,
    F2,
)

cert(
    "F3 from root orbit",
    F3_from_orbit,
    F3,
)


# ==============================================================================
# [25] SAFE PURE ORBIT LINEAR-COMBINATION SEARCH
# ==============================================================================

print("[25] PURE TRANSLATED-PRODUCT COMBINATION SEARCH")
print("-" * 78)

c0, c1, c2_sym, c3 = sp.symbols(
    "c0 c1 c2_sym c3"
)

pure_combo = (
    c0 * L0
    + c1 * L1
    + c2_sym * L2
    + c3
)

target_combo = F2

#
# IMPORTANT:
#
#   pure_combo - target_combo
#
# is rational in r.
#
# Clear denominators FIRST.
#
# The original script failed here because Poly() received r**(-2).
#

diff25 = sp.together(
    pure_combo - target_combo
)

num25, den25 = diff25.as_numer_denom()

poly25 = sp.Poly(
    sp.expand(num25),
    A,
    B,
    r,
)

solutions25 = sp.solve(
    poly25.coeffs(),
    [c0, c1, c2_sym, c3],
    dict=True,
)

print("  common denominator =")
sp.pprint(den25)
print()

print("  solutions =")
sp.pprint(solutions25)
print()

for sol in solutions25:
    reconstructed = simp(
        pure_combo.subs(sol)
    )

    show("reconstructed", reconstructed)

    cert(
        "pure orbit reconstruction",
        reconstructed,
        target_combo,
    )


# ==============================================================================
# [26] SAFE ORBIT + INVARIANT SEARCH
# ==============================================================================

print("[26] ORBIT + INVARIANT MINIMALITY SEARCH")
print("-" * 78)

a0, a1, aD, ac = sp.symbols(
    "a0 a1 aD ac"
)

combo2 = (
    a0 * L0
    + a1 * L1
    + aD * Delta
    + ac
)

diff26 = sp.together(
    combo2 - F2
)

num26, den26 = diff26.as_numer_denom()

poly26 = sp.Poly(
    sp.expand(num26),
    A,
    B,
    r,
)

solutions26 = sp.solve(
    poly26.coeffs(),
    [a0, a1, aD, ac],
    dict=True,
)

print("  common denominator =")
sp.pprint(den26)
print()

print("  solutions =")
sp.pprint(solutions26)
print()

for sol in solutions26:
    reconstructed = simp(
        combo2.subs(sol)
    )

    show("reconstructed", reconstructed)

    cert(
        "L0,L1,Delta reconstruction",
        reconstructed,
        F2,
    )


# ==============================================================================
# [27] GENERAL DISCRETE TRANSLATION ORBIT
# ==============================================================================

print("[27] GENERAL DISCRETE TRANSLATION ORBIT")
print("-" * 78)

def L(j):
    return simp(
        (A + j * r) * (B + j * r) / r**2
    )


Lm = L(-m)
Lp = L(m)

show("L(m)", Lp)
show("L(-m)", Lm)

cert(
    "L(m)=m^2+Sm+N",
    Lp,
    m**2 + S * m + N,
)

cert(
    "L(-m)=m^2-Sm+N",
    Lm,
    m**2 - S * m + N,
)


# ==============================================================================
# [28] CENTRAL DIFFERENCE
# ==============================================================================

print("[28] CENTRAL DIFFERENCE")
print("-" * 78)

central = simp(
    L(h) - L(-h)
)

show("L(h)-L(-h)", central)

cert(
    "central difference = 2hS",
    central,
    2 * h * S,
)


# ==============================================================================
# [29] DISCRETE S RECOVERY
# ==============================================================================

print("[29] DISCRETE S RECOVERY")
print("-" * 78)

S_from_orbit = simp(
    L1 - L0 - 1
)

show("S from L1-L0-1", S_from_orbit)

cert(
    "S=L1-L0-1",
    S_from_orbit,
    S,
)


# ==============================================================================
# [30] GAP RECOVERY
# ==============================================================================

print("[30] GAP RECOVERY")
print("-" * 78)

cert(
    "Delta=(A-B)^2/r^2",
    Delta,
    (A - B)**2 / r**2,
)


# ==============================================================================
# [31] HIGHER TRANSLATED PRODUCTS
# ==============================================================================

print("[31] HIGHER TRANSLATED PRODUCTS")
print("-" * 78)

L3 = L(3)
L4 = L(4)
L5 = L(5)

show("L3", L3)
show("L4", L4)
show("L5", L5)

cert(
    "L3=N+3S+9",
    L3,
    N + 3*S + 9,
)

cert(
    "L4=N+4S+16",
    L4,
    N + 4*S + 16,
)

cert(
    "L5=N+5S+25",
    L5,
    N + 5*S + 25,
)


# ==============================================================================
# [32] SECOND DIFFERENCE AT HIGHER LEVELS
# ==============================================================================

print("[32] HIGHER ORBIT CURVATURE")
print("-" * 78)

cert(
    "L3-2L2+L1=2",
    L3 - 2*L2 + L1,
    2,
)

cert(
    "L4-2L3+L2=2",
    L4 - 2*L3 + L2,
    2,
)

cert(
    "L5-2L4+L3=2",
    L5 - 2*L4 + L3,
    2,
)


# ==============================================================================
# [33] FIRST DIFFERENCES
# ==============================================================================

print("[33] ORBIT FIRST DIFFERENCES")
print("-" * 78)

cert(
    "L1-L0=S+1",
    L1 - L0,
    S + 1,
)

cert(
    "L2-L1=S+3",
    L2 - L1,
    S + 3,
)

cert(
    "L3-L2=S+5",
    L3 - L2,
    S + 5,
)

cert(
    "L4-L3=S+7",
    L4 - L3,
    S + 7,
)


# ==============================================================================
# [34] RATIO-LADDER ANALOGUE
# ==============================================================================

print("[34] ROOT-ORBIT RATIO LADDER ANALOGUE")
print("-" * 78)

#
# The root-coordinate orbit gives an explicit analogue of the
# proposed homogeneous ratio ladder:
#
#     L_{m+1}-L_m = S+2m+1.
#
# This is additive rather than multiplicative.
#

for j in range(0, 5):
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
# [35] M1 RELATION
# ==============================================================================

print("[35] M1 RELATION")
print("-" * 78)

cert(
    "M1=L1",
    M1,
    L1,
)

cert(
    "M1-N=S+1",
    M1-N,
    S+1,
)


# ==============================================================================
# [36] STRUCTURAL SUMMARY
# ==============================================================================

print("[36] STRUCTURAL SUMMARY")
print("-" * 78)

print(
    r"""
  ROOT COORDINATES:

      A = A0
      B = B0

  NORMALIZATION:

      P = A/r
      Q = B/r

  SYMMETRIC CHANNELS:

      N     = AB/r^2
      S     = (A+B)/r
      Delta = (A-B)^2/r^2

  TRANSLATED PRODUCT ORBIT:

      L_m
        = (A+mr)(B+mr)/r^2
        = N + mS + m^2.

  FIRST DIFFERENCE:

      L_{m+1}-L_m
        = S + 2m + 1.

  SECOND DIFFERENCE:

      L_{m+2}-2L_{m+1}+L_m
        = 2.

  SPECTRAL REFLECTION:

      L_{-z}
        = z^2-Sz+N
        = (z-A/r)(z-B/r).

  KAPPA +1 TRANSLATION:

      Q(z-1)
        = (z-(A+r)/r)(z-(B+r)/r).

  HOMOGENEOUS LAYER:

      F2
        = 2N+S-Delta
        = L0+L1-1-Delta.

  COMPANION LAYER:

      F3
        = (S+1)F2
        = (L1-L0)F2.

  The remaining question is now specifically whether the
  higher homogeneous layers can be expressed using the same
  translated-product orbit L_m, rather than being independent
  objects.
"""
)


# ==============================================================================
# [37] FINAL AUDIT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 554 FINISHED")
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
The important result is:

    F2 = L0 + L1 - 1 - Delta

with

    L_m = (A0+mr)(B0+mr)/r^2.

The root-coordinate orbit itself has the exact affine ladder

    L_{m+1}-L_m = S+2m+1.

This is structurally related to the previously proposed
homogeneous-layer ratio ladder, but it is not the same object.

Therefore the next experiment should compare actual higher
homogeneous layers against the translated-product basis:

    L0, L1, L2, L3, ...

Specifically test F4/F5/F6, when their genuine formulas are
available, against low-degree combinations such as:

    L_n
    L_n + L_{n+1}
    L_{n+1}-L_n
    L_n + L_{n+1} - Delta
    products of consecutive L_n.

Do not infer F4/F5 from the proposed ratio ladder.
Use their actual homogeneous-layer formulas if available.
"""
)

# Explicit reminder that all polynomial searches are denominator-safe.
print()
print("NOTE: rational identities were converted to numerator")
print("      polynomial identities before every sp.Poly() search.")