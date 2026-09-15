#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 553
# ==============================================================================
# 2020 FACTOR COORDINATES -> HOMOGENEOUS LAYER FACTORIZATION
#
# Goal:
#
#   Work directly with
#
#       A0 = 2y - 2x + 3
#       B0 = 2y + 2x - 3
#
#   instead of treating S as an unknown.
#
# Known normalized coordinates:
#
#       P = A0/r
#       Q = B0/r
#
#       N     = P Q
#       S     = P + Q
#       Delta = (P-Q)^2
#
# Known homogeneous layers:
#
#       F2 = 6N - S^2 + S
#       F3 = (S+1)F2
#
# The experiment searches for exact representations of F2 and F3
# using:
#
#       A0, B0
#       A0+B0
#       A0-B0
#       A0*B0
#       A0+r
#       B0+r
#       (A0+r)(B0+r)
#
# and their low-degree products.
#
# No numerical fitting.
# No divisor enumeration.
# No factor recovery from N.
#
# ==============================================================================


print("=" * 78)
print("EXPERIMENT 553 START")
print("=" * 78)
print("2020 FACTOR COORDINATES -> HOMOGENEOUS LAYER FACTORIZATION")
print()


# ==============================================================================
# SYMBOLS
# ==============================================================================

A, B, r, z = sp.symbols("A B r z", nonzero=True)
x, y = sp.symbols("x y")

failures = 0


# ==============================================================================
# HELPERS
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
# [1] HISTORICAL COORDINATES
# ==============================================================================

print("[1] HISTORICAL FACTOR COORDINATES")
print("-" * 78)

A0 = 2*y - 2*x + 3
B0 = 2*y + 2*x - 3

print("  A0 =")
sp.pprint(A0)
print()
print("  B0 =")
sp.pprint(B0)
print()

cert("A0 definition", A0, 2*y - 2*x + 3)
cert("B0 definition", B0, 2*y + 2*x - 3)


# ==============================================================================
# [2] NORMALIZED ROOT COORDINATES
# ==============================================================================

print("[2] NORMALIZED ROOT COORDINATES")
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
# [3] BASIC SYMMETRIC IDENTITIES
# ==============================================================================

print("[3] BASIC A,B SYMMETRIC IDENTITIES")
print("-" * 78)

cert(
    "S = (A+B)/r",
    S,
    (A + B) / r,
)

cert(
    "N = AB/r^2",
    N,
    A * B / r**2,
)

cert(
    "Delta = (A-B)^2/r^2",
    Delta,
    (A - B)**2 / r**2,
)

cert(
    "M1 = (A+r)(B+r)/r^2",
    M1,
    (A + r) * (B + r) / r**2,
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

show("F2", F2)
show("F3", F3)

cert(
    "F3=(S+1)F2",
    F3,
    (S + 1) * F2,
)


# ==============================================================================
# [5] F2 PURE A,B FORM
# ==============================================================================

print("[5] F2 PURE A,B FORM")
print("-" * 78)

F2_AB = simp(
    6*(A*B/r**2)
    - ((A+B)/r)**2
    + (A+B)/r
)

show("F2(A,B)", F2_AB)

expected_F2_AB = simp(
    (
        4*A*B
        - A**2
        - B**2
        + r*A
        + r*B
    ) / r**2
)

show("expected expanded form", expected_F2_AB)

cert(
    "F2 = [4AB-A^2-B^2+r(A+B)]/r^2",
    F2_AB,
    expected_F2_AB,
)


# ==============================================================================
# [6] F2 IN SUM / GAP COORDINATES
# ==============================================================================

print("[6] F2 IN A+B AND A-B")
print("-" * 78)

sum_ab = A + B
gap_ab = A - B

F2_sum_gap = simp(
    (
        2*A*B
        - (A-B)**2
        + r*(A+B)
    ) / r**2
)

show("F2 via sum/gap", F2_sum_gap)

cert(
    "F2 = [2AB-(A-B)^2+r(A+B)]/r^2",
    F2_AB,
    F2_sum_gap,
)


# ==============================================================================
# [7] F2 AS N, S, DELTA
# ==============================================================================

print("[7] F2 = 2N+S-DELTA")
print("-" * 78)

F2_NSD = simp(
    2*N + S - Delta
)

show("2N+S-Delta", F2_NSD)

cert(
    "F2 = 2N+S-Delta",
    F2,
    F2_NSD,
)


# ==============================================================================
# [8] F2 AS N AND M1
# ==============================================================================

print("[8] F2 USING N AND M1")
print("-" * 78)

#
# Since
#
#     M1 = N+S+1
#
# and
#
#     F2 = 2N+S-Delta,
#
# one gets
#
#     F2 = N + M1 - 1 - Delta.
#

candidate_N_M1 = simp(
    N + M1 - 1 - Delta
)

show("N+M1-1-Delta", candidate_N_M1)

cert(
    "F2 = N+M1-1-Delta",
    F2,
    candidate_N_M1,
)


# ==============================================================================
# [9] SHIFTED FACTOR FORM
# ==============================================================================

print("[9] SHIFTED FACTOR COORDINATES")
print("-" * 78)

Ap = A + r
Bp = B + r

show("A+r", Ap)
show("B+r", Bp)

cert(
    "M1",
    M1,
    Ap*Bp/r**2,
)


# ==============================================================================
# [10] SEARCH SIMPLE PRODUCTS
# ==============================================================================

print("[10] SIMPLE PRODUCT SEARCH")
print("-" * 78)

products = {
    "AB": A*B,
    "A(B+r)": A*(B+r),
    "B(A+r)": B*(A+r),
    "(A+r)(B+r)": (A+r)*(B+r),
    "A(A-B)": A*(A-B),
    "B(A-B)": B*(A-B),
    "(A+B)(A-B)": (A+B)*(A-B),
    "(A+r)(A-B)": (A+r)*(A-B),
    "(B+r)(A-B)": (B+r)*(A-B),
}

target = simp(
    r**2 * F2
)

print("  target = r^2 F2 =")
sp.pprint(target)
print()

for name, expr in products.items():
    print(f"  {name} =")
    sp.pprint(simp(expr))
    print()


# ==============================================================================
# [11] LINEAR COMBINATION SEARCH
# ==============================================================================

print("[11] LOW-DEGREE SYMMETRIC COMBINATION SEARCH")
print("-" * 78)

alpha, beta, gamma, delta_c = sp.symbols(
    "alpha beta gamma delta"
)

basis = [
    A*B,
    (A-B)**2,
    r*(A+B),
    r**2,
]

combo = (
    alpha*basis[0]
    + beta*basis[1]
    + gamma*basis[2]
    + delta_c*basis[3]
)

eq_poly = sp.Poly(
    sp.expand(
        combo - target
    ),
    A,
    B,
    r,
)

solution = sp.solve(
    eq_poly.coeffs(),
    [alpha, beta, gamma, delta_c],
    dict=True,
)

print("  basis:")
for item in basis:
    sp.pprint(item)

print()
print("  coefficient solution =")
sp.pprint(solution)
print()

if solution:
    combo_solution = simp(
        combo.subs(solution[0])
    )

    show(
        "reconstructed target",
        combo_solution,
    )

    cert(
        "linear symmetric reconstruction",
        combo_solution,
        target,
    )


# ==============================================================================
# [12] SEARCH SHIFTED-PRODUCT REPRESENTATIONS
# ==============================================================================

print("[12] SHIFTED-PRODUCT REPRESENTATION SEARCH")
print("-" * 78)

# Look for constants u,v,w such that
#
#   target = (A+u r)(B+v r) + w(A-B)^2
#
# or with the symmetric orientation.

u, v, w = sp.symbols("u v w")

candidate_shifted = (
    (A + u*r)*(B + v*r)
    + w*(A-B)**2
)

poly_shifted = sp.Poly(
    sp.expand(
        candidate_shifted - target
    ),
    A,
    B,
    r,
)

sol_shifted = sp.solve(
    poly_shifted.coeffs(),
    [u, v, w],
    dict=True,
)

print("  equation:")
sp.pprint(candidate_shifted)
print()

print("  solutions =")
sp.pprint(sol_shifted)
print()

for sol in sol_shifted:
    reconstructed = simp(
        candidate_shifted.subs(sol)
    )

    print("  reconstructed:")
    sp.pprint(reconstructed)
    print()

    cert(
        "shifted-product reconstruction",
        reconstructed,
        target,
    )


# ==============================================================================
# [13] F3 PURE A,B FORM
# ==============================================================================

print("[13] F3 PURE A,B FORM")
print("-" * 78)

F3_AB = simp(
    (A+B+r)/r * F2_AB
)

show("F3(A,B)", F3_AB)

cert(
    "F3 = (A+B+r)F2/r",
    F3_AB,
    F3,
)


# ==============================================================================
# [14] F3/F2
# ==============================================================================

print("[14] F3/F2")
print("-" * 78)

ratio = simp(
    F3_AB / F2_AB
)

show("F3/F2", ratio)

cert(
    "F3/F2 = (A+B+r)/r",
    ratio,
    (A+B+r)/r,
)


# ==============================================================================
# [15] DIRECT 2020 S CHANNEL
# ==============================================================================

print("[15] DIRECT 2020 S CHANNEL")
print("-" * 78)

S_direct = simp(
    ratio - 1
)

show("S = F3/F2-1", S_direct)

cert(
    "F3/F2-1 = (A+B)/r",
    S_direct,
    (A+B)/r,
)


# ==============================================================================
# [16] KAPPA QUADRATIC IN A,B
# ==============================================================================

print("[16] KAPPA QUADRATIC")
print("-" * 78)

Qz = simp(
    z**2 - S*z + N
)

Qz_AB = simp(
    (z - A/r)*(z - B/r)
)

show("Q(z)", Qz)
show("(z-A/r)(z-B/r)", Qz_AB)

cert(
    "Q factorization",
    Qz,
    Qz_AB,
)


# ==============================================================================
# [17] Q(z-1)
# ==============================================================================

print("[17] TRANSLATED KAPPA QUADRATIC")
print("-" * 78)

Q_shift = simp(
    Qz.subs(z, z-1)
)

Q_shift_fact = simp(
    (z - (A+r)/r) * (z - (B+r)/r)
)

show("Q(z-1)", Q_shift)
show("translated roots", Q_shift_fact)

cert(
    "Q(z-1) translated roots",
    Q_shift,
    Q_shift_fact,
)


# ==============================================================================
# [18] FOUR-BASE PRODUCT
# ==============================================================================

print("[18] FOUR-BASE CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

chi = simp(
    Qz * Q_shift
)

show("chi(z)", chi)

cert(
    "chi=Q(z)Q(z-1)",
    chi,
    Qz * Q_shift,
)


# ==============================================================================
# [19] HISTORICAL TRANSLATION
# ==============================================================================

print("[19] HISTORICAL FACTOR TRANSLATION")
print("-" * 78)

h = sp.symbols("h")

Ah = A + h*r
Bh = B + h*r

Lh = simp(
    Ah*Bh/r**2
)

show("L(h)", Lh)

cert(
    "L(h)=N+hS+h^2",
    Lh,
    N + h*S + h**2,
)

cert(
    "L(-z)=Q(z)",
    Lh.subs(h, -z),
    Qz,
)


# ==============================================================================
# [20] TRANSLATION OF F2
# ==============================================================================

print("[20] TRANSLATED F2 IN A,B")
print("-" * 78)

Ph = simp(Ah/r)
Qh = simp(Bh/r)

Nh = simp(Ph*Qh)
Sh = simp(Ph+Qh)

F2h = simp(
    6*Nh
    - Sh**2
    + Sh
)

F2_difference = simp(
    F2h - F2_AB
)

show("F2(h)", F2h)
show("F2(h)-F2", F2_difference)

cert(
    "F2(h)-F2 = 2h(S+h+1)",
    F2_difference,
    2*h*(S+h+1),
)


# ==============================================================================
# [21] F2 TRANSLATION COEFFICIENTS
# ==============================================================================

print("[21] F2 TRANSLATION COEFFICIENTS")
print("-" * 78)

hpoly = sp.Poly(
    sp.expand(F2h),
    h,
)

c2 = simp(hpoly.coeff_monomial(h**2))
c1 = simp(hpoly.coeff_monomial(h))
c0 = simp(hpoly.coeff_monomial(1))

show("h^2 coefficient", c2)
show("h coefficient", c1)
show("constant", c0)

cert(
    "F2 h^2 coefficient",
    c2,
    2,
)

cert(
    "F2 h coefficient",
    c1,
    2*(S+1),
)

cert(
    "F2 constant",
    c0,
    F2_AB,
)


# ==============================================================================
# [22] F2 ORBIT DISCRIMINANT
# ==============================================================================

print("[22] F2 ORBIT DISCRIMINANT")
print("-" * 78)

D_F2 = simp(
    c1**2 - 4*c2*c0
)

show("disc_h(F2(h))", D_F2)

cert(
    "disc_h(F2(h)) = 12 Delta + 4",
    D_F2,
    12*Delta + 4,
)


# ==============================================================================
# [23] HISTORICAL DELTA CHANNEL
# ==============================================================================

print("[23] HISTORICAL DISCRIMINANT CHANNEL")
print("-" * 78)

Delta_hist = simp(
    (A-B)**2/r**2
)

show("Delta", Delta_hist)

cert(
    "Delta=(A-B)^2/r^2",
    Delta,
    Delta_hist,
)

Delta_from_orbit = simp(
    (D_F2 - 4)/12
)

cert(
    "Delta=(disc(F2 orbit)-4)/12",
    Delta_from_orbit,
    Delta,
)


# ==============================================================================
# [24] M1 TRANSLATION RELATION
# ==============================================================================

print("[24] M1 TRANSLATION")
print("-" * 78)

M1_AB = simp(
    (A+r)*(B+r)/r**2
)

cert(
    "M1 = N+S+1",
    M1_AB,
    N+S+1,
)

cert(
    "M1-N = S+1",
    M1_AB-N,
    S+1,
)


# ==============================================================================
# [25] COMPARISON OF ALL NATURAL CHANNELS
# ==============================================================================

print("[25] CHANNEL COMPARISON")
print("-" * 78)

channels = {
    "product": A*B,
    "sum": A+B,
    "gap": A-B,
    "shifted product": (A+r)*(B+r),
    "F2 numerator": target,
}

for name, expr in channels.items():
    print(f"  {name}:")
    sp.pprint(sp.factor(expr))
    print()


# ==============================================================================
# [26] DIRECT STRUCTURAL CERTIFICATE
# ==============================================================================

print("[26] DIRECT STRUCTURAL CERTIFICATE")
print("-" * 78)

print(
    """
  The 2020 factor coordinates already provide:

      A0 = 2y-2x+3
      B0 = 2y+2x-3

  and after normalization:

      P = A0/r
      Q = B0/r.

  Consequently:

      N     = A0 B0/r^2
      S     = (A0+B0)/r
      Delta = (A0-B0)^2/r^2

  and:

      Q(z)
        = (z-A0/r)(z-B0/r).

  The known F2 layer is:

      F2
        = 4PQ-P^2-Q^2+P+Q

        = 2N+S-Delta.

  In A0,B0 coordinates:

      r^2 F2
        = 4A0B0
          - (A0-B0)^2
          + r(A0+B0).

  The translated factor orbit is:

      L(h)
        = (A0+hr)(B0+hr)/r^2

        = N+hS+h^2.

  Its spectral reflection is exactly:

      L(-z)=Q(z).

  Therefore the historical coordinates do not merely
  reconstruct S.

  They directly supply the ROOT COORDINATES from which the
  modern KAPPA quadratic is formed.
"""
)


# ==============================================================================
# [27] FINAL AUDIT
# ==============================================================================

print("=" * 78)
print("EXPERIMENT 553 FINISHED")
print("=" * 78)
print()

print(f"SYMBOLIC FAILURES = {failures}")
print(f"OVERALL EXACT AUDIT = {failures == 0}")

print()
print("=" * 78)
print("NEXT RESEARCH TARGET")
print("=" * 78)

print(
    """
The important upstream objects are now explicitly:

    A0 = 2y-2x+3
    B0 = 2y+2x-3

The next experiment should inspect the ORIGINAL 2020 algebra before
any modern KAPPA substitution.

Specifically search whether the homogeneous-layer expressions can
be written directly as low-degree combinations of:

    A0
    B0
    A0+B0
    A0-B0
    A0B0
    (A0+r)(B0+r)
    A0^2
    B0^2.

The strongest possible result would be an exact layer identity
whose natural factorization is already visible in A0 and B0.

That would establish that the modern KAPPA layers are not merely
functions of S,N,Delta, but are inherited from the original
2020 factor-coordinate algebra.
"""
)
