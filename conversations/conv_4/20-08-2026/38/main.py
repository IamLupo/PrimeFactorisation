#!/usr/bin/env python3

import sympy as sp

print("=" * 78)
print("EXPERIMENT 514 START")
print("=" * 78)
print("CENTERED KAPPA SEQUENCE / EVEN-ODD RECURRENCE / CONIC INVARIANT")
print("=" * 78)


# ============================================================================
# SYMBOLS
# ============================================================================

U, V, r = sp.symbols("U V r", nonzero=True)
n = sp.symbols("n", integer=True, nonnegative=True)
N, S, Delta, K = sp.symbols("N S Delta K")
x = sp.symbols("x")


def fact(expr):
    return sp.factor(sp.expand(expr))


def simp(expr):
    return sp.simplify(sp.expand(expr))


def zero(expr):
    return sp.expand(expr) == 0


def show_cert(name, expr):
    d = fact(expr)
    print(f"  {name}")
    print(f"    difference = {d}")
    print(f"    PASS = {d == 0}")
    print()


# ============================================================================
# [1] HISTORICAL / NORMALIZED COORDINATES
# ============================================================================

print()
print("[1] NORMALIZED HISTORICAL COORDINATES")
print("-" * 78)

P = fact((U - V) / r)
Q = fact((U + V) / r)

N_expr = fact(P * Q)
S_expr = fact(P + Q)
Delta_expr = fact((P - Q) ** 2)

print(f"  P = {P}")
print(f"  Q = {Q}")
print(f"  N = {N_expr}")
print(f"  S = {S_expr}")
print(f"  Delta = {Delta_expr}")
print()


# ============================================================================
# [2] ORIGINAL FOUR-EXPONENTIAL KAPPA SEQUENCE
# ============================================================================

print()
print("[2] ORIGINAL FOUR-EXPONENTIAL SEQUENCE")
print("-" * 78)

# The normalized F_n kernel:
#
#   F_n = P*(Q+1)^n + Q*(P+1)^n
#         - (Q+1)*P^n - (P+1)*Q^n
#
# This agrees with the established KAPPA kernel.

c = [
    -(Q + 1),
    -(P + 1),
    P,
    Q,
]

lam = [
    P,
    Q,
    P + 1,
    Q + 1,
]


def F(k):
    return fact(sum(c[i] * lam[i] ** k for i in range(4)))


for k in range(0, 9):
    print(f"  F_{k} = {F(k)}")

print()


# ============================================================================
# [3] FOUR-BASE CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[3] FOUR-BASE CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

z = sp.symbols("z")

chi = fact(sp.prod(z - a for a in lam))

print("  chi(z) =")
print(f"    {chi}")
print()


# ============================================================================
# [4] CENTER THE SPECTRUM
# ============================================================================

print()
print("[4] SPECTRAL CENTER")
print("-" * 78)

center = fact((P + Q + 1) / 2)

print(f"  center = {center}")
print()

centered_lam = [fact(a - center) for a in lam]

for i, a in enumerate(centered_lam, start=1):
    print(f"  mu_{i} = {a}")

print()


# ============================================================================
# [5] CENTERED ROOT STRUCTURE
# ============================================================================

print()
print("[5] CENTERED ROOT PAIRING")
print("-" * 78)

print("  Testing:")
print("    mu_4 = -mu_1")
print("    mu_3 = -mu_2")
print()

show_cert(
    "mu_4 + mu_1",
    centered_lam[3] + centered_lam[0],
)

show_cert(
    "mu_3 + mu_2",
    centered_lam[2] + centered_lam[1],
)


# ============================================================================
# [6] CENTERED CHARACTERISTIC POLYNOMIAL
# ============================================================================

print()
print("[6] CENTERED CHARACTERISTIC POLYNOMIAL")
print("-" * 78)

xi = sp.symbols("xi")

chi_centered = fact(
    chi.subs(z, xi + center)
)

print("  chi_centered(xi) =")
print(f"    {chi_centered}")
print()


# Expected centered form:
#
#   [ (xi+1/2)^2 - V^2/r^2 ]
#   [ (xi-1/2)^2 - V^2/r^2 ]


centered_target = fact(
    (
        (xi + sp.Rational(1, 2)) ** 2
        - V ** 2 / r ** 2
    )
    *
    (
        (xi - sp.Rational(1, 2)) ** 2
        - V ** 2 / r ** 2
    )
)

show_cert(
    "centered quartic",
    chi_centered - centered_target,
)


# ============================================================================
# [7] EVEN POLYNOMIAL FORM
# ============================================================================

print()
print("[7] EVEN CENTERED POLYNOMIAL")
print("-" * 78)

even_centered = fact(
    xi**4
    - (
        2 * V**2 / r**2
        + sp.Rational(1, 2)
    ) * xi**2
    + (
        V**2 / r**2
        - sp.Rational(1, 4)
    )**2
)

print("  chi_centered(xi) =")
print(f"    {even_centered}")
print()

show_cert(
    "even quartic expansion",
    chi_centered - even_centered,
)


# ============================================================================
# [8] RECURRENCE FOR CENTERED EXPONENTIALS
# ============================================================================

print()
print("[8] CENTERED FOURTH-ORDER RECURRENCE")
print("-" * 78)

mu = centered_lam

# Elementary symmetric coefficients of the centered roots.
e1 = fact(sum(mu))
e2 = fact(sum(mu[i] * mu[j] for i in range(4) for j in range(i + 1, 4)))
e3 = fact(
    sum(
        mu[i] * mu[j] * mu[k]
        for i in range(4)
        for j in range(i + 1, 4)
        for k in range(j + 1, 4)
    )
)
e4 = fact(sp.prod(mu))

print(f"  e1 = {e1}")
print(f"  e2 = {e2}")
print(f"  e3 = {e3}")
print(f"  e4 = {e4}")
print()

# Recurrence:
#
# C[n+4] - e1 C[n+3] + e2 C[n+2] - e3 C[n+1] + e4 C[n] = 0
#
# Because e1=e3=0, it must reduce to:
#
# C[n+4] + e2 C[n+2] + e4 C[n] = 0


show_cert("e1 = 0", e1)
show_cert("e3 = 0", e3)


# ============================================================================
# [9] RECURRENCE COEFFICIENTS IN V^2
# ============================================================================

print()
print("[9] RECURRENCE COEFFICIENTS")
print("-" * 78)

A_rec = fact(e2)
B_rec = fact(e4)

print(f"  recurrence:")
print("    C[n+4] + A*C[n+2] + B*C[n] = 0")
print()
print(f"  A = {A_rec}")
print(f"  B = {B_rec}")
print()

A_target = fact(
    -(
        2 * V**2 / r**2
        + sp.Rational(1, 2)
    )
)

B_target = fact(
    (
        V**2 / r**2
        - sp.Rational(1, 4)
    )**2
)

show_cert(
    "A depends only on V^2",
    A_rec - A_target,
)

show_cert(
    "B depends only on V^2",
    B_rec - B_target,
)


# ============================================================================
# [10] U-ELIMINATION OF THE RECURRENCE
# ============================================================================

print()
print("[10] U-ELIMINATION")
print("-" * 78)

print(f"  A contains U = {A_rec.has(U)}")
print(f"  B contains U = {B_rec.has(U)}")
print()

u_elimination = not A_rec.has(U) and not B_rec.has(U)

print(f"  recurrence U-free = {u_elimination}")
print()


# ============================================================================
# [11] CONIC SUBSTITUTION
# ============================================================================

print()
print("[11] CONIC INVARIANT SUBSTITUTION")
print("-" * 78)

# K = U^2 - V^2
#
# Hence
#
#   V^2 = U^2 - K.
#
# We test whether the recurrence remains structurally
# separable into the conic invariant K and the centered
# coordinate U.

K_expr = fact(U**2 - V**2)

V2_from_K = fact(U**2 - K)

A_K = fact(A_rec.subs(V**2, V2_from_K))
B_K = fact(B_rec.subs(V**2, V2_from_K))

print(f"  K = {K_expr}")
print(f"  A(V^2=U^2-K) = {A_K}")
print(f"  B(V^2=U^2-K) = {B_K}")
print()


# ============================================================================
# [12] ORIGINAL/SHIFTED PRODUCT RELATION
# ============================================================================

print()
print("[12] SHIFTED PRODUCT CHANNEL")
print("-" * 78)

M0 = fact(P * Q)
M1 = fact((P + 1) * (Q + 1))

print(f"  M0 = {M0}")
print(f"  M1 = {M1}")
print()

show_cert(
    "M0 = (U^2-V^2)/r^2",
    M0 - K_expr / r**2,
)

show_cert(
    "M1 = ((U+r)^2-V^2)/r^2",
    M1 - ((U + r)**2 - V**2) / r**2,
)


# ============================================================================
# [13] DIFFERENCE OF SHIFTED PRODUCTS
# ============================================================================

print()
print("[13] TRANSLATION DIFFERENCE")
print("-" * 78)

M_difference = fact(M1 - M0)

print(f"  M1-M0 = {M_difference}")
print()

show_cert(
    "M1-M0 = 2U/r + 1",
    M_difference - (2 * U / r + 1),
)


# ============================================================================
# [14] CENTERED ROOT SQUARES
# ============================================================================

print()
print("[14] CENTERED ROOT SQUARES")
print("-" * 78)

mu_squared = [fact(a**2) for a in mu]

for i, a in enumerate(mu_squared, start=1):
    print(f"  mu_{i}^2 = {a}")

print()


# The two distinct squared roots must be
#
#   (V/r + 1/2)^2
#   (V/r - 1/2)^2


sq_target_1 = fact((V / r + sp.Rational(1, 2))**2)
sq_target_2 = fact((V / r - sp.Rational(1, 2))**2)

show_cert(
    "mu_1^2 - (V/r+1/2)^2",
    mu[0]**2 - sq_target_1,
)

show_cert(
    "mu_2^2 - (V/r-1/2)^2",
    mu[1]**2 - sq_target_2,
)


# ============================================================================
# [15] EVEN / ODD SEQUENCE SPLIT
# ============================================================================

print()
print("[15] EVEN / ODD CENTERED SPECTRAL MODES")
print("-" * 78)

# Define abstract centered spectral moments:
#
#   E_m = sum c_i * mu_i^m
#
# The characteristic polynomial is even, so the roots occur
# in ± pairs. We inspect parity structure symbolically.

weights = c

C = lambda m: fact(
    sum(weights[i] * centered_lam[i]**m for i in range(4))
)

for m in range(0, 9):
    print(f"  C_{m} = {C(m)}")

print()


# ============================================================================
# [16] PARITY CHECK
# ============================================================================

print()
print("[16] PARITY CHECK")
print("-" * 78)

# A centered spectrum with ± roots admits a natural decomposition
# into even and odd powers. We test the direct symbolic relation
# against the roots rather than assuming the weights are symmetric.

parity_rows = []

for m in range(0, 9):
    cpos = fact(C(m))
    cneg = fact(
        sum(
            weights[i] * (-centered_lam[i])**m
            for i in range(4)
        )
    )

    even_part = fact((cpos + cneg) / 2)
    odd_part = fact((cpos - cneg) / 2)

    parity_rows.append((m, even_part, odd_part))

for m, even_part, odd_part in parity_rows:
    print(f"  m={m}")
    print(f"    even = {even_part}")
    print(f"    odd  = {odd_part}")

print()


# ============================================================================
# [17] WEIGHT TRANSFORMATION
# ============================================================================

print()
print("[17] WEIGHT SYMMETRY UNDER CENTERING")
print("-" * 78)

for i in range(4):
    print(f"  c_{i+1} = {fact(weights[i])}")

print()

# Test pair differences:
#
# c(P) vs c(Q)
# c(P+1) vs c(Q+1)

w1, w2, w3, w4 = weights

print(f"  c(P)-c(Q) = {fact(w1-w2)}")
print(f"  c(P+1)-c(Q+1) = {fact(w3-w4)}")
print(f"  c(P)+c(Q) = {fact(w1+w2)}")
print(f"  c(P+1)+c(Q+1) = {fact(w3+w4)}")
print()


# ============================================================================
# [18] POSSIBLE U-FREE OBSERVABLES
# ============================================================================

print()
print("[18] U-FREE OBSERVABLE SEARCH")
print("-" * 78)

observables = {
    "e1": e1,
    "e2": e2,
    "e3": e3,
    "e4": e4,
    "M0": M0,
    "M1": M1,
    "M1-M0": M_difference,
    "Delta": Delta_expr,
}

for name, expr in observables.items():
    print(f"  {name}")
    print(f"    = {expr}")
    print(f"    contains U = {expr.has(U)}")

print()


# ============================================================================
# [19] KEY SYMBOLIC CANDIDATE
# ============================================================================

print()
print("[19] KEY RECURRENCE COMBINATION")
print("-" * 78)

# Since
#
#   A = -(2V^2/r^2 + 1/2)
#
# we obtain directly
#
#   V^2/r^2 = -A/2 - 1/4.
#
# This is a candidate for extracting the gap square
# directly from a centered recurrence coefficient.

V2_from_A = fact(-A_rec / 2 - sp.Rational(1, 4))

print("  reconstructed V^2/r^2 =")
print(f"    {V2_from_A}")
print()

show_cert(
    "V^2/r^2 reconstruction",
    V2_from_A - V**2 / r**2,
)


# ============================================================================
# [20] GAP SQUARE FROM RECURRENCE
# ============================================================================

print()
print("[20] GAP-SQUARE FROM RECURRENCE")
print("-" * 78)

gap_candidate = fact(
    4 * V2_from_A
)

gap_target = fact(
    (P - Q)**2
)

print("  candidate =")
print(f"    {gap_candidate}")

print()
print("  target =")
print(f"    {gap_target}")

print()

show_cert(
    "4*(V^2/r^2) = Delta",
    gap_candidate - gap_target,
)


# ============================================================================
# [21] CONIC RELATION
# ============================================================================

print()
print("[21] CONIC RELATION")
print("-" * 78)

# In historical coordinates:
#
#   U^2 - V^2 = K.
#
# After obtaining V^2, we have:
#
#   U^2 = K + V^2.
#
# The experiment tests the reconstruction algebraically.

U2_candidate = fact(
    K_expr + V**2
)

show_cert(
    "U^2 = K + V^2",
    U2_candidate - U**2,
)


# ============================================================================
# [22] INFORMATION FLOW
# ============================================================================

print()
print("[22] INFORMATION FLOW")
print("-" * 78)

print(
r"""
The centered spectrum has the exact form

    mu =
      { -(V/r+1/2),
        -(V/r-1/2),
          V/r-1/2,
          V/r+1/2 }.

Therefore its characteristic polynomial is even:

    chi_c(xi)
      =
      ((xi+1/2)^2 - V^2/r^2)
      ((xi-1/2)^2 - V^2/r^2).

The fourth-order recurrence consequently contains only
the coefficients

    A = -(2 V^2/r^2 + 1/2)
    B = (V^2/r^2 - 1/4)^2.

Hence

    V^2/r^2 = -A/2 - 1/4.

This is significant because the gap information appears
as a recurrence coefficient rather than as an explicitly
inserted factor difference.

The remaining question is now very specific:

    Can A or B be computed from the original KAPPA sequence
    without explicitly constructing P and Q?

If yes, then

    recurrence coefficient
        -> V^2/r^2
        -> Delta
        -> factor gap

and the historical conic supplies

    U^2-V^2 = K

to reconstruct the symmetric coordinate.
"""
)


# ============================================================================
# [23] FINAL AUDIT
# ============================================================================

print()
print("[23] FINAL SYMBOLIC AUDIT")
print("-" * 78)

audits = {
    "centered roots pair": zero(centered_lam[3] + centered_lam[0])
    and zero(centered_lam[2] + centered_lam[1]),

    "centered characteristic": zero(
        chi_centered - centered_target
    ),

    "even quartic": zero(
        chi_centered - even_centered
    ),

    "e1 zero": zero(e1),

    "e3 zero": zero(e3),

    "A coefficient": zero(
        A_rec - A_target
    ),

    "B coefficient": zero(
        B_rec - B_target
    ),

    "V2 reconstruction": zero(
        V2_from_A - V**2 / r**2
    ),

    "gap reconstruction": zero(
        gap_candidate - gap_target
    ),

    "conic reconstruction": zero(
        U2_candidate - U**2
    ),
}

for name, value in audits.items():
    print(f"  {name:<28} = {value}")

overall = all(audits.values())

print()
print(f"  OVERALL EXACT AUDIT = {overall}")


# ============================================================================
# [24] NEXT TARGET
# ============================================================================

print()
print("[24] NEXT RESEARCH TARGET")
print("-" * 78)

print(
r"""
The important result is not another expression for Delta.

The centered four-base spectrum produces an even quartic,
and its recurrence has the form

    C[n+4] + A C[n+2] + B C[n] = 0

with

    A = -(2 V^2/r^2 + 1/2).

Therefore the gap-square is encoded directly in a recurrence
coefficient:

    Delta = 4 V^2/r^2 = -2A - 1.

The next experiment should therefore construct the recurrence
DIRECTLY from the observable KAPPA sequence F_n and ask:

    Can its recurrence coefficient A be extracted from F_n
    alone?

Specifically:

    F_n
      -> recurrence coefficients
      -> A
      -> Delta
      -> historical conic K
      -> U^2
      -> U
      -> p,q.

This is a substantially sharper upstream target than searching
arbitrary combinations of F_n, H2(x), continued fractions,
or divisor pairs.
"""
)

print()
print("=" * 78)
print("EXPERIMENT 514 FINISHED")
print("=" * 78)
