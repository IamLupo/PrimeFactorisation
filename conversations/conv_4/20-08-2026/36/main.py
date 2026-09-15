#!/usr/bin/env python3

"""
EXPERIMENT 512
SCALE-NORMALIZED HISTORICAL COORDINATES / KAPPA OPERATOR ALGEBRA

Pure symbolic experiment.
No data analysis.
No numerical fitting.
No divisor enumeration.
No continued fractions.

Definitions:

    A0 = 2*y - 2*x + 3
    B0 = 2*y + 2*x - 3

For the historical scale r:

    r = 1  when 4*v-1 = N
    r = 3  when 4*v-1 = 3*N

Define normalized coordinates

    P = A0/r
    Q = B0/r

so universally

    N = P*Q
    S = P+Q
    Delta = (P-Q)^2.

The KAPPA shifted kernel is

    H2(t) = 6*N - S^2 + S*t.

The experiment searches for exact symbolic identities relating
historical coordinates to the KAPPA operator.
"""

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

P, Q = sp.symbols("P Q")
A, B = sp.symbols("A B")
r = sp.symbols("r", nonzero=True)
t = sp.symbols("t")
u, v, w = sp.symbols("u v w")

X, Y = sp.symbols("X Y")


# =============================================================================
# HELPERS
# =============================================================================

def fact(expr):
    """Expanded symbolic expression followed by factorization."""
    return sp.factor(sp.expand(expr))


def diff(expr1, expr2):
    """Exact symbolic difference."""
    return sp.expand(expr1 - expr2)


def passes(expr1, expr2):
    """True iff symbolic difference is exactly zero."""
    return sp.expand(expr1 - expr2) == 0


def print_identity(label, lhs, rhs):
    d = fact(lhs - rhs)
    print(f"  {label}")
    print(f"    difference = {d}")
    print(f"    PASS = {d == 0}")
    print()


# =============================================================================
# MODERN COORDINATES
# =============================================================================

Npq = P * Q
Spq = P + Q
Deltapq = (P - Q) ** 2
M1pq = (P + 1) * (Q + 1)


def H2(shift):
    return sp.expand(
        6 * Npq
        - Spq ** 2
        + Spq * shift
    )


# =============================================================================
# START
# =============================================================================

print("=" * 78)
print("EXPERIMENT 512 START")
print("=" * 78)
print("SCALE-NORMALIZED HISTORICAL COORDINATES / KAPPA OPERATOR ALGEBRA")
print("=" * 78)


# =============================================================================
# [1] NORMALIZED COORDINATES
# =============================================================================

print()
print("[1] NORMALIZED HISTORICAL COORDINATES")
print("-" * 78)

print("  P = A0/r")
print("  Q = B0/r")
print()
print("  N     = P*Q")
print("  S     = P+Q")
print("  Delta = (P-Q)^2")
print("  M1    = (P+1)(Q+1)")
print()

print_identity("N", Npq, P * Q)
print_identity("S", Spq, P + Q)
print_identity("Delta", Deltapq, (P - Q) ** 2)
print_identity("M1", M1pq, (P + 1) * (Q + 1))


# =============================================================================
# [2] KAPPA KERNEL
# =============================================================================

print()
print("[2] GENERAL KAPPA SHIFTED KERNEL")
print("-" * 78)

H_general = fact(H2(t))

print("  H2(t) =")
print(f"    {H_general}")
print()


# =============================================================================
# [3] GENERAL LINEAR HISTORICAL SHIFT
# =============================================================================

print()
print("[3] GENERAL LINEAR SHIFT")
print("-" * 78)

T = u * P + v * Q + w
H_linear = fact(H2(T))

print("  t = u*P + v*Q + w")
print()
print("  H2(t) =")
print(f"    {H_linear}")
print()

poly = sp.Poly(sp.expand(H_linear), P, Q)

print(f"  coefficient(P^2) = {poly.coeff_monomial(P**2)}")
print(f"  coefficient(P*Q) = {poly.coeff_monomial(P*Q)}")
print(f"  coefficient(Q^2) = {poly.coeff_monomial(Q**2)}")
print(f"  coefficient(P)   = {poly.coeff_monomial(P)}")
print(f"  coefficient(Q)   = {poly.coeff_monomial(Q)}")
print()


# =============================================================================
# [4] UNIQUE SUM SHIFT
# =============================================================================

print()
print("[4] SUM SHIFT")
print("-" * 78)

H_sum = fact(H2(P + Q))

print("  t = P+Q")
print(f"  H2(P+Q) = {H_sum}")
print()

print_identity(
    "H2(P+Q)=6N",
    H_sum,
    6 * Npq,
)


# =============================================================================
# [5] DIFFERENCE SHIFTS
# =============================================================================

print()
print("[5] DIFFERENCE SHIFTS")
print("-" * 78)

H_qp = fact(H2(Q - P))
H_pq = fact(H2(P - Q))

print(f"  H2(Q-P) = {H_qp}")
print(f"  H2(P-Q) = {H_pq}")
print()

print_identity(
    "H2(Q-P)",
    H_qp,
    -2 * P * (P - 2 * Q),
)

print_identity(
    "H2(P-Q)",
    H_pq,
    2 * Q * (2 * P - Q),
)


# =============================================================================
# [6] OPPOSITE SHIFT SUM AND DIFFERENCE
# =============================================================================

print()
print("[6] OPPOSITE SHIFT OPERATOR")
print("-" * 78)

Hsum_diff = fact(H_qp + H_pq)
Hdiff_diff = fact(H_qp - H_pq)

print(f"  H2(Q-P)+H2(P-Q) = {Hsum_diff}")
print(f"  H2(Q-P)-H2(P-Q) = {Hdiff_diff}")
print()

print_identity(
    "sum identity",
    Hsum_diff,
    2 * (6 * Npq - Spq ** 2),
)

print_identity(
    "difference identity",
    Hdiff_diff,
    -2 * Spq * (P - Q),
)


# =============================================================================
# [7] GAP SQUARE
# =============================================================================

print()
print("[7] GAP-SQUARE EXTRACTION")
print("-" * 78)

gap_from_operator = fact(
    Hdiff_diff ** 2 / (4 * Spq ** 2)
)

print("  [H2(Q-P)-H2(P-Q)]^2 / [4*S^2]")
print(f"    = {gap_from_operator}")
print()

print_identity(
    "gap square",
    gap_from_operator,
    Deltapq,
)


# =============================================================================
# [8] OPERATOR PRODUCT
# =============================================================================

print()
print("[8] OPPOSITE-SHIFT PRODUCT")
print("-" * 78)

product_op = fact(H_qp * H_pq)

print("  H2(Q-P)*H2(P-Q) =")
print(f"    {product_op}")
print()


# =============================================================================
# [9] TARGET SEARCH
# =============================================================================

print()
print("[9] EXACT TARGET SEARCH")
print("-" * 78)

candidate_shifts = {
    "P+Q": P + Q,
    "Q-P": Q - P,
    "P-Q": P - Q,
    "P": P,
    "Q": Q,
    "P+1": P + 1,
    "Q+1": Q + 1,
    "P+Q+1": P + Q + 1,
    "P+Q-1": P + Q - 1,
    "Q-P+1": Q - P + 1,
    "P-Q+1": P - Q + 1,
}

targets = {
    "N": Npq,
    "S": Spq,
    "Delta": Deltapq,
    "M1": M1pq,
    "6N": 6 * Npq,
    "S^2": Spq ** 2,
    "2N-Delta": 2 * Npq - Deltapq,
    "N+S+1": Npq + Spq + 1,
}

for name, shift in candidate_shifts.items():
    expr = fact(H2(shift))

    hits = []
    for target_name, target_expr in targets.items():
        if passes(expr, target_expr):
            hits.append(target_name)

    print(f"  t={name}")
    print(f"    H2 = {expr}")
    print(f"    exact target hits = {hits}")
    print()


# =============================================================================
# [10] QUADRATIC FORM ANALYSIS
# =============================================================================

print()
print("[10] QUADRATIC FORM OF H2(uP+vQ+w)")
print("-" * 78)

quad_matrix = sp.Matrix([
    [u - 1, (5 + u + v) / 2],
    [(5 + u + v) / 2, v - 1],
])

quad_det = fact(quad_matrix.det())

print("  quadratic matrix =")
print(quad_matrix)
print()

print(f"  determinant = {quad_det}")
print()

print("  rank-one condition:")
print(f"    {quad_det} = 0")
print()


# =============================================================================
# [11] SOLVE RANK-ONE CONDITION
# =============================================================================

print()
print("[11] RANK-ONE SHIFT FAMILIES")
print("-" * 78)

rank_solutions = sp.solve(
    sp.Eq(quad_det, 0),
    v,
    dict=True,
)

if not rank_solutions:
    print("  no symbolic solutions returned")
else:
    for sol in rank_solutions:
        expression = fact(H_linear.subs(sol))

        print(f"  relation: {sol}")
        print(f"  H2 = {expression}")
        print()


# =============================================================================
# [12] NORMALIZED UNNORMALIZED HISTORICAL COORDINATES
# =============================================================================

print()
print("[12] UNNORMALIZED HISTORICAL COORDINATES")
print("-" * 78)

P_ab = A / r
Q_ab = B / r

N_ab = fact(P_ab * Q_ab)
S_ab = fact(P_ab + Q_ab)

def H2_ab(shift):
    return fact(
        6 * N_ab
        - S_ab ** 2
        + S_ab * shift
    )


print("  N(A,B) =")
print(f"    {N_ab}")

print()
print("  S(A,B) =")
print(f"    {S_ab}")

print()

H_ab_sum = H2_ab((A + B) / r)
H_ab_diff = H2_ab((B - A) / r)
H_ab_rev = H2_ab((A - B) / r)

print("  H2((A+B)/r) =")
print(f"    {H_ab_sum}")
print()

print("  H2((B-A)/r) =")
print(f"    {H_ab_diff}")
print()

print("  H2((A-B)/r) =")
print(f"    {H_ab_rev}")
print()

print_identity(
    "normalized sum shift",
    H_ab_sum,
    6 * A * B / r ** 2,
)


# =============================================================================
# [13] SCALE 1 AND SCALE 3
# =============================================================================

print()
print("[13] BRANCH SCALE CONSEQUENCES")
print("-" * 78)

for scale in (1, 3):
    sum_expr = fact(H_ab_sum.subs(r, scale))
    diff_expr = fact(H_ab_diff.subs(r, scale))
    rev_expr = fact(H_ab_rev.subs(r, scale))

    print(f"  r={scale}")
    print(f"    H2((A+B)/r) = {sum_expr}")
    print(f"    H2((B-A)/r) = {diff_expr}")
    print(f"    H2((A-B)/r) = {rev_expr}")
    print()


# =============================================================================
# [14] HISTORICAL x,y SUBSTITUTION
# =============================================================================

print()
print("[14] HISTORICAL x,y SUBSTITUTION")
print("-" * 78)

A0_xy = 2 * Y - 2 * X + 3
B0_xy = 2 * Y + 2 * X - 3

print("  A0 = 2Y-2X+3")
print("  B0 = 2Y+2X-3")
print()

for scale in (1, 3):
    expr_sum = fact(
        H2_ab(
            (A + B) / r
        ).subs({
            A: A0_xy,
            B: B0_xy,
            r: scale,
        })
    )

    expr_diff = fact(
        H2_ab(
            (B - A) / r
        ).subs({
            A: A0_xy,
            B: B0_xy,
            r: scale,
        })
    )

    print(f"  r={scale}")
    print(f"    sum shift  = {expr_sum}")
    print(f"    diff shift = {expr_diff}")
    print()


# =============================================================================
# [15] HISTORICAL GAP OPERATOR
# =============================================================================

print()
print("[15] HISTORICAL GAP OPERATOR")
print("-" * 78)

gap_ab = fact(
    ((A - B) / r) ** 2
)

operator_gap_ab = fact(
    (
        H_ab_diff - H_ab_rev
    ) ** 2
    /
    (
        4 * S_ab ** 2
    )
)

print("  Delta =")
print(f"    {gap_ab}")
print()

print("  operator reconstruction =")
print(f"    {operator_gap_ab}")
print()

print_identity(
    "operator reconstruction of Delta",
    operator_gap_ab,
    gap_ab,
)


# =============================================================================
# [16] M1
# =============================================================================

print()
print("[16] M1 IN NORMALIZED HISTORICAL COORDINATES")
print("-" * 78)

M1_ab = fact(
    (A / r + 1) * (B / r + 1)
)

print(f"  M1 = {M1_ab}")
print()

print_identity(
    "M1=N+S+1",
    M1_ab,
    N_ab + S_ab + 1,
)


# =============================================================================
# [17] CENTRAL SYMBOLIC IDENTITY
# =============================================================================

print()
print("[17] CENTRAL SYMBOLIC IDENTITY")
print("-" * 78)

print(
    r"""
After normalization, define

    P = A0/r
    Q = B0/r.

Then

    N = P Q
    S = P+Q
    Delta = (P-Q)^2.

The shifted KAPPA operator is

    H2(t)
      = 6PQ-(P+Q)^2+(P+Q)t.

The distinguished shifts are

    t = P+Q,
    t = Q-P,
    t = P-Q.

They give

    H2(P+Q) = 6PQ,

    H2(Q-P) = -2P(P-2Q),

    H2(P-Q) = 2Q(2P-Q).

Furthermore,

    H2(Q-P)-H2(P-Q)
      = -2(P+Q)(P-Q),

and consequently

    [H2(Q-P)-H2(P-Q)]^2
    -------------------------------- = (P-Q)^2.
              4(P+Q)^2

Thus the historical factor-gap coordinate is precisely
the normalized antisymmetric component of the KAPPA
shift operator.
"""
)


# =============================================================================
# [18] FINAL AUDIT
# =============================================================================

print()
print("[18] FINAL SYMBOLIC AUDIT")
print("-" * 78)

basic_ok = (
    passes(Npq, P * Q)
    and passes(Spq, P + Q)
    and passes(Deltapq, (P - Q) ** 2)
    and passes(M1pq, (P + 1) * (Q + 1))
)

operator_ok = (
    passes(H_sum, 6 * Npq)
    and passes(
        H_qp,
        -2 * P * (P - 2 * Q),
    )
    and passes(
        H_pq,
        2 * Q * (2 * P - Q),
    )
)

gap_ok = passes(
    gap_from_operator,
    Deltapq,
)

print(f"  normalized coordinates = {basic_ok}")
print(f"  operator identities    = {operator_ok}")
print(f"  gap extraction         = {gap_ok}")

overall = basic_ok and operator_ok and gap_ok

print()
print(f"  OVERALL SYMBOLIC AUDIT = {overall}")

print()
print("  NEXT SYMBOLIC TARGET:")
print("    Find whether the historical conic equation")
print("    generates a recurrence/operator in P,Q whose")
print("    characteristic structure matches the KAPPA")
print("    exponential bases.")
print()
print("    In particular, investigate whether")
print()
print("        P, Q, P+1, Q+1")
print()
print("    arise naturally from the historical conic")
print("    transformations rather than being inserted")
print("    after the factor pair is already known.")


print()
print("=" * 78)
print("EXPERIMENT 512 FINISHED")
print("=" * 78)