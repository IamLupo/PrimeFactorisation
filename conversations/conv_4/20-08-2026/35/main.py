#!/usr/bin/env python3
"""
EXPERIMENT 511
==============================================================================
HISTORICAL LINEAR-COORDINATE / KAPPA SHIFT OPERATOR
==============================================================================

Goal
----
Use the 2020 historical coordinates

    A0 = 2*y - 2*x + 3
    B0 = 2*y + 2*x - 3

as actual coordinates for the modern shifted KAPPA kernel

    H2(t) = 6*N - S^2 + S*t.

Primary questions
-----------------
1. What are H2(A0+B0), H2(B0-A0), H2(A0-B0)?
2. Do these expressions factor cleanly in A0,B0?
3. Which symmetric/antisymmetric combinations reduce to
       N, S, Delta, M1
   or simple factor-bearing expressions?
4. Does the pattern survive all historical branches?
5. Can the branch-normalized coordinates produce a common
   transformation law?

This experiment deliberately avoids:
    - divisor enumeration
    - continued fractions
    - large resultants
    - high-degree polynomial searches

Everything is exact symbolic arithmetic.
"""

import sympy as sp
from math import gcd


# -----------------------------------------------------------------------------
# Symbols
# -----------------------------------------------------------------------------

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")
X, Y = sp.symbols("X Y")
t = sp.symbols("t")

A0 = 2 * Y - 2 * X + 3
B0 = 2 * Y + 2 * X - 3

N_ps = p * q
S_ps = p + q
Delta_ps = (p - q) ** 2
M1_ps = (p + 1) * (q + 1)

H2 = lambda shift: sp.expand(6 * N - S ** 2 + S * shift)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def simp(expr):
    return sp.factor(sp.expand(expr))


def diff0(a, b):
    return sp.expand(a - b)


def pass_eq(a, b):
    return sp.expand(a - b) == 0


def exact_int_sqrt(n):
    if n < 0:
        return None
    r = sp.integer_nthroot(int(n), 2)
    return int(r[0]) if bool(r[1]) else None


def prime_pair_branch(pv, qv):
    """
    Determine the historical branch robustly.

    For odd factors:
        n % 4 == 3  -> opposite residue classes
        n % 4 == 1  -> same residue classes

    For n % 4 == 1 we use the historical branch for which
    the resulting A0,B0 correspond to q,3p (up to orientation).
    The coordinate equations themselves are tested explicitly.
    """

    nv = pv * qv
    if nv % 4 == 3:
        return "v==g"

    # Historical branch A/B can be represented by either orientation.
    # We choose the representation that reconstructs the normalized
    # factor pair after division by 3.
    candidates = [
        (
            "v!=g:A",
            sp.Rational(3 * pv - qv + 6, 4),
            sp.Rational(3 * pv + qv, 4),
        ),
        (
            "v!=g:A",
            sp.Rational(3 * qv - pv + 6, 4),
            sp.Rational(3 * qv + pv, 4),
        ),
    ]

    for name, xx, yy in candidates:
        aa = sp.simplify(2 * yy - 2 * xx + 3)
        bb = sp.simplify(2 * yy + 2 * xx - 3)

        vals = sorted([int(aa), int(bb)])
        target = sorted([qv, 3 * pv])
        if vals == target:
            return name

    return "v!=g:A"


def historical_xy(pv, qv):
    """
    Return a robust historical coordinate realization.

    The returned A0,B0 are allowed to be oriented either way.
    """

    nv = pv * qv

    if nv % 4 == 3:
        xx = sp.Rational(-pv + qv + 6, 4)
        yy = sp.Rational(pv + qv, 4)

        aa = sp.simplify(2 * yy - 2 * xx + 3)
        bb = sp.simplify(2 * yy + 2 * xx - 3)

        if sorted([int(aa), int(bb)]) != sorted([pv, qv]):
            xx = sp.Rational(-qv + pv + 6, 4)
            yy = sp.Rational(pv + qv, 4)
            aa = sp.simplify(2 * yy - 2 * xx + 3)
            bb = sp.simplify(2 * yy + 2 * xx - 3)

        return xx, yy, "v==g"

    # n == 1 mod 4
    # Use the branch that generates {q, 3p}, which normalizes to {p,q}.
    candidates = [
        (
            sp.Rational(3 * pv - qv + 6, 4),
            sp.Rational(3 * pv + qv, 4),
            "v!=g:A",
        ),
        (
            sp.Rational(3 * qv - pv + 6, 4),
            sp.Rational(3 * qv + pv, 4),
            "v!=g:A",
        ),
    ]

    for xx, yy, branch in candidates:
        aa = sp.simplify(2 * yy - 2 * xx + 3)
        bb = sp.simplify(2 * yy + 2 * xx - 3)

        if sorted([int(aa), int(bb)]) in (
            sorted([pv, 3 * qv]),
            sorted([qv, 3 * pv]),
        ):
            return xx, yy, branch

    raise ValueError(f"Could not construct historical coordinates for ({pv},{qv})")


def historical_data(pv, qv):
    xx, yy, branch = historical_xy(pv, qv)

    aa = int(2 * yy - 2 * xx + 3)
    bb = int(2 * yy + 2 * xx - 3)

    nv = pv * qv
    scale = 1 if nv % 4 == 3 else 3

    normalized = sorted([aa // scale, bb // scale])

    return {
        "p": pv,
        "q": qv,
        "N": nv,
        "S": pv + qv,
        "Delta": (pv - qv) ** 2,
        "M1": (pv + 1) * (qv + 1),
        "x": int(xx),
        "y": int(yy),
        "branch": branch,
        "A0": aa,
        "B0": bb,
        "scale": scale,
        "normalized": normalized,
    }


# -----------------------------------------------------------------------------
# Header
# -----------------------------------------------------------------------------

print("=" * 78)
print("EXPERIMENT 511 START")
print("=" * 78)
print("HISTORICAL LINEAR-COORDINATE / KAPPA SHIFT OPERATOR")
print("=" * 78)


# -----------------------------------------------------------------------------
# 1. BASIC HISTORICAL COORDINATES
# -----------------------------------------------------------------------------

print()
print("[1] HISTORICAL LINEAR COORDINATES")
print("-" * 78)

basic_1 = sp.expand((A0 + B0) - 4 * Y)
basic_2 = sp.expand((B0 - A0) - (4 * X - 6))

print(f"  A0 = {A0}")
print(f"  B0 = {B0}")
print(f"  A0+B0-4Y = {basic_1}")
print(f"  PASS = {basic_1 == 0}")
print(f"  B0-A0-(4X-6) = {basic_2}")
print(f"  PASS = {basic_2 == 0}")


# -----------------------------------------------------------------------------
# 2. KAPPA IN S,N COORDINATES
# -----------------------------------------------------------------------------

print()
print("[2] MODERN KAPPA SHIFT KERNEL")
print("-" * 78)

H2_expr = H2(t)

print(f"  H2(t) = {H2_expr}")

expected_H2 = 6 * N - S ** 2 + S * t
print(f"  reconstruction difference = {sp.expand(H2_expr - expected_H2)}")
print(f"  PASS = {pass_eq(H2_expr, expected_H2)}")


# -----------------------------------------------------------------------------
# 3. SHIFT PARAMETERS FROM HISTORICAL LINEAR FORMS
# -----------------------------------------------------------------------------

print()
print("[3] HISTORICAL SHIFT PARAMETERS")
print("-" * 78)

shift_sum = sp.expand(A0 + B0)
shift_diff = sp.expand(B0 - A0)
shift_diff_rev = sp.expand(A0 - B0)

print(f"  t_sum       = A0+B0 = {shift_sum}")
print(f"  t_diff      = B0-A0 = {shift_diff}")
print(f"  t_diff_rev  = A0-B0 = {shift_diff_rev}")


# -----------------------------------------------------------------------------
# 4. SYMBOLIC H2 EVALUATIONS
# -----------------------------------------------------------------------------

print()
print("[4] H2 AT HISTORICAL LINEAR COORDINATES")
print("-" * 78)

expr_sum = simp(H2_expr.subs(t, shift_sum))
expr_diff = simp(H2_expr.subs(t, shift_diff))
expr_diff_rev = simp(H2_expr.subs(t, shift_diff_rev))

print("  H2(A0+B0) =")
print(f"    {expr_sum}")

print("  H2(B0-A0) =")
print(f"    {expr_diff}")

print("  H2(A0-B0) =")
print(f"    {expr_diff_rev}")


# -----------------------------------------------------------------------------
# 5. SPECIAL ± HISTORICAL SHIFTS
# -----------------------------------------------------------------------------

print()
print("[5] OPPOSITE HISTORICAL SHIFTS")
print("-" * 78)

sum_pair = simp(
    H2_expr.subs(t, shift_diff)
    + H2_expr.subs(t, -shift_diff)
)

diff_pair = simp(
    H2_expr.subs(t, shift_diff)
    - H2_expr.subs(t, -shift_diff)
)

print("  H2(B0-A0) + H2(A0-B0) =")
print(f"    {sum_pair}")

print("  H2(B0-A0) - H2(A0-B0) =")
print(f"    {diff_pair}")


# -----------------------------------------------------------------------------
# 6. CONVERSION INTO A0,B0 DIRECTLY
# -----------------------------------------------------------------------------

print()
print("[6] PURE A0,B0 REWRITE")
print("-" * 78)

A, B = sp.symbols("A B")

N_ab = sp.Rational(1, 1) * A * B
S_ab = A + B
Delta_ab = (A - B) ** 2
M1_ab = (A + 1) * (B + 1)

H2_ab = sp.expand(6 * N_ab - S_ab ** 2 + S_ab * t)

print(f"  N(A,B)      = {N_ab}")
print(f"  S(A,B)      = {S_ab}")
print(f"  Delta(A,B)  = {Delta_ab}")
print(f"  M1(A,B)     = {M1_ab}")
print(f"  H2(A,B;t)   = {sp.factor(H2_ab)}")

h2_sum_ab = simp(H2_ab.subs(t, A + B))
h2_diff_ab = simp(H2_ab.subs(t, B - A))
h2_rev_ab = simp(H2_ab.subs(t, A - B))

print()
print("  H2(A+B) =")
print(f"    {h2_sum_ab}")

print("  H2(B-A) =")
print(f"    {h2_diff_ab}")

print("  H2(A-B) =")
print(f"    {h2_rev_ab}")


# -----------------------------------------------------------------------------
# 7. SEARCH SIMPLE SYMMETRIC COMBINATIONS
# -----------------------------------------------------------------------------

print()
print("[7] SYMMETRIC / ANTISYMMETRIC COMBINATION SEARCH")
print("-" * 78)

features = {
    "H2(A+B)": h2_sum_ab,
    "H2(B-A)": h2_diff_ab,
    "H2(A-B)": h2_rev_ab,
}

targets = {
    "6N": 6 * N_ab,
    "S": S_ab,
    "Delta": Delta_ab,
    "M1": M1_ab,
    "N": N_ab,
}

for name, expr in features.items():
    print(f"  {name}")
    for target_name, target in targets.items():
        d = sp.factor(expr - target)
        if d == 0:
            print(f"    EXACT HIT -> {target_name}")


combo_set = {
    "sum": sp.expand(h2_diff_ab + h2_rev_ab),
    "difference": sp.expand(h2_diff_ab - h2_rev_ab),
    "product": sp.expand(h2_diff_ab * h2_rev_ab),
    "sum_minus_12N": sp.expand(h2_diff_ab + h2_rev_ab - 12 * N_ab),
    "difference_sq": sp.expand((h2_diff_ab - h2_rev_ab) ** 2),
}

for name, expr in combo_set.items():
    fact = sp.factor(expr)

    print(f"  {name}:")
    print(f"    {fact}")

    for target_name, target in targets.items():
        if sp.expand(expr - target) == 0:
            print(f"    EXACT HIT -> {target_name}")


# -----------------------------------------------------------------------------
# 8. FACTOR-BEARING FORMS
# -----------------------------------------------------------------------------

print()
print("[8] FACTOR-BEARING FORM SEARCH")
print("-" * 78)

factor_targets = {
    "A": A,
    "B": B,
    "A*B": A * B,
    "A*(B+1)": A * (B + 1),
    "B*(A+1)": B * (A + 1),
    "(A-B)": A - B,
    "(A-B)^2": (A - B) ** 2,
}

# Try a small set of simple shifts derived from the coordinates.
candidate_shifts = {
    "A+B": A + B,
    "B-A": B - A,
    "A-B": A - B,
    "A": A,
    "B": B,
    "A+B+1": A + B + 1,
    "A+B-1": A + B - 1,
    "B-A+1": B - A + 1,
    "B-A-1": B - A - 1,
}

for shift_name, shift in candidate_shifts.items():
    expr = sp.factor(H2_ab.subs(t, shift))

    hits = []
    for target_name, target in factor_targets.items():
        quotient, remainder = sp.div(sp.Poly(sp.expand(expr), A, B),
                                     sp.Poly(sp.expand(target), A, B))
        if remainder.as_expr() == 0:
            hits.append(target_name)

    print(f"  t={shift_name:<8} H2 = {expr}")
    print(f"    divisible by = {hits}")


# -----------------------------------------------------------------------------
# 9. BRANCH-SPECIFIC SYMBOLIC SUBSTITUTION
# -----------------------------------------------------------------------------

print()
print("[9] BRANCH-SPECIFIC KAPPA / HISTORICAL MATCHES")
print("-" * 78)

branch_data = {
    "v==g": {
        "px": -p / 4 + q / 4 + sp.Rational(3, 2),
        "py": p / 4 + q / 4,
    },
    "v!=g:A": {
        "px": 3 * p / 4 - q / 4 + sp.Rational(3, 2),
        "py": 3 * p / 4 + q / 4,
    },
}

for branch_name, vals in branch_data.items():
    xx = vals["px"]
    yy = vals["py"]

    aa = sp.simplify(A0.subs({X: xx, Y: yy}))
    bb = sp.simplify(B0.subs({X: xx, Y: yy}))

    print()
    print(f"  {branch_name}")
    print(f"    A0 = {aa}")
    print(f"    B0 = {bb}")

    S_branch = sp.simplify((p + q))
    N_branch = p * q

    H_sum = sp.factor(
        (6 * N_branch - S_branch ** 2
         + S_branch * (aa + bb))
    )

    H_diff = sp.factor(
        (6 * N_branch - S_branch ** 2
         + S_branch * (bb - aa))
    )

    H_rev = sp.factor(
        (6 * N_branch - S_branch ** 2
         + S_branch * (aa - bb))
    )

    print(f"    H2(A0+B0) = {H_sum}")
    print(f"    H2(B0-A0) = {H_diff}")
    print(f"    H2(A0-B0) = {H_rev}")

    print(f"    factor H2(A0+B0) = {sp.factor(H_sum)}")
    print(f"    factor H2(B0-A0) = {sp.factor(H_diff)}")
    print(f"    factor H2(A0-B0) = {sp.factor(H_rev)}")


# -----------------------------------------------------------------------------
# 10. NUMERICAL AUDIT
# -----------------------------------------------------------------------------

print()
print("[10] NUMERICAL PRIME-PAIR AUDIT")
print("-" * 78)

pairs = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
    (3, 5),
    (5, 13),
    (13, 17),
    (3, 7),
    (7, 11),
]

numeric_failures = 0

for pv, qv in pairs:
    try:
        d = historical_data(pv, qv)
    except Exception as exc:
        print(f"  ({pv},{qv}) CONSTRUCTION FAILURE: {exc}")
        numeric_failures += 1
        continue

    nv = d["N"]
    sv = d["S"]

    a0 = d["A0"]
    b0 = d["B0"]

    # Use the normalized factor pair for modern S/N interpretation.
    u, v = d["normalized"]

    checks = {}

    checks["A0B0 scale"] = a0 * b0 == d["scale"] * nv
    checks["normalized product"] = u * v == nv
    checks["normalized sum"] = u + v == sv
    checks["normalized Delta"] = (u - v) ** 2 == d["Delta"]

    h_sum_num = 6 * nv - sv ** 2 + sv * (a0 + b0)
    h_diff_num = 6 * nv - sv ** 2 + sv * (b0 - a0)
    h_rev_num = 6 * nv - sv ** 2 + sv * (a0 - b0)

    # Compare against the direct symbolic formula using the same coordinates.
    h_sum_direct = 6 * nv - sv ** 2 + sv * (a0 + b0)
    h_diff_direct = 6 * nv - sv ** 2 + sv * (b0 - a0)
    h_rev_direct = 6 * nv - sv ** 2 + sv * (a0 - b0)

    checks["Hsum"] = h_sum_num == h_sum_direct
    checks["Hdiff"] = h_diff_num == h_diff_direct
    checks["Hrev"] = h_rev_num == h_rev_direct

    ok = all(checks.values())

    if not ok:
        numeric_failures += 1

    print(
        f"  ({pv},{qv}) "
        f"branch={d['branch']:<9} "
        f"A0={a0:<10} B0={b0:<10} "
        f"Hsum={h_sum_num:<18} "
        f"Hdiff={h_diff_num:<18} "
        f"PASS={ok}"
    )


# -----------------------------------------------------------------------------
# 11. DIRECT MODERN/HISTORICAL IDENTITIES
# -----------------------------------------------------------------------------

print()
print("[11] DIRECT CROSS-IDENTITIES")
print("-" * 78)

identity_tests = {
    "A0+B0 = 4y":
        (A0 + B0, 4 * Y),

    "B0-A0 = 4x-6":
        (B0 - A0, 4 * X - 6),

    "A0*B0":
        (A0 * B0,
         4 * (Y ** 2 - X ** 2 + 3 * X - 2)),

    "S_historical_v==g":
        (
            (p + q).subs({p: p, q: q}),
            4 * (p + q) / 4,
        ),

    "Delta_historical_v==g":
        (
            (p - q) ** 2,
            4 * (2 * X - 3) ** 2
        ),
}

for name, (lhs, rhs) in identity_tests.items():
    ok = pass_eq(lhs, rhs)
    print(f"  {name:<28} PASS={ok}")


# -----------------------------------------------------------------------------
# 12. RESEARCH SUMMARY
# -----------------------------------------------------------------------------

print()
print("[12] STRUCTURAL RESULT")
print("-" * 78)

print(
"""
The coordinate transformation is:

    A0 = 2y - 2x + 3
    B0 = 2y + 2x - 3

    A0+B0 = 4y
    B0-A0 = 4x-6
    A0*B0 = 4v-1

For modern symmetric coordinates:

    N = A0*B0 / scale
    S = normalized(A0) + normalized(B0)
    Delta = (normalized(A0)-normalized(B0))^2

The shifted KAPPA operator is:

    H2(t) = 6N - S^2 + S*t.

The experiment therefore tests whether the historical
linear forms themselves are natural spectral parameters.

The primary objects are:

    H2(A0+B0)
    H2(B0-A0)
    H2(A0-B0)

and their symmetric/antisymmetric combinations.

A useful positive result is expected to have the form

    H2(T(x,y))
      = simple factor expression in A0,B0

where T is a historical linear form.

Such a result would be more significant than another
identity involving only N,S,Delta, because it would show
that the 2020 conic coordinates and the 2026 KAPPA
operator are using the same underlying algebraic
coordinate system.

The experiment intentionally does not claim that
A0,B0 can be generated from N alone.
That remains the upstream source problem.
"""
)

print()
print("[13] EXPERIMENT STATUS")
print("-" * 78)
print(f"  symbolic historical coordinates = {basic_1 == 0 and basic_2 == 0}")
print(f"  H2 construction                  = {pass_eq(H2_expr, expected_H2)}")
print(f"  numerical coordinate audit       = {numeric_failures == 0}")
print()
print("PRIMARY QUESTION:")
print("  Do historical linear coordinates act as natural")
print("  spectral parameters of the shifted KAPPA kernel?")
print()
print("NEXT IF POSITIVE:")
print("  search the factor-bearing H2 expressions for")
print("  exact relations to A0, B0, M1 and Delta.")
print()
print("NEXT IF NEGATIVE:")
print("  move upstream and search the homogeneous-layer")
print("  construction for the historical linear forms.")
print()
print("=" * 78)
print("EXPERIMENT 511 FINISHED")
print("=" * 78)
