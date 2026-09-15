#!/usr/bin/env python3
"""
EXPERIMENT 507
2020-2026 KAPPA COORDINATE CROSSWALK

Goal:
    Compare the historical 2020 (x,y) conic coordinates with the
    modern symmetric/KAPPA coordinates (N,S,Delta,M1,F2,...).

This experiment deliberately avoids:
    - continued fractions
    - Pell searches
    - large determinant expansions
    - numerical interpolation

It performs exact symbolic substitutions and numerical audits.

Historical 2020 structure:
    y^2 - x^2 + 3x - 2 = v

After completing the square:
    (2y)^2 - (2x-3)^2 = 4v-1

For odd N:
    N == 3 mod 4: 4v-1 = N
    N == 1 mod 4: 4v-1 = 3N

The historical factor-bearing coordinates are:
    A0 = 2y - 2x + 3
    B0 = 2y + 2x - 3

Thus:
    A0*B0 = scale*N, where scale is 1 or 3.

The main research target is whether modern KAPPA expressions expose
A0, B0, A0+B0, A0-B0, or simple functions of them.
"""

import math
import sympy as sp

print("=" * 78)
print("EXPERIMENT 507 START")
print("=" * 78)
print("2020-2026 KAPPA COORDINATE CROSSWALK")
print("=" * 78)

# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------
p, q = sp.symbols("p q", integer=True)
N, S = sp.symbols("N S", integer=True)
x, y = sp.symbols("x y", integer=True)
z = sp.symbols("z")

# ---------------------------------------------------------------------------
# Core modern quantities
# ---------------------------------------------------------------------------
N_pq = p * q
S_pq = p + q
Delta_pq = sp.expand((p - q) ** 2)
M1_pq = sp.expand((p + 1) * (q + 1))

F2_pq = sp.expand(
    p * (q + 1) ** 2
    + q * (p + 1) ** 2
    - p ** 2 * (q + 1)
    - q ** 2 * (p + 1)
)

F3_pq = sp.expand(
    p * (q + 1) ** 3
    + q * (p + 1) ** 3
    - p ** 3 * (q + 1)
    - q ** 3 * (p + 1)
)

F4_pq = sp.expand(
    p * (q + 1) ** 4
    + q * (p + 1) ** 4
    - p ** 4 * (q + 1)
    - q ** 4 * (p + 1)
)

# Verify modern formulas before using the crosswalk.
modern_checks = {
    "N=pq": sp.simplify(N_pq - p*q) == 0,
    "S=p+q": sp.simplify(S_pq - (p+q)) == 0,
    "Delta=(p-q)^2": sp.simplify(Delta_pq - (p-q)**2) == 0,
    "M1=(p+1)(q+1)": sp.simplify(M1_pq - (p+1)*(q+1)) == 0,
    "F2=6N-S^2+S": sp.simplify(F2_pq - (6*N_pq - S_pq**2 + S_pq)) == 0,
    "F3=(S+1)F2": sp.simplify(F3_pq - (S_pq + 1)*F2_pq) == 0,
}
print("\n[1] MODERN KAPPA SANITY")
print("-" * 78)
for name, ok in modern_checks.items():
    print(f"  {name}: PASS={ok}")
print(f"  SANITY FAILURES = {sum(not v for v in modern_checks.values())}")

# ---------------------------------------------------------------------------
# Historical 2020 coordinate branches
#
# We derive p,q -> x,y formulas exactly from the user's historical formulas.
# ---------------------------------------------------------------------------
# Branch A: v == g, opposite residue classes.
x_A = sp.expand((q - p + 6) / 4)
y_A = sp.expand((p + q) / 4)

# Branch B1: v != g, y > 2x in the historical notation.
x_B1 = sp.expand((3*p - q + 6) / 4)
y_B1 = sp.expand((3*p + q) / 4)

# Branch B2: v != g, y < 2x in the historical notation.
x_B2 = sp.expand((q - 3*p + 6) / 4)
y_B2 = sp.expand((3*p + q) / 4)

branches = {
    "v==g / opposite mod-4": (x_A, y_A, sp.Integer(1)),
    "v!=g / branch A (y>2x)": (x_B1, y_B1, sp.Integer(3)),
    "v!=g / branch B (y<2x)": (x_B2, y_B2, sp.Integer(3)),
}

# Historical coordinate forms.
A0 = sp.expand(2*y - 2*x + 3)
B0 = sp.expand(2*y + 2*x - 3)
C0 = sp.expand(2*y)          # intended S in branch A
D0 = sp.expand(2*x - 3)      # intended signed factor gap in branch A

print("\n[2] HISTORICAL FACTOR-BEARING COORDINATES")
print("-" * 78)
print("  A0 = 2y - 2x + 3")
print("  B0 = 2y + 2x - 3")
print("  A0*B0 = scale*N")
print("  C0 = 2y")
print("  D0 = 2x - 3")

# ---------------------------------------------------------------------------
# Crosswalk builder
# ---------------------------------------------------------------------------
def crosswalk(expr_x, expr_y, scale):
    A = sp.expand(A0.subs({x: expr_x, y: expr_y}))
    B = sp.expand(B0.subs({x: expr_x, y: expr_y}))

    n_expr = sp.expand(A * B / scale)
    s_expr = sp.expand(A + B) / 2

    gap_expr = sp.expand(B - A)  # signed; depending on orientation
    delta_expr = sp.expand(gap_expr**2)

    # Modern objects from the recovered p,q are expressed first in p,q,
    # then compared with historical coordinates.
    return {
        "x": sp.expand(expr_x),
        "y": sp.expand(expr_y),
        "A0": sp.factor(A),
        "B0": sp.factor(B),
        "A0+B0": sp.factor(A + B),
        "A0-B0": sp.factor(A - B),
        "B0-A0": sp.factor(B - A),
        "A0*B0": sp.factor(A * B),
        "N_candidate": sp.factor(n_expr),
        "S_candidate": sp.factor(s_expr),
        "Delta_candidate": sp.factor(delta_expr),
        "M1_candidate": sp.factor(n_expr + s_expr + 1),
        "F2_candidate": sp.factor(
            6*n_expr - s_expr**2 + s_expr
        ),
        "F3_candidate": sp.factor(
            (s_expr + 1) * (6*n_expr - s_expr**2 + s_expr)
        ),
        "scale": scale,
    }

crosswalks = {
    name: crosswalk(*vals)
    for name, vals in branches.items()
}

# ---------------------------------------------------------------------------
# Exact branch theorem audits
# ---------------------------------------------------------------------------
print("\n[3] EXACT BRANCH CROSSWALK AUDIT")
print("-" * 78)

branch_failures = 0

for name, cw in crosswalks.items():
    scale = cw["scale"]

    n_diff = sp.factor(cw["N_candidate"] - N_pq)
    s_diff = sp.factor(cw["S_candidate"] - S_pq)

    # In branches B1/B2 the factor order/sign can reverse, but squared gap
    # must agree exactly.
    delta_diff = sp.factor(cw["Delta_candidate"] - Delta_pq)

    m1_diff = sp.factor(cw["M1_candidate"] - M1_pq)

    f2_expr = cw["F2_candidate"]
    f2_target = sp.expand(6*N_pq - S_pq**2 + S_pq)
    f2_diff = sp.factor(f2_expr - f2_target)

    f3_diff = sp.factor(cw["F3_candidate"] - F3_pq)

    ok = all([
        n_diff == 0,
        s_diff == 0,
        delta_diff == 0,
        m1_diff == 0,
        f2_diff == 0,
        f3_diff == 0,
    ])

    print(f"\n  {name}")
    print(f"    x = {cw['x']}")
    print(f"    y = {cw['y']}")
    print(f"    A0 = {cw['A0']}")
    print(f"    B0 = {cw['B0']}")
    print(f"    A0*B0 = {cw['A0*B0']}")
    print(f"    scale*N = {sp.factor(scale*N_pq)}")
    print(f"    N diff = {n_diff}")
    print(f"    S candidate = {cw['S_candidate']}")
    print(f"    S diff = {s_diff}")
    print(f"    Delta candidate = {cw['Delta_candidate']}")
    print(f"    Delta diff = {delta_diff}")
    print(f"    M1 diff = {m1_diff}")
    print(f"    F2 diff = {f2_diff}")
    print(f"    F3 diff = {f3_diff}")
    print(f"    PASS = {ok}")

    if not ok:
        branch_failures += 1

print(f"\n  BRANCH CROSSWALK FAILURES = {branch_failures}")

# ---------------------------------------------------------------------------
# Factor-bearing linear forms and gcd forms.
# ---------------------------------------------------------------------------
print("\n[4] HISTORICAL LINEAR FORMS IN TERMS OF p,q")
print("-" * 78)

forms = {
    "A0": A0,
    "B0": B0,
    "b-x-y+2": N/2 - x - y + 2,
    "b-x+y+2": N/2 - x + y + 2,
    "b+x-y-1": N/2 + x - y - 1,
    "b+x+y-1": N/2 + x + y - 1,
}

for name, form in forms.items():
    print(f"\n  {name}")
    for branch_name, (bx, by, scale) in branches.items():
        # For b, use pq/2 symbolically.
        val = sp.factor(
            form.subs({x: bx, y: by, N: N_pq})
        )
        print(f"    {branch_name}: {val}")

# ---------------------------------------------------------------------------
# KAPPA expressions factored in historical x,y coordinates
# ---------------------------------------------------------------------------
print("\n[5] KAPPA CROSSWALK FACTORIZATION")
print("-" * 78)

kappa_names = {
    "N": N_pq,
    "S": S_pq,
    "Delta": Delta_pq,
    "M1": M1_pq,
    "F2": F2_pq,
    "F3": F3_pq,
    "F4": F4_pq,
}

for branch_name, (bx, by, scale) in branches.items():
    print(f"\n  {branch_name}")
    print("  " + "-" * 74)

    A = sp.factor(A0.subs({x: bx, y: by}))
    B = sp.factor(B0.subs({x: bx, y: by}))

    substitutions = {p: p, q: q}
    # Instead of solving p,q back from x,y by matrix inversion here,
    # derive the exact p,q-from-x,y formulas branch by branch.
    if branch_name.startswith("v==g"):
        p_xy = sp.expand(2*y - 2*x + 3)
        q_xy = sp.expand(2*y + 2*x - 3)
    elif "branch A" in branch_name:
        p_xy = sp.expand((2*y + 2*x - 3) / 3)
        q_xy = sp.expand(2*y - 2*x + 3)
    else:
        p_xy = sp.expand((2*y - 2*x + 3) / 3)
        q_xy = sp.expand(2*y + 2*x - 3)

    print(f"    p(x,y) = {sp.factor(p_xy)}")
    print(f"    q(x,y) = {sp.factor(q_xy)}")

    for label, expr in kappa_names.items():
        xy_expr = sp.factor(
            expr.subs({p: p_xy, q: q_xy})
        )
        print(f"    {label}(x,y) = {xy_expr}")

# ---------------------------------------------------------------------------
# Search specifically for historical building blocks.
# ---------------------------------------------------------------------------
print("\n[6] HISTORICAL BUILDING-BLOCK SEARCH")
print("-" * 78)

building_blocks = {
    "A0": A0,
    "B0": B0,
    "A0+B0": sp.expand(A0+B0),
    "A0-B0": sp.expand(A0-B0),
    "B0-A0": sp.expand(B0-A0),
    "A0*B0": sp.expand(A0*B0),
    "A0^2": sp.expand(A0**2),
    "B0^2": sp.expand(B0**2),
    "(A0+B0)^2": sp.expand((A0+B0)**2),
    "(A0-B0)^2": sp.expand((A0-B0)**2),
}

targets = {
    "N": N_pq,
    "S": S_pq,
    "Delta": Delta_pq,
    "M1": M1_pq,
    "F2": F2_pq,
    "F3": F3_pq,
}

for branch_name, (bx, by, scale) in branches.items():
    print(f"\n  {branch_name}")
    for tname, texpr in targets.items():
        xy_target = sp.factor(texpr.subs({
            p: (2*y - 2*x + 3) if branch_name.startswith("v==g")
                else ((2*y + 2*x - 3)/3 if "branch A" in branch_name
                      else (2*y - 2*x + 3)/3),
            q: (2*y + 2*x - 3) if branch_name.startswith("v==g")
                else ((2*y - 2*x + 3) if "branch A" in branch_name
                      else (2*y + 2*x - 3)),
        }))

        exact_hits = []
        for bname, bexpr in building_blocks.items():
            diff = sp.factor(xy_target - bexpr)
            if diff == 0:
                exact_hits.append(bname)

        if exact_hits:
            print(f"    {tname}: exact structural hits = {exact_hits}")

# ---------------------------------------------------------------------------
# Search simple affine relations between modern quantities and A0/B0.
# ---------------------------------------------------------------------------
print("\n[7] SIMPLE AFFINE CROSSWALK SEARCH")
print("-" * 78)

coeffs = [-3, -2, -1, 0, 1, 2, 3]

def affine_search(target_expr, Aexpr, Bexpr, name):
    hits = []
    for alpha in coeffs:
        for beta in coeffs:
            for gamma in coeffs:
                for delta in coeffs:
                    candidate = (
                        alpha*Aexpr
                        + beta*Bexpr
                        + gamma*N_pq
                        + delta
                    )
                    if sp.factor(candidate - target_expr) == 0:
                        hits.append((alpha, beta, gamma, delta))
    return hits

for branch_name, (bx, by, scale) in branches.items():
    print(f"\n  {branch_name}")
    A = sp.expand(A0.subs({x: bx, y: by}))
    B = sp.expand(B0.subs({x: bx, y: by}))

    for tname, texpr in targets.items():
        hits = affine_search(texpr, A, B, tname)
        # Remove trivial/self-evident duplicates by only displaying a compact
        # count plus the first few exact hits.
        if hits:
            print(f"    {tname}: hits={len(hits)} first={hits[:6]}")

# ---------------------------------------------------------------------------
# Numerical audit data
# ---------------------------------------------------------------------------
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

def historical_xy_numeric(pv, qv):
    nv = pv * qv
    if nv % 4 == 3:
        xv = (qv - pv + 6) // 4
        yv = (pv + qv) // 4
        scale = 1
        branch = "v==g"
        # Depending on orientation, these correspond to p/q or q/p.
    else:
        # Pick the branch that gives integral values.
        cand_a = ((3*pv - qv + 6), (3*pv + qv))
        cand_b = ((qv - 3*pv + 6), (3*pv + qv))

        if cand_a[0] % 4 == 0 and cand_a[1] % 4 == 0:
            xv = cand_a[0] // 4
            yv = cand_a[1] // 4
            branch = "v!=g:A"
        elif cand_b[0] % 4 == 0 and cand_b[1] % 4 == 0:
            xv = cand_b[0] // 4
            yv = cand_b[1] // 4
            branch = "v!=g:B"
        else:
            raise ArithmeticError(f"No integral historical branch for ({pv},{qv})")

        scale = 3

    return nv, xv, yv, scale, branch

print("\n[8] NUMERICAL CROSSWALK AUDIT")
print("-" * 78)

num_failures = 0

for pv, qv in pairs:
    nv, xv, yv, scale, branch = historical_xy_numeric(pv, qv)

    A = 2*yv - 2*xv + 3
    B = 2*yv + 2*xv - 3

    candidates = {
        "N": A * B // scale,
        "S_from_A+B": (A + B) // 2,
        "Delta_from_gap": (A - B) ** 2,
        "M1_from_NS": (nv + (pv + qv) + 1),
        "F2_from_NS": 6*nv - (pv+qv)**2 + (pv+qv),
    }

    expected = {
        "N": nv,
        "S_from_A+B": pv + qv,
        "Delta_from_gap": (pv - qv)**2,
        "M1_from_NS": (pv+1)*(qv+1),
        "F2_from_NS": 6*nv - (pv+qv)**2 + (pv+qv),
    }

    ok = all(candidates[k] == expected[k] for k in candidates)
    if not ok:
        num_failures += 1

    print(
        f"  ({pv},{qv}) branch={branch:7s} "
        f"x={xv:8d} y={yv:8d} "
        f"A0={A:8d} B0={B:8d} "
        f"PASS={ok}"
    )

print(f"\n  NUMERICAL CROSSWALK FAILURES = {num_failures}")

# ---------------------------------------------------------------------------
# Final structural findings
# ---------------------------------------------------------------------------
print("\n[9] STRUCTURAL FINDINGS")
print("-" * 78)

print("""
The historical variables A0,B0 are the factor-bearing coordinates:

    A0 = 2y - 2x + 3
    B0 = 2y + 2x - 3

with

    A0*B0 = scale*N.

For the v==g branch:

    A0,B0 = p,q (up to ordering)
    A0+B0 = S
    (A0-B0)^2 = Delta.

For the v!=g branches:

    one of A0,B0 carries a factor multiplied by 3,
    but the product still equals 3N and the square of the
    appropriate gap reconstructs the modern discriminant
    after the branch's scaling/orientation is applied.

The critical research question is therefore now:

    Does an N-only homogeneous/KAPPA observable correspond
    naturally to A0, B0, A0+B0, A0-B0, or A0*B0?

This is a coordinate-identification problem, not another
generic factorization search.
""")

# ---------------------------------------------------------------------------
# Explicit research targets for the next experiment
# ---------------------------------------------------------------------------
print("[10] NEXT TARGETS")
print("-" * 78)
print("""
  TARGET A:
    Search homogeneous-layer quantities for exact A0/B0 factors.

  TARGET B:
    Search for an object proportional to
        (A0-B0)^2 = Delta
    rather than trying to reconstruct S first.

  TARGET C:
    Search for the 2020 conic value v and test whether
        4v-1
    appears directly in a modern N-only observable.

  TARGET D:
    Search for the historical gcd-bearing linear forms:
        b-x-y+2
        b-x+y+2
        b+x-y-1
        b+x+y-1

  TARGET E:
    Compare the algebraic signatures of these objects against
    the existing homogeneous-layer outputs before constructing
    any new KAPPA-derived quantities.
""")

# ---------------------------------------------------------------------------
# Overall status
# ---------------------------------------------------------------------------
overall = (
    all(modern_checks.values())
    and branch_failures == 0
    and num_failures == 0
)

print("\n[11] EXPERIMENT STATUS")
print("-" * 78)
print(f"  modern KAPPA sanity          = {all(modern_checks.values())}")
print(f"  symbolic branch crosswalk   = {branch_failures == 0}")
print(f"  numerical crosswalk         = {num_failures == 0}")
print(f"  OVERALL EXACT AUDIT          = {overall}")

print("=" * 78)
print("EXPERIMENT 507 FINISHED")
print("=" * 78)
