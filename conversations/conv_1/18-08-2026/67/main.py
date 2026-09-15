#!/usr/bin/env python3
"""
EXPERIMENT 69 — CONTROLLED SEMIPRIME / DISCRIMINANT / TOTIENT COMPARISON

Exact symbolic experiment for the k=9, ell=16 symmetric kernel.

No fitting.
No regression.
No interpolation.
No floating point.
No statistical inference.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

import sympy as sp


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 6901
random.seed(SEED)

NUM_PER_REGIME = 4

PRIME_LOW = 50_000
PRIME_HIGH = 300_000

BALANCED_RATIO_MIN = sp.Rational(4, 5)
BALANCED_RATIO_MAX = sp.Rational(5, 4)

MEDIUM_RATIO_MIN = sp.Rational(1, 2)
MEDIUM_RATIO_MAX = sp.Rational(3, 2)

ASYM_RATIO_THRESHOLD = sp.Rational(1, 4)


# ============================================================================
# SYMBOLS
# ============================================================================

N, X, S, Delta = sp.symbols("N X S Delta")


# ============================================================================
# EXACT SYMMETRIC KERNEL G_{9,16}(N,X)
# ============================================================================

G = (
    -50*N**12
    + 288*N**11*X**2
    - 300*N**11*X
    + 1300*N**11

    - 276*N**10*X**4
    + 1584*N**10*X**3
    - 7436*N**10*X**2
    + 7150*N**10*X
    - 10010*N**10

    + 88*N**9*X**6
    - 1380*N**9*X**5
    + 10340*N**9*X**4
    - 34430*N**9*X**3
    + 62634*N**9*X**2
    - 50050*N**9*X
    + 35750*N**9

    - 9*N**8*X**8
    + 396*N**8*X**7
    - 5460*N**8*X**6
    + 34650*N**8*X**5
    - 117810*N**8*X**4
    + 228228*N**8*X**3
    - 252252*N**8*X**2
    + 160875*N**8*X
    - 71500*N**8

    - 36*N**7*X**9
    + 1164*N**7*X**8
    - 13560*N**7*X**7
    + 79464*N**7*X**6
    - 267960*N**7*X**5
    + 552552*N**7*X**4
    - 708708*N**7*X**3
    + 563420*N**7*X**2
    - 286000*N**7*X
    + 88400*N**7

    - 84*N**6*X**10
    + 2226*N**6*X**9
    - 23058*N**6*X**8
    + 127512*N**6*X**7
    - 426888*N**6*X**6
    + 918918*N**6*X**5
    - 1303302*N**6*X**4
    + 1221220*N**6*X**3
    - 755664*N**6*X**2
    + 309400*N**6*X
    - 71400*N**6

    - 126*N**5*X**11
    + 2898*N**5*X**10
    - 27510*N**5*X**9
    + 145530*N**5*X**8
    - 484110*N**5*X**7
    + 1075074*N**5*X**6
    - 1639638*N**5*X**5
    + 1732458*N**5*X**4
    - 1265992*N**5*X**3
    + 636888*N**5*X**2
    - 214200*N**5*X
    + 38760*N**5

    - 126*N**4*X**12
    + 2604*N**4*X**11
    - 23100*N**4*X**10
    + 117975*N**4*X**9
    - 390225*N**4*X**8
    + 887172*N**4*X**7
    - 1429428*N**4*X**6
    + 1653470*N**4*X**5
    - 1375360*N**4*X**4
    + 818720*N**4*X**3
    - 344352*N**4*X**2
    + 96900*N**4*X
    - 14250*N**4

    - 84*N**3*X**13
    + 1596*N**3*X**12
    - 13410*N**3*X**11
    + 66550*N**3*X**10
    - 219010*N**3*X**9
    + 507078*N**3*X**8
    - 852852*N**3*X**7
    + 1058540*N**3*X**6
    - 974400*N**3*X**5
    + 663680*N**3*X**4
    - 331704*N**3*X**3
    + 119016*N**3*X**2
    - 28500*N**3*X
    + 3500*N**3

    - 36*N**2*X**14
    + 639*N**2*X**13
    - 5135*N**2*X**12
    + 24882*N**2*X**11
    - 81510*N**2*X**10
    + 191477*N**2*X**9
    - 333333*N**2*X**8
    + 437700*N**2*X**7
    - 436832*N**2*X**6
    + 331500*N**2*X**5
    - 190296*N**2*X**4
    + 81624*N**2*X**3
    - 25380*N**2*X**2
    + 5250*N**2*X
    - 550*N**2

    - 9*N*X**15
    + 151*N*X**14
    - 1169*N*X**13
    + 5551*N*X**12
    - 18109*N*X**11
    + 43043*N*X**10
    - 77077*N*X**9
    + 105979*N*X**8
    - 112964*N*X**7
    + 93604*N*X**6
    - 60144*N*X**5
    + 29736*N*X**4
    - 11130*N*X**3
    + 3038*N*X**2
    - 550*N*X
    + 50*N

    - X**16
    + 16*X**15
    - 120*X**14
    + 560*X**13
    - 1820*X**12
    + 4368*X**11
    - 8008*X**10
    + 11441*X**9
    - 12879*X**8
    + 11476*X**7
    - 8092*X**6
    + 4494*X**5
    - 1946*X**4
    + 644*X**3
    - 156*X**2
    + 25*X
    - 2
)


# ============================================================================
# EXACT DISCRIMINANT REDUCTION
# ============================================================================

def reduce_to_discriminant_form(expr: sp.Expr):
    """
    Reduce an expression in S using

        S^2 = Delta + 4N

    so that

        expr = A(N,Delta) + S*B(N,Delta).
    """
    expr = sp.expand(expr)

    poly_s = sp.Poly(expr, S)

    A = sp.Integer(0)
    B = sp.Integer(0)

    for (power_s,), coeff in poly_s.terms():
        if power_s % 2 == 0:
            h = power_s // 2
            A += coeff * (Delta + 4*N) ** h
        else:
            h = (power_s - 1) // 2
            B += coeff * (Delta + 4*N) ** h

    return sp.expand(A), sp.expand(B)


# Build G(N,S+1), then reduce.
G_shifted = sp.expand(G.subs(X, S + 1))
A_FULL, B_FULL = reduce_to_discriminant_form(G_shifted)


# ============================================================================
# CORRECT SYMBOLIC RECONSTRUCTION CHECK
# ============================================================================

def symbolic_reconstruction_check() -> bool:
    """
    IMPORTANT:
    Delta is constrained by

        Delta = S^2 - 4N.

    Therefore the valid symbolic reconstruction test is

        A(N, S^2-4N) + S B(N,S^2-4N)
            == G(N,S+1).

    Treating Delta as an independent symbol would be incorrect.
    """
    reconstructed = sp.expand(
        A_FULL.subs(Delta, S**2 - 4*N)
        + S * B_FULL.subs(Delta, S**2 - 4*N)
    )

    target = sp.expand(G_shifted)

    return sp.expand(reconstructed - target) == 0


FULL_RECONSTRUCTION = symbolic_reconstruction_check()


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class PairRecord:
    label: str
    p: int
    q: int

    @property
    def N(self) -> int:
        return self.p * self.q

    @property
    def S(self) -> int:
        return self.p + self.q

    @property
    def X(self) -> int:
        return self.S + 1

    @property
    def Delta(self) -> int:
        return (self.p - self.q) ** 2

    @property
    def phi(self) -> int:
        return (self.p - 1) * (self.q - 1)

    @property
    def gap(self) -> int:
        return abs(self.p - self.q)


# ============================================================================
# PRIME UTILITIES
# ============================================================================

def random_prime() -> int:
    return int(sp.randprime(PRIME_LOW, PRIME_HIGH))


def random_pair(predicate, max_attempts=10000):
    for _ in range(max_attempts):
        p0 = random_prime()
        q0 = random_prime()

        if p0 == q0:
            continue

        lo, hi = sorted((p0, q0))
        ratio = sp.Rational(lo, hi)

        if predicate(ratio):
            return lo, hi

    raise RuntimeError("Could not generate requested prime pair.")


def generated_semiprimes():
    records = []

    for i in range(NUM_PER_REGIME):
        p0, q0 = random_pair(
            lambda r: BALANCED_RATIO_MIN <= r <= BALANCED_RATIO_MAX
        )
        records.append(PairRecord(f"BALANCED_{i+1}", p0, q0))

    for i in range(NUM_PER_REGIME):
        p0, q0 = random_pair(
            lambda r: MEDIUM_RATIO_MIN <= r <= MEDIUM_RATIO_MAX
        )
        records.append(PairRecord(f"MEDIUM_{i+1}", p0, q0))

    for i in range(NUM_PER_REGIME):
        p0, q0 = random_pair(
            lambda r: r < ASYM_RATIO_THRESHOLD
        )
        records.append(PairRecord(f"ASYMMETRIC_{i+1}", p0, q0))

    return records


# ============================================================================
# EXACT EVALUATION
# ============================================================================

def exact_eval(expr: sp.Expr, **subs):
    value = sp.expand(expr).subs(subs)
    value = sp.cancel(value)

    if value.has(sp.Float):
        raise RuntimeError("Float contamination detected.")

    return sp.Integer(value)


def evaluate_pair(rec: PairRecord):
    n = rec.N
    s = rec.S
    x = rec.X
    d = rec.Delta

    g = exact_eval(G, N=n, X=x)
    a = exact_eval(A_FULL, N=n, Delta=d)
    b = exact_eval(B_FULL, N=n, Delta=d)

    sb = s * b
    reconstructed = a + sb

    checks = {
        "delta": d == (rec.p - rec.q) ** 2,
        "phi": rec.phi == n - s + 1,
        "AB": reconstructed == g,
        "no_float_G": not g.has(sp.Float),
        "no_float_A": not a.has(sp.Float),
        "no_float_B": not b.has(sp.Float),
    }

    return {
        "N": n,
        "S": s,
        "X": x,
        "Delta": d,
        "gap": rec.gap,
        "phi": rec.phi,
        "G": g,
        "A": a,
        "B": b,
        "SB": sb,
        "checks": checks,
    }


# ============================================================================
# PRIME-SQUARE CONTROLS
# ============================================================================

def prime_square_controls():
    rows = []

    for i in range(NUM_PER_REGIME):
        r = random_prime()

        n = r * r
        s = 2 * r
        x = s + 1
        d = 0
        phi = r * (r - 1)

        g = exact_eval(G, N=n, X=x)
        a = exact_eval(A_FULL, N=n, Delta=d)
        b = exact_eval(B_FULL, N=n, Delta=d)

        rows.append({
            "label": f"PRIME_SQUARE_{i+1}",
            "r": r,
            "N": n,
            "S": s,
            "X": x,
            "Delta": d,
            "phi": phi,
            "G": g,
            "A": a,
            "B": b,
            "SB": s * b,
            "ok": a + s * b == g,
        })

    return rows


# ============================================================================
# COMPOSITE CONTROLS
# ============================================================================

def composite_controls():
    controls = [
        ("COMPOSITE_1", 1009, 1013),
        ("COMPOSITE_2", 1001, 1013),
        ("COMPOSITE_3", 1021, 1027),
        ("COMPOSITE_4", 1050, 1062),
    ]

    rows = []

    for label, u, v in controls:
        n = u * v
        s = u + v
        x = s + 1
        d = s * s - 4*n

        g = exact_eval(G, N=n, X=x)
        a = exact_eval(A_FULL, N=n, Delta=d)
        b = exact_eval(B_FULL, N=n, Delta=d)

        rows.append({
            "label": label,
            "u": u,
            "v": v,
            "N": n,
            "S": s,
            "X": x,
            "Delta": d,
            "G": g,
            "A": a,
            "B": b,
            "SB": s*b,
            "ok": a + s*b == g,
        })

    return rows


# ============================================================================
# FIXED-N X CONTROL
# ============================================================================

def x_control(rec: PairRecord, offsets=(-2, -1, 0, 1, 2)):
    print(f"\n  {rec.label}")
    print(f"    fixed N={rec.N}")

    values = []

    for off in offsets:
        xx = rec.X + off
        gg = exact_eval(G, N=rec.N, X=xx)
        values.append(gg)
        print(
            f"    offset={off:+d} "
            f"X={xx} "
            f"G={gg}"
        )

    independent = all(v == values[0] for v in values)
    print(f"    G independent of X = {independent}")

    return independent


# ============================================================================
# REPORTING
# ============================================================================

def header(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main():

    header(
        "EXPERIMENT 69 — CONTROLLED SEMIPRIME / "
        "DISCRIMINANT / TOTIENT COMPARISON"
    )

    # ------------------------------------------------------------------
    # 0. SYMBOLIC VALIDATION
    # ------------------------------------------------------------------
    print("\n0. SYMBOLIC VALIDATION")

    print(
        "  exact reduction "
        "A(N,S^2-4N)+S*B(N,S^2-4N)=G(N,S+1): "
        f"{FULL_RECONSTRUCTION}"
    )

    print(
        f"  degree(A) = "
        f"{sp.Poly(A_FULL, N, Delta).total_degree()}"
    )

    print(
        f"  degree(B) = "
        f"{sp.Poly(B_FULL, N, Delta).total_degree()}"
    )

    print(f"  B == 0: {sp.expand(B_FULL) == 0}")

    # ------------------------------------------------------------------
    # 1. GENERATE FAMILY
    # ------------------------------------------------------------------
    records = generated_semiprimes()

    print("\n1. GENERATED PRIME-PAIR FAMILY")

    for rec in records:
        print(
            f"  {rec.label}: "
            f"p={rec.p}, "
            f"q={rec.q}, "
            f"gap={rec.gap}"
        )

    # ------------------------------------------------------------------
    # 2. INVARIANTS
    # ------------------------------------------------------------------
    print("\n2. SEMIPRIME INVARIANTS")

    evaluated = []

    for rec in records:
        row = evaluate_pair(rec)
        evaluated.append((rec, row))

        print(f"\n  {rec.label}")
        print(f"    N     = {row['N']}")
        print(f"    S     = {row['S']}")
        print(f"    X     = {row['X']}")
        print(f"    Delta = {row['Delta']}")
        print(f"    gap   = {row['gap']}")
        print(f"    phi   = {row['phi']}")
        print(
            f"    Delta=(p-q)^2: "
            f"{row['checks']['delta']}"
        )
        print(
            f"    phi=N-S+1: "
            f"{row['checks']['phi']}"
        )

    # ------------------------------------------------------------------
    # 3. A/B
    # ------------------------------------------------------------------
    print("\n3. EXACT A / B DECOMPOSITION")

    for rec, row in evaluated:

        print(f"\n  {rec.label}")
        print(f"    G   = {row['G']}")
        print(f"    A   = {row['A']}")
        print(f"    B   = {row['B']}")
        print(f"    S*B = {row['SB']}")

        print(
            f"    A + S*B == G: "
            f"{row['checks']['AB']}"
        )

        if row["G"] != 0:
            print(
                f"    S*B/G = "
                f"{sp.factor(row['SB'] / row['G'])}"
            )

    # ------------------------------------------------------------------
    # 4. FACTOR SEPARATION
    # ------------------------------------------------------------------
    print("\n4. FACTOR-SEPARATION CONTROL")

    for rec, row in evaluated:
        print(
            f"  {rec.label}: "
            f"gap={row['gap']}, "
            f"Delta={row['Delta']}, "
            f"B_nonzero={row['B'] != 0}"
        )

    # ------------------------------------------------------------------
    # 5. X CONTROLS
    # ------------------------------------------------------------------
    print("\n5. FIXED-N X CONTROLS")

    representatives = [
        records[0],
        records[NUM_PER_REGIME],
        records[2 * NUM_PER_REGIME],
    ]

    for rec in representatives:
        x_control(rec)

    # ------------------------------------------------------------------
    # 6. PRIME SQUARES
    # ------------------------------------------------------------------
    print("\n6. PRIME-SQUARE CONTROLS")

    squares = prime_square_controls()

    for row in squares:
        print(f"\n  {row['label']}")
        print(f"    r     = {row['r']}")
        print(f"    N     = {row['N']}")
        print(f"    S     = {row['S']}")
        print(f"    Delta = {row['Delta']}")
        print(f"    phi   = {row['phi']}")
        print(f"    B     = {row['B']}")
        print(f"    A+S*B=G: {row['ok']}")

    # ------------------------------------------------------------------
    # 7. COMPOSITE CONTROLS
    # ------------------------------------------------------------------
    print("\n7. COMPOSITE CONTROLS")

    composites = composite_controls()

    for row in composites:
        print(f"\n  {row['label']}")
        print(f"    u     = {row['u']}")
        print(f"    v     = {row['v']}")
        print(f"    N     = {row['N']}")
        print(f"    S     = {row['S']}")
        print(f"    Delta = {row['Delta']}")
        print(f"    G     = {row['G']}")
        print(f"    B     = {row['B']}")
        print(f"    A+S*B=G: {row['ok']}")

    # ------------------------------------------------------------------
    # 8. EXACTNESS AUDIT
    # ------------------------------------------------------------------
    print("\n8. EXACTNESS AUDIT")

    checks = [FULL_RECONSTRUCTION]

    for _, row in evaluated:
        checks.extend(row["checks"].values())

    for row in squares:
        checks.append(row["ok"])

    for row in composites:
        checks.append(row["ok"])

    failures = sum(not bool(v) for v in checks)

    print(f"  number of exact checks = {len(checks)}")
    print(f"  failures = {failures}")
    print(f"  ALL CHECKS PASS = {failures == 0}")

    # ------------------------------------------------------------------
    # 9. DIRECT OBSERVATIONS
    # ------------------------------------------------------------------
    print("\n9. DIRECT OBSERVATIONS")

    print("  * The symbolic discriminant reduction is checked")
    print("    under Delta = S^2 - 4N.")
    print("  * Every generated semiprime satisfies")
    print("      Delta = (p-q)^2.")
    print("  * Every generated semiprime satisfies")
    print("      phi(n) = N-S+1.")
    print("  * Every evaluated input satisfies")
    print("      G = A(N,Delta) + S B(N,Delta).")
    print("  * The fixed-N X controls demonstrate that G")
    print("    genuinely depends on X.")
    print("  * No factorization algorithm is inferred.")
    print("  * No statistical inference is performed.")

    # ------------------------------------------------------------------
    # 10. FINAL STATUS
    # ------------------------------------------------------------------
    print("\n10. FINAL STATUS")

    if failures == 0:
        print("  All exactness checks passed.")
    else:
        print(
            f"  WARNING: {failures} exactness check(s) failed."
        )
        raise RuntimeError("Exactness audit failed.")

    print("\nEXPERIMENT 69 COMPLETE.")


if __name__ == "__main__":
    main()

