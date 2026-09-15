#!/usr/bin/env python3

import sympy as sp


print("=" * 78)
print("EXPERIMENT 510 START")
print("=" * 78)
print("2020 CONIC -> KAPPA ALGEBRAIC CROSS-BRIDGE SEARCH")
print("=" * 78)


# =============================================================================
# SYMBOLS
# =============================================================================

p, q = sp.symbols("p q")
x, y = sp.symbols("x y")
N, S = sp.symbols("N S")

A0 = 2 * y - 2 * x + 3
B0 = 2 * y + 2 * x - 3

v = y**2 - x**2 + 3*x - 2


# =============================================================================
# BRANCHES
# =============================================================================

branches = {
    "v==g": {
        "p": -2*x + 2*y + 3,
        "q": 2*x + 2*y - 3,
        "scale": 1,
    },
    "v!=g:A": {
        "p": (2*x + 2*y - 3) / 3,
        "q": -2*x + 2*y + 3,
        "scale": 3,
    },
    "v!=g:B": {
        "p": (-2*x + 2*y + 3) / 3,
        "q": 2*x + 2*y - 3,
        "scale": 3,
    },
}


# =============================================================================
# HELPERS
# =============================================================================

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def half_floor(n):
    return (n - 1) // 2


def compute_vg(n):
    b = half_floor(n)
    vv = (b * b) % n
    g = b // 2 + 1
    return b, vv, g


def historical_candidates(pv, qv):
    return [
        (
            "v==g",
            sp.Rational(qv - pv + 6, 4),
            sp.Rational(pv + qv, 4),
        ),
        (
            "v!=g:A",
            sp.Rational(3 * pv - qv + 6, 4),
            sp.Rational(3 * pv + qv, 4),
        ),
        (
            "v!=g:B",
            sp.Rational(qv - 3 * pv + 6, 4),
            sp.Rational(3 * pv + qv, 4),
        ),
    ]


def is_integer_rational(z):
    return isinstance(z, sp.Rational) and z.q == 1


# =============================================================================
# DATASET
# =============================================================================

TEST_CASES = [
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


# =============================================================================
# [1] BASIC HISTORICAL IDENTITIES
# =============================================================================

print("\n[1] BASIC 2020 / 2026 COORDINATE SYSTEM")
print("-" * 78)

basic_checks = {
    "A0+B0 = 4y": A0 + B0 - 4*y,
    "B0-A0 = 4x-6": B0 - A0 - (4*x - 6),
    "A0*B0 = 4v-1": A0 * B0 - (4*v - 1),
}

basic_failures = 0

for name, expr in basic_checks.items():
    passed = sp.expand(expr) == 0
    print(f"  {name}: PASS={passed}")
    if not passed:
        print(f"    difference = {sp.factor(expr)}")
        basic_failures += 1

print(f"  BASIC FAILURES = {basic_failures}")


# =============================================================================
# [2] MODERN KAPPA QUANTITIES IN HISTORICAL COORDINATES
# =============================================================================

print("\n[2] MODERN KAPPA QUANTITIES IN (x,y)")
print("-" * 78)

for branch_name, info in branches.items():

    pp = info["p"]
    qq = info["q"]

    NN = simp(pp * qq)
    SS = simp(pp + qq)
    Delta = simp((pp - qq)**2)
    M1 = simp((pp + 1) * (qq + 1))

    F2 = simp(6 * NN - SS**2 + SS)
    F3 = simp((SS + 1) * F2)

    print(f"\n  {branch_name}")
    print(f"    N     = {NN}")
    print(f"    S     = {SS}")
    print(f"    Delta = {Delta}")
    print(f"    M1    = {M1}")
    print(f"    F2    = {F2}")
    print(f"    F3    = {F3}")


# =============================================================================
# [3] CONIC CERTIFICATE
# =============================================================================

print("\n[3] HISTORICAL CONIC CERTIFICATE")
print("-" * 78)

conic_difference = simp(
    y**2 - x**2 + 3*x - 2 - v
)

print(f"  difference = {conic_difference}")
print(f"  PASS = {conic_difference == 0}")


# =============================================================================
# [4] BRANCH RELATIONS
# =============================================================================

print("\n[4] v / g / N BRANCH RELATIONS")
print("-" * 78)

branch_failures = 0

for branch_name, info in branches.items():

    pp = info["p"]
    qq = info["q"]
    scale = info["scale"]

    NN = simp(pp * qq)
    K = simp(4*v - 1)

    if branch_name == "v==g":
        checks = [
            simp(K - NN),
            simp(v - (NN + 1) / 4),
        ]
    else:
        checks = [
            simp(K - 3*NN),
            simp(v - (3*NN + 1) / 4),
        ]

    passed = all(expr == 0 for expr in checks)

    print(f"  {branch_name}: PASS={passed}")

    if not passed:
        for expr in checks:
            if expr != 0:
                print(f"    difference = {expr}")
        branch_failures += 1

print(f"  BRANCH FAILURES = {branch_failures}")


# =============================================================================
# [5] KAPPA -> DISCRIMINANT
# =============================================================================

print("\n[5] KAPPA F2 -> DISCRIMINANT BRIDGE")
print("-" * 78)

F2_symbolic = 6*N - S**2 + S
Delta_symbolic = S**2 - 4*N

delta_identity = simp(
    (2*N - F2_symbolic + S) - Delta_symbolic
)

print("  Delta = 2N - F2 + S")
print(f"  difference = {delta_identity}")
print(f"  PASS = {delta_identity == 0}")


# =============================================================================
# [6] F2/F3 -> DISCRIMINANT
# =============================================================================

print("\n[6] F2/F3 -> DISCRIMINANT")
print("-" * 78)

F3_symbolic = (S + 1) * F2_symbolic

delta_from_f = simp(
    2*N
    - F2_symbolic
    + F3_symbolic / F2_symbolic
    - 1
)

delta_from_f_difference = simp(
    delta_from_f - Delta_symbolic
)

print(f"  reconstructed Delta = {delta_from_f}")
print(f"  target Delta        = {Delta_symbolic}")
print(f"  difference           = {delta_from_f_difference}")
print(f"  PASS = {delta_from_f_difference == 0}")


# =============================================================================
# [7] SMALL LINEAR COMBINATION SEARCH
# =============================================================================

print("\n[7] SMALL LINEAR COMBINATION SEARCH")
print("-" * 78)

targets = {
    "x": x,
    "y": y,
    "2x-3": 2*x - 3,
    "2y": 2*y,
    "A0": A0,
    "B0": B0,
    "A0+B0": A0 + B0,
    "B0-A0": B0 - A0,
    "(B0-A0)^2": (B0 - A0)**2,
    "A0*B0": A0 * B0,
}

a, b, c, d, e = sp.symbols("a b c d e")

for branch_name, info in branches.items():

    pp = info["p"]
    qq = info["q"]

    NN = simp(pp * qq)

    if branch_name == "v==g":
        gg = v
    else:
        gg = simp((NN + 3) / 4)

    FF2 = simp(6*NN - (pp + qq)**2 + (pp + qq))

    basis = [FF2, NN, v, gg, sp.Integer(1)]

    print(f"\n  {branch_name}")

    total_hits = 0

    for target_name, target in targets.items():

        expr = sp.Poly(
            sp.expand(
                a*basis[0]
                + b*basis[1]
                + c*basis[2]
                + d*basis[3]
                + e
                - target
            ),
            x,
            y,
        )

        equations = [
            coeff
            for coeff in expr.coeffs()
        ]

        if not equations:
            continue

        sol = sp.solve(
            equations,
            (a, b, c, d, e),
            dict=True
        )

        if sol:
            valid = []

            for candidate in sol:
                # Reject arbitrary-parameter families.
                unresolved = any(
                    symbol not in candidate
                    for symbol in (a, b, c, d, e)
                )

                if not unresolved:
                    valid.append(candidate)

            if valid:
                total_hits += len(valid)
                print(f"    {target_name}: {valid[:3]}")

    print(f"    total exact target hits = {total_hits}")


# =============================================================================
# [8] GAP-SQUARE SEARCH
# =============================================================================

print("\n[8] GAP-SQUARE SIGNATURE SEARCH")
print("-" * 78)

for branch_name, info in branches.items():

    pp = info["p"]
    qq = info["q"]

    NN = simp(pp * qq)
    SS = simp(pp + qq)
    FF2 = simp(6*NN - SS**2 + SS)

    Delta = simp((pp - qq)**2)
    historical_gap = simp((B0 - A0)**2)

    print(f"\n  {branch_name}")
    print(f"    2N-F2+S       = {simp(2*NN - FF2 + SS)}")
    print(f"    Delta         = {Delta}")
    print(f"    (B0-A0)^2     = {historical_gap}")
    print(
        f"    Delta gap     = "
        f"{simp(Delta - historical_gap)}"
    )


# =============================================================================
# [9] SHIFTED H2 SEARCH
# =============================================================================

print("\n[9] SHIFTED H2 PARAMETER SEARCH")
print("-" * 78)

t_candidates = {
    "0": sp.Integer(0),
    "1": sp.Integer(1),
    "-1": sp.Integer(-1),
    "x": x,
    "-x": -x,
    "2x-3": 2*x - 3,
    "3-2x": 3 - 2*x,
    "4x-6": 4*x - 6,
    "6-4x": 6 - 4*x,
    "2y": 2*y,
    "-2y": -2*y,
    "4y": 4*y,
    "-4y": -4*y,
}

for branch_name, info in branches.items():

    pp = info["p"]
    qq = info["q"]

    NN = simp(pp * qq)
    SS = simp(pp + qq)

    print(f"\n  {branch_name}")

    for tname, tv in t_candidates.items():

        H = simp(6*NN - SS**2 + SS*tv)

        numerator, factors = sp.factor_list(
            sp.expand(H)
        )

        factor_names = []

        for factor, exponent in factors:
            if simp(factor - A0) == 0:
                factor_names.append("A0")
            elif simp(factor - B0) == 0:
                factor_names.append("B0")
            elif simp(factor - (B0-A0)) == 0:
                factor_names.append("B0-A0")

        interesting = (
            len(factors) > 1
            or bool(factor_names)
        )

        if interesting:
            print(
                f"    t={tname:<8} "
                f"H2(t)={sp.factor(H)}"
            )


# =============================================================================
# [10] HISTORICAL FACTOR DIVISIBILITY SEARCH
# =============================================================================

print("\n[10] HISTORICAL FACTOR DIVISIBILITY SEARCH")
print("-" * 78)

historical_forms = {
    "A0": A0,
    "B0": B0,
    "A0-B0": A0 - B0,
    "B0-A0": B0 - A0,
    "A0*B0": A0 * B0,
}

for branch_name, info in branches.items():

    pp = info["p"]
    qq = info["q"]

    NN = simp(pp * qq)
    SS = simp(pp + qq)

    print(f"\n  {branch_name}")

    for tname, tv in t_candidates.items():

        H = sp.expand(
            6*NN - SS**2 + SS*tv
        )

        hits = []

        for form_name, form in historical_forms.items():

            divisor_poly = sp.Poly(
                sp.expand(form),
                x,
                y
            )

            if divisor_poly.is_zero:
                continue

            quotient, remainder = sp.div(
                sp.Poly(H, x, y),
                divisor_poly
            )

            if remainder.as_expr() == 0:
                hits.append(form_name)

        if hits:
            print(
                f"    t={tname:<8} "
                f"divisible by {hits}"
            )


# =============================================================================
# [11] NUMERICAL CROSSWALK
# =============================================================================

print("\n[11] NUMERICAL CROSS-CHECK")
print("-" * 78)

numeric_failures = 0

for pv, qv in TEST_CASES:

    nn = pv * qv

    selected = None

    for branch_name, xx, yy in historical_candidates(
        pv, qv
    ):

        if not (
            is_integer_rational(xx)
            and is_integer_rational(yy)
        ):
            continue

        xx_i = int(xx)
        yy_i = int(yy)

        aa = 2*yy_i - 2*xx_i + 3
        bb0 = 2*yy_i + 2*xx_i - 3

        product = aa * bb0

        if product not in (nn, 3*nn):
            continue

        scale = 1 if product == nn else 3

        # FIX:
        # only construct normalized candidates when they are
        # actually integral; never insert None into sorted().
        normalized_values = []

        if aa % scale == 0:
            normalized_values.append(aa // scale)

        if bb0 % scale == 0:
            normalized_values.append(bb0 // scale)

        if len(normalized_values) != 2:
            continue

        normalized = sorted(normalized_values)

        if normalized == sorted([pv, qv]):
            selected = (
                branch_name,
                xx_i,
                yy_i,
                aa,
                bb0,
                scale,
            )
            break

    if selected is None:
        print(f"  ({pv},{qv}) SELECTION FAILURE")
        numeric_failures += 1
        continue

    branch_name, xx_i, yy_i, aa, bb0, scale = selected

    normalized_a = aa // scale
    normalized_b = bb0 // scale

    S_true = pv + qv
    Delta_true = (pv - qv)**2
    M1_true = (pv + 1) * (qv + 1)

    S_hist = normalized_a + normalized_b
    Delta_hist = (normalized_a - normalized_b)**2
    M1_hist = (normalized_a + 1) * (normalized_b + 1)

    passed = (
        S_hist == S_true
        and Delta_hist == Delta_true
        and M1_hist == M1_true
    )

    print(
        f"  ({pv},{qv}) "
        f"branch={branch_name:<10} "
        f"x={xx_i:<8} "
        f"y={yy_i:<8} "
        f"S={S_hist:<10} "
        f"Delta={Delta_hist:<14} "
        f"PASS={passed}"
    )

    if not passed:
        numeric_failures += 1

print(f"\n  NUMERICAL FAILURES = {numeric_failures}")


# =============================================================================
# [12] CENTRAL CROSS-BRIDGE
# =============================================================================

print("\n[12] CENTRAL CROSS-BRIDGE")
print("-" * 78)

print("  A0+B0       =", sp.factor(A0 + B0))
print("  B0-A0       =", sp.factor(B0 - A0))
print("  A0*B0       =", sp.factor(A0 * B0))
print("  4v-1        =", sp.factor(4*v - 1))
print("  (B0-A0)^2   =", sp.factor((B0-A0)**2))

print()
print("  Modern coordinates:")
print("    N     = pq")
print("    S     = p+q")
print("    Delta = (p-q)^2")
print("    M1    = (p+1)(q+1)")
print()
print("  Historical coordinates:")
print("    A0 = 2y-2x+3")
print("    B0 = 2y+2x-3")
print("    A0+B0 = S       after normalization")
print("    (A0-B0)^2 = Delta after normalization")


# =============================================================================
# [13] RESEARCH INTERPRETATION
# =============================================================================

print("\n[13] RESEARCH INTERPRETATION")
print("-" * 78)

print(
"""
The experiment treats the 2020 variables as a genuine
coordinate system rather than as a factor-search procedure.

Exact identities:

    A0+B0 = 4y
    B0-A0 = 4x-6
    A0*B0 = 4v-1

After branch normalization:

    A0,B0 -> factor pair
    A0+B0 -> S
    (A0-B0)^2 -> Delta
    (A0+1)(B0+1) -> M1

The modern KAPPA sequence also satisfies:

    F2 = 6N-S^2+S

and therefore:

    Delta = 2N-F2+S.

The experiment searches for representations of
these same quantities in the historical coordinates
that have a simpler factor structure.

The main objective is NOT to rediscover the fact that
4v-1 factors.

The objective is to find a new identity where a modern
KAPPA observable factors through one of the historical
linear coordinates:

    A0
    B0
    A0-B0
    2x-3
    2y

or through a particularly simple shifted H2(t).

Such an identity would constitute an actual cross-era
structural bridge rather than another factorization
reformulation.
"""
)


# =============================================================================
# [14] FINAL STATUS
# =============================================================================

print("\n[14] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  historical basic identities = {basic_failures == 0}"
)
print(
    f"  branch relations           = {branch_failures == 0}"
)
print(
    f"  numerical crosswalk        = {numeric_failures == 0}"
)
print("  divisor enumeration        = False")
print("  continued fractions       = False")
print("  primary focus              = (x,y) <-> KAPPA")
print()
print("  PRIMARY QUESTION:")
print("    Does the historical conic coordinate system expose")
print("    a genuinely new low-degree representation of")
print("    F2, Delta, M1, or shifted H2(t)?")

print()
print("=" * 79)
print("EXPERIMENT 510 FINISHED")
print("=" * 79)