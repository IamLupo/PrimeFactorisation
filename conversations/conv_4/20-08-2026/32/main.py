#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 508
==============================================================================
2020-2026 HISTORICAL CONIC / KAPPA COORDINATE CROSSWALK
==============================================================================

Goals:

1. Verify the N-only historical quantities b, v, g.
2. Verify the exact conic

       y^2 - x^2 + 3x - 2 = v

3. Verify its difference-of-squares form

       (2y - 2x + 3)(2y + 2x - 3) = 4v - 1.

4. Determine valid historical (x,y) coordinates robustly without assuming
   the old y > 2x / y < 2x branch condition.

5. Recover the factor-bearing coordinates

       A0 = 2y - 2x + 3
       B0 = 2y + 2x - 3.

6. Normalize the scale-3 branch.

7. Recover the modern coordinates

       S     = p+q
       Delta = (p-q)^2.

8. Audit the historical gcd-bearing linear forms.

9. Search the integer conic solution space for the given examples.

This experiment does NOT assume that every historically stated branch
inequality is globally valid. It validates the actual algebraic identities.
==============================================================================

"""

from __future__ import annotations

import math
from typing import Optional


# ============================================================================
# DATA
# ============================================================================

PRIME_CASES = [
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


# ============================================================================
# BASIC HISTORICAL N-ONLY FUNCTIONS
# ============================================================================

def b_value(n: int) -> int:
    return n // 2


def v_value(n: int) -> int:
    b = b_value(n)
    return (b * b) % n


def g_value(n: int) -> int:
    b = b_value(n)
    return b // 2 + 1


def closed_v(n: int) -> int:
    r = n % 4

    if r == 1:
        return (3 * n + 1) // 4

    if r == 3:
        return (n + 1) // 4

    raise ValueError("closed_v requires odd n")


def closed_g(n: int) -> int:
    r = n % 4

    if r == 1:
        return (n + 3) // 4

    if r == 3:
        return (n + 1) // 4

    raise ValueError("closed_g requires odd n")


def scale_from_mod4(n: int) -> int:
    """
    4v - 1 = scale * N

    N == 3 mod 4 -> scale = 1
    N == 1 mod 4 -> scale = 3
    """
    r = n % 4

    if r == 1:
        return 3

    if r == 3:
        return 1

    raise ValueError("scale_from_mod4 requires odd n")


# ============================================================================
# HISTORICAL COORDINATE CANDIDATES
# ============================================================================

def candidate_coordinate_systems(p: int, q: int):
    """
    Generate all historically relevant coordinate forms.

    The old formulas were:

        opposite:
            x = (-p + q + 6)/4
            y = ( p + q    )/4

        branch A:
            x = ( 3p - q + 6)/4
            y = ( 3p + q    )/4

        branch B:
            x = (-3p + q + 6)/4
            y = ( 3p + q    )/4

    Because the old branch inequality can fail for some valid examples,
    we test both prime orientations p,q and all three formulas.
    """

    raw = []

    orientations = [
        ("pq", p, q),
        ("qp", q, p),
    ]

    for orientation, P, Q in orientations:

        formulas = [
            (
                "v==g",
                Q - P + 6,
                P + Q,
            ),
            (
                "v!=g:A",
                3 * P - Q + 6,
                3 * P + Q,
            ),
            (
                "v!=g:B",
                Q - 3 * P + 6,
                3 * P + Q,
            ),
        ]

        for branch, xn, yn in formulas:

            if xn % 4 != 0 or yn % 4 != 0:
                continue

            x = xn // 4
            y = yn // 4

            raw.append(
                {
                    "orientation": orientation,
                    "branch": branch,
                    "x": x,
                    "y": y,
                }
            )

    return raw


# ============================================================================
# HISTORICAL COORDINATE MAP
# ============================================================================

def A0(x: int, y: int) -> int:
    return 2 * y - 2 * x + 3


def B0(x: int, y: int) -> int:
    return 2 * y + 2 * x - 3


# ============================================================================
# VALIDATE A CANDIDATE COORDINATE SYSTEM
# ============================================================================

def validate_coordinate_candidate(
    n: int,
    p: int,
    q: int,
    x: int,
    y: int,
):
    v = v_value(n)
    scale = scale_from_mod4(n)

    conic = y * y - x * x + 3 * x - 2
    a = A0(x, y)
    b = B0(x, y)

    return {
        "conic_ok": conic == v,
        "factor_ok": a * b == scale * n,
        "product": a * b,
        "scale_n": scale * n,
        "a": a,
        "b": b,
        "conic_value": conic,
    }


def choose_historical_coordinates(
    p: int,
    q: int,
) -> Optional[dict]:

    n = p * q

    candidates = candidate_coordinate_systems(p, q)

    valid = []

    for candidate in candidates:
        result = validate_coordinate_candidate(
            n,
            p,
            q,
            candidate["x"],
            candidate["y"],
        )

        if result["conic_ok"] and result["factor_ok"]:
            row = dict(candidate)
            row.update(result)
            valid.append(row)

    if not valid:
        return None

    # Prefer the historical branch naming where possible.
    priority = {
        "v==g": 0,
        "v!=g:A": 1,
        "v!=g:B": 2,
    }

    valid.sort(
        key=lambda item: (
            priority.get(item["branch"], 99),
            item["orientation"],
        )
    )

    return valid[0]


# ============================================================================
# NORMALIZATION
# ============================================================================

def normalize_factor_coordinates(
    n: int,
    a: int,
    b: int,
):
    """
    Convert A0,B0 into the actual factor pair.

    scale=1:
        A0*B0=N

    scale=3:
        A0*B0=3N

    In the scale-3 case, one coordinate must carry the factor 3.
    """

    scale = scale_from_mod4(n)

    if scale == 1:
        if a * b != n:
            return None

        return tuple(sorted((a, b)))

    if a * b != 3 * n:
        return None

    if a % 3 == 0:
        return tuple(sorted((a // 3, b)))

    if b % 3 == 0:
        return tuple(sorted((a, b // 3)))

    return None


# ============================================================================
# HISTORICAL GCD FORMS
# ============================================================================

def gf1(n: int, x: int, y: int) -> int:
    return n // 2 - x - y + 2


def gf2(n: int, x: int, y: int) -> int:
    return n // 2 - x + y + 2


def gf3(n: int, x: int, y: int) -> int:
    return n // 2 + x - y - 1


def gf4(n: int, x: int, y: int) -> int:
    return n // 2 + x + y - 1


# ============================================================================
# SMALL CONIC ENUMERATOR
# ============================================================================

def enumerate_conic_solutions(
    n: int,
    max_abs_x: int = 100,
):
    """
    Enumerate integer x with bounded |x| and recover y if integral:

        y^2 = x^2 - 3x + 2 + v.
    """

    v = v_value(n)

    solutions = []

    for x in range(-max_abs_x, max_abs_x + 1):
        rhs = x * x - 3 * x + 2 + v

        if rhs < 0:
            continue

        y = math.isqrt(rhs)

        if y * y != rhs:
            continue

        solutions.append((x, y))

        if y != 0:
            solutions.append((x, -y))

    return sorted(set(solutions))


# ============================================================================
# START
# ============================================================================

print("=" * 78)
print("EXPERIMENT 508 START")
print("=" * 78)
print("2020-2026 HISTORICAL CONIC / KAPPA COORDINATE CROSSWALK")
print("=" * 78)


# ============================================================================
# [1]
# ============================================================================

print()
print("[1] N-ONLY HISTORICAL CHANNEL")
print("-" * 78)

failures = 0

for p, q in PRIME_CASES:
    n = p * q

    b = b_value(n)
    v = v_value(n)
    g = g_value(n)

    cv = closed_v(n)
    cg = closed_g(n)

    ok = (v == cv and g == cg)

    print(
        f"  N={n:<14} "
        f"b={b:<14} "
        f"v={v:<14} "
        f"g={g:<14} "
        f"PASS={ok}"
    )

    if not ok:
        failures += 1

print(f"  FAILURES = {failures}")


# ============================================================================
# [2]
# ============================================================================

print()
print("[2] 4v-1 SCALE LAW")
print("-" * 78)

failures = 0

for p, q in PRIME_CASES:
    n = p * q
    v = v_value(n)

    scale = scale_from_mod4(n)

    lhs = 4 * v - 1
    rhs = scale * n

    ok = lhs == rhs

    print(
        f"  N={n:<14} "
        f"4v-1={lhs:<18} "
        f"scale*N={rhs:<18} "
        f"scale={scale} "
        f"PASS={ok}"
    )

    if not ok:
        failures += 1

print(f"  SCALE FAILURES = {failures}")


# ============================================================================
# [3]
# ============================================================================

print()
print("[3] ROBUST HISTORICAL COORDINATE SELECTION")
print("-" * 78)

selection_failures = 0
selected = []

for p, q in PRIME_CASES:
    row = choose_historical_coordinates(p, q)

    if row is None:
        print(
            f"  ({p},{q}) "
            f"NO VALID HISTORICAL COORDINATE SYSTEM"
        )
        selection_failures += 1
        continue

    selected.append((p, q, row))

    print(
        f"  ({p},{q}) "
        f"orientation={row['orientation']:<3} "
        f"branch={row['branch']:<8} "
        f"x={row['x']:<10} "
        f"y={row['y']:<10} "
        f"A0={row['a']:<12} "
        f"B0={row['b']:<12} "
        f"PASS=True"
    )

print(
    f"  COORDINATE SELECTION FAILURES = "
    f"{selection_failures}"
)


# ============================================================================
# [4]
# ============================================================================

print()
print("[4] HISTORICAL COORDINATE IDENTITIES")
print("-" * 78)

identity_failures = 0

for p, q, row in selected:

    x = row["x"]
    y = row["y"]

    n = p * q
    v = v_value(n)
    scale = scale_from_mod4(n)

    a = A0(x, y)
    b = B0(x, y)

    lhs_conic = y * y - x * x + 3 * x - 2
    lhs_factor = a * b

    ok = (
        lhs_conic == v
        and lhs_factor == scale * n
        and b - a == 4 * x - 6
        and a + b == 4 * y
    )

    print(
        f"  ({p},{q}) "
        f"conic={lhs_conic == v} "
        f"factor={lhs_factor == scale*n} "
        f"inverse={b-a == 4*x-6 and a+b == 4*y} "
        f"PASS={ok}"
    )

    if not ok:
        identity_failures += 1

print(f"  IDENTITY FAILURES = {identity_failures}")


# ============================================================================
# [5]
# ============================================================================

print()
print("[5] NORMALIZED FACTOR CROSSWALK")
print("-" * 78)

crosswalk_failures = 0

for p, q, row in selected:

    n = p * q

    normalized = normalize_factor_coordinates(
        n,
        row["a"],
        row["b"],
    )

    if normalized is None:
        print(
            f"  ({p},{q}) "
            f"NORMALIZATION FAILED"
        )
        crosswalk_failures += 1
        continue

    u, vfactor = normalized

    S = u + vfactor
    Delta = (u - vfactor) ** 2

    ok = (
        {u, vfactor} == {p, q}
        and S == p + q
        and Delta == (p - q) ** 2
    )

    print(
        f"  ({p},{q}) "
        f"normalized=({u},{vfactor}) "
        f"S={S} "
        f"Delta={Delta} "
        f"PASS={ok}"
    )

    if not ok:
        crosswalk_failures += 1

print(
    f"  CROSSWALK FAILURES = "
    f"{crosswalk_failures}"
)


# ============================================================================
# [6]
# ============================================================================

print()
print("[6] HISTORICAL GCD-BEARING FORMS")
print("-" * 78)

gcd_failures = 0

for p, q, row in selected:

    n = p * q
    x = row["x"]
    y = row["y"]

    values = [
        gf1(n, x, y),
        gf2(n, x, y),
        gf3(n, x, y),
        gf4(n, x, y),
    ]

    gcds = [
        math.gcd(value, n)
        for value in values
    ]

    proper = [
        g
        for g in gcds
        if 1 < g < n
    ]

    ok = len(proper) > 0

    print(
        f"  ({p},{q}) "
        f"forms={values} "
        f"gcds={gcds} "
        f"proper={proper} "
        f"PASS={ok}"
    )

    if not ok:
        gcd_failures += 1

print(f"  GCD EXPOSURE FAILURES = {gcd_failures}")


# ============================================================================
# [7]
# ============================================================================

print()
print("[7] EXACT GCD REWRITE THROUGH A0/B0")
print("-" * 78)

rewrite_failures = 0

for p, q, row in selected:

    n = p * q
    x = row["x"]
    y = row["y"]

    a = A0(x, y)
    b = B0(x, y)

    lhs = [
        gf1(n, x, y),
        gf2(n, x, y),
        gf3(n, x, y),
        gf4(n, x, y),
    ]

    rhs = [
        (n - b) / 2,
        (n + a) / 2,
        (n - a) / 2,
        (n + b) / 2,
    ]

    ok = all(
        2 * lhs[i] == 2 * rhs[i]
        for i in range(4)
    )

    print(
        f"  ({p},{q}) "
        f"PASS={ok}"
    )

    if not ok:
        rewrite_failures += 1

print(
    f"  GCD REWRITE FAILURES = "
    f"{rewrite_failures}"
)


# ============================================================================
# [8]
# ============================================================================

print()
print("[8] HISTORICAL COORDINATES -> MODERN KAPPA TARGETS")
print("-" * 78)

modern_failures = 0

for p, q, row in selected:

    n = p * q

    normalized = normalize_factor_coordinates(
        n,
        row["a"],
        row["b"],
    )

    if normalized is None:
        modern_failures += 1
        continue

    u, w = normalized

    S0 = u + w
    Delta0 = (u - w) ** 2
    M1 = (u + 1) * (w + 1)

    target_S = p + q
    target_Delta = (p - q) ** 2
    target_M1 = (p + 1) * (q + 1)

    ok = (
        S0 == target_S
        and Delta0 == target_Delta
        and M1 == target_M1
    )

    print(
        f"  ({p},{q}) "
        f"S={S0} "
        f"Delta={Delta0} "
        f"M1={M1} "
        f"PASS={ok}"
    )

    if not ok:
        modern_failures += 1

print(
    f"  MODERN TARGET FAILURES = "
    f"{modern_failures}"
)


# ============================================================================
# [9]
# ============================================================================

print()
print("[9] SAME-N CONTROL")
print("-" * 78)

same_n_data = {
    12: [(2, 6), (3, 4)],
    18: [(2, 9), (3, 6)],
    20: [(2, 10), (4, 5)],
    30: [(3, 10), (5, 6)],
}

same_n_failures = 0

for n, pairs in same_n_data.items():

    b = b_value(n)
    v = v_value(n)
    g = g_value(n)

    sums = {p + q for p, q in pairs}
    gaps = {(p - q) ** 2 for p, q in pairs}

    print(
        f"  N={n:<4} "
        f"(b,v,g)=({b},{v},{g}) "
        f"S-values={sorted(sums)} "
        f"Delta-values={sorted(gaps)}"
    )

    ok = len(sums) > 1 and len(gaps) > 1

    print(
        f"      same N-only data, different S/Delta = "
        f"{ok}"
    )

    if not ok:
        same_n_failures += 1

print(
    f"  SAME-N CONTROL FAILURES = "
    f"{same_n_failures}"
)


# ============================================================================
# [10]
# ============================================================================

print()
print("[10] SMALL CONIC SOLUTION SEARCH")
print("-" * 78)

small_n_values = [15, 21, 35, 77, 105, 165]

for n in small_n_values:

    solutions = enumerate_conic_solutions(
        n,
        max_abs_x=100,
    )

    factor_solutions = []

    scale = scale_from_mod4(n)

    for x, y in solutions:

        a = A0(x, y)
        b = B0(x, y)

        if a * b != scale * n:
            continue

        normalized = normalize_factor_coordinates(
            n,
            a,
            b,
        )

        factor_solutions.append(
            (x, y, a, b, normalized)
        )

    print(
        f"  N={n:<4} "
        f"conic solutions={len(solutions):<4} "
        f"factor-bearing={len(factor_solutions)}"
    )

    for item in factor_solutions[:10]:
        x, y, a, b, normalized = item

        print(
            f"      x={x:<5} "
            f"y={y:<5} "
            f"A0={a:<5} "
            f"B0={b:<5} "
            f"normalized={normalized}"
        )


# ============================================================================
# [11]
# ============================================================================

print()
print("[11] CENTRAL ALGEBRAIC CROSSWALK")
print("-" * 78)

print(
    r"""
Define

    A0 = 2y - 2x + 3
    B0 = 2y + 2x - 3.

Then exactly:

    A0 + B0 = 4y
    B0 - A0 = 4x - 6

and

    A0 B0
      = 4(y^2 - x^2 + 3x - 2)
      = 4v - 1.

For odd N:

    N == 3 mod 4:
        4v - 1 = N

    N == 1 mod 4:
        4v - 1 = 3N.

After normalization of the scale-3 coordinate:

    factor_1 = p
    factor_2 = q

and therefore:

    S     = factor_1 + factor_2
    Delta = (factor_1-factor_2)^2
    M1    = (factor_1+1)(factor_2+1).

The key distinction is:

    N -> v

is completely source-free, while

    v -> (A0,B0)

requires choosing a particular integer factorization of
4v-1.

That selection step is exactly where the factor-bearing
information enters.
"""
)


# ============================================================================
# [12]
# ============================================================================

print()
print("[12] RESEARCH TARGET")
print("-" * 78)

print(
    r"""
The useful object is NOT another rewrite of v.

The real missing interface is:

        N
        |
        v
      4v-1
        |
        |  select integer factors
        v
      A0, B0
        |
        +----> A0+B0      -> S
        |
        +----> A0-B0      -> gap
        |
        +----> (A0-B0)^2  -> Delta
        |
        +----> (A0+1)(B0+1) -> M1.

Therefore the next experiments should compare the
existing historical / homogeneous-layer observables
against quantities attached to the FACTOR SELECTION
problem:

    divisor pairs of 4v-1
    near-square divisor pairs
    factor-pair sum
    factor-pair difference
    gcd-bearing linear forms.

The important question is whether any existing N-only
observable naturally selects the same divisor pair,
rather than whether it algebraically reproduces v.
"""
)


# ============================================================================
# [13]
# ============================================================================

overall = (
    selection_failures == 0
    and identity_failures == 0
    and crosswalk_failures == 0
    and gcd_failures == 0
    and rewrite_failures == 0
    and modern_failures == 0
    and same_n_failures == 0
)

print()
print("[13] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  coordinate selection          = "
    f"{selection_failures == 0}"
)

print(
    f"  historical identities         = "
    f"{identity_failures == 0}"
)

print(
    f"  normalized factor crosswalk   = "
    f"{crosswalk_failures == 0}"
)

print(
    f"  gcd exposure                  = "
    f"{gcd_failures == 0}"
)

print(
    f"  gcd rewrite                   = "
    f"{rewrite_failures == 0}"
)

print(
    f"  modern S/Delta/M1 crosswalk   = "
    f"{modern_failures == 0}"
)

print(
    f"  same-N control                = "
    f"{same_n_failures == 0}"
)

print(
    f"  OVERALL EXACT AUDIT           = "
    f"{overall}"
)


print()
print("=" * 78)
print("EXPERIMENT 508 FINISHED")
print("=" * 78)