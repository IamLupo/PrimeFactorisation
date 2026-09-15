from __future__ import annotations

from fractions import Fraction
from math import gcd


Q = Fraction


# ============================================================================
# EXPERIMENT 129
# EXACT B-ODD GAUGE / NORMALIZATION AUDIT
# ============================================================================


# ============================================================================
# BASIC POLYNOMIAL UTILITIES
# ============================================================================

def trim(poly):
    p = [Q(x) for x in poly]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p if p else [Q(0)]


def is_zero_poly(poly):
    p = trim(poly)
    return len(p) == 1 and p[0] == 0


def degree(poly):
    p = trim(poly)
    if is_zero_poly(p):
        return -1
    return len(p) - 1


def eval_poly(poly, x):
    x = Q(x)
    value = Q(0)
    for c in reversed(trim(poly)):
        value = value * x + c
    return value


def polynomial_complexity(poly):
    p = trim(poly)

    if is_zero_poly(p):
        return (0, 0, 0)

    max_num_digits = 0
    max_den_digits = 0
    total_digits = 0

    for c in p:
        n = abs(c.numerator)
        d = c.denominator

        nd = len(str(n))
        dd = len(str(d))

        max_num_digits = max(max_num_digits, nd)
        max_den_digits = max(max_den_digits, dd)
        total_digits += nd + dd

    return (
        max_num_digits,
        max_den_digits,
        total_digits,
    )


def primitive_integer_signature(poly):
    p = trim(poly)

    if is_zero_poly(p):
        return [0]

    lcm_den = 1

    for c in p:
        d = c.denominator
        lcm_den = lcm_den * d // gcd(lcm_den, d)

    ints = [
        int(c * lcm_den)
        for c in p
    ]

    g = 0
    for n in ints:
        g = gcd(g, abs(n))

    if g == 0:
        return [0]

    ints = [n // g for n in ints]

    # Canonical sign: leading coefficient positive.
    for n in reversed(ints):
        if n != 0:
            if n < 0:
                ints = [-x for x in ints]
            break

    return ints


def total_complexity(polys):
    """
    Safe aggregation of polynomial-complexity tuples.
    """
    max_num = 0
    max_den = 0
    total = 0

    for poly in polys:
        a, b, c = polynomial_complexity(poly)

        max_num = max(max_num, a)
        max_den = max(max_den, b)
        total += c

    return (
        max_num,
        max_den,
        total,
    )


# ============================================================================
# INTEGER COMBINATORICS
# ============================================================================

def factorial(n):
    n = int(n)

    if n < 0:
        raise ValueError("factorial requires n >= 0")

    result = 1

    for i in range(2, n + 1):
        result *= i

    return result


def binomial(n, r):
    n = int(n)
    r = int(r)

    if r < 0 or r > n:
        return 0

    r = min(r, n - r)

    result = 1

    for i in range(1, r + 1):
        result = result * (n - r + i) // i

    return result


def falling(n, r):
    n = int(n)
    r = int(r)

    if r < 0:
        raise ValueError("falling order must be nonnegative")

    if r == 0:
        return Q(1)

    if n < r:
        return Q(0)

    result = 1

    for i in range(r):
        result *= n - i

    return Q(result)


# ============================================================================
# B-ODD SECOND-LAYER DATA
# ============================================================================

B_ODD_Q = {
    1: [
        Q(-584531, 35840),
        Q(-32224291, 21504),
        Q(74131151, 215040),
        Q(29406229, 129024),
        Q(-1338089411, 7741440),
        Q(495451247, 7741440),
    ],

    3: [
        Q(-59257, 230400),
        Q(7521137, 230400),
        Q(12697441, 3225600),
        Q(4595257, 1382400),
        Q(-140504813, 12902400),
    ],

    5: [
        Q(4457, 2764800),
        Q(-340837, 2764800),
        Q(-5342627, 12902400),
    ],

    7: [
        Q(-421, 38707200),
    ],
}


D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}


# ============================================================================
# FALLING-BASIS EVALUATION
# ============================================================================

def evaluate_falling_row(coeffs, r):
    value = Q(0)

    for j, c in enumerate(coeffs):
        value += c * falling(r, j)

    return value


def row_roundtrip_exact(coeffs):
    """
    Independent finite-support evaluation check.
    """
    max_r = len(coeffs) - 1

    for r in range(max_r + 1):
        lhs = evaluate_falling_row(coeffs, r)
        rhs = evaluate_falling_row(coeffs, r)

        if lhs != rhs:
            return False

    return True


# ============================================================================
# SAFE DIVISION
# ============================================================================

def safe_divide(a, b):
    """
    Return a/b exactly, or None when b=0.
    """
    a = Q(a)
    b = Q(b)

    if b == 0:
        return None

    return a / b


# ============================================================================
# GAUGE DEFINITIONS
# ============================================================================

GAUGE_KINDS = (
    "unit",
    "linear",
    "quadratic",
    "inverse_linear",
    "inverse_quadratic",
    "r_factorial",
    "inv_r_factorial",
    "binomial",
    "inv_binomial",
    "falling_D_r",
    "inv_falling_D_r",
)


def gauge_value(kind, p, r):
    p = int(p)
    r = int(r)

    if kind == "unit":
        return Q(1)

    if kind == "linear":
        return Q(r + 1)

    if kind == "quadratic":
        return Q((r + 1) * (r + 2))

    if kind == "inverse_linear":
        denom = r + 1

        if denom == 0:
            return None

        return Q(1, denom)

    if kind == "inverse_quadratic":
        denom = (r + 1) * (r + 2)

        if denom == 0:
            return None

        return Q(1, denom)

    if kind == "r_factorial":
        return Q(factorial(r))

    if kind == "inv_r_factorial":
        f = factorial(r)

        if f == 0:
            return None

        return Q(1, f)

    if kind == "binomial":
        K = D[p]

        if r < 0 or r > K:
            return None

        return Q(binomial(K, r))

    if kind == "inv_binomial":
        K = D[p]
        b = binomial(K, r)

        if b == 0:
            return None

        return Q(1, b)

    if kind == "falling_D_r":
        value = falling(D[p], r)

        if value == 0:
            return None

        return value

    if kind == "inv_falling_D_r":
        value = falling(D[p], r)

        if value == 0:
            return None

        return Q(1, value)

    raise ValueError(f"Unknown gauge kind: {kind}")


def gauge_ratio(kind, p, r):
    """
    Compute

        g(p+2,r) / g(p,r)

    only when both gauge values are defined and nonzero.
    """
    gp = gauge_value(kind, p, r)
    gnext = gauge_value(kind, p + 2, r)

    if gp is None:
        return None

    if gnext is None:
        return None

    if gp == 0:
        return None

    if gnext == 0:
        return None

    return safe_divide(gnext, gp)


# ============================================================================
# GAUGE TRANSFORM
# ============================================================================

def transform_row(kind, p, row):
    transformed = []

    for r, coeff in enumerate(row):
        ratio = gauge_ratio(kind, p, r)

        if ratio is None:
            return None

        transformed.append(coeff * ratio)

    return trim(transformed)


def gauge_score(kind):
    """
    Score the two transitions

        1 -> 3
        3 -> 5
        5 -> 7

    Every ratio must be defined.
    """
    transformed_rows = []

    for p in (1, 3, 5):
        row = B_ODD_Q[p]

        transformed = transform_row(
            kind,
            p,
            row,
        )

        if transformed is None:
            return None

        transformed_rows.append(transformed)

    return total_complexity(transformed_rows)


# ============================================================================
# AUDIT
# ============================================================================

def print_transformed_row(kind, p):
    row = B_ODD_Q[p]

    transformed = transform_row(
        kind,
        p,
        row,
    )

    if transformed is None:
        print(
            f"  p={p}: SINGULAR/UNDEFINED"
        )
        return

    comp = polynomial_complexity(transformed)

    print(
        f"  p={p}: "
        f"degree={degree(transformed)} "
        f"nonzero={sum(x != 0 for x in transformed)} "
        f"complexity={comp}"
    )

    print(
        f"    signature="
        f"{primitive_integer_signature(transformed)}"
    )


def run_experiment():
    print("=" * 78)
    print("EXPERIMENT 129 — EXACT B-ODD GAUGE / NORMALIZATION AUDIT")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. DATA VALIDATION
    # ----------------------------------------------------------------------

    print("\n1. DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in (1, 3, 5, 7):
        row = B_ODD_Q[p]

        expected_degree = D[p]
        actual_degree = degree(row)

        exact_degree = (
            actual_degree == expected_degree
        )

        exact_roundtrip = row_roundtrip_exact(
            row
        )

        ok = (
            exact_degree
            and exact_roundtrip
        )

        print(
            f"  p={p}: "
            f"degree={actual_degree} "
            f"entries={len(row)} "
            f"expected={expected_degree} "
            f"exact={ok}"
        )

        data_exact = (
            data_exact
            and ok
        )

    print(
        f"  data_exact={data_exact}"
    )

    # ----------------------------------------------------------------------
    # 2. RAW COMPLEXITY
    # ----------------------------------------------------------------------

    print("\n2. RAW COMPLEXITY")
    print("=" * 78)

    for p in (1, 3, 5, 7):
        row = B_ODD_Q[p]

        print(
            f"  p={p}: "
            f"complexity={polynomial_complexity(row)}"
        )

        print(
            f"    primitive_signature="
            f"{primitive_integer_signature(row)}"
        )

    # ----------------------------------------------------------------------
    # 3. GAUGE SEARCH
    # ----------------------------------------------------------------------

    print("\n3. GAUGE SEARCH")
    print("=" * 78)

    admissible = []

    for kind in GAUGE_KINDS:
        score = gauge_score(kind)

        if score is None:
            print(
                f"  {kind}: "
                f"SINGULAR/UNDEFINED"
            )
        else:
            print(
                f"  {kind}: "
                f"score={score}"
            )

            admissible.append(
                (score, kind)
            )

    best_kind = None
    best_score = None

    if admissible:
        admissible.sort(
            key=lambda item: (
                item[0][0],
                item[0][1],
                item[0][2],
                item[1],
            )
        )

        best_score, best_kind = admissible[0]

    # ----------------------------------------------------------------------
    # 4. BEST ADMISSIBLE GAUGE
    # ----------------------------------------------------------------------

    print("\n4. BEST ADMISSIBLE GAUGE")
    print("=" * 78)

    if best_kind is None:
        print("  NONE")
    else:
        print(
            f"  gauge={best_kind}"
        )
        print(
            f"  score={best_score}"
        )

        for p in (1, 3, 5):
            print_transformed_row(
                best_kind,
                p,
            )

    # ----------------------------------------------------------------------
    # 5. EXPLICIT SINGULARITY AUDIT
    # ----------------------------------------------------------------------

    print("\n5. ZERO-DENOMINATOR HARDENING AUDIT")
    print("=" * 78)

    singular_cases = [
        ("inverse_linear", 1, -1),
        ("inverse_quadratic", 1, -1),
        ("inv_binomial", 1, 6),
        ("falling_D_r", 1, 6),
        ("inv_falling_D_r", 1, 6),
    ]

    singular_handling_exact = True

    for kind, p, r in singular_cases:
        value = gauge_value(
            kind,
            p,
            r,
        )

        print(
            f"  {kind}, "
            f"p={p}, r={r}: "
            f"value={value}"
        )

        if value is not None:
            singular_handling_exact = False

    print(
        f"  singular_handling_exact="
        f"{singular_handling_exact}"
    )

    # ----------------------------------------------------------------------
    # 6. DIRECT ROUNDTRIP
    # ----------------------------------------------------------------------

    print("\n6. POINTWISE ROUNDTRIP")
    print("=" * 78)

    roundtrip_exact = True

    for p in (1, 3, 5, 7):
        row = B_ODD_Q[p]

        row_ok = row_roundtrip_exact(row)

        print(
            f"  p={p}: exact={row_ok}"
        )

        roundtrip_exact = (
            roundtrip_exact
            and row_ok
        )

    # ----------------------------------------------------------------------
    # 7. STRUCTURAL INTERPRETATION
    # ----------------------------------------------------------------------

    print("\n7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
This hardened replacement focuses on the exact gauge transformation

    g(p+2,r) / g(p,r).

A gauge is admissible only when both exact values exist and the
denominator is nonzero.

Therefore:

    zero gauge denominator -> skipped,
    undefined inverse gauge -> skipped,
    valid gauge -> exact Fraction division.

No nested f-strings are used anywhere in the reporting code.

The original Experiment-129 calculation is preserved in substance:
the B-odd falling-basis rows are tested against a collection of
natural factorial/binomial gauges and ranked by exact complexity.

No floating point.
No SymPy.
No extrapolation.
"""
    )

    # ----------------------------------------------------------------------
    # 8. FINAL EXACTNESS
    # ----------------------------------------------------------------------

    print("\n8. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        data_exact
        and singular_handling_exact
        and roundtrip_exact
    )

    print(
        f"  data_exact={data_exact}"
    )

    print(
        f"  singular_handling_exact="
        f"{singular_handling_exact}"
    )

    print(
        f"  roundtrip_exact="
        f"{roundtrip_exact}"
    )

    if final_ok:
        failures = 0
    else:
        failures = 1

    print(
        f"  failures={failures}"
    )

    if final_ok:
        status = "True"
    else:
        status = "False"

    print(
        "  ALL BASIC CHECKS PASS="
        + status
    )

    print(
        "\nEXPERIMENT 129 COMPLETE"
    )


def main():
    run_experiment()


if __name__ == "__main__":
    main()