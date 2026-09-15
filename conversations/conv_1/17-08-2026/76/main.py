#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 140
POST-BOUNDARY RECURRENCE + EXPLICIT LAURENT BOUNDARY SOURCE

PURPOSE
-------
Experiment 139 tested the recurrence on the FULL Laurent support and failed
because the recurrence crosses the finite Laurent boundary.

This experiment separates:

    1. the admissible post-boundary j-strip;
    2. the boundary-crossing step.

For fixed (k,ell),

    j = ell-k-d-1,

and d changes by 2, so admissible j-values have the same parity.

ODD d
-----
    A_j = (-1)^((d-1)/2) (k+ell)

hence, INSIDE the admissible strip,

    A_(j+2) + A_j = 0.

At the outer boundary this recurrence is forced by the missing next term.
We explicitly measure that boundary source.

EVEN d
------
    A_j = (-1)^(d/2+1) (k+ell) j/2

hence, inside the admissible strip,

    A_(j+4) + 2 A_(j+2) + A_j = 0.

Again we separately measure the first boundary crossing.

GOAL
----
Show:

    admissible recurrence = exact,
    boundary source = simple closed form.

This directly repairs the over-strong global recurrence tested in 139.

NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO FITTING
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")
x, t, u = sp.symbols("x t u")


# ============================================================================
# Exact quotient
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q)
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q)
    )

    Q, R = sp.div(PF, PD)

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = {sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric conversion
# ============================================================================

def symmetric_to_NS(expr):
    expr = sp.expand(
        expr.subs(q, S - p)
    )

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N)
    )

    poly = sp.Poly(
        expr,
        p,
        domain=sp.QQ.frac_field(S, N)
    )

    rem = sp.expand(
        sp.rem(poly, modulus).as_expr()
    )

    residual_p = sp.expand(
        rem.coeff(p, 1)
    )

    if residual_p != 0:
        raise ArithmeticError(
            f"non-symmetric remainder: {sp.factor(residual_p)}"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


# ============================================================================
# Weighted layers
# ============================================================================

def weighted_degree(term):
    powers = term.as_powers_dict()

    s_exp = sp.sympify(
        powers.get(S, 0)
    )
    n_exp = sp.sympify(
        powers.get(N, 0)
    )

    if not s_exp.is_Integer or not n_exp.is_Integer:
        raise ArithmeticError(
            f"non-integer weighted exponent in {term}"
        )

    return int(s_exp) + 2 * int(n_exp)


def weighted_layers(expr):
    layers = {}

    for term in sp.Add.make_args(
        sp.expand(expr)
    ):
        w = weighted_degree(term)
        layers[w] = sp.expand(
            layers.get(w, 0) + term
        )

    return {
        w: sp.expand(v)
        for w, v in layers.items()
        if sp.expand(v) != 0
    }


def layer_to_x(layer):
    scaled = sp.expand(
        layer.subs(
            {
                S: t*x,
                N: t**2
            }
        )
    )

    poly = sp.Poly(
        scaled,
        t,
        domain=sp.QQ.frac_field(x)
    )

    return sp.expand(
        poly.LC()
    )


def get_layers(k, ell):
    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    layers = weighted_layers(Q)

    top = max(layers)

    top_H = layer_to_x(
        layers[top]
    )

    lower_H = sp.Integer(0)

    if top - 1 in layers:
        lower_H = layer_to_x(
            layers[top - 1]
        )

    return top, top_H, lower_H


# ============================================================================
# Laurent coefficient
# ============================================================================

def laurent_coeff(H, j):
    """
    Exact [u^j] H(u+u^-1) using the binomial identity.
    """

    poly = sp.Poly(
        sp.expand(H),
        x,
        domain=sp.QQ
    )

    total = sp.Integer(0)

    for (m,), coeff in poly.terms():

        m = int(m)

        if m < j:
            continue

        if (m - j) % 2:
            continue

        r = (m - j) // 2

        total += coeff * sp.binomial(
            m,
            r
        )

    return sp.expand(total)


# ============================================================================
# Admissible j values
# ============================================================================

def admissible_rows(k, ell):
    rows = []

    for d in range(
        1,
        ell - k - 1
    ):
        j = ell - k - d - 1

        if j < 1:
            continue

        rows.append(
            {
                "d": d,
                "j": j
            }
        )

    return rows


# ============================================================================
# Exact target formulas
# ============================================================================

def expected_A(k, ell, d):
    j = ell - k - d - 1

    if d % 2 == 1:
        return sp.expand(
            (-1) ** ((d - 1) // 2)
            * (k + ell)
        )

    return sp.expand(
        (-1) ** (d // 2 + 1)
        * sp.Rational(k + ell, 2)
        * j
    )


# ============================================================================
# Analyze one pair
# ============================================================================

def analyze_pair(k, ell):

    top_weight, top_H, lower_H = get_layers(
        k,
        ell
    )

    rows = admissible_rows(
        k,
        ell
    )

    odd_rows = [
        r for r in rows
        if r["d"] % 2 == 1
    ]

    even_rows = [
        r for r in rows
        if r["d"] % 2 == 0
    ]

    # ------------------------------------------------------------
    # Exact admissible coefficients
    # ------------------------------------------------------------

    for r in rows:

        d = r["d"]
        j = r["j"]

        H = (
            top_H
            if d % 2 == 1
            else lower_H
        )

        r["A"] = laurent_coeff(
            H,
            j
        )

        r["expected"] = expected_A(
            k,
            ell,
            d
        )

        r["law_ok"] = (
            sp.expand(
                r["A"] - r["expected"]
            ) == 0
        )

    # ------------------------------------------------------------
    # Odd branch: A_(j+2)+A_j = 0
    # ONLY where both indices are admissible.
    # ------------------------------------------------------------

    odd_rows_sorted = sorted(
        odd_rows,
        key=lambda z: z["j"]
    )

    odd_interior = []
    odd_boundary = None

    for a, b in zip(
        odd_rows_sorted,
        odd_rows_sorted[1:]
    ):

        if b["j"] == a["j"] + 2:
            residual = sp.expand(
                b["A"] + a["A"]
            )

            odd_interior.append(
                {
                    "j": a["j"],
                    "next_j": b["j"],
                    "residual": residual
                }
            )

    # First step that crosses beyond the admissible maximum j.
    if odd_rows_sorted:

        max_row = max(
            odd_rows_sorted,
            key=lambda z: z["j"]
        )

        j0 = max_row["j"]
        boundary_next = j0 + 2

        boundary_value = laurent_coeff(
            top_H,
            boundary_next
        )

        odd_boundary = {
            "j": j0,
            "next_j": boundary_next,
            "A_j": max_row["A"],
            "A_next": boundary_value,
            "residual": sp.expand(
                boundary_value + max_row["A"]
            )
        }

    # ------------------------------------------------------------
    # Even branch: second difference with step 2
    # ONLY where all three j-values are admissible.
    # ------------------------------------------------------------

    even_rows_sorted = sorted(
        even_rows,
        key=lambda z: z["j"]
    )

    even_interior = []

    lookup_even = {
        r["j"]: r
        for r in even_rows_sorted
    }

    even_js = sorted(
        lookup_even.keys()
    )

    for j0 in even_js:

        if (
            j0 + 2 in lookup_even
            and j0 + 4 in lookup_even
        ):

            A0 = lookup_even[j0]["A"]
            A2 = lookup_even[j0 + 2]["A"]
            A4 = lookup_even[j0 + 4]["A"]

            residual = sp.expand(
                A4 + 2 * A2 + A0
            )

            even_interior.append(
                {
                    "j": j0,
                    "residual": residual
                }
            )

    even_boundary = None

    if even_rows_sorted:

        max_j = max(
            r["j"]
            for r in even_rows_sorted
        )

        if max_j + 4 <= top_weight:
            A0 = laurent_coeff(
                lower_H,
                max_j
            )
            A2 = laurent_coeff(
                lower_H,
                max_j + 2
            )
            A4 = laurent_coeff(
                lower_H,
                max_j + 4
            )

            even_boundary = {
                "j": max_j,
                "residual": sp.expand(
                    A4 + 2 * A2 + A0
                )
            }

    return {
        "k": k,
        "ell": ell,
        "top_weight": top_weight,
        "top_H": top_H,
        "lower_H": lower_H,
        "rows": rows,
        "odd_interior": odd_interior,
        "odd_boundary": odd_boundary,
        "even_interior": even_interior,
        "even_boundary": even_boundary
    }


# ============================================================================
# Print pair
# ============================================================================

def print_pair(result):

    k = result["k"]
    ell = result["ell"]

    print(
        f"({k},{ell})"
    )

    print(
        f"  top H(x)   = {sp.factor(result['top_H'])}"
    )

    print(
        f"  lower H(x) = {sp.factor(result['lower_H'])}"
    )

    print()

    odd_fail = [
        z for z in result["odd_interior"]
        if z["residual"] != 0
    ]

    even_fail = [
        z for z in result["even_interior"]
        if z["residual"] != 0
    ]

    law_fail = [
        r for r in result["rows"]
        if not r["law_ok"]
    ]

    print(
        f"  admissible-law failures = {len(law_fail)}"
    )

    print(
        f"  odd interior recurrence failures = {len(odd_fail)}"
    )

    print(
        f"  even interior recurrence failures = {len(even_fail)}"
    )

    if result["odd_interior"]:
        print(
            "  odd interior residuals:"
        )
        for z in result["odd_interior"]:
            print(
                f"    j={z['j']} -> {z['next_j']}: "
                f"{sp.factor(z['residual'])}"
            )

    if result["even_interior"]:
        print(
            "  even interior residuals:"
        )
        for z in result["even_interior"]:
            print(
                f"    j={z['j']}: "
                f"{sp.factor(z['residual'])}"
            )

    if result["odd_boundary"] is not None:
        z = result["odd_boundary"]

        print(
            "  odd boundary crossing:"
        )

        print(
            f"    j={z['j']} -> {z['next_j']}"
        )

        print(
            f"    A_j = {sp.factor(z['A_j'])}"
        )

        print(
            f"    A_next = {sp.factor(z['A_next'])}"
        )

        print(
            f"    boundary residual = "
            f"{sp.factor(z['residual'])}"
        )

    print()


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 140")
    print("POST-BOUNDARY RECURRENCE + EXPLICIT LAURENT BOUNDARY")
    print("ADMISSIBLE j STRIP ONLY")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO POLYNOMIAL FITTING")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    PAIRS = [
        (1, 7),
        (1, 9),
        (1, 11),
        (3, 9),
        (3, 11),
        (5, 11),
        (1, 13),
        (3, 13)
    ]

    results = []

    total_law_failures = 0
    total_odd_failures = 0
    total_even_failures = 0

    boundary_values = []

    print("=" * 78)
    print("1. ADMISSIBLE-STRIP RECURRENCE TEST")
    print("=" * 78)

    for k, ell in PAIRS:

        result = analyze_pair(
            k,
            ell
        )

        results.append(
            result
        )

        print_pair(
            result
        )

        total_law_failures += sum(
            not r["law_ok"]
            for r in result["rows"]
        )

        total_odd_failures += sum(
            z["residual"] != 0
            for z in result["odd_interior"]
        )

        total_even_failures += sum(
            z["residual"] != 0
            for z in result["even_interior"]
        )

        if result["odd_boundary"] is not None:
            boundary_values.append(
                (
                    k,
                    ell,
                    result["odd_boundary"]["residual"]
                )
            )

    print("=" * 78)
    print("2. BOUNDARY SOURCE SUMMARY")
    print("=" * 78)

    for k, ell, value in boundary_values:
        print(
            f"({k},{ell}) boundary source = "
            f"{sp.factor(value)}"
        )

    print()

    print("=" * 78)
    print("3. GLOBAL SUMMARY")
    print("=" * 78)

    print(
        f"pairs tested = {len(results)}"
    )

    print(
        f"admissible law failures = "
        f"{total_law_failures}"
    )

    print(
        f"odd interior failures = "
        f"{total_odd_failures}"
    )

    print(
        f"even interior failures = "
        f"{total_even_failures}"
    )

    print()

    # Check the observed boundary source.
    # The first boundary crossing for the top layer is expected
    # to produce +(k+ell) for the tested orientation.
    boundary_failures = []

    for k, ell, residual in boundary_values:

        expected = sp.Integer(
            k + ell
        )

        if sp.expand(
            residual - expected
        ) != 0:

            boundary_failures.append(
                (
                    k,
                    ell,
                    residual,
                    expected
                )
            )

    print(
        f"boundary-source failures = "
        f"{len(boundary_failures)}"
    )

    if boundary_failures:

        print(
            "FIRST BOUNDARY FAILURES"
        )

        for item in boundary_failures[:8]:
            print(
                f"  ({item[0]},{item[1]}) "
                f"got={item[2]} "
                f"expected={item[3]}"
            )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    total_failures = (
        total_law_failures
        + total_odd_failures
        + total_even_failures
        + len(boundary_failures)
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The recurrence holds on the admissible post-boundary"
        )
        print(
            "j-strip, while the failure observed in Experiment 139"
        )
        print(
            "is isolated to the first boundary crossing."
        )
        print()
        print(
            "The next target is therefore the exact symbolic"
        )
        print(
            "boundary coefficient, followed by the recurrence"
        )
        print(
            "induction."
        )

    else:

        print("STATUS = FAIL")
        print(
            f"total failures = {total_failures}"
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

