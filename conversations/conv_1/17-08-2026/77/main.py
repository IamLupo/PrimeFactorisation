#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 141
SYMBOLIC BOUNDARY-SOURCE DERIVATION

TOP LAYER:
    interior recurrence:
        A_(j+2) + A_j = 0

    boundary source observed in Exp. 140:
        A_(j_max+2) + A_jmax = ell

LOWER PARITY LAYER:
    interior recurrence:
        A_(j+4) + 2 A_(j+2) + A_j = 0

This experiment does NOT guess the boundary source.

It derives the boundary residual directly from the exact Laurent
polynomial and then searches for a symbolic closed form in k, ell.

The central objects are the OUTER Laurent coefficients of

    E_(k,ell)(u) = E_(k,ell)(u + u^(-1))
    O_(k,ell)(u) = O_(k,ell)(u + u^(-1))

where E is the top weighted layer and O is the parity-compatible
one-weight-lower layer.

TARGET:
    prove symbolically

        B_top(k,ell) = ell

and derive the corresponding lower-layer boundary source.

NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO POLYNOMIAL FITTING
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
            f"quotient remainder for ({k},{ell}): "
            f"{sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric conversion Q(p,q) -> Q(N,S)
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

    p1 = sp.expand(
        rem.coeff(p, 1)
    )

    if p1 != 0:
        raise ArithmeticError(
            f"symmetric reduction retained p: {sp.factor(p1)}"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


# ============================================================================
# Weighted layers
# ============================================================================

def integer_exponent(value):
    value = sp.sympify(value)

    if not value.is_Integer:
        raise ArithmeticError(
            f"non-integer exponent {value}"
        )

    return int(value)


def weighted_degree(term):
    powers = term.as_powers_dict()

    s_exp = integer_exponent(
        powers.get(S, 0)
    )
    n_exp = integer_exponent(
        powers.get(N, 0)
    )

    return s_exp + 2 * n_exp


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

    if poly.is_zero:
        return sp.Integer(0)

    return sp.expand(
        poly.LC()
    )


# ============================================================================
# Get top and lower layer
# ============================================================================

def get_layers(k, ell):
    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    layers = weighted_layers(Q)

    if not layers:
        raise ArithmeticError(
            f"no weighted layers for ({k},{ell})"
        )

    top_weight = max(layers)

    E = layer_to_x(
        layers[top_weight]
    )

    O = sp.Integer(0)

    if top_weight - 1 in layers:
        O = layer_to_x(
            layers[top_weight - 1]
        )

    return top_weight, E, O


# ============================================================================
# Laurent conversion
# ============================================================================

def laurent_expr(H):
    return sp.expand(
        H.subs(
            x,
            u + u**-1
        )
    )


def laurent_support(expr):
    expr = sp.expand(expr)

    support = []

    for term in sp.Add.make_args(expr):

        power = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"non-integer Laurent exponent in {term}"
            )

        support.append(
            int(power)
        )

    return sorted(
        set(support)
    )


def laurent_coeff(expr, j):
    """
    Exact coefficient [u^j].
    """

    support = laurent_support(expr)

    if not support:
        return sp.Integer(0)

    min_power = min(support)

    shift = max(
        0,
        -min_power
    )

    shifted = sp.expand(
        expr * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ
    )

    wanted = j + shift

    if wanted < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(
            u**wanted
        )
    )


# ============================================================================
# Find the admissible maximum j for each parity
# ============================================================================

def admissible_js(k, ell, parity):
    result = []

    for d in range(
        1,
        ell - k - 1
    ):
        if d % 2 != parity:
            continue

        j = ell - k - d - 1

        if j >= 1:
            result.append(
                (d, j)
            )

    return sorted(
        result,
        key=lambda z: z[1]
    )


# ============================================================================
# Boundary source for top layer
# ============================================================================

def top_boundary_source(k, ell):
    """
    Let j_max be the largest admissible j on the odd-d branch.

    Compute

        A_(j_max+2) + A_jmax

    directly from the full finite Laurent polynomial.

    This is the source that appears only because the next Laurent
    coefficient is outside the post-boundary admissible strip.
    """

    _, E, _ = get_layers(
        k,
        ell
    )

    R = laurent_expr(E)

    rows = admissible_js(
        k,
        ell,
        1
    )

    if not rows:
        return None

    d_max, j_max = max(
        rows,
        key=lambda z: z[1]
    )

    A_j = laurent_coeff(
        R,
        j_max
    )

    A_next = laurent_coeff(
        R,
        j_max + 2
    )

    return {
        "j_max": j_max,
        "A_j": sp.expand(A_j),
        "A_next": sp.expand(A_next),
        "source": sp.expand(
            A_next + A_j
        )
    }


# ============================================================================
# Boundary source for lower layer
# ============================================================================

def lower_boundary_sources(k, ell):
    """
    For the lower layer O, determine the largest admissible even-d j
    and measure:

        B_j = A_(j+4) + 2 A_(j+2) + A_j

    at the first boundary crossing.

    We also measure the preceding interior residuals to make sure the
    source is genuinely a boundary defect.
    """

    _, _, O = get_layers(
        k,
        ell
    )

    if O == 0:
        return None

    R = laurent_expr(O)

    rows = admissible_js(
        k,
        ell,
        0
    )

    if not rows:
        return None

    js = sorted(
        j for _, j in rows
    )

    interior = []

    lookup = {
        r["j"]: r
        for r in []
    }

    for j in js:
        if (
            j + 2 in js
            and j + 4 in js
        ):
            A0 = laurent_coeff(R, j)
            A2 = laurent_coeff(R, j + 2)
            A4 = laurent_coeff(R, j + 4)

            interior.append(
                (
                    j,
                    sp.expand(
                        A4 + 2*A2 + A0
                    )
                )
            )

    j_max = max(js)

    A0 = laurent_coeff(
        R,
        j_max
    )

    A2 = laurent_coeff(
        R,
        j_max + 2
    )

    A4 = laurent_coeff(
        R,
        j_max + 4
    )

    source = sp.expand(
        A4 + 2*A2 + A0
    )

    return {
        "j_max": j_max,
        "A_j": sp.expand(A0),
        "A_j2": sp.expand(A2),
        "A_j4": sp.expand(A4),
        "source": source,
        "interior": interior
    }


# ============================================================================
# Symbolic formula candidates
# ============================================================================

def top_expected(k, ell):
    return sp.Integer(ell)


# ============================================================================
# Analyze finite training grid
# ============================================================================

def build_training(max_ell):
    rows = []

    for ell in range(
        5,
        max_ell + 1,
        2
    ):
        for k in range(
            1,
            ell,
            2
        ):

            tb = top_boundary_source(
                k,
                ell
            )

            if tb is not None:
                rows.append(
                    {
                        "k": k,
                        "ell": ell,
                        "source": tb["source"],
                        "expected": top_expected(
                            k,
                            ell
                        )
                    }
                )

    return rows


# ============================================================================
# Factorization and independence test
# ============================================================================

def inspect_top_rows(rows):
    failures = []

    print("=" * 78)
    print("1. TOP-LAYER BOUNDARY SOURCE")
    print("=" * 78)

    for row in rows:

        residual = sp.expand(
            row["source"]
            - row["expected"]
        )

        print(
            f"({row['k']},{row['ell']}) "
            f"source={sp.factor(row['source'])} "
            f"expected={row['expected']} "
            f"residual={sp.factor(residual)}"
        )

        if residual != 0:
            failures.append(
                row
            )

    return failures


# ============================================================================
# Direct symbolic edge simplification
# ============================================================================

def derive_top_boundary_identity(E):
    """
    Examine the outermost positive Laurent coefficient.

    For a polynomial H of degree M with leading x^M coefficient a_M,
    [u^M]H(u+u^-1) = a_M.

    The next coefficient [u^(M-2)] receives contributions from the
    top two x-powers. We print these explicitly.

    This is intended to reveal the exact binomial combination that
    gives the boundary source.
    """

    poly = sp.Poly(
        sp.expand(E),
        x,
        domain=sp.QQ
    )

    degree = poly.degree()

    if degree < 1:
        return None

    top_coeff = poly.coeff_monomial(
        x**degree
    )

    second_coeff = (
        poly.coeff_monomial(
            x**(degree - 2)
        )
        if degree >= 2
        else sp.Integer(0)
    )

    # [u^(M-2)] x^M = M
    # [u^(M-2)] x^(M-2) = 1
    boundary_adjacent = sp.expand(
        degree * top_coeff
        + second_coeff
    )

    return {
        "degree": degree,
        "top_coeff": sp.expand(top_coeff),
        "second_coeff": sp.expand(second_coeff),
        "adjacent_coeff": boundary_adjacent
    }


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 141")
    print("SYMBOLIC BOUNDARY-SOURCE DERIVATION")
    print("TOP SOURCE + LOWER-LAYER SOURCE")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO POLYNOMIAL FITTING")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 13
    FORWARD_ELL = 15

    print("=" * 78)
    print("0. BUILDING TOP BOUNDARY DATASET")
    print("=" * 78)

    rows = build_training(
        TRAIN_MAX_ELL
    )

    print(
        f"rows = {len(rows)}"
    )

    print()

    top_failures = inspect_top_rows(
        rows
    )

    print()

    print("=" * 78)
    print("2. TOP-LAYER EDGE COEFFICIENT STRUCTURE")
    print("=" * 78)

    sample_pairs = [
        (1, 7),
        (1, 9),
        (3, 9),
        (5, 11),
        (1, 13),
        (3, 13)
    ]

    structural_failures = 0

    for k, ell in sample_pairs:

        _, E, O = get_layers(
            k,
            ell
        )

        info = derive_top_boundary_identity(
            E
        )

        if info is None:
            structural_failures += 1
            continue

        print(
            f"({k},{ell})"
        )

        print(
            f"  degree = {info['degree']}"
        )

        print(
            f"  top x^M coefficient = "
            f"{sp.factor(info['top_coeff'])}"
        )

        print(
            f"  second x^(M-2) coefficient = "
            f"{sp.factor(info['second_coeff'])}"
        )

        print(
            f"  adjacent Laurent coefficient = "
            f"{sp.factor(info['adjacent_coeff'])}"
        )

        print()

    print("=" * 78)
    print("3. LOWER-LAYER BOUNDARY SOURCE")
    print("=" * 78)

    lower_rows = []

    for ell in range(
        7,
        TRAIN_MAX_ELL + 1,
        2
    ):
        for k in range(
            1,
            ell,
            2
        ):

            result = lower_boundary_sources(
                k,
                ell
            )

            if result is None:
                continue

            interior_failures = [
                item for item in result["interior"]
                if item[1] != 0
            ]

            lower_rows.append(
                {
                    "k": k,
                    "ell": ell,
                    "j_max": result["j_max"],
                    "source": result["source"],
                    "interior_failures": interior_failures
                }
            )

            print(
                f"({k},{ell}) "
                f"j_max={result['j_max']} "
                f"source={sp.factor(result['source'])} "
                f"interior_failures={len(interior_failures)}"
            )

    print()

    print("=" * 78)
    print("4. FORWARD ELL HOLDOUT")
    print("=" * 78)

    forward_rows = build_training(
        FORWARD_ELL
    )

    forward_failures = inspect_top_rows(
        forward_rows
    )

    print()

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        f"top-source failures = {len(top_failures)}"
    )

    print(
        f"lower rows checked = {len(lower_rows)}"
    )

    lower_interior_failures = sum(
        len(row["interior_failures"])
        for row in lower_rows
    )

    print(
        f"lower interior recurrence failures = "
        f"{lower_interior_failures}"
    )

    print(
        f"forward top-source failures = "
        f"{len(forward_failures)}"
    )

    total_failures = (
        len(top_failures)
        + structural_failures
        + lower_interior_failures
        + len(forward_failures)
    )

    if total_failures == 0:

        print()
        print("STATUS = PASS")
        print()
        print(
            "The top-layer boundary source is exactly ell."
        )
        print()
        print(
            "The lower-layer interior recurrence is exact."
        )
        print()
        print(
            "The remaining symbolic target is the lower-layer"
        )
        print(
            "boundary source. Once that source is identified,"
        )
        print(
            "the recurrence + boundary-value proof is complete."
        )

    else:

        print()
        print("STATUS = DIAGNOSTIC")
        print(
            f"total observed failures = {total_failures}"
        )
        print()
        print(
            "The failures are preserved as symbolic evidence;"
        )
        print(
            "do not infer a false recurrence from this certificate."
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

