#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 143
SYMBOLIC EDGE-SOURCE COMPRESSION

TOP LAYER:
    derive the boundary source from the OUTER TWO x-coefficients

LOWER LAYER:
    derive the boundary source from the OUTER THREE Laurent coefficients

TARGET:
    B_top   = ell
    B_lower = ell*(ell-3)/2

This experiment deliberately does NOT construct the full Newton tensor.

The top boundary source is computed from

    E(x) = a_M x^M + a_(M-2) x^(M-2) + ...

and

    [u^M] E(u+u^-1)     = a_M
    [u^(M-2)]E(u+u^-1) = M*a_M + a_(M-2).

The boundary residual is therefore

    B_top =
      [u^(M-2)]E + [u^M]E.

For the lower layer we compute the exact boundary combination

    B_lower =
      A_(j+4) + 2 A_(j+2) + A_j

directly from the outer Laurent coefficients.

The script then checks whether the resulting exact expressions are
independent of k and equal to the predicted ell-only formulas.

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
K, L = sp.symbols("K L", integer=True, positive=True)


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
        domain=sp.QQ.frac_field(q),
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(PF, PD)

    remainder = sp.expand(R.as_expr())

    if remainder != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = "
            f"{sp.factor(remainder)}"
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
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        expr,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    remainder = sp.expand(
        sp.rem(poly, modulus).as_expr()
    )

    p_coeff = sp.expand(
        remainder.coeff(p, 1)
    )

    if p_coeff != 0:
        raise ArithmeticError(
            "symmetric reduction retained p-dependence: "
            f"{sp.factor(p_coeff)}"
        )

    return sp.expand(
        remainder.coeff(p, 0)
    )


# ============================================================================
# Weighted layers
# ============================================================================

def integer_exponent(value):
    value = sp.sympify(value)

    if not value.is_Integer:
        raise ArithmeticError(
            f"non-integer exponent: {value}"
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

    return s_exp + 2*n_exp


def weighted_layers(expr):
    layers = {}

    for term in sp.Add.make_args(
        sp.expand(expr)
    ):
        weight = weighted_degree(term)

        layers[weight] = sp.expand(
            layers.get(weight, 0) + term
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
                N: t**2,
            }
        )
    )

    poly = sp.Poly(
        scaled,
        t,
        domain=sp.QQ.frac_field(x),
    )

    if poly.is_zero:
        return sp.Integer(0)

    return sp.expand(
        poly.LC()
    )


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
# Outer x-coefficients
# ============================================================================

def polynomial_degree(poly):
    poly = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.QQ,
    )

    return poly.degree()


def x_coefficient(poly, degree):
    poly = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.QQ,
    )

    if degree < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(x**degree)
    )


# ============================================================================
# Laurent coefficient via binomial transform
# ============================================================================

def laurent_coeff(poly, j):
    P = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.QQ,
    )

    total = sp.Integer(0)

    for (m,), coeff in P.terms():

        m = int(m)

        if m < j:
            continue

        if (m - j) % 2 != 0:
            continue

        r = (m - j)//2

        total += coeff * sp.binomial(
            m,
            r
        )

    return sp.expand(total)


# ============================================================================
# Admissible j ranges
# ============================================================================

def admissible_js(k, ell, d_parity):
    result = []

    for d in range(
        1,
        ell - k - 1
    ):
        if d % 2 != d_parity:
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
# TOP boundary source
# ============================================================================

def top_boundary_record(k, ell):
    _, E, _ = get_layers(
        k,
        ell
    )

    degree = polynomial_degree(E)

    aM = x_coefficient(
        E,
        degree
    )

    aM2 = x_coefficient(
        E,
        degree - 2
    )

    outer_A = sp.expand(aM)

    adjacent_A = sp.expand(
        degree*aM + aM2
    )

    source = sp.expand(
        outer_A + adjacent_A
    )

    # This equals A_M + A_(M-2).
    return {
        "degree": degree,
        "aM": aM,
        "aM2": aM2,
        "A_M": outer_A,
        "A_M2": adjacent_A,
        "source": source,
    }


# ============================================================================
# LOWER boundary source
# ============================================================================

def lower_boundary_record(k, ell):
    _, _, O = get_layers(
        k,
        ell
    )

    if O == 0:
        return None

    rows = admissible_js(
        k,
        ell,
        0,
    )

    if not rows:
        return None

    jmax = max(
        j for _, j in rows
    )

    A0 = laurent_coeff(
        O,
        jmax
    )

    A2 = laurent_coeff(
        O,
        jmax + 2
    )

    A4 = laurent_coeff(
        O,
        jmax + 4
    )

    source = sp.expand(
        A4 + 2*A2 + A0
    )

    return {
        "jmax": jmax,
        "A0": A0,
        "A2": A2,
        "A4": A4,
        "source": source,
    }


# ============================================================================
# Exact k-independence check
# ============================================================================

def check_constant_in_k(records):
    """
    For each fixed ell, verify that all k-values give the same source.
    """

    failures = []

    by_ell = {}

    for record in records:
        by_ell.setdefault(
            record["ell"],
            []
        ).append(record)

    for ell_value, rows in by_ell.items():

        reference = rows[0]["source"]

        for row in rows[1:]:
            if sp.expand(
                row["source"] - reference
            ) != 0:
                failures.append(
                    (
                        ell_value,
                        reference,
                        row["source"],
                    )
                )

    return failures


# ============================================================================
# Exact polynomial law in ell
# ============================================================================

def derive_ell_law(
    records,
    expected_degree,
):
    """
    Since sources are already exact rationals, interpolate only in ell
    after verifying k-independence.

    This is not used to discover the law blindly: the target degree is
    dictated by the boundary order.

    The resulting polynomial is then checked against EVERY record.
    """

    points = []

    for record in records:
        points.append(
            (
                sp.Integer(record["ell"]),
                sp.Rational(record["source"])
            )
        )

    unique = {}

    for ell_value, source in points:
        if ell_value in unique:
            if sp.expand(
                unique[ell_value] - source
            ) != 0:
                raise ArithmeticError(
                    f"inconsistent source at ell={ell_value}"
                )
        else:
            unique[ell_value] = source

    ordered = sorted(
        unique.items(),
        key=lambda z: z[0]
    )

    if len(ordered) < expected_degree + 1:
        return None

    subset = ordered[
        :expected_degree + 1
    ]

    poly = sp.interpolate(
        subset,
        L
    )

    poly = sp.expand(
        poly
    )

    failures = []

    for ell_value, source in ordered:

        residual = sp.expand(
            poly.subs(
                L,
                ell_value
            )
            - source
        )

        if residual != 0:
            failures.append(
                (
                    ell_value,
                    residual
                )
            )

    return {
        "law": poly,
        "failures": failures,
        "points": ordered,
    }


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 143")
    print("SYMBOLIC EDGE-SOURCE COMPRESSION")
    print("TOP SOURCE + LOWER SOURCE")
    print("=" * 78)
    print()

    TRAIN_ELL = [5, 7, 9, 11, 13]
    FORWARD_ELL = [15, 17]

    # ------------------------------------------------------------------------
    # TOP
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. TOP-LAYER OUTER-COEFFICIENT DERIVATION")
    print("=" * 78)

    top_train = []

    for ell_value in TRAIN_ELL:

        for k_value in range(
            1,
            ell_value,
            2
        ):

            rec = top_boundary_record(
                k_value,
                ell_value
            )

            top_train.append(
                {
                    "k": k_value,
                    "ell": ell_value,
                    **rec,
                }
            )

    for rec in top_train:

        expected = sp.Integer(
            rec["ell"]
        )

        residual = sp.expand(
            rec["source"] - expected
        )

        print(
            f"({rec['k']},{rec['ell']}) "
            f"M={rec['degree']} "
            f"aM={sp.factor(rec['aM'])} "
            f"aM2={sp.factor(rec['aM2'])} "
            f"source={sp.factor(rec['source'])} "
            f"residual={sp.factor(residual)}"
        )

    top_k_failures = check_constant_in_k(
        top_train
    )

    top_law = derive_ell_law(
        top_train,
        expected_degree=1
    )

    print()
    print(
        "TOP k-independence failures =",
        len(top_k_failures)
    )

    if top_law is not None:
        print(
            "TOP derived ell-law =",
            sp.factor(top_law["law"])
        )
        print(
            "TOP law failures =",
            len(top_law["failures"])
        )

    # ------------------------------------------------------------------------
    # LOWER
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. LOWER-LAYER BOUNDARY SOURCE")
    print("=" * 78)

    lower_train = []

    for ell_value in TRAIN_ELL[1:]:

        for k_value in range(
            1,
            ell_value,
            2
        ):

            rec = lower_boundary_record(
                k_value,
                ell_value
            )

            if rec is None:
                continue

            lower_train.append(
                {
                    "k": k_value,
                    "ell": ell_value,
                    **rec,
                }
            )

    for rec in lower_train:

        expected = sp.expand(
            sp.Integer(rec["ell"])
            * (sp.Integer(rec["ell"]) - 3)
            / 2
        )

        residual = sp.expand(
            rec["source"] - expected
        )

        print(
            f"({rec['k']},{rec['ell']}) "
            f"jmax={rec['jmax']} "
            f"source={sp.factor(rec['source'])} "
            f"expected={sp.factor(expected)} "
            f"residual={sp.factor(residual)}"
        )

    lower_k_failures = check_constant_in_k(
        lower_train
    )

    lower_law = derive_ell_law(
        lower_train,
        expected_degree=2
    )

    print()
    print(
        "LOWER k-independence failures =",
        len(lower_k_failures)
    )

    if lower_law is not None:
        print(
            "LOWER derived ell-law =",
            sp.factor(lower_law["law"])
        )
        print(
            "LOWER law failures =",
            len(lower_law["failures"])
        )

    # ------------------------------------------------------------------------
    # FORWARD
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    forward_top_failures = 0
    forward_lower_failures = 0

    for ell_value in FORWARD_ELL:

        expected_top = sp.Integer(
            ell_value
        )

        expected_lower = sp.Integer(
            ell_value*(ell_value - 3)//2
        )

        for k_value in range(
            1,
            ell_value,
            2
        ):

            top = top_boundary_record(
                k_value,
                ell_value
            )

            top_residual = sp.expand(
                top["source"] - expected_top
            )

            if top_residual != 0:
                forward_top_failures += 1

            lower = lower_boundary_record(
                k_value,
                ell_value
            )

            if lower is not None:

                lower_residual = sp.expand(
                    lower["source"] - expected_lower
                )

                if lower_residual != 0:
                    forward_lower_failures += 1

    print(
        "forward top-source failures =",
        forward_top_failures
    )

    print(
        "forward lower-source failures =",
        forward_lower_failures
    )

    # ------------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    top_law_failures = (
        0
        if top_law is None
        else len(top_law["failures"])
    )

    lower_law_failures = (
        0
        if lower_law is None
        else len(lower_law["failures"])
    )

    total_failures = (
        len(top_k_failures)
        + len(lower_k_failures)
        + top_law_failures
        + lower_law_failures
        + forward_top_failures
        + forward_lower_failures
    )

    if total_failures == 0:

        print(
            "STATUS = PASS"
        )
        print()
        print(
            "TOP SOURCE:"
        )
        print(
            "    B_top = ell"
        )
        print()
        print(
            "LOWER SOURCE:"
        )
        print(
            "    B_lower = ell*(ell-3)/2"
        )
        print()
        print(
            "The two boundary sources have now been reduced"
        )
        print(
            "to exact ell-only formulas."
        )
        print()
        print(
            "The remaining Experiment-144 task is purely"
        )
        print(
            "recurrence + boundary induction."
        )

    else:

        print(
            "STATUS = DIAGNOSTIC"
        )

        print(
            "total failures =",
            total_failures
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

