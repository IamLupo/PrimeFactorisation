#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 134
SYMBOLIC TWO-LAYER TOP-EDGE DERIVATION

GOAL
----
Turn Experiment 133 into a symbolic derivation.

For

    j = ell - k - d - 1

and the weighted quotient

    deg(S)=1
    deg(N)=2,

Experiment 133 found:

    d odd:
        leading coefficient comes from the absolute top layer;

    d even:
        absolute top layer has the wrong parity,
        and the leading coefficient comes from exactly one lower layer.

This experiment does NOT build the full Newton tensor.

Instead it symbolically extracts only:

    H_top
    H_next

and compares their homogeneous Newton coefficient at Pi_j.

TARGET:
    derive the closed forms

        d odd:
            LC = (-1)^((d-1)/2) (k+ell)

        d even:
            LC = (-1)^(d/2)
                 (k+ell)/2
                 (k-ell+d+1)

The experiment also checks whether the relevant layer itself has
a compact factorized form in terms of k, ell and x.

NO FULL NEWTON TENSOR
NO C/D TENSOR
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
x, t = sp.symbols("x t")


# ============================================================================
# Exact finite kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact quotient
# ============================================================================

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

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric reduction
# ============================================================================

def symmetric_to_NS(expr):
    expr_sub = sp.expand(
        expr.subs(q, S - p)
    )

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        expr_sub,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.expand(
        sp.rem(poly, modulus).as_expr()
    )

    p_coeff = sp.expand(
        rem.coeff(p, 1)
    )

    if p_coeff != 0:
        raise ArithmeticError(
            "Symmetric conversion retained p-dependence: "
            f"{sp.factor(p_coeff)}"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


# ============================================================================
# Weighted layers
# ============================================================================

def safe_int(value):
    value = sp.sympify(value)

    if not value.is_Integer:
        raise ArithmeticError(
            f"Expected integer exponent, got {value}"
        )

    return int(value)


def weighted_degree_term(term):
    powers = sp.sympify(term).as_powers_dict()

    return (
        safe_int(powers.get(S, 0))
        + 2 * safe_int(powers.get(N, 0))
    )


def weighted_layers(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return {}

    layers = {}

    for term in sp.Add.make_args(expr):
        w = weighted_degree_term(term)

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
                N: t**2,
            }
        )
    )

    pt = sp.Poly(
        scaled,
        t,
        domain=sp.QQ.frac_field(x),
    )

    if pt.is_zero:
        return sp.Integer(0)

    return sp.expand(
        pt.LC()
    )


# ============================================================================
# Homogeneous Newton basis
# ============================================================================

def build_pi(max_j):
    Pi = [sp.Integer(2)]

    if max_j >= 1:
        Pi.append(x)

    for j in range(2, max_j + 1):
        Pi.append(
            sp.expand(
                x * Pi[j - 1] - Pi[j - 2]
            )
        )

    return Pi


def pi_coefficients(H, max_j):
    Pi = build_pi(max_j)

    work = sp.Poly(
        sp.expand(H),
        x,
        domain=sp.QQ,
    )

    coeffs = {}

    for j in range(max_j, 0, -1):
        c = sp.expand(
            work.coeff_monomial(x**j)
        )

        coeffs[j] = c

        if c != 0:
            work -= sp.Poly(
                sp.expand(c * Pi[j]),
                x,
                domain=sp.QQ,
            )

    residual = sp.expand(work.as_expr())

    if residual != 0:
        if residual.is_Rational:
            coeffs[0] = sp.simplify(
                residual / 2
            )
        else:
            raise ArithmeticError(
                f"Nonconstant Pi residual: {sp.factor(residual)}"
            )

    return coeffs


# ============================================================================
# Established target law
# ============================================================================

def expected_lead(k, ell, d):
    if d % 2 == 1:
        return sp.expand(
            (-1) ** ((d - 1) // 2)
            * (k + ell)
        )

    return sp.expand(
        (-1) ** (d // 2)
        * sp.Rational(k + ell, 2)
        * (k - ell + d + 1)
    )


# ============================================================================
# Symbolic layer extraction for one case
# ============================================================================

def analyze_case(k, ell, d):
    j = ell - k - d - 1

    if j < 1:
        return None

    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    layers = weighted_layers(Q)

    weights = sorted(
        layers,
        reverse=True
    )

    top_weight = weights[0]
    top_x = layer_to_x(
        layers[top_weight]
    )

    top_coeff = sp.expand(
        pi_coefficients(
            top_x,
            ell - 1
        ).get(j, 0)
    )

    if d % 2 == 1:
        selected_weight = top_weight
    else:
        selected_weight = top_weight - 1

    if selected_weight not in layers:
        raise ArithmeticError(
            f"Expected lower layer missing for "
            f"({k},{ell}), d={d}"
        )

    selected_x = layer_to_x(
        layers[selected_weight]
    )

    selected_coeff = sp.expand(
        pi_coefficients(
            selected_x,
            ell - 1
        ).get(j, 0)
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "j": j,
        "top_weight": top_weight,
        "top_x": top_x,
        "top_coeff": top_coeff,
        "selected_weight": selected_weight,
        "selected_x": selected_x,
        "selected_coeff": selected_coeff,
        "expected": expected_lead(
            k,
            ell,
            d
        ),
    }


# ============================================================================
# Small exact grid
# ============================================================================

def build_records(max_ell):
    records = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            for d in range(1, ell - k - 1):
                rec = analyze_case(
                    k,
                    ell,
                    d
                )

                if rec is not None:
                    records.append(rec)

    return records


# ============================================================================
# Empirical symbolic-form search
# ============================================================================

def factor_selected_examples(records):
    """
    This does not fit laws.
    It only prints factorized selected homogeneous layers so we can
    identify a common closed form.
    """

    print("=" * 78)
    print("1. SELECTED HOMOGENEOUS LAYER FACTORIZATION")
    print("=" * 78)

    for rec in records[:]:
        if rec["d"] > 6:
            continue

        print(
            f"({rec['k']},{rec['ell']}) "
            f"d={rec['d']} j={rec['j']}"
        )

        print(
            f"  selected weight = "
            f"{rec['selected_weight']}"
        )

        print(
            f"  selected H(x) = "
            f"{sp.factor(rec['selected_x'])}"
        )

        print(
            f"  Pi_j coefficient = "
            f"{sp.factor(rec['selected_coeff'])}"
        )

        print(
            f"  expected = "
            f"{sp.factor(rec['expected'])}"
        )

        print()


# ============================================================================
# Validation
# ============================================================================

def validate(records):
    top_fail = []
    selected_fail = []

    for rec in records:

        if rec["d"] % 2 == 1:
            if sp.expand(
                rec["top_coeff"]
                - rec["expected"]
            ) != 0:
                top_fail.append(rec)

        if sp.expand(
            rec["selected_coeff"]
            - rec["expected"]
        ) != 0:
            selected_fail.append(rec)

    print("=" * 78)
    print("2. TWO-LAYER LEADING-COEFFICIENT CERTIFICATE")
    print("=" * 78)

    print(
        f"records = {len(records)}"
    )

    print(
        f"odd-d top-layer failures = "
        f"{len(top_fail)}"
    )

    print(
        f"selected-layer failures = "
        f"{len(selected_fail)}"
    )

    if top_fail:
        print()
        print("FIRST ODD-d FAILURES")
        for rec in top_fail[:10]:
            print(
                f"({rec['k']},{rec['ell']}) "
                f"d={rec['d']} "
                f"j={rec['j']}"
            )
            print(
                f"  top coeff = "
                f"{rec['top_coeff']}"
            )
            print(
                f"  expected = "
                f"{rec['expected']}"
            )

    if selected_fail:
        print()
        print("FIRST SELECTED-LAYER FAILURES")
        for rec in selected_fail[:10]:
            print(
                f"({rec['k']},{rec['ell']}) "
                f"d={rec['d']} "
                f"j={rec['j']}"
            )
            print(
                f"  selected coeff = "
                f"{rec['selected_coeff']}"
            )
            print(
                f"  expected = "
                f"{rec['expected']}"
            )

    return top_fail, selected_fail


# ============================================================================
# Forward ell
# ============================================================================

def forward_holdout(ell):
    records = build_records(ell)

    failures = [
        rec for rec in records
        if sp.expand(
            rec["selected_coeff"]
            - rec["expected"]
        ) != 0
    ]

    print("=" * 78)
    print("3. FORWARD ELL HOLDOUT")
    print("=" * 78)

    print(
        f"ell={ell}"
    )

    print(
        f"rows = {len(records)}"
    )

    print(
        f"selected-layer failures = "
        f"{len(failures)}"
    )

    return failures


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 134")
    print("SYMBOLIC TWO-LAYER TOP-EDGE DERIVATION")
    print("TARGET j = ell-k-d-1")
    print("ODD d -> TOP LAYER")
    print("EVEN d -> ONE LOWER LAYER")
    print("NO FULL NEWTON TENSOR")
    print("NO C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 13
    FORWARD_ELL = 15

    print("=" * 78)
    print("0. BUILDING EXACT TWO-LAYER DATASET")
    print("=" * 78)

    records = build_records(
        TRAIN_MAX_ELL
    )

    print(
        f"records = {len(records)}"
    )

    print()

    factor_selected_examples(
        records
    )

    top_fail, selected_fail = validate(
        records
    )

    print()

    forward_fail = forward_holdout(
        FORWARD_ELL
    )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    if not top_fail and not selected_fail and not forward_fail:

        print("STATUS = PASS")
        print()
        print(
            "The observed Experiment-133 parity mechanism"
        )
        print(
            "is reproduced exactly by the two relevant"
        )
        print(
            "weighted homogeneous layers."
        )
        print()
        print(
            "The remaining task is now symbolic factorization"
        )
        print(
            "of those two layers themselves."
        )

    else:

        print("STATUS = FAIL")

        print(
            f"odd-d top failures = {len(top_fail)}"
        )

        print(
            f"selected-layer failures = {len(selected_fail)}"
        )

        print(
            f"forward failures = {len(forward_fail)}"
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

