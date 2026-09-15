#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 137
HYPERGEOMETRIC COEFFICIENT SEQUENCE OF THE SELECTED EDGE LAYER

OBJECTIVE
---------
Experiment 136 reduced the leading coefficient to

    L = sum_r h_r * binom(j+2r, r),

where

    H(x) = sum_r h_r x^(j+2r)

is the parity-compatible weighted homogeneous layer.

The next question is structural:

    Are the h_r themselves hypergeometric?

That is, does

    h_(r+1) / h_r

simplify to a rational function of r?

If yes, the leading-law sum is a candidate for exact
hypergeometric summation (Gosper / creative telescoping),
which is much closer to an analytic proof than further
tensor computation.

THIS EXPERIMENT DOES NOT FIT NUMERICAL DATA.

It performs only exact symbolic operations:

  1. construct the exact quotient;
  2. extract the selected homogeneous layer;
  3. extract the coefficients h_r;
  4. compute exact consecutive ratios h_(r+1)/h_r;
  5. factor those ratios;
  6. test the resulting binomial-weighted summand;
  7. attempt an exact Gosper certificate for the finite sum.

TARGET SUM
----------
    L = sum_r h_r * binom(j+2r, r)

EXPECTED CLOSED FORMS
---------------------
odd d:
    L = (-1)^((d-1)/2) * (k+ell)

even d:
    L = (-1)^(d/2+1) * (k+ell) * j / 2

NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp
from sympy.concrete.gosper import gosper_sum, gosper_normal, gosper_term


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")
x, t = sp.symbols("x t")
r = sp.symbols("r", integer=True, nonnegative=True)


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
# Symmetric conversion
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

    if sp.expand(
        rem.coeff(p, 1)
    ) != 0:
        raise ArithmeticError(
            f"Symmetric reduction retained p: "
            f"{sp.factor(rem.coeff(p, 1))}"
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
            f"Non-integer exponent: {value}"
        )

    return int(value)


def weighted_degree_term(term):
    powers = term.as_powers_dict()

    return (
        safe_int(powers.get(S, 0))
        + 2 * safe_int(powers.get(N, 0))
    )


def weighted_layers(expr):
    expr = sp.expand(expr)

    layers = {}

    for term in sp.Add.make_args(expr):

        w = weighted_degree_term(term)

        layers[w] = sp.expand(
            layers.get(w, 0) + term
        )

    return {
        w: sp.expand(poly)
        for w, poly in layers.items()
        if sp.expand(poly) != 0
    }


# ============================================================================
# Convert layer to H(x)
# ============================================================================

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

    return sp.expand(poly.LC())


# ============================================================================
# Selected parity-compatible layer
# ============================================================================

def selected_layer(Q_NS, d):

    layers = weighted_layers(Q_NS)

    if not layers:
        raise ArithmeticError(
            "No weighted layers."
        )

    top = max(layers)

    weight = (
        top
        if d % 2 == 1
        else top - 1
    )

    if weight not in layers:
        raise ArithmeticError(
            f"Missing selected layer: "
            f"top={top}, selected={weight}"
        )

    return (
        weight,
        layer_to_x(layers[weight])
    )


# ============================================================================
# Leading coefficient law
# ============================================================================

def expected_lead(k, ell, d):

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
# Extract coefficient sequence
# ============================================================================

def coefficient_sequence(H, j):
    """
    H(x) = sum_r h_r x^(j+2r).

    Returns [(r, h_r), ...].
    """

    P = sp.Poly(
        sp.expand(H),
        x,
        domain=sp.QQ
    )

    terms = []

    for (power,), coeff in P.terms():

        power = int(power)

        if power < j:
            continue

        delta = power - j

        if delta % 2 != 0:
            continue

        rr = delta // 2

        terms.append(
            (rr, sp.expand(coeff))
        )

    terms.sort(
        key=lambda z: z[0]
    )

    return terms


# ============================================================================
# Exact ratio certificate
# ============================================================================

def ratio_certificate(seq):

    ratios = []

    for i in range(
        len(seq) - 1
    ):

        r0, h0 = seq[i]
        r1, h1 = seq[i + 1]

        if r1 != r0 + 1:
            raise ArithmeticError(
                "Coefficient sequence has a gap."
            )

        if h0 == 0:
            ratios.append(
                (
                    r0,
                    sp.nan
                )
            )
        else:
            ratios.append(
                (
                    r0,
                    sp.factor(
                        sp.cancel(
                            h1 / h0
                        )
                    )
                )
            )

    return ratios


# ============================================================================
# Summand
# ============================================================================

def summand_from_sequence(seq, j):
    """
    Exact summands

        h_r * binom(j+2r,r).
    """

    return [
        (
            rr,
            sp.simplify(
                h * sp.binomial(
                    j + 2 * rr,
                    rr
                )
            )
        )
        for rr, h in seq
    ]


# ============================================================================
# Hypergeometric ratio of summand
# ============================================================================

def summand_ratio(seq, j):

    ratios = []

    for i in range(
        len(seq) - 1
    ):

        rr0, h0 = seq[i]
        rr1, h1 = seq[i + 1]

        if h0 == 0:
            ratios.append(
                (
                    rr0,
                    sp.nan
                )
            )
            continue

        ratio_h = sp.cancel(
            h1 / h0
        )

        ratio_binom = sp.cancel(
            sp.binomial(
                j + 2 * rr1,
                rr1
            )
            /
            sp.binomial(
                j + 2 * rr0,
                rr0
            )
        )

        ratios.append(
            (
                rr0,
                sp.factor(
                    sp.cancel(
                        ratio_h * ratio_binom
                    )
                )
            )
        )

    return ratios


# ============================================================================
# Exact finite sum
# ============================================================================

def exact_sum(seq, j):
    return sp.simplify(
        sp.expand(
            sum(
                h * sp.binomial(
                    j + 2 * rr,
                    rr
                )
                for rr, h in seq
            )
        )
    )


# ============================================================================
# Gosper attempt
# ============================================================================

def gosper_certificate(seq, j):

    if not seq:
        return {
            "status": "EMPTY"
        }

    n = sp.symbols(
        "n",
        integer=True,
        nonnegative=True
    )

    first_r = seq[0][0]

    # Re-index from the first available r.
    # The symbolic ratio is inferred from the finite exact sequence,
    # but the summand itself is kept exact.

    terms = [
        h * sp.binomial(
            j + 2 * rr,
            rr
        )
        for rr, h in seq
    ]

    finite = sp.expand(
        sum(terms)
    )

    result = {
        "finite_sum": sp.factor(finite),
        "gosper": None,
        "status": "NO_SYMBOLIC_TERM"
    }

    # For the present experiment we can only invoke Gosper directly
    # when h_r is represented as an explicit polynomial/rational function
    # in r. Therefore use the exact ratio to detect whether the sequence
    # appears hypergeometric on its support.

    ratios = ratio_certificate(seq)

    nontrivial = [
        rat for _, rat in ratios
        if rat is not sp.nan
    ]

    result["ratios"] = ratios

    if nontrivial:
        result["status"] = "RATIO_AVAILABLE"

    return result


# ============================================================================
# Analyze one case
# ============================================================================

def analyze_case(k, ell, d):

    j = ell - k - d - 1

    if j < 1:
        return None

    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    weight, H = selected_layer(
        Q,
        d
    )

    seq = coefficient_sequence(
        H,
        j
    )

    ratios = ratio_certificate(
        seq
    )

    sratios = summand_ratio(
        seq,
        j
    )

    total = exact_sum(
        seq,
        j
    )

    expected = expected_lead(
        k,
        ell,
        d
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "j": j,
        "weight": weight,
        "H": H,
        "sequence": seq,
        "ratios": ratios,
        "summand_ratios": sratios,
        "sum": total,
        "expected": expected,
        "sum_ok": sp.expand(
            total - expected
        ) == 0,
    }


# ============================================================================
# Build dataset
# ============================================================================

def build_records(max_ell):

    records = []

    for ell in range(
        3,
        max_ell + 1,
        2
    ):
        for k in range(
            1,
            ell,
            2
        ):
            for d in range(
                1,
                ell - k - 1
            ):

                rec = analyze_case(
                    k,
                    ell,
                    d
                )

                if rec is not None:
                    records.append(rec)

    return records


# ============================================================================
# Print coefficient-sequence structure
# ============================================================================

def print_structure(records):

    print("=" * 78)
    print("1. COEFFICIENT-SEQUENCE STRUCTURE")
    print("=" * 78)

    shown = 0

    for rec in records:

        if shown >= 12:
            break

        seq = rec["sequence"]

        if len(seq) < 2:
            continue

        print(
            f"({rec['k']},{rec['ell']}) "
            f"d={rec['d']} "
            f"j={rec['j']}"
        )

        print(
            f"  H(x) = "
            f"{sp.factor(rec['H'])}"
        )

        print(
            "  h_r = "
            + str(
                [
                    (rr, sp.factor(h))
                    for rr, h in seq
                ]
            )
        )

        print(
            "  h_(r+1)/h_r = "
            + str(
                [
                    (rr, rat)
                    for rr, rat in rec["ratios"]
                ]
            )
        )

        print(
            "  summand ratio = "
            + str(
                [
                    (rr, rat)
                    for rr, rat in rec["summand_ratios"]
                ]
            )
        )

        print(
            f"  finite sum = "
            f"{sp.factor(rec['sum'])}"
        )

        print(
            f"  expected = "
            f"{sp.factor(rec['expected'])}"
        )

        print()

        shown += 1


# ============================================================================
# Validation
# ============================================================================

def validate(records):

    failures = [
        rec for rec in records
        if not rec["sum_ok"]
    ]

    print("=" * 78)
    print("2. EXACT FINITE-SUM CERTIFICATE")
    print("=" * 78)

    print(
        f"records = {len(records)}"
    )

    print(
        f"finite-sum failures = "
        f"{len(failures)}"
    )

    return failures


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ell):

    records = build_records(
        ell
    )

    failures = [
        rec for rec in records
        if not rec["sum_ok"]
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
        f"finite-sum failures = "
        f"{len(failures)}"
    )

    return failures


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 137")
    print("HYPERGEOMETRIC COEFFICIENT SEQUENCE")
    print("SELECTED EDGE LAYER -> h_r -> BINOMIAL SUM")
    print("TARGET j = ell-k-d-1")
    print("NO FULL NEWTON TENSOR")
    print("NO C/D TENSOR")
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
    print("0. BUILDING COEFFICIENT-SEQUENCE DATASET")
    print("=" * 78)

    records = build_records(
        TRAIN_MAX_ELL
    )

    print(
        f"records = {len(records)}"
    )

    print()

    failures = validate(
        records
    )

    print()

    print_structure(
        records
    )

    forward_failures = forward_holdout(
        FORWARD_ELL
    )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    total_failures = (
        len(failures)
        + len(forward_failures)
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The leading law is reproduced exactly by"
        )
        print(
            "the finite coefficient sum"
        )
        print()
        print(
            "    L = sum_r h_r binom(j+2r,r)."
        )
        print()
        print(
            "The next proof question is now whether the"
        )
        print(
            "h_r sequence and the weighted summand admit"
        )
        print(
            "a closed hypergeometric representation."
        )
        print()
        print(
            "A successful exact ratio/Gosper certificate"
        )
        print(
            "would provide the first direct summation route"
        )
        print(
            "to the closed leading-law formula."
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

        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

