#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 136
DIRECT LAURENT EDGE FROM THE FINITE QUOTIENT
NO NEWTON BASIS DECOMPOSITION

CORE IDEA
---------
Experiment 135 established

    [Pi_j] H(x) = [u^j] H(u + u^(-1)).

For

    H(x) = sum_m h_m x^m,

we have

    [u^j] (u + u^(-1))^m
      = binom(m, (m-j)/2)

when

    m >= j
    m-j is even,

and zero otherwise.

Therefore

    [Pi_j] H
      =
    sum_{m >= j, m == j (mod 2)}
        h_m * binom(m, (m-j)/2).

This experiment removes BOTH:

    1. recursive Newton decomposition;
    2. explicit Laurent substitution.

The leading coefficient is computed directly from the ordinary
x-coefficients of the selected weighted homogeneous layer.

TARGET
------
j = ell-k-d-1

odd d:
    selected layer = top layer

even d:
    selected layer = one weighted layer below the top

EXPECTED LAW
------------
odd d:
    L = (-1)^((d-1)/2) (k+ell)

even d:
    L = (-1)^(d/2+1) * (k+ell) * j / 2

The experiment also independently compares:

    A) direct binomial-transform Laurent coefficient
    B) literal Laurent substitution
    C) exact Newton coefficient from Q(N,S)

No fitting is performed.

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


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")
x, u, t = sp.symbols("x u t")


# ============================================================================
# Original finite kernel
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
        expr.subs(
            q,
            S - p
        )
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
        sp.rem(
            poly,
            modulus
        ).as_expr()
    )

    p_coeff = sp.expand(
        rem.coeff(
            p,
            1
        )
    )

    if p_coeff != 0:
        raise ArithmeticError(
            f"Symmetric conversion retained p-dependence: "
            f"{sp.factor(p_coeff)}"
        )

    return sp.expand(
        rem.coeff(
            p,
            0
        )
    )


# ============================================================================
# Weighted degree
# ============================================================================

def safe_int(value):
    value = sp.sympify(value)

    if not value.is_Integer:
        raise ArithmeticError(
            f"Expected integer exponent, got {value}"
        )

    return int(value)


def weighted_degree_term(term):
    powers = sp.sympify(
        term
    ).as_powers_dict()

    return (
        safe_int(
            powers.get(
                S,
                0
            )
        )
        + 2 * safe_int(
            powers.get(
                N,
                0
            )
        )
    )


def weighted_layers(expr):
    expr = sp.expand(expr)

    if expr == 0:
        return {}

    layers = {}

    for term in sp.Add.make_args(expr):

        w = weighted_degree_term(
            term
        )

        layers[w] = sp.expand(
            layers.get(
                w,
                0
            ) + term
        )

    return {
        w: sp.expand(poly)
        for w, poly in layers.items()
        if sp.expand(poly) != 0
    }


# ============================================================================
# Convert a homogeneous weighted layer to H(x)
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

    return sp.expand(
        poly.LC()
    )


# ============================================================================
# Select the correct parity layer
# ============================================================================

def selected_layer(Q_NS, d):
    layers = weighted_layers(
        Q_NS
    )

    if not layers:
        raise ArithmeticError(
            "No weighted layers found."
        )

    top_weight = max(
        layers
    )

    if d % 2 == 1:
        selected_weight = top_weight
    else:
        selected_weight = top_weight - 1

    if selected_weight not in layers:
        raise ArithmeticError(
            f"Selected weighted layer missing: "
            f"top={top_weight}, selected={selected_weight}"
        )

    return (
        selected_weight,
        layer_to_x(
            layers[selected_weight]
        ),
        top_weight,
        layers
    )


# ============================================================================
# Direct Laurent coefficient from ordinary powers
# ============================================================================

def direct_binomial_laurent_coefficient(Hx, j):
    """
    If

        H(x) = sum_m h_m x^m,

    then

        [u^j] H(u+u^-1)
          =
        sum h_m * binom(m, (m-j)/2).

    This avoids constructing any negative powers of u.
    """

    poly = sp.Poly(
        sp.expand(Hx),
        x,
        domain=sp.QQ,
    )

    total = sp.Integer(0)

    for (m,), coeff in poly.terms():

        m = int(m)

        if m < j:
            continue

        if (m - j) % 2 != 0:
            continue

        r = (m - j) // 2

        total += coeff * sp.binomial(
            m,
            r
        )

    return sp.simplify(
        sp.expand(total)
    )


# ============================================================================
# Literal Laurent substitution
# ============================================================================

def literal_laurent_coefficient(Hx, j):
    expr_u = sp.expand(
        Hx.subs(
            x,
            u + 1/u
        )
    )

    terms = sp.Add.make_args(
        expr_u
    )

    min_power = 0

    for term in terms:

        power = sp.sympify(
            term.as_powers_dict().get(
                u,
                0
            )
        )

        if not power.is_Integer:
            raise ArithmeticError(
                f"Unexpected u exponent: {term}"
            )

        min_power = min(
            min_power,
            int(power)
        )

    shift = max(
        0,
        -min_power
    )

    shifted = sp.expand(
        expr_u * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ,
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
# Exact Newton basis, independent check only
# ============================================================================

def build_newton_basis(max_j):
    P = [
        sp.Integer(2)
    ]

    if max_j >= 1:
        P.append(S)

    for j in range(
        2,
        max_j + 1
    ):
        P.append(
            sp.expand(
                S * P[j - 1]
                - N * P[j - 2]
            )
        )

    return P


def exact_newton_coefficient(Q_NS, ell, j):
    P = build_newton_basis(
        ell - 1
    )

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    for idx in range(
        ell - 1,
        j - 1,
        -1
    ):

        c = sp.expand(
            work.coeff_monomial(
                S**idx
            )
        )

        if idx == j:
            return sp.expand(
                c
            )

        if c != 0:
            work -= sp.Poly(
                sp.expand(
                    c * P[idx]
                ),
                S,
                domain=sp.QQ.frac_field(N),
            )

    raise ArithmeticError(
        f"Could not extract exact Newton coefficient j={j}"
    )


# ============================================================================
# Established leading law
# ============================================================================

def expected_lead(k, ell, d):
    j = ell - k - d - 1

    if d % 2 == 1:
        return sp.expand(
            (-1) ** (
                (d - 1) // 2
            )
            * (k + ell)
        )

    return sp.expand(
        (-1) ** (
            d // 2 + 1
        )
        * sp.Rational(
            k + ell,
            2
        )
        * j
    )


# ============================================================================
# Analyze one case
# ============================================================================

def analyze_case(k, ell, d):
    j = ell - k - d - 1

    if j < 1:
        return None

    Q_NS = symmetric_to_NS(
        quotient_Q(
            k,
            ell
        )
    )

    selected_weight, Hx, top_weight, layers = selected_layer(
        Q_NS,
        d
    )

    direct_coeff = direct_binomial_laurent_coefficient(
        Hx,
        j
    )

    literal_coeff = literal_laurent_coefficient(
        Hx,
        j
    )

    exact_coeff = exact_newton_coefficient(
        Q_NS,
        ell,
        j
    )

    lead_poly = sp.Poly(
        exact_coeff,
        N,
        domain=sp.QQ,
    )

    exact_lead = sp.expand(
        lead_poly.LC()
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
        "top_weight": top_weight,
        "selected_weight": selected_weight,
        "weight_gap": (
            top_weight - selected_weight
        ),
        "H": Hx,
        "direct_coeff": direct_coeff,
        "literal_coeff": literal_coeff,
        "exact_coeff": exact_coeff,
        "exact_lead": exact_lead,
        "expected": expected,
        "direct_ok": sp.expand(
            direct_coeff - expected
        ) == 0,
        "literal_ok": sp.expand(
            literal_coeff - direct_coeff
        ) == 0,
        "exact_ok": sp.expand(
            exact_lead - direct_coeff
        ) == 0,
    }


# ============================================================================
# Build grid
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
                    records.append(
                        rec
                    )

    return records


# ============================================================================
# Summary by parity
# ============================================================================

def parity_summary(records):

    print("=" * 78)
    print("2. PARITY / WEIGHT SUMMARY")
    print("=" * 78)

    for parity in (
        "odd",
        "even"
    ):

        subset = [
            r for r in records
            if (
                r["d"] % 2 == 1
                if parity == "odd"
                else r["d"] % 2 == 0
            )
        ]

        gaps = sorted(
            set(
                r["weight_gap"]
                for r in subset
            )
        )

        print(
            f"{parity:>4} d: "
            f"rows={len(subset):3d} "
            f"weight_gaps={gaps}"
        )


# ============================================================================
# Examples
# ============================================================================

def print_examples(records):

    wanted = {
        (1, 5, 1),
        (1, 5, 2),
        (1, 7, 2),
        (1, 9, 3),
        (1, 9, 4),
        (3, 11, 4),
        (1, 13, 6),
        (5, 13, 5),
    }

    print("=" * 78)
    print("3. DIRECT LAURENT-EDGE EXAMPLES")
    print("=" * 78)

    for rec in records:

        key = (
            rec["k"],
            rec["ell"],
            rec["d"]
        )

        if key not in wanted:
            continue

        print(
            f"({rec['k']},{rec['ell']}) "
            f"d={rec['d']} "
            f"j={rec['j']}"
        )

        print(
            f"  selected weight = "
            f"{rec['selected_weight']}"
        )

        print(
            f"  H(x) = "
            f"{sp.factor(rec['H'])}"
        )

        print(
            f"  direct binomial coefficient = "
            f"{sp.factor(rec['direct_coeff'])}"
        )

        print(
            f"  literal Laurent coefficient = "
            f"{sp.factor(rec['literal_coeff'])}"
        )

        print(
            f"  exact Newton leading coefficient = "
            f"{sp.factor(rec['exact_lead'])}"
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

    direct_fail = [
        r for r in records
        if not r["direct_ok"]
    ]

    literal_fail = [
        r for r in records
        if not r["literal_ok"]
    ]

    exact_fail = [
        r for r in records
        if not r["exact_ok"]
    ]

    print("=" * 78)
    print("1. DIRECT BINOMIAL-LAURENT CERTIFICATE")
    print("=" * 78)

    print(
        f"records = {len(records)}"
    )

    print(
        f"direct-law failures = "
        f"{len(direct_fail)}"
    )

    print(
        f"literal-vs-direct failures = "
        f"{len(literal_fail)}"
    )

    print(
        f"exact-vs-direct failures = "
        f"{len(exact_fail)}"
    )

    if direct_fail:

        print()
        print(
            "FIRST DIRECT-LAW FAILURES"
        )
        print("-" * 78)

        for r in direct_fail[:10]:

            print(
                f"({r['k']},{r['ell']}) "
                f"d={r['d']} "
                f"j={r['j']}"
            )

            print(
                f"  direct = {r['direct_coeff']}"
            )

            print(
                f"  expected = {r['expected']}"
            )

            print()

    return (
        direct_fail,
        literal_fail,
        exact_fail
    )


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ell):

    records = build_records(
        ell
    )

    direct_fail = sum(
        not r["direct_ok"]
        for r in records
    )

    literal_fail = sum(
        not r["literal_ok"]
        for r in records
    )

    exact_fail = sum(
        not r["exact_ok"]
        for r in records
    )

    print("=" * 78)
    print("4. FORWARD ELL HOLDOUT")
    print("=" * 78)

    print(
        f"ell={ell}"
    )

    print(
        f"rows = {len(records)}"
    )

    print(
        f"direct-law failures = "
        f"{direct_fail}"
    )

    print(
        f"literal-vs-direct failures = "
        f"{literal_fail}"
    )

    print(
        f"exact-vs-direct failures = "
        f"{exact_fail}"
    )

    return (
        direct_fail,
        literal_fail,
        exact_fail
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 136")
    print("DIRECT LAURENT EDGE FROM THE FINITE QUOTIENT")
    print("NO NEWTON BASIS DECOMPOSITION")
    print("TARGET j = ell-k-d-1")
    print("ODD d -> TOP LAYER")
    print("EVEN d -> ONE LOWER LAYER")
    print("x = u + u^(-1)")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 13
    FORWARD_ELL = 15

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    print("=" * 78)
    print("0. BUILDING DIRECT LAURENT-EDGE DATASET")
    print("=" * 78)

    records = build_records(
        TRAIN_MAX_ELL
    )

    print(
        f"records = {len(records)}"
    )

    print()

    # ------------------------------------------------------------------
    # Validate
    # ------------------------------------------------------------------

    direct_fail, literal_fail, exact_fail = validate(
        records
    )

    print()

    parity_summary(
        records
    )

    print()

    print_examples(
        records
    )

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    forward = forward_holdout(
        FORWARD_ELL
    )

    # ------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    total_failures = (
        len(direct_fail)
        + len(literal_fail)
        + len(exact_fail)
        + sum(forward)
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The leading coefficient is obtained directly from"
        )
        print(
            "the ordinary x-coefficients of the selected homogeneous"
        )
        print(
            "layer using the binomial Laurent transform:"
        )
        print()
        print(
            "  [u^j](u+u^(-1))^m"
        )
        print(
            "      = binom(m,(m-j)/2)."
        )
        print()
        print(
            "Thus the leading-law calculation no longer requires"
        )
        print(
            "either Newton-basis recursion or explicit Laurent"
        )
        print(
            "polynomial expansion."
        )
        print()
        print(
            "The remaining symbolic problem is now narrowed to"
        )
        print(
            "a closed-form evaluation of the finite binomial sum"
        )
        print(
            "for the selected homogeneous layer."
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

