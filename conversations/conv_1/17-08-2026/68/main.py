#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 133
PARITY-ADJUSTED WEIGHTED NEWTON EDGE
FIRST HOMOGENEOUS LAYER WITH THE CORRECT NEWTON PARITY

TARGET:
    j = ell-k-d-1

WEIGHT:
    deg(S)=1
    deg(N)=2

KEY OBSERVATION FROM 132R
-------------------------
The absolute top weighted quotient is even in x.

Therefore:

    d odd  -> j even -> top layer can contribute directly

    d even -> j odd  -> top layer is parity-obstructed

For even d, the leading coefficient must come from the first lower
weighted layer containing the required odd parity.

GOAL
----
For each exact quotient Q(N,S):

  1. compute all weighted homogeneous layers;
  2. determine the parity of each layer in x;
  3. select the highest layer whose x-parity matches j;
  4. decompose that layer in the homogeneous Newton basis Pi_j;
  5. compare the Pi_j coefficient with the established Experiment-131
     leading coefficient.

This experiment tests the parity mechanism directly.

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
t, x = sp.symbols("t x")


# ============================================================================
# Kernel
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

    p_coeff = sp.expand(
        rem.coeff(p, 1)
    )

    if p_coeff != 0:
        raise ArithmeticError(
            f"Symmetric reduction retained p: {sp.factor(p_coeff)}"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


# ============================================================================
# Safe exponent handling
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

    s_pow = safe_int(
        powers.get(S, 0)
    )

    n_pow = safe_int(
        powers.get(N, 0)
    )

    return s_pow + 2 * n_pow


# ============================================================================
# Weighted layers
# ============================================================================

def weighted_layers(expr):
    """
    Return a dictionary

        weight -> homogeneous polynomial in S,N

    with all zero layers omitted.
    """

    expr = sp.expand(expr)

    if expr == 0:
        return {}

    layers = {}

    for term in sp.Add.make_args(expr):
        w = weighted_degree_term(term)
        layers[w] = sp.expand(
            layers.get(
                w,
                sp.Integer(0)
            )
            + term
        )

    return {
        w: sp.expand(poly)
        for w, poly in layers.items()
        if sp.expand(poly) != 0
    }


# ============================================================================
# Convert one layer to x under S=t*x, N=t^2
# ============================================================================

def layer_to_x(layer):
    """
    Every monomial in a weighted-homogeneous layer of weight W
    becomes t^W times a polynomial in x.
    """

    scaled = sp.expand(
        layer.subs(
            {
                S: t*x,
                N: t**2,
            }
        )
    )

    poly_t = sp.Poly(
        scaled,
        t,
        domain=sp.QQ.frac_field(x),
    )

    if poly_t.is_zero:
        return sp.Integer(0)

    return sp.expand(
        poly_t.LC()
    )


# ============================================================================
# Parity
# ============================================================================

def x_parity(poly):
    """
    Return:
        +1  even
        -1  odd
         0  mixed
    """

    poly = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.QQ
    )

    even = True
    odd = True

    for (power,), coeff in poly.terms():
        if coeff == 0:
            continue

        if power % 2 == 0:
            odd = False
        else:
            even = False

    if even:
        return +1

    if odd:
        return -1

    return 0


# ============================================================================
# Homogeneous Newton basis
# ============================================================================

def build_pi(max_j):
    Pi = [sp.Integer(2)]

    if max_j >= 1:
        Pi.append(x)

    for j in range(
        2,
        max_j + 1
    ):
        Pi.append(
            sp.expand(
                x * Pi[j - 1]
                - Pi[j - 2]
            )
        )

    return Pi


# ============================================================================
# Homogeneous Newton decomposition
# ============================================================================

def homogeneous_newton_coefficients(H, max_j):
    Pi = build_pi(
        max_j
    )

    work = sp.Poly(
        sp.expand(H),
        x,
        domain=sp.QQ,
    )

    coeffs = {}

    for j in range(
        max_j,
        0,
        -1
    ):

        c = sp.expand(
            work.coeff_monomial(
                x**j
            )
        )

        coeffs[j] = c

        if c != 0:
            work -= sp.Poly(
                sp.expand(
                    c * Pi[j]
                ),
                x,
                domain=sp.QQ,
            )

    residual = sp.expand(
        work.as_expr()
    )

    if residual != 0:

        if residual.is_Rational:
            coeffs[0] = sp.simplify(
                residual / 2
            )

            residual = sp.Integer(0)

        else:
            raise ArithmeticError(
                "Nonconstant residual after homogeneous Newton "
                f"decomposition: {sp.factor(residual)}"
            )

    return coeffs


# ============================================================================
# Established Experiment-131 laws
# ============================================================================

def expected_degree(k, d):
    return k + (d - 1) // 2


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
# Exact Newton coefficient
# ============================================================================

def build_exact_newton_basis(max_j):
    P = [sp.Integer(2)]

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


def exact_target(Q_NS, ell, j):
    P = build_exact_newton_basis(
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
            return c

        if c != 0:
            work -= sp.Poly(
                sp.expand(
                    c * P[idx]
                ),
                S,
                domain=sp.QQ.frac_field(N),
            )

    raise ArithmeticError(
        f"Could not extract Newton coefficient j={j}"
    )


# ============================================================================
# Analyze one case
# ============================================================================

def analyze_case(k, ell, d):

    j = ell - k - d - 1

    if j < 1:
        return None

    Qpq = quotient_Q(
        k,
        ell
    )

    QNS = symmetric_to_NS(
        Qpq
    )

    c = exact_target(
        QNS,
        ell,
        j
    )

    poly_N = sp.Poly(
        sp.expand(c),
        N,
        domain=sp.QQ,
    )

    degree = int(
        poly_N.degree()
    )

    lead = sp.expand(
        poly_N.LC()
    )

    expected_D = expected_degree(
        k,
        d
    )

    expected_L = expected_lead(
        k,
        ell,
        d
    )

    layers = weighted_layers(
        QNS
    )

    expected_parity = (
        +1 if j % 2 == 0
        else -1
    )

    compatible_layers = []

    for weight in sorted(
        layers.keys(),
        reverse=True
    ):

        Hx = layer_to_x(
            layers[weight]
        )

        parity = x_parity(
            Hx
        )

        if parity == expected_parity:
            compatible_layers.append(
                (
                    weight,
                    Hx
                )
            )

    if not compatible_layers:
        raise ArithmeticError(
            f"No parity-compatible layer found "
            f"for ({k},{ell}), d={d}, j={j}"
        )

    selected_weight, selected_H = (
        compatible_layers[0]
    )

    coeffs = homogeneous_newton_coefficients(
        selected_H,
        ell - 1
    )

    selected_coeff = sp.expand(
        coeffs.get(
            j,
            0
        )
    )

    # Also record the absolute top layer for comparison.
    absolute_weight = max(
        layers
    )

    absolute_H = layer_to_x(
        layers[absolute_weight]
    )

    absolute_coeff = sp.expand(
        homogeneous_newton_coefficients(
            absolute_H,
            ell - 1
        ).get(
            j,
            0
        )
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "j": j,
        "degree": degree,
        "lead": lead,
        "expected_D": expected_D,
        "expected_L": expected_L,
        "degree_ok": degree == expected_D,
        "direct_ok": sp.expand(
            lead - expected_L
        ) == 0,
        "absolute_weight": absolute_weight,
        "absolute_H": absolute_H,
        "absolute_coeff": absolute_coeff,
        "selected_weight": selected_weight,
        "selected_H": selected_H,
        "selected_coeff": selected_coeff,
        "selected_ok": sp.expand(
            selected_coeff - expected_L
        ) == 0,
        "parity_gap": (
            absolute_weight
            - selected_weight
        ),
        "expected_parity": expected_parity,
    }


# ============================================================================
# Grid
# ============================================================================

def build_grid(max_ell):

    rows = []

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

                row = analyze_case(
                    k,
                    ell,
                    d
                )

                if row is not None:
                    rows.append(
                        row
                    )

    return rows


# ============================================================================
# Validation
# ============================================================================

def validate(rows):

    degree_fail = [
        r for r in rows
        if not r["degree_ok"]
    ]

    direct_fail = [
        r for r in rows
        if not r["direct_ok"]
    ]

    selected_fail = [
        r for r in rows
        if not r["selected_ok"]
    ]

    print("=" * 78)
    print("1. PARITY-ADJUSTED TOP-LAYER CERTIFICATE")
    print("=" * 78)

    print(
        f"rows = {len(rows)}"
    )

    print(
        f"degree failures = "
        f"{len(degree_fail)}"
    )

    print(
        f"direct leading failures = "
        f"{len(direct_fail)}"
    )

    print(
        f"parity-adjusted layer failures = "
        f"{len(selected_fail)}"
    )

    gaps = sorted(
        {
            r["parity_gap"]
            for r in rows
        }
    )

    print(
        f"observed weighted-layer gaps = "
        f"{gaps}"
    )

    by_d = {}

    for r in rows:
        d = r["d"]
        by_d.setdefault(
            d,
            []
        ).append(
            r["parity_gap"]
        )

    print()
    print("GAP BY d")
    print("-" * 78)

    for d in sorted(
        by_d
    ):
        vals = sorted(
            set(by_d[d])
        )

        print(
            f"d={d:2d}: gaps={vals}"
        )

    if selected_fail:

        print()
        print(
            "FIRST PARITY-ADJUSTED FAILURES"
        )
        print("-" * 78)

        for r in selected_fail[:12]:

            print(
                f"({r['k']},{r['ell']}) "
                f"d={r['d']} "
                f"j={r['j']}"
            )

            print(
                f"  selected weight = "
                f"{r['selected_weight']}"
            )

            print(
                f"  selected H(x) = "
                f"{sp.factor(r['selected_H'])}"
            )

            print(
                f"  selected Pi_j coeff = "
                f"{r['selected_coeff']}"
            )

            print(
                f"  expected = "
                f"{r['expected_L']}"
            )

            print()

    return (
        degree_fail,
        direct_fail,
        selected_fail
    )


# ============================================================================
# Examples
# ============================================================================

def print_examples(rows):

    wanted = {
        (1, 5, 2),
        (1, 7, 2),
        (1, 7, 4),
        (3, 7, 2),
        (1, 9, 2),
        (1, 9, 4),
        (1, 9, 6),
        (3, 11, 4),
        (1, 13, 6),
        (3, 7, 1),
        (1, 9, 3),
    }

    print("=" * 78)
    print(
        "2. PARITY-ADJUSTED LAYER EXAMPLES"
    )
    print("=" * 78)

    for r in rows:

        key = (
            r["k"],
            r["ell"],
            r["d"]
        )

        if key not in wanted:
            continue

        print(
            f"({r['k']},{r['ell']}) "
            f"d={r['d']} "
            f"j={r['j']}"
        )

        print(
            f"  absolute weight = "
            f"{r['absolute_weight']}"
        )

        print(
            f"  absolute Pi_j coeff = "
            f"{r['absolute_coeff']}"
        )

        print(
            f"  selected weight = "
            f"{r['selected_weight']}"
        )

        print(
            f"  parity gap = "
            f"{r['parity_gap']}"
        )

        print(
            f"  selected H(x) = "
            f"{sp.factor(r['selected_H'])}"
        )

        print(
            f"  selected Pi_j coeff = "
            f"{sp.factor(r['selected_coeff'])}"
        )

        print(
            f"  expected LC = "
            f"{r['expected_L']}"
        )

        print(
            f"  exact LC = "
            f"{r['lead']}"
        )

        print()


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ells):

    rows = build_grid_for_ells(
        ells
    )

    degree_fail = sum(
        not r["degree_ok"]
        for r in rows
    )

    direct_fail = sum(
        not r["direct_ok"]
        for r in rows
    )

    selected_fail = sum(
        not r["selected_ok"]
        for r in rows
    )

    print("=" * 78)
    print(
        "3. FORWARD ELL HOLDOUT"
    )
    print("=" * 78)

    print(
        f"forward rows = {len(rows)}"
    )

    print(
        f"degree failures = {degree_fail}"
    )

    print(
        f"direct leading failures = {direct_fail}"
    )

    print(
        f"parity-adjusted layer failures = "
        f"{selected_fail}"
    )

    return rows


def build_grid_for_ells(ells):

    rows = []

    for ell in ells:

        for k in range(
            1,
            ell,
            2
        ):

            for d in range(
                1,
                ell - k - 1
            ):

                row = analyze_case(
                    k,
                    ell,
                    d
                )

                if row is not None:
                    rows.append(
                        row
                    )

    return rows


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 133")
    print(
        "PARITY-ADJUSTED WEIGHTED NEWTON EDGE"
    )
    print(
        "FIRST HOMOGENEOUS LAYER WITH CORRECT NEWTON PARITY"
    )
    print(
        "TARGET j = ell-k-d-1"
    )
    print("NO FULL NEWTON TENSOR")
    print("NO C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 13
    FORWARD_ELL = [15]

    # ------------------------------------------------------------------
    # Quotient sanity
    # ------------------------------------------------------------------

    print("=" * 78)
    print("0. QUOTIENT SANITY")
    print("=" * 78)

    sanity = [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]

    sanity_fail = 0

    for k, ell in sanity:

        try:
            quotient_Q(
                k,
                ell
            )

            print(
                f"({k},{ell}) quotient = PASS"
            )

        except Exception as exc:

            sanity_fail += 1

            print(
                f"({k},{ell}) quotient = FAIL: "
                f"{type(exc).__name__}: {exc}"
            )

    print(
        f"quotient failures = "
        f"{sanity_fail}/{len(sanity)}"
    )

    if sanity_fail:
        raise ArithmeticError(
            "Quotient sanity failed."
        )

    print()

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    print("=" * 78)
    print("BUILDING PARITY-ADJUSTED DATASET")
    print("=" * 78)

    rows = build_grid(
        TRAIN_MAX_ELL
    )

    print(
        f"training rows = {len(rows)}"
    )

    print()

    degree_fail, direct_fail, selected_fail = validate(
        rows
    )

    print()

    print_examples(
        rows
    )

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------

    forward = forward_holdout(
        FORWARD_ELL
    )

    forward_degree_fail = sum(
        not r["degree_ok"]
        for r in forward
    )

    forward_direct_fail = sum(
        not r["direct_ok"]
        for r in forward
    )

    forward_selected_fail = sum(
        not r["selected_ok"]
        for r in forward
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
        len(degree_fail)
        + len(direct_fail)
        + len(selected_fail)
        + forward_degree_fail
        + forward_direct_fail
        + forward_selected_fail
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The parity obstruction identified in Experiment 132R"
        )
        print(
            "is resolved by selecting the highest weighted layer"
        )
        print(
            "whose x-parity matches the Newton index j."
        )
        print()
        print(
            "This supports the mechanism:"
        )
        print(
            "  odd d  -> absolute top layer contributes"
        )
        print(
            "  even d -> first lower parity-compatible layer contributes"
        )
        print()
        print(
            "The selected homogeneous Pi_j coefficient reproduces"
        )
        print(
            "the exact Experiment-131 leading coefficient."
        )

    else:

        print("STATUS = FAIL")
        print(
            f"total failures = "
            f"{total_failures}"
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
