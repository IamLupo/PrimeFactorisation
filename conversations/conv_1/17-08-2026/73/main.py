#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 138
J-RECURRENCE CERTIFICATE FOR THE SELECTED LAURENT EDGE

PURPOSE
-------
Experiment 137 reduced the leading coefficient to

    L_j = [u^j] H(u + u^(-1))

where H is one of the two parity-compatible homogeneous layers.

For fixed (k,ell), the layer H does NOT change when d changes by 2.
Only

    j = ell-k-d-1

moves by 2.

Therefore we can eliminate d and study the coefficient sequence

    A_j = [u^j] H(u + u^(-1))

directly.

The experimentally predicted laws imply simple recurrences in j.

ODD-d BRANCH
------------
    A_j = sign(j) * (k+ell)

Hence, along the admissible parity class,

    A_(j+2) + A_j = 0.

EVEN-d BRANCH
-------------
    A_j = sign(j) * (k+ell) * j / 2

with alternating sign, hence

    A_(j+2) + A_j
        = -2 * sign-at-j * (k+ell)/2
        = alternating constant.

Equivalently, after removing the alternating sign,

    B_j = (-1)^((j-j0)/2) A_j

should be linear in j.

This experiment does NOT search for a closed h_r formula yet.
Instead it asks a narrower symbolic question:

    Can the final Laurent coefficient be characterized by a
    first-order recurrence in j?

If yes, the final proof may reduce to:
    1. one boundary coefficient;
    2. one exact j-recurrence.

This is potentially much cheaper than deriving the complete h_r sequence.

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
x, t = sp.symbols("x t")


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

    if sp.expand(rem.coeff(p, 1)) != 0:
        raise ArithmeticError(
            f"Symmetric reduction retained p-dependence: "
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
        w: sp.expand(v)
        for w, v in layers.items()
        if sp.expand(v) != 0
    }


# ============================================================================
# Layer -> x
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
# Get both parity layers for fixed (k,ell)
# ============================================================================

def get_pair_layers(k, ell):
    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    layers = weighted_layers(Q)

    if not layers:
        raise ArithmeticError(
            f"No weighted layers for ({k},{ell})"
        )

    top_weight = max(layers)

    even_layer = layer_to_x(
        layers[top_weight]
    )

    odd_layer = sp.Integer(0)

    if top_weight - 1 in layers:
        odd_layer = layer_to_x(
            layers[top_weight - 1]
        )

    return {
        "Q": Q,
        "top_weight": top_weight,
        "even_layer": even_layer,
        "odd_layer": odd_layer,
    }


# ============================================================================
# Direct coefficient [u^j] H(u+u^-1)
# ============================================================================

def laurent_edge(Hx, j):

    P = sp.Poly(
        sp.expand(Hx),
        x,
        domain=sp.QQ,
    )

    total = sp.Integer(0)

    for (m,), coeff in P.terms():

        m = int(m)

        if m < j:
            continue

        delta = m - j

        if delta % 2 != 0:
            continue

        r = delta // 2

        total += coeff * sp.binomial(
            m,
            r
        )

    return sp.expand(total)


# ============================================================================
# Admissible j values
# ============================================================================

def admissible_js(k, ell, parity):
    """
    parity = 1 for odd d
    parity = 0 for even d
    """

    values = []

    for d in range(
        1,
        ell - k - 1
    ):

        if d % 2 != parity:
            continue

        j = ell - k - d - 1

        if j >= 1:
            values.append(
                (d, j)
            )

    return values


# ============================================================================
# Predicted coefficient
# ============================================================================

def expected_coefficient(k, ell, d, j):

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
# Analyze one fixed pair
# ============================================================================

def analyze_pair(k, ell):

    data = get_pair_layers(
        k,
        ell
    )

    even_layer = data["even_layer"]
    odd_layer = data["odd_layer"]

    odd_rows = []
    even_rows = []

    # ------------------------------------------------------------
    # odd d -> top/even layer
    # ------------------------------------------------------------

    odd_js = admissible_js(
        k,
        ell,
        1
    )

    for d, j in odd_js:

        A = laurent_edge(
            even_layer,
            j
        )

        expected = expected_coefficient(
            k,
            ell,
            d,
            j
        )

        odd_rows.append(
            {
                "d": d,
                "j": j,
                "A": A,
                "expected": expected,
                "ok": sp.expand(
                    A - expected
                ) == 0,
            }
        )

    # ------------------------------------------------------------
    # even d -> lower/odd layer
    # ------------------------------------------------------------

    if odd_layer != 0:

        even_js = admissible_js(
            k,
            ell,
            0
        )

        for d, j in even_js:

            A = laurent_edge(
                odd_layer,
                j
            )

            expected = expected_coefficient(
                k,
                ell,
                d,
                j
            )

            even_rows.append(
                {
                    "d": d,
                    "j": j,
                    "A": A,
                    "expected": expected,
                    "ok": sp.expand(
                        A - expected
                    ) == 0,
                }
            )

    return {
        "k": k,
        "ell": ell,
        "even_layer": even_layer,
        "odd_layer": odd_layer,
        "odd_rows": odd_rows,
        "even_rows": even_rows,
    }


# ============================================================================
# Recurrence tests
# ============================================================================

def test_odd_recurrence(rows):
    """
    Ordered by increasing j.

    Expected:

        A_(j+2) + A_j = 0
    """

    rows = sorted(
        rows,
        key=lambda r: r["j"]
    )

    failures = []

    for a, b in zip(
        rows,
        rows[1:]
    ):

        if b["j"] != a["j"] + 2:
            continue

        residual = sp.expand(
            b["A"] + a["A"]
        )

        if residual != 0:
            failures.append(
                (
                    a["j"],
                    b["j"],
                    residual
                )
            )

    return failures


def test_even_linearized(rows):
    """
    For even d, define

        B_j = (-1)^((j-j_min)/2) A_j.

    Since j changes by 2, the sign is unambiguous.

    B_j should be linear in j.

    We test constant first differences.
    """

    rows = sorted(
        rows,
        key=lambda r: r["j"]
    )

    if len(rows) < 3:
        return [], []

    j0 = rows[0]["j"]

    transformed = []

    for row in rows:

        step = (row["j"] - j0) // 2

        sign = (-1) ** step

        B = sp.expand(
            sign * row["A"]
        )

        transformed.append(
            (
                row["j"],
                B
            )
        )

    diff_rows = []

    for a, b in zip(
        transformed,
        transformed[1:]
    ):

        diff_rows.append(
            (
                a[0],
                b[0],
                sp.expand(
                    b[1] - a[1]
                )
            )
        )

    failures = []

    if diff_rows:

        reference = diff_rows[0][2]

        for item in diff_rows[1:]:

            if sp.expand(
                item[2] - reference
            ) != 0:
                failures.append(
                    item
                )

    return failures, diff_rows


# ============================================================================
# Print pair summary
# ============================================================================

def print_pair_summary(result):

    k = result["k"]
    ell = result["ell"]

    print(
        f"({k},{ell})"
    )

    print(
        f"  even layer = "
        f"{sp.factor(result['even_layer'])}"
    )

    if result["odd_layer"] != 0:
        print(
            f"  odd layer  = "
            f"{sp.factor(result['odd_layer'])}"
        )
    else:
        print(
            "  odd layer  = 0"
        )

    odd_fail = test_odd_recurrence(
        result["odd_rows"]
    )

    even_fail, diff_rows = test_even_linearized(
        result["even_rows"]
    )

    print(
        f"  odd-d rows = "
        f"{len(result['odd_rows'])} "
        f"recurrence_failures={len(odd_fail)}"
    )

    print(
        f"  even-d rows = "
        f"{len(result['even_rows'])} "
        f"linear_recurrence_failures={len(even_fail)}"
    )

    if diff_rows:
        print(
            "  even branch transformed differences:"
        )

        for item in diff_rows:
            print(
                f"    j={item[0]}->{item[1]} "
                f"delta={sp.factor(item[2])}"
            )

    print()


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 138")
    print("J-RECURRENCE CERTIFICATE FOR THE SELECTED LAURENT EDGE")
    print("FIXED (k,ell), VARY j BY STEPS OF 2")
    print("NO FULL NEWTON TENSOR")
    print("NO C/D TENSOR")
    print("NO POLYNOMIAL FITTING")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    # Keep this deliberately small.
    # The point is symbolic recurrence structure, not mass enumeration.
    PAIRS = [
        (1, 7),
        (1, 9),
        (1, 11),
        (3, 9),
        (3, 11),
        (5, 11),
        (1, 13),
        (3, 13),
    ]

    print("=" * 78)
    print("1. FIXED-PAIR RECURRENCE TEST")
    print("=" * 78)

    results = []

    total_identity_failures = 0
    total_odd_recurrence_failures = 0
    total_even_recurrence_failures = 0

    for k, ell in PAIRS:

        result = analyze_pair(
            k,
            ell
        )

        results.append(
            result
        )

        identity_failures = sum(
            not row["ok"]
            for row in result["odd_rows"]
        ) + sum(
            not row["ok"]
            for row in result["even_rows"]
        )

        odd_fail = test_odd_recurrence(
            result["odd_rows"]
        )

        even_fail, _ = test_even_linearized(
            result["even_rows"]
        )

        total_identity_failures += (
            identity_failures
        )

        total_odd_recurrence_failures += (
            len(odd_fail)
        )

        total_even_recurrence_failures += (
            len(even_fail)
        )

        print_pair_summary(
            result
        )

    print("=" * 78)
    print("2. GLOBAL RECURRENCE SUMMARY")
    print("=" * 78)

    print(
        f"pairs tested = {len(results)}"
    )

    print(
        f"direct coefficient-law failures = "
        f"{total_identity_failures}"
    )

    print(
        f"odd-d recurrence failures = "
        f"{total_odd_recurrence_failures}"
    )

    print(
        f"even-d linear recurrence failures = "
        f"{total_even_recurrence_failures}"
    )

    print()

    print("=" * 78)
    print("3. BOUNDARY-VALUE CHECK")
    print("=" * 78)

    boundary_failures = 0

    for result in results:

        for rows in (
            result["odd_rows"],
            result["even_rows"],
        ):

            if not rows:
                continue

            first = min(
                rows,
                key=lambda z: z["j"]
            )

            if not first["ok"]:
                boundary_failures += 1

    print(
        f"boundary failures = "
        f"{boundary_failures}"
    )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    total_failures = (
        total_identity_failures
        + total_odd_recurrence_failures
        + total_even_recurrence_failures
        + boundary_failures
    )

    if total_failures == 0:

        print(
            "STATUS = PASS"
        )

        print()
        print(
            "For every tested fixed (k,ell), the selected"
        )
        print(
            "Laurent coefficient satisfies the expected"
        )
        print(
            "j-step recurrence."
        )

        print()
        print(
            "This reduces the symbolic proof target to:"
        )

        print(
            "  1. derive the recurrence in j;"
        )

        print(
            "  2. establish one boundary coefficient;"
        )

        print(
            "  3. conclude the closed leading law."
        )

        print()
        print(
            "This is intentionally narrower than another"
        )
        print(
            "hypergeometric search: it tests whether the"
        )
        print(
            "final coefficient can be characterized by"
        )
        print(
            "a recurrence without deriving every h_r."
        )

    else:

        print(
            "STATUS = FAIL"
        )

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

