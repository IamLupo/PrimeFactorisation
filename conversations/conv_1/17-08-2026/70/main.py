#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 135
LAURENT / GENERATING-FUNCTION NEWTON COEFFICIENT EXTRACTION

GOAL
----
Exploit the exact parametrization

    x = u + u^(-1)

for the homogeneous Newton basis

    Pi_0 = 2
    Pi_j = x Pi_{j-1} - Pi_{j-2}

which gives

    Pi_j(u + u^(-1))
        = u^j + u^(-j),     j >= 1.

Therefore, if

    H(x) = sum_j a_j Pi_j(x),

then

    H(u + u^(-1))
        = a_0 + sum_{j>=1} a_j (u^j + u^(-j)),

and for j >= 1,

    a_j = [u^j] H(u + u^(-1)).

This experiment tests whether the Experiment-134 selected homogeneous
layers can be extracted directly by Laurent coefficients, eliminating
the recursive Newton decomposition entirely.

TARGET
------
For

    j = ell-k-d-1,

verify

    [u^j] E_{k,ell}(u + u^(-1))
        = (-1)^((d-1)/2) (k+ell)

for odd d,

and

    [u^j] O_{k,ell}(u + u^(-1))
        = (-1)^((d/2)+1) (k+ell) j / 2

for even d.

The even formula is the Experiment-129/130 law rewritten using

    k-ell+d+1 = -j.

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
x, u, t = sp.symbols("x u t")


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
# Symmetric conversion to S,N
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
            f"Symmetric reduction retained p-dependence: "
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


# ============================================================================
# Convert one weighted layer to x
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
# Pi basis sanity check
# ============================================================================

def build_pi(max_j):
    Pi = [sp.Integer(2)]

    if max_j >= 1:
        Pi.append(x)

    for j in range(2, max_j + 1):
        Pi.append(
            sp.expand(
                x * Pi[j - 1]
                - Pi[j - 2]
            )
        )

    return Pi


def verify_laurent_pi_identity(max_j):
    failures = []

    for j in range(1, max_j + 1):
        pi = build_pi(j)[j]

        lhs = sp.expand(
            pi.subs(
                x,
                u + 1/u
            )
        )

        rhs = sp.expand(
            u**j + u**(-j)
        )

        if sp.expand(lhs - rhs) != 0:
            failures.append(
                (
                    j,
                    sp.expand(lhs),
                    sp.expand(rhs)
                )
            )

    return failures


# ============================================================================
# Laurent coefficient extractor
# ============================================================================

def laurent_coefficient(expr, power):
    """
    Exact coefficient of u^power in a Laurent polynomial.

    We use Poly after shifting by a sufficiently large power of u.
    """

    expr = sp.expand(expr)

    terms = sp.Add.make_args(expr)

    min_pow = 0

    for term in terms:
        powers = term.as_powers_dict()

        exponent = sp.sympify(
            powers.get(u, 0)
        )

        if not exponent.is_Integer:
            raise ArithmeticError(
                f"Non-integer u exponent: {term}"
            )

        min_pow = min(
            min_pow,
            int(exponent)
        )

    shift = -min_pow if min_pow < 0 else 0

    shifted = sp.expand(
        expr * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ
    )

    wanted = power + shift

    if wanted < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(
            u**wanted
        )
    )


# ============================================================================
# Laurent extraction of Newton coefficient
# ============================================================================

def coefficient_via_laurent(Hx, j):
    transformed = sp.expand(
        Hx.subs(
            x,
            u + 1/u
        )
    )

    return sp.expand(
        laurent_coefficient(
            transformed,
            j
        )
    )


# ============================================================================
# Established leading law
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
# Analyze one case
# ============================================================================

def analyze_case(k, ell, d):

    j = ell - k - d - 1

    if j < 1:
        return None

    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    layers = weighted_layers(Q)

    top_weight = max(layers)

    if d % 2 == 1:
        selected_weight = top_weight
    else:
        selected_weight = top_weight - 1

    if selected_weight not in layers:
        raise ArithmeticError(
            f"Missing selected layer for "
            f"({k},{ell}) d={d}"
        )

    Hx = layer_to_x(
        layers[selected_weight]
    )

    coeff = coefficient_via_laurent(
        Hx,
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
        "weight": selected_weight,
        "H": Hx,
        "coefficient": coeff,
        "expected": expected,
        "ok": sp.expand(
            coeff - expected
        ) == 0,
    }


# ============================================================================
# Training grid
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
# Grouped sign / j-law test
# ============================================================================

def verify_closed_forms(records):

    odd_fail = []
    even_fail = []

    for rec in records:

        if rec["d"] % 2 == 1:
            if rec["coefficient"] != (
                (-1) ** ((rec["d"] - 1) // 2)
                * (rec["k"] + rec["ell"])
            ):
                odd_fail.append(rec)

        else:
            expected = (
                (-1) ** (rec["d"] // 2 + 1)
                * sp.Rational(
                    rec["k"] + rec["ell"],
                    2
                )
                * rec["j"]
            )

            if sp.expand(
                rec["coefficient"] - expected
            ) != 0:
                even_fail.append(rec)

    return odd_fail, even_fail


# ============================================================================
# Examples
# ============================================================================

def print_examples(records):

    wanted = {
        (1,5,1),
        (1,5,2),
        (1,7,2),
        (1,9,3),
        (1,9,4),
        (3,11,4),
        (1,13,6),
        (5,13,5),
    }

    print("=" * 78)
    print("1. LAURENT NEWTON-COEFFICIENT EXAMPLES")
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
            f"{rec['weight']}"
        )

        print(
            f"  H(x) = "
            f"{sp.factor(rec['H'])}"
        )

        print(
            f"  [u^{rec['j']}] H(u+u^-1) = "
            f"{sp.factor(rec['coefficient'])}"
        )

        print(
            f"  expected = "
            f"{sp.factor(rec['expected'])}"
        )

        print()


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ell):

    records = build_records(
        ell
    )

    failures = [
        rec for rec in records
        if not rec["ok"]
    ]

    print("=" * 78)
    print("2. FORWARD ELL HOLDOUT")
    print("=" * 78)

    print(
        f"ell={ell}"
    )

    print(
        f"rows = {len(records)}"
    )

    print(
        f"failures = {len(failures)}"
    )

    return failures


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 135")
    print("LAURENT / GENERATING-FUNCTION NEWTON EXTRACTION")
    print("TARGET j = ell-k-d-1")
    print("ODD d -> EVEN LAYER")
    print("EVEN d -> ODD LAYER")
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
    # Basis identity
    # ------------------------------------------------------------------

    print("=" * 78)
    print("0. LAURENT NEWTON BASIS SANITY")
    print("=" * 78)

    basis_failures = verify_laurent_pi_identity(
        TRAIN_MAX_ELL
    )

    print(
        f"basis checks = {TRAIN_MAX_ELL}"
    )

    print(
        f"basis failures = "
        f"{len(basis_failures)}"
    )

    if basis_failures:
        for item in basis_failures[:5]:
            print(
                f"j={item[0]}: "
                f"{item[1]} != {item[2]}"
            )

        raise ArithmeticError(
            "Laurent Newton basis identity failed."
        )

    print(
        "STATUS = PASS"
    )
    print()

    # ------------------------------------------------------------------
    # Dataset
    # ------------------------------------------------------------------

    print("=" * 78)
    print("BUILDING LAURENT EXTRACTION DATASET")
    print("=" * 78)

    records = build_records(
        TRAIN_MAX_ELL
    )

    print(
        f"records = {len(records)}"
    )

    print()

    # ------------------------------------------------------------------
    # Closed-form validation
    # ------------------------------------------------------------------

    odd_fail, even_fail = verify_closed_forms(
        records
    )

    print("=" * 78)
    print("1. EXACT LAURENT EXTRACTION CERTIFICATE")
    print("=" * 78)

    print(
        f"records = {len(records)}"
    )

    print(
        f"odd-d failures = {len(odd_fail)}"
    )

    print(
        f"even-d failures = {len(even_fail)}"
    )

    print(
        f"total failures = "
        f"{len(odd_fail) + len(even_fail)}"
    )

    print()

    print_examples(
        records
    )

    # ------------------------------------------------------------------
    # Forward holdout
    # ------------------------------------------------------------------

    forward_failures = forward_holdout(
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
        len(odd_fail)
        + len(even_fail)
        + len(forward_failures)
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The Newton coefficient extraction is reproduced"
        )
        print(
            "directly as a Laurent coefficient:"
        )
        print()
        print(
            "    [Pi_j] H(x)"
        )
        print(
            "        ="
        )
        print(
            "    [u^j] H(u + u^(-1))."
        )
        print()
        print(
            "This removes recursive Newton-basis decomposition"
        )
        print(
            "from the leading-edge calculation."
        )
        print()
        print(
            "The next symbolic target is to derive closed forms"
        )
        print(
            "for the selected layer H(x) itself."
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

