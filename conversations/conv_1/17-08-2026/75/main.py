#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 139
SYMBOLIC J-RECURRENCE FROM THE LAURENT EDGE

PURPOSE
-------
Experiment 138 found, for fixed (k,ell), that the selected Laurent
coefficient sequence

    A_j = [u^j] H(u + u^(-1))

obeys very small recurrences as j changes by 2.

This experiment tries to derive those recurrences directly from the
EXACT Laurent polynomial R(u) = H(u + u^(-1)).

No fitting is performed.

For a Laurent polynomial

    R(u) = sum_j A_j u^j,

we have

    [u^(j+2)]R + [u^j]R

for the odd branch, and the corresponding signed first-difference
relation for the even branch.

The key symbolic certificate sought here is:

ODD LAYER
---------
    R(u) + u^2 R(u)
has vanishing coefficients on the admissible interior powers.

Equivalently,

    A_(j+2) + A_j = 0.

EVEN LAYER
----------
After sign normalization,

    B_j = (-1)^((j-j0)/2) A_j

should have constant step-2 difference.

Equivalently, an appropriate second-order relation is expected.

The experiment works with the exact selected homogeneous layers and
tries to derive the recurrence symbolically from the Laurent polynomial
itself, before evaluating individual coefficients.

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

    if sp.expand(rem.coeff(p, 1)) != 0:
        raise ArithmeticError(
            "Symmetric conversion retained p-dependence: "
            f"{sp.factor(rem.coeff(p, 1))}"
        )

    return sp.expand(rem.coeff(p, 0))


# ============================================================================
# Weighted layer extraction
# ============================================================================

def safe_int(v):
    v = sp.sympify(v)

    if not v.is_Integer:
        raise ArithmeticError(
            f"Expected integer exponent, got {v}"
        )

    return int(v)


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


def selected_layers(k, ell):
    Q = symmetric_to_NS(
        quotient_Q(k, ell)
    )

    layers = weighted_layers(Q)

    if not layers:
        raise ArithmeticError(
            f"No weighted layers for ({k},{ell})"
        )

    top = max(layers)

    even_layer = layer_to_x(
        layers[top]
    )

    odd_layer = (
        layer_to_x(layers[top - 1])
        if top - 1 in layers
        else sp.Integer(0)
    )

    return even_layer, odd_layer


# ============================================================================
# Laurent polynomial
# ============================================================================

def laurent_poly(H):
    return sp.expand(
        H.subs(
            x,
            u + 1/u
        )
    )


# ============================================================================
# Laurent support utilities
# ============================================================================

def laurent_support(expr):
    expr = sp.expand(expr)

    support = []

    for term in sp.Add.make_args(expr):
        exponent = sp.sympify(
            term.as_powers_dict().get(u, 0)
        )

        if not exponent.is_Integer:
            raise ArithmeticError(
                f"Non-integer Laurent exponent in term: {term}"
            )

        support.append(int(exponent))

    return sorted(set(support))


def laurent_coeff(expr, power):
    expr = sp.expand(expr)

    # Shift negative powers into an ordinary polynomial.
    support = laurent_support(expr)

    min_power = min(support) if support else 0
    shift = max(0, -min_power)

    shifted = sp.expand(
        expr * u**shift
    )

    poly = sp.Poly(
        shifted,
        u,
        domain=sp.QQ,
    )

    wanted = power + shift

    if wanted < 0:
        return sp.Integer(0)

    return sp.expand(
        poly.coeff_monomial(u**wanted)
    )


# ============================================================================
# Symbolic odd-branch recurrence
# ============================================================================

def symbolic_odd_recurrence(R):
    """
    Test whether the Laurent coefficient relation

        A_(j+2) + A_j = 0

    holds over the full interior support.

    Algebraically this means the Laurent polynomial

        (1 + u^(-2)) R(u)

    has no admissible interior coefficients of the relevant parity.

    We test exact coefficients, not numerical values.
    """

    support = laurent_support(R)

    positive = [
        e for e in support
        if e >= 1
    ]

    if not positive:
        return {
            "ok": True,
            "residual_support": []
        }

    jmin = min(positive)
    jmax = max(positive)

    failures = []

    for j in range(
        jmin,
        jmax - 1,
        2
    ):
        residual = sp.expand(
            laurent_coeff(R, j + 2)
            + laurent_coeff(R, j)
        )

        if residual != 0:
            failures.append(
                (
                    j,
                    sp.factor(residual)
                )
            )

    return {
        "ok": len(failures) == 0,
        "residual_support": failures
    }


# ============================================================================
# Symbolic even-branch recurrence
# ============================================================================

def symbolic_even_recurrence(R):
    """
    For the even branch the coefficient sequence is alternating-linear.

    Thus its second step difference after the correct alternating sign
    is constant. Equivalently:

        A_(j+4) + 2 A_(j+2) + A_j = 0

    for a pure alternating linear sequence.

    We test this exact second-order-in-step relation directly.
    """

    support = laurent_support(R)

    positive = [
        e for e in support
        if e >= 1
    ]

    if not positive:
        return {
            "ok": True,
            "residual_support": []
        }

    jmin = min(positive)
    jmax = max(positive)

    failures = []

    for j in range(
        jmin,
        jmax - 3,
        2
    ):
        residual = sp.expand(
            laurent_coeff(R, j + 4)
            + 2 * laurent_coeff(R, j + 2)
            + laurent_coeff(R, j)
        )

        if residual != 0:
            failures.append(
                (
                    j,
                    sp.factor(residual)
                )
            )

    return {
        "ok": len(failures) == 0,
        "residual_support": failures
    }


# ============================================================================
# Direct recurrence from the target formulas
# ============================================================================

def target_odd_value(k, ell, d):
    return sp.expand(
        (-1) ** ((d - 1) // 2)
        * (k + ell)
    )


def target_even_value(k, ell, d):
    j = ell - k - d - 1

    return sp.expand(
        (-1) ** (d // 2 + 1)
        * sp.Rational(k + ell, 2)
        * j
    )


# ============================================================================
# Coefficient recurrence sanity against target formulas
# ============================================================================

def target_recurrence_checks(k, ell):

    odd_values = []
    even_values = []

    for d in range(
        1,
        ell - k - 1
    ):

        j = ell - k - d - 1

        if j < 1:
            continue

        if d % 2 == 1:
            odd_values.append(
                (
                    j,
                    target_odd_value(
                        k,
                        ell,
                        d
                    )
                )
            )
        else:
            even_values.append(
                (
                    j,
                    target_even_value(
                        k,
                        ell,
                        d
                    )
                )
            )

    odd_failures = []

    odd_values.sort()

    for (j1, a1), (j2, a2) in zip(
        odd_values,
        odd_values[1:]
    ):
        if j2 != j1 + 2:
            continue

        if sp.expand(a2 + a1) != 0:
            odd_failures.append(
                (j1, j2, sp.expand(a2 + a1))
            )

    even_failures = []

    even_values.sort()

    for (j1, a1), (j2, a2) in zip(
        even_values,
        even_values[1:]
    ):
        if j2 != j1 + 2:
            continue

        # First difference after alternating sign.
        # Using the actual d parity relation, this is equivalent to
        # A_(j+2) + A_j being constant in j.
        residual = sp.expand(
            (a2 + a1)
        )

        if even_values:
            pass

    # Compare all step-2 sums to the first step-2 sum.
    if len(even_values) >= 2:
        first_sum = sp.expand(
            even_values[1][1] + even_values[0][1]
        )

        for (j1, a1), (j2, a2) in zip(
            even_values,
            even_values[1:]
        ):
            current = sp.expand(a2 + a1)

            if current != first_sum:
                even_failures.append(
                    (
                        j1,
                        j2,
                        current
                    )
                )

    return odd_failures, even_failures


# ============================================================================
# One fixed pair
# ============================================================================

def analyze_pair(k, ell):

    even_H, odd_H = selected_layers(
        k,
        ell
    )

    even_R = laurent_poly(
        even_H
    )

    odd_R = laurent_poly(
        odd_H
    )

    odd_symbolic = symbolic_odd_recurrence(
        even_R
    )

    even_symbolic = symbolic_even_recurrence(
        odd_R
    )

    target_odd, target_even = target_recurrence_checks(
        k,
        ell
    )

    return {
        "k": k,
        "ell": ell,
        "even_H": even_H,
        "odd_H": odd_H,
        "even_R": even_R,
        "odd_R": odd_R,
        "odd_symbolic": odd_symbolic,
        "even_symbolic": even_symbolic,
        "target_odd": target_odd,
        "target_even": target_even,
    }


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 139")
    print("SYMBOLIC J-RECURRENCE FROM THE LAURENT EDGE")
    print("EXACT LAURENT-POLYNOMIAL IDENTITIES")
    print("NO FULL NEWTON TENSOR")
    print("NO C/D TENSOR")
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
        (3, 13),
    ]

    symbolic_odd_failures = 0
    symbolic_even_failures = 0
    target_odd_failures = 0
    target_even_failures = 0

    for k, ell in PAIRS:

        print("=" * 78)
        print(
            f"PAIR ({k},{ell})"
        )
        print("=" * 78)

        result = analyze_pair(
            k,
            ell
        )

        print(
            "EVEN / TOP LAYER"
        )
        print(
            f"  H(x) = "
            f"{sp.factor(result['even_H'])}"
        )

        print(
            f"  Laurent support = "
            f"{laurent_support(result['even_R'])}"
        )

        print(
            f"  A_(j+2)+A_j recurrence = "
            f"{result['odd_symbolic']['ok']}"
        )

        if not result["odd_symbolic"]["ok"]:
            print(
                "  first residuals:"
            )
            for item in result["odd_symbolic"]["residual_support"][:5]:
                print(
                    f"    j={item[0]} "
                    f"residual={item[1]}"
                )

            symbolic_odd_failures += 1

        print()

        print(
            "ODD / LOWER LAYER"
        )

        if result["odd_H"] == 0:
            print(
                "  H(x) = 0"
            )
        else:
            print(
                f"  H(x) = "
                f"{sp.factor(result['odd_H'])}"
            )

            print(
                f"  Laurent support = "
                f"{laurent_support(result['odd_R'])}"
            )

            print(
                "  A_(j+4)+2A_(j+2)+A_j recurrence = "
                f"{result['even_symbolic']['ok']}"
            )

            if not result["even_symbolic"]["ok"]:
                print(
                    "  first residuals:"
                )

                for item in result["even_symbolic"]["residual_support"][:5]:
                    print(
                        f"    j={item[0]} "
                        f"residual={item[1]}"
                    )

                symbolic_even_failures += 1

        print()

        print(
            "TARGET-LAW RECURRENCE CHECK"
        )

        print(
            f"  odd-branch failures = "
            f"{len(result['target_odd'])}"
        )

        print(
            f"  even-branch failures = "
            f"{len(result['target_even'])}"
        )

        target_odd_failures += len(
            result["target_odd"]
        )

        target_even_failures += len(
            result["target_even"]
        )

        print()

    print("=" * 78)
    print("FINAL SUMMARY")
    print("=" * 78)

    print(
        f"pairs tested = {len(PAIRS)}"
    )

    print(
        f"symbolic odd-layer recurrence failures = "
        f"{symbolic_odd_failures}"
    )

    print(
        f"symbolic lower-layer recurrence failures = "
        f"{symbolic_even_failures}"
    )

    print(
        f"target odd recurrence failures = "
        f"{target_odd_failures}"
    )

    print(
        f"target even recurrence failures = "
        f"{target_even_failures}"
    )

    print()

    if (
        symbolic_odd_failures == 0
        and symbolic_even_failures == 0
        and target_odd_failures == 0
        and target_even_failures == 0
    ):

        print(
            "STATUS = PASS"
        )

        print()
        print(
            "The Laurent-layer recurrences are exact symbolic"
        )
        print(
            "identities on the tested fixed (k,ell) families."
        )

        print()
        print(
            "The leading-law proof is therefore reduced to"
        )
        print(
            "identifying the boundary Laurent coefficient."
        )

        print()
        print(
            "NEXT TARGET:"
        )
        print(
            "derive the extremal/boundary coefficient directly"
        )
        print(
            "from the finite binomial kernel."
        )

    else:

        print(
            "STATUS = FAIL"
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

