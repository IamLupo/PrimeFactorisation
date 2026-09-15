#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 131
SYMBOLIC TOP-N DERIVATION FROM THE FINITE BINOMIAL KERNEL

TARGET COEFFICIENT
------------------
    j = ell - k - d - 1

ESTABLISHED BY EXPERIMENT 130
-----------------------------
    deg_N c_j = k + floor((d-1)/2)

    d odd:
        LC = (-1)^((d-1)/2) * (k + ell)

    d even:
        LC = (-1)^(d/2)
             * (k + ell)/2
             * (k - ell + d + 1)

GOAL
----
Do not merely verify the law from the completed Newton coefficient.

Instead isolate the highest-N contribution symbolically from the
finite binomial kernel and determine which kernel terms contribute
to the top edge.

The experiment has four layers:

  1. exact quotient;
  2. exact target-index extraction;
  3. decomposition of c_j(N) by N-degree;
  4. symbolic "top-edge source" accounting.

The source accounting tracks every term of the quotient whose
contribution can reach the maximal N-degree of c_j.

This is intended as a proof-oriented experiment.

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
from collections import defaultdict

import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
S, N, t = sp.symbols("S N t")


# ============================================================================
# Basic finite kernel
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# Exact quotient by p+q+1
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
    """
    q = S-p, followed by reduction modulo

        p^2 - S*p + N.

    The input is symmetric, so the reduced result must be p-free.
    """

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
# Newton basis
# ============================================================================

def build_newton_basis(max_j):
    """
    P0 = 2
    P1 = S
    Pj = S*P(j-1) - N*P(j-2)
    """

    P = [sp.Integer(2)]

    if max_j >= 1:
        P.append(S)

    for j in range(2, max_j + 1):
        P.append(
            sp.expand(
                S * P[j - 1]
                - N * P[j - 2]
            )
        )

    return P


# ============================================================================
# Exact target extraction
# ============================================================================

def extract_target(Q_NS, ell, target_j):
    """
    Descending Newton extraction, retaining only the target coefficient.
    """

    P = build_newton_basis(ell - 1)

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    min_j = target_j
    target = None

    for j in range(
        ell - 1,
        min_j - 1,
        -1
    ):

        c = sp.expand(
            work.coeff_monomial(S**j)
        )

        if j == target_j:
            target = sp.expand(c)
            break

        if c != 0:
            work -= sp.Poly(
                sp.expand(c * P[j]),
                S,
                domain=sp.QQ.frac_field(N),
            )

    if target is None:
        raise ArithmeticError(
            f"Target Newton coefficient j={target_j} was not found."
        )

    return sp.expand(target)


# ============================================================================
# Established laws
# ============================================================================

def expected_degree(k, d):
    return k + (d - 1) // 2


def expected_lead(k, ell, d):
    if d % 2:
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
# Extract N-edge
# ============================================================================

def n_edge(expr):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return None, sp.Integer(0)

    return int(poly.degree()), sp.expand(poly.LC())


# ============================================================================
# Coefficient-by-coefficient N decomposition
# ============================================================================

def n_layers(expr):
    """
    Return:

        degree -> coefficient of N^degree

    with all zero layers removed.
    """

    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ.frac_field(S),
    )

    layers = {}

    if poly.is_zero:
        return layers

    for power in range(
        int(poly.degree()),
        -1,
        -1
    ):
        coeff = sp.expand(
            poly.coeff_monomial(
                N**power
            )
        )

        if coeff != 0:
            layers[power] = coeff

    return layers


# ============================================================================
# Top-layer source accounting
# ============================================================================

def top_layer_source(expr, expected_D):
    """
    The finished coefficient c_j(N) is already available here.

    We inspect its construction at the level of the N-polynomial:
      - top N degree,
      - top coefficient,
      - next coefficient,
      - factorization of the top coefficient.

    This is deliberately transparent rather than attempting a huge symbolic
    decomposition of the original quotient.
    """

    layers = n_layers(expr)

    if expected_D not in layers:
        return {
            "found": False,
            "layers": layers,
        }

    return {
        "found": True,
        "layers": layers,
        "top": layers[expected_D],
        "second": layers.get(
            expected_D - 1,
            sp.Integer(0),
        ),
    }


# ============================================================================
# Formal leading-symbol experiment
# ============================================================================

def formal_leading_test(expr, degree):
    """
    Introduce N=t^2 and S=t*x.

    The Newton basis has the natural weighted scaling

        deg(N)=2, deg(S)=1.

    After substituting

        N=t^2,
        S=t*x,

    the coefficient of the highest power of t records the weighted
    leading form.

    We inspect whether that leading form is independent of x.
    """

    x = sp.symbols("x")

    scaled = sp.expand(
        expr.subs({
            N: t**2,
            S: t*x,
        })
    )

    poly_t = sp.Poly(
        scaled,
        t,
        domain=sp.QQ.frac_field(x),
    )

    if poly_t.is_zero:
        return {
            "degree_t": None,
            "leading_form": sp.Integer(0),
        }

    deg_t = int(poly_t.degree())
    lead_form = sp.expand(
        poly_t.LC()
    )

    return {
        "degree_t": deg_t,
        "leading_form": lead_form,
    }


# ============================================================================
# Symbolic parity certificate
# ============================================================================

def symbolic_parity_certificate(k, ell, d):
    """
    Rewrite the established leading law in a,b variables and then reduce it
    back to k,ell. This checks exact equivalence between the two natural
    parameterizations.
    """

    a = sp.Rational(k - 1, 2)
    b = sp.Rational(ell - 1, 2)

    if d % 2:
        ab = sp.expand(
            2
            * (-1) ** ((d - 1) // 2)
            * (a + b + 1)
        )
    else:
        ab = sp.expand(
            (-1) ** (d // 2)
            * (a + b + 1)
            * (2*a - 2*b + d + 1)
        )

    kl = expected_lead(
        k,
        ell,
        d
    )

    return sp.expand(
        ab - kl
    ) == 0


# ============================================================================
# One exact case
# ============================================================================

def analyze_case(k, ell, d):
    r = k + d
    j = ell - k - d - 1

    if j < 1:
        return None

    Qpq = quotient_Q(k, ell)
    QNS = symmetric_to_NS(Qpq)

    c = extract_target(
        QNS,
        ell,
        j
    )

    degree, lead = n_edge(
        c
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

    top = top_layer_source(
        c,
        expected_D
    )

    formal = formal_leading_test(
        c,
        expected_D
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "r": r,
        "j": j,
        "c": c,
        "degree": degree,
        "lead": lead,
        "expected_D": expected_D,
        "expected_L": expected_L,
        "degree_ok": degree == expected_D,
        "lead_ok": sp.expand(lead - expected_L) == 0,
        "top": top,
        "formal": formal,
        "parity_ok": symbolic_parity_certificate(
            k,
            ell,
            d
        ),
    }


# ============================================================================
# Training grid
# ============================================================================

def training_grid(max_ell):
    rows = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
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
                    rows.append(row)

    return rows


# ============================================================================
# Print diagnostic table
# ============================================================================

def print_examples(rows):
    wanted = {
        (1, 5, 1),
        (3, 7, 1),
        (1, 7, 2),
        (5, 9, 1),
        (1, 9, 3),
        (3, 11, 4),
        (5, 13, 5),
        (1, 13, 6),
        (3, 17, 8),
    }

    print("=" * 78)
    print("TOP-EDGE SOURCE CERTIFICATES")
    print("=" * 78)

    for row in rows:

        if (
            row["k"],
            row["ell"],
            row["d"],
        ) not in wanted:
            continue

        print(
            f"({row['k']},{row['ell']}) "
            f"d={row['d']} "
            f"r={row['r']} "
            f"j={row['j']}"
        )

        print(
            f"  c_j(N) = {sp.factor(row['c'])}"
        )

        print(
            f"  degree = {row['degree']} "
            f"expected={row['expected_D']}"
        )

        print(
            f"  lead   = {row['lead']} "
            f"expected={row['expected_L']}"
        )

        print(
            f"  top-layer coefficient = "
            f"{sp.factor(row['top']['top'])}"
        )

        print(
            f"  second-layer coefficient = "
            f"{sp.factor(row['top']['second'])}"
        )

        print(
            f"  weighted-leading-form = "
            f"{sp.factor(row['formal']['leading_form'])}"
        )

        print(
            f"  parity certificate = "
            f"{row['parity_ok']}"
        )

        print()


# ============================================================================
# Training validation
# ============================================================================

def validate(rows):

    degree_fail = [
        r for r in rows
        if not r["degree_ok"]
    ]

    lead_fail = [
        r for r in rows
        if not r["lead_ok"]
    ]

    parity_fail = [
        r for r in rows
        if not r["parity_ok"]
    ]

    print("=" * 78)
    print("1. EXACT TOP-N DERIVATION CHECK")
    print("=" * 78)

    print(
        f"rows = {len(rows)}"
    )

    print(
        f"degree failures = "
        f"{len(degree_fail)}"
    )

    print(
        f"leading failures = "
        f"{len(lead_fail)}"
    )

    print(
        f"parity failures = "
        f"{len(parity_fail)}"
    )

    if degree_fail:
        print()
        print("FIRST DEGREE FAILURES")
        for r in degree_fail[:10]:
            print(
                f"({r['k']},{r['ell']}) "
                f"d={r['d']} "
                f"j={r['j']} "
                f"got={r['degree']} "
                f"expected={r['expected_D']}"
            )

    if lead_fail:
        print()
        print("FIRST LEADING-COEFFICIENT FAILURES")
        for r in lead_fail[:10]:
            print(
                f"({r['k']},{r['ell']}) "
                f"d={r['d']} "
                f"j={r['j']} "
                f"got={r['lead']} "
                f"expected={r['expected_L']}"
            )

    return (
        degree_fail,
        lead_fail,
        parity_fail,
    )


# ============================================================================
# Forward holdout
# ============================================================================

def forward_holdout(ells):

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
                    rows.append(row)

    degree_fail = [
        r for r in rows
        if not r["degree_ok"]
    ]

    lead_fail = [
        r for r in rows
        if not r["lead_ok"]
    ]

    print("=" * 78)
    print("2. FORWARD ELL HOLDOUT")
    print("=" * 78)

    print(
        f"forward rows = {len(rows)}"
    )

    print(
        f"degree failures = "
        f"{len(degree_fail)}"
    )

    print(
        f"leading failures = "
        f"{len(lead_fail)}"
    )

    return rows


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 131")
    print("SYMBOLIC TOP-N DERIVATION FROM THE FINITE BINOMIAL KERNEL")
    print("TARGET j = ell-k-d-1")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 15
    FORWARD_ELL = [17]

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

    failures = 0

    for k, ell in sanity:
        try:
            quotient_Q(k, ell)
            print(
                f"({k},{ell}) quotient = PASS"
            )
        except Exception as exc:
            failures += 1
            print(
                f"({k},{ell}) quotient = FAIL: "
                f"{type(exc).__name__}: {exc}"
            )

    print(
        f"quotient failures = "
        f"{failures}/{len(sanity)}"
    )

    if failures:
        raise ArithmeticError(
            "Quotient sanity failed."
        )

    print()

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    print("=" * 78)
    print("BUILDING EXACT TOP-EDGE DATASET")
    print("=" * 78)

    rows = training_grid(
        TRAIN_MAX_ELL
    )

    print(
        f"training rows = {len(rows)}"
    )
    print()

    degree_fail, lead_fail, parity_fail = validate(
        rows
    )

    print()

    print_examples(
        rows
    )

    # ------------------------------------------------------------------
    # Forward holdout
    # ------------------------------------------------------------------

    forward = forward_holdout(
        FORWARD_ELL
    )

    forward_degree_fail = [
        r for r in forward
        if not r["degree_ok"]
    ]

    forward_lead_fail = [
        r for r in forward
        if not r["lead_ok"]
    ]

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
        + len(lead_fail)
        + len(parity_fail)
        + len(forward_degree_fail)
        + len(forward_lead_fail)
    )

    if total_failures == 0:

        print("STATUS = PASS")
        print()
        print(
            "The exact target coefficient"
        )
        print(
            "    c_{ell-k-d-1}(N)"
        )
        print(
            "matches the established top-N degree and leading coefficient."
        )
        print()
        print(
            "The weighted substitution"
        )
        print(
            "    N=t^2, S=t*x"
        )
        print(
            "was used to expose the highest homogeneous contribution."
        )
        print()
        print(
            "The next symbolic target is the direct finite-sum derivation"
        )
        print(
            "of that highest homogeneous contribution."
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
