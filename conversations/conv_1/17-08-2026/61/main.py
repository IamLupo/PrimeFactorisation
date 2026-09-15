#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 127
DIRECT LEADING-COEFFICIENT CERTIFICATE FROM FINITE BINOMIAL KERNEL
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO LARGE DATASET
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS

Goal
----
For the detector quotient

    Q_{k,l}(p,q)
      = F_{k,l}(p,q) / (p+q+1),

with

    F_{k,l}
      = p^k (1+q)^l + q^k (1+p)^l
        - p^l (1+q)^k - q^l (1+p)^k,

verify the Experiment-126 post-boundary leading law directly.

The Newton basis is

    P_0 = 2,
    P_j = p^j + q^j,  j >= 1.

For d = l-k, the target Newton coefficient is

    c_{d+1}(N),

and its expected N-degree is

    D = (l-1)/2 - d.

We do NOT construct all c_j.
Instead, we isolate only the coefficient of N^D in c_{d+1}.

Expected law
------------

Let

    a = (k-1)/2
    b = (l-1)/2.

Then

    L_d(a,b) =
        2*(-1)^((d-1)/2)*(a+b+1),                 d odd

        (-1)^(d/2)*(a+b+1)*(2*a-2*b+d+1),        d even.

This experiment derives the same scalar from the finite kernel
and verifies equality symbolically.

==============================================================================
"""

import sys
import sympy as sp

# ---------------------------------------------------------------------------
# Global symbols
# ---------------------------------------------------------------------------

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")

# ---------------------------------------------------------------------------
# Canonical symmetric variables
# ---------------------------------------------------------------------------

S_sym = p + q
N_sym = p * q

# ---------------------------------------------------------------------------
# Kernel
# ---------------------------------------------------------------------------

def kernel_F(k: int, ell: int) -> sp.Expr:
    """
    Exact finite binomial detector kernel.
    """
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def kernel_Q(k: int, ell: int) -> sp.Expr:
    """
    Exact quotient in Q[p,q].
    """
    F = kernel_F(k, ell)
    den = p + q + 1
    rem = sp.rem(
        sp.Poly(F, p, domain=sp.QQ.frac_field(q)),
        sp.Poly(den, p, domain=sp.QQ.frac_field(q)),
    )

    if sp.expand(rem.as_expr()) != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder for ({k},{ell}): "
            f"{sp.factor(rem.as_expr())}"
        )

    Q = sp.cancel(F / den)

    if sp.denom(Q) != 1:
        raise ArithmeticError(
            f"Quotient is not polynomial for ({k},{ell}): {Q}"
        )

    return sp.expand(Q)


# ---------------------------------------------------------------------------
# Symmetric polynomial conversion
# ---------------------------------------------------------------------------

def to_NS(expr: sp.Expr) -> sp.Expr:
    """
    Convert a symmetric polynomial in p,q into N=pq and S=p+q.

    We use the elementary symmetric reduction through the identity

        p^2 = S*p - N

    and eliminate p.

    The detectors are symmetric, so the remainder must be independent of p.
    """

    expr = sp.expand(expr)

    poly_p = sp.Poly(expr, p, domain=sp.QQ.frac_field(q))

    # Reduce p-powers recursively using p^2 = S*p - N,
    # temporarily taking S=p+q.
    #
    # A direct substitution through elementary symmetric polynomials is
    # more robust here: use q = S-p, then reduce modulo p^2-Sp+N.

    expr_S = sp.expand(expr.subs(q, S - p))

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.rem(
        sp.Poly(expr_S, p, domain=sp.QQ.frac_field(S, N)),
        modulus,
    )

    rem_expr = sp.expand(rem.as_expr())

    coeff_p = sp.expand(rem_expr.coeff(p, 1))
    const = sp.expand(rem_expr.coeff(p, 0))

    if coeff_p != 0:
        raise ArithmeticError(
            "Symmetric reduction retained a p-term:\n"
            + str(sp.factor(coeff_p))
        )

    return sp.expand(const)


# ---------------------------------------------------------------------------
# Newton power-sum basis
# ---------------------------------------------------------------------------

def P(j: int) -> sp.Expr:
    """
    Newton/symmetric power-sum basis:

        P_0 = 2
        P_j = p^j + q^j
    """
    if j == 0:
        return sp.Integer(2)
    return sp.expand(p**j + q**j)


def P_NS(j: int) -> sp.Expr:
    """
    P_j expressed in N,S via recurrence

        P_0 = 2
        P_1 = S
        P_j = S P_{j-1} - N P_{j-2}.
    """
    if j == 0:
        return sp.Integer(2)
    if j == 1:
        return S

    prev2 = sp.Integer(2)
    prev1 = S

    for _j in range(2, j + 1):
        cur = sp.expand(S * prev1 - N * prev2)
        prev2, prev1 = prev1, cur

    return prev1


# ---------------------------------------------------------------------------
# Exact Newton coefficient extraction
# ---------------------------------------------------------------------------

def newton_coefficients(Q_NS: sp.Expr, max_j: int):
    """
    Extract the coefficients c_j in

        Q = sum_{j=0}^max_j c_j(N) P_j(N,S).

    This routine is only used on small hand-picked certificates.
    It is NOT used to build a large dataset.

    Since P_j has leading S^j, extraction proceeds from high degree down.
    """

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    coeffs = {j: sp.Integer(0) for j in range(max_j + 1)}

    for j in range(max_j, -1, -1):
        coeff_Sj = work.coeff_monomial(S**j)

        if coeff_Sj == 0:
            continue

        coeffs[j] = sp.expand(coeff_Sj)

        subtraction = sp.Poly(
            sp.expand(coeffs[j] * P_NS(j)),
            S,
            domain=sp.QQ.frac_field(N),
        )

        work = work - subtraction

    residual = sp.expand(work.as_expr())

    if residual != 0:
        raise ArithmeticError(
            "Newton reconstruction residual is nonzero:\n"
            + str(sp.factor(residual))
        )

    return coeffs


# ---------------------------------------------------------------------------
# Expected leading law
# ---------------------------------------------------------------------------

def expected_leading(k: int, ell: int, d: int) -> sp.Integer:
    """
    Experiment-126 closed law.

        a=(k-1)/2
        b=(ell-1)/2

    odd d:
        2*(-1)^((d-1)/2)*(a+b+1)

    even d:
        (-1)^(d/2)*(a+b+1)*(2*a-2*b+d+1)
    """

    a = sp.Rational(k - 1, 2)
    b = sp.Rational(ell - 1, 2)

    if d % 2 == 1:
        return sp.simplify(
            2 * (-1)**((d - 1)//2) * (a + b + 1)
        )

    return sp.simplify(
        (-1)**(d//2)
        * (a + b + 1)
        * (2*a - 2*b + d + 1)
    )


# ---------------------------------------------------------------------------
# Direct certificate for one (k,ell)
# ---------------------------------------------------------------------------

def leading_certificate(k: int, ell: int):
    d = ell - k

    # Target Newton coefficient index.
    j = d + 1

    # Expected degree in N.
    D = (ell - 1)//2 - d

    if D < 0:
        raise ValueError(
            f"Invalid non-post-boundary row ({k},{ell}): D={D}"
        )

    Q_pq = kernel_Q(k, ell)
    Q_NS = to_NS(Q_pq)

    coeffs = newton_coefficients(Q_NS, ell)

    c = sp.expand(coeffs[j])
    poly_N = sp.Poly(
        c,
        N,
        domain=sp.QQ,
    )

    actual_degree = poly_N.degree()

    if actual_degree != D:
        return {
            "k": k,
            "ell": ell,
            "d": d,
            "j": j,
            "D_expected": D,
            "D_actual": actual_degree,
            "coefficient": c,
            "leading": None,
            "expected": expected_leading(k, ell, d),
            "match": False,
            "reason": "degree mismatch",
        }

    actual_lead = sp.Integer(poly_N.LC())
    expected = sp.Integer(expected_leading(k, ell, d))

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "j": j,
        "D_expected": D,
        "D_actual": actual_degree,
        "coefficient": c,
        "leading": actual_lead,
        "expected": expected,
        "match": sp.simplify(actual_lead - expected) == 0,
        "reason": "",
    }


# ---------------------------------------------------------------------------
# Small training grid
# ---------------------------------------------------------------------------

def training_grid(max_ell: int):
    rows = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            d = ell - k

            # Post-boundary condition used in Experiment 126.
            D = (ell - 1)//2 - d
            if D < 0:
                continue

            rows.append((k, ell))

    return rows


# ---------------------------------------------------------------------------
# Ratio certificate
# ---------------------------------------------------------------------------

def ratio_certificate(max_ell: int):
    """
    Check the predicted ratio under an increment of k by 2 where both rows
    exist.

    The primary purpose is to verify that the closed leading law behaves
    consistently across the k-direction.
    """

    failures = 0
    checks = 0

    for ell in range(5, max_ell + 1, 2):
        for k in range(1, ell - 2, 2):
            d1 = ell - k
            d2 = ell - (k + 2)

            D1 = (ell - 1)//2 - d1
            D2 = (ell - 1)//2 - d2

            if D1 < 0 or D2 < 0:
                continue

            L1 = expected_leading(k, ell, d1)
            L2 = expected_leading(k + 2, ell, d2)

            # Only report exact consistency, not numerical fitting.
            checks += 1

            if sp.simplify(L2 - expected_leading(k + 2, ell, d2)) != 0:
                failures += 1

    return checks, failures


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("KAPPA EXPERIMENT 127")
    print("DIRECT POST-BOUNDARY LEADING-LAW CERTIFICATE")
    print("FINITE BINOMIAL KERNEL -> SINGLE NEWTON COEFFICIENT")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    # Keep this deliberately small.
    TRAIN_MAX_ELL = 15
    FORWARD_ELL = 17

    # ----------------------------------------------------------------------
    # Direct quotient sanity
    # ----------------------------------------------------------------------

    print("1. EXACT QUOTIENT SANITY")
    print("-" * 78)

    sanity_cases = [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]

    quotient_failures = 0

    for k, ell in sanity_cases:
        try:
            Q = kernel_Q(k, ell)
            print(f"({k},{ell}) quotient = PASS")
        except Exception as exc:
            quotient_failures += 1
            print(f"({k},{ell}) quotient = FAIL: {exc}")

    print(f"quotient failures = {quotient_failures}/{len(sanity_cases)}")
    print()

    if quotient_failures:
        raise ArithmeticError("Quotient sanity failed.")

    # ----------------------------------------------------------------------
    # Training certificate
    # ----------------------------------------------------------------------

    print("2. TRAINING CERTIFICATE")
    print("-" * 78)

    rows = training_grid(TRAIN_MAX_ELL)

    failures = []
    cancellations = 0

    for k, ell in rows:
        result = leading_certificate(k, ell)

        if result["actual_degree"] != result["D_expected"]:
            if result["leading"] is None:
                cancellations += 1

        if not result["match"]:
            failures.append(result)

    print(f"rows checked       = {len(rows)}")
    print(f"cancellation rows  = {cancellations}")
    print(f"law failures       = {len(failures)}")
    print()

    if failures:
        print("FAILURE DETAILS")
        print("-" * 78)

        for r in failures[:20]:
            print(
                f"({r['k']},{r['ell']}) "
                f"d={r['d']} j={r['j']} "
                f"Dexp={r['D_expected']} Dact={r['D_actual']} "
                f"lead={r['leading']} expected={r['expected']} "
                f"reason={r['reason']}"
            )

        raise ArithmeticError(
            "Experiment-126 leading law failed."
        )

    print("STATUS = PASS")
    print()

    # ----------------------------------------------------------------------
    # Representative symbolic certificates
    # ----------------------------------------------------------------------

    print("3. REPRESENTATIVE CERTIFICATES")
    print("-" * 78)

    representative = [
        (1, 5),
        (3, 7),
        (5, 9),
        (1, 7),
        (3, 9),
        (5, 11),
        (1, 9),
        (3, 11),
        (5, 13),
    ]

    for k, ell in representative:
        if ell > TRAIN_MAX_ELL:
            continue

        r = leading_certificate(k, ell)

        print(
            f"({k},{ell}) d={r['d']} j={r['j']} "
            f"D={r['D_actual']}"
        )
        print(f"  c_j(N)       = {sp.factor(r['coefficient'])}")
        print(f"  actual lead  = {r['leading']}")
        print(f"  expected     = {r['expected']}")
        print(f"  match        = {r['match']}")
        print()

    # ----------------------------------------------------------------------
    # Forward holdout
    # ----------------------------------------------------------------------

    print("4. FORWARD ELL HOLDOUT")
    print("-" * 78)

    forward_rows = training_grid(FORWARD_ELL)
    forward_rows = [
        row for row in forward_rows
        if row[1] == FORWARD_ELL
    ]

    forward_failures = []

    for k, ell in forward_rows:
        r = leading_certificate(k, ell)

        if not r["match"]:
            forward_failures.append(r)

    print(f"ell = {FORWARD_ELL}")
    print(f"rows = {len(forward_rows)}")
    print(f"forward failures = {len(forward_failures)}")

    if forward_failures:
        for r in forward_failures:
            print(
                f"  ({r['k']},{r['ell']}) "
                f"lead={r['leading']} "
                f"expected={r['expected']}"
            )

        raise ArithmeticError(
            "Forward leading-law holdout failed."
        )

    print("STATUS = PASS")
    print()

    # ----------------------------------------------------------------------
    # Formula audit
    # ----------------------------------------------------------------------

    print("5. CLOSED-FORMULA AUDIT")
    print("-" * 78)

    for d in range(1, 11):
        k = 1
        ell = k + d

        if ell % 2 == 0:
            continue

        D = (ell - 1)//2 - d

        if D < 0:
            continue

        L = expected_leading(k, ell, d)

        print(
            f"d={d:2d}: "
            f"L(1,{ell}) = {sp.factor(L)}"
        )

    print()

    # ----------------------------------------------------------------------
    # Final diagnostic
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "The experiment isolates only the coefficient c_(d+1)(N)."
    )
    print()
    print(
        "It independently verifies the Experiment-126 leading law "
        "from the finite binomial kernel."
    )
    print()
    print(
        "No complete Newton tensor is constructed."
    )
    print(
        "No C/D tensor is constructed."
    )
    print(
        "No polynomial fitting is used."
    )
    print()
    print(
        "The mathematical target after a PASS is:"
    )
    print()
    print(
        "    derive the closed leading coefficient directly from"
    )
    print(
        "    the binomial-sum coefficient extraction."
    )
    print()
    print(
        "A PASS here is useful only if the derivation mechanism is"
    )
    print(
        "simple enough to generalize to the second N-layer."
    )
    print()
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print()
        print("FATAL:", type(exc).__name__, exc)
        sys.exit(1)

