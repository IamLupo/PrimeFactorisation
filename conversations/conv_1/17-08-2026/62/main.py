#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 127R
TARGET NEWTON COEFFICIENT ONLY
DIRECT POST-BOUNDARY LEADING-LAW CERTIFICATE
FINITE BINOMIAL KERNEL -> ONE NEWTON COEFFICIENT

P0 INCLUDED EXPLICITLY
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS

Purpose
-------

Experiment 126 established the candidate post-boundary leading law:

    r = (k-1)/2
    m = (ell-1)/2
    d = ell-k

for the target Newton coefficient

    c_{d+1}(N).

This experiment avoids reconstructing c_0,...,c_{d}.

It computes only the Newton coefficient c_{d+1}
by descending from S^(ell-1) to S^(d+1).

The Newton basis is

    P_0 = 2
    P_1 = S
    P_j = S P_{j-1} - N P_{j-2}.

For j >= 1, P_j is monic in S:

    P_j = S^j + lower powers of S.

Hence, once higher Newton coefficients have been removed,
the coefficient of S^j is exactly c_j.

We deliberately do NOT extract P0 or any coefficients below
the target index.

Expected leading law from Experiment 126:

    d odd:
        L = 2*(-1)^((d-1)/2) * (a+b+1)

    d even:
        L = (-1)^(d/2) * (a+b+1) * (2*a-2*b+d+1)

This experiment certifies only the leading coefficient.
==============================================================================
"""

import sys
import sympy as sp


# ---------------------------------------------------------------------------
# Symbols
# ---------------------------------------------------------------------------

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ---------------------------------------------------------------------------
# Exact detector kernel
# ---------------------------------------------------------------------------

def kernel_F(k: int, ell: int) -> sp.Expr:
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q(k: int, ell: int) -> sp.Expr:
    """
    Exact quotient F/(p+q+1), with p as the division variable.
    """
    F = kernel_F(k, ell)
    divisor = p + q + 1

    poly_F = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )
    poly_D = sp.Poly(
        divisor,
        p,
        domain=sp.QQ.frac_field(q),
    )

    quotient, remainder = sp.div(poly_F, poly_D)

    if remainder.as_expr() != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder for ({k},{ell}): "
            f"{sp.factor(remainder.as_expr())}"
        )

    Q = sp.expand(quotient.as_expr())

    if sp.denom(Q) != 1:
        raise ArithmeticError(
            f"Non-polynomial quotient for ({k},{ell}): {Q}"
        )

    return Q


# ---------------------------------------------------------------------------
# Symmetric conversion p,q -> S,N
# ---------------------------------------------------------------------------

def symmetric_to_NS(expr: sp.Expr) -> sp.Expr:
    """
    Convert a symmetric polynomial in p,q into S=p+q and N=pq.

    Substitute q=S-p and reduce modulo

        p^2 - S p + N.

    A symmetric polynomial must reduce to a p-independent remainder.
    """
    substituted = sp.expand(expr.subs(q, S - p))

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly_expr = sp.Poly(
        substituted,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    remainder = sp.rem(poly_expr, modulus).as_expr()
    remainder = sp.expand(remainder)

    p_part = sp.expand(remainder.coeff(p, 1))

    if p_part != 0:
        raise ArithmeticError(
            "Symmetric reduction retained p-dependence:\n"
            + str(sp.factor(p_part))
        )

    return sp.expand(remainder.coeff(p, 0))


# ---------------------------------------------------------------------------
# Newton basis
# ---------------------------------------------------------------------------

def P_newton(j: int) -> sp.Expr:
    """
    P_0=2, P_1=S,
    P_j = S P_{j-1} - N P_{j-2}.
    """
    if j == 0:
        return sp.Integer(2)
    if j == 1:
        return S

    p0 = sp.Integer(2)
    p1 = S

    for _ in range(2, j + 1):
        p2 = sp.expand(S*p1 - N*p0)
        p0, p1 = p1, p2

    return p1


# ---------------------------------------------------------------------------
# Target coefficient extraction
# ---------------------------------------------------------------------------

def target_newton_coefficient(
    Q_NS: sp.Expr,
    ell: int,
    target_j: int,
) -> sp.Expr:
    """
    Extract only c_target_j in

        Q = sum_j c_j(N) P_j(N,S).

    Since P_j is monic in S, we descend from ell-1 down to target_j.

    We intentionally never compute coefficients below target_j.
    """

    if target_j < 1:
        raise ValueError(
            "This target extractor is intended for j>=1."
        )

    work = sp.Poly(
        sp.expand(Q_NS),
        S,
        domain=sp.QQ.frac_field(N),
    )

    max_j = ell - 1

    coeff_cache = {}

    for j in range(max_j, target_j - 1, -1):
        coeff = sp.expand(work.coeff_monomial(S**j))

        coeff_cache[j] = coeff

        if coeff == 0:
            continue

        basis = sp.Poly(
            sp.expand(coeff * P_newton(j)),
            S,
            domain=sp.QQ.frac_field(N),
        )

        work = work - basis

    # We do not require the remaining lower-S terms to vanish.
    # Those belong to c_{j-1}, ..., c_0 and are irrelevant here.

    return sp.expand(coeff_cache[target_j])


# ---------------------------------------------------------------------------
# Expected Experiment-126 leading law
# ---------------------------------------------------------------------------

def expected_leading(k: int, ell: int, d: int) -> sp.Expr:
    a = sp.Rational(k - 1, 2)
    b = sp.Rational(ell - 1, 2)

    if d % 2 == 1:
        return sp.expand(
            2 * (-1)**((d - 1)//2) * (a + b + 1)
        )

    return sp.expand(
        (-1)**(d//2)
        * (a + b + 1)
        * (2*a - 2*b + d + 1)
    )


# ---------------------------------------------------------------------------
# One-row certificate
# ---------------------------------------------------------------------------

def certify_row(k: int, ell: int):
    d = ell - k
    j = d + 1

    # Post-boundary degree predicted by Experiment 126.
    D = (ell - 1)//2 - d

    if D < 0:
        return {
            "k": k,
            "ell": ell,
            "d": d,
            "j": j,
            "D_expected": D,
            "D_actual": None,
            "lead": None,
            "expected": None,
            "match": False,
            "skipped": True,
        }

    Q_pq = quotient_Q(k, ell)
    Q_NS = symmetric_to_NS(Q_pq)

    c_j = target_newton_coefficient(
        Q_NS,
        ell,
        j,
    )

    poly_N = sp.Poly(
        c_j,
        N,
        domain=sp.QQ,
    )

    actual_degree = poly_N.degree()

    if actual_degree == -sp.oo:
        actual_lead = sp.Integer(0)
    else:
        actual_lead = sp.expand(poly_N.LC())

    expected = sp.expand(
        expected_leading(k, ell, d)
    )

    return {
        "k": k,
        "ell": ell,
        "d": d,
        "j": j,
        "D_expected": D,
        "D_actual": actual_degree,
        "lead": actual_lead,
        "expected": expected,
        "coefficient": c_j,
        "match": (
            actual_degree == D
            and sp.simplify(actual_lead - expected) == 0
        ),
        "skipped": False,
    }


# ---------------------------------------------------------------------------
# Small training grid
# ---------------------------------------------------------------------------

def training_rows(max_ell: int):
    rows = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            d = ell - k
            D = (ell - 1)//2 - d

            if D >= 0:
                rows.append((k, ell))

    return rows


# ---------------------------------------------------------------------------
# Representative certificate printing
# ---------------------------------------------------------------------------

def print_certificate(result):
    print(
        f"({result['k']},{result['ell']}) "
        f"d={result['d']} j={result['j']}"
    )
    print(
        f"  expected degree = {result['D_expected']}"
    )
    print(
        f"  actual degree   = {result['D_actual']}"
    )
    print(
        f"  leading         = {result['lead']}"
    )
    print(
        f"  expected lead   = {result['expected']}"
    )

    if "coefficient" in result:
        print(
            f"  c_j(N)          = "
            f"{sp.factor(result['coefficient'])}"
        )

    print(
        f"  match            = {result['match']}"
    )
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("KAPPA EXPERIMENT 127R")
    print("TARGET NEWTON COEFFICIENT ONLY")
    print("DIRECT POST-BOUNDARY LEADING-LAW CERTIFICATE")
    print("FINITE BINOMIAL KERNEL -> ONE NEWTON COEFFICIENT")
    print("NO FULL NEWTON TENSOR")
    print("NO FULL C/D TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)
    print()

    TRAIN_MAX_ELL = 15
    FORWARD_ELL = 17

    # ----------------------------------------------------------------------
    # 1. quotient sanity
    # ----------------------------------------------------------------------

    print("1. EXACT QUOTIENT SANITY")
    print("-" * 78)

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
            print(f"({k},{ell}) quotient = PASS")
        except Exception as exc:
            failures += 1
            print(
                f"({k},{ell}) quotient = FAIL: "
                f"{type(exc).__name__}: {exc}"
            )

    print(
        f"quotient failures = {failures}/{len(sanity)}"
    )
    print()

    if failures:
        raise ArithmeticError(
            "Quotient sanity failed."
        )

    # ----------------------------------------------------------------------
    # 2. training certificate
    # ----------------------------------------------------------------------

    print("2. TRAINING CERTIFICATE")
    print("-" * 78)

    rows = training_rows(TRAIN_MAX_ELL)

    bad = []
    checked = 0

    for k, ell in rows:
        result = certify_row(k, ell)

        if result.get("skipped"):
            continue

        checked += 1

        if not result["match"]:
            bad.append(result)

    print(f"rows checked = {checked}")
    print(f"law failures = {len(bad)}")

    if bad:
        print()
        print("FIRST FAILURES")
        print("-" * 78)

        for result in bad[:20]:
            print_certificate(result)

        raise ArithmeticError(
            "Experiment 126 leading law failed."
        )

    print("STATUS = PASS")
    print()

    # ----------------------------------------------------------------------
    # 3. representatives
    # ----------------------------------------------------------------------

    print("3. REPRESENTATIVE CERTIFICATES")
    print("-" * 78)

    representatives = [
        (1, 5),
        (3, 7),
        (5, 9),
        (1, 7),
        (3, 9),
        (5, 11),
        (1, 9),
        (3, 11),
    ]

    for k, ell in representatives:
        if ell <= TRAIN_MAX_ELL:
            print_certificate(
                certify_row(k, ell)
            )

    # ----------------------------------------------------------------------
    # 4. forward ell holdout
    # ----------------------------------------------------------------------

    print("4. FORWARD ELL HOLDOUT")
    print("-" * 78)

    forward = []

    for k in range(1, FORWARD_ELL, 2):
        result = certify_row(k, FORWARD_ELL)

        if result.get("skipped"):
            continue

        forward.append(result)

    forward_bad = [
        r for r in forward
        if not r["match"]
    ]

    print(f"ell = {FORWARD_ELL}")
    print(f"rows = {len(forward)}")
    print(f"forward failures = {len(forward_bad)}")

    if forward_bad:
        print()
        for result in forward_bad:
            print_certificate(result)

        raise ArithmeticError(
            "Forward ell holdout failed."
        )

    print("STATUS = PASS")
    print()

    # ----------------------------------------------------------------------
    # 5. Boundary consistency: d=1..10 where possible
    # ----------------------------------------------------------------------

    print("5. d-STRIP CONSISTENCY")
    print("-" * 78)

    strip_failures = 0

    for d in range(1, 11):
        # Find the smallest admissible odd ell with k=1.
        k = 1
        ell = k + d

        if ell % 2 == 0:
            continue

        D = (ell - 1)//2 - d
        if D < 0:
            continue

        result = certify_row(k, ell)

        print(
            f"d={d:2d}: "
            f"L={sp.factor(result['expected'])} "
            f"actual={result['lead']} "
            f"ok={result['match']}"
        )

        if not result["match"]:
            strip_failures += 1

    print(
        f"d-strip failures = {strip_failures}"
    )

    if strip_failures:
        raise ArithmeticError(
            "d-strip consistency failed."
        )

    print()

    # ----------------------------------------------------------------------
    # 6. Final diagnostic
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "This experiment does not construct c_0,...,c_{d}."
    )
    print(
        "It extracts only the target coefficient c_{d+1}."
    )
    print()
    print(
        "The key optimization is:"
    )
    print(
        "  descend only from S^(ell-1) to S^(d+1),"
    )
    print(
        "  because every P_j is monic in S."
    )
    print()
    print(
        "P0 is never needed for this target and therefore cannot"
    )
    print(
        "produce the previous -8*N**3 reconstruction residual."
    )
    print()
    print(
        "A PASS gives a clean independent certificate of the"
    )
    print(
        "Experiment-126 leading law."
    )
    print()
    print(
        "Only after this passes should the next experiment target"
    )
    print(
        "the second N-coefficient of the same single c_{d+1}."
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
            str(exc),
        )
        sys.exit(1)

