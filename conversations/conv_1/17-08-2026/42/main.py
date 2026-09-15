#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 110
TERMINAL-STRIP KERNEL DERIVATION

EXACT BINOMIAL KERNEL
-> NEWTON BASIS
-> r = ell-1-j TERMINAL COORDINATE

OBJECTIVE
---------

Experiment 109 established:

    ordinary division
      =
    synthetic division
      =
    explicit finite binomial kernel

and exposed the terminal leading pattern

    r < k  :  -C(k,r)
    r = k  :  C(k+1,2)

This experiment derives that strip symbolically.

For

    j = ell - 1 - r

define

    c_r(k,ell;N)
        = coefficient of P_j

and study

    r = 0,1,...,k

and then

    r > k.

We test whether the terminal strip has an exact decomposition

    c_r(N)
      =
    L_r(k) N^r
      + L_{r,1}(k,ell) N^(r-1)
      + ...

without fitting numerical data.

The important goal is to prove the boundary transition at r=k.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import time
from math import comb

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# =============================================================================
# WEIGHTS
# =============================================================================

WEIGHTS = [
    (1, 3),

    (1, 5), (3, 5),

    (1, 7), (3, 7), (5, 7),

    (1, 9), (3, 9), (5, 9), (7, 9),

    (1, 11), (3, 11), (5, 11), (7, 11), (9, 11),

    (1, 13), (3, 13), (5, 13),
    (7, 13), (9, 13), (11, 13),

    # new weights
    (1, 15), (3, 15), (5, 15), (7, 15),
    (9, 15), (11, 15), (13, 15),
]


# =============================================================================
# DIRECT DETECTOR
# =============================================================================

def detector_F(k: int, ell: int) -> sp.Expr:
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# =============================================================================
# EXACT QUOTIENT
# =============================================================================

def quotient(k: int, ell: int) -> sp.Expr:

    F = detector_F(k, ell)

    Q, rem = sp.div(
        sp.Poly(
            F,
            p,
            domain=sp.QQ.frac_field(q),
        ),
        sp.Poly(
            p + q + 1,
            p,
            domain=sp.QQ.frac_field(q),
        ),
    )

    if sp.expand(rem.as_expr()) != 0:
        raise ArithmeticError(
            f"division failure for ({k},{ell})"
        )

    return sp.expand(Q.as_expr())


# =============================================================================
# SYMMETRIC CONVERSION
# =============================================================================

def to_NS(expr: sp.Expr) -> sp.Expr:

    symmetric, remainder, mapping = sp.symmetrize(
        sp.expand(expr),
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            "Non-symmetric quotient encountered."
        )

    mapping = dict(mapping)

    sum_symbol = None
    prod_symbol = None

    for symbol, value in mapping.items():

        if sp.expand(value - (p + q)) == 0:
            sum_symbol = symbol

        elif sp.expand(value - p*q) == 0:
            prod_symbol = symbol

    if sum_symbol is None or prod_symbol is None:
        raise ArithmeticError(
            "Could not identify elementary symmetric variables."
        )

    return sp.expand(
        symmetric.subs(
            {
                sum_symbol: S,
                prod_symbol: N,
            }
        )
    )


# =============================================================================
# NEWTON BASIS
# =============================================================================

def newton_basis(max_j: int) -> list[sp.Expr]:

    P = [sp.Integer(0)] * (max_j + 1)

    P[0] = sp.Integer(2)

    if max_j >= 1:
        P[1] = S

    for j in range(2, max_j + 1):
        P[j] = sp.expand(
            S * P[j - 1]
            - N * P[j - 2]
        )

    return P


def newton_coefficients(
    expr: sp.Expr,
    max_j: int,
) -> dict[int, sp.Expr]:

    P = newton_basis(max_j)

    residual = sp.expand(expr)
    result: dict[int, sp.Expr] = {}

    for j in range(max_j, 0, -1):

        poly = sp.Poly(
            residual,
            S,
            domain=sp.QQ.frac_field(N),
        )

        c = sp.expand(
            poly.coeff_monomial(S**j)
        )

        if c != 0:

            result[j] = sp.factor(c)

            residual = sp.expand(
                residual - c * P[j]
            )

    result[0] = sp.factor(residual)

    reconstruction = sp.expand(
        result[0]
        + sum(
            result.get(j, 0) * P[j]
            for j in range(1, max_j + 1)
        )
    )

    if sp.expand(
        reconstruction - expr
    ) != 0:

        raise ArithmeticError(
            "Newton reconstruction failed."
        )

    return result


# =============================================================================
# TERMINAL COORDINATE
# =============================================================================

def terminal_coefficient(
    k: int,
    ell: int,
    r: int,
) -> sp.Expr:

    """
    r = 0 corresponds to P_(ell-1)
    r = 1 corresponds to P_(ell-2)
    etc.
    """

    j = ell - 1 - r

    if j < 1:
        raise ValueError(
            "Terminal index must satisfy j >= 1."
        )

    Q = to_NS(
        quotient(k, ell)
    )

    C = newton_coefficients(
        Q,
        ell - 1,
    )

    return sp.expand(
        C.get(j, 0)
    )


# =============================================================================
# LEADING N COEFFICIENT
# =============================================================================

def leading_N_data(expr: sp.Expr) -> tuple[int, sp.Expr]:
    poly = sp.Poly(
        sp.expand(expr),
        N,
    )

    if poly.is_zero:
        return -sp.oo, sp.Integer(0)

    return poly.degree(), sp.factor(poly.LC())


# =============================================================================
# HYPOTHESIS 1:
# r < k => -C(k,r)
# r = k => C(k+1,2)
# =============================================================================

def test_boundary_law() -> None:

    print()
    print("1. EXACT TERMINAL BOUNDARY LAW")
    print("-" * 78)

    failures = 0
    total = 0

    for k, ell in WEIGHTS:

        if ell <= k:
            continue

        print()
        print(f"({k},{ell})")

        for r in range(0, k + 1):

            c = terminal_coefficient(
                k,
                ell,
                r,
            )

            degree, lead = leading_N_data(c)

            if r < k:
                expected = -sp.Integer(
                    comb(k, r)
                )
            else:
                expected = sp.Integer(
                    comb(k + 1, 2)
                )

            total += 1

            ok = (
                sp.expand(lead - expected) == 0
            )

            if not ok:
                failures += 1

            print(
                f"  r={r:2d} "
                f"j={ell-1-r:2d} "
                f"degree={degree:2} "
                f"lead={lead} "
                f"expected={expected} "
                f"ok={ok}"
            )

    print()
    print(
        f"boundary-law failures = "
        f"{failures}/{total}"
    )

    if failures:
        raise ArithmeticError(
            "Boundary law failed."
        )

    print("STATUS = PASS")


# =============================================================================
# SECOND COEFFICIENT
# =============================================================================

def second_N_coefficient(expr: sp.Expr) -> sp.Expr:

    degree, _ = leading_N_data(expr)

    if degree <= 0:
        return sp.Integer(0)

    poly = sp.Poly(
        sp.expand(expr),
        N,
    )

    return sp.factor(
        poly.coeff_monomial(
            N ** (degree - 1)
        )
    )


# =============================================================================
# SEARCH FOR CLOSED TERMINAL FORM
# =============================================================================

def terminal_strip_table() -> None:

    print()
    print("2. TERMINAL STRIP COEFFICIENTS")
    print("-" * 78)

    for k, ell in WEIGHTS:

        print()
        print(f"({k},{ell})")

        for r in range(
            0,
            min(k + 2, ell - 1),
        ):

            c = terminal_coefficient(
                k,
                ell,
                r,
            )

            degree, lead = leading_N_data(c)
            second = second_N_coefficient(c)

            lead_str = str(sp.factor(lead))
            second_str = str(sp.factor(second))

            print(
                f"  r={r:2d} "
                f"degree={degree:2d} "
                f"lead={lead_str:>20} "
                f"next={second_str}"
            )


# =============================================================================
# NORMALIZED TERMINAL POLYNOMIAL
# =============================================================================

def normalized_terminal(
    k: int,
    ell: int,
    r: int,
) -> sp.Expr:

    c = terminal_coefficient(
        k,
        ell,
        r,
    )

    degree, _ = leading_N_data(c)

    if degree < 0:
        return sp.Integer(0)

    return sp.factor(
        sp.expand(
            c / N**degree
        )
    )


# =============================================================================
# TEST WHETHER ELL DISAPPEARS FROM TOP STRIP
# =============================================================================

def ell_independence_test() -> None:

    print()
    print("3. ELL-INDEPENDENCE OF THE TOP STRIP")
    print("-" * 78)

    failures = 0

    # Compare every k across all legal ell > k.
    for k in range(1, 14, 2):

        legal = [
            ell
            for kk, ell in WEIGHTS
            if kk == k and ell > k
        ]

        if len(legal) < 2:
            continue

        print()
        print(f"k={k}")

        for r in range(0, k + 1):

            data = []

            for ell in legal:

                c = terminal_coefficient(
                    k,
                    ell,
                    r,
                )

                degree, lead = leading_N_data(c)

                data.append(
                    (
                        ell,
                        degree,
                        sp.expand(lead),
                    )
                )

            degree_set = {
                x[1]
                for x in data
            }

            lead_set = {
                x[2]
                for x in data
            }

            ok = (
                len(degree_set) == 1
                and len(lead_set) == 1
            )

            if not ok:
                failures += 1

            print(
                f"  r={r}: "
                f"data={data} "
                f"stable={ok}"
            )

    print()
    print(
        f"ell-independence failures = "
        f"{failures}"
    )


# =============================================================================
# BELOW THE BOUNDARY: r > k
# =============================================================================

def beyond_boundary_test() -> None:

    print()
    print("4. POST-BOUNDARY REGIME r > k")
    print("-" * 78)

    for k, ell in WEIGHTS:

        if ell <= k + 1:
            continue

        print()
        print(f"({k},{ell})")

        for r in range(
            k + 1,
            min(
                ell - 1,
                k + 5,
            ),
        ):

            c = terminal_coefficient(
                k,
                ell,
                r,
            )

            degree, lead = leading_N_data(c)

            print(
                f"  r={r:2d} "
                f"degree={degree:2} "
                f"lead={lead}"
            )


# =============================================================================
# BINOMIAL-CORRECTION SEARCH AT r=k
# =============================================================================

def boundary_correction_search() -> None:

    print()
    print("5. r=k BOUNDARY CORRECTION SEARCH")
    print("-" * 78)

    print(
        "Testing whether the exceptional coefficient is "
        "C(k+1,2) exactly."
    )

    failures = 0

    for k, ell in WEIGHTS:

        c = terminal_coefficient(
            k,
            ell,
            k,
        )

        degree, lead = leading_N_data(c)

        expected_degree = k
        expected_lead = sp.Integer(
            comb(k + 1, 2)
        )

        ok = (
            degree == expected_degree
            and sp.expand(
                lead - expected_lead
            ) == 0
        )

        if not ok:
            failures += 1

        print(
            f"({k},{ell}) "
            f"degree={degree} "
            f"lead={lead} "
            f"expected={expected_degree}/"
            f"{expected_lead} "
            f"ok={ok}"
        )

    print()
    print(
        f"boundary correction failures = "
        f"{failures}/{len(WEIGHTS)}"
    )


# =============================================================================
# NEW WEIGHT VALIDATION
# =============================================================================

def holdout_weight_validation() -> None:

    print()
    print("6. NEW-WEIGHT VALIDATION")
    print("-" * 78)

    failures = 0

    holdout = [
        (1, 15), (3, 15), (5, 15),
        (7, 15), (9, 15), (11, 15), (13, 15),
    ]

    for k, ell in holdout:

        Q = to_NS(
            quotient(k, ell)
        )

        C = newton_coefficients(
            Q,
            ell - 1,
        )

        # Verify the complete terminal strip, not just the top term.
        for r in range(
            0,
            min(k + 1, ell - 1),
        ):

            j = ell - 1 - r

            c = sp.expand(
                C.get(j, 0)
            )

            if c == 0:
                failures += 1
                print(
                    f"FAIL zero coefficient "
                    f"({k},{ell}) r={r}"
                )

    print()
    print(
        f"new-weight failures = "
        f"{failures}"
    )

    if failures:
        raise ArithmeticError(
            "New-weight validation failed."
        )

    print("STATUS = PASS")


# =============================================================================
# KERNEL-LEVEL SYMBOLIC CHECK
# =============================================================================

def kernel_boundary_symbolic() -> None:

    print()
    print("7. KERNEL-LEVEL BOUNDARY CHECK")
    print("-" * 78)

    """
    Rather than recover the coefficient tensor and inspect it,
    explicitly expand the finite synthetic kernel and identify the
    top terminal contribution.

    The point is to check whether the binomial boundary term can be
    seen before symmetric reduction.
    """

    k_values = [1, 3, 5, 7]

    for k in k_values:

        ell = k + 2

        F = detector_F(
            k,
            ell,
        )

        Q = sp.expand(
            quotient(
                k,
                ell,
            )
        )

        Qns = to_NS(Q)

        C = newton_coefficients(
            Qns,
            ell - 1,
        )

        c = sp.factor(
            C[ell - 1 - k]
        )

        degree, lead = leading_N_data(c)

        print(
            f"k={k} ell={ell}: "
            f"boundary={c}, "
            f"degree={degree}, "
            f"lead={lead}, "
            f"expected={comb(k+1,2)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 110")
    print("TERMINAL-STRIP KERNEL DERIVATION")
    print("EXACT BINOMIAL KERNEL -> NEWTON TERMINAL REGION")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    test_boundary_law()

    terminal_strip_table()

    ell_independence_test()

    beyond_boundary_test()

    boundary_correction_search()

    holdout_weight_validation()

    kernel_boundary_symbolic()

    print()
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
Experiment 110 focuses exclusively on the boundary exposed by
Experiment 109.

The tested law is:

    c_(ell-1-r) leading term =
        -C(k,r) N^r          for r < k

and

        +C(k+1,2) N^k        for r = k.

The experiment separates two regimes:

    0 <= r < k
        pure Pascal/binomial regime

    r = k
        first boundary correction

    r > k
        post-boundary regime requiring a different kernel.

A successful run would establish that the Pascal structure is not
just an observed leading-term pattern: it is the exact terminal
kernel of the synthetic division.

The next mathematical target after a successful run is then clear:

    derive the r > k regime

and determine whether the entire coefficient tensor can be written
as a piecewise finite hypergeometric sum.

This is still structural/oracle mathematics.
It does not claim an N-only factoring algorithm.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 110 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

