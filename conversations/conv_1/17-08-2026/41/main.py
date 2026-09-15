#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 109
EXACT SYNTHETIC-DIVISION BINOMIAL KERNEL

DIRECT DETECTOR
    ->
BINOMIAL EXPANSION
    ->
SYNTHETIC DIVISION BY (p+q+1)
    ->
EXACT FINITE-SUM QUOTIENT KERNEL
    ->
SYMMETRIC (N,S)
    ->
NEWTON MOMENT COEFFICIENTS

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN

Purpose
-------
Experiment 108 established:

  * direct detector identity = exact;
  * independent Newton coefficient engine = exact;
  * corrected terminal coefficients = exact;
  * ell=13 symbolic holdout = exact.

But the attempted degree law was wrong.

This experiment derives the quotient coefficients directly from the
binomial expansion using synthetic division in p:

    F(p,q) = (p + q + 1) Q(p,q).

If

    F(p,q) = sum_i f_i(q) p^i

and

    Q(p,q) = sum_i g_i(q) p^i,

then

    g_{d-1} = f_d
    g_{i-1} = f_i - (q+1) g_i.

Therefore

    g_i =
        sum_{t=0}^{d-1-i}
            (-1)^t (q+1)^t f_{i+1+t}.

This gives an exact finite binomial kernel.

The experiment compares:

  A. ordinary polynomial division;
  B. synthetic-division recurrence;
  C. explicit finite-sum formula;
  D. symmetric Newton decomposition.

It then searches for:

  * terminal coefficient laws;
  * binomial leading-coefficient laws;
  * factorization of the finite-sum coefficients;
  * stability under new weights.

The objective is to discover the actual kernel, not fit it.
==============================================================================
"""

from __future__ import annotations

from math import comb
import time

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# =============================================================================
# WEIGHTS
# =============================================================================

DISCOVERY = (
    (1, 3),
    (1, 5), (3, 5),
    (1, 7), (3, 7), (5, 7),
    (1, 9), (3, 9), (5, 9), (7, 9),
    (1, 11), (3, 11), (5, 11),
    (7, 11), (9, 11),
)

HOLDOUT = (
    (1, 13), (3, 13), (5, 13),
    (7, 13), (9, 13), (11, 13),
)

ALL = DISCOVERY + HOLDOUT


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
# ORDINARY DIVISION
# =============================================================================

def ordinary_quotient(k: int, ell: int) -> sp.Expr:

    F = detector_F(k, ell)

    Q, rem = sp.div(
        sp.Poly(F, p, domain=sp.QQ.frac_field(q)),
        sp.Poly(p + q + 1, p, domain=sp.QQ.frac_field(q)),
    )

    if sp.expand(rem.as_expr()) != 0:
        raise ArithmeticError(
            f"ordinary division remainder != 0 for ({k},{ell})"
        )

    return sp.expand(Q.as_expr())


# =============================================================================
# COEFFICIENT ARRAY f_i(q)
# =============================================================================

def p_coefficients(expr: sp.Expr) -> list[sp.Expr]:

    poly = sp.Poly(
        sp.expand(expr),
        p,
        domain=sp.QQ.frac_field(q),
    )

    degree = poly.degree()

    coeffs = [
        sp.expand(
            poly.coeff_monomial(p**i)
        )
        for i in range(degree + 1)
    ]

    return coeffs


# =============================================================================
# SYNTHETIC DIVISION
# =============================================================================

def synthetic_quotient(
    k: int,
    ell: int,
) -> sp.Expr:

    F = detector_F(k, ell)

    f = p_coefficients(F)

    d = len(f) - 1

    # Quotient degree is d-1.
    g = [sp.Integer(0)] * d

    g[d - 1] = f[d]

    for i in range(d - 1, 0, -1):
        g[i - 1] = sp.expand(
            f[i] - (q + 1) * g[i]
        )

    Q = sum(
        g[i] * p**i
        for i in range(d)
    )

    Q = sp.expand(Q)

    # Exact reconstruction.
    if sp.expand(
        (p + q + 1) * Q - F
    ) != 0:
        raise ArithmeticError(
            f"synthetic quotient reconstruction failed "
            f"for ({k},{ell})"
        )

    return Q


# =============================================================================
# EXPLICIT FINITE SUM
# =============================================================================

def synthetic_kernel_coeff(
    f: list[sp.Expr],
    i: int,
) -> sp.Expr:
    """
    Explicit finite kernel:

        g_i =
          sum_{t=0}^{d-1-i}
            (-1)^t (q+1)^t f_{i+1+t}.
    """

    d = len(f) - 1

    result = sp.Integer(0)

    for t in range(d - 1 - i + 1):

        result += (
            (-1)**t
            * (q + 1)**t
            * f[i + 1 + t]
        )

    return sp.expand(result)


def explicit_kernel_quotient(
    k: int,
    ell: int,
) -> sp.Expr:

    F = detector_F(k, ell)

    f = p_coefficients(F)

    d = len(f) - 1

    g = [
        synthetic_kernel_coeff(
            f,
            i,
        )
        for i in range(d)
    ]

    Q = sp.expand(
        sum(
            g[i] * p**i
            for i in range(d)
        )
    )

    if sp.expand(
        (p + q + 1) * Q - F
    ) != 0:
        raise ArithmeticError(
            f"explicit kernel reconstruction failed "
            f"for ({k},{ell})"
        )

    return Q


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
            "Expression is not symmetric."
        )

    mapping = dict(mapping)

    s_symbol = None
    n_symbol = None

    for symbol, value in mapping.items():

        if sp.expand(
            value - (p + q)
        ) == 0:
            s_symbol = symbol

        if sp.expand(
            value - p*q
        ) == 0:
            n_symbol = symbol

    if s_symbol is None or n_symbol is None:
        raise ArithmeticError(
            "Could not identify symmetric generators."
        )

    return sp.expand(
        symmetric.subs(
            {
                s_symbol: S,
                n_symbol: N,
            }
        )
    )


# =============================================================================
# NEWTON BASIS
# =============================================================================

def newton_basis(max_j: int) -> list[sp.Expr]:

    P = [
        sp.Integer(0)
        for _ in range(max_j + 1)
    ]

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
    Q: sp.Expr,
    max_j: int,
) -> dict[int, sp.Expr]:

    P = newton_basis(max_j)

    residual = sp.expand(Q)

    coefficients: dict[int, sp.Expr] = {}

    for j in range(max_j, 0, -1):

        poly = sp.Poly(
            residual,
            S,
            domain=sp.QQ.frac_field(N),
        )

        c = sp.factor(
            poly.coeff_monomial(S**j)
        )

        if c != 0:

            coefficients[j] = c

            residual = sp.expand(
                residual - c * P[j]
            )

    coefficients[0] = sp.factor(
        residual
    )

    reconstruction = sp.expand(
        coefficients[0]
        + sum(
            coefficients.get(j, 0) * P[j]
            for j in range(1, max_j + 1)
        )
    )

    if sp.expand(
        reconstruction - Q
    ) != 0:
        raise ArithmeticError(
            "Newton reconstruction failed."
        )

    return coefficients


# =============================================================================
# DIRECT BINOMIAL COEFFICIENTS f_i(q)
# =============================================================================

def raw_binomial_coefficients(
    k: int,
    ell: int,
) -> list[sp.Expr]:

    """
    Build the p-coefficient array directly from the binomial theorem.

    Contributions:

      p^k (1+q)^ell
      q^k (1+p)^ell
      -p^ell (1+q)^k
      -q^ell (1+p)^k
    """

    d = max(k, ell)

    f = [
        sp.Integer(0)
        for _ in range(d + 1)
    ]

    # p^k (1+q)^ell
    for r in range(ell + 1):

        power_p = k

        if power_p <= d:
            f[power_p] += (
                sp.Integer(
                    comb(ell, r)
                )
                * q**r
            )

    # q^k (1+p)^ell
    for r in range(ell + 1):

        f[r] += (
            q**k
            * comb(ell, r)
        )

    # -p^ell (1+q)^k
    for r in range(k + 1):

        f[ell] -= (
            comb(k, r)
            * q**r
        )

    # -q^ell (1+p)^k
    for r in range(k + 1):

        f[r] -= (
            q**ell
            * comb(k, r)
        )

    return [
        sp.expand(x)
        for x in f
    ]


# =============================================================================
# VALIDATE RAW BINOMIAL ARRAY
# =============================================================================

def validate_raw_array() -> None:

    print()
    print("1. BINOMIAL RAW-ARRAY VALIDATION")
    print("-" * 78)

    failures = 0

    for k, ell in ALL:

        F = detector_F(
            k,
            ell,
        )

        f1 = p_coefficients(F)
        f2 = raw_binomial_coefficients(
            k,
            ell,
        )

        if len(f1) != len(f2):
            failures += 1
            print(
                f"FAIL length ({k},{ell})"
            )
            continue

        for a, b in zip(f1, f2):

            if sp.expand(a - b) != 0:
                failures += 1
                print(
                    f"FAIL raw coefficient ({k},{ell})"
                )
                break

    print(
        f"raw coefficient failures = "
        f"{failures}/{len(ALL)}"
    )

    if failures:
        raise ArithmeticError(
            "Raw binomial coefficient construction failed."
        )


# =============================================================================
# COMPARE ALL THREE QUOTIENT ENGINES
# =============================================================================

def compare_quotient_engines() -> None:

    print()
    print("2. THREE-WAY QUOTIENT VALIDATION")
    print("-" * 78)

    failures = 0

    for k, ell in ALL:

        Q1 = ordinary_quotient(
            k,
            ell,
        )

        Q2 = synthetic_quotient(
            k,
            ell,
        )

        Q3 = explicit_kernel_quotient(
            k,
            ell,
        )

        if sp.expand(Q1 - Q2) != 0:
            failures += 1
            print(
                f"FAIL ordinary/synthetic ({k},{ell})"
            )

        if sp.expand(Q1 - Q3) != 0:
            failures += 1
            print(
                f"FAIL ordinary/kernel ({k},{ell})"
            )

    print(
        f"quotient-engine failures = "
        f"{failures}/{len(ALL) * 2}"
    )

    if failures:
        raise ArithmeticError(
            "Quotient engines disagree."
        )


# =============================================================================
# COEFFICIENT TENSOR FROM EXPLICIT KERNEL
# =============================================================================

def compare_newton_tensor() -> None:

    print()
    print("3. NEWTON TENSOR FROM EXPLICIT BINOMIAL KERNEL")
    print("-" * 78)

    failures = 0

    for k, ell in ALL:

        Q_kernel = to_NS(
            explicit_kernel_quotient(
                k,
                ell,
            )
        )

        Q_direct = to_NS(
            ordinary_quotient(
                k,
                ell,
            )
        )

        C_kernel = newton_coefficients(
            Q_kernel,
            ell - 1,
        )

        C_direct = newton_coefficients(
            Q_direct,
            ell - 1,
        )

        for j in range(0, ell):

            a = sp.factor(
                C_kernel.get(j, 0)
            )

            b = sp.factor(
                C_direct.get(j, 0)
            )

            if sp.expand(a - b) != 0:
                failures += 1
                print(
                    f"FAIL ({k},{ell}) j={j}"
                )

    print(
        f"Newton coefficient failures = "
        f"{failures}"
    )

    if failures:
        raise ArithmeticError(
            "Newton tensor reconstruction failed."
        )


# =============================================================================
# TERMINAL LAW
# =============================================================================

def terminal_law() -> None:

    print()
    print("4. TERMINAL LAWS")
    print("-" * 78)

    failures = 0
    total = 0

    for k, ell in ALL:

        Q = to_NS(
            explicit_kernel_quotient(
                k,
                ell,
            )
        )

        C = newton_coefficients(
            Q,
            ell - 1,
        )

        top = sp.expand(
            C.get(ell - 1, 0)
        )

        prev = sp.expand(
            C.get(ell - 2, 0)
        )

        expected_top = sp.Integer(-1)

        expected_prev = (
            sp.Integer(1)
            if k == 1
            else sp.expand(1 - k*N)
        )

        total += 2

        if sp.expand(
            top - expected_top
        ) != 0:

            failures += 1
            print(
                f"FAIL top ({k},{ell})"
            )

        if sp.expand(
            prev - expected_prev
        ) != 0:

            failures += 1
            print(
                f"FAIL prev ({k},{ell}): "
                f"{prev} != {expected_prev}"
            )

    print(
        f"terminal failures = "
        f"{failures}/{total}"
    )

    if failures:
        raise ArithmeticError(
            "Terminal law failed."
        )

    print("STATUS = PASS")


# =============================================================================
# LEADING DEGREE TABLE
# =============================================================================

def leading_table() -> None:

    print()
    print("5. LEADING N-DEGREE / COEFFICIENT TABLE")
    print("-" * 78)

    for ell in sorted(
        {ell for _, ell in ALL}
    ):

        print()
        print(
            f"ell={ell}"
        )

        for k, e in ALL:

            if e != ell:
                continue

            Q = to_NS(
                explicit_kernel_quotient(
                    k,
                    ell,
                )
            )

            C = newton_coefficients(
                Q,
                ell - 1,
            )

            row = []

            for j in range(1, ell):

                c = sp.expand(
                    C.get(j, 0)
                )

                poly = sp.Poly(
                    c,
                    N,
                )

                row.append(
                    (
                        j,
                        poly.degree(),
                        sp.factor(
                            poly.LC()
                        ),
                    )
                )

            print(
                f"  k={k}: {row}"
            )


# =============================================================================
# BINOMIAL LEADING-COEFFICIENT SEARCH
# =============================================================================

def binomial_leading_search() -> None:

    print()
    print("6. BINOMIAL LEADING-COEFFICIENT SEARCH")
    print("-" * 78)

    """
    The data in Experiment 108 showed rows such as

      9,36,84,126,126,84,36,9,1

    which are binomial coefficients.

    We test a small collection of exact hypotheses.

    For each coefficient, let:

        d = deg_N(c_j)

    and L = leading coefficient.

    Candidate families tested:

      A. +/- C(k, r)
      B. +/- k*C(k-1, r)
      C. +/- C(k, r) times a small polynomial in k
      D. +/- C(ell-j-1, r)

    We do NOT fit arbitrary polynomials.
    """

    successes = []

    for k, ell in ALL:

        Q = to_NS(
            explicit_kernel_quotient(
                k,
                ell,
            )
        )

        C = newton_coefficients(
            Q,
            ell - 1,
        )

        for j in range(1, ell):

            c = sp.expand(
                C.get(j, 0)
            )

            if c == 0:
                continue

            poly = sp.Poly(
                c,
                N,
            )

            d = poly.degree()
            L = sp.expand(
                poly.LC()
            )

            # Search small integer r values.
            candidates = []

            for r in range(k + 2):

                if r <= k:
                    candidates.extend(
                        [
                            ("C(k,r)", sp.Integer(comb(k, r))),
                            (
                                "-C(k,r)",
                                -sp.Integer(comb(k, r)),
                            ),
                        ]
                    )

            for r in range(ell + 2):

                if r <= ell:
                    candidates.extend(
                        [
                            (
                                "C(ell,r)",
                                sp.Integer(comb(ell, r)),
                            ),
                            (
                                "-C(ell,r)",
                                -sp.Integer(comb(ell, r)),
                            ),
                        ]
                    )

            found = [
                name
                for name, value in candidates
                if sp.expand(L - value) == 0
            ]

            if found:
                successes.append(
                    (
                        k,
                        ell,
                        j,
                        d,
                        found[:4],
                    )
                )

    print(
        f"leading coefficients matching raw binomials = "
        f"{len(successes)}"
    )

    for item in successes[:80]:
        print(
            f"  ({item[0]},{item[1]}) "
            f"j={item[2]} degree={item[3]} "
            f"matches={item[4]}"
        )


# =============================================================================
# TERMINAL BINOMIAL STRIP
# =============================================================================

def terminal_binomial_strip() -> None:

    print()
    print("7. TERMINAL BINOMIAL STRIP")
    print("-" * 78)

    """
    Inspect

        c_(ell-1-r)

    for small r.

    We test whether the top N-degree coefficient follows

        +/- C(k,r)

    for the first several r.

    This is exactly the visible Pascal-row phenomenon from Experiment 108.
    """

    for k, ell in ALL:

        Q = to_NS(
            explicit_kernel_quotient(
                k,
                ell,
            )
        )

        C = newton_coefficients(
            Q,
            ell - 1,
        )

        print()
        print(
            f"({k},{ell})"
        )

        max_r = min(
            6,
            ell - 2,
            k,
        )

        for r in range(
            0,
            max_r + 1,
        ):

            j = ell - 1 - r

            c = sp.expand(
                C.get(j, 0)
            )

            poly = sp.Poly(
                c,
                N,
            )

            degree = poly.degree()
            leading = sp.factor(
                poly.LC()
            )

            expected = -comb(
                k,
                r,
            )

            print(
                f"  r={r:2d} "
                f"j={j:2d} "
                f"degree={degree} "
                f"lead={leading} "
                f"expected_binomial={expected} "
                f"match={sp.expand(leading-expected)==0}"
            )


# =============================================================================
# HOLDOUT
# =============================================================================

def holdout_validation() -> None:

    print()
    print("8. NEW-WEIGHT HOLDOUT")
    print("-" * 78)

    failures = 0

    for k, ell in HOLDOUT:

        Q_kernel = to_NS(
            explicit_kernel_quotient(
                k,
                ell,
            )
        )

        Q_direct = to_NS(
            ordinary_quotient(
                k,
                ell,
            )
        )

        if sp.expand(
            Q_kernel - Q_direct
        ) != 0:

            failures += 1

            print(
                f"FAIL ({k},{ell})"
            )

    print(
        f"holdout failures = "
        f"{failures}/{len(HOLDOUT)}"
    )

    if failures:
        raise ArithmeticError(
            "Holdout symbolic validation failed."
        )

    print("STATUS = PASS")


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 109")
    print("EXACT SYNTHETIC-DIVISION BINOMIAL KERNEL")
    print("DIRECT DETECTOR -> FINITE BINOMIAL SUM -> NEWTON TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    validate_raw_array()

    compare_quotient_engines()

    compare_newton_tensor()

    terminal_law()

    leading_table()

    binomial_leading_search()

    terminal_binomial_strip()

    holdout_validation()

    print()
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
Experiment 109 changes the question.

We no longer ask whether the Newton coefficients can be fitted.

We construct them exactly from:

    F_(k,l)
      =
    p^k(1+q)^l + q^k(1+p)^l
    - p^l(1+q)^k - q^l(1+p)^k

followed by exact synthetic division by

    p+q+1.

This produces the finite kernel

    g_i =
      sum_t (-1)^t (q+1)^t f_(i+1+t)(q),

where every f_i(q) is itself a finite binomial sum.

The important outcomes are:

  A. THREE-WAY EXACT AGREEMENT
     ordinary division =
     synthetic recurrence =
     explicit finite binomial kernel.

  B. NEWTON TENSOR DERIVATION
     The Newton coefficient tensor follows from an explicit
     finite combinatorial construction.

  C. BINOMIAL LEADING TERMS
     Test whether the observed Pascal-row coefficients are
     exact consequences of the kernel.

  D. HOLDOUT WEIGHTS
     The construction must remain exact for ell=13.

If the terminal/leading coefficients reduce to binomial sums,
the next step is to simplify those sums analytically into a
closed hypergeometric or cyclotomic expression.

If that succeeds, we will have explained the detector tensor
from first principles.

This still does not solve the N-only evaluation problem.
It is a structural derivation experiment.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 109 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

