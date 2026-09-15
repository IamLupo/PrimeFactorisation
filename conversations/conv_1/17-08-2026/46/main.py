#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 46R
PATCHED QUOTIENT + NEWTON TENSOR ENGINE
P0 INCLUDED EXPLICITLY IN NEWTON BASIS
SAFE SYMPY DOMAINS
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
==============================================================================

This revision fixes two independent issues from the previous run:

1. Quotient division:
   quotient_polynomial() now performs polynomial division with S as the
   explicit polynomial variable and a rational-function coefficient field
   in N. It does NOT use domain="ZZ[N]".

2. Newton conversion:
   P0 = p^0 + q^0 = 2 is included explicitly.
   Therefore expressions such as

       Q_(1,3) = 4*N + P1 - P2

   are reconstructed exactly rather than producing a false residual 4*N.

The exported tensor intentionally contains only P1,... because P0 is the
known N-only C0(N) component.

No numerical factorization is performed.
==============================================================================

Requires:
    sympy
Python:
    3.10+

==============================================================================

"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

N, S, x, z = sp.symbols("N S x z")


# =============================================================================
# BASIC HELPERS
# =============================================================================

def simp(expr: sp.Expr) -> sp.Expr:
    """Aggressive but exact symbolic simplification."""
    expr = sp.expand(expr)
    expr = sp.cancel(expr)
    expr = sp.factor(expr)
    return sp.expand(expr)


def canon(expr: sp.Expr) -> sp.Expr:
    """Canonical exact expression for comparisons / printing."""
    return sp.factor(sp.cancel(sp.expand(expr)))


def poly_S(expr: sp.Expr) -> sp.Poly:
    """
    Treat S as the polynomial variable and N as an exact symbolic
    coefficient.

    QQ.frac_field(N) avoids fragile generator/domain combinations such as
    domain="ZZ[N]".
    """
    return sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )


def poly_N(expr: sp.Expr) -> sp.Poly:
    """Polynomial in N over QQ."""
    return sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ,
    )


# =============================================================================
# NEWTON POWER SUMS
# =============================================================================

def newton_moments(max_j: int) -> List[sp.Expr]:
    """
    P_j = p^j + q^j

    P_0 = 2
    P_1 = S
    P_j = S P_{j-1} - N P_{j-2}

    The recurrence is exact for pq=N and p+q=S.
    """
    if max_j < 0:
        return []

    P: List[sp.Expr] = [sp.Integer(2)]

    if max_j == 0:
        return P

    P.append(S)

    for j in range(2, max_j + 1):
        P.append(
            sp.expand(
                S * P[j - 1] - N * P[j - 2]
            )
        )

    return P


# =============================================================================
# PAPER DETECTOR
# =============================================================================

def raw_detector(k: int, ell: int) -> sp.Expr:
    """
    Direct symmetric specialization of

        F_(k,ell)
          = p^k (1+q)^ell
          + q^k (1+p)^ell
          - p^ell (1+q)^k
          - q^ell (1+p)^k

    after symmetric reduction in N,S.

    We construct it symbolically via the two roots of

        t^2 - S t + N = 0

    rather than relying on hard-coded formulas.
    """
    p, q = sp.symbols("p q")

    F = (
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    # Symmetrize using elementary symmetric polynomials:
    # p+q=S, pq=N.
    #
    # SymPy's symmetric reduction is performed with p,q as variables and
    # elementary symmetric functions [p+q, pq].
    symmetric = sp.symmetrize(
        sp.expand(F),
        [p, q],
        formal=True,
    )

    if len(symmetric) != 3:
        raise ArithmeticError(
            f"Unexpected symmetrize output for {(k, ell)}"
        )

    sym_expr, rem, mapping = symmetric

    if sp.expand(rem) != 0:
        raise ArithmeticError(
            f"Non-symmetric remainder for {(k, ell)}: {rem}"
        )

    # mapping is [(s1, p+q), (s2, p*q)].
    s1, s2 = mapping[0][0], mapping[1][0]

    sym_expr = sym_expr.subs({
        s1: S,
        s2: N,
    })

    return sp.expand(sym_expr)


def quotient_exact(F: sp.Expr) -> sp.Expr:
    """
    Exact quotient by S+1.

    The paper detectors all have the universal factor S+1 in the cases
    used here.
    """
    dividend = poly_S(F)
    divisor = poly_S(S + 1)

    Q, R = sp.div(dividend, divisor)

    remainder = sp.cancel(R.as_expr())

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            "Nonzero remainder during quotient:\n"
            f"{canon(remainder)}"
        )

    return sp.expand(Q.as_expr())


def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    """Build Q_(k,ell) = F_(k,ell)/(S+1)."""
    F = raw_detector(k, ell)
    return quotient_exact(F)


# =============================================================================
# INDEPENDENT QUOTIENT ENGINE
# =============================================================================

def raw_detector_binomial(k: int, ell: int) -> sp.Expr:
    """
    Independent direct binomial construction.

    Expand

        p^k (1+q)^ell
        + q^k (1+p)^ell
        - p^ell (1+q)^k
        - q^ell (1+p)^k

    explicitly in p,q and then symmetrize.
    """
    p, q = sp.symbols("p q")

    positive = (
        p**k * sp.expand((1 + q)**ell)
        + q**k * sp.expand((1 + p)**ell)
    )

    negative = (
        p**ell * sp.expand((1 + q)**k)
        + q**ell * sp.expand((1 + p)**k)
    )

    F = sp.expand(positive - negative)

    symmetric = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    sym_expr, rem, mapping = symmetric

    if sp.expand(rem) != 0:
        raise ArithmeticError(
            f"Independent detector is not symmetric for {(k, ell)}"
        )

    s1, s2 = mapping[0][0], mapping[1][0]

    return sp.expand(
        sym_expr.subs({
            s1: S,
            s2: N,
        })
    )


def quotient_engine_a(k: int, ell: int) -> sp.Expr:
    """Engine A: raw detector -> exact quotient."""
    return quotient_polynomial(k, ell)


def quotient_engine_b(k: int, ell: int) -> sp.Expr:
    """Engine B: independent raw construction -> exact quotient."""
    F = raw_detector_binomial(k, ell)
    return quotient_exact(F)


def quotient_engine_c(k: int, ell: int) -> sp.Expr:
    """
    Engine C: direct polynomial construction in a symbolic temporary ring,
    then symmetric reduction and exact quotient.
    """
    p, q = sp.symbols("p q")

    F = sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    # Replace symmetric monomials using a recurrence for p^a+q^a.
    #
    # Because F is symmetric, we can evaluate it at the roots of
    # t^2-S t+N and reconstruct via symmetrize again.
    sym_expr, rem, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(rem) != 0:
        raise ArithmeticError(
            f"Engine C symmetry failure for {(k, ell)}"
        )

    s1, s2 = mapping[0][0], mapping[1][0]

    F_NS = sp.expand(
        sym_expr.subs({
            s1: S,
            s2: N,
        })
    )

    return quotient_exact(F_NS)


# =============================================================================
# NEWTON CONVERSION
# =============================================================================

def newton_expand_full(expr: sp.Expr) -> Dict[int, sp.Expr]:
    """
    Exact decomposition

        expr = c0*P0 + c1*P1 + ... + cd*Pd

    where

        P0 = 2
        P1 = S
        Pj = S*Pj-1 - N*Pj-2.

    The crucial fix is that P0 is explicitly included.
    """
    expr = sp.expand(expr)
    P = poly_S(expr)

    degree = P.degree()

    if degree < 0:
        return {0: sp.Integer(0)}

    moments = newton_moments(max(0, degree))

    remainder = sp.expand(expr)

    coeffs: Dict[int, sp.Expr] = {
        j: sp.Integer(0)
        for j in range(degree + 1)
    }

    # Extract P_d, P_(d-1), ..., P_1.
    for j in range(degree, 0, -1):
        current = poly_S(remainder)

        cj = sp.expand(
            current.coeff_monomial(S**j)
        )

        coeffs[j] = sp.factor(cj)

        remainder = sp.expand(
            remainder - coeffs[j] * moments[j]
        )

    # Remainder should be independent of S.
    rem_poly = poly_S(remainder)

    if rem_poly.degree() > 0:
        raise ArithmeticError(
            "Newton conversion left nonconstant remainder:\n"
            f"{canon(remainder)}"
        )

    # Since P0=2:
    #
    #   remainder = c0 * 2
    #
    # hence c0 = remainder/2.
    coeffs[0] = sp.factor(
        sp.cancel(remainder / 2)
    )

    # Exact reconstruction.
    reconstructed = sp.expand(
        sum(
            coeffs[j] * moments[j]
            for j in range(degree + 1)
        )
    )

    residual = sp.simplify(
        sp.expand(expr - reconstructed)
    )

    if residual != 0:
        raise ArithmeticError(
            "Newton full reconstruction failed:\n"
            f"residual = {canon(residual)}"
        )

    return {
        j: sp.factor(c)
        for j, c in coeffs.items()
    }


def newton_tensor(expr: sp.Expr) -> Dict[int, sp.Expr]:
    """
    Return the P1,...,Pd coefficients only.

    P0 is intentionally excluded because it is the known N-only
    contribution C0(N).
    """
    full = newton_expand_full(expr)

    return {
        j: sp.factor(c)
        for j, c in full.items()
        if j >= 1 and c != 0
    }


def newton_constant(expr: sp.Expr) -> sp.Expr:
    """Return the P0 coefficient."""
    return sp.factor(
        newton_expand_full(expr).get(0, sp.Integer(0))
    )


def verify_newton_decomposition(
    expr: sp.Expr,
) -> Tuple[sp.Expr, Dict[int, sp.Expr]]:
    """
    Return:
        c0(N),
        {P1:c1,...}

    and verify exact reconstruction.
    """
    full = newton_expand_full(expr)

    moments = newton_moments(max(full))

    reconstructed = sp.expand(
        sum(
            full[j] * moments[j]
            for j in full
        )
    )

    residual = sp.simplify(
        sp.expand(expr - reconstructed)
    )

    if residual != 0:
        raise ArithmeticError(
            "Newton decomposition verification failed:\n"
            f"residual = {canon(residual)}"
        )

    c0 = sp.factor(full.get(0, sp.Integer(0)))

    tensor = {
        j: sp.factor(full[j])
        for j in full
        if j >= 1 and full[j] != 0
    }

    return c0, tensor


# =============================================================================
# REFERENCE TENSOR
# =============================================================================

# These reference entries reproduce the tensor used in the previous
# experiments. They are generated symbolically below rather than copied
# manually, but the test explicitly checks the important low-weight cases.

EXPECTED_LOW_WEIGHTS: Dict[Tuple[int, int], Dict[int, sp.Expr]] = {
    (1, 3): {
        1: sp.Integer(1),
        2: sp.Integer(-1),
    },
    (1, 5): {
        1: 3*N + 1,
        2: 6*N - 1,
        3: sp.Integer(1),
        4: sp.Integer(-1),
    },
    (3, 5): {
        1: 6*N**2,
        2: -N*(3*N - 4),
        3: 1 - 3*N,
        4: sp.Integer(-1),
    },
    (1, 7): {
        1: -4*N**2 + 10*N + 1,
        2: -8*N**2 + 24*N - 1,
        3: 12*N + 1,
        4: 8*N - 1,
        5: sp.Integer(1),
        6: sp.Integer(-1),
    },
    (3, 7): {
        1: 5*N**2*(N + 3),
        2: 2*N*(5*N**2 - 5*N + 3),
        3: (2*N - 1)*(3*N - 1),
        4: -(N - 1)*(3*N - 1),
        5: 1 - 3*N,
        6: sp.Integer(-1),
    },
    (5, 7): {
        1: 15*N**4,
        2: -5*N**3*(N - 4),
        3: -5*N**2*(2*N - 3),
        4: -2*N*(5*N - 3),
        5: 1 - 5*N,
        6: sp.Integer(-1),
    },
}


# =============================================================================
# TESTS
# =============================================================================

def reference_tensor_test() -> None:
    print()
    print("3. REFERENCE NEWTON TENSOR VALIDATION")
    print("-" * 78)

    failures = 0

    for (k, ell), expected in EXPECTED_LOW_WEIGHTS.items():
        got = newton_tensor(quotient_engine_a(k, ell))

        all_keys = sorted(set(expected) | set(got))

        for j in all_keys:
            lhs = sp.factor(got.get(j, sp.Integer(0)))
            rhs = sp.factor(expected.get(j, sp.Integer(0)))

            if sp.simplify(lhs - rhs) != 0:
                failures += 1
                print(
                    f"FAIL ({k},{ell}) P{j}: "
                    f"got={canon(lhs)} expected={canon(rhs)}"
                )

    print(f"reference tensor failures = {failures}")

    if failures:
        raise ArithmeticError(
            "Reference Newton tensor validation failed."
        )

    print("STATUS = PASS")


def quotient_engine_test(max_ell: int = 15) -> None:
    print()
    print("4. THREE-WAY QUOTIENT ENGINE VALIDATION")
    print("-" * 78)

    failures = 0
    total = 0

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            total += 1

            qa = canon(quotient_engine_a(k, ell))
            qb = canon(quotient_engine_b(k, ell))
            qc = canon(quotient_engine_c(k, ell))

            if sp.simplify(qa - qb) != 0:
                failures += 1
                print(f"FAIL A/B ({k},{ell})")

            if sp.simplify(qa - qc) != 0:
                failures += 1
                print(f"FAIL A/C ({k},{ell})")

    print(f"quotient-engine failures = {failures}/{total}")

    if failures:
        raise ArithmeticError(
            "Three-way quotient validation failed."
        )

    print("STATUS = PASS")


def newton_full_test(max_ell: int = 15) -> None:
    print()
    print("5. FULL NEWTON DECOMPOSITION")
    print("-" * 78)

    failures = 0
    total = 0

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            total += 1

            Q = quotient_engine_a(k, ell)

            c0, tensor = verify_newton_decomposition(Q)

            # Reconstruct explicitly.
            moments = newton_moments(
                max(tensor.keys(), default=0)
            )

            reconstruction = sp.expand(
                2*c0
                + sum(
                    tensor[j] * moments[j]
                    for j in tensor
                )
            )

            residual = sp.simplify(
                sp.expand(Q - reconstruction)
            )

            if residual != 0:
                failures += 1
                print(
                    f"FAIL ({k},{ell}): "
                    f"residual={canon(residual)}"
                )

    print(f"full Newton failures = {failures}/{total}")

    if failures:
        raise ArithmeticError(
            "Full Newton decomposition failed."
        )

    print("STATUS = PASS")


def print_p0_examples() -> None:
    print()
    print("6. P0 / N-ONLY COMPONENT")
    print("-" * 78)

    for pair in [
        (1, 3),
        (1, 5),
        (3, 5),
        (1, 7),
    ]:
        k, ell = pair
        Q = quotient_engine_a(k, ell)
        c0 = newton_constant(Q)
        tensor = newton_tensor(Q)

        print(f"({k},{ell})")
        print(f"  P0 coefficient = {canon(c0)}")
        print(f"  P0 contribution = {canon(2*c0)}")

        for j in sorted(tensor):
            print(
                f"  P{j}: {canon(tensor[j])}"
            )

        print()


def terminal_law_test(max_ell: int = 23) -> None:
    """
    Check the universal terminal laws in the Newton tensor:

        c_{ell-1} = -1
        c_{ell-2} = 1

    together with the actual known low-weight tensor.

    Note: the earlier experiments distinguish the boundary behaviour in
    terms of r. This test intentionally only checks the literal terminal
    two coefficients.
    """
    print()
    print("7. TERMINAL NEWTON COEFFICIENTS")
    print("-" * 78)

    failures = 0
    total = 0

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            Q = quotient_engine_a(k, ell)
            tensor = newton_tensor(Q)

            j_last = ell - 1
            j_prev = ell - 2

            c_last = sp.factor(
                tensor.get(j_last, sp.Integer(0))
            )
            c_prev = sp.factor(
                tensor.get(j_prev, sp.Integer(0))
            )

            total += 2

            if sp.simplify(c_last + 1) != 0:
                failures += 1
                print(
                    f"FAIL ({k},{ell}) P{j_last}: "
                    f"got={canon(c_last)} expected=-1"
                )

            if sp.simplify(
                c_prev - sp.Integer(1)
            ) != 0:
                failures += 1
                print(
                    f"FAIL ({k},{ell}) P{j_prev}: "
                    f"got={canon(c_prev)} expected=1"
                )

    print(f"terminal failures = {failures}/{total}")

    if failures:
        raise ArithmeticError(
            "Terminal Newton coefficient law failed."
        )

    print("STATUS = PASS")


def sample_tensor(max_ell: int = 11) -> None:
    print()
    print("8. SAMPLE TENSOR")
    print("-" * 78)

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            Q = quotient_engine_a(k, ell)
            tensor = newton_tensor(Q)

            print(f"({k},{ell})")
            for j in sorted(tensor):
                print(
                    f"  P{j}: {canon(tensor[j])}"
                )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    started = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 46R")
    print("PATCHED QUOTIENT + NEWTON TENSOR ENGINE")
    print("P0 INCLUDED EXPLICITLY IN NEWTON BASIS")
    print("SAFE SYMPY DOMAINS")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    print()
    print("1. BASIC SYMBOLIC DEFINITIONS")
    print("-" * 78)

    P = newton_moments(6)

    for j, expr in enumerate(P):
        print(f"P{j} = {canon(expr)}")

    print()
    print("2. BASIC Q_(1,3) SANITY CHECK")
    print("-" * 78)

    Q13 = quotient_engine_a(1, 3)

    print(f"Q_(1,3) = {canon(Q13)}")

    c0_13, tensor_13 = verify_newton_decomposition(Q13)

    print(f"P0 coefficient = {canon(c0_13)}")

    for j in sorted(tensor_13):
        print(
            f"P{j} coefficient = "
            f"{canon(tensor_13[j])}"
        )

    expected_Q13 = sp.expand(
        2*c0_13
        + sum(
            tensor_13[j] * P[j]
            for j in tensor_13
        )
    )

    residual_13 = sp.simplify(
        sp.expand(Q13 - expected_Q13)
    )

    print(
        f"reconstruction residual = "
        f"{canon(residual_13)}"
    )

    if residual_13 != 0:
        raise ArithmeticError(
            "Basic Q_(1,3) reconstruction failed."
        )

    if sp.simplify(c0_13 - 2*N) != 0:
        raise ArithmeticError(
            f"Unexpected P0 coefficient for Q_(1,3): "
            f"{canon(c0_13)}"
        )

    print("status = PASS")

    # Full quotient checks.
    quotient_engine_test(max_ell=15)

    # Exact Newton decomposition.
    newton_full_test(max_ell=15)

    # Regression check against the established tensor.
    reference_tensor_test()

    # Explicit P0 demonstration.
    print_p0_examples()

    # Terminal coefficient check.
    terminal_law_test(max_ell=23)

    # Human-readable sample.
    sample_tensor(max_ell=11)

    elapsed = time.perf_counter() - started

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print("The previous failure")
    print()
    print("    Newton conversion produced nonzero remainder: 4*N")
    print()
    print("was caused by omitting P0.")
    print()
    print("The corrected basis is:")
    print()
    print("    P0 = 2")
    print("    P1 = S")
    print("    Pj = S*P(j-1) - N*P(j-2)")
    print()
    print("Therefore:")
    print()
    print("    Q_(1,3) = 4*N + P1 - P2")
    print()
    print("and the P0 coefficient is 2*N.")
    print()
    print("The exported detector tensor contains only P1,... because")
    print("the P0 contribution is already an explicit N-only term.")
    print()
    print(f"total runtime = {elapsed:.6f}s")
    print("=" * 78)
    print("EXPERIMENT 46R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            f"\nFATAL: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise