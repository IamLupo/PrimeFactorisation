#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 103R
PATCHED INVERSE-RESOLVENT ARITHMETIC COLLAPSE

FOLLOW-UP TO 102R4-P
------------------------------------------------------------------------------
Goal:

    exact rational P1 resolvent
        ->
    rewrite in sigma_r(N) coordinates
        ->
    test whether the result collapses arithmetically.

PATCH:
    SymPy symbolic matrix equality is checked entry-by-entry with simplify().

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
STRICT TARGET HOLDOUT
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
import random
import time

import sympy as sp


# =============================================================================
# PARAMETERS
# =============================================================================

P_MIN = 2_000_000
P_MAX = 4_200_000

TARGET_COUNT = 50
TRAIN_COUNT = 35

RANDOM_SEED = 103

DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

MAX_MOMENT = 6


# =============================================================================
# TARGET
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> list[int]:

    if hi < 2 or lo > hi:
        return []

    flags = bytearray(
        b"\x01"
    ) * (hi + 1)

    flags[0] = 0
    flags[1] = 0

    root = isqrt(hi)

    for p in range(
        2,
        root + 1,
    ):
        if flags[p]:

            start = p * p
            count = (
                (hi - start) // p
                + 1
            )

            flags[
                start:hi + 1:p
            ] = (
                b"\x00"
            ) * count

    return [
        n
        for n in range(
            max(2, lo),
            hi + 1,
        )
        if flags[n]
    ]


# =============================================================================
# TARGET GENERATION
# =============================================================================

def generate_targets(
    primes: list[int],
    count: int,
) -> list[Target]:

    rng = random.Random(
        RANDOM_SEED
    )

    out: list[Target] = []
    seen: set[int] = set()

    while len(out) < count:

        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        out.append(
            Target(
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return out


# =============================================================================
# POWER SUM POLYNOMIALS
# =============================================================================

def power_sum_polynomials(
    N: sp.Symbol,
    S: sp.Symbol,
    max_j: int,
) -> list[sp.Expr]:

    P = [
        sp.Integer(0)
    ] * (max_j + 1)

    P[0] = sp.Integer(2)

    if max_j >= 1:
        P[1] = S

    for j in range(
        2,
        max_j + 1,
    ):
        P[j] = sp.expand(
            S * P[j - 1]
            - N * P[j - 2]
        )

    return P


# =============================================================================
# NUMERIC SIGMA_r
# =============================================================================

def sigma_from_pq(
    p: int,
    q: int,
    r: int,
) -> int:

    return (
        (1 + p ** r)
        * (1 + q ** r)
    )


# =============================================================================
# PAPER DETECTOR
# =============================================================================

def detector_direct(
    k: int,
    ell: int,
    p: int,
    q: int,
) -> int:

    return (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )


# =============================================================================
# QUOTIENT POLYNOMIAL
# =============================================================================

def quotient_polynomial(
    k: int,
    ell: int,
) -> tuple[
    sp.Symbol,
    sp.Symbol,
    sp.Expr,
]:

    N, S = sp.symbols(
        "N S"
    )

    p, q = sp.symbols(
        "p q"
    )

    f = (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )

    symmetric, remainder, mapping = (
        sp.symmetrize(
            sp.expand(f),
            [p, q],
            formal=True,
        )
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            f"Symmetric reduction failed "
            f"for f_({k},{ell})."
        )

    symbols = list(
        dict(mapping).keys()
    )

    if len(symbols) != 2:
        raise ArithmeticError(
            "Unexpected symmetric mapping."
        )

    s1, s2 = symbols

    symmetric = sp.expand(
        symmetric.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )

    field = sp.QQ.frac_field(N)

    numerator = sp.Poly(
        symmetric,
        S,
        domain=field,
    )

    divisor = sp.Poly(
        S + 1,
        S,
        domain=field,
    )

    quotient, remainder_poly = (
        numerator.div(
            divisor
        )
    )

    if sp.expand(
        remainder_poly.as_expr()
    ) != 0:
        raise ArithmeticError(
            f"f_({k},{ell}) is not divisible "
            f"by S+1."
        )

    Q = sp.cancel(
        sp.expand(
            quotient.as_expr()
        )
    )

    if sp.expand(
        (S + 1) * Q - symmetric
    ) != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed "
            f"for f_({k},{ell})."
        )

    return (
        N,
        S,
        sp.factor(Q),
    )


# =============================================================================
# MOMENT DECOMPOSITION
# =============================================================================

def decompose_in_power_sums(
    Q: sp.Expr,
    N: sp.Symbol,
    S: sp.Symbol,
    max_moment: int,
) -> dict[int, sp.Expr]:

    field = sp.QQ.frac_field(N)

    poly = sp.Poly(
        sp.expand(Q),
        S,
        domain=field,
    )

    degree = poly.degree()

    if degree > max_moment:
        raise ArithmeticError(
            f"S degree {degree} exceeds "
            f"MAX_MOMENT={max_moment}."
        )

    P = power_sum_polynomials(
        N,
        S,
        max_moment,
    )

    remainder = sp.expand(Q)
    coefficients: dict[int, sp.Expr] = {}

    for d in range(
        degree,
        0,
        -1,
    ):

        current = sp.Poly(
            sp.expand(remainder),
            S,
            domain=field,
        )

        lead = sp.expand(
            current.coeff_monomial(
                S ** d
            )
        )

        coefficients[d] = sp.factor(
            lead
        )

        remainder = sp.expand(
            remainder
            - lead * P[d]
        )

    remainder = sp.expand(
        remainder
    )

    if sp.degree(
        remainder,
        S,
    ) > 0:
        raise ArithmeticError(
            "Residual S dependence."
        )

    if remainder != 0:
        coefficients[0] = sp.factor(
            remainder
        )

    return {
        j: sp.factor(c)
        for j, c in coefficients.items()
    }


# =============================================================================
# BUILD MOMENT SYSTEM
# =============================================================================

def build_system():

    N = sp.Symbol("N")

    decompositions = {}

    for detector in DETECTORS:

        k, ell = detector

        N2, S, Q = quotient_polynomial(
            k,
            ell,
        )

        d = decompose_in_power_sums(
            Q,
            N2,
            S,
            MAX_MOMENT,
        )

        decompositions[
            detector
        ] = d

    M_rows = []
    C0 = []

    for detector in DETECTORS:

        d = decompositions[
            detector
        ]

        C0.append(
            sp.factor(
                d.get(
                    0,
                    sp.Integer(0),
                )
            )
        )

        M_rows.append(
            [
                sp.factor(
                    d.get(
                        j,
                        sp.Integer(0),
                    )
                )
                for j in range(
                    1,
                    MAX_MOMENT + 1,
                )
            ]
        )

    M = sp.Matrix(
        M_rows
    )

    C0 = sp.Matrix(
        C0
    )

    return (
        N,
        M,
        C0,
        decompositions,
    )


# =============================================================================
# ROBUST SYMBOLIC MATRIX CHECK
# =============================================================================

def matrix_identity_failures(
    A: sp.Matrix,
    B: sp.Matrix,
) -> int:

    if A.shape != B.shape:
        return max(
            A.rows * A.cols,
            B.rows * B.cols,
        )

    failures = 0

    for i in range(A.rows):
        for j in range(A.cols):

            diff = sp.cancel(
                sp.together(
                    A[i, j]
                    - B[i, j]
                )
            )

            if sp.simplify(
                diff
            ) != 0:
                failures += 1

    return failures


# =============================================================================
# ROBUST EXACT INVERSE
# =============================================================================

def exact_inverse(
    M: sp.Matrix,
) -> sp.Matrix:

    if M.rows != M.cols:
        raise ArithmeticError(
            "Moment matrix must be square."
        )

    determinant = sp.factor(
        M.det()
    )

    if sp.simplify(
        determinant
    ) == 0:
        raise ArithmeticError(
            "Moment matrix is identically singular."
        )

    # Compute once.
    inverse = M.inv()

    # Simplify every entry separately.
    inverse = inverse.applyfunc(
        lambda x: sp.factor(
            sp.cancel(x)
        )
    )

    # Robust symbolic validation.
    product = (
        M * inverse
    )

    failures = matrix_identity_failures(
        product,
        sp.eye(M.rows),
    )

    if failures:
        raise ArithmeticError(
            "M inverse verification failed "
            f"after simplification: "
            f"{failures} nonzero entries."
        )

    return inverse


# =============================================================================
# ORACLE Q VECTOR
# =============================================================================

def oracle_q_vector(
    t: Target,
) -> sp.Matrix:

    values = []

    for k, ell in DETECTORS:

        f = detector_direct(
            k,
            ell,
            t.p,
            t.q,
        )

        values.append(
            sp.Rational(
                f,
                t.s + 1,
            )
        )

    return sp.Matrix(
        values
    )


# =============================================================================
# P1 RESOLVENT
# =============================================================================

def build_p1_resolvent(
    N: sp.Symbol,
    M: sp.Matrix,
    C0: sp.Matrix,
):

    inverse = exact_inverse(
        M
    )

    Qs = sp.symbols(
        "Q13 Q15 Q17 Q35 Q37 Q57"
    )

    P = (
        inverse
        * (
            sp.Matrix(Qs)
            - C0
        )
    )

    P1 = sp.factor(
        sp.cancel(
            P[0]
        )
    )

    return (
        inverse,
        Qs,
        P1,
    )


# =============================================================================
# SIGMA-BASIS REWRITE
# =============================================================================

def derive_p1_sigma_expression(
    N: sp.Symbol,
    inverse: sp.Matrix,
    C0: sp.Matrix,
    decompositions,
):

    Sigma = sp.symbols(
        "Sigma1:7"
    )

    sigma_symbols = {
        j: Sigma[j - 1]
        for j in range(
            1,
            7,
        )
    }

    Q_sigma = []

    for detector in DETECTORS:

        d = decompositions[
            detector
        ]

        expr = sp.Integer(0)

        for j, coefficient in d.items():

            if j == 0:

                expr += coefficient

            else:

                expr += (
                    coefficient
                    * (
                        sigma_symbols[j]
                        - 1
                        - N ** j
                    )
                )

        Q_sigma.append(
            sp.expand(expr)
        )

    P_vector = (
        inverse
        * (
            sp.Matrix(Q_sigma)
            - C0
        )
    )

    P1_sigma = sp.factor(
        sp.cancel(
            sp.expand(
                P_vector[0]
            )
        )
    )

    return (
        Sigma,
        P1_sigma,
        P_vector,
    )


# =============================================================================
# SIGMA DEPENDENCY
# =============================================================================

def sigma_dependency(
    expr: sp.Expr,
    Sigma,
) -> list[int]:

    deps = []

    for j, symbol in enumerate(
        Sigma,
        start=1,
    ):

        derivative = sp.diff(
            expr,
            symbol,
        )

        if sp.simplify(
            derivative
        ) != 0:
            deps.append(j)

    return deps


# =============================================================================
# NUMERICAL VALIDATION
# =============================================================================

def validate_sigma_formula(
    expr: sp.Expr,
    Sigma,
    N: sp.Symbol,
    targets: list[Target],
) -> int:

    failures = 0

    for t in targets:

        substitutions = {
            N: t.n,
        }

        for j, symbol in enumerate(
            Sigma,
            start=1,
        ):

            substitutions[
                symbol
            ] = sigma_from_pq(
                t.p,
                t.q,
                j,
            )

        value = sp.cancel(
            sp.together(
                expr.subs(
                    substitutions
                )
            )
        )

        if sp.simplify(
            value - t.s
        ) != 0:
            failures += 1

    return failures


# =============================================================================
# MAIN
# =============================================================================

def main():

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 103R")
    print("PATCHED INVERSE-RESOLVENT ARITHMETIC COLLAPSE")
    print("EXACT PAPER DETECTOR -> SIGMA BASIS")
    print("STRICT TARGET HOLDOUT")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # =========================================================================
    # PRIME POPULATION
    # =========================================================================

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)

    print(
        f"prime population = "
        f"{len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter()-t0:.6f}s"
    )

    targets = generate_targets(
        primes,
        TARGET_COUNT,
    )

    for i, t in enumerate(
        targets[:24],
        start=1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if TARGET_COUNT > 24:
        print(
            "... remaining targets omitted"
        )

    train = targets[
        :TRAIN_COUNT
    ]

    test = targets[
        TRAIN_COUNT:
    ]

    print()
    print("2. TARGET HOLDOUT")
    print("-" * 78)

    print(
        f"training targets = "
        f"{len(train)}"
    )

    print(
        f"test targets = "
        f"{len(test)}"
    )

    # =========================================================================
    # BUILD SYSTEM
    # =========================================================================

    print()
    print("3. BUILDING PAPER MOMENT SYSTEM")
    print("-" * 78)

    N, M, C0, decompositions = (
        build_system()
    )

    for detector in DETECTORS:

        print()
        print(
            f"Q_{detector}"
        )

        for j in sorted(
            decompositions[
                detector
            ]
        ):

            print(
                f"  P{j}: "
                f"{sp.factor(
                    decompositions[
                        detector
                    ][j]
                )}"
            )

    print()
    print("status = PASS")

    # =========================================================================
    # DETERMIMANT + INVERSE
    # =========================================================================

    print()
    print("4. EXACT RATIONAL INVERSE")
    print("-" * 78)

    determinant = sp.factor(
        M.det()
    )

    print(
        f"det(M) = "
        f"{determinant}"
    )

    start = time.perf_counter()

    inverse, Qs, P1_direct = (
        build_p1_resolvent(
            N,
            M,
            C0,
        )
    )

    print(
        f"inverse time = "
        f"{time.perf_counter()-start:.6f}s"
    )

    print()
    print("P1 in Q basis:")
    print(
        P1_direct
    )

    # =========================================================================
    # SIGMA REWRITE
    # =========================================================================

    print()
    print("5. REWRITE P1 IN SIGMA_r BASIS")
    print("-" * 78)

    Sigma, P1_sigma, P_vector = (
        derive_p1_sigma_expression(
            N,
            inverse,
            C0,
            decompositions,
        )
    )

    print(
        "P1_sigma ="
    )

    print(
        P1_sigma
    )

    deps = sigma_dependency(
        P1_sigma,
        Sigma,
    )

    print()
    print(
        f"sigma dependencies = "
        f"{deps}"
    )

    # =========================================================================
    # VALIDATION
    # =========================================================================

    print()
    print("6. SIGMA-BASIS VALIDATION")
    print("-" * 78)

    all_failures = (
        validate_sigma_formula(
            P1_sigma,
            Sigma,
            N,
            targets,
        )
    )

    train_failures = (
        validate_sigma_formula(
            P1_sigma,
            Sigma,
            N,
            train,
        )
    )

    test_failures = (
        validate_sigma_formula(
            P1_sigma,
            Sigma,
            N,
            test,
        )
    )

    print(
        f"all failures  = "
        f"{all_failures}/{len(targets)}"
    )

    print(
        f"train failures = "
        f"{train_failures}/{len(train)}"
    )

    print(
        f"test failures  = "
        f"{test_failures}/{len(test)}"
    )

    # =========================================================================
    # SIGMA COEFFICIENTS
    # =========================================================================

    print()
    print("7. SIGMA COEFFICIENT SPECTRUM")
    print("-" * 78)

    for j, symbol in enumerate(
        Sigma,
        start=1,
    ):

        coefficient = sp.factor(
            sp.cancel(
                sp.diff(
                    P1_sigma,
                    symbol,
                )
            )
        )

        if coefficient == 0:
            continue

        num, den = sp.fraction(
            sp.together(
                coefficient
            )
        )

        num_poly = sp.Poly(
            sp.expand(num),
            N,
        )

        den_poly = sp.Poly(
            sp.expand(den),
            N,
        )

        print()
        print(
            f"sigma_{j} coefficient = "
            f"{coefficient}"
        )

        print(
            f"  numerator degree = "
            f"{num_poly.degree()}"
        )

        print(
            f"  denominator degree = "
            f"{den_poly.degree()}"
        )

    # =========================================================================
    # TAUTOLOGY TEST
    # =========================================================================

    print()
    print("8. TAUTOLOGY TEST")
    print("-" * 78)

    sigma1 = Sigma[0]

    difference = sp.factor(
        sp.cancel(
            P1_sigma
            - (
                sigma1
                - N
                - 1
            )
        )
    )

    print(
        "P1 - (sigma1 - N - 1) ="
    )

    print(
        difference
    )

    if sp.simplify(
        difference
    ) == 0:

        print(
            "RESULT = EXACT SIGMA1 TAUTOLOGY"
        )

    else:

        print(
            "RESULT = NONTRIVIAL SIGMA COMBINATION"
        )

    # =========================================================================
    # LOWER-WEIGHT COLLAPSE
    # =========================================================================

    print()
    print("9. LOWER-WEIGHT COLLAPSE")
    print("-" * 78)

    for cutoff in range(
        1,
        7,
    ):

        substitutions = {
            Sigma[j - 1]: 0
            for j in range(
                cutoff + 1,
                7,
            )
        }

        truncated = sp.factor(
            sp.cancel(
                P1_sigma.subs(
                    substitutions
                )
            )
        )

        residual = sp.factor(
            sp.cancel(
                truncated
                - (
                    sigma1
                    - N
                    - 1
                )
            )
        )

        print()
        print(
            f"cutoff sigma_{cutoff}:"
        )

        print(
            f"residual = {residual}"
        )

    # =========================================================================
    # ORACLE TEST
    # =========================================================================

    print()
    print("10. HELD-OUT ORACLE SANITY")
    print("-" * 78)

    for i, t in enumerate(
        test[:5],
        start=1,
    ):

        substitutions = {
            N: t.n,
        }

        for j, symbol in enumerate(
            Sigma,
            start=1,
        ):

            substitutions[
                symbol
            ] = sigma_from_pq(
                t.p,
                t.q,
                j,
            )

        recovered = sp.simplify(
            P1_sigma.subs(
                substitutions
            )
        )

        print(
            f"target "
            f"{TRAIN_COUNT+i:2d}: "
            f"recovered={recovered} "
            f"true_s={t.s} "
            f"ok={recovered == t.s}"
        )

    # =========================================================================
    # FINAL
    # =========================================================================

    print()
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
103R is the patched continuation of 102R4-P.

The previous failure was only due to fragile symbolic matrix equality
checking. The inverse is now verified entry-by-entry using simplify().

The experiment asks whether the exact P1 inverse has an arithmetic
collapse in the divisor-sum basis:

    P_r = sigma_r(N) - N^r - 1.

Interpretation:

A. NONTRIVIAL LOW-WEIGHT COLLAPSE
   A compact combination of several sigma_r appears.

B. EXACT SIGMA1 TAUTOLOGY
   The whole construction reduces to

       P1 = sigma1(N) - N - 1.

   This confirms consistency but gives no new N-only algorithm.

C. HIGHER-WEIGHT EXPRESSION
   The inverse remains algebraically exact but computationally no
   simpler than the original detector system.

The crucial distinction remains:

    algebraic identity
        !=
    efficient evaluation from N alone.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter()-total_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 103R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()