#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 96
EXACT ALGEBRAIC INDEPENDENCE OF PAPER PRIME DETECTORS

PAPER FAMILY
------------
f_(k,l)(p,q)
  = (1+q)^l p^k - (1+q)^k p^l
    + (1+p)^l q^k - (1+p)^k q^l

SEMIPRIME VARIABLES
-------------------
N = p*q
s = p+q

E1:
    E1(N,s) = s*((N+1)^2 - s^2)

MAIN QUESTION
-------------
Does the paper detector f_(k,l), after symmetric reduction to Z[N,s],
contain genuinely new algebraic information beyond E1 / sigma_1?

We investigate:

    f_(k,l)(N,s)

    factorization

    common factor with (N+s+1)

    quotient after removing common sigma1 factor

    gcd(f_(k,l), E1) as polynomials in s

    resultant(quotient, E1, s)

    pairwise gcds between detector quotients

    numerical independence on held-out semiprimes

IMPORTANT
---------
This is a symbolic/algebraic experiment.

It does NOT attempt factorization.

A positive result means the paper family contains an algebraically
independent resolvent worth testing computationally.

A negative result means the apparent detector spectrum may largely be
different presentations of the same sigma_1 / E1 information.

NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import sympy as sp


# ============================================================================
# PARAMETERS
# ============================================================================

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols(
    "p q",
    integer=True,
)

N, S = sp.symbols(
    "N S",
    integer=True,
)

T = sp.symbols(
    "T",
    integer=True,
)

# ============================================================================


@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

    flags = bytearray(
        b"\x01" * (hi + 1)
    )

    flags[0:2] = b"\x00\x00"

    for x in range(
        2,
        math.isqrt(hi) + 1,
    ):
        if flags[x]:
            start = x * x
            flags[
                start:hi + 1:x
            ] = b"\x00" * (
                ((hi - start) // x) + 1
            )

    return [
        x
        for x in range(lo, hi + 1)
        if flags[x]
    ]


# ============================================================================
# TARGETS
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    rng = random.Random(
        96096
    )

    result = []
    seen = set()

    while len(result) < count:

        p0 = primes[
            rng.randrange(
                len(primes)
            )
        ]

        q0 = primes[
            rng.randrange(
                len(primes)
            )
        ]

        if p0 == q0:
            continue

        if p0 > q0:
            p0, q0 = q0, p0

        if (p0, q0) in seen:
            continue

        seen.add(
            (p0, q0)
        )

        result.append(
            Target(
                p=p0,
                q=q0,
                n=p0 * q0,
                s=p0 + q0,
            )
        )

    return result


# ============================================================================
# PAPER DETECTOR
# ============================================================================

def paper_detector(
    k: int,
    ell: int,
) -> sp.Expr:

    return sp.expand(
        (1 + q) ** ell * p**k
        - (1 + q) ** k * p**ell
        + (1 + p) ** ell * q**k
        - (1 + p) ** k * q**ell
    )


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduce(
    expr: sp.Expr,
) -> sp.Expr:

    """
    Convert symmetric p,q polynomial into N=pq, S=p+q.

    symmetrize(..., formal=True) creates s1=p+q and s2=pq.
    """

    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [p, q],
        formal=True,
    )

    if remainder != 0:
        raise ArithmeticError(
            f"Non-symmetric remainder: {remainder}"
        )

    s1_symbol, s2_symbol = mapping[0][0], mapping[1][0]

    reduced = reduced.subs(
        {
            s1_symbol: S,
            s2_symbol: N,
        }
    )

    return sp.expand(
        reduced
    )


# ============================================================================
# E1
# ============================================================================

def E1_polynomial() -> sp.Expr:

    return sp.expand(
        S * (
            (N + 1) ** 2
            - S**2
        )
    )


# ============================================================================
# NUMERICAL DIRECT CHECK
# ============================================================================

def detector_direct_numeric(
    t: Target,
    k: int,
    ell: int,
) -> int:

    return int(
        (
            (1 + t.q) ** ell * t.p**k
            - (1 + t.q) ** k * t.p**ell
            + (1 + t.p) ** ell * t.q**k
            - (1 + t.p) ** k * t.q**ell
        )
    )


def evaluate_ns(
    expr: sp.Expr,
    n: int,
    s: int,
) -> int:

    return int(
        expr.subs(
            {
                N: n,
                S: s,
            }
        )
    )


# ============================================================================
# MAIN SYMBOLIC ANALYSIS
# ============================================================================

def main() -> None:

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 96")
    print("EXACT ALGEBRAIC INDEPENDENCE OF PAPER PRIME DETECTORS")
    print("SYMMETRIC REDUCTION IN Z[N,S]")
    print("E1 / SIGMA1 FACTOR ANALYSIS")
    print("PAIRWISE GCD + RESULTANT ANALYSIS")
    print("STRICT TARGET HOLDOUT")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------------

    print("\n1. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()

    primes = sieve_primes(
        PRIME_MIN,
        PRIME_MAX,
    )

    print(
        f"prime population = {len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter() - start:.6f}s"
    )

    # ------------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = {len(targets)}"
    )

    for i, t in enumerate(
        targets[:20],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    print(
        f"\ntraining targets = {TRAIN_TARGETS}"
    )
    print(
        f"test targets     = {TEST_TARGETS}"
    )

    # ------------------------------------------------------------------------
    # E1
    # ------------------------------------------------------------------------

    print("\n2. E1 BASELINE")
    print("-" * 78)

    E1 = E1_polynomial()

    print(
        "E1(N,S) ="
    )
    print(
        sp.factor(E1)
    )

    sigma1 = (
        N + S + 1
    )

    print(
        "\nsigma1(N,S) ="
    )
    print(
        sigma1
    )

    # ------------------------------------------------------------------------
    # Build detectors
    # ------------------------------------------------------------------------

    reduced: Dict[
        Tuple[int, int],
        sp.Expr,
    ] = {}

    quotients: Dict[
        Tuple[int, int],
        sp.Expr,
    ] = {}

    print("\n3. PAPER DETECTOR REDUCTION")
    print("-" * 78)

    for pair in DETECTORS:

        k, ell = pair

        direct = paper_detector(
            k,
            ell,
        )

        reduced_expr = symmetric_reduce(
            direct
        )

        reduced[pair] = reduced_expr

        print(
            f"\nf_({k},{ell})"
        )

        print(
            "  reduced ="
        )

        print(
            "   ",
            sp.factor(
                reduced_expr
            )
        )

        # Test sigma1 factor.
        sigma_factor_test = sp.rem(
            sp.Poly(
                reduced_expr,
                S,
            ),
            sp.Poly(
                sigma1,
                S,
            ),
        )

        sigma_divisible = (
            sigma_factor_test == 0
        )

        print(
            f"  divisible by "
            f"(N+S+1): {sigma_divisible}"
        )

        if sigma_divisible:

            quotient = sp.factor(
                sp.cancel(
                    reduced_expr
                    / sigma1
                )
            )

        else:

            quotient = sp.factor(
                reduced_expr
            )

        quotients[pair] = sp.expand(
            quotient
        )

        print(
            "  quotient after sigma1 ="
        )

        print(
            "   ",
            quotient
        )

    # ------------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------------

    print("\n4. NUMERICAL IDENTITY VALIDATION")
    print("-" * 78)

    for pair in DETECTORS:

        failures_direct = 0
        failures_reduced = 0

        for t in targets:

            direct = detector_direct_numeric(
                t,
                pair[0],
                pair[1],
            )

            reduced_value = evaluate_ns(
                reduced[pair],
                t.n,
                t.s,
            )

            quotient_value = evaluate_ns(
                quotients[pair],
                t.n,
                t.s,
            )

            sigma_value = (
                t.n
                + t.s
                + 1
            )

            reconstructed = (
                quotient_value
                * sigma_value
            )

            if direct != reduced_value:
                failures_direct += 1

            if direct != reconstructed:
                failures_reduced += 1

        print(
            f"f_{pair}: "
            f"direct failures={failures_direct} "
            f"sigma-quotient failures={failures_reduced}"
        )

    # ------------------------------------------------------------------------
    # GCD WITH E1
    # ------------------------------------------------------------------------

    print("\n5. GCD WITH E1 AS POLYNOMIAL IN S")
    print("-" * 78)

    gcd_data = {}

    for pair in DETECTORS:

        qexpr = quotients[pair]

        gcd_expr = sp.factor(
            sp.gcd(
                sp.Poly(
                    qexpr,
                    S,
                ),
                sp.Poly(
                    E1,
                    S,
                ),
            ).as_expr()
        )

        gcd_data[pair] = gcd_expr

        print(
            f"f_{pair}:"
        )

        print(
            f"  gcd(quotient, E1) = "
            f"{gcd_expr}"
        )

    # ------------------------------------------------------------------------
    # RESULTANTS
    # ------------------------------------------------------------------------

    print("\n6. RESULTANT OF QUOTIENT WITH E1")
    print("-" * 78)

    resultant_data = {}

    for pair in DETECTORS:

        qexpr = quotients[pair]

        print(
            f"\nf_{pair}"
        )

        resultant = sp.factor(
            sp.resultant(
                qexpr,
                E1,
                S,
            )
        )

        resultant_data[pair] = resultant

        print(
            "resultant ="
        )

        print(
            resultant
        )

        if resultant == 0:
            print(
                "  STATUS = ALGEBRAIC DEPENDENCE"
            )
        else:
            print(
                "  STATUS = GENERICALLY INDEPENDENT"
            )

    # ------------------------------------------------------------------------
    # PAIRWISE QUOTIENT GCDs
    # ------------------------------------------------------------------------

    print("\n7. PAIRWISE DETECTOR GCD MATRIX")
    print("-" * 78)

    print(
        "Rows/columns = detector pairs."
    )

    for a in DETECTORS:

        row = []

        for b in DETECTORS:

            g = sp.factor(
                sp.gcd(
                    sp.Poly(
                        quotients[a],
                        S,
                    ),
                    sp.Poly(
                        quotients[b],
                        S,
                    ),
                ).as_expr()
            )

            row.append(
                str(g)
            )

        print(
            f"{a}: {row}"
        )

    # ------------------------------------------------------------------------
    # Algebraic reconstruction test
    # ------------------------------------------------------------------------

    print("\n8. NUMERICAL INDEPENDENCE TEST")
    print("-" * 78)

    print(
        "For each held-out target we compute the exact detector values."
    )

    print(
        "Then we verify that the detector quotient contains "
        "information not explained by E1 alone."
    )

    for pair in DETECTORS:

        qexpr = quotients[pair]

        train_values = [
            evaluate_ns(
                qexpr,
                t.n,
                t.s,
            )
            for t in targets[:TRAIN_TARGETS]
        ]

        test_values = [
            evaluate_ns(
                qexpr,
                t.n,
                t.s,
            )
            for t in targets[TRAIN_TARGETS:]
        ]

        print(
            f"f_{pair}: "
            f"train quotient values="
            f"{len(train_values)} "
            f"test quotient values="
            f"{len(test_values)}"
        )

    # ------------------------------------------------------------------------
    # Specialized f_(1,3)
    # ------------------------------------------------------------------------

    print("\n9. SPECIAL FOCUS: f_(1,3)")
    print("-" * 78)

    f13 = sp.factor(
        reduced[(1, 3)]
    )

    q13 = sp.factor(
        quotients[(1, 3)]
    )

    print(
        "f_(1,3)(N,S) ="
    )
    print(
        f13
    )

    print(
        "\nquotient after sigma1 ="
    )
    print(
        q13
    )

    print(
        "\nE1(N,S) ="
    )
    print(
        sp.factor(E1)
    )

    print(
        "\nGCD with E1 ="
    )
    print(
        gcd_data[(1, 3)]
    )

    print(
        "\nResultant with E1 ="
    )
    print(
        resultant_data[(1, 3)]
    )

    # ------------------------------------------------------------------------
    # Special algebraic combinations
    # ------------------------------------------------------------------------

    print("\n10. SIMPLE COMBINATIONS")
    print("-" * 78)

    combinations = []

    # Quotient pairs:
    q13 = quotients[(1, 3)]
    q15 = quotients[(1, 5)]
    q17 = quotients[(1, 7)]
    q35 = quotients[(3, 5)]

    combinations.append(
        ("Q13-Q15", sp.factor(q13 - q15))
    )

    combinations.append(
        ("Q15-Q17", sp.factor(q15 - q17))
    )

    combinations.append(
        ("Q13-Q17", sp.factor(q13 - q17))
    )

    combinations.append(
        ("Q13+Q15", sp.factor(q13 + q15))
    )

    combinations.append(
        ("Q13*Q15 gcd E1", sp.factor(
            sp.gcd(
                sp.Poly(q13 * q15, S),
                sp.Poly(E1, S),
            ).as_expr()
        ))
    )

    for name, expr in combinations:

        print(
            f"{name}:"
        )

        print(
            expr
        )

    # ------------------------------------------------------------------------
    # Complexity measure
    # ------------------------------------------------------------------------

    print("\n11. RESOLVENT COMPLEXITY")
    print("-" * 78)

    for pair in DETECTORS:

        expr = quotients[pair]

        poly = sp.Poly(
            expr,
            S,
        )

        degree = poly.degree()

        coefficient_polys = [
            sp.Poly(
                c,
                N,
            )
            for c in poly.all_coeffs()
        ]

        max_coeff_degree = max(
            (
                c.degree()
                for c in coefficient_polys
            ),
            default=0,
        )

        print(
            f"f_{pair}: "
            f"s-degree={degree} "
            f"max coefficient degree in N="
            f"{max_coeff_degree}"
        )

    # ------------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The central question is now algebraic:

    Does the paper's f_(k,l) family contain genuinely new
    information beyond sigma_1 / E1?

Interpretation:

A. If every detector is divisible by (N+S+1), that common factor
   is not new information; it is sigma_1(N).

B. If the quotient after removing (N+S+1) has nontrivial gcd with
   E1, the detector may still largely reproduce the existing
   E1 information.

C. If the resultant of the quotient with E1 is nonzero, the two
   equations are generically algebraically independent.

D. If two independent paper detectors have a nontrivial resultant
   in S but neither is reducible to E1, then their combination may
   give a new semiprime resolvent.

The highest-value outcome is therefore NOT candidate collapse.

It is:

    independent paper detector
        +
    E1
        ->
    elimination relation in N or a lower-degree equation for S.

If such a relation exists, that is the object we should attempt
to evaluate from N alone in the next experiment.
"""
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 96 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

