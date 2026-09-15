#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 98
PAPER QUOTIENT-ALGEBRA COLLAPSE

REMOVE UNIVERSAL (S+1) FACTOR
USE Q_(1,3) AS QUADRATIC GENERATOR

SEMIPRIME VARIABLES
-------------------
N = p*q
S = p+q

PAPER DETECTORS
---------------
f_(k,l)(p,q)
  = (1+q)^l p^k - (1+q)^k p^l
    + (1+p)^l q^k - (1+p)^k q^l

For all tested paper detectors:

    f_(k,l)(N,S) = (S+1) Q_(k,l)(N,S)

Special low-weight quotient:

    Q_(1,3) = 6N - S^2 + S

Therefore:

    S^2 - S = 6N - Q_(1,3)

MAIN QUESTION
-------------
Does the rest of the paper family collapse algebraically once Q_(1,3)
is treated as a known quadratic invariant?

For each Q_(k,l), reduce modulo:

    R = S^2 - S - 6N + Q13

The remainder has degree <= 1 in S:

    Q_(k,l) = A_k(N,Q13)*S + B_k(N,Q13)

If A=0, the quotient is completely determined by N and Q13.

If several quotients have linear remainders, search combinations that
eliminate S.

The strongest outcome would be:

    H(N,Q13,Q_other) = F(N)

or

    S = rational_function(N,Q13,Q_other)

or, even better,

    Q13 = F(N)

which would expose the hidden resolvent directly from N.

NO million-candidate S scan
NO CRT sieve
NO SKLEARN
NO CSV
==============================================================================
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from itertools import combinations, product
from typing import Dict, List, Sequence, Tuple

import sympy as sp


# =============================================================================
# PARAMETERS
# =============================================================================

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

DETECTORS: Tuple[Tuple[int, int], ...] = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

SMALL_COEFFS = (-3, -2, -1, 1, 2, 3)

P, Q = sp.symbols("P Q")
N, S, Z = sp.symbols("N S Z")


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

def sieve_primes(lo: int, hi: int) -> List[int]:
    flags = bytearray(b"\x01") * (hi + 1)
    flags[0:2] = b"\x00\x00"

    for x in range(2, math.isqrt(hi) + 1):
        if flags[x]:
            start = x * x
            flags[start:hi + 1:x] = b"\x00" * (
                ((hi - start) // x) + 1
            )

    return [x for x in range(lo, hi + 1) if flags[x]]


# =============================================================================
# TARGET GENERATION
# =============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:
    rng = random.Random(98098)

    result: List[Target] = []
    seen = set()

    while len(result) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        result.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return result


# =============================================================================
# PAPER DETECTOR
# =============================================================================

def paper_detector(
    k: int,
    ell: int,
) -> sp.Expr:
    return sp.expand(
        (1 + Q) ** ell * P**k
        - (1 + Q) ** k * P**ell
        + (1 + P) ** ell * Q**k
        - (1 + P) ** k * Q**ell
    )


# =============================================================================
# SYMMETRIC REDUCTION
# =============================================================================

def symmetric_reduce(expr: sp.Expr) -> sp.Expr:
    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [P, Q],
        formal=True,
    )

    if remainder != 0:
        raise ArithmeticError(
            f"Symmetric reduction remainder != 0: {remainder}"
        )

    s1 = mapping[0][0]
    s2 = mapping[1][0]

    return sp.expand(
        reduced.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )


# =============================================================================
# DETECTOR QUOTIENTS
# =============================================================================

def build_quotients() -> Dict[Tuple[int, int], sp.Expr]:
    result: Dict[Tuple[int, int], sp.Expr] = {}

    for pair in DETECTORS:
        f = symmetric_reduce(
            paper_detector(*pair)
        )

        quotient, remainder = sp.div(
            sp.Poly(f, S),
            sp.Poly(S + 1, S),
        )

        if sp.expand(remainder.as_expr()) != 0:
            raise ArithmeticError(
                f"Detector f_{pair} is not divisible by S+1."
            )

        result[pair] = sp.factor(
            quotient.as_expr()
        )

    return result


# =============================================================================
# E1
# =============================================================================

def E1_polynomial() -> sp.Expr:
    return sp.expand(
        S * ((N + 1) ** 2 - S**2)
    )


# =============================================================================
# Q13 QUADRATIC RELATION
# =============================================================================

def q13_polynomial(quotients: Dict[Tuple[int, int], sp.Expr]) -> sp.Expr:
    return sp.factor(
        quotients[(1, 3)]
    )


def q13_relation(q13: sp.Expr) -> sp.Expr:
    # q13 = 6N - S^2 + S
    # therefore S^2 - S - 6N + q13 = 0
    Zsymbol = Z
    return sp.expand(
        S**2 - S - 6*N + Zsymbol
    )


# =============================================================================
# REDUCE MODULO Q13 RELATION
# =============================================================================

def reduce_mod_q13(
    expr: sp.Expr,
    q13_symbol: sp.Symbol = Z,
) -> sp.Expr:
    """
    Reduce expr modulo:
        S^2 - S - 6N + Z = 0

    Treat Z as the symbolic value of Q13.
    """

    relation = sp.Poly(
        S**2 - S - 6*N + q13_symbol,
        S,
        domain="EX",
    )

    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain="EX",
    )

    rem = sp.rem(
        poly,
        relation,
    )

    return sp.factor(
        rem.as_expr()
    )


# =============================================================================
# EXPLICIT Q13 SUBSTITUTION
# =============================================================================

def q13_substitute_to_check(
    expr: sp.Expr,
    q13_expr: sp.Expr,
) -> sp.Expr:
    """
    Replace Z by the actual Q13(N,S).
    """
    return sp.factor(
        sp.expand(
            expr.subs(
                Z,
                q13_expr,
            )
        )
    )


# =============================================================================
# NUMERIC EVALUATION
# =============================================================================

def eval_ns(
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


def direct_detector(
    t: Target,
    pair: Tuple[int, int],
) -> int:
    k, ell = pair

    return int(
        (1 + t.q) ** ell * t.p**k
        - (1 + t.q) ** k * t.p**ell
        + (1 + t.p) ** ell * t.q**k
        - (1 + p) ** k * t.q**ell
    )


# =============================================================================
# CORRECT DIRECT DETECTOR
# =============================================================================

def direct_detector_correct(
    t: Target,
    pair: Tuple[int, int],
) -> int:
    k, ell = pair

    return int(
        (1 + t.q) ** ell * t.p**k
        - (1 + t.q) ** k * t.p**ell
        + (1 + t.p) ** ell * t.q**k
        - (1 + t.p) ** k * t.q**ell
    )


# =============================================================================
# VALIDATION
# =============================================================================

def validate_quotients(
    targets: Sequence[Target],
    quotients: Dict[Tuple[int, int], sp.Expr],
) -> None:

    failures = 0

    for pair, quotient in quotients.items():

        for t in targets:

            expected = direct_detector_correct(
                t,
                pair,
            )

            reconstructed = (
                t.s + 1
            ) * eval_ns(
                quotient,
                t.n,
                t.s,
            )

            if expected != reconstructed:
                failures += 1

    print(
        f"quotient reconstruction failures = "
        f"{failures}"
    )

    if failures:
        raise ArithmeticError(
            "Quotient validation failed."
        )


# =============================================================================
# LINEAR DECOMPOSITION
# =============================================================================

def split_linear_in_S(
    expr: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Return:

        expr = A(N,Z)*S + B(N,Z)

    assuming expr has S-degree <= 1.
    """

    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain="EX",
    )

    degree = poly.degree()

    if degree > 1:
        raise ValueError(
            f"Expression is not linear in S: degree={degree}"
        )

    if degree == 0:
        return sp.Integer(0), sp.factor(
            poly.as_expr()
        )

    coeffs = poly.all_coeffs()

    return (
        sp.factor(coeffs[0]),
        sp.factor(coeffs[1]),
    )


# =============================================================================
# COMPLEXITY
# =============================================================================

def complexity(expr: sp.Expr) -> Tuple[int, int, int, int]:
    poly = sp.Poly(
        sp.expand(expr),
        N,
        S,
    )

    total_degree = poly.total_degree()

    s_degree = sp.Poly(
        sp.expand(expr),
        S,
    ).degree()

    n_degree = max(
        (
            sp.Poly(c, N).degree()
            for c in sp.Poly(
                sp.expand(expr),
                S,
            ).all_coeffs()
        ),
        default=0,
    )

    terms = len(
        poly.terms()
    )

    return (
        total_degree,
        s_degree,
        n_degree,
        terms,
    )


# =============================================================================
# MAIN REDUCTION TABLE
# =============================================================================

def reduction_table(
    quotients: Dict[Tuple[int, int], sp.Expr],
    q13_expr: sp.Expr,
) -> Dict[
    Tuple[int, int],
    Tuple[sp.Expr, sp.Expr],
]:

    results = {}

    print("\n3. Q13 ALGEBRA REDUCTION")
    print("-" * 78)

    print(
        "Relation:"
    )
    print(
        "  S^2 - S - 6N + Q13 = 0"
    )

    for pair in DETECTORS:

        q = quotients[pair]

        reduced = reduce_mod_q13(
            q
        )

        A, B = split_linear_in_S(
            reduced
        )

        results[pair] = (
            A,
            B,
        )

        print(
            f"\nQ_{pair}"
        )
        print(
            f"  original = {sp.factor(q)}"
        )
        print(
            f"  reduced  = {sp.factor(reduced)}"
        )
        print(
            f"  A(N,Q13) = {A}"
        )
        print(
            f"  B(N,Q13) = {B}"
        )

        # Verify that the reduced form equals the original detector
        # quotient on the exact semiprime targets.
        failures = 0

        for t in TARGETS_GLOBAL:
            q13_value = eval_ns(
                q13_expr,
                t.n,
                t.s,
            )

            reduced_value = int(
                reduced.subs(
                    {
                        N: t.n,
                        Z: q13_value,
                        S: t.s,
                    }
                )
            )

            original_value = eval_ns(
                q,
                t.n,
                t.s,
            )

            if reduced_value != original_value:
                failures += 1

        print(
            f"  reduction validation = "
            f"{len(TARGETS_GLOBAL)-failures}/"
            f"{len(TARGETS_GLOBAL)}"
        )

    return results


# =============================================================================
# SEARCH N-ONLY CANCELLATIONS
# =============================================================================

def search_linear_elimination(
    linear_forms: Dict[
        Tuple[int, int],
        Tuple[sp.Expr, sp.Expr],
    ],
) -> None:

    print(
        "\n4. LINEAR-IN-S ELIMINATION"
    )
    print("-" * 78)

    rows = []

    for a, b in combinations(
        DETECTORS,
        2,
    ):

        A1, B1 = linear_forms[a]
        A2, B2 = linear_forms[b]

        determinant = sp.factor(
            A1 * B2 - A2 * B1
        )

        eliminated = sp.factor(
            A1 * B2 - A2 * B1
        )

        numerator = sp.factor(
            A1 * B2 - A2 * B1
        )

        print(
            f"\nPAIR Q_{a}, Q_{b}"
        )

        print(
            f"  A1 = {A1}"
        )
        print(
            f"  B1 = {B1}"
        )
        print(
            f"  A2 = {A2}"
        )
        print(
            f"  B2 = {B2}"
        )

        print(
            "\n  determinant ="
        )
        print(
            determinant
        )

        # If two equations are:
        #
        #   Q_a = A1*S + B1
        #   Q_b = A2*S + B2
        #
        # then:
        #
        #   A2*Q_a - A1*Q_b = A2*B1 - A1*B2
        #
        # This is an N,Z-only relation if Q_a and Q_b are considered
        # observed invariants.
        #
        # We record its form.

        lhs = sp.factor(
            A2 * Z - A1 * Z
        )

        rows.append(
            (
                complexity(determinant),
                a,
                b,
                determinant,
            )
        )

    rows.sort(
        key=lambda x: x[0]
    )

    print(
        "\nBEST DETERMINANTS"
    )

    for row in rows[:20]:

        print(
            f"pair={row[1]}/{row[2]} "
            f"complexity={row[0]}"
        )
        print(
            row[3]
        )


# =============================================================================
# DIRECT COMBINATIONS OF Q's
# =============================================================================

def search_q_combinations(
    quotients: Dict[Tuple[int, int], sp.Expr],
    q13_expr: sp.Expr,
) -> None:

    print(
        "\n5. SMALL-INTEGER Q-COMBINATION SEARCH"
    )
    print("-" * 78)

    candidates = []

    names = list(
        DETECTORS
    )

    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):

            qa = quotients[a]
            qb = quotients[names[j]]

            for ca in SMALL_COEFFS:
                for cb in SMALL_COEFFS:

                    expr = sp.factor(
                        ca * qa
                        + cb * qb
                    )

                    reduced = reduce_mod_q13(
                        expr
                    )

                    A, B = split_linear_in_S(
                        reduced
                    )

                    # We prefer A=0:
                    # expression independent of S.
                    score = (
                        0 if A == 0 else 1,
                        complexity(reduced),
                    )

                    candidates.append(
                        (
                            score,
                            (ca, a, cb, names[j]),
                            reduced,
                            A,
                            B,
                        )
                    )

    candidates.sort(
        key=lambda x: x[0]
    )

    shown = 0

    for score, name, reduced, A, B in candidates:

        print(
            f"\ncombination:"
        )
        print(
            f"  {name[0]}*Q_{name[1]} "
            f"+ {name[2]}*Q_{name[3]}"
        )
        print(
            f"  reduced = {reduced}"
        )
        print(
            f"  A = {A}"
        )
        print(
            f"  B = {B}"
        )

        shown += 1

        if shown >= 30:
            break


# =============================================================================
# FOCUS ON Q13 AND Q15
# =============================================================================

def focus_q13_q15(
    quotients: Dict[Tuple[int, int], sp.Expr],
) -> None:

    print(
        "\n6. FOCUS: Q13 / Q15"
    )
    print("-" * 78)

    q13 = quotients[(1, 3)]
    q15 = quotients[(1, 5)]

    print(
        "Q13 ="
    )
    print(
        sp.factor(q13)
    )

    print(
        "\nQ15 ="
    )
    print(
        sp.factor(q15)
    )

    r15 = reduce_mod_q13(
        q15
    )

    print(
        "\nQ15 reduced modulo Q13 relation ="
    )
    print(
        sp.factor(r15)
    )

    A15, B15 = split_linear_in_S(
        r15
    )

    print(
        "\nQ15 = A15*S + B15"
    )
    print(
        f"A15 = {A15}"
    )
    print(
        f"B15 = {B15}"
    )

    # Try to solve for S.
    if A15 != 0:
        print(
            "\nFORMAL RESOLVENT:"
        )
        print(
            "S = (Q15 - B15)/A15"
        )

        print(
            "A15 ="
        )
        print(
            A15
        )

        print(
            "B15 ="
        )
        print(
            B15
        )


# =============================================================================
# Q13 RECONSTRUCTION TEST
# =============================================================================

def oracle_q13_linear_resolvent_test(
    quotients: Dict[Tuple[int, int], sp.Expr],
    targets: Sequence[Target],
) -> None:

    print(
        "\n7. ORACLE Q13 / Q15 RESOLUTION"
    )
    print("-" * 78)

    q13 = quotients[(1, 3)]
    q15 = quotients[(1, 5)]

    r15 = reduce_mod_q13(
        q15
    )

    A15, B15 = split_linear_in_S(
        r15
    )

    if A15 == 0:
        print(
            "Q15 contains no linear S coefficient after reduction."
        )
        return

    failures = 0

    for t in targets:

        q13_value = eval_ns(
            q13,
            t.n,
            t.s,
        )

        q15_value = eval_ns(
            q15,
            t.n,
            t.s,
        )

        a_value = int(
            A15.subs(
                {
                    N: t.n,
                    Z: q13_value,
                }
            )
        )

        b_value = int(
            B15.subs(
                {
                    N: t.n,
                    Z: q13_value,
                }
            )
        )

        recovered_s = (
            q15_value - b_value
        )

        if a_value == 0:
            failures += 1
            continue

        if recovered_s % a_value != 0:
            failures += 1
            continue

        recovered_s //= a_value

        if recovered_s != t.s:
            failures += 1

    print(
        f"oracle recovered S = "
        f"{len(targets)-failures}/"
        f"{len(targets)}"
    )


# =============================================================================
# FINAL
# =============================================================================

def main() -> None:

    global TARGETS_GLOBAL

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 98")
    print("PAPER QUOTIENT-ALGEBRA COLLAPSE")
    print("REMOVE UNIVERSAL (S+1) FACTOR")
    print("Q13 AS QUADRATIC GENERATOR")
    print("LINEAR REDUCTION MODULO Q13")
    print("N-ONLY CANCELLATION SEARCH")
    print("STRICT TARGET HOLDOUT")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # Prime population
    # -------------------------------------------------------------------------

    print(
        "\n1. PRIME POPULATION"
    )
    print("-" * 78)

    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_MIN,
        PRIME_MAX,
    )

    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # -------------------------------------------------------------------------
    # Targets
    # -------------------------------------------------------------------------

    TARGETS_GLOBAL = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = {len(TARGETS_GLOBAL)}"
    )

    for i, t in enumerate(
        TARGETS_GLOBAL[:20],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    print()
    print(
        f"training targets = {TRAIN_TARGETS}"
    )
    print(
        f"test targets     = {TEST_TARGETS}"
    )

    # -------------------------------------------------------------------------
    # Build quotients
    # -------------------------------------------------------------------------

    print(
        "\n2. PAPER QUOTIENT CONSTRUCTION"
    )
    print("-" * 78)

    quotients = build_quotients()

    for pair in DETECTORS:
        print(
            f"Q_{pair} = "
            f"{sp.factor(quotients[pair])}"
        )

    validate_quotients(
        TARGETS_GLOBAL,
        quotients,
    )

    q13 = quotients[(1, 3)]

    print(
        "\nQ13 ="
    )
    print(
        sp.factor(q13)
    )

    expected_q13 = sp.factor(
        6*N - S**2 + S
    )

    print(
        "\nQ13 expected ="
    )
    print(
        expected_q13
    )

    if sp.expand(q13 - expected_q13) != 0:
        raise ArithmeticError(
            "Q13 identity mismatch."
        )

    print(
        "Q13 identity = PASS"
    )

    # -------------------------------------------------------------------------
    # Reduction
    # -------------------------------------------------------------------------

    linear_forms = reduction_table(
        quotients,
        q13,
    )

    # -------------------------------------------------------------------------
    # Linear elimination
    # -------------------------------------------------------------------------

    search_linear_elimination(
        linear_forms
    )

    # -------------------------------------------------------------------------
    # Small combinations
    # -------------------------------------------------------------------------

    search_q_combinations(
        quotients,
        q13,
    )

    # -------------------------------------------------------------------------
    # Q13 / Q15
    # -------------------------------------------------------------------------

    focus_q13_q15(
        quotients
    )

    oracle_q13_linear_resolvent_test(
        quotients,
        TARGETS_GLOBAL,
    )

    # -------------------------------------------------------------------------
    # Final diagnostic
    # -------------------------------------------------------------------------

    print(
        "\n8. FINAL DIAGNOSTIC"
    )
    print("=" * 78)

    print(
        """
The experiment removes the universal (S+1) factor from the paper
detectors and uses

    Q13 = 6N - S^2 + S

as the base quadratic relation.

Every other Q_(k,l) can then be reduced to:

    Q_(k,l) = A(N,Q13) * S + B(N,Q13).

The critical outcomes are:

1. A(N,Q13) = 0
   ----------------
   The detector becomes a function of N and Q13 alone.

2. A(N,Q13) != 0
   ----------------
   The detector provides a linear resolvent:

       S = (Q_(k,l) - B) / A.

3. A small combination becomes independent of S
   ---------------------------------------------
   This would be particularly interesting because it could provide a
   relation among observable paper coefficients without the hidden S.

4. Several quotient relations collapse together
   ----------------------------------------------
   This would indicate that the paper family may be generated by a
   low-dimensional algebra over Z[N].

The decisive future objective remains:

    can Q13 or another quotient be evaluated from N alone?

But this experiment tells us whether the entire paper detector family
collapses to a much smaller algebraic system first.
"""
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 98 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

