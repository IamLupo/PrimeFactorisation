#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 97
MULTI-DETECTOR LINEAR RESOLVENT SEARCH

PAPER PRIME-DETECTOR FAMILY
---------------------------
f_(k,l)(p,q)
  = (1+q)^l p^k - (1+q)^k p^l
    + (1+p)^l q^k - (1+p)^k q^l

SEMIPRIME VARIABLES
-------------------
N = p*q
S = p+q

KNOWN:
    E1 = S*((N+1)^2 - S^2)

OBJECTIVE
---------
Search combinations of paper-derived coefficients for the simplest
algebraic relation that determines S.

Especially test whether:

    a(N) * S = b(N)

or

    a(N) * S + b(N) = 0

can be obtained from two or more known coefficient values.

This is different from the earlier experiments:

    - no million-candidate S sweep
    - no local character classifier
    - no CRT candidate collapse
    - no factor-pair test

The main output is an algebraic complexity ranking.

A successful result would identify a very low-complexity resolvent
which can then be targeted in the next experiment for N-only access.

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
from itertools import combinations
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

DETECTORS: Tuple[Tuple[int, int], ...] = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

# Limit coefficient multiplier search.
SMALL_COEFFS = (-3, -2, -1, 1, 2, 3)


# ============================================================================
# SYMBOLS
# ============================================================================

P, Q = sp.symbols("P Q")
N, S = sp.symbols("N S")


# ============================================================================
# TARGET
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

def sieve_primes(lo: int, hi: int) -> List[int]:
    flags = bytearray(b"\x01") * (hi + 1)
    flags[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(hi) + 1):
        if flags[p]:
            start = p * p
            flags[start : hi + 1 : p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x for x in range(lo, hi + 1)
        if flags[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:
    rng = random.Random(97097)

    out: List[Target] = []
    seen = set()

    while len(out) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        out.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return out


# ============================================================================
# PAPER DETECTOR
# ============================================================================

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


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduce(expr: sp.Expr) -> sp.Expr:
    reduced, remainder, mapping = sp.symmetrize(
        expr,
        [P, Q],
        formal=True,
    )

    if remainder != 0:
        raise ArithmeticError(
            f"Non-symmetric remainder: {remainder}"
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


# ============================================================================
# E1
# ============================================================================

E1 = sp.expand(
    S * ((N + 1) ** 2 - S**2)
)


# ============================================================================
# DETECTOR POLYNOMIALS
# ============================================================================

def build_detector_polynomials() -> Dict[Tuple[int, int], sp.Expr]:
    result = {}

    for pair in DETECTORS:
        result[pair] = sp.factor(
            symmetric_reduce(
                paper_detector(*pair)
            )
        )

    return result


# ============================================================================
# VALIDATION
# ============================================================================

def eval_direct(
    target: Target,
    pair: Tuple[int, int],
) -> int:
    k, ell = pair

    return int(
        (1 + target.q) ** ell * target.p**k
        - (1 + target.q) ** k * target.p**ell
        + (1 + target.p) ** ell * target.q**k
        - (1 + target.p) ** k * target.q**ell
    )


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


def validate_polynomials(
    targets: Sequence[Target],
    detector_polys: Dict[Tuple[int, int], sp.Expr],
) -> None:

    failures = 0

    for pair, expr in detector_polys.items():
        for t in targets:
            direct = eval_direct(t, pair)
            reduced = eval_ns(expr, t.n, t.s)

            if direct != reduced:
                failures += 1

    print(
        f"detector identity failures = {failures}"
    )

    if failures:
        raise ArithmeticError(
            "Symmetric reduction validation failed."
        )


# ============================================================================
# AFFINE REDUCTION IN S
# ============================================================================

def reduce_mod_candidate_relation(
    expr: sp.Expr,
) -> sp.Expr:
    """
    Polynomial cleanup over Z[N,S].
    """
    return sp.factor(
        sp.cancel(
            sp.expand(expr)
        )
    )


# ============================================================================
# RESULTANT ELIMINATION
# ============================================================================

def eliminate_s(
    f: sp.Expr,
    g: sp.Expr,
) -> sp.Expr:
    """
    Eliminate S between f=0 and g=0.

    This is an algebraic diagnostic only.
    """
    return sp.factor(
        sp.resultant(
            f,
            g,
            S,
        )
    )


# ============================================================================
# SUBRESULTANT SEARCH
# ============================================================================

def subresultant_linear_candidate(
    f: sp.Expr,
    g: sp.Expr,
) -> Tuple[sp.Expr | None, int | None]:
    """
    Search the subresultant sequence for a polynomial of degree 1 in S.

    Returns:
        (linear polynomial, degree)
    """
    pf = sp.Poly(f, S)
    pg = sp.Poly(g, S)

    try:
        seq = sp.subresultants(
            pf,
            pg,
            S,
        )
    except Exception:
        return None, None

    best = None

    for sub in seq:
        expr = sp.factor(
            sub.as_expr()
        )

        degree = sp.Poly(
            expr,
            S,
        ).degree()

        if degree == 1:
            best = expr

    if best is None:
        return None, None

    return best, 1


# ============================================================================
# NORMALIZE LINEAR POLYNOMIAL
# ============================================================================

def normalize_linear(
    expr: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr] | None:
    """
    Normalize:

        a(N) S + b(N) = 0

    into a(N), b(N).
    """
    poly = sp.Poly(
        expr,
        S,
    )

    if poly.degree() != 1:
        return None

    coeffs = poly.all_coeffs()

    a = sp.factor(coeffs[0])
    b = sp.factor(coeffs[1])

    # Remove common polynomial factor if possible.
    common = sp.gcd(
        sp.Poly(a, N),
        sp.Poly(b, N),
    ).as_expr()

    if common not in (0, 1, -1):
        a = sp.factor(a / common)
        b = sp.factor(b / common)

    return (
        sp.factor(a),
        sp.factor(b),
    )


# ============================================================================
# COMPLEXITY SCORE
# ============================================================================

def expression_complexity(
    expr: sp.Expr,
) -> Tuple[int, int, int, int]:
    """
    Lower is better.

    Returns:
        total degree
        S degree
        N degree of coefficients
        expanded term count
    """
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

    coeff_degrees = []

    for c in sp.Poly(
        sp.expand(expr),
        S,
    ).all_coeffs():
        coeff_degrees.append(
            sp.Poly(
                c,
                N,
            ).degree()
        )

    max_n_degree = max(
        coeff_degrees,
        default=0,
    )

    terms = len(
        sp.Poly(
            sp.expand(expr),
            N,
            S,
        ).terms()
    )

    return (
        total_degree,
        s_degree,
        max_n_degree,
        terms,
    )


# ============================================================================
# VERIFY LINEAR RESOLVENT NUMERICALLY
# ============================================================================

def verify_linear_resolvent(
    expr: sp.Expr,
    targets: Sequence[Target],
) -> Tuple[int, int]:

    normalized = normalize_linear(
        expr
    )

    if normalized is None:
        return 0, len(targets)

    a, b = normalized

    ok = 0
    fail = 0

    for t in targets:
        lhs = (
            int(
                a.subs(N, t.n)
            )
            * t.s
            + int(
                b.subs(N, t.n)
            )
        )

        if lhs == 0:
            ok += 1
        else:
            fail += 1

    return ok, fail


# ============================================================================
# SPECIAL TWO-EQUATION RESOLVENT
# ============================================================================

def special_resolvent_check(
    targets: Sequence[Target],
    detector_expr: sp.Expr,
    pair: Tuple[int, int],
) -> None:

    print(
        f"\nSPECIAL CHECK: E1 + f_{pair}"
    )

    f = detector_expr

    # IMPORTANT:
    #
    # The actual detector values A and B are not set to zero.
    #
    # For every target:
    #
    #   E1(N,S) - A = 0
    #   f(N,S)  - B = 0
    #
    # We therefore construct the symbolic resultant with symbolic
    # constants A and B.

    A, B = sp.symbols(
        "A B"
    )

    F = sp.expand(
        E1 - A
    )

    G = sp.expand(
        f - B
    )

    print(
        "searching subresultants..."
    )

    linear_found, _ = (
        subresultant_linear_candidate(
            F,
            G,
        )
    )

    if linear_found is None:
        print(
            "No linear subresultant found."
        )
        return

    print(
        "linear subresultant:"
    )
    print(
        sp.factor(
            linear_found
        )
    )

    normalized = normalize_linear(
        linear_found
    )

    if normalized is not None:
        a, b = normalized

        print(
            "\nnormalized:"
        )

        print(
            f"a(N,A,B) = {a}"
        )

        print(
            f"b(N,A,B) = {b}"
        )

    # Numerical verification.
    failures = 0

    for t in targets:

        A_value = eval_ns(
            E1,
            t.n,
            t.s,
        )

        B_value = eval_ns(
            f,
            t.n,
            t.s,
        )

        relation = sp.expand(
            linear_found
            .subs(
                {
                    A: A_value,
                    B: B_value,
                }
            )
        )

        relation_value = int(
            relation.subs(
                {
                    N: t.n,
                    S: t.s,
                }
            )
        )

        if relation_value != 0:
            failures += 1

    print(
        f"held-out identity failures = "
        f"{failures}/{len(targets)}"
    )


# ============================================================================
# SEARCH PAIRS
# ============================================================================

def search_detector_pairs(
    detector_polys: Dict[Tuple[int, int], sp.Expr],
    targets: Sequence[Target],
) -> None:

    print(
        "\n3. TWO-DETECTOR SUBRESULTANT SEARCH"
    )
    print("-" * 78)

    rows = []

    for a, b in combinations(
        DETECTORS,
        2,
    ):

        f = detector_polys[a]
        g = detector_polys[b]

        print(
            f"\nPAIR f_{a} / f_{b}"
        )

        sub, _ = (
            subresultant_linear_candidate(
                f,
                g,
            )
        )

        if sub is not None:

            print(
                "  linear subresultant:"
            )
            print(
                "   ",
                sp.factor(sub)
            )

            norm = normalize_linear(
                sub
            )

            if norm is not None:

                aa, bb = norm

                print(
                    f"  a(N) = {aa}"
                )
                print(
                    f"  b(N) = {bb}"
                )

                ok, fail = (
                    verify_linear_resolvent(
                        sub,
                        targets,
                    )
                )

                print(
                    f"  direct target verification = "
                    f"{ok}/{len(targets)}"
                )

                rows.append(
                    {
                        "pair": (a, b),
                        "expr": sub,
                        "complexity": expression_complexity(sub),
                        "verified": ok,
                    }
                )

        else:
            print(
                "  no linear subresultant"
            )

        # Also compute polynomial gcd.
        gcd_expr = sp.factor(
            sp.gcd(
                sp.Poly(f, S),
                sp.Poly(g, S),
            ).as_expr()
        )

        print(
            f"  gcd in S = {gcd_expr}"
        )

    print(
        "\n4. BEST DIRECT DETECTOR-PAIR RESOLVENTS"
    )
    print("-" * 78)

    rows.sort(
        key=lambda x: (
            x["verified"] != len(targets),
            x["complexity"],
        )
    )

    for row in rows:
        print(
            f"pair={row['pair']} "
            f"complexity={row['complexity']} "
            f"verified={row['verified']}/{len(targets)}"
        )
        print(
            " ",
            sp.factor(
                row["expr"]
            )
        )


# ============================================================================
# E1 + DETECTOR SEARCH
# ============================================================================

def run_e1_detector_search(
    detector_polys: Dict[Tuple[int, int], sp.Expr],
    targets: Sequence[Target],
) -> None:

    print(
        "\n5. E1 + PAPER DETECTOR RESOLVENTS"
    )
    print("-" * 78)

    for pair in DETECTORS:
        special_resolvent_check(
            targets,
            detector_polys[pair],
            pair,
        )


# ============================================================================
# INTEGER LINEAR COMBINATION SEARCH
# ============================================================================

def search_simple_linear_combinations(
    detector_polys: Dict[Tuple[int, int], sp.Expr],
    targets: Sequence[Target],
) -> None:

    print(
        "\n6. SMALL-INTEGER LINEAR COMBINATION SEARCH"
    )
    print("-" * 78)

    basis_names = [
        "E1"
    ] + [
        f"f_{pair}"
        for pair in DETECTORS
    ]

    basis_exprs = [
        E1
    ] + [
        detector_polys[pair]
        for pair in DETECTORS
    ]

    candidates = []

    # Pair combinations only, to keep this experiment interpretable.
    for i in range(
        len(basis_exprs)
    ):
        for j in range(
            i + 1,
            len(basis_exprs)
        ):

            for c1 in SMALL_COEFFS:
                for c2 in SMALL_COEFFS:

                    expr = sp.factor(
                        c1 * basis_exprs[i]
                        + c2 * basis_exprs[j]
                    )

                    # Ignore zero.
                    if expr == 0:
                        continue

                    degree = sp.Poly(
                        expr,
                        S,
                    ).degree()

                    # We are looking for low S degree.
                    if degree > 3:
                        continue

                    complexity = (
                        expression_complexity(
                            expr
                        )
                    )

                    candidates.append(
                        (
                            complexity,
                            f"{c1}*{basis_names[i]} "
                            f"+ {c2}*{basis_names[j]}",
                            expr,
                        )
                    )

    candidates.sort(
        key=lambda x: x[0]
    )

    print(
        f"candidate combinations = "
        f"{len(candidates)}"
    )

    shown = 0

    for complexity, name, expr in candidates:

        print(
            f"\n{name}"
        )
        print(
            f"complexity={complexity}"
        )
        print(
            sp.factor(expr)
        )

        shown += 1

        if shown >= 30:
            break


# ============================================================================
# TEST THE KNOWN f13 - E1 RELATION
# ============================================================================

def explicit_f13_difference_check(
    detector_polys: Dict[Tuple[int, int], sp.Expr],
    targets: Sequence[Target],
) -> None:

    print(
        "\n7. EXPLICIT f_(1,3) - E1 IDENTITY"
    )
    print("-" * 78)

    f13 = detector_polys[(1, 3)]

    difference = sp.factor(
        E1 - f13
    )

    print(
        "E1 - f_(1,3) ="
    )
    print(
        difference
    )

    print(
        "\nnormalized relation:"
    )

    relation = sp.factor(
        (
            E1
            - f13
            - N * (
                S * (N - 4)
                - 6
            )
        )
    )

    print(
        relation
    )

    if relation == 0:
        print(
            "SYMBOLIC CHECK = PASS"
        )
    else:
        print(
            "SYMBOLIC CHECK = FAIL"
        )
        raise ArithmeticError(
            "f13-E1 relation failed."
        )

    failures = 0

    for t in targets:

        value = eval_ns(
            difference,
            t.n,
            t.s,
        )

        expected = (
            t.n
            * (
                t.s * (
                    t.n - 4
                )
                - 6
            )
        )

        if value != expected:
            failures += 1

    print(
        f"numeric validation = "
        f"{len(targets)-failures}/{len(targets)}"
    )


# ============================================================================
# FINAL
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 97")
    print("MULTI-DETECTOR LINEAR RESOLVENT SEARCH")
    print("E1 + PAPER PRIME-DETECTOR FAMILY")
    print("SUBRESULTANTS / RESULTANTS / LOW-DEGREE COMBINATIONS")
    print("STRICT TARGET HOLDOUT")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

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
        f"prime population = "
        f"{len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = "
        f"{len(targets)}"
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

    print()
    print(
        f"training targets = "
        f"{TRAIN_TARGETS}"
    )
    print(
        f"test targets     = "
        f"{TEST_TARGETS}"
    )

    print(
        "\n2. PAPER DETECTOR CONSTRUCTION"
    )
    print("-" * 78)

    detector_polys = (
        build_detector_polynomials()
    )

    for pair in DETECTORS:
        print(
            f"f_{pair} = "
            f"{sp.factor(detector_polys[pair])}"
        )

    validate_polynomials(
        targets,
        detector_polys,
    )

    print(
        "\nidentity validation = PASS"
    )

    # ----------------------------------------------------------------------
    # SEARCH
    # ----------------------------------------------------------------------

    search_detector_pairs(
        detector_polys,
        targets,
    )

    run_e1_detector_search(
        detector_polys,
        targets,
    )

    search_simple_linear_combinations(
        detector_polys,
        targets,
    )

    explicit_f13_difference_check(
        detector_polys,
        targets,
    )

    # ----------------------------------------------------------------------
    # Final diagnostic
    # ----------------------------------------------------------------------

    print(
        "\n8. FINAL DIAGNOSTIC"
    )
    print("=" * 78)

    print(
        """
The objective is to find the lowest-complexity equation involving
N and S that can be constructed from the paper-derived coefficients.

The most valuable outcome is:

    A(N) * S + B(N) = 0

because then:

    S = -B(N) / A(N)

and the semiprime factorization follows from:

    x^2 - S*x + N = 0.

Three outcomes are distinguished:

1. DIRECT LINEAR RESOLVENT
   A combination of actual paper coefficients produces a linear
   polynomial in S.

2. ALGEBRAIC RESOLVENT AFTER TWO KNOWN COEFFICIENT VALUES
   The subresultant of
       E1(N,S)-A
       f(N,S)-B
   contains a linear polynomial in S.
   This proves that two paper coefficients are sufficient to recover S
   algebraically once their values are accessible.

3. NO LOW-COMPLEXITY RESOLVENT
   The family remains high-degree and does not simplify.

IMPORTANT:
A relation involving oracle coefficient values is NOT an N-only
factorization algorithm.

The next experiment would only become computationally decisive if one
of these coefficient combinations can itself be evaluated from N
without first knowing p or q.
"""
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 97 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

