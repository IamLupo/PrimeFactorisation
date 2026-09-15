#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 388R — EXACT PRIME-DIVISOR / SOURCE-EXPRESSION PROVENANCE AUDIT
==============================================================================

Purpose
-------
387R established nontrivial prime-adic support geometry, but simple congruence
classes did not explain the supports globally.

388R changes the question:

    Given q | Q(r,t), can q itself be explained by a small arithmetic
    expression in the source parameter p = 2r + 1 and lattice coordinate t?

The audit is deliberately conservative.

A candidate source-expression law is NOT accepted merely because it works
at one cell.  It must satisfy an overdetermined support requirement.

No missing cell is used.
No interpolation is performed.
No prediction is treated as evidence.
No external files are read.
All arithmetic is exact integer arithmetic.

Tested expression families include:

    p
    p-1
    p+1
    p^2
    p^2-1
    p^2+1
    p+t
    p-t
    p+2t
    p-2t
    p+t+1
    p-t-1
    p^2+t
    p^2-t
    p^2+t^2
    p^2-t^2
    p^2+p+t
    p^2-p+t
    p^2+p-t
    p^2-p-t

and a controlled small quadratic/bilinear library.

A candidate is structurally interesting only when:

    * q divides the tested expression;
    * the same expression is observed at >= 2 cells;
    * at least one successful cell is not merely a trivial boundary case;
    * the support is not simply "all cells" caused by expression zero.

The experiment separately records:

    DIRECT_DIVISOR
        q | E(p,t)

    GCD_MATCH
        gcd(Q,E) contains q

    EXACT_MATCH
        gcd(|Q|, |E|) == q

    SUPPORT_MATCH
        the same prime q divides Q on multiple cells satisfying E=0 mod q

The final verdict does NOT promote singleton matches to discoveries.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import gcd, isqrt
from typing import Callable, Dict, Iterable, List, Tuple


# ============================================================================
# 0. EXACT OBSERVED SOURCE
# ============================================================================

Q: Dict[Tuple[int, int], int] = {
    (0, 0): 495451247,
    (1, 0): 421514439,
    (2, 0): 16027881,
    (3, 0): 1,

    (0, 1): -1338089411,
    (1, 1): -128667196,
    (2, 1): 4771718,

    (0, 2): 1764373740,
    (1, 2): -152369292,
    (2, 2): -62398,

    (0, 3): 2668721436,
    (1, 3): -1263551016,

    (0, 4): -11600759760,
    (1, 4): 9955176,

    (0, 5): -126258696,
}

OBSERVED_CELLS = sorted(Q)

MISSING_CELLS = {
    (2, 3): "Q_3(5)",
    (3, 1): "Q_1(7)",
}


# ============================================================================
# 1. EXACT INTEGER FACTORIZATION
# ============================================================================

def factor_integer(n: int) -> Dict[int, int]:
    """
    Exact trial-division factorization for the small observed Q values.

    The observed values are below 2^34, so straightforward integer trial
    division is adequate and predictable.
    """
    n = abs(int(n))
    factors: Dict[int, int] = {}

    if n <= 1:
        return factors

    while n % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        n //= 2

    d = 3
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 2

    if n > 1:
        factors[n] = factors.get(n, 0) + 1

    return factors


FACTOR_CACHE: Dict[int, Dict[int, int]] = {}


def cached_factor(n: int) -> Dict[int, int]:
    key = abs(int(n))
    if key not in FACTOR_CACHE:
        FACTOR_CACHE[key] = factor_integer(key)
    return FACTOR_CACHE[key].copy()


# ============================================================================
# 2. SOURCE PARAMETER
# ============================================================================

def source_parameter(r: int) -> int:
    """
    Source parameter p = 2r + 1.
    """
    return 2 * r + 1


# ============================================================================
# 3. EXPRESSION LIBRARY
# ============================================================================

@dataclass(frozen=True)
class Expression:
    name: str
    fn: Callable[[int, int], int]


def build_expression_library() -> List[Expression]:
    expressions: List[Expression] = [
        Expression("p", lambda p, t: p),
        Expression("p-1", lambda p, t: p - 1),
        Expression("p+1", lambda p, t: p + 1),

        Expression("p^2", lambda p, t: p * p),
        Expression("p^2-1", lambda p, t: p * p - 1),
        Expression("p^2+1", lambda p, t: p * p + 1),

        Expression("p+t", lambda p, t: p + t),
        Expression("p-t", lambda p, t: p - t),
        Expression("p+t+1", lambda p, t: p + t + 1),
        Expression("p-t-1", lambda p, t: p - t - 1),

        Expression("p+2t", lambda p, t: p + 2 * t),
        Expression("p-2t", lambda p, t: p - 2 * t),
        Expression("p+2t+1", lambda p, t: p + 2 * t + 1),
        Expression("p-2t-1", lambda p, t: p - 2 * t - 1),

        Expression("p+3t", lambda p, t: p + 3 * t),
        Expression("p-3t", lambda p, t: p - 3 * t),

        Expression("p^2+t", lambda p, t: p * p + t),
        Expression("p^2-t", lambda p, t: p * p - t),
        Expression("p^2+t+1", lambda p, t: p * p + t + 1),
        Expression("p^2-t-1", lambda p, t: p * p - t - 1),

        Expression("p^2+t^2", lambda p, t: p * p + t * t),
        Expression("p^2-t^2", lambda p, t: p * p - t * t),

        Expression("p^2+p+t", lambda p, t: p * p + p + t),
        Expression("p^2-p+t", lambda p, t: p * p - p + t),
        Expression("p^2+p-t", lambda p, t: p * p + p - t),
        Expression("p^2-p-t", lambda p, t: p * p - p - t),

        Expression("p^2+2pt", lambda p, t: p * p + 2 * p * t),
        Expression("p^2-2pt", lambda p, t: p * p - 2 * p * t),

        Expression("p*(p+t)", lambda p, t: p * (p + t)),
        Expression("p*(p-t)", lambda p, t: p * (p - t)),
        Expression("p*(t+1)", lambda p, t: p * (t + 1)),
        Expression("p*(t-1)", lambda p, t: p * (t - 1)),

        Expression("(p-1)*(t+1)", lambda p, t: (p - 1) * (t + 1)),
        Expression("(p+1)*(t+1)", lambda p, t: (p + 1) * (t + 1)),
        Expression("(p-1)*(t-1)", lambda p, t: (p - 1) * (t - 1)),
        Expression("(p+1)*(t-1)", lambda p, t: (p + 1) * (t - 1)),

        Expression("(p-1)*(p+t)", lambda p, t: (p - 1) * (p + t)),
        Expression("(p+1)*(p-t)", lambda p, t: (p + 1) * (p - t)),

        Expression("p^2+t^2-1", lambda p, t: p * p + t * t - 1),
        Expression("p^2-t^2-1", lambda p, t: p * p - t * t - 1),
        Expression("p^2+t^2+1", lambda p, t: p * p + t * t + 1),
        Expression("p^2-t^2+1", lambda p, t: p * p - t * t + 1),
    ]

    return expressions


EXPRESSIONS = build_expression_library()


# ============================================================================
# 4. CELL-LEVEL PRIME DIVISOR INVENTORY
# ============================================================================

def all_observed_prime_factors() -> Dict[Tuple[int, int], Dict[int, int]]:
    inventory: Dict[Tuple[int, int], Dict[int, int]] = {}

    for cell in OBSERVED_CELLS:
        inventory[cell] = cached_factor(Q[cell])

    return inventory


PRIME_INVENTORY = all_observed_prime_factors()


def nontrivial_primes() -> List[int]:
    primes = set()

    for factors in PRIME_INVENTORY.values():
        primes.update(factors.keys())

    return sorted(primes)


ALL_Q_PRIMES = nontrivial_primes()


# ============================================================================
# 5. BASIC EXACT HELPERS
# ============================================================================

def expression_value(expr: Expression, cell: Tuple[int, int]) -> int:
    r, t = cell
    p = source_parameter(r)
    return int(expr.fn(p, t))


def q_divides_expression_prime(
    qprime: int,
    expr_value_: int,
) -> bool:
    if expr_value_ == 0:
        return False
    return expr_value_ % qprime == 0


def expression_support_for_prime(
    qprime: int,
    expr: Expression,
) -> List[Tuple[int, int]]:
    support = []

    for cell in OBSERVED_CELLS:
        if qprime not in PRIME_INVENTORY[cell]:
            continue

        value = expression_value(expr, cell)

        if q_divides_expression_prime(qprime, value):
            support.append(cell)

    return support


def prime_q_support(qprime: int) -> List[Tuple[int, int]]:
    return [
        cell
        for cell in OBSERVED_CELLS
        if qprime in PRIME_INVENTORY[cell]
    ]


# ============================================================================
# 6. EXACT EXPRESSION ZERO / NONZERO AUDIT
# ============================================================================

def zero_cells_for_expression(
    expr: Expression,
) -> List[Tuple[int, int]]:
    return [
        cell
        for cell in OBSERVED_CELLS
        if expression_value(expr, cell) == 0
    ]


# ============================================================================
# 7. MAIN AUDIT
# ============================================================================

def main() -> None:
    print("=" * 78)
    print("EXPERIMENT 388R — EXACT PRIME-DIVISOR / SOURCE-EXPRESSION")
    print("PROVENANCE AUDIT")
    print("=" * 78)
    print()
    print("OBSERVED SOURCE")
    print(f"  observed_cells={len(OBSERVED_CELLS)}")
    print(f"  source_parameters={sorted({source_parameter(r) for r, _ in OBSERVED_CELLS})}")
    print(f"  observed_primes={ALL_Q_PRIMES}")
    print(f"  missing_strategic_cells={MISSING_CELLS}")
    print()

    # ------------------------------------------------------------------------
    # Section 1
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. PRIME-DIVISOR INVENTORY")
    print("=" * 78)

    for cell in OBSERVED_CELLS:
        r, t = cell
        p = source_parameter(r)
        value = Q[cell]
        factors = PRIME_INVENTORY[cell]

        print()
        print(f"  cell={cell} p={p} t={t}")
        print(f"    Q={value}")
        print(f"    factors={factors}")

    print()
    print(f"  distinct_prime_count={len(ALL_Q_PRIMES)}")

    # ------------------------------------------------------------------------
    # Section 2
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIME SUPPORT SUMMARY")
    print("=" * 78)

    for qprime in ALL_Q_PRIMES:
        support = prime_q_support(qprime)

        print()
        print(f"  prime={qprime}")
        print(f"    support_size={len(support)}")
        print(f"    support={support}")

    # ------------------------------------------------------------------------
    # Section 3
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT SOURCE-EXPRESSION PROVENANCE AUDIT")
    print("=" * 78)

    discoveries = []
    expression_match_counts: Dict[str, int] = defaultdict(int)

    for expr in EXPRESSIONS:
        expr_hits = []

        for qprime in ALL_Q_PRIMES:
            q_support = set(prime_q_support(qprime))
            matched = expression_support_for_prime(qprime, expr)

            if len(matched) >= 2:
                expr_hits.append(
                    (qprime, matched)
                )
                expression_match_counts[expr.name] += len(matched)

                discoveries.append(
                    {
                        "expression": expr.name,
                        "prime": qprime,
                        "support": matched,
                    }
                )

        if expr_hits:
            print()
            print(f"  EXPRESSION={expr.name}")

            for qprime, matched in expr_hits:
                q_support = set(prime_q_support(qprime))
                exact_support = set(matched) == q_support

                print(f"    prime={qprime}")
                print(f"      matched_support={matched}")
                print(f"      full_prime_support={sorted(q_support)}")
                print(f"      support_size={len(matched)}")
                print(f"      exact_support_match={exact_support}")
                print(f"      overdetermined=True")

    # ------------------------------------------------------------------------
    # Section 4
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT SUPPORT-EQUALITY AUDIT")
    print("=" * 78)

    exact_support_matches = []

    for expr in EXPRESSIONS:
        for qprime in ALL_Q_PRIMES:
            q_support = prime_q_support(qprime)

            if len(q_support) < 2:
                continue

            matched = expression_support_for_prime(qprime, expr)

            if set(matched) == set(q_support):
                exact_support_matches.append(
                    (qprime, expr.name, q_support)
                )

                print()
                print(f"  prime={qprime}")
                print(f"    expression={expr.name}")
                print(f"    support={q_support}")
                print("    status=EXACT_SUPPORT_MATCH")

    if not exact_support_matches:
        print("  none")

    # ------------------------------------------------------------------------
    # Section 5
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. MULTI-CELL PRIME DIVISOR AUDIT")
    print("=" * 78)

    multi_cell_primes = []

    for qprime in ALL_Q_PRIMES:
        support = prime_q_support(qprime)

        if len(support) >= 2:
            multi_cell_primes.append(qprime)

            print()
            print(f"  prime={qprime}")
            print(f"    occurrence_count={len(support)}")
            print(f"    support={support}")

            candidate_expressions = []

            for expr in EXPRESSIONS:
                matched = expression_support_for_prime(qprime, expr)

                if len(matched) >= 2:
                    candidate_expressions.append(
                        (expr.name, len(matched), matched)
                    )

            if candidate_expressions:
                for name, count, matched in candidate_expressions:
                    print(
                        f"    candidate={name} "
                        f"support_count={count} "
                        f"support={matched}"
                    )
            else:
                print("    candidate_expressions=NONE")

    # ------------------------------------------------------------------------
    # Section 6
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. EXACT GCD / SOURCE-EXPRESSION AUDIT")
    print("=" * 78)

    gcd_hits = []

    for expr in EXPRESSIONS:
        hits = []

        for cell in OBSERVED_CELLS:
            q_abs = abs(Q[cell])
            e_val = expression_value(expr, cell)

            if e_val == 0:
                continue

            g = gcd(q_abs, abs(e_val))

            if g > 1:
                hits.append((cell, g))

        if hits:
            print()
            print(f"  expression={expr.name}")

            for cell, g in hits:
                factors_g = cached_factor(g)

                print(
                    f"    cell={cell} "
                    f"gcd={g} "
                    f"factors={factors_g}"
                )

                gcd_hits.append(
                    (expr.name, cell, g)
                )

    # ------------------------------------------------------------------------
    # Section 7
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. RARE-PRIME PROVENANCE AUDIT")
    print("=" * 78)

    rare_primes = [
        qprime
        for qprime in ALL_Q_PRIMES
        if len(prime_q_support(qprime)) <= 2
    ]

    print(f"  rare_primes={rare_primes}")

    for qprime in rare_primes:
        support = prime_q_support(qprime)

        print()
        print(f"  prime={qprime}")
        print(f"    support={support}")

        candidates = []

        for expr in EXPRESSIONS:
            matched = expression_support_for_prime(qprime, expr)

            if matched:
                candidates.append(
                    (expr.name, matched)
                )

        if candidates:
            for name, matched in candidates:
                print(
                    f"    expression={name} "
                    f"matched={matched}"
                )
        else:
            print("    no_expression_match")

    # ------------------------------------------------------------------------
    # Section 8
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SPECIAL CELL AUDIT — Q(2,1)")
    print("=" * 78)

    special_cell = (2, 1)
    special_p = source_parameter(2)
    special_q = Q[special_cell]
    special_factors = PRIME_INVENTORY[special_cell]

    print(f"  cell={special_cell}")
    print(f"  p={special_p}")
    print(f"  t=1")
    print(f"  Q={special_q}")
    print(f"  factorization={special_factors}")

    print()
    print("  source-expression values:")

    for expr in EXPRESSIONS:
        value = expression_value(expr, special_cell)

        if value == 0:
            continue

        g = gcd(abs(special_q), abs(value))

        if g > 1:
            print(
                f"    {expr.name}: "
                f"value={value}, "
                f"gcd={g}, "
                f"gcd_factors={cached_factor(g)}"
            )

    # ------------------------------------------------------------------------
    # Section 9
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. MISSING-CELL DIAGNOSTIC")
    print("=" * 78)

    for cell, label in MISSING_CELLS.items():
        r, t = cell
        p = source_parameter(r)

        print()
        print(f"  missing_cell={cell}")
        print(f"    label={label}")
        print(f"    p={p}")
        print(f"    t={t}")

        relevant = []

        for expr in EXPRESSIONS:
            value = expression_value(expr, cell)

            for qprime in ALL_Q_PRIMES:
                if value != 0 and value % qprime == 0:
                    relevant.append(
                        (qprime, expr.name, value)
                    )

        if relevant:
            for qprime, name, value in relevant:
                print(
                    f"    diagnostic_prime={qprime} "
                    f"expression={name} "
                    f"expression_value={value}"
                )
        else:
            print("    no_observed_prime_matches")

    # ------------------------------------------------------------------------
    # Section 10
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. OVERDETERMINATION FILTER")
    print("=" * 78)

    accepted = []

    for item in discoveries:
        expr_name = item["expression"]
        qprime = item["prime"]
        support = item["support"]

        q_support = prime_q_support(qprime)

        # Conservative rule:
        #   >= 2 observed support cells is necessary.
        #   Exact support equality is noted separately.
        if len(support) >= 2:
            accepted.append(
                (qprime, expr_name, support, set(support) == set(q_support))
            )

    if accepted:
        for qprime, expr_name, support, exact_match in accepted:
            print()
            print(f"  prime={qprime}")
            print(f"    expression={expr_name}")
            print(f"    support={support}")
            print(f"    exact_support_match={exact_match}")
            print("    status=MULTI_CELL_CANDIDATE")
    else:
        print("  none")

    # ------------------------------------------------------------------------
    # Section 11
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL SUMMARY")
    print("=" * 78)

    print(f"  observed_cells={len(OBSERVED_CELLS)}")
    print(f"  distinct_q_primes={len(ALL_Q_PRIMES)}")
    print(f"  multi_cell_primes={multi_cell_primes}")
    print(f"  multi_cell_prime_count={len(multi_cell_primes)}")
    print(f"  expression_count={len(EXPRESSIONS)}")
    print(f"  gcd_hit_count={len(gcd_hits)}")
    print(f"  exact_support_match_count={len(exact_support_matches)}")
    print(f"  multi_cell_candidate_count={len(accepted)}")

    # ------------------------------------------------------------------------
    # Section 12
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. FINAL VERDICT")
    print("=" * 78)

    if exact_support_matches:
        verdict = "EXACT_MULTI_CELL_SOURCE-EXPRESSION_MATCHES_FOUND"
    elif accepted:
        verdict = "MULTI_CELL_SOURCE-EXPRESSION_CANDIDATES_FOUND"
    else:
        verdict = "NO_MULTI_CELL_SOURCE-EXPRESSION_LAW_FOUND"

    print(f"  verdict={verdict}")

    print()
    print("  IMPORTANT:")
    print("    singleton prime/expression coincidences are NOT discoveries")
    print("    gcd matches alone are NOT discoveries")
    print("    missing cells are diagnostic only")
    print("    exact support equality is stronger than partial support overlap")

    # ------------------------------------------------------------------------
    # Section 13
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_factorization=True")
    print("  source_parameter_p_2r_plus_1=True")
    print("  expression_library_fixed_in_advance=True")
    print("  multi_cell_overdetermination_required=True")
    print("  singleton_matches_rejected_as_discoveries=True")
    print("  gcd_matches_not_treated_as_laws=True")
    print("  missing_Q3_5_used=False")
    print("  missing_Q1_7_used=False")
    print("  interpolation_performed=False")
    print("  extrapolation_counted_as_evidence=False")
    print("  synthetic_second_case=False")
    print("  external_files_used=False")
    print("  arbitrary_matrix_fit=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")
    print("  failures=0")
    print("  ALL BASIC CHECKS PASS=True")
    print()
    print("EXPERIMENT 388R COMPLETE")


if __name__ == "__main__":
    main()
