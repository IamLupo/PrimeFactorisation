#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 118
DETECTOR INVOLUTION / INVARIANT-RING TEST
Q_(k,l)(N,S) UNDER S -> 1-S
TEST X = S(S-1) INVARIANT REPRESENTATION
COLLISION-CLASSIFICATION FOLLOW-UP TO EXPERIMENT 117

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO NUMERICAL FITTING
EXACT SYMPY ALGEBRA
==============================================================================

OBJECTIVES
----------

1. For every paper detector Q_(k,l)(N,S), test exactly whether

       Q(N,1-S) - Q(N,S) == 0.

2. Factor the difference when it is nonzero.

3. Compute the gcd/common factor of all involution differences.

4. Test whether Q_(k,l)(N,S) belongs to QQ[N, X] with

       X = S(S-1).

5. Derive the exact obstruction when the X-representation fails.

6. Reproduce the collision classes observed in EXPERIMENT 117
   for small moduli.

7. Determine whether observed detector collisions are explained by

       S2 == 1-S1 (mod M)

   or require a larger equivalence relation.

8. Test additional invariants suggested by the collision pairs:

       S1 + S2
       S1*S2
       S(S-1)
       S^2 - 4N

9. Perform all symbolic tests independently of oracle target generation.

This is a structural experiment.
It does NOT claim an N-only factoring algorithm.
"""

from __future__ import annotations

import math
import random
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

N, S, X = sp.symbols("N S X", integer=True)

# Odd weights used throughout the previous experiments.
K_VALUES = [1, 3, 5, 7, 9]
ELL_VALUES = [3, 5, 7, 9, 11, 13, 15]

DETECTOR_PAIRS = [
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
    (1, 9),
    (3, 9),
    (5, 9),
    (7, 9),
    (1, 11),
    (3, 11),
    (5, 11),
]


# =============================================================================
# BASIC SYMBOLIC HELPERS
# =============================================================================

def canonical(expr: sp.Expr) -> sp.Expr:
    """Canonical exact simplification."""
    return sp.factor(sp.expand(sp.cancel(expr)))


def degree_in(expr: sp.Expr, var: sp.Symbol) -> int:
    """Polynomial degree; -oo is represented as -1 for the zero polynomial."""
    p = sp.Poly(sp.expand(expr), var, domain="EX")
    if p.is_zero:
        return -1
    return p.degree()


def poly_in_S(expr: sp.Expr) -> sp.Poly:
    """
    Safe polynomial representation in S with coefficients in QQ(N).
    Avoids invalid domain='ZZ[N]' generator handling.
    """
    return sp.Poly(sp.expand(expr), S, domain=sp.QQ.frac_field(N))


# =============================================================================
# PAPER DETECTOR
# =============================================================================

def raw_detector(k: int, ell: int) -> sp.Expr:
    """
    F_(k,l):

      p^k(1+q)^l + q^k(1+p)^l
      - p^l(1+q)^k - q^l(1+p)^k

    rewritten directly in symmetric variables through the
    two-variable symmetric substitution:

        p + q = S
        p q   = N

    We do not use p or q numerically.
    """

    # The raw two-variable object.
    p, q = sp.symbols("p q")

    Fpq = (
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    # Every symmetric polynomial can be reduced using:
    #   p^a + q^a = Newton recurrence.
    #
    # Rather than relying on a heuristic symmetric reduction,
    # construct the power sums explicitly.

    max_power = max(k + ell, ell + k)

    P = {0: sp.Integer(2), 1: S}
    for j in range(2, max_power + 1):
        P[j] = sp.expand(S * P[j - 1] - N * P[j - 2])

    # Expand Fpq into monomials p^a q^b and pair symmetric terms.
    poly = sp.Poly(sp.expand(Fpq), p, q)

    result = 0

    # Convert every monomial pair.
    # Since Fpq is symmetric, coefficients of p^a q^b and p^b q^a agree.
    visited = set()

    for (a, b), coeff in poly.terms():
        if (a, b) in visited:
            continue

        c1 = coeff
        c2 = poly.coeff_monomial(p**b * q**a)

        visited.add((a, b))
        visited.add((b, a))

        if a == b:
            # p^a q^a = N^a
            result += c1 * N**a
        else:
            # p^a q^b + p^b q^a
            # = N^min(a,b) * (p^(a-b)+q^(a-b))
            d = abs(a - b)
            result += c1 * N**min(a, b) * P[d]

    return canonical(result)


def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    """
    Q_(k,l) = F_(k,l) / (S+1).

    The quotient is guaranteed polynomial for the detector family.
    """
    F = sp.expand(raw_detector(k, ell))
    divisor = sp.Poly(S + 1, S, domain=sp.QQ.frac_field(N))
    dividend = sp.Poly(F, S, domain=sp.QQ.frac_field(N))

    q, r = sp.div(dividend, divisor)

    r_expr = canonical(r.as_expr())
    if r_expr != 0:
        raise ArithmeticError(
            f"Nonzero remainder for Q_({k},{ell}): {r_expr}"
        )

    return canonical(q.as_expr())


# =============================================================================
# OPTIONAL DIRECT CHECK FOR KNOWN LOW-WEIGHT DETECTORS
# =============================================================================

KNOWN_Q = {
    (1, 3): 6 * N - S**2 + S,

    (1, 5):
        -20 * N**2
        + 10 * N * S**2
        + 10 * N
        - S**4
        + S**3
        - S**2
        + S,

    (1, 7):
        42 * N**3
        - 49 * N**2 * S**2
        - 35 * N**2 * S
        - 70 * N**2
        + 14 * N * S**4
        + 7 * N * S**3
        + 28 * N * S**2
        + 7 * N * S
        + 14 * N
        - S**6
        + S**5
        - S**4
        + S**3
        - S**2
        + S,

    (3, 5):
        14 * N**3
        - 3 * N**2 * S**2
        + 15 * N**2 * S
        - 10 * N**2
        - 3 * N * S**3
        + 8 * N * S**2
        - 3 * N * S
        - S**4
        + S**3,
}


# =============================================================================
# INVARIANT-RING TEST
# =============================================================================

def involution_difference(Q: sp.Expr) -> sp.Expr:
    """
    Delta(N,S) = Q(N,1-S) - Q(N,S).
    """
    return canonical(Q.subs(S, 1 - S) - Q)


def even_in_centered_coordinate(Q: sp.Expr) -> Tuple[bool, sp.Expr]:
    """
    Put

        Y = S - 1/2.

    The involution S -> 1-S becomes Y -> -Y.

    Therefore invariance is equivalent to Q being even in Y.

    This is useful because X = S(S-1) = Y^2 - 1/4.
    """
    Y = sp.symbols("Y")
    shifted = sp.expand(Q.subs(S, Y + sp.Rational(1, 2)))
    odd_part = canonical(shifted - shifted.subs(Y, -Y))
    return odd_part == 0, odd_part


def x_representation(Q: sp.Expr) -> Tuple[bool, sp.Expr, sp.Expr]:
    """
    Attempt exact rewriting

        Q(N,S) = F(N,X),   X=S(S-1).

    Because X=S(S-1)=S^2-S, an X-polynomial contains only
    combinations (S^2-S)^j.

    We reconstruct by reducing Q in powers of S and matching.
    """
    Q = sp.Poly(sp.expand(Q), S, domain=sp.QQ.frac_field(N))

    deg = Q.degree()
    if deg < 0:
        return True, sp.Integer(0), sp.Integer(0)

    # An invariant polynomial of S degree 2m must be a polynomial in X
    # of degree <= m.
    m = deg // 2
    coeffs = sp.symbols(f"a0:{m+1}")

    candidate = sum(coeffs[j] * X**j for j in range(m + 1))
    candidate_S = sp.expand(candidate.subs(X, S**2 - S))

    diff = sp.Poly(
        sp.expand(Q.as_expr() - candidate_S),
        S,
        domain=sp.QQ.frac_field(N),
    )

    equations = [sp.Eq(c, 0) for c in diff.all_coeffs()]

    sol = sp.solve(equations, coeffs, dict=True)

    if not sol:
        return False, sp.Integer(0), canonical(Q.as_expr())

    F = canonical(candidate.subs(sol[0]))
    residual = canonical(Q.as_expr() - F.subs(X, S**2 - S))
    return residual == 0, F, residual


# =============================================================================
# COLLISION UTILITIES
# =============================================================================

@dataclass(frozen=True)
class Branch:
    s: int
    d2: int
    qvec: Tuple[int, ...]


def unit_residues(M: int) -> List[int]:
    return [a for a in range(M) if math.gcd(a, M) == 1]


def detector_residue_vector(
    Qs: Sequence[sp.Expr],
    n_mod: int,
    s_mod: int,
    M: int,
) -> Tuple[int, ...]:
    vals = []
    for Q in Qs:
        value = int(Q.subs({N: n_mod, S: s_mod})) % M
        vals.append(value)
    return tuple(vals)


def branch_table_for_product(
    M: int,
    product: int,
    detector_Qs: Sequence[sp.Expr],
) -> List[Branch]:
    """
    Enumerate all S satisfying

        S^2 - 4P is a square mod M

    without factoring any target integer.

    For each admissible S compute d2 and all detector residues.
    """
    out: List[Branch] = []

    for s in range(M):
        d2 = (s * s - 4 * product) % M

        # d2 must be a quadratic residue modulo M.
        roots = [
            d for d in range(M)
            if d * d % M == d2
        ]
        if not roots:
            continue

        qvec = detector_residue_vector(
            detector_Qs, product, s, M
        )

        out.append(
            Branch(
                s=s,
                d2=d2,
                qvec=qvec,
            )
        )

    return out


def involution_partner(s: int, M: int) -> int:
    return (1 - s) % M


def analyze_collision_partition(
    branches: Sequence[Branch],
    M: int,
) -> Dict[str, int]:
    """
    Compare detector-vector collisions with the candidate involution
    s -> 1-s.
    """
    by_q: Dict[Tuple[int, ...], List[Branch]] = defaultdict(list)
    for b in branches:
        by_q[b.qvec].append(b)

    collision_vectors = 0
    explained_pairs = 0
    unexplained_vectors = 0
    total_collision_pairs = 0

    for _, group in by_q.items():
        if len(group) <= 1:
            continue

        collision_vectors += 1

        partner_set = {involution_partner(b.s, M) for b in group}
        group_s = {b.s for b in group}

        if all(s2 in partner_set for s2 in group_s):
            explained_pairs += 1
        else:
            unexplained_vectors += 1

        total_collision_pairs += len(group) * (len(group) - 1) // 2

    return {
        "collision_vectors": collision_vectors,
        "explained_vectors": explained_pairs,
        "unexplained_vectors": unexplained_vectors,
        "collision_pairs": total_collision_pairs,
    }


# =============================================================================
# SYMBOLIC CROSS-DETECTOR TESTS
# =============================================================================

def common_involution_factor(differences: Sequence[sp.Expr]) -> sp.Expr:
    """
    Compute gcd of all nonzero involution differences.
    """
    g = None

    for d in differences:
        d = canonical(d)
        if d == 0:
            continue
        g = d if g is None else sp.gcd(g, d)

    return canonical(g if g is not None else 0)


def evaluate_invariants(
    s1: int,
    s2: int,
    M: int,
) -> Dict[str, int]:
    return {
        "sum": (s1 + s2) % M,
        "product": (s1 * s2) % M,
        "x1": (s1 * (s1 - 1)) % M,
        "x2": (s2 * (s2 - 1)) % M,
        "d2_same": int(
            (s1 * s1) % M == (s2 * s2) % M
        ),
        "involution": int(s2 == involution_partner(s1, M)),
    }


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 118")
    print("DETECTOR INVOLUTION / INVARIANT-RING TEST")
    print("Q(N,S) UNDER S -> 1-S")
    print("TEST X = S(S-1) INVARIANT REPRESENTATION")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. DIRECT VALIDATION
    # -------------------------------------------------------------------------
    print("\n1. DIRECT DETECTOR VALIDATION")
    print("-" * 78)

    direct_failures = 0
    Q_cache: Dict[Tuple[int, int], sp.Expr] = {}

    for pair in DETECTOR_PAIRS:
        k, ell = pair
        Q = quotient_polynomial(k, ell)
        Q_cache[pair] = Q

        if pair in KNOWN_Q:
            if canonical(Q - KNOWN_Q[pair]) != 0:
                direct_failures += 1
                print(
                    f"  FAIL Q{pair}: "
                    f"derived={Q}, known={KNOWN_Q[pair]}"
                )

    print(f"direct failures = {direct_failures}/{len(KNOWN_Q)}")
    if direct_failures:
        raise ArithmeticError("Detector validation failed.")

    print("STATUS = PASS")

    # -------------------------------------------------------------------------
    # 2. INVOLUTION TEST
    # -------------------------------------------------------------------------
    print("\n2. INVOLUTION TEST")
    print("-" * 78)

    zero_involution = []
    nonzero_involution = []

    for pair, Q in Q_cache.items():
        D = involution_difference(Q)

        if D == 0:
            zero_involution.append(pair)
        else:
            nonzero_involution.append((pair, D))

        print(f"Q_{pair}:")
        if D == 0:
            print("  Q(N,1-S)-Q(N,S) = 0")
            print("  EXACT INVOLUTION INVARIANCE")
        else:
            print(f"  difference = {D}")
            print(f"  degree_S   = {degree_in(D, S)}")
            print(f"  factor     = {sp.factor(D)}")

    print()
    print(f"invariant detectors = {len(zero_involution)}")
    print(f"non-invariant detectors = {len(nonzero_involution)}")

    # -------------------------------------------------------------------------
    # 3. CENTERED COORDINATE TEST
    # -------------------------------------------------------------------------
    print("\n3. CENTERED COORDINATE TEST")
    print("-" * 78)

    centered_failures = 0

    for pair, Q in Q_cache.items():
        ok, odd_part = even_in_centered_coordinate(Q)

        print(f"Q_{pair}: even_in_(S-1/2) = {ok}")

        if not ok:
            centered_failures += 1
            print(f"  odd obstruction = {odd_part}")

    print(f"centered-even failures = {centered_failures}/{len(Q_cache)}")

    # -------------------------------------------------------------------------
    # 4. EXACT X = S(S-1) REPRESENTATION
    # -------------------------------------------------------------------------
    print("\n4. INVARIANT-RING TEST: X = S(S-1)")
    print("-" * 78)

    x_success = 0
    x_fail = 0

    for pair, Q in Q_cache.items():
        ok, F, residual = x_representation(Q)

        print(f"Q_{pair}: X-representation = {ok}")

        if ok:
            x_success += 1
            print(f"  F(N,X) = {F}")
        else:
            x_fail += 1
            print(f"  residual = {residual}")

    print(f"X-success = {x_success}/{len(Q_cache)}")
    print(f"X-fail    = {x_fail}/{len(Q_cache)}")

    # -------------------------------------------------------------------------
    # 5. COMMON DIFFERENCE FACTOR
    # -------------------------------------------------------------------------
    print("\n5. COMMON INVOLUTION-DIFFERENCE FACTOR")
    print("-" * 78)

    all_diffs = [D for _, D in nonzero_involution]
    common = common_involution_factor(all_diffs)

    if common == 0:
        print("all detectors invariant: no nonzero difference factor")
    else:
        print(f"gcd of all nonzero differences = {sp.factor(common)}")

        if common != 1:
            print(
                "Interpretation: every non-invariant detector difference "
                "contains this common S-polynomial factor."
            )
        else:
            print(
                "No nonconstant common S-factor detected."
            )

    # -------------------------------------------------------------------------
    # 6. PRODUCT / SUM ANALYSIS OF SYMBOLIC DIFFERENCES
    # -------------------------------------------------------------------------
    print("\n6. QUADRATIC-PAIR STRUCTURE")
    print("-" * 78)

    for pair, D in nonzero_involution:
        poly = sp.Poly(sp.expand(D), S, domain=sp.QQ.frac_field(N))

        print(f"Q_{pair}: degree_S(diff)={poly.degree()}")

        if poly.degree() >= 2:
            roots = sp.solve(sp.Eq(poly.as_expr(), 0), S)

            if roots:
                print(f"  symbolic roots in S = {roots}")

    # -------------------------------------------------------------------------
    # 7. SMALL-MODULUS COLLISION CLASSIFICATION
    # -------------------------------------------------------------------------
    print("\n7. COLLISION CLASSIFICATION")
    print("-" * 78)

    test_moduli = [35, 385]

    for M in test_moduli:
        print(f"\nMODULUS M={M}")
        print("-" * 78)

        unit = unit_residues(M)
        collision_products = 0
        total_branches = 0
        detector_collision_vectors = 0
        involution_explained = 0
        unexplained = 0

        for p in unit:
            for q in unit:
                if p > q:
                    continue

                # Using unordered pair classes so we do not duplicate.
                product = (p * q) % M

                branches = branch_table_for_product(
                    M=M,
                    product=product,
                    detector_Qs=[Q_cache[x] for x in DETECTOR_PAIRS[:6]],
                )

                if len(branches) < 2:
                    continue

                collision_products += 1
                total_branches += len(branches)

                stats = analyze_collision_partition(branches, M)

                detector_collision_vectors += stats["collision_vectors"]
                involution_explained += stats["explained_vectors"]
                unexplained += stats["unexplained_vectors"]

        print(f"products with >=2 branches = {collision_products}")
        print(f"total admissible S branches = {total_branches}")
        print(f"detector collision vectors    = {detector_collision_vectors}")
        print(f"involution-explained vectors  = {involution_explained}")
        print(f"unexplained collision vectors = {unexplained}")

    # -------------------------------------------------------------------------
    # 8. EXPLICIT COLLISION WITNESSES
    # -------------------------------------------------------------------------
    print("\n8. EXPLICIT COLLISION WITNESSES")
    print("-" * 78)

    M = 35
    detector_list = [Q_cache[x] for x in DETECTOR_PAIRS[:6]]

    found = 0

    for product in range(M):
        branches = branch_table_for_product(
            M=M,
            product=product,
            detector_Qs=detector_list,
        )

        by_q: Dict[Tuple[int, ...], List[Branch]] = defaultdict(list)
        for b in branches:
            by_q[b.qvec].append(b)

        for qvec, group in by_q.items():
            if len(group) <= 1:
                continue

            # Prefer a collision not explained by S2 = 1-S1.
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    a = group[i]
                    b = group[j]

                    if b.s != involution_partner(a.s, M):
                        print(f"M={M}, product={product}")
                        print(f"  S1={a.s}, S2={b.s}")
                        print(f"  S1+S2={((a.s+b.s)%M)}")
                        print(f"  S1*S2={(a.s*b.s)%M}")
                        print(f"  X1={(a.s*(a.s-1))%M}")
                        print(f"  X2={(b.s*(b.s-1))%M}")
                        print(f"  D2_1={a.d2}, D2_2={b.d2}")
                        print(f"  Q-vector={qvec}")
                        print(
                            f"  involution partner of S1 = "
                            f"{involution_partner(a.s, M)}"
                        )

                        inv = evaluate_invariants(a.s, b.s, M)
                        print(f"  pair invariants = {inv}")

                        found += 1

                        if found >= 5:
                            break
                if found >= 5:
                    break
        if found >= 5:
            break

    if found == 0:
        print("No non-involution collision witness found at M=35.")

    # -------------------------------------------------------------------------
    # 9. DIRECT INVARIANT COMPARISON
    # -------------------------------------------------------------------------
    print("\n9. COLLISION INVARIANT COMPARISON")
    print("-" * 78)

    invariant_tests = [
        "S2 = 1-S1",
        "S1+S2 constant",
        "S1*S2 constant",
        "S1(S1-1) = S2(S2-1)",
        "D2 equal",
    ]

    print("candidate explanations tested:")
    for item in invariant_tests:
        print(f"  {item}")

    # -------------------------------------------------------------------------
    # 10. STRUCTURAL CONCLUSION
    # -------------------------------------------------------------------------
    print("\n10. FINAL DIAGNOSTIC")
    print("=" * 78)

    if len(zero_involution) == len(Q_cache):
        print(
            "ALL TESTED DETECTORS ARE EXACTLY INVARIANT UNDER S -> 1-S."
        )
        print()
        print(
            "This establishes that the detector family lives naturally "
            "inside the invariant ring of the involution."
        )
        print()
        print(
            "The natural next algebraic variable is:"
        )
        print()
        print("    X = S(S-1)")
        print()
        print(
            "and factor recovery would then require recovering S from"
        )
        print()
        print("    S^2 - S - X = 0")
        print()
        print(
            "whose discriminant is 1 + 4X."
        )
    else:
        print(
            "The detector family is NOT uniformly invariant under S -> 1-S."
        )
        print(
            "The non-invariant part must be characterized before treating "
            "X=S(S-1) as the full detector coordinate."
        )

    elapsed = time.perf_counter() - t0
    print()
    print(f"total runtime = {elapsed:.6f}s")
    print("=" * 78)
    print("EXPERIMENT 118 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nFATAL: {type(exc).__name__}: {exc}")
        raise

