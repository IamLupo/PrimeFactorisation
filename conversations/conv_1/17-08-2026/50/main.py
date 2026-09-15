#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 117
QUADRATIC-BRANCH COMPRESSION / DETECTOR INFORMATION MATRIX
EXACT RESIDUE-LEVEL ANALYSIS

GOAL
----
Experiment 116R showed that, for fixed N mod M, different trace branches
S mod M can produce different paper-detector vectors.

This experiment asks the sharper question:

    For fixed N mod M,

        S  --->  (Q13,Q15,Q17,Q35,Q37,Q57)

    how much information is retained?

We measure:

  1. Number of distinct trace branches S.
  2. Number of distinct D^2=(S^2-4N) branches.
  3. Number of distinct full detector vectors.
  4. Trace collisions under the full detector vector.
  5. D^2 collisions under the full detector vector.
  6. Whether the detector vector determines S.
  7. Whether the detector vector determines D^2.
  8. Minimal detector subsets that determine S.
  9. Minimal detector subsets that determine D^2.
 10. Explicit branch-collision witnesses.

IMPORTANT
---------
This is a residue-algebra experiment.

It does NOT:
  - factor N;
  - search factor pairs;
  - use SKLEARN;
  - use CSV;
  - depend on a prime population;
  - claim an N-only algorithm.

For each modulus M, all unit residue pairs are analyzed directly.

Supported moduli are intentionally limited so the computation stays
memory-safe. M=55055 is sampled rather than exhaustively enumerated.
==============================================================================
"""

from __future__ import annotations

from itertools import combinations
from math import gcd
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy as sp


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

EXHAUSTIVE_MODULI = [35, 385, 5005]
SAMPLED_MODULI = [55055]

# For large sampled moduli:
SAMPLE_PRODUCTS = 250
SAMPLE_P_PER_PRODUCT = 400

DETECTOR_NAMES = (
    "Q13",
    "Q15",
    "Q17",
    "Q35",
    "Q37",
    "Q57",
)

# Six exact paper quotient formulas from the established experiments.
N, S = sp.symbols("N S")

Q_EXPR = {
    "Q13": 6 * N - S**2 + S,

    "Q15": (
        -20 * N**2
        + 10 * N * S**2
        + 10 * N
        - S**4
        + S**3
        - S**2
        + S
    ),

    "Q17": (
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
        + S
    ),

    "Q35": (
        14 * N**3
        - 3 * N**2 * S**2
        + 15 * N**2 * S
        - 10 * N**2
        - 3 * N * S**3
        + 8 * N * S**2
        - 3 * N * S
        - S**4
        + S**3
    ),

    "Q37": (
        -36 * N**4
        + 22 * N**3 * S**2
        - 28 * N**3 * S
        + 70 * N**3
        - 3 * N**2 * S**4
        + 21 * N**2 * S**3
        - 35 * N**2 * S**2
        + 35 * N**2 * S
        - 14 * N**2
        - 3 * N * S**5
        + 10 * N * S**4
        - 10 * N * S**3
        + 10 * N * S**2
        - 3 * N * S
        - S**6
        + S**5
        - S**4
        + S**3
    ),

    "Q57": (
        22 * N**5
        - 5 * N**4 * S**2
        + 45 * N**4 * S
        - 60 * N**4
        - 10 * N**3 * S**3
        + 60 * N**3 * S**2
        - 70 * N**3 * S
        + 14 * N**3
        - 10 * N**2 * S**4
        + 40 * N**2 * S**3
        - 33 * N**2 * S**2
        + 5 * N**2 * S
        - 5 * N * S**5
        + 12 * N * S**4
        - 5 * N * S**3
        - S**6
        + S**5
    ),
}


# ---------------------------------------------------------------------------
# Basic modular helpers
# ---------------------------------------------------------------------------

def unit_residues(M: int) -> List[int]:
    return [x for x in range(1, M) if gcd(x, M) == 1]


def mod_expr(expr: sp.Expr, n: int, s: int, M: int) -> int:
    """
    Exact integer evaluation followed by reduction mod M.
    SymPy is only used for correctness; no floating point is involved.
    """
    value = expr.subs({N: n, S: s})
    return int(value) % M


def detector_vector(n: int, s: int, M: int) -> Tuple[int, ...]:
    return tuple(
        mod_expr(Q_EXPR[name], n, s, M)
        for name in DETECTOR_NAMES
    )


def d2_value(n: int, s: int, M: int) -> int:
    return (s * s - 4 * n) % M


# ---------------------------------------------------------------------------
# Residue-pair generation
# ---------------------------------------------------------------------------

def multiplicative_collisions(
    M: int,
    *,
    sample_products: int | None = None,
    sample_p_per_product: int | None = None,
) -> Iterable[Tuple[int, List[Tuple[int, int]]]]:
    """
    Yield (product residue, [(p,q), ...]) for unit residue pairs.

    Exhaustive mode:
        every unit p for every product P.

    Sample mode:
        a deterministic subset of products and p-values.
    """
    units = unit_residues(M)
    if not units:
        return

    inv = {x: pow(x, -1, M) for x in units}

    if sample_products is None:
        products = units
    else:
        step = max(1, len(units) // sample_products)
        products = units[::step][:sample_products]

    for P in products:
        if sample_p_per_product is None:
            ps = units
        else:
            step = max(1, len(units) // sample_p_per_product)
            ps = units[::step][:sample_p_per_product]

        pairs = []
        for p in ps:
            q = (P * inv[p]) % M
            pairs.append((p, q))

        yield P, pairs


# ---------------------------------------------------------------------------
# Per-product branch analysis
# ---------------------------------------------------------------------------

def analyze_product(
    M: int,
    P: int,
    pairs: Sequence[Tuple[int, int]],
) -> dict:
    """
    Analyze all sampled/exhaustive pairs with pq = P mod M.

    Since (p,q) and (q,p) produce the same S and D^2, we collapse at the
    trace level for information analysis.
    """
    branch_by_s: Dict[int, dict] = {}

    for p, q in pairs:
        s = (p + q) % M
        d2 = d2_value(P, s, M)
        qv = detector_vector(P, s, M)

        branch_by_s[s] = {
            "d2": d2,
            "qv": qv,
            "pair": (p, q),
        }

    if not branch_by_s:
        return {
            "P": P,
            "branches": 0,
            "distinct_d2": 0,
            "distinct_qv": 0,
            "trace_injective": True,
            "d2_injective": True,
            "branch_by_s": {},
        }

    s_values = sorted(branch_by_s)
    q_vectors = [branch_by_s[s]["qv"] for s in s_values]
    d2_values = [branch_by_s[s]["d2"] for s in s_values]

    qv_to_s: Dict[Tuple[int, ...], List[int]] = {}
    d2_to_s: Dict[int, List[int]] = {}

    for s in s_values:
        qv = branch_by_s[s]["qv"]
        d2 = branch_by_s[s]["d2"]

        qv_to_s.setdefault(qv, []).append(s)
        d2_to_s.setdefault(d2, []).append(s)

    trace_injective = len(qv_to_s) == len(s_values)
    d2_injective = len(qv_to_s) == len(set(d2_values))

    return {
        "P": P,
        "branches": len(s_values),
        "distinct_d2": len(set(d2_values)),
        "distinct_qv": len(qv_to_s),
        "trace_injective": trace_injective,
        "d2_injective": d2_injective,
        "qv_to_s": qv_to_s,
        "d2_to_s": d2_to_s,
        "branch_by_s": branch_by_s,
    }


# ---------------------------------------------------------------------------
# Detector subset analysis
# ---------------------------------------------------------------------------

def subset_signature(
    qv: Tuple[int, ...],
    indices: Sequence[int],
) -> Tuple[int, ...]:
    return tuple(qv[i] for i in indices)


def subset_is_injective(
    branch_by_s: Dict[int, dict],
    indices: Sequence[int],
) -> bool:
    seen = set()

    for info in branch_by_s.values():
        sig = subset_signature(info["qv"], indices)
        if sig in seen:
            return False
        seen.add(sig)

    return True


def minimal_injective_subsets(
    branch_by_s: Dict[int, dict],
) -> Tuple[List[Tuple[str, ...]], List[Tuple[str, ...]]]:
    """
    Returns:
        subsets determining S
        subsets determining D^2

    We test all 2^6 detector subsets, including singleton subsets.
    """
    trace_minimal: List[Tuple[str, ...]] = []
    d2_minimal: List[Tuple[str, ...]] = []

    for r in range(1, len(DETECTOR_NAMES) + 1):
        for indices in combinations(range(len(DETECTOR_NAMES)), r):
            names = tuple(DETECTOR_NAMES[i] for i in indices)

            # Q-vector determines S iff it is injective over S.
            if subset_is_injective(branch_by_s, indices):
                trace_minimal.append(names)

            # Q-vector determines D^2 iff equal signatures imply equal d2.
            sig_to_d2: Dict[Tuple[int, ...], int] = {}
            ok = True

            for info in branch_by_s.values():
                sig = subset_signature(info["qv"], indices)
                d2 = info["d2"]

                if sig in sig_to_d2 and sig_to_d2[sig] != d2:
                    ok = False
                    break

                sig_to_d2[sig] = d2

            if ok:
                d2_minimal.append(names)

        # Once the smallest size is found, stop searching larger sizes.
        if trace_minimal and d2_minimal:
            min_trace_size = len(trace_minimal[0])
            min_d2_size = len(d2_minimal[0])

            trace_minimal = [
                x for x in trace_minimal if len(x) == min_trace_size
            ]
            d2_minimal = [
                x for x in d2_minimal if len(x) == min_d2_size
            ]
            break

    return trace_minimal, d2_minimal


# ---------------------------------------------------------------------------
# Global analysis
# ---------------------------------------------------------------------------

def analyze_modulus(
    M: int,
    *,
    sample_products: int | None = None,
    sample_p_per_product: int | None = None,
) -> None:
    print("=" * 78)
    print(f"MODULUS M={M}")
    print("=" * 78)

    units = unit_residues(M)
    print(f"unit residues = {len(units)}")

    total_products = 0
    analyzed_products = 0
    total_branches = 0
    total_q_vectors = 0
    total_d2 = 0

    trace_collisions = 0
    d2_collisions = 0

    best_trace_collision = None
    best_d2_collision = None

    subset_trace_counts: Dict[Tuple[str, ...], int] = {}
    subset_d2_counts: Dict[Tuple[str, ...], int] = {}

    for P, pairs in multiplicative_collisions(
        M,
        sample_products=sample_products,
        sample_p_per_product=sample_p_per_product,
    ):
        total_products += 1

        result = analyze_product(M, P, pairs)

        if result["branches"] <= 1:
            continue

        analyzed_products += 1
        total_branches += result["branches"]
        total_q_vectors += result["distinct_qv"]
        total_d2 += result["distinct_d2"]

        if not result["trace_injective"]:
            trace_collisions += 1

            collision_groups = [
                (qv, ss)
                for qv, ss in result["qv_to_s"].items()
                if len(ss) > 1
            ]

            collision_groups.sort(key=lambda x: (-len(x[1]), x[1]))

            if collision_groups:
                qv, ss = collision_groups[0]
                candidate = (len(ss), P, ss, qv)

                if (
                    best_trace_collision is None
                    or candidate[0] > best_trace_collision[0]
                ):
                    best_trace_collision = candidate

        # D^2 distinguishability is weaker than S distinguishability.
        # Here we test whether equal detector vectors can correspond to
        # different D^2 values.
        bad_d2 = False
        for qv, ss in result["qv_to_s"].items():
            values = {
                result["branch_by_s"][s]["d2"]
                for s in ss
            }
            if len(values) > 1:
                bad_d2 = True
                break

        if bad_d2:
            d2_collisions += 1

        # Find minimal useful detector subsets for a representative branch
        # class when it has enough distinct traces to make the question
        # meaningful.
        if analyzed_products <= 12 and result["branches"] >= 3:
            tr_min, d2_min = minimal_injective_subsets(result["branch_by_s"])

            for x in tr_min:
                subset_trace_counts[x] = subset_trace_counts.get(x, 0) + 1

            for x in d2_min:
                subset_d2_counts[x] = subset_d2_counts.get(x, 0) + 1

    print()
    print("GLOBAL BRANCH STATISTICS")
    print("-" * 78)
    print(f"products examined             = {total_products}")
    print(f"products with >=2 S branches  = {analyzed_products}")
    print(f"trace branches                = {total_branches}")
    print(f"distinct detector vectors     = {total_q_vectors}")
    print(f"distinct D^2 branches         = {total_d2}")
    print(f"products with S-collisions    = {trace_collisions}")
    print(f"products with D^2 collisions  = {d2_collisions}")

    if analyzed_products:
        print(
            "average detector compression = "
            f"{total_q_vectors / total_branches:.6f}"
        )
        print(
            "average D2 compression       = "
            f"{total_d2 / total_branches:.6f}"
        )

    print()
    print("DETECTOR VECTOR VS TRACE")
    print("-" * 78)

    if trace_collisions == 0:
        print("full six-detector vector is injective on every analyzed trace set")
    else:
        print(
            "full six-detector vector is NOT injective on all analyzed "
            f"trace sets ({trace_collisions} collision products)"
        )

    print()
    print("DETECTOR VECTOR VS D^2")
    print("-" * 78)

    if d2_collisions == 0:
        print("full detector vector determines D^2 on every analyzed branch set")
    else:
        print(
            "full detector vector has D^2 ambiguity on "
            f"{d2_collisions} analyzed products"
        )

    print()
    print("MINIMAL DETECTOR SUBSETS")
    print("-" * 78)

    if subset_trace_counts:
        ordered = sorted(
            subset_trace_counts.items(),
            key=lambda kv: (-kv[1], len(kv[0]), kv[0]),
        )[:12]

        print("most frequent smallest subsets determining S:")
        for names, count in ordered:
            print(f"  {names}: {count} representative products")

    else:
        print("no trace-determining subset found in representative products")

    if subset_d2_counts:
        ordered = sorted(
            subset_d2_counts.items(),
            key=lambda kv: (-kv[1], len(kv[0]), kv[0]),
        )[:12]

        print("most frequent smallest subsets determining D^2:")
        for names, count in ordered:
            print(f"  {names}: {count} representative products")

    else:
        print("no D^2-determining subset found in representative products")

    print()
    print("EXPLICIT COLLISION WITNESS")
    print("-" * 78)

    if best_trace_collision is None:
        print("no full-vector trace collision found")
    else:
        multiplicity, P, ss, qv = best_trace_collision
        print(f"product residue P = {P}")
        print(f"trace collision multiplicity = {multiplicity}")
        print(f"S values = {ss}")
        print(f"Q vector = {qv}")

        # Print D^2 for the colliding traces.
        for s in ss:
            info = analyze_product(
                M,
                P,
                [(s, (P * pow(s, -1, M)) % M)]
                if gcd(s, M) == 1
                else [(0, 0)],
            )

        print("colliding branches:")
        # Pull branch information directly by rebuilding the branch table.
        inv = {
            x: pow(x, -1, M)
            for x in unit_residues(M)
        }

        for s in ss:
            # Find p satisfying p+q=s with pq=P.
            found = None
            for p in unit_residues(M):
                q = (P * inv[p]) % M
                if (p + q) % M == s:
                    found = (p, q)
                    break

            if found is None:
                print(f"  S={s}: no representative pair found")
                continue

            d2 = d2_value(P, s, M)
            print(f"  S={s:>6} D2={d2:>6} pair={found}")

    print()


# ---------------------------------------------------------------------------
# Direct algebra sanity checks
# ---------------------------------------------------------------------------

def sanity_check() -> None:
    print("=" * 78)
    print("DIRECT ALGEBRA SANITY CHECK")
    print("=" * 78)

    samples = [
        (3, 7, 35),
        (8, 13, 385),
        (11, 17, 5005),
    ]

    for p, q, M0 in samples:
        n = p * q
        s = p + q
        t = s * s - 3 * n

        assert (t - n) % M0 == ((p - q) ** 2) % M0

        for name, expr in Q_EXPR.items():
            direct = mod_expr(expr, n, s, M0)

            # For the basic detector, independently verify Q13 from its
            # simple closed form.
            if name == "Q13":
                independent = (6 * n - s * s + s) % M0
                assert direct == independent

    print("sanity checks = PASS")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 78)
    print("KAPPA EXPERIMENT 117")
    print("QUADRATIC-BRANCH COMPRESSION / DETECTOR INFORMATION MATRIX")
    print("EXACT RESIDUE-LEVEL ANALYSIS")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    sanity_check()

    for M in EXHAUSTIVE_MODULI:
        analyze_modulus(M)

    for M in SAMPLED_MODULI:
        print(
            f"Sampling large modulus M={M}: "
            f"products={SAMPLE_PRODUCTS}, "
            f"p-per-product={SAMPLE_P_PER_PRODUCT}"
        )
        analyze_modulus(
            M,
            sample_products=SAMPLE_PRODUCTS,
            sample_p_per_product=SAMPLE_P_PER_PRODUCT,
        )

    print("=" * 78)
    print("EXPERIMENT 117 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

