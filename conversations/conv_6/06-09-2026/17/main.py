#!/usr/bin/env python3

"""
START EXPERIMENT 174

MRV + FORWARD CHECKING ON THE COMPLETE 5x5 RESIDUE CSP

Variables:

    A_i = p mod R1[i]
    B_j = q mod R2[j]

Constraints:

    (A_i, B_j) must belong to the observed C relation R_ij.

Previous experiments:

    162:
        C depends only on z_ij = A_i * B_j.

    167:
        rarest-cell anchoring is useful.

    168-169:
        CRT parameterization reduces Cartesian explosion.

    170:
        pairwise/path consistency did not materially improve
        the surviving domains.

    171-173:
        2x2 / product rank-1 closure removes ambiguity but
        costs a lot and still leaves broad domains.

Experiment 174 changes the search architecture.

Instead of:

    anchor
      -> AC3
      -> many CRT searches
      -> repeat

we use:

    MRV
      -> assign one residue
      -> forward-check all neighbors
      -> propagate singleton domains
      -> choose next smallest domain
      -> recurse
      -> only at a COMPLETE residue assignment:
             CRT p
             CRT q
             test p*q=n

This measures whether the CSP itself contains enough higher-order
information to collapse with ordinary search, without paying for
repeated CRT enumeration.

Important:

    The true solution must survive.

This experiment reports:

    - initial relation sizes
    - AC3 reduction
    - MRV node count
    - backtracks
    - maximum depth
    - complete residue assignments examined
    - CRT tests
    - exact recovery time

FINISHED EXPERIMENT 174
"""

from __future__ import annotations

import itertools
import math
import random
import time

from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

from sympy import randprime


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS_LIST = [30, 36, 42, 48, 54]

SEED = 1742026

MAX_NODES = 5_000_000

MAX_COMPLETE_ASSIGNMENTS = 1_000_000

STOP_ON_EXACT = True


# ============================================================
# CRT
# ============================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:

    t = (
        (a2 - a1)
        * pow(m1, -1, m2)
    ) % m2

    x = a1 + m1 * t
    M = m1 * m2

    return x % M, M


def crt_many(
    residues: List[int],
    moduli: List[int],
) -> Tuple[int, int]:

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli,
    ):

        x, M = crt_pair(
            x,
            M,
            a,
            m,
        )

    return x, M


# ============================================================
# C
# ============================================================

def exact_C(
    n: int,
    r: int,
    s: int,
    a: int,
    b: int,
) -> int:

    z = a * b

    beta = (
        (n - z)
        * pow(r, -1, s)
    ) % s

    alpha = (
        (n - z)
        * pow(s, -1, r)
    ) % r

    return (
        r * beta
        + s * alpha
        + z
    ) // (r * s)


# ============================================================
# C GRID
# ============================================================

def make_C_grid(
    n: int,
    p: int,
    q: int,
) -> List[List[int]]:

    A = [
        p % r
        for r in R1
    ]

    B = [
        q % s
        for s in R2
    ]

    return [
        [
            exact_C(
                n,
                r,
                s,
                A[i],
                B[j],
            )
            for j, s in enumerate(R2)
        ]
        for i, r in enumerate(R1)
    ]


# ============================================================
# RELATIONS
# ============================================================

@dataclass
class Relation:

    pairs: Set[Tuple[int, int]]

    by_a: Dict[int, Set[int]]
    by_b: Dict[int, Set[int]]


def build_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> Relation:

    pairs = set()

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}

    for a in range(r):

        for b in range(s):

            if exact_C(
                n,
                r,
                s,
                a,
                b,
            ) != target_c:

                continue

            pairs.add(
                (a, b)
            )

            by_a.setdefault(
                a,
                set(),
            ).add(b)

            by_b.setdefault(
                b,
                set(),
            ).add(a)

    return Relation(
        pairs=pairs,
        by_a=by_a,
        by_b=by_b,
    )


def build_relations(
    n: int,
    c_grid: List[List[int]],
) -> Dict[
    Tuple[int, int],
    Relation,
]:

    result = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            result[(i, j)] = build_relation(
                n,
                r,
                s,
                c_grid[i][j],
            )

    return result


# ============================================================
# VARIABLES
# ============================================================

# A variables:
#     0 .. 4
#
# B variables:
#     5 .. 9

def is_A(
    v: int,
) -> bool:

    return v < len(R1)


def A_index(
    v: int,
) -> int:

    return v


def B_index(
    v: int,
) -> int:

    return v - len(R1)


def variable_name(
    v: int,
) -> str:

    if is_A(v):
        return f"A{v}"

    return f"B{B_index(v)}"


def variable_modulus(
    v: int,
) -> int:

    if is_A(v):
        return R1[A_index(v)]

    return R2[B_index(v)]


def initial_domain(
    v: int,
) -> Set[int]:

    return set(
        range(
            variable_modulus(v)
        )
    )


# ============================================================
# NEIGHBORS
# ============================================================

def neighbors(
    v: int,
) -> List[int]:

    if is_A(v):

        i = A_index(v)

        return [
            len(R1) + j
            for j in range(len(R2))
        ]

    j = B_index(v)

    return list(
        range(len(R1))
    )


def relation_for(
    u: int,
    v: int,
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> Relation:

    if is_A(u) and not is_A(v):

        return relations[
            (
                A_index(u),
                B_index(v),
            )
        ]

    if not is_A(u) and is_A(v):

        return relations[
            (
                A_index(v),
                B_index(u),
            )
        ]

    raise ValueError(
        "A-B relation expected"
    )


def supported_values(
    u: int,
    candidate_u: int,
    v: int,
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> Set[int]:

    rel = relation_for(
        u,
        v,
        relations,
    )

    if is_A(u):

        return set(
            rel.by_a.get(
                candidate_u,
                set(),
            )
        )

    return set(
        rel.by_b.get(
            candidate_u,
            set(),
        )
    )


# ============================================================
# AC3 / FORWARD CHECKING
# ============================================================

@dataclass
class SearchStats:

    nodes: int = 0

    assignments: int = 0

    backtracks: int = 0

    complete_assignments: int = 0

    crt_tests_p: int = 0

    crt_tests_q: int = 0

    exact_hits: int = 0

    max_depth: int = 0

    propagation_rounds: int = 0

    values_removed: int = 0


def propagate_domains(
    domains: Dict[int, Set[int]],
    assigned: Dict[int, int],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
    stats: SearchStats,
) -> bool:

    """
    Full AC-3 style propagation over the current domains.

    Assigned variables are represented as singleton domains.

    Every time an edge removes values, the affected neighbors
    are placed back into the queue.
    """

    queue = []

    for u in range(10):

        for v in neighbors(u):

            if u < v:
                queue.append(
                    (u, v)
                )
                queue.append(
                    (v, u)
                )

    while queue:

        u, v = queue.pop()

        stats.propagation_rounds += 1

        if not domains[u] or not domains[v]:
            return False

        rel = relation_for(
            u,
            v,
            relations,
        )

        new_u = set()

        if is_A(u):

            for a in domains[u]:

                if (
                    rel.by_a.get(
                        a,
                        set(),
                    )
                    & domains[v]
                ):

                    new_u.add(a)

        else:

            for b in domains[u]:

                if (
                    rel.by_b.get(
                        b,
                        set(),
                    )
                    & domains[v]
                ):

                    new_u.add(b)

        if new_u != domains[u]:

            stats.values_removed += (
                len(domains[u])
                - len(new_u)
            )

            domains[u] = new_u

            if not new_u:
                return False

            for w in neighbors(u):

                if w != v:

                    queue.append(
                        (w, u)
                    )

    return True


# ============================================================
# INITIAL AC3
# ============================================================

def initial_ac3(
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> Tuple[
    Dict[int, Set[int]],
    SearchStats,
]:

    stats = SearchStats()

    domains = {
        v: initial_domain(v)
        for v in range(10)
    }

    ok = propagate_domains(
        domains,
        {},
        relations,
        stats,
    )

    if not ok:

        raise RuntimeError(
            "Initial CSP is inconsistent"
        )

    return domains, stats


# ============================================================
# MRV SELECTION
# ============================================================

def choose_mrv_variable(
    domains: Dict[int, Set[int]],
    assigned: Dict[int, int],
) -> int:

    candidates = [
        v
        for v in range(10)
        if v not in assigned
    ]

    # MRV first.
    #
    # Tie-break:
    # choose variable with the largest number
    # of still-unassigned neighbors.
    #
    # This tends to expose more constraints sooner.

    return min(
        candidates,
        key=lambda v: (
            len(domains[v]),
            -sum(
                w not in assigned
                for w in neighbors(v)
            ),
            v,
        ),
    )


# ============================================================
# VALUE ORDERING
# ============================================================

def value_score(
    v: int,
    value: int,
    domains: Dict[int, Set[int]],
    assigned: Dict[int, int],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> int:

    """
    Least-constraining value score.

    Larger score = leaves more neighboring possibilities.

    We search values with the smallest score first only because
    this tends to find contradictions quickly.
    """

    score = 0

    for w in neighbors(v):

        if w in assigned:
            continue

        score += len(
            supported_values(
                v,
                value,
                w,
                relations,
            )
            & domains[w]
        )

    return score


# ============================================================
# EXACT ASSIGNMENT CHECK
# ============================================================

def complete_assignment_to_factors(
    assignment: Dict[int, int],
    n: int,
    stats: SearchStats,
) -> Tuple[int | None, int | None]:

    A = [
        assignment[i]
        for i in range(len(R1))
    ]

    B = [
        assignment[
            len(R1) + j
        ]
        for j in range(len(R2))
    ]

    p_candidate, M_p = crt_many(
        A,
        R1,
    )

    q_candidate, M_q = crt_many(
        B,
        R2,
    )

    stats.crt_tests_p += 1
    stats.crt_tests_q += 1

    # The complete CRT residues determine p modulo M_p and
    # q modulo M_q. We do not yet know which representative
    # is the actual factor.
    #
    # Search the smaller representative side first.

    # p = p0 + M_p*k
    #
    # q = q0 + M_q*l
    #
    # Because p*q=n, candidates are bounded by sqrt(n)
    # for the smaller factor.

    sqrt_n = math.isqrt(n)

    # Search whichever CRT modulus is larger first because
    # that produces fewer possible representatives.

    if M_p >= M_q:

        if M_p > sqrt_n:

            if (
                p_candidate > 1
                and n % p_candidate == 0
            ):

                q = n // p_candidate

                if q <= 1:
                    return None, None

                return (
                    p_candidate,
                    q,
                )

        else:

            kmax = (
                sqrt_n
                - p_candidate
            ) // M_p

            for k in range(
                max(0, kmax + 1)
            ):

                p = (
                    p_candidate
                    + k * M_p
                )

                if p <= 1:
                    continue

                if n % p != 0:
                    continue

                return (
                    p,
                    n // p,
                )

    else:

        if M_q > sqrt_n:

            if (
                q_candidate > 1
                and n % q_candidate == 0
            ):

                p = n // q_candidate

                if p <= 1:
                    return None, None

                return (
                    p,
                    q_candidate,
                )

        else:

            kmax = (
                sqrt_n
                - q_candidate
            ) // M_q

            for k in range(
                max(0, kmax + 1)
            ):

                q = (
                    q_candidate
                    + k * M_q
                )

                if q <= 1:
                    continue

                if n % q != 0:
                    continue

                return (
                    n // q,
                    q,
                )

    return None, None


# ============================================================
# FULL VERIFICATION
# ============================================================

def verify_factor(
    n: int,
    p: int,
    q: int,
    c_grid: List[List[int]],
) -> bool:

    if p <= 1 or q <= 1:
        return False

    if p * q != n:
        return False

    for i, r in enumerate(R1):

        a = p % r

        for j, s in enumerate(R2):

            b = q % s

            if exact_C(
                n,
                r,
                s,
                a,
                b,
            ) != c_grid[i][j]:

                return False

    return True


# ============================================================
# DFS
# ============================================================

def dfs(
    *,
    domains: Dict[int, Set[int]],
    assigned: Dict[int, int],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
    n: int,
    c_grid: List[List[int]],
    stats: SearchStats,
) -> Tuple[int | None, int | None]:

    stats.nodes += 1

    depth = len(
        assigned
    )

    stats.max_depth = max(
        stats.max_depth,
        depth,
    )

    if stats.nodes > MAX_NODES:

        raise RuntimeError(
            "MAX_NODES exceeded"
        )

    # --------------------------------------------------------
    # Complete assignment.
    # --------------------------------------------------------

    if depth == 10:

        stats.complete_assignments += 1

        if (
            stats.complete_assignments
            > MAX_COMPLETE_ASSIGNMENTS
        ):

            raise RuntimeError(
                "MAX_COMPLETE_ASSIGNMENTS exceeded"
            )

        p, q = (
            complete_assignment_to_factors(
                assigned,
                n,
                stats,
            )
        )

        if p is None:
            return None, None

        stats.exact_hits += 1

        if verify_factor(
            n,
            p,
            q,
            c_grid,
        ):

            return p, q

        return None, None

    # --------------------------------------------------------
    # MRV.
    # --------------------------------------------------------

    v = choose_mrv_variable(
        domains,
        assigned,
    )

    values = sorted(
        domains[v],
        key=lambda value: (
            value_score(
                v,
                value,
                domains,
                assigned,
                relations,
            ),
            value,
        ),
    )

    for value in values:

        stats.assignments += 1

        # ----------------------------------------------------
        # Copy domains.
        # ----------------------------------------------------

        child_domains = {
            u: set(d)
            for u, d in domains.items()
        }

        child_domains[v] = {
            value
        }

        child_assigned = dict(
            assigned
        )

        child_assigned[v] = value

        # ----------------------------------------------------
        # Forward checking / AC3.
        # ----------------------------------------------------

        ok = propagate_domains(
            child_domains,
            child_assigned,
            relations,
            stats,
        )

        if not ok:

            stats.backtracks += 1
            continue

        # ----------------------------------------------------
        # Synchronize singleton domains into assignment.
        # ----------------------------------------------------

        consistent = True

        changed = True

        while changed:

            changed = False

            for u in range(10):

                if (
                    u in child_assigned
                ):
                    continue

                if len(
                    child_domains[u]
                ) != 1:

                    continue

                only = next(
                    iter(
                        child_domains[u]
                    )
                )

                # Check against already assigned neighbors.
                for w in neighbors(u):

                    if (
                        w not in child_assigned
                    ):
                        continue

                    rel = relation_for(
                        u,
                        w,
                        relations,
                    )

                    if is_A(u):

                        if only not in (
                            rel.by_a.get(
                                child_assigned[w],
                                set(),
                            )
                        ):

                            consistent = False
                            break

                    else:

                        if only not in (
                            rel.by_b.get(
                                child_assigned[w],
                                set(),
                            )
                        ):

                            consistent = False
                            break

                if not consistent:
                    break

                child_assigned[u] = only
                changed = True

            if not consistent:
                break

            if changed:

                if not propagate_domains(
                    child_domains,
                    child_assigned,
                    relations,
                    stats,
                ):

                    consistent = False
                    break

        if not consistent:

            stats.backtracks += 1
            continue

        # ----------------------------------------------------
        # Recurse.
        # ----------------------------------------------------

        result = dfs(
            domains=child_domains,
            assigned=child_assigned,
            relations=relations,
            n=n,
            c_grid=c_grid,
            stats=stats,
        )

        if result[0] is not None:

            return result

    return None, None


# ============================================================
# SAMPLE GENERATION
# ============================================================

def make_semiprime(
    bits: int,
    rng: random.Random,
) -> Tuple[int, int, int]:

    lo = 1 << (
        bits // 2 - 1
    )

    hi = 1 << (
        bits // 2 + 1
    )

    p = int(
        randprime(
            lo,
            hi,
        )
    )

    q = int(
        randprime(
            lo,
            hi,
        )
    )

    while q == p:

        q = int(
            randprime(
                lo,
                hi,
            )
        )

    return (
        p,
        q,
        p * q,
    )


# ============================================================
# SAMPLE
# ============================================================

def run_sample(
    bits: int,
    sample_id: int,
    rng: random.Random,
) -> dict:

    p, q, n = make_semiprime(
        bits,
        rng,
    )

    true_A = [
        p % r
        for r in R1
    ]

    true_B = [
        q % s
        for s in R2
    ]

    c_grid = make_C_grid(
        n,
        p,
        q,
    )

    relations = build_relations(
        n,
        c_grid,
    )

    original_pairs = sum(
        len(rel.pairs)
        for rel in relations.values()
    )

    start_ac3 = time.perf_counter()

    domains, stats = initial_ac3(
        relations
    )

    ac3_time = (
        time.perf_counter()
        - start_ac3
    )

    initial_widths = [
        len(domains[v])
        for v in range(10)
    ]

    # --------------------------------------------------------
    # Check true residues survive initial propagation.
    # --------------------------------------------------------

    true_survive = all(
        true_A[i]
        in domains[i]
        for i in range(len(R1))
    ) and all(
        true_B[j]
        in domains[
            len(R1) + j
        ]
        for j in range(len(R2))
    )

    # --------------------------------------------------------
    # DFS.
    # --------------------------------------------------------

    search_start = time.perf_counter()

    try:

        recovered_p, recovered_q = dfs(
            domains=domains,
            assigned={},
            relations=relations,
            n=n,
            c_grid=c_grid,
            stats=stats,
        )

    except RuntimeError as exc:

        recovered_p = None
        recovered_q = None
        search_error = str(exc)

    else:

        search_error = None

    search_time = (
        time.perf_counter()
        - search_start
    )

    total_time = (
        ac3_time
        + search_time
    )

    recovered = (
        recovered_p is not None
        and recovered_q is not None
        and {
            recovered_p,
            recovered_q,
        }
        == {
            p,
            q,
        }
    )

    # --------------------------------------------------------
    # Output.
    # --------------------------------------------------------

    print()
    print("=" * 88)
    print(
        f"EXPERIMENT 174 SAMPLE {sample_id}"
    )
    print("=" * 88)

    print(
        f"bits              = {bits}"
    )

    print(
        f"p                 = {p}"
    )

    print(
        f"q                 = {q}"
    )

    print(
        f"n                 = {n}"
    )

    print(
        f"true A            = {true_A}"
    )

    print(
        f"true B            = {true_B}"
    )

    print()
    print(
        "C-grid:"
    )

    for row in c_grid:

        print(
            "   ",
            row,
        )

    print()
    print(
        "--- RELATIONS ---"
    )

    print(
        f"original relation pairs       = "
        f"{original_pairs}"
    )

    print()
    print(
        "--- INITIAL AC-3 ---"
    )

    print(
        f"domain widths                 = "
        f"{initial_widths}"
    )

    print(
        f"total domain values           = "
        f"{sum(initial_widths)}"
    )

    print(
        f"true residues survive         = "
        f"{true_survive}"
    )

    print(
        f"AC-3 propagation rounds       = "
        f"{stats.propagation_rounds}"
    )

    print(
        f"values removed                = "
        f"{stats.values_removed}"
    )

    print(
        f"AC-3 time                     = "
        f"{ac3_time:.6f}s"
    )

    print()
    print(
        "--- MRV SEARCH ---"
    )

    print(
        f"nodes                         = "
        f"{stats.nodes}"
    )

    print(
        f"assignments                   = "
        f"{stats.assignments}"
    )

    print(
        f"backtracks                    = "
        f"{stats.backtracks}"
    )

    print(
        f"max depth                     = "
        f"{stats.max_depth}"
    )

    print(
        f"complete assignments          = "
        f"{stats.complete_assignments}"
    )

    print(
        f"CRT p tests                   = "
        f"{stats.crt_tests_p}"
    )

    print(
        f"CRT q tests                   = "
        f"{stats.crt_tests_q}"
    )

    print(
        f"exact hits                    = "
        f"{stats.exact_hits}"
    )

    if search_error is not None:

        print(
            f"search status                 = "
            f"STOPPED: {search_error}"
        )

    else:

        print(
            f"search status                 = "
            f"completed"
        )

    print()
    print(
        f"recovered                     = "
        f"{recovered}"
    )

    print(
        f"recovered factors             = "
        f"{(recovered_p, recovered_q)}"
    )

    print()
    print(
        f"search time                   = "
        f"{search_time:.6f}s"
    )

    print(
        f"total time                    = "
        f"{total_time:.6f}s"
    )

    print("=" * 88)

    return {
        "bits": bits,
        "original_pairs": original_pairs,
        "domain_widths": initial_widths,
        "total_domain_values": sum(
            initial_widths
        ),
        "true_survive": true_survive,
        "values_removed": stats.values_removed,
        "nodes": stats.nodes,
        "assignments": stats.assignments,
        "backtracks": stats.backtracks,
        "max_depth": stats.max_depth,
        "complete_assignments": stats.complete_assignments,
        "crt_tests": (
            stats.crt_tests_p
            + stats.crt_tests_q
        ),
        "exact_hits": stats.exact_hits,
        "recovered": recovered,
        "total_time": total_time,
        "search_error": search_error,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(
        SEED
    )

    print(
        "START EXPERIMENT 174"
    )

    print()
    print(
        "R1 =",
        R1,
    )

    print(
        "R2 =",
        R2,
    )

    print(
        "BITS =",
        BITS_LIST,
    )

    print(
        "MAX_NODES =",
        MAX_NODES,
    )

    print(
        "MAX_COMPLETE_ASSIGNMENTS =",
        MAX_COMPLETE_ASSIGNMENTS,
    )

    print()

    results = []

    for sample_id, bits in enumerate(
        BITS_LIST,
        start=1,
    ):

        results.append(
            run_sample(
                bits,
                sample_id,
                rng,
            )
        )

    print()
    print("=" * 120)
    print(
        "EXPERIMENT 174 SUMMARY"
    )
    print("=" * 120)

    print()
    print(
        "bits | widths | removed | "
        "nodes | backtracks | complete | "
        "CRT | exact | recovered | time"
    )

    print("-" * 120)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{str(r['domain_widths']):>35} | "
            f"{r['values_removed']:>7} | "
            f"{r['nodes']:>5} | "
            f"{r['backtracks']:>10} | "
            f"{r['complete_assignments']:>8} | "
            f"{r['crt_tests']:>3} | "
            f"{r['exact_hits']:>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['total_time']:.6f}s"
        )

    print()
    print(
        "FINISHED EXPERIMENT 174"
    )


if __name__ == "__main__":
    main()
