#!/usr/bin/env python3

"""
START EXPERIMENT 170

PAIRWISE / PATH CONSISTENCY CLOSURE
+ ADAPTIVE CRT-PARAMETER SEARCH

Experiment 169 showed that:

    - rarest-cell anchoring works;
    - AC-3 preserves the true solution;
    - CRT-parameterization avoids the large Cartesian product;
    - but AC-3 often leaves many globally incompatible residue
      assignments alive.

The weakness is that ordinary AC-3 only checks:

    A_i value -> some B_j value
    B_j value -> some A_i value

independently for each edge.

Experiment 170 adds induced same-side compatibility.

For two A variables A_i and A_k, values (a_i,a_k) are compatible
iff for EVERY column j there exists at least one b_j such that:

    (a_i,b_j) is allowed in cell (i,j)
    (a_k,b_j) is allowed in cell (k,j)

Likewise, for B_j and B_l, values (b_j,b_l) are compatible iff for
EVERY row i there exists at least one a_i supporting both values.

This gives a path-consistency-style closure on the complete
5x5 bipartite C-grid.

After this stronger closure:

    1. choose the rarest C-cell;
    2. anchor each surviving pair;
    3. run pairwise consistency closure;
    4. choose the cheaper CRT direction;
    5. perform incremental CRT parameter search;
    6. verify n = p*q and the entire C-grid.

The experiment measures how much stronger pairwise consistency is
than ordinary AC-3.

FINISHED EXPERIMENT 170
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

BITS_LIST = [30, 36, 42, 48, 54, 60, 66, 72, 78, 84]

SEED = 1702026

MAX_X_STATES = 2_000_000
MAX_T_STATES = 2_000_000

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
# EXACT C
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
# RELATION
# ============================================================

@dataclass
class Relation:

    r: int
    s: int
    target_c: int

    pairs: Set[Tuple[int, int]]

    by_a: Dict[int, Set[int]]
    by_b: Dict[int, Set[int]]

    products: Set[int]


def build_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> Relation:

    pairs = set()

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}

    products = set()

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

            products.add(
                a * b
            )

    return Relation(
        r=r,
        s=s,
        target_c=target_c,
        pairs=pairs,
        by_a=by_a,
        by_b=by_b,
        products=products,
    )


def build_relations(
    n: int,
    c_grid: List[List[int]],
) -> Dict[Tuple[int, int], Relation]:

    relations = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            relations[(i, j)] = (
                build_relation(
                    n,
                    r,
                    s,
                    c_grid[i][j],
                )
            )

    return relations


# ============================================================
# ORDINARY AC-3
# ============================================================

def revise_A(
    i: int,
    j: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[Tuple[int, int], Relation],
) -> bool:

    rel = relations[(i, j)]

    old = A[i]

    new = {
        a
        for a in old
        if rel.by_a.get(
            a,
            set(),
        ) & B[j]
    }

    if new == old:
        return False

    A[i] = new
    return True


def revise_B(
    i: int,
    j: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[Tuple[int, int], Relation],
) -> bool:

    rel = relations[(i, j)]

    old = B[j]

    new = {
        b
        for b in old
        if rel.by_b.get(
            b,
            set(),
        ) & A[i]
    }

    if new == old:
        return False

    B[j] = new
    return True


def ac3_original(
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[Tuple[int, int], Relation],
) -> bool:

    queue = []

    for i in range(len(R1)):

        for j in range(len(R2)):

            queue.append(
                ("A", i, j)
            )

            queue.append(
                ("B", i, j)
            )

    while queue:

        kind, i, j = queue.pop()

        if kind == "A":

            changed = revise_A(
                i,
                j,
                A,
                B,
                relations,
            )

            if changed:

                if not A[i]:
                    return False

                for jj in range(len(R2)):

                    if jj != j:
                        queue.append(
                            ("B", i, jj)
                        )

        else:

            changed = revise_B(
                i,
                j,
                A,
                B,
                relations,
            )

            if changed:

                if not B[j]:
                    return False

                for ii in range(len(R1)):

                    if ii != i:
                        queue.append(
                            ("A", ii, j)
                        )

    return True


# ============================================================
# INDUCED A-A COMPATIBILITY
# ============================================================

def build_AA_compatibility(
    relations: Dict[Tuple[int, int], Relation],
) -> Dict[Tuple[int, int], Dict[int, Set[int]]]:

    """
    compat[(i,k)][a_i] = set of a_k values compatible with a_i.

    Two A values are compatible iff for EVERY column j they have
    at least one common B_j support.

    This is stronger than ordinary AC.
    """

    compat = {}

    nr = len(R1)
    ns = len(R2)

    for i in range(nr):

        for k in range(i + 1, nr):

            mapping: Dict[int, Set[int]] = {}

            ri_values = range(R1[i])
            rk_values = range(R1[k])

            for a in ri_values:

                allowed_k = set()

                for ak in rk_values:

                    ok = True

                    for j in range(ns):

                        support_i = (
                            relations[(i, j)]
                            .by_a
                            .get(a, set())
                        )

                        support_k = (
                            relations[(k, j)]
                            .by_a
                            .get(ak, set())
                        )

                        if not (
                            support_i
                            & support_k
                        ):

                            ok = False
                            break

                    if ok:
                        allowed_k.add(ak)

                mapping[a] = allowed_k

            compat[(i, k)] = mapping

    return compat


# ============================================================
# INDUCED B-B COMPATIBILITY
# ============================================================

def build_BB_compatibility(
    relations: Dict[Tuple[int, int], Relation],
) -> Dict[Tuple[int, int], Dict[int, Set[int]]]:

    """
    compat[(j,l)][b_j] = set of b_l values compatible with b_j.

    Two B values are compatible iff for EVERY row i they have
    at least one common A_i support.
    """

    compat = {}

    nr = len(R1)
    ns = len(R2)

    for j in range(ns):

        for l in range(j + 1, ns):

            mapping: Dict[int, Set[int]] = {}

            for b in range(R2[j]):

                allowed_l = set()

                for bl in range(R2[l]):

                    ok = True

                    for i in range(nr):

                        support_j = (
                            relations[(i, j)]
                            .by_b
                            .get(b, set())
                        )

                        support_l = (
                            relations[(i, l)]
                            .by_b
                            .get(bl, set())
                        )

                        if not (
                            support_j
                            & support_l
                        ):

                            ok = False
                            break

                    if ok:
                        allowed_l.add(bl)

                mapping[b] = allowed_l

            compat[(j, l)] = mapping

    return compat


# ============================================================
# PATH-CONSISTENCY CLOSURE
# ============================================================

def revise_AA(
    i: int,
    k: int,
    A: List[Set[int]],
    aa: Dict[Tuple[int, int], Dict[int, Set[int]]],
) -> bool:

    if i < k:

        mapping = aa[(i, k)]

        old = A[i]

        new = {
            a
            for a in old
            if mapping.get(
                a,
                set(),
            ) & A[k]
        }

    else:

        mapping = aa[(k, i)]

        old = A[i]

        new = {
            a
            for a in old
            if any(
                a in mapping.get(
                    ak,
                    set(),
                )
                for ak in A[k]
            )
        }

    if new == old:
        return False

    A[i] = new
    return True


def revise_BB(
    j: int,
    l: int,
    B: List[Set[int]],
    bb: Dict[Tuple[int, int], Dict[int, Set[int]]],
) -> bool:

    if j < l:

        mapping = bb[(j, l)]

        old = B[j]

        new = {
            b
            for b in old
            if mapping.get(
                b,
                set(),
            ) & B[l]
        }

    else:

        mapping = bb[(l, j)]

        old = B[j]

        new = {
            b
            for b in old
            if any(
                b in mapping.get(
                    bl,
                    set(),
                )
                for bl in B[l]
            )
        }

    if new == old:
        return False

    B[j] = new
    return True


def path_consistency(
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[Tuple[int, int], Relation],
    aa: Dict[Tuple[int, int], Dict[int, Set[int]]],
    bb: Dict[Tuple[int, int], Dict[int, Set[int]]],
) -> bool:

    """
    Alternate:

        original A-B arc consistency
        induced A-A consistency
        induced B-B consistency

    until stable.

    This is a domain-consistency approximation to path consistency.
    """

    changed = True

    while changed:

        changed = False

        # -----------------------------------------------
        # Original bipartite relations.
        # -----------------------------------------------

        if not ac3_original(
            A,
            B,
            relations,
        ):
            return False

        # -----------------------------------------------
        # A-A relations.
        # -----------------------------------------------

        for i in range(len(R1)):

            for k in range(
                i + 1,
                len(R1),
            ):

                if revise_AA(
                    i,
                    k,
                    A,
                    aa,
                ):

                    changed = True

                    if not A[i]:
                        return False

                if revise_AA(
                    k,
                    i,
                    A,
                    aa,
                ):

                    changed = True

                    if not A[k]:
                        return False

        # -----------------------------------------------
        # B-B relations.
        # -----------------------------------------------

        for j in range(len(R2)):

            for l in range(
                j + 1,
                len(R2),
            ):

                if revise_BB(
                    j,
                    l,
                    B,
                    bb,
                ):

                    changed = True

                    if not B[j]:
                        return False

                if revise_BB(
                    l,
                    j,
                    B,
                    bb,
                ):

                    changed = True

                    if not B[l]:
                        return False

    return True


# ============================================================
# ANCHOR PROPAGATION
# ============================================================

def propagate_anchor(
    anchor_i: int,
    anchor_j: int,
    a0: int,
    b0: int,
    relations: Dict[Tuple[int, int], Relation],
    aa: Dict[Tuple[int, int], Dict[int, Set[int]]],
    bb: Dict[Tuple[int, int], Dict[int, Set[int]]],
) -> Tuple[
    bool,
    List[Set[int]],
    List[Set[int]],
    List[Set[int]],
    List[Set[int]],
]:

    # -----------------------------------------------
    # Initial domains.
    # -----------------------------------------------

    A = [
        set(range(r))
        for r in R1
    ]

    B = [
        set(range(s))
        for s in R2
    ]

    # -----------------------------------------------
    # Anchor.
    # -----------------------------------------------

    A[anchor_i] = {a0}
    B[anchor_j] = {b0}

    # -----------------------------------------------
    # First ordinary AC-3 only.
    # -----------------------------------------------

    initial_A = [
        set(x)
        for x in A
    ]

    initial_B = [
        set(x)
        for x in B
    ]

    if not ac3_original(
        A,
        B,
        relations,
    ):

        return (
            False,
            initial_A,
            initial_B,
            A,
            B,
        )

    # -----------------------------------------------
    # Full pairwise closure.
    # -----------------------------------------------

    if not path_consistency(
        A,
        B,
        relations,
        aa,
        bb,
    ):

        return (
            False,
            initial_A,
            initial_B,
            A,
            B,
        )

    return (
        True,
        initial_A,
        initial_B,
        A,
        B,
    )


# ============================================================
# FULL FACTOR VERIFICATION
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
# CRT PLAN
# ============================================================

@dataclass
class Plan:

    x_indices: Tuple[int, ...]
    y_indices: Tuple[int, ...]

    M: int
    t_modulus: int

    x_cost: int
    t_cost: int

    estimated: int


def subset_product(
    domains: List[Set[int]],
    indices: Tuple[int, ...],
) -> int:

    result = 1

    for i in indices:

        result *= len(
            domains[i]
        )

    return result


def subset_modulus(
    moduli: List[int],
    indices: Tuple[int, ...],
) -> int:

    result = 1

    for i in indices:

        result *= moduli[i]

    return result


def choose_plan(
    X: List[Set[int]],
    X_moduli: List[int],
    Y: List[Set[int]],
    Y_moduli: List[int],
    sqrt_n: int,
) -> Plan:

    best = None

    nx = len(X_moduli)
    ny = len(Y_moduli)

    # -----------------------------------------------
    # Only subsets which make all remaining t unique
    # are useful.
    # -----------------------------------------------

    for xm in range(
        1,
        1 << nx,
    ):

        x_indices = tuple(
            i
            for i in range(nx)
            if xm & (1 << i)
        )

        x_cost = subset_product(
            X,
            x_indices,
        )

        if x_cost > MAX_X_STATES:
            continue

        M = subset_modulus(
            X_moduli,
            x_indices,
        )

        t_max = (
            sqrt_n // M
        )

        for ym in range(
            1,
            1 << ny,
        ):

            y_indices = tuple(
                j
                for j in range(ny)
                if ym & (1 << j)
            )

            t_modulus = subset_modulus(
                Y_moduli,
                y_indices,
            )

            if t_modulus <= t_max:
                continue

            t_cost = subset_product(
                Y,
                y_indices,
            )

            if t_cost > MAX_T_STATES:
                continue

            estimated = (
                x_cost
                * t_cost
            )

            candidate = Plan(
                x_indices=x_indices,
                y_indices=y_indices,
                M=M,
                t_modulus=t_modulus,
                x_cost=x_cost,
                t_cost=t_cost,
                estimated=estimated,
            )

            if (
                best is None
                or candidate.estimated
                < best.estimated
            ):

                best = candidate

    if best is None:
        raise RuntimeError(
            "No feasible CRT plan"
        )

    return best


# ============================================================
# t RESIDUES
# ============================================================

def t_residues(
    n: int,
    P: int,
    M: int,
    s: int,
    B_domain: Set[int],
) -> Set[int]:

    result = set()

    n_mod = n % s
    P_mod = P % s
    M_mod = M % s

    M_inv = pow(
        M_mod,
        -1,
        s,
    )

    for b in B_domain:

        b_mod = b % s

        if b_mod == 0:

            if n_mod == 0:

                result.update(
                    range(s)
                )

            continue

        b_inv = pow(
            b_mod,
            -1,
            s,
        )

        t = (
            (
                n_mod * b_inv
                - P_mod
            )
            * M_inv
        ) % s

        result.add(t)

    return result


# ============================================================
# INCREMENTAL T CRT
# ============================================================

@dataclass
class TStats:

    states_created: int = 0

    levels: List[int] = None

    final_candidates: int = 0

    max_level_states: int = 0

    def __post_init__(self):
        if self.levels is None:
            self.levels = []


def incremental_t_crt(
    constraints: List[Tuple[int, Set[int]]],
    t_max: int,
    stats: TStats,
) -> List[int]:

    if not constraints:

        return list(
            range(
                t_max + 1
            )
        )

    # Strongest constraint first.
    ordered = sorted(
        constraints,
        key=lambda x: (
            len(x[1]),
            x[0],
        ),
    )

    # (residue, modulus)
    states = {
        (0, 1)
    }

    for modulus, allowed in ordered:

        next_states = set()

        for current, current_M in states:

            for residue in allowed:

                new_t, new_M = crt_pair(
                    current,
                    current_M,
                    residue,
                    modulus,
                )

                stats.states_created += 1

                if new_M > t_max:

                    if new_t <= t_max:

                        next_states.add(
                            (
                                new_t,
                                new_M,
                            )
                        )

                else:

                    next_states.add(
                        (
                            new_t,
                            new_M,
                        )
                    )

                if len(next_states) > MAX_T_STATES:

                    raise RuntimeError(
                        "T state limit exceeded"
                    )

        states = next_states

        size = len(states)

        stats.levels.append(size)

        stats.max_level_states = max(
            stats.max_level_states,
            size,
        )

        if not states:
            return []

    candidates = set()

    for residue, modulus in states:

        if modulus > t_max:

            if residue <= t_max:
                candidates.add(
                    residue
                )

        else:

            kmax = (
                t_max - residue
            ) // modulus

            for k in range(
                kmax + 1
            ):

                candidates.add(
                    residue
                    + k * modulus
                )

    result = sorted(
        candidates
    )

    stats.final_candidates += len(
        result
    )

    return result


# ============================================================
# SIDE SEARCH
# ============================================================

@dataclass
class SideStats:

    plans: int = 0

    x_states: int = 0

    t_constraints: int = 0

    t_states_created: int = 0

    t_final_candidates: int = 0

    divisibility_tests: int = 0

    c_tests: int = 0

    max_t_level_states: int = 0


def search_side(
    *,
    n: int,
    sqrt_n: int,
    X: List[Set[int]],
    X_moduli: List[int],
    Y: List[Set[int]],
    Y_moduli: List[int],
    c_grid: List[List[int]],
    X_is_A: bool,
    stats: SideStats,
) -> Tuple[
    int | None,
    int | None,
    Plan | None,
]:

    try:

        plan = choose_plan(
            X,
            X_moduli,
            Y,
            Y_moduli,
            sqrt_n,
        )

    except RuntimeError:

        return None, None, None

    stats.plans += 1

    x_lists = [
        sorted(
            X[i]
        )
        for i in plan.x_indices
    ]

    for x_values in itertools.product(
        *x_lists
    ):

        stats.x_states += 1

        if (
            stats.x_states
            > MAX_X_STATES
        ):
            raise RuntimeError(
                "X state limit exceeded"
            )

        P, M = crt_many(
            list(x_values),
            [
                X_moduli[i]
                for i in plan.x_indices
            ],
        )

        if P > sqrt_n:
            continue

        t_max = (
            sqrt_n - P
        ) // M

        if t_max < 0:
            continue

        constraints = []

        impossible = False

        for j in plan.y_indices:

            allowed = t_residues(
                n,
                P,
                M,
                Y_moduli[j],
                Y[j],
            )

            stats.t_constraints += 1

            if not allowed:

                impossible = True
                break

            constraints.append(
                (
                    Y_moduli[j],
                    allowed,
                )
            )

        if impossible:
            continue

        tstats = TStats()

        try:

            candidates = incremental_t_crt(
                constraints,
                t_max,
                tstats,
            )

        except RuntimeError:

            continue

        stats.t_states_created += (
            tstats.states_created
        )

        stats.t_final_candidates += (
            len(candidates)
        )

        stats.max_t_level_states = max(
            stats.max_t_level_states,
            tstats.max_level_states,
        )

        for t in candidates:

            x = P + M * t

            if (
                x <= 1
                or x > sqrt_n
            ):
                continue

            stats.divisibility_tests += 1

            if n % x != 0:
                continue

            y = n // x

            stats.c_tests += 1

            if X_is_A:

                pp = x
                qq = y

            else:

                pp = y
                qq = x

            if not verify_factor(
                n,
                pp,
                qq,
                c_grid,
            ):
                continue

            return pp, qq, plan

    return None, None, plan


# ============================================================
# SEMIPRIME GENERATION
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

    sqrt_n = math.isqrt(n)

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

    # --------------------------------------------------------
    # Strongest cell.
    # --------------------------------------------------------

    ranked = []

    for (i, j), rel in relations.items():

        ranked.append(
            (
                len(rel.pairs),
                len(rel.products),
                i,
                j,
                rel,
            )
        )

    ranked.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    (
        rare_pairs,
        rare_products,
        anchor_i,
        anchor_j,
        anchor_rel,
    ) = ranked[0]

    true_anchor = (
        true_A[anchor_i],
        true_B[anchor_j],
    )

    true_rank = None

    for rank, pair in enumerate(
        sorted(anchor_rel.pairs),
        start=1,
    ):

        if pair == true_anchor:

            true_rank = rank
            break

    # --------------------------------------------------------
    # Build induced pairwise relations.
    # --------------------------------------------------------

    closure_start = time.perf_counter()

    aa = build_AA_compatibility(
        relations
    )

    bb = build_BB_compatibility(
        relations
    )

    closure_build_time = (
        time.perf_counter()
        - closure_start
    )

    print()
    print("=" * 82)
    print(
        f"EXPERIMENT 170 SAMPLE {sample_id}"
    )
    print("=" * 82)

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
        f"sqrt(n)           = {sqrt_n}"
    )

    print()
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
        f"rarest cell       = "
        f"({anchor_i},{anchor_j}) "
        f"r={R1[anchor_i]} "
        f"s={R2[anchor_j]}"
    )

    print(
        f"target C          = "
        f"{anchor_rel.target_c}"
    )

    print(
        f"allowed pairs     = "
        f"{rare_pairs}"
    )

    print(
        f"allowed products  = "
        f"{rare_products}"
    )

    print(
        f"true anchor pair  = "
        f"{true_anchor}"
    )

    print(
        f"true anchor rank  = "
        f"{true_rank}"
    )

    print()
    print(
        f"pairwise build time = "
        f"{closure_build_time:.6f}s"
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    start = time.perf_counter()

    anchor_attempts = 0

    ordinary_ac3_survivors = 0
    pairwise_survivors = 0

    true_anchor_survived = False

    true_A_domains = None
    true_B_domains = None

    A_stats = SideStats()
    B_stats = SideStats()

    exact = None

    for a0, b0 in sorted(
        anchor_rel.pairs
    ):

        anchor_attempts += 1

        (
            ok,
            before_A,
            before_B,
            A,
            B,
        ) = propagate_anchor(
            anchor_i,
            anchor_j,
            a0,
            b0,
            relations,
            aa,
            bb,
        )

        if not ok:
            continue

        ordinary_ac3_survivors += 1

        # ----------------------------------------------------
        # Check whether path closure actually survives.
        # ----------------------------------------------------

        # propagate_anchor already executes path closure,
        # therefore reaching this point means pairwise survived.
        pairwise_survivors += 1

        if (
            a0 == true_anchor[0]
            and b0 == true_anchor[1]
        ):

            true_anchor_survived = True

            true_A_domains = [
                set(x)
                for x in A
            ]

            true_B_domains = [
                set(x)
                for x in B
            ]

        # ----------------------------------------------------
        # Estimate both directions.
        # ----------------------------------------------------

        plan_A = None
        plan_B = None

        try:

            plan_A = choose_plan(
                A,
                R1,
                B,
                R2,
                sqrt_n,
            )

        except RuntimeError:
            pass

        try:

            plan_B = choose_plan(
                B,
                R2,
                A,
                R1,
                sqrt_n,
            )

        except RuntimeError:
            pass

        if (
            plan_A is None
            and plan_B is None
        ):
            continue

        if (
            plan_B is None
            or (
                plan_A is not None
                and (
                    plan_A.estimated
                    <= plan_B.estimated
                )
            )
        ):

            first = "A"

        else:

            first = "B"

        # ----------------------------------------------------
        # Search first direction.
        # ----------------------------------------------------

        if first == "A":

            pp, qq, _ = search_side(
                n=n,
                sqrt_n=sqrt_n,
                X=A,
                X_moduli=R1,
                Y=B,
                Y_moduli=R2,
                c_grid=c_grid,
                X_is_A=True,
                stats=A_stats,
            )

        else:

            qq, pp, _ = search_side(
                n=n,
                sqrt_n=sqrt_n,
                X=B,
                X_moduli=R2,
                Y=A,
                Y_moduli=R1,
                c_grid=c_grid,
                X_is_A=False,
                stats=B_stats,
            )

        if pp is not None:

            exact = (
                pp,
                qq,
                f"{first}->factor",
                a0,
                b0,
            )

            break

        # ----------------------------------------------------
        # Opposite direction.
        # ----------------------------------------------------

        if first == "A":

            qq, pp, _ = search_side(
                n=n,
                sqrt_n=sqrt_n,
                X=B,
                X_moduli=R2,
                Y=A,
                Y_moduli=R1,
                c_grid=c_grid,
                X_is_A=False,
                stats=B_stats,
            )

        else:

            pp, qq, _ = search_side(
                n=n,
                sqrt_n=sqrt_n,
                X=A,
                X_moduli=R1,
                Y=B,
                Y_moduli=R2,
                c_grid=c_grid,
                X_is_A=True,
                stats=A_stats,
            )

        if pp is not None:

            exact = (
                pp,
                qq,
                f"{'B' if first == 'A' else 'A'}->factor",
                a0,
                b0,
            )

            break

    elapsed = (
        time.perf_counter()
        - start
    )

    # --------------------------------------------------------
    # Output.
    # --------------------------------------------------------

    print()
    print(
        "--- SEARCH RESULT ---"
    )

    print(
        f"anchor attempts             = "
        f"{anchor_attempts}"
    )

    print(
        f"AC-3 survivors              = "
        f"{ordinary_ac3_survivors}"
    )

    print(
        f"pairwise survivors          = "
        f"{pairwise_survivors}"
    )

    print(
        f"true anchor survived        = "
        f"{true_anchor_survived}"
    )

    if true_A_domains is not None:

        print()
        print(
            "TRUE ANCHOR DOMAINS"
        )

        print(
            "A widths                     =",
            [
                len(x)
                for x in true_A_domains
            ],
        )

        print(
            "B widths                     =",
            [
                len(x)
                for x in true_B_domains
            ],
        )

        print(
            "A domains                    =",
            [
                sorted(x)
                for x in true_A_domains
            ],
        )

        print(
            "B domains                    =",
            [
                sorted(x)
                for x in true_B_domains
            ],
        )

        print(
            "true A still present         =",
            all(
                true_A[i]
                in true_A_domains[i]
                for i in range(len(R1))
            ),
        )

        print(
            "true B still present         =",
            all(
                true_B[j]
                in true_B_domains[j]
                for j in range(len(R2))
            ),
        )

    print()
    print(
        "--- A -> p STATISTICS ---"
    )

    print(
        f"plans                        = "
        f"{A_stats.plans}"
    )

    print(
        f"X CRT states                 = "
        f"{A_stats.x_states}"
    )

    print(
        f"t constraint calls           = "
        f"{A_stats.t_constraints}"
    )

    print(
        f"t CRT states created         = "
        f"{A_stats.t_states_created}"
    )

    print(
        f"final t candidates           = "
        f"{A_stats.t_final_candidates}"
    )

    print(
        f"divisibility tests           = "
        f"{A_stats.divisibility_tests}"
    )

    print(
        f"C-grid tests                 = "
        f"{A_stats.c_tests}"
    )

    print(
        f"max t-level states           = "
        f"{A_stats.max_t_level_states}"
    )

    print()
    print(
        "--- B -> q STATISTICS ---"
    )

    print(
        f"plans                        = "
        f"{B_stats.plans}"
    )

    print(
        f"X CRT states                 = "
        f"{B_stats.x_states}"
    )

    print(
        f"t constraint calls           = "
        f"{B_stats.t_constraints}"
    )

    print(
        f"t CRT states created         = "
        f"{B_stats.t_states_created}"
    )

    print(
        f"final t candidates           = "
        f"{B_stats.t_final_candidates}"
    )

    print(
        f"divisibility tests           = "
        f"{B_stats.divisibility_tests}"
    )

    print(
        f"C-grid tests                 = "
        f"{B_stats.c_tests}"
    )

    print(
        f"max t-level states           = "
        f"{B_stats.max_t_level_states}"
    )

    print()
    print(
        f"exact                         = "
        f"{exact}"
    )

    recovered = False

    if exact is not None:

        fp, fq, direction, aa0, bb0 = exact

        recovered = (
            fp * fq == n
            and {
                fp,
                fq,
            }
            == {
                p,
                q,
            }
        )

        print(
            f"exact original factors      = "
            f"{recovered}"
        )

    print()
    print(
        f"runtime                      = "
        f"{elapsed:.6f}s"
    )

    print("=" * 82)

    return {
        "bits": bits,
        "rarest_pairs": rare_pairs,
        "rarest_products": rare_products,
        "true_anchor_rank": true_rank,
        "anchor_attempts": anchor_attempts,
        "ac3_survivors": ordinary_ac3_survivors,
        "pairwise_survivors": pairwise_survivors,
        "true_anchor_survived": true_anchor_survived,
        "A_x_states": A_stats.x_states,
        "A_t_states": A_stats.t_states_created,
        "A_candidates": A_stats.t_final_candidates,
        "B_x_states": B_stats.x_states,
        "B_t_states": B_stats.t_states_created,
        "B_candidates": B_stats.t_final_candidates,
        "exact": exact,
        "recovered": recovered,
        "runtime": elapsed,
        "pairwise_build_time": closure_build_time,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(
        SEED
    )

    print(
        "START EXPERIMENT 170"
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

    print()

    results = []

    for sample_id, bits in enumerate(
        BITS_LIST,
        start=1,
    ):

        result = run_sample(
            bits,
            sample_id,
            rng,
        )

        results.append(
            result
        )

    print()
    print("=" * 115)
    print(
        "EXPERIMENT 170 SUMMARY"
    )
    print("=" * 115)

    print()
    print(
        "bits | rare | attempts | AC3 | "
        "PCons | A-X | A-t | B-X | B-t | "
        "exact | recovered | time"
    )

    print("-" * 115)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['rarest_pairs']:>4} | "
            f"{r['anchor_attempts']:>8} | "
            f"{r['ac3_survivors']:>3} | "
            f"{r['pairwise_survivors']:>5} | "
            f"{r['A_x_states']:>4} | "
            f"{r['A_t_states']:>5} | "
            f"{r['B_x_states']:>4} | "
            f"{r['B_t_states']:>5} | "
            f"{str(r['exact'] is not None):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['runtime']:.6f}s"
        )

    print()
    print(
        "FINISHED EXPERIMENT 170"
    )


if __name__ == "__main__":
    main()
