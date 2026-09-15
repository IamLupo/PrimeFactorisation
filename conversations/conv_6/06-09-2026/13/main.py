#!/usr/bin/env python3

"""
START EXPERIMENT 171

2x2 RECTANGLE CONSISTENCY
+ RAREST RECTANGLE SEEDING
+ ADAPTIVE CRT PARAMETER SEARCH

Background
----------

Experiments 162-169 established:

    C(i,j) = F(n, r_i, s_j, a_i*b_j)

and therefore every cell relation is a constraint on:

    z_ij = a_i*b_j

with the global multiplicative rank-1 condition:

    z_ij*z_kl = z_il*z_kj.

Experiment 170 added A-A and B-B pairwise consistency.

However, Experiment 170 showed:

    AC3 survivors == pairwise survivors

for every tested sample.

Therefore ordinary pairwise consistency is not seeing the
correlation between FOUR cells simultaneously.

Experiment 171 introduces 2x2 rectangle consistency.

For rows i,k and columns j,l:

    (a_i,b_j)
    (a_i,b_l)
    (a_k,b_j)
    (a_k,b_l)

must all be valid simultaneously.

For a pivot pair (a_i,b_j), we require some b_l and some a_k
such that all four corresponding cell relations exist.

The same test is performed with every cell of a selected
2x2 rectangle as pivot.

We:

    1. Build the exact C-cell relations.
    2. Rank all 2x2 rectangles by relation size.
    3. Select the rarest rectangles.
    4. Iteratively remove cell pairs which have no compatible
       completion inside those rectangles.
    5. Select the smallest surviving cell relation.
    6. Anchor on its surviving pairs.
    7. Run ordinary AC-3 as a final local propagation step.
    8. Run the incremental CRT-parameter search.
    9. Verify n = p*q and the complete C-grid.

The main measurements are:

    original relation sizes
    rectangle seed size
    number removed by rectangle consistency
    true pair survival
    final anchor count
    CRT states
    exact recovery time

FINISHED EXPERIMENT 171
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

SEED = 1712026

# Number of rarest 2x2 rectangles to use.
RECTANGLE_LIMIT = 12

MAX_X_STATES = 2_000_000
MAX_T_STATES = 2_000_000


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


def rebuild_relation(
    old: Relation,
    pairs: Set[Tuple[int, int]],
) -> Relation:

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}
    products: Set[int] = set()

    for a, b in pairs:

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
        r=old.r,
        s=old.s,
        target_c=old.target_c,
        pairs=pairs,
        by_a=by_a,
        by_b=by_b,
        products=products,
    )


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

    result = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            result[(i, j)] = (
                build_relation(
                    n,
                    r,
                    s,
                    c_grid[i][j],
                )
            )

    return result


# ============================================================
# ORDINARY AC-3
# ============================================================

def ac3(
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[Tuple[int, int], Relation],
) -> bool:

    changed = True

    while changed:

        changed = False

        for i in range(len(R1)):

            for j in range(len(R2)):

                rel = relations[(i, j)]

                # A -> B
                old_A = A[i]

                new_A = {
                    a
                    for a in old_A
                    if rel.by_a.get(
                        a,
                        set(),
                    ) & B[j]
                }

                if new_A != old_A:

                    A[i] = new_A
                    changed = True

                    if not new_A:
                        return False

                # B -> A
                old_B = B[j]

                new_B = {
                    b
                    for b in old_B
                    if rel.by_b.get(
                        b,
                        set(),
                    ) & A[i]
                }

                if new_B != old_B:

                    B[j] = new_B
                    changed = True

                    if not new_B:
                        return False

    return True


# ============================================================
# RECTANGLE DEFINITIONS
# ============================================================

Rectangle = Tuple[int, int, int, int]


def all_rectangles() -> List[Rectangle]:

    result = []

    for i in range(len(R1)):

        for k in range(
            i + 1,
            len(R1),
        ):

            for j in range(len(R2)):

                for l in range(
                    j + 1,
                    len(R2),
                ):

                    result.append(
                        (
                            i,
                            k,
                            j,
                            l,
                        )
                    )

    return result


# ============================================================
# RECTANGLE SCORE
# ============================================================

def rectangle_score(
    rect: Rectangle,
    relations: Dict[Tuple[int, int], Relation],
) -> int:

    i, k, j, l = rect

    return (
        len(relations[(i, j)].pairs)
        * len(relations[(i, l)].pairs)
        * len(relations[(k, j)].pairs)
        * len(relations[(k, l)].pairs)
    )


# ============================================================
# PIVOT RECTANGLE SUPPORT
# ============================================================

def pivot_supported(
    pivot_i: int,
    pivot_j: int,
    a: int,
    b: int,
    other_i: int,
    other_j: int,
    relations: Dict[Tuple[int, int], Relation],
) -> bool:

    """
    Test whether:

        (a,b)

    in cell (pivot_i,pivot_j) can be completed to a valid
    2x2 rectangle using:

        other_i
        other_j

    We need some d such that:

        (a,d) in (pivot_i,other_j)

    and some c such that:

        (c,b) in (other_i,pivot_j)
        (c,d) in (other_i,other_j)

    """

    top_other = relations[
        (pivot_i, other_j)
    ]

    left_other = relations[
        (other_i, pivot_j)
    ]

    opposite = relations[
        (other_i, other_j)
    ]

    possible_d = top_other.by_a.get(
        a,
        set(),
    )

    if not possible_d:
        return False

    for d in possible_d:

        possible_c = (
            left_other.by_b.get(
                b,
                set(),
            )
            &
            opposite.by_b.get(
                d,
                set(),
            )
        )

        if possible_c:
            return True

    return False


# ============================================================
# ONE RECTANGLE PRUNING PASS
# ============================================================

def prune_rectangle(
    rect: Rectangle,
    relations: Dict[Tuple[int, int], Relation],
) -> Tuple[int, bool]:

    i, k, j, l = rect

    cells = [
        (i, j, k, l),
        (i, l, k, j),
        (k, j, i, l),
        (k, l, i, j),
    ]

    removed = 0
    any_changed = False

    for pi, pj, oi, oj in cells:

        rel = relations[(pi, pj)]

        if not rel.pairs:
            continue

        kept = set()

        for a, b in rel.pairs:

            if pivot_supported(
                pi,
                pj,
                a,
                b,
                oi,
                oj,
                relations,
            ):

                kept.add(
                    (a, b)
                )

        if len(kept) != len(
            rel.pairs
        ):

            removed += (
                len(rel.pairs)
                - len(kept)
            )

            relations[(pi, pj)] = (
                rebuild_relation(
                    rel,
                    kept,
                )
            )

            any_changed = True

    return removed, any_changed


# ============================================================
# MULTI-RECTANGLE CLOSURE
# ============================================================

@dataclass
class RectangleClosureStats:

    passes: int = 0
    removed: int = 0
    rectangles_processed: int = 0


def rectangle_closure(
    relations: Dict[Tuple[int, int], Relation],
    selected_rectangles: List[Rectangle],
) -> RectangleClosureStats:

    stats = RectangleClosureStats()

    while True:

        stats.passes += 1

        changed = False

        for rect in selected_rectangles:

            stats.rectangles_processed += 1

            removed, did_change = (
                prune_rectangle(
                    rect,
                    relations,
                )
            )

            stats.removed += removed

            if did_change:
                changed = True

            # Any empty relation means contradiction.
            i, k, j, l = rect

            if (
                not relations[(i, j)].pairs
                or
                not relations[(i, l)].pairs
                or
                not relations[(k, j)].pairs
                or
                not relations[(k, l)].pairs
            ):

                return stats

        if not changed:
            break

    return stats


# ============================================================
# RAREST CELL
# ============================================================

def find_rarest_cell(
    relations: Dict[Tuple[int, int], Relation],
) -> Tuple[int, int, Relation]:

    candidates = []

    for (i, j), rel in relations.items():

        candidates.append(
            (
                len(rel.pairs),
                len(rel.products),
                i,
                j,
                rel,
            )
        )

    candidates.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    _, _, i, j, rel = candidates[0]

    return i, j, rel


# ============================================================
# ANCHOR AC-3
# ============================================================

def anchor_propagation(
    anchor_i: int,
    anchor_j: int,
    a0: int,
    b0: int,
    relations: Dict[Tuple[int, int], Relation],
) -> Tuple[
    bool,
    List[Set[int]],
    List[Set[int]],
]:

    A = [
        set(range(r))
        for r in R1
    ]

    B = [
        set(range(s))
        for s in R2
    ]

    A[anchor_i] = {a0}
    B[anchor_j] = {b0}

    ok = ac3(
        A,
        B,
        relations,
    )

    return ok, A, B


# ============================================================
# C-GRID VERIFICATION
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
    y_cost: int

    estimated: int


def domain_product(
    domains: List[Set[int]],
    indices: Tuple[int, ...],
) -> int:

    value = 1

    for i in indices:
        value *= len(domains[i])

    return value


def modulus_product(
    moduli: List[int],
    indices: Tuple[int, ...],
) -> int:

    value = 1

    for i in indices:
        value *= moduli[i]

    return value


def choose_plan(
    X: List[Set[int]],
    X_moduli: List[int],
    Y: List[Set[int]],
    Y_moduli: List[int],
    sqrt_n: int,
) -> Plan:

    best = None

    for xmask in range(
        1,
        1 << len(X_moduli),
    ):

        x_indices = tuple(
            i
            for i in range(
                len(X_moduli)
            )
            if xmask & (1 << i)
        )

        x_cost = domain_product(
            X,
            x_indices,
        )

        if x_cost > MAX_X_STATES:
            continue

        M = modulus_product(
            X_moduli,
            x_indices,
        )

        t_max = (
            sqrt_n // M
        )

        for ymask in range(
            1,
            1 << len(Y_moduli),
        ):

            y_indices = tuple(
                j
                for j in range(
                    len(Y_moduli)
                )
                if ymask & (1 << j)
            )

            T = modulus_product(
                Y_moduli,
                y_indices,
            )

            if T <= t_max:
                continue

            y_cost = domain_product(
                Y,
                y_indices,
            )

            if y_cost > MAX_T_STATES:
                continue

            estimated = (
                x_cost
                * y_cost
            )

            candidate = Plan(
                x_indices=x_indices,
                y_indices=y_indices,
                M=M,
                t_modulus=T,
                x_cost=x_cost,
                y_cost=y_cost,
                estimated=estimated,
            )

            if (
                best is None
                or estimated
                < best.estimated
            ):

                best = candidate

    if best is None:

        raise RuntimeError(
            "No CRT plan"
        )

    return best


# ============================================================
# T RESIDUES
# ============================================================

def t_residue_set(
    n: int,
    P: int,
    M: int,
    s: int,
    B_domain: Set[int],
) -> Set[int]:

    result = set()

    n_mod = n % s
    P_mod = P % s

    M_inv = pow(
        M % s,
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
# INCREMENTAL T SEARCH
# ============================================================

@dataclass
class TStats:

    states_created: int = 0
    max_states: int = 0
    final_candidates: int = 0


def incremental_t_search(
    constraints: List[
        Tuple[int, Set[int]]
    ],
    t_max: int,
    stats: TStats,
) -> List[int]:

    ordered = sorted(
        constraints,
        key=lambda x: (
            len(x[1]),
            x[0],
        ),
    )

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

                if (
                    len(next_states)
                    > MAX_T_STATES
                ):

                    raise RuntimeError(
                        "T state limit exceeded"
                    )

        states = next_states

        stats.max_states = max(
            stats.max_states,
            len(states),
        )

        if not states:
            break

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
                    residue + k * modulus
                )

    stats.final_candidates += len(
        candidates
    )

    return sorted(candidates)


# ============================================================
# SIDE SEARCH
# ============================================================

@dataclass
class SideStats:

    plans: int = 0
    x_states: int = 0
    t_constraints: int = 0
    t_states: int = 0
    candidates: int = 0
    divisions: int = 0
    c_tests: int = 0


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

        return None, None

    stats.plans += 1

    x_lists = [
        sorted(X[i])
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

            allowed = t_residue_set(
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

            candidates = (
                incremental_t_search(
                    constraints,
                    t_max,
                    tstats,
                )
            )

        except RuntimeError:

            continue

        stats.t_states += (
            tstats.states_created
        )

        stats.candidates += (
            len(candidates)
        )

        for t in candidates:

            x = P + M * t

            if (
                x <= 1
                or x > sqrt_n
            ):
                continue

            stats.divisions += 1

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

            if verify_factor(
                n,
                pp,
                qq,
                c_grid,
            ):

                return pp, qq

    return None, None


# ============================================================
# SEMIPRIME
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

    return p, q, p * q


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

    original_sizes = {
        cell: len(rel.pairs)
        for cell, rel in relations.items()
    }

    # --------------------------------------------------------
    # Rank rectangles.
    # --------------------------------------------------------

    rectangles = all_rectangles()

    rectangles.sort(
        key=lambda rect: (
            rectangle_score(
                rect,
                relations,
            ),
            rect,
        )
    )

    selected = rectangles[
        :RECTANGLE_LIMIT
    ]

    # --------------------------------------------------------
    # Measure true pairs before pruning.
    # --------------------------------------------------------

    true_cells = {
        (i, j):
        (
            true_A[i],
            true_B[j],
        )
        for i in range(len(R1))
        for j in range(len(R2))
    }

    # --------------------------------------------------------
    # Rectangle closure.
    # --------------------------------------------------------

    closure_start = time.perf_counter()

    closure_stats = (
        rectangle_closure(
            relations,
            selected,
        )
    )

    closure_time = (
        time.perf_counter()
        - closure_start
    )

    # --------------------------------------------------------
    # Find rarest surviving cell.
    # --------------------------------------------------------

    anchor_i, anchor_j, anchor_rel = (
        find_rarest_cell(
            relations
        )
    )

    true_anchor = true_cells[
        (anchor_i, anchor_j)
    ]

    true_anchor_present = (
        true_anchor
        in anchor_rel.pairs
    )

    true_anchor_rank = None

    for rank, pair in enumerate(
        sorted(anchor_rel.pairs),
        start=1,
    ):

        if pair == true_anchor:

            true_anchor_rank = rank
            break

    # --------------------------------------------------------
    # Count surviving relation pairs.
    # --------------------------------------------------------

    surviving_total = sum(
        len(rel.pairs)
        for rel in relations.values()
    )

    original_total = sum(
        original_sizes.values()
    )

    removed_total = (
        original_total
        - surviving_total
    )

    # --------------------------------------------------------
    # Search.
    # --------------------------------------------------------

    A_stats = SideStats()
    B_stats = SideStats()

    anchor_attempts = 0
    ac3_survivors = 0

    exact = None

    true_anchor_domains = None
    true_B_domains = None

    start = time.perf_counter()

    for a0, b0 in sorted(
        anchor_rel.pairs
    ):

        anchor_attempts += 1

        ok, A, B = (
            anchor_propagation(
                anchor_i,
                anchor_j,
                a0,
                b0,
                relations,
            )
        )

        if not ok:
            continue

        ac3_survivors += 1

        if (
            a0 == true_anchor[0]
            and b0 == true_anchor[1]
        ):

            true_anchor_domains = [
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

        try:

            plan_A = choose_plan(
                A,
                R1,
                B,
                R2,
                sqrt_n,
            )

        except RuntimeError:

            plan_A = None

        try:

            plan_B = choose_plan(
                B,
                R2,
                A,
                R1,
                sqrt_n,
            )

        except RuntimeError:

            plan_B = None

        if (
            plan_A is None
            and plan_B is None
        ):

            continue

        if (
            plan_B is None
            or (
                plan_A is not None
                and plan_A.estimated
                <= plan_B.estimated
            )
        ):

            first = "A"

        else:

            first = "B"

        if first == "A":

            pp, qq = search_side(
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

            qq, pp = search_side(
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
                first,
                a0,
                b0,
            )

            break

        # ----------------------------------------------------
        # Opposite direction.
        # ----------------------------------------------------

        if first == "A":

            qq, pp = search_side(
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

            pp, qq = search_side(
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
                (
                    "B"
                    if first == "A"
                    else "A"
                ),
                a0,
                b0,
            )

            break

    search_time = (
        time.perf_counter()
        - start
    )

    total_time = (
        closure_time
        + search_time
    )

    recovered = False

    if exact is not None:

        pp, qq, direction, a0, b0 = (
            exact
        )

        recovered = (
            pp * qq == n
            and {
                pp,
                qq,
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
    print("=" * 84)
    print(
        f"EXPERIMENT 171 SAMPLE {sample_id}"
    )
    print("=" * 84)

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
        "--- RECTANGLE SELECTION ---"
    )

    print(
        f"rectangles total              = "
        f"{len(rectangles)}"
    )

    print(
        f"rectangles selected           = "
        f"{len(selected)}"
    )

    for rank, rect in enumerate(
        selected[:5],
        start=1,
    ):

        i, k, j, l = rect

        print(
            f"  #{rank}: "
            f"rows=({i},{k}) "
            f"cols=({j},{l}) "
            f"score={rectangle_score(rect, {
                cell: relations[cell]
                for cell in relations
            })}"
        )

    print()
    print(
        "--- RECTANGLE CLOSURE ---"
    )

    print(
        f"original relation pairs       = "
        f"{original_total}"
    )

    print(
        f"surviving relation pairs      = "
        f"{surviving_total}"
    )

    print(
        f"pairs removed                 = "
        f"{removed_total}"
    )

    print(
        f"closure passes                = "
        f"{closure_stats.passes}"
    )

    print(
        f"closure time                  = "
        f"{closure_time:.6f}s"
    )

    print()
    print(
        f"rarest surviving cell         = "
        f"({anchor_i},{anchor_j}) "
        f"r={R1[anchor_i]} "
        f"s={R2[anchor_j]}"
    )

    print(
        f"target C                      = "
        f"{anchor_rel.target_c}"
    )

    print(
        f"surviving pairs               = "
        f"{len(anchor_rel.pairs)}"
    )

    print(
        f"surviving products            = "
        f"{len(anchor_rel.products)}"
    )

    print(
        f"true anchor                   = "
        f"{true_anchor}"
    )

    print(
        f"true anchor present           = "
        f"{true_anchor_present}"
    )

    print(
        f"true anchor rank              = "
        f"{true_anchor_rank}"
    )

    print()
    print(
        "--- SEARCH RESULT ---"
    )

    print(
        f"anchor attempts               = "
        f"{anchor_attempts}"
    )

    print(
        f"AC-3 survivors                = "
        f"{ac3_survivors}"
    )

    if true_anchor_domains is not None:

        print()
        print(
            "TRUE ANCHOR DOMAINS"
        )

        print(
            "A widths                       =",
            [
                len(x)
                for x in true_anchor_domains
            ],
        )

        print(
            "B widths                       =",
            [
                len(x)
                for x in true_B_domains
            ],
        )

        print(
            "true A still present           =",
            all(
                true_A[i]
                in true_anchor_domains[i]
                for i in range(len(R1))
            ),
        )

        print(
            "true B still present           =",
            all(
                true_B[j]
                in true_B_domains[j]
                for j in range(len(R2))
            ),
        )

    print()
    print(
        "--- A -> p ---"
    )

    print(
        f"plans                          = "
        f"{A_stats.plans}"
    )

    print(
        f"X states                       = "
        f"{A_stats.x_states}"
    )

    print(
        f"t constraints                  = "
        f"{A_stats.t_constraints}"
    )

    print(
        f"t states                       = "
        f"{A_stats.t_states}"
    )

    print(
        f"candidates                     = "
        f"{A_stats.candidates}"
    )

    print(
        f"divisibility tests             = "
        f"{A_stats.divisions}"
    )

    print(
        f"C-grid tests                   = "
        f"{A_stats.c_tests}"
    )

    print()
    print(
        "--- B -> q ---"
    )

    print(
        f"plans                          = "
        f"{B_stats.plans}"
    )

    print(
        f"X states                       = "
        f"{B_stats.x_states}"
    )

    print(
        f"t constraints                  = "
        f"{B_stats.t_constraints}"
    )

    print(
        f"t states                       = "
        f"{B_stats.t_states}"
    )

    print(
        f"candidates                     = "
        f"{B_stats.candidates}"
    )

    print(
        f"divisibility tests             = "
        f"{B_stats.divisions}"
    )

    print(
        f"C-grid tests                   = "
        f"{B_stats.c_tests}"
    )

    print()
    print(
        f"exact                          = "
        f"{exact}"
    )

    print(
        f"exact original factors        = "
        f"{recovered}"
    )

    print()
    print(
        f"search time                    = "
        f"{search_time:.6f}s"
    )

    print(
        f"total time                     = "
        f"{total_time:.6f}s"
    )

    print("=" * 84)

    return {
        "bits": bits,
        "original_pairs": original_total,
        "surviving_pairs": surviving_total,
        "removed_pairs": removed_total,
        "closure_passes": closure_stats.passes,
        "closure_time": closure_time,
        "rarest_pairs": len(anchor_rel.pairs),
        "true_anchor_present": true_anchor_present,
        "true_anchor_rank": true_anchor_rank,
        "anchor_attempts": anchor_attempts,
        "ac3_survivors": ac3_survivors,
        "A_x_states": A_stats.x_states,
        "A_t_states": A_stats.t_states,
        "B_x_states": B_stats.x_states,
        "B_t_states": B_stats.t_states,
        "exact": exact,
        "recovered": recovered,
        "total_time": total_time,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(
        SEED
    )

    print(
        "START EXPERIMENT 171"
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
        "RECTANGLE_LIMIT =",
        RECTANGLE_LIMIT,
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
    print("=" * 115)
    print(
        "EXPERIMENT 171 SUMMARY"
    )
    print("=" * 115)

    print()
    print(
        "bits | orig-pairs | kept | removed | "
        "rarest | attempts | AC3 | "
        "A-X | A-t | B-X | B-t | "
        "exact | recovered | time"
    )

    print("-" * 115)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['original_pairs']:>10} | "
            f"{r['surviving_pairs']:>4} | "
            f"{r['removed_pairs']:>7} | "
            f"{r['rarest_pairs']:>6} | "
            f"{r['anchor_attempts']:>8} | "
            f"{r['ac3_survivors']:>3} | "
            f"{r['A_x_states']:>4} | "
            f"{r['A_t_states']:>5} | "
            f"{r['B_x_states']:>4} | "
            f"{r['B_t_states']:>5} | "
            f"{str(r['exact'] is not None):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['total_time']:.6f}s"
        )

    print()
    print(
        "FINISHED EXPERIMENT 171"
    )


if __name__ == "__main__":
    main()
