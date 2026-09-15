#!/usr/bin/env python3

"""
START EXPERIMENT 172

FULL 2x2 RECTANGLE CONSISTENCY
USING INTEGER BITSETS

Experiment 171:
    - selected only 12 of the 100 possible 2x2 rectangles;
    - removed many relation pairs;
    - but the remaining domains were still broad.

Experiment 172:
    - uses ALL 100 possible 2x2 rectangles;
    - iterates rectangle closure until stable;
    - accelerates support tests with Python integer bitsets.

For a cell relation R_ij we have pairs (a,b).

For a 2x2 rectangle:

        (i,j)   (i,l)
        (k,j)   (k,l)

a pair (a,b) in R_ij survives iff there exist c,d such that:

        (a,d) in R_i,l
        (c,b) in R_k,j
        (c,d) in R_k,l

This is exactly the 2x2 multiplicative consistency condition.

Bitset representation:
    by_a[a] = bitset of compatible b values
    by_b[b] = bitset of compatible a values

For the opposite relation, we precompute:

    neighbor_union[d] = bitset of c values adjacent to d

Then a pivot pair survives iff:

    OR over d in top_support(a)
        opposite.by_b[d]

intersects left_support(b).

The OR operation is performed using Python integers.

After complete rectangle closure:
    - choose the rarest surviving cell;
    - anchor each surviving pair;
    - run ordinary AC-3;
    - use adaptive incremental CRT search;
    - verify exact factorization.

The primary scientific question:

    Does FULL 2x2 consistency materially reduce the
    remaining residue domains beyond Experiment 171?

FINISHED EXPERIMENT 172
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

SEED = 1722026

MAX_X_STATES = 2_000_000
MAX_T_STATES = 2_000_000

MAX_ANCHOR_ATTEMPTS = 2_000


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

    # Integer bitsets.
    bits_by_a: Dict[int, int]
    bits_by_b: Dict[int, int]

    # opposite_by_b[b] gives a-bitset.
    products: Set[int]


def build_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> Relation:

    pairs: Set[Tuple[int, int]] = set()

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}

    bits_by_a: Dict[int, int] = {}
    bits_by_b: Dict[int, int] = {}

    products: Set[int] = set()

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

            bits_by_a[a] = (
                bits_by_a.get(a, 0)
                | (1 << b)
            )

            bits_by_b[b] = (
                bits_by_b.get(b, 0)
                | (1 << a)
            )

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
        bits_by_a=bits_by_a,
        bits_by_b=bits_by_b,
        products=products,
    )


def rebuild_relation(
    old: Relation,
    pairs: Set[Tuple[int, int]],
) -> Relation:

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}

    bits_by_a: Dict[int, int] = {}
    bits_by_b: Dict[int, int] = {}

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

        bits_by_a[a] = (
            bits_by_a.get(a, 0)
            | (1 << b)
        )

        bits_by_b[b] = (
            bits_by_b.get(b, 0)
            | (1 << a)
        )

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
        bits_by_a=bits_by_a,
        bits_by_b=bits_by_b,
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
# RECTANGLES
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
# BIT ITERATION
# ============================================================

def iter_bits(
    x: int,
):

    while x:

        low = (
            x
            & -x
        )

        bit = (
            low.bit_length()
            - 1
        )

        yield bit

        x -= low


# ============================================================
# RECTANGLE SUPPORT
# ============================================================

def pair_has_rectangle_support(
    pivot_a: int,
    pivot_b: int,
    top: Relation,
    left: Relation,
    opposite: Relation,
) -> bool:

    """
    Pivot relation:

        pivot = (i,j)

    top:
        (i,l)

    left:
        (k,j)

    opposite:
        (k,l)

    Need:

        exists d:
            (pivot_a,d) in top

        and

        exists c:
            (c,pivot_b) in left
            (c,d) in opposite
    """

    top_bits = top.bits_by_a.get(
        pivot_a,
        0,
    )

    if not top_bits:
        return False

    left_bits = left.bits_by_b.get(
        pivot_b,
        0,
    )

    if not left_bits:
        return False

    # Walk over whichever side has fewer candidates.
    top_count = top_bits.bit_count()
    left_count = left_bits.bit_count()

    if top_count <= left_count:

        for d in iter_bits(
            top_bits
        ):

            opposite_c = (
                opposite.bits_by_b.get(
                    d,
                    0,
                )
            )

            if (
                opposite_c
                & left_bits
            ):
                return True

    else:

        # Build the opposite-side support
        # for all d supported by pivot_a.
        union_c = 0

        bits = top_bits

        while bits:

            low = (
                bits
                & -bits
            )

            d = (
                low.bit_length()
                - 1
            )

            union_c |= (
                opposite.bits_by_b.get(
                    d,
                    0,
                )
            )

            bits -= low

        if union_c & left_bits:
            return True

    return False


# ============================================================
# ONE RECTANGLE PASS
# ============================================================

def prune_rectangle(
    rect: Rectangle,
    relations: Dict[Tuple[int, int], Relation],
) -> Tuple[int, bool]:

    i, k, j, l = rect

    orientations = [
        (
            (i, j),
            (i, l),
            (k, j),
            (k, l),
        ),
        (
            (i, l),
            (i, j),
            (k, l),
            (k, j),
        ),
        (
            (k, j),
            (k, l),
            (i, j),
            (i, l),
        ),
        (
            (k, l),
            (k, j),
            (i, l),
            (i, j),
        ),
    ]

    removed = 0
    changed = False

    for pivot_cell, top_cell, left_cell, opposite_cell in orientations:

        pivot = relations[
            pivot_cell
        ]

        if not pivot.pairs:
            continue

        top = relations[
            top_cell
        ]

        left = relations[
            left_cell
        ]

        opposite = relations[
            opposite_cell
        ]

        kept = set()

        for a, b in pivot.pairs:

            if pair_has_rectangle_support(
                a,
                b,
                top,
                left,
                opposite,
            ):

                kept.add(
                    (a, b)
                )

        if len(kept) != len(
            pivot.pairs
        ):

            removed += (
                len(pivot.pairs)
                - len(kept)
            )

            relations[pivot_cell] = (
                rebuild_relation(
                    pivot,
                    kept,
                )
            )

            changed = True

    return removed, changed


# ============================================================
# FULL RECTANGLE CLOSURE
# ============================================================

@dataclass
class ClosureStats:

    passes: int = 0
    rectangles_processed: int = 0
    removed: int = 0
    max_removed_one_pass: int = 0


def rectangle_closure(
    relations: Dict[Tuple[int, int], Relation],
    rectangles: List[Rectangle],
) -> ClosureStats:

    stats = ClosureStats()

    while True:

        stats.passes += 1

        removed_this_pass = 0
        changed = False

        for rect in rectangles:

            stats.rectangles_processed += 1

            removed, did_change = (
                prune_rectangle(
                    rect,
                    relations,
                )
            )

            removed_this_pass += removed

            if did_change:
                changed = True

            i, k, j, l = rect

            if (
                not relations[(i, j)].pairs
                or not relations[(i, l)].pairs
                or not relations[(k, j)].pairs
                or not relations[(k, l)].pairs
            ):

                return stats

        stats.removed += (
            removed_this_pass
        )

        stats.max_removed_one_pass = max(
            stats.max_removed_one_pass,
            removed_this_pass,
        )

        if not changed:
            break

    return stats


# ============================================================
# AC3
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
# ANCHOR PROPAGATION
# ============================================================

def propagate_anchor(
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
# FACTOR VERIFICATION
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

    estimated: int


def domain_product(
    domains: List[Set[int]],
    indices: Tuple[int, ...],
) -> int:

    result = 1

    for i in indices:
        result *= len(
            domains[i]
        )

    return result


def modulus_product(
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

            t_modulus = modulus_product(
                Y_moduli,
                y_indices,
            )

            if t_modulus <= t_max:
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
                t_modulus=t_modulus,
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
            "No CRT plan available"
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

        if b % s == 0:

            if n_mod == 0:
                result.update(
                    range(s)
                )

            continue

        b_inv = pow(
            b,
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

        result.add(
            t
        )

    return result


# ============================================================
# INCREMENTAL T CRT
# ============================================================

@dataclass
class TStats:

    states_created: int = 0
    max_states: int = 0
    candidates: int = 0


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
            return []

    result = set()

    for residue, modulus in states:

        if modulus > t_max:

            if residue <= t_max:
                result.add(
                    residue
                )

        else:

            kmax = (
                t_max - residue
            ) // modulus

            for k in range(
                kmax + 1
            ):

                result.add(
                    residue
                    + k * modulus
                )

    stats.candidates += len(
        result
    )

    return sorted(result)


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
    max_t_states: int = 0


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

            candidates = incremental_t_search(
                constraints,
                t_max,
                tstats,
            )

        except RuntimeError:

            continue

        stats.t_states += (
            tstats.states_created
        )

        stats.candidates += (
            len(candidates)
        )

        stats.max_t_states = max(
            stats.max_t_states,
            tstats.max_states,
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

    original_total = sum(
        len(rel.pairs)
        for rel in relations.values()
    )

    rectangles = all_rectangles()

    # --------------------------------------------------------
    # FULL 100-rectangle closure.
    # --------------------------------------------------------

    closure_start = time.perf_counter()

    closure_stats = rectangle_closure(
        relations,
        rectangles,
    )

    closure_time = (
        time.perf_counter()
        - closure_start
    )

    surviving_total = sum(
        len(rel.pairs)
        for rel in relations.values()
    )

    removed_total = (
        original_total
        - surviving_total
    )

    # --------------------------------------------------------
    # Strongest surviving cell.
    # --------------------------------------------------------

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

    (
        rare_pairs,
        rare_products,
        anchor_i,
        anchor_j,
        anchor_rel,
    ) = candidates[0]

    true_A = [
        p % r
        for r in R1
    ]

    true_B = [
        q % s
        for s in R2
    ]

    true_anchor = (
        true_A[anchor_i],
        true_B[anchor_j],
    )

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
    # Search.
    # --------------------------------------------------------

    start = time.perf_counter()

    anchor_attempts = 0
    ac3_survivors = 0

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

        if anchor_attempts > MAX_ANCHOR_ATTEMPTS:

            break

        ok, A, B = (
            propagate_anchor(
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
        # Pick cheaper CRT direction.
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
    # OUTPUT
    # --------------------------------------------------------

    print()
    print("=" * 86)
    print(
        f"EXPERIMENT 172 SAMPLE {sample_id}"
    )
    print("=" * 86)

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
        "--- FULL 2x2 CLOSURE ---"
    )

    print(
        f"rectangles                  = "
        f"{len(rectangles)}"
    )

    print(
        f"closure passes              = "
        f"{closure_stats.passes}"
    )

    print(
        f"rectangles processed        = "
        f"{closure_stats.rectangles_processed}"
    )

    print(
        f"original relation pairs     = "
        f"{original_total}"
    )

    print(
        f"surviving relation pairs    = "
        f"{surviving_total}"
    )

    print(
        f"pairs removed               = "
        f"{removed_total}"
    )

    if original_total:

        print(
            f"survival ratio              = "
            f"{surviving_total / original_total:.6f}"
        )

    print(
        f"max removed/pass            = "
        f"{closure_stats.max_removed_one_pass}"
    )

    print(
        f"closure time                = "
        f"{closure_time:.6f}s"
    )

    print()
    print(
        "--- RAREST SURVIVING CELL ---"
    )

    print(
        f"cell                        = "
        f"({anchor_i},{anchor_j})"
    )

    print(
        f"r                           = "
        f"{R1[anchor_i]}"
    )

    print(
        f"s                           = "
        f"{R2[anchor_j]}"
    )

    print(
        f"C                           = "
        f"{anchor_rel.target_c}"
    )

    print(
        f"pairs                       = "
        f"{rare_pairs}"
    )

    print(
        f"products                    = "
        f"{rare_products}"
    )

    print(
        f"true anchor                 = "
        f"{true_anchor}"
    )

    print(
        f"true anchor present         = "
        f"{true_anchor_present}"
    )

    print(
        f"true anchor rank            = "
        f"{true_anchor_rank}"
    )

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
        f"{ac3_survivors}"
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
        "--- A -> p ---"
    )

    print(
        f"plans                        = "
        f"{A_stats.plans}"
    )

    print(
        f"X states                     = "
        f"{A_stats.x_states}"
    )

    print(
        f"t constraints                = "
        f"{A_stats.t_constraints}"
    )

    print(
        f"t states                     = "
        f"{A_stats.t_states}"
    )

    print(
        f"candidates                   = "
        f"{A_stats.candidates}"
    )

    print(
        f"divisibility tests           = "
        f"{A_stats.divisions}"
    )

    print(
        f"C-grid tests                 = "
        f"{A_stats.c_tests}"
    )

    print(
        f"max t states                 = "
        f"{A_stats.max_t_states}"
    )

    print()
    print(
        "--- B -> q ---"
    )

    print(
        f"plans                        = "
        f"{B_stats.plans}"
    )

    print(
        f"X states                     = "
        f"{B_stats.x_states}"
    )

    print(
        f"t constraints                = "
        f"{B_stats.t_constraints}"
    )

    print(
        f"t states                     = "
        f"{B_stats.t_states}"
    )

    print(
        f"candidates                   = "
        f"{B_stats.candidates}"
    )

    print(
        f"divisibility tests           = "
        f"{B_stats.divisions}"
    )

    print(
        f"C-grid tests                 = "
        f"{B_stats.c_tests}"
    )

    print(
        f"max t states                 = "
        f"{B_stats.max_t_states}"
    )

    print()
    print(
        f"exact                       = "
        f"{exact}"
    )

    print(
        f"exact original factors      = "
        f"{recovered}"
    )

    print()
    print(
        f"search time                  = "
        f"{search_time:.6f}s"
    )

    print(
        f"total time                   = "
        f"{total_time:.6f}s"
    )

    print("=" * 86)

    return {
        "bits": bits,
        "original_pairs": original_total,
        "surviving_pairs": surviving_total,
        "removed_pairs": removed_total,
        "closure_passes": closure_stats.passes,
        "closure_time": closure_time,
        "rarest_pairs": rare_pairs,
        "true_anchor_present": true_anchor_present,
        "true_anchor_rank": true_anchor_rank,
        "anchor_attempts": anchor_attempts,
        "ac3_survivors": ac3_survivors,
        "true_anchor_survived": true_anchor_survived,
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
        "START EXPERIMENT 172"
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
        "EXPERIMENT 172 SUMMARY"
    )
    print("=" * 120)

    print()

    print(
        "bits | orig | kept | removed | "
        "passes | rarest | attempts | AC3 | "
        "A-X | A-t | B-X | B-t | "
        "exact | recovered | time"
    )

    print("-" * 120)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['original_pairs']:>4} | "
            f"{r['surviving_pairs']:>4} | "
            f"{r['removed_pairs']:>7} | "
            f"{r['closure_passes']:>6} | "
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
        "FINISHED EXPERIMENT 172"
    )


if __name__ == "__main__":
    main()
