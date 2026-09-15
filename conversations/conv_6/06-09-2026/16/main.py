#!/usr/bin/env python3

"""
START EXPERIMENT 173

INTEGER PRODUCT-MATRIX RANK-1 CLOSURE

Core identity
-------------

From Experiment 162:

    C(i,j) depends only on

        z_ij = a_i * b_j

where

        a_i = p mod r_i
        b_j = q mod s_j.

Therefore the complete product matrix is rank 1 over Z:

        z_ij = a_i b_j

and every 2x2 minor vanishes:

        z_ij * z_kl = z_il * z_kj.

Experiment 173 works directly with the allowed product sets:

        Z_ij = { a*b :
                 C(n,r_i,s_j,a,b) = observed C_ij }.

For every 2x2 rectangle:

        Z00 * Z11 = Z01 * Z10

is required.

A product value z survives in a cell only if it can participate in
at least one rank-1 completion of the rectangle.

This is deliberately different from Experiment 172:

    172:
        works with complete (a,b) pairs.

    173:
        works with integer products z=a*b.

The product representation is much smaller in many cells.

After product-level closure:

    1. Rebuild the compatible (a,b) pairs by requiring a*b in Z_ij.
    2. Run ordinary AC-3.
    3. Select the rarest surviving cell.
    4. Anchor on surviving (a,b) pairs.
    5. Use the incremental CRT-parameter search.
    6. Verify exact factorization and complete C-grid.

The experiment specifically measures:

    - original product-set size;
    - product-set reduction;
    - whether the true product matrix survives;
    - resulting residue-domain widths;
    - CRT workload;
    - exact recovery.

FINISHED EXPERIMENT 173
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

SEED = 1732026

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

    return (
        x % (m1 * m2),
        m1 * m2,
    )


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
# PRODUCT RELATION
# ============================================================

@dataclass
class ProductRelation:

    r: int
    s: int
    target_c: int

    products: Set[int]


def build_product_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> ProductRelation:

    products = set()

    for a in range(r):

        for b in range(s):

            if exact_C(
                n,
                r,
                s,
                a,
                b,
            ) == target_c:

                products.add(
                    a * b
                )

    return ProductRelation(
        r=r,
        s=s,
        target_c=target_c,
        products=products,
    )


def build_product_relations(
    n: int,
    c_grid: List[List[int]],
) -> Dict[Tuple[int, int], ProductRelation]:

    result = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            result[(i, j)] = (
                build_product_relation(
                    n,
                    r,
                    s,
                    c_grid[i][j],
                )
            )

    return result


# ============================================================
# PAIR RELATIONS
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


def products_to_relation(
    original: ProductRelation,
    products: Set[int],
) -> Relation:

    pairs = set()

    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}

    for a in range(
        original.r
    ):

        for b in range(
            original.s
        ):

            if (
                a * b
                not in products
            ):
                continue

            # The zero product is technically possible in the
            # residue relation, so do not discard a=0 or b=0 here.
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
        r=original.r,
        s=original.s,
        target_c=original.target_c,
        pairs=pairs,
        by_a=by_a,
        by_b=by_b,
        products=set(products),
    )


def build_pair_relations(
    product_relations: Dict[
        Tuple[int, int],
        ProductRelation,
    ],
) -> Dict[
    Tuple[int, int],
    Relation,
]:

    return {
        cell: products_to_relation(
            rel,
            rel.products,
        )
        for cell, rel
        in product_relations.items()
    }


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
                        (i, k, j, l)
                    )

    return result


# ============================================================
# PRODUCT-MINOR SUPPORT
# ============================================================

def build_cross_product_set(
    X: Set[int],
    Y: Set[int],
) -> Set[int]:

    """
    Return all integer products x*y.

    Sets are small enough here that direct construction is practical.
    """

    result = set()

    for x in X:

        for y in Y:

            result.add(
                x * y
            )

    return result


# ============================================================
# RECTANGLE PRODUCT CLOSURE
# ============================================================

def prune_product_rectangle(
    rect: Rectangle,
    relations: Dict[
        Tuple[int, int],
        ProductRelation,
    ],
) -> Tuple[int, bool]:

    """
    For:

        Z00 * Z11 = Z01 * Z10

    remove values which cannot participate in ANY equality.

    We construct the possible cross-product values for both
    diagonal pairs.

    For each z in Z00 we need some w in Z11 such that:

        z*w

    belongs to the cross-product set:

        Z01*Z10.
    """

    i, k, j, l = rect

    c00 = (i, j)
    c01 = (i, l)
    c10 = (k, j)
    c11 = (k, l)

    z00 = relations[c00].products
    z01 = relations[c01].products
    z10 = relations[c10].products
    z11 = relations[c11].products

    if not z00 or not z01 or not z10 or not z11:
        return 0, False

    # --------------------------------------------------------
    # Left diagonal target products.
    # --------------------------------------------------------

    right_products = build_cross_product_set(
        z01,
        z10,
    )

    # --------------------------------------------------------
    # Support for values in Z00 / Z11.
    # --------------------------------------------------------

    keep00 = set()

    # Simple modular divisibility shortcut:
    # if z00 == 0 then every product is 0.
    for z in z00:

        if z == 0:

            if 0 in right_products:
                keep00.add(z)

            continue

        supported = False

        for w in z11:

            if (
                z * w
                in right_products
            ):

                supported = True
                break

        if supported:
            keep00.add(z)

    keep11 = set()

    for w in z11:

        if w == 0:

            if 0 in right_products:
                keep11.add(w)

            continue

        supported = False

        for z in z00:

            if (
                z * w
                in right_products
            ):

                supported = True
                break

        if supported:
            keep11.add(w)

    # --------------------------------------------------------
    # Reverse diagonal.
    #
    # Need:
    #
    #       z01*z10 = z00*z11
    #
    # so use left diagonal products.
    # --------------------------------------------------------

    left_products = build_cross_product_set(
        keep00,
        keep11,
    )

    keep01 = set()

    for z in z01:

        if z == 0:

            if 0 in left_products:
                keep01.add(z)

            continue

        supported = False

        for w in z10:

            if (
                z * w
                in left_products
            ):

                supported = True
                break

        if supported:
            keep01.add(z)

    keep10 = set()

    for w in z10:

        if w == 0:

            if 0 in left_products:
                keep10.add(w)

            continue

        supported = False

        for z in z01:

            if (
                z * w
                in left_products
            ):

                supported = True
                break

        if supported:
            keep10.add(w)

    changed = (
        keep00 != z00
        or keep01 != z01
        or keep10 != z10
        or keep11 != z11
    )

    removed = (
        len(z00)
        - len(keep00)
        + len(z01)
        - len(keep01)
        + len(z10)
        - len(keep10)
        + len(z11)
        - len(keep11)
    )

    if not changed:
        return 0, False

    relations[c00].products = keep00
    relations[c01].products = keep01
    relations[c10].products = keep10
    relations[c11].products = keep11

    return removed, True


# ============================================================
# GLOBAL PRODUCT CLOSURE
# ============================================================

@dataclass
class ProductClosureStats:

    passes: int = 0
    rectangles_processed: int = 0
    removed: int = 0
    max_removed_pass: int = 0


def product_closure(
    relations: Dict[
        Tuple[int, int],
        ProductRelation,
    ],
    rectangles: List[Rectangle],
) -> ProductClosureStats:

    stats = ProductClosureStats()

    while True:

        stats.passes += 1

        removed_pass = 0
        changed = False

        for rect in rectangles:

            stats.rectangles_processed += 1

            removed, did_change = (
                prune_product_rectangle(
                    rect,
                    relations,
                )
            )

            removed_pass += removed

            if did_change:
                changed = True

            i, k, j, l = rect

            if (
                not relations[(i, j)].products
                or not relations[(i, l)].products
                or not relations[(k, j)].products
                or not relations[(k, l)].products
            ):

                return stats

        stats.removed += removed_pass

        stats.max_removed_pass = max(
            stats.max_removed_pass,
            removed_pass,
        )

        if not changed:
            break

    return stats


# ============================================================
# PAIR AC3
# ============================================================

def ac3(
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[
        Tuple[int, int],
        Relation,
    ],
) -> bool:

    changed = True

    while changed:

        changed = False

        for i in range(len(R1)):

            for j in range(len(R2)):

                rel = relations[
                    (i, j)
                ]

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

            estimate = (
                x_cost
                * y_cost
            )

            candidate = Plan(
                x_indices=x_indices,
                y_indices=y_indices,
                M=M,
                t_modulus=T,
                estimated=estimate,
            )

            if (
                best is None
                or estimate
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

    if not constraints:

        return list(
            range(
                t_max + 1
            )
        )

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

    stats.candidates += len(
        candidates
    )

    return sorted(
        candidates
    )


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
        sorted(
            X[i]
        )
        for i in plan.x_indices
    ]

    for values in itertools.product(
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
            list(values),
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

    true_Z = [
        [
            true_A[i] * true_B[j]
            for j in range(len(R2))
        ]
        for i in range(len(R1))
    ]

    c_grid = make_C_grid(
        n,
        p,
        q,
    )

    product_relations = (
        build_product_relations(
            n,
            c_grid,
        )
    )

    original_product_total = sum(
        len(rel.products)
        for rel
        in product_relations.values()
    )

    original_product_sizes = {
        cell: len(rel.products)
        for cell, rel
        in product_relations.items()
    }

    # --------------------------------------------------------
    # Product rank-1 closure.
    # --------------------------------------------------------

    rectangles = all_rectangles()

    closure_start = time.perf_counter()

    closure_stats = product_closure(
        product_relations,
        rectangles,
    )

    closure_time = (
        time.perf_counter()
        - closure_start
    )

    surviving_product_total = sum(
        len(rel.products)
        for rel
        in product_relations.values()
    )

    removed_product_total = (
        original_product_total
        - surviving_product_total
    )

    # --------------------------------------------------------
    # Check true Z matrix survival.
    # --------------------------------------------------------

    true_products_survive = all(
        true_Z[i][j]
        in product_relations[(i, j)].products
        for i in range(len(R1))
        for j in range(len(R2))
    )

    # --------------------------------------------------------
    # Convert products back into pair relations.
    # --------------------------------------------------------

    relations = build_pair_relations(
        product_relations
    )

    # --------------------------------------------------------
    # Find rarest surviving pair relation.
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

    A_stats = SideStats()
    B_stats = SideStats()

    anchor_attempts = 0
    ac3_survivors = 0

    true_anchor_survived = False

    true_A_domains = None
    true_B_domains = None

    exact = None

    search_start = time.perf_counter()

    for a0, b0 in sorted(
        anchor_rel.pairs
    ):

        anchor_attempts += 1

        if (
            anchor_attempts
            > MAX_ANCHOR_ATTEMPTS
        ):
            break

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

        # ----------------------------------------------------
        # First side.
        # ----------------------------------------------------

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
        # Opposite side.
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
        - search_start
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
    print("=" * 88)
    print(
        f"EXPERIMENT 173 SAMPLE {sample_id}"
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
        "--- PRODUCT RANK-1 CLOSURE ---"
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
        f"original product values     = "
        f"{original_product_total}"
    )

    print(
        f"surviving product values    = "
        f"{surviving_product_total}"
    )

    print(
        f"product values removed      = "
        f"{removed_product_total}"
    )

    if original_product_total:

        print(
            f"product survival ratio      = "
            f"{surviving_product_total / original_product_total:.6f}"
        )

    print(
        f"max removed/pass            = "
        f"{closure_stats.max_removed_pass}"
    )

    print(
        f"true product matrix survives= "
        f"{true_products_survive}"
    )

    print(
        f"closure time                = "
        f"{closure_time:.6f}s"
    )

    # --------------------------------------------------------
    # Product statistics.
    # --------------------------------------------------------

    print()
    print(
        "PRODUCT SET SIZES"
    )

    for i in range(len(R1)):

        row = []

        for j in range(len(R2)):

            row.append(
                len(
                    product_relations[
                        (i, j)
                    ].products
                )
            )

        print(
            "   ",
            row,
        )

    # --------------------------------------------------------
    # Rarest relation.
    # --------------------------------------------------------

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

    # --------------------------------------------------------
    # Search.
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

    # --------------------------------------------------------
    # A side.
    # --------------------------------------------------------

    print()
    print(
        "--- A -> p ---"
    )

    print(
        f"plans                       = "
        f"{A_stats.plans}"
    )

    print(
        f"X states                    = "
        f"{A_stats.x_states}"
    )

    print(
        f"t constraints               = "
        f"{A_stats.t_constraints}"
    )

    print(
        f"t states                    = "
        f"{A_stats.t_states}"
    )

    print(
        f"candidates                  = "
        f"{A_stats.candidates}"
    )

    print(
        f"divisibility tests          = "
        f"{A_stats.divisions}"
    )

    print(
        f"C-grid tests                = "
        f"{A_stats.c_tests}"
    )

    print(
        f"max t states                = "
        f"{A_stats.max_t_states}"
    )

    # --------------------------------------------------------
    # B side.
    # --------------------------------------------------------

    print()
    print(
        "--- B -> q ---"
    )

    print(
        f"plans                       = "
        f"{B_stats.plans}"
    )

    print(
        f"X states                    = "
        f"{B_stats.x_states}"
    )

    print(
        f"t constraints               = "
        f"{B_stats.t_constraints}"
    )

    print(
        f"t states                    = "
        f"{B_stats.t_states}"
    )

    print(
        f"candidates                  = "
        f"{B_stats.candidates}"
    )

    print(
        f"divisibility tests          = "
        f"{B_stats.divisions}"
    )

    print(
        f"C-grid tests                = "
        f"{B_stats.c_tests}"
    )

    print(
        f"max t states                = "
        f"{B_stats.max_t_states}"
    )

    # --------------------------------------------------------
    # Final.
    # --------------------------------------------------------

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
        f"search time                 = "
        f"{search_time:.6f}s"
    )

    print(
        f"total time                  = "
        f"{total_time:.6f}s"
    )

    print("=" * 88)

    return {
        "bits": bits,
        "original_products": original_product_total,
        "surviving_products": surviving_product_total,
        "removed_products": removed_product_total,
        "closure_passes": closure_stats.passes,
        "closure_time": closure_time,
        "true_products_survive": true_products_survive,
        "rarest_pairs": rare_pairs,
        "rarest_products": rare_products,
        "true_anchor_present": true_anchor_present,
        "true_anchor_rank": true_anchor_rank,
        "anchor_attempts": anchor_attempts,
        "ac3_survivors": ac3_survivors,
        "true_anchor_survived": true_anchor_survived,
        "A_x_states": A_stats.x_states,
        "A_t_states": A_stats.t_states,
        "A_candidates": A_stats.candidates,
        "B_x_states": B_stats.x_states,
        "B_t_states": B_stats.t_states,
        "B_candidates": B_stats.candidates,
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
        "START EXPERIMENT 173"
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
        "EXPERIMENT 173 SUMMARY"
    )
    print("=" * 120)

    print()
    print(
        "bits | orig-Z | kept-Z | removed | "
        "ratio | rarest | attempts | AC3 | "
        "A-X | A-t | B-X | B-t | "
        "exact | recovered | time"
    )

    print("-" * 120)

    for r in results:

        ratio = (
            r["surviving_products"]
            / r["original_products"]
            if r["original_products"]
            else 0
        )

        print(
            f"{r['bits']:>4} | "
            f"{r['original_products']:>6} | "
            f"{r['surviving_products']:>6} | "
            f"{r['removed_products']:>7} | "
            f"{ratio:>5.3f} | "
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
        "FINISHED EXPERIMENT 173"
    )


if __name__ == "__main__":
    main()
