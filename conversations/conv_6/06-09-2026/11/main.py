#!/usr/bin/env python3

"""
START EXPERIMENT 169

ADAPTIVE BIDIRECTIONAL CRT-TREE PRUNING

Experiment 168 showed:

    x = P + M*t

is much better than enumerating the complete Cartesian product of
all residue domains.

Experiment 169 goes one step further.

For every propagated anchor:

1. Estimate both directions:

       A -> p
       B -> q

2. Select the cheaper direction.

3. Choose an X-side CRT subset.

4. Build the t constraints from the opposite-side residue domains.

5. Instead of computing the full Cartesian product of all t residue
   sets, construct the t CRT state incrementally.

6. After each CRT level, once the current modulus exceeds t_max,
   every remaining CRT class has at most one representative inside
   [0, t_max].

7. Deduplicate CRT states modulo the current modulus.

8. Order t constraints by estimated branching factor.

9. Test exact divisibility as soon as a complete candidate appears.

This experiment specifically measures whether the remaining
t-CRT explosion from Experiment 168 can be compressed further.

FINISHED EXPERIMENT 169
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

SEED = 1692026

STOP_ON_EXACT = True

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

    t = ((a2 - a1) * pow(m1, -1, m2)) % m2

    x = a1 + m1 * t
    M = m1 * m2

    return x % M, M


def crt_many(
    residues: List[int],
    moduli: List[int],
) -> Tuple[int, int]:

    x = 0
    M = 1

    for a, m in zip(residues, moduli):
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
# RELATIONS
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
    by_a = {}
    by_b = {}
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

            pairs.add((a, b))
            products.add(a * b)

            by_a.setdefault(
                a,
                set(),
            ).add(b)

            by_b.setdefault(
                b,
                set(),
            ).add(a)

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

            result[(i, j)] = build_relation(
                n,
                r,
                s,
                c_grid[i][j],
            )

    return result


# ============================================================
# AC-3
# ============================================================

def revise_A(
    i: int,
    j: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relation: Relation,
) -> bool:

    old = A[i]

    new = {
        a
        for a in old
        if relation.by_a.get(a, set()) & B[j]
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
    relation: Relation,
) -> bool:

    old = B[j]

    new = {
        b
        for b in old
        if relation.by_b.get(b, set()) & A[i]
    }

    if new == old:
        return False

    B[j] = new
    return True


def ac3(
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

        rel = relations[(i, j)]

        if kind == "A":

            changed = revise_A(
                i,
                j,
                A,
                B,
                rel,
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
                rel,
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
# C CHECK
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
# SUBSET COST
# ============================================================

def subset_product(
    domains: List[Set[int]],
    indices: Tuple[int, ...],
) -> int:

    value = 1

    for i in indices:
        value *= len(domains[i])

    return value


def subset_modulus(
    moduli: List[int],
    indices: Tuple[int, ...],
) -> int:

    value = 1

    for i in indices:
        value *= moduli[i]

    return value


@dataclass
class Plan:

    x_indices: Tuple[int, ...]
    y_indices: Tuple[int, ...]

    M: int
    t_modulus: int

    estimated_x: int
    estimated_t: int
    estimated_total: int


def choose_best_plan(
    X: List[Set[int]],
    X_moduli: List[int],
    Y: List[Set[int]],
    Y_moduli: List[int],
    sqrt_n: int,
) -> Plan:

    best = None

    nx = len(X_moduli)
    ny = len(Y_moduli)

    for xmask in range(
        1,
        1 << nx,
    ):

        x_indices = tuple(
            i
            for i in range(nx)
            if xmask & (1 << i)
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

        t_max_estimate = sqrt_n // M

        for ymask in range(
            1,
            1 << ny,
        ):

            y_indices = tuple(
                j
                for j in range(ny)
                if ymask & (1 << j)
            )

            Y_modulus = subset_modulus(
                Y_moduli,
                y_indices,
            )

            if Y_modulus <= t_max_estimate:
                continue

            t_cost = subset_product(
                Y,
                y_indices,
            )

            if t_cost > MAX_T_STATES:
                continue

            total = x_cost * t_cost

            candidate = Plan(
                x_indices=x_indices,
                y_indices=y_indices,
                M=M,
                t_modulus=Y_modulus,
                estimated_x=x_cost,
                estimated_t=t_cost,
                estimated_total=total,
            )

            if (
                best is None
                or candidate.estimated_total
                < best.estimated_total
            ):

                best = candidate

    if best is None:
        raise RuntimeError(
            "No CRT search plan found"
        )

    return best


# ============================================================
# T RESIDUE GENERATION
# ============================================================

def t_residues_for_modulus(
    n: int,
    P: int,
    M: int,
    s: int,
    B_domain: Set[int],
) -> Set[int]:

    """
    x = P + M*t

    x*b == n mod s

    =>
        t == (n*b^-1 - P) M^-1 mod s
    """

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
class TSearchResult:

    candidates: List[int]

    residue_sets: List[Tuple[int, int]]

    states_created: int
    states_after_each_level: List[int]


def incremental_t_crt(
    t_constraints: List[Tuple[int, Set[int]]],
    t_max: int,
    stats: TSearchResult | None = None,
) -> TSearchResult:

    """
    Build:

        t == residue_i mod modulus_i

    incrementally.

    States are represented as:

        (t0, M)

    with:

        t == t0 mod M.

    Once M > t_max, each residue class has at most one
    representative inside [0, t_max].

    We then materialize that representative immediately.

    Constraints are sorted by the number of allowed residues.
    """

    if not t_constraints:

        return TSearchResult(
            candidates=list(range(t_max + 1)),
            residue_sets=[],
            states_created=0,
            states_after_each_level=[],
        )

    ordered = sorted(
        t_constraints,
        key=lambda item: (
            len(item[1]),
            item[0],
        ),
    )

    states = {
        (0, 1)
    }

    states_created = 0
    states_after_each_level = []

    used = []

    for modulus, residue_set in ordered:

        used.append(
            (modulus, len(residue_set))
        )

        next_states = set()

        for current, current_M in states:

            for residue in residue_set:

                new_t, new_M = crt_pair(
                    current,
                    current_M,
                    residue,
                    modulus,
                )

                states_created += 1

                if new_M > t_max:

                    # There can be at most one integer
                    # representative <= t_max.

                    candidate = new_t

                    if candidate <= t_max:

                        next_states.add(
                            (
                                candidate,
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
                        "T CRT state limit exceeded"
                    )

        states = next_states

        states_after_each_level.append(
            len(states)
        )

        if not states:
            break

    candidates = []

    # Materialize final representatives.
    for residue, modulus in states:

        if modulus > t_max:

            if residue <= t_max:
                candidates.append(residue)

        else:

            # The final modulus is not always > t_max.
            # Enumerate all representatives in range.
            kmax = (
                t_max - residue
            ) // modulus

            for k in range(
                kmax + 1
            ):

                candidates.append(
                    residue + k * modulus
                )

    candidates = sorted(
        set(candidates)
    )

    return TSearchResult(
        candidates=candidates,
        residue_sets=used,
        states_created=states_created,
        states_after_each_level=states_after_each_level,
    )


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


def side_cost_estimate(
    X: List[Set[int]],
    X_moduli: List[int],
    Y: List[Set[int]],
    Y_moduli: List[int],
    sqrt_n: int,
) -> Plan:

    return choose_best_plan(
        X,
        X_moduli,
        Y,
        Y_moduli,
        sqrt_n,
    )


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
) -> Tuple[int | None, int | None, Plan | None]:

    try:

        plan = side_cost_estimate(
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
        sorted(X[i])
        for i in plan.x_indices
    ]

    for residues in itertools.product(
        *x_lists
    ):

        stats.x_states += 1

        P, M = crt_many(
            list(residues),
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

            opposite = Y[j]

            T = t_residues_for_modulus(
                n,
                P,
                M,
                Y_moduli[j],
                opposite,
            )

            stats.t_constraints += 1

            if not T:

                impossible = True
                break

            constraints.append(
                (
                    Y_moduli[j],
                    T,
                )
            )

        if impossible:
            continue

        try:

            result = incremental_t_crt(
                constraints,
                t_max,
            )

        except RuntimeError:

            continue

        stats.t_states_created += (
            result.states_created
        )

        stats.t_final_candidates += (
            len(result.candidates)
        )

        if result.states_after_each_level:

            stats.max_t_level_states = max(
                stats.max_t_level_states,
                max(
                    result.states_after_each_level
                ),
            )

        for t in result.candidates:

            x = P + M * t

            if x <= 1:
                continue

            if x > sqrt_n:
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

                qq = x
                pp = y

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
# SAMPLE
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
# RUN SAMPLE
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
    # Rarest cell.
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

    rare_pairs, rare_products, anchor_i, anchor_j, anchor_rel = ranked[0]

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

    print()
    print("=" * 80)
    print(
        f"EXPERIMENT 169 SAMPLE {sample_id}"
    )
    print("=" * 80)

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
    print("C-grid:")

    for row in c_grid:
        print("   ", row)

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
    print("five rarest cells:")

    for rank, item in enumerate(
        ranked[:5],
        start=1,
    ):

        pc, zc, i, j, rel = item

        print(
            f"  #{rank}: "
            f"cell=({i},{j}) "
            f"r={R1[i]} "
            f"s={R2[j]} "
            f"C={rel.target_c} "
            f"pairs={pc} "
            f"products={zc}"
        )

    # --------------------------------------------------------
    # Anchor search.
    # --------------------------------------------------------

    anchor_attempts = 0
    ac3_survivors = 0

    A_stats = SideStats()
    B_stats = SideStats()

    true_anchor_survived = False

    true_A_domains = None
    true_B_domains = None

    exact = None

    start = time.perf_counter()

    # --------------------------------------------------------
    # We estimate both directions for the true domain sizes
    # only after propagation.  Then choose the cheaper one.
    # --------------------------------------------------------

    for anchor_a, anchor_b in sorted(
        anchor_rel.pairs
    ):

        anchor_attempts += 1

        ok, A, B = propagate_anchor(
            anchor_i,
            anchor_j,
            anchor_a,
            anchor_b,
            relations,
        )

        if not ok:
            continue

        ac3_survivors += 1

        if (
            anchor_a == true_anchor[0]
            and anchor_b == true_anchor[1]
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
        # Estimate both sides.
        # ----------------------------------------------------

        plan_A = None
        plan_B = None

        try:

            plan_A = side_cost_estimate(
                A,
                R1,
                B,
                R2,
                sqrt_n,
            )

        except RuntimeError:
            pass

        try:

            plan_B = side_cost_estimate(
                B,
                R2,
                A,
                R1,
                sqrt_n,
            )

        except RuntimeError:
            pass

        if plan_A is None and plan_B is None:
            continue

        # Search cheapest side first.
        if (
            plan_B is None
            or (
                plan_A is not None
                and plan_A.estimated_total
                <= plan_B.estimated_total
            )
        ):

            first = "A"
            first_plan = plan_A

        else:

            first = "B"
            first_plan = plan_B

        print(
            f"\nanchor {anchor_attempts}: "
            f"({anchor_a},{anchor_b}) "
            f"survived"
        )

        print(
            f"  A widths = "
            f"{[len(x) for x in A]}"
        )

        print(
            f"  B widths = "
            f"{[len(x) for x in B]}"
        )

        print(
            f"  estimated A cost = "
            f"{plan_A.estimated_total if plan_A else None}"
        )

        print(
            f"  estimated B cost = "
            f"{plan_B.estimated_total if plan_B else None}"
        )

        print(
            f"  selected direction = "
            f"{first}"
        )

        # ----------------------------------------------------
        # First direction.
        # ----------------------------------------------------

        if first == "A":

            pp, qq, used_plan = search_side(
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

            qq, pp, used_plan = search_side(
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
                anchor_a,
                anchor_b,
            )

            break

        # ----------------------------------------------------
        # Try the opposite direction only if necessary.
        # ----------------------------------------------------

        if first == "A":

            qq, pp, used_plan = search_side(
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
                    "B->factor",
                    anchor_a,
                    anchor_b,
                )

                break

        else:

            pp, qq, used_plan = search_side(
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
                    "A->factor",
                    anchor_a,
                    anchor_b,
                )

                break

    elapsed = time.perf_counter() - start

    # --------------------------------------------------------
    # Results.
    # --------------------------------------------------------

    print()
    print("--- SEARCH RESULT ---")

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
        print("TRUE ANCHOR DOMAINS")

        print(
            "A widths                     =",
            [len(x) for x in true_A_domains],
        )

        print(
            "B widths                     =",
            [len(x) for x in true_B_domains],
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
    print("--- A -> p STATISTICS ---")

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
    print("--- B -> q STATISTICS ---")

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

        fp, fq, direction, aa, bb = exact

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

    print("=" * 80)

    return {
        "bits": bits,
        "rarest_pairs": rare_pairs,
        "rarest_products": rare_products,
        "true_anchor_rank": true_rank,
        "anchor_attempts": anchor_attempts,
        "ac3_survivors": ac3_survivors,
        "true_anchor_survived": true_anchor_survived,
        "A_x_states": A_stats.x_states,
        "A_t_states": A_stats.t_states_created,
        "A_candidates": A_stats.t_final_candidates,
        "A_divisions": A_stats.divisibility_tests,
        "B_x_states": B_stats.x_states,
        "B_t_states": B_stats.t_states_created,
        "B_candidates": B_stats.t_final_candidates,
        "B_divisions": B_stats.divisibility_tests,
        "exact": exact,
        "recovered": recovered,
        "runtime": elapsed,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 169")
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

        results.append(result)

    print()
    print("=" * 110)
    print("EXPERIMENT 169 SUMMARY")
    print("=" * 110)

    print()
    print(
        "bits | rare | attempts | AC3 | "
        "A-X | A-t | A-cand | "
        "B-X | B-t | B-cand | "
        "exact | recovered | time"
    )

    print("-" * 110)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['rarest_pairs']:>4} | "
            f"{r['anchor_attempts']:>8} | "
            f"{r['ac3_survivors']:>3} | "
            f"{r['A_x_states']:>4} | "
            f"{r['A_t_states']:>5} | "
            f"{r['A_candidates']:>7} | "
            f"{r['B_x_states']:>4} | "
            f"{r['B_t_states']:>5} | "
            f"{r['B_candidates']:>7} | "
            f"{str(r['exact'] is not None):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['runtime']:.6f}s"
        )

    print()
    print("FINISHED EXPERIMENT 169")


if __name__ == "__main__":
    main()
