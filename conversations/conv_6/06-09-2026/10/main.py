#!/usr/bin/env python3

"""
START EXPERIMENT 168

INCREMENTAL CRT-PARAMETER PRUNING

Experiment 167 showed that:
    - the rarest C-cell is a useful anchor;
    - AC-3 preserves the true residue assignment;
    - the bottleneck is Cartesian CRT enumeration.

Experiment 168 changes the CRT stage.

Suppose we choose some factor-side congruences:

    x == P (mod M)

so every candidate has the form

    x = P + M*t.

For a constraint-side modulus s, with

    y == b (mod s)

and

    x*y == n (mod s),

we obtain

    (P + M*t)*b == n (mod s)

and, for b != 0,

    t == (n*b^(-1) - P) * M^(-1) (mod s).

Therefore each Y-domain induces a small set of allowed
residues for t.

We choose enough Y moduli that their CRT modulus exceeds

    floor(sqrt(n) / M).

Then every t-range contains at most one representative
per CRT solution.

This avoids enumerating the Cartesian product of all
complete A residue vectors.

The search is performed from both sides:

    A -> recover p
    B -> recover q

so no assumption about which prime is smaller is required.

FINISHED EXPERIMENT 168
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

SEED = 1682026

STOP_ON_EXACT = True

# Hard safety limits.
MAX_A_BRANCHES = 2_000_000
MAX_T_BRANCHES = 2_000_000


# ============================================================
# BASIC CRT
# ============================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:

    t = ((a2 - a1) * pow(m1, -1, m2)) % m2

    x = a1 + m1 * t
    m = m1 * m2

    return x % m, m


def crt_many(
    residues: List[int],
    moduli: List[int],
) -> Tuple[int, int]:

    x = 0
    M = 1

    for a, m in zip(residues, moduli):
        x, M = crt_pair(x, M, a, m)

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
    r1: List[int],
    r2: List[int],
) -> List[List[int]]:

    A = [p % r for r in r1]
    B = [q % s for s in r2]

    return [
        [
            exact_C(
                n,
                r,
                s,
                A[i],
                B[j],
            )
            for j, s in enumerate(r2)
        ]
        for i, r in enumerate(r1)
    ]


# ============================================================
# RELATIONS
# ============================================================

@dataclass
class CellRelation:

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
) -> CellRelation:

    pairs: Set[Tuple[int, int]] = set()
    by_a: Dict[int, Set[int]] = {}
    by_b: Dict[int, Set[int]] = {}
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

            pairs.add((a, b))
            products.add(a * b)

            by_a.setdefault(a, set()).add(b)
            by_b.setdefault(b, set()).add(a)

    return CellRelation(
        r=r,
        s=s,
        target_c=target_c,
        pairs=pairs,
        by_a=by_a,
        by_b=by_b,
        products=products,
    )


def build_all_relations(
    n: int,
    c_grid: List[List[int]],
    r1: List[int],
    r2: List[int],
) -> Dict[Tuple[int, int], CellRelation]:

    relations = {}

    for i, r in enumerate(r1):
        for j, s in enumerate(r2):

            relations[(i, j)] = build_relation(
                n,
                r,
                s,
                c_grid[i][j],
            )

    return relations


# ============================================================
# AC-3
# ============================================================

def revise_A(
    i: int,
    j: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relation: CellRelation,
) -> bool:

    old = A[i]
    new = set()

    for a in old:

        supported_b = relation.by_a.get(
            a,
            set(),
        )

        if supported_b & B[j]:
            new.add(a)

    changed = new != old

    if changed:
        A[i] = new

    return changed


def revise_B(
    i: int,
    j: int,
    A: List[Set[int]],
    B: List[Set[int]],
    relation: CellRelation,
) -> bool:

    old = B[j]
    new = set()

    for b in old:

        supported_a = relation.by_b.get(
            b,
            set(),
        )

        if supported_a & A[i]:
            new.add(b)

    changed = new != old

    if changed:
        B[j] = new

    return changed


def ac3(
    A: List[Set[int]],
    B: List[Set[int]],
    relations: Dict[Tuple[int, int], CellRelation],
    nr: int,
    ns: int,
) -> bool:

    queue = []

    for i in range(nr):
        for j in range(ns):

            queue.append(
                ("A", i, j)
            )

            queue.append(
                ("B", i, j)
            )

    while queue:

        kind, i, j = queue.pop()

        relation = relations[(i, j)]

        if kind == "A":

            changed = revise_A(
                i,
                j,
                A,
                B,
                relation,
            )

            if changed:

                if not A[i]:
                    return False

                for jj in range(ns):

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
                relation,
            )

            if changed:

                if not B[j]:
                    return False

                for ii in range(nr):

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
    anchor_a: int,
    anchor_b: int,
    relations: Dict[Tuple[int, int], CellRelation],
    nr: int,
    ns: int,
) -> Tuple[
    bool,
    List[Set[int]],
    List[Set[int]],
]:

    A = [
        set(range(relations[(i, 0)].r))
        for i in range(nr)
    ]

    B = [
        set(range(relations[(0, j)].s))
        for j in range(ns)
    ]

    A[anchor_i] = {anchor_a}
    B[anchor_j] = {anchor_b}

    ok = ac3(
        A,
        B,
        relations,
        nr,
        ns,
    )

    return ok, A, B


# ============================================================
# FULL C CHECK
# ============================================================

def c_grid_compatible(
    n: int,
    p: int,
    q: int,
    c_grid: List[List[int]],
    r1: List[int],
    r2: List[int],
) -> bool:

    for i, r in enumerate(r1):

        a = p % r

        for j, s in enumerate(r2):

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
# SUBSET SELECTION
# ============================================================

def product_lengths(
    domains: List[Set[int]],
    indices: Tuple[int, ...],
) -> int:

    result = 1

    for i in indices:
        result *= len(domains[i])

    return result


def modulus_product(
    moduli: List[int],
    indices: Tuple[int, ...],
) -> int:

    result = 1

    for i in indices:
        result *= moduli[i]

    return result


def choose_search_plan(
    X_domains: List[Set[int]],
    X_moduli: List[int],
    Y_domains: List[Set[int]],
    Y_moduli: List[int],
    sqrt_n: int,
) -> Tuple[
    Tuple[int, ...],
    Tuple[int, ...],
    int,
    int,
    int,
]:
    """
    Choose:

        X subset -> x = P + M*t

        Y subset -> t mod product(Y)

    such that:

        product(Y moduli) > sqrt(n) / M.

    The cost estimate is:

        X Cartesian width * Y Cartesian width.

    """

    nx = len(X_moduli)
    ny = len(Y_moduli)

    best = None

    # Try all nonempty X subsets.
    for xmask in range(1, 1 << nx):

        x_indices = tuple(
            i
            for i in range(nx)
            if xmask & (1 << i)
        )

        x_branches = product_lengths(
            X_domains,
            x_indices,
        )

        if x_branches > MAX_A_BRANCHES:
            continue

        M = modulus_product(
            X_moduli,
            x_indices,
        )

        t_bound = sqrt_n // M

        # We need Y modulus > maximum possible t.
        required = t_bound

        for ymask in range(1, 1 << ny):

            y_indices = tuple(
                j
                for j in range(ny)
                if ymask & (1 << j)
            )

            y_M = modulus_product(
                Y_moduli,
                y_indices,
            )

            if y_M <= required:
                continue

            y_branches = product_lengths(
                Y_domains,
                y_indices,
            )

            if y_branches > MAX_T_BRANCHES:
                continue

            estimated = (
                x_branches
                * y_branches
            )

            candidate = (
                estimated,
                x_branches,
                y_branches,
                x_indices,
                y_indices,
                M,
                y_M,
            )

            if best is None or candidate < best:
                best = candidate

    if best is None:
        raise RuntimeError(
            "No feasible CRT-parameter search plan"
        )

    (
        estimated,
        x_branches,
        y_branches,
        x_indices,
        y_indices,
        M,
        y_M,
    ) = best

    return (
        x_indices,
        y_indices,
        M,
        y_M,
        estimated,
    )


# ============================================================
# BUILD t RESIDUE SET
# ============================================================

def allowed_t_residues(
    n: int,
    P: int,
    M: int,
    s: int,
    B_domain: Set[int],
) -> Set[int]:
    """
    x = P + M*t

    y = b mod s

    x*y == n mod s

    For b != 0:

        P*b + M*t*b == n

        t == (n*b^(-1)-P)*M^(-1) mod s
    """

    result: Set[int] = set()

    P_mod = P % s
    M_mod = M % s

    # Since all R1/R2 radices are distinct primes,
    # M is coprime to s.
    M_inv = pow(M_mod, -1, s)

    n_mod = n % s

    for b in B_domain:

        b_mod = b % s

        if b_mod == 0:

            # Then x*b == 0 mod s.
            # This is possible only when n == 0 mod s.
            if n_mod == 0:
                result.update(range(s))

            continue

        b_inv = pow(
            b_mod,
            -1,
            s,
        )

        t = (
            (n_mod * b_inv - P_mod)
            * M_inv
        ) % s

        result.add(t)

    return result


# ============================================================
# SEARCH ONE SIDE
# ============================================================

@dataclass
class SideSearchStats:

    plans: int = 0
    x_states: int = 0
    t_constraint_calls: int = 0
    t_combinations: int = 0
    candidate_values: int = 0
    divisibility_tests: int = 0
    c_tests: int = 0


def search_factor_side(
    *,
    n: int,
    sqrt_n: int,
    X_domains: List[Set[int]],
    X_moduli: List[int],
    Y_domains: List[Set[int]],
    Y_moduli: List[int],
    c_grid: List[List[int]],
    r1: List[int],
    r2: List[int],
    x_is_p: bool,
    stats: SideSearchStats,
) -> Tuple[int | None, int | None]:

    """
    Recover either p or q.

    X side:
        x = P + M*t

    Y side:
        y residue domains constrain t.

    Returns:
        (x, y)
    """

    (
        x_indices,
        y_indices,
        M,
        T_modulus,
        estimated,
    ) = choose_search_plan(
        X_domains,
        X_moduli,
        Y_domains,
        Y_moduli,
        sqrt_n,
    )

    stats.plans += 1

    # --------------------------------------------------------
    # Enumerate X residue combinations.
    # --------------------------------------------------------

    x_lists = [
        sorted(X_domains[i])
        for i in x_indices
    ]

    for x_residues in itertools.product(*x_lists):

        stats.x_states += 1

        P, actual_M = crt_many(
            list(x_residues),
            [X_moduli[i] for i in x_indices],
        )

        if actual_M != M:
            raise AssertionError(
                "CRT modulus mismatch"
            )

        if P > sqrt_n:
            continue

        t_max = (
            sqrt_n - P
        ) // M

        if t_max < 0:
            continue

        # Because T_modulus > sqrt_n/M,
        # at most one t representative per CRT class
        # lies in [0,t_max].
        if T_modulus <= t_max:
            # The plan theoretically prevents this.
            continue

        # ----------------------------------------------------
        # Construct allowed t residues for every selected Y
        # modulus.
        # ----------------------------------------------------

        t_sets = []

        impossible = False

        for j in y_indices:

            stats.t_constraint_calls += 1

            residues = allowed_t_residues(
                n,
                P,
                M,
                Y_moduli[j],
                Y_domains[j],
            )

            if not residues:
                impossible = True
                break

            t_sets.append(
                sorted(residues)
            )

        if impossible:
            continue

        branch_count = 1

        for values in t_sets:
            branch_count *= len(values)

        if branch_count > MAX_T_BRANCHES:
            continue

        stats.t_combinations += branch_count

        # ----------------------------------------------------
        # CRT the t residues.
        # ----------------------------------------------------

        for t_residues in itertools.product(*t_sets):

            t, t_M = crt_many(
                list(t_residues),
                [
                    Y_moduli[j]
                    for j in y_indices
                ],
            )

            if t_M != T_modulus:
                raise AssertionError(
                    "t CRT modulus mismatch"
                )

            if t > t_max:
                continue

            x = P + M * t

            if x <= 1:
                continue

            if x > sqrt_n:
                continue

            stats.candidate_values += 1
            stats.divisibility_tests += 1

            if n % x != 0:
                continue

            y = n // x

            stats.c_tests += 1

            if x_is_p:

                if not c_grid_compatible(
                    n,
                    x,
                    y,
                    c_grid,
                    r1,
                    r2,
                ):
                    continue

            else:

                if not c_grid_compatible(
                    n,
                    y,
                    x,
                    c_grid,
                    r1,
                    r2,
                ):
                    continue

            return x, y

    return None, None


# ============================================================
# SAMPLE GENERATION
# ============================================================

def make_semiprime(
    bits: int,
    rng: random.Random,
) -> Tuple[int, int, int]:

    lo = 1 << (bits // 2 - 1)
    hi = 1 << (bits // 2 + 1)

    p = int(randprime(lo, hi))

    q = int(randprime(lo, hi))

    while q == p:
        q = int(randprime(lo, hi))

    return p, q, p * q


# ============================================================
# ONE SAMPLE
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
        R1,
        R2,
    )

    relations = build_all_relations(
        n,
        c_grid,
        R1,
        R2,
    )

    # --------------------------------------------------------
    # Find rarest cell.
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
        key=lambda z: (
            z[0],
            z[1],
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

    true_anchor_rank = None

    for rank, pair in enumerate(
        sorted(anchor_rel.pairs),
        start=1,
    ):

        if pair == true_anchor:
            true_anchor_rank = rank
            break

    print()
    print("=" * 78)
    print(f"EXPERIMENT 168 SAMPLE {sample_id}")
    print("=" * 78)

    print(f"bits              = {bits}")
    print(f"p                 = {p}")
    print(f"q                 = {q}")
    print(f"n                 = {n}")
    print(f"sqrt(n)           = {sqrt_n}")

    print()
    print(f"true A            = {true_A}")
    print(f"true B            = {true_B}")

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

    print(f"target C          = {anchor_rel.target_c}")
    print(f"allowed pairs     = {rare_pairs}")
    print(f"allowed products  = {rare_products}")
    print(f"true anchor pair  = {true_anchor}")
    print(f"true anchor rank  = {true_anchor_rank}")

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

    start = time.perf_counter()

    anchor_attempts = 0
    ac3_survivors = 0

    true_anchor_survived = False

    true_A_after = None
    true_B_after = None

    A_stats = SideSearchStats()
    B_stats = SideSearchStats()

    exact = None

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
            len(R1),
            len(R2),
        )

        if not ok:
            continue

        ac3_survivors += 1

        if (
            anchor_a == true_anchor[0]
            and anchor_b == true_anchor[1]
        ):

            true_anchor_survived = True

            true_A_after = [
                set(x)
                for x in A
            ]

            true_B_after = [
                set(x)
                for x in B
            ]

        # ----------------------------------------------------
        # Search p side.
        # ----------------------------------------------------

        pp, qq = search_factor_side(
            n=n,
            sqrt_n=sqrt_n,
            X_domains=A,
            X_moduli=R1,
            Y_domains=B,
            Y_moduli=R2,
            c_grid=c_grid,
            r1=R1,
            r2=R2,
            x_is_p=True,
            stats=A_stats,
        )

        if pp is not None:

            exact = (
                pp,
                qq,
                "A->p",
                anchor_a,
                anchor_b,
            )

            break

        # ----------------------------------------------------
        # Search q side.
        # ----------------------------------------------------

        qq, pp = search_factor_side(
            n=n,
            sqrt_n=sqrt_n,
            X_domains=B,
            X_moduli=R2,
            Y_domains=A,
            Y_moduli=R1,
            c_grid=c_grid,
            r1=R1,
            r2=R2,
            x_is_p=False,
            stats=B_stats,
        )

        if qq is not None:

            exact = (
                pp,
                qq,
                "B->q",
                anchor_a,
                anchor_b,
            )

            break

    elapsed = time.perf_counter() - start

    # --------------------------------------------------------
    # Search result.
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

    if true_A_after is not None:

        print()
        print("TRUE ANCHOR DOMAINS")

        print(
            "A widths                     =",
            [len(x) for x in true_A_after],
        )

        print(
            "B widths                     =",
            [len(x) for x in true_B_after],
        )

        print(
            "true A still present         =",
            all(
                true_A[i] in true_A_after[i]
                for i in range(len(R1))
            ),
        )

        print(
            "true B still present         =",
            all(
                true_B[j] in true_B_after[j]
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
        f"{A_stats.t_constraint_calls}"
    )

    print(
        f"t CRT combinations           = "
        f"{A_stats.t_combinations}"
    )

    print(
        f"candidate values             = "
        f"{A_stats.candidate_values}"
    )

    print(
        f"divisibility tests           = "
        f"{A_stats.divisibility_tests}"
    )

    print(
        f"C-grid tests                 = "
        f"{A_stats.c_tests}"
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
        f"{B_stats.t_constraint_calls}"
    )

    print(
        f"t CRT combinations           = "
        f"{B_stats.t_combinations}"
    )

    print(
        f"candidate values             = "
        f"{B_stats.candidate_values}"
    )

    print(
        f"divisibility tests           = "
        f"{B_stats.divisibility_tests}"
    )

    print(
        f"C-grid tests                 = "
        f"{B_stats.c_tests}"
    )

    print()
    print(f"exact                         = {exact}")

    recovered = False

    if exact is not None:

        fp, fq, direction, aa, bb = exact

        recovered = (
            fp * fq == n
            and {
                fp,
                fq,
            } == {
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

    print("=" * 78)

    return {
        "bits": bits,
        "p": p,
        "q": q,
        "n": n,
        "rarest_pairs": rare_pairs,
        "rarest_products": rare_products,
        "true_anchor_rank": true_anchor_rank,
        "anchor_attempts": anchor_attempts,
        "ac3_survivors": ac3_survivors,
        "true_anchor_survived": true_anchor_survived,
        "A_x_states": A_stats.x_states,
        "A_t_combinations": A_stats.t_combinations,
        "A_candidates": A_stats.candidate_values,
        "B_x_states": B_stats.x_states,
        "B_t_combinations": B_stats.t_combinations,
        "B_candidates": B_stats.candidate_values,
        "exact": exact,
        "recovered": recovered,
        "runtime": elapsed,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 168")
    print()
    print("R1 =", R1)
    print("R2 =", R2)
    print("BITS =", BITS_LIST)
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
    print("=" * 90)
    print("EXPERIMENT 168 SUMMARY")
    print("=" * 90)

    print()
    print(
        "bits | rare-pairs | attempts | AC3 | "
        "A-states | A-t | B-states | B-t | "
        "exact | recovered | time"
    )

    print("-" * 90)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['rarest_pairs']:>10} | "
            f"{r['anchor_attempts']:>8} | "
            f"{r['ac3_survivors']:>3} | "
            f"{r['A_x_states']:>8} | "
            f"{r['A_t_combinations']:>3} | "
            f"{r['B_x_states']:>8} | "
            f"{r['B_t_combinations']:>3} | "
            f"{str(r['exact'] is not None):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['runtime']:.6f}s"
        )

    print()
    print("FINISHED EXPERIMENT 168")


if __name__ == "__main__":
    main()
