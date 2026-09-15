#!/usr/bin/env python3

"""
START EXPERIMENT 167

RAREST-CELL MULTIPLICATIVE C-ANCHOR

Goal:
    Replace the previous sqrt(n) candidate scan.

    For each sample:
      1. Build the exact C-grid from the known factorization.
      2. For every cell (r_i, s_j), build the relation
             C(n,r_i,s_j,a,b) = target_C[i][j]
         over all residue pairs (a,b).
      3. Select the cell having the fewest valid (a,b) pairs.
      4. Anchor on each pair in that rarest relation.
      5. Run AC-3 propagation over all A_i / B_j domains.
      6. CRT-lift the surviving domains.
      7. Test exact divisibility.
      8. Verify the recovered factor pair against the complete C-grid.

This is an oracle experiment: the target C-grid is obtained from
the true p,q solely so that we can measure how much information
the C-grid itself contains.

The important difference from the previous experiment is that we
never scan x = 1..sqrt(n).

FINISHED EXPERIMENT 167
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

SEED = 1672026

# Stop after finding the exact factorization.
STOP_ON_EXACT = True

# Safety limit for CRT enumeration.
MAX_CRT_BRANCHES = 5_000_000


# ============================================================
# BASIC MATH
# ============================================================

def crt_pair(a1: int, m1: int, a2: int, m2: int) -> Tuple[int, int]:
    """
    Combine:
        x = a1 mod m1
        x = a2 mod m2

    Assumes gcd(m1,m2)=1.
    """
    t = ((a2 - a1) * pow(m1, -1, m2)) % m2
    x = a1 + m1 * t
    m = m1 * m2
    return x % m, m


def crt_many(residues: List[int], moduli: List[int]) -> Tuple[int, int]:
    x = 0
    m = 1

    for a, mod in zip(residues, moduli):
        x, m = crt_pair(x, m, a, mod)

    return x, m


def exact_C(n: int, r: int, s: int, a: int, b: int) -> int:
    """
    Exact carry C for one radix cell.

    beta = (k*b) mod s
    alpha = (ell*a) mod r

    Using:
        n = (rk+a)(sl+b)

    we have:
        rk*b = n - ab - sl*a
        mod s:
            r*k*b == n-ab mod s

    therefore:
        beta = (n-ab) * r^{-1} mod s

    Similarly:
        alpha = (n-ab) * s^{-1} mod r
    """

    if math.gcd(r, s) != 1:
        raise ValueError(f"r={r}, s={s} are not coprime")

    z = a * b

    beta = ((n - z) * pow(r, -1, s)) % s
    alpha = ((n - z) * pow(s, -1, r)) % r

    return (r * beta + s * alpha + z) // (r * s)


def make_C_grid(
    n: int,
    p: int,
    q: int,
    r1: List[int],
    r2: List[int],
) -> List[List[int]]:
    """
    Build the oracle C-grid from the true factors.
    """

    A = [p % r for r in r1]
    B = [q % s for s in r2]

    grid = []

    for i, r in enumerate(r1):
        row = []

        for j, s in enumerate(r2):
            row.append(exact_C(n, r, s, A[i], B[j]))

        grid.append(row)

    return grid


# ============================================================
# RELATION BUILDING
# ============================================================

@dataclass
class CellRelation:
    r: int
    s: int
    target_c: int

    pairs: Set[Tuple[int, int]]

    # For fast AC-3 support checks.
    by_a: Dict[int, Set[int]]
    by_b: Dict[int, Set[int]]

    # Distinct ab products.
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
            c = exact_C(n, r, s, a, b)

            if c != target_c:
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


# ============================================================
# GRID CONSTRUCTION
# ============================================================

def build_all_relations(
    n: int,
    c_grid: List[List[int]],
    r1: List[int],
    r2: List[int],
) -> Dict[Tuple[int, int], CellRelation]:

    relations = {}

    for i, r in enumerate(r1):
        for j, s in enumerate(r2):
            rel = build_relation(
                n,
                r,
                s,
                c_grid[i][j],
            )

            relations[(i, j)] = rel

    return relations


# ============================================================
# AC-3
# ============================================================

def revise_A(
    i: int,
    j: int,
    A_domains: List[Set[int]],
    B_domains: List[Set[int]],
    relation: CellRelation,
) -> bool:

    old = A_domains[i]

    new_domain = set()

    for a in old:
        # IMPORTANT:
        # .get(a, set()) prevents the old KeyError: 0.
        supported_b = relation.by_a.get(a, set())

        if supported_b & B_domains[j]:
            new_domain.add(a)

    changed = new_domain != old

    if changed:
        A_domains[i] = new_domain

    return changed


def revise_B(
    i: int,
    j: int,
    A_domains: List[Set[int]],
    B_domains: List[Set[int]],
    relation: CellRelation,
) -> bool:

    old = B_domains[j]

    new_domain = set()

    for b in old:
        supported_a = relation.by_b.get(b, set())

        if supported_a & A_domains[i]:
            new_domain.add(b)

    changed = new_domain != old

    if changed:
        B_domains[j] = new_domain

    return changed


def ac3(
    A_domains: List[Set[int]],
    B_domains: List[Set[int]],
    relations: Dict[Tuple[int, int], CellRelation],
    nr: int,
    ns: int,
) -> bool:

    queue = []

    for i in range(nr):
        for j in range(ns):
            queue.append(("A", i, j))
            queue.append(("B", i, j))

    while queue:

        kind, i, j = queue.pop()

        rel = relations[(i, j)]

        if kind == "A":
            changed = revise_A(
                i,
                j,
                A_domains,
                B_domains,
                rel,
            )

            if changed:

                if not A_domains[i]:
                    return False

                for jj in range(ns):
                    if jj != j:
                        queue.append(("B", i, jj))

        else:
            changed = revise_B(
                i,
                j,
                A_domains,
                B_domains,
                rel,
            )

            if changed:

                if not B_domains[j]:
                    return False

                for ii in range(nr):
                    if ii != i:
                        queue.append(("A", ii, j))

    return True


# ============================================================
# FULL C-GRID CHECK
# ============================================================

def c_grid_compatible(
    n: int,
    A: List[int],
    B: List[int],
    c_grid: List[List[int]],
    r1: List[int],
    r2: List[int],
) -> bool:

    for i, r in enumerate(r1):
        for j, s in enumerate(r2):

            if exact_C(
                n,
                r,
                s,
                A[i],
                B[j],
            ) != c_grid[i][j]:

                return False

    return True


# ============================================================
# CRT SUBSET SELECTION
# ============================================================

def choose_best_crt_subset(
    domains: List[Set[int]],
    moduli: List[int],
    target_modulus: int,
) -> Tuple[List[int], int]:

    """
    Find the subset whose modulus product exceeds target_modulus
    while minimizing the estimated number of CRT branches.

    Since we only have 5 radices this exhaustive subset search is cheap.
    """

    best_indices = None
    best_cost = None

    n = len(moduli)

    for mask in range(1, 1 << n):

        indices = [
            i for i in range(n)
            if mask & (1 << i)
        ]

        modulus_product = 1
        branch_cost = 1

        for i in indices:
            modulus_product *= moduli[i]
            branch_cost *= len(domains[i])

            if branch_cost > MAX_CRT_BRANCHES:
                break

        if modulus_product <= target_modulus:
            continue

        if branch_cost > MAX_CRT_BRANCHES:
            continue

        if best_cost is None or branch_cost < best_cost:
            best_cost = branch_cost
            best_indices = indices

    if best_indices is None:
        raise RuntimeError(
            "Could not find a CRT subset exceeding target modulus"
        )

    return best_indices, best_cost


# ============================================================
# CRT ENUMERATION
# ============================================================

def enumerate_crt_candidates(
    domains: List[Set[int]],
    moduli: List[int],
    target_bound: int,
) -> Tuple[List[int], int, int]:

    """
    Enumerate CRT combinations from a chosen subset.

    Returns:
        candidates
        number of branches
        modulus used
    """

    target_modulus = target_bound

    indices, estimated = choose_best_crt_subset(
        domains,
        moduli,
        target_modulus,
    )

    selected_moduli = [moduli[i] for i in indices]

    candidate_domains = [
        sorted(domains[i])
        for i in indices
    ]

    modulus_product = math.prod(selected_moduli)

    if modulus_product <= target_bound:
        raise AssertionError("CRT modulus is not large enough")

    branches = 0
    candidates = []

    for residues in itertools.product(*candidate_domains):

        branches += 1

        if branches > MAX_CRT_BRANCHES:
            raise RuntimeError(
                f"CRT branch limit exceeded: {MAX_CRT_BRANCHES}"
            )

        x, M = crt_many(
            list(residues),
            selected_moduli,
        )

        if 0 < x <= target_bound:
            candidates.append(x)

    return candidates, branches, modulus_product


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
) -> Tuple[bool, List[Set[int]], List[Set[int]]]:

    A_domains = [
        set(range(rel.r))
        for rel in (
            relations[(i, 0)]
            for i in range(nr)
        )
    ]

    B_domains = [
        set(range(rel.s))
        for rel in (
            relations[(0, j)]
            for j in range(ns)
        )
    ]

    # Force anchor.
    A_domains[anchor_i] = {anchor_a}
    B_domains[anchor_j] = {anchor_b}

    ok = ac3(
        A_domains,
        B_domains,
        relations,
        nr,
        ns,
    )

    return ok, A_domains, B_domains


# ============================================================
# SAMPLE GENERATION
# ============================================================

def make_semiprime(bits: int, rng: random.Random) -> Tuple[int, int, int]:

    lo = 1 << (bits // 2 - 1)
    hi = 1 << (bits // 2 + 1)

    p = int(randprime(lo, hi))
    q = int(randprime(lo, hi))

    while q == p:
        q = int(randprime(lo, hi))

    n = p * q

    return p, q, n


# ============================================================
# ONE SAMPLE
# ============================================================

def run_sample(
    bits: int,
    sample_id: int,
    rng: random.Random,
) -> dict:

    p, q, n = make_semiprime(bits, rng)

    sqrt_n = math.isqrt(n)

    true_A = [p % r for r in R1]
    true_B = [q % s for s in R2]

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

    ranked_cells = []

    for (i, j), rel in relations.items():

        ranked_cells.append(
            (
                len(rel.pairs),
                len(rel.products),
                i,
                j,
                rel,
            )
        )

    ranked_cells.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    pair_count, product_count, anchor_i, anchor_j, anchor_rel = ranked_cells[0]

    anchor_true_pair = (
        true_A[anchor_i],
        true_B[anchor_j],
    )

    anchor_true_rank = None

    for rank, pair in enumerate(
        sorted(anchor_rel.pairs),
        start=1,
    ):
        if pair == anchor_true_pair:
            anchor_true_rank = rank
            break

    print()
    print("=" * 72)
    print(f"EXPERIMENT 167 SAMPLE {sample_id}")
    print("=" * 72)

    print(f"bits              = {bits}")
    print(f"p                 = {p}")
    print(f"q                 = {q}")
    print(f"n                 = {n}")
    print(f"sqrt(n)           = {sqrt_n}")

    print()
    print(f"R1                = {R1}")
    print(f"R2                = {R2}")

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
    print(f"allowed pairs     = {pair_count}")
    print(f"allowed products  = {product_count}")
    print(f"true anchor pair  = {anchor_true_pair}")
    print(f"true anchor rank  = {anchor_true_rank}")

    # Show the five rarest cells.
    print()
    print("five rarest cells:")

    for rank, item in enumerate(ranked_cells[:5], start=1):

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
    # Anchor enumeration.
    # --------------------------------------------------------

    start = time.perf_counter()

    anchor_attempts = 0
    propagation_survivors = 0
    exact_hits = []

    total_crt_branches_A = 0
    total_crt_branches_B = 0

    true_anchor_survived = False
    true_anchor_propagation_domains = None

    successful_anchor = None

    for anchor_a, anchor_b in sorted(anchor_rel.pairs):

        anchor_attempts += 1

        ok, A_domains, B_domains = propagate_anchor(
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

        propagation_survivors += 1

        if (anchor_a, anchor_b) == anchor_true_pair:
            true_anchor_survived = True
            true_anchor_propagation_domains = (
                [set(x) for x in A_domains],
                [set(x) for x in B_domains],
            )

        # ----------------------------------------------------
        # Full singleton case first.
        # ----------------------------------------------------

        if all(len(x) == 1 for x in A_domains + B_domains):

            A = [next(iter(x)) for x in A_domains]
            B = [next(iter(x)) for x in B_domains]

            if not c_grid_compatible(
                n,
                A,
                B,
                c_grid,
                R1,
                R2,
            ):
                continue

            # CRT-derived factor candidates.
            p_candidate = None
            q_candidate = None

            try:
                p_candidates, branches_a, modulus_a = enumerate_crt_candidates(
                    A_domains,
                    R1,
                    sqrt_n,
                )

                q_candidates, branches_b, modulus_b = enumerate_crt_candidates(
                    B_domains,
                    R2,
                    sqrt_n,
                )

                total_crt_branches_A += branches_a
                total_crt_branches_B += branches_b

                for pp in p_candidates:
                    if pp != 0 and n % pp == 0:
                        qq = n // pp

                        if qq == q or qq == p:
                            exact_hits.append(
                                (pp, qq)
                            )

                        if (
                            pp * qq == n
                            and pp > 1
                            and qq > 1
                        ):
                            successful_anchor = (
                                anchor_a,
                                anchor_b,
                                pp,
                                qq,
                            )

                            if STOP_ON_EXACT:
                                break

                if successful_anchor and STOP_ON_EXACT:
                    break

            except RuntimeError:
                pass

            continue

        # ----------------------------------------------------
        # CRT lift A side.
        # ----------------------------------------------------

        try:

            p_candidates, branches_a, modulus_a = enumerate_crt_candidates(
                A_domains,
                R1,
                sqrt_n,
            )

            total_crt_branches_A += branches_a

            for candidate_p in p_candidates:

                if candidate_p <= 1:
                    continue

                if n % candidate_p != 0:
                    continue

                candidate_q = n // candidate_p

                if candidate_q <= 1:
                    continue

                # Verify exact C-grid.
                AA = [candidate_p % r for r in R1]
                BB = [candidate_q % s for s in R2]

                if not c_grid_compatible(
                    n,
                    AA,
                    BB,
                    c_grid,
                    R1,
                    R2,
                ):
                    continue

                exact_hits.append(
                    (
                        candidate_p,
                        candidate_q,
                    )
                )

                successful_anchor = (
                    anchor_a,
                    anchor_b,
                    candidate_p,
                    candidate_q,
                )

                if STOP_ON_EXACT:
                    break

            if successful_anchor and STOP_ON_EXACT:
                break

        except RuntimeError:
            pass

        # ----------------------------------------------------
        # CRT lift B side.
        # ----------------------------------------------------

        if successful_anchor and STOP_ON_EXACT:
            break

        try:

            q_candidates, branches_b, modulus_b = enumerate_crt_candidates(
                B_domains,
                R2,
                sqrt_n,
            )

            total_crt_branches_B += branches_b

            for candidate_q in q_candidates:

                if candidate_q <= 1:
                    continue

                if n % candidate_q != 0:
                    continue

                candidate_p = n // candidate_q

                if candidate_p <= 1:
                    continue

                AA = [candidate_p % r for r in R1]
                BB = [candidate_q % s for s in R2]

                if not c_grid_compatible(
                    n,
                    AA,
                    BB,
                    c_grid,
                    R1,
                    R2,
                ):
                    continue

                exact_hits.append(
                    (
                        candidate_p,
                        candidate_q,
                    )
                )

                successful_anchor = (
                    anchor_a,
                    anchor_b,
                    candidate_p,
                    candidate_q,
                )

                if STOP_ON_EXACT:
                    break

        except RuntimeError:
            pass

        if successful_anchor and STOP_ON_EXACT:
            break

    elapsed = time.perf_counter() - start

    # --------------------------------------------------------
    # Output propagation result for true anchor.
    # --------------------------------------------------------

    print()
    print("--- SEARCH RESULT ---")

    print(f"anchor attempts             = {anchor_attempts}")
    print(f"AC-3 survivors              = {propagation_survivors}")
    print(f"true anchor survived        = {true_anchor_survived}")

    if true_anchor_propagation_domains is not None:

        A_dom, B_dom = true_anchor_propagation_domains

        print()
        print("TRUE ANCHOR DOMAINS")
        print(
            "A widths                     =",
            [len(x) for x in A_dom],
        )
        print(
            "B widths                     =",
            [len(x) for x in B_dom],
        )

        print(
            "A domains                    =",
            [sorted(x) for x in A_dom],
        )

        print(
            "B domains                    =",
            [sorted(x) for x in B_dom],
        )

        true_A_survives = all(
            true_A[i] in A_dom[i]
            for i in range(len(R1))
        )

        true_B_survives = all(
            true_B[j] in B_dom[j]
            for j in range(len(R2))
        )

        print(
            f"true A still present         = {true_A_survives}"
        )

        print(
            f"true B still present         = {true_B_survives}"
        )

    print()
    print(
        f"CRT branches A               = "
        f"{total_crt_branches_A}"
    )

    print(
        f"CRT branches B               = "
        f"{total_crt_branches_B}"
    )

    print(
        f"exact hits                   = "
        f"{exact_hits}"
    )

    print(
        f"successful anchor            = "
        f"{successful_anchor}"
    )

    if successful_anchor:

        _, _, fp, fq = successful_anchor

        recovered = (
            fp * fq == n
            and (
                (fp == p and fq == q)
                or
                (fp == q and fq == p)
            )
        )

        print(
            f"exact original factors      = "
            f"{recovered}"
        )

    else:
        recovered = False

    print()
    print(f"runtime                      = {elapsed:.6f}s")

    print("=" * 72)

    return {
        "bits": bits,
        "p": p,
        "q": q,
        "n": n,
        "rarest_cell": (anchor_i, anchor_j),
        "rarest_pairs": pair_count,
        "rarest_products": product_count,
        "true_anchor_rank": anchor_true_rank,
        "anchor_attempts": anchor_attempts,
        "propagation_survivors": propagation_survivors,
        "true_anchor_survived": true_anchor_survived,
        "crt_branches_A": total_crt_branches_A,
        "crt_branches_B": total_crt_branches_B,
        "exact_hits": exact_hits,
        "successful_anchor": successful_anchor,
        "recovered": recovered,
        "runtime": elapsed,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 167")
    print()
    print("R1 =", R1)
    print("R2 =", R2)
    print("BITS =", BITS_LIST)
    print()

    results = []

    for sample_id, bits in enumerate(BITS_LIST, start=1):

        result = run_sample(
            bits,
            sample_id,
            rng,
        )

        results.append(result)

    print()
    print("=" * 72)
    print("EXPERIMENT 167 SUMMARY")
    print("=" * 72)

    print()
    print(
        "bits | rarest | attempts | AC3 | CRT-A | CRT-B | "
        "exact | recovered | time"
    )
    print("-" * 72)

    for r in results:

        print(
            f"{r['bits']:>4} | "
            f"{r['rarest_pairs']:>6} | "
            f"{r['anchor_attempts']:>8} | "
            f"{r['propagation_survivors']:>3} | "
            f"{r['crt_branches_A']:>5} | "
            f"{r['crt_branches_B']:>5} | "
            f"{len(r['exact_hits']):>5} | "
            f"{str(r['recovered']):>9} | "
            f"{r['runtime']:.6f}s"
        )

    print()
    print("FINISHED EXPERIMENT 167")


if __name__ == "__main__":
    main()
