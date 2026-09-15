#!/usr/bin/env python3

# ============================================================
# START EXPERIMENT 163
# Global Multiplicative Rank-1 Recovery from C-Product Sets
# ============================================================

from __future__ import annotations

import math
import random
import time
from collections import deque


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48, 54]
SAMPLES_PER_SIZE = 1

SEED = 163

# Limit when explicitly enumerating complete solutions.
DFS_SOLUTION_CAP = 5000

# How many of the most constrained anchors to inspect deeply,
# in addition to the true anchor.
TOP_ANCHORS_TO_SEARCH = 12


# ------------------------------------------------------------
# Primality / semiprime generation
# ------------------------------------------------------------

def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37,
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    bases = [2, 3, 5, 7, 11, 13, 17]

    for a in bases:
        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(bits: int, rng: random.Random) -> int:
    while True:
        x = rng.getrandbits(bits)
        x |= 1 << (bits - 1)
        x |= 1

        if is_probable_prime(x):
            return x


def make_semiprime(
    bits: int,
    rng: random.Random,
) -> tuple[int, int, int]:

    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:
        p = random_prime(p_bits, rng)
        q = random_prime(q_bits, rng)

        if p != q:
            return p, q, p * q


# ------------------------------------------------------------
# CRT
# ------------------------------------------------------------

def crt_pair(
    a: int,
    m: int,
    b: int,
    n: int,
) -> int:

    t = ((b - a) * pow(m, -1, n)) % n
    return a + m * t


def crt_vector(
    residues: list[int],
    moduli: list[int],
) -> int:

    x = 0
    M = 1

    for a, m in zip(residues, moduli):

        t = ((a - x) * pow(M, -1, m)) % m

        x += M * t
        M *= m

    return x


# ------------------------------------------------------------
# Exact C
# ------------------------------------------------------------

def exact_C_from_ab(
    n: int,
    a: int,
    b: int,
    r: int,
    s: int,
) -> int:

    D = n - a * b

    beta = (D * pow(r, -1, s)) % s
    alpha = (D * pow(s, -1, r)) % r

    return (r * beta + s * alpha + a * b) // (r * s)


# ------------------------------------------------------------
# Build the product relation
# ------------------------------------------------------------

def build_product_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
) -> tuple[set[int], list[set[int]], list[set[int]]]:
    """
    Returns:

        product_set
            z values satisfying F(z)=target_c

        support_a[a]
            b values compatible with a

        support_b[b]
            a values compatible with b
    """

    product_set: set[int] = set()

    support_a = [
        set()
        for _ in range(r)
    ]

    support_b = [
        set()
        for _ in range(s)
    ]

    for a in range(1, r):

        for b in range(1, s):

            c = exact_C_from_ab(
                n,
                a,
                b,
                r,
                s,
            )

            if c != target_c:
                continue

            z = a * b

            product_set.add(z)
            support_a[a].add(b)
            support_b[b].add(a)

    return (
        product_set,
        support_a,
        support_b,
    )


# ------------------------------------------------------------
# Build the complete grid
# ------------------------------------------------------------

def build_grid(
    n: int,
    p: int,
    q: int,
):
    """
    Grid entry:

        grid[i][j] = {
            r,
            s,
            observed C,
            allowed products Z,
            support tables
        }
    """

    grid = []

    for i, r in enumerate(R1):

        row = []

        for j, s in enumerate(R2):

            a_true = p % r
            b_true = q % s

            true_c = exact_C_from_ab(
                n,
                a_true,
                b_true,
                r,
                s,
            )

            (
                product_set,
                support_a,
                support_b,
            ) = build_product_relation(
                n,
                r,
                s,
                true_c,
            )

            row.append(
                {
                    "r": r,
                    "s": s,
                    "C": true_c,
                    "products": product_set,
                    "support_a": support_a,
                    "support_b": support_b,
                    "true_a": a_true,
                    "true_b": b_true,
                }
            )

        grid.append(row)

    return grid


# ------------------------------------------------------------
# Initial domains from an anchor
# ------------------------------------------------------------

def anchored_domains(
    grid,
    anchor_a: int,
    anchor_b: int,
):
    """
    Fix:

        A_0 = {anchor_a}
        B_0 = {anchor_b}

    and immediately propagate through row 0
    and column 0.
    """

    A = [
        set(range(1, r))
        for r in R1
    ]

    B = [
        set(range(1, s))
        for s in R2
    ]

    # Anchor.
    A[0] = {anchor_a}
    B[0] = {anchor_b}

    # Column 0 -> all A_i.
    for i in range(len(R1)):

        relation = grid[i][0]

        if i == 0:
            continue

        new_domain = set()

        for a in A[i]:

            if anchor_b in relation["support_a"][a]:
                new_domain.add(a)

        A[i] = new_domain

        if not A[i]:
            return None

    # Row 0 -> all B_j.
    for j in range(len(R2)):

        relation = grid[0][j]

        if j == 0:
            continue

        new_domain = set()

        for b in B[j]:

            if anchor_a in relation["support_b"][b]:
                new_domain.add(b)

        B[j] = new_domain

        if not B[j]:
            return None

    return A, B


# ------------------------------------------------------------
# Arc-consistency propagation
# ------------------------------------------------------------

def propagate(
    grid,
    A,
    B,
) -> tuple[list[set[int]], list[set[int]]] | None:

    queue = deque()

    for i in range(len(R1)):
        for j in range(len(R2)):
            queue.append(("A", i, j))
            queue.append(("B", i, j))

    while queue:

        side, i, j = queue.popleft()

        relation = grid[i][j]

        if side == "A":

            old_domain = A[i]

            new_domain = set()

            for a in old_domain:

                possible_b = relation["support_a"][a]

                if possible_b & B[j]:
                    new_domain.add(a)

            if new_domain == old_domain:
                continue

            if not new_domain:
                return None

            A[i] = new_domain

            # A_i changed, so every B neighbor must be checked.
            for jj in range(len(R2)):
                if jj != j:
                    queue.append(("B", i, jj))

        else:

            old_domain = B[j]

            new_domain = set()

            for b in old_domain:

                possible_a = relation["support_b"][b]

                if possible_a & A[i]:
                    new_domain.add(b)

            if new_domain == old_domain:
                continue

            if not new_domain:
                return None

            B[j] = new_domain

            # B_j changed, so every A neighbor must be checked.
            for ii in range(len(R1)):
                if ii != i:
                    queue.append(("A", ii, j))

    return A, B


# ------------------------------------------------------------
# Domain signature
# ------------------------------------------------------------

def domain_signature(A, B):
    widths_a = [len(x) for x in A]
    widths_b = [len(x) for x in B]

    total = sum(widths_a) + sum(widths_b)

    product_space = 1

    for x in A:
        product_space *= len(x)

    for x in B:
        product_space *= len(x)

    return {
        "widths_a": widths_a,
        "widths_b": widths_b,
        "sum_widths": total,
        "cartesian_space": product_space,
    }


# ------------------------------------------------------------
# Full assignment verification
# ------------------------------------------------------------

def assignment_is_valid(
    grid,
    A_values,
    B_values,
) -> bool:

    for i in range(len(R1)):
        for j in range(len(R2)):

            a = A_values[i]
            b = B_values[j]

            if (
                exact_C_from_ab(
                    CURRENT_N,
                    a,
                    b,
                    R1[i],
                    R2[j],
                )
                != grid[i][j]["C"]
            ):
                return False

    return True


# ------------------------------------------------------------
# DFS after propagation
# ------------------------------------------------------------

def search_assignments(
    grid,
    A,
    B,
    true_A,
    true_B,
    cap: int,
):
    """
    Enumerate compatible full assignments after AC propagation.

    Uses MRV.

    Returns:
        number_found
        cap_hit
        true_found
    """

    solutions = 0
    cap_hit = False
    true_found = False

    def recurse(A_cur, B_cur):

        nonlocal solutions
        nonlocal cap_hit
        nonlocal true_found

        if solutions >= cap:
            cap_hit = True
            return

        propagated = propagate(
            grid,
            [set(x) for x in A_cur],
            [set(x) for x in B_cur],
        )

        if propagated is None:
            return

        A2, B2 = propagated

        # Fully determined?
        complete = (
            all(len(x) == 1 for x in A2)
            and all(len(x) == 1 for x in B2)
        )

        if complete:

            av = [next(iter(x)) for x in A2]
            bv = [next(iter(x)) for x in B2]

            solutions += 1

            if av == true_A and bv == true_B:
                true_found = True

            return

        # MRV.
        best_side = None
        best_index = None
        best_values = None
        best_size = 10**9

        for i, domain in enumerate(A2):

            if 1 < len(domain) < best_size:
                best_side = "A"
                best_index = i
                best_values = sorted(domain)
                best_size = len(domain)

        for j, domain in enumerate(B2):

            if 1 < len(domain) < best_size:
                best_side = "B"
                best_index = j
                best_values = sorted(domain)
                best_size = len(domain)

        if best_values is None:
            return

        for value in best_values:

            A3 = [set(x) for x in A2]
            B3 = [set(x) for x in B2]

            if best_side == "A":
                A3[best_index] = {value}
            else:
                B3[best_index] = {value}

            recurse(A3, B3)

            if solutions >= cap:
                cap_hit = True
                return

    recurse(A, B)

    return solutions, cap_hit, true_found


# ------------------------------------------------------------
# Analyze one sample
# ------------------------------------------------------------

def analyze_sample(
    n: int,
    p: int,
    q: int,
):

    global CURRENT_N
    CURRENT_N = n

    grid = build_grid(n, p, q)

    true_A = [p % r for r in R1]
    true_B = [q % s for s in R2]

    true_anchor = (
        true_A[0],
        true_B[0],
    )

    anchors = []

    # --------------------------------------------------------
    # All corner-compatible anchors
    # --------------------------------------------------------

    corner = grid[0][0]

    corner_anchors = []

    for a in range(1, R1[0]):
        if not corner["support_a"][a]:
            continue

        for b in corner["support_a"][a]:

            if b in corner["support_a"][a]:
                corner_anchors.append((a, b))

    # --------------------------------------------------------
    # Test every anchor with row/column filtering + AC-3
    # --------------------------------------------------------

    start = time.perf_counter()

    feasible_anchors = []

    for anchor_a, anchor_b in corner_anchors:

        result = anchored_domains(
            grid,
            anchor_a,
            anchor_b,
        )

        if result is None:
            continue

        A, B = result

        propagated = propagate(
            grid,
            A,
            B,
        )

        if propagated is None:
            continue

        A, B = propagated

        sig = domain_signature(A, B)

        is_true_anchor = (
            anchor_a == true_anchor[0]
            and anchor_b == true_anchor[1]
        )

        feasible_anchors.append(
            (
                sig["sum_widths"],
                sig["cartesian_space"],
                anchor_a,
                anchor_b,
                A,
                B,
                is_true_anchor,
            )
        )

    elapsed = time.perf_counter() - start

    # Sort by total remaining domain width.
    feasible_anchors.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    # --------------------------------------------------------
    # Deep-search the most constrained anchors.
    # Always include the true anchor.
    # --------------------------------------------------------

    selected = []

    seen = set()

    for item in feasible_anchors:

        key = (item[2], item[3])

        if key not in seen:
            selected.append(item)
            seen.add(key)

        if len(selected) >= TOP_ANCHORS_TO_SEARCH:
            break

    for item in feasible_anchors:

        if item[6]:

            key = (item[2], item[3])

            if key not in seen:
                selected.append(item)
                seen.add(key)

            break

    # --------------------------------------------------------
    # Deep search
    # --------------------------------------------------------

    deep_results = []

    for item in selected:

        (
            sum_widths,
            cartesian_space,
            anchor_a,
            anchor_b,
            A,
            B,
            is_true_anchor,
        ) = item

        solutions, cap_hit, true_found = \
            search_assignments(
                grid,
                A,
                B,
                true_A,
                true_B,
                DFS_SOLUTION_CAP,
            )

        deep_results.append(
            {
                "anchor": (anchor_a, anchor_b),
                "sum_widths": sum_widths,
                "cartesian_space": cartesian_space,
                "solutions": solutions,
                "cap_hit": cap_hit,
                "true_anchor": is_true_anchor,
                "true_found": true_found,
                "A": A,
                "B": B,
            }
        )

    # --------------------------------------------------------
    # Locate true anchor result
    # --------------------------------------------------------

    true_result = None

    for item in deep_results:
        if item["true_anchor"]:
            true_result = item
            break

    # --------------------------------------------------------
    # Exact CRT reconstruction
    # --------------------------------------------------------

    recovered = []

    if true_result is not None:

        A = true_result["A"]
        B = true_result["B"]

        if (
            all(len(x) == 1 for x in A)
            and all(len(x) == 1 for x in B)
        ):

            rec_A = [next(iter(x)) for x in A]
            rec_B = [next(iter(x)) for x in B]

            recovered_p = crt_vector(
                rec_A,
                R1,
            )

            recovered_q = crt_vector(
                rec_B,
                R2,
            )

            recovered = [
                recovered_p == p,
                recovered_q == q,
                recovered_p,
                recovered_q,
            ]

    return {
        "grid": grid,
        "true_A": true_A,
        "true_B": true_B,
        "true_anchor": true_anchor,
        "corner_anchors": len(corner_anchors),
        "feasible_anchors": feasible_anchors,
        "deep_results": deep_results,
        "true_result": true_result,
        "elapsed": elapsed,
        "recovered": recovered,
    }


# ------------------------------------------------------------
# Display
# ------------------------------------------------------------

def print_sample_result(
    bits: int,
    sample: int,
    p: int,
    q: int,
    n: int,
    result,
):

    feasible = result["feasible_anchors"]
    deep = result["deep_results"]

    print()
    print(f"Sample {sample}:")
    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  n = {n}")
    print()

    print("True residue vectors:")

    print(
        "  A = ["
        + ", ".join(
            str(x)
            for x in result["true_A"]
        )
        + "]"
    )

    print(
        "  B = ["
        + ", ".join(
            str(x)
            for x in result["true_B"]
        )
        + "]"
    )

    print(
        f"  true anchor = {result['true_anchor']}"
    )

    print()
    print("Anchor filtering:")
    print(
        f"  corner-compatible anchors = "
        f"{result['corner_anchors']}"
    )
    print(
        f"  globally feasible anchors = "
        f"{len(feasible)}"
    )

    print(
        f"  anchor propagation time = "
        f"{result['elapsed']:.6f}s"
    )

    print()
    print("Most constrained surviving anchors:")

    for item in deep:

        print(
            f"  anchor={item['anchor']} "
            f"sum_width={item['sum_widths']:3d} "
            f"cartesian={item['cartesian_space']}"
        )

        print(
            f"      A="
            f"{[sorted(x) for x in item['A']]}"
        )

        print(
            f"      B="
            f"{[sorted(x) for x in item['B']]}"
        )

        print(
            f"      DFS solutions={item['solutions']} "
            f"cap_hit={item['cap_hit']} "
            f"true_found={item['true_found']}"
        )

    print()

    if result["true_result"] is None:

        print(
            "TRUE ANCHOR: not included in deep-search set"
        )

    else:

        tr = result["true_result"]

        print("TRUE ANCHOR RESULT:")

        print(
            f"  anchor              = "
            f"{tr['anchor']}"
        )

        print(
            f"  remaining sum width = "
            f"{tr['sum_widths']}"
        )

        print(
            f"  remaining Cartesian "
            f"space               = "
            f"{tr['cartesian_space']}"
        )

        print(
            f"  DFS solutions       = "
            f"{tr['solutions']}"
        )

        print(
            f"  DFS cap hit        = "
            f"{tr['cap_hit']}"
        )

        print(
            f"  true assignment found = "
            f"{tr['true_found']}"
        )

    if result["recovered"]:

        ok_p, ok_q, rec_p, rec_q = result["recovered"]

        print()
        print("CRT RECONSTRUCTION:")
        print(f"  recovered p = {rec_p}")
        print(f"  recovered q = {rec_q}")
        print(f"  p exact     = {ok_p}")
        print(f"  q exact     = {ok_q}")


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 163")
    print("=" * 72)

    print(
        "Global Multiplicative Rank-1 Recovery "
        "from C-Product Sets"
    )

    print()
    print("Experiment 162 established:")
    print()
    print("    C_ij = F_ij(a_i * b_j)")
    print()
    print("with")
    print()
    print("    a_i = p mod r_i")
    print("    b_j = q mod s_j")
    print()
    print("Therefore every cell produces an allowed")
    print("product set Z_ij, while globally")
    print()
    print("    z_ij = a_i * b_j")
    print()
    print("This experiment anchors one (a_0,b_0),")
    print("propagates divisibility/product constraints,")
    print("then applies full bipartite arc consistency.")
    print("=" * 72)

    global_best = []

    for bits in BIT_SIZES:

        print()
        print("-" * 72)
        print(f"BIT SIZE = {bits}")
        print("-" * 72)

        for sample in range(1, SAMPLES_PER_SIZE + 1):

            p, q, n = make_semiprime(bits, rng)

            result = analyze_sample(
                n,
                p,
                q,
            )

            print_sample_result(
                bits,
                sample,
                p,
                q,
                n,
                result,
            )

            feasible = result["feasible_anchors"]

            if feasible:

                best = feasible[0]

                global_best.append(
                    (
                        bits,
                        best[0],
                        best[1],
                        best[2],
                        best[3],
                    )
                )

    # --------------------------------------------------------
    # Final summary
    # --------------------------------------------------------

    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    print()

    for bits, width, space, a, b in global_best:

        print(
            f"{bits:2d}-bit: "
            f"best anchor=({a},{b}), "
            f"sum_width={width}, "
            f"Cartesian={space}"
        )

    print()
    print("Interpretation:")
    print()
    print("1. Corner anchoring supplies candidate A_i and B_j domains.")
    print("2. Cross-cell product constraints are then propagated.")
    print("3. The key measurement is whether the true anchor collapses")
    print("   to singleton residue vectors.")
    print("4. A singleton 5x5 residue vector can be CRT-reconstructed.")
    print()
    print("M1 = product(R1) =", math.prod(R1))
    print("M2 = product(R2) =", math.prod(R2))
    print()
    print("If the recovered CRT residues are smaller than their")
    print("corresponding moduli products, they directly give p and q.")
    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 163")
    print("=" * 72)


if __name__ == "__main__":
    main()
