#!/usr/bin/env python3

# ============================================================
# START EXPERIMENT 164
# CRT Interval Lift After Multiplicative C-Propagation
# ============================================================

from __future__ import annotations

import itertools
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

SEED = 164

# Maximum number of explicit CRT branch combinations.
# Larger searches are reported as capped.
MAX_BRANCHES = 5_000_000


# ------------------------------------------------------------
# Primality
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

    for a in [2, 3, 5, 7, 11, 13, 17]:
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

    p = random_prime(bits // 2, rng)
    q = random_prime(bits - bits // 2, rng)

    while p == q:
        q = random_prime(bits - bits // 2, rng)

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


def crt_assignments(
    assignments: list[tuple[int, int]],
) -> tuple[int, int]:

    """
    CRT for a list [(residue, modulus), ...].

    Returns:
        x, M
    """

    x = 0
    M = 1

    for a, m in assignments:

        t = ((a - x) * pow(M, -1, m)) % m

        x += M * t
        M *= m

    return x, M


# ------------------------------------------------------------
# C formula
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

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)


# ------------------------------------------------------------
# Build one C-product relation
# ------------------------------------------------------------

def build_relation(
    n: int,
    r: int,
    s: int,
    target_c: int,
):

    support_a = [
        set()
        for _ in range(r)
    ]

    support_b = [
        set()
        for _ in range(s)
    ]

    products = set()

    for a in range(1, r):

        for b in range(1, s):

            if exact_C_from_ab(
                n,
                a,
                b,
                r,
                s,
            ) != target_c:
                continue

            support_a[a].add(b)
            support_b[b].add(a)

            products.add(a * b)

    return {
        "r": r,
        "s": s,
        "C": target_c,
        "support_a": support_a,
        "support_b": support_b,
        "products": products,
    }


# ------------------------------------------------------------
# Construct entire grid
# ------------------------------------------------------------

def build_grid(
    n: int,
    p: int,
    q: int,
):

    grid = []

    for r in R1:

        row = []

        for s in R2:

            a = p % r
            b = q % s

            c = exact_C_from_ab(
                n,
                a,
                b,
                r,
                s,
            )

            row.append(
                build_relation(
                    n,
                    r,
                    s,
                    c,
                )
            )

        grid.append(row)

    return grid


# ------------------------------------------------------------
# AC propagation
# ------------------------------------------------------------

def propagate(
    grid,
    A,
    B,
):

    queue = deque()

    for i in range(len(R1)):
        for j in range(len(R2)):
            queue.append(("A", i, j))
            queue.append(("B", i, j))

    while queue:

        side, i, j = queue.popleft()

        rel = grid[i][j]

        if side == "A":

            old = A[i]

            new = {
                a
                for a in old
                if rel["support_a"][a] & B[j]
            }

            if not new:
                return None

            if new == old:
                continue

            A[i] = new

            for jj in range(len(R2)):
                queue.append(("B", i, jj))

        else:

            old = B[j]

            new = {
                b
                for b in old
                if rel["support_b"][b] & A[i]
            }

            if not new:
                return None

            if new == old:
                continue

            B[j] = new

            for ii in range(len(R1)):
                queue.append(("A", ii, j))

    return A, B


# ------------------------------------------------------------
# Anchor initialization
# ------------------------------------------------------------

def anchored_domains(
    grid,
    anchor_a: int,
    anchor_b: int,
):

    A = [
        set(range(1, r))
        for r in R1
    ]

    B = [
        set(range(1, s))
        for s in R2
    ]

    A[0] = {anchor_a}
    B[0] = {anchor_b}

    result = propagate(
        grid,
        A,
        B,
    )

    return result


# ------------------------------------------------------------
# CRT subset selection
# ------------------------------------------------------------

def find_best_subset(
    moduli: list[int],
    domains: list[set[int]],
    upper: int,
):
    """
    Find the subset of residue variables which:

        product(moduli) > upper

    while minimizing the Cartesian branch count

        product(len(domain_i)).

    Since there are only five variables, exhaustive subset
    search is trivial.
    """

    best = None

    n = len(moduli)

    for mask in range(1, 1 << n):

        indices = [
            i
            for i in range(n)
            if mask & (1 << i)
        ]

        M = 1
        branches = 1

        for i in indices:
            M *= moduli[i]
            branches *= len(domains[i])

        if M <= upper:
            continue

        candidate = (
            branches,
            len(indices),
            -M,
            indices,
            M,
        )

        if best is None or candidate < best:
            best = candidate

    return best


# ------------------------------------------------------------
# Enumerate CRT representatives
# ------------------------------------------------------------

def enumerate_crt_candidates(
    moduli: list[int],
    domains: list[set[int]],
    upper: int,
    remaining_domains: list[set[int]],
):
    """
    Choose a subset whose modulus product exceeds 'upper'.

    Enumerate assignments only on that subset.

    Because M > upper, each CRT residue corresponds to at most
    one positive integer <= upper.

    Then verify all remaining residue constraints.
    """

    subset_info = find_best_subset(
        moduli,
        domains,
        upper,
    )

    if subset_info is None:
        return {
            "status": "NO_SUBSET",
            "branches": 0,
            "candidates": [],
            "subset": None,
            "M": None,
        }

    branches, _, _, indices, M = subset_info

    if branches > MAX_BRANCHES:

        return {
            "status": "CAP",
            "branches": branches,
            "candidates": [],
            "subset": indices,
            "M": M,
        }

    value_lists = [
        sorted(domains[i])
        for i in indices
    ]

    candidates = []

    tested = 0

    for values in itertools.product(*value_lists):

        tested += 1

        assignments = [
            (value, moduli[i])
            for value, i in zip(values, indices)
        ]

        x, crt_M = crt_assignments(assignments)

        # x is the canonical representative [0,M).
        if x < 2 or x > upper:
            continue

        # Verify all residue domains.
        valid = True

        for i, m in enumerate(moduli):

            if x % m not in domains[i]:
                valid = False
                break

        if not valid:
            continue

        candidates.append(x)

        if len(candidates) > 1000:
            break

    return {
        "status": "OK",
        "branches": branches,
        "tested": tested,
        "candidates": candidates,
        "subset": indices,
        "M": M,
    }


# ------------------------------------------------------------
# Try one side
# ------------------------------------------------------------

def solve_side(
    name: str,
    moduli: list[int],
    domains: list[set[int]],
    upper: int,
    other_moduli: list[int],
    other_domains: list[set[int]],
    n: int,
):

    result = enumerate_crt_candidates(
        moduli,
        domains,
        upper,
        other_domains,
    )

    if result["status"] != "OK":
        return {
            "name": name,
            "result": result,
            "factor_candidates": [],
            "exact_factors": [],
        }

    candidates = result["candidates"]

    exact = []

    for x in candidates:

        if n % x != 0:
            continue

        y = n // x

        if y < 2:
            continue

        # Verify the other-side residue vector.
        good_other = True

        for m, domain in zip(
            other_moduli,
            other_domains,
        ):

            if y % m not in domain:
                good_other = False
                break

        if good_other:
            exact.append(
                (x, y)
            )

    return {
        "name": name,
        "result": result,
        "factor_candidates": candidates,
        "exact_factors": exact,
    }


# ------------------------------------------------------------
# Analyze one anchored state
# ------------------------------------------------------------

def analyze_anchor(
    grid,
    n: int,
    p: int,
    q: int,
    anchor_a: int,
    anchor_b: int,
):

    domains = anchored_domains(
        grid,
        anchor_a,
        anchor_b,
    )

    if domains is None:
        return None

    A, B = domains

    sqrt_n = math.isqrt(n)

    # We assume p <= q for solving.
    # The generated primes are approximately balanced, and
    # exact verification below handles the actual factor.
    upper_p = sqrt_n
    upper_q = n // 2

    A_domains = A
    B_domains = B

    a_solution = solve_side(
        "A/p",
        R1,
        A_domains,
        upper_p,
        R2,
        B_domains,
        n,
    )

    b_solution = solve_side(
        "B/q",
        R2,
        B_domains,
        upper_p,
        R1,
        A_domains,
        n,
    )

    true_anchor = (
        anchor_a == p % R1[0]
        and anchor_b == q % R2[0]
    )

    return {
        "A": A,
        "B": B,
        "true_anchor": true_anchor,
        "A_solution": a_solution,
        "B_solution": b_solution,
    }


# ------------------------------------------------------------
# Find corner anchors
# ------------------------------------------------------------

def get_corner_anchors(
    grid,
):

    rel = grid[0][0]

    anchors = []

    for a in range(1, R1[0]):

        for b in rel["support_a"][a]:

            anchors.append(
                (a, b)
            )

    return anchors


# ------------------------------------------------------------
# Main sample
# ------------------------------------------------------------

def analyze_sample(
    p: int,
    q: int,
    n: int,
):

    grid = build_grid(
        n,
        p,
        q,
    )

    anchors = get_corner_anchors(grid)

    true_anchor = (
        p % R1[0],
        q % R2[0],
    )

    print()
    print(
        f"true p={p}"
    )
    print(
        f"true q={q}"
    )
    print(
        f"n={n}"
    )

    print(
        f"sqrt(n)={math.isqrt(n)}"
    )

    print(
        f"corner anchors={len(anchors)}"
    )

    # --------------------------------------------------------
    # First cheaply rank anchors by AC domain volume.
    # --------------------------------------------------------

    ranked = []

    for anchor in anchors:

        result = anchored_domains(
            grid,
            anchor[0],
            anchor[1],
        )

        if result is None:
            continue

        A, B = result

        total_width = (
            sum(len(x) for x in A)
            + sum(len(x) for x in B)
        )

        # Estimate cheapest CRT bridge on both sides.
        a_subset = find_best_subset(
            R1,
            A,
            math.isqrt(n),
        )

        b_subset = find_best_subset(
            R2,
            B,
            math.isqrt(n),
        )

        a_cost = (
            a_subset[0]
            if a_subset is not None
            else 10**100
        )

        b_cost = (
            b_subset[0]
            if b_subset is not None
            else 10**100
        )

        crt_cost = min(
            a_cost,
            b_cost,
        )

        ranked.append(
            (
                crt_cost,
                total_width,
                anchor,
                A,
                B,
            )
        )

    ranked.sort(
        key=lambda x: (
            x[0],
            x[1],
        )
    )

    print()
    print("Best anchors by CRT bridge cost:")

    # Always include true anchor.
    selected = []
    seen = set()

    for item in ranked[:12]:

        if item[2] not in seen:
            selected.append(item)
            seen.add(item[2])

    for item in ranked:

        if item[2] == true_anchor:
            if item[2] not in seen:
                selected.append(item)
                seen.add(item[2])
            break

    for item in selected:

        crt_cost, width, anchor, A, B = item

        a_sub = find_best_subset(
            R1,
            A,
            math.isqrt(n),
        )

        b_sub = find_best_subset(
            R2,
            B,
            math.isqrt(n),
        )

        print(
            f"  anchor={anchor} "
            f"CRTcost={crt_cost:,} "
            f"width={width}"
        )

        if a_sub is not None:

            print(
                f"      A subset="
                f"{a_sub[3]} "
                f"M={a_sub[4]:,} "
                f"branches={a_sub[0]:,}"
            )

        if b_sub is not None:

            print(
                f"      B subset="
                f"{b_sub[3]} "
                f"M={b_sub[4]:,} "
                f"branches={b_sub[0]:,}"
            )

    # --------------------------------------------------------
    # Actually solve selected anchors.
    # --------------------------------------------------------

    print()
    print("CRT lifting:")

    total_exact = 0
    true_recovered = False

    for item in selected:

        anchor = item[2]

        t0 = time.perf_counter()

        result = analyze_anchor(
            grid,
            n,
            p,
            q,
            anchor[0],
            anchor[1],
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        if result is None:
            print(
                f"  anchor={anchor}: "
                f"AC contradiction"
            )
            continue

        ares = result["A_solution"]
        bres = result["B_solution"]

        print()
        print(
            f"  anchor={anchor} "
            f"time={elapsed:.6f}s"
        )

        for side_result in (
            ares,
            bres,
        ):

            side = side_result["name"]
            rr = side_result["result"]

            print(
                f"    {side}: "
                f"status={rr['status']} "
                f"branches={rr.get('branches', 0):,} "
                f"candidates="
                f"{len(side_result['factor_candidates'])}"
            )

            if rr.get("subset") is not None:
                print(
                    f"        subset="
                    f"{rr['subset']} "
                    f"M={rr['M']:,}"
                )

            if side_result["exact_factors"]:

                print(
                    f"        EXACT="
                    f"{side_result['exact_factors']}"
                )

                total_exact += len(
                    side_result["exact_factors"]
                )

                for pair in side_result["exact_factors"]:

                    if pair == (p, q):
                        true_recovered = True

        if result["true_anchor"]:

            print(
                "    TRUE ANCHOR"
            )

            print(
                f"      true p={p}"
            )

            print(
                f"      true q={q}"
            )

            print(
                f"      recovered="
                f"{true_recovered}"
            )

    return {
        "anchors": len(anchors),
        "true_recovered": true_recovered,
        "exact_total": total_exact,
    }


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 164")
    print("=" * 72)
    print(
        "CRT Interval Lift After "
        "Multiplicative C-Propagation"
    )
    print()
    print("Experiment 163 showed that C-product propagation")
    print("does not determine the residue vectors by itself.")
    print()
    print("But A=[p mod r_i] and B=[q mod s_j] are not")
    print("independent variables: each complete vector must")
    print("be the CRT representation of one actual integer.")
    print()
    print("Since p < sqrt(n), a CRT modulus M > sqrt(n)")
    print("makes a residue vector correspond to at most one p.")
    print()
    print("This experiment searches only the cheapest subset")
    print("of residue variables whose modulus product exceeds")
    print("sqrt(n), then checks the resulting integer against")
    print("all remaining residue domains and n = p*q.")
    print("=" * 72)

    summaries = []

    for bits in BIT_SIZES:

        print()
        print("-" * 72)
        print(f"BIT SIZE = {bits}")
        print("-" * 72)

        for sample in range(
            1,
            SAMPLES_PER_SIZE + 1,
        ):

            p, q, n = make_semiprime(
                bits,
                rng,
            )

            result = analyze_sample(
                p,
                q,
                n,
            )

            summaries.append(
                (
                    bits,
                    result["true_recovered"],
                    result["exact_total"],
                )
            )

    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    for bits, recovered, exact_total in summaries:

        print(
            f"{bits:2d}-bit: "
            f"true_recovered={recovered} "
            f"exact_hits={exact_total}"
        )

    print()
    print("FINISHED EXPERIMENT 164")


if __name__ == "__main__":
    main()
