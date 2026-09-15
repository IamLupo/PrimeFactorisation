#!/usr/bin/env python3

# ============================================================
# START EXPERIMENT 166
# Global CRT Candidate Intersection
# ============================================================

from __future__ import annotations

import itertools
import math
import random
import time
from collections import defaultdict, deque


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48, 54]
SAMPLES_PER_SIZE = 1

SEED = 166


# ------------------------------------------------------------
# Primality
# ------------------------------------------------------------

def is_probable_prime(n: int) -> bool:

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37,
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in (
        2, 3, 5, 7, 11, 13, 17
    ):

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = x * x % n

            if x == n - 1:
                break

        else:
            return False

    return True


def random_prime(
    bits: int,
    rng: random.Random,
) -> int:

    while True:

        x = rng.getrandbits(bits)

        x |= 1 << (bits - 1)
        x |= 1

        if is_probable_prime(x):
            return x


def make_semiprime(
    bits: int,
    rng: random.Random,
):

    while True:

        p = random_prime(bits // 2, rng)
        q = random_prime(bits - bits // 2, rng)

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
):

    t = ((b - a) * pow(m, -1, n)) % n

    return a + m * t


def crt_list(
    residues,
    moduli,
):

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli,
    ):

        t = ((a - x) * pow(M, -1, m)) % m

        x += M * t
        M *= m

    return x, M


# ------------------------------------------------------------
# C
# ------------------------------------------------------------

def C_value(
    n: int,
    a: int,
    b: int,
    r: int,
    s: int,
) -> int:

    D = n - a * b

    beta = (
        D
        * pow(r, -1, s)
    ) % s

    alpha = (
        D
        * pow(s, -1, r)
    ) % r

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)


# ------------------------------------------------------------
# Build cell
# ------------------------------------------------------------

def build_cell(
    n: int,
    r: int,
    s: int,
    C: int,
):

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

            if C_value(
                n,
                a,
                b,
                r,
                s,
            ) != C:
                continue

            support_a[a].add(b)
            support_b[b].add(a)

    return {
        "r": r,
        "s": s,
        "C": C,
        "support_a": support_a,
        "support_b": support_b,
    }


# ------------------------------------------------------------
# Build grid
# ------------------------------------------------------------

def build_grid(
    n: int,
    p: int,
    q: int,
):

    grid = []

    for r in R1:

        row = []

        a = p % r

        for s in R2:

            b = q % s

            c = C_value(
                n,
                a,
                b,
                r,
                s,
            )

            row.append(
                build_cell(
                    n,
                    r,
                    s,
                    c,
                )
            )

        grid.append(row)

    return grid


# ------------------------------------------------------------
# AC propagation for one anchor
# ------------------------------------------------------------

def propagate(
    grid,
    A,
    B,
):

    queue = deque()

    for i in range(5):
        for j in range(5):

            queue.append(
                ("A", i, j)
            )

            queue.append(
                ("B", i, j)
            )

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

            if new != old:

                A[i] = new

                for jj in range(5):
                    queue.append(
                        ("B", i, jj)
                    )

        else:

            old = B[j]

            new = {
                b
                for b in old
                if rel["support_b"][b] & A[i]
            }

            if not new:
                return None

            if new != old:

                B[j] = new

                for ii in range(5):
                    queue.append(
                        ("A", ii, j)
                    )

    return A, B


# ------------------------------------------------------------
# Enumerate corner anchors and propagated domains
# ------------------------------------------------------------

def all_anchor_states(
    n: int,
    grid,
):

    corner = grid[0][0]

    states = []

    for a in range(1, R1[0]):

        for b in corner["support_a"][a]:

            A = [
                set(range(1, r))
                for r in R1
            ]

            B = [
                set(range(1, s))
                for s in R2
            ]

            A[0] = {a}
            B[0] = {b}

            result = propagate(
                grid,
                A,
                B,
            )

            if result is None:
                continue

            A, B = result

            states.append(
                {
                    "anchor": (a, b),
                    "A": A,
                    "B": B,
                }
            )

    return states


# ------------------------------------------------------------
# Choose a CRT subset
# ------------------------------------------------------------

def cheapest_subset(
    moduli,
    domains,
    upper,
):

    best = None

    for mask in range(
        1,
        1 << len(moduli),
    ):

        idx = [
            i
            for i in range(len(moduli))
            if mask & (1 << i)
        ]

        M = 1
        branches = 1

        for i in idx:

            M *= moduli[i]
            branches *= len(domains[i])

        if M <= upper:
            continue

        candidate = (
            branches,
            -M,
            idx,
            M,
        )

        if best is None or candidate < best:
            best = candidate

    return best


# ------------------------------------------------------------
# Build direct candidate lookup
# ------------------------------------------------------------

def build_candidate_index(
    n: int,
    moduli,
    domains,
    upper,
):

    info = cheapest_subset(
        moduli,
        domains,
        upper,
    )

    if info is None:
        return {
            "status": "NO_SUBSET",
            "candidates": set(),
            "branches": 0,
            "M": None,
            "indices": None,
        }

    branches, neg_M, indices, M = info

    value_lists = [
        sorted(domains[i])
        for i in indices
    ]

    candidates = set()

    for values in itertools.product(
        *value_lists
    ):

        residues = [
            v
            for v in values
        ]

        subset_moduli = [
            moduli[i]
            for i in indices
        ]

        x, crt_M = crt_list(
            residues,
            subset_moduli,
        )

        if x < 2:
            continue

        if x > upper:
            continue

        # Verify every currently available residue.
        valid = True

        for m, domain in zip(
            moduli,
            domains,
        ):

            if x % m not in domain:
                valid = False
                break

        if valid:
            candidates.add(x)

    return {
        "status": "OK",
        "candidates": candidates,
        "branches": branches,
        "M": M,
        "indices": indices,
    }


# ------------------------------------------------------------
# Global candidate intersection
# ------------------------------------------------------------

def analyze_sample(
    p: int,
    q: int,
    n: int,
):

    sqrt_n = math.isqrt(n)

    t0 = time.perf_counter()

    grid = build_grid(
        n,
        p,
        q,
    )

    states = all_anchor_states(
        n,
        grid,
    )

    anchor_time = (
        time.perf_counter()
        - t0
    )

    true_anchor = (
        p % R1[0],
        q % R2[0],
    )

    true_state = None

    for state in states:

        if state["anchor"] == true_anchor:
            true_state = state
            break

    print()
    print(
        f"p={p}"
    )
    print(
        f"q={q}"
    )
    print(
        f"n={n}"
    )
    print(
        f"sqrt(n)={sqrt_n}"
    )
    print(
        f"corner anchors={len(states)}"
    )
    print(
        f"anchor propagation time="
        f"{anchor_time:.6f}s"
    )

    # --------------------------------------------------------
    # Global candidate sets.
    #
    # Rather than solving each anchor independently, collect
    # every possible small factor generated by every anchor.
    # --------------------------------------------------------

    global_A = set()
    global_B = set()

    cheapest_A = None
    cheapest_B = None

    for state in states:

        A = state["A"]
        B = state["B"]

        info_A = build_candidate_index(
            n,
            R1,
            A,
            sqrt_n,
        )

        info_B = build_candidate_index(
            n,
            R2,
            B,
            sqrt_n,
        )

        if info_A["status"] == "OK":

            global_A.update(
                info_A["candidates"]
            )

            if cheapest_A is None:
                cheapest_A = info_A
            else:
                if (
                    info_A["branches"]
                    < cheapest_A["branches"]
                ):
                    cheapest_A = info_A

        if info_B["status"] == "OK":

            global_B.update(
                info_B["candidates"]
            )

            if cheapest_B is None:
                cheapest_B = info_B
            else:
                if (
                    info_B["branches"]
                    < cheapest_B["branches"]
                ):
                    cheapest_B = info_B

    # --------------------------------------------------------
    # Test all generated integer candidates against n.
    # --------------------------------------------------------

    exact_pairs = set()

    factor_candidates = (
        global_A
        | global_B
    )

    for x in factor_candidates:

        if x < 2:
            continue

        if n % x != 0:
            continue

        y = n // x

        if (
            x == p
            and y == q
        ) or (
            x == q
            and y == p
        ):
            exact_pairs.add(
                tuple(sorted((x, y)))
            )

    # --------------------------------------------------------
    # True-anchor direct diagnostic.
    # --------------------------------------------------------

    true_A_candidates = set()
    true_B_candidates = set()

    if true_state is not None:

        A = true_state["A"]
        B = true_state["B"]

        info_A = build_candidate_index(
            n,
            R1,
            A,
            sqrt_n,
        )

        info_B = build_candidate_index(
            n,
            R2,
            B,
            sqrt_n,
        )

        true_A_candidates = info_A.get(
            "candidates",
            set(),
        )

        true_B_candidates = info_B.get(
            "candidates",
            set(),
        )

        print()
        print("TRUE ANCHOR DIAGNOSTIC")

        print(
            f"  A widths="
            f"{[len(x) for x in A]}"
        )

        print(
            f"  B widths="
            f"{[len(x) for x in B]}"
        )

        print(
            f"  A CRT branches="
            f"{info_A.get('branches', 0):,}"
        )

        print(
            f"  B CRT branches="
            f"{info_B.get('branches', 0):,}"
        )

        print(
            f"  A candidates="
            f"{len(true_A_candidates)}"
        )

        print(
            f"  B candidates="
            f"{len(true_B_candidates)}"
        )

        print(
            f"  p in A candidates="
            f"{p in true_A_candidates}"
        )

        print(
            f"  q in B candidates="
            f"{q in true_B_candidates}"
        )

    else:

        print()
        print(
            "WARNING: true anchor was not "
            "among propagated states."
        )

    return {
        "states": len(states),
        "global_A": len(global_A),
        "global_B": len(global_B),
        "exact": exact_pairs,
        "true_found": (
            (p, q) in exact_pairs
            or (q, p) in exact_pairs
        ),
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 166")
    print("=" * 72)

    print(
        "Global CRT Candidate Intersection"
    )

    print()
    print(
        "Experiment 165 became computationally expensive because"
    )
    print(
        "it repeatedly performed CRT enumeration separately for"
    )
    print(
        "every anchor and both factor directions."
    )

    print()
    print(
        "Experiment 166 keeps the C-propagation stage, but treats"
    )
    print(
        "the CRT lift as a candidate-generation problem."
    )

    print()
    print(
        "For each surviving anchor:"
    )

    print(
        "    A-domains -> CRT candidates for x <= sqrt(n)"
    )

    print(
        "    B-domains -> CRT candidates for x <= sqrt(n)"
    )

    print()
    print(
        "All integer candidates are then merged globally and"
    )

    print(
        "tested against n."
    )

    print("=" * 72)

    summaries = []

    for bits in BIT_SIZES:

        print()
        print("-" * 72)
        print(
            f"BIT SIZE = {bits}"
        )
        print("-" * 72)

        for sample in range(
            1,
            SAMPLES_PER_SIZE + 1,
        ):

            p, q, n = make_semiprime(
                bits,
                rng,
            )

            t0 = time.perf_counter()

            result = analyze_sample(
                p,
                q,
                n,
            )

            elapsed = (
                time.perf_counter()
                - t0
            )

            print()
            print(
                "SAMPLE SUMMARY"
            )

            print(
                f"  global A candidates = "
                f"{result['global_A']:,}"
            )

            print(
                f"  global B candidates = "
                f"{result['global_B']:,}"
            )

            print(
                f"  exact factorizations = "
                f"{sorted(result['exact'])}"
            )

            print(
                f"  TRUE FOUND = "
                f"{result['true_found']}"
            )

            print(
                f"  total time = "
                f"{elapsed:.6f}s"
            )

            summaries.append(
                (
                    bits,
                    result["true_found"],
                    result["states"],
                    result["global_A"],
                    result["global_B"],
                    elapsed,
                )
            )

    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    for (
        bits,
        found,
        states,
        ca,
        cb,
        elapsed,
    ) in summaries:

        print(
            f"{bits:2d}-bit: "
            f"anchors={states} "
            f"A_candidates={ca:,} "
            f"B_candidates={cb:,} "
            f"found={found} "
            f"time={elapsed:.6f}s"
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 166")
    print("=" * 72)


if __name__ == "__main__":
    main()