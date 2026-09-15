#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 582
# ==============================================================================
# HIERARCHICAL DOMAIN REDUCTION
#
# User-specified construction
# ---------------------------
#
#     k = 2^(z-1)
#
#     m = 2^z
#
# CASE A:
#
#     n mod m == m-1
#
#     p = k*y - k*x + 3
#     q = k*y + k*x - 3
#
# CASE B:
#
#     otherwise
#
#     p = (k*y + k*x - 3) / 3
#     q = k*y - k*x + 3
#
#
# This experiment asks:
#
#     How strongly does each modular level reduce the admissible
#     (x,y) domain?
#
# We measure:
#
#     - number of semiprimes
#     - number of factor-pair states
#     - number of valid (x,y) states
#     - unique x
#     - unique y
#     - unique (x,y)
#     - x/y ranges
#     - compression ratio
#     - branch populations
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 8

# Only odd semiprimes p*q.
ONLY_SEMIPRIMES = True

# Print individual valid states for small z.
SHOW_EXAMPLES = True
MAX_EXAMPLES = 10


# ==============================================================================
# DATA STRUCTURE
# ==============================================================================

@dataclass(frozen=True)
class State:

    n: int
    p: int
    q: int

    z: int
    k: int
    modulus: int

    branch: str

    x: int
    y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    if limit >= 0:
        sieve[0] = 0

    if limit >= 1:
        sieve[1] = 0

    root = isqrt(limit)

    for p in range(
        2,
        root + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        count = (
            (limit - start) // p
        ) + 1

        sieve[
            start:limit + 1:p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIME GENERATOR
# ==============================================================================

def generate_semiprimes(
    max_n,
    sieve,
):

    primes = [
        p
        for p in range(
            3,
            max_n + 1,
            2,
        )
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > max_n:
            break

        max_q = max_n // p

        for q in primes[i:]:

            if q > max_q:
                break

            result.append(
                (
                    p * q,
                    p,
                    q,
                )
            )

    return result


# ==============================================================================
# CASE A INVERSION
# ==============================================================================

def solve_case_a(
    p,
    q,
    k,
):
    """
    Case A:

        p = ky-kx+3
        q = ky+kx-3

    Therefore:

        p-3 = k(y-x)
        q+3 = k(y+x)

    Hence:

        y-x = (p-3)/k
        y+x = (q+3)/k
    """

    a = p - 3
    b = q + 3

    if a % k != 0:
        return None

    if b % k != 0:
        return None

    ymx = a // k
    ypx = b // k

    if (ypx + ymx) % 2 != 0:
        return None

    if (ypx - ymx) % 2 != 0:
        return None

    x = (
        ypx - ymx
    ) // 2

    y = (
        ypx + ymx
    ) // 2

    return x, y


# ==============================================================================
# CASE B INVERSION
# ==============================================================================

def solve_case_b(
    p,
    q,
    k,
):
    """
    Case B:

        p = (ky+kx-3)/3
        q = ky-kx+3

    Therefore:

        3p+3 = k(y+x)
        q-3  = k(y-x)

    Hence:

        y+x = (3p+3)/k
        y-x = (q-3)/k
    """

    a = 3 * p + 3
    b = q - 3

    if a % k != 0:
        return None

    if b % k != 0:
        return None

    ypx = a // k
    ymx = b // k

    if (ypx + ymx) % 2 != 0:
        return None

    if (ypx - ymx) % 2 != 0:
        return None

    x = (
        ypx - ymx
    ) // 2

    y = (
        ypx + ymx
    ) // 2

    return x, y


# ==============================================================================
# FORWARD VALIDATION
# ==============================================================================

def forward_case_a(
    x,
    y,
    k,
):
    p = k * y - k * x + 3
    q = k * y + k * x - 3

    return p, q


def forward_case_b(
    x,
    y,
    k,
):
    numerator = (
        k * y
        + k * x
        - 3
    )

    if numerator % 3 != 0:
        return None

    p = numerator // 3

    q = (
        k * y
        - k * x
        + 3
    )

    return p, q


# ==============================================================================
# PROCESS ONE FACTOR PAIR
# ==============================================================================

def process_pair(
    n,
    p,
    q,
    z,
):
    """
    Return the State if the pair satisfies the user's branch construction.
    """

    k = 1 << (z - 1)
    modulus = 1 << z

    residue = n % modulus

    if residue == modulus - 1:

        solution = solve_case_a(
            p,
            q,
            k,
        )

        if solution is None:
            return None

        x, y = solution

        check_p, check_q = forward_case_a(
            x,
            y,
            k,
        )

        if (
            check_p != p
            or check_q != q
        ):
            raise RuntimeError(
                "Case A reconstruction failed"
            )

        branch = "A"

    else:

        solution = solve_case_b(
            p,
            q,
            k,
        )

        if solution is None:
            return None

        x, y = solution

        result = forward_case_b(
            x,
            y,
            k,
        )

        if result is None:
            return None

        check_p, check_q = result

        if (
            check_p != p
            or check_q != q
        ):
            raise RuntimeError(
                "Case B reconstruction failed"
            )

        branch = "B"

    if x < 0 or y < 0:
        return None

    if p * q != n:
        raise RuntimeError(
            "Factorization mismatch"
        )

    return State(
        n=n,
        p=p,
        q=q,
        z=z,
        k=k,
        modulus=modulus,
        branch=branch,
        x=x,
        y=y,
    )


# ==============================================================================
# STATISTICS
# ==============================================================================

def state_statistics(
    states,
):

    if not states:

        return {
            "count": 0,
        }

    xs = [
        s.x
        for s in states
    ]

    ys = [
        s.y
        for s in states
    ]

    xy = {
        (
            s.x,
            s.y,
        )
        for s in states
    }

    ns = {
        s.n
        for s in states
    }

    pq = {
        (
            s.p,
            s.q,
        )
        for s in states
    }

    return {
        "count": len(states),

        "unique_n": len(ns),

        "unique_pq": len(pq),

        "unique_x": len(
            set(xs)
        ),

        "unique_y": len(
            set(ys)
        ),

        "unique_xy": len(xy),

        "xmin": min(xs),
        "xmax": max(xs),

        "ymin": min(ys),
        "ymax": max(ys),
    }


# ==============================================================================
# BRANCH STATISTICS
# ==============================================================================

def branch_statistics(
    states,
):

    result = {}

    for branch in (
        "A",
        "B",
    ):

        subset = [
            s
            for s in states
            if s.branch == branch
        ]

        result[branch] = (
            state_statistics(
                subset
            )
        )

    return result


# ==============================================================================
# MAIN EXPERIMENT
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 582 START"
    )
    print("=" * 90)

    print()
    print(
        "HIERARCHICAL DOMAIN REDUCTION"
    )

    print()
    print(
        "Testing:"
    )

    print(
        "    k = 2^(z-1)"
    )

    print(
        "    Case A: residue == 2^z - 1"
    )

    print(
        "    Case B: all other residues"
    )

    # --------------------------------------------------------------------------
    # Generate semiprimes.
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    print(
        "[2] SEMIPRIME GENERATION"
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    total semiprimes = "
        f"{len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Baseline.
    # --------------------------------------------------------------------------

    baseline_n = {
        n
        for n, _, _
        in semiprimes
    }

    baseline_pq = {
        (
            p,
            q,
        )
        for _, p, q
        in semiprimes
    }

    print()
    print(
        "[3] BASELINE"
    )

    print(
        f"    unique n  = {len(baseline_n)}"
    )

    print(
        f"    factor states = {len(baseline_pq)}"
    )

    # --------------------------------------------------------------------------
    # Per-level experiment.
    # --------------------------------------------------------------------------

    all_levels = {}

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        k = 1 << (z - 1)
        modulus = 1 << z

        print()
        print("=" * 90)
        print(
            f"LEVEL z={z}"
        )
        print(
            f"    modulus = {modulus}"
        )
        print(
            f"    k       = {k}"
        )
        print("=" * 90)

        states = []

        for n, p, q in semiprimes:

            state = process_pair(
                n,
                p,
                q,
                z,
            )

            if state is not None:
                states.append(
                    state
                )

        all_levels[z] = states

        stats = state_statistics(
            states
        )

        branches = branch_statistics(
            states
        )

        # ----------------------------------------------------------------------
        # Overall.
        # ----------------------------------------------------------------------

        print()
        print(
            "OVERALL"
        )

        print(
            f"    valid states = "
            f"{stats.get('count', 0)}"
        )

        print(
            f"    unique n    = "
            f"{stats.get('unique_n', 0)}"
        )

        print(
            f"    unique p,q  = "
            f"{stats.get('unique_pq', 0)}"
        )

        print(
            f"    unique x    = "
            f"{stats.get('unique_x', 0)}"
        )

        print(
            f"    unique y    = "
            f"{stats.get('unique_y', 0)}"
        )

        print(
            f"    unique x,y  = "
            f"{stats.get('unique_xy', 0)}"
        )

        if stats.get("count", 0):

            print(
                f"    x range     = "
                f"[{stats['xmin']}, "
                f"{stats['xmax']}]"
            )

            print(
                f"    y range     = "
                f"[{stats['ymin']}, "
                f"{stats['ymax']}]"
            )

        # ----------------------------------------------------------------------
        # Compression relative to baseline.
        # ----------------------------------------------------------------------

        valid_count = stats.get(
            "count",
            0,
        )

        print()
        print(
            "REDUCTION"
        )

        if len(semiprimes):

            ratio = (
                valid_count
                / len(semiprimes)
            )

            reduction = (
                1.0 - ratio
            )

            print(
                f"    surviving fraction = "
                f"{ratio:.8f}"
            )

            print(
                f"    reduction          = "
                f"{reduction:.8f}"
            )

        # ----------------------------------------------------------------------
        # Branch counts.
        # ----------------------------------------------------------------------

        print()
        print(
            "BRANCHES"
        )

        for branch in (
            "A",
            "B",
        ):

            branch_states = [
                s
                for s in states
                if s.branch == branch
            ]

            branch_stats = branches[
                branch
            ]

            print()
            print(
                f"    CASE {branch}"
            )

            print(
                f"        states = "
                f"{len(branch_states)}"
            )

            if branch_states:

                print(
                    f"        unique x = "
                    f"{branch_stats['unique_x']}"
                )

                print(
                    f"        unique y = "
                    f"{branch_stats['unique_y']}"
                )

                print(
                    f"        unique xy = "
                    f"{branch_stats['unique_xy']}"
                )

                print(
                    f"        x range = "
                    f"[{branch_stats['xmin']}, "
                    f"{branch_stats['xmax']}]"
                )

                print(
                    f"        y range = "
                    f"[{branch_stats['ymin']}, "
                    f"{branch_stats['ymax']}]"
                )

        # ----------------------------------------------------------------------
        # Sample states.
        # ----------------------------------------------------------------------

        if SHOW_EXAMPLES:

            print()
            print(
                "EXAMPLES"
            )

            for state in states[
                :MAX_EXAMPLES
            ]:

                print(
                    f"    {state.n} "
                    f"{state.p} "
                    f"{state.q} "
                    f"{state.x} "
                    f"{state.y} "
                    f"branch={state.branch}"
                )

    # ==========================================================================
    # CROSS-LEVEL COMPARISON
    # ==========================================================================

    print()
    print("=" * 90)
    print(
        "CROSS-LEVEL DOMAIN SHRINKAGE"
    )
    print("=" * 90)

    previous = None

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        states = all_levels[z]

        current_xy = {
            (
                s.x,
                s.y,
            )
            for s in states
        }

        if previous is None:

            print(
                f"z={z}: "
                f"unique(x,y)={len(current_xy)}"
            )

        else:

            previous_xy = previous

            intersection = (
                current_xy
                & previous_xy
            )

            print(
                f"z={z}: "
                f"unique(x,y)={len(current_xy):<8} "
                f"previous={len(previous_xy):<8} "
                f"intersection={len(intersection):<8}"
            )

        previous = current_xy

    # ==========================================================================
    # MODULAR RESIDUE DISTRIBUTION
    # ==========================================================================

    print()
    print("=" * 90)
    print(
        "RESIDUE DISTRIBUTION"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        modulus = 1 << z

        counts = defaultdict(int)

        for n, _, _ in semiprimes:

            counts[
                n % modulus
            ] += 1

        print()
        print(
            f"MOD {modulus}"
        )

        for residue in sorted(
            counts
        ):

            print(
                f"    {residue:<4}: "
                f"{counts[residue]}"
            )

    # ==========================================================================
    # FINAL
    # ==========================================================================

    print()
    print("=" * 90)
    print(
        "FINAL AUDIT"
    )
    print("=" * 90)

    print(
        """
For each z we measured the surviving domain after applying the
user-defined branch equations.

The quantities of greatest interest are:

    valid states
    unique x
    unique y
    unique (x,y)
    x range
    y range

If these shrink systematically as z increases, then the modular
condition is acting as a domain restriction on the x,y function
inputs.

The next step after this experiment is to recover the actual
function on each surviving child domain and compare it to the
parent function after the corresponding coordinate transformation.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 582 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
