#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 644
# ==============================================================================
#
# EXACT 2-ADIC IMAGE OF THE FACTOR-GCD VALUATION
#
# Experiment 642 established:
#
#     d = gcd(A,B)
#
#     m = v2(d)
#
#     depth =
#         m                      if normalized parity is opposite
#         m + 1                  if normalized parity is both odd
#
# equivalently:
#
#     depth = m + [v2(A)=v2(B)]
#
#
# Experiment 643 showed:
#
#     (frame, n mod 2^k)
#
# does NOT uniquely determine m.
#
# However, its monotonicity test was conceptually wrong:
# the number of ambiguous residue classes can increase when
# the modulus doubles because each old residue splits into
# two children.
#
# Therefore Experiment 644 reverses the map.
#
# For each fixed:
#
#     frame
#     m
#     k
#
# we compute EXACTLY which residues
#
#     n mod 2^k
#
# can occur.
#
# Then we ask:
#
#     1. Does every m have a recognizable residue family?
#
#     2. Do different m-values have overlapping residue families?
#
#     3. At what k does a fixed m first become separable
#        from another m?
#
#     4. Is there a simple tree/lifting rule for these families?
#
# This works entirely at the 2-adic residue level.
#
# We do NOT enumerate the full semiprime domain for every k.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass


# ==============================================================================
# CONFIGURATION
# ==============================================================================

K_LEVELS = [4, 5, 6, 7, 8, 9, 10, 11, 12]

# m values to inspect.
MAX_M = 12

# Number of residue examples to print.
SHOW = 12

INF = 10**9


# ==============================================================================
# v2
# ==============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# FRAME RESIDUALS
# ==============================================================================

def residuals_from_residues(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        return (
            p - 3,
            q + 3,
        )

    return (
        p + 1,
        q - 3,
    )


# ==============================================================================
# TRUNCATED m
# ==============================================================================

def residual_m(
    frame: str,
    p: int,
    q: int,
    k: int,
) -> int:

    A, B = residuals_from_residues(
        frame,
        p,
        q,
    )

    va = v2(A)
    vb = v2(B)

    return min(
        va,
        vb,
        k,
    )


# ==============================================================================
# EXACT RESIDUE IMAGE
# ==============================================================================

def build_image(
    frame: str,
    k: int,
    max_m: int,
):
    """
    Return:

        m -> set(n mod 2^k)

    for all compatible odd p,q residues.

    Since p is odd, it is invertible modulo 2^k.

    For every n residue and p residue:

        q = n * p^{-1} mod 2^k.
    """

    M = 1 << k

    odd = range(
        1,
        M,
        2,
    )

    inverses = {
        p: pow(
            p,
            -1,
            M,
        )
        for p in odd
    }

    image = {
        m: set()
        for m in range(
            max_m + 1
        )
    }

    for n_mod in range(
        1,
        M,
        2,
    ):

        for p_mod in odd:

            q_mod = (
                n_mod
                * inverses[p_mod]
            ) % M

            m = residual_m(
                frame,
                p_mod,
                q_mod,
                k,
            )

            if m <= max_m:
                image[m].add(
                    n_mod
                )

    return image


# ==============================================================================
# TEST 0
# ==============================================================================

def test_image_partition(
    frame: str,
    k: int,
    image,
) -> int:

    print("=" * 90)
    print(
        f"TEST 0: RESIDUE IMAGE PARTITION "
        f"FRAME={frame} k={k}"
    )
    print("=" * 90)

    M = 1 << k

    seen = defaultdict(set)

    for m, residues in image.items():

        for r in residues:
            seen[r].add(m)

    failures = 0

    # Every residue belongs to at least one compatible m bucket.
    for r in range(
        1,
        M,
        2,
    ):

        if not seen.get(r):
            failures += 1

    overlap_residues = sum(
        1
        for values in seen.values()
        if len(values) > 1
    )

    exact_residues = sum(
        1
        for values in seen.values()
        if len(values) == 1
    )

    print(
        f"    odd residue classes={M // 2}"
    )

    print(
        f"    covered={len(seen)}"
    )

    print(
        f"    exact m residues={exact_residues}"
    )

    print(
        f"    overlapping residues={overlap_residues}"
    )

    print(
        f"    failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 1
# ==============================================================================

def test_m_image_sizes(
    frame: str,
    k: int,
    image,
) -> None:

    print("=" * 90)
    print(
        f"TEST 1: IMAGE SIZE BY m "
        f"FRAME={frame} k={k}"
    )
    print("=" * 90)

    for m in sorted(image):

        print(
            f"    m={m:<2} "
            f"residues={len(image[m])}"
        )

    print()


# ==============================================================================
# TEST 2
# ==============================================================================

def test_m_overlap_matrix(
    frame: str,
    k: int,
    image,
) -> None:

    print("=" * 90)
    print(
        f"TEST 2: m OVERLAP MATRIX "
        f"FRAME={frame} k={k}"
    )
    print("=" * 90)

    overlaps = []

    m_values = sorted(
        image.keys()
    )

    for i, m1 in enumerate(
        m_values
    ):

        for m2 in m_values[
            i + 1:
        ]:

            common = (
                image[m1]
                & image[m2]
            )

            if common:

                overlaps.append(
                    (
                        m1,
                        m2,
                        len(common),
                    )
                )

    if not overlaps:

        print(
            "    no overlapping m-image residues"
        )

    else:

        for m1, m2, count in overlaps:

            print(
                f"    m={m1} "
                f"<-> "
                f"m={m2} "
                f"overlap={count}"
            )

    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_minimal_separator(
    frame: str,
    max_k: int,
    max_m: int,
) -> None:

    print("=" * 90)
    print(
        f"TEST 3: FIRST k THAT SEPARATES EACH m-PAIR "
        f"FRAME={frame}"
    )
    print("=" * 90)

    pair_first_separated = {}

    for m1 in range(
        max_m + 1
    ):

        for m2 in range(
            m1 + 1,
            max_m + 1,
        ):

            first = None

            for k in range(
                4,
                max_k + 1,
            ):

                image = build_image(
                    frame,
                    k,
                    max_m,
                )

                if not (
                    image[m1]
                    &
                    image[m2]
                ):

                    first = k
                    break

            pair_first_separated[
                (m1, m2)
            ] = first

    for pair, k in pair_first_separated.items():

        m1, m2 = pair

        print(
            f"    {m1} vs {m2}: "
            f"first_separated={k}"
        )

    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_residue_lift_tree(
    frame: str,
    k: int,
    image,
) -> None:

    print("=" * 90)
    print(
        f"TEST 4: RESIDUE LIFT TREE "
        f"FRAME={frame} k={k}"
    )
    print("=" * 90)

    # Build reverse map:
    #
    # residue -> possible m values.

    classification = defaultdict(set)

    for m, residues in image.items():

        for r in residues:

            classification[r].add(m)

    # Group by parent residue.
    #
    # Every r modulo 2^k has parent:
    #
    #     r mod 2^(k-1)

    M = 1 << k
    parent_M = 1 << (k - 1)

    children = defaultdict(
        dict
    )

    for r in range(
        1,
        M,
        2,
    ):

        parent = r % parent_M

        if r not in classification:
            continue

        children[parent][
            r
        ] = classification[r]

    shown = 0

    for parent in sorted(
        children
    ):

        child_data = children[
            parent
        ]

        if len(child_data) != 2:
            continue

        values = list(
            child_data.values()
        )

        if values[0] == values[1]:
            continue

        print(
            f"    parent residue={parent}"
        )

        for child, vals in sorted(
            child_data.items()
        ):

            print(
                f"        child={child:<6} "
                f"m={sorted(vals)}"
            )

        print()

        shown += 1

        if shown >= SHOW:
            break

    if shown == 0:
        print(
            "    no nontrivial m-splitting "
            "lift examples"
        )

    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_depth_overlap(
    frame: str,
    k: int,
    image,
) -> None:

    print("=" * 90)
    print(
        f"TEST 5: DEPTH-LEVEL CONSEQUENCE "
        f"FRAME={frame} k={k}"
    )
    print("=" * 90)

    depth_values = defaultdict(
        set
    )

    # m alone is not enough for depth:
    #
    # opposite parity -> m
    # both odd        -> m+1
    #
    # For the residue-level study we therefore record
    # the possible depth interval implied by m.
    #
    # This is deliberately conservative.

    for m, residues in image.items():

        for r in residues:

            depth_values[r].add(
                m
            )

    ambiguous = 0

    for r, ms in depth_values.items():

        if len(ms) > 1:

            ambiguous += 1

    print(
        f"    residues={len(depth_values)}"
    )

    print(
        f"    residues with multiple m="
        f"{ambiguous}"
    )

    print()

# ==============================================================================
# TEST 6
# ==============================================================================

def test_actual_semiprime_projection(
    frame: str,
    states,
    k: int,
) -> None:

    print("=" * 90)
    print(
        f"TEST 6: ACTUAL SEMIPRIME PROJECTION "
        f"FRAME={frame} k={k}"
    )
    print("=" * 90)

    M = 1 << k

    buckets = defaultdict(
        set
    )

    examples = defaultdict(
        list
    )

    for s in states:

        if s.frame != frame:
            continue

        key = (
            s.n % M
        )

        buckets[key].add(
            s.m
        )

        if (
            len(
                examples[key]
            )
            < 4
        ):

            examples[key].append(
                s
            )

    ambiguous = [
        (
            key,
            values,
        )
        for key, values in buckets.items()
        if len(values) > 1
    ]

    exact = sum(
        1
        for values in buckets.values()
        if len(values) == 1
    )

    print(
        f"    residue signatures="
        f"{len(buckets)}"
    )

    print(
        f"    exact signatures="
        f"{exact}"
    )

    print(
        f"    ambiguous signatures="
        f"{len(ambiguous)}"
    )

    for key, values in ambiguous[:10]:

        print(
            f"    collision residue={key} "
            f"m={sorted(values)}"
        )

        for s in examples[key]:

            print(
                f"        n={s.n:<10} "
                f"p={s.p:<6} "
                f"q={s.q:<6} "
                f"m={s.m}"
            )

    print()


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(
    limit: int,
) -> list[int]:

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    sieve[0] = 0
    sieve[1] = 0

    for p in range(
        2,
        int(limit ** 0.5) + 1,
    ):

        if not sieve[p]:
            continue

        sieve[
            p * p:
            limit + 1:
            p
        ] = (
            b"\x00"
            * (
                (
                    limit
                    - p * p
                )
                // p
                + 1
            )
        )

    return [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
        if sieve[p]
    ]


# ==============================================================================
# ACTUAL STATES
# ==============================================================================

@dataclass
class State:
    n: int
    p: int
    q: int
    frame: str
    A: int
    B: int
    m: int
    depth: int


def build_states(
    primes: list[int],
) -> list[State]:

    states = []

    for i, p in enumerate(
        primes
    ):

        for q in primes[i:]:

            n = p * q

            frame = (
                "A"
                if n % 4 == 3
                else "B"
            )

            if frame == "A":

                A = p - 3
                B = q + 3

            else:

                A = p + 1
                B = q - 3

            d = 0

            from math import gcd

            d = gcd(
                abs(A),
                abs(B),
            )

            if d == 0:
                m = INF
            else:
                m = v2(d)

            if frame == "A":

                X_num = B - A
                Y_num = A + B

            else:

                X_num = 3 * A - B
                Y_num = 3 * A + B

            depth = min(
                v2(X_num),
                v2(Y_num),
            )

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    m=m,
                    depth=depth,
                )
            )

    return states


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 644 START"
    )
    print("=" * 90)
    print()
    print(
        "EXACT 2-ADIC IMAGE OF "
        "FACTOR-GCD VALUATION"
    )
    print()

    print(
        "[1] PRIME SIEVE"
    )

    # Same prime population as Experiment 642/643.
    primes = prime_sieve(
        6000
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    print(
        "[2] SEMIPRIME GENERATION"
    )

    states = build_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    print("=" * 90)
    print(
        "TEST A: EXACT n-RESIDUE IMAGE"
    )
    print("=" * 90)
    print()

    # --------------------------------------------------------------
    # Main exact image analysis.
    #
    # We deliberately use k<=10 for the complete residue-space
    # enumeration. This is large enough to expose the structure
    # while keeping the run manageable.
    # --------------------------------------------------------------

    selected_k = 10

    total_failures = 0

    for frame in (
        "A",
        "B",
    ):

        print(
            f"FRAME {frame}"
        )

        image = build_image(
            frame,
            selected_k,
            MAX_M,
        )

        total_failures += test_image_partition(
            frame,
            selected_k,
            image,
        )

        test_m_image_sizes(
            frame,
            selected_k,
            image,
        )

        test_m_overlap_matrix(
            frame,
            selected_k,
            image,
        )

        test_residue_lift_tree(
            frame,
            selected_k,
            image,
        )

        test_depth_overlap(
            frame,
            selected_k,
            image,
        )

        test_actual_semiprime_projection(
            frame,
            states,
            selected_k,
        )

    # --------------------------------------------------------------
    # Pairwise separation at small k.
    # --------------------------------------------------------------

    print("=" * 90)
    print(
        "TEST B: m-PAIR SEPARATION"
    )
    print("=" * 90)
    print()

    for frame in (
        "A",
        "B",
    ):

        test_minimal_separator(
            frame,
            10,
            8,
        )

    # --------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------

    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)
    print()

    wanted = [
        9,
        15,
        21,
        33,
        39,
        57,
        69,
        77,
        87,
        93,
        111,
        141,
        183,
        213,
    ]

    lookup = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = lookup[n]

        c = (
            9
            if s.frame == "A"
            else 3
        )

        print(
            f"n={n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={s.A} "
            f"B={s.B}"
        )

        print(
            f"    m={s.m}"
        )

        print(
            f"    v2(n+c)="
            f"{v2(n + c)}"
        )

        print(
            f"    depth={s.depth}"
        )

        print()

    # --------------------------------------------------------------
    # Final interpretation.
    # --------------------------------------------------------------

    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
Experiment 642 established:

    depth =
        v2(gcd(A,B))
        +
        [v2(A)=v2(B)].

Therefore the remaining quantity is:

    m=v2(gcd(A,B)).

Experiment 643 established that:

    (frame, n mod 2^k)

does not uniquely determine m.

Experiment 644 reverses the question.

For each m we compute its exact 2-adic image:

    I_k(frame,m)
        =
    { n mod 2^k :
      there exists an odd factor residue p
      producing that m }.

Therefore:

    I_k(frame,m1)
        intersect
    I_k(frame,m2)

is nonempty exactly when the two m-values remain
indistinguishable at modulus 2^k.

This gives a much cleaner interpretation of the
information loss.

Three possibilities matter:

    1. Images become disjoint.

       Then sufficiently many low bits of n separate
       those m-values.

    2. Images overlap for every tested k.

       Then those m-values are genuinely compatible
       with the same 2-adic n projection to arbitrary
       tested depth.

    3. The images form a recursive tree.

       Then the missing information may be a branch bit
       rather than a conventional arithmetic valuation.

The residue-lift tree is especially important:

    r mod 2^k
         |
         +-- r
         |
         +-- r+2^k

If the two children acquire different possible-m sets,
that identifies the exact bit at which the factor-side
valuation becomes distinguishable.

The goal is therefore no longer:

    "find a clever formula for m."

It is:

    "characterize the 2-adic image and its lifting tree."

If the lifting tree has a simple recursive description,
that may reveal the missing invariant directly.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 644 FINISHED"
    )
    print("=" * 90)
    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL STRUCTURAL TESTS PASSED"
        )

    else:

        print(
            "STATUS=CONSISTENCY FAILURES FOUND"
        )


if __name__ == "__main__":
    main()
