#!/usr/bin/env python3
# =============================================================================
# EXPERIMENT 674
# =============================================================================
#
# DISCRIMINANT / 2-ADIC ROOT AMBIGUITY
#
# Established:
#
#   S = p + q
#   D = p - q
#
#   D^2 = S^2 - 4n
#
# Frame A:
#   C = q + 3
#   2C = S - D + 6
#   c = 9
#
# Frame B:
#   C = p + 1
#   2C = S + D + 2
#   c = 3
#
#   depth = v2(gcd(2C, n+c))
#
# The previous experiment showed:
#
#   n + S
#
# eventually separates the finite semiprime sample, but only
# at a relatively high bit depth.
#
# This experiment asks why.
#
# Since:
#
#   D^2 = S^2 - 4n,
#
# the missing information may be exactly the choice of a
# 2-adic square root of the discriminant.
#
# We measure:
#
#   1. discriminant consistency
#   2. number of square-root residues modulo 2^k
#   3. whether those root branches produce different C valuations
#   4. whether the ambiguity of (n,S) is explained entirely
#      by the root branches
#   5. whether adding the sign of D resolves the ambiguity
#
# No O(N^2) pairwise state comparison is used.
#
# =============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import gcd


PRIME_LIMIT = 6000
MAX_K = 16
INF = 10**9


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF
    x = abs(x)
    return (x & -x).bit_length() - 1


def sieve(limit: int) -> list[int]:
    if limit < 3:
        return []

    flags = bytearray(b"\x01") * (limit + 1)
    flags[0] = 0
    flags[1] = 0

    p = 2
    while p * p <= limit:
        if flags[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            flags[start : limit + 1 : p] = b"\x00" * count
        p += 1

    return [
        p
        for p in range(3, limit + 1, 2)
        if flags[p]
    ]


# =============================================================================
# STATE
# =============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str

    S: int
    D: int

    C: int
    c: int

    depth: int
    vc: int


def make_state(p: int, q: int) -> State:
    n = p * q
    S = p + q
    D = p - q

    if n % 4 == 3:
        # FRAME A
        frame = "A"
        C = q + 3
        c = 9
    else:
        # FRAME B
        frame = "B"
        C = p + 1
        c = 3

    depth = v2(gcd(2 * C, n + c))

    return State(
        n=n,
        p=p,
        q=q,
        frame=frame,
        S=S,
        D=D,
        C=C,
        c=c,
        depth=depth,
        vc=v2(C),
    )


def build_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(make_state(p, q))

    return states


# =============================================================================
# BASELINE
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = v2(gcd(2 * s.C, s.n + s.c))
        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1: DISCRIMINANT IDENTITY
# =============================================================================

def test_discriminant_identity(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: EXACT DISCRIMINANT IDENTITY")
    print("=" * 90)

    failures = 0

    for s in states:
        delta = s.S * s.S - 4 * s.n

        if delta != s.D * s.D:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# ROOT ENUMERATION
# =============================================================================

def square_roots_mod_2k(delta: int, k: int) -> list[int]:
    """
    Enumerate x with x^2 == delta (mod 2^k).

    This is intentionally used only for small k.
    At k <= 12 the search space is tiny and is useful
    as an exact diagnostic.

    For larger k we use the observed state residues
    rather than brute-forcing all residues.
    """
    mod = 1 << k
    target = delta % mod

    return [
        x
        for x in range(mod)
        if (x * x) % mod == target
    ]


# =============================================================================
# TEST 2
# =============================================================================

def test_root_branch_structure(states: list[State]) -> None:
    print("=" * 90)
    print("TEST 2: 2-ADIC DISCRIMINANT ROOT BRANCHES")
    print("=" * 90)

    for k in range(1, 13):
        mod = 1 << k

        signatures: dict[
            tuple,
            set[int]
        ] = defaultdict(set)

        root_counts: dict[
            tuple,
            set[int]
        ] = defaultdict(set)

        for s in states:
            delta = (s.S * s.S - 4 * s.n) % mod

            roots = square_roots_mod_2k(delta, k)

            sig = (
                s.frame,
                s.n % mod,
                s.S % mod,
            )

            signatures[sig].add(s.depth)
            root_counts[sig].add(len(roots))

        ambiguous_depth = sum(
            1
            for values in signatures.values()
            if len(values) > 1
        )

        ambiguous_roots = sum(
            1
            for values in root_counts.values()
            if len(values) > 1
        )

        root_hist: dict[int, int] = defaultdict(int)

        for values in root_counts.values():
            for count in values:
                root_hist[count] += 1

        print(
            f"k={k:2d} "
            f"signatures={len(signatures):7d} "
            f"depth_ambiguous={ambiguous_depth:6d} "
            f"root-count-ambiguous={ambiguous_roots:6d} "
            f"root-counts={sorted(root_hist)}"
        )

    print()


# =============================================================================
# TEST 3
# =============================================================================

def test_observed_root_branches(states: list[State]) -> int:
    """
    For each actual state compare the observed D residue with
    all square-root branches of the discriminant.

    The important question is whether D is just one of the
    available roots modulo 2^k.
    """

    print("=" * 90)
    print("TEST 3: ACTUAL D INSIDE THE DISCRIMINANT ROOT FIBER")
    print("=" * 90)

    failures = 0

    for k in range(1, 13):
        mod = 1 << k

        for s in states:
            delta = (s.S * s.S - 4 * s.n) % mod
            roots = square_roots_mod_2k(delta, k)

            if s.D % mod not in roots:
                failures += 1

    print(f"checks={len(states) * 12}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_root_branch_depths(states: list[State]) -> None:
    """
    For each (n,S) residue, enumerate all discriminant roots
    and ask what C valuations those roots could generate.

    This directly tests whether the ambiguity in n+S is the
    ambiguity of the 2-adic factor branch.
    """

    print("=" * 90)
    print("TEST 4: ROOT BRANCH -> C VALUATION")
    print("=" * 90)

    for k in range(2, 11):

        mod = 1 << k

        buckets: dict[
            tuple,
            set[int]
        ] = defaultdict(set)

        for s in states:
            sig = (
                s.frame,
                s.n % mod,
                s.S % mod,
            )

            delta = (s.S * s.S - 4 * s.n) % mod
            roots = square_roots_mod_2k(delta, k)

            # Every D root induces a candidate C.
            for d in roots:
                if s.frame == "A":
                    # 2C = S - D + 6
                    numerator = (
                        (s.S % mod)
                        - d
                        + 6
                    ) % mod
                else:
                    # 2C = S + D + 2
                    numerator = (
                        (s.S % mod)
                        + d
                        + 2
                    ) % mod

                # numerator is 2C.
                #
                # If numerator is zero modulo 2^k,
                # then v2(C) >= k-1.
                #
                # Otherwise:
                vc = v2(numerator) - 1

                if vc < 0:
                    # Odd numerator would contradict 2C.
                    continue

                buckets[sig].add(min(vc, k - 1))

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        print(
            f"k={k:2d} "
            f"signatures={len(buckets):7d} "
            f"candidate-vC ambiguous={ambiguous:6d}"
        )

    print()


# =============================================================================
# TEST 5
# =============================================================================

def test_sign_of_D(states: list[State]) -> None:
    """
    Compare:

        n + S

    against:

        n + S + sign(D)

    The hypothesis is that a large portion of the remaining
    ambiguity comes from root-branch selection.
    """

    print("=" * 90)
    print("TEST 5: DOES THE SIGN OF D RESOLVE n+S AMBIGUITY?")
    print("=" * 90)

    for k in range(2, MAX_K + 1):
        mod = 1 << k

        plain: dict[tuple, set[int]] = defaultdict(set)
        signed: dict[tuple, set[int]] = defaultdict(set)

        for s in states:

            plain_sig = (
                s.frame,
                s.n % mod,
                s.S % mod,
            )

            sign = 0
            if s.D > 0:
                sign = 1
            elif s.D < 0:
                sign = -1

            signed_sig = (
                s.frame,
                s.n % mod,
                s.S % mod,
                sign,
            )

            plain[plain_sig].add(s.depth)
            signed[signed_sig].add(s.depth)

        plain_ambiguous = sum(
            1
            for values in plain.values()
            if len(values) > 1
        )

        signed_ambiguous = sum(
            1
            for values in signed.values()
            if len(values) > 1
        )

        print(
            f"k={k:2d} "
            f"plain={plain_ambiguous:6d} "
            f"+sign(D)={signed_ambiguous:6d}"
        )

    print()


# =============================================================================
# TEST 6
# =============================================================================

def test_root_valuation_compression(states: list[State]) -> None:
    """
    Instead of storing D itself, store how close D is to zero
    2-adically.

    Compare:

        v2(D)

    and the actual branch-required quantity:

        v2(C).

    This is diagnostic only; previous experiments already
    showed neither simple valuation alone is sufficient.
    """

    print("=" * 90)
    print("TEST 6: ROOT VALUATION VS C VALUATION")
    print("=" * 90)

    buckets: dict[
        tuple,
        set[int]
    ] = defaultdict(set)

    for s in states:
        sig = (
            s.frame,
            v2(s.D),
            v2(s.S),
            v2(s.n + s.c),
        )

        buckets[sig].add(s.vc)

    ambiguous = {
        sig: values
        for sig, values in buckets.items()
        if len(values) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first collisions:")

        for i, (sig, values) in enumerate(
            ambiguous.items()
        ):
            if i >= 15:
                break

            print(
                f"    {sig} -> "
                f"{sorted(values)}"
            )

    print()


# =============================================================================
# TEST 7
# =============================================================================

def test_discriminant_residue_signature(states: list[State]) -> None:
    """
    Replace S by the discriminant:

        Delta = S^2 - 4n.

    Test whether:

        (frame, n, Delta)

    performs identically to:

        (frame, n, S)

    for depth reconstruction.
    """

    print("=" * 90)
    print("TEST 7: n + DISCRIMINANT RESIDUE")
    print("=" * 90)

    for k in range(2, MAX_K + 1):
        mod = 1 << k

        s_buckets: dict[
            tuple,
            set[int]
        ] = defaultdict(set)

        delta_buckets: dict[
            tuple,
            set[int]
        ] = defaultdict(set)

        for s in states:
            delta = (
                s.S * s.S
                - 4 * s.n
            )

            s_buckets[
                (
                    s.frame,
                    s.n % mod,
                    s.S % mod,
                )
            ].add(s.depth)

            delta_buckets[
                (
                    s.frame,
                    s.n % mod,
                    delta % mod,
                )
            ].add(s.depth)

        s_ambiguous = sum(
            1 for values in s_buckets.values()
            if len(values) > 1
        )

        delta_ambiguous = sum(
            1 for values in delta_buckets.values()
            if len(values) > 1
        )

        print(
            f"k={k:2d} "
            f"S-ambiguous={s_ambiguous:6d} "
            f"Delta-ambiguous={delta_ambiguous:6d}"
        )

    print()


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states: list[State]) -> None:
    print("=" * 90)
    print("TEST 8: EXAMPLES")
    print("=" * 90)

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
        485879,
        5579767,
    ]

    by_n = {s.n: s for s in states}

    for n in wanted:
        s = by_n.get(n)

        if s is None:
            continue

        delta = s.S * s.S - 4 * s.n

        print(
            f"n={s.n} p={s.p} q={s.q} "
            f"frame={s.frame}"
        )

        print(f"    S={s.S}")
        print(f"    D={s.D}")
        print(f"    D^2={s.D * s.D}")
        print(f"    Delta=S^2-4n={delta}")
        print(f"    v2(Delta)={v2(delta)}")
        print(f"    C={s.C}")
        print(f"    v2(C)={s.vc}")
        print(f"    v2(n+c)={v2(s.n + s.c)}")
        print(f"    depth={s.depth}")

        for k in (4, 8, 10, 12):
            mod = 1 << k

            roots = square_roots_mod_2k(
                delta,
                k,
            )

            print(
                f"    k={k:2d} "
                f"D mod 2^k={s.D % mod:6d} "
                f"roots={roots[:16]}"
            )

        print()

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 674 START")
    print("=" * 90)
    print()

    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve(PRIME_LIMIT)

    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)

    total_failures += test_discriminant_identity(states)

    test_root_branch_structure(states)

    total_failures += test_observed_root_branches(states)

    test_root_branch_depths(states)

    test_sign_of_D(states)

    test_root_valuation_compression(states)

    test_discriminant_residue_signature(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("The exact algebraic relation is:")
    print()
    print("    D^2 = S^2 - 4n")
    print()

    print("Therefore the previous n+S information boundary")
    print("can be interpreted as a 2-adic square-root problem:")
    print()
    print("    (n,S)")
    print("       -> Delta = S^2 - 4n")
    print("       -> D^2 = Delta")
    print("       -> D branch")
    print("       -> C")
    print("       -> depth")
    print()

    print("The key question is now:")
    print()
    print("    Is the remaining factor-side information")
    print("    exactly the choice of a 2-adic square root")
    print("    of the discriminant?")
    print()

    print("If adding the root branch removes the ambiguity,")
    print("then the missing information is not an arbitrary")
    print("factor invariant. It is the 2-adic root branch")
    print("of the quadratic x^2-Sx+n=0.")
    print()

    print(f"TOTAL TEST FAILURES={total_failures}")
    print()
    print("=" * 90)
    print("EXPERIMENT 674 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
