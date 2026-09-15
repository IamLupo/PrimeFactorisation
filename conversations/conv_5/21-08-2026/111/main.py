#!/usr/bin/env python3
# =============================================================================
# EXPERIMENT 673
# =============================================================================
#
# ONE SYMMETRIC COORDINATE + n : INFORMATION BOUNDARY
#
# Established:
#
#   FRAME A:
#       C = q + 3
#       c = 9
#
#   FRAME B:
#       C = p + 1
#       c = 3
#
#   depth = v2(gcd(2C, n+c))
#
# Equivalent:
#
#   depth = min(v2(C)+1, v2(n+c))
#
# Also:
#
#   S = p+q
#   D = p-q
#
# and:
#
#   FRAME A:
#       2C = S - D + 6
#
#   FRAME B:
#       2C = S + D + 2
#
# This experiment asks:
#
#   Given n and ONLY ONE symmetric coordinate,
#
#       (frame, n mod 2^k, S mod 2^k)
#
#   or
#
#       (frame, n mod 2^k, D mod 2^k)
#
#   when does the signature determine the exact depth?
#
# We also test:
#
#   (frame, n mod 2^k, S mod 2^(k+1))
#   (frame, n mod 2^k, D mod 2^(k+1))
#
# because 2C is reconstructed from S,D with a one-bit
# division by 2.
#
# No O(N^2) pairwise comparison is performed.
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
# NUMBER THEORY
# =============================================================================

def v2(x: int) -> int:
    """2-adic valuation, with v2(0)=INF."""
    if x == 0:
        return INF
    x = abs(x)
    return (x & -x).bit_length() - 1


def sieve(limit: int) -> list[int]:
    if limit < 3:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = 0
    is_prime[1] = 0

    p = 2
    while p * p <= limit:
        if is_prime[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            is_prime[start : limit + 1 : p] = b"\x00" * count
        p += 1

    return [
        p
        for p in range(3, limit + 1, 2)
        if is_prime[p]
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

    A: int
    B: int
    C: int

    S: int
    D: int

    c: int
    depth: int

    vc: int
    wn: int


# =============================================================================
# STATE CONSTRUCTION
# =============================================================================

def make_state(p: int, q: int) -> State:
    n = p * q

    S = p + q
    D = p - q

    if n % 4 == 3:
        # FRAME A
        #
        # C = q+3
        # c = 9
        A = p - 3
        B = q + 3
        C = B
        c = 9

    else:
        # FRAME B
        #
        # C = p+1
        # c = 3
        A = p + 1
        B = q - 3
        C = A
        c = 3

    depth = v2(gcd(2 * C, n + c))

    return State(
        n=n,
        p=p,
        q=q,
        frame="A" if n % 4 == 3 else "B",
        A=A,
        B=B,
        C=C,
        S=S,
        D=D,
        c=c,
        depth=depth,
        vc=v2(C),
        wn=v2(n + c),
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
# GENERIC BUCKET TEST
# =============================================================================

def bucket_test(
    states: list[State],
    name: str,
    signature_builder,
    value_builder,
    show_examples: bool = False,
) -> tuple[int, int]:
    buckets: dict[tuple, set[int]] = defaultdict(set)
    examples: dict[tuple, State] = {}

    for s in states:
        sig = signature_builder(s)

        buckets[sig].add(value_builder(s))
        examples.setdefault(sig, s)

    ambiguous = {
        sig: values
        for sig, values in buckets.items()
        if len(values) > 1
    }

    print(name)
    print(
        f"    signatures={len(buckets):7d} "
        f"ambiguous={len(ambiguous):6d}"
    )

    if show_examples and ambiguous:
        print("    first collisions:")

        for i, (sig, values) in enumerate(ambiguous.items()):
            if i >= 10:
                break

            s = examples[sig]

            print(
                f"        signature={sig} "
                f"depths={sorted(values)}"
            )

            print(
                f"            n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame}"
            )

    return len(buckets), len(ambiguous)


# =============================================================================
# TEST 1
# =============================================================================

def test_single_coordinate(
    states: list[State],
    coordinate: str,
) -> None:

    print("=" * 90)
    print(
        "TEST 1: n + "
        + coordinate
        + " RESIDUE -> EXACT DEPTH"
    )
    print("=" * 90)

    first_exact = None

    for k in range(1, MAX_K + 1):
        mod = 1 << k

        if coordinate == "S":
            build = lambda s, mod=mod: (
                s.frame,
                s.n % mod,
                s.S % mod,
            )

        else:
            build = lambda s, mod=mod: (
                s.frame,
                s.n % mod,
                s.D % mod,
            )

        _, ambiguous = bucket_test(
            states,
            f"    k={k:2d}",
            build,
            lambda s: s.depth,
            show_examples=(k in (2, 4, 8)),
        )

        if ambiguous == 0 and first_exact is None:
            first_exact = k

    print()

    if first_exact is None:
        print("first exact k=None")
    else:
        print(f"first exact k={first_exact}")

    print()


# =============================================================================
# TEST 2
# =============================================================================

def test_shifted_coordinate(
    states: list[State],
    coordinate: str,
) -> None:

    print("=" * 90)
    print(
        "TEST 2: n + "
        + coordinate
        + " WITH ONE EXTRA BIT"
    )
    print("=" * 90)

    first_exact = None

    for k in range(1, MAX_K + 1):
        nmod = 1 << k
        coordmod = 1 << (k + 1)

        if coordinate == "S":
            build = lambda s, nmod=nmod, coordmod=coordmod: (
                s.frame,
                s.n % nmod,
                s.S % coordmod,
            )

        else:
            build = lambda s, nmod=nmod, coordmod=coordmod: (
                s.frame,
                s.n % nmod,
                s.D % coordmod,
            )

        _, ambiguous = bucket_test(
            states,
            f"    k={k:2d}",
            build,
            lambda s: s.depth,
        )

        if ambiguous == 0 and first_exact is None:
            first_exact = k

    print()

    if first_exact is None:
        print("first exact k=None")
    else:
        print(f"first exact k={first_exact}")

    print()


# =============================================================================
# TEST 3
# =============================================================================

def test_symmetric_pair_reconstruction(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: RECONSTRUCT 2C FROM ONE SYMMETRIC COORDINATE")
    print("=" * 90)

    failures = 0

    for s in states:
        if s.frame == "A":
            reconstructed = s.S - s.D + 6
        else:
            reconstructed = s.S + s.D + 2

        if reconstructed != 2 * s.C:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_valuation_from_one_coordinate(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: v2(C) FROM S/D + n MODULAR DATA")
    print("=" * 90)

    failures = 0

    # The theoretical relation is:
    #
    # FRAME A:
    #   2C = S-D+6
    #
    # FRAME B:
    #   2C = S+D+2
    #
    # We verify that sufficiently large residues reproduce
    # the actual C valuation.

    for s in states:
        k = MAX_K

        mod = 1 << (k + 1)

        S = s.S % mod
        D = s.D % mod

        if s.frame == "A":
            two_c = (S - D + 6) % mod
        else:
            two_c = (S + D + 2) % mod

        # Since two_c is only known modulo 2^(k+1),
        # v2(two_c) is exact whenever it is <= k.
        recovered_vc = v2(two_c) - 1

        if recovered_vc != s.vc:
            # For the rare case where C's valuation is at
            # the modulus boundary, inspect whether the
            # observed residue is actually zero.
            if not (
                s.vc >= k
                and two_c == 0
            ):
                failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_minimal_information(states: list[State]) -> None:
    print("=" * 90)
    print("TEST 5: MINIMAL INFORMATION BOUNDARY")
    print("=" * 90)

    print()
    print("S-only:")
    test_single_coordinate(states, "S")

    print("D-only:")
    test_single_coordinate(states, "D")

    print("S with one extra bit:")
    test_shifted_coordinate(states, "S")

    print("D with one extra bit:")
    test_shifted_coordinate(states, "D")


# =============================================================================
# TEST 6
# =============================================================================

def test_depth_formula_from_coordinate(
    states: list[State],
    coordinate: str,
) -> None:

    print("=" * 90)
    print(
        "TEST 6: "
        + coordinate
        + " + n -> C VALUATION -> DEPTH"
    )
    print("=" * 90)

    first_exact = None

    for k in range(1, MAX_K + 1):
        nmod = 1 << k
        coordmod = 1 << (k + 1)

        buckets: dict[tuple, set[int]] = defaultdict(set)

        for s in states:
            if coordinate == "S":
                coord = s.S % coordmod
            else:
                coord = s.D % coordmod

            sig = (
                s.frame,
                s.n % nmod,
                coord,
            )

            # Store depth directly.
            buckets[sig].add(s.depth)

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        print(
            f"    k={k:2d} "
            f"signatures={len(buckets):7d} "
            f"ambiguous={ambiguous:6d}"
        )

        if ambiguous == 0 and first_exact is None:
            first_exact = k

    print()

    if first_exact is None:
        print("first exact depth k=None")
    else:
        print(f"first exact depth k={first_exact}")

    print()


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states: list[State]) -> None:
    print("=" * 90)
    print("TEST 7: EXAMPLES")
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

        print(
            f"n={s.n} p={s.p} q={s.q} "
            f"frame={s.frame}"
        )

        print(f"    S={s.S}")
        print(f"    D={s.D}")
        print(f"    C={s.C}")
        print(f"    v2(S)={v2(s.S)}")
        print(f"    v2(D)={v2(s.D)}")
        print(f"    v2(C)={s.vc}")
        print(f"    v2(n+c)={s.wn}")
        print(f"    depth={s.depth}")

        if s.frame == "A":
            print(
                f"    2C=S-D+6 -> "
                f"{s.S}-{s.D}+6={2 * s.C}"
            )
        else:
            print(
                f"    2C=S+D+2 -> "
                f"{s.S}+{s.D}+2={2 * s.C}"
            )

        print()

    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 673 START")
    print("=" * 90)
    print()

    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve(PRIME_LIMIT)

    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    # -------------------------------------------------------------------------
    # Baseline
    # -------------------------------------------------------------------------

    total_failures += test_baseline(states)

    # -------------------------------------------------------------------------
    # Exact algebra
    # -------------------------------------------------------------------------

    total_failures += test_symmetric_pair_reconstruction(states)

    # -------------------------------------------------------------------------
    # Information boundary
    # -------------------------------------------------------------------------

    test_single_coordinate(states, "S")
    test_single_coordinate(states, "D")

    test_shifted_coordinate(states, "S")
    test_shifted_coordinate(states, "D")

    # -------------------------------------------------------------------------
    # Direct valuation reconstruction
    # -------------------------------------------------------------------------

    total_failures += test_valuation_from_one_coordinate(states)

    # -------------------------------------------------------------------------
    # Full depth reconstruction
    # -------------------------------------------------------------------------

    test_depth_formula_from_coordinate(states, "S")
    test_depth_formula_from_coordinate(states, "D")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Established:")
    print()
    print("    FRAME A:")
    print("        2C = S - D + 6")
    print()
    print("    FRAME B:")
    print("        2C = S + D + 2")
    print()

    print("Therefore:")
    print()
    print("    (S,D) -> C")
    print("    (S,D) -> depth")
    print()

    print("This experiment asks the stronger question:")
    print()
    print("    n + S")
    print("or:")
    print("    n + D")
    print()
    print("Is one symmetric coordinate enough?")
    print()

    print("The important distinction is:")
    print()
    print("    n alone")
    print("        loses the factor branch.")
    print()
    print("    n + one symmetric coordinate")
    print("        may retain enough branch information")
    print("        to reconstruct the 2-adic depth.")
    print()

    print(
        "A persistent ambiguity means both n and that "
        "single symmetric coordinate still leave multiple "
        "possible depths."
    )

    print()

    print(f"TOTAL TEST FAILURES={total_failures}")
    print()

    print("=" * 90)
    print("EXPERIMENT 673 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
