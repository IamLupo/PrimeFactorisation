#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 585
# ==============================================================================
# HIERARCHICAL CHILD-BRANCH / COMPRESSED-x TEST
#
# No files.
# No web.
#
# Definitions:
#
#     k = 2^(z-1)
#
#     u = y-x
#     v = y+x
#
# Case A:
#
#     p = ku + 3
#     q = kv - 3
#
# Case B:
#
#     3p = kv - 3
#     q  = ku + 3
#
#
# Since:
#
#     v-u = 2x
#
# the next modular branch can potentially depend on:
#
#     x mod 2
#
# rather than on the full x,y state.
#
# After consuming that bit we test:
#
#     x_next = floor(x/2)
#
# and ask whether the same process repeats.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 9

SHOW_EXAMPLES = True
MAX_EXAMPLES = 10


# ==============================================================================
# SAMPLE
# ==============================================================================

@dataclass(frozen=True)
class Sample:
    n: int
    p: int
    q: int

    z: int
    k: int
    modulus: int

    branch: str

    x: int
    y: int

    u: int
    v: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> bytearray:

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = isqrt(limit)

    for p in range(2, root + 1):

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
# SEMIPRIMES
# ==============================================================================

def generate_semiprimes(
    max_n: int,
    sieve: bytearray,
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
# CASE A INVERSE
# ==============================================================================

def solve_case_a(
    p: int,
    q: int,
    k: int,
):

    # p = ku + 3
    # q = kv - 3

    if (p - 3) % k != 0:
        return None

    if (q + 3) % k != 0:
        return None

    u = (
        p - 3
    ) // k

    v = (
        q + 3
    ) // k

    if (v - u) % 2 != 0:
        return None

    if (v + u) % 2 != 0:
        return None

    x = (
        v - u
    ) // 2

    y = (
        v + u
    ) // 2

    return x, y, u, v


# ==============================================================================
# CASE B INVERSE
# ==============================================================================

def solve_case_b(
    p: int,
    q: int,
    k: int,
):

    # 3p = kv - 3
    # q  = ku + 3

    if (3 * p + 3) % k != 0:
        return None

    if (q - 3) % k != 0:
        return None

    v = (
        3 * p + 3
    ) // k

    u = (
        q - 3
    ) // k

    if (v - u) % 2 != 0:
        return None

    if (v + u) % 2 != 0:
        return None

    x = (
        v - u
    ) // 2

    y = (
        v + u
    ) // 2

    return x, y, u, v


# ==============================================================================
# LEVEL SAMPLE
# ==============================================================================

def make_sample(
    n: int,
    p: int,
    q: int,
    z: int,
):

    k = 1 << (z - 1)
    modulus = 1 << z
    residue = n % modulus

    # --------------------------------------------------------------------------
    # Case A
    # --------------------------------------------------------------------------

    if residue == modulus - 1:

        solved = solve_case_a(
            p,
            q,
            k,
        )

        if solved is None:
            return None

        x, y, u, v = solved

        check_p = (
            k * u + 3
        )

        check_q = (
            k * v - 3
        )

        branch = "A"

    # --------------------------------------------------------------------------
    # Case B
    # --------------------------------------------------------------------------

    else:

        solved = solve_case_b(
            p,
            q,
            k,
        )

        if solved is None:
            return None

        x, y, u, v = solved

        numerator = (
            k * v - 3
        )

        if numerator % 3 != 0:
            return None

        check_p = numerator // 3

        check_q = (
            k * u + 3
        )

        branch = "B"

    # --------------------------------------------------------------------------
    # Exact validation
    # --------------------------------------------------------------------------

    if check_p != p:
        raise RuntimeError(
            "p reconstruction failed"
        )

    if check_q != q:
        raise RuntimeError(
            "q reconstruction failed"
        )

    if p * q != n:
        raise RuntimeError(
            "n reconstruction failed"
        )

    if v - u != 2 * x:
        raise RuntimeError(
            "v-u != 2x"
        )

    if v + u != 2 * y:
        raise RuntimeError(
            "v+u != 2y"
        )

    return Sample(
        n=n,
        p=p,
        q=q,
        z=z,
        k=k,
        modulus=modulus,
        branch=branch,
        x=x,
        y=y,
        u=u,
        v=v,
    )


# ==============================================================================
# BUILD ALL LEVEL DATA
# ==============================================================================

def generate_level_data(
    semiprimes,
):

    data = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        for n, p, q in semiprimes:

            sample = make_sample(
                n,
                p,
                q,
                z,
            )

            if sample is not None:
                data[z].append(sample)

    return data


# ==============================================================================
# NEXT BIT
# ==============================================================================

def next_bit(
    n: int,
    z: int,
) -> int:

    return (
        n >> z
    ) & 1


# ==============================================================================
# CHILD RESIDUE
# ==============================================================================

def child_residue(
    n: int,
    z: int,
) -> int:

    return n % (
        1 << (z + 1)
    )


# ==============================================================================
# TEST x MOD 2
# ==============================================================================

def test_x_parity(
    data,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: CHILD BRANCH FROM x MOD 2"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"LEVEL z={z} "
            f"(mod {1 << z} -> mod {1 << (z + 1)})"
        )

        for branch in (
            "A",
            "B",
        ):

            subset = [
                s
                for s in data[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            mapping = defaultdict(set)

            for sample in subset:

                parity = sample.x & 1
                bit = next_bit(
                    sample.n,
                    z,
                )

                mapping[
                    parity
                ].add(bit)

            exact = all(
                len(values) == 1
                for values in mapping.values()
            )

            print(
                f"    branch={branch} "
                f"exact={exact}"
            )

            for parity in sorted(mapping):

                print(
                    f"        x mod 2 = {parity} "
                    f"-> bits = "
                    f"{sorted(mapping[parity])}"
                )


# ==============================================================================
# TEST FULL x VALUE
# ==============================================================================

def test_x_resolution(
    data,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: MINIMUM x RESOLUTION"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"z={z}"
        )

        for branch in (
            "A",
            "B",
        ):

            subset = [
                s
                for s in data[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            found = None

            for level in range(
                1,
                z + 2,
            ):

                modulus = 1 << level

                mapping = defaultdict(set)

                for sample in subset:

                    key = (
                        sample.x
                        % modulus
                    )

                    mapping[key].add(
                        next_bit(
                            sample.n,
                            z,
                        )
                    )

                exact = all(
                    len(values) == 1
                    for values in mapping.values()
                )

                if exact:

                    found = (
                        level,
                        modulus,
                        mapping,
                    )

                    break

            if found is None:

                print(
                    f"    branch={branch}: "
                    "NO EXACT x RULE"
                )

            else:

                level, modulus, mapping = found

                print(
                    f"    branch={branch}: "
                    f"x mod {modulus} "
                    f"(2^{level}) is exact"
                )

                print(
                    f"        states={len(mapping)}"
                )


# ==============================================================================
# SYMBOLIC FORM
# ==============================================================================

def symbolic_transition():

    print()
    print("=" * 90)
    print(
        "SYMBOLIC CHILD-BIT TRANSITION"
    )
    print("=" * 90)

    for z in range(
        3,
        MAX_Z,
    ):

        modulus_child = 1 << (
            z + 1
        )

        power = 1 << z

        print()
        print(
            f"z={z}"
        )

        # ----------------------------------------------------------------------
        # A
        # ----------------------------------------------------------------------

        a_even = (
            -9
        ) % modulus_child

        a_odd = (
            -9
            + 3 * power
        ) % modulus_child

        print(
            "    CASE A:"
        )

        print(
            f"        x even -> "
            f"n mod 2^{z+1} = "
            f"{a_even}"
        )

        print(
            f"        x odd  -> "
            f"n mod 2^{z+1} = "
            f"{a_odd}"
        )

        print(
            f"        next bit even = "
            f"{(a_even >> z) & 1}"
        )

        print(
            f"        next bit odd  = "
            f"{(a_odd >> z) & 1}"
        )

        # ----------------------------------------------------------------------
        # B
        # ----------------------------------------------------------------------

        b_even = (
            -3
        ) % modulus_child

        b_odd = (
            -3
            + power
        ) % modulus_child

        print(
            "    CASE B:"
        )

        print(
            f"        x even -> "
            f"n mod 2^{z+1} = "
            f"{b_even}"
        )

        print(
            f"        x odd  -> "
            f"n mod 2^{z+1} = "
            f"{b_odd}"
        )

        print(
            f"        next bit even = "
            f"{(b_even >> z) & 1}"
        )

        print(
            f"        next bit odd  = "
            f"{(b_odd >> z) & 1}"
        )


# ==============================================================================
# RECURSIVE x COMPRESSION
# ==============================================================================

def test_recursive_halving(
    data,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: RECURSIVE x -> floor(x/2)"
    )
    print("=" * 90)

    print()
    print(
        "The current bit consumes:"
    )

    print(
        "    x mod 2"
    )

    print()
    print(
        "Then define:"
    )

    print(
        "    x_next = floor(x/2)"
    )

    print()
    print(
        "We test whether x_next mod 2 "
        "controls the subsequent structure."
    )

    for z in range(
        MIN_Z,
        MAX_Z - 1,
    ):

        print()
        print(
            f"LEVEL z={z}"
        )

        for branch in (
            "A",
            "B",
        ):

            subset = [
                s
                for s in data[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            mapping = defaultdict(set)

            for sample in subset:

                x_next = sample.x // 2

                mapping[
                    x_next & 1
                ].add(
                    next_bit(
                        sample.n,
                        z,
                    )
                )

            exact = all(
                len(values) == 1
                for values in mapping.values()
            )

            print(
                f"    branch={branch} "
                f"x_next=x//2 "
                f"exact={exact}"
            )

            for key in sorted(mapping):

                print(
                    f"        x_next mod 2={key} "
                    f"bits={sorted(mapping[key])}"
                )


# ==============================================================================
# PARENT / CHILD SIGNATURES
# ==============================================================================

def branch_signature(
    data,
):

    print()
    print("=" * 90)
    print(
        "PARENT -> CHILD SIGNATURES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"z={z}"
        )

        for branch in (
            "A",
            "B",
        ):

            subset = [
                s
                for s in data[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            signature = defaultdict(
                set
            )

            for sample in subset:

                state = (
                    sample.x & 1
                )

                signature[state].add(
                    child_residue(
                        sample.n,
                        z,
                    )
                )

            print(
                f"    branch={branch}"
            )

            for state in sorted(signature):

                print(
                    f"        x%2={state} "
                    f"children="
                    f"{sorted(signature[state])}"
                )


# ==============================================================================
# DOMAIN REPORT
# ==============================================================================

def domain_report(
    data,
):

    print()
    print("=" * 90)
    print(
        "NORMALIZED x DOMAIN"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print()
        print(
            f"z={z}"
        )

        for branch in (
            "A",
            "B",
        ):

            subset = [
                s
                for s in data[z]
                if s.branch == branch
            ]

            if not subset:
                continue

            xs = {
                s.x
                for s in subset
            }

            even = sum(
                (s.x & 1) == 0
                for s in subset
            )

            odd = len(subset) - even

            print(
                f"    branch={branch} "
                f"states={len(subset):<7} "
                f"unique_x={len(xs):<7} "
                f"even={even:<7} "
                f"odd={odd:<7} "
                f"range=[{min(xs)}, {max(xs)}]"
            )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    data,
):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print()
        print(
            f"z={z}"
        )

        for sample in data[z][
            :MAX_EXAMPLES
        ]:

            print(
                f"    n={sample.n} "
                f"p={sample.p} "
                f"q={sample.q} "
                f"x={sample.x} "
                f"y={sample.y} "
                f"x%2={sample.x & 1} "
                f"next_bit="
                f"{next_bit(sample.n, z)} "
                f"branch={sample.branch}"
            )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 585 START"
    )
    print("=" * 90)

    print()
    print(
        "HIERARCHICAL CHILD-BRANCH / COMPRESSED-x TEST"
    )

    # --------------------------------------------------------------------------
    # Generate.
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
        f"    semiprimes={len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Build all levels.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] LEVEL DATA"
    )

    data = generate_level_data(
        semiprimes
    )

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print(
            f"    z={z}: "
            f"{len(data[z])} states"
        )

    # --------------------------------------------------------------------------
    # Symbolic derivation.
    # --------------------------------------------------------------------------

    symbolic_transition()

    # --------------------------------------------------------------------------
    # x parity.
    # --------------------------------------------------------------------------

    test_x_parity(
        data
    )

    # --------------------------------------------------------------------------
    # Minimum x resolution.
    # --------------------------------------------------------------------------

    test_x_resolution(
        data
    )

    # --------------------------------------------------------------------------
    # Recursive halving.
    # --------------------------------------------------------------------------

    test_recursive_halving(
        data
    )

    # --------------------------------------------------------------------------
    # Signatures.
    # --------------------------------------------------------------------------

    branch_signature(
        data
    )

    # --------------------------------------------------------------------------
    # Domain.
    # --------------------------------------------------------------------------

    domain_report(
        data
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    examples(
        data
    )

    # --------------------------------------------------------------------------
    # Final.
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "FINAL AUDIT"
    )
    print("=" * 90)

    print(
        """
The previous parity experiment used:

    (v-u) mod 2

but:

    v-u = 2x

so that test discarded the important bit.

The corrected hierarchy is:

    v-u
      =
    2x
      ->
    x mod 2
      ->
    child bit

After consuming x mod 2, the remaining coordinate is:

    x_next = floor(x/2)

which is tested for recursive behavior.

A strong result would therefore look like:

    parent branch
        +
    x mod 2
        ->
    child branch

followed by:

    x -> floor(x/2)

and the same mechanism repeating at the next level.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 585 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()