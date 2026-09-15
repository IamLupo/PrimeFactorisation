#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict
from math import isqrt


# ==============================================================================
# EXPERIMENT 608
# ==============================================================================
# FACTOR-SIDE 2-ADIC SURVIVAL EQUIVALENCE
#
# Existing invariant system:
#
#   FRAME A:
#       X = (q-p+6)/2
#       Y = (p+q)/2
#
#   FRAME B:
#       X = (3p-q+6)/2
#       Y = (3p+q)/2
#
# Proposed factor-side survival laws:
#
#   FRAME A:
#       alive at level z iff
#           p ==  3 mod 2^z
#           q == -3 mod 2^z
#
#   FRAME B:
#       alive at level z iff
#           p == -1 mod 2^z
#           q ==  3 mod 2^z
#
# Deepest-level predictions:
#
#   FRAME A:
#       1 + min(v2(p-3), v2(q+3))
#
#   FRAME B:
#       1 + min(v2(p+1), v2(q-3))
#
# This experiment compares those expressions against the
# previously established GLOBAL (FRAME,X,Y) representation.
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

LIMIT = 2_000_000
MAX_Z = 24
EXAMPLE_LIMIT = 20


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str
    X: int
    Y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def sieve(limit: int) -> bytearray:
    flags = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        flags[0] = 0
    if limit >= 1:
        flags[1] = 0

    root = isqrt(limit)

    for p in range(2, root + 1):
        if not flags[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        flags[start : limit + 1 : p] = b"\x00" * count

    return flags


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(limit: int, flags: bytearray) -> list[State]:
    primes = [
        p
        for p in range(3, limit + 1, 2)
        if flags[p]
    ]

    states: list[State] = []

    for i, p in enumerate(primes):
        max_q = limit // p

        for q in primes[i:]:
            if q > max_q:
                break

            n = p * q

            if n & 3 == 3:
                frame = "A"

                x_num = q - p + 6
                y_num = p + q
            else:
                frame = "B"

                x_num = 3 * p - q + 6
                y_num = 3 * p + q

            if x_num % 2 != 0:
                raise AssertionError(
                    f"non-integral X for n={n}: {x_num}"
                )

            if y_num % 2 != 0:
                raise AssertionError(
                    f"non-integral Y for n={n}: {y_num}"
                )

            X = x_num // 2
            Y = y_num // 2

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    X=X,
                    Y=Y,
                )
            )

    return states


# ==============================================================================
# 2-ADIC VALUATION
# ==============================================================================

def v2(value: int) -> int:
    if value == 0:
        return 10**9

    value = abs(value)
    count = 0

    while (value & 1) == 0:
        value >>= 1
        count += 1

    return count


# ==============================================================================
# BIT EXTRACTION
# ==============================================================================
# For this experiment we mostly use congruences, but this helper makes the
# intended bit position explicit.
#
# IMPORTANT:
# Python's right shift on negative integers is arithmetic. That is not what
# we want for a residue bit. Therefore extract the bit using abs(value).
#
# The factor-side congruence tests themselves remain exact modulo arithmetic.
# ==============================================================================

def abs_bit(value: int, bit: int) -> int:
    return (abs(value) >> bit) & 1


# ==============================================================================
# GLOBAL X,Y SURVIVAL
# ==============================================================================
# Level z is represented when:
#
#     2^(z-1) | X
#     2^(z-1) | Y
# ==============================================================================

def xy_alive(X: int, Y: int, z: int) -> bool:
    modulus = 1 << (z - 1)

    return (
        X % modulus == 0
        and Y % modulus == 0
    )


# ==============================================================================
# FACTOR-SIDE SURVIVAL
# ==============================================================================

def factor_alive(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> bool:

    modulus = 1 << z

    if frame == "A":
        return (
            (p - 3) % modulus == 0
            and
            (q + 3) % modulus == 0
        )

    if frame == "B":
        return (
            (p + 1) % modulus == 0
            and
            (q - 3) % modulus == 0
        )

    raise ValueError(
        f"unknown frame {frame!r}"
    )


# ==============================================================================
# FACTOR-SIDE LINEAR FORM
# ==============================================================================

def linear_factor_alive(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> bool:

    if frame == "A":

        c1 = p + q
        c2 = q - p + 6

    elif frame == "B":

        c1 = 3 * p + q
        c2 = 3 * p - q + 6

    else:
        raise ValueError(
            f"unknown frame {frame!r}"
        )

    return (
        c1 % (1 << z) == 0
        and
        c2 % (1 << (z + 1)) == 0
    )


# ==============================================================================
# FACTOR-SIDE DEPTH
# ==============================================================================

def factor_depth_formula(
    frame: str,
    p: int,
    q: int,
) -> int:

    if frame == "A":
        return (
            min(
                v2(p - 3),
                v2(q + 3),
            )
            + 1
        )

    if frame == "B":
        return (
            min(
                v2(p + 1),
                v2(q - 3),
            )
            + 1
        )

    raise ValueError(
        f"unknown frame {frame!r}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_factor_survival_equivalence(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 1: X,Y SURVIVAL == FACTOR-SIDE SURVIVAL"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for state in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            xy = xy_alive(
                state.X,
                state.Y,
                z,
            )

            factor = factor_alive(
                state.frame,
                state.p,
                state.q,
                z,
            )

            checked += 1

            if xy != factor:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"frame={state.frame} "
                        f"p={state.p} "
                        f"q={state.q} "
                        f"X={state.X} "
                        f"Y={state.Y} "
                        f"z={z} "
                        f"XY={xy} "
                        f"factor={factor}"
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_linear_equivalence(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 2: FACTOR CONGRUENCE == LINEAR-FORM CONDITION"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            direct = factor_alive(
                state.frame,
                state.p,
                state.q,
                z,
            )

            linear = linear_factor_alive(
                state.frame,
                state.p,
                state.q,
                z,
            )

            checked += 1

            if direct != linear:
                failures += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_deepest_factor_formula(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 3: FACTOR-SIDE DEEPEST LEVEL"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for state in states:

        xy_depth = (
            min(
                v2(state.X),
                v2(state.Y),
            )
            + 1
        )

        factor_depth = factor_depth_formula(
            state.frame,
            state.p,
            state.q,
        )

        checked += 1

        if xy_depth != factor_depth:

            failures += 1

            if shown < EXAMPLE_LIMIT:

                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"frame={state.frame} "
                    f"X,Y=({state.X},{state.Y}) "
                    f"xy_depth={xy_depth} "
                    f"factor_depth={factor_depth}"
                )

                shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================
# Report the exact active factor residues for all observed levels.
# ==============================================================================

def test_factor_residue_stream(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 4: ACTIVE FACTOR RESIDUE STREAM"
    )
    print("=" * 90)

    failures = 0

    active: dict[
        tuple[int, str],
        dict[str, set[int]],
    ] = defaultdict(
        lambda: {
            "p": set(),
            "q": set(),
        }
    )

    for state in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            if not xy_alive(
                state.X,
                state.Y,
                z,
            ):
                continue

            modulus = 1 << z

            actual_p = state.p % modulus
            actual_q = state.q % modulus

            if state.frame == "A":

                expected_p = 3 % modulus
                expected_q = (-3) % modulus

            else:

                expected_p = (-1) % modulus
                expected_q = 3 % modulus

            active[
                (z, state.frame)
            ]["p"].add(
                actual_p
            )

            active[
                (z, state.frame)
            ]["q"].add(
                actual_q
            )

            if (
                actual_p != expected_p
                or
                actual_q != expected_q
            ):

                failures += 1

    for key in sorted(active):

        values = active[key]

        print(
            f"    z={key[0]} "
            f"frame={key[1]} "
            f"p={sorted(values['p'])} "
            f"q={sorted(values['q'])}"
        )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================
# For an already-live parent at level z, test the NEXT level directly.
# ==============================================================================

def test_next_factor_bits(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: NEXT FACTOR-LEVEL SURVIVAL"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not xy_alive(
                state.X,
                state.Y,
                z,
            ):
                continue

            next_xy = xy_alive(
                state.X,
                state.Y,
                z + 1,
            )

            next_factor = factor_alive(
                state.frame,
                state.p,
                state.q,
                z + 1,
            )

            checked += 1

            if next_xy != next_factor:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={state.frame}"
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================
# Determine whether one factor's next bit is enough to predict survival,
# and whether both are needed.
# ==============================================================================

def test_minimal_factor_bits(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 6: MINIMAL FACTOR-SIDE BIT STATE"
    )
    print("=" * 90)

    maps: dict[
        str,
        defaultdict,
    ] = {
        "pbit": defaultdict(set),
        "qbit": defaultdict(set),
        "pqbit": defaultdict(set),
        "frame_pbit": defaultdict(set),
        "frame_qbit": defaultdict(set),
    }

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not xy_alive(
                state.X,
                state.Y,
                z,
            ):
                continue

            survive = xy_alive(
                state.X,
                state.Y,
                z + 1,
            )

            # We are looking at the factor information newly
            # exposed when moving from modulus 2^z to 2^(z+1).
            #
            # Relative to the required active residue, define
            # the residual quotient parity.

            if state.frame == "A":

                p_residual = (
                    state.p - 3
                )

                q_residual = (
                    state.q + 3
                )

            else:

                p_residual = (
                    state.p + 1
                )

                q_residual = (
                    state.q - 3
                )

            pbit = (
                abs_bit(
                    p_residual,
                    z,
                )
            )

            qbit = (
                abs_bit(
                    q_residual,
                    z,
                )
            )

            maps[
                "pbit"
            ][
                (z, pbit)
            ].add(
                survive
            )

            maps[
                "qbit"
            ][
                (z, qbit)
            ].add(
                survive
            )

            maps[
                "pqbit"
            ][
                (z, pbit, qbit)
            ].add(
                survive
            )

            maps[
                "frame_pbit"
            ][
                (z, state.frame, pbit)
            ].add(
                survive
            )

            maps[
                "frame_qbit"
            ][
                (z, state.frame, qbit)
            ].add(
                survive
            )

    for name, mapping in maps.items():

        ambiguous = sum(
            1
            for values in mapping.values()
            if len(values) > 1
        )

        print(
            f"    {name:<18} "
            f"ambiguous_states={ambiguous}"
        )


# ==============================================================================
# TEST 7
# ==============================================================================
# Determine the child residue from factor residual bits.
#
# A live parent satisfies:
#
#   A:
#       p = 3 + 2^z a
#       q = -3 + 2^z b
#
#   B:
#       p = -1 + 2^z a
#       q = 3 + 2^z b
#
# The residue modulo 2^(z+1) is determined by a+b modulo 2.
# ==============================================================================

def test_factor_child_residue_bits(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 7: FACTOR-SIDE CHILD RESIDUE BIT LAW"
    )
    print("=" * 90)

    mapping: defaultdict[
        tuple,
        set[int],
    ] = defaultdict(set)

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not xy_alive(
                state.X,
                state.Y,
                z,
            ):
                continue

            if state.frame == "A":

                a = (
                    (
                        state.p - 3
                    )
                    // (1 << z)
                ) & 1

                b = (
                    (
                        state.q + 3
                    )
                    // (1 << z)
                ) & 1

            else:

                a = (
                    (
                        state.p + 1
                    )
                    // (1 << z)
                ) & 1

                b = (
                    (
                        state.q - 3
                    )
                    // (1 << z)
                ) & 1

            residue = state.n % (
                1 << (z + 1)
            )

            mapping[
                (
                    z,
                    state.frame,
                    a,
                    b,
                )
            ].add(
                residue
            )

    ambiguous = 0

    for key, values in sorted(
        mapping.items()
    ):

        if len(values) > 1:
            ambiguous += 1

        if key[0] <= 8:

            print(
                f"    {key} -> "
                f"{sorted(values)}"
            )

    print(
        f"ambiguous states={ambiguous}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================
# Direct cross-check against the X/Y bit model.
# ==============================================================================

def test_xy_factor_crosswalk(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 8: X/Y <-> FACTOR RESIDUAL CROSSWALK"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not xy_alive(
                state.X,
                state.Y,
                z,
            ):
                continue

            xb = abs_bit(
                state.X,
                z - 1,
            )

            yb = abs_bit(
                state.Y,
                z - 1,
            )

            next_xy = xy_alive(
                state.X,
                state.Y,
                z + 1,
            )

            if state.frame == "A":

                a = (
                    (
                        state.p - 3
                    )
                    // (1 << z)
                ) & 1

                b = (
                    (
                        state.q + 3
                    )
                    // (1 << z)
                ) & 1

            else:

                a = (
                    (
                        state.p + 1
                    )
                    // (1 << z)
                ) & 1

                b = (
                    (
                        state.q - 3
                    )
                    // (1 << z)
                ) & 1

            factor_next = (
                a == 0
                and
                b == 0
            )

            checked += 1

            if next_xy != factor_next:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={state.frame} "
                        f"XYbits=({xb},{yb}) "
                        f"factorbits=({a},{b}) "
                        f"xy_next={next_xy} "
                        f"factor_next={factor_next}"
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# SYMBOLIC SUMMARY
# ==============================================================================

def print_symbolic_summary() -> None:

    print()
    print("=" * 90)
    print(
        "SYMBOLIC FACTOR-SIDE HIERARCHY"
    )
    print("=" * 90)

    print(
r"""
FRAME A

    X = (q-p+6)/2
    Y = (p+q)/2

    Therefore:

        2^(z-1) | X
        2^(z-1) | Y

    iff:

        p =  3 mod 2^z
        q = -3 mod 2^z.

    Hence:

        depth_A
          =
        1 + min(
            v2(p-3),
            v2(q+3)
        ).


FRAME B

    X = (3p-q+6)/2
    Y = (3p+q)/2

    Therefore:

        2^(z-1) | X
        2^(z-1) | Y

    iff:

        p = -1 mod 2^z
        q =  3 mod 2^z.

    Hence:

        depth_B
          =
        1 + min(
            v2(p+1),
            v2(q-3)
        ).


ACTIVE FACTOR RESIDUES

    FRAME A:

        p ≡  3 (mod 2^z)
        q ≡ -3 (mod 2^z)

    FRAME B:

        p ≡ -1 (mod 2^z)
        q ≡  3 (mod 2^z)


NEXT-LEVEL SURVIVAL

    At an already-live level z:

    FRAME A survives to z+1 iff:

        bit_z(p-3) = 0
        bit_z(q+3) = 0

    FRAME B survives to z+1 iff:

        bit_z(p+1) = 0
        bit_z(q-3) = 0


CHILD RESIDUE

    FRAME A:

        p = 3  + 2^z a
        q = -3 + 2^z b

        n = pq

        n mod 2^(z+1)
          =
        -9 + 2^z(a+b)

    FRAME B:

        p = -1 + 2^z a
        q =  3 + 2^z b

        n mod 2^(z+1)
          =
        -3 + 2^z(a+b)

    Therefore the child residue depends on:

        a XOR b

    because only the parity of a+b matters.

    This gives a second factor-side finite-state layer:

        factor residual bit pair
            ->
        child residue bit.
"""
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def show_examples(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    shown = 0

    for state in states:

        if shown >= EXAMPLE_LIMIT:
            break

        if state.frame == "A":

            f1 = state.p - 3
            f2 = state.q + 3

        else:

            f1 = state.p + 1
            f2 = state.q - 3

        depth_xy = (
            min(
                v2(state.X),
                v2(state.Y),
            )
            + 1
        )

        depth_factor = (
            min(
                v2(f1),
                v2(f2),
            )
            + 1
        )

        print()
        print(
            f"n={state.n} "
            f"p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    frame={state.frame}"
        )

        print(
            f"    X={state.X} "
            f"Y={state.Y}"
        )

        print(
            f"    factor residuals="
            f"({f1},{f2})"
        )

        print(
            f"    v2(X)={v2(state.X)} "
            f"v2(Y)={v2(state.Y)}"
        )

        print(
            f"    factor v2="
            f"({v2(f1)},{v2(f2)})"
        )

        print(
            f"    depth_xy={depth_xy} "
            f"depth_factor={depth_factor}"
        )

        for z in range(
            2,
            min(
                depth_factor + 1,
                MAX_Z + 1,
            ),
        ):

            modulus = 1 << z

            if state.frame == "A":

                expected_p = 3 % modulus
                expected_q = (-3) % modulus

            else:

                expected_p = (-1) % modulus
                expected_q = 3 % modulus

            print(
                f"    z={z:<2} "
                f"(p mod 2^z,q mod 2^z)="
                f"({state.p % modulus},"
                f"{state.q % modulus}) "
                f"expected="
                f"({expected_p},{expected_q})"
            )

        shown += 1


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print(
        "EXPERIMENT 608 START"
    )
    print("=" * 90)

    print()
    print(
        "[1] PRIME SIEVE"
    )

    flags = sieve(
        LIMIT
    )

    print()
    print(
        "[2] SEMIPRIME GENERATION"
    )

    states = generate_states(
        LIMIT,
        flags,
    )

    print(
        f"    semiprimes={len(states)}"
    )

    print()
    print(
        "[3] BUILD GLOBAL STATES"
    )

    print(
        f"    global states={len(states)}"
    )

    test_factor_survival_equivalence(
        states
    )

    test_linear_equivalence(
        states
    )

    test_deepest_factor_formula(
        states
    )

    test_factor_residue_stream(
        states
    )

    test_next_factor_bits(
        states
    )

    test_minimal_factor_bits(
        states
    )

    test_factor_child_residue_bits(
        states
    )

    test_xy_factor_crosswalk(
        states
    )

    print_symbolic_summary()

    show_examples(
        states
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 608 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()