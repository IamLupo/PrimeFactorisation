#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 609
# ==============================================================================
# EXACT FACTOR-SIDE 2-ADIC RECONSTRUCTION
#
# CENTRAL CORRECTION
#
# The global coordinates are:
#
#   FRAME A:
#       X = (q-p+6)/2
#       Y = (p+q)/2
#
#   FRAME B:
#       X = (3p-q+6)/2
#       Y = (3p+q)/2
#
# Level z means:
#
#       2^(z-1) | X
#       2^(z-1) | Y
#
# Therefore the equivalent factor conditions are:
#
#   FRAME A:
#       2^(z-1) | (p-3)
#       2^(z-1) | (q+3)
#
#   FRAME B:
#       2^(z-1) | (p+1)
#       2^(z-1) | (q-3)
#
# NOT 2^z.
#
# This experiment verifies:
#
#   1. Global X,Y survival
#      ==
#      factor-side congruences.
#
#   2. Exact deepest-level equality.
#
#   3. Exact same-N transition equality.
#
#   4. Factor-side next-bit survival.
#
#   5. Correct child-residue formulas.
#
#   6. The relation between X/Y bits and factor residual bits.
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

        flags[start:limit + 1:p] = b"\x00" * count

    return flags


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(
    limit: int,
    flags: bytearray,
) -> list[State]:

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

            if n % 4 == 3:

                frame = "A"

                X_num = q - p + 6
                Y_num = p + q

            else:

                frame = "B"

                X_num = 3 * p - q + 6
                Y_num = 3 * p + q

            if X_num % 2 != 0:
                raise AssertionError(
                    f"X not integral for n={n}"
                )

            if Y_num % 2 != 0:
                raise AssertionError(
                    f"Y not integral for n={n}"
                )

            X = X_num // 2
            Y = Y_num // 2

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

    out = 0

    while (value & 1) == 0:
        value >>= 1
        out += 1

    return out


# ==============================================================================
# LEVEL SURVIVAL FROM GLOBAL X,Y
# ==============================================================================

def xy_alive(
    X: int,
    Y: int,
    z: int,
) -> bool:

    modulus = 1 << (z - 1)

    return (
        X % modulus == 0
        and
        Y % modulus == 0
    )


# ==============================================================================
# CORRECTED FACTOR-SIDE SURVIVAL
# ==============================================================================
#
# IMPORTANT:
#
#       level z
#           =>
#       modulus = 2^(z-1)
#
# ==============================================================================

def factor_alive(
    frame: str,
    p: int,
    q: int,
    z: int,
) -> bool:

    modulus = 1 << (z - 1)

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
        f"unknown frame: {frame}"
    )


# ==============================================================================
# FACTOR DEPTH
# ==============================================================================

def factor_depth(
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
        f"unknown frame: {frame}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_survival_equivalence(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 1: X,Y SURVIVAL == CORRECTED FACTOR SURVIVAL"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for s in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            xy = xy_alive(
                s.X,
                s.Y,
                z,
            )

            factor = factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            )

            checked += 1

            if xy != factor:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"frame={s.frame}",
                        f"p={s.p}",
                        f"q={s.q}",
                        f"X={s.X}",
                        f"Y={s.Y}",
                        f"z={z}",
                        f"XY={xy}",
                        f"factor={factor}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_depth(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 2: EXACT DEPTH EQUALITY"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for s in states:

        xy_depth = (
            min(
                v2(s.X),
                v2(s.Y),
            )
            + 1
        )

        f_depth = factor_depth(
            s.frame,
            s.p,
            s.q,
        )

        checked += 1

        if xy_depth != f_depth:

            failures += 1

            if shown < EXAMPLE_LIMIT:

                print(
                    "    mismatch",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"X={s.X}",
                    f"Y={s.Y}",
                    f"xy_depth={xy_depth}",
                    f"factor_depth={f_depth}",
                )

                shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================
# Exact factor-side transition:
#
# alive(z+1)
# iff
# the z-th residual bit of both factor-side expressions is zero.
# ==============================================================================

def test_factor_next_bit(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 3: FACTOR-SIDE NEXT BIT SURVIVAL"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for s in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            ):
                continue

            next_actual = factor_alive(
                s.frame,
                s.p,
                s.q,
                z + 1,
            )

            if s.frame == "A":

                u = (
                    s.p - 3
                ) // (1 << (z - 1))

                v = (
                    s.q + 3
                ) // (1 << (z - 1))

            else:

                u = (
                    s.p + 1
                ) // (1 << (z - 1))

                v = (
                    s.q - 3
                ) // (1 << (z - 1))

            predicted = (
                (u & 1) == 0
                and
                (v & 1) == 0
            )

            checked += 1

            if predicted != next_actual:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                        f"frame={s.frame}",
                        f"u={u}",
                        f"v={v}",
                        f"pred={predicted}",
                        f"actual={next_actual}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================
# Same-N transition:
#
# if z and z+1 both exist, compare local coordinates.
# ==============================================================================

def test_same_n_halving(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 4: SAME-N COORDINATE HALVING"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for s in states:

        for z in range(
            2,
            MAX_Z,
        ):

            current = factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            )

            nxt = factor_alive(
                s.frame,
                s.p,
                s.q,
                z + 1,
            )

            if not (
                current
                and
                nxt
            ):
                continue

            k = 1 << (z - 1)
            k_next = 2 * k

            x = s.X // k
            y = s.Y // k

            x_next = s.X // k_next
            y_next = s.Y // k_next

            checked += 1

            if (
                x_next != x // 2
                or
                y_next != y // 2
                or
                x != 2 * x_next
                or
                y != 2 * y_next
            ):

                failures += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================
# Exact child residue from the global state.
#
# For z >= 3:
#
#   A:
#       r' = -9 + 2^z * Xbit
#
#   B:
#       r' = -3 + 2^z * Xbit
#
# For z = 2 the quadratic term survives:
#
#   A:
#       r' = 7 - 4*Ybit
#
#   B:
#       r' = 5 - 4*Ybit
# ==============================================================================

def global_child_residue(
    frame: str,
    Xbit: int,
    Ybit: int,
    z: int,
) -> int:

    modulus = 1 << (z + 1)

    if z == 2:

        if frame == "A":

            return (
                7 - 4 * Ybit
            ) % modulus

        return (
            5 - 4 * Ybit
        ) % modulus

    if frame == "A":

        return (
            -9
            + (1 << z) * Xbit
        ) % modulus

    return (
        -3
        + (1 << z) * Xbit
    ) % modulus


def test_child_residue(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: EXACT CHILD RESIDUE LAW"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for s in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            ):
                continue

            Xbit = (
                abs(s.X)
                >> (z - 1)
            ) & 1

            Ybit = (
                abs(s.Y)
                >> (z - 1)
            ) & 1

            predicted = global_child_residue(
                s.frame,
                Xbit,
                Ybit,
                z,
            )

            actual = s.n % (
                1 << (z + 1)
            )

            checked += 1

            if predicted != actual:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"frame={s.frame}",
                        f"z={z}",
                        f"Xbit={Xbit}",
                        f"Ybit={Ybit}",
                        f"pred={predicted}",
                        f"actual={actual}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================
# Factor-side child residue using live factor-side residuals.
#
# At level z, define:
#
#   A:
#       a = (p-3)/2^(z-1)
#       b = (q+3)/2^(z-1)
#
#   B:
#       a = (p+1)/2^(z-1)
#       b = (q-3)/2^(z-1)
#
# Then modulo 2^(z+1):
#
#   A:
#       n = -9 + 2^z(a+b)
#
#   B:
#       n = -3 + 2^z(a+b)
#
# ==============================================================================

def factor_child_residue(
    frame: str,
    a: int,
    b: int,
    z: int,
) -> int:

    modulus = 1 << (z + 1)

    if frame == "A":

        return (
            -9
            + (1 << z) * (
                (a + b) & 1
            )
        ) % modulus

    return (
        -3
        + (1 << z) * (
            (a + b) & 1
        )
    ) % modulus


def test_factor_child_residue(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 6: FACTOR RESIDUAL -> CHILD RESIDUE"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for s in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            ):
                continue

            denominator = 1 << (
                z - 1
            )

            if s.frame == "A":

                a = (
                    s.p - 3
                ) // denominator

                b = (
                    s.q + 3
                ) // denominator

            else:

                a = (
                    s.p + 1
                ) // denominator

                b = (
                    s.q - 3
                ) // denominator

            predicted = factor_child_residue(
                s.frame,
                a,
                b,
                z,
            )

            actual = s.n % (
                1 << (z + 1)
            )

            checked += 1

            if predicted != actual:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                        f"frame={s.frame}",
                        f"a={a}",
                        f"b={b}",
                        f"pred={predicted}",
                        f"actual={actual}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================
# Compare the factor residual bits to the X,Y bits.
#
# For A:
#
#   p-3 = 2^(z-1) * (Y-X)/? ...
#
# The experiment derives the exact parity map empirically and prints
# the deterministic tables.
# ==============================================================================

def test_bit_crosswalk(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 7: FACTOR RESIDUAL BIT <-> X,Y BIT CROSSWALK"
    )
    print("=" * 90)

    tables: dict[
        tuple,
        set[tuple[int, int]],
    ] = defaultdict(set)

    for s in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            ):
                continue

            Xbit = (
                abs(s.X)
                >> (z - 1)
            ) & 1

            Ybit = (
                abs(s.Y)
                >> (z - 1)
            ) & 1

            denominator = 1 << (
                z - 1
            )

            if s.frame == "A":

                a = (
                    (
                        s.p - 3
                    ) // denominator
                ) & 1

                b = (
                    (
                        s.q + 3
                    ) // denominator
                ) & 1

            else:

                a = (
                    (
                        s.p + 1
                    ) // denominator
                ) & 1

                b = (
                    (
                        s.q - 3
                    ) // denominator
                ) & 1

            tables[
                (
                    z,
                    s.frame,
                    Xbit,
                    Ybit,
                )
            ].add(
                (a, b)
            )

    ambiguous = 0

    for key in sorted(tables):

        values = tables[key]

        if len(values) > 1:
            ambiguous += 1

        if key[0] <= 8:

            print(
                f"    {key} -> "
                f"{sorted(values)}"
            )

    print(
        f"ambiguous={ambiguous}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================
# Active residue itself.
# ==============================================================================

def test_active_residue(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 8: ACTIVE RESIDUE"
    )
    print("=" * 90)

    values: defaultdict[
        tuple,
        set[int],
    ] = defaultdict(set)

    for s in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            if not factor_alive(
                s.frame,
                s.p,
                s.q,
                z,
            ):
                continue

            values[
                (z, s.frame)
            ].add(
                s.n % (1 << z)
            )

    failures = 0

    for key in sorted(values):

        vals = values[key]

        z, frame = key

        if frame == "A":
            expected = (-9) % (1 << z)
        else:
            expected = (-3) % (1 << z)

        if vals != {expected}:

            failures += 1

        print(
            f"    z={z:<2} "
            f"frame={frame} "
            f"observed={sorted(vals)} "
            f"expected={expected}"
        )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# TEST 9
# ==============================================================================

def test_complete_reconstruction(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 9: COMPLETE RECONSTRUCTION"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for s in states:

        if s.frame == "A":

            p2 = (
                s.Y
                - s.X
                + 3
            )

            q2 = (
                s.Y
                + s.X
                - 3
            )

        else:

            numerator = (
                s.Y
                + s.X
                - 3
            )

            if numerator % 3 != 0:
                failures += 1
                continue

            p2 = numerator // 3

            q2 = (
                s.Y
                - s.X
                + 3
            )

        n2 = p2 * q2

        checked += 1

        if (
            p2 != s.p
            or
            q2 != s.q
            or
            n2 != s.n
        ):

            failures += 1

    print(
        f"checked={checked} failures={failures}"
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

    count = 0

    for s in states:

        if count >= EXAMPLE_LIMIT:
            break

        xy_depth = (
            min(
                v2(s.X),
                v2(s.Y),
            )
            + 1
        )

        factor_depth_value = factor_depth(
            s.frame,
            s.p,
            s.q,
        )

        print()
        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q}"
        )

        print(
            f"    frame={s.frame}"
        )

        print(
            f"    X={s.X} "
            f"Y={s.Y}"
        )

        print(
            f"    v2(X)={v2(s.X)} "
            f"v2(Y)={v2(s.Y)}"
        )

        print(
            f"    xy_depth={xy_depth}"
        )

        print(
            f"    factor_depth={factor_depth_value}"
        )

        for z in range(
            2,
            min(
                xy_depth + 1,
                MAX_Z + 1,
            ),
        ):

            modulus = 1 << (z - 1)

            if s.frame == "A":

                rp = (
                    (s.p - 3)
                    // modulus
                )

                rq = (
                    (s.q + 3)
                    // modulus
                )

            else:

                rp = (
                    (s.p + 1)
                    // modulus
                )

                rq = (
                    (s.q - 3)
                    // modulus
                )

            print(
                f"    z={z:<2} "
                f"modulus=2^{z-1:<2} "
                f"residuals=({rp},{rq}) "
                f"bits=({rp & 1},{rq & 1}) "
                f"n mod 2^z="
                f"{s.n % (1 << z)}"
            )

        count += 1


# ==============================================================================
# SYMBOLIC SUMMARY
# ==============================================================================

def symbolic_summary() -> None:

    print()
    print("=" * 90)
    print(
        "SYMBOLIC SUMMARY"
    )
    print("=" * 90)

    print(
r"""
The correct correspondence is:

FRAME A

    X = (q-p+6)/2
    Y = (p+q)/2

    2^(z-1) | X,Y

    iff

    p-3 = 0 mod 2^(z-1)
    q+3 = 0 mod 2^(z-1).


FRAME B

    X = (3p-q+6)/2
    Y = (3p+q)/2

    2^(z-1) | X,Y

    iff

    p+1 = 0 mod 2^(z-1)
    q-3 = 0 mod 2^(z-1).


The previous factor-side experiment used 2^z here,
which shifted the predicted depth by exactly one level.

For an already-live level z define:

FRAME A:

    a = (p-3)/2^(z-1)
    b = (q+3)/2^(z-1)

FRAME B:

    a = (p+1)/2^(z-1)
    b = (q-3)/2^(z-1)

Then:

    z -> z+1 survives
        iff
    a even AND b even.

The next residue is:

FRAME A:

    n mod 2^(z+1)
       =
    -9 + 2^z(a+b)

FRAME B:

    n mod 2^(z+1)
       =
    -3 + 2^z(a+b)

so the new residue bit is:

    a XOR b.

Thus the corrected hierarchy contains:

    GLOBAL STATE:
        (FRAME, X, Y)

    LEVEL SURVIVAL:
        Xbit, Ybit

    FACTOR-SIDE SURVIVAL:
        residual bit a, residual bit b

    ACTIVE RESIDUE:
        FRAME

    CHILD RESIDUE BIT:
        a XOR b

The key question now is whether the factor residual bits
are themselves a simple linear transformation of the
current X,Y bits.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print(
        "EXPERIMENT 609 START"
    )
    print("=" * 90)

    print()
    print("[1] PRIME SIEVE")

    flags = sieve(
        LIMIT
    )

    print()
    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        LIMIT,
        flags,
    )

    print(
        f"    semiprimes={len(states)}"
    )

    print()
    print("[3] GLOBAL STATES")

    print(
        f"    global states={len(states)}"
    )

    test_survival_equivalence(
        states
    )

    test_depth(
        states
    )

    test_factor_next_bit(
        states
    )

    test_same_n_halving(
        states
    )

    test_child_residue(
        states
    )

    test_factor_child_residue(
        states
    )

    test_bit_crosswalk(
        states
    )

    test_active_residue(
        states
    )

    test_complete_reconstruction(
        states
    )

    symbolic_summary()

    show_examples(
        states
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 609 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
