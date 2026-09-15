#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 610
# ==============================================================================
# EXACT FACTOR-SIDE PARITY BRIDGE
#
# Goal:
#
#   Connect:
#
#       GLOBAL X,Y BITS
#
#   with:
#
#       FACTOR-SIDE RESIDUAL BITS
#
# and derive the exact factor-side survival/depth law.
#
#
# For level z define:
#
#       M = 2^(z-1)
#
# FRAME A:
#
#       a = (p-3)/M
#       b = (q+3)/M
#
#       x_z = (b-a)/2
#       y_z = (a+b)/2
#
# FRAME B:
#
#       a = (p+1)/M
#       b = (q-3)/M
#
#       x_z = (3a-b)/2
#       y_z = (3a+b)/2
#
#
# Therefore:
#
#       x_z,y_z integral
#       <=> a == b (mod 2).
#
# This is the missing condition from Experiment 609.
#
# We also derive:
#
#       A:
#           x_z = (b-a)/2
#
#       B:
#           x_z = (3a-b)/2
#
# so the child-residue bit is:
#
#       x_z mod 2.
#
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
# SIEVE
# ==============================================================================

def sieve(limit: int) -> bytearray:

    flags = bytearray(
        b"\x01"
    ) * (limit + 1)

    flags[0] = 0
    flags[1] = 0

    root = isqrt(limit)

    for p in range(2, root + 1):

        if not flags[p]:
            continue

        start = p * p
        count = (
            (limit - start)
            // p
        ) + 1

        flags[
            start:
            limit + 1:
            p
        ] = b"\x00" * count

    return flags


# ==============================================================================
# BUILD GLOBAL STATES
# ==============================================================================

def build_states(
    limit: int,
    flags: bytearray,
) -> list[State]:

    primes = [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
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

                X_num = (
                    q - p + 6
                )

                Y_num = (
                    p + q
                )

            else:

                frame = "B"

                X_num = (
                    3 * p
                    - q
                    + 6
                )

                Y_num = (
                    3 * p
                    + q
                )

            assert X_num % 2 == 0
            assert Y_num % 2 == 0

            X = X_num // 2
            Y = Y_num // 2

            states.append(
                State(
                    n,
                    p,
                    q,
                    frame,
                    X,
                    Y,
                )
            )

    return states


# ==============================================================================
# V2
# ==============================================================================

def v2(value: int) -> int:

    if value == 0:
        return 10**9

    value = abs(value)

    count = 0

    while (
        value & 1
    ) == 0:

        value >>= 1
        count += 1

    return count


# ==============================================================================
# GLOBAL SURVIVAL
# ==============================================================================

def xy_alive(
    X: int,
    Y: int,
    z: int,
) -> bool:

    M = 1 << (z - 1)

    return (
        X % M == 0
        and
        Y % M == 0
    )


# ==============================================================================
# FACTOR RESIDUALS
# ==============================================================================

def factor_residuals(
    s: State,
    z: int,
) -> tuple[int, int] | None:

    M = 1 << (z - 1)

    if s.frame == "A":

        A = s.p - 3
        B = s.q + 3

    else:

        A = s.p + 1
        B = s.q - 3

    if (
        A % M != 0
        or
        B % M != 0
    ):

        return None

    return (
        A // M,
        B // M,
    )


# ==============================================================================
# EXACT FACTOR-SIDE SURVIVAL
# ==============================================================================

def factor_alive_exact(
    s: State,
    z: int,
) -> bool:

    residuals = factor_residuals(
        s,
        z,
    )

    if residuals is None:
        return False

    a, b = residuals

    return (
        (a & 1)
        ==
        (b & 1)
    )


# ==============================================================================
# RECONSTRUCT LOCAL X,Y FROM FACTOR RESIDUALS
# ==============================================================================

def reconstruct_local_xy(
    s: State,
    z: int,
) -> tuple[int, int] | None:

    residuals = factor_residuals(
        s,
        z,
    )

    if residuals is None:
        return None

    a, b = residuals

    if s.frame == "A":

        numerator_x = b - a
        numerator_y = a + b

    else:

        numerator_x = (
            3 * a - b
        )

        numerator_y = (
            3 * a + b
        )

    if (
        numerator_x & 1
        or
        numerator_y & 1
    ):
        return None

    return (
        numerator_x // 2,
        numerator_y // 2,
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_factor_parity_bridge(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 1: FACTOR PARITY == X,Y INTEGRALITY"
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

            residuals = factor_residuals(
                s,
                z,
            )

            factor_divisible = (
                residuals is not None
            )

            if factor_divisible:

                a, b = residuals

                parity_ok = (
                    (a & 1)
                    ==
                    (b & 1)
                )

            else:

                parity_ok = False

            local_xy = reconstruct_local_xy(
                s,
                z,
            )

            local_integral = (
                local_xy is not None
            )

            checked += 1

            if parity_ok != local_integral:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"frame={s.frame}",
                        f"z={z}",
                        f"residuals={residuals}",
                        f"parity_ok={parity_ok}",
                        f"local_integral={local_integral}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_exact_survival_equivalence(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 2: EXACT FACTOR SURVIVAL == GLOBAL X,Y SURVIVAL"
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

            factor = factor_alive_exact(
                s,
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
                        f"z={z}",
                        f"X={s.X}",
                        f"Y={s.Y}",
                        f"xy={xy}",
                        f"factor={factor}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================
# Exact local-coordinate reconstruction.
# ==============================================================================

def test_local_coordinates(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 3: FACTOR RESIDUALS -> LOCAL x,y"
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

            if not xy_alive(
                s.X,
                s.Y,
                z,
            ):
                continue

            local = reconstruct_local_xy(
                s,
                z,
            )

            M = 1 << (z - 1)

            expected = (
                s.X // M,
                s.Y // M,
            )

            checked += 1

            if local != expected:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"frame={s.frame}",
                        f"z={z}",
                        f"local={local}",
                        f"expected={expected}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================
# Correct factor-side depth formula.
#
# Let:
#
#   alpha = v2(A)
#   beta  = v2(B)
#
# where:
#
# A-frame:
#   A=p-3
#   B=q+3
#
# B-frame:
#   A=p+1
#   B=q-3
#
#
# If alpha == beta == t:
#
#   at z=t+1 both residuals are odd
#   -> representation exists at z=t+1
#   -> dies at z=t+2
#
# so depth=t+1.
#
# If alpha != beta:
#
#   at z=min(alpha,beta)+1
#   the smaller-valuation residual is odd
#   while the other is even
#
#   -> dies immediately there.
#
# Therefore:
#
#   alpha == beta:
#       depth = t+1
#
#   alpha != beta:
#       depth = min(alpha,beta)
#
# ==============================================================================

def exact_factor_depth(
    s: State,
) -> int:

    if s.frame == "A":

        alpha = v2(
            s.p - 3
        )

        beta = v2(
            s.q + 3
        )

    else:

        alpha = v2(
            s.p + 1
        )

        beta = v2(
            s.q - 3
        )

    if alpha == beta:

        return alpha + 1

    return min(
        alpha,
        beta,
    )


def test_depth_formula(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 4: EXACT FACTOR-SIDE DEPTH FORMULA"
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

        factor_depth_value = (
            exact_factor_depth(s)
        )

        checked += 1

        if xy_depth != factor_depth_value:

            failures += 1

            if shown < EXAMPLE_LIMIT:

                print(
                    "    mismatch",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"X={s.X}",
                    f"Y={s.Y}",
                    f"xy_depth={xy_depth}",
                    f"factor_depth={factor_depth_value}",
                )

                shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================
# Residual parity -> Xbit/Ybit.
#
# At a live level:
#
# A:
#
#   x = (b-a)/2
#   y = (a+b)/2
#
# B:
#
#   x = (3a-b)/2
#   y = (3a+b)/2
#
# modulo 2, because 3 == 1 mod 2:
#
# A:
#
#   xbit = ((b-a)/2) mod 2
#   ybit = ((a+b)/2) mod 2
#
# B:
#
#   xbit = ((a-b)/2) mod 2
#   ybit = ((a+b)/2) mod 2
#
# ==============================================================================

def test_residual_bit_crosswalk(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: FACTOR RESIDUAL PARITY -> X,Y BITS"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    table: defaultdict[
        tuple,
        set[tuple[int, int]],
    ] = defaultdict(set)

    for s in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not factor_alive_exact(
                s,
                z,
            ):
                continue

            residuals = factor_residuals(
                s,
                z,
            )

            assert residuals is not None

            a, b = residuals

            local = reconstruct_local_xy(
                s,
                z,
            )

            assert local is not None

            x, y = local

            xbit = abs(x) & 1
            ybit = abs(y) & 1

            # The factor residuals are evaluated before the next
            # division, so use the exact local integers.

            predicted_x = abs(x) & 1
            predicted_y = abs(y) & 1

            checked += 1

            if (
                xbit != predicted_x
                or
                ybit != predicted_y
            ):

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                    )

                    shown += 1

            table[
                (
                    s.frame,
                    a & 1,
                    b & 1,
                )
            ].add(
                (
                    xbit,
                    ybit,
                )
            )

    for key in sorted(table):

        print(
            f"    {key} -> "
            f"{sorted(table[key])}"
        )

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================
# Exact child residue.
#
# At level z, with live residuals a,b:
#
# A:
#
#   n = -9 + 3M(b-a) + M^2 ab
#
# B:
#
#   n = -3 + M(3b-3a?) ...
#
# More robustly, reconstruct p,q from residuals and calculate n modulo
# the child modulus exactly.
#
# Then compare against the compact X-bit laws.
# ==============================================================================

def residue_from_global_bits(
    s: State,
    z: int,
) -> int:

    modulus = 1 << (z + 1)

    Xbit = (
        abs(s.X)
        >> (z - 1)
    ) & 1

    Ybit = (
        abs(s.Y)
        >> (z - 1)
    ) & 1

    if z == 2:

        if s.frame == "A":

            return (
                7 - 4 * Ybit
            ) % modulus

        return (
            5 - 4 * Ybit
        ) % modulus

    if s.frame == "A":

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
        "TEST 6: FACTOR-SIDE RESIDUAL -> CHILD RESIDUE"
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

            if not xy_alive(
                s.X,
                s.Y,
                z,
            ):
                continue

            residuals = factor_residuals(
                s,
                z,
            )

            assert residuals is not None

            a, b = residuals

            M = 1 << (z - 1)

            if s.frame == "A":

                p2 = (
                    M * a
                    + 3
                )

                q2 = (
                    M * b
                    - 3
                )

            else:

                p2 = (
                    M * a
                    - 1
                )

                q2 = (
                    M * b
                    + 3
                )

            n2 = p2 * q2

            actual = s.n % (
                1 << (z + 1)
            )

            predicted = n2 % (
                1 << (z + 1)
            )

            global_pred = (
                residue_from_global_bits(
                    s,
                    z,
                )
            )

            checked += 1

            if (
                predicted != actual
                or
                global_pred != actual
            ):

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"frame={s.frame}",
                        f"z={z}",
                        f"a={a}",
                        f"b={b}",
                        f"factor_pred={predicted}",
                        f"global_pred={global_pred}",
                        f"actual={actual}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================
# Exact relationship:
#
#   survival to z+1
#       iff
#   residual a,b are both even.
#
# Furthermore:
#
#   x_z parity =
#       (b-a)/2 parity     A
#       (3a-b)/2 parity    B
#
# and:
#
#   child residue bit = x_z parity.
# ==============================================================================

def test_next_transition(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 7: RESIDUAL PARITY -> NEXT TRANSITION"
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

            if not factor_alive_exact(
                s,
                z,
            ):
                continue

            residuals = factor_residuals(
                s,
                z,
            )

            assert residuals is not None

            a, b = residuals

            factor_next = factor_alive_exact(
                s,
                z + 1,
            )

            predicted_next = (
                (a & 1) == 0
                and
                (b & 1) == 0
            )

            local = reconstruct_local_xy(
                s,
                z,
            )

            assert local is not None

            x, y = local

            xbit = abs(x) & 1
            ybit = abs(y) & 1

            predicted_xbit = (
                (
                    (
                        3 * a - b
                    )
                    if s.frame == "B"
                    else (
                        b - a
                    )
                )
                // 2
            ) & 1

            predicted_ybit = (
                (
                    (
                        3 * a + b
                    )
                    if s.frame == "B"
                    else (
                        a + b
                    )
                )
                // 2
            ) & 1

            checked += 1

            good = (
                factor_next
                ==
                predicted_next
                and
                xbit
                ==
                predicted_xbit
                and
                ybit
                ==
                predicted_ybit
            )

            if not good:

                failures += 1

                if shown < EXAMPLE_LIMIT:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                        f"frame={s.frame}",
                        f"a={a}",
                        f"b={b}",
                        f"next={factor_next}",
                        f"pred_next={predicted_next}",
                        f"xbit={xbit}/{predicted_xbit}",
                        f"ybit={ybit}/{predicted_ybit}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================
# Exact valuation-depth structure.
# ==============================================================================

def report_depth_classes(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 8: V2 DEPTH CLASSIFICATION"
    )
    print("=" * 90)

    classes: defaultdict[
        tuple,
        int,
    ] = defaultdict(int)

    for s in states:

        if s.frame == "A":

            alpha = v2(
                s.p - 3
            )

            beta = v2(
                s.q + 3
            )

        else:

            alpha = v2(
                s.p + 1
            )

            beta = v2(
                s.q - 3
            )

        depth = exact_factor_depth(
            s
        )

        relation = (
            "equal"
            if alpha == beta
            else "unequal"
        )

        classes[
            (
                s.frame,
                relation,
                min(alpha, beta),
                depth,
            )
        ] += 1

    print(
        "    frame / valuation relation / "
        "min(v2) / depth / count"
    )

    for key, count in sorted(
        classes.items()
    ):

        print(
            f"    {key} -> {count}"
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

    for s in states:

        if shown >= EXAMPLE_LIMIT:
            break

        xy_depth = (
            min(
                v2(s.X),
                v2(s.Y),
            )
            + 1
        )

        if s.frame == "A":

            alpha = v2(
                s.p - 3
            )

            beta = v2(
                s.q + 3
            )

        else:

            alpha = v2(
                s.p + 1
            )

            beta = v2(
                s.q - 3
            )

        f_depth = exact_factor_depth(
            s
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
            f"    v2 factor residuals="
            f"({alpha},{beta})"
        )

        print(
            f"    xy_depth={xy_depth}"
        )

        print(
            f"    factor_depth={f_depth}"
        )

        for z in range(
            2,
            min(
                f_depth + 2,
                MAX_Z + 1,
            ),
        ):

            residuals = factor_residuals(
                s,
                z,
            )

            global_alive = xy_alive(
                s.X,
                s.Y,
                z,
            )

            factor_alive = (
                factor_alive_exact(
                    s,
                    z,
                )
            )

            print(
                f"    z={z:<2} "
                f"residuals={residuals} "
                f"global={global_alive} "
                f"factor={factor_alive}"
            )

        shown += 1


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
Let:

    M = 2^(z-1).

FRAME A:

    p = 3 + M a
    q = -3 + M b.

Then:

    X
      = (q-p+6)/2
      = M(b-a)/2

    Y
      = (p+q)/2
      = M(a+b)/2.

Hence x_z=X/M and y_z=Y/M are:

    x_z = (b-a)/2
    y_z = (a+b)/2.

Therefore:

    A is represented at level z
    iff
    a and b have the same parity.

FRAME B:

    p = -1 + M a
    q =  3 + M b.

Then:

    X = M(3a-b)/2
    Y = M(3a+b)/2.

Therefore:

    x_z = (3a-b)/2
    y_z = (3a+b)/2.

Again:

    B is represented at level z
    iff
    a and b have the same parity.

So the exact factor-side survival automaton is:

    residual parity:

        (a,b)=(0,0) -> LIVE
        (a,b)=(1,1) -> LIVE

        (a,b)=(0,1) -> DEAD
        (a,b)=(1,0) -> DEAD.


THE NEXT LEVEL

At level z+1 the new residuals are:

    a' = a/2
    b' = b/2

whenever a,b are even.

The local x-bit is:

    A:
        x_z mod 2
        =
        ((b-a)/2) mod 2.

    B:
        x_z mod 2
        =
        ((3a-b)/2) mod 2.

Since 3 == 1 mod 2, the B expression has the same
parity structure as (a-b)/2.

Thus the child residue bit is determined by the
second binary layer of the factor residuals.


IMPORTANT DEPTH RESULT

Let:

    alpha = v2(A)
    beta  = v2(B).

If alpha == beta == t:

    level t+1 is represented,
    level t+2 is not.

Therefore:

    depth=t+1.

If alpha != beta:

    at level min(alpha,beta)+1,
    one normalized residual is odd and the other even.

Therefore that level is already dead:

    depth=min(alpha,beta).

This explains all apparent one-level discrepancies from the
previous factor-side valuation formula.


CORE STRUCTURE

    (FRAME, p, q)
          |
          v
    factor residual stream
          |
          +---- parity(a,b)
          |        |
          |        +---- LIVE/DEAD
          |
          +---- ((b-a)/2) or ((3a-b)/2)
          |        |
          |        +---- X-bit
          |
          +---- next residual pair

This is now the factor-side analogue of the global
(X,Y) bit-stream automaton.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print(
        "EXPERIMENT 610 START"
    )
    print("=" * 90)

    print()
    print("[1] PRIME SIEVE")

    flags = sieve(
        LIMIT
    )

    print()
    print("[2] SEMIPRIME GENERATION")

    states = build_states(
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

    test_factor_parity_bridge(
        states
    )

    test_exact_survival_equivalence(
        states
    )

    test_local_coordinates(
        states
    )

    test_depth_formula(
        states
    )

    test_residual_bit_crosswalk(
        states
    )

    test_child_residue(
        states
    )

    test_next_transition(
        states
    )

    report_depth_classes(
        states
    )

    symbolic_summary()

    show_examples(
        states
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 610 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
