#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 611
# ==============================================================================
#
# EXACT 2-ADIC RESIDUAL CONGRUENCE LADDER
#
# Experiment 610 established:
#
#   level-z representation
#       <=> factor residuals a,b have equal parity.
#
# But Test 7 showed that:
#
#   a,b even
#
# is NOT enough for survival to z+1.
#
# The exact next-level condition is expected to be:
#
#   a ≡ b ≡ 0 (mod 4)
#
# for the residuals at level z.
#
# More generally, this experiment tests the complete ladder:
#
#   live at level z:
#       a ≡ b (mod 2)
#
#   live at z+1:
#       a ≡ b ≡ 0 (mod 4)
#
#   live at z+2:
#       a ≡ b ≡ 0 (mod 8)
#
# etc.
#
# We will derive:
#
#   A:
#       a = (p-3) / 2^(z-1)
#       b = (q+3) / 2^(z-1)
#
#   B:
#       a = (p+1) / 2^(z-1)
#       b = (q-3) / 2^(z-1)
#
# and test whether the complete hierarchy is equivalent to:
#
#       a ≡ b ≡ 0 mod 2^m
#
# for the appropriate number of future descents.
#
# ==============================================================================


LIMIT = 2_000_000
MAX_Z = 24
EXAMPLES = 20


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

    for p in range(
        2,
        root + 1,
    ):

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

            if n % 4 == 3:

                frame = "A"

                X = (
                    q - p + 6
                ) // 2

                Y = (
                    p + q
                ) // 2

            else:

                frame = "B"

                X = (
                    3 * p
                    - q
                    + 6
                ) // 2

                Y = (
                    3 * p
                    + q
                ) // 2

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
# V2
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return 10**9

    x = abs(x)

    result = 0

    while (x & 1) == 0:

        x >>= 1
        result += 1

    return result


# ==============================================================================
# GLOBAL DEPTH
# ==============================================================================

def xy_depth(
    s: State,
) -> int:

    return (
        min(
            v2(s.X),
            v2(s.Y),
        )
        + 1
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
# EXACT LEVEL REPRESENTATION
# ==============================================================================

def factor_live(
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
        (a - b) % 2 == 0
    )


# ==============================================================================
# EXPECTED RESIDUAL CONDITION
# ==============================================================================
#
# At level z:
#
#   LIVE
#       <=> a == b (mod 2)
#
# For survival through k additional transitions:
#
#   a and b must have enough common 2-adic divisibility.
#
# The conjectured exact condition is:
#
#   survive k transitions
#       <=>
#   a ≡ b ≡ 0 (mod 2^k)
#
# once z is itself a valid level.
#
# This is tested below.
# ==============================================================================

def residual_survival_prediction(
    a: int,
    b: int,
    transitions: int,
) -> bool:

    if transitions == 0:

        return (
            (a - b) % 2 == 0
        )

    modulus = 1 << transitions

    return (
        a % modulus == 0
        and
        b % modulus == 0
    )


# ==============================================================================
# ACTUAL SURVIVAL k LEVELS FORWARD
# ==============================================================================

def actual_survives_k(
    s: State,
    z: int,
    k: int,
) -> bool:

    for j in range(
        k + 1
    ):

        level = z + j

        if not factor_live(
            s,
            level,
        ):

            return False

    return True


# ==============================================================================
# TEST 1
# ==============================================================================

def test_one_step_condition(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 1: EXACT ONE-STEP RESIDUAL CONDITION"
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

            residuals = factor_residuals(
                s,
                z,
            )

            if residuals is None:
                continue

            a, b = residuals

            actual = factor_live(
                s,
                z + 1,
            )

            # Derived condition:
            #
            # x_z even and y_z even.
            #
            # A:
            #     x=(b-a)/2
            #     y=(a+b)/2
            #
            # B:
            #     x=(3a-b)/2
            #     y=(3a+b)/2
            #
            # In both frames the exact condition is:
            #
            #     a ≡ b ≡ 0 mod 4
            #
            predicted = (
                a % 4 == 0
                and
                b % 4 == 0
            )

            checked += 1

            if actual != predicted:

                failures += 1

                if shown < EXAMPLES:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                        f"frame={s.frame}",
                        f"a={a}",
                        f"b={b}",
                        f"actual={actual}",
                        f"predicted={predicted}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================
# k-step residual congruence ladder.
#
# At a live level z, test:
#
#     actual_survives_k(z,k)
#
# against:
#
#     a,b divisible by 2^(k+1)?
#
# Careful:
#
# k=0 means current level only:
#
#     a == b mod 2
#
# k=1 means current + child:
#
#     a,b divisible by 4.
#
# k=2 means three consecutive live levels:
#
#     a,b divisible by 8.
#
# So:
#
#     k future transitions
#         ->
#     divisibility by 2^(k+1).
# ==============================================================================

def test_congruence_ladder(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 2: MULTI-LEVEL 2-ADIC CONGRUENCE LADDER"
    )
    print("=" * 90)

    for k in range(
        0,
        6,
    ):

        checked = 0
        failures = 0
        shown = 0

        for s in states:

            for z in range(
                2,
                MAX_Z - k,
            ):

                residuals = factor_residuals(
                    s,
                    z,
                )

                if residuals is None:
                    continue

                a, b = residuals

                actual = actual_survives_k(
                    s,
                    z,
                    k,
                )

                modulus = 1 << (
                    k + 1
                )

                predicted = (
                    a % modulus == 0
                    and
                    b % modulus == 0
                )

                checked += 1

                if actual != predicted:

                    failures += 1

                    if shown < 3:

                        print(
                            "    mismatch",
                            f"k={k}",
                            f"n={s.n}",
                            f"z={z}",
                            f"a={a}",
                            f"b={b}",
                            f"actual={actual}",
                            f"predicted={predicted}",
                        )

                        shown += 1

        print(
            f"k={k} "
            f"checked={checked} "
            f"failures={failures}"
        )


# ==============================================================================
# TEST 3
# ==============================================================================
# Connect the residual ladder directly to X,Y valuation.
#
# From Experiment 610:
#
# A:
#
#     X = M(b-a)/2
#     Y = M(a+b)/2
#
# B:
#
#     X = M(3a-b)/2
#     Y = M(3a+b)/2.
#
# We investigate whether:
#
#     min(v2(X),v2(Y))
#
# is exactly determined by:
#
#     z-2 + min(v2(a-b), v2(a+b))
#
# which must hold for both frames.
#
# ==============================================================================

def test_valuation_bridge(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 3: V2(X,Y) <-> FACTOR RESIDUAL V2"
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

            if residuals is None:
                continue

            a, b = residuals

            if s.frame == "A":

                alpha = v2(
                    b - a
                )

                beta = v2(
                    a + b
                )

            else:

                alpha = v2(
                    3 * a - b
                )

                beta = v2(
                    3 * a + b
                )

            local_min = min(
                alpha,
                beta,
            )

            predicted = (
                (z - 2)
                + local_min
            )

            actual = min(
                v2(s.X),
                v2(s.Y),
            )

            checked += 1

            if predicted != actual:

                failures += 1

                if shown < EXAMPLES:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                        f"frame={s.frame}",
                        f"a={a}",
                        f"b={b}",
                        f"predicted={predicted}",
                        f"actual={actual}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================
# Residual-depth formula.
#
# Instead of comparing valuations of the shifted factors directly,
# calculate:
#
#     d = z + number of successful future transitions.
#
# The conjecture:
#
#     if a,b are both divisible by 2^t but not 2^(t+1),
#
#       deepest = z + t - 1
#
# provided the current level z is live.
#
# More precisely, for residuals at level z:
#
#     t = min(v2(a),v2(b))
#
#     if t is finite:
#
#         deepest = z + t - 1
#
# This should also handle the equal-valuation cases.
#
# ==============================================================================

def residual_depth_prediction(
    s: State,
    z: int,
) -> int | None:

    residuals = factor_residuals(
        s,
        z,
    )

    if residuals is None:
        return None

    a, b = residuals

    if (a - b) % 2 != 0:
        return None

    t = min(
        v2(a),
        v2(b),
    )

    if t >= 10**9:
        return MAX_Z

    return (
        z + t - 1
    )


def test_residual_depth(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 4: RESIDUAL DEPTH FROM CURRENT LEVEL"
    )
    print("=" * 90)

    checked = 0
    failures = 0
    shown = 0

    for s in states:

        true_depth = xy_depth(
            s
        )

        for z in range(
            2,
            min(
                true_depth + 1,
                MAX_Z + 1,
            ),
        ):

            prediction = residual_depth_prediction(
                s,
                z,
            )

            if prediction is None:
                continue

            checked += 1

            # Actual maximum represented level
            # from this point is simply true_depth.
            actual = true_depth

            if prediction != actual:

                failures += 1

                if shown < EXAMPLES:

                    print(
                        "    mismatch",
                        f"n={s.n}",
                        f"z={z}",
                        f"prediction={prediction}",
                        f"actual={actual}",
                    )

                    shown += 1

    print(
        f"checked={checked} failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================
# Residual automaton modulo 2^m.
#
# Determine the smallest modulus needed to decide survival
# for m future transitions.
# ==============================================================================

def test_minimal_residual_modulus(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: MINIMAL RESIDUAL MODULUS"
    )
    print("=" * 90)

    for future in range(
        0,
        6,
    ):

        modulus = 1 << (
            future + 1
        )

        ambiguous = set()

        mapping: defaultdict[
            tuple[int, int],
            set[bool],
        ] = defaultdict(set)

        for s in states:

            for z in range(
                2,
                MAX_Z - future,
            ):

                residuals = factor_residuals(
                    s,
                    z,
                )

                if residuals is None:
                    continue

                a, b = residuals

                key = (
                    a % modulus,
                    b % modulus,
                )

                value = actual_survives_k(
                    s,
                    z,
                    future,
                )

                mapping[
                    key
                ].add(
                    value
                )

        for key, values in mapping.items():

            if len(values) > 1:

                ambiguous.add(
                    key
                )

        print(
            f"future={future} "
            f"modulus={modulus} "
            f"states={len(mapping)} "
            f"ambiguous={len(ambiguous)}"
        )


# ==============================================================================
# TEST 6
# ==============================================================================
# Show the actual residual recursion.
#
# If:
#
#     a,b divisible by 2,
#
# then:
#
#     a' = a/2
#     b' = b/2
#
# But a future representation requires:
#
#     a' == b' (mod 2)
#
# i.e.
#
#     a == b (mod 4).
#
# At another step:
#
#     a'' = a/4
#     b'' = b/4
#
# and so forth.
#
# ==============================================================================

def show_recursion(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 6: RESIDUAL RECURSION EXAMPLES"
    )
    print("=" * 90)

    count = 0

    for s in states:

        if count >= EXAMPLES:
            break

        depth = xy_depth(
            s
        )

        residuals = factor_residuals(
            s,
            2,
        )

        if residuals is None:
            continue

        a, b = residuals

        print()
        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    depth={depth}"
        )

        print(
            f"    level 2 residuals="
            f"({a},{b})"
        )

        for j in range(
            0,
            min(
                6,
                MAX_Z - 2,
            ),
        ):

            modulus = 1 << (
                j + 1
            )

            print(
                f"    future={j} "
                f"need a=b=0 mod {modulus}: "
                f"{a % modulus == 0 and b % modulus == 0}"
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
THE NEW RESULT TO TEST

Experiment 610 established:

    current representation:

        a == b (mod 2).

But Experiment 610 Test 7 showed that:

    a,b even

does not imply another surviving level.

For example:

    n=9
    frame=B
    level z=2:

        a=2
        b=0

    both are even,

    but:

        x_2 = (3a-b)/2
            = 3

        y_2 = (3a+b)/2
            = 3

    so the next level requires both x_2 and y_2
    to be even, which they are not.

Therefore the exact one-step condition is:

    a ≡ b ≡ 0 (mod 4).


GENERALIZATION

At level z let:

    a_z
    b_z

be the normalized factor residuals.

Then:

    representation at z:
        a_z ≡ b_z (mod 2)

    representation at z+1:
        a_z ≡ b_z ≡ 0 (mod 4)

    representation at z+2:
        a_z ≡ b_z ≡ 0 (mod 8)

and generally:

    survive k additional levels
        iff

    a_z ≡ b_z ≡ 0
        (mod 2^(k+1)).

Equivalently the residuals are consumed by a common
2-adic shift:

    (a,b)
        ->
    (a/2,b/2)
        ->
    (a/4,b/4)
        ->
    ...

and every step requires the new pair to have equal parity.

Thus the factor-side automaton is itself a binary
divisibility stream.


GLOBAL/FECTOR BRIDGE

The global variables satisfy:

    A-frame:
        X = 2^(z-2)(b-a)
        Y = 2^(z-2)(a+b)

    B-frame:
        X = 2^(z-2)(3a-b)
        Y = 2^(z-2)(3a+b).

Consequently each additional successful level consumes
one more common binary digit of the factor residual pair.


EXPECTED DEPTH

At level z, define:

    t = min(v2(a), v2(b)).

For a currently represented state:

    deepest = z + t - 1

should follow.

The experiment checks this against the exact global depth.


IMPORTANT DISTINCTION

There are now three related but distinct notions:

    1. representation at current level

        a == b (mod 2)

    2. survival to next level

        a == b == 0 (mod 4)

    3. survival k levels

        a == b == 0 (mod 2^(k+1)).

This is the factor-side version of the global
binary-bit descent.


TARGET

If all tests pass, the hierarchy can be expressed
without enumerating levels:

    STATE = (FRAME, residual pair)

    consume one common 2-adic bit layer per descent.

The factor residuals then provide a direct bridge between:

    prime factors
        <-->
    invariant X,Y
        <-->
    level hierarchy
        <-->
    SAT constraints.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print(
        "EXPERIMENT 611 START"
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

    test_one_step_condition(
        states
    )

    test_congruence_ladder(
        states
    )

    test_valuation_bridge(
        states
    )

    test_residual_depth(
        states
    )

    test_minimal_residual_modulus(
        states
    )

    show_recursion(
        states
    )

    symbolic_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 611 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
