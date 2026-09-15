#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 602
# ==============================================================================
# GLOBAL INVARIANT / 2-ADIC DESCENT COLLAPSE
#
# Hypothesis
# ---------
#
# At level z:
#
#     k_z = 2^(z-1)
#
#     x_z = X / k_z
#     y_z = Y / k_z
#
# where X,Y are level-independent invariant coordinates.
#
# Experiment 601 established empirically:
#
#     child exists
#         iff
#     x_z even AND y_z even
#
# and when it survives:
#
#     x_(z+1) = x_z / 2
#     y_(z+1) = y_z / 2
#
#     frame_(z+1) = frame_z
#
# This experiment asks whether the COMPLETE hierarchy can therefore
# be described directly from the initial X,Y.
#
# Main tests
# ----------
#
# 1. Recover X,Y from z=2.
# 2. Predict every x_z,y_z directly from X,Y.
# 3. Predict survival directly from divisibility:
#
#        2^z | X
#        2^z | Y
#
# 4. Predict deepest level:
#
#        deepest = 1 + min(v2(X), v2(Y))
#
#    with the measured hierarchy starting at z=2.
#
# 5. Predict every active residue directly from X,Y parity/divisibility.
# 6. Check whether frame is completely determined by the initial
#    n mod 4 residue.
# 7. Derive minimal Boolean state:
#
#        frame
#        X bit z-1
#        Y bit z-1
#
# 8. Verify the CORRECT CNF encoding of:
#
#        S <-> (!X & !Y)
#
#    and compare it against the observed data.
#
# No external files.
# No external sources.
# No output files.
# ==============================================================================

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

SHOW_EXAMPLES = 20


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:

    n: int
    p: int
    q: int

    z: int
    residue: int

    frame: str

    x: int
    y: int

    X: int
    Y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
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
            start:
            limit + 1:
            p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_semiprimes(
    limit,
    sieve,
):

    primes = [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > limit:
            break

        max_q = limit // p

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
# FRAME A
# ==============================================================================

def frame_A(
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    denominator = 2 * k

    x_num = (
        q - p + 6
    )

    y_num = (
        p + q
    )

    if x_num % denominator:
        return None

    if y_num % denominator:
        return None

    return (
        x_num // denominator,
        y_num // denominator,
    )


# ==============================================================================
# FRAME B
# ==============================================================================

def frame_B(
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    denominator = 2 * k

    x_num = (
        3 * p
        - q
        + 6
    )

    y_num = (
        3 * p
        + q
    )

    if x_num % denominator:
        return None

    if y_num % denominator:
        return None

    return (
        x_num // denominator,
        y_num // denominator,
    )


# ==============================================================================
# INVARIANT COORDINATES
# ==============================================================================

def invariant_from_state(
    state,
):

    k = 1 << (
        state.z - 1
    )

    return (
        k * state.x,
        k * state.y,
    )


# ==============================================================================
# LEVEL BUILD
# ==============================================================================

def build_levels(
    semiprimes,
):

    levels = {
        z: {}
        for z in range(
            MIN_Z,
            MAX_Z + 1,
        )
    }

    for n, p, q in semiprimes:

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            a = frame_A(
                p,
                q,
                z,
            )

            if a is not None:

                k = 1 << (
                    z - 1
                )

                X = k * a[0]
                Y = k * a[1]

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=n % (
                        1 << z
                    ),
                    frame="A",
                    x=a[0],
                    y=a[1],
                    X=X,
                    Y=Y,
                )

                continue

            b = frame_B(
                p,
                q,
                z,
            )

            if b is not None:

                k = 1 << (
                    z - 1
                )

                X = k * b[0]
                Y = k * b[1]

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=n % (
                        1 << z
                    ),
                    frame="B",
                    x=b[0],
                    y=b[1],
                    X=X,
                    Y=Y,
                )

    return levels


# ==============================================================================
# v2
# ==============================================================================

def v2(value):

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
# CHILD
# ==============================================================================

def child_state(
    levels,
    z,
    n,
):

    return levels[
        z + 1
    ].get(n)


# ==============================================================================
# FRAME FROM INITIAL MOD 4
# ==============================================================================

def frame_from_initial_residue(
    residue,
):

    if residue == 3:
        return "A"

    if residue == 1:
        return "B"

    return None


# ==============================================================================
# DIRECT LEVEL COORDINATES
# ==============================================================================

def predicted_xy_from_XY(
    X,
    Y,
    z,
):

    k = 1 << (
        z - 1
    )

    if X % k:
        return None

    if Y % k:
        return None

    return (
        X // k,
        Y // k,
    )


# ==============================================================================
# DIRECT SURVIVAL PREDICTION
# ==============================================================================

def predicted_survival(
    X,
    Y,
    z,
):

    required = 1 << z

    return (
        X % required == 0
        and
        Y % required == 0
    )


# ==============================================================================
# DIRECT FRAME PREDICTION
# ==============================================================================

def predicted_frame(
    n,
):

    return frame_from_initial_residue(
        n % 4
    )


# ==============================================================================
# DIRECT RESIDUE PREDICTION
# ==============================================================================

def predicted_residue(
    X,
    Y,
    frame,
    z,
):

    modulus = 1 << z

    k = 1 << (
        z - 1
    )

    x = X // k

    if frame == "A":

        value = (
            (
                k * (Y // k - X // k)
                + 3
            )
            *
            (
                k * (Y // k + X // k)
                - 3
            )
        )

        return value % modulus

    if frame == "B":

        numerator = (
            (
                k * (X // k + Y // k)
                - 3
            )
            *
            (
                k * (Y // k - X // k)
                + 3
            )
        )

        if numerator % 3:
            return None

        return (
            numerator // 3
        ) % modulus

    return None


# ==============================================================================
# [1] INVARIANT CONSISTENCY
# ==============================================================================

def test_invariant_consistency(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: GLOBAL X,Y INVARIANCE"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    by_n = {}

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        for n, state in levels[z].items():

            current = (
                state.X,
                state.Y,
            )

            if n in by_n:

                if current != by_n[n]:

                    failures += 1

                    if failures <= SHOW_EXAMPLES:

                        print(
                            f"    mismatch "
                            f"n={n} "
                            f"z={z} "
                            f"previous={by_n[n]} "
                            f"current={current}"
                        )

            else:

                by_n[n] = current

            checked += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# [2] DIRECT x,y RECONSTRUCTION
# ==============================================================================

def test_direct_coordinates(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: DIRECT X,Y -> x_z,y_z"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        failures = 0
        checked = 0

        for state in levels[z].values():

            predicted = (
                predicted_xy_from_XY(
                    state.X,
                    state.Y,
                    z,
                )
            )

            checked += 1

            if predicted != (
                state.x,
                state.y,
            ):

                failures += 1

        print(
            f"z={z} "
            f"checked={checked} "
            f"failures={failures}"
        )


# ==============================================================================
# [3] DIRECT SURVIVAL
# ==============================================================================

def test_direct_survival(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: DIRECT X,Y SURVIVAL"
    )
    print("=" * 90)

    global_failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        checked = 0
        failures = 0

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            actual = (
                child is not None
            )

            predicted = (
                predicted_survival(
                    state.X,
                    state.Y,
                    z,
                )
            )

            checked += 1

            if predicted != actual:

                failures += 1

                if failures <= SHOW_EXAMPLES:

                    print(
                        f"    "
                        f"z={z} "
                        f"n={state.n} "
                        f"X={state.X} "
                        f"Y={state.Y} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

        global_failures += failures

        print(
            f"z={z} -> z={z+1} "
            f"checked={checked} "
            f"failures={failures}"
        )

    print()
    print(
        f"GLOBAL failures="
        f"{global_failures}"
    )


# ==============================================================================
# [4] DIRECT FRAME
# ==============================================================================

def test_direct_frame(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: INITIAL RESIDUE -> FRAME"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        for state in levels[z].values():

            predicted = (
                predicted_frame(
                    state.n
                )
            )

            checked += 1

            if predicted != state.frame:

                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# [5] DEEPEST LEVEL
# ==============================================================================

def test_deepest_level(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: EXACT DEEPEST LEVEL FORMULA"
    )
    print("=" * 90)

    # Every n exists at z=2.

    failures = 0
    checked = 0

    depth_distribution = Counter()

    examples = []

    for n, state2 in levels[
        MIN_Z
    ].items():

        vx = v2(
            state2.X
        )

        vy = v2(
            state2.Y
        )

        predicted = (
            min(vx, vy)
            + 1
        )

        actual = MIN_Z - 1

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            if n in levels[z]:

                actual = z

        # Limit is only up to MAX_Z,
        # so do not claim a mismatch for
        # predicted depths beyond the
        # observed range.

        if predicted <= MAX_Z:

            checked += 1

            if predicted != actual:

                failures += 1

                if len(examples) < SHOW_EXAMPLES:

                    examples.append(
                        (
                            n,
                            state2.p,
                            state2.q,
                            state2.X,
                            state2.Y,
                            vx,
                            vy,
                            predicted,
                            actual,
                        )
                    )

            depth_distribution[
                actual
            ] += 1

    print(
        f"bounded checked={checked} "
        f"failures={failures}"
    )

    if examples:

        print()
        print(
            "COUNTEREXAMPLES"
        )

        for item in examples:

            print(
                "    "
                f"n={item[0]} "
                f"p={item[1]} "
                f"q={item[2]} "
                f"X={item[3]} "
                f"Y={item[4]} "
                f"v2X={item[5]} "
                f"v2Y={item[6]} "
                f"predicted={item[7]} "
                f"actual={item[8]}"
            )

    print()
    print(
        "OBSERVED DEPTH DISTRIBUTION"
    )

    for depth, count in sorted(
        depth_distribution.items()
    ):

        print(
            f"    depth={depth}: "
            f"{count}"
        )


# ==============================================================================
# [6] DIRECT RESIDUE FROM X,Y
# ==============================================================================

def test_residue_prediction(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: DIRECT RESIDUE FROM X,Y"
    )
    print("=" * 90)

    global_failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        failures = 0
        checked = 0

        for state in levels[z].values():

            predicted = predicted_residue(
                state.X,
                state.Y,
                state.frame,
                z,
            )

            checked += 1

            if predicted != state.residue:

                failures += 1

                if failures <= SHOW_EXAMPLES:

                    print(
                        f"    "
                        f"z={z} "
                        f"n={state.n} "
                        f"frame={state.frame} "
                        f"pred={predicted} "
                        f"actual={state.residue}"
                    )

        global_failures += failures

        print(
            f"z={z} "
            f"checked={checked} "
            f"failures={failures}"
        )

    print()
    print(
        f"GLOBAL failures="
        f"{global_failures}"
    )


# ==============================================================================
# [7] BOOLEAN MINIMALITY
# ==============================================================================

def test_boolean_minimality(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: MINIMAL BOOLEAN SURVIVAL STATE"
    )
    print("=" * 90)

    # --------------------------------------------------------------
    # frame + xbit + ybit
    # --------------------------------------------------------------

    groups_full = defaultdict(
        set
    )

    # frame + xbit
    groups_no_y = defaultdict(
        set
    )

    # frame + ybit
    groups_no_x = defaultdict(
        set
    )

    # xbit + ybit
    groups_no_frame = defaultdict(
        set
    )

    # xbit only
    groups_x = defaultdict(
        set
    )

    # ybit only
    groups_y = defaultdict(
        set
    )

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            outcome = (
                child is not None
            )

            xb = state.x & 1
            yb = state.y & 1

            groups_full[
                (
                    state.frame,
                    xb,
                    yb,
                )
            ].add(
                outcome
            )

            groups_no_y[
                (
                    state.frame,
                    xb,
                )
            ].add(
                outcome
            )

            groups_no_x[
                (
                    state.frame,
                    yb,
                )
            ].add(
                outcome
            )

            groups_no_frame[
                (
                    xb,
                    yb,
                )
            ].add(
                outcome
            )

            groups_x[
                xb
            ].add(
                outcome
            )

            groups_y[
                yb
            ].add(
                outcome
            )

    tests = [
        (
            "frame+xbit+ybit",
            groups_full,
        ),
        (
            "frame+xbit",
            groups_no_y,
        ),
        (
            "frame+ybit",
            groups_no_x,
        ),
        (
            "xbit+ybit",
            groups_no_frame,
        ),
        (
            "xbit only",
            groups_x,
        ),
        (
            "ybit only",
            groups_y,
        ),
    ]

    for name, groups in tests:

        deterministic = all(
            len(values) == 1
            for values in groups.values()
        )

        ambiguous = sum(
            len(values) > 1
            for values in groups.values()
        )

        print(
            f"    {name:<20} "
            f"deterministic="
            f"{deterministic} "
            f"ambiguous_states="
            f"{ambiguous}"
        )


# ==============================================================================
# [8] CORRECT CNF TEST
# ==============================================================================

def test_correct_cnf(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: CORRECT CNF FOR SURVIVAL"
    )
    print("=" * 90)

    print(
        """
We encode:

    S <-> (!X & !Y)

Correct CNF:

    (!S v !X)
    (!S v !Y)
    ( S v  X v  Y)

where:

    X = x mod 2
    Y = y mod 2
    S = child survives.
"""
    )

    failures = 0
    checked = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            X = state.x & 1
            Y = state.y & 1
            S = 1 if (
                child is not None
            ) else 0

            clauses = [

                (
                    S == 0
                    or X == 0
                ),

                (
                    S == 0
                    or Y == 0
                ),

                (
                    S == 1
                    or X == 1
                    or Y == 1
                ),
            ]

            checked += 1

            if not all(clauses):

                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# [9] BITWISE INTERPRETATION
# ==============================================================================

def bitwise_interpretation(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 9: BITWISE DESCENT INTERPRETATION"
    )
    print("=" * 90)

    print(
        """
For a fixed invariant pair X,Y:

    x_z = X >> (z-1)
    y_z = Y >> (z-1)

provided the corresponding representation exists.

The next descent requires:

    bit_(z-1)(X) = 0
    bit_(z-1)(Y) = 0.

Therefore the hierarchy should be equivalent to scanning the
binary expansions of X and Y simultaneously.

The first nonzero bit in either coordinate terminates descent.
"""
    )

    failures = 0

    for n, state in levels[
        MIN_Z
    ].items():

        X = state.X
        Y = state.Y

        for z in range(
            MIN_Z,
            MAX_Z,
        ):

            child = child_state(
                levels,
                z,
                n,
            )

            required_bit = z

            xbit = (
                (abs(X) >> required_bit)
                & 1
            )

            ybit = (
                (abs(Y) >> required_bit)
                & 1
            )

            predicted = (
                xbit == 0
                and
                ybit == 0
            )

            actual = (
                child is not None
            )

            if predicted != actual:

                failures += 1

                if failures <= SHOW_EXAMPLES:

                    print(
                        f"    mismatch "
                        f"n={n} "
                        f"z={z} "
                        f"X={X} "
                        f"Y={Y} "
                        f"xbit={xbit} "
                        f"ybit={ybit} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

    print(
        f"GLOBAL failures={failures}"
    )


# ==============================================================================
# [10] COMPACT RULE TABLE
# ==============================================================================

def print_rule_table():

    print()
    print("=" * 90)
    print(
        "TEST 10: MINIMAL TRANSITION TABLE"
    )
    print("=" * 90)

    print(
        """
FRAME A:

    x y | result
    ----+--------
    0 0 | A
    0 1 | DEAD
    1 0 | DEAD
    1 1 | DEAD

FRAME B:

    x y | result
    ----+--------
    0 0 | B
    0 1 | DEAD
    1 0 | DEAD
    1 1 | DEAD

Residue of the surviving child:

    A:
        r' = -9 + 2^z*x
             mod 2^(z+1)

    B:
        r' = -3 + 2^z*x
             mod 2^(z+1)
"""
    )


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary():

    print()
    print("=" * 90)
    print(
        "FINAL AUDIT"
    )
    print("=" * 90)

    print(
        r"""
The strongest possible structure is now:

    INITIAL STATE
        =
    (frame_2, X, Y)

where:

    frame_2 = A  iff n mod 4 = 3
    frame_2 = B  iff n mod 4 = 1.

At level z:

    x_z = X / 2^(z-1)
    y_z = Y / 2^(z-1)

and descent occurs iff:

    2^z | X
    2^z | Y.

Equivalently:

    bit_z(X) = 0
    bit_z(Y) = 0.

When descent occurs:

    X,Y remain unchanged,
    frame remains unchanged,
    x and y are halved.

The deepest reachable level is therefore:

    deepest =
        1 + min(v2(X), v2(Y))

subject to the measured level range.

This would replace the whole level-by-level construction by a
single 2-adic divisibility process on two invariant integers.

The remaining nontrivial component is then the factor/frame
mapping:

    FRAME A:
        p = Y - X + 3
        q = Y + X - 3

    FRAME B:
        3p = Y + X - 3
        q  = Y - X + 3

with X,Y interpreted as the invariant coordinates at the chosen
normalization.

If every test above is exact, the hierarchical system has collapsed
to:

    (frame, X, Y)
          |
          +-- inspect binary bits of X,Y
          |
          +-- divide X,Y by 2 at every surviving level
          |
          +-- retain the same frame.

That is substantially smaller than maintaining independent
A_z/B_z functions at every level.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 602 START"
    )
    print("=" * 90)

    print()
    print(
        "GLOBAL INVARIANT / 2-ADIC DESCENT COLLAPSE"
    )

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    print()
    print(
        "[2] SEMIPRIME GENERATION"
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes="
        f"{len(semiprimes)}"
    )

    print()
    print(
        "[3] BUILD LEVELS"
    )

    levels = build_levels(
        semiprimes
    )

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print(
            f"    z={z} "
            f"states={len(levels[z])}"
        )

    test_invariant_consistency(
        levels
    )

    test_direct_coordinates(
        levels
    )

    test_direct_survival(
        levels
    )

    test_direct_frame(
        levels
    )

    test_deepest_level(
        levels
    )

    test_residue_prediction(
        levels
    )

    test_boolean_minimality(
        levels
    )

    test_correct_cnf(
        levels
    )

    bitwise_interpretation(
        levels
    )

    print_rule_table()

    final_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 602 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
