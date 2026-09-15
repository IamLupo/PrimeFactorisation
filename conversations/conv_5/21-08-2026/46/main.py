#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 605
# ==============================================================================
# COMPLETE GLOBAL COLLAPSE:
#
#     (p,q)
#       |
#       v
#   FRAME, X, Y
#       |
#       +--> p,q
#       +--> n
#       +--> n mod 2^z
#       +--> survival depth
#       +--> complete residue chain
#
# No external files.
# No external sources.
#
# Main hypotheses:
#
#   FRAME A iff n mod 4 = 3
#   FRAME B iff n mod 4 = 1
#
# FRAME A:
#
#   X = (q-p+6)/2
#   Y = (p+q)/2
#
#   p = Y-X+3
#   q = Y+X-3
#
# FRAME B:
#
#   X = (3p-q+6)/2
#   Y = (3p+q)/2
#
#   p = (Y+X-3)/3
#   q = Y-X+3
#
# The level representation is:
#
#   x_z = X / 2^(z-1)
#   y_z = Y / 2^(z-1)
#
# and z -> z+1 survives iff:
#
#   bit_(z-1)(X) = bit_(z-1)(Y) = 0.
#
# Hence:
#
#   deepest = 1 + min(v2(X), v2(Y))
#
# subject to z >= 2.
#
# ==============================================================================


from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 2_000_000
MAX_Z = 20
EXAMPLE_COUNT = 20


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class Pair:
    n: int
    p: int
    q: int


@dataclass(frozen=True)
class GlobalState:
    n: int
    p: int
    q: int
    frame: str
    X: int
    Y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def sieve_primes(limit: int):

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

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

        sieve[
            start:
            limit + 1:
            p
        ] = b"\x00" * (
            ((limit - start) // p) + 1
        )

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
                Pair(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return result


# ==============================================================================
# FRAME SELECTION
# ==============================================================================

def frame_from_n(n: int):

    r = n % 4

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Odd semiprime has impossible residue n mod 4={r}"
    )


# ==============================================================================
# DIRECT GLOBAL COORDINATES FROM FACTORS
# ==============================================================================

def global_XY_from_factors(
    p: int,
    q: int,
    frame: str,
):

    if frame == "A":

        x_num = (
            q - p + 6
        )

        y_num = (
            p + q
        )

    elif frame == "B":

        x_num = (
            3 * p
            - q
            + 6
        )

        y_num = (
            3 * p
            + q
        )

    else:

        raise ValueError(
            f"Unknown frame {frame}"
        )

    if x_num % 2:
        return None

    if y_num % 2:
        return None

    return (
        x_num // 2,
        y_num // 2,
    )


# ==============================================================================
# FACTORS FROM GLOBAL COORDINATES
# ==============================================================================

def factors_from_XY(
    frame,
    X,
    Y,
):

    if frame == "A":

        p = (
            Y - X + 3
        )

        q = (
            Y + X - 3
        )

        return (
            p,
            q,
        )

    if frame == "B":

        numerator = (
            Y + X - 3
        )

        if numerator % 3:
            return None

        p = (
            numerator // 3
        )

        q = (
            Y - X + 3
        )

        return (
            p,
            q,
        )

    return None


# ==============================================================================
# N FROM GLOBAL STATE
# ==============================================================================

def n_from_state(
    frame,
    X,
    Y,
):

    factors = factors_from_XY(
        frame,
        X,
        Y,
    )

    if factors is None:
        return None

    p, q = factors

    return p * q


# ==============================================================================
# 2-ADIC VALUATION
# ==============================================================================

def v2(
    value: int,
):

    if value == 0:
        return 10**9

    value = abs(value)

    result = 0

    while (
        value & 1
    ) == 0:

        value >>= 1
        result += 1

    return result


# ==============================================================================
# BIT
# ==============================================================================

def bit(
    value: int,
    position: int,
):

    return (
        abs(value) >> position
    ) & 1


# ==============================================================================
# LEVEL COORDINATES
# ==============================================================================

def level_xy(
    X,
    Y,
    z,
):

    divisor = 1 << (
        z - 1
    )

    if X % divisor:
        return None

    if Y % divisor:
        return None

    return (
        X // divisor,
        Y // divisor,
    )


# ==============================================================================
# LEVEL EXISTENCE
# ==============================================================================

def level_exists(
    X,
    Y,
    z,
):

    modulus = 1 << (
        z - 1
    )

    return (
        X % modulus == 0
        and
        Y % modulus == 0
    )


# ==============================================================================
# LEVEL RESIDUE
# ==============================================================================

def residue_from_state(
    frame,
    X,
    Y,
    z,
):

    n = n_from_state(
        frame,
        X,
        Y,
    )

    if n is None:
        return None

    return n % (
        1 << z
    )


# ==============================================================================
# SYMBOLIC ACTIVE RESIDUE
# ==============================================================================

def active_residue_formula(
    frame,
    z,
):

    modulus = 1 << z

    if frame == "A":
        return (-9) % modulus

    if frame == "B":
        return (-3) % modulus

    raise ValueError(
        frame
    )


# ==============================================================================
# SYMBOLIC CHILD RESIDUE
# ==============================================================================

def child_residue_formula(
    frame,
    X,
    z,
):

    x_bit = bit(
        X,
        z - 1,
    )

    modulus = 1 << (
        z + 1
    )

    if frame == "A":
        base = -9
    else:
        base = -3

    return (
        base
        + (1 << z) * x_bit
    ) % modulus


# ==============================================================================
# GLOBAL STATE BUILD
# ==============================================================================

def build_global_states(
    semiprimes,
):

    states = []

    for item in semiprimes:

        frame = frame_from_n(
            item.n
        )

        XY = global_XY_from_factors(
            item.p,
            item.q,
            frame,
        )

        if XY is None:
            raise AssertionError(
                "Could not construct X,Y"
            )

        X, Y = XY

        states.append(
            GlobalState(
                n=item.n,
                p=item.p,
                q=item.q,
                frame=frame,
                X=X,
                Y=Y,
            )
        )

    return states


# ==============================================================================
# TEST 1
# ==============================================================================

def test_frame_from_n(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: INITIAL FRAME FROM n mod 4"
    )
    print("=" * 90)

    failures = 0

    for state in states:

        expected = (
            "A"
            if state.n % 4 == 3
            else "B"
        )

        if state.frame != expected:

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_XY_forward_inverse(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: FACTORS <-> GLOBAL X,Y"
    )
    print("=" * 90)

    failures = 0

    for state in states:

        factors = factors_from_XY(
            state.frame,
            state.X,
            state.Y,
        )

        if factors != (
            state.p,
            state.q,
        ):

            failures += 1

            if failures <= EXAMPLE_COUNT:

                print(
                    f"    mismatch n={state.n} "
                    f"frame={state.frame} "
                    f"X={state.X} "
                    f"Y={state.Y} "
                    f"got={factors} "
                    f"expected={(state.p,state.q)}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_global_n(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: GLOBAL (FRAME,X,Y) -> n"
    )
    print("=" * 90)

    failures = 0

    for state in states:

        reconstructed = n_from_state(
            state.frame,
            state.X,
            state.Y,
        )

        if reconstructed != state.n:

            failures += 1

            if failures <= EXAMPLE_COUNT:

                print(
                    f"    mismatch n={state.n} "
                    f"reconstructed={reconstructed}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_global_active_residue(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: ACTIVE RESIDUE DIRECTLY FROM FRAME"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            if not level_exists(
                state.X,
                state.Y,
                z,
            ):
                continue

            residue = (
                residue_from_state(
                    state.frame,
                    state.X,
                    state.Y,
                    z,
                )
            )

            formula = active_residue_formula(
                state.frame,
                z,
            )

            checked += 1

            if residue != formula:

                failures += 1

                if failures <= EXAMPLE_COUNT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={state.frame} "
                        f"residue={residue} "
                        f"formula={formula}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_survival(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: GLOBAL SURVIVAL FROM ONE BIT PAIR"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            actual = (
                level_exists(
                    state.X,
                    state.Y,
                    z + 1,
                )
            )

            predicted = (
                bit(
                    state.X,
                    z - 1,
                ) == 0
                and
                bit(
                    state.Y,
                    z - 1,
                ) == 0
            )

            checked += 1

            if actual != predicted:

                failures += 1

                if failures <= EXAMPLE_COUNT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"X={state.X} "
                        f"Y={state.Y} "
                        f"bits=("
                        f"{bit(state.X,z-1)},"
                        f"{bit(state.Y,z-1)}"
                        f") "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_deepest(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: GLOBAL DEEPEST LEVEL"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    distribution = Counter()

    for state in states:

        predicted = (
            1
            + min(
                v2(state.X),
                v2(state.Y),
            )
        )

        actual = 1

        for z in range(
            2,
            MAX_Z + 1,
        ):

            if level_exists(
                state.X,
                state.Y,
                z,
            ):

                actual = z

        if predicted <= MAX_Z:

            checked += 1

            if predicted != actual:

                failures += 1

                if failures <= EXAMPLE_COUNT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"X={state.X} "
                        f"Y={state.Y} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

        distribution[
            actual
        ] += 1

    print(
        f"bounded checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "DEPTH DISTRIBUTION"
    )

    for depth in sorted(
        distribution
    ):

        print(
            f"    depth={depth}: "
            f"{distribution[depth]}"
        )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_child_residues(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: CHILD RESIDUE FROM SINGLE X BIT"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not level_exists(
                state.X,
                state.Y,
                z + 1,
            ):

                continue

            actual = residue_from_state(
                state.frame,
                state.X,
                state.Y,
                z + 1,
            )

            predicted = child_residue_formula(
                state.frame,
                state.X,
                z,
            )

            checked += 1

            if actual != predicted:

                failures += 1

                if failures <= EXAMPLE_COUNT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={state.frame} "
                        f"Xbit={bit(state.X,z-1)} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

    print(
        f"surviving transitions checked="
        f"{checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================

def test_frame_necessary_only_for_factor_mapping(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: FRAME ROLE SEPARATION"
    )
    print("=" * 90)

    print(
        """
SURVIVAL:

    depends only on:
        bit_(z-1)(X)
        bit_(z-1)(Y)

RESIDUE:

    depends on:
        FRAME
        bit_(z-1)(X)

FACTOR RECONSTRUCTION:

    requires:
        FRAME
        X
        Y
"""
    )

    survival_failures = 0
    residue_failures = 0

    for state in states:

        for z in range(
            2,
            min(
                MAX_Z,
                12,
            ),
        ):

            survival_direct = (
                bit(
                    state.X,
                    z - 1,
                ) == 0
                and
                bit(
                    state.Y,
                    z - 1,
                ) == 0
            )

            survival_actual = (
                level_exists(
                    state.X,
                    state.Y,
                    z + 1,
                )
            )

            if (
                survival_direct
                != survival_actual
            ):

                survival_failures += 1

            if not survival_actual:
                continue

            predicted = child_residue_formula(
                state.frame,
                state.X,
                z,
            )

            actual = residue_from_state(
                state.frame,
                state.X,
                state.Y,
                z + 1,
            )

            if predicted != actual:

                residue_failures += 1

    print(
        f"survival failures="
        f"{survival_failures}"
    )

    print(
        f"residue failures="
        f"{residue_failures}"
    )


# ==============================================================================
# SAT FORM
# ==============================================================================

def print_sat_model():

    print()
    print("=" * 90)
    print(
        "SAT MODEL"
    )
    print("=" * 90)

    print(
        """
GLOBAL VARIABLES:

    F = initial frame
        F=1 -> A
        F=0 -> B

    X_0, X_1, X_2, ...
        bits of |X|

    Y_0, Y_1, Y_2, ...
        bits of |Y|

At transition z -> z+1:

    use X_(z-1)
    use Y_(z-1)

Survival variable:

    S_z

Exact equation:

    S_z <-> (!X_(z-1) & !Y_(z-1))

CNF:

    (!S_z v !X_(z-1))
    (!S_z v !Y_(z-1))
    ( S_z v  X_(z-1) v  Y_(z-1))

Frame:

    F_(z+1) = F_z

whenever S_z.

No x_z or y_z integer variables are necessary.

The entire coordinate hierarchy can therefore be represented
by the bit streams of X and Y.
"""
    )


# ==============================================================================
# BIT STREAM REPORT
# ==============================================================================

def show_bit_streams(
    states,
):

    print()
    print("=" * 90)
    print(
        "GLOBAL BIT-STREAM REPRESENTATION"
    )
    print("=" * 90)

    shown = 0

    for state in states:

        if shown >= EXAMPLE_COUNT:
            break

        depth = (
            1
            + min(
                v2(state.X),
                v2(state.Y),
            )
        )

        print()
        print(
            f"n={state.n} "
            f"frame={state.frame} "
            f"X={state.X} "
            f"Y={state.Y}"
        )

        print(
            f"    X binary="
            f"{abs(state.X):b}"
        )

        print(
            f"    Y binary="
            f"{abs(state.Y):b}"
        )

        print(
            f"    v2(X)={v2(state.X)} "
            f"v2(Y)={v2(state.Y)} "
            f"deepest={depth}"
        )

        print(
            "    transition bits:"
        )

        for z in range(
            2,
            min(
                depth + 1,
                MAX_Z,
            ),
        ):

            xb = bit(
                state.X,
                z - 1,
            )

            yb = bit(
                state.Y,
                z - 1,
            )

            print(
                f"        "
                f"z={z}->{z+1} "
                f"bit=({xb},{yb}) "
                f"survive={xb == 0 and yb == 0}"
            )

        shown += 1


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 605 START"
    )
    print("=" * 90)

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = sieve_primes(
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
        f"    semiprimes="
        f"{len(semiprimes)}"
    )

    print()
    print(
        "[3] GLOBAL STATE"
    )

    states = build_global_states(
        semiprimes
    )

    print(
        f"    global states="
        f"{len(states)}"
    )

    test_frame_from_n(
        states
    )

    test_XY_forward_inverse(
        states
    )

    test_global_n(
        states
    )

    test_global_active_residue(
        states
    )

    test_survival(
        states
    )

    test_deepest(
        states
    )

    test_child_residues(
        states
    )

    test_frame_necessary_only_for_factor_mapping(
        states
    )

    print_sat_model()

    show_bit_streams(
        states
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 605 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
