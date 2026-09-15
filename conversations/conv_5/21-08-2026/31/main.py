#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 588
# ==============================================================================
# 2-ADIC INVARIANT COORDINATE SURVIVAL
#
# NO FILES
# NO WEB
#
# From Experiment 587:
#
#     X = 2^(z-1) * x_z
#     Y = 2^(z-1) * y_z
#
# and X,Y are invariant across levels.
#
# Therefore:
#
#     x_z = X / 2^(z-1)
#     y_z = Y / 2^(z-1)
#
# For the SAME factor pair to survive from level z to z+1,
# both child coordinates must remain integral:
#
#     x_(z+1) = X / 2^z
#     y_(z+1) = Y / 2^z
#
# hence:
#
#     2^z | X
#     2^z | Y
#
# or equivalently:
#
#     v2(X) >= z
#     v2(Y) >= z
#
#
# MAIN HYPOTHESIS
# ---------------
#
# The maximum surviving level should therefore be controlled by:
#
#     min(v2(X), v2(Y)).
#
# This experiment:
#
#   1. generates semiprimes;
#   2. computes invariant X,Y at level z=2;
#   3. computes v2(X), v2(Y);
#   4. predicts the maximum possible level;
#   5. compares prediction with the actual highest observed level;
#   6. tests exact transition survival;
#   7. analyzes branch/residue information after the divisibility filter;
#   8. reports how many SAT states remain after each additional bit.
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
MAX_Z = 12

SHOW_EXAMPLES = True
MAX_EXAMPLES = 20


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class FactorState:
    n: int
    p: int
    q: int
    X: int
    Y: int


@dataclass(frozen=True)
class LevelState:
    n: int
    p: int
    q: int

    z: int
    k: int
    branch: str

    x: int
    y: int

    X: int
    Y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> bytearray:

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = isqrt(limit)

    for p in range(2, root + 1):

        if not sieve[p]:
            continue

        start = p * p
        count = (limit - start) // p + 1

        sieve[start:limit + 1:p] = (
            b"\x00" * count
        )

    return sieve


# ==============================================================================
# SEMIPRIME GENERATION
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
# 2-ADIC VALUATION
# ==============================================================================

def v2(n: int) -> int:
    """
    Return the exponent of 2 dividing n.

    Convention:
        v2(0) = infinity
    """

    if n == 0:
        return 10**9

    n = abs(n)

    count = 0

    while (
        n & 1
    ) == 0:

        n >>= 1
        count += 1

    return count


# ==============================================================================
# INVARIANT COORDINATES FROM FACTORS
# ==============================================================================
#
# At z=2:
#
#     k=2
#
# Case A:
#
#     p = 2(y-x)+3
#     q = 2(y+x)-3
#
# Therefore:
#
#     X = kx
#     Y = ky
#
# gives:
#
#     p = Y-X+3
#     q = Y+X-3.
#
#
# Case B:
#
#     3p = Y+X-3
#     q  = Y-X+3.
#
# We recover X,Y directly from p,q.
#
# ==============================================================================

def invariant_from_pq(
    p: int,
    q: int,
):

    # --------------------------------------------------------------------------
    # Candidate A
    # --------------------------------------------------------------------------

    # p = Y-X+3
    # q = Y+X-3
    #
    # Therefore:
    #
    # Y-X = p-3
    # Y+X = q+3
    #
    # Hence:
    #
    # 2Y = p+q
    # 2X = q-p+6

    if (
        (p + q) % 2 == 0
        and
        (q - p + 6) % 2 == 0
    ):

        Y_a = (
            p + q
        ) // 2

        X_a = (
            q - p + 6
        ) // 2

        # Validate.

        if (
            Y_a - X_a + 3 == p
            and
            Y_a + X_a - 3 == q
        ):

            return "A", X_a, Y_a

    # --------------------------------------------------------------------------
    # Candidate B
    # --------------------------------------------------------------------------

    # 3p = Y+X-3
    # q  = Y-X+3
    #
    # Therefore:
    #
    # Y+X = 3p+3
    # Y-X = q-3
    #
    # Hence:
    #
    # 2Y = 3p+q
    # 2X = 3p-q+6

    if (
        (3 * p + q) % 2 == 0
        and
        (3 * p - q + 6) % 2 == 0
    ):

        Y_b = (
            3 * p + q
        ) // 2

        X_b = (
            3 * p - q + 6
        ) // 2

        if (
            3 * p == Y_b + X_b - 3
            and
            q == Y_b - X_b + 3
        ):

            return "B", X_b, Y_b

    return None


# ==============================================================================
# BUILD INVARIANT STATES
# ==============================================================================

def build_invariant_states(
    semiprimes,
):

    states = []

    for n, p, q in semiprimes:

        result = invariant_from_pq(
            p,
            q,
        )

        if result is None:
            continue

        _, X, Y = result

        states.append(
            FactorState(
                n=n,
                p=p,
                q=q,
                X=X,
                Y=Y,
            )
        )

    return states


# ==============================================================================
# PREDICT MAX LEVEL
# ==============================================================================

def predicted_max_level(
    X: int,
    Y: int,
):

    """
    Level z requires:

        2^(z-1) | X
        2^(z-1) | Y

    So the largest possible z is:

        min(v2(X), v2(Y)) + 1

    assuming the coordinates themselves are defined at z.
    """

    vx = v2(X)
    vy = v2(Y)

    finite = min(
        vx,
        vy,
    )

    if finite >= 10**8:
        return 10**9

    return finite + 1


# ==============================================================================
# ACTUAL LEVEL CONSTRUCTION
# ==============================================================================

def actual_level_state(
    state: FactorState,
    z: int,
):

    k = 1 << (
        z - 1
    )

    if (
        state.X % k != 0
        or
        state.Y % k != 0
    ):
        return None

    x = state.X // k
    y = state.Y // k

    modulus = 1 << z
    residue = state.n % modulus

    # --------------------------------------------------------------------------
    # Determine which factor function applies.
    # --------------------------------------------------------------------------

    if residue == modulus - 1:

        # A:
        #
        # p = Y-X+3
        # q = Y+X-3

        if (
            state.p
            != state.Y - state.X + 3
        ):
            return None

        if (
            state.q
            != state.Y + state.X - 3
        ):
            return None

        branch = "A"

    else:

        # B:
        #
        # 3p = Y+X-3
        # q  = Y-X+3

        if (
            3 * state.p
            != state.Y + state.X - 3
        ):
            return None

        if (
            state.q
            != state.Y - state.X + 3
        ):
            return None

        branch = "B"

    return LevelState(
        n=state.n,
        p=state.p,
        q=state.q,
        z=z,
        k=k,
        branch=branch,
        x=x,
        y=y,
        X=state.X,
        Y=state.Y,
    )


# ==============================================================================
# ACTUAL HIGHEST LEVEL
# ==============================================================================

def actual_highest_level(
    state: FactorState,
):

    highest = None

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        if actual_level_state(
            state,
            z,
        ) is not None:

            highest = z

    return highest


# ==============================================================================
# PREDICTION AUDIT
# ==============================================================================

def prediction_audit(
    invariant_states,
):

    print()
    print("=" * 90)
    print(
        "PREDICTED VS ACTUAL MAXIMUM LEVEL"
    )
    print("=" * 90)

    exact = 0
    mismatch = 0

    mismatch_examples = []

    for state in invariant_states:

        predicted = predicted_max_level(
            state.X,
            state.Y,
        )

        actual = actual_highest_level(
            state
        )

        if actual == predicted:

            exact += 1

        else:

            mismatch += 1

            if len(
                mismatch_examples
            ) < MAX_EXAMPLES:

                mismatch_examples.append(
                    (
                        state,
                        predicted,
                        actual,
                    )
                )

    print(
        f"    total states = "
        f"{len(invariant_states)}"
    )

    print(
        f"    exact        = "
        f"{exact}"
    )

    print(
        f"    mismatches   = "
        f"{mismatch}"
    )

    if mismatch_examples:

        print()
        print(
            "    FIRST MISMATCHES"
        )

        for state, predicted, actual in mismatch_examples:

            print(
                f"        n={state.n} "
                f"p={state.p} "
                f"q={state.q} "
                f"X={state.X} "
                f"Y={state.Y} "
                f"v2X={v2(state.X)} "
                f"v2Y={v2(state.Y)} "
                f"predicted={predicted} "
                f"actual={actual}"
            )


# ==============================================================================
# SURVIVAL BY DIVISIBILITY
# ==============================================================================

def divisibility_survival(
    invariant_states,
):

    print()
    print("=" * 90)
    print(
        "2-ADIC SURVIVAL COUNTS"
    )
    print("=" * 90)

    print()
    print(
        "A state survives level z iff:"
    )

    print(
        "    2^(z-1) | X"
    )

    print(
        "    2^(z-1) | Y"
    )

    print()

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        required = 1 << (
            z - 1
        )

        count = 0

        for state in invariant_states:

            if (
                state.X % required == 0
                and
                state.Y % required == 0
            ):

                count += 1

        print(
            f"    z={z:<2} "
            f"required=2^{z-1:<2} "
            f"survivors={count}"
        )


# ==============================================================================
# TRANSITION AUDIT
# ==============================================================================

def transition_audit(
    invariant_states,
):

    print()
    print("=" * 90)
    print(
        "EXACT z -> z+1 SURVIVAL TEST"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        required_parent = 1 << (
            z - 1
        )

        required_child = 1 << z

        predicted_survivors = 0
        actual_survivors = 0
        disagreement = 0

        for state in invariant_states:

            parent_exists = (
                state.X
                % required_parent
                == 0
                and
                state.Y
                % required_parent
                == 0
            )

            child_exists = (
                state.X
                % required_child
                == 0
                and
                state.Y
                % required_child
                == 0
            )

            actual_parent = (
                actual_level_state(
                    state,
                    z,
                )
                is not None
            )

            actual_child = (
                actual_level_state(
                    state,
                    z + 1,
                )
                is not None
            )

            predicted_transition = (
                parent_exists
                and
                child_exists
            )

            actual_transition = (
                actual_parent
                and
                actual_child
            )

            if predicted_transition:
                predicted_survivors += 1

            if actual_transition:
                actual_survivors += 1

            if (
                predicted_transition
                !=
                actual_transition
            ):

                disagreement += 1

        print()
        print(
            f"    z={z} -> z={z+1}"
        )

        print(
            f"        predicted={predicted_survivors}"
        )

        print(
            f"        actual={actual_survivors}"
        )

        print(
            f"        disagreement={disagreement}"
        )


# ==============================================================================
# v2 DISTRIBUTION
# ==============================================================================

def valuation_distribution(
    invariant_states,
):

    print()
    print("=" * 90)
    print(
        "v2(X), v2(Y) DISTRIBUTION"
    )
    print("=" * 90)

    pairs = defaultdict(int)

    for state in invariant_states:

        vx = v2(state.X)
        vy = v2(state.Y)

        pairs[
            (
                vx,
                vy,
            )
        ] += 1

    for (vx, vy), count in sorted(
        pairs.items()
    ):

        print(
            f"    v2(X)={vx:<3} "
            f"v2(Y)={vy:<3} "
            f"count={count}"
        )


# ==============================================================================
# BRANCH VS LOW BITS
# ==============================================================================

def low_bit_branch_test(
    invariant_states,
):

    print()
    print("=" * 90)
    print(
        "BRANCH / LOW-BIT STRUCTURE OF X,Y"
    )
    print("=" * 90)

    for bits in range(
        1,
        7,
    ):

        modulus = 1 << bits

        table = defaultdict(set)

        for state in invariant_states:

            Xr = state.X % modulus
            Yr = state.Y % modulus

            level = actual_highest_level(
                state
            )

            branch = "NONE"

            if level is not None:

                current = actual_level_state(
                    state,
                    min(
                        level,
                        MAX_Z,
                    ),
                )

                if current is not None:
                    branch = current.branch

            table[
                (
                    Xr,
                    Yr,
                )
            ].add(
                branch
            )

        deterministic = all(
            len(values) == 1
            for values in table.values()
        )

        print()
        print(
            f"    bits={bits} "
            f"modulus={modulus} "
            f"states={len(table)} "
            f"deterministic={deterministic}"
        )


# ==============================================================================
# SAT REDUCTION
# ==============================================================================

def sat_reduction_report(
    invariant_states,
):

    print()
    print("=" * 90)
    print(
        "POTENTIAL SAT STATE REDUCTION"
    )
    print("=" * 90)

    total = len(
        invariant_states
    )

    print(
        f"    base invariant states = {total}"
    )

    print()

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        required = 1 << (
            z - 1
        )

        survivors = [
            state
            for state in invariant_states
            if (
                state.X % required == 0
                and
                state.Y % required == 0
            )
        ]

        print(
            f"    level z={z:<2} "
            f"survivors={len(survivors):<8} "
            f"fraction="
            f"{len(survivors) / total:.10f}"
        )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    invariant_states,
):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLES WITH 2-ADIC VALUATION"
    )
    print("=" * 90)

    shown = 0

    for state in invariant_states:

        vx = v2(
            state.X
        )

        vy = v2(
            state.Y
        )

        predicted = predicted_max_level(
            state.X,
            state.Y,
        )

        actual = actual_highest_level(
            state
        )

        print(
            f"    n={state.n} "
            f"p={state.p} "
            f"q={state.q} "
            f"X={state.X} "
            f"Y={state.Y} "
            f"v2X={vx} "
            f"v2Y={vy} "
            f"predicted={predicted} "
            f"actual={actual}"
        )

        shown += 1

        if shown >= MAX_EXAMPLES:
            break


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 588 START"
    )
    print("=" * 90)

    print()
    print(
        "2-ADIC INVARIANT COORDINATE SURVIVAL"
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
        f"    semiprimes = "
        f"{len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Invariant coordinates.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD INVARIANT COORDINATES"
    )

    invariant_states = build_invariant_states(
        semiprimes
    )

    print(
        f"    invariant states = "
        f"{len(invariant_states)}"
    )

    # --------------------------------------------------------------------------
    # Prediction audit.
    # --------------------------------------------------------------------------

    prediction_audit(
        invariant_states
    )

    # --------------------------------------------------------------------------
    # Divisibility survival.
    # --------------------------------------------------------------------------

    divisibility_survival(
        invariant_states
    )

    # --------------------------------------------------------------------------
    # Transition audit.
    # --------------------------------------------------------------------------

    transition_audit(
        invariant_states
    )

    # --------------------------------------------------------------------------
    # Valuations.
    # --------------------------------------------------------------------------

    valuation_distribution(
        invariant_states
    )

    # --------------------------------------------------------------------------
    # Branch / low-bit structure.
    # --------------------------------------------------------------------------

    low_bit_branch_test(
        invariant_states
    )

    # --------------------------------------------------------------------------
    # SAT reduction.
    # --------------------------------------------------------------------------

    sat_reduction_report(
        invariant_states
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    examples(
        invariant_states
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
The proposed recursive coordinate system is:

    X = 2^(z-1) * x_z
    Y = 2^(z-1) * y_z

A transition from z to z+1 requires:

    x_(z+1) = X / 2^z
    y_(z+1) = Y / 2^z

so the exact survival condition is:

    2^z | X
    2^z | Y.

Equivalently:

    v2(X) >= z
    v2(Y) >= z.

Therefore the predicted deepest level is:

    max_z =
        min(v2(X), v2(Y)) + 1.

Experiment 588 tests this prediction against the actual
level-dependent algebra.

If the prediction has zero mismatches, then the entire
shrinking-coordinate hierarchy can be interpreted as a
2-adic divisibility filtration of the invariant coordinates.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 588 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
