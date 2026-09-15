#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 600
# ==============================================================================
# ODD-BIT DESCENDANTS + HALF-OFFSET FRAME TRANSFORMATION
#
# Purpose
# -------
# Experiment 599 established:
#
#   for z >= 3:
#
#       A:
#           n mod 2^(z+1)
#             = -9 + 2^z*(x mod 2)
#
#       B:
#           n mod 2^(z+1)
#             = -3 + 2^z*(x mod 2)
#
# and that the active child exists exactly when x,y are even.
#
# This experiment investigates the discarded odd-x states.
#
# Instead of assuming that
#
#     x' = x/2
#     y' = y/2
#
# is the only possible coordinate transition, search:
#
#     x' = (x - ex)/2
#     y' = (y - ey)/2
#
# where:
#
#     ex, ey in {0,1}
#
# and test every possible parent-frame -> child-frame combination:
#
#     A -> A
#     A -> B
#     B -> A
#     B -> B
#
# The goal is to determine whether the odd branch is:
#
#     1. simply discarded,
#     2. represented by another frame,
#     3. represented after a half-offset transformation,
#     4. or produces a residue that has no valid normalized representation.
#
# No external data.
# No files.
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

SHOW_MISMATCHES = 12
SHOW_EXAMPLES = 12


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
# SEMIPRIMES
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
# FRAME A COORDINATES
# ==============================================================================

def frame_A(
    p,
    q,
    z,
):

    k = 1 << (z - 1)

    denominator = 2 * k

    x_num = (
        q - p + 6
    )

    y_num = (
        p + q
    )

    if x_num % denominator != 0:
        return None

    if y_num % denominator != 0:
        return None

    return (
        x_num // denominator,
        y_num // denominator,
    )


# ==============================================================================
# FRAME B COORDINATES
# ==============================================================================

def frame_B(
    p,
    q,
    z,
):

    k = 1 << (z - 1)

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

    if x_num % denominator != 0:
        return None

    if y_num % denominator != 0:
        return None

    return (
        x_num // denominator,
        y_num // denominator,
    )


# ==============================================================================
# BUILD LEVELS
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

            residue = (
                n % (1 << z)
            )

            a = frame_A(
                p,
                q,
                z,
            )

            b = frame_B(
                p,
                q,
                z,
            )

            if a is not None:

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=residue,
                    frame="A",
                    x=a[0],
                    y=a[1],
                )

            elif b is not None:

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=residue,
                    frame="B",
                    x=b[0],
                    y=b[1],
                )

    return levels


# ==============================================================================
# CHILD STATE
# ==============================================================================

def get_child(
    levels,
    z,
    n,
):

    return levels[
        z + 1
    ].get(n)


# ==============================================================================
# EXACT RESIDUE FORMULA
# ==============================================================================

def predicted_child_residue(
    state,
):

    z = state.z

    modulus = 1 << (
        z + 1
    )

    k = 1 << (
        z - 1
    )

    x = state.x
    y = state.y

    # Exact calculation for z=2 and above.
    #
    # This avoids making the z>=3 simplification at z=2.

    if state.frame == "A":

        value = (
            (k * (y - x) + 3)
            *
            (k * (y + x) - 3)
        )

        return value % modulus

    if state.frame == "B":

        numerator = (
            (k * (x + y) - 3)
            *
            (k * (y - x) + 3)
        )

        if numerator % 3 != 0:
            return None

        return (
            numerator // 3
        ) % modulus

    raise ValueError(
        state.frame
    )


# ==============================================================================
# CLOSED FORM FOR z >= 3
# ==============================================================================

def predicted_child_residue_closed(
    state,
):

    z = state.z

    modulus = 1 << (
        z + 1
    )

    parity = state.x & 1

    if state.frame == "A":

        return (
            -9
            + (1 << z) * parity
        ) % modulus

    if state.frame == "B":

        return (
            -3
            + (1 << z) * parity
        ) % modulus

    raise ValueError(
        state.frame
    )


# ==============================================================================
# HALF-OFFSET TRANSFORM
# ==============================================================================

def transformed_coordinate(
    x,
    y,
    ex,
    ey,
):

    x_num = x - ex
    y_num = y - ey

    if x_num % 2 != 0:
        return None

    if y_num % 2 != 0:
        return None

    return (
        x_num // 2,
        y_num // 2,
    )


# ==============================================================================
# FRAME EQUATION EVALUATION
# ==============================================================================

def evaluate_frame(
    frame,
    x,
    y,
    z,
):

    k = 1 << (
        z - 1
    )

    if frame == "A":

        p = (
            k * (y - x)
            + 3
        )

        q = (
            k * (y + x)
            - 3
        )

        return p, q

    if frame == "B":

        numerator = (
            k * (x + y)
            - 3
        )

        if numerator % 3 != 0:
            return None

        p = numerator // 3

        q = (
            k * (y - x)
            + 3
        )

        return p, q

    raise ValueError(
        frame
    )


# ==============================================================================
# [1] CLOSED FORM AUDIT
# ==============================================================================

def closed_form_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "CLOSED-FORM CHILD RESIDUE AUDIT"
    )
    print("=" * 90)

    for z in range(
        2,
        MAX_Z,
    ):

        exact_fail = 0
        closed_fail = 0

        for state in levels[z].values():

            actual = (
                state.n
                % (1 << (z + 1))
            )

            exact = (
                predicted_child_residue(
                    state
                )
            )

            if exact != actual:

                exact_fail += 1

            if z >= 3:

                closed = (
                    predicted_child_residue_closed(
                        state
                    )
                )

                if closed != actual:

                    closed_fail += 1

        print(
            f"z={z} -> z={z+1}"
            f"  exact_fail={exact_fail}"
            f"  closed_fail={closed_fail}"
        )


# ==============================================================================
# [2] EVEN / ODD RESIDUE LAW
# ==============================================================================

def parity_residue_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "PARITY -> CHILD RESIDUE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        groups = defaultdict(
            set
        )

        for state in levels[z].values():

            bit = state.x & 1

            residue = (
                state.n
                % (1 << (z + 1))
            )

            groups[
                (
                    state.frame,
                    bit,
                )
            ].add(
                residue
            )

        print()
        print(
            f"z={z}"
        )

        for key in sorted(groups):

            print(
                f"    "
                f"{key}"
                f" -> "
                f"{sorted(groups[key])}"
            )


# ==============================================================================
# [3] ACTIVE CHILD TEST
# ==============================================================================

def active_child_test(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ACTIVE CHILD VS PARITY"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        counts = Counter()

        for n, state in levels[z].items():

            child = get_child(
                levels,
                z,
                n,
            )

            key = (
                state.frame,
                state.x & 1,
                state.y & 1,
                child is not None,
            )

            counts[key] += 1

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for key, count in sorted(
            counts.items()
        ):

            print(
                f"    "
                f"{key} = {count}"
            )


# ==============================================================================
# [4] HALF-OFFSET SEARCH
# ==============================================================================

def half_offset_search(
    levels,
):

    print()
    print("=" * 90)
    print(
        "HALF-OFFSET FRAME TRANSFORMATION SEARCH"
    )
    print("=" * 90)

    transformations = []

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent_level = levels[z]
        child_level = levels[z + 1]

        for parent_frame in (
            "A",
            "B",
        ):

            # We test all four coordinate offsets.

            for ex in (
                0,
                1,
            ):

                for ey in (
                    0,
                    1,
                ):

                    for child_frame in (
                        "A",
                        "B",
                    ):

                        total = 0
                        coordinate_valid = 0
                        exact_factor = 0

                        for n, parent in (
                            parent_level.items()
                        ):

                            if (
                                parent.frame
                                != parent_frame
                            ):
                                continue

                            total += 1

                            transformed = (
                                transformed_coordinate(
                                    parent.x,
                                    parent.y,
                                    ex,
                                    ey,
                                )
                            )

                            if (
                                transformed
                                is None
                            ):
                                continue

                            coordinate_valid += 1

                            x2, y2 = (
                                transformed
                            )

                            factor_pair = (
                                evaluate_frame(
                                    child_frame,
                                    x2,
                                    y2,
                                    z + 1,
                                )
                            )

                            if factor_pair is None:
                                continue

                            p2, q2 = (
                                factor_pair
                            )

                            if (
                                p2 == parent.p
                                and
                                q2 == parent.q
                            ):

                                exact_factor += 1

                        transformations.append(
                            (
                                z,
                                parent_frame,
                                ex,
                                ey,
                                child_frame,
                                total,
                                coordinate_valid,
                                exact_factor,
                            )
                        )

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        rows = [
            row
            for row in transformations
            if row[0] == z
        ]

        rows.sort(
            key=lambda row:
            (
                -row[7],
                -row[6],
                row[1],
                row[2],
                row[3],
                row[4],
            )
        )

        for row in rows[:8]:

            (
                _z,
                pf,
                ex,
                ey,
                cf,
                total,
                coordinate_valid,
                exact_factor,
            ) = row

            print(
                f"    "
                f"{pf} -> {cf} "
                f"x'=(x-{ex})/2 "
                f"y'=(y-{ey})/2 "
                f"eligible={coordinate_valid}/"
                f"{total} "
                f"factor_exact={exact_factor}"
            )


# ==============================================================================
# [5] ODD BRANCH DETAILED ANALYSIS
# ==============================================================================

def odd_branch_analysis(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ODD BRANCH ANALYSIS"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for frame in (
            "A",
            "B",
        ):

            odd_states = [
                state
                for state in levels[z].values()
                if (
                    state.frame == frame
                    and (
                        state.x & 1
                    )
                )
            ]

            child_exists = sum(
                get_child(
                    levels,
                    z,
                    state.n,
                )
                is not None
                for state in odd_states
            )

            print(
                f"    frame={frame}"
                f" odd_x={len(odd_states)}"
                f" active_child={child_exists}"
            )

            residues = Counter()

            for state in odd_states:

                residue = (
                    state.n
                    % (1 << (z + 1))
                )

                residues[
                    residue
                ] += 1

            if residues:

                top = sorted(
                    residues.items()
                )

                print(
                    f"        child residues="
                    f"{top[:12]}"
                )


# ==============================================================================
# [6] DETERMINE WHETHER ODD CHILD IS OTHER FRAME
# ==============================================================================

def odd_other_frame_test(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ODD PARENT -> OTHER FRAME TEST"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for parent_frame in (
            "A",
            "B",
        ):

            candidates = [
                state
                for state in levels[z].values()
                if (
                    state.frame
                    == parent_frame
                    and (
                        state.x & 1
                    )
                )
            ]

            target_counts = Counter()

            for state in candidates:

                child = get_child(
                    levels,
                    z,
                    state.n,
                )

                if child is None:

                    target_counts[
                        "NONE"
                    ] += 1

                else:

                    target_counts[
                        child.frame
                    ] += 1

            print(
                f"    "
                f"{parent_frame}: "
                f"{dict(target_counts)}"
            )


# ==============================================================================
# [7] INACTIVE RESIDUE VERIFICATION
# ==============================================================================

def inactive_residue_test(
    levels,
):

    print()
    print("=" * 90)
    print(
        "INACTIVE CHILD RESIDUE TEST"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        active_residues = {
            state.residue
            for state
            in levels[z + 1].values()
        }

        observed_odd = defaultdict(
            int
        )

        for state in levels[z].values():

            if (
                state.x & 1
            ) == 0:
                continue

            residue = (
                state.n
                % (1 << (z + 1))
            )

            observed_odd[
                (
                    state.frame,
                    residue,
                )
            ] += 1

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for (
            frame,
            residue,
        ), count in sorted(
            observed_odd.items()
        )[:16]:

            status = (
                "ACTIVE"
                if residue in active_residues
                else "INACTIVE"
            )

            print(
                f"    "
                f"{frame} "
                f"residue={residue} "
                f"count={count} "
                f"{status}"
            )


# ==============================================================================
# [8] SAME-N REPRESENTATION CHECK
# ==============================================================================

def representation_count(
    levels,
):

    print()
    print("=" * 90)
    print(
        "NUMBER OF VALID CHILD REPRESENTATIONS"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        counts = Counter()

        for n, state in levels[z].items():

            a = frame_A(
                state.p,
                state.q,
                z + 1,
            )

            b = frame_B(
                state.p,
                state.q,
                z + 1,
            )

            count = (
                int(a is not None)
                +
                int(b is not None)
            )

            counts[
                (
                    state.frame,
                    state.x & 1,
                    count,
                )
            ] += 1

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for key, count in sorted(
            counts.items()
        ):

            print(
                f"    "
                f"{key} = {count}"
            )


# ==============================================================================
# [9] EXAMPLES
# ==============================================================================

def examples(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ODD-BRANCH EXAMPLES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        shown = 0

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for n in sorted(
            levels[z]
        ):

            if shown >= SHOW_EXAMPLES:
                break

            state = levels[z][n]

            if (
                state.x & 1
            ) == 0:
                continue

            residue = (
                n
                % (1 << (z + 1))
            )

            child = get_child(
                levels,
                z,
                n,
            )

            print(
                f"    "
                f"n={n} "
                f"{state.frame}"
                f"{state.residue} "
                f"x={state.x} "
                f"y={state.y} "
                f"bits="
                f"({state.x & 1},"
                f"{state.y & 1}) "
                f"child_residue="
                f"{residue} "
                f"child="
                f"{child.frame + str(child.residue) if child else 'NONE'}"
            )

            shown += 1


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary():

    print()
    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
        r"""
Experiment 600 distinguishes two separate phenomena.

1. MODULAR DESCENT

For z >= 3:

    A:
        r_(z+1)
          =
        -9 + 2^z(x mod 2)
        mod 2^(z+1)

    B:
        r_(z+1)
          =
        -3 + 2^z(x mod 2)
        mod 2^(z+1).

Thus the discarded x-bit determines the child residue.

2. REPRESENTATION SURVIVAL

The normalized child coordinates normally require:

    x' = x/2
    y' = y/2.

Therefore an active child representation requires:

    x even
    y even.

The experiment now asks whether the odd branch can instead use:

    x'=(x-1)/2
    y'=(y-1)/2

or other half-offset combinations.

The four possibilities tested are:

    (x/2, y/2)
    ((x-1)/2, y/2)
    (x/2, (y-1)/2)
    ((x-1)/2, (y-1)/2)

and each is tested against BOTH frames.

This is the direct test of whether the apparently discarded branch
contains another A/B representation.

The expected possibilities are therefore:

    odd branch -> other frame
    odd branch -> shifted frame
    odd branch -> inactive residue
    odd branch -> no representation at this level.

The experiment does not assume which one is correct.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 600 START"
    )
    print("=" * 90)

    print()
    print(
        "ODD-BIT DESCENDANTS + HALF-OFFSET FRAME TRANSFORMATION"
    )

    # --------------------------------------------------------------------------
    # Sieve
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    # --------------------------------------------------------------------------
    # Semiprimes
    # --------------------------------------------------------------------------

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

    # --------------------------------------------------------------------------
    # Levels
    # --------------------------------------------------------------------------

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
            f"    z={z}"
            f" states={len(levels[z])}"
        )

    # --------------------------------------------------------------------------
    # Audits
    # --------------------------------------------------------------------------

    closed_form_audit(
        levels
    )

    parity_residue_audit(
        levels
    )

    active_child_test(
        levels
    )

    half_offset_search(
        levels
    )

    odd_branch_analysis(
        levels
    )

    odd_other_frame_test(
        levels
    )

    inactive_residue_test(
        levels
    )

    representation_count(
        levels
    )

    examples(
        levels
    )

    final_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 600 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
