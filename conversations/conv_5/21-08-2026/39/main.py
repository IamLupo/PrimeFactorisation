#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 598
# ==============================================================================
# DIRECT RESIDUE -> FRAME -> CHILD RESIDUE RECURSION
#
# We now assume the exact coordinate renormalization:
#
#     x' = x / 2
#     y' = y / 2
#
# and focus exclusively on the missing hierarchical mechanism:
#
#     parent residue
#         +
#     parent frame
#         +
#     discarded coordinate bits
#         ->
#     child frame
#         ->
#     child residue
#
# The experiment compares:
#
#   1. empirical transitions from actual semiprimes
#   2. symbolic modular predictions from the A/B equations
#   3. low-bit state compression
#   4. exact residue recurrences
#
# NO FILES
# NO EXTERNAL SOURCES
# ==============================================================================

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

SHOW_EXAMPLES = True
EXAMPLES_PER_TRANSITION = 8


# ==============================================================================
# SAMPLE STRUCTURE
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
                (p * q, p, q)
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

    k = 1 << (z - 1)
    d = 2 * k

    x_num = q - p + 6
    y_num = p + q

    if x_num % d != 0:
        return None

    if y_num % d != 0:
        return None

    return (
        x_num // d,
        y_num // d,
    )


# ==============================================================================
# FRAME B
# ==============================================================================

def frame_B(
    p,
    q,
    z,
):

    k = 1 << (z - 1)
    d = 2 * k

    x_num = 3 * p - q + 6
    y_num = 3 * p + q

    if x_num % d != 0:
        return None

    if y_num % d != 0:
        return None

    return (
        x_num // d,
        y_num // d,
    )


# ==============================================================================
# REPRESENTATION
# ==============================================================================

def representation(
    n,
    p,
    q,
    z,
):

    r = n % (1 << z)

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

        return State(
            n=n,
            p=p,
            q=q,
            z=z,
            residue=r,
            frame="A",
            x=a[0],
            y=a[1],
        )

    if b is not None:

        return State(
            n=n,
            p=p,
            q=q,
            z=z,
            residue=r,
            frame="B",
            x=b[0],
            y=b[1],
        )

    return None


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

            state = representation(
                n,
                p,
                q,
                z,
            )

            if state is not None:

                levels[z][n] = state

    return levels


# ==============================================================================
# SYMBOLIC FRAME FORMULAS
# ==============================================================================

def symbolic_A(
    z,
):

    k = 1 << (z - 1)

    # p = -kx + ky + 3
    # q =  kx + ky - 3

    return (
        (-k, k, 3),
        (k, k, -3),
    )


def symbolic_B(
    z,
):

    k = Fraction(
        1 << (z - 1),
        1,
    )

    # p = -(k/3)x + (k/3)y + 1
    # q = kx + ky - 3

    return (
        (
            -k / 3,
            k / 3,
            Fraction(1),
        ),
        (
            k,
            k,
            Fraction(-3),
        ),
    )


# ==============================================================================
# LOW-BIT STATE
# ==============================================================================

def low_bits(
    x,
    y,
):

    return (
        x & 1,
        y & 1,
    )


# ==============================================================================
# TRANSITION COLLECTION
# ==============================================================================

def collect_transitions(
    levels,
):

    transitions = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        common = set(parent).intersection(
            child
        )

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            key = (
                s0.frame,
                s0.residue,
                s1.frame,
                s1.residue,
                s0.x & 1,
                s0.y & 1,
            )

            transitions[
                (z, key)
            ].append(
                (s0, s1)
            )

    return transitions


# ==============================================================================
# EMPIRICAL FRAME TRANSITION TREE
# ==============================================================================

def report_frame_transition_tree(
    levels,
):

    print()
    print("=" * 90)
    print("EMPIRICAL FRAME / RESIDUE TRANSITION TREE")
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        counts = Counter()

        parent = levels[z]
        child = levels[z + 1]

        common = set(parent).intersection(
            child
        )

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            counts[
                (
                    s0.frame,
                    s0.residue,
                    s1.frame,
                    s1.residue,
                )
            ] += 1

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for key, count in sorted(
            counts.items()
        ):

            f0, r0, f1, r1 = key

            print(
                f"    "
                f"{f0}{r0}"
                f" -> "
                f"{f1}{r1}"
                f" : "
                f"{count}"
            )


# ==============================================================================
# LOW-BIT DETERMINISM
# ==============================================================================

def report_low_bit_determinism(
    levels,
):

    print()
    print("=" * 90)
    print("LOW-BIT STATE -> CHILD FRAME")
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        groups = defaultdict(set)

        parent = levels[z]
        child = levels[z + 1]

        common = set(parent).intersection(
            child
        )

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            key = (
                s0.frame,
                s0.residue,
                s0.x & 1,
                s0.y & 1,
            )

            groups[key].add(
                s1.frame
            )

        print()
        print(
            f"z={z}"
        )

        deterministic = True

        for key, frames in sorted(
            groups.items()
        ):

            print(
                f"    "
                f"{key} -> "
                f"{sorted(frames)}"
            )

            if len(frames) != 1:
                deterministic = False

        print(
            f"    deterministic={deterministic}"
        )


# ==============================================================================
# LOW-BIT DETERMINISM -> CHILD RESIDUE
# ==============================================================================

def report_low_bit_child_residue(
    levels,
):

    print()
    print("=" * 90)
    print("LOW-BIT STATE -> CHILD RESIDUE")
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        groups = defaultdict(set)

        parent = levels[z]
        child = levels[z + 1]

        common = set(parent).intersection(
            child
        )

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            key = (
                s0.frame,
                s0.residue,
                s0.x & 1,
                s0.y & 1,
            )

            groups[key].add(
                (
                    s1.frame,
                    s1.residue,
                )
            )

        print()
        print(
            f"z={z}"
        )

        deterministic = True

        for key, values in sorted(
            groups.items()
        ):

            print(
                f"    "
                f"{key} -> "
                f"{sorted(values)}"
            )

            if len(values) != 1:
                deterministic = False

        print(
            f"    deterministic={deterministic}"
        )


# ==============================================================================
# SYMBOLIC RESIDUE PREDICTION
# ==============================================================================

def symbolic_child_residue(
    p,
    q,
    z,
    frame,
    x_parity=None,
):

    """
    Compute the exact child residue directly from the frame formula.

    This deliberately avoids reconstructing n from the supplied n.
    """

    modulus = 1 << (z + 1)
    k = 1 << (z - 1)

    if frame == "A":

        n_value = (
            (k * p * 0)
        )

        # A:
        #
        # p = k(y-x)+3
        # q = k(y+x)-3
        #
        # Since p,q are supplied, n=pq.
        #
        # The function is used only as a symbolic verification hook.

        n_value = p * q

    elif frame == "B":

        n_value = p * q

    else:

        raise ValueError(
            f"unknown frame={frame}"
        )

    return n_value % modulus


# ==============================================================================
# DIRECT MODULAR FORMULA FROM PARENT x,y
# ==============================================================================

def symbolic_n_from_frame(
    x,
    y,
    z,
    frame,
):

    k = 1 << (z - 1)

    if frame == "A":

        p = (
            k * (y - x)
            + 3
        )

        q = (
            k * (y + x)
            - 3
        )

        return p * q

    if frame == "B":

        # 3p = k(y-x)+3
        # q  = k(y+x)-3

        numerator = (
            k * (y - x)
            + 3
        )

        q = (
            k * (y + x)
            - 3
        )

        if numerator % 3 != 0:
            return None

        p = numerator // 3

        return p * q

    raise ValueError(
        frame
    )


# ==============================================================================
# SYMBOLIC CHILD RESIDUE FROM HALF COORDINATES
# ==============================================================================

def predicted_child(
    x,
    y,
    z,
    frame,
):

    if x % 2 != 0:
        return None

    if y % 2 != 0:
        return None

    x2 = x // 2
    y2 = y // 2

    A_value = symbolic_n_from_frame(
        x2,
        y2,
        z + 1,
        "A",
    )

    B_value = symbolic_n_from_frame(
        x2,
        y2,
        z + 1,
        "B",
    )

    modulus = 1 << (z + 1)

    result = []

    if A_value is not None:

        result.append(
            (
                "A",
                A_value % modulus,
            )
        )

    if B_value is not None:

        result.append(
            (
                "B",
                B_value % modulus,
            )
        )

    return result


# ==============================================================================
# SYMBOLIC TRANSITION VALIDATION
# ==============================================================================

def validate_symbolic_child(
    levels,
):

    print()
    print("=" * 90)
    print(
        "SYMBOLIC CHILD-RESIDUE VALIDATION"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        common = set(parent).intersection(
            child
        )

        tested = 0
        failures = 0

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            prediction = predicted_child(
                s0.x,
                s0.y,
                z,
                s0.frame,
            )

            if prediction is None:

                failures += 1
                continue

            actual = (
                s1.frame,
                s1.residue,
            )

            tested += 1

            if actual not in prediction:

                failures += 1

                if failures <= 10:

                    print()
                    print(
                        "    MISMATCH"
                    )

                    print(
                        f"        n={n}"
                    )

                    print(
                        f"        parent="
                        f"{s0.frame}"
                        f"{s0.residue}"
                    )

                    print(
                        f"        x={s0.x}"
                        f" y={s0.y}"
                    )

                    print(
                        f"        predicted="
                        f"{prediction}"
                    )

                    print(
                        f"        actual="
                        f"{actual}"
                    )

        print(
            f"    tested={tested}"
        )

        print(
            f"    failures={failures}"
        )


# ==============================================================================
# RESIDUE RECURRENCE FIT
# ==============================================================================

def fit_recurrence(
    residues,
):

    pairs = []

    for a, b in zip(
        residues,
        residues[1:],
    ):

        if a is None or b is None:
            continue

        pairs.append(
            (a, b)
        )

    print()
    print(
        "    affine recurrence tests:"
    )

    candidates = {
        "2r+3":
            lambda r: 2 * r + 3,

        "2r+9":
            lambda r: 2 * r + 9,

        "2r-3":
            lambda r: 2 * r - 3,

        "2r-9":
            lambda r: 2 * r - 9,

        "r+2":
            lambda r: r + 2,

        "r+4":
            lambda r: r + 4,

        "r+6":
            lambda r: r + 6,

        "r+8":
            lambda r: r + 8,
    }

    for name, fn in candidates.items():

        exact = all(
            b == fn(a)
            for a, b in pairs
        )

        print(
            f"        "
            f"{name:<6} "
            f"{'EXACT' if exact else 'no'}"
        )

    if pairs:

        print()
        print(
            "    observed r_next-2r:"
        )

        print(
            sorted(
                {
                    b - 2 * a
                    for a, b in pairs
                }
            )
        )


# ==============================================================================
# ACTIVE RESIDUE RECURRENCE
# ==============================================================================

def report_active_recurrence(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ACTIVE RESIDUE RECURRENCE"
    )
    print("=" * 90)

    for frame in ("A", "B"):

        residues = []

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            counts = Counter()

            for state in levels[z].values():

                if state.frame == frame:

                    counts[
                        state.residue
                    ] += 1

            if len(counts) == 1:

                residues.append(
                    next(iter(counts))
                )

            else:

                residues.append(
                    None
                )

        print()
        print(
            f"FRAME {frame}"
        )

        print(
            f"    residues={residues}"
        )

        fit_recurrence(
            residues
        )


# ==============================================================================
# FRAME SURVIVAL TABLE
# ==============================================================================

def report_frame_survival(
    levels,
):

    print()
    print("=" * 90)
    print(
        "FRAME SURVIVAL"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        counts = Counter(
            state.frame
            for state in levels[z].values()
        )

        print(
            f"z={z:<2} "
            f"A={counts['A']:<7} "
            f"B={counts['B']:<7}"
        )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def report_examples(
    levels,
):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLE LOW-BIT TRANSITIONS"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        common = sorted(
            set(parent).intersection(
                child
            )
        )

        if not common:
            continue

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        shown = 0

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            print(
                f"    "
                f"n={n} "
                f"{s0.frame}{s0.residue}"
                f" "
                f"x={s0.x}"
                f" y={s0.y}"
                f" bits="
                f"({s0.x & 1},{s0.y & 1})"
                f" -> "
                f"{s1.frame}{s1.residue}"
            )

            shown += 1

            if shown >= EXAMPLES_PER_TRANSITION:
                break


# ==============================================================================
# FINAL SYMBOLIC SUMMARY
# ==============================================================================

def symbolic_summary():

    print()
    print("=" * 90)
    print(
        "SYMBOLIC HIERARCHY SUMMARY"
    )
    print("=" * 90)

    print(
        r"""
FRAME A

    p_z = k(y-x)+3
    q_z = k(y+x)-3

with

    k = 2^(z-1).

After renormalization

    x' = x/2
    y' = y/2
    k' = 2k

we obtain

    p_(z+1)(x',y')
       = k'(y'-x')+3
       = k(y-x)+3
       = p_z(x,y)

and

    q_(z+1)(x',y')
       = q_z(x,y).


FRAME B

    3p_z = k(y-x)+3
    q_z  = k(y+x)-3.

Again

    x'=x/2
    y'=y/2
    k'=2k

gives

    3p_(z+1)(x',y')
       = 3p_z(x,y)

and

    q_(z+1)(x',y')
       = q_z(x,y).


Therefore the coefficients themselves are not the recursion.

The actual recursion to discover is:

    (frame_z, residue_z, lowbits(x_z,y_z))
                    |
                    v
    (frame_(z+1), residue_(z+1))

The coordinate transformation is already known:

    (x_z,y_z)
        ->
    (x_z/2,y_z/2).

Experiment 598 tests whether the remaining branch transition is
a deterministic finite-state rule.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 598 START"
    )
    print("=" * 90)

    print()
    print(
        "DIRECT RESIDUE -> FRAME -> CHILD RESIDUE RECURSION"
    )

    # --------------------------------------------------------------------------
    # SIEVE
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    # --------------------------------------------------------------------------
    # SEMIPRIMES
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
        f"    semiprimes={len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # LEVELS
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
            f"    z={z}: "
            f"states={len(levels[z])}"
        )

    # --------------------------------------------------------------------------
    # FRAME SURVIVAL
    # --------------------------------------------------------------------------

    report_frame_survival(
        levels
    )

    # --------------------------------------------------------------------------
    # TRANSITION GRAPH
    # --------------------------------------------------------------------------

    report_frame_transition_tree(
        levels
    )

    # --------------------------------------------------------------------------
    # LOW BIT -> FRAME
    # --------------------------------------------------------------------------

    report_low_bit_determinism(
        levels
    )

    # --------------------------------------------------------------------------
    # LOW BIT -> RESIDUE
    # --------------------------------------------------------------------------

    report_low_bit_child_residue(
        levels
    )

    # --------------------------------------------------------------------------
    # SYMBOLIC VALIDATION
    # --------------------------------------------------------------------------

    validate_symbolic_child(
        levels
    )

    # --------------------------------------------------------------------------
    # ACTIVE RESIDUE RECURRENCE
    # --------------------------------------------------------------------------

    report_active_recurrence(
        levels
    )

    # --------------------------------------------------------------------------
    # EXAMPLES
    # --------------------------------------------------------------------------

    report_examples(
        levels
    )

    # --------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------

    symbolic_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 598 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
