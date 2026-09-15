#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 589
# ==============================================================================
# BRANCH-COMPATIBILITY OF THE LEVEL-INVARIANT COORDINATES
#
# NO FILES
# NO WEB
#
# Experiment 587 established:
#
#     X = 2^(z-1) * x_z
#     Y = 2^(z-1) * y_z
#
# and:
#
#     X,Y are invariant across levels.
#
# Experiment 588 showed:
#
#     2^(z-1) | X,Y
#
# is necessary for a coordinate representation,
# but NOT sufficient for the A/B function at that level.
#
# Example:
#
#     n=39
#     X=8
#     Y=8
#
# has enough 2-adic divisibility for z=4,
# but the actual construction is only valid through z=3.
#
#
# Therefore this experiment isolates the missing condition:
#
#     BRANCH COMPATIBILITY.
#
#
# At any level z:
#
#     A:
#         p_A = Y-X+3
#         q_A = Y+X-3
#
#     B:
#         3*p_B = Y+X-3
#         q_B  = Y-X+3
#
#
# The experiment tests both candidate branches independently.
#
# It also computes:
#
#     n_A = p_A*q_A
#     n_B = p_B*q_B
#
# and compares their residues modulo 2^z and 2^(z+1).
#
#
# MAIN QUESTION:
#
#     Can the valid branch be determined from low bits of X,Y?
#
# If yes, the deepest valid level should be expressible as:
#
#     divisibility condition
#         +
#     branch compatibility condition.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

MAX_BITS_SEARCH = 10

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
            start:limit + 1:p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_semiprimes(
    max_n,
    sieve,
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
# INVARIANT COORDINATES
# ==============================================================================

def invariant_from_pq(
    p,
    q,
):

    # --------------------------------------------------------------------------
    # A:
    #
    # p = Y-X+3
    # q = Y+X-3
    #
    # --------------------------------------------------------------------------

    if (
        (p + q) % 2 == 0
        and
        (q - p + 6) % 2 == 0
    ):

        Y = (
            p + q
        ) // 2

        X = (
            q - p + 6
        ) // 2

        if (
            Y - X + 3 == p
            and
            Y + X - 3 == q
        ):

            return FactorState(
                n=p * q,
                p=p,
                q=q,
                X=X,
                Y=Y,
            )

    # --------------------------------------------------------------------------
    # B:
    #
    # 3p = Y+X-3
    # q  = Y-X+3
    #
    # --------------------------------------------------------------------------

    if (
        (3 * p + q) % 2 == 0
        and
        (3 * p - q + 6) % 2 == 0
    ):

        Y = (
            3 * p + q
        ) // 2

        X = (
            3 * p - q + 6
        ) // 2

        if (
            3 * p
            == Y + X - 3
            and
            q
            == Y - X + 3
        ):

            return FactorState(
                n=p * q,
                p=p,
                q=q,
                X=X,
                Y=Y,
            )

    return None


# ==============================================================================
# BUILD INVARIANT STATES
# ==============================================================================

def build_states(
    semiprimes,
):

    result = []

    for n, p, q in semiprimes:

        state = invariant_from_pq(
            p,
            q,
        )

        if state is not None:
            result.append(state)

    return result


# ==============================================================================
# CANDIDATE BRANCH A
# ==============================================================================

def branch_a(
    X,
    Y,
):

    p = (
        Y - X + 3
    )

    q = (
        Y + X - 3
    )

    return p, q


# ==============================================================================
# CANDIDATE BRANCH B
# ==============================================================================

def branch_b(
    X,
    Y,
):

    numerator = (
        Y + X - 3
    )

    if numerator % 3 != 0:
        return None

    p = numerator // 3

    q = (
        Y - X + 3
    )

    return p, q


# ==============================================================================
# BRANCH MATCH
# ==============================================================================

def branch_match(
    state,
    branch,
):

    if branch == "A":

        result = branch_a(
            state.X,
            state.Y,
        )

    else:

        result = branch_b(
            state.X,
            state.Y,
        )

    if result is None:
        return False

    p, q = result

    return (
        p == state.p
        and
        q == state.q
    )


# ==============================================================================
# BRANCH RESIDUE
# ==============================================================================

def branch_residue(
    state,
    branch,
    z,
):

    modulus = 1 << z

    if branch == "A":

        result = branch_a(
            state.X,
            state.Y,
        )

    else:

        result = branch_b(
            state.X,
            state.Y,
        )

    if result is None:
        return None

    p, q = result

    return (
        p * q
    ) % modulus


# ==============================================================================
# ACTUAL LEVEL VALIDITY
# ==============================================================================

def actual_level(
    state,
    z,
):

    k = 1 << (
        z - 1
    )

    # Coordinate divisibility.
    if state.X % k != 0:
        return None

    if state.Y % k != 0:
        return None

    # Determine actual residue branch.
    modulus = 1 << z
    residue = state.n % modulus

    if residue == modulus - 1:

        if branch_match(
            state,
            "A",
        ):

            return "A"

        return None

    if branch_match(
        state,
        "B",
    ):

        return "B"

    return None


# ==============================================================================
# BRANCH COMPATIBILITY TABLE
# ==============================================================================

def compatibility_report(
    states,
):

    print()
    print("=" * 90)
    print(
        "BRANCH COMPATIBILITY BY LEVEL"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        counts = {
            "A": 0,
            "B": 0,
            "NONE": 0,
            "BOTH": 0,
        }

        for state in states:

            a_match = branch_match(
                state,
                "A",
            )

            b_match = branch_match(
                state,
                "B",
            )

            actual = actual_level(
                state,
                z,
            )

            if (
                a_match
                and
                b_match
            ):

                counts["BOTH"] += 1

            elif a_match:

                counts["A"] += 1

            elif b_match:

                counts["B"] += 1

            else:

                counts["NONE"] += 1

            # Verify actual construction.
            if (
                actual == "A"
                and
                not a_match
            ):
                raise RuntimeError(
                    "Actual A branch failed"
                )

            if (
                actual == "B"
                and
                not b_match
            ):
                raise RuntimeError(
                    "Actual B branch failed"
                )

        print()
        print(
            f"z={z}"
        )

        print(
            f"    A-compatible = "
            f"{counts['A']}"
        )

        print(
            f"    B-compatible = "
            f"{counts['B']}"
        )

        print(
            f"    BOTH          = "
            f"{counts['BOTH']}"
        )

        print(
            f"    NONE          = "
            f"{counts['NONE']}"
        )


# ==============================================================================
# LOW-BIT BRANCH PREDICTOR
# ==============================================================================

def low_bit_branch_predictor(
    states,
):

    print()
    print("=" * 90)
    print(
        "LOW-BIT BRANCH PREDICTOR"
    )
    print("=" * 90)

    for bits in range(
        1,
        MAX_BITS_SEARCH + 1,
    ):

        modulus = 1 << bits

        table = defaultdict(set)

        for state in states:

            key = (
                state.X % modulus,
                state.Y % modulus,
            )

            if branch_match(
                state,
                "A",
            ):

                label = "A"

            elif branch_match(
                state,
                "B",
            ):

                label = "B"

            else:

                label = "NONE"

            table[key].add(
                label
            )

        deterministic = all(
            len(values) == 1
            for values
            in table.values()
        )

        print()
        print(
            f"bits={bits:<2} "
            f"modulus={modulus:<5} "
            f"states={len(table):<7} "
            f"deterministic={deterministic}"
        )

        if deterministic:

            print(
                "    EXACT low-bit branch rule found."
            )

            # Do not dump enormous tables.
            if bits <= 4:

                for key in sorted(table):

                    print(
                        f"        {key} -> "
                        f"{sorted(table[key])}"
                    )

            print()


# ==============================================================================
# LEVEL COMPATIBILITY PREDICTOR
# ==============================================================================

def level_compatibility_table(
    states,
):

    print()
    print("=" * 90)
    print(
        "LEVEL COMPATIBILITY FROM X,Y LOW BITS"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        # At level z, only the low z-1 bits are relevant to whether
        # X,Y can be divided by 2^(z-1), plus branch compatibility.

        bits = z - 1

        modulus = 1 << bits

        table = defaultdict(set)

        for state in states:

            key = (
                state.X % modulus,
                state.Y % modulus,
            )

            valid = (
                actual_level(
                    state,
                    z,
                )
                is not None
            )

            table[key].add(
                valid
            )

        deterministic = all(
            len(values) == 1
            for values
            in table.values()
        )

        print(
            f"z={z:<2} "
            f"lowbits=2^{bits:<2} "
            f"states={len(table):<7} "
            f"deterministic={deterministic}"
        )


# ==============================================================================
# MODULAR FACTOR FORMULAS
# ==============================================================================

def modular_formula_report():

    print()
    print("=" * 90)
    print(
        "SYMBOLIC BRANCH RESIDUES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        modulus = 1 << z

        print()
        print(
            f"z={z}, modulus={modulus}"
        )

        print(
            "    A:"
        )

        print(
            "        n_A = "
            "(Y-X+3)(Y+X-3)"
        )

        print(
            "    B:"
        )

        print(
            "        n_B = "
            "(Y+X-3)(Y-X+3)/3"
        )

        # ----------------------------------------------------------------------
        # Difference between A and B.
        # ----------------------------------------------------------------------

        print(
            "    Difference:"
        )

        print(
            "        n_A = 3*n_B"
        )

        print(
            "        whenever B is integral."
        )

        print(
            "    Therefore the branch distinction is not simply"
        )

        print(
            "        divisibility of X,Y."
        )

        print(
            "    It also depends on which scaled factorization"
        )

        print(
            "    matches the actual n."
        )


# ==============================================================================
# DEEPEST VALID LEVEL
# ==============================================================================

def deepest_levels(
    states,
):

    print()
    print("=" * 90)
    print(
        "DEEPEST VALID LEVEL"
    )
    print("=" * 90)

    distribution = defaultdict(int)

    examples = []

    for state in states:

        deepest = None

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            if actual_level(
                state,
                z,
            ) is not None:

                deepest = z

        distribution[
            deepest
        ] += 1

        if (
            deepest is not None
            and
            len(examples) < MAX_EXAMPLES
        ):

            examples.append(
                (
                    state,
                    deepest,
                )
            )

    for level in sorted(
        distribution,
        key=lambda value: (
            value is None,
            value if value is not None else -1,
        ),
    ):

        print(
            f"    deepest={level}: "
            f"{distribution[level]}"
        )

    print()

    print(
        "    EXAMPLES"
    )

    for state, deepest in examples:

        print(
            f"        n={state.n} "
            f"p={state.p} "
            f"q={state.q} "
            f"X={state.X} "
            f"Y={state.Y} "
            f"deepest={deepest}"
        )


# ==============================================================================
# EXPLAIN MISMATCHES
# ==============================================================================

def explain_valuation_mismatch(
    states,
):

    print()
    print("=" * 90)
    print(
        "DIVISIBILITY-VALID BUT BRANCH-INVALID"
    )
    print("=" * 90)

    shown = 0

    for state in states:

        vx = 0
        value = abs(state.X)

        if value == 0:
            vx = 10**9
        else:
            while (
                value & 1
            ) == 0:
                value >>= 1
                vx += 1

        vy = 0
        value = abs(state.Y)

        if value == 0:
            vy = 10**9
        else:
            while (
                value & 1
            ) == 0:
                value >>= 1
                vy += 1

        candidate_level = min(
            vx,
            vy,
        ) + 1

        if candidate_level < MIN_Z:
            continue

        if candidate_level > MAX_Z:
            candidate_level = MAX_Z

        actual = actual_level(
            state,
            candidate_level,
        )

        if actual is None:

            print()
            print(
                f"    n={state.n} "
                f"p={state.p} "
                f"q={state.q}"
            )

            print(
                f"        X={state.X} "
                f"Y={state.Y}"
            )

            print(
                f"        v2X={vx} "
                f"v2Y={vy}"
            )

            print(
                f"        candidate level="
                f"{candidate_level}"
            )

            print(
                f"        actual branch at "
                f"candidate = NONE"
            )

            a_residue = branch_residue(
                state,
                "A",
                candidate_level,
            )

            b_residue = branch_residue(
                state,
                "B",
                candidate_level,
            )

            actual_residue = (
                state.n
                % (1 << candidate_level)
            )

            print(
                f"        n residue="
                f"{actual_residue}"
            )

            print(
                f"        A residue="
                f"{a_residue}"
            )

            print(
                f"        B residue="
                f"{b_residue}"
            )

            shown += 1

            if shown >= MAX_EXAMPLES:
                break


# ==============================================================================
# SAME-N LEVEL CHAIN
# ==============================================================================

def same_n_chain_report(
    states,
):

    print()
    print("=" * 90)
    print(
        "LEVEL CHAINS IN INVARIANT COORDINATES"
    )
    print("=" * 90)

    shown = 0

    for state in states:

        chain = []

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            branch = actual_level(
                state,
                z,
            )

            if branch is None:
                break

            chain.append(
                (
                    z,
                    branch,
                )
            )

        if not chain:
            continue

        print()
        print(
            f"n={state.n} "
            f"X={state.X} "
            f"Y={state.Y}"
        )

        print(
            "    chain = "
            + " -> ".join(
                f"{z}:{branch}"
                for z, branch
                in chain
            )
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
        "EXPERIMENT 589 START"
    )
    print("=" * 90)

    print()
    print(
        "BRANCH-COMPATIBILITY OF INVARIANT X,Y"
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
    # Invariant states.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] INVARIANT STATES"
    )

    states = build_states(
        semiprimes
    )

    print(
        f"    states={len(states)}"
    )

    # --------------------------------------------------------------------------
    # Branch compatibility.
    # --------------------------------------------------------------------------

    compatibility_report(
        states
    )

    # --------------------------------------------------------------------------
    # Low-bit rule.
    # --------------------------------------------------------------------------

    low_bit_branch_predictor(
        states
    )

    # --------------------------------------------------------------------------
    # Level compatibility.
    # --------------------------------------------------------------------------

    level_compatibility_table(
        states
    )

    # --------------------------------------------------------------------------
    # Symbolic formulas.
    # --------------------------------------------------------------------------

    modular_formula_report()

    # --------------------------------------------------------------------------
    # Deepest level.
    # --------------------------------------------------------------------------

    deepest_levels(
        states
    )

    # --------------------------------------------------------------------------
    # Explain mismatches.
    # --------------------------------------------------------------------------

    explain_valuation_mismatch(
        states
    )

    # --------------------------------------------------------------------------
    # Chains.
    # --------------------------------------------------------------------------

    same_n_chain_report(
        states
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
Experiment 588 showed:

    divisibility:
        2^(z-1) | X,Y

is not sufficient for a level representation.

Experiment 589 isolates the second condition:

    BRANCH COMPATIBILITY.

For every level we independently evaluate:

    A:
        p = Y-X+3
        q = Y+X-3

    B:
        3p = Y+X-3
        q  = Y-X+3

and compare them against the actual factor pair.

The desired result is a compact rule of the form:

    level z
        +
    low bits of X,Y
        ->
    valid branch

If such a rule exists, the complete hierarchy becomes:

    invariant X,Y
        |
        +-- divisibility filter
        |
        +-- branch-compatibility filter
        |
        +-- next modular level

That would provide a much more precise SAT reduction model than
using 2-adic valuation alone.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 589 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
