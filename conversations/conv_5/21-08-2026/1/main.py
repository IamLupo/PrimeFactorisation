#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 571
# ==============================================================================
# HIERARCHICAL MODULAR x,y STRUCTURE
#
# Goal
# ----
#
# Test the hypothesis:
#
#     mod 4 rules
#          |
#          +---- determine mod 8 rules
#                    |
#                    +---- determine mod 16 rules
#                              |
#                              +---- determine mod 32 rules
#
#
# Odd residue tree:
#
#          mod 4
#         /     \
#        1       3
#       / \     / \
#      1   5   3   7        mod 8
#
# Each node is then split again at the next power of two.
#
#
# 2020 coordinate system used
# ---------------------------
#
#     A0 = 2y - 2x + 3
#     B0 = 2y + 2x - 3
#
#     P = A0 / r
#     Q = B0 / r
#
# Therefore:
#
#     P + Q = 4y / r
#
#     P - Q = (6 - 4x) / r
#
# and
#
#     H_P = P - 2Q - 1/2
#
#     C_P = 2P - 4Q - 1.
#
#
# IMPORTANT
# ---------
#
# This experiment does NOT claim that a particular x,y,n relationship from
# the original 2020 experiment has been reconstructed.
#
# Instead, several controlled constructions are tested and their modular
# behavior is compared.
#
# The experiment can therefore reveal whether the hierarchical phenomenon
# survives different natural interpretations of the coordinates.
#
# ==============================================================================


from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 100_000

# Only odd n.
ONLY_PRIMES = False

# For composite n, only use semiprime n=p*q if True.
ONLY_SEMIPRIMES = True

# Maximum samples displayed per node.
SHOW_SAMPLES = 5

# Levels:
#
#   2 -> mod 4
#   3 -> mod 8
#   4 -> mod 16
#   5 -> mod 32
#
MAX_LEVEL = 5


# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass(frozen=True)
class Sample:
    n: int
    p: int
    q: int
    r: int
    x: int
    y: int

    P: Fraction
    Q: Fraction

    S: Fraction
    gap: Fraction
    Delta: Fraction
    H_P: Fraction
    C_P: Fraction


# ==============================================================================
# BASIC NUMBER THEORY
# ==============================================================================

def is_prime(n: int) -> bool:

    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:

        if n % d == 0:
            return False

        d += 2

    return True


def smallest_factor(n: int) -> int:

    if n % 2 == 0:
        return 2

    d = 3

    while d * d <= n:

        if n % d == 0:
            return d

        d += 2

    return n


def semiprime_factor_pair(n: int):

    p = smallest_factor(n)

    if p == n:
        return None

    q = n // p

    if is_prime(q):
        return p, q

    return None


# ==============================================================================
# 2020 CROSSWALK
# ==============================================================================

def xy_from_pq(p: int, q: int):
    """
    Controlled reconstruction of x,y from a factor pair.

    We deliberately expose the normalization r.

    For the simplest normalization r=1:

        2y - 2x + 3 = p
        2y + 2x - 3 = q

    Adding:

        4y = p+q

        y = (p+q)/4

    Subtracting:

        4x = q-p+6

        x = (q-p+6)/4

    These need not be integers for arbitrary p,q.

    Therefore this function only accepts pairs for which the coordinates
    are integral.

    The experiment additionally tests r=2 and r=4.
    """

    candidates = []

    for r in (1, 2, 4):

        y_num = p + q

        x_num = q - p + 6

        if y_num % (4 * r) != 0:
            continue

        if x_num % (4 * r) != 0:
            continue

        y = y_num // (4 * r)
        x = x_num // (4 * r)

        candidates.append((r, x, y))

    return candidates


# ==============================================================================
# ROOT COORDINATES
# ==============================================================================

def make_sample(n: int, p: int, q: int, r: int, x: int, y: int):

    P = Fraction(
        2 * y - 2 * x + 3,
        r,
    )

    Q = Fraction(
        2 * y + 2 * x - 3,
        r,
    )

    S = P + Q

    gap = P - Q

    Delta = gap * gap

    H_P = P - 2 * Q - Fraction(1, 2)

    C_P = 2 * P - 4 * Q - 1

    return Sample(
        n=n,
        p=p,
        q=q,
        r=r,
        x=x,
        y=y,
        P=P,
        Q=Q,
        S=S,
        gap=gap,
        Delta=Delta,
        H_P=H_P,
        C_P=C_P,
    )


# ==============================================================================
# DATA GENERATION
# ==============================================================================

def generate_samples():

    samples = []

    for n in range(3, MAX_N + 1, 2):

        if ONLY_PRIMES:

            if not is_prime(n):
                continue

            # No factor pair exists for a prime.
            continue

        if ONLY_SEMIPRIMES:

            pair = semiprime_factor_pair(n)

            if pair is None:
                continue

            p, q = pair

        else:

            # Generic factor pair.
            p = smallest_factor(n)

            if p == n:
                continue

            q = n // p

            if p * q != n:
                continue

        candidates = xy_from_pq(
            p,
            q,
        )

        for r, x, y in candidates:

            sample = make_sample(
                n,
                p,
                q,
                r,
                x,
                y,
            )

            samples.append(sample)

    return samples


# ==============================================================================
# RESIDUE HELPERS
# ==============================================================================

def modulus(level):

    return 2 ** level


def odd_residues(m):

    return range(
        1,
        m,
        2,
    )


def parent_of(residue, level):

    m = modulus(level)

    return residue % (m // 2)


def children_of(residue, level):

    m = modulus(level)

    return (
        residue,
        residue + m,
    )


# ==============================================================================
# GROUPING
# ==============================================================================

def group_samples(samples, level):

    m = modulus(level)

    groups = defaultdict(list)

    for s in samples:

        groups[
            s.n % m
        ].append(s)

    return groups


# ==============================================================================
# SIGNATURES
# ==============================================================================

def unique_values(samples, attribute):

    return {
        getattr(s, attribute)
        for s in samples
    }


def signature(samples):

    keys = [
        "n",
        "p",
        "q",
        "r",
        "x",
        "y",
        "P",
        "Q",
        "S",
        "gap",
        "Delta",
        "H_P",
        "C_P",
    ]

    result = {}

    for key in keys:

        values = unique_values(
            samples,
            key,
        )

        result[key] = {
            "count": len(values),
            "values": values,
        }

    return result


# ==============================================================================
# AFFINE RELATION
# ==============================================================================

def exact_affine_relation(
    samples,
    source,
    target,
):

    if len(samples) < 2:
        return None

    pairs = [
        (
            getattr(s, source),
            getattr(s, target),
        )
        for s in samples
    ]

    # Need two different source values.
    first = None

    for i in range(len(pairs)):

        x1, y1 = pairs[i]

        for j in range(i + 1, len(pairs)):

            x2, y2 = pairs[j]

            if x1 == x2:
                continue

            a = Fraction(
                y2 - y1,
                x2 - x1,
            )

            b = y1 - a * x1

            first = (a, b)

            break

        if first is not None:
            break

    if first is None:
        return None

    a, b = first

    for x, y in pairs:

        if a * x + b != y:
            return None

    return a, b


# ==============================================================================
# CONDITIONAL MODULAR RULE
# ==============================================================================

def modular_signature(samples, modulus):

    """
    For every observable determine its distribution modulo modulus.

    This is more useful than simply looking at n mod modulus.
    """

    result = {}

    keys = [
        "x",
        "y",
        "r",
        "P",
        "Q",
        "S",
        "gap",
        "Delta",
        "H_P",
        "C_P",
    ]

    for key in keys:

        counter = Counter()

        for s in samples:

            value = getattr(s, key)

            if isinstance(value, Fraction):

                if value.denominator == 1:

                    value = value.numerator

                else:
                    # Keep numerator/denominator structure separate.
                    value = (
                        value.numerator % modulus,
                        value.denominator % modulus,
                    )

            else:

                value %= modulus

            counter[value] += 1

        result[key] = counter

    return result


# ==============================================================================
# LEVEL REPORT
# ==============================================================================

def report_level(samples, level):

    m = modulus(level)

    print()
    print("=" * 90)
    print(f"LEVEL {level}: MOD {m}")
    print("=" * 90)

    groups = group_samples(
        samples,
        level,
    )

    for residue in odd_residues(m):

        subset = groups.get(
            residue,
            [],
        )

        print()
        print(
            f"[{residue:>2} mod {m:<2}] "
            f"samples={len(subset)}"
        )

        if not subset:
            continue

        for key in [
            "x",
            "y",
            "r",
            "P",
            "Q",
            "S",
            "gap",
            "Delta",
            "H_P",
            "C_P",
        ]:

            values = unique_values(
                subset,
                key,
            )

            print(
                f"    {key:<7}: "
                f"{len(values):>5} unique"
            )

        for s in subset[:SHOW_SAMPLES]:

            print(
                "    "
                f"n={s.n:<8} "
                f"p={s.p:<6} "
                f"q={s.q:<6} "
                f"r={s.r:<2} "
                f"x={s.x:<6} "
                f"y={s.y:<6}"
            )


# ==============================================================================
# PARENT / CHILD ANALYSIS
# ==============================================================================

def analyze_transition(
    samples,
    level,
):

    child_level = level + 1

    parent_m = modulus(level)
    child_m = modulus(child_level)

    parents = group_samples(
        samples,
        level,
    )

    children = group_samples(
        samples,
        child_level,
    )

    print()
    print("=" * 90)
    print(
        f"TRANSITION: MOD {parent_m} -> MOD {child_m}"
    )
    print("=" * 90)

    for parent in odd_residues(parent_m):

        parent_samples = parents.get(
            parent,
            [],
        )

        child_a, child_b = children_of(
            parent,
            level,
        )

        A = children.get(
            child_a,
            [],
        )

        B = children.get(
            child_b,
            [],
        )

        print()
        print(
            f"PARENT {parent} mod {parent_m}"
        )

        print(
            f"    child A = {child_a} mod {child_m}: "
            f"{len(A)} samples"
        )

        print(
            f"    child B = {child_b} mod {child_m}: "
            f"{len(B)} samples"
        )

        if not A or not B:
            print(
                "    insufficient child data"
            )
            continue

        # ----------------------------------------------------------------------
        # Compare exact observable sets.
        # ----------------------------------------------------------------------

        for key in [
            "x",
            "y",
            "P",
            "Q",
            "S",
            "gap",
            "Delta",
            "H_P",
            "C_P",
        ]:

            values_A = unique_values(
                A,
                key,
            )

            values_B = unique_values(
                B,
                key,
            )

            same = values_A == values_B

            overlap = len(
                values_A & values_B
            )

            print(
                f"    {key:<7} "
                f"A={len(values_A):<5} "
                f"B={len(values_B):<5} "
                f"overlap={overlap:<5} "
                f"same={same}"
            )


# ==============================================================================
# CHILD DIFFERENCE TEST
# ==============================================================================

def child_difference_test(
    samples,
    level,
):

    parent_m = modulus(level)
    child_m = modulus(level + 1)

    groups = group_samples(
        samples,
        level,
    )

    child_groups = group_samples(
        samples,
        level + 1,
    )

    print()
    print("=" * 90)
    print(
        f"CHILD DIFFERENCE TEST "
        f"MOD {parent_m} -> MOD {child_m}"
    )
    print("=" * 90)

    for parent in odd_residues(parent_m):

        A_id, B_id = children_of(
            parent,
            level,
        )

        A = child_groups.get(
            A_id,
            [],
        )

        B = child_groups.get(
            B_id,
            [],
        )

        if len(A) < 2 or len(B) < 2:
            continue

        print()
        print(
            f"Parent {parent} mod {parent_m}"
        )

        # --------------------------------------------------------------
        # Find observables whose distributions differ between children.
        # --------------------------------------------------------------

        discriminators = []

        for key in [
            "x",
            "y",
            "r",
            "P",
            "Q",
            "S",
            "gap",
            "Delta",
            "H_P",
            "C_P",
        ]:

            set_A = unique_values(
                A,
                key,
            )

            set_B = unique_values(
                B,
                key,
            )

            if set_A != set_B:

                discriminators.append(
                    key
                )

        if discriminators:

            print(
                "    CHILD-DISCRIMINATING OBSERVABLES:"
            )

            print(
                "    "
                + ", ".join(
                    discriminators
                )
            )

        else:

            print(
                "    No observable separates the two children."
            )


# ==============================================================================
# INHERITED RULE TEST
# ==============================================================================

def inherited_rule_test(
    samples,
    level,
):

    parent_m = modulus(level)

    child_m = modulus(level + 1)

    parent_groups = group_samples(
        samples,
        level,
    )

    child_groups = group_samples(
        samples,
        level + 1,
    )

    print()
    print("=" * 90)
    print(
        f"INHERITED-RULE TEST "
        f"MOD {parent_m} -> MOD {child_m}"
    )
    print("=" * 90)

    observables = [
        "x",
        "y",
        "P",
        "Q",
        "S",
        "gap",
        "Delta",
        "H_P",
        "C_P",
    ]

    for parent in odd_residues(parent_m):

        child_A_id, child_B_id = children_of(
            parent,
            level,
        )

        A = child_groups.get(
            child_A_id,
            [],
        )

        B = child_groups.get(
            child_B_id,
            [],
        )

        if len(A) < 3 or len(B) < 3:
            continue

        print()
        print(
            f"Parent {parent} mod {parent_m}"
        )

        found = False

        for source in observables:

            for target in observables:

                rel_A = exact_affine_relation(
                    A,
                    source,
                    target,
                )

                rel_B = exact_affine_relation(
                    B,
                    source,
                    target,
                )

                if (
                    rel_A is not None
                    and rel_B is not None
                    and rel_A == rel_B
                ):

                    a, b = rel_A

                    print(
                        f"    SHARED CHILD LAW:"
                        f" {target} = "
                        f"({a})*{source} + ({b})"
                    )

                    found = True

        if not found:

            print(
                "    No shared exact affine child law found."
            )


# ==============================================================================
# CONDITIONAL RESIDUE TABLE
# ==============================================================================

def conditional_table(
    samples,
    parent_level,
    observable,
):

    parent_m = modulus(parent_level)

    child_m = modulus(parent_level + 1)

    groups = group_samples(
        samples,
        parent_level + 1,
    )

    print()
    print("=" * 90)
    print(
        f"CONDITIONAL TABLE: "
        f"{observable}"
    )
    print(
        f"MOD {parent_m} parent -> MOD {child_m} child"
    )
    print("=" * 90)

    for parent in odd_residues(parent_m):

        A_id, B_id = children_of(
            parent,
            parent_level,
        )

        A = groups.get(
            A_id,
            [],
        )

        B = groups.get(
            B_id,
            [],
        )

        print()

        print(
            f"parent {parent} mod {parent_m}:"
        )

        for label, subset in [
            (A_id, A),
            (B_id, B),
        ]:

            if not subset:
                continue

            residues = Counter()

            for s in subset:

                value = getattr(
                    s,
                    observable,
                )

                if isinstance(value, Fraction):

                    if value.denominator == 1:

                        value = value.numerator

                    else:

                        residues[
                            (
                                value.numerator,
                                value.denominator,
                            )
                        ] += 1

                        continue

                residues[
                    value % child_m
                ] += 1

            print(
                f"    child {label}: "
                f"{dict(residues)}"
            )


# ==============================================================================
# TREE DISPLAY
# ==============================================================================

def show_tree(max_level):

    print()
    print("=" * 90)
    print("ODD RESIDUE TREE")
    print("=" * 90)

    for level in range(
        2,
        max_level + 1,
    ):

        m = modulus(level)
        pm = m // 2

        print()
        print(
            f"mod {pm} -> mod {m}"
        )

        for parent in odd_residues(pm):

            a, b = children_of(
                parent,
                level - 1,
            )

            print(
                f"    {parent:>2} "
                f"-> {a:>2}, {b:>2}"
            )


# ==============================================================================
# BIT INTERPRETATION
# ==============================================================================

def show_bits(max_level):

    print()
    print("=" * 90)
    print("2-ADIC / BINARY INTERPRETATION")
    print("=" * 90)

    for level in range(
        2,
        max_level + 1,
    ):

        m = modulus(level)

        print()
        print(
            f"MOD {m}"
        )

        for r in odd_residues(m):

            bits = format(
                r,
                f"0{level}b",
            )

            print(
                f"    {r:>3} = {bits}"
            )


# ==============================================================================
# SUMMARY
# ==============================================================================

def final_summary(samples):

    print()
    print("=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    print(
        f"Total generated samples: {len(samples)}"
    )

    if not samples:
        return

    print()

    print(
        "The experiment is testing the following hypothesis:"
    )

    print()
    print(
        "    R_4"
    )

    print(
        "     / \\"
    )

    print(
        "    R_8a R_8b"
    )

    print(
        "      /     \\"
    )

    print(
        "     R_16 children"
    )

    print(
        "         ..."
    )

    print()

    print(
        "A genuine hierarchical law requires more than different"
    )

    print(
        "frequencies in different residue classes."
    )

    print()

    print(
        "Evidence becomes interesting when:"
    )

    print(
        "    1. the two children of a parent are distinguishable,"
    )

    print(
        "    2. the distinction is expressed through x/y/P/Q observables,"
    )

    print(
        "    3. the child rules depend on the parent branch,"
    )

    print(
        "    4. the same transformation repeats at the next level."
    )

    print()

    print(
        "That would support a recursive:"
    )

    print()
    print(
        "    mod 4 -> mod 8 -> mod 16 -> mod 32"
    )

    print(
        "    rule hierarchy."
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 571 START")
    print("=" * 90)
    print()
    print(
        "HIERARCHICAL MODULAR x,y STRUCTURE"
    )
    print()

    # --------------------------------------------------------------------------
    # Generate.
    # --------------------------------------------------------------------------

    print("[1] GENERATING EXACT INTEGER DATA")
    print("-" * 90)

    samples = generate_samples()

    print(
        f"Generated {len(samples)} samples."
    )

    if not samples:

        print(
            "No valid samples."
        )

        return

    # --------------------------------------------------------------------------
    # Levels.
    # --------------------------------------------------------------------------

    print(
        "[2] LEVEL REPORTS"
    )

    for level in range(
        2,
        MAX_LEVEL + 1,
    ):

        report_level(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Parent / child.
    # --------------------------------------------------------------------------

    print(
        "[3] PARENT -> CHILD TESTS"
    )

    for level in range(
        2,
        MAX_LEVEL,
    ):

        analyze_transition(
            samples,
            level,
        )

        child_difference_test(
            samples,
            level,
        )

        inherited_rule_test(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Conditional tables.
    # --------------------------------------------------------------------------

    print(
        "[4] CONDITIONAL OBSERVABLE TABLES"
    )

    for observable in [
        "x",
        "y",
        "S",
        "gap",
        "H_P",
        "C_P",
    ]:

        conditional_table(
            samples,
            2,
            observable,
        )

    # --------------------------------------------------------------------------
    # Tree.
    # --------------------------------------------------------------------------

    print(
        "[5] RESIDUE TREE"
    )

    show_tree(
        MAX_LEVEL
    )

    # --------------------------------------------------------------------------
    # Binary interpretation.
    # --------------------------------------------------------------------------

    print(
        "[6] BINARY / 2-ADIC BRANCHES"
    )

    show_bits(
        MAX_LEVEL
    )

    # --------------------------------------------------------------------------
    # Final.
    # --------------------------------------------------------------------------

    final_summary(
        samples
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 571 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()