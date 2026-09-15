#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 572
# ==============================================================================
# HIERARCHICAL MODULAR BRANCH-PREDICTION EXPERIMENT
#
# Main question
# ------------
#
# Does a modular parent branch determine the next modular branch?
#
# For odd n:
#
#     n mod 4
#         |
#         +--> n mod 8
#                 |
#                 +--> n mod 16
#                         |
#                         +--> n mod 32
#
# At every level there is exactly one newly introduced binary bit.
#
# For example:
#
#     n mod 4 = 3
#
# gives two possibilities:
#
#     n mod 8 = 3
#     n mod 8 = 7
#
# Write:
#
#     n = r + 2^k b
#
# where:
#
#     r = n mod 2^k
#     b = 0 or 1
#
# Then b is the "next branch bit".
#
#
# This experiment asks:
#
#     Can b be predicted exactly from:
#
#         x
#         y
#         r
#         P
#         Q
#         S
#         P-Q
#         Delta
#         H_P
#         C_P
#
# and does the prediction depend on the parent residue?
#
#
# IMPORTANT
# ---------
#
# We retain the same controlled 2020-style coordinate construction used in
# Experiment 571:
#
#     A0 = 2y - 2x + 3
#     B0 = 2y + 2x - 3
#
#     P = A0/r
#     Q = B0/r
#
# with integer factor-pair construction.
#
# No tautological tests such as:
#
#     x = 1*x + 0
#
# are used here.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 100_000

# Levels tested:
#
#   2 -> mod 4
#   3 -> mod 8
#   4 -> mod 16
#   5 -> mod 32
#
MIN_LEVEL = 2
MAX_LEVEL = 5

# Maximum number of examples printed per rule.
SAMPLE_LIMIT = 8


# ==============================================================================
# DATA STRUCTURE
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
# NUMBER THEORY
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


def semiprime_pair(n: int):

    p = smallest_factor(n)

    if p == n:
        return None

    q = n // p

    if is_prime(q):
        return p, q

    return None


# ==============================================================================
# 2020 x,y CROSSWALK
# ==============================================================================

def xy_from_pq(p: int, q: int):

    """
    Solve:

        2y - 2x + 3 = r P
        2y + 2x - 3 = r Q

    using the same controlled normalization family as Experiment 571.

    For the raw factor coordinates:

        y = (p+q)/(4r)
        x = (q-p+6)/(4r)

    We keep r in {1,2,4}.
    """

    result = []

    for r in (1, 2, 4):

        y_num = p + q
        x_num = q - p + 6

        denominator = 4 * r

        if y_num % denominator != 0:
            continue

        if x_num % denominator != 0:
            continue

        y = y_num // denominator
        x = x_num // denominator

        result.append(
            (r, x, y)
        )

    return result


# ==============================================================================
# ROOT COORDINATES
# ==============================================================================

def build_sample(
    n,
    p,
    q,
    r,
    x,
    y,
):

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
# GENERATE DATA
# ==============================================================================

def generate_samples():

    samples = []

    for n in range(
        3,
        MAX_N + 1,
        2,
    ):

        pair = semiprime_pair(n)

        if pair is None:
            continue

        p, q = pair

        for r, x, y in xy_from_pq(p, q):

            samples.append(
                build_sample(
                    n,
                    p,
                    q,
                    r,
                    x,
                    y,
                )
            )

    return samples


# ==============================================================================
# MODULAR HELPERS
# ==============================================================================

def modulus(level):

    return 2 ** level


def parent_residue(
    n,
    level,
):

    return n % modulus(level)


def next_branch_bit(
    n,
    level,
):

    """
    Given parent level k:

        parent = n mod 2^k

    next level is:

        n mod 2^(k+1)

    The new bit is:

        floor(n / 2^k) mod 2
    """

    return (
        n // modulus(level)
    ) & 1


def child_residue(
    n,
    level,
):

    """
    The child residue at level k+1.
    """

    return n % modulus(level + 1)


# ==============================================================================
# FRACTION MODULAR REPRESENTATION
# ==============================================================================

def mod_value(
    value,
    m,
):

    """
    Represent exact rational values canonically.

    For integer values use value mod m.

    For non-integral rationals preserve:
        numerator mod m
        denominator mod m
    """

    if isinstance(value, Fraction):

        if value.denominator == 1:

            return value.numerator % m

        return (
            value.numerator % m,
            value.denominator % m,
        )

    return value % m


# ==============================================================================
# OBSERVABLE ACCESS
# ==============================================================================

OBSERVABLES = (
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
)


def observable_value(
    sample,
    name,
):

    return getattr(
        sample,
        name,
    )


# ==============================================================================
# GROUP BY PARENT
# ==============================================================================

def group_by_parent(
    samples,
    level,
):

    groups = defaultdict(list)

    for sample in samples:

        parent = parent_residue(
            sample.n,
            level,
        )

        groups[parent].append(sample)

    return groups


# ==============================================================================
# BASIC PARENT REPORT
# ==============================================================================

def report_parent_bits(
    samples,
    level,
):

    m = modulus(level)

    groups = group_by_parent(
        samples,
        level,
    )

    print()
    print("=" * 90)
    print(
        f"PARENT LEVEL: MOD {m}"
    )
    print(
        f"NEXT LEVEL:   MOD {2*m}"
    )
    print("=" * 90)

    for parent in sorted(groups):

        subset = groups[parent]

        bit_counts = Counter(
            next_branch_bit(
                s.n,
                level,
            )
            for s in subset
        )

        print()
        print(
            f"PARENT {parent} mod {m}"
        )

        print(
            f"    total = {len(subset)}"
        )

        print(
            f"    next bit 0 = {bit_counts[0]}"
        )

        print(
            f"    next bit 1 = {bit_counts[1]}"
        )

        if bit_counts[0] and bit_counts[1]:

            print(
                "    STATUS = both children exist"
            )

        elif bit_counts[0]:

            print(
                "    STATUS = only child bit 0"
            )

        elif bit_counts[1]:

            print(
                "    STATUS = only child bit 1"
            )


# ==============================================================================
# PURE OBSERVABLE MOD-2 PREDICTION
# ==============================================================================

def observable_mod2_partition(
    samples,
    parent_level,
):

    """
    For each observable v, test whether:

        v mod 2

    determines the next branch bit.
    """

    parent_groups = group_by_parent(
        samples,
        parent_level,
    )

    print()
    print("=" * 90)
    print(
        f"ONE-OBSERVABLE BIT PREDICTION"
    )
    print(
        f"Parent mod {modulus(parent_level)}"
        f" -> child mod {2 * modulus(parent_level)}"
    )
    print("=" * 90)

    for parent in sorted(parent_groups):

        subset = parent_groups[parent]

        if len(subset) < 4:
            continue

        print()
        print(
            f"PARENT {parent} mod "
            f"{modulus(parent_level)}"
        )

        for observable in OBSERVABLES:

            table = defaultdict(set)

            for sample in subset:

                value = observable_value(
                    sample,
                    observable,
                )

                if isinstance(
                    value,
                    Fraction,
                ):

                    if value.denominator == 1:
                        key = value.numerator & 1
                    else:
                        key = (
                            value.numerator & 1,
                            value.denominator & 1,
                        )

                else:

                    key = value & 1

                bit = next_branch_bit(
                    sample.n,
                    parent_level,
                )

                table[key].add(bit)

            deterministic = all(
                len(bits) == 1
                for bits in table.values()
            )

            if deterministic:

                rules = {
                    key: next(iter(bits))
                    for key, bits in table.items()
                }

                print(
                    f"    EXACT: "
                    f"next_bit = f({observable} mod 2)"
                    f" -> {rules}"
                )


# ==============================================================================
# MODULAR SIGNATURE PREDICTION
# ==============================================================================

def search_single_feature_rules(
    samples,
    parent_level,
):

    """
    Search:

        next_bit = f(observable mod m)

    for several small moduli.

    This is the main non-trivial test.
    """

    parent_m = modulus(parent_level)

    parent_groups = group_by_parent(
        samples,
        parent_level,
    )

    test_moduli = (
        2,
        4,
        8,
        16,
    )

    print()
    print("=" * 90)
    print(
        "SINGLE-FEATURE MODULAR RULE SEARCH"
    )
    print(
        f"Parent level = mod {parent_m}"
    )
    print(
        f"Next bit = child bit at mod {2 * parent_m}"
    )
    print("=" * 90)

    for parent in sorted(parent_groups):

        subset = parent_groups[parent]

        if len(subset) < 8:
            continue

        print()
        print(
            f"PARENT {parent} mod {parent_m}"
        )

        found = False

        for observable in OBSERVABLES:

            for test_modulus in test_moduli:

                table = defaultdict(set)

                for sample in subset:

                    raw = observable_value(
                        sample,
                        observable,
                    )

                    key = mod_value(
                        raw,
                        test_modulus,
                    )

                    bit = next_branch_bit(
                        sample.n,
                        parent_level,
                    )

                    table[key].add(bit)

                if not table:
                    continue

                deterministic = all(
                    len(bits) == 1
                    for bits in table.values()
                )

                if not deterministic:
                    continue

                rule = {
                    key: next(iter(bits))
                    for key, bits in sorted(
                        table.items(),
                        key=str,
                    )
                }

                print()
                print(
                    f"    EXACT RULE:"
                )

                print(
                    f"        parent = {parent} mod {parent_m}"
                )

                print(
                    f"        feature = {observable}"
                )

                print(
                    f"        modulus = {test_modulus}"
                )

                print(
                    f"        rule = {rule}"
                )

                found = True

        if not found:

            print(
                "    No exact single-feature rule found."
            )


# ==============================================================================
# TWO-FEATURE RULE SEARCH
# ==============================================================================

def search_two_feature_rules(
    samples,
    parent_level,
):

    """
    Search exact deterministic rules of the form:

        next_bit = f(A mod m, B mod m)

    This is important because the x,y system may require two coordinates.
    """

    parent_m = modulus(parent_level)

    groups = group_by_parent(
        samples,
        parent_level,
    )

    candidate_features = list(
        OBSERVABLES
    )

    test_moduli = (
        2,
        4,
        8,
    )

    print()
    print("=" * 90)
    print(
        "TWO-FEATURE MODULAR RULE SEARCH"
    )
    print(
        f"Parent mod {parent_m}"
    )
    print("=" * 90)

    for parent in sorted(groups):

        subset = groups[parent]

        if len(subset) < 8:
            continue

        print()
        print(
            f"PARENT {parent} mod {parent_m}"
        )

        found = False

        for feature_a, feature_b in combinations(
            candidate_features,
            2,
        ):

            for test_modulus in test_moduli:

                table = defaultdict(set)

                for sample in subset:

                    value_a = observable_value(
                        sample,
                        feature_a,
                    )

                    value_b = observable_value(
                        sample,
                        feature_b,
                    )

                    key_a = mod_value(
                        value_a,
                        test_modulus,
                    )

                    key_b = mod_value(
                        value_b,
                        test_modulus,
                    )

                    key = (
                        key_a,
                        key_b,
                    )

                    bit = next_branch_bit(
                        sample.n,
                        parent_level,
                    )

                    table[key].add(bit)

                if not table:
                    continue

                if not all(
                    len(bits) == 1
                    for bits in table.values()
                ):
                    continue

                print()
                print(
                    "    EXACT TWO-FEATURE RULE"
                )

                print(
                    f"        features = "
                    f"({feature_a}, {feature_b})"
                )

                print(
                    f"        modulus = "
                    f"{test_modulus}"
                )

                print(
                    f"        states = "
                    f"{len(table)}"
                )

                rules = {
                    key: next(iter(bits))
                    for key, bits in table.items()
                }

                for key, bit in list(
                    rules.items()
                )[:16]:

                    print(
                        f"            {key} -> {bit}"
                    )

                found = True

        if not found:

            print(
                "    No exact two-feature rule found."
            )


# ==============================================================================
# CHILD SIGNATURE COMPARISON
# ==============================================================================

def child_signature(
    subset,
    observable,
    m,
):

    counter = Counter()

    for sample in subset:

        value = observable_value(
            sample,
            observable,
        )

        counter[
            mod_value(
                value,
                m,
            )
        ] += 1

    return counter


def compare_children(
    samples,
    level,
):

    parent_m = modulus(level)
    child_m = modulus(level + 1)

    groups = group_by_parent(
        samples,
        level,
    )

    child_groups = group_by_parent(
        samples,
        level + 1,
    )

    print()
    print("=" * 90)
    print(
        f"CHILD SIGNATURE DIFFERENCE"
    )
    print(
        f"mod {parent_m} -> mod {child_m}"
    )
    print("=" * 90)

    for parent in sorted(groups):

        child_a = parent
        child_b = parent + parent_m

        A = child_groups.get(
            child_a,
            [],
        )

        B = child_groups.get(
            child_b,
            [],
        )

        if not A or not B:
            continue

        print()
        print(
            f"PARENT {parent} mod {parent_m}"
        )

        for observable in OBSERVABLES:

            sig_A = child_signature(
                A,
                observable,
                child_m,
            )

            sig_B = child_signature(
                B,
                observable,
                child_m,
            )

            # Compare normalized probability distributions.
            total_A = sum(sig_A.values())
            total_B = sum(sig_B.values())

            keys = set(
                sig_A
            ) | set(
                sig_B
            )

            distance = Fraction(
                0,
                1,
            )

            for key in keys:

                a = Fraction(
                    sig_A.get(key, 0),
                    total_A,
                )

                b = Fraction(
                    sig_B.get(key, 0),
                    total_B,
                )

                distance += abs(
                    a - b
                )

            if distance != 0:

                print(
                    f"    {observable:<7} "
                    f"L1={float(distance):.6f}"
                )


# ==============================================================================
# PARENT MEMORY TEST
# ==============================================================================

def parent_memory_test(
    samples,
    level,
):

    """
    Compare the same feature across different parents.

    If a rule is truly inherited, we expect the child-bit mapping
    to be conditioned by the parent branch.

    We therefore look for:

        feature -> bit

    mappings that are different for different parent residues.
    """

    parent_m = modulus(level)

    groups = group_by_parent(
        samples,
        level,
    )

    print()
    print("=" * 90)
    print(
        "PARENT-MEMORY TEST"
    )
    print(
        f"Parent modulus = {parent_m}"
    )
    print("=" * 90)

    for observable in OBSERVABLES:

        parent_rules = {}

        for parent, subset in groups.items():

            table = defaultdict(set)

            for sample in subset:

                value = observable_value(
                    sample,
                    observable,
                )

                key = mod_value(
                    value,
                    2,
                )

                bit = next_branch_bit(
                    sample.n,
                    level,
                )

                table[key].add(bit)

            if not table:
                continue

            if all(
                len(bits) == 1
                for bits in table.values()
            ):

                rule = tuple(
                    sorted(
                        (
                            key,
                            next(iter(bits)),
                        )
                        for key, bits in table.items()
                    )
                )

                parent_rules[parent] = rule

        if len(set(parent_rules.values())) > 1:

            print()
            print(
                f"FEATURE {observable}"
            )

            for parent, rule in sorted(
                parent_rules.items()
            ):

                print(
                    f"    parent {parent}: "
                    f"{rule}"
                )


# ==============================================================================
# EXAMPLE TRAJECTORIES
# ==============================================================================

def show_trajectories(
    samples,
):

    print()
    print("=" * 90)
    print(
        "EXACT BRANCH TRAJECTORIES"
    )
    print("=" * 90)

    # Pick first sample for each first-level branch
    # and then show its complete binary residue path.

    seen = set()

    for sample in samples:

        parent4 = sample.n % 4

        if parent4 in seen:
            continue

        seen.add(parent4)

        print()
        print(
            f"n = {sample.n}"
        )

        print(
            f"x = {sample.x}, "
            f"y = {sample.y}, "
            f"r = {sample.r}"
        )

        for level in range(
            MIN_LEVEL,
            MAX_LEVEL + 1,
        ):

            m = modulus(level)

            residue = sample.n % m

            print(
                f"    mod {m:<3} "
                f"residue={residue:<3}"
            )

            if level < MAX_LEVEL:

                bit = next_branch_bit(
                    sample.n,
                    level,
                )

                print(
                    f"        next bit = {bit}"
                )

        print()


# ==============================================================================
# SUMMARY
# ==============================================================================

def summary(
    samples,
):

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 572 SUMMARY"
    )
    print("=" * 90)

    print(
        f"Samples = {len(samples)}"
    )

    print()

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_m = modulus(level)
        child_m = modulus(level + 1)

        groups = group_by_parent(
            samples,
            level,
        )

        both = 0
        only_zero = 0
        only_one = 0

        for subset in groups.values():

            bits = {
                next_branch_bit(
                    s.n,
                    level,
                )
                for s in subset
            }

            if bits == {0, 1}:
                both += 1

            elif bits == {0}:
                only_zero += 1

            elif bits == {1}:
                only_one += 1

        print(
            f"mod {parent_m:<3} -> "
            f"mod {child_m:<3}: "
            f"both={both}, "
            f"only0={only_zero}, "
            f"only1={only_one}"
        )

    print()
    print(
        "INTERPRETATION:"
    )

    print()
    print(
        "The strongest evidence would be an EXACT RULE in which"
    )

    print(
        "the next branch bit is determined by x/y/P/Q data."
    )

    print()
    print(
        "Even stronger would be a rule whose form changes when the"
    )

    print(
        "parent residue changes."
    )

    print()
    print(
        "That would support:"
    )

    print()
    print(
        "    parent rule"
    )

    print(
        "        -> child-bit rule"
    )

    print(
        "            -> next parent rule"
    )

    print()
    print(
        "which is the hierarchical mechanism being investigated."
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 572 START")
    print("=" * 90)

    print()
    print(
        "HIERARCHICAL MODULAR BRANCH-PREDICTION"
    )

    print()

    # --------------------------------------------------------------------------
    # Generate.
    # --------------------------------------------------------------------------

    print(
        "[1] GENERATING DATA"
    )

    print(
        "-" * 90
    )

    samples = generate_samples()

    print(
        f"Generated {len(samples)} samples."
    )

    if not samples:

        print(
            "No samples."
        )

        return

    # --------------------------------------------------------------------------
    # Parent bits.
    # --------------------------------------------------------------------------

    print(
        "[2] PARENT / CHILD BIT COUNTS"
    )

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        report_parent_bits(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Mod-2 tests.
    # --------------------------------------------------------------------------

    print(
        "[3] MOD-2 SINGLE FEATURE TEST"
    )

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        observable_mod2_partition(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Modular feature rules.
    # --------------------------------------------------------------------------

    print(
        "[4] SINGLE-FEATURE MODULAR RULE SEARCH"
    )

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        search_single_feature_rules(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Two-feature rules.
    # --------------------------------------------------------------------------

    print(
        "[5] TWO-FEATURE RULE SEARCH"
    )

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        search_two_feature_rules(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Child distributions.
    # --------------------------------------------------------------------------

    print(
        "[6] CHILD SIGNATURE DIFFERENCES"
    )

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        compare_children(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Parent memory.
    # --------------------------------------------------------------------------

    print(
        "[7] PARENT MEMORY"
    )

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_memory_test(
            samples,
            level,
        )

    # --------------------------------------------------------------------------
    # Trajectories.
    # --------------------------------------------------------------------------

    print(
        "[8] EXAMPLE TRAJECTORIES"
    )

    show_trajectories(
        samples,
    )

    # --------------------------------------------------------------------------
    # Final summary.
    # --------------------------------------------------------------------------

    summary(
        samples,
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 572 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
