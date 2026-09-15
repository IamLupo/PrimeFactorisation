#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 575
# ==============================================================================
# HIERARCHICAL MODULAR RULE RECURSION
#
# Main hypothesis
# ---------------
#
# Let
#
#     r_k = n mod 2^k
#
# and let the next binary branch bit be
#
#     b_k = floor(n / 2^k) mod 2.
#
# Thus
#
#     n mod 2^(k+1) = r_k + b_k * 2^k.
#
# We test whether, conditioned on the parent residue r_k,
#
#     b_k = f_(r_k)(x mod 2^(k-1)).
#
# Expected resolution:
#
#     mod 4   -> x mod 2
#     mod 8   -> x mod 4
#     mod 16  -> x mod 8
#     mod 32  -> x mod 16
#     mod 64  -> x mod 32
#     mod 128 -> x mod 64
#
# We also test:
#
#     gap = P-Q
#
# and:
#
#     y
#
# using the same 2020-style crosswalk used in Experiments 571-574:
#
#     P = (2y - 2x + 3) / r
#     Q = (2y + 2x - 3) / r
#
# with
#
#     y = (p+q)/(4r)
#     x = (q-p+6)/(4r)
#
# for r in {1,2,4}.
#
# No floating point is used.
#
# ==============================================================================

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_LEVEL = 2       # 2^2 = 4
MAX_LEVEL = 8       # 2^8 = 256

SHOW_EXACT_RULES = True
SHOW_EXAMPLE_PATHS = True

MAX_RULE_ENTRIES = 32
MAX_EXAMPLE_PATHS = 12


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
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> bytearray:
    """
    sieve[n] == 1 exactly when n is prime.
    """

    if limit < 2:
        return bytearray(limit + 1)

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = isqrt(limit)

    for p in range(2, root + 1):

        if not sieve[p]:
            continue

        start = p * p
        count = (
            (limit - start) // p
        ) + 1

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
    """
    Generate unique odd semiprimes

        n = p*q

    with p <= q and n <= max_n.
    """

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
# 2020 x,y CROSSWALK
# ==============================================================================

def xy_from_pq(
    p: int,
    q: int,
):
    """
    Reconstruct x,y using:

        y = (p+q)/(4r)
        x = (q-p+6)/(4r)

    for r in {1,2,4}.
    """

    result = []

    for r in (
        1,
        2,
        4,
    ):

        denominator = 4 * r

        y_num = p + q
        x_num = q - p + 6

        if y_num % denominator != 0:
            continue

        if x_num % denominator != 0:
            continue

        y = y_num // denominator
        x = x_num // denominator

        result.append(
            (
                r,
                x,
                y,
            )
        )

    return result


# ==============================================================================
# ROOT OBSERVABLES
# ==============================================================================

def make_sample(
    n: int,
    p: int,
    q: int,
    r: int,
    x: int,
    y: int,
):
    """
    Construct exact root observables.
    """

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

    H_P = (
        P
        - 2 * Q
        - Fraction(1, 2)
    )

    C_P = (
        2 * P
        - 4 * Q
        - 1
    )

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

    print(
        "[1] Building prime sieve..."
    )

    sieve = prime_sieve(
        MAX_N
    )

    print(
        "[2] Generating semiprimes..."
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        "    semiprimes =",
        len(semiprimes),
    )

    samples = []

    for n, p, q in semiprimes:

        candidates = xy_from_pq(
            p,
            q,
        )

        for r, x, y in candidates:

            samples.append(
                make_sample(
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

def modulus(level: int) -> int:
    return 1 << level


def parent_residue(
    sample: Sample,
    level: int,
) -> int:
    return sample.n % modulus(level)


def next_bit(
    sample: Sample,
    level: int,
) -> int:
    return (
        sample.n >> level
    ) & 1


def group_by_parent(
    samples,
    level: int,
):
    groups = defaultdict(list)

    for sample in samples:

        parent = parent_residue(
            sample,
            level,
        )

        groups[parent].append(
            sample
        )

    return groups


# ==============================================================================
# FEATURE FUNCTIONS
# ==============================================================================

def x_feature(
    sample: Sample,
    feature_modulus: int,
):
    return sample.x % feature_modulus


def y_feature(
    sample: Sample,
    feature_modulus: int,
):
    return sample.y % feature_modulus


def gap_feature(
    sample: Sample,
    feature_modulus: int,
):
    value = sample.gap

    if value.denominator == 1:
        return value.numerator % feature_modulus

    return (
        value.numerator % feature_modulus,
        value.denominator % feature_modulus,
    )


# ==============================================================================
# EXACT RULE BUILDER
# ==============================================================================

def build_rule(
    subset,
    feature_function,
    feature_modulus: int,
    parent_level: int,
):
    """
    Return an exact mapping

        feature -> next bit

    or None when some feature state produces both bits.
    """

    mapping = defaultdict(set)

    for sample in subset:

        key = feature_function(
            sample,
            feature_modulus,
        )

        bit = next_bit(
            sample,
            parent_level,
        )

        mapping[key].add(
            bit
        )

    if not mapping:
        return None

    for values in mapping.values():

        if len(values) != 1:
            return None

    result = {}

    for key, values in mapping.items():
        result[key] = next(
            iter(values)
        )

    return result


# ==============================================================================
# RULE STATISTICS
# ==============================================================================

def rule_stats(rule):

    if rule is None:
        return None

    zero = 0
    one = 0

    for bit in rule.values():

        if bit == 0:
            zero += 1
        else:
            one += 1

    return {
        "states": len(rule),
        "zero": zero,
        "one": one,
    }


# ==============================================================================
# RULE COMPRESSION
# ==============================================================================

def compress_rule(
    rule,
    full_modulus: int,
):
    """
    Try smaller power-of-two moduli.

    Example:

        x mod 16

    might actually only require

        x mod 4.
    """

    if rule is None:
        return None

    candidate = 2

    while candidate <= full_modulus:

        reduced = {}
        valid = True

        for residue, bit in rule.items():

            if not isinstance(
                residue,
                int,
            ):
                valid = False
                break

            key = residue % candidate

            if key in reduced:

                if reduced[key] != bit:
                    valid = False
                    break

            else:

                reduced[key] = bit

        if valid:
            return (
                candidate,
                reduced,
            )

        candidate *= 2

    return None


# ==============================================================================
# BRANCH STRUCTURE
# ==============================================================================

def analyze_branch_structure(
    samples,
):

    print()
    print("=" * 90)
    print(
        "BRANCH STRUCTURE"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        child_modulus = modulus(
            level + 1
        )

        groups = group_by_parent(
            samples,
            level,
        )

        total_parents = 0
        both_children = 0

        print()
        print(
            "MOD "
            f"{parent_modulus}"
            " -> MOD "
            f"{child_modulus}"
        )

        for parent in sorted(groups):

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if not bits:
                continue

            total_parents += 1

            if bits == {0, 1}:
                both_children += 1

            print(
                f"    parent={parent:<4} "
                f"count={len(subset):<6} "
                f"children={sorted(bits)}"
            )

        print(
            "    total parents =",
            total_parents,
        )

        print(
            "    both children =",
            both_children,
        )


# ==============================================================================
# EXPECTED x RULE
# ==============================================================================

def analyze_expected_x_rule(
    samples,
):
    """
    Main test:

        parent n mod 2^k

            -> x mod 2^(k-1)

            -> next bit
    """

    print()
    print("=" * 90)
    print(
        "EXPECTED x-RESOLUTION TEST"
    )
    print("=" * 90)

    results = []

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        expected_x_modulus = modulus(
            level - 1
        )

        groups = group_by_parent(
            samples,
            level,
        )

        exact_count = 0
        two_child_count = 0

        print()
        print(
            "MOD "
            f"{parent_modulus}"
            " -> MOD "
            f"{parent_modulus * 2}"
        )

        print(
            "    expected feature = x mod "
            f"{expected_x_modulus}"
        )

        for parent in sorted(groups):

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if bits != {0, 1}:
                continue

            two_child_count += 1

            rule = build_rule(
                subset,
                x_feature,
                expected_x_modulus,
                level,
            )

            if rule is None:

                print(
                    f"    parent={parent:<4} "
                    "FAIL"
                )

                results.append(
                    (
                        level,
                        parent,
                        False,
                        None,
                    )
                )

                continue

            exact_count += 1

            stats = rule_stats(
                rule
            )

            compressed = compress_rule(
                rule,
                expected_x_modulus,
            )

            if compressed is not None:

                compressed_modulus, _ = (
                    compressed
                )

                print(
                    f"    parent={parent:<4} "
                    "EXACT "
                    f"states={stats['states']:<3} "
                    "compressed=x mod "
                    f"{compressed_modulus}"
                )

            else:

                print(
                    f"    parent={parent:<4} "
                    "EXACT "
                    f"states={stats['states']}"
                )

            results.append(
                (
                    level,
                    parent,
                    True,
                    rule,
                )
            )

        print(
            "    EXACT = "
            f"{exact_count}/{two_child_count}"
        )

    return results


# ==============================================================================
# MINIMUM x RESOLUTION
# ==============================================================================

def minimum_x_resolution(
    samples,
):
    """
    Find the smallest modulus of x that exactly predicts the next bit.
    """

    print()
    print("=" * 90)
    print(
        "MINIMUM x-RESOLUTION"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        groups = group_by_parent(
            samples,
            level,
        )

        print()
        print(
            f"PARENT MOD {parent_modulus}"
        )

        for parent in sorted(groups):

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if bits != {0, 1}:
                continue

            found = None

            feature_level = 1

            while feature_level <= (
                MAX_LEVEL + 2
            ):

                feature_modulus = modulus(
                    feature_level
                )

                rule = build_rule(
                    subset,
                    x_feature,
                    feature_modulus,
                    level,
                )

                if rule is not None:

                    found = (
                        feature_modulus,
                        feature_level,
                        rule,
                    )

                    break

                feature_level += 1

            if found is None:

                print(
                    f"    parent={parent:<4} "
                    "NO EXACT RULE"
                )

                continue

            feature_modulus, feature_level, _ = (
                found
            )

            expected = modulus(
                level - 1
            )

            if feature_modulus == expected:
                status = "MATCH"
            elif feature_modulus < expected:
                status = "LESS"
            else:
                status = "MORE"

            print(
                f"    parent={parent:<4} "
                f"minimum=x mod {feature_modulus:<4} "
                f"expected={expected:<4} "
                f"{status}"
            )


# ==============================================================================
# GAP TEST
# ==============================================================================

def analyze_gap(
    samples,
):

    print()
    print("=" * 90)
    print(
        "GAP CROSS-CHECK"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        gap_modulus = modulus(
            level - 1
        )

        groups = group_by_parent(
            samples,
            level,
        )

        total = 0
        exact = 0

        for parent in sorted(groups):

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if bits != {0, 1}:
                continue

            total += 1

            rule = build_rule(
                subset,
                gap_feature,
                gap_modulus,
                level,
            )

            if rule is not None:
                exact += 1

        print(
            f"MOD {parent_modulus:<4} -> "
            f"MOD {parent_modulus * 2:<4} "
            f"| gap mod {gap_modulus:<4} "
            f"| EXACT {exact}/{total}"
        )


# ==============================================================================
# Y TEST
# ==============================================================================

def analyze_y(
    samples,
):

    print()
    print("=" * 90)
    print(
        "y CROSS-CHECK"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        y_modulus = modulus(
            level - 1
        )

        groups = group_by_parent(
            samples,
            level,
        )

        total = 0
        exact = 0

        for parent in sorted(groups):

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if bits != {0, 1}:
                continue

            total += 1

            rule = build_rule(
                subset,
                y_feature,
                y_modulus,
                level,
            )

            if rule is not None:
                exact += 1

        print(
            f"MOD {parent_modulus:<4} -> "
            f"MOD {parent_modulus * 2:<4} "
            f"| y mod {y_modulus:<4} "
            f"| EXACT {exact}/{total}"
        )


# ==============================================================================
# PARENT RULE TABLES
# ==============================================================================

def print_compact_rules(
    samples,
):

    print()
    print("=" * 90)
    print(
        "COMPACT EXACT x-RULES"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        feature_modulus = modulus(
            level - 1
        )

        groups = group_by_parent(
            samples,
            level,
        )

        print()
        print(
            f"MOD {parent_modulus} "
            f"-> MOD {parent_modulus * 2}"
        )

        for parent in sorted(groups):

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if bits != {0, 1}:
                continue

            rule = build_rule(
                subset,
                x_feature,
                feature_modulus,
                level,
            )

            if rule is None:
                continue

            compressed = compress_rule(
                rule,
                feature_modulus,
            )

            if compressed is not None:

                compressed_modulus, reduced_rule = (
                    compressed
                )

                ones = []

                for residue, bit in sorted(
                    reduced_rule.items()
                ):

                    if bit == 1:
                        ones.append(
                            residue
                        )

                print(
                    f"    parent={parent:<4} "
                    f"x_mod={compressed_modulus:<4} "
                    f"bit1={ones}"
                )

            else:

                print(
                    f"    parent={parent:<4} "
                    f"x_mod={feature_modulus:<4}"
                )

                if SHOW_EXACT_RULES:

                    entries = sorted(
                        rule.items(),
                        key=lambda item: (
                            str(item[0])
                        ),
                    )

                    limited = entries[
                        :MAX_RULE_ENTRIES
                    ]

                    print(
                        "        rule =",
                        limited,
                    )

                    if len(entries) > MAX_RULE_ENTRIES:

                        print(
                            "        ... "
                            f"{len(entries) - MAX_RULE_ENTRIES} "
                            "more states"
                        )


# ==============================================================================
# RULE STABILITY
# ==============================================================================

def analyze_rule_stability(
    samples,
):

    print()
    print("=" * 90)
    print(
        "RULE STABILITY"
    )
    print("=" * 90)

    for level in range(
        MIN_LEVEL,
        MAX_LEVEL,
    ):

        parent_modulus = modulus(
            level
        )

        feature_modulus = modulus(
            level - 1
        )

        groups = group_by_parent(
            samples,
            level,
        )

        signatures = Counter()

        for parent in groups:

            subset = groups[parent]

            bits = set()

            for sample in subset:
                bits.add(
                    next_bit(
                        sample,
                        level,
                    )
                )

            if bits != {0, 1}:
                continue

            rule = build_rule(
                subset,
                x_feature,
                feature_modulus,
                level,
            )

            if rule is None:

                signatures[
                    "FAIL"
                ] += 1

                continue

            ones = tuple(
                sorted(
                    residue
                    for residue, bit in rule.items()
                    if bit == 1
                )
            )

            signature = (
                len(rule),
                ones,
            )

            signatures[
                signature
            ] += 1

        print()
        print(
            f"MOD {parent_modulus} "
            f"(x mod {feature_modulus})"
        )

        if not signatures:

            print(
                "    no two-child nodes"
            )

            continue

        for signature, count in (
            signatures.most_common(8)
        ):

            if signature == "FAIL":

                print(
                    f"    FAIL nodes = {count}"
                )

                continue

            states, ones = signature

            print(
                f"    count={count:<3} "
                f"states={states:<3} "
                f"bit1={list(ones)}"
            )


# ==============================================================================
# SAMPLE PATHS
# ==============================================================================

def print_example_paths(
    samples,
):

    if not SHOW_EXAMPLE_PATHS:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLE 2-ADIC PATHS"
    )
    print("=" * 90)

    seen_n = set()

    for sample in samples:

        if sample.n in seen_n:
            continue

        seen_n.add(
            sample.n
        )

        print()
        print(
            f"n={sample.n} "
            f"p={sample.p} "
            f"q={sample.q} "
            f"x={sample.x} "
            f"y={sample.y} "
            f"r={sample.r}"
        )

        for level in range(
            MIN_LEVEL,
            MAX_LEVEL + 1,
        ):

            m = modulus(level)
            residue = sample.n % m

            line = (
                f"    mod {m:<4} "
                f"residue={residue:<4}"
            )

            if level < MAX_LEVEL:

                bit = next_bit(
                    sample,
                    level,
                )

                line += (
                    f" next_bit={bit}"
                )

            print(line)

        if len(seen_n) >= MAX_EXAMPLE_PATHS:
            break


# ==============================================================================
# FINAL AUDIT
# ==============================================================================

def final_audit(
    samples,
):

    print()
    print("=" * 90)
    print(
        "FINAL AUDIT"
    )
    print("=" * 90)

    print(
        "samples =",
        len(samples),
    )

    print()
    print(
        "Target:"
    )

    print(
        "    r_k = n mod 2^k"
    )

    print(
        "    b_k = floor(n / 2^k) mod 2"
    )

    print(
        "    b_k = f_(r_k)(x mod 2^(k-1))"
    )

    print()
    print(
        "A strong result requires:"
    )

    print(
        "    - both children at many consecutive levels;"
    )

    print(
        "    - exact x-rule on every observed parent;"
    )

    print(
        "    - the required x resolution scaling with the level;"
    )

    print(
        "    - a compact or repeating rule family;"
    )

    print(
        "    - persistence when the range of n is increased."
    )

    print()
    print(
        "An exact lookup table by itself is not sufficient to"
    )

    print(
        "establish a recursive law."
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 575 START"
    )
    print("=" * 90)

    print()
    print(
        "HIERARCHICAL MODULAR RULE RECURSION"
    )

    print()
    print(
        "Testing:"
    )

    print(
        "    b_k = f_(r_k)(x mod 2^(k-1))"
    )

    # --------------------------------------------------------------------------
    # Generate
    # --------------------------------------------------------------------------

    samples = generate_samples()

    print()
    print(
        "Generated samples =",
        len(samples),
    )

    if not samples:

        print(
            "No samples generated."
        )

        return

    # --------------------------------------------------------------------------
    # Structure
    # --------------------------------------------------------------------------

    analyze_branch_structure(
        samples
    )

    # --------------------------------------------------------------------------
    # Main hypothesis
    # --------------------------------------------------------------------------

    analyze_expected_x_rule(
        samples
    )

    # --------------------------------------------------------------------------
    # Minimum x information
    # --------------------------------------------------------------------------

    minimum_x_resolution(
        samples
    )

    # --------------------------------------------------------------------------
    # Alternative observables
    # --------------------------------------------------------------------------

    analyze_gap(
        samples
    )

    analyze_y(
        samples
    )

    # --------------------------------------------------------------------------
    # Compact rule family
    # --------------------------------------------------------------------------

    analyze_rule_stability(
        samples
    )

    if SHOW_EXACT_RULES:

        print_compact_rules(
            samples
        )

    # --------------------------------------------------------------------------
    # Paths
    # --------------------------------------------------------------------------

    print_example_paths(
        samples
    )

    # --------------------------------------------------------------------------
    # Final
    # --------------------------------------------------------------------------

    final_audit(
        samples
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 575 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()