#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 114
N-ONLY ACCESSIBILITY BARRIER
RESIDUE FINGERPRINTS VS PAPER QUOTIENTS
NO FACTOR-PAIR SEARCH IN INVERSE STAGE
NO CSV
NO SKLEARN
==============================================================================

PURPOSE
-------

Experiments 46R / 102R / 103R established that the paper quotient family has
an exact algebraic description in the hidden symmetric variable

    S = p + q

and that the Newton/P0 bookkeeping is internally consistent.

The question now changes:

    Can a bounded, genuinely N-only residue fingerprint determine a paper
    quotient without recovering S, p, q, or an equivalent divisor-sum object?

This experiment does NOT try to "prove impossibility". It measures a concrete
barrier:

    N mod m  ->  Q_(k,l) mod m

and then

    (N mod m1, ..., N mod mr)
        -> Q_(k,l) mod L

for progressively stronger fingerprints.

A useful negative result is:

    many targets have identical N-only fingerprints but different oracle Q.

That means the corresponding bounded residue information cannot determine the
detector quotient.

The experiment also checks the same phenomenon for S itself.

IMPORTANT
---------
* p,q are used only to build oracle targets and validate Q.
* The inverse fingerprint tests use only N residues.
* No factor-pair search is performed.
* No CSV.
* No sklearn.
==============================================================================

"""

from __future__ import annotations

import math
import statistics
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import sympy as sp


# =============================================================================
# PARAMETERS
# =============================================================================

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

TARGETS = 6000
TRAIN_FRACTION = 0.75

# Small odd moduli. Keep these deliberately modest: the point is to test
# bounded residue information, not to encode N itself.
MODULI = (
    5, 7, 11, 13, 17, 19,
    23, 29, 31, 37,
    41, 43, 47,
)

# Fingerprint stages.
FINGERPRINT_STAGES = (
    (5,),
    (5, 7),
    (5, 7, 11),
    (5, 7, 11, 13),
    (5, 7, 11, 13, 17),
    (5, 7, 11, 13, 17, 19),
    (5, 7, 11, 13, 17, 19, 23),
    (5, 7, 11, 13, 17, 19, 23, 29),
)

# Quotients to test.
DETECTORS = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

# We test quotient residues modulo these output moduli.
OUTPUT_MODULI = (
    5, 7, 11, 13, 17, 19, 23, 29, 31,
)

# Deterministic seed.
SEED = 114


# =============================================================================
# SYMBOLIC DEFINITIONS
# =============================================================================

N, S = sp.symbols("N S")


# =============================================================================
# TARGET DATA
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int

    @property
    def N(self) -> int:
        return self.p * self.q

    @property
    def S(self) -> int:
        return self.p + self.q


# =============================================================================
# PRIME POPULATION
# =============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    """
    Return primes in [lo, hi] using a segmented sieve.

    The experiment only needs prime generation. No factor-pair search is
    performed on the resulting semiprimes.
    """
    if hi < lo:
        return []

    base_limit = int(math.isqrt(hi)) + 1

    base = bytearray(b"\x01") * (base_limit + 1)
    base[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(base_limit)) + 1):
        if base[p]:
            start = p * p
            base[start:base_limit + 1:p] = (
                b"\x00"
                * (((base_limit - start) // p) + 1)
            )

    base_primes = [
        p for p in range(2, base_limit + 1)
        if base[p]
    ]

    segment = max(100_000, hi - lo + 1)

    out: List[int] = []

    for left in range(lo, hi + 1, segment):
        right = min(left + segment - 1, hi)
        arr = bytearray(b"\x01") * (right - left + 1)

        for p in base_primes:
            first = max(
                p * p,
                ((left + p - 1) // p) * p
            )

            if first > right:
                continue

            arr[first - left:right - left + 1:p] = (
                b"\x00"
                * (((right - first) // p) + 1)
            )

        if left == 0:
            arr[0] = 0
            if right >= 1:
                arr[1] = 0
        elif left == 1:
            arr[0] = 0

        out.extend(
            left + i
            for i, flag in enumerate(arr)
            if flag
        )

    return out


def build_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:
    """
    Deterministic semiprime target construction.

    We deliberately generate a broad distribution of p,q rather than making
    q close to p. This is useful because a residue-fingerprint barrier should
    not depend on near-square semiprimes.
    """
    if len(primes) < 2:
        raise ValueError("Need at least two primes.")

    targets: List[Target] = []

    # Deterministic modular walk through the prime population.
    npr = len(primes)

    state = SEED

    used = set()

    while len(targets) < count:
        state = (state * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)

        i = state % npr

        state = (state * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)

        # Use a relatively independent second index.
        j = (state % npr)

        if i == j:
            continue

        p = int(primes[i])
        q = int(primes[j])

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in used:
            continue

        used.add(key)
        targets.append(Target(p, q))

    return targets


# =============================================================================
# PAPER DETECTOR
# =============================================================================

def raw_detector(k: int, ell: int, p: int, q: int) -> int:
    """
    Direct paper numerator for semiprime p*q.
    """
    return (
        p**k * (q + 1)**ell
        + q**k * (p + 1)**ell
        - p**ell * (q + 1)**k
        - q**ell * (p + 1)**k
    )


def quotient_value(k: int, ell: int, p: int, q: int) -> int:
    """
    Exact oracle quotient

        Q_(k,l) = F_(k,l)/(p+q+1).

    This uses p,q only to construct the oracle value.
    """
    F = raw_detector(k, ell, p, q)
    divisor = p + q + 1

    if F % divisor != 0:
        raise ArithmeticError(
            f"Non-exact detector quotient for ({k},{ell}), "
            f"p={p}, q={q}"
        )

    return F // divisor


# =============================================================================
# RESIDUE FINGERPRINTS
# =============================================================================

def fingerprint(N_value: int, moduli: Sequence[int]) -> Tuple[int, ...]:
    return tuple(
        N_value % m
        for m in moduli
    )


def output_residue(
    value: int,
    modulus: int,
) -> int:
    return value % modulus


# =============================================================================
# COLLISION STATISTICS
# =============================================================================

@dataclass
class CollisionSummary:
    groups: int
    repeated_groups: int
    repeated_targets: int
    ambiguous_groups: int
    deterministic_groups: int
    max_group_size: int
    total_targets: int


def summarize_groups(
    groups: Dict[Tuple[int, ...], List[int]],
) -> CollisionSummary:
    sizes = [len(v) for v in groups.values()]

    repeated = [
        size for size in sizes
        if size > 1
    ]

    return CollisionSummary(
        groups=len(groups),
        repeated_groups=len(repeated),
        repeated_targets=sum(repeated),
        ambiguous_groups=len(repeated),
        deterministic_groups=len(groups) - len(repeated),
        max_group_size=max(sizes) if sizes else 0,
        total_targets=sum(sizes),
    )


# =============================================================================
# ACCESSIBILITY TEST
# =============================================================================

def build_feature_groups(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> Dict[Tuple[int, ...], List[int]]:
    groups: Dict[Tuple[int, ...], List[int]] = defaultdict(list)

    for idx, target in enumerate(targets):
        fp = fingerprint(target.N, moduli)
        groups[fp].append(idx)

    return dict(groups)


def detect_output_ambiguity(
    targets: Sequence[Target],
    values: Sequence[int],
    feature_groups: Dict[Tuple[int, ...], List[int]],
    output_modulus: int,
) -> Tuple[int, int, int]:
    """
    Returns:
        ambiguous_groups
        deterministic_groups
        differing_pairs

    A feature group is ambiguous if the same N-only fingerprint maps to more
    than one output residue.
    """
    ambiguous = 0
    deterministic = 0
    differing_pairs = 0

    for indices in feature_groups.values():
        residues = {
            values[i] % output_modulus
            for i in indices
        }

        if len(residues) <= 1:
            deterministic += 1
        else:
            ambiguous += 1

            # Count unordered residue disagreements.
            residues_list = list(
                values[i] % output_modulus
                for i in indices
            )

            for i in range(len(residues_list)):
                for j in range(i + 1, len(residues_list)):
                    if residues_list[i] != residues_list[j]:
                        differing_pairs += 1

    return ambiguous, deterministic, differing_pairs


def fingerprint_entropy_like(
    targets: Sequence[Target],
    values: Sequence[int],
    moduli: Sequence[int],
    output_modulus: int,
) -> Tuple[float, float]:
    """
    Two simple descriptive quantities:

    ambiguity_rate:
        fraction of fingerprint groups that admit multiple output residues.

    target_ambiguity_rate:
        fraction of targets lying inside an ambiguous fingerprint group.
    """
    groups: Dict[Tuple[int, ...], set] = defaultdict(set)
    counts: Dict[Tuple[int, ...], int] = defaultdict(int)

    for target, value in zip(targets, values):
        fp = fingerprint(target.N, moduli)
        groups[fp].add(value % output_modulus)
        counts[fp] += 1

    ambiguous_groups = sum(
        1 for residues in groups.values()
        if len(residues) > 1
    )

    ambiguous_targets = sum(
        counts[fp]
        for fp, residues in groups.items()
        if len(residues) > 1
    )

    total_groups = len(groups)
    total_targets = len(targets)

    return (
        ambiguous_groups / total_groups
        if total_groups else 0.0,
        ambiguous_targets / total_targets
        if total_targets else 0.0,
    )


# =============================================================================
# EXACT N-ONLY LOW-DEGREE TEST
# =============================================================================

def low_degree_polynomial_fit_failure(
    train_targets: Sequence[Target],
    train_values: Sequence[int],
    test_targets: Sequence[Target],
    test_values: Sequence[int],
    degree: int,
) -> Tuple[int, int]:
    """
    Fit an exact polynomial Q(N) over QQ using training points.

    This is deliberately a very restricted N-only model.

    Returns:
        train_failures, test_failures
    """
    if len(train_targets) <= degree:
        return (
            len(train_targets),
            len(test_targets),
        )

    x = sp.symbols("x")

    pairs = [
        (sp.Integer(t.N), sp.Integer(v))
        for t, v in zip(train_targets, train_values)
    ]

    poly = sp.interpolate(
        pairs[:degree + 1],
        x,
    )

    poly = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.QQ,
    )

    train_fail = 0
    for t, v in zip(train_targets, train_values):
        got = poly.eval(t.N)
        if got != v:
            train_fail += 1

    test_fail = 0
    for t, v in zip(test_targets, test_values):
        got = poly.eval(t.N)
        if got != v:
            test_fail += 1

    return train_fail, test_fail


# =============================================================================
# SAME-N RESIDUE VS SAME-N FINGERPRINT DIAGNOSTIC
# =============================================================================

def residue_progression(
    targets: Sequence[Target],
    value_fn,
    label: str,
) -> None:
    print()
    print(f"  {label}")
    print("  " + "-" * 72)

    values = [
        value_fn(t)
        for t in targets
    ]

    for mods in FINGERPRINT_STAGES:
        groups = build_feature_groups(
            targets,
            mods,
        )

        summary = summarize_groups(groups)

        amb_groups, amb_targets = fingerprint_entropy_like(
            targets,
            values,
            mods,
            5,
        )

        print(
            f"  mods={mods!s:<42} "
            f"groups={summary.groups:4d} "
            f"repeated={summary.repeated_groups:4d} "
            f"ambig_groups={amb_groups:7.3f} "
            f"ambig_targets={amb_targets:7.3f}"
        )


# =============================================================================
# CROSS-OUTPUT-COMPATIBILITY TEST
# =============================================================================

def cross_output_test(
    targets: Sequence[Target],
    detector: Tuple[int, int],
) -> None:
    """
    For one detector, test whether increasingly rich N fingerprints determine
    Q modulo several small output moduli.

    This is intentionally stronger than testing just one output modulus.
    """
    k, ell = detector

    values = [
        quotient_value(k, ell, t.p, t.q)
        for t in targets
    ]

    print()
    print(f"DETECTOR Q_({k},{ell})")
    print("-" * 78)

    for mods in FINGERPRINT_STAGES:
        groups = build_feature_groups(
            targets,
            mods,
        )

        summary = summarize_groups(groups)

        worst_ambiguous = 0
        worst_target_ambiguity = 0.0

        for out_mod in OUTPUT_MODULI:
            amb_groups, amb_targets, _ = detect_output_ambiguity(
                targets,
                values,
                groups,
                out_mod,
            )

            worst_ambiguous = max(
                worst_ambiguous,
                amb_groups,
            )

            if summary.groups:
                target_rate = 0.0

                ambiguous_target_count = 0

                for fp, indices in groups.items():
                    residues = {
                        values[i] % out_mod
                        for i in indices
                    }

                    if len(residues) > 1:
                        ambiguous_target_count += len(indices)

                target_rate = (
                    ambiguous_target_count / len(targets)
                )

                worst_target_ambiguity = max(
                    worst_target_ambiguity,
                    target_rate,
                )

        print(
            f"mods={mods!s:<42} "
            f"groups={summary.groups:4d} "
            f"max_group={summary.max_group_size:3d} "
            f"worst_amb_groups={worst_ambiguous:4d} "
            f"worst_amb_targets={worst_target_ambiguity:7.3f}"
        )


# =============================================================================
# HOLDOUT TEST FOR FINGERPRINT DETERMINISM
# =============================================================================

def holdout_fingerprint_test(
    targets: Sequence[Target],
    detector: Tuple[int, int],
) -> None:
    """
    Train/test are used only to determine whether a fingerprint seen in
    training has a unique oracle residue.

    This avoids any interpolation and avoids using p,q as features.
    """
    k, ell = detector

    values = [
        quotient_value(k, ell, t.p, t.q)
        for t in targets
    ]

    split = int(len(targets) * TRAIN_FRACTION)

    train_targets = targets[:split]
    test_targets = targets[split:]

    train_values = values[:split]
    test_values = values[split:]

    print()
    print(
        f"HOLDOUT FINGERPRINT TEST Q_({k},{ell})"
    )
    print("-" * 78)

    for mods in FINGERPRINT_STAGES:
        lookup: Dict[Tuple[int, ...], set] = defaultdict(set)

        for t, v in zip(train_targets, train_values):
            lookup[
                fingerprint(t.N, mods)
            ].add(v % 11)

        covered = 0
        predictable = 0
        contradictions = 0

        for t, v in zip(test_targets, test_values):
            fp = fingerprint(t.N, mods)

            if fp not in lookup:
                continue

            covered += 1

            possible = lookup[fp]

            if len(possible) == 1:
                predictable += 1

                if next(iter(possible)) != v % 11:
                    contradictions += 1
            else:
                contradictions += 1

        print(
            f"mods={mods!s:<42} "
            f"covered={covered:3d} "
            f"predictable={predictable:3d} "
            f"contradictions={contradictions:3d}"
        )


# =============================================================================
# SAME-FINGERPRINT EXPLICIT COUNTEREXAMPLES
# =============================================================================

def print_counterexamples(
    targets: Sequence[Target],
    detector: Tuple[int, int],
    max_examples: int = 10,
) -> None:
    """
    Print a few explicit collisions where:

        N fingerprints are identical
        but Q residues differ.

    This is the most concrete output of the experiment.
    """
    k, ell = detector

    values = [
        quotient_value(k, ell, t.p, t.q)
        for t in targets
    ]

    print()
    print(
        f"COUNTEREXAMPLES Q_({k},{ell})"
    )
    print("-" * 78)

    shown = 0

    for mods in FINGERPRINT_STAGES:
        groups = build_feature_groups(
            targets,
            mods,
        )

        for fp, indices in groups.items():
            for i in range(len(indices)):
                for j in range(i + 1, len(indices)):
                    a = indices[i]
                    b = indices[j]

                    if values[a] % 11 == values[b] % 11:
                        continue

                    print(
                        f"mods={mods}"
                    )
                    print(
                        f"  fingerprint = {fp}"
                    )
                    print(
                        f"  A: N={targets[a].N} "
                        f"S={targets[a].S} "
                        f"Q mod 11={values[a] % 11}"
                    )
                    print(
                        f"  B: N={targets[b].N} "
                        f"S={targets[b].S} "
                        f"Q mod 11={values[b] % 11}"
                    )

                    shown += 1

                    if shown >= max_examples:
                        return


# =============================================================================
# LOW-DEGREE N-ONLY BASELINE
# =============================================================================

def low_degree_baseline(
    targets: Sequence[Target],
    detector: Tuple[int, int],
) -> None:
    k, ell = detector

    values = [
        quotient_value(k, ell, t.p, t.q)
        for t in targets
    ]

    split = int(len(targets) * TRAIN_FRACTION)

    train_targets = targets[:split]
    test_targets = targets[split:]

    train_values = values[:split]
    test_values = values[split:]

    print()
    print(
        f"LOW-DEGREE N-ONLY BASELINE Q_({k},{ell})"
    )
    print("-" * 78)

    for degree in (0, 1, 2, 3, 4, 5):
        train_fail, test_fail = low_degree_polynomial_fit_failure(
            train_targets,
            train_values,
            test_targets,
            test_values,
            degree,
        )

        print(
            f"degree={degree}: "
            f"train_fail={train_fail}/{len(train_targets)} "
            f"test_fail={test_fail}/{len(test_targets)}"
        )


# =============================================================================
# GLOBAL SUMMARY
# =============================================================================

def global_summary(
    targets: Sequence[Target],
) -> None:
    print()
    print("=" * 78)
    print("GLOBAL N-ONLY ACCESSIBILITY SUMMARY")
    print("=" * 78)

    for detector in DETECTORS:
        k, ell = detector

        values = [
            quotient_value(k, ell, t.p, t.q)
            for t in targets
        ]

        print()
        print(
            f"Q_({k},{ell})"
        )

        for mods in (
            (5,),
            (5, 7, 11),
            (5, 7, 11, 13, 17, 19),
        ):
            groups = build_feature_groups(
                targets,
                mods,
            )

            ambiguous_by_mod11 = 0
            ambiguous_targets = 0

            for fp, indices in groups.items():
                residues = {
                    values[i] % 11
                    for i in indices
                }

                if len(residues) > 1:
                    ambiguous_by_mod11 += 1
                    ambiguous_targets += len(indices)

            print(
                f"  mods={mods}: "
                f"groups={len(groups)} "
                f"ambiguous_groups={ambiguous_by_mod11} "
                f"ambiguous_targets={ambiguous_targets}"
            )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 114")
    print("N-ONLY ACCESSIBILITY BARRIER")
    print("RESIDUE FINGERPRINTS VS PAPER QUOTIENTS")
    print("NO FACTOR-PAIR SEARCH IN INVERSE STAGE")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. PRIME POPULATION
    # -------------------------------------------------------------------------
    print()
    print("1. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    elapsed = time.perf_counter() - start

    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = {elapsed:.6f}s"
    )

    if len(primes) < 1000:
        raise ArithmeticError(
            "Prime population unexpectedly small."
        )

    # -------------------------------------------------------------------------
    # 2. TARGETS
    # -------------------------------------------------------------------------
    print()
    print("2. TARGET POPULATION")
    print("-" * 78)

    targets = build_targets(
        primes,
        TARGETS,
    )

    print(
        f"targets = {len(targets)}"
    )

    if len(targets) != TARGETS:
        raise ArithmeticError(
            "Target generation failed."
        )

    # -------------------------------------------------------------------------
    # 3. BASIC ORACLE VALIDATION
    # -------------------------------------------------------------------------
    print()
    print("3. ORACLE DETECTOR VALIDATION")
    print("-" * 78)

    failures = 0

    for target in targets[: min(100, len(targets))]:
        for k, ell in DETECTORS:
            F = raw_detector(
                k,
                ell,
                target.p,
                target.q,
            )

            d = target.S + 1

            if F % d != 0:
                failures += 1

    print(
        f"oracle quotient failures = "
        f"{failures}"
    )

    if failures:
        raise ArithmeticError(
            "Oracle quotient validation failed."
        )

    print("STATUS = PASS")

    # -------------------------------------------------------------------------
    # 4. N -> S IS NOT LOCALLY RESIDUE-DETERMINISTIC
    # -------------------------------------------------------------------------
    print()
    print("4. N-RESIDUE FINGERPRINT VS S")
    print("-" * 78)

    residue_progression(
        targets,
        lambda t: t.S,
        "hidden trace S=p+q",
    )

    # -------------------------------------------------------------------------
    # 5. N-RESIDUE FINGERPRINT VS EACH DETECTOR
    # -------------------------------------------------------------------------
    print()
    print("5. N-RESIDUE FINGERPRINT VS PAPER QUOTIENTS")
    print("-" * 78)

    for detector in DETECTORS:
        cross_output_test(
            targets,
            detector,
        )

    # -------------------------------------------------------------------------
    # 6. EXPLICIT COUNTEREXAMPLES
    # -------------------------------------------------------------------------
    print()
    print("6. EXPLICIT RESIDUE COLLISIONS")
    print("-" * 78)

    print_counterexamples(
        targets,
        (1, 3),
        max_examples=8,
    )

    # -------------------------------------------------------------------------
    # 7. HOLDOUT ACCESSIBILITY
    # -------------------------------------------------------------------------
    print()
    print("7. STRICT HOLDOUT ACCESSIBILITY")
    print("-" * 78)

    for detector in (
        (1, 3),
        (1, 5),
        (3, 5),
    ):
        holdout_fingerprint_test(
            targets,
            detector,
        )

    # -------------------------------------------------------------------------
    # 8. LOW-DEGREE N-ONLY BASELINE
    # -------------------------------------------------------------------------
    print()
    print("8. LOW-DEGREE N-ONLY BASELINE")
    print("-" * 78)

    for detector in (
        (1, 3),
        (1, 5),
        (3, 5),
    ):
        low_degree_baseline(
            targets,
            detector,
        )

    # -------------------------------------------------------------------------
    # 9. GLOBAL SUMMARY
    # -------------------------------------------------------------------------
    global_summary(targets)

    elapsed = time.perf_counter() - t0

    # -------------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # -------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)
    print()
    print(
        "This experiment deliberately stops deriving more symbolic identities."
    )
    print()
    print(
        "The question is whether bounded N-only residue information can "
        "determine the detector quotient."
    )
    print()
    print(
        "A strong barrier result is:"
    )
    print(
        "    same N fingerprint"
    )
    print(
        "        +"
    )
    print(
        "    different Q residue"
    )
    print(
        "        ->"
    )
    print(
        "    fingerprint is insufficient."
    )
    print()
    print(
        "The strict holdout section checks that this is not merely an "
        "in-sample collision pattern."
    )
    print()
    print(
        "The low-degree polynomial section is only a baseline. It should "
        "not be interpreted as a proof of impossibility."
    )
    print()
    print(
        "The important distinction remains:"
    )
    print()
    print(
        "    Q = R(N) exists mathematically"
    )
    print(
        "    != "
    )
    print(
        "    R(N) is efficiently computable without recovering equivalent "
        "divisor information."
    )
    print()
    print(
        "If the residue barrier is strong across several detector outputs "
        "and moduli, the next experiments should move away from the paper "
        "detector tensor and toward genuinely different N-only information."
    )
    print()
    print(
        f"total runtime = {elapsed:.6f}s"
    )
    print("=" * 78)
    print("EXPERIMENT 114 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            f"\nFATAL: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise

