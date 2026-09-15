#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 89
EXACT FINITE-FIELD RESOLVENT NULL

FOCUS:
    T = s^2 - 3n
      = p^2 - p q + q^2
      = n + (p-q)^2

QUESTION:
    Is the apparent N -> T modular signal at small ell anything more
    than unavoidable finite-field algebra plus the empirical distribution
    of prime residues?

TESTS:
    1. EXACT F_l enumeration
    2. ACTUAL PRIME-PAIR empirical distribution
    3. EMPIRICAL-RESIDUE SYNTHETIC NULL
    4. UNIFORM NONZERO-RESIDUE NULL
    5. OOS lookup using ONLY n mod ell
    6. Exact Bayes-optimal accuracy from the finite-field model
    7. Comparison of actual-vs-theoretical conditional distributions

NO CSV
NO SKLEARN
==============================================================================

IMPORTANT:

For a fixed ell and residue pair (a,b):

    n = a*b mod ell
    T = a^2 - a*b + b^2 mod ell

The exact finite-field enumeration gives the mathematically allowed
conditional distribution P(T | N) before any target data are considered.

The synthetic empirical null samples residue pairs from the observed
prime-residue frequencies modulo ell. This preserves the finite-population
prime-residue structure while destroying any higher-order target structure.

If the observed OOS performance is essentially identical to this null,
then the earlier "signal" is explained by finite-field residue algebra.

==============================================================================
"""

from __future__ import annotations

import time
import math
import random
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 4000
TRAIN_TARGETS = 3000

P_MIN = 2_000_000
P_MAX = 4_200_000

MODULI = [
    5, 7, 11, 13, 17, 19, 23, 31, 37, 61, 67
]

RNG_SEED = 89089

SYNTHETIC_SAMPLES = 100_000

DETAIL_MODULI = {
    5, 7, 11, 13, 17, 19, 23, 31, 37, 61, 67
}


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

    lo = max(2, lo)

    sieve = bytearray(
        b"\x01" * (hi + 1)
    )

    sieve[0:2] = b"\x00\x00"

    for p in range(
        2,
        math.isqrt(hi) + 1,
    ):
        if sieve[p]:
            start = p * p
            sieve[
                start:hi + 1:p
            ] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x
        for x in range(lo, hi + 1)
        if sieve[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    rng = random.Random(
        RNG_SEED
    )

    out: List[Target] = []
    seen = set()

    while len(out) < count:

        p = primes[
            rng.randrange(
                len(primes)
            )
        ]

        q = primes[
            rng.randrange(
                len(primes)
            )
        ]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        out.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return out


# ============================================================================
# FINITE-FIELD MAP
# ============================================================================

def exact_pair_map(
    ell: int,
) -> Dict[int, Counter]:

    """
    Exact mathematical map

        n -> distribution of T

    over all ordered nonzero residue pairs.

    We include both orientations (a,b) because the unordered semiprime
    pair is symmetric anyway.
    """

    result = {
        n: Counter()
        for n in range(ell)
    }

    for a in range(
        1,
        ell,
    ):
        for b in range(
            1,
            ell,
        ):

            n = (
                a * b
            ) % ell

            T = (
                a * a
                - a * b
                + b * b
            ) % ell

            result[n][T] += 1

    return result


def theoretical_bayes_accuracy(
    pair_map: Dict[int, Counter],
) -> float:

    """
    Optimal prediction accuracy if the only information available is n.
    """

    total = 0
    correct = 0

    for counter in pair_map.values():
        n = sum(counter.values())

        if n == 0:
            continue

        correct += max(
            counter.values()
        )

        total += n

    if total == 0:
        return float("nan")

    return correct / total


# ============================================================================
# PRIME-RESIDUE DISTRIBUTION
# ============================================================================

def prime_residue_distribution(
    primes: Sequence[int],
    ell: int,
) -> Counter:

    c = Counter()

    for p in primes:
        r = p % ell

        if r != 0:
            c[r] += 1

    return c


# ============================================================================
# EMPIRICAL SYNTHETIC NULL
# ============================================================================

def synthetic_empirical_pairs(
    residue_counts: Counter,
    ell: int,
    samples: int,
    rng: random.Random,
) -> List[Tuple[int, int]]:

    residues = list(
        residue_counts.keys()
    )

    weights = [
        residue_counts[r]
        for r in residues
    ]

    pairs = []

    for _ in range(samples):

        a = rng.choices(
            residues,
            weights=weights,
            k=1,
        )[0]

        b = rng.choices(
            residues,
            weights=weights,
            k=1,
        )[0]

        pairs.append(
            (a, b)
        )

    return pairs


def synthetic_uniform_pairs(
    ell: int,
    samples: int,
    rng: random.Random,
) -> List[Tuple[int, int]]:

    residues = list(
        range(1, ell)
    )

    return [
        (
            rng.choice(residues),
            rng.choice(residues),
        )
        for _ in range(samples)
    ]


# ============================================================================
# TRANSFORM
# ============================================================================

def pair_to_nt(
    a: int,
    b: int,
    ell: int,
) -> Tuple[int, int]:

    n = (
        a * b
    ) % ell

    T = (
        a * a
        - a * b
        + b * b
    ) % ell

    return n, T


# ============================================================================
# LOOKUP MODEL
# ============================================================================

def build_lookup(
    pairs: Sequence[Tuple[int, int]],
    ell: int,
) -> Dict[int, int]:

    buckets = defaultdict(
        Counter
    )

    for a, b in pairs:

        n, T = pair_to_nt(
            a,
            b,
            ell,
        )

        buckets[n][T] += 1

    return {
        n: counter.most_common(1)[0][0]
        for n, counter in buckets.items()
    }


def evaluate_lookup_on_targets(
    train: Sequence[Target],
    test: Sequence[Target],
    ell: int,
) -> Tuple[
    float,
    float,
    float,
]:

    """
    Returns:

        accuracy
        balanced accuracy
        coverage
    """

    model = defaultdict(
        Counter
    )

    for t in train:

        n = t.n % ell
        T = (
            t.s * t.s
            - 3 * t.n
        ) % ell

        model[n][T] += 1

    predictor = {
        n: c.most_common(1)[0][0]
        for n, c in model.items()
    }

    truth = []
    pred = []

    for t in test:

        n = t.n % ell

        if n not in predictor:
            continue

        T = (
            t.s * t.s
            - 3 * t.n
        ) % ell

        truth.append(T)
        pred.append(
            predictor[n]
        )

    if not truth:
        return (
            float("nan"),
            float("nan"),
            0.0,
        )

    acc = sum(
        a == b
        for a, b in zip(
            truth,
            pred,
        )
    ) / len(truth)

    labels = sorted(
        set(truth)
    )

    recalls = []

    for label in labels:

        idx = [
            i
            for i, y in enumerate(truth)
            if y == label
        ]

        if not idx:
            continue

        recalls.append(
            sum(
                pred[i] == label
                for i in idx
            ) / len(idx)
        )

    bal = (
        statistics.fmean(recalls)
        if recalls
        else float("nan")
    )

    return (
        acc,
        bal,
        len(truth) / len(test),
    )


# ============================================================================
# DISTRIBUTION COMPARISON
# ============================================================================

def distribution_error(
    observed: Counter,
    expected: Counter,
    support: Iterable[int],
) -> float:

    """
    Total variation distance.
    """

    total_o = sum(
        observed.values()
    )

    total_e = sum(
        expected.values()
    )

    if total_o == 0 or total_e == 0:
        return float("nan")

    tv = 0.0

    for x in support:

        po = (
            observed[x]
            / total_o
        )

        pe = (
            expected[x]
            / total_e
        )

        tv += abs(
            po - pe
        )

    return 0.5 * tv


# ============================================================================
# CONDITIONAL TABLES
# ============================================================================

def conditional_counter(
    pairs: Sequence[Tuple[int, int]],
    ell: int,
) -> Dict[int, Counter]:

    table = {
        n: Counter()
        for n in range(ell)
    }

    for a, b in pairs:

        n, T = pair_to_nt(
            a,
            b,
            ell,
        )

        table[n][T] += 1

    return table


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("KAPPA EXPERIMENT 89")
    print("EXACT FINITE-FIELD RESOLVENT NULL")
    print("T = s^2 - 3n = p^2 - pq + q^2")
    print("N -> T LOCAL SIGNAL DECOMPOSITION")
    print("EXACT THEORY VS PRIME EMPIRICAL VS MATCHED NULL")
    print("STRICT TARGET HOLDOUT")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    overall_start = time.perf_counter()

    # ------------------------------------------------------------------
    # Population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = {len(targets)}"
    )

    for i, t in enumerate(
        targets[:24],
        1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if NUM_TARGETS > 24:
        print(
            "... remaining targets omitted"
        )

    train = targets[
        :TRAIN_TARGETS
    ]

    test = targets[
        TRAIN_TARGETS:
    ]

    print("\n2. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train)}"
    )
    print(
        f"test targets = {len(test)}"
    )

    # ------------------------------------------------------------------
    # Main finite-field analysis
    # ------------------------------------------------------------------

    summary = []

    for ell in MODULI:

        print("\n")
        print("=" * 78)
        print(f"MODULUS ell = {ell}")
        print("=" * 78)

        # --------------------------------------------------------------
        # Exact field map
        # --------------------------------------------------------------

        field_map = exact_pair_map(
            ell
        )

        bayes_acc = (
            theoretical_bayes_accuracy(
                field_map
            )
        )

        print("\n3. EXACT FINITE-FIELD MAP")
        print("-" * 78)

        for n in range(ell):

            counter = field_map[n]

            if not counter:
                continue

            total = sum(
                counter.values()
            )

            top = counter.most_common()

            print(
                f"n={n:3d} "
                f"pairs={total:4d} "
                f"states={len(counter):2d} "
                f"top={top[:5]}"
            )

        print(
            f"theoretical Bayes accuracy = "
            f"{bayes_acc:.6f}"
        )

        # --------------------------------------------------------------
        # Prime residue distribution
        # --------------------------------------------------------------

        residue_counts = (
            prime_residue_distribution(
                primes,
                ell,
            )
        )

        print("\n4. PRIME RESIDUE DISTRIBUTION")
        print("-" * 78)

        print(
            sorted(
                residue_counts.items()
            )
        )

        # --------------------------------------------------------------
        # Actual target distribution
        # --------------------------------------------------------------

        actual_pairs = [
            (
                t.p % ell,
                t.q % ell,
            )
            for t in targets
        ]

        actual_table = conditional_counter(
            actual_pairs,
            ell,
        )

        # --------------------------------------------------------------
        # Synthetic matched null
        # --------------------------------------------------------------

        rng = random.Random(
            RNG_SEED + ell
        )

        empirical_pairs = (
            synthetic_empirical_pairs(
                residue_counts,
                ell,
                SYNTHETIC_SAMPLES,
                rng,
            )
        )

        uniform_pairs = (
            synthetic_uniform_pairs(
                ell,
                SYNTHETIC_SAMPLES,
                rng,
            )
        )

        empirical_table = (
            conditional_counter(
                empirical_pairs,
                ell,
            )
        )

        uniform_table = (
            conditional_counter(
                uniform_pairs,
                ell,
            )
        )

        # --------------------------------------------------------------
        # Conditional distribution distance
        # --------------------------------------------------------------

        actual_tv = []
        empirical_tv = []
        uniform_tv = []

        for n in range(ell):

            if not actual_table[n]:
                continue

            actual_tv.append(
                distribution_error(
                    actual_table[n],
                    field_map[n],
                    range(ell),
                )
            )

            empirical_tv.append(
                distribution_error(
                    actual_table[n],
                    empirical_table[n],
                    range(ell),
                )
            )

            uniform_tv.append(
                distribution_error(
                    actual_table[n],
                    uniform_table[n],
                    range(ell),
                )
            )

        # --------------------------------------------------------------
        # OOS actual lookup
        # --------------------------------------------------------------

        acc, bal, coverage = (
            evaluate_lookup_on_targets(
                train,
                test,
                ell,
            )
        )

        # --------------------------------------------------------------
        # Synthetic lookup Bayes accuracy
        # --------------------------------------------------------------

        empirical_model = build_lookup(
            empirical_pairs,
            ell,
        )

        empirical_correct = 0
        empirical_total = 0

        for t in test:

            n = t.n % ell

            if n not in empirical_model:
                continue

            actual_T = (
                t.s * t.s
                - 3 * t.n
            ) % ell

            empirical_correct += (
                empirical_model[n]
                == actual_T
            )

            empirical_total += 1

        empirical_acc = (
            empirical_correct
            / empirical_total
            if empirical_total
            else float("nan")
        )

        # --------------------------------------------------------------
        # Print
        # --------------------------------------------------------------

        mean_actual_tv = (
            statistics.fmean(actual_tv)
            if actual_tv
            else float("nan")
        )

        mean_empirical_tv = (
            statistics.fmean(empirical_tv)
            if empirical_tv
            else float("nan")
        )

        mean_uniform_tv = (
            statistics.fmean(uniform_tv)
            if uniform_tv
            else float("nan")
        )

        print("\n5. CONDITIONAL DISTRIBUTION DISTANCE")
        print("-" * 78)

        print(
            f"actual vs exact-field TV = "
            f"{mean_actual_tv:.6f}"
        )

        print(
            f"actual vs empirical-null TV = "
            f"{mean_empirical_tv:.6f}"
        )

        print(
            f"actual vs uniform-null TV = "
            f"{mean_uniform_tv:.6f}"
        )

        print("\n6. OOS LOOKUP")
        print("-" * 78)

        print(
            f"actual train-state lookup "
            f"accuracy = {acc:.6f}"
        )

        print(
            f"actual train-state lookup "
            f"balanced = {bal:.6f}"
        )

        print(
            f"actual coverage = {coverage:.6f}"
        )

        print(
            f"matched empirical-residue "
            f"synthetic accuracy = "
            f"{empirical_acc:.6f}"
        )

        print(
            f"exact finite-field Bayes "
            f"accuracy = {bayes_acc:.6f}"
        )

        summary.append(
            (
                ell,
                bayes_acc,
                acc,
                bal,
                coverage,
                empirical_acc,
                mean_actual_tv,
                mean_empirical_tv,
                mean_uniform_tv,
            )
        )

        # --------------------------------------------------------------
        # Detailed interesting small modulus
        # --------------------------------------------------------------

        if ell in DETAIL_MODULI:

            print("\n7. DETAILED CONDITIONAL STRUCTURE")
            print("-" * 78)

            for n in range(ell):

                theory = field_map[n]

                if not theory:
                    continue

                actual = actual_table[n]

                if not actual:
                    continue

                print(
                    f"N={n:3d}: "
                    f"theory={theory.most_common()} "
                    f"actual={actual.most_common()}"
                )

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("8. CROSS-MODULUS SUMMARY")
    print("=" * 78)

    print(
        "ell | theory | actual_OOS | balanced | "
        "coverage | empirical_null | TV(field) | TV(emp)"
    )

    for row in summary:

        (
            ell,
            theory,
            acc,
            bal,
            cov,
            emp,
            tvf,
            tve,
            _tvu,
        ) = row

        print(
            f"{ell:3d} | "
            f"{theory:.6f} | "
            f"{acc:.6f} | "
            f"{bal:.6f} | "
            f"{cov:.6f} | "
            f"{emp:.6f} | "
            f"{tvf:.6f} | "
            f"{tve:.6f}"
        )

    # ------------------------------------------------------------------
    # Key diagnostic
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The key question is whether the observed N -> T signal "
        "exceeds what finite-field algebra itself predicts."
    )

    print()
    print(
        "Interpretation A:"
    )
    print(
        "    actual OOS ~= empirical synthetic null"
    )
    print(
        "    and actual conditional tables resemble exact F_l theory"
    )
    print(
        "    => the earlier signal is explained by ordinary "
        "finite-field residue structure."
    )

    print()
    print(
        "Interpretation B:"
    )
    print(
        "    actual OOS substantially exceeds empirical null"
    )
    print(
        "    across several ell"
    )
    print(
        "    => there is additional structure beyond the basic "
        "residue-pair mechanism."
    )

    print()
    print(
        "Interpretation C:"
    )
    print(
        "    ell=7 is special but larger ell behave normally"
    )
    print(
        "    => investigate a modulus-specific cyclotomic identity "
        "rather than a universal resolvent law."
    )

    print()
    print(
        "Most important:"
    )
    print(
        "This experiment does NOT use H6/H8 as a predictor."
    )
    print(
        "It directly tests the finite-field structure of"
    )
    print(
        "    T = p^2 - pq + q^2"
    )
    print(
        "conditioned on"
    )
    print(
        "    N = pq."
    )

    print()
    print(
        "That isolates the most plausible confounder from "
        "Experiment 88."
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - overall_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 89 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

