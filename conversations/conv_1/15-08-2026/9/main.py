#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 67
JOINT N-ONLY SIGNATURE / CONDITIONAL INFORMATION TEST
CROSS-MODULUS INFORMATION FOR RELATIVE LEGENDRE ORIENTATION
TARGET + MODULUS + 2D HOLDOUT
GLOBAL + WITHIN-TARGET PERMUTATION NULL
NAIVE LOGISTIC + INFORMATION DIAGNOSTICS
NO CSV OUTPUT
NO SKLEARN
==============================================================================

QUESTION
--------
Does the joint collection of N-only residue information across several
moduli contain predictive information about relative Legendre orientation
at a held-out modulus?

For target n and modulus ell:

    y(ell) = epsilon_ell * epsilon_reference

where epsilon is the common Legendre sign of the two conjugate
discriminants.

IMPORTANT
---------
The true factorization is used ONLY to construct the ground-truth label.

ALL predictor features are computed from n alone.

For a response modulus ell, features belonging to ell itself are excluded.
This makes the main test genuinely cross-modulus.

TESTS
-----
1. Target holdout
2. Leave-one-modulus-out
3. Two-dimensional unseen-target + unseen-modulus holdout
4. Global permutation null
5. Within-target permutation null
6. Simple N-only baselines
7. Joint categorical information diagnostics
8. Per-target test accuracy
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Sequence


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 67067

NUM_TARGETS = 120

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

# FIXED: correct variable name.
CYCLOTOMIC = [
    7,
    13,
    19,
    31,
    37,
    61,
    67,
    79,
    127,
    307,
    331,
    631,
    1723,
]

RAW_CONTROLS = [
    673,
    1109,
    2309,
    2351,
    4421,
    4561,
    4759,
    6211,
    7879,
    7951,
    8273,
    8689,
    9781,
]

VALID_CONTROLS = [
    ell for ell in RAW_CONTROLS
    if ell % 3 == 1
]

CYC_REFERENCE = 7
CONTROL_REFERENCE = 673

TRAIN_TARGET_FRACTION = 0.75

PERMUTATIONS = 400


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


@dataclass(frozen=True)
class Sample:
    target_id: int
    ell: int
    y: int
    x: tuple[float, ...]


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> list[int]:
    """Return all primes in [lo, hi]."""

    if hi < 2:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    root = int(math.isqrt(hi))

    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            count = ((hi - start) // p) + 1
            sieve[start:hi + 1:p] = b"\x00" * count

    return [
        x for x in range(lo, hi + 1)
        if sieve[x]
    ]


def generate_targets(
    primes: Sequence[int],
    count: int,
    seed: int,
) -> list[Target]:
    """
    Deterministic fresh target population.

    p and q are distinct primes in the requested range.
    """

    rng = random.Random(seed)

    targets: list[Target] = []
    used: set[tuple[int, int]] = set()

    while len(targets) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in used:
            continue

        used.add(pair)

        targets.append(
            Target(
                idx=len(targets) + 1,
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return targets


# ============================================================================
# NUMBER THEORY
# ============================================================================

def legendre(a: int, ell: int) -> int:
    """
    Legendre symbol (a / ell).

    Returns:
        +1 quadratic residue
        -1 non-residue
         0 zero
    """

    a %= ell

    if a == 0:
        return 0

    value = pow(a, (ell - 1) // 2, ell)

    if value == 1:
        return 1

    if value == ell - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre result: a={a}, ell={ell}, value={value}"
    )


def roots_cyclotomic(ell: int) -> tuple[int, int]:
    """
    Find the two roots of

        x^2 + x + 1 = 0 mod ell

    for primes ell == 1 mod 3.
    """

    if ell % 3 != 1:
        raise ValueError(
            f"{ell} is not valid for x^2+x+1: ell % 3 != 1"
        )

    exponent = (ell - 1) // 3

    for g in range(2, ell):
        z = pow(g, exponent, ell)

        if z == 1:
            continue

        z2 = (z * z) % ell

        if (
            (z * z + z + 1) % ell == 0
            and (z2 * z2 + z2 + 1) % ell == 0
            and z != z2
        ):
            return min(z, z2), max(z, z2)

    raise RuntimeError(
        f"Could not find roots of x^2+x+1 modulo {ell}"
    )


def discriminant_pair(
    n: int,
    ell: int,
) -> tuple[int, int, int, int]:
    """
    Return:
        w1, w2, d1, d2
    """

    w1, w2 = roots_cyclotomic(ell)

    d1 = (w1 * w1 - 4 * n) % ell
    d2 = (w2 * w2 - 4 * n) % ell

    return w1, w2, d1, d2


def common_orientation(
    n: int,
    ell: int,
) -> int | None:
    """
    Return common Legendre sign if both discriminants are nonzero and
    have the same sign.

    Returns None for:
        zero symbol
        opposite signs
    """

    _, _, d1, d2 = discriminant_pair(n, ell)

    c1 = legendre(d1, ell)
    c2 = legendre(d2, ell)

    if c1 == 0 or c2 == 0:
        return None

    if c1 != c2:
        return None

    if c1 not in (-1, 1):
        raise RuntimeError(
            f"Invalid common orientation: "
            f"n={n}, ell={ell}, c1={c1}, c2={c2}"
        )

    return c1


def n_only_scalar_features(
    n: int,
    ell: int,
) -> tuple[int, int, int, int, int]:
    """
    N-only feature block:

        chi(n)
        chi(n+1)
        chi(n-1)
        chi(4n+1)
        chi(G(n))

    where

        G(n) = 16n^2 + 4n + 1.
    """

    g = 16 * n * n + 4 * n + 1

    return (
        legendre(n, ell),
        legendre(n + 1, ell),
        legendre(n - 1, ell),
        legendre(4 * n + 1, ell),
        legendre(g, ell),
    )


# ============================================================================
# LABEL CONSTRUCTION
# ============================================================================

def build_relative_labels(
    targets: Sequence[Target],
    moduli: Sequence[int],
    reference: int,
) -> dict[tuple[int, int], int]:
    """
    Compute

        y(ell) = epsilon_ell * epsilon_reference.

    Only nondegenerate common-sign cells are retained.
    """

    labels: dict[tuple[int, int], int] = {}

    for target in targets:
        e_ref = common_orientation(target.n, reference)

        if e_ref is None:
            continue

        for ell in moduli:
            if ell == reference:
                continue

            e = common_orientation(target.n, ell)

            if e is None:
                continue

            labels[(target.idx, ell)] = e * e_ref

    return labels


# ============================================================================
# CROSS-MODULUS FEATURE VECTOR
# ============================================================================

def feature_vector(
    n: int,
    target_ell: int,
    panel: Sequence[int],
) -> tuple[float, ...]:
    """
    Build a joint N-only feature vector.

    The response modulus target_ell is excluded.

    For each remaining modulus we include five scalar Legendre features.

    Then add panel aggregates:
        sum of each feature
        positive-minus-negative count of each feature
    """

    blocks: list[tuple[int, ...]] = []

    for ell in panel:
        if ell == target_ell:
            continue

        blocks.append(
            n_only_scalar_features(n, ell)
        )

    out: list[float] = []

    # Raw cross-modulus values.
    for block in blocks:
        out.extend(float(v) for v in block)

    if blocks:
        width = len(blocks[0])

        # Sum of each feature across the panel.
        for j in range(width):
            values = [b[j] for b in blocks]
            out.append(float(sum(values)))

        # Positive-minus-negative balance.
        for j in range(width):
            values = [b[j] for b in blocks]
            pos = sum(v == 1 for v in values)
            neg = sum(v == -1 for v in values)
            out.append(float(pos - neg))

    return tuple(out)


# ============================================================================
# DATASET CONSTRUCTION
# ============================================================================

def make_samples(
    targets: Sequence[Target],
    panel: Sequence[int],
    reference: int,
) -> tuple[list[Sample], int]:
    """
    Build all usable relative-orientation samples.
    """

    labels = build_relative_labels(
        targets,
        panel,
        reference,
    )

    samples: list[Sample] = []

    total_cells = 0

    target_by_id = {
        t.idx: t
        for t in targets
    }

    for target in targets:
        for ell in panel:
            if ell == reference:
                continue

            total_cells += 1

            key = (target.idx, ell)

            if key not in labels:
                continue

            x = feature_vector(
                target.n,
                ell,
                panel,
            )

            samples.append(
                Sample(
                    target_id=target.idx,
                    ell=ell,
                    y=labels[key],
                    x=x,
                )
            )

    degenerate = total_cells - len(samples)

    return samples, degenerate


# ============================================================================
# PURE PYTHON LOGISTIC MODEL
# ============================================================================

class LogisticModel:
    """
    Small deterministic L2-regularized logistic regression.

    Implemented directly so there are no sklearn compatibility warnings.
    """

    def __init__(
        self,
        learning_rate: float = 0.07,
        epochs: int = 700,
        l2: float = 1.5,
    ):
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.l2 = l2

        self.w: list[float] = []
        self.b: float = 0.0

        self.means: list[float] = []
        self.scales: list[float] = []

    @staticmethod
    def sigmoid(z: float) -> float:
        if z >= 0:
            e = math.exp(-z)
            return 1.0 / (1.0 + e)

        e = math.exp(z)
        return e / (1.0 + e)

    def fit(
        self,
        X: Sequence[Sequence[float]],
        y: Sequence[int],
    ) -> None:
        if not X:
            raise ValueError("Cannot train on empty data.")

        if len(X) != len(y):
            raise ValueError(
                f"X/y mismatch: {len(X)} != {len(y)}"
            )

        feature_count = len(X[0])

        self.means = []
        self.scales = []

        for j in range(feature_count):
            column = [
                float(row[j])
                for row in X
            ]

            mu = statistics.fmean(column)

            variance = statistics.fmean(
                (v - mu) ** 2
                for v in column
            )

            scale = math.sqrt(variance)

            if scale < 1e-12:
                scale = 1.0

            self.means.append(mu)
            self.scales.append(scale)

        Z = self.transform(X)

        self.w = [0.0] * feature_count
        self.b = 0.0

        n = float(len(Z))

        for _ in range(self.epochs):
            grad_w = [0.0] * feature_count
            grad_b = 0.0

            for row, label in zip(Z, y):
                target = 1.0 if label == 1 else 0.0

                z = (
                    self.b
                    + sum(
                        wi * xi
                        for wi, xi in zip(self.w, row)
                    )
                )

                p = self.sigmoid(z)
                error = p - target

                grad_b += error

                for j, value in enumerate(row):
                    grad_w[j] += error * value

            grad_b /= n

            for j in range(feature_count):
                grad_w[j] = (
                    grad_w[j] / n
                    + self.l2 * self.w[j]
                )

            self.b -= self.learning_rate * grad_b

            for j in range(feature_count):
                self.w[j] -= (
                    self.learning_rate * grad_w[j]
                )

    def transform(
        self,
        X: Sequence[Sequence[float]],
    ) -> list[list[float]]:
        return [
            [
                (
                    float(row[j]) - self.means[j]
                ) / self.scales[j]
                for j in range(len(self.means))
            ]
            for row in X
        ]

    def predict_one(
        self,
        x: Sequence[float],
    ) -> int:
        if not self.w:
            raise RuntimeError("Model is not fitted.")

        zrow = [
            (
                float(x[j]) - self.means[j]
            ) / self.scales[j]
            for j in range(len(self.w))
        ]

        z = (
            self.b
            + sum(
                wi * xi
                for wi, xi in zip(self.w, zrow)
            )
        )

        return 1 if z >= 0 else -1

    def predict(
        self,
        X: Sequence[Sequence[float]],
    ) -> list[int]:
        return [
            self.predict_one(row)
            for row in X
        ]


# ============================================================================
# METRICS
# ============================================================================

def accuracy(
    samples: Sequence[Sample],
    predictions: Sequence[int],
) -> float:
    if not samples:
        return float("nan")

    return sum(
        pred == sample.y
        for sample, pred in zip(samples, predictions)
    ) / len(samples)


def balanced_accuracy(
    samples: Sequence[Sample],
    predictions: Sequence[int],
) -> float:
    if not samples:
        return float("nan")

    positives = 0
    negatives = 0
    true_positive = 0
    true_negative = 0

    for sample, pred in zip(samples, predictions):
        if sample.y == 1:
            positives += 1
            if pred == 1:
                true_positive += 1
        else:
            negatives += 1
            if pred == -1:
                true_negative += 1

    if positives == 0 or negatives == 0:
        return float("nan")

    return 0.5 * (
        true_positive / positives
        + true_negative / negatives
    )


def majority_accuracy(
    samples: Sequence[Sample],
) -> float:
    if not samples:
        return float("nan")

    counts = Counter(
        sample.y
        for sample in samples
    )

    return max(counts.values()) / len(samples)


# ============================================================================
# TARGET SPLIT
# ============================================================================

def split_targets(
    targets: Sequence[Target],
    seed: int,
) -> tuple[set[int], set[int]]:

    ids = [
        target.idx
        for target in targets
    ]

    rng = random.Random(seed)
    rng.shuffle(ids)

    train_count = int(
        round(
            len(ids)
            * TRAIN_TARGET_FRACTION
        )
    )

    train_ids = set(
        ids[:train_count]
    )

    test_ids = set(
        ids[train_count:]
    )

    return train_ids, test_ids


def filter_targets(
    samples: Sequence[Sample],
    ids: set[int],
) -> list[Sample]:
    return [
        sample
        for sample in samples
        if sample.target_id in ids
    ]


# ============================================================================
# MODEL SCORING
# ============================================================================

def train_test_model(
    train_samples: Sequence[Sample],
    test_samples: Sequence[Sample],
) -> tuple[float, float, float]:
    """
    Returns:
        train accuracy
        test accuracy
        test balanced accuracy
    """

    if not train_samples or not test_samples:
        return (
            float("nan"),
            float("nan"),
            float("nan"),
        )

    model = LogisticModel()

    model.fit(
        [sample.x for sample in train_samples],
        [sample.y for sample in train_samples],
    )

    train_pred = model.predict(
        [sample.x for sample in train_samples]
    )

    test_pred = model.predict(
        [sample.x for sample in test_samples]
    )

    return (
        accuracy(train_samples, train_pred),
        accuracy(test_samples, test_pred),
        balanced_accuracy(test_samples, test_pred),
    )


# ============================================================================
# INFORMATION THEORY
# ============================================================================

def binary_entropy(
    values: Sequence[int],
) -> float:
    if not values:
        return float("nan")

    counts = Counter(values)
    n = len(values)

    entropy = 0.0

    for count in counts.values():
        p = count / n

        if p > 0:
            entropy -= p * math.log2(p)

    return entropy


def empirical_conditional_entropy(
    x_values: Sequence[tuple],
    y_values: Sequence[int],
    min_group: int = 3,
) -> float:
    """
    Estimate H(Y|X) by grouping identical joint signatures.

    Tiny groups are pooled into one bucket so singleton groups do not
    trivially force conditional entropy to zero.
    """

    if not x_values:
        return float("nan")

    groups: dict[tuple, list[int]] = defaultdict(list)

    for x, y in zip(x_values, y_values):
        groups[x].append(y)

    pooled: dict[tuple, list[int]] = defaultdict(list)

    for key, labels in groups.items():
        if len(labels) < min_group:
            pooled[("__OTHER__",)].extend(labels)
        else:
            pooled[key].extend(labels)

    n = len(y_values)

    result = 0.0

    for labels in pooled.values():
        weight = len(labels) / n
        result += (
            weight
            * binary_entropy(labels)
        )

    return result


# ============================================================================
# EXACT JOINT SIGNATURE
# ============================================================================

def exact_joint_signature(
    n: int,
    target_ell: int,
    panel: Sequence[int],
) -> tuple[int, ...]:
    """
    Exact discrete joint N-only signature.

    Excludes target_ell.
    """

    signature: list[int] = []

    for ell in panel:
        if ell == target_ell:
            continue

        signature.extend(
            n_only_scalar_features(
                n,
                ell,
            )
        )

    return tuple(signature)


def information_summary(
    targets: Sequence[Target],
    samples: Sequence[Sample],
    panel: Sequence[int],
) -> tuple[float, float, float]:
    """
    Returns:
        H(Y)
        H(Y|joint signature)
        I(Y; signature)
    """

    target_map = {
        target.idx: target
        for target in targets
    }

    xs: list[tuple[int, ...]] = []
    ys: list[int] = []

    for sample in samples:
        target = target_map[
            sample.target_id
        ]

        xs.append(
            exact_joint_signature(
                target.n,
                sample.ell,
                panel,
            )
        )

        ys.append(sample.y)

    h_y = binary_entropy(ys)

    h_y_given_x = (
        empirical_conditional_entropy(
            xs,
            ys,
            min_group=3,
        )
    )

    if (
        math.isnan(h_y)
        or math.isnan(h_y_given_x)
    ):
        mutual_information = float("nan")
    else:
        mutual_information = max(
            0.0,
            h_y - h_y_given_x,
        )

    return (
        h_y,
        h_y_given_x,
        mutual_information,
    )


# ============================================================================
# LEAVE-ONE-MODULUS-OUT
# ============================================================================

def leave_one_modulus_out(
    samples: Sequence[Sample],
    moduli: Sequence[int],
    train_ids: set[int],
    test_ids: set[int],
    reference: int,
) -> list[tuple[int, float, float, int, int]]:

    results = []

    for heldout in moduli:
        if heldout == reference:
            continue

        train = [
            sample
            for sample in samples
            if (
                sample.ell == heldout
                and sample.target_id in train_ids
            )
        ]

        test = [
            sample
            for sample in samples
            if (
                sample.ell == heldout
                and sample.target_id in test_ids
            )
        ]

        if not train or not test:
            results.append(
                (
                    heldout,
                    float("nan"),
                    float("nan"),
                    len(train),
                    len(test),
                )
            )
            continue

        model = LogisticModel()

        model.fit(
            [sample.x for sample in train],
            [sample.y for sample in train],
        )

        predictions = model.predict(
            [sample.x for sample in test]
        )

        results.append(
            (
                heldout,
                accuracy(test, predictions),
                balanced_accuracy(test, predictions),
                len(train),
                len(test),
            )
        )

    return results


# ============================================================================
# TWO-DIMENSIONAL HOLDOUT
# ============================================================================

def two_dimensional_holdout(
    samples: Sequence[Sample],
    heldout_moduli: Sequence[int],
    train_ids: set[int],
    test_ids: set[int],
) -> tuple[float, float, int]:

    train = [
        sample
        for sample in samples
        if (
            sample.ell not in heldout_moduli
            and sample.target_id in train_ids
        )
    ]

    test = [
        sample
        for sample in samples
        if (
            sample.ell in heldout_moduli
            and sample.target_id in test_ids
        )
    ]

    if not train or not test:
        return (
            float("nan"),
            float("nan"),
            len(test),
        )

    model = LogisticModel()

    model.fit(
        [sample.x for sample in train],
        [sample.y for sample in train],
    )

    predictions = model.predict(
        [sample.x for sample in test]
    )

    return (
        accuracy(test, predictions),
        balanced_accuracy(test, predictions),
        len(test),
    )


# ============================================================================
# PERMUTATION NULLS
# ============================================================================

def permute_global(
    samples: Sequence[Sample],
    seed: int,
) -> list[Sample]:

    rng = random.Random(seed)

    labels = [
        sample.y
        for sample in samples
    ]

    rng.shuffle(labels)

    return [
        Sample(
            target_id=sample.target_id,
            ell=sample.ell,
            y=label,
            x=sample.x,
        )
        for sample, label in zip(
            samples,
            labels,
        )
    ]


def permute_within_target(
    samples: Sequence[Sample],
    seed: int,
) -> list[Sample]:

    rng = random.Random(seed)

    groups: dict[int, list[Sample]] = defaultdict(list)

    for sample in samples:
        groups[sample.target_id].append(sample)

    output: list[Sample] = []

    for group in groups.values():
        labels = [
            sample.y
            for sample in group
        ]

        rng.shuffle(labels)

        for sample, label in zip(
            group,
            labels,
        ):
            output.append(
                Sample(
                    target_id=sample.target_id,
                    ell=sample.ell,
                    y=label,
                    x=sample.x,
                )
            )

    return output


def permutation_score(
    samples: Sequence[Sample],
    train_ids: set[int],
    test_ids: set[int],
    within_target: bool,
    seed: int,
) -> float:

    if within_target:
        randomized = permute_within_target(
            samples,
            seed,
        )
    else:
        randomized = permute_global(
            samples,
            seed,
        )

    train = filter_targets(
        randomized,
        train_ids,
    )

    test = filter_targets(
        randomized,
        test_ids,
    )

    _, _, score = train_test_model(
        train,
        test,
    )

    return score


def permutation_pvalue(
    samples: Sequence[Sample],
    train_ids: set[int],
    test_ids: set[int],
    observed: float,
    within_target: bool,
    permutations: int,
    seed: int,
) -> float:

    if math.isnan(observed):
        return float("nan")

    rng = random.Random(seed)

    exceed = 0

    for _ in range(permutations):
        random_seed = rng.randrange(
            1 << 62
        )

        score = permutation_score(
            samples,
            train_ids,
            test_ids,
            within_target,
            random_seed,
        )

        if (
            not math.isnan(score)
            and score >= observed
        ):
            exceed += 1

    return (
        exceed + 1
    ) / (
        permutations + 1
    )


# ============================================================================
# VALIDATION
# ============================================================================

def validate_moduli(
    moduli: Sequence[int],
) -> None:

    print("4. MODULUS VALIDATION")
    print("-" * 78)

    failures = 0

    for ell in moduli:
        r1, r2 = roots_cyclotomic(ell)

        ok = (
            (r1 * r1 + r1 + 1) % ell == 0
            and
            (r2 * r2 + r2 + 1) % ell == 0
            and
            r1 != r2
            and
            (r1 * r2) % ell == 1
        )

        print(
            f"ell={ell:5d} "
            f"roots=({r1:5d},{r2:5d}) "
            f"identity_ok={ok}"
        )

        if not ok:
            failures += 1

    print(
        f"validation failures = {failures}"
    )

    if failures:
        raise RuntimeError(
            "Modulus validation failed."
        )

    print("status = PASS")
    print()


# ============================================================================
# SIMPLE BASELINES
# ============================================================================

def simple_baseline(
    samples: Sequence[Sample],
    targets: Sequence[Target],
    panel: Sequence[int],
    feature_index: int,
    train_ids: set[int],
    test_ids: set[int],
) -> tuple[float, float]:

    target_map = {
        target.idx: target
        for target in targets
    }

    transformed: list[Sample] = []

    for sample in samples:
        target = target_map[
            sample.target_id
        ]

        values = []

        for ell in panel:
            if ell == sample.ell:
                continue

            block = n_only_scalar_features(
                target.n,
                ell,
            )

            values.append(
                block[feature_index]
            )

        if not values:
            continue

        transformed.append(
            Sample(
                target_id=sample.target_id,
                ell=sample.ell,
                y=sample.y,
                x=(statistics.fmean(values),),
            )
        )

    train = filter_targets(
        transformed,
        train_ids,
    )

    test = filter_targets(
        transformed,
        test_ids,
    )

    train_acc, test_acc, _ = train_test_model(
        train,
        test,
    )

    return train_acc, test_acc


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start_total = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 67")
    print("JOINT N-ONLY SIGNATURE / CONDITIONAL INFORMATION TEST")
    print("CROSS-MODULUS INFORMATION FOR RELATIVE LEGENDRE ORIENTATION")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("GLOBAL + WITHIN-TARGET PERMUTATION NULL")
    print("NAIVE LOGISTIC + INFORMATION DIAGNOSTICS")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    # ----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ----------------------------------------------------------------------

    prime_start = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    print("1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time  = "
        f"{time.perf_counter() - prime_start:.6f}s"
    )
    print()

    # ----------------------------------------------------------------------
    # 2. TARGETS
    # ----------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
        SEED,
    )

    print("2. TARGETS")
    print("-" * 78)
    print(
        f"total targets = {len(targets)}"
    )

    for target in targets[:24]:
        print(
            f"target {target.idx:3d}: "
            f"p={target.p} "
            f"q={target.q} "
            f"n={target.n} "
            f"s={target.s}"
        )

    print(
        "... remaining targets generated but omitted"
    )
    print()

    # ----------------------------------------------------------------------
    # 3. MODULUS FAMILIES
    # ----------------------------------------------------------------------

    cyclotomic = list(CYCLOTOMIC)
    controls = list(VALID_CONTROLS)

    excluded_controls = [
        ell
        for ell in RAW_CONTROLS
        if ell not in controls
    ]

    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(
        f"cyclotomic        = {cyclotomic}"
    )
    print(
        f"valid controls    = {controls}"
    )
    print(
        f"excluded controls = {excluded_controls}"
    )
    print()

    validate_moduli(
        cyclotomic + controls
    )

    # ----------------------------------------------------------------------
    # 5. TARGET HOLDOUT
    # ----------------------------------------------------------------------

    train_ids, test_ids = split_targets(
        targets,
        SEED + 1,
    )

    print("5. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train_ids)}"
    )
    print(
        f"test targets     = {len(test_ids)}"
    )
    print(
        f"test ids         = {sorted(test_ids)}"
    )
    print()

    # ----------------------------------------------------------------------
    # 6. DATASET CONSTRUCTION
    # ----------------------------------------------------------------------

    print("6. CROSS-MODULUS DATASET")
    print("-" * 78)

    dataset_start = time.perf_counter()

    cyc_samples, cyc_degenerate = make_samples(
        targets,
        cyclotomic,
        CYC_REFERENCE,
    )

    ctrl_samples, ctrl_degenerate = make_samples(
        targets,
        controls,
        CONTROL_REFERENCE,
    )

    print(
        f"cyclotomic samples   = {len(cyc_samples)}"
    )
    print(
        f"cyclotomic degenerate= {cyc_degenerate}"
    )
    print(
        f"control samples      = {len(ctrl_samples)}"
    )
    print(
        f"control degenerate   = {ctrl_degenerate}"
    )
    print(
        f"construction time    = "
        f"{time.perf_counter() - dataset_start:.6f}s"
    )
    print()

    # ----------------------------------------------------------------------
    # 7. LABEL DISTRIBUTIONS
    # ----------------------------------------------------------------------

    print("7. LABEL DISTRIBUTIONS")
    print("-" * 78)

    for name, samples in (
        ("cyclotomic", cyc_samples),
        ("control", ctrl_samples),
    ):
        positive = sum(
            sample.y == 1
            for sample in samples
        )

        negative = sum(
            sample.y == -1
            for sample in samples
        )

        print(
            f"{name:12s} "
            f"positive={positive:5d} "
            f"negative={negative:5d} "
            f"majority={majority_accuracy(samples):.6f}"
        )

    print()

    # ----------------------------------------------------------------------
    # 8. TARGET HOLDOUT
    # ----------------------------------------------------------------------

    print("8. TARGET HOLDOUT")
    print("-" * 78)

    cyc_train = filter_targets(
        cyc_samples,
        train_ids,
    )

    cyc_test = filter_targets(
        cyc_samples,
        test_ids,
    )

    ctrl_train = filter_targets(
        ctrl_samples,
        train_ids,
    )

    ctrl_test = filter_targets(
        ctrl_samples,
        test_ids,
    )

    cyc_train_acc, cyc_test_acc, cyc_test_bal = (
        train_test_model(
            cyc_train,
            cyc_test,
        )
    )

    ctrl_train_acc, ctrl_test_acc, ctrl_test_bal = (
        train_test_model(
            ctrl_train,
            ctrl_test,
        )
    )

    print(
        f"cyclotomic train={cyc_train_acc:.6f} "
        f"test={cyc_test_acc:.6f} "
        f"balanced={cyc_test_bal:.6f}"
    )

    print(
        f"control    train={ctrl_train_acc:.6f} "
        f"test={ctrl_test_acc:.6f} "
        f"balanced={ctrl_test_bal:.6f}"
    )

    print(
        f"test delta = "
        f"{cyc_test_acc - ctrl_test_acc:+.6f}"
    )

    print()

    # ----------------------------------------------------------------------
    # 9. INFORMATION THEORY
    # ----------------------------------------------------------------------

    print("9. JOINT SIGNATURE INFORMATION")
    print("-" * 78)

    for name, samples, panel in (
        (
            "cyclotomic",
            cyc_test,
            cyclotomic,
        ),
        (
            "control",
            ctrl_test,
            controls,
        ),
    ):
        h_y, h_cond, mutual_info = information_summary(
            targets,
            samples,
            panel,
        )

        print(
            f"{name:12s} "
            f"H(Y)={h_y:.6f} "
            f"H(Y|X)={h_cond:.6f} "
            f"I={mutual_info:.6f} bits"
        )

    print()

    # ----------------------------------------------------------------------
    # 10. LEAVE-ONE-MODULUS-OUT
    # ----------------------------------------------------------------------

    print("10. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    cyc_loo = leave_one_modulus_out(
        cyc_samples,
        cyclotomic,
        train_ids,
        test_ids,
        CYC_REFERENCE,
    )

    ctrl_loo = leave_one_modulus_out(
        ctrl_samples,
        controls,
        train_ids,
        test_ids,
        CONTROL_REFERENCE,
    )

    for (
        ell,
        acc,
        bal,
        train_n,
        test_n,
    ) in cyc_loo:
        if math.isnan(acc):
            print(
                f"cyclo ell={ell:5d} "
                f"accuracy=nan "
                f"train={train_n} test={test_n}"
            )
        else:
            print(
                f"cyclo ell={ell:5d} "
                f"accuracy={acc:.6f} "
                f"balanced={bal:.6f} "
                f"train={train_n} test={test_n}"
            )

    print()

    for (
        ell,
        acc,
        bal,
        train_n,
        test_n,
    ) in ctrl_loo:
        if math.isnan(acc):
            print(
                f"control ell={ell:5d} "
                f"accuracy=nan "
                f"train={train_n} test={test_n}"
            )
        else:
            print(
                f"control ell={ell:5d} "
                f"accuracy={acc:.6f} "
                f"balanced={bal:.6f} "
                f"train={train_n} test={test_n}"
            )

    cyc_loo_valid = [
        acc
        for _, acc, _, _, _ in cyc_loo
        if not math.isnan(acc)
    ]

    ctrl_loo_valid = [
        acc
        for _, acc, _, _, _ in ctrl_loo
        if not math.isnan(acc)
    ]

    if cyc_loo_valid:
        print(
            f"cyclotomic LOO mean   = "
            f"{statistics.fmean(cyc_loo_valid):.6f}"
        )
        print(
            f"cyclotomic LOO median = "
            f"{statistics.median(cyc_loo_valid):.6f}"
        )
    else:
        print(
            "cyclotomic LOO mean   = nan"
        )
        print(
            "cyclotomic LOO median = nan"
        )

    if ctrl_loo_valid:
        print(
            f"control LOO mean      = "
            f"{statistics.fmean(ctrl_loo_valid):.6f}"
        )
        print(
            f"control LOO median    = "
            f"{statistics.median(ctrl_loo_valid):.6f}"
        )
    else:
        print(
            "control LOO mean      = nan"
        )
        print(
            "control LOO median    = nan"
        )

    print()

    # ----------------------------------------------------------------------
    # 11. TWO-DIMENSIONAL HOLDOUT
    # ----------------------------------------------------------------------

    print("11. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)

    heldout_cyc_moduli = [
        61,
        127,
        307,
    ]

    heldout_ctrl_moduli = [
        6211,
        7951,
    ]

    cyc_2d_acc, cyc_2d_bal, cyc_2d_n = (
        two_dimensional_holdout(
            cyc_samples,
            heldout_cyc_moduli,
            train_ids,
            test_ids,
        )
    )

    ctrl_2d_acc, ctrl_2d_bal, ctrl_2d_n = (
        two_dimensional_holdout(
            ctrl_samples,
            heldout_ctrl_moduli,
            train_ids,
            test_ids,
        )
    )

    print(
        f"cyclotomic 2D = {cyc_2d_acc:.6f} "
        f"balanced={cyc_2d_bal:.6f} "
        f"n={cyc_2d_n} "
        f"heldout={heldout_cyc_moduli}"
    )

    print(
        f"control 2D    = {ctrl_2d_acc:.6f} "
        f"balanced={ctrl_2d_bal:.6f} "
        f"n={ctrl_2d_n} "
        f"heldout={heldout_ctrl_moduli}"
    )

    if (
        not math.isnan(cyc_2d_acc)
        and not math.isnan(ctrl_2d_acc)
    ):
        print(
            f"2D delta = "
            f"{cyc_2d_acc - ctrl_2d_acc:+.6f}"
        )
    else:
        print(
            "2D delta = nan"
        )

    print()

    # ----------------------------------------------------------------------
    # 12. PERMUTATION NULLS
    # ----------------------------------------------------------------------

    print("12. PERMUTATION NULL")
    print("-" * 78)

    observed = cyc_test_bal

    p_global = permutation_pvalue(
        cyc_samples,
        train_ids,
        test_ids,
        observed,
        within_target=False,
        permutations=PERMUTATIONS,
        seed=SEED + 100,
    )

    p_within = permutation_pvalue(
        cyc_samples,
        train_ids,
        test_ids,
        observed,
        within_target=True,
        permutations=PERMUTATIONS,
        seed=SEED + 101,
    )

    print(
        f"observed cyclotomic balanced accuracy "
        f"= {observed:.6f}"
    )

    print(
        f"permutations = {PERMUTATIONS}"
    )

    print(
        f"global-label permutation p "
        f"= {p_global:.6f}"
    )

    print(
        f"within-target permutation p "
        f"= {p_within:.6f}"
    )

    print()

    # ----------------------------------------------------------------------
    # 13. SIMPLE N-ONLY BASELINES
    # ----------------------------------------------------------------------

    print("13. SIMPLE N-ONLY BASELINES")
    print("-" * 78)

    baseline_names = [
        "mean_chi(n)",
        "mean_chi(n+1)",
        "mean_chi(n-1)",
        "mean_chi(4n+1)",
        "mean_chi(G)",
    ]

    for feature_index, name in enumerate(
        baseline_names
    ):
        train_acc, test_acc = simple_baseline(
            cyc_samples,
            targets,
            cyclotomic,
            feature_index,
            train_ids,
            test_ids,
        )

        print(
            f"{name:18s} "
            f"train={train_acc:.6f} "
            f"test={test_acc:.6f}"
        )

    print()

    # ----------------------------------------------------------------------
    # 14. PER-TARGET TEST ACCURACY
    # ----------------------------------------------------------------------

    print("14. PER-TARGET TEST ACCURACY")
    print("-" * 78)

    model = LogisticModel()

    if cyc_train and cyc_test:
        model.fit(
            [sample.x for sample in cyc_train],
            [sample.y for sample in cyc_train],
        )

        for target_id in sorted(
            test_ids
        ):
            subset = [
                sample
                for sample in cyc_test
                if sample.target_id == target_id
            ]

            if not subset:
                print(
                    f"target {target_id:3d}: "
                    f"no usable cells"
                )
                continue

            predictions = model.predict(
                [sample.x for sample in subset]
            )

            acc = accuracy(
                subset,
                predictions,
            )

            bal = balanced_accuracy(
                subset,
                predictions,
            )

            print(
                f"target {target_id:3d}: "
                f"n={len(subset):2d} "
                f"accuracy={acc:.6f} "
                f"balanced={bal:.6f}"
            )

    print()

    # ----------------------------------------------------------------------
    # 15. FINAL DIAGNOSTIC
    # ----------------------------------------------------------------------

    print("=" * 78)
    print("15. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        f"cyclotomic target holdout = "
        f"{cyc_test_acc:.6f}"
    )

    print(
        f"control target holdout    = "
        f"{ctrl_test_acc:.6f}"
    )

    print(
        f"target delta              = "
        f"{cyc_test_acc - ctrl_test_acc:+.6f}"
    )

    if cyc_loo_valid:
        print(
            f"cyclotomic LOO mean       = "
            f"{statistics.fmean(cyc_loo_valid):.6f}"
        )
    else:
        print(
            "cyclotomic LOO mean       = nan"
        )

    if ctrl_loo_valid:
        print(
            f"control LOO mean          = "
            f"{statistics.fmean(ctrl_loo_valid):.6f}"
        )
    else:
        print(
            "control LOO mean          = nan"
        )

    print(
        f"cyclotomic 2D             = "
        f"{cyc_2d_acc:.6f}"
    )

    print(
        f"control 2D                = "
        f"{ctrl_2d_acc:.6f}"
    )

    if (
        not math.isnan(cyc_2d_acc)
        and not math.isnan(ctrl_2d_acc)
    ):
        print(
            f"2D delta                  = "
            f"{cyc_2d_acc - ctrl_2d_acc:+.6f}"
        )
    else:
        print(
            "2D delta                  = nan"
        )

    print(
        f"global permutation p      = "
        f"{p_global:.6f}"
    )

    print(
        f"within-target permutation p = "
        f"{p_within:.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)

    print(
        "The key test is whether joint N-only information survives "
        "both target and modulus holdout."
    )

    print(
        "A training improvement alone is insufficient."
    )

    print(
        "A useful positive result should show:"
    )

    print(
        "  * target holdout above baseline and controls;"
    )
    print(
        "  * leave-one-modulus-out above chance;"
    )
    print(
        "  * 2D unseen-target + unseen-modulus above controls;"
    )
    print(
        "  * permutation p-values substantially below 0.05;"
    )
    print(
        "  * positive reduction in H(Y|joint signature)."
    )

    print()
    print(
        "If target holdout improves but 2D holdout and permutation "
        "tests fail, treat the effect as non-transferable."
    )

    print()
    print(
        "If all tests fail, this is strong evidence that combining "
        "the existing N-only local character information does not "
        "recover the missing orientation."
    )

    elapsed = (
        time.perf_counter()
        - start_total
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 67 COMPLETE")
    print(
        f"total runtime = {elapsed:.6f}s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()