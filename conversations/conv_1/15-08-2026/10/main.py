#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 68R2
HIGH-POWER BALANCED CROSS-MODULUS ORIENTATION TEST
N-ONLY RELATIVE LEGENDRE ORIENTATION
ROBUST ROOT CONSTRUCTION / CACHED LOCAL FEATURES
TARGET + MODULUS + 2D HOLDOUT
GLOBAL + WITHIN-TARGET PERMUTATION NULL
EXACT MUTUAL INFORMATION
NO CSV OUTPUT
NO SKLEARN
==============================================================================

PATCHES FROM 68R
-----------------
1. Removed invalid Legendre evaluation modulo ell_a * ell_b.
   That modulus is composite and is not a Legendre-symbol modulus.

2. Cross-modulus features now consist only of:
      - legitimate local Legendre features at ell_a
      - legitimate local Legendre features at ell_b
      - products/differences of those local features
      - cubic residue-class one-hot features
      - modulus-order indicator

3. Root discovery remains direct and independently validated.

4. No factor p, q, or true sum s is exposed to the predictor.

5. All calculations use Python's standard library only.
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

SEED = 6802

NUM_TARGETS = 2000

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

CONTROL = [
    673, 4561, 4759, 6211,
    7879, 7951, 8689, 9781
]

TRAIN_TARGET_FRACTION = 0.75

PERMUTATIONS_GLOBAL = 300
PERMUTATIONS_WITHIN_TARGET = 300

HELDOUT_CYCLO = {61, 127, 307}
HELDOUT_CONTROL = {6211, 7951}


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
class PairSample:
    target_id: int
    ell_a: int
    ell_b: int
    label: int
    feature: tuple[float, ...]


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> list[int]:
    if hi < 2 or lo > hi:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            count = ((hi - start) // p) + 1
            sieve[start:hi + 1:p] = b"\x00" * count

    return [
        x
        for x in range(max(lo, 2), hi + 1)
        if sieve[x]
    ]


# ============================================================================
# CYCLOTOMIC ROOTS
# ============================================================================

_ROOT_CACHE: dict[int, tuple[int, int]] = {}


def roots_cyclotomic(ell: int) -> tuple[int, int]:
    """
    Directly enumerate the roots of

        x^2 + x + 1 = 0 mod ell.

    This is deliberately simple because all experimental moduli are <= 9781.
    """

    if ell < 2:
        raise ValueError(f"invalid ell={ell}")

    if ell % 3 != 1:
        raise ValueError(
            f"ell={ell} does not satisfy ell == 1 mod 3"
        )

    cached = _ROOT_CACHE.get(ell)
    if cached is not None:
        return cached

    roots = [
        x
        for x in range(ell)
        if (x * x + x + 1) % ell == 0
    ]

    if len(roots) != 2:
        raise RuntimeError(
            f"Expected two roots for ell={ell}, found {roots}"
        )

    a, b = roots

    if a == b:
        raise RuntimeError(
            f"Repeated root for prime ell={ell}"
        )

    if (a + b) % ell != ell - 1:
        raise RuntimeError(
            f"Vieta sum failure for ell={ell}"
        )

    if (a * b) % ell != 1:
        raise RuntimeError(
            f"Vieta product failure for ell={ell}"
        )

    cached = (a, b)
    _ROOT_CACHE[ell] = cached

    return cached


def prepare_roots(moduli: Sequence[int]) -> None:
    for ell in sorted(set(moduli)):
        roots_cyclotomic(ell)


def validate_roots(moduli: Sequence[int]) -> None:
    print("4. MODULUS VALIDATION")
    print("-" * 78)

    failures = 0

    for ell in moduli:
        try:
            a, b = roots_cyclotomic(ell)

            ok = (
                (a * a + a + 1) % ell == 0
                and
                (b * b + b + 1) % ell == 0
                and
                (a + b) % ell == ell - 1
                and
                (a * b) % ell == 1
            )

            print(
                f"ell={ell:5d} "
                f"roots=({a:5d},{b:5d}) "
                f"identity_ok={ok}"
            )

            if not ok:
                failures += 1

        except Exception as exc:
            print(
                f"ell={ell:5d} ERROR={exc}"
            )
            failures += 1

    print(f"validation failures = {failures}")
    print(
        "status = PASS"
        if failures == 0
        else "status = FAIL"
    )
    print()

    if failures:
        raise RuntimeError("Modulus validation failed")


# ============================================================================
# LEGENDRE SYMBOL
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    """
    Legendre symbol for odd prime p.
    """

    if p <= 2:
        raise ValueError(
            f"Legendre modulus must be odd prime; got p={p}"
        )

    a %= p

    if a == 0:
        return 0

    v = pow(a, (p - 1) // 2, p)

    if v == 1:
        return 1

    if v == p - 1:
        return -1

    # This should never occur for prime p.
    raise RuntimeError(
        f"Invalid Legendre evaluation: a={a}, p={p}, result={v}"
    )


def chi(a: int, ell: int) -> int:
    return legendre_symbol(a % ell, ell)


# ============================================================================
# TRUE LOCAL ORIENTATION
# ============================================================================

def local_sigma(n: int, ell: int) -> int | None:
    """
    Construct the true common-sign orientation at ell.

    Returns:
        +1 / -1 if both local discriminants have the same nonzero
        Legendre sign.

        None if either is zero or the two signs differ.
    """

    w1, w2 = roots_cyclotomic(ell)

    nmod = n % ell

    d1 = (w1 * w1 - 4 * nmod) % ell
    d2 = (w2 * w2 - 4 * nmod) % ell

    c1 = legendre_symbol(d1, ell)
    c2 = legendre_symbol(d2, ell)

    if c1 == 0 or c2 == 0:
        return None

    if c1 != c2:
        return None

    return c1


# ============================================================================
# CUBIC CHARACTER SUPPORT
# ============================================================================

_FACTOR_CACHE: dict[int, set[int]] = {}
_PRIMITIVE_ROOT_CACHE: dict[int, int] = {}
_LOG_CACHE: dict[int, dict[int, int]] = {}


def factor_set(n: int) -> set[int]:
    cached = _FACTOR_CACHE.get(n)
    if cached is not None:
        return cached

    x = n
    factors: set[int] = set()

    if x % 2 == 0:
        factors.add(2)
        while x % 2 == 0:
            x //= 2

    d = 3

    while d * d <= x:
        if x % d == 0:
            factors.add(d)

            while x % d == 0:
                x //= d

        d += 2

    if x > 1:
        factors.add(x)

    _FACTOR_CACHE[n] = factors
    return factors


def primitive_root(p: int) -> int:
    cached = _PRIMITIVE_ROOT_CACHE.get(p)

    if cached is not None:
        return cached

    factors = factor_set(p - 1)

    for g in range(2, p):
        if all(
            pow(g, (p - 1) // q, p) != 1
            for q in factors
        ):
            _PRIMITIVE_ROOT_CACHE[p] = g
            return g

    raise RuntimeError(
        f"No primitive root found modulo {p}"
    )


def discrete_log_table(p: int) -> dict[int, int]:
    cached = _LOG_CACHE.get(p)

    if cached is not None:
        return cached

    g = primitive_root(p)

    table: dict[int, int] = {}
    x = 1

    for k in range(p - 1):
        table[x] = k
        x = (x * g) % p

    _LOG_CACHE[p] = table
    return table


def cubic_class(value: int, ell: int) -> int:
    """
    Return exponent class modulo 3.

    0 = divisible by ell or exponent == 0 mod 3
    1 = exponent == 1 mod 3
    2 = exponent == 2 mod 3
    """

    v = value % ell

    if v == 0:
        return 0

    table = discrete_log_table(ell)
    exponent = table[v]

    return exponent % 3


# ============================================================================
# N-ONLY LOCAL FEATURES
# ============================================================================

def local_features(n: int, ell: int) -> tuple[float, ...]:
    """
    Every quantity here is computable from n and ell only.
    """

    G = 16 * n * n + 4 * n + 1

    values = [
        n,
        n + 1,
        n - 1,
        4 * n + 1,
        G,
    ]

    result: list[float] = []

    # Quadratic characters.
    for v in values:
        result.append(
            float(chi(v, ell))
        )

    # Cubic classes, one-hot encoded.
    for v in (n, n + 1, 4 * n + 1):
        c = cubic_class(v, ell)

        result.extend([
            1.0 if c == 0 else 0.0,
            1.0 if c == 1 else 0.0,
            1.0 if c == 2 else 0.0,
        ])

    return tuple(result)


# ============================================================================
# CROSS-MODULUS FEATURES
# ============================================================================

def pair_features(
    n: int,
    ell_a: int,
    ell_b: int,
) -> tuple[float, ...]:
    """
    LEGAL N-ONLY cross-modulus feature construction.

    Crucially, every Legendre evaluation is still performed modulo a
    PRIME modulus ell_a or ell_b.

    We never evaluate a Legendre symbol modulo ell_a * ell_b.
    """

    fa = local_features(n, ell_a)
    fb = local_features(n, ell_b)

    out: list[float] = []

    # Original local features.
    out.extend(fa)
    out.extend(fb)

    # Pointwise products.
    for a, b in zip(fa, fb):
        out.append(a * b)

    # Pointwise differences.
    for a, b in zip(fa, fb):
        out.append(a - b)

    # Absolute differences.
    for a, b in zip(fa, fb):
        out.append(abs(a - b))

    # Explicit N-only quadratic-character relationships.
    # These remain valid because each chi is evaluated at its own prime.
    common_values = (
        n,
        n + 1,
        n - 1,
        4 * n + 1,
        16 * n * n + 4 * n + 1,
    )

    for value in common_values:
        ca = chi(value, ell_a)
        cb = chi(value, ell_b)

        out.append(float(ca))
        out.append(float(cb))
        out.append(float(ca * cb))
        out.append(float(ca - cb))
        out.append(float(ca == cb))

    # Modulus-independent comparison features.
    out.extend([
        1.0 if ell_a < ell_b else 0.0,
        math.log(ell_a),
        math.log(ell_b),
        math.log(ell_a * ell_b),
    ])

    return tuple(out)


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    rng: random.Random,
) -> list[Target]:

    targets: list[Target] = []
    seen_n: set[int] = set()

    while len(targets) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen_n:
            continue

        seen_n.add(n)

        targets.append(
            Target(
                idx=len(targets) + 1,
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return targets


# ============================================================================
# LOGISTIC MODEL
# ============================================================================

class LogisticModel:
    """
    Small standard-library logistic regression.

    This avoids sklearn warnings/deprecations entirely.
    """

    def __init__(
        self,
        learning_rate: float = 0.05,
        epochs: int = 500,
        l2: float = 0.8,
    ):
        self.learning_rate = learning_rate
        self.epochs = epochs
        self.l2 = l2

        self.means: list[float] = []
        self.scales: list[float] = []
        self.weights: list[float] = []
        self.bias = 0.0

    @staticmethod
    def sigmoid(z: float) -> float:
        if z >= 0.0:
            e = math.exp(-z)
            return 1.0 / (1.0 + e)

        e = math.exp(z)
        return e / (1.0 + e)

    def fit(
        self,
        X: Sequence[Sequence[float]],
        y: Sequence[int],
    ) -> None:

        if len(X) != len(y):
            raise ValueError("X/y mismatch")

        if not X:
            raise ValueError("empty training set")

        m = len(X)
        d = len(X[0])

        self.means = []
        self.scales = []

        for j in range(d):
            col = [
                float(X[i][j])
                for i in range(m)
            ]

            mean = statistics.fmean(col)

            var = statistics.fmean(
                (x - mean) ** 2
                for x in col
            )

            scale = math.sqrt(var)

            self.means.append(mean)
            self.scales.append(
                scale if scale > 1e-12 else 1.0
            )

        Z = [
            [
                (
                    float(X[i][j]) - self.means[j]
                ) / self.scales[j]
                for j in range(d)
            ]
            for i in range(m)
        ]

        self.weights = [0.0] * d
        self.bias = 0.0

        for _ in range(self.epochs):
            grad_w = [0.0] * d
            grad_b = 0.0

            for i in range(m):
                z = self.bias

                row = Z[i]

                for j in range(d):
                    z += self.weights[j] * row[j]

                p = self.sigmoid(z)
                err = p - float(y[i])

                grad_b += err

                for j in range(d):
                    grad_w[j] += err * row[j]

            inv_m = 1.0 / m

            grad_b *= inv_m

            for j in range(d):
                grad_w[j] = (
                    grad_w[j] * inv_m
                    + self.l2 * self.weights[j]
                )

            self.bias -= (
                self.learning_rate * grad_b
            )

            for j in range(d):
                self.weights[j] -= (
                    self.learning_rate * grad_w[j]
                )

    def predict_one(
        self,
        x: Sequence[float],
    ) -> int:

        z = self.bias

        for j, value in enumerate(x):
            z += self.weights[j] * (
                (float(value) - self.means[j])
                / self.scales[j]
            )

        return 1 if self.sigmoid(z) >= 0.5 else 0

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
    y_true: Sequence[int],
    y_pred: Sequence[int],
) -> float:

    if not y_true:
        return float("nan")

    return sum(
        a == b
        for a, b in zip(y_true, y_pred)
    ) / len(y_true)


def balanced_accuracy(
    y_true: Sequence[int],
    y_pred: Sequence[int],
) -> float:

    pos = sum(
        1
        for y in y_true
        if y == 1
    )

    neg = sum(
        1
        for y in y_true
        if y == 0
    )

    if pos == 0 or neg == 0:
        return float("nan")

    tp = sum(
        1
        for a, b in zip(y_true, y_pred)
        if a == 1 and b == 1
    )

    tn = sum(
        1
        for a, b in zip(y_true, y_pred)
        if a == 0 and b == 0
    )

    return 0.5 * (
        tp / pos
        + tn / neg
    )


def entropy(labels: Sequence[int]) -> float:
    if not labels:
        return 0.0

    n = len(labels)
    c = Counter(labels)

    h = 0.0

    for count in c.values():
        p = count / n
        h -= p * math.log2(p)

    return h


def mutual_information(
    X: Sequence[tuple],
    Y: Sequence[int],
) -> float:

    if not X:
        return 0.0

    n = len(X)

    joint = Counter(zip(X, Y))
    px = Counter(X)
    py = Counter(Y)

    mi = 0.0

    for (x, y), count in joint.items():
        pxy = count / n
        pxv = px[x] / n
        pyv = py[y] / n

        mi += pxy * math.log2(
            pxy / (pxv * pyv)
        )

    return mi


# ============================================================================
# DATASET
# ============================================================================

def build_pair_dataset(
    targets: Sequence[Target],
    moduli: Sequence[int],
    rng: random.Random,
) -> tuple[list[PairSample], int, int]:

    all_samples: list[PairSample] = []

    raw_pos = 0
    raw_neg = 0

    # Precompute every local sigma and every N-only feature exactly once
    # for each target/modulus pair.
    for t in targets:

        sigma: dict[int, int | None] = {}
        feature_cache: dict[int, tuple[float, ...]] = {}

        for ell in moduli:
            sigma[ell] = local_sigma(t.n, ell)
            feature_cache[ell] = local_features(
                t.n,
                ell,
            )

        for i in range(len(moduli)):
            ell_a = moduli[i]
            sa = sigma[ell_a]

            if sa is None:
                continue

            for j in range(i + 1, len(moduli)):
                ell_b = moduli[j]
                sb = sigma[ell_b]

                if sb is None:
                    continue

                label = 1 if sa == sb else 0

                if label:
                    raw_pos += 1
                else:
                    raw_neg += 1

                # Build the feature vector without recomputing local
                # Legendre/cubic features.
                fa = feature_cache[ell_a]
                fb = feature_cache[ell_b]

                f: list[float] = []

                f.extend(fa)
                f.extend(fb)

                f.extend(
                    a * b
                    for a, b in zip(fa, fb)
                )

                f.extend(
                    a - b
                    for a, b in zip(fa, fb)
                )

                f.extend(
                    abs(a - b)
                    for a, b in zip(fa, fb)
                )

                values = (
                    t.n,
                    t.n + 1,
                    t.n - 1,
                    4 * t.n + 1,
                    16 * t.n * t.n
                    + 4 * t.n
                    + 1,
                )

                for value in values:
                    ca = chi(value, ell_a)
                    cb = chi(value, ell_b)

                    f.extend([
                        float(ca),
                        float(cb),
                        float(ca * cb),
                        float(ca - cb),
                        float(ca == cb),
                    ])

                f.extend([
                    1.0 if ell_a < ell_b else 0.0,
                    math.log(ell_a),
                    math.log(ell_b),
                    math.log(ell_a * ell_b),
                ])

                all_samples.append(
                    PairSample(
                        target_id=t.idx,
                        ell_a=ell_a,
                        ell_b=ell_b,
                        label=label,
                        feature=tuple(f),
                    )
                )

    # Exact global class balance.
    positives = [
        s
        for s in all_samples
        if s.label == 1
    ]

    negatives = [
        s
        for s in all_samples
        if s.label == 0
    ]

    take = min(
        len(positives),
        len(negatives),
    )

    positives = rng.sample(
        positives,
        take,
    )

    negatives = rng.sample(
        negatives,
        take,
    )

    balanced = positives + negatives
    rng.shuffle(balanced)

    return (
        balanced,
        len(balanced),
        raw_pos + raw_neg,
    )


# ============================================================================
# TARGET SPLIT
# ============================================================================

def split_targets(
    targets: Sequence[Target],
    rng: random.Random,
) -> tuple[set[int], set[int]]:

    ids = [
        t.idx
        for t in targets
    ]

    rng.shuffle(ids)

    cut = int(
        len(ids) * TRAIN_TARGET_FRACTION
    )

    return (
        set(ids[:cut]),
        set(ids[cut:]),
    )


# ============================================================================
# SCORING
# ============================================================================

def train_and_score(
    train: Sequence[PairSample],
    test: Sequence[PairSample],
) -> tuple[float, float, float]:

    if not train or not test:
        return (
            float("nan"),
            float("nan"),
            float("nan"),
        )

    model = LogisticModel()

    X_train = [
        s.feature
        for s in train
    ]

    y_train = [
        s.label
        for s in train
    ]

    X_test = [
        s.feature
        for s in test
    ]

    y_test = [
        s.label
        for s in test
    ]

    model.fit(
        X_train,
        y_train,
    )

    tr_pred = model.predict(
        X_train
    )

    te_pred = model.predict(
        X_test
    )

    return (
        accuracy(
            y_train,
            tr_pred,
        ),
        accuracy(
            y_test,
            te_pred,
        ),
        balanced_accuracy(
            y_test,
            te_pred,
        ),
    )


# ============================================================================
# TARGET HOLDOUT
# ============================================================================

def target_holdout(
    samples: Sequence[PairSample],
    train_ids: set[int],
    test_ids: set[int],
) -> tuple[float, float, float]:

    train = [
        s
        for s in samples
        if s.target_id in train_ids
    ]

    test = [
        s
        for s in samples
        if s.target_id in test_ids
    ]

    _, test_acc, test_bal = train_and_score(
        train,
        test,
    )

    return (
        test_acc,
        test_bal,
        len(test),
    )


# ============================================================================
# LOO MODULUS HOLDOUT
# ============================================================================

def loo_modulus(
    samples: Sequence[PairSample],
    ell: int,
    train_ids: set[int],
    test_ids: set[int],
) -> tuple[float, float, int]:

    train = [
        s
        for s in samples
        if (
            s.target_id in train_ids
            and s.ell_a != ell
            and s.ell_b != ell
        )
    ]

    test = [
        s
        for s in samples
        if (
            s.target_id in test_ids
            and (
                s.ell_a == ell
                or s.ell_b == ell
            )
        )
    ]

    if not train or not test:
        return (
            float("nan"),
            float("nan"),
            len(test),
        )

    _, acc, bal = train_and_score(
        train,
        test,
    )

    return (
        acc,
        bal,
        len(test),
    )


# ============================================================================
# 2D HOLDOUT
# ============================================================================

def holdout_2d(
    samples: Sequence[PairSample],
    heldout: set[int],
    train_ids: set[int],
    test_ids: set[int],
) -> tuple[float, float, int]:

    train = [
        s
        for s in samples
        if (
            s.target_id in train_ids
            and s.ell_a not in heldout
            and s.ell_b not in heldout
        )
    ]

    test = [
        s
        for s in samples
        if (
            s.target_id in test_ids
            and s.ell_a in heldout
            and s.ell_b in heldout
        )
    ]

    if not train or not test:
        return (
            float("nan"),
            float("nan"),
            len(test),
        )

    _, acc, bal = train_and_score(
        train,
        test,
    )

    return (
        acc,
        bal,
        len(test),
    )


# ============================================================================
# PERMUTATION NULL
# ============================================================================

def permutation_pvalue(
    samples: Sequence[PairSample],
    train_ids: set[int],
    test_ids: set[int],
    observed: float,
    rng: random.Random,
    repetitions: int,
    within_target: bool,
) -> float:

    train = [
        s
        for s in samples
        if s.target_id in train_ids
    ]

    test = [
        s
        for s in samples
        if s.target_id in test_ids
    ]

    if not train or not test:
        return float("nan")

    if math.isnan(observed):
        return float("nan")

    exceed = 0

    # Precompute group indices if doing within-target permutation.
    groups: dict[int, list[int]] = defaultdict(list)

    if within_target:
        for i, s in enumerate(train):
            groups[s.target_id].append(i)

    original_labels = [
        s.label
        for s in train
    ]

    for _ in range(repetitions):

        labels = original_labels[:]

        if within_target:
            for indices in groups.values():
                vals = [
                    labels[i]
                    for i in indices
                ]

                rng.shuffle(vals)

                for i, v in zip(indices, vals):
                    labels[i] = v
        else:
            rng.shuffle(labels)

        perm_train = [
            PairSample(
                target_id=s.target_id,
                ell_a=s.ell_a,
                ell_b=s.ell_b,
                label=labels[i],
                feature=s.feature,
            )
            for i, s in enumerate(train)
        ]

        _, acc, _ = train_and_score(
            perm_train,
            test,
        )

        if not math.isnan(acc) and acc >= observed:
            exceed += 1

    return (
        exceed + 1
    ) / (
        repetitions + 1
    )


# ============================================================================
# BASELINE
# ============================================================================

def baseline_from_feature(
    samples: Sequence[PairSample],
    train_ids: set[int],
    test_ids: set[int],
    index: int,
) -> float:

    train = [
        s
        for s in samples
        if s.target_id in train_ids
    ]

    test = [
        s
        for s in samples
        if s.target_id in test_ids
    ]

    if not train or not test:
        return float("nan")

    score = sum(
        (
            1
            if s.label == 1
            else -1
        )
        * s.feature[index]
        for s in train
    )

    direction = (
        1
        if score >= 0
        else -1
    )

    pred = [
        1
        if direction * s.feature[index] >= 0
        else 0
        for s in test
    ]

    return accuracy(
        [s.label for s in test],
        pred,
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    total_start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 68R2")
    print("HIGH-POWER BALANCED CROSS-MODULUS ORIENTATION TEST")
    print("N-ONLY RELATIVE LEGENDRE ORIENTATION")
    print("ROBUST ROOT CONSTRUCTION / CACHED LOCAL FEATURES")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("GLOBAL + WITHIN-TARGET PERMUTATION NULL")
    print("EXACT MUTUAL INFORMATION")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    # ----------------------------------------------------------------------
    # PRIME POPULATION
    # ----------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_MIN,
        PRIME_MAX,
    )

    print("1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )
    print()

    # ----------------------------------------------------------------------
    # TARGETS
    # ----------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
        rng,
    )

    print("2. TARGET SUMMARY")
    print("-" * 78)
    print(
        f"total targets = {len(targets)}"
    )

    for t in targets[:24]:
        print(
            f"target {t.idx:4d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    print(
        "... remaining generated targets omitted"
    )
    print()

    # ----------------------------------------------------------------------
    # MODULI
    # ----------------------------------------------------------------------

    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(
        f"cyclotomic        = {CYCLOTOMIC}"
    )
    print(
        f"valid controls    = {CONTROL}"
    )
    print(
        "excluded controls = "
        "[1109, 2309, 2351, 4421, 8273]"
    )
    print()

    # ----------------------------------------------------------------------
    # ROOT VALIDATION
    # ----------------------------------------------------------------------

    prepare_roots(
        CYCLOTOMIC + CONTROL
    )

    validate_roots(
        CYCLOTOMIC + CONTROL
    )

    # ----------------------------------------------------------------------
    # HOLDOUT
    # ----------------------------------------------------------------------

    train_ids, test_ids = split_targets(
        targets,
        rng,
    )

    print("5. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train_ids)}"
    )
    print(
        f"test targets     = {len(test_ids)}"
    )
    print()

    # ----------------------------------------------------------------------
    # DATASETS
    # ----------------------------------------------------------------------

    t0 = time.perf_counter()

    cyclo_samples, cyclo_count, cyclo_raw = (
        build_pair_dataset(
            targets,
            CYCLOTOMIC,
            rng,
        )
    )

    control_samples, control_count, control_raw = (
        build_pair_dataset(
            targets,
            CONTROL,
            rng,
        )
    )

    print("6. CROSS-MODULUS DATASET")
    print("-" * 78)
    print(
        f"cyclotomic raw samples = "
        f"{cyclo_raw}"
    )
    print(
        f"cyclotomic balanced   = "
        f"{cyclo_count}"
    )
    print(
        f"control raw samples   = "
        f"{control_raw}"
    )
    print(
        f"control balanced      = "
        f"{control_count}"
    )
    print(
        f"construction time = "
        f"{time.perf_counter() - t0:.6f}s"
    )
    print()

    # ----------------------------------------------------------------------
    # LABEL BALANCE
    # ----------------------------------------------------------------------

    def report_balance(
        name: str,
        samples: Sequence[PairSample],
    ) -> None:

        pos = sum(
            s.label
            for s in samples
        )

        neg = len(samples) - pos

        if samples:
            majority = (
                max(pos, neg)
                / len(samples)
            )
        else:
            majority = float("nan")

        print(
            f"{name:12s} "
            f"positive={pos:7d} "
            f"negative={neg:7d} "
            f"majority={majority:.6f}"
        )

    print("7. LABEL DISTRIBUTIONS")
    print("-" * 78)

    report_balance(
        "cyclotomic",
        cyclo_samples,
    )

    report_balance(
        "control",
        control_samples,
    )

    print()

    # ----------------------------------------------------------------------
    # TARGET HOLDOUT
    # ----------------------------------------------------------------------

    print("8. TARGET HOLDOUT")
    print("-" * 78)

    ctest, cbal, cn = target_holdout(
        cyclo_samples,
        train_ids,
        test_ids,
    )

    rtest, rbal, rn = target_holdout(
        control_samples,
        train_ids,
        test_ids,
    )

    print(
        f"cyclotomic test = {ctest:.6f} "
        f"balanced={cbal:.6f} "
        f"n={cn}"
    )

    print(
        f"control    test = {rtest:.6f} "
        f"balanced={rbal:.6f} "
        f"n={rn}"
    )

    if not math.isnan(ctest) and not math.isnan(rtest):
        print(
            f"test delta = "
            f"{ctest-rtest:+.6f}"
        )

    print()

    # ----------------------------------------------------------------------
    # EXACT INFORMATION
    # ----------------------------------------------------------------------

    print("9. EXACT JOINT SIGNATURE INFORMATION")
    print("-" * 78)

    def signature(s: PairSample) -> tuple:
        # Include local chi(n), chi(n+1), chi(G)
        # from both moduli.
        return (
            s.feature[0],
            s.feature[1],
            s.feature[4],
            s.feature[13],
            s.feature[14],
            s.feature[17],
        )

    for name, samples in (
        ("cyclotomic", cyclo_samples),
        ("control", control_samples),
    ):

        X = [
            signature(s)
            for s in samples
        ]

        Y = [
            s.label
            for s in samples
        ]

        mi = mutual_information(
            X,
            Y,
        )

        hy = entropy(Y)

        hcond = hy - mi

        print(
            f"{name:12s} "
            f"H(Y)={hy:.6f} "
            f"H(Y|X)={hcond:.6f} "
            f"I={mi:.6f} bits"
        )

    print()

    # ----------------------------------------------------------------------
    # SIMPLE BASELINES
    # ----------------------------------------------------------------------

    print("10. SIMPLE N-ONLY BASELINES")
    print("-" * 78)

    baseline_map = {
        0: "chi_a(n)",
        1: "chi_a(n+1)",
        2: "chi_a(n-1)",
        3: "chi_a(4n+1)",
        4: "chi_a(G)",
    }

    for idx, name in baseline_map.items():

        c = baseline_from_feature(
            cyclo_samples,
            train_ids,
            test_ids,
            idx,
        )

        r = baseline_from_feature(
            control_samples,
            train_ids,
            test_ids,
            idx,
        )

        print(
            f"{name:15s} "
            f"C={c:.6f} "
            f"R={r:.6f} "
            f"delta={c-r:+.6f}"
        )

    print()

    # ----------------------------------------------------------------------
    # LOO
    # ----------------------------------------------------------------------

    print("11. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    cyclo_loo: list[float] = []
    control_loo: list[float] = []

    for ell in CYCLOTOMIC[1:]:

        acc, bal, n = loo_modulus(
            cyclo_samples,
            ell,
            train_ids,
            test_ids,
        )

        if not math.isnan(acc):
            cyclo_loo.append(acc)

        print(
            f"cyclo ell={ell:5d} "
            f"accuracy={acc:.6f} "
            f"balanced={bal:.6f} "
            f"n={n}"
        )

    for ell in CONTROL[1:]:

        acc, bal, n = loo_modulus(
            control_samples,
            ell,
            train_ids,
            test_ids,
        )

        if not math.isnan(acc):
            control_loo.append(acc)

        print(
            f"control ell={ell:5d} "
            f"accuracy={acc:.6f} "
            f"balanced={bal:.6f} "
            f"n={n}"
        )

    cmean = (
        statistics.fmean(cyclo_loo)
        if cyclo_loo
        else float("nan")
    )

    rmean = (
        statistics.fmean(control_loo)
        if control_loo
        else float("nan")
    )

    print(
        f"cyclotomic LOO mean = "
        f"{cmean:.6f}"
    )

    print(
        f"control LOO mean    = "
        f"{rmean:.6f}"
    )

    print()

    # ----------------------------------------------------------------------
    # 2D HOLDOUT
    # ----------------------------------------------------------------------

    c2d, c2db, c2dn = holdout_2d(
        cyclo_samples,
        HELDOUT_CYCLO,
        train_ids,
        test_ids,
    )

    r2d, r2db, r2dn = holdout_2d(
        control_samples,
        HELDOUT_CONTROL,
        train_ids,
        test_ids,
    )

    print("12. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)

    print(
        f"cyclotomic 2D = {c2d:.6f} "
        f"balanced={c2db:.6f} "
        f"n={c2dn} "
        f"heldout={sorted(HELDOUT_CYCLO)}"
    )

    print(
        f"control 2D    = {r2d:.6f} "
        f"balanced={r2db:.6f} "
        f"n={r2dn} "
        f"heldout={sorted(HELDOUT_CONTROL)}"
    )

    if not math.isnan(c2d) and not math.isnan(r2d):
        print(
            f"2D delta = "
            f"{c2d-r2d:+.6f}"
        )

    print()

    # ----------------------------------------------------------------------
    # PERMUTATIONS
    # ----------------------------------------------------------------------

    print("13. PERMUTATION NULL")
    print("-" * 78)

    pg = permutation_pvalue(
        cyclo_samples,
        train_ids,
        test_ids,
        ctest,
        rng,
        PERMUTATIONS_GLOBAL,
        False,
    )

    pw = permutation_pvalue(
        cyclo_samples,
        train_ids,
        test_ids,
        ctest,
        rng,
        PERMUTATIONS_WITHIN_TARGET,
        True,
    )

    print(
        f"observed cyclotomic accuracy = "
        f"{ctest:.6f}"
    )

    print(
        f"global-label permutation p = "
        f"{pg:.6f}"
    )

    print(
        f"within-target permutation p = "
        f"{pw:.6f}"
    )

    print()

    # ----------------------------------------------------------------------
    # FINAL
    # ----------------------------------------------------------------------

    print("14. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        f"cyclotomic target holdout = "
        f"{ctest:.6f}"
    )

    print(
        f"control target holdout    = "
        f"{rtest:.6f}"
    )

    if not math.isnan(ctest) and not math.isnan(rtest):
        print(
            f"target delta = "
            f"{ctest-rtest:+.6f}"
        )

    print(
        f"cyclotomic LOO mean = "
        f"{cmean:.6f}"
    )

    print(
        f"control LOO mean = "
        f"{rmean:.6f}"
    )

    print(
        f"cyclotomic 2D = "
        f"{c2d:.6f}"
    )

    print(
        f"control 2D = "
        f"{r2d:.6f}"
    )

    print(
        f"global permutation p = "
        f"{pg:.6f}"
    )

    print(
        f"within-target permutation p = "
        f"{pw:.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)
    print(
        """
The previous crash was caused by an invalid composite modulus entering
the Legendre-symbol routine.

This version keeps every Legendre evaluation local to an odd prime ell.

A useful result should require simultaneous evidence from:

  1. target holdout;
  2. balanced target holdout;
  3. leave-one-modulus-out transfer;
  4. two-dimensional unseen-target/unseen-modulus transfer;
  5. low permutation p-values;
  6. positive exact mutual information.

A target-holdout advantage alone is not sufficient.

The 68R2 feature construction contains no p, q, s, local roots,
discriminant values, or true orientation labels.
"""
    )

    elapsed = time.perf_counter() - total_start

    print("=" * 78)
    print("EXPERIMENT 68R2 COMPLETE")
    print(
        f"total runtime = {elapsed:.6f}s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()