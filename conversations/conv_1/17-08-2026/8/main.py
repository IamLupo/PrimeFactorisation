#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 78
COMPRESSED MODULUS-NORMALIZED CHARACTER INVARIANTS
E1 CHARACTER PREDICTION WITHOUT EXACT STATE LOOKUP
TARGET + MODULUS + 2D HOLDOUT
3-CLASS E1 CHARACTER + NONZERO SIGN TASK
LEAVE-ONE-MODULUS-OUT
WITHIN-TARGET PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Main question
-------------
Do small, reusable, modulus-normalized summaries of N-only character data
carry transferable information about the TRUE E1 character?

This experiment deliberately avoids raw tuples such as

    (chi_7(n), chi_13(n), ..., cubic_1723(...))

and instead compresses them into counts / aggregate statistics:

    Legendre class counts over OTHER moduli
    cubic class counts over OTHER moduli
    aggregate character sums

The goal is to distinguish genuine reusable structure from exact-state
collision / memorisation effects seen in Experiments 76-77.

Tasks
-----
A. Full 3-class prediction:
       chi_ell(E1) in {-1, 0, +1}

B. Nonzero sign prediction:
       chi_ell(E1) in {-1, +1}, restricted to chi != 0

C. Target holdout
D. Leave-one-modulus-out
E. Two-dimensional target + modulus holdout
F. Within-target permutation null
G. Occupancy / class-balance diagnostics

Expected runtime: seconds, not minutes.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 78001

TARGETS_TOTAL = 120
TRAIN_TARGET_FRACTION = 0.75

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

CONTROLS = [
    673, 4561, 4759, 6211, 7879, 7951, 8689, 9781
]

# Polynomials evaluated on N for the compressed signature.
POLYS = (
    "n-1",
    "n",
    "n+1",
    "n2+n+1",
)

# Cubic-character inputs.
CUBIC_POLYS = (
    "n",
    "n+1",
    "4n+1",
)

PERMUTATIONS = 200

# Two-dimensional holdout.
HOLDOUT_CYCLO_MODULI = [61, 127, 307]
HOLDOUT_CONTROL_MODULI = [6211, 7951]


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q

    @property
    def e1(self) -> int:
        # E1 = (p^2 - 1)(q^2 - 1)(p + q)
        #    = s * ((n + 1)^2 - s^2)
        return self.s * ((self.n + 1) ** 2 - self.s ** 2)


@dataclass(frozen=True)
class Sample:
    target_idx: int
    ell: int
    x: Tuple[float, ...]
    y3: int
    ysign: int | None


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi < 2 or lo > hi:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))
    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start : hi + 1 : p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [p for p in range(max(2, lo), hi + 1) if sieve[p]]


def make_targets(primes: Sequence[int], count: int, rng: random.Random) -> List[Target]:
    out: List[Target] = []
    used: set[Tuple[int, int]] = set()

    while len(out) < count:
        p, q = rng.sample(primes, 2)
        if p > q:
            p, q = q, p

        key = (p, q)
        if key in used:
            continue

        # Keep the same general semiprime scale as the earlier experiments.
        if not (p > 2_000_000 and q > 2_000_000):
            continue

        used.add(key)
        out.append(Target(len(out) + 1, p, q))

    return out


# ============================================================================
# FINITE-FIELD / CHARACTER UTILITIES
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    r = pow(a, (p - 1) // 2, p)

    if r == 1:
        return 1
    if r == p - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre result a={a}, p={p}, pow={r}"
    )


def cube_root_pair_prime(p: int) -> Tuple[int, int]:
    """
    Roots of x^2 + x + 1 = 0 mod p.

    Valid for p == 1 mod 3.
    """
    if p == 3 or p % 3 != 1:
        raise ValueError(f"{p} is not valid for cubic roots")

    # Search a non-trivial cube root.
    # The moduli are tiny, so this is intentionally simple and robust.
    for z in range(2, p):
        if (z * z + z + 1) % p == 0:
            z2 = (-1 - z) % p

            if z2 != z and (z2 * z2 + z2 + 1) % p == 0:
                return min(z, z2), max(z, z2)

    raise RuntimeError(f"No cubic roots found for p={p}")


def cubic_class(a: int, p: int) -> int:
    """
    Multiplicative cubic-character class.

    Returns:
        0 if a == 0 mod p
        1 or 2 otherwise, according to a^(p-1)/3.

    This deliberately uses a class label rather than attempting to
    identify an absolute cubic root orientation.
    """
    a %= p

    if a == 0:
        return 0

    r = pow(a, (p - 1) // 3, p)

    # Because p == 1 mod 3, nonzero cubic character values are
    # 1, zeta, zeta^2. We collapse them to class 0/1/2 relative
    # to a canonical primitive cube root.
    if r == 1:
        return 0

    z, z2 = cube_root_pair_prime(p)

    if r == z:
        return 1

    if r == z2:
        return 2

    # Depending on representative convention r may be the other
    # nontrivial root, but that is precisely z or z^2.
    raise RuntimeError(
        f"Unexpected cubic character a={a}, p={p}, r={r}, roots=({z},{z2})"
    )


def eval_poly(name: str, n: int) -> int:
    if name == "n-1":
        return n - 1
    if name == "n":
        return n
    if name == "n+1":
        return n + 1
    if name == "n2+n+1":
        return n * n + n + 1

    if name == "4n+1":
        return 4 * n + 1

    raise ValueError(name)


# ============================================================================
# CHARACTER TABLES
# ============================================================================

def prepare_root_table(moduli: Sequence[int]) -> Dict[int, Tuple[int, int]]:
    roots: Dict[int, Tuple[int, int]] = {}

    for ell in moduli:
        r1, r2 = cube_root_pair_prime(ell)

        if (r1 * r1 + r1 + 1) % ell != 0:
            raise RuntimeError(f"Root validation failed for ell={ell}")

        if (r2 * r2 + r2 + 1) % ell != 0:
            raise RuntimeError(f"Root validation failed for ell={ell}")

        roots[ell] = (r1, r2)

    return roots


def character_signature(
    n: int,
    moduli: Sequence[int],
    exclude_ell: int,
) -> Tuple[int, ...]:
    """
    Compressed, modulus-normalized N-only Legendre signature.

    For every OTHER modulus:
      four polynomials -> Legendre classes

    Then append aggregate sums / zero counts.

    This is intentionally much smaller than the exact raw state.
    """
    values: List[int] = []

    # Class counts for each polynomial.
    counts = {
        poly: {-1: 0, 0: 0, 1: 0}
        for poly in POLYS
    }

    sums = {
        poly: 0
        for poly in POLYS
    }

    for m in moduli:
        if m == exclude_ell:
            continue

        for poly in POLYS:
            c = legendre_symbol(eval_poly(poly, n), m)
            counts[poly][c] += 1
            sums[poly] += c

    # Ordered compressed counts.
    for poly in POLYS:
        values.extend([
            counts[poly][-1],
            counts[poly][0],
            counts[poly][1],
        ])

    # Aggregate character sums.
    for poly in POLYS:
        values.append(sums[poly])

    return tuple(values)


def cubic_signature(
    n: int,
    moduli: Sequence[int],
    exclude_ell: int,
) -> Tuple[int, ...]:
    values: List[int] = []

    counts = {
        poly: {0: 0, 1: 0, 2: 0}
        for poly in CUBIC_POLYS
    }

    for m in moduli:
        if m == exclude_ell:
            continue

        for poly in CUBIC_POLYS:
            c = cubic_class(eval_poly(poly, n), m)
            counts[poly][c] += 1

    for poly in CUBIC_POLYS:
        values.extend([
            counts[poly][0],
            counts[poly][1],
            counts[poly][2],
        ])

    return tuple(values)


def combined_feature(
    n: int,
    family: Sequence[int],
    ell: int,
) -> Tuple[float, ...]:
    """
    Small reusable N-only invariant vector.

    Includes:
      * compressed Legendre counts
      * Legendre aggregate sums
      * compressed cubic counts
    """
    a = character_signature(n, family, ell)
    b = cubic_signature(n, family, ell)

    return tuple(float(x) for x in (*a, *b))


# ============================================================================
# LABELS
# ============================================================================

def e1_label(target: Target, ell: int) -> int:
    return legendre_symbol(target.e1, ell)


def e1_sign_label(target: Target, ell: int) -> int | None:
    y = e1_label(target, ell)
    if y == 0:
        return None
    return y


# ============================================================================
# DATASET CONSTRUCTION
# ============================================================================

def build_dataset(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> List[Sample]:
    out: List[Sample] = []

    for t in targets:
        for ell in moduli:
            y3 = e1_label(t, ell)
            ysign = e1_sign_label(t, ell)

            # Require at least two other moduli so the compressed
            # signature represents a genuine cross-modulus state.
            if len(moduli) <= 2:
                raise RuntimeError("Need at least 3 moduli")

            x = combined_feature(t.n, moduli, ell)

            out.append(
                Sample(
                    target_idx=t.idx,
                    ell=ell,
                    x=x,
                    y3=y3,
                    ysign=ysign,
                )
            )

    return out


# ============================================================================
# SIMPLE MODULUS-INVARIANT CLASSIFIER
# ============================================================================

class NearestCentroid:
    """
    Simple class centroid model.

    No sklearn, no hyperparameter tuning.
    """

    def __init__(self) -> None:
        self.centroids: Dict[int, Tuple[float, ...]] = {}
        self.mean: List[float] = []
        self.scale: List[float] = []

    def fit(self, samples: Sequence[Sample], label_attr: str) -> None:
        if not samples:
            raise ValueError("No training samples")

        xs = [s.x for s in samples]
        ys = [getattr(s, label_attr) for s in samples]

        if any(y is None for y in ys):
            raise ValueError("None labels not allowed")

        d = len(xs[0])

        self.mean = [
            statistics.fmean(x[i] for x in xs)
            for i in range(d)
        ]

        self.scale = []

        for i in range(d):
            vals = [x[i] for x in xs]
            mu = self.mean[i]

            var = statistics.fmean((v - mu) ** 2 for v in vals)
            sd = math.sqrt(var)

            self.scale.append(sd if sd > 1e-12 else 1.0)

        groups: Dict[int, List[List[float]]] = defaultdict(list)

        for x, y in zip(xs, ys):
            z = [
                (x[i] - self.mean[i]) / self.scale[i]
                for i in range(d)
            ]
            groups[int(y)].append(z)

        self.centroids = {}

        for y, rows in groups.items():
            self.centroids[y] = tuple(
                statistics.fmean(row[i] for row in rows)
                for i in range(d)
            )

    def predict_one(self, x: Tuple[float, ...]) -> int:
        if not self.centroids:
            raise RuntimeError("Model not fitted")

        z = [
            (x[i] - self.mean[i]) / self.scale[i]
            for i in range(len(x))
        ]

        best_y = None
        best_dist = float("inf")

        for y, c in self.centroids.items():
            dist = sum(
                (z[i] - c[i]) ** 2
                for i in range(len(z))
            )

            if dist < best_dist:
                best_dist = dist
                best_y = y

        assert best_y is not None
        return best_y

    def predict(self, samples: Sequence[Sample]) -> List[int]:
        return [self.predict_one(s.x) for s in samples]


# ============================================================================
# METRICS
# ============================================================================

def accuracy(y_true: Sequence[int], y_pred: Sequence[int]) -> float:
    if not y_true:
        return float("nan")

    return sum(a == b for a, b in zip(y_true, y_pred)) / len(y_true)


def balanced_accuracy(
    y_true: Sequence[int],
    y_pred: Sequence[int],
) -> float:
    classes = sorted(set(y_true))

    recalls = []

    for c in classes:
        idx = [i for i, y in enumerate(y_true) if y == c]

        if not idx:
            continue

        recalls.append(
            sum(y_pred[i] == c for i in idx) / len(idx)
        )

    if not recalls:
        return float("nan")

    return statistics.fmean(recalls)


def majority_baseline(y_true: Sequence[int]) -> float:
    if not y_true:
        return float("nan")

    counts = Counter(y_true)
    return max(counts.values()) / len(y_true)


def mutual_information(
    x: Sequence[Tuple[float, ...]],
    y: Sequence[int],
) -> float:
    """
    Exact empirical mutual information for a discrete feature tuple.

    This is diagnostic only; it is susceptible to sparse states.
    The experiment's main metric is OOS prediction.
    """
    if not x:
        return 0.0

    joint = Counter(zip(x, y))
    px = Counter(x)
    py = Counter(y)
    n = len(y)

    mi = 0.0

    for (xx, yy), cxy in joint.items():
        pxy = cxy / n
        pxv = px[xx] / n
        pyv = py[yy] / n

        mi += pxy * math.log2(pxy / (pxv * pyv))

    return mi


# ============================================================================
# DATA SPLITS
# ============================================================================

def split_targets(
    targets: Sequence[Target],
    rng: random.Random,
) -> Tuple[List[Target], List[Target]]:
    idxs = list(range(len(targets)))
    rng.shuffle(idxs)

    cut = int(len(idxs) * TRAIN_TARGET_FRACTION)

    train_ids = set(idxs[:cut])
    train = [t for i, t in enumerate(targets) if i in train_ids]
    test = [t for i, t in enumerate(targets) if i not in train_ids]

    return train, test


def filter_by_targets(
    samples: Sequence[Sample],
    targets: Sequence[Target],
) -> List[Sample]:
    wanted = {t.idx for t in targets}
    return [s for s in samples if s.target_idx in wanted]


# ============================================================================
# WITHIN-TARGET PERMUTATION
# ============================================================================

def permutation_test(
    train: Sequence[Sample],
    test: Sequence[Sample],
    label_attr: str,
    observed: float,
    rng: random.Random,
    permutations: int,
) -> float:
    grouped: Dict[int, List[int]] = defaultdict(list)

    y = [getattr(s, label_attr) for s in train]

    for i, s in enumerate(train):
        grouped[s.target_idx].append(i)

    exceed = 0

    for _ in range(permutations):
        yp = list(y)

        for idxs in grouped.values():
            vals = [yp[i] for i in idxs]
            rng.shuffle(vals)

            for j, i in enumerate(idxs):
                yp[i] = vals[j]

        perm_samples: List[Sample] = []

        for s, yy in zip(train, yp):
            # Preserve only features; replace the label.
            perm_samples.append(
                Sample(
                    target_idx=s.target_idx,
                    ell=s.ell,
                    x=s.x,
                    y3=int(yy),
                    ysign=s.ysign,
                )
            )

        model = NearestCentroid()
        model.fit(perm_samples, label_attr)

        preds = model.predict(test)
        yt = [getattr(s, label_attr) for s in test]

        score = balanced_accuracy(yt, preds)

        if score >= observed:
            exceed += 1

    return (exceed + 1) / (permutations + 1)


# ============================================================================
# MODULUS HOLDOUT
# ============================================================================

def leave_one_modulus_out(
    samples: Sequence[Sample],
    family: Sequence[int],
    label_attr: str,
) -> List[Tuple[int, float, int]]:
    results = []

    for ell in family:
        train = [s for s in samples if s.ell != ell]
        test = [s for s in samples if s.ell == ell]

        if not train or not test:
            results.append((ell, float("nan"), len(test)))
            continue

        # Only evaluate classes that actually exist in training.
        train_labels = {getattr(s, label_attr) for s in train}
        test = [
            s for s in test
            if getattr(s, label_attr) in train_labels
        ]

        if not test:
            results.append((ell, float("nan"), 0))
            continue

        model = NearestCentroid()
        model.fit(train, label_attr)

        preds = model.predict(test)
        yt = [getattr(s, label_attr) for s in test]

        results.append(
            (ell, balanced_accuracy(yt, preds), len(test))
        )

    return results


# ============================================================================
# 2D HOLDOUT
# ============================================================================

def two_dimensional_holdout(
    samples: Sequence[Sample],
    train_targets: Sequence[Target],
    holdout_moduli: Sequence[int],
    label_attr: str,
) -> Tuple[float, int]:
    train_target_ids = {t.idx for t in train_targets}

    train = [
        s for s in samples
        if s.target_idx in train_target_ids
        and s.ell not in holdout_moduli
    ]

    test = [
        s for s in samples
        if s.target_idx not in train_target_ids
        and s.ell in holdout_moduli
    ]

    if not train or not test:
        return float("nan"), 0

    train_labels = {getattr(s, label_attr) for s in train}

    test = [
        s for s in test
        if getattr(s, label_attr) in train_labels
    ]

    if not test:
        return float("nan"), 0

    model = NearestCentroid()
    model.fit(train, label_attr)

    pred = model.predict(test)
    yt = [getattr(s, label_attr) for s in test]

    return balanced_accuracy(yt, pred), len(test)


# ============================================================================
# REPORTING
# ============================================================================

def report_family(
    name: str,
    samples: Sequence[Sample],
    train_targets: Sequence[Target],
    test_targets: Sequence[Target],
    family: Sequence[int],
    rng: random.Random,
) -> None:
    print()
    print("=" * 78)
    print(name)
    print("=" * 78)

    train = filter_by_targets(samples, train_targets)
    test = filter_by_targets(samples, test_targets)

    # ----------------------------------------------------------------------
    # Full 3-class E1 character
    # ----------------------------------------------------------------------

    model = NearestCentroid()
    model.fit(train, "y3")

    pred = model.predict(test)

    yt = [s.y3 for s in test]

    acc = accuracy(yt, pred)
    bal = balanced_accuracy(yt, pred)
    base = majority_baseline(yt)

    mi_train = mutual_information(
        [s.x for s in train],
        [s.y3 for s in train],
    )

    print("3-CLASS E1 CHARACTER")
    print(f"  train samples        = {len(train)}")
    print(f"  test samples         = {len(test)}")
    print(f"  train states         = {len(set(s.x for s in train))}")
    print(f"  test states          = {len(set(s.x for s in test))}")
    print(f"  majority baseline    = {base:.6f}")
    print(f"  OOS accuracy         = {acc:.6f}")
    print(f"  OOS balanced         = {bal:.6f}")
    print(f"  train state MI       = {mi_train:.6f} bits")

    p3 = permutation_test(
        train,
        test,
        "y3",
        bal,
        rng,
        PERMUTATIONS,
    )

    print(f"  within-target permutation p = {p3:.6f}")

    # ----------------------------------------------------------------------
    # Nonzero sign task
    # ----------------------------------------------------------------------

    train_sign = [s for s in train if s.ysign is not None]
    test_sign = [s for s in test if s.ysign is not None]

    sign_model = None

    if train_sign and test_sign:
        sign_model = NearestCentroid()
        sign_model.fit(train_sign, "ysign")

        pred_s = sign_model.predict(test_sign)
        yt_s = [s.ysign for s in test_sign]

        s_acc = accuracy(yt_s, pred_s)
        s_bal = balanced_accuracy(yt_s, pred_s)
        s_base = majority_baseline(yt_s)

        p_sign = permutation_test(
            train_sign,
            test_sign,
            "ysign",
            s_bal,
            rng,
            PERMUTATIONS,
        )
    else:
        s_acc = s_bal = s_base = p_sign = float("nan")

    print()
    print("NONZERO SIGN TASK")
    print(f"  train nonzero = {len(train_sign)}")
    print(f"  test nonzero  = {len(test_sign)}")
    print(f"  majority      = {s_base:.6f}")
    print(f"  OOS accuracy  = {s_acc:.6f}")
    print(f"  OOS balanced  = {s_bal:.6f}")
    print(f"  permutation p = {p_sign:.6f}")

    # ----------------------------------------------------------------------
    # LOO modulus
    # ----------------------------------------------------------------------

    loo = leave_one_modulus_out(samples, family, "y3")
    valid = [a for _, a, n in loo if not math.isnan(a) and n > 0]

    print()
    print("LEAVE-ONE-MODULUS-OUT")
    for ell, score, n in loo:
        print(
            f"  ell={ell:5d} balanced={score!s:>10} n={n:4d}"
        )

    loo_mean = statistics.fmean(valid) if valid else float("nan")
    print(f"  mean LOO balanced = {loo_mean:.6f}")

    # ----------------------------------------------------------------------
    # 2D holdout
    # ----------------------------------------------------------------------

    if name.upper().startswith("CYC"):
        holdout = HOLDOUT_CYCLO_MODULI
    else:
        holdout = HOLDOUT_CONTROL_MODULI

    score_2d, n2d = two_dimensional_holdout(
        samples,
        train_targets,
        holdout,
        "y3",
    )

    print()
    print("2D TARGET + MODULUS HOLDOUT")
    print(f"  held-out moduli = {holdout}")
    print(f"  balanced accuracy = {score_2d:.6f}")
    print(f"  evaluated samples = {n2d}")

    # ----------------------------------------------------------------------
    # State occupancy
    # ----------------------------------------------------------------------

    occ = Counter(s.x for s in train)

    print()
    print("COMPRESSED STATE OCCUPANCY")
    print(f"  unique states     = {len(occ)}")
    print(f"  mean occupancy    = {statistics.fmean(occ.values()):.3f}")
    print(f"  median occupancy  = {statistics.median(occ.values()):.3f}")
    print(f"  max occupancy     = {max(occ.values())}")

    for threshold in (1, 2, 3, 5):
        states = [v for v in occ.values() if v >= threshold]
        print(
            f"  occupancy >= {threshold}: "
            f"states={len(states):4d} samples={sum(states):5d}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    start = time.perf_counter()
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 78")
    print("COMPRESSED MODULUS-NORMALIZED CHARACTER INVARIANTS")
    print("E1 CHARACTER PREDICTION WITHOUT EXACT STATE LOOKUP")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("LEAVE-ONE-MODULUS-OUT")
    print("WITHIN-TARGET PERMUTATION NULL")
    print("3-CLASS + NONZERO SIGN")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # Prime population
    # ----------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {time.perf_counter() - t0:.6f}s")

    # ----------------------------------------------------------------------
    # Targets
    # ----------------------------------------------------------------------

    targets = make_targets(primes, TARGETS_TOTAL, rng)

    print()
    print("2. TARGET SUMMARY")
    print("-" * 78)

    for t in targets[:24]:
        print(
            f"target {t.idx:4d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining targets generated but omitted")

    # ----------------------------------------------------------------------
    # Modulus validation
    # ----------------------------------------------------------------------

    all_moduli = list(dict.fromkeys(CYCLOTOMIC + CONTROLS))
    roots = prepare_root_table(all_moduli)

    print()
    print("3. MODULUS VALIDATION")
    print("-" * 78)

    failures = 0

    for ell in all_moduli:
        r1, r2 = roots[ell]
        ok = (
            (r1 * r1 + r1 + 1) % ell == 0
            and
            (r2 * r2 + r2 + 1) % ell == 0
            and
            r1 != r2
        )

        print(
            f"ell={ell:5d} roots=({r1:5d},{r2:5d}) "
            f"identity_ok={ok}"
        )

        if not ok:
            failures += 1

    print(f"validation failures = {failures}")

    if failures:
        raise RuntimeError("Modulus validation failed")

    print("status = PASS")

    # ----------------------------------------------------------------------
    # E1 identity validation
    # ----------------------------------------------------------------------

    print()
    print("4. E1 IDENTITY VALIDATION")
    print("-" * 78)

    identity_failures = 0

    for t in targets:
        e1_factor = (t.p * t.p - 1) * (t.q * t.q - 1) * t.s
        e1_s = t.s * ((t.n + 1) ** 2 - t.s ** 2)
        cubic = (
            t.s ** 3
            - (t.n + 1) ** 2 * t.s
            + e1_factor
        )

        ok = (
            e1_factor == e1_s
            and cubic == 0
        )

        if not ok:
            identity_failures += 1

    print(f"identity failures = {identity_failures}")

    if identity_failures:
        raise RuntimeError("E1 identity validation failed")

    print("status = PASS")

    # ----------------------------------------------------------------------
    # Target split
    # ----------------------------------------------------------------------

    train_targets, test_targets = split_targets(targets, rng)

    print()
    print("5. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_targets)}")
    print(f"test targets     = {len(test_targets)}")

    # ----------------------------------------------------------------------
    # Dataset construction
    # ----------------------------------------------------------------------

    t_dataset = time.perf_counter()

    cyclo_samples = build_dataset(
        targets,
        CYCLOTOMIC,
    )

    control_samples = build_dataset(
        targets,
        CONTROLS,
    )

    print()
    print("6. COMPRESSED DATASET")
    print("-" * 78)
    print(f"cyclotomic samples = {len(cyclo_samples)}")
    print(f"control samples    = {len(control_samples)}")
    print(
        f"construction time  = "
        f"{time.perf_counter() - t_dataset:.6f}s"
    )

    # ----------------------------------------------------------------------
    # Reports
    # ----------------------------------------------------------------------

    report_family(
        "CYCLOTOMIC",
        cyclo_samples,
        train_targets,
        test_targets,
        CYCLOTOMIC,
        rng,
    )

    report_family(
        "CONTROL",
        control_samples,
        train_targets,
        test_targets,
        CONTROLS,
        rng,
    )

    # ----------------------------------------------------------------------
    # Final comparison
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment is intentionally different from Experiments 76-77.

It does NOT use exact raw character-state lookup.

Instead it asks whether a SMALL CROSS-MODULUS SUMMARY of N-only
characters predicts the E1 character on unseen data.

A genuinely interesting result requires:

  1. compressed states to repeat frequently,
  2. target-holdout performance above majority baseline,
  3. cyclotomic performance to beat controls,
  4. leave-one-modulus-out performance to remain above chance,
  5. 2D target+modulus holdout to remain above chance,
  6. permutation p-values to remain small.

A particularly important outcome is:

    good cyclotomic OOS
    + good cyclotomic LOO
    + good cyclotomic 2D

while controls remain near chance.

That would be qualitatively stronger evidence than memorized
exact-state collisions.

A failure is also useful: it would strongly indicate that the
E1 correlations seen previously are not captured by any small
modulus-normalized character invariant of this form.
"""
    )

    print(f"total runtime = {time.perf_counter() - start:.6f}s")
    print("=" * 78)
    print("EXPERIMENT 78 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

