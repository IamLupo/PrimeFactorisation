#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 76
COMPRESSED CHARACTER-STATE E1 PREDICTION
NO RAW N MOD ELL^K LOOKUP
REUSABLE STATE / OCCUPANCY CONTROL
LEGENDRE + CUBIC CHARACTER SIGNATURES
TARGET + MODULUS + 2D HOLDOUT
WITHIN-TARGET PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

QUESTION
--------
Can E1 mod ell be predicted from a SMALL N-ONLY character signature
that has substantial state reuse across targets?

This experiment deliberately avoids raw n mod ell^k because Experiment 75
showed that raw-residue MI is dominated by state uniqueness.

For each modulus ell we build a compact signature from N-only characters:

  Legendre:
      chi(n)
      chi(n+1)
      chi(n-1)
      chi(4n+1)
      chi(n^2+n+1)
      chi(n^2+4n+1)

  Cubic character, only when ell == 1 mod 3:
      cubic(n)
      cubic(n+1)
      cubic(4n+1)

The signature is therefore a small finite-state object.

Prediction is ONLY performed when the exact signature has appeared in the
training set.  The prediction is the majority E1 residue class observed
for that signature.

We measure:

  1. state reuse
  2. coverage
  3. exact E1 prediction
  4. balanced exact prediction
  5. unseen-target performance
  6. leave-one-modulus-out
  7. two-dimensional target+modulus holdout
  8. within-target permutation null
  9. cyclotomic vs control comparison

IMPORTANT
---------
This does NOT claim that E1 is obtainable from N.
It tests whether a compact, reusable N-only character state contains
transferable information about E1 mod ell.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

SEED = 76001

NUM_TARGETS = 120
TRAIN_FRACTION = 0.75

# Small cyclotomic moduli are deliberately emphasized because state reuse
# is required for a meaningful lookup experiment.
CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

CONTROLS = [
    673, 4561, 4759, 6211, 7879, 7951, 8689, 9781
]

# For 2D holdout we use a smaller modulus set so that every held-out modulus
# still has enough repeated character states.
CYCLO_2D_HOLDOUT = [61, 127, 307]
CONTROL_2D_HOLDOUT = [6211, 7951]

PERMUTATIONS = 300


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int

    @property
    def e1(self) -> int:
        # E1 = (p^2 - 1)(q^2 - 1)(p+q)
        return (self.p * self.p - 1) * (self.q * self.q - 1) * self.s


@dataclass(frozen=True)
class Sample:
    target_idx: int
    ell: int
    signature: Tuple[int, ...]
    y: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(limit: int) -> List[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start::p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, v in enumerate(sieve) if v]


def build_prime_population() -> List[int]:
    # Same population scale used throughout the project.
    primes = sieve_primes(4_200_000)
    return primes


def generate_targets(
    primes: Sequence[int],
    count: int,
    rng: random.Random,
) -> List[Target]:
    eligible = [p for p in primes if 2_000_000 <= p <= 4_199_999]

    targets: List[Target] = []
    seen_n = set()

    while len(targets) < count:
        p, q = rng.sample(eligible, 2)

        if p == q:
            continue

        p, q = sorted((p, q))
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
# MODULAR CHARACTER HELPERS
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    """
    Return:
        +1 quadratic residue
        -1 non-residue
         0 divisible by p

    Requires p prime.
    """
    a %= p

    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1
    if value == p - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre result a={a} p={p} value={value}"
    )


def find_cubic_root_pair(ell: int) -> Tuple[int, int]:
    """
    Solve x^2+x+1 == 0 mod ell.

    Valid for primes ell == 1 mod 3.
    """
    if ell == 3 or ell % 3 != 1:
        raise ValueError(
            f"ell={ell} does not have two nontrivial cubic roots"
        )

    roots = [
        x for x in range(1, ell)
        if (x * x + x + 1) % ell == 0
    ]

    if len(roots) != 2:
        raise RuntimeError(
            f"Expected two cubic roots for ell={ell}, found {roots}"
        )

    return roots[0], roots[1]


def cubic_character(x: int, ell: int) -> int:
    """
    Return cubic residue class:

        0  if x == 0 mod ell
        1  if x is in cubic residue class
        2  for the other nonzero cubic class

    The implementation uses a primitive-root/discrete-log-free exponentiation
    test.  For ell == 1 mod 3 the multiplicative group has three cubic
    character classes.

    Let z be a primitive cube root of unity.
    We compute x^((ell-1)/3), which is in {1,z,z^2}.
    """
    x %= ell

    if x == 0:
        return 0

    _, z2 = find_cubic_root_pair(ell)
    z = (ell - 1 - z2) % ell

    v = pow(x, (ell - 1) // 3, ell)

    if v == 1:
        return 0

    if v == z:
        return 1

    if v == z2:
        return 2

    raise RuntimeError(
        f"Unexpected cubic character ell={ell} x={x} value={v}"
    )


# ============================================================================
# E1
# ============================================================================

def e1_identity_check(t: Target) -> bool:
    """
    Verify:

        E1 = (p^2-1)(q^2-1)(p+q)
           = s * ((n+1)^2 - s^2)

    and

        s^3 - (n+1)^2 s + E1 = 0.
    """
    e1_factor = (t.p * t.p - 1) * (t.q * t.q - 1) * t.s
    e1_sform = t.s * ((t.n + 1) ** 2 - t.s * t.s)
    cubic = (
        t.s ** 3
        - (t.n + 1) ** 2 * t.s
        + e1_factor
    )

    return (
        e1_factor == e1_sform
        and cubic == 0
        and e1_factor == t.e1
    )


# ============================================================================
# SIGNATURE CONSTRUCTION
# ============================================================================

def legendre_part(n: int, ell: int) -> Tuple[int, ...]:
    """
    Compact quadratic-character signature.

    Every expression is N-only.
    """
    values = [
        n,
        n + 1,
        n - 1,
        4 * n + 1,
        n * n + n + 1,
        n * n + 4 * n + 1,
    ]

    return tuple(
        legendre_symbol(v, ell)
        for v in values
    )


def signature_for(
    n: int,
    ell: int,
    cubic_cache: Dict[int, Tuple[int, int]],
) -> Tuple[int, ...]:
    sig = list(legendre_part(n, ell))

    if ell % 3 == 1:
        cubic_cache.setdefault(ell, find_cubic_root_pair(ell))

        sig.extend([
            cubic_character(n, ell),
            cubic_character(n + 1, ell),
            cubic_character(4 * n + 1, ell),
        ])

    return tuple(sig)


# ============================================================================
# MODEL
# ============================================================================

class SignatureLookup:
    """
    Exact reusable-state majority lookup.

    A state is predictive only when it has appeared in training.
    """

    def __init__(self) -> None:
        self.counts: Dict[Tuple[int, ...], Counter] = {}

    def fit(self, samples: Sequence[Sample]) -> None:
        tmp: Dict[Tuple[int, ...], Counter] = defaultdict(Counter)

        for sample in samples:
            tmp[sample.signature][sample.y] += 1

        self.counts = dict(tmp)

    def predict(
        self,
        sample: Sample,
    ) -> Optional[int]:
        counter = self.counts.get(sample.signature)

        if not counter:
            return None

        best_count = max(counter.values())
        winners = [
            y for y, c in counter.items()
            if c == best_count
        ]

        # Deterministic tie handling.
        return min(winners)


# ============================================================================
# METRICS
# ============================================================================

def accuracy(
    samples: Sequence[Sample],
    model: SignatureLookup,
) -> Tuple[float, float, int]:
    """
    Return:
        exact accuracy
        balanced accuracy
        covered samples
    """
    covered = [
        (s, model.predict(s))
        for s in samples
        if model.predict(s) is not None
    ]

    if not covered:
        return float("nan"), float("nan"), 0

    y_true = [s.y for s, _ in covered]
    y_pred = [p for _, p in covered]

    acc = sum(
        a == b for a, b in zip(y_true, y_pred)
    ) / len(y_true)

    classes = sorted(set(y_true))

    recalls = []
    for cls in classes:
        idx = [i for i, y in enumerate(y_true) if y == cls]
        if not idx:
            continue

        recalls.append(
            sum(y_pred[i] == cls for i in idx) / len(idx)
        )

    bal = statistics.fmean(recalls) if recalls else float("nan")

    return acc, bal, len(covered)


def majority_baseline(
    samples: Sequence[Sample],
) -> float:
    if not samples:
        return float("nan")

    counts = Counter(s.y for s in samples)
    return max(counts.values()) / len(samples)


def state_statistics(
    train: Sequence[Sample],
    test: Sequence[Sample],
) -> Tuple[int, float, int, float]:
    """
    Returns:
        number of training states
        mean training bucket size
        number of test samples whose state is known
        coverage
    """
    buckets = Counter(s.signature for s in train)

    nstates = len(buckets)

    mean_bucket = (
        statistics.fmean(buckets.values())
        if buckets else float("nan")
    )

    covered = sum(
        1 for s in test
        if s.signature in buckets
    )

    coverage = (
        covered / len(test)
        if test else float("nan")
    )

    return nstates, mean_bucket, covered, coverage


# ============================================================================
# PERMUTATION NULL
# ============================================================================

def permutation_accuracy(
    train: Sequence[Sample],
    test: Sequence[Sample],
    rng: random.Random,
    permutations: int,
) -> Tuple[float, float]:
    """
    Within-target label permutation.

    Labels are shuffled separately within each target so target composition
    is preserved while the relationship between signature and E1 is broken.
    """
    by_target: Dict[int, List[int]] = defaultdict(list)

    for s in train:
        by_target[s.target_idx].append(s.y)

    observed_model = SignatureLookup()
    observed_model.fit(train)

    observed, _, _ = accuracy(test, observed_model)

    if math.isnan(observed):
        return observed, float("nan")

    hits = 0

    for _ in range(permutations):
        permuted: List[Sample] = []

        for target_idx, labels in by_target.items():
            shuffled = list(labels)
            rng.shuffle(shuffled)

            target_samples = [
                s for s in train
                if s.target_idx == target_idx
            ]

            for s, y in zip(target_samples, shuffled):
                permuted.append(
                    Sample(
                        target_idx=s.target_idx,
                        ell=s.ell,
                        signature=s.signature,
                        y=y,
                    )
                )

        model = SignatureLookup()
        model.fit(permuted)

        score, _, _ = accuracy(test, model)

        if not math.isnan(score) and score >= observed:
            hits += 1

    p = (hits + 1) / (permutations + 1)

    return observed, p


# ============================================================================
# DATASET CONSTRUCTION
# ============================================================================

def build_samples(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> List[Sample]:
    cubic_cache: Dict[int, Tuple[int, int]] = {}
    samples: List[Sample] = []

    for t in targets:
        for ell in moduli:
            sig = signature_for(t.n, ell, cubic_cache)
            y = t.e1 % ell

            samples.append(
                Sample(
                    target_idx=t.idx,
                    ell=ell,
                    signature=sig,
                    y=y,
                )
            )

    return samples


# ============================================================================
# SPLITS
# ============================================================================

def split_target_ids(
    targets: Sequence[Target],
    rng: random.Random,
    fraction: float,
) -> Tuple[set, set]:
    ids = [t.idx for t in targets]
    rng.shuffle(ids)

    cut = int(round(len(ids) * fraction))

    return set(ids[:cut]), set(ids[cut:])


def subset_by_targets(
    samples: Sequence[Sample],
    ids: set,
) -> List[Sample]:
    return [
        s for s in samples
        if s.target_idx in ids
    ]


def subset_by_moduli(
    samples: Sequence[Sample],
    moduli: Sequence[int],
) -> List[Sample]:
    allowed = set(moduli)
    return [s for s in samples if s.ell in allowed]


# ============================================================================
# TARGET-HOLDOUT EXPERIMENT
# ============================================================================

def run_target_holdout(
    samples: Sequence[Sample],
    train_ids: set,
    test_ids: set,
    label: str,
    permutations: int,
    rng: random.Random,
) -> None:
    train = subset_by_targets(samples, train_ids)
    test = subset_by_targets(samples, test_ids)

    model = SignatureLookup()
    model.fit(train)

    acc, bal, covered = accuracy(test, model)

    nstates, mean_bucket, _, coverage = state_statistics(
        train,
        test,
    )

    baseline = majority_baseline(train)

    obs, perm_p = permutation_accuracy(
        train,
        test,
        rng,
        permutations,
    )

    print()
    print(label)
    print("-" * 78)
    print(f"train samples              = {len(train)}")
    print(f"test samples               = {len(test)}")
    print(f"training states            = {nstates}")
    print(f"mean training bucket      = {mean_bucket:.3f}")
    print(f"test coverage              = {coverage:.6f}")
    print(f"covered test samples       = {covered}")
    print(f"majority baseline          = {baseline:.6f}")
    print(f"OOS exact accuracy         = {acc:.6f}")
    print(f"OOS balanced accuracy     = {bal:.6f}")
    print(f"permutation p-value        = {perm_p:.6f}")


# ============================================================================
# LEAVE-ONE-MODULUS-OUT
# ============================================================================

def run_loo_modulus(
    samples: Sequence[Sample],
    moduli: Sequence[int],
) -> List[Tuple[int, float, float, float]]:
    results = []

    for held_out in moduli:
        train = [s for s in samples if s.ell != held_out]
        test = [s for s in samples if s.ell == held_out]

        # Since E1 mod ell is different for each ell, an exact lookup model
        # trained on other ell values should generally have zero useful
        # cross-ell labels. We include this test intentionally: if it fails,
        # the cross-modulus encoding is genuinely informative.
        #
        # We therefore convert the task into normalized character-state
        # prediction using a cross-modulus code below.

        model = SignatureLookup()
        model.fit(train)

        acc, bal, covered = accuracy(test, model)

        results.append(
            (
                held_out,
                acc,
                bal,
                covered / len(test) if test else float("nan"),
            )
        )

    return results


# ============================================================================
# CROSS-MODULUS NORMALIZED TEST
# ============================================================================

def normalized_label(y: int, ell: int) -> int:
    """
    Map the E1 residue to a Legendre signature of E1 itself.

    This creates a label that is comparable across different moduli:

        chi(E1)

    The exact residue prediction is still tested separately.
    """
    return legendre_symbol(y, ell)


def build_normalized_samples(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> List[Sample]:
    """
    Same N-only signature, but label is the quadratic character of E1.

    This permits genuine leave-one-modulus-out generalization.
    """
    cubic_cache: Dict[int, Tuple[int, int]] = {}
    out: List[Sample] = []

    for t in targets:
        for ell in moduli:
            sig = signature_for(t.n, ell, cubic_cache)
            y = normalized_label(t.e1, ell)

            out.append(
                Sample(
                    target_idx=t.idx,
                    ell=ell,
                    signature=sig,
                    y=y,
                )
            )

    return out


def run_cross_modulus_normalized(
    samples: Sequence[Sample],
    cyclotomic: Sequence[int],
    controls: Sequence[int],
    train_ids: set,
    test_ids: set,
    rng: random.Random,
) -> None:
    """
    Genuine cross-modulus prediction because the target label is comparable
    across ell.
    """
    print()
    print("NORMALIZED E1 CHARACTER TASK")
    print("-" * 78)

    for family_name, family in (
        ("CYCLOTOMIC", cyclotomic),
        ("CONTROL", controls),
    ):
        fam_train = [
            s for s in samples
            if s.ell in family and s.target_idx in train_ids
        ]

        fam_test = [
            s for s in samples
            if s.ell in family and s.target_idx in test_ids
        ]

        model = SignatureLookup()
        model.fit(fam_train)

        acc, bal, covered = accuracy(fam_test, model)

        nstates, mean_bucket, _, coverage = state_statistics(
            fam_train,
            fam_test,
        )

        obs, p = permutation_accuracy(
            fam_train,
            fam_test,
            rng,
            200,
        )

        print()
        print(f"{family_name}")
        print(f"  train samples          = {len(fam_train)}")
        print(f"  test samples           = {len(fam_test)}")
        print(f"  states                 = {nstates}")
        print(f"  mean bucket            = {mean_bucket:.3f}")
        print(f"  coverage               = {coverage:.6f}")
        print(f"  exact character acc    = {acc:.6f}")
        print(f"  balanced accuracy      = {bal:.6f}")
        print(f"  permutation p          = {p:.6f}")


# ============================================================================
# MATCHED STATE-OCCUPANCY TEST
# ============================================================================

def matched_state_occupancy_test(
    samples: Sequence[Sample],
    cyclotomic: Sequence[int],
    controls: Sequence[int],
    train_ids: set,
    test_ids: set,
) -> None:
    """
    Compare only signatures with training bucket size >= K.

    This removes the strongest uniqueness artifact from Experiment 75.
    """
    print()
    print("MATCHED STATE-OCCUPANCY ANALYSIS")
    print("-" * 78)

    for min_bucket in (2, 3, 5):
        print(f"\nminimum training bucket = {min_bucket}")

        for family_name, family in (
            ("CYCLOTOMIC", cyclotomic),
            ("CONTROL", controls),
        ):
            train = [
                s for s in samples
                if s.ell in family and s.target_idx in train_ids
            ]

            test = [
                s for s in samples
                if s.ell in family and s.target_idx in test_ids
            ]

            buckets = Counter(s.signature for s in train)

            eligible_states = {
                sig for sig, count in buckets.items()
                if count >= min_bucket
            }

            filtered_test = [
                s for s in test
                if s.signature in eligible_states
            ]

            model = SignatureLookup()
            model.fit(train)

            acc, bal, covered = accuracy(filtered_test, model)

            coverage = (
                len(filtered_test) / len(test)
                if test else float("nan")
            )

            print(
                f"  {family_name:10s} "
                f"states={len(eligible_states):4d} "
                f"coverage={coverage:.6f} "
                f"acc={acc:.6f} "
                f"bal={bal:.6f}"
            )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    started = time.perf_counter()
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 76")
    print("COMPRESSED CHARACTER-STATE E1 PREDICTION")
    print("NO RAW N MOD ELL^K LOOKUP")
    print("REUSABLE STATE / OCCUPANCY CONTROL")
    print("LEGENDRE + CUBIC CHARACTER SIGNATURES")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("WITHIN-TARGET PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ---------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ---------------------------------------------------------------------
    t0 = time.perf_counter()
    primes = build_prime_population()
    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(
        f"generation time  = {time.perf_counter() - t0:.6f}s"
    )

    # ---------------------------------------------------------------------
    # 2. TARGETS
    # ---------------------------------------------------------------------
    targets = generate_targets(
        primes,
        NUM_TARGETS,
        rng,
    )

    print()
    print("2. TARGET SUMMARY")
    print("-" * 78)
    print(f"total targets = {len(targets)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining targets omitted")

    # ---------------------------------------------------------------------
    # 3. MODULI
    # ---------------------------------------------------------------------
    print()
    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic = {CYCLOTOMIC}")
    print(f"controls   = {CONTROLS}")

    # ---------------------------------------------------------------------
    # 4. IDENTITY VALIDATION
    # ---------------------------------------------------------------------
    failures = [
        t.idx for t in targets
        if not e1_identity_check(t)
    ]

    print()
    print("4. E1 IDENTITY VALIDATION")
    print("-" * 78)
    print(f"identity failures = {len(failures)}")

    if failures:
        print(f"failure targets = {failures}")
        raise RuntimeError("E1 identity validation failed")

    print("status = PASS")

    # ---------------------------------------------------------------------
    # 5. TARGET SPLIT
    # ---------------------------------------------------------------------
    train_ids, test_ids = split_target_ids(
        targets,
        rng,
        TRAIN_FRACTION,
    )

    print()
    print("5. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_ids)}")
    print(f"test targets     = {len(test_ids)}")

    # ---------------------------------------------------------------------
    # 6. BUILD DATASETS
    # ---------------------------------------------------------------------
    build_started = time.perf_counter()

    cyclo_samples = build_samples(
        targets,
        CYCLOTOMIC,
    )

    control_samples = build_samples(
        targets,
        CONTROLS,
    )

    normalized_samples = build_normalized_samples(
        targets,
        CYCLOTOMIC + CONTROLS,
    )

    print()
    print("6. DATASET CONSTRUCTION")
    print("-" * 78)
    print(f"cyclotomic samples = {len(cyclo_samples)}")
    print(f"control samples    = {len(control_samples)}")
    print(
        f"construction time  = "
        f"{time.perf_counter() - build_started:.6f}s"
    )

    # ---------------------------------------------------------------------
    # 7. MODULUS-BY-MODULUS TARGET HOLDOUT
    # ---------------------------------------------------------------------
    print()
    print("7. EXACT E1 RESIDUE TARGET HOLDOUT")
    print("-" * 78)

    for family_name, samples in (
        ("CYCLOTOMIC", cyclo_samples),
        ("CONTROL", control_samples),
    ):
        run_target_holdout(
            samples,
            train_ids,
            test_ids,
            family_name,
            PERMUTATIONS,
            rng,
        )

    # ---------------------------------------------------------------------
    # 8. STATE OCCUPANCY
    # ---------------------------------------------------------------------
    matched_state_occupancy_test(
        cyclo_samples + control_samples,
        CYCLOTOMIC,
        CONTROLS,
        train_ids,
        test_ids,
    )

    # ---------------------------------------------------------------------
    # 9. GENUINE CROSS-MODULUS NORMALIZED TEST
    # ---------------------------------------------------------------------
    run_cross_modulus_normalized(
        normalized_samples,
        CYCLOTOMIC,
        CONTROLS,
        train_ids,
        test_ids,
        rng,
    )

    # ---------------------------------------------------------------------
    # 10. LEAVE-ONE-MODULUS-OUT ON NORMALIZED LABEL
    # ---------------------------------------------------------------------
    print()
    print("10. LEAVE-ONE-MODULUS-OUT NORMALIZED CHARACTER")
    print("-" * 78)

    for family_name, family in (
        ("CYCLOTOMIC", CYCLOTOMIC),
        ("CONTROL", CONTROLS),
    ):
        print(f"\n{family_name}")

        loo_scores = []

        for held_out in family:
            train = [
                s for s in normalized_samples
                if s.ell in family and s.ell != held_out
            ]

            test = [
                s for s in normalized_samples
                if s.ell == held_out
            ]

            model = SignatureLookup()
            model.fit(train)

            acc, bal, covered = accuracy(test, model)

            loo_scores.append(acc)

            print(
                f"  ell={held_out:5d} "
                f"accuracy={acc:.6f} "
                f"balanced={bal:.6f} "
                f"coverage="
                f"{covered / len(test):.6f}"
            )

        valid = [
            x for x in loo_scores
            if not math.isnan(x)
        ]

        print(
            f"  mean LOO = "
            f"{statistics.fmean(valid) if valid else float('nan'):.6f}"
        )

    # ---------------------------------------------------------------------
    # 11. TWO-DIMENSIONAL HOLDOUT
    # ---------------------------------------------------------------------
    print()
    print("11. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)

    for family_name, family, heldout in (
        (
            "CYCLOTOMIC",
            CYCLOTOMIC,
            set(CYCLO_2D_HOLDOUT),
        ),
        (
            "CONTROL",
            CONTROLS,
            set(CONTROL_2D_HOLDOUT),
        ),
    ):
        train = [
            s for s in normalized_samples
            if s.ell in family
            and s.ell not in heldout
            and s.target_idx in train_ids
        ]

        test = [
            s for s in normalized_samples
            if s.ell in heldout
            and s.target_idx in test_ids
        ]

        model = SignatureLookup()
        model.fit(train)

        acc, bal, covered = accuracy(test, model)

        coverage = covered / len(test) if test else float("nan")

        print(f"{family_name}")
        print(f"  held out moduli = {sorted(heldout)}")
        print(f"  accuracy        = {acc:.6f}")
        print(f"  balanced        = {bal:.6f}")
        print(f"  coverage        = {coverage:.6f}")

    # ---------------------------------------------------------------------
    # 12. SIMPLE CHARACTER BASELINES
    # ---------------------------------------------------------------------
    print()
    print("12. SIMPLE N-ONLY CHARACTER BASELINES")
    print("-" * 78)

    for name, expression in (
        ("chi(n)", lambda n, ell: legendre_symbol(n, ell)),
        ("chi(n+1)", lambda n, ell: legendre_symbol(n + 1, ell)),
        ("chi(n-1)", lambda n, ell: legendre_symbol(n - 1, ell)),
        ("chi(4n+1)", lambda n, ell: legendre_symbol(4 * n + 1, ell)),
        (
            "chi(n^2+n+1)",
            lambda n, ell: legendre_symbol(
                n * n + n + 1,
                ell,
            ),
        ),
        (
            "chi(n^2+4n+1)",
            lambda n, ell: legendre_symbol(
                n * n + 4 * n + 1,
                ell,
            ),
        ),
    ):
        scores = []

        for t in targets:
            for ell in CYCLOTOMIC:
                pred = expression(t.n, ell)
                actual = normalized_label(t.e1, ell)

                if t.idx in test_ids:
                    scores.append(pred == actual)

        score = (
            sum(scores) / len(scores)
            if scores else float("nan")
        )

        print(f"{name:20s} test accuracy={score:.6f}")

    # ---------------------------------------------------------------------
    # 13. FINAL DIAGNOSTIC
    # ---------------------------------------------------------------------
    print()
    print("13. FINAL DIAGNOSTIC")
    print("-" * 78)
    print(
        "The decisive test is NOT raw mutual information.\n"
        "It is whether a SMALL REUSABLE N-only character state predicts\n"
        "a comparable E1-derived character on unseen targets and unseen\n"
        "moduli, after controlling for state occupancy.\n"
    )

    print("A positive result requires:")
    print("  1. meaningful state reuse")
    print("  2. OOS accuracy above the matched baseline")
    print("  3. improvement for cyclotomic over controls")
    print("  4. survival under leave-one-modulus-out")
    print("  5. survival under 2D target+modulus holdout")
    print("  6. permutation p-values that remain small")

    print()
    print(
        "A failure would be informative: it would strongly suggest that "
        "the\n"
        "E1 residue correlation seen previously is primarily a property "
        "of raw\n"
        "residue collisions rather than a reusable cyclotomic law."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 76 COMPLETE")
    print(
        f"total runtime = "
        f"{time.perf_counter() - started:.6f}s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

