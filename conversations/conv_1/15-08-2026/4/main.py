"""
==============================================================================
KAPPA EXPERIMENT 62
CROSS-MODULUS / CROSS-TARGET ORIENTATION GENERALIZATION
UNSEEN ELL TEST + UNSEEN TARGET TEST + PERMUTATION NULL
NO CSV OUTPUT
==============================================================================

Purpose
-------
Experiment 61R showed a preliminary N-only orientation signal:

    cyclotomic holdout  ~= 0.627
    control holdout     ~= 0.472

But the test still had only 51 cyclotomic observations, and the
orientation model had not been tested rigorously on completely unseen
cyclotomic moduli.

This experiment asks three stricter questions:

A. TARGET HOLDOUT
   Train on some targets and test on completely unseen targets.

B. MODULUS HOLDOUT
   Train on 12 cyclotomic moduli and test on the omitted modulus.
   Repeat for every cyclotomic ell.

C. TWO-DIMENSIONAL HOLDOUT
   Test simultaneously on unseen targets and an unseen ell.

Controls are subjected to the same modulus-holdout procedure.

A permutation null estimates how often the observed target-holdout
accuracy could arise after destroying the orientation labels.

Important
---------
This is an information test only.

No factorization is performed.
No candidate sieve is claimed.
No CSV files are written.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

try:
    import numpy as np
    from sklearn.linear_model import LogisticRegression
except ImportError as exc:
    raise RuntimeError(
        "This experiment requires numpy and scikit-learn.\n"
        "Install with:\n"
        "  pip install numpy scikit-learn"
    ) from exc


# ============================================================================
# CONSTANTS
# ============================================================================

SEED = 62062

PRIME_LIMIT = 4_200_000

NUM_TARGETS = 120
LEGACY_TARGETS = 24
EXTRA_TARGETS = NUM_TARGETS - LEGACY_TARGETS

TRAIN_TARGET_FRACTION = 0.75

PERMUTATIONS = 500

MIN_S = 4_000_000
MAX_S = 8_399_998

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723,
]

RAW_CONTROLS = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781,
]

# Only primes ell == 1 mod 3 have two nontrivial roots of
#
#     x^2 + x + 1 == 0 mod ell.
#
# This is the orientation-valid control subset.
CONTROLS = [p for p in RAW_CONTROLS if p % 3 == 1]


# ============================================================================
# LEGACY TARGETS
# ============================================================================

LEGACY_PAIRS = [
    (3318013, 4042603),
    (2129167, 3402323),
    (2224517, 3978749),
    (3685051, 4020281),
    (2399627, 2452649),
    (2593039, 2996527),
    (2149859, 2772097),
    (2060543, 2514401),
    (2675423, 2883973),
    (2828887, 3960137),
    (3497381, 3793241),
    (2193509, 4011353),
    (3429689, 3983927),
    (2515871, 3050581),
    (2261297, 3374827),
    (2345537, 3673349),
    (2212039, 3452809),
    (2175373, 2476921),
    (2438509, 4188577),
    (2367553, 2399407),
    (2072897, 2087077),
    (3552023, 3707453),
    (3112909, 3210167),
    (2965819, 3239963),
]


@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


@dataclass(frozen=True)
class LocalSample:
    target_idx: int
    ell: int
    label: int
    features: Tuple[float, ...]


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(limit: int) -> List[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, is_prime in enumerate(sieve) if is_prime]


# ============================================================================
# BASIC MODULAR FUNCTIONS
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    v = pow(a, (p - 1) // 2, p)

    if v == 1:
        return 1

    if v == p - 1:
        return -1

    raise AssertionError(
        f"Unexpected Legendre value {v} for a={a}, p={p}"
    )


def roots_cyclotomic(ell: int) -> Tuple[int, int]:
    """
    Return the two nontrivial roots of x^2+x+1 modulo ell.

    Roots exist exactly for primes ell == 1 mod 3.
    """
    if ell % 3 != 1:
        raise ValueError(
            f"{ell} is not orientation-valid: ell % 3 != 1"
        )

    roots = [
        x for x in range(1, ell)
        if (x * x + x + 1) % ell == 0
    ]

    if len(roots) != 2:
        raise RuntimeError(
            f"{ell} expected two roots, found {roots}"
        )

    return tuple(sorted(roots))


def cubic_root_of_unity(ell: int) -> int:
    """
    Find a primitive cube root of unity modulo ell.
    """
    r1, r2 = roots_cyclotomic(ell)

    # Deterministic choice.
    return r1


def cubic_class(a: int, ell: int) -> int:
    """
    Cubic character class for a nonzero residue.

    Returns:
        0 for residue class 1 under the cubic character,
        1 for zeta,
        2 for zeta^2.

    Zero is encoded as -1 because it is degenerate for orientation.
    """
    a %= ell

    if a == 0:
        return -1

    zeta = cubic_root_of_unity(ell)
    zeta2 = (zeta * zeta) % ell

    value = pow(a, (ell - 1) // 3, ell)

    if value == 1:
        return 0
    if value == zeta:
        return 1
    if value == zeta2:
        return 2

    raise AssertionError(
        f"Unexpected cubic-character value {value} for ell={ell}"
    )


# ============================================================================
# TARGET GENERATION
# ============================================================================

def build_targets(primes: List[int]) -> List[Target]:
    rng = random.Random(SEED)

    targets: List[Target] = []

    # Preserve the previous 24 targets exactly.
    for idx, (p, q) in enumerate(LEGACY_PAIRS, start=1):
        n = p * q
        s = p + q

        if p == q:
            raise RuntimeError("Legacy target has p == q.")

        if not (MIN_S <= s <= MAX_S):
            raise RuntimeError(
                f"Legacy target {idx} has sum outside the scan domain."
            )

        targets.append(Target(idx, p, q, n, s))

    prime_set = set(primes)

    # Generate additional independent targets.
    seen_pairs = {
        tuple(sorted((t.p, t.q)))
        for t in targets
    }

    while len(targets) < NUM_TARGETS:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        pair = tuple(sorted((p, q)))

        if pair in seen_pairs:
            continue

        s = p + q

        # Keep the same factor-sum domain as previous experiments.
        if not (MIN_S <= s <= MAX_S):
            continue

        # Both primes should remain in the intended population.
        if p not in prime_set or q not in prime_set:
            continue

        n = p * q

        idx = len(targets) + 1
        targets.append(Target(idx, p, q, n, s))
        seen_pairs.add(pair)

    return targets


# ============================================================================
# N-ONLY FEATURES
# ============================================================================

def build_features(n: int, ell: int) -> Tuple[float, ...]:
    """
    Feature vector intentionally depends ONLY on n and ell.

    No p, q, s, discriminant, or true local root orientation is exposed.
    """

    chi_n = legendre_symbol(n, ell)
    chi_np1 = legendre_symbol(n + 1, ell)
    chi_nm1 = legendre_symbol(n - 1, ell)
    chi_4np1 = legendre_symbol(4 * n + 1, ell)

    # N-only conjugate-product quantity from Experiment 56.
    G = 16 * n * n + 4 * n + 1
    chi_G = legendre_symbol(G, ell)

    c_n = cubic_class(n, ell)
    c_np1 = cubic_class(n + 1, ell)
    c_4np1 = cubic_class(4 * n + 1, ell)

    # One-hot cubic features.
    return (
        float(chi_n),
        float(chi_np1),
        float(chi_nm1),
        float(chi_4np1),
        float(chi_G),

        float(c_n == 0),
        float(c_n == 1),
        float(c_n == 2),

        float(c_np1 == 0),
        float(c_np1 == 1),
        float(c_np1 == 2),

        float(c_4np1 == 0),
        float(c_4np1 == 1),
        float(c_4np1 == 2),
    )


# ============================================================================
# LOCAL ORIENTATION
# ============================================================================

def local_orientation(n: int, ell: int):
    """
    Determine the ordered local orientation.

    Let r1 < r2 be the two roots of x^2+x+1.

    Define
        d_i = r_i^2 - 4n mod ell.

    The cell is CLEAN exactly when:
        chi(d1), chi(d2) = (+1,-1) or (-1,+1).

    Label:
        1  => first root is QR, second is non-QR
        0  => first root is non-QR, second is QR

    Degenerate cells are returned as None.
    """

    r1, r2 = roots_cyclotomic(ell)

    d1 = (r1 * r1 - 4 * n) % ell
    d2 = (r2 * r2 - 4 * n) % ell

    c1 = legendre_symbol(d1, ell)
    c2 = legendre_symbol(d2, ell)

    if (c1, c2) == (1, -1):
        return 1

    if (c1, c2) == (-1, 1):
        return 0

    return None


# ============================================================================
# DATASET BUILDING
# ============================================================================

def build_dataset(
    targets: List[Target],
    moduli: List[int],
) -> Tuple[List[LocalSample], int]:
    samples: List[LocalSample] = []
    degenerate = 0

    for t in targets:
        for ell in moduli:
            label = local_orientation(t.n, ell)

            if label is None:
                degenerate += 1
                continue

            features = build_features(t.n, ell)

            samples.append(
                LocalSample(
                    target_idx=t.idx,
                    ell=ell,
                    label=label,
                    features=features,
                )
            )

    return samples, degenerate


# ============================================================================
# MODEL
# ============================================================================

def make_model() -> LogisticRegression:
    return LogisticRegression(
        solver="liblinear",
        penalty="l2",
        C=1.0,
        max_iter=2000,
        random_state=SEED,
    )


def fit_predict(
    train_samples: List[LocalSample],
    test_samples: List[LocalSample],
):
    if not train_samples or not test_samples:
        return []

    X_train = np.array(
        [s.features for s in train_samples],
        dtype=float,
    )
    y_train = np.array(
        [s.label for s in train_samples],
        dtype=int,
    )

    X_test = np.array(
        [s.features for s in test_samples],
        dtype=float,
    )

    # Guard against pathological training partitions.
    if len(set(y_train.tolist())) < 2:
        majority = int(round(float(np.mean(y_train))))
        return [majority] * len(test_samples)

    model = make_model()
    model.fit(X_train, y_train)

    return model.predict(X_test).astype(int).tolist()


def accuracy(
    predictions: List[int],
    samples: List[LocalSample],
) -> float:
    if not samples:
        return float("nan")

    correct = sum(
        pred == sample.label
        for pred, sample in zip(predictions, samples)
    )

    return correct / len(samples)


# ============================================================================
# TARGET HOLDOUT
# ============================================================================

def target_holdout(
    samples: List[LocalSample],
    train_targets: set[int],
    test_targets: set[int],
) -> float:
    train = [
        s for s in samples
        if s.target_idx in train_targets
    ]

    test = [
        s for s in samples
        if s.target_idx in test_targets
    ]

    pred = fit_predict(train, test)
    return accuracy(pred, test)


# ============================================================================
# MODULUS HOLDOUT
# ============================================================================

def leave_one_modulus_out(
    samples: List[LocalSample],
    moduli: List[int],
) -> Dict[int, float]:

    results: Dict[int, float] = {}

    for ell in moduli:
        train = [
            s for s in samples
            if s.ell != ell
        ]

        test = [
            s for s in samples
            if s.ell == ell
        ]

        pred = fit_predict(train, test)
        results[ell] = accuracy(pred, test)

    return results


# ============================================================================
# TWO-DIMENSIONAL HOLDOUT
# ============================================================================

def two_dimensional_holdout(
    samples: List[LocalSample],
    test_targets: set[int],
    heldout_ell: int,
) -> float:

    train = [
        s for s in samples
        if s.target_idx not in test_targets
        and s.ell != heldout_ell
    ]

    test = [
        s for s in samples
        if s.target_idx in test_targets
        and s.ell == heldout_ell
    ]

    pred = fit_predict(train, test)
    return accuracy(pred, test)


# ============================================================================
# CONTROL TESTS
# ============================================================================

def evaluate_controls(
    targets: List[Target],
    controls: List[int],
):
    control_samples, degenerate = build_dataset(
        targets,
        controls,
    )

    loo = leave_one_modulus_out(
        control_samples,
        controls,
    )

    return control_samples, degenerate, loo


# ============================================================================
# FIXED RULE BASELINES
# ============================================================================

def baseline_accuracy(
    samples: List[LocalSample],
    feature_index: int,
) -> float:
    """
    Majority-like sign classifier:

        feature >= 0 => predicts 1
        feature < 0  => predicts 0

    Used only as a simple transparent baseline.
    """

    if not samples:
        return float("nan")

    correct = 0

    for s in samples:
        pred = 1 if s.features[feature_index] >= 0 else 0

        if pred == s.label:
            correct += 1

    return correct / len(samples)


# ============================================================================
# PERMUTATION NULL
# ============================================================================

def permutation_target_holdout(
    samples: List[LocalSample],
    train_targets: set[int],
    test_targets: set[int],
    observed_accuracy: float,
    rng: random.Random,
) -> float:

    train = [
        s for s in samples
        if s.target_idx in train_targets
    ]

    test = [
        s for s in samples
        if s.target_idx in test_targets
    ]

    if not train or not test:
        return float("nan")

    labels = [s.label for s in train]
    ge_count = 0

    for _ in range(PERMUTATIONS):
        shuffled = labels[:]
        rng.shuffle(shuffled)

        perm_train = [
            LocalSample(
                target_idx=s.target_idx,
                ell=s.ell,
                label=shuffled[i],
                features=s.features,
            )
            for i, s in enumerate(train)
        ]

        pred = fit_predict(perm_train, test)
        acc = accuracy(pred, test)

        if acc >= observed_accuracy:
            ge_count += 1

    return (ge_count + 1) / (PERMUTATIONS + 1)


# ============================================================================
# TARGET-LEVEL REPEATED SPLITS
# ============================================================================

def repeated_target_holdout(
    samples: List[LocalSample],
    target_indices: List[int],
    repeats: int,
) -> List[float]:

    rng = random.Random(SEED + 1000)

    scores = []

    target_indices = list(target_indices)

    for _ in range(repeats):
        shuffled = target_indices[:]
        rng.shuffle(shuffled)

        split = int(len(shuffled) * TRAIN_TARGET_FRACTION)

        train_targets = set(shuffled[:split])
        test_targets = set(shuffled[split:])

        score = target_holdout(
            samples,
            train_targets,
            test_targets,
        )

        if not math.isnan(score):
            scores.append(score)

    return scores


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 62")
    print("CROSS-MODULUS / CROSS-TARGET ORIENTATION GENERALIZATION")
    print("UNSEEN ELL + UNSEEN TARGET + PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # ---------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ---------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LIMIT)
    prime_time = time.perf_counter() - t0

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.6f}s")

    # ---------------------------------------------------------------------
    # 2. TARGETS
    # ---------------------------------------------------------------------

    targets = build_targets(primes)

    print()
    print("2. TARGETS")
    print("-" * 78)
    print(f"total targets = {len(targets)}")
    print(f"legacy targets = {LEGACY_TARGETS}")
    print(f"new targets    = {EXTRA_TARGETS}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    print("...")
    print("remaining generated targets omitted from detailed listing")

    # ---------------------------------------------------------------------
    # 3. MODULI
    # ---------------------------------------------------------------------

    print()
    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic = {CYCLOTOMIC}")
    print(f"raw controls = {RAW_CONTROLS}")
    print(f"valid controls = {CONTROLS}")
    print(
        "excluded controls = "
        + str([x for x in RAW_CONTROLS if x not in CONTROLS])
    )

    # ---------------------------------------------------------------------
    # 4. LOCAL VALIDATION
    # ---------------------------------------------------------------------

    print()
    print("4. ROOT VALIDATION")
    print("-" * 78)

    for ell in CYCLOTOMIC + CONTROLS:
        r1, r2 = roots_cyclotomic(ell)

        ok = (
            (r1 * r1 + r1 + 1) % ell == 0
            and (r2 * r2 + r2 + 1) % ell == 0
            and r1 != r2
        )

        print(
            f"ell={ell:5d} roots=({r1:5d},{r2:5d}) "
            f"identity_ok={ok}"
        )

    # ---------------------------------------------------------------------
    # 5. BUILD DATASETS
    # ---------------------------------------------------------------------

    t0 = time.perf_counter()

    cyclo_samples, cyclo_deg = build_dataset(
        targets,
        CYCLOTOMIC,
    )

    control_samples, control_deg = build_dataset(
        targets,
        CONTROLS,
    )

    dataset_time = time.perf_counter() - t0

    print()
    print("5. DATASET CONSTRUCTION")
    print("-" * 78)
    print(f"cyclotomic samples = {len(cyclo_samples):,}")
    print(f"cyclotomic degenerate = {cyclo_deg:,}")
    print(f"control samples     = {len(control_samples):,}")
    print(f"control degenerate  = {control_deg:,}")
    print(f"time = {dataset_time:.6f}s")

    # ---------------------------------------------------------------------
    # 6. TRAIN / TEST TARGET SPLIT
    # ---------------------------------------------------------------------

    all_target_ids = [t.idx for t in targets]

    rng = random.Random(SEED)

    shuffled = all_target_ids[:]
    rng.shuffle(shuffled)

    split = int(len(shuffled) * TRAIN_TARGET_FRACTION)

    train_targets = set(shuffled[:split])
    test_targets = set(shuffled[split:])

    print()
    print("6. TARGET-LEVEL HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_targets)}")
    print(f"test targets     = {len(test_targets)}")
    print(f"train ids = {sorted(train_targets)}")
    print(f"test ids  = {sorted(test_targets)}")

    # ---------------------------------------------------------------------
    # 7. TARGET HOLDOUT
    # ---------------------------------------------------------------------

    cyclo_target_acc = target_holdout(
        cyclo_samples,
        train_targets,
        test_targets,
    )

    control_target_acc = target_holdout(
        control_samples,
        train_targets,
        test_targets,
    )

    print()
    print("7. TARGET HOLDOUT ACCURACY")
    print("-" * 78)
    print(f"cyclotomic accuracy = {cyclo_target_acc:.6f}")
    print(f"control accuracy    = {control_target_acc:.6f}")
    print(
        f"difference          = "
        f"{cyclo_target_acc - control_target_acc:+.6f}"
    )

    # ---------------------------------------------------------------------
    # 8. REPEATED TARGET HOLDOUT
    # ---------------------------------------------------------------------

    cyclo_repeated = repeated_target_holdout(
        cyclo_samples,
        all_target_ids,
        repeats=25,
    )

    control_repeated = repeated_target_holdout(
        control_samples,
        all_target_ids,
        repeats=25,
    )

    print()
    print("8. REPEATED TARGET HOLDOUT")
    print("-" * 78)

    print(
        f"cyclotomic mean   = {statistics.mean(cyclo_repeated):.6f}"
    )
    print(
        f"cyclotomic median = {statistics.median(cyclo_repeated):.6f}"
    )
    print(
        f"cyclotomic min    = {min(cyclo_repeated):.6f}"
    )
    print(
        f"cyclotomic max    = {max(cyclo_repeated):.6f}"
    )

    print(
        f"control mean      = {statistics.mean(control_repeated):.6f}"
    )
    print(
        f"control median    = {statistics.median(control_repeated):.6f}"
    )
    print(
        f"control min       = {min(control_repeated):.6f}"
    )
    print(
        f"control max       = {max(control_repeated):.6f}"
    )

    # ---------------------------------------------------------------------
    # 9. LEAVE-ONE-MODULUS-OUT
    # ---------------------------------------------------------------------

    cyclo_loo = leave_one_modulus_out(
        cyclo_samples,
        CYCLOTOMIC,
    )

    control_loo = leave_one_modulus_out(
        control_samples,
        CONTROLS,
    )

    print()
    print("9. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    print("CYCLOTOMIC")
    for ell in CYCLOTOMIC:
        print(
            f"  ell={ell:5d} accuracy={cyclo_loo[ell]:.6f}"
        )

    print("CONTROL")
    for ell in CONTROLS:
        print(
            f"  ell={ell:5d} accuracy={control_loo[ell]:.6f}"
        )

    print()
    print(
        f"cyclotomic LOO mean   = "
        f"{statistics.mean(cyclo_loo.values()):.6f}"
    )
    print(
        f"cyclotomic LOO median = "
        f"{statistics.median(cyclo_loo.values()):.6f}"
    )
    print(
        f"control LOO mean      = "
        f"{statistics.mean(control_loo.values()):.6f}"
    )
    print(
        f"control LOO median    = "
        f"{statistics.median(control_loo.values()):.6f}"
    )

    # ---------------------------------------------------------------------
    # 10. TWO-DIMENSIONAL HOLDOUT
    # ---------------------------------------------------------------------

    # Hold out the same target test set and each ell in turn.
    cyclo_2d = {}

    for ell in CYCLOTOMIC:
        cyclo_2d[ell] = two_dimensional_holdout(
            cyclo_samples,
            test_targets,
            ell,
        )

    control_2d = {}

    for ell in CONTROLS:
        control_2d[ell] = two_dimensional_holdout(
            control_samples,
            test_targets,
            ell,
        )

    print()
    print("10. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)
    print(
        "Both the target group and the modulus are unseen during training."
    )

    print("CYCLOTOMIC")
    for ell in CYCLOTOMIC:
        print(
            f"  ell={ell:5d} accuracy={cyclo_2d[ell]:.6f}"
        )

    print("CONTROL")
    for ell in CONTROLS:
        print(
            f"  ell={ell:5d} accuracy={control_2d[ell]:.6f}"
        )

    print()
    print(
        f"cyclotomic 2D mean   = "
        f"{statistics.mean(cyclo_2d.values()):.6f}"
    )
    print(
        f"cyclotomic 2D median = "
        f"{statistics.median(cyclo_2d.values()):.6f}"
    )
    print(
        f"control 2D mean      = "
        f"{statistics.mean(control_2d.values()):.6f}"
    )
    print(
        f"control 2D median    = "
        f"{statistics.median(control_2d.values()):.6f}"
    )

    # ---------------------------------------------------------------------
    # 11. SIMPLE N-ONLY BASELINES
    # ---------------------------------------------------------------------

    baseline_names = {
        0: "chi(n)",
        1: "chi(n+1)",
        2: "chi(n-1)",
        3: "chi(4n+1)",
        4: "chi(G)",
    }

    print()
    print("11. SIMPLE N-ONLY BASELINES")
    print("-" * 78)

    for idx, name in baseline_names.items():
        score = baseline_accuracy(
            [
                s for s in cyclo_samples
                if s.target_idx in test_targets
            ],
            idx,
        )

        print(
            f"{name:12s} test_accuracy={score:.6f}"
        )

    # ---------------------------------------------------------------------
    # 12. PERMUTATION NULL
    # ---------------------------------------------------------------------

    perm_rng = random.Random(SEED + 999)

    observed = cyclo_target_acc

    p_value = permutation_target_holdout(
        cyclo_samples,
        train_targets,
        test_targets,
        observed,
        perm_rng,
    )

    print()
    print("12. PERMUTATION NULL")
    print("-" * 78)
    print(f"observed cyclotomic target accuracy = {observed:.6f}")
    print(f"permutations                        = {PERMUTATIONS}")
    print(
        f"empirical p-value (>= observed)     = "
        f"{p_value:.6f}"
    )

    # ---------------------------------------------------------------------
    # 13. FEATURE-LEVEL MODEL
    # ---------------------------------------------------------------------

    train = [
        s for s in cyclo_samples
        if s.target_idx in train_targets
    ]

    test = [
        s for s in cyclo_samples
        if s.target_idx in test_targets
    ]

    if len(set(s.label for s in train)) >= 2:

        X_train = np.array(
            [s.features for s in train],
            dtype=float,
        )
        y_train = np.array(
            [s.label for s in train],
            dtype=int,
        )

        model = make_model()
        model.fit(X_train, y_train)

        names = [
            "chi(n)",
            "chi(n+1)",
            "chi(n-1)",
            "chi(4n+1)",
            "chi(G)",

            "cubic_n==0",
            "cubic_n==1",
            "cubic_n==2",

            "cubic_n1==0",
            "cubic_n1==1",
            "cubic_n1==2",

            "cubic_4n1==0",
            "cubic_4n1==1",
            "cubic_4n1==2",
        ]

        ranking = sorted(
            zip(names, model.coef_[0]),
            key=lambda x: abs(x[1]),
            reverse=True,
        )

        print()
        print("13. FEATURE COEFFICIENT RANKING")
        print("-" * 78)

        for name, coeff in ranking:
            print(
                f"{name:18s} {coeff:+.8f}"
            )

    # ---------------------------------------------------------------------
    # 14. TARGET-LEVEL ACCURACY
    # ---------------------------------------------------------------------

    train = [
        s for s in cyclo_samples
        if s.target_idx in train_targets
    ]

    test_targets_sorted = sorted(test_targets)

    print()
    print("14. PER-TARGET HOLDOUT ACCURACY")
    print("-" * 78)

    for target_idx in test_targets_sorted:
        local_train = train

        local_test = [
            s for s in cyclo_samples
            if s.target_idx == target_idx
        ]

        pred = fit_predict(
            local_train,
            local_test,
        )

        acc = accuracy(
            pred,
            local_test,
        )

        print(
            f"target {target_idx:3d}: "
            f"clean={len(local_test):2d} "
            f"accuracy={acc:.6f}"
        )

    # ---------------------------------------------------------------------
    # 15. FINAL DIAGNOSTIC
    # ---------------------------------------------------------------------

    cyclo_loo_mean = statistics.mean(cyclo_loo.values())
    control_loo_mean = statistics.mean(control_loo.values())

    cyclo_2d_mean = statistics.mean(cyclo_2d.values())
    control_2d_mean = statistics.mean(control_2d.values())

    print()
    print("=" * 78)
    print("15. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print("TARGET HOLDOUT")
    print(
        f"  cyclotomic = {cyclo_target_acc:.6f}"
    )
    print(
        f"  control    = {control_target_acc:.6f}"
    )
    print(
        f"  delta      = "
        f"{cyclo_target_acc - control_target_acc:+.6f}"
    )

    print()
    print("REPEATED TARGET HOLDOUT")
    print(
        f"  cyclotomic mean = "
        f"{statistics.mean(cyclo_repeated):.6f}"
    )
    print(
        f"  control mean    = "
        f"{statistics.mean(control_repeated):.6f}"
    )

    print()
    print("UNSEEN-MODULUS HOLDOUT")
    print(
        f"  cyclotomic mean = {cyclo_loo_mean:.6f}"
    )
    print(
        f"  control mean    = {control_loo_mean:.6f}"
    )
    print(
        f"  delta           = "
        f"{cyclo_loo_mean - control_loo_mean:+.6f}"
    )

    print()
    print("UNSEEN-TARGET + UNSEEN-MODULUS")
    print(
        f"  cyclotomic mean = {cyclo_2d_mean:.6f}"
    )
    print(
        f"  control mean    = {control_2d_mean:.6f}"
    )
    print(
        f"  delta           = "
        f"{cyclo_2d_mean - control_2d_mean:+.6f}"
    )

    print()
    print("PERMUTATION")
    print(
        f"  empirical p-value = {p_value:.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)

    print(
        """
NO SIGNAL
---------
Cyclotomic accuracy remains near 0.50 under target and modulus holdouts,
with no stable advantage over controls.

MODULUS-SPECIFIC EFFECT
-----------------------
Target holdout is above 0.50, but leave-one-ell-out accuracy collapses.
This means the model likely learned properties tied to the observed
moduli rather than a transferable orientation law.

TARGET-SPECIFIC EFFECT
----------------------
Leave-one-ell-out remains weak/moderate, but target holdout is substantially
better. This suggests a possible correlation with the target distribution,
not necessarily a universal cyclotomic rule.

CYCLOTOMIC GENERALIZATION
-------------------------
The interesting case is when ALL of the following survive:

    target holdout > 0.5
    modulus holdout > control modulus holdout
    two-dimensional holdout > 0.5
    permutation p-value is small

That would be materially stronger evidence than Experiment 61R.

STRONGEST RESULT
----------------
The strongest outcome would be a cyclotomic advantage that survives both
unseen TARGETS and unseen ELL values.

Only after that should we attempt another CRT reconstruction experiment.

IMPORTANT
---------
This experiment still does not establish a factorization algorithm.
It tests whether an N-only orientation predictor has genuine
out-of-sample generalization.
        """
    )

    total_time = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 62 COMPLETE")
    print(f"total runtime = {total_time:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

