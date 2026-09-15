#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 69
PAIRWISE N-ONLY CHARACTER COMPRESSION / STRICT OUT-OF-SAMPLE TEST
RELATIVE LEGENDRE ORIENTATION FROM CROSS-MODULUS CHARACTER FEATURES
TARGET + MODULUS + 2D HOLDOUT
EXACT SINGLE-FEATURE INFORMATION + TOP-K LOOKUP MODEL
LIGHTWEIGHT WITHIN-TARGET PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Goal
----
Test whether the relative orientation label

    Y(a,b) = sigma_a * sigma_b

contains predictable information in a richer but still strictly N-only
character family evaluated independently at the two moduli.

The experiment deliberately avoids:
  * p or q in any feature
  * s in any feature
  * factor-dependent quantities
  * an oracle root orientation

Local sigma is used only as the LABEL.

A local cell is usable only when

    chi(d1) == chi(d2) != 0

with roots ordered canonically as the two integer roots r1 < r2.

Features are generated from N-only functions:
    n + c
    4n + c
    G(n) + c
    n^2 + n + c
    n^2 + 4n + c

for small integer offsets c.

For each pair (ell_a, ell_b), the feature vector contains:
  * each local character separately
  * the product of the two local characters

The classifier is intentionally simple:
  1. rank single features by exact mutual information on TRAIN targets
  2. select the top K
  3. fit an exact lookup table on their ternary feature tuple
  4. evaluate on unseen targets
  5. repeat for unseen moduli
  6. repeat for unseen targets + unseen moduli

This is a structural information test, not a factorization algorithm.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

SEED = 690069
N_TARGETS = 1200

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

TEST_TARGET_FRACTION = 0.25

# Keep the feature family rich enough to be meaningful, but still cheap.
OFFSETS = (-2, -1, 0, 1, 2)

# Top features used by the exact lookup model.
TOP_K = 4

# Pair sampling cap per family. None = all valid pair samples.
# 1200 targets * 78 pairs is still manageable, but this cap prevents
# accidental runtime explosion on unusually dense data.
MAX_PAIR_SAMPLES = 120_000

# Lightweight permutation null.
N_PERMUTATIONS = 40
PERM_SAMPLE_SIZE = 15_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127,
    307, 331, 631, 1723
]

RAW_CONTROLS = [
    673, 1109, 2309, 2351, 4421, 4561, 4759, 6211,
    7879, 7951, 8273, 8689, 9781
]

# x^2+x+1 has two nontrivial roots iff ell == 1 mod 3.
VALID_CONTROLS = [p for p in RAW_CONTROLS if p % 3 == 1]
EXCLUDED_CONTROLS = [p for p in RAW_CONTROLS if p % 3 != 1]


# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

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


@dataclass(frozen=True)
class LocalState:
    ell: int
    r1: int
    r2: int
    chi_d1: int
    chi_d2: int
    sigma: int | None

    @property
    def usable(self) -> bool:
        return self.sigma in (-1, 1)


@dataclass(frozen=True)
class Sample:
    target_id: int
    ell_a: int
    ell_b: int
    y: int
    features: Tuple[int, ...]


# ---------------------------------------------------------------------------
# PRIME GENERATION
# ---------------------------------------------------------------------------

def prime_sieve(lo: int, hi: int) -> List[int]:
    """Return primes in [lo, hi]."""
    limit = hi
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    root = int(math.isqrt(limit))
    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(lo, hi + 1) if sieve[p]]


def generate_targets(primes: Sequence[int], count: int, rng: random.Random) -> List[Target]:
    """Generate distinct semiprimes from distinct primes."""
    out: List[Target] = []
    seen_n = set()

    attempts = 0
    while len(out) < count:
        attempts += 1
        if attempts > count * 100:
            raise RuntimeError("Could not generate enough distinct targets.")

        p, q = rng.sample(primes, 2)
        if p > q:
            p, q = q, p

        n = p * q
        if n in seen_n:
            continue

        seen_n.add(n)
        out.append(Target(len(out) + 1, p, q))

    return out


# ---------------------------------------------------------------------------
# ROOTS / LEGENDRE
# ---------------------------------------------------------------------------

def roots_cyclotomic(ell: int) -> Tuple[int, int]:
    """
    Return the two roots of x^2+x+1 mod ell.

    This is deliberately brute-forced because the modulus list is tiny
    and small enough that validation is essentially free.
    """
    roots = [x for x in range(1, ell) if (x * x + x + 1) % ell == 0]
    roots = sorted(set(roots))

    if len(roots) != 2:
        raise RuntimeError(
            f"{ell} does not have exactly two roots of x^2+x+1"
        )

    return roots[0], roots[1]


def validate_moduli(moduli: Sequence[int]) -> Dict[int, Tuple[int, int]]:
    roots: Dict[int, Tuple[int, int]] = {}
    failures = []

    print("4. MODULUS VALIDATION")
    print("-" * 78)

    for ell in moduli:
        try:
            r1, r2 = roots_cyclotomic(ell)
            ok = (
                (r1 * r1 + r1 + 1) % ell == 0
                and (r2 * r2 + r2 + 1) % ell == 0
                and (r1 * r2) % ell == 1
            )
        except Exception:
            ok = False
            r1 = r2 = -1

        print(
            f"ell={ell:5d} roots=({r1:5d}, {r2:5d}) "
            f"identity_ok={ok}"
        )

        if ok:
            roots[ell] = (r1, r2)
        else:
            failures.append(ell)

    print(f"validation failures = {len(failures)}")
    print(f"status = {'PASS' if not failures else 'FAIL'}")
    print()

    if failures:
        raise RuntimeError(f"Invalid moduli: {failures}")

    return roots


def legendre_symbol(a: int, p: int) -> int:
    """Return -1, 0, +1 for odd prime p."""
    a %= p
    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1
    if value == p - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre result for a={a}, p={p}: {value}"
    )


# ---------------------------------------------------------------------------
# LOCAL STATES
# ---------------------------------------------------------------------------

def G(n: int) -> int:
    return 16 * n * n + 4 * n + 1


def local_state(
    n: int,
    ell: int,
    roots: Dict[int, Tuple[int, int]],
) -> LocalState:
    r1, r2 = roots[ell]

    d1 = (r1 * r1 - 4 * n) % ell
    d2 = (r2 * r2 - 4 * n) % ell

    c1 = legendre_symbol(d1, ell)
    c2 = legendre_symbol(d2, ell)

    if c1 != 0 and c1 == c2:
        sigma = c1
    else:
        sigma = None

    return LocalState(
        ell=ell,
        r1=r1,
        r2=r2,
        chi_d1=c1,
        chi_d2=c2,
        sigma=sigma,
    )


# ---------------------------------------------------------------------------
# N-ONLY FEATURE FAMILY
# ---------------------------------------------------------------------------

def feature_functions(n: int) -> List[Tuple[str, int]]:
    """
    Return (feature_name, integer_value) functions.

    Values are not themselves the final feature values; they are passed
    through the Legendre symbol at each modulus.
    """
    out: List[Tuple[str, int]] = []

    for c in OFFSETS:
        out.append((f"n{c:+d}", n + c))

    for c in OFFSETS:
        out.append((f"4n{c:+d}", 4 * n + c))

    g = G(n)
    for c in OFFSETS:
        out.append((f"G{c:+d}", g + c))

    for c in OFFSETS:
        out.append((f"n2+n{c:+d}", n * n + n + c))

    for c in OFFSETS:
        out.append((f"n2+4n{c:+d}", n * n + 4 * n + c))

    return out


def local_feature_vector(
    n: int,
    ell: int,
) -> Tuple[Tuple[str, ...], Tuple[int, ...]]:
    names: List[str] = []
    values: List[int] = []

    for name, raw in feature_functions(n):
        names.append(name)
        values.append(legendre_symbol(raw, ell))

    return tuple(names), tuple(values)


# ---------------------------------------------------------------------------
# SAMPLE CONSTRUCTION
# ---------------------------------------------------------------------------

def build_pair_samples(
    targets: Sequence[Target],
    moduli: Sequence[int],
    roots: Dict[int, Tuple[int, int]],
    rng: random.Random,
    max_samples: int | None,
) -> List[Sample]:
    """
    Construct pair samples.

    Feature vector contains:
      local features for ell_a
      local features for ell_b
      pairwise products of corresponding features

    The label is sigma_a * sigma_b.
    """
    samples: List[Sample] = []

    feature_names: Tuple[str, ...] | None = None

    for t in targets:
        states = {
            ell: local_state(t.n, ell, roots)
            for ell in moduli
        }

        usable = [ell for ell in moduli if states[ell].usable]

        for i in range(len(usable)):
            ell_a = usable[i]

            names_a, fa = local_feature_vector(t.n, ell_a)

            for j in range(i + 1, len(usable)):
                ell_b = usable[j]

                names_b, fb = local_feature_vector(t.n, ell_b)

                if feature_names is None:
                    pair_names = (
                        tuple(f"{ell_a}:A:{x}" for x in names_a)
                        + tuple(f"{ell_b}:B:{x}" for x in names_b)
                        + tuple(f"{ell_a}:{ell_b}:P:{x}" for x in names_a)
                    )
                    feature_names = pair_names

                products = tuple(
                    a * b for a, b in zip(fa, fb)
                )

                x = fa + fb + products

                y = states[ell_a].sigma * states[ell_b].sigma

                samples.append(
                    Sample(
                        target_id=t.idx,
                        ell_a=ell_a,
                        ell_b=ell_b,
                        y=y,
                        features=x,
                    )
                )

                if max_samples is not None and len(samples) >= max_samples:
                    rng.shuffle(samples)
                    return samples

    rng.shuffle(samples)
    return samples


# ---------------------------------------------------------------------------
# INFORMATION THEORY
# ---------------------------------------------------------------------------

def entropy_binary(labels: Sequence[int]) -> float:
    if not labels:
        return 0.0

    n = len(labels)
    counts = Counter(labels)

    h = 0.0
    for c in counts.values():
        p = c / n
        if p > 0:
            h -= p * math.log2(p)

    return h


def mutual_information_single(
    values: Sequence[int],
    labels: Sequence[int],
) -> float:
    if not values:
        return 0.0

    n = len(labels)

    joint = Counter(zip(values, labels))
    xv = Counter(values)
    yv = Counter(labels)

    mi = 0.0

    for (x, y), count in joint.items():
        pxy = count / n
        px = xv[x] / n
        py = yv[y] / n

        if pxy > 0:
            mi += pxy * math.log2(pxy / (px * py))

    return mi


# ---------------------------------------------------------------------------
# LOOKUP MODEL
# ---------------------------------------------------------------------------

def select_top_features(
    samples: Sequence[Sample],
    top_k: int,
) -> List[int]:
    if not samples:
        return []

    width = len(samples[0].features)
    labels = [s.y for s in samples]

    scored = []

    for j in range(width):
        values = [s.features[j] for s in samples]
        mi = mutual_information_single(values, labels)
        scored.append((mi, j))

    scored.sort(reverse=True)
    return [j for _, j in scored[:top_k]]


def fit_lookup(
    samples: Sequence[Sample],
    feature_ids: Sequence[int],
) -> Dict[Tuple[int, ...], int]:
    """
    Exact majority lookup on the selected ternary feature tuple.
    """
    counts: Dict[Tuple[int, ...], Counter] = defaultdict(Counter)

    for s in samples:
        key = tuple(s.features[j] for j in feature_ids)
        counts[key][s.y] += 1

    model: Dict[Tuple[int, ...], int] = {}

    for key, counter in counts.items():
        if counter[1] > counter[-1]:
            model[key] = 1
        elif counter[-1] > counter[1]:
            model[key] = -1
        else:
            model[key] = 1

    return model


def predict_lookup(
    model: Dict[Tuple[int, ...], int],
    samples: Sequence[Sample],
    feature_ids: Sequence[int],
) -> List[int]:
    out = []

    # Global fallback based on training majority.
    fallback = 1
    if model:
        total_pos = sum(1 for v in model.values() if v == 1)
        total_neg = sum(1 for v in model.values() if v == -1)
        if total_neg > total_pos:
            fallback = -1

    for s in samples:
        key = tuple(s.features[j] for j in feature_ids)
        out.append(model.get(key, fallback))

    return out


def accuracy(
    samples: Sequence[Sample],
    pred: Sequence[int],
) -> float:
    if not samples:
        return float("nan")

    return sum(
        int(s.y == yhat)
        for s, yhat in zip(samples, pred)
    ) / len(samples)


def balanced_accuracy(
    samples: Sequence[Sample],
    pred: Sequence[int],
) -> float:
    if not samples:
        return float("nan")

    tp = tn = pos = neg = 0

    for s, yhat in zip(samples, pred):
        if s.y == 1:
            pos += 1
            if yhat == 1:
                tp += 1
        else:
            neg += 1
            if yhat == -1:
                tn += 1

    if pos == 0 or neg == 0:
        return float("nan")

    return 0.5 * (tp / pos + tn / neg)


# ---------------------------------------------------------------------------
# SPLITS
# ---------------------------------------------------------------------------

def split_targets(
    targets: Sequence[Target],
    rng: random.Random,
) -> Tuple[set[int], set[int]]:
    ids = [t.idx for t in targets]
    rng.shuffle(ids)

    cut = int(len(ids) * (1.0 - TEST_TARGET_FRACTION))

    train_ids = set(ids[:cut])
    test_ids = set(ids[cut:])

    return train_ids, test_ids


def filter_target_ids(
    samples: Sequence[Sample],
    ids: set[int],
) -> List[Sample]:
    return [s for s in samples if s.target_id in ids]


def filter_modulus_holdout(
    samples: Sequence[Sample],
    train_moduli: set[int],
    test_moduli: set[int],
) -> Tuple[List[Sample], List[Sample]]:
    train = []
    test = []

    for s in samples:
        if s.ell_a in test_moduli or s.ell_b in test_moduli:
            test.append(s)
        elif s.ell_a in train_moduli and s.ell_b in train_moduli:
            train.append(s)

    return train, test


# ---------------------------------------------------------------------------
# PERMUTATION NULL
# ---------------------------------------------------------------------------

def within_target_permutation(
    samples: Sequence[Sample],
    feature_ids: Sequence[int],
    base_train_ids: set[int],
    n_perm: int,
    rng: random.Random,
) -> Tuple[float, float]:
    """
    Lightweight permutation null.

    We sample at most PERM_SAMPLE_SIZE training samples.
    Labels are shuffled within target, preserving each target's label mix.
    """
    train = [
        s for s in samples
        if s.target_id in base_train_ids
    ]

    if len(train) > PERM_SAMPLE_SIZE:
        train = rng.sample(train, PERM_SAMPLE_SIZE)

    grouped: Dict[int, List[int]] = defaultdict(list)
    for i, s in enumerate(train):
        grouped[s.target_id].append(i)

    observed_train_model = fit_lookup(train, feature_ids)

    # Observed is calculated by cross-validation-free lookup accuracy.
    observed = accuracy(
        train,
        predict_lookup(observed_train_model, train, feature_ids),
    )

    null_values = []

    for _ in range(n_perm):
        ys = [s.y for s in train]

        for idxs in grouped.values():
            vals = [ys[i] for i in idxs]
            rng.shuffle(vals)
            for i, value in zip(idxs, vals):
                ys[i] = value

        permuted = [
            Sample(
                target_id=s.target_id,
                ell_a=s.ell_a,
                ell_b=s.ell_b,
                y=ys[i],
                features=s.features,
            )
            for i, s in enumerate(train)
        ]

        model = fit_lookup(permuted, feature_ids)
        pred = predict_lookup(model, permuted, feature_ids)
        null_values.append(accuracy(permuted, pred))

    ge = sum(v >= observed for v in null_values)
    p = (ge + 1) / (len(null_values) + 1)

    return observed, p


# ---------------------------------------------------------------------------
# REPORT HELPERS
# ---------------------------------------------------------------------------

def print_feature_ranking(
    samples: Sequence[Sample],
    feature_count: int,
    feature_labels: Sequence[str],
    top: int = 15,
) -> List[int]:
    labels = [s.y for s in samples]
    ranked = []

    for j in range(feature_count):
        vals = [s.features[j] for s in samples]
        mi = mutual_information_single(vals, labels)
        ranked.append((mi, j))

    ranked.sort(reverse=True)

    print("10. TOP SINGLE N-ONLY FEATURES")
    print("-" * 78)

    for rank, (mi, j) in enumerate(ranked[:top], 1):
        print(
            f"{rank:2d}. {feature_labels[j]:35s} "
            f"MI={mi:.8f}"
        )

    print()

    return [j for _, j in ranked]


def make_pair_feature_names(
    moduli: Sequence[int],
) -> List[str]:
    base_names = [name for name, _ in feature_functions(1)]

    names = []
    for ell in moduli:
        names.extend([f"A ell={ell}:{x}" for x in base_names])
        names.extend([f"B ell={ell}:{x}" for x in base_names])

    # Actual pairwise feature labels depend on the two moduli. The dataset
    # uses fixed column positions only because the first pair establishes
    # the feature layout. For interpretation we use generic positional names.
    half = len(names) // 2
    local_a = names[:half]
    local_b = names[half:]

    return (
        local_a
        + local_b
        + [f"PAIR:{base_names[i]}" for i in range(len(base_names))]
    )


def family_stats(
    samples: Sequence[Sample],
) -> str:
    if not samples:
        return "empty"

    pos = sum(s.y == 1 for s in samples)
    neg = len(samples) - pos
    return (
        f"n={len(samples)} +={pos} -={neg} "
        f"majority={max(pos, neg) / len(samples):.6f}"
    )


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.perf_counter()
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 69")
    print("PAIRWISE N-ONLY CHARACTER COMPRESSION / STRICT OOS TEST")
    print("RELATIVE LEGENDRE ORIENTATION")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("LIGHTWEIGHT WITHIN-TARGET PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    # ----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ----------------------------------------------------------------------
    t = time.perf_counter()
    primes = prime_sieve(PRIME_LO, PRIME_HI)

    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time = {time.perf_counter() - t:.6f}s")
    print()

    # ----------------------------------------------------------------------
    # 2. TARGETS
    # ----------------------------------------------------------------------
    targets = generate_targets(primes, N_TARGETS, rng)

    print("2. TARGET SUMMARY")
    print("-" * 78)
    print(f"total targets = {len(targets)}")
    for target in targets[:24]:
        print(
            f"target {target.idx:4d}: "
            f"p={target.p} q={target.q} "
            f"n={target.n} s={target.s}"
        )
    if len(targets) > 24:
        print("... remaining generated targets omitted")
    print()

    # ----------------------------------------------------------------------
    # 3. MODULI
    # ----------------------------------------------------------------------
    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic       = {CYCLOTOMIC}")
    print(f"valid controls   = {VALID_CONTROLS}")
    print(f"excluded controls= {EXCLUDED_CONTROLS}")
    print()

    all_moduli = CYCLOTOMIC + VALID_CONTROLS
    roots = validate_moduli(all_moduli)

    # ----------------------------------------------------------------------
    # 5. TARGET HOLDOUT
    # ----------------------------------------------------------------------
    train_ids, test_ids = split_targets(targets, rng)

    print("5. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_ids)}")
    print(f"test targets     = {len(test_ids)}")
    print()

    # ----------------------------------------------------------------------
    # 6. DATASET
    # ----------------------------------------------------------------------
    t = time.perf_counter()

    cyclo_samples = build_pair_samples(
        targets,
        CYCLOTOMIC,
        roots,
        rng,
        MAX_PAIR_SAMPLES,
    )

    control_samples = build_pair_samples(
        targets,
        VALID_CONTROLS,
        roots,
        rng,
        MAX_PAIR_SAMPLES,
    )

    print("6. CROSS-MODULUS DATASET")
    print("-" * 78)
    print(f"cyclotomic samples = {len(cyclo_samples)}")
    print(f"control samples    = {len(control_samples)}")
    print(f"construction time  = {time.perf_counter() - t:.6f}s")
    print()

    # Feature layout is fixed in the sample representation because the
    # corresponding local feature positions are identical for every pair.
    base_feature_names = [
        name for name, _ in feature_functions(targets[0].n)
    ]

    feature_count = len(base_feature_names) * 3

    positional_feature_names = (
        [f"A:{x}" for x in base_feature_names]
        + [f"B:{x}" for x in base_feature_names]
        + [f"PAIR:{x}" for x in base_feature_names]
    )

    # ----------------------------------------------------------------------
    # 7. LABEL DISTRIBUTIONS
    # ----------------------------------------------------------------------
    print("7. LABEL DISTRIBUTIONS")
    print("-" * 78)

    for name, samples in (
        ("cyclotomic", cyclo_samples),
        ("control", control_samples),
    ):
        pos = sum(s.y == 1 for s in samples)
        neg = sum(s.y == -1 for s in samples)

        print(
            f"{name:12s} += {pos:7d} -= {neg:7d} "
            f"majority={max(pos, neg)/len(samples):.6f}"
        )

    print()

    # ----------------------------------------------------------------------
    # 8. SINGLE FEATURE INFORMATION
    # ----------------------------------------------------------------------
    print("8. INFORMATION CONTENT")
    print("-" * 78)

    cyclo_train = filter_target_ids(cyclo_samples, train_ids)
    cyclo_test = filter_target_ids(cyclo_samples, test_ids)

    control_train = filter_target_ids(control_samples, train_ids)
    control_test = filter_target_ids(control_samples, test_ids)

    print(f"cyclotomic train = {len(cyclo_train)}")
    print(f"cyclotomic test  = {len(cyclo_test)}")
    print(f"control train    = {len(control_train)}")
    print(f"control test     = {len(control_test)}")
    print()

    cyclo_ranked = print_feature_ranking(
        cyclo_train,
        feature_count,
        positional_feature_names,
    )

    # ----------------------------------------------------------------------
    # 9. TOP-K EXACT LOOKUP
    # ----------------------------------------------------------------------
    print("9. TOP-K EXACT LOOKUP MODEL")
    print("-" * 78)

    selected = cyclo_ranked[:TOP_K]
    print("selected features:")
    for j in selected:
        print(f"  {j:3d}: {positional_feature_names[j]}")

    cyclo_model = fit_lookup(cyclo_train, selected)
    cyclo_pred = predict_lookup(cyclo_model, cyclo_test, selected)

    control_selected = (
        print_feature_ranking(
            control_train,
            feature_count,
            positional_feature_names,
            top=10,
        )[:TOP_K]
    )

    control_model = fit_lookup(control_train, control_selected)
    control_pred = predict_lookup(
        control_model,
        control_test,
        control_selected,
    )

    c_train_acc = accuracy(
        cyclo_train,
        predict_lookup(cyclo_model, cyclo_train, selected),
    )
    c_test_acc = accuracy(cyclo_test, cyclo_pred)
    c_test_bal = balanced_accuracy(cyclo_test, cyclo_pred)

    r_train_acc = accuracy(
        control_train,
        predict_lookup(control_model, control_train, control_selected),
    )
    r_test_acc = accuracy(control_test, control_pred)
    r_test_bal = balanced_accuracy(control_test, control_pred)

    print()
    print(
        f"cyclotomic train={c_train_acc:.6f} "
        f"test={c_test_acc:.6f} "
        f"balanced={c_test_bal:.6f}"
    )
    print(
        f"control    train={r_train_acc:.6f} "
        f"test={r_test_acc:.6f} "
        f"balanced={r_test_bal:.6f}"
    )
    print(
        f"test delta={c_test_acc-r_test_acc:+.6f}"
    )
    print()

    # ----------------------------------------------------------------------
    # 10. UNSEEN MODULUS
    # ----------------------------------------------------------------------
    print("10. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    cyclo_loo = []

    for ell in CYCLOTOMIC:
        train_moduli = set(CYCLOTOMIC) - {ell}
        test_moduli = {ell}

        train_s, test_s = filter_modulus_holdout(
            cyclo_samples,
            train_moduli,
            test_moduli,
        )

        train_s = [s for s in train_s if s.target_id in train_ids]
        test_s = [s for s in test_s if s.target_id in test_ids]

        if not train_s or not test_s:
            acc = float("nan")
        else:
            ranked = select_top_features(train_s, TOP_K)
            model = fit_lookup(train_s, ranked)
            pred = predict_lookup(model, test_s, ranked)
            acc = balanced_accuracy(test_s, pred)

        cyclo_loo.append((ell, acc))

        print(
            f"cyclo ell={ell:5d} "
            f"accuracy={acc:.6f} "
            f"train={len(train_s):6d} "
            f"test={len(test_s):6d}"
        )

    valid = [a for _, a in cyclo_loo if not math.isnan(a)]
    cyclo_loo_mean = statistics.fmean(valid) if valid else float("nan")

    print(f"cyclotomic LOO mean = {cyclo_loo_mean:.6f}")
    print()

    control_loo = []

    for ell in VALID_CONTROLS:
        train_moduli = set(VALID_CONTROLS) - {ell}
        test_moduli = {ell}

        train_s, test_s = filter_modulus_holdout(
            control_samples,
            train_moduli,
            test_moduli,
        )

        train_s = [s for s in train_s if s.target_id in train_ids]
        test_s = [s for s in test_s if s.target_id in test_ids]

        if not train_s or not test_s:
            acc = float("nan")
        else:
            ranked = select_top_features(train_s, TOP_K)
            model = fit_lookup(train_s, ranked)
            pred = predict_lookup(model, test_s, ranked)
            acc = balanced_accuracy(test_s, pred)

        control_loo.append((ell, acc))

        print(
            f"control ell={ell:5d} "
            f"accuracy={acc:.6f} "
            f"train={len(train_s):6d} "
            f"test={len(test_s):6d}"
        )

    valid = [a for _, a in control_loo if not math.isnan(a)]
    control_loo_mean = statistics.fmean(valid) if valid else float("nan")

    print(f"control LOO mean = {control_loo_mean:.6f}")
    print()

    # ----------------------------------------------------------------------
    # 11. TWO-DIMENSIONAL HOLDOUT
    # ----------------------------------------------------------------------
    print("11. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)

    heldout_cyclo = set(sorted(CYCLOTOMIC)[::4])
    heldout_control = set(sorted(VALID_CONTROLS)[::3])

    c_train_2d, c_test_2d = filter_modulus_holdout(
        cyclo_samples,
        set(CYCLOTOMIC) - heldout_cyclo,
        heldout_cyclo,
    )

    c_train_2d = [s for s in c_train_2d if s.target_id in train_ids]
    c_test_2d = [s for s in c_test_2d if s.target_id in test_ids]

    r_train_2d, r_test_2d = filter_modulus_holdout(
        control_samples,
        set(VALID_CONTROLS) - heldout_control,
        heldout_control,
    )

    r_train_2d = [s for s in r_train_2d if s.target_id in train_ids]
    r_test_2d = [s for s in r_test_2d if s.target_id in test_ids]

    if c_train_2d and c_test_2d:
        ranked = select_top_features(c_train_2d, TOP_K)
        model = fit_lookup(c_train_2d, ranked)
        pred = predict_lookup(model, c_test_2d, ranked)
        c2d = balanced_accuracy(c_test_2d, pred)
    else:
        c2d = float("nan")

    if r_train_2d and r_test_2d:
        ranked = select_top_features(r_train_2d, TOP_K)
        model = fit_lookup(r_train_2d, ranked)
        pred = predict_lookup(model, r_test_2d, ranked)
        r2d = balanced_accuracy(r_test_2d, pred)
    else:
        r2d = float("nan")

    print(
        f"cyclotomic 2D={c2d:.6f} "
        f"heldout={sorted(heldout_cyclo)} "
        f"n={len(c_test_2d)}"
    )
    print(
        f"control    2D={r2d:.6f} "
        f"heldout={sorted(heldout_control)} "
        f"n={len(r_test_2d)}"
    )

    if not math.isnan(c2d) and not math.isnan(r2d):
        print(f"2D delta={c2d-r2d:+.6f}")

    print()

    # ----------------------------------------------------------------------
    # 12. EXACT JOINT INFORMATION
    # ----------------------------------------------------------------------
    print("12. EXACT JOINT INFORMATION")
    print("-" * 78)

    def joint_mi(samples: Sequence[Sample], feature_ids: Sequence[int]) -> float:
        if not samples:
            return 0.0

        keys = [
            tuple(s.features[j] for j in feature_ids)
            for s in samples
        ]
        ys = [s.y for s in samples]
        return mutual_information_single(keys, ys)

    c_mi = joint_mi(cyclo_train, selected)
    r_mi = joint_mi(control_train, control_selected)

    print(f"cyclotomic top-{TOP_K} MI = {c_mi:.8f} bits")
    print(f"control    top-{TOP_K} MI = {r_mi:.8f} bits")
    print(
        f"MI delta = {c_mi-r_mi:+.8f} bits"
    )
    print()

    # ----------------------------------------------------------------------
    # 13. LIGHTWEIGHT WITHIN-TARGET PERMUTATION NULL
    # ----------------------------------------------------------------------
    print("13. WITHIN-TARGET PERMUTATION NULL")
    print("-" * 78)

    # Use the selected cyclotomic features.
    observed, pvalue = within_target_permutation(
        cyclo_samples,
        selected,
        train_ids,
        N_PERMUTATIONS,
        rng,
    )

    print(f"observed training lookup accuracy = {observed:.6f}")
    print(f"permutations = {N_PERMUTATIONS}")
    print(f"empirical p-value = {pvalue:.6f}")
    print()

    # ----------------------------------------------------------------------
    # 14. FINAL DIAGNOSTIC
    # ----------------------------------------------------------------------
    print("14. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        f"cyclotomic target holdout = {c_test_acc:.6f}"
    )
    print(
        f"control target holdout    = {r_test_acc:.6f}"
    )
    print(
        f"target delta              = {c_test_acc-r_test_acc:+.6f}"
    )
    print(
        f"cyclotomic LOO mean       = {cyclo_loo_mean:.6f}"
    )
    print(
        f"control LOO mean          = {control_loo_mean:.6f}"
    )
    print(
        f"cyclotomic 2D             = {c2d:.6f}"
    )
    print(
        f"control 2D                = {r2d:.6f}"
    )
    print(
        f"cyclotomic joint MI       = {c_mi:.8f}"
    )
    print(
        f"control joint MI          = {r_mi:.8f}"
    )
    print(
        f"permutation p             = {pvalue:.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)
    print(
        "POSITIVE SIGNAL requires the same direction to survive "
        "target holdout, unseen-modulus holdout, and 2D holdout."
    )
    print(
        "A large training accuracy with weak out-of-sample accuracy "
        "is classified as overfitting."
    )
    print(
        "Cyclotomic advantage is interesting only when it exceeds "
        "the corresponding control behavior."
    )
    print(
        "This experiment still does not establish a factorization "
        "algorithm; it tests whether N-only cross-modulus characters "
        "contain transferable information about relative orientation."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 69 COMPLETE")
    print(f"total runtime = {time.perf_counter() - t0:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

