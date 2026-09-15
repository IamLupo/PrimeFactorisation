#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 75
LOCAL E1 INFORMATION / N-MOD-ELL^K TEST
ZERO-EVENT DEBIASING
STRICT UNSEEN-TARGET VALIDATION
CYCLOTOMIC VS MATCHED CONTROL MODULI
PERMUTATION-ADJUSTED MUTUAL INFORMATION
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Research question
-----------------
Experiment 74 found an apparent N-only E1 advantage, but the best formula was
simply:

    E1_pred(n) = 0

That can be explained by the unusually frequent event E1 == 0 (mod ell) for
small ell. Therefore this experiment does NOT search polynomials.

Instead, for each modulus ell we test:

    n mod ell^k  ->  E1 mod ell

for k = 1, 2, 3.

We measure:

1. global mutual information I(X;Y)
2. zero-event-only information
3. conditional mutual information after removing Y=0
4. exact residue prediction on unseen targets
5. prediction accuracy after a strict target holdout
6. permutation-adjusted MI
7. collision rate of the N-derived state X
8. cyclotomic vs control comparison

The key question is:

    Does n carry information about E1 beyond the trivial E1 == 0 event?

A useful result must survive:
    - unseen targets
    - zero-event removal
    - permutation null
    - comparison with controls
    - sparse-state diagnostics

Important:
This is still an information experiment. It does not claim a factorization
algorithm.

==============================================================================

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
# Configuration
# ---------------------------------------------------------------------------

SEED = 75001

NUM_TARGETS = 120
TRAIN_FRACTION = 0.75

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

# Local information depth.
K_VALUES = (1, 2, 3)

# Permutation count. Keep moderate so the experiment stays fast.
PERMUTATIONS = 200

# We only trust sparse-state MI cautiously when there are enough repeated X
# states. This is diagnostic, not a hard rejection of a result.
MIN_MEAN_BUCKET_SIZE = 1.50

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

CONTROL = [
    673, 4561, 4759, 6211, 7879, 7951, 8689, 9781
]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int
    e1: int


@dataclass(frozen=True)
class Sample:
    target_idx: int
    n: int
    ell: int
    k: int
    x: int
    y: int


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> List[int]:
    """Return all primes in [lo, hi]."""
    if hi < 2 or hi < lo:
        return []

    limit = hi + 1
    sieve = bytearray(b"\x01") * limit
    sieve[:2] = b"\x00\x00"

    root = int(math.isqrt(hi))
    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit:p] = b"\x00" * (
                ((limit - 1 - start) // p) + 1
            )

    return [i for i in range(lo, hi + 1) if sieve[i]]


# ---------------------------------------------------------------------------
# E1
# ---------------------------------------------------------------------------

def compute_e1(p: int, q: int) -> int:
    """
    E1 = (p^2 - 1)(q^2 - 1)(p + q)
    """
    return (p * p - 1) * (q * q - 1) * (p + q)


def verify_e1_identity(p: int, q: int, n: int, s: int, e1: int) -> bool:
    """
    E1 = s * ((n + 1)^2 - s^2)
    """
    rhs = s * ((n + 1) * (n + 1) - s * s)
    return e1 == rhs


# ---------------------------------------------------------------------------
# Target generation
# ---------------------------------------------------------------------------

def generate_targets(primes: Sequence[int], count: int, rng: random.Random) -> List[Target]:
    if len(primes) < 2 * count:
        raise RuntimeError("Not enough primes to generate targets.")

    chosen = rng.sample(list(primes), 2 * count)

    out: List[Target] = []

    for i in range(count):
        p = chosen[2 * i]
        q = chosen[2 * i + 1]

        if p == q:
            raise RuntimeError("Generated duplicate prime pair.")

        if p > q:
            p, q = q, p

        n = p * q
        s = p + q
        e1 = compute_e1(p, q)

        if not verify_e1_identity(p, q, n, s, e1):
            raise RuntimeError(f"E1 identity failure at target {i + 1}")

        out.append(
            Target(
                idx=i + 1,
                p=p,
                q=q,
                n=n,
                s=s,
                e1=e1,
            )
        )

    return out


# ---------------------------------------------------------------------------
# Entropy / mutual information
# ---------------------------------------------------------------------------

def entropy(values: Sequence[int]) -> float:
    if not values:
        return 0.0

    counts = Counter(values)
    total = len(values)

    h = 0.0
    for c in counts.values():
        p = c / total
        h -= p * math.log2(p)

    return h


def mutual_information(x: Sequence[int], y: Sequence[int]) -> float:
    """
    Plug-in empirical mutual information in bits.
    """
    if len(x) != len(y):
        raise ValueError("x/y length mismatch.")
    if not x:
        return 0.0

    n = len(x)

    joint = Counter(zip(x, y))
    cx = Counter(x)
    cy = Counter(y)

    mi = 0.0

    for (vx, vy), cxy in joint.items():
        pxy = cxy / n
        px = cx[vx] / n
        py = cy[vy] / n

        mi += pxy * math.log2(pxy / (px * py))

    return mi


def conditional_mi_given_nonzero_y(
    x: Sequence[int],
    y: Sequence[int],
) -> float:
    """
    I(X;Y | Y != 0).
    """
    filtered = [
        (vx, vy)
        for vx, vy in zip(x, y)
        if vy != 0
    ]

    if not filtered:
        return 0.0

    fx, fy = zip(*filtered)
    return mutual_information(fx, fy)


def zero_indicator_mi(
    x: Sequence[int],
    y: Sequence[int],
) -> float:
    """
    I(X; 1[Y == 0]).
    """
    z = [1 if v == 0 else 0 for v in y]
    return mutual_information(x, z)


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

def majority_prediction(values: Sequence[int]) -> int:
    if not values:
        raise RuntimeError("Cannot compute majority of empty sequence.")

    counts = Counter(values)

    # Deterministic tie-break.
    best_count = max(counts.values())
    best = min(v for v, c in counts.items() if c == best_count)
    return best


def train_lookup_model(samples: Sequence[Sample]) -> Tuple[Dict[int, int], int]:
    """
    State -> most common Y.

    Global fallback is the global training-mode Y.
    """
    buckets: Dict[int, List[int]] = defaultdict(list)

    for s in samples:
        buckets[s.x].append(s.y)

    lookup = {
        x: majority_prediction(ys)
        for x, ys in buckets.items()
    }

    fallback = majority_prediction([s.y for s in samples])

    return lookup, fallback


def score_lookup_model(
    lookup: Dict[int, int],
    fallback: int,
    samples: Sequence[Sample],
) -> Tuple[float, float]:
    """
    Returns:
        exact accuracy
        balanced accuracy (when both classes exist)
    """
    if not samples:
        return float("nan"), float("nan")

    correct = 0
    pred = []

    for s in samples:
        yhat = lookup.get(s.x, fallback)
        pred.append(yhat)
        if yhat == s.y:
            correct += 1

    accuracy = correct / len(samples)

    labels = sorted(set(s.y for s in samples))
    if len(labels) < 2:
        return accuracy, float("nan")

    recalls = []
    for label in labels:
        idx = [i for i, s in enumerate(samples) if s.y == label]
        if not idx:
            continue

        hit = sum(pred[i] == label for i in idx)
        recalls.append(hit / len(idx))

    balanced = statistics.fmean(recalls) if recalls else float("nan")
    return accuracy, balanced


# ---------------------------------------------------------------------------
# State-collision diagnostics
# ---------------------------------------------------------------------------

def state_diagnostics(samples: Sequence[Sample]) -> Dict[str, float]:
    if not samples:
        return {
            "states": 0.0,
            "collision_rate": float("nan"),
            "mean_bucket": float("nan"),
            "median_bucket": float("nan"),
            "singleton_fraction": float("nan"),
        }

    counts = Counter(s.x for s in samples)
    sizes = list(counts.values())

    singleton_fraction = sum(1 for s in sizes if s == 1) / len(sizes)

    repeated = sum(1 for s in sizes if s > 1)
    collision_rate = repeated / len(sizes)

    return {
        "states": float(len(counts)),
        "collision_rate": collision_rate,
        "mean_bucket": statistics.fmean(sizes),
        "median_bucket": statistics.median(sizes),
        "singleton_fraction": singleton_fraction,
    }


# ---------------------------------------------------------------------------
# Permutation null
# ---------------------------------------------------------------------------

def permuted_mi_values(
    samples: Sequence[Sample],
    rng: random.Random,
    count: int,
) -> List[float]:
    """
    Permute Y within the same modulus.

    Since every sample set passed here is already a single modulus/depth,
    this preserves the marginal Y distribution while destroying the
    N-to-E1 association.
    """
    if not samples:
        return []

    x = [s.x for s in samples]
    y = [s.y for s in samples]

    out = []

    for _ in range(count):
        yp = list(y)
        rng.shuffle(yp)
        out.append(mutual_information(x, yp))

    return out


def empirical_p_value(observed: float, null_values: Sequence[float]) -> float:
    if not null_values:
        return float("nan")

    extreme = sum(v >= observed - 1e-15 for v in null_values)

    return (extreme + 1) / (len(null_values) + 1)


# ---------------------------------------------------------------------------
# Dataset construction
# ---------------------------------------------------------------------------

def build_samples(
    targets: Sequence[Target],
    moduli: Sequence[int],
    k: int,
) -> Dict[int, List[Sample]]:
    """
    Returns samples grouped by modulus.
    """
    out: Dict[int, List[Sample]] = {ell: [] for ell in moduli}

    for t in targets:
        for ell in moduli:
            modulus = ell ** k

            x = t.n % modulus
            y = t.e1 % ell

            out[ell].append(
                Sample(
                    target_idx=t.idx,
                    n=t.n,
                    ell=ell,
                    k=k,
                    x=x,
                    y=y,
                )
            )

    return out


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def summarize_family(
    name: str,
    samples_by_ell: Dict[int, List[Sample]],
    train_ids: set[int],
    test_ids: set[int],
    rng: random.Random,
) -> None:

    print()
    print("=" * 78)
    print(f"{name}")
    print("=" * 78)

    all_raw_mi = []
    all_zero_mi = []
    all_cond_mi = []
    all_test_acc = []
    all_test_balanced = []
    all_p = []

    for ell, samples in samples_by_ell.items():

        train = [s for s in samples if s.target_idx in train_ids]
        test = [s for s in samples if s.target_idx in test_ids]

        if not train or not test:
            continue

        x_train = [s.x for s in train]
        y_train = [s.y for s in train]

        raw_mi = mutual_information(x_train, y_train)
        zero_mi = zero_indicator_mi(x_train, y_train)
        cond_mi = conditional_mi_given_nonzero_y(x_train, y_train)

        lookup, fallback = train_lookup_model(train)
        test_acc, test_bal = score_lookup_model(
            lookup,
            fallback,
            test,
        )

        null = permuted_mi_values(
            train,
            rng,
            PERMUTATIONS,
        )

        p = empirical_p_value(raw_mi, null)

        diag = state_diagnostics(train)

        zero_rate = sum(v == 0 for v in y_train) / len(y_train)

        print(
            f"ell={ell:5d} "
            f"MI={raw_mi:.6f} "
            f"zeroMI={zero_mi:.6f} "
            f"condMI={cond_mi:.6f} "
            f"zero_rate={zero_rate:.4f} "
            f"test_acc={test_acc:.6f} "
            f"test_bal={test_bal:.6f} "
            f"p={p:.6f} "
            f"states={int(diag['states']):3d} "
            f"mean_bucket={diag['mean_bucket']:.2f}"
        )

        all_raw_mi.append(raw_mi)
        all_zero_mi.append(zero_mi)
        all_cond_mi.append(cond_mi)
        all_test_acc.append(test_acc)

        if not math.isnan(test_bal):
            all_test_balanced.append(test_bal)

        all_p.append(p)

    if all_raw_mi:
        print()
        print(f"{name} SUMMARY")
        print(f"  mean raw MI          = {statistics.fmean(all_raw_mi):.6f} bits")
        print(f"  mean zero-event MI   = {statistics.fmean(all_zero_mi):.6f} bits")
        print(f"  mean nonzero-cond MI = {statistics.fmean(all_cond_mi):.6f} bits")
        print(f"  mean OOS accuracy    = {statistics.fmean(all_test_acc):.6f}")

        if all_test_balanced:
            print(
                f"  mean OOS balanced    = "
                f"{statistics.fmean(all_test_balanced):.6f}"
            )

        print(
            f"  fraction p<=0.05    = "
            f"{sum(p <= 0.05 for p in all_p) / len(all_p):.6f}"
        )


def family_aggregate(
    samples_by_ell: Dict[int, List[Sample]],
    train_ids: set[int],
    test_ids: set[int],
    rng: random.Random,
) -> Dict[str, float]:

    raw = []
    zero = []
    cond = []
    acc = []
    bal = []
    pvals = []

    for ell, samples in samples_by_ell.items():
        train = [s for s in samples if s.target_idx in train_ids]
        test = [s for s in samples if s.target_idx in test_ids]

        if not train or not test:
            continue

        x = [s.x for s in train]
        y = [s.y for s in train]

        raw.append(mutual_information(x, y))
        zero.append(zero_indicator_mi(x, y))
        cond.append(conditional_mi_given_nonzero_y(x, y))

        lookup, fallback = train_lookup_model(train)
        a, b = score_lookup_model(lookup, fallback, test)
        acc.append(a)

        if not math.isnan(b):
            bal.append(b)

        null = permuted_mi_values(train, rng, PERMUTATIONS)
        pvals.append(empirical_p_value(raw[-1], null))

    return {
        "raw_mi": statistics.fmean(raw) if raw else float("nan"),
        "zero_mi": statistics.fmean(zero) if zero else float("nan"),
        "cond_mi": statistics.fmean(cond) if cond else float("nan"),
        "acc": statistics.fmean(acc) if acc else float("nan"),
        "balanced": statistics.fmean(bal) if bal else float("nan"),
        "p_fraction": (
            sum(p <= 0.05 for p in pvals) / len(pvals)
            if pvals else float("nan")
        ),
    }


def main() -> None:
    rng = random.Random(SEED)

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 75")
    print("LOCAL E1 INFORMATION / N-MOD-ELL^K TEST")
    print("ZERO-EVENT DEBIASING")
    print("STRICT UNSEEN-TARGET VALIDATION")
    print("CYCLOTOMIC VS MATCHED CONTROL MODULI")
    print("PERMUTATION-ADJUSTED MUTUAL INFORMATION")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # 1. Prime population
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(PRIME_LO, PRIME_HI)

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {time.perf_counter() - t0:.6f}s")

    if len(primes) < 2 * NUM_TARGETS:
        raise RuntimeError("Prime population too small.")

    # -----------------------------------------------------------------------
    # 2. Targets
    # -----------------------------------------------------------------------

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
            f"target {t.idx:4d}: "
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining generated targets omitted")

    # -----------------------------------------------------------------------
    # 3. Train/test split by target
    # -----------------------------------------------------------------------

    ids = [t.idx for t in targets]
    rng.shuffle(ids)

    split = int(NUM_TARGETS * TRAIN_FRACTION)

    train_ids = set(ids[:split])
    test_ids = set(ids[split:])

    print()
    print("3. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_ids)}")
    print(f"test targets     = {len(test_ids)}")

    # -----------------------------------------------------------------------
    # 4. Identity validation
    # -----------------------------------------------------------------------

    failures = 0

    for t in targets:
        if not verify_e1_identity(t.p, t.q, t.n, t.s, t.e1):
            failures += 1

    print()
    print("4. E1 IDENTITY VALIDATION")
    print("-" * 78)
    print(f"identity failures = {failures}")
    print(f"status = {'PASS' if failures == 0 else 'FAIL'}")

    if failures:
        raise RuntimeError("E1 identity validation failed.")

    # -----------------------------------------------------------------------
    # 5. Run each k
    # -----------------------------------------------------------------------

    family_results = {}

    for k in K_VALUES:

        print()
        print("=" * 78)
        print(f"5. LOCAL INFORMATION DEPTH k={k}")
        print("-" * 78)

        print(
            "X = n mod ell^k"
        )
        print(
            "Y = E1 mod ell"
        )

        cyc_samples = build_samples(
            targets,
            CYCLOTOMIC,
            k,
        )

        ctrl_samples = build_samples(
            targets,
            CONTROL,
            k,
        )

        print()
        print("CYCLOTOMIC")
        summarize_family(
            "CYCLOTOMIC",
            cyc_samples,
            train_ids,
            test_ids,
            rng,
        )

        print()
        print("CONTROL")
        summarize_family(
            "CONTROL",
            ctrl_samples,
            train_ids,
            test_ids,
            rng,
        )

        cyc_agg = family_aggregate(
            cyc_samples,
            train_ids,
            test_ids,
            rng,
        )

        ctrl_agg = family_aggregate(
            ctrl_samples,
            train_ids,
            test_ids,
            rng,
        )

        mi_delta = cyc_agg["raw_mi"] - ctrl_agg["raw_mi"]
        cond_delta = cyc_agg["cond_mi"] - ctrl_agg["cond_mi"]
        acc_delta = cyc_agg["acc"] - ctrl_agg["acc"]

        print()
        print("DEPTH COMPARISON")
        print(f"  cyclotomic mean MI          = {cyc_agg['raw_mi']:.6f}")
        print(f"  control mean MI             = {ctrl_agg['raw_mi']:.6f}")
        print(f"  MI delta                    = {mi_delta:+.6f}")
        print()
        print(
            f"  cyclotomic zero-event MI    = "
            f"{cyc_agg['zero_mi']:.6f}"
        )
        print(
            f"  control zero-event MI       = "
            f"{ctrl_agg['zero_mi']:.6f}"
        )
        print(
            f"  zero-event MI delta         = "
            f"{cyc_agg['zero_mi'] - ctrl_agg['zero_mi']:+.6f}"
        )
        print()
        print(
            f"  cyclotomic conditional MI   = "
            f"{cyc_agg['cond_mi']:.6f}"
        )
        print(
            f"  control conditional MI      = "
            f"{ctrl_agg['cond_mi']:.6f}"
        )
        print(
            f"  conditional MI delta        = "
            f"{cond_delta:+.6f}"
        )
        print()
        print(
            f"  cyclotomic OOS accuracy     = "
            f"{cyc_agg['acc']:.6f}"
        )
        print(
            f"  control OOS accuracy        = "
            f"{ctrl_agg['acc']:.6f}"
        )
        print(
            f"  OOS accuracy delta          = "
            f"{acc_delta:+.6f}"
        )

        family_results[k] = {
            "cyc": cyc_agg,
            "ctrl": ctrl_agg,
        }

    # -----------------------------------------------------------------------
    # 6. Cross-depth comparison
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CROSS-DEPTH SUMMARY")
    print("-" * 78)

    print(
        "k | C_MI       R_MI       delta_MI   "
        "C_condMI    R_condMI    delta_cond   "
        "C_acc      R_acc"
    )

    for k in K_VALUES:
        c = family_results[k]["cyc"]
        r = family_results[k]["ctrl"]

        print(
            f"{k:1d} | "
            f"{c['raw_mi']:.8f} "
            f"{r['raw_mi']:.8f} "
            f"{c['raw_mi'] - r['raw_mi']:+.8f}   "
            f"{c['cond_mi']:.8f} "
            f"{r['cond_mi']:.8f} "
            f"{c['cond_mi'] - r['cond_mi']:+.8f}   "
            f"{c['acc']:.6f} "
            f"{r['acc']:.6f}"
        )

    # -----------------------------------------------------------------------
    # 7. Principal diagnostic
    # -----------------------------------------------------------------------

    best_k = max(
        K_VALUES,
        key=lambda kk: (
            family_results[kk]["cyc"]["cond_mi"]
            - family_results[kk]["ctrl"]["cond_mi"]
        ),
    )

    best_c = family_results[best_k]["cyc"]
    best_r = family_results[best_k]["ctrl"]

    print()
    print("=" * 78)
    print("7. PRINCIPAL DIAGNOSTIC")
    print("-" * 78)

    print(f"best conditional-information depth = k={best_k}")
    print(
        f"cyclotomic conditional MI = "
        f"{best_c['cond_mi']:.8f} bits"
    )
    print(
        f"control conditional MI    = "
        f"{best_r['cond_mi']:.8f} bits"
    )
    print(
        f"conditional MI advantage   = "
        f"{best_c['cond_mi'] - best_r['cond_mi']:+.8f} bits"
    )
    print()
    print(
        f"cyclotomic OOS exact rate   = "
        f"{best_c['acc']:.6f}"
    )
    print(
        f"control OOS exact rate      = "
        f"{best_r['acc']:.6f}"
    )
    print(
        f"OOS accuracy advantage      = "
        f"{best_c['acc'] - best_r['acc']:+.6f}"
    )

    print()
    print("INTERPRETATION")
    print("-" * 78)

    cond_adv = (
        best_c["cond_mi"] - best_r["cond_mi"]
    )
    acc_adv = (
        best_c["acc"] - best_r["acc"]
    )

    if cond_adv <= 0 and acc_adv <= 0:
        print(
            "NO NONZERO E1 SIGNAL:"
        )
        print(
            "After removing E1==0, n mod ell^k provides no "
            "clear cyclotomic advantage."
        )
    elif cond_adv > 0 and acc_adv <= 0:
        print(
            "INFORMATION WITHOUT PREDICTIVE POWER:"
        )
        print(
            "Cyclotomic conditional information is higher, but "
            "the effect does not produce useful unseen-target residue prediction."
        )
    elif cond_adv <= 0 and acc_adv > 0:
        print(
            "PREDICTION EFFECT WITHOUT INFORMATION SEPARATION:"
        )
        print(
            "This requires deeper scrutiny; it may reflect class imbalance "
            "or finite-sample behavior."
        )
    else:
        print(
            "POTENTIAL LOCAL E1 SIGNAL:"
        )
        print(
            "Both conditional information and unseen-target prediction "
            "favor the cyclotomic family."
        )

    print()
    print(
        "CRITICAL CHECK:"
    )
    print(
        "Compare conditional MI, not just raw MI."
    )
    print(
        "A positive raw MI caused entirely by E1==0 is not considered a "
        "successful discovery."
    )

    print()
    print(
        "A genuinely interesting result would require:"
    )
    print(
        "  1. positive conditional MI"
    )
    print(
        "  2. cyclotomic > control"
    )
    print(
        "  3. unseen-target exact prediction > matched baseline"
    )
    print(
        "  4. permutation p-values small for multiple moduli"
    )
    print(
        "  5. the result to remain visible for k=1 before using larger k"
    )

    total = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 75 COMPLETE")
    print(f"total runtime = {total:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

