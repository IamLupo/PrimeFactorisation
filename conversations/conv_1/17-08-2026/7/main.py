#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 77
ANTI-PREDICTION / SIGN-INVERSION TEST
DIRECT VS INVERTED N-ONLY CHARACTER PREDICTION
REUSABLE SIGNATURE STATES
TARGET + MODULUS + 2D HOLDOUT
WITHIN-TARGET PERMUTATION NULL
OCCUPANCY-MATCHED DIAGNOSTICS
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Research question
------------------
Experiment 76 produced a cyclotomic normalized-character OOS accuracy below
50%, while its permutation p-value was small.

Experiment 77 asks a very specific question:

    Is that apparently reproducible effect actually an INVERTED relation?

For every predictor we compare:

    DIRECT      : predict Y_hat = +Y_model
    INVERTED    : predict Y_hat = -Y_model

where the target label is:

    Y = Legendre(E1 mod ell) in {-1, +1}

Zero E1-character cells are excluded from the prediction task.

The experiment uses a fixed N-only signature. No target factor information
is used in the predictor.

Outputs
-------
1. direct accuracy
2. inverted accuracy
3. balanced accuracy
4. direct - inverted comparison
5. target holdout
6. leave-one-modulus-out
7. 2D target+modulus holdout
8. within-target permutation null
9. occupancy-matched state analysis
10. per-modulus diagnostics

Interpretation
--------------
A genuine inversion effect would require:

    inverted > direct
    inverted > 0.5
    advantage survives unseen targets
    advantage survives unseen moduli
    advantage survives 2D holdout
    controls do NOT show the same effect

This is a diagnostic experiment, not a factorization algorithm.
"""

from __future__ import annotations

import math
import random
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

SEED = 77001

NUM_TARGETS = 120
TRAIN_TARGETS = 90
TEST_TARGETS = NUM_TARGETS - TRAIN_TARGETS

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

CONTROLS = [
    673, 4561, 4759, 6211,
    7879, 7951, 8689, 9781
]

# Two-dimensional holdout: unseen targets AND unseen moduli.
CYCLO_2D_HELDOUT = [61, 127, 307]
CONTROL_2D_HELDOUT = [6211, 7951]

PERMUTATIONS = 250

# Occupancy thresholds.
OCCUPANCY_THRESHOLDS = [1, 2, 3, 5]

# Signature definition.
# All entries are N-only and modulus-local.
FEATURE_OFFSETS = (
    0,       # n
    1,       # n+1
    -1,      # n-1
    2,       # n+2
    -2,      # n-2
    3,       # n+3
    -3,      # n-3
)

# Polynomial characters whose values are fed through Legendre symbols.
POLYNOMIAL_NAMES = (
    "n",
    "n+1",
    "n-1",
    "n2+n+1",
    "n2+4n+1",
    "4n+1",
    "4n-1",
)

# Maximum number of distinct signature states retained for sanity checks.
# We do not hash/collide states; this only limits diagnostic output.
MAX_PRINT_STATES = 12


# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

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
    family: str
    signature: Tuple[int, ...]
    y: int


@dataclass
class EvalResult:
    total: int
    covered: int
    direct_acc: float
    inverted_acc: float
    direct_bal: float
    inverted_bal: float
    direct_better: bool
    inversion_gain: float
    direct_gain_vs_chance: float
    inverted_gain_vs_chance: float


# ---------------------------------------------------------------------------
# PRIME GENERATION
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> List[int]:
    """
    Inclusive prime sieve.
    """
    if hi < 2 or hi < lo:
        return []

    size = hi + 1
    is_prime = bytearray(b"\x01") * size
    is_prime[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))
    for p in range(2, limit + 1):
        if is_prime[p]:
            start = p * p
            step_count = ((hi - start) // p) + 1
            is_prime[start: hi + 1: p] = b"\x00" * step_count

    return [x for x in range(lo, hi + 1) if is_prime[x]]


def generate_targets(
    prime_pool: Sequence[int],
    count: int,
    rng: random.Random,
) -> List[Target]:
    """
    Generate semiprimes from distinct primes.

    Keeps p,q in the same broad size region as prior Kappa experiments.
    """
    targets: List[Target] = []
    used_pairs = set()

    if len(prime_pool) < 2 * count:
        raise RuntimeError("Prime population too small for requested targets.")

    max_attempts = count * 100
    attempts = 0

    while len(targets) < count and attempts < max_attempts:
        attempts += 1

        p = prime_pool[rng.randrange(len(prime_pool))]
        q = prime_pool[rng.randrange(len(prime_pool))]

        if p == q:
            continue

        a, b = sorted((p, q))
        if (a, b) in used_pairs:
            continue

        used_pairs.add((a, b))
        s = a + b
        n = a * b

        targets.append(
            Target(
                idx=len(targets) + 1,
                p=a,
                q=b,
                n=n,
                s=s,
            )
        )

    if len(targets) != count:
        raise RuntimeError(
            f"Could only generate {len(targets)} targets; required {count}."
        )

    return targets


# ---------------------------------------------------------------------------
# NUMBER THEORY
# ---------------------------------------------------------------------------

def legendre_symbol(a: int, p: int) -> int:
    """
    Legendre symbol for odd prime p.

    Returns:
        -1, 0, +1
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
        f"Unexpected Legendre result for a={a}, p={p}: {value}"
    )


def cubic_character_class(a: int, p: int) -> int:
    """
    For p == 1 mod 3, return the cubic character class:

        0  if a == 0 mod p
        1  if a is a cubic residue
        2  if a belongs to one of the two nontrivial cubic classes

    This is only a compact 3-state diagnostic, not a complex-valued cubic
    character.
    """
    a %= p

    if a == 0:
        return 0

    if (p - 1) % 3 != 0:
        raise ValueError(
            f"p={p} does not support cubic character class."
        )

    exponent = (p - 1) // 3
    z = pow(a, exponent, p)

    if z == 1:
        return 0

    # For a chosen primitive cube-root omega:
    # z = omega or omega^2 => classes 1 or 2.
    omega = find_cube_root(p)
    if z == omega:
        return 1
    if z == (omega * omega) % p:
        return 2

    # Defensive fallback; should not occur for p == 1 mod 3.
    raise RuntimeError(
        f"Unexpected cubic class for a={a}, p={p}, z={z}"
    )


def find_cube_root(p: int) -> int:
    """
    Return a nontrivial cube root of unity modulo prime p == 1 mod 3.

    Deterministic small search is sufficient because all experiment moduli
    are tiny.
    """
    if p % 3 != 1:
        raise ValueError(f"p={p} is not 1 mod 3")

    for z in range(2, p):
        if pow(z, 3, p) == 1 and z != 1:
            return z

    raise RuntimeError(f"No nontrivial cube root found for p={p}.")


# ---------------------------------------------------------------------------
# E1
# ---------------------------------------------------------------------------

def e1_from_factors(p: int, q: int) -> int:
    return (p * p - 1) * (q * q - 1) * (p + q)


def e1_from_n_s(n: int, s: int) -> int:
    return s * ((n + 1) ** 2 - s * s)


def validate_e1(t: Target) -> bool:
    e_factor = e1_from_factors(t.p, t.q)
    e_s = e1_from_n_s(t.n, t.s)

    cubic = (
        t.s ** 3
        - (t.n + 1) ** 2 * t.s
        + e_factor
    )

    return e_factor == e_s and cubic == 0


# ---------------------------------------------------------------------------
# ROOTS OF x^2 + x + 1
# ---------------------------------------------------------------------------

def cyclotomic_roots(ell: int) -> Tuple[int, int]:
    """
    Find the two nontrivial roots of x^2+x+1 mod ell.
    """
    roots = [
        x for x in range(1, ell)
        if (x * x + x + 1) % ell == 0
    ]

    if len(roots) != 2:
        raise RuntimeError(
            f"{ell} does not have exactly two roots of x^2+x+1"
        )

    return roots[0], roots[1]


def validate_moduli(moduli: Sequence[int]) -> None:
    for ell in moduli:
        r1, r2 = cyclotomic_roots(ell)
        if r1 == r2:
            raise RuntimeError(f"Repeated roots at ell={ell}")

        if (r1 * r1 + r1 + 1) % ell != 0:
            raise RuntimeError(f"Bad first root for ell={ell}")

        if (r2 * r2 + r2 + 1) % ell != 0:
            raise RuntimeError(f"Bad second root for ell={ell}")


# ---------------------------------------------------------------------------
# N-ONLY SIGNATURE
# ---------------------------------------------------------------------------

def normalized_poly_values(n: int) -> Tuple[int, ...]:
    """
    Polynomial family deliberately fixed before evaluation.

    These are cheap N-only transforms. The resulting values are fed through
    Legendre symbols modulo ell.
    """
    values = (
        n,
        n + 1,
        n - 1,
        n * n + n + 1,
        n * n + 4 * n + 1,
        4 * n + 1,
        4 * n - 1,
    )
    return values


def n_only_signature(n: int, ell: int) -> Tuple[int, ...]:
    """
    Signature consists of:

      7 Legendre-character states
      + 7 cubic-character classes when ell == 1 mod 3

    The cubic part is included because Experiment 76 used Legendre + cubic
    state compression.

    For controls/cyclotomic primes here, ell == 1 mod 3 holds.
    """
    values = normalized_poly_values(n)

    leg = tuple(
        legendre_symbol(v, ell)
        for v in values
    )

    cub = tuple(
        cubic_character_class(v, ell)
        for v in values
    )

    return leg + cub


# ---------------------------------------------------------------------------
# LABEL
# ---------------------------------------------------------------------------

def e1_character(e1: int, ell: int) -> int:
    y = legendre_symbol(e1, ell)

    # Prediction task requires a binary label.
    if y == 0:
        raise ValueError("E1 character is zero; sample is degenerate.")

    return y


# ---------------------------------------------------------------------------
# DATASET
# ---------------------------------------------------------------------------

def build_samples(
    targets: Sequence[Target],
    moduli: Sequence[int],
    family: str,
) -> List[Sample]:
    samples: List[Sample] = []

    for t in targets:
        e1 = e1_from_factors(t.p, t.q)

        for ell in moduli:
            y = legendre_symbol(e1, ell)

            if y == 0:
                # Exclude zero-character cells from the binary prediction task.
                continue

            sig = n_only_signature(t.n, ell)

            samples.append(
                Sample(
                    target_id=t.idx,
                    ell=ell,
                    family=family,
                    signature=sig,
                    y=y,
                )
            )

    return samples


# ---------------------------------------------------------------------------
# STATE MODEL
# ---------------------------------------------------------------------------

class StateMajorityModel:
    """
    Simple state lookup model.

    Training:
        state -> majority target label

    Prediction:
        known state -> majority label
        unknown state -> None
    """

    def __init__(self) -> None:
        self.counts: Dict[Tuple[int, ...], Counter] = {}

    def fit(self, samples: Sequence[Sample]) -> None:
        buckets: Dict[Tuple[int, ...], Counter] = defaultdict(Counter)

        for s in samples:
            buckets[s.signature][s.y] += 1

        self.counts = dict(buckets)

    def predict(self, sample: Sample) -> Optional[int]:
        c = self.counts.get(sample.signature)

        if not c:
            return None

        if c[1] > c[-1]:
            return 1
        if c[-1] > c[1]:
            return -1

        # Fixed deterministic tie-break.
        return 1

    def occupancy(self, signature: Tuple[int, ...]) -> int:
        c = self.counts.get(signature)
        return 0 if c is None else sum(c.values())


# ---------------------------------------------------------------------------
# METRICS
# ---------------------------------------------------------------------------

def accuracy(y_true: Sequence[int], y_pred: Sequence[int]) -> float:
    if not y_true:
        return float("nan")
    return sum(a == b for a, b in zip(y_true, y_pred)) / len(y_true)


def balanced_accuracy(y_true: Sequence[int], y_pred: Sequence[int]) -> float:
    positives = [(a, b) for a, b in zip(y_true, y_pred) if a == 1]
    negatives = [(a, b) for a, b in zip(y_true, y_pred) if a == -1]

    if not positives or not negatives:
        return float("nan")

    tpr = sum(a == b for a, b in positives) / len(positives)
    tnr = sum(a == b for a, b in negatives) / len(negatives)

    return 0.5 * (tpr + tnr)


def mean_or_nan(values: Iterable[float]) -> float:
    vals = [x for x in values if not math.isnan(x)]
    return statistics.fmean(vals) if vals else float("nan")


def sign_inverse(y: int) -> int:
    return -y


def collect_predictions(
    model: StateMajorityModel,
    samples: Sequence[Sample],
    min_occupancy: int = 1,
) -> Tuple[List[int], List[int], List[int]]:
    """
    Returns:
        y_true
        direct_predictions
        inverse_predictions

    Unknown or insufficient-occupancy states are omitted.
    """
    yt: List[int] = []
    pd: List[int] = []
    pi: List[int] = []

    for s in samples:
        if model.occupancy(s.signature) < min_occupancy:
            continue

        pred = model.predict(s)
        if pred is None:
            continue

        yt.append(s.y)
        pd.append(pred)
        pi.append(-pred)

    return yt, pd, pi


def evaluate_model(
    model: StateMajorityModel,
    samples: Sequence[Sample],
    min_occupancy: int = 1,
) -> EvalResult:
    yt, pd, pi = collect_predictions(
        model,
        samples,
        min_occupancy=min_occupancy,
    )

    total = len(samples)
    covered = len(yt)

    if covered == 0:
        return EvalResult(
            total=total,
            covered=0,
            direct_acc=float("nan"),
            inverted_acc=float("nan"),
            direct_bal=float("nan"),
            inverted_bal=float("nan"),
            direct_better=False,
            inversion_gain=float("nan"),
            direct_gain_vs_chance=float("nan"),
            inverted_gain_vs_chance=float("nan"),
        )

    da = accuracy(yt, pd)
    ia = accuracy(yt, pi)
    db = balanced_accuracy(yt, pd)
    ib = balanced_accuracy(yt, pi)

    return EvalResult(
        total=total,
        covered=covered,
        direct_acc=da,
        inverted_acc=ia,
        direct_bal=db,
        inverted_bal=ib,
        direct_better=da >= ia,
        inversion_gain=ia - da,
        direct_gain_vs_chance=da - 0.5,
        inverted_gain_vs_chance=ia - 0.5,
    )


# ---------------------------------------------------------------------------
# TARGET / MODULUS SPLITS
# ---------------------------------------------------------------------------

def split_targets(
    targets: Sequence[Target],
    rng: random.Random,
) -> Tuple[List[int], List[int]]:
    ids = [t.idx for t in targets]
    rng.shuffle(ids)

    train_ids = ids[:TRAIN_TARGETS]
    test_ids = ids[TRAIN_TARGETS:]

    return sorted(train_ids), sorted(test_ids)


def by_target(
    samples: Sequence[Sample],
    ids: Sequence[int],
) -> List[Sample]:
    allowed = set(ids)
    return [s for s in samples if s.target_id in allowed]


def by_modulus(
    samples: Sequence[Sample],
    moduli: Sequence[int],
) -> List[Sample]:
    allowed = set(moduli)
    return [s for s in samples if s.ell in allowed]


# ---------------------------------------------------------------------------
# PERMUTATION NULL
# ---------------------------------------------------------------------------

def permutation_pvalue(
    train_samples: Sequence[Sample],
    test_samples: Sequence[Sample],
    permutations: int,
    rng: random.Random,
    min_occupancy: int = 1,
    use_inverted: bool = False,
) -> float:
    """
    Within-target permutation null.

    We preserve each target's label multiset but randomly permute labels among
    its modulus cells.

    This is more appropriate than globally shuffling all labels because target
    difficulty / label balance can differ.
    """
    observed_model = StateMajorityModel()
    observed_model.fit(train_samples)

    observed = evaluate_model(
        observed_model,
        test_samples,
        min_occupancy=min_occupancy,
    )

    observed_score = (
        observed.inverted_acc if use_inverted
        else observed.direct_acc
    )

    if math.isnan(observed_score):
        return float("nan")

    # Group train samples by target.
    train_groups: Dict[int, List[Sample]] = defaultdict(list)
    for s in train_samples:
        train_groups[s.target_id].append(s)

    exceed = 0

    for _ in range(permutations):
        shuffled: List[Sample] = []

        for tid, group in train_groups.items():
            labels = [s.y for s in group]
            rng.shuffle(labels)

            for s, y in zip(group, labels):
                shuffled.append(
                    Sample(
                        target_id=s.target_id,
                        ell=s.ell,
                        family=s.family,
                        signature=s.signature,
                        y=y,
                    )
                )

        model = StateMajorityModel()
        model.fit(shuffled)

        score = evaluate_model(
            model,
            test_samples,
            min_occupancy=min_occupancy,
        )

        null_score = (
            score.inverted_acc if use_inverted
            else score.direct_acc
        )

        if not math.isnan(null_score) and null_score >= observed_score:
            exceed += 1

    return (exceed + 1) / (permutations + 1)


# ---------------------------------------------------------------------------
# LOO MODULUS
# ---------------------------------------------------------------------------

def loo_modulus(
    all_samples: Sequence[Sample],
    moduli: Sequence[int],
    rng: random.Random,
    min_occupancy: int = 1,
) -> Dict[int, EvalResult]:
    out: Dict[int, EvalResult] = {}

    for ell in moduli:
        train = [s for s in all_samples if s.ell != ell]
        test = [s for s in all_samples if s.ell == ell]

        model = StateMajorityModel()
        model.fit(train)

        out[ell] = evaluate_model(
            model,
            test,
            min_occupancy=min_occupancy,
        )

    return out


# ---------------------------------------------------------------------------
# 2D HOLDOUT
# ---------------------------------------------------------------------------

def two_dimensional_holdout(
    all_samples: Sequence[Sample],
    heldout_moduli: Sequence[int],
    heldout_targets: Sequence[int],
    min_occupancy: int = 1,
) -> EvalResult:
    heldout_moduli = set(heldout_moduli)
    heldout_targets = set(heldout_targets)

    train = [
        s for s in all_samples
        if s.ell not in heldout_moduli
        and s.target_id not in heldout_targets
    ]

    test = [
        s for s in all_samples
        if s.ell in heldout_moduli
        and s.target_id in heldout_targets
    ]

    model = StateMajorityModel()
    model.fit(train)

    return evaluate_model(
        model,
        test,
        min_occupancy=min_occupancy,
    )


# ---------------------------------------------------------------------------
# REPORTING
# ---------------------------------------------------------------------------

def print_result(label: str, result: EvalResult) -> None:
    print(
        f"{label:<20} "
        f"n={result.total:4d} "
        f"covered={result.covered:4d} "
        f"direct={result.direct_acc:.6f} "
        f"inverted={result.inverted_acc:.6f} "
        f"direct_bal={result.direct_bal:.6f} "
        f"inverted_bal={result.inverted_bal:.6f} "
        f"inv-direct={result.inversion_gain:+.6f}"
    )


def print_state_stats(
    train_samples: Sequence[Sample],
    family: str,
) -> None:
    counts = Counter(s.signature for s in train_samples)

    occupancies = list(counts.values())

    if not occupancies:
        print(f"{family}: no states")
        return

    print(f"{family}")
    print(f"  unique states       = {len(counts)}")
    print(f"  mean occupancy      = {statistics.fmean(occupancies):.3f}")
    print(f"  median occupancy    = {statistics.median(occupancies):.3f}")
    print(f"  max occupancy       = {max(occupancies)}")

    for k in OCCUPANCY_THRESHOLDS:
        covered = sum(v >= k for v in occupancies)
        samples = sum(v for v in occupancies if v >= k)
        print(
            f"  occupancy >= {k}: "
            f"states={covered:4d} "
            f"samples={samples:5d}"
        )


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 77")
    print("ANTI-PREDICTION / SIGN-INVERSION TEST")
    print("DIRECT VS INVERTED N-ONLY CHARACTER PREDICTION")
    print("REUSABLE SIGNATURE STATES")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("WITHIN-TARGET PERMUTATION NULL")
    print("OCCUPANCY-MATCHED DIAGNOSTICS")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # -----------------------------------------------------------------------
    primes = sieve_primes(PRIME_LO, PRIME_HI)

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")

    # -----------------------------------------------------------------------
    # 2. TARGETS
    # -----------------------------------------------------------------------
    targets = generate_targets(primes, NUM_TARGETS, rng)

    print("\n2. TARGET SUMMARY")
    print("-" * 78)

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if NUM_TARGETS > 24:
        print("... remaining generated targets omitted")

    # -----------------------------------------------------------------------
    # 3. VALIDATE IDENTITIES
    # -----------------------------------------------------------------------
    failures = sum(not validate_e1(t) for t in targets)

    print("\n3. E1 IDENTITY VALIDATION")
    print("-" * 78)
    print(f"identity failures = {failures}")

    if failures:
        raise RuntimeError("E1 identity validation failed")

    print("status = PASS")

    # -----------------------------------------------------------------------
    # 4. MODULUS VALIDATION
    # -----------------------------------------------------------------------
    validate_moduli(CYCLOTOMIC)
    validate_moduli(CONTROLS)

    print("\n4. MODULUS VALIDATION")
    print("-" * 78)

    for ell in CYCLOTOMIC:
        r1, r2 = cyclotomic_roots(ell)
        print(
            f"cyclo ell={ell:5d} roots=({r1:5d},{r2:5d}) "
            f"identity_ok=True"
        )

    for ell in CONTROLS:
        r1, r2 = cyclotomic_roots(ell)
        print(
            f"control ell={ell:5d} roots=({r1:5d},{r2:5d}) "
            f"identity_ok=True"
        )

    # -----------------------------------------------------------------------
    # 5. TARGET HOLDOUT
    # -----------------------------------------------------------------------
    train_ids, test_ids = split_targets(targets, rng)

    print("\n5. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_ids)}")
    print(f"test targets     = {len(test_ids)}")

    # -----------------------------------------------------------------------
    # 6. BUILD DATASETS
    # -----------------------------------------------------------------------
    cyclo_samples = build_samples(
        targets,
        CYCLOTOMIC,
        "cyclotomic",
    )

    control_samples = build_samples(
        targets,
        CONTROLS,
        "control",
    )

    print("\n6. DATASET CONSTRUCTION")
    print("-" * 78)
    print(f"cyclotomic samples = {len(cyclo_samples)}")
    print(f"control samples    = {len(control_samples)}")

    # -----------------------------------------------------------------------
    # 7. GLOBAL STATE ANALYSIS
    # -----------------------------------------------------------------------
    print("\n7. TRAINING STATE OCCUPANCY")
    print("-" * 78)

    cyclo_train = by_target(cyclo_samples, train_ids)
    control_train = by_target(control_samples, train_ids)

    print_state_stats(cyclo_train, "CYCLOTOMIC")
    print()
    print_state_stats(control_train, "CONTROL")

    # -----------------------------------------------------------------------
    # 8. TARGET HOLDOUT
    # -----------------------------------------------------------------------
    cyclo_test = by_target(cyclo_samples, test_ids)
    control_test = by_target(control_samples, test_ids)

    cyclo_model = StateMajorityModel()
    cyclo_model.fit(cyclo_train)

    control_model = StateMajorityModel()
    control_model.fit(control_train)

    print("\n8. TARGET HOLDOUT")
    print("-" * 78)

    c_result = evaluate_model(
        cyclo_model,
        cyclo_test,
        min_occupancy=1,
    )

    r_result = evaluate_model(
        control_model,
        control_test,
        min_occupancy=1,
    )

    print_result("cyclotomic", c_result)
    print_result("control", r_result)

    print(
        f"target delta direct       = "
        f"{c_result.direct_acc - r_result.direct_acc:+.6f}"
    )
    print(
        f"target delta inverted     = "
        f"{c_result.inverted_acc - r_result.inverted_acc:+.6f}"
    )
    print(
        f"cyclotomic inversion gain = "
        f"{c_result.inversion_gain:+.6f}"
    )

    # -----------------------------------------------------------------------
    # 9. OCCUPANCY THRESHOLDS
    # -----------------------------------------------------------------------
    print("\n9. OCCUPANCY-MATCHED TARGET HOLDOUT")
    print("-" * 78)

    for k in OCCUPANCY_THRESHOLDS:
        print(f"\nminimum training occupancy = {k}")

        cr = evaluate_model(
            cyclo_model,
            cyclo_test,
            min_occupancy=k,
        )

        rr = evaluate_model(
            control_model,
            control_test,
            min_occupancy=k,
        )

        print_result("cyclotomic", cr)
        print_result("control", rr)

    # -----------------------------------------------------------------------
    # 10. LEAVE-ONE-MODULUS-OUT
    # -----------------------------------------------------------------------
    print("\n10. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    cyclo_loo = loo_modulus(
        cyclo_samples,
        CYCLOTOMIC,
        rng,
        min_occupancy=1,
    )

    control_loo = loo_modulus(
        control_samples,
        CONTROLS,
        rng,
        min_occupancy=1,
    )

    print("CYCLOTOMIC")
    for ell in CYCLOTOMIC:
        res = cyclo_loo[ell]
        print(
            f"ell={ell:5d} "
            f"covered={res.covered:4d}/{res.total:4d} "
            f"direct={res.direct_acc:.6f} "
            f"inverted={res.inverted_acc:.6f} "
            f"inv-direct={res.inversion_gain:+.6f}"
        )

    print(
        f"cyclotomic LOO direct mean   = "
        f"{mean_or_nan(r.direct_acc for r in cyclo_loo.values()):.6f}"
    )
    print(
        f"cyclotomic LOO inverted mean = "
        f"{mean_or_nan(r.inverted_acc for r in cyclo_loo.values()):.6f}"
    )

    print("\nCONTROL")
    for ell in CONTROLS:
        res = control_loo[ell]
        print(
            f"ell={ell:5d} "
            f"covered={res.covered:4d}/{res.total:4d} "
            f"direct={res.direct_acc:.6f} "
            f"inverted={res.inverted_acc:.6f} "
            f"inv-direct={res.inversion_gain:+.6f}"
        )

    print(
        f"control LOO direct mean      = "
        f"{mean_or_nan(r.direct_acc for r in control_loo.values()):.6f}"
    )
    print(
        f"control LOO inverted mean    = "
        f"{mean_or_nan(r.inverted_acc for r in control_loo.values()):.6f}"
    )

    # -----------------------------------------------------------------------
    # 11. TWO-DIMENSIONAL HOLDOUT
    # -----------------------------------------------------------------------
    print("\n11. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)

    heldout_targets = test_ids

    c2d = two_dimensional_holdout(
        cyclo_samples,
        CYCLO_2D_HELDOUT,
        heldout_targets,
        min_occupancy=1,
    )

    r2d = two_dimensional_holdout(
        control_samples,
        CONTROL_2D_HELDOUT,
        heldout_targets,
        min_occupancy=1,
    )

    print_result("cyclotomic", c2d)
    print_result("control", r2d)

    # -----------------------------------------------------------------------
    # 12. WITHIN-TARGET PERMUTATION NULL
    # -----------------------------------------------------------------------
    print("\n12. WITHIN-TARGET PERMUTATION NULL")
    print("-" * 78)

    # Keep permutations moderate so this remains lightweight.
    c_p_direct = permutation_pvalue(
        cyclo_train,
        cyclo_test,
        PERMUTATIONS,
        rng,
        min_occupancy=1,
        use_inverted=False,
    )

    c_p_inverted = permutation_pvalue(
        cyclo_train,
        cyclo_test,
        PERMUTATIONS,
        rng,
        min_occupancy=1,
        use_inverted=True,
    )

    r_p_direct = permutation_pvalue(
        control_train,
        control_test,
        PERMUTATIONS,
        rng,
        min_occupancy=1,
        use_inverted=False,
    )

    r_p_inverted = permutation_pvalue(
        control_train,
        control_test,
        PERMUTATIONS,
        rng,
        min_occupancy=1,
        use_inverted=True,
    )

    print(f"cyclotomic direct   p = {c_p_direct:.6f}")
    print(f"cyclotomic inverted p = {c_p_inverted:.6f}")
    print(f"control direct      p = {r_p_direct:.6f}")
    print(f"control inverted    p = {r_p_inverted:.6f}")

    # -----------------------------------------------------------------------
    # 13. CROSS-MODULUS INVERSION SUMMARY
    # -----------------------------------------------------------------------
    print("\n13. SIGN-INVERSION SUMMARY")
    print("-" * 78)

    c_inv_oos = c_result.inverted_acc
    c_dir_oos = c_result.direct_acc
    r_inv_oos = r_result.inverted_acc
    r_dir_oos = r_result.direct_acc

    print(f"cyclotomic direct OOS   = {c_dir_oos:.6f}")
    print(f"cyclotomic inverted OOS = {c_inv_oos:.6f}")
    print(f"control direct OOS      = {r_dir_oos:.6f}")
    print(f"control inverted OOS    = {r_inv_oos:.6f}")

    print()
    if (
        not math.isnan(c_inv_oos)
        and c_inv_oos > c_dir_oos
        and c_inv_oos > 0.5
        and c_inv_oos > r_inv_oos
    ):
        print("STATUS = INVERSION SIGNAL WORTH FOLLOW-UP")
    else:
        print("STATUS = NO TRANSFERABLE INVERSION SIGNAL")

    # -----------------------------------------------------------------------
    # 14. FINAL DIAGNOSTIC
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("14. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        f"cyclotomic target direct   = {c_result.direct_acc:.6f}"
    )
    print(
        f"cyclotomic target inverted = {c_result.inverted_acc:.6f}"
    )
    print(
        f"control target direct      = {r_result.direct_acc:.6f}"
    )
    print(
        f"control target inverted    = {r_result.inverted_acc:.6f}"
    )

    print(
        f"cyclotomic LOO direct      = "
        f"{mean_or_nan(r.direct_acc for r in cyclo_loo.values()):.6f}"
    )
    print(
        f"cyclotomic LOO inverted    = "
        f"{mean_or_nan(r.inverted_acc for r in cyclo_loo.values()):.6f}"
    )
    print(
        f"control LOO direct         = "
        f"{mean_or_nan(r.direct_acc for r in control_loo.values()):.6f}"
    )
    print(
        f"control LOO inverted       = "
        f"{mean_or_nan(r.inverted_acc for r in control_loo.values()):.6f}"
    )

    print(
        f"cyclotomic 2D direct      = {c2d.direct_acc:.6f}"
    )
    print(
        f"cyclotomic 2D inverted    = {c2d.inverted_acc:.6f}"
    )
    print(
        f"control 2D direct         = {r2d.direct_acc:.6f}"
    )
    print(
        f"control 2D inverted       = {r2d.inverted_acc:.6f}"
    )

    print()
    print("Interpretation")
    print("--------------")
    print(
        "A genuine sign-inversion effect requires the inverted predictor "
        "to outperform the direct predictor on unseen targets AND unseen "
        "moduli, while the same inversion does not appear in controls."
    )
    print(
        "An inverted score below 50% is not evidence by itself; its complement "
        "must also survive the same out-of-sample tests."
    )
    print(
        "This experiment is diagnostic only. It does not establish a "
        "factorization algorithm or an E1-from-N identity."
    )

    print("\n" + "=" * 78)
    print("EXPERIMENT 77 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

