#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 66
TRUE LEGENDRE RELATIVE ORIENTATION
N-ONLY PREDICTION OF CHI_ell(d1) * CHI_ref(d1_ref)
TARGET + MODULUS + 2D HOLDOUT
LEAVE-ONE-MODULUS-OUT
GLOBAL + WITHIN-TARGET PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

IMPORTANT
---------
Experiment 65 used

    (d1-d2)/(zeta-zeta2)

which is identically +1.

This experiment fixes that by using the actual quadratic-character
orientation:

    sigma_ell = chi_ell(d1)

only when

    chi_ell(d1) == chi_ell(d2) != 0.

For a fixed reference modulus ell0:

    rho_ell = sigma_ell * sigma_ell0

where rho is genuinely +1 or -1.

The label is therefore NOT algebraically forced.

No p, q, or s is used in the feature vector.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 6601

NUM_TARGETS = 120
TRAIN_TARGETS = 90

PERMUTATIONS = 500

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

RAW_CONTROLS = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

CYCLO_REFERENCE = 7
CONTROL_REFERENCE = 673


# ============================================================================
# TARGET DATA
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


def generate_prime_population(
    lo: int = 2_000_000,
    hi: int = 4_200_000,
) -> list[int]:
    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(hi)) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [i for i in range(lo, hi + 1) if sieve[i]]


def generate_targets(
    primes: list[int],
    count: int,
    seed: int,
) -> list[Target]:
    rng = random.Random(seed)
    used: set[tuple[int, int]] = set()
    out: list[Target] = []

    while len(out) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pair = (p, q)
        if pair in used:
            continue

        used.add(pair)
        out.append(Target(len(out) + 1, p, q))

    return out


# ============================================================================
# MODULAR ALGEBRA
# ============================================================================

def roots_x2_x1(ell: int) -> tuple[int, int] | None:
    """
    Roots of x^2+x+1 = 0 mod ell.

    For prime ell, these exist exactly when ell == 1 mod 3.
    """
    if ell == 2 or ell % 3 != 1:
        return None

    roots = [
        x for x in range(1, ell)
        if (x * x + x + 1) % ell == 0
    ]

    if len(roots) != 2:
        raise RuntimeError(
            f"{ell}: expected two roots, found {roots}"
        )

    return tuple(sorted(roots))


def legendre(a: int, ell: int) -> int:
    a %= ell

    if a == 0:
        return 0

    y = pow(a, (ell - 1) // 2, ell)

    if y == 1:
        return 1

    if y == ell - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre value ell={ell}, a={a}, y={y}"
    )


@dataclass(frozen=True)
class LocalState:
    ell: int
    zeta: int
    zeta2: int
    d1: int
    d2: int
    chi1: int
    chi2: int
    same_sign: bool
    sigma: int | None


def local_state(n: int, ell: int) -> LocalState:
    roots = roots_x2_x1(ell)

    if roots is None:
        raise RuntimeError(
            f"ell={ell} has no x^2+x+1 roots"
        )

    zeta, zeta2 = roots

    d1 = (zeta - 4 * n) % ell
    d2 = (zeta2 - 4 * n) % ell

    chi1 = legendre(d1, ell)
    chi2 = legendre(d2, ell)

    same = (
        chi1 != 0
        and chi2 != 0
        and chi1 == chi2
    )

    sigma = chi1 if same else None

    return LocalState(
        ell=ell,
        zeta=zeta,
        zeta2=zeta2,
        d1=d1,
        d2=d2,
        chi1=chi1,
        chi2=chi2,
        same_sign=same,
        sigma=sigma,
    )


# ============================================================================
# N-ONLY FEATURES
# ============================================================================

def cubic_class(x: int, ell: int) -> int:
    x %= ell

    if x == 0:
        return 0

    y = pow(x, (ell - 1) // 3, ell)

    if y == 1:
        return 0

    roots = roots_x2_x1(ell)
    assert roots is not None

    zeta, zeta2 = roots

    if y == zeta:
        return 1

    if y == zeta2:
        return 2

    raise RuntimeError(
        f"Unexpected cubic class ell={ell}, x={x}, y={y}"
    )


def n_features(n: int, ell: int) -> list[float]:
    """
    STRICTLY N-ONLY.

    No p, q, s, d1, d2, or sigma enters here.
    """
    G = 16 * n * n + 4 * n + 1

    values = [
        n,
        n + 1,
        n - 1,
        4 * n + 1,
        G,
    ]

    out: list[float] = []

    for x in values:
        out.append(float(legendre(x, ell)))

    for x in [n, n + 1, 4 * n + 1]:
        c = cubic_class(x, ell)

        out.extend([
            1.0 if c == 0 else 0.0,
            1.0 if c == 1 else 0.0,
            1.0 if c == 2 else 0.0,
        ])

    return out


# ============================================================================
# DATASET
# ============================================================================

@dataclass
class Sample:
    target_id: int
    ell: int
    x: list[float]
    y: int


def build_relative_dataset(
    targets: list[Target],
    moduli: list[int],
    reference_ell: int,
) -> tuple[
    list[Sample],
    int,
    int,
    int,
]:
    """
    Label:

        rho = sigma_ell * sigma_reference

    only if both local cells are same-sign and nonzero.
    """
    samples: list[Sample] = []

    total_cells = 0
    degenerate = 0

    for t in targets:
        ref = local_state(t.n, reference_ell)

        for ell in moduli:
            if ell == reference_ell:
                continue

            total_cells += 1

            cur = local_state(t.n, ell)

            if not ref.same_sign or not cur.same_sign:
                degenerate += 1
                continue

            assert ref.sigma in (-1, 1)
            assert cur.sigma in (-1, 1)

            rho = ref.sigma * cur.sigma

            if rho not in (-1, 1):
                raise RuntimeError("Invalid relative orientation")

            samples.append(
                Sample(
                    target_id=t.idx,
                    ell=ell,
                    x=n_features(t.n, ell),
                    y=rho,
                )
            )

    positive = sum(s.y == 1 for s in samples)
    negative = sum(s.y == -1 for s in samples)

    return samples, positive, negative, total_cells - len(samples)


# ============================================================================
# PURE PYTHON RIDGE CLASSIFIER
# ============================================================================

def transpose(A: list[list[float]]) -> list[list[float]]:
    return [list(c) for c in zip(*A)]


def matmul(
    A: list[list[float]],
    B: list[list[float]],
) -> list[list[float]]:
    rows = len(A)
    inner = len(B)
    cols = len(B[0])

    C = [[0.0] * cols for _ in range(rows)]

    for i in range(rows):
        for k in range(inner):
            if A[i][k] == 0:
                continue

            for j in range(cols):
                C[i][j] += A[i][k] * B[k][j]

    return C


def solve(A: list[list[float]], b: list[float]) -> list[float]:
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for col in range(n):
        pivot = max(
            range(col, n),
            key=lambda r: abs(M[r][col]),
        )

        if abs(M[pivot][col]) < 1e-14:
            raise RuntimeError("Singular system")

        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]

        pv = M[col][col]

        for j in range(col, n + 1):
            M[col][j] /= pv

        for r in range(n):
            if r == col:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, n + 1):
                M[r][j] -= factor * M[col][j]

    return [M[i][n] for i in range(n)]


class RidgeBinary:
    def __init__(self, alpha: float = 2.0):
        self.alpha = alpha
        self.mean: list[float] | None = None
        self.scale: list[float] | None = None
        self.w: list[float] | None = None

    def fit(self, X: list[list[float]], y: list[int]) -> None:
        if not X:
            raise ValueError("X is empty")

        if len(X) != len(y):
            raise ValueError(
                f"X/y length mismatch: len(X)={len(X)}, len(y)={len(y)}"
            )

        d = len(X[0])

        if any(len(row) != d for row in X):
            raise ValueError("Inconsistent feature dimensions in X")

        # Standardize features.
        self.mean = [
            statistics.fmean(row[j] for row in X)
            for j in range(d)
        ]

        self.scale = []

        for j in range(d):
            vals = [row[j] for row in X]
            m = self.mean[j]

            var = statistics.fmean(
                (v - m) ** 2
                for v in vals
            )

            s = math.sqrt(var)

            self.scale.append(
                1.0 if s < 1e-12 else s
            )

        Xs = [
            [
                (row[j] - self.mean[j]) / self.scale[j]
                for j in range(d)
            ]
            for row in X
        ]

        # Add intercept column.
        Z = [
            [1.0] + row
            for row in Xs
        ]

        Zt = transpose(Z)

        # Z^T Z
        ZtZ = matmul(Zt, Z)

        # Ridge penalty on feature coefficients,
        # but NOT on the intercept.
        for i in range(1, len(ZtZ)):
            ZtZ[i][i] += self.alpha

        # Correct calculation of Z^T y.
        #
        # Z has shape:
        #     samples x (features + 1)
        #
        # Therefore:
        #     (Z^T y)[i] = sum_r Z[r][i] * y[r]
        Zty = [
            sum(
                Z[r][i] * float(y[r])
                for r in range(len(y))
            )
            for i in range(len(Zt))
        ]

        self.w = solve(ZtZ, Zty)

    def predict(self, X: list[list[float]]) -> list[int]:
        if self.mean is None:
            raise RuntimeError("Model has not been fitted")

        if self.scale is None:
            raise RuntimeError("Model has not been fitted")

        if self.w is None:
            raise RuntimeError("Model has not been fitted")

        out = []

        for row in X:
            if len(row) != len(self.mean):
                raise ValueError(
                    "Prediction feature dimension does not match training"
                )

            z = [1.0]

            for j in range(len(row)):
                z.append(
                    (row[j] - self.mean[j]) / self.scale[j]
                )

            score = sum(
                a * b
                for a, b in zip(self.w, z)
            )

            out.append(
                1 if score >= 0 else -1
            )

        return out


# ============================================================================
# EVALUATION
# ============================================================================

def split_targets(
    samples: list[Sample],
    train_ids: set[int],
) -> tuple[list[Sample], list[Sample]]:
    train = [
        s for s in samples
        if s.target_id in train_ids
    ]

    test = [
        s for s in samples
        if s.target_id not in train_ids
    ]

    return train, test


def score_model(
    train: list[Sample],
    test: list[Sample],
) -> tuple[float, float]:
    if not train or not test:
        return float("nan"), float("nan")

    labels = set(s.y for s in train)

    if len(labels) < 2:
        majority = max(
            sum(s.y == 1 for s in test),
            sum(s.y == -1 for s in test),
        ) / len(test)

        return 1.0, majority

    model = RidgeBinary(alpha=2.0)

    model.fit(
        [s.x for s in train],
        [s.y for s in train],
    )

    train_pred = model.predict(
        [s.x for s in train]
    )

    test_pred = model.predict(
        [s.x for s in test]
    )

    train_acc = sum(
        p == s.y
        for p, s in zip(train_pred, train)
    ) / len(train)

    test_acc = sum(
        p == s.y
        for p, s in zip(test_pred, test)
    ) / len(test)

    return train_acc, test_acc


def majority_baseline(samples: list[Sample]) -> float:
    if not samples:
        return float("nan")

    p = sum(s.y == 1 for s in samples)
    m = len(samples) - p

    return max(p, m) / len(samples)


# ============================================================================
# PERMUTATION NULLS
# ============================================================================

def permutation_test(
    train: list[Sample],
    test: list[Sample],
    observed: float,
    permutations: int,
    seed: int,
    within_target: bool = False,
) -> float:
    rng = random.Random(seed)

    groups: dict[int, list[Sample]] = {}

    if within_target:
        for s in train:
            groups.setdefault(
                s.target_id, []
            ).append(s)

    exceed = 0

    for _ in range(permutations):
        if within_target:
            shuffled_y: dict[int, list[int]] = {}

            for tid, group in groups.items():
                ys = [s.y for s in group]
                rng.shuffle(ys)
                shuffled_y[tid] = ys

            idxs = {
                tid: 0
                for tid in groups
            }

            perm_train = []

            for s in train:
                y = shuffled_y[s.target_id][idxs[s.target_id]]
                idxs[s.target_id] += 1

                perm_train.append(
                    Sample(
                        target_id=s.target_id,
                        ell=s.ell,
                        x=s.x,
                        y=y,
                    )
                )

        else:
            ys = [s.y for s in train]
            rng.shuffle(ys)

            perm_train = [
                Sample(
                    target_id=s.target_id,
                    ell=s.ell,
                    x=s.x,
                    y=y,
                )
                for s, y in zip(train, ys)
            ]

        _, score = score_model(
            perm_train,
            test,
        )

        if score >= observed:
            exceed += 1

    return (exceed + 1) / (permutations + 1)


# ============================================================================
# MODULUS HOLDOUT
# ============================================================================

def leave_one_modulus_out(
    samples: list[Sample],
    train_ids: set[int],
    moduli: list[int],
) -> list[tuple[int, float, int, int]]:
    results = []

    for ell in moduli:
        train = [
            s for s in samples
            if s.target_id in train_ids
            and s.ell != ell
        ]

        test = [
            s for s in samples
            if s.target_id in train_ids
            and s.ell == ell
        ]

        if (
            len(train) == 0
            or len(test) == 0
            or len(set(s.y for s in train)) < 2
        ):
            results.append(
                (ell, float("nan"), len(train), len(test))
            )
            continue

        model = RidgeBinary(alpha=2.0)

        model.fit(
            [s.x for s in train],
            [s.y for s in train],
        )

        pred = model.predict(
            [s.x for s in test]
        )

        acc = sum(
            p == s.y
            for p, s in zip(pred, test)
        ) / len(test)

        results.append(
            (ell, acc, len(train), len(test))
        )

    return results


# ============================================================================
# TWO DIMENSIONAL HOLDOUT
# ============================================================================

def two_dimensional_holdout(
    samples: list[Sample],
    train_ids: set[int],
    heldout_moduli: set[int],
) -> float:
    train = [
        s for s in samples
        if s.target_id in train_ids
        and s.ell not in heldout_moduli
    ]

    test = [
        s for s in samples
        if s.target_id not in train_ids
        and s.ell in heldout_moduli
    ]

    if not train or not test:
        return float("nan")

    if len(set(s.y for s in train)) < 2:
        return majority_baseline(test)

    model = RidgeBinary(alpha=2.0)

    model.fit(
        [s.x for s in train],
        [s.y for s in train],
    )

    pred = model.predict(
        [s.x for s in test]
    )

    return sum(
        p == s.y
        for p, s in zip(pred, test)
    ) / len(test)


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 66")
    print("TRUE LEGENDRE RELATIVE ORIENTATION")
    print("N-ONLY PREDICTION")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("GLOBAL + WITHIN-TARGET PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. Population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = generate_prime_population()

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(
        f"generation time  = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    targets = generate_targets(
        primes,
        NUM_TARGETS,
        SEED,
    )

    print("\n2. TARGETS")
    print("-" * 78)
    print(f"total targets = {len(targets)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    print("... remaining targets omitted")

    # ------------------------------------------------------------------
    # 3. Moduli
    # ------------------------------------------------------------------

    valid_controls = [
        e for e in RAW_CONTROLS
        if roots_x2_x1(e) is not None
    ]

    excluded_controls = [
        e for e in RAW_CONTROLS
        if roots_x2_x1(e) is None
    ]

    print("\n3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic       = {CYCLOTOMIC}")
    print(f"valid controls   = {valid_controls}")
    print(f"excluded controls= {excluded_controls}")

    # ------------------------------------------------------------------
    # 4. Reference validation
    # ------------------------------------------------------------------

    for ref in [CYCLO_REFERENCE, CONTROL_REFERENCE]:
        roots = roots_x2_x1(ref)

        if roots is None:
            raise RuntimeError(
                f"Reference ell={ref} invalid"
            )

    print("\n4. REFERENCE VALIDATION")
    print("-" * 78)
    print(
        f"cyclotomic reference = {CYCLO_REFERENCE} "
        f"roots={roots_x2_x1(CYCLO_REFERENCE)}"
    )
    print(
        f"control reference    = {CONTROL_REFERENCE} "
        f"roots={roots_x2_x1(CONTROL_REFERENCE)}"
    )

    # ------------------------------------------------------------------
    # 5. Direct label validation
    # ------------------------------------------------------------------

    direct_failures = 0

    for t in targets:
        for ell in CYCLOTOMIC:
            st = local_state(t.n, ell)

            if st.same_sign:
                if st.sigma not in (-1, 1):
                    direct_failures += 1

    print("\n5. LOCAL LEGENDRE ORIENTATION VALIDATION")
    print("-" * 78)
    print(
        "same-sign sigma values must be only +/-1"
    )
    print(f"validation failures = {direct_failures}")

    if direct_failures:
        raise RuntimeError(
            "Local Legendre orientation validation failed"
        )

    # ------------------------------------------------------------------
    # 6. Dataset
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    cyclo_samples, cpos, cneg, cdeg = build_relative_dataset(
        targets,
        CYCLOTOMIC,
        CYCLO_REFERENCE,
    )

    control_samples, rpos, rneg, rdeg = build_relative_dataset(
        targets,
        valid_controls,
        CONTROL_REFERENCE,
    )

    print("\n6. RELATIVE-ORIENTATION DATASET")
    print("-" * 78)

    print(
        f"cyclotomic samples = {len(cyclo_samples)} "
        f"positive={cpos} negative={cneg} "
        f"degenerate={cdeg}"
    )

    print(
        f"control samples    = {len(control_samples)} "
        f"positive={rpos} negative={rneg} "
        f"degenerate={rdeg}"
    )

    print(
        f"cyclotomic majority baseline = "
        f"{majority_baseline(cyclo_samples):.6f}"
    )

    print(
        f"control majority baseline = "
        f"{majority_baseline(control_samples):.6f}"
    )

    print(
        f"dataset construction time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # Immediate sanity guard.
    if cpos == 0 or cneg == 0:
        raise RuntimeError(
            "Cyclotomic relative label is constant; "
            "this experiment is trivial."
        )

    if rpos == 0 or rneg == 0:
        print(
            "WARNING: control relative label is constant."
        )

    # ------------------------------------------------------------------
    # 7. Target split
    # ------------------------------------------------------------------

    ids = list(range(1, NUM_TARGETS + 1))

    rng = random.Random(SEED + 1)
    rng.shuffle(ids)

    train_ids = set(ids[:TRAIN_TARGETS])
    test_ids = set(ids[TRAIN_TARGETS:])

    print("\n7. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_ids)}")
    print(f"test targets     = {len(test_ids)}")
    print(f"test ids         = {sorted(test_ids)}")

    ctrain, ctest = split_targets(
        cyclo_samples,
        train_ids,
    )

    rtrain, rtest = split_targets(
        control_samples,
        train_ids,
    )

    # ------------------------------------------------------------------
    # 8. Target holdout
    # ------------------------------------------------------------------

    ctrain_acc, ctest_acc = score_model(
        ctrain,
        ctest,
    )

    rtrain_acc, rtest_acc = score_model(
        rtrain,
        rtest,
    )

    print("\n8. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"cyclotomic train={ctrain_acc:.6f} "
        f"test={ctest_acc:.6f}"
    )
    print(
        f"control    train={rtrain_acc:.6f} "
        f"test={rtest_acc:.6f}"
    )
    print(
        f"test delta = "
        f"{ctest_acc-rtest_acc:+.6f}"
    )

    # ------------------------------------------------------------------
    # 9. Repeated target holdout
    # ------------------------------------------------------------------

    print("\n9. REPEATED TARGET HOLDOUT")
    print("-" * 78)

    c_scores = []
    r_scores = []

    for rep in range(10):
        rrng = random.Random(SEED + 100 + rep)
        rid = list(range(1, NUM_TARGETS + 1))
        rrng.shuffle(rid)

        tr_ids = set(rid[:TRAIN_TARGETS])

        ct, cv = split_targets(
            cyclo_samples,
            tr_ids,
        )

        rt, rv = split_targets(
            control_samples,
            tr_ids,
        )

        _, ca = score_model(ct, cv)
        _, ra = score_model(rt, rv)

        c_scores.append(ca)
        r_scores.append(ra)

    print(
        f"cyclotomic mean   = "
        f"{statistics.fmean(c_scores):.6f}"
    )
    print(
        f"cyclotomic median = "
        f"{statistics.median(c_scores):.6f}"
    )
    print(
        f"cyclotomic min    = "
        f"{min(c_scores):.6f}"
    )
    print(
        f"cyclotomic max    = "
        f"{max(c_scores):.6f}"
    )

    print(
        f"control mean      = "
        f"{statistics.fmean(r_scores):.6f}"
    )
    print(
        f"control median    = "
        f"{statistics.median(r_scores):.6f}"
    )

    # ------------------------------------------------------------------
    # 10. Leave-one-modulus-out
    # ------------------------------------------------------------------

    print("\n10. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    c_moduli = [
        e for e in CYCLOTOMIC
        if e != CYCLO_REFERENCE
    ]

    r_moduli = [
        e for e in valid_controls
        if e != CONTROL_REFERENCE
    ]

    c_loo = leave_one_modulus_out(
        cyclo_samples,
        train_ids,
        c_moduli,
    )

    r_loo = leave_one_modulus_out(
        control_samples,
        train_ids,
        r_moduli,
    )

    for ell, acc, ntrain, ntest in c_loo:
        print(
            f"cyclo ell={ell:5d} "
            f"accuracy={acc if not math.isnan(acc) else float('nan'):.6f} "
            f"train={ntrain} test={ntest}"
        )

    cvals = [
        acc for _, acc, _, _ in c_loo
        if not math.isnan(acc)
    ]

    if cvals:
        print(
            f"cyclotomic LOO mean = "
            f"{statistics.fmean(cvals):.6f}"
        )
        print(
            f"cyclotomic LOO median = "
            f"{statistics.median(cvals):.6f}"
        )
    else:
        print(
            "cyclotomic LOO mean = nan "
            "(no valid LOO evaluations)"
        )

    for ell, acc, ntrain, ntest in r_loo:
        print(
            f"control ell={ell:5d} "
            f"accuracy={acc if not math.isnan(acc) else float('nan'):.6f} "
            f"train={ntrain} test={ntest}"
        )

    rvals = [
        acc for _, acc, _, _ in r_loo
        if not math.isnan(acc)
    ]

    if rvals:
        print(
            f"control LOO mean = "
            f"{statistics.fmean(rvals):.6f}"
        )
        print(
            f"control LOO median = "
            f"{statistics.median(rvals):.6f}"
        )
    else:
        print(
            "control LOO mean = nan "
            "(no valid LOO evaluations)"
        )

    # ------------------------------------------------------------------
    # 11. Two-dimensional holdout
    # ------------------------------------------------------------------

    held_cyclo = {61, 127, 307}
    held_control = {6211, 7951}

    c2d = two_dimensional_holdout(
        cyclo_samples,
        train_ids,
        held_cyclo,
    )

    r2d = two_dimensional_holdout(
        control_samples,
        train_ids,
        held_control,
    )

    print("\n11. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)
    print(
        f"cyclotomic 2D = {c2d:.6f} "
        f"heldout={sorted(held_cyclo)}"
    )
    print(
        f"control 2D    = {r2d:.6f} "
        f"heldout={sorted(held_control)}"
    )
    print(
        f"2D delta = {c2d-r2d:+.6f}"
    )

    # ------------------------------------------------------------------
    # 12. Baselines
    # ------------------------------------------------------------------

    print("\n12. SIMPLE N-ONLY BASELINES")
    print("-" * 78)

    for name, idx in [
        ("chi(n)", 0),
        ("chi(n+1)", 1),
        ("chi(n-1)", 2),
        ("chi(4n+1)", 3),
        ("chi(G)", 4),
    ]:
        usable_c = [
            s for s in ctest
            if s.x[idx] != 0
        ]

        usable_r = [
            s for s in rtest
            if s.x[idx] != 0
        ]

        if usable_c:
            cpred = [
                1 if s.x[idx] > 0 else -1
                for s in usable_c
            ]
            ca = sum(
                p == s.y
                for p, s in zip(cpred, usable_c)
            ) / len(usable_c)
        else:
            ca = float("nan")

        if usable_r:
            rpred = [
                1 if s.x[idx] > 0 else -1
                for s in usable_r
            ]
            ra = sum(
                p == s.y
                for p, s in zip(rpred, usable_r)
            ) / len(usable_r)
        else:
            ra = float("nan")

        print(
            f"{name:12s} "
            f"C={ca:.6f} R={ra:.6f} "
            f"delta={ca-ra:+.6f}"
        )

    # ------------------------------------------------------------------
    # 13. Permutation tests
    # ------------------------------------------------------------------

    print("\n13. PERMUTATION NULL")
    print("-" * 78)

    global_p = permutation_test(
        ctrain,
        ctest,
        ctest_acc,
        PERMUTATIONS,
        SEED + 5000,
        within_target=False,
    )

    within_target_p = permutation_test(
        ctrain,
        ctest,
        ctest_acc,
        PERMUTATIONS,
        SEED + 6000,
        within_target=True,
    )

    print(
        f"observed cyclotomic accuracy = "
        f"{ctest_acc:.6f}"
    )
    print(
        f"global-label permutation p = "
        f"{global_p:.6f}"
    )
    print(
        f"within-target permutation p = "
        f"{within_target_p:.6f}"
    )

    # ------------------------------------------------------------------
    # 14. Relative-label balance by modulus
    # ------------------------------------------------------------------

    print("\n14. RELATIVE-LABEL BALANCE BY MODULUS")
    print("-" * 78)

    for ell in c_moduli:
        rows = [
            s for s in cyclo_samples
            if s.ell == ell
        ]

        pos = sum(s.y == 1 for s in rows)
        neg = sum(s.y == -1 for s in rows)

        total = len(rows)

        print(
            f"ell={ell:5d} "
            f"n={total:3d} "
            f"+={pos:3d} "
            f"-={neg:3d} "
            f"positive_fraction="
            f"{pos/total if total else float('nan'):.6f}"
        )

    # ------------------------------------------------------------------
    # 15. Final diagnostic
    # ------------------------------------------------------------------

    c_loo_mean = (
        statistics.fmean(cvals)
        if cvals else float("nan")
    )

    r_loo_mean = (
        statistics.fmean(rvals)
        if rvals else float("nan")
    )

    print("\n15. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        f"cyclotomic target holdout = {ctest_acc:.6f}"
    )
    print(
        f"control target holdout    = {rtest_acc:.6f}"
    )
    print(
        f"target delta              = "
        f"{ctest_acc-rtest_acc:+.6f}"
    )

    print(
        f"cyclotomic repeated mean  = "
        f"{statistics.fmean(c_scores):.6f}"
    )
    print(
        f"control repeated mean     = "
        f"{statistics.fmean(r_scores):.6f}"
    )

    print(
        f"cyclotomic LOO mean       = "
        f"{c_loo_mean:.6f}"
    )
    print(
        f"control LOO mean          = "
        f"{r_loo_mean:.6f}"
    )

    print(
        f"cyclotomic 2D             = "
        f"{c2d:.6f}"
    )
    print(
        f"control 2D                = "
        f"{r2d:.6f}"
    )

    print(
        f"global permutation p      = "
        f"{global_p:.6f}"
    )

    print(
        f"within-target permutation p = "
        f"{within_target_p:.6f}"
    )

    print("\nINTERPRETATION")
    print("-" * 78)

    if (
        ctest_acc > 0.55
        and ctest_acc > rtest_acc + 0.05
        and (
            math.isnan(c_loo_mean)
            or math.isnan(r_loo_mean)
            or c_loo_mean > r_loo_mean + 0.03
        )
        and c2d > 0.55
        and within_target_p < 0.05
    ):
        print(
            "POTENTIAL RELATIVE-ORIENTATION SIGNAL"
        )
        print(
            "The nontrivial Legendre relative orientation "
            "appears to generalize beyond the training targets."
        )
    elif (
        ctest_acc > 0.55
        and ctest_acc > rtest_acc + 0.05
    ):
        print(
            "TARGET-LEVEL EFFECT ONLY"
        )
        print(
            "There is some target-holdout separation, "
            "but cross-modulus/generalization evidence is insufficient."
        )
    else:
        print(
            "NO ROBUST N-ONLY RELATIVE-ORIENTATION SIGNAL"
        )
        print(
            "The observed accuracy is not sufficiently different "
            "from the controls/null."
        )

    print("\nCRITICAL DISTINCTION")
    print("-" * 78)
    print(
        "Unlike Experiment 65, this label is NOT algebraically "
        "identically +1."
    )
    print(
        "It is the actual relative Legendre branch orientation."
    )
    print(
        "Therefore a nontrivial result here would be meaningful."
    )

    print("\n" + "=" * 78)
    print("EXPERIMENT 66 COMPLETE")
    print(
        f"total runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

