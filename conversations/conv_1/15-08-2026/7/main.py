#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 65
REFERENCE-ANCHORED RELATIVE ORIENTATION
N-ONLY PREDICTION OF EPSILON_ELL * EPSILON_REFERENCE
TARGET + MODULUS + 2D HOLDOUT
LEAVE-ONE-MODULUS-OUT
PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

CORE IDEA
---------
For each valid modulus ell choose an ordered root pair:

    zeta, zeta^2

and define

    d1 = zeta^2 - 4n
    d2 = zeta^2_conjugate - 4n

with

    epsilon_ell = (d1 - d2) / omega,   omega = zeta - zeta^2

so epsilon_ell is expected to be +/-1 whenever both branches are
non-degenerate.

Instead of predicting epsilon_ell directly, anchor everything to a
fixed reference modulus ell0:

    rho_ell = epsilon_ell * epsilon_ell0

This removes the arbitrary global sign ambiguity.

The experiment asks:

    Can N-only information predict rho_ell for unseen targets
    and unseen moduli?

A useful result requires:
    * balanced nontrivial labels
    * target-holdout > baseline
    * leave-one-modulus-out > baseline
    * 2D holdout > baseline
    * cyclotomic advantage over controls
    * small permutation p-value

If rho is constant, the experiment automatically reports that and
does NOT pretend that 100% accuracy is meaningful.
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 6501

NUM_TARGETS = 120
TRAIN_TARGETS = 90
TEST_TARGETS = NUM_TARGETS - TRAIN_TARGETS

PERMUTATIONS = 500

# Fixed reference moduli.
CYCLO_REFERENCE = 7
CONTROL_REFERENCE = 673

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

RAW_CONTROLS = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]


# ---------------------------------------------------------------------------
# Target generation
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


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def generate_prime_population(lo: int = 2_000_000,
                              hi: int = 4_200_000) -> list[int]:
    """
    Simple sieve of Eratosthenes.
    """
    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(hi)) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [i for i in range(lo, hi + 1) if sieve[i]]


def generate_targets(primes: list[int],
                     count: int,
                     seed: int) -> list[Target]:
    rng = random.Random(seed)

    # Deterministic but reasonably separated pairs.
    targets: list[Target] = []
    used_pairs: set[tuple[int, int]] = set()

    while len(targets) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pair = (p, q)
        if pair in used_pairs:
            continue

        # Keep the same general regime as the earlier experiments.
        if not (2_000_000 <= p <= 4_200_000):
            continue
        if not (2_000_000 <= q <= 4_200_000):
            continue

        used_pairs.add(pair)
        targets.append(Target(len(targets) + 1, p, q))

    return targets


# ---------------------------------------------------------------------------
# Modular helpers
# ---------------------------------------------------------------------------

def legendre_symbol(a: int, p: int) -> int:
    """
    Returns:
        +1 quadratic residue
         0 zero
        -1 non-residue
    """
    a %= p
    if a == 0:
        return 0

    x = pow(a, (p - 1) // 2, p)
    if x == 1:
        return 1
    if x == p - 1:
        return -1
    raise RuntimeError("Unexpected Legendre result")


def roots_x2_x1(ell: int) -> tuple[int, int] | None:
    """
    Find the two nontrivial roots of x^2+x+1 mod ell.

    For odd prime ell this exists exactly when ell == 1 mod 3.
    """
    if ell == 2 or ell % 3 != 1:
        return None

    roots = []
    for x in range(1, ell):
        if (x * x + x + 1) % ell == 0:
            roots.append(x)

    if len(roots) != 2:
        raise RuntimeError(
            f"{ell} was expected to have exactly two roots, got {roots}"
        )

    return tuple(sorted(roots))


def omega_for_roots(roots: tuple[int, int], ell: int) -> int:
    zeta, zeta2 = roots
    return (zeta - zeta2) % ell


def orientation_epsilon(n: int,
                        ell: int,
                        roots: tuple[int, int]) -> int | None:
    """
    epsilon = (d1-d2) / omega in F_ell.

    If either discriminant is zero, return None (degenerate).
    """
    zeta, zeta2 = roots
    omega = (zeta - zeta2) % ell

    if omega == 0:
        raise RuntimeError("omega must be nonzero")

    d1 = (zeta - 4 * n) % ell
    d2 = (zeta2 - 4 * n) % ell

    if d1 == 0 or d2 == 0:
        return None

    eps = ((d1 - d2) * pow(omega, -1, ell)) % ell

    if eps == 1:
        return 1
    if eps == ell - 1:
        return -1

    raise RuntimeError(
        f"Orientation failure ell={ell}, n={n}, "
        f"d1={d1}, d2={d2}, omega={omega}, eps={eps}"
    )


# ---------------------------------------------------------------------------
# N-only feature extraction
# ---------------------------------------------------------------------------

def cubic_class(x: int, ell: int) -> int:
    """
    Classifies a nonzero element using x^((ell-1)/3).

    0 is returned for zero.
    1/2 are the two nontrivial cubic-character classes.
    """
    x %= ell

    if x == 0:
        return 0

    y = pow(x, (ell - 1) // 3, ell)

    if y == 1:
        return 0

    roots = roots_x2_x1(ell)
    if roots is None:
        raise RuntimeError(
            f"cubic_class called on invalid ell={ell}"
        )

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
    All features are N-only.

    No factor p/q/s enters here.
    """
    values = [
        n,
        n + 1,
        n - 1,
        4 * n + 1,
        16 * n * n + 4 * n + 1,  # P
    ]

    features: list[float] = []

    # Legendre features.
    for x in values:
        chi = legendre_symbol(x, ell)
        features.append(float(chi))

    # Cubic classes as one-hot variables.
    cubic_values = [
        n,
        n + 1,
        4 * n + 1,
    ]

    for x in cubic_values:
        c = cubic_class(x, ell)
        features.extend([
            1.0 if c == 0 else 0.0,
            1.0 if c == 1 else 0.0,
            1.0 if c == 2 else 0.0,
        ])

    return features


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

@dataclass
class Sample:
    target_id: int
    ell: int
    x: list[float]
    y: int


def build_orientation_tables(moduli: Iterable[int]) -> dict[int, tuple[int, int]]:
    tables: dict[int, tuple[int, int]] = {}

    for ell in moduli:
        roots = roots_x2_x1(ell)
        if roots is None:
            continue
        tables[ell] = roots

    return tables


def build_relative_dataset(
    targets: list[Target],
    moduli: list[int],
    reference_ell: int,
    root_tables: dict[int, tuple[int, int]],
) -> tuple[list[Sample], int, int]:
    """
    For each target and ell != reference:
        rho = epsilon_ell * epsilon_reference.

    Samples where either orientation is degenerate are discarded.
    """
    samples: list[Sample] = []
    degenerate = 0

    if reference_ell not in root_tables:
        raise RuntimeError("Reference modulus has no valid roots")

    for t in targets:
        eps_ref = orientation_epsilon(
            t.n,
            reference_ell,
            root_tables[reference_ell],
        )

        for ell in moduli:
            if ell == reference_ell:
                continue

            eps_ell = orientation_epsilon(
                t.n,
                ell,
                root_tables[ell],
            )

            if eps_ref is None or eps_ell is None:
                degenerate += 1
                continue

            y = eps_ell * eps_ref
            assert y in (-1, 1)

            samples.append(
                Sample(
                    target_id=t.idx,
                    ell=ell,
                    x=n_features(t.n, ell),
                    y=y,
                )
            )

    positive = sum(s.y == 1 for s in samples)
    negative = sum(s.y == -1 for s in samples)

    return samples, positive, negative


# ---------------------------------------------------------------------------
# Tiny pure-Python ridge classifier
# ---------------------------------------------------------------------------

def mat_transpose(A: list[list[float]]) -> list[list[float]]:
    return [list(col) for col in zip(*A)]


def mat_mul(A: list[list[float]],
            B: list[list[float]]) -> list[list[float]]:
    rows = len(A)
    cols = len(B[0])
    inner = len(B)

    out = [[0.0] * cols for _ in range(rows)]

    for i in range(rows):
        for k in range(inner):
            aik = A[i][k]
            if aik == 0.0:
                continue
            for j in range(cols):
                out[i][j] += aik * B[k][j]

    return out


def solve_linear_system(A: list[list[float]],
                        b: list[float]) -> list[float]:
    """
    Gaussian elimination with partial pivoting.
    """
    n = len(A)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for col in range(n):
        pivot = max(
            range(col, n),
            key=lambda r: abs(M[r][col])
        )

        if abs(M[pivot][col]) < 1e-14:
            raise RuntimeError("Singular linear system")

        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]

        pivot_value = M[col][col]

        for j in range(col, n + 1):
            M[col][j] /= pivot_value

        for r in range(n):
            if r == col:
                continue

            factor = M[r][col]
            if factor == 0.0:
                continue

            for j in range(col, n + 1):
                M[r][j] -= factor * M[col][j]

    return [M[i][n] for i in range(n)]


class RidgeBinary:
    def __init__(self, alpha: float = 1.0):
        self.alpha = alpha
        self.mean: list[float] | None = None
        self.scale: list[float] | None = None
        self.weights: list[float] | None = None

    def _standardize_fit(
        self,
        X: list[list[float]]
    ) -> list[list[float]]:
        d = len(X[0])
        self.mean = [
            statistics.fmean(row[j] for row in X)
            for j in range(d)
        ]

        scale: list[float] = []
        for j in range(d):
            vals = [row[j] for row in X]
            m = self.mean[j]
            var = statistics.fmean((v - m) ** 2 for v in vals)
            s = math.sqrt(var)

            if s < 1e-12:
                s = 1.0

            scale.append(s)

        self.scale = scale

        return [
            [
                (row[j] - self.mean[j]) / self.scale[j]
                for j in range(d)
            ]
            for row in X
        ]

    def _standardize_apply(
        self,
        X: list[list[float]]
    ) -> list[list[float]]:
        assert self.mean is not None
        assert self.scale is not None

        return [
            [
                (row[j] - self.mean[j]) / self.scale[j]
                for j in range(len(row))
            ]
            for row in X
        ]

    def fit(self, X: list[list[float]], y: list[int]) -> None:
        Xs = self._standardize_fit(X)

        # Add intercept.
        Z = [[1.0] + row for row in Xs]

        p = len(Z[0])

        Xt = mat_transpose(Z)
        XtX = mat_mul(Xt, Z)

        # Ridge regularization except intercept.
        for i in range(1, p):
            XtX[i][i] += self.alpha

        Xty = [0.0] * p
        for i in range(p):
            Xty[i] = sum(
                Xt[r][i] * float(y[r])
                for r in range(len(y))
            )

        self.weights = solve_linear_system(XtX, Xty)

    def predict_score(self, X: list[list[float]]) -> list[float]:
        assert self.weights is not None

        Xs = self._standardize_apply(X)
        Z = [[1.0] + row for row in Xs]

        return [
            sum(w * x for w, x in zip(self.weights, row))
            for row in Z
        ]

    def predict(self, X: list[list[float]]) -> list[int]:
        return [
            1 if s >= 0.0 else -1
            for s in self.predict_score(X)
        ]


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def accuracy(model: RidgeBinary,
             samples: list[Sample]) -> float:
    if not samples:
        return float("nan")

    X = [s.x for s in samples]
    y = [s.y for s in samples]

    pred = model.predict(X)

    return sum(a == b for a, b in zip(y, pred)) / len(y)


def constant_accuracy(samples: list[Sample]) -> float:
    if not samples:
        return float("nan")

    pos = sum(s.y == 1 for s in samples)
    neg = len(samples) - pos

    return max(pos, neg) / len(samples)


def split_by_targets(
    samples: list[Sample],
    train_ids: set[int],
) -> tuple[list[Sample], list[Sample]]:
    train = [s for s in samples if s.target_id in train_ids]
    test = [s for s in samples if s.target_id not in train_ids]
    return train, test


def fit_and_score(
    train_samples: list[Sample],
    test_samples: list[Sample],
    alpha: float = 2.0,
) -> tuple[float, float]:
    model = RidgeBinary(alpha=alpha)

    X = [s.x for s in train_samples]
    y = [s.y for s in train_samples]

    if len(set(y)) < 2:
        return constant_accuracy(train_samples), constant_accuracy(test_samples)

    model.fit(X, y)

    return (
        accuracy(model, train_samples),
        accuracy(model, test_samples),
    )


# ---------------------------------------------------------------------------
# Permutation null
# ---------------------------------------------------------------------------

def permutation_pvalue(
    train_samples: list[Sample],
    test_samples: list[Sample],
    observed: float,
    permutations: int,
    seed: int,
) -> float:
    rng = random.Random(seed)

    y_original = [s.y for s in train_samples]

    exceed = 0

    for _ in range(permutations):
        shuffled = y_original[:]
        rng.shuffle(shuffled)

        perm_train = [
            Sample(
                target_id=s.target_id,
                ell=s.ell,
                x=s.x,
                y=y,
            )
            for s, y in zip(train_samples, shuffled)
        ]

        _, score = fit_and_score(
            perm_train,
            test_samples,
        )

        if score >= observed:
            exceed += 1

    return (exceed + 1) / (permutations + 1)


# ---------------------------------------------------------------------------
# Target and modulus splitting
# ---------------------------------------------------------------------------

def deterministic_target_split(
    targets: list[Target],
    seed: int,
    train_count: int,
) -> tuple[set[int], set[int]]:
    ids = [t.idx for t in targets]

    rng = random.Random(seed)
    rng.shuffle(ids)

    train_ids = set(ids[:train_count])
    test_ids = set(ids[train_count:])

    return train_ids, test_ids


def loo_modulus(
    samples: list[Sample],
    moduli: list[int],
    train_ids: set[int],
    alpha: float = 2.0,
) -> list[tuple[int, float, int]]:
    results = []

    for ell in moduli:
        train = [
            s for s in samples
            if s.target_id in train_ids and s.ell != ell
        ]

        test = [
            s for s in samples
            if s.target_id in train_ids and s.ell == ell
        ]

        if not test or len(set(s.y for s in train)) < 2:
            results.append((ell, float("nan"), len(test)))
            continue

        model = RidgeBinary(alpha=alpha)
        model.fit([s.x for s in train], [s.y for s in train])

        results.append((ell, accuracy(model, test), len(test)))

    return results


def two_dimensional_holdout(
    cyclo_samples: list[Sample],
    control_samples: list[Sample],
    train_ids: set[int],
    heldout_cyclo: set[int],
    heldout_control: set[int],
    alpha: float = 2.0,
) -> tuple[float, float]:

    def run(
        samples: list[Sample],
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
            return constant_accuracy(test)

        model = RidgeBinary(alpha=alpha)
        model.fit([s.x for s in train], [s.y for s in train])

        return accuracy(model, test)

    return (
        run(cyclo_samples, heldout_cyclo),
        run(control_samples, heldout_control),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 65")
    print("REFERENCE-ANCHORED RELATIVE ORIENTATION")
    print("N-ONLY PREDICTION OF EPSILON_ELL * EPSILON_REFERENCE")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("LEAVE-ONE-MODULUS-OUT")
    print("PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. Prime population + targets
    # ------------------------------------------------------------------

    p0 = time.perf_counter()

    primes = generate_prime_population()

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {time.perf_counter() - p0:.6f}s")

    targets = generate_targets(
        primes,
        NUM_TARGETS,
        seed=SEED,
    )

    print("\n2. TARGET SUMMARY")
    print("-" * 78)
    print(f"total targets = {len(targets)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    print("... remaining targets omitted")

    # ------------------------------------------------------------------
    # 3. Moduli
    # ------------------------------------------------------------------

    valid_controls = [
        ell for ell in RAW_CONTROLS
        if roots_x2_x1(ell) is not None
    ]

    excluded_controls = [
        ell for ell in RAW_CONTROLS
        if roots_x2_x1(ell) is None
    ]

    print("\n3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic      = {CYCLOTOMIC}")
    print(f"valid controls  = {valid_controls}")
    print(f"excluded controls= {excluded_controls}")

    if CYCLO_REFERENCE not in CYCLOTOMIC:
        raise RuntimeError("Cyclotomic reference missing")

    if CONTROL_REFERENCE not in valid_controls:
        raise RuntimeError("Control reference invalid")

    # ------------------------------------------------------------------
    # 4. Root validation
    # ------------------------------------------------------------------

    all_moduli = CYCLOTOMIC + valid_controls
    root_tables = build_orientation_tables(all_moduli)

    failures = 0

    print("\n4. ROOT / OMEGA VALIDATION")
    print("-" * 78)

    for ell, roots in root_tables.items():
        zeta, zeta2 = roots
        omega = omega_for_roots(roots, ell)

        if omega == 0:
            failures += 1

        ok = (
            (zeta * zeta + zeta + 1) % ell == 0
            and (zeta2 * zeta2 + zeta2 + 1) % ell == 0
            and (zeta != zeta2)
        )

        if not ok:
            failures += 1

        print(
            f"ell={ell:5d} roots={roots} "
            f"omega={omega:5d} identity_ok={ok}"
        )

    print(f"validation failures = {failures}")

    if failures:
        raise RuntimeError("Root validation failed")

    # ------------------------------------------------------------------
    # 5. Build relative orientation datasets
    # ------------------------------------------------------------------

    print("\n5. REFERENCE-ANCHORED DATASET")
    print("-" * 78)
    print(f"cyclotomic reference = {CYCLO_REFERENCE}")
    print(f"control reference    = {CONTROL_REFERENCE}")

    cyclo_samples, cpos, cneg = build_relative_dataset(
        targets,
        CYCLOTOMIC,
        CYCLO_REFERENCE,
        root_tables,
    )

    control_samples, rpos, rneg = build_relative_dataset(
        targets,
        valid_controls,
        CONTROL_REFERENCE,
        root_tables,
    )

    print(
        f"cyclotomic samples = {len(cyclo_samples)} "
        f"positive={cpos} negative={cneg}"
    )
    print(
        f"control samples    = {len(control_samples)} "
        f"positive={rpos} negative={rneg}"
    )

    print(
        f"cyclotomic constant baseline = "
        f"{constant_accuracy(cyclo_samples):.6f}"
    )
    print(
        f"control constant baseline     = "
        f"{constant_accuracy(control_samples):.6f}"
    )

    # ------------------------------------------------------------------
    # 6. Target split
    # ------------------------------------------------------------------

    train_ids, test_ids = deterministic_target_split(
        targets,
        seed=SEED + 1,
        train_count=TRAIN_TARGETS,
    )

    print("\n6. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_ids)}")
    print(f"test targets     = {len(test_ids)}")
    print(f"test ids         = {sorted(test_ids)}")

    ctrain, ctest = split_by_targets(cyclo_samples, train_ids)
    rtrain, rtest = split_by_targets(control_samples, train_ids)

    # ------------------------------------------------------------------
    # 7. Target holdout
    # ------------------------------------------------------------------

    c_train_acc, c_test_acc = fit_and_score(
        ctrain, ctest, alpha=2.0
    )
    r_train_acc, r_test_acc = fit_and_score(
        rtrain, rtest, alpha=2.0
    )

    print("\n7. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"cyclotomic train = {c_train_acc:.6f} "
        f"test = {c_test_acc:.6f}"
    )
    print(
        f"control    train = {r_train_acc:.6f} "
        f"test = {r_test_acc:.6f}"
    )
    print(
        f"test delta = {c_test_acc - r_test_acc:+.6f}"
    )

    # ------------------------------------------------------------------
    # 8. Repeated target holdout
    # ------------------------------------------------------------------

    print("\n8. REPEATED TARGET HOLDOUT")
    print("-" * 78)

    cyclo_scores = []
    control_scores = []

    for repeat in range(10):
        tr_ids, te_ids = deterministic_target_split(
            targets,
            seed=SEED + 100 + repeat,
            train_count=TRAIN_TARGETS,
        )

        ct, cv = split_by_targets(cyclo_samples, tr_ids)
        rt, rv = split_by_targets(control_samples, tr_ids)

        _, cs = fit_and_score(ct, cv, alpha=2.0)
        _, rs = fit_and_score(rt, rv, alpha=2.0)

        cyclo_scores.append(cs)
        control_scores.append(rs)

    print(
        f"cyclotomic mean   = {statistics.fmean(cyclo_scores):.6f}"
    )
    print(
        f"cyclotomic median = {statistics.median(cyclo_scores):.6f}"
    )
    print(
        f"cyclotomic min    = {min(cyclo_scores):.6f}"
    )
    print(
        f"cyclotomic max    = {max(cyclo_scores):.6f}"
    )
    print(
        f"control mean      = {statistics.fmean(control_scores):.6f}"
    )
    print(
        f"control median    = {statistics.median(control_scores):.6f}"
    )

    # ------------------------------------------------------------------
    # 9. Leave-one-modulus-out
    # ------------------------------------------------------------------

    print("\n9. LEAVE-ONE-MODULUS-OUT")
    print("-" * 78)

    cyclo_nonref = [e for e in CYCLOTOMIC if e != CYCLO_REFERENCE]
    control_nonref = [e for e in valid_controls if e != CONTROL_REFERENCE]

    cyclo_loo = loo_modulus(
        cyclo_samples,
        cyclo_nonref,
        train_ids,
        alpha=2.0,
    )

    control_loo = loo_modulus(
        control_samples,
        control_nonref,
        train_ids,
        alpha=2.0,
    )

    for ell, acc, count in cyclo_loo:
        print(
            f"cyclo ell={ell:5d} "
            f"accuracy={acc:.6f} n={count}"
        )

    print(
        f"cyclotomic LOO mean = "
        f"{statistics.fmean(a for _, a, _ in cyclo_loo if not math.isnan(a)):.6f}"
    )

    for ell, acc, count in control_loo:
        print(
            f"control ell={ell:5d} "
            f"accuracy={acc:.6f} n={count}"
        )

    print(
        f"control LOO mean = "
        f"{statistics.fmean(a for _, a, _ in control_loo if not math.isnan(a)):.6f}"
    )

    # ------------------------------------------------------------------
    # 10. Two-dimensional holdout
    # ------------------------------------------------------------------

    # Hold out several whole moduli, while simultaneously holding out
    # unseen targets.
    cyclo_heldout_moduli = {61, 127, 307}
    control_heldout_moduli = {6211, 7951}

    c2d, r2d = two_dimensional_holdout(
        cyclo_samples,
        control_samples,
        train_ids,
        cyclo_heldout_moduli,
        control_heldout_moduli,
        alpha=2.0,
    )

    print("\n10. TWO-DIMENSIONAL HOLDOUT")
    print("-" * 78)
    print(
        f"cyclotomic 2D accuracy = {c2d:.6f} "
        f"heldout_moduli={sorted(cyclo_heldout_moduli)}"
    )
    print(
        f"control 2D accuracy    = {r2d:.6f} "
        f"heldout_moduli={sorted(control_heldout_moduli)}"
    )
    print(f"2D delta = {c2d - r2d:+.6f}")

    # ------------------------------------------------------------------
    # 11. Fixed N-only baselines
    # ------------------------------------------------------------------

    print("\n11. FIXED N-ONLY BASELINES")
    print("-" * 78)

    def baseline_feature(
        samples: list[Sample],
        feature_index: int,
    ) -> float:
        usable = []

        for s in samples:
            v = s.x[feature_index]
            if v > 0:
                pred = 1
            elif v < 0:
                pred = -1
            else:
                continue

            usable.append(pred == s.y)

        if not usable:
            return float("nan")

        return sum(usable) / len(usable)

    for name, idx in [
        ("chi(n)", 0),
        ("chi(n+1)", 1),
        ("chi(n-1)", 2),
        ("chi(4n+1)", 3),
        ("chi(P)", 4),
    ]:
        ca = baseline_feature(ctest, idx)
        ra = baseline_feature(rtest, idx)

        print(
            f"{name:12s} C={ca:.6f} R={ra:.6f} "
            f"delta={ca-ra:+.6f}"
        )

    # ------------------------------------------------------------------
    # 12. Permutation null
    # ------------------------------------------------------------------

    print("\n12. PERMUTATION NULL")
    print("-" * 78)

    p_value = permutation_pvalue(
        ctrain,
        ctest,
        c_test_acc,
        permutations=PERMUTATIONS,
        seed=SEED + 9000,
    )

    print(f"observed cyclotomic test accuracy = {c_test_acc:.6f}")
    print(f"permutations = {PERMUTATIONS}")
    print(f"empirical p-value = {p_value:.6f}")

    # ------------------------------------------------------------------
    # 13. Per-modulus target holdout
    # ------------------------------------------------------------------

    print("\n13. PER-MODULUS TARGET HOLDOUT")
    print("-" * 78)

    for ell in cyclo_nonref:
        tr = [
            s for s in ctrain
            if s.ell == ell
        ]
        te = [
            s for s in ctest
            if s.ell == ell
        ]

        if len(tr) < 10 or len(set(s.y for s in tr)) < 2:
            print(
                f"ell={ell:5d}: insufficient "
                f"train={len(tr)}"
            )
            continue

        model = RidgeBinary(alpha=2.0)
        model.fit([s.x for s in tr], [s.y for s in tr])

        acc = accuracy(model, te)

        print(
            f"ell={ell:5d} "
            f"train={len(tr):3d} "
            f"test={len(te):3d} "
            f"accuracy={acc:.6f}"
        )

    # ------------------------------------------------------------------
    # 14. Direct sanity check:
    #     does relative orientation contain actual variation?
    # ------------------------------------------------------------------

    print("\n14. RELATIVE-ORIENTATION SANITY CHECK")
    print("-" * 78)

    if cpos == 0 or cneg == 0:
        print(
            "WARNING: cyclotomic relative-orientation label is constant."
        )
        print(
            "This means 100% classification would be trivial."
        )
    else:
        print(
            "Cyclotomic relative orientation has both signs."
        )

    if rpos == 0 or rneg == 0:
        print(
            "WARNING: control relative-orientation label is constant."
        )
        print(
            "Control 100% accuracy would be trivial."
        )
    else:
        print(
            "Control relative orientation has both signs."
        )

    # ------------------------------------------------------------------
    # 15. Final diagnostic
    # ------------------------------------------------------------------

    cyclo_repeat = statistics.fmean(cyclo_scores)
    control_repeat = statistics.fmean(control_scores)

    cyclo_loo_vals = [
        a for _, a, _ in cyclo_loo
        if not math.isnan(a)
    ]
    control_loo_vals = [
        a for _, a, _ in control_loo
        if not math.isnan(a)
    ]

    cyclo_loo_mean = statistics.fmean(cyclo_loo_vals)
    control_loo_mean = statistics.fmean(control_loo_vals)

    print("\n15. FINAL DIAGNOSTIC")
    print("-" * 78)
    print(
        f"target-holdout cyclotomic = {c_test_acc:.6f}"
    )
    print(
        f"target-holdout control    = {r_test_acc:.6f}"
    )
    print(
        f"target-holdout delta      = "
        f"{c_test_acc-r_test_acc:+.6f}"
    )
    print(
        f"repeated cyclotomic mean  = {cyclo_repeat:.6f}"
    )
    print(
        f"repeated control mean     = {control_repeat:.6f}"
    )
    print(
        f"cyclotomic LOO mean       = {cyclo_loo_mean:.6f}"
    )
    print(
        f"control LOO mean          = {control_loo_mean:.6f}"
    )
    print(
        f"2D cyclotomic             = {c2d:.6f}"
    )
    print(
        f"2D control                = {r2d:.6f}"
    )
    print(
        f"permutation p-value       = {p_value:.6f}"
    )

    print("\nINTERPRETATION")
    print("-" * 78)

    if cpos == 0 or cneg == 0:
        print(
            "TRIVIAL LABEL:"
            " relative orientation is constant."
        )
        print(
            "The experiment does not test orientation recovery."
        )
    elif (
        c_test_acc > 0.55
        and c_test_acc > r_test_acc + 0.05
        and cyclo_loo_mean > control_loo_mean + 0.03
        and c2d > 0.55
        and p_value < 0.05
    ):
        print(
            "POTENTIAL CYCLOTOMIC RELATIVE-ORIENTATION SIGNAL:"
        )
        print(
            "The relative sign survives target holdout, "
            "modulus holdout, 2D holdout, and permutation testing."
        )
        print(
            "This would justify a subsequent CRT-reconstruction experiment."
        )
    elif c_test_acc > 0.55 and c_test_acc > r_test_acc + 0.05:
        print(
            "TARGET-LEVEL SIGNAL ONLY:"
        )
        print(
            "There may be target-distribution correlation, "
            "but cross-modulus transfer is not established."
        )
    else:
        print(
            "NO ROBUST RELATIVE-ORIENTATION SIGNAL:"
        )
        print(
            "The N-only features do not demonstrate stable "
            "cross-modulus orientation recovery."
        )

    print("\nIMPORTANT")
    print("-" * 78)
    print(
        "This experiment deliberately predicts epsilon_ell * epsilon_reference."
    )
    print(
        "It does NOT use p, q, or s in the feature vector."
    )
    print(
        "The reference multiplication removes global branch ambiguity "
        "without collapsing the label to a trivial constant."
    )
    print(
        "A successful result would be considerably stronger than "
        "Experiment 64."
    )

    print("\n" + "=" * 78)
    print("EXPERIMENT 65 COMPLETE")
    print(
        f"total runtime = {time.perf_counter() - t0:.6f}s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()
