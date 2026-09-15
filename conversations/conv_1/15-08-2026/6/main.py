#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 64
EISENSTEIN ORIENTATION / CROSS-MODULUS SIGN-CONSISTENCY TEST
N-ONLY FEATURES VS LOCAL ORIENTATION LABEL
TARGET + MODULUS + 2D HOLDOUT
PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Purpose
-------
Experiment 63 showed that low-degree symmetric N-only expressions in

    S = d1 + d2
    P = d1*d2

do not generalize reliably across unseen moduli.

The exact identity is

    d1 + d2       = -8n - 1
    d1*d2         = 16n^2 + 4n + 1
    (d1-d2)^2     = -3.

For every ell == 1 mod 3, let zeta be a fixed nontrivial cube root
of unity:

    zeta^2 + zeta + 1 == 0.

Then

    omega = zeta - zeta^2

satisfies

    omega^2 = -3.

Therefore the local ordered pair has an orientation sign

    epsilon_ell = (d1 - d2) / omega in {+1,-1}.

The important question is NOT whether epsilon_ell itself can be read
from a symmetric polynomial in n.

Instead we test whether the COLLECTION of orientation signs across
different ell has a transferable structure that can be predicted from
N-only modular characters.

In particular we examine:

  * individual epsilon_ell prediction;
  * pairwise orientation products epsilon_a * epsilon_b;
  * normalized Eisenstein signs;
  * N-only Legendre/cubic character features;
  * target holdout;
  * leave-one-modulus-out;
  * two-dimensional holdout;
  * permutation null.

A positive result must survive unseen targets AND unseen ell.

IMPORTANT
---------
The local epsilon label depends on an ordering of the two roots. That
ordering is only a diagnostic/oracle label. The predictor itself gets
N-only features.

No factorization algorithm is claimed by this experiment.
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import isqrt
from random import Random
from statistics import mean, median
from typing import Dict, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 64001
TARGET_COUNT = 120

PRIME_LIMIT = 4_200_000

CYCLotomic = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

RAW_CONTROLS = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

VALID_CONTROLS = [p for p in RAW_CONTROLS if p % 3 == 1]

TEST_TARGET_FRACTION = 0.25
PERMUTATIONS = 500

# Candidate local N-only shifts.
SHIFTS = (-3, -2, -1, 0, 1, 2, 3)

# Candidate pair feature names.
PAIR_FEATURE_LIMIT = 40


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


@dataclass
class LocalState:
    ell: int
    r1: int
    r2: int
    d1: int
    d2: int
    epsilon: int
    chi_d1: int
    chi_d2: int
    chi_prod: int
    cubic_n: int
    cubic_np1: int
    cubic_nm1: int
    chi_n: int
    chi_np1: int
    chi_nm1: int
    chi_4n1: int


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def sieve_primes(limit: int) -> List[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, v in enumerate(sieve) if v]


# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

def legendre_symbol(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    x = pow(a, (p - 1) // 2, p)
    return 1 if x == 1 else -1


def roots_cyclotomic(ell: int) -> Tuple[int, int]:
    """
    Return the two roots of x^2 + x + 1 mod ell.

    Only primes ell == 1 mod 3 have two nontrivial roots.
    """
    if ell % 3 != 1:
        raise ValueError(
            f"{ell} is not orientation-valid: ell must be 1 mod 3"
        )

    roots = [x for x in range(2, ell) if (x * x + x + 1) % ell == 0]

    if len(roots) != 2:
        raise RuntimeError(
            f"{ell}: expected two nontrivial roots, got {roots}"
        )

    return tuple(roots)  # type: ignore[return-value]


def cubic_class(a: int, ell: int) -> int:
    """
    For ell == 1 mod 3, classify a nonzero residue into the three
    cubic-character classes using a fixed primitive cube root zeta.

    We use the discrete-log-free criterion:

        a^((ell-1)/3)

    which lands in {1, zeta, zeta^2}.
    """
    a %= ell

    if a == 0:
        return -1

    zeta, zeta2 = roots_cyclotomic(ell)
    value = pow(a, (ell - 1) // 3, ell)

    if value == 1:
        return 0
    if value == zeta:
        return 1
    if value == zeta2:
        return 2

    raise RuntimeError(
        f"unexpected cubic value ell={ell} a={a}: {value}"
    )


def chi(a: int, ell: int) -> int:
    return legendre_symbol(a, ell)


# ---------------------------------------------------------------------------
# Target generation
# ---------------------------------------------------------------------------

def generate_targets(
    primes: Sequence[int],
    count: int,
    rng: Random,
) -> List[Target]:
    """
    Generate semantically matched semiprimes with even factor sums.

    Both odd primes imply even s automatically.
    """
    out: List[Target] = []

    # Avoid tiny/huge imbalance.
    eligible = [
        p for p in primes
        if 2_000_000 <= p <= 4_200_000
    ]

    seen = set()

    while len(out) < count:
        p = eligible[rng.randrange(len(eligible))]
        q = eligible[rng.randrange(len(eligible))]

        if p == q:
            continue

        key = tuple(sorted((p, q)))
        if key in seen:
            continue

        seen.add(key)
        out.append(Target(len(out) + 1, p, q))

    return out


# ---------------------------------------------------------------------------
# Local orientation
# ---------------------------------------------------------------------------

def normalized_eisenstein_orientation(
    d1: int,
    d2: int,
    zeta: int,
    zeta2: int,
    ell: int,
) -> int:
    """
    omega = zeta - zeta^2, omega^2 = -3.

    Since (d1-d2)^2 = -3, the quotient

        (d1-d2) / omega

    must be +1 or -1.

    This uses the inverse of omega modulo ell.
    """
    omega = (zeta - zeta2) % ell

    if omega == 0:
        raise RuntimeError(f"omega vanished mod {ell}")

    inv_omega = pow(omega, -1, ell)

    eps = ((d1 - d2) % ell) * inv_omega % ell

    if eps == 1:
        return 1
    if eps == ell - 1:
        return -1

    raise RuntimeError(
        f"orientation quotient was not ±1: ell={ell}, eps={eps}"
    )


def local_state(n: int, ell: int) -> LocalState:
    r1, r2 = roots_cyclotomic(ell)

    # F(w)=w^2+w+1=0
    # d(w)=w^2 - 4n = -w - 1 - 4n
    d1 = (-r1 - 1 - 4 * n) % ell
    d2 = (-r2 - 1 - 4 * n) % ell

    zeta, zeta2 = r1, r2

    epsilon = normalized_eisenstein_orientation(
        d1, d2, zeta, zeta2, ell
    )

    # Verify the structural identity.
    if ((d1 - d2) * (d1 - d2) + 3) % ell != 0:
        raise AssertionError(
            f"(d1-d2)^2 != -3 mod {ell}"
        )

    cp_n = cubic_class(n, ell)
    cp_np1 = cubic_class(n + 1, ell)
    cp_nm1 = cubic_class(n - 1, ell)

    return LocalState(
        ell=ell,
        r1=r1,
        r2=r2,
        d1=d1,
        d2=d2,
        epsilon=epsilon,
        chi_d1=chi(d1, ell),
        chi_d2=chi(d2, ell),
        chi_prod=chi(d1, ell) * chi(d2, ell),
        cubic_n=cp_n,
        cubic_np1=cp_np1,
        cubic_nm1=cp_nm1,
        chi_n=chi(n, ell),
        chi_np1=chi(n + 1, ell),
        chi_nm1=chi(n - 1, ell),
        chi_4n1=chi(4 * n + 1, ell),
    )


# ---------------------------------------------------------------------------
# Feature engineering
# ---------------------------------------------------------------------------

def local_feature_vector(st: LocalState) -> Tuple[int, ...]:
    """
    N-only features.

    No d1/d2, no p/q, no s.
    """
    vals: List[int] = []

    # Quadratic characters.
    for x in (
        st.chi_n,
        st.chi_np1,
        st.chi_nm1,
        st.chi_4n1,
    ):
        vals.extend([0 if x == -1 else 1])

    # One-hot cubic classes.
    for c in (
        st.cubic_n,
        st.cubic_np1,
        st.cubic_nm1,
    ):
        vals.extend([
            int(c == 0),
            int(c == 1),
            int(c == 2),
        ])

    return tuple(vals)


def pair_feature_vector(
    a: LocalState,
    b: LocalState,
) -> Tuple[int, ...]:
    """
    N-only pair features.

    The target label is epsilon_a * epsilon_b.

    Features intentionally avoid local orientation labels.
    """
    fa = local_feature_vector(a)
    fb = local_feature_vector(b)

    out: List[int] = list(fa)
    out.extend(fb)

    # Cross products of the N-only scalar characters.
    scalar_a = (
        a.chi_n,
        a.chi_np1,
        a.chi_nm1,
        a.chi_4n1,
    )

    scalar_b = (
        b.chi_n,
        b.chi_np1,
        b.chi_nm1,
        b.chi_4n1,
    )

    for xa in scalar_a:
        for xb in scalar_b:
            out.append(xa * xb)

    # Modulus-only indicator, but NOT modulus identity.
    out.append(int(a.ell % 3 == 1))
    out.append(int(b.ell % 3 == 1))

    return tuple(out)


# ---------------------------------------------------------------------------
# Simple regularized linear classifier
# ---------------------------------------------------------------------------

def sigmoid(z: float) -> float:
    if z >= 0:
        e = pow(2.718281828459045, -z)
        return 1.0 / (1.0 + e)
    e = pow(2.718281828459045, z)
    return e / (1.0 + e)


def standardize(
    X: Sequence[Sequence[float]],
) -> Tuple[List[List[float]], List[float], List[float]]:
    if not X:
        return [], [], []

    m = len(X[0])
    means = [0.0] * m
    scales = [1.0] * m

    n = len(X)

    for j in range(m):
        means[j] = sum(row[j] for row in X) / n

    for j in range(m):
        var = sum((row[j] - means[j]) ** 2 for row in X) / n
        scales[j] = var ** 0.5
        if scales[j] < 1e-12:
            scales[j] = 1.0

    Z = [
        [(row[j] - means[j]) / scales[j] for j in range(m)]
        for row in X
    ]

    return Z, means, scales


def apply_standardize(
    X: Sequence[Sequence[float]],
    means: Sequence[float],
    scales: Sequence[float],
) -> List[List[float]]:
    return [
        [(row[j] - means[j]) / scales[j] for j in range(len(row))]
        for row in X
    ]


def fit_logistic(
    X: Sequence[Sequence[float]],
    y: Sequence[int],
    l2: float = 5.0,
    steps: int = 700,
    lr: float = 0.03,
) -> Tuple[List[float], float]:
    """
    Small deterministic logistic regression implementation.

    No sklearn dependency.
    """
    if not X:
        raise ValueError("empty training matrix")

    n = len(X)
    m = len(X[0])

    w = [0.0] * m
    b = 0.0

    for _ in range(steps):
        gw = [0.0] * m
        gb = 0.0

        for i in range(n):
            z = b + sum(w[j] * X[i][j] for j in range(m))
            pr = sigmoid(z)
            err = pr - y[i]

            gb += err

            for j in range(m):
                gw[j] += err * X[i][j]

        # Average gradient + L2.
        gb /= n

        for j in range(m):
            gw[j] = gw[j] / n + l2 * w[j]

        b -= lr * gb

        for j in range(m):
            w[j] -= lr * gw[j]

    return w, b


def predict_prob(
    X: Sequence[Sequence[float]],
    w: Sequence[float],
    b: float,
) -> List[float]:
    out = []

    for row in X:
        z = b + sum(w[j] * row[j] for j in range(len(row)))
        out.append(sigmoid(z))

    return out


def predict_label(
    X: Sequence[Sequence[float]],
    w: Sequence[float],
    b: float,
) -> List[int]:
    return [1 if p >= 0.5 else -1 for p in predict_prob(X, w, b)]


def accuracy(
    truth: Sequence[int],
    pred: Sequence[int],
) -> float:
    if not truth:
        return float("nan")

    return sum(a == b for a, b in zip(truth, pred)) / len(truth)


# ---------------------------------------------------------------------------
# Sample construction
# ---------------------------------------------------------------------------

@dataclass
class Sample:
    target_id: int
    ell_a: int
    ell_b: int
    X: Tuple[int, ...]
    y: int


def build_pair_samples(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> List[Sample]:
    out: List[Sample] = []

    for t in targets:
        states = []

        for ell in moduli:
            st = local_state(t.n, ell)

            # Pair data are only interesting where both branches are
            # non-degenerate in the ordinary Legendre sense.
            if st.chi_d1 == 0 or st.chi_d2 == 0:
                continue

            states.append(st)

        for a, b in combinations(states, 2):
            y = a.epsilon * b.epsilon

            out.append(
                Sample(
                    target_id=t.idx,
                    ell_a=a.ell,
                    ell_b=b.ell,
                    X=pair_feature_vector(a, b),
                    y=y,
                )
            )

    return out


# ---------------------------------------------------------------------------
# Split utilities
# ---------------------------------------------------------------------------

def split_target_ids(
    target_ids: Sequence[int],
    fraction: float,
    rng: Random,
) -> Tuple[set[int], set[int]]:
    ids = list(target_ids)
    rng.shuffle(ids)

    n_test = max(1, int(round(len(ids) * fraction)))

    test_ids = set(ids[:n_test])
    train_ids = set(ids[n_test:])

    return train_ids, test_ids


def rows_for_targets(
    samples: Sequence[Sample],
    ids: set[int],
) -> List[Sample]:
    return [s for s in samples if s.target_id in ids]


def prepare_matrix(
    samples: Sequence[Sample],
) -> Tuple[List[List[float]], List[int]]:
    X = [list(map(float, s.X)) for s in samples]
    y = [1 if s.y == 1 else 0 for s in samples]
    return X, y


# ---------------------------------------------------------------------------
# Training / evaluation
# ---------------------------------------------------------------------------

def train_and_test(
    train_samples: Sequence[Sample],
    test_samples: Sequence[Sample],
) -> float:
    if not train_samples or not test_samples:
        return float("nan")

    Xtr, ytr01 = prepare_matrix(train_samples)
    Xte, yte01 = prepare_matrix(test_samples)

    ytr = [1 if y == 1 else -1 for y in ytr01]
    yte = [1 if y == 1 else -1 for y in yte01]

    Xtr_z, means, scales = standardize(Xtr)
    Xte_z = apply_standardize(Xte, means, scales)

    w, b = fit_logistic(
        Xtr_z,
        ytr01,
        l2=4.0,
        steps=650,
        lr=0.035,
    )

    pred01 = predict_label(Xte_z, w, b)
    pred = [1 if x == 1 else -1 for x in pred01]

    return accuracy(yte, pred)


# ---------------------------------------------------------------------------
# Direct baseline rules
# ---------------------------------------------------------------------------

def baseline_pair_rule(
    samples: Sequence[Sample],
    mode: str,
) -> float:
    """
    Test simple N-only cross-modulus product rules.
    """

    if not samples:
        return float("nan")

    correct = 0
    total = 0

    for s in samples:
        # Decode the first handful of scalar features.
        # feature layout:
        # [chi_n, chi_np1, chi_nm1, chi_4n1,
        #  3 one-hot cubic groups ...]
        x = s.X

        if mode == "chi_n":
            # Product of first scalar feature encoded as {-1,+1}.
            a = -1 if x[0] == 0 else 1
            b = -1 if x[4 + 0] == 0 else 1
            pred = a * b

        elif mode == "chi_np1":
            a = -1 if x[1] == 0 else 1
            b = -1 if x[4 + 1] == 0 else 1
            pred = a * b

        else:
            raise ValueError(mode)

        if pred == s.y:
            correct += 1

        total += 1

    return correct / total if total else float("nan")


# ---------------------------------------------------------------------------
# Pair permutation null
# ---------------------------------------------------------------------------

def permutation_test(
    train_samples: Sequence[Sample],
    test_samples: Sequence[Sample],
    observed: float,
    rng: Random,
    permutations: int,
) -> float:
    """
    Shuffle training labels while preserving all feature structure.
    """
    if not train_samples or not test_samples:
        return float("nan")

    Xtr, ytr = prepare_matrix(train_samples)
    Xte, yte = prepare_matrix(test_samples)

    Xtr_z, means, scales = standardize(Xtr)
    Xte_z = apply_standardize(Xte, means, scales)

    count = 0

    labels = list(ytr)

    for _ in range(permutations):
        rng.shuffle(labels)

        w, b = fit_logistic(
            Xtr_z,
            labels,
            l2=4.0,
            steps=450,
            lr=0.035,
        )

        pred = predict_label(Xte_z, w, b)
        hit = accuracy(
            [1 if y == 1 else -1 for y in yte],
            pred,
        )

        if hit >= observed - 1e-12:
            count += 1

    return (count + 1) / (permutations + 1)


# ---------------------------------------------------------------------------
# Modulus holdout
# ---------------------------------------------------------------------------

def leave_one_modulus_out(
    samples: Sequence[Sample],
    moduli: Sequence[int],
) -> Dict[int, float]:
    out = {}

    for ell in moduli:
        train = [
            s for s in samples
            if s.ell_a != ell and s.ell_b != ell
        ]

        test = [
            s for s in samples
            if s.ell_a == ell or s.ell_b == ell
        ]

        out[ell] = train_and_test(train, test)

    return out


# ---------------------------------------------------------------------------
# Structural diagnostics
# ---------------------------------------------------------------------------

def validate_orientation_identity(
    targets: Sequence[Target],
    moduli: Sequence[int],
) -> None:
    failures = 0

    for t in targets:
        for ell in moduli:
            st = local_state(t.n, ell)

            if st.epsilon not in (-1, 1):
                failures += 1
                continue

            if (
                ((st.d1 - st.d2) % ell)
                !=
                (st.epsilon * ((st.r1 - st.r2) % ell)) % ell
            ):
                failures += 1

    print(f"orientation identity failures = {failures}")
    print(f"status = {'PASS' if failures == 0 else 'FAIL'}")


def summarize_pair_balance(
    samples: Sequence[Sample],
) -> None:
    pos = sum(s.y == 1 for s in samples)
    neg = sum(s.y == -1 for s in samples)

    total = pos + neg

    print(f"pair samples = {total}")
    print(f"epsilon_pair +1 = {pos}")
    print(f"epsilon_pair -1 = {neg}")

    if total:
        print(f"positive fraction = {pos / total:.6f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    rng = Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 64")
    print("EISENSTEIN ORIENTATION / CROSS-MODULUS SIGN-CONSISTENCY TEST")
    print("N-ONLY FEATURES VS LOCAL ORIENTATION LABEL")
    print("TARGET + MODULUS + 2D HOLDOUT")
    print("PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. Prime population
    # ----------------------------------------------------------------------
    print("\n1. PRIME POPULATION")
    primes = sieve_primes(PRIME_LIMIT)
    print(f"prime population = {len(primes)}")

    # ----------------------------------------------------------------------
    # 2. Targets
    # ----------------------------------------------------------------------
    targets = generate_targets(primes, TARGET_COUNT, rng)

    print("\n2. TARGET SUMMARY")
    print(f"total targets = {len(targets)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining generated targets omitted")

    # ----------------------------------------------------------------------
    # 3. Modulus families
    # ----------------------------------------------------------------------
    print("\n3. MODULUS FAMILIES")
    print(f"cyclotomic = {CYCLotomic}")
    print(f"raw controls = {RAW_CONTROLS}")
    print(f"valid controls = {VALID_CONTROLS}")

    print(
        f"excluded controls = "
        f"{[x for x in RAW_CONTROLS if x not in VALID_CONTROLS]}"
    )

    # ----------------------------------------------------------------------
    # 4. Root/orientation validation
    # ----------------------------------------------------------------------
    print("\n4. ORIENTATION IDENTITY VALIDATION")
    validate_orientation_identity(targets[:10], CYCLotomic)

    # ----------------------------------------------------------------------
    # 5. Root table
    # ----------------------------------------------------------------------
    print("\n5. CYCLOTOMIC ROOTS")

    for ell in CYCLotomic:
        r1, r2 = roots_cyclotomic(ell)
        omega = (r1 - r2) % ell

        print(
            f"ell={ell:5d} "
            f"roots=({r1:5d},{r2:5d}) "
            f"omega=zeta-zeta^2={omega:5d} "
            f"omega^2 mod ell={omega * omega % ell:5d}"
        )

    # ----------------------------------------------------------------------
    # 6. Build local pair dataset
    # ----------------------------------------------------------------------
    print("\n6. CROSS-MODULUS PAIR DATASET")

    cyclo_samples = build_pair_samples(
        targets,
        CYCLotomic,
    )

    control_samples = build_pair_samples(
        targets,
        VALID_CONTROLS,
    )

    summarize_pair_balance(cyclo_samples)
    print(f"control pair samples = {len(control_samples)}")

    # ----------------------------------------------------------------------
    # 7. Target split
    # ----------------------------------------------------------------------
    target_ids = [t.idx for t in targets]

    train_ids, test_ids = split_target_ids(
        target_ids,
        TEST_TARGET_FRACTION,
        rng,
    )

    print("\n7. TARGET HOLDOUT")
    print(f"training targets = {len(train_ids)}")
    print(f"test targets = {len(test_ids)}")

    train_c = rows_for_targets(cyclo_samples, train_ids)
    test_c = rows_for_targets(cyclo_samples, test_ids)

    train_r = rows_for_targets(control_samples, train_ids)
    test_r = rows_for_targets(control_samples, test_ids)

    # ----------------------------------------------------------------------
    # 8. Target holdout
    # ----------------------------------------------------------------------
    print("\n8. TARGET HOLDOUT")

    acc_c = train_and_test(train_c, test_c)
    acc_r = train_and_test(train_r, test_r)

    print(f"cyclotomic accuracy = {acc_c:.6f}")
    print(f"control accuracy    = {acc_r:.6f}")
    print(f"difference          = {acc_c - acc_r:+.6f}")

    # ----------------------------------------------------------------------
    # 9. Simple baselines
    # ----------------------------------------------------------------------
    print("\n9. N-ONLY CROSS-MODULUS BASELINES")

    for name in ("chi_n", "chi_np1"):
        c = baseline_pair_rule(test_c, name)
        r = baseline_pair_rule(test_r, name)

        print(
            f"{name:12s} "
            f"C={c:.6f} R={r:.6f} delta={c-r:+.6f}"
        )

    # ----------------------------------------------------------------------
    # 10. Leave-one-modulus-out
    # ----------------------------------------------------------------------
    print("\n10. LEAVE-ONE-MODULUS-OUT")

    loo_c = leave_one_modulus_out(
        cyclo_samples,
        CYCLotomic,
    )

    loo_r = leave_one_modulus_out(
        control_samples,
        VALID_CONTROLS,
    )

    for ell in CYCLotomic:
        print(
            f"cyclo ell={ell:5d} accuracy={loo_c[ell]:.6f}"
        )

    for ell in VALID_CONTROLS:
        print(
            f"control ell={ell:5d} accuracy={loo_r[ell]:.6f}"
        )

    print(
        f"cyclotomic LOO mean   = {mean(loo_c.values()):.6f}"
    )
    print(
        f"cyclotomic LOO median = {median(loo_c.values()):.6f}"
    )
    print(
        f"control LOO mean      = {mean(loo_r.values()):.6f}"
    )
    print(
        f"control LOO median    = {median(loo_r.values()):.6f}"
    )

    # ----------------------------------------------------------------------
    # 11. Two-dimensional holdout
    # ----------------------------------------------------------------------
    heldout_c_mods = set(
        CYCLotomic[::3]
    )
    heldout_r_mods = set(
        VALID_CONTROLS[::3]
    )

    train2_c = [
        s for s in train_c
        if s.ell_a not in heldout_c_mods
        and s.ell_b not in heldout_c_mods
    ]

    test2_c = [
        s for s in test_c
        if s.ell_a in heldout_c_mods
        or s.ell_b in heldout_c_mods
    ]

    train2_r = [
        s for s in train_r
        if s.ell_a not in heldout_r_mods
        and s.ell_b not in heldout_r_mods
    ]

    test2_r = [
        s for s in test_r
        if s.ell_a in heldout_r_mods
        or s.ell_b in heldout_r_mods
    ]

    print("\n11. TWO-DIMENSIONAL HOLDOUT")
    print(f"heldout cyclotomic ell = {sorted(heldout_c_mods)}")
    print(f"heldout control ell    = {sorted(heldout_r_mods)}")

    acc2_c = train_and_test(train2_c, test2_c)
    acc2_r = train_and_test(train2_r, test2_r)

    print(f"cyclotomic 2D = {acc2_c:.6f}")
    print(f"control 2D    = {acc2_r:.6f}")
    print(f"delta         = {acc2_c - acc2_r:+.6f}")

    # ----------------------------------------------------------------------
    # 12. Pair-sign consistency / internal algebra
    # ----------------------------------------------------------------------
    print("\n12. PAIR SIGN CONSISTENCY")

    # For a,b,c:
    # (eps_a eps_b)(eps_b eps_c) = eps_a eps_c.
    # This relation is algebraically exact. We verify it on observed data.
    consistency_checks = 0
    consistency_failures = 0

    for t in targets[: min(20, len(targets))]:
        states = {
            ell: local_state(t.n, ell)
            for ell in CYCLotomic
        }

        for a, b, c in combinations(CYCLotomic, 3):
            lhs = states[a].epsilon * states[b].epsilon
            rhs = states[b].epsilon * states[c].epsilon
            predicted = lhs * rhs

            direct = states[a].epsilon * states[c].epsilon

            consistency_checks += 1

            if predicted != direct:
                consistency_failures += 1

    print(f"checks   = {consistency_checks}")
    print(f"failures = {consistency_failures}")
    print(
        f"status   = "
        f"{'PASS' if consistency_failures == 0 else 'FAIL'}"
    )

    # ----------------------------------------------------------------------
    # 13. Permutation null
    # ----------------------------------------------------------------------
    print("\n13. PERMUTATION NULL")

    observed = acc_c
    p_value = permutation_test(
        train_c,
        test_c,
        observed,
        rng,
        PERMUTATIONS,
    )

    print(f"observed cyclotomic accuracy = {observed:.6f}")
    print(f"permutations = {PERMUTATIONS}")
    print(f"empirical p-value = {p_value:.6f}")

    # ----------------------------------------------------------------------
    # 14. Final diagnostic
    # ----------------------------------------------------------------------
    print("\n14. FINAL DIAGNOSTIC")

    print(
        f"target holdout cyclotomic = {acc_c:.6f}"
    )
    print(
        f"target holdout control    = {acc_r:.6f}"
    )
    print(
        f"target delta              = {acc_c - acc_r:+.6f}"
    )

    print(
        f"LOO cyclotomic mean       = {mean(loo_c.values()):.6f}"
    )
    print(
        f"LOO control mean          = {mean(loo_r.values()):.6f}"
    )

    print(
        f"2D cyclotomic             = {acc2_c:.6f}"
    )
    print(
        f"2D control                = {acc2_r:.6f}"
    )

    print(
        f"permutation p-value       = {p_value:.6f}"
    )

    print("\nINTERPRETATION")
    print("-" * 78)

    if (
        acc_c > 0.55
        and acc_c > acc_r + 0.03
        and acc2_c > 0.55
        and p_value < 0.05
        and mean(loo_c.values()) > mean(loo_r.values())
    ):
        print(
            "POTENTIAL CYCLOTOMIC ORIENTATION SIGNAL\n"
            "The cross-modulus sign structure survives target and modulus\n"
            "holdout better than the control family.\n"
            "This would justify a subsequent CRT reconstruction experiment."
        )

    elif acc_c > 0.55 and acc_c > acc_r + 0.03:
        print(
            "TARGET-LEVEL EFFECT WITHOUT ROBUST CROSS-MODULUS GENERALIZATION\n"
            "There is a target-level separation, but the modulus-holdout or\n"
            "permutation tests do not support a universal orientation law."
        )

    else:
        print(
            "NO ROBUST CROSS-MODULUS ORIENTATION SIGNAL\n"
            "The Eisenstein orientation appears not to be predictably\n"
            "recoverable from the tested N-only features."
        )

    print("\nIMPORTANT")
    print("-" * 78)
    print(
        "The pair-sign construction itself is exact algebraically."
    )
    print(
        "The experiment asks whether N-only information predicts those signs."
    )
    print(
        "A successful classifier alone is not a factorization algorithm."
    )
    print(
        "Only a result surviving unseen targets, unseen moduli, and the"
        " permutation null should motivate CRT reconstruction."
    )

    print("\n" + "=" * 78)
    print("EXPERIMENT 64 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

