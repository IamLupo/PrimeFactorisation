#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 88
LOCAL RESOLVENT ACCESSIBILITY TEST

TARGET:
    R(N) = H8 / H6
         = (s^2 - 3N) / 4

QUESTION:
    Is the H6/H8 resolvent locally determined by N alone?

TESTS:
    X = N mod ell^k
    Y = R(N) mod ell

Depths:
    k = 1,2,3

For each modulus:
    - exact train lookup
    - strict unseen-target test
    - majority baseline
    - permutation null
    - balanced accuracy
    - unresolved H6 == 0 cases
    - residual conditional entropy

IMPORTANT:
    The oracle target R(N) is computed from p,q.
    The predictor state uses ONLY N mod ell^k.

This is NOT factorization.
It tests whether the new algebraic resolvent has a reusable
N-only local representation.

NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 2000
TRAIN_TARGETS = 1500

P_MIN = 2_000_000
P_MAX = 4_200_000

MODULI = [
    7, 13, 19, 31, 37,
    61, 67, 79, 127,
]

CONTROL = [
    673, 4561, 4759,
    6211, 7879,
]

DEPTHS = [1, 2, 3]

RNG_SEED = 88088

PERMUTATIONS = 200

# Only print detailed rows for these moduli.
DETAIL_MODULI = {
    7, 13, 19, 31, 37, 61, 67
}


# ============================================================================
# TARGET
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    lo = max(2, lo)

    sieve = bytearray(
        b"\x01" * (hi + 1)
    )

    sieve[0:2] = b"\x00\x00"

    for p in range(
        2,
        math.isqrt(hi) + 1,
    ):
        if sieve[p]:
            start = p * p
            sieve[
                start:hi + 1:p
            ] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x for x in range(lo, hi + 1)
        if sieve[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    rng = random.Random(RNG_SEED)

    out: List[Target] = []
    seen = set()

    while len(out) < count:

        p = primes[
            rng.randrange(len(primes))
        ]

        q = primes[
            rng.randrange(len(primes))
        ]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)

        out.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return out


# ============================================================================
# PAPER RESOLVENT
# ============================================================================

def sigma_semiprime(
    n: int,
    s: int,
    j: int,
) -> int:

    R = [0] * (j + 1)

    R[0] = 2

    if j >= 1:
        R[1] = s

    for r in range(
        2,
        j + 1,
    ):
        R[r] = (
            s * R[r - 1]
            - n * R[r - 2]
        )

    return (
        1
        + R[j]
        + n ** j
    )


def H6(
    n: int,
    s: int,
) -> int:

    num = (
        (n * n - n + 1)
        * sigma_semiprime(n, s, 1)
        - sigma_semiprime(n, s, 3)
    )

    if num % 6 != 0:
        raise ArithmeticError(
            "H6 integrality failure"
        )

    return num // 6


def H8(
    n: int,
    s: int,
) -> int:

    num = (
        -n * n
        * sigma_semiprime(n, s, 1)
        + (n * n + 1)
        * sigma_semiprime(n, s, 3)
        - sigma_semiprime(n, s, 5)
    )

    if num % 24 != 0:
        raise ArithmeticError(
            "H8 integrality failure"
        )

    return num // 24


# ============================================================================
# RESOLVENT
# ============================================================================

def resolvent_exact(
    n: int,
    s: int,
) -> Tuple[int | None, int | None]:

    """
    Return:

        R = H8/H6
        T = s^2 - 3n

    R may be non-integral, while T is always integral.

    For modular work we use:
        T = s^2 - 3n
    because it avoids division by 4.
    """

    h6 = H6(n, s)
    h8 = H8(n, s)

    if h6 == 0:
        return None, None

    numerator = h8
    denominator = h6

    if numerator % denominator == 0:
        ratio = numerator // denominator
    else:
        ratio = None

    T = s * s - 3 * n

    return ratio, T


# ============================================================================
# MODULAR RESOLVENT
# ============================================================================

def resolvent_mod(
    n: int,
    s: int,
    ell: int,
) -> int | None:

    """
    Primary target:

        T = s^2 - 3n mod ell

    This is equivalent to the H6/H8 resolvent:

        T = 4 H8/H6.

    Using T avoids H6 inversion degeneracy and tests the underlying
    algebraic quantity directly.

    The target is still derived from the true factor sum s.
    """

    return (
        (s % ell) ** 2
        - 3 * (n % ell)
    ) % ell


# ============================================================================
# STATE
# ============================================================================

def state_value(
    n: int,
    ell: int,
    k: int,
) -> int:

    modulus = ell ** k

    return n % modulus


# ============================================================================
# ENTROPY
# ============================================================================

def entropy_from_counts(
    counts: Counter,
) -> float:

    total = sum(
        counts.values()
    )

    if total == 0:
        return 0.0

    h = 0.0

    for c in counts.values():

        p = c / total

        if p > 0:
            h -= p * math.log2(p)

    return h


def conditional_entropy(
    pairs: Sequence[Tuple[int, int]],
) -> float:

    buckets = defaultdict(
        Counter
    )

    total = len(pairs)

    if total == 0:
        return 0.0

    for x, y in pairs:
        buckets[x][y] += 1

    h = 0.0

    for bucket in buckets.values():

        bucket_total = sum(
            bucket.values()
        )

        weight = (
            bucket_total / total
        )

        h += (
            weight
            * entropy_from_counts(bucket)
        )

    return h


# ============================================================================
# MAJORITY LOOKUP MODEL
# ============================================================================

def build_lookup(
    X: Sequence[int],
    Y: Sequence[int],
) -> Dict[int, int]:

    buckets = defaultdict(
        Counter
    )

    for x, y in zip(X, Y):
        buckets[x][y] += 1

    model = {}

    for x, counter in buckets.items():

        model[x] = counter.most_common(
            1
        )[0][0]

    return model


def accuracy(
    y_true: Sequence[int],
    y_pred: Sequence[int],
) -> float:

    if not y_true:
        return float("nan")

    return sum(
        a == b
        for a, b in zip(
            y_true,
            y_pred,
        )
    ) / len(y_true)


def balanced_accuracy(
    y_true: Sequence[int],
    y_pred: Sequence[int],
) -> float:

    labels = sorted(
        set(y_true)
    )

    recalls = []

    for label in labels:

        idx = [
            i for i, y in enumerate(
                y_true
            )
            if y == label
        ]

        if not idx:
            continue

        recalls.append(
            sum(
                y_pred[i] == label
                for i in idx
            )
            / len(idx)
        )

    if not recalls:
        return float("nan")

    return statistics.fmean(
        recalls
    )


# ============================================================================
# PERMUTATION NULL
# ============================================================================

def permutation_pvalue(
    train_x: Sequence[int],
    train_y: Sequence[int],
    observed_acc: float,
    permutations: int,
    rng: random.Random,
) -> float:

    if not train_y:
        return float("nan")

    exceed = 0

    model = None

    shuffled = list(train_y)

    for _ in range(permutations):

        rng.shuffle(shuffled)

        model = build_lookup(
            train_x,
            shuffled,
        )

        pred = [
            model.get(
                x,
                None,
            )
            for x in train_x
        ]

        # Unseen state inside training is impossible,
        # but keep the check.
        valid = [
            i
            for i, p in enumerate(pred)
            if p is not None
        ]

        if not valid:
            continue

        acc = sum(
            pred[i] == shuffled[i]
            for i in valid
        ) / len(valid)

        if acc >= observed_acc:
            exceed += 1

    return (
        exceed + 1
    ) / (
        permutations + 1
    )


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def run_local_test(
    targets: Sequence[Target],
    ell: int,
    k: int,
    within_nonzero: bool,
    rng: random.Random,
) -> Dict[str, float]:

    train = targets[
        :TRAIN_TARGETS
    ]

    test = targets[
        TRAIN_TARGETS:
    ]

    train_x = []
    train_y = []

    test_x = []
    test_y = []

    degenerate_train = 0
    degenerate_test = 0

    for t in train:

        h6 = H6(
            t.n,
            t.s,
        )

        if h6 == 0:
            degenerate_train += 1
            continue

        x = state_value(
            t.n,
            ell,
            k,
        )

        y = resolvent_mod(
            t.n,
            t.s,
            ell,
        )

        if y is None:
            continue

        # Optional alternative:
        # target only nonzero resolvent residues.
        if within_nonzero and y == 0:
            continue

        train_x.append(x)
        train_y.append(y)

    for t in test:

        h6 = H6(
            t.n,
            t.s,
        )

        if h6 == 0:
            degenerate_test += 1
            continue

        x = state_value(
            t.n,
            ell,
            k,
        )

        y = resolvent_mod(
            t.n,
            t.s,
            ell,
        )

        if y is None:
            continue

        if within_nonzero and y == 0:
            continue

        test_x.append(x)
        test_y.append(y)

    if not train_x:
        return {
            "train": 0,
            "test": 0,
            "coverage": 0.0,
            "accuracy": float("nan"),
            "balanced": float("nan"),
            "baseline": float("nan"),
            "cond_h": float("nan"),
            "p": float("nan"),
            "deg_train": degenerate_train,
            "deg_test": degenerate_test,
        }

    model = build_lookup(
        train_x,
        train_y,
    )

    pred = [
        model.get(
            x,
            None,
        )
        for x in test_x
    ]

    covered = [
        i
        for i, p in enumerate(pred)
        if p is not None
    ]

    if not covered:
        return {
            "train": len(train_x),
            "test": len(test_x),
            "coverage": 0.0,
            "accuracy": float("nan"),
            "balanced": float("nan"),
            "baseline": max(
                Counter(test_y).values()
            ) / len(test_y),
            "cond_h": conditional_entropy(
                list(zip(train_x, train_y))
            ),
            "p": float("nan"),
            "deg_train": degenerate_train,
            "deg_test": degenerate_test,
        }

    y_cov = [
        test_y[i]
        for i in covered
    ]

    p_cov = [
        pred[i]
        for i in covered
    ]

    acc = accuracy(
        y_cov,
        p_cov,
    )

    bal = balanced_accuracy(
        y_cov,
        p_cov,
    )

    baseline = (
        max(
            Counter(y_cov).values()
        )
        / len(y_cov)
    )

    pval = permutation_pvalue(
        train_x,
        train_y,
        accuracy(
            train_y,
            [
                model[x]
                for x in train_x
            ],
        ),
        PERMUTATIONS,
        rng,
    )

    return {
        "train": len(train_x),
        "test": len(test_x),
        "states": len(set(train_x)),
        "coverage": len(covered)
        / len(test_x),
        "accuracy": acc,
        "balanced": bal,
        "baseline": baseline,
        "cond_h": conditional_entropy(
            list(
                zip(
                    train_x,
                    train_y,
                )
            )
        ),
        "p": pval,
        "deg_train": degenerate_train,
        "deg_test": degenerate_test,
    }


# ============================================================================
# DIRECT N-ONLY BASELINES
# ============================================================================

def baseline_mod_formula(
    n: int,
    ell: int,
    mode: str,
) -> int:

    if mode == "n":
        return n % ell

    if mode == "n2":
        return (n * n) % ell

    if mode == "n2+n+1":
        return (
            n * n + n + 1
        ) % ell

    if mode == "n2+4n+1":
        return (
            n * n + 4*n + 1
        ) % ell

    raise ValueError(mode)


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 88")
    print("LOCAL RESOLVENT ACCESSIBILITY TEST")
    print("TARGET = H8/H6 = (s^2 - 3n)/4")
    print("STATE = n mod ell^k")
    print("STRICT UNSEEN-TARGET VALIDATION")
    print("PERMUTATION NULL")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)

    print(
        f"prime population = "
        f"{len(primes)}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    print(
        f"total targets = "
        f"{len(targets)}"
    )

    for i, t in enumerate(
        targets[:24],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if NUM_TARGETS > 24:
        print(
            "... remaining generated targets omitted"
        )

    # ------------------------------------------------------------------
    # Exact identity
    # ------------------------------------------------------------------

    print("\n2. RESOLVENT VALIDATION")
    print("-" * 78)

    failures = 0

    for i, t in enumerate(
        targets,
        1,
    ):

        h6 = H6(
            t.n,
            t.s,
        )

        h8 = H8(
            t.n,
            t.s,
        )

        if h6 == 0:
            failures += 1
            ok = False
        else:
            lhs = (
                4 * h8
            )

            rhs = (
                h6
                * (
                    t.s * t.s
                    - 3 * t.n
                )
            )

            ok = (
                lhs == rhs
            )

            if not ok:
                failures += 1

        if i <= 24:
            print(
                f"target {i:3d}: "
                f"identity_ok={ok}"
            )

    print(
        f"identity failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "H6/H8 resolvent validation failed"
        )

    print("status = PASS")

    # ------------------------------------------------------------------
    # Run local information tests
    # ------------------------------------------------------------------

    print("\n3. LOCAL RESOLVENT TEST")
    print("-" * 78)

    rng = random.Random(
        RNG_SEED + 1
    )

    all_results = []

    for ell in MODULI:

        if ell not in DETAIL_MODULI:
            continue

        print("\n")
        print(
            f"ell = {ell}"
        )
        print("-" * 78)

        for k in DEPTHS:

            result = run_local_test(
                targets,
                ell,
                k,
                False,
                rng,
            )

            all_results.append(
                (
                    ell,
                    k,
                    result,
                )
            )

            print(
                f"k={k} "
                f"train={result.get('train',0):4d} "
                f"test={result.get('test',0):4d} "
                f"states={result.get('states',0):4d} "
                f"coverage={result['coverage']:.4f} "
                f"acc={result['accuracy']:.6f} "
                f"bal={result['balanced']:.6f} "
                f"baseline={result['baseline']:.6f} "
                f"H(Y|X)={result['cond_h']:.6f} "
                f"p={result['p']:.6f}"
            )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("4. CROSS-DEPTH SUMMARY")
    print("=" * 78)

    print(
        "ell | k | coverage | accuracy | balanced | baseline | H(Y|X)"
    )

    for ell, k, r in all_results:

        print(
            f"{ell:3d} | "
            f"{k:1d} | "
            f"{r['coverage']:.6f} | "
            f"{r['accuracy']:.6f} | "
            f"{r['balanced']:.6f} | "
            f"{r['baseline']:.6f} | "
            f"{r['cond_h']:.6f}"
        )

    # ------------------------------------------------------------------
    # N-only state baselines
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("5. DIRECT N-ONLY MODULAR BASELINES")
    print("=" * 78)

    for ell in [
        7, 13, 19, 31, 37, 61, 67
    ]:

        train = targets[
            :TRAIN_TARGETS
        ]

        test = targets[
            TRAIN_TARGETS:
        ]

        for mode in [
            "n",
            "n2",
            "n2+n+1",
            "n2+4n+1",
        ]:

            buckets = defaultdict(
                Counter
            )

            for t in train:

                y = resolvent_mod(
                    t.n,
                    t.s,
                    ell,
                )

                x = baseline_mod_formula(
                    t.n,
                    ell,
                    mode,
                )

                buckets[x][y] += 1

            model = {
                x: counter.most_common(1)[0][0]
                for x, counter in buckets.items()
            }

            ys = []
            ps = []

            for t in test:

                y = resolvent_mod(
                    t.n,
                    t.s,
                    ell,
                )

                x = baseline_mod_formula(
                    t.n,
                    ell,
                    mode,
                )

                if x not in model:
                    continue

                ys.append(y)
                ps.append(model[x])

            if ys:

                print(
                    f"ell={ell:3d} "
                    f"{mode:10s} "
                    f"acc={accuracy(ys,ps):.6f} "
                    f"bal={balanced_accuracy(ys,ps):.6f}"
                )

    # ------------------------------------------------------------------
    # Nonzero-only test
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("6. NONZERO-RESOLVENT CONDITIONAL TEST")
    print("=" * 78)

    for ell in [
        7, 13, 19, 31, 37, 61, 67
    ]:

        r = run_local_test(
            targets,
            ell,
            1,
            True,
            rng,
        )

        print(
            f"ell={ell:3d} "
            f"k=1 "
            f"coverage={r['coverage']:.6f} "
            f"acc={r['accuracy']:.6f} "
            f"bal={r['balanced']:.6f} "
            f"baseline={r['baseline']:.6f} "
            f"p={r['p']:.6f}"
        )

    # ------------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The new target is not E1 or a hidden sign."
    )

    print(
        "It is the concrete algebraic quantity:"
    )

    print(
        "    T(N) = s^2 - 3N"
    )

    print(
        "equivalently:"
    )

    print(
        "    4 H8/H6 = T(N)"
    )

    print()
    print(
        "A successful result would show that T(N) is "
        "predictable from a small reusable state"
    )

    print(
        "    N mod ell^k"
    )

    print(
        "on unseen targets."
    )

    print()
    print(
        "THAT would be materially stronger than our earlier "
        "E1-character experiments because the target is now "
        "the exact quantity required to reconstruct s^2."
    )

    print()
    print(
        "Interpretation:"
    )

    print(
        "  OOS near baseline:"
    )
    print(
        "      no evidence for a local N-only resolvent."
    )

    print(
        "  OOS > baseline but poor LOO / permutation:"
    )
    print(
        "      likely residue collisions / overfitting."
    )

    print(
        "  OOS > baseline across k and multiple ell,"
    )
    print(
        "  with small permutation p-values:"
    )
    print(
        "      strong evidence for a reusable local "
        "N-only representation of the resolvent."
    )

    print()
    print(
        "If this fails, the next question should be whether "
        "the divisor-sum representation of H8/H6 admits a "
        "different arithmetic shortcut rather than another "
        "character classifier."
    )

    runtime = (
        time.perf_counter()
        - start
    )

    print(
        f"\ntotal runtime = "
        f"{runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 88 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

