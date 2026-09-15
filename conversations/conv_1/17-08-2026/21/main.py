#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 91
DIRECT SIGMA_1 / E2 ACCESSIBILITY TEST

CORE OBSERVATION
----------------
For n = p*q:

    sigma_1(n) = (p+1)(q+1)
               = n + s + 1

Therefore:

    s = sigma_1(n) - n - 1

and once s is known:

    p,q = roots of x^2 - s*x + n.

This experiment asks whether sigma_1(n) itself has a reusable,
low-complexity modular representation.

PAPER / QUASIMODULAR CONNECTION
--------------------------------
E2(q) = 1 - 24 * sum_{n>=1} sigma_1(n) q^n

so the coefficient of E2 directly contains sigma_1(n).

TESTS
-----
1. Exact sigma_1 sequence modulo ell.
2. Simple period search.
3. Small-state prediction:
       n mod ell^k -> sigma_1(n) mod ell
4. Strict unseen-target validation.
5. Matched synthetic residue null.
6. CRT reconstruction of s from predicted sigma_1 residues.
7. Direct factor recovery from reconstructed s.

The main question is:

    Can sigma_1(n) mod ell be obtained from a small N-only state?

This is different from Experiments 88-90:
we are now targeting the arithmetic function that would DIRECTLY
give the factor sum.

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
    5, 7, 11, 13, 17,
    19, 23, 31, 37, 61, 67
]

DEPTHS = [1, 2, 3]

SEQ_LIMIT = 50000

PERIOD_MAX = 2000

SYNTHETIC_SAMPLES = 100000

RNG_SEED = 91091

PRINT_TARGETS = 24


# ============================================================================
# DATA
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

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

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
        x
        for x in range(lo, hi + 1)
        if sieve[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    rng = random.Random(
        RNG_SEED
    )

    targets: List[Target] = []
    seen = set()

    while len(targets) < count:

        p = primes[
            rng.randrange(
                len(primes)
            )
        ]

        q = primes[
            rng.randrange(
                len(primes)
            )
        ]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add(
            (p, q)
        )

        targets.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return targets


# ============================================================================
# SIGMA_1
# ============================================================================

def sigma1_semiprime(
    t: Target,
) -> int:

    return (
        t.p + 1
    ) * (
        t.q + 1
    )


def sigma1_from_ns(
    n: int,
    s: int,
) -> int:

    return (
        n + s + 1
    )


# ============================================================================
# SEQUENCE GENERATION
# ============================================================================

def sigma1_sequence(
    limit: int,
) -> List[int]:

    """
    sigma1(n) for every 0 <= n <= limit.

    Sieve:
        sigma1(j) += d
        for every divisor d | j.
    """

    sigma = [
        0
    ] * (
        limit + 1
    )

    for d in range(
        1,
        limit + 1,
    ):

        for j in range(
            d,
            limit + 1,
            d,
        ):
            sigma[j] += d

    return sigma


# ============================================================================
# E2 COEFFICIENT
# ============================================================================

def e2_coefficient_mod(
    sigma1: int,
    ell: int,
) -> int:

    """
    q^n coefficient of E2:

        -24 sigma1(n)

    The constant term of E2 is irrelevant here.
    """

    return (
        -24 * sigma1
    ) % ell


# ============================================================================
# PERIOD TEST
# ============================================================================

def period_errors(
    sequence: Sequence[int],
    period: int,
    warmup: int,
) -> int:

    """
    Count:

        a[n] != a[n-period]

    for n >= warmup + period.
    """

    start = max(
        period,
        warmup,
    )

    errors = 0

    for n in range(
        start,
        len(sequence),
    ):

        if (
            sequence[n]
            != sequence[n - period]
        ):
            errors += 1

    return errors


def search_best_period(
    sequence: Sequence[int],
    max_period: int,
    warmup: int,
) -> Tuple[int, int, float]:

    best_period = None
    best_errors = None

    comparisons = max(
        1,
        len(sequence) - warmup - 1,
    )

    for period in range(
        1,
        min(
            max_period,
            len(sequence) - warmup,
        ) + 1,
    ):

        errors = period_errors(
            sequence,
            period,
            warmup,
        )

        if (
            best_errors is None
            or errors < best_errors
        ):
            best_errors = errors
            best_period = period

    return (
        int(best_period),
        int(best_errors),
        best_errors / comparisons,
    )


# ============================================================================
# STATE LOOKUP MODEL
# ============================================================================

def build_state_model(
    targets: Sequence[Target],
    ell: int,
    k: int,
) -> Dict[int, int]:

    buckets = defaultdict(
        Counter
    )

    modulus = ell ** k

    for t in targets:

        x = t.n % modulus

        y = (
            sigma1_semiprime(t)
            % ell
        )

        buckets[x][y] += 1

    return {
        x: counter.most_common(
            1
        )[0][0]
        for x, counter in buckets.items()
    }


def evaluate_state_model(
    train: Sequence[Target],
    test: Sequence[Target],
    ell: int,
    k: int,
) -> Dict[str, float]:

    model = build_state_model(
        train,
        ell,
        k,
    )

    modulus = ell ** k

    truth = []
    predictions = []

    for t in test:

        x = t.n % modulus

        if x not in model:
            continue

        truth.append(
            sigma1_semiprime(t)
            % ell
        )

        predictions.append(
            model[x]
        )

    if not truth:

        return {
            "states": len(model),
            "coverage": 0.0,
            "accuracy": float("nan"),
            "balanced": float("nan"),
            "baseline": float("nan"),
        }

    accuracy = sum(
        a == b
        for a, b in zip(
            truth,
            predictions,
        )
    ) / len(truth)

    labels = sorted(
        set(truth)
    )

    recalls = []

    for label in labels:

        idx = [
            i
            for i, y in enumerate(truth)
            if y == label
        ]

        if not idx:
            continue

        recalls.append(
            sum(
                predictions[i]
                == label
                for i in idx
            )
            / len(idx)
        )

    balanced = (
        statistics.fmean(
            recalls
        )
        if recalls
        else float("nan")
    )

    baseline = (
        max(
            Counter(truth).values()
        )
        / len(truth)
    )

    return {
        "states": len(model),
        "coverage": (
            len(truth)
            / len(test)
        ),
        "accuracy": accuracy,
        "balanced": balanced,
        "baseline": baseline,
    }


# ============================================================================
# EMPIRICAL NULL
# ============================================================================

def synthetic_prime_residue_null(
    targets: Sequence[Target],
    ell: int,
    samples: int,
) -> List[Tuple[int, int, int]]:

    """
    Preserve the empirical marginal distribution of p mod ell and q mod ell.

    Return:
        (n residue, sigma1 residue, s residue)

    This is the finite-field matched null.
    """

    counts = Counter()

    for t in targets:
        counts[t.p % ell] += 1
        counts[t.q % ell] += 1

    residues = sorted(
        counts.keys()
    )

    weights = [
        counts[r]
        for r in residues
    ]

    rng = random.Random(
        RNG_SEED + ell
    )

    result = []

    for _ in range(samples):

        p = rng.choices(
            residues,
            weights=weights,
            k=1,
        )[0]

        q = rng.choices(
            residues,
            weights=weights,
            k=1,
        )[0]

        n = (
            p * q
        ) % ell

        sigma1 = (
            (p + 1)
            * (q + 1)
        ) % ell

        s = (
            p + q
        ) % ell

        result.append(
            (
                n,
                sigma1,
                s,
            )
        )

    return result


# ============================================================================
# DIRECT RECOVERY FROM SIGMA1
# ============================================================================

def recover_s_from_sigma1(
    n: int,
    sigma1_value: int,
) -> int:

    """
    Exact identity:

        sigma1(n) = n + s + 1

    Therefore:

        s = sigma1(n) - n - 1.
    """

    return (
        sigma1_value
        - n
        - 1
    )


def factor_from_ns(
    n: int,
    s: int,
) -> Tuple[int, int] | None:

    discriminant = (
        s * s
        - 4 * n
    )

    if discriminant < 0:
        return None

    d = math.isqrt(
        discriminant
    )

    if d * d != discriminant:
        return None

    if (
        (s + d) % 2
        != 0
    ):
        return None

    p = (
        s + d
    ) // 2

    q = (
        s - d
    ) // 2

    if p * q != n:
        return None

    return (
        min(p, q),
        max(p, q),
    )


# ============================================================================
# MODULAR CRT
# ============================================================================

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:

    inv = pow(
        m1,
        -1,
        m2,
    )

    k = (
        (a2 - a1)
        * inv
    ) % m2

    x = (
        a1
        + m1 * k
    )

    return (
        x % (m1 * m2),
        m1 * m2,
    )


def reconstruct_residue(
    residues: Sequence[Tuple[int, int]],
) -> Tuple[int, int]:

    x = 0
    m = 1

    for r, mod in residues:

        x, m = crt_pair(
            x,
            m,
            r,
            mod,
        )

    return (
        x,
        m,
    )


# ============================================================================
# MODULAR SIGMA1 STATE EXPERIMENT
# ============================================================================

def modular_sigma1_recovery(
    train: Sequence[Target],
    test: Sequence[Target],
    moduli: Sequence[int],
    k: int,
) -> Dict[str, float]:

    models = {}

    for ell in moduli:

        models[ell] = build_state_model(
            train,
            ell,
            k,
        )

    recovered = 0
    unique_s = 0
    factor_hits = 0
    usable = 0

    for t in test:

        residues = []
        possible = True

        for ell in moduli:

            state = (
                t.n
                % (
                    ell ** k
                )
            )

            model = models[ell]

            if state not in model:
                possible = False
                break

            sigma_residue = (
                model[state]
            ) % ell

            residues.append(
                (
                    sigma_residue,
                    ell,
                )
            )

        if not possible:
            continue

        usable += 1

        sigma_mod_M, M = reconstruct_residue(
            residues
        )

        # Need sigma1(n) itself to reconstruct s exactly.
        # Only accept if CRT modulus covers the possible sigma1 range.
        #
        # sigma1(n) = n+s+1 and s <= 2*sqrt(n) for positive p,q.
        upper = (
            t.n
            + 2 * math.isqrt(t.n)
            + 1
        )

        if M <= upper:
            continue

        # Candidate sigma1 is sigma_mod_M itself, since it is now
        # the unique representative in [0,M).
        sigma_candidate = sigma_mod_M

        if not (
            0
            <= sigma_candidate
            <= upper
        ):
            continue

        s_candidate = (
            sigma_candidate
            - t.n
            - 1
        )

        recovered += (
            s_candidate
            == t.s
        )

        factors = factor_from_ns(
            t.n,
            s_candidate,
        )

        if factors == (
            min(t.p, t.q),
            max(t.p, t.q),
        ):
            factor_hits += 1

    return {
        "usable": usable,
        "s_recovery": (
            recovered
            / len(test)
        ),
        "factor_recovery": (
            factor_hits
            / len(test)
        ),
    }


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 91")
    print("DIRECT SIGMA_1 / E2 ACCESSIBILITY TEST")
    print("SIGMA_1(N) = N + S + 1")
    print("E2 COEFFICIENT = -24 SIGMA_1(N)")
    print("LOCAL MODULAR ACCESS + CRT")
    print("STRICT TARGET HOLDOUT")
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
        f"{time.perf_counter()-t0:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    train = targets[
        :TRAIN_TARGETS
    ]

    test = targets[
        TRAIN_TARGETS:
    ]

    print(
        f"total targets = "
        f"{len(targets)}"
    )

    for i, t in enumerate(
        targets[:PRINT_TARGETS],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s} "
            f"sigma1={sigma1_semiprime(t)}"
        )

    if NUM_TARGETS > PRINT_TARGETS:
        print(
            "... remaining targets omitted"
        )

    print("\n2. TARGET HOLDOUT")
    print("-" * 78)

    print(
        f"training targets = "
        f"{len(train)}"
    )

    print(
        f"test targets = "
        f"{len(test)}"
    )

    # ------------------------------------------------------------------
    # Exact identity validation
    # ------------------------------------------------------------------

    print("\n3. SIGMA1 IDENTITY VALIDATION")
    print("-" * 78)

    failures = 0

    for i, t in enumerate(
        targets,
        1,
    ):

        sigma_direct = (
            sigma1_semiprime(t)
        )

        sigma_ns = (
            sigma1_from_ns(
                t.n,
                t.s,
            )
        )

        ok = (
            sigma_direct
            == sigma_ns
        )

        if not ok:
            failures += 1

        if i <= PRINT_TARGETS:
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
            "sigma1 identity validation failed"
        )

    print("status = PASS")

    # ------------------------------------------------------------------
    # E2 relationship
    # ------------------------------------------------------------------

    print("\n4. E2 COEFFICIENT CHECK")
    print("-" * 78)

    for ell in MODULI[:7]:

        mismatches = 0

        for t in test:

            sigma = (
                sigma1_semiprime(t)
            )

            e2 = e2_coefficient_mod(
                sigma,
                ell,
            )

            if (
                e2
                != (
                    -24 * sigma
                ) % ell
            ):
                mismatches += 1

        print(
            f"ell={ell:3d} "
            f"E2 coefficient mismatches="
            f"{mismatches}/{len(test)}"
        )

    # ------------------------------------------------------------------
    # Sequence / period tests
    # ------------------------------------------------------------------

    print("\n5. SIGMA1 / E2 SEQUENCE PERIOD SEARCH")
    print("-" * 78)

    seq_start = time.perf_counter()

    sigma_seq = sigma1_sequence(
        SEQ_LIMIT
    )

    print(
        f"sigma1 sequence limit = "
        f"{SEQ_LIMIT}"
    )

    print(
        f"sequence generation = "
        f"{time.perf_counter()-seq_start:.6f}s"
    )

    for ell in MODULI[:7]:

        seq_mod = [
            sigma_seq[n] % ell
            for n in range(
                1,
                SEQ_LIMIT + 1,
            )
        ]

        period, errors, rate = (
            search_best_period(
                seq_mod,
                PERIOD_MAX,
                500,
            )
        )

        print(
            f"ell={ell:3d} "
            f"best_period={period:5d} "
            f"errors={errors:6d} "
            f"error_rate={rate:.6f}"
        )

    # ------------------------------------------------------------------
    # Local state prediction
    # ------------------------------------------------------------------

    print("\n6. LOCAL SIGMA1 ACCESSIBILITY")
    print("-" * 78)

    print(
        "ell | k | states | coverage | "
        "accuracy | balanced | baseline"
    )

    local_results = []

    for ell in MODULI:

        for k in DEPTHS:

            result = evaluate_state_model(
                train,
                test,
                ell,
                k,
            )

            local_results.append(
                (
                    ell,
                    k,
                    result,
                )
            )

            print(
                f"{ell:3d} | "
                f"{k:1d} | "
                f"{result['states']:6d} | "
                f"{result['coverage']:.6f} | "
                f"{result['accuracy']:.6f} | "
                f"{result['balanced']:.6f} | "
                f"{result['baseline']:.6f}"
            )

    # ------------------------------------------------------------------
    # Matched finite-field null
    # ------------------------------------------------------------------

    print("\n7. EMPIRICAL PRIME-RESIDUE NULL")
    print("-" * 78)

    for ell in MODULI:

        synthetic = (
            synthetic_prime_residue_null(
                targets,
                ell,
                SYNTHETIC_SAMPLES,
            )
        )

        table = defaultdict(
            Counter
        )

        for n, sigma, _s in synthetic:
            table[n][sigma] += 1

        correct = 0
        total = 0

        for n, counter in table.items():

            if not counter:
                continue

            best = (
                counter.most_common(
                    1
                )[0][0]
            )

            correct += max(
                counter.values()
            )

            total += sum(
                counter.values()
            )

        synthetic_acc = (
            correct / total
            if total
            else float("nan")
        )

        print(
            f"ell={ell:3d} "
            f"synthetic finite-field "
            f"Bayes accuracy={synthetic_acc:.6f}"
        )

    # ------------------------------------------------------------------
    # CRT recovery
    # ------------------------------------------------------------------

    print("\n8. CRT SIGMA1 -> S RECONSTRUCTION")
    print("-" * 78)

    families = [
        (
            "C3",
            [7, 13, 19],
        ),
        (
            "C5",
            [7, 13, 19, 31, 37],
        ),
        (
            "C7",
            [7, 13, 19, 31, 37, 61, 67],
        ),
        (
            "R3",
            [673, 4561, 4759],
        ),
        (
            "R5",
            [673, 4561, 4759, 6211, 7879],
        ),
    ]

    for name, moduli in families:

        result = modular_sigma1_recovery(
            train,
            test,
            moduli,
            1,
        )

        print(
            f"{name:3s} "
            f"usable={result['usable']:3d}/"
            f"{len(test)} "
            f"s_recovery={result['s_recovery']:.6f} "
            f"factor_recovery={result['factor_recovery']:.6f}"
        )

    # ------------------------------------------------------------------
    # Direct comparison
    # ------------------------------------------------------------------

    print("\n9. DIRECT COMPARISON")
    print("-" * 78)

    print(
        "The exact identity says:"
    )

    print(
        "    sigma1(N) = N + s + 1"
    )

    print(
        "Therefore any successful exact sigma1 reconstruction "
        "immediately gives:"
    )

    print(
        "    s = sigma1(N) - N - 1"
    )

    print(
        "and then:"
    )

    print(
        "    p,q = roots of x^2 - s*x + N."
    )

    print()
    print(
        "The key distinction is:"
    )

    print(
        "    sigma1 accessibility"
    )

    print(
        "versus"
    )

    print(
        "    H6/H8 resolvent accessibility."
    )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Experiment 91 changes the target from a derived tower ratio"
    )

    print(
        "to the fundamental arithmetic function sigma_1(N)."
    )

    print()
    print(
        "This is the critical identity:"
    )

    print(
        "    sigma_1(N) = N + p + q + 1."
    )

    print()
    print(
        "Therefore:"
    )

    print(
        "    sigma_1(N)  ->  s  ->  p,q."
    )

    print()
    print(
        "A convincing positive result requires:"
    )

    print(
        "  1. high OOS sigma1 prediction;"
    )

    print(
        "  2. survival for several ell;"
    )

    print(
        "  3. survival at k=1;"
    )

    print(
        "  4. CRT reconstruction of sigma1 beyond its "
        "possible range;"
    )

    print(
        "  5. actual held-out factor recovery."
    )

    print()
    print(
        "A negative result would be equally useful:"
    )

    print(
        "it would strongly suggest that the partition/quasimodular "
        "route does not provide a simple local N-only evaluation of "
        "sigma_1(N)."
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter()-start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 91 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

