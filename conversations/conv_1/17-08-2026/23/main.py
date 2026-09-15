#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 93
PARTIAL SIGMA_1 / EARLY-CERTIFICATE BARRIER

QUESTION
--------
Can sigma_1(N) reveal the hidden factor sum

    s = p + q

before we encounter a nontrivial divisor of N?

For semiprime N = p*q:

    sigma_1(N) = 1 + p + q + N
               = N + s + 1

Therefore, if a partial divisor sum

    S(B) = sum_{d|N, d<=B} d

contains enough information to force the value of s before B reaches
min(p,q), we may have a genuinely new computational shortcut.

If nothing useful happens until B reaches the first factor, then the
sigma_1 representation is information-rich but computationally equivalent
to divisor/factor search.

METHOD
------
For each target:

    N = p*q

we increase B through several checkpoints.

At each B we compute the exact N-only partial divisor sum:

    S(B) = sum_{d<=B, d|N} d

WITHOUT using p or q.

We then ask:

1. Have we encountered a nontrivial divisor?
2. How much of sigma_1(N) is known?
3. What interval for s remains possible?
4. Can s be uniquely reconstructed?
5. How many candidate s values remain?
6. What is the first B at which exact recovery occurs?

The experiment explicitly tracks the "barrier":

    first_nontrivial_divisor

versus

    first_sigma1_certificate.

If these coincide systematically, the divisor-sum route does not bypass
factor discovery.

NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import List, Optional, Tuple


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 16

P_MIN = 2_000_000
P_MAX = 4_200_000

RNG_SEED = 93093

# B checkpoints.
B_CHECKPOINTS = [
    10,
    100,
    1_000,
    10_000,
    100_000,
    250_000,
    500_000,
    1_000_000,
    1_500_000,
    2_000_000,
    2_500_000,
    3_000_000,
    3_500_000,
    4_000_000,
]


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
# PRIME SIEVE
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
    primes: List[int],
    count: int,
) -> List[Target]:

    rng = random.Random(
        RNG_SEED
    )

    targets = []
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
# PARTIAL DIVISOR SUM
# ============================================================================

def partial_sigma1(
    n: int,
    B: int,
) -> Tuple[int, Optional[int]]:
    """
    Compute

        S(B) = sum_{d|n, d<=B} d

    using only N and trial division.

    Returns:

        partial_sum
        first_nontrivial_divisor_seen

    The function intentionally does NOT use p or q.
    """

    limit = min(
        B,
        math.isqrt(n),
    )

    total = 1
    first_factor = None

    if limit < 2:
        return (
            total,
            None,
        )

    for d in range(
        2,
        limit + 1,
    ):

        if n % d != 0:
            continue

        total += d

        if first_factor is None:
            first_factor = d

        # We do not stop here because the purpose is to compute
        # the entire partial divisor sum through B.
        #
        # There can be only one nontrivial divisor below sqrt(N)
        # for a semiprime, but keeping the general logic explicit
        # makes the experiment reusable.

    return (
        total,
        first_factor,
    )


# ============================================================================
# PARTIAL SIGMA CERTIFICATE
# ============================================================================

def candidate_s_interval(
    n: int,
    partial_sum: int,
    B: int,
) -> Tuple[int, int]:

    """
    If no factor <= B has appeared, and the target construction assumes
    both prime factors lie in [P_MIN, P_MAX], derive a conservative
    interval for s.

    The identity is:

        sigma1(n) = n + s + 1

    and

        sigma1(n) = partial_sum + sum_of_unseen_divisors.

    Without knowing the unseen divisors, the partial sum alone gives
    only a broad interval.

    We use the guaranteed minimum contribution 0 and the maximal possible
    remaining contribution implied by the target scale.
    """

    # Very conservative lower bound:
    # sigma1 >= n + 1 + s,
    # but partial_sum <= sigma1.
    #
    # Thus:
    #
    # s >= partial_sum - n - 1.
    #
    lower = max(
        0,
        partial_sum - n - 1,
    )

    # From the configured prime range:
    upper = (
        2 * P_MAX
    )

    return (
        lower,
        upper,
    )


# ============================================================================
# EXACT RECOVERY CHECK
# ============================================================================

def recover_from_full_sigma(
    n: int,
    sigma1: int,
) -> Tuple[
    Optional[int],
    Optional[Tuple[int, int]],
]:

    s = (
        sigma1
        - n
        - 1
    )

    if s <= 0:
        return (
            None,
            None,
        )

    d2 = (
        s * s
        - 4 * n
    )

    if d2 < 0:
        return (
            None,
            None,
        )

    d = math.isqrt(
        d2
    )

    if d * d != d2:
        return (
            None,
            None,
        )

    if (
        (s + d) % 2
        != 0
    ):
        return (
            None,
            None,
        )

    p = (
        s + d
    ) // 2

    q = (
        s - d
    ) // 2

    if p * q != n:
        return (
            None,
            None,
        )

    return (
        s,
        (
            min(p, q),
            max(p, q),
        ),
    )


# ============================================================================
# ONE TARGET
# ============================================================================

def analyze_target(
    t: Target,
) -> dict:

    results = []

    first_factor_seen = None
    first_exact_sigma = None
    first_exact_s = None

    sigma_full = (
        t.n
        + t.s
        + 1
    )

    for B in B_CHECKPOINTS:

        if B >= math.isqrt(
            t.n
        ):

            # Full divisor range below sqrt(N).
            # For the semiprime this is enough to discover p.
            actual_B = math.isqrt(
                t.n
            )

        else:
            actual_B = B

        t0 = time.perf_counter()

        partial, found = partial_sigma1(
            t.n,
            actual_B,
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        if (
            found is not None
            and first_factor_seen is None
        ):
            first_factor_seen = (
                actual_B
            )

        lower, upper = (
            candidate_s_interval(
                t.n,
                partial,
                actual_B,
            )
        )

        exact_partial = (
            partial
            == sigma_full
        )

        if (
            exact_partial
            and first_exact_sigma is None
        ):
            first_exact_sigma = (
                actual_B
            )

        # A deliberately conservative certificate:
        # If partial sum is already the complete sigma1, exact s follows.
        s_candidate = None
        factors = None

        if exact_partial:

            (
                s_candidate,
                factors,
            ) = recover_from_full_sigma(
                t.n,
                partial,
            )

            if (
                s_candidate is not None
                and first_exact_s is None
            ):
                first_exact_s = (
                    actual_B
                )

        results.append(
            {
                "B": actual_B,
                "partial": partial,
                "fraction": (
                    partial
                    / sigma_full
                ),
                "first_factor": found,
                "lower_s": lower,
                "upper_s": upper,
                "candidate_width": (
                    max(
                        0,
                        upper - lower,
                    )
                ),
                "exact_sigma": exact_partial,
                "s_candidate": s_candidate,
                "factors": factors,
                "time": elapsed,
            }
        )

        # Avoid repeatedly doing essentially the same sqrt(N) scan.
        if (
            actual_B
            >= math.isqrt(t.n)
        ):
            break

    return {
        "target": t,
        "results": results,
        "first_factor_B": first_factor_seen,
        "first_exact_sigma_B": first_exact_sigma,
        "first_exact_s_B": first_exact_s,
        "sigma_full": sigma_full,
    }


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 93")
    print("PARTIAL SIGMA_1 / EARLY-CERTIFICATE BARRIER")
    print("N-ONLY DIVISOR-SUM ACCUMULATION")
    print("SIGMA_1 -> S -> FACTORS")
    print("NO LOCAL CHARACTER CLASSIFIER")
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

    print(
        f"total targets = "
        f"{len(targets)}"
    )

    for i, t in enumerate(
        targets[:20],
        1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if NUM_TARGETS > 20:
        print(
            "... remaining targets omitted"
        )

    # ------------------------------------------------------------------
    # Analyze targets
    # ------------------------------------------------------------------

    analyses = []

    print("\n2. PARTIAL SIGMA ANALYSIS")
    print("-" * 78)

    for i, t in enumerate(
        targets,
        1,
    ):

        analysis = analyze_target(
            t
        )

        analyses.append(
            analysis
        )

        if (
            i <= 8
            or i % 4 == 0
        ):

            first_factor = (
                analysis[
                    "first_factor_B"
                ]
            )

            first_sigma = (
                analysis[
                    "first_exact_sigma_B"
                ]
            )

            first_s = (
                analysis[
                    "first_exact_s_B"
                ]
            )

            print(
                f"processed {i:2d}/{NUM_TARGETS} "
                f"first_factor_B="
                f"{first_factor} "
                f"first_sigma_B="
                f"{first_sigma} "
                f"first_s_B="
                f"{first_s}"
            )

    # ------------------------------------------------------------------
    # Per-target checkpoint analysis
    # ------------------------------------------------------------------

    print("\n3. CHECKPOINT SUMMARY")
    print("-" * 78)

    checkpoint_rows = []

    for analysis in analyses:

        t = analysis["target"]

        for row in analysis["results"]:

            checkpoint_rows.append(
                (
                    t,
                    row,
                )
            )

    for B in B_CHECKPOINTS:

        rows = [
            row
            for t, row
            in checkpoint_rows
            if row["B"] == B
        ]

        if not rows:
            continue

        mean_fraction = statistics.fmean(
            r["fraction"]
            for r in rows
        )

        mean_width = statistics.fmean(
            r["candidate_width"]
            for r in rows
        )

        factor_seen = sum(
            r["first_factor"] is not None
            for r in rows
        )

        exact_sigma = sum(
            r["exact_sigma"]
            for r in rows
        )

        print(
            f"B={B:9d} "
            f"mean_sigma_fraction="
            f"{mean_fraction:.8f} "
            f"mean_s_width="
            f"{mean_width:.2f} "
            f"factor_seen="
            f"{factor_seen:2d}/{len(rows)} "
            f"exact_sigma="
            f"{exact_sigma:2d}/{len(rows)}"
        )

    # ------------------------------------------------------------------
    # Barrier comparison
    # ------------------------------------------------------------------

    print("\n4. FACTOR BARRIER VS SIGMA1 CERTIFICATE")
    print("-" * 78)

    factor_Bs = [
        a["first_factor_B"]
        for a in analyses
        if a["first_factor_B"] is not None
    ]

    sigma_Bs = [
        a["first_exact_sigma_B"]
        for a in analyses
        if a["first_exact_sigma_B"] is not None
    ]

    s_Bs = [
        a["first_exact_s_B"]
        for a in analyses
        if a["first_exact_s_B"] is not None
    ]

    print(
        f"targets with first nontrivial divisor found = "
        f"{len(factor_Bs)}/{NUM_TARGETS}"
    )

    print(
        f"targets with exact sigma1 certificate = "
        f"{len(sigma_Bs)}/{NUM_TARGETS}"
    )

    print(
        f"targets with exact s recovery = "
        f"{len(s_Bs)}/{NUM_TARGETS}"
    )

    if factor_Bs:

        print(
            f"mean first-factor B = "
            f"{statistics.fmean(factor_Bs):.2f}"
        )

        print(
            f"median first-factor B = "
            f"{statistics.median(factor_Bs):.2f}"
        )

    # ------------------------------------------------------------------
    # Partial information before first factor
    # ------------------------------------------------------------------

    print("\n5. PRE-FACTOR INFORMATION")
    print("-" * 78)

    fractions_before_factor = []

    for analysis in analyses:

        first_factor = (
            analysis[
                "first_factor_B"
            ]
        )

        if first_factor is None:
            continue

        candidate_rows = [
            r
            for r in analysis[
                "results"
            ]
            if r["B"] < first_factor
        ]

        if not candidate_rows:
            continue

        best = max(
            candidate_rows,
            key=lambda r:
                r["fraction"],
        )

        fractions_before_factor.append(
            best["fraction"]
        )

    if fractions_before_factor:

        print(
            "Largest partial sigma1 fraction "
            "observed before first factor:"
        )

        print(
            f"mean = "
            f"{statistics.fmean(fractions_before_factor):.8f}"
        )

        print(
            f"median = "
            f"{statistics.median(fractions_before_factor):.8f}"
        )

        print(
            f"max = "
            f"{max(fractions_before_factor):.8f}"
        )

    else:

        print(
            "No pre-factor checkpoint available."
        )

    # ------------------------------------------------------------------
    # Detailed examples
    # ------------------------------------------------------------------

    print("\n6. SAMPLE TARGET DETAILS")
    print("-" * 78)

    for i, analysis in enumerate(
        analyses[:5],
        1,
    ):

        t = analysis["target"]

        print(
            f"\nTarget {i}:"
        )

        print(
            f"  p={t.p}"
        )

        print(
            f"  q={t.q}"
        )

        print(
            f"  n={t.n}"
        )

        print(
            f"  s={t.s}"
        )

        print(
            f"  sqrt(n)="
            f"{math.isqrt(t.n)}"
        )

        print(
            f"  true sigma1="
            f"{analysis['sigma_full']}"
        )

        for row in analysis[
            "results"
        ]:

            print(
                f"  B={row['B']:9d} "
                f"partial={row['partial']} "
                f"fraction={row['fraction']:.8f} "
                f"first_factor="
                f"{row['first_factor']} "
                f"exact_sigma="
                f"{row['exact_sigma']}"
            )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The exact identity is:"
    )

    print(
        "    sigma1(N) = N + s + 1."
    )

    print()
    print(
        "This experiment asks whether partial divisor information "
        "can determine sigma1, and therefore s, BEFORE the first "
        "nontrivial divisor is encountered."
    )

    print()
    print(
        "POSITIVE OUTCOME:"
    )

    print(
        "    exact or nearly exact s certificate appears at B "
        "well below min(p,q)."
    )

    print()
    print(
        "NEGATIVE OUTCOME:"
    )

    print(
        "    sigma1 remains essentially N+1 until the first factor "
        "appears, and exact recovery starts only after that point."
    )

    print()
    print(
        "For the current semiprime construction, the latter result "
        "would demonstrate a structural barrier:"
    )

    print(
        "    the divisor-sum representation does not reveal the "
        "missing factor sum before divisor discovery."
    )

    print()
    print(
        "That would rule out another large class of sigma1-based "
        "shortcuts and justify moving away from direct coefficient "
        "evaluation."
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter()-start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 93 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

