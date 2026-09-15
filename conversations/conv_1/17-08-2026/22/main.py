#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 92
FOUR-SQUARE GLOBAL REPRESENTATION OF SIGMA_1(N)

KEY IDENTITY
------------
For odd N:

    r_4(N) = 8 sigma_1(N)

where r_4(N) counts ordered signed integer solutions of

    x1^2 + x2^2 + x3^2 + x4^2 = N.

For semiprime N = p q:

    sigma_1(N) = N + s + 1

therefore:

    s = r_4(N)/8 - N - 1.

Once s is known:

    p,q are roots of x^2 - s*x + N.

PURPOSE
-------
Earlier experiments attacked sigma_1 locally through N mod ell^k.

Those attempts failed.

This experiment asks whether sigma_1(N) has a completely different,
GLOBAL representation that can be computed without explicitly factoring N.

We compare several methods:

A. Direct 4-square enumeration
B. Meet-in-the-middle pair-square counting
C. Cached square multiplicities
D. Exact sigma_1 oracle for verification

IMPORTANT
---------
A positive identity result does NOT establish an efficient factorization
algorithm. The purpose is to measure whether the representation admits a
computational path that scales better than explicit divisor/factor search.

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
from typing import Dict, List, Tuple


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 40

P_MIN = 2_000_000
P_MAX = 4_200_000

RNG_SEED = 92092

# To keep the experiment runnable, use smaller benchmark targets as well.
BENCHMARK_BITS = [
    20,
    22,
    24,
]

PRINT_TARGETS = 20


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

    out = []
    seen = set()

    while len(out) < count:

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

        seen.add((p, q))

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
# EXACT ORACLE
# ============================================================================

def sigma1_oracle(
    t: Target,
) -> int:

    return (
        t.p + 1
    ) * (
        t.q + 1
    )


# ============================================================================
# FOUR-SQUARE MATHEMATICS
# ============================================================================

def square_count_map(
    limit: int,
) -> Dict[int, int]:

    """
    Count ordered signed square representations of one variable:

        x^2

    For each square value <= limit:

        0  -> 1 representation
        a^2 -> 2 representations for a != 0

    This allows fast convolution-style counting.
    """

    result = {}

    root = math.isqrt(
        limit
    )

    result[0] = 1

    for x in range(
        1,
        root + 1,
    ):

        result[x * x] = 2

    return result


def r4_meet_in_middle(
    n: int,
) -> int:

    """
    Exact meet-in-the-middle computation:

        r4(n)
          =
        sum_m r2(m) * r2(n-m)

    where r2(k) counts signed ordered representations

        x^2 + y^2 = k.

    This is O(n) in the naive implementation and therefore intended
    primarily as a scaling diagnostic.

    We deliberately do not use factorization information.
    """

    if n < 0:
        return 0

    root = math.isqrt(n)

    pair_counts = Counter()

    for x in range(
        -root,
        root + 1,
    ):

        remaining = (
            n - x * x
        )

        if remaining < 0:
            continue

        y = math.isqrt(
            remaining
        )

        if y * y == remaining:

            # Count ±y.
            if y == 0:
                pair_counts[
                    x * x + y * y
                ] += 1
            else:
                pair_counts[
                    x * x + y * y
                ] += 2

    total = 0

    # Direct 4-square enumeration through pair sums.
    #
    # This keeps the implementation exact and transparent.
    #
    # For demonstration targets, we enumerate possible x,y pairs
    # and count compatible z,w pairs.

    pair_counts.clear()

    for x in range(
        -root,
        root + 1,
    ):

        x2 = x * x

        if x2 > n:
            continue

        yroot = math.isqrt(
            n - x2
        )

        for y in range(
            -yroot,
            yroot + 1,
        ):

            s2 = (
                x2
                + y * y
            )

            if s2 > n:
                continue

            pair_counts[s2] += 1

    # Now count the complementary pair.
    #
    # Since pair_counts contains ALL ordered signed (x,y),
    # convolution gives r4.

    for a, count_a in pair_counts.items():

        b = n - a

        if b in pair_counts:

            total += (
                count_a
                * pair_counts[b]
            )

    return total


# ============================================================================
# OPTIMIZED TWO-PAIR ENUMERATION
# ============================================================================

def r4_pair_hash(
    n: int,
) -> int:

    """
    Hash-based two-pair method.

    Build all ordered signed pairs (x,y), grouped by x^2+y^2.
    Then match m and n-m.

    This is faster than explicitly enumerating four coordinates,
    but still scales roughly linearly with n.
    """

    root = math.isqrt(n)

    counts = Counter()

    # Instead of storing every pair, count them immediately.
    for x in range(
        -root,
        root + 1,
    ):

        xx = x * x

        if xx > n:
            continue

        yr = math.isqrt(
            n - xx
        )

        # We need all y in [-yr,yr].
        # Each y produces one ordered signed pair.
        for y in range(
            -yr,
            yr + 1,
        ):

            value = (
                xx
                + y * y
            )

            counts[value] += 1

    total = 0

    for a, c in counts.items():

        b = n - a

        if b in counts:

            total += (
                c
                * counts[b]
            )

    return total


# ============================================================================
# SMALLER BENCHMARK: DIRECT FOUR SQUARES
# ============================================================================

def r4_direct_small(
    n: int,
) -> int:

    """
    Brute-force exact four-square count for small n.

    Used only for validating the optimized implementation.
    """

    root = math.isqrt(n)

    total = 0

    for x in range(
        -root,
        root + 1,
    ):

        x2 = x * x

        if x2 > n:
            continue

        yroot = math.isqrt(
            n - x2
        )

        for y in range(
            -yroot,
            yroot + 1,
        ):

            xy = (
                x2
                + y * y
            )

            if xy > n:
                continue

            zroot = math.isqrt(
                n - xy
            )

            for z in range(
                -zroot,
                zroot + 1,
            ):

                xyz = (
                    xy
                    + z * z
                )

                if xyz > n:
                    continue

                w2 = (
                    n - xyz
                )

                w = math.isqrt(
                    w2
                )

                if w * w != w2:
                    continue

                if w == 0:
                    total += 1
                else:
                    total += 2

    return total


# ============================================================================
# RECONSTRUCTION
# ============================================================================

def reconstruct_from_r4(
    n: int,
    r4: int,
) -> Tuple[
    int | None,
    Tuple[int, int] | None,
]:

    if r4 % 8 != 0:
        return (
            None,
            None,
        )

    sigma1 = (
        r4 // 8
    )

    s = (
        sigma1
        - n
        - 1
    )

    factors = None

    if s > 0:

        d2 = (
            s * s
            - 4 * n
        )

        if d2 >= 0:

            d = math.isqrt(
                d2
            )

            if (
                d * d == d2
                and (s + d) % 2 == 0
            ):

                p = (
                    s + d
                ) // 2

                q = (
                    s - d
                ) // 2

                if p * q == n:
                    factors = (
                        min(p, q),
                        max(p, q),
                    )

    return (
        s,
        factors,
    )


# ============================================================================
# SMALL-N IDENTITY VALIDATION
# ============================================================================

def validate_four_square_identity() -> None:

    print("\n2. FOUR-SQUARE IDENTITY VALIDATION")
    print("-" * 78)

    failures = 0

    # Small odd semiprimes where brute force is feasible.
    tests = [
        3 * 5,
        5 * 7,
        7 * 11,
        11 * 13,
        13 * 17,
        17 * 19,
        19 * 23,
    ]

    for n in tests:

        expected = None

        # sigma1 for semiprime from factoring only for validation.
        # This section is deliberately a mathematical verification.
        #
        # Recover p,q from the known test construction.
        for p in range(
            2,
            math.isqrt(n) + 1,
        ):
            if (
                n % p == 0
                and n // p > p
            ):
                q = n // p

                expected = (
                    p + 1
                ) * (
                    q + 1
                )

                break

        if expected is None:
            failures += 1
            continue

        exact = r4_direct_small(
            n
        )

        ok = (
            exact
            == 8 * expected
        )

        print(
            f"n={n:6d} "
            f"r4={exact:6d} "
            f"8*sigma1={8*expected:6d} "
            f"ok={ok}"
        )

        if not ok:
            failures += 1

    print(
        f"identity failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Four-square identity validation failed"
        )

    print("status = PASS")


# ============================================================================
# PERFORMANCE BENCHMARK
# ============================================================================

def benchmark_method(
    method_name: str,
    method,
    ns: Sequence[int],
) -> None:

    print(
        f"\n{method_name}"
    )

    for n in ns:

        t0 = time.perf_counter()

        value = method(
            n
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        print(
            f"n={n:9d} "
            f"time={elapsed:.6f}s "
            f"r4={value}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 92")
    print("FOUR-SQUARE GLOBAL REPRESENTATION OF SIGMA_1(N)")
    print("r4(N) = 8 sigma_1(N) FOR ODD N")
    print("SIGMA1 -> S -> FACTORS")
    print("GLOBAL REPRESENTATION / COMPLEXITY TEST")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Identity validation
    # ------------------------------------------------------------------

    validate_four_square_identity()

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print("\n3. PRIME POPULATION")
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
    # Real targets
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
        targets[:PRINT_TARGETS],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if NUM_TARGETS > PRINT_TARGETS:
        print(
            "... remaining targets omitted"
        )

    # ------------------------------------------------------------------
    # Small benchmark scaling
    # ------------------------------------------------------------------

    print("\n4. GLOBAL R4 PERFORMANCE BENCHMARK")
    print("-" * 78)

    benchmark_ns = [
        (1 << b) - 1
        for b in BENCHMARK_BITS
    ]

    # Only benchmark the hash-pair method on sizes that remain realistic.
    for n in benchmark_ns:

        t0 = time.perf_counter()

        r4 = r4_pair_hash(
            n
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        print(
            f"n={n:9d} "
            f"time={elapsed:.6f}s "
            f"r4={r4}"
        )

    # ------------------------------------------------------------------
    # Real target feasibility sample
    # ------------------------------------------------------------------

    print("\n5. REAL TARGET FEASIBILITY")
    print("-" * 78)

    # Do NOT attempt all 40 targets.
    # Try only tiny representative samples and stop if the cost explodes.
    samples = [
        targets[0],
        targets[1],
    ]

    for i, t in enumerate(
        samples,
        1,
    ):

        print(
            f"\nTarget {i}: n={t.n}"
        )

        print(
            "Attempting global four-square evaluation..."
        )

        t0 = time.perf_counter()

        # Guard against accidental multi-minute runs.
        # We run the method directly but report runtime prominently.
        r4 = r4_pair_hash(
            t.n
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        sigma = (
            r4 // 8
            if r4 % 8 == 0
            else None
        )

        s_recovered = None
        factors = None

        if sigma is not None:

            (
                s_recovered,
                factors,
            ) = reconstruct_from_r4(
                t.n,
                r4,
            )

        print(
            f"elapsed={elapsed:.6f}s"
        )

        print(
            f"r4={r4}"
        )

        print(
            f"sigma1_candidate={sigma}"
        )

        print(
            f"recovered_s={s_recovered}"
        )

        print(
            f"true_s={t.s}"
        )

        print(
            f"recovered_factors={factors}"
        )

        print(
            f"true_factors="
            f"{(t.p, t.q)}"
        )

        if elapsed > 30:
            print(
                "WARNING: global four-square method "
                "already exceeds 30 seconds."
            )

    # ------------------------------------------------------------------
    # Exact oracle comparison
    # ------------------------------------------------------------------

    print("\n6. ORACLE RELATION")
    print("-" * 78)

    checks = 0

    for t in targets:

        sigma = sigma1_oracle(
            t
        )

        expected_s = (
            sigma
            - t.n
            - 1
        )

        factors = reconstruct_from_r4(
            t.n,
            8 * sigma,
        )

        ok = (
            expected_s == t.s
            and factors[1]
            == (
                min(t.p, t.q),
                max(t.p, t.q),
            )
        )

        checks += ok

    print(
        f"oracle sigma1 -> s -> factor "
        f"checks = {checks}/{len(targets)}"
    )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The identity itself is exact:"
    )

    print(
        "    r4(N) = 8 sigma1(N)"
    )

    print(
        "and for semiprime N=pq:"
    )

    print(
        "    sigma1(N) = N+s+1."
    )

    print()
    print(
        "Therefore a computationally efficient r4(N) evaluator "
        "would immediately solve the factor-sum problem."
    )

    print()
    print(
        "The decisive observation is computational:"
    )

    print(
        "If the four-square representation takes time scaling "
        "like O(N) or worse, it is not a practical shortcut."
    )

    print()
    print(
        "The experiment is therefore successful in either case:"
    )

    print(
        "  POSITIVE:"
    )

    print(
        "    find an unexpectedly fast global representation."
    )

    print(
        "  NEGATIVE:"
    )

    print(
        "    establish that this attractive modular-form identity "
        "does not translate into a computational shortcut."
    )

    print()
    print(
        "That is the question we need answered before continuing "
        "deeper into the quasimodular tower."
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter()-start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 92 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

