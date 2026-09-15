#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 90R
H6/H8 GCD / DIVISIBILITY FACTOR SIGNATURE SEARCH

PATCHED:
    No unvalidated H10/H12 formulas.
    Uses only the experimentally validated H6/H8 identities.

KNOWN EXACT IDENTITIES
----------------------
    6 H6 = A

    24 H8 = A * (s^2 - 3n)

    A = s * ((n+1)^2 - s^2)

Therefore:

    H8/H6 = (s^2 - 3n)/4

QUESTION
--------
Can a simple expression computed from N alone have a nontrivial gcd
with H6 or H8 that exposes p or q?

We search:

    gcd(F(N), H6)
    gcd(F(N), H8)

for a family of low-degree N-only polynomials F.

The decisive metric is unseen-target factor recovery.

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
from typing import Callable, List, Tuple


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 120
TRAIN_TARGETS = 90

P_MIN = 2_000_000
P_MAX = 4_200_000

RNG_SEED = 90090

COEFFS = range(-4, 5)

PRINT_TOP = 50


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
    primes: List[int],
    count: int,
) -> List[Target]:

    rng = random.Random(
        RNG_SEED
    )

    out: List[Target] = []
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
# POWER SUM
# ============================================================================

def power_sum(
    n: int,
    s: int,
    k: int,
) -> int:

    if k == 0:
        return 2

    if k == 1:
        return s

    r0 = 2
    r1 = s

    for _ in range(2, k + 1):
        r0, r1 = (
            r1,
            s * r1 - n * r0,
        )

    return r1


# ============================================================================
# DIVISOR SUMS FOR pq
# ============================================================================

def sigma_semiprime(
    n: int,
    s: int,
    k: int,
) -> int:

    return (
        1
        + power_sum(n, s, k)
        + n ** k
    )


# ============================================================================
# VALIDATED H6
# ============================================================================

def h6(
    n: int,
    s: int,
) -> int:

    """
    Exact validated formula:

        H6 = E1/6

        E1 =
            (n^2-n+1) sigma_1(n)
            - sigma_3(n)
    """

    numerator = (
        (n * n - n + 1)
        * sigma_semiprime(
            n,
            s,
            1,
        )
        - sigma_semiprime(
            n,
            s,
            3,
        )
    )

    if numerator % 6 != 0:
        raise ArithmeticError(
            "H6 integrality failure"
        )

    return numerator // 6


# ============================================================================
# VALIDATED H8
# ============================================================================

def h8(
    n: int,
    s: int,
) -> int:

    """
    Exact validated formula:

        24 H8 =
            -n^2 sigma_1(n)
            +(n^2+1) sigma_3(n)
            -sigma_5(n)
    """

    numerator = (
        -n * n
        * sigma_semiprime(
            n,
            s,
            1,
        )
        + (n * n + 1)
        * sigma_semiprime(
            n,
            s,
            3,
        )
        - sigma_semiprime(
            n,
            s,
            5,
        )
    )

    if numerator % 24 != 0:
        raise ArithmeticError(
            "H8 integrality failure"
        )

    return numerator // 24


# ============================================================================
# TOWER VALIDATION
# ============================================================================

def validate_target(
    t: Target,
) -> bool:

    H6 = h6(
        t.n,
        t.s,
    )

    H8 = h8(
        t.n,
        t.s,
    )

    A = (
        t.s
        * (
            (t.n + 1) ** 2
            - t.s * t.s
        )
    )

    if (
        6 * H6
        != A
    ):
        return False

    if (
        24 * H8
        != A
        * (
            t.s * t.s
            - 3 * t.n
        )
    ):
        return False

    # Equivalent resolvent check.
    if (
        4 * H8
        != H6
        * (
            t.s * t.s
            - 3 * t.n
        )
    ):
        return False

    return True


# ============================================================================
# N-ONLY PROBE
# ============================================================================

@dataclass(frozen=True)
class Probe:
    name: str
    fn: Callable[[int], int]


def make_probes() -> List[Probe]:

    probes = [
        Probe(
            "n",
            lambda n: n,
        ),
        Probe(
            "n-1",
            lambda n: n - 1,
        ),
        Probe(
            "n+1",
            lambda n: n + 1,
        ),
        Probe(
            "n^2-1",
            lambda n: n * n - 1,
        ),
        Probe(
            "n^2+n+1",
            lambda n: n * n + n + 1,
        ),
        Probe(
            "n^2-n+1",
            lambda n: n * n - n + 1,
        ),
        Probe(
            "n^2+4n+1",
            lambda n: n * n + 4 * n + 1,
        ),
        Probe(
            "n^3-1",
            lambda n: n ** 3 - 1,
        ),
        Probe(
            "n^3+1",
            lambda n: n ** 3 + 1,
        ),
    ]

    # Bounded cubic family.
    for a in COEFFS:
        for b in COEFFS:
            for c in COEFFS:
                for d in COEFFS:

                    if (
                        a == 0
                        and b == 0
                        and c == 0
                        and d == 0
                    ):
                        continue

                    complexity = (
                        abs(a)
                        + abs(b)
                        + abs(c)
                        + abs(d)
                    )

                    if complexity > 8:
                        continue

                    name = (
                        f"{a}n^3+"
                        f"{b}n^2+"
                        f"{c}n+"
                        f"{d}"
                    )

                    probes.append(
                        Probe(
                            name=name,
                            fn=lambda n,
                                a=a,
                                b=b,
                                c=c,
                                d=d:
                                (
                                    a * n ** 3
                                    + b * n ** 2
                                    + c * n
                                    + d
                                ),
                        )
                    )

    return probes


# ============================================================================
# GCD METRICS
# ============================================================================

def gcd_metrics(
    probe: Probe,
    targets: List[Target],
    h6_values: List[int],
    h8_values: List[int],
) -> Tuple[
    float,
    float,
    float,
    float,
]:

    """
    Returns:

        nontrivial_rate
        exact_factor_rate
        mean_max_gcd_digits
        mean_gcd_over_n

    A nontrivial gcd means:

        1 < g < n

    An exact factor hit means:

        g == p or g == q
    """

    nontrivial = 0
    exact = 0

    log_values = []
    fractions = []

    for t, H6, H8 in zip(
        targets,
        h6_values,
        h8_values,
    ):

        x = probe.fn(
            t.n
        )

        best_g = 1
        factor_hit = False

        for H in (
            H6,
            H8,
        ):

            g = math.gcd(
                abs(x),
                abs(H),
            )

            if g > best_g:
                best_g = g

            if (
                g > 1
                and g < t.n
            ):
                nontrivial += 1
                # Only count each target once.
                break

        for H in (
            H6,
            H8,
        ):

            g = math.gcd(
                abs(x),
                abs(H),
            )

            if (
                g == t.p
                or g == t.q
            ):
                factor_hit = True
                break

        if factor_hit:
            exact += 1

        log_values.append(
            math.log10(
                best_g + 1
            )
        )

        fractions.append(
            best_g / t.n
        )

    count = len(targets)

    return (
        nontrivial / count,
        exact / count,
        statistics.fmean(
            log_values
        ),
        statistics.fmean(
            fractions
        ),
    )


# ============================================================================
# BASELINE RANDOM PROBES
# ============================================================================

def random_probe_result(
    targets: List[Target],
    h6_values: List[int],
    h8_values: List[int],
    count: int,
) -> Tuple[float, float]:

    rng = random.Random(
        RNG_SEED + 12345
    )

    nontrivial = 0
    exact = 0

    for t, H6, H8 in zip(
        targets,
        h6_values,
        h8_values,
    ):

        # Random cubic coefficients.
        coeffs = [
            rng.randint(
                -8,
                8,
            )
            for _ in range(4)
        ]

        x = (
            coeffs[0] * t.n ** 3
            + coeffs[1] * t.n ** 2
            + coeffs[2] * t.n
            + coeffs[3]
        )

        hit = False

        for H in (
            H6,
            H8,
        ):

            g = math.gcd(
                abs(x),
                abs(H),
            )

            if (
                g > 1
                and g < t.n
            ):
                hit = True

            if (
                g == t.p
                or g == t.q
            ):
                exact += 1
                break

        if hit:
            nontrivial += 1

    return (
        nontrivial / len(targets),
        exact / len(targets),
    )


# ============================================================================
# DIRECT FACTOR DIVISIBILITY TABLE
# ============================================================================

def direct_divisibility(
    label: str,
    targets: List[Target],
    values: List[int],
) -> None:

    p_hits = sum(
        value % t.p == 0
        for t, value in zip(
            targets,
            values,
        )
    )

    q_hits = sum(
        value % t.q == 0
        for t, value in zip(
            targets,
            values,
        )
    )

    print(
        f"{label:8s} "
        f"divisible_by_p="
        f"{p_hits}/{len(targets)} "
        f"divisible_by_q="
        f"{q_hits}/{len(targets)}"
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 90R")
    print("PATCHED H6/H8 GCD / DIVISIBILITY FACTOR SIGNATURE SEARCH")
    print("NO UNVALIDATED H10/H12")
    print("STRICT TARGET HOLDOUT")
    print("N-ONLY PROBE FAMILY")
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
            "... remaining targets omitted"
        )

    train = targets[
        :TRAIN_TARGETS
    ]

    test = targets[
        TRAIN_TARGETS:
    ]

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
    # H6/H8
    # ------------------------------------------------------------------

    print("\n3. H6/H8 CONSTRUCTION")
    print("-" * 78)

    all_h6 = []
    all_h8 = []

    failures = 0

    for i, t in enumerate(
        targets,
        1,
    ):

        H6 = h6(
            t.n,
            t.s,
        )

        H8 = h8(
            t.n,
            t.s,
        )

        ok = validate_target(
            t
        )

        if not ok:
            failures += 1

        all_h6.append(H6)
        all_h8.append(H8)

        if i <= 24:
            print(
                f"target {i:3d}: "
                f"identity_ok={ok}"
            )

    if failures:
        raise RuntimeError(
            "H6/H8 validation failed"
        )

    print(
        "identity failures = 0"
    )

    print(
        "status = PASS"
    )

    train_h6 = all_h6[
        :TRAIN_TARGETS
    ]

    train_h8 = all_h8[
        :TRAIN_TARGETS
    ]

    test_h6 = all_h6[
        TRAIN_TARGETS:
    ]

    test_h8 = all_h8[
        TRAIN_TARGETS:
    ]

    # ------------------------------------------------------------------
    # Structural facts
    # ------------------------------------------------------------------

    print("\n4. STRUCTURAL IDENTITIES")
    print("-" * 78)

    print(
        "6H6 = s*((n+1)^2-s^2)"
    )

    print(
        "24H8 = "
        "6H6*(s^2-3n)"
    )

    print(
        "4H8/H6 = s^2-3n"
    )

    # ------------------------------------------------------------------
    # Direct gcds
    # ------------------------------------------------------------------

    print("\n5. DIRECT H6/H8 GCD STATISTICS")
    print("-" * 78)

    for label, ts, h6s, h8s in (
        (
            "TRAIN",
            train,
            train_h6,
            train_h8,
        ),
        (
            "TEST",
            test,
            test_h6,
            test_h8,
        ),
    ):

        print(
            f"\n{label}"
        )

        pairs = [
            (
                "gcd(H6,H8)",
                [
                    math.gcd(
                        a,
                        b,
                    )
                    for a, b in zip(
                        h6s,
                        h8s,
                    )
                ],
            )
        ]

        for name, values in pairs:

            print(
                f"{name:12s} "
                f"mean_log10="
                f"{statistics.fmean(math.log10(v+1) for v in values): .6f}"
            )

            factor_hits = sum(
                (
                    v == t.p
                    or v == t.q
                )
                for t, v in zip(
                    ts,
                    values,
                )
            )

            print(
                f"{name:12s} "
                f"exact_factor_gcd="
                f"{factor_hits}/{len(ts)}"
            )

    # ------------------------------------------------------------------
    # N-only probes
    # ------------------------------------------------------------------

    print("\n6. N-ONLY GCD PROBE SEARCH")
    print("-" * 78)

    probes = make_probes()

    print(
        f"candidate probes = "
        f"{len(probes)}"
    )

    rows = []

    for probe in probes:

        train_m = gcd_metrics(
            probe,
            train,
            train_h6,
            train_h8,
        )

        test_m = gcd_metrics(
            probe,
            test,
            test_h6,
            test_h8,
        )

        rows.append(
            (
                probe.name,
                train_m,
                test_m,
            )
        )

    rows.sort(
        key=lambda row: (
            row[2][1],
            row[2][0],
            row[2][2],
        ),
        reverse=True,
    )

    print(
        "\nTOP TEST PROBES"
    )

    print(
        "probe | "
        "train_nontriv | "
        "test_nontriv | "
        "train_exact | "
        "test_exact | "
        "test_log10(g)"
    )

    for name, train_m, test_m in (
        rows[:PRINT_TOP]
    ):

        print(
            f"{name:28s} "
            f"{train_m[0]:.4f} "
            f"{test_m[0]:.4f} "
            f"{train_m[1]:.4f} "
            f"{test_m[1]:.4f} "
            f"{test_m[2]:.4f}"
        )

    # ------------------------------------------------------------------
    # Explicit probes
    # ------------------------------------------------------------------

    print("\n7. EXPLICIT STRUCTURAL PROBES")
    print("-" * 78)

    for probe_name in [
        "n",
        "n-1",
        "n+1",
        "n^2-1",
        "n^2+n+1",
        "n^2-n+1",
        "n^2+4n+1",
        "n^3-1",
        "n^3+1",
    ]:

        for name, train_m, test_m in rows:

            if name != probe_name:
                continue

            print(
                f"{name:15s} "
                f"train_nontriv={train_m[0]:.4f} "
                f"test_nontriv={test_m[0]:.4f} "
                f"train_exact={train_m[1]:.4f} "
                f"test_exact={test_m[1]:.4f}"
            )

            break

    # ------------------------------------------------------------------
    # Random null
    # ------------------------------------------------------------------

    print("\n8. RANDOM CUBIC NULL")
    print("-" * 78)

    null_train = random_probe_result(
        train,
        train_h6,
        train_h8,
        1000,
    )

    null_test = random_probe_result(
        test,
        test_h6,
        test_h8,
        1000,
    )

    print(
        f"train nontrivial="
        f"{null_train[0]:.4f} "
        f"exact="
        f"{null_train[1]:.4f}"
    )

    print(
        f"test  nontrivial="
        f"{null_test[0]:.4f} "
        f"exact="
        f"{null_test[1]:.4f}"
    )

    # ------------------------------------------------------------------
    # Direct algebraic divisibility checks
    # ------------------------------------------------------------------

    print("\n9. DIRECT FACTOR DIVISIBILITY")
    print("-" * 78)

    direct_divisibility(
        "H6",
        train,
        train_h6,
    )

    direct_divisibility(
        "H8",
        train,
        train_h8,
    )

    direct_divisibility(
        "H6",
        test,
        test_h6,
    )

    direct_divisibility(
        "H8",
        test,
        test_h8,
    )

    # ------------------------------------------------------------------
    # Candidate suspicious results
    # ------------------------------------------------------------------

    print("\n10. SUSPICIOUS TEST PROBES")
    print("-" * 78)

    suspicious = []

    for name, train_m, test_m in rows:

        if (
            test_m[1] > 0
            or (
                test_m[0]
                > null_test[0] + 0.10
            )
        ):
            suspicious.append(
                (
                    name,
                    train_m,
                    test_m,
                )
            )

    if suspicious:

        for name, train_m, test_m in (
            suspicious[:30]
        ):
            print(
                f"{name:28s} "
                f"test_nontriv={test_m[0]:.4f} "
                f"test_exact={test_m[1]:.4f}"
            )

    else:

        print(
            "No probe showed a convincing "
            "held-out factor signature."
        )

    # ------------------------------------------------------------------
    # Final
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "This patched experiment deliberately stops at H6/H8."
    )

    print()
    print(
        "The central question is:"
    )

    print(
        "Can some low-degree N-only expression F(N) "
        "produce a nontrivial gcd with H6 or H8?"
    )

    print()
    print(
        "A genuinely interesting result would be:"
    )

    print(
        "    gcd(F(N), H6/H8) = p or q"
    )

    print(
        "on unseen targets, at a rate far above the "
        "random-polynomial null."
    )

    print()
    print(
        "If that does NOT happen, then the common factor"
    )

    print(
        "    A = s*((n+1)^2-s^2)"
    )

    print(
        "is not exposed by these simple N-only gcd probes."
    )

    print()
    print(
        "That would leave the real bottleneck:"
    )

    print(
        "    evaluate H6/H8 or their ratio from N "
        "without knowing s."
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 90R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()