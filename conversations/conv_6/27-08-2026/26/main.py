#!/usr/bin/env python3
"""
====================================================================================================
PERMUTATION NULL / MATCHED ANCHOR FACTOR-GEOMETRY EXPERIMENT
====================================================================================================

Purpose
-------
Test whether the factor geometry of the TRUE semiprime anchors is distinguishable
from randomly permuted factor-pair assignments on the exact same collision landscape.

For each of N fixed actual semiprimes:

    n_i = p_i * q_i

For every detected collision:

    N_t = n_i + 2*M*t = a*b

The collision event (t, N_t, a, b) is kept fixed.

NULL MODEL:
    randomly permute the original (p_i,q_i) pairs across the 300 anchors.

This preserves:
    * the exact n_i values
    * the exact collision events
    * the exact t distribution
    * the exact shifted semiprimes
    * the exact factorization geometry (a,b)

Only the association between anchor n_i and factor pair (p_i,q_i) changes.

Statistics tested:
    mean normalized dp
    mean normalized dq
    mean dSum
    mean dDiff
    corr(t, dSum)
    corr(t, dDiff)
    same-sign fraction
    opposite-sign fraction

The actual assignment is compared against 1000 permutation realizations.

====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
from dataclasses import dataclass
from collections import Counter, defaultdict


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435
STEP = 2 * M

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

T_MIN = -25
T_MAX = 25

TRIALS = 300
PERMUTATIONS = 1000

SEED = 1511464998


# ================================================================================================
# DATA STRUCTURES
# ================================================================================================

@dataclass(frozen=True)
class Event:
    trial: int
    n: int
    p: int
    q: int
    a: int
    b: int
    t: int


# ================================================================================================
# PRIME / SEMIPRIME GENERATION
# ================================================================================================

def sieve(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []

    composite = bytearray(limit + 1)
    primes = []

    for i in range(2, limit + 1):
        if composite[i]:
            continue

        primes.append(i)

        if i * i <= limit:
            composite[i * i : limit + 1 : i] = b"\x01" * (
                ((limit - i * i) // i) + 1
            )

    return primes


def generate_unique_semiprimes(
    count: int,
    lo: int,
    hi: int,
    rng: random.Random,
) -> list[tuple[int, int, int]]:
    """
    Generate count unique semiprimes p*q where both p,q are primes in [lo,hi].
    """

    primes = sieve(hi)
    primes = [p for p in primes if p >= lo]

    if len(primes) < 2:
        raise RuntimeError("Not enough primes in requested range.")

    seen: set[int] = set()
    result: list[tuple[int, int, int]] = []

    while len(result) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)
        result.append((p, q, n))

    return result


# ================================================================================================
# FACTORIZATION OF SHIFTED VALUES
# ================================================================================================

def factor_in_range(n: int, lo: int, hi: int) -> list[tuple[int, int]]:
    """
    Find factor pairs a*b=n with:

        lo <= a <= hi
        lo <= b <= hi
        a <= b
    """

    out: list[tuple[int, int]] = []

    a0 = max(lo, math.isqrt(n))

    # We only need divisors <= sqrt(n).
    start = lo

    for a in range(start, min(hi, math.isqrt(n)) + 1):
        if n % a != 0:
            continue

        b = n // a

        if b < a:
            continue

        if b > hi:
            continue

        out.append((a, b))

    return out


# ================================================================================================
# COLLISION LANDSCAPE
# ================================================================================================

def build_collision_landscape(
    anchors: list[tuple[int, int, int]]
) -> list[Event]:
    """
    Build the complete collision landscape.

    IMPORTANT:
    The collision landscape is tied only to n and t.
    """

    events: list[Event] = []

    for trial, (p, q, n) in enumerate(anchors, start=1):

        for t in range(T_MIN, T_MAX + 1):
            if t == 0:
                continue

            shifted = n + STEP * t

            if shifted <= 0:
                continue

            pairs = factor_in_range(
                shifted,
                FACTOR_MIN,
                FACTOR_MAX
            )

            for a, b in pairs:
                events.append(
                    Event(
                        trial=trial,
                        n=n,
                        p=p,
                        q=q,
                        a=a,
                        b=b,
                        t=t,
                    )
                )

    return events


# ================================================================================================
# CORRELATION
# ================================================================================================

def pearson(xs: list[float], ys: list[float]) -> float:
    """Pearson correlation without external dependencies."""

    if len(xs) != len(ys) or len(xs) < 2:
        return float("nan")

    mx = statistics.fmean(xs)
    my = statistics.fmean(ys)

    num = 0.0
    dx2 = 0.0
    dy2 = 0.0

    for x, y in zip(xs, ys):
        dx = x - mx
        dy = y - my

        num += dx * dy
        dx2 += dx * dx
        dy2 += dy * dy

    den = math.sqrt(dx2 * dy2)

    if den == 0:
        return float("nan")

    return num / den


# ================================================================================================
# STATISTICS
# ================================================================================================

STAT_NAMES = (
    "mean_ndp",
    "mean_ndq",
    "mean_dsum",
    "mean_ddiff",
    "corr_t_dsum",
    "corr_t_ddiff",
    "same_sign_fraction",
    "opposite_sign_fraction",
)


def calculate_statistics(
    events: list[Event],
    factor_assignment: dict[int, tuple[int, int]],
) -> dict[str, float]:

    if not events:
        raise RuntimeError("No collision events found.")

    ndp: list[float] = []
    ndq: list[float] = []

    dsum: list[float] = []
    ddiff: list[float] = []

    ts: list[float] = []

    same = 0
    opposite = 0

    for e in events:

        p, q = factor_assignment[e.trial]

        dp = e.a - p
        dq = e.b - q

        s = dp + dq
        d = dp - dq

        ndp.append(dp / 100_000.0)
        ndq.append(dq / 100_000.0)

        dsum.append(float(s))
        ddiff.append(float(d))

        ts.append(float(e.t))

        if dp == 0 or dq == 0:
            # Do not force zero into either sign category.
            continue

        if (dp > 0 and dq > 0) or (dp < 0 and dq < 0):
            same += 1
        else:
            opposite += 1

    classified = same + opposite

    return {
        "mean_ndp": statistics.fmean(ndp),
        "mean_ndq": statistics.fmean(ndq),
        "mean_dsum": statistics.fmean(dsum),
        "mean_ddiff": statistics.fmean(ddiff),
        "corr_t_dsum": pearson(ts, dsum),
        "corr_t_ddiff": pearson(ts, ddiff),
        "same_sign_fraction": (
            same / classified if classified else float("nan")
        ),
        "opposite_sign_fraction": (
            opposite / classified if classified else float("nan")
        ),
    }


# ================================================================================================
# PERMUTATION
# ================================================================================================

def permutation_assignment(
    anchors: list[tuple[int, int, int]],
    rng: random.Random,
) -> dict[int, tuple[int, int]]:

    pairs = [(p, q) for p, q, _ in anchors]

    rng.shuffle(pairs)

    return {
        trial: pairs[trial - 1]
        for trial in range(1, len(anchors) + 1)
    }


# ================================================================================================
# EMPIRICAL P VALUE
# ================================================================================================

def two_sided_empirical_p(
    actual: float,
    null_values: list[float],
) -> float:

    null_values = [
        x for x in null_values
        if math.isfinite(x)
    ]

    if not null_values or not math.isfinite(actual):
        return float("nan")

    null_mean = statistics.fmean(null_values)

    actual_distance = abs(actual - null_mean)

    extreme = sum(
        abs(x - null_mean) >= actual_distance
        for x in null_values
    )

    return (extreme + 1) / (len(null_values) + 1)


# ================================================================================================
# MAIN
# ================================================================================================

def main() -> None:

    print("=" * 100)
    print("PERMUTATION NULL / MATCHED ANCHOR FACTOR-GEOMETRY EXPERIMENT")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"2M                   = {STEP:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range              = [{T_MIN}, {T_MAX:+d}]")
    print(f"actual trials        = {TRIALS}")
    print(f"permutations         = {PERMUTATIONS}")
    print(f"random seed          = {SEED}")
    print()

    rng = random.Random(SEED)

    # --------------------------------------------------------------------------------------------
    # Generate one fixed actual dataset
    # --------------------------------------------------------------------------------------------

    anchors = generate_unique_semiprimes(
        TRIALS,
        FACTOR_MIN,
        FACTOR_MAX,
        rng
    )

    # --------------------------------------------------------------------------------------------
    # Build the collision landscape ONCE.
    #
    # This means every permutation sees exactly the same events.
    # --------------------------------------------------------------------------------------------

    events = build_collision_landscape(anchors)

    print("=" * 100)
    print("FIXED COLLISION LANDSCAPE")
    print("=" * 100)
    print(f"events = {len(events)}")
    print()

    if not events:
        print("No events found.")
        return

    # Actual assignment.
    actual_assignment = {
        trial: (p, q)
        for trial, (p, q, n) in enumerate(anchors, start=1)
    }

    actual_stats = calculate_statistics(
        events,
        actual_assignment
    )

    # --------------------------------------------------------------------------------------------
    # Permutation null distribution
    # --------------------------------------------------------------------------------------------

    null_values = {
        name: []
        for name in STAT_NAMES
    }

    for permutation in range(1, PERMUTATIONS + 1):

        assignment = permutation_assignment(
            anchors,
            rng
        )

        stats = calculate_statistics(
            events,
            assignment
        )

        for name in STAT_NAMES:
            value = stats[name]

            if math.isfinite(value):
                null_values[name].append(value)

        if permutation % 100 == 0:
            print(
                f"permutation {permutation:4d}/{PERMUTATIONS}"
            )

    # --------------------------------------------------------------------------------------------
    # Results
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL vs PERMUTATION NULL")
    print("=" * 100)

    print(
        f"{'STATISTIC':28s}"
        f"{'ACTUAL':>16s}"
        f"{'NULL MEAN':>16s}"
        f"{'NULL STD':>16s}"
        f"{'Z':>14s}"
        f"{'P(2-sided)':>16s}"
    )

    print("-" * 100)

    for name in STAT_NAMES:

        actual = actual_stats[name]
        values = null_values[name]

        if not values or not math.isfinite(actual):
            print(
                f"{name:28s}"
                f"{str(actual):>16s}"
                f"{'N/A':>16s}"
                f"{'N/A':>16s}"
                f"{'N/A':>14s}"
                f"{'N/A':>16s}"
            )
            continue

        mu = statistics.fmean(values)

        if len(values) >= 2:
            sigma = statistics.stdev(values)
        else:
            sigma = 0.0

        if sigma > 0:
            z = (actual - mu) / sigma
        else:
            z = float("nan")

        p = two_sided_empirical_p(
            actual,
            values
        )

        print(
            f"{name:28s}"
            f"{actual:16.8f}"
            f"{mu:16.8f}"
            f"{sigma:16.8f}"
            f"{z:14.4f}"
            f"{p:16.6f}"
        )

    # --------------------------------------------------------------------------------------------
    # Quantiles
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL POSITION INSIDE NULL DISTRIBUTION")
    print("=" * 100)

    for name in STAT_NAMES:

        actual = actual_stats[name]
        values = sorted(null_values[name])

        if not values or not math.isfinite(actual):
            continue

        below = sum(x <= actual for x in values)
        percentile = below / len(values)

        q025 = values[int(0.025 * (len(values) - 1))]
        q500 = values[int(0.500 * (len(values) - 1))]
        q975 = values[int(0.975 * (len(values) - 1))]

        print(f"{name}")
        print(f"    actual percentile = {percentile:.6f}")
        print(f"    null 2.5%         = {q025:.8f}")
        print(f"    null 50%          = {q500:.8f}")
        print(f"    null 97.5%        = {q975:.8f}")

    # --------------------------------------------------------------------------------------------
    # t-conditioned analysis
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("t-CONDITIONED ACTUAL vs PERMUTATION NULL")
    print("=" * 100)

    events_by_t: dict[int, list[Event]] = defaultdict(list)

    for e in events:
        events_by_t[e.t].append(e)

    print(
        f"{'t':>4s}"
        f"{'N':>6s}"
        f"{'ACT dSum':>14s}"
        f"{'NULL dSum':>14s}"
        f"{'NULL SD':>14s}"
        f"{'Z':>12s}"
    )

    print("-" * 100)

    for t in range(T_MIN, T_MAX + 1):

        if t == 0:
            continue

        tev = events_by_t.get(t)

        if not tev:
            continue

        actual_values = []

        for e in tev:
            p, q = actual_assignment[e.trial]
            actual_values.append(
                float((e.a - p) + (e.b - q))
            )

        actual_mean = statistics.fmean(actual_values)

        permutation_means = []

        for _ in range(250):

            assignment = permutation_assignment(
                anchors,
                rng
            )

            vals = []

            for e in tev:
                p, q = assignment[e.trial]
                vals.append(
                    float((e.a - p) + (e.b - q))
                )

            permutation_means.append(
                statistics.fmean(vals)
            )

        null_mean = statistics.fmean(permutation_means)

        if len(permutation_means) >= 2:
            null_sd = statistics.stdev(permutation_means)
        else:
            null_sd = 0.0

        if null_sd > 0:
            z = (actual_mean - null_mean) / null_sd
        else:
            z = float("nan")

        print(
            f"{t:4d}"
            f"{len(tev):6d}"
            f"{actual_mean:14.2f}"
            f"{null_mean:14.2f}"
            f"{null_sd:14.2f}"
            f"{z:12.4f}"
        )

    # --------------------------------------------------------------------------------------------
    # Identity sanity check
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("LANDSCAPE SANITY CHECK")
    print("=" * 100)

    landscape = {
        (e.trial, e.t, e.n, e.a, e.b)
        for e in events
    }

    print(f"unique collision states = {len(landscape)}")
    print(f"event list              = {len(events)}")

    if len(landscape) == len(events):
        print("duplicate collision states = 0")
    else:
        print(
            "duplicate collision states = "
            f"{len(events) - len(landscape)}"
        )

    # --------------------------------------------------------------------------------------------
    # Final interpretation helper
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("INTERPRETATION GUIDE")
    print("=" * 100)
    print(
        "A statistic is interesting only if the actual value lies unusually far "
        "from the permutation null distribution."
    )
    print(
        "The strongest evidence would be a large |Z| together with a small "
        "empirical two-sided p-value across a preselected statistic."
    )
    print(
        "The t-conditioned table tests whether any apparent geometry is still "
        "unusual after fixing t."
    )
    print()
    print("Experiment complete.")


if __name__ == "__main__":
    main()
