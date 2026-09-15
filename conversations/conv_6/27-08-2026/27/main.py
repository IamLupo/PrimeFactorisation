#!/usr/bin/env python3

import math
import random
import statistics
from dataclasses import dataclass
from collections import defaultdict


# ================================================================================================
# CONFIG
# ================================================================================================

M = 111_546_435
STEP = 2 * M

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

T_MIN = -25
T_MAX = 25

TRIALS = 300

# Number of approximately equal-size n groups.
# With 300 anchors, 20 means ~15 anchors/block.
BLOCKS = 20

PERMUTATIONS = 5000

SEED = 1511464998


# ================================================================================================
# DATA
# ================================================================================================

@dataclass(frozen=True)
class Anchor:
    trial: int
    p: int
    q: int
    n: int


@dataclass(frozen=True)
class Event:
    trial: int
    n: int
    a: int
    b: int
    t: int


# ================================================================================================
# PRIME SIEVE
# ================================================================================================

def sieve(limit):

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


# ================================================================================================
# GENERATE UNIQUE SEMIPRIME ANCHORS
# ================================================================================================

def generate_anchors(count, rng):

    primes = [
        p
        for p in sieve(FACTOR_MAX)
        if p >= FACTOR_MIN
    ]

    seen = set()
    anchors = []

    while len(anchors) < count:

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

        anchors.append(
            Anchor(
                trial=len(anchors) + 1,
                p=p,
                q=q,
                n=n
            )
        )

    return anchors


# ================================================================================================
# FACTOR PAIRS
# ================================================================================================

def factor_pairs(n):

    result = []

    r = math.isqrt(n)

    start = FACTOR_MIN
    stop = min(r, FACTOR_MAX)

    for a in range(start, stop + 1):

        if n % a != 0:
            continue

        b = n // a

        if b > FACTOR_MAX:
            continue

        if a <= b:
            result.append((a, b))

    return result


# ================================================================================================
# BUILD FIXED COLLISION LANDSCAPE
# ================================================================================================

def build_events(anchors):

    events = []

    for anchor in anchors:

        for t in range(T_MIN, T_MAX + 1):

            if t == 0:
                continue

            Nt = anchor.n + STEP * t

            if Nt <= 0:
                continue

            for a, b in factor_pairs(Nt):

                events.append(
                    Event(
                        trial=anchor.trial,
                        n=anchor.n,
                        a=a,
                        b=b,
                        t=t
                    )
                )

    return events


# ================================================================================================
# PEARSON
# ================================================================================================

def pearson(x, y):

    n = len(x)

    if n < 2:
        return float("nan")

    sx = sum(x)
    sy = sum(y)

    mx = sx / n
    my = sy / n

    num = 0.0
    xx = 0.0
    yy = 0.0

    for a, b in zip(x, y):

        dx = a - mx
        dy = b - my

        num += dx * dy
        xx += dx * dx
        yy += dy * dy

    den = math.sqrt(xx * yy)

    if den == 0:
        return float("nan")

    return num / den


# ================================================================================================
# FAST CORRELATION OF t WITH dSum
# ================================================================================================

def correlation_from_assignment(events, assignment):

    # We only need:
    #
    # t
    # dSum = (a-p) + (b-q)
    #
    # Since:
    #
    # dSum = (a+b) - (p+q)
    #
    # precompute event-side a+b.

    n = len(events)

    sum_t = 0.0
    sum_d = 0.0

    for e in events:

        p, q = assignment[e.trial]

        d = (e.a + e.b) - (p + q)

        sum_t += e.t
        sum_d += d

    mt = sum_t / n
    md = sum_d / n

    numerator = 0.0
    tt = 0.0
    dd = 0.0

    for e in events:

        p, q = assignment[e.trial]

        d = (e.a + e.b) - (p + q)

        dt = e.t - mt
        ddv = d - md

        numerator += dt * ddv
        tt += dt * dt
        dd += ddv * ddv

    den = math.sqrt(tt * dd)

    if den == 0:
        return float("nan")

    return numerator / den


# ================================================================================================
# MORE STATISTICS
# ================================================================================================

def assignment_statistics(events, assignment):

    n = len(events)

    sum_t = 0.0
    sum_d = 0.0

    same = 0

    dsum_values = []

    for e in events:

        p, q = assignment[e.trial]

        dp = e.a - p
        dq = e.b - q

        d = dp + dq

        sum_t += e.t
        sum_d += d

        dsum_values.append(d)

        if (dp >= 0 and dq >= 0) or (dp < 0 and dq < 0):
            same += 1

    mean_t = sum_t / n
    mean_d = sum_d / n

    num = 0.0
    den_t = 0.0
    den_d = 0.0

    for e in events:

        p, q = assignment[e.trial]

        d = (e.a - p) + (e.b - q)

        dt = e.t - mean_t
        dd = d - mean_d

        num += dt * dd
        den_t += dt * dt
        den_d += dd * dd

    corr = num / math.sqrt(den_t * den_d)

    return {
        "corr_t_dsum": corr,
        "mean_dsum": mean_d,
        "same_sign_fraction": same / n,
    }


# ================================================================================================
# SIZE BLOCKS
# ================================================================================================

def make_size_blocks(anchors):

    ordered = sorted(
        anchors,
        key=lambda x: x.n
    )

    blocks = []

    # Exactly BLOCKS blocks as evenly as possible.
    base = len(ordered) // BLOCKS
    remainder = len(ordered) % BLOCKS

    pos = 0

    for i in range(BLOCKS):

        size = base + (1 if i < remainder else 0)

        block = ordered[pos:pos + size]

        blocks.append(block)

        pos += size

    return blocks


# ================================================================================================
# BLOCK PERMUTATION
# ================================================================================================

def permute_factor_pairs_within_blocks(blocks, rng):

    assignment = {}

    for block in blocks:

        pairs = [
            (a.p, a.q)
            for a in block
        ]

        rng.shuffle(pairs)

        for anchor, pair in zip(block, pairs):

            assignment[anchor.trial] = pair

    return assignment


# ================================================================================================
# EMPIRICAL P VALUE
# ================================================================================================

def empirical_p(actual, null):

    mu = statistics.fmean(null)

    distance = abs(actual - mu)

    extreme = sum(
        abs(x - mu) >= distance
        for x in null
    )

    return (extreme + 1) / (len(null) + 1)


# ================================================================================================
# QUANTILE
# ================================================================================================

def quantile(values, q):

    values = sorted(values)

    if not values:
        return float("nan")

    pos = q * (len(values) - 1)

    lo = math.floor(pos)
    hi = math.ceil(pos)

    if lo == hi:
        return values[lo]

    f = pos - lo

    return (
        values[lo] * (1.0 - f)
        + values[hi] * f
    )


# ================================================================================================
# T-CONDITIONED DATA
# ================================================================================================

def actual_t_curve(events, assignment):

    by_t = defaultdict(list)

    for e in events:
        by_t[e.t].append(e)

    result = {}

    for t in sorted(by_t):

        vals = []

        for e in by_t[t]:

            p, q = assignment[e.trial]

            vals.append(
                (e.a - p) + (e.b - q)
            )

        result[t] = {
            "N": len(vals),
            "mean": statistics.fmean(vals)
        }

    return result


# ================================================================================================
# MAIN
# ================================================================================================

def main():

    print("=" * 100)
    print("SIZE-CONTROLLED BLOCK-PERMUTATION EXPERIMENT")
    print("=" * 100)

    print(f"M                    = {M:,}")
    print(f"2M                   = {STEP:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range              = [{T_MIN}, {T_MAX:+d}]")
    print(f"actual trials        = {TRIALS}")
    print(f"size blocks          = {BLOCKS}")
    print(f"permutations         = {PERMUTATIONS}")
    print(f"random seed          = {SEED}")

    rng = random.Random(SEED)

    # --------------------------------------------------------------------------------------------
    # Anchors
    # --------------------------------------------------------------------------------------------

    anchors = generate_anchors(
        TRIALS,
        rng
    )

    # --------------------------------------------------------------------------------------------
    # Fixed landscape
    # --------------------------------------------------------------------------------------------

    events = build_events(anchors)

    print()
    print("=" * 100)
    print("FIXED COLLISION LANDSCAPE")
    print("=" * 100)

    print(f"events = {len(events)}")

    # --------------------------------------------------------------------------------------------
    # Actual assignment
    # --------------------------------------------------------------------------------------------

    actual_assignment = {
        a.trial: (a.p, a.q)
        for a in anchors
    }

    actual = assignment_statistics(
        events,
        actual_assignment
    )

    # --------------------------------------------------------------------------------------------
    # Size blocks
    # --------------------------------------------------------------------------------------------

    blocks = make_size_blocks(anchors)

    print()
    print("=" * 100)
    print("SIZE BLOCKS")
    print("=" * 100)

    for i, block in enumerate(blocks):

        print(
            f"block={i:02d} "
            f"N={len(block):3d} "
            f"n=[{block[0].n:,} ... {block[-1].n:,}]"
        )

    # --------------------------------------------------------------------------------------------
    # Actual t curve
    # --------------------------------------------------------------------------------------------

    curve = actual_t_curve(
        events,
        actual_assignment
    )

    print()
    print("=" * 100)
    print("ACTUAL t-CURVE")
    print("=" * 100)

    print(
        f"{'t':>4s}"
        f"{'N':>6s}"
        f"{'mean dSum':>18s}"
    )

    print("-" * 100)

    for t, row in curve.items():

        print(
            f"{t:4d}"
            f"{row['N']:6d}"
            f"{row['mean']:18.3f}"
        )

    # --------------------------------------------------------------------------------------------
    # Null
    # --------------------------------------------------------------------------------------------

    null_corr = []
    null_mean_dsum = []
    null_same = []

    print()
    print("=" * 100)
    print("BLOCK PERMUTATIONS")
    print("=" * 100)

    for i in range(1, PERMUTATIONS + 1):

        assignment = permute_factor_pairs_within_blocks(
            blocks,
            rng
        )

        stats = assignment_statistics(
            events,
            assignment
        )

        null_corr.append(
            stats["corr_t_dsum"]
        )

        null_mean_dsum.append(
            stats["mean_dsum"]
        )

        null_same.append(
            stats["same_sign_fraction"]
        )

        if i % 500 == 0:
            print(
                f"permutation {i:5d}/{PERMUTATIONS}"
            )

    # --------------------------------------------------------------------------------------------
    # Results
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL vs SIZE-CONTROLLED NULL")
    print("=" * 100)

    print(
        f"{'STATISTIC':28s}"
        f"{'ACTUAL':>16s}"
        f"{'NULL MEAN':>16s}"
        f"{'NULL STD':>16s}"
        f"{'Z':>12s}"
        f"{'P(2-sided)':>14s}"
    )

    print("-" * 100)

    datasets = [
        (
            "corr_t_dsum",
            actual["corr_t_dsum"],
            null_corr
        ),
        (
            "mean_dsum",
            actual["mean_dsum"],
            null_mean_dsum
        ),
        (
            "same_sign_fraction",
            actual["same_sign_fraction"],
            null_same
        ),
    ]

    for name, actual_value, values in datasets:

        mu = statistics.fmean(values)
        sd = statistics.stdev(values)

        z = (
            (actual_value - mu) / sd
            if sd != 0
            else float("nan")
        )

        p = empirical_p(
            actual_value,
            values
        )

        print(
            f"{name:28s}"
            f"{actual_value:16.8f}"
            f"{mu:16.8f}"
            f"{sd:16.8f}"
            f"{z:12.4f}"
            f"{p:14.6f}"
        )

    # --------------------------------------------------------------------------------------------
    # Percentiles
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("NULL PERCENTILES")
    print("=" * 100)

    percentile_data = [
        (
            "corr_t_dsum",
            actual["corr_t_dsum"],
            null_corr
        ),
        (
            "mean_dsum",
            actual["mean_dsum"],
            null_mean_dsum
        ),
        (
            "same_sign_fraction",
            actual["same_sign_fraction"],
            null_same
        ),
    ]

    for name, actual_value, values in percentile_data:

        q025 = quantile(values, 0.025)
        q50 = quantile(values, 0.50)
        q975 = quantile(values, 0.975)

        below = sum(
            x <= actual_value
            for x in values
        )

        percentile = below / len(values)

        print(name)
        print(f"    actual percentile = {percentile:.6f}")
        print(f"    null 2.5%         = {q025:.8f}")
        print(f"    null 50%          = {q50:.8f}")
        print(f"    null 97.5%        = {q975:.8f}")

    # --------------------------------------------------------------------------------------------
    # Compare with unrestricted-null result
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("COMPARISON WITH PREVIOUS UNRESTRICTED NULL")
    print("=" * 100)

    print()
    print("Previous unrestricted permutation:")
    print("    corr(t,dSum) actual = 0.94108869")
    print("    null mean           = 0.42395873")
    print("    Z                   = 13.9167")
    print("    p                   = 0.000999")

    print()
    print("This experiment:")
    print("    shuffles factor pairs ONLY within similar-n blocks.")
    print()
    print("Therefore a surviving large correlation is evidence that")
    print("the observed t -> dSum relationship is not merely an")
    print("artifact of assigning different factor-pair scales to")
    print("different n scales.")

    # --------------------------------------------------------------------------------------------
    # Checks
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("LANDSCAPE SANITY")
    print("=" * 100)

    states = {
        (
            e.trial,
            e.n,
            e.t,
            e.a,
            e.b
        )
        for e in events
    }

    print(
        f"event list              = {len(events)}"
    )

    print(
        f"unique states           = {len(states)}"
    )

    print(
        f"duplicate states        = "
        f"{len(events) - len(states)}"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
