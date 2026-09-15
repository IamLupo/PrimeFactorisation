#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict
from statistics import mean

# ==============================================================================
# CONFIG
# ==============================================================================

SEED = 0x5A17_2026
random.seed(SEED)

M = 111_546_435
STEP = 2 * M

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
TMAX = 25

RANDOM_ANCHORS_PER_TRIAL = 20
ANCHOR_DELTA = 5_000_000

TS = [t for t in range(-TMAX, TMAX + 1) if t != 0]


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def sieve(limit):
    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(limit) + 1):
        if s[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            s[start:limit + 1:p] = b"\x00" * count

    return [i for i in range(FACTOR_MIN, limit + 1)
            if s[i]]


FACTOR_PRIMES = sieve(FACTOR_MAX)
FACTOR_SET = set(FACTOR_PRIMES)

print("=" * 100)
print("THREE-WAY ANCHOR / FACTOR-GEOMETRY EXPERIMENT")
print("=" * 100)
print(f"M                    = {M:,}")
print(f"2M                   = {STEP:,}")
print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
print(f"T range              = [-{TMAX}, +{TMAX}]")
print(f"actual trials        = {TRIALS}")
print(f"random anchors/trial = {RANDOM_ANCHORS_PER_TRIAL}")
print(f"anchor delta         = +/-{ANCHOR_DELTA:,}")
print()


# ==============================================================================
# RANDOM ACTUAL SEMIPRIME
# ==============================================================================

def generate_actual():
    """
    Generate a unique semiprime p*q where both p,q are primes in [10k,100k].
    """

    while True:
        p = random.choice(FACTOR_PRIMES)
        q = random.choice(FACTOR_PRIMES)

        if p > q:
            p, q = q, p

        n = p * q
        return n, p, q


# ==============================================================================
# DIRECT FACTORIZATION
# ==============================================================================

def find_factor_pair(n):
    """
    Find a factorization n=a*b where:
        a,b are prime
        10,000 <= a <= b <= 100,000

    Returns (a,b), or None.
    """

    limit = min(math.isqrt(n), FACTOR_MAX)

    for a in FACTOR_PRIMES:

        if a > limit:
            break

        if n % a != 0:
            continue

        b = n // a

        if (
            a >= FACTOR_MIN
            and a <= b <= FACTOR_MAX
            and b in FACTOR_SET
        ):
            return a, b

    return None


# ==============================================================================
# RANDOM CONTROL
# ==============================================================================

def random_anchor_near(n):
    while True:

        delta = random.randint(-ANCHOR_DELTA, ANCHOR_DELTA)

        if delta == 0:
            continue

        x = n + delta

        if x <= 0:
            continue

        if math.gcd(x, STEP) != 1:
            continue

        return x


# ==============================================================================
# BUILD UNIQUE ACTUAL DATASET
# ==============================================================================

actual_trials = []
used_n = set()

while len(actual_trials) < TRIALS:

    n, p, q = generate_actual()

    if n in used_n:
        continue

    used_n.add(n)

    actual_trials.append({
        "n": n,
        "p": p,
        "q": q,
    })


# ==============================================================================
# SHUFFLED FACTOR ATTACHMENTS
# ==============================================================================

shuffled_pairs = [
    (x["p"], x["q"])
    for x in actual_trials
]

random.shuffle(shuffled_pairs)

shuffled_trials = []

for x, (p, q) in zip(actual_trials, shuffled_pairs):

    shuffled_trials.append({
        "n": x["n"],
        "p": p,
        "q": q,
    })


# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

populations = ("actual", "shuffled", "random")

per_t = {
    name: Counter()
    for name in populations
}

event_records = {
    name: []
    for name in populations
}

sign_counts = {
    name: Counter()
    for name in populations
}

# Actual and shuffled MUST have identical collision landscapes.
landscape = {
    "actual": set(),
    "shuffled": set(),
}

random_anchor_count = 0


# ==============================================================================
# EVENT SCANNER
# ==============================================================================

def scan(anchor_n, p=None, q=None):

    events = []

    for t in TS:

        Nt = anchor_n + STEP * t

        if Nt <= 0:
            continue

        pair = find_factor_pair(Nt)

        if pair is None:
            continue

        a, b = pair

        event = {
            "t": t,
            "anchor_n": anchor_n,
            "Nt": Nt,
            "a": a,
            "b": b,
        }

        if p is not None:

            event["p"] = p
            event["q"] = q

            event["dp"] = a - p
            event["dq"] = b - q

            event["dSum"] = (a + b) - (p + q)
            event["dDiff"] = (b - a) - (q - p)

            scale = math.sqrt(abs(STEP * t))

            event["ndp"] = event["dp"] / scale
            event["ndq"] = event["dq"] / scale

        events.append(event)

    return events


# ==============================================================================
# ACTUAL
# ==============================================================================

for trial_id, x in enumerate(actual_trials, 1):

    events = scan(
        x["n"],
        x["p"],
        x["q"],
    )

    for e in events:

        key = (e["t"], e["Nt"])

        landscape["actual"].add(key)

        per_t["actual"][e["t"]] += 1

        sign = (
            "same"
            if e["dp"] * e["dq"] > 0
            else "opposite"
            if e["dp"] * e["dq"] < 0
            else "zero"
        )

        sign_counts["actual"][sign] += 1

        event_records["actual"].append(
            (trial_id, e)
        )


# ==============================================================================
# SHUFFLED
# ==============================================================================

for trial_id, x in enumerate(shuffled_trials, 1):

    events = scan(
        x["n"],
        x["p"],
        x["q"],
    )

    for e in events:

        key = (e["t"], e["Nt"])

        landscape["shuffled"].add(key)

        per_t["shuffled"][e["t"]] += 1

        sign = (
            "same"
            if e["dp"] * e["dq"] > 0
            else "opposite"
            if e["dp"] * e["dq"] < 0
            else "zero"
        )

        sign_counts["shuffled"][sign] += 1

        event_records["shuffled"].append(
            (trial_id, e)
        )


# ==============================================================================
# RANDOM CONTROLS
# ==============================================================================

for actual in actual_trials:

    for _ in range(RANDOM_ANCHORS_PER_TRIAL):

        random_anchor_count += 1

        n0 = random_anchor_near(actual["n"])

        events = scan(n0)

        for e in events:

            per_t["random"][e["t"]] += 1

            event_records["random"].append(
                (random_anchor_count, e)
            )


# ==============================================================================
# SUMMARY
# ==============================================================================

print("=" * 100)
print("SUMMARY")
print("=" * 100)

for name in populations:

    events = sum(per_t[name].values())

    if name == "random":
        anchors = random_anchor_count
    else:
        anchors = TRIALS

    print(
        f"{name:10s} "
        f"anchors={anchors:5d} "
        f"events={events:5d} "
        f"events/anchor={events/anchors:.8f}"
    )

print()


# ==============================================================================
# LANDSCAPE INVARIANT
# ==============================================================================

print("=" * 100)
print("LANDSCAPE INVARIANT")
print("=" * 100)

print(
    f"actual (t,Nt) events   = "
    f"{len(landscape['actual']):5d}"
)

print(
    f"shuffled (t,Nt) events = "
    f"{len(landscape['shuffled']):5d}"
)

print(
    "identical collision landscape =",
    landscape["actual"] == landscape["shuffled"]
)

print()


# ==============================================================================
# PER-t
# ==============================================================================

print("=" * 100)
print("t DISTRIBUTION")
print("=" * 100)

print(
    f"{'t':>4} "
    f"{'ACT':>7} "
    f"{'SHUF':>7} "
    f"{'RANDOM':>7}"
)

for t in TS:

    print(
        f"{t:+4d} "
        f"{per_t['actual'][t]:7d} "
        f"{per_t['shuffled'][t]:7d} "
        f"{per_t['random'][t]:7d}"
    )

print()


# ==============================================================================
# SYMMETRIC COMPARISON
# ==============================================================================

print("=" * 100)
print("SYMMETRIC |t| COMPARISON")
print("=" * 100)

print(
    f"{'|t|':>4} "
    f"{'A-':>6} {'A+':>6} "
    f"{'S-':>6} {'S+':>6} "
    f"{'R-':>6} {'R+':>6}"
)

for k in range(1, TMAX + 1):

    print(
        f"{k:4d} "
        f"{per_t['actual'][-k]:6d} "
        f"{per_t['actual'][+k]:6d} "
        f"{per_t['shuffled'][-k]:6d} "
        f"{per_t['shuffled'][+k]:6d} "
        f"{per_t['random'][-k]:6d} "
        f"{per_t['random'][+k]:6d}"
    )

print()


# ==============================================================================
# GEOMETRY SUMMARY
# ==============================================================================

def geometry_summary(name):

    out = {}

    for t in TS:

        events = [
            e
            for _, e in event_records[name]
            if e["t"] == t and "dp" in e
        ]

        if not events:
            continue

        dps = [e["dp"] for e in events]
        dqs = [e["dq"] for e in events]

        dsum = [e["dSum"] for e in events]
        ddiff = [e["dDiff"] for e in events]

        ndp = [e["ndp"] for e in events]
        ndq = [e["ndq"] for e in events]

        out[t] = {
            "count": len(events),

            "dp_min": min(dps),
            "dp_max": max(dps),

            "dq_min": min(dqs),
            "dq_max": max(dqs),

            "dSum_mean": mean(dsum),
            "dDiff_mean": mean(ddiff),

            "ndp_mean": mean(ndp),
            "ndq_mean": mean(ndq),
        }

    return out


actual_geometry = geometry_summary("actual")
shuffled_geometry = geometry_summary("shuffled")


print("=" * 100)
print("ACTUAL vs SHUFFLED GEOMETRY")
print("=" * 100)

print(
    f"{'t':>4} "
    f"{'A_nDP':>10} {'S_nDP':>10} "
    f"{'A_nDQ':>10} {'S_nDQ':>10} "
    f"{'A_dSum':>12} {'S_dSum':>12}"
)

for t in TS:

    if t not in actual_geometry:
        continue

    if t not in shuffled_geometry:
        continue

    a = actual_geometry[t]
    s = shuffled_geometry[t]

    print(
        f"{t:+4d} "
        f"{a['ndp_mean']:+10.6f} "
        f"{s['ndp_mean']:+10.6f} "
        f"{a['ndq_mean']:+10.6f} "
        f"{s['ndq_mean']:+10.6f} "
        f"{a['dSum_mean']:+12.1f} "
        f"{s['dSum_mean']:+12.1f}"
    )

print()


# ==============================================================================
# GLOBAL GEOMETRY
# ==============================================================================

print("=" * 100)
print("GLOBAL GEOMETRY")
print("=" * 100)

for name in ("actual", "shuffled"):

    events = [
        e
        for _, e in event_records[name]
        if "dp" in e
    ]

    if not events:
        continue

    dp = [e["dp"] for e in events]
    dq = [e["dq"] for e in events]

    dsum = [e["dSum"] for e in events]
    ddiff = [e["dDiff"] for e in events]

    ndp = [e["ndp"] for e in events]
    ndq = [e["ndq"] for e in events]

    print(name)

    print(f"  events       = {len(events)}")

    print(
        f"  dp range     = "
        f"{min(dp):+d} ... {max(dp):+d}"
    )

    print(
        f"  dq range     = "
        f"{min(dq):+d} ... {max(dq):+d}"
    )

    print(
        f"  dSum range   = "
        f"{min(dsum):+d} ... {max(dsum):+d}"
    )

    print(
        f"  dDiff range  = "
        f"{min(ddiff):+d} ... {max(ddiff):+d}"
    )

    print(
        f"  mean ndp     = {mean(ndp):+.8f}"
    )

    print(
        f"  mean ndq     = {mean(ndq):+.8f}"
    )

    print(
        f"  same-sign    = {sign_counts[name]['same']}"
    )

    print(
        f"  opposite     = {sign_counts[name]['opposite']}"
    )

    print()


# ==============================================================================
# REPEATED EXACT GEOMETRY
# ==============================================================================

print("=" * 100)
print("REPEATED EXACT (t,dp,dq)")
print("=" * 100)

for name in ("actual", "shuffled"):

    counts = Counter()

    for _, e in event_records[name]:

        if "dp" not in e:
            continue

        key = (
            e["t"],
            e["dp"],
            e["dq"],
        )

        counts[key] += 1

    repeated = [
        (count, key)
        for key, count in counts.items()
        if count > 1
    ]

    repeated.sort(reverse=True)

    print(name)

    if not repeated:
        print("  none")
    else:
        for count, (t, dp, dq) in repeated[:25]:

            print(
                f"  count={count:3d} "
                f"t={t:+3d} "
                f"dp={dp:+7d} "
                f"dq={dq:+7d}"
            )

    print()


# ==============================================================================
# EXTREMES
# ==============================================================================

print("=" * 100)
print("EXTREME ACTUAL EVENTS")
print("=" * 100)

actual_events = [
    e
    for _, e in event_records["actual"]
    if "dp" in e
]

for e in sorted(
    actual_events,
    key=lambda x: abs(x["dp"]) + abs(x["dq"]),
    reverse=True
)[:25]:

    print(
        f"t={e['t']:+3d} "
        f"actual={e['p']:,}*{e['q']:,} "
        f"shifted={e['a']:,}*{e['b']:,} "
        f"dp={e['dp']:+7d} "
        f"dq={e['dq']:+7d} "
        f"dSum={e['dSum']:+8d} "
        f"dDiff={e['dDiff']:+8d}"
    )

print()


# ==============================================================================
# FINAL
# ==============================================================================

print("=" * 100)
print("FINAL")
print("=" * 100)

print(
    "Actual/shuffled (t,Nt) landscape identical:",
    landscape["actual"] == landscape["shuffled"],
)

print(
    "The shuffled experiment changes only the association between"
)
print(
    "the original anchor factor pair (p,q) and the same shifted Nt."
)

print()
print("Experiment complete.")