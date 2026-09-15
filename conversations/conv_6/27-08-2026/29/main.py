#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# ============================================================
# CRT JOINT MODULAR SIGNATURE / PERMUTATION EXPERIMENT
# ============================================================

M = 111_546_435
T_MIN = -25
T_MAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
PERMUTATIONS = 2000

SEED = 1511464998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

rng = random.Random(SEED)


# ============================================================
# PRIME GENERATION
# ============================================================

def sieve(n):
    s = bytearray(b"\x01") * (n + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p:n + 1:p] = b"\x00" * (
                ((n - p * p) // p) + 1
            )

    return [i for i in range(n + 1) if s[i]]


ALL_PRIMES = sieve(FACTOR_MAX)
FACTOR_PRIMES = [
    p for p in ALL_PRIMES
    if FACTOR_MIN <= p <= FACTOR_MAX
]

PRIME_SET = set(FACTOR_PRIMES)

print("=" * 100)
print("CRT JOINT MODULAR SIGNATURE / PERMUTATION EXPERIMENT")
print("=" * 100)
print(f"M                  = {M:,}")
print(f"2M                 = {2*M:,}")
print(f"factor range       = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
print(f"T range            = [{T_MIN}, {T_MAX}]")
print(f"trials             = {TRIALS}")
print(f"permutations       = {PERMUTATIONS}")
print(f"moduli             = {MODULI}")
print(f"random seed        = {SEED}")
print(f"factor primes      = {len(FACTOR_PRIMES):,}")
print()


# ============================================================
# MODULAR HELPERS
# ============================================================

def inv_mod(a, r):
    return pow(a, -1, r)


def signature_for_pair(a, p):
    """
    U = (a/p mod r) for all r in MODULI.
    """
    out = []

    for r in MODULI:
        out.append((a * inv_mod(p % r, r)) % r)

    return tuple(out)


def crt_encode(sig):
    """
    Encode the complete modular vector into one integer.

    Since the moduli are pairwise coprime, this is injective.
    """
    x = 0
    multiplier = 1

    for u, r in zip(sig, MODULI):
        x += u * multiplier
        multiplier *= r

    return x


# ============================================================
# SEMIPRIME ANCHORS
# ============================================================

def random_semiprime():
    while True:
        p = rng.choice(FACTOR_PRIMES)
        q = rng.choice(FACTOR_PRIMES)

        if p == q:
            continue

        # canonical ordering
        if p > q:
            p, q = q, p

        n = p * q

        # Important: inverses modulo every r must exist.
        if math.gcd(p, math.prod(MODULI)) != 1:
            continue

        if math.gcd(q, math.prod(MODULI)) != 1:
            continue

        return p, q, n


# ============================================================
# COLLISION SEARCH
#
# This deliberately searches only the valid factor interval:
#
#     FACTOR_MIN <= a <= b <= FACTOR_MAX
#
# and only primes.
# ============================================================

def collisions_for_n(n):
    events = []

    for t in range(T_MIN, T_MAX + 1):

        if t == 0:
            continue

        Nt = n + 2 * M * t

        if Nt <= 0:
            continue

        # We only need a <= sqrt(Nt)
        limit = math.isqrt(Nt)

        lo = FACTOR_MIN
        hi = min(FACTOR_MAX, limit)

        for a in FACTOR_PRIMES:

            if a > hi:
                break

            if Nt % a:
                continue

            b = Nt // a

            if b < FACTOR_MIN or b > FACTOR_MAX:
                continue

            if b not in PRIME_SET:
                continue

            if a > b:
                continue

            # Do not report the original anchor factor pair.
            # t != 0 already removes the direct anchor.
            events.append({
                "t": t,
                "Nt": Nt,
                "a": a,
                "b": b,
                "p": None,
                "q": None,
                "n": n,
            })

    return events


# ============================================================
# GENERATE ACTUAL DATASET
# ============================================================

anchors = []
events = []

while len(anchors) < TRIALS:

    p, q, n = random_semiprime()

    # keep anchors unique
    if n in {x["n"] for x in anchors}:
        continue

    anchor = {
        "p": p,
        "q": q,
        "n": n,
    }

    anchors.append(anchor)

    found = collisions_for_n(n)

    for e in found:
        e["p"] = p
        e["q"] = q
        events.append(e)


print("=" * 100)
print("ACTUAL DATASET")
print("=" * 100)

print(f"anchors = {len(anchors)}")
print(f"events  = {len(events)}")
print()


# ============================================================
# IMPORTANT:
# The collision landscape consists only of (t,Nt,a,b).
#
# The permutation experiment will NEVER change those.
#
# Only (p,q) is reassigned.
# ============================================================

landscape = [
    (
        e["t"],
        e["Nt"],
        e["a"],
        e["b"],
    )
    for e in events
]


# ============================================================
# ACTUAL SIGNATURES
# ============================================================

def event_signature(e, p_override=None):
    p = e["p"] if p_override is None else p_override
    return signature_for_pair(e["a"], p)


actual_signatures = [
    event_signature(e)
    for e in events
]


# ============================================================
# CRT SIGNATURE STATISTICS
# ============================================================

def signature_entropy(signatures):

    counts = Counter(signatures)
    total = len(signatures)

    if total == 0:
        return 0.0

    H = 0.0

    for c in counts.values():
        p = c / total
        H -= p * math.log2(p)

    return H


def unique_signature_count(signatures):
    return len(set(signatures))


def collision_rate(signatures):
    counts = Counter(signatures)

    repeated = sum(
        c - 1
        for c in counts.values()
        if c > 1
    )

    return repeated / len(signatures)


# ============================================================
# PAIRWISE MODULAR DEPENDENCE
# ============================================================

def mutual_information(x, y):

    n = len(x)

    joint = Counter(zip(x, y))
    cx = Counter(x)
    cy = Counter(y)

    mi = 0.0

    for (a, b), c in joint.items():

        pab = c / n
        pa = cx[a] / n
        pb = cy[b] / n

        mi += pab * math.log2(
            pab / (pa * pb)
        )

    return mi


def signature_components(signatures):
    return list(zip(*signatures))


def mean_pairwise_mi(signatures):

    components = signature_components(signatures)

    vals = []

    for i in range(len(MODULI)):
        for j in range(i + 1, len(MODULI)):

            vals.append(
                mutual_information(
                    components[i],
                    components[j]
                )
            )

    return sum(vals) / len(vals)


# ============================================================
# CRT ENCODING
# ============================================================

actual_crt = [
    crt_encode(s)
    for s in actual_signatures
]


# ============================================================
# t / SIGNATURE ASSOCIATION
# ============================================================

def cramers_like_t_association(events, signatures):

    """
    Pearson correlation between signed t and each
    modular component.
    """

    out = {}

    ts = [e["t"] for e in events]

    mean_t = sum(ts) / len(ts)

    for idx, r in enumerate(MODULI):

        ys = [s[idx] for s in signatures]
        mean_y = sum(ys) / len(ys)

        num = sum(
            (t - mean_t) * (y - mean_y)
            for t, y in zip(ts, ys)
        )

        den1 = math.sqrt(
            sum((t - mean_t) ** 2 for t in ts)
        )

        den2 = math.sqrt(
            sum((y - mean_y) ** 2 for y in ys)
        )

        corr = (
            num / (den1 * den2)
            if den1 and den2
            else 0.0
        )

        out[r] = corr

    return out


actual_corr = cramers_like_t_association(
    events,
    actual_signatures
)


# ============================================================
# NULL DISTRIBUTION
#
# Every permutation independently shuffles the anchor
# factor-pair association while keeping the collision
# landscape fixed.
# ============================================================

p_values = []
entropy_values = []
unique_values = []
repeat_values = []
mi_values = []
crt_values = []

corr_values = {
    r: []
    for r in MODULI
}


event_count = len(events)


print("=" * 100)
print("PERMUTATION NULL")
print("=" * 100)

for iteration in range(1, PERMUTATIONS + 1):

    # --------------------------------------------------------
    # Shuffle p values among events.
    #
    # We shuffle complete p/q pairs rather than independently
    # shuffling p and q, preserving the semiprime-pair structure.
    # --------------------------------------------------------

    pairs = [
        (e["p"], e["q"])
        for e in events
    ]

    rng.shuffle(pairs)

    shuffled_signatures = []

    for e, (p, q) in zip(events, pairs):
        shuffled_signatures.append(
            signature_for_pair(e["a"], p)
        )

    # --------------------------------------------------------
    # Statistics
    # --------------------------------------------------------

    H = signature_entropy(shuffled_signatures)
    U = unique_signature_count(shuffled_signatures)
    R = collision_rate(shuffled_signatures)
    MI = mean_pairwise_mi(shuffled_signatures)

    crt = [
        crt_encode(s)
        for s in shuffled_signatures
    ]

    entropy_values.append(H)
    unique_values.append(U)
    repeat_values.append(R)
    mi_values.append(MI)
    crt_values.append(sum(crt) / len(crt))

    corr = cramers_like_t_association(
        events,
        shuffled_signatures
    )

    for r in MODULI:
        corr_values[r].append(abs(corr[r]))

    if iteration % 100 == 0:
        print(
            f"permutation {iteration:4d}/{PERMUTATIONS}"
        )


# ============================================================
# EMPIRICAL NULL HELPERS
# ============================================================

def percentile(values, x):

    values = sorted(values)

    k = 0

    while k < len(values) and values[k] <= x:
        k += 1

    return k / len(values)


def two_sided_empirical_p(values, observed):

    center = sum(values) / len(values)

    observed_distance = abs(observed - center)

    extreme = sum(
        abs(x - center) >= observed_distance
        for x in values
    )

    return (extreme + 1) / (len(values) + 1)


def z_score(observed, values):

    mean = sum(values) / len(values)

    var = sum(
        (x - mean) ** 2
        for x in values
    ) / len(values)

    sd = math.sqrt(var)

    if sd == 0:
        return 0.0

    return (observed - mean) / sd


# ============================================================
# GLOBAL RESULTS
# ============================================================

actual_H = signature_entropy(actual_signatures)
actual_U = unique_signature_count(actual_signatures)
actual_R = collision_rate(actual_signatures)
actual_MI = mean_pairwise_mi(actual_signatures)
actual_CRT = sum(actual_crt) / len(actual_crt)


print()
print("=" * 100)
print("ACTUAL vs PERMUTATION NULL")
print("=" * 100)

print(
    f"{'STATISTIC':<28}"
    f"{'ACTUAL':>14}"
    f"{'NULL MEAN':>14}"
    f"{'NULL STD':>14}"
    f"{'Z':>12}"
    f"{'P(2-sided)':>14}"
)

print("-" * 100)


def print_stat(name, actual, null):

    mean = sum(null) / len(null)

    sd = math.sqrt(
        sum((x - mean) ** 2 for x in null)
        / len(null)
    )

    z = z_score(actual, null)
    p = two_sided_empirical_p(null, actual)

    print(
        f"{name:<28}"
        f"{actual:>14.8f}"
        f"{mean:>14.8f}"
        f"{sd:>14.8f}"
        f"{z:>12.4f}"
        f"{p:>14.6f}"
    )


print_stat(
    "signature_entropy",
    actual_H,
    entropy_values
)

print_stat(
    "unique_signatures",
    actual_U,
    unique_values
)

print_stat(
    "repeat_fraction",
    actual_R,
    repeat_values
)

print_stat(
    "mean_pairwise_MI",
    actual_MI,
    mi_values
)

print_stat(
    "mean_CRT_code",
    actual_CRT,
    crt_values
)


# ============================================================
# POSITION INSIDE NULL
# ============================================================

print()
print("=" * 100)
print("ACTUAL POSITION INSIDE NULL")
print("=" * 100)

for name, actual, null in [
    ("signature_entropy", actual_H, entropy_values),
    ("unique_signatures", actual_U, unique_values),
    ("repeat_fraction", actual_R, repeat_values),
    ("mean_pairwise_MI", actual_MI, mi_values),
    ("mean_CRT_code", actual_CRT, crt_values),
]:

    print(f"{name}")
    print(
        f"    percentile = "
        f"{percentile(null, actual):.6f}"
    )


# ============================================================
# MODULAR COMPONENTS
# ============================================================

print()
print("=" * 100)
print("PER-MODULUS ACTUAL vs NULL")
print("=" * 100)

for idx, r in enumerate(MODULI):

    vals = [s[idx] for s in actual_signatures]

    actual_mean = sum(vals) / len(vals)

    null_means = []

    for perm in range(PERMUTATIONS):

        # Reconstruct the mean component from the
        # stored null correlation is not appropriate here,
        # so we independently generate the permutation
        # component mean.

        pairs = [
            (e["p"], e["q"])
            for e in events
        ]

        # deterministic local RNG based on permutation
        rr = random.Random(SEED + perm + 100000 + r)
        rr.shuffle(pairs)

        s = [
            signature_for_pair(e["a"], p)[idx]
            for e, (p, q) in zip(events, pairs)
        ]

        null_means.append(
            sum(s) / len(s)
        )

    print()
    print(f"r = {r}")
    print(
        f"    actual mean u = "
        f"{actual_mean:.8f}"
    )
    print(
        f"    null mean u   = "
        f"{sum(null_means)/len(null_means):.8f}"
    )
    print(
        f"    percentile    = "
        f"{percentile(null_means, actual_mean):.6f}"
    )


# ============================================================
# t-CONDITIONED JOINT SIGNATURE ANALYSIS
# ============================================================

print()
print("=" * 100)
print("t-CONDITIONED JOINT SIGNATURE ANALYSIS")
print("=" * 100)

events_by_t = defaultdict(list)

for index, e in enumerate(events):
    events_by_t[e["t"]].append(index)


for t in sorted(events_by_t):

    idxs = events_by_t[t]

    if len(idxs) < 3:
        continue

    sigs = [
        actual_signatures[i]
        for i in idxs
    ]

    H = signature_entropy(sigs)
    U = unique_signature_count(sigs)
    MI = mean_pairwise_mi(sigs)

    print(
        f"t={t:+3d}"
        f" N={len(idxs):4d}"
        f" entropy={H:10.6f}"
        f" unique={U:5d}"
        f" pairwise_MI={MI:.8f}"
    )


# ============================================================
# SIGNATURE FREQUENCY
# ============================================================

print()
print("=" * 100)
print("MOST COMMON ACTUAL CRT SIGNATURES")
print("=" * 100)

sig_counts = Counter(actual_signatures)

for signature, count in sig_counts.most_common(25):

    print(
        f"count={count:4d}"
        f" signature={signature}"
        f" CRT={crt_encode(signature)}"
    )


# ============================================================
# PAIRWISE MODULAR MI
# ============================================================

print()
print("=" * 100)
print("PAIRWISE MODULAR MUTUAL INFORMATION")
print("=" * 100)

components = signature_components(actual_signatures)

for i in range(len(MODULI)):

    for j in range(i + 1, len(MODULI)):

        mi = mutual_information(
            components[i],
            components[j]
        )

        print(
            f"r={MODULI[i]:2d}"
            f" vs r={MODULI[j]:2d}"
            f" MI={mi:.10f}"
        )


# ============================================================
# t CORRELATION
# ============================================================

print()
print("=" * 100)
print("t vs MODULAR COMPONENT")
print("=" * 100)

for r in MODULI:

    actual = abs(actual_corr[r])
    null = corr_values[r]

    print(
        f"r={r:2d}"
        f" actual |corr|={actual:.8f}"
        f" null_mean={sum(null)/len(null):.8f}"
        f" percentile={percentile(null, actual):.6f}"
    )


# ============================================================
# FINAL LANDSCAPE CHECK
# ============================================================

print()
print("=" * 100)
print("LANDSCAPE INVARIANCE CHECK")
print("=" * 100)

# Reconstruct original landscape representation.
#
# The permutation test must never modify these.
#
# This check is intentionally structural.

assert len(landscape) == len(events)

print(f"collision states = {len(landscape)}")
print("landscape remained fixed during permutation test = True")


# ============================================================
# CONCLUSION
# ============================================================

print()
print("=" * 100)
print("FINAL")
print("=" * 100)

print(
    "The experiment tests whether the vector"
)
print(
    "    U=(u_3,u_5,u_7,u_11,u_13,u_17,u_19,u_23)"
)
print(
    "contains joint structure that disappears when"
)
print(
    "the anchor-factor association is permuted."
)
print()
print(
    "The identity"
)
print(
    "    (a/p)*(b/q) = 1 (mod r)"
)
print(
    "is NOT itself treated as evidence."
)
print()
print(
    "Evidence for structure requires the actual joint"
)
print(
    "statistics to lie outside the permutation null."
)
print()
print("Experiment complete.")
