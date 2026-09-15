#!/usr/bin/env python3
"""
====================================================================================================
SHIFTED SEMIPRIME RESIDUE-MATRIX / MATCHED-NONCOLLISION EXPERIMENT
====================================================================================================

Goal
----
Test whether the factorization of shifted values

    N_t = n + 2*M*t

has unusual modular residue geometry.

For every collision event

    N_t = a*b

we record

    a mod r
    b mod r
    N_t mod r

for several small prime moduli r.

The important comparison is NOT (p,q) vs (a,b).

Instead, for each observed collision N_t = a*b we construct a matched
null population of factor pairs drawn from the SAME factor range and
conditioned on:

    x*y = N_t

so that the comparison is about the residue structure of the factorization
of the SAME integer.

We also construct a second null:

    random semiprime pairs with the same factor-size distribution.

Statistics
----------
1. Residue-pair entropy
2. Number of distinct residue pairs
3. Concentration of residue pairs
4. Mutual information I(a mod r ; b mod r)
5. Chi-square distance from the expected multiplicative residue relation
6. Residue-pair enrichment
7. Cross-modulus joint-signature entropy

Important
---------
The identity

    (a*b) mod r = N_t mod r

is tautological.

The experiment therefore does NOT count that identity as evidence.

Evidence would require the observed residue distribution to differ from
matched factorizations / matched random semiprimes.

====================================================================================================
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass

# ================================================================================================
# CONFIG
# ================================================================================================

M = 111_546_435
T_MIN = -25
T_MAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
RANDOM_SEED = 1511464998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

# Number of matched residue-null samples per observed collision.
MATCHED_SAMPLES = 200

# Number of generic random semiprime samples for the global control.
RANDOM_SAMPLES = 200_000

# ================================================================================================
# PRIMES
# ================================================================================================

def sieve(limit: int) -> list[int]:
    bs = bytearray(b"\x01") * (limit + 1)
    bs[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if bs[p]:
            bs[p * p : limit + 1 : p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [i for i, v in enumerate(bs) if v]


PRIMES = sieve(FACTOR_MAX)

FACTOR_PRIMES = [
    p for p in PRIMES
    if FACTOR_MIN <= p <= FACTOR_MAX
]

FACTOR_PRIME_SET = set(FACTOR_PRIMES)

# ================================================================================================
# HELPERS
# ================================================================================================

def is_semiprime_from_factor_pair(a: int, b: int) -> bool:
    return (
        a in FACTOR_PRIME_SET
        and b in FACTOR_PRIME_SET
        and a * b > 0
    )


def factor_pair(n: int) -> tuple[int, int] | None:
    """
    Find the two prime factors when n is a semiprime whose factors lie in
    [FACTOR_MIN, FACTOR_MAX].

    This is deliberately simple because n is only ~1e8..1e10.
    """
    limit = math.isqrt(n)

    for p in FACTOR_PRIMES:
        if p > limit:
            break

        if n % p == 0:
            q = n // p

            if (
                q in FACTOR_PRIME_SET
                and p != q
            ):
                return (p, q)

    return None


def entropy(counter: Counter) -> float:
    total = sum(counter.values())

    if total == 0:
        return 0.0

    h = 0.0

    for c in counter.values():
        p = c / total
        h -= p * math.log2(p)

    return h


def mutual_information(
    pairs: list[tuple[int, int]]
) -> float:
    if not pairs:
        return 0.0

    xy = Counter(pairs)
    x = Counter(a for a, _ in pairs)
    y = Counter(b for _, b in pairs)

    n = len(pairs)
    mi = 0.0

    for (a, b), c in xy.items():
        pxy = c / n
        px = x[a] / n
        py = y[b] / n

        mi += pxy * math.log2(pxy / (px * py))

    return mi


def concentration(counter: Counter) -> float:
    """
    Simpson concentration:

        sum p_i^2

    Larger = more concentrated.
    """
    n = sum(counter.values())

    if n == 0:
        return 0.0

    return sum((c / n) ** 2 for c in counter.values())


def residue_signature(
    a: int,
    b: int,
    moduli: list[int],
) -> tuple[int, ...]:
    """
    Canonical unordered representation.

    Since a*b is symmetric, ordering the factor pair should not create
    artificial differences.
    """
    out = []

    for r in moduli:
        x = a % r
        y = b % r

        if x <= y:
            out.extend((x, y))
        else:
            out.extend((y, x))

    return tuple(out)


def multiplicative_residue_ok(
    a: int,
    b: int,
    n: int,
    r: int,
) -> bool:
    return (a * b) % r == n % r


def residue_pair_distribution(
    pairs: list[tuple[int, int]],
    r: int,
) -> Counter:
    c = Counter()

    for a, b in pairs:
        x = a % r
        y = b % r

        if x <= y:
            c[(x, y)] += 1
        else:
            c[(y, x)] += 1

    return c


# ================================================================================================
# GENERATE ACTUAL COLLISIONS
# ================================================================================================

def generate_actual_dataset(
    trials: int,
    seed: int,
) -> tuple[list[tuple[int, int, int, int, int]], list[int]]:
    """
    Returns:

        events:
            (trial, n, t, a, b)

        anchors:
            original n values

    The event generator is based entirely on shifted semiprimes.
    """
    rng = random.Random(seed)

    events = []
    anchors = []

    for trial in range(1, trials + 1):

        # Generate an actual semiprime anchor.
        p = rng.choice(FACTOR_PRIMES)
        q = rng.choice(FACTOR_PRIMES)

        while p == q:
            q = rng.choice(FACTOR_PRIMES)

        n = p * q
        anchors.append(n)

        for t in range(T_MIN, T_MAX + 1):
            if t == 0:
                continue

            Nt = n + 2 * M * t

            if Nt <= 0:
                continue

            pair = factor_pair(Nt)

            if pair is None:
                continue

            a, b = pair

            events.append((trial, n, t, a, b))

    return events, anchors


# ================================================================================================
# MATCHED NULL FOR THE SAME Nt
# ================================================================================================

def enumerate_factorizations_in_range(n: int) -> list[tuple[int, int]]:
    """
    Enumerate ALL prime factor pairs a*b=n with both factors inside the
    prescribed range.

    For a genuine semiprime this normally produces one unordered pair.

    The routine is kept explicit so the experiment is easy to audit.
    """
    out = []

    limit = math.isqrt(n)

    for p in FACTOR_PRIMES:
        if p > limit:
            break

        if n % p != 0:
            continue

        q = n // p

        if q not in FACTOR_PRIME_SET:
            continue

        if p <= q:
            out.append((p, q))

    return out


def generate_matched_null(
    actual_events: list[tuple[int, int, int, int, int]],
    samples_per_event: int,
    seed: int,
) -> list[tuple[int, int]]:
    """
    Null population based on the SAME Nt values.

    Each actual event contributes factor pairs from the same Nt.
    Because a semiprime normally has only one factorization, we also add
    randomized label/orientation draws to generate a legitimate null
    distribution over the factor roles.

    The null therefore tests whether the residue statistics can be
    explained entirely by the arithmetic of Nt itself.
    """
    rng = random.Random(seed)

    null_pairs = []

    for _, n, t, a, b in actual_events:

        Nt = n + 2 * M * t

        factorizations = enumerate_factorizations_in_range(Nt)

        if not factorizations:
            continue

        for _ in range(samples_per_event):

            x, y = rng.choice(factorizations)

            if rng.random() < 0.5:
                x, y = y, x

            null_pairs.append((x, y))

    return null_pairs


# ================================================================================================
# GENERIC RANDOM SEMIPRIME CONTROL
# ================================================================================================

def generate_random_semiprimes(
    count: int,
    seed: int,
) -> list[tuple[int, int]]:
    rng = random.Random(seed)

    out = []

    for _ in range(count):
        p = rng.choice(FACTOR_PRIMES)
        q = rng.choice(FACTOR_PRIMES)

        while p == q:
            q = rng.choice(FACTOR_PRIMES)

        out.append((p, q))

    return out


# ================================================================================================
# MODULAR ANALYSIS
# ================================================================================================

@dataclass
class ModStats:
    r: int
    entropy_actual: float
    entropy_matched: float
    entropy_random: float

    concentration_actual: float
    concentration_matched: float
    concentration_random: float

    mi_actual: float
    mi_matched: float
    mi_random: float

    observed_pairs: int
    random_pairs: int


def analyse_modulus(
    events: list[tuple[int, int, int, int, int]],
    matched: list[tuple[int, int]],
    random_pairs: list[tuple[int, int]],
    r: int,
) -> ModStats:

    actual_pairs = [
        (a, b)
        for _, _, _, a, b in events
    ]

    actual_dist = residue_pair_distribution(actual_pairs, r)
    matched_dist = residue_pair_distribution(matched, r)
    random_dist = residue_pair_distribution(random_pairs, r)

    return ModStats(
        r=r,

        entropy_actual=entropy(actual_dist),
        entropy_matched=entropy(matched_dist),
        entropy_random=entropy(random_dist),

        concentration_actual=concentration(actual_dist),
        concentration_matched=concentration(matched_dist),
        concentration_random=concentration(random_dist),

        mi_actual=mutual_information(
            [(a % r, b % r) for a, b in actual_pairs]
        ),
        mi_matched=mutual_information(
            [(a % r, b % r) for a, b in matched
        ]),
        mi_random=mutual_information(
            [(a % r, b % r) for a, b in random_pairs]
        ),

        observed_pairs=len(actual_pairs),
        random_pairs=len(random_pairs),
    )


# ================================================================================================
# JOINT SIGNATURE ANALYSIS
# ================================================================================================

def joint_signature_distribution(
    pairs: list[tuple[int, int]],
    moduli: list[int],
) -> Counter:

    return Counter(
        residue_signature(a, b, moduli)
        for a, b in pairs
    )


def print_joint_summary(
    name: str,
    pairs: list[tuple[int, int]],
    moduli: list[int],
):
    sigs = joint_signature_distribution(pairs, moduli)

    print(f"{name}")
    print(f"  observations        = {len(pairs):8d}")
    print(f"  unique signatures   = {len(sigs):8d}")
    print(f"  entropy             = {entropy(sigs):12.8f}")
    print(f"  concentration       = {concentration(sigs):12.8f}")


# ================================================================================================
# RESIDUE ENRICHMENT
# ================================================================================================

def compare_residue_enrichment(
    actual: list[tuple[int, int]],
    matched: list[tuple[int, int]],
    r: int,
    top_k: int = 15,
):
    a = residue_pair_distribution(actual, r)
    b = residue_pair_distribution(matched, r)

    na = sum(a.values())
    nb = sum(b.values())

    rows = []

    keys = set(a) | set(b)

    for key in keys:
        pa = a[key] / na if na else 0.0
        pb = b[key] / nb if nb else 0.0

        if pb == 0:
            ratio = float("inf")
        else:
            ratio = pa / pb

        rows.append(
            (
                ratio,
                key,
                a[key],
                b[key],
                pa,
                pb,
            )
        )

    rows.sort(reverse=True)

    print()
    print(f"r={r} TOP RESIDUE ENRICHMENT")
    print(
        "  residue-pair     ACT      NULL     "
        "ACT_FREQ       NULL_FREQ       RATIO"
    )

    for ratio, key, ca, cb, pa, pb in rows[:top_k]:
        print(
            f"  {str(key):14s} "
            f"{ca:7d} {cb:8d} "
            f"{pa:11.6f} {pb:13.6f} "
            f"{ratio:10.4f}"
        )


# ================================================================================================
# ACTUAL SHIFTED RESIDUES
# ================================================================================================

def shifted_residue_analysis(
    events: list[tuple[int, int, int, int, int]],
    moduli: list[int],
):
    print()
    print("=" * 100)
    print("SHIFTED N_t RESIDUE ANALYSIS")
    print("=" * 100)

    for r in moduli:

        counts = Counter()

        for _, n, t, a, b in events:
            Nt = n + 2 * M * t

            counts[(Nt % r, a % r, b % r)] += 1

        print()
        print(f"r={r}")

        # Verify arithmetic identity but do not treat it as evidence.
        failures = 0

        for (nr, ar, br), count in counts.items():
            if (ar * br) % r != nr:
                failures += count

        print(f"  residue identity failures = {failures}")

        # Show most common triples.
        for triple, count in counts.most_common(12):
            print(
                f"  Nt={triple[0]:2d} "
                f"a={triple[1]:2d} "
                f"b={triple[2]:2d} "
                f"count={count:4d}"
            )


# ================================================================================================
# MAIN
# ================================================================================================

def main():

    print("=" * 100)
    print("SHIFTED SEMIPRIME RESIDUE-MATRIX / MATCHED-NONCOLLISION EXPERIMENT")
    print("=" * 100)
    print(f"M                 = {M:,}")
    print(f"2M                = {2*M:,}")
    print(f"factor range      = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range            = [{T_MIN}, {T_MAX}]")
    print(f"trials              = {TRIALS}")
    print(f"matched/event      = {MATCHED_SAMPLES}")
    print(f"random controls    = {RANDOM_SAMPLES:,}")
    print(f"moduli             = {MODULI}")
    print(f"factor primes      = {len(FACTOR_PRIMES)}")
    print(f"random seed        = {RANDOM_SEED:,}")

    # --------------------------------------------------------------------------------------------
    # ACTUAL
    # --------------------------------------------------------------------------------------------

    events, anchors = generate_actual_dataset(
        TRIALS,
        RANDOM_SEED,
    )

    print()
    print("=" * 100)
    print("ACTUAL DATASET")
    print("=" * 100)
    print(f"anchors = {len(anchors)}")
    print(f"events  = {len(events)}")

    # --------------------------------------------------------------------------------------------
    # MATCHED NULL
    # --------------------------------------------------------------------------------------------

    matched = generate_matched_null(
        events,
        MATCHED_SAMPLES,
        RANDOM_SEED ^ 0xA5A5A5A5,
    )

    # --------------------------------------------------------------------------------------------
    # GENERIC RANDOM CONTROL
    # --------------------------------------------------------------------------------------------

    random_pairs = generate_random_semiprimes(
        RANDOM_SAMPLES,
        RANDOM_SEED ^ 0x5A5A5A5A,
    )

    # --------------------------------------------------------------------------------------------
    # MODULAR RESULTS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("PER-MODULUS RESULTS")
    print("=" * 100)

    print(
        "   r   "
        "ACT_H      NULL_H     RAND_H    "
        "ACT_C      NULL_C     RAND_C    "
        "ACT_MI     NULL_MI    RAND_MI"
    )

    stats = []

    for r in MODULI:

        s = analyse_modulus(
            events,
            matched,
            random_pairs,
            r,
        )

        stats.append(s)

        print(
            f"{r:4d} "
            f"{s.entropy_actual:10.6f} "
            f"{s.entropy_matched:10.6f} "
            f"{s.entropy_random:10.6f} "
            f"{s.concentration_actual:10.6f} "
            f"{s.concentration_matched:10.6f} "
            f"{s.concentration_random:10.6f} "
            f"{s.mi_actual:10.6f} "
            f"{s.mi_matched:10.6f} "
            f"{s.mi_random:10.6f}"
        )

    # --------------------------------------------------------------------------------------------
    # JOINT MODULAR SIGNATURE
    # --------------------------------------------------------------------------------------------

    actual_pairs = [
        (a, b)
        for _, _, _, a, b in events
    ]

    print()
    print("=" * 100)
    print("JOINT MODULAR SIGNATURE")
    print("=" * 100)

    print_joint_summary(
        "ACTUAL",
        actual_pairs,
        MODULI,
    )

    print_joint_summary(
        "MATCHED SAME-N_t",
        matched,
        MODULI,
    )

    print_joint_summary(
        "RANDOM SEMIPRIME",
        random_pairs,
        MODULI,
    )

    # --------------------------------------------------------------------------------------------
    # RESIDUE ENRICHMENT
    # --------------------------------------------------------------------------------------------

    for r in MODULI:
        compare_residue_enrichment(
            actual_pairs,
            matched,
            r,
            top_k=10,
        )

    # --------------------------------------------------------------------------------------------
    # SHIFTED N_t ANALYSIS
    # --------------------------------------------------------------------------------------------

    shifted_residue_analysis(
        events,
        MODULI,
    )

    # --------------------------------------------------------------------------------------------
    # CROSS-MODULUS SIGNATURE
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CROSS-MODULUS SIGNATURE REPETITION")
    print("=" * 100)

    actual_sig = joint_signature_distribution(
        actual_pairs,
        MODULI,
    )

    matched_sig = joint_signature_distribution(
        matched,
        MODULI,
    )

    repeated_actual = [
        (sig, count)
        for sig, count in actual_sig.items()
        if count > 1
    ]

    repeated_matched = [
        (sig, count)
        for sig, count in matched_sig.items()
        if count > 1
    ]

    print(
        f"actual repeated signatures = "
        f"{len(repeated_actual)}"
    )

    print(
        f"matched repeated signatures = "
        f"{len(repeated_matched)}"
    )

    if repeated_actual:
        print()
        print("Most repeated ACTUAL signatures:")
        for sig, count in sorted(
            repeated_actual,
            key=lambda x: x[1],
            reverse=True
        )[:20]:
            print(
                f"  count={count:4d} "
                f"signature={sig}"
            )

    # --------------------------------------------------------------------------------------------
    # IMPORTANT CONTROL:
    # TEST FOR FACTOR-SIZE CONFUSION
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FACTOR-SIZE NORMALIZATION")
    print("=" * 100)

    actual_norm = []
    matched_norm = []

    for _, _, _, a, b in events:
        actual_norm.append(
            (
                a / FACTOR_MAX,
                b / FACTOR_MAX,
            )
        )

    for a, b in matched:
        matched_norm.append(
            (
                a / FACTOR_MAX,
                b / FACTOR_MAX,
            )
        )

    mean_actual_a = (
        sum(x for x, _ in actual_norm)
        / len(actual_norm)
    )

    mean_actual_b = (
        sum(y for _, y in actual_norm)
        / len(actual_norm)
    )

    mean_matched_a = (
        sum(x for x, _ in matched_norm)
        / len(matched_norm)
    )

    mean_matched_b = (
        sum(y for _, y in matched_norm)
        / len(matched_norm)
    )

    print(
        f"actual mean normalized factor 1 = "
        f"{mean_actual_a:.8f}"
    )

    print(
        f"actual mean normalized factor 2 = "
        f"{mean_actual_b:.8f}"
    )

    print(
        f"matched mean normalized factor 1 = "
        f"{mean_matched_a:.8f}"
    )

    print(
        f"matched mean normalized factor 2 = "
        f"{mean_matched_b:.8f}"
    )

    # --------------------------------------------------------------------------------------------
    # FINAL
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FINAL")
    print("=" * 100)

    print(
        "The arithmetic identity "
        "(a*b) mod r = N_t mod r was verified, "
        "but is not counted as evidence."
    )

    print()
    print(
        "The primary question is whether the observed residue distribution "
        "of the shifted semiprime factors differs from:"
    )

    print(
        "  1. factorization-conditioned matched controls, and"
    )

    print(
        "  2. generic random semiprimes from the same factor range."
    )

    print()
    print(
        "A strong result would require the ACTUAL residue statistics "
        "to remain anomalous against the SAME-N_t matched null."
    )

    print()
    print("Experiment complete.")


if __name__ == "__main__":
    main()
