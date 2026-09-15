#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# ============================================================
# CRT-HYPERBOLA / MOD-M FACTOR-SELECTION EXPERIMENT
# ============================================================

M = 111_546_435
T_MIN = -25
T_MAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
RANDOM_NULL_SAMPLES = 200_000
HYPERBOLA_SAMPLES = 200_000

SEED = 1_511_464_998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]


# ============================================================
# PRIME GENERATION
# ============================================================

def sieve(n: int):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(n)) + 1):
        if sieve[p]:
            sieve[p * p:n + 1:p] = b"\x00" * (
                ((n - p * p) // p) + 1
            )

    return [i for i in range(2, n + 1) if sieve[i]]


ALL_PRIMES = [
    p for p in sieve(FACTOR_MAX)
    if p >= FACTOR_MIN
]

PRIME_SET = set(ALL_PRIMES)


# ============================================================
# BASIC HELPERS
# ============================================================

def pair_key(a: int, b: int):
    return (min(a, b), max(a, b))


def residue_vector(x: int):
    return tuple(x % r for r in MODULI)


def modular_inverse(a: int, m: int):
    return pow(a, -1, m)


def normalize_factor(x: int):
    return (x - FACTOR_MIN) / (FACTOR_MAX - FACTOR_MIN)


def mutual_information(xs, ys):
    n = len(xs)

    if n == 0:
        return 0.0

    joint = Counter(zip(xs, ys))
    cx = Counter(xs)
    cy = Counter(ys)

    mi = 0.0

    for (x, y), count in joint.items():
        pxy = count / n
        px = cx[x] / n
        py = cy[y] / n

        mi += pxy * math.log2(pxy / (px * py))

    return mi


def entropy(values):
    n = len(values)

    if n == 0:
        return 0.0

    c = Counter(values)

    return -sum(
        (v / n) * math.log2(v / n)
        for v in c.values()
    )


# ============================================================
# GENERATE ACTUAL COLLISION EVENTS
# ============================================================

def generate_actual_events():

    events = []

    # We deliberately generate new semiprimes until we have
    # exactly TRIALS usable anchors.

    while len(events) < TRIALS:

        p = random.choice(ALL_PRIMES)
        q = random.choice(ALL_PRIMES)

        if p == q:
            continue

        n = p * q

        for t in range(T_MIN, T_MAX + 1):

            if t == 0:
                continue

            Nt = n + 2 * M * t

            # We only care about positive Nt.
            if Nt <= 0:
                continue

            # Find all factor pairs in our prime interval.
            found = None

            for a in ALL_PRIMES:

                if a * a > Nt:
                    break

                if Nt % a != 0:
                    continue

                b = Nt // a

                if b in PRIME_SET and FACTOR_MIN <= b <= FACTOR_MAX:
                    found = pair_key(a, b)
                    break

            if found is None:
                continue

            a, b = found

            events.append({
                "p": min(p, q),
                "q": max(p, q),
                "n": n,
                "t": t,
                "Nt": Nt,
                "a": a,
                "b": b,
            })

    return events


# ============================================================
# NOTE:
# In your previous experiments the anchor set was already
# generated externally. If you want exactly the same 300
# anchors as before, replace generate_actual_events() with
# your existing anchor loader.
#
# The analysis below is independent of how the anchors were
# obtained.
# ============================================================


# ============================================================
# NULL 1
# RANDOM PRIME PAIRS
# ============================================================

def generate_random_prime_pairs(count):

    result = []

    for _ in range(count):

        a = random.choice(ALL_PRIMES)
        b = random.choice(ALL_PRIMES)

        if a == b:
            continue

        result.append(pair_key(a, b))

    return result


# ============================================================
# NULL 2
# MOD-M HYPERBOLA NULL
#
# Given C = p*q mod M:
#
#       a*b = C mod M
#
# Therefore:
#
#       b = C * a^(-1) mod M
#
# We randomly choose prime a and ask whether the modularly
# determined b is ALSO a prime in the allowed range.
# ============================================================

def generate_hyperbola_pairs(events, target_count):

    result = []

    attempts = 0
    max_attempts = target_count * 500

    while len(result) < target_count and attempts < max_attempts:

        attempts += 1

        event = random.choice(events)

        p = event["p"]
        q = event["q"]

        C = (p * q) % M

        a = random.choice(ALL_PRIMES)

        # Since gcd(a, M)=1 for the primes in this range,
        # inverse exists.
        if math.gcd(a, M) != 1:
            continue

        inv_a = modular_inverse(a, M)

        b_residue = (C * inv_a) % M

        # b must be an actual integer in our allowed interval.
        # Since M > FACTOR_MAX, there is at most one candidate.
        if not (FACTOR_MIN <= b_residue <= FACTOR_MAX):
            continue

        if b_residue not in PRIME_SET:
            continue

        if b_residue == a:
            continue

        result.append(pair_key(a, b_residue))

    return result


# ============================================================
# NULL 3
# CONDITIONED SAME Nt
#
# For each actual event, its Nt is kept fixed.
# The factor pair is therefore the unique prime factor pair
# for that Nt.
#
# This is useful as a sanity check and MUST reproduce the
# actual residues.
# ============================================================

def build_same_Nt_null(events):

    result = []

    for e in events:
        result.append((e["a"], e["b"]))

    return result


# ============================================================
# STATISTICS
# ============================================================

def analyze_dataset(name, pairs):

    if not pairs:
        print(name, "EMPTY")
        return None

    a_values = [a for a, b in pairs]
    b_values = [b for a, b in pairs]

    stats = {}

    stats["mean_a_norm"] = sum(
        normalize_factor(a) for a in a_values
    ) / len(a_values)

    stats["mean_b_norm"] = sum(
        normalize_factor(b) for b in b_values
    ) / len(b_values)

    stats["entropy_a"] = entropy(
        residue_vector(a) for a in a_values
    )

    stats["entropy_b"] = entropy(
        residue_vector(b) for b in b_values
    )

    for r in MODULI:

        xr = [a % r for a in a_values]
        yr = [b % r for b in b_values]

        stats[f"mi_{r}"] = mutual_information(xr, yr)

        stats[f"joint_entropy_{r}"] = entropy(
            list(zip(xr, yr))
        )

    # Full CRT signature
    signatures = [
        (residue_vector(a), residue_vector(b))
        for a, b in pairs
    ]

    stats["full_entropy"] = entropy(signatures)
    stats["full_unique"] = len(set(signatures))
    stats["full_concentration"] = (
        max(Counter(signatures).values()) / len(signatures)
    )

    # modular hyperbola residual
    failures = 0

    for a, b in pairs:

        # This is the basic relation tested by the experiment.
        # For arbitrary random pairs this will generally fail.
        pass

    return stats


# ============================================================
# PER-MODULUS HYPERBOLA OCCUPANCY
# ============================================================

def hyperbola_occupancy(events):

    output = {}

    for r in MODULI:

        actual = Counter()
        expected = Counter()

        for e in events:

            C = (e["p"] * e["q"]) % r

            a = e["a"] % r
            b = e["b"] % r

            actual[(a, b)] += 1

            expected[
                (a, (C * pow(a, -1, r)) % r)
            ] += 1

        output[r] = {
            "actual": actual,
            "expected": expected,
        }

    return output


# ============================================================
# GLOBAL CRT SIGNATURE
# ============================================================

def print_crt_table(name, pairs):

    print()
    print("=" * 100)
    print(name)
    print("=" * 100)

    print(
        f"{'r':>4} "
        f"{'MI':>14} "
        f"{'JointH':>14} "
        f"{'UniquePairs':>14}"
    )

    for r in MODULI:

        xs = [a % r for a, b in pairs]
        ys = [b % r for a, b in pairs]

        joint = list(zip(xs, ys))

        print(
            f"{r:4d} "
            f"{mutual_information(xs, ys):14.8f} "
            f"{entropy(joint):14.8f} "
            f"{len(set(joint)):14d}"
        )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(SEED)

    print("=" * 100)
    print("CRT-HYPERBOLA / MOD-M FACTOR-SELECTION EXPERIMENT")
    print("=" * 100)

    print(f"M                 = {M:,}")
    print(f"factor range      = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range           = [{T_MIN}, {T_MAX}]")
    print(f"actual anchors    = {TRIALS:,}")
    print(f"random null       = {RANDOM_NULL_SAMPLES:,}")
    print(f"hyperbola null    = {HYPERBOLA_SAMPLES:,}")
    print(f"seed              = {SEED}")
    print()

    # --------------------------------------------------------
    # ACTUAL
    # --------------------------------------------------------

    actual_events = generate_actual_events()

    actual_pairs = [
        (e["a"], e["b"])
        for e in actual_events
    ]

    print(
        f"actual collision events = {len(actual_events):,}"
    )

    # --------------------------------------------------------
    # NULL 1
    # --------------------------------------------------------

    random_pairs = generate_random_prime_pairs(
        RANDOM_NULL_SAMPLES
    )

    print(
        f"random prime pairs      = {len(random_pairs):,}"
    )

    # --------------------------------------------------------
    # NULL 2
    # --------------------------------------------------------

    hyperbola_pairs = generate_hyperbola_pairs(
        actual_events,
        HYPERBOLA_SAMPLES
    )

    print(
        f"hyperbola pairs         = {len(hyperbola_pairs):,}"
    )

    # --------------------------------------------------------
    # NULL 3
    # --------------------------------------------------------

    same_nt_pairs = build_same_Nt_null(
        actual_events
    )

    # --------------------------------------------------------
    # STATISTICS
    # --------------------------------------------------------

    stats_actual = analyze_dataset(
        "ACTUAL",
        actual_pairs
    )

    stats_random = analyze_dataset(
        "RANDOM",
        random_pairs
    )

    stats_hyper = analyze_dataset(
        "HYPERBOLA",
        hyperbola_pairs
    )

    stats_same = analyze_dataset(
        "SAME-Nt",
        same_nt_pairs
    )

    # --------------------------------------------------------
    # PER-MODULUS
    # --------------------------------------------------------

    print_crt_table(
        "PER-MODULUS RESIDUE STRUCTURE",
        actual_pairs
    )

    print_crt_table(
        "RANDOM PRIME-PAIR NULL",
        random_pairs
    )

    print_crt_table(
        "MOD-M HYPERBOLA NULL",
        hyperbola_pairs
    )

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL vs NULL")
    print("=" * 100)

    print(
        f"{'STATISTIC':<28}"
        f"{'ACTUAL':>16}"
        f"{'RANDOM':>16}"
        f"{'HYPERBOLA':>16}"
    )

    print("-" * 100)

    keys = [
        "mean_a_norm",
        "mean_b_norm",
        "entropy_a",
        "entropy_b",
        "full_entropy",
        "full_unique",
        "full_concentration",
    ]

    for key in keys:

        print(
            f"{key:<28}"
            f"{stats_actual[key]:16.8f}"
            f"{stats_random[key]:16.8f}"
            f"{stats_hyper[key]:16.8f}"
        )

    print()
    print("=" * 100)
    print("MUTUAL INFORMATION")
    print("=" * 100)

    print(
        f"{'r':>4}"
        f"{'ACTUAL':>16}"
        f"{'RANDOM':>16}"
        f"{'HYPERBOLA':>16}"
    )

    for r in MODULI:

        print(
            f"{r:4d}"
            f"{stats_actual[f'mi_{r}']:16.8f}"
            f"{stats_random[f'mi_{r}']:16.8f}"
            f"{stats_hyper[f'mi_{r}']:16.8f}"
        )

    # --------------------------------------------------------
    # EXACT MODULAR RELATION
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL MOD-M IDENTITY")
    print("=" * 100)

    failures = 0

    for e in actual_events:

        a = e["a"]
        b = e["b"]
        C = e["p"] * e["q"]

        if (a * b - C) % M != 0:
            failures += 1

    print(
        f"(a*b - p*q) mod M failures = {failures}"
    )

    # --------------------------------------------------------
    # HYPERBOLA CONDITION
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("HYPERBOLA NULL CONDITION")
    print("=" * 100)

    hyper_failures = 0

    for a, b in hyperbola_pairs:

        # The hyperbola null is constructed using a selected
        # event's C. Because the pair itself does not retain C,
        # this global check is intentionally not performed here.
        #
        # The important comparison is the distribution.
        pass

    print(
        "NULL 2 explicitly enforces b = C*a^(-1) mod M"
    )

    # --------------------------------------------------------
    # T / FACTOR SIZE
    # --------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL FACTOR GEOMETRY BY t")
    print("=" * 100)

    grouped = defaultdict(list)

    for e in actual_events:
        grouped[e["t"]].append(e)

    for t in sorted(grouped):

        events = grouped[t]

        ndp = []
        ndq = []

        for e in events:

            p = e["p"]
            q = e["q"]
            a = e["a"]
            b = e["b"]

            # Keep factor ordering aligned by magnitude.
            if p <= q:
                dp = a - p
                dq = b - q
            else:
                dp = b - q
                dq = a - p

            ndp.append(dp / 100000.0)
            ndq.append(dq / 100000.0)

        print(
            f"t={t:+3d}"
            f" count={len(events):4d}"
            f" mean_ndp={sum(ndp)/len(ndp):+.8f}"
            f" mean_ndq={sum(ndq)/len(ndq):+.8f}"
        )

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print("""
The critical comparison is:

    ACTUAL
        observed factors of N_t = p*q + 2*M*t

    RANDOM
        unrestricted random prime pairs

    HYPERBOLA
        random prime a, with
        b = (p*q) * a^(-1) mod M,
        requiring b itself to be an allowed prime

The HYPERBOLA null already satisfies the exact modular
constraint. Therefore:

    ACTUAL ~= HYPERBOLA

would mean that the observed modular structure is largely
explained by the congruence plus prime-range selection.

But:

    ACTUAL differs significantly from HYPERBOLA

would be evidence for an additional structure beyond the
obvious modulo-M identity.

The most important outputs are therefore the differences
in MI for r=11,13,17,19,23 and the full CRT signature
statistics.
""")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
