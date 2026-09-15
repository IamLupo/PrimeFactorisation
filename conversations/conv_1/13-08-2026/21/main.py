#!/usr/bin/env python3

import math
import random
from collections import Counter

# ============================================================
# KAPPA EXPERIMENT 21
# PROGRESSIVE MODULUS DISCRIMINATION / SIGNATURE COLLISION DEPTH
#
# NO CSV OUTPUT
# ============================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

TARGETS = 12

# Deliberately smaller than Experiment 20 because we test
# many prefix lengths.
WRONG_PAIRS = 20_000

# Number of initial moduli used in each progressive test.
PREFIXES = [3, 5, 7, 9, 11, 13, 15]

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000


# ============================================================
# BASIC NUMBER THEORY
# ============================================================

def is_prime(n):
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def random_prime():
    while True:
        n = random.randrange(PRIME_LOW, PRIME_HIGH)

        if n % 2 == 0:
            n += 1

        if is_prime(n):
            return n


def factor_pair():
    while True:
        p = random_prime()
        q = random_prime()

        if p != q:
            return p, q


# ============================================================
# MODULUS INVENTORY
#
# EXACTLY THE EXPERIMENT-20 FORMULA
# ============================================================

MODULI = [r * r + 3 for r in R_VALUES]


def units_mod(m):
    return [
        x for x in range(1, m)
        if math.gcd(x, m) == 1
    ]


def multiplicative_order(a, m):
    if math.gcd(a, m) != 1:
        return None

    x = 1

    for order in range(1, m + 1):
        x = (x * a) % m

        if x == 1:
            return order

    raise RuntimeError(
        f"multiplicative order failed for a={a}, m={m}"
    )


def build_inventory():

    inventory = {}

    for m in MODULI:

        units = units_mod(m)

        orders = {
            u: multiplicative_order(u, m)
            for u in units
        }

        inventory[m] = {
            "units": units,
            "orbit": len(units),
            "orders": orders,
        }

    return inventory


INVENTORY = build_inventory()


# ============================================================
# LOCAL SIGNATURE
# ============================================================

def local_signature(p, q, m):

    a = p % m
    b = q % m

    oa = INVENTORY[m]["orders"].get(a)
    ob = INVENTORY[m]["orders"].get(b)

    return (
        (a, b),                         # Fpair
        tuple(sorted((a, b))),          # Fsorted
        (a + b) % m,                    # Fsum
        (a * b) % m,                    # Fprod
        (a - b) % m,                    # Fdiff
        (a ** 3 + b ** 3) % m,          # cube_sum
        (oa, ob),                       # order_pair
        tuple(sorted((oa, ob))),        # order_sorted
    )


SIGNATURES = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "order_pair",
    "order_sorted",
]

SIG_INDEX = {
    name: i
    for i, name in enumerate(SIGNATURES)
}


# ============================================================
# COMPLETE LOCAL DATA
#
# Computed ONCE per pair.
# This is the main speed improvement over Experiment 20.
# ============================================================

def pair_signature_data(p, q):

    result = [[] for _ in SIGNATURES]

    for m in MODULI:

        local = local_signature(p, q, m)

        for i in range(len(SIGNATURES)):
            result[i].append(local[i])

    return tuple(
        tuple(values)
        for values in result
    )


# ============================================================
# PREFIX SIGNATURE
# ============================================================

def prefix_signature(data, signature_index, prefix):

    return data[signature_index][:prefix]


# ============================================================
# PROGRESSIVE COLLISION TEST
# ============================================================

def progressive_collision_test(
    p,
    q,
    wrong_pairs,
    prefixes
):

    actual = pair_signature_data(p, q)

    # collisions[signature][prefix]
    collisions = {
        name: {
            prefix: 0
            for prefix in prefixes
        }
        for name in SIGNATURES
    }

    # We generate each wrong pair only ONCE.
    # Its complete signature data is then reused for every prefix
    # and every signature.
    for _ in range(wrong_pairs):

        a, b = factor_pair()

        candidate = pair_signature_data(a, b)

        for name in SIGNATURES:

            i = SIG_INDEX[name]

            actual_values = actual[i]
            candidate_values = candidate[i]

            for prefix in prefixes:

                if (
                    actual_values[:prefix]
                    == candidate_values[:prefix]
                ):
                    collisions[name][prefix] += 1

    return collisions


# ============================================================
# FPAIR THEORETICAL BASELINE
# ============================================================

def fpair_prefix_baseline(prefix):

    return -sum(
        math.log10(
            INVENTORY[MODULI[i]]["orbit"]
        )
        for i in range(prefix)
    )


# ============================================================
# FORMATTING
# ============================================================

def fmt_rate(hits, trials):

    if hits == 0:
        return f"<1/{trials:,}"

    return f"{hits / trials:.8g}"


# ============================================================
# MAIN
# ============================================================

print("=" * 78)
print("KAPPA EXPERIMENT 21")
print("PROGRESSIVE MODULUS DISCRIMINATION / SIGNATURE COLLISION DEPTH")
print("NO CSV OUTPUT")
print("=" * 78)

print(f"random seed   = {SEED}")
print(f"R values      = {R_VALUES}")
print(f"targets       = {TARGETS}")
print(f"wrong pairs   = {WRONG_PAIRS:,}")
print(f"prefixes      = {PREFIXES}")
print()


# ============================================================
# MODULUS INVENTORY
# ============================================================

print("=" * 78)
print("1. MODULUS INVENTORY")
print("=" * 78)

for r, m in zip(R_VALUES, MODULI):

    orbit = INVENTORY[m]["orbit"]

    print(
        f"r={r:3d} "
        f"m={m:6d} "
        f"units={orbit:7d}"
    )

print()


# ============================================================
# THEORETICAL FPAIR PREFIX BASELINES
# ============================================================

print("=" * 78)
print("2. THEORETICAL FPAIR PREFIX BASELINES")
print("=" * 78)

print(
    f"{'prefix':>8s}"
    f"{'last m':>10s}"
    f"{'log10 generic freq':>24s}"
)

for prefix in PREFIXES:

    m = MODULI[prefix - 1]

    baseline = fpair_prefix_baseline(prefix)

    print(
        f"{prefix:8d}"
        f"{m:10d}"
        f"{baseline:24.6f}"
    )

print()


# ============================================================
# AGGREGATE STORAGE
# ============================================================

aggregate_collisions = {
    name: {
        prefix: []
        for prefix in PREFIXES
    }
    for name in SIGNATURES
}


# ============================================================
# TARGET LOOP
# ============================================================

for target in range(1, TARGETS + 1):

    p, q = factor_pair()

    n = p * q
    s = p + q
    delta = (p - q) ** 2

    print("=" * 78)
    print(f"TARGET {target}")
    print("=" * 78)

    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"s       = {s}")
    print(f"Delta   = {delta}")
    print(f"n bits  = {n.bit_length()}")
    print()

    collisions = progressive_collision_test(
        p,
        q,
        WRONG_PAIRS,
        PREFIXES
    )

    # --------------------------------------------------------
    # Collision depth table
    # --------------------------------------------------------

    print("PROGRESSIVE COMPLETE-SIGNATURE COLLISION RATES")
    print("-" * 78)

    print(
        f"{'prefix':>7s}"
        f"{'m':>8s}"
        + "".join(
            f"{name:>17s}"
            for name in SIGNATURES
        )
    )

    for prefix in PREFIXES:

        m = MODULI[prefix - 1]

        row = [
            f"{prefix:7d}",
            f"{m:8d}",
        ]

        for name in SIGNATURES:

            hits = collisions[name][prefix]

            aggregate_collisions[name][prefix].append(
                hits
            )

            row.append(
                f"{fmt_rate(hits, WRONG_PAIRS):>17s}"
            )

        print("".join(row))

    print()

    # --------------------------------------------------------
    # First-zero depth
    # --------------------------------------------------------

    print("FIRST OBSERVED ZERO-COLLISION PREFIX")
    print("-" * 78)

    for name in SIGNATURES:

        zero_prefix = None

        for prefix in PREFIXES:

            if collisions[name][prefix] == 0:
                zero_prefix = prefix
                break

        if zero_prefix is None:
            text = "none"
        else:
            text = str(zero_prefix)

        print(
            f"{name:18s} {text}"
        )

    print()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print()
print("MEDIAN COLLISION COUNTS BY PREFIX")
print("-" * 78)

print(
    f"{'prefix':>7s}"
    + "".join(
        f"{name:>17s}"
        for name in SIGNATURES
    )
)

for prefix in PREFIXES:

    row = [f"{prefix:7d}"]

    for name in SIGNATURES:

        values = sorted(
            aggregate_collisions[name][prefix]
        )

        if len(values) % 2:

            median = values[len(values) // 2]

        else:

            median = (
                values[len(values) // 2 - 1]
                + values[len(values) // 2]
            ) / 2

        row.append(
            f"{median:17.1f}"
        )

    print("".join(row))


# ============================================================
# ZERO-COLLISION TARGET COUNT
# ============================================================

print()
print("TARGETS WITH ZERO OBSERVED COLLISIONS")
print("-" * 78)

print(
    f"{'signature':18s}"
    + "".join(
        f"{prefix:>10d}"
        for prefix in PREFIXES
    )
)

for name in SIGNATURES:

    row = [f"{name:18s}"]

    for prefix in PREFIXES:

        values = aggregate_collisions[name][prefix]

        zero_count = sum(
            1
            for x in values
            if x == 0
        )

        row.append(
            f"{zero_count:10d}"
        )

    print("".join(row))


# ============================================================
# FPAIR OBSERVED VS THEORETICAL
# ============================================================

print()
print("FPAIR: OBSERVED MEDIAN VS THEORETICAL BASELINE")
print("-" * 78)

print(
    f"{'prefix':>8s}"
    f"{'theory log10':>18s}"
    f"{'median empirical':>20s}"
)

for prefix in PREFIXES:

    theory = fpair_prefix_baseline(prefix)

    values = sorted(
        aggregate_collisions["Fpair"][prefix]
    )

    if len(values) % 2:

        median_hits = values[len(values) // 2]

    else:

        median_hits = (
            values[len(values) // 2 - 1]
            + values[len(values) // 2]
        ) / 2

    if median_hits == 0:

        empirical = (
            f"<1/{WRONG_PAIRS:,}"
        )

    else:

        empirical_rate = median_hits / WRONG_PAIRS
        empirical = f"{math.log10(empirical_rate):.6f}"

    print(
        f"{prefix:8d}"
        f"{theory:18.6f}"
        f"{empirical:>20s}"
    )


# ============================================================
# WIN COUNTS
# ============================================================

print()
print("DISCRIMINATION WIN COUNTS")
print("-" * 78)

print(
    "For each prefix, the winner is the signature with the "
    "fewest total collisions across all targets."
)
print()

for prefix in PREFIXES:

    totals = {}

    for name in SIGNATURES:

        totals[name] = sum(
            aggregate_collisions[name][prefix]
        )

    best = min(
        totals.values()
    )

    winners = [
        name
        for name in SIGNATURES
        if totals[name] == best
    ]

    print(
        f"prefix={prefix:2d} "
        f"winner={', '.join(winners):20s} "
        f"total collisions={best}"
    )


# ============================================================
# INTERPRETATION
# ============================================================

print()
print("=" * 78)
print("INTERPRETATION")
print("=" * 78)

print("""
Experiment 21 measures COLLISION DEPTH.

Instead of asking only whether the complete 15-modulus signature
collides, the experiment asks:

    How many moduli are required before wrong factor pairs
    become statistically indistinguishable from the target?

For every signature, prefixes of length:

    3, 5, 7, 9, 11, 13, 15

are tested.

Fpair is the strongest reference because its generic probability
has the direct theoretical baseline:

    product(1 / unit-orbit-size)

over the selected prefix.

The important quantity is therefore the transition from:

    many collisions
        ->
    few collisions
        ->
    zero observed collisions.

If Fpair reaches zero collisions substantially earlier than the
compressed signatures, that supports the hypothesis that retaining
the complete ordered residue pair gives substantially more
discriminating information than scalar invariants.

Fsorted provides the ordering-control experiment.

Fsum/Fprod/Fdiff/cube_sum test compressed algebraic information.

order_pair/order_sorted test multiplicative-order information.

IMPORTANT:

Zero observed collisions does NOT prove mathematical uniqueness.
With 20,000 random wrong pairs, zero collisions only establishes
an empirical upper bound on the collision probability.

The particularly useful result is the PREFIX at which each
signature consistently reaches zero collisions across targets.
""")


print()
print("=" * 78)
print("DONE")
print("=" * 78)