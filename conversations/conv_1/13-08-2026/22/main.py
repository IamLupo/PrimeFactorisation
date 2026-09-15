#!/usr/bin/env python3

import math
import random
from collections import defaultdict

# ============================================================
# KAPPA EXPERIMENT 22
# COLLISION THRESHOLD / ADVERSARIAL COLLISION SEARCH
#
# NO CSV OUTPUT
# ============================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

TARGETS = 12

# Random wrong pairs used to establish the baseline.
WRONG_PAIRS = 20_000

# Number of additional candidates searched specifically for
# collisions at the critical prefixes.
ADVERSARIAL_PAIRS = 250_000

# Prefixes around the observed transition.
PREFIXES = [3, 4, 5, 6, 7]

# Prime generation interval.
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


def random_wrong_pair(p, q):
    while True:
        a = random_prime()
        b = random_prime()

        if a == b:
            continue

        if {a, b} != {p, q}:
            return a, b


# ============================================================
# MODULUS / UNIT INVENTORY
# ============================================================

MODULI = [r * r + 3 for r in R_VALUES]


def units_mod(m):
    return [
        x
        for x in range(1, m)
        if math.gcd(x, m) == 1
    ]


def multiplicative_order(a, m):
    if math.gcd(a, m) != 1:
        return None

    x = 1
    order = 0

    while True:
        order += 1
        x = (x * a) % m

        if x == 1:
            return order


def build_inventory():
    inventory = {}

    for m in MODULI:

        units = units_mod(m)

        orders = {}

        for u in units:
            orders[u] = multiplicative_order(u, m)

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

    return {
        "pair": (a, b),

        "sorted": tuple(
            sorted((a, b))
        ),

        "sum": (a + b) % m,

        "prod": (a * b) % m,

        "diff": (a - b) % m,

        "order_pair": (oa, ob),

        "order_sorted": tuple(
            sorted((oa, ob))
        ),

        "cube_sum":
            (a ** 3 + b ** 3) % m,
    }


# ============================================================
# SIGNATURE TYPES
# ============================================================

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


# ============================================================
# PREFIX SIGNATURE
# ============================================================

def prefix_signature(
    p,
    q,
    signature_name,
    prefix
):

    out = []

    for i in range(prefix):

        m = MODULI[i]

        s = local_signature(
            p,
            q,
            m
        )

        if signature_name == "Fpair":
            value = s["pair"]

        elif signature_name == "Fsorted":
            value = s["sorted"]

        elif signature_name == "Fsum":
            value = s["sum"]

        elif signature_name == "Fprod":
            value = s["prod"]

        elif signature_name == "Fdiff":
            value = s["diff"]

        elif signature_name == "cube_sum":
            value = s["cube_sum"]

        elif signature_name == "order_pair":
            value = s["order_pair"]

        elif signature_name == "order_sorted":
            value = s["order_sorted"]

        else:
            raise ValueError(signature_name)

        out.append(value)

    return tuple(out)


# ============================================================
# COMPLETE PREFIX SIGNATURE
# ============================================================

def all_prefix_signatures(
    p,
    q,
    prefix
):

    result = {}

    for name in SIGNATURES:

        result[name] = prefix_signature(
            p,
            q,
            name,
            prefix
        )

    return result


# ============================================================
# GENERIC FPAIR BASELINE
# ============================================================

def theoretical_log10(prefix):

    return -sum(
        math.log10(
            INVENTORY[
                MODULI[i]
            ]["orbit"]
        )
        for i in range(prefix)
    )


# ============================================================
# RANDOM COLLISION BASELINE
# ============================================================

def random_collision_test(
    p,
    q,
    prefix
):

    actual = all_prefix_signatures(
        p,
        q,
        prefix
    )

    hits = {
        name: 0
        for name in SIGNATURES
    }

    for _ in range(WRONG_PAIRS):

        a, b = random_wrong_pair(
            p,
            q
        )

        candidate = all_prefix_signatures(
            a,
            b,
            prefix
        )

        for name in SIGNATURES:

            if candidate[name] == actual[name]:
                hits[name] += 1

    return hits


# ============================================================
# ADVERSARIAL COLLISION SEARCH
# ============================================================

def adversarial_search(
    p,
    q,
    prefix
):
    """
    Search specifically for actual wrong pairs that collide.

    Unlike the random experiment, this records the first
    concrete collision and the corresponding factor pair.

    This is especially important at prefixes 5 and 7,
    where Experiment 21 observed mostly zero collisions.
    """

    actual = all_prefix_signatures(
        p,
        q,
        prefix
    )

    found = {
        name: []
        for name in SIGNATURES
    }

    for _ in range(ADVERSARIAL_PAIRS):

        a, b = random_wrong_pair(
            p,
            q
        )

        candidate = all_prefix_signatures(
            a,
            b,
            prefix
        )

        for name in SIGNATURES:

            if len(found[name]) > 0:
                continue

            if candidate[name] == actual[name]:

                found[name].append(
                    (a, b)
                )

        if all(
            found[name]
            for name in SIGNATURES
        ):
            break

    return found


# ============================================================
# FIRST ZERO PREFIX
# ============================================================

def first_zero_prefix(
    collision_counts
):

    for prefix in PREFIXES:

        if collision_counts[prefix] == 0:
            return prefix

    return None


# ============================================================
# FORMAT
# ============================================================

def fmt_prob(
    hits,
    trials
):

    if hits == 0:
        return f"<1/{trials:,}"

    return f"{hits / trials:.8g}"


# ============================================================
# MAIN HEADER
# ============================================================

print("=" * 78)
print("KAPPA EXPERIMENT 22")
print("COLLISION THRESHOLD / ADVERSARIAL COLLISION SEARCH")
print("NO CSV OUTPUT")
print("=" * 78)

print(
    f"random seed       = {SEED}"
)

print(
    f"R values          = {R_VALUES}"
)

print(
    f"targets           = {TARGETS}"
)

print(
    f"random pairs      = {WRONG_PAIRS:,}"
)

print(
    f"adversarial pairs = {ADVERSARIAL_PAIRS:,}"
)

print(
    f"prefixes          = {PREFIXES}"
)

print()


# ============================================================
# 1. MODULUS INVENTORY
# ============================================================

print("=" * 78)
print("1. MODULUS INVENTORY")
print("=" * 78)

for r in R_VALUES:

    m = r * r + 3

    orbit = INVENTORY[m]["orbit"]

    print(
        f"r={r:3d} "
        f"m={m:6d} "
        f"units={orbit:7d}"
    )

print()


# ============================================================
# 2. THEORETICAL PREFIX BASELINES
# ============================================================

print("=" * 78)
print("2. THEORETICAL FPAIR PREFIX BASELINES")
print("=" * 78)

print(
    f"{'prefix':>8s}"
    f"{'last m':>12s}"
    f"{'log10 generic freq':>24s}"
)

for prefix in PREFIXES:

    print(
        f"{prefix:8d}"
        f"{MODULI[prefix - 1]:12d}"
        f"{theoretical_log10(prefix):24.6f}"
    )

print()


# ============================================================
# STORAGE FOR FINAL SUMMARY
# ============================================================

all_random_hits = {
    prefix: {
        name: []
        for name in SIGNATURES
    }
    for prefix in PREFIXES
}

all_adversarial_found = {
    prefix: {
        name: 0
        for name in SIGNATURES
    }
    for prefix in PREFIXES
}

zero_counts = {
    prefix: {
        name: 0
        for name in SIGNATURES
    }
    for prefix in PREFIXES
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


    # --------------------------------------------------------
    # RANDOM COLLISION DEPTH
    # --------------------------------------------------------

    print(
        "RANDOM COLLISION DEPTH"
    )

    print("-" * 78)

    print(
        f"{'prefix':>7s} "
        f"{'m':>7s} "
        + " ".join(
            f"{name:>13s}"
            for name in SIGNATURES
        )
    )

    target_random = {}

    for prefix in PREFIXES:

        hits = random_collision_test(
            p,
            q,
            prefix
        )

        target_random[prefix] = hits

        row = [
            f"{prefix:7d}",
            f"{MODULI[prefix - 1]:7d}",
        ]

        for name in SIGNATURES:

            h = hits[name]

            all_random_hits[
                prefix
            ][name].append(h)

            if h == 0:
                zero_counts[
                    prefix
                ][name] += 1

            row.append(
                f"{fmt_prob(h, WRONG_PAIRS):>13s}"
            )

        print(" ".join(row))

    print()


    # --------------------------------------------------------
    # FIRST ZERO PREFIX
    # --------------------------------------------------------

    print(
        "FIRST OBSERVED ZERO-COLLISION PREFIX"
    )

    print("-" * 78)

    for name in SIGNATURES:

        counts = [
            target_random[prefix][name]
            for prefix in PREFIXES
        ]

        mapping = {
            prefix: counts[i]
            for i, prefix in enumerate(PREFIXES)
        }

        first = first_zero_prefix(
            mapping
        )

        if first is None:
            text = "none"

        else:
            text = str(first)

        print(
            f"{name:18s}{text:>8s}"
        )

    print()


    # --------------------------------------------------------
    # ADVERSARIAL SEARCH
    # --------------------------------------------------------

    print(
        "ADVERSARIAL COLLISION SEARCH"
    )

    print("-" * 78)

    for prefix in [5, 6, 7]:

        print(
            f"PREFIX {prefix} "
            f"(m={MODULI[prefix - 1]})"
        )

        found = adversarial_search(
            p,
            q,
            prefix
        )

        for name in SIGNATURES:

            pair_list = found[name]

            if pair_list:

                a, b = pair_list[0]

                all_adversarial_found[
                    prefix
                ][name] += 1

                print(
                    f"  {name:18s}"
                    f"COLLISION "
                    f"({a}, {b})"
                )

            else:

                print(
                    f"  {name:18s}"
                    f"NONE FOUND"
                )

        print()

    # --------------------------------------------------------
    # TARGET CONCLUSION
    # --------------------------------------------------------

    print(
        "TARGET CONCLUSION"
    )

    print("-" * 78)

    for name in SIGNATURES:

        first = None

        for prefix in PREFIXES:

            if target_random[prefix][name] == 0:

                first = prefix
                break

        if first is None:

            print(
                f"{name:18s}"
                f"no zero prefix"
            )

        else:

            print(
                f"{name:18s}"
                f"first zero = {first}"
            )

    print()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print()


# ============================================================
# MEDIAN RANDOM COLLISION COUNTS
# ============================================================

print(
    "MEDIAN RANDOM COLLISION COUNTS"
)

print("-" * 78)

print(
    f"{'prefix':>7s} "
    + " ".join(
        f"{name:>13s}"
        for name in SIGNATURES
    )
)

for prefix in PREFIXES:

    row = [
        f"{prefix:7d}"
    ]

    for name in SIGNATURES:

        values = sorted(
            all_random_hits[
                prefix
            ][name]
        )

        middle = len(values) // 2

        if len(values) % 2:

            median = values[middle]

        else:

            median = (
                values[middle - 1]
                + values[middle]
            ) / 2

        row.append(
            f"{median:13.1f}"
        )

    print(" ".join(row))

print()


# ============================================================
# ZERO-COUNT TARGETS
# ============================================================

print(
    "TARGETS WITH ZERO RANDOM COLLISIONS"
)

print("-" * 78)

print(
    f"{'signature':18s}"
    + "".join(
        f"{prefix:>10d}"
        for prefix in PREFIXES
    )
)

for name in SIGNATURES:

    row = [
        f"{name:18s}"
    ]

    for prefix in PREFIXES:

        row.append(
            f"{zero_counts[prefix][name]:10d}"
        )

    print("".join(row))

print()


# ============================================================
# ADVERSARIAL SEARCH SUMMARY
# ============================================================

print(
    "ADVERSARIAL COLLISION SUMMARY"
)

print("-" * 78)

print(
    "Number of targets for which at least one"
)

print(
    "concrete wrong-pair collision was found."
)

print()

print(
    f"{'signature':18s}"
    f"{'prefix 5':>12s}"
    f"{'prefix 6':>12s}"
    f"{'prefix 7':>12s}"
)

for name in SIGNATURES:

    print(
        f"{name:18s}"
        f"{all_adversarial_found[5][name]:12d}"
        f"{all_adversarial_found[6][name]:12d}"
        f"{all_adversarial_found[7][name]:12d}"
    )

print()


# ============================================================
# EMPIRICAL VS THEORY
# ============================================================

print(
    "FPAIR: EMPIRICAL VS THEORETICAL"
)

print("-" * 78)

print(
    f"{'prefix':>8s}"
    f"{'theory log10':>18s}"
    f"{'median empirical':>22s}"
)

for prefix in PREFIXES:

    values = sorted(
        all_random_hits[
            prefix
        ]["Fpair"]
    )

    middle = len(values) // 2

    if len(values) % 2:

        median_hits = values[middle]

    else:

        median_hits = (
            values[middle - 1]
            + values[middle]
        ) / 2

    if median_hits == 0:

        empirical = (
            f"<1/{WRONG_PAIRS:,}"
        )

    else:

        empirical = (
            f"{median_hits / WRONG_PAIRS:.8g}"
        )

    print(
        f"{prefix:8d}"
        f"{theoretical_log10(prefix):18.6f}"
        f"{empirical:>22s}"
    )

print()


# ============================================================
# INTERPRETATION
# ============================================================

print("=" * 78)
print("INTERPRETATION")
print("=" * 78)

print("""
Experiment 22 is designed around the transition observed in
Experiment 21.

Experiment 21 showed:

    prefix 3:
        substantial collisions for every signature

    prefix 5:
        Fpair/Fsorted consistently reached zero observed collisions

    prefix 7:
        essentially every signature reached zero observed collisions

The key question is now whether prefix 5 represents a genuine
collision-depth threshold or merely an artifact of sampling
20,000 random wrong pairs.

This experiment therefore performs TWO tests.

1. RANDOM COLLISION DEPTH

The same random-pair methodology is applied at prefixes:

    3, 4, 5, 6, 7

This resolves the transition more finely.

2. ADVERSARIAL COLLISION SEARCH

For prefixes 5, 6 and 7, substantially more wrong pairs are
searched and the first concrete collision is recorded.

This is important because:

    zero random collisions

does not mean:

    no collision exists.

A concrete collision found by the adversarial search would establish
that the corresponding signature is not injective at that prefix.

Conversely, failure to find a collision is still empirical evidence,
not a proof of uniqueness.

The especially important comparison is:

    Fpair vs Fsorted
    Fpair vs Fsum
    Fpair vs Fprod
    Fpair vs Fdiff
    Fpair vs cube_sum
    Fpair vs order_pair
    Fpair vs order_sorted

If Fpair repeatedly reaches zero at prefix 5 while compressed
signatures require prefix 6 or 7, that strengthens the conclusion
that retaining the complete residue pair carries materially more
discriminating information.

If adversarial search finds collisions for Fpair at prefix 5,
then the apparent Experiment-21 threshold was only a sampling
threshold and the next experiment should move toward exhaustive
local-state enumeration rather than larger random samples.

The most important output is therefore not the absolute number
of moduli, but the first prefix at which collisions become
empirically absent AND remain absent under the stronger search.
""")

print()
print("=" * 78)
print("DONE")
print("=" * 78)

