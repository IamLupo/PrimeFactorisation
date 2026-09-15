#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# ============================================================
# KAPPA EXPERIMENT 20
# COMPLETE-SIGNATURE COLLISION / DISCRIMINATION TEST
#
# NO CSV OUTPUT
# ============================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
            31, 37, 41, 43, 47]

TARGETS = 12

# Number of wrong candidate pairs tested per target.
# Increase for a stronger collision experiment.
WRONG_PAIRS = 100_000

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


# ============================================================
# MODULUS / UNIT ORBIT INVENTORY
# ============================================================

MODULI = [r * r + 3 for r in R_VALUES]


def units_mod(m):
    return [x for x in range(1, m) if math.gcd(x, m) == 1]


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
    """
    Return the complete local information used by the experiment.

    Pair is ordered.
    Sorted removes ordering.
    Sum/product/difference are scalar signatures.
    Order pair records multiplicative orders.
    """

    a = p % m
    b = q % m

    oa = INVENTORY[m]["orders"].get(a)
    ob = INVENTORY[m]["orders"].get(b)

    return {
        "pair": (a, b),
        "sorted": tuple(sorted((a, b))),
        "sum": (a + b) % m,
        "prod": (a * b) % m,
        "diff": (a - b) % m,
        "order_pair": (oa, ob),
        "order_sorted": tuple(sorted((oa, ob))),
        "cube_sum": (a ** 3 + b ** 3) % m,
    }


# ============================================================
# COMPLETE SIGNATURE
# ============================================================

def complete_signature(p, q):
    result = {
        "Fpair": [],
        "Fsorted": [],
        "Fsum": [],
        "Fprod": [],
        "Fdiff": [],
        "cube_sum": [],
        "order_pair": [],
        "order_sorted": [],
    }

    for m in MODULI:
        s = local_signature(p, q, m)

        result["Fpair"].append(s["pair"])
        result["Fsorted"].append(s["sorted"])
        result["Fsum"].append(s["sum"])
        result["Fprod"].append(s["prod"])
        result["Fdiff"].append(s["diff"])
        result["cube_sum"].append(s["cube_sum"])
        result["order_pair"].append(s["order_pair"])
        result["order_sorted"].append(s["order_sorted"])

    return {k: tuple(v) for k, v in result.items()}


# ============================================================
# PARTIAL SIGNATURE
# ============================================================

def partial_signature(p, q, signature_name, modulus_indices):
    out = []

    for i in modulus_indices:
        m = MODULI[i]
        s = local_signature(p, q, m)

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
# JOINT EMPIRICAL COLLISION TEST
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


def random_wrong_pair(p, q):
    while True:
        a = random_prime()
        b = random_prime()

        if a == b:
            continue

        if {a, b} != {p, q}:
            return a, b


def collision_test(p, q, wrong_pairs):
    actual = complete_signature(p, q)

    counts = Counter()

    exact_collisions = defaultdict(int)

    for _ in range(wrong_pairs):
        a, b = random_wrong_pair(p, q)
        candidate = complete_signature(a, b)

        for name in SIGNATURES:
            if candidate[name] == actual[name]:
                counts[name] += 1

    return counts


# ============================================================
# MODULUS-BY-MODULUS MATCH TEST
# ============================================================

def local_match_rates(p, q, wrong_pairs):
    """
    For each signature and modulus, estimate:

        P(candidate local signature == actual local signature)

    using random wrong factor pairs.
    """

    hits = {
        name: [0] * len(MODULI)
        for name in SIGNATURES
    }

    for _ in range(wrong_pairs):
        a, b = random_wrong_pair(p, q)

        for i, m in enumerate(MODULI):
            actual = local_signature(p, q, m)
            cand = local_signature(a, b, m)

            for name in SIGNATURES:
                if name == "Fpair":
                    x = actual["pair"]
                    y = cand["pair"]
                elif name == "Fsorted":
                    x = actual["sorted"]
                    y = cand["sorted"]
                elif name == "Fsum":
                    x = actual["sum"]
                    y = cand["sum"]
                elif name == "Fprod":
                    x = actual["prod"]
                    y = cand["prod"]
                elif name == "Fdiff":
                    x = actual["diff"]
                    y = cand["diff"]
                elif name == "cube_sum":
                    x = actual["cube_sum"]
                    y = cand["cube_sum"]
                elif name == "order_pair":
                    x = actual["order_pair"]
                    y = cand["order_pair"]
                elif name == "order_sorted":
                    x = actual["order_sorted"]
                    y = cand["order_sorted"]

                if x == y:
                    hits[name][i] += 1

    return hits


# ============================================================
# LEAVE-ONE-MODULUS-OUT ANALYSIS
# ============================================================

def leave_one_out(p, q):
    actual = complete_signature(p, q)

    result = {}

    for name in SIGNATURES:
        result[name] = []

        for omitted in range(len(MODULI)):
            indices = [
                i for i in range(len(MODULI))
                if i != omitted
            ]

            sig = tuple(actual[name][i] for i in indices)
            result[name].append(sig)

    return result


# ============================================================
# FORMATTING
# ============================================================

def fmt_prob(x, n):
    if n == 0:
        return "0"

    p = x / n

    if p == 0:
        return "<1/" + f"{n:,}"

    return f"{p:.8g}"


def log10_product_orbit():
    return -sum(
        math.log10(INVENTORY[m]["orbit"])
        for m in MODULI
    )


# ============================================================
# MAIN
# ============================================================

print("=" * 78)
print("KAPPA EXPERIMENT 20")
print("COMPLETE-SIGNATURE COLLISION / DISCRIMINATION TEST")
print("NO CSV OUTPUT")
print("=" * 78)

print(f"random seed = {SEED}")
print(f"R values    = {R_VALUES}")
print(f"targets     = {TARGETS}")
print(f"wrong pairs = {WRONG_PAIRS}")
print()

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

baseline = log10_product_orbit()

print()
print("GENERIC FPAIR BASELINE")
print("-" * 78)
print(f"log10(product 1/orbit) = {baseline:.6f}")
print()


# ============================================================
# TARGET LOOP
# ============================================================

all_collision_counts = {
    name: []
    for name in SIGNATURES
}

all_local_hits = {
    name: [0] * len(MODULI)
    for name in SIGNATURES
}

all_local_trials = WRONG_PAIRS * TARGETS

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

    actual = complete_signature(p, q)

    # --------------------------------------------------------
    # Exact complete-signature collision test
    # --------------------------------------------------------

    counts = collision_test(
        p,
        q,
        WRONG_PAIRS
    )

    print("COMPLETE SIGNATURE COLLISION TEST")
    print("-" * 78)
    print(
        f"{'signature':18s}"
        f"{'wrong matches':>18s}"
        f"{'empirical freq':>20s}"
    )

    for name in SIGNATURES:
        c = counts[name]

        print(
            f"{name:18s}"
            f"{c:18d}"
            f"{fmt_prob(c, WRONG_PAIRS):>20s}"
        )

        all_collision_counts[name].append(c)

    print()

    # --------------------------------------------------------
    # Local match test
    # --------------------------------------------------------

    hits = local_match_rates(
        p,
        q,
        max(1000, WRONG_PAIRS // 10)
    )

    local_trials = max(1000, WRONG_PAIRS // 10)

    print("LOCAL MATCH RATES")
    print("-" * 78)

    header = (
        f"{'m':>6s} "
        + " ".join(f"{name:>13s}" for name in SIGNATURES)
    )

    print(header)

    for i, m in enumerate(MODULI):
        row = [f"{m:6d}"]

        for name in SIGNATURES:
            h = hits[name][i]

            all_local_hits[name][i] += h

            row.append(
                f"{fmt_prob(h, local_trials):>13s}"
            )

        print(" ".join(row))

    print()

    # --------------------------------------------------------
    # Leave-one-out structural test
    # --------------------------------------------------------

    print("LEAVE-ONE-MODULUS-OUT FPAIR COLLISION STRUCTURE")
    print("-" * 78)

    for omitted, m in enumerate(MODULI):

        indices = [
            i for i in range(len(MODULI))
            if i != omitted
        ]

        sig = tuple(actual["Fpair"][i] for i in indices)

        # The theoretical generic frequency is obtained by
        # removing the omitted orbit factor.
        remaining_log10 = -sum(
            math.log10(INVENTORY[MODULI[i]]["orbit"])
            for i in indices
        )

        print(
            f"omit m={m:6d} "
            f"log10 generic Fpair = {remaining_log10: .6f}"
        )

    print()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print()
print("COMPLETE-SIGNATURE COLLISION RATES")
print("-" * 78)

print(
    f"{'signature':18s}"
    f"{'median hits':>16s}"
    f"{'median rate':>18s}"
)

for name in SIGNATURES:

    values = all_collision_counts[name]

    ordered = sorted(values)

    if len(ordered) % 2:
        median_hits = ordered[len(ordered) // 2]
    else:
        median_hits = (
            ordered[len(ordered) // 2 - 1]
            + ordered[len(ordered) // 2]
        ) / 2

    median_rate = median_hits / WRONG_PAIRS

    print(
        f"{name:18s}"
        f"{median_hits:16.1f}"
        f"{median_rate:18.8g}"
    )


print()
print("WIN COUNTS: FEWEST WRONG-PAIR COLLISIONS")
print("-" * 78)

wins = Counter()

for t in range(TARGETS):

    best_name = min(
        SIGNATURES,
        key=lambda name:
        all_collision_counts[name][t]
    )

    wins[best_name] += 1

for name in SIGNATURES:
    print(
        f"{name:18s}"
        f"{wins[name]:3d}/{TARGETS}"
    )


print()
print("AGGREGATE LOCAL MATCH RATES")
print("-" * 78)

print(
    f"{'m':>6s} "
    + " ".join(f"{name:>13s}" for name in SIGNATURES)
)

for i, m in enumerate(MODULI):

    row = [f"{m:6d}"]

    for name in SIGNATURES:
        rate = (
            all_local_hits[name][i]
            / all_local_trials
        )

        row.append(
            f"{rate:13.8g}"
        )

    print(" ".join(row))


print()
print("=" * 78)
print("INTERPRETATION")
print("=" * 78)

print("""
The critical statistic is COMPLETE-SIGNATURE COLLISION RATE.

For each target, random incorrect factor pairs are generated.
A collision means that the incorrect pair produces exactly the
same complete signature as the actual factor pair.

Fpair should have extremely few or zero collisions if the complete
local residue pair really identifies the factor pair.

Fsorted tests whether ordering is essential.

Fsum/Fprod/Fdiff/cube_sum test whether compressed scalar invariants
retain enough information to reproduce the same signature.

order_pair/order_sorted test multiplicative-order discrimination.

The leave-one-modulus-out section shows how much of the Fpair
discrimination survives after removing each modulus.

IMPORTANT:
The theoretical Fpair baseline is determined largely by the chosen
unit-orbit sizes. Therefore empirical collision testing against
wrong candidate pairs is the decisive next test.
""")

print()
print("=" * 78)
print("DONE")
print("=" * 78)

