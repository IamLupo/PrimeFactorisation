#!/usr/bin/env python3

import math
import random
from collections import defaultdict

# ============================================================
# KAPPA EXPERIMENT 25
# CRT-CONSTRUCTIVE FPAIR COLLISION SEARCH
#
# NO CSV OUTPUT
#
# Goal:
#   Avoid O(P^2) pair enumeration.
#
#   For each prefix of moduli, construct the exact residue
#   requirements for Fpair and then look for distinct primes
#   in the prime pool occupying those residue classes.
#
#   This tests whether concrete Fpair collisions survive when
#   the residue constraints are attacked directly.
# ============================================================

SEED = 20260814
random.seed(SEED)

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

TARGETS = 12

PRIME_POOL_SIZE = 5_000

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]


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


def generate_prime_pool():

    pool = set()

    while len(pool) < PRIME_POOL_SIZE:

        n = random.randrange(
            PRIME_LOW,
            PRIME_HIGH
        )

        if n % 2 == 0:
            n += 1

        if is_prime(n):
            pool.add(n)

    return sorted(pool)


def factor_pair(prime_pool):

    while True:

        p = random.choice(prime_pool)
        q = random.choice(prime_pool)

        if p != q:
            return p, q


# ============================================================
# MODULI
# ============================================================

MODULI = [
    r * r + 3
    for r in R_VALUES
]


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


# ============================================================
# INVENTORY
# ============================================================

INVENTORY = {}

for m in MODULI:

    units = units_mod(m)

    orders = {
        u: multiplicative_order(u, m)
        for u in units
    }

    INVENTORY[m] = {
        "units": units,
        "orbit": len(units),
        "orders": orders,
    }


# ============================================================
# CRT
# ============================================================

def gcd_extended(a, b):

    if b == 0:
        return a, 1, 0

    g, x1, y1 = gcd_extended(
        b,
        a % b
    )

    return (
        g,
        y1,
        x1 - (a // b) * y1
    )


def crt_pair(residues, moduli):

    """
    General CRT for pairwise coprime moduli.

    Returns:

        x modulo product(moduli)

    or None if incompatible.
    """

    x = 0
    M = 1

    for a, m in zip(residues, moduli):

        g = math.gcd(M, m)

        if (a - x) % g != 0:
            return None

        M1 = M // g
        m1 = m // g

        _, inv, _ = gcd_extended(
            M1,
            m1
        )

        k = (
            ((a - x) // g)
            * inv
        ) % m1

        x += M * k
        M *= m1

        x %= M

    return x


# ============================================================
# PREFIX MODULI
# ============================================================

def prefix_moduli(prefix):

    return MODULI[:prefix]


def prefix_modulus_product(prefix):

    value = 1

    for m in prefix_moduli(prefix):
        value *= m

    return value


# ============================================================
# TARGET RESIDUES
# ============================================================

def target_residues(p, q, prefix):

    mods = prefix_moduli(prefix)

    rp = [
        p % m
        for m in mods
    ]

    rq = [
        q % m
        for m in mods
    ]

    return rp, rq


# ============================================================
# DIRECT PRIME CLASS INDEX
# ============================================================

def build_residue_index(prime_pool, modulus):

    index = defaultdict(list)

    for p in prime_pool:

        index[p % modulus].append(p)

    return index


# ============================================================
# CONSTRUCTIVE FPAIR SEARCH
# ============================================================

def constructive_fpair_search(
    target_p,
    target_q,
    prefix,
    residue_indexes
):

    mods = prefix_moduli(prefix)

    target_rp, target_rq = target_residues(
        target_p,
        target_q,
        prefix
    )

    # --------------------------------------------------------
    # Search orientation 1:
    #
    # a has target_p residues
    # b has target_q residues
    # --------------------------------------------------------

    key_a = tuple(target_rp)
    key_b = tuple(target_rq)

    candidates_a = residue_indexes[prefix].get(
        key_a,
        []
    )

    candidates_b = residue_indexes[prefix].get(
        key_b,
        []
    )

    for a in candidates_a:

        for b in candidates_b:

            if a == b:
                continue

            if {a, b} == {target_p, target_q}:
                continue

            return a, b

    # --------------------------------------------------------
    # Search reversed orientation.
    # --------------------------------------------------------

    candidates_a = residue_indexes[prefix].get(
        key_b,
        []
    )

    candidates_b = residue_indexes[prefix].get(
        key_a,
        []
    )

    for a in candidates_a:

        for b in candidates_b:

            if a == b:
                continue

            if {a, b} == {target_p, target_q}:
                continue

            return a, b

    return None


# ============================================================
# COMPLETE FPAIR SIGNATURE
# ============================================================

def fpair_signature(p, q, prefix):

    return tuple(
        (
            p % m,
            q % m
        )
        for m in prefix_moduli(prefix)
    )


# ============================================================
# FSORTED SIGNATURE
# ============================================================

def fsorted_signature(p, q, prefix):

    return tuple(
        tuple(sorted((p % m, q % m)))
        for m in prefix_moduli(prefix)
    )


# ============================================================
# VERIFY COLLISION
# ============================================================

def verify_collision(
    target_p,
    target_q,
    candidate_p,
    candidate_q,
    prefix
):

    return (
        fpair_signature(
            target_p,
            target_q,
            prefix
        )
        ==
        fpair_signature(
            candidate_p,
            candidate_q,
            prefix
        )
    )


# ============================================================
# THEORETICAL BASELINE
# ============================================================

def theoretical_log10(prefix):

    return -sum(
        math.log10(
            INVENTORY[m]["orbit"]
        )
        for m in prefix_moduli(prefix)
    )


# ============================================================
# MAIN
# ============================================================

print("=" * 78)
print("KAPPA EXPERIMENT 25")
print("CRT-CONSTRUCTIVE FPAIR COLLISION SEARCH")
print("NO CSV OUTPUT")
print("=" * 78)

print(
    f"random seed      = {SEED}"
)

print(
    f"R values         = {R_VALUES}"
)

print(
    f"targets          = {TARGETS}"
)

print(
    f"prime pool       = {PRIME_POOL_SIZE:,}"
)

print(
    f"prime interval   = "
    f"[{PRIME_LOW:,}, {PRIME_HIGH:,}]"
)

print(
    f"prefixes         = {PREFIXES}"
)

print()


# ============================================================
# MODULUS INVENTORY
# ============================================================

print("=" * 78)
print("1. MODULUS INVENTORY")
print("=" * 78)

for r, m in zip(R_VALUES, MODULI):

    print(
        f"r={r:3d} "
        f"m={m:6d} "
        f"units={INVENTORY[m]['orbit']:7d}"
    )

print()


# ============================================================
# BASELINES
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
# PRIME POOL
# ============================================================

print("=" * 78)
print("3. PRIME POOL GENERATION")
print("=" * 78)

print("Generating unique prime pool...")

prime_pool = generate_prime_pool()

print(
    f"generated primes = {len(prime_pool):,}"
)

print()


# ============================================================
# RESIDUE INDEXES
# ============================================================

print("=" * 78)
print("4. BUILDING PREFIX RESIDUE INDEXES")
print("=" * 78)

"""
Instead of indexing every prime-pair combination, we index
individual primes by their COMPLETE residue vector.

For example, prefix=5 stores:

    (
        p mod 7,
        p mod 12,
        p mod 28,
        p mod 52,
        p mod 124
    )

The target pair therefore requires another prime with exactly
the target's residue vector.

This reduces the search from approximately:

    P^2

to approximately:

    P

per prefix.
"""

residue_indexes = {}

for prefix in PREFIXES:

    print(
        f"building prefix {prefix} "
        f"(last m={MODULI[prefix - 1]})..."
    )

    index = defaultdict(list)

    mods = prefix_moduli(prefix)

    for p in prime_pool:

        key = tuple(
            p % m
            for m in mods
        )

        index[key].append(p)

    residue_indexes[prefix] = index

    occupied = len(index)

    print(
        f"  unique residue classes = "
        f"{occupied:,}"
    )

print()

print("residue indexes complete")
print()


# ============================================================
# TARGET LOOP
# ============================================================

results = {
    prefix: 0
    for prefix in PREFIXES
}

print("=" * 78)
print("5. CONSTRUCTIVE COLLISION SEARCH")
print("=" * 78)

for target in range(1, TARGETS + 1):

    p, q = factor_pair(prime_pool)

    n = p * q
    s = p + q
    delta = (p - q) ** 2

    print()
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

    print(
        "CONSTRUCTIVE FPAIR COLLISION SEARCH"
    )

    print("-" * 78)

    print(
        f"{'prefix':>8s}"
        f"{'last m':>10s}"
        f"{'candidate':>34s}"
        f"{'verified':>12s}"
    )

    for prefix in PREFIXES:

        candidate = constructive_fpair_search(
            p,
            q,
            prefix,
            residue_indexes
        )

        if candidate is None:

            print(
                f"{prefix:8d}"
                f"{MODULI[prefix - 1]:10d}"
                f"{'NONE FOUND':>34s}"
                f"{'--':>12s}"
            )

            continue

        a, b = candidate

        verified = verify_collision(
            p,
            q,
            a,
            b,
            prefix
        )

        if verified:
            results[prefix] += 1

        print(
            f"{prefix:8d}"
            f"{MODULI[prefix - 1]:10d}"
            f"{str(candidate):>34s}"
            f"{str(verified):>12s}"
        )

    print()


# ============================================================
# PREFIX SURVIVAL
# ============================================================

print("=" * 78)
print("6. COLLISION SURVIVAL SUMMARY")
print("=" * 78)

print()

print(
    f"{'prefix':>8s}"
    f"{'last m':>10s}"
    f"{'targets with concrete collision':>36s}"
)

for prefix in PREFIXES:

    print(
        f"{prefix:8d}"
        f"{MODULI[prefix - 1]:10d}"
        f"{results[prefix]:36d}/{TARGETS}"
    )

print()


# ============================================================
# INTERPRETATION
# ============================================================

print("=" * 78)
print("INTERPRETATION")
print("=" * 78)

print("""
Experiment 25 changes the search strategy.

Experiments 23 and 24 attempted to discover collisions by
examining large numbers of candidate prime pairs.

This experiment does NOT enumerate prime pairs.

Instead, every prime is indexed by its complete residue vector
over the selected prefix.

For example, prefix 5 represents:

    p mod 7
    p mod 12
    p mod 28
    p mod 52
    p mod 124

A concrete Fpair collision therefore exists in the pool exactly
when another prime occupies the same complete residue vector as
one of the target factors.

The computational cost is approximately O(P * number_of_prefixes),
rather than O(P^2).

IMPORTANT:

A collision found here is much stronger evidence than a random
collision estimate because it provides an explicit wrong factor
pair.

Conversely, NONE FOUND does not prove uniqueness. It only means
that the selected finite prime pool contains no second prime in
the required residue class.

The most important transition is:

    collision at prefix 3
        ->
    collision at prefix 4
        ->
    collision at prefix 5
        ->
    collision disappears

If the disappearance consistently occurs around prefix 5 or 6,
that gives a much sharper empirical estimate of the constructive
collision depth than random sampling alone.

Also note that this experiment specifically tests Fpair.

That is intentional: Fpair is the reference signature against
which the compressed signatures can later be compared.
""")

print()
print("=" * 78)
print("DONE")
print("=" * 78)

