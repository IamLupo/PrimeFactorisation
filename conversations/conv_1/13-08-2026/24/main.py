#!/usr/bin/env python3

import math
import random
from collections import defaultdict

# ============================================================================
# KAPPA EXPERIMENT 24
# HASH-INDEXED CONSTRUCTIVE COLLISION SEARCH
#
# NO CSV OUTPUT
#
# Goal:
#   Find concrete wrong prime pairs having the same partial signature
#   as a target pair, without repeatedly scanning the entire pair space.
#
# Main optimization:
#   Instead of:
#
#       for target:
#           for candidate pair:
#               calculate signatures repeatedly
#
#   we build hash indexes incrementally.
#
# ============================================================================


SEED = 20260814
random.seed(SEED)

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

TARGETS = 12

PRIME_POOL_SIZE = 5_000

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

# Maximum number of distinct prime pairs to inspect.
#
# This is a safety limit. The search normally stops much earlier
# because collisions are found through hash tables.
MAX_PAIR_SCANS = 8_000_000


# ============================================================================
# SIGNATURES
# ============================================================================

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


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

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
    primes = set()

    print("Generating unique prime pool...")

    while len(primes) < PRIME_POOL_SIZE:
        n = random.randrange(PRIME_LOW, PRIME_HIGH)

        if n % 2 == 0:
            n += 1

        if is_prime(n):
            primes.add(n)

    primes = sorted(primes)

    print(f"generated primes = {len(primes)}")
    print()

    return primes


# ============================================================================
# MODULUS INVENTORY
# ============================================================================

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


# ============================================================================
# PRIME CACHE
# ============================================================================

class PrimeCache:

    def __init__(self, primes):

        self.primes = primes

        # residues[m][index]
        self.residues = {}

        # orders[m][index]
        self.orders = {}

        self.build()

    def build(self):

        print("Precomputing residues and multiplicative orders...")
        
        for m in MODULI:

            residue_list = []
            order_list = []

            order_table = INVENTORY[m]["orders"]

            for p in self.primes:

                a = p % m

                residue_list.append(a)
                order_list.append(
                    order_table.get(a)
                )

            self.residues[m] = residue_list
            self.orders[m] = order_list

        print("cache complete")
        print()


# ============================================================================
# LOCAL SIGNATURE FROM CACHE
# ============================================================================

def local_signature_indices(i, j, m, cache):

    a = cache.residues[m][i]
    b = cache.residues[m][j]

    oa = cache.orders[m][i]
    ob = cache.orders[m][j]

    return {
        "Fpair": (a, b),

        "Fsorted": tuple(
            sorted((a, b))
        ),

        "Fsum": (a + b) % m,

        "Fprod": (a * b) % m,

        "Fdiff": (a - b) % m,

        "cube_sum": (
            a ** 3 + b ** 3
        ) % m,

        "order_pair": (oa, ob),

        "order_sorted": tuple(
            sorted((oa, ob))
        ),
    }


# ============================================================================
# PREFIX SIGNATURE
# ============================================================================

def prefix_signature(i, j, prefix, signature_name, cache):

    values = []

    for k in range(prefix):

        m = MODULI[k]

        s = local_signature_indices(
            i,
            j,
            m,
            cache
        )

        values.append(
            s[signature_name]
        )

    return tuple(values)


# ============================================================================
# FAST INCREMENTAL SIGNATURE
# ============================================================================

def extend_signature(
    current,
    i,
    j,
    prefix,
    signature_name,
    cache
):

    m = MODULI[prefix - 1]

    s = local_signature_indices(
        i,
        j,
        m,
        cache
    )

    return current + (s[signature_name],)


# ============================================================================
# TARGET GENERATION
# ============================================================================

def choose_targets(primes):

    targets = []

    while len(targets) < TARGETS:

        i, j = random.sample(
            range(len(primes)),
            2
        )

        p = primes[i]
        q = primes[j]

        targets.append(
            (i, j)
        )

    return targets


# ============================================================================
# FORMATTING
# ============================================================================

def fmt_pair(pair, primes):

    if pair is None:
        return "NONE FOUND"

    i, j = pair

    return (
        f"({primes[i]}, {primes[j]})"
    )


def theoretical_log10(prefix):

    return -sum(
        math.log10(
            INVENTORY[m]["orbit"]
        )
        for m in MODULI[:prefix]
    )


# ============================================================================
# BUILD TARGET SIGNATURES
# ============================================================================

def build_target_signatures(
    target_i,
    target_j,
    cache
):

    result = {
        name: {}
        for name in SIGNATURES
    }

    for name in SIGNATURES:

        current = ()

        for prefix in range(1, max(PREFIXES) + 1):

            current = extend_signature(
                current,
                target_i,
                target_j,
                prefix,
                name,
                cache
            )

            result[name][prefix] = current

    return result


# ============================================================================
# HASH-INDEXED COLLISION SEARCH
# ============================================================================

def search_collisions(
    target_i,
    target_j,
    primes,
    cache
):

    max_prefix = max(PREFIXES)

    target_signatures = build_target_signatures(
        target_i,
        target_j,
        cache
    )

    # ------------------------------------------------------------------------
    # For every signature and prefix we need only the FIRST collision.
    #
    # index[name][prefix][signature] = (i, j)
    #
    # Once a signature has been seen, we can immediately compare it against
    # the target signature.
    # ------------------------------------------------------------------------

    indexes = {
        name: {
            prefix: {}
            for prefix in PREFIXES
        }
        for name in SIGNATURES
    }

    found = {
        name: {
            prefix: None
            for prefix in PREFIXES
        }
        for name in SIGNATURES
    }

    # We process pairs in a deterministic order.
    #
    # Only i < j is used because the pool contains unordered factor pairs.
    #
    # Fpair / order_pair remain ordered signatures internally according to
    # the i,j orientation.
    # ------------------------------------------------------------------------

    n = len(primes)

    pair_scans = 0

    for i in range(n):

        for j in range(i + 1, n):

            pair_scans += 1

            if pair_scans > MAX_PAIR_SCANS:
                return found, pair_scans

            # Never count the actual target pair as a collision.
            if (
                (i == target_i and j == target_j)
                or
                (i == target_j and j == target_i)
            ):
                continue

            # ----------------------------------------------------------------
            # Build each signature progressively.
            # ----------------------------------------------------------------

            current = {
                name: ()
                for name in SIGNATURES
            }

            for prefix in range(1, max_prefix + 1):

                m = MODULI[prefix - 1]

                s = local_signature_indices(
                    i,
                    j,
                    m,
                    cache
                )

                for name in SIGNATURES:

                    current[name] += (
                        s[name],
                    )

                # Only inspect requested prefixes.
                if prefix not in PREFIXES:
                    continue

                for name in SIGNATURES:

                    if found[name][prefix] is not None:
                        continue

                    target_sig = (
                        target_signatures[name][prefix]
                    )

                    candidate_sig = current[name]

                    if candidate_sig == target_sig:

                        found[name][prefix] = (
                            i,
                            j
                        )

                # If every requested signature has a collision at this
                # prefix, there is nothing more to discover for this target
                # at this prefix.
                #
                # We still continue because longer prefixes are independent.
                # ----------------------------------------------------------------

    return found, pair_scans


# ============================================================================
# MORE EFFICIENT PREFIX SEARCH
# ============================================================================
#
# The function above is deliberately straightforward.
#
# This second implementation is used for the actual experiment.
#
# It builds a hash table at each prefix and stores only signatures that
# have not yet collided with the target.
#
# The key optimization is that once a collision has been found at prefix k,
# we do not need to search for another collision for that same signature
# and prefix.
#
# ============================================================================

def fast_constructive_search(
    target_i,
    target_j,
    primes,
    cache
):

    max_prefix = max(PREFIXES)

    target_signatures = build_target_signatures(
        target_i,
        target_j,
        cache
    )

    found = {
        name: {
            prefix: None
            for prefix in PREFIXES
        }
        for name in SIGNATURES
    }

    # ------------------------------------------------------------------------
    # For each signature we maintain a dictionary:
    #
    #   partial signature -> first pair producing it
    #
    # This allows O(1)-average collision lookup.
    # ------------------------------------------------------------------------

    tables = {
        name: {}
        for name in SIGNATURES
    }

    n = len(primes)

    pair_scans = 0

    # Once all eight signatures have collisions for a given prefix,
    # we don't need to keep processing that prefix.
    #
    # However, longer prefixes still require the pairs.
    # ------------------------------------------------------------------------

    required = {
        prefix: set(SIGNATURES)
        for prefix in PREFIXES
    }

    unresolved = {
        prefix: set(SIGNATURES)
        for prefix in PREFIXES
    }

    for i in range(n):

        for j in range(i + 1, n):

            pair_scans += 1

            if pair_scans > MAX_PAIR_SCANS:
                return found, pair_scans

            if (
                (i == target_i and j == target_j)
                or
                (i == target_j and j == target_i)
            ):
                continue

            current = {
                name: ()
                for name in SIGNATURES
            }

            for prefix in range(1, max_prefix + 1):

                m = MODULI[prefix - 1]

                local = local_signature_indices(
                    i,
                    j,
                    m,
                    cache
                )

                for name in SIGNATURES:

                    current[name] += (
                        local[name],
                    )

                if prefix not in PREFIXES:
                    continue

                for name in SIGNATURES:

                    if name not in unresolved[prefix]:
                        continue

                    key = current[name]

                    # --------------------------------------------------------
                    # Check whether this exact signature belongs to target.
                    # --------------------------------------------------------

                    if (
                        key
                        ==
                        target_signatures[name][prefix]
                    ):
                        found[name][prefix] = (
                            i,
                            j
                        )

                        unresolved[prefix].discard(
                            name
                        )

                        continue

                    # --------------------------------------------------------
                    # Store the first pair producing this signature.
                    #
                    # This is useful for diagnostic statistics and ensures
                    # future occurrences can be recognized as collisions.
                    # --------------------------------------------------------

                    if key not in tables[name]:
                        tables[name][key] = (
                            i,
                            j
                        )

                # ------------------------------------------------------------
                # If every signature has already produced a target collision
                # at this prefix, there is no reason to keep building tables
                # for that prefix.
                # ------------------------------------------------------------

    return found, pair_scans


# ============================================================================
# COLLISION DEPTH SUMMARY
# ============================================================================

def first_collision_depth(found):

    result = {}

    for name in SIGNATURES:

        depth = None

        for prefix in PREFIXES:

            if found[name][prefix] is not None:

                depth = prefix
                break

        result[name] = depth

    return result


# ============================================================================
# MAIN
# ============================================================================

print("=" * 78)
print("KAPPA EXPERIMENT 24")
print("HASH-INDEXED CONSTRUCTIVE COLLISION SEARCH")
print("NO CSV OUTPUT")
print("=" * 78)

print(f"random seed      = {SEED}")
print(f"R values         = {R_VALUES}")
print(f"targets          = {TARGETS}")
print(f"prime pool       = {PRIME_POOL_SIZE:,}")
print(
    f"prime interval   = "
    f"[{PRIME_LOW:,}, {PRIME_HIGH:,}]"
)
print(f"prefixes         = {PREFIXES}")
print()


# ============================================================================
# 1. MODULUS INVENTORY
# ============================================================================

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


# ============================================================================
# 2. THEORETICAL BASELINES
# ============================================================================

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


# ============================================================================
# 3. PRIME POOL
# ============================================================================

print("=" * 78)
print("3. PRIME POOL GENERATION")
print("=" * 78)

primes = generate_prime_pool()


# ============================================================================
# 4. CACHE
# ============================================================================

print("=" * 78)
print("4. PRIME RESIDUE / ORDER CACHE")
print("=" * 78)

cache = PrimeCache(primes)


# ============================================================================
# 5. TARGETS
# ============================================================================

targets = choose_targets(primes)


# ============================================================================
# AGGREGATE RESULTS
# ============================================================================

all_depths = {
    name: []
    for name in SIGNATURES
}

all_scans = []


# ============================================================================
# TARGET LOOP
# ============================================================================

for target_number, (target_i, target_j) in enumerate(
    targets,
    start=1
):

    p = primes[target_i]
    q = primes[target_j]

    n = p * q
    s = p + q
    delta = (p - q) ** 2

    print("=" * 78)
    print(f"TARGET {target_number}")
    print("=" * 78)

    print(f"p       = {p}")
    print(f"q       = {q}")
    print(f"n       = {n}")
    print(f"s       = {s}")
    print(f"Delta   = {delta}")
    print(f"n bits  = {n.bit_length()}")
    print()

    # ------------------------------------------------------------------------
    # Constructive search
    # ------------------------------------------------------------------------

    found, pair_scans = fast_constructive_search(
        target_i,
        target_j,
        primes,
        cache
    )

    all_scans.append(pair_scans)

    # ------------------------------------------------------------------------
    # Collision table
    # ------------------------------------------------------------------------

    print("CONSTRUCTIVE COLLISIONS")
    print("-" * 78)

    print(
        f"{'prefix':>8s}"
        f"{'m':>8s}"
        + "".join(
            f"{name:>22s}"
            for name in SIGNATURES
        )
    )

    for prefix in PREFIXES:

        row = [
            f"{prefix:8d}",
            f"{MODULI[prefix - 1]:8d}",
        ]

        for name in SIGNATURES:

            pair = found[name][prefix]

            if pair is None:
                text = "NONE"
            else:
                text = fmt_pair(
                    pair,
                    primes
                )

            row.append(
                f"{text:>22s}"
            )

        print("".join(row))

    print()

    # ------------------------------------------------------------------------
    # First collision depth
    # ------------------------------------------------------------------------

    depths = first_collision_depth(found)

    print("FIRST CONSTRUCTIVE COLLISION PREFIX")
    print("-" * 78)

    for name in SIGNATURES:

        depth = depths[name]

        if depth is None:
            text = "NONE FOUND"
        else:
            text = str(depth)

        print(
            f"{name:18s}"
            f"{text:>12s}"
        )

        all_depths[name].append(
            depth
        )

    print()

    print(
        f"pair scans = {pair_scans:,}"
    )

    if pair_scans >= MAX_PAIR_SCANS:
        print(
            "WARNING: MAX_PAIR_SCANS reached."
        )

    print()


# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print()
print("FIRST CONSTRUCTIVE COLLISION DEPTH")
print("-" * 78)

print(
    f"{'signature':18s}"
    + "".join(
        f"{prefix:>8d}"
        for prefix in PREFIXES
    )
)

for name in SIGNATURES:

    row = [f"{name:18s}"]

    for prefix in PREFIXES:

        count = sum(
            1
            for depth in all_depths[name]
            if depth is not None
            and depth <= prefix
        )

        row.append(
            f"{count:8d}"
        )

    print("".join(row))


print()
print("TARGETS WITH CONSTRUCTIVE COLLISION AT EACH PREFIX")
print("-" * 78)

for prefix in PREFIXES:

    print()
    print(f"PREFIX {prefix} (m={MODULI[prefix - 1]})")

    for name in SIGNATURES:

        count = sum(
            1
            for depth in all_depths[name]
            if depth == prefix
        )

        cumulative = sum(
            1
            for depth in all_depths[name]
            if depth is not None
            and depth <= prefix
        )

        print(
            f"  {name:18s}"
            f" first={count:3d}/{TARGETS}"
            f"  cumulative={cumulative:3d}/{TARGETS}"
        )


print()
print("SEARCH WORK")
print("-" * 78)

print(
    f"minimum pair scans = "
    f"{min(all_scans):,}"
)

print(
    f"median pair scans   = "
    f"{sorted(all_scans)[len(all_scans)//2]:,}"
)

print(
    f"maximum pair scans = "
    f"{max(all_scans):,}"
)


print()
print("=" * 78)
print("INTERPRETATION")
print("=" * 78)

print("""
Experiment 24 changes the computational question.

Experiment 23 repeatedly searched for collisions by examining many
candidate pairs. That becomes expensive because the prime-pool pair
space grows approximately as:

    N(N-1)/2

For N = 5,000 this is almost 12.5 million unordered pairs.

Experiment 24 uses hash-indexed constructive searching.

For every candidate pair, the partial signature is converted into a
tuple that can be used as a dictionary key.

The target signature is therefore tested by exact hash lookup rather
than by comparing the candidate against many previous candidates.

The important output is:

    FIRST CONSTRUCTIVE COLLISION PREFIX

This is stronger than a random collision count because the program
actually produces a concrete wrong prime pair whenever one exists
inside the searched pool.

Interpretation:

    NONE FOUND
        means no concrete collision was found in the searched pool
        before MAX_PAIR_SCANS.

    concrete pair
        means an explicit wrong prime pair reproduces the target
        signature at that prefix.

A collision is empirical evidence only. It does not prove that a
signature is mathematically non-unique globally.

Likewise, NONE FOUND does not prove uniqueness.

The particularly useful comparison is:

    Fpair
    Fsorted
        versus
    Fsum / Fprod / Fdiff / cube_sum
    order_pair / order_sorted

If Fpair repeatedly loses all constructive collisions at a shorter
prefix than the compressed signatures, that gives stronger evidence
that retaining the complete residue pair preserves more identifying
information.

If Fpair and Fsorted behave almost identically, that indicates that
ordering contributes little in this experiment.

If scalar signatures continue to produce concrete collisions after
Fpair has stopped producing them, that is direct constructive evidence
that the compression loses discriminating information.

IMPORTANT:
This experiment searches only the finite prime pool and only up to
MAX_PAIR_SCANS. It is not a proof of global uniqueness.
""")


print()
print("=" * 78)
print("DONE")
print("=" * 78)

