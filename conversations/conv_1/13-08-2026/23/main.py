#!/usr/bin/env python3

import math
import random
from collections import defaultdict

# ============================================================
# KAPPA EXPERIMENT 23
# FAST SIGNATURE COLLISION DATABASE / CONSTRUCTIVE SEARCH
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

# Prime generation interval.
PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

# Number of primes used for the collision database.
# Increase if the machine has enough RAM.
PRIME_POOL_SIZE = 5_000

# Prefixes around the observed transition.
PREFIXES = [3, 4, 5, 6, 7]

# Only these signatures need explicit collision databases.
SEARCH_SIGNATURES = [
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
# MODULI
# ============================================================

MODULI = [
    r * r + 3
    for r in R_VALUES
]


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


def generate_prime_pool(size):
    """
    Generate a unique pool of primes.

    The pool is generated once and then reused for every target.
    """

    primes = set()

    while len(primes) < size:

        n = random.randrange(
            PRIME_LOW,
            PRIME_HIGH
        )

        if n % 2 == 0:
            n += 1

        if is_prime(n):
            primes.add(n)

    return sorted(primes)


def random_prime_from_pool(pool):
    return random.choice(pool)


def factor_pair(pool):
    while True:

        p = random_prime_from_pool(pool)
        q = random_prime_from_pool(pool)

        if p != q:
            return p, q


# ============================================================
# INVENTORY
# ============================================================

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
            orders[u] = multiplicative_order(
                u,
                m
            )

        inventory[m] = {
            "units": units,
            "orbit": len(units),
            "orders": orders,
        }

    return inventory


INVENTORY = build_inventory()


# ============================================================
# PRIME INFORMATION CACHE
# ============================================================

def build_prime_info(primes):
    """
    Precompute all residue/order information needed by the
    experiment.

    This is the main performance optimization over Experiment 22.
    """

    info = {}

    for p in primes:

        residues = []
        orders = []

        for m in MODULI:

            r = p % m

            residues.append(r)

            orders.append(
                INVENTORY[m]["orders"].get(r)
            )

        info[p] = {
            "residues": tuple(residues),
            "orders": tuple(orders),
        }

    return info


# ============================================================
# SIGNATURE VALUE FROM CACHED PRIME DATA
# ============================================================

def signature_value(
    p,
    q,
    signature,
    prefix,
    prime_info
):

    ip = prime_info[p]
    iq = prime_info[q]

    values = []

    for i in range(prefix):

        m = MODULI[i]

        a = ip["residues"][i]
        b = iq["residues"][i]

        if signature == "Fpair":

            value = (
                a,
                b
            )

        elif signature == "Fsorted":

            value = tuple(
                sorted((a, b))
            )

        elif signature == "Fsum":

            value = (
                (a + b) % m
            )

        elif signature == "Fprod":

            value = (
                (a * b) % m
            )

        elif signature == "Fdiff":

            value = (
                (a - b) % m
            )

        elif signature == "cube_sum":

            value = (
                (a ** 3 + b ** 3) % m
            )

        elif signature == "order_pair":

            value = (
                ip["orders"][i],
                iq["orders"][i]
            )

        elif signature == "order_sorted":

            value = tuple(
                sorted(
                    (
                        ip["orders"][i],
                        iq["orders"][i]
                    )
                )
            )

        else:
            raise ValueError(signature)

        values.append(value)

    return tuple(values)


# ============================================================
# ALL SIGNATURES FOR ONE PAIR
# ============================================================

def all_signatures(
    p,
    q,
    prefix,
    prime_info
):

    return {
        name: signature_value(
            p,
            q,
            name,
            prefix,
            prime_info
        )
        for name in SEARCH_SIGNATURES
    }


# ============================================================
# BUILD COLLISION DATABASE
# ============================================================

def build_collision_database(
    primes,
    prime_info,
    signature,
    prefix
):
    """
    Build:

        signature -> first prime pair

    for the selected prefix.

    Only one representative pair is retained for each signature.
    """

    database = {}

    count = 0

    total = len(primes)

    for i in range(total):

        p = primes[i]

        for j in range(i + 1, total):

            q = primes[j]

            key = signature_value(
                p,
                q,
                signature,
                prefix,
                prime_info
            )

            if key not in database:
                database[key] = (
                    p,
                    q
                )

            count += 1

    return database


# ============================================================
# TARGET COLLISION SEARCH
# ============================================================

def find_collision(
    target_p,
    target_q,
    primes,
    prime_info,
    signature,
    prefix
):
    """
    Find a different prime pair with exactly the same signature
    as the target.

    Uses the database generated from the prime pool.
    """

    target_key = signature_value(
        target_p,
        target_q,
        signature,
        prefix,
        prime_info
    )

    database = {}

    total = len(primes)

    for i in range(total):

        p = primes[i]

        for j in range(i + 1, total):

            q = primes[j]

            pair = (
                p,
                q
            )

            if {
                p,
                q
            } == {
                target_p,
                target_q
            }:
                continue

            key = signature_value(
                p,
                q,
                signature,
                prefix,
                prime_info
            )

            if key == target_key:
                return (
                    p,
                    q
                )

            # Keep the database bounded by only storing the
            # signatures encountered during the search.
            database[key] = pair

    return None


# ============================================================
# FASTER COLLISION SEARCH USING HASH BUCKETS
# ============================================================

def collision_bucket_search(
    target_p,
    target_q,
    primes,
    prime_info,
    signature,
    prefix
):
    """
    Construct a hash bucket for the target signature.

    Instead of storing every signature globally, candidates are
    compared through a direct dictionary lookup.

    The target signature is calculated once.
    """

    target_key = signature_value(
        target_p,
        target_q,
        signature,
        prefix,
        prime_info
    )

    total = len(primes)

    for i in range(total):

        p = primes[i]

        for j in range(i + 1, total):

            q = primes[j]

            if {
                p,
                q
            } == {
                target_p,
                target_q
            }:
                continue

            key = signature_value(
                p,
                q,
                signature,
                prefix,
                prime_info
            )

            if key == target_key:

                return (
                    p,
                    q
                )

    return None


# ============================================================
# CONSTRUCTIVE RESIDUE SEARCH
# ============================================================

def residue_candidates(
    target_p,
    target_q,
    primes,
    prime_info,
    prefix
):
    """
    Search for primes having the same individual residue vector
    as target_p and target_q.

    This is especially useful for Fpair.

    Because Fpair is ordered, we search:

        p' == p residue vector
        q' == q residue vector

    separately.

    For Fsorted, the two residue vectors may also be exchanged.
    """

    target_p_residues = (
        prime_info[target_p]["residues"][:prefix]
    )

    target_q_residues = (
        prime_info[target_q]["residues"][:prefix]
    )

    p_candidates = []
    q_candidates = []

    for x in primes:

        residues = (
            prime_info[x]["residues"][:prefix]
        )

        if residues == target_p_residues:
            p_candidates.append(x)

        if residues == target_q_residues:
            q_candidates.append(x)

    return (
        p_candidates,
        q_candidates
    )


# ============================================================
# CONSTRUCTIVE FPAIR COLLISION
# ============================================================

def constructive_fpair_collision(
    target_p,
    target_q,
    primes,
    prime_info,
    prefix
):
    """
    Look for another ordered prime pair with exactly the same
    Fpair residue vector.

    This is much more targeted than random search.
    """

    p_candidates, q_candidates = residue_candidates(
        target_p,
        target_q,
        primes,
        prime_info,
        prefix
    )

    for p in p_candidates:

        for q in q_candidates:

            if p == q:
                continue

            if (
                p == target_p
                and q == target_q
            ):
                continue

            return (
                p,
                q
            )

    return None


# ============================================================
# CONSTRUCTIVE FSORTED COLLISION
# ============================================================

def constructive_fsorted_collision(
    target_p,
    target_q,
    primes,
    prime_info,
    prefix
):
    """
    Search for another unordered residue pair with the same
    Fsorted signature.
    """

    p_candidates, q_candidates = residue_candidates(
        target_p,
        target_q,
        primes,
        prime_info,
        prefix
    )

    # Ordered orientation.
    for p in p_candidates:

        for q in q_candidates:

            if p == q:
                continue

            if {
                p,
                q
            } == {
                target_p,
                target_q
            }:
                continue

            return (
                p,
                q
            )

    # Reversed orientation.
    for p in q_candidates:

        for q in p_candidates:

            if p == q:
                continue

            if {
                p,
                q
            } == {
                target_p,
                target_q
            }:
                continue

            return (
                p,
                q
            )

    return None


# ============================================================
# FORMAT HELPERS
# ============================================================

def pair_text(pair):

    if pair is None:
        return "NONE FOUND"

    return f"({pair[0]}, {pair[1]})"


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
# MAIN
# ============================================================

print("=" * 78)
print("KAPPA EXPERIMENT 23")
print("FAST SIGNATURE COLLISION DATABASE / CONSTRUCTIVE SEARCH")
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
# INVENTORY
# ============================================================

print("=" * 78)
print("1. MODULUS INVENTORY")
print("=" * 78)

for r in R_VALUES:

    m = r * r + 3

    print(
        f"r={r:3d} "
        f"m={m:6d} "
        f"units={INVENTORY[m]['orbit']:7d}"
    )

print()


# ============================================================
# THEORETICAL BASELINES
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

print(
    "Generating unique prime pool..."
)

PRIMES = generate_prime_pool(
    PRIME_POOL_SIZE
)

print(
    f"generated primes = {len(PRIMES):,}"
)

print()


# ============================================================
# PRIME INFORMATION
# ============================================================

print("=" * 78)
print("4. PRIME RESIDUE / ORDER CACHE")
print("=" * 78)

print(
    "Precomputing residues and multiplicative orders..."
)

PRIME_INFO = build_prime_info(
    PRIMES
)

print("cache complete")
print()


# ============================================================
# TARGETS
# ============================================================

all_results = defaultdict(
    lambda: defaultdict(int)
)


for target in range(
    1,
    TARGETS + 1
):

    p, q = factor_pair(
        PRIMES
    )

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
    # CONSTRUCTIVE FPAIR / FSORTED SEARCH
    # --------------------------------------------------------

    print(
        "CONSTRUCTIVE FPAIR / FSORTED SEARCH"
    )

    print("-" * 78)

    print(
        f"{'prefix':>8s}"
        f"{'m':>8s}"
        f"{'Fpair collision':>28s}"
        f"{'Fsorted collision':>28s}"
    )

    for prefix in PREFIXES:

        fp = constructive_fpair_collision(
            p,
            q,
            PRIMES,
            PRIME_INFO,
            prefix
        )

        fs = constructive_fsorted_collision(
            p,
            q,
            PRIMES,
            PRIME_INFO,
            prefix
        )

        if fp is not None:
            all_results[
                "Fpair"
            ][prefix] += 1

        if fs is not None:
            all_results[
                "Fsorted"
            ][prefix] += 1

        print(
            f"{prefix:8d}"
            f"{MODULI[prefix - 1]:8d}"
            f"{pair_text(fp):>28s}"
            f"{pair_text(fs):>28s}"
        )

    print()


    # --------------------------------------------------------
    # GENERAL SIGNATURE SEARCH
    # --------------------------------------------------------

    print(
        "GENERAL SIGNATURE COLLISION SEARCH"
    )

    print("-" * 78)

    print(
        "Searching the prime pool for concrete wrong-pair"
    )

    print(
        "collisions at prefixes 4, 5 and 6."
    )

    print()

    for prefix in [4, 5, 6]:

        print(
            f"PREFIX {prefix} "
            f"(m={MODULI[prefix - 1]})"
        )

        for name in SEARCH_SIGNATURES:

            # Fpair and Fsorted have already been handled by
            # the constructive search.
            if name == "Fpair":

                pair = constructive_fpair_collision(
                    p,
                    q,
                    PRIMES,
                    PRIME_INFO,
                    prefix
                )

            elif name == "Fsorted":

                pair = constructive_fsorted_collision(
                    p,
                    q,
                    PRIMES,
                    PRIME_INFO,
                    prefix
                )

            else:

                pair = collision_bucket_search(
                    p,
                    q,
                    PRIMES,
                    PRIME_INFO,
                    name,
                    prefix
                )

            if pair is None:

                print(
                    f"  {name:18s}"
                    f"NONE FOUND"
                )

            else:

                all_results[
                    name
                ][prefix] += 1

                print(
                    f"  {name:18s}"
                    f"{pair_text(pair)}"
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
# CONCRETE COLLISION COUNTS
# ============================================================

print(
    "TARGETS WITH CONCRETE COLLISIONS"
)

print("-" * 78)

print(
    f"{'signature':18s}"
    f"{'prefix 3':>12s}"
    f"{'prefix 4':>12s}"
    f"{'prefix 5':>12s}"
    f"{'prefix 6':>12s}"
    f"{'prefix 7':>12s}"
)

for name in SEARCH_SIGNATURES:

    print(
        f"{name:18s}"
        f"{all_results[name][3]:12d}"
        f"{all_results[name][4]:12d}"
        f"{all_results[name][5]:12d}"
        f"{all_results[name][6]:12d}"
        f"{all_results[name][7]:12d}"
    )

print()


# ============================================================
# THRESHOLD SUMMARY
# ============================================================

print(
    "EMPIRICAL COLLISION THRESHOLD"
)

print("-" * 78)

for name in SEARCH_SIGNATURES:

    first_zero = None

    for prefix in PREFIXES:

        if all_results[name][prefix] == 0:

            first_zero = prefix

            break

    if first_zero is None:

        print(
            f"{name:18s}"
            f"no zero-collision prefix"
        )

    else:

        print(
            f"{name:18s}"
            f"first zero observed = {first_zero}"
        )

print()


# ============================================================
# INTERPRETATION
# ============================================================

print("=" * 78)
print("INTERPRETATION")
print("=" * 78)

print("""
Experiment 23 changes the collision experiment from repeated
Monte-Carlo signature calculation to a cached prime-pool search.

The prime pool is generated once.

For every prime, all residues modulo the selected moduli and all
multiplicative orders are precomputed once.

This makes subsequent signature comparisons substantially cheaper.

The most important experiment is the CONSTRUCTIVE FPAIR search.

For a target (p,q), another pair (p',q') is a genuine Fpair
collision at prefix k when:

    p' mod m_i = p mod m_i
    q' mod m_i = q mod m_i

for every modulus in the prefix.

The search therefore directly tests whether the residue fingerprint
can occur for another pair of primes in the selected prime population.

Fsorted additionally allows the two residue vectors to be exchanged.

The general signature search then asks the same question for:

    Fsum
    Fprod
    Fdiff
    cube_sum
    order_pair
    order_sorted

A CONCRETE COLLISION is much stronger evidence than zero collisions
in a random sample: it gives an explicit wrong prime pair that
survives the tested signature.

However, failure to find a collision is still NOT a mathematical
proof of uniqueness. It only establishes that no collision was
found inside the selected prime pool.

The critical comparison is therefore:

    Fpair prefix 4
    versus
    compressed-signature prefixes 5 and 6.

If Fpair has no concrete collisions at prefix 4 while compressed
signatures still have many concrete collisions, that strengthens
the collision-depth separation observed in Experiments 21 and 22.

If a concrete Fpair collision appears at prefix 4 or 5, that gives
us an explicit counterexample and should become the target of the
next experiment.
""")

print()
print("=" * 78)
print("DONE")
print("=" * 78)

