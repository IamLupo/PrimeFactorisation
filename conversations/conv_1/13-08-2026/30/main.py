import random
import time
import sympy


# ==============================================================================
# KAPPA EXPERIMENT 30
# MINIMAL SIGNATURE SEPARATION — COMPACT
# NO CSV OUTPUT
# ==============================================================================

SEED = 20260814
TARGETS = 12
POOL_SIZE = 5000

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

SIGNATURES = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "order_pair",
]


# ==============================================================================
# HELPERS
# ==============================================================================

def banner(text):
    print()
    print("=" * 78)
    print(text)
    print("=" * 78)


def section(text):
    print()
    print("-" * 78)
    print(text)
    print("-" * 78)


def generate_prime_pool():
    """
    Generate exactly POOL_SIZE unique primes in the requested interval.
    """
    rng = random.Random(SEED)

    primes = set()

    while len(primes) < POOL_SIZE:
        x = rng.randint(PRIME_LO, PRIME_HI)

        if x % 2 == 0:
            x += 1

        if x >= PRIME_HI:
            continue

        if sympy.isprime(x):
            primes.add(x)

    return sorted(primes)


def generate_targets(primes):
    """
    Reproduce the fixed target construction from the previous experiment.

    The target pairs are selected deterministically from the same seed.
    """
    rng = random.Random(SEED)

    targets = []

    while len(targets) < TARGETS:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        pair = tuple(sorted((p, q)))

        if pair not in targets:
            targets.append(pair)

    return targets


def build_modulus_inventory():
    inventory = []

    for r in R_VALUES:
        m = 3 * r + 1

        units = sum(
            1
            for x in range(1, m)
            if sympy.gcd(x, m) == 1
        )

        inventory.append((r, m, units))

    return inventory


# ==============================================================================
# SIGNATURE FUNCTIONS
# ==============================================================================

def signature_values(p, q, m):
    """
    Compute the seven signatures for a pair.

    This is intentionally kept in one place so the experiment can be
    changed without changing the search machinery.
    """

    rp = p % m
    rq = q % m

    # Pair of residues.
    Fpair = (rp, rq)

    # Sorted pair.
    Fsorted = tuple(sorted((rp, rq)))

    # Sum modulo m.
    Fsum = (rp + rq) % m

    # Product modulo m.
    Fprod = (rp * rq) % m

    # Difference modulo m.
    Fdiff = (rp - rq) % m

    # Cubes modulo m.
    cube_sum = (pow(rp, 3, m) + pow(rq, 3, m)) % m

    # Multiplicative orders.
    if sympy.gcd(rp, m) == 1:
        op = sympy.n_order(rp, m)
    else:
        op = 0

    if sympy.gcd(rq, m) == 1:
        oq = sympy.n_order(rq, m)
    else:
        oq = 0

    order_pair = tuple(sorted((op, oq)))

    return {
        "Fpair": Fpair,
        "Fsorted": Fsorted,
        "Fsum": Fsum,
        "Fprod": Fprod,
        "Fdiff": Fdiff,
        "cube_sum": cube_sum,
        "order_pair": order_pair,
    }


# ==============================================================================
# PREFIX CACHE
# ==============================================================================

def build_pair_signature_cache(primes, m):
    """
    Compute all signatures for every unordered prime pair.

    Only the signature tuple is retained.
    """

    cache = []

    for i in range(len(primes)):
        p = primes[i]

        for j in range(i + 1, len(primes)):
            q = primes[j]

            values = signature_values(p, q, m)

            cache.append(
                (
                    p,
                    q,
                    values,
                )
            )

    return cache


# ==============================================================================
# CANDIDATE SETS
# ==============================================================================

def build_signature_sets(pair_cache):
    """
    Convert the pair cache into:

        signature -> value -> set(pair index)

    This allows very cheap intersections later.
    """

    indexes = {
        name: {}
        for name in SIGNATURES
    }

    for idx, (_, _, values) in enumerate(pair_cache):

        for name in SIGNATURES:

            value = values[name]

            bucket = indexes[name].setdefault(value, set())

            bucket.add(idx)

    return indexes


def target_candidate_set(pair_cache, indexes, target, signature):
    """
    Return all pair indexes matching the target under one signature.
    """

    p, q = target

    for idx, (a, b, _) in enumerate(pair_cache):
        if a == p and b == q:
            target_idx = idx
            break
    else:
        raise RuntimeError("Target pair not found in pair cache.")

    target_value = pair_cache[target_idx][2][signature]

    return set(indexes[signature][target_value])


# ==============================================================================
# GREEDY MINIMAL SEPARATION
# ==============================================================================

def greedy_separation(pair_cache, indexes, target):
    """
    Greedily choose the signature which produces the smallest surviving
    candidate set at every step.

    Important:
        The target itself must remain in the set.

    Returns:
        selected signatures
        final survivor count
        target survives
        step counts
    """

    remaining = set(range(len(pair_cache)))

    selected = []
    unused = set(SIGNATURES)

    target_idx = None

    for idx, (p, q, _) in enumerate(pair_cache):
        if (p, q) == target:
            target_idx = idx
            break

    if target_idx is None:
        raise RuntimeError("Target missing from candidate pool.")

    steps = []

    while unused:

        best_signature = None
        best_set = None
        best_size = None

        for signature in unused:

            target_matches = target_candidate_set(
                pair_cache,
                indexes,
                target,
                signature,
            )

            candidate = remaining & target_matches

            if target_idx not in candidate:
                continue

            size = len(candidate)

            if best_size is None or size < best_size:
                best_size = size
                best_signature = signature
                best_set = candidate

        if best_signature is None:
            break

        remaining = best_set
        selected.append(best_signature)
        unused.remove(best_signature)

        steps.append(
            (
                best_signature,
                len(remaining),
            )
        )

        if len(remaining) <= 1:
            break

    return selected, remaining, steps


# ==============================================================================
# COMPACT PREFIX EXPERIMENT
# ==============================================================================

def run_prefix(primes, targets, m, prefix):
    """
    Run one prefix and print only useful information.
    """

    start = time.perf_counter()

    pair_cache = build_pair_signature_cache(primes, m)

    cache_time = time.perf_counter() - start

    start = time.perf_counter()

    indexes = build_signature_sets(pair_cache)

    index_time = time.perf_counter() - start

    print(
        f"prefix {prefix}: "
        f"pairs={len(pair_cache):,} "
        f"cache={cache_time:.2f}s "
        f"index={index_time:.2f}s"
    )

    zero_count = 0
    nonzero_count = 0

    for target_no, target in enumerate(targets, 1):

        selected, remaining, steps = greedy_separation(
            pair_cache,
            indexes,
            target,
        )

        final_count = len(remaining)

        if final_count == 1:
            zero_count += 1
            status = "UNIQUE"
        else:
            nonzero_count += 1
            status = f"{final_count} survivors"

        print(
            f"  target {target_no:2d}: "
            f"{status:<18} "
            f"steps={len(selected)} "
            f"order={' -> '.join(selected)}"
        )

    print(
        f"  summary: "
        f"unique={zero_count}/{len(targets)} "
        f"nonzero={nonzero_count}/{len(targets)}"
    )

    return zero_count, nonzero_count


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    experiment_start = time.perf_counter()

    random.seed(SEED)

    banner(
        "KAPPA EXPERIMENT 30\n"
        "MINIMAL SIGNATURE SEPARATION — COMPACT\n"
        "NO CSV OUTPUT"
    )

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {POOL_SIZE:,}")
    print(f"prime interval   = [{PRIME_LO:,}, {PRIME_HI:,}]")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # --------------------------------------------------------------------------
    # MODULUS INVENTORY
    # --------------------------------------------------------------------------

    section("MODULUS INVENTORY")

    inventory = build_modulus_inventory()

    for r, m, units in inventory:
        print(
            f"r={r:3d} "
            f"m={m:7d} "
            f"units={units:7d}"
        )

    # --------------------------------------------------------------------------
    # PRIME POOL
    # --------------------------------------------------------------------------

    section("PRIME POOL")

    start = time.perf_counter()

    print("Generating prime pool...")

    primes = generate_prime_pool()

    print(
        f"generated primes = {len(primes):,}"
    )

    print(
        f"generation time = "
        f"{time.perf_counter() - start:.2f}s"
    )

    # --------------------------------------------------------------------------
    # TARGETS
    # --------------------------------------------------------------------------

    section("TARGETS")

    targets = generate_targets(primes)

    for i, (p, q) in enumerate(targets, 1):

        n = p * q

        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    # --------------------------------------------------------------------------
    # EXPERIMENT
    # --------------------------------------------------------------------------

    banner("MINIMAL SIGNATURE SEPARATION")

    global_results = {}

    for prefix, (_, m, _) in zip(PREFIXES, inventory):

        print()

        zero, nonzero = run_prefix(
            primes,
            targets,
            m,
            prefix,
        )

        global_results[prefix] = (zero, nonzero)

    # --------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # --------------------------------------------------------------------------

    banner("GLOBAL SUMMARY")

    print(
        "prefix     unique       nonzero"
    )
    print(
        "--------------------------------"
    )

    for prefix in PREFIXES:

        unique, nonzero = global_results[prefix]

        print(
            f"{prefix:6d} "
            f"{unique:10d} "
            f"{nonzero:12d}"
        )

    # --------------------------------------------------------------------------
    # END
    # --------------------------------------------------------------------------

    elapsed = time.perf_counter() - experiment_start

    banner("EXPERIMENT 30 COMPLETE")

    print(
        "unique = target isolated to itself within the finite searchable pool"
    )
    print(
        "nonzero = one or more wrong candidates remain"
    )
    print(
        "neither result proves global uniqueness"
    )

    print()
    print(
        f"total runtime = {elapsed:.2f}s "
        f"({elapsed / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

