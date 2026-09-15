import random
import time
import sympy


# ==============================================================================
# KAPPA EXPERIMENT 31
# PRUNED MINIMAL SIGNATURE SEPARATION
# NO CSV OUTPUT
# ==============================================================================

SEED = 20260814

TARGETS = 12
POOL_SIZE = 5_000

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

# Cheap signatures first.
#
# order_pair is deliberately last because multiplicative order is
# substantially more expensive than the residue-based signatures.
SIGNATURES = [
    "Fpair",
    "Fdiff",
    "cube_sum",
    "Fsum",
    "Fprod",
    "Fsorted",
    "order_pair",
]


# ==============================================================================
# OUTPUT
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


# ==============================================================================
# PRIME GENERATION
# ==============================================================================

def generate_prime_pool():
    rng = random.Random(SEED)

    primes = set()

    while len(primes) < POOL_SIZE:

        x = rng.randint(PRIME_LO, PRIME_HI)

        if x < 2:
            continue

        if x != 2 and x % 2 == 0:
            x += 1

        if x >= PRIME_HI:
            continue

        if sympy.isprime(x):
            primes.add(x)

    return sorted(primes)


# ==============================================================================
# TARGET GENERATION
# ==============================================================================

def generate_targets(primes):
    """
    Deterministic target generation.

    The target pairs are selected from the same generated prime pool.
    """

    rng = random.Random(SEED)

    targets = []
    used = set()

    while len(targets) < TARGETS:

        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        pair = tuple(sorted((p, q)))

        if pair in used:
            continue

        used.add(pair)
        targets.append(pair)

    return targets


# ==============================================================================
# MODULUS
# ==============================================================================

def build_modulus_inventory():
    """
    Same modulus construction used by the compact experiment.

    m = 3r + 1
    """

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
# SIGNATURE CALCULATIONS
# ==============================================================================

def cheap_signatures(p, q, m):
    """
    Calculate the cheap residue signatures.

    No multiplicative orders here.
    """

    rp = p % m
    rq = q % m

    return {
        "Fpair": (rp, rq),

        "Fsorted": (
            min(rp, rq),
            max(rp, rq),
        ),

        "Fdiff": (rp - rq) % m,

        "Fsum": (rp + rq) % m,

        "Fprod": (rp * rq) % m,

        "cube_sum": (
            pow(rp, 3, m) +
            pow(rq, 3, m)
        ) % m,
    }


def order_signature(p, q, m):
    """
    Calculate multiplicative-order signature.

    Only called after the cheaper signatures have survived.
    """

    rp = p % m
    rq = q % m

    if sympy.gcd(rp, m) != 1:
        op = 0
    else:
        op = sympy.n_order(rp, m)

    if sympy.gcd(rq, m) != 1:
        oq = 0
    else:
        oq = sympy.n_order(rq, m)

    return tuple(sorted((op, oq)))


# ==============================================================================
# TARGET SIGNATURES
# ==============================================================================

def get_target_signatures(p, q, m):
    """
    Calculate all target signature values once.
    """

    values = cheap_signatures(p, q, m)

    values["order_pair"] = order_signature(p, q, m)

    return values


# ==============================================================================
# CANDIDATE TESTING
# ==============================================================================

def matches_signature(
    p,
    q,
    m,
    signature,
    target_values,
):
    """
    Test one candidate against one target signature.

    This intentionally calculates only the requested signature.
    """

    rp = p % m
    rq = q % m

    if signature == "Fpair":
        return (rp, rq) == target_values["Fpair"]

    if signature == "Fsorted":
        return (
            min(rp, rq),
            max(rp, rq)
        ) == target_values["Fsorted"]

    if signature == "Fdiff":
        return (
            (rp - rq) % m
            ==
            target_values["Fdiff"]
        )

    if signature == "Fsum":
        return (
            (rp + rq) % m
            ==
            target_values["Fsum"]
        )

    if signature == "Fprod":
        return (
            (rp * rq) % m
            ==
            target_values["Fprod"]
        )

    if signature == "cube_sum":
        return (
            (
                pow(rp, 3, m)
                +
                pow(rq, 3, m)
            ) % m
            ==
            target_values["cube_sum"]
        )

    if signature == "order_pair":
        return (
            order_signature(p, q, m)
            ==
            target_values["order_pair"]
        )

    raise ValueError(f"Unknown signature: {signature}")


# ==============================================================================
# TARGET INDEX
# ==============================================================================

def find_target_index(primes, target):
    """
    Locate target components in the prime pool.
    """

    p, q = target

    try:
        pi = primes.index(p)
        qi = primes.index(q)
    except ValueError:
        raise RuntimeError(
            f"Target {target} is not present in prime pool."
        )

    return pi, qi


# ==============================================================================
# INITIAL CANDIDATES
# ==============================================================================

def initial_candidates(primes, target):
    """
    Generate unordered candidate pairs.

    We keep the pair as indexes into the prime list.

    Important optimization:
    the pair is represented by (i, j), not by copying prime values.
    """

    n = len(primes)

    return [
        (i, j)
        for i in range(n)
        for j in range(i + 1, n)
    ]


# ==============================================================================
# PRUNED SEARCH
# ==============================================================================

def run_target_prefix(
    primes,
    target,
    m,
):
    """
    Search one target at one prefix.

    The search is intentionally sequential.

    First use cheap signatures to eliminate candidates.
    Only survivors reach expensive signatures.
    """

    p_target, q_target = target

    target_values = get_target_signatures(
        p_target,
        q_target,
        m,
    )

    target_i, target_j = find_target_index(
        primes,
        target,
    )

    # --------------------------------------------------------------------------
    # Instead of constructing all 12.5 million pair tuples,
    # maintain candidates as pairs of prime indexes.
    #
    # We start with the full triangular search space.
    # --------------------------------------------------------------------------

    candidates = None

    # We benchmark the cheap signatures on a small sample first.
    #
    # This lets us determine which signature gives the strongest pruning
    # without assuming the best order.
    sample_limit = 10_000

    sample = []

    count = 0

    for i in range(len(primes)):

        p = primes[i]

        for j in range(i + 1, len(primes)):

            q = primes[j]

            sample.append((p, q))

            count += 1

            if count >= sample_limit:
                break

        if count >= sample_limit:
            break

    # --------------------------------------------------------------------------
    # Benchmark cheap signatures.
    # --------------------------------------------------------------------------

    benchmark = []

    cheap = [
        "Fpair",
        "Fdiff",
        "cube_sum",
        "Fsum",
        "Fprod",
        "Fsorted",
    ]

    for signature in cheap:

        start = time.perf_counter()

        survivors = 0

        for p, q in sample:

            if matches_signature(
                p,
                q,
                m,
                signature,
                target_values,
            ):
                survivors += 1

        elapsed = time.perf_counter() - start

        benchmark.append(
            (
                signature,
                survivors,
                elapsed,
            )
        )

    # --------------------------------------------------------------------------
    # Rank by pruning efficiency first, then speed.
    #
    # A lower survivor ratio is better.
    # --------------------------------------------------------------------------

    benchmark.sort(
        key=lambda x: (
            x[1] / len(sample),
            x[2],
        )
    )

    signature_order = [
        x[0]
        for x in benchmark
    ]

    # order_pair stays last.
    signature_order.append("order_pair")

    # --------------------------------------------------------------------------
    # Actual candidate scan.
    #
    # Instead of generating every pair and then calculating seven signatures,
    # calculate signatures one at a time and immediately reject failures.
    # --------------------------------------------------------------------------

    survivors_by_step = []

    # For the first signature, scan the whole pool.
    #
    # Subsequent signatures only see survivors.
    #
    # To avoid storing millions of tuples, use integer encoded pairs:
    #
    # encoded = i * n + j
    #
    # This is considerably smaller than storing two Python integers in a tuple.
    n = len(primes)

    current = None

    for step, signature in enumerate(signature_order, 1):

        start = time.perf_counter()

        next_candidates = []

        if current is None:

            # First pass over all unordered prime pairs.
            for i in range(n):

                p = primes[i]

                for j in range(i + 1, n):

                    q = primes[j]

                    if matches_signature(
                        p,
                        q,
                        m,
                        signature,
                        target_values,
                    ):
                        next_candidates.append(
                            i * n + j
                        )

        else:

            for encoded in current:

                i = encoded // n
                j = encoded - i * n

                p = primes[i]
                q = primes[j]

                if matches_signature(
                    p,
                    q,
                    m,
                    signature,
                    target_values,
                ):
                    next_candidates.append(encoded)

        elapsed = time.perf_counter() - start

        current = next_candidates

        # Target must always survive.
        target_encoded = target_i * n + target_j

        target_survives = target_encoded in current

        survivors_by_step.append(
            (
                signature,
                len(current),
                target_survives,
                elapsed,
            )
        )

        if len(current) <= 1:
            break

        if not target_survives:
            break

    return survivors_by_step


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    experiment_start = time.perf_counter()

    banner(
        "KAPPA EXPERIMENT 31\n"
        "PRUNED MINIMAL SIGNATURE SEPARATION\n"
        "NO CSV OUTPUT"
    )

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {POOL_SIZE:,}")
    print(
        f"prime interval   = "
        f"[{PRIME_LO:,}, {PRIME_HI:,}]"
    )
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # --------------------------------------------------------------------------
    # MODULI
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

        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={p * q}"
        )

    # --------------------------------------------------------------------------
    # EXPERIMENT
    # --------------------------------------------------------------------------

    banner("PRUNED MINIMAL SIGNATURE SEPARATION")

    global_stats = {
        prefix: {
            "unique": 0,
            "nonunique": 0,
            "failed": 0,
        }
        for prefix in PREFIXES
    }

    for target_no, target in enumerate(targets, 1):

        section(
            f"TARGET {target_no:2d} "
            f"({target[0]}, {target[1]})"
        )

        for prefix, (_, m, _) in zip(
            PREFIXES,
            inventory,
        ):

            start = time.perf_counter()

            result = run_target_prefix(
                primes,
                target,
                m,
            )

            elapsed = time.perf_counter() - start

            final_count = result[-1][1]

            order = [
                x[0]
                for x in result
            ]

            target_survives = result[-1][2]

            if not target_survives:

                status = "TARGET LOST"
                global_stats[prefix]["failed"] += 1

            elif final_count == 1:

                status = "UNIQUE"
                global_stats[prefix]["unique"] += 1

            else:

                status = f"{final_count} survivors"
                global_stats[prefix]["nonunique"] += 1

            counts = [
                str(x[1])
                for x in result
            ]

            print(
                f"prefix {prefix}: "
                f"{' -> '.join(counts):<35} "
                f"{status:<18} "
                f"{' -> '.join(order):<50} "
                f"[{elapsed:.2f}s]"
            )

    # --------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # --------------------------------------------------------------------------

    banner("GLOBAL SUMMARY")

    print(
        "prefix      unique      nonunique      target-lost"
    )
    print(
        "---------------------------------------------------"
    )

    for prefix in PREFIXES:

        stats = global_stats[prefix]

        print(
            f"{prefix:6d}"
            f"{stats['unique']:12d}"
            f"{stats['nonunique']:14d}"
            f"{stats['failed']:16d}"
        )

    # --------------------------------------------------------------------------
    # FINAL
    # --------------------------------------------------------------------------

    elapsed = time.perf_counter() - experiment_start

    banner("EXPERIMENT 31 COMPLETE")

    print(
        "UNIQUE     = exactly one surviving candidate in the finite pool"
    )
    print(
        "NONUNIQUE  = two or more candidates remain"
    )
    print(
        "TARGET LOST = diagnostic failure; signature direction/order"
        " excluded the target"
    )
    print(
        "NONE of these results proves global uniqueness."
    )

    print()
    print(
        f"total runtime = {elapsed:.2f}s "
        f"({elapsed / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

