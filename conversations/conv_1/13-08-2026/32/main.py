import random
import time
from collections import defaultdict

from sympy import isprime


# ==============================================================================
# KAPPA EXPERIMENT 32
# FAST PRUNED MINIMAL SIGNATURE SEPARATION
# NO CSV OUTPUT
# ==============================================================================
#
# Main optimization:
#
#   OLD:
#       scan all prime pairs for every signature
#
#   NEW:
#       build candidate set from ONE strong signature,
#       then test ONLY those candidates against the remaining signatures.
#
# Additional pruning:
#   - Fsorted is not independently computed; it is derived from Fpair.
#   - no pairwise intersection matrix
#   - no complete collision lists unless survivor count is small
#   - stop immediately when target is isolated
#
# ==============================================================================


SEED = 20260814

TARGETS = 12
POOL_SIZE = 5_000

LOW = 2_000_000
HIGH = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]


# ==============================================================================
# HELPERS
# ==============================================================================

def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def section(title):
    print()
    print("-" * 78)
    print(title)
    print("-" * 78)


def elapsed(t0):
    return f"{time.perf_counter() - t0:.3f}s"


# ==============================================================================
# MODULUS CONSTRUCTION
# ==============================================================================

def build_modulus_inventory():
    inventory = []

    for r in R_VALUES:
        m = 3 * r + 1

        units = sum(
            1
            for x in range(1, m)
            if __import__("math").gcd(x, m) == 1
        )

        inventory.append((r, m, units))

    return inventory


# ==============================================================================
# PRIME POOL
# ==============================================================================

def generate_prime_pool():
    rng = random.Random(SEED)

    primes = []
    seen = set()

    while len(primes) < POOL_SIZE:
        x = rng.randint(LOW, HIGH)

        if x in seen:
            continue

        if isprime(x):
            seen.add(x)
            primes.append(x)

    primes.sort()

    return primes


# ==============================================================================
# FIXED TARGET GENERATION
# ==============================================================================

def generate_targets(primes):
    #
    # Keep target generation deterministic.
    #
    rng = random.Random(SEED)

    targets = []

    while len(targets) < TARGETS:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair not in targets:
            targets.append(pair)

    return targets


# ==============================================================================
# SIGNATURE CORE
# ==============================================================================
#
# IMPORTANT:
#
# These functions preserve the signature definitions used by the previous
# experiments. If your previous script has slightly different formulas for
# these functions, replace ONLY these functions. The pruning engine below
# remains unchanged.
#
# ==============================================================================


def residue(x, m):
    return x % m


def Fpair(p, q, m):
    return (
        residue(p, m),
        residue(q, m),
    )


def Fsorted(p, q, m):
    a = residue(p, m)
    b = residue(q, m)
    return tuple(sorted((a, b)))


def Fsum(p, q, m):
    return (p + q) % m


def Fprod(p, q, m):
    return (p * q) % m


def Fdiff(p, q, m):
    return (p - q) % m


def cube_sum(p, q, m):
    return (pow(p, 3, m) + pow(q, 3, m)) % m


def multiplicative_order(a, m):
    #
    # Lazy order computation.
    #
    # If gcd(a,m) != 1, there is no multiplicative order.
    #
    import math

    a %= m

    if math.gcd(a, m) != 1:
        return 0

    order = 1
    value = a % m

    while value != 1:
        value = (value * a) % m
        order += 1

    return order


def order_pair(p, q, m):
    return (
        multiplicative_order(p, m),
        multiplicative_order(q, m),
    )


# ==============================================================================
# TARGET SIGNATURES
# ==============================================================================

def target_signature(target, m, signature):
    p, q = target

    if signature == "Fpair":
        return Fpair(p, q, m)

    if signature == "Fsorted":
        return Fsorted(p, q, m)

    if signature == "Fsum":
        return Fsum(p, q, m)

    if signature == "Fprod":
        return Fprod(p, q, m)

    if signature == "Fdiff":
        return Fdiff(p, q, m)

    if signature == "cube_sum":
        return cube_sum(p, q, m)

    if signature == "order_pair":
        return order_pair(p, q, m)

    raise ValueError(signature)


# ==============================================================================
# SIGNATURE EVALUATION
# ==============================================================================

def signature_value(pair, m, signature):
    p, q = pair

    if signature == "Fpair":
        return Fpair(p, q, m)

    if signature == "Fsorted":
        return Fsorted(p, q, m)

    if signature == "Fsum":
        return Fsum(p, q, m)

    if signature == "Fprod":
        return Fprod(p, q, m)

    if signature == "Fdiff":
        return Fdiff(p, q, m)

    if signature == "cube_sum":
        return cube_sum(p, q, m)

    if signature == "order_pair":
        return order_pair(p, q, m)

    raise ValueError(signature)


# ==============================================================================
# FAST CANDIDATE INDEX
# ==============================================================================

def build_signature_index(primes, m, signature):
    """
    Build:
        signature value -> list of prime pairs

    IMPORTANT:
    We store only unordered pairs p < q.
    """

    index = defaultdict(list)

    n = len(primes)

    for i in range(n):
        p = primes[i]

        for j in range(i + 1, n):
            q = primes[j]

            sig = signature_value((p, q), m, signature)

            index[sig].append((p, q))

    return index


# ==============================================================================
# FASTER SPECIAL INDEXES
# ==============================================================================

def build_Fpair_index(primes, m):
    """
    Fpair is especially cheap.

    We first group primes by residue and then create pairs only inside the
    required residue classes.

    This avoids evaluating Fpair for every pair individually.
    """

    buckets = defaultdict(list)

    for p in primes:
        buckets[p % m].append(p)

    index = defaultdict(list)

    for r1, values1 in buckets.items():

        # Same residue.
        for i in range(len(values1)):
            for j in range(i + 1, len(values1)):
                index[(r1, r1)].append(
                    (values1[i], values1[j])
                )

        # Different residues.
        for r2 in buckets:
            if r2 <= r1:
                continue

            values2 = buckets[r2]

            key = (r1, r2)

            for p in values1:
                for q in values2:
                    if p < q:
                        index[key].append((p, q))
                    else:
                        index[key].append((q, p))

    return index


# ==============================================================================
# FILTER
# ==============================================================================

def filter_candidates(candidates, target, m, signature):
    """
    Keep candidates matching the target's signature.

    The target itself is always in candidates when the initial candidate
    generation is correct.
    """

    wanted = target_signature(target, m, signature)

    result = []

    for pair in candidates:
        if signature_value(pair, m, signature) == wanted:
            result.append(pair)

    return result


# ==============================================================================
# SPECIAL FILTERS
# ==============================================================================

def filter_Fsorted_from_Fpair(candidates, target, m):
    """
    Fsorted contains no information beyond Fpair.

    If Fpair was already applied, Fsorted can never reduce the set.
    """

    return candidates


# ==============================================================================
# CHOOSE INITIAL SIGNATURE
# ==============================================================================

def choose_initial_signature(prefix):
    """
    Strong signatures are preferred.

    At small prefixes Fdiff is cheap and frequently useful.

    At larger prefixes Fpair is extremely strong.

    The important point is that we NEVER start with an expensive signature
    such as order_pair unless necessary.
    """

    if prefix >= 5:
        return "Fpair"

    if prefix == 4:
        return "Fpair"

    return "Fpair"


# ==============================================================================
# GREEDY ORDER
# ==============================================================================

def signature_order(prefix):
    """
    Cheap/strong first.

    Fsorted is deliberately omitted because it is redundant with Fpair.
    """

    if prefix >= 5:
        return [
            "Fpair",
            "Fdiff",
            "cube_sum",
            "order_pair",
            "Fsum",
            "Fprod",
        ]

    if prefix == 4:
        return [
            "Fpair",
            "Fdiff",
            "cube_sum",
            "order_pair",
            "Fsum",
            "Fprod",
        ]

    return [
        "Fpair",
        "Fdiff",
        "cube_sum",
        "order_pair",
        "Fsum",
        "Fprod",
    ]


# ==============================================================================
# TARGET SEPARATION
# ==============================================================================

def separate_target(primes, target, prefix, modulus):
    """
    Fast greedy separation.

    We deliberately do NOT build every signature index.

    First Fpair creates the initial candidate set.
    Every later signature scans only that candidate set.
    """

    t0 = time.perf_counter()

    m = modulus

    order = signature_order(prefix)

    print()
    print(
        f"prefix {prefix}: target={target} m={m}"
    )

    # ------------------------------------------------------------------
    # FIRST FILTER
    # ------------------------------------------------------------------

    first = order[0]

    if first == "Fpair":
        index = build_Fpair_index(primes, m)
        wanted = target_signature(target, m, first)

        candidates = list(index.get(wanted, []))

        # Release large dictionary as soon as possible.
        del index

    else:
        index = build_signature_index(primes, m, first)
        wanted = target_signature(target, m, first)
        candidates = list(index.get(wanted, []))
        del index

    print(
        f"  {first:12s}: {len(candidates):8d}"
    )

    if target not in candidates:
        print("  WARNING: target missing from initial candidate set")
        return {
            "survivors": candidates,
            "order": [first],
            "isolated": False,
            "target_survives": False,
            "time": time.perf_counter() - t0,
        }

    # ------------------------------------------------------------------
    # GREEDY FILTERS
    # ------------------------------------------------------------------

    used = [first]

    for signature in order[1:]:

        # Fsorted is redundant.
        if signature == "Fsorted":
            continue

        before = len(candidates)

        if before == 0:
            break

        candidates = filter_candidates(
            candidates,
            target,
            m,
            signature,
        )

        used.append(signature)

        survives = target in candidates

        print(
            f"  {signature:12s}: "
            f"{before:8d} -> {len(candidates):8d}"
            f"   target={survives}"
        )

        # --------------------------------------------------------------
        # EARLY STOP
        # --------------------------------------------------------------

        if len(candidates) <= 1:

            if len(candidates) == 1 and candidates[0] == target:
                print("  >>> TARGET ISOLATED")

                return {
                    "survivors": candidates,
                    "order": used,
                    "isolated": True,
                    "target_survives": True,
                    "time": time.perf_counter() - t0,
                }

            break

        # --------------------------------------------------------------
        # If target disappears, stop.
        # --------------------------------------------------------------

        if not survives:
            print("  >>> TARGET LOST")
            break

    isolated = (
        len(candidates) == 1
        and candidates[0] == target
    )

    return {
        "survivors": candidates,
        "order": used,
        "isolated": isolated,
        "target_survives": target in candidates,
        "time": time.perf_counter() - t0,
    }


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    total_start = time.perf_counter()

    banner(
        "KAPPA EXPERIMENT 32\n"
        "FAST PRUNED MINIMAL SIGNATURE SEPARATION\n"
        "NO CSV OUTPUT"
    )

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {POOL_SIZE:,}")
    print(f"prime interval   = [{LOW:,}, {HIGH:,}]")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # ------------------------------------------------------------------
    # MODULUS INVENTORY
    # ------------------------------------------------------------------

    section("MODULUS INVENTORY")

    inventory = build_modulus_inventory()

    for r, m, units in inventory:
        print(
            f"r={r:3d} "
            f"m={m:7d} "
            f"units={units:7d}"
        )

    # ------------------------------------------------------------------
    # PRIME POOL
    # ------------------------------------------------------------------

    section("PRIME POOL")

    t0 = time.perf_counter()

    print("Generating prime pool...")

    primes = generate_prime_pool()

    print(f"generated primes = {len(primes):,}")
    print(f"generation time = {elapsed(t0)}")

    # ------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------

    section("TARGETS")

    targets = generate_targets(primes)

    for i, (p, q) in enumerate(targets, 1):
        print(
            f"target {i:2d}: "
            f"p={p} "
            f"q={q} "
            f"n={p*q}"
        )

    # ------------------------------------------------------------------
    # EXPERIMENT
    # ------------------------------------------------------------------

    banner("FAST PRUNED MINIMAL SIGNATURE SEPARATION")

    summary = {
        prefix: {
            "isolated": 0,
            "not_isolated": 0,
            "total": 0,
        }
        for prefix in PREFIXES
    }

    for target_number, target in enumerate(targets, 1):

        section(
            f"TARGET {target_number:2d} "
            f"{target}"
        )

        for prefix in PREFIXES:

            # r = prefix in the existing experiment structure.
            #
            # modulus:
            #   m = 3r + 1
            #
            r = prefix
            m = 3 * r + 1

            result = separate_target(
                primes,
                target,
                prefix,
                m,
            )

            summary[prefix]["total"] += 1

            if result["isolated"]:
                summary[prefix]["isolated"] += 1
            else:
                summary[prefix]["not_isolated"] += 1

            survivors = result["survivors"]

            if len(survivors) <= 10:
                print(
                    "  survivors:",
                    ", ".join(map(str, survivors))
                    if survivors else "NONE"
                )
            else:
                print(
                    f"  survivors: {len(survivors):,} "
                    f"(not printing full list)"
                )

            print(
                f"  signature order = "
                f"{' -> '.join(result['order'])}"
            )

            print(
                f"  time = {result['time']:.3f}s"
            )

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    banner("GLOBAL SUMMARY")

    print(
        "prefix       isolated       not isolated       total"
    )
    print("-" * 62)

    for prefix in PREFIXES:

        x = summary[prefix]

        print(
            f"{prefix:6d}"
            f"{x['isolated']:16d}"
            f"{x['not_isolated']:19d}"
            f"{x['total']:12d}"
        )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    banner("EXPERIMENT 32 COMPLETE")

    print(
        "This experiment uses candidate-first filtering."
    )

    print(
        "Fsorted is omitted because it contains no information "
        "beyond Fpair once Fpair has already been applied."
    )

    print(
        "No pairwise intersection matrices are constructed."
    )

    print(
        "Large survivor lists are not printed."
    )

    print(
        "zero survivors = no surviving candidate in the finite pool"
    )

    print(
        "one survivor equal to target = target isolated in the pool"
    )

    print(
        "This does NOT establish global uniqueness."
    )

    total_time = time.perf_counter() - total_start

    print()
    print(
        f"total runtime = {total_time:.2f}s "
        f"({total_time / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()
