#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 29
MINIMAL SIGNATURE SEPARATION / INTERSECTION LADDER
NO CSV OUTPUT
==============================================================================

Purpose
-------
Experiment 28 attempted to construct large hash indexes independently for
every signature. That became unnecessarily expensive.

Experiment 29 instead asks:

    Starting from a target pair (p,q), how many candidate pairs survive
    progressively stronger intersections of signatures?

For every prefix:

    1. Find collisions for individual signatures.
    2. Measure pairwise intersections.
    3. Find the strongest useful 3-way intersections.
    4. Continue until the intersection becomes empty.
    5. Record the minimum number of signatures needed to separate the target
       from the searchable prime pool.

Important
---------
This is a finite pool experiment.

"NONE FOUND" means no collision was found in the current 5,000-prime
population. It does NOT prove global uniqueness.

No CSV output.
Only terminal printing.

Uses:
    from sympy import isprime
"""

import random
import time
from collections import defaultdict
from itertools import combinations

from sympy import isprime
from sympy.ntheory import n_order


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

TARGETS = 12
PRIME_POOL_SIZE = 5000

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29,
    31, 37, 41, 43, 47,
]

# Maximum number of candidates retained for verbose examples.
DISPLAY_LIMIT = 8

# Maximum intersection size we explicitly print.
INTERSECTION_DISPLAY_LIMIT = 8

SIGNATURE_ORDER = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "order_pair",
]


# ============================================================================
# MODULUS FUNCTIONS
# ============================================================================

def modulus_for_r(r):
    """
    m = r^2 + r + 1
    """
    return r * r + r + 1


def unit_count(m):
    count = 0
    for x in range(1, m):
        if __import__("math").gcd(x, m) == 1:
            count += 1
    return count


def build_modulus_inventory():
    inventory = {}

    for r in R_VALUES:
        m = modulus_for_r(r)
        inventory[r] = {
            "m": m,
            "units": unit_count(m),
        }

    return inventory


# ============================================================================
# PRIME POOL
# ============================================================================

def generate_prime_pool(rng):
    """
    Generate exactly PRIME_POOL_SIZE unique primes in the requested interval.

    sympy.isprime is deliberately used here, as requested.
    """

    primes = set()

    while len(primes) < PRIME_POOL_SIZE:
        x = rng.randrange(PRIME_MIN, PRIME_MAX)

        if isprime(x):
            primes.add(x)

    return sorted(primes)


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(rng, prime_pool):
    """
    Select 12 distinct factor pairs from the prime pool.

    Targets are selected deterministically from the seeded pool.

    The target factors are added to the searchable population automatically,
    so the target itself is always represented in the candidate universe.
    """

    if len(prime_pool) < 2 * TARGETS:
        raise RuntimeError("Prime pool too small.")

    shuffled = list(prime_pool)
    rng.shuffle(shuffled)

    targets = []

    used = set()

    for x in shuffled:
        if x in used:
            continue

        for y in shuffled:
            if y == x or y in used:
                continue

            targets.append((x, y))
            used.add(x)
            used.add(y)
            break

        if len(targets) >= TARGETS:
            break

    return targets


# ============================================================================
# RESIDUE CACHE
# ============================================================================

def build_residue_cache(primes, inventory):
    """
    residue_cache[p][r] = p mod m_r
    """

    cache = {}

    for p in primes:
        row = {}

        for r in R_VALUES:
            m = inventory[r]["m"]
            row[r] = p % m

        cache[p] = row

    return cache


# ============================================================================
# PER-PRIME POWER SIGNATURES
# ============================================================================

def build_power_vectors(primes, inventory):
    """
    power_vectors[p][prefix] =

        (
            p^2 mod m_2,
            p^3 mod m_3,
            ...
        )

    where prefix k means the first k entries of R_VALUES.

    More precisely, for each r we store:

        p^r mod m_r

    This is the fundamental information used by Fpair, Fsorted, Fsum,
    Fprod and Fdiff.
    """

    result = {}

    for p in primes:
        prefix_vectors = {}

        row = []

        for r in R_VALUES:
            m = inventory[r]["m"]
            row.append(pow(p, r, m))

        for prefix in PREFIXES:
            prefix_len = prefix

            prefix_vectors[prefix] = tuple(
                row[:prefix_len]
            )

        result[p] = prefix_vectors

    return result


# ============================================================================
# ORDER CACHE
# ============================================================================

def build_order_vectors(primes, inventory):
    """
    Lazy multiplicative-order cache.

    We intentionally calculate orders only when needed.

    n_order(a,m) requires gcd(a,m)=1. For the prime pool used here, this is
    normally satisfied, but we explicitly handle non-units.
    """

    cache = {}

    def get(p, prefix):
        if p not in cache:
            cache[p] = {}

        if prefix in cache[p]:
            return cache[p][prefix]

        values = []

        for r in R_VALUES[:prefix]:
            m = inventory[r]["m"]
            a = p % m

            import math

            if math.gcd(a, m) != 1:
                values.append(0)
            else:
                values.append(int(n_order(a, m)))

        value = tuple(values)

        cache[p][prefix] = value
        return value

    return cache, get


# ============================================================================
# SIGNATURE CONSTRUCTION
# ============================================================================

def signature_from_vectors(
    p,
    q,
    prefix,
    power_vectors,
    order_getter,
    name,
):
    """
    Construct one signature for a candidate pair.

    Definitions
    -----------
    Fpair:
        ordered pair of power vectors

    Fsorted:
        lexicographically sorted power vectors

    Fsum:
        coordinate-wise sum modulo m

    Fprod:
        coordinate-wise product modulo m

    Fdiff:
        coordinate-wise absolute difference modulo m

    cube_sum:
        p^3 + q^3 modulo m_r for every r

    order_pair:
        ordered pair of multiplicative-order vectors
    """

    vp = power_vectors[p][prefix]
    vq = power_vectors[q][prefix]

    mods = [
        modulus_for_r(r)
        for r in R_VALUES[:prefix]
    ]

    if name == "Fpair":
        return (vp, vq)

    if name == "Fsorted":
        return tuple(sorted((vp, vq)))

    if name == "Fsum":
        return tuple(
            (a + b) % m
            for a, b, m in zip(vp, vq, mods)
        )

    if name == "Fprod":
        return tuple(
            (a * b) % m
            for a, b, m in zip(vp, vq, mods)
        )

    if name == "Fdiff":
        return tuple(
            abs(a - b) % m
            for a, b, m in zip(vp, vq, mods)
        )

    if name == "cube_sum":
        return tuple(
            (pow(p, 3, m) + pow(q, 3, m)) % m
            for m in mods
        )

    if name == "order_pair":
        op = order_getter(p, prefix)
        oq = order_getter(q, prefix)

        return (op, oq)

    raise ValueError(f"Unknown signature: {name}")


# ============================================================================
# FAST COMPLEMENT INDEXES
# ============================================================================

def build_vector_index(primes, power_vectors, prefix):
    """
    Index:

        power vector -> primes

    This is enough to solve Fsum/Fprod/Fdiff without constructing O(P^2)
    candidate pairs.
    """

    index = defaultdict(list)

    for p in primes:
        index[power_vectors[p][prefix]].append(p)

    return index


def find_sum_collisions(
    target_p,
    target_q,
    primes,
    power_vectors,
    prefix,
):
    """
    Fsum collision search.

    Given target vector T:

        x + y = T (mod m)

    For each x, calculate the required y-vector and look it up.
    """

    vp = power_vectors[target_p][prefix]
    vq = power_vectors[target_q][prefix]

    target_sum = tuple(
        (a + b) % modulus_for_r(r)
        for a, b, r in zip(vp, vq, R_VALUES[:prefix])
    )

    index = build_vector_index(
        primes,
        power_vectors,
        prefix,
    )

    result = set()

    for p in primes:
        vp_candidate = power_vectors[p][prefix]

        needed = tuple(
            (target_sum[i] - vp_candidate[i])
            % modulus_for_r(R_VALUES[i])
            for i in range(prefix)
        )

        for q in index.get(needed, ()):
            if p != q:
                result.add(
                    tuple(sorted((p, q)))
                )

    return result


def find_prod_collisions(
    target_p,
    target_q,
    primes,
    power_vectors,
    prefix,
):
    """
    Fprod collision search.

    Because every target prime is coprime to m_r, calculate:

        q^r = target_prod_r * inverse(p^r) mod m_r

    If p^r is not invertible, skip that candidate.
    """

    import math

    vp = power_vectors[target_p][prefix]
    vq = power_vectors[target_q][prefix]

    target_prod = tuple(
        (a * b) % modulus_for_r(r)
        for a, b, r in zip(vp, vq, R_VALUES[:prefix])
    )

    index = build_vector_index(
        primes,
        power_vectors,
        prefix,
    )

    result = set()

    for p in primes:
        vp_candidate = power_vectors[p][prefix]

        needed = []

        valid = True

        for i, r in enumerate(R_VALUES[:prefix]):
            m = modulus_for_r(r)
            a = vp_candidate[i]

            if math.gcd(a, m) != 1:
                valid = False
                break

            inv = pow(a, -1, m)

            needed.append(
                (target_prod[i] * inv) % m
            )

        if not valid:
            continue

        needed = tuple(needed)

        for q in index.get(needed, ()):
            if p != q:
                result.add(
                    tuple(sorted((p, q)))
                )

    return result


def find_diff_collisions(
    target_p,
    target_q,
    primes,
    power_vectors,
    prefix,
):
    """
    Fdiff collision search.

    The absolute difference gives two possible coordinate vectors:

        q = p + d
        q = p - d

    modulo each m_r.

    We therefore create both possible required vectors and look them up.
    """

    vp = power_vectors[target_p][prefix]
    vq = power_vectors[target_q][prefix]

    diff = tuple(
        abs(a - b) % modulus_for_r(r)
        for a, b, r in zip(vp, vq, R_VALUES[:prefix])
    )

    index = build_vector_index(
        primes,
        power_vectors,
        prefix,
    )

    result = set()

    for p in primes:
        vp_candidate = power_vectors[p][prefix]

        needed_plus = tuple(
            (vp_candidate[i] + diff[i])
            % modulus_for_r(R_VALUES[i])
            for i in range(prefix)
        )

        needed_minus = tuple(
            (vp_candidate[i] - diff[i])
            % modulus_for_r(R_VALUES[i])
            for i in range(prefix)
        )

        for needed in (needed_plus, needed_minus):
            for q in index.get(needed, ()):
                if p != q:
                    result.add(
                        tuple(sorted((p, q)))
                    )

    return result


# ============================================================================
# GENERIC SIGNATURE INDEX
# ============================================================================

def build_pair_signature_index(
    primes,
    prefix,
    power_vectors,
    order_getter,
    signature_name,
):
    """
    Build a signature index by enumerating unordered prime pairs.

    This function is intentionally NOT used for Fsum/Fprod/Fdiff.

    It is reserved for signatures where a complement lookup is not available.

    For the 5,000-prime pool, this is approximately 12.5 million pairs, so
    we use it only for signatures where absolutely necessary.
    """

    index = defaultdict(list)

    n = len(primes)

    for i in range(n):
        p = primes[i]

        for j in range(i + 1, n):
            q = primes[j]

            sig = signature_from_vectors(
                p,
                q,
                prefix,
                power_vectors,
                order_getter,
                signature_name,
            )

            index[sig].append((p, q))

    return index


# ============================================================================
# TARGET SIGNATURE
# ============================================================================

def target_signature(
    p,
    q,
    prefix,
    name,
    power_vectors,
    order_getter,
):
    return signature_from_vectors(
        p,
        q,
        prefix,
        power_vectors,
        order_getter,
        name,
    )


# ============================================================================
# DIRECT COLLISION SEARCH
# ============================================================================

def direct_signature_collisions(
    target_p,
    target_q,
    primes,
    prefix,
    name,
    power_vectors,
    order_getter,
):
    """
    Search collisions for a single signature.

    For Fsum/Fprod/Fdiff use the O(P)-style complement search.

    For Fpair/Fsorted/cube_sum/order_pair we use a two-stage strategy:

        1. calculate target signature;
        2. scan candidate p;
        3. use structural pruning before calculating the expensive full pair
           signature.

    Fpair/Fsorted can be solved directly using the individual vectors.

    cube_sum is also reducible to an individual vector lookup because
    p^3 mod m_r is simply another per-prime vector.

    order_pair is handled by order-vector complement indexing.
    """

    target_sig = target_signature(
        target_p,
        target_q,
        prefix,
        name,
        power_vectors,
        order_getter,
    )

    result = set()

    # ------------------------------------------------------------------
    # Fpair
    # ------------------------------------------------------------------

    if name == "Fpair":
        vp_target = power_vectors[target_p][prefix]
        vq_target = power_vectors[target_q][prefix]

        index = build_vector_index(
            primes,
            power_vectors,
            prefix,
        )

        for p in primes:
            if p == target_p:
                pass

            vp = power_vectors[p][prefix]

            for q in index.get(vq_target, ()):
                if p == q:
                    continue

                if vp == vp_target:
                    result.add(tuple(sorted((p, q))))

        return result

    # ------------------------------------------------------------------
    # Fsorted
    # ------------------------------------------------------------------

    if name == "Fsorted":
        vp_target = power_vectors[target_p][prefix]
        vq_target = power_vectors[target_q][prefix]

        index = build_vector_index(
            primes,
            power_vectors,
            prefix,
        )

        for first_target, second_target in (
            (vp_target, vq_target),
            (vq_target, vp_target),
        ):
            for p in index.get(first_target, ()):
                for q in index.get(second_target, ()):
                    if p != q:
                        result.add(tuple(sorted((p, q))))

        return result

    # ------------------------------------------------------------------
    # Fsum
    # ------------------------------------------------------------------

    if name == "Fsum":
        return find_sum_collisions(
            target_p,
            target_q,
            primes,
            power_vectors,
            prefix,
        )

    # ------------------------------------------------------------------
    # Fprod
    # ------------------------------------------------------------------

    if name == "Fprod":
        return find_prod_collisions(
            target_p,
            target_q,
            primes,
            power_vectors,
            prefix,
        )

    # ------------------------------------------------------------------
    # Fdiff
    # ------------------------------------------------------------------

    if name == "Fdiff":
        return find_diff_collisions(
            target_p,
            target_q,
            primes,
            power_vectors,
            prefix,
        )

    # ------------------------------------------------------------------
    # cube_sum
    # ------------------------------------------------------------------

    if name == "cube_sum":
        mods = [
            modulus_for_r(r)
            for r in R_VALUES[:prefix]
        ]

        target_cube = tuple(
            (
                pow(target_p, 3, m)
                + pow(target_q, 3, m)
            ) % m
            for m in mods
        )

        cube_vectors = {}

        for p in primes:
            cube_vectors[p] = tuple(
                pow(p, 3, m) for m in mods
            )

        index = defaultdict(list)

        for p in primes:
            index[cube_vectors[p]].append(p)

        for p in primes:
            vp = cube_vectors[p]

            needed = tuple(
                (target_cube[i] - vp[i]) % mods[i]
                for i in range(prefix)
            )

            for q in index.get(needed, ()):
                if p != q:
                    result.add(tuple(sorted((p, q))))

        return result

    # ------------------------------------------------------------------
    # order_pair
    # ------------------------------------------------------------------

    if name == "order_pair":

        op_target = order_getter(target_p, prefix)
        oq_target = order_getter(target_q, prefix)

        order_index = defaultdict(list)

        for p in primes:
            order_index[
                order_getter(p, prefix)
            ].append(p)

        for a, b in (
            (op_target, oq_target),
            (oq_target, op_target),
        ):
            for p in order_index.get(a, ()):
                for q in order_index.get(b, ()):
                    if p != q:
                        result.add(tuple(sorted((p, q))))

        return result

    raise ValueError(name)


# ============================================================================
# INTERSECTION HELPERS
# ============================================================================

def intersection_all(sets):
    if not sets:
        return set()

    result = set(sets[0])

    for s in sets[1:]:
        result.intersection_update(s)

        if not result:
            break

    return result


def format_pair(pair):
    return f"({pair[0]}, {pair[1]})"


def format_candidates(candidates):
    if not candidates:
        return "NONE"

    ordered = sorted(candidates)

    shown = ordered[:INTERSECTION_DISPLAY_LIMIT]

    text = ", ".join(
        format_pair(x)
        for x in shown
    )

    if len(ordered) > len(shown):
        text += f", ... +{len(ordered) - len(shown)} more"

    return text


# ============================================================================
# ANALYSIS
# ============================================================================

def analyze_intersections(
    target_p,
    target_q,
    prefix,
    collision_sets,
):
    """
    Calculate:

        self survival
        pair intersections
        triple intersections
        ...

    We use a greedy minimal-separation ladder:

        Start with the signature having the fewest collisions.

        Add the signature that produces the smallest remaining intersection.

    This is NOT claimed to be the globally optimal set cover.

    It is a deterministic greedy diagnostic.
    """

    target_pair = tuple(sorted((target_p, target_q)))

    names = [
        name
        for name in SIGNATURE_ORDER
        if name in collision_sets
    ]

    print()
    print("INTERSECTION SUMMARY")
    print("-" * 76)

    print(
        f"target = {format_pair(target_pair)}"
    )

    print()

    print(
        f"{'signature':<16}"
        f"{'self collisions':>17}"
        f"{'contains target':>17}"
    )

    print("-" * 76)

    for name in names:
        s = collision_sets[name]

        print(
            f"{name:<16}"
            f"{len(s):>17}"
            f"{str(target_pair in s):>17}"
        )

    # --------------------------------------------------------------
    # Pairwise intersection matrix
    # --------------------------------------------------------------

    print()
    print("PAIRWISE INTERSECTION MATRIX")
    print("-" * 76)

    print(
        "Each cell = number of candidate pairs satisfying BOTH signatures."
    )

    print()

    header = "source".ljust(16)

    for name in names:
        header += f"{name[:12]:>13}"

    print(header)
    print("-" * len(header))

    for a in names:
        row = a.ljust(16)

        for b in names:
            value = len(
                collision_sets[a]
                & collision_sets[b]
            )

            row += f"{value:>13}"

        print(row)

    # --------------------------------------------------------------
    # Greedy minimal-separation ladder
    # --------------------------------------------------------------

    print()
    print("GREEDY MINIMAL-SEPARATION LADDER")
    print("-" * 76)

    remaining_names = set(names)

    current = None
    selected = []

    step = 0

    while remaining_names:

        step += 1

        best_name = None
        best_set = None

        for name in sorted(remaining_names):
            if current is None:
                candidate = set(
                    collision_sets[name]
                )
            else:
                candidate = (
                    current
                    & collision_sets[name]
                )

            if best_set is None:
                best_name = name
                best_set = candidate
                continue

            # Prefer fewer candidates.
            # Break ties lexicographically.
            if len(candidate) < len(best_set):
                best_name = name
                best_set = candidate
            elif (
                len(candidate) == len(best_set)
                and name < best_name
            ):
                best_name = name
                best_set = candidate

        selected.append(best_name)
        remaining_names.remove(best_name)

        current = best_set

        contains_target = target_pair in current

        print(
            f"step {step}: "
            f"+ {best_name:<12} "
            f"survivors = {len(current):<8} "
            f"target_survives = {contains_target}"
        )

        if len(current) <= INTERSECTION_DISPLAY_LIMIT:
            print(
                "           "
                + format_candidates(current)
            )

        if not current:
            print()
            print(
                "           ZERO-CANDIDATE INTERSECTION REACHED."
            )
            print(
                "           Remaining signatures are not needed for"
                " this finite pool."
            )
            break

    # --------------------------------------------------------------
    # First separating prefix within greedy ladder
    # --------------------------------------------------------------

    print()
    print("MINIMAL-SEPARATION RESULT")
    print("-" * 76)

    if current is not None and len(current) == 0:
        print(
            "A zero-candidate intersection was reached."
        )
        print(
            "This means the selected signature combination separates"
            " the target from every other searchable pair in this pool."
        )
    else:
        print(
            "No zero-candidate intersection was reached."
        )
        print(
            "The finite pool still contains at least one surviving"
            " candidate under the full tested signature family."
        )

    print()
    print(
        "Greedy signature order:"
    )
    print(
        "  " + " -> ".join(selected)
    )

    return selected


# ============================================================================
# MAIN
# ============================================================================

def main():

    total_start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 29")
    print("MINIMAL SIGNATURE SEPARATION / INTERSECTION LADDER")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {PRIME_POOL_SIZE:,}")
    print(
        f"prime interval   = "
        f"[{PRIME_MIN:,}, {PRIME_MAX:,}]"
    )
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # ------------------------------------------------------------------
    # MODULUS INVENTORY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    inventory = build_modulus_inventory()

    for r in R_VALUES:
        m = inventory[r]["m"]
        units = inventory[r]["units"]

        print(
            f"r={r:>3} "
            f"m={m:>7} "
            f"units={units:>7}"
        )

    # ------------------------------------------------------------------
    # PRIME POOL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIME POOL GENERATION")
    print("=" * 78)

    t0 = time.perf_counter()

    print("Generating unique prime pool...")

    prime_pool = generate_prime_pool(rng)

    print(
        f"generated primes = {len(prime_pool):,}"
    )

    print(
        f"pool generation time = "
        f"{time.perf_counter() - t0:.2f}s"
    )

    # ------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TARGET GENERATION")
    print("=" * 78)

    targets = generate_targets(
        rng,
        prime_pool,
    )

    target_primes = {
        p
        for pair in targets
        for p in pair
    }

    searchable_primes = sorted(
        set(prime_pool)
        | target_primes
    )

    print("Generating fixed target set...")

    for i, (p, q) in enumerate(targets, 1):
        n = p * q

        print(
            f"target {i:>2}: "
            f"p={p} "
            f"q={q} "
            f"n={n}"
        )

    if len(searchable_primes) != len(prime_pool):
        print()
        print(
            "NOTE: fixed target factors were not all in the original"
            " generated pool."
        )
        print(
            "They were added to the searchable population."
        )

    # ------------------------------------------------------------------
    # RESIDUE CACHE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PRIME RESIDUE CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    residue_cache = build_residue_cache(
        searchable_primes,
        inventory,
    )

    print("residue cache complete")

    print(
        f"cache time = "
        f"{time.perf_counter() - t0:.2f}s"
    )

    # ------------------------------------------------------------------
    # POWER VECTORS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. POWER VECTOR CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    power_vectors = build_power_vectors(
        searchable_primes,
        inventory,
    )

    print("power-vector cache complete")

    print(
        f"cache time = "
        f"{time.perf_counter() - t0:.2f}s"
    )

    # ------------------------------------------------------------------
    # ORDER CACHE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. MULTIPLICATIVE ORDER CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    order_cache, order_getter = build_order_vectors(
        searchable_primes,
        inventory,
    )

    print(
        "using lazy multiplicative-order cache"
    )

    print(
        f"cache setup time = "
        f"{time.perf_counter() - t0:.2f}s"
    )

    # ------------------------------------------------------------------
    # EXPERIMENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. SIGNATURE INTERSECTION EXPERIMENT")
    print("=" * 78)

    all_results = []

    for target_index, (target_p, target_q) in enumerate(
        targets,
        1,
    ):

        print()
        print("=" * 78)
        print(f"TARGET {target_index}")
        print("=" * 78)

        n = target_p * target_q
        s = target_p + target_q
        delta = (
            target_p * target_p
            + target_q * target_q
        )

        print(f"p       = {target_p}")
        print(f"q       = {target_q}")
        print(f"n       = {n}")
        print(f"s       = {s}")
        print(f"Delta   = {delta}")
        print(f"n bits  = {n.bit_length()}")

        target_result = {}

        for prefix in PREFIXES:

            print()
            print("-" * 78)
            print(
                f"PREFIX {prefix} "
                f"(m={modulus_for_r(R_VALUES[prefix - 1])})"
            )
            print("-" * 78)

            collision_sets = {}

            for name in SIGNATURE_ORDER:

                t0 = time.perf_counter()

                collisions = direct_signature_collisions(
                    target_p,
                    target_q,
                    searchable_primes,
                    prefix,
                    name,
                    power_vectors,
                    order_getter,
                )

                elapsed = (
                    time.perf_counter()
                    - t0
                )

                collision_sets[name] = collisions

                target_pair = tuple(
                    sorted((target_p, target_q))
                )

                self_collision_count = len(
                    collisions
                )

                wrong_collisions = (
                    collisions
                    - {target_pair}
                )

                print(
                    f"{name:<16}"
                    f"self={self_collision_count:<7}"
                    f"wrong={len(wrong_collisions):<7}"
                    f"[{elapsed:.4f}s]"
                )

                if (
                    wrong_collisions
                    and len(wrong_collisions)
                    <= DISPLAY_LIMIT
                ):
                    print(
                        "    wrong candidates: "
                        + format_candidates(
                            wrong_collisions
                        )
                    )

            # ----------------------------------------------------------
            # INTERSECTION ANALYSIS
            # ----------------------------------------------------------

            selected = analyze_intersections(
                target_p,
                target_q,
                prefix,
                collision_sets,
            )

            target_result[prefix] = {
                "collisions": collision_sets,
                "greedy_order": selected,
            }

        all_results.append(
            target_result
        )

    # ------------------------------------------------------------------
    # GLOBAL SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. GLOBAL SUMMARY")
    print("=" * 78)

    print()
    print(
        "For each prefix, the following reports the number of targets"
        " for which each signature had at least one WRONG collision."
    )

    print()

    print(
        f"{'signature':<16}"
        + "".join(
            f"{prefix:>10}"
            for prefix in PREFIXES
        )
    )

    print("-" * 76)

    for name in SIGNATURE_ORDER:

        row = f"{name:<16}"

        for prefix in PREFIXES:

            count = 0

            for target_result in all_results:

                collisions = target_result[prefix][
                    "collisions"
                ][name]

                # Remove the genuine target.
                # Any remaining pair is a wrong collision.
                if len(collisions) > 1:
                    count += 1

            row += f"{count:>10}"

        print(row)

    # ------------------------------------------------------------------
    # ZERO SURVIVAL SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. ZERO-SURVIVAL INTERSECTION SUMMARY")
    print("=" * 78)

    for prefix in PREFIXES:

        print()
        print(f"PREFIX {prefix}")
        print("-" * 50)

        for target_index, target_result in enumerate(
            all_results,
            1,
        ):

            selected = target_result[prefix][
                "greedy_order"
            ]

            collisions = target_result[prefix][
                "collisions"
            ]

            current = None
            zero_step = None

            for step, name in enumerate(
                selected,
                1,
            ):

                if current is None:
                    current = set(
                        collisions[name]
                    )
                else:
                    current &= collisions[name]

                if (
                    not current
                    and zero_step is None
                ):
                    zero_step = step

            if zero_step is None:
                final_count = (
                    len(current)
                    if current is not None
                    else 0
                )

                print(
                    f"target {target_index:>2}: "
                    f"no zero intersection "
                    f"(final={final_count})"
                )

            else:
                print(
                    f"target {target_index:>2}: "
                    f"zero intersection at "
                    f"step {zero_step}"
                )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    elapsed = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 29 COMPLETE")
    print("=" * 78)

    print()
    print(
        "Experiment 29 measures the minimum practical signature"
        " intersection needed to eliminate the finite-pool"
        " collision set."
    )

    print()
    print(
        "Important:"
    )
    print(
        "  zero survivors = uniqueness within this searchable pool"
    )
    print(
        "  NONE FOUND     = no collision observed in this pool"
    )
    print(
        "  neither result proves global uniqueness"
    )

    print()
    print(
        f"total runtime = {elapsed:.2f}s "
        f"({elapsed / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

