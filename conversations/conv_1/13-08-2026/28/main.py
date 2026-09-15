#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 28
SIGNATURE INTERSECTION / MINIMAL SEPARATION
NO CSV OUTPUT
==============================================================================

Question:

    Experiment 27 measured whether a collision for signature X
    also fools signature Y.

    Experiment 28 asks:

        Which COMBINATIONS of signatures are actually independent?

For every target, prefix and source signature:

    1. Find a constructive collision for the source signature.
    2. Test that candidate against every signature.
    3. Record the complete survival mask.
    4. Aggregate survival over all targets.
    5. Evaluate every non-empty signature subset.
    6. Find the smallest subsets that eliminate collisions.

Important:

    NONE FOUND does NOT imply uniqueness.

The experiment only measures collisions inside the finite prime pool.

No CSV files are written.
All results are printed.
==============================================================================

random seed      = 20260814
targets          = 12
prime pool       = 5,000
prime interval   = [2,000,000, 4,200,000]
prefixes         = [3, 4, 5, 6, 7]

R values:
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

Signatures:
    Fpair
    Fsorted
    Fsum
    Fprod
    Fdiff
    cube_sum
    order_pair
    order_sorted
"""

import random
import time
from itertools import combinations

import sympy as sp


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

TARGETS = 12
POOL_SIZE = 5_000

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

SIGNATURE_NAMES = [
    "Fpair",
    "Fsorted",
    "Fsum",
    "Fprod",
    "Fdiff",
    "cube_sum",
    "order_pair",
    "order_sorted",
]

random.seed(SEED)


# ============================================================================
# MODULUS CONSTRUCTION
# ============================================================================

def modulus_for_r(r):
    """
    m = r^2 + r + 1

    This reproduces:

        r=2  -> 7
        r=3  -> 13
        r=5  -> 31
        r=7  -> 57
        ...
    """
    return r * r + r + 1


def unit_count(m):
    return sum(1 for x in range(1, m) if sp.gcd(x, m) == 1)


# ============================================================================
# PRIME GENERATION
# ============================================================================

def generate_prime_pool():
    print("Generating unique prime pool...")

    pool = set()

    while len(pool) < POOL_SIZE:
        x = random.randint(PRIME_LOW, PRIME_HIGH)

        if sp.isprime(x):
            pool.add(x)

    return sorted(pool)


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(pool):
    """
    Construct 12 semiprime targets from distinct primes.

    The targets are deliberately chosen from the same searchable
    population so every target factor can be evaluated by the
    signature machinery.

    We keep p != q.
    """

    print("Generating fixed target set...")

    targets = []

    # Shuffle a copy so the target pairs are reproducible but not simply
    # the first twelve pool entries.
    candidates = pool[:]
    random.shuffle(candidates)

    i = 0

    while len(targets) < TARGETS:
        p = candidates[i]
        q = candidates[i + 1]
        i += 2

        if p == q:
            continue

        targets.append((p, q))

    return targets


# ============================================================================
# PRIME RESIDUE CACHE
# ============================================================================

def build_residue_cache(primes, moduli):
    cache = {}

    for p in primes:
        cache[p] = {
            m: p % m
            for m in moduli
        }

    return cache


# ============================================================================
# MULTIPLICATIVE ORDER
# ============================================================================

def multiplicative_order(a, m):
    """
    Return ord_m(a).

    If gcd(a,m) != 1, the multiplicative order is undefined.
    We represent it by 0.
    """

    if sp.gcd(a, m) != 1:
        return 0

    return int(sp.n_order(a, m))


def build_lazy_order_cache(primes, moduli, residue_cache):
    """
    Lazy cache.

    We only calculate an order when a signature actually requests it.
    This is considerably cheaper than calculating every order for every
    prime and every modulus up front.
    """

    cache = {}

    def get(p, m):
        if p not in cache:
            cache[p] = {}

        if m not in cache[p]:
            a = residue_cache[p][m]
            cache[p][m] = multiplicative_order(a, m)

        return cache[p][m]

    return get


# ============================================================================
# SIGNATURE COMPONENTS
# ============================================================================

def pair_residue(p, q, m, residue_cache):
    return (
        residue_cache[p][m],
        residue_cache[q][m],
    )


def sorted_residue(p, q, m, residue_cache):
    return tuple(sorted((
        residue_cache[p][m],
        residue_cache[q][m],
    )))


def sum_residue(p, q, m, residue_cache):
    return (
        residue_cache[p][m] +
        residue_cache[q][m]
    ) % m


def product_residue(p, q, m, residue_cache):
    return (
        residue_cache[p][m] *
        residue_cache[q][m]
    ) % m


def difference_residue(p, q, m, residue_cache):
    return (
        residue_cache[p][m] -
        residue_cache[q][m]
    ) % m


def cube_sum_residue(p, q, m, residue_cache):
    a = residue_cache[p][m]
    b = residue_cache[q][m]

    return (a ** 3 + b ** 3) % m


def order_pair(p, q, m, residue_cache, order_get):
    return (
        order_get(p, m),
        order_get(q, m),
    )


def order_sorted(p, q, m, residue_cache, order_get):
    return tuple(sorted((
        order_get(p, m),
        order_get(q, m),
    )))


# ============================================================================
# SIGNATURE
# ============================================================================

def signature(
    name,
    p,
    q,
    prefix,
    residue_cache,
    order_get,
):
    """
    A signature is the tuple of values over the first `prefix`
    R-values.

    prefix=3 means:

        R = [2,3,5]

    prefix=4 means:

        R = [2,3,5,7]

    etc.
    """

    values = []

    for r in R_VALUES[:prefix]:
        m = modulus_for_r(r)

        if name == "Fpair":
            v = pair_residue(p, q, m, residue_cache)

        elif name == "Fsorted":
            v = sorted_residue(p, q, m, residue_cache)

        elif name == "Fsum":
            v = sum_residue(p, q, m, residue_cache)

        elif name == "Fprod":
            v = product_residue(p, q, m, residue_cache)

        elif name == "Fdiff":
            v = difference_residue(p, q, m, residue_cache)

        elif name == "cube_sum":
            v = cube_sum_residue(p, q, m, residue_cache)

        elif name == "order_pair":
            v = order_pair(
                p,
                q,
                m,
                residue_cache,
                order_get,
            )

        elif name == "order_sorted":
            v = order_sorted(
                p,
                q,
                m,
                residue_cache,
                order_get,
            )

        else:
            raise ValueError(f"Unknown signature: {name}")

        values.append(v)

    return tuple(values)


# ============================================================================
# HASH INDEX
# ============================================================================

def build_index(pool, prefix, name, residue_cache, order_get):
    """
    Maps:

        signature -> first prime pair

    We only need one candidate because the experiment asks whether
    a constructive collision exists.

    The index stores ordered pairs p < q.
    """

    index = {}

    for i, p in enumerate(pool):
        for q in pool[i + 1:]:
            sig = signature(
                name,
                p,
                q,
                prefix,
                residue_cache,
                order_get,
            )

            if sig not in index:
                index[sig] = (p, q)

    return index


# ============================================================================
# TARGET SIGNATURE
# ============================================================================

def target_signature_map(
    target,
    prefix,
    residue_cache,
    order_get,
):
    p, q = target

    result = {}

    for name in SIGNATURE_NAMES:
        result[name] = signature(
            name,
            p,
            q,
            prefix,
            residue_cache,
            order_get,
        )

    return result


# ============================================================================
# FIND SOURCE COLLISION
# ============================================================================

def source_collision(
    target,
    source_name,
    prefix,
    indexes,
    residue_cache,
    order_get,
):
    """
    Return the first wrong candidate matching source_name.

    If the target itself is present in the index, skip it and continue
    searching for a genuine collision.

    Because the index stores only one pair per signature, we additionally
    build a collision search when that first pair is the target.
    """

    target_sig = signature(
        source_name,
        target[0],
        target[1],
        prefix,
        residue_cache,
        order_get,
    )

    candidate = indexes[source_name].get(target_sig)

    if candidate is not None:
        if set(candidate) != set(target):
            return candidate

    # The first indexed pair can be the target itself.
    # Search explicitly for another pair with the same signature.
    found = None

    for i, p in enumerate(PRIME_POOL):
        for q in PRIME_POOL[i + 1:]:
            if {p, q} == set(target):
                continue

            sig = signature(
                source_name,
                p,
                q,
                prefix,
                residue_cache,
                order_get,
            )

            if sig == target_sig:
                found = (p, q)
                break

        if found is not None:
            break

    return found


# ============================================================================
# COMBINATION SURVIVAL
# ============================================================================

def candidate_matches(
    candidate,
    target_sigs,
    subset,
    prefix,
    residue_cache,
    order_get,
):
    """
    Does candidate satisfy every signature in subset?
    """

    p, q = candidate

    for name in subset:
        candidate_sig = signature(
            name,
            p,
            q,
            prefix,
            residue_cache,
            order_get,
        )

        if candidate_sig != target_sigs[name]:
            return False

    return True


# ============================================================================
# BITMASK UTILITIES
# ============================================================================

def signature_mask(candidate, target_sigs, prefix, residue_cache, order_get):
    mask = 0

    for i, name in enumerate(SIGNATURE_NAMES):
        candidate_sig = signature(
            name,
            candidate[0],
            candidate[1],
            prefix,
            residue_cache,
            order_get,
        )

        if candidate_sig == target_sigs[name]:
            mask |= (1 << i)

    return mask


def mask_names(mask):
    result = []

    for i, name in enumerate(SIGNATURE_NAMES):
        if mask & (1 << i):
            result.append(name)

    return result


# ============================================================================
# EXPERIMENT
# ============================================================================

def run_experiment(
    pool,
    targets,
    residue_cache,
    order_get,
):
    global PRIME_POOL
    PRIME_POOL = pool

    # ------------------------------------------------------------------------
    # Build indexes
    # ------------------------------------------------------------------------

    indexes_by_prefix = {}

    print()
    print("=" * 78)
    print("HASH INDEX CONSTRUCTION")
    print("=" * 78)

    for prefix in PREFIXES:
        indexes_by_prefix[prefix] = {}

        for name in SIGNATURE_NAMES:
            t0 = time.perf_counter()

            idx = build_index(
                pool,
                prefix,
                name,
                residue_cache,
                order_get,
            )

            indexes_by_prefix[prefix][name] = idx

            dt = time.perf_counter() - t0

            print(
                f"prefix {prefix}: "
                f"{name:<14} "
                f"{len(idx):>7} classes "
                f"[{dt:.3f}s]"
            )

    # ------------------------------------------------------------------------
    # Main experiment
    # ------------------------------------------------------------------------

    all_results = {}

    for prefix in PREFIXES:
        print()
        print("=" * 78)
        print(f"PREFIX {prefix}")
        print("=" * 78)

        prefix_results = []

        for target_number, target in enumerate(targets, 1):
            p, q = target

            print()
            print(f"TARGET {target_number}")
            print("-" * 78)
            print(f"target = ({p}, {q})")
            print(f"n      = {p * q}")
            print(f"s      = {p + q}")
            print(f"Delta  = {(p - q) ** 2}")
            print(f"bits   = {(p * q).bit_length()}")

            target_sigs = target_signature_map(
                target,
                prefix,
                residue_cache,
                order_get,
            )

            target_result = {
                "target": target,
                "sources": {},
            }

            # ---------------------------------------------------------------
            # Find source collisions
            # ---------------------------------------------------------------

            for source_name in SIGNATURE_NAMES:

                t0 = time.perf_counter()

                candidate = source_collision(
                    target,
                    source_name,
                    prefix,
                    indexes_by_prefix[prefix],
                    residue_cache,
                    order_get,
                )

                dt = time.perf_counter() - t0

                if candidate is None:
                    print(
                        f"{source_name:<14} "
                        f"NONE FOUND [{dt:.3f}s]"
                    )

                    target_result["sources"][source_name] = None
                    continue

                mask = signature_mask(
                    candidate,
                    target_sigs,
                    prefix,
                    residue_cache,
                    order_get,
                )

                matching = mask_names(mask)

                print(
                    f"{source_name:<14} "
                    f"candidate=({candidate[0]}, {candidate[1]}) "
                    f"[{dt:.3f}s]"
                )

                print(
                    f"{'':14} "
                    f"survives {len(matching)}/{len(SIGNATURE_NAMES)}:"
                )

                print(
                    " " * 14 +
                    ", ".join(matching)
                )

                target_result["sources"][source_name] = {
                    "candidate": candidate,
                    "mask": mask,
                    "matching": matching,
                }

            prefix_results.append(target_result)

        all_results[prefix] = prefix_results

    return all_results


# ============================================================================
# AGGREGATION
# ============================================================================

def analyze_intersections(all_results):
    print()
    print("=" * 78)
    print("INTERSECTION ANALYSIS")
    print("=" * 78)

    # ------------------------------------------------------------------------
    # Per-prefix source survival
    # ------------------------------------------------------------------------

    for prefix in PREFIXES:
        results = all_results[prefix]

        print()
        print(f"PREFIX {prefix}")
        print("-" * 78)

        print(
            f"{'source':<15}"
            f"{'collisions':>12}"
            f"{'cross':>12}"
            f"{'all':>12}"
            f"{'mean survivors':>18}"
        )

        for source in SIGNATURE_NAMES:

            collision_count = 0
            cross_count = 0
            all_count = 0
            total_survivors = 0

            for result in results:
                item = result["sources"].get(source)

                if item is None:
                    continue

                collision_count += 1

                survivors = len(item["matching"])
                total_survivors += survivors

                if survivors >= 2:
                    cross_count += 1

                if survivors == len(SIGNATURE_NAMES):
                    all_count += 1

            if collision_count:
                mean_survivors = (
                    total_survivors / collision_count
                )
            else:
                mean_survivors = 0.0

            print(
                f"{source:<15}"
                f"{collision_count:>12}"
                f"{cross_count:>12}"
                f"{all_count:>12}"
                f"{mean_survivors:>18.2f}"
            )

    # ------------------------------------------------------------------------
    # Pairwise conditional survival
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("PAIRWISE CROSS-SURVIVAL")
    print("=" * 78)

    print(
        "Entry X -> Y = number of X-collisions that also satisfy Y."
    )
    print()

    for prefix in PREFIXES:

        print()
        print(f"PREFIX {prefix}")
        print("-" * 78)

        header = f"{'source':<15}"
        for name in SIGNATURE_NAMES:
            header += f"{name[:8]:>10}"

        print(header)

        for source in SIGNATURE_NAMES:

            row = f"{source:<15}"

            for destination in SIGNATURE_NAMES:

                if source == destination:
                    value = "-"
                else:
                    total = 0
                    survive = 0

                    for result in all_results[prefix]:
                        item = result["sources"].get(source)

                        if item is None:
                            continue

                        total += 1

                        if destination in item["matching"]:
                            survive += 1

                    value = (
                        f"{survive}/{total}"
                        if total
                        else "0/0"
                    )

                row += f"{value:>10}"

            print(row)

    # ------------------------------------------------------------------------
    # Minimal separating subsets
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("MINIMAL SEPARATING SUBSETS")
    print("=" * 78)

    print(
        """
A subset is called separating for a particular source collision when
the candidate satisfies the source signature but fails at least one
signature in the subset.

The important quantity is:

    survivors = collisions satisfying ALL signatures in the subset

Smaller survivor counts mean stronger separation.
"""
    )

    for prefix in PREFIXES:

        print()
        print(f"PREFIX {prefix}")
        print("-" * 78)

        # Collect all source collision masks.
        masks = []

        for result in all_results[prefix]:
            for source in SIGNATURE_NAMES:
                item = result["sources"].get(source)

                if item is not None:
                    masks.append(item["mask"])

        if not masks:
            print("No collisions.")
            continue

        total_collisions = len(masks)

        print(
            f"total source collisions analyzed = {total_collisions}"
        )
        print()

        separating_subsets = []

        for size in range(1, len(SIGNATURE_NAMES) + 1):

            best_survivors = None
            best_subsets = []

            for subset in combinations(SIGNATURE_NAMES, size):

                subset_bits = 0

                for name in subset:
                    idx = SIGNATURE_NAMES.index(name)
                    subset_bits |= (1 << idx)

                survivors = sum(
                    1
                    for mask in masks
                    if (mask & subset_bits) == subset_bits
                )

                if best_survivors is None or survivors < best_survivors:
                    best_survivors = survivors
                    best_subsets = [subset]

                elif survivors == best_survivors:
                    best_subsets.append(subset)

            print(
                f"size {size}: "
                f"minimum survivors = "
                f"{best_survivors}/{total_collisions}"
            )

            for subset in best_subsets[:10]:
                print(
                    " " * 9 +
                    " + ".join(subset)
                )

            if best_survivors == 0:
                print()
                print(
                    f"FIRST ZERO-SURVIVAL SUBSET SIZE = {size}"
                )
                break

    # ------------------------------------------------------------------------
    # Full signature-mask distribution
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("SURVIVAL MASK DISTRIBUTION")
    print("=" * 78)

    for prefix in PREFIXES:

        counts = {}

        for result in all_results[prefix]:
            for source in SIGNATURE_NAMES:
                item = result["sources"].get(source)

                if item is None:
                    continue

                mask = item["mask"]
                counts[mask] = counts.get(mask, 0) + 1

        print()
        print(f"PREFIX {prefix}")
        print("-" * 78)

        for mask, count in sorted(
            counts.items(),
            key=lambda x: (-x[1], x[0]),
        ):
            names = mask_names(mask)

            print(
                f"{count:>5}  "
                f"{len(names):>2}/{len(SIGNATURE_NAMES)}  "
                f"{', '.join(names)}"
            )

    # ------------------------------------------------------------------------
    # Particularly interesting collision classes
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("HIGH-CORRELATION COLLISIONS")
    print("=" * 78)

    for prefix in PREFIXES:

        print()
        print(f"PREFIX {prefix}")
        print("-" * 78)

        for target_number, result in enumerate(
            all_results[prefix],
            1,
        ):
            for source in SIGNATURE_NAMES:

                item = result["sources"].get(source)

                if item is None:
                    continue

                if len(item["matching"]) >= 6:

                    candidate = item["candidate"]

                    print(
                        f"target {target_number:<2} "
                        f"source={source:<14} "
                        f"candidate=({candidate[0]}, {candidate[1]})"
                    )

                    print(
                        f"    survives "
                        f"{len(item['matching'])}/"
                        f"{len(SIGNATURE_NAMES)}: "
                        f"{', '.join(item['matching'])}"
                    )


# ============================================================================
# FINAL SUMMARY
# ============================================================================

def final_summary(all_results):

    print()
    print("=" * 78)
    print("FINAL SUMMARY")
    print("=" * 78)

    print()
    print(
        "FIRST ZERO-SURVIVAL PREFIX BY SOURCE"
    )
    print()

    for source in SIGNATURE_NAMES:

        zero_prefix = None

        for prefix in PREFIXES:

            found_collision = False
            surviving = 0

            for result in all_results[prefix]:

                item = result["sources"].get(source)

                if item is None:
                    continue

                found_collision = True

                if len(item["matching"]) == 1:
                    surviving += 1

            # This source has reached zero constructive collisions
            # only when there were no source collisions at all.
            if not found_collision:
                zero_prefix = prefix
                break

        if zero_prefix is None:
            print(
                f"{source:<15} "
                f"first zero-source-collision prefix = not reached"
            )
        else:
            print(
                f"{source:<15} "
                f"first zero-source-collision prefix = "
                f"{zero_prefix}"
            )

    print()
    print(
        "INTERPRETATION"
    )
    print(
        """
Experiment 28 is primarily about redundancy.

If:

    {Fpair, Fsum}

performs almost as well as:

    {Fpair, Fsorted, Fsum, Fprod, Fdiff, cube_sum,
     order_pair, order_sorted}

then most of those signatures are contributing little additional
information.

If adding:

    order_pair

causes a large drop in surviving collisions, then multiplicative
order is providing information that the residue-polynomial signatures
do not contain.

Likewise, if:

    cube_sum

produces a large additional reduction, it is carrying an independent
constraint.

The strongest result would be a very small subset with zero survivors.

That would identify a compact separating signature family for the
finite experiment.

However:

    zero survivors != mathematical uniqueness.

It only means that no collision was found inside the tested prime pool.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    t_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 28")
    print("SIGNATURE INTERSECTION / MINIMAL SEPARATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {POOL_SIZE:,}")
    print(f"prime interval   = [{PRIME_LOW:,}, {PRIME_HIGH:,}]")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # ------------------------------------------------------------------------
    # 1. MODULUS INVENTORY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. MODULUS INVENTORY")
    print("=" * 78)

    moduli = []

    for r in R_VALUES:
        m = modulus_for_r(r)
        u = unit_count(m)

        moduli.append(m)

        print(
            f"r={r:>3} "
            f"m={m:>7} "
            f"units={u:>7}"
        )

    # ------------------------------------------------------------------------
    # 2. PRIME POOL
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIME POOL GENERATION")
    print("=" * 78)

    t0 = time.perf_counter()

    pool = generate_prime_pool()

    dt = time.perf_counter() - t0

    print(f"generated primes = {len(pool):,}")
    print(f"pool generation time = {dt:.2f}s")

    # ------------------------------------------------------------------------
    # 3. TARGETS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. TARGET GENERATION")
    print("=" * 78)

    targets = generate_targets(pool)

    for i, (p, q) in enumerate(targets, 1):
        print(
            f"target {i:>2}: "
            f"p={p} "
            f"q={q} "
            f"n={p*q}"
        )

    # ------------------------------------------------------------------------
    # 4. RESIDUE CACHE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PRIME RESIDUE CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    residue_cache = build_residue_cache(
        pool,
        moduli,
    )

    # Target factors are already in the pool, but keep this robust if
    # target generation is changed later.
    for p, q in targets:
        if p not in residue_cache:
            residue_cache[p] = {
                m: p % m
                for m in moduli
            }

        if q not in residue_cache:
            residue_cache[q] = {
                m: q % m
                for m in moduli
            }

    dt = time.perf_counter() - t0

    print("residue cache complete")
    print(f"cache time = {dt:.2f}s")

    # ------------------------------------------------------------------------
    # 5. ORDER CACHE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. MULTIPLICATIVE ORDER CACHE")
    print("=" * 78)

    t0 = time.perf_counter()

    all_primes = set(pool)

    for p, q in targets:
        all_primes.add(p)
        all_primes.add(q)

    order_get = build_lazy_order_cache(
        all_primes,
        moduli,
        residue_cache,
    )

    print("using lazy order cache")
    print("order cache ready")

    dt = time.perf_counter() - t0

    print(f"cache setup time = {dt:.2f}s")

    # ------------------------------------------------------------------------
    # 6. EXPERIMENT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SIGNATURE INTERSECTION EXPERIMENT")
    print("=" * 78)

    all_results = run_experiment(
        pool,
        targets,
        residue_cache,
        order_get,
    )

    # ------------------------------------------------------------------------
    # 7. ANALYSIS
    # ------------------------------------------------------------------------

    analyze_intersections(all_results)

    # ------------------------------------------------------------------------
    # 8. SUMMARY
    # ------------------------------------------------------------------------

    final_summary(all_results)

    # ------------------------------------------------------------------------
    # RUNTIME
    # ------------------------------------------------------------------------

    elapsed = time.perf_counter() - t_start

    print()
    print("=" * 78)
    print(
        f"total runtime = {elapsed:.2f}s "
        f"({elapsed / 60:.2f} min)"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

