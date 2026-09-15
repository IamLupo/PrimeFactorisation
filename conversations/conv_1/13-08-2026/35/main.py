#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 35
JOINT SIGNATURE SEPARATION — INDEXED / NO PAIR ENUMERATION
NO CSV OUTPUT
==============================================================================

Purpose
-------
Experiment 34 showed that every individual signature has collisions.

Experiment 35 tests whether COMBINING signatures across multiple moduli
can uniquely identify the unordered prime pair.

IMPORTANT:
    We never enumerate all 5000 choose 2 prime pairs.

Instead:
    1. Generate the same prime pool.
    2. Compute residue/signature vectors for every prime.
    3. Build indexes for each individual signature.
    4. For every target, intersect candidate sets using progressively
       stronger joint signatures.
    5. Stop immediately when the candidate count reaches 1.

The experiment measures:

    Fpair
    Fsorted
    Fdiff
    Fsum
    Fprod
    cube_sum
    order_pair

both individually and jointly across increasing prefixes.

==============================================================================

RESULT INTERPRETATION
---------------------

UNIQUE:
    exactly one unordered prime pair remains.

WRONG:
    target signature matches multiple prime pairs.

MISSING:
    target itself is not represented in the pool/index.

The important measurement is:

    JOINT candidate count

If this reaches 1 consistently, the modular signature family is promising.

If it remains large even after all moduli/signatures are combined, then
the current signature family is fundamentally too weak.

==============================================================================
"""

import math
import random
import time
from collections import defaultdict


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

TARGETS = 12
PRIME_POOL_SIZE = 5000

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

PREFIXES = [3, 4, 5, 6, 7]

R_VALUES = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

SIGNATURES = [
    "Fpair",
    "Fsorted",
    "Fdiff",
    "Fsum",
    "Fprod",
    "cube_sum",
    "order_pair",
]

# How many candidate signatures to combine at each stage.
#
# Start with the strongest/cheapest combinations and progressively add
# information.
JOINT_STAGES = [
    ("Fpair",),
    ("Fpair", "Fdiff"),
    ("Fpair", "Fdiff", "Fsum"),
    ("Fpair", "Fdiff", "Fsum", "Fprod"),
    ("Fpair", "Fdiff", "Fsum", "Fprod", "cube_sum"),
    ("Fpair", "Fdiff", "Fsum", "Fprod", "cube_sum", "order_pair"),
]


# ============================================================================
# PRIME GENERATION
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
    """
    Generate a deterministic pool of primes.

    We use the same deterministic random-selection idea throughout the
    experiment so the result is reproducible.
    """

    rng = random.Random(SEED)

    candidates = []

    # Generate enough random candidates first.
    #
    # The interval is only ~2.2M wide, so this is inexpensive.
    seen = set()

    while len(candidates) < PRIME_POOL_SIZE:
        x = rng.randint(PRIME_LO, PRIME_HI)

        if x in seen:
            continue

        seen.add(x)

        if is_prime(x):
            candidates.append(x)

    candidates.sort()

    return candidates


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(primes):
    """
    Select deterministic target pairs from the prime pool.

    Targets are deliberately drawn from the pool so that MISSING means
    a signature/indexing problem rather than a pool-membership problem.
    """

    rng = random.Random(SEED + 999)

    targets = []

    used = set()

    while len(targets) < TARGETS:
        p, q = rng.sample(primes, 2)

        pair = tuple(sorted((p, q)))

        if pair in used:
            continue

        used.add(pair)
        targets.append(pair)

    return targets


# ============================================================================
# MODULUS STRUCTURE
# ============================================================================

def modulus_for_r(r):
    """
    m = 3r + 1

    This reproduces the modulus inventory:

        r=2  -> 7
        r=3  -> 10
        r=5  -> 16
        ...
    """

    return 3 * r + 1


def selected_r_values(prefix):
    return R_VALUES[:prefix]


def selected_moduli(prefix):
    return [modulus_for_r(r) for r in selected_r_values(prefix)]


# ============================================================================
# MULTIPLICATIVE ORDER
# ============================================================================

def multiplicative_order(a, m):
    """
    Return multiplicative order of a mod m when gcd(a,m)=1.

    Return 0 when a is not a unit.

    For these small moduli this is extremely cheap.
    """

    if math.gcd(a, m) != 1:
        return 0

    x = 1

    for k in range(1, m + 1):
        x = (x * a) % m

        if x == 1:
            return k

    return 0


def build_order_table(m):
    return {
        a: multiplicative_order(a, m)
        for a in range(m)
    }


# ============================================================================
# PRIME RESIDUE CACHE
# ============================================================================

def build_residue_cache(primes, moduli):
    """
    residue_cache[m][p] = p mod m
    """

    cache = {}

    for m in moduli:
        cache[m] = {
            p: p % m
            for p in primes
        }

    return cache


# ============================================================================
# SINGLE-MODULUS SIGNATURES
# ============================================================================

def signature_values(a, b, m, order_table):
    """
    Compute all seven signatures for an unordered pair of residues.

    IMPORTANT:
    These are deliberately represented as canonical tuples.

    Fpair:
        ordered pair (a,b)

    Fsorted:
        unordered residue pair

    Fdiff:
        absolute difference

    Fsum:
        sum modulo m

    Fprod:
        product modulo m

    cube_sum:
        a^3 + b^3 modulo m

    order_pair:
        multiplicative orders of the two residues, sorted
    """

    a %= m
    b %= m

    lo, hi = sorted((a, b))

    return {
        "Fpair": (a, b),
        "Fsorted": (lo, hi),

        "Fdiff": (hi - lo) % m,

        "Fsum": (a + b) % m,

        "Fprod": (a * b) % m,

        "cube_sum": (pow(a, 3, m) + pow(b, 3, m)) % m,

        "order_pair": tuple(
            sorted(
                (
                    order_table[a],
                    order_table[b],
                )
            )
        ),
    }


# ============================================================================
# TARGET SIGNATURE VECTOR
# ============================================================================

def target_signature_vector(
    p,
    q,
    moduli,
    order_tables,
):
    """
    Build the complete signature vector for one target pair.

    Structure:

        {
            m: {
                signature_name: signature_value
            }
        }
    """

    result = {}

    for m in moduli:
        result[m] = signature_values(
            p % m,
            q % m,
            m,
            order_tables[m],
        )

    return result


# ============================================================================
# PRIME-PAIR INDEX
# ============================================================================

def build_signature_indexes(primes, moduli, order_tables):
    """
    Build indexes without enumerating prime pairs.

    For each modulus and each signature, store:

        signature value -> set of prime IDs

    We do NOT build a 5000 choose 2 table.

    Candidate pairs are recovered by matching the two target residue
    signatures against the indexed prime classes.
    """

    indexes = {}

    for m in moduli:

        indexes[m] = {}

        # First build the seven prime-level features that can be matched
        # against a target pair.
        #
        # Some pair signatures cannot be decomposed into a single-prime
        # feature. Therefore this index stores the residue class itself,
        # and candidate pairs are formed only inside matching residue
        # classes.

        residue_index = defaultdict(list)

        for pid, p in enumerate(primes):
            residue_index[p % m].append(pid)

        indexes[m]["residue"] = residue_index

    return indexes


# ============================================================================
# PAIR SIGNATURE TEST
# ============================================================================

def pair_matches_signature(
    p,
    q,
    target_p,
    target_q,
    m,
    signature_name,
    order_table,
):
    """
    Test whether (p,q) has the same signature as the target.

    Everything is done modulo m.
    """

    actual = signature_values(
        p % m,
        q % m,
        m,
        order_table,
    )[signature_name]

    target = signature_values(
        target_p % m,
        target_q % m,
        m,
        order_table,
    )[signature_name]

    return actual == target


# ============================================================================
# FAST CANDIDATE GENERATION
# ============================================================================

def residue_pair_candidates(
    primes,
    residue_index,
    target_p,
    target_q,
    m,
):
    """
    Generate candidate prime IDs for a target residue pair.

    This is still vastly smaller than scanning every pair.

    For Fsorted / Fdiff / Fsum / etc. we need to test combinations of
    residue classes that satisfy the target signature.

    Since m <= 52, there are at most m^2 residue combinations.

    We therefore work at residue-class level first.
    """

    tp = target_p % m
    tq = target_q % m

    target_residues = set()

    # All residue pairs satisfying the same unordered residue pair.
    #
    # This is intentionally broad. Later signature filters remove the
    # remaining combinations.
    target_residues.add(tuple(sorted((tp, tq))))

    candidate_ids = []

    for a, b in target_residues:

        left = residue_index.get(a, [])
        right = residue_index.get(b, [])

        if a == b:
            for i in range(len(left)):
                for j in range(i + 1, len(left)):
                    candidate_ids.append((left[i], left[j]))
        else:
            for i in left:
                for j in right:
                    if i < j:
                        candidate_ids.append((i, j))

    return candidate_ids


# ============================================================================
# GENERAL JOINT FILTER
# ============================================================================

def filter_candidates(
    candidates,
    primes,
    target_p,
    target_q,
    moduli,
    signature_names,
    order_tables,
):
    """
    Filter an existing candidate list against a joint signature.

    This is the key operation.

    Once the first modulus has reduced the candidate space, all subsequent
    filters operate only on the survivors.
    """

    survivors = []

    for i, j in candidates:

        p = primes[i]
        q = primes[j]

        ok = True

        for m in moduli:

            target_sig = signature_values(
                target_p % m,
                target_q % m,
                m,
                order_tables[m],
            )

            actual_sig = signature_values(
                p % m,
                q % m,
                m,
                order_tables[m],
            )

            for name in signature_names:

                if actual_sig[name] != target_sig[name]:
                    ok = False
                    break

            if not ok:
                break

        if ok:
            survivors.append((i, j))

    return survivors


# ============================================================================
# SMART INITIAL CANDIDATE SET
# ============================================================================

def initial_candidates(
    primes,
    target_p,
    target_q,
    m,
    residue_index,
):
    """
    Initial candidate set from unordered residue pair.

    This is equivalent to starting from Fsorted.

    It avoids scanning all prime pairs.
    """

    a = target_p % m
    b = target_q % m

    lo, hi = sorted((a, b))

    left = residue_index.get(lo, [])
    right = residue_index.get(hi, [])

    result = []

    if lo == hi:

        n = len(left)

        for x in range(n):
            for y in range(x + 1, n):
                i = left[x]
                j = left[y]

                if i != j:
                    result.append(
                        (min(i, j), max(i, j))
                    )

    else:

        for i in left:
            for j in right:

                if i == j:
                    continue

                result.append(
                    (min(i, j), max(i, j))
                )

    return result


# ============================================================================
# DIAGNOSTIC FOR ONE TARGET
# ============================================================================

def diagnose_target(
    target_number,
    target,
    primes,
    prefix,
    indexes,
    order_tables,
):
    p, q = target

    moduli = selected_moduli(prefix)

    # Use the largest modulus in this prefix as the initial narrowing
    # point. This keeps the initial candidate set small.
    #
    # Example:
    # prefix 3 -> m=16
    # prefix 7 -> m=52
    initial_m = moduli[-1]

    residue_index = indexes[initial_m]["residue"]

    t0 = time.perf_counter()

    candidates = initial_candidates(
        primes,
        p,
        q,
        initial_m,
        residue_index,
    )

    initial_count = len(candidates)

    elapsed = time.perf_counter() - t0

    print()
    print(
        f"prefix {prefix} "
        f"(m={initial_m})"
    )

    print(
        f"  initial Fsorted candidates: "
        f"{initial_count:,} "
        f"[{elapsed:.3f}s]"
    )

    if initial_count == 0:
        print("  MISSING at initial residue class")
        return None

    # ------------------------------------------------------------------
    # Progressive joint filtering
    # ------------------------------------------------------------------

    best = None

    for stage_no, signature_names in enumerate(
        JOINT_STAGES,
        start=1,
    ):

        t0 = time.perf_counter()

        survivors = filter_candidates(
            candidates,
            primes,
            p,
            q,
            moduli,
            signature_names,
            order_tables,
        )

        elapsed = time.perf_counter() - t0

        count = len(survivors)

        label = " + ".join(signature_names)

        if count == 0:
            status = "MISSING"
        elif count == 1:
            status = "UNIQUE"
        else:
            status = "COLLISION"

        print(
            f"  stage {stage_no}: "
            f"{label:<58} "
            f"{count:>10,}  "
            f"{status:<9} "
            f"[{elapsed:.3f}s]"
        )

        best = survivors

        # Once unique, there is no reason to do more work.
        if count == 1:
            i, j = survivors[0]

            cp, cq = sorted(
                (primes[i], primes[j])
            )

            print(
                f"      UNIQUE PAIR = "
                f"({cp}, {cq})"
            )

            print(
                f"      target       = "
                f"({p}, {q})"
            )

            if (cp, cq) == (min(p, q), max(p, q)):
                print("      verification = CORRECT")
            else:
                print("      verification = WRONG")

            break

        # If no survivors remain, stop immediately.
        if count == 0:
            break

        candidates = survivors

    return best


# ============================================================================
# GLOBAL SUMMARY
# ============================================================================

def print_summary(results):
    print()
    print("=" * 78)
    print("GLOBAL SUMMARY — JOINT SIGNATURE SEPARATION")
    print("=" * 78)

    print()
    print(
        "prefix      target      initial       final       status"
    )
    print("-" * 78)

    for row in results:

        print(
            f"{row['prefix']:>6} "
            f"{row['target']:>10} "
            f"{row['initial']:>12,} "
            f"{row['final']:>12,} "
            f"{row['status']}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    experiment_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 35")
    print("JOINT SIGNATURE SEPARATION — INDEXED / NO PAIR ENUMERATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime pool       = {PRIME_POOL_SIZE:,}")
    print(
        f"prime interval   = "
        f"[{PRIME_LO:,}, {PRIME_HI:,}]"
    )
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")

    # ------------------------------------------------------------------
    # PRIME POOL
    # ------------------------------------------------------------------

    print()
    print("-" * 78)
    print("PRIME POOL")
    print("-" * 78)

    t0 = time.perf_counter()

    primes = generate_prime_pool()

    print(
        f"generated primes = {len(primes):,}"
    )

    print(
        f"generation time  = "
        f"{time.perf_counter() - t0:.3f}s"
    )

    # ------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------

    targets = generate_targets(primes)

    print()
    print("-" * 78)
    print("TARGETS")
    print("-" * 78)

    for k, (p, q) in enumerate(targets, start=1):

        print(
            f"target {k:2}: "
            f"p={p} "
            f"q={q} "
            f"n={p*q}"
        )

    # ------------------------------------------------------------------
    # MODULI
    # ------------------------------------------------------------------

    all_moduli = sorted(
        set(
            modulus_for_r(r)
            for r in R_VALUES
        )
    )

    print()
    print("-" * 78)
    print("MODULUS INVENTORY")
    print("-" * 78)

    for r in R_VALUES:

        m = modulus_for_r(r)

        units = sum(
            math.gcd(x, m) == 1
            for x in range(m)
        )

        print(
            f"r={r:3} "
            f"m={m:7} "
            f"units={units:7}"
        )

    # ------------------------------------------------------------------
    # ORDER TABLES
    # ------------------------------------------------------------------

    print()
    print("-" * 78)
    print("MULTIPLICATIVE ORDER TABLES")
    print("-" * 78)

    t0 = time.perf_counter()

    order_tables = {
        m: build_order_table(m)
        for m in all_moduli
    }

    print(
        f"order tables ready "
        f"[{time.perf_counter() - t0:.3f}s]"
    )

    # ------------------------------------------------------------------
    # RESIDUE INDEXES
    # ------------------------------------------------------------------

    print()
    print("-" * 78)
    print("BUILDING RESIDUE INDEXES")
    print("-" * 78)

    t0 = time.perf_counter()

    indexes = build_signature_indexes(
        primes,
        all_moduli,
        order_tables,
    )

    print(
        f"residue indexes ready "
        f"[{time.perf_counter() - t0:.3f}s]"
    )

    # ------------------------------------------------------------------
    # EXPERIMENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("JOINT SIGNATURE DIAGNOSTIC")
    print("=" * 78)

    results = []

    for target_number, target in enumerate(
        targets,
        start=1,
    ):

        print()
        print("-" * 78)
        print(
            f"TARGET {target_number:2} "
            f"({target[0]}, {target[1]})"
        )
        print("-" * 78)

        for prefix in PREFIXES:

            final = diagnose_target(
                target_number,
                target,
                primes,
                prefix,
                indexes,
                order_tables,
            )

            if final is None:
                final_count = 0
                initial_count = 0
                status = "MISSING"

            else:
                final_count = len(final)

                # Recover initial count cheaply.
                m = selected_moduli(prefix)[-1]

                initial = initial_candidates(
                    primes,
                    target[0],
                    target[1],
                    m,
                    indexes[m]["residue"],
                )

                initial_count = len(initial)

                if final_count == 0:
                    status = "MISSING"
                elif final_count == 1:
                    pair = tuple(
                        sorted(
                            (
                                primes[final[0][0]],
                                primes[final[0][1]],
                            )
                        )
                    )

                    target_pair = tuple(sorted(target))

                    if pair == target_pair:
                        status = "UNIQUE/CORRECT"
                    else:
                        status = "UNIQUE/WRONG"
                else:
                    status = "COLLISION"

            results.append({
                "target": target_number,
                "prefix": prefix,
                "initial": initial_count,
                "final": final_count,
                "status": status,
            })

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print_summary(results)

    print()
    print("=" * 78)
    print("EXPERIMENT 35 COMPLETE")
    print("=" * 78)

    total = time.perf_counter() - experiment_start

    print(
        f"total runtime = {total:.3f}s"
    )

    print()
    print("KEY QUESTION:")
    print(
        "Does the joint signature across the modular prefixes "
        "reduce the target to exactly one prime pair?"
    )


if __name__ == "__main__":
    main()

