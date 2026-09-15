#!/usr/bin/env python3

import math
import random
import time
from collections import Counter

# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 1_511_464_998

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

# Recursive descendant radices.
# These are deliberately much smaller than r1/r2 so that the quotient
# coordinates can actually descend several levels.
INNER_MIN = 3
INNER_MAX = 31

DEPTH = 4

PRINT_ANCHOR_EVERY = 10

# Limit candidate roots so a pathological anchor cannot explode the tree.
MAX_ROOT_CANDIDATES = 20_000


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def is_prime(n: int) -> bool:
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


def prime_list(lo: int, hi: int):
    return [x for x in range(lo, hi + 1) if is_prime(x)]


def factor_integer(n: int):
    """
    Exact factorization of the small quotient products occurring here.

    Returns:
        [(prime, exponent), ...]
    """
    if n < 1:
        return []

    factors = []

    d = 2

    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            factors.append((d, e))

        d = 3 if d == 2 else d + 2

    if n > 1:
        factors.append((n, 1))

    return factors


def divisors_from_factorization(factors):
    divisors = [1]

    for p, e in factors:
        old = divisors[:]
        mul = 1

        for _ in range(e):
            mul *= p
            for d in old:
                divisors.append(d * mul)

    return sorted(divisors)


def divisor_pairs(n: int):
    """
    All ordered positive factor pairs (x,y) with x*y=n.
    """
    if n <= 0:
        return []

    out = []

    for d in divisors_from_factorization(factor_integer(n)):
        if d * d > n:
            break

        if n % d == 0:
            q = n // d

            out.append((d, q))

            if d != q:
                out.append((q, d))

    return sorted(out)


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def build_factor_pool():
    return prime_list(FACTOR_MIN, FACTOR_MAX)


def build_modulus_pool():
    return prime_list(MOD_MIN, MOD_MAX)


def choose_close_pair(mods, rng):
    """
    Pick a random close prime pair.
    """
    candidates = []

    for r1 in mods:
        lo = r1
        hi = int(math.floor(r1 * (1.0 + CLOSE_RATIO)))

        for r2 in mods:
            if r2 == r1:
                continue

            if lo <= r2 <= hi:
                candidates.append((r1, r2))

    if not candidates:
        raise RuntimeError("No close modulus pairs available.")

    return rng.choice(candidates)


def build_anchors(factor_primes, modulus_primes, count, seed):
    rng = random.Random(seed)

    anchors = []

    for _ in range(count):
        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        n = p * q

        r1, r2 = choose_close_pair(modulus_primes, rng)

        if r1 > r2:
            r1, r2 = r2, r1

        k = p // r1
        l = q // r2

        if k <= 0 or l <= 0:
            continue

        anchors.append({
            "p": p,
            "q": q,
            "n": n,
            "r1": r1,
            "r2": r2,
            "k": k,
            "l": l,
            "K": k * l,
        })

    if not anchors:
        raise RuntimeError("No anchors created.")

    return anchors


# =============================================================================
# ROOT QUOTIENT PRODUCT
# =============================================================================

def observable_root_window(anchor):
    n = anchor["n"]
    r1 = anchor["r1"]
    r2 = anchor["r2"]

    pmax = FACTOR_MAX
    qmax = FACTOR_MAX

    # From:
    #
    # n/(r1*r2) = kl + positive cross terms
    #
    # and the previous bound E <= k+l,
    # obtain a safe observable upper bound from factor-range limits.

    T = n // (r1 * r2)

    kmax = pmax // r1
    lmax = qmax // r2

    emax = kmax + lmax

    low = max(1, T - emax)
    high = T

    return T, low, high, emax


def root_candidates(anchor):
    """
    Generate candidate (k,l) pairs from the observable T-E window.

    We don't use the true K here.
    """
    T, low, high, emax = observable_root_window(anchor)

    r1 = anchor["r1"]
    r2 = anchor["r2"]

    kmin = max(1, math.ceil(FACTOR_MIN / r1))
    kmax = FACTOR_MAX // r1

    lmin = max(1, math.ceil(FACTOR_MIN / r2))
    lmax = FACTOR_MAX // r2

    candidates = set()

    for K in range(low, high + 1):
        for k, l in divisor_pairs(K):

            if kmin <= k <= kmax and lmin <= l <= lmax:
                candidates.add((k, l))

                if len(candidates) >= MAX_ROOT_CANDIDATES:
                    return sorted(candidates)

    return sorted(candidates)


# =============================================================================
# RECURSIVE DESCENT
# =============================================================================

def deterministic_inner_radices(anchor_id, level):
    """
    Deterministic per-anchor / per-level radices.

    These are deliberately small. They form a recursive coordinate system:

        x = u*s1 + a
        y = v*s2 + b
    """
    seed = (
        SEED
        ^ (anchor_id * 0x9E3779B1)
        ^ (level * 0x85EBCA77)
    ) & 0xFFFFFFFF

    rng = random.Random(seed)

    primes = prime_list(INNER_MIN, INNER_MAX)

    s1 = rng.choice(primes)
    s2 = rng.choice(primes)

    return s1, s2


def recursive_step(x, y, s1, s2):
    """
    Given a coordinate pair (x,y), construct its next quotient coordinates.

        x = u*s1 + a
        y = v*s2 + b

    Then:

        U = u*v

    The returned data describes the next-generation quotient product.
    """
    if x <= 0 or y <= 0:
        return None

    u, ax = divmod(x, s1)
    v, ay = divmod(y, s2)

    if u <= 0 or v <= 0:
        return {
            "u": u,
            "v": v,
            "a": ax,
            "b": ay,
            "product": u * v,
            "factor_pairs": [],
        }

    U = u * v
    pairs = divisor_pairs(U)

    return {
        "u": u,
        "v": v,
        "a": ax,
        "b": ay,
        "product": U,
        "factor_pairs": pairs,
    }


def genealogy(anchor_id, k, l):
    """
    Follow one candidate branch recursively.

    We intentionally do NOT use n after the root stage.
    The question is whether the quotient-product genealogy itself
    contains a distinctive path.
    """
    x = k
    y = l

    levels = []

    for depth in range(1, DEPTH + 1):
        s1, s2 = deterministic_inner_radices(anchor_id, depth)

        step = recursive_step(x, y, s1, s2)

        if step is None:
            break

        level = {
            "depth": depth,
            "x": x,
            "y": y,
            "s1": s1,
            "s2": s2,
            "u": step["u"],
            "v": step["v"],
            "a": step["a"],
            "b": step["b"],
            "product": step["product"],
            "factor_pair_count": len(step["factor_pairs"]),
        }

        levels.append(level)

        if step["u"] <= 0 or step["v"] <= 0:
            break

        x = step["u"]
        y = step["v"]

    return levels


def genealogy_signature(levels):
    """
    A compact signature of the recursive descendant structure.

    This deliberately ignores p and q.
    """
    return tuple(
        (
            level["factor_pair_count"],
            level["product"],
            level["u"],
            level["v"],
        )
        for level in levels
    )


# =============================================================================
# ROOT ANALYSIS
# =============================================================================

def analyze_anchor(anchor_id, anchor):
    p = anchor["p"]
    q = anchor["q"]
    n = anchor["n"]

    r1 = anchor["r1"]
    r2 = anchor["r2"]

    true_k = anchor["k"]
    true_l = anchor["l"]
    true_K = anchor["K"]

    # -------------------------------------------------------------------------
    # Observable root quantity
    # -------------------------------------------------------------------------

    T, low, high, emax = observable_root_window(anchor)

    # -------------------------------------------------------------------------
    # Exact/oracle factorization of true K = k*l
    # -------------------------------------------------------------------------

    true_pairs = divisor_pairs(true_K)

    admissible_true_pairs = []

    kmin = max(1, math.ceil(FACTOR_MIN / r1))
    kmax = FACTOR_MAX // r1

    lmin = max(1, math.ceil(FACTOR_MIN / r2))
    lmax = FACTOR_MAX // r2

    for x, y in true_pairs:
        if kmin <= x <= kmax and lmin <= y <= lmax:
            admissible_true_pairs.append((x, y))

    true_pair_unique = len(admissible_true_pairs) == 1

    # -------------------------------------------------------------------------
    # Observable candidate generation
    # -------------------------------------------------------------------------

    candidates = root_candidates(anchor)

    true_present = (true_k, true_l) in candidates

    # -------------------------------------------------------------------------
    # Genealogy signatures for observable candidates
    # -------------------------------------------------------------------------

    signature_counts = Counter()

    candidate_signatures = {}

    for x, y in candidates:
        levels = genealogy(anchor_id, x, y)
        sig = genealogy_signature(levels)

        candidate_signatures[(x, y)] = sig
        signature_counts[sig] += 1

    true_signature = None

    if true_present:
        true_signature = candidate_signatures[(true_k, true_l)]

    # Is true genealogy unique among observable candidates?
    true_signature_unique = (
        true_present
        and
        signature_counts[true_signature] == 1
    )

    # -------------------------------------------------------------------------
    # Recursive branching on the TRUE branch
    # -------------------------------------------------------------------------

    true_genealogy = genealogy(anchor_id, true_k, true_l)

    branch_counts = [
        level["factor_pair_count"]
        for level in true_genealogy
    ]

    # -------------------------------------------------------------------------
    # Competing genealogy complexity
    # -------------------------------------------------------------------------

    candidate_depth_lengths = []

    for levels in candidate_signatures.values():
        candidate_depth_lengths.append(len(levels))

    avg_signature_depth = (
        sum(candidate_depth_lengths) / len(candidate_depth_lengths)
        if candidate_depth_lengths
        else 0.0
    )

    return {
        "p": p,
        "q": q,
        "n": n,
        "r1": r1,
        "r2": r2,
        "true_k": true_k,
        "true_l": true_l,
        "true_K": true_K,
        "T": T,
        "E_true": T - true_K,
        "E_max": emax,
        "root_factor_pairs": len(true_pairs),
        "root_admissible_pairs": len(admissible_true_pairs),
        "root_true_unique": true_pair_unique,
        "observable_candidates": len(candidates),
        "true_present": true_present,
        "true_signature_unique": true_signature_unique,
        "branch_counts": branch_counts,
        "true_genealogy": true_genealogy,
        "signature_count": len(signature_counts),
        "avg_signature_depth": avg_signature_depth,
    }


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def run():
    start_time = time.perf_counter()

    print("=" * 100)
    print("RECURSIVE QUOTIENT-PRODUCT GENEALOGY / FACTOR-DESCENT EXPERIMENT")
    print("=" * 100)

    print(f"N anchors                 = {N_ANCHORS:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"recursive depth           = {DEPTH}")
    print(f"inner radix range         = {INNER_MIN} - {INNER_MAX}")
    print(f"seed                      = {SEED:,}")

    print()
    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = build_factor_pool()
    modulus_primes = build_modulus_pool()

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    print()
    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        modulus_primes,
        N_ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")

    print()
    print("=" * 100)
    print("RUNNING RECURSIVE QUOTIENT GENEALOGY")
    print("=" * 100)

    results = []

    for i, anchor in enumerate(anchors, start=1):
        result = analyze_anchor(i, anchor)
        results.append(result)

        if (
            i % PRINT_ANCHOR_EVERY == 0
            or i == len(anchors)
        ):
            print(f"anchor {i:3d}/{len(anchors)}")

    # -------------------------------------------------------------------------
    # Aggregate statistics
    # -------------------------------------------------------------------------

    avg_observable_candidates = (
        sum(r["observable_candidates"] for r in results)
        / len(results)
    )

    avg_root_pairs = (
        sum(r["root_factor_pairs"] for r in results)
        / len(results)
    )

    avg_admissible_root_pairs = (
        sum(r["root_admissible_pairs"] for r in results)
        / len(results)
    )

    true_root_unique = sum(
        r["root_true_unique"] for r in results
    )

    true_present = sum(
        r["true_present"] for r in results
    )

    true_signature_unique = sum(
        r["true_signature_unique"] for r in results
    )

    avg_E = (
        sum(r["E_true"] for r in results)
        / len(results)
    )

    avg_Emax = (
        sum(r["E_max"] for r in results)
        / len(results)
    )

    all_branch_counts = []

    for r in results:
        all_branch_counts.extend(r["branch_counts"])

    avg_branch = (
        sum(all_branch_counts) / len(all_branch_counts)
        if all_branch_counts
        else 0.0
    )

    max_branch = max(all_branch_counts) if all_branch_counts else 0

    # -------------------------------------------------------------------------
    # Signature statistics
    # -------------------------------------------------------------------------

    signature_histogram = Counter(
        r["signature_count"]
        for r in results
    )

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(f"anchors analyzed                 = {len(results):,}")
    print(f"average true E                   = {avg_E:.3f}")
    print(f"average E-window upper bound     = {avg_Emax:.3f}")

    print()
    print("ORACLE TRUE-K ANALYSIS")
    print(f"average factor pairs of k*l      = {avg_root_pairs:.3f}")
    print(
        f"average admissible factor pairs  = "
        f"{avg_admissible_root_pairs:.3f}"
    )
    print(
        f"true (k,l) uniquely determined   = "
        f"{true_root_unique}/{len(results)}"
    )

    print()
    print("OBSERVABLE T-E ANALYSIS")
    print(
        f"average observable root pairs    = "
        f"{avg_observable_candidates:.3f}"
    )
    print(
        f"true (k,l) contained             = "
        f"{true_present}/{len(results)}"
    )

    print()
    print("RECURSIVE GENEALOGY")
    print(
        f"true genealogy branch avg        = "
        f"{avg_branch:.3f}"
    )
    print(f"maximum descendant branching     = {max_branch}")
    print(
        f"true genealogy signature unique  = "
        f"{true_signature_unique}/{len(results)}"
    )

    print()
    print("=" * 100)
    print("KEY TEST")
    print("=" * 100)

    if true_present:
        print(
            "The root product factorization preserves the true quotient "
            "pair for a subset of anchors."
        )
    else:
        print(
            "The defect window fails to recover the true quotient pair "
            "for at least some anchors."
        )

    if true_root_unique > 0:
        print(
            "Some anchors have a uniquely determined admissible "
            "(k,l) from factorization of the true k*l product."
        )

    if true_signature_unique > 0:
        print(
            "Some anchors produce a recursively distinctive genealogy."
        )

    print()
    print("=" * 100)
    print("ANCHOR EXAMPLES")
    print("=" * 100)

    examples = results[:20]

    for r in examples:
        branch = ",".join(map(str, r["branch_counts"]))

        print(
            f"n={r['n']:,} "
            f"p={r['p']:,} "
            f"q={r['q']:,} "
            f"mods=({r['r1']},{r['r2']}) "
            f"k={r['true_k']} "
            f"l={r['true_l']} "
            f"K={r['true_K']} "
            f"T={r['T']} "
            f"E={r['E_true']} "
            f"Emax={r['E_max']} "
            f"rootPairs={r['root_factor_pairs']} "
            f"admPairs={r['root_admissible_pairs']} "
            f"obsPairs={r['observable_candidates']} "
            f"present={r['true_present']} "
            f"branches=[{branch}] "
            f"sigUnique={r['true_signature_unique']}"
        )

    print()
    print("=" * 100)
    print("GENEALOGY DEPTH EXAMPLES")
    print("=" * 100)

    shown = 0

    for r in results:
        if shown >= 10:
            break

        print()
        print(
            f"n={r['n']:,} "
            f"true=(k,l)=({r['true_k']},{r['true_l']}) "
            f"K={r['true_K']}"
        )

        for level in r["true_genealogy"]:
            print(
                f"  depth={level['depth']} "
                f"x={level['x']} "
                f"y={level['y']} "
                f"radix=({level['s1']},{level['s2']}) "
                f"next=({level['u']},{level['v']}) "
                f"product={level['product']} "
                f"factorPairs={level['factor_pair_count']} "
                f"lowDigits=({level['a']},{level['b']})"
            )

        shown += 1

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
The experiment separates three different questions.

1. ORACLE PRODUCT INFORMATION

   If the exact quantity

       K = k*l

   were known, factor K and enumerate all divisor pairs.

   If the true (k,l) pair is frequently unique among the admissible
   divisor pairs, then K really does encode substantial information
   about the quotient coordinates.

2. OBSERVABLE PRODUCT INFORMATION

   In an actual attack we do not know K.

   We only know

       T = floor(n/(r1*r2))

   and

       K = T-E.

   Therefore the experiment factors every value in the finite
   defect window and collects all admissible quotient pairs.

   This measures how much information is lost by the unknown defect E.

3. RECURSIVE GENEALOGY

   For a quotient pair

       (x,y)

   choose smaller radices s1,s2 and write

       x = u*s1 + a
       y = v*s2 + b.

   The descendant quotient product is

       U = u*v.

   We then factor U and count its possible divisor-pair descendants.

   The important question is whether the true branch becomes
   structurally special as the recursion continues.

A genuinely interesting positive result would look like:

       n
        |
        v
       K = k*l
        |
        v
     unique/small (k,l)
        |
        v
     unique/small (u,v)
        |
        v
     unique/small descendants
        |
        v
     eventually tiny coordinates

while the competing divisor genealogies die out.

A negative result would be:

       K
       |
       +-- many divisor pairs
       |
       +-- many equally plausible descendants
       |
       +-- no unique genealogy.

That would show that recursive factoring of k*l does not
automatically reveal the parent coordinates.

The most important distinction is:

       "K is smaller than n"

is not enough.

We need:

       "factor(K) determines useful information about the
        ancestor coordinates."

That is exactly what this experiment measures.

No knowledge of p or q is used during candidate generation.
The true p and q are used only afterward to score whether the
genealogy retained the actual branch.
"""
    )

    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    elapsed = time.perf_counter() - start_time

    print(f"total runtime                  = {elapsed:.3f} seconds")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
