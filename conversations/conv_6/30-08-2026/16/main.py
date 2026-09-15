#!/usr/bin/env python3

import math
import random
from collections import defaultdict, Counter

# =================================================================================================
# CONFIG
# =================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

ACTUAL_ANCHORS = 300

SEED = 1_511_464_998

PROGRESS_EVERY = 25

# =================================================================================================
# PRIME TEST
# =================================================================================================

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


# =================================================================================================
# PRIME POOLS
# =================================================================================================

def build_primes(lo: int, hi: int):
    return [n for n in range(lo, hi + 1) if is_prime(n)]


# =================================================================================================
# DISTINCT PRIME PAIRS
# =================================================================================================

def build_actual_anchors(primes, count, rng):
    anchors = []

    seen = set()

    while len(anchors) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))
        anchors.append((p, q))

    return anchors


# =================================================================================================
# MODULUS PRIME POOL
# =================================================================================================

def build_modulus_primes(lo, hi):
    return build_primes(lo, hi)


# =================================================================================================
# MAXIMUM-PRODUCT CLOSE TRIPLE
# =================================================================================================

def choose_max_product_triple(n, modulus_primes, close_ratio):
    """
    Find r1 < r2 < r3 such that:

        r1*r2*r3 < n

    and the three primes are mutually close:

        max(r_i)/min(r_i) <= 1 + close_ratio

    Among all such triples, choose the one with maximum product.
    """

    best = None
    best_product = -1

    L = len(modulus_primes)

    for i in range(L):
        r1 = modulus_primes[i]

        # Lower bound on the largest allowed prime.
        max_allowed = int(r1 * (1.0 + close_ratio))

        for j in range(i + 1, L):
            r2 = modulus_primes[j]

            if r2 > max_allowed:
                break

            partial = r1 * r2

            if partial >= n:
                break

            # We only care about r3 > r2.
            for k in range(j + 1, L):
                r3 = modulus_primes[k]

                if r3 > max_allowed:
                    break

                product = partial * r3

                if product >= n:
                    break

                if product > best_product:
                    best_product = product
                    best = (r1, r2, r3)

    return best


# =================================================================================================
# RESIDUE / FINGERPRINT FUNCTIONS
# =================================================================================================

def raw_fingerprint(p, mods):
    return tuple(p % r for r in mods)


def sum_product_fingerprint(p, mods):
    a = p % mods[0]
    b = p % mods[1]
    c = p % mods[2]

    e1 = a + b + c
    e3 = a * b * c

    return e1, e3


def symmetric_fingerprint(p, mods):
    a = p % mods[0]
    b = p % mods[1]
    c = p % mods[2]

    e1 = a + b + c
    e2 = a * b + a * c + b * c
    e3 = a * b * c

    return e1, e2, e3


# =================================================================================================
# BUILD GLOBAL PRIME INDEX FOR ONE MODULUS TRIPLE
# =================================================================================================

def build_fingerprint_maps(primes, mods):
    raw_map = defaultdict(list)
    sumprod_map = defaultdict(list)
    sym_map = defaultdict(list)

    for p in primes:
        raw_map[raw_fingerprint(p, mods)].append(p)
        sumprod_map[sum_product_fingerprint(p, mods)].append(p)
        sym_map[symmetric_fingerprint(p, mods)].append(p)

    return raw_map, sumprod_map, sym_map


# =================================================================================================
# FIND FACTOR-PAIR CANDIDATES
# =================================================================================================

def factor_pair_candidates(n, primes, mods):
    """
    Find unordered prime pairs (a,b) such that

        (a*b) mod r == n mod r

    for all three moduli.

    The search uses the product fingerprint.

    For each a, the required residue of b is computed modulo each r.
    """

    n_residues = tuple(n % r for r in mods)

    # Index primes by their residue triple.
    residue_index = defaultdict(list)

    for p in primes:
        key = (
            p % mods[0],
            p % mods[1],
            p % mods[2],
        )
        residue_index[key].append(p)

    result = set()

    for a in primes:
        a1 = a % mods[0]
        a2 = a % mods[1]
        a3 = a % mods[2]

        # a is prime and the moduli are > a? Not necessarily.
        # In our range they are much smaller than a, so inverses exist
        # provided a is not divisible by r. Since a > 10000 and r <= 3000,
        # this is equivalent to a % r != 0.
        if any(a1 == 0 for _ in [0]):
            pass

        if a1 == 0 or a2 == 0 or a3 == 0:
            continue

        try:
            b1 = (n_residues[0] * pow(a1, -1, mods[0])) % mods[0]
            b2 = (n_residues[1] * pow(a2, -1, mods[1])) % mods[1]
            b3 = (n_residues[2] * pow(a3, -1, mods[2])) % mods[2]
        except ValueError:
            continue

        required = (b1, b2, b3)

        for b in residue_index.get(required, ()):
            if a > b:
                pair = (b, a)
            else:
                pair = (a, b)

            # Exact verification against the original integer.
            if pair[0] * pair[1] == n:
                result.add(pair)

    return sorted(result)


# =================================================================================================
# SINGLE-ANCHOR ANALYSIS
# =================================================================================================

def analyze_anchor(p, q, primes, mods, fingerprint_cache):
    n = p * q

    raw_map, sumprod_map, sym_map = fingerprint_cache

    raw_key = raw_fingerprint(p, mods)
    sp_key = sum_product_fingerprint(p, mods)
    sym_key = symmetric_fingerprint(p, mods)

    raw_bucket = raw_map[raw_key]
    sp_bucket = sumprod_map[sp_key]
    sym_bucket = sym_map[sym_key]

    pair_candidates = factor_pair_candidates(n, primes, mods)

    return {
        "p": p,
        "q": q,
        "n": n,
        "mods": mods,
        "R": mods[0] * mods[1] * mods[2],

        "raw_bucket": len(raw_bucket),
        "sumprod_bucket": len(sp_bucket),
        "sym_bucket": len(sym_bucket),

        "raw_values": raw_bucket,
        "sumprod_values": sp_bucket,
        "sym_values": sym_bucket,

        "pair_candidates": pair_candidates,
    }


# =================================================================================================
# MAIN
# =================================================================================================

def run():
    rng = random.Random(SEED)

    print("=" * 100)
    print("FULL PRIME-POPULATION / THREE-MODULUS COMPRESSED FINGERPRINT EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors            = {ACTUAL_ANCHORS}")
    print(f"modulus prime range      = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio              = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # ---------------------------------------------------------------------------------------------
    # PRIME POOLS
    # ---------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = build_primes(FACTOR_MIN, FACTOR_MAX)
    modulus_primes = build_modulus_primes(MODULUS_MIN, MODULUS_MAX)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    # ---------------------------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # ---------------------------------------------------------------------------------------------

    anchors = build_actual_anchors(
        factor_primes,
        ACTUAL_ANCHORS,
        rng,
    )

    print(f"actual anchors             = {len(anchors):,}")
    print()

    # ---------------------------------------------------------------------------------------------
    # CHOOSE r1,r2,r3 FOR EVERY ANCHOR
    # ---------------------------------------------------------------------------------------------

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE PRIME TRIPLES")
    print("=" * 100)

    selected = []

    for i, (p, q) in enumerate(anchors, 1):
        n = p * q

        mods = choose_max_product_triple(
            n,
            modulus_primes,
            CLOSE_RATIO,
        )

        if mods is None:
            print(
                f"WARNING anchor {i}: no suitable modulus triple "
                f"for n={n:,}"
            )
            continue

        selected.append((p, q, n, mods))

        if i % PROGRESS_EVERY == 0:
            print(f"anchor {i:3d}/{len(anchors)}")

    print()
    print(f"usable anchors = {len(selected):,}")

    # ---------------------------------------------------------------------------------------------
    # CACHE FINGERPRINT MAPS
    #
    # Many anchors reuse the same modulus triple. Build maps once per distinct triple.
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING GLOBAL FINGERPRINT MAPS")
    print("=" * 100)

    unique_mod_sets = sorted({
        mods
        for _, _, _, mods in selected
    })

    print(f"unique modulus triples = {len(unique_mod_sets):,}")

    maps_by_mods = {}

    for i, mods in enumerate(unique_mod_sets, 1):
        maps_by_mods[mods] = build_fingerprint_maps(
            factor_primes,
            mods,
        )

        if i % PROGRESS_EVERY == 0 or i == len(unique_mod_sets):
            print(f"triple {i:4d}/{len(unique_mod_sets)}")

    # ---------------------------------------------------------------------------------------------
    # ANALYZE ACTUAL ANCHORS
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANALYZING ACTUAL ANCHORS AGAINST ENTIRE PRIME POOL")
    print("=" * 100)

    results = []

    for i, (p, q, n, mods) in enumerate(selected, 1):

        maps = maps_by_mods[mods]

        result = analyze_anchor(
            p,
            q,
            factor_primes,
            mods,
            maps,
        )

        results.append(result)

        if i % PROGRESS_EVERY == 0:
            print(f"anchor {i:3d}/{len(selected)}")

    # ---------------------------------------------------------------------------------------------
    # SANITY
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SANITY CHECK")
    print("=" * 100)

    identity_failures = 0
    pair_failures = 0

    for r in results:
        p, q, n = r["p"], r["q"], r["n"]
        r1, r2, r3 = r["mods"]

        if n % r1 != (p * q) % r1:
            identity_failures += 1

        if n % r2 != (p * q) % r2:
            identity_failures += 1

        if n % r3 != (p * q) % r3:
            identity_failures += 1

        if (p, q) not in r["pair_candidates"]:
            pair_failures += 1

    print(f"modular identity failures = {identity_failures}")
    print(f"actual pair missing        = {pair_failures}")

    # ---------------------------------------------------------------------------------------------
    # FULL PRIME-POPULATION FINGERPRINT RESULTS
    # ---------------------------------------------------------------------------------------------

    raw_unique = 0
    raw_ambiguous = 0

    sp_unique = 0
    sp_ambiguous = 0

    sym_unique = 0
    sym_ambiguous = 0

    pair_unique = 0
    pair_ambiguous = 0
    pair_zero = 0

    raw_collision_prime_count = 0
    sp_collision_prime_count = 0
    sym_collision_prime_count = 0

    raw_max = 0
    sp_max = 0
    sym_max = 0

    for r in results:
        rb = r["raw_bucket"]
        sb = r["sumprod_bucket"]
        yb = r["sym_bucket"]

        raw_max = max(raw_max, rb)
        sp_max = max(sp_max, sb)
        sym_max = max(sym_max, yb)

        if rb == 1:
            raw_unique += 1
        else:
            raw_ambiguous += 1
            raw_collision_prime_count += rb

        if sb == 1:
            sp_unique += 1
        else:
            sp_ambiguous += 1
            sp_collision_prime_count += sb

        if yb == 1:
            sym_unique += 1
        else:
            sym_ambiguous += 1
            sym_collision_prime_count += yb

        pc = len(r["pair_candidates"])

        if pc == 0:
            pair_zero += 1
        elif pc == 1:
            pair_unique += 1
        else:
            pair_ambiguous += 1

    # ---------------------------------------------------------------------------------------------
    # SUMMARY
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FULL PRIME-POPULATION FINGERPRINT UNIQUENESS")
    print("=" * 100)

    print(
        f"{'FINGERPRINT':30s}"
        f"{'UNIQUE':>12s}"
        f"{'AMBIGUOUS':>14s}"
        f"{'MAX BUCKET':>14s}"
    )
    print("-" * 75)

    print(
        f"{'RAW (a1,a2,a3)':30s}"
        f"{raw_unique:12d}"
        f"{raw_ambiguous:14d}"
        f"{raw_max:14d}"
    )

    print(
        f"{'SUM+PRODUCT (e1,e3)':30s}"
        f"{sp_unique:12d}"
        f"{sp_ambiguous:14d}"
        f"{sp_max:14d}"
    )

    print(
        f"{'FULL SYMMETRIC (e1,e2,e3)':30s}"
        f"{sym_unique:12d}"
        f"{sym_ambiguous:14d}"
        f"{sym_max:14d}"
    )

    # ---------------------------------------------------------------------------------------------
    # FACTOR PAIR SUMMARY
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FACTOR-PAIR RECOVERY")
    print("=" * 100)

    print(f"unique factor pair       = {pair_unique:3d}/{len(results)}")
    print(f"ambiguous factor pair    = {pair_ambiguous:3d}/{len(results)}")
    print(f"zero compatible pairs    = {pair_zero:3d}/{len(results)}")

    pair_counts = Counter(
        len(r["pair_candidates"])
        for r in results
    )

    print()
    print("candidate distribution:")

    for count in sorted(pair_counts):
        print(
            f"    {count:2d} candidates : "
            f"{pair_counts[count]:3d} anchors"
        )

    # ---------------------------------------------------------------------------------------------
    # COLLISION RATIO
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("COMPRESSED-FINGERPRINT COLLISION LOAD")
    print("=" * 100)

    print(f"raw collision prime population       = {raw_collision_prime_count:,}")
    print(f"sum+product collision population     = {sp_collision_prime_count:,}")
    print(f"symmetric collision population       = {sym_collision_prime_count:,}")

    # ---------------------------------------------------------------------------------------------
    # EXAMPLES WHERE SUM+PRODUCT COLLIDES
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUM+PRODUCT COLLISION EXAMPLES")
    print("=" * 100)

    shown = 0

    for r in results:
        if r["sumprod_bucket"] <= 1:
            continue

        print()
        print(
            f"anchor p={r['p']:,} q={r['q']:,} "
            f"mods={r['mods']}"
        )

        print(
            f"  SUM+PRODUCT = {sum_product_fingerprint(r['p'], r['mods'])}"
        )

        print(
            f"  collision count = {r['sumprod_bucket']}"
        )

        values = r["sumprod_values"]

        for value in values[:20]:
            mark = "*" if value == r["p"] else " "
            print(f"  {mark} p={value:,}")

        if len(values) > 20:
            print(
                f"  ... {len(values) - 20} more"
            )

        shown += 1

        if shown >= 20:
            break

    if shown == 0:
        print("none")

    # ---------------------------------------------------------------------------------------------
    # ANCHOR TABLE
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        f"{'ID':>3s} "
        f"{'p':>7s} "
        f"{'q':>7s} "
        f"{'r1':>5s} "
        f"{'r2':>5s} "
        f"{'r3':>5s} "
        f"{'R/n':>12s} "
        f"{'RAW':>5s} "
        f"{'S+P':>5s} "
        f"{'SYM':>5s} "
        f"{'PAIR':>5s}"
    )
    print("-" * 90)

    for i, r in enumerate(results, 1):
        R = r["R"]
        ratio = R / r["n"]

        print(
            f"{i:3d} "
            f"{r['p']:7,d} "
            f"{r['q']:7,d} "
            f"{r['mods'][0]:5d} "
            f"{r['mods'][1]:5d} "
            f"{r['mods'][2]:5d} "
            f"{ratio:12.9f} "
            f"{r['raw_bucket']:5d} "
            f"{r['sumprod_bucket']:5d} "
            f"{r['sym_bucket']:5d} "
            f"{len(r['pair_candidates']):5d}"
        )

    # ---------------------------------------------------------------------------------------------
    # GLOBAL COLLISION TEST ACROSS ALL DISTINCT MODULUS TRIPLES
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("GLOBAL TEST: DOES SUM+PRODUCT EVER COLLIDE?")
    print("=" * 100)

    global_sp_collision_triples = 0
    global_sp_collision_pairs = 0
    global_max_sp_bucket = 0

    for mods in unique_mod_sets:
        _, sp_map, _ = maps_by_mods[mods]

        local_max = max(
            (len(v) for v in sp_map.values()),
            default=0
        )

        collisions = sum(
            1
            for v in sp_map.values()
            if len(v) > 1
        )

        collision_members = sum(
            len(v)
            for v in sp_map.values()
            if len(v) > 1
        )

        if collisions:
            global_sp_collision_triples += 1
            global_sp_collision_pairs += collision_members

        global_max_sp_bucket = max(
            global_max_sp_bucket,
            local_max
        )

    print(
        f"modulus triples containing SUM+PRODUCT collisions = "
        f"{global_sp_collision_triples:,}/{len(unique_mod_sets):,}"
    )

    print(
        f"total prime members in collision buckets          = "
        f"{global_sp_collision_pairs:,}"
    )

    print(
        f"largest SUM+PRODUCT collision bucket              = "
        f"{global_max_sp_bucket:,}"
    )

    # ---------------------------------------------------------------------------------------------
    # IMPORTANT COMPARISON
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("KEY COMPARISON")
    print("=" * 100)

    print(
        "For the actual anchors:"
    )
    print()

    print(
        f"RAW unique              = "
        f"{raw_unique}/{len(results)} "
        f"({100.0 * raw_unique / len(results):.2f}%)"
    )

    print(
        f"SUM+PRODUCT unique      = "
        f"{sp_unique}/{len(results)} "
        f"({100.0 * sp_unique / len(results):.2f}%)"
    )

    print(
        f"FULL SYMMETRIC unique   = "
        f"{sym_unique}/{len(results)} "
        f"({100.0 * sym_unique / len(results):.2f}%)"
    )

    print(
        f"PAIR unique             = "
        f"{pair_unique}/{len(results)} "
        f"({100.0 * pair_unique / len(results):.2f}%)"
    )

    # ---------------------------------------------------------------------------------------------
    # INTERPRETATION
    # ---------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
This experiment tests the proposed compressed p-fingerprint against the
ENTIRE prime population in the factor interval.

For each actual anchor we first choose its maximum-product close triple

    r1 < r2 < r3

with

    r1*r2*r3 < n.

For every prime p in the complete factor interval we calculate

    a1 = p mod r1
    a2 = p mod r2
    a3 = p mod r3.

The three tested fingerprints are:

    RAW:
        (a1,a2,a3)

    SUM+PRODUCT:
        (a1+a2+a3, a1*a2*a3)

    FULL SYMMETRIC:
        (e1,e2,e3)

The important difference from the previous experiment is that the bucket
is now measured against all 8,363 factor primes.

Therefore:

    SUM+PRODUCT bucket = 1

means that no other prime in the entire tested factor range shares the
same two-number fingerprint.

That is much stronger than observing uniqueness among only the 300 actual
factors.

The final factor-pair test remains separate.

A unique p fingerprint does NOT automatically imply unique factorization,
because q must also satisfy the product constraints.

The strongest observed result would therefore be:

    SUM+PRODUCT bucket = 1
    AND
    factor-pair candidates = 1.

If this happens repeatedly for the selected maximum-product triples, the
two quantities

    sum(p mod r_i)
    product(p mod r_i)

are providing a highly discriminating representation of p in this
finite population.

However, this still remains a finite-range information experiment.

It does NOT establish a universal factoring method.

The next mathematically important question after this is whether the same
SUM+PRODUCT uniqueness survives when the r1,r2,r3 selection is changed,
especially when R/n is deliberately matched against alternative triples.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
