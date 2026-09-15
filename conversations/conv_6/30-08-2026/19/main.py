#!/usr/bin/env python3

import math
import random
from dataclasses import dataclass


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

ACTUAL_ANCHORS = 300

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_EVERY = 25


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass(frozen=True)
class Triple:
    r1: int
    r2: int
    r3: int

    @property
    def R(self) -> int:
        return self.r1 * self.r2 * self.r3


# =============================================================================
# SAFE SIEVE
# =============================================================================

def sieve(limit: int) -> list[int]:
    """
    Return all primes <= limit.

    Deliberately uses a straightforward boolean bytearray implementation
    instead of extended-slice assignment so there is no slice-length
    off-by-one hazard.
    """
    if limit < 2:
        return []

    composite = bytearray(limit + 1)

    composite[0] = 1
    composite[1] = 1

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if composite[p]:
            continue

        start = p * p

        for multiple in range(start, limit + 1, p):
            composite[multiple] = 1

    return [
        n
        for n in range(2, limit + 1)
        if not composite[n]
    ]


# =============================================================================
# BINARY SEARCH
# =============================================================================

def bisect_right(a: list[int], x: int) -> int:
    lo = 0
    hi = len(a)

    while lo < hi:
        mid = (lo + hi) // 2

        if a[mid] <= x:
            lo = mid + 1
        else:
            hi = mid

    return lo


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def generate_anchors(
    primes: list[int],
    count: int,
    rng: random.Random,
) -> list[Anchor]:

    result: list[Anchor] = []
    seen: set[tuple[int, int]] = set()

    while len(result) < count:

        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in seen:
            continue

        seen.add(pair)
        result.append(Anchor(p, q))

    return result


# =============================================================================
# MAXIMUM-PRODUCT CLOSE TRIPLE
# =============================================================================

def choose_triple(
    n: int,
    primes: list[int],
    close_ratio: float,
) -> Triple:

    best: Triple | None = None
    best_R = -1

    prime_count = len(primes)

    for i in range(prime_count - 2):

        r1 = primes[i]

        # Optional closeness restriction relative to r1.
        max_allowed_r3 = int(r1 * (1.0 + close_ratio))

        for j in range(i + 1, prime_count - 1):

            r2 = primes[j]

            r12 = r1 * r2

            if r12 >= n:
                break

            # Need:
            #
            #     r1*r2*r3 < n
            #
            # therefore:
            #
            #     r3 <= (n-1)/(r1*r2)
            #
            max_r3 = (n - 1) // r12

            if max_r3 <= r2:
                continue

            max_r3 = min(max_r3, max_allowed_r3)

            k = bisect_right(primes, max_r3) - 1

            if k <= j:
                continue

            r3 = primes[k]

            R = r1 * r2 * r3

            if R >= n:
                continue

            if R > best_R:
                best_R = R
                best = Triple(r1, r2, r3)

    if best is None:
        raise RuntimeError(
            f"No suitable modulus triple found for n={n:,}"
        )

    return best


# =============================================================================
# BOUNDED VALUES WITH A GIVEN RESIDUE
# =============================================================================

def bounded_values_with_residue(
    residue: int,
    modulus: int,
) -> list[int]:

    first = residue

    while first < FACTOR_MIN:
        first += modulus

    if first > FACTOR_MAX:
        return []

    count = ((FACTOR_MAX - first) // modulus) + 1

    return [
        first + i * modulus
        for i in range(count)
    ]


# =============================================================================
# TRIAL PRIMALITY
# =============================================================================

def is_prime(n: int) -> bool:

    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    if n % 3 == 0:
        return n == 3

    limit = math.isqrt(n)

    d = 5

    while d <= limit:
        if n % d == 0:
            return False

        if d + 2 <= limit and n % (d + 2) == 0:
            return False

        d += 6

    return True


# =============================================================================
# DIRECT ANALYSIS
# =============================================================================

def analyze_anchor(
    anchor: Anchor,
    triple: Triple,
) -> dict:

    n = anchor.n

    r1 = triple.r1
    r2 = triple.r2
    r3 = triple.r3

    # -------------------------------------------------------------------------
    # BASELINE
    # -------------------------------------------------------------------------

    baseline = FACTOR_MAX - FACTOR_MIN + 1

    # -------------------------------------------------------------------------
    # STAGE 1
    #
    # Represent every bounded p by p mod r1.
    #
    # For each nonzero residue a:
    #
    #     q == n * a^-1 (mod r1)
    #
    # -------------------------------------------------------------------------

    step1: list[tuple[int, int]] = []

    nr1 = n % r1

    for a in range(1, r1):

        b = (nr1 * pow(a, -1, r1)) % r1

        values = bounded_values_with_residue(
            a,
            r1,
        )

        for p in values:
            step1.append((p, b))

    # -------------------------------------------------------------------------
    # STAGE 2
    #
    # Exact divisibility.
    #
    # This is NOT the intended modular speedup; it is included as a baseline
    # comparison.
    # -------------------------------------------------------------------------

    divisible: list[tuple[int, int]] = []

    for p, _ in step1:

        if n % p != 0:
            continue

        q = n // p

        if FACTOR_MIN <= q <= FACTOR_MAX:
            divisible.append((p, q))

    # -------------------------------------------------------------------------
    # STAGE 3
    #
    # Verify the second modulus.
    # -------------------------------------------------------------------------

    mod2: list[tuple[int, int]] = []

    nr2 = n % r2

    for p, q in divisible:

        if (p * q) % r2 == nr2:
            mod2.append((p, q))

    # -------------------------------------------------------------------------
    # STAGE 4
    #
    # Verify the third modulus.
    # -------------------------------------------------------------------------

    mod3: list[tuple[int, int]] = []

    nr3 = n % r3

    for p, q in mod2:

        if (p * q) % r3 == nr3:
            mod3.append((p, q))

    # -------------------------------------------------------------------------
    # DIRECT TWO-MODULUS CRT RECONSTRUCTION
    #
    # We reconstruct p itself from:
    #
    #     p == a1 (mod r1)
    #     p == a2 (mod r2)
    #
    # yielding
    #
    #     p == x (mod r1*r2).
    #
    # -------------------------------------------------------------------------

    R12 = r1 * r2

    inv_r1_r2 = pow(r1, -1, r2)

    crt12_set: set[int] = set()

    # p mod r1 can be any nonzero residue
    # p mod r2 can be any nonzero residue.
    #
    # We do NOT test divisibility by n here.

    for a1 in range(1, r1):

        for a2 in range(1, r2):

            t = ((a2 - a1) * inv_r1_r2) % r2

            x = a1 + r1 * t

            # x is in [0,R12).
            #
            # Find all representatives in the bounded factor interval.

            if x < FACTOR_MIN:

                x += (
                    ((FACTOR_MIN - x + R12 - 1) // R12)
                    * R12
                )

            while x <= FACTOR_MAX:

                if x >= FACTOR_MIN:
                    crt12_set.add(x)

                x += R12

    crt12 = sorted(crt12_set)

    # -------------------------------------------------------------------------
    # THIRD MODULUS RESIDUE FILTER
    #
    # Now we filter CRT candidates by whether there exists a compatible
    # bounded q residue modulo r3.
    #
    # IMPORTANT:
    # this stage still does not use n % p.
    #
    # Given p:
    #
    #     q == n * p^-1 (mod r3)
    #
    # We only ask whether that q residue has a bounded representative.
    # -------------------------------------------------------------------------

    crt3: list[int] = []

    for p in crt12:

        if p % r3 == 0:
            if nr3 != 0:
                continue

            # If n == 0 mod r3, p may be zero mod r3.
            # The general q-residue treatment is handled below.

        else:

            q_residue = (
                nr3 * pow(p, -1, r3)
            ) % r3

            first_q = q_residue

            while first_q < FACTOR_MIN:
                first_q += r3

            if first_q > FACTOR_MAX:
                continue

            crt3.append(p)

            continue

        # Special case p == 0 mod r3.
        #
        # We need q such that:
        #
        #     p*q == n == 0 (mod r3)
        #
        # so any q is compatible.
        #

        crt3.append(p)

    # -------------------------------------------------------------------------
    # THIRD MODULUS + EXACT n TEST
    #
    # This is now deliberately delayed until after CRT candidate collapse.
    # -------------------------------------------------------------------------

    exact_candidates: list[tuple[int, int]] = []

    for p in crt3:

        if n % p != 0:
            continue

        q = n // p

        if not (FACTOR_MIN <= q <= FACTOR_MAX):
            continue

        if q % r3 != (
            nr3 * pow(p % r3, -1, r3) % r3
            if p % r3
            else 0
        ):
            continue

        exact_candidates.append((p, q))

    # -------------------------------------------------------------------------
    # PRIME PAIR RECOVERY
    # -------------------------------------------------------------------------

    prime_pairs = []

    for p, q in exact_candidates:

        if not is_prime(p):
            continue

        if not is_prime(q):
            continue

        pair = tuple(sorted((p, q)))

        if pair not in prime_pairs:
            prime_pairs.append(pair)

    prime_pairs.sort()

    return {
        "baseline": baseline,
        "step1": len(step1),
        "divisible": len(divisible),
        "mod2": len(mod2),
        "mod3": len(mod3),
        "crt12": len(crt12),
        "crt3": len(crt3),
        "exact": len(exact_candidates),
        "prime_pairs": prime_pairs,
    }


# =============================================================================
# MAIN
# =============================================================================

def run():

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-MODULUS DIRECT CRT FACTOR-RESIDUE RECONSTRUCTION EXPERIMENT")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(
        f"factor range         = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"anchors              = {ACTUAL_ANCHORS}")
    print(
        f"modulus prime range  = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )
    print(f"close ratio           = {CLOSE_RATIO:.0%}")
    print(f"seed                  = {SEED:,}")

    # -------------------------------------------------------------------------
    # PRIME POOLS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    all_primes = sieve(FACTOR_MAX)

    factor_primes = [
        p
        for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in sieve(MODULUS_MAX)
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(
        f"factor primes        = {len(factor_primes):,}"
    )

    print(
        f"modulus primes       = {len(modulus_primes):,}"
    )

    # -------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # -------------------------------------------------------------------------

    anchors = generate_anchors(
        factor_primes,
        ACTUAL_ANCHORS,
        rng,
    )

    print(
        f"actual anchors       = {len(anchors)}"
    )

    # -------------------------------------------------------------------------
    # SELECT MODULUS TRIPLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE PRIME TRIPLES")
    print("=" * 100)

    triples: list[Triple] = []

    for i, anchor in enumerate(anchors, 1):

        triple = choose_triple(
            anchor.n,
            modulus_primes,
            CLOSE_RATIO,
        )

        triples.append(triple)

        if i % PROGRESS_EVERY == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    # -------------------------------------------------------------------------
    # ANALYZE
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("DIRECT CRT ANALYSIS")
    print("=" * 100)

    results = []

    for i, (anchor, triple) in enumerate(
        zip(anchors, triples),
        1,
    ):

        result = analyze_anchor(
            anchor,
            triple,
        )

        results.append(result)

        if i % PROGRESS_EVERY == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    # -------------------------------------------------------------------------
    # MEAN
    # -------------------------------------------------------------------------

    def mean(key: str) -> float:
        return (
            sum(result[key] for result in results)
            / len(results)
        )

    # -------------------------------------------------------------------------
    # SEARCH SPACE
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SEARCH SPACE")
    print("=" * 100)

    print(
        f"baseline p values            = "
        f"{mean('baseline'):,.2f}"
    )

    print(
        f"after r1 representation      = "
        f"{mean('step1'):,.2f}"
    )

    print(
        f"after exact divisibility     = "
        f"{mean('divisible'):,.2f}"
    )

    print(
        f"after r2                    = "
        f"{mean('mod2'):,.2f}"
    )

    print(
        f"after r3                    = "
        f"{mean('mod3'):,.2f}"
    )

    print(
        f"CRT12 bounded p candidates  = "
        f"{mean('crt12'):,.2f}"
    )

    print(
        f"CRT12 + r3 candidates       = "
        f"{mean('crt3'):,.2f}"
    )

    print(
        f"exact n=pq candidates       = "
        f"{mean('exact'):,.2f}"
    )

    # -------------------------------------------------------------------------
    # REDUCTIONS
    # -------------------------------------------------------------------------

    baseline = mean("baseline")

    print()
    print("=" * 100)
    print("REDUCTION RATIOS")
    print("=" * 100)

    for label, key in [
        ("r1 representation", "step1"),
        ("CRT12", "crt12"),
        ("CRT12 + r3", "crt3"),
        ("exact", "exact"),
    ]:

        value = mean(key)

        reduction = (
            100.0 * (1.0 - value / baseline)
        )

        print(
            f"{label:20s} = "
            f"{reduction:10.6f}%"
        )

    # -------------------------------------------------------------------------
    # RECOVERY
    # -------------------------------------------------------------------------

    recovered = 0
    unique_crt = 0
    unique_exact = 0

    for anchor, result in zip(anchors, results):

        target = tuple(
            sorted((anchor.p, anchor.q))
        )

        if target in result["prime_pairs"]:
            recovered += 1

        if len(result["crt3"]) == 1:
            unique_crt += 1

        if len(result["prime_pairs"]) == 1:
            unique_exact += 1

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"correct prime pair recovered = "
        f"{recovered}/{len(anchors)}"
    )

    print(
        f"unique CRT3 candidate        = "
        f"{unique_crt}/{len(anchors)}"
    )

    print(
        f"unique final prime pair      = "
        f"{unique_exact}/{len(anchors)}"
    )

    # -------------------------------------------------------------------------
    # CRT12 DISTRIBUTION
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CRT12 CANDIDATE DISTRIBUTION")
    print("=" * 100)

    dist = {}

    for result in results:

        value = result["crt12"]

        dist[value] = dist.get(value, 0) + 1

    for value in sorted(dist):

        print(
            f"CRT12 candidates = "
            f"{value:6d} "
            f"anchors = {dist[value]:4d}"
        )

    # -------------------------------------------------------------------------
    # CRT3 DISTRIBUTION
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CRT3 CANDIDATE DISTRIBUTION")
    print("=" * 100)

    dist3 = {}

    for result in results:

        value = result["crt3"]

        dist3[value] = dist3.get(value, 0) + 1

    for value in sorted(dist3):

        print(
            f"CRT3 candidates = "
            f"{value:6d} "
            f"anchors = {dist3[value]:4d}"
        )

    # -------------------------------------------------------------------------
    # ANCHOR TABLE
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q      "
        "r1     r2     r3       R/n     "
        "BASE   CRT12   CRT3   EXACT"
    )

    print("-" * 100)

    for i, (
        anchor,
        triple,
        result,
    ) in enumerate(
        zip(anchors, triples, results),
        1,
    ):

        print(
            f"{i:3d} "
            f"{anchor.p:8,d} "
            f"{anchor.q:8,d} "
            f"{triple.r1:6,d} "
            f"{triple.r2:6,d} "
            f"{triple.r3:6,d} "
            f"{triple.R / anchor.n:9.7f} "
            f"{result['baseline']:6,d} "
            f"{result['crt12']:7,d} "
            f"{result['crt3']:6,d} "
            f"{result['exact']:6,d}"
        )

    # -------------------------------------------------------------------------
    # BEST RESULTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BEST CRT COLLAPSES")
    print("=" * 100)

    ranked = sorted(
        zip(anchors, triples, results),
        key=lambda item: (
            item[2]["crt3"],
            item[2]["crt12"],
        ),
    )

    for anchor, triple, result in ranked[:20]:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"mods=({triple.r1},{triple.r2},{triple.r3}) "
            f"R/n={triple.R / anchor.n:.10f} "
            f"CRT12={result['crt12']} "
            f"CRT3={result['crt3']} "
            f"exact={result['exact']} "
            f"pairs={result['prime_pairs']}"
        )

    # -------------------------------------------------------------------------
    # VERIFY ACTUAL FACTOR RESIDUES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL RESIDUE CHECK")
    print("=" * 100)

    failures = 0

    for anchor, triple in zip(anchors, triples):

        n = anchor.n

        for r in (
            triple.r1,
            triple.r2,
            triple.r3,
        ):

            if (
                anchor.p * anchor.q
            ) % r != n % r:

                failures += 1

    print(
        f"residue identity failures = {failures}"
    )

    # -------------------------------------------------------------------------
    # THEORY
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
This experiment changes the attack surface.

The previous experiment searched over:

    s = p + q

and used:

    D = s^2 - 4n.

This experiment instead asks:

    Can p itself be reconstructed from modular coordinates?

For two primes r1 and r2:

    p == a1 (mod r1)
    p == a2 (mod r2)

CRT gives:

    p == x (mod r1*r2).

Because p is restricted to:

    10,000 <= p <= 100,000,

the number of bounded representatives of x can become very small
once:

    r1*r2 >> 100,000.

The third modulus is then used as another coordinate:

    p == a3 (mod r3).

The experiment therefore measures a direct pipeline:

    bounded p population
             |
             v
       two-modulus CRT
             |
             v
       third-modulus filter
             |
             v
       exact n = p*q
             |
             v
       prime-pair recovery.

The critical result is CRT12.

If CRT12 is already tiny before checking whether p divides n,
then the modular information itself is collapsing the bounded
factor search.

That is substantially different from merely observing that a
known factor has a unique fingerprint.

The experiment also keeps the third modulus separate so we can
see whether it contributes genuine additional reduction.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()