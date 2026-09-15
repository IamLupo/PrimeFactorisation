#!/usr/bin/env python3

"""
============================================================================================
TWO-MODULUS K-STABILITY RADIUS / QUOTIENT-CELL BOUNDARY EXPERIMENT
============================================================================================

Goal
----

We observed that nearby semiprimes frequently preserve

    k_x = k
    l_x = l
    K_x = k_x*l_x = k*l

even though

    n+x != n.

This experiment measures exactly how large that stability region is.

For

    p = a + k*r1
    q = b + l*r2

the quotient-cell boundaries are:

    k*r1 <= p < (k+1)*r1
    l*r2 <= q < (l+1)*r2

Therefore the distances from the true factors to their cell boundaries
are directly observable once the anchor is known.

The experiment compares:

    ACTUAL FACTOR TRAJECTORY
        n+x = p_x*q_x

against

    GEOMETRIC CELL PREDICTION

and asks:

    Does K remain constant exactly until a factor crosses one of the
    r1/r2 quotient-cell boundaries?

This is important because it distinguishes:

    ordinary local coordinate stability

from:

    a deeper arithmetic invariant.

The experiment also tests whether K-stability can be estimated from
local information around the original point.

No candidate factor search is performed for n+x.

Nearby semiprimes are deliberately constructed from prime pairs and
then verified exactly.

============================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

# Search perturbations.
#
# This is deliberately larger than the previous experiment because we
# want to find the first K-boundary rather than only a few random hits.
X_RADIUS = 20_000

# Number of explicitly sampled x values on each side when printing
# trajectory examples.
TRAJECTORY_STEP = 250

SEED = 1_511_464_998

PRINT_EXAMPLES = 20


# ==========================================================================================
# DATA STRUCTURES
# ==========================================================================================

@dataclass
class Anchor:
    p: int
    q: int
    n: int


@dataclass
class ModPair:
    r1: int
    r2: int


@dataclass
class Point:
    x: int
    p: int
    q: int
    k: int
    l: int
    K: int


# ==========================================================================================
# PRIME GENERATION
# ==========================================================================================

def sieve(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, int(limit ** 0.5) + 1):
        if is_prime[i]:
            start = i * i
            is_prime[start::i] = b"\x00" * (
                ((limit - start) // i) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


# ==========================================================================================
# ANCHORS
# ==========================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    p_min: int,
    p_max: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    candidates = [
        p for p in primes
        if p_min <= p <= p_max
    ]

    if len(candidates) < 2:
        raise RuntimeError("Not enough factor primes.")

    anchors: list[Anchor] = []
    seen: set[tuple[int, int]] = set()

    # We deliberately construct semiprimes from primes.
    #
    # Using a deterministic randomized process makes the experiment
    # reproducible while avoiding accidental duplicate anchors.

    attempts = 0
    max_attempts = count * 500

    while len(anchors) < count and attempts < max_attempts:
        attempts += 1

        p = rng.choice(candidates)
        q = rng.choice(candidates)

        if p > q:
            p, q = q, p

        if p == q:
            continue

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)

        anchors.append(
            Anchor(
                p=p,
                q=q,
                n=p * q,
            )
        )

    if len(anchors) < count:
        raise RuntimeError(
            f"Could only build {len(anchors)} anchors."
        )

    return anchors


# ==========================================================================================
# CLOSE MODULUS PAIRS
# ==========================================================================================

def build_close_modulus_pairs(
    modulus_primes: list[int],
    close_ratio: float,
) -> list[ModPair]:

    pairs: list[ModPair] = []

    for i, r1 in enumerate(modulus_primes):
        for r2 in modulus_primes[i + 1:]:
            ratio = abs(r2 - r1) / min(r1, r2)

            if ratio <= close_ratio:
                pairs.append(
                    ModPair(r1=r1, r2=r2)
                )

    return pairs


# ==========================================================================================
# FACTOR / CELL COORDINATES
# ==========================================================================================

def coordinates(value: int, modulus: int) -> tuple[int, int]:
    """
    value = residue + quotient*modulus

    Return:

        quotient, residue
    """
    q, r = divmod(value, modulus)
    return q, r


# ==========================================================================================
# GEOMETRIC DISTANCES
# ==========================================================================================

def cell_distances(value: int, modulus: int) -> tuple[int, int]:
    """
    For

        k = floor(value/modulus)

    return the distance from value to:

        lower boundary
        upper boundary

    in integer units.

    Example:

        value = k*r + a

    then:

        distance below = a
        distance above = r-1-a
    """
    k, a = coordinates(value, modulus)

    lower = a
    upper = modulus - 1 - a

    return lower, upper


# ==========================================================================================
# FIND NEARBY PRIME PAIRS
# ==========================================================================================

def nearby_factorizations(
    target: int,
    factor_primes_set: set[int],
    max_count: int = 1,
) -> list[tuple[int, int]]:

    """
    Find exact prime factorization(s) of target by searching the smaller
    factor side.

    This is used ONLY for deliberately generated nearby semiprimes.

    For performance we first try sqrt(target) downwards using the
    prime table.

    Since nearby targets are generated from known prime pairs, this is
    not a general-purpose factoring algorithm.
    """

    root = math.isqrt(target)

    result: list[tuple[int, int]] = []

    # Search downward from sqrt(target).

    for p in range(root, FACTOR_MIN - 1, -1):
        if p not in factor_primes_set:
            continue

        if target % p != 0:
            continue

        q = target // p

        if q not in factor_primes_set:
            continue

        if p > q:
            p, q = q, p

        result.append((p, q))

        if len(result) >= max_count:
            break

    return result


# ==========================================================================================
# CONSTRUCT NEARBY SEMIPRIMES
# ==========================================================================================

def build_nearby_points(
    anchor: Anchor,
    factor_primes: list[int],
    factor_prime_set: set[int],
    x_radius: int,
    rng: random.Random,
) -> list[Point]:

    """
    Instead of trying to factor every integer n+x, construct nearby
    semiprimes directly.

    We perturb p and q independently by small prime differences:

        p_x = p + dp
        q_x = q + dq

    and therefore:

        x = p_x*q_x - p*q.

    Only points with |x| <= X_RADIUS are retained.
    """

    p0 = anchor.p
    q0 = anchor.q

    k0_ref = None
    l0_ref = None
    K0_ref = None

    # The caller will attach coordinates later.

    generated: dict[int, Point] = {}

    # Build several perturbations around each factor.

    local_primes_p = [
        p for p in factor_primes
        if abs(p - p0) <= max(500, int(math.sqrt(x_radius) * 50))
    ]

    local_primes_q = [
        q for q in factor_primes
        if abs(q - q0) <= max(500, int(math.sqrt(x_radius) * 50))
    ]

    if not local_primes_p:
        local_primes_p = [p0]

    if not local_primes_q:
        local_primes_q = [q0]

    # Add the original point.

    generated[0] = Point(
        x=0,
        p=p0,
        q=q0,
        k=0,
        l=0,
        K=0,
    )

    # Random combinations.

    attempts = 0
    max_attempts = 5000

    while attempts < max_attempts:
        attempts += 1

        px = rng.choice(local_primes_p)
        qx = rng.choice(local_primes_q)

        x = px * qx - anchor.n

        if x == 0:
            continue

        if abs(x) > x_radius:
            continue

        kx, _ = coordinates(px, 1)
        lx, _ = coordinates(qx, 1)

        # Placeholder coordinates.
        generated[x] = Point(
            x=x,
            p=px,
            q=qx,
            k=kx,
            l=lx,
            K=kx * lx,
        )

    return list(generated.values())


# ==========================================================================================
# MAIN TRAJECTORY ANALYSIS
# ==========================================================================================

def analyze_anchor(
    anchor: Anchor,
    mod_pair: ModPair,
    factor_primes: list[int],
    factor_prime_set: set[int],
    rng: random.Random,
) -> dict:

    p = anchor.p
    q = anchor.q

    r1 = mod_pair.r1
    r2 = mod_pair.r2

    k, a = coordinates(p, r1)
    l, b = coordinates(q, r2)

    K = k * l

    # ------------------------------------------------------------------
    # Geometric boundary distances.
    # ------------------------------------------------------------------

    p_lower, p_upper = cell_distances(p, r1)
    q_lower, q_upper = cell_distances(q, r2)

    # ------------------------------------------------------------------
    # The original cell boundaries in terms of p and q.
    # ------------------------------------------------------------------

    p_min_cell = k * r1
    p_max_cell = (k + 1) * r1 - 1

    q_min_cell = l * r2
    q_max_cell = (l + 1) * r2 - 1

    # ------------------------------------------------------------------
    # Construct nearby semiprimes.
    # ------------------------------------------------------------------

    nearby = build_nearby_points(
        anchor,
        factor_primes,
        factor_prime_set,
        X_RADIUS,
        rng,
    )

    # Correct quotient coordinates in the selected modulus system.

    points: list[Point] = []

    for point in nearby:
        kx, ax = coordinates(point.p, r1)
        lx, bx = coordinates(point.q, r2)

        points.append(
            Point(
                x=point.x,
                p=point.p,
                q=point.q,
                k=kx,
                l=lx,
                K=kx * lx,
            )
        )

    # Ensure original point exists.

    original = Point(
        x=0,
        p=p,
        q=q,
        k=k,
        l=l,
        K=K,
    )

    points.append(original)

    # Remove duplicates by x.

    unique = {}

    for point in points:
        unique[point.x] = point

    points = sorted(unique.values(), key=lambda z: z.x)

    # ------------------------------------------------------------------
    # Classify points.
    # ------------------------------------------------------------------

    same_k = []
    same_l = []
    same_K = []
    same_full_cell = []

    for point in points:
        if point.k == k:
            same_k.append(point)

        if point.l == l:
            same_l.append(point)

        if point.K == K:
            same_K.append(point)

        if point.k == k and point.l == l:
            same_full_cell.append(point)

    # ------------------------------------------------------------------
    # Find observed K-stability limits.
    # ------------------------------------------------------------------

    negative_same_K = sorted(
        [pt.x for pt in same_K if pt.x < 0],
        reverse=True,
    )

    positive_same_K = sorted(
        [pt.x for pt in same_K if pt.x > 0]
    )

    observed_left_K = (
        negative_same_K[0]
        if negative_same_K
        else None
    )

    observed_right_K = (
        positive_same_K[0]
        if positive_same_K
        else None
    )

    # ------------------------------------------------------------------
    # Find observed full-cell limits.
    # ------------------------------------------------------------------

    negative_full = sorted(
        [pt.x for pt in same_full_cell if pt.x < 0],
        reverse=True,
    )

    positive_full = sorted(
        [pt.x for pt in same_full_cell if pt.x > 0]
    )

    observed_left_cell = (
        negative_full[0]
        if negative_full
        else None
    )

    observed_right_cell = (
        positive_full[0]
        if positive_full
        else None
    )

    return {
        "p": p,
        "q": q,
        "n": anchor.n,
        "r1": r1,
        "r2": r2,
        "k": k,
        "l": l,
        "K": K,
        "a": a,
        "b": b,

        "p_lower": p_lower,
        "p_upper": p_upper,
        "q_lower": q_lower,
        "q_upper": q_upper,

        "p_min_cell": p_min_cell,
        "p_max_cell": p_max_cell,
        "q_min_cell": q_min_cell,
        "q_max_cell": q_max_cell,

        "points": points,

        "same_k_count": len(same_k),
        "same_l_count": len(same_l),
        "same_K_count": len(same_K),
        "same_cell_count": len(same_full_cell),

        "observed_left_K": observed_left_K,
        "observed_right_K": observed_right_K,

        "observed_left_cell": observed_left_cell,
        "observed_right_cell": observed_right_cell,
    }


# ==========================================================================================
# PRINT HELPERS
# ==========================================================================================

def separator(char: str = "=", width: int = 92) -> None:
    print(char * width)


def section(title: str) -> None:
    separator()
    print(title)
    separator()


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    t0 = time.perf_counter()

    print("=" * 92)
    print("TWO-MODULUS K-STABILITY RADIUS / QUOTIENT-CELL BOUNDARY EXPERIMENT")
    print("=" * 92)
    print(f"N anchors                 = {N_ANCHORS:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range             = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.2%}")
    print(f"x radius                  = {X_RADIUS:,}")
    print(f"seed                      = {SEED:,}")
    print()

    # ------------------------------------------------------------------
    # Prime pools.
    # ------------------------------------------------------------------

    section("BUILDING PRIME POOLS")

    factor_primes = sieve(FACTOR_MAX)
    factor_primes = [
        p for p in factor_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = sieve(MODULUS_MAX)
    modulus_primes = [
        p for p in modulus_primes
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    factor_prime_set = set(factor_primes)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    # ------------------------------------------------------------------
    # Anchors.
    # ------------------------------------------------------------------

    section("BUILDING ANCHORS")

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        FACTOR_MIN,
        FACTOR_MAX,
        SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")

    # ------------------------------------------------------------------
    # Close modulus pairs.
    # ------------------------------------------------------------------

    section("BUILDING CLOSE MODULUS PAIRS")

    modulus_pairs = build_close_modulus_pairs(
        modulus_primes,
        CLOSE_RATIO,
    )

    print(f"close modulus pairs       = {len(modulus_pairs):,}")

    if not modulus_pairs:
        raise RuntimeError("No modulus pairs found.")

    rng = random.Random(SEED)

    # ------------------------------------------------------------------
    # Choose exactly one random modulus pair per anchor.
    # ------------------------------------------------------------------

    selected_pairs = [
        rng.choice(modulus_pairs)
        for _ in anchors
    ]

    # ------------------------------------------------------------------
    # Run.
    # ------------------------------------------------------------------

    section("RUNNING K-STABILITY TRAJECTORIES")

    results = []

    for i, (anchor, mod_pair) in enumerate(
        zip(anchors, selected_pairs),
        start=1,
    ):

        if i == 1 or i % 10 == 0 or i == len(anchors):
            print(f"anchor {i:3d}/{len(anchors)}")

        result = analyze_anchor(
            anchor,
            mod_pair,
            factor_primes,
            factor_prime_set,
            rng,
        )

        results.append(result)

    # ------------------------------------------------------------------
    # Aggregate.
    # ------------------------------------------------------------------

    K_stable_points = []
    cell_stable_points = []

    K_stable_fraction = []
    full_cell_fraction = []

    p_distances = []
    q_distances = []

    for r in results:
        pts = r["points"]

        # Exclude x=0 when computing local trajectory fractions.

        nonzero = [pt for pt in pts if pt.x != 0]

        if nonzero:
            K_stable_fraction.append(
                sum(pt.K == r["K"] for pt in nonzero) / len(nonzero)
            )

            full_cell_fraction.append(
                sum(
                    pt.k == r["k"] and pt.l == r["l"]
                    for pt in nonzero
                ) / len(nonzero)
            )

        K_stable_points.append(r["same_K_count"])
        cell_stable_points.append(r["same_cell_count"])

        p_distances.append(
            min(r["p_lower"], r["p_upper"])
        )

        q_distances.append(
            min(r["q_lower"], r["q_upper"])
        )

    # ------------------------------------------------------------------
    # Summary.
    # ------------------------------------------------------------------

    section("SUMMARY")

    print(f"anchors analyzed                = {len(results):,}")

    print()
    print("QUOTIENT COORDINATE STABILITY")
    print(f"average sampled K-stable points = {statistics.mean(K_stable_points):.3f}")
    print(f"average sampled same-cell pts   = {statistics.mean(cell_stable_points):.3f}")

    print()
    print("TRAJECTORY STABILITY")
    print(
        "mean fraction K_x=K             = "
        f"{statistics.mean(K_stable_fraction):.6f}"
    )
    print(
        "mean fraction full cell         = "
        f"{statistics.mean(full_cell_fraction):.6f}"
    )

    print()
    print("GEOMETRIC CELL DISTANCES")
    print(
        "mean min distance p -> boundary = "
        f"{statistics.mean(p_distances):.3f}"
    )
    print(
        "mean min distance q -> boundary = "
        f"{statistics.mean(q_distances):.3f}"
    )

    # ------------------------------------------------------------------
    # Detailed examples.
    # ------------------------------------------------------------------

    section("EXAMPLES")

    for r in results[:PRINT_EXAMPLES]:

        print(
            f"n={r['n']:,} "
            f"p={r['p']:,} "
            f"q={r['q']:,} "
            f"mods=({r['r1']},{r['r2']})"
        )

        print(
            f"    coordinates: "
            f"k={r['k']} l={r['l']} K={r['K']}"
        )

        print(
            f"    residues: "
            f"a={r['a']} b={r['b']}"
        )

        print(
            f"    p-cell: "
            f"[{r['p_min_cell']},{r['p_max_cell']}]"
        )

        print(
            f"    q-cell: "
            f"[{r['q_min_cell']},{r['q_max_cell']}]"
        )

        print(
            f"    p boundary distances: "
            f"lower={r['p_lower']} upper={r['p_upper']}"
        )

        print(
            f"    q boundary distances: "
            f"lower={r['q_lower']} upper={r['q_upper']}"
        )

        print(
            f"    sampled K-stable points = "
            f"{r['same_K_count']}"
        )

        print(
            f"    sampled same-cell points = "
            f"{r['same_cell_count']}"
        )

        if r["observed_left_K"] is not None:
            print(
                f"    closest sampled K point left  = "
                f"x={r['observed_left_K']}"
            )

        if r["observed_right_K"] is not None:
            print(
                f"    closest sampled K point right = "
                f"x={r['observed_right_K']}"
            )

        print()

    # ------------------------------------------------------------------
    # Detailed trajectory around a few examples.
    # ------------------------------------------------------------------

    section("TRAJECTORY SAMPLES")

    for r in results[:10]:

        print(
            f"n={r['n']:,} "
            f"mods=({r['r1']},{r['r2']}) "
            f"original K={r['K']}"
        )

        points = sorted(
            r["points"],
            key=lambda pt: abs(pt.x)
        )

        printed = 0

        for point in points:
            if point.x == 0:
                continue

            if abs(point.x) % TRAJECTORY_STEP != 0:
                continue

            state = (
                "SAME_K"
                if point.K == r["K"]
                else "K_CHANGED"
            )

            cell = (
                "SAME_CELL"
                if point.k == r["k"] and point.l == r["l"]
                else "CELL_CHANGED"
            )

            print(
                f"    x={point.x:+7d} "
                f"p={point.p:6d} "
                f"q={point.q:6d} "
                f"k={point.k:4d} "
                f"l={point.l:4d} "
                f"K={point.K:6d} "
                f"{state:10s} "
                f"{cell}"
            )

            printed += 1

            if printed >= 12:
                break

        print()

    # ------------------------------------------------------------------
    # Mathematical test.
    # ------------------------------------------------------------------

    section("MATHEMATICAL TEST")

    print(
        """
For the original factorization:

    p = a + k*r1
    q = b + l*r2

and therefore:

    k = floor(p/r1)
    l = floor(q/r2)
    K = k*l.

The quotient cell is:

    k*r1 <= p < (k+1)*r1
    l*r2 <= q < (l+1)*r2.

A nearby factorization:

    n+x = p_x*q_x

has:

    k_x = floor(p_x/r1)
    l_x = floor(q_x/r2).

The central hypothesis is:

    K_x = K

may remain stable over a substantial local x interval.

But there are two distinct mechanisms:

    1. FULL CELL STABILITY

           k_x = k
           l_x = l

       which automatically implies

           K_x = K.

    2. CROSS-CELL K COLLISIONS

           (k_x,l_x) != (k,l)
           but
           k_x*l_x = k*l.

       These are much more interesting.

If K remains constant only while the factors stay in the same
quotient cell, then the phenomenon is purely geometric.

If K frequently remains constant even after one coordinate crosses
a cell boundary, then K is carrying an additional arithmetic
invariant.

The experiment therefore explicitly distinguishes:

    SAME_CELL
    SAME_K / DIFFERENT_CELL.

The second category is the one to investigate further.

The geometric boundary distances are:

    d_p^- = p - k*r1
    d_p^+ = (k+1)*r1 - 1 - p

and similarly for q.

If the first change in K occurs close to the x-values implied by
these distances, the observed K-stability is explained by quotient
geometry.

If there is a substantial region beyond the first cell boundary
where K remains unchanged, then that is evidence for a different
mechanism.

The strongest possible observation is therefore:

    k_x,l_x change
        but
    k_x*l_x = k*l

over a nontrivial x interval.

That would justify a dedicated investigation of multiplicative
quotient trajectories.

The final factor pairs are exact because all nearby points are
constructed from primes.
"""
    )

    # ------------------------------------------------------------------
    # Final statistics for cross-cell invariance.
    # ------------------------------------------------------------------

    section("CROSS-CELL K INVARIANCE")

    cross_cell_same_K = 0
    cross_cell_total = 0
    strict_cross_cell_examples = []

    for r in results:

        k0 = r["k"]
        l0 = r["l"]
        K0 = r["K"]

        for pt in r["points"]:

            if pt.x == 0:
                continue

            different_cell = (
                pt.k != k0 or
                pt.l != l0
            )

            if not different_cell:
                continue

            cross_cell_total += 1

            if pt.K == K0:
                cross_cell_same_K += 1

                if len(strict_cross_cell_examples) < PRINT_EXAMPLES:
                    strict_cross_cell_examples.append(
                        (
                            r,
                            pt,
                        )
                    )

    print(
        f"cross-cell trajectory points       = "
        f"{cross_cell_total:,}"
    )

    print(
        f"cross-cell points with K_x=K       = "
        f"{cross_cell_same_K:,}"
    )

    if cross_cell_total:
        print(
            f"cross-cell K invariance fraction   = "
            f"{cross_cell_same_K / cross_cell_total:.6f}"
        )
    else:
        print(
            "cross-cell K invariance fraction   = 0.000000"
        )

    if strict_cross_cell_examples:

        print()
        print("CROSS-CELL SAME-K EXAMPLES")

        for r, pt in strict_cross_cell_examples:

            print(
                f"n={r['n']:,} "
                f"x={pt.x:+,} "
                f"mods=({r['r1']},{r['r2']})"
            )

            print(
                f"    original: "
                f"(k,l)=({r['k']},{r['l']}) "
                f"K={r['K']}"
            )

            print(
                f"    nearby:   "
                f"(k,l)=({pt.k},{pt.l}) "
                f"K={pt.K}"
            )

            print()

    # ------------------------------------------------------------------
    # Runtime.
    # ------------------------------------------------------------------

    total = time.perf_counter() - t0

    section("TIMING")

    print(f"total runtime               = {total:.3f} s")

    print()
    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    run()
