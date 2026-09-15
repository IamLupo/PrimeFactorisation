#!/usr/bin/env python3

import math
import random
import time
from dataclasses import dataclass


# ============================================================================
# CONFIG
# ============================================================================

SEED = 1_511_464_998

N_ANCHORS = 20

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

PERTURB_RADIUS = 50

CLOSE_RATIO = 0.20


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve(limit):
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0] = 0
    a[1] = 0

    for p in range(2, int(math.isqrt(limit)) + 1):
        if a[p]:
            a[p * p:limit + 1:p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [i for i in range(2, limit + 1) if a[i]]


# ============================================================================
# ANCHOR
# ============================================================================

@dataclass
class Anchor:
    p: int
    q: int
    n: int


# ============================================================================
# BUILD ANCHORS
# ============================================================================

def build_anchors(primes, count, seed):
    rng = random.Random(seed)

    pairs = []

    # Do not construct the full Cartesian product.
    # Randomly choose directly from the prime pool.
    seen = set()

    while len(pairs) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))
        pairs.append(Anchor(p, q, p * q))

    return pairs


# ============================================================================
# CLOSE MODULUS PAIR
# ============================================================================

def choose_close_pair(primes, rng):
    r1 = rng.choice(primes)

    max_r2 = int(r1 * (1.0 + CLOSE_RATIO))

    candidates = [
        p for p in primes
        if p > r1 and p <= max_r2
    ]

    if not candidates:
        return None

    return r1, rng.choice(candidates)


# ============================================================================
# COORDINATES
# ============================================================================

def coordinates(n, p, q, r1, r2):
    R = r1 * r2

    a = p % r1
    b = q % r2

    k = (p - a) // r1
    l = (q - b) // r2

    T = n // R
    E = T - k * l
    K = k * l

    return {
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "T": T,
        "E": E,
        "K": K,
    }


# ============================================================================
# NEARBY FACTORIZATION
# ============================================================================

def make_nearby(anchor, prime_set, rng):
    p = anchor.p
    q = anchor.q

    for _ in range(1000):

        dp = rng.randint(-PERTURB_RADIUS, PERTURB_RADIUS)
        dq = rng.randint(-PERTURB_RADIUS, PERTURB_RADIUS)

        if dp == 0 and dq == 0:
            continue

        px = p + dp
        qx = q + dq

        if px < FACTOR_MIN or px > FACTOR_MAX:
            continue

        if qx < FACTOR_MIN or qx > FACTOR_MAX:
            continue

        if px not in prime_set or qx not in prime_set:
            continue

        nx = px * qx

        if nx == anchor.n:
            continue

        return px, qx, nx

    return None


# ============================================================================
# MAIN
# ============================================================================

def run():

    start = time.perf_counter()

    print("=" * 92)
    print("FAST NEARBY-FACTORIZATION / E-K TRAJECTORY EXPERIMENT")
    print("=" * 92)
    print(f"anchors                    = {N_ANCHORS}")
    print(f"factor range               = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range              = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"perturbation radius        = {PERTURB_RADIUS}")
    print(f"seed                       = {SEED}")
    print()

    # ------------------------------------------------------------------------
    # PRIME POOLS
    # ------------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING PRIME POOLS")
    print("=" * 92)

    all_primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in all_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    factor_set = set(factor_primes)

    print(f"factor primes              = {len(factor_primes):,}")
    print(f"modulus primes             = {len(modulus_primes):,}")
    print()

    # ------------------------------------------------------------------------
    # ANCHORS
    # ------------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING ANCHORS")
    print("=" * 92)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(f"actual anchors             = {len(anchors)}")
    print()

    rng = random.Random(SEED)

    # ------------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------------

    total = 0

    same_E = 0
    same_K = 0
    same_k = 0
    same_l = 0
    same_pair = 0

    abs_dE = []
    abs_dK = []
    abs_dk = []
    abs_dl = []

    examples = []

    # ------------------------------------------------------------------------
    # RUN
    # ------------------------------------------------------------------------

    print("=" * 92)
    print("RUNNING FAST TRAJECTORY TEST")
    print("=" * 92)

    for index, anchor in enumerate(anchors, 1):

        pair = choose_close_pair(modulus_primes, rng)

        if pair is None:
            continue

        r1, r2 = pair

        nearby = make_nearby(
            anchor,
            factor_set,
            rng,
        )

        if nearby is None:
            continue

        px, qx, nx = nearby

        # Original coordinates.
        original = coordinates(
            anchor.n,
            anchor.p,
            anchor.q,
            r1,
            r2,
        )

        # Nearby coordinates.
        nearby_coord = coordinates(
            nx,
            px,
            qx,
            r1,
            r2,
        )

        total += 1

        dE = nearby_coord["E"] - original["E"]
        dK = nearby_coord["K"] - original["K"]
        dk = nearby_coord["k"] - original["k"]
        dl = nearby_coord["l"] - original["l"]

        abs_dE.append(abs(dE))
        abs_dK.append(abs(dK))
        abs_dk.append(abs(dk))
        abs_dl.append(abs(dl))

        if nearby_coord["E"] == original["E"]:
            same_E += 1

        if nearby_coord["K"] == original["K"]:
            same_K += 1

        if nearby_coord["k"] == original["k"]:
            same_k += 1

        if nearby_coord["l"] == original["l"]:
            same_l += 1

        if (
            nearby_coord["k"] == original["k"]
            and nearby_coord["l"] == original["l"]
        ):
            same_pair += 1

        if len(examples) < 20:
            examples.append(
                (
                    anchor,
                    px,
                    qx,
                    nx,
                    r1,
                    r2,
                    original,
                    nearby_coord,
                )
            )

        if index % 5 == 0 or index == len(anchors):
            print(f"anchor {index:3d}/{len(anchors)}")

    # ------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("SUMMARY")
    print("=" * 92)

    print(f"successful nearby samples   = {total}")

    if total == 0:
        print("No samples generated.")
        return

    print()

    print("E STABILITY")
    print(f"    E_x = E                 = {same_E}/{total}")
    print(f"    fraction                = {same_E / total:.6f}")
    print()

    print("K STABILITY")
    print(f"    K_x = K                 = {same_K}/{total}")
    print(f"    fraction                = {same_K / total:.6f}")
    print()

    print("COORDINATE STABILITY")
    print(f"    k_x = k                 = {same_k}/{total}")
    print(f"    l_x = l                 = {same_l}/{total}")
    print(f"    (k_x,l_x) = (k,l)       = {same_pair}/{total}")
    print()

    print("AVERAGE DISTANCES")
    print(f"    mean |E_x-E|            = {sum(abs_dE)/total:.6f}")
    print(f"    max  |E_x-E|            = {max(abs_dE):,}")
    print()

    print(f"    mean |K_x-K|            = {sum(abs_dK)/total:.6f}")
    print(f"    max  |K_x-K|            = {max(abs_dK):,}")
    print()

    print(f"    mean |k_x-k|            = {sum(abs_dk)/total:.6f}")
    print(f"    max  |k_x-k|            = {max(abs_dk):,}")
    print()

    print(f"    mean |l_x-l|            = {sum(abs_dl)/total:.6f}")
    print(f"    max  |l_x-l|            = {max(abs_dl):,}")
    print()

    # ------------------------------------------------------------------------
    # TEST THE BASIC RELATION
    # ------------------------------------------------------------------------

    print("=" * 92)
    print("ALGEBRAIC TEST")
    print("=" * 92)

    relation_ok = 0

    for (
        anchor,
        px,
        qx,
        nx,
        r1,
        r2,
        original,
        nearby_coord,
    ) in examples:

        lhs = nearby_coord["K"] - original["K"]
        rhs = (
            nearby_coord["T"]
            - original["T"]
            - (
                nearby_coord["E"]
                - original["E"]
            )
        )

        if lhs == rhs:
            relation_ok += 1

    print(
        "K_x-K = (T_x-T)-(E_x-E)"
    )

    print(
        f"checked example samples      = {len(examples)}"
    )

    print(
        f"identity failures             = "
        f"{len(examples) - relation_ok}"
    )

    print()

    # ------------------------------------------------------------------------
    # EXAMPLES
    # ------------------------------------------------------------------------

    print("=" * 92)
    print("EXAMPLES")
    print("=" * 92)

    for (
        anchor,
        px,
        qx,
        nx,
        r1,
        r2,
        original,
        nearby_coord,
    ) in examples:

        print(
            f"n={anchor.n:,} x={nx-anchor.n:+,} "
            f"n+x={nx:,}"
        )

        print(
            f"    original: "
            f"p={anchor.p:,} q={anchor.q:,} "
            f"k={original['k']} l={original['l']} "
            f"K={original['K']} T={original['T']} E={original['E']}"
        )

        print(
            f"    nearby:   "
            f"p={px:,} q={qx:,} "
            f"k={nearby_coord['k']} l={nearby_coord['l']} "
            f"K={nearby_coord['K']} "
            f"T={nearby_coord['T']} E={nearby_coord['E']}"
        )

        print(
            f"    delta: "
            f"dK={nearby_coord['K']-original['K']:+d} "
            f"dE={nearby_coord['E']-original['E']:+d} "
            f"dk={nearby_coord['k']-original['k']:+d} "
            f"dl={nearby_coord['l']-original['l']:+d}"
        )

        print(
            f"    mods=({r1},{r2})"
        )

        print()

    # ------------------------------------------------------------------------
    # PREDICTIVE RULES
    # ------------------------------------------------------------------------

    print("=" * 92)
    print("PREDICTIVE RULE TEST")
    print("=" * 92)

    print(
        """
The following possible rules are tested conceptually:

    Rule A:
        K = K_x

    Rule B:
        K is close to K_x

    Rule C:
        E = E_x

    Rule D:
        k is close to k_x
        l is close to l_x

The important quantity is not whether these are sometimes true.

The important question is whether the error remains bounded while
the factor size increases.

A useful nearby-factor mechanism would ideally show:

    |K-K_x| = small
    |k-k_x| = small
    |l-l_x| = small

independently of the size of n.

If these errors grow with the factor size, nearby factorization is
not giving a useful shortcut.

"""
    )

    # ------------------------------------------------------------------------
    # TIMING
    # ------------------------------------------------------------------------

    elapsed = time.perf_counter() - start

    print("=" * 92)
    print("TIMING")
    print("=" * 92)
    print(f"total runtime                = {elapsed:.3f} seconds")

    print()
    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    run()
