#!/usr/bin/env python3

import math
import random
import statistics
import time
from collections import Counter


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ANCHORS = 300
MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20
SEED = 1_511_464_998

PROGRESS_EVERY = 25


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    sieve_ = bytearray(b"\x01") * (limit + 1)
    sieve_[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve_[p]:
            start = p * p
            sieve_[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, v in enumerate(sieve_) if v]


# =============================================================================
# PRIME POOLS
# =============================================================================

def build_prime_pools():
    primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in primes
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    return factor_primes, modulus_primes


# =============================================================================
# ANCHORS
# =============================================================================

def build_anchors(factor_primes: list[int]) -> list[tuple[int, int, int]]:
    """
    Deterministically generate distinct prime factor pairs.

    The experiment is interested in structural behavior rather than
    reproducing a particular previous random list byte-for-byte.
    """

    rng = random.Random(SEED)

    pairs = set()

    while len(pairs) < ANCHORS:
        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        if p >= q:
            continue

        n = p * q

        if not (FACTOR_MIN * FACTOR_MIN <= n <= FACTOR_MAX * FACTOR_MAX):
            continue

        pairs.add((p, q, n))

    return sorted(pairs, key=lambda x: x[2])


# =============================================================================
# CLOSE MODULUS TRIPLE
# =============================================================================

def select_close_triple(
    n: int,
    modulus_primes: list[int],
) -> tuple[int, int, int]:

    target = n ** (1.0 / 3.0)

    lo = target * (1.0 - CLOSE_RATIO)
    hi = target * (1.0 + CLOSE_RATIO)

    candidates = [
        r for r in modulus_primes
        if lo <= r <= hi
    ]

    best = None
    best_product = -1

    # Search all triples in the relatively small local neighborhood.
    for i in range(len(candidates)):
        r1 = candidates[i]

        for j in range(i + 1, len(candidates)):
            r2 = candidates[j]

            product12 = r1 * r2

            if product12 >= n:
                break

            for k in range(j + 1, len(candidates)):
                r3 = candidates[k]

                R = product12 * r3

                if R > n:
                    break

                if R > best_product:
                    best_product = R
                    best = (r1, r2, r3)

    if best is None:
        raise RuntimeError(f"No valid modulus triple for n={n}")

    return best


# =============================================================================
# MODULAR INVERSE TABLE
# =============================================================================

def inverse_table(prime: int) -> list[int]:
    """
    inv[x] = x^-1 mod prime for x != 0.

    Since the prime is at most 3000 this is tiny.
    """

    inv = [0] * prime
    inv[1] = 1

    for x in range(2, prime):
        inv[x] = prime - (prime // x) * inv[prime % x] % prime

    return inv


# =============================================================================
# CRT HELPERS
# =============================================================================

def crt_three(
    a1: int,
    a2: int,
    a3: int,
    r1: int,
    r2: int,
    r3: int,
    inv_r1_mod_r2: int,
    inv_r1r2_mod_r3: int,
) -> tuple[int, int]:
    """
    Combine:

        x = a1 (mod r1)
        x = a2 (mod r2)
        x = a3 (mod r3)

    Returns:

        x0, modulus

    with:

        0 <= x0 < modulus
    """

    # First combine r1 and r2.
    t = ((a2 - a1) * inv_r1_mod_r2) % r2
    x12 = a1 + r1 * t
    R12 = r1 * r2

    # Then combine with r3.
    u = ((a3 - x12) * inv_r1r2_mod_r3) % r3
    x = x12 + R12 * u

    return x, R12 * r3


# =============================================================================
# HYPERBOLA INTERVAL
# =============================================================================

def q_interval(
    n: int,
    p: int,
    factor_min: int,
    factor_max: int,
) -> tuple[int, int] | None:

    # Integer q satisfying:

    #     p*q <= n < p*(q+1)

    # equivalently:

    #     floor(n/(p+1)) + 1 <= q <= floor(n/p)

    q_lo = n // (p + 1) + 1
    q_hi = n // p

    q_lo = max(q_lo, factor_min)
    q_hi = min(q_hi, factor_max)

    if q_lo > q_hi:
        return None

    return q_lo, q_hi


# =============================================================================
# FIRST-MODULUS-ONLY TEST
# =============================================================================

def residue_intersects_interval(
    q_lo: int,
    q_hi: int,
    residue: int,
    modulus: int,
) -> bool:

    if residue == 0:
        first = modulus
    else:
        first = residue

    if first < q_lo:
        first += ((q_lo - first + modulus - 1) // modulus) * modulus

    return first <= q_hi


# =============================================================================
# MAIN PER-ANCHOR SEARCH
# =============================================================================

def run_anchor(
    p_true: int,
    q_true: int,
    n: int,
    mods: tuple[int, int, int],
):
    r1, r2, r3 = mods

    R = r1 * r2 * r3

    # CRT precomputation.
    inv_r1_mod_r2 = pow(r1, -1, r2)
    inv_r1r2_mod_r3 = pow((r1 * r2) % r3, -1, r3)

    # Inverse tables.
    inv1 = inverse_table(r1)
    inv2 = inverse_table(r2)
    inv3 = inverse_table(r3)

    p_limit = min(
        FACTOR_MAX,
        int(math.isqrt(n)),
    )

    if p_limit < FACTOR_MIN:
        return {
            "p_states": 0,
            "r1_survivors": 0,
            "r12_survivors": 0,
            "crt_survivors": 0,
            "exact_solutions": 0,
            "true_seen": False,
        }

    p_states = 0
    r1_survivors = 0
    r12_survivors = 0
    crt_survivors = 0
    exact_solutions = 0

    true_seen = False

    survivors = []

    for p in range(FACTOR_MIN, p_limit + 1):

        p_states += 1

        # p must be invertible modulo all three primes.
        a1 = p % r1
        a2 = p % r2
        a3 = p % r3

        if a1 == 0 or a2 == 0 or a3 == 0:
            continue

        q_range = q_interval(
            n,
            p,
            FACTOR_MIN,
            FACTOR_MAX,
        )

        if q_range is None:
            continue

        q_lo, q_hi = q_range

        # q = n / p (mod ri)
        b1 = (n % r1) * inv1[a1] % r1

        if not residue_intersects_interval(
            q_lo,
            q_hi,
            b1,
            r1,
        ):
            continue

        r1_survivors += 1

        b2 = (n % r2) * inv2[a2] % r2

        if not residue_intersects_interval(
            q_lo,
            q_hi,
            b2,
            r2,
        ):
            continue

        r12_survivors += 1

        b3 = (n % r3) * inv3[a3] % r3

        if not residue_intersects_interval(
            q_lo,
            q_hi,
            b3,
            r3,
        ):
            continue

        # Combine all three q congruences.
        q0, modulus = crt_three(
            b1,
            b2,
            b3,
            r1,
            r2,
            r3,
            inv_r1_mod_r2,
            inv_r1r2_mod_r3,
        )

        # Find the first q in [q_lo, q_hi] in the CRT class.
        if q0 == 0:
            first_q = modulus
        else:
            first_q = q0

        if first_q < q_lo:
            first_q += (
                (q_lo - first_q + modulus - 1)
                // modulus
            ) * modulus

        if first_q > q_hi:
            continue

        crt_survivors += 1

        # There can be at most one such q because modulus is normally
        # larger than the factor interval. Still handle generality.
        q = first_q

        while q <= q_hi:

            survivors.append((p, q))

            if p * q == n:
                exact_solutions += 1

                if (
                    (p == p_true and q == q_true)
                    or
                    (p == q_true and q == p_true)
                ):
                    true_seen = True

            q += modulus

    return {
        "p_states": p_states,
        "r1_survivors": r1_survivors,
        "r12_survivors": r12_survivors,
        "crt_survivors": crt_survivors,
        "exact_solutions": exact_solutions,
        "true_seen": true_seen,
        "survivors": survivors,
    }


# =============================================================================
# REPORTING
# =============================================================================

def pct(a: float, b: float) -> float:
    if b == 0:
        return 0.0
    return 100.0 * a / b


def run():

    total_start = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME DIRECT P-RESIDUE / CRT-HYPERBOLA FEASIBILITY EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes, modulus_primes = build_prime_pools()

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(factor_primes)

    print(f"actual anchors             = {len(anchors)}")
    print()

    print("=" * 100)
    print("SELECTING ANCHOR-SPECIFIC CLOSE TRIPLES")
    print("=" * 100)

    selected = []

    for idx, (p, q, n) in enumerate(anchors, 1):

        mods = select_close_triple(
            n,
            modulus_primes,
        )

        selected.append(
            (p, q, n, mods)
        )

        if idx % PROGRESS_EVERY == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    print()

    print("=" * 100)
    print("VALIDATING MODULUS TRIPLES")
    print("=" * 100)

    validation_failures = 0
    unique_triples = set()

    for p, q, n, mods in selected:

        r1, r2, r3 = mods

        if not (
            MODULUS_MIN <= r1 <= MODULUS_MAX
            and MODULUS_MIN <= r2 <= MODULUS_MAX
            and MODULUS_MIN <= r3 <= MODULUS_MAX
            and r1 < r2 < r3
            and all(
                math.gcd(a, b) == 1
                for a, b in [
                    (r1, r2),
                    (r1, r3),
                    (r2, r3),
                ]
            )
        ):
            validation_failures += 1

        unique_triples.add(mods)

    print(f"triple validation failures = {validation_failures}")
    print(f"unique modulus triples      = {len(unique_triples)}")
    print()

    if validation_failures:
        raise RuntimeError("Invalid modulus triple selection.")

    print("=" * 100)
    print("RUNNING DIRECT CRT-HYPERBOLA SEARCH")
    print("=" * 100)

    search_start = time.perf_counter()

    totals = Counter()

    exact_recovery = 0

    r1_list = []
    r12_list = []
    crt_list = []
    pstate_list = []

    strongest = []
    weakest = []

    for idx, (p_true, q_true, n, mods) in enumerate(selected, 1):

        result = run_anchor(
            p_true,
            q_true,
            n,
            mods,
        )

        p_states = result["p_states"]
        r1_survivors = result["r1_survivors"]
        r12_survivors = result["r12_survivors"]
        crt_survivors = result["crt_survivors"]
        exact_solutions = result["exact_solutions"]

        totals["p_states"] += p_states
        totals["r1_survivors"] += r1_survivors
        totals["r12_survivors"] += r12_survivors
        totals["crt_survivors"] += crt_survivors
        totals["exact_solutions"] += exact_solutions

        pstate_list.append(p_states)
        r1_list.append(r1_survivors)
        r12_list.append(r12_survivors)
        crt_list.append(crt_survivors)

        if result["true_seen"]:
            exact_recovery += 1

        record = (
            crt_survivors,
            p_true,
            q_true,
            n,
            mods,
            p_states,
            r1_survivors,
            r12_survivors,
            exact_solutions,
        )

        strongest.append(record)
        strongest.sort(key=lambda x: x[0])

        if len(strongest) > 20:
            strongest.pop()

        weakest.append(record)
        weakest.sort(key=lambda x: x[0], reverse=True)

        if len(weakest) > 20:
            weakest.pop()

        if idx % PROGRESS_EVERY == 0:
            print(f"anchor {idx:3d}/{len(selected)}")

    search_time = time.perf_counter() - search_start
    total_time = time.perf_counter() - total_start

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------

    avg_p_states = statistics.mean(pstate_list)
    avg_r1 = statistics.mean(r1_list)
    avg_r12 = statistics.mean(r12_list)
    avg_crt = statistics.mean(crt_list)
    avg_exact = totals["exact_solutions"] / len(selected)

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(f"anchors analyzed            = {len(selected)}")
    print(f"average integer p states    = {avg_p_states:,.3f}")
    print(f"average r1 survivors        = {avg_r1:,.3f}")
    print(f"average r1+r2 survivors     = {avg_r12:,.3f}")
    print(f"average three-modulus CRT    = {avg_crt:,.3f}")
    print(f"average exact solutions      = {avg_exact:,.3f}")
    print()

    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    print(
        f"r1 / p states               = "
        f"{avg_r1 / avg_p_states:.12f}"
    )

    print(
        f"r1+r2 / p states            = "
        f"{avg_r12 / avg_p_states:.12f}"
    )

    print(
        f"CRT / p states              = "
        f"{avg_crt / avg_p_states:.12f}"
    )

    print(
        f"CRT reduction               = "
        f"{100.0 * (1.0 - avg_crt / avg_p_states):.6f}%"
    )

    print()

    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(f"correctly recovered         = {exact_recovery}/{len(selected)}")
    print(
        f"recovery rate               = "
        f"{pct(exact_recovery, len(selected)):.4f}%"
    )

    print()

    print("=" * 100)
    print("STRONGEST CRT COLLAPSES")
    print("=" * 100)

    for (
        crt_survivors,
        p_true,
        q_true,
        n,
        mods,
        p_states,
        r1_survivors,
        r12_survivors,
        exact_solutions,
    ) in strongest:

        print(
            f"n={n:,} "
            f"p={p_true:,} "
            f"q={q_true:,} "
            f"mods={mods} "
            f"p={p_states:,} "
            f"r1={r1_survivors:,} "
            f"r12={r12_survivors:,} "
            f"CRT={crt_survivors:,} "
            f"exact={exact_solutions}"
        )

    print()

    print("=" * 100)
    print("WEAKEST CRT COLLAPSES")
    print("=" * 100)

    for (
        crt_survivors,
        p_true,
        q_true,
        n,
        mods,
        p_states,
        r1_survivors,
        r12_survivors,
        exact_solutions,
    ) in weakest:

        print(
            f"n={n:,} "
            f"p={p_true:,} "
            f"q={q_true:,} "
            f"mods={mods} "
            f"p={p_states:,} "
            f"r1={r1_survivors:,} "
            f"r12={r12_survivors:,} "
            f"CRT={crt_survivors:,} "
            f"exact={exact_solutions}"
        )

    print()

    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(f"search runtime               = {search_time:.3f} seconds")
    print(f"total runtime                = {total_time:.3f} seconds")

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print()
    print(
        """
For every possible integer p in the factor interval:

    q_lo <= q <= q_hi

is obtained directly from the hyperbola:

    q_lo = floor(n/(p+1)) + 1
    q_hi = floor(n/p).

For a genuine factor pair:

    p*q = n

and therefore for every modulus r_i:

    q = n * p^(-1) (mod r_i).

The three conditions are combined into one CRT class:

    q = Q_p (mod r1*r2*r3).

The experiment then asks:

    Does this CRT class intersect the hyperbola interval?

No enumeration of:

    a
    a3
    b
    b3
    k
    l
    q

is performed.

Therefore this is a clean test of whether the three close moduli
can prune the p-domain before constructing candidate q values.

There are two fundamentally different outcomes.

CASE A:

    CRT survivors << integer p states

This proves that the modular/hyperbolic combination is an extremely
strong necessary-condition filter.

But it does NOT yet give a fast factorization algorithm, because the
experiment still scanned every p.

CASE B:

    CRT survivors is tiny, but obtaining them requires scanning
    essentially all p.

Then the modular information is useful for candidate verification,
but it has not removed the fundamental p-search.

The important next target is therefore NOT merely:

    CRT survivors ~= 1.

The real target is:

    Can the surviving p values themselves be generated directly,
    without scanning the complete p interval?

That is the point at which a genuinely new factor-search mechanism
would begin to emerge.

Every survivor is finally tested with:

    p*q == n

so exact recovery is separate from modular survival.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
