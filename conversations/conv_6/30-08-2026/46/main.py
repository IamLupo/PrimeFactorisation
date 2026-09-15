#!/usr/bin/env python3

"""
================================================================================
RECURSIVE DEFECT-CONTRACTION / OBSERVABLE E-SET EXPERIMENT
================================================================================

Question:

    Does the E-structure recurse in a useful way?

Level 1:

    T1 = floor(n / (r1*r2))

    T1 = k*l + E1

Level 2:

    k = u*s1 + alpha
    l = v*s2 + beta

    T2 = floor(T1 / (s1*s2))

    T2 = u*v + E2

where

    E2 =
        floor(
            u*s1*beta
            + v*s2*alpha
            + alpha*beta
            + E1
        / (s1*s2)
        )

Level 3 repeats the same construction.

IMPORTANT:

This experiment does NOT follow only the true factor branch.

It constructs the OBSERVABLE candidate set:

    E1 candidates
        ->
    candidate (k,l) pairs
        ->
    candidate E2 values
        ->
    candidate (u,v) pairs
        ->
    candidate E3 values
        -> ...

The question is whether the uncertainty itself contracts.

================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict


# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

# Number of recursive levels after the initial E1 level.
DEPTH = 4

# Small radix pool.
RADIX_MIN = 3
RADIX_MAX = 31

# Number of radix pairs to test per anchor/candidate.
RADIX_PAIR_COUNT = 3


# ==============================================================================
# PRIME GENERATION
# ==============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


# ==============================================================================
# FACTORIZATION
# ==============================================================================

def factorize(n: int) -> list[tuple[int, int]]:
    """
    Return prime factors as (prime, exponent).
    """
    out: list[tuple[int, int]] = []

    d = 2
    x = n

    while d * d <= x:
        if x % d == 0:
            e = 0
            while x % d == 0:
                x //= d
                e += 1
            out.append((d, e))
        d = 3 if d == 2 else d + 2

    if x > 1:
        out.append((x, 1))

    return out


def divisors(n: int) -> list[int]:
    fs = factorize(n)

    ds = [1]

    for p, e in fs:
        old = ds[:]
        mult = 1

        expanded = []
        for power in range(e + 1):
            expanded.extend(d * mult for d in old)
            mult *= p

        ds = expanded

    return sorted(set(ds))


# ==============================================================================
# CLOSE MODULUS PAIRS
# ==============================================================================

def choose_close_pair(
    modulus_primes: list[int],
    rng: random.Random,
) -> tuple[int, int]:

    pairs = []

    for i, r1 in enumerate(modulus_primes):
        for r2 in modulus_primes[i + 1:]:
            if r2 / r1 <= 1.0 + CLOSE_RATIO:
                pairs.append((r1, r2))

    if not pairs:
        raise RuntimeError("No close modulus pairs found.")

    return rng.choice(pairs)


# ==============================================================================
# BUILD ANCHORS
# ==============================================================================

def build_anchors(
    factor_primes: list[int],
    modulus_primes: list[int],
    count: int,
    rng: random.Random,
):
    anchors = []

    seen = set()

    while len(anchors) < count:

        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        if p == q:
            continue

        n = p * q

        if n in seen:
            continue

        if not (FACTOR_MIN * FACTOR_MIN <= n <= FACTOR_MAX * FACTOR_MAX):
            continue

        r1, r2 = choose_close_pair(modulus_primes, rng)

        # We want r1*r2 to be below n so T1 is meaningful.
        if r1 * r2 >= n:
            continue

        anchors.append({
            "n": n,
            "p": min(p, q),
            "q": max(p, q),
            "r1": r1,
            "r2": r2,
        })

        seen.add(n)

    return anchors


# ==============================================================================
# LEVEL-1 CANDIDATE GENERATION
# ==============================================================================

def level1_candidates(
    n: int,
    r1: int,
    r2: int,
    p_min: int,
    p_max: int,
    q_min: int,
    q_max: int,
):
    """
    Observable candidates.

    T1 = floor(n/(r1*r2))

    E1 = T1 - k*l

    Use a rigorous upper bound:

        E1 <= floor(p_max/r1) + floor(q_max/r2)

    Then factor every possible K = T1-E1 and retain divisor pairs
    that fit the quotient-coordinate bounds.
    """

    T1 = n // (r1 * r2)

    k_min = math.ceil(p_min / r1)
    k_max = p_max // r1

    l_min = math.ceil(q_min / r2)
    l_max = q_max // r2

    e_max = k_max + l_max

    candidates = []

    for e1 in range(e_max + 1):

        K = T1 - e1

        if K <= 0:
            continue

        for k in divisors(K):

            if k < k_min or k > k_max:
                continue

            if K % k:
                continue

            l = K // k

            if l < l_min or l > l_max:
                continue

            candidates.append((k, l, e1))

    return T1, e_max, candidates


# ==============================================================================
# RECURSIVE DEFECT TRANSITION
# ==============================================================================

def descend_pair(
    x: int,
    y: int,
    E: int,
    s1: int,
    s2: int,
):
    """
    Given:

        x = u*s1 + alpha
        y = v*s2 + beta

    and:

        current quantity = x*y + E

    compute:

        T_next = floor((x*y + E)/(s1*s2))
        T_next = u*v + E_next

    using the exact recursive identity.
    """

    u, alpha = divmod(x, s1)
    v, beta = divmod(y, s2)

    current = x * y + E

    T_next = current // (s1 * s2)

    E_next = (
        u * s1 * beta
        + v * s2 * alpha
        + alpha * beta
        + E
    ) // (s1 * s2)

    assert T_next == u * v + E_next

    return {
        "x": x,
        "y": y,
        "u": u,
        "v": v,
        "alpha": alpha,
        "beta": beta,
        "E_in": E,
        "T_next": T_next,
        "E_next": E_next,
    }


# ==============================================================================
# RADIX SELECTION
# ==============================================================================

def choose_radix_pairs(
    radix_primes: list[int],
    count: int,
    rng: random.Random,
):
    all_pairs = []

    for s1 in radix_primes:
        for s2 in radix_primes:
            if s1 * s2 <= 0:
                continue
            all_pairs.append((s1, s2))

    if len(all_pairs) <= count:
        return all_pairs

    return rng.sample(all_pairs, count)


# ==============================================================================
# VALIDATE TRUE BRANCH
# ==============================================================================

def validate_true_branch(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    radix_pairs,
):
    T1 = n // (r1 * r2)

    k, a = divmod(p, r1)
    l, b = divmod(q, r2)

    E1 = T1 - k * l

    assert n == (k * r1 + a) * (l * r2 + b)
    assert T1 == k * l + E1

    out = []

    x = k
    y = l
    E = E1

    for depth, (s1, s2) in enumerate(radix_pairs, start=2):

        state = descend_pair(
            x,
            y,
            E,
            s1,
            s2,
        )

        out.append(state)

        x = state["u"]
        y = state["v"]
        E = state["E_next"]

        if x == 0 or y == 0:
            break

    return E1, out


# ==============================================================================
# OBSERVABLE RECURSIVE ANALYSIS
# ==============================================================================

def analyze_anchor(
    anchor,
    radix_pairs,
):
    n = anchor["n"]
    r1 = anchor["r1"]
    r2 = anchor["r2"]

    p = anchor["p"]
    q = anchor["q"]

    T1, emax, level1 = level1_candidates(
        n=n,
        r1=r1,
        r2=r2,
        p_min=FACTOR_MIN,
        p_max=FACTOR_MAX,
        q_min=FACTOR_MIN,
        q_max=FACTOR_MAX,
    )

    # --------------------------------------------------------------------------
    # Level-1 defect set
    # --------------------------------------------------------------------------

    E1_set = {e for _, _, e in level1}

    # --------------------------------------------------------------------------
    # Recursive observable states
    # --------------------------------------------------------------------------

    states = [
        {
            "x": k,
            "y": l,
            "E": e1,
        }
        for k, l, e1 in level1
    ]

    level_stats = []

    # Level 1 statistics.
    level_stats.append({
        "level": 1,
        "states": len(states),
        "unique_E": len(E1_set),
        "E_min": min(E1_set) if E1_set else None,
        "E_max": max(E1_set) if E1_set else None,
        "products": len({k * l for k, l, _ in level1}),
    })

    # --------------------------------------------------------------------------
    # Recursive descent
    # --------------------------------------------------------------------------

    for depth, radix_pair in enumerate(radix_pairs, start=2):

        s1, s2 = radix_pair

        next_states = []

        for state in states:

            if state["x"] <= 0 or state["y"] <= 0:
                continue

            descended = descend_pair(
                state["x"],
                state["y"],
                state["E"],
                s1,
                s2,
            )

            u = descended["u"]
            v = descended["v"]
            e_next = descended["E_next"]

            if u < 0 or v < 0:
                continue

            next_states.append({
                "x": u,
                "y": v,
                "E": e_next,
            })

        states = next_states

        Es = {s["E"] for s in states}
        products = {s["x"] * s["y"] for s in states}

        level_stats.append({
            "level": depth,
            "states": len(states),
            "unique_E": len(Es),
            "E_min": min(Es) if Es else None,
            "E_max": max(Es) if Es else None,
            "products": len(products),
        })

        if not states:
            break

    # --------------------------------------------------------------------------
    # True branch validation
    # --------------------------------------------------------------------------

    k_true, a_true = divmod(p, r1)
    l_true, b_true = divmod(q, r2)

    E1_true = T1 - k_true * l_true

    true_branch_present = any(
        s["x"] == k_true and
        s["y"] == l_true and
        s["E"] == E1_true
        for s in [
            {
                "x": k,
                "y": l,
                "E": e,
            }
            for k, l, e in level1
        ]
    )

    true_chain = []

    x = k_true
    y = l_true
    E = E1_true

    for s1, s2 in radix_pairs:

        if x <= 0 or y <= 0:
            break

        st = descend_pair(
            x,
            y,
            E,
            s1,
            s2,
        )

        true_chain.append(st)

        x = st["u"]
        y = st["v"]
        E = st["E_next"]

    return {
        "T1": T1,
        "Emax": emax,
        "level1_count": len(level1),
        "E1_set": E1_set,
        "level_stats": level_stats,
        "true_E1": E1_true,
        "true_k": k_true,
        "true_l": l_true,
        "true_branch_present": true_branch_present,
        "true_chain": true_chain,
        "radix_pairs": radix_pairs,
    }


# ==============================================================================
# MAIN
# ==============================================================================

def run():
    start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("RECURSIVE DEFECT-CONTRACTION / OBSERVABLE E-SET EXPERIMENT")
    print("=" * 100)
    print(f"N anchors                    = {N_ANCHORS:,}")
    print(f"factor range                 = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range                = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio                  = {CLOSE_RATIO:.0%}")
    print(f"recursive depth              = {DEPTH}")
    print(f"radix range                  = {RADIX_MIN} - {RADIX_MAX}")
    print(f"radix pairs / anchor         = {RADIX_PAIR_COUNT}")
    print(f"seed                         = {SEED:,}")
    print()

    # --------------------------------------------------------------------------
    # Prime pools
    # --------------------------------------------------------------------------

    factor_primes = [
        p for p in sieve(FACTOR_MAX)
        if p >= FACTOR_MIN
    ]

    modulus_primes = [
        p for p in sieve(MOD_MAX)
        if p >= MOD_MIN
    ]

    radix_primes = [
        p for p in sieve(RADIX_MAX)
        if p >= RADIX_MIN
    ]

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)
    print(f"factor primes                 = {len(factor_primes):,}")
    print(f"modulus primes                = {len(modulus_primes):,}")
    print(f"radix primes                  = {len(radix_primes):,}")
    print()

    # --------------------------------------------------------------------------
    # Anchors
    # --------------------------------------------------------------------------

    anchors = build_anchors(
        factor_primes,
        modulus_primes,
        N_ANCHORS,
        rng,
    )

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)
    print(f"actual anchors                = {len(anchors)}")
    print()

    # --------------------------------------------------------------------------
    # Analyze
    # --------------------------------------------------------------------------

    results = []

    for idx, anchor in enumerate(anchors, start=1):

        radix_pairs = choose_radix_pairs(
            radix_primes,
            RADIX_PAIR_COUNT,
            rng,
        )

        result = analyze_anchor(
            anchor,
            radix_pairs,
        )

        results.append((anchor, result))

        if idx % 10 == 0 or idx == len(anchors):
            print(f"anchor {idx:3d}/{len(anchors)}")

    # ==========================================================================
    # Aggregate
    # ==========================================================================

    def avg(values):
        values = [v for v in values if v is not None]
        return statistics.mean(values) if values else 0.0

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(f"anchors analyzed              = {len(results)}")
    print()

    print("LEVEL-1 OBSERVABLE SET")
    print(
        f"average E1 candidates         = "
        f"{avg([r['level1_count'] for _, r in results]):.3f}"
    )
    print(
        f"average distinct E1 values    = "
        f"{avg([len(r['E1_set']) for _, r in results]):.3f}"
    )
    print(
        f"average Emax                  = "
        f"{avg([r['Emax'] for _, r in results]):.3f}"
    )

    print()

    # --------------------------------------------------------------------------
    # Per-depth aggregate statistics
    # --------------------------------------------------------------------------

    for depth in range(1, DEPTH + 1):

        rows = []

        for _, r in results:
            for st in r["level_stats"]:
                if st["level"] == depth:
                    rows.append(st)

        if not rows:
            continue

        print(f"LEVEL {depth}")
        print(
            f"  average states              = "
            f"{avg([x['states'] for x in rows]):.3f}"
        )
        print(
            f"  average distinct E          = "
            f"{avg([x['unique_E'] for x in rows]):.3f}"
        )
        print(
            f"  average distinct products   = "
            f"{avg([x['products'] for x in rows]):.3f}"
        )
        print(
            f"  average E width             = "
            f"{avg([x['E_max'] - x['E_min'] if x['E_min'] is not None else 0 for x in rows]):.3f}"
        )

    # --------------------------------------------------------------------------
    # Contraction ratios
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("DEFECT-SET CONTRACTION")
    print("=" * 100)

    for depth in range(2, DEPTH + 1):

        ratios = []
        product_ratios = []
        state_ratios = []

        for _, r in results:

            stats = {
                st["level"]: st
                for st in r["level_stats"]
            }

            if depth not in stats or (depth - 1) not in stats:
                continue

            prev = stats[depth - 1]
            cur = stats[depth]

            if prev["unique_E"] > 0:
                ratios.append(
                    cur["unique_E"] / prev["unique_E"]
                )

            if prev["products"] > 0:
                product_ratios.append(
                    cur["products"] / prev["products"]
                )

            if prev["states"] > 0:
                state_ratios.append(
                    cur["states"] / prev["states"]
                )

        print(
            f"E-set level {depth-1} -> {depth}   = "
            f"{avg(ratios):.6f}"
        )
        print(
            f"product set {depth-1} -> {depth} = "
            f"{avg(product_ratios):.6f}"
        )
        print(
            f"state set {depth-1} -> {depth}   = "
            f"{avg(state_ratios):.6f}"
        )

    # --------------------------------------------------------------------------
    # True branch retention
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("TRUE BRANCH RETENTION")
    print("=" * 100)

    retained = sum(
        1
        for _, r in results
        if r["true_branch_present"]
    )

    print(
        f"true (k,l) present at level 1 = "
        f"{retained}/{len(results)}"
    )

    # --------------------------------------------------------------------------
    # E1 -> E2 examples
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("DEFECT CHAIN EXAMPLES")
    print("=" * 100)

    for anchor, result in results[:20]:

        print(
            f"n={anchor['n']:,} "
            f"p={anchor['p']:,} "
            f"q={anchor['q']:,} "
            f"mods=({anchor['r1']},{anchor['r2']})"
        )

        print(
            f"    T1={result['T1']} "
            f"E1_true={result['true_E1']} "
            f"E1_set={len(result['E1_set'])}"
        )

        for i, st in enumerate(result["true_chain"], start=2):

            print(
                f"    depth={i} "
                f"radix={result['radix_pairs'][i-2]} "
                f"x={st['x']} "
                f"y={st['y']} "
                f"u={st['u']} "
                f"v={st['v']} "
                f"E_in={st['E_in']} "
                f"E_next={st['E_next']}"
            )

        print()

    # --------------------------------------------------------------------------
    # Strongest contractions
    # --------------------------------------------------------------------------

    scored = []

    for anchor, result in results:

        stats = {
            st["level"]: st
            for st in result["level_stats"]
        }

        if 1 not in stats or 2 not in stats:
            continue

        s1 = stats[1]
        s2 = stats[2]

        if s1["unique_E"] == 0:
            continue

        ratio = s2["unique_E"] / s1["unique_E"]

        scored.append((ratio, anchor, result))

    scored.sort(key=lambda x: x[0])

    print("=" * 100)
    print("STRONGEST E-SET CONTRACTIONS")
    print("=" * 100)

    for ratio, anchor, result in scored[:20]:

        print(
            f"n={anchor['n']:,} "
            f"T1={result['T1']} "
            f"E1={len(result['E1_set'])} "
            f"E2={result['level_stats'][1]['unique_E'] if len(result['level_stats']) > 1 else 0} "
            f"ratio={ratio:.6f}"
        )

    # --------------------------------------------------------------------------
    # Weakest contractions
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("WEAKEST E-SET CONTRACTIONS")
    print("=" * 100)

    for ratio, anchor, result in scored[-20:]:

        print(
            f"n={anchor['n']:,} "
            f"T1={result['T1']} "
            f"E1={len(result['E1_set'])} "
            f"E2={result['level_stats'][1]['unique_E'] if len(result['level_stats']) > 1 else 0} "
            f"ratio={ratio:.6f}"
        )

    # --------------------------------------------------------------------------
    # Mathematical interpretation
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)
    print(
        """
The experiment tests whether the unknown correction term itself has
a recursive contraction law.

LEVEL 1:

    T1 = k*l + E1.

The observable uncertainty is represented by all admissible (k,l)
pairs generated from:

    K = T1-E1.

LEVEL 2:

    k = u*s1 + alpha
    l = v*s2 + beta

and therefore:

    T2 = floor((k*l+E1)/(s1*s2))
       = u*v + E2.

with:

    E2 =
        floor(
            u*s1*beta
          + v*s2*alpha
          + alpha*beta
          + E1
        / (s1*s2)
        ).

The old defect E1 is therefore explicitly carried into E2.

LEVEL 3 repeats the same operation.

The experiment measures the observable sets:

    S1 = possible E1 values
    S2 = possible E2 values
    S3 = possible E3 values
    ...

and also:

    P1 = possible k*l products
    P2 = possible u*v products
    ...

A positive result would look like:

    |S2| << |S1|
    |S3| << |S2|
    |S4| << |S3|

while the true branch remains present.

That would indicate genuine recursive contraction of the uncertainty.

An even stronger result would be:

    possible products
        ->
    very few descendant products
        ->
    tiny terminal coordinates.

A negative result is:

    defect-set size remains large,
    or shrinks only because the coordinates themselves were already
    enumerated,
    or the number of possible coordinate states explodes while the
    defect set becomes small.

The distinction is critical.

A small E value by itself is not enough.

We need observable contraction:

    UNKNOWN INPUT
        |
        v
       E1 set
        |
        v
       E2 set
        |
        v
       E3 set
        |
        v
    tiny terminal state

without inserting the true factor coordinates.

The true p,q values are used only to verify that the actual branch
survives the observable construction.

Every recursive identity is checked exactly.
"""
    )

    elapsed = time.perf_counter() - start

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime                 = {elapsed:.3f} seconds")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
