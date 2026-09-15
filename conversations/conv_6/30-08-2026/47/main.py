#!/usr/bin/env python3

import math
import random
import statistics
from dataclasses import dataclass
from collections import Counter


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 1_511_464_998

N_ANCHORS = 300

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

PRINT_EXAMPLES = 20


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Anchor:
    p: int
    q: int
    n: int

    r1: int
    r2: int

    a: int
    b: int

    k: int
    l: int

    T: int
    E: int

    c1: int
    d1: int

    c2: int
    d2: int

    c3: int

    g: int

    # Reversed orientation
    a2: int
    b2: int
    k2: int
    l2: int
    E2: int


# =============================================================================
# PRIME GENERATION
# =============================================================================

def build_primes(lo: int, hi: int):
    """
    Sieve of Eratosthenes.
    Much faster and simpler than trial-dividing every integer.
    """

    if hi < 2:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = math.isqrt(hi)

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x for x in range(lo, hi + 1)
        if sieve[x]
    ]


# =============================================================================
# CLOSE MODULUS PAIRS
# =============================================================================

def build_close_pairs(mod_primes):
    """
    Precompute only valid close modulus pairs.

    This is at most a few hundred thousand small tuples and is safe.
    """

    pairs = []

    for i, r1 in enumerate(mod_primes):
        max_r2 = int(math.floor(r1 * (1.0 + CLOSE_RATIO)))

        for j in range(i + 1, len(mod_primes)):

            r2 = mod_primes[j]

            if r2 > max_r2:
                break

            if r2 <= r1:
                continue

            if math.gcd(r1, r2) != 1:
                continue

            pairs.append((r1, r2))

    return pairs


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def build_anchors(factor_primes, close_pairs, count, rng):

    """
    Generate anchors directly.

    IMPORTANT:
    We do NOT construct all possible p*q combinations.

    This avoids the previous O(P^2) memory explosion.
    """

    anchors = []

    if not close_pairs:
        raise RuntimeError("No valid close modulus pairs found.")

    attempts = 0
    max_attempts = count * 100

    while len(anchors) < count and attempts < max_attempts:

        attempts += 1

        # Random factor primes.
        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        # Keep p <= q only for a stable orientation.
        if p > q:
            p, q = q, p

        n = p * q

        r1, r2 = rng.choice(close_pairs)

        # We want positive quotient coordinates.
        if p <= r1:
            continue

        if q <= r2:
            continue

        a = p % r1
        b = q % r2

        k = p // r1
        l = q // r2

        if k <= 0 or l <= 0:
            continue

        R = r1 * r2

        T = n // R
        E = T - k * l

        # ---------------------------------------------------------------------
        # Carry decomposition
        # ---------------------------------------------------------------------

        # a*l = c1*r1 + d1
        c1, d1 = divmod(a * l, r1)

        # b*k = c2*r2 + d2
        c2, d2 = divmod(b * k, r2)

        # Final carry.
        #
        # After replacing:
        #
        #   a*l = c1*r1 + d1
        #   b*k = c2*r2 + d2
        #
        # the remaining contribution is:
        #
        #   d1*r2 + d2*r1 + a*b
        #
        # divided by r1*r2.
        #
        c3 = (
            d1 * r2
            + d2 * r1
            + a * b
        ) // R

        g = n % R

        # ---------------------------------------------------------------------
        # Reversed orientation
        # ---------------------------------------------------------------------

        a2 = p % r2
        b2 = q % r1

        k2 = p // r2
        l2 = q // r1

        E2 = (n // R) - k2 * l2

        anchor = Anchor(
            p=p,
            q=q,
            n=n,

            r1=r1,
            r2=r2,

            a=a,
            b=b,

            k=k,
            l=l,

            T=T,
            E=E,

            c1=c1,
            d1=d1,

            c2=c2,
            d2=d2,

            c3=c3,

            g=g,

            a2=a2,
            b2=b2,

            k2=k2,
            l2=l2,

            E2=E2,
        )

        anchors.append(anchor)

    if len(anchors) < count:
        raise RuntimeError(
            f"Only generated {len(anchors)}/{count} valid anchors."
        )

    return anchors


# =============================================================================
# VALIDATION
# =============================================================================

def validate_anchor(x: Anchor):

    r1 = x.r1
    r2 = x.r2
    R = r1 * r2

    # Original coordinates.
    assert x.p == x.a + x.k * r1
    assert x.q == x.b + x.l * r2

    assert 0 <= x.a < r1
    assert 0 <= x.b < r2

    # Quotient-product relation.
    assert x.T == x.k * x.l + x.E

    # Original expansion.
    expanded = (
        x.k * x.l * R
        + x.a * x.l * r2
        + x.b * x.k * r1
        + x.a * x.b
    )

    assert expanded == x.n

    # Carry equations.
    assert x.a * x.l == x.c1 * r1 + x.d1
    assert x.b * x.k == x.c2 * r2 + x.d2

    assert 0 <= x.d1 < r1
    assert 0 <= x.d2 < r2

    # Exact E decomposition.
    assert x.E == x.c1 + x.c2 + x.c3

    # Final carry should be very small.
    assert x.c3 in (0, 1, 2)

    # Original E bound.
    assert 0 <= x.E <= x.k + x.l

    # Carry bounds.
    assert 0 <= x.c1 <= x.l - 1
    assert 0 <= x.c2 <= x.k - 1

    # Remainder.
    assert x.g == x.n % R


# =============================================================================
# CORRELATION
# =============================================================================

def correlation(xs, ys):

    if len(xs) != len(ys) or len(xs) < 2:
        return 0.0

    mx = statistics.mean(xs)
    my = statistics.mean(ys)

    numerator = sum(
        (x - mx) * (y - my)
        for x, y in zip(xs, ys)
    )

    dx = math.sqrt(
        sum((x - mx) ** 2 for x in xs)
    )

    dy = math.sqrt(
        sum((y - my) ** 2 for y in ys)
    )

    if dx == 0 or dy == 0:
        return 0.0

    return numerator / (dx * dy)


# =============================================================================
# SMALL INTEGER FACTORIZATION
# =============================================================================

def factor_integer(n):

    if n < 2:
        return []

    factors = []

    while n % 2 == 0:
        factors.append(2)
        n //= 2

    d = 3

    while d * d <= n:

        while n % d == 0:
            factors.append(d)
            n //= d

        d += 2

    if n > 1:
        factors.append(n)

    return factors


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def main():

    rng = random.Random(SEED)

    print("=" * 100)
    print("TWO-MODULUS E1 CARRY-GENEALOGY / CROSS-TERM DECOMPOSITION EXPERIMENT")
    print("=" * 100)

    print(f"N anchors                 = {N_ANCHORS:,}")
    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(
        f"modulus range             = "
        f"{MOD_MIN:,} - {MOD_MAX:,}"
    )
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # =========================================================================
    # PRIME POOLS
    # =========================================================================

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = build_primes(
        FACTOR_MIN,
        FACTOR_MAX,
    )

    mod_primes = build_primes(
        MOD_MIN,
        MOD_MAX,
    )

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(mod_primes):,}"
    )

    # =========================================================================
    # CLOSE PAIRS
    # =========================================================================

    print()

    print("=" * 100)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 100)

    close_pairs = build_close_pairs(mod_primes)

    print(
        f"close modulus pairs       = "
        f"{len(close_pairs):,}"
    )

    if not close_pairs:
        raise RuntimeError(
            "No close modulus pairs available."
        )

    # =========================================================================
    # ANCHORS
    # =========================================================================

    print()

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        close_pairs,
        N_ANCHORS,
        rng,
    )

    print(
        f"actual anchors            = "
        f"{len(anchors):,}"
    )

    # =========================================================================
    # VALIDATION
    # =========================================================================

    print()

    print("=" * 100)
    print("VALIDATING CARRY DECOMPOSITION")
    print("=" * 100)

    failures = 0

    for x in anchors:

        try:
            validate_anchor(x)
        except AssertionError:
            failures += 1

    print(
        f"carry identity failures   = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Carry decomposition validation failed."
        )

    # =========================================================================
    # BASIC ARRAYS
    # =========================================================================

    E = [x.E for x in anchors]

    c1 = [x.c1 for x in anchors]
    c2 = [x.c2 for x in anchors]
    c3 = [x.c3 for x in anchors]

    k = [x.k for x in anchors]
    l = [x.l for x in anchors]

    E2 = [x.E2 for x in anchors]

    # =========================================================================
    # E STRUCTURE
    # =========================================================================

    print()

    print("=" * 100)
    print("E1 CARRY STRUCTURE")
    print("=" * 100)

    print(
        f"average E                 = "
        f"{statistics.mean(E):.3f}"
    )

    print(
        f"minimum E                 = "
        f"{min(E)}"
    )

    print(
        f"maximum E                 = "
        f"{max(E)}"
    )

    print()

    print(
        f"average c1                = "
        f"{statistics.mean(c1):.3f}"
    )

    print(
        f"average c2                = "
        f"{statistics.mean(c2):.3f}"
    )

    print(
        f"average c3                = "
        f"{statistics.mean(c3):.3f}"
    )

    print(
        f"maximum c3                = "
        f"{max(c3)}"
    )

    # =========================================================================
    # EXACT DECOMPOSITION
    # =========================================================================

    print()

    print("=" * 100)
    print("EXACT E DECOMPOSITION")
    print("=" * 100)

    exact_count = sum(
        x.E == x.c1 + x.c2 + x.c3
        for x in anchors
    )

    print(
        f"E = c1+c2+c3             = "
        f"{exact_count}/{len(anchors)}"
    )

    # =========================================================================
    # C3
    # =========================================================================

    print()

    print("=" * 100)
    print("FINAL CARRY DISTRIBUTION")
    print("=" * 100)

    c3_counts = Counter(c3)

    for value in sorted(c3_counts):

        print(
            f"c3 = {value:2d} : "
            f"{c3_counts[value]:5d}"
        )

    # =========================================================================
    # PARTIAL E
    # =========================================================================

    print()

    print("=" * 100)
    print("PARTIAL E RECONSTRUCTION")
    print("=" * 100)

    partial = [
        x.E - x.c1 - x.c2
        for x in anchors
    ]

    partial_counts = Counter(partial)

    for value in sorted(partial_counts):

        print(
            f"E-(c1+c2) = {value:3d} : "
            f"{partial_counts[value]:5d}"
        )

    # =========================================================================
    # CARRY RATIOS
    # =========================================================================

    print()

    print("=" * 100)
    print("E COMPONENT RATIOS")
    print("=" * 100)

    nonzero_E = [
        x for x in anchors
        if x.E > 0
    ]

    if nonzero_E:

        c1_ratio = statistics.mean(
            x.c1 / x.E
            for x in nonzero_E
        )

        c2_ratio = statistics.mean(
            x.c2 / x.E
            for x in nonzero_E
        )

        c3_ratio = statistics.mean(
            x.c3 / x.E
            for x in nonzero_E
        )

        print(
            f"mean c1/E               = "
            f"{c1_ratio:.6f}"
        )

        print(
            f"mean c2/E               = "
            f"{c2_ratio:.6f}"
        )

        print(
            f"mean c3/E               = "
            f"{c3_ratio:.6f}"
        )

    # =========================================================================
    # SCALE
    # =========================================================================

    print()

    print("=" * 100)
    print("CARRY SCALE")
    print("=" * 100)

    print(
        f"average k                = "
        f"{statistics.mean(k):.3f}"
    )

    print(
        f"average l                = "
        f"{statistics.mean(l):.3f}"
    )

    print(
        f"average c1/l             = "
        f"{statistics.mean(x.c1 / max(1, x.l) for x in anchors):.6f}"
    )

    print(
        f"average c2/k             = "
        f"{statistics.mean(x.c2 / max(1, x.k) for x in anchors):.6f}"
    )

    # =========================================================================
    # CORRELATIONS
    # =========================================================================

    print()

    print("=" * 100)
    print("CORRELATION ANALYSIS")
    print("=" * 100)

    print(
        f"corr(E,c1)               = "
        f"{correlation(E, c1):.6f}"
    )

    print(
        f"corr(E,c2)               = "
        f"{correlation(E, c2):.6f}"
    )

    print(
        f"corr(E,c3)               = "
        f"{correlation(E, c3):.6f}"
    )

    print(
        f"corr(E,k)                = "
        f"{correlation(E, k):.6f}"
    )

    print(
        f"corr(E,l)                = "
        f"{correlation(E, l):.6f}"
    )

    # =========================================================================
    # NEW:
    # DOES E FOLLOW ONE OF THE CROSS TERMS?
    # =========================================================================

    print()

    print("=" * 100)
    print("CROSS-TERM RELATIONS")
    print("=" * 100)

    cross1 = [
        x.a * x.l
        for x in anchors
    ]

    cross2 = [
        x.b * x.k
        for x in anchors
    ]

    print(
        f"average a*l              = "
        f"{statistics.mean(cross1):.3f}"
    )

    print(
        f"average b*k              = "
        f"{statistics.mean(cross2):.3f}"
    )

    print(
        f"corr(E,a*l)              = "
        f"{correlation(E, cross1):.6f}"
    )

    print(
        f"corr(E,b*k)              = "
        f"{correlation(E, cross2):.6f}"
    )

    # =========================================================================
    # GCD STRUCTURE
    # =========================================================================

    print()

    print("=" * 100)
    print("GCD STRUCTURE")
    print("=" * 100)

    gcd_c1_k = [
        math.gcd(x.c1, x.k)
        for x in anchors
    ]

    gcd_c2_l = [
        math.gcd(x.c2, x.l)
        for x in anchors
    ]

    print(
        f"average gcd(c1,k)       = "
        f"{statistics.mean(gcd_c1_k):.3f}"
    )

    print(
        f"average gcd(c2,l)       = "
        f"{statistics.mean(gcd_c2_l):.3f}"
    )

    print(
        f"c1 shares factor with k  = "
        f"{sum(g > 1 for g in gcd_c1_k)}/{len(anchors)}"
    )

    print(
        f"c2 shares factor with l  = "
        f"{sum(g > 1 for g in gcd_c2_l)}/{len(anchors)}"
    )

    # =========================================================================
    # INNER CARRY FACTORIZATION
    # =========================================================================

    print()

    print("=" * 100)
    print("INNER CARRY FACTORIZATION")
    print("=" * 100)

    c1_factor_counts = []
    c2_factor_counts = []

    for x in anchors:

        if x.c1 > 1:
            c1_factor_counts.append(
                len(factor_integer(x.c1))
            )

        if x.c2 > 1:
            c2_factor_counts.append(
                len(factor_integer(x.c2))
            )

    if c1_factor_counts:

        print(
            f"mean prime factors(c1) = "
            f"{statistics.mean(c1_factor_counts):.3f}"
        )

    else:

        print(
            "mean prime factors(c1) = 0"
        )

    if c2_factor_counts:

        print(
            f"mean prime factors(c2) = "
            f"{statistics.mean(c2_factor_counts):.3f}"
        )

    else:

        print(
            "mean prime factors(c2) = 0"
        )

    # =========================================================================
    # REVERSED ORIENTATION
    # =========================================================================

    print()

    print("=" * 100)
    print("REVERSED r1/r2 ORIENTATION")
    print("=" * 100)

    print(
        f"average E'               = "
        f"{statistics.mean(E2):.3f}"
    )

    delta_E = [
        x.E - x.E2
        for x in anchors
    ]

    print(
        f"average |E-E'|           = "
        f"{statistics.mean(abs(v) for v in delta_E):.3f}"
    )

    print(
        f"maximum |E-E'|           = "
        f"{max(abs(v) for v in delta_E)}"
    )

    print(
        f"corr(E,E')               = "
        f"{correlation(E, E2):.6f}"
    )

    # =========================================================================
    # EXAMPLES
    # =========================================================================

    print()

    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for x in anchors[:PRINT_EXAMPLES]:

        print(
            f"n={x.n:,} "
            f"p={x.p:,} "
            f"q={x.q:,} "
            f"mods=({x.r1},{x.r2})"
        )

        print(
            f"    a={x.a} "
            f"b={x.b} "
            f"k={x.k} "
            f"l={x.l}"
        )

        print(
            f"    T={x.T} "
            f"E={x.E}"
        )

        print(
            f"    a*l={x.a*x.l} "
            f"= {x.c1}*{x.r1}+{x.d1}"
        )

        print(
            f"    b*k={x.b*x.k} "
            f"= {x.c2}*{x.r2}+{x.d2}"
        )

        print(
            f"    E={x.E} "
            f"= {x.c1}+{x.c2}+{x.c3}"
        )

        print(
            f"    c3={x.c3} "
            f"g=n mod (r1*r2)={x.g}"
        )

        print(
            f"    reversed: "
            f"k'={x.k2} "
            f"l'={x.l2} "
            f"E'={x.E2}"
        )

        print()

    # =========================================================================
    # MATHEMATICAL INTERPRETATION
    # =========================================================================

    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
For

    p = a + k*r1
    q = b + l*r2

the exact product expansion is

    n =
        k*l*r1*r2
        + a*l*r2
        + b*k*r1
        + a*b.

Let

    T = floor(n/(r1*r2))

and

    E = T-k*l.

The cross terms can be decomposed as

    a*l = c1*r1 + d1
    b*k = c2*r2 + d2.

This produces the exact identity

    E = c1 + c2 + c3,

where

    c3 =
        floor(
            (d1*r2 + d2*r1 + a*b)
            /(r1*r2)
        ).

Because

    0 <= d1 < r1
    0 <= d2 < r2
    0 <= a < r1
    0 <= b < r2,

the final carry c3 is very small.

The research question is therefore no longer merely:

    "Does E have structure?"

It definitely does.

The interesting question is:

    "Does the structure of E recursively expose information about
     k and l?"

In particular, c1 and c2 satisfy

    c1 = floor(a*l/r1)
    c2 = floor(b*k/r2).

These are themselves quotient variables generated by products.

That creates a possible second-level recursive structure:

    E1
     |
     +---- c1 = floor(a*l/r1)
     |
     +---- c2 = floor(b*k/r2)
     |
     +---- c3 = final carry

The experiment therefore measures whether c1 and c2 have additional
relationships with k,l, their factorization, gcd structure, or the
reversed modulus orientation.

The key point is that the identity itself is not a discovery.

The discovery would be a stronger invariant such as:

    c1 = f(k,l,r1,r2)

with a substantially narrower range than the trivial bound,

or a recursive relationship in which c1/c2 can be constrained
without knowing a,b.

That would move the research from:

    exact decomposition

toward:

    predictive decomposition.

This distinction is essential.

The script never constructs all possible factor pairs. Anchors are
generated directly, so memory usage remains small.
"""
    )

    print()

    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()