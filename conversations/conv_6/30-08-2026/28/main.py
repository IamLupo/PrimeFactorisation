#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME QUOTIENT-RESULTANT / COORDINATE-ELIMINATION EXPERIMENT
====================================================================================================

Goal
----
Test whether the closeness

    r2 = r1 + d2
    r3 = r1 + d3

creates an exploitable algebraic relation between the quotient coordinates of
the unknown factors, without enumerating the full factor-prime population.

For

    p = a + k*r1
    q = b + l*r1

we have

    a*b == g (mod r1)

where

    g = n mod r1.

For any second modulus r = r1 + d:

    p mod r = a - k*d (mod r)
    q mod r = b - l*d (mod r)

and because pq == g (mod r),

    (a-kd)(b-ld) == g (mod r).

Subtracting the r1 equation gives:

    d * (a*l + b*k - d*k*l) == 0 (mod r)

because r is prime and d != 0 mod r,

    a*l + b*k - d*k*l == 0 (mod r).

Define

    T = a*l + b*k
    K = k*l.

Then for r2,r3:

    T == d2*K (mod r2)
    T == d3*K (mod r3).

Thus the pair of close moduli imposes two linear congruences on the SAME
quantity T.

This experiment:

    1. enumerates quotient coordinates only,
    2. derives T from K using CRT,
    3. exploits the small physical range of T,
    4. reconstructs possible (a,b),
    5. performs exact verification.

The important measurement is whether the CRT-derived T state collapses
the coordinate search substantially BEFORE exact factor verification.

This is an experimental search, not a claim of a general factoring algorithm.
====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ANCHORS = 300

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

# Only report these many detailed anchors.
DETAIL_ANCHORS = 20

# Progress interval.
PROGRESS = 25


# ================================================================================================
# BASIC MATH
# ================================================================================================

def sieve(limit: int) -> list[int]:
    """
    Reliable bytearray sieve.

    Avoids extended-slice size issues by assigning exact slice lengths.
    """
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start::p] = b"\x00" * count

    return [i for i, v in enumerate(a) if v]


def egcd(a: int, b: int):
    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def invmod(a: int, m: int) -> int:
    a %= m
    g, x, _ = egcd(a, m)

    if g != 1:
        raise ValueError(f"{a} has no inverse modulo {m}")

    return x % m


def crt_pair(a1: int, m1: int, a2: int, m2: int) -> int:
    """
    Return x in [0,m1*m2) satisfying

        x == a1 mod m1
        x == a2 mod m2

    Assumes gcd(m1,m2)=1.
    """
    t = ((a2 - a1) * invmod(m1, m2)) % m2
    return a1 + m1 * t


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


# ================================================================================================
# CLOSE TRIPLE SELECTION
# ================================================================================================

def select_close_triples(
    factor_primes: list[int],
    modulus_primes: list[int],
    anchors: list[tuple[int, int]],
    rng: random.Random,
) -> list[tuple[int, tuple[int, int, int]]]:
    """
    Select one maximum-product close triple for every anchor.

    The anchor controls how close the three moduli need to be relative to
    the candidate modulus population.
    """
    selected = []

    for idx, (p, q) in enumerate(anchors, 1):
        n = p * q

        usable = []

        for i in range(len(modulus_primes)):
            r1 = modulus_primes[i]

            upper = int(r1 * (1.0 + CLOSE_RATIO))

            # Search candidates after r1.
            for j in range(i + 1, len(modulus_primes)):
                r2 = modulus_primes[j]

                if r2 > upper:
                    break

                for k in range(j + 1, len(modulus_primes)):
                    r3 = modulus_primes[k]

                    if r3 > upper:
                        break

                    R = r1 * r2 * r3

                    if R < n:
                        usable.append((R, r1, r2, r3))

        if not usable:
            continue

        usable.sort(reverse=True)
        _, r1, r2, r3 = usable[0]

        selected.append((n, (r1, r2, r3)))

        if idx % PROGRESS == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    return selected


# ================================================================================================
# ANCHOR GENERATION
# ================================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    rng: random.Random,
) -> list[tuple[int, int]]:
    """
    Random distinct prime pairs from the factor interval.
    """
    pairs = []

    for _ in range(count):
        p, q = rng.sample(factor_primes, 2)
        if p > q:
            p, q = q, p

        pairs.append((p, q))

    return pairs


# ================================================================================================
# QUOTIENT DOMAIN
# ================================================================================================

def quotient_range(value_min: int, value_max: int, r: int) -> range:
    """
    k such that

        value_min <= a + k*r <= value_max

    for some residue a in [0,r).
    """
    # Conservative quotient range.
    k_min = max(0, value_min // r - 1)
    k_max = value_max // r + 1

    return range(k_min, k_max + 1)


# ================================================================================================
# CORE RESULTANT EXPERIMENT
# ================================================================================================

def solve_anchor(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
):
    r1, r2, r3 = mods

    assert r1 < r2 < r3

    d2 = r2 - r1
    d3 = r3 - r1

    g = n % r1

    # --------------------------------------------------------------------------------------------
    # Precompute b = g * a^{-1} mod r1.
    #
    # For prime r1 and gcd(g,r1)=1 this exists for every nonzero a.
    # a==0 can occur only when r1 divides p. Since p is prime and factor-prime
    # interval is much larger than r1, this is possible only if p == r1,
    # which cannot happen for the selected triples here. Still handled safely.
    # --------------------------------------------------------------------------------------------

    residue_pairs = []

    for a in range(r1):
        if a == 0:
            continue

        b = (g * invmod(a, r1)) % r1

        # Verify the base modular hyperbola.
        if (a * b - g) % r1 != 0:
            raise AssertionError("r1 residue relation failed")

        residue_pairs.append((a, b))

    # --------------------------------------------------------------------------------------------
    # Build possible quotient coordinates.
    #
    # p = a + k*r1
    # q = b + l*r1
    #
    # We exploit:
    #
    #   T = a*l + b*k
    #   K = k*l
    #
    # and
    #
    #   T == d2*K (mod r2)
    #   T == d3*K (mod r3).
    #
    # Therefore for every K we obtain one CRT class for T.
    #
    # Rather than enumerate every (a,k,l), we enumerate feasible K values
    # and inspect only T representatives that can physically fit inside
    # the quotient-coordinate bounds.
    # --------------------------------------------------------------------------------------------

    k_values = list(quotient_range(FACTOR_MIN, FACTOR_MAX, r1))
    l_values = k_values[:]  # same factor interval

    max_abs_T = 2 * r1 * (max(k_values) + max(l_values) + 2)

    R23 = r2 * r3

    # Because T is small while R23 is usually much larger, only a few CRT
    # representatives can lie in the physical T interval.
    #
    # We use a moderate K domain, derived from k*l.
    max_K = max(k_values) * max(l_values)

    candidate_T_by_K = {}

    states_K = 0
    states_T = 0

    for k in k_values:
        for l in l_values:
            K = k * l

            if K in candidate_T_by_K:
                continue

            states_K += 1

            t2 = (d2 * K) % r2
            t3 = (d3 * K) % r3

            T0 = crt_pair(t2, r2, t3, r3)

            # T can differ by R23.
            #
            # For our coordinate representation the physically useful T is
            # approximately bounded by |a*l+b*k| < 2*r1*(k+l).
            #
            # Usually R23 is so large that there is at most one possible
            # representative.
            m_min = math.ceil((-max_abs_T - T0) / R23)
            m_max = math.floor((max_abs_T - T0) / R23)

            reps = []

            for m in range(m_min, m_max + 1):
                T = T0 + m * R23

                if -max_abs_T <= T <= max_abs_T:
                    reps.append(T)

            candidate_T_by_K[K] = reps
            states_T += len(reps)

    # --------------------------------------------------------------------------------------------
    # Actual coordinate test.
    #
    # Now examine (a,b) states, but eliminate l using:
    #
    #   T = a*l + b*k
    #
    # Therefore
    #
    #   l = (T - b*k) / a
    #
    # whenever a != 0.
    #
    # This avoids explicitly enumerating l for every (a,k).
    # --------------------------------------------------------------------------------------------

    exact_pairs = set()
    prime_pairs = set()

    surviving_coordinate_states = 0
    divisibility_states = 0
    physical_states = 0

    actual_state = None

    for a, b in residue_pairs:

        for k in k_values:
            p = a + k * r1

            if not (FACTOR_MIN <= p <= FACTOR_MAX):
                continue

            # p need not be prime at this stage.
            # We intentionally delay the primality test.

            for T_values in candidate_T_by_K.values():
                # We do not know K directly here.
                # This branch is intentionally replaced below by a targeted
                # K loop. Keeping the algebra explicit is useful for auditing.
                pass

            # Targeted enumeration over l is avoided by using K buckets.
            #
            # For fixed k we have K = k*l.
            if k == 0:
                continue

            for l in l_values:
                q = b + l * r1

                if not (FACTOR_MIN <= q <= FACTOR_MAX):
                    continue

                K = k * l
                Ts = candidate_T_by_K.get(K)

                if not Ts:
                    continue

                T_actual = a * l + b * k

                if T_actual not in Ts:
                    continue

                surviving_coordinate_states += 1

                # Third congruence is already implicit in the CRT T construction,
                # but verify it independently.
                if (T_actual - d2 * K) % r2 != 0:
                    continue

                if (T_actual - d3 * K) % r3 != 0:
                    continue

                divisibility_states += 1

                # Exact factor pair.
                if p * q != n:
                    continue

                physical_states += 1

                pp, qq = sorted((p, q))

                exact_pairs.add((pp, qq))

                if is_prime(pp) and is_prime(qq):
                    prime_pairs.add((pp, qq))

                if {pp, qq} == {p_true, q_true}:
                    actual_state = (a, b, k, l, T_actual, K)

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "g": g,
        "r23": R23,
        "k_domain": len(k_values),
        "a_domain": len(residue_pairs),
        "K_states": states_K,
        "T_states": states_T,
        "surviving_coordinate_states": surviving_coordinate_states,
        "divisibility_states": divisibility_states,
        "physical_states": physical_states,
        "exact_pairs": sorted(exact_pairs),
        "prime_pairs": sorted(prime_pairs),
        "actual_state": actual_state,
    }


# ================================================================================================
# FAST AUDIT VERSION
# ================================================================================================

def fast_anchor_audit(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
):
    """
    Much faster independent audit.

    Instead of rebuilding all coordinate combinations, directly verify the
    resultant identity for the true factor pair and measure the size of the
    quotient-coordinate domain.
    """

    r1, r2, r3 = mods

    d2 = r2 - r1
    d3 = r3 - r1

    p = p_true
    q = q_true

    a = p % r1
    b = q % r1

    k = (p - a) // r1
    l = (q - b) // r1

    K = k * l
    T = a * l + b * k

    ok_r1 = (a * b) % r1 == n % r1
    ok_r2 = (T - d2 * K) % r2 == 0
    ok_r3 = (T - d3 * K) % r3 == 0

    assert ok_r1
    assert ok_r2
    assert ok_r3

    return {
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "ok": ok_r1 and ok_r2 and ok_r3,
        "R23": r2 * r3,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def run():
    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME QUOTIENT-RESULTANT / COORDINATE-ELIMINATION EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # --------------------------------------------------------------------------------------------
    # PRIME POOLS
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    all_primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in all_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # --------------------------------------------------------------------------------------------
    # ANCHORS
    # --------------------------------------------------------------------------------------------

    anchors = build_anchors(factor_primes, ANCHORS, rng)

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    selected = select_close_triples(
        factor_primes,
        modulus_primes,
        anchors,
        rng,
    )

    print()
    print(f"usable anchors            = {len(selected)}")
    print()

    # --------------------------------------------------------------------------------------------
    # RESULT COLLECTION
    # --------------------------------------------------------------------------------------------

    results = []

    K_counts = []
    T_counts = []
    surviving_counts = []
    divisibility_counts = []
    exact_counts = []

    recovery = 0
    prime_recovery = 0

    t0 = time.perf_counter()

    print("=" * 100)
    print("RUNNING QUOTIENT-RESULTANT ANALYSIS")
    print("=" * 100)

    for idx, (n, mods) in enumerate(selected, 1):

        p_true, q_true = anchors[idx - 1]

        # Make sure the selected anchor actually matches n.
        if p_true * q_true != n:
            raise AssertionError("anchor mismatch")

        audit = fast_anchor_audit(
            n,
            p_true,
            q_true,
            mods,
        )

        if not audit["ok"]:
            raise AssertionError("true resultant identity failed")

        # --------------------------------------------------------------------
        # Full coordinate experiment.
        #
        # To keep runtime bounded, execute the expensive search only for the
        # first DETAIL_ANCHORS anchors.
        #
        # The remaining anchors still undergo the independent identity audit.
        # --------------------------------------------------------------------

        if idx <= DETAIL_ANCHORS:
            result = solve_anchor(
                n,
                p_true,
                q_true,
                mods,
            )

            results.append(result)

            K_counts.append(result["K_states"])
            T_counts.append(result["T_states"])
            surviving_counts.append(
                result["surviving_coordinate_states"]
            )
            divisibility_counts.append(
                result["divisibility_states"]
            )
            exact_counts.append(
                result["physical_states"]
            )

            if (p_true, q_true) in result["prime_pairs"]:
                prime_recovery += 1

            if result["actual_state"] is not None:
                recovery += 1

        if idx % PROGRESS == 0:
            print(f"anchor {idx:3d}/{len(selected)}")

    elapsed = time.perf_counter() - t0

    # --------------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("IDENTITY AUDIT")
    print("=" * 100)

    print(f"anchors identity-verified  = {len(selected)}/{len(selected)}")
    print()

    if results:
        print("=" * 100)
        print("RESULTANT SEARCH SUMMARY")
        print("=" * 100)

        print(
            f"full coordinate anchors    = {len(results)}"
        )
        print(
            f"mean K states              = {statistics.mean(K_counts):.3f}"
        )
        print(
            f"mean T states              = {statistics.mean(T_counts):.3f}"
        )
        print(
            f"mean surviving states      = "
            f"{statistics.mean(surviving_counts):.3f}"
        )
        print(
            f"mean divisor-valid states  = "
            f"{statistics.mean(divisibility_counts):.3f}"
        )
        print(
            f"mean exact pair states     = "
            f"{statistics.mean(exact_counts):.3f}"
        )

        print()
        print(
            f"minimum T states           = {min(T_counts):,}"
        )
        print(
            f"maximum T states           = {max(T_counts):,}"
        )
        print(
            f"minimum surviving states   = {min(surviving_counts):,}"
        )
        print(
            f"maximum surviving states   = {max(surviving_counts):,}"
        )

        print()
        print("=" * 100)
        print("RECOVERY")
        print("=" * 100)

        print(
            f"resultant state recovery   = "
            f"{recovery}/{len(results)}"
        )
        print(
            f"prime-pair recovery        = "
            f"{prime_recovery}/{len(results)}"
        )

        print()
        print("=" * 100)
        print("CRT T-STATE DISTRIBUTION")
        print("=" * 100)

        dist = Counter(T_counts)

        for value, count in sorted(dist.items()):
            print(
                f"T states = {value:6d} anchors = {count:3d}"
            )

        print()
        print("=" * 100)
        print("STRONGEST COLLAPSES")
        print("=" * 100)

        ranked = sorted(
            results,
            key=lambda x: (
                x["surviving_coordinate_states"],
                x["T_states"],
            ),
        )

        for r in ranked[:20]:
            print(
                f"mods=({r['r1']},{r['r2']},{r['r3']}) "
                f"R23={r['r23']:,} "
                f"K={r['K_states']:,} "
                f"T={r['T_states']:,} "
                f"survive={r['surviving_coordinate_states']:,} "
                f"div={r['divisibility_states']:,} "
                f"exact={r['physical_states']:,}"
            )

        print()
        print("=" * 100)
        print("WEAKEST COLLAPSES")
        print("=" * 100)

        for r in ranked[-20:]:
            print(
                f"mods=({r['r1']},{r['r2']},{r['r3']}) "
                f"R23={r['r23']:,} "
                f"K={r['K_states']:,} "
                f"T={r['T_states']:,} "
                f"survive={r['surviving_coordinate_states']:,} "
                f"div={r['divisibility_states']:,} "
                f"exact={r['physical_states']:,}"
            )

        print()
        print("=" * 100)
        print("ACTUAL QUOTIENT-RESULTANT STATES")
        print("=" * 100)

        for i, r in enumerate(results[:20], 1):
            state = r["actual_state"]

            print(
                f"{i:3d}: "
                f"mods=({r['r1']},{r['r2']},{r['r3']}) "
                f"a={state[0]} "
                f"b={state[1]} "
                f"k={state[2]} "
                f"l={state[3]} "
                f"K={state[5]} "
                f"T={state[4]}"
            )

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"full-search timing           = {elapsed:.3f} seconds")

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
For

    p = a + k*r1
    q = b + l*r1

the first modulus gives

    a*b == g (mod r1).

Let

    r2 = r1 + d2
    r3 = r1 + d3.

The second and third moduli imply

    T = a*l + b*k
    K = k*l

and therefore

    T == d2*K (mod r2)
    T == d3*K (mod r3).

So for each K we can CRT-reconstruct T.

Because T is constrained by the physical quotient ranges, the CRT class
may have very few representatives.

The hoped-for chain is:

    quotient coordinates
          |
          v
        K = kl
          |
          v
      CRT(T)
          |
          v
    coordinate consistency
          |
          v
       p,q interval
          |
          v
      exact p*q = n.

The critical measurements are:

    K states
    T states
    surviving coordinate states
    exact pair states.

The important comparison is not whether the identity works — it necessarily
does for the true factors.

The question is whether the close-modulus resultant condition reduces the
unknown coordinate search BEFORE exact multiplication/divisibility testing.

A strong result would look like:

    huge coordinate domain
          ->
    very small T-state domain
          ->
    very few coordinate states
          ->
    one exact factor pair.

A weak result means the resultant condition is simply another encoding of
the original factor equation.

NOTE
----
Only DETAIL_ANCHORS receive the full coordinate search so that the script
finishes quickly. The identity itself is independently verified on every
anchor.

To scale the experiment, increase DETAIL_ANCHORS after measuring its runtime.
"""
    )

    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
