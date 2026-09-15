#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME EXACT QUOTIENT / LINEARIZED SECOND-MODULUS EXPERIMENT
====================================================================================================

This is the corrected version of the previous quotient-resultant experiment.

IMPORTANT CORRECTION
--------------------

We know:

    a*b == g (mod r1)

but NOT:

    a*b == g (mod r2)

Instead define the exact integer quotient

    h = (a*b - g) / r1.

With

    p = a + k*r1
    q = b + l*r1
    r2 = r1 + Delta

we have

    p mod r2 = a - k*Delta
    q mod r2 = b - l*Delta.

The valid r2 condition is therefore:

    (a-k*Delta)(b-l*Delta) == g (mod r2).

Expanding:

    a*b
    - a*l*Delta
    - b*k*Delta
    + k*l*Delta^2
    - g
    == 0 (mod r2).

Using

    a*b - g = h*r1,

we obtain:

    h*r1
    - a*l*Delta
    - b*k*Delta
    + k*l*Delta^2
    == 0 (mod r2).

For fixed a,b,k:

    l * (k*Delta^2 - a*Delta)
    ==
    b*k*Delta - h*r1
    (mod r2).

Thus l can be solved directly modulo r2 whenever the coefficient
is invertible.

The third modulus is then used as an independent verification.

This experiment measures whether this exact quotient treatment gives
a useful reduction without making the invalid substitution from the
previous experiment.

====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter


# ================================================================================================
# CONFIG
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

ANCHORS = 300

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS = 25

# Number of anchors receiving the full exact quotient search.
# Raise after checking runtime.
FULL_SEARCH_ANCHORS = 20


# ================================================================================================
# PRIME SIEVE
# ================================================================================================

def sieve(limit: int) -> list[int]:
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

    return [i for i, x in enumerate(a) if x]


# ================================================================================================
# MODULAR INVERSE
# ================================================================================================

def egcd(a: int, b: int):
    if b == 0:
        return b + a, 1, 0

    g, x1, y1 = egcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def invmod(a: int, m: int) -> int:
    a %= m

    g, x, _ = egcd(a, m)

    if g != 1:
        raise ValueError("not invertible")

    return x % m


# ================================================================================================
# PRIME TEST
# ================================================================================================

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
# ANCHOR GENERATION
# ================================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    rng: random.Random,
) -> list[tuple[int, int]]:

    result = []

    for _ in range(count):
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        result.append((p, q))

    return result


# ================================================================================================
# CLOSE MODULUS TRIPLE
# ================================================================================================

def choose_close_triple(
    n: int,
    modulus_primes: list[int],
) -> tuple[int, int, int]:

    best = None
    best_R = -1

    for i, r1 in enumerate(modulus_primes):

        upper = int(r1 * (1.0 + CLOSE_RATIO))

        for j in range(i + 1, len(modulus_primes)):

            r2 = modulus_primes[j]

            if r2 > upper:
                break

            for k in range(j + 1, len(modulus_primes)):

                r3 = modulus_primes[k]

                if r3 > upper:
                    break

                R = r1 * r2 * r3

                if R < n and R > best_R:
                    best_R = R
                    best = (r1, r2, r3)

    if best is None:
        raise RuntimeError("no usable modulus triple")

    return best


# ================================================================================================
# VALIDATE TRUE FACTOR STATE
# ================================================================================================

def validate_true_state(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    r3: int,
):
    g1 = n % r1
    g2 = n % r2
    g3 = n % r3

    assert g1 == n % r1
    assert g2 == n % r2
    assert g3 == n % r3

    a = p % r1
    b = q % r1

    k = (p - a) // r1
    l = (q - b) // r1

    h_num = a * b - g1

    assert h_num % r1 == 0

    h = h_num // r1

    delta2 = r2 - r1
    delta3 = r3 - r1

    # r1 relation
    assert (a * b - g1) % r1 == 0

    # Exact r2 equation.
    lhs2 = (
        h * r1
        - a * l * delta2
        - b * k * delta2
        + k * l * delta2 * delta2
    )

    assert lhs2 % r2 == 0

    # Exact r3 equation.
    lhs3 = (
        h * r1
        - a * l * delta3
        - b * k * delta3
        + k * l * delta3 * delta3
    )

    assert lhs3 % r3 == 0

    # Direct modular check.
    assert (
        ((a - k * delta2) * (b - l * delta2) - n)
        % r2
        == 0
    )

    assert (
        ((a - k * delta3) * (b - l * delta3) - n)
        % r3
        == 0
    )

    return a, b, k, l, h


# ================================================================================================
# SOLVE ONE ANCHOR
# ================================================================================================

def solve_anchor(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
):
    r1, r2, r3 = mods

    delta2 = r2 - r1
    delta3 = r3 - r1

    g = n % r1

    # ------------------------------------------------------------------------
    # Quotient domains.
    # ------------------------------------------------------------------------

    k_min = max(0, (FACTOR_MIN - (r1 - 1) + r1 - 1) // r1)
    k_max = FACTOR_MAX // r1

    k_values = range(k_min, k_max + 1)

    # ------------------------------------------------------------------------
    # Residue hyperbola at r1:
    #
    #     a*b == g (mod r1)
    #
    # Only invertible residues are needed here.
    # ------------------------------------------------------------------------

    residue_pairs = []

    for a in range(1, r1):

        b = (g * pow(a, -1, r1)) % r1

        # Exact quotient.
        h_num = a * b - g

        if h_num % r1 != 0:
            raise AssertionError(
                "r1 residue quotient is not integral"
            )

        h = h_num // r1

        residue_pairs.append((a, b, h))

    candidate_count = 0
    invertible_count = 0
    interval_count = 0
    r3_count = 0
    prime_count = 0
    exact_count = 0

    exact_pairs = set()
    prime_pairs = set()

    actual_found = False
    actual_l = None

    # ------------------------------------------------------------------------
    # Search.
    #
    # For each (a,b,k) solve l modulo r2 directly.
    #
    #     l * C == RHS (mod r2)
    #
    # with
    #
    #     C = k*delta2^2 - a*delta2
    #
    #     RHS = b*k*delta2 - h*r1.
    #
    # ------------------------------------------------------------------------

    for a, b, h in residue_pairs:

        for k in k_values:

            p = a + k * r1

            if not (FACTOR_MIN <= p <= FACTOR_MAX):
                continue

            C = (
                k * delta2 * delta2
                - a * delta2
            ) % r2

            RHS = (
                b * k * delta2
                - h * r1
            ) % r2

            # r2 is prime, therefore C has either one inverse or is zero.
            if C == 0:

                # Degenerate linear congruence.
                if RHS != 0:
                    continue

                # Every l modulo r2 is possible.
                #
                # We cannot exploit this case directly, so enumerate the
                # small quotient domain.
                for l in range(k_min, k_max + 1):

                    candidate_count += 1

                    q = b + l * r1

                    if not (FACTOR_MIN <= q <= FACTOR_MAX):
                        continue

                    interval_count += 1

                    # Independent r3 verification.
                    x3 = (a - k * delta3) % r3
                    y3 = (b - l * delta3) % r3

                    if (x3 * y3 - g) % r3 != 0:
                        continue

                    r3_count += 1

                    if p * q != n:
                        continue

                    exact_count += 1

                    pp, qq = sorted((p, q))
                    exact_pairs.add((pp, qq))

                    if is_prime(pp) and is_prime(qq):
                        prime_count += 1
                        prime_pairs.add((pp, qq))

                    if {pp, qq} == {p_true, q_true}:
                        actual_found = True
                        actual_l = l

                continue

            invertible_count += 1

            l0 = (RHS * pow(C, -1, r2)) % r2

            # Because r2 is around 300..3000 and the quotient k/l range is
            # tiny, there are normally only a few l = l0 + t*r2 to inspect.
            #
            # Find the allowed quotient range.
            #
            # q = b + l*r1
            #
            # FACTOR_MIN <= q <= FACTOR_MAX
            #

            l_low = max(
                0,
                (FACTOR_MIN - b + r1 - 1) // r1,
            )

            l_high = (
                FACTOR_MAX - b
            ) // r1

            # Enumerate only the arithmetic progression l0 mod r2.
            #
            # This is the main elimination step.
            if l0 < l_low:
                t0 = (l_low - l0 + r2 - 1) // r2
                l = l0 + t0 * r2
            else:
                l = l0

            while l <= l_high:

                candidate_count += 1

                q = b + l * r1

                if not (FACTOR_MIN <= q <= FACTOR_MAX):
                    l += r2
                    continue

                interval_count += 1

                # ------------------------------------------------------------
                # Independent third-modulus consistency check.
                # ------------------------------------------------------------

                x3 = (a - k * delta3) % r3
                y3 = (b - l * delta3) % r3

                if (x3 * y3 - g) % r3 != 0:
                    l += r2
                    continue

                r3_count += 1

                # ------------------------------------------------------------
                # Exact product.
                # ------------------------------------------------------------

                if p * q == n:

                    exact_count += 1

                    pp, qq = sorted((p, q))

                    exact_pairs.add((pp, qq))

                    if is_prime(pp) and is_prime(qq):
                        prime_count += 1
                        prime_pairs.add((pp, qq))

                    if {pp, qq} == {p_true, q_true}:
                        actual_found = True
                        actual_l = l

                l += r2

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "delta2": delta2,
        "delta3": delta3,
        "r2r3": r2 * r3,
        "residue_states": len(residue_pairs),
        "k_states": len(k_values),
        "candidate_count": candidate_count,
        "invertible_count": invertible_count,
        "interval_count": interval_count,
        "r3_count": r3_count,
        "prime_count": prime_count,
        "exact_count": exact_count,
        "exact_pairs": sorted(exact_pairs),
        "prime_pairs": sorted(prime_pairs),
        "actual_found": actual_found,
        "actual_l": actual_l,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME EXACT QUOTIENT / LINEARIZED SECOND-MODULUS EXPERIMENT")
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

    primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # --------------------------------------------------------------------------------------------
    # ANCHORS
    # --------------------------------------------------------------------------------------------

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        rng,
    )

    # --------------------------------------------------------------------------------------------
    # MODULUS TRIPLES
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    triples = []

    for i, (p, q) in enumerate(anchors, 1):

        n = p * q

        mods = choose_close_triple(
            n,
            modulus_primes,
        )

        triples.append((n, p, q, mods))

        if i % PROGRESS == 0:
            print(f"anchor {i:3d}/{len(anchors)}")

    print()
    print(f"usable anchors            = {len(triples)}")
    print()

    # --------------------------------------------------------------------------------------------
    # TRUE-STATE ALGEBRA CHECK
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("VALIDATING EXACT QUOTIENT IDENTITY")
    print("=" * 100)

    identity_failures = 0

    for i, (n, p, q, mods) in enumerate(triples, 1):

        try:
            r1, r2, r3 = mods

            validate_true_state(
                n,
                p,
                q,
                r1,
                r2,
                r3,
            )

        except AssertionError:
            identity_failures += 1

        if i % PROGRESS == 0:
            print(f"anchor {i:3d}/{len(triples)}")

    print()
    print(
        f"identity failures         = {identity_failures}"
    )
    print()

    if identity_failures:
        raise RuntimeError(
            "Exact quotient identity failed. "
            "The implementation must not continue."
        )

    # --------------------------------------------------------------------------------------------
    # FULL EXPERIMENT
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("RUNNING LINEARIZED r2 SEARCH")
    print("=" * 100)

    t0 = time.perf_counter()

    results = []

    for i, (n, p, q, mods) in enumerate(
        triples[:FULL_SEARCH_ANCHORS],
        1,
    ):

        result = solve_anchor(
            n,
            p,
            q,
            mods,
        )

        results.append(
            (n, p, q, result)
        )

        if i % PROGRESS == 0:
            print(
                f"full-search anchor "
                f"{i:3d}/{FULL_SEARCH_ANCHORS}"
            )

    elapsed = time.perf_counter() - t0

    # --------------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    if not results:
        print("No full-search anchors.")
        return

    residue_states = [
        x[3]["residue_states"]
        for x in results
    ]

    k_states = [
        x[3]["k_states"]
        for x in results
    ]

    candidate_counts = [
        x[3]["candidate_count"]
        for x in results
    ]

    interval_counts = [
        x[3]["interval_count"]
        for x in results
    ]

    r3_counts = [
        x[3]["r3_count"]
        for x in results
    ]

    prime_counts = [
        x[3]["prime_count"]
        for x in results
    ]

    exact_counts = [
        x[3]["exact_count"]
        for x in results
    ]

    print(
        f"full-search anchors       = {len(results)}"
    )

    print(
        f"average r1 residue states = "
        f"{statistics.mean(residue_states):.3f}"
    )

    print(
        f"average k states          = "
        f"{statistics.mean(k_states):.3f}"
    )

    print(
        f"average (a,k) states      = "
        f"{statistics.mean(residue_states) * statistics.mean(k_states):.3f}"
    )

    print(
        f"average solved l states   = "
        f"{statistics.mean(candidate_counts):.3f}"
    )

    print(
        f"average interval states   = "
        f"{statistics.mean(interval_counts):.3f}"
    )

    print(
        f"average r3 survivors      = "
        f"{statistics.mean(r3_counts):.3f}"
    )

    print(
        f"average prime candidates  = "
        f"{statistics.mean(prime_counts):.3f}"
    )

    print(
        f"average exact solutions   = "
        f"{statistics.mean(exact_counts):.3f}"
    )

    print()
    print("=" * 100)
    print("BASELINE COMPARISON")
    print("=" * 100)

    baseline = len(factor_primes)

    mean_candidates = statistics.mean(candidate_counts)
    mean_r3 = statistics.mean(r3_counts)

    print(
        f"ordinary factor-prime enumeration = {baseline:,}"
    )

    print(
        f"average solved-l candidates        = "
        f"{mean_candidates:.3f}"
    )

    print(
        f"average third-modulus survivors     = "
        f"{mean_r3:.3f}"
    )

    print(
        f"candidate ratio                     = "
        f"{mean_candidates / baseline:.9f}"
    )

    print(
        f"r3 ratio                            = "
        f"{mean_r3 / baseline:.9f}"
    )

    if mean_candidates:
        print(
            f"candidate reduction                = "
            f"{100.0 * (1.0 - mean_candidates / baseline):.6f}%"
        )

    if mean_r3:
        print(
            f"r3 reduction                        = "
            f"{100.0 * (1.0 - mean_r3 / baseline):.6f}%"
        )

    # --------------------------------------------------------------------------------------------
    # RECOVERY
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    recovered = sum(
        1 for _, _, _, r in results
        if r["actual_found"]
    )

    print(
        f"correctly recovered       = "
        f"{recovered}/{len(results)}"
    )

    # --------------------------------------------------------------------------------------------
    # DISTRIBUTIONS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SOLVED-l CANDIDATE DISTRIBUTION")
    print("=" * 100)

    for value, count in sorted(
        Counter(candidate_counts).items()
    ):
        print(
            f"candidates = {value:6d} "
            f"anchors = {count:3d}"
        )

    print()
    print("=" * 100)
    print("R3 SURVIVOR DISTRIBUTION")
    print("=" * 100)

    for value, count in sorted(
        Counter(r3_counts).items()
    ):
        print(
            f"r3 survivors = {value:6d} "
            f"anchors = {count:3d}"
        )

    # --------------------------------------------------------------------------------------------
    # EXAMPLES
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for n, p, q, r in results[:20]:

        mods = (
            r["r1"],
            r["r2"],
            r["r3"],
        )

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods={mods} "
            f"candidates={r['candidate_count']} "
            f"interval={r['interval_count']} "
            f"r3={r['r3_count']} "
            f"prime={r['prime_count']} "
            f"exact={r['exact_count']} "
            f"recovered={r['actual_found']}"
        )

    # --------------------------------------------------------------------------------------------
    # MOST AGGRESSIVE
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOST AGGRESSIVE COLLAPSES")
    print("=" * 100)

    ranked = sorted(
        results,
        key=lambda x: (
            x[3]["r3_count"],
            x[3]["candidate_count"],
        ),
    )

    for n, p, q, r in ranked[:20]:

        R = r["r1"] * r["r2"] * r["r3"]

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods=({r['r1']},{r['r2']},{r['r3']}) "
            f"R/n={R/n:.10f} "
            f"candidate={r['candidate_count']} "
            f"r3={r['r3_count']} "
            f"exact={r['exact_count']}"
        )

    # --------------------------------------------------------------------------------------------
    # ANCHOR TABLE
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q      r1     r2     r3"
        "     RESIDUE    K       L      R3    EXACT"
    )

    print("-" * 100)

    for i, (n, p, q, r) in enumerate(results, 1):

        print(
            f"{i:3d} "
            f"{p:8,d} "
            f"{q:8,d} "
            f"{r['r1']:6d} "
            f"{r['r2']:6d} "
            f"{r['r3']:6d} "
            f"{r['residue_states']:9,d} "
            f"{r['k_states']:6,d} "
            f"{r['candidate_count']:7,d} "
            f"{r['r3_count']:5,d} "
            f"{r['exact_count']:6,d}"
        )

    # --------------------------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"full-search runtime      = {elapsed:.3f} seconds"
    )

    # --------------------------------------------------------------------------------------------
    # MATHEMATICAL INTERPRETATION
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
The previous experiment failed because it silently replaced:

    a*b == g (mod r1)

with an equation modulo r2.

That is invalid.

The corrected construction retains the exact quotient:

    h = (a*b-g)/r1.

For

    p = a+k*r1
    q = b+l*r1
    r2 = r1+Delta,

the exact second-modulus condition becomes:

    h*r1
      - a*l*Delta
      - b*k*Delta
      + k*l*Delta^2
      == 0 (mod r2).

For fixed a,b,k this is linear in l:

    l*(k*Delta^2-a*Delta)
      ==
    b*k*Delta-h*r1
      (mod r2).

Therefore:

    (a,b,k)
          |
          v
    solve l modulo r2
          |
          v
    only a tiny arithmetic progression
    of l values can fit the factor interval
          |
          v
    verify r3
          |
          v
    exact p*q=n.

This is the actual close-modulus coordinate elimination.

The critical number is NOT merely the percentage reduction.

The useful quantity is:

    number of solved l candidates

per (a,k) state.

If that is consistently O(1), and the r3 stage also becomes O(1),
the next question is whether we can eliminate the remaining a/k
enumeration itself.

The script deliberately reports all stages separately so that a
coordinate transformation cannot be mistaken for a factoring speedup.

The identity is independently checked against every true factor pair
before the search begins.
"""
    )

    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
