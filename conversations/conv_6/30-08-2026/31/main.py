#!/usr/bin/env python3

from __future__ import annotations

import math
import random
import time
from collections import Counter
from statistics import mean


# ================================================================================================
# CONFIGURATION
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

# Run 20 first. Set to ANCHORS after validation.
SEARCH_ANCHORS = 20


# ================================================================================================
# OUTPUT
# ================================================================================================

def banner(s: str) -> None:
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)


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

    return [i for i, v in enumerate(a) if v]


def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


# ================================================================================================
# CLOSE MODULUS TRIPLE
# ================================================================================================

def choose_close_triple(
    n: int,
    primes: list[int],
) -> tuple[int, int, int]:

    best = None
    best_R = -1

    for i, r1 in enumerate(primes):

        limit = int(r1 * (1.0 + CLOSE_RATIO))

        for j in range(i + 1, len(primes)):

            r2 = primes[j]

            if r2 > limit:
                break

            for k in range(j + 1, len(primes)):

                r3 = primes[k]

                if r3 > limit:
                    break

                R = r1 * r2 * r3

                if R < n and R > best_R:
                    best_R = R
                    best = (r1, r2, r3)

    if best is None:
        raise RuntimeError(
            f"No usable triple for n={n}"
        )

    return best


# ================================================================================================
# ANCHORS
# ================================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    rng: random.Random,
) -> list[tuple[int, int]]:

    out = []

    for _ in range(count):
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        out.append((p, q))

    return out


# ================================================================================================
# EXACT IDENTITY VALIDATION
# ================================================================================================

def validate_anchor(
    n: int,
    p: int,
    q: int,
    mods: tuple[int, int, int],
) -> None:

    r1, r2, r3 = mods

    g = n % r1
    t = (n - g) // r1

    a = p % r1
    b = q % r1

    k = (p - a) // r1
    ell = (q - b) // r1

    h_num = a * b - g

    if h_num % r1 != 0:
        raise AssertionError(
            "r1 quotient failure"
        )

    h = h_num // r1

    assert n == g + t * r1
    assert a * b == g + h * r1

    d2 = r2 - r1
    d3 = r3 - r1

    # Exact r2 identity.
    e2 = (
        (h - t) * r1
        - a * ell * d2
        - b * k * d2
        + k * ell * d2 * d2
    )

    if e2 % r2 != 0:
        raise AssertionError(
            f"r2 identity failed: {e2}"
        )

    # Exact r3 identity.
    e3 = (
        (h - t) * r1
        - a * ell * d3
        - b * k * d3
        + k * ell * d3 * d3
    )

    if e3 % r3 != 0:
        raise AssertionError(
            f"r3 identity failed: {e3}"
        )


# ================================================================================================
# QUADRATIC UTILITIES
# ================================================================================================

def solve_linear_mod(
    A: int,
    B: int,
    mod: int,
) -> list[int]:

    A %= mod
    B %= mod

    if A == 0:
        return list(range(mod)) if B == 0 else []

    return [(B * pow(A, -1, mod)) % mod]


def solve_quadratic_prime(
    A: int,
    B: int,
    C: int,
    p: int,
) -> list[int]:

    A %= p
    B %= p
    C %= p

    if A == 0:
        return solve_linear_mod(B, (-C) % p, p)

    # Discriminant.
    D = (B * B - 4 * A * C) % p

    if D == 0:
        x = (-B * pow(2 * A, -1, p)) % p
        return [x]

    # Euler criterion.
    if pow(D, (p - 1) // 2, p) != 1:
        return []

    # Python pow(..., (p+1)//4, p) works for p == 3 mod 4.
    if p % 4 == 3:
        sqrt_D = pow(D, (p + 1) // 4, p)
        roots = {
            sqrt_D,
            (-sqrt_D) % p,
        }
    else:
        # Small prime range: use a direct square-root table.
        roots = {
            x
            for x in range(p)
            if x * x % p == D
        }

    inv_2A = pow(2 * A, -1, p)

    out = set()

    for s in roots:
        x = ((-B + s) * inv_2A) % p
        out.add(x)

    return sorted(out)


# ================================================================================================
# DERIVE THE r2 LINEAR EQUATION
# ================================================================================================

def r2_linear_coefficients(
    a: int,
    b: int,
    h: int,
    t: int,
    k: int,
    r1: int,
    r2: int,
) -> tuple[int, int]:

    delta = r2 - r1

    # l*(k*d²-a*d)
    #       =
    # b*k*d-(h-t)*r1
    #
    A = (
        k * delta * delta
        - a * delta
    )

    B = (
        b * k * delta
        - (h - t) * r1
    )

    return A % r2, B % r2


# ================================================================================================
# SEARCH
# ================================================================================================

def solve_anchor(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
) -> dict:

    r1, r2, r3 = mods

    g = n % r1
    t = (n - g) // r1

    d2 = r2 - r1
    d3 = r3 - r1

    # p = a + k*r1
    k_min = max(
        0,
        (FACTOR_MIN - (r1 - 1) + r1 - 1) // r1,
    )

    k_max = FACTOR_MAX // r1

    residue_states = 0
    k_states = 0

    r2_linear_solutions = 0
    r3_quadratic_roots = 0

    interval_candidates = 0
    prime_candidates = 0
    exact_pairs: set[tuple[int, int]] = set()

    actual_found = False

    # How many quadratic k roots are produced for each a.
    quadratic_distribution = Counter()

    # --------------------------------------------------------------------------------------------
    # Enumerate only a.
    #
    # b is fixed by:
    #
    #     a*b = g (mod r1)
    #
    # --------------------------------------------------------------------------------------------

    for a in range(1, r1):

        b = (
            g * pow(a, -1, r1)
        ) % r1

        h_num = a * b - g

        if h_num % r1 != 0:
            raise AssertionError(
                "Bad r1 state"
            )

        h = h_num // r1

        residue_states += 1

        # ----------------------------------------------------------------------------------------
        # Instead of scanning k and then using r3 as a filter, derive the r3 equation after
        # substituting the r2 solution.
        #
        # The r2 equation is:
        #
        #     l * A2(k) = B2(k)   mod r2
        #
        # Since l is the bounded quotient coordinate, when r2 > l-range it is normally the
        # unique representative l0. The resulting r3 condition is then tested as a function
        # of k.
        #
        # We construct a polynomial by evaluating the resulting expression at several k values
        # and interpolating modulo r3.
        #
        # This is intentionally an experimental "quotient resultant" construction.
        # ----------------------------------------------------------------------------------------

        k_samples = [
            k_min,
            min(k_min + 1, k_max),
            min(k_min + 2, k_max),
            min(k_min + 3, k_max),
        ]

        sample_values: list[tuple[int, int]] = []

        for k in k_samples:

            p = a + k * r1

            if not (
                FACTOR_MIN <= p <= FACTOR_MAX
            ):
                continue

            A2, B2 = r2_linear_coefficients(
                a,
                b,
                h,
                t,
                k,
                r1,
                r2,
            )

            if A2 == 0:
                continue

            l = (
                B2 *
                pow(A2, -1, r2)
            ) % r2

            # Only use the canonical r2 solution when it lies in the actual quotient range.
            l_low = max(
                0,
                (FACTOR_MIN - b + r1 - 1) // r1,
            )

            l_high = (
                FACTOR_MAX - b
            ) // r1

            if l_low <= l <= l_high:
                x3 = (
                    a - k * d3
                ) % r3

                y3 = (
                    b - l * d3
                ) % r3

                value = (
                    x3 * y3 - n
                ) % r3

                sample_values.append(
                    (k % r3, value)
                )

        # ----------------------------------------------------------------------------------------
        # Because the canonical l(k) includes an inverse modulo r2, it is not globally a
        # polynomial over r3. Therefore we do NOT pretend the interpolation is exact.
        #
        # Instead, use the experimentally constructed quadratic resultant only as a candidate
        # generator and always verify every candidate against the original congruences.
        #
        # The fallback candidate domain is k modulo r3.
        # ----------------------------------------------------------------------------------------

        candidates_k: set[int] = set()

        # Direct quadratic approximation from freezing the r2 inverse relation.
        #
        # This is deliberately treated as a hypothesis generator, never as proof.
        #
        # Use all k residues modulo r3 and obtain candidate roots from the actual r3 product
        # relation with the corresponding r2-derived l.
        #
        # This gives a measurable "coordinate elimination" baseline.

        for k_residue in range(r3):

            k = k_residue

            if k < k_min:
                k += (
                    (k_min - k + r3 - 1)
                    // r3
                ) * r3

            if k > k_max:
                continue

            p = a + k * r1

            if not (
                FACTOR_MIN <= p <= FACTOR_MAX
            ):
                continue

            k_states += 1

            A2, B2 = r2_linear_coefficients(
                a,
                b,
                h,
                t,
                k,
                r1,
                r2,
            )

            if A2 == 0:
                continue

            l = (
                B2 *
                pow(A2, -1, r2)
            ) % r2

            l_low = max(
                0,
                (FACTOR_MIN - b + r1 - 1) // r1,
            )

            l_high = (
                FACTOR_MAX - b
            ) // r1

            if not (
                l_low <= l <= l_high
            ):
                continue

            r2_linear_solutions += 1

            x3 = (
                a - k * d3
            ) % r3

            y3 = (
                b - l * d3
            ) % r3

            if (
                x3 * y3 - n
            ) % r3 != 0:
                continue

            r3_quadratic_roots += 1

            candidates_k.add(k)

            q = b + l * r1

            interval_candidates += 1

            if p * q != n:
                continue

            pair = tuple(
                sorted((p, q))
            )

            exact_pairs.add(pair)

            if (
                is_prime(p)
                and is_prime(q)
            ):
                prime_candidates += 1

                if pair == tuple(
                    sorted((p_true, q_true))
                ):
                    actual_found = True

        quadratic_distribution[len(candidates_k)] += 1

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "R": r1 * r2 * r3,
        "R_over_n": (r1 * r2 * r3) / n,
        "g": g,
        "t": t,
        "d2": d2,
        "d3": d3,
        "residue_states": residue_states,
        "k_states": k_states,
        "r2_linear_solutions": r2_linear_solutions,
        "r3_quadratic_roots": r3_quadratic_roots,
        "interval_candidates": interval_candidates,
        "prime_candidates": prime_candidates,
        "exact_pairs": sorted(exact_pairs),
        "actual_found": actual_found,
        "quadratic_distribution": quadratic_distribution,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    rng = random.Random(SEED)

    banner(
        "THREE-CLOSE-PRIME r3 COORDINATE-ELIMINATION / "
        "QUOTIENT-RESULTANT EXPERIMENT"
    )

    print(
        f"M                         = {M:,}"
    )
    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(
        f"anchors                   = {ANCHORS}"
    )
    print(
        f"modulus prime range       = "
        f"{MOD_MIN:,} - {MOD_MAX:,}"
    )
    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.0%}"
    )
    print(
        f"seed                      = {SEED:,}"
    )

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    banner("BUILDING PRIME POOLS")

    all_primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in all_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    # ============================================================================================
    # ANCHORS
    # ============================================================================================

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        rng,
    )

    banner(
        "SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES"
    )

    triples = []

    for i, (p, q) in enumerate(
        anchors,
        1,
    ):

        n = p * q

        mods = choose_close_triple(
            n,
            modulus_primes,
        )

        triples.append(
            (n, p, q, mods)
        )

        if i % PROGRESS == 0:
            print(
                f"anchor {i:3d}/{ANCHORS}"
            )

    print()
    print(
        f"usable anchors            = "
        f"{len(triples)}"
    )

    # ============================================================================================
    # EXACT ALGEBRA VALIDATION
    # ============================================================================================

    banner(
        "VALIDATING THREE-MODULUS ALGEBRA"
    )

    failures = 0

    for i, (
        n,
        p,
        q,
        mods,
    ) in enumerate(
        triples,
        1,
    ):

        try:

            validate_anchor(
                n,
                p,
                q,
                mods,
            )

        except AssertionError as exc:

            failures += 1

            print()
            print(
                f"FIRST FAILURE: anchor {i}"
            )
            print(
                f"n={n:,} p={p:,} q={q:,}"
            )
            print(
                f"mods={mods}"
            )
            print(
                f"error={exc}"
            )

            break

        if i % PROGRESS == 0:
            print(
                f"anchor {i:3d}/{ANCHORS}"
            )

    print()
    print(
        f"identity failures         = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Three-modulus algebra validation failed."
        )

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    search_data = triples[
        :SEARCH_ANCHORS
    ]

    banner(
        "RUNNING r3 COORDINATE-ELIMINATION SEARCH"
    )

    print(
        f"search anchors             = "
        f"{len(search_data)}"
    )

    start = time.perf_counter()

    results = []

    for i, (
        n,
        p,
        q,
        mods,
    ) in enumerate(
        search_data,
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
                f"anchor {i:3d}/{len(search_data)}"
            )

    elapsed = time.perf_counter() - start

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    banner("SUMMARY")

    print(
        f"anchors analyzed              = "
        f"{len(results)}"
    )

    avg_a = mean(
        r["residue_states"]
        for _, _, _, r in results
    )

    avg_k = mean(
        r["k_states"]
        for _, _, _, r in results
    )

    avg_r2 = mean(
        r["r2_linear_solutions"]
        for _, _, _, r in results
    )

    avg_r3 = mean(
        r["r3_quadratic_roots"]
        for _, _, _, r in results
    )

    avg_exact = mean(
        r["interval_candidates"]
        for _, _, _, r in results
    )

    avg_prime = mean(
        r["prime_candidates"]
        for _, _, _, r in results
    )

    print(
        f"average a states              = "
        f"{avg_a:.3f}"
    )

    print(
        f"average k residues tested     = "
        f"{avg_k:.3f}"
    )

    print(
        f"average r2-compatible states  = "
        f"{avg_r2:.3f}"
    )

    print(
        f"average r3 survivors          = "
        f"{avg_r3:.3f}"
    )

    print(
        f"average final candidates      = "
        f"{avg_exact:.3f}"
    )

    print(
        f"average prime candidates      = "
        f"{avg_prime:.3f}"
    )

    # ============================================================================================
    # BASELINE
    # ============================================================================================

    banner("SEARCH REDUCTION")

    baseline = len(factor_primes)

    print(
        f"ordinary factor-prime tests   = "
        f"{baseline:,}"
    )

    print(
        f"r2-compatible / baseline      = "
        f"{avg_r2 / baseline:.9f}"
    )

    print(
        f"r3 survivors / baseline       = "
        f"{avg_r3 / baseline:.9f}"
    )

    print(
        f"prime candidates / baseline   = "
        f"{avg_prime / baseline:.9f}"
    )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    banner("RECOVERY")

    recovered = sum(
        1
        for _, _, _, r in results
        if r["actual_found"]
    )

    print(
        f"correctly recovered            = "
        f"{recovered}/{len(results)}"
    )

    # ============================================================================================
    # R3 DISTRIBUTION
    # ============================================================================================

    banner(
        "r3 SURVIVOR DISTRIBUTION"
    )

    dist = Counter(
        r["r3_quadratic_roots"]
        for _, _, _, r in results
    )

    for value, count in sorted(
        dist.items()
    ):
        print(
            f"r3 survivors = "
            f"{value:6d} "
            f"anchors = "
            f"{count:3d}"
        )

    # ============================================================================================
    # STRONGEST COLLAPSES
    # ============================================================================================

    banner(
        "STRONGEST COORDINATE COLLAPSES"
    )

    ranked = sorted(
        results,
        key=lambda x: (
            x[3]["r3_quadratic_roots"],
            x[3]["r2_linear_solutions"],
        ),
    )

    for n, p, q, r in ranked[:20]:

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods=({r['r1']},{r['r2']},{r['r3']}) "
            f"R/n={r['R_over_n']:.10f} "
            f"r2={r['r2_linear_solutions']} "
            f"r3={r['r3_quadratic_roots']} "
            f"prime={r['prime_candidates']} "
            f"exact={len(r['exact_pairs'])}"
        )

    # ============================================================================================
    # EXAMPLES
    # ============================================================================================

    banner("EXAMPLES")

    for n, p, q, r in results[:20]:

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods=({r['r1']},{r['r2']},{r['r3']}) "
            f"g={r['g']} "
            f"t={r['t']} "
            f"a={r['residue_states']} "
            f"r2={r['r2_linear_solutions']} "
            f"r3={r['r3_quadratic_roots']} "
            f"recovered={r['actual_found']}"
        )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    banner("TIMING")

    print(
        f"search runtime               = "
        f"{elapsed:.3f} seconds"
    )

    # ============================================================================================
    # INTERPRETATION
    # ============================================================================================

    banner("MATHEMATICAL INTERPRETATION")

    print(
        r"""
The corrected experiment established the exact identity:

    n = g + t*r1

and:

    a*b = g + h*r1.

Therefore:

    a*b - n = (h-t)*r1.

With:

    p = a + k*r1
    q = b + l*r1
    Delta_i = r_i-r1,

the second modulus gives:

    l*(k*Delta_2^2-a*Delta_2)
       =
    b*k*Delta_2-(h-t)*r1
       (mod r2).

The previous experiment used this relation to solve l for each
(a,k).

This experiment asks a different question:

    Can the third-modulus relation eliminate k itself?

The third equation is:

    l*(k*Delta_3^2-a*Delta_3)
      =
    b*k*Delta_3-(h-t)*r1
      (mod r3).

Thus r2 and r3 provide two different linear constraints involving
the same quotient coordinates k and l.

The hoped-for structure is:

    r2 equation
         |
         v
    eliminate l
         |
         v
    resultant in k
         |
         v
    solve candidate k
         |
         v
    reconstruct l
         |
         v
    exact product

The critical statistic is therefore not merely the number of r3
survivors. It is whether the third modulus can replace an enumeration
over the approximately 90,000 possible quotient positions.

There is an important caveat:

    solving l modulo r2 does not automatically produce a polynomial
    function of k modulo r3 because the inverse coefficient is taken
    modulo r2.

Therefore every generated k candidate is independently verified
against the original r2 and r3 equations and finally against:

    p*q == n.

A successful result would require both:

    very few k candidates

and:

    100% exact recovery.

If r3 still behaves merely as a filter over a large k domain, then
the close-modulus equations are primarily performing coordinate
compression rather than eliminating the factor search.
"""
    )

    banner("EXPERIMENT COMPLETE")


if __name__ == "__main__":
    run()
