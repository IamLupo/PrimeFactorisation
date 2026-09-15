#!/usr/bin/env python3

"""
============================================================================================
E(n+x) LOCAL STABILITY / FACTOR-TRAJECTORY EXPERIMENT
============================================================================================

Question:

    Does

        E(n+1) = E(n) = E(n-1)

    actually hold when n+1 and n-1 are independently factorized?

We distinguish:

    A) FIXED-(k,l) E:
       Keep the original k,l and vary n -> n+x.
       This tests the exact quotient/remainder staircase.

    B) ACTUAL FACTORIZATION E:
       Factor n+x independently.
       Compute its own residues, k_x, l_x and E_x.

The interesting question is whether the equality survives
the change in factorization.

All results are exact.

============================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter

try:
    import sympy as sp
except ImportError:
    raise SystemExit(
        "sympy is required.\n"
        "Install with:\n"
        "    pip install sympy"
    )


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

X_RADIUS = 10

SEED = 1_511_464_998

# Only count x values whose number has a two-factor decomposition
# inside the same factor range.
REQUIRE_FACTOR_RANGE = True


# ==========================================================================================
# PRIME POOLS
# ==========================================================================================

def build_primes(lo: int, hi: int) -> list[int]:
    return list(sp.primerange(lo, hi + 1))


# ==========================================================================================
# ANCHORS
# ==========================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[tuple[int, int, int]]:
    rng = random.Random(seed)

    anchors = []
    seen = set()

    while len(anchors) < count:
        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        seen.add(n)
        anchors.append((p, q, n))

    return anchors


# ==========================================================================================
# CLOSE MODULUS PAIR
# ==========================================================================================

def choose_close_pair(
    modulus_primes: list[int],
    rng: random.Random,
) -> tuple[int, int]:

    candidates = []

    for r1 in modulus_primes:
        for r2 in modulus_primes:
            if r1 >= r2:
                continue

            ratio = (r2 - r1) / r1

            if ratio <= CLOSE_RATIO:
                candidates.append((r1, r2))

    if not candidates:
        raise RuntimeError("No close modulus pairs found.")

    return rng.choice(candidates)


# ==========================================================================================
# COORDINATES
# ==========================================================================================

def coordinates(
    p: int,
    q: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int]:
    a = p % r1
    b = q % r2

    k = (p - a) // r1
    l = (q - b) // r2

    return a, b, k, l


# ==========================================================================================
# E COMPUTATION
# ==========================================================================================

def compute_E(
    n: int,
    k: int,
    l: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int]:
    R = r1 * r2

    T = n // R
    K = k * l
    E = T - K
    g = n % R

    return T, K, E, g


# ==========================================================================================
# FIXED-(k,l) MODEL
# ==========================================================================================

def fixed_coordinate_E(
    n: int,
    x: int,
    k: int,
    l: int,
    R: int,
) -> tuple[int, int, int]:
    """
    Keep k,l fixed and calculate

        E_x = floor((n+x)/R) - k*l

    together with the quotient and remainder.
    """
    m = n + x
    T = m // R
    E = T - k * l
    g = m % R

    return T, E, g


# ==========================================================================================
# FACTOR n+x
# ==========================================================================================

def factor_nearby(
    m: int,
    factor_min: int,
    factor_max: int,
) -> tuple[int, int] | None:

    if m <= 0:
        return None

    # Exact factorization.
    fac = sp.factorint(m)

    # We specifically want two prime factors.
    factors = []

    for p, exponent in fac.items():
        factors.extend([int(p)] * exponent)

    if len(factors) != 2:
        return None

    p, q = sorted(factors)

    if REQUIRE_FACTOR_RANGE:
        if not (factor_min <= p <= factor_max):
            return None
        if not (factor_min <= q <= factor_max):
            return None

    return p, q


# ==========================================================================================
# ACTUAL E FOR n+x
# ==========================================================================================

def actual_E_for_factorization(
    m: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
) -> dict:

    a, b, k, l = coordinates(p, q, r1, r2)
    T, K, E, g = compute_E(m, k, l, r1, r2)

    return {
        "m": m,
        "p": p,
        "q": q,
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "E": E,
        "g": g,
    }


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    t0 = time.perf_counter()

    print("=" * 100)
    print("E(n+x) LOCAL STABILITY / FACTOR-TRAJECTORY EXPERIMENT")
    print("=" * 100)
    print(f"N anchors                 = {N_ANCHORS}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"x radius                  = {X_RADIUS}")
    print(f"seed                      = {SEED}")
    print()

    # --------------------------------------------------------------------------------------
    # PRIME POOLS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = build_primes(FACTOR_MIN, FACTOR_MAX)
    modulus_primes = build_primes(MOD_MIN, MOD_MAX)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # --------------------------------------------------------------------------------------
    # ANCHORS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")
    print()

    # --------------------------------------------------------------------------------------
    # RESULT STORAGE
    # --------------------------------------------------------------------------------------

    rng = random.Random(SEED + 1234567)

    fixed_equal_n_pm1 = 0
    actual_equal_n_p1 = 0
    actual_equal_n_m1 = 0
    actual_equal_triplet = 0

    actual_available_p1 = 0
    actual_available_m1 = 0

    actual_x_values = 0
    factorization_failures = 0

    fixed_delta_counter = Counter()

    actual_delta_plus = []
    actual_delta_minus = []

    examples = []

    # --------------------------------------------------------------------------------------
    # MAIN LOOP
    # --------------------------------------------------------------------------------------

    for idx, (p0, q0, n) in enumerate(anchors, 1):

        if idx % 10 == 0 or idx == 1:
            print(f"anchor {idx:3d}/{len(anchors)}")

        r1, r2 = choose_close_pair(modulus_primes, rng)

        R = r1 * r2

        # Original coordinates
        a0, b0, k0, l0 = coordinates(
            p0,
            q0,
            r1,
            r2,
        )

        T0, K0, E0, g0 = compute_E(
            n,
            k0,
            l0,
            r1,
            r2,
        )

        # ----------------------------------------------------------------------
        # FIXED COORDINATE TEST
        # ----------------------------------------------------------------------

        fixed = {}

        for x in range(-X_RADIUS, X_RADIUS + 1):

            T_x, E_x, g_x = fixed_coordinate_E(
                n,
                x,
                k0,
                l0,
                R,
            )

            fixed[x] = {
                "T": T_x,
                "E": E_x,
                "g": g_x,
            }

            fixed_delta_counter[E_x - E0] += 1

        # The x = +1 / -1 fixed-coordinate test.
        fixed_pm1_ok = (
            fixed[1]["E"] == E0
            and fixed[-1]["E"] == E0
        )

        if fixed_pm1_ok:
            fixed_equal_n_pm1 += 1

        # ----------------------------------------------------------------------
        # ACTUAL FACTORIZATION TEST
        # ----------------------------------------------------------------------

        actual = {}

        for x in range(-X_RADIUS, X_RADIUS + 1):

            if x == 0:
                actual[x] = {
                    "m": n,
                    "p": p0,
                    "q": q0,
                    "E": E0,
                    "K": K0,
                    "k": k0,
                    "l": l0,
                    "g": g0,
                }
                actual_x_values += 1
                continue

            m = n + x

            pair = factor_nearby(
                m,
                FACTOR_MIN,
                FACTOR_MAX,
            )

            if pair is None:
                factorization_failures += 1
                continue

            px, qx = pair

            data = actual_E_for_factorization(
                m,
                px,
                qx,
                r1,
                r2,
            )

            actual[x] = data
            actual_x_values += 1

        # ----------------------------------------------------------------------
        # ACTUAL +/-1
        # ----------------------------------------------------------------------

        plus_ok = False
        minus_ok = False

        if 1 in actual:
            actual_available_p1 += 1
            plus_ok = actual[1]["E"] == E0

            actual_delta_plus.append(
                actual[1]["E"] - E0
            )

            if plus_ok:
                actual_equal_n_p1 += 1

        if -1 in actual:
            actual_available_m1 += 1
            minus_ok = actual[-1]["E"] == E0

            actual_delta_minus.append(
                actual[-1]["E"] - E0
            )

            if minus_ok:
                actual_equal_n_m1 += 1

        # Both actual factorizations available and both equal.
        if plus_ok and minus_ok:
            actual_equal_triplet += 1

        # ----------------------------------------------------------------------
        # EXAMPLES
        # ----------------------------------------------------------------------

        if len(examples) < 20:
            examples.append(
                {
                    "n": n,
                    "p": p0,
                    "q": q0,
                    "r1": r1,
                    "r2": r2,
                    "R": R,
                    "k": k0,
                    "l": l0,
                    "E0": E0,
                    "g0": g0,
                    "fixed_m1": fixed[-1]["E"],
                    "fixed_p1": fixed[1]["E"],
                    "actual_m1": actual.get(-1),
                    "actual_p1": actual.get(1),
                }
            )

    # ======================================================================================
    # SUMMARY
    # ======================================================================================

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(f"anchors analyzed             = {len(anchors):,}")
    print(f"nearby factorizations tested = {actual_x_values:,}")
    print(f"factorization unavailable    = {factorization_failures:,}")
    print()

    print("FIXED-(k,l) TEST")
    print(
        "    E(n-1)=E(n)=E(n+1)      = "
        f"{fixed_equal_n_pm1}/{len(anchors)}"
    )
    print()

    print("ACTUAL FACTORIZATION TEST")
    print(
        "    n+1 factorizations       = "
        f"{actual_available_p1}/{len(anchors)}"
    )
    print(
        "    E(n+1)=E(n)              = "
        f"{actual_equal_n_p1}/{actual_available_p1}"
        if actual_available_p1
        else
        "    E(n+1)=E(n)              = N/A"
    )

    print(
        "    n-1 factorizations       = "
        f"{actual_available_m1}/{len(anchors)}"
    )
    print(
        "    E(n-1)=E(n)              = "
        f"{actual_equal_n_m1}/{actual_available_m1}"
        if actual_available_m1
        else
        "    E(n-1)=E(n)              = N/A"
    )

    both_available = sum(
        1
        for _ in range(len(anchors))
    )

    print(
        "    E(n-1)=E(n)=E(n+1)       = "
        f"{actual_equal_triplet}"
    )

    # ======================================================================================
    # DELTA DISTRIBUTION
    # ======================================================================================

    print()
    print("=" * 100)
    print("FIXED-COORDINATE E DELTA DISTRIBUTION")
    print("=" * 100)

    for delta, count in sorted(fixed_delta_counter.items()):
        print(
            f"E(n+x)-E(n) = {delta:3d} : {count:5d}"
        )

    # ======================================================================================
    # ACTUAL DELTAS
    # ======================================================================================

    print()
    print("=" * 100)
    print("ACTUAL FACTORIZATION DELTAS")
    print("=" * 100)

    if actual_delta_plus:
        print(
            f"mean  E(n+1)-E(n) = "
            f"{statistics.mean(actual_delta_plus):.6f}"
        )
        print(
            f"median E(n+1)-E(n) = "
            f"{statistics.median(actual_delta_plus):.6f}"
        )
        print(
            f"min   E(n+1)-E(n) = "
            f"{min(actual_delta_plus):d}"
        )
        print(
            f"max   E(n+1)-E(n) = "
            f"{max(actual_delta_plus):d}"
        )

    if actual_delta_minus:
        print(
            f"mean  E(n-1)-E(n) = "
            f"{statistics.mean(actual_delta_minus):.6f}"
        )
        print(
            f"median E(n-1)-E(n) = "
            f"{statistics.median(actual_delta_minus):.6f}"
        )
        print(
            f"min   E(n-1)-E(n) = "
            f"{min(actual_delta_minus):d}"
        )
        print(
            f"max   E(n-1)-E(n) = "
            f"{max(actual_delta_minus):d}"
        )

    # ======================================================================================
    # EXAMPLES
    # ======================================================================================

    print()
    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for ex in examples:

        print(
            f"n={ex['n']:,} "
            f"p={ex['p']:,} q={ex['q']:,} "
            f"mods=({ex['r1']},{ex['r2']}) "
            f"R={ex['R']:,}"
        )

        print(
            f"    original: k={ex['k']} l={ex['l']} "
            f"E={ex['E0']} g={ex['g0']}"
        )

        print(
            f"    fixed: "
            f"E(n-1)={ex['fixed_m1']} "
            f"E(n)={ex['E0']} "
            f"E(n+1)={ex['fixed_p1']}"
        )

        if ex["actual_m1"] is not None:
            d = ex["actual_m1"]
            print(
                f"    actual n-1: "
                f"p={d['p']:,} q={d['q']:,} "
                f"k={d['k']} l={d['l']} "
                f"E={d['E']} "
                f"delta={d['E']-ex['E0']:+d}"
            )
        else:
            print("    actual n-1: unavailable")

        if ex["actual_p1"] is not None:
            d = ex["actual_p1"]
            print(
                f"    actual n+1: "
                f"p={d['p']:,} q={d['q']:,} "
                f"k={d['k']} l={d['l']} "
                f"E={d['E']} "
                f"delta={d['E']-ex['E0']:+d}"
            )
        else:
            print("    actual n+1: unavailable")

        print()

    # ======================================================================================
    # MATHEMATICAL INTERPRETATION
    # ======================================================================================

    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)
    print(
        """
Write

    R = r1*r2
    n = T*R + g
    E = T-k*l.

For FIXED k,l:

    E(n+x)
      = floor((n+x)/R) - k*l
      = E(n) + floor((g+x)/R).

Therefore, whenever

    0 <= g+x < R,

we necessarily have

    E(n+x) = E(n).

In particular, if

    0 < g < R-1,

then

    E(n-1) = E(n) = E(n+1)

for the SAME k,l.

That statement is exact and contains no experimental uncertainty.

But this is not yet the interesting question.

For the ACTUAL FACTORIZATION of n+x we instead obtain

    n+x = p_x*q_x

with

    p_x = a_x + k_x*r1
    q_x = b_x + l_x*r2

and therefore

    E_x
      = floor((n+x)/R) - k_x*l_x.

Now k_x*l_x generally changes because the factorization of
n+x is different.

So:

    fixed-coordinate equality
        !=
    factor-trajectory equality.

The experiment therefore tests whether the latter happens
empirically.

If

    E(n+1) = E(n)

also survives under independent refactorization of n+1,
that would be much more interesting.

If instead the equality disappears immediately, then the
simple equality was only a consequence of holding k,l fixed.

The latter is actually what we should expect in general.

The strongest possible result would therefore be a nontrivial
trajectory such as

    E(n-1) ~= E(n) ~= E(n+1)

despite

    (k,l)_(n-1)
        !=
    (k,l)_n
        !=
    (k,l)_(n+1).

That would indicate that E is capturing some property of
nearby factorizations rather than merely the original quotient
coordinates.

The experiment also records the full fixed-coordinate staircase
over:

    x = -X_RADIUS ... +X_RADIUS.

This lets us compare the exact theoretical staircase against
the independently factorized trajectory.

All actual factor pairs are verified through exact
factorization, and only two-prime factorizations in the selected
range are considered.
"""
    )

    print("=" * 100)

    elapsed = time.perf_counter() - t0
    print("TIMING")
    print("=" * 100)
    print(
        f"total runtime               = {elapsed:.3f} seconds"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
