#!/usr/bin/env python3

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter


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

# Full expensive search can be limited while debugging.
# Set to ANCHORS after the algebra is confirmed.
FULL_SEARCH_ANCHORS = 20


# ================================================================================================
# OUTPUT
# ================================================================================================

def banner(title: str) -> None:
    print()
    print("=" * 100)
    print(title)
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

    return [i for i, x in enumerate(a) if x]


# ================================================================================================
# PRIME TEST
# ================================================================================================

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
# ANCHORS
# ================================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    rng: random.Random,
) -> list[tuple[int, int]]:

    anchors = []

    for _ in range(count):
        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        anchors.append((p, q))

    return anchors


# ================================================================================================
# CLOSE PRIME TRIPLE
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
        raise RuntimeError(
            f"No usable modulus triple for n={n}"
        )

    return best


# ================================================================================================
# EXACT ALGEBRA VALIDATION
# ================================================================================================

def validate_identity(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    r3: int,
) -> dict:

    # --------------------------------------------------------------------------------------------
    # First modulus.
    # --------------------------------------------------------------------------------------------

    g = n % r1

    a = p % r1
    b = q % r1

    k = (p - a) // r1
    ell = (q - b) // r1

    assert p == a + k * r1
    assert q == b + ell * r1

    # ab-g must be an exact multiple of r1.
    h_num = a * b - g

    if h_num % r1 != 0:
        raise AssertionError(
            f"r1 quotient failed: a={a}, b={b}, g={g}"
        )

    h = h_num // r1

    # --------------------------------------------------------------------------------------------
    # THIS WAS THE MISSING TERM IN THE PREVIOUS SCRIPT.
    #
    # n = g + t*r1
    # --------------------------------------------------------------------------------------------

    n_num = n - g

    if n_num % r1 != 0:
        raise AssertionError(
            "n-g is not divisible by r1"
        )

    t = n_num // r1

    assert n == g + t * r1

    # --------------------------------------------------------------------------------------------
    # SECOND MODULUS
    # --------------------------------------------------------------------------------------------

    delta2 = r2 - r1

    x2 = (a - k * delta2) % r2
    y2 = (b - ell * delta2) % r2

    # Direct truth check.
    if (x2 * y2 - n) % r2 != 0:
        raise AssertionError(
            f"direct r2 failure: "
            f"x2={x2}, y2={y2}, n={n}, r2={r2}"
        )

    # Correct expanded expression:
    #
    # h*r1
    # - t*r1
    # - a*ell*delta
    # - b*k*delta
    # + k*ell*delta^2
    #
    # must be 0 mod r2.

    exact_r2 = (
        (h - t) * r1
        - a * ell * delta2
        - b * k * delta2
        + k * ell * delta2 * delta2
    )

    if exact_r2 % r2 != 0:
        raise AssertionError(
            f"expanded r2 failure: "
            f"exact_r2={exact_r2}, "
            f"r2={r2}"
        )

    # --------------------------------------------------------------------------------------------
    # THIRD MODULUS
    # --------------------------------------------------------------------------------------------

    delta3 = r3 - r1

    x3 = (a - k * delta3) % r3
    y3 = (b - ell * delta3) % r3

    if (x3 * y3 - n) % r3 != 0:
        raise AssertionError(
            f"direct r3 failure: "
            f"x3={x3}, y3={y3}, n={n}, r3={r3}"
        )

    exact_r3 = (
        (h - t) * r1
        - a * ell * delta3
        - b * k * delta3
        + k * ell * delta3 * delta3
    )

    if exact_r3 % r3 != 0:
        raise AssertionError(
            f"expanded r3 failure: "
            f"exact_r3={exact_r3}, "
            f"r3={r3}"
        )

    return {
        "g": g,
        "a": a,
        "b": b,
        "k": k,
        "l": ell,
        "h": h,
        "t": t,
        "h_minus_t": h - t,
        "delta2": delta2,
        "delta3": delta3,
    }


# ================================================================================================
# SOLVE ONE ANCHOR
# ================================================================================================

def solve_anchor(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
) -> dict:

    r1, r2, r3 = mods

    delta2 = r2 - r1
    delta3 = r3 - r1

    # --------------------------------------------------------------------------------------------
    # n = g + t*r1
    # --------------------------------------------------------------------------------------------

    g = n % r1
    t = (n - g) // r1

    # --------------------------------------------------------------------------------------------
    # k domain.
    # --------------------------------------------------------------------------------------------

    k_min = max(
        0,
        (FACTOR_MIN - (r1 - 1) + r1 - 1) // r1,
    )

    k_max = FACTOR_MAX // r1

    residue_states = 0
    k_states = 0
    solved_l_states = 0
    interval_states = 0
    r3_states = 0
    prime_candidates = 0
    exact_solutions = 0

    exact_pairs = set()

    actual_found = False

    # Distribution of l residues produced modulo r2.
    l_residue_counter = Counter()

    # --------------------------------------------------------------------------------------------
    # r1 hyperbola:
    #
    #       a*b == g (mod r1)
    #
    # For prime r1 and g != 0, b is uniquely determined for every
    # nonzero a.
    #
    # --------------------------------------------------------------------------------------------

    for a in range(1, r1):

        b = (g * pow(a, -1, r1)) % r1

        h_num = a * b - g

        if h_num % r1 != 0:
            raise AssertionError(
                "Internal r1 quotient error"
            )

        h = h_num // r1

        residue_states += 1

        for k in range(k_min, k_max + 1):

            p = a + k * r1

            if not (FACTOR_MIN <= p <= FACTOR_MAX):
                continue

            k_states += 1

            # --------------------------------------------------------------------------------
            # Correct linearized equation.
            #
            # l*(k*delta²-a*delta)
            #     ==
            # b*k*delta - (h-t)*r1
            #                  (mod r2)
            # --------------------------------------------------------------------------------

            C = (
                k * delta2 * delta2
                - a * delta2
            ) % r2

            RHS = (
                b * k * delta2
                - (h - t) * r1
            ) % r2

            # Non-degenerate case.
            if C != 0:

                l0 = (
                    RHS *
                    pow(C, -1, r2)
                ) % r2

                l_residue_counter[l0] += 1

                # Allowed integer l range:
                #
                # q = b + l*r1
                #
                l_low = max(
                    0,
                    (FACTOR_MIN - b + r1 - 1) // r1,
                )

                l_high = (
                    FACTOR_MAX - b
                ) // r1

                # First l >= l_low satisfying:
                #
                # l == l0 (mod r2)
                #

                if l0 < l_low:
                    step = (
                        (l_low - l0 + r2 - 1)
                        // r2
                    )

                    ell = l0 + step * r2
                else:
                    ell = l0

                while ell <= l_high:

                    solved_l_states += 1

                    q = b + ell * r1

                    if not (
                        FACTOR_MIN <= q <= FACTOR_MAX
                    ):
                        ell += r2
                        continue

                    interval_states += 1

                    # ------------------------------------------------------------------------
                    # Independent third-modulus check.
                    # ------------------------------------------------------------------------

                    x3 = (
                        a - k * delta3
                    ) % r3

                    y3 = (
                        b - ell * delta3
                    ) % r3

                    if (x3 * y3 - n) % r3 != 0:
                        ell += r2
                        continue

                    r3_states += 1

                    # ------------------------------------------------------------------------
                    # Exact product.
                    # ------------------------------------------------------------------------

                    if p * q == n:

                        exact_solutions += 1

                        pair = tuple(sorted((p, q)))
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

                    ell += r2

                continue

            # --------------------------------------------------------------------------------
            # Degenerate case:
            #
            # C == 0.
            #
            # Then either no l works, or the congruence places no restriction
            # on l modulo r2.
            #
            # --------------------------------------------------------------------------------

            if RHS != 0:
                continue

            # Rare/degenerate branch.
            #
            # Enumerate only the actual bounded l interval.
            #
            l_low = max(
                0,
                (FACTOR_MIN - b + r1 - 1) // r1,
            )

            l_high = (
                FACTOR_MAX - b
            ) // r1

            for ell in range(l_low, l_high + 1):

                solved_l_states += 1

                q = b + ell * r1

                interval_states += 1

                x3 = (
                    a - k * delta3
                ) % r3

                y3 = (
                    b - ell * delta3
                ) % r3

                if (x3 * y3 - n) % r3 != 0:
                    continue

                r3_states += 1

                if p * q != n:
                    continue

                exact_solutions += 1

                pair = tuple(sorted((p, q)))
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

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "R": r1 * r2 * r3,
        "R_over_n": (r1 * r2 * r3) / n,
        "g": g,
        "t": t,
        "delta2": delta2,
        "delta3": delta3,
        "residue_states": residue_states,
        "k_states": k_states,
        "solved_l_states": solved_l_states,
        "interval_states": interval_states,
        "r3_states": r3_states,
        "prime_candidates": prime_candidates,
        "exact_solutions": exact_solutions,
        "exact_pairs": sorted(exact_pairs),
        "actual_found": actual_found,
        "l_residue_counter": l_residue_counter,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    rng = random.Random(SEED)

    banner(
        "THREE-CLOSE-PRIME CORRECTED QUOTIENT / "
        "LINEARIZED SECOND-MODULUS EXPERIMENT"
    )

    print(f"M                         = {M:,}")
    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"anchors                   = {ANCHORS}")
    print(
        f"modulus prime range       = "
        f"{MOD_MIN:,} - {MOD_MAX:,}"
    )
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")

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

    for i, (p, q) in enumerate(anchors, 1):

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
    # ALGEBRA VALIDATION
    # ============================================================================================

    banner(
        "VALIDATING CORRECT QUOTIENT IDENTITY"
    )

    identity_failures = 0

    for i, (n, p, q, mods) in enumerate(
        triples,
        1,
    ):

        r1, r2, r3 = mods

        try:

            state = validate_identity(
                n,
                p,
                q,
                r1,
                r2,
                r3,
            )

        except AssertionError as exc:

            identity_failures += 1

            print()
            print(
                f"FIRST FAILURE: anchor {i}"
            )

            print(
                f"n={n:,} p={p:,} q={q:,}"
            )

            print(
                f"mods=({r1},{r2},{r3})"
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
        f"{identity_failures}"
    )

    if identity_failures:
        raise RuntimeError(
            "Corrected identity still failed."
        )

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    search_triples = triples[
        :FULL_SEARCH_ANCHORS
    ]

    banner(
        "RUNNING CORRECTED LINEARIZED r2 SEARCH"
    )

    print(
        f"search anchors             = "
        f"{len(search_triples)}"
    )

    start = time.perf_counter()

    results = []

    for i, (n, p, q, mods) in enumerate(
        search_triples,
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
                f"anchor {i:3d}/{len(search_triples)}"
            )

    elapsed = time.perf_counter() - start

    if not results:
        raise RuntimeError(
            "No results were produced."
        )

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    banner("SUMMARY")

    mean_residue = statistics.mean(
        r["residue_states"]
        for _, _, _, r in results
    )

    mean_k = statistics.mean(
        r["k_states"]
        for _, _, _, r in results
    )

    mean_l = statistics.mean(
        r["solved_l_states"]
        for _, _, _, r in results
    )

    mean_interval = statistics.mean(
        r["interval_states"]
        for _, _, _, r in results
    )

    mean_r3 = statistics.mean(
        r["r3_states"]
        for _, _, _, r in results
    )

    mean_prime = statistics.mean(
        r["prime_candidates"]
        for _, _, _, r in results
    )

    mean_exact = statistics.mean(
        r["exact_solutions"]
        for _, _, _, r in results
    )

    print(
        f"average r1 residue states = "
        f"{mean_residue:.3f}"
    )

    print(
        f"average k states          = "
        f"{mean_k:.3f}"
    )

    print(
        f"average (a,k) states      = "
        f"{mean_residue * mean_k:.3f}"
    )

    print(
        f"average solved-l states   = "
        f"{mean_l:.3f}"
    )

    print(
        f"average interval states   = "
        f"{mean_interval:.3f}"
    )

    print(
        f"average r3 survivors      = "
        f"{mean_r3:.3f}"
    )

    print(
        f"average prime candidates  = "
        f"{mean_prime:.3f}"
    )

    print(
        f"average exact solutions   = "
        f"{mean_exact:.3f}"
    )

    # ============================================================================================
    # BASELINE
    # ============================================================================================

    banner("BASELINE COMPARISON")

    baseline = len(factor_primes)

    print(
        f"ordinary factor-prime enumeration = "
        f"{baseline:,}"
    )

    print(
        f"average solved-l states            = "
        f"{mean_l:.3f}"
    )

    print(
        f"average r3 states                  = "
        f"{mean_r3:.3f}"
    )

    print(
        f"average prime candidates           = "
        f"{mean_prime:.3f}"
    )

    print(
        f"solved-l / baseline                = "
        f"{mean_l / baseline:.9f}"
    )

    print(
        f"r3 / baseline                      = "
        f"{mean_r3 / baseline:.9f}"
    )

    print(
        f"prime / baseline                   = "
        f"{mean_prime / baseline:.9f}"
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
        f"correctly recovered       = "
        f"{recovered}/{len(results)}"
    )

    # ============================================================================================
    # DISTRIBUTIONS
    # ============================================================================================

    banner("SOLVED-l DISTRIBUTION")

    l_distribution = Counter(
        r["solved_l_states"]
        for _, _, _, r in results
    )

    for value, count in sorted(
        l_distribution.items()
    ):
        print(
            f"l candidates = {value:6d} "
            f"anchors = {count:3d}"
        )

    banner("R3 SURVIVOR DISTRIBUTION")

    r3_distribution = Counter(
        r["r3_states"]
        for _, _, _, r in results
    )

    for value, count in sorted(
        r3_distribution.items()
    ):
        print(
            f"r3 survivors = {value:6d} "
            f"anchors = {count:3d}"
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
            f"R/n={r['R_over_n']:.10f} "
            f"g={r['g']} "
            f"t={r['t']} "
            f"residue={r['residue_states']} "
            f"k={r['k_states']} "
            f"solved_l={r['solved_l_states']} "
            f"r3={r['r3_states']} "
            f"prime={r['prime_candidates']} "
            f"exact={r['exact_solutions']} "
            f"recovered={r['actual_found']}"
        )

    # ============================================================================================
    # STRONGEST COLLAPSES
    # ============================================================================================

    banner("STRONGEST LINEARIZED COLLAPSES")

    ranked = sorted(
        results,
        key=lambda x: (
            x[3]["r3_states"],
            x[3]["solved_l_states"],
        ),
    )

    for n, p, q, r in ranked[:20]:

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods=({r['r1']},{r['r2']},{r['r3']}) "
            f"R/n={r['R_over_n']:.10f} "
            f"solved_l={r['solved_l_states']} "
            f"r3={r['r3_states']} "
            f"prime={r['prime_candidates']} "
            f"exact={r['exact_solutions']}"
        )

    # ============================================================================================
    # ANCHOR TABLE
    # ============================================================================================

    banner("ANCHOR RESULTS")

    print(
        " ID        p        q      r1    r2    r3"
        "     (a,k)     L      R3    EXACT"
    )

    print("-" * 100)

    for i, (n, p, q, r) in enumerate(
        results,
        1,
    ):

        print(
            f"{i:3d} "
            f"{p:8,d} "
            f"{q:8,d} "
            f"{r['r1']:6d} "
            f"{r['r2']:6d} "
            f"{r['r3']:6d} "
            f"{r['residue_states'] * r['k_states']:9,d} "
            f"{r['solved_l_states']:6,d} "
            f"{r['r3_states']:6,d} "
            f"{r['exact_solutions']:6,d}"
        )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    banner("TIMING")

    print(
        f"search runtime            = "
        f"{elapsed:.3f} seconds"
    )

    # ============================================================================================
    # INTERPRETATION
    # ============================================================================================

    banner("MATHEMATICAL INTERPRETATION")

    print(
        """
The previous experiment contained a second algebra mistake.

It correctly defined:

    g = n mod r1

but then treated g as though it were also the residue of n modulo r2.

That is false in general.

The exact relation is:

    n = g + t*r1

where:

    t = (n-g)/r1.

At r1 we have:

    a*b = g + h*r1.

Therefore:

    a*b - n = (h-t)*r1.

With:

    p = a + k*r1
    q = b + l*r1
    Delta = r2-r1,

we obtain:

    p mod r2 = a-k*Delta
    q mod r2 = b-l*Delta.

Expanding the real r2 condition gives:

    (a-k*Delta)(b-l*Delta) = n (mod r2)

and therefore:

    (h-t)*r1
      - a*l*Delta
      - b*k*Delta
      + k*l*Delta^2
      = 0 (mod r2).

For fixed a,b,k this becomes:

    l*(k*Delta^2-a*Delta)
      =
    b*k*Delta-(h-t)*r1
      (mod r2).

So the close-modulus elimination survives, but only when the
quotient information of BOTH n and ab is retained.

The important measurements remain:

    (a,k) states
         ->
    solved l states
         ->
    third-modulus survivors
         ->
    exact factor pairs.

The experiment should now first pass the exact algebra validation
for all anchors. Only after that should the linearized search be
trusted.

The eventual goal is still to determine whether the remaining
(a,k) enumeration can itself be eliminated. If it cannot, then this
is a coordinate reduction rather than a fundamentally new factoring
algorithm.
"""
    )

    banner("EXPERIMENT COMPLETE")


if __name__ == "__main__":
    run()
