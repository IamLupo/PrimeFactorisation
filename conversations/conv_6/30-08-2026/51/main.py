#!/usr/bin/env python3

"""
============================================================================================
NEARBY SEMIPRIME FACTOR-TRAJECTORY / E-INVARIANT EXPERIMENT
============================================================================================

Purpose
-------

Test the hypothesis:

    n   = p*q
    n+x = p_x*q_x

with small x, where (p_x,q_x) is a DIFFERENT prime factor pair.

For fixed moduli r1,r2 define:

    p  = a  + k*r1
    q  = b  + l*r2

and therefore:

    E = floor(n/(r1*r2)) - k*l.

For every deliberately constructed nearby semiprime n+x we measure:

    E_x - E_0
    k_x*l_x - k_0*l_0
    g_x - g_0

The central hypothesis is:

    if |x| < r1*r2 - g,

then floor((n+x)/(r1*r2)) remains constant, so:

    E_x - E_0 = -(k_x*l_x-k_0*l_0).

Therefore:

    E_x = E_0
        <=>
    k_x*l_x = k_0*l_0.

This experiment asks whether such equality occurs for genuinely
different nearby factorizations.

No factorization oracle is used for n+x.

Nearby semiprimes are deliberately constructed from prime
perturbations of the original factors.
"""

import math
import random
import statistics
import time


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

SEED = 1_511_464_998

N_ANCHORS = 100

P_MIN = 10_000
P_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

# We deliberately construct nearby factors:
DELTA_RADIUS = 300

# Desired |x| range:
X_RADIUS = 1_000

# Number of nearby semiprime trajectories retained per anchor:
TARGET_NEIGHBORS = 20


# ==========================================================================================
# PRIME GENERATION
# ==========================================================================================

def sieve(limit: int):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, is_prime in enumerate(sieve) if is_prime]


# ==========================================================================================
# ANCHORS
# ==========================================================================================

def build_anchors(primes, count, seed):
    rng = random.Random(seed)

    eligible = [
        p for p in primes
        if P_MIN <= p <= P_MAX
    ]

    anchors = []

    # Construct semiprimes with both factors in the requested range.
    while len(anchors) < count:
        p = rng.choice(eligible)
        q = rng.choice(eligible)

        if p > q:
            p, q = q, p

        if p == q:
            continue

        n = p * q

        anchors.append((p, q, n))

    return anchors


# ==========================================================================================
# MODULUS PAIRS
# ==========================================================================================

def build_close_pairs(modulus_primes):
    pairs = []

    for i, r1 in enumerate(modulus_primes):
        for r2 in modulus_primes[i + 1:]:
            ratio = abs(r2 - r1) / min(r1, r2)

            if ratio <= CLOSE_RATIO:
                pairs.append((r1, r2))

    return pairs


def choose_modulus_pair(rng, pairs):
    return rng.choice(pairs)


# ==========================================================================================
# BASIC UTILITIES
# ==========================================================================================

def is_prime_lookup(prime_set, x):
    return x in prime_set


def ceil_div(a, b):
    return -((-a) // b)


# ==========================================================================================
# E / COORDINATE CALCULATION
# ==========================================================================================

def coordinates(p, q, r1, r2):
    a = p % r1
    b = q % r2

    k = (p - a) // r1
    l = (q - b) // r2

    R = r1 * r2
    T = (p * q) // R

    E = T - k * l
    g = (p * q) % R

    return {
        "p": p,
        "q": q,
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "K": k * l,
        "T": T,
        "E": E,
        "g": g,
        "R": R,
    }


# ==========================================================================================
# DELIBERATELY CONSTRUCT NEARBY SEMIPRIMES
# ==========================================================================================

def build_nearby_factorizations(
    p,
    q,
    prime_set,
    delta_radius,
    x_radius,
    target_count,
):
    """
    Generate nearby semiprimes by replacing

        p -> p + dp
        q -> q + dq

    with both new factors prime.

    We retain cases where:

        |(p+dp)(q+dq) - p*q| <= x_radius

    and exclude the original pair.

    The construction is deliberately factor-aware so we can test
    nearby factor trajectories without needing to factor n+x.
    """

    n = p * q

    candidates = {}

    deltas_p = [
        dp for dp in range(-delta_radius, delta_radius + 1)
        if dp != 0 or True
    ]

    deltas_q = [
        dq for dq in range(-delta_radius, delta_radius + 1)
    ]

    # First build prime perturbation lists.
    p_variants = []
    q_variants = []

    for dp in deltas_p:
        px = p + dp

        if px < P_MIN or px > P_MAX:
            continue

        if px in prime_set:
            p_variants.append((dp, px))

    for dq in deltas_q:
        qx = q + dq

        if qx < P_MIN or qx > P_MAX:
            continue

        if qx in prime_set:
            q_variants.append((dq, qx))

    # Exact nearby product test.
    #
    # We use the identity:
    #
    # x = p*dq + q*dp + dp*dq
    #
    # which avoids a full multiplication until needed.
    for dp, px in p_variants:
        for dq, qx in q_variants:

            if px == p and qx == q:
                continue

            x = p * dq + q * dp + dp * dq

            if abs(x) > x_radius:
                continue

            nx = px * qx

            if nx != n + x:
                raise AssertionError("Nearby factor construction failed.")

            # Canonicalize factor order.
            f1, f2 = sorted((px, qx))

            key = (x, f1, f2)

            candidates[key] = {
                "x": x,
                "p": f1,
                "q": f2,
            }

    # Prefer x close to zero, then factor distance.
    result = sorted(
        candidates.values(),
        key=lambda z: (
            abs(z["x"]),
            abs(z["p"] - min(p, q)),
            abs(z["q"] - max(p, q)),
        )
    )

    return result[:target_count]


# ==========================================================================================
# SINGLE ANCHOR EXPERIMENT
# ==========================================================================================

def run_anchor(p, q, r1, r2):

    original = coordinates(p, q, r1, r2)

    nearby = build_nearby_factorizations(
        p,
        q,
        PRIME_SET,
        DELTA_RADIUS,
        X_RADIUS,
        TARGET_NEIGHBORS,
    )

    trajectory = []

    for item in nearby:

        x = item["x"]
        px = item["p"]
        qx = item["q"]

        current = coordinates(px, qx, r1, r2)

        trajectory.append({
            "x": x,
            "p": px,
            "q": qx,
            "E": current["E"],
            "K": current["K"],
            "k": current["k"],
            "l": current["l"],
            "g": current["g"],
            "R": current["R"],
            "dE": current["E"] - original["E"],
            "dK": current["K"] - original["K"],
            "same_E": current["E"] == original["E"],
            "same_K": current["K"] == original["K"],
            "different_factorization": (
                px != p or qx != q
            ),
        })

    return original, trajectory


# ==========================================================================================
# MAIN EXPERIMENT
# ==========================================================================================

def run():

    global PRIME_SET

    start_total = time.perf_counter()

    print("=" * 92)
    print("NEARBY SEMIPRIME FACTOR-TRAJECTORY / E-INVARIANT EXPERIMENT")
    print("=" * 92)
    print(f"N anchors                 = {N_ANCHORS}")
    print(f"factor range              = {P_MIN:,} - {P_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"delta radius              = {DELTA_RADIUS}")
    print(f"x radius                  = {X_RADIUS}")
    print(f"neighbors / anchor        = {TARGET_NEIGHBORS}")
    print(f"seed                      = {SEED}")
    print()

    # --------------------------------------------------------------------------------------
    # Prime pools
    # --------------------------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING PRIME POOLS")
    print("=" * 92)

    factor_primes = sieve(P_MAX)
    modulus_primes = [
        p for p in factor_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    PRIME_SET = set(factor_primes)

    print(f"factor primes             = {len([p for p in factor_primes if P_MIN <= p <= P_MAX]):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # --------------------------------------------------------------------------------------
    # Anchors
    # --------------------------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING ANCHORS")
    print("=" * 92)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")
    print()

    # --------------------------------------------------------------------------------------
    # Modulus pairs
    # --------------------------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 92)

    modulus_pairs = build_close_pairs(modulus_primes)

    print(f"close modulus pairs       = {len(modulus_pairs):,}")
    print()

    rng = random.Random(SEED ^ 0xC0FFEE)

    # --------------------------------------------------------------------------------------
    # Statistics
    # --------------------------------------------------------------------------------------

    all_dE = []
    all_dK = []

    same_E = 0
    same_K = 0
    different_factorizations = 0

    useful_same_E = []
    nontrivial_same_E = []

    total_neighbors = 0

    trajectory_rows = []

    # --------------------------------------------------------------------------------------
    # Run
    # --------------------------------------------------------------------------------------

    print("=" * 92)
    print("RUNNING NEARBY FACTOR TRAJECTORIES")
    print("=" * 92)

    for index, (p, q, n) in enumerate(anchors, 1):

        if index % 10 == 0 or index == 1:
            print(f"anchor {index:3d}/{len(anchors):3d}")

        r1, r2 = choose_modulus_pair(rng, modulus_pairs)

        original, trajectory = run_anchor(
            p,
            q,
            r1,
            r2,
        )

        total_neighbors += len(trajectory)

        for row in trajectory:

            all_dE.append(row["dE"])
            all_dK.append(row["dK"])

            if row["same_E"]:
                same_E += 1
                useful_same_E.append(
                    (
                        n,
                        p,
                        q,
                        row["x"],
                        row["p"],
                        row["q"],
                        r1,
                        r2,
                        original["k"],
                        original["l"],
                        row["k"],
                        row["l"],
                        original["E"],
                        row["E"],
                        original["K"],
                        row["K"],
                    )
                )

                if row["different_factorization"]:
                    nontrivial_same_E.append(useful_same_E[-1])

            if row["same_K"]:
                same_K += 1

            if row["different_factorization"]:
                different_factorizations += 1

            trajectory_rows.append(
                (
                    n,
                    p,
                    q,
                    r1,
                    r2,
                    row,
                )
            )

    # --------------------------------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("SUMMARY")
    print("=" * 92)

    print(f"anchors analyzed               = {len(anchors):,}")
    print(f"nearby semiprimes tested       = {total_neighbors:,}")
    print(f"different factorizations       = {different_factorizations:,}")

    if all_dE:
        print(f"average |E_x-E|               = {statistics.mean(abs(x) for x in all_dE):.6f}")
        print(f"maximum |E_x-E|               = {max(abs(x) for x in all_dE):,}")

    if all_dK:
        print(f"average |K_x-K|               = {statistics.mean(abs(x) for x in all_dK):.6f}")
        print(f"maximum |K_x-K|               = {max(abs(x) for x in all_dK):,}")

    print()

    print("=" * 92)
    print("E STABILITY")
    print("=" * 92)

    print(f"E_x = E                           = {same_E:,}/{total_neighbors:,}")
    print(f"fraction                         = {same_E / total_neighbors if total_neighbors else 0:.6%}")

    print(f"K_x = K                           = {same_K:,}/{total_neighbors:,}")
    print(f"fraction                         = {same_K / total_neighbors if total_neighbors else 0:.6%}")

    print()

    print("=" * 92)
    print("NONTRIVIAL E-INVARIANT CASES")
    print("=" * 92)

    print(
        "Cases where E_x = E while the factorization (p_x,q_x) differs "
        "from (p,q):"
    )

    print()

    for row in nontrivial_same_E[:20]:

        (
            n,
            p,
            q,
            x,
            px,
            qx,
            r1,
            r2,
            k0,
            l0,
            kx,
            lx,
            E0,
            Ex,
            K0,
            Kx,
        ) = row

        print(
            f"n={n:,} x={x:+d} "
            f"orig=({p:,},{q:,}) "
            f"new=({px:,},{qx:,}) "
            f"mods=({r1},{r2})"
        )

        print(
            f"    (k,l)=({k0},{l0}) -> ({kx},{lx}) "
            f"K={K0}->{Kx} "
            f"E={E0}->{Ex}"
        )

        print(
            f"    SAME_K={K0 == Kx} "
            f"SAME_E={E0 == Ex}"
        )

    if not nontrivial_same_E:
        print("No nontrivial E-invariant trajectories were found.")

    print()

    # --------------------------------------------------------------------------------------
    # By x
    # --------------------------------------------------------------------------------------

    print("=" * 92)
    print("E DELTA BY x")
    print("=" * 92)

    by_x = {}

    for (
        n,
        p,
        q,
        r1,
        r2,
        row,
    ) in trajectory_rows:

        x = row["x"]

        by_x.setdefault(
            x,
            []
        ).append(
            row["dE"]
        )

    for x in sorted(by_x):
        vals = by_x[x]

        print(
            f"x={x:+6d} "
            f"count={len(vals):4d} "
            f"mean_dE={statistics.mean(vals):9.4f} "
            f"same={sum(v == 0 for v in vals):4d}"
        )

    # --------------------------------------------------------------------------------------
    # Exact algebraic check
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("ALGEBRAIC CONSISTENCY CHECK")
    print("=" * 92)

    failures = 0
    trivial_fixed_t = 0

    for (
        n,
        p,
        q,
        r1,
        r2,
        row,
    ) in trajectory_rows:

        R = r1 * r2

        # Recompute directly.
        nx = n + row["x"]

        Tx = nx // R

        expected_dE = (
            Tx - row["K"]
        ) - (
            n // R - (
                (p % r1 and (p - p % r1) // r1)
                or (p - p % r1) // r1
            ) * (
                (q - q % r2) // r2
            )
        )

        # Simpler exact comparison using stored original coordinates
        # can be obtained from the trajectory input itself.
        #
        # The direct relationship for the same modulus pair is:
        #
        # E_x - E_0 = (T_x-T_0) - (K_x-K_0).
        #
        original = coordinates(p, q, r1, r2)

        lhs = row["dE"]
        rhs = (Tx - original["T"]) - (row["K"] - original["K"])

        if lhs != rhs:
            failures += 1

        if Tx == original["T"]:
            if row["dE"] == -row["dK"]:
                trivial_fixed_t += 1

    print(f"identity failures               = {failures:,}")
    print(f"fixed-T samples satisfying")
    print(f"    dE = -dK                    = {trivial_fixed_t:,}/{len(trajectory_rows):,}")

    # --------------------------------------------------------------------------------------
    # Most interesting trajectories
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("MOST INTERESTING TRAJECTORIES")
    print("=" * 92)

    ranked = sorted(
        trajectory_rows,
        key=lambda item: (
            abs(item[5]["dE"]),
            abs(item[5]["x"]),
        )
    )

    for (
        n,
        p,
        q,
        r1,
        r2,
        row,
    ) in ranked[:20]:

        print(
            f"n={n:,} x={row['x']:+d} "
            f"({p:,},{q:,})->({row['p']:,},{row['q']:,}) "
            f"mods=({r1},{r2}) "
            f"K={row['K']} "
            f"E={row['E']} "
            f"dE={row['dE']:+d}"
        )

    # --------------------------------------------------------------------------------------
    # Important interpretation
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 92)

    print(
        """
For the original factorization:

    n = p*q

write:

    R = r1*r2
    T = floor(n/R)
    E = T-k*l.

For a nearby deliberately constructed semiprime:

    n+x = p_x*q_x

we similarly have:

    E_x = floor((n+x)/R) - k_x*l_x.

Therefore exactly:

    E_x-E
      =
    [floor((n+x)/R)-floor(n/R)]
      -
    [k_x*l_x-k*l].

If:

    floor((n+x)/R) = floor(n/R),

then:

    E_x-E = -(k_x*l_x-k*l).

Hence:

    E_x = E
        <=> 
    k_x*l_x = k*l.

This experiment directly tests that equivalence on genuinely
different nearby factorizations.

The crucial distinction is:

    fixed p,q coordinates
        versus
    independently changed prime factors.

A nontrivial observation would therefore be:

    n+x = p_x*q_x
    (p_x,q_x) != (p,q)
    k_x*l_x = k*l
    E_x = E

for many nearby x.

That would indicate that E is not merely an artifact of one
particular coordinate representation and may reflect a local
invariant of nearby semiprime factorizations.

Conversely, if:

    E_x != E

almost whenever the factorization changes, then the apparent
stability of E is mostly the consequence of holding the original
k,l fixed.

One especially important quantity is:

    dK = k_x*l_x - k*l

because in the constant-T regime:

    dE = -dK.

That tells us exactly what E is tracking.

No factorization oracle is needed for n+x because the nearby
semiprimes are deliberately constructed from prime perturbations.
The final products are nevertheless checked exactly.
"""
    )

    elapsed = time.perf_counter() - start_total

    print()
    print("=" * 92)
    print("TIMING")
    print("=" * 92)
    print(f"total runtime                  = {elapsed:.3f} s")

    print()
    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    run()
