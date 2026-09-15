#!/usr/bin/env python3

import math
import random
import time
from collections import defaultdict

import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 1_511_464_998

SCALES = (
    10**9,
    10**12,
    10**16,
)

ANCHORS_PER_SCALE = 8

# Theoretical R regimes.
R_EXPONENTS = (
    0.25,
    1.0 / 3.0,
    0.50,
    2.0 / 3.0,
)

# Small multiplicative offsets.
R_OFFSETS = (
    0.90,
    1.00,
    1.10,
)

# Deterministic auxiliary constructions.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Modulus bounds.
MOD_MIN = 300
MOD_MAX = 3000

# r1/r2 relative closeness.
MAX_MOD_RATIO = 0.20

# Don't let one pathological integer create huge output/work.
MAX_AUX_PAIRS = 100_000


# =============================================================================
# MODULUS PRIME POOL
# =============================================================================

def build_modulus_primes():
    return list(
        sp.primerange(
            MOD_MIN,
            MOD_MAX + 1,
        )
    )


def build_modulus_pairs(primes):
    """
    Precompute all valid prime pairs once.

    We deliberately keep r1 <= r2 here.

    Every pair satisfies:

        |r2-r1| / r1 <= MAX_MOD_RATIO
    """

    pairs = []

    for i, r1 in enumerate(primes):

        for r2 in primes[i:]:

            ratio = (r2 - r1) / r1

            if ratio <= MAX_MOD_RATIO:
                R = r1 * r2

                pairs.append(
                    (
                        R,
                        r1,
                        r2,
                    )
                )

    pairs.sort(key=lambda x: x[0])

    return pairs


# =============================================================================
# FEASIBLE R SELECTION
# =============================================================================

def choose_modulus_pair(target_R, modulus_pairs):
    """
    Find the feasible modulus pair whose product is closest to target_R.

    This function always returns a valid pair.

    Important:
    The requested R may be impossible because r1,r2 are constrained to
    [MOD_MIN, MOD_MAX]. Therefore we clamp the target to the nearest
    feasible interval before searching.
    """

    min_R = modulus_pairs[0][0]
    max_R = modulus_pairs[-1][0]

    effective_target = min(
        max(
            target_R,
            min_R,
        ),
        max_R,
    )

    # Binary search for insertion position.
    lo = 0
    hi = len(modulus_pairs)

    while lo < hi:

        mid = (lo + hi) // 2

        if modulus_pairs[mid][0] < effective_target:
            lo = mid + 1
        else:
            hi = mid

    candidates = []

    if lo < len(modulus_pairs):
        candidates.append(modulus_pairs[lo])

    if lo > 0:
        candidates.append(modulus_pairs[lo - 1])

    best = min(
        candidates,
        key=lambda item: (
            abs(item[0] - effective_target),
            abs(item[1] - item[2]),
        ),
    )

    return {
        "target_R": target_R,
        "effective_target_R": effective_target,
        "R": best[0],
        "r1": best[1],
        "r2": best[2],
        "clamped": (
            target_R != effective_target
        ),
    }


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def generate_anchor(scale, rng):
    """
    Generate a semiprime near the requested scale.

    The factors are approximately sqrt(scale).
    """

    center = math.isqrt(scale)

    lo = max(
        100,
        int(center * 0.70),
    )

    hi = int(
        center * 1.30
    )

    p = int(
        sp.randprime(
            lo,
            hi,
        )
    )

    q = int(
        sp.randprime(
            lo,
            hi,
        )
    )

    while q == p:
        q = int(
            sp.randprime(
                lo,
                hi,
            )
        )

    n = p * q

    return p, q, n


# =============================================================================
# DIVISOR PAIRS
# =============================================================================

def divisor_pairs(value):
    """
    Exact factorization with SymPy.

    Returns:

        (d,e)

    with:

        d*e=value
        d <= e
    """

    factors = sp.factorint(value)

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(divisors)

        powers = [1]
        current = 1

        for _ in range(exponent):
            current *= prime
            powers.append(current)

        divisors = []

        for d in old:
            for power in powers:
                divisors.append(
                    d * power
                )

    divisors = sorted(set(divisors))

    root = math.isqrt(value)

    pairs = []

    for d in divisors:

        if d > root:
            break

        if value % d != 0:
            continue

        pairs.append(
            (
                d,
                value // d,
            )
        )

        if len(pairs) >= MAX_AUX_PAIRS:
            break

    return pairs


# =============================================================================
# DETERMINISTIC AUXILIARY VALUE
# =============================================================================

def deterministic_x(n, S):
    """
    Smallest x >= 0 such that:

        n+x == 0 mod S.
    """
    return (-n) % S


# =============================================================================
# AUXILIARY TEST
# =============================================================================

def test_auxiliary(
    n,
    p,
    q,
    r1,
    r2,
    S,
):
    """
    Factor n+x exactly and test Kx=K.
    """

    x = deterministic_x(
        n,
        S,
    )

    # x=0 would merely rediscover n.
    if x == 0:
        return {
            "x": 0,
            "nx": n,
            "pairs": 0,
            "matches": 0,
            "same": 0,
            "cross": 0,
            "min_delta": None,
        }

    nx = n + x

    if nx <= 0:
        return {
            "x": x,
            "nx": nx,
            "pairs": 0,
            "matches": 0,
            "same": 0,
            "cross": 0,
            "min_delta": None,
        }

    pairs = divisor_pairs(nx)

    k = p // r1
    l = q // r2

    K = k * l

    matches = 0
    same = 0
    cross = 0

    min_delta = None

    for px, qx in pairs:

        # Both orientations.
        for pa, qa in (
            (px, qx),
            (qx, px),
        ):

            kx = pa // r1
            lx = qa // r2

            Kx = kx * lx

            delta = abs(
                Kx - K
            )

            if (
                min_delta is None
                or delta < min_delta
            ):
                min_delta = delta

            if Kx == K:

                matches += 1

                if (
                    kx == k
                    and lx == l
                ):
                    same += 1
                else:
                    cross += 1

    return {
        "x": x,
        "nx": nx,
        "pairs": len(pairs),
        "matches": matches,
        "same": same,
        "cross": cross,
        "min_delta": min_delta,
    }


# =============================================================================
# ONE MODULUS REGIME
# =============================================================================

def run_regime(
    n,
    p,
    q,
    r1,
    r2,
):
    """
    Run all S values for one selected R.
    """

    k = p // r1
    l = q // r2
    K = k * l

    T = n // (r1 * r2)
    E = T - K

    total_pairs = 0
    total_matches = 0
    total_same = 0
    total_cross = 0

    S_hits = 0

    min_delta = None

    for S in S_VALUES:

        result = test_auxiliary(
            n,
            p,
            q,
            r1,
            r2,
            S,
        )

        total_pairs += result["pairs"]
        total_matches += result["matches"]
        total_same += result["same"]
        total_cross += result["cross"]

        if result["matches"] > 0:
            S_hits += 1

        if result["min_delta"] is not None:

            if (
                min_delta is None
                or result["min_delta"] < min_delta
            ):
                min_delta = result["min_delta"]

    return {
        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "E": E,

        "aux_pairs": total_pairs,
        "matches": total_matches,
        "same": total_same,
        "cross": total_cross,
        "S_hits": S_hits,
        "min_delta": min_delta,
    }


# =============================================================================
# SCALE EXECUTION
# =============================================================================

def run_scale(
    scale,
    modulus_pairs,
    rng,
):
    """
    Run all anchors and all R regimes for one scale.
    """

    rows = []

    print()
    print("=" * 100)
    print(f"SCALE {scale:.0e}")
    print("=" * 100)

    for anchor_id in range(
        1,
        ANCHORS_PER_SCALE + 1,
    ):

        p, q, n = generate_anchor(
            scale,
            rng,
        )

        print(
            f"anchor {anchor_id:2d}/{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        for exponent in R_EXPONENTS:

            base_R = int(
                round(
                    n ** exponent
                )
            )

            for offset in R_OFFSETS:

                target_R = int(
                    round(
                        base_R * offset
                    )
                )

                modulus = choose_modulus_pair(
                    target_R,
                    modulus_pairs,
                )

                regime = run_regime(
                    n,
                    p,
                    q,
                    modulus["r1"],
                    modulus["r2"],
                )

                rows.append(
                    {
                        "scale": scale,
                        "n": n,
                        "p": p,
                        "q": q,

                        "exponent": exponent,
                        "offset": offset,

                        "target_R": modulus[
                            "target_R"
                        ],

                        "effective_target_R": modulus[
                            "effective_target_R"
                        ],

                        "R": modulus["R"],
                        "r1": modulus["r1"],
                        "r2": modulus["r2"],

                        "clamped": modulus[
                            "clamped"
                        ],

                        **regime,
                    }
                )

    return rows


# =============================================================================
# SUMMARY BY R EXPONENT
# =============================================================================

def summarize_by_exponent(rows):

    grouped = defaultdict(list)

    for row in rows:
        grouped[
            round(
                row["exponent"],
                6,
            )
        ].append(row)

    print()
    print("-" * 100)
    print("R-SCALE SUMMARY")
    print("-" * 100)

    print(
        "exp       mean-R         mean-K       "
        "auxPairs        Kx=K      same     cross    "
        "S-hit%    clamp%"
    )

    print("-" * 100)

    for exponent in sorted(grouped):

        group = grouped[
            exponent
        ]

        mean_R = (
            sum(
                r["R"]
                for r in group
            )
            / len(group)
        )

        mean_K = (
            sum(
                r["K"]
                for r in group
            )
            / len(group)
        )

        aux_pairs = sum(
            r["aux_pairs"]
            for r in group
        )

        matches = sum(
            r["matches"]
            for r in group
        )

        same = sum(
            r["same"]
            for r in group
        )

        cross = sum(
            r["cross"]
            for r in group
        )

        S_hits = sum(
            r["S_hits"]
            for r in group
        )

        clamp_count = sum(
            1
            for r in group
            if r["clamped"]
        )

        total_S = (
            len(group)
            * len(S_VALUES)
        )

        s_hit_pct = (
            100.0
            * S_hits
            / total_S
            if total_S
            else 0.0
        )

        clamp_pct = (
            100.0
            * clamp_count
            / len(group)
            if group
            else 0.0
        )

        print(
            f"{exponent:0.3f} "
            f"{mean_R:14.3f} "
            f"{mean_K:14.3f} "
            f"{aux_pairs:12,d} "
            f"{matches:10,d} "
            f"{same:9,d} "
            f"{cross:9,d} "
            f"{s_hit_pct:8.2f}% "
            f"{clamp_pct:8.2f}%"
        )


# =============================================================================
# CROSS-SCALE SUMMARY
# =============================================================================

def summarize_scale(rows):

    if not rows:
        return

    scale = rows[0]["scale"]

    total_aux = sum(
        r["aux_pairs"]
        for r in rows
    )

    total_matches = sum(
        r["matches"]
        for r in rows
    )

    total_same = sum(
        r["same"]
        for r in rows
    )

    total_cross = sum(
        r["cross"]
        for r in rows
    )

    rate = (
        total_matches
        / total_aux
        if total_aux
        else 0.0
    )

    print()
    print(
        f"SCALE {scale:.0e} TOTAL"
    )
    print(
        f"    regimes              = {len(rows)}"
    )
    print(
        f"    auxiliary pairs      = {total_aux:,}"
    )
    print(
        f"    Kx == K              = {total_matches:,}"
    )
    print(
        f"    same-cell            = {total_same:,}"
    )
    print(
        f"    cross-cell           = {total_cross:,}"
    )
    print(
        f"    match / aux-pair     = {rate:.10f}"
    )


# =============================================================================
# REPRESENTATIVE MATCHES
# =============================================================================

def show_examples(rows):

    print()
    print("=" * 100)
    print("REPRESENTATIVE Kx = K MATCHES")
    print("=" * 100)

    shown = 0

    for row in rows:

        if row["matches"] <= 0:
            continue

        print(
            f"scale={row['scale']:.0e} "
            f"n={row['n']:,} "
            f"exp={row['exponent']:.3f} "
            f"R={row['R']:,} "
            f"(r1,r2)=({row['r1']},{row['r2']}) "
            f"K={row['K']:,} "
            f"matches={row['matches']} "
            f"same={row['same']} "
            f"cross={row['cross']}"
        )

        shown += 1

        if shown >= 10:
            break

    if shown == 0:
        print("none")


# =============================================================================
# FINAL INTERPRETATION
# =============================================================================

def print_interpretation():

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        r"""
The purpose of this experiment is to determine whether the
modulus product

    R = r1*r2

has a preferred scale for the recursive construction

    n
     |
     v
    n+x
     |
     v
   (px,qx)
     |
     v
   (kx,lx)
     |
     v
    Kx.

For the original factorization:

    p = a + k*r1
    q = b + l*r2

we have

    K = k*l

and

    T = floor(n/R)

with

    E = T-K.

Thus:

    K approximately n/R.

The experiment compares four theoretical regimes:

    R ~ n^(1/4)
    R ~ n^(1/3)
    R ~ n^(1/2)
    R ~ n^(2/3).

Their corresponding K scales are:

    n^(3/4)
    n^(2/3)
    n^(1/2)
    n^(1/3).

The last regime therefore produces the smallest K object.

For every auxiliary value n+x:

    n+x = px*qx

is factored exactly.

Then:

    kx = floor(px/r1)
    lx = floor(qx/r2)

and:

    Kx = kx*lx.

The principal event is:

    Kx = K.

This is split into:

    same-cell:
        (kx,lx) = (k,l)

and:

    cross-cell:
        (kx,lx) != (k,l)
        but
        kx*lx = k*l.

Only the second is a genuine cross-cell multiplicative
collision.

IMPORTANT:

The modulus range is no longer fixed in this experiment.

r1 and r2 are allowed to grow with n.

This is necessary because, for example, at n=10^16:

    sqrt(n) = 10^8,

which cannot be represented by r1,r2 <= 3000.

The actual prime pair is therefore chosen near each theoretical
R target.

If a theoretical target lies outside the realizable range of the
selected prime-pair construction, it is explicitly clamped.

The experiment therefore does not silently treat an impossible
R as if it were realizable.

The central question is:

    Does some R scale simultaneously give

        small K

        and

        non-negligible Kx=K frequency?

If the collision probability collapses as R grows, then making
K smaller may destroy the useful auxiliary relationship.

If a particular scaling regime preserves cross-cell collisions,
that regime becomes the next candidate for recursive investigation.

SymPy factorization is used only as an experimental oracle.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 100)
    print("TRUE MODULUS-SCALING / K-RECURSION EXPERIMENT")
    print("=" * 100)

    print(
        f"scales                  = "
        f"{[f'{x:.0e}' for x in SCALES]}"
    )

    print(
        f"anchors / scale         = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"R exponents             = "
        f"{R_EXPONENTS}"
    )

    print(
        f"R offsets               = "
        f"{R_OFFSETS}"
    )

    print(
        f"S values                = "
        f"{S_VALUES}"
    )

    print(
        f"modulus range           = "
        f"{MOD_MIN} .. {MOD_MAX}"
    )

    print(
        f"seed                    = "
        f"{SEED}"
    )

    print()
    print("=" * 100)
    print("BUILDING MODULUS PRIME POOL")
    print("=" * 100)

    primes = build_modulus_primes()

    print(
        f"modulus primes          = "
        f"{len(primes)}"
    )

    print()
    print("=" * 100)
    print("BUILDING VALID MODULUS PAIRS")
    print("=" * 100)

    modulus_pairs = build_modulus_pairs(
        primes
    )

    print(
        f"valid prime pairs       = "
        f"{len(modulus_pairs)}"
    )

    print(
        f"minimum feasible R      = "
        f"{modulus_pairs[0][0]:,}"
    )

    print(
        f"maximum feasible R      = "
        f"{modulus_pairs[-1][0]:,}"
    )

    rng = random.Random(
        SEED
    )

    start = time.perf_counter()

    all_rows = []

    for scale in SCALES:

        rows = run_scale(
            scale,
            modulus_pairs,
            rng,
        )

        all_rows.extend(rows)

        summarize_by_exponent(
            rows
        )

        summarize_scale(
            rows
        )

    show_examples(
        all_rows
    )

    # -------------------------------------------------------------------------
    # FINAL COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CROSS-SCALE COMPARISON")
    print("=" * 100)

    grouped = defaultdict(list)

    for row in all_rows:
        grouped[
            row["scale"]
        ].append(row)

    print(
        "scale       auxPairs       Kx=K       "
        "same       cross       match/aux"
    )

    print("-" * 100)

    for scale in sorted(grouped):

        rows = grouped[scale]

        aux = sum(
            r["aux_pairs"]
            for r in rows
        )

        matches = sum(
            r["matches"]
            for r in rows
        )

        same = sum(
            r["same"]
            for r in rows
        )

        cross = sum(
            r["cross"]
            for r in rows
        )

        rate = (
            matches / aux
            if aux
            else 0.0
        )

        print(
            f"{scale:>5.0e} "
            f"{aux:14,d} "
            f"{matches:10,d} "
            f"{same:10,d} "
            f"{cross:10,d} "
            f"{rate:16.10f}"
        )

    # -------------------------------------------------------------------------
    # TIMING
    # -------------------------------------------------------------------------

    elapsed = (
        time.perf_counter()
        - start
    )

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime              = "
        f"{elapsed:.3f}s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()