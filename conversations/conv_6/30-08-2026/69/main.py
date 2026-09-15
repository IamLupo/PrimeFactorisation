#!/usr/bin/env python3

import math
import random
import time

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

ANCHORS_PER_SCALE = 12

# Desired recursive object.
K_TARGET = 1000

# Search around n/K.
R_MULTIPLIERS = (
    0.80,
    0.85,
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
    1.15,
    1.20,
)

# Deterministic auxiliary constructions.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Number of candidate R values examined around n/K.
R_CANDIDATES = 31

# Prime-search radius around sqrt(R).
PAIR_WINDOW = 2500

# Keep the output small.
MAX_MATCH_EXAMPLES = 12


# =============================================================================
# PRIME UTILITIES
# =============================================================================

def next_prime(x: int) -> int:
    return int(sp.nextprime(max(2, x - 1)))


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def generate_anchor(
    scale: int,
    rng: random.Random,
):
    """
    Generate a balanced semiprime near the requested scale.
    """

    root = math.isqrt(scale)

    lo = int(root * 0.80)
    hi = int(root * 1.20)

    lo = max(lo, 100)

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

    return p, q, p * q


# =============================================================================
# BALANCED PRIME PAIR
# =============================================================================

def choose_balanced_prime_pair(
    target_R: int,
):
    """
    Find a balanced prime pair with product near target_R.
    """

    center = math.isqrt(
        max(target_R, 4)
    )

    best = None

    for delta in range(
        PAIR_WINDOW + 1
    ):

        for c in (
            center - delta,
            center + delta,
        ):

            if c < 2:
                continue

            r1 = next_prime(c)

            desired = max(
                2,
                target_R // r1,
            )

            for probe in (
                desired - 2,
                desired - 1,
                desired,
                desired + 1,
                desired + 2,
            ):

                if probe < 2:
                    continue

                r2 = next_prime(
                    probe
                )

                lo = min(r1, r2)
                hi = max(r1, r2)

                # Balanced moduli.
                if hi / lo > 1.30:
                    continue

                R = r1 * r2

                score = (
                    abs(R - target_R),
                    abs(r1 - r2),
                )

                if (
                    best is None
                    or score < best[0]
                ):
                    best = (
                        score,
                        r1,
                        r2,
                        R,
                    )

        if (
            best is not None
            and best[0][0] == 0
        ):
            break

    if best is None:
        raise RuntimeError(
            f"Could not find balanced prime pair "
            f"near R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# FIND UNIQUE R REGIMES
# =============================================================================

def build_R_targets(
    n: int,
):
    """
    Build several target R values around:

        R = n / K_TARGET

    plus a small local integer neighborhood.

    This tests whether the fixed-K target is sensitive to R.
    """

    base = n / K_TARGET

    targets = set()

    for multiplier in R_MULTIPLIERS:

        center = int(
            round(
                base * multiplier
            )
        )

        # Add a small deterministic local neighborhood.
        step = max(
            1,
            center // 500
        )

        for j in range(
            -R_CANDIDATES // 2,
            R_CANDIDATES // 2 + 1,
        ):

            targets.add(
                center + j * step
            )

    return sorted(
        x for x in targets
        if x > 3
    )


# =============================================================================
# FACTOR PAIRS
# =============================================================================

def factor_pairs(
    value: int,
):
    """
    Exact SymPy factorization and divisor-pair generation.
    """

    factors = sp.factorint(
        value
    )

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(
            divisors
        )

        powers = [1]
        cur = 1

        for _ in range(exponent):
            cur *= prime
            powers.append(cur)

        divisors = [
            d * p
            for d in old
            for p in powers
        ]

    divisors = sorted(
        set(divisors)
    )

    root = math.isqrt(
        value
    )

    pairs = []

    for d in divisors:

        if d > root:
            break

        if value % d == 0:
            pairs.append(
                (
                    d,
                    value // d,
                )
            )

    return factors, pairs


# =============================================================================
# AUXILIARY n+x
# =============================================================================

def auxiliary_test(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    K: int,
):
    """
    Construct deterministic n+x values and test Kx=K.
    """

    true_k = p // r1
    true_l = q // r2

    matches = 0
    same = 0
    cross = 0
    aux_pairs = 0
    S_hits = 0

    examples = []

    for S in S_VALUES:

        x = (-n) % S

        if x == 0:
            continue

        nx = n + x

        _, pairs = factor_pairs(
            nx
        )

        aux_pairs += len(
            pairs
        )

        local_match = False

        for px, qx in pairs:

            for a, b in (
                (px, qx),
                (qx, px),
            ):

                kx = a // r1
                lx = b // r2

                Kx = kx * lx

                if Kx == K:

                    matches += 1
                    local_match = True

                    if (
                        kx == true_k
                        and lx == true_l
                    ):
                        same += 1
                    else:
                        cross += 1

                    if len(examples) < 4:
                        examples.append(
                            (
                                S,
                                x,
                                a,
                                b,
                                kx,
                                lx,
                            )
                        )

        if local_match:
            S_hits += 1

    return {
        "aux_pairs": aux_pairs,
        "matches": matches,
        "same": same,
        "cross": cross,
        "S_hits": S_hits,
        "examples": examples,
    }


# =============================================================================
# R REGIME
# =============================================================================

def evaluate_R(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    R: int,
):
    """
    Evaluate one R choice.
    """

    k = p // r1
    l = q // r2

    if k <= 0 or l <= 0:
        return None

    K = k * l

    T = n // R
    E = T - K

    aux = auxiliary_test(
        n,
        p,
        q,
        r1,
        r2,
        K,
    )

    return {
        "r1": r1,
        "r2": r2,
        "R": R,

        "k": k,
        "l": l,
        "K": K,

        "T": T,
        "E": E,

        "target_error": abs(
            K - K_TARGET
        ),

        "aux_pairs":
            aux["aux_pairs"],

        "matches":
            aux["matches"],

        "same":
            aux["same"],

        "cross":
            aux["cross"],

        "S_hits":
            aux["S_hits"],

        "examples":
            aux["examples"],
    }


# =============================================================================
# ONE ANCHOR
# =============================================================================

def run_anchor(
    n: int,
    p: int,
    q: int,
):
    """
    Search several R values and retain the R regimes whose K is
    closest to K_TARGET.
    """

    candidates = []

    targets = build_R_targets(
        n
    )

    # Avoid repeating identical R values.
    seen_R = set()

    for target in targets:

        r1, r2, R = (
            choose_balanced_prime_pair(
                target
            )
        )

        if R in seen_R:
            continue

        seen_R.add(R)

        result = evaluate_R(
            n,
            p,
            q,
            r1,
            r2,
            R,
        )

        if result is not None:
            candidates.append(
                result
            )

    candidates.sort(
        key=lambda r: (
            r["target_error"],
            -r["matches"],
            abs(
                r["R"]
                - n / K_TARGET
            ),
        )
    )

    return candidates


# =============================================================================
# SCALE RUN
# =============================================================================

def run_scale(
    scale: int,
    rng: random.Random,
):
    rows = []

    print()
    print("=" * 100)
    print(
        f"SCALE {scale:.0e}"
    )
    print("=" * 100)

    for i in range(
        ANCHORS_PER_SCALE
    ):

        p, q, n = (
            generate_anchor(
                scale,
                rng,
            )
        )

        print(
            f"anchor {i + 1:2d}/"
            f"{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        regimes = run_anchor(
            n,
            p,
            q,
        )

        # Keep only the best three R choices for each anchor.
        for regime in regimes[:3]:

            rows.append(
                {
                    "scale": scale,
                    "p": p,
                    "q": q,
                    "n": n,
                    **regime,
                }
            )

    return rows


# =============================================================================
# SCALE SUMMARY
# =============================================================================

def print_scale_summary(
    rows,
):
    print()
    print("=" * 100)
    print(
        "BEST R REGIMES BY SCALE"
    )
    print("=" * 100)

    print(
        "scale       anchors   "
        "mean K       mean |K-1000|   "
        "mean E      auxPairs       "
        "Kx=K       same      cross"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r for r in rows
            if r["scale"] == scale
        ]

        if not group:
            continue

        mean_K = (
            sum(
                r["K"]
                for r in group
            )
            / len(group)
        )

        mean_error = (
            sum(
                r["target_error"]
                for r in group
            )
            / len(group)
        )

        mean_E = (
            sum(
                r["E"]
                for r in group
            )
            / len(group)
        )

        aux = sum(
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

        print(
            f"{scale:>5.0e} "
            f"{len(group):9d} "
            f"{mean_K:13.3f} "
            f"{mean_error:17.3f} "
            f"{mean_E:11.3f} "
            f"{aux:13,d} "
            f"{matches:10,d} "
            f"{same:9,d} "
            f"{cross:9,d}"
        )


# =============================================================================
# K TARGET QUALITY
# =============================================================================

def print_target_quality(
    rows,
):
    print()
    print("=" * 100)
    print(
        "K TARGET QUALITY"
    )
    print("=" * 100)

    print(
        "scale       K within 5    K within 10   "
        "K within 25   K within 50   "
        "K within 100"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r for r in rows
            if r["scale"] == scale
        ]

        if not group:
            continue

        def count(bound):
            return sum(
                abs(
                    r["K"]
                    - K_TARGET
                ) <= bound
                for r in group
            )

        print(
            f"{scale:>5.0e} "
            f"{count(5):13d} "
            f"{count(10):14d} "
            f"{count(25):14d} "
            f"{count(50):14d} "
            f"{count(100):15d}"
        )


# =============================================================================
# COLLISION SUMMARY
# =============================================================================

def print_collision_summary(
    rows,
):
    print()
    print("=" * 100)
    print(
        "Kx = K AUXILIARY COLLISIONS"
    )
    print("=" * 100)

    print(
        "scale       auxPairs       "
        "Kx=K        same       cross      "
        "match/aux"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r for r in rows
            if r["scale"] == scale
        ]

        aux = sum(
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

        ratio = (
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
            f"{ratio:14.10f}"
        )


# =============================================================================
# MATCH EXAMPLES
# =============================================================================

def print_matches(
    rows,
):
    print()
    print("=" * 100)
    print(
        "REPRESENTATIVE Kx = K MATCHES"
    )
    print("=" * 100)

    shown = 0

    for row in rows:

        if row["matches"] <= 0:
            continue

        ex = row[
            "examples"
        ]

        print(
            f"scale={row['scale']:.0e} "
            f"n={row['n']:,}"
        )

        print(
            f"    R={row['R']:,} "
            f"(r1,r2)=("
            f"{row['r1']},"
            f"{row['r2']})"
        )

        print(
            f"    true (k,l)=("
            f"{row['k']},"
            f"{row['l']}) "
            f"K={row['K']:,} "
            f"T={row['T']:,} "
            f"E={row['E']:,}"
        )

        print(
            f"    K-target error="
            f"{row['target_error']:,} "
            f"Kx=K={row['matches']} "
            f"same={row['same']} "
            f"cross={row['cross']}"
        )

        for (
            S,
            x,
            px,
            qx,
            kx,
            lx,
        ) in ex:

            print(
                f"    S={S} x={x} "
                f"aux=({px:,},{qx:,}) "
                f"(kx,lx)=({kx},{lx})"
            )

        print()

        shown += 1

        if shown >= MAX_MATCH_EXAMPLES:
            break

    if shown == 0:
        print("none")


# =============================================================================
# INTERPRETATION
# =============================================================================

def print_interpretation():

    print()
    print("=" * 100)
    print(
        "MATHEMATICAL INTERPRETATION"
    )
    print("=" * 100)

    print(
        r"""
The experiment tests a fixed recursive object size:

    K_target = 1000.

Starting from:

    n = p*q

and:

    R = r1*r2

define:

    k = floor(p/r1)
    l = floor(q/r2)

and therefore:

    K = k*l.

The quotient identity is:

    T = floor(n/R)
    E = T-K.

Ignoring the bounded floor/carry effects:

    K ~= n/R.

Therefore the natural modulus product is:

    R ~= n/K_target.

The experiment searches several nearby R values instead of assuming
that one exact modulus product is optimal.

For every R it computes:

    K = k*l.

The first question is:

    How accurately can the desired recursive size

        K ~= 1000

    actually be realized?

The second question is:

    Does a good R also improve the auxiliary trajectory?

For deterministic auxiliary values:

    x = (-n) mod S

we have:

    n+x == 0 mod S.

After exact factorization:

    n+x = px*qx

we calculate:

    kx = floor(px/r1)
    lx = floor(qx/r2)

and:

    Kx = kx*lx.

The stronger recursive event is:

    Kx = K.

It is separated into:

    same-cell:
        (kx,lx)=(k,l)

and:

    cross-cell:
        (kx,lx)!=(k,l)
        but
        kx*lx=K.

The main object of interest is therefore the combined chain:

    n
     |
     v
    R ~= n/1000
     |
     v
    T ~= 1000
     |
     v
    K = T-E ~= 1000
     |
     v
    factor(K)
     |
     v
    candidate (k,l)

together with:

    n+x
      |
      v
    factor(n+x)
      |
      v
    (kx,lx)
      |
      v
    Kx.

This experiment deliberately does not assume that the best R is
exactly n/1000. It searches a neighborhood and asks whether there
is a reproducible optimum.

The final p,q values are used only to score the resulting
quotient coordinates and auxiliary Kx collisions.

SymPy is used as the exact factorization oracle.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    print("=" * 100)
    print(
        "ADAPTIVE R = n/K FIXED-SIZE "
        "RECURSIVE EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"scales                 = "
        f"{[f'{s:.0e}' for s in SCALES]}"
    )

    print(
        f"anchors / scale        = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"K target               = "
        f"{K_TARGET}"
    )

    print(
        f"R multipliers          = "
        f"{R_MULTIPLIERS}"
    )

    print(
        f"S values               = "
        f"{S_VALUES}"
    )

    print(
        f"seed                   = "
        f"{SEED}"
    )

    rng = random.Random(
        SEED
    )

    all_rows = []

    for scale in SCALES:

        rows = run_scale(
            scale,
            rng,
        )

        all_rows.extend(
            rows
        )

    print_scale_summary(
        all_rows
    )

    print_target_quality(
        all_rows
    )

    print_collision_summary(
        all_rows
    )

    print_matches(
        all_rows
    )

    print_interpretation()

    elapsed = (
        time.perf_counter()
        - start
    )

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime = "
        f"{elapsed:.3f}s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
