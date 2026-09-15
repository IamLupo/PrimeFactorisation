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

ANCHORS_PER_SCALE = 10

# We deliberately target:
#
#       K_target ~= 1000
#
# Since:
#
#       K ~= n / R
#
# choose:
#
#       R ~= n / K_target.
#
K_TARGET = 1000

# Test small variation around the target R.
R_OFFSETS = (
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
)

# Deterministic auxiliary constructions.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Balanced modulus requirement:
#
#       r1 ~= r2 ~= sqrt(R)
#
# This is important because p and q are also approximately balanced.
MAX_MOD_RATIO = 0.20

# Local search around sqrt(target_R).
PAIR_SEARCH_WINDOW = 10000

MAX_AUX_PAIRS = 100_000


# =============================================================================
# PRIME HELPERS
# =============================================================================

def next_prime(n: int) -> int:
    return int(sp.nextprime(max(2, n - 1)))


# =============================================================================
# BALANCED MODULUS PAIR
# =============================================================================

def choose_balanced_prime_pair(target_R: int):
    """
    Find primes r1,r2 satisfying approximately

        r1*r2 ~= target_R

    with r1 and r2 close to one another.
    """

    if target_R < 4:
        return 2, 2, 4

    center = math.isqrt(target_R)

    best = None

    for delta in range(
        PAIR_SEARCH_WINDOW + 1
    ):

        for candidate in (
            center - delta,
            center + delta,
        ):

            if candidate < 2:
                continue

            r1 = next_prime(candidate)

            desired_r2 = max(
                2,
                target_R // r1,
            )

            for probe in (
                desired_r2 - 2,
                desired_r2 - 1,
                desired_r2,
                desired_r2 + 1,
                desired_r2 + 2,
            ):

                if probe < 2:
                    continue

                r2 = next_prime(probe)

                lo = min(r1, r2)
                hi = max(r1, r2)

                ratio = (
                    (hi - lo) / lo
                )

                if ratio > MAX_MOD_RATIO:
                    continue

                R = r1 * r2

                abs_error = abs(
                    R - target_R
                )

                rel_error = (
                    abs_error / target_R
                )

                score = (
                    rel_error,
                    abs_error,
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
            and best[0][0] < 1e-12
        ):
            break

    if best is None:
        raise RuntimeError(
            f"Could not construct balanced prime pair "
            f"for target R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def generate_anchor(
    scale: int,
    rng: random.Random,
):
    """
    Generate a balanced semiprime around the requested scale.
    """

    center = math.isqrt(scale)

    lo = max(
        100,
        int(center * 0.75),
    )

    hi = int(
        center * 1.25
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

    return p, q, p * q


# =============================================================================
# FACTOR PAIRS
# =============================================================================

def factor_pairs(value: int):
    """
    Exact factorization using SymPy.
    """

    fac = sp.factorint(
        value
    )

    divisors = [1]

    for prime, exponent in fac.items():

        old = list(
            divisors
        )

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
# DETERMINISTIC AUXILIARY
# =============================================================================

def deterministic_x(
    n: int,
    S: int,
) -> int:

    return (-n) % S


# =============================================================================
# ONE AUXILIARY TEST
# =============================================================================

def test_auxiliary(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    K: int,
    S: int,
):
    """
    Construct n+x deterministically with

        n+x == 0 mod S

    then factor it and calculate Kx.
    """

    x = deterministic_x(
        n,
        S,
    )

    if x == 0:
        return {
            "x": x,
            "pairs": 0,
            "matches": 0,
            "same": 0,
            "cross": 0,
            "min_delta": None,
        }

    nx = n + x

    pairs = factor_pairs(
        nx
    )

    k = p // r1
    l = q // r2

    matches = 0
    same = 0
    cross = 0
    min_delta = None

    for px, qx in pairs:

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
        "pairs": len(pairs),
        "matches": matches,
        "same": same,
        "cross": cross,
        "min_delta": min_delta,
    }


# =============================================================================
# ONE R TEST
# =============================================================================

def run_R_test(
    p: int,
    q: int,
    n: int,
    r1: int,
    r2: int,
):
    """
    Run all deterministic S constructions for one R.
    """

    R = r1 * r2

    k = p // r1
    l = q // r2

    K = k * l

    T = n // R
    E = T - K

    result = {
        "R": R,
        "r1": r1,
        "r2": r2,

        "k": k,
        "l": l,
        "K": K,
        "T": T,
        "E": E,

        "K_error": K - K_TARGET,
        "K_rel_error": (
            abs(K - K_TARGET)
            / K_TARGET
        ),

        "aux_tests": 0,
        "aux_pairs": 0,

        "matches": 0,
        "same": 0,
        "cross": 0,

        "S_hits": 0,

        "min_delta": None,
    }

    for S in S_VALUES:

        aux = test_auxiliary(
            n,
            p,
            q,
            r1,
            r2,
            K,
            S,
        )

        result[
            "aux_tests"
        ] += 1

        result[
            "aux_pairs"
        ] += aux[
            "pairs"
        ]

        result[
            "matches"
        ] += aux[
            "matches"
        ]

        result[
            "same"
        ] += aux[
            "same"
        ]

        result[
            "cross"
        ] += aux[
            "cross"
        ]

        if aux[
            "matches"
        ] > 0:
            result[
                "S_hits"
            ] += 1

        if (
            aux[
                "min_delta"
            ] is not None
        ):

            if (
                result[
                    "min_delta"
                ] is None
                or aux[
                    "min_delta"
                ]
                < result[
                    "min_delta"
                ]
            ):
                result[
                    "min_delta"
                ] = aux[
                    "min_delta"
                ]

    return result


# =============================================================================
# ONE SCALE
# =============================================================================

def run_scale(
    scale: int,
    rng: random.Random,
):

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
            f"anchor {anchor_id:2d}/"
            f"{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        # The central target:
        #
        #     R = n / 1000
        #
        target_base = (
            n / K_TARGET
        )

        for offset in R_OFFSETS:

            target_R = int(
                round(
                    target_base * offset
                )
            )

            r1, r2, R = (
                choose_balanced_prime_pair(
                    target_R
                )
            )

            result = run_R_test(
                p,
                q,
                n,
                r1,
                r2,
            )

            rows.append(
                {
                    "scale": scale,
                    "anchor": anchor_id,

                    "p": p,
                    "q": q,
                    "n": n,

                    "target_R": target_R,
                    "offset": offset,

                    **result,
                }
            )

    return rows


# =============================================================================
# SUMMARY
# =============================================================================

def print_summary(
    rows,
):

    print()
    print("-" * 100)
    print("K-TARGET SUMMARY")
    print("-" * 100)

    grouped = defaultdict(list)

    for row in rows:
        grouped[
            row["offset"]
        ].append(row)

    print(
        "offset      mean R          mean K        "
        "mean |K-1000|    auxPairs      Kx=K      "
        "same      cross"
    )

    print("-" * 100)

    for offset in sorted(
        grouped
    ):

        group = grouped[
            offset
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

        mean_error = (
            sum(
                abs(
                    r["K"]
                    - K_TARGET
                )
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
            f"{offset:0.2f} "
            f"{mean_R:15.3f} "
            f"{mean_K:15.3f} "
            f"{mean_error:16.3f} "
            f"{aux:12,d} "
            f"{matches:10,d} "
            f"{same:9,d} "
            f"{cross:9,d}"
        )


# =============================================================================
# SCALE SUMMARY
# =============================================================================

def print_scale_summary(
    rows,
):

    print()
    print("=" * 100)
    print("PER-SCALE SUMMARY")
    print("=" * 100)

    grouped = defaultdict(list)

    for row in rows:
        grouped[
            row["scale"]
        ].append(row)

    print(
        "scale       mean K        mean E       "
        "mean |K-1000|      Kx=K      cross       "
        "match/aux"
    )

    print("-" * 100)

    for scale in sorted(
        grouped
    ):

        group = grouped[
            scale
        ]

        mean_K = (
            sum(
                r["K"]
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

        mean_error = (
            sum(
                abs(
                    r["K"]
                    - K_TARGET
                )
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
            f"{mean_K:14.3f} "
            f"{mean_E:13.3f} "
            f"{mean_error:17.3f} "
            f"{matches:10,d} "
            f"{cross:10,d} "
            f"{ratio:14.10f}"
        )


# =============================================================================
# BEST R VALUES
# =============================================================================

def print_best_R(rows):

    print()
    print("=" * 100)
    print("CLOSEST R TO THE K=1000 TARGET")
    print("=" * 100)

    grouped = defaultdict(list)

    for row in rows:

        grouped[
            (
                row["scale"],
                row["anchor"],
            )
        ].append(row)

    shown = 0

    for key in sorted(
        grouped
    ):

        group = grouped[
            key
        ]

        best = min(
            group,
            key=lambda r: abs(
                r["K"]
                - K_TARGET
            ),
        )

        print(
            f"scale={best['scale']:.0e} "
            f"n={best['n']:,} "
            f"R={best['R']:,} "
            f"(r1,r2)=({best['r1']},{best['r2']}) "
            f"K={best['K']:,} "
            f"E={best['E']:,}"
        )

        shown += 1

        if shown >= 12:
            break


# =============================================================================
# REPRESENTATIVE COLLISIONS
# =============================================================================

def print_matches(rows):

    print()
    print("=" * 100)
    print("REPRESENTATIVE Kx = K COLLISIONS")
    print("=" * 100)

    count = 0

    for row in rows:

        if row[
            "matches"
        ] == 0:
            continue

        print(
            f"scale={row['scale']:.0e} "
            f"n={row['n']:,} "
            f"R={row['R']:,} "
            f"(r1,r2)=({row['r1']},{row['r2']}) "
            f"K={row['K']:,} "
            f"matches={row['matches']} "
            f"same={row['same']} "
            f"cross={row['cross']}"
        )

        count += 1

        if count >= 12:
            break

    if count == 0:
        print("none")


# =============================================================================
# MATHEMATICAL INTERPRETATION
# =============================================================================

def print_interpretation():

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        r"""
This experiment fixes the recursive target size rather than fixing
the modulus range.

The target is:

    K_target = 1000.

Since:

    K = k*l

and:

    T = floor(n/R)
    E = T-K,

we have approximately:

    K ~= n/R.

Therefore choose:

    R ~= n/1000.

This gives the intended hierarchy:

    n = 10^9
        R ~= 10^6
        K ~= 10^3

    n = 10^12
        R ~= 10^9
        K ~= 10^3

    n = 10^16
        R ~= 10^13
        K ~= 10^3.

Balanced moduli satisfy approximately:

    r1 ~= r2 ~= sqrt(R).

For balanced semiprimes:

    p ~= q ~= sqrt(n).

Therefore:

    k ~= p/r1
    l ~= q/r2

and:

    k*l ~= n/(r1*r2)
        ~= 1000.

This is the central construction being tested.

The experiment does NOT require the same exponent alpha at every
scale.

Instead it holds the recursive object size fixed:

    K ~= 1000.

That means:

    R ~= n/K.

This is a different hypothesis from:

    R = n^(5/6).

At n=10^16, for example:

    R ~= 10^13

produces:

    K ~= 1000.

The auxiliary process is:

    n+x
       |
       v
    factor(n+x)
       |
       v
    (kx,lx)
       |
       v
    Kx = kx*lx.

The key event is:

    Kx = K.

It is divided into:

    same-cell:
        (kx,lx)=(k,l)

and:

    cross-cell:
        (kx,lx)!=(k,l)
        but
        kx*lx=K.

The experiment therefore asks:

    Can the original factorization problem

        factor(n)

    be represented by an auxiliary family whose recursive
    quotient-product object remains only about 1000 in size?

A particularly interesting outcome would be that K remains near
1000 over all three scales while Kx=K collisions remain observable.

The next conceptual step would then be to determine whether the
correct K can be identified from observable quantities

    T = floor(n/R)

and the structure of the defect E,

without already knowing k,l.

SymPy is used here only as an exact factorization oracle.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 100)
    print("FIXED K-TARGET / R = n/K RECURSIVE EXPERIMENT")
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
        f"K target                = "
        f"{K_TARGET}"
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
        f"modulus selection       = "
        f"balanced primes, no fixed upper bound"
    )

    print(
        f"seed                    = "
        f"{SEED}"
    )

    rng = random.Random(
        SEED
    )

    start = time.perf_counter()

    all_rows = []

    for scale in SCALES:

        rows = run_scale(
            scale,
            rng,
        )

        all_rows.extend(
            rows
        )

    print_summary(
        all_rows
    )

    print_scale_summary(
        all_rows
    )

    print_best_R(
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
        f"total runtime              = "
        f"{elapsed:.3f}s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
