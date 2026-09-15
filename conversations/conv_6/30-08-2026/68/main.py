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

ANCHORS_PER_SCALE = 10

K_TARGET = 1000

R_OFFSETS = (
    0.95,
    1.00,
    1.05,
)

S_VALUES = (
    30,
    210,
    2310,
    30030,
)

PAIR_SEARCH_WINDOW = 5000

MAX_AUX_PAIRS = 100_000

MAX_EXAMPLES = 6


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
    Generate a balanced semiprime approximately equal to `scale`.

    IMPORTANT:
    The factor range scales with sqrt(scale).

    This allows:

        scale = 1e9
        scale = 1e12
        scale = 1e16

    without using the old 10,000..100,000 restriction.
    """

    root = math.isqrt(scale)

    # Approximately balanced factors.
    lo = max(
        101,
        int(root * 0.80),
    )

    hi = max(
        lo + 100,
        int(root * 1.20),
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
# BALANCED PRIME MODULUS PAIR
# =============================================================================

def choose_balanced_prime_pair(
    target_R: int,
):
    """
    Find r1,r2 such that:

        r1*r2 ~= target_R

    with r1 ~= r2.

    No fixed modulus ceiling is imposed.
    """

    if target_R < 4:
        return 2, 2, 4

    center = math.isqrt(
        target_R
    )

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

            r1 = next_prime(
                candidate
            )

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

                lo = min(
                    r1,
                    r2,
                )

                hi = max(
                    r1,
                    r2,
                )

                # Keep the pair reasonably balanced.
                if (
                    hi / lo
                    > 1.25
                ):
                    continue

                R = r1 * r2

                score = (
                    abs(
                        R - target_R
                    ),
                    abs(
                        r1 - r2
                    ),
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
            "Could not construct balanced prime pair "
            f"near R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# FACTOR PAIRS
# =============================================================================

def factor_pairs(
    value: int,
):
    """
    Exact factorization with SymPy.
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
        current = 1

        for _ in range(exponent):
            current *= prime
            powers.append(
                current
            )

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

    result = []

    for d in divisors:

        if d > root:
            break

        if value % d == 0:

            result.append(
                (
                    d,
                    value // d,
                )
            )

            if (
                len(result)
                >= MAX_AUX_PAIRS
            ):
                break

    return factors, result


# =============================================================================
# AUXILIARY CONSTRUCTION
# =============================================================================

def deterministic_x(
    n: int,
    S: int,
) -> int:
    """
        x = (-n) mod S

    guarantees:

        n+x == 0 mod S.
    """

    return (-n) % S


# =============================================================================
# AUXILIARY Kx TEST
# =============================================================================

def test_auxiliary(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    K_true: int,
    S: int,
):
    """
    Factor deterministic n+x and inspect all quotient products Kx.
    """

    x = deterministic_x(
        n,
        S,
    )

    if x == 0:
        return {
            "x": x,
            "nx": n,
            "pairs": 0,
            "matches": 0,
            "same": 0,
            "cross": 0,
            "min_delta": 0,
        }

    nx = n + x

    _, pairs = factor_pairs(
        nx
    )

    k_true = p // r1
    l_true = q // r2

    matches = 0
    same = 0
    cross = 0

    min_delta = None

    for px, qx in pairs:

        # Both orientations are valid for the quotient analysis.
        for a, b in (
            (px, qx),
            (qx, px),
        ):

            kx = a // r1
            lx = b // r2

            Kx = kx * lx

            delta = abs(
                Kx - K_true
            )

            if (
                min_delta is None
                or delta < min_delta
            ):
                min_delta = delta

            if Kx == K_true:

                matches += 1

                if (
                    kx == k_true
                    and lx == l_true
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
# ONE REGIME
# =============================================================================

def run_regime(
    p: int,
    q: int,
    n: int,
    target_R: int,
):
    """

    Calculate:

        R
        T
        E
        K
        k
        l

    and test deterministic n+x trajectories.
    """

    r1, r2, R = (
        choose_balanced_prime_pair(
            target_R
        )
    )

    k = p // r1
    l = q // r2

    # Important:
    #
    # If r1 > p or r2 > q, k/l can become zero.
    #
    # Such a regime is not useful for the recursive construction.
    if k <= 0 or l <= 0:
        return {
            "r1": r1,
            "r2": r2,
            "R": R,
            "k": k,
            "l": l,
            "K": 0,
            "T": n // R,
            "E": n // R,
            "candidate_count": 0,
            "aux_pairs": 0,
            "matches": 0,
            "same": 0,
            "cross": 0,
            "S_hits": 0,
            "min_delta": None,
            "invalid": True,
        }

    K = k * l

    T = n // R

    E = T - K

    # For this structural experiment we can use the exact carry bound
    # implied by the actual quotient coordinates to expose how small
    # the recursive object is.
    #
    # This is reported separately from the observable quantities.
    E_bound_structural = (
        k
        + l
        + 2
    )

    candidate_count = min(
        E_bound_structural,
        T,
    ) + 1

    aux_pairs = 0
    matches = 0
    same = 0
    cross = 0
    S_hits = 0

    min_delta = None

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

        aux_pairs += aux[
            "pairs"
        ]

        matches += aux[
            "matches"
        ]

        same += aux[
            "same"
        ]

        cross += aux[
            "cross"
        ]

        if aux[
            "matches"
        ] > 0:
            S_hits += 1

        if (
            aux["min_delta"]
            is not None
        ):

            if (
                min_delta is None
                or aux["min_delta"]
                < min_delta
            ):
                min_delta = (
                    aux["min_delta"]
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

        "candidate_count": candidate_count,

        "E_bound_structural":
            E_bound_structural,

        "aux_pairs": aux_pairs,
        "matches": matches,
        "same": same,
        "cross": cross,
        "S_hits": S_hits,

        "min_delta": min_delta,

        "invalid": False,
    }


# =============================================================================
# SCALE
# =============================================================================

def run_scale(
    scale: int,
    rng: random.Random,
):
    """

    Generate anchors and test:

        R ~= 0.95 n/K_target
        R ~= 1.00 n/K_target
        R ~= 1.05 n/K_target
    """

    rows = []

    print()
    print("=" * 100)
    print(
        f"SCALE {scale:.0e}"
    )
    print("=" * 100)

    for anchor_id in range(
        1,
        ANCHORS_PER_SCALE + 1,
    ):

        p, q, n = (
            generate_anchor(
                scale,
                rng,
            )
        )

        print(
            f"anchor {anchor_id:2d}/"
            f"{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        for offset in R_OFFSETS:

            target_R = int(
                round(
                    n
                    * offset
                    / K_TARGET
                )
            )

            regime = run_regime(
                p,
                q,
                n,
                target_R,
            )

            rows.append(
                {
                    "scale": scale,
                    "anchor": anchor_id,
                    "p": p,
                    "q": q,
                    "n": n,
                    "offset": offset,
                    "target_R":
                        target_R,
                    **regime,
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
    print(
        "FIXED-K TARGET SUMMARY"
    )
    print("-" * 100)

    print(
        "scale       cases   mean R        "
        "mean K      mean E      "
        "Kstates      Kx=K      same    cross"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r
            for r in rows
            if r["scale"] == scale
        ]

        if not group:
            continue

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

        mean_E = (
            sum(
                r["E"]
                for r in group
            )
            / len(group)
        )

        mean_states = (
            sum(
                r["candidate_count"]
                for r in group
            )
            / len(group)
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
            f"{len(group):8d} "
            f"{mean_R:15.3f} "
            f"{mean_K:12.3f} "
            f"{mean_E:11.3f} "
            f"{mean_states:12.3f} "
            f"{matches:10,d} "
            f"{same:8,d} "
            f"{cross:8,d}"
        )


# =============================================================================
# E / K SCALING
# =============================================================================

def print_scaling(
    rows,
):
    print()
    print("=" * 100)
    print(
        "E / K SCALING"
    )
    print("=" * 100)

    print(
        "scale       mean K        mean E       "
        "mean E/K       mean Kstates"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r
            for r in rows
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

        mean_E = (
            sum(
                r["E"]
                for r in group
            )
            / len(group)
        )

        mean_states = (
            sum(
                r["candidate_count"]
                for r in group
            )
            / len(group)
        )

        ratio = (
            mean_E / mean_K
            if mean_K
            else 0
        )

        print(
            f"{scale:>5.0e} "
            f"{mean_K:14.3f} "
            f"{mean_E:14.3f} "
            f"{ratio:14.8f} "
            f"{mean_states:14.3f}"
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
        "Kx = K COLLISION SUMMARY"
    )
    print("=" * 100)

    print(
        "scale       aux pairs       "
        "Kx=K        same       cross       "
        "match/aux"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r
            for r in rows
            if r["scale"] == scale
        ]

        if not group:
            continue

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
# REPRESENTATIVE EXAMPLES
# =============================================================================

def print_examples(
    rows,
):
    print()
    print("=" * 100)
    print(
        "REPRESENTATIVE CASES"
    )
    print("=" * 100)

    shown = 0

    for row in rows:

        if row["invalid"]:
            continue

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
            f"    (k,l)=("
            f"{row['k']},"
            f"{row['l']}) "
            f"K={row['K']:,} "
            f"T={row['T']:,} "
            f"E={row['E']:,}"
        )

        print(
            f"    target K={K_TARGET} "
            f"|K-target|="
            f"{abs(row['K'] - K_TARGET):,}"
        )

        print(
            f"    E-bound="
            f"{row['E_bound_structural']:,} "
            f"Kstates="
            f"{row['candidate_count']:,}"
        )

        print(
            f"    Kx=K="
            f"{row['matches']} "
            f"same="
            f"{row['same']} "
            f"cross="
            f"{row['cross']}"
        )

        print()

        shown += 1

        if shown >= MAX_EXAMPLES:
            break


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
The experiment fixes the recursive object size rather than fixing
the exponent alpha.

The target is:

    K_target = 1000.

Since:

    K = k*l

and:

    T = floor(n/R)
    E = T-K,

we approximately have:

    K ~= n/R.

Therefore choose:

    R ~= n/1000.

This gives:

    n=10^9
        R ~= 10^6
        K ~= 10^3

    n=10^12
        R ~= 10^9
        K ~= 10^3

    n=10^16
        R ~= 10^13
        K ~= 10^3.

For balanced factors and balanced moduli:

    p ~= q ~= sqrt(n)

and:

    r1 ~= r2 ~= sqrt(R).

Therefore:

    k ~= p/r1
    l ~= q/r2

and:

    k*l ~= n/R ~= 1000.

This is the key scaling construction.

The recursive chain being tested is:

    n
     |
     v
    choose R ~= n/1000
     |
     v
    T = floor(n/R) ~= 1000
     |
     v
    K = T-E
     |
     v
    factor(K)
     |
     v
    divisor pairs (k,l).

The experiment also constructs deterministic auxiliary values:

    x = (-n) mod S

so:

    n+x == 0 mod S.

The auxiliary value is exactly factored:

    n+x = px*qx

giving:

    kx = floor(px/r1)
    lx = floor(qx/r2)

and:

    Kx = kx*lx.

The stronger event remains:

    Kx = K.

This is separated into:

    same-cell:
        (kx,lx)=(k,l)

and:

    cross-cell:
        (kx,lx)!=(k,l)
        but
        kx*lx=K.

A particularly interesting result would be that the recursive
object remains around 1000 over all three scales and that the
observable K candidate space remains small.

Important distinction:

    R ~= n/1000

is a fixed-K construction.

It is not the same as choosing a fixed power:

    R=n^(5/6).

At different n, the corresponding alpha therefore changes.

The present experiment is designed specifically to test whether
holding the recursive problem size approximately constant is
more informative than holding its exponent constant.

The actual E and k,l values are reported for validation, but
the auxiliary n+x construction itself does not enumerate the
original p interval.

SymPy is used as an exact factorization oracle.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    print("=" * 100)
    print(
        "FIXED K=1000 / R=n/K "
        "RECURSIVE SCALING EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"scales             = "
        f"{[f'{s:.0e}' for s in SCALES]}"
    )

    print(
        f"anchors / scale    = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"K target           = "
        f"{K_TARGET}"
    )

    print(
        f"R offsets          = "
        f"{R_OFFSETS}"
    )

    print(
        f"S values           = "
        f"{S_VALUES}"
    )

    print(
        f"seed               = "
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

    print_summary(
        all_rows
    )

    print_scaling(
        all_rows
    )

    print_collision_summary(
        all_rows
    )

    print_examples(
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
