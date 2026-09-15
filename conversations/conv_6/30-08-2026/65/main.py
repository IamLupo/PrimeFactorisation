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

# The interesting high-R regimes.
R_EXPONENTS = (
    0.50,
    2.0 / 3.0,
    0.75,
    5.0 / 6.0,
    11.0 / 12.0,
)

# Small variation around each theoretical R.
R_OFFSETS = (
    0.95,
    1.00,
    1.05,
)

# Deterministic auxiliary constructions.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# r1 and r2 should be relatively close.
MAX_MOD_RATIO = 0.20

# Number of nearby integer values to inspect around the theoretical
# sqrt(R) when searching for balanced prime moduli.
PAIR_SEARCH_WINDOW = 3000

# Limit auxiliary divisor-pair enumeration.
MAX_AUX_PAIRS = 100_000

# For each R regime, only keep a tiny number of representative matches.
MAX_MATCH_EXAMPLES = 2


# =============================================================================
# PRIME HELPERS
# =============================================================================

def next_prime(n: int) -> int:
    return int(sp.nextprime(max(2, n - 1)))


def previous_prime(n: int) -> int:
    if n <= 2:
        return 2
    return int(sp.prevprime(n + 1))


# =============================================================================
# MODULUS PAIR SELECTION
# =============================================================================

def choose_balanced_prime_pair(target_R: int):
    """
    Choose r1,r2 primes such that:

        r1*r2 ~= target_R

    while keeping r1 and r2 close.

    The search is deliberately local around sqrt(target_R), which
    is the balanced configuration.
    """

    if target_R < 4:
        return 2, 2

    center = math.isqrt(target_R)

    best = None

    for delta in range(PAIR_SEARCH_WINDOW + 1):

        candidates = []

        a = center - delta
        b = center + delta

        if a >= 2:
            candidates.append(a)

        candidates.append(b)

        for candidate in candidates:

            r1 = next_prime(candidate)

            desired = max(
                2,
                target_R // r1,
            )

            probes = (
                desired - 2,
                desired - 1,
                desired,
                desired + 1,
                desired + 2,
            )

            for probe in probes:

                if probe < 2:
                    continue

                r2 = next_prime(probe)

                lo = min(r1, r2)
                hi = max(r1, r2)

                ratio = (hi - lo) / lo

                if ratio > MAX_MOD_RATIO:
                    continue

                R = r1 * r2

                absolute_error = abs(R - target_R)

                relative_error = (
                    absolute_error / target_R
                )

                key = (
                    relative_error,
                    absolute_error,
                    abs(r1 - r2),
                )

                if best is None or key < best[0]:

                    best = (
                        key,
                        r1,
                        r2,
                        R,
                    )

        if best is not None:

            # Extremely good match.
            if best[0][0] < 1e-10:
                break

    if best is None:
        raise RuntimeError(
            f"Could not construct balanced prime pair near "
            f"target R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# SEMIPRIME GENERATION
# =============================================================================

def generate_anchor(scale: int, rng: random.Random):

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

    n = p * q

    return p, q, n


# =============================================================================
# DIVISOR PAIRS
# =============================================================================

def factor_pairs(value: int):

    factors = sp.factorint(value)

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(divisors)

        powers = [1]
        current = 1

        for _ in range(exponent):
            current *= prime
            powers.append(current)

        new_divisors = []

        for d in old:
            for power in powers:
                new_divisors.append(
                    d * power
                )

        divisors = new_divisors

    divisors = sorted(set(divisors))

    root = math.isqrt(value)

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

            if len(pairs) >= MAX_AUX_PAIRS:
                break

    return pairs


# =============================================================================
# DETERMINISTIC AUXILIARY
# =============================================================================

def make_x(n: int, S: int) -> int:
    return (-n) % S


# =============================================================================
# E / K CANDIDATE WINDOW
# =============================================================================

def compute_E_bound(k: int, l: int) -> int:
    """
    Safe bound from the carry structure:

        E = c1+c2+c3

    with approximately:

        c1 < l
        c2 < k
        c3 <= 2.

    Therefore E <= k+l (conservative enough here).
    """

    return k + l + 2


def observable_K_candidates(n: int, R: int, k: int, l: int):
    """
    Generate the observable T-E candidate window.

    TRUE E is deliberately not used to construct the list.
    The bound comes from k,l as a controlled experimental oracle.

    This is useful for measuring the size of the recursive K problem.
    """

    T = n // R

    E_max = compute_E_bound(
        k,
        l,
    )

    candidates = []

    for E in range(
        0,
        min(
            E_max,
            T,
        ) + 1,
    ):

        K = T - E

        if K <= 0:
            break

        candidates.append(
            K
        )

    return T, E_max, candidates


# =============================================================================
# AUXILIARY K TEST
# =============================================================================

def evaluate_auxiliary(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    K: int,
    S: int,
):

    x = make_x(
        n,
        S,
    )

    # Skip the trivial n+x=n case.
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

    if nx <= 0:
        return {
            "x": x,
            "pairs": 0,
            "matches": 0,
            "same": 0,
            "cross": 0,
            "min_delta": None,
        }

    pairs = factor_pairs(nx)

    k = p // r1
    l = q // r2

    matches = 0
    same = 0
    cross = 0

    min_delta = None

    for px, qx in pairs:

        # Test both orientations.
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
# ONE REGIME
# =============================================================================

def run_regime(
    p: int,
    q: int,
    n: int,
    r1: int,
    r2: int,
):

    R = r1 * r2

    k = p // r1
    l = q // r2

    K = k * l

    T = n // R
    E = T - K

    T_obs, E_bound, K_candidates = (
        observable_K_candidates(
            n,
            R,
            k,
            l,
        )
    )

    result = {
        "R": R,
        "r1": r1,
        "r2": r2,

        "k": k,
        "l": l,

        "K": K,
        "T": T,
        "E": E,

        "E_bound": E_bound,

        "K_candidate_count": len(
            K_candidates
        ),

        "aux_tests": 0,
        "aux_pairs": 0,

        "matches": 0,
        "same": 0,
        "cross": 0,

        "S_hits": 0,

        "min_delta": None,

        "examples": [],
    }

    for S in S_VALUES:

        result["aux_tests"] += 1

        aux = evaluate_auxiliary(
            n,
            p,
            q,
            r1,
            r2,
            K,
            S,
        )

        result["aux_pairs"] += aux[
            "pairs"
        ]

        result["matches"] += aux[
            "matches"
        ]

        result["same"] += aux[
            "same"
        ]

        result["cross"] += aux[
            "cross"
        ]

        if aux["matches"] > 0:

            result["S_hits"] += 1

            if len(
                result["examples"]
            ) < MAX_MATCH_EXAMPLES:

                result["examples"].append(
                    {
                        "S": S,
                        "x": aux["x"],
                        "pairs": aux["pairs"],
                        "matches": aux["matches"],
                        "same": aux["same"],
                        "cross": aux["cross"],
                    }
                )

        if (
            aux["min_delta"] is not None
        ):

            if (
                result["min_delta"] is None
                or aux["min_delta"]
                < result["min_delta"]
            ):

                result["min_delta"] = (
                    aux["min_delta"]
                )

    return result


# =============================================================================
# SCALE RUN
# =============================================================================

def run_scale(
    scale: int,
    modulus_pairs,
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

                r1, r2, R = (
                    choose_balanced_prime_pair(
                        target_R
                    )
                )

                regime = run_regime(
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
                        "n": n,
                        "p": p,
                        "q": q,

                        "exponent": exponent,
                        "offset": offset,

                        "target_R": target_R,

                        **regime,
                    }
                )

        print(
            f"    anchor "
            f"{anchor_id:2d}/{ANCHORS_PER_SCALE}"
        )

    return rows


# =============================================================================
# SUMMARY
# =============================================================================

def summarize(rows):

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
        "exp       mean R           mean K       "
        "mean E     mean Ebound     Kstates    "
        "auxPairs       Kx=K      same    cross"
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

        mean_E = (
            sum(
                r["E"]
                for r in group
            )
            / len(group)
        )

        mean_Ebound = (
            sum(
                r["E_bound"]
                for r in group
            )
            / len(group)
        )

        mean_Kstates = (
            sum(
                r["K_candidate_count"]
                for r in group
            )
            / len(group)
        )

        auxPairs = sum(
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
            f"{exponent:0.3f} "
            f"{mean_R:16.3f} "
            f"{mean_K:16.3f} "
            f"{mean_E:10.3f} "
            f"{mean_Ebound:14.3f} "
            f"{mean_Kstates:10.3f} "
            f"{auxPairs:12,d} "
            f"{matches:10,d} "
            f"{same:8,d} "
            f"{cross:8,d}"
        )


# =============================================================================
# CROSS-SCALE
# =============================================================================

def cross_scale(rows):

    grouped = defaultdict(list)

    for row in rows:
        grouped[
            row["scale"]
        ].append(row)

    print()
    print("=" * 100)
    print("CROSS-SCALE COMPARISON")
    print("=" * 100)

    print(
        "scale       regimes      mean K       "
        "mean E       Kstates      Kx=K      "
        "same     cross     match/aux"
    )

    print("-" * 100)

    for scale in sorted(grouped):

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

        mean_states = (
            sum(
                r["K_candidate_count"]
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

        ratio = (
            matches / aux
            if aux
            else 0.0
        )

        print(
            f"{scale:>5.0e} "
            f"{len(group):12,d} "
            f"{mean_K:14.3f} "
            f"{mean_E:12.3f} "
            f"{mean_states:12.3f} "
            f"{matches:10,d} "
            f"{same:8,d} "
            f"{cross:8,d} "
            f"{ratio:14.10f}"
        )


# =============================================================================
# REPRESENTATIVE MATCHES
# =============================================================================

def show_matches(rows):

    print()
    print("=" * 100)
    print("REPRESENTATIVE Kx = K EVENTS")
    print("=" * 100)

    shown = 0

    for row in rows:

        if row["matches"] <= 0:
            continue

        for example in row["examples"]:

            print(
                f"scale={row['scale']:.0e} "
                f"n={row['n']:,} "
                f"exp={row['exponent']:.3f} "
                f"R={row['R']:,} "
                f"(r1,r2)=({row['r1']},{row['r2']})"
            )

            print(
                f"    K={row['K']:,} "
                f"T={row['T']:,} "
                f"E={row['E']:,} "
                f"Ebound={row['E_bound']:,} "
                f"Kstates={row['K_candidate_count']}"
            )

            print(
                f"    S={example['S']} "
                f"x={example['x']:,} "
                f"auxPairs={example['pairs']} "
                f"matches={example['matches']} "
                f"same={example['same']} "
                f"cross={example['cross']}"
            )

            shown += 1

            if shown >= 12:
                return

    if shown == 0:
        print("none")


# =============================================================================
# E-WINDOW SCALING
# =============================================================================

def summarize_E_scaling(rows):

    print()
    print("=" * 100)
    print("E / K SCALING")
    print("=" * 100)

    grouped = defaultdict(list)

    for row in rows:
        grouped[
            (
                row["scale"],
                round(
                    row["exponent"],
                    6,
                ),
            )
        ].append(row)

    print(
        "scale       exp      mean K       "
        "mean E       E/K          Kstates"
    )

    print("-" * 100)

    for (
        scale,
        exponent,
    ), group in sorted(grouped.items()):

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
                r["K_candidate_count"]
                for r in group
            )
            / len(group)
        )

        ratio = (
            mean_E / mean_K
            if mean_K
            else 0.0
        )

        print(
            f"{scale:>5.0e} "
            f"{exponent:8.3f} "
            f"{mean_K:14.3f} "
            f"{mean_E:14.3f} "
            f"{ratio:12.8f} "
            f"{mean_states:12.3f}"
        )


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
This experiment specifically investigates the high-R hypothesis.

Set:

    R = r1*r2.

For:

    p = a + k*r1
    q = b + l*r2

we have:

    K = k*l
    T = floor(n/R)
    E = T-K.

Ignoring the small floor/carry effects:

    K ~ n/R.

Therefore if:

    R ~ n^alpha,

then:

    K ~ n^(1-alpha).

The tested regimes are:

    alpha = 1/2
        K ~ n^(1/2)

    alpha = 2/3
        K ~ n^(1/3)

    alpha = 3/4
        K ~ n^(1/4)

    alpha = 5/6
        K ~ n^(1/6)

    alpha = 11/12
        K ~ n^(1/12).

For balanced factors and balanced moduli:

    p ~ q ~ n^(1/2)

and:

    r1 ~ r2 ~ n^(alpha/2),

so:

    k,l ~ n^((1-alpha)/2).

At alpha=5/6:

    K ~ n^(1/6)

and:

    k,l ~ n^(1/12).

This is the regime motivating the experiment.

The defect satisfies:

    E = T-K.

The carry decomposition gives a scale of roughly:

    E = O(k+l),

so for balanced cases at alpha=5/6:

    E = O(n^(1/12))

while:

    K = O(n^(1/6)).

Thus:

    E/K -> 0

heuristically as n grows.

This experiment therefore measures three related quantities:

    1. size of K,

    2. size of the observable E/K uncertainty,

    3. frequency of auxiliary Kx=K collisions.

The auxiliary construction is deterministic:

    x = -n mod S.

Hence:

    n+x == 0 mod S.

The auxiliary number is then factored exactly:

    n+x = px*qx

and:

    kx = floor(px/r1)
    lx = floor(qx/r2)

giving:

    Kx = kx*lx.

The event:

    Kx = K

is split into:

    SAME CELL:
        (kx,lx)=(k,l)

and:

    CROSS CELL:
        (kx,lx)!=(k,l)
        but
        kx*lx=K.

The experiment also explicitly constructs the observable
candidate window:

    K_candidate = T-E

for all E within the carry-derived upper bound.

This measures the size of the Level-2 recursive problem.

The central hypothesis is therefore:

    larger R
        ->
    smaller K
        ->
    smaller relative E uncertainty
        ->
    potentially easier recursive factorization.

The experiment does NOT assume that this implication is true.

In particular, a very small K is not sufficient.

The critical question remains whether the auxiliary trajectory
continues to produce useful Kx information when R becomes large.

SymPy is used as an exact factorization oracle for the auxiliary
numbers only.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 100)
    print(
        "HIGH-R MODULUS SCALING / "
        "n^(5/6) K-RECURSION EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"scales                 = "
        f"{[f'{x:.0e}' for x in SCALES]}"
    )

    print(
        f"anchors / scale        = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"R exponents            = "
        f"{R_EXPONENTS}"
    )

    print(
        f"R offsets              = "
        f"{R_OFFSETS}"
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

    start = time.perf_counter()

    all_rows = []

    for scale in SCALES:

        rows = run_scale(
            scale,
            None,
            rng,
        )

        all_rows.extend(rows)

    # The scale runner needs modulus_pairs only for historical API
    # compatibility. The actual pair selector is scale-free.
    #
    # All calculations above are complete.

    summarize(
        all_rows
    )

    summarize_E_scaling(
        all_rows
    )

    cross_scale(
        all_rows
    )

    show_matches(
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
