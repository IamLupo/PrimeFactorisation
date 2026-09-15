#!/usr/bin/env python3

import math
import random
import time
import sympy as sp


# =============================================================================
# CONFIG
# =============================================================================

SEED = 1_511_464_998

SCALES = (
    10**9,
    10**12,
    10**16,
)

ANCHORS_PER_SCALE = 12

K_TARGET = 1000

# Small perturbations around R = n / K_TARGET.
R_OFFSETS = (
    0.95,
    1.00,
    1.05,
)

# Deterministic auxiliary moduli.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Prime modulus search.
PAIR_WINDOW = 5000

# Experimental factor range used to generate anchors.
FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

# Only print a small number of examples.
MAX_EXAMPLES = 4


# =============================================================================
# PRIME UTILITIES
# =============================================================================

def next_prime(x: int) -> int:
    return int(sp.nextprime(max(2, x - 1)))


# =============================================================================
# ANCHORS
# =============================================================================

def generate_anchor(
    scale: int,
    rng: random.Random,
):
    """
    Generate p,q near sqrt(scale).
    """

    root = math.isqrt(scale)

    lo = max(
        FACTOR_MIN,
        int(root * 0.75),
    )

    hi = min(
        FACTOR_MAX,
        int(root * 1.25),
    )

    if lo >= hi:
        raise RuntimeError(
            f"Invalid anchor interval for scale={scale}"
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
# BALANCED R
# =============================================================================

def choose_balanced_pair(
    target_R: int,
):
    """
    Find balanced primes r1,r2 with product close to target_R.
    """

    center = math.isqrt(target_R)

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
                desired - 1,
                desired,
                desired + 1,
            ):

                if probe < 2:
                    continue

                r2 = next_prime(
                    probe
                )

                R = r1 * r2

                # Balancedness.
                lo = min(r1, r2)
                hi = max(r1, r2)

                if (
                    hi / lo
                    > 1.25
                ):
                    continue

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

    if best is None:
        raise RuntimeError(
            f"No balanced prime pair near R={target_R}"
        )

    return (
        best[1],
        best[2],
        best[3],
    )


# =============================================================================
# FACTORIZATION / DIVISOR PAIRS
# =============================================================================

def divisor_pairs(
    value: int,
):
    """
    Exact factorization with SymPy.
    """

    fac = sp.factorint(
        value
    )

    divisors = [1]

    for prime, exp in fac.items():

        old = list(divisors)

        powers = [1]
        cur = 1

        for _ in range(exp):
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

    out = []

    for d in divisors:

        if d > root:
            break

        if value % d == 0:
            out.append(
                (
                    d,
                    value // d,
                )
            )

    return fac, out


# =============================================================================
# CARRY BOUND
# =============================================================================

def observable_E_bound(
    r1: int,
    r2: int,
):
    """
    Build a conservative bound on E using the known experimental
    factor-generation range.

    Since:

        k <= p_max / r1
        l <= q_max / r2

    and:

        E = c1 + c2 + c3
        c1 < l
        c2 < k
        c3 <= 2

    use:

        E <= floor((P_MAX-1)/r1)
           + floor((P_MAX-1)/r2)
           + 2

    This is deliberately conservative.
    """

    k_max = (
        FACTOR_MAX - 1
    ) // r1

    l_max = (
        FACTOR_MAX - 1
    ) // r2

    return k_max + l_max + 2


# =============================================================================
# OBSERVABLE K SEARCH
# =============================================================================

def observable_K_search(
    n: int,
    r1: int,
    r2: int,
):
    """
    We know:

        T = floor(n/R)

    but not E.

    Generate all K=T-E values allowed by the conservative
    experimental E bound and factor each K.

    The TRUE k,l are not used here.
    """

    R = r1 * r2

    T = n // R

    E_bound = min(
        observable_E_bound(
            r1,
            r2,
        ),
        T,
    )

    candidates = []

    for E in range(
        E_bound + 1
    ):

        K = T - E

        if K <= 0:
            break

        candidates.append(
            (
                E,
                K,
            )
        )

    return T, E_bound, candidates


# =============================================================================
# AUXILIARY n+x
# =============================================================================

def auxiliary_hits(
    n: int,
    K_true: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
):
    """
    Construct a few deterministic n+x values and factor them.

    Measure whether their Kx values recover K_true.
    """

    k_true = p // r1
    l_true = q // r2

    results = []

    for S in S_VALUES:

        x = (-n) % S

        if x == 0:
            continue

        nx = n + x

        fac, pairs = divisor_pairs(
            nx
        )

        local_hits = []

        for px, qx in pairs:

            for a, b in (
                (px, qx),
                (qx, px),
            ):

                kx = a // r1
                lx = b // r2

                Kx = kx * lx

                if Kx == K_true:

                    local_hits.append(
                        (
                            a,
                            b,
                            kx,
                            lx,
                        )
                    )

        results.append(
            {
                "S": S,
                "x": x,
                "nx": nx,
                "factor_count": len(
                    pairs
                ),
                "hits": local_hits,
                "factorization": fac,
            }
        )

    return results


# =============================================================================
# ONE CASE
# =============================================================================

def run_case(
    p: int,
    q: int,
    n: int,
    target_R: int,
):
    """
    Full observable K experiment for one modulus scale.
    """

    r1, r2, R = choose_balanced_pair(
        target_R
    )

    k_true = p // r1
    l_true = q // r2

    K_true = k_true * l_true

    T, E_bound, candidates = (
        observable_K_search(
            n,
            r1,
            r2,
        )
    )

    E_true = T - K_true

    # Factor every candidate K.
    candidate_records = []

    true_K_found = False
    true_pair_found = False

    for E, K in candidates:

        fac, pairs = divisor_pairs(
            K
        )

        contains_true = (
            any(
                (
                    a == k_true
                    and b == l_true
                )
                or (
                    a == l_true
                    and b == k_true
                )
                for a, b in pairs
            )
        )

        if K == K_true:
            true_K_found = True

        if contains_true:
            true_pair_found = True

        candidate_records.append(
            {
                "E": E,
                "K": K,
                "factors": fac,
                "pairs": pairs,
                "contains_true": contains_true,
            }
        )

    aux = auxiliary_hits(
        n,
        K_true,
        p,
        q,
        r1,
        r2,
    )

    aux_hits = sum(
        len(a["hits"])
        for a in aux
    )

    return {
        "p": p,
        "q": q,
        "n": n,
        "target_R": target_R,

        "r1": r1,
        "r2": r2,
        "R": R,

        "k": k_true,
        "l": l_true,
        "K": K_true,

        "T": T,
        "E": E_true,
        "E_bound": E_bound,

        "candidate_count": len(
            candidates
        ),

        "true_K_found": true_K_found,
        "true_pair_found": true_pair_found,

        "candidate_records": candidate_records,

        "aux": aux,
        "aux_hits": aux_hits,
    }


# =============================================================================
# MAIN SCALE LOOP
# =============================================================================

def run_scale(
    scale: int,
    rng: random.Random,
):
    """

    Run the fixed-K experiment over several anchors and R offsets.
    """

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

        p, q, n = generate_anchor(
            scale,
            rng,
        )

        print(
            f"anchor {i + 1:2d}/"
            f"{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        for offset in R_OFFSETS:

            target_R = int(
                round(
                    n
                    / K_TARGET
                    * offset
                )
            )

            result = run_case(
                p,
                q,
                n,
                target_R,
            )

            result[
                "offset"
            ] = offset

            rows.append(
                result
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
        "FIXED-K / OBSERVABLE-E SUMMARY"
    )
    print("-" * 100)

    print(
        "scale       cases   mean K      "
        "mean E   mean Ebound   Kstates   "
        "trueK    true(k,l)"
    )

    print("-" * 100)

    scales = sorted(
        set(
            r["n"]
            for r in rows
        )
    )

    # Group by the requested scale indirectly
    groups = {}

    for row in rows:

        scale_key = min(
            SCALES,
            key=lambda s:
                abs(
                    math.log10(
                        s
                    )
                    -
                    math.log10(
                        row["n"]
                    )
                )
        )

        groups.setdefault(
            scale_key,
            [],
        ).append(
            row
        )

    for scale in SCALES:

        group = groups.get(
            scale,
            [],
        )

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

        mean_bound = (
            sum(
                r["E_bound"]
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

        trueK = sum(
            r["true_K_found"]
            for r in group
        )

        truepair = sum(
            r["true_pair_found"]
            for r in group
        )

        print(
            f"{scale:>5.0e} "
            f"{len(group):8d} "
            f"{mean_K:11.3f} "
            f"{mean_E:10.3f} "
            f"{mean_bound:13.3f} "
            f"{mean_states:10.3f} "
            f"{trueK:7d}/{len(group):<3d} "
            f"{truepair:7d}/{len(group):<3d}"
        )


# =============================================================================
# AUXILIARY SUMMARY
# =============================================================================

def print_aux_summary(
    rows,
):

    print()
    print("=" * 100)
    print(
        "AUXILIARY n+x -> Kx TEST"
    )
    print("=" * 100)

    print(
        "scale       cases   aux K-matches   "
        "cases with hit   exact/orientation"
    )

    print("-" * 100)

    for scale in SCALES:

        group = [
            r
            for r in rows
            if min(
                SCALES,
                key=lambda s:
                    abs(
                        math.log10(s)
                        -
                        math.log10(
                            r["n"]
                        )
                    )
            ) == scale
        ]

        hits = sum(
            r["aux_hits"]
            for r in group
        )

        cases = sum(
            r["aux_hits"] > 0
            for r in group
        )

        print(
            f"{scale:>5.0e} "
            f"{len(group):10d} "
            f"{hits:15d} "
            f"{cases:16d}"
        )


# =============================================================================
# REPRESENTATIVE CASES
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

        print()
        print(
            f"n={row['n']:,}"
        )

        print(
            f"    R={row['R']:,} "
            f"(r1,r2)=({row['r1']},{row['r2']})"
        )

        print(
            f"    true (k,l)=("
            f"{row['k']},{row['l']}) "
            f"K={row['K']:,} "
            f"T={row['T']:,} "
            f"E={row['E']:,}"
        )

        print(
            f"    observable E bound="
            f"{row['E_bound']:,} "
            f"K candidates="
            f"{row['candidate_count']:,}"
        )

        print(
            f"    true K in candidate set="
            f"{row['true_K_found']} "
            f"true pair reconstructed="
            f"{row['true_pair_found']}"
        )

        print(
            f"    auxiliary Kx=K hits="
            f"{row['aux_hits']}"
        )

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
The experiment fixes the desired recursive object size:

    K_target = 1000.

For:

    R = r1*r2

we have:

    T = floor(n/R)
    E = T-K
    K = k*l.

The construction chooses:

    R ~= n/1000.

Therefore:

    T ~= 1000

and consequently the true K should also be close to 1000.

The important change is that the candidate search is now
OBSERVABLE.

The experiment knows only:

    n
    r1
    r2
    T

and constructs:

    K_candidate = T-E

over an independently derived conservative E range.

For each candidate K:

    factor(K)
        ->
    divisor pairs
        ->
    candidate (k,l).

The true k,l are used only for validation.

The experiment therefore distinguishes:

    1. TRUE K IDENTIFICATION

       Is the correct K present in the observable K candidate set?

    2. TRUE (k,l) IDENTIFICATION

       Does factoring one of those K values produce the correct
       quotient coordinates?

The auxiliary recursion is independently tested through:

    x = (-n) mod S

so that:

    n+x == 0 mod S.

Then:

    n+x = px*qx

and:

    kx = floor(px/r1)
    lx = floor(qx/r2)
    Kx = kx*lx.

The stronger recursive event is:

    Kx = K.

The central fixed-K hypothesis is:

    R ~= n/K_target

rather than:

    R = n^alpha.

Thus the intended scaling is:

    n=10^9
        R ~= 10^6
        K ~= 10^3

    n=10^12
        R ~= 10^9
        K ~= 10^3

    n=10^16
        R ~= 10^13
        K ~= 10^3.

The experiment specifically asks whether the recursive object K
can remain approximately constant while n increases by seven
orders of magnitude.

A positive structural result would be:

    K remains small
    +
    observable E range remains small
    +
    factor(K) yields few candidate (k,l)
    +
    auxiliary Kx=K events remain common.

This still would not by itself constitute a factorization
algorithm, but it would establish the computational bridge that
the recursive construction requires.

SymPy is used only as an exact experimental factorization oracle.
"""
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    print("=" * 100)
    print(
        "FIXED-K RECURSIVE FACTORIZATION EXPERIMENT"
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

    print_aux_summary(
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
