#!/usr/bin/env python3

"""
====================================================================================================
START EXPERIMENT 73
STATIC K -> (k,l) LOOKUP / RECURSIVE RECONSTRUCTION EXPERIMENT
====================================================================================================

Purpose
-------
Build a static lookup table:

    K -> all divisor pairs (k,l) with k*l = K

and use it in the recursive reconstruction pipeline:

    n+x
      -> factor(n+x)
      -> (kx,lx)
      -> Kx
      -> STATIC_LOOKUP[Kx]
      -> candidate (k,l)
      -> quotient-cell reconstruction

The experiment compares:

    MODE A: dynamically factor Kx
    MODE B: static lookup[Kx]

The mathematical path is identical after obtaining the candidate
divisor pairs. The lookup table only removes repeated factorization
of small K values.

The experiment does NOT use the true K to select a candidate.
The true factorization is used only for validation.

Output is intentionally compact.
====================================================================================================
"""

from __future__ import annotations

import math
import random
import time
from collections import defaultdict

from sympy import factorint, nextprime


# ================================================================================================
# CONFIGURATION
# ================================================================================================

SEED = 1511464998

SCALES = {
    "1e+09": 10**9,
    "1e+12": 10**12,
    "1e+16": 10**16,
}

ANCHORS_PER_SCALE = 3

K_TARGET = 1000
R_OFFSETS = (0.95, 1.00, 1.05)

S_VALUES = (30, 210, 2310, 30030)

# Static table range.
#
# We deliberately make this comfortably larger than the observed
# K range around the K=1000 construction.
K_TABLE_MAX = 5000

# Number of auxiliary factorizations retained per S.
MAX_AUX_FACTOR_PAIRS = 250

# Quotient-cell reconstruction limits.
MAX_GEOMETRY_SURVIVORS = 5000

# Progress output.
PROGRESS_EVERY = 1


# ================================================================================================
# FORMATTING
# ================================================================================================

def fmt(n: int) -> str:
    return f"{n:,}"


def fmt_s(x: float) -> str:
    return f"{x:.4f}s"


# ================================================================================================
# PRIME / ANCHOR GENERATION
# ================================================================================================

def generate_balanced_prime_pair(
    target_r: float,
    rng: random.Random,
    search_radius: int = 5000,
) -> tuple[int, int]:
    """
    Find two reasonably balanced primes whose product is close to target_r.

    This deliberately does not require exact equality.
    """
    if target_r <= 0:
        raise ValueError("target_r must be positive")

    root = max(2, int(math.sqrt(target_r)))

    best = None

    # Search around sqrt(target_r).
    for d1 in range(-search_radius, search_radius + 1):
        p = root + d1
        if p < 2:
            continue

        p = int(nextprime(p - 1))

        q_target = max(2, int(round(target_r / p)))
        q = int(nextprime(q_target - 1))

        R = p * q
        error = abs(R - target_r)

        if best is None or error < best[0]:
            best = (error, p, q)

            if error == 0:
                break

    if best is None:
        raise RuntimeError(f"Could not construct prime pair near R={target_r}")

    _, p, q = best
    return p, q


def generate_anchor(
    scale: int,
    rng: random.Random,
) -> tuple[int, int, int]:
    """
    Generate a balanced semiprime n=p*q near the requested scale.
    """
    root = int(math.isqrt(scale))

    # Keep factors roughly balanced.
    low = max(10_000, int(root * 0.72))
    high = max(low + 1000, int(root * 1.28))

    p = int(nextprime(rng.randint(low, high)))
    q = int(nextprime(rng.randint(low, high)))

    # Avoid trivial equality for the usual semiprime case.
    if p == q:
        q = int(nextprime(q))

    n = p * q

    return p, q, n


# ================================================================================================
# STATIC K LOOKUP TABLE
# ================================================================================================

def build_static_k_lookup(max_k: int) -> dict[int, tuple[tuple[int, int], ...]]:
    """
    Build:

        K -> ((k1,l1), (k2,l2), ...)

    including both orientations.

    Example:

        12 -> ((1,12), (2,6), (3,4),
               (4,3), (6,2), (12,1))
    """
    table: dict[int, list[tuple[int, int]]] = defaultdict(list)

    for k in range(1, max_k + 1):
        limit = math.isqrt(k)
        for d in range(1, limit + 1):
            if k % d == 0:
                q = k // d

                table[k].append((d, q))

                if d != q:
                    table[k].append((q, d))

        table[k].sort()

    return {
        K: tuple(pairs)
        for K, pairs in table.items()
    }


# ================================================================================================
# AUXILIARY CONSTRUCTION
# ================================================================================================

def deterministic_x(n: int, S: int) -> int:
    """
    Smallest nonnegative x such that:

        n+x == 0 mod S
    """
    return (-n) % S


def factor_auxiliary(n_plus_x: int) -> list[tuple[int, int]]:
    """
    Exact factorization oracle using SymPy.

    Returns all divisor pairs (p,q) with p <= q.
    """
    fac = factorint(n_plus_x)

    divisors = [1]
    for prime, exponent in fac.items():
        current = list(divisors)
        power = 1

        for _ in range(exponent):
            power *= prime
            for d in current:
                divisors.append(d * power)

    pairs = []

    for d in divisors:
        if d > n_plus_x // d:
            continue

        q = n_plus_x // d

        # We are interested in nontrivial factor pairs.
        if d > 1 and q > 1:
            pairs.append((d, q))

    pairs.sort()

    return pairs


# ================================================================================================
# QUOTIENT GEOMETRY
# ================================================================================================

def quotient_coordinates(
    p: int,
    q: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int, int]:
    """
    Returns:

        a = p mod r1
        b = q mod r2
        k = floor(p/r1)
        l = floor(q/r2)
        K = k*l
    """
    k = p // r1
    l = q // r2
    a = p % r1
    b = q % r2
    K = k * l

    return a, b, k, l, K


def quotient_cell_contains_factor(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
) -> bool:
    """
    Exact cell test.

    p must be in:

        [k*r1, (k+1)*r1 - 1]

    q must be in:

        [l*r2, (l+1)*r2 - 1]

    and p*q must equal n.
    """
    if k < 0 or l < 0:
        return False

    lo_p = k * r1
    hi_p = (k + 1) * r1 - 1

    lo_q = l * r2
    hi_q = (l + 1) * r2 - 1

    return (
        lo_p <= p <= hi_p
        and lo_q <= q <= hi_q
        and p * q == n
    )


def reconstruct_from_k(
    n: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
) -> tuple[int, int] | None:
    """
    Direct exact reconstruction.

    We need:

        p = k*r1 + a
        q = l*r2 + b

    with:

        0 <= a < r1
        0 <= b < r2
        p*q = n.

    Rather than enumerate all a,b, derive q from p:

        q = n/p.

    The candidate p interval is:

        k*r1 <= p < (k+1)*r1.

    We enumerate the SHORTER quotient-cell interval.

    This is the remaining reconstruction cost of this experiment.
    """
    p_lo = k * r1
    p_hi = (k + 1) * r1 - 1

    q_lo = l * r2
    q_hi = (l + 1) * r2 - 1

    if p_lo < 2:
        p_lo = 2

    if p_hi < p_lo:
        return None

    # Use the shorter side where possible.
    p_span = p_hi - p_lo + 1
    q_span = q_hi - q_lo + 1

    if p_span <= q_span:
        span = min(p_span, MAX_GEOMETRY_SURVIVORS + 1)

        # For very large cells do not perform an accidental huge scan.
        if p_span > MAX_GEOMETRY_SURVIVORS:
            return None

        for p in range(p_lo, p_hi + 1):
            if p == 0 or n % p:
                continue

            q = n // p

            if q_lo <= q <= q_hi and p * q == n:
                return p, q

        return None

    else:
        if q_span > MAX_GEOMETRY_SURVIVORS:
            return None

        for q in range(q_lo, q_hi + 1):
            if q == 0 or n % q:
                continue

            p = n // q

            if p_lo <= p <= p_hi and p * q == n:
                return p, q

        return None


# ================================================================================================
# ONE AUXILIARY FACTORIZATION
# ================================================================================================

def collect_auxiliary_kx(
    n: int,
    r1: int,
    r2: int,
    S: int,
) -> tuple[int, list[tuple[int, int, int, int]]]:
    """
    Returns:

        x
        records = [(px,qx,kx,lx,Kx), ...]

    The full auxiliary n+x is factored exactly.
    """
    x = deterministic_x(n, S)
    nx = n + x

    factor_pairs = factor_auxiliary(nx)

    if len(factor_pairs) > MAX_AUX_FACTOR_PAIRS:
        factor_pairs = factor_pairs[:MAX_AUX_FACTOR_PAIRS]

    records = []

    for px, qx in factor_pairs:
        kx = px // r1
        lx = qx // r2
        Kx = kx * lx

        records.append((px, qx, kx, lx, Kx))

    return x, records


# ================================================================================================
# DYNAMIC MODE
# ================================================================================================

def dynamic_pairs_for_k(K: int) -> tuple[tuple[int, int], ...]:
    """
    Old-style runtime factorization of K using SymPy.
    """
    if K <= 0:
        return ()

    fac = factorint(K)

    divisors = [1]

    for prime, exponent in fac.items():
        current = list(divisors)
        power = 1

        for _ in range(exponent):
            power *= prime
            for d in current:
                divisors.append(d * power)

    pairs = []

    for d in divisors:
        q = K // d
        pairs.append((d, q))

    pairs = sorted(set(pairs))

    return tuple(pairs)


# ================================================================================================
# RECONSTRUCTION FROM K
# ================================================================================================

def reconstruct_candidates(
    n: int,
    r1: int,
    r2: int,
    candidate_K: int,
    pair_provider,
) -> tuple[int, bool]:
    """
    Try all (k,l) pairs for K and apply quotient-cell reconstruction.

    Returns:

        number_of_pairs_tested
        solved
    """
    pairs = pair_provider(candidate_K)

    if len(pairs) == 0:
        return 0, False

    tested = 0

    for k, l in pairs:
        tested += 1

        result = reconstruct_from_k(
            n=n,
            r1=r1,
            r2=r2,
            k=k,
            l=l,
        )

        if result is not None:
            return tested, True

    return tested, False


# ================================================================================================
# ONE CASE
# ================================================================================================

def run_case(
    n: int,
    p_true: int,
    q_true: int,
    r1: int,
    r2: int,
    lookup: dict[int, tuple[tuple[int, int], ...]],
) -> dict:
    R = r1 * r2

    _, _, k_true, l_true, K_true = quotient_coordinates(
        p_true,
        q_true,
        r1,
        r2,
    )

    T = n // R
    E = T - K_true

    # ------------------------------------------------------------
    # Auxiliary stage
    # ------------------------------------------------------------

    aux_seen = set()

    auxiliary_records = []

    for S in S_VALUES:
        x, records = collect_auxiliary_kx(
            n=n,
            r1=r1,
            r2=r2,
            S=S,
        )

        for px, qx, kx, lx, Kx in records:
            if Kx <= 0:
                continue

            aux_seen.add(Kx)

            auxiliary_records.append(
                {
                    "S": S,
                    "x": x,
                    "px": px,
                    "qx": qx,
                    "kx": kx,
                    "lx": lx,
                    "Kx": Kx,
                }
            )

    # ------------------------------------------------------------
    # Compare dynamic factorization versus static lookup.
    # ------------------------------------------------------------

    dynamic_time = 0.0
    static_time = 0.0

    dynamic_tests = 0
    static_tests = 0

    dynamic_solved = False
    static_solved = False

    dynamic_candidates = 0
    static_candidates = 0

    dynamic_seen = set()
    static_seen = set()

    # Important:
    # The order is deterministic and based only on observed Kx values.
    ordered_Kx = sorted(aux_seen)

    for Kx in ordered_Kx:
        # --------------------------------------------
        # Dynamic
        # --------------------------------------------
        t0 = time.perf_counter()

        pairs = dynamic_pairs_for_k(Kx)

        dynamic_time += time.perf_counter() - t0

        dynamic_candidates += len(pairs)

        if not dynamic_solved:
            t1 = time.perf_counter()

            tested, solved = reconstruct_candidates(
                n=n,
                r1=r1,
                r2=r2,
                candidate_K=Kx,
                pair_provider=dynamic_pairs_for_k,
            )

            dynamic_time += time.perf_counter() - t1

            dynamic_tests += tested

            if solved:
                dynamic_solved = True

        dynamic_seen.add(Kx)

        # --------------------------------------------
        # Static lookup
        # --------------------------------------------
        t2 = time.perf_counter()

        pairs_static = lookup.get(Kx, ())

        static_time += time.perf_counter() - t2

        static_candidates += len(pairs_static)

        if not static_solved:
            t3 = time.perf_counter()

            tested, solved = reconstruct_candidates(
                n=n,
                r1=r1,
                r2=r2,
                candidate_K=Kx,
                pair_provider=lambda value, table=lookup: table.get(value, ()),
            )

            static_time += time.perf_counter() - t3

            static_tests += tested

            if solved:
                static_solved = True

        static_seen.add(Kx)

    # ------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------

    true_kx_exposed = K_true in aux_seen
    true_lookup_exists = K_true in lookup

    true_lookup_pairs = lookup.get(K_true, ())

    true_pair_present = (k_true, l_true) in true_lookup_pairs

    return {
        "n": n,
        "p_true": p_true,
        "q_true": q_true,
        "r1": r1,
        "r2": r2,
        "R": R,
        "k_true": k_true,
        "l_true": l_true,
        "K_true": K_true,
        "T": T,
        "E": E,
        "distinct_Kx": len(aux_seen),
        "aux_records": len(auxiliary_records),
        "true_Kx_exposed": true_kx_exposed,
        "true_lookup_exists": true_lookup_exists,
        "true_pair_present": true_pair_present,
        "dynamic_candidates": dynamic_candidates,
        "static_candidates": static_candidates,
        "dynamic_tests": dynamic_tests,
        "static_tests": static_tests,
        "dynamic_solved": dynamic_solved,
        "static_solved": static_solved,
        "dynamic_time": dynamic_time,
        "static_time": static_time,
    }


# ================================================================================================
# SCALE
# ================================================================================================

def run_scale(
    label: str,
    scale: int,
    rng: random.Random,
    lookup: dict[int, tuple[tuple[int, int], ...]],
) -> list[dict]:
    print()
    print("=" * 100)
    print(f"SCALE {label}")
    print("=" * 100)

    rows = []

    for anchor_idx in range(1, ANCHORS_PER_SCALE + 1):
        p_true, q_true, n = generate_anchor(scale, rng)

        print(
            f"anchor {anchor_idx}/{ANCHORS_PER_SCALE} "
            f"n={fmt(n)}"
        )

        # Deduplicate R values.
        R_values = []

        for offset in R_OFFSETS:
            target_R = n / K_TARGET
            effective_target = target_R * offset

            r1, r2 = generate_balanced_prime_pair(
                effective_target,
                rng,
            )

            R = r1 * r2

            candidate = (R, r1, r2)

            if candidate not in R_values:
                R_values.append(candidate)

        for R, r1, r2 in R_values:
            row = run_case(
                n=n,
                p_true=p_true,
                q_true=q_true,
                r1=r1,
                r2=r2,
                lookup=lookup,
            )

            row["scale"] = label

            rows.append(row)

            status = "YES" if row["static_solved"] else "NO"

            print(
                f"    R={fmt(R):>15} "
                f"(r1,r2)=({fmt(r1)},{fmt(r2)}) "
                f"K={fmt(row['K_true']):>6} "
                f"Kx={row['distinct_Kx']:>4} "
                f"K-exposed={'YES' if row['true_Kx_exposed'] else 'NO ':>3} "
                f"static-solved={status} "
                f"static={fmt_s(row['static_time'])}"
            )

    return rows


# ================================================================================================
# MAIN
# ================================================================================================

def main() -> None:
    total_start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("START EXPERIMENT 73")
    print("STATIC K -> (k,l) LOOKUP / RECURSIVE RECONSTRUCTION EXPERIMENT")
    print("=" * 100)

    print()
    print("configuration")
    print(f"    scales                  = {list(SCALES.keys())}")
    print(f"    anchors / scale        = {ANCHORS_PER_SCALE}")
    print(f"    K target               = {K_TARGET}")
    print(f"    R offsets              = {R_OFFSETS}")
    print(f"    S values               = {S_VALUES}")
    print(f"    static K table max     = {K_TABLE_MAX}")
    print(f"    seed                   = {SEED}")

    # ============================================================================================
    # BUILD STATIC TABLE
    # ============================================================================================

    print()
    print("=" * 100)
    print("BUILDING STATIC K LOOKUP TABLE")
    print("=" * 100)

    t_lookup = time.perf_counter()

    lookup = build_static_k_lookup(K_TABLE_MAX)

    lookup_time = time.perf_counter() - t_lookup

    total_pairs = sum(len(v) for v in lookup.values())

    print(f"    K entries              = {len(lookup):,}")
    print(f"    total divisor pairs    = {total_pairs:,}")
    print(f"    build time              = {fmt_s(lookup_time)}")

    # Some validation of the table itself.
    for K in (1, 48, 864, 900, 1000, 1080, 1089):
        pairs = lookup.get(K, ())
        if K <= K_TABLE_MAX:
            assert pairs
            assert all(a * b == K for a, b in pairs)

    # ============================================================================================
    # RUN SCALES
    # ============================================================================================

    all_rows = []

    for label, scale in SCALES.items():
        rows = run_scale(
            label=label,
            scale=scale,
            rng=rng,
            lookup=lookup,
        )

        all_rows.extend(rows)

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"{'scale':<8} "
        f"{'cases':>6} "
        f"{'K-exposed':>10} "
        f"{'dynamic solved':>15} "
        f"{'static solved':>14} "
        f"{'mean Kx':>10} "
        f"{'dyn tests':>12} "
        f"{'stat tests':>12}"
    )

    print("-" * 100)

    for label in SCALES:
        rows = [r for r in all_rows if r["scale"] == label]

        cases = len(rows)
        exposed = sum(r["true_Kx_exposed"] for r in rows)
        dynamic_solved = sum(r["dynamic_solved"] for r in rows)
        static_solved = sum(r["static_solved"] for r in rows)

        mean_Kx = (
            sum(r["distinct_Kx"] for r in rows) / cases
            if cases else 0.0
        )

        dyn_tests = sum(r["dynamic_tests"] for r in rows)
        stat_tests = sum(r["static_tests"] for r in rows)

        print(
            f"{label:<8} "
            f"{cases:>6} "
            f"{exposed:>10} "
            f"{dynamic_solved:>15} "
            f"{static_solved:>14} "
            f"{mean_Kx:>10.2f} "
            f"{dyn_tests:>12,} "
            f"{stat_tests:>12,}"
        )

    # ============================================================================================
    # LOOKUP VS DYNAMIC
    # ============================================================================================

    print()
    print("=" * 100)
    print("STATIC LOOKUP VS DYNAMIC K FACTORIZATION")
    print("=" * 100)

    total_dynamic_time = sum(r["dynamic_time"] for r in all_rows)
    total_static_time = sum(r["static_time"] for r in all_rows)

    total_dynamic_candidates = sum(
        r["dynamic_candidates"]
        for r in all_rows
    )

    total_static_candidates = sum(
        r["static_candidates"]
        for r in all_rows
    )

    print(f"    dynamic K factorization time = {fmt_s(total_dynamic_time)}")
    print(f"    static lookup time           = {fmt_s(total_static_time)}")
    print(
        f"    dynamic divisor pairs seen  = "
        f"{total_dynamic_candidates:,}"
    )
    print(
        f"    static divisor pairs seen   = "
        f"{total_static_candidates:,}"
    )

    if total_static_time > 0:
        print(
            f"    lookup/dynamic time ratio    = "
            f"{total_static_time / total_dynamic_time:.6f}"
        )

    # ============================================================================================
    # K COVERAGE
    # ============================================================================================

    print()
    print("=" * 100)
    print("K COVERAGE")
    print("=" * 100)

    outside_table = [
        r for r in all_rows
        if r["K_true"] > K_TABLE_MAX
    ]

    exposed = sum(r["true_Kx_exposed"] for r in all_rows)
    valid_pairs = sum(r["true_pair_present"] for r in all_rows)

    print(f"    total cases               = {len(all_rows)}")
    print(f"    true K exposed by Kx     = {exposed}/{len(all_rows)}")
    print(f"    true pair in static table= {valid_pairs}/{len(all_rows)}")
    print(
        f"    cases outside table      = "
        f"{len(outside_table)}"
    )

    # ============================================================================================
    # REPRESENTATIVE MATCHES
    # ============================================================================================

    print()
    print("=" * 100)
    print("REPRESENTATIVE K MATCHES")
    print("=" * 100)

    shown = 0

    for row in all_rows:
        if not row["true_Kx_exposed"]:
            continue

        print(
            f"scale={row['scale']} "
            f"n={fmt(row['n'])} "
            f"R={fmt(row['R'])} "
            f"K={fmt(row['K_true'])} "
            f"(k,l)=({row['k_true']},{row['l_true']}) "
            f"Kx-distinct={row['distinct_Kx']} "
            f"static-solved={row['static_solved']}"
        )

        pairs = lookup.get(row["K_true"], ())

        print(
            f"    static lookup pairs={len(pairs)} "
            f"contains true pair={row['true_pair_present']}"
        )

        # Print only a few pairs.
        preview = pairs[:10]
        if len(pairs) > 10:
            suffix = " ..."
        else:
            suffix = ""

        print(f"    pairs={preview}{suffix}")

        shown += 1

        if shown >= 8:
            break

    # ============================================================================================
    # CHECKSUM / CONSISTENCY
    # ============================================================================================

    print()
    print("=" * 100)
    print("CONSISTENCY CHECK")
    print("=" * 100)

    lookup_failures = 0
    identity_failures = 0

    for K, pairs in lookup.items():
        seen = set()

        for k, l in pairs:
            if k * l != K:
                lookup_failures += 1

            if (k, l) in seen:
                lookup_failures += 1

            seen.add((k, l))

    for row in all_rows:
        if row["T"] - row["E"] != row["K_true"]:
            identity_failures += 1

    print(f"    lookup failures          = {lookup_failures}")
    print(f"    K=T-E identity failures  = {identity_failures}")

    # ============================================================================================
    # TIMING
    # ============================================================================================

    total_time = time.perf_counter() - total_start

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"    lookup build             = {fmt_s(lookup_time)}")
    print(f"    dynamic K factorization  = {fmt_s(total_dynamic_time)}")
    print(f"    static lookup            = {fmt_s(total_static_time)}")
    print(f"    total runtime            = {fmt_s(total_time)}")

    print()
    print("=" * 100)
    print("FINISHED EXPERIMENT 73")
    print("=" * 100)


if __name__ == "__main__":
    main()
