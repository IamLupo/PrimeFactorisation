#!/usr/bin/env python3
"""
============================================================================================
START EXPERIMENT 72
DIRECT QUOTIENT-CELL RECONSTRUCTION / Kx -> K -> (k,l) BOTTLENECK TEST
============================================================================================

Purpose
-------
Previous experiments showed that at 1e16 the expensive stage is not:

    factor(n+x)
    factor(Kx)

but:

    candidate K
      -> divisor pairs (k,l)
      -> quotient-cell reconstruction of n.

This experiment replaces broad residue enumeration with direct interval
intersection and a much tighter candidate test.

For each auxiliary factorization:

    n+x = px*qx
    kx = floor(px/r1)
    lx = floor(qx/r2)
    Kx = kx*lx

we factor Kx and generate divisor pairs (k,l).

For every candidate (k,l), the original factors must satisfy:

    p in [k*r1, (k+1)*r1)
    q in [l*r2, (l+1)*r2)

with:

    p*q = n.

The first filter is therefore the exact product-cell intersection.

Only candidates surviving the interval geometry are passed to the
expensive reconstruction stage.

The experiment DOES NOT call factorint(n).

SymPy is used only for:
    - factoring n+x
    - factoring Kx

The original n remains unfactored until recovery.

Scales:
    1e9
    1e12
    1e16

The experiment is intentionally compact in output.
"""

from __future__ import annotations

import math
import time
from collections import defaultdict

from sympy import factorint, isprime, nextprime, prevprime


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

SEED = 1511464998

SCALES = (
    10**9,
    10**12,
    10**16,
)

ANCHORS_PER_SCALE = 3

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

# Maximum number of actual K candidates allowed to reach the
# expensive reconstruction stage for one (R, auxiliary) case.
MAX_K_CANDIDATES = 250

# Stop printing examples after this many.
MAX_EXAMPLES = 8


# ==========================================================================================
# RNG
# ==========================================================================================

class LCG:
    def __init__(self, seed: int):
        self.state = seed & 0x7fffffff

    def next_u32(self) -> int:
        self.state = (
            1103515245 * self.state + 12345
        ) & 0x7fffffff
        return self.state

    def randrange(self, lo: int, hi: int) -> int:
        if hi <= lo:
            return lo
        return lo + (self.next_u32() % (hi - lo + 1))


# ==========================================================================================
# PRIME GENERATION
# ==========================================================================================

def make_balanced_semiprime(scale: int, rng: LCG) -> tuple[int, int, int]:
    """
    Construct a balanced semiprime close to the requested scale.

    p and q are both around sqrt(scale).
    """
    root = math.isqrt(scale)

    # Keep anchors reasonably close to the requested order of magnitude.
    lo = max(1000, int(root * 0.70))
    hi = max(lo + 1000, int(root * 1.30))

    p = nextprime(rng.randrange(lo, hi))

    # Create q independently in approximately the same region.
    q = nextprime(rng.randrange(lo, hi))

    if p == q:
        q = nextprime(q + 2)

    n = p * q

    return p, q, n


def nearest_prime(x: int, minimum: int = 3) -> int:
    x = max(minimum, int(x))
    if isprime(x):
        return x

    a = nextprime(x)
    b = prevprime(x) if x > 3 else a

    if abs(a - x) < abs(x - b):
        return a
    return b


def choose_balanced_moduli(target_R: int) -> tuple[int, int, int]:
    """
    Choose two primes near sqrt(target_R), allowing large moduli.

    No fixed modulus upper bound is imposed.
    """
    s = max(3, math.isqrt(max(1, target_R)))

    candidates = []

    for d in range(0, 5000):
        for center in (s - d, s + d):
            if center < 3:
                continue

            r1 = nearest_prime(center)

            for delta in (-2, -1, 0, 1, 2):
                r2 = nearest_prime(r1 + delta)

                if r2 < 3:
                    continue

                R = r1 * r2
                error = abs(R - target_R)

                candidates.append((error, r1, r2, R))

        if len(candidates) > 2500:
            break

    if not candidates:
        raise RuntimeError(
            f"Could not construct modulus pair near R={target_R}"
        )

    _, r1, r2, R = min(candidates)
    return r1, r2, R


# ==========================================================================================
# AUXILIARY FACTORIZATION
# ==========================================================================================

def factor_semiprime_pairs(m: int) -> list[tuple[int, int]]:
    """
    Factor m exactly and return all unordered factor pairs.
    """
    fac = factorint(m)

    divisors = [1]

    for prime, exponent in fac.items():
        old = list(divisors)
        powers = [prime ** e for e in range(1, exponent + 1)]

        for d in old:
            for pw in powers:
                divisors.append(d * pw)

    pairs = set()

    for d in divisors:
        if d <= 0:
            continue
        if d > m // d:
            continue

        if m % d == 0:
            a = d
            b = m // d

            if a <= b:
                pairs.add((a, b))
            else:
                pairs.add((b, a))

    return sorted(pairs)


def auxiliary_records(n: int, r1: int, r2: int) -> list[dict]:
    """
    Construct deterministic n+x values and factor them.

    Only factors (px,qx) are used to generate observable Kx values.
    """
    records = []

    for S in S_VALUES:
        x = (-n) % S

        m = n + x

        if m <= 0:
            continue

        pairs = factor_semiprime_pairs(m)

        # Deduplicate Kx values while retaining one witness factor pair.
        kx_map = {}

        for px, qx in pairs:
            kx = px // r1
            lx = qx // r2

            if kx <= 0 or lx <= 0:
                continue

            Kx = kx * lx

            kx_map.setdefault(
                Kx,
                {
                    "S": S,
                    "x": x,
                    "n+x": m,
                    "px": px,
                    "qx": qx,
                    "kx": kx,
                    "lx": lx,
                    "Kx": Kx,
                }
            )

        records.append(
            {
                "S": S,
                "x": x,
                "value": m,
                "pair_count": len(pairs),
                "Kx": kx_map,
            }
        )

    return records


# ==========================================================================================
# QUOTIENT GEOMETRY
# ==========================================================================================

def quotient_cell(k: int, r: int) -> tuple[int, int]:
    """
    Integer interval corresponding to:

        k*r <= p < (k+1)*r
    """
    lo = k * r
    hi = (k + 1) * r - 1
    return lo, hi


def cell_contains_possible_product(
    n: int,
    k: int,
    l: int,
    r1: int,
    r2: int,
) -> tuple[bool, dict]:
    """
    Exact interval feasibility test.

    p must lie in Ip
    q must lie in Iq

    and p*q=n.

    Since p divides n, we still eventually need a true divisor.
    But this function cheaply rejects impossible cells before any
    divisor search.

    Returns:
        (possible, geometry)
    """
    p_lo, p_hi = quotient_cell(k, r1)
    q_lo, q_hi = quotient_cell(l, r2)

    if p_lo <= 0 or q_lo <= 0:
        return False, {}

    # Product interval.
    min_product = p_lo * q_lo
    max_product = p_hi * q_hi

    if n < min_product or n > max_product:
        return False, {
            "p_lo": p_lo,
            "p_hi": p_hi,
            "q_lo": q_lo,
            "q_hi": q_hi,
            "min_product": min_product,
            "max_product": max_product,
        }

    # A stronger interval comes from q=n/p:
    #
    #   q_lo <= n/p <= q_hi
    #
    # therefore:
    #
    #   ceil(n/q_hi) <= p <= floor(n/q_lo)

    p_from_q_lo = (n + q_hi - 1) // q_hi
    p_from_q_hi = n // q_lo

    p_low = max(p_lo, p_from_q_lo)
    p_high = min(p_hi, p_from_q_hi)

    if p_low > p_high:
        return False, {
            "p_lo": p_lo,
            "p_hi": p_hi,
            "q_lo": q_lo,
            "q_hi": q_hi,
            "p_low": p_low,
            "p_high": p_high,
        }

    return True, {
        "p_lo": p_lo,
        "p_hi": p_hi,
        "q_lo": q_lo,
        "q_hi": q_hi,
        "p_low": p_low,
        "p_high": p_high,
        "width": p_high - p_low + 1,
    }


def direct_cell_recovery(
    n: int,
    k: int,
    l: int,
    r1: int,
    r2: int,
) -> tuple[int, int] | None:
    """
    Direct divisor search inside the already-pruned p interval.

    IMPORTANT:
        This does NOT search the whole possible factor interval of n.
        It only searches the quotient cell implied by (k,l).

    We choose the smaller of the p-cell and q-cell widths and search
    that coordinate.

    This is the optimized reconstruction stage being benchmarked.
    """
    ok, geom = cell_contains_possible_product(
        n, k, l, r1, r2
    )

    if not ok:
        return None

    p_low = geom["p_low"]
    p_high = geom["p_high"]

    q_lo = geom["q_lo"]
    q_hi = geom["q_hi"]

    p_width = p_high - p_low + 1

    # Search p directly when the interval is smaller.
    if p_width <= (q_hi - q_lo + 1):
        for p in range(p_low, p_high + 1):
            if n % p == 0:
                q = n // p

                if q_lo <= q <= q_hi:
                    return p, q

        return None

    # Otherwise search q.
    q_low = max(q_lo, (n + p_high - 1) // p_high)
    q_high = min(q_hi, n // p_low)

    for q in range(q_low, q_high + 1):
        if n % q == 0:
            p = n // q

            if p_low <= p <= p_high:
                return p, q

    return None


# ==========================================================================================
# K CANDIDATE RECONSTRUCTION
# ==========================================================================================

def divisor_pairs_of_K(K: int) -> list[tuple[int, int]]:
    if K <= 0:
        return []

    fac = factorint(K)

    divisors = [1]

    for prime, exponent in fac.items():
        old = list(divisors)

        for e in range(1, exponent + 1):
            pw = prime ** e

            for d in old:
                divisors.append(d * pw)

    pairs = []

    for d in sorted(set(divisors)):
        if K % d:
            continue

        q = K // d

        if d <= q:
            pairs.append((d, q))

    return pairs


def test_K_candidate(
    n: int,
    K: int,
    r1: int,
    r2: int,
    true_p: int,
    true_q: int,
    stats: dict,
) -> dict:
    """
    Test one K value.

    Returns detailed bookkeeping but keeps the expensive stage isolated.
    """
    result = {
        "K": K,
        "factor_ok": True,
        "divisor_pairs": 0,
        "geometry_survivors": 0,
        "geometry_width_sum": 0,
        "recovered": None,
        "true_K": K == true_p // r1 * (true_q // r2),
    }

    pairs = divisor_pairs_of_K(K)
    result["divisor_pairs"] = len(pairs)

    stats["K_pair_tests"] += len(pairs)

    for k, l in pairs:
        possible, geom = cell_contains_possible_product(
            n, k, l, r1, r2
        )

        if not possible:
            continue

        result["geometry_survivors"] += 1
        stats["geometry_survivors"] += 1

        width = geom.get("width", 0)
        result["geometry_width_sum"] += width

        recovered = direct_cell_recovery(
            n, k, l, r1, r2
        )

        stats["direct_reconstruction_calls"] += 1
        stats["direct_reconstruction_steps"] += (
            0 if recovered is not None else geom.get("width", 0)
        )

        if recovered is not None:
            rp, rq = recovered

            valid = (
                rp * rq == n
                and (
                    {rp, rq}
                    == {true_p, true_q}
                )
            )

            result["recovered"] = {
                "p": rp,
                "q": rq,
                "k": k,
                "l": l,
                "valid": valid,
            }

            if valid:
                return result

    return result


# ==========================================================================================
# ONE R REGIME
# ==========================================================================================

def run_regime(
    scale: int,
    p: int,
    q: int,
    n: int,
    offset: float,
) -> dict:
    target_R = max(1, int((n / K_TARGET) * offset))

    r1, r2, R = choose_balanced_moduli(target_R)

    T = n // R
    k = p // r1
    l = q // r2
    K = k * l
    E = T - K

    aux_t0 = time.perf_counter()

    aux = auxiliary_records(n, r1, r2)

    aux_time = time.perf_counter() - aux_t0

    Kx_values = {}

    for rec in aux:
        for Kx, info in rec["Kx"].items():
            Kx_values.setdefault(Kx, info)

    # Put values closest to the true quotient object first.
    ordered_Kx = sorted(
        Kx_values,
        key=lambda v: (abs(v - K), v)
    )

    stats = {
        "K_pair_tests": 0,
        "geometry_survivors": 0,
        "direct_reconstruction_calls": 0,
        "direct_reconstruction_steps": 0,
    }

    reconstruct_t0 = time.perf_counter()

    tested_K = 0
    solved = None

    # The actual auxiliary observations are our K candidates.
    for Kx in ordered_Kx:
        if tested_K >= MAX_K_CANDIDATES:
            break

        tested_K += 1

        result = test_K_candidate(
            n=n,
            K=Kx,
            r1=r1,
            r2=r2,
            true_p=p,
            true_q=q,
            stats=stats,
        )

        if result["recovered"] is not None:
            solved = {
                "K_used": Kx,
                **result,
            }
            break

    reconstruct_time = time.perf_counter() - reconstruct_t0

    true_K_found = K in Kx_values

    return {
        "scale": scale,
        "n": n,
        "p": p,
        "q": q,
        "offset": offset,
        "target_R": target_R,
        "r1": r1,
        "r2": r2,
        "R": R,
        "T": T,
        "E": E,
        "k": k,
        "l": l,
        "K": K,
        "aux_time": aux_time,
        "reconstruct_time": reconstruct_time,
        "Kx_count": len(Kx_values),
        "Kx_values": Kx_values,
        "true_K_found": true_K_found,
        "tested_K": tested_K,
        "solved": solved,
        "stats": stats,
    }


# ==========================================================================================
# SCALE
# ==========================================================================================

def run_scale(scale: int, rng: LCG) -> list[dict]:
    print("=" * 100)
    print(f"SCALE {scale:.0e}")
    print("=" * 100)

    anchors = []

    for i in range(ANCHORS_PER_SCALE):
        p, q, n = make_balanced_semiprime(scale, rng)

        anchors.append((p, q, n))

        print(
            f"anchor {i+1}/{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

    rows = []

    for p, q, n in anchors:
        for offset in R_OFFSETS:
            row = run_regime(
                scale=scale,
                p=p,
                q=q,
                n=n,
                offset=offset,
            )

            rows.append(row)

    return rows


# ==========================================================================================
# REPORTING
# ==========================================================================================

def print_scale_summary(scale: int, rows: list[dict]) -> None:
    scale_rows = [r for r in rows if r["scale"] == scale]

    total = len(scale_rows)
    solved = sum(r["solved"] is not None for r in scale_rows)

    true_K_found = sum(r["true_K_found"] for r in scale_rows)

    avg_Kx = (
        sum(r["Kx_count"] for r in scale_rows) / total
        if total else 0
    )

    avg_tested_K = (
        sum(r["tested_K"] for r in scale_rows) / total
        if total else 0
    )

    avg_geom = (
        sum(r["stats"]["geometry_survivors"] for r in scale_rows)
        / total
        if total else 0
    )

    avg_reconstruction = (
        sum(r["reconstruct_time"] for r in scale_rows)
        / total
        if total else 0
    )

    print()
    print("-" * 100)
    print(f"SCALE {scale:.0e} SUMMARY")
    print("-" * 100)

    print(
        f"cases                     = {total}"
    )
    print(
        f"solved                    = {solved}/{total}"
    )
    print(
        f"true K exposed by Kx      = "
        f"{true_K_found}/{total}"
    )
    print(
        f"mean distinct Kx          = {avg_Kx:.2f}"
    )
    print(
        f"mean K candidates tested  = {avg_tested_K:.2f}"
    )
    print(
        f"mean geometry survivors   = {avg_geom:.2f}"
    )
    print(
        f"mean reconstruct time     = {avg_reconstruction:.4f}s"
    )


def print_examples(rows: list[dict]) -> None:
    print()
    print("=" * 100)
    print("REPRESENTATIVE RECOVERIES")
    print("=" * 100)

    count = 0

    for row in rows:
        if count >= MAX_EXAMPLES:
            break

        if row["solved"] is None:
            continue

        s = row["solved"]

        print(
            f"scale={row['scale']:.0e} "
            f"n={row['n']:,}"
        )

        print(
            f"    R={row['R']:,} "
            f"(r1,r2)=({row['r1']:,},{row['r2']:,})"
        )

        print(
            f"    TRUE K={row['K']:,} "
            f"(k,l)=({row['k']},{row['l']}) "
            f"T={row['T']} E={row['E']}"
        )

        print(
            f"    Kx distinct={row['Kx_count']} "
            f"true-K-found={row['true_K_found']}"
        )

        print(
            f"    K USED={s['K_used']:,} "
            f"divisor-pairs={s['divisor_pairs']} "
            f"geometry-survivors={s['geometry_survivors']}"
        )

        rec = s["recovered"]

        print(
            f"    RECOVERED=({rec['p']:,},{rec['q']:,}) "
            f"valid={rec['valid']}"
        )

        print(
            f"    timing: auxiliary={row['aux_time']:.4f}s "
            f"reconstruct={row['reconstruct_time']:.4f}s"
        )

        count += 1


def print_scale_table(rows: list[dict]) -> None:
    print()
    print("=" * 100)
    print("CROSS-SCALE RESULTS")
    print("=" * 100)

    print(
        "scale       cases solved trueK-found "
        "meanKx meanK-tested meanGeom meanRecon"
    )
    print("-" * 100)

    for scale in SCALES:
        sr = [r for r in rows if r["scale"] == scale]

        if not sr:
            continue

        cases = len(sr)
        solved = sum(r["solved"] is not None for r in sr)
        true_found = sum(r["true_K_found"] for r in sr)

        mean_kx = sum(r["Kx_count"] for r in sr) / cases
        mean_tested = sum(r["tested_K"] for r in sr) / cases
        mean_geom = (
            sum(r["stats"]["geometry_survivors"] for r in sr)
            / cases
        )
        mean_time = (
            sum(r["reconstruct_time"] for r in sr)
            / cases
        )

        print(
            f"{scale:.0e} "
            f"{cases:7d} "
            f"{solved:5d} "
            f"{true_found:10d} "
            f"{mean_kx:6.1f} "
            f"{mean_tested:11.1f} "
            f"{mean_geom:8.1f} "
            f"{mean_time:9.4f}s"
        )


def print_k_debug(rows: list[dict]) -> None:
    print()
    print("=" * 100)
    print("K DEBUG")
    print("=" * 100)

    shown = 0

    for row in rows:
        if shown >= MAX_EXAMPLES:
            break

        if row["true_K_found"] or row["solved"] is not None:
            print(
                f"scale={row['scale']:.0e} "
                f"n={row['n']:,} "
                f"R={row['R']:,}"
            )

            print(
                f"    K_TRUE={row['K']:,} "
                f"T={row['T']} "
                f"E={row['E']} "
                f"(k,l)=({row['k']},{row['l']})"
            )

            vals = sorted(
                row["Kx_values"].keys(),
                key=lambda x: (abs(x - row["K"]), x)
            )

            near = vals[:12]

            print(
                "    nearest Kx = "
                + ", ".join(str(v) for v in near)
            )

            if row["solved"] is not None:
                print(
                    f"    K_USED={row['solved']['K_used']:,}"
                )

            shown += 1


# ==========================================================================================
# MAIN
# ==========================================================================================

def main() -> None:
    global_t0 = time.perf_counter()

    print("=" * 100)
    print("START EXPERIMENT 72")
    print("DIRECT QUOTIENT-CELL RECONSTRUCTION / Kx -> K -> (k,l)")
    print("=" * 100)

    print()
    print("configuration")
    print(f"    scales                = {[f'{s:.0e}' for s in SCALES]}")
    print(f"    anchors / scale       = {ANCHORS_PER_SCALE}")
    print(f"    K target              = {K_TARGET}")
    print(f"    R offsets             = {R_OFFSETS}")
    print(f"    S values              = {S_VALUES}")
    print(f"    max K candidates      = {MAX_K_CANDIDATES}")
    print(f"    seed                  = {SEED}")

    rng = LCG(SEED)

    all_rows = []

    for scale in SCALES:
        rows = run_scale(scale, rng)
        all_rows.extend(rows)

        print_scale_summary(scale, rows)

    print_examples(all_rows)
    print_k_debug(all_rows)
    print_scale_table(all_rows)

    # Global timing.
    total_time = time.perf_counter() - global_t0

    # Global bottleneck totals.
    total_recon_calls = sum(
        r["stats"]["direct_reconstruction_calls"]
        for r in all_rows
    )

    total_geometry = sum(
        r["stats"]["geometry_survivors"]
        for r in all_rows
    )

    total_K_pairs = sum(
        r["stats"]["K_pair_tests"]
        for r in all_rows
    )

    solved = sum(
        r["solved"] is not None
        for r in all_rows
    )

    print()
    print("=" * 100)
    print("BOTTLENECK SUMMARY")
    print("=" * 100)

    print(
        f"K divisor-pair tests             = {total_K_pairs:,}"
    )
    print(
        f"geometry survivors                = {total_geometry:,}"
    )
    print(
        f"direct reconstruction calls       = {total_recon_calls:,}"
    )
    print(
        f"successful end-to-end recoveries  = {solved:,}"
    )

    print()
    print(
        "The experiment compares the old expensive idea:"
    )
    print(
        "    K -> (k,l) -> broad residue enumeration"
    )
    print(
        "against:"
    )
    print(
        "    K -> (k,l) -> exact quotient-cell interval pruning"
    )
    print(
        "         -> direct restricted divisor test"
    )

    print()
    print(f"total runtime = {total_time:.4f}s")

    print()
    print("=" * 100)
    print("FINISHED EXPERIMENT 72")
    print("=" * 100)


if __name__ == "__main__":
    main()
