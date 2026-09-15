#!/usr/bin/env python3

import bisect
import math
import random
import statistics
import time
from collections import Counter


# ============================================================================
# CONFIGURATION
# ============================================================================

N_MIN = 10_000
N_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

ANCHORS = 300
CLOSE_RATIO = 0.20

SEED = 1_511_464_998


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if s[p]:
            start = p * p
            s[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, v in enumerate(s) if v]


# ============================================================================
# ANCHOR CONSTRUCTION
# ============================================================================

def build_factor_anchors(
    primes: list[int],
    count: int,
    seed: int,
) -> list[tuple[int, int, int]]:
    """
    Returns:
        (n, p, q)

    with p <= q and both inside the factor range.
    """

    rng = random.Random(seed)

    candidates = []

    for p in primes:
        if p < N_MIN:
            continue
        if p > N_MAX:
            break

        for q in primes:
            if q < p:
                continue
            if q > N_MAX:
                break

            n = p * q

            # Keep the experiment in a useful semiprime range.
            candidates.append((n, p, q))

            if len(candidates) > count * 80:
                break

        if len(candidates) > count * 80:
            break

    rng.shuffle(candidates)

    # Prefer a broad spread of n values.
    candidates.sort(key=lambda x: x[0])

    selected = []
    used = set()

    if not candidates:
        raise RuntimeError("No factor anchors found.")

    # Quantile-like selection.
    for i in range(count):
        idx = round(i * (len(candidates) - 1) / max(1, count - 1))
        item = candidates[idx]

        if item[0] not in used:
            selected.append(item)
            used.add(item[0])

    # Fill any holes.
    if len(selected) < count:
        for item in candidates:
            if item[0] not in used:
                selected.append(item)
                used.add(item[0])

            if len(selected) == count:
                break

    return selected[:count]


# ============================================================================
# CLOSE MODULUS TRIPLE
# ============================================================================

def select_close_triple(
    n: int,
    mods: list[int],
    ratio: float,
) -> tuple[int, int, int]:
    """
    Choose r1 < r2 < r3 whose product is <= n and as large as possible.

    This is deliberately deterministic.

    We exploit the sorted modulus list and binary search for r3.
    """

    best = None
    best_product = -1

    for i, r1 in enumerate(mods):
        for j in range(i + 1, len(mods)):
            r2 = mods[j]

            partial = r1 * r2

            if partial >= n:
                break

            limit = n // partial

            pos = bisect.bisect_right(mods, limit) - 1

            if pos <= j:
                continue

            r3 = mods[pos]
            R = partial * r3

            if R > best_product:
                best_product = R
                best = (r1, r2, r3)

    if best is None:
        raise RuntimeError(f"No modulus triple found for n={n}")

    R = math.prod(best)
    ratio_ok = R >= int((1.0 - ratio) * n)

    if not ratio_ok:
        raise RuntimeError(
            f"Selected triple not sufficiently close for n={n}: "
            f"R/n={R/n:.6f}"
        )

    return best


# ============================================================================
# MATHEMATICS
# ============================================================================

def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def t_range_from_rectangle(
    n: int,
    R: int,
    p_min: int,
    p_max: int,
) -> tuple[int, int]:
    """
    Broad rectangle:

        p_min <= p <= p_max
        p_min <= q <= p_max

    Therefore:

        p_min^2 <= pq <= p_max^2

    and:

        pq = n + tR.

    This is deliberately a VERY broad t-domain.
    """

    product_min = p_min * p_min
    product_max = p_max * p_max

    t_min = ceil_div(product_min - n, R)
    t_max = (product_max - n) // R

    return t_min, t_max


def hyperbola_t_bound(
    n: int,
    R: int,
    p_max: int,
) -> tuple[int, int]:
    """
    Along the integer hyperbola quotient:

        q = floor(n/p)

    we have:

        -p < pq-n <= 0

    hence:

        -(p_max) < pq-n <= 0.

    If R > p_max, the only possible integer t is zero.
    """

    t_min = ceil_div(-(p_max - 1), R)
    t_max = 0

    return t_min, t_max


def t_states(t_min: int, t_max: int) -> int:
    if t_max < t_min:
        return 0
    return t_max - t_min + 1


# ============================================================================
# DIRECT HYPERBOLA CONTROL
# ============================================================================

def hyperbola_control(
    n: int,
    R: int,
    p_min: int,
    p_max: int,
) -> tuple[int, int, int]:
    """
    A small control calculation.

    We deliberately do NOT search p values.

    Instead we verify the universal bound:

        -p < pq-n <= 0

    for p <= p_max.

    Returns:

        max_abs_error_bound
        max_possible_t_magnitude
        forced_zero
    """

    max_abs_error = p_max - 1

    possible_t = max_abs_error // R

    forced_zero = int(R > max_abs_error)

    return max_abs_error, possible_t, forced_zero


# ============================================================================
# OPTIONAL CONTROL: SAMPLE INTEGER HYPERBOLA
# ============================================================================

def sampled_hyperbola_t(
    n: int,
    R: int,
    p_min: int,
    p_max: int,
    samples: int,
    rng: random.Random,
) -> Counter:
    """
    Sanity check only.

    Samples p values, computes q=floor(n/p), and records:

        t = floor((p*q-n)/R)

    This should always be zero when R > p_max.

    This does NOT form part of the proposed factoring mechanism.
    """

    counts = Counter()

    if p_max <= p_min:
        return counts

    for _ in range(samples):
        p = rng.randint(p_min, p_max)

        q = n // p
        defect = p * q - n

        # defect is <= 0.
        t = defect // R

        counts[t] += 1

    return counts


# ============================================================================
# RUN
# ============================================================================

def run() -> None:
    start_total = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME t-QUOTIENT / HYPERBOLA EXACTNESS EXPERIMENT")
    print("=" * 100)
    print(f"factor range              = {N_MIN:,} - {N_MAX:,}")
    print(f"modulus prime range       = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"anchors                   = {ANCHORS}")
    print(f"seed                      = {SEED}")
    print()

    # ------------------------------------------------------------------------
    # PRIME POOLS
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    all_primes = sieve(N_MAX)
    factor_primes = [
        p for p in all_primes
        if N_MIN <= p <= N_MAX
    ]

    modulus_primes = [
        p for p in all_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # ------------------------------------------------------------------------
    # ANCHORS
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_factor_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    print(f"actual anchors             = {len(anchors)}")
    print()

    # ------------------------------------------------------------------------
    # MODULI
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("SELECTING CLOSE MODULUS TRIPLES")
    print("=" * 100)

    records = []

    for idx, (n, p, q) in enumerate(anchors, start=1):
        mods = select_close_triple(
            n,
            modulus_primes,
            CLOSE_RATIO,
        )

        r1, r2, r3 = mods
        R = r1 * r2 * r3

        records.append({
            "id": idx,
            "n": n,
            "p": p,
            "q": q,
            "r1": r1,
            "r2": r2,
            "r3": r3,
            "R": R,
        })

        if idx % 25 == 0 or idx == len(anchors):
            print(f"anchor {idx:3d}/{len(anchors)}")

    print()

    # ------------------------------------------------------------------------
    # ANALYSIS
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("RUNNING t-DOMAIN ANALYSIS")
    print("=" * 100)

    rectangle_states = []
    hyperbola_states = []

    forced_zero_count = 0
    R_gt_pmax_count = 0

    ratios = []
    gaps = []

    for rec in records:
        n = rec["n"]
        R = rec["R"]

        t_rect_min, t_rect_max = t_range_from_rectangle(
            n,
            R,
            N_MIN,
            N_MAX,
        )

        t_hyp_min, t_hyp_max = hyperbola_t_bound(
            n,
            R,
            N_MAX,
        )

        rect_count = t_states(
            t_rect_min,
            t_rect_max,
        )

        hyp_count = t_states(
            t_hyp_min,
            t_hyp_max,
        )

        _, possible_t_mag, forced_zero = hyperbola_control(
            n,
            R,
            N_MIN,
            N_MAX,
        )

        rectangle_states.append(rect_count)
        hyperbola_states.append(hyp_count)

        if forced_zero:
            forced_zero_count += 1

        if R > N_MAX - 1:
            R_gt_pmax_count += 1

        ratios.append(R / n)
        gaps.append(n - R)

        rec["t_rect_min"] = t_rect_min
        rec["t_rect_max"] = t_rect_max
        rec["t_rect_count"] = rect_count

        rec["t_hyp_min"] = t_hyp_min
        rec["t_hyp_max"] = t_hyp_max
        rec["t_hyp_count"] = hyp_count

        rec["forced_zero"] = forced_zero

    # ------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"average rectangle t states = "
        f"{statistics.mean(rectangle_states):.3f}"
    )

    print(
        f"minimum rectangle t states = "
        f"{min(rectangle_states)}"
    )

    print(
        f"maximum rectangle t states = "
        f"{max(rectangle_states)}"
    )

    print(
        f"average hyperbola t states = "
        f"{statistics.mean(hyperbola_states):.3f}"
    )

    print(
        f"minimum hyperbola t states = "
        f"{min(hyperbola_states)}"
    )

    print(
        f"maximum hyperbola t states = "
        f"{max(hyperbola_states)}"
    )

    print(
        f"R > p_max-1               = "
        f"{R_gt_pmax_count}/{len(records)}"
    )

    print(
        f"hyperbola forces t=0      = "
        f"{forced_zero_count}/{len(records)}"
    )

    print()

    # ------------------------------------------------------------------------
    # REDUCTION
    # ------------------------------------------------------------------------

    avg_rect = statistics.mean(rectangle_states)
    avg_hyp = statistics.mean(hyperbola_states)

    print("=" * 100)
    print("t-SPACE REDUCTION")
    print("=" * 100)

    print(
        f"hyperbola / rectangle      = "
        f"{avg_hyp / avg_rect:.12f}"
    )

    print(
        f"rectangle -> hyperbola     = "
        f"{100.0 * (1.0 - avg_hyp / avg_rect):.6f}%"
    )

    print()

    # ------------------------------------------------------------------------
    # R/N
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("MODULUS SCALE")
    print("=" * 100)

    print(
        f"mean R/n                   = "
        f"{statistics.mean(ratios):.12f}"
    )

    print(
        f"minimum R/n                = "
        f"{min(ratios):.12f}"
    )

    print(
        f"maximum R/n                = "
        f"{max(ratios):.12f}"
    )

    print(
        f"mean gap n-R               = "
        f"{statistics.mean(gaps):.3f}"
    )

    print(
        f"maximum gap n-R            = "
        f"{max(gaps):,}"
    )

    print()

    # ------------------------------------------------------------------------
    # SAMPLE HYPERBOLA CONTROL
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("HYPERBOLA t CONTROL")
    print("=" * 100)

    rng = random.Random(SEED ^ 0xA5A5A5A5)

    sampled_nonzero = 0
    sampled_total = 0
    sampled_distribution = Counter()

    for rec in records[:20]:
        sample = sampled_hyperbola_t(
            rec["n"],
            rec["R"],
            N_MIN,
            N_MAX,
            samples=500,
            rng=rng,
        )

        sampled_distribution.update(sample)
        sampled_total += sum(sample.values())

        for t, count in sample.items():
            if t != 0:
                sampled_nonzero += count

    print(f"sampled hyperbola states   = {sampled_total:,}")
    print(f"sampled nonzero t states   = {sampled_nonzero:,}")

    print()
    print("sampled t distribution:")
    for t in sorted(sampled_distribution):
        print(
            f"    t = {t:4d} : "
            f"{sampled_distribution[t]:,}"
        )

    print()

    # ------------------------------------------------------------------------
    # EXAMPLES
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for rec in records[:20]:
        print(
            f"n={rec['n']:,} "
            f"p={rec['p']:,} "
            f"q={rec['q']:,} "
            f"R={rec['R']:,} "
            f"R/n={rec['R']/rec['n']:.10f} "
            f"t_rect=[{rec['t_rect_min']},{rec['t_rect_max']}] "
            f"rect={rec['t_rect_count']} "
            f"t_hyp=[{rec['t_hyp_min']},{rec['t_hyp_max']}] "
            f"hyp={rec['t_hyp_count']} "
            f"forced_zero={bool(rec['forced_zero'])}"
        )

    # ------------------------------------------------------------------------
    # MATHEMATICAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)
    print()
    print(
        "The candidate congruence is:"
    )
    print()
    print(
        "    p*q == n (mod R)"
    )
    print()
    print(
        "which is equivalent to:"
    )
    print()
    print(
        "    p*q = n + t*R."
    )
    print()
    print(
        "The first measurement uses only the broad factor box:"
    )
    print()
    print(
        "    p_min <= p,q <= p_max."
    )
    print()
    print(
        "This gives a potentially nontrivial t interval."
    )
    print()
    print(
        "But the actual hyperbola relation is much stronger:"
    )
    print()
    print(
        "    q = floor(n/p)."
    )
    print()
    print(
        "Therefore:"
    )
    print()
    print(
        "    p*q <= n"
    )
    print(
        "and:"
    )
    print()
    print(
        "    n - p*q < p <= p_max."
    )
    print()
    print(
        "Hence:"
    )
    print()
    print(
        "    -(p_max-1) <= p*q-n <= 0."
    )
    print()
    print(
        "If:"
    )
    print()
    print(
        "    R > p_max-1,"
    )
    print()
    print(
        "then the only multiple of R in this entire defect interval"
    )
    print(
        "is zero."
    )
    print()
    print(
        "Therefore:"
    )
    print()
    print(
        "    p*q == n (mod R)"
    )
    print(
        "    + hyperbola quotient"
    )
    print(
        "        =>"
    )
    print(
        "    t = 0"
    )
    print(
        "        =>"
    )
    print(
        "    p*q = n."
    )
    print()
    print(
        "This is an important diagnostic."
    )
    print()
    print(
        "If the experiment reports:"
    )
    print()
    print(
        "    hyperbola t states = 1"
    )
    print()
    print(
        "for essentially every anchor, then t itself cannot provide"
    )
    print(
        "a new search dimension."
    )
    print()
    print(
        "The modular condition has become an exact divisibility"
    )
    print(
        "condition once the hyperbola bound is imposed."
    )
    print()
    print(
        "The remaining problem is therefore:"
    )
    print()
    print(
        "    generate p (or q)"
    )
    print(
        "without enumerating the whole factor interval."
    )
    print()
    print(
        "That means the next useful direction is not another"
    )
    print(
        "reparameterization of t."
    )
    print()
    print(
        "The next serious target is a direct solver for the bilinear"
    )
    print(
        "constraint:"
    )
    print()
    print(
        "    p*q = n"
    )
    print(
        "    p ≡ A (mod R)"
    )
    print()
    print(
        "or an equivalent low-dimensional lattice problem, where the"
    )
    print(
        "candidate p/q values themselves are generated rather than"
    )
    print(
        "enumerated."
    )
    print()
    print(
        "Every eventual candidate must still satisfy:"
    )
    print()
    print(
        "    p*q == n."
    )

    # ------------------------------------------------------------------------
    # TIMING
    # ------------------------------------------------------------------------

    total_time = time.perf_counter() - start_total

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime               = {total_time:.3f} seconds")
    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
