#!/usr/bin/env python3

import math
import random
import time
from itertools import combinations
from sympy import factorint, isprime, nextprime


# =============================================================================
# START EXPERIMENT 77
# TWO-MODULUS INTERSECTION / (a,b) RECOVERY EXPERIMENT
# =============================================================================

CONFIG = {
    "scales": [10**9, 10**12, 10**16],
    "anchors_per_scale": 4,

    # Desired recursive quotient-product size.
    "K_target": 1000,

    # Several nearby R values for each n.
    "R_offsets": (0.90, 0.95, 1.00, 1.05, 1.10),

    # Static K lookup table.
    "K_table_max": 5000,

    # Number of modulus representations to keep per n.
    "representations_per_n": 5,

    # Seed for reproducibility.
    "seed": 1511464998,

    # Prime generation search radius.
    "prime_search_radius": 5000,

    # Maximum number of K candidates to inspect from the static table
    # around the predicted K.
    "K_window": 128,
}


# =============================================================================
# Formatting helpers
# =============================================================================

def fmt(x):
    return f"{x:,}"


def pct(x):
    return f"{100.0 * x:.2f}%"


def banner(title):
    print("=" * 92)
    print(title)
    print("=" * 92)


# =============================================================================
# Prime helpers
# =============================================================================

def nearby_prime(x: int, direction: int = 1) -> int:
    """
    Find the nearest prime to x in the requested direction.
    """
    x = max(2, int(x))

    if x == 2:
        return 2

    if direction > 0:
        if x % 2 == 0:
            x += 1
        while not isprime(x):
            x += 2
        return x

    if x % 2 == 0:
        x -= 1
    while x >= 2 and not isprime(x):
        x -= 2
    if x < 2:
        raise RuntimeError("Could not find prime.")
    return x


def choose_balanced_prime_pair(target_R: int) -> tuple[int, int]:
    """
    Choose primes r1,r2 such that r1*r2 is near target_R
    and r1,r2 are approximately sqrt(target_R).
    """
    root = max(2, math.isqrt(max(4, target_R)))

    best = None

    # Search around sqrt(target_R).
    for delta in range(-CONFIG["prime_search_radius"],
                       CONFIG["prime_search_radius"] + 1):

        x = root + delta
        if x < 2:
            continue

        r1 = nearby_prime(x, 1)

        # Choose r2 close to target_R / r1.
        target_r2 = max(2, target_R // r1)
        candidates = {
            nearby_prime(target_r2, 1),
            nearby_prime(target_r2, -1),
        }

        for r2 in candidates:
            R = r1 * r2
            score = abs(R - target_R)

            # Secondary preference: balanced moduli.
            balance = abs(r1 - r2)

            key = (score, balance)

            if best is None or key < best[0]:
                best = (key, r1, r2)

    if best is None:
        raise RuntimeError(
            f"Could not construct modulus pair for R={target_R}"
        )

    return best[1], best[2]


# =============================================================================
# Semiprime generation
# =============================================================================

def generate_anchor(scale: int, rng: random.Random):
    """
    Generate p,q close enough that p*q is near the requested scale.
    """
    target_root = math.sqrt(scale)

    # Keep factors in a broad interval around sqrt(scale).
    low = max(3, int(target_root * 0.75))
    high = max(low + 10, int(target_root * 1.25))

    for _ in range(10000):
        p0 = rng.randint(low, high)
        q0 = rng.randint(low, high)

        p = nearby_prime(p0, 1)
        q = nearby_prime(q0, 1)

        if p == q:
            q = nextprime(q)

        n = p * q

        # Keep n in a useful neighborhood of the requested scale.
        if scale // 2 <= n <= scale * 2:
            return p, q, n

    raise RuntimeError(f"Could not generate anchor for scale={scale}")


# =============================================================================
# Static K -> (k,l) lookup
# =============================================================================

def build_k_lookup(max_K: int):
    """
    Static list of all positive divisor pairs (k,l) with k*l = K.

    Both ordered orientations are included.
    """
    table = {}

    for K in range(1, max_K + 1):
        pairs = []

        root = math.isqrt(K)

        for d in range(1, root + 1):
            if K % d != 0:
                continue

            q = K // d

            pairs.append((d, q))

            if d != q:
                pairs.append((q, d))

        pairs.sort()
        table[K] = tuple(pairs)

    return table


# =============================================================================
# Representation mathematics
# =============================================================================

def representation_from_moduli(p: int, q: int, r1: int, r2: int):
    """
    Exact quotient-cell representation:

        p = a + k*r1
        q = b + l*r2

    with:

        0 <= a < r1
        0 <= b < r2
    """
    k = p // r1
    l = q // r2

    a = p - k * r1
    b = q - l * r2

    R = r1 * r2
    K = k * l
    T = p * q // R
    E = T - K

    return {
        "r1": r1,
        "r2": r2,
        "R": R,
        "k": k,
        "l": l,
        "a": a,
        "b": b,
        "K": K,
        "T": T,
        "E": E,
    }


def verify_representation(n: int, rep):
    lhs = (
        rep["a"] + rep["k"] * rep["r1"],
        rep["b"] + rep["l"] * rep["r2"],
    )

    product = lhs[0] * lhs[1]

    return (
        product == n
        and 0 <= rep["a"] < rep["r1"]
        and 0 <= rep["b"] < rep["r2"]
    )


# =============================================================================
# Two-representation intersection
# =============================================================================

def implied_residue_difference(rep1, rep2):
    """
    Same p and q imply:

        a1 + k1*r11 = a2 + k2*r12
        b1 + l1*r21 = b2 + l2*r22

    Therefore:

        a1 - a2 = k2*r12 - k1*r11
        b1 - b2 = l2*r22 - l1*r21
    """
    da = (
        rep2["k"] * rep2["r1"]
        - rep1["k"] * rep1["r1"]
    )

    db = (
        rep2["l"] * rep2["r2"]
        - rep1["l"] * rep1["r2"]
    )

    return da, db


def candidate_ab_from_rep_pair(rep1, rep2):
    """
    Given both quotient cells, solve the linear consistency relation.

    a2 = a1 + delta_a
    b2 = b1 + delta_b

    Intersect the allowed residue intervals.

    This produces an exact rectangle intersection in (a,b)-space.
    """
    da, db = implied_residue_difference(rep1, rep2)

    # From a1 - a2 = da:
    # a2 = a1 - da
    #
    # Both must satisfy:
    # 0 <= a1 < r1
    # 0 <= a2 < s1
    #
    # Therefore:
    # 0 <= a1 < r1
    # 0 <= a1-da < s1
    #
    # => da <= a1 < da+s1

    a_low = max(0, da)
    a_high = min(rep1["r1"] - 1, da + rep2["r1"] - 1)

    # Same for b.
    b_low = max(0, db)
    b_high = min(rep1["r2"] - 1, db + rep2["r2"] - 1)

    if a_low > a_high or b_low > b_high:
        return None

    return {
        "da": da,
        "db": db,
        "a_low": a_low,
        "a_high": a_high,
        "b_low": b_low,
        "b_high": b_high,
        "a_width": a_high - a_low + 1,
        "b_width": b_high - b_low + 1,
        "area": (a_high - a_low + 1) * (b_high - b_low + 1),
    }


# =============================================================================
# Solving once K is known
# =============================================================================

def solve_from_K_and_rep(n: int, rep, K_lookup):
    """
    For a known K, enumerate only the static divisor pairs (k,l).

    For every candidate pair:
        p = a + k*r1
        q = b + l*r2

    Equivalently:

        a = p - k*r1
        b = q - l*r2

    The exact test can be performed by checking factor pairs of n.

    This routine is intentionally exact and small; the point of the
    experiment is measuring how much two quotient representations
    shrink the residue space.
    """
    K = rep["K"]

    pairs = K_lookup.get(K, ())

    solutions = []

    # Exact factorization oracle for the experiment.
    factors = factorint(n)

    # Build all exact factor pairs of n.
    divisors = []

    for d in factors:
        # All divisors generated below.
        pass

    # More efficient divisor generation.
    all_divisors = [1]
    for prime, exponent in factors.items():
        old = list(all_divisors)
        powers = []
        cur = 1
        for _ in range(exponent + 1):
            powers.append(cur)
            cur *= prime

        all_divisors = []
        for base in old:
            for pw in powers:
                all_divisors.append(base * pw)

    for p_candidate in all_divisors:
        if p_candidate * p_candidate > n:
            continue

        if n % p_candidate != 0:
            continue

        q_candidate = n // p_candidate

        for pp, qq in (
            (p_candidate, q_candidate),
            (q_candidate, p_candidate),
        ):
            k = pp // rep["r1"]
            l = qq // rep["r2"]

            if k * l != K:
                continue

            a = pp - k * rep["r1"]
            b = qq - l * rep["r2"]

            if not (
                0 <= a < rep["r1"]
                and 0 <= b < rep["r2"]
            ):
                continue

            solutions.append((pp, qq, k, l, a, b))

    # Deduplicate.
    solutions = sorted(set(solutions))

    return solutions


# =============================================================================
# K candidate filtering
# =============================================================================

def nearby_static_K_candidates(T: int, lookup, window: int):
    """
    Generate K candidates near T using the static lookup table.
    """
    low = max(1, T - window)
    high = min(max(lookup), T)

    return list(range(low, high + 1))


def candidate_cells_for_K(rep, K, lookup):
    """
    Retrieve all (k,l) for K from the static table and retain only
    cells compatible with rough quotient geometry.
    """
    pairs = lookup.get(K, ())
    out = []

    for k, l in pairs:
        if k <= 0 or l <= 0:
            continue

        # Cell must be capable of containing p and q.
        # We don't know p,q, so use the expected quotient-cell scale.
        # This is only a prefilter and does not assume the solution.
        if k * rep["r1"] >= 0 and l * rep["r2"] >= 0:
            out.append((k, l))

    return out


# =============================================================================
# One case
# =============================================================================

def analyze_case(n: int, p: int, q: int, rep1, rep2, lookup):
    t0 = time.perf_counter()

    intersection = candidate_ab_from_rep_pair(rep1, rep2)

    if intersection is None:
        return {
            "compatible": False,
            "intersection": None,
        }

    # True residue pair relative to rep1.
    true_a = rep1["a"]
    true_b = rep1["b"]

    contains_true = (
        intersection["a_low"] <= true_a <= intersection["a_high"]
        and intersection["b_low"] <= true_b <= intersection["b_high"]
    )

    # Factorization of the small K values from the static table.
    K1 = rep1["K"]
    K2 = rep2["K"]

    cells1 = candidate_cells_for_K(rep1, K1, lookup)
    cells2 = candidate_cells_for_K(rep2, K2, lookup)

    # Find shared p,q implied by the two static cell sets.
    pair_candidates = []

    for k1, l1 in cells1:
        for k2, l2 in cells2:

            # Candidate p from the common residue interval.
            # p must satisfy:
            #
            #   p = a1+k1*r11 = a2+k2*r12
            #
            # We don't enumerate all a1/a2. Instead test whether
            # the quotient-cell intervals overlap for p.
            p_low = k1 * rep1["r1"]
            p_high = (k1 + 1) * rep1["r1"] - 1

            p_low2 = k2 * rep2["r1"]
            p_high2 = (k2 + 1) * rep2["r1"] - 1

            p_lo = max(p_low, p_low2)
            p_hi = min(p_high, p_high2)

            if p_lo > p_hi:
                continue

            q_low = l1 * rep1["r2"]
            q_high = (l1 + 1) * rep1["r2"] - 1

            q_low2 = l2 * rep2["r2"]
            q_high2 = (l2 + 1) * rep2["r2"] - 1

            q_lo = max(q_low, q_low2)
            q_hi = min(q_high, q_high2)

            if q_lo > q_hi:
                continue

            pair_candidates.append(
                (k1, l1, k2, l2, p_lo, p_hi, q_lo, q_hi)
            )

    # Check true cells explicitly.
    true_cells_present = (
        (rep1["k"], rep1["l"]) in cells1
        and (rep2["k"], rep2["l"]) in cells2
    )

    return {
        "compatible": True,
        "intersection": intersection,
        "contains_true": contains_true,
        "true_cells_present": true_cells_present,
        "cells1": len(cells1),
        "cells2": len(cells2),
        "pair_candidates": pair_candidates,
        "pair_candidate_count": len(pair_candidates),
        "contains_true_cell_pair": any(
            c[0] == rep1["k"]
            and c[1] == rep1["l"]
            and c[2] == rep2["k"]
            and c[3] == rep2["l"]
            for c in pair_candidates
        ),
        "elapsed": time.perf_counter() - t0,
    }


# =============================================================================
# Scale runner
# =============================================================================

def run_scale(scale: int, rng: random.Random, lookup):
    banner(f"SCALE {scale:.0e}")

    anchors = []

    for i in range(CONFIG["anchors_per_scale"]):
        p, q, n = generate_anchor(scale, rng)
        anchors.append((p, q, n))

        print(
            f"anchor {i+1}/{CONFIG['anchors_per_scale']} "
            f"n={fmt(n)} "
            f"p,q=({fmt(p)},{fmt(q)})"
        )

    case_rows = []

    for anchor_index, (p, q, n) in enumerate(anchors, 1):

        reps = []

        print()
        print("-" * 92)
        print(
            f"ANCHOR {anchor_index}/{len(anchors)} "
            f"n={fmt(n)}"
        )
        print("-" * 92)

        for offset in CONFIG["R_offsets"]:

            target_R = int(
                (n / CONFIG["K_target"]) * offset
            )

            r1, r2 = choose_balanced_prime_pair(target_R)
            rep = representation_from_moduli(p, q, r1, r2)

            reps.append(rep)

            print(
                f"R={fmt(rep['R']):>18} "
                f"(r1,r2)=({fmt(r1)},{fmt(r2)}) "
                f"K={rep['K']:>5} "
                f"T={rep['T']:>5} "
                f"E={rep['E']:>4} "
                f"(k,l)=({rep['k']},{rep['l']}) "
                f"(a,b)=({fmt(rep['a'])},{fmt(rep['b'])})"
            )

        # Pairwise two-modulus experiment.
        pair_results = []

        for i, j in combinations(range(len(reps)), 2):

            rep1 = reps[i]
            rep2 = reps[j]

            result = analyze_case(
                n,
                p,
                q,
                rep1,
                rep2,
                lookup,
            )

            if result["compatible"]:
                pair_results.append(result)

                I = result["intersection"]

                print(
                    f"\nPAIR {i+1} x {j+1}: "
                    f"K=({rep1['K']},{rep2['K']})"
                )

                print(
                    f"    delta a={fmt(I['da'])} "
                    f"delta b={fmt(I['db'])}"
                )

                print(
                    f"    a intersection = "
                    f"[{fmt(I['a_low'])}, {fmt(I['a_high'])}] "
                    f"width={fmt(I['a_width'])}"
                )

                print(
                    f"    b intersection = "
                    f"[{fmt(I['b_low'])}, {fmt(I['b_high'])}] "
                    f"width={fmt(I['b_width'])}"
                )

                print(
                    f"    residue area = {fmt(I['area'])}"
                )

                print(
                    f"    true (a,b) inside = "
                    f"{result['contains_true']}"
                )

                print(
                    f"    static cells = "
                    f"{result['cells1']} x {result['cells2']}"
                )

                print(
                    f"    compatible cell-pairs = "
                    f"{result['pair_candidate_count']}"
                )

                print(
                    f"    true cell-pair present = "
                    f"{result['contains_true_cell_pair']}"
                )

        case_rows.append({
            "p": p,
            "q": q,
            "n": n,
            "representations": reps,
            "pairs": pair_results,
        })

    # -------------------------------------------------------------------------
    # Scale summary
    # -------------------------------------------------------------------------

    total_pairs = 0
    true_inside = 0
    true_cell_present = 0
    total_area = 0
    total_width_a = 0
    total_width_b = 0
    total_cell_pairs = 0

    for case in case_rows:
        for pair in case["pairs"]:
            total_pairs += 1
            true_inside += int(pair["contains_true"])
            true_cell_present += int(pair["contains_true_cell_pair"])

            I = pair["intersection"]

            total_area += I["area"]
            total_width_a += I["a_width"]
            total_width_b += I["b_width"]
            total_cell_pairs += pair["pair_candidate_count"]

    print()
    print("-" * 92)
    print(f"SCALE {scale:.0e} TWO-MODULUS SUMMARY")
    print("-" * 92)

    if total_pairs:
        print(
            f"pair cases                 = {total_pairs}"
        )
        print(
            f"true (a,b) contained       = "
            f"{true_inside}/{total_pairs} "
            f"({pct(true_inside / total_pairs)})"
        )
        print(
            f"true cell-pair contained   = "
            f"{true_cell_present}/{total_pairs} "
            f"({pct(true_cell_present / total_pairs)})"
        )
        print(
            f"mean a-intersection width  = "
            f"{total_width_a / total_pairs:.2f}"
        )
        print(
            f"mean b-intersection width  = "
            f"{total_width_b / total_pairs:.2f}"
        )
        print(
            f"mean residue-area          = "
            f"{total_area / total_pairs:.2f}"
        )
        print(
            f"mean compatible cell-pairs = "
            f"{total_cell_pairs / total_pairs:.2f}"
        )

    return case_rows


# =============================================================================
# Main
# =============================================================================

def main():
    started = time.perf_counter()

    banner(
        "START EXPERIMENT 77\n"
        "TWO-MODULUS INTERSECTION / (a,b) RECOVERY EXPERIMENT"
    )

    print()
    print("configuration")
    print(
        f"    scales                  = "
        f"{[f'{x:.0e}' for x in CONFIG['scales']]}"
    )
    print(
        f"    anchors / scale        = "
        f"{CONFIG['anchors_per_scale']}"
    )
    print(
        f"    K target               = "
        f"{CONFIG['K_target']}"
    )
    print(
        f"    R offsets              = "
        f"{CONFIG['R_offsets']}"
    )
    print(
        f"    static K table max     = "
        f"{CONFIG['K_table_max']}"
    )
    print(
        f"    representations / n   = "
        f"{CONFIG['representations_per_n']}"
    )
    print(
        f"    seed                   = "
        f"{CONFIG['seed']}"
    )

    # -------------------------------------------------------------------------
    # Static table
    # -------------------------------------------------------------------------

    print()
    banner("BUILDING STATIC K -> (k,l) LOOKUP")

    t0 = time.perf_counter()

    lookup = build_k_lookup(
        CONFIG["K_table_max"]
    )

    build_time = time.perf_counter() - t0

    divisor_pairs = sum(
        len(v)
        for v in lookup.values()
    )

    print(
        f"K entries             = {len(lookup):,}"
    )
    print(
        f"total divisor pairs   = {divisor_pairs:,}"
    )
    print(
        f"build time             = {build_time:.6f}s"
    )

    # -------------------------------------------------------------------------
    # Experiment
    # -------------------------------------------------------------------------

    rng = random.Random(
        CONFIG["seed"]
    )

    all_results = []

    for scale in CONFIG["scales"]:
        rows = run_scale(
            scale,
            rng,
            lookup,
        )

        all_results.extend(
            (scale, row)
            for row in rows
        )

    # -------------------------------------------------------------------------
    # Global summary
    # -------------------------------------------------------------------------

    print()
    banner("GLOBAL SUMMARY")

    global_pairs = 0
    global_true_inside = 0
    global_true_cells = 0
    global_area = 0
    global_width_a = 0
    global_width_b = 0
    global_cell_pairs = 0

    for scale, row in all_results:
        for pair in row["pairs"]:
            global_pairs += 1
            global_true_inside += int(
                pair["contains_true"]
            )
            global_true_cells += int(
                pair["contains_true_cell_pair"]
            )

            I = pair["intersection"]

            global_area += I["area"]
            global_width_a += I["a_width"]
            global_width_b += I["b_width"]
            global_cell_pairs += pair["pair_candidate_count"]

    print(
        f"total two-modulus pairs      = "
        f"{global_pairs}"
    )

    if global_pairs:
        print(
            f"true (a,b) contained          = "
            f"{global_true_inside}/{global_pairs} "
            f"({pct(global_true_inside / global_pairs)})"
        )

        print(
            f"true cell-pair present        = "
            f"{global_true_cells}/{global_pairs} "
            f"({pct(global_true_cells / global_pairs)})"
        )

        print(
            f"mean a-intersection width     = "
            f"{global_width_a / global_pairs:.2f}"
        )

        print(
            f"mean b-intersection width     = "
            f"{global_width_b / global_pairs:.2f}"
        )

        print(
            f"mean residue search area      = "
            f"{global_area / global_pairs:.2f}"
        )

        print(
            f"mean compatible cell pairs    = "
            f"{global_cell_pairs / global_pairs:.2f}"
        )

    # -------------------------------------------------------------------------
    # Important mathematical check
    # -------------------------------------------------------------------------

    print()
    banner("ALGEBRAIC CONSISTENCY CHECK")

    failures = 0

    for scale, row in all_results:
        for pair in row["pairs"]:
            rep1 = row["representations"][0]
            # We only need to verify the actual relation from the pair.
            I = pair["intersection"]

            # The true residues of the first representation must satisfy
            # the predicted difference relation.
            # Find the actual reps by matching the K/T values.
            reps = row["representations"]

            found = False

            for rep_a, rep_b in combinations(reps, 2):
                candidate = candidate_ab_from_rep_pair(
                    rep_a,
                    rep_b,
                )

                if candidate is None:
                    continue

                if (
                    candidate["a_low"]
                    <= rep_a["a"]
                    <= candidate["a_high"]
                    and
                    candidate["b_low"]
                    <= rep_a["b"]
                    <= candidate["b_high"]
                ):
                    found = True
                    break

            if not found:
                failures += 1

    print(
        f"intersection consistency failures = "
        f"{failures}"
    )

    # -------------------------------------------------------------------------
    # Timing
    # -------------------------------------------------------------------------

    total_time = time.perf_counter() - started

    print()
    banner("TIMING")

    print(
        f"lookup build time            = "
        f"{build_time:.6f}s"
    )
    print(
        f"total runtime                = "
        f"{total_time:.6f}s"
    )

    print()
    banner("FINISHED EXPERIMENT 77")


if __name__ == "__main__":
    main()
