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

ANCHORS_PER_SCALE = 12

MOD_MIN = 300
MOD_MAX = 3000

CLOSE_RATIO = 0.20

# Test these theoretical R scales.
R_EXPONENTS = (
    0.25,
    1.0 / 3.0,
    0.50,
    2.0 / 3.0,
)

# Deterministic auxiliary numbers.
S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Maximum number of displayed matches.
MAX_MATCH_EXAMPLES = 3


# =============================================================================
# PRIME POOL
# =============================================================================

def build_prime_pool():
    return list(sp.primerange(MOD_MIN, MOD_MAX + 1))


# =============================================================================
# MODULUS PAIRS
# =============================================================================

def build_modulus_pairs(primes):
    """
    Construct all prime pairs satisfying the close-ratio condition.

    The pair may be ordered later, but storage is canonical.
    """
    pairs = []

    for i, r1 in enumerate(primes):
        for r2 in primes[i:]:
            lo = min(r1, r2)
            hi = max(r1, r2)

            if (hi - lo) / lo <= CLOSE_RATIO:
                pairs.append((r1, r2))

    return pairs


def choose_best_pair(target_R, modulus_pairs):
    """
    Find the feasible prime pair whose product is closest to target_R.

    Because the modulus range is fixed, the smallest possible product is:

        MOD_MIN^2

    and the largest is:

        MOD_MAX^2.

    Therefore theoretical R targets outside that interval are clamped.
    """
    min_R = MOD_MIN * MOD_MIN
    max_R = MOD_MAX * MOD_MAX

    clamped_R = max(min_R, min(max_R, target_R))

    best = None

    for r1, r2 in modulus_pairs:
        R = r1 * r2

        # Primary objective: absolute distance.
        distance = abs(R - clamped_R)

        # Secondary objective: relative distance.
        rel = distance / clamped_R

        key = (
            distance,
            rel,
            abs(r1 - r2),
            r1,
            r2,
        )

        if best is None or key < best[0]:
            best = (key, r1, r2, R)

    if best is None:
        raise RuntimeError("No feasible modulus pairs were generated.")

    return {
        "target_R": target_R,
        "effective_R": clamped_R,
        "r1": best[1],
        "r2": best[2],
        "R": best[3],
        "clamped": clamped_R != target_R,
    }


# =============================================================================
# ANCHORS
# =============================================================================

def generate_anchor(scale, rng):
    """
    Generate two reasonably balanced primes whose product is near scale.
    """
    center = math.isqrt(scale)

    lo = max(10, int(center * 0.70))
    hi = int(center * 1.30)

    p = int(sp.randprime(lo, hi))
    q = int(sp.randprime(lo, hi))

    if p == q:
        q = int(sp.nextprime(q))

    return p, q, p * q


# =============================================================================
# FACTORIZATION
# =============================================================================

def factor_pairs(value):
    """
    Factor a value with SymPy and return all divisor pairs d,e such that:

        d*e=value

    and d <= e.

    This is intentionally exact.
    """
    fac = sp.factorint(value)

    divisors = [1]

    for prime, exponent in fac.items():
        old = list(divisors)

        powers = [1]
        current = 1

        for _ in range(exponent):
            current *= prime
            powers.append(current)

        divisors = []

        for d in old:
            for pw in powers:
                divisors.append(d * pw)

    divisors = sorted(set(divisors))

    pairs = []

    root = math.isqrt(value)

    for d in divisors:
        if d > root:
            break

        if value % d == 0:
            pairs.append((d, value // d))

    return pairs


# =============================================================================
# AUXILIARY TEST
# =============================================================================

def deterministic_x(n, S):
    return (-n) % S


def test_auxiliary(
    n,
    p,
    q,
    r1,
    r2,
    S,
):
    """
    Construct n+x divisible by S and test every exact factor pair.
    """
    x = deterministic_x(n, S)
    nx = n + x

    if nx <= 0:
        return {
            "x": x,
            "nx": nx,
            "pairs": 0,
            "matches": 0,
            "same_cell": 0,
            "cross_cell": 0,
            "min_delta": None,
            "examples": [],
        }

    pairs = factor_pairs(nx)

    k = p // r1
    l = q // r2
    K = k * l

    matches = 0
    same_cell = 0
    cross_cell = 0

    min_delta = None
    examples = []

    for px, qx in pairs:

        # Test both orientations.
        for pax, qax in ((px, qx), (qx, px)):

            kx = pax // r1
            lx = qax // r2
            Kx = kx * lx

            delta = abs(Kx - K)

            if min_delta is None or delta < min_delta:
                min_delta = delta

            if Kx == K:
                matches += 1

                if kx == k and lx == l:
                    same_cell += 1
                else:
                    cross_cell += 1

                if len(examples) < MAX_MATCH_EXAMPLES:
                    examples.append(
                        {
                            "px": pax,
                            "qx": qax,
                            "kx": kx,
                            "lx": lx,
                            "Kx": Kx,
                        }
                    )

    return {
        "x": x,
        "nx": nx,
        "pairs": len(pairs),
        "matches": matches,
        "same_cell": same_cell,
        "cross_cell": cross_cell,
        "min_delta": min_delta,
        "examples": examples,
    }


# =============================================================================
# ONE ANCHOR / ONE R SCALE
# =============================================================================

def run_R_case(
    n,
    p,
    q,
    r1,
    r2,
):
    K = (p // r1) * (q // r2)

    total_pairs = 0
    total_matches = 0
    total_same = 0
    total_cross = 0
    S_hits = 0

    min_delta = None

    examples = []

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
        total_same += result["same_cell"]
        total_cross += result["cross_cell"]

        if result["matches"] > 0:
            S_hits += 1

            if len(examples) < MAX_MATCH_EXAMPLES:
                examples.append(
                    {
                        "S": S,
                        **result,
                    }
                )

        if result["min_delta"] is not None:
            if min_delta is None or result["min_delta"] < min_delta:
                min_delta = result["min_delta"]

    return {
        "K": K,
        "pairs": total_pairs,
        "matches": total_matches,
        "same": total_same,
        "cross": total_cross,
        "S_hits": S_hits,
        "min_delta": min_delta,
        "examples": examples,
    }


# =============================================================================
# SCALE
# =============================================================================

def run_scale(scale, modulus_pairs, rng):
    print()
    print("=" * 100)
    print(f"SCALE {scale:.0e}")
    print("=" * 100)

    scale_results = []

    for anchor_id in range(1, ANCHORS_PER_SCALE + 1):

        p, q, n = generate_anchor(scale, rng)

        print(
            f"anchor {anchor_id:2d}/{ANCHORS_PER_SCALE} "
            f"n={n:,}"
        )

        for exponent in R_EXPONENTS:

            theoretical_R = max(
                1,
                int(round(n ** exponent))
            )

            modulus = choose_best_pair(
                theoretical_R,
                modulus_pairs,
            )

            result = run_R_case(
                n,
                p,
                q,
                modulus["r1"],
                modulus["r2"],
            )

            record = {
                "scale": scale,
                "anchor": anchor_id,
                "n": n,
                "p": p,
                "q": q,
                "exponent": exponent,
                "theoretical_R": theoretical_R,
                **modulus,
                **result,
            }

            scale_results.append(record)

    return scale_results


# =============================================================================
# SUMMARY
# =============================================================================

def summarize(results):
    grouped = defaultdict(list)

    for rec in results:
        grouped[round(rec["exponent"], 6)].append(rec)

    print()
    print("-" * 100)
    print("R-SCALE SUMMARY")
    print("-" * 100)

    print(
        "exp       target-R        effective-R       "
        "auxPairs       Kx=K       same      cross      S-hit%"
    )

    print("-" * 100)

    for exponent in sorted(grouped):

        rows = grouped[exponent]

        total_pairs = sum(x["pairs"] for x in rows)
        total_matches = sum(x["matches"] for x in rows)
        total_same = sum(x["same"] for x in rows)
        total_cross = sum(x["cross"] for x in rows)

        hit_count = sum(
            1
            for x in rows
            if x["S_hits"] > 0
        )

        target_mean = (
            sum(x["target_R"] for x in rows)
            / len(rows)
        )

        effective_mean = (
            sum(x["effective_R"] for x in rows)
            / len(rows)
        )

        hit_rate = hit_count / len(rows)

        print(
            f"{exponent:0.3f}   "
            f"{target_mean:14.3f} "
            f"{effective_mean:15.3f} "
            f"{total_pairs:12,d} "
            f"{total_matches:11,d} "
            f"{total_same:9,d} "
            f"{total_cross:9,d} "
            f"{100.0 * hit_rate:8.2f}%"
        )


def print_match_examples(results):
    examples = []

    for rec in results:
        for item in rec["examples"]:
            examples.append(
                (
                    rec,
                    item,
                )
            )

    print()
    print("-" * 100)
    print("Kx = K EXAMPLES")
    print("-" * 100)

    if not examples:
        print("none")
        return

    for rec, item in examples[:12]:

        print(
            f"scale={rec['scale']:.0e} "
            f"n={rec['n']:,} "
            f"exp={rec['exponent']:.3f} "
            f"R={rec['R']:,} "
            f"(r1,r2)=({rec['r1']},{rec['r2']}) "
            f"K={rec['K']:,} "
            f"S={item['S']} "
            f"x={item['x']:,}"
        )

        for e in item["examples"][:3]:
            print(
                f"    ({e['px']:,},{e['qx']:,}) "
                f"-> ({e['kx']},{e['lx']}) "
                f"Kx={e['Kx']:,}"
            )


# =============================================================================
# FINAL CROSS-SCALE TABLE
# =============================================================================

def cross_scale(results):
    grouped = defaultdict(list)

    for rec in results:
        grouped[rec["scale"]].append(rec)

    print()
    print("=" * 100)
    print("CROSS-SCALE COMPARISON")
    print("=" * 100)

    print(
        "scale       tests       auxPairs       Kx=K       "
        "same       cross       Kx=K/auxPair"
    )

    print("-" * 100)

    for scale in sorted(grouped):

        rows = grouped[scale]

        tests = len(rows)
        aux_pairs = sum(x["pairs"] for x in rows)
        matches = sum(x["matches"] for x in rows)
        same = sum(x["same"] for x in rows)
        cross = sum(x["cross"] for x in rows)

        ratio = matches / aux_pairs if aux_pairs else 0.0

        print(
            f"{scale:>5.0e} "
            f"{tests:10,d} "
            f"{aux_pairs:14,d} "
            f"{matches:11,d} "
            f"{same:10,d} "
            f"{cross:10,d} "
            f"{ratio:17.10f}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 100)
    print("MODULUS-SCALE / Kx=K EXPERIMENT")
    print("=" * 100)

    print(f"scales                 = {[f'{x:.0e}' for x in SCALES]}")
    print(f"anchors / scale        = {ANCHORS_PER_SCALE}")
    print(f"modulus range          = {MOD_MIN} .. {MOD_MAX}")
    print(f"close ratio             = {CLOSE_RATIO:.2%}")
    print(f"R exponents             = {R_EXPONENTS}")
    print(f"S values                = {S_VALUES}")
    print(f"seed                    = {SEED}")

    print()
    print("=" * 100)
    print("BUILDING MODULUS PRIME POOL")
    print("=" * 100)

    primes = build_prime_pool()

    print(f"modulus primes          = {len(primes)}")

    print()
    print("=" * 100)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 100)

    modulus_pairs = build_modulus_pairs(primes)

    print(f"close modulus pairs     = {len(modulus_pairs)}")

    min_R = MOD_MIN * MOD_MIN
    max_R = MOD_MAX * MOD_MAX

    print(f"feasible R range        = {min_R:,} .. {max_R:,}")

    print()
    print(
        "Targets outside this interval are clamped to the nearest "
        "feasible R."
    )

    rng = random.Random(SEED)

    all_results = []

    start = time.perf_counter()

    for scale in SCALES:

        scale_results = run_scale(
            scale,
            modulus_pairs,
            rng,
        )

        all_results.extend(scale_results)

        summarize(scale_results)
        print_match_examples(scale_results)

    elapsed = time.perf_counter() - start

    cross_scale(all_results)

    # =========================================================================
    # CLAMPING REPORT
    # =========================================================================

    print()
    print("=" * 100)
    print("R-TARGET CLAMPING")
    print("=" * 100)

    for scale in SCALES:

        print(f"{scale:.0e}:")

        for exponent in R_EXPONENTS:

            target = int(round(scale ** exponent))

            if target < min_R:
                print(
                    f"    n^{exponent:.3f}: "
                    f"{target:,} -> clamped upward to >= {min_R:,}"
                )

            elif target > max_R:
                print(
                    f"    n^{exponent:.3f}: "
                    f"{target:,} -> clamped downward to <= {max_R:,}"
                )

            else:
                print(
                    f"    n^{exponent:.3f}: "
                    f"{target:,} -> feasible"
                )

    # =========================================================================
    # MATHEMATICAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        r"""
For

    p = a + k*r1
    q = b + l*r2

define

    R = r1*r2
    K = k*l
    T = floor(n/R)
    E = T-K.

The hypothesis being tested is that the density of auxiliary
collisions

    Kx = K

depends on the scale of R.

The theoretical target scales are

    R ~ n^(1/4)
    R ~ n^(1/3)
    R ~ n^(1/2)
    R ~ n^(2/3).

However, r1 and r2 are restricted to the actual modulus range

    MOD_MIN <= r1,r2 <= MOD_MAX,

so not every theoretical R is realizable.

The experiment therefore compares each theoretical target against
the nearest feasible prime product

    R = r1*r2.

For an auxiliary value

    n+x,

the exact factorization gives

    n+x = px*qx

and therefore

    kx = floor(px/r1)
    lx = floor(qx/r2)

and

    Kx = kx*lx.

The principal event is

    Kx = K.

It is divided into:

    same-cell:
        (kx,lx) = (k,l)

and:

    cross-cell:
        (kx,lx) != (k,l)
        but
        kx*lx = k*l.

The important comparison is therefore not merely whether
Kx=K occurs, but whether its frequency changes systematically
with the modulus-product scale.

No factorization shortcut is assumed.

SymPy is used only as an exact experimental oracle for
factorizing the auxiliary n+x values.

The final experiment does not claim that a high collision rate
implies an algorithm. It only tests whether a preferred
modulus scale exists.
"""
    )

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(f"total runtime              = {elapsed:.3f}s")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()