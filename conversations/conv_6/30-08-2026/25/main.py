#!/usr/bin/env python3

import math
import random
import time
from collections import Counter


# ================================================================================================
# CONFIG
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

ANCHORS = 300

SEED = 1_511_464_998

MAX_D = FACTOR_MAX - FACTOR_MIN

# Both factors are odd primes, so:
#
#   p+q  is even
#   q-p  is even
#
# Therefore d only needs even values.
D_MIN = 0
D_MAX = MAX_D


# ================================================================================================
# SIEVE
# ================================================================================================

def sieve(limit):
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    root = math.isqrt(limit)

    for p in range(2, root + 1):

        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start::p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if a[i]]


# ================================================================================================
# ANCHORS
# ================================================================================================

def generate_anchors(primes, count, rng):

    result = []
    seen = set()

    while len(result) < count:

        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in seen:
            continue

        seen.add(pair)

        result.append(
            (
                p,
                q,
                p * q,
            )
        )

    return result


# ================================================================================================
# MAXIMUM-PRODUCT CLOSE PRIME TRIPLE
# ================================================================================================

def select_maximum_product_triple(n, modulus_primes):

    best = None
    best_R = -1

    L = len(modulus_primes)

    for i in range(L):

        r1 = modulus_primes[i]

        for j in range(i + 1, L):

            r2 = modulus_primes[j]

            if r2 > r1 * (1.0 + CLOSE_RATIO):
                break

            for k in range(j + 1, L):

                r3 = modulus_primes[k]

                if r3 > r1 * (1.0 + CLOSE_RATIO):
                    break

                R = r1 * r2 * r3

                if R < n and R > best_R:

                    best_R = R
                    best = (r1, r2, r3)

    return best


# ================================================================================================
# QUADRATIC RESIDUE TABLE
# ================================================================================================

def build_qr_table(p):

    qr = bytearray(p)

    for x in range(p):
        qr[(x * x) % p] = 1

    return qr


# ================================================================================================
# EXACT RECOVERY
# ================================================================================================

def recover_from_d(n, d):

    D = 4 * n + d * d

    s = math.isqrt(D)

    if s * s != D:
        return None

    if (s - d) & 1:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if not (
        FACTOR_MIN <= p <= FACTOR_MAX
        and
        FACTOR_MIN <= q <= FACTOR_MAX
    ):
        return None

    if p * q != n:
        return None

    return p, q, s


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    t0 = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME MODULAR DIFFERENCE DOMAIN EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in primes
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    print()

    # ============================================================================================
    # ANCHORS
    # ============================================================================================

    anchors = generate_anchors(
        factor_primes,
        ANCHORS,
        rng,
    )

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    records = []

    for i, (p, q, n) in enumerate(anchors, 1):

        mods = select_maximum_product_triple(
            n,
            modulus_primes,
        )

        if mods is None:
            continue

        r1, r2, r3 = mods

        R = r1 * r2 * r3
        gap = n - R
        d_actual = q - p

        records.append(
            {
                "id": i,
                "p": p,
                "q": q,
                "n": n,
                "mods": mods,
                "R": R,
                "gap": gap,
                "d": d_actual,
            }
        )

        if i % 25 == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    print()
    print(
        f"usable anchors = {len(records)}"
    )
    print()

    # ============================================================================================
    # QR TABLES
    # ============================================================================================

    print("=" * 100)
    print("BUILDING QR TABLES")
    print("=" * 100)

    used_moduli = sorted(
        {
            r
            for rec in records
            for r in rec["mods"]
        }
    )

    qr_tables = {}

    for i, r in enumerate(used_moduli, 1):

        qr_tables[r] = build_qr_table(r)

        print(
            f"QR table {i:3d}/{len(used_moduli)} "
            f"r={r}"
        )

    print()

    # ============================================================================================
    # GLOBAL COUNTERS
    # ============================================================================================

    total_d = 0

    total_qr1 = 0
    total_qr2 = 0
    total_qr3 = 0

    total_exact_after_qr1 = 0
    total_exact_after_qr2 = 0
    total_exact_after_qr3 = 0

    total_plain_exact = 0

    recovered = 0

    d_positions = []
    d_rankings = []

    per_anchor = []

    # ============================================================================================
    # MAIN DIFFERENCE SEARCH
    # ============================================================================================

    print("=" * 100)
    print("RUNNING DIFFERENCE-DOMAIN SEARCH")
    print("=" * 100)

    for anchor_index, rec in enumerate(records, 1):

        n = rec["n"]
        actual_d = rec["d"]

        r1, r2, r3 = rec["mods"]

        qr1 = qr_tables[r1]
        qr2 = qr_tables[r2]
        qr3 = qr_tables[r3]

        # ----------------------------------------------------------------------------------------
        # Plain exact difference search
        # ----------------------------------------------------------------------------------------

        plain_hits = []

        for d in range(D_MIN, D_MAX + 1, 2):

            total_d += 1

            D = 4 * n + d * d
            s = math.isqrt(D)

            if s * s != D:
                continue

            pair = recover_from_d(n, d)

            if pair is not None:
                plain_hits.append(
                    (
                        d,
                        pair[0],
                        pair[1],
                    )
                )

        total_plain_exact += len(plain_hits)

        # ----------------------------------------------------------------------------------------
        # Modular filtering
        # ----------------------------------------------------------------------------------------

        count1 = 0
        count2 = 0
        count3 = 0

        exact1 = []
        exact2 = []
        exact3 = []

        for d in range(D_MIN, D_MAX + 1, 2):

            D = 4 * n + d * d

            # First modulus.
            if not qr1[D % r1]:
                continue

            count1 += 1

            # Exact square after one QR filter.
            s = math.isqrt(D)

            if s * s == D:

                pair = recover_from_d(n, d)

                if pair is not None:
                    exact1.append(
                        (
                            d,
                            pair[0],
                            pair[1],
                        )
                    )

            # Second modulus.
            if not qr2[D % r2]:
                continue

            count2 += 1

            # Exact square after two QR filters.
            s = math.isqrt(D)

            if s * s == D:

                pair = recover_from_d(n, d)

                if pair is not None:
                    exact2.append(
                        (
                            d,
                            pair[0],
                            pair[1],
                        )
                    )

            # Third modulus.
            if not qr3[D % r3]:
                continue

            count3 += 1

            # Exact square after three QR filters.
            s = math.isqrt(D)

            if s * s == D:

                pair = recover_from_d(n, d)

                if pair is not None:
                    exact3.append(
                        (
                            d,
                            pair[0],
                            pair[1],
                        )
                    )

        total_qr1 += count1
        total_qr2 += count2
        total_qr3 += count3

        total_exact_after_qr1 += len(exact1)
        total_exact_after_qr2 += len(exact2)
        total_exact_after_qr3 += len(exact3)

        # ----------------------------------------------------------------------------------------
        # Position of actual d inside the modular candidate stream
        # ----------------------------------------------------------------------------------------

        qr3_values = []

        for d in range(D_MIN, D_MAX + 1, 2):

            D = 4 * n + d * d

            if not qr1[D % r1]:
                continue

            if not qr2[D % r2]:
                continue

            if not qr3[D % r3]:
                continue

            qr3_values.append(d)

        try:
            d_rank = qr3_values.index(actual_d) + 1
        except ValueError:
            d_rank = None

        if d_rank is not None:

            d_positions.append(
                (
                    d_rank,
                    len(qr3_values),
                )
            )

            d_rankings.append(
                d_rank / len(qr3_values)
            )

        # ----------------------------------------------------------------------------------------
        # Recovery
        # ----------------------------------------------------------------------------------------

        recovered_this_anchor = False

        for d, p2, q2 in exact3:

            if (
                d == actual_d
                and
                p2 == rec["p"]
                and
                q2 == rec["q"]
            ):
                recovered_this_anchor = True
                break

        if recovered_this_anchor:
            recovered += 1

        per_anchor.append(
            {
                "id": rec["id"],
                "p": rec["p"],
                "q": rec["q"],
                "n": n,
                "mods": rec["mods"],
                "R": rec["R"],
                "gap": rec["gap"],
                "d": actual_d,
                "qr1": count1,
                "qr2": count2,
                "qr3": count3,
                "exact1": len(exact1),
                "exact2": len(exact2),
                "exact3": len(exact3),
                "d_rank": d_rank,
                "d_total": len(qr3_values),
            }
        )

        if anchor_index % 25 == 0:
            print(
                f"anchor {anchor_index:3d}/{len(records)}"
            )

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    N = len(records)

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"anchors analyzed             = {N}"
    )

    print(
        f"even d values / anchor       = "
        f"{total_d / N:,.3f}"
    )

    print(
        f"after QR(r1) / anchor        = "
        f"{total_qr1 / N:,.3f}"
    )

    print(
        f"after QR(r1,r2) / anchor     = "
        f"{total_qr2 / N:,.3f}"
    )

    print(
        f"after QR(r1,r2,r3) / anchor  = "
        f"{total_qr3 / N:,.3f}"
    )

    print()

    print(
        f"plain exact solutions         = "
        f"{total_plain_exact}"
    )

    print(
        f"exact after QR1              = "
        f"{total_exact_after_qr1}"
    )

    print(
        f"exact after QR2              = "
        f"{total_exact_after_qr2}"
    )

    print(
        f"exact after QR3              = "
        f"{total_exact_after_qr3}"
    )

    print()

    # ============================================================================================
    # REDUCTION
    # ============================================================================================

    print("=" * 100)
    print("DIFFERENCE-DOMAIN REDUCTION")
    print("=" * 100)

    base = total_d

    print(
        f"QR1 reduction                = "
        f"{100.0 * (1.0 - total_qr1 / base):.6f}%"
    )

    print(
        f"QR1+QR2 reduction            = "
        f"{100.0 * (1.0 - total_qr2 / base):.6f}%"
    )

    print(
        f"QR1+QR2+QR3 reduction        = "
        f"{100.0 * (1.0 - total_qr3 / base):.6f}%"
    )

    print()

    # ============================================================================================
    # EXPECTED RANDOM BASELINE
    # ============================================================================================

    print("=" * 100)
    print("RANDOM-RESIDUE EXPECTATION")
    print("=" * 100)

    print()
    print(
        "For a random nonzero residue modulo a large prime:"
    )

    print()
    print(
        "    P(QR) ≈ 1/2"
    )

    print()

    print(
        "Therefore three independent QR tests should retain "
        "approximately 1/8 of the difference domain."
    )

    expected = (
        total_d / 8.0
        if total_d
        else 0.0
    )

    print(
        f"expected after 3 QR tests      ≈ "
        f"{expected:,.3f}"
    )

    print(
        f"observed after 3 QR tests      = "
        f"{total_qr3:,.0f}"
    )

    print()

    # ============================================================================================
    # ACTUAL d POSITION
    # ============================================================================================

    print("=" * 100)
    print("ACTUAL d POSITION INSIDE MODULAR CANDIDATES")
    print("=" * 100)

    if d_rankings:

        print(
            f"mean percentile               = "
            f"{sum(d_rankings) / len(d_rankings):.6f}"
        )

        ordered = sorted(d_rankings)

        mid = len(ordered) // 2

        if len(ordered) % 2:
            median = ordered[mid]
        else:
            median = (
                ordered[mid - 1] +
                ordered[mid]
            ) / 2

        print(
            f"median percentile             = "
            f"{median:.6f}"
        )

        print(
            f"best percentile               = "
            f"{min(d_rankings):.6f}"
        )

        print(
            f"worst percentile              = "
            f"{max(d_rankings):.6f}"
        )

    print()

    # ============================================================================================
    # R/n RELATIONSHIP
    # ============================================================================================

    print("=" * 100)
    print("R/n VS DIFFERENCE SEARCH LOAD")
    print("=" * 100)

    # Simple Pearson correlation.
    xs = []
    ys = []

    for row in per_anchor:

        xs.append(
            row["R"] / row["n"]
        )

        ys.append(
            row["qr3"]
        )

    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)

    num = sum(
        (x - mx) * (y - my)
        for x, y in zip(xs, ys)
    )

    denx = math.sqrt(
        sum(
            (x - mx) ** 2
            for x in xs
        )
    )

    deny = math.sqrt(
        sum(
            (y - my) ** 2
            for y in ys
        )
    )

    corr = (
        num / (denx * deny)
        if denx and deny
        else 0.0
    )

    print(
        f"corr(R/n, QR3 candidates) = "
        f"{corr:.6f}"
    )

    print()

    # ============================================================================================
    # GAP CORRELATION
    # ============================================================================================

    gaps = [
        row["gap"]
        for row in per_anchor
    ]

    loads = [
        row["qr3"]
        for row in per_anchor
    ]

    mg = sum(gaps) / len(gaps)
    ml = sum(loads) / len(loads)

    num = sum(
        (g - mg) * (l - ml)
        for g, l in zip(gaps, loads)
    )

    deng = math.sqrt(
        sum(
            (g - mg) ** 2
            for g in gaps
        )
    )

    denl = math.sqrt(
        sum(
            (l - ml) ** 2
            for l in loads
        )
    )

    gap_corr = (
        num / (deng * denl)
        if deng and denl
        else 0.0
    )

    print("=" * 100)
    print("GAP VS DIFFERENCE-DOMAIN LOAD")
    print("=" * 100)

    print(
        f"corr(gap, QR3 candidates) = "
        f"{gap_corr:.6f}"
    )

    print()

    # ============================================================================================
    # MOST AGGRESSIVE CASES
    # ============================================================================================

    print("=" * 100)
    print("MOST AGGRESSIVE DIFFERENCE COLLAPSES")
    print("=" * 100)

    strongest = sorted(
        per_anchor,
        key=lambda x: (
            x["qr3"],
            x["qr2"],
            x["qr1"],
        )
    )

    for row in strongest[:20]:

        r1, r2, r3 = row["mods"]

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"d={row['d']:,} "
            f"mods=({r1},{r2},{r3}) "
            f"R/n={row['R']/row['n']:.10f} "
            f"gap={row['gap']:,} "
            f"QR1={row['qr1']} "
            f"QR2={row['qr2']} "
            f"QR3={row['qr3']} "
            f"dRank={row['d_rank']}/"
            f"{row['d_total']}"
        )

    # ============================================================================================
    # WORST CASES
    # ============================================================================================

    print()
    print("=" * 100)
    print("LEAST AGGRESSIVE DIFFERENCE COLLAPSES")
    print("=" * 100)

    weakest = sorted(
        per_anchor,
        key=lambda x: (
            -x["qr3"],
            -x["qr2"],
        )
    )

    for row in weakest[:20]:

        r1, r2, r3 = row["mods"]

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"d={row['d']:,} "
            f"mods=({r1},{r2},{r3}) "
            f"R/n={row['R']/row['n']:.10f} "
            f"gap={row['gap']:,} "
            f"QR3={row['qr3']} "
            f"dRank={row['d_rank']}/"
            f"{row['d_total']}"
        )

    # ============================================================================================
    # ANCHOR TABLE
    # ============================================================================================

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q       d       R/n       "
        "QR1   QR2   QR3   DRANK"
    )

    print("-" * 100)

    for row in per_anchor:

        if row["d_rank"] is None:
            rank_text = "NONE"
        else:
            rank_text = (
                f"{row['d_rank']}/"
                f"{row['d_total']}"
            )

        print(
            f"{row['id']:3d} "
            f"{row['p']:8,d} "
            f"{row['q']:8,d} "
            f"{row['d']:7,d} "
            f"{row['R']/row['n']:.7f} "
            f"{row['qr1']:5d} "
            f"{row['qr2']:5d} "
            f"{row['qr3']:5d} "
            f"{rank_text:>12}"
        )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"correct actual factors recovered = "
        f"{recovered}/{N}"
    )

    print(
        f"recovery rate                    = "
        f"{100.0 * recovered / N:.4f}%"
    )

    # ============================================================================================
    # MATHEMATICAL INTERPRETATION
    # ============================================================================================

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print()
    print("We now work with:")
    print()
    print("    d = |p-q|")
    print()
    print("and:")
    print()
    print("    s^2 = 4n + d^2")
    print()

    print(
        "For the true factor pair, s is therefore exactly "
        "recoverable from d whenever 4n+d² is a square."
    )

    print()
    print("Modulo each close prime:")
    print()
    print("    s² ≡ 4n + d² (mod r_i)")
    print()
    print(
        "so every true d must satisfy three quadratic-residue "
        "conditions."
    )

    print()
    print("The search therefore becomes:")
    print()
    print("    all possible d")
    print("          |")
    print("          v")
    print("    QR(r1)")
    print("          |")
    print("          v")
    print("    QR(r2)")
    print("          |")
    print("          v")
    print("    QR(r3)")
    print("          |")
    print("          v")
    print("    exact square")
    print("          |")
    print("          v")
    print("    p=(s-d)/2, q=(s+d)/2")
    print()

    print(
        "The crucial comparison with the previous experiment is that "
        "there is no p enumeration and no s enumeration."
    )

    print()
    print(
        "The difference domain contains only 45,001 even candidates "
        "for this factor interval."
    )

    print()
    print(
        "The strongest possible observation would be that the "
        "three-modulus QR intersection is substantially smaller "
        "than the expected 1/8 random survival rate, or that the "
        "actual d repeatedly appears unusually early in the surviving "
        "candidate ordering."
    )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    elapsed = time.perf_counter() - t0

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime = {elapsed:.3f} seconds"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
