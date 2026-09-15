#!/usr/bin/env python3

import math
import random
import time
from collections import Counter


# ==================================================================================================
# CONFIG
# ==================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

ANCHORS = 300

SEED = 1_511_464_998


# ==================================================================================================
# SIEVE
# ==================================================================================================

def sieve(limit):
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start::p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if a[i]]


# ==================================================================================================
# TRIPLE SELECTION
# ==================================================================================================

def select_triple(n, modulus_primes):
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


# ==================================================================================================
# ANCHORS
# ==================================================================================================

def generate_anchors(primes, count, rng):
    seen = set()
    out = []

    while len(out) < count:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))
        out.append((p, q, p * q))

    return out


# ==================================================================================================
# BASIC NUMBER THEORY
# ==================================================================================================

def factor_pair_from_s(n, s):
    """
    If s = p+q then

        x^2 - s*x + n = 0

    and

        D = s^2 - 4n = (p-q)^2.
    """

    D = s * s - 4 * n

    if D < 0:
        return None

    root = math.isqrt(D)

    if root * root != D:
        return None

    if (s + root) % 2:
        return None

    p = (s - root) // 2
    q = (s + root) // 2

    if p * q != n:
        return None

    return p, q


# ==================================================================================================
# CRT
# ==================================================================================================

def crt3(residues, moduli):
    x = 0
    R = math.prod(moduli)

    for a, r in zip(residues, moduli):
        m = R // r
        x += a * m * pow(m, -1, r)

    return x % R


# ==================================================================================================
# DIRECT CRT SEARCH
# ==================================================================================================

def direct_crt_candidates(n, mods, factor_primes):
    r1, r2, r3 = mods

    R = r1 * r2 * r3

    hits = []

    for p in factor_primes:

        if math.gcd(p, R) != 1:
            continue

        q1 = (n % r1) * pow(p, -1, r1) % r1
        q2 = (n % r2) * pow(p, -1, r2) % r2
        q3 = (n % r3) * pow(p, -1, r3) % r3

        Q = crt3(
            (q1, q2, q3),
            mods
        )

        if FACTOR_MIN <= Q <= FACTOR_MAX:
            hits.append((p, Q))

    return hits


# ==================================================================================================
# GAP DECOMPOSITION
# ==================================================================================================

def gap_decomposition(n, p, q, mods):
    r1, r2, r3 = mods

    R = r1 * r2 * r3
    g = n - R

    # Basic symmetric quantities.
    s = p + q
    d = abs(p - q)

    # Compare factors against modulus scale.
    rp = [p % r for r in mods]
    rq = [q % r for r in mods]

    # Quotients.
    kp = [(p - a) // r for p, a, r in zip([p] * 3, rp, mods)]
    kq = [(q - a) // r for q, a, r in zip([q] * 3, rq, mods)]

    # Product decomposition:
    #
    # p = k_i r_i + a_i
    # q = l_i r_i + b_i
    #
    # n - R can be compared against the lower-order terms.
    terms = []

    for r, a, b, k, l in zip(mods, rp, rq, kp, kq):

        terms.append(
            {
                "r": r,
                "a": a,
                "b": b,
                "k": k,
                "l": l,
                "ab": a * b,
                "k_b": k * b,
                "l_a": l * a,
            }
        )

    # Distance from the three-modulus product.
    #
    # This is intentionally expressed in several ways to look for
    # small or repeated components.
    #
    # For each i:
    #
    # n = (k_i r_i + a_i)(l_i r_i + b_i)
    #
    # = k_i l_i r_i^2 + (k_i b_i+l_i a_i)r_i + a_i b_i
    #
    # so:
    #
    # g = n - R
    #
    # may have useful representations involving r_i.

    local_values = []

    for t in terms:
        r = t["r"]
        k = t["k"]
        l = t["l"]
        a = t["a"]
        b = t["b"]

        value = (
            k * l * r * r
            + (k * b + l * a) * r
            + a * b
            - R
        )

        local_values.append(value)

    # CRT residual of the gap.
    gap_residues = tuple(g % r for r in mods)

    # gcds of the gap with each modulus and factor-derived quantities.
    gap_gcds = tuple(math.gcd(g, r) for r in mods)

    return {
        "R": R,
        "g": g,
        "s": s,
        "d": d,
        "p_res": tuple(rp),
        "q_res": tuple(rq),
        "p_quot": tuple(kp),
        "q_quot": tuple(kq),
        "gap_res": gap_residues,
        "gap_gcds": gap_gcds,
        "local_values": tuple(local_values),
    }


# ==================================================================================================
# GAP-DERIVED DIVISIBILITY SEARCH
# ==================================================================================================

def gap_divisibility_search(n, mods, factor_primes):
    """
    Explore consequences of

        pq = R + g

    where R = r1*r2*r3 and g = n-R.

    For every divisor d of g we examine whether assigning a
    factor difference structure to d can produce a plausible
    factor pair.

    This is deliberately exploratory: it does not assume that
    p or q must divide g.
    """

    R = math.prod(mods)
    g = n - R

    abs_g = abs(g)

    divisors = []

    if abs_g > 0:
        root = math.isqrt(abs_g)

        for d in range(1, root + 1):
            if abs_g % d == 0:
                divisors.append(d)

                if d * d != abs_g:
                    divisors.append(abs_g // d)

    # Only small divisors are interesting for direct relations.
    divisors.sort()

    pset = set(factor_primes)

    direct_divisor_hits = []

    for d in divisors:

        # Test simple possibilities:
        #
        # p = d
        # q = n/d
        #
        # and
        #
        # p-q = d
        #
        # and
        #
        # p+q = d.
        if d in pset and n % d == 0:
            q = n // d

            if q in pset:
                direct_divisor_hits.append(
                    ("factor_divides_gap", d, d, q)
                )

        # Difference candidate:
        # q-p = d  =>
        #
        # p(p+d)=n
        #
        # p solves p^2 + d p - n = 0.
        D = d * d + 4 * n
        rootD = math.isqrt(D)

        if rootD * rootD == D:
            num = rootD - d

            if num > 0 and num % 2 == 0:
                p = num // 2
                q = p + d

                if p in pset and q in pset and p * q == n:
                    direct_divisor_hits.append(
                        ("factor_difference_is_gap_divisor", d, p, q)
                    )

    return {
        "gap": g,
        "divisor_count": len(divisors),
        "direct_hits": direct_divisor_hits,
        "small_divisors": divisors[:50],
    }


# ==================================================================================================
# MODULUS-OFFSET SEARCH
# ==================================================================================================

def modulus_offset_search(n, p, q, mods):
    """
    Look at the relationship between r_i and the factor residues.

    Define

        a_i = p mod r_i
        b_i = q mod r_i

    Then

        a_i*b_i == n mod r_i.

    We inspect several low-dimensional combinations involving:

        a_i
        b_i
        a_i+b_i
        a_i*b_i
        p//r_i
        q//r_i

    The point is to find quantities which are anomalously small or
    tightly related to g = n-R.
    """

    R = math.prod(mods)
    g = n - R

    rows = []

    for r in mods:

        a = p % r
        b = q % r

        kp = p // r
        kq = q // r

        rows.append(
            {
                "r": r,
                "a": a,
                "b": b,
                "a+b": a + b,
                "a*b": a * b,
                "kp": kp,
                "kq": kq,
                "kp+kq": kp + kq,
                "kp*kq": kp * kq,
                "r-(a+b)": r - a - b,
                "a*b-g": a * b - g,
            }
        )

    return rows


# ==================================================================================================
# MAIN
# ==================================================================================================

def run():

    t0 = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME GAP EXPLOITATION / R - N DEFECT EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                  = {ANCHORS}")
    print(f"modulus prime range      = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio              = {CLOSE_RATIO:.0%}")
    print(f"seed                     = {SEED:,}")
    print()

    # ==============================================================================================
    # PRIME POOLS
    # ==============================================================================================

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    all_primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in all_primes
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # ==============================================================================================
    # ANCHORS
    # ==============================================================================================

    anchors = generate_anchors(
        factor_primes,
        ANCHORS,
        rng
    )

    # ==============================================================================================
    # TRIPLES
    # ==============================================================================================

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    records = []

    for idx, (p, q, n) in enumerate(anchors, 1):

        mods = select_triple(
            n,
            modulus_primes
        )

        if mods is None:
            continue

        R = math.prod(mods)

        records.append(
            {
                "id": idx,
                "p": p,
                "q": q,
                "n": n,
                "mods": mods,
                "R": R,
                "gap": n - R,
            }
        )

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    print()
    print(f"usable anchors = {len(records)}")
    print()

    # ==============================================================================================
    # DIRECT GAP ANALYSIS
    # ==============================================================================================

    print("=" * 100)
    print("GAP STATISTICS")
    print("=" * 100)

    gaps = [r["gap"] for r in records]

    print(f"mean gap             = {sum(gaps) / len(gaps):,.3f}")
    print(f"median gap           = {sorted(gaps)[len(gaps)//2]:,}")
    print(f"minimum gap          = {min(gaps):,}")
    print(f"maximum gap          = {max(gaps):,}")
    print()

    # ==============================================================================================
    # TEST WHETHER GAP HAS SPECIAL ARITHMETIC STRUCTURE
    # ==============================================================================================

    gap_mod_classes = Counter()

    for rec in records:

        g = rec["gap"]
        mods = rec["mods"]

        residues = tuple(g % r for r in mods)

        gap_mod_classes[residues] += 1

    print("=" * 100)
    print("GAP RESIDUE STRUCTURE")
    print("=" * 100)

    print(f"unique gap residue vectors = {len(gap_mod_classes):,}")
    print(
        f"largest repeated vector   = "
        f"{max(gap_mod_classes.values())}"
    )

    print()

    # ==============================================================================================
    # MAIN SEARCH
    # ==============================================================================================

    total_divisors = 0
    total_direct_gap_hits = 0

    all_rows = []

    print("=" * 100)
    print("RUNNING GAP EXPLOITATION")
    print("=" * 100)

    for idx, rec in enumerate(records, 1):

        n = rec["n"]
        p = rec["p"]
        q = rec["q"]
        mods = rec["mods"]

        stats = gap_decomposition(
            n,
            p,
            q,
            mods
        )

        div = gap_divisibility_search(
            n,
            mods,
            factor_primes
        )

        total_divisors += div["divisor_count"]
        total_direct_gap_hits += len(div["direct_hits"])

        all_rows.append(
            {
                "id": rec["id"],
                "p": p,
                "q": q,
                "n": n,
                "mods": mods,
                "R": stats["R"],
                "gap": stats["g"],
                "ratio": stats["R"] / n,
                "sum": stats["s"],
                "diff": stats["d"],
                "gap_gcds": stats["gap_gcds"],
                "p_res": stats["p_res"],
                "q_res": stats["q_res"],
                "p_quot": stats["p_quot"],
                "q_quot": stats["q_quot"],
                "divisors": div["divisor_count"],
                "direct_hits": div["direct_hits"],
            }
        )

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{len(records)}")

    # ==============================================================================================
    # GAP / FACTOR RELATIONSHIPS
    # ==============================================================================================

    print()
    print("=" * 100)
    print("GAP VS FACTOR RELATIONSHIPS")
    print("=" * 100)

    def correlation(xs, ys):

        if len(xs) < 2:
            return 0.0

        xm = sum(xs) / len(xs)
        ym = sum(ys) / len(ys)

        num = sum(
            (x - xm) * (y - ym)
            for x, y in zip(xs, ys)
        )

        dx = math.sqrt(
            sum((x - xm) ** 2 for x in xs)
        )

        dy = math.sqrt(
            sum((y - ym) ** 2 for y in ys)
        )

        if dx == 0 or dy == 0:
            return 0.0

        return num / (dx * dy)

    gaps = [x["gap"] for x in all_rows]
    sums = [x["sum"] for x in all_rows]
    diffs = [x["diff"] for x in all_rows]
    products = [x["p"] * x["q"] for x in all_rows]
    pvals = [x["p"] for x in all_rows]
    qvals = [x["q"] for x in all_rows]

    print(
        f"corr(gap, p+q)       = "
        f"{correlation(gaps, sums):.6f}"
    )

    print(
        f"corr(gap, |p-q|)     = "
        f"{correlation(gaps, diffs):.6f}"
    )

    print(
        f"corr(gap, p)         = "
        f"{correlation(gaps, pvals):.6f}"
    )

    print(
        f"corr(gap, q)         = "
        f"{correlation(gaps, qvals):.6f}"
    )

    # ==============================================================================================
    # SPECIAL TEST:
    # CAN THE GAP PREDICT THE FACTOR DIFFERENCE?
    # ==============================================================================================

    print()
    print("=" * 100)
    print("GAP -> FACTOR DIFFERENCE TEST")
    print("=" * 100)

    exact_difference_hits = 0
    approximate_difference_hits = 0

    for row in all_rows:

        n = row["n"]
        g = row["gap"]
        d = row["diff"]

        # Try simple transformations.
        candidates = {
            g,
            math.isqrt(g),
            2 * math.isqrt(g),
            g // 2,
            g // 3,
            g // 4,
            g + 1,
            g - 1 if g > 0 else 0,
        }

        if d in candidates:
            exact_difference_hits += 1

        # Approximate scaling test:
        #
        # d^2 = (p-q)^2
        #
        # and
        #
        # n = ((p+q)^2 - d^2)/4.
        #
        # Check whether a low-degree polynomial in g gives d.
        for scale in range(1, 11):

            if abs(d * d - scale * g) <= 2 * scale:
                approximate_difference_hits += 1
                break

    print(
        f"simple exact hits           = "
        f"{exact_difference_hits}/{len(all_rows)}"
    )

    print(
        f"low-scale d²≈k·gap hits     = "
        f"{approximate_difference_hits}/{len(all_rows)}"
    )

    # ==============================================================================================
    # MODULUS-OFFSET STRUCTURE
    # ==============================================================================================

    print()
    print("=" * 100)
    print("MODULUS-OFFSET STRUCTURE")
    print("=" * 100)

    offset_counter = Counter()

    for row in all_rows:

        p = row["p"]
        q = row["q"]
        mods = row["mods"]

        rows = modulus_offset_search(
            row["n"],
            p,
            q,
            mods
        )

        for item in rows:

            signature = (
                item["a+b"],
                item["r-(a+b)"],
                item["kp+kq"],
            )

            offset_counter[signature] += 1

    print(
        f"unique offset signatures = "
        f"{len(offset_counter):,}"
    )

    print(
        f"largest repeated offset signature = "
        f"{max(offset_counter.values())}"
    )

    # ==============================================================================================
    # DIRECT GAP DIVISOR RESULTS
    # ==============================================================================================

    print()
    print("=" * 100)
    print("GAP DIVISOR SEARCH")
    print("=" * 100)

    print(
        f"total gap divisors examined = "
        f"{total_divisors:,}"
    )

    print(
        f"direct factor/divisor hits   = "
        f"{total_direct_gap_hits:,}"
    )

    # ==============================================================================================
    # MOST INTERESTING CASES
    # ==============================================================================================

    print()
    print("=" * 100)
    print("SMALLEST GAP CASES")
    print("=" * 100)

    smallest = sorted(
        all_rows,
        key=lambda x: x["gap"]
    )

    for row in smallest[:25]:

        p = row["p"]
        q = row["q"]
        r1, r2, r3 = row["mods"]

        print(
            f"n={row['n']:,} "
            f"p={p:,} q={q:,} "
            f"mods=({r1},{r2},{r3}) "
            f"R={row['R']:,} "
            f"gap={row['gap']:,} "
            f"R/n={row['ratio']:.12f} "
            f"|p-q|={row['diff']:,}"
        )

    # ==============================================================================================
    # RESIDUE DETAILS FOR SMALLEST GAPS
    # ==============================================================================================

    print()
    print("=" * 100)
    print("SMALLEST GAP RESIDUE DETAILS")
    print("=" * 100)

    for row in smallest[:15]:

        print(
            f"n={row['n']:,} "
            f"gap={row['gap']:,} "
            f"mods={row['mods']} "
            f"p_res={row['p_res']} "
            f"q_res={row['q_res']} "
            f"gap_mods={tuple(row['gap'] % r for r in row['mods'])} "
            f"gap_gcds={row['gap_gcds']}"
        )

    # ==============================================================================================
    # COMPLETE ANCHOR SUMMARY
    # ==============================================================================================

    print()
    print("=" * 100)
    print("ANCHOR SUMMARY")
    print("=" * 100)

    print(
        " ID        p        q       GAP       R/n       "
        "|p-q|      p+r       q+r"
    )

    print("-" * 100)

    for row in all_rows:

        p = row["p"]
        q = row["q"]
        r1, r2, r3 = row["mods"]

        print(
            f"{row['id']:3d} "
            f"{p:8,d} "
            f"{q:8,d} "
            f"{row['gap']:9,d} "
            f"{row['ratio']:.7f} "
            f"{row['diff']:9,d} "
            f"{p+r1:9,d} "
            f"{q+r3:9,d}"
        )

    # ==============================================================================================
    # FINAL EXPLOITATION TEST
    # ==============================================================================================

    print()
    print("=" * 100)
    print("FINAL EXPLOITATION TEST")
    print("=" * 100)

    print()
    print("Known:")
    print()
    print("    n = R + g")
    print()
    print("Therefore:")
    print()
    print("    p*q = R + g")
    print()
    print("and hence:")
    print()
    print("    p*q - R = g")
    print()

    print(
        "The experiment searched for direct arithmetic relationships "
        "between g and:"
    )

    print()
    print("    p")
    print("    q")
    print("    p+q")
    print("    |p-q|")
    print("    p mod r_i")
    print("    q mod r_i")
    print("    quotient coordinates")
    print("    divisors of g")
    print()

    print(
        "The decisive outcome is whether any of these quantities "
        "produces a substantially smaller search space than the "
        "8,363-prime baseline."
    )

    # ==============================================================================================
    # MOST IMPORTANT OBSERVATION
    # ==============================================================================================

    print()
    print("=" * 100)
    print("MOST IMPORTANT ALGEBRAIC OBSERVATION")
    print("=" * 100)

    print()
    print("Because R < n:")
    print()
    print("    n mod R = n - R = g")
    print()

    print("But also:")
    print()
    print("    p*q mod R = g")
    print()

    print(
        "So the entire factorization problem can be viewed as:"
    )

    print()
    print("    find p,q in the factor interval such that")
    print()
    print("        p*q ≡ g (mod R)")
    print()
    print("    AND")
    print()
    print("        p*q = R+g"
    )
    print()

    print(
        "The next step is therefore NOT another generic fingerprint "
        "experiment. It is to investigate whether the congruence"
    )

    print()
    print("    p*q ≡ g (mod r1*r2*r3)")
    print()
    print(
        "can be decomposed into a smaller bilinear search using the "
        "fact that r1,r2,r3 are close and R is just below n."
    )

    # ==============================================================================================
    # TIMING
    # ==============================================================================================

    elapsed = time.perf_counter() - t0

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime = {elapsed:.3f} seconds")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
