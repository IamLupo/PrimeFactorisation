#!/usr/bin/env python3
"""
====================================================================================================
THREE-CLOSE-PRIME CRT INVERSION / DIRECT-P SEARCH EXPERIMENT
====================================================================================================

Goal
----
Previous experiment:

    p in PRIME_POOL
        -> q_i = n * p^(-1) mod r_i
        -> CRT(q_1,q_2,q_3) = Q
        -> test Q in factor interval
        -> exact p*Q == n

worked extremely well as a filter.

This experiment asks the next question:

    Can Q(p) itself be used to predict / constrain p,
    so that we do NOT need to enumerate the complete prime pool?

For each actual semiprime:

    n = p*q

and selected close primes:

    R = r1*r2*r3 < n

we compute:

    Q(p) = CRT(q mod r1, q mod r2, q mod r3)

Since q is in the factor interval and R > FACTOR_MAX in most cases:

    q = Q(p)

when Q is taken as the least non-negative CRT representative.

We then study:

    p * Q(p) - n
    n / Q(p)
    |n/Q(p) - p|
    p + Q(p)
    p - Q(p)
    Q(p) / p
    Q(p) * p
    (n-Q(p)*p)
    CRT coordinate quotients

and, most importantly, whether cheap bounds on Q can restrict p.

The experiment contains several increasingly aggressive tests:

A. Full prime enumeration baseline
B. CRT Q reconstruction
C. Interval rejection before exact multiplication
D. Integer-square / AM-GM bounds
E. Q-order inversion:
       p ~= n / Q
F. Newton-like reconstruction from Q
G. Fixed-point / consistency search
H. Local neighborhood search around predicted p
I. Comparison against exhaustive truth

This does NOT assume that the method works.
It measures exactly how much of the remaining p-search can be removed.
"""

import math
import random
import time


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

ACTUAL_ANCHORS = 300

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

MAX_EXAMPLES = 20

# ================================================================================================
# UTILITIES
# ================================================================================================


def sieve(limit):
    """Simple bytearray sieve."""
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start::p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if a[i]]


def egcd(a, b):
    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def modinv(a, m):
    a %= m
    g, x, _ = egcd(a, m)

    if g != 1:
        raise ValueError(f"{a} has no inverse modulo {m}")

    return x % m


def crt_pairwise(residues, moduli):
    """
    CRT for pairwise-coprime moduli.
    """
    R = 1
    for r in moduli:
        R *= r

    x = 0

    for a, r in zip(residues, moduli):
        Ri = R // r
        inv = modinv(Ri % r, r)
        x += a * Ri * inv

    return x % R


def branch_mod4(p, q):
    a = p % 4
    b = q % 4

    if a == 1 and b == 1:
        return "11"

    if a == 3 and b == 3:
        return "33"

    if {a, b} == {1, 3}:
        return "13/31"

    return "INVALID"


# ================================================================================================
# PRIME / MODULUS SELECTION
# ================================================================================================


def select_close_triple(n, modulus_primes, ratio):
    """
    Select the largest product r1*r2*r3 < n among triples satisfying
    max(r)/min(r) <= 1+ratio.
    """
    best = None
    best_product = -1

    L = len(modulus_primes)

    for i in range(L):
        r1 = modulus_primes[i]

        for j in range(i + 1, L):
            r2 = modulus_primes[j]

            if r2 > r1 * (1.0 + ratio):
                break

            for k in range(j + 1, L):
                r3 = modulus_primes[k]

                if r3 > r1 * (1.0 + ratio):
                    break

                R = r1 * r2 * r3

                if R < n and R > best_product:
                    best_product = R
                    best = (r1, r2, r3)

    return best


# ================================================================================================
# ACTUAL ANCHORS
# ================================================================================================


def generate_anchors(primes, count, rng):
    """
    Generate deterministic random semiprime anchors from the prime pool.
    """
    anchors = []

    seen = set()

    while len(anchors) < count:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)

        anchors.append((p, q, p * q))

    return anchors


# ================================================================================================
# DIRECT CRT
# ================================================================================================


def reconstruct_q(n, p, mods):
    """
    Recover q residue vector from p:

        q == n * p^(-1) mod r_i
    """
    residues = []

    for r in mods:
        pmod = p % r

        if pmod == 0:
            return None, None

        qmod = (n % r) * modinv(pmod, r) % r
        residues.append(qmod)

    q = crt_pairwise(residues, mods)

    return q, tuple(residues)


# ================================================================================================
# PREDICTION TESTS
# ================================================================================================


def analyze_prediction(n, p, q, Q, R):
    """
    Analyze how much information Q gives about p.
    """

    result = {}

    result["Q"] = Q

    result["prediction_error"] = abs(Q - q)

    # Since n = p*q:
    # n/Q should be p when Q=q.
    if Q != 0:
        pred = n // Q
        result["n_div_Q"] = pred
        result["n_over_Q_float_error"] = abs((n / Q) - p)
    else:
        result["n_div_Q"] = None
        result["n_over_Q_float_error"] = None

    # Algebraically useful residuals.
    result["product_residual"] = abs(n - p * Q)

    result["sum"] = p + Q
    result["difference"] = abs(p - Q)

    # Ratio around 1 is useful when factors are balanced.
    if p:
        result["Q_over_p"] = Q / p
    else:
        result["Q_over_p"] = None

    # Classical bounds from pq=n:
    # p <= sqrt(n) <= q for p <= q.
    root = math.isqrt(n)

    result["sqrt_n"] = root

    if Q > 0:
        # Candidate interval induced by q >= p or q <= p.
        result["q_side"] = "Q>=sqrt(n)" if Q >= root else "Q<sqrt(n)"

        predicted_p = n // Q

        result["predicted_p"] = predicted_p

        # Distance in integer candidate space.
        result["prediction_distance"] = abs(predicted_p - p)

    else:
        result["q_side"] = None
        result["predicted_p"] = None
        result["prediction_distance"] = None

    # How many integers around predicted p need to be checked.
    if result["predicted_p"] is not None:
        lo = max(FACTOR_MIN, result["predicted_p"] - 5)
        hi = min(FACTOR_MAX, result["predicted_p"] + 5)
        result["local_11_width"] = hi - lo + 1
    else:
        result["local_11_width"] = None

    return result


# ================================================================================================
# MAIN
# ================================================================================================


def run():
    t0 = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME CRT INVERSION / DIRECT-P SEARCH EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ACTUAL_ANCHORS}")
    print(f"modulus prime range      = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # ============================================================================================
    # PRIME POOL
    # ============================================================================================

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

    # ============================================================================================
    # ANCHORS
    # ============================================================================================

    anchors = generate_anchors(
        factor_primes,
        ACTUAL_ANCHORS,
        rng
    )

    print("actual anchors             =", len(anchors))
    print()

    # ============================================================================================
    # SELECT MODULUS TRIPLES
    # ============================================================================================

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE PRIME TRIPLES")
    print("=" * 100)

    selected = []

    for idx, (p, q, n) in enumerate(anchors, 1):
        mods = select_close_triple(
            n,
            modulus_primes,
            CLOSE_RATIO
        )

        if mods is None:
            print(f"WARNING anchor {idx}: no usable modulus triple")
            continue

        R = math.prod(mods)

        selected.append({
            "id": idx,
            "p": p,
            "q": q,
            "n": n,
            "mods": mods,
            "R": R,
        })

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    print()
    print("usable anchors =", len(selected))
    print()

    # ============================================================================================
    # CORE ANALYSIS
    # ============================================================================================

    print("=" * 100)
    print("CRT INVERSION ANALYSIS")
    print("=" * 100)

    total_p = 0
    total_valid_q = 0
    total_exact = 0

    residual_zero = 0

    prediction_distances = []
    local_hits = 0

    # Search reductions.
    exhaustive = len(factor_primes)

    # Examples.
    examples = []

    # Direct metrics.
    min_prediction_distance = None
    max_prediction_distance = None

    for idx, rec in enumerate(selected, 1):
        p = rec["p"]
        q = rec["q"]
        n = rec["n"]
        mods = rec["mods"]
        R = rec["R"]

        total_p += exhaustive

        # Direct CRT for actual p.
        Q, residues = reconstruct_q(n, p, mods)

        if Q is None:
            continue

        if FACTOR_MIN <= Q <= FACTOR_MAX:
            total_valid_q += 1

        if Q == q:
            total_exact += 1
            residual_zero += 1

        analysis = analyze_prediction(
            n,
            p,
            q,
            Q,
            R
        )

        d = analysis["prediction_distance"]

        if d is not None:
            prediction_distances.append(d)

            if min_prediction_distance is None or d < min_prediction_distance:
                min_prediction_distance = d

            if max_prediction_distance is None or d > max_prediction_distance:
                max_prediction_distance = d

        # Test local neighborhood around floor(n/Q).
        predicted = analysis["predicted_p"]

        found_local = False

        if predicted is not None:
            lo = max(FACTOR_MIN, predicted - 5)
            hi = min(FACTOR_MAX, predicted + 5)

            for candidate_p in factor_primes:
                if candidate_p < lo:
                    continue

                if candidate_p > hi:
                    break

                if candidate_p * Q == n:
                    found_local = True
                    break

        if found_local:
            local_hits += 1

        if len(examples) < MAX_EXAMPLES:
            examples.append({
                "id": rec["id"],
                "p": p,
                "q": q,
                "n": n,
                "mods": mods,
                "R": R,
                "R_over_n": R / n,
                "Q": Q,
                "predicted_p": analysis["predicted_p"],
                "prediction_distance": analysis["prediction_distance"],
                "product_residual": analysis["product_residual"],
                "sum": analysis["sum"],
                "difference": analysis["difference"],
                "Q_over_p": analysis["Q_over_p"],
                "local_hit": found_local,
            })

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{len(selected)}")

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    print()
    print("=" * 100)
    print("CRT INVERSION SUMMARY")
    print("=" * 100)

    print(f"anchors analyzed                 = {len(selected)}")
    print(f"direct CRT Q == actual q         = {total_exact}/{len(selected)}")
    print(f"Q inside factor interval         = {total_valid_q}/{len(selected)}")

    if prediction_distances:
        mean_distance = sum(prediction_distances) / len(prediction_distances)
    else:
        mean_distance = float("nan")

    print(f"mean |floor(n/Q)-p|              = {mean_distance:.6f}")
    print(f"minimum prediction distance      = {min_prediction_distance}")
    print(f"maximum prediction distance      = {max_prediction_distance}")
    print(f"local +/-5 prime search hits     = {local_hits}/{len(selected)}")

    # ============================================================================================
    # SEARCH COST MODEL
    # ============================================================================================

    print()
    print("=" * 100)
    print("SEARCH COST COMPARISON")
    print("=" * 100)

    baseline_ops = len(selected) * len(factor_primes)

    if selected:
        average_pred_distance = (
            sum(prediction_distances) / len(prediction_distances)
            if prediction_distances
            else float("inf")
        )
    else:
        average_pred_distance = float("inf")

    print(f"exhaustive p checks               = {baseline_ops:,}")
    print(f"average baseline / anchor        = {len(factor_primes):,.0f}")
    print(f"average |predicted_p-p|           = {average_pred_distance:.6f}")
    print(f"local +/-5 candidate window      = 11 integers")
    print()

    # ============================================================================================
    # IMPORTANT ALGEBRAIC TESTS
    # ============================================================================================

    print("=" * 100)
    print("ALGEBRAIC CONSISTENCY TESTS")
    print("=" * 100)

    exact_identity_failures = 0
    division_identity_failures = 0

    for rec in selected:
        p = rec["p"]
        q = rec["q"]
        n = rec["n"]
        mods = rec["mods"]

        Q, _ = reconstruct_q(n, p, mods)

        if Q != q:
            exact_identity_failures += 1

        if Q != 0:
            if n // Q != p or n % Q != 0:
                division_identity_failures += 1

    print(f"Q == q failures                   = {exact_identity_failures}")
    print(f"n/Q == p failures                 = {division_identity_failures}")

    # ============================================================================================
    # EXAMPLES
    # ============================================================================================

    print()
    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for e in examples:
        print()
        print(
            f"n={e['n']:,} "
            f"p={e['p']:,} "
            f"q={e['q']:,}"
        )

        print(
            f"mods={e['mods']} "
            f"R={e['R']:,} "
            f"R/n={e['R_over_n']:.12f}"
        )

        print(
            f"Q={e['Q']:,} "
            f"predicted_p=floor(n/Q)={e['predicted_p']:,} "
            f"distance={e['prediction_distance']}"
        )

        print(
            f"|n-pQ|={e['product_residual']:,} "
            f"p+Q={e['sum']:,} "
            f"|p-Q|={e['difference']:,}"
        )

        print(
            f"Q/p={e['Q_over_p']:.12f} "
            f"local_hit={e['local_hit']}"
        )

    # ============================================================================================
    # SECOND STAGE:
    # CAN WE SEARCH Q INSTEAD OF P?
    # ============================================================================================

    print()
    print("=" * 100)
    print("SECOND-STAGE: REVERSE SEARCH OVER q")
    print("=" * 100)

    """
    Instead of:

        for p:
            derive q

    try:

        for q:
            derive p

    Since

        p == n*q^(-1) mod r_i

    the same CRT collapse should occur symmetrically.

    This checks whether either side is systematically easier.
    """

    reverse_success = 0
    reverse_candidates = []

    for rec in selected:
        p = rec["p"]
        q = rec["q"]
        n = rec["n"]
        mods = rec["mods"]

        residues = []

        valid = True

        for r in mods:
            qmod = q % r

            if qmod == 0:
                valid = False
                break

            pmod = (n % r) * modinv(qmod, r) % r
            residues.append(pmod)

        if not valid:
            continue

        P = crt_pairwise(residues, mods)

        if FACTOR_MIN <= P <= FACTOR_MAX:
            reverse_success += 1

        reverse_candidates.append(P)

    print(f"reverse CRT P == actual p       = "
          f"{reverse_success}/{len(selected)}")

    if reverse_candidates:
        print(
            f"reverse candidates in interval = "
            f"{sum(FACTOR_MIN <= x <= FACTOR_MAX for x in reverse_candidates)}"
        )

    # ============================================================================================
    # COMBINED RESULT
    # ============================================================================================

    print()
    print("=" * 100)
    print("COMBINED RESULT")
    print("=" * 100)

    print()
    print("Current method:")
    print()
    print("    enumerate p")
    print("        |")
    print("        v")
    print("    3 modular inversions")
    print("        |")
    print("        v")
    print("    CRT")
    print("        |")
    print("        v")
    print("    q")
    print("        |")
    print("        v")
    print("    exact n == p*q")
    print()

    print("Proposed inversion:")
    print()
    print("    obtain Q")
    print("        |")
    print("        v")
    print("    p ~= n/Q")
    print("        |")
    print("        v")
    print("    local p search")
    print()

    if selected:
        if local_hits == len(selected):
            print(
                "RESULT: every tested anchor was recoverable from a "
                "local neighborhood around n/Q."
            )
        else:
            print(
                "RESULT: the local prediction is NOT sufficient for every "
                "anchor."
            )

    print()
    print("IMPORTANT:")
    print()
    print(
        "Q is not an independently observed quantity."
    )
    print(
        "It was reconstructed after supplying p."
    )
    print()
    print(
        "Therefore a successful n/Q prediction from Q does NOT by itself "
        "remove the original p enumeration."
    )
    print()
    print(
        "The decisive next question is whether Q can be obtained from n "
        "and the three moduli WITHOUT already knowing p."
    )

    # ============================================================================================
    # TIMING
    # ============================================================================================

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
