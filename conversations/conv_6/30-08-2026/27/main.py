#!/usr/bin/env python3

import math
import random
import time
from bisect import bisect_right


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

ANCHORS = 300
CLOSE_RATIO = 0.20
SEED = 1_511_464_998

PROGRESS_EVERY = 25

# If True, print only the important summary.
COMPACT = False


# =============================================================================
# BASIC UTILITIES
# =============================================================================

def banner(title):
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def ceil_div(a, b):
    return -((-a) // b)


def isqrt_exact(n):
    if n < 0:
        return None

    r = math.isqrt(n)
    if r * r == n:
        return r

    return None


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(limit):
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    stop = math.isqrt(limit)

    for p in range(2, stop + 1):
        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start::p] = b"\x00" * count

    return [i for i, v in enumerate(a) if v]


# =============================================================================
# CRT
# =============================================================================

def crt_two(a1, m1, a2, m2):
    """
    Solve:
        x = a1 (mod m1)
        x = a2 (mod m2)

    m1 and m2 must be coprime.
    """

    k = ((a2 - a1) * pow(m1, -1, m2)) % m2
    x = a1 + m1 * k

    return x % (m1 * m2)


def crt_coefficient(c1, c2, m1, m2):
    """
    CRT-combine one coefficient modulo m1 and m2.
    """

    return crt_two(c1 % m1, m1, c2 % m2, m2)


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def build_anchors(factor_primes, count, seed):
    """
    Deterministic random factor-prime pairs.

    The pairs are chosen independently, matching the repeated-p examples
    seen in the previous experiments.
    """

    rng = random.Random(seed)

    anchors = []

    for _ in range(count):
        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        while q == p:
            q = rng.choice(factor_primes)

        if p > q:
            p, q = q, p

        n = p * q

        anchors.append((p, q, n))

    return anchors


# =============================================================================
# CLOSE MAXIMUM-PRODUCT MODULUS TRIPLE
# =============================================================================

def select_best_triple(n, modulus_primes):
    """
    Select:

        r1 < r2 < r3

    subject to:

        r3 / r1 <= 1 + CLOSE_RATIO
        r1*r2*r3 < n

    and maximize:

        r1*r2*r3

    This is the same structural condition used throughout the experiments:
    make R as close to n from below as possible while keeping the three
    primes close together.
    """

    best = None
    best_product = -1

    L = len(modulus_primes)

    ratio_limit = 1.0 + CLOSE_RATIO

    for i in range(L - 2):

        r1 = modulus_primes[i]

        # r3/r1 <= ratio_limit
        max_r3 = int(r1 * ratio_limit)

        j_limit = bisect_right(modulus_primes, max_r3)

        if j_limit <= i + 1:
            continue

        for j in range(i + 1, j_limit - 1):

            r2 = modulus_primes[j]

            # Need r3 > r2.
            if r1 * r2 * r2 >= n:
                break

            max_r3_by_product = (n - 1) // (r1 * r2)

            max_r3 = min(max_r3, max_r3_by_product)

            k = bisect_right(modulus_primes, max_r3) - 1

            if k <= j:
                continue

            r3 = modulus_primes[k]

            R = r1 * r2 * r3

            if R >= n:
                continue

            if R > best_product:
                best_product = R
                best = (r1, r2, r3)

    return best


# =============================================================================
# GAP-BILINEAR COEFFICIENT CONSTRUCTION
# =============================================================================

def build_bilinear_coefficients(a, b, g, r1, r2, r3):
    """
    Define:

        d2 = r2-r1
        d3 = r3-r1

    Since

        p = a + k*r1
        q = b + l*r1

    we have modulo ri:

        p = a-k*di
        q = b-l*di

    and therefore:

        (a-k*di)(b-l*di) - g == 0 (mod ri)

    Expanding:

        di^2*k*l
        - b*di*k
        - a*di*l
        + (a*b-g)

    For i=2 and i=3 we CRT-combine every coefficient into one equation:

        A*k*l + B*k + C*l + D == 0 (mod r2*r3)
    """

    d2 = r2 - r1
    d3 = r3 - r1

    # r2 equation
    A2 = d2 * d2
    B2 = -b * d2
    C2 = -a * d2
    D2 = a * b - g

    # r3 equation
    A3 = d3 * d3
    B3 = -b * d3
    C3 = -a * d3
    D3 = a * b - g

    R23 = r2 * r3

    A = crt_coefficient(A2, A3, r2, r3)
    B = crt_coefficient(B2, B3, r2, r3)
    C = crt_coefficient(C2, C3, r2, r3)
    D = crt_coefficient(D2, D3, r2, r3)

    return A, B, C, D, R23


# =============================================================================
# HYPERBOLA TRANSFORMATION
# =============================================================================

def hyperbola_rhs(A, B, C, D, t, R23):
    """
    From

        A*k*l + B*k + C*l + D = t*R23

    multiply by A and complete the product:

        (A*k+C)(A*l+B)
            =
        A*t*R23 + B*C - A*D
    """

    return A * t * R23 + B * C - A * D


# =============================================================================
# SEARCH RANGE
# =============================================================================

def quotient_range(residue, modulus, low, high):
    """
    Find all integer k for which:

        low <= residue + k*modulus <= high
    """

    k_min = ceil_div(low - residue, modulus)
    k_max = (high - residue) // modulus

    return k_min, k_max


# =============================================================================
# HYPERBOLA SEARCH FOR ONE RESIDUE STATE
# =============================================================================

def solve_a_state(
    a,
    b,
    n,
    g,
    r1,
    r2,
    r3,
):
    """
    Search one a = p mod r1 state.

    Returns statistics and recovered pairs.
    """

    A, B, C, D, R23 = build_bilinear_coefficients(
        a, b, g, r1, r2, r3
    )

    if A == 0:
        return {
            "valid": False,
            "reason": "A=0",
        }

    k_min, k_max = quotient_range(
        a,
        r1,
        FACTOR_MIN,
        FACTOR_MAX,
    )

    l_min, l_max = quotient_range(
        b,
        r1,
        FACTOR_MIN,
        FACTOR_MAX,
    )

    if k_min > k_max or l_min > l_max:
        return {
            "valid": False,
            "reason": "empty quotient range",
        }

    # ---------------------------------------------------------
    # The canonical CRT equation is:
    #
    #   F = Akl + Bk + Cl + D
    #
    # All coefficients are non-negative residues modulo R23.
    # Therefore F is monotonic in k and l.
    # ---------------------------------------------------------

    F_min = (
        A * k_min * l_min
        + B * k_min
        + C * l_min
        + D
    )

    F_max = (
        A * k_max * l_max
        + B * k_max
        + C * l_max
        + D
    )

    t_min = ceil_div(F_min, R23)
    t_max = F_max // R23

    if t_min > t_max:
        return {
            "valid": True,
            "k_count": k_max - k_min + 1,
            "l_count": l_max - l_min + 1,
            "t_count": 0,
            "t_survivors": 0,
            "exact_candidates": 0,
            "solutions": [],
        }

    t_count = t_max - t_min + 1

    # Product-form intervals.
    x_min = A * k_min + C
    x_max = A * k_max + C

    y_min = A * l_min + B
    y_max = A * l_max + B

    product_min = x_min * y_min
    product_max = x_max * y_max

    t_survivors = 0
    exact_candidates = 0
    solutions = []

    # ---------------------------------------------------------
    # Search quotient parameter t.
    #
    # Only t values whose hyperbola can physically intersect
    # the (k,l) rectangle survive the interval test.
    # ---------------------------------------------------------

    for t in range(t_min, t_max + 1):

        N = hyperbola_rhs(
            A,
            B,
            C,
            D,
            t,
            R23,
        )

        if N <= 0:
            continue

        if N < product_min or N > product_max:
            continue

        t_survivors += 1

        # -----------------------------------------------------
        # Exact lattice recovery.
        #
        # We use the shorter quotient dimension.
        # -----------------------------------------------------

        k_span = k_max - k_min + 1
        l_span = l_max - l_min + 1

        if k_span <= l_span:

            for k in range(k_min, k_max + 1):

                x = A * k + C

                if x == 0:
                    continue

                if N % x != 0:
                    continue

                y = N // x

                if (y - B) % A != 0:
                    continue

                l = (y - B) // A

                if not (l_min <= l <= l_max):
                    continue

                p = a + k * r1
                q = b + l * r1

                if not (
                    FACTOR_MIN <= p <= FACTOR_MAX
                    and FACTOR_MIN <= q <= FACTOR_MAX
                ):
                    continue

                exact_candidates += 1

                if p * q == n:
                    pair = tuple(sorted((p, q)))

                    if pair not in solutions:
                        solutions.append(pair)

        else:

            for l in range(l_min, l_max + 1):

                y = A * l + B

                if y == 0:
                    continue

                if N % y != 0:
                    continue

                x = N // y

                if (x - C) % A != 0:
                    continue

                k = (x - C) // A

                if not (k_min <= k <= k_max):
                    continue

                p = a + k * r1
                q = b + l * r1

                if not (
                    FACTOR_MIN <= p <= FACTOR_MAX
                    and FACTOR_MIN <= q <= FACTOR_MAX
                ):
                    continue

                exact_candidates += 1

                if p * q == n:
                    pair = tuple(sorted((p, q)))

                    if pair not in solutions:
                        solutions.append(pair)

    return {
        "valid": True,
        "k_count": k_max - k_min + 1,
        "l_count": l_max - l_min + 1,
        "t_count": t_count,
        "t_survivors": t_survivors,
        "exact_candidates": exact_candidates,
        "solutions": solutions,
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "R23": R23,
        "k_min": k_min,
        "k_max": k_max,
        "l_min": l_min,
        "l_max": l_max,
        "t_min": t_min,
        "t_max": t_max,
    }


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def run():

    start_total = time.perf_counter()

    banner("THREE-CLOSE-PRIME GAP-RESULTANT / HYPERBOLA ELIMINATION EXPERIMENT")

    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                  = {ANCHORS}")
    print(f"modulus prime range      = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio              = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")

    # -------------------------------------------------------------------------
    # PRIME POOLS
    # -------------------------------------------------------------------------

    banner("BUILDING PRIME POOLS")

    factor_primes = sieve(FACTOR_MAX)
    factor_primes = [
        p for p in factor_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = sieve(MOD_MAX)
    modulus_primes = [
        p for p in modulus_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    # -------------------------------------------------------------------------
    # ANCHORS
    # -------------------------------------------------------------------------

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    # -------------------------------------------------------------------------
    # SELECT TRIPLES
    # -------------------------------------------------------------------------

    banner("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")

    usable = []

    for idx, (p, q, n) in enumerate(anchors, 1):

        triple = select_best_triple(
            n,
            modulus_primes,
        )

        if triple is None:
            continue

        r1, r2, r3 = triple

        R = r1 * r2 * r3

        usable.append({
            "id": idx,
            "p": p,
            "q": q,
            "n": n,
            "r1": r1,
            "r2": r2,
            "r3": r3,
            "R": R,
            "gap": n - R,
        })

        if idx % PROGRESS_EVERY == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    print()
    print(f"usable anchors = {len(usable)}")

    # -------------------------------------------------------------------------
    # MAIN SEARCH
    # -------------------------------------------------------------------------

    banner("RUNNING GAP-RESULTANT / HYPERBOLA SEARCH")

    total_a_states = 0
    total_original_ak_states = 0
    total_t_states = 0
    total_t_survivors = 0
    total_exact_candidates = 0

    recovery_count = 0
    unique_recovery_count = 0

    t_ratio_values = []
    survivor_ratio_values = []

    strongest = []
    weakest = []

    rows = []

    search_start = time.perf_counter()

    for idx, anchor in enumerate(usable, 1):

        p_actual = anchor["p"]
        q_actual = anchor["q"]
        n = anchor["n"]

        r1 = anchor["r1"]
        r2 = anchor["r2"]
        r3 = anchor["r3"]

        R = anchor["R"]
        g = anchor["gap"]

        # p mod r1
        actual_a = p_actual % r1

        if actual_a == 0:
            continue

        # From:
        #
        #   p*q == g (mod r1)
        #
        # and:
        #
        #   a*b == g (mod r1)
        #
        # recover b.
        actual_b = (g * pow(actual_a, -1, r1)) % r1

        # ---------------------------------------------------------------------
        # Enumerate all possible residues a.
        # ---------------------------------------------------------------------

        anchor_a_states = 0
        anchor_ak_states = 0
        anchor_t_states = 0
        anchor_t_survivors = 0
        anchor_exact = 0

        recovered_pairs = set()

        actual_state = None

        for a in range(1, r1):

            if math.gcd(a, r1) != 1:
                continue

            b = (g * pow(a, -1, r1)) % r1

            k_min, k_max = quotient_range(
                a,
                r1,
                FACTOR_MIN,
                FACTOR_MAX,
            )

            l_min, l_max = quotient_range(
                b,
                r1,
                FACTOR_MIN,
                FACTOR_MAX,
            )

            if k_min > k_max or l_min > l_max:
                continue

            anchor_a_states += 1

            k_count = k_max - k_min + 1
            l_count = l_max - l_min + 1

            anchor_ak_states += k_count

            result = solve_a_state(
                a,
                b,
                n,
                g,
                r1,
                r2,
                r3,
            )

            if not result["valid"]:
                continue

            anchor_t_states += result["t_count"]
            anchor_t_survivors += result["t_survivors"]
            anchor_exact += result["exact_candidates"]

            for pair in result["solutions"]:
                recovered_pairs.add(pair)

            if a == actual_a and b == actual_b:
                actual_state = result

        # ---------------------------------------------------------------------
        # Recovery checks
        # ---------------------------------------------------------------------

        actual_pair = tuple(sorted((p_actual, q_actual)))

        recovered = actual_pair in recovered_pairs

        if recovered:
            recovery_count += 1

        if len(recovered_pairs) == 1 and actual_pair in recovered_pairs:
            unique_recovery_count += 1

        # ---------------------------------------------------------------------
        # Statistics
        # ---------------------------------------------------------------------

        t_ratio = (
            anchor_t_states / anchor_ak_states
            if anchor_ak_states
            else 0.0
        )

        survivor_ratio = (
            anchor_t_survivors / anchor_t_states
            if anchor_t_states
            else 0.0
        )

        t_ratio_values.append(t_ratio)
        survivor_ratio_values.append(survivor_ratio)

        total_a_states += anchor_a_states
        total_original_ak_states += anchor_ak_states
        total_t_states += anchor_t_states
        total_t_survivors += anchor_t_survivors
        total_exact_candidates += anchor_exact

        rows.append({
            "id": anchor["id"],
            "p": p_actual,
            "q": q_actual,
            "n": n,
            "r1": r1,
            "r2": r2,
            "r3": r3,
            "R": R,
            "gap": g,
            "R_over_n": R / n,
            "a_states": anchor_a_states,
            "ak_states": anchor_ak_states,
            "t_states": anchor_t_states,
            "t_survivors": anchor_t_survivors,
            "exact": anchor_exact,
            "recovered": recovered,
            "unique": (
                len(recovered_pairs) == 1
                and actual_pair in recovered_pairs
            ),
        })

        if idx % PROGRESS_EVERY == 0:
            print(f"anchor {idx:3d}/{len(usable)}")

    search_time = time.perf_counter() - search_start

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------

    banner("SUMMARY")

    N = len(rows)

    if N == 0:
        print("No usable anchors.")
        return

    avg_a = total_a_states / N
    avg_ak = total_original_ak_states / N
    avg_t = total_t_states / N
    avg_survivors = total_t_survivors / N
    avg_exact = total_exact_candidates / N

    print(f"anchors analyzed               = {N:,}")
    print()
    print(f"average a states               = {avg_a:,.3f}")
    print(f"average (a,k) states           = {avg_ak:,.3f}")
    print(f"average t states               = {avg_t:,.3f}")
    print(f"average surviving t states     = {avg_survivors:,.3f}")
    print(f"average exact candidates       = {avg_exact:,.3f}")

    # -------------------------------------------------------------------------
    # REDUCTION
    # -------------------------------------------------------------------------

    banner("STATE-SPACE REDUCTION")

    baseline = len(factor_primes)

    print(f"factor-prime baseline          = {baseline:,}")
    print(f"average (a,k) states           = {avg_ak:,.3f}")
    print(f"average t states               = {avg_t:,.3f}")
    print(f"average surviving t states     = {avg_survivors:,.3f}")

    if avg_ak:
        print()
        print(
            "t vs (a,k) reduction           = "
            f"{(1.0 - avg_t / avg_ak) * 100:.6f}%"
        )

    if avg_t:
        print(
            "hyperbola interval reduction   = "
            f"{(1.0 - avg_survivors / avg_t) * 100:.6f}%"
        )

    if baseline:
        print()
        print(
            "t / prime baseline             = "
            f"{avg_t / baseline:.6f}"
        )

        print(
            "surviving t / prime baseline   = "
            f"{avg_survivors / baseline:.6f}"
        )

    # -------------------------------------------------------------------------
    # RECOVERY
    # -------------------------------------------------------------------------

    banner("RECOVERY")

    print(
        f"correct factor pair recovered  = "
        f"{recovery_count}/{N}"
    )

    print(
        f"unique recovered pair          = "
        f"{unique_recovery_count}/{N}"
    )

    print(
        f"recovery rate                  = "
        f"{100.0 * recovery_count / N:.4f}%"
    )

    # -------------------------------------------------------------------------
    # RATIO DISTRIBUTION
    # -------------------------------------------------------------------------

    banner("QUOTIENT-PARAMETER DISTRIBUTION")

    print(
        f"mean t/(a,k)                   = "
        f"{sum(t_ratio_values) / len(t_ratio_values):.8f}"
    )

    print(
        f"minimum t/(a,k)                = "
        f"{min(t_ratio_values):.8f}"
    )

    print(
        f"maximum t/(a,k)                = "
        f"{max(t_ratio_values):.8f}"
    )

    print()
    print(
        f"mean surviving-t fraction      = "
        f"{sum(survivor_ratio_values) / len(survivor_ratio_values):.8f}"
    )

    # -------------------------------------------------------------------------
    # BEST/WORST ANCHORS
    # -------------------------------------------------------------------------

    rows_sorted = sorted(
        rows,
        key=lambda x: x["t_survivors"]
    )

    banner("STRONGEST HYPERBOLA COLLAPSES")

    for row in rows_sorted[:20]:

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods=({row['r1']},{row['r2']},{row['r3']}) "
            f"R/n={row['R_over_n']:.10f} "
            f"gap={row['gap']:,} "
            f"(a,k)={row['ak_states']:,} "
            f"t={row['t_states']:,} "
            f"survive={row['t_survivors']:,} "
            f"exact={row['exact']:,} "
            f"recovered={row['recovered']}"
        )

    banner("WEAKEST HYPERBOLA COLLAPSES")

    for row in rows_sorted[-20:]:

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods=({row['r1']},{row['r2']},{row['r3']}) "
            f"R/n={row['R_over_n']:.10f} "
            f"gap={row['gap']:,} "
            f"(a,k)={row['ak_states']:,} "
            f"t={row['t_states']:,} "
            f"survive={row['t_survivors']:,} "
            f"exact={row['exact']:,} "
            f"recovered={row['recovered']}"
        )

    # -------------------------------------------------------------------------
    # SELECTED ANCHOR TABLE
    # -------------------------------------------------------------------------

    banner("ANCHOR RESULTS")

    print(
        f"{'ID':>4} "
        f"{'p':>8} "
        f"{'q':>8} "
        f"{'r1':>5} "
        f"{'r2':>5} "
        f"{'r3':>5} "
        f"{'R/n':>10} "
        f"{'A,K':>8} "
        f"{'T':>8} "
        f"{'SURV':>7} "
        f"{'EXACT':>6}"
    )

    print("-" * 100)

    for row in rows[:100]:

        print(
            f"{row['id']:4d} "
            f"{row['p']:8,d} "
            f"{row['q']:8,d} "
            f"{row['r1']:5d} "
            f"{row['r2']:5d} "
            f"{row['r3']:5d} "
            f"{row['R_over_n']:10.7f} "
            f"{row['ak_states']:8,d} "
            f"{row['t_states']:8,d} "
            f"{row['t_survivors']:7,d} "
            f"{row['exact']:6,d}"
        )

    # -------------------------------------------------------------------------
    # FINAL INTERPRETATION
    # -------------------------------------------------------------------------

    banner("INTERPRETATION")

    print(
        """
The previous bilinear formulation searched the quotient lattice:

    (a,k,l)

with the second and third close moduli imposing bilinear congruences.

Here those two congruences are CRT-combined into:

    A*k*l + B*k + C*l + D == 0 (mod r2*r3)

and then written as:

    A*k*l + B*k + C*l + D = t*(r2*r3).

Multiplying by A gives the hyperbola form:

    (A*k + C)(A*l + B)
        =
    A*t*(r2*r3) + B*C - A*D.

Therefore the experiment asks whether the small modulus gaps allow us
to replace the two-dimensional quotient search by a substantially
smaller search over t.

The important measurements are:

    average (a,k) states
                vs
    average t states
                vs
    surviving t states.

A genuinely interesting result would be:

    t states << (a,k) states

AND

    surviving t states << t states

while still recovering every factor pair.

If t is approximately the same size as the original quotient lattice,
then the transformation is only a coordinate change.

If surviving t states become very small, the next question becomes
whether the hyperbola equation itself can be solved without enumerating
the remaining quotient coordinate.

The experiment deliberately separates those two effects.
"""
    )

    # -------------------------------------------------------------------------
    # TIMING
    # -------------------------------------------------------------------------

    total_time = time.perf_counter() - start_total

    banner("TIMING")

    print(f"search time                   = {search_time:.3f} seconds")
    print(f"total runtime                 = {total_time:.3f} seconds")

    banner("EXPERIMENT COMPLETE")


if __name__ == "__main__":
    run()
