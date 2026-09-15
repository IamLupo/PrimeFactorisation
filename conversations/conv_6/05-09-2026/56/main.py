#!/usr/bin/env python3

import math
import random
from fractions import Fraction


# ========================================================================
# START EXPERIMENT 139
# Multi-column rational intersection for quotient ratios
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 139")
print("Multi-column rational intersection for quotient ratios")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43, 59]
R2 = [19, 37, 61]

print()
print("R1 =", R1)
print("R2 =", R2)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = x * x % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        p, q = sorted((p, q))

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ------------------------------------------------------------------------
# CELL
# ------------------------------------------------------------------------

def cell(n, p, q, r, s):

    k, a = divmod(p, r)
    ell, b = divmod(q, s)

    c1 = (k * b) // s
    c2 = (ell * a) // r

    beta = (k * b) % s
    alpha = (ell * a) % r

    c3 = (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)

    E = c1 + c2 + c3

    Q = n // (r * s)
    K = Q - E

    assert K == k * ell

    return {
        "Q": Q,
        "K": K,
        "E": E,
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c3": c3,
    }


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

def build_grid(n, p, q):

    g = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            g[(i, j)] = cell(
                n,
                p,
                q,
                r,
                s,
            )

    return g


# ------------------------------------------------------------------------
# K BOUNDS
#
# E = c1+c2+c3
#
# c1 < k
# c2 < ell
# c3 <= 2
#
# Hence:
#
#     E <= k + ell + 2.
#
# Use observable conservative bounds.
# ------------------------------------------------------------------------

def E_bound(n, r, s):

    root = math.isqrt(n)

    kmax = root // r + 1
    lmax = root // s + 1

    return kmax + lmax + 2


# ------------------------------------------------------------------------
# RATIO INTERVAL
#
# K_i/K_0
#
# with:
#
#     Q_i - B_i <= K_i <= Q_i
#
#     Q_0 - B_0 <= K_0 <= Q_0
#
# all positive.
#
# Therefore:
#
#     (Q_i-B_i)/Q_0
#       <= K_i/K_0 <=
#
#     Q_i/(Q_0-B_0).
# ------------------------------------------------------------------------

def ratio_interval(
    Qi,
    Bi,
    Q0,
    B0,
):

    low = Fraction(
        Qi - Bi,
        Q0,
    )

    high = Fraction(
        Qi,
        Q0 - B0,
    )

    return low, high


# ------------------------------------------------------------------------
# INTERSECT RATIONAL INTERVALS
# ------------------------------------------------------------------------

def intersect_intervals(intervals):

    if not intervals:
        return None

    low = max(
        x[0]
        for x in intervals
    )

    high = min(
        x[1]
        for x in intervals
    )

    if low > high:
        return None

    return low, high


# ------------------------------------------------------------------------
# CONTINUED FRACTION
# ------------------------------------------------------------------------

def continued_fraction(x):

    num = x.numerator
    den = x.denominator

    result = []

    while den:

        a = num // den

        result.append(a)

        num, den = den, num - a * den

    return result


# ------------------------------------------------------------------------
# CONVERGENTS
# ------------------------------------------------------------------------

def convergents(cf):

    pm2, pm1 = 0, 1
    qm2, qm1 = 1, 0

    out = []

    for a in cf:

        p = a * pm1 + pm2
        q = a * qm1 + qm2

        out.append(
            Fraction(p, q)
        )

        pm2, pm1 = pm1, p
        qm2, qm1 = qm1, q

    return out


# ------------------------------------------------------------------------
# SEMICONVERGENTS
# ------------------------------------------------------------------------

def semiconvergents(cf):

    out = []

    pm2, pm1 = 0, 1
    qm2, qm1 = 1, 0

    for a in cf:

        for t in range(1, a + 1):

            p = t * pm1 + pm2
            q = t * qm1 + qm2

            out.append(
                Fraction(p, q)
            )

        pm2, pm1 = pm1, a * pm1 + pm2
        qm2, qm1 = qm1, a * qm1 + qm2

    return out


# ------------------------------------------------------------------------
# FRACTIONS IN INTERVAL
#
# We don't enumerate denominators.
#
# We generate continued-fraction candidates from the midpoint.
# ------------------------------------------------------------------------

def rational_candidates_in_interval(
    low,
    high,
    denominator_bound,
):

    if low > high:
        return []

    midpoint = (
        low + high
    ) / 2

    cf = continued_fraction(
        midpoint
    )

    candidates = set()

    for x in convergents(cf):
        if x.denominator <= denominator_bound:
            if low <= x <= high:
                candidates.add(x)

    for x in semiconvergents(cf):
        if x.denominator <= denominator_bound:
            if low <= x <= high:
                candidates.add(x)

    return sorted(candidates)


# ------------------------------------------------------------------------
# RADIX COMPATIBILITY
#
# Suppose:
#
#     k_i/k_0 = u/v
#
# in lowest terms.
#
# Then:
#
#     k_0 = v*t
#     k_i = u*t
#
# and:
#
#     r0*k0 - ri*ki = a_i-a_0
#
# hence:
#
#     t*(r0*v-ri*u) = delta
#
# with small delta.
# ------------------------------------------------------------------------

def compatible_scales(
    u,
    v,
    r0,
    ri,
):

    if u <= 0 or v <= 0:
        return []

    if math.gcd(u, v) != 1:
        return []

    d = (
        r0 * v
        - ri * u
    )

    max_delta = max(
        r0 - 1,
        ri - 1,
    )

    results = []

    if d == 0:

        results.append(
            {
                "u": u,
                "v": v,
                "t": None,
                "k0": None,
                "ki": None,
                "delta": 0,
            }
        )

        return results

    for delta in range(
        -max_delta,
        max_delta + 1,
    ):

        if delta % d != 0:
            continue

        t = delta // d

        if t <= 0:
            continue

        results.append(
            {
                "u": u,
                "v": v,
                "t": t,
                "k0": v * t,
                "ki": u * t,
                "delta": delta,
            }
        )

    return results


# ------------------------------------------------------------------------
# MULTI-COLUMN ANALYSIS FOR ONE i
# ------------------------------------------------------------------------

def analyze_ratio(
    n,
    g,
    i,
):

    intervals = []

    for j in range(3):

        Qi = g[(i, j)]["Q"]
        Q0 = g[(0, j)]["Q"]

        Bi = E_bound(
            n,
            R1[i],
            R2[j],
        )

        B0 = E_bound(
            n,
            R1[0],
            R2[j],
        )

        lo, hi = ratio_interval(
            Qi,
            Bi,
            Q0,
            B0,
        )

        intervals.append(
            (
                lo,
                hi,
                j,
            )
        )

    intersection = intersect_intervals(
        [
            (x[0], x[1])
            for x in intervals
        ]
    )

    return intervals, intersection


# ------------------------------------------------------------------------
# DIRECT CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    for p in range(
        2,
        math.isqrt(n) + 1,
    ):

        if n % p == 0:
            return p, n // p

    return None, None


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(139)

TEST_BITS = [
    30,
    36,
    42,
    48,
    54,
]


for bits in TEST_BITS:

    print()
    print("=" * 72)
    print(
        "GENERATING",
        bits,
        "-BIT SEMIPRIME"
    )
    print("=" * 72)

    n, p, q = random_semiprime(bits)

    g = build_grid(
        n,
        p,
        q,
    )

    true_k = [
        g[(i, 0)]["k"]
        for i in range(3)
    ]

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print(
        "n      =",
        n
    )

    print(
        "true p =",
        p
    )

    print(
        "true q =",
        q
    )

    print()
    print(
        "TRUE k =",
        true_k
    )

    # ------------------------------------------------------------
    # Direct scan count.
    # ------------------------------------------------------------

    direct_count = math.isqrt(n) - 1

    print()
    print(
        "DIRECT SEARCH SCALE ~",
        direct_count
    )

    # ------------------------------------------------------------
    # Analyze k1/k0 and k2/k0.
    # ------------------------------------------------------------

    recovered_sets = []

    for i in (1, 2):

        intervals, intersection = analyze_ratio(
            n,
            g,
            i,
        )

        true_ratio = Fraction(
            true_k[i],
            true_k[0],
        )

        print()
        print(
            f"k[{i}]/k[0]"
        )

        print(
            "true ratio =",
            true_ratio,
            "=",
            float(true_ratio),
        )

        for lo, hi, j in intervals:

            print(
                "column",
                j,
                ":"
            )

            print(
                "  low  =",
                lo,
                "=",
                float(lo),
            )

            print(
                "  high =",
                hi,
                "=",
                float(hi),
            )

            print(
                "  width =",
                float(hi - lo),
            )

            print(
                "  true inside =",
                lo <= true_ratio <= hi,
            )

        if intersection is None:

            print()
            print(
                "INTERSECTION EMPTY"
            )

            recovered_sets.append([])

            continue

        low, high = intersection

        print()
        print(
            "INTERSECTION"
        )

        print(
            "low  =",
            low,
            "=",
            float(low),
        )

        print(
            "high =",
            high,
            "=",
            float(high),
        )

        print(
            "width =",
            float(high - low),
        )

        print(
            "true inside =",
            low <= true_ratio <= high,
        )

        # --------------------------------------------------------
        # CF reconstruction.
        # --------------------------------------------------------

        denominator_bound = (
            math.isqrt(n)
            // R1[0]
            + 2
        )

        candidates = rational_candidates_in_interval(
            low,
            high,
            denominator_bound,
        )

        print()
        print(
            "CF candidates =",
            len(candidates)
        )

        found_true = (
            true_ratio in candidates
        )

        print(
            "TRUE RATIO RECOVERED =",
            found_true
        )

        compatible = []

        for frac in candidates:

            u = frac.numerator
            v = frac.denominator

            states = compatible_scales(
                u,
                v,
                R1[0],
                R1[i],
            )

            compatible.extend(
                (
                    frac,
                    state,
                )
                for state in states
            )

        print(
            "radix-compatible states =",
            len(compatible)
        )

        if compatible:

            print()
            print(
                "FIRST COMPATIBLE STATES"
            )

            for frac, state in compatible[:15]:

                exact_ratio = (
                    frac
                    == true_ratio
                )

                exact_k = (
                    state["k0"]
                    == true_k[0]
                    and
                    state["ki"]
                    == true_k[i]
                )

                print(
                    frac,
                    "t =",
                    state["t"],
                    "k0 =",
                    state["k0"],
                    "ki =",
                    state["ki"],
                    "delta =",
                    state["delta"],
                    "TRUE_RATIO =",
                    exact_ratio,
                    "TRUE_SCALE =",
                    exact_k,
                )

        recovered_sets.append(
            compatible
        )

    # ------------------------------------------------------------
    # Combine i=1 and i=2.
    # ------------------------------------------------------------

    print()
    print(
        "CROSS-RATIO INTERSECTION"
    )

    combined = []

    if recovered_sets[0] and recovered_sets[1]:

        for frac1, s1 in recovered_sets[0]:

            if s1["k0"] is None:
                continue

            for frac2, s2 in recovered_sets[1]:

                if s2["k0"] is None:
                    continue

                if (
                    s1["k0"]
                    != s2["k0"]
                ):
                    continue

                combined.append(
                    (
                        frac1,
                        frac2,
                        s1,
                        s2,
                    )
                )

    print(
        "combined states =",
        len(combined)
    )

    true_vector_found = False

    for frac1, frac2, s1, s2 in combined:

        if (
            s1["k0"]
            == true_k[0]
            and s1["ki"]
            == true_k[1]
            and s2["ki"]
            == true_k[2]
        ):

            true_vector_found = True

            print()
            print(
                ">>> TRUE k VECTOR RECOVERED <<<"
            )

            print(
                "k =",
                (
                    s1["k0"],
                    s1["ki"],
                    s2["ki"],
                )
            )

            print(
                "ratio1 =",
                frac1
            )

            print(
                "ratio2 =",
                frac2
            )

            print(
                "delta1 =",
                s1["delta"]
            )

            print(
                "delta2 =",
                s2["delta"]
            )

            break

    print(
        "TRUE k VECTOR FOUND =",
        true_vector_found
    )

    # ------------------------------------------------------------
    # Show first combined states.
    # ------------------------------------------------------------

    if combined:

        print()
        print(
            "FIRST COMBINED STATES"
        )

        for frac1, frac2, s1, s2 in combined[:10]:

            print(
                "k =",
                (
                    s1["k0"],
                    s1["ki"],
                    s2["ki"],
                ),
                "ratios =",
                (
                    frac1,
                    frac2,
                ),
            )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 139
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 139")
print("=" * 72)
