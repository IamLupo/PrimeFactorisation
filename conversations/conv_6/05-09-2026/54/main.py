#!/usr/bin/env python3

import math
import random
from fractions import Fraction


# ========================================================================
# START EXPERIMENT 137
# Continued-fraction recovery of quotient ratios
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 137")
print("Continued-fraction recovery of quotient ratios")
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
# CARRY CELL
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

    return {
        "Q": Q,
        "K": K,
        "E": E,
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
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
# CONTINUED FRACTIONS
# ------------------------------------------------------------------------

def continued_fraction(num, den):

    out = []

    while den != 0:

        a = num // den

        out.append(a)

        num, den = den, num - a * den

    return out


def convergents(cf):

    p_minus2 = 0
    p_minus1 = 1

    q_minus2 = 1
    q_minus1 = 0

    result = []

    for a in cf:

        p = (
            a * p_minus1
            + p_minus2
        )

        q = (
            a * q_minus1
            + q_minus2
        )

        result.append(
            (p, q)
        )

        p_minus2, p_minus1 = (
            p_minus1,
            p,
        )

        q_minus2, q_minus1 = (
            q_minus1,
            q,
        )

    return result


# ------------------------------------------------------------------------
# SEMICONVERGENTS
#
# Include intermediate fractions around the convergents.
# ------------------------------------------------------------------------

def semiconvergents(cf):

    result = []

    p0, p1 = 0, 1
    q0, q1 = 1, 0

    for a in cf:

        for t in range(1, a + 1):

            p = t * p1 + p0
            q = t * q1 + q0

            result.append(
                (p, q)
            )

        p0, p1 = p1, a * p1 + p0
        q0, q1 = q1, a * q1 + q0

    return result


# ------------------------------------------------------------------------
# TRUE RATIO TEST
# ------------------------------------------------------------------------

def ratio_error(p, q, x):

    return abs(
        Fraction(p, q)
        - x
    )


# ------------------------------------------------------------------------
# RECOVER RATIONAL APPROXIMATION
#
# Given:
#
#     x ~= numerator / denominator
#
# seek rational candidates with denominator <= BOUND.
#
# We use:
#
#     convergents
#     semiconvergents
#
# but do not scan denominators.
# ------------------------------------------------------------------------

def rational_recovery(
    numerator,
    denominator,
    denominator_bound,
):

    x = Fraction(
        numerator,
        denominator,
    )

    cf = continued_fraction(
        numerator,
        denominator,
    )

    conv = convergents(cf)
    semi = semiconvergents(cf)

    candidates = set()

    for p, q in conv:

        if q <= denominator_bound:
            candidates.add(
                (p, q)
            )

    for p, q in semi:

        if q <= denominator_bound:
            candidates.add(
                (p, q)
            )

    # ------------------------------------------------------------
    # Sort by approximation quality.
    # ------------------------------------------------------------

    candidates = sorted(
        candidates,
        key=lambda z: (
            abs(
                Fraction(z[0], z[1])
                - x
            ),
            z[1],
        ),
    )

    return x, candidates


# ------------------------------------------------------------------------
# TEST RADIX-COMPATIBLE RATIO
#
# A candidate k_i/k0 must satisfy:
#
#     r0*k0 - ri*ki = ai-a0
#
# with:
#
#     -(r0-1) <= ai-a0 <= ri-1.
#
# We don't know k0 yet.
#
# Given candidate ratio ki/k0 = u/v, solve:
#
#     r0*v*t - ri*u*t = delta
#
# where k0=v*t and ki=u*t.
#
# Thus:
#
#     t*(r0*v-ri*u) = delta.
#
# Since delta is tiny, this is a very strong arithmetic test.
# ------------------------------------------------------------------------

def ratio_radix_compatibility(
    u,
    v,
    r0,
    ri,
):

    if u <= 0 or v <= 0:
        return []

    g = (
        r0 * v
        - ri * u
    )

    results = []

    # delta = ai-a0
    #
    # |delta| <= max(radices)-1
    max_delta = max(
        r0 - 1,
        ri - 1,
    )

    if g == 0:

        # Then delta must be zero.
        results.append(
            {
                "u": u,
                "v": v,
                "t": None,
                "delta": 0,
                "k0": None,
                "ki": None,
                "exact_scale": False,
            }
        )

        return results

    for delta in range(
        -max_delta,
        max_delta + 1,
    ):

        if delta % g != 0:
            continue

        t = delta // g

        if t <= 0:
            continue

        k0 = v * t
        ki = u * t

        results.append(
            {
                "u": u,
                "v": v,
                "t": t,
                "delta": delta,
                "k0": k0,
                "ki": ki,
                "exact_scale": True,
            }
        )

    return results


# ------------------------------------------------------------------------
# VERIFY FULL QUOTIENT VECTOR
# ------------------------------------------------------------------------

def verify_k_candidate(
    n,
    p,
    q,
    k0,
    expected_i,
):

    true_g = build_grid(
        n,
        p,
        q,
    )

    true_k = [
        true_g[(i, 0)]["k"]
        for i in range(3)
    ]

    return (
        k0 == true_k[0]
        and expected_i == true_k[1:]
    )


# ------------------------------------------------------------------------
# DIRECT CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    tested = 0

    for p in range(
        2,
        math.isqrt(n) + 1,
    ):

        tested += 1

        if n % p == 0:
            return (
                p,
                n // p,
                tested,
            )

    return None, None, tested


# ------------------------------------------------------------------------
# MAIN EXPERIMENT
# ------------------------------------------------------------------------

random.seed(137)

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
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, p, q = random_semiprime(bits)

    g = build_grid(
        n,
        p,
        q,
    )

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    true_k = [
        g[(i, 0)]["k"]
        for i in range(3)
    ]

    print()
    print("TRUE k VECTOR")
    print(
        "k =",
        true_k
    )

    # ------------------------------------------------------------
    # Direct control.
    # ------------------------------------------------------------

    _, _, direct_tests = direct_factor(n)

    print()
    print("DIRECT SQRT(n) CONTROL")
    print(
        "p tested =",
        direct_tests
    )

    # ------------------------------------------------------------
    # k1/k0 and k2/k0 from observable Q ratios.
    # ------------------------------------------------------------

    Q00 = g[(0, 0)]["Q"]
    Q10 = g[(1, 0)]["Q"]
    Q20 = g[(2, 0)]["Q"]

    print()
    print("OBSERVABLE k RATIOS")

    ratios = [
        (
            1,
            Q10,
            Q00,
        ),
        (
            2,
            Q20,
            Q00,
        ),
    ]

    recovered_ratio_sets = []

    k0_bound = (
        math.isqrt(n)
        // R1[0]
    ) + 1

    for i, num, den in ratios:

        observed = Fraction(
            num,
            den,
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
            "observed =",
            observed,
        )

        print(
            "true     =",
            true_ratio,
        )

        print(
            "absolute error =",
            float(
                abs(
                    observed
                    - true_ratio
                )
            ),
        )

        # --------------------------------------------------------
        # Rational reconstruction.
        # --------------------------------------------------------

        x, candidates = rational_recovery(
            num,
            den,
            k0_bound,
        )

        print(
            "CF candidates =",
            len(candidates)
        )

        found_true_ratio = False

        compatible = []

        for u, v in candidates:

            if math.gcd(u, v) != 1:
                continue

            # Candidate fraction may be >1 etc.
            # Keep physically sensible positive ratios.
            if u <= 0 or v <= 0:
                continue

            if Fraction(u, v) > 1:
                continue

            compat = ratio_radix_compatibility(
                u,
                v,
                R1[0],
                R1[i],
            )

            if compat:

                compatible.extend(
                    compat
                )

            if (
                Fraction(u, v)
                == true_ratio
            ):

                found_true_ratio = True

        print(
            "TRUE RATIO FOUND BY CF =",
            found_true_ratio
        )

        print(
            "radix-compatible candidates =",
            len(compatible)
        )

        # Show first few candidates.
        if candidates:

            print()
            print("BEST CF CANDIDATES")

            for u, v in candidates[:10]:

                print(
                    f"{u}/{v}",
                    "=",
                    float(
                        Fraction(u, v)
                    ),
                    "error =",
                    float(
                        abs(
                            Fraction(u, v)
                            - observed
                        )
                    ),
                )

        recovered_ratio_sets.append(
            compatible
        )

    # ------------------------------------------------------------
    # Cross-ratio consistency.
    #
    # If we recover:
    #
    #     k1/k0 = u1/v1
    #     k2/k0 = u2/v2
    #
    # then the same k0 scale must satisfy both radix equations.
    # ------------------------------------------------------------

    print()
    print("CROSS-RATIO SCALE CONSISTENCY")

    if all(recovered_ratio_sets):

        combined = []

        for c1 in recovered_ratio_sets[0]:

            for c2 in recovered_ratio_sets[1]:

                if (
                    c1["k0"] is None
                    or c2["k0"] is None
                ):
                    continue

                if c1["k0"] == c2["k0"]:

                    combined.append(
                        (
                            c1,
                            c2,
                        )
                    )

        print(
            "combined consistent states =",
            len(combined)
        )

        true_scale_found = False

        for c1, c2 in combined:

            if (
                c1["k0"]
                == true_k[0]
                and c1["ki"]
                == true_k[1]
                and c2["ki"]
                == true_k[2]
            ):

                true_scale_found = True

                print()
                print(
                    ">>> TRUE k VECTOR RECOVERED <<<"
                )

                print(
                    "k0 =",
                    c1["k0"]
                )

                print(
                    "k1 =",
                    c1["ki"]
                )

                print(
                    "k2 =",
                    c2["ki"]
                )

                print(
                    "delta1 =",
                    c1["delta"]
                )

                print(
                    "delta2 =",
                    c2["delta"]
                )

                break

        print(
            "TRUE k VECTOR FOUND =",
            true_scale_found
        )

        if combined:

            print()
            print("FIRST COMBINED STATES")

            for c1, c2 in combined[:10]:

                print(
                    "k =",
                    (
                        c1["k0"],
                        c1["ki"],
                        c2["ki"],
                    ),
                    "deltas =",
                    (
                        c1["delta"],
                        c2["delta"],
                    )
                )

    else:

        print(
            "No combined states."
        )

    # ------------------------------------------------------------
    # Direct exact ratio comparison.
    # ------------------------------------------------------------

    print()
    print("TRUE RATIO DATA")

    for i in (1, 2):

        true_ratio = Fraction(
            true_k[i],
            true_k[0],
        )

        observed = Fraction(
            g[(i, 0)]["Q"],
            g[(0, 0)]["Q"],
        )

        print(
            f"k[{i}]/k[0] =",
            true_ratio,
            "observed =",
            observed,
        )

    print()
    print("=" * 72)


# ========================================================================
# FINISHED EXPERIMENT 137
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 137")
print("=" * 72)
