#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 127
# Bilinear factor transform -> constrained divisor search
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 127")
print("Bilinear factor transform -> constrained divisor search")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
# ------------------------------------------------------------------------

R1 = [17, 43]
R2 = [19, 37]

r = R1[0]
s = R2[0]

R = r * s

print()
print("r =", r)
print("s =", s)
print("R = r*s =", R)


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
    s0 = 0

    while d % 2 == 0:
        d //= 2
        s0 += 1

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

        for _ in range(s0 - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# RANDOM PRIME
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

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() == bits:
            return n, p, q


# ------------------------------------------------------------------------
# DIRECT FACTOR CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    limit = math.isqrt(n)

    tested = 0

    for p in range(2, limit + 1):

        tested += 1

        if n % p == 0:

            return p, n // p, tested

    return None, None, tested


# ------------------------------------------------------------------------
# CARRY STATE
# ------------------------------------------------------------------------

def carry_state(p, q, r1, r2):

    k, a = divmod(p, r1)
    ell, b = divmod(q, r2)

    c1 = (k * b) // r2
    c2 = (ell * a) // r1

    beta = (k * b) % r2
    alpha = (ell * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    E = c1 + c2 + c3

    Q = (p * q) // (r1 * r2)

    return {
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "E": E,
        "Q": Q,
        "K": Q - E,
    }


# ------------------------------------------------------------------------
# NEW TRANSFORMATION
#
# Starting with:
#
#     p = r*x + a
#     q = s*y + b
#
# we have:
#
#     n = (r*x+a)(s*y+b)
#
# Multiply by r*s:
#
#     r*s*n = (s*p)(r*q)
#
# Define:
#
#     U = s*p
#     V = r*q
#
# Then:
#
#     U*V = r*s*n
#
# and because:
#
#     p == a (mod r)
#
# we have:
#
#     U = s*p == s*a (mod r*s)
#
# Similarly:
#
#     V == r*b (mod r*s)
#
# Therefore this is a divisor search for:
#
#     N' = r*s*n
#
# restricted to one pair of residue classes modulo r*s.
#
# The question:
#
#     Does this representation actually remove work?
#
# ------------------------------------------------------------------------


def constrained_divisor_search(n, true_p, true_q):

    N = R * n

    limit = math.isqrt(N)

    residue_classes = 0
    U_tested = 0
    divisibility_hits = 0
    exact_hits = []

    # ------------------------------------------------------------
    # Enumerate possible a,b for the base radix pair.
    # ------------------------------------------------------------

    for a in range(r):

        for b in range(s):

            residue_classes += 1

            U0 = (s * a) % R

            # ----------------------------------------------------
            # We need:
            #
            #     U == s*a mod R
            #
            # Search U <= sqrt(N).
            # ----------------------------------------------------

            if U0 == 0:

                U = R

            else:

                U = U0

                if U < 2:
                    U += (
                        ((2 - U + R - 1) // R)
                        * R
                    )

            if U > limit:
                continue

            while U <= limit:

                U_tested += 1

                if N % U == 0:

                    divisibility_hits += 1

                    V = N // U

                    # ------------------------------------------------
                    # Both factors must have the expected residue.
                    # ------------------------------------------------

                    if V % R == (r * b) % R:

                        p = U // s
                        q = V // r

                        # U must actually be divisible by s.
                        # V must actually be divisible by r.
                        if (
                            U % s == 0
                            and V % r == 0
                        ):

                            if p * q == n:

                                if (
                                    {p, q}
                                    == {true_p, true_q}
                                ):
                                    exact_hits.append(
                                        (
                                            p,
                                            q,
                                            a,
                                            b,
                                            U,
                                            V,
                                        )
                                    )

                U += R

    return {
        "N": N,
        "limit": limit,
        "residue_classes": residue_classes,
        "U_tested": U_tested,
        "divisibility_hits": divisibility_hits,
        "exact_hits": exact_hits,
    }


# ------------------------------------------------------------------------
# SAME SEARCH, BUT IMPLEMENTED DIRECTLY IN p-SPACE
#
# This is included as an algebraic control.
#
# Since U = s*p and U advances by r*s,
#
#     U = s*p
#
# means p itself advances by r.
#
# Therefore the constrained-U search may simply be the ordinary
# p search reorganized by residue class.
# ------------------------------------------------------------------------

def direct_residue_class_scan(
    n,
    true_p,
    true_q
):

    limit = math.isqrt(n)

    tested = 0
    residue_classes = 0
    exact_hits = []

    for a in range(r):

        residue_classes += 1

        p = a

        if p < 2:
            p += (
                ((2 - p + r - 1) // r)
                * r
            )

        while p <= limit:

            tested += 1

            if n % p == 0:

                q = n // p

                if (
                    {p, q}
                    == {true_p, true_q}
                ):
                    exact_hits.append(
                        (
                            p,
                            q,
                            a,
                        )
                    )

            p += r

    return {
        "residue_classes": residue_classes,
        "p_tested": tested,
        "exact_hits": exact_hits,
    }


# ------------------------------------------------------------------------
# FULL GRID VERIFICATION
# ------------------------------------------------------------------------

def verify_full_grid(p, q):

    expected = build_grid(p, q)

    for key, item in expected.items():

        r1, r2 = key

        C = carry_state(
            p,
            q,
            r1,
            r2,
        )

        if (
            C["Q"] != item["Q"]
            or C["E"] != item["E"]
            or C["K"] != item["K"]
        ):
            return False

    return True


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

def build_grid(p, q):

    out = {}

    for r1 in R1:

        for r2 in R2:

            C = carry_state(
                p,
                q,
                r1,
                r2,
            )

            out[(r1, r2)] = C

    return out


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(127)

TEST_BITS = [
    30,
    36,
    42,
    48,
]

for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, true_p, true_q = random_semiprime(bits)

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("TRUE RESIDUES")

    print(
        "a = p mod r =",
        true_p % r
    )

    print(
        "b = q mod s =",
        true_q % s
    )

    print()
    print("TRANSFORMED FACTORS")

    U_true = s * true_p
    V_true = r * true_q

    print("U = s*p =", U_true)
    print("V = r*q =", V_true)

    print()
    print("CHECK")

    print(
        "U*V == r*s*n =",
        U_true * V_true == R * n
    )

    print(
        "U mod R =",
        U_true % R
    )

    print(
        "s*a mod R =",
        (s * (true_p % r)) % R
    )

    print(
        "V mod R =",
        V_true % R
    )

    print(
        "r*b mod R =",
        (r * (true_q % s)) % R
    )

    # ------------------------------------------------------------
    # Direct control
    # ------------------------------------------------------------

    print()
    print("DIRECT SQRT(n) SCAN")

    t0 = time.perf_counter()

    pd, qd, tested_direct = direct_factor(n)

    t1 = time.perf_counter()

    direct_time = t1 - t0

    print(
        "p tested =",
        tested_direct
    )

    print(
        "runtime =",
        f"{direct_time:.6f}",
        "s"
    )

    # ------------------------------------------------------------
    # Direct residue-class control
    # ------------------------------------------------------------

    print()
    print("DIRECT RESIDUE-CLASS SCAN")

    t0 = time.perf_counter()

    residue_result = direct_residue_class_scan(
        n,
        true_p,
        true_q,
    )

    t1 = time.perf_counter()

    residue_time = t1 - t0

    print(
        "residue classes =",
        residue_result["residue_classes"]
    )

    print(
        "p tested =",
        residue_result["p_tested"]
    )

    print(
        "runtime =",
        f"{residue_time:.6f}",
        "s"
    )

    # ------------------------------------------------------------
    # New transformed search
    # ------------------------------------------------------------

    print()
    print("CONSTRAINED U-V DIVISOR SEARCH")

    t0 = time.perf_counter()

    result = constrained_divisor_search(
        n,
        true_p,
        true_q,
    )

    t1 = time.perf_counter()

    transformed_time = t1 - t0

    print()
    print(
        "N' = r*s*n =",
        result["N"]
    )

    print(
        "sqrt(N') =",
        result["limit"]
    )

    print(
        "residue classes =",
        result["residue_classes"]
    )

    print(
        "U tested =",
        result["U_tested"]
    )

    print(
        "divisibility hits =",
        result["divisibility_hits"]
    )

    print(
        "exact hits =",
        len(result["exact_hits"])
    )

    print(
        "runtime =",
        f"{transformed_time:.6f}",
        "s"
    )

    recovered = any(
        {x[0], x[1]}
        == {true_p, true_q}
        for x in result["exact_hits"]
    )

    print(
        "TRUE FACTORIZATION FOUND =",
        recovered
    )

    # ------------------------------------------------------------
    # Candidate ratios
    # ------------------------------------------------------------

    print()
    print("SEARCH RATIOS")

    if residue_result["p_tested"] > 0:

        print(
            "transformed / residue-class =",
            f"{transformed_time / residue_time:.3f}x"
        )

    if tested_direct > 0:

        print(
            "transformed / direct =",
            f"{transformed_time / direct_time:.3f}x"
        )

    # ------------------------------------------------------------
    # Print exact result
    # ------------------------------------------------------------

    if result["exact_hits"]:

        print()
        print("EXACT RECOVERY")

        for item in result["exact_hits"]:

            p, q, a, b, U, V = item

            print("p =", p)
            print("q =", q)
            print("a =", a)
            print("b =", b)
            print("U =", U)
            print("V =", V)

            print(
                "U*V == R*n =",
                U * V == R * n
            )

            print(
                "FULL GRID VERIFIED =",
                verify_full_grid(p, q)
            )


# ========================================================================
# FINISHED EXPERIMENT 127
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 127")
print("=" * 72)
