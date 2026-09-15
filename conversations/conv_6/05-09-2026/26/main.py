#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 109
# Divisibility elimination
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20

RANDOM_SEED = 109


# ============================================================
# PRIME TEST
# ============================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    )

    for p in small:

        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in (
        2, 3, 5, 7,
        11, 13, 17
    ):

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


def random_prime(bits):

    while True:

        p = random.getrandbits(bits)

        p |= 1 << (bits - 1)
        p |= 1

        if is_probable_prime(p):
            return p


def make_semiprime(bits):

    fb = bits // 2

    while True:

        p = random_prime(fb)
        q = random_prime(fb)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() >= bits - 1:
            return n, p, q


# ============================================================
# CARRY
# ============================================================

def carries(r1, r2, k, l, a, b):

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    return c1, c2, c3, c1 + c2 + c3


# ============================================================
# FORCE q MOD R00
# ============================================================

def q_residue(n, p_res):

    if math.gcd(p_res, R00) != 1:
        return None

    return (
        n
        * pow(p_res, -1, R00)
    ) % R00


def decode_q(q0):

    # q = 19*l + b
    #
    # q0 = 19*t + b
    #
    # where t = l mod 17.

    b = q0 % r20
    t = q0 // r20

    return t, b


# ============================================================
# BUILD THE DIVISIBILITY CONSTANT
# ============================================================

def build_constant(n, a, s):

    """
    Fix

        a = p mod 17
        s = k mod 19

    so

        k = 19u + s
        p = 323u + (17s+a)

    The modular equation determines q mod 323 and therefore
    fixes b and t.

    The first-cell equation becomes

        A - B*u = p*m

    and therefore

        p | 323*A + B*C

    where

        C = 17*s + a.

    Returns all useful quantities.
    """

    Q = n // R00

    C = 17 * s + a

    p_res = C % R00

    q0 = q_residue(
        n,
        p_res
    )

    if q0 is None:
        return None

    t, b = decode_q(q0)

    # --------------------------------------------------------
    # beta = k*b mod 19.
    #
    # Since k = 19u+s, beta is fixed by s.
    # --------------------------------------------------------

    beta = (
        s * b
    ) % r20

    c1_const = (
        s * b
    ) // r20

    # --------------------------------------------------------
    # alpha = t*a mod 17.
    # --------------------------------------------------------

    alpha = (
        t * a
    ) % r10

    c2_const = (
        t * a
    ) // r10

    c3 = (
        r10 * beta
        + r20 * alpha
        + a * b
    ) // R00

    # --------------------------------------------------------
    # First-cell equation:
    #
    # Q - (19u+s)t
    # - (u*b + floor(s*b/19))
    # - c2_const
    # - c3
    #
    # =
    #
    # (323u + C)m
    #
    # Hence:
    #
    # A - B*u = (323u+C)m
    # --------------------------------------------------------

    A = (
        Q
        - s * t
        - c1_const
        - c2_const
        - c3
    )

    B = (
        19 * t + b
    )

    H = (
        323 * A
        + B * C
    )

    return {
        "a": a,
        "s": s,
        "C": C,
        "q0": q0,
        "t": t,
        "b": b,
        "A": A,
        "B": B,
        "H": H,
    }


# ============================================================
# VERIFY A GCD CANDIDATE
# ============================================================

def verify_candidate(n, p, info):

    C = info["C"]

    # C + 323u = p
    if p < C:
        return None

    if (p - C) % 323 != 0:
        return None

    u = (p - C) // 323

    if u < 0:
        return None

    k = 19 * u + info["s"]

    # Reconstruct q residue.

    q0 = info["q0"]

    b = info["b"]
    t = info["t"]

    # Solve actual l.

    numerator = (
        info["A"]
        - info["B"] * u
    )

    if numerator < 0:
        return None

    if numerator % p != 0:
        return None

    m = numerator // p

    # Here m is the quotient from the linearized equation.
    # Recover l directly:

    l = t + r10 * m

    if l <= 0:
        return None

    a = info["a"]

    reconstructed_p = r10 * k + a

    if reconstructed_p != p:
        return None

    q = r20 * l + b

    # Full first-cell verification.

    _, _, _, Ecarry = carries(
        r10,
        r20,
        k,
        l,
        a,
        b
    )

    Equot = (
        n // R00
        - k * l
    )

    if Ecarry != Equot:
        return None

    return {
        "p": p,
        "q": q,
        "k": k,
        "l": l,
        "a": a,
        "b": b,
        "t": t,
        "m": m,
        "E": Equot,
    }


# ============================================================
# MAIN DIVISIBILITY SEARCH
# ============================================================

def search(n):

    gcd_tests = 0
    nontrivial_gcds = 0
    reconstructed = []
    verified = []

    start = time.perf_counter()

    # --------------------------------------------------------
    # 17 possible a values
    # 19 possible k residues
    # --------------------------------------------------------

    for a in range(r10):

        for s in range(r20):

            info = build_constant(
                n,
                a,
                s
            )

            if info is None:
                continue

            H = abs(info["H"])

            if H == 0:
                continue

            g = math.gcd(
                n,
                H
            )

            gcd_tests += 1

            # Ignore trivial gcds.

            if g == 1 or g == n:
                continue

            nontrivial_gcds += 1

            # The gcd may contain p, q, or both.
            #
            # Try the gcd and complementary factor.

            candidates = {
                g,
                n // g,
            }

            for candidate in candidates:

                result = verify_candidate(
                    n,
                    candidate,
                    info
                )

                if result is None:
                    continue

                reconstructed.append(result)

                if result["p"] * result["q"] == n:

                    verified.append(result)

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "gcd_tests": gcd_tests,
        "nontrivial_gcds": nontrivial_gcds,
        "reconstructed": reconstructed,
        "verified": verified,
        "runtime": elapsed,
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    result
):

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print()
    print(
        "n      =",
        n
    )

    print(
        "true p =",
        true_p
    )

    print(
        "true q =",
        true_q
    )

    print()
    print("GRID")

    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)

    print()
    print("SEARCH RESULTS")

    print(
        "gcd tests            =",
        result["gcd_tests"]
    )

    print(
        "nontrivial gcds      =",
        result["nontrivial_gcds"]
    )

    print(
        "reconstructed        =",
        len(result["reconstructed"])
    )

    print(
        "verified             =",
        len(result["verified"])
    )

    print(
        "runtime              = %.9f s"
        % result["runtime"]
    )

    print()
    print("RECONSTRUCTED")

    for r in result["reconstructed"][:20]:

        print(
            "  "
            f"p={r['p']} "
            f"q={r['q']} "
            f"k={r['k']} "
            f"l={r['l']} "
            f"a={r['a']} "
            f"b={r['b']} "
            f"t={r['t']} "
            f"E={r['E']} "
            f"exact={r['p'] * r['q'] == n}"
        )

    if len(result["reconstructed"]) > 20:

        print(
            "  ..."
            f"{len(result['reconstructed']) - 20}"
            " more"
        )

    print()
    print("VERIFIED FACTORIZATIONS")

    for r in result["verified"]:

        print(
            "  "
            f"{r['p']} * {r['q']} = "
            f"{r['p'] * r['q']}"
        )

    print()
    print("TRUE FACTORIZATION FOUND =")

    print(
        any(
            (
                r["p"] == true_p
                and
                r["q"] == true_q
            )
            or
            (
                r["p"] == true_q
                and
                r["q"] == true_p
            )
            for r in result["verified"]
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 109")
    print("Divisibility elimination")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)

    for bits in (30, 36, 42, 48):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        result = search(n)

        report(
            n,
            p,
            q,
            result
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 109")
    print("=" * 72)


if __name__ == "__main__":
    main()