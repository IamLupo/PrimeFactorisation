#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 118
# Quotient-variable bilinear search
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

R00 = R1[0] * R2[0]
R11 = R1[1] * R2[1]

M = math.lcm(R00, R11)

RANDOM_SEED = 118


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
# CRT HELPERS
# ============================================================

def solve_linear_congruence(a, b, m):
    """
    Solve

        a*x == b (mod m)

    """

    g = math.gcd(a, m)

    if b % g != 0:
        return []

    aa = a // g
    bb = b // g
    mm = m // g

    if mm == 1:
        return [0]

    inv = pow(
        aa % mm,
        -1,
        mm
    )

    x0 = (
        bb * inv
    ) % mm

    return [
        x0 + t * mm
        for t in range(g)
    ]


def q_residues(n, p):

    """
    Solve

        p*q == n (mod M)
    """

    return solve_linear_congruence(
        p,
        n,
        M
    )


# ============================================================
# RESIDUE PAIR GENERATION
# ============================================================

def residue_pairs(n):

    """
    Find all invertible residue pairs

        a*b == n (mod M)

    with

        0 <= a,b < M.

    Instead of M^2 enumeration, enumerate a and solve for b.
    """

    pairs = []

    n_mod = n % M

    for a in range(M):

        if math.gcd(a, M) != 1:
            continue

        b = (
            n_mod
            * pow(a, -1, M)
        ) % M

        pairs.append(
            (a, b)
        )

    return pairs


# ============================================================
# BILINEAR SEARCH FOR FIXED RESIDUE PAIR
# ============================================================

def solve_residue_pair(n, a, b):

    """
    p = a + M*u
    q = b + M*v

    with

        p*q = n.

    Expansion:

        M^2*u*v
        + M*a*v
        + M*b*u
        + ab
        = n

    Divide by M:

        M*u*v + a*v + b*u
        = (n-ab)/M

    Define

        C = (n-ab)/M.

    Then

        v(Mu+a) = C-bu

    so for each u,

        v = (C-bu)/(Mu+a).

    The purpose of the experiment is to determine whether
    the u-space is significantly smaller than the original
    p-space.
    """

    if (n - a * b) % M != 0:
        return []

    C = (
        n - a * b
    ) // M

    # p <= sqrt(n)
    p_limit = math.isqrt(n)

    if a > p_limit:
        return []

    u_max = (
        p_limit - a
    ) // M

    candidates = []

    for u in range(
        u_max + 1
    ):

        p = a + M * u

        if p < 2:
            continue

        denominator = (
            M * u + a
        )

        numerator = (
            C - b * u
        )

        if numerator < 0:
            continue

        if numerator % denominator != 0:
            continue

        v = (
            numerator
            // denominator
        )

        if v < 0:
            continue

        q = b + M * v

        if p > q:
            continue

        candidates.append(
            (p, q, u, v)
        )

    return candidates


# ============================================================
# RESIDUE-CLASS SEARCH
# ============================================================

def bilinear_search(n):

    pairs = residue_pairs(n)

    pair_tests = 0
    u_tests = 0
    candidates = []

    start = time.perf_counter()

    p_limit = math.isqrt(n)

    for a, b in pairs:

        pair_tests += 1

        if a > p_limit:
            continue

        u_max = (
            p_limit - a
        ) // M

        u_tests += (
            u_max + 1
        )

        found = solve_residue_pair(
            n,
            a,
            b
        )

        candidates.extend(
            found
        )

    elapsed = (
        time.perf_counter()
        - start
    )

    exact = [
        x
        for x in candidates
        if x[0] * x[1] == n
    ]

    return {
        "pairs": len(pairs),
        "pair_tests": pair_tests,
        "u_tests": u_tests,
        "candidates": candidates,
        "exact": exact,
        "time": elapsed,
    }


# ============================================================
# DIRECT P BASELINE
# ============================================================

def direct_scan(n):

    limit = math.isqrt(n)

    tested = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        tested += 1

        if n % p == 0:

            q = n // p

            return {
                "tested": tested,
                "p": p,
                "q": q,
                "time": (
                    time.perf_counter()
                    - start
                ),
            }

    return {
        "tested": tested,
        "p": None,
        "q": None,
        "time": (
            time.perf_counter()
            - start
        ),
    }


# ============================================================
# MODULAR BASELINE
# ============================================================

def modular_scan(n):

    limit = math.isqrt(n)

    tested = 0
    q_candidates = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        tested += 1

        if math.gcd(p, M) != 1:
            continue

        q_res = (
            n
            * pow(
                p,
                -1,
                M
            )
        ) % M

        q_low = 1

        # Since p*q=n exactly at the solution:
        #
        # q = n/p
        #
        # The experiment deliberately does NOT use that
        # as a rejection shortcut.
        #
        # Instead use the second-cell quotient bucket.

        Q11 = n // R11

        q_low = (
            R11 * Q11 + p - 1
        ) // p

        q_high = (
            R11 * (Q11 + 1) - 1
        ) // p

        q = (
            q_low
            + (q_res - q_low) % M
        )

        while q <= q_high:

            q_candidates += 1

            if p <= q and p * q == n:

                return {
                    "tested": tested,
                    "q_candidates": q_candidates,
                    "p": p,
                    "q": q,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                }

            q += M

    return {
        "tested": tested,
        "q_candidates": q_candidates,
        "p": None,
        "q": None,
        "time": (
            time.perf_counter()
            - start
        ),
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    direct,
    modular,
    bilinear
):

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("MODULUS")

    print("R00 =", R00)
    print("R11 =", R11)
    print("M   =", M)

    print()
    print("TRUE FACTOR RESIDUES")

    print(
        "p mod M =",
        true_p % M
    )

    print(
        "q mod M =",
        true_q % M
    )

    print()
    print("DIRECT P-SCAN")

    print(
        "p tested =",
        direct["tested"]
    )

    print(
        "time     = %.6f s"
        % direct["time"]
    )

    print()
    print("TWO-CELL MODULAR SCAN")

    print(
        "p tested       =",
        modular["tested"]
    )

    print(
        "q candidates   =",
        modular["q_candidates"]
    )

    print(
        "time           = %.6f s"
        % modular["time"]
    )

    print()
    print("BILINEAR RESIDUE SEARCH")

    print(
        "residue pairs  =",
        bilinear["pairs"]
    )

    print(
        "pair tests     =",
        bilinear["pair_tests"]
    )

    print(
        "u tests        =",
        bilinear["u_tests"]
    )

    print(
        "candidates     =",
        len(bilinear["candidates"])
    )

    print(
        "exact          =",
        len(bilinear["exact"])
    )

    print(
        "time           = %.6f s"
        % bilinear["time"]
    )

    print()
    print("BILINEAR RESULTS")

    for p, q, u, v in bilinear["candidates"][:20]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"u={u} "
            f"v={v} "
            f"exact={p*q == n}"
        )

    if len(bilinear["candidates"]) > 20:

        print(
            "  ..."
            f"{len(bilinear['candidates']) - 20}"
            " more"
        )

    print()
    print("TRUE FOUND =")

    print(
        any(
            p == true_p
            and q == true_q
            for p, q, _, _
            in bilinear["exact"]
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 118")
    print("Quotient-variable bilinear search")
    print("=" * 72)

    print()

    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)
    print("R11 =", R11)
    print("M =", M)

    for bits in (
        30,
        36,
        42,
        48,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(
            bits
        )

        direct = direct_scan(
            n
        )

        modular = modular_scan(
            n
        )

        bilinear = bilinear_search(
            n
        )

        report(
            n,
            p,
            q,
            direct,
            modular,
            bilinear
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 118")
    print("=" * 72)


if __name__ == "__main__":
    main()
