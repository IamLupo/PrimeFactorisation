#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 115
# Difference constraint from the carry state
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20
R11 = r11 * r21

RANDOM_SEED = 115


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
# DIRECT BASELINE
# ============================================================

def direct_factor(n):

    limit = math.isqrt(n)
    tested = 0

    start = time.perf_counter()

    for p in range(3, limit + 1, 2):

        tested += 1

        if n % p == 0:

            q = n // p

            return {
                "p": p,
                "q": q,
                "tested": tested,
                "time": time.perf_counter() - start,
            }

    return {
        "p": None,
        "q": None,
        "tested": tested,
        "time": time.perf_counter() - start,
    }


# ============================================================
# BUILD FIRST-CELL RESIDUE STATE
# ============================================================

def residue_state(n, p):

    if math.gcd(p, R00) != 1:
        return None

    q0 = (
        n
        * pow(p, -1, R00)
    ) % R00

    k = p // r10
    a = p % r10

    l_mod = q0 // r20
    b = q0 % r20

    return k, a, l_mod, b


# ============================================================
# SEARCH USING d = q-p
# ============================================================

def difference_search(n):

    """
    Exploratory search.

    We still begin from the first-cell residue structure,
    but instead of treating p as the principal unknown,
    derive the possible difference

        d = q-p.

    For fixed p-state:

        p = r10*k+a
        q = r20*l+b

    hence

        d = r20*l-r10*k+(b-a).

    The first-cell quotient equation gives

        n / R00 = k*l + E.

    Therefore

        l = (Q-E)/k.

    The experiment enumerates possible E values and asks
    whether the resulting d is integral.

    The goal is to measure whether this produces a sparse
    difference space.
    """

    Q = n // R00

    limit = math.isqrt(n)

    p_tested = 0
    states = 0
    d_candidates = 0
    exact = []

    start = time.perf_counter()

    for p in range(3, limit + 1, 2):

        p_tested += 1

        state = residue_state(
            n,
            p
        )

        if state is None:
            continue

        states += 1

        k, a, l_mod, b = state

        # ----------------------------------------------------
        # Conservative E range.
        #
        # We don't know l yet.
        #
        # Use a broad local bound. This is deliberately not
        # the old algebraic l solver.
        # ----------------------------------------------------

        # Approximate l around sqrt(n)/r20.
        l_est = math.isqrt(n) // r20

        # E = Q-k*l.
        #
        # Around the actual solution E is much smaller than Q.
        #
        # Search a symmetric local region whose width is
        # controlled by the carry scale.
        #
        # This is exploratory rather than an assumed theorem.
        E_center = max(
            0,
            Q - k * l_est
        )

        E_span = (
            r10
            + r20
            + 20
        )

        for E in range(
            max(0, E_center - E_span),
            E_center + E_span + 1
        ):

            numerator = Q - E

            if numerator <= 0:
                continue

            if numerator % k != 0:
                continue

            l = numerator // k

            if l <= 0:
                continue

            # Force residue condition.
            if l % r10 != l_mod:
                continue

            q = r20 * l + b

            d = q - p

            if d < 0:
                continue

            d_candidates += 1

            if p * q == n:

                exact.append(
                    (
                        p,
                        q,
                        d,
                        E
                    )
                )

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "p_tested": p_tested,
        "states": states,
        "d_candidates": d_candidates,
        "exact": exact,
        "time": elapsed,
    }


# ============================================================
# DIFFERENCE-ONLY FERMAT
# ============================================================

def fermat(n):

    a = math.isqrt(n)

    if a * a < n:
        a += 1

    iterations = 0

    start = time.perf_counter()

    while True:

        iterations += 1

        b2 = a * a - n
        b = math.isqrt(b2)

        if b * b == b2:

            p = a - b
            q = a + b

            if p > 1 and p * q == n:

                return {
                    "p": p,
                    "q": q,
                    "iterations": iterations,
                    "time": (
                        time.perf_counter()
                        - start
                    ),
                }

        a += 1


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    direct,
    diff,
    fermat_result,
):

    gap = true_q - true_p

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
    print("FACTOR GAP")

    print(
        "q-p =",
        gap
    )

    print(
        "sqrt(n) =",
        math.isqrt(n)
    )

    print()
    print("DIRECT")

    print(
        "p tested =",
        direct["tested"]
    )

    print(
        "runtime  = %.6f s"
        % direct["time"]
    )

    print()
    print("DIFFERENCE SEARCH")

    print(
        "p tested       =",
        diff["p_tested"]
    )

    print(
        "residue states =",
        diff["states"]
    )

    print(
        "d candidates   =",
        diff["d_candidates"]
    )

    print(
        "exact          =",
        len(diff["exact"])
    )

    print(
        "runtime        = %.6f s"
        % diff["time"]
    )

    print()
    print("FERMAT CONTROL")

    print(
        "iterations     =",
        fermat_result["iterations"]
    )

    print(
        "runtime        = %.6f s"
        % fermat_result["time"]
    )

    print(
        "p              =",
        fermat_result["p"]
    )

    print(
        "q              =",
        fermat_result["q"]
    )

    print()
    print("EXACT DIFFERENCE RESULTS")

    for p, q, d, E in diff["exact"]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"d={d} "
            f"E={E}"
        )

    print()
    print(
        "TRUE FOUND =",
        any(
            p == true_p
            and q == true_q
            for p, q, _, _ in diff["exact"]
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
    print("START EXPERIMENT 115")
    print("Difference constraint from the carry state")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)

    print(
        "R00 =",
        R00
    )

    print(
        "R11 =",
        R11
    )

    for bits in (
        30,
        36,
        42,
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

        direct = direct_factor(
            n
        )

        diff = difference_search(
            n
        )

        fermat_result = fermat(
            n
        )

        report(
            n,
            p,
            q,
            direct,
            diff,
            fermat_result
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 115")
    print("=" * 72)


if __name__ == "__main__":
    main()

