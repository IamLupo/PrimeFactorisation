#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 108
# Algebraic k-congruence search
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20

RANDOM_SEED = 108


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

        p |= (1 << (bits - 1))
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
# FIRST-CELL ALGEBRA
# ============================================================

def solve_l0_fixed_t(Q, r1, r2, k, a, b, t):

    c1 = (k * b) // r2
    beta = (k * b) % r2

    alpha = (t * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    floor_t = (t * a) // r1

    numerator = (
        Q
        - k * t
        - c1
        - floor_t
        - c3
    )

    denominator = k * r1 + a

    if numerator < 0:
        return None

    if numerator % denominator != 0:
        return None

    m = numerator // denominator

    if m < 0:
        return None

    l = t + r1 * m

    if l <= 0:
        return None

    # Exact validation.

    _, _, _, Ecarry = carries(
        r1,
        r2,
        k,
        l,
        a,
        b
    )

    Equot = Q - k * l

    if Ecarry != Equot:
        return None

    return l


# ============================================================
# DIRECT BASELINE
# ============================================================

def baseline_scan(n):

    Q00 = n // R00

    limit = math.isqrt(n)

    tested = 0
    candidates = 0
    exact = []

    start = time.perf_counter()

    for p in range(3, limit + 1, 2):

        tested += 1

        if math.gcd(p, R00) != 1:
            continue

        k = p // r10
        a = p % r10

        q0 = (
            n * pow(p, -1, R00)
        ) % R00

        b = q0 % r20
        t = q0 // r20

        l = solve_l0_fixed_t(
            Q00,
            r10,
            r20,
            k,
            a,
            b,
            t
        )

        if l is None:
            continue

        candidates += 1

        q = r20 * l + b

        if p > q:
            continue

        if p * q == n:
            exact.append((p, q))

    elapsed = (
        time.perf_counter() - start
    )

    return {
        "tested": tested,
        "candidates": candidates,
        "exact": exact,
        "time": elapsed,
    }


# ============================================================
# CONGRUENCE CLASS ANALYSIS
# ============================================================

def derive_k_classes(n):
    """
    Analyze k modulo R00.

    Since

        p = 17k + a

    and a in [0,16],

    every k/a pair defines p.

    We determine which k residue classes can possibly satisfy
    the modular inverse relation.

    This is intentionally an exploratory stage:
    it does NOT assume a closed-form solution exists.

    Returns allowed (k_mod, a) pairs.
    """

    allowed = []

    for a in range(r10):

        for kres in range(R00):

            p_res = (
                r10 * kres + a
            ) % R00

            # p must be invertible modulo R00.
            if math.gcd(
                p_res,
                R00
            ) != 1:
                continue

            # Corresponding q residue.

            q_res = (
                n
                * pow(p_res, -1, R00)
            ) % R00

            # Decode q residue.

            b = q_res % r20
            t = q_res // r20

            # We don't know the large k yet, but the residue
            # pair is arithmetically admissible.

            allowed.append(
                (
                    kres,
                    a,
                    b,
                    t
                )
            )

    return allowed


# ============================================================
# SEARCH BY RESIDUE CLASSES
# ============================================================

def residue_class_search(n, allowed):

    Q00 = n // R00

    p_limit = math.isqrt(n)

    # Maximum k.
    k_max = p_limit // r10

    tested_k = 0
    generated_k = 0
    l_candidates = 0

    solutions = []

    start = time.perf_counter()

    for (
        kres,
        a,
        b,
        t
    ) in allowed:

        # k = kres + m*R00

        first = kres

        if first == 0:
            first = R00

        if first > k_max:
            continue

        for k in range(
            first,
            k_max + 1,
            R00
        ):

            tested_k += 1

            p = r10 * k + a

            if p > p_limit:
                continue

            generated_k += 1

            l = solve_l0_fixed_t(
                Q00,
                r10,
                r20,
                k,
                a,
                b,
                t
            )

            if l is None:
                continue

            l_candidates += 1

            q = r20 * l + b

            if p > q:
                continue

            product = p * q

            if product // R00 != n // R00:
                continue

            if product % R00 != n % R00:
                continue

            solutions.append(
                (p, q)
            )

    elapsed = (
        time.perf_counter() - start
    )

    return {
        "tested_k": tested_k,
        "generated_k": generated_k,
        "l_candidates": l_candidates,
        "solutions": solutions,
        "time": elapsed,
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    p,
    q,
    baseline,
    allowed,
    modular
):

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())

    print()
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    print()
    print("GRID")

    print("r1  =", R1)
    print("r2  =", R2)

    print("R00 =", R00)

    print()
    print("BASELINE DIRECT-p SEARCH")

    print(
        "p tested       =",
        baseline["tested"]
    )

    print(
        "first-cell candidates =",
        baseline["candidates"]
    )

    print(
        "exact          =",
        len(baseline["exact"])
    )

    print(
        "runtime        = %.6f s"
        % baseline["time"]
    )

    print()
    print("CONGRUENCE CLASS SPACE")

    print(
        "allowed (k mod R00, a) states =",
        len(allowed)
    )

    print()
    print("RESIDUE-CLASS SEARCH")

    print(
        "k values tested =",
        modular["tested_k"]
    )

    print(
        "generated k     =",
        modular["generated_k"]
    )

    print(
        "l candidates    =",
        modular["l_candidates"]
    )

    print(
        "solutions       =",
        len(modular["solutions"])
    )

    print(
        "runtime         = %.6f s"
        % modular["time"]
    )

    print()
    print("SOLUTIONS")

    for x in modular["solutions"][:20]:

        pp, qq = x

        print(
            "  "
            f"p={pp} "
            f"q={qq} "
            f"exact={pp * qq == n}"
        )

    print()
    print("TRUE FACTORIZATION")

    print(
        "found =",
        (p, q) in modular["solutions"]
    )

    if baseline["time"] > 0:

        print()
        print("SPEEDUP")

        print(
            "baseline / residue-class = %.3fx"
            % (
                baseline["time"]
                / max(
                    modular["time"],
                    1e-12
                )
            )
        )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 108")
    print("Algebraic k-congruence search")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)

    for bits in (30, 36, 42):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        # ----------------------------------------------------
        # Baseline.
        # ----------------------------------------------------

        baseline = baseline_scan(n)

        # ----------------------------------------------------
        # Build modular state space.
        # ----------------------------------------------------

        allowed = derive_k_classes(n)

        # ----------------------------------------------------
        # Search each class.
        # ----------------------------------------------------

        modular = residue_class_search(
            n,
            allowed
        )

        report(
            n,
            p,
            q,
            baseline,
            allowed,
            modular
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 108")
    print("=" * 72)


if __name__ == "__main__":
    main()

