#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 110
# Corrected divisibility recovery
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20

RANDOM_SEED = 110


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
# BUILD H
# ============================================================

def build_H(n, a, s):

    Q = n // R00

    # p = 323*u + C
    C = r10 * s + a

    # p mod 323
    p_res = C % R00

    if math.gcd(p_res, R00) != 1:
        return None

    # q = n * p^(-1) mod 323
    q0 = (
        n
        * pow(p_res, -1, R00)
    ) % R00

    # q0 = 19*t + b
    t = q0 // r20
    b = q0 % r20

    # k = 19*u + s

    # floor(k*b/19)
    beta = (s * b) % r20
    c1_const = (s * b) // r20

    # floor(t*a/17)
    alpha = (t * a) % r10
    c2_const = (t * a) // r10

    # c3 is fixed for the residue class
    c3 = (
        r10 * beta
        + r20 * alpha
        + a * b
    ) // R00

    # First-cell equation:
    #
    # A - B*u = p*m
    #
    # where
    #
    # p = 323*u + C

    A = (
        Q
        - s * t
        - c1_const
        - c2_const
        - c3
    )

    B = (
        r20 * t + b
    )

    H = (
        R00 * A
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
# VERIFY RESIDUE CLASS
# ============================================================

def verify_class(n, info):

    H = info["H"]

    # gcd may be n itself.
    g = math.gcd(n, abs(H))

    candidates = set()

    if g > 1:
        candidates.add(g)

    if g != 0 and n % g == 0:
        candidates.add(n // g)

    verified = []

    for p in candidates:

        # Need p to belong to this residue class.

        C = info["C"]

        if p < C:
            continue

        if (p - C) % R00 != 0:
            continue

        u = (p - C) // R00

        k = r20 * u + info["s"]

        a = info["a"]
        b = info["b"]
        t = info["t"]

        # ----------------------------------------------------
        # Recover m from:
        #
        # A - B*u = p*m
        # ----------------------------------------------------

        numerator = (
            info["A"]
            - info["B"] * u
        )

        if numerator < 0:
            continue

        if numerator % p != 0:
            continue

        m = numerator // p

        l = t + r10 * m

        if l <= 0:
            continue

        # Recover q.

        q = r20 * l + b

        # Basic structural checks.

        if p > q:
            continue

        # ----------------------------------------------------
        # Exact first-cell carry validation.
        # ----------------------------------------------------

        _, _, _, Ec = carries(
            r10,
            r20,
            k,
            l,
            a,
            b
        )

        Eq = (
            n // R00
            - k * l
        )

        if Ec != Eq:
            continue

        verified.append(
            {
                "p": p,
                "q": q,
                "k": k,
                "l": l,
                "a": a,
                "b": b,
                "t": t,
                "m": m,
                "E": Eq,
                "gcd": g,
                "H": H,
            }
        )

    return g, verified


# ============================================================
# SEARCH
# ============================================================

def search(n):

    counters = {
        "classes": 0,
        "gcd_1": 0,
        "gcd_nontrivial": 0,
        "gcd_n": 0,
        "verified": 0,
        "exact": 0,
    }

    results = []

    start = time.perf_counter()

    for a in range(r10):

        for s in range(r20):

            info = build_H(
                n,
                a,
                s
            )

            if info is None:
                continue

            counters["classes"] += 1

            g, verified = verify_class(
                n,
                info
            )

            # ------------------------------------------------
            # Classify gcd.
            # ------------------------------------------------

            if g == 1:
                counters["gcd_1"] += 1

            elif g == n:
                counters["gcd_n"] += 1

            elif 1 < g < n:
                counters["gcd_nontrivial"] += 1

            # ------------------------------------------------
            # Record verified structural solutions.
            # ------------------------------------------------

            for result in verified:

                counters["verified"] += 1

                result["a_class"] = a
                result["s_class"] = s

                result["H_equals_n"] = (
                    result["H"] == n
                )

                if (
                    result["p"] * result["q"]
                    == n
                ):
                    counters["exact"] += 1

                results.append(result)

    elapsed = (
        time.perf_counter()
        - start
    )

    return counters, results, elapsed


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    counters,
    results,
    elapsed
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

    print("r1  =", R1)
    print("r2  =", R2)
    print("R00 =", R00)

    print()
    print("GCD CLASSIFICATION")

    print(
        "classes tested       =",
        counters["classes"]
    )

    print(
        "gcd = 1              =",
        counters["gcd_1"]
    )

    print(
        "1 < gcd < n          =",
        counters["gcd_nontrivial"]
    )

    print(
        "gcd = n              =",
        counters["gcd_n"]
    )

    print()
    print("RECOVERY")

    print(
        "verified candidates  =",
        counters["verified"]
    )

    print(
        "exact factorizations =",
        counters["exact"]
    )

    print(
        "runtime              = %.9f s"
        % elapsed
    )

    print()
    print("VERIFIED RESULTS")

    for r in results[:30]:

        print(
            "  "
            f"a={r['a_class']} "
            f"s={r['s_class']} | "
            f"p={r['p']} "
            f"q={r['q']} | "
            f"gcd={r['gcd']} | "
            f"H=n:{r['H_equals_n']} | "
            f"exact:{r['p'] * r['q'] == n}"
        )

    if len(results) > 30:

        print(
            "  ..."
            f"{len(results) - 30}"
            " more"
        )

    # --------------------------------------------------------
    # Locate true class.
    # --------------------------------------------------------

    true_k = true_p // r10
    true_a = true_p % r10

    true_s = true_k % r20

    true_info = build_H(
        n,
        true_a,
        true_s
    )

    print()
    print("TRUE RESIDUE CLASS")

    print(
        "a =", true_a
    )

    print(
        "s =", true_s
    )

    if true_info is None:

        print(
            "true class unexpectedly non-invertible"
        )

    else:

        print(
            "C =",
            true_info["C"]
        )

        print(
            "q0 =",
            true_info["q0"]
        )

        print(
            "t =",
            true_info["t"]
        )

        print(
            "b =",
            true_info["b"]
        )

        print(
            "H =",
            true_info["H"]
        )

        print(
            "H == n =",
            true_info["H"] == n
        )

    print()
    print(
        "TRUE FACTORIZATION FOUND =",
        any(
            (
                r["p"] == true_p
                and
                r["q"] == true_q
            )
            for r in results
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 110")
    print("Corrected divisibility recovery")
    print("=" * 72)

    print()
    print(
        "r1 =",
        R1
    )

    print(
        "r2 =",
        R2
    )

    print(
        "R00 =",
        R00
    )

    # --------------------------------------------------------
    # Run several sizes.
    # --------------------------------------------------------

    for bits in (
        30,
        36,
        42,
        48,
        54,
        60,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        counters, results, elapsed = search(n)

        report(
            n,
            p,
            q,
            counters,
            results,
            elapsed
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 110")
    print("=" * 72)


if __name__ == "__main__":
    main()
