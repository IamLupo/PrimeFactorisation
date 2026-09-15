#!/usr/bin/env python3

import math
import random


# ============================================================
# START EXPERIMENT 111
# Test whether H is identically n
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20

RANDOM_SEED = 111


# ============================================================
# PRIME GENERATION
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
# BUILD H
# ============================================================

def build_H(n, a, s):

    Q = n // R00

    # p residue represented by:
    #
    # p = R00*u + C
    #
    # where
    #
    # C = 17*s + a

    C = r10 * s + a

    p_res = C % R00

    if math.gcd(p_res, R00) != 1:
        return None

    # --------------------------------------------------------
    # Force q modulo R00 from
    #
    # p*q == n (mod R00)
    # --------------------------------------------------------

    q0 = (
        n
        * pow(
            p_res,
            -1,
            R00
        )
    ) % R00

    # q0 = 19*t + b
    t = q0 // r20
    b = q0 % r20

    # --------------------------------------------------------
    # Fixed residue carry pieces.
    # --------------------------------------------------------

    beta = (
        s * b
    ) % r20

    c1_const = (
        s * b
    ) // r20

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
    # Linearized first-cell equation.
    # --------------------------------------------------------

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
        "p_res": p_res,
        "q0": q0,
        "t": t,
        "b": b,
        "A": A,
        "B": B,
        "H": H,
    }


# ============================================================
# ANALYZE ALL STATES
# ============================================================

def analyze(n):

    states = []

    count = 0

    H_equal_n = 0
    H_divides_n = 0
    H_multiple_n = 0

    differences = {}

    gcd_values = {}

    for a in range(r10):

        for s in range(r20):

            info = build_H(
                n,
                a,
                s
            )

            if info is None:
                continue

            count += 1

            H = info["H"]

            delta = H - n

            differences[delta] = (
                differences.get(delta, 0) + 1
            )

            if H == n:
                H_equal_n += 1

            if H != 0 and n % H == 0:
                H_divides_n += 1

            if H != 0 and H % n == 0:
                H_multiple_n += 1

            g = math.gcd(
                n,
                abs(H)
            )

            gcd_values[g] = (
                gcd_values.get(g, 0) + 1
            )

            states.append(
                (
                    a,
                    s,
                    info
                )
            )

    return {
        "states": states,
        "count": count,
        "H_equal_n": H_equal_n,
        "H_divides_n": H_divides_n,
        "H_multiple_n": H_multiple_n,
        "differences": differences,
        "gcd_values": gcd_values,
    }


# ============================================================
# REPORT
# ============================================================

def report(n, p, q, result):

    print()
    print("-" * 72)

    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    print()
    print("GRID")
    print("r1  =", R1)
    print("r2  =", R2)
    print("R00 =", R00)

    print()
    print("STATE COUNT")

    print(
        "admissible states =",
        result["count"]
    )

    print()
    print("H RELATION TO n")

    print(
        "H == n             =",
        result["H_equal_n"]
    )

    print(
        "H divides n        =",
        result["H_divides_n"]
    )

    print(
        "H is multiple of n =",
        result["H_multiple_n"]
    )

    print()
    print("H - n VALUES")

    sorted_diffs = sorted(
        result["differences"].items(),
        key=lambda x: (
            abs(x[0]),
            x[0]
        )
    )

    for delta, count in sorted_diffs[:30]:

        print(
            f"  H-n={delta} "
            f"count={count}"
        )

    if len(sorted_diffs) > 30:

        print(
            "  ..."
            f"{len(sorted_diffs) - 30}"
            " more distinct values"
        )

    print()
    print("GCD VALUES")

    gcd_items = sorted(
        result["gcd_values"].items(),
        key=lambda x: x[1],
        reverse=True
    )

    for g, count in gcd_items[:20]:

        print(
            f"  gcd={g} "
            f"count={count}"
        )

    print()
    print("SAMPLE STATES")

    for a, s, info in result["states"][:20]:

        print(
            "  "
            f"a={a} "
            f"s={s} | "
            f"C={info['C']} | "
            f"q0={info['q0']} | "
            f"t={info['t']} "
            f"b={info['b']} | "
            f"A={info['A']} "
            f"B={info['B']} | "
            f"H={info['H']} | "
            f"H-n={info['H']-n}"
        )

    # --------------------------------------------------------
    # True state.
    # --------------------------------------------------------

    true_k = p // r10
    true_a = p % r10
    true_s = true_k % r20

    true_info = build_H(
        n,
        true_a,
        true_s
    )

    print()
    print("TRUE STATE")

    print(
        "a =",
        true_a
    )

    print(
        "s =",
        true_s
    )

    if true_info is not None:

        print(
            "H =",
            true_info["H"]
        )

        print(
            "H-n =",
            true_info["H"] - n
        )

        print(
            "H == n =",
            true_info["H"] == n
        )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 111")
    print("Test whether H is identically n")
    print("=" * 72)

    print()

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

        result = analyze(n)

        report(
            n,
            p,
            q,
            result
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 111")
    print("=" * 72)


if __name__ == "__main__":
    main()
