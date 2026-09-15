#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 113
# Corrected bidirectional product reconstruction
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20       # 15
R11 = r11 * r21       # 899
M = math.lcm(R00, R11)  # 13485

RANDOM_SEED = 113


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
    factor_bits = bits // 2

    while True:
        p = random_prime(factor_bits)
        q = random_prime(factor_bits)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() >= bits - 1:
            return n, p, q


# ============================================================
# EXACT CARRY DECOMPOSITION
# ============================================================

def carries(r1, r2, k, l, a, b):
    """
    p = r1*k + a
    q = r2*l + b

    E = floor(p*q/(r1*r2)) - k*l

    and

        E = c1 + c2 + c3
    """

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
# GENERAL LINEAR CONGRUENCE
# ============================================================

def solve_linear_congruence(a, b, m):
    """
    Solve

        a*x = b (mod m)

    and return all x modulo m.
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


# ============================================================
# q RESIDUES FOR FIXED p
# ============================================================

def q_residues_for_p(n, p, modulus):
    """
    Solve

        p*q = n (mod modulus).
    """

    return solve_linear_congruence(
        p % modulus,
        n % modulus,
        modulus
    )


# ============================================================
# FIRST CELL STATE
# ============================================================

def first_cell_states(n, p):
    """
    For the first modulus:

        q = 5*l0 + b0
        q mod 15 is known.
    """

    states = []

    q0_values = q_residues_for_p(
        n,
        p,
        R00
    )

    if not q0_values:
        return states

    k0 = p // r10
    a0 = p % r10

    for q0 in q0_values:

        # q0 = 5*t + b
        t0 = q0 // r20
        b0 = q0 % r20

        states.append(
            {
                "q0": q0,
                "t0": t0,
                "b0": b0,
                "k0": k0,
                "a0": a0,
            }
        )

    return states


# ============================================================
# CRT FOR q
# ============================================================

def combine_q_congruences(
    q1,
    m1,
    q2,
    m2
):
    """
    Combine

        q = q1 (mod m1)
        q = q2 (mod m2)

    even when m1 and m2 are not coprime.
    """

    g = math.gcd(m1, m2)

    if (q2 - q1) % g != 0:
        return []

    m1r = m1 // g
    m2r = m2 // g

    diff = (q2 - q1) // g

    if m2r == 1:
        t = 0
    else:
        inv = pow(
            m1r % m2r,
            -1,
            m2r
        )

        t = (
            diff * inv
        ) % m2r

    x = q1 + m1 * t

    lcm = m1 * m2r

    return [x % lcm]


# ============================================================
# SECOND-CELL QUOTIENT BUCKET
# ============================================================

def second_bucket_q_interval(n, p):
    """
    floor(p*q / R11) = floor(n / R11)

    implies

        R11*Q11 <= p*q < R11*(Q11+1).

    """

    Q11 = n // R11

    q_low = (
        R11 * Q11 + p - 1
    ) // p

    q_high = (
        R11 * (Q11 + 1) - 1
    ) // p

    return q_low, q_high


# ============================================================
# FULL 2x2 CARRY CHECK
# ============================================================

def full_check(n, p, q, state):
    """
    Validate the complete 2x2 carry system.
    """

    if p > q:
        return None

    k0 = state["k0"]
    a0 = state["a0"]

    l0 = q // r20
    b0 = q % r20

    # --------------------------------------------------------
    # Cell (0,0)
    # --------------------------------------------------------

    Q00 = n // R00

    E00 = (
        Q00
        - k0 * l0
    )

    got00 = carries(
        r10,
        r20,
        k0,
        l0,
        a0,
        b0
    )[3]

    if got00 != E00:
        return None

    # --------------------------------------------------------
    # Second scale.
    # --------------------------------------------------------

    k1 = p // r11
    a1 = p % r11

    l1 = q // r21
    b1 = q % r21

    # --------------------------------------------------------
    # E values from n.
    # --------------------------------------------------------

    E10 = (
        n // (r11 * r20)
        - k1 * l0
    )

    E01 = (
        n // (r10 * r21)
        - k0 * l1
    )

    E11 = (
        n // R11
        - k1 * l1
    )

    # --------------------------------------------------------
    # Cell (1,0)
    # --------------------------------------------------------

    got10 = carries(
        r11,
        r20,
        k1,
        l0,
        a1,
        b0
    )[3]

    if got10 != E10:
        return None

    # --------------------------------------------------------
    # Cell (0,1)
    # --------------------------------------------------------

    got01 = carries(
        r10,
        r21,
        k0,
        l1,
        a0,
        b1
    )[3]

    if got01 != E01:
        return None

    # --------------------------------------------------------
    # Cell (1,1)
    # --------------------------------------------------------

    got11 = carries(
        r11,
        r21,
        k1,
        l1,
        a1,
        b1
    )[3]

    if got11 != E11:
        return None

    return {
        "p": p,
        "q": q,

        "k": (k0, k1),
        "l": (l0, l1),

        "a": (a0, a1),
        "b": (b0, b1),

        "E": (
            E00,
            E10,
            E01,
            E11,
        ),
    }


# ============================================================
# MAIN P-DIRECTION SEARCH
# ============================================================

def search_p_direction(n):

    limit = math.isqrt(n)

    counters = {
        "p_tested": 0,
        "p_with_R00_solution": 0,
        "p_with_R11_solution": 0,

        "combined_residues": 0,
        "q_bucket_candidates": 0,
        "carry_candidates": 0,
        "exact": 0,
    }

    bucket = []
    survivors = []
    exact = []

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):
        counters["p_tested"] += 1

        # ----------------------------------------------------
        # First congruence
        # ----------------------------------------------------

        q1_values = q_residues_for_p(
            n,
            p,
            R00
        )

        if q1_values:
            counters["p_with_R00_solution"] += 1

        # ----------------------------------------------------
        # Second congruence
        # ----------------------------------------------------

        q2_values = q_residues_for_p(
            n,
            p,
            R11
        )

        if q2_values:
            counters["p_with_R11_solution"] += 1

        if not q1_values or not q2_values:
            continue

        # State information from first cell.

        states = first_cell_states(
            n,
            p
        )

        for state in states:

            q1 = state["q0"]

            for q2 in q2_values:

                combined = combine_q_congruences(
                    q1,
                    R00,
                    q2,
                    R11
                )

                counters["combined_residues"] += len(
                    combined
                )

                if not combined:
                    continue

                q_low, q_high = (
                    second_bucket_q_interval(
                        n,
                        p
                    )
                )

                for q_res in combined:

                    # ------------------------------------------------
                    # First q >= q_low with q == q_res (mod M)
                    # ------------------------------------------------

                    delta = (
                        q_res - q_low
                    ) % M

                    q = q_low + delta

                    while q <= q_high:

                        counters[
                            "q_bucket_candidates"
                        ] += 1

                        bucket.append(
                            (
                                p,
                                q,
                                p * q - n
                            )
                        )

                        checked = full_check(
                            n,
                            p,
                            q,
                            state
                        )

                        if checked is not None:

                            counters[
                                "carry_candidates"
                            ] += 1

                            survivors.append(
                                checked
                            )

                            if (
                                checked["p"]
                                * checked["q"]
                                == n
                            ):

                                counters["exact"] += 1

                                exact.append(
                                    (
                                        checked["p"],
                                        checked["q"]
                                    )
                                )

                        q += M

    counters["runtime"] = (
        time.perf_counter()
        - start
    )

    return (
        counters,
        bucket,
        survivors,
        exact
    )


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    counters,
    bucket,
    survivors,
    exact
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
    print("MODULI")

    print(
        "R00 =",
        R00
    )

    print(
        "R11 =",
        R11
    )

    print(
        "M   =",
        M
    )

    print()
    print("SEARCH RESULTS")

    print(
        "p tested                 =",
        counters["p_tested"]
    )

    print(
        "p with R00 solution      =",
        counters["p_with_R00_solution"]
    )

    print(
        "p with R11 solution      =",
        counters["p_with_R11_solution"]
    )

    print(
        "combined q residues      =",
        counters["combined_residues"]
    )

    print(
        "q bucket candidates      =",
        counters["q_bucket_candidates"]
    )

    print(
        "carry survivors          =",
        counters["carry_candidates"]
    )

    print(
        "exact                    =",
        counters["exact"]
    )

    print(
        "runtime                  = %.6f s"
        % counters["runtime"]
    )

    # --------------------------------------------------------
    # Bucket error
    # --------------------------------------------------------

    if bucket:

        max_error = max(
            abs(error)
            for _, _, error in bucket
        )

        print()
        print("BUCKET ERROR")

        print(
            "max |pq-n| =",
            max_error
        )

        print(
            "R11-1      =",
            R11 - 1
        )

    # --------------------------------------------------------
    # Carry survivors
    # --------------------------------------------------------

    print()
    print("CARRY SURVIVORS")

    for s in survivors[:20]:

        print(
            "  "
            f"p={s['p']} "
            f"q={s['q']} | "
            f"k={s['k']} "
            f"l={s['l']} | "
            f"a={s['a']} "
            f"b={s['b']} | "
            f"E={s['E']} | "
            f"exact={s['p'] * s['q'] == n}"
        )

    if len(survivors) > 20:

        print(
            "  ..."
            f"{len(survivors) - 20}"
            " more"
        )

    # --------------------------------------------------------
    # Exact
    # --------------------------------------------------------

    print()
    print("EXACT")

    for p, q in exact:

        print(
            "  "
            f"{p} * {q} = {p*q}"
        )

    print()
    print(
        "TRUE FACTORIZATION FOUND =",
        (
            true_p,
            true_q
        ) in exact
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 113")
    print("Corrected bidirectional product reconstruction")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)
    print("R11 =", R11)
    print("M   =", M)

    for bits in (
        30,
        36,
        42,
        48,
        54,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        (
            counters,
            bucket,
            survivors,
            exact
        ) = search_p_direction(n)

        report(
            n,
            p,
            q,
            counters,
            bucket,
            survivors,
            exact
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 113")
    print("=" * 72)


if __name__ == "__main__":
    main()
