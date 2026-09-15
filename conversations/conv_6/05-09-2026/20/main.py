#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 103
# Semiprime carry recovery
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

RANDOM_SEED = 103

# Test several semiprime sizes.
# Set this to fewer sizes if you want a quick run.
TEST_BITS = [30, 36, 42]


# ============================================================
# PRIME GENERATION
# ============================================================

def is_probable_prime(n):
    """
    Deterministic Miller-Rabin for the integer sizes used here.
    The bases below are sufficient for our experimental range.
    """

    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    )

    for p in small_primes:
        if n % p == 0:
            return n == p

    # n-1 = d * 2^s
    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    # Good deterministic set for the sizes used here.
    bases = (
        2, 3, 5, 7, 11,
        13, 17
    )

    for a in bases:

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        composite = True

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True


def random_prime(bits):
    """
    Generate a random odd prime with exactly `bits` bits.
    """

    while True:

        p = random.getrandbits(bits)

        # Force top bit.
        p |= (1 << (bits - 1))

        # Force odd.
        p |= 1

        if is_probable_prime(p):
            return p


def make_semiprime(bits):
    """
    Generate p < q, both prime, with approximately bits/2
    bits each.

    n therefore has approximately `bits` bits.
    """

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
# CARRY EQUATIONS
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
# ALGEBRAIC SOLVER FOR l0
# ============================================================

def solve_l0(Q, r1, r2, k, a, b):
    """
    Solve exactly

        Q - k*l = E(k,l,a,b)

    by splitting

        l = t + r1*m

    for t = 0,...,r1-1.
    """

    solutions = []

    c1 = (k * b) // r2
    beta = (k * b) % r2

    for t in range(r1):

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
            continue

        if numerator % denominator != 0:
            continue

        m = numerator // denominator

        if m < 0:
            continue

        l = t + r1 * m

        if l <= 0:
            continue

        # Final exact check of the algebraic solution.
        _, _, _, E1 = carries(
            r1, r2,
            k, l,
            a, b
        )

        E2 = Q - k * l

        if E1 != E2:
            continue

        solutions.append(l)

    return solutions


# ============================================================
# FULL 2x2 SEARCH
# ============================================================

def search(n):

    r10, r11 = R1
    r20, r21 = R2

    Q00 = n // (r10 * r20)
    Q10 = n // (r11 * r20)
    Q01 = n // (r10 * r21)
    Q11 = n // (r11 * r21)

    # p <= q by construction.
    p_limit = math.isqrt(n)

    k0_max = p_limit // r10

    counters = {
        "k0_tested": 0,
        "residue_tests": 0,
        "first_cell_solutions": 0,
        "cross_tests": 0,
        "carry_tests": 0,
        "survivors": 0,
    }

    survivors = []

    # --------------------------------------------------------
    # Scan k0 only.
    # --------------------------------------------------------

    for k0 in range(1, k0_max + 1):

        counters["k0_tested"] += 1

        # ----------------------------------------------------
        # The first cell has only r10*r20 residue combinations.
        # Here: 3*5 = 15.
        # ----------------------------------------------------

        for a0 in range(r10):

            for b0 in range(r20):

                counters["residue_tests"] += 1

                # Algebraically recover l0.
                l0_values = solve_l0(
                    Q00,
                    r10,
                    r20,
                    k0,
                    a0,
                    b0,
                )

                counters["first_cell_solutions"] += len(
                    l0_values
                )

                for l0 in l0_values:

                    # ------------------------------------------------
                    # First-scale values define p and q.
                    # ------------------------------------------------

                    p = r10 * k0 + a0
                    q = r20 * l0 + b0

                    # Search only p <= q.
                    if p > q:
                        continue

                    # ------------------------------------------------
                    # Second scale is completely forced.
                    # ------------------------------------------------

                    k1 = p // r11
                    a1 = p % r11

                    l1 = q // r21
                    b1 = q % r21

                    # ------------------------------------------------
                    # Derive E values from n.
                    # ------------------------------------------------

                    E00 = Q00 - k0 * l0
                    E10 = Q10 - k1 * l0
                    E01 = Q01 - k0 * l1
                    E11 = Q11 - k1 * l1

                    if (
                        E00 < 0
                        or E10 < 0
                        or E01 < 0
                        or E11 < 0
                    ):
                        continue

                    counters["cross_tests"] += 1

                    # ------------------------------------------------
                    # Full exact carry system.
                    # ------------------------------------------------

                    _, _, _, got00 = carries(
                        r10,
                        r20,
                        k0,
                        l0,
                        a0,
                        b0
                    )

                    counters["carry_tests"] += 1

                    if got00 != E00:
                        continue

                    _, _, _, got10 = carries(
                        r11,
                        r20,
                        k1,
                        l0,
                        a1,
                        b0
                    )

                    counters["carry_tests"] += 1

                    if got10 != E10:
                        continue

                    _, _, _, got01 = carries(
                        r10,
                        r21,
                        k0,
                        l1,
                        a0,
                        b1
                    )

                    counters["carry_tests"] += 1

                    if got01 != E01:
                        continue

                    _, _, _, got11 = carries(
                        r11,
                        r21,
                        k1,
                        l1,
                        a1,
                        b1
                    )

                    counters["carry_tests"] += 1

                    if got11 != E11:
                        continue

                    # ------------------------------------------------
                    # Carry-consistent candidate.
                    #
                    # NOTE:
                    # We do NOT reject with p*q != n here.
                    # ------------------------------------------------

                    survivor = {
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

                    survivors.append(survivor)

                    counters["survivors"] += 1

    return counters, survivors


# ============================================================
# REPORT
# ============================================================

def print_case(n, p, q, counters, survivors):

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())

    print()
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)

    print()
    print("FACTOR CHECK")
    print("p*q == n =", p * q == n)

    print()
    print("GROUND TRUTH")

    for i, r in enumerate(R1):

        k = p // r
        a = p % r

        print(
            f"r1[{i}]={r} -> "
            f"k={k} a={a}"
        )

    for j, r in enumerate(R2):

        l = q // r
        b = q % r

        print(
            f"r2[{j}]={r} -> "
            f"l={l} b={b}"
        )

    print()
    print("TRUE E MATRIX")

    for r1 in R1:

        row = []

        for r2 in R2:

            k = p // r1
            l = q // r2

            E = n // (r1 * r2) - k * l

            row.append(E)

        print(row)

    print()
    print("SEARCH RESULTS")

    for key, value in counters.items():

        print(
            f"{key:<22} = {value}"
        )

    # --------------------------------------------------------
    # Exact verification AFTER the carry search.
    # --------------------------------------------------------

    exact = []

    for s in survivors:

        if s["p"] * s["q"] == n:
            exact.append(
                (s["p"], s["q"])
            )

    print()
    print("POST-SEARCH VERIFICATION")

    print(
        "exact factorizations =",
        len(exact)
    )

    true_found = (
        (p, q) in exact
    )

    print(
        "true factorization found =",
        true_found
    )

    print()
    print("CARRY SURVIVORS")

    for s in survivors[:25]:

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

    if len(survivors) > 25:
        print(
            f"  ... {len(survivors) - 25} more"
        )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 103")
    print("Semiprime carry recovery")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("test sizes =", TEST_BITS)

    total_start = time.perf_counter()

    for bits in TEST_BITS:

        print()
        print("=" * 72)
        print(f"GENERATING {bits}-BIT SEMIPRIME")
        print("=" * 72)

        gen_start = time.perf_counter()

        n, p, q = make_semiprime(bits)

        gen_time = time.perf_counter() - gen_start

        print()
        print(
            "prime generation time = %.6f s"
            % gen_time
        )

        search_start = time.perf_counter()

        counters, survivors = search(n)

        search_time = (
            time.perf_counter()
            - search_start
        )

        print_case(
            n,
            p,
            q,
            counters,
            survivors
        )

        print()
        print(
            "search runtime = %.6f s"
            % search_time
        )

    total_time = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 103")
    print("=" * 72)

    print(
        "total runtime = %.6f s"
        % total_time
    )


if __name__ == "__main__":
    main()
