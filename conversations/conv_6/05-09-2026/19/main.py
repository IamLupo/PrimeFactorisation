#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 102
# Algebraic elimination of l0
# ============================================================

R1 = [3, 29]
R2 = [5, 31]


# ------------------------------------------------------------
# Exact carry calculation
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# First-cell algebraic solver
# ------------------------------------------------------------

def solve_l0(Q, r1, r2, k, a, b):
    """
    Solve exactly:

        Q - k*l = carry_E(k,l,a,b)

    for l.

    For fixed k,a,b:

        c1 = floor(k*b/r2)

    is constant.

    Split l by residue modulo r1:

        l = t + r1*m

    where

        t in [0, r1-1].

    Then

        floor(l*a/r1)
        =
        floor(t*a/r1) + a*m

    and alpha = (l*a) mod r1 is constant.

    Therefore the entire equation becomes linear in m.

    Returns all exact integer solutions l.
    """

    solutions = []

    c1 = (k * b) // r2
    beta = (k * b) % r2

    for t in range(r1):

        # alpha is fixed inside this residue class.
        alpha = (t * a) % r1

        c3 = (
            r1 * beta
            + r2 * alpha
            + a * b
        ) // (r1 * r2)

        floor_t = (t * a) // r1

        # Equation:
        #
        # Q - k*(t+r1*m)
        #
        # =
        #
        # c1 + floor_t + a*m + c3
        #
        # => Q-k*t-c1-floor_t-c3
        #    =
        #    (k*r1+a)*m

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

        # Final exact verification.
        _, _, _, E_carry = carries(
            r1,
            r2,
            k,
            l,
            a,
            b,
        )

        E_quotient = Q - k * l

        if E_carry != E_quotient:
            continue

        solutions.append(l)

    return solutions


# ------------------------------------------------------------
# Full 2x2 search
# ------------------------------------------------------------

def search(n, r1s, r2s):

    r10, r11 = r1s
    r20, r21 = r2s

    R00 = r10 * r20
    R10 = r11 * r20
    R01 = r10 * r21
    R11 = r11 * r21

    Q00 = n // R00
    Q10 = n // R10
    Q01 = n // R01
    Q11 = n // R11

    # --------------------------------------------------------
    # We orient the search so p <= q.
    #
    # Then p <= sqrt(n), giving a finite tight k0 range.
    # --------------------------------------------------------

    p_limit = math.isqrt(n)

    k0_max = p_limit // r10

    k0_tested = 0
    first_cell_solutions = 0
    cross_tests = 0
    carry_tests = 0

    survivors = []
    exact = set()

    # --------------------------------------------------------
    # Scan only k0.
    # --------------------------------------------------------

    for k0 in range(1, k0_max + 1):

        k0_tested += 1

        # ----------------------------------------------------
        # Only 15 residue combinations for the first cell.
        # ----------------------------------------------------

        for a0 in range(r10):
            for b0 in range(r20):

                # ------------------------------------------------
                # Algebraically solve l0.
                # ------------------------------------------------

                l0_values = solve_l0(
                    Q00,
                    r10,
                    r20,
                    k0,
                    a0,
                    b0,
                )

                first_cell_solutions += len(l0_values)

                for l0 in l0_values:

                    # ------------------------------------------------
                    # Construct p and q immediately.
                    # ------------------------------------------------

                    p = r10 * k0 + a0
                    q = r20 * l0 + b0

                    # Orientation condition.
                    if p > q:
                        continue

                    # ------------------------------------------------
                    # Force second-scale quotients/residues.
                    # ------------------------------------------------

                    k1 = p // r11
                    a1 = p % r11

                    l1 = q // r21
                    b1 = q % r21

                    # ------------------------------------------------
                    # Compute all E values from n.
                    # ------------------------------------------------

                    E00 = Q00 - k0 * l0
                    E10 = Q10 - k1 * l0
                    E01 = Q01 - k0 * l1
                    E11 = Q11 - k1 * l1

                    if E00 < 0:
                        continue

                    if E10 < 0:
                        continue

                    if E01 < 0:
                        continue

                    if E11 < 0:
                        continue

                    cross_tests += 1

                    # ------------------------------------------------
                    # Full 2x2 carry check.
                    # ------------------------------------------------

                    _, _, _, got00 = carries(
                        r10,
                        r20,
                        k0,
                        l0,
                        a0,
                        b0,
                    )

                    carry_tests += 1

                    if got00 != E00:
                        continue

                    _, _, _, got10 = carries(
                        r11,
                        r20,
                        k1,
                        l0,
                        a1,
                        b0,
                    )

                    carry_tests += 1

                    if got10 != E10:
                        continue

                    _, _, _, got01 = carries(
                        r10,
                        r21,
                        k0,
                        l1,
                        a0,
                        b1,
                    )

                    carry_tests += 1

                    if got01 != E01:
                        continue

                    _, _, _, got11 = carries(
                        r11,
                        r21,
                        k1,
                        l1,
                        a1,
                        b1,
                    )

                    carry_tests += 1

                    if got11 != E11:
                        continue

                    # ------------------------------------------------
                    # Full carry-consistent survivor.
                    # ------------------------------------------------

                    item = (
                        p,
                        q,
                        k0,
                        k1,
                        l0,
                        l1,
                        a0,
                        a1,
                        b0,
                        b1,
                        E00,
                        E10,
                        E01,
                        E11,
                    )

                    survivors.append(item)

                    # Exact factorization.
                    if p * q == n:
                        exact.add((p, q))

    return {
        "k0_tested": k0_tested,
        "first_cell_solutions": first_cell_solutions,
        "cross_tests": cross_tests,
        "carry_tests": carry_tests,
        "survivors": survivors,
        "exact": sorted(exact),
    }


# ------------------------------------------------------------
# Test case
# ------------------------------------------------------------

def make_test_case(bits=30):

    root = 1 << (bits // 2)

    p = root + random.randint(1000, 15000)
    q = root + random.randint(1000, 15000)

    if p > q:
        p, q = q, p

    return p * q, p, q


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    random.seed(102)

    n, true_p, true_q = make_test_case()

    print("=" * 72)
    print("START EXPERIMENT 102")
    print("Algebraic elimination of l0")
    print("=" * 72)

    print()
    print("n       =", n)
    print("true p  =", true_p)
    print("true q  =", true_q)
    print("r1      =", R1)
    print("r2      =", R2)

    # --------------------------------------------------------
    # Ground truth
    # --------------------------------------------------------

    print()
    print("GROUND TRUTH")

    for i, r in enumerate(R1):

        k = true_p // r
        a = true_p % r

        print(
            f"r1[{i}]={r} -> k={k} a={a}"
        )

    for j, r in enumerate(R2):

        l = true_q // r
        b = true_q % r

        print(
            f"r2[{j}]={r} -> l={l} b={b}"
        )

    # --------------------------------------------------------
    # E matrix
    # --------------------------------------------------------

    print()
    print("TRUE E MATRIX")

    for r1 in R1:

        row = []

        for r2 in R2:

            k = true_p // r1
            l = true_q // r2

            E = n // (r1 * r2) - k * l

            row.append(E)

        print(row)

    # --------------------------------------------------------
    # Search
    # --------------------------------------------------------

    t0 = time.perf_counter()

    result = search(
        n,
        R1,
        R2,
    )

    t1 = time.perf_counter()

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------

    print()
    print("SEARCH RESULTS")

    print(
        "k0 tested             =",
        result["k0_tested"]
    )

    print(
        "first-cell solutions  =",
        result["first_cell_solutions"]
    )

    print(
        "cross tests           =",
        result["cross_tests"]
    )

    print(
        "carry tests           =",
        result["carry_tests"]
    )

    print(
        "survivors             =",
        len(result["survivors"])
    )

    print(
        "exact                 =",
        len(result["exact"])
    )

    print(
        "runtime               = %.6f s"
        % (t1 - t0)
    )

    # --------------------------------------------------------
    # True factorization
    # --------------------------------------------------------

    true_found = (
        (true_p, true_q)
        in result["exact"]
    )

    print()
    print(
        "TRUE FACTORIZATION FOUND =",
        true_found
    )

    # --------------------------------------------------------
    # Print survivors
    # --------------------------------------------------------

    print()

    for s in result["survivors"][:30]:

        (
            p,
            q,
            k0,
            k1,
            l0,
            l1,
            a0,
            a1,
            b0,
            b1,
            E00,
            E10,
            E01,
            E11,
        ) = s

        print(
            "  "
            f"p={p} q={q} | "
            f"k=({k0},{k1}) "
            f"l=({l0},{l1}) | "
            f"a=({a0},{a1}) "
            f"b=({b0},{b1}) | "
            f"E=({E00},{E10},{E01},{E11})"
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 102")
    print("=" * 72)


if __name__ == "__main__":
    main()