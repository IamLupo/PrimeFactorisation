#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 101
# Forced cross-scale carry collapse
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

CARRY_SLACK = 10


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
# Conservative E bound
# ------------------------------------------------------------

def E_valid(k, l, E):
    return 0 <= E <= k + l + CARRY_SLACK


# ------------------------------------------------------------
# Generate test case
# ------------------------------------------------------------

def make_test_case(bits=30):

    root = 1 << (bits // 2)

    p = root + random.randint(1000, 15000)
    q = root + random.randint(1000, 15000)

    return p * q, p, q


# ------------------------------------------------------------
# Main search
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
    # We only scan k0.
    #
    # l0 follows from:
    #
    #   E00 = Q00 - k0*l0
    #
    # and
    #
    #   0 <= E00 <= k0+l0+C
    # --------------------------------------------------------

    k0_max = int(math.isqrt(n) / r10) + 5000

    structures = 0
    residue_tests = 0

    range_reject = 0
    E_reject = 0
    cross_reject = 0
    carry_reject = 0

    survivors = []
    exact = []

    for k0 in range(1, k0_max + 1):

        # Conservative l0 interval.

        l0_lo = math.ceil(
            (Q00 - k0 - CARRY_SLACK)
            / (k0 + 1)
        )

        l0_hi = Q00 // k0

        if l0_lo < 1:
            l0_lo = 1

        if l0_lo > l0_hi:
            continue

        for l0 in range(l0_lo, l0_hi + 1):

            E00 = Q00 - k0 * l0

            if not E_valid(k0, l0, E00):
                E_reject += 1
                continue

            # ------------------------------------------------
            # ONLY 15 residue possibilities.
            # ------------------------------------------------

            for a0 in range(r10):
                for b0 in range(r20):

                    residue_tests += 1

                    # Force p and q.

                    p = r10 * k0 + a0
                    q = r20 * l0 + b0

                    # ------------------------------------------------
                    # Force the second-scale quotients.
                    # ------------------------------------------------

                    k1 = p // r11
                    l1 = q // r21

                    # Forced residues.

                    a1 = p % r11
                    b1 = q % r21

                    # ------------------------------------------------
                    # Compute all E values.
                    # ------------------------------------------------

                    E10 = Q10 - k1 * l0
                    E01 = Q01 - k0 * l1
                    E11 = Q11 - k1 * l1

                    # Negative E is impossible.

                    if E10 < 0 or E01 < 0 or E11 < 0:
                        cross_reject += 1
                        continue

                    # ------------------------------------------------
                    # E bounds.
                    # ------------------------------------------------

                    if not E_valid(k1, l0, E10):
                        E_reject += 1
                        continue

                    if not E_valid(k0, l1, E01):
                        E_reject += 1
                        continue

                    if not E_valid(k1, l1, E11):
                        E_reject += 1
                        continue

                    structures += 1

                    # ------------------------------------------------
                    # Full carry equations.
                    # ------------------------------------------------

                    _, _, _, got00 = carries(
                        r10,
                        r20,
                        k0,
                        l0,
                        a0,
                        b0,
                    )

                    if got00 != E00:
                        carry_reject += 1
                        continue

                    _, _, _, got10 = carries(
                        r11,
                        r20,
                        k1,
                        l0,
                        a1,
                        b0,
                    )

                    if got10 != E10:
                        carry_reject += 1
                        continue

                    _, _, _, got01 = carries(
                        r10,
                        r21,
                        k0,
                        l1,
                        a0,
                        b1,
                    )

                    if got01 != E01:
                        carry_reject += 1
                        continue

                    _, _, _, got11 = carries(
                        r11,
                        r21,
                        k1,
                        l1,
                        a1,
                        b1,
                    )

                    if got11 != E11:
                        carry_reject += 1
                        continue

                    # ------------------------------------------------
                    # Full carry-consistent survivor.
                    # ------------------------------------------------

                    item = {
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

                    survivors.append(item)

                    if p * q == n:
                        exact.append((p, q))

    return {
        "structures": structures,
        "residue_tests": residue_tests,
        "range_reject": range_reject,
        "E_reject": E_reject,
        "cross_reject": cross_reject,
        "carry_reject": carry_reject,
        "survivors": survivors,
        "exact": exact,
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    random.seed(101)

    n, true_p, true_q = make_test_case()

    print("=" * 72)
    print("START EXPERIMENT 101")
    print("Forced cross-scale carry collapse")
    print("=" * 72)

    print()
    print("n       =", n)
    print("true p  =", true_p)
    print("true q  =", true_q)
    print("r1      =", R1)
    print("r2      =", R2)

    # Ground truth

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
    # True E matrix
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
        "structures      =",
        result["structures"]
    )

    print(
        "residue tests   =",
        result["residue_tests"]
    )

    print(
        "range rejects   =",
        result["range_reject"]
    )

    print(
        "E rejects       =",
        result["E_reject"]
    )

    print(
        "cross rejects   =",
        result["cross_reject"]
    )

    print(
        "carry rejects   =",
        result["carry_reject"]
    )

    print(
        "survivors       =",
        len(result["survivors"])
    )

    print(
        "exact           =",
        len(result["exact"])
    )

    print(
        "runtime         = %.6f s"
        % (t1 - t0)
    )

    # --------------------------------------------------------
    # True structure
    # --------------------------------------------------------

    true_found = False

    for s in result["survivors"]:

        if (
            (s["p"] == true_p and s["q"] == true_q)
            or
            (s["p"] == true_q and s["q"] == true_p)
        ):
            true_found = True
            break

    print()
    print(
        "TRUE STRUCTURE SURVIVED =",
        true_found
    )

    print(
        "TRUE FACTORIZATION      =",
        true_p * true_q == n
    )

    # --------------------------------------------------------
    # Survivors
    # --------------------------------------------------------

    print()

    if result["survivors"]:

        print("SURVIVORS")

        for s in result["survivors"][:30]:

            print(
                "  "
                f"p={s['p']} "
                f"q={s['q']} | "
                f"k={s['k']} "
                f"l={s['l']} | "
                f"a={s['a']} "
                f"b={s['b']} | "
                f"E={s['E']}"
            )

    print()
    print(
        "=" * 72
    )
    print(
        "FINISHED EXPERIMENT 101"
    )
    print(
        "=" * 72
    )


if __name__ == "__main__":
    main()