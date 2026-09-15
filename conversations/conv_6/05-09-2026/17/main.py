#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 100
# Carry-aware 2x2 transition collapse
# ============================================================

CARRY_SLACK = 10

R1 = [3, 29]
R2 = [5, 31]


# ------------------------------------------------------------
# Mathematical helpers
# ------------------------------------------------------------

def carries(r1, r2, k, l, a, b):
    """
    Exact carry decomposition.

        p = r1*k + a
        q = r2*l + b

        E = floor(p*q/(r1*r2)) - k*l

    with

        c1 = floor(k*b/r2)
        c2 = floor(l*a/r1)
        c3 = floor((r1*beta + r2*alpha + a*b)/(r1*r2))

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


def exact_E(n, r1, r2, k, l):
    R = r1 * r2
    Q = n // R
    return Q - k * l


def valid_E_range(k, l):
    """
    Conservative carry bound:

        E <= k + l + CARRY_SLACK

    and E >= 0.
    """
    lo = 0
    hi = k + l + CARRY_SLACK
    return lo, hi


# ------------------------------------------------------------
# Generate a balanced test factorization
# ------------------------------------------------------------

def make_test_case(bits=30):
    """
    Produce a semiprime with factors near sqrt(n).

    This keeps the experiment focused on the regime used by
    the previous experiments.
    """

    root = 1 << (bits // 2)

    # Slightly perturb around sqrt(n)
    p = root + random.randint(1000, 15000)
    q = root + random.randint(1000, 15000)

    n = p * q

    return n, p, q


# ------------------------------------------------------------
# Build quotient-transition structures
# ------------------------------------------------------------

def transition_structures(n, r1s, r2s):
    """
    Construct candidate

        (k0, k1, l0, l1)

    structures using only quotient/carry bounds.

    We deliberately do NOT use p*q=n here.
    """

    r10, r11 = r1s
    r20, r21 = r2s

    Q00 = n // (r10 * r20)
    Q10 = n // (r11 * r20)
    Q01 = n // (r10 * r21)
    Q11 = n // (r11 * r21)

    # Balanced-factor regime:
    # k0 is roughly sqrt(n)/r10.
    #
    # Give ourselves some room above that estimate.
    k0_max = int(math.isqrt(n) / r10) + 5000

    structures = []

    for k0 in range(1, k0_max + 1):

        # From
        #
        # 0 <= Q00-k0*l0 <= k0+l0+C
        #
        # obtain the conservative l0 interval.

        l0_lo = math.ceil(
            (Q00 - k0 - CARRY_SLACK) / (k0 + 1)
        )

        l0_hi = Q00 // k0

        if l0_lo < 1:
            l0_lo = 1

        if l0_lo > l0_hi:
            continue

        for l0 in range(l0_lo, l0_hi + 1):

            E00 = Q00 - k0 * l0

            loE, hiE = valid_E_range(k0, l0)

            if not (loE <= E00 <= hiE):
                continue

            # ------------------------------------------------
            # Infer k1 interval from:
            #
            # 0 <= Q10-k1*l0 <= k1+l0+C
            # ------------------------------------------------

            k1_lo = math.ceil(
                (Q10 - l0 - CARRY_SLACK) / (l0 + 1)
            )

            k1_hi = Q10 // l0

            if k1_lo < 1:
                k1_lo = 1

            if k1_lo > k1_hi:
                continue

            for k1 in range(k1_lo, k1_hi + 1):

                E10 = Q10 - k1 * l0

                loE10, hiE10 = valid_E_range(k1, l0)

                if not (loE10 <= E10 <= hiE10):
                    continue

                # ------------------------------------------------
                # Infer l1 interval from:
                #
                # 0 <= Q01-k0*l1 <= k0+l1+C
                # ------------------------------------------------

                l1_lo = math.ceil(
                    (Q01 - k0 - CARRY_SLACK) / (k0 + 1)
                )

                l1_hi = Q01 // k0

                if l1_lo < 1:
                    l1_lo = 1

                if l1_lo > l1_hi:
                    continue

                for l1 in range(l1_lo, l1_hi + 1):

                    E01 = Q01 - k0 * l1

                    loE01, hiE01 = valid_E_range(k0, l1)

                    if not (loE01 <= E01 <= hiE01):
                        continue

                    # Final cell
                    E11 = Q11 - k1 * l1

                    if E11 < 0:
                        continue

                    loE11, hiE11 = valid_E_range(k1, l1)

                    if not (loE11 <= E11 <= hiE11):
                        continue

                    structures.append(
                        (k0, k1, l0, l1,
                         E00, E10, E01, E11)
                    )

    return structures


# ------------------------------------------------------------
# Carry-aware residue collapse
# ------------------------------------------------------------

def residue_collapse(
    n,
    r1s,
    r2s,
    structures
):
    """
    For every transition structure:

      1. enumerate only (a0,b0) for the cheapest cell
      2. construct p,q
      3. force all remaining residues
      4. check residue ranges
      5. check every carry equation
      6. only THEN test p*q == n
    """

    r10, r11 = r1s
    r20, r21 = r2s

    survivors = []
    exact = []

    residue_attempts = 0
    range_rejections = 0
    carry_rejections = 0

    for (
        k0, k1,
        l0, l1,
        E00, E10,
        E01, E11
    ) in structures:

        # Cheapest possible residue search.
        #
        # r10=3, r20=5 -> only 15 possibilities.
        for a0 in range(r10):
            for b0 in range(r20):

                residue_attempts += 1

                p = r10 * k0 + a0
                q = r20 * l0 + b0

                # ------------------------------------------------
                # Force the remaining residues.
                # ------------------------------------------------

                a1 = p - r11 * k1
                b1 = q - r21 * l1

                if not (0 <= a1 < r11):
                    range_rejections += 1
                    continue

                if not (0 <= b1 < r21):
                    range_rejections += 1
                    continue

                # ------------------------------------------------
                # Check the full 2x2 carry system.
                # ------------------------------------------------

                _, _, _, got_E00 = carries(
                    r10, r20,
                    k0, l0,
                    a0, b0
                )

                if got_E00 != E00:
                    carry_rejections += 1
                    continue

                _, _, _, got_E10 = carries(
                    r11, r20,
                    k1, l0,
                    a1, b0
                )

                if got_E10 != E10:
                    carry_rejections += 1
                    continue

                _, _, _, got_E01 = carries(
                    r10, r21,
                    k0, l1,
                    a0, b1
                )

                if got_E01 != E01:
                    carry_rejections += 1
                    continue

                _, _, _, got_E11 = carries(
                    r11, r21,
                    k1, l1,
                    a1, b1
                )

                if got_E11 != E11:
                    carry_rejections += 1
                    continue

                survivors.append(
                    {
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
                )

                if p * q == n:
                    exact.append((p, q))

    return (
        survivors,
        exact,
        residue_attempts,
        range_rejections,
        carry_rejections,
    )


# ------------------------------------------------------------
# Pretty printer
# ------------------------------------------------------------

def print_survivors(survivors, exact, limit=20):

    print()
    print("SURVIVORS:", len(survivors))
    print("EXACT    :", len(exact))

    if survivors:
        print()
        print("First survivors:")

        for s in survivors[:limit]:
            print(
                "  p=%d q=%d | "
                "k=%s l=%s | "
                "a=%s b=%s | "
                "E=%s"
                % (
                    s["p"],
                    s["q"],
                    s["k"],
                    s["l"],
                    s["a"],
                    s["b"],
                    s["E"],
                )
            )

    if exact:
        print()
        print("EXACT FACTORIZATIONS:")

        for p, q in exact:
            print(
                "  %d * %d = %d"
                % (p, q, p * q)
            )


# ============================================================
# Main
# ============================================================

def main():

    random.seed(100)

    # --------------------------------------------------------
    # Generate test semiprime
    # --------------------------------------------------------

    n, true_p, true_q = make_test_case(bits=30)

    print("=" * 72)
    print("START EXPERIMENT 100")
    print("Carry-aware 2x2 transition collapse")
    print("=" * 72)

    print()
    print("n       =", n)
    print("true p  =", true_p)
    print("true q  =", true_q)
    print("r1      =", R1)
    print("r2      =", R2)

    # --------------------------------------------------------
    # Ground-truth structural values
    # --------------------------------------------------------

    print()
    print("GROUND TRUTH")

    for i, r1 in enumerate(R1):
        k = true_p // r1
        a = true_p % r1

        print(
            "r1[%d]=%d -> k=%d a=%d"
            % (i, r1, k, a)
        )

    for j, r2 in enumerate(R2):
        l = true_q // r2
        b = true_q % r2

        print(
            "r2[%d]=%d -> l=%d b=%d"
            % (j, r2, l, b)
        )

    print()
    print("TRUE E MATRIX")

    for r1 in R1:
        row = []

        for r2 in R2:
            k = true_p // r1
            l = true_q // r2

            E = exact_E(
                n,
                r1,
                r2,
                k,
                l
            )

            row.append(E)

        print(row)

    # --------------------------------------------------------
    # Transition stage
    # --------------------------------------------------------

    t0 = time.perf_counter()

    structures = transition_structures(
        n,
        R1,
        R2
    )

    t1 = time.perf_counter()

    print()
    print("TRANSITION STAGE")
    print("structures =", len(structures))
    print(
        "time       = %.6f s"
        % (t1 - t0)
    )

    # --------------------------------------------------------
    # Carry-aware residue collapse
    # --------------------------------------------------------

    t2 = time.perf_counter()

    (
        survivors,
        exact,
        residue_attempts,
        range_rejections,
        carry_rejections,
    ) = residue_collapse(
        n,
        R1,
        R2,
        structures
    )

    t3 = time.perf_counter()

    print()
    print("CARRY-AWARE RESIDUE STAGE")
    print("residue attempts   =", residue_attempts)
    print("range rejections   =", range_rejections)
    print("carry rejections   =", carry_rejections)
    print("survivors          =", len(survivors))
    print("exact              =", len(exact))
    print(
        "time               = %.6f s"
        % (t3 - t2)
    )

    # --------------------------------------------------------
    # Check whether the true pair survived
    # --------------------------------------------------------

    true_found = False

    for s in survivors:
        if (
            (s["p"] == true_p and s["q"] == true_q)
            or
            (s["p"] == true_q and s["q"] == true_p)
        ):
            true_found = True
            break

    print()
    print("TRUE STRUCTURE SURVIVED =", true_found)
    print("TRUE FACTORIZATION      =", (true_p * true_q == n))

    print_survivors(
        survivors,
        exact
    )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 100")
    print("=" * 72)


if __name__ == "__main__":
    main()
