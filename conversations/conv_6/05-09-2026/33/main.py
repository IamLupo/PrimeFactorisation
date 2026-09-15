#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 116
# Exact 2x2 E-matrix closure
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20
R10 = r11 * r20
R01 = r10 * r21
R11 = r11 * r21

RANDOM_SEED = 116


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
# GROUND TRUTH E MATRIX
# ============================================================

def true_E_matrix(n, p, q):

    out = []

    for r1 in R1:

        row = []

        k = p // r1

        for r2 in R2:

            l = q // r2

            E = (
                n // (r1 * r2)
                - k * l
            )

            row.append(E)

        out.append(row)

    return out


# ============================================================
# EXACT RANK-1 TEST
# ============================================================

def rank1(E00, E10, E01, E11, n):

    """
    Let

        Aij = Qij - Eij = ki*lj.

    Exact multiplicative rank-1 requires

        A00*A11 = A01*A10.
    """

    Q00 = n // R00
    Q10 = n // R10
    Q01 = n // R01
    Q11 = n // R11

    A00 = Q00 - E00
    A10 = Q10 - E10
    A01 = Q01 - E01
    A11 = Q11 - E11

    if min(
        A00,
        A10,
        A01,
        A11
    ) <= 0:
        return False

    return (
        A00 * A11
        ==
        A10 * A01
    )


# ============================================================
# BUILD CANDIDATE K/L VALUES
# ============================================================

def recover_kl(
    n,
    E00,
    E10,
    E01,
    E11
):

    Q00 = n // R00
    Q10 = n // R10
    Q01 = n // R01
    Q11 = n // R11

    A00 = Q00 - E00
    A10 = Q10 - E10
    A01 = Q01 - E01
    A11 = Q11 - E11

    if min(
        A00,
        A10,
        A01,
        A11
    ) <= 0:
        return None

    if (
        A00 * A11
        !=
        A10 * A01
    ):
        return None

    # --------------------------------------------------------
    # Factor A00 = k0*l0.
    #
    # For this experiment we enumerate divisors of A00.
    # This is intentionally a diagnostic experiment:
    # the important question is whether the exact E-matrix
    # contains enough information to recover k/l.
    # --------------------------------------------------------

    limit = math.isqrt(A00)

    for d in range(1, limit + 1):

        if A00 % d != 0:
            continue

        k0 = d
        l0 = A00 // d

        # Derive k1 and l1 from the other cells.

        if A10 % l0 != 0:
            continue

        if A01 % k0 != 0:
            continue

        k1 = A10 // l0
        l1 = A01 // k0

        if k1 * l1 != A11:
            continue

        yield (
            k0,
            k1,
            l0,
            l1
        )


# ============================================================
# E-STATE ENUMERATION
# ============================================================

def E_bounds(n):

    """
    Conservative carry bounds.

    E is roughly the sum of the three carry components.
    For this experiment we use a broad but finite range
    based on the quotient magnitudes.
    """

    Q00 = n // R00
    Q10 = n // R10
    Q01 = n // R01
    Q11 = n // R11

    # We don't enumerate to Q itself.
    #
    # Use the fact that the carry is dominated by the quotient
    # variables. The scale factors make the following safe
    # for our experimental range.

    # Return broad ranges to be tightened from the actual
    # residue states below.

    return (
        Q00,
        Q10,
        Q01,
        Q11
    )


# ============================================================
# RESIDUE STATE -> E VALUES
# ============================================================

def state_E(
    r1,
    r2,
    k,
    l,
    a,
    b
):

    return carries(
        r1,
        r2,
        k,
        l,
        a,
        b
    )[3]


# ============================================================
# DIRECT FINITE RESIDUE EXPERIMENT
# ============================================================

def search(n):

    Q00 = n // R00
    Q10 = n // R10
    Q01 = n // R01
    Q11 = n // R11

    counters = {
        "residue_states": 0,
        "rank1_tests": 0,
        "rank1_pass": 0,
        "kl_recoveries": 0,
        "residue_matches": 0,
        "exact": 0,
    }

    survivors = []

    # --------------------------------------------------------
    # We use the residue state at BOTH scales.
    #
    # a0 in [0,2]
    # b0 in [0,4]
    # a1 in [0,28]
    # b1 in [0,30]
    #
    # Total raw state space:
    #
    # 3*5*29*31 = 13485.
    #
    # This is deliberately finite.
    # --------------------------------------------------------

    start = time.perf_counter()

    for a0 in range(r10):

        for b0 in range(r20):

            # ------------------------------------------------
            # k0 and l0 are still unknown here.
            #
            # We derive possible E00 carry forms symbolically.
            # ------------------------------------------------

            for a1 in range(r11):

                for b1 in range(r21):

                    counters["residue_states"] += 1

                    # ------------------------------------------------
                    # A consistency condition between residues.
                    #
                    # Since the same p is represented at both r1
                    # scales:
                    #
                    # r10*k0+a0 = r11*k1+a1.
                    #
                    # Therefore
                    #
                    # 17*k0 + a0
                    # =
                    # 29*k1 + a1.
                    #
                    # Likewise q.
                    #
                    # We solve these congruences.
                    # ------------------------------------------------

                    # 17*k0 + a0 == 29*k1 + a1
                    #
                    # mod 29:
                    #
                    # 17*k0 == a1-a0 mod 29.

                    inv17 = pow(
                        r10,
                        -1,
                        r11
                    )

                    k0_res = (
                        (a1 - a0)
                        * inv17
                    ) % r11

                    # 5*l0+b0 == 31*l1+b1
                    #
                    # mod 31:
                    #
                    # 5*l0 == b1-b0 mod31.

                    inv5 = pow(
                        r20,
                        -1,
                        r21
                    )

                    l0_res = (
                        (b1 - b0)
                        * inv5
                    ) % r21

                    # ------------------------------------------------
                    # At this stage:
                    #
                    # k0 = k0_res mod 29
                    # l0 = l0_res mod 31
                    #
                    # Therefore:
                    #
                    # k0 = k0_res + 29*x
                    # l0 = l0_res + 31*y
                    #
                    # Use A00 = k0*l0 and E00 = Q00-A00.
                    #
                    # The key experiment is whether the rank-1
                    # equations can determine x,y without a full
                    # sqrt(n) scan.
                    # ------------------------------------------------

                    # Instead of enumerating x,y indefinitely,
                    # derive approximate k0 and l0 from sqrt(Q00).

                    root = math.isqrt(
                        Q00
                    )

                    # Test a narrow window around the balanced
                    # scale, but deliberately record failures.
                    #
                    # This is a diagnostic step, not an assumed
                    # theorem.

                    x0 = max(
                        0,
                        (
                            root
                            - k0_res
                        ) // r11
                        - 100
                    )

                    x1 = (
                        root
                        - k0_res
                    ) // r11 + 100

                    y0 = max(
                        0,
                        (
                            root
                            - l0_res
                        ) // r21
                        - 100
                    )

                    y1 = (
                        root
                        - l0_res
                    ) // r21 + 100

                    for x in range(
                        x0,
                        x1 + 1
                    ):

                        k0 = (
                            k0_res
                            + r11 * x
                        )

                        if k0 <= 0:
                            continue

                        for y in range(
                            y0,
                            y1 + 1
                        ):

                            l0 = (
                                l0_res
                                + r21 * y
                            )

                            if l0 <= 0:
                                continue

                            E00 = (
                                Q00
                                - k0 * l0
                            )

                            if E00 < 0:
                                continue

                            # ------------------------------------------------
                            # Derive p and q.
                            # ------------------------------------------------

                            p = (
                                r10 * k0
                                + a0
                            )

                            q = (
                                r20 * l0
                                + b0
                            )

                            # Verify second-scale residues.
                            if p % r11 != a1:
                                continue

                            if q % r21 != b1:
                                continue

                            k1 = p // r11
                            l1 = q // r21

                            E10 = (
                                Q10
                                - k1 * l0
                            )

                            E01 = (
                                Q01
                                - k0 * l1
                            )

                            E11 = (
                                Q11
                                - k1 * l1
                            )

                            if min(
                                E10,
                                E01,
                                E11
                            ) < 0:
                                continue

                            counters[
                                "rank1_tests"
                            ] += 1

                            if not rank1(
                                n,
                                E00,
                                E10,
                                E01,
                                E11
                            ):
                                continue

                            counters[
                                "rank1_pass"
                            ] += 1

                            # ------------------------------------------------
                            # Exact residue carry equations.
                            # ------------------------------------------------

                            got00 = state_E(
                                r10,
                                r20,
                                k0,
                                l0,
                                a0,
                                b0
                            )

                            got10 = state_E(
                                r11,
                                r20,
                                k1,
                                l0,
                                a1,
                                b0
                            )

                            got01 = state_E(
                                r10,
                                r21,
                                k0,
                                l1,
                                a0,
                                b1
                            )

                            got11 = state_E(
                                r11,
                                r21,
                                k1,
                                l1,
                                a1,
                                b1
                            )

                            if (
                                got00 != E00
                                or got10 != E10
                                or got01 != E01
                                or got11 != E11
                            ):
                                continue

                            counters[
                                "residue_matches"
                            ] += 1

                            item = {
                                "p": p,
                                "q": q,
                                "k": (
                                    k0,
                                    k1
                                ),
                                "l": (
                                    l0,
                                    l1
                                ),
                                "a": (
                                    a0,
                                    a1
                                ),
                                "b": (
                                    b0,
                                    b1
                                ),
                                "E": (
                                    E00,
                                    E10,
                                    E01,
                                    E11
                                ),
                            }

                            survivors.append(
                                item
                            )

                            if p * q == n:
                                counters[
                                    "exact"
                                ] += 1

    counters["runtime"] = (
        time.perf_counter()
        - start
    )

    return counters, survivors


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    counters,
    survivors
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

    print(
        "r1 =",
        R1
    )

    print(
        "r2 =",
        R2
    )

    print()
    print("SEARCH RESULTS")

    for key, value in counters.items():

        if key == "runtime":

            print(
                "runtime = %.6f s"
                % value
            )

        else:

            print(
                f"{key:<22} = {value}"
            )

    print()
    print("SURVIVORS")

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
            "  ..."
            f"{len(survivors)-25}"
            " more"
        )

    print()
    print(
        "TRUE FOUND =",
        any(
            s["p"] == true_p
            and
            s["q"] == true_q
            for s in survivors
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
    print("START EXPERIMENT 116")
    print("Exact 2x2 E-matrix closure")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)
    print("R10 =", R10)
    print("R01 =", R01)
    print("R11 =", R11)

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

        counters, survivors = search(
            n
        )

        report(
            n,
            p,
            q,
            counters,
            survivors
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 116")
    print("=" * 72)


if __name__ == "__main__":
    main()
