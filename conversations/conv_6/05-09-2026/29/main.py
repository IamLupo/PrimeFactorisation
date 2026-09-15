#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 112
# Two-cell early quotient filter
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20
R11 = r11 * r21

RANDOM_SEED = 112


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
# FIRST-CELL EXACT SOLVER
# ============================================================

def solve_l0_fixed_t(
    Q,
    r1,
    r2,
    k,
    a,
    b,
    t,
):

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

    # Exact first-cell verification.

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
# FORCE q MOD 15
# ============================================================

def force_q_mod_R00(n, p):

    if math.gcd(p, R00) != 1:
        return None

    return (
        n
        * pow(
            p,
            -1,
            R00
        )
    ) % R00


def decode_q0(q0):

    # q0 = 5*t + b
    #
    # t = l0 mod 3
    # b = q mod 5

    t = q0 // r20
    b = q0 % r20

    return t, b


# ============================================================
# SEARCH
# ============================================================

def search(n):

    Q00 = n // R00
    Q11 = n // R11

    p_limit = math.isqrt(n)

    counters = {
        "p_tested": 0,
        "modular_forced": 0,
        "modular_failed": 0,

        "second_cell_bucket_tests": 0,
        "second_cell_bucket_pass": 0,

        "first_cell_solver_tests": 0,
        "first_cell_candidates": 0,

        "final_candidates": 0,
        "exact": 0,
    }

    second_bucket = []
    first_candidates = []
    final_candidates = []

    start = time.perf_counter()

    # --------------------------------------------------------
    # Only odd p for odd semiprimes.
    # --------------------------------------------------------

    for p in range(3, p_limit + 1, 2):

        counters["p_tested"] += 1

        if math.gcd(p, R00) != 1:

            counters["modular_failed"] += 1

            continue

        # ----------------------------------------------------
        # First modular constraint.
        # ----------------------------------------------------

        q0 = force_q_mod_R00(
            n,
            p
        )

        if q0 is None:
            counters["modular_failed"] += 1
            continue

        counters["modular_forced"] += 1

        t, b0 = decode_q0(q0)

        k0 = p // r10
        a0 = p % r10

        # ----------------------------------------------------
        # ====================================================
        # SECOND-CELL EARLY FILTER
        # ====================================================
        #
        # Instead of solving l0 first, derive the possible
        # q values from the second quotient bucket.
        #
        # q must satisfy
        #
        # floor(p*q / R11) = floor(n / R11)
        #
        # so
        #
        # R11*Q11 <= p*q < R11*(Q11+1)
        #
        # giving
        #
        # ceil(R11*Q11/p)
        # <= q <=
        # floor((R11*(Q11+1)-1)/p).
        #
        # We only need q ≡ q0 (mod 5) here.
        # ----------------------------------------------------

        counters["second_cell_bucket_tests"] += 1

        q_lo = (
            R11 * Q11 + p - 1
        ) // p

        q_hi = (
            R11 * (Q11 + 1) - 1
        ) // p

        # Restrict q to:
        #
        # q = 5*l + b0
        #
        # Find first q >= q_lo with this residue.

        first_q = q_lo

        rem = (
            first_q - b0
        ) % r20

        if rem:
            first_q += r20 - rem

        if first_q > q_hi:

            continue

        # There can be more than one q in this tiny interval.
        # Enumerate only this arithmetic progression.

        q = first_q

        second_passed = False

        while q <= q_hi:

            # Must satisfy the first-cell modular relation.

            if q % R00 == q0:

                second_passed = True

                second_bucket.append(
                    (p, q)
                )

            q += r20

        if not second_passed:
            continue

        counters["second_cell_bucket_pass"] += 1

        # ----------------------------------------------------
        # Now solve the first-cell equation.
        # ----------------------------------------------------

        l_candidates_here = []

        for _, q_candidate in [
            x for x in second_bucket
            if x[0] == p
        ]:

            l0 = (
                q_candidate - b0
            ) // r20

            counters["first_cell_solver_tests"] += 1

            # Verify first-cell carry equation directly.

            _, _, _, Ecarry = carries(
                r10,
                r20,
                k0,
                l0,
                a0,
                b0
            )

            Equot = Q00 - k0 * l0

            if Ecarry != Equot:
                continue

            l_candidates_here.append(
                (
                    l0,
                    q_candidate
                )
            )

        if not l_candidates_here:
            continue

        counters["first_cell_candidates"] += len(
            l_candidates_here
        )

        # ----------------------------------------------------
        # Full structural verification.
        # ----------------------------------------------------

        for l0, q in l_candidates_here:

            k1 = p // r11
            a1 = p % r11

            l1 = q // r21
            b1 = q % r21

            E00 = Q00 - k0 * l0
            E10 = (
                n // (r11 * r20)
                - k1 * l0
            )
            E01 = (
                n // (r10 * r21)
                - k0 * l1
            )
            E11 = Q11 - k1 * l1

            # ------------------------------------------------
            # Full four-cell carry check.
            # ------------------------------------------------

            if carries(
                r10, r20,
                k0, l0,
                a0, b0
            )[3] != E00:
                continue

            if carries(
                r11, r20,
                k1, l0,
                a1, b0
            )[3] != E10:
                continue

            if carries(
                r10, r21,
                k0, l1,
                a0, b1
            )[3] != E01:
                continue

            if carries(
                r11, r21,
                k1, l1,
                a1, b1
            )[3] != E11:
                continue

            candidate = {
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

            final_candidates.append(candidate)

            if p * q == n:
                counters["exact"] += 1

    counters["final_candidates"] = len(
        final_candidates
    )

    counters["runtime"] = (
        time.perf_counter() - start
    )

    return (
        counters,
        second_bucket,
        final_candidates,
    )


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    counters,
    second_bucket,
    final_candidates,
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

    print(
        "R00 =",
        R00
    )

    print(
        "R11 =",
        R11
    )

    print()
    print("SEARCH RESULTS")

    print(
        "p tested                 =",
        counters["p_tested"]
    )

    print(
        "modular forced           =",
        counters["modular_forced"]
    )

    print(
        "modular failed           =",
        counters["modular_failed"]
    )

    print(
        "second-cell bucket tests =",
        counters["second_cell_bucket_tests"]
    )

    print(
        "second-cell bucket pass  =",
        counters["second_cell_bucket_pass"]
    )

    print(
        "first-cell solver tests  =",
        counters["first_cell_solver_tests"]
    )

    print(
        "first-cell candidates    =",
        counters["first_cell_candidates"]
    )

    print(
        "final candidates         =",
        counters["final_candidates"]
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
    # Second-cell survivors.
    # --------------------------------------------------------

    print()
    print("SECOND-CELL BUCKET SURVIVORS")

    for p, q in second_bucket[:20]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"error={p*q-n}"
        )

    if len(second_bucket) > 20:

        print(
            "  ..."
            f"{len(second_bucket)-20}"
            " more"
        )

    # --------------------------------------------------------
    # Final candidates.
    # --------------------------------------------------------

    print()
    print("FINAL CANDIDATES")

    for s in final_candidates[:20]:

        print(
            "  "
            f"p={s['p']} "
            f"q={s['q']} | "
            f"k={s['k']} "
            f"l={s['l']} | "
            f"a={s['a']} "
            f"b={s['b']} | "
            f"E={s['E']} | "
            f"exact={s['p']*s['q']==n}"
        )

    if len(final_candidates) > 20:

        print(
            "  ..."
            f"{len(final_candidates)-20}"
            " more"
        )

    print()
    print(
        "TRUE FACTORIZATION FOUND =",
        any(
            s["p"] == true_p
            and s["q"] == true_q
            for s in final_candidates
        )
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 112")
    print("Two-cell early quotient filter")
    print("=" * 72)

    print()

    for bits in (
        30,
        36,
        42,
        48,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        counters, second_bucket, final_candidates = search(n)

        report(
            n,
            p,
            q,
            counters,
            second_bucket,
            final_candidates,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 112")
    print("=" * 72)


if __name__ == "__main__":
    main()