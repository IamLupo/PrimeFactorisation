#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 119
# Modular sum / discriminant search
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

R00 = R1[0] * R2[0]
R11 = R1[1] * R2[1]

M = math.lcm(R00, R11)

RANDOM_SEED = 119


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
# BUILD ALL POSSIBLE SUM RESIDUES
# ============================================================

def build_sum_residues(n):

    """
    We know

        pq == n (mod M).

    For every invertible p residue a:

        q == n*a^(-1) (mod M).

    Therefore

        S = p+q
          == a + n*a^(-1) (mod M).

    Many different a can produce the same S residue, so
    store only unique S residues while keeping the associated
    p/q residue pairs for diagnostics.
    """

    n_mod = n % M

    sums = {}
    pairs = 0

    for a in range(M):

        if math.gcd(a, M) != 1:
            continue

        b = (
            n_mod
            * pow(a, -1, M)
        ) % M

        s = (
            a + b
        ) % M

        pairs += 1

        sums.setdefault(
            s,
            []
        ).append(
            (a, b)
        )

    return sums, pairs


# ============================================================
# DISCRIMINANT TEST
# ============================================================

def check_sum(n, S):

    """
    For

        S = p+q

    we need

        D = S^2 - 4n

    to be a perfect square.

    Then

        p = (S - sqrt(D))/2
        q = (S + sqrt(D))/2.
    """

    D = S * S - 4 * n

    if D < 0:
        return None

    d = math.isqrt(D)

    if d * d != D:
        return None

    if (S - d) % 2 != 0:
        return None

    p = (
        S - d
    ) // 2

    q = (
        S + d
    ) // 2

    if p < 2:
        return None

    if p * q != n:
        return None

    return p, q, d


# ============================================================
# MODULAR-SUM SEARCH
# ============================================================

def modular_sum_search(
    n,
    sum_residues,
    max_gap_multiplier=1.0,
):
    """
    Search S values satisfying

        S == s (mod M).

    For balanced semiprimes we have

        S = 2*sqrt(n) + O((q-p)^2/sqrt(n)).

    The caller chooses the amount of S-space to inspect
    using max_gap_multiplier.

    Specifically we estimate a search range from the factor
    scale and inspect a configurable number of M-steps.
    """

    root = math.isqrt(n)

    if root * root < n:
        root += 1

    # --------------------------------------------------------
    # Start at ceil(2*sqrt(n)).
    # --------------------------------------------------------

    S_min = 2 * root

    # --------------------------------------------------------
    # We do NOT know the factor gap in advance.
    #
    # For the generated experimental semiprimes the factors
    # have equal bit length, so the sum remains on the sqrt(n)
    # scale.
    #
    # Search a configurable multiple of sqrt(n).
    # --------------------------------------------------------

    S_max = int(
        4 * root * max_gap_multiplier
    )

    S_residue_tests = 0
    S_values_tested = 0
    discriminants = 0

    results = []

    start = time.perf_counter()

    for s in sorted(sum_residues):

        # First S >= S_min with S == s mod M.

        delta = (
            s - S_min
        ) % M

        S = S_min + delta

        while S <= S_max:

            S_residue_tests += 1

            D = (
                S * S
                - 4 * n
            )

            discriminants += 1

            if D >= 0:

                d = math.isqrt(D)

                if d * d == D:

                    S_values_tested += 1

                    result = check_sum(
                        n,
                        S
                    )

                    if result is not None:

                        if result not in results:

                            results.append(
                                result
                            )

            S += M

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "sum_residue_count":
            len(sum_residues),

        "sum_residue_tests":
            S_residue_tests,

        "discriminants":
            discriminants,

        "square_hits":
            S_values_tested,

        "results":
            results,

        "runtime":
            elapsed,
    }


# ============================================================
# FERMAT BASELINE
# ============================================================

def fermat(n):

    root = math.isqrt(n)

    if root * root < n:
        root += 1

    a = root

    iterations = 0

    start = time.perf_counter()

    while True:

        iterations += 1

        b2 = (
            a * a
            - n
        )

        b = math.isqrt(
            b2
        )

        if b * b == b2:

            p = a - b
            q = a + b

            if p > 1 and p * q == n:

                return {
                    "p": p,
                    "q": q,
                    "iterations":
                        iterations,
                    "time":
                        time.perf_counter()
                        - start,
                }

        a += 1


# ============================================================
# MODULAR p SCAN BASELINE
# ============================================================

def modular_p_scan(n):

    limit = math.isqrt(n)

    tested = 0
    candidates = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        tested += 1

        if math.gcd(
            p,
            M
        ) != 1:

            continue

        q_res = (
            n
            * pow(
                p,
                -1,
                M
            )
        ) % M

        Q11 = n // R11

        q_low = (
            R11 * Q11
            + p - 1
        ) // p

        q_high = (
            R11 * (Q11 + 1)
            - 1
        ) // p

        q = (
            q_low
            + (q_res - q_low) % M
        )

        while q <= q_high:

            candidates += 1

            if (
                p <= q
                and
                p * q == n
            ):

                return {
                    "p": p,
                    "q": q,
                    "tested": tested,
                    "candidates":
                        candidates,
                    "time":
                        time.perf_counter()
                        - start,
                }

            q += M

    return {
        "p": None,
        "q": None,
        "tested": tested,
        "candidates":
            candidates,
        "time":
            time.perf_counter()
            - start,
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    sum_info,
    fermat_info,
    p_scan_info,
):

    gap = true_q - true_p

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
    print("TRUE SUM / GAP")

    print(
        "p+q =",
        true_p + true_q
    )

    print(
        "q-p =",
        gap
    )

    print(
        "sqrt(n) =",
        math.isqrt(n)
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
        "M =",
        M
    )

    print()
    print("SUM RESIDUE SPACE")

    print(
        "raw invertible p residues =",
        sum_info["raw_pairs"]
    )

    print(
        "unique sum residues       =",
        sum_info["unique_sums"]
    )

    print()
    print("MODULAR SUM SEARCH")

    print(
        "sum residues tested =",
        sum_info["sum_residue_count"]
    )

    print(
        "S values tested     =",
        sum_info["sum_residue_tests"]
    )

    print(
        "discriminants       =",
        sum_info["discriminants"]
    )

    print(
        "square hits         =",
        sum_info["square_hits"]
    )

    print(
        "results             =",
        len(sum_info["results"])
    )

    print(
        "runtime             = %.6f s"
        % sum_info["runtime"]
    )

    print()
    print("FERMAT")

    print(
        "iterations =",
        fermat_info["iterations"]
    )

    print(
        "runtime    = %.6f s"
        % fermat_info["time"]
    )

    print()
    print("MODULAR p SCAN")

    print(
        "p tested    =",
        p_scan_info["tested"]
    )

    print(
        "q candidates =",
        p_scan_info["candidates"]
    )

    print(
        "runtime      = %.6f s"
        % p_scan_info["time"]
    )

    print()
    print("SUM SEARCH RESULTS")

    for p, q, d in sum_info["results"]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"gap={d} "
            f"exact={p*q == n}"
        )

    print()
    print(
        "TRUE FACTORIZATION FOUND =",
        any(
            p == true_p
            and q == true_q
            for p, q, _ in sum_info["results"]
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
    print("START EXPERIMENT 119")
    print("Modular sum / discriminant search")
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

        # ----------------------------------------------------
        # Build modular sum states.
        # ----------------------------------------------------

        sums, raw_pairs = (
            build_sum_residues(n)
        )

        sum_info = {
            "unique_sums":
                len(sums),
            "raw_pairs":
                raw_pairs,
        }

        # ----------------------------------------------------
        # Search sums.
        #
        # This is deliberately bounded because our generator
        # creates equal-bit-length factors.
        # ----------------------------------------------------

        sum_search = modular_sum_search(
            n,
            sums,
            max_gap_multiplier=1.0,
        )

        sum_info.update(
            sum_search
        )

        # ----------------------------------------------------
        # Fermat control.
        # ----------------------------------------------------

        fermat_info = fermat(
            n
        )

        # ----------------------------------------------------
        # Experiment-113 style p scan.
        # ----------------------------------------------------

        p_scan_info = modular_p_scan(
            n
        )

        report(
            n,
            p,
            q,
            sum_info,
            fermat_info,
            p_scan_info,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 119")
    print("=" * 72)


if __name__ == "__main__":
    main()
