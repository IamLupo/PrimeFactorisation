#!/usr/bin/env python3

# ============================================================
# START EXPERIMENT 167
# Global Factor Candidate -> C-Grid Compatibility
# ============================================================

from __future__ import annotations

import math
import random
import time

from collections import defaultdict


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BIT_SIZES = [30, 36, 42, 48, 54]
SAMPLES_PER_SIZE = 1

SEED = 167


# ------------------------------------------------------------
# Primality
# ------------------------------------------------------------

def is_probable_prime(n: int) -> bool:

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37,
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in (
        2, 3, 5, 7, 11, 13, 17
    ):

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = x * x % n

            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(
    bits: int,
    rng: random.Random,
) -> int:

    while True:

        x = rng.getrandbits(bits)

        x |= 1 << (bits - 1)
        x |= 1

        if is_probable_prime(x):
            return x


def make_semiprime(
    bits: int,
    rng: random.Random,
):

    while True:

        p = random_prime(
            bits // 2,
            rng,
        )

        q = random_prime(
            bits - bits // 2,
            rng,
        )

        if p != q:
            return p, q, p * q


# ------------------------------------------------------------
# C identity
# ------------------------------------------------------------

def C_value(
    n: int,
    a: int,
    b: int,
    r: int,
    s: int,
) -> int:

    D = n - a * b

    beta = (
        D
        * pow(r, -1, s)
    ) % s

    alpha = (
        D
        * pow(s, -1, r)
    ) % r

    return (
        r * beta
        + s * alpha
        + a * b
    ) // (r * s)


# ------------------------------------------------------------
# Build observed C matrix
# ------------------------------------------------------------

def observed_C_grid(
    n: int,
    p: int,
    q: int,
):

    C = [
        [0] * len(R2)
        for _ in R1
    ]

    for i, r in enumerate(R1):

        a = p % r

        for j, s in enumerate(R2):

            b = q % s

            C[i][j] = C_value(
                n,
                a,
                b,
                r,
                s,
            )

    return C


# ------------------------------------------------------------
# Allowed b values for (r,s,C,a)
# ------------------------------------------------------------

def allowed_b(
    n: int,
    r: int,
    s: int,
    c: int,
    a: int,
):

    result = set()

    for b in range(1, s):

        if C_value(
            n,
            a,
            b,
            r,
            s,
        ) == c:

            result.add(b)

    return result


# ------------------------------------------------------------
# Precompute cell compatibility
# ------------------------------------------------------------

def build_support_tables(
    n: int,
    C,
):

    tables = {}

    for i, r in enumerate(R1):

        for j, s in enumerate(R2):

            c = C[i][j]

            by_a = {}

            for a in range(1, r):

                by_a[a] = allowed_b(
                    n,
                    r,
                    s,
                    c,
                    a,
                )

            tables[(i, j)] = by_a

    return tables


# ------------------------------------------------------------
# Candidate compatibility test
# ------------------------------------------------------------

def candidate_is_C_compatible(
    n: int,
    x: int,
    is_p_side: bool,
    C,
    support,
):

    """
    Test whether x can be one factor of n.

    We do NOT require x | n here.

    Instead, derive its residue vector and determine whether
    there EXISTS a compatible residue vector for the opposite
    factor.

    This is useful because the divisibility test can be delayed.
    """

    if x < 2:
        return False, None

    if is_p_side:

        A = [
            x % r
            for r in R1
        ]

        # For each column j, determine which b_j values
        # are compatible with every row i.
        B_domains = []

        for j, s in enumerate(R2):

            possible = set(
                range(1, s)
            )

            for i, _r in enumerate(R1):

                possible &= support[(i, j)][A[i]]

                if not possible:
                    break

            if not possible:
                return False, None

            B_domains.append(possible)

        return True, (
            tuple(A),
            tuple(
                frozenset(x)
                for x in B_domains
            ),
        )

    else:

        B = [
            x % s
            for s in R2
        ]

        A_domains = []

        for i, r in enumerate(R1):

            possible = set(
                range(1, r)
            )

            for j, _s in enumerate(R2):

                # Reverse lookup: find a values that support
                # this fixed b value.
                good = {
                    a
                    for a in range(1, r)
                    if B[j]
                    in support[(i, j)].get(
                        a,
                        set(),
                    )
                }

                possible &= good

                if not possible:
                    break

            if not possible:
                return False, None

            A_domains.append(
                possible
            )

        return True, (
            tuple(
                frozenset(x)
                for x in A_domains
            ),
            tuple(B),
        )


# ------------------------------------------------------------
# CRT helper
# ------------------------------------------------------------

def crt_list(
    residues,
    moduli,
):

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli,
    ):

        t = (
            (a - x)
            * pow(M, -1, m)
        ) % m

        x += M * t
        M *= m

    return x, M


# ------------------------------------------------------------
# Generate all x <= sqrt(n) from CRT coordinates
# ------------------------------------------------------------

def generate_crt_candidates(
    moduli,
    upper,
):

    """
    Since the CRT modulus exceeds upper, every x <= upper
    has a unique residue vector.

    Rather than looping over residue tuples, directly loop
    over x.

    This is deliberately a control measurement:
    how much of the 1..sqrt(n) interval survives the
    C-compatibility test?
    """

    for x in range(
        2,
        upper + 1,
    ):

        yield x


# ------------------------------------------------------------
# Factor-side sieve
# ------------------------------------------------------------

def sieve_factor_side(
    n: int,
    upper: int,
    C,
    support,
    is_p_side: bool,
):

    compatible = []
    exact = []

    tested = 0

    t0 = time.perf_counter()

    for x in generate_crt_candidates(
        None,
        upper,
    ):

        tested += 1

        ok, state = candidate_is_C_compatible(
            n,
            x,
            is_p_side,
            C,
            support,
        )

        if not ok:
            continue

        compatible.append(
            (x, state)
        )

        # Now perform the expensive exact divisor test.
        if n % x == 0:

            y = n // x

            exact.append(
                (x, y)
            )

    elapsed = (
        time.perf_counter()
        - t0
    )

    return {
        "tested": tested,
        "compatible": compatible,
        "exact": exact,
        "time": elapsed,
    }


# ------------------------------------------------------------
# Improved cached sieve
# ------------------------------------------------------------

def cached_factor_sieve(
    n: int,
    upper: int,
    C,
    support,
    is_p_side: bool,
):

    """
    Important optimization.

    The C-grid only depends on x through the residue vector
    on the chosen side.

    Cache compatibility by residue vector.

    For x <= sqrt(n), many values have repeated low-radix
    residue patterns when only a subset of the grid is used.
    """

    cache = {}

    compatible = []
    exact = []

    tested = 0
    cache_hits = 0

    t0 = time.perf_counter()

    for x in range(
        2,
        upper + 1,
    ):

        tested += 1

        if is_p_side:

            key = tuple(
                x % r
                for r in R1
            )

        else:

            key = tuple(
                x % s
                for s in R2
            )

        if key in cache:

            ok, state = cache[key]

            cache_hits += 1

        else:

            ok, state = candidate_is_C_compatible(
                n,
                x,
                is_p_side,
                C,
                support,
            )

            cache[key] = (
                ok,
                state,
            )

        if not ok:
            continue

        compatible.append(
            x
        )

        if n % x == 0:

            y = n // x

            exact.append(
                (x, y)
            )

    elapsed = (
        time.perf_counter()
        - t0
    )

    return {
        "tested": tested,
        "compatible": compatible,
        "exact": exact,
        "cache_hits": cache_hits,
        "cache_size": len(cache),
        "time": elapsed,
    }


# ------------------------------------------------------------
# Main experiment
# ------------------------------------------------------------

def run_sample(
    bits: int,
    sample: int,
    p: int,
    q: int,
    n: int,
):

    upper = math.isqrt(n)

    print()
    print(
        f"p={p}"
    )
    print(
        f"q={q}"
    )
    print(
        f"n={n}"
    )
    print(
        f"sqrt(n)={upper}"
    )

    t0 = time.perf_counter()

    C = observed_C_grid(
        n,
        p,
        q,
    )

    support = build_support_tables(
        n,
        C,
    )

    prep_time = (
        time.perf_counter()
        - t0
    )

    print()
    print("Observed C matrix:")

    for row in C:
        print(
            "  "
            + " ".join(
                str(x)
                for x in row
            )
        )

    print(
        f"support preparation="
        f"{prep_time:.6f}s"
    )

    # --------------------------------------------------------
    # Fast C-only sieve on both directions
    # --------------------------------------------------------

    print()
    print("C-only factor sieve:")

    pside = cached_factor_sieve(
        n,
        upper,
        C,
        support,
        True,
    )

    print()
    print("P SIDE")

    print(
        f"  tested={pside['tested']:,}"
    )

    print(
        f"  C-compatible="
        f"{len(pside['compatible']):,}"
    )

    print(
        f"  cache size="
        f"{pside['cache_size']:,}"
    )

    print(
        f"  cache hits="
        f"{pside['cache_hits']:,}"
    )

    print(
        f"  exact="
        f"{pside['exact']}"
    )

    print(
        f"  time="
        f"{pside['time']:.6f}s"
    )

    qside = cached_factor_sieve(
        n,
        upper,
        C,
        support,
        False,
    )

    print()
    print("Q SIDE")

    print(
        f"  tested={qside['tested']:,}"
    )

    print(
        f"  C-compatible="
        f"{len(qside['compatible']):,}"
    )

    print(
        f"  cache size="
        f"{qside['cache_size']:,}"
    )

    print(
        f"  cache hits="
        f"{qside['cache_hits']:,}"
    )

    print(
        f"  exact="
        f"{qside['exact']}"
    )

    print(
        f"  time="
        f"{qside['time']:.6f}s"
    )

    # --------------------------------------------------------
    # Compare to truth
    # --------------------------------------------------------

    p_survives = p in pside["compatible"]
    q_survives = q in qside["compatible"]

    found = (
        (p, q) in pside["exact"]
        or (q, p) in qside["exact"]
    )

    print()
    print("TRUTH CHECK")

    print(
        f"  p C-compatible="
        f"{p_survives}"
    )

    print(
        f"  q C-compatible="
        f"{q_survives}"
    )

    print(
        f"  p divides n="
        f"{n % p == 0}"
    )

    print(
        f"  q divides n="
        f"{n % q == 0}"
    )

    print(
        f"  exact recovery="
        f"{found}"
    )

    return {
        "bits": bits,
        "p_compatible": p_survives,
        "q_compatible": q_survives,
        "found": found,
        "p_compatible_count": len(
            pside["compatible"]
        ),
        "q_compatible_count": len(
            qside["compatible"]
        ),
        "p_time": pside["time"],
        "q_time": qside["time"],
    }


# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------

def main():

    rng = random.Random(SEED)

    print("START EXPERIMENT 167")
    print("=" * 72)

    print(
        "Global Factor Candidate -> "
        "C-Grid Compatibility"
    )

    print()
    print(
        "Experiment 166 searched from the anchor outward:"
    )

    print(
        "    anchor -> domains -> CRT candidates"
    )

    print()
    print(
        "This experiment reverses the direction:"
    )

    print(
        "    candidate factor -> residue vector -> C compatibility"
    )

    print()
    print(
        "For every x <= sqrt(n), determine whether its complete"
    )

    print(
        "residue vector admits ANY opposite residue vector that"
    )

    print(
        "satisfies the entire 5x5 C matrix."
    )

    print()
    print(
        "Only after a candidate survives all C constraints do we"
    )

    print(
        "perform the exact test n % x == 0."
    )

    print("=" * 72)

    summaries = []

    for bits in BIT_SIZES:

        print()
        print("-" * 72)
        print(
            f"BIT SIZE = {bits}"
        )
        print("-" * 72)

        for sample in range(
            1,
            SAMPLES_PER_SIZE + 1,
        ):

            p, q, n = make_semiprime(
                bits,
                rng,
            )

            result = run_sample(
                bits,
                sample,
                p,
                q,
                n,
            )

            summaries.append(
                result
            )

    print()
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)

    for result in summaries:

        print(
            f"{result['bits']:2d}-bit: "
            f"P-compatible="
            f"{result['p_compatible_count']:,} "
            f"Q-compatible="
            f"{result['q_compatible_count']:,} "
            f"found="
            f"{result['found']} "
            f"Ptime="
            f"{result['p_time']:.4f}s "
            f"Qtime="
            f"{result['q_time']:.4f}s"
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 167")
    print("=" * 72)


if __name__ == "__main__":
    main()
