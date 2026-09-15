#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 125
# Multi-radix CRT residue compression
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 125")
print("Multi-radix CRT residue compression")
print("=" * 72)


# ------------------------------------------------------------------------
# RADICES
#
# All radices within each side are pairwise coprime.
#
# Product of R1:
#
#   17 * 43 * 47 = 34357
#
# Product of R2:
#
#   19 * 37 * 53 = 37271
#
# Combined modulus:
#
#   34357 * 37271 ~= 1.28e9
#
# This is intentionally much larger than the moduli used previously.
# ------------------------------------------------------------------------

R1 = [17, 43, 47]
R2 = [19, 37, 53]

M1 = math.prod(R1)
M2 = math.prod(R2)

print()
print("R1 =", R1)
print("R2 =", R2)
print()
print("M1 =", M1)
print("M2 =", M2)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
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


# ------------------------------------------------------------------------
# RANDOM PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= (1 << (bits - 1))
        x |= 1

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        n = p * q

        if n.bit_length() == bits:
            return n, min(p, q), max(p, q)


# ------------------------------------------------------------------------
# CRT
# ------------------------------------------------------------------------

def crt_pairwise(residues, moduli):

    """
    Solve

        x == residues[i] mod moduli[i]

    for pairwise-coprime moduli.
    """

    x = 0
    M = 1

    for a, m in zip(residues, moduli):

        inv = pow(M, -1, m)

        t = ((a - x) * inv) % m

        x += M * t
        M *= m

        x %= M

    return x, M


# ------------------------------------------------------------------------
# DIRECT SQRT SCAN
# ------------------------------------------------------------------------

def direct_factor(n):

    limit = math.isqrt(n)

    tested = 0

    for p in range(2, limit + 1):

        tested += 1

        if n % p == 0:

            return p, n // p, tested

    return None, None, tested


# ------------------------------------------------------------------------
# CRT-RESIDUE SEARCH
#
# We enumerate the possible residue vector:
#
#   a_i = p mod r1_i
#
# giving one combined A modulo M1.
#
# Then:
#
#   q = n * A^(-1) mod M1.
#
# This gives a simultaneous residue condition on q.
#
# We then search only values of p inside the sqrt(n) range that satisfy
#
#   p == A mod M1.
#
# The important statistic is how many p values remain.
# ------------------------------------------------------------------------

def crt_residue_search(n, true_p, true_q):

    sqrt_n = math.isqrt(n)

    tested_residue_vectors = 0
    invertible_vectors = 0
    valid_q_residue_vectors = 0

    total_p_candidates = 0

    exact_found = []

    # ------------------------------------------------------------
    # Enumerate all possible p residue vectors.
    # ------------------------------------------------------------

    for a0 in range(R1[0]):

        for a1 in range(R1[1]):

            for a2 in range(R1[2]):

                residues = [a0, a1, a2]

                tested_residue_vectors += 1

                A, modulus = crt_pairwise(
                    residues,
                    R1
                )

                # ------------------------------------------------
                # p must be coprime to modulus in order to invert.
                # ------------------------------------------------

                if math.gcd(A, modulus) != 1:
                    continue

                invertible_vectors += 1

                qres = (
                    n * pow(A, -1, modulus)
                ) % modulus

                # ------------------------------------------------
                # q must have some valid residue modulo every
                # R2[j].
                #
                # We don't know b_j yet, so any residue is allowed.
                # The test is therefore primarily the fact that q
                # shares the same CRT-derived modulus.
                # ------------------------------------------------

                valid_q_residue_vectors += 1

                # ------------------------------------------------
                # Find the first positive p in this residue class.
                #
                # p = A + t*M1
                # ------------------------------------------------

                if A == 0:
                    first = modulus
                else:
                    first = A

                if first < 2:
                    first += (
                        ((2 - first + modulus - 1) // modulus)
                        * modulus
                    )

                if first > sqrt_n:
                    continue

                # Number of p values in this CRT class below sqrt(n).
                count = (
                    (sqrt_n - first) // modulus
                ) + 1

                total_p_candidates += count

                # ------------------------------------------------
                # Test candidates from this residue class.
                # ------------------------------------------------

                p = first

                for _ in range(count):

                    if n % p == 0:

                        q = n // p

                        if (
                            {p, q}
                            == {true_p, true_q}
                        ):
                            exact_found.append(
                                (
                                    p,
                                    q,
                                    A,
                                    modulus,
                                    qres,
                                )
                            )

                    p += modulus

    return {
        "tested_residue_vectors":
            tested_residue_vectors,

        "invertible_vectors":
            invertible_vectors,

        "valid_q_residue_vectors":
            valid_q_residue_vectors,

        "total_p_candidates":
            total_p_candidates,

        "exact_found":
            exact_found,
    }


# ------------------------------------------------------------------------
# TRUE RESIDUE DISPLAY
# ------------------------------------------------------------------------

def print_true_residues(p, q):

    a = [p % r for r in R1]
    b = [q % r for r in R2]

    A, M = crt_pairwise(a, R1)
    B, N = crt_pairwise(b, R2)

    print()
    print("TRUE p residues")
    print("a =", a)

    print()
    print("TRUE q residues")
    print("b =", b)

    print()
    print("TRUE CRT")
    print("p mod M1 =", A)
    print("M1       =", M)

    print("q mod M2 =", B)
    print("M2       =", N)


# ------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------

random.seed(125)

TEST_BITS = [
    30,
    36,
    42,
    48,
]

for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, true_p, true_q = random_semiprime(bits)

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())
    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print_true_residues(
        true_p,
        true_q
    )

    # ------------------------------------------------------------
    # Direct control
    # ------------------------------------------------------------

    print()
    print("DIRECT SQRT(n) SCAN")

    t0 = time.perf_counter()

    p0, q0, tested0 = direct_factor(n)

    t1 = time.perf_counter()

    direct_time = t1 - t0

    direct_ok = (
        {p0, q0}
        == {true_p, true_q}
    )

    print()
    print("p tested =", tested0)
    print("runtime  =", f"{direct_time:.6f}", "s")
    print("correct  =", direct_ok)

    # ------------------------------------------------------------
    # CRT search
    # ------------------------------------------------------------

    print()
    print("CRT RESIDUE SEARCH")

    t0 = time.perf_counter()

    result = crt_residue_search(
        n,
        true_p,
        true_q,
    )

    t1 = time.perf_counter()

    crt_time = t1 - t0

    print()
    print(
        "residue vectors tested =",
        result["tested_residue_vectors"]
    )

    print(
        "invertible vectors    =",
        result["invertible_vectors"]
    )

    print(
        "valid q vectors       =",
        result["valid_q_residue_vectors"]
    )

    print(
        "remaining p candidates =",
        result["total_p_candidates"]
    )

    print(
        "exact candidates found =",
        len(result["exact_found"])
    )

    print(
        "runtime =",
        f"{crt_time:.6f}",
        "s"
    )

    recovered = any(
        {x[0], x[1]}
        == {true_p, true_q}
        for x in result["exact_found"]
    )

    print(
        "TRUE FACTORIZATION FOUND =",
        recovered
    )

    # ------------------------------------------------------------
    # Reduction factor
    # ------------------------------------------------------------

    if result["total_p_candidates"] > 0:

        reduction = (
            tested0
            / result["total_p_candidates"]
        )

        print()
        print(
            "sqrt-search candidate reduction =",
            f"{reduction:.3f}x"
        )

    # ------------------------------------------------------------
    # Timing comparison
    # ------------------------------------------------------------

    print()
    print("TIMING")

    if crt_time > 0:

        print(
            "direct / CRT =",
            f"{direct_time / crt_time:.3f}x"
        )

    # ------------------------------------------------------------
    # Exact result
    # ------------------------------------------------------------

    if result["exact_found"]:

        print()
        print("EXACT RECOVERY")

        for item in result["exact_found"]:

            p, q, A, M, qres = item

            print("p =", p)
            print("q =", q)
            print("A =", A)
            print("M =", M)
            print("q residue from n*A^-1 =", qres)


# ========================================================================
# FINISHED EXPERIMENT 125
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 125")
print("=" * 72)
