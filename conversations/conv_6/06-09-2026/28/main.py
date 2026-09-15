#!/usr/bin/env python3

"""
START EXPERIMENT 184

MODULAR CO-FACTOR SIEVE

Experiments 180-183 established that:

    * C-values derived from n are tautological.
    * p-residue classes alone do not select the factor.
    * q-residue classes are determined by p through pq=n mod M.
    * CRT therefore partitions the ordinary factor search rather than
      eliminating its fundamental ambiguity.

Experiment 184 changes direction.

We keep the useful modular identity

    q ≡ n*p^{-1} (mod t)

but use it as a SIEVE.

For every small prime t not dividing n, a valid factor p must satisfy:

    p != 0 mod t

and

    q = n/p != 0 mod t.

Therefore:

    p mod t != 0
    n * inverse(p mod t, t) mod t != 0.

This removes candidates for which either factor is divisible by t.

This is essentially a modular wheel for the factorization problem.

The experiment compares:

    A) plain odd p scanning
    B) a wheel based on divisibility of p
    C) a stronger co-factor-aware sieve that rejects candidates when q
       is divisible by a small prime
    D) the same sieve restricted to a CRT progression.

The purpose is NOT to claim a new asymptotic factorization algorithm.

Instead we measure:

    * candidate reduction
    * actual n % p tests
    * wall-clock time
    * scaling with bit size.

The important question is whether the q-side modular information provides
any substantial gain beyond ordinary wheel factorization.

No C-values are used.

No factor residues are assumed.

The hidden factors are used only for verification.

"""


import math
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

# Small primes used by the modular sieve.
SIEVE_PRIMES = [
    3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43,
    47, 53, 59, 61, 67, 71,
    73, 79, 83, 89,
]

# Baseline hard cap so a bad experiment does not run forever.
MAX_TESTS = 20_000_000

# Only test these CRT moduli.
CRT_MODULI = [
    [17],
    [17, 43],
    [17, 43, 59],
]


# ============================================================
# PRIME TEST
# ============================================================

def is_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37,
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = 41
    step = 2

    while d * d <= n:

        if n % d == 0:
            return False

        d += step
        step = 6 - step

    return True


# ============================================================
# EGCD / INVERSE
# ============================================================

def egcd(a, b):

    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(
        b,
        a % b,
    )

    return (
        g,
        y1,
        x1 - (a // b) * y1,
    )


def inv_mod(a, m):

    g, x, _ = egcd(
        a,
        m,
    )

    if g != 1:

        raise ValueError(
            f"{a} not invertible mod {m}"
        )

    return x % m


# ============================================================
# CRT
# ============================================================

def crt_pair(
    a1,
    m1,
    a2,
    m2,
):

    g = math.gcd(
        m1,
        m2,
    )

    if (a2 - a1) % g != 0:
        return None

    m1r = m1 // g
    m2r = m2 // g

    rhs = (
        (a2 - a1)
        // g
    ) % m2r

    if m2r == 1:

        t = 0

    else:

        t = (
            rhs
            * inv_mod(
                m1r % m2r,
                m2r,
            )
        ) % m2r

    x = a1 + m1 * t

    modulus = m1 * m2r

    return (
        x % modulus,
        modulus,
    )


# ============================================================
# SEMIPRIME
# ============================================================

def make_semiprime(bits):

    low = 1 << (
        bits // 2 - 1
    )

    high = 1 << (
        bits // 2 + 1
    )

    for p in range(
        low | 1,
        high,
        2,
    ):

        if not is_prime(p):
            continue

        target = (
            (1 << bits)
            // p
        )

        for delta in range(
            -1000,
            1001,
            2,
        ):

            q = target + delta

            if q <= p:
                continue

            if not is_prime(q):
                continue

            n = p * q

            if n.bit_length() == bits:

                return p, q, n

    raise RuntimeError(
        f"cannot generate {bits}-bit semiprime"
    )


# ============================================================
# BASELINE ODD SCAN
# ============================================================

def baseline_scan(n):

    limit = math.isqrt(n)

    tests = 0

    start = time.perf_counter()

    # p must be odd for the test semiprimes.
    for p in range(
        3,
        limit + 1,
        2,
    ):

        tests += 1

        if tests > MAX_TESTS:

            return {
                "found": None,
                "tests": tests,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_TESTS",
            }

        if n % p == 0:

            return {
                "found": p,
                "tests": tests,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "tests": tests,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# P-ONLY WHEEL
# ============================================================

def p_only_sieve_scan(
    n,
    primes,
):

    limit = math.isqrt(n)

    tests = 0

    candidates = 0

    rejected = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2,
    ):

        bad = False

        for t in primes:

            if p == t:
                continue

            if p % t == 0:

                bad = True
                rejected += 1
                break

        if bad:
            continue

        candidates += 1

        tests += 1

        if tests > MAX_TESTS:

            return {
                "found": None,
                "tests": tests,
                "candidates":
                    candidates,
                "rejected": rejected,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_TESTS",
            }

        if n % p == 0:

            return {
                "found": p,
                "tests": tests,
                "candidates":
                    candidates,
                "rejected": rejected,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "tests": tests,
        "candidates":
            candidates,
        "rejected": rejected,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# CO-FACTOR SIEVE
# ============================================================

def cofactor_sieve_scan(
    n,
    primes,
):

    limit = math.isqrt(n)

    tests = 0

    candidates = 0

    p_rejected = 0

    q_rejected = 0

    start = time.perf_counter()

    # --------------------------------------------------------
    # For each small t, precompute whether it divides n.
    # If t | n we cannot use inverse modulo t.
    # --------------------------------------------------------

    usable_primes = [
        t
        for t in primes
        if n % t != 0
    ]

    for p in range(
        3,
        limit + 1,
        2,
    ):

        bad = False

        for t in usable_primes:

            # p itself divisible by t.
            if p % t == 0:

                # p=t is a legitimate possible factor,
                # so don't reject that exact small prime.
                if p != t:

                    bad = True
                    p_rejected += 1
                    break

                continue

            # q = n/p mod t.
            #
            # If this residue is 0, then t divides q.
            #
            # n mod t is nonzero and p is invertible.
            inv_p = inv_mod(
                p % t,
                t,
            )

            q_mod_t = (
                (n % t)
                * inv_p
            ) % t

            if q_mod_t == 0:

                bad = True
                q_rejected += 1
                break

        if bad:
            continue

        candidates += 1

        tests += 1

        if tests > MAX_TESTS:

            return {
                "found": None,
                "tests": tests,
                "candidates":
                    candidates,
                "p_rejected":
                    p_rejected,
                "q_rejected":
                    q_rejected,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_TESTS",
            }

        if n % p == 0:

            return {
                "found": p,
                "tests": tests,
                "candidates":
                    candidates,
                "p_rejected":
                    p_rejected,
                "q_rejected":
                    q_rejected,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "tests": tests,
        "candidates":
            candidates,
        "p_rejected":
            p_rejected,
        "q_rejected":
            q_rejected,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# PRECOMPUTED MODULAR REJECTION TABLE
# ============================================================

def build_rejection_table(
    n,
    primes,
):

    table = {}

    for t in primes:

        if n % t == 0:
            continue

        allowed = set()

        n_mod = n % t

        for x in range(
            1,
            t,
        ):

            inv = inv_mod(
                x,
                t,
            )

            q_mod = (
                n_mod
                * inv
            ) % t

            # x=p mod t is valid iff neither p nor q
            # is zero modulo t.

            if q_mod != 0:

                allowed.add(x)

        table[t] = allowed

    return table


# ============================================================
# CRT-PROGRESSION + COFACTOR SIEVE
# ============================================================

def crt_progression_scan(
    n,
    radices,
    sieve_primes,
):

    limit = math.isqrt(n)

    # --------------------------------------------------------
    # Build the full CRT modulus.
    #
    # Instead of enumerating residue vectors, consider EVERY
    # residue P modulo M that is admissible.
    #
    # This is intentionally a pure factor-blind search.
    # --------------------------------------------------------

    M = math.prod(radices)

    tests = 0

    candidates = 0

    rejected = 0

    start = time.perf_counter()

    # We can enumerate p directly in [1,sqrt(n)] and only
    # retain values that are units / satisfy the modular sieve.
    #
    # This isolates the effect of the CRT modulus from the
    # residue-vector generation problem.

    usable = [
        t
        for t in sieve_primes
        if n % t != 0
    ]

    # Cache allowed residues for each sieve prime.
    table = build_rejection_table(
        n,
        usable,
    )

    for p in range(
        3,
        limit + 1,
        2,
    ):

        # p modulo the CRT modulus.
        #
        # Every p automatically belongs to one CRT residue
        # class. This makes the point explicit: the CRT modulus
        # itself does not remove p values.
        p_crt = p % M

        if p_crt == 0:
            # p is divisible by every radix, impossible for
            # these balanced semiprimes.
            rejected += 1
            continue

        bad = False

        for t in usable:

            residue = p % t

            # Exact p residue check.
            if residue == 0:

                if p != t:

                    bad = True
                    rejected += 1
                    break

                continue

            if residue not in table[t]:

                bad = True
                rejected += 1
                break

        if bad:
            continue

        candidates += 1

        tests += 1

        if tests > MAX_TESTS:

            return {
                "found": None,
                "M": M,
                "tests": tests,
                "candidates":
                    candidates,
                "rejected":
                    rejected,
                "time":
                    time.perf_counter()
                    - start,
                "aborted":
                    "MAX_TESTS",
            }

        if n % p == 0:

            return {
                "found": p,
                "M": M,
                "tests": tests,
                "candidates":
                    candidates,
                "rejected":
                    rejected,
                "time":
                    time.perf_counter()
                    - start,
                "aborted": None,
            }

    return {
        "found": None,
        "M": M,
        "tests": tests,
        "candidates":
            candidates,
        "rejected":
            rejected,
        "time":
            time.perf_counter()
            - start,
        "aborted": None,
    }


# ============================================================
# RESIDUE DENSITY
# ============================================================

def theoretical_density(
    n,
    primes,
):

    density = 1.0

    for t in primes:

        if n % t == 0:
            continue

        # Among residues modulo t:
        #
        # p != 0
        # q != 0
        #
        # both conditions remove two residue classes.
        #
        # Since q = n/p and n is a unit:
        #
        # p != 0 and q != 0 are actually the same condition:
        #
        # q = 0 iff p = 0.
        #
        # Therefore only ONE class is removed.
        #
        density *= (
            (t - 1) / t
        )

    return density


# ============================================================
# RUN INSTANCE
# ============================================================

def run_instance(bits):

    print()
    print("=" * 72)
    print(
        f"START INSTANCE {bits}-BIT"
    )
    print("=" * 72)

    p, q, n = make_semiprime(
        bits
    )

    limit = math.isqrt(n)

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = {limit}"
    )

    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    # --------------------------------------------------------
    # Baseline.
    # --------------------------------------------------------

    print()
    print(
        "BASELINE ODD SCAN"
    )

    baseline = baseline_scan(
        n
    )

    print(
        f"    tests = "
        f"{baseline['tests']}"
    )

    print(
        f"    time = "
        f"{baseline['time']:.6f}s"
    )

    if baseline.get("found"):

        print(
            f"    found = "
            f"{baseline['found']}"
        )

    if baseline.get("aborted"):

        print(
            f"    status = "
            f"{baseline['aborted']}"
        )

    # --------------------------------------------------------
    # P-only wheel.
    # --------------------------------------------------------

    print()
    print(
        "P-ONLY MODULAR SIEVE"
    )

    p_only = p_only_sieve_scan(
        n,
        SIEVE_PRIMES,
    )

    print(
        f"    candidates = "
        f"{p_only['candidates']}"
    )

    print(
        f"    tests = "
        f"{p_only['tests']}"
    )

    print(
        f"    time = "
        f"{p_only['time']:.6f}s"
    )

    if p_only.get("found"):

        print(
            f"    found = "
            f"{p_only['found']}"
        )

    if p_only.get("aborted"):

        print(
            f"    status = "
            f"{p_only['aborted']}"
        )

    # --------------------------------------------------------
    # Co-factor-aware sieve.
    # --------------------------------------------------------

    print()
    print(
        "CO-FACTOR-AWARE SIEVE"
    )

    cof = cofactor_sieve_scan(
        n,
        SIEVE_PRIMES,
    )

    print(
        f"    candidates = "
        f"{cof['candidates']}"
    )

    print(
        f"    tests = "
        f"{cof['tests']}"
    )

    print(
        f"    p rejected = "
        f"{cof['p_rejected']}"
    )

    print(
        f"    q rejected = "
        f"{cof['q_rejected']}"
    )

    print(
        f"    time = "
        f"{cof['time']:.6f}s"
    )

    if cof.get("found"):

        print(
            f"    found = "
            f"{cof['found']}"
        )

    if cof.get("aborted"):

        print(
            f"    status = "
            f"{cof['aborted']}"
        )

    # --------------------------------------------------------
    # Theoretical density.
    # --------------------------------------------------------

    density = theoretical_density(
        n,
        SIEVE_PRIMES,
    )

    print()
    print(
        "THEORETICAL SIEVE DENSITY"
    )

    print(
        f"    density = "
        f"{density:.12e}"
    )

    print(
        f"    expected candidates <= sqrt(n) "
        f"≈ {int(limit * density)}"
    )

    # --------------------------------------------------------
    # CRT progression experiments.
    # --------------------------------------------------------

    print()
    print(
        "CRT + CO-FACTOR SIEVE"
    )

    for radices in CRT_MODULI:

        print()
        print(
            f"    radices = "
            f"{radices}"
        )

        result = crt_progression_scan(
            n,
            radices,
            SIEVE_PRIMES,
        )

        print(
            f"        M = "
            f"{result['M']}"
        )

        print(
            f"        candidates = "
            f"{result['candidates']}"
        )

        print(
            f"        tests = "
            f"{result['tests']}"
        )

        print(
            f"        rejected = "
            f"{result['rejected']}"
        )

        print(
            f"        time = "
            f"{result['time']:.6f}s"
        )

        if result.get("found"):

            print(
                f"        FOUND = "
                f"{result['found']}"
            )

        if result.get("aborted"):

            print(
                f"        status = "
                f"{result['aborted']}"
            )

    # --------------------------------------------------------
    # Verification.
    # --------------------------------------------------------

    print()
    print(
        "VERIFICATION"
    )

    print(
        f"    p is prime = "
        f"{is_prime(p)}"
    )

    print(
        f"    q is prime = "
        f"{is_prime(q)}"
    )

    for t in SIEVE_PRIMES:

        if n % t == 0:
            continue

        p_ok = (
            p % t != 0
        )

        q_ok = (
            q % t != 0
        )

        if not (
            p_ok and q_ok
        ):

            print(
                f"    WARNING t={t}: "
                f"factor violates sieve"
            )

    print()
    print(
        f"FINISHED INSTANCE {bits}-BIT"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print(
        "START EXPERIMENT 184"
    )

    print()

    print(
        f"R1 = {R1}"
    )

    print(
        f"R2 = {R2}"
    )

    print(
        f"BITS = {BITS}"
    )

    print(
        f"SIEVE_PRIMES = {SIEVE_PRIMES}"
    )

    print(
        f"MAX_TESTS = {MAX_TESTS}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 184"
    )


if __name__ == "__main__":
    main()
