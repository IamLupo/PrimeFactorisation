#!/usr/bin/env python3

import math
import random
import time


# =============================================================================
# KAPPA EXPERIMENT 72
# RANDOM PRODUCT / CRT RESIDUE EXPERIMENT
#
# Demonstrates:
#
#       n = p*q
#       M = p1*p2
#       r = n mod M
#       n = r + k*M
#
# Also verifies:
#
#       r mod p1 = n mod p1
#       r mod p2 = n mod p2
#
# and reconstructs r from those two residues using CRT.
#
# No CSV output
# No external dependencies
# =============================================================================


SEED = 20260828

TRIALS = 20

# Factor range for n = p*q
FACTOR_LOW = 100_000
FACTOR_HIGH = 500_000

# Range for p1 and p2
MOD_LOW = 100
MOD_HIGH = 10_000


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def gcd(a, b):
    return math.gcd(a, b)


def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)


def crt_two(a1, m1, a2, m2):
    """
    Solve:

        x = a1 (mod m1)
        x = a2 (mod m2)

    Requires gcd(m1, m2) = 1.

    Returns the least non-negative solution modulo m1*m2.
    """

    g = math.gcd(m1, m2)

    if g != 1:
        raise ValueError(
            f"CRT requires coprime moduli: gcd({m1}, {m2})={g}"
        )

    # x = a1 + m1*t
    #
    # a1 + m1*t = a2 (mod m2)
    #
    # m1*t = a2-a1 (mod m2)
    #
    inv = pow(m1, -1, m2)

    t = ((a2 - a1) * inv) % m2

    x = a1 + m1 * t

    return x % (m1 * m2)


# =============================================================================
# PRIME GENERATION
# =============================================================================

def is_prime(n):
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False

        d += 2

    return True


def generate_primes(low, high):
    return [
        x
        for x in range(low, high)
        if is_prime(x)
    ]


# =============================================================================
# RANDOM COPRIME MODULI
# =============================================================================

def random_coprime_pair(rng, low, high):
    """
    Pick p1,p2 with:

        p1 != p2
        gcd(p1,p2) = 1
    """

    while True:

        p1 = rng.randint(low, high)
        p2 = rng.randint(low, high)

        if p1 == p2:
            continue

        if math.gcd(p1, p2) != 1:
            continue

        return p1, p2


# =============================================================================
# SINGLE TRIAL
# =============================================================================

def run_trial(trial_id, p, q, p1, p2):
    n = p * q
    M = p1 * p2

    # ---------------------------------------------------------
    # Direct residue
    # ---------------------------------------------------------

    r = n % M

    # ---------------------------------------------------------
    # Recover k
    # ---------------------------------------------------------

    k = (n - r) // M

    # Verify exact identity
    identity_ok = (
        r + k * M == n
    )

    # ---------------------------------------------------------
    # Local residues
    # ---------------------------------------------------------

    a1 = n % p1
    a2 = n % p2

    # ---------------------------------------------------------
    # CRT reconstruction
    # ---------------------------------------------------------

    r_crt = crt_two(
        a1,
        p1,
        a2,
        p2
    )

    crt_ok = (
        r_crt == r
        and r_crt % p1 == a1
        and r_crt % p2 == a2
    )

    # ---------------------------------------------------------
    # Product represented by r+kM
    # ---------------------------------------------------------

    recovered_n = r + k * M

    # ---------------------------------------------------------
    # Number of possible k values if p,q are restricted
    # to [FACTOR_LOW, FACTOR_HIGH)
    #
    # L^2 <= r+kM <= U^2
    # ---------------------------------------------------------

    L = FACTOR_LOW
    U = FACTOR_HIGH - 1

    lower_product = L * L
    upper_product = U * U

    k_min = max(
        0,
        math.ceil(
            (lower_product - r) / M
        )
    )

    k_max = math.floor(
        (upper_product - r) / M
    )

    if k_min <= k_max:
        k_count = k_max - k_min + 1
    else:
        k_count = 0

    # ---------------------------------------------------------
    # Relative scales
    # ---------------------------------------------------------

    M_over_n = M / n
    n_over_M = n / M

    return {
        "trial": trial_id,
        "p": p,
        "q": q,
        "n": n,
        "p1": p1,
        "p2": p2,
        "M": M,
        "r": r,
        "k": k,
        "a1": a1,
        "a2": a2,
        "r_crt": r_crt,
        "identity_ok": identity_ok,
        "crt_ok": crt_ok,
        "recovered_n": recovered_n,
        "k_min": k_min,
        "k_max": k_max,
        "k_count": k_count,
        "M_over_n": M_over_n,
        "n_over_M": n_over_M,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    rng = random.Random(SEED)

    start_total = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 72")
    print("RANDOM PRODUCT / CRT RESIDUE EXPERIMENT")
    print("=" * 78)
    print()

    print(f"seed              = {SEED}")
    print(f"trials            = {TRIALS}")
    print(
        f"factor range      = "
        f"[{FACTOR_LOW:,}, {FACTOR_HIGH:,})"
    )
    print(
        f"modulus range     = "
        f"[{MOD_LOW:,}, {MOD_HIGH:,}]"
    )
    print()

    # =========================================================================
    # 1. PRIME POPULATION
    # =========================================================================

    print("=" * 78)
    print("1. PRIME POPULATION")
    print("=" * 78)

    start = time.perf_counter()

    # For this demonstration we generate a smaller candidate population
    # containing enough primes for random sampling.
    primes = generate_primes(
        FACTOR_LOW,
        FACTOR_HIGH
    )

    generation_time = time.perf_counter() - start

    print(
        f"prime population = "
        f"{len(primes):,}"
    )

    print(
        f"generation time  = "
        f"{generation_time:.6f}s"
    )

    print()

    # =========================================================================
    # 2. RANDOM TRIALS
    # =========================================================================

    print("=" * 78)
    print("2. RANDOM n = p*q / M = p1*p2 TRIALS")
    print("=" * 78)
    print()

    results = []

    for trial in range(1, TRIALS + 1):

        # Random semiprime factors
        p, q = rng.sample(primes, 2)

        # Random coprime moduli
        p1, p2 = random_coprime_pair(
            rng,
            MOD_LOW,
            MOD_HIGH
        )

        result = run_trial(
            trial,
            p,
            q,
            p1,
            p2
        )

        results.append(result)

        print(
            f"TRIAL {trial:2d}"
        )

        print(
            f"  p,q        = "
            f"{result['p']:,} * {result['q']:,}"
        )

        print(
            f"  n          = "
            f"{result['n']:,}"
        )

        print(
            f"  p1,p2      = "
            f"{result['p1']:,} * {result['p2']:,}"
        )

        print(
            f"  M          = "
            f"{result['M']:,}"
        )

        print(
            f"  r = n mod M= "
            f"{result['r']:,}"
        )

        print(
            f"  k          = "
            f"{result['k']:,}"
        )

        print(
            f"  check      = "
            f"{result['r']:,} + "
            f"{result['k']:,} * "
            f"{result['M']:,} "
            f"== {result['n']:,}"
        )

        print(
            f"  identity_ok = "
            f"{result['identity_ok']}"
        )

        print()

        print(
            f"  local residue mod p1 = "
            f"{result['a1']:,}"
        )

        print(
            f"  local residue mod p2 = "
            f"{result['a2']:,}"
        )

        print(
            f"  CRT reconstructed r   = "
            f"{result['r_crt']:,}"
        )

        print(
            f"  CRT check             = "
            f"{result['crt_ok']}"
        )

        print()

        print(
            f"  n/M          = "
            f"{result['n_over_M']:.6f}"
        )

        print(
            f"  M/n          = "
            f"{result['M_over_n']:.12f}"
        )

        print(
            f"  k interval in "
            f"[L^2,U^2] = "
            f"[{result['k_min']:,}, "
            f"{result['k_max']:,}]"
        )

        print(
            f"  number of possible k = "
            f"{result['k_count']:,}"
        )

        print("-" * 78)

    # =========================================================================
    # 3. GLOBAL IDENTITY CHECK
    # =========================================================================

    print()
    print("=" * 78)
    print("3. GLOBAL CORRECTNESS")
    print("=" * 78)

    identity_failures = sum(
        not x["identity_ok"]
        for x in results
    )

    crt_failures = sum(
        not x["crt_ok"]
        for x in results
    )

    recovered_failures = sum(
        x["recovered_n"] != x["n"]
        for x in results
    )

    print(
        f"identity failures      = "
        f"{identity_failures}"
    )

    print(
        f"CRT failures           = "
        f"{crt_failures}"
    )

    print(
        f"n reconstruction fails = "
        f"{recovered_failures}"
    )

    print(
        "status = "
        + (
            "PASS"
            if identity_failures == 0
            and crt_failures == 0
            and recovered_failures == 0
            else "FAIL"
        )
    )

    # =========================================================================
    # 4. HOW MANY k VALUES?
    # =========================================================================

    print()
    print("=" * 78)
    print("4. k-SPACE ANALYSIS")
    print("=" * 78)

    k_counts = [
        x["k_count"]
        for x in results
    ]

    zero = sum(
        k == 0
        for k in k_counts
    )

    one = sum(
        k == 1
        for k in k_counts
    )

    multiple = sum(
        k > 1
        for k in k_counts
    )

    print(
        f"trials with no possible k = "
        f"{zero}"
    )

    print(
        f"trials with exactly one k = "
        f"{one}"
    )

    print(
        f"trials with multiple k     = "
        f"{multiple}"
    )

    print()

    print(
        "This measures how sparse the arithmetic progression"
    )

    print(
        "    r, r+M, r+2M, ..."
    )

    print(
        "is inside the possible product range."
    )

    # =========================================================================
    # 5. EXTREME M CASE
    # =========================================================================

    print()
    print("=" * 78)
    print("5. LARGE-M LIMIT")
    print("=" * 78)

    print("""
For

    n = r + kM

the spacing between possible values of n having the same residue r
is exactly M.

Therefore:

    larger M  -> fewer possible k values
    M > product-range width -> at most one candidate product

However, this does NOT mean r becomes easier to discover.

The residue remains:

    r = n mod M.

The experiment therefore separates two questions:

    A. Given n, how sparse is r+kM?
    B. Without knowing n's residue, can r itself be predicted?

Question A is purely arithmetic.
Question B is the potentially difficult research problem.
""")

    # =========================================================================
    # 6. SUMMARY TABLE
    # =========================================================================

    print()
    print("=" * 78)
    print("6. COMPACT SUMMARY")
    print("=" * 78)

    print(
        f"{'trial':>5s} "
        f"{'M':>12s} "
        f"{'r':>12s} "
        f"{'k':>12s} "
        f"{'k-count':>12s} "
        f"{'CRT':>6s}"
    )

    print("-" * 78)

    for x in results:

        print(
            f"{x['trial']:5d} "
            f"{x['M']:12,d} "
            f"{x['r']:12,d} "
            f"{x['k']:12,d} "
            f"{x['k_count']:12,d} "
            f"{str(x['crt_ok']):>6s}"
        )

    # =========================================================================
    # 7. FINAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("7. FINAL INTERPRETATION")
    print("=" * 78)

    print("""
The experiment confirms three facts.

1. Every integer n has a unique residue

       r = n mod M

   with

       n = r + kM.

2. For M=p1*p2 with gcd(p1,p2)=1, the same r is completely
   determined by the pair

       n mod p1
       n mod p2

   through CRT.

3. Increasing M makes the progression

       r + kM

   sparser, but does not by itself provide a method for discovering r.

Therefore the interesting factoring question is not whether the
representation

       n = r + kM

exists.

It always does.

The interesting question is whether some Kappa identity can constrain
or predict r modulo a large M using information that does not already
require testing candidate factors.

If such an identity existed, the resulting candidate set

       r + kM

could potentially become extremely small.
""")

    # =========================================================================
    # FINAL
    # =========================================================================

    total_time = (
        time.perf_counter() - start_total
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 72 COMPLETE")
    print("=" * 78)

    print(
        f"total runtime = "
        f"{total_time:.6f}s"
    )


if __name__ == "__main__":
    main()
