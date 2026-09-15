#!/usr/bin/env python3

"""
==============================================================================
KAPPA NEXT EXPERIMENT
ADVERSARIAL N-ONLY / S-RECOVERY COLLISION SEARCH
==============================================================================

Goal
----

The previous experiment showed that, among only 12 random targets, a large
multi-modulus n-signature sometimes uniquely identified the observed s-signature.

That is not yet evidence of an n -> s relation, because the sample is tiny.

This experiment attacks the proposed relation directly.

For each auxiliary modulus

    m = F(r) = r^2-r+1

we search for residue pairs

    (x,y)
    (u,v)

such that

    x*y == u*v (mod m)

but

    x+y != u+v (mod m).

If such pairs exist, then n mod m does NOT determine s mod m.

We then try to lift the residue pairs to actual prime factors:

    p == x (mod m)
    q == y (mod m)

and

    P == u (mod m)
    Q == v (mod m).

This produces actual semiprimes n1=pq and n2=PQ satisfying

    n1 == n2 (mod m)

while

    p+q != P+Q (mod m).

That is an adversarial counterexample to any n-only S-recovery rule modulo m.

No CSV output.
Only stdout.
"""

import math
import random
from collections import defaultdict


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47
]

# Number of residue collisions to inspect per modulus.
MAX_COLLISIONS_PER_MODULUS = 8

# Number of actual prime lifts attempted for each residue pair.
PRIME_LIFT_ATTEMPTS = 3000

# Prime size for lifted factors.
# Keep factors in roughly this range so the resulting n is nontrivial.
LIFT_BITS = 24

# For combined-modulus experiments we use selected moduli.
COMBINED_R_SETS = [
    [2, 3],
    [3, 7],
    [7, 11],
    [13, 17],
    [19, 23],
    [29, 31],
]

random.seed(SEED)


# ============================================================================
# BASIC POLYNOMIAL
# ============================================================================

def F(x):
    return x * x - x + 1


# ============================================================================
# DETERMINISTIC MILLER-RABIN FOR 64-BIT INTEGERS
# ============================================================================

def is_prime(n):
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

    # Deterministic for all 64-bit integers.
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

        witness = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                witness = False
                break

        if witness:
            return False

    return True


# ============================================================================
# PRIME LIFT
# ============================================================================

def random_prime_with_residue(residue, modulus, bits=LIFT_BITS):
    """
    Find a prime p satisfying

        p == residue (mod modulus)

    using random lifts.
    """

    residue %= modulus

    low = 1 << (bits - 1)
    high = (1 << bits) - 1

    first_k = max(0, (low - residue + modulus - 1) // modulus)
    last_k = (high - residue) // modulus

    if first_k > last_k:
        return None

    for _ in range(PRIME_LIFT_ATTEMPTS):
        k = random.randint(first_k, last_k)
        candidate = residue + k * modulus

        if candidate >= 2 and is_prime(candidate):
            return candidate

    return None


# ============================================================================
# RESIDUE COLLISION SEARCH
# ============================================================================

def find_product_collisions(m):
    """
    Build

        product -> [(x,y,sum), ...]

    over units modulo m.

    We deliberately use units because actual prime factors are normally
    invertible modulo the auxiliary modulus unless they divide m.
    """

    buckets = defaultdict(list)

    for x in range(1, m):
        if math.gcd(x, m) != 1:
            continue

        inv_x = pow(x, -1, m)

        for y in range(1, m):
            if math.gcd(y, m) != 1:
                continue

            product = (x * y) % m
            summation = (x + y) % m

            buckets[product].append((x, y, summation))

    collisions = []

    for product, entries in buckets.items():

        by_sum = {}

        for x, y, summation in entries:
            by_sum.setdefault(summation, (x, y))

        if len(by_sum) <= 1:
            continue

        sums = list(by_sum.keys())

        # Take distinct sums.
        for i in range(len(sums)):
            for j in range(i + 1, len(sums)):

                s1 = sums[i]
                s2 = sums[j]

                x, y = by_sum[s1]
                u, v = by_sum[s2]

                if s1 != s2:
                    collisions.append(
                        (product, x, y, s1, u, v, s2)
                    )

                if len(collisions) >= MAX_COLLISIONS_PER_MODULUS:
                    return collisions

    return collisions


# ============================================================================
# ACTUAL PRIME COUNTEREXAMPLE
# ============================================================================

def lift_collision(m, collision):
    """
    Try to turn

        (x,y) and (u,v)

    into actual prime pairs.
    """

    product, x, y, s1, u, v, s2 = collision

    p = random_prime_with_residue(x, m)
    q = random_prime_with_residue(y, m)

    if p is None or q is None:
        return None

    P = random_prime_with_residue(u, m)
    Q = random_prime_with_residue(v, m)

    if P is None or Q is None:
        return None

    if p == q:
        return None

    if P == Q:
        return None

    n1 = p * q
    n2 = P * Q

    actual_s1 = p + q
    actual_s2 = P + Q

    if n1 % m != n2 % m:
        return None

    if actual_s1 % m == actual_s2 % m:
        return None

    return {
        "p": p,
        "q": q,
        "P": P,
        "Q": Q,
        "n1": n1,
        "n2": n2,
        "s1": actual_s1,
        "s2": actual_s2,
        "n_residue": n1 % m,
        "s1_residue": actual_s1 % m,
        "s2_residue": actual_s2 % m,
    }


# ============================================================================
# SINGLE-MODULUS EXPERIMENT
# ============================================================================

def experiment_single_modulus(r):
    m = F(r)

    print()
    print("-" * 78)
    print("r =", r, "m = F(r) =", m)
    print("-" * 78)

    units = sum(
        1 for x in range(m)
        if math.gcd(x, m) == 1
    )

    print("unit residues =", units)

    collisions = find_product_collisions(m)

    print(
        "adversarial product collisions found =",
        len(collisions)
    )

    if not collisions:
        print(
            "RESULT: no collision found in unit search."
        )
        return 0, 0

    print()
    print("FIRST RESIDUE COLLISIONS")
    print()

    for i, c in enumerate(collisions[:MAX_COLLISIONS_PER_MODULUS], 1):
        product, x, y, s1, u, v, s2 = c

        print(
            f"[{i}] product={product:>5}  "
            f"({x},{y}) -> sum={s1:>5}   "
            f"({u},{v}) -> sum={s2:>5}"
        )

    print()
    print("TRYING ACTUAL PRIME LIFTS")
    print()

    success = 0

    for i, collision in enumerate(
        collisions[:MAX_COLLISIONS_PER_MODULUS],
        1
    ):
        example = lift_collision(m, collision)

        if example is None:
            print(
                f"[{i}] prime lift: FAILED"
            )
            continue

        success += 1

        print(
            f"[{i}] prime lift: SUCCESS"
        )

        print(
            "    p,q       =",
            example["p"],
            example["q"]
        )

        print(
            "    P,Q       =",
            example["P"],
            example["Q"]
        )

        print(
            "    n1 mod m  =",
            example["n1"] % m
        )

        print(
            "    n2 mod m  =",
            example["n2"] % m
        )

        print(
            "    s1 mod m  =",
            example["s1"] % m
        )

        print(
            "    s2 mod m  =",
            example["s2"] % m
        )

        print(
            "    n residues equal =",
            example["n1"] % m == example["n2"] % m
        )

        print(
            "    s residues differ =",
            example["s1"] % m != example["s2"] % m
        )

        print(
            "    n1 bits =",
            example["n1"].bit_length(),
            "n2 bits =",
            example["n2"].bit_length()
        )

    return len(collisions), success


# ============================================================================
# GENERALIZED CRT
# ============================================================================

def crt_pair(a, m, b, n):
    """
    Solve

        x == a mod m
        x == b mod n

    returning (x, lcm(m,n)), or None if inconsistent.
    """

    g = math.gcd(m, n)

    if (b - a) % g != 0:
        return None

    m1 = m // g
    n1 = n // g

    # x = a + m*k
    #
    # m*k == b-a (mod n)
    #
    # m1*k == (b-a)/g (mod n1)

    rhs = (b - a) // g

    if n1 == 1:
        k = 0
    else:
        inv = pow(m1, -1, n1)
        k = (rhs * inv) % n1

    lcm = m * n1
    x = (a + m * k) % lcm

    return x, lcm


# ============================================================================
# COMBINED MODULUS CONSTRUCTION
# ============================================================================

def combined_modulus(rs):
    M = 1

    for r in rs:
        M = math.lcm(M, F(r))

    return M


def combined_residue_collision(rs):
    """
    Search directly modulo the combined modulus.

    We use unit residues x and choose y=1 for the first pair,
    then search another representation of the same product.

    This is deliberately small-modulus and adversarial rather than
    target-sampling based.
    """

    M = combined_modulus(rs)

    print()
    print("=" * 78)
    print("COMBINED MODULUS")
    print("=" * 78)

    print("r values       =", rs)
    print("F(r) values    =", [F(r) for r in rs])
    print("combined M     =", M)
    print("M bits         =", M.bit_length())

    # Do not attempt enormous exhaustive searches.
    if M > 2_000_000:
        print(
            "M too large for exhaustive residue search."
        )
        print(
            "Using constructive unit collision search."
        )

    # Simple universal construction:
    #
    # (x, y) and (x*t, y*t^{-1})
    #
    # have the same product modulo M.
    #
    # Their sums need not be equal.
    #
    # This is exactly the symmetry that destroys n-only recovery.

    candidates = []

    for _ in range(10000):

        # Choose small units when possible.
        x = random.randrange(1, min(M, 100000))

        if math.gcd(x, M) != 1:
            continue

        y = random.randrange(1, min(M, 100000))

        if math.gcd(y, M) != 1:
            continue

        t = random.randrange(2, min(M, 100000))

        if math.gcd(t, M) != 1:
            continue

        inv_t = pow(t, -1, M)

        u = (x * t) % M
        v = (y * inv_t) % M

        n1 = (x * y) % M
        n2 = (u * v) % M

        s1 = (x + y) % M
        s2 = (u + v) % M

        if n1 == n2 and s1 != s2:
            candidates.append(
                (x, y, u, v, n1, s1, s2)
            )

            break

    if not candidates:
        print("No constructive collision found.")
        return

    x, y, u, v, nres, s1, s2 = candidates[0]

    print()
    print("CONSTRUCTIVE COLLISION")
    print()
    print("pair A:")
    print("    x =", x)
    print("    y =", y)
    print("    xy mod M =", nres)
    print("    x+y mod M =", s1)

    print()
    print("pair B:")
    print("    u =", u)
    print("    v =", v)
    print("    uv mod M =", nres)
    print("    u+v mod M =", s2)

    print()
    print("same n residue =", (x * y) % M == (u * v) % M)
    print("different s residue =", s1 != s2)

    print()
    print(
        "INTERPRETATION:"
    )
    print(
        "  Product information alone does not determine"
    )
    print(
        "  the corresponding sum modulo the combined modulus."
    )


# ============================================================================
# INFORMATION-THEORETIC LOCAL TEST
# ============================================================================

def exact_local_ambiguity(r):
    """
    For each n residue, calculate how many different s residues can occur
    among unit factor pairs.

    This is stronger than testing a handful of actual targets.
    """

    m = F(r)

    by_product = defaultdict(set)

    units = [
        x for x in range(m)
        if math.gcd(x, m) == 1
    ]

    for x in units:
        for y in units:
            product = (x * y) % m
            summation = (x + y) % m
            by_product[product].add(summation)

    ambiguity = [len(v) for v in by_product.values()]

    max_ambiguity = max(ambiguity) if ambiguity else 0
    avg_ambiguity = (
        sum(ambiguity) / len(ambiguity)
        if ambiguity else 0
    )

    deterministic_products = sum(
        1 for v in by_product.values()
        if len(v) == 1
    )

    print(
        f"r={r:>2} "
        f"m={m:>6} "
        f"unit_products={len(by_product):>6} "
        f"deterministic={deterministic_products:>6} "
        f"max_s_ambiguity={max_ambiguity:>5} "
        f"avg_s_ambiguity={avg_ambiguity:>8.2f}"
    )

    return by_product


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA NEXT EXPERIMENT")
    print("ADVERSARIAL N-ONLY / S-RECOVERY COLLISION SEARCH")
    print("=" * 78)

    print()
    print("random seed =", SEED)
    print("R values    =", R_VALUES)
    print("targets     = NONE")
    print("CSV output  = NONE")

    # ----------------------------------------------------------------------
    # 1. Exact local information-theoretic test
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT LOCAL N -> S AMBIGUITY")
    print("=" * 78)

    print()
    print(
        "For each m=F(r), enumerate unit factor pairs."
    )

    print(
        "For each product residue n, count possible sum residues s."
    )

    print()

    local_results = {}

    for r in R_VALUES:
        local_results[r] = exact_local_ambiguity(r)

    # ----------------------------------------------------------------------
    # 2. Adversarial residue collisions
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ADVERSARIAL RESIDUE COLLISION SEARCH")
    print("=" * 78)

    total_collisions = 0
    total_prime_counterexamples = 0

    for r in R_VALUES:
        collisions, successes = experiment_single_modulus(r)

        total_collisions += collisions
        total_prime_counterexamples += successes

    # ----------------------------------------------------------------------
    # 3. Combined modulus constructions
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. COMBINED-MODULUS ADVERSARIAL TEST")
    print("=" * 78)

    for rs in COMBINED_R_SETS:
        combined_residue_collision(rs)

    # ----------------------------------------------------------------------
    # 4. Direct symmetry explanation
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ALGEBRAIC CONTROL")
    print("=" * 78)

    print()
    print("For any unit t modulo M:")
    print()
    print("    (x,y) -> (x*t, y*t^(-1))")
    print()
    print("preserves")
    print()
    print("    xy mod M")
    print()
    print("but generally changes")
    print()
    print("    x+y mod M.")
    print()

    print(
        "Therefore an n-only recovery of s would require"
    )

    print(
        "additional structure that forbids these alternative"
    )

    print(
        "factor residue pairs."
    )

    # ----------------------------------------------------------------------
    # 5. Final classification
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FINAL CLASSIFICATION")
    print("=" * 78)

    print()

    print(
        "Residue collisions found          =",
        total_collisions
    )

    print(
        "Actual prime counterexamples       =",
        total_prime_counterexamples
    )

    print()

    if total_prime_counterexamples > 0:
        print("RESULT: STRONG NEGATIVE SIGNAL")
        print()
        print(
            "We have actual semiprimes with equal n residue"
        )
        print(
            "but different s residue modulo the same F(r)."
        )
        print()
        print(
            "Therefore n mod F(r) cannot determine s mod F(r)"
        )
        print(
            "for that auxiliary modulus."
        )
        print()
        print(
            "This also explains why a low-degree polynomial"
        )
        print(
            "n -> s search is not the right next direction."
        )

    else:
        print("RESULT: INCONCLUSIVE")
        print()
        print(
            "Residue-level collisions may exist, but no actual"
        )
        print(
            "prime lift was found within the configured search."
        )
        print(
            "Increase PRIME_LIFT_ATTEMPTS before drawing conclusions."
        )

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

