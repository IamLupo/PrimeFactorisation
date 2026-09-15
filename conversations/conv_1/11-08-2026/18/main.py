#!/usr/bin/env python3

"""
KAPPA NEXT EXPERIMENTS
======================

Research target
---------------

For

    F(x) = x^2 - x + 1

and

    K3(p,q,r) = 1 - F(p)F(q)F(r),

we know that if r is controlled,

    A = F(p)F(q) = (1-K3)/F(r).

For n = p*q,

    A = u^2 - (n+1)u + n^2 - n + 1,

where

    u = p+q.

Hence

    D_u = 4A - 3(n-1)^2.

The central open question is:

    Can the required information be recovered from n alone?

This script investigates:

1. Exact identities.
2. Prime-pair datasets.
3. kappa values for controlled r.
4. Whether kappa contains information beyond n.
5. Collisions for fixed n.
6. Modular behavior.
7. gcd relationships.
8. Cubic-root-of-unity / F(x) structure.
9. Whether F(p)F(q) correlates with simple functions of n.
10. Reconstruction from an oracle.
11. n = p*q*r with known r.
12. n = p*q*r*s with known r,s.
13. Approximate/generalized controlled variables.
14. Dataset generation.

No claims of polynomial-time factoring are made.
"""

from math import gcd, isqrt
from itertools import combinations
from collections import defaultdict
import csv


# =============================================================================
# BASIC MATHEMATICS
# =============================================================================

def F(x):
    return x*x - x + 1


def kappa(xs):
    """
    K_N = 1 - product F(x_i)
    """
    product = 1
    for x in xs:
        product *= F(x)
    return 1 - product


def kappa3(p, q, r):
    return 1 - F(p) * F(q) * F(r)


def kappa4(p, q, r, s):
    return 1 - F(p) * F(q) * F(r) * F(s)


def primes_up_to(limit):
    sieve = [True] * (limit + 1)

    if limit >= 0:
        sieve[0] = False
    if limit >= 1:
        sieve[1] = False

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            for j in range(p*p, limit + 1, p):
                sieve[j] = False

    return [i for i in range(2, limit + 1) if sieve[i]]


# =============================================================================
# SECTION 1
# EXACT STRUCTURE
# =============================================================================

def section_1_exact_structure():
    print("=" * 80)
    print("1. EXACT ALGEBRAIC STRUCTURE")
    print("=" * 80)

    print("""
F(x) = x^2 - x + 1

K_N = 1 - product F(x_i)

Therefore:

K3 = 1 - F(p)F(q)F(r)

K4 = 1 - F(p)F(q)F(r)F(s)

The exact transition is

K4 - K3
    = F(p)F(q)F(r) * (1-F(s))
    = s(1-s)F(p)F(q)F(r)

and therefore

K4 - K3 = s(1-s)(1-K3).

This section is exact, not experimental.
""")


# =============================================================================
# SECTION 2
# ORIGINAL DEFINITION CHECK
# =============================================================================

def direct_kappa(xs):
    """
    Original ratio:

        P_N = product(1+x_i^3) / product(1+x_i)

    K_N = P_N - 1

    Since

        1+x^3 = (1+x)(x^2-x+1),

    this equals

        product F(x_i) - 1.

    Note:
        The research convention used above is

        K = 1 - product F(x_i).

    So the sign depends on which convention is used.
    """

    numerator = 1
    denominator = 1

    for x in xs:
        numerator *= 1 + x**3
        denominator *= 1 + x

    return numerator // denominator - 1


def section_2_definition_check():
    print("\n" + "=" * 80)
    print("2. ORIGINAL DEFINITION VS F-PRODUCT")
    print("=" * 80)

    tests = [
        [2],
        [2, 3],
        [2, 3, 5],
        [2, 3, 5, 7],
        [2, 3, 5, 7, 11],
    ]

    for xs in tests:
        direct = direct_kappa(xs)

        product_form = 1
        for x in xs:
            product_form *= F(x)

        product_form -= 1

        print(f"xs={xs}")
        print(f"  direct       = {direct}")
        print(f"  product form = {product_form}")

        assert direct == product_form

    print("\nPASS: original ratio agrees with product F(x).")


# =============================================================================
# SECTION 3
# SYMMETRIC REDUCTION
# =============================================================================

def A_from_uv(u, v):
    return u*u - u*v - u + v*v - v + 1


def section_3_symmetric_formula():
    print("\n" + "=" * 80)
    print("3. SYMMETRIC REDUCTION")
    print("=" * 80)

    print("""
For

    u = p+q
    v = pq

we have

    F(p)F(q)
      = u^2 - uv - u + v^2 - v + 1.
""")

    for p, q in [(2, 3), (5, 11), (17, 23), (31, 47)]:
        u = p + q
        v = p * q

        direct = F(p) * F(q)
        reduced = A_from_uv(u, v)

        print(
            f"p={p:3d}, q={q:3d}: "
            f"direct={direct:10d}, reduced={reduced:10d}"
        )

        assert direct == reduced

    print("\nPASS: symmetric reduction confirmed.")


# =============================================================================
# SECTION 4
# ORACLE RECOVERY
# =============================================================================

def recover_from_kappa3(n, r, K3):
    """
    Conditional factor recovery.

    Assumes:

        n = p*q
        r is known
        K3 is available.

    Then:

        A = (1-K3)/F(r)

    and solve for u=p+q.
    """

    fr = F(r)

    numerator = 1 - K3

    if numerator % fr != 0:
        return {
            "success": False,
            "reason": "K3 is inconsistent with integer A"
        }

    A = numerator // fr

    D = 4*A - 3*(n-1)**2

    if D < 0:
        return {
            "success": False,
            "reason": "negative discriminant",
            "A": A,
            "D_u": D
        }

    root = isqrt(D)

    if root * root != D:
        return {
            "success": False,
            "reason": "D_u is not a square",
            "A": A,
            "D_u": D
        }

    candidates = []

    for sign in [1, -1]:
        numerator_u = n + 1 + sign * root

        if numerator_u % 2 != 0:
            continue

        u = numerator_u // 2

        D_factor = u*u - 4*n

        if D_factor < 0:
            continue

        rf = isqrt(D_factor)

        if rf * rf != D_factor:
            continue

        a = (u - rf) // 2
        b = (u + rf) // 2

        if a * b == n:
            candidates.append(tuple(sorted((a, b))))

    candidates = sorted(set(candidates))

    return {
        "success": len(candidates) > 0,
        "A": A,
        "D_u": D,
        "factor_candidates": candidates,
        "reason": (
            "factorization recovered"
            if candidates
            else "integer discriminant but no factor pair"
        )
    }


def section_4_oracle_test():
    print("\n" + "=" * 80)
    print("4. CONDITIONAL ORACLE RECOVERY")
    print("=" * 80)

    examples = [
        (17, 23, 5),
        (31, 47, 7),
        (101, 103, 11),
        (127, 131, 13),
    ]

    for p, q, r in examples:
        n = p*q
        K = kappa3(p, q, r)

        result = recover_from_kappa3(n, r, K)

        print(
            f"p={p}, q={q}, n={n}, r={r}, "
            f"K3={K}"
        )
        print("  ", result)

        assert result["success"]

    print("\nPASS: oracle recovery works on examples.")


# =============================================================================
# SECTION 5
# EXHAUSTIVE PRIME PAIR TEST
# =============================================================================

def section_5_exhaustive(limit=100):
    print("\n" + "=" * 80)
    print("5. EXHAUSTIVE PRIME-PAIR ORACLE TEST")
    print("=" * 80)

    ps = primes_up_to(limit)

    total = 0
    success = 0

    r_values = [0, 1, 2, 3, 5, 7]

    for p, q in combinations(ps, 2):
        n = p*q

        for r in r_values:
            total += 1

            K = kappa3(p, q, r)

            result = recover_from_kappa3(n, r, K)

            if result["success"]:
                success += 1

    print(f"Prime limit : {limit}")
    print(f"Tests       : {total}")
    print(f"Successes   : {success}")
    print(f"Failures    : {total-success}")

    if total:
        print(f"Success rate: {success/total:.6f}")


# =============================================================================
# SECTION 6
# KAPPA COLLISIONS FOR FIXED n
# =============================================================================

def section_6_fixed_n_collisions():
    print("\n" + "=" * 80)
    print("6. FIXED-n FACTORIZATION / KAPPA COLLISION TEST")
    print("=" * 80)

    examples = [60, 72, 84, 90, 120, 180, 210, 360]

    for n in examples:
        print(f"\nn = {n}")

        factor_pairs = []

        for a in range(2, isqrt(n) + 1):
            if n % a == 0:
                b = n // a
                factor_pairs.append((a, b))

        values = defaultdict(list)

        r = 5

        for p, q in factor_pairs:
            K = kappa3(p, q, r)
            values[K].append((p, q))

            print(
                f"  {p:4d} * {q:4d}"
                f" -> K3={K:12d}"
            )

        collisions = {
            k: v for k, v in values.items()
            if len(v) > 1
        }

        if collisions:
            print("  COLLISIONS:")
            for K, pairs in collisions.items():
                print(f"    K={K}: {pairs}")
        else:
            print("  No K3 collisions for this n.")


# =============================================================================
# SECTION 7
# DOES KAPPA DISTINGUISH FACTORIZATIONS?
# =============================================================================

def section_7_collision_statistics(max_n=1000, r=5):
    print("\n" + "=" * 80)
    print("7. GLOBAL KAPPA COLLISION STATISTICS")
    print("=" * 80)

    total_n = 0
    n_with_collision = 0
    total_factorizations = 0

    for n in range(4, max_n + 1):

        values = defaultdict(list)

        for p in range(2, isqrt(n) + 1):
            if n % p != 0:
                continue

            q = n // p

            K = kappa3(p, q, r)

            values[K].append((p, q))

            total_factorizations += 1

        if not values:
            continue

        total_n += 1

        collisions = [
            pairs for pairs in values.values()
            if len(pairs) > 1
        ]

        if collisions:
            n_with_collision += 1

    print(f"Range: 4 .. {max_n}")
    print(f"Composite n with factor pairs: {total_n}")
    print(f"n having K3 collisions: {n_with_collision}")
    print(
        "Collision fraction: "
        f"{n_with_collision/total_n if total_n else 0:.6f}"
    )


# =============================================================================
# SECTION 8
# MODULAR STRUCTURE OF F
# =============================================================================

def roots_mod(m):
    return [x for x in range(m) if F(x) % m == 0]


def section_8_modular_F(limit=50):
    print("\n" + "=" * 80)
    print("8. MODULAR ROOT STRUCTURE OF F(x)")
    print("=" * 80)

    for m in range(2, limit + 1):
        roots = roots_mod(m)

        if roots:
            print(f"mod {m:2d}: roots = {roots}")


# =============================================================================
# SECTION 9
# CUBIC ROOT OF UNITY CONNECTION
# =============================================================================

def section_9_cubic_structure():
    print("\n" + "=" * 80)
    print("9. CUBIC-ROOT-OF-UNITY STRUCTURE")
    print("=" * 80)

    print("""
F(x) = x^2 - x + 1.

Let omega satisfy

    omega^2 + omega + 1 = 0.

Then

    F(x) = (1 + omega*x)(1 + omega^2*x).

Therefore

    product F(x_i)

is the product of two conjugate cubic-root-of-unity evaluations.

For an integer x:

    F(x) = x^2-x+1.

Modulo a prime l, F(x)=0 is equivalent to

    x^3 = -1

with x != -1.

Thus roots exist only when the relevant cubic-root-of-unity
structure exists modulo l.

We test this computationally below.
""")

    for prime in primes_up_to(100):
        roots = roots_mod(prime)

        if roots:
            print(
                f"prime={prime:3d}: "
                f"roots={roots}"
            )


# =============================================================================
# SECTION 10
# GCD STRUCTURE
# =============================================================================

def section_10_gcd_structure(limit=500):
    print("\n" + "=" * 80)
    print("10. GCD STRUCTURE")
    print("=" * 80)

    print("""
Test:

    gcd(F(p), F(q))

for distinct primes p,q.

Also test relationships involving:

    gcd(F(p)F(q), n)
    gcd(F(p)F(q), n-1)
    gcd(F(p)F(q), n+1)
""")

    interesting = []

    ps = primes_up_to(limit)

    for p, q in combinations(ps, 2):
        n = p*q

        A = F(p)*F(q)

        values = {
            "gcd(A,n)": gcd(A, n),
            "gcd(A,n-1)": gcd(A, n-1),
            "gcd(A,n+1)": gcd(A, n+1),
            "gcd(A, n*n-1)": gcd(A, n*n-1),
        }

        if any(v > 1 for v in values.values()):
            interesting.append((p, q, n, A, values))

    print(f"Prime pairs checked: {len(list(combinations(ps, 2)))}")
    print(f"Interesting cases: {len(interesting)}")

    for row in interesting[:30]:
        p, q, n, A, values = row
        print(
            f"p={p}, q={q}, n={n}, "
            f"A={A}, gcds={values}"
        )


# =============================================================================
# SECTION 11
# SIMPLE FUNCTIONS OF n
# =============================================================================

def section_11_compare_simple_functions(limit=300):
    print("\n" + "=" * 80)
    print("11. SEARCH FOR SIMPLE n-ONLY RELATIONS")
    print("=" * 80)

    print("""
For each prime pair, calculate

    A = F(p)F(q)

and compare against simple polynomial expressions in n.

We look for exact divisibility or gcd patterns rather than
assuming a relation exists.
""")

    tests = {
        "A mod n": lambda n, A: A % n,
        "A mod (n-1)": lambda n, A: A % (n-1),
        "A mod (n+1)": lambda n, A: A % (n+1),
        "A mod (n^2-1)": lambda n, A: A % (n*n-1),
        "gcd(A,n)": lambda n, A: gcd(A, n),
        "gcd(A,n-1)": lambda n, A: gcd(A, n-1),
        "gcd(A,n+1)": lambda n, A: gcd(A, n+1),
    }

    ps = primes_up_to(limit)

    counters = {name: defaultdict(int) for name in tests}

    total = 0

    for p, q in combinations(ps, 2):
        n = p*q
        A = F(p)*F(q)

        total += 1

        for name, fn in tests.items():
            value = fn(n, A)
            counters[name][value] += 1

    print(f"Prime pairs: {total}")

    for name, counter in counters.items():

        most_common = sorted(
            counter.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]

        print(f"\n{name}")
        print("  most common:")
        for value, count in most_common:
            print(f"    {value}: {count}")


# =============================================================================
# SECTION 12
# DISCRIMINANT TEST WITHOUT KAPPA
# =============================================================================

def section_12_discriminant_structure(limit=300):
    print("\n" + "=" * 80)
    print("12. DISCRIMINANT STRUCTURE")
    print("=" * 80)

    print("""
For A=F(p)F(q):

    D_u = 4A - 3(n-1)^2.

For a genuine factor pair this is exactly

    D_u = (p+q-(n+1))^2.

We verify this directly.
""")

    ps = primes_up_to(limit)

    failures = 0

    for p, q in combinations(ps, 2):
        n = p*q
        u = p+q
        A = F(p)*F(q)

        D = 4*A - 3*(n-1)**2

        expected = (u-(n+1))**2

        if D != expected:
            failures += 1

    print(f"Prime pairs checked: {len(list(combinations(ps, 2)))}")
    print(f"Failures: {failures}")

    if failures == 0:
        print("PASS: discriminant identity exact.")


# =============================================================================
# SECTION 13
# CONTROLLED r EXPERIMENT
# =============================================================================

def section_13_controlled_r():
    print("\n" + "=" * 80)
    print("13. CONTROLLED-r EXPERIMENT")
    print("=" * 80)

    p = 17
    q = 23

    n = p*q

    print(f"Hidden p={p}, q={q}, n={n}")

    for r in range(-10, 11):

        K = kappa3(p, q, r)

        A = (1-K) // F(r)

        print(
            f"r={r:3d} "
            f"F(r)={F(r):5d} "
            f"K3={K:15d} "
            f"A={A:10d}"
        )

        assert A == F(p)*F(q)


# =============================================================================
# SECTION 14
# CONTROLLED r,s FOR FOUR FACTORS
# =============================================================================

def section_14_controlled_rs():
    print("\n" + "=" * 80)
    print("14. n = p*q*r*s WITH CONTROLLED r,s")
    print("=" * 80)

    p = 11
    q = 17
    r = 5
    s = 7

    n = p*q*r*s

    K4 = kappa4(p, q, r, s)

    known = F(r)*F(s)

    A = (1-K4) // known

    reduced_n = n // (r*s)

    print(f"n={n}")
    print(f"known r={r}, s={s}")
    print(f"K4={K4}")
    print(f"F(r)F(s)={known}")
    print(f"Recovered F(p)F(q)={A}")
    print(f"Reduced pq={reduced_n}")

    assert A == F(p)*F(q)
    assert reduced_n == p*q

    print("PASS.")


# =============================================================================
# SECTION 15
# GENERAL n = p*q*r WITH KNOWN r
# =============================================================================

def section_15_pqr():
    print("\n" + "=" * 80)
    print("15. n = p*q*r WITH KNOWN r")
    print("=" * 80)

    examples = [
        (11, 17, 23),
        (13, 19, 29),
        (23, 31, 37),
        (41, 43, 47),
    ]

    for p, q, r in examples:

        n = p*q*r
        K = kappa3(p, q, r)

        reduced_n = n // r
        A = (1-K) // F(r)

        result = recover_from_kappa3(
            reduced_n,
            r,
            K
        )

        print(
            f"p={p}, q={q}, r={r}, "
            f"n={n}"
        )
        print(f"  reduced pq={reduced_n}")
        print(f"  recovered A={A}")
        print(f"  result={result}")

        assert A == F(p)*F(q)
        assert result["success"]


# =============================================================================
# SECTION 16
# MODULAR KAPPA EXPERIMENT
# =============================================================================

def section_16_modular_kappa(moduli=(2, 3, 4, 5, 7, 9, 11, 13, 19)):
    print("\n" + "=" * 80)
    print("16. MODULAR KAPPA EXPERIMENT")
    print("=" * 80)

    ps = primes_up_to(50)

    for m in moduli:

        print(f"\nmod {m}")

        for r in range(m):

            residues = set()

            for p, q in combinations(ps, 2):

                if p % m == 0 or q % m == 0:
                    continue

                K = kappa3(p, q, r)

                residues.add(K % m)

            if len(residues) <= 1:
                print(
                    f"  r={r:2d}, "
                    f"F(r)={F(r)%m:2d}, "
                    f"K residues={sorted(residues)}"
                )


# =============================================================================
# SECTION 17
# CAN MODULAR INFORMATION RECOVER FACTORS?
# =============================================================================

def section_17_modular_factor_search(limit=500, moduli=(3, 5, 7, 11, 13, 19)):
    print("\n" + "=" * 80)
    print("17. MODULAR FACTOR-RECOVERY SEARCH")
    print("=" * 80)

    print("""
For each n=pq we know n mod m.

We compare this with the hidden information:

    p mod m
    q mod m
    F(p) mod m
    F(q) mod m
    F(p)F(q) mod m.

We test whether any small modulus gives a unique
signature for the unordered factor pair.
""")

    ps = primes_up_to(limit)

    for m in moduli:

        signatures = defaultdict(set)

        for p, q in combinations(ps, 2):

            n = p*q
            A = F(p)*F(q)

            signature = (
                n % m,
                A % m
            )

            pair_signature = tuple(sorted((p % m, q % m)))

            signatures[signature].add(pair_signature)

        unique = sum(
            1 for values in signatures.values()
            if len(values) == 1
        )

        total = len(signatures)

        print(
            f"mod {m:2d}: "
            f"signatures={total:6d}, "
            f"unique={unique:6d}"
        )


# =============================================================================
# SECTION 18
# SEARCH FOR KAPPA RELATIONS BETWEEN DIFFERENT r
# =============================================================================

def section_18_multi_r_relations(limit=100):
    print("\n" + "=" * 80)
    print("18. RELATIONS BETWEEN KAPPA VALUES FOR DIFFERENT r")
    print("=" * 80)

    ps = primes_up_to(limit)

    r_values = [-2, -1, 0, 1, 2, 3, 5]

    print("""
For fixed p,q:

    1-K3(r) = A*F(r)

Therefore ratios satisfy

    (1-K3(r1))/(1-K3(r2))
        = F(r1)/F(r2).

This is exact.

We test whether eliminating A produces
new n-only information.
""")

    for p, q in list(combinations(ps, 2))[:10]:

        values = {
            r: 1-kappa3(p, q, r)
            for r in r_values
        }

        print(f"\np={p}, q={q}")

        for r in r_values:
            print(
                f"  r={r:2d}: "
                f"1-K3={values[r]}"
            )

        for r1, r2 in combinations(r_values, 2):

            lhs_num = values[r1] * F(r2)
            lhs_den = values[r2] * F(r1)

            assert lhs_num == lhs_den

    print("\nPASS: all tested cross-r relations exact.")


# =============================================================================
# SECTION 19
# SEARCH FOR SIMPLE POLYNOMIAL RELATIONS
# =============================================================================

def section_19_polynomial_residuals(limit=100):
    print("\n" + "=" * 80)
    print("19. POLYNOMIAL RESIDUAL SEARCH")
    print("=" * 80)

    print("""
We test simple candidates involving

    A = F(p)F(q)
    n = pq
    u = p+q

The goal is NOT to assume a formula,
but to identify exact low-degree residuals.
""")

    candidates = {
        "A - n^2": lambda n, u, A: A - n*n,
        "A - n(n-1)": lambda n, u, A: A - n*(n-1),
        "A - (n^2-n+1)": lambda n, u, A: A-(n*n-n+1),
        "A - (n+1)^2": lambda n, u, A: A-(n+1)**2,
        "A - u^2": lambda n, u, A: A-u*u,
        "A - n*u": lambda n, u, A: A-n*u,
    }

    ps = primes_up_to(limit)

    for name, fn in candidates.items():

        zero_count = 0
        total = 0

        for p, q in combinations(ps, 2):

            n = p*q
            u = p+q
            A = F(p)*F(q)

            total += 1

            if fn(n, u, A) == 0:
                zero_count += 1

        print(
            f"{name:25s}: "
            f"{zero_count}/{total} exact zeros"
        )


# =============================================================================
# SECTION 20
# DATASET GENERATION
# =============================================================================

def section_20_generate_dataset(
    filename="kappa_experiment_dataset.csv",
    prime_limit=200,
    r_values=(-2, -1, 0, 1, 2, 3, 5, 7, 11)
):
    print("\n" + "=" * 80)
    print("20. GENERATING DATASET")
    print("=" * 80)

    ps = primes_up_to(prime_limit)

    rows = []

    for p, q in combinations(ps, 2):

        n = p*q
        u = p+q
        A = F(p)*F(q)

        for r in r_values:

            K = kappa3(p, q, r)

            D_u = 4*A - 3*(n-1)**2

            rows.append({
                "p": p,
                "q": q,
                "n": n,
                "u": u,
                "A": A,
                "r": r,
                "F_r": F(r),
                "K3": K,
                "D_u": D_u,
                "D_u_mod_n": D_u % n,
                "A_mod_n": A % n,
                "A_mod_n_minus_1": A % (n-1),
                "A_mod_n_plus_1": A % (n+1),
                "gcd_A_n": gcd(A, n),
                "gcd_A_n_minus_1": gcd(A, n-1),
                "gcd_A_n_plus_1": gcd(A, n+1),
                "p_mod_3": p % 3,
                "q_mod_3": q % 3,
                "p_mod_7": p % 7,
                "q_mod_7": q % 7,
                "F_p_mod_7": F(p) % 7,
                "F_q_mod_7": F(q) % 7,
                "K3_mod_7": K % 7,
            })

    fieldnames = list(rows[0].keys())

    with open(filename, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )

        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {filename}")
    print(f"Rows : {len(rows)}")


# =============================================================================
# SECTION 21
# IMPORTANT COMPLEXITY CHECK
# =============================================================================

def section_21_complexity_warning():
    print("\n" + "=" * 80)
    print("21. COMPLEXITY / INFORMATION CHECK")
    print("=" * 80)

    print("""
The experiments above establish conditional algebraic recovery:

    n + r + K3  --->  p,q

or

    n + r+s + K4 ---> p,q

But the key question is whether K3/K4 can itself be computed
from n efficiently.

If computing K3 requires p and q, then the identity does not
immediately provide a new factoring algorithm.

Therefore future experiments should specifically attack:

    n ---> K3

rather than repeatedly testing:

    K3 ---> p,q

The latter direction is already solved.
""")


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 80)
    print("KAPPA CONTROLLED-VARIABLE / FACTOR-RECOVERY ANALYSIS")
    print("=" * 80)

    section_1_exact_structure()
    section_2_definition_check()
    section_3_symmetric_formula()
    section_4_oracle_test()
    section_5_exhaustive(limit=100)

    section_6_fixed_n_collisions()
    section_7_collision_statistics(max_n=1000)

    section_8_modular_F(limit=50)
    section_9_cubic_structure()

    section_10_gcd_structure(limit=300)
    section_11_compare_simple_functions(limit=200)
    section_12_discriminant_structure(limit=200)

    section_13_controlled_r()
    section_14_controlled_rs()
    section_15_pqr()

    section_16_modular_kappa()
    section_17_modular_factor_search(limit=300)

    section_18_multi_r_relations(limit=50)
    section_19_polynomial_residuals(limit=100)

    section_20_generate_dataset(
        filename="kappa_next_experiment_dataset.csv",
        prime_limit=200
    )

    section_21_complexity_warning()

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()

