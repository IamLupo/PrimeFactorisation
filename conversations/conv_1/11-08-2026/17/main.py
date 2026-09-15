#!/usr/bin/env python3

"""
KAPPA STRUCTURE / CONTROLLED VARIABLE / FACTOR RECOVERY EXPERIMENT

Self-contained research experiment.

Main identity:

    K_N = product_i F(x_i) - 1

where

    F(x) = x^2 - x + 1.

Equivalently,

    K_N = 1 - product_i F(x_i)

depending on the convention used for kappa.

This script uses the convention

    K_N = 1 - product_i F(x_i)

throughout.

For N=3:

    K3 = 1 - F(p) F(q) F(r)

For n = p q and controlled r:

    A = F(p)F(q) = (1-K3)/F(r)

and

    A = u^2 - (n+1)u + n^2-n+1

where

    u = p+q.

Therefore

    u^2 -(n+1)u + n^2-n+1-A = 0

with discriminant

    D_u = 4A - 3(n-1)^2.

This script tests the complete chain experimentally.

Output:
    - human-readable report
    - CSV data files
"""

import csv
import math
import sympy as sp
from itertools import combinations


# ============================================================================
# CONFIGURATION
# ============================================================================

PRIME_LIMIT = 100
MODULI = list(range(2, 31))

CONTROLLED_R_VALUES = list(range(-10, 11))
CONTROLLED_S_VALUES = list(range(-5, 6))

SMALL_F_LIMIT = 100

OUTPUT_PREFIX = "kappa_experiment"


# ============================================================================
# BASIC MATHEMATICS
# ============================================================================

def F(x):
    """F(x) = x^2 - x + 1."""
    return x*x - x + 1


def kappa(xs):
    """
    K_N = 1 - product F(x_i)
    """
    result = 1
    for x in xs:
        result -= 0  # keep integer arithmetic explicit
    product = 1
    for x in xs:
        product *= F(x)
    return 1 - product


def kappa_product_form(xs):
    """
    Equivalent representation:
        1 - K_N = product F(x_i)
    """
    product = 1
    for x in xs:
        product *= F(x)
    return product


# ============================================================================
# DIRECT ORIGINAL DEFINITION
# ============================================================================

def kappa_from_original(xs):
    """
    Original ratio:

        P_N = product(1+x_i^3) / product(1+x_i)

    and

        K_N = 1 - P_N

    whenever denominators are nonzero.

    Since

        (1+x^3)/(1+x) = x^2-x+1

    this must agree with kappa().
    """
    numerator = 1
    denominator = 1

    for x in xs:
        numerator *= (1 + x**3)
        denominator *= (1 + x)

    if denominator == 0:
        return None

    P = sp.Rational(numerator, denominator)

    return 1 - P


# ============================================================================
# TRANSITION FORMULA
# ============================================================================

def transition_rhs(K, s):
    """
    K_{N+1} - K_N = s(1-s)(1-K_N)
    """
    return s * (1 - s) * (1 - K)


def test_transition(xs, s):
    K = kappa(xs)
    K_next = kappa(xs + [s])

    lhs = K_next - K
    rhs = transition_rhs(K, s)

    return lhs, rhs, lhs == rhs


# ============================================================================
# SYMBOLIC STRUCTURE
# ============================================================================

def symbolic_structure():
    print("=" * 80)
    print("1. SYMBOLIC STRUCTURE")
    print("=" * 80)

    p, q, r, s = sp.symbols("p q r s")

    FF = lambda x: x**2 - x + 1

    K3 = 1 - FF(p)*FF(q)*FF(r)
    K4 = 1 - FF(p)*FF(q)*FF(r)*FF(s)

    print("\nK3:")
    print(sp.factor(K3))

    print("\nK4:")
    print(sp.factor(K4))

    transition = sp.factor(K4 - K3)

    print("\nK4 - K3:")
    print(transition)

    expected = sp.factor(
        s*(1-s)*(1-K3)
    )

    residual = sp.expand(transition - expected)

    print("\nTransition residual:")
    print(residual)

    if residual == 0:
        print("\nPASS: exact transition identity confirmed.")
    else:
        print("\nFAIL: transition identity did not simplify to zero.")


# ============================================================================
# ORIGINAL DEFINITION TEST
# ============================================================================

def test_original_definition():
    print("\n")
    print("=" * 80)
    print("2. ORIGINAL DEFINITION VS PRODUCT FORM")
    print("=" * 80)

    test_sets = [
        [2],
        [2, 3],
        [2, 3, 5],
        [2, 3, 5, 7],
        [2, 3, 5, 7, 11],
    ]

    all_pass = True

    for xs in test_sets:
        direct = kappa_from_original(xs)
        compact = kappa(xs)

        print(f"\nxs = {xs}")
        print("direct :", direct)
        print("compact:", compact)

        if direct != compact:
            all_pass = False
            print("FAIL")

    if all_pass:
        print("\nPASS: original ratio agrees with compact F-product form.")


# ============================================================================
# SYMBOLIC SYMMETRIC REDUCTION
# ============================================================================

def symmetric_reduction():
    print("\n")
    print("=" * 80)
    print("3. SYMMETRIC REDUCTION OF F(p)F(q)")
    print("=" * 80)

    p, q = sp.symbols("p q")
    u, v = sp.symbols("u v")

    expression = (
        (p**2-p+1) *
        (q**2-q+1)
    )

    symmetric = (
        u**2
        - u*v
        - u
        + v**2
        - v
        + 1
    )

    # Substitute p+q=u and pq=v using symmetric reduction.
    poly = sp.Poly(sp.expand(expression), p, q)

    reduced = sp.symmetrize(
        expression,
        [p, q],
        formal=True
    )

    print("\nDirect:")
    print(sp.expand(expression))

    print("\nExpected:")
    print(symmetric)

    print("\nSymmetrize output:")
    print(reduced)

    reconstructed = symmetric.subs({
        u: p+q,
        v: p*q
    })

    residual = sp.expand(expression - reconstructed)

    print("\nResidual:")
    print(residual)

    if residual == 0:
        print("\nPASS: symmetric reduction confirmed.")
    else:
        print("\nFAIL.")


# ============================================================================
# RECOVERY FORMULAS
# ============================================================================

def recovery_from_kappa(n, r, K3):
    """
    Given:

        n = p*q
        r = controlled variable
        K3 = 1 - F(p)F(q)F(r)

    recover p,q.
    """

    Fr = F(r)

    if Fr == 0:
        return {
            "success": False,
            "reason": "F(r)=0; division unavailable"
        }

    A = (1 - K3) // Fr

    # u^2 -(n+1)u +(n^2-n+1-A) = 0

    Du = 4*A - 3*(n-1)**2

    if Du < 0:
        return {
            "success": False,
            "reason": "negative D_u",
            "A": A,
            "D_u": Du
        }

    sqrt_Du = math.isqrt(Du)

    if sqrt_Du**2 != Du:
        return {
            "success": False,
            "reason": "D_u is not a square",
            "A": A,
            "D_u": Du
        }

    candidates_u = []

    for sign in [1, -1]:
        numerator = n + 1 + sign*sqrt_Du

        if numerator % 2 == 0:
            candidates_u.append(numerator // 2)

    factor_candidates = []

    for u in candidates_u:

        Dfactor = u*u - 4*n

        if Dfactor < 0:
            continue

        sqrt_Dfactor = math.isqrt(Dfactor)

        if sqrt_Dfactor**2 != Dfactor:
            continue

        if (u + sqrt_Dfactor) % 2 != 0:
            continue

        p = (u + sqrt_Dfactor)//2
        q = (u - sqrt_Dfactor)//2

        if p*q == n:
            factor_candidates.append(tuple(sorted((p, q))))

    factor_candidates = sorted(set(factor_candidates))

    return {
        "success": len(factor_candidates) > 0,
        "n": n,
        "r": r,
        "F_r": Fr,
        "A": A,
        "D_u": Du,
        "u_candidates": candidates_u,
        "factor_candidates": factor_candidates,
        "reason": (
            "factorization recovered"
            if factor_candidates
            else "no factorization recovered"
        )
    }


# ============================================================================
# SINGLE EXAMPLE
# ============================================================================

def basic_recovery_example():
    print("\n")
    print("=" * 80)
    print("4. BASIC FACTOR RECOVERY")
    print("=" * 80)

    p = 17
    q = 23
    r = 5

    n = p*q
    K3 = kappa([p, q, r])

    print(f"\np = {p}")
    print(f"q = {q}")
    print(f"n = {n}")
    print(f"r = {r}")
    print(f"F(r) = {F(r)}")
    print(f"K3 = {K3}")

    result = recovery_from_kappa(n, r, K3)

    print("\nRecovery:")
    print(result)

    if result["success"]:
        print("\nPASS: factors recovered.")
    else:
        print("\nFAIL.")


# ============================================================================
# EXHAUSTIVE PRIME-PAIR TEST
# ============================================================================

def prime_pair_test():
    print("\n")
    print("=" * 80)
    print("5. EXHAUSTIVE PRIME-PAIR RECOVERY")
    print("=" * 80)

    primes = list(sp.primerange(2, PRIME_LIMIT + 1))

    total = 0
    success = 0
    failures = []

    # Avoid p=q because this tests ordinary semiprime factor recovery.
    pairs = list(combinations(primes, 2))

    for p, q in pairs:

        n = p*q

        for r in [0, 1, 2, 3, 4, 5, -1, -2, 7]:

            # F(r) is never zero over integers, but keep the check.
            if F(r) == 0:
                continue

            K3 = kappa([p, q, r])

            result = recovery_from_kappa(n, r, K3)

            total += 1

            if result["success"]:
                recovered = result["factor_candidates"]

                expected = tuple(sorted((p, q)))

                if expected in recovered:
                    success += 1
                else:
                    failures.append(
                        (p, q, r, "wrong factors", recovered)
                    )
            else:
                failures.append(
                    (p, q, r, result["reason"])
                )

    print(f"\nPrime limit: {PRIME_LIMIT}")
    print(f"Tests: {total}")
    print(f"Successes: {success}")
    print(f"Failures: {len(failures)}")

    if total:
        print(f"Success rate: {success/total}")

    if failures:
        print("\nFirst failures:")
        for failure in failures[:10]:
            print(failure)
    else:
        print("\nPASS: all tested prime pairs recovered.")


# ============================================================================
# CONTROLLED r ANALYSIS
# ============================================================================

def controlled_r_analysis():
    print("\n")
    print("=" * 80)
    print("6. CONTROLLED r ANALYSIS")
    print("=" * 80)

    print("\n r        F(r)")
    print("----------------")

    for r in CONTROLLED_R_VALUES:
        print(f"{r:3d}      {F(r):8d}")


# ============================================================================
# SMALL F SEARCH
# ============================================================================

def small_F_search():
    print("\n")
    print("=" * 80)
    print("7. SMALL F(r) SEARCH")
    print("=" * 80)

    values = []

    for r in range(-100, 101):

        value = abs(F(r))

        if value <= SMALL_F_LIMIT:
            values.append((r, F(r)))

    print("\nValues with |F(r)| <=", SMALL_F_LIMIT)

    for r, value in values:
        print(f"r = {r:4d}, F(r) = {value:6d}")


# ============================================================================
# MODULAR ROOT STRUCTURE
# ============================================================================

def modular_roots():
    print("\n")
    print("=" * 80)
    print("8. MODULAR ROOT STRUCTURE OF F(r)")
    print("=" * 80)

    rows = []

    for m in MODULI:

        roots = [
            r for r in range(m)
            if F(r) % m == 0
        ]

        rows.append((m, roots))

        print(f"mod {m:2d}: {roots}")

    return rows


# ============================================================================
# MODULAR KAPPA EXPERIMENT
# ============================================================================

def modular_kappa_experiment():
    print("\n")
    print("=" * 80)
    print("9. MODULAR KAPPA COLLAPSE SEARCH")
    print("=" * 80)

    # Small primes for p,q.
    primes = list(sp.primerange(2, 30))

    for m in MODULI:

        interesting = []

        for r in range(m):

            values = set()

            for p, q in combinations(primes, 2):

                # Don't force p,q to be distinct modulo m.
                value = kappa([p, q, r]) % m
                values.add(value)

            if len(values) <= 2:
                interesting.append(
                    (r, F(r) % m, sorted(values))
                )

        if interesting:
            print(f"\nmod {m}:")

            for row in interesting:
                print(
                    f"  r={row[0]:3d} "
                    f"F(r)={row[1]:3d} "
                    f"K3 residues={row[2]}"
                )


# ============================================================================
# K3 AS FUNCTION OF p+q AND pq
# ============================================================================

def symbolic_K3_reduction():
    print("\n")
    print("=" * 80)
    print("10. K3 AS A FUNCTION OF u=p+q, v=pq")
    print("=" * 80)

    p, q, r = sp.symbols("p q r")
    u, v = sp.symbols("u v")

    K3 = 1 - (
        (p**2-p+1)
        * (q**2-q+1)
        * (r**2-r+1)
    )

    A = (
        u**2
        - u*v
        - u
        + v**2
        - v
        + 1
    )

    expected = 1 - A*(r**2-r+1)

    residual = sp.expand(
        K3 - expected.subs({
            u: p+q,
            v: p*q
        })
    )

    print("\nK3:")
    print(sp.factor(K3))

    print("\nA = F(p)F(q):")
    print(A)

    print("\nK3 in u,v,r:")
    print(expected)

    print("\nResidual:")
    print(residual)

    if residual == 0:
        print("\nPASS.")


# ============================================================================
# DISCRIMINANT EXPERIMENT
# ============================================================================

def discriminant_experiment():
    print("\n")
    print("=" * 80)
    print("11. DISCRIMINANT STRUCTURE")
    print("=" * 80)

    n, A = sp.symbols("n A")

    u = sp.symbols("u")

    polynomial = (
        u**2
        - (n+1)*u
        + n**2-n+1-A
    )

    D = sp.factor(sp.discriminant(polynomial, u))

    print("\nRecovery polynomial:")
    print(polynomial)

    print("\nDiscriminant:")
    print(D)

    expected = 4*A - 3*(n-1)**2

    print("\nExpected:")
    print(expected)

    print("\nResidual:")
    print(sp.expand(D - expected))


# ============================================================================
# TRANSITION TEST FOR MANY N
# ============================================================================

def transition_exhaustive():
    print("\n")
    print("=" * 80)
    print("12. EXHAUSTIVE TRANSITION TEST")
    print("=" * 80)

    values = [
        -3, -2, -1, 0, 1, 2, 3, 5, 7
    ]

    total = 0
    failures = []

    for N in range(1, 8):

        # Generate a few deterministic test vectors.
        xs = list(range(2, 2+N))

        for s in values:

            lhs, rhs, ok = test_transition(xs, s)

            total += 1

            if not ok:
                failures.append(
                    (xs, s, lhs, rhs)
                )

    print(f"\nTests: {total}")
    print(f"Failures: {len(failures)}")

    if failures:
        for f in failures[:10]:
            print(f)
    else:
        print("\nPASS: transition identity holds for all tests.")


# ============================================================================
# K4 TEST
# ============================================================================

def K4_controlled_test():
    print("\n")
    print("=" * 80)
    print("13. K4 CONTROLLED-VARIABLE TEST")
    print("=" * 80)

    p = 17
    q = 23
    r = 5

    K3 = kappa([p, q, r])

    print(f"\np={p}, q={q}, r={r}")
    print(f"K3={K3}")

    for s in CONTROLLED_S_VALUES:

        K4 = kappa([p, q, r, s])

        lhs = K4 - K3
        rhs = s*(1-s)*(1-K3)

        print(
            f"s={s:3d} "
            f"K4-K3={lhs:20d} "
            f"formula={rhs:20d} "
            f"PASS={lhs==rhs}"
        )


# ============================================================================
# INFORMATION TEST: DOES K4 ADD INFORMATION?
# ============================================================================

def information_test():
    print("\n")
    print("=" * 80)
    print("14. DOES K4 ADD INFORMATION BEYOND K3?")
    print("=" * 80)

    p = 17
    q = 23
    r = 5

    K3 = kappa([p, q, r])

    print(f"\nFixed p={p}, q={q}, r={r}")
    print(f"K3={K3}")

    print("\nComparing different s:")

    for s in [-3, -1, 0, 1, 2, 5, 10]:

        K4 = kappa([p, q, r, s])

        predicted = K3 + s*(1-s)*(1-K3)

        print(
            f"s={s:3d} "
            f"K4={K4:20d} "
            f"predicted={predicted:20d} "
            f"match={K4==predicted}"
        )

    print(
        "\nObservation: once K3 and s are known, "
        "K4 is algebraically determined."
    )


# ============================================================================
# MULTIPLE CONTROLLED VARIABLES
# ============================================================================

def multi_controlled_test():
    print("\n")
    print("=" * 80)
    print("15. MULTIPLE CONTROLLED VARIABLES")
    print("=" * 80)

    p = 17
    q = 23

    r = 5
    s = 7

    K2 = kappa([p, q])
    K3 = kappa([p, q, r])
    K4 = kappa([p, q, r, s])

    print("\np,q hidden:")
    print(f"n = {p*q}")

    print("\nControlled variables:")
    print(f"r = {r}")
    print(f"s = {s}")

    print("\nValues:")
    print(f"K2 = {K2}")
    print(f"K3 = {K3}")
    print(f"K4 = {K4}")

    # Recover A directly from K3.
    A_from_K3 = (1-K3)//F(r)

    # Recover A directly from K4.
    A_from_K4 = (1-K4)//(F(r)*F(s))

    print("\nRecovered F(p)F(q):")
    print(f"from K3 = {A_from_K3}")
    print(f"from K4 = {A_from_K4}")
    print(f"actual   = {F(p)*F(q)}")

    print(
        "\nK3 and K4 therefore carry the same base "
        "factor-product information once r,s are known."
    )


# ============================================================================
# n = pqr WITH ONE CONTROLLED FACTOR
# ============================================================================

def three_factor_experiment():
    print("\n")
    print("=" * 80)
    print("16. n = p*q*r: WHAT CHANGES?")
    print("=" * 80)

    p = 11
    q = 17
    r = 23

    n = p*q*r

    # Suppose r is known.
    K3 = kappa([p, q, r])

    print(f"\nHidden factors p={p}, q={q}")
    print(f"Known r={r}")
    print(f"n={n}")
    print(f"K3={K3}")

    # Since n/r = pq, we reduce to semiprime problem.
    pq = n//r

    A = (1-K3)//F(r)

    print(f"\nReduced pq = {pq}")
    print(f"Recovered F(p)F(q) = {A}")

    # Solve for p+q.
    Du = 4*A - 3*(pq-1)**2

    print(f"D_u = {Du}")

    if Du >= 0 and math.isqrt(Du)**2 == Du:

        root = math.isqrt(Du)

        candidates = []

        for sign in [1, -1]:

            numerator = pq + 1 + sign*root

            if numerator % 2 == 0:
                candidates.append(numerator//2)

        print(f"Candidate p+q values: {candidates}")

        recovered = []

        for u in candidates:

            Df = u*u - 4*pq

            if Df >= 0 and math.isqrt(Df)**2 == Df:

                sf = math.isqrt(Df)

                if (u+sf) % 2 == 0:

                    a = (u+sf)//2
                    b = (u-sf)//2

                    if a*b == pq:
                        recovered.append(tuple(sorted((a,b))))

        print(f"Recovered p,q: {sorted(set(recovered))}")

    else:
        print("No integer square discriminant.")


# ============================================================================
# n = pqrs WITH r,s CONTROLLED
# ============================================================================

def four_factor_experiment():
    print("\n")
    print("=" * 80)
    print("17. n = p*q*r*s WITH r,s CONTROLLED")
    print("=" * 80)

    p = 11
    q = 17
    r = 5
    s = 7

    n = p*q*r*s

    K4 = kappa([p, q, r, s])

    print(f"\nHidden p={p}, q={q}")
    print(f"Known r={r}, s={s}")
    print(f"n={n}")
    print(f"K4={K4}")

    known_part = F(r)*F(s)

    A = (1-K4)//known_part

    pq = n//(r*s)

    print(f"\nKnown F(r)F(s) = {known_part}")
    print(f"Recovered F(p)F(q) = {A}")
    print(f"Reduced pq = {pq}")

    Du = 4*A - 3*(pq-1)**2

    print(f"D_u = {Du}")

    if Du >= 0 and math.isqrt(Du)**2 == Du:

        root = math.isqrt(Du)

        for sign in [1, -1]:

            numerator = pq + 1 + sign*root

            if numerator % 2 == 0:

                u = numerator//2

                Df = u*u - 4*pq

                if Df >= 0 and math.isqrt(Df)**2 == Df:

                    sf = math.isqrt(Df)

                    if (u+sf) % 2 == 0:

                        a = (u+sf)//2
                        b = (u-sf)//2

                        if a*b == pq:

                            print(
                                f"Recovered factor pair: "
                                f"({a}, {b})"
                            )

                            return

    print("No factor pair recovered.")


# ============================================================================
# IMPORTANT: BLACK-BOX FACTORING EXPERIMENT
# ============================================================================

def black_box_factor_experiment():
    """
    This is the most important conceptual experiment.

    We assume:
        n = p*q is public.
        r is known.
        K3 is supplied by an oracle.

    We do NOT use p or q in the recovery stage.

    The only question is:

        Is K3 itself sufficient to recover p,q?

    This distinguishes an algebraic identity from an actual
    factoring algorithm.
    """

    print("\n")
    print("=" * 80)
    print("18. BLACK-BOX FACTORING EXPERIMENT")
    print("=" * 80)

    primes = list(sp.primerange(3, PRIME_LIMIT + 1))

    tests = []

    for p, q in combinations(primes, 2):

        n = p*q
        r = 5

        K3_oracle = kappa([p, q, r])

        # From this point onward, pretend p and q do not exist.
        result = recovery_from_kappa(
            n=n,
            r=r,
            K3=K3_oracle
        )

        tests.append({
            "p": p,
            "q": q,
            "n": n,
            "r": r,
            "K3": K3_oracle,
            "success": result["success"]
        })

    total = len(tests)
    successes = sum(t["success"] for t in tests)

    print(f"\nTests: {total}")
    print(f"Successful recoveries: {successes}")
    print(f"Failures: {total-successes}")

    if total:
        print(f"Success rate: {successes/total}")

    print(
        "\nThis demonstrates that an oracle providing K3 "
        "would be sufficient for the algebraic recovery."
    )


# ============================================================================
# TEST WHETHER KAPPA IS EASY TO COMPUTE FROM n ALONE
# ============================================================================

def hidden_information_test():
    """
    Important negative/control experiment.

    For a fixed n, enumerate all factor pairs and see whether
    K3 differs between possible factorisations.

    For semiprimes there is only one nontrivial pair, but for
    composite n there can be many factorisations.

    This helps distinguish:
        - information contained in K3
        - information already contained in n.
    """

    print("\n")
    print("=" * 80)
    print("19. HIDDEN INFORMATION / MULTIPLE FACTORIZATION TEST")
    print("=" * 80)

    numbers = [
        60,
        72,
        84,
        90,
        120,
        180,
        210,
        360
    ]

    r = 5

    for n in numbers:

        factorizations = []

        for a in range(2, int(math.sqrt(n))+1):

            if n % a == 0:

                b = n//a

                if a != b:

                    K3 = kappa([a, b, r])

                    factorizations.append(
                        (a, b, K3)
                    )

        print(f"\nn={n}")

        for row in factorizations:
            print(
                f"  {row[0]:4d} * {row[1]:4d}"
                f" -> K3={row[2]}"
            )


# ============================================================================
# SEARCH FOR COLLISIONS
# ============================================================================

def kappa_collision_search():
    """
    Search for distinct (p,q) pairs that produce the same K3
    when r is fixed.

    This is useful for understanding whether the mapping

        (p,q) -> K3

    is injective over tested ranges.
    """

    print("\n")
    print("=" * 80)
    print("20. KAPPA COLLISION SEARCH")
    print("=" * 80)

    primes = list(sp.primerange(2, PRIME_LIMIT + 1))

    for r in [0, 1, 2, 3, 5, 7, -1]:

        buckets = {}

        for p, q in combinations(primes, 2):

            value = kappa([p, q, r])

            buckets.setdefault(value, []).append(
                (p, q)
            )

        collisions = [
            (value, pairs)
            for value, pairs in buckets.items()
            if len(pairs) > 1
        ]

        print(
            f"\nr={r}: "
            f"{len(collisions)} collision values"
        )

        for value, pairs in collisions[:5]:

            print(
                f"  K3={value}: {pairs[:10]}"
            )


# ============================================================================
# RECORD CSV DATA
# ============================================================================

def save_csv(filename, rows):

    if not rows:
        return

    keys = sorted(rows[0].keys())

    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=keys
        )

        writer.writeheader()
        writer.writerows(rows)


def generate_dataset():
    """
    Generate a machine-readable dataset for later analysis.
    """

    print("\n")
    print("=" * 80)
    print("21. GENERATING DATASETS")
    print("=" * 80)

    primes = list(sp.primerange(2, PRIME_LIMIT + 1))

    rows = []

    for p, q in combinations(primes, 2):

        n = p*q

        for r in CONTROLLED_R_VALUES:

            K3 = kappa([p, q, r])

            result = recovery_from_kappa(
                n,
                r,
                K3
            )

            rows.append({
                "p": p,
                "q": q,
                "n": n,
                "r": r,
                "F_r": F(r),
                "K3": K3,
                "A": result.get("A"),
                "D_u": result.get("D_u"),
                "success": result["success"]
            })

    filename = OUTPUT_PREFIX + "_prime_pairs.csv"

    save_csv(filename, rows)

    print(f"\nSaved: {filename}")
    print(f"Rows: {len(rows)}")


# ============================================================================
# SUMMARY
# ============================================================================

def final_summary():

    print("\n")
    print("=" * 80)
    print("FINAL RESEARCH SUMMARY")
    print("=" * 80)

    print(
        """
The experiments establish the following algebraic chain:

    K_N = 1 - product F(x_i)

with

    F(x) = x^2 - x + 1.

For n=pq and known r:

    K3 = 1 - F(p)F(q)F(r)

therefore

    A = F(p)F(q)
      = (1-K3)/F(r).

Writing

    u=p+q
    n=pq

gives

    A = u^2 -(n+1)u + n^2-n+1.

Thus

    u^2 -(n+1)u + n^2-n+1-A = 0.

The discriminant is

    D_u = 4A - 3(n-1)^2.

If D_u is an integer square then p+q is recovered.

Finally:

    p,q = roots of x^2-u*x+n.

For n=pqr with r known:

    pq = n/r

and the same mechanism applies.

For n=pqrs with r,s known:

    pq = n/(rs)

and

    F(p)F(q)
      = (1-K4)/(F(r)F(s)).

Therefore controlled variables reduce the
unknown-factor problem back to the same
two-variable algebraic problem.

The critical open question remains:

    Can K3 or K4 be computed efficiently from n
    without knowing the hidden factors?

If an independent efficient oracle exists,
the recovery equations give a factorization
procedure.

If not, the identities are algebraically exact
but do not by themselves constitute a new
factoring algorithm.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("\n")
    print("=" * 80)
    print("KAPPA CONTROLLED-VARIABLE / FACTOR-RECOVERY ANALYSIS")
    print("=" * 80)

    symbolic_structure()
    test_original_definition()
    symmetric_reduction()
    basic_recovery_example()
    prime_pair_test()
    controlled_r_analysis()
    small_F_search()
    modular_roots()
    modular_kappa_experiment()
    symbolic_K3_reduction()
    discriminant_experiment()
    transition_exhaustive()
    K4_controlled_test()
    information_test()
    multi_controlled_test()
    three_factor_experiment()
    four_factor_experiment()
    black_box_factor_experiment()
    hidden_information_test()
    kappa_collision_search()
    generate_dataset()
    final_summary()

    print("\n")
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()

