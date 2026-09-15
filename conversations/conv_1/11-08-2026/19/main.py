#!/usr/bin/env python3
"""
KAPPA CONTROLLED-VARIABLE / FACTOR-RECOVERY EXPERIMENT

Core definitions
----------------
    F(x) = x^2 - x + 1

    K_N = 1 - product(F(x_i))

For n = p*q and known r:

    K3 = 1 - F(p)F(q)F(r)

For n = p*q*r*s and known r,s:

    K4 = 1 - F(p)F(q)F(r)F(s)

The script tests the algebraic identities, conditional oracle
factor recovery, collisions, modular structure, and searches for
simple n-only relations.

IMPORTANT:
This script does NOT assume that K3 can be efficiently computed
from n. The central open computational question is whether there
is an efficient n -> K3 mechanism that does not already require
factoring n.
"""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from itertools import combinations, product
from pathlib import Path


# ============================================================================
# CONFIGURATION
# ============================================================================

PRIME_LIMIT = 300
COLLISION_LIMIT = 5000
MODULUS_LIMIT = 150
CONTROLLED_R_MIN = -25
CONTROLLED_R_MAX = 25

OUTPUT_DIR = Path("kappa_results")
OUTPUT_DIR.mkdir(exist_ok=True)


# ============================================================================
# BASIC DEFINITIONS
# ============================================================================

def F(x: int) -> int:
    """F(x) = x^2 - x + 1."""
    return x * x - x + 1


def kappa(xs) -> int:
    """K_N = 1 - product F(x_i)."""
    prod_value = 1
    for x in xs:
        prod_value *= F(x)
    return 1 - prod_value


def primes_upto(limit: int) -> list[int]:
    """Simple deterministic prime generator."""
    primes = []

    for x in range(2, limit + 1):
        is_prime = True

        for p in primes:
            if p * p > x:
                break
            if x % p == 0:
                is_prime = False
                break

        if is_prime:
            primes.append(x)

    return primes


def factor_pairs(n: int):
    """Return unordered positive factor pairs excluding 1*n."""
    pairs = []

    for d in range(2, math.isqrt(n) + 1):
        if n % d == 0:
            pairs.append((d, n // d))

    return pairs


def write_csv(filename: str, rows: list[dict]):
    """Write rows to the result directory."""
    if not rows:
        return

    path = OUTPUT_DIR / filename

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {path} ({len(rows)} rows)")


# ============================================================================
# 1. EXACT PRODUCT IDENTITY
# ============================================================================

def direct_product(xs):
    result = 1
    for x in xs:
        result *= F(x)
    return result


def test_product_identity():
    print("\n" + "=" * 80)
    print("1. EXACT PRODUCT IDENTITY")
    print("=" * 80)

    test_sets = [
        [2],
        [2, 3],
        [2, 3, 5],
        [2, 3, 5, 7],
        [2, 3, 5, 7, 11],
    ]

    for xs in test_sets:
        direct = direct_product(xs)
        K = kappa(xs)

        assert K == 1 - direct

        print(
            f"xs={str(xs):<20} "
            f"product={direct:<15} "
            f"K={K}"
        )

    print("PASS: K_N = 1 - product F(x_i)")


# ============================================================================
# 2. EXACT K4-K3 TRANSITION
# ============================================================================

def test_transition_identity():
    print("\n" + "=" * 80)
    print("2. EXACT K4-K3 TRANSITION")
    print("=" * 80)

    tests = [
        (17, 23, 5, -3),
        (17, 23, 5, 0),
        (17, 23, 5, 1),
        (17, 23, 5, 7),
        (11, 17, 23, 13),
    ]

    for p, q, r, s in tests:
        K3 = kappa([p, q, r])
        K4 = kappa([p, q, r, s])

        lhs = K4 - K3
        rhs = s * (1 - s) * (1 - K3)

        assert lhs == rhs

        print(
            f"p={p}, q={q}, r={r}, s={s}: "
            f"K4-K3={lhs}, formula={rhs}"
        )

    print("PASS: K4-K3 = s(1-s)(1-K3)")


# ============================================================================
# 3. SYMMETRIC REDUCTION
# ============================================================================

def symmetric_A(u: int, v: int) -> int:
    """
    F(p)F(q) expressed using

        u = p+q
        v = pq
    """
    return u*u - u*v - u + v*v - v + 1


def test_symmetric_reduction():
    print("\n" + "=" * 80)
    print("3. SYMMETRIC REDUCTION")
    print("=" * 80)

    tests = [
        (2, 3),
        (5, 11),
        (17, 23),
        (31, 47),
        (101, 103),
    ]

    for p, q in tests:
        u = p + q
        v = p * q

        direct = F(p) * F(q)
        reduced = symmetric_A(u, v)

        assert direct == reduced

        print(
            f"p={p:4}, q={q:4}: "
            f"direct={direct:12}, reduced={reduced:12}"
        )

    print("PASS: symmetric reduction confirmed")


# ============================================================================
# 4. DISCRIMINANT IDENTITY
# ============================================================================

def discriminant_from_A(A: int, n: int) -> int:
    """
    Recovery polynomial:

        u^2 - (n+1)u + n^2-n+1-A = 0

    Discriminant:

        D = 4A - 3(n-1)^2
    """
    return 4 * A - 3 * (n - 1) ** 2


def test_discriminant():
    print("\n" + "=" * 80)
    print("4. DISCRIMINANT STRUCTURE")
    print("=" * 80)

    rows = []
    failures = 0

    ps = primes_upto(PRIME_LIMIT)

    for p, q in combinations(ps, 2):
        n = p * q
        u = p + q
        A = F(p) * F(q)

        D = discriminant_from_A(A, n)

        # IMPORTANT:
        # D = (2u - (n+1))^2
        expected = (2*u - (n + 1)) ** 2

        ok = D == expected

        if not ok:
            failures += 1

        rows.append({
            "p": p,
            "q": q,
            "n": n,
            "u": u,
            "A": A,
            "D": D,
            "expected_D": expected,
            "pass": ok,
        })

    print(f"Prime pairs checked: {len(rows)}")
    print(f"Failures: {failures}")

    assert failures == 0

    write_csv("discriminant_identity.csv", rows)

    print("PASS: D = (2u-(n+1))^2")


# ============================================================================
# 5. CONDITIONAL ORACLE RECOVERY
# ============================================================================

def recover_from_oracle(n: int, r: int, K3: int) -> dict:
    """
    Given n=pq, known r, and an oracle K3:

        A = (1-K3)/F(r)

    Then solve for u=p+q and finally p,q.
    """

    fr = F(r)

    if fr == 0:
        return {
            "success": False,
            "reason": "F(r)=0",
        }

    numerator = 1 - K3

    if numerator % fr != 0:
        return {
            "success": False,
            "reason": "A is not integral",
        }

    A = numerator // fr

    D = discriminant_from_A(A, n)

    if D < 0:
        return {
            "success": False,
            "A": A,
            "D_u": D,
            "factor_candidates": [],
            "reason": "negative discriminant",
        }

    sqrt_D = math.isqrt(D)

    if sqrt_D * sqrt_D != D:
        return {
            "success": False,
            "A": A,
            "D_u": D,
            "factor_candidates": [],
            "reason": "discriminant not square",
        }

    candidates = []

    # u = ((n+1) +/- sqrt(D))/2
    for sign in (+1, -1):
        numerator_u = n + 1 + sign * sqrt_D

        if numerator_u % 2 != 0:
            continue

        u = numerator_u // 2

        disc_pq = u*u - 4*n

        if disc_pq < 0:
            continue

        sqrt_pq = math.isqrt(disc_pq)

        if sqrt_pq * sqrt_pq != disc_pq:
            continue

        if (u + sqrt_pq) % 2 != 0:
            continue

        p = (u + sqrt_pq) // 2
        q = (u - sqrt_pq) // 2

        if p > 0 and q > 0 and p*q == n:
            candidates.append(tuple(sorted((p, q))))

    candidates = sorted(set(candidates))

    return {
        "success": bool(candidates),
        "A": A,
        "D_u": D,
        "factor_candidates": candidates,
        "reason": (
            "factorization recovered"
            if candidates
            else "no valid factor pair"
        ),
    }


def test_oracle_recovery():
    print("\n" + "=" * 80)
    print("5. CONDITIONAL ORACLE RECOVERY")
    print("=" * 80)

    examples = [
        (17, 23, 5),
        (31, 47, 7),
        (101, 103, 11),
        (127, 131, 13),
    ]

    for p, q, r in examples:
        n = p*q
        K3 = kappa([p, q, r])

        result = recover_from_oracle(n, r, K3)

        print(
            f"p={p}, q={q}, n={n}, r={r}, "
            f"K3={K3}"
        )
        print("   ", result)

        assert result["success"]
        assert (p, q) in result["factor_candidates"]

    print("PASS: oracle recovery works")


# ============================================================================
# 6. EXHAUSTIVE PRIME-PAIR ORACLE TEST
# ============================================================================

def exhaustive_oracle_test():
    print("\n" + "=" * 80)
    print("6. EXHAUSTIVE PRIME-PAIR ORACLE TEST")
    print("=" * 80)

    ps = primes_upto(PRIME_LIMIT)

    rows = []
    failures = 0

    # r=0 is especially interesting because F(0)=1.
    controlled_values = [0, 1, 2, 5, 7]

    for p, q in combinations(ps, 2):
        for r in controlled_values:
            n = p*q
            K3 = kappa([p, q, r])

            result = recover_from_oracle(n, r, K3)

            ok = (
                result["success"]
                and (p, q) in result["factor_candidates"]
            )

            if not ok:
                failures += 1

            rows.append({
                "p": p,
                "q": q,
                "n": n,
                "r": r,
                "K3": K3,
                "A": result.get("A"),
                "D_u": result.get("D_u"),
                "recovered": ok,
            })

    print(f"Prime limit : {PRIME_LIMIT}")
    print(f"Tests       : {len(rows)}")
    print(f"Failures    : {failures}")
    print(
        f"Success rate: "
        f"{(len(rows)-failures)/len(rows):.6f}"
    )

    assert failures == 0

    write_csv("oracle_recovery.csv", rows)

    print("PASS: all tested oracle recoveries")


# ============================================================================
# 7. K3 COLLISION SEARCH
# ============================================================================

def collision_search():
    print("\n" + "=" * 80)
    print("7. GLOBAL KAPPA COLLISION SEARCH")
    print("=" * 80)

    collision_rows = []
    summary_rows = []

    collision_ns = []

    for n in range(4, COLLISION_LIMIT + 1):
        pairs = factor_pairs(n)

        if len(pairs) < 2:
            continue

        groups = defaultdict(list)

        # r=0 => F(r)=1 and K3 = 1-F(p)F(q)
        for p, q in pairs:
            K3 = kappa([p, q, 0])
            groups[K3].append((p, q))

        collisions = {
            value: pair_list
            for value, pair_list in groups.items()
            if len(pair_list) > 1
        }

        summary_rows.append({
            "n": n,
            "factor_pair_count": len(pairs),
            "collision": bool(collisions),
            "collision_group_count": len(collisions),
        })

        if collisions:
            collision_ns.append(n)

            for K3, pair_list in collisions.items():
                collision_rows.append({
                    "n": n,
                    "K3": K3,
                    "pairs": str(pair_list),
                    "pair_count": len(pair_list),
                })

    composite_count = len(summary_rows)

    print(f"Range: 4 .. {COLLISION_LIMIT}")
    print(f"Composite n with >=2 factor pairs: {composite_count}")
    print(f"n having K3 collisions: {len(collision_ns)}")

    if composite_count:
        print(
            f"Collision fraction: "
            f"{len(collision_ns)/composite_count:.8f}"
        )

    print("Collision n:", collision_ns[:100])

    write_csv("collision_summary.csv", summary_rows)
    write_csv("k3_collisions.csv", collision_rows)


# ============================================================================
# 8. MODULAR ROOT STRUCTURE
# ============================================================================

def modular_root_search():
    print("\n" + "=" * 80)
    print("8. MODULAR ROOT STRUCTURE OF F")
    print("=" * 80)

    rows = []

    for m in range(2, MODULUS_LIMIT + 1):
        roots = [
            x for x in range(m)
            if F(x) % m == 0
        ]

        if roots:
            rows.append({
                "modulus": m,
                "roots": str(roots),
                "root_count": len(roots),
            })

            print(
                f"mod {m:3}: roots={roots}"
            )

    write_csv("modular_roots.csv", rows)


# ============================================================================
# 9. PRIME MODULI / CUBIC ROOT STRUCTURE
# ============================================================================

def prime_modular_roots():
    print("\n" + "=" * 80)
    print("9. PRIME MODULI AND CUBIC-ROOT STRUCTURE")
    print("=" * 80)

    rows = []

    for p in primes_upto(MODULUS_LIMIT):
        roots = [
            x for x in range(p)
            if F(x) % p == 0
        ]

        rows.append({
            "prime": p,
            "roots": str(roots),
            "root_count": len(roots),
        })

        if roots:
            print(
                f"prime={p:3}: roots={roots}"
            )

    write_csv("prime_modular_roots.csv", rows)


# ============================================================================
# 10. CONTROLLED-r EXPERIMENT
# ============================================================================

def controlled_r_experiment():
    print("\n" + "=" * 80)
    print("10. CONTROLLED-r EXPERIMENT")
    print("=" * 80)

    p, q = 17, 23
    A_actual = F(p)*F(q)

    rows = []

    for r in range(CONTROLLED_R_MIN, CONTROLLED_R_MAX + 1):
        K3 = kappa([p, q, r])

        A_recovered = (1-K3)//F(r)

        assert A_recovered == A_actual

        rows.append({
            "p": p,
            "q": q,
            "r": r,
            "F_r": F(r),
            "K3": K3,
            "A_actual": A_actual,
            "A_recovered": A_recovered,
            "match": A_actual == A_recovered,
        })

        print(
            f"r={r:4}, "
            f"F(r)={F(r):8}, "
            f"K3={K3:15}, "
            f"A={A_recovered}"
        )

    write_csv("controlled_r.csv", rows)

    print("PASS: A recovered for every tested r")


# ============================================================================
# 11. CONTROLLED s / K4 EXPERIMENT
# ============================================================================

def controlled_s_experiment():
    print("\n" + "=" * 80)
    print("11. CONTROLLED-s / K4 EXPERIMENT")
    print("=" * 80)

    p, q, r = 17, 23, 5
    K3 = kappa([p, q, r])

    rows = []

    for s in range(-15, 16):
        K4 = kappa([p, q, r, s])

        predicted_difference = s*(1-s)*(1-K3)
        predicted_K4 = K3 + predicted_difference

        ok = K4 == predicted_K4

        assert ok

        rows.append({
            "p": p,
            "q": q,
            "r": r,
            "s": s,
            "K3": K3,
            "K4": K4,
            "predicted_K4": predicted_K4,
            "difference": K4-K3,
            "formula_difference": predicted_difference,
            "pass": ok,
        })

        print(
            f"s={s:4}: "
            f"K4={K4:16}, "
            f"predicted={predicted_K4:16}, "
            f"PASS={ok}"
        )

    write_csv("controlled_s.csv", rows)


# ============================================================================
# 12. MULTIPLE CONTROLLED VARIABLES
# ============================================================================

def controlled_rs_recovery():
    print("\n" + "=" * 80)
    print("12. n = p*q*r*s WITH CONTROLLED r,s")
    print("=" * 80)

    examples = [
        (11, 17, 5, 7),
        (13, 19, 7, 11),
        (23, 31, 5, 13),
        (41, 43, 7, 17),
    ]

    rows = []

    for p, q, r, s in examples:
        n = p*q*r*s
        K4 = kappa([p, q, r, s])

        known = F(r)*F(s)
        numerator = 1-K4

        A = numerator // known

        expected = F(p)*F(q)

        pq = n // (r*s)

        ok = (
            n % (r*s) == 0
            and A == expected
            and pq == p*q
        )

        assert ok

        rows.append({
            "p": p,
            "q": q,
            "r": r,
            "s": s,
            "n": n,
            "K4": K4,
            "F(r)F(s)": known,
            "A_recovered": A,
            "A_expected": expected,
            "pq_recovered": pq,
            "pass": ok,
        })

        print(
            f"n={n}, r={r}, s={s}, "
            f"K4={K4}, A={A}, PASS={ok}"
        )

    write_csv("controlled_rs_recovery.csv", rows)


# ============================================================================
# 13. n = p*q*r WITH KNOWN r
# ============================================================================

def known_r_three_factor_test():
    print("\n" + "=" * 80)
    print("13. n = p*q*r WITH KNOWN r")
    print("=" * 80)

    examples = [
        (11, 17, 23),
        (13, 19, 29),
        (23, 31, 37),
        (41, 43, 47),
    ]

    rows = []

    for p, q, r in examples:
        n = p*q*r
        K3 = kappa([p, q, r])

        reduced_n = n//r
        result = recover_from_oracle(reduced_n, r, K3)

        ok = (
            result["success"]
            and (p,q) in result["factor_candidates"]
        )

        assert ok

        rows.append({
            "p": p,
            "q": q,
            "r": r,
            "n": n,
            "K3": K3,
            "pq": reduced_n,
            "A": result["A"],
            "D_u": result["D_u"],
            "recovered": ok,
        })

        print(
            f"p={p}, q={q}, r={r}, n={n}, "
            f"pq={reduced_n}, "
            f"A={result['A']}, "
            f"factors={result['factor_candidates']}"
        )

    write_csv("known_r_three_factor.csv", rows)


# ============================================================================
# 14. CROSS-r ELIMINATION
# ============================================================================

def cross_r_test():
    print("\n" + "=" * 80)
    print("14. CROSS-r ELIMINATION")
    print("=" * 80)

    pairs = [
        (2,3),
        (2,5),
        (17,23),
        (31,47),
        (101,103),
    ]

    r_values = [-3, -2, -1, 0, 1, 2, 3, 5, 7]

    rows = []

    for p, q in pairs:
        A = F(p)*F(q)

        for r1, r2 in combinations(r_values, 2):
            K1 = kappa([p,q,r1])
            K2 = kappa([p,q,r2])

            lhs = (1-K1)*F(r2)
            rhs = (1-K2)*F(r1)

            ok = lhs == rhs

            assert ok

            rows.append({
                "p": p,
                "q": q,
                "r1": r1,
                "r2": r2,
                "K3_r1": K1,
                "K3_r2": K2,
                "lhs": lhs,
                "rhs": rhs,
                "A": A,
                "pass": ok,
            })

    print(f"Cross-r tests: {len(rows)}")
    print("PASS: (1-K3(r1))F(r2) = (1-K3(r2))F(r1)")

    write_csv("cross_r_relations.csv", rows)


# ============================================================================
# 15. GCD STRUCTURE
# ============================================================================

def gcd_structure_test():
    print("\n" + "=" * 80)
    print("15. GCD STRUCTURE")
    print("=" * 80)

    ps = primes_upto(PRIME_LIMIT)

    rows = []

    for p, q in combinations(ps, 2):
        n = p*q
        A = F(p)*F(q)

        values = {
            "gcd_A_n": math.gcd(A, n),
            "gcd_A_n_minus_1": math.gcd(A, n-1),
            "gcd_A_n_plus_1": math.gcd(A, n+1),
            "gcd_A_n2_minus_1": math.gcd(A, n*n-1),
        }

        rows.append({
            "p": p,
            "q": q,
            "n": n,
            "A": A,
            **values,
        })

    print(f"Prime pairs checked: {len(rows)}")

    for key in [
        "gcd_A_n",
        "gcd_A_n_minus_1",
        "gcd_A_n_plus_1",
        "gcd_A_n2_minus_1",
    ]:
        counts = Counter(row[key] for row in rows)
        print(f"\n{key}")
        print(counts.most_common(15))

    write_csv("gcd_structure.csv", rows)


# ============================================================================
# 16. SIMPLE n-ONLY RELATIONS
# ============================================================================

def simple_relation_test():
    print("\n" + "=" * 80)
    print("16. SIMPLE n-ONLY RELATION SEARCH")
    print("=" * 80)

    ps = primes_upto(PRIME_LIMIT)

    candidates = {
        "A-n^2": lambda n,u,A: A-n*n,
        "A-n(n-1)": lambda n,u,A: A-n*(n-1),
        "A-(n^2-n+1)": lambda n,u,A: A-(n*n-n+1),
        "A-(n+1)^2": lambda n,u,A: A-(n+1)**2,
        "A-u^2": lambda n,u,A: A-u*u,
        "A-nu": lambda n,u,A: A-n*u,
    }

    rows = []

    for p,q in combinations(ps, 2):
        n = p*q
        u = p+q
        A = F(p)*F(q)

        row = {
            "p": p,
            "q": q,
            "n": n,
            "u": u,
            "A": A,
        }

        for name, func in candidates.items():
            residual = func(n,u,A)
            row[name] = residual

        rows.append(row)

    for name in candidates:
        zeros = sum(row[name] == 0 for row in rows)
        print(
            f"{name:25}: "
            f"{zeros}/{len(rows)} exact zeros"
        )

    write_csv("simple_n_relations.csv", rows)


# ============================================================================
# 17. LOW-DEGREE POLYNOMIAL RESIDUAL SEARCH
# ============================================================================

def polynomial_residual_search():
    """
    Search a small family of expressions of the form

        c * A^i * n^j * u^k

    and identify exact integer linear relations among monomials.

    This is exploratory rather than a proof engine.
    """

    print("\n" + "=" * 80)
    print("17. LOW-DEGREE POLYNOMIAL RESIDUAL SEARCH")
    print("=" * 80)

    ps = primes_upto(min(PRIME_LIMIT, 100))

    data = []

    for p,q in combinations(ps,2):
        n = p*q
        u = p+q
        A = F(p)*F(q)

        # Known relation should emerge:
        # A - u^2 + (n+1)u - n^2 + n - 1 = 0
        residual = (
            A
            - u*u
            + (n+1)*u
            - n*n
            + n
            - 1
        )

        data.append({
            "p": p,
            "q": q,
            "n": n,
            "u": u,
            "A": A,
            "known_residual": residual,
        })

        assert residual == 0

    print(f"Rows tested: {len(data)}")
    print(
        "Exact relation found/verified:\n"
        "A - u^2 + (n+1)u - n^2 + n - 1 = 0"
    )

    write_csv("polynomial_residuals.csv", data)


# ============================================================================
# 18. MODULAR FACTOR SIGNATURES
# ============================================================================

def modular_signature_test():
    print("\n" + "=" * 80)
    print("18. MODULAR FACTOR-RECOVERY SIGNATURE SEARCH")
    print("=" * 80)

    ps = primes_upto(100)

    moduli = [3,5,7,11,13,17,19,23,29,31]

    summary = []
    detailed = []

    for m in moduli:
        signatures = defaultdict(set)

        for p,q in combinations(ps,2):
            n = p*q
            A = F(p)*F(q)

            # Signature uses only residue-level quantities.
            signature = (
                n % m,
                A % m,
                F(p) % m,
                F(q) % m,
            )

            pair = tuple(sorted((p % m, q % m)))
            signatures[signature].add(pair)

            detailed.append({
                "modulus": m,
                "p": p,
                "q": q,
                "n_mod_m": n % m,
                "A_mod_m": A % m,
                "Fp_mod_m": F(p) % m,
                "Fq_mod_m": F(q) % m,
                "signature": str(signature),
            })

        unique = sum(
            len(v) == 1
            for v in signatures.values()
        )

        summary.append({
            "modulus": m,
            "signature_count": len(signatures),
            "unique_signature_count": unique,
        })

        print(
            f"mod {m:3}: "
            f"signatures={len(signatures):6}, "
            f"unique={unique:6}"
        )

    write_csv("modular_signatures_summary.csv", summary)
    write_csv("modular_signatures_detail.csv", detailed)


# ============================================================================
# 19. DATASET GENERATION
# ============================================================================

def generate_master_dataset():
    print("\n" + "=" * 80)
    print("19. MASTER DATASET")
    print("=" * 80)

    ps = primes_upto(PRIME_LIMIT)

    rows = []

    for p,q in combinations(ps,2):
        n = p*q
        u = p+q
        A = F(p)*F(q)

        for r in [0,1,2,3,5,7]:
            K3 = kappa([p,q,r])

            rows.append({
                "p": p,
                "q": q,
                "n": n,
                "u": u,
                "A": A,
                "r": r,
                "F_r": F(r),
                "K3": K3,
                "one_minus_K3": 1-K3,
                "A_recovered": (
                    (1-K3)//F(r)
                    if F(r) != 0 else None
                ),
                "D_u": discriminant_from_A(A,n),
                "sqrt_D": math.isqrt(
                    discriminant_from_A(A,n)
                ),
            })

    write_csv("kappa_master_dataset.csv", rows)


# ============================================================================
# 20. FINAL REPORT
# ============================================================================

def final_report():
    print("\n" + "=" * 80)
    print("FINAL RESEARCH REPORT")
    print("=" * 80)

    print(
        """
EXACT RESULTS
-------------

F(x) = x^2 - x + 1

K_N = 1 - product F(x_i)

For K3:

    K3 = 1 - F(p)F(q)F(r)

For K4:

    K4 = 1 - F(p)F(q)F(r)F(s)

Therefore:

    K4-K3 = s(1-s)(1-K3)

This is an exact algebraic identity.


TWO-FACTOR REDUCTION
--------------------

Let

    n = pq
    u = p+q
    A = F(p)F(q)

Then

    A = u^2 -(n+1)u+n^2-n+1.

Hence u satisfies

    u^2 -(n+1)u+n^2-n+1-A = 0.

Its discriminant is

    D = 4A - 3(n-1)^2

and, for a genuine factor pair,

    D = (2u-(n+1))^2.

Therefore an oracle supplying A, equivalently K3 when r is
known, permits algebraic recovery of p and q.


CRITICAL COMPLEXITY QUESTION
----------------------------

The experiments establish:

    n + r + K3  -> p,q

and

    n + r+s + K4 -> p,q.

They do NOT establish:

    n -> K3

efficiently.

If K3 requires knowledge of p and q to calculate, then the
factor-recovery identity is conditional and does not by itself
constitute a new factoring algorithm.

The next meaningful direction is therefore to search for an
efficient n-only method of obtaining K3, A, or an equivalent
quantity that encodes p+q.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 80)
    print("KAPPA CONTROLLED-VARIABLE / FACTOR-RECOVERY ANALYSIS")
    print("=" * 80)

    test_product_identity()
    test_transition_identity()
    test_symmetric_reduction()
    test_discriminant()
    test_oracle_recovery()
    exhaustive_oracle_test()
    collision_search()
    modular_root_search()
    prime_modular_roots()
    controlled_r_experiment()
    controlled_s_experiment()
    controlled_rs_recovery()
    known_r_three_factor_test()
    cross_r_test()
    gcd_structure_test()
    simple_relation_test()
    polynomial_residual_search()
    modular_signature_test()
    generate_master_dataset()
    final_report()

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"All CSV files are in: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()

