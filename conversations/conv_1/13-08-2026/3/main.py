#!/usr/bin/env python3

"""
==============================================================================
DUAL CYCLOTOMIC INVARIANT EXPERIMENT
==============================================================================

Goal
----

We now study two related cyclotomic factors simultaneously:

    Phi_3(x) = x^2 + x + 1
    Phi_6(x) = x^2 - x + 1

For

    n = p*q
    u = p+q

define

    A3 = Phi_3(p) Phi_3(q)
    A6 = Phi_6(p) Phi_6(q)

The key identities to test are expected to be

    A3 = u^2 + (n+1)u + n^2+n+1

    A6 = u^2 - (n+1)u + n^2-n+1

and therefore

    A3 - A6 = 2(n+1)u + 2n

so that

    u = (A3-A6-2n)/(2(n+1)).

This experiment asks whether the two-invariant system gives genuinely
new leverage over the one-invariant Kappa/Phi_6 system.

We test:

1. Exact symbolic reductions of A3 and A6.
2. Linear recovery of u from A3-A6.
3. Recovery of p,q from the two invariants.
4. Controlled auxiliary primes for BOTH cyclotomic families.
5. Whether auxiliary values remain one-dimensional separately,
   but become useful when two families are combined.
6. Modular recovery of u from partial A3 and A6 residues.
7. Whether n-only arithmetic can constrain A3 and A6 residues.
8. Whether the pair (A3,A6) has collisions for fixed n.
9. Resultant/elimination searches involving the two families.
10. Multi-target numerical testing.

NO CSV FILES.
EVERYTHING IS PRINTED.
"""


from __future__ import annotations

import random
from itertools import combinations
from math import gcd, isqrt

import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 500

AUX_LIMIT = 50

MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47
]

TEST_PRIME_LIMIT = 200

RANDOM_SEED = 20260813


# =============================================================================
# SYMBOLS
# =============================================================================

p, q = sp.symbols("p q")
n, u = sp.symbols("n u")
r, s = sp.symbols("r s")

A3, A6 = sp.symbols("A3 A6")

x = sp.symbols("x")


# =============================================================================
# BASIC FUNCTIONS
# =============================================================================

def phi3(x):
    return x**2 + x + 1


def phi6(x):
    return x**2 - x + 1


def F3_int(x: int) -> int:
    return x * x + x + 1


def F6_int(x: int) -> int:
    return x * x - x + 1


def A3_int(p: int, q: int) -> int:
    return F3_int(p) * F3_int(q)


def A6_int(p: int, q: int) -> int:
    return F6_int(p) * F6_int(q)


def bitlen(x: int) -> int:
    return abs(x).bit_length() if x else 1


def is_prime(n_: int) -> bool:
    if n_ < 2:
        return False

    if n_ == 2:
        return True

    if n_ % 2 == 0:
        return False

    d = 3

    while d * d <= n_:
        if n_ % d == 0:
            return False
        d += 2

    return True


def primes_up_to(limit: int) -> list[int]:
    return [
        x
        for x in range(2, limit + 1)
        if is_prime(x)
    ]


PRIMES = primes_up_to(PRIME_LIMIT)


def exact_sqrt(value: int):
    if value < 0:
        return None

    root = isqrt(value)

    if root * root == value:
        return root

    return None


def section(title: str):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line():
    print("-" * 78)


# =============================================================================
# EXACT TWO-BODY FORMULAS
# =============================================================================

def A3_nu(nv, uv):
    return (
        uv**2
        + (nv + 1) * uv
        + nv**2
        + nv
        + 1
    )


def A6_nu(nv, uv):
    return (
        uv**2
        - (nv + 1) * uv
        + nv**2
        - nv
        + 1
    )


def recover_u_from_A_pair(nv: int, a3: int, a6: int):
    numerator = a3 - a6 - 2 * nv
    denominator = 2 * (nv + 1)

    if numerator % denominator != 0:
        return None

    return numerator // denominator


# =============================================================================
# SECTION 1
# EXACT SYMBOLIC REDUCTION
# =============================================================================

def section_symbolic_reduction():
    section("1. EXACT SYMBOLIC REDUCTION")

    a3_raw = sp.expand(
        phi3(p) * phi3(q)
    )

    a6_raw = sp.expand(
        phi6(p) * phi6(q)
    )

    print()
    print("A3 = Phi_3(p) Phi_3(q)")
    print(a3_raw)

    print()
    print("A6 = Phi_6(p) Phi_6(q)")
    print(a6_raw)

    print()
    print("Expected symmetric variables:")
    print("    n = p*q")
    print("    u = p+q")

    a3_expected = sp.expand(
        u**2
        + (n + 1) * u
        + n**2
        + n
        + 1
    )

    a6_expected = sp.expand(
        u**2
        - (n + 1) * u
        + n**2
        - n
        + 1
    )

    a3_sym = sp.symmetrize(
        a3_raw,
        [p, q],
        formal=True
    )

    a6_sym = sp.symmetrize(
        a6_raw,
        [p, q],
        formal=True
    )

    print()
    print("A3 symmetric:")
    print(a3_sym)

    print()
    print("A6 symmetric:")
    print(a6_sym)

    e1, e2 = sp.symbols("e1 e2")

    a3_e = (
        e1**2
        + e1
        + e2**2
        + e2
        + 1
        + e1 * e2
    )

    a6_e = (
        e1**2
        - e1
        + e2**2
        - e2
        + 1
        - e1 * e2
    )

    print()
    print("A3 in e1,e2:")
    print(a3_e)

    print()
    print("A6 in e1,e2:")
    print(a6_e)

    print()
    print("A3-A6:")
    print(
        sp.factor(
            a3_e - a6_e
        )
    )


# =============================================================================
# SECTION 2
# DIFFERENCE / SUM STRUCTURE
# =============================================================================

def section_sum_difference():
    section("2. SUM / DIFFERENCE STRUCTURE")

    a3 = A3_nu(n, u)
    a6 = A6_nu(n, u)

    difference = sp.factor(
        sp.expand(a3 - a6)
    )

    summation = sp.factor(
        sp.expand(a3 + a6)
    )

    product = sp.factor(
        sp.expand(a3 * a6)
    )

    print()
    print("A3-A6 =")
    print(difference)

    print()
    print("A3+A6 =")
    print(summation)

    print()
    print("A3*A6 =")
    print(product)

    recovered_u = sp.factor(
        (a3 - a6 - 2*n)
        / (2*(n+1))
    )

    print()
    print("u reconstructed from A3-A6:")
    print(recovered_u)

    print()
    print(
        "Expected:"
    )
    print("u")


# =============================================================================
# SECTION 3
# NUMERICAL EXACT RECOVERY
# =============================================================================

def section_exact_recovery():
    section("3. EXACT RECOVERY FROM (A3,A6)")

    examples = [
        (2, 3),
        (5, 11),
        (17, 23),
        (31, 47),
        (101, 103),
        (127, 131),
        (179, 181),
    ]

    failures = 0

    print()
    print(
        "p       q       n          u          A3            A6"
    )
    line()

    for pv, qv in examples:

        nv = pv * qv
        uv = pv + qv

        a3v = A3_int(pv, qv)
        a6v = A6_int(pv, qv)

        recovered_u = recover_u_from_A_pair(
            nv,
            a3v,
            a6v
        )

        print(
            f"{pv:5d}"
            f"{qv:8d}"
            f"{nv:11d}"
            f"{uv:11d}"
            f"{a3v:15d}"
            f"{a6v:15d}"
        )

        if recovered_u != uv:
            failures += 1

    print()
    print(f"Recovery failures = {failures}")

    if failures == 0:
        print("PASS: A3 and A6 recover u exactly.")


# =============================================================================
# SECTION 4
# DIRECT POLYNOMIAL FACTORIZATION
# =============================================================================

def section_factor_recovery():
    section("4. FACTOR RECOVERY FROM THE DUAL INVARIANT PAIR")

    examples = [
        (17, 23),
        (31, 47),
        (101, 103),
        (127, 131),
        (179, 181),
    ]

    for pv, qv in examples:

        nv = pv * qv
        uv = pv + qv

        a3v = A3_int(pv, qv)
        a6v = A6_int(pv, qv)

        recovered_u = recover_u_from_A_pair(
            nv,
            a3v,
            a6v
        )

        print()
        print(
            f"p={pv}, q={qv}, n={nv}"
        )

        print(
            f"A3={a3v}"
        )

        print(
            f"A6={a6v}"
        )

        print(
            f"Recovered u={recovered_u}"
        )

        if recovered_u is None:
            print("FAIL")
            continue

        discriminant = (
            recovered_u**2
            - 4 * nv
        )

        root = exact_sqrt(discriminant)

        print(
            f"factor discriminant={discriminant}"
        )

        if root is not None:

            p1 = (recovered_u + root) // 2
            q1 = (recovered_u - root) // 2

            print(
                f"Recovered factors=({p1},{q1})"
            )


# =============================================================================
# SECTION 5
# AUXILIARY K3/K6 FAMILIES
# =============================================================================

def section_auxiliary_families():
    section("5. DUAL AUXILIARY FAMILIES")

    print()
    print(
        "For hidden A3,A6 define:"
    )

    print(
        "K3(r) = 1 - A3*Phi_3(r)"
    )

    print(
        "K6(r) = 1 - A6*Phi_6(r)"
    )

    print()
    print(
        "Representative hidden pair:"
    )

    pv = 22612043
    qv = 26706517

    a3v = A3_int(pv, qv)
    a6v = A6_int(pv, qv)

    print(f"p={pv}")
    print(f"q={qv}")
    print(f"A3={a3v}")
    print(f"A6={a6v}")

    print()
    print(
        "r       Phi3(r)       K3(r)"
        "                 Phi6(r)       K6(r)"
    )
    line()

    for rv in PRIMES:

        if rv > AUX_LIMIT:
            break

        f3 = F3_int(rv)
        f6 = F6_int(rv)

        k3 = 1 - a3v * f3
        k6 = 1 - a6v * f6

        print(
            f"{rv:3d}"
            f"{f3:14d}"
            f"{k3:24d}"
            f"{f6:14d}"
            f"{k6:24d}"
        )


# =============================================================================
# SECTION 6
# PAIRWISE ELIMINATION INSIDE EACH FAMILY
# =============================================================================

def section_family_elimination():
    section("6. ELIMINATION WITHIN EACH CYCLOTOMIC FAMILY")

    print()
    print("For Phi_3 family:")
    print(
        "K3(r)=1-A3*Phi3(r)"
    )

    print(
        "F3(s)*K3(r)-F3(r)*K3(s)"
    )

    expression3 = sp.factor(
        phi3(s) * (
            1 - A3 * phi3(r)
        )
        - phi3(r) * (
            1 - A3 * phi3(s)
        )
    )

    print(
        " ="
    )
    print(expression3)

    print()
    print("For Phi_6 family:")

    expression6 = sp.factor(
        phi6(s) * (
            1 - A6 * phi6(r)
        )
        - phi6(r) * (
            1 - A6 * phi6(s)
        )
    )

    print(
        "F6(s)*K6(r)-F6(r)*K6(s) ="
    )
    print(expression6)

    print()
    print(
        "These should eliminate A3/A6 completely."
    )


# =============================================================================
# SECTION 7
# CROSS-FAMILY ELIMINATION
# =============================================================================

def section_cross_family():
    section("7. CROSS-FAMILY COMBINATIONS")

    print()
    print(
        "Now test whether combining K3 and K6 produces"
    )
    print(
        "information about u rather than merely A3 or A6."
    )

    # Generic K's.
    K3r = 1 - A3 * phi3(r)
    K6r = 1 - A6 * phi6(r)

    expressions = {
        "K3(r)-K6(r)":
            K3r - K6r,

        "Phi6(r) K3(r) - Phi3(r) K6(r)":
            phi6(r) * K3r
            - phi3(r) * K6r,

        "Phi3(r) K6(r) - Phi6(r) K3(r)":
            phi3(r) * K6r
            - phi6(r) * K3r,

        "K3(r)*Phi6(s)-K6(s)*Phi3(r)":
            K3r * phi6(s)
            - K6r.subs(r, s) * phi3(r),
    }

    for name, expr in expressions.items():

        print()
        print(name)
        print("=")
        print(
            sp.factor(
                sp.expand(expr)
            )
        )


# =============================================================================
# SECTION 8
# SUM / DIFFERENCE FROM AUXILIARY K VALUES
# =============================================================================

def section_recover_A_pair_from_auxiliary():
    section("8. RECOVER (A3,A6) FROM CONTROLLED AUXILIARY VALUES")

    pv = 17
    qv = 23

    a3v = A3_int(pv, qv)
    a6v = A6_int(pv, qv)

    print()
    print(
        f"Hidden p={pv}, q={qv}"
    )

    print(
        f"True A3={a3v}"
    )

    print(
        f"True A6={a6v}"
    )

    test_r = [2, 3, 5, 7, 11]

    failures = 0

    for rv in test_r:

        f3 = F3_int(rv)
        f6 = F6_int(rv)

        k3 = 1 - a3v * f3
        k6 = 1 - a6v * f6

        recovered_a3 = (1 - k3) // f3
        recovered_a6 = (1 - k6) // f6

        nvalue = pv * qv

        recovered_u = recover_u_from_A_pair(
            nvalue,
            recovered_a3,
            recovered_a6
        )

        print()
        print(
            f"r={rv}"
        )

        print(
            f"  recovered A3={recovered_a3}"
        )

        print(
            f"  recovered A6={recovered_a6}"
        )

        print(
            f"  recovered u={recovered_u}"
        )

        if (
            recovered_a3 != a3v
            or recovered_a6 != a6v
            or recovered_u != pv + qv
        ):
            failures += 1

    print()
    print(
        f"Failures={failures}"
    )

    if failures == 0:
        print(
            "PASS: one controlled r gives both invariant channels."
        )


# =============================================================================
# SECTION 9
# MODULAR DUAL-INVARIANT RECOVERY
# =============================================================================

def section_modular_dual_recovery():
    section("9. MODULAR DUAL-INVARIANT RECOVERY")

    pv = 22612043
    qv = 26706517

    nv = pv * qv
    uv = pv + qv

    a3v = A3_int(pv, qv)
    a6v = A6_int(pv, qv)

    print()
    print(
        f"n={nv}"
    )

    print(
        f"true u={uv}"
    )

    print(
        f"A3 bits={bitlen(a3v)}"
    )

    print(
        f"A6 bits={bitlen(a6v)}"
    )

    print()
    print(
        "m   A3 mod m   A6 mod m   possible u residues"
    )
    line()

    for m in MODULI:

        a3mod = a3v % m
        a6mod = a6v % m

        survivors = []

        for ur in range(m):

            c3 = (
                ur**2
                + (nv % m + 1) * ur
                + (nv % m)**2
                + nv % m
                + 1
            ) % m

            c6 = (
                ur**2
                - (nv % m + 1) * ur
                + (nv % m)**2
                - nv % m
                + 1
            ) % m

            if c3 == a3mod and c6 == a6mod:
                survivors.append(ur)

        print(
            f"{m:3d}"
            f"{a3mod:12d}"
            f"{a6mod:12d}"
            f"{str(survivors):>24s}"
        )


# =============================================================================
# SECTION 10
# N-ONLY DUAL-INVARIANT RESIDUE SEARCH
# =============================================================================

def section_n_only_dual_search():
    section("10. n-ONLY DUAL-INVARIANT SEARCH")

    print()
    print(
        "Critical experiment:"
    )

    print(
        "For each modulus m, compute the possible pair"
    )

    print(
        "    (A3(u), A6(u)) mod m"
    )

    print(
        "as u ranges over all residues."
    )

    print(
        "If this pair collapses to one value, then both"
    )

    print(
        "invariants would be n-only modulo m."
    )

    nv = 22612043 * 26706517

    print()
    print(
        "m      distinct (A3,A6) pairs     n-only?"
    )
    line()

    hits = []

    for m in MODULI:

        pairs = set()

        nmod = nv % m

        for ur in range(m):

            a3 = (
                ur**2
                + (nmod + 1) * ur
                + nmod**2
                + nmod
                + 1
            ) % m

            a6 = (
                ur**2
                - (nmod + 1) * ur
                + nmod**2
                - nmod
                + 1
            ) % m

            pairs.add(
                (a3, a6)
            )

        unique = len(pairs) == 1

        print(
            f"{m:5d}"
            f"{len(pairs):32d}"
            f"{str(unique):>12s}"
        )

        if unique:
            hits.append(
                (m, pairs)
            )

    print()

    if hits:
        print(
            "IMPORTANT DUAL-INVARIANT n-only residues:"
        )

        for m, pairs in hits:
            print(
                m,
                pairs
            )

    else:
        print(
            "No tested modulus made the dual invariant pair"
        )
        print(
            "uniquely determined by n alone."
        )


# =============================================================================
# SECTION 11
# RESULTANT / ELIMINATION SEARCH
# =============================================================================

def section_resultant_search():
    section("11. CROSS-CYCLOTOMIC RESULTANT SEARCH")

    print()
    print(
        "Search for identities involving the two cyclotomic families"
    )

    print(
        "that eliminate p,q and possibly leave useful relations in n,u."
    )

    e1, e2 = sp.symbols("e1 e2")

    a3e = (
        e1**2
        + e1 * e2
        + e1
        + e2**2
        + e2
        + 1
    )

    a6e = (
        e1**2
        - e1 * e2
        - e1
        + e2**2
        - e2
        + 1
    )

    print()
    print(
        "A3(e1,e2) ="
    )
    print(a3e)

    print()
    print(
        "A6(e1,e2) ="
    )
    print(a6e)

    print()
    print(
        "A3-A6 ="
    )
    print(
        sp.factor(
            a3e - a6e
        )
    )

    print()
    print(
        "Set e1=u, e2=n:"
    )

    print(
        "A3-A6 ="
    )
    print(
        sp.factor(
            (a3e-a6e).subs({
                e1: u,
                e2: n
            })
        )
    )

    print()
    print(
        "Search resultant eliminating u from the difference equation."
    )

    difference_eq = sp.expand(
        A3 - A6
        - 2*(n+1)*u
        - 2*n
    )

    result = sp.resultant(
        difference_eq,
        u - (
            A3-A6-2*n
        )/(2*(n+1)),
        u
    )

    print(
        "Result:"
    )
    print(
        sp.factor(result)
    )

    print()
    print(
        "The main question is whether adding another independent"
    )

    print(
        "cyclotomic family creates a second equation rather than"
    )

    print(
        "another reparameterization of the same u."
    )


# =============================================================================
# SECTION 12
# MULTI-TARGET TEST
# =============================================================================

def section_multi_target():
    section("12. MULTI-TARGET DUAL-INVARIANT TEST")

    examples = [
        (2, 3),
        (5, 11),
        (17, 23),
        (31, 47),
        (41, 43),
        (53, 59),
        (61, 67),
        (71, 73),
        (79, 83),
        (89, 97),
        (101, 103),
        (127, 131),
        (137, 139),
        (149, 151),
        (157, 163),
        (167, 173),
        (179, 181),
        (191, 193),
    ]

    failures = 0

    print()
    print(
        "p       q       n             u          recovered"
    )
    line()

    for pv, qv in examples:

        nv = pv * qv
        uv = pv + qv

        a3v = A3_int(pv, qv)
        a6v = A6_int(pv, qv)

        recovered = recover_u_from_A_pair(
            nv,
            a3v,
            a6v
        )

        ok = recovered == uv

        print(
            f"{pv:5d}"
            f"{qv:8d}"
            f"{nv:14d}"
            f"{uv:12d}"
            f"{str(ok):>12s}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"Failures={failures}"
    )

    if failures == 0:
        print(
            "PASS: dual invariant recovery works on all targets."
        )


# =============================================================================
# SECTION 13
# FINAL REPORT
# =============================================================================

def section_final():
    section("13. FINAL RESEARCH REPORT")

    print()
    print(
        "The new invariant pair is:"
    )

    print()
    print(
        "    A3 = Phi_3(p)Phi_3(q)"
    )

    print(
        "    A6 = Phi_6(p)Phi_6(q)"
    )

    print()
    print(
        "For n=pq and u=p+q:"
    )

    print()
    print(
        "    A3 = u^2 +(n+1)u+n^2+n+1"
    )

    print(
        "    A6 = u^2 -(n+1)u+n^2-n+1"
    )

    print()
    print(
        "Therefore:"
    )

    print()
    print(
        "    A3-A6 = 2(n+1)u+2n"
    )

    print()
    print(
        "and:"
    )

    print()
    print(
        "    u = (A3-A6-2n)/(2(n+1))"
    )

    print()
    print(
        "This is stronger than the single-A6 quadratic recovery."
    )

    print()
    print(
        "The key unresolved problem is:"
    )

    print()
    print(
        "Can BOTH A3 and A6, or enough information about their"
    )

    print(
        "difference, be obtained from n alone?"
    )

    print()
    print(
        "A successful future result would therefore be:"
    )

    print()
    print(
        "    n"
    )
    print(
        "      ↓"
    )
    print(
        "    dual cyclotomic residue structure"
    )
    print(
        "      ↓"
    )
    print(
        "    A3-A6 modulo useful moduli"
    )
    print(
        "      ↓"
    )
    print(
        "    u=p+q"
    )
    print(
        "      ↓"
    )
    print(
        "    p,q"
    )

    print()
    print(
        "NO CSV FILES WERE CREATED."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    random.seed(RANDOM_SEED)

    print()
    print("=" * 78)
    print(
        "DUAL CYCLOTOMIC INVARIANT / FACTOR-RECOVERY EXPERIMENT"
    )
    print("=" * 78)

    section_symbolic_reduction()

    section_sum_difference()

    section_exact_recovery()

    section_factor_recovery()

    section_auxiliary_families()

    section_family_elimination()

    section_cross_family()

    section_recover_A_pair_from_auxiliary()

    section_modular_dual_recovery()

    section_n_only_dual_search()

    section_resultant_search()

    section_multi_target()

    section_final()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

