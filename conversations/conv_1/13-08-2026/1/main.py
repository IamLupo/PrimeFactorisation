#!/usr/bin/env python3
"""
KAPPA / CYCLOTOMIC RESULTANT / DIVISOR-SUM STRUCTURE SEARCH
============================================================

Purpose
-------
This script attacks the next high-potential questions directly.

We work with

    F(x) = x^2 - x + 1 = Phi_6(x)

and, for squarefree

    n = p*q,

define

    A = F(p)F(q).

The already-established identities are

    A = u^2 -(n+1)u + n^2-n+1
    where u = p+q,

and

    A = sigma_3(n) / sigma_1(n)

for squarefree n.

The new questions tested here are:

1. Can A be expressed through divisor-sum quantities in useful ways?
2. What relations arise modulo F(r)?
3. What happens if TWO auxiliary variables satisfy
       F(r)=0, F(s)=0?
4. Can resultants eliminate r,s and leave a relation involving
       n and u?
5. Are there nontrivial low-degree polynomial relations
   beyond identities already implied by A?
6. Does the cyclotomic structure produce useful n-only
   information?
7. Does the root-of-unity representation reveal additional
   exact structure?
8. Can partial residue information about A translate into
   strong restrictions on u?

NO CSV FILES ARE CREATED.
Everything is printed.
"""

from __future__ import annotations

import math
import random
from itertools import combinations
from collections import defaultdict

import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 200

TEST_PRIME_PAIRS = 120

TARGET_BITS = [20, 30, 40, 50]

AUX = [2, 3, 5, 7, 11, 13, 17, 19, 23]

MODULI = [
    3, 5, 7, 13, 19, 31, 37, 43, 61, 67
]

RANDOM_SEED = 20260813


# =============================================================================
# SYMBOLS
# =============================================================================

p, q = sp.symbols("p q", integer=True)
r, s, t = sp.symbols("r s t", integer=True)

n, u = sp.symbols("n u", integer=True)

x, y, z = sp.symbols("x y z")

A = sp.symbols("A")

e1, e2, e3 = sp.symbols("e1 e2 e3")

omega = sp.symbols("omega")


# =============================================================================
# BASIC FUNCTIONS
# =============================================================================

def F(x):
    return x**2 - x + 1


def F_int(x: int) -> int:
    return x * x - x + 1


def is_prime(n_: int) -> bool:
    if n_ < 2:
        return False

    if n_ % 2 == 0:
        return n_ == 2

    d = 3

    while d * d <= n_:
        if n_ % d == 0:
            return False
        d += 2

    return True


def primes_up_to(limit: int) -> list[int]:
    return [k for k in range(2, limit + 1) if is_prime(k)]


PRIMES = primes_up_to(PRIME_LIMIT)


def A_pq(pv, qv):
    return F(pv) * F(qv)


def A_nu(nv, uv):
    return (
        uv**2
        - (nv + 1) * uv
        + nv**2
        - nv
        + 1
    )


def K_from_A(A_value, rv):
    return 1 - A_value * F(rv)


def line(ch="-", width=78):
    print(ch * width)


# =============================================================================
# DIVISOR SUMS
# =============================================================================

def sigma_k_from_factorization(factors, k):
    """
    factors = [(p,a), ...]
    """
    result = 1

    for prime, exponent in factors:
        numerator = prime ** (k * (exponent + 1)) - 1
        denominator = prime**k - 1
        result *= numerator // denominator

    return result


def sigma_k_squarefree_primes(primes, k):
    result = 1

    for pv in primes:
        result *= 1 + pv**k

    return result


# =============================================================================
# SECTION 1
# DIVISOR-SUM STRUCTURE
# =============================================================================

def section_divisor_sum_structure():
    print()
    print("=" * 78)
    print("1. DIVISOR-SUM / CYCLOTOMIC STRUCTURE")
    print("=" * 78)

    e1_local = p + q
    e2_local = p * q

    sigma1 = (p + 1) * (q + 1)
    sigma3 = (p**3 + 1) * (q**3 + 1)

    A_expr = F(p) * F(q)

    print()
    print("F(x) =")
    print(sp.expand(F(x)))

    print()
    print("sigma_1(pq) =")
    print(sp.expand(sigma1))

    print()
    print("sigma_3(pq) =")
    print(sp.expand(sigma3))

    print()
    print("sigma_3 / sigma_1 =")
    print(sp.factor(A_expr))

    identity = sp.factor(
        sigma3 - sigma1 * A_expr
    )

    print()
    print("Check sigma_3 - sigma_1*A:")
    print(identity)

    print()
    print("PASS" if identity == 0 else "FAIL")

    print()
    print("For n=pq and u=p+q:")
    print()
    print("sigma_1(n) = n + u + 1")
    print("A = sigma_3(n)/sigma_1(n)")

    sigma1_nu = n + u + 1

    sigma3_nu = sp.expand(
        (p**3 + 1) * (q**3 + 1)
    )

    sigma3_nu_reduced = sp.expand(
        u**3
        - 3*n*u
        + 3*u
        + 1
        + n**3
        + 3*n**2
        + 3*n
    )

    print()
    print("sigma_3 in terms of n,u:")
    print(
        "sigma_3 = "
        "n^3 + 3n^2 + 3n + "
        "u^3 - 3nu + 3u + 1"
    )

    direct = sp.expand(
        sigma3_nu.subs({
            p**2: sp.Symbol("p2")
        })
    )

    # Direct symmetric reduction.
    sym_sigma3 = sp.symmetrize(
        sigma3_nu,
        [p, q],
        formal=True
    )

    print()
    print("Symmetric sigma_3 result:")
    print(sym_sigma3)


# =============================================================================
# SECTION 2
# A - F(n) FACTORIZATION
# =============================================================================

def section_A_minus_Fn():
    print()
    print("=" * 78)
    print("2. A - F(n): HIDDEN FACTORIZATION STRUCTURE")
    print("=" * 78)

    A_pq_expr = sp.expand(F(p) * F(q))
    Fn_expr = F(p * q)

    difference = sp.factor(
        A_pq_expr - Fn_expr
    )

    print()
    print("A = F(p)F(q)")
    print("F(n) = F(pq)")
    print()
    print("A - F(n) =")
    print(difference)

    print()
    print("Expected:")
    print("(p+q-1)(p+q-pq)")

    expected = sp.expand(
        (p + q - 1) * (p + q - p*q)
    )

    print("Residual:")
    print(
        sp.expand(difference - expected)
    )

    print()
    print("In n,u:")
    print(
        sp.factor(
            A_nu(n, u) - F(n)
        )
    )

    print()
    print(
        "Therefore:"
    )
    print()
    print(
        "A - F(n) = (u-1)(u-n)"
    )

    print()
    print(
        "This is one of the main objects for the new search."
    )


# =============================================================================
# SECTION 3
# DISCRIMINANT AS A SQUARE
# =============================================================================

def section_discriminant():
    print()
    print("=" * 78)
    print("3. DISCRIMINANT / SQUARE STRUCTURE")
    print("=" * 78)

    A_nu_expr = A_nu(n, u)

    polynomial = sp.expand(
        u**2
        - (n + 1) * u
        + n**2
        - n
        + 1
        - A
    )

    D = sp.factor(
        sp.discriminant(
            polynomial,
            u
        )
    )

    print()
    print("Polynomial:")
    print(polynomial)

    print()
    print("Discriminant:")
    print(D)

    print()
    print("Expected:")
    print(
        "4A - 3(n-1)^2"
    )

    print()
    print(
        "Residual:"
    )
    print(
        sp.expand(
            D - (
                4*A
                - 3*(n - 1)**2
            )
        )
    )

    print()
    print(
        "For the genuine factor pair:"
    )
    print()
    print(
        "D = (2u-(n+1))^2"
    )


# =============================================================================
# SECTION 4
# ONE AUXILIARY CYCLOTOMIC QUOTIENT
# =============================================================================

def section_one_auxiliary():
    print()
    print("=" * 78)
    print("4. ONE AUXILIARY VARIABLE: QUOTIENT BY F(r)")
    print("=" * 78)

    relation_r = r**2 - r + 1

    print()
    print("Cyclotomic relation:")
    print(
        "F(r) = r^2-r+1 = 0"
    )

    print()
    print(
        "Reduce low-degree powers of r modulo F(r)."
    )

    for power in range(0, 9):
        rem = sp.rem(
            r**power,
            relation_r,
            r
        )

        print(
            f"r^{power:<2d} -> {sp.expand(rem)}"
        )

    print()
    print(
        "Now reduce generic expressions involving n and u."
    )

    candidates = [
        u - r,
        u + r,
        u**2 - u*r + r**2,
        F(u),
        F(n),
        A_nu(n, u) - F(r),
        A_nu(n, u) - F(r)*n,
        A_nu(n, u) - F(r)*(n+1),
    ]

    for expr in candidates:
        reduced = sp.rem(
            sp.expand(expr),
            relation_r,
            r
        )

        print()
        print("Expression:")
        print(expr)
        print("Reduction modulo F(r):")
        print(sp.expand(reduced))


# =============================================================================
# SECTION 5
# TWO AUXILIARY VARIABLES
# =============================================================================

def section_two_auxiliary():
    print()
    print("=" * 78)
    print("5. TWO AUXILIARY VARIABLES / RESULTANT SEARCH")
    print("=" * 78)

    Fr = F(r)
    Fs = F(s)

    print()
    print("Relations:")
    print("Fr =", Fr)
    print("Fs =", Fs)

    # Symmetric combinations of r,s.
    v1 = sp.expand(r + s)
    v2 = sp.expand(r*s)

    print()
    print("Elementary symmetric auxiliary variables:")
    print("a = r+s")
    print("b = rs")

    # Derive relations when both r,s satisfy F=0.
    # Since r^2-r+1=0 and s^2-s+1=0.
    print()
    print("Relations among a=r+s and b=rs:")

    rel1 = sp.symmetrize(
        Fr + Fs,
        [r, s],
        formal=True
    )

    rel2 = sp.symmetrize(
        Fr * Fs,
        [r, s],
        formal=True
    )

    print("Fr+Fs symmetric:")
    print(rel1)

    print("Fr*Fs symmetric:")
    print(rel2)

    # Direct substitution using:
    # r^2 = r-1
    # s^2 = s-1
    print()
    print("Reduction table for monomials r^i s^j:")

    for i in range(0, 5):
        row = []

        for j in range(0, 5):
            expr = r**i * s**j

            expr1 = sp.rem(
                expr,
                Fr,
                r
            )

            expr2 = sp.rem(
                expr1,
                Fs,
                s
            )

            row.append(
                sp.expand(expr2)
            )

        print(f"i={i}: {row}")

    # Resultant search.
    print()
    print("RESULTANT SEARCH")
    line()

    G_candidates = [
        u - r,
        u + r,
        u**2 - (n+1)*u + n**2-n+1,
        F(u) - F(r),
        F(n) - F(r),
        A_nu(n, u) - F(r)*F(s),
        A_nu(n, u) - F(r) - F(s),
    ]

    for G in G_candidates:

        Rr = sp.resultant(
            G,
            Fr,
            r
        )

        if r not in Rr.free_symbols:
            print()
            print("G =")
            print(G)
            print()
            print("Resultant eliminating r:")
            print(
                sp.factor(Rr)
            )

    print()
    print(
        "Now eliminate r and s from simple combinations."
    )

    H_candidates = [
        r + s - 1,
        r*s - 1,
        F(r) + F(s),
        F(r) - F(s),
        F(r)*F(s),
        (r-s)**2,
    ]

    for H in H_candidates:

        R1 = sp.resultant(
            H,
            Fr,
            r
        )

        R2 = sp.resultant(
            R1,
            Fs,
            s
        )

        print()
        print("H =", H)
        print("double resultant =")
        print(sp.factor(R2))


# =============================================================================
# SECTION 6
# EXACT AUXILIARY ROOTS AND CUBIC ROOTS OF UNITY
# =============================================================================

def section_root_of_unity():
    print()
    print("=" * 78)
    print("6. ROOT-OF-UNITY REPRESENTATION")
    print("=" * 78)

    relation = omega**2 + omega + 1

    print()
    print("omega^2 + omega + 1 = 0")

    factor_identity = sp.expand(
        (1 + omega*x)
        * (1 + omega**2*x)
    )

    print()
    print("(1+omega*x)(1+omega^2*x) =")
    print(factor_identity)

    reduced = sp.rem(
        factor_identity,
        relation,
        omega
    )

    print()
    print("Modulo omega^2+omega+1:")
    print(sp.expand(reduced))

    print()
    print(
        "Therefore:"
    )
    print()
    print(
        "F(x) = (1+omega*x)(1+omega^2*x)"
    )

    # Generic polynomial E(t) for 2-body and 3-body.
    tvar = sp.symbols("t")

    E2 = (
        1
        + e1*tvar
        + e2*tvar**2
    )

    print()
    print("For two factors:")
    print(
        "E(t) =", E2
    )

    product_root = sp.expand(
        E2.subs(tvar, omega)
        * E2.subs(tvar, omega**2)
    )

    reduced_root = sp.rem(
        product_root,
        relation,
        omega
    )

    print()
    print(
        "E(omega)E(omega^2) reduced:"
    )
    print(
        sp.expand(reduced_root)
    )

    print()
    print(
        "Expected A in e1,e2:"
    )

    A_e = sp.expand(
        e1**2
        - e1*e2
        - e1
        + e2**2
        - e2
        + 1
    )

    print(A_e)

    print()
    print("Residual:")
    print(
        sp.expand(
            reduced_root - A_e
        )
    )


# =============================================================================
# SECTION 7
# THREE AUXILIARY VARIABLES
# =============================================================================

def section_three_auxiliary():
    print()
    print("=" * 78)
    print("7. THREE AUXILIARY VARIABLES")
    print("=" * 78)

    Ft = F(t)

    print()
    print(
        "All three satisfy:"
    )
    print(
        "F(r)=F(s)=F(t)=0"
    )

    # Symmetric sums of three roots of F(x).
    a, b, c = sp.symbols("a b c")

    print()
    print(
        "Search for possible values of:"
    )
    print(
        "a = r+s+t"
    )
    print(
        "b = rs+rt+st"
    )
    print(
        "c = rst"
    )

    # Since r,s,t are each roots of same quadratic,
    # each lies in a 2-dimensional algebra.
    # Search numerically over a modest set of actual roots modulo primes.
    print()
    print(
        "Modular enumeration:"
    )

    for mod in [3, 7, 13, 19, 31]:

        roots = [
            z
            for z in range(mod)
            if (z*z-z+1) % mod == 0
        ]

        triples = set()

        for rv in roots:
            for sv in roots:
                for tv in roots:

                    triples.add(
                        (
                            (rv + sv + tv) % mod,
                            (rv*sv + rv*tv + sv*tv) % mod,
                            (rv*sv*tv) % mod,
                        )
                    )

        print()
        print(
            f"mod {mod}: roots={roots}"
        )
        print(
            f"  distinct (a,b,c) triples = {len(triples)}"
        )
        print(
            f"  triples = {sorted(triples)}"
        )


# =============================================================================
# SECTION 8
# SEARCH FOR n-u RELATIONS FROM AUXILIARY RESULTANTS
# =============================================================================

def section_resultant_candidate_search():
    print()
    print("=" * 78)
    print("8. LOW-DEGREE RESULTANT CANDIDATE SEARCH")
    print("=" * 78)

    print()
    print(
        "We construct a small family of polynomials G(n,u,r,s)"
    )
    print(
        "and eliminate r,s under F(r)=F(s)=0."
    )

    Fr = F(r)
    Fs = F(s)

    candidates = []

    # Basic linear/quadratic expressions.
    for alpha in range(-2, 3):
        for beta in range(-2, 3):
            for gamma in range(-2, 3):

                if alpha == beta == gamma == 0:
                    continue

                candidates.append(
                    sp.expand(
                        alpha*u
                        + beta*r
                        + gamma*s
                    )
                )

    # Add structurally more interesting expressions.
    candidates += [
        u - (r+s),
        u - r*s,
        u - F(r),
        u - F(s),
        u - (r+s-r*s),

        A_nu(n,u) - F(r)*F(s),
        A_nu(n,u) - F(r),
        A_nu(n,u) - F(s),

        F(u) - F(r)*F(s),

        F(n) - F(r)*F(s),

        A_nu(n,u) - F(n),
    ]

    seen = set()
    hits = []

    print()
    print(
        "Testing candidates..."
    )

    for G in candidates:

        key = str(sp.expand(G))

        if key in seen:
            continue

        seen.add(key)

        # Eliminate r.
        R1 = sp.resultant(
            G,
            Fr,
            r
        )

        # Eliminate s.
        R2 = sp.resultant(
            R1,
            Fs,
            s
        )

        R2 = sp.factor(
            sp.expand(R2)
        )

        # Interesting if result depends on n/u
        # and isn't identically zero or a pure constant.
        free = R2.free_symbols

        if R2 != 0 and (
            n in free or u in free
        ):
            hits.append(
                (G, R2)
            )

    print()
    print(
        f"Candidates tested: {len(seen)}"
    )
    print(
        f"n/u-dependent resultant hits: {len(hits)}"
    )

    for G, R in hits[:40]:

        print()
        print("G =")
        print(G)

        print()
        print("Eliminated relation =")
        print(R)

        print()
        print(
            "Factorized:"
        )
        print(
            sp.factor(R)
        )

    print()
    print(
        "Important:"
    )
    print(
        "A resultant is only genuinely interesting if the resulting"
    )
    print(
        "n,u relation is not merely a disguised restatement of the"
    )
    print(
        "known A(u) identity."
    )


# =============================================================================
# SECTION 9
# NUMERICAL TEST OF ALL CANDIDATE RELATIONS
# =============================================================================

def section_numerical_relations():
    print()
    print("=" * 78)
    print("9. NUMERICAL CROSS-TARGET RESULTANT TEST")
    print("=" * 78)

    random.seed(RANDOM_SEED)

    # Fixed moderate semiprime examples.
    test_pairs = [
        (17, 23),
        (31, 47),
        (41, 43),
        (53, 59),
        (61, 67),
        (71, 73),
        (79, 83),
        (89, 97),
    ]

    print()
    print(
        "Testing the main exact identities across independent p,q."
    )

    failures = 0

    for pv, qv in test_pairs:

        nv = pv * qv
        uv = pv + qv
        Av = A_pq(pv, qv)

        check1 = (
            Av
            - (
                uv**2
                - (nv+1)*uv
                + nv**2
                - nv
                + 1
            )
        )

        check2 = (
            Av
            - F(nv)
            - (
                (uv-1)
                * (uv-nv)
            )
        )

        check3 = (
            4*Av
            - 3*(nv-1)**2
            - (2*uv-(nv+1))**2
        )

        if check1 != 0 or check2 != 0 or check3 != 0:
            failures += 1

        print()
        print(
            f"p={pv:3d} q={qv:3d} "
            f"n={nv:6d} u={uv:4d}"
        )
        print(
            f"  A(u) identity : {check1}"
        )
        print(
            f"  A-F(n) identity: {check2}"
        )
        print(
            f"  discriminant  : {check3}"
        )

    print()
    print(
        f"Failures: {failures}"
    )

    if failures == 0:
        print("PASS")


# =============================================================================
# SECTION 10
# PARTIAL A RESIDUE -> U SURVIVAL
# =============================================================================

def section_partial_A():
    print()
    print("=" * 78)
    print("10. PARTIAL A INFORMATION -> u SURVIVAL")
    print("=" * 78)

    pv, qv = 22612043, 26706517

    nv = pv * qv
    uv = pv + qv
    Av = A_pq(pv, qv)

    print()
    print(f"p={pv}")
    print(f"q={qv}")
    print(f"n={nv}")
    print(f"true u={uv}")
    print(f"A bits={Av.bit_length()}")

    print()
    print(
        "For each modulus m, give ONLY A mod m."
    )
    print(
        "Count how many u residues survive."
    )

    print()
    print(
        "m       A mod m        surviving u residues       true u mod m"
    )
    line()

    for m in MODULI:

        ares = Av % m

        survivors = []

        nm = nv % m

        for ur in range(m):

            candidate_A = (
                ur**2
                - (nm + 1)*ur
                + nm**2
                - nm
                + 1
            ) % m

            if candidate_A == ares:
                survivors.append(ur)

        print(
            f"{m:3d}"
            f"{ares:14d}"
            f"{str(survivors):>28s}"
            f"{uv % m:18d}"
        )


# =============================================================================
# SECTION 11
# SEARCH FOR A SECOND INVARIANT FROM TWO AUXILIARIES
# =============================================================================

def section_second_invariant():
    print()
    print("=" * 78)
    print("11. SEARCH FOR A SECOND AUXILIARY INVARIANT")
    print("=" * 78)

    print()
    print(
        "If all auxiliary K-values are generated by"
    )
    print(
        "    K_r = 1 - A F(r),"
    )
    print(
        "then every polynomial expression in K_r values should"
    )
    print(
        "factor through the single variable A."
    )

    # Symbolic demonstration.
    A_sym = sp.symbols("A_sym")

    kr = 1 - A_sym*F(r)
    ks = 1 - A_sym*F(s)

    expressions = {
        "K_r + K_s":
            kr + ks,

        "K_r - K_s":
            kr - ks,

        "K_r*K_s":
            kr*ks,

        "F(s)K_r - F(r)K_s":
            F(s)*kr - F(r)*ks,

        "(K_r-K_s)^2":
            (kr-ks)**2,
    }

    for name, expr in expressions.items():

        print()
        print(name)
        print("= ")
        print(
            sp.factor(
                sp.expand(expr)
            )
        )

    print()
    print(
        "This makes the structural obstruction explicit:"
    )
    print()
    print(
        "the entire auxiliary K-family has only one hidden scalar A."
    )
    print(
        "A genuinely new invariant must therefore come from some"
    )
    print(
        "structure outside the K_r affine family."
    )


# =============================================================================
# SECTION 12
# FINAL SYNTHESIS
# =============================================================================

def section_final():
    print()
    print("=" * 78)
    print("12. FINAL SYNTHESIS")
    print("=" * 78)

    print()
    print("ESTABLISHED:")
    print()

    print(
        "1. F(x) = x^2-x+1 = Phi_6(x)"
    )

    print(
        "2. A(n) = product_{p|n} F(p)"
    )

    print(
        "3. For squarefree n:"
    )

    print(
        "       A(n) = sigma_3(n)/sigma_1(n)"
    )

    print(
        "4. For n=pq and u=p+q:"
    )

    print(
        "       A = u^2-(n+1)u+n^2-n+1"
    )

    print(
        "5. Also:"
    )

    print(
        "       A-F(n) = (u-1)(u-n)"
    )

    print(
        "6. And:"
    )

    print(
        "       4A-3(n-1)^2 = (2u-(n+1))^2"
    )

    print()
    print("AUXILIARY STRUCTURE:")
    print()

    print(
        "7. K_r = 1-AF(r)"
    )

    print(
        "8. Therefore all auxiliary K-values are affine evaluations"
    )

    print(
        "   of the same hidden scalar A."
    )

    print()
    print("NEW TARGET:")
    print()

    print(
        "9. The next useful object is NOT another combination of K_r."
    )

    print(
        "10. Search for a second invariant arising from the"
    )

    print(
        "    cyclotomic condition F(r)=0,"
    )

    print(
        "    especially after eliminating r,s by resultants."
    )

    print()
    print("HIGH-VALUE OUTCOMES:")
    print()

    print(
        "A. Find a nontrivial G(n,u)=0 mod M(r,s,...)"
    )

    print(
        "B. Find a congruence for u or p/q computable from n"
    )

    print(
        "C. Find an efficiently computable divisor-sum residue"
    )

    print(
        "D. Prove the auxiliary family contributes no second invariant"
    )

    print()
    print(
        "The script deliberately treats A and K as ORACLE quantities"
    )

    print(
        "only in the control sections, so a discovered relation in"
    )

    print(
        "the elimination sections would be genuinely new."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print()
    print("=" * 78)
    print("KAPPA / CYCLOTOMIC / RESULTANT RESEARCH SCRIPT")
    print("=" * 78)

    section_divisor_sum_structure()

    section_A_minus_Fn()

    section_discriminant()

    section_one_auxiliary()

    section_two_auxiliary()

    section_root_of_unity()

    section_three_auxiliary()

    section_resultant_candidate_search()

    section_numerical_relations()

    section_partial_A()

    section_second_invariant()

    section_final()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

