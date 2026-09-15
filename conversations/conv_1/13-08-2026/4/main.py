#!/usr/bin/env python3

"""
==============================================================================
DUAL CYCLOTOMIC DIFFERENCE / N-ONLY RESIDUE SEARCH
==============================================================================

Main target
-----------

Let

    Phi_3(x) = x^2 + x + 1
    Phi_6(x) = x^2 - x + 1

For

    n = p*q
    u = p+q

define

    A3 = Phi_3(p) Phi_3(q)
    A6 = Phi_6(p) Phi_6(q)

Then the exact identities are

    A3 = u^2 +(n+1)u+n^2+n+1

    A6 = u^2 -(n+1)u+n^2-n+1

and therefore

    Dphi = A3-A6
         = 2(n+1)u.

Hence

    u = Dphi / (2(n+1)).

This experiment focuses ONLY on Dphi.

The central question is:

    Can we obtain useful information about

        Dphi mod m

    from n and auxiliary cyclotomic structure alone?

If yes, then

    Dphi mod m
        ->
    u mod m
        ->
    CRT
        ->
    u
        ->
    p,q.

NO CSV FILES.
EVERYTHING IS PRINTED.

Important separation
--------------------

N-ONLY sections:
    do not use p,q,u,A3,A6,Dphi to choose constraints.

CONTROL sections:
    use the true values only to measure what would happen
    if a Dphi residue were available.

==============================================================================
"""

from __future__ import annotations

from itertools import combinations
from math import gcd, isqrt
import random
import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 500
AUX_LIMIT = 60

MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

RANDOM_SEED = 20260813


# =============================================================================
# BASIC FUNCTIONS
# =============================================================================

def phi3_int(x: int) -> int:
    return x * x + x + 1


def phi6_int(x: int) -> int:
    return x * x - x + 1


def Dphi_true(p: int, q: int) -> int:
    return (
        phi3_int(p) * phi3_int(q)
        - phi6_int(p) * phi6_int(q)
    )


def A3_true(p: int, q: int) -> int:
    return phi3_int(p) * phi3_int(q)


def A6_true(p: int, q: int) -> int:
    return phi6_int(p) * phi6_int(q)


def bitlen(x: int) -> int:
    return abs(x).bit_length() if x else 1


def exact_sqrt(x: int):
    if x < 0:
        return None

    r = isqrt(x)

    return r if r * r == x else None


def is_prime(n: int) -> bool:
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def primes_up_to(limit: int):
    return [
        x for x in range(2, limit + 1)
        if is_prime(x)
    ]


PRIMES = primes_up_to(PRIME_LIMIT)


def section(title: str):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line():
    print("-" * 78)


# =============================================================================
# SYMBOLIC DEFINITIONS
# =============================================================================

p, q = sp.symbols("p q")
n, u = sp.symbols("n u")
r, s = sp.symbols("r s")


# =============================================================================
# SECTION 1
# EXACT SYMBOLIC DIFFERENCE
# =============================================================================

def section_exact_symbolic():
    section("1. EXACT DUAL-CYCLOTOMIC DIFFERENCE")

    P3 = p**2 + p + 1
    Q3 = q**2 + q + 1

    P6 = p**2 - p + 1
    Q6 = q**2 - q + 1

    A3 = sp.expand(P3 * Q3)
    A6 = sp.expand(P6 * Q6)

    D = sp.factor(A3 - A6)

    print()
    print("A3 =")
    print(A3)

    print()
    print("A6 =")
    print(A6)

    print()
    print("Dphi = A3 - A6 =")
    print(D)

    print()
    print("Expected:")
    print("2*(p+q)*(p*q+1)")

    expected = 2 * (p + q) * (p * q + 1)

    print()
    print("Residual:")
    print(
        sp.expand(D - expected)
    )

    e1, e2 = sp.symbols("e1 e2")

    D_e = sp.expand(
        2 * (e1 * e2 + e1)
    )

    print()
    print("Elementary-symmetric form:")
    print(D_e)

    print()
    print(
        "With e1=u and e2=n:"
    )
    print(
        "Dphi = 2(n+1)u"
    )


# =============================================================================
# SECTION 2
# EXACT NUMERICAL CHECKS
# =============================================================================

def section_exact_numeric():
    section("2. NUMERICAL DPHI IDENTITY")

    examples = [
        (2, 3),
        (5, 11),
        (17, 23),
        (31, 47),
        (101, 103),
        (127, 131),
        (179, 181),
        (191, 193),
    ]

    failures = 0

    print()
    print(
        "p      q       n          u          Dphi"
    )
    line()

    for pv, qv in examples:

        nv = pv * qv
        uv = pv + qv

        d = Dphi_true(pv, qv)

        expected = 2 * (nv + 1) * uv

        ok = d == expected

        print(
            f"{pv:5d}"
            f"{qv:8d}"
            f"{nv:11d}"
            f"{uv:11d}"
            f"{d:18d}"
            f"   {ok}"
        )

        if not ok:
            failures += 1

    print()
    print(f"Failures = {failures}")

    if failures == 0:
        print("PASS")


# =============================================================================
# SECTION 3
# RECOVER u EXACTLY FROM DPHI
# =============================================================================

def section_exact_recovery():
    section("3. EXACT u RECOVERY FROM DPHI")

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
        "p       q       n          true u      recovered u"
    )
    line()

    for pv, qv in examples:

        nv = pv * qv
        uv = pv + qv

        d = Dphi_true(pv, qv)

        denominator = 2 * (nv + 1)

        if d % denominator != 0:
            recovered = None
        else:
            recovered = d // denominator

        print(
            f"{pv:6d}"
            f"{qv:8d}"
            f"{nv:12d}"
            f"{uv:14d}"
            f"{str(recovered):>16s}"
        )

        if recovered != uv:
            failures += 1

    print()
    print(f"Failures = {failures}")

    if failures == 0:
        print(
            "PASS: Dphi recovers u linearly."
        )


# =============================================================================
# SECTION 4
# DPHI MOD m -> u MOD m
# =============================================================================

def section_modular_recovery():
    section("4. DPHI MOD m -> u MOD m")

    pv = 22612043
    qv = 26706517

    nv = pv * qv
    uv = pv + qv

    d = Dphi_true(pv, qv)

    print()
    print(f"n = {nv}")
    print(f"u = {uv}")
    print(f"Dphi = {d}")
    print()

    print(
        "m     Dphi mod m   gcd(2(n+1),m)   possible u residues"
    )
    line()

    for m in MODULI:

        dmod = d % m

        coeff = (2 * (nv + 1)) % m

        residues = []

        for ur in range(m):

            if (
                coeff * ur
            ) % m == dmod:
                residues.append(ur)

        print(
            f"{m:3d}"
            f"{dmod:15d}"
            f"{gcd(2*(nv+1),m):18d}"
            f"{str(residues):>30s}"
        )


# =============================================================================
# SECTION 5
# N-ONLY SEARCH:
# CAN n DETERMINE DPHI MOD m?
# =============================================================================

def section_n_only_dphi():
    section("5. N-ONLY DPHI RESIDUE SEARCH")

    pv = 22612043
    qv = 26706517

    nv = pv * qv

    print()
    print(
        "This section does NOT use the true p,q,u."
    )

    print(
        "For each modulus m, compute"
    )

    print(
        "    Dphi(u) = 2(n+1)u mod m"
    )

    print(
        "for every u residue."
    )

    print(
        "If only one Dphi residue occurs, then Dphi mod m"
    )

    print(
        "would be determined by n alone."
    )

    print()
    print(
        "m     distinct Dphi residues     n-only?"
    )
    line()

    hits = []

    for m in MODULI:

        coeff = (2 * (nv + 1)) % m

        values = {
            (coeff * ur) % m
            for ur in range(m)
        }

        unique = len(values) == 1

        print(
            f"{m:3d}"
            f"{len(values):30d}"
            f"{str(unique):>12s}"
        )

        if unique:
            hits.append(
                (m, values)
            )

    print()

    if hits:
        print(
            "IMPORTANT: n-only Dphi residues found:"
        )

        for item in hits:
            print(item)

    else:
        print(
            "No tested modulus made Dphi mod m unique from n alone."
        )


# =============================================================================
# SECTION 6
# N-ONLY AUXILIARY MODULI
# =============================================================================

def n_only_polynomials(nv: int):
    return {
        "n-1": nv - 1,
        "n+1": nv + 1,
        "n^2-1": nv**2 - 1,
        "n^2+n+1": nv**2 + nv + 1,
        "n^2-n+1": nv**2 - nv + 1,
        "n^3-1": nv**3 - 1,
        "n^3+1": nv**3 + 1,
        "n^3-n": nv**3 - nv,
        "n^4-1": nv**4 - 1,
    }


def section_adaptive_moduli():
    section("6. ADAPTIVE AUXILIARY MODULUS SEARCH")

    nv = 22612043 * 26706517

    print()
    print(
        "Only n and auxiliary r are used."
    )

    print(
        "For m=Phi_3(r) and m=Phi_6(r), inspect gcds with"
    )

    print(
        "natural n-only polynomials."
    )

    polys = n_only_polynomials(nv)

    for family_name, family in [
        ("Phi_3", phi3_int),
        ("Phi_6", phi6_int),
    ]:

        print()
        print(
            f"FAMILY = {family_name}"
        )

        print(
            "r    m=F(r)      bits   strongest gcd      expression"
        )
        line()

        candidates = []

        for rv in PRIMES:

            if rv > AUX_LIMIT:
                break

            m = family(rv)

            best_g = 1
            best_name = None

            for name, value in polys.items():

                g = gcd(abs(value), m)

                if g > best_g:
                    best_g = g
                    best_name = name

            if best_g > 1:

                candidates.append(
                    (
                        best_g,
                        rv,
                        m,
                        best_name
                    )
                )

        candidates.sort(
            reverse=True
        )

        for g, rv, m, name in candidates[:30]:

            print(
                f"{rv:3d}"
                f"{m:12d}"
                f"{bitlen(m):8d}"
                f"{g:18d}"
                f"   {name}"
            )


# =============================================================================
# SECTION 7
# DOES AN AUXILIARY MODULUS PROVIDE DPHI INFORMATION?
# =============================================================================

def section_auxiliary_dphi_constraints():
    section("7. AUXILIARY-CYCLOTOMIC DPHI CONSTRAINT SEARCH")

    nv = 22612043 * 26706517

    print()
    print(
        "We test whether congruences involving n and Phi_j(r)"
    )

    print(
        "can force Dphi modulo a divisor of Phi_j(r)."
    )

    print()
    print(
        "There is no A/K oracle in this section."
    )

    for family_name, family in [
        ("Phi_3", phi3_int),
        ("Phi_6", phi6_int),
    ]:

        print()
        print(
            f"FAMILY = {family_name}"
        )

        for rv in PRIMES:

            if rv > AUX_LIMIT:
                break

            m = family(rv)

            # Only manageable moduli for exhaustive residue search.
            if m > 5000:
                continue

            coeff = (2 * (nv + 1)) % m

            values = {
                (coeff * ur) % m
                for ur in range(m)
            }

            # Number of possible Dphi residues.
            print(
                f"r={rv:2d} "
                f"m={m:5d} "
                f"distinct Dphi residues={len(values):5d}"
            )


# =============================================================================
# SECTION 8
# COMBINE TWO AUXILIARY FAMILIES
# =============================================================================

def section_dual_auxiliary_search():
    section("8. DUAL AUXILIARY FAMILY SEARCH")

    nv = 22612043 * 26706517

    print()
    print(
        "Use both"
    )

    print(
        "    Phi_3(r)"
    )

    print(
        "and"
    )

    print(
        "    Phi_6(s)"
    )

    print(
        "to construct n-only moduli."
    )

    print()
    print(
        "Search for shared factors:"
    )

    print(
        "gcd(Phi_3(r), Phi_6(s))."
    )

    print()
    print(
        "r   s     Phi3(r)    Phi6(s)    gcd"
    )
    line()

    hits = []

    small_primes = [
        z
        for z in PRIMES
        if z <= AUX_LIMIT
    ]

    for rv in small_primes:

        for sv in small_primes:

            g = gcd(
                phi3_int(rv),
                phi6_int(sv)
            )

            if g > 1:

                hits.append(
                    (
                        g,
                        rv,
                        sv,
                        phi3_int(rv),
                        phi6_int(sv)
                    )
                )

    hits.sort(
        reverse=True
    )

    for g, rv, sv, m3, m6 in hits[:50]:

        print(
            f"{rv:3d}"
            f"{sv:4d}"
            f"{m3:12d}"
            f"{m6:12d}"
            f"{g:8d}"
        )

    if not hits:
        print(
            "No nontrivial shared factors found."
        )

    print()
    print(
        "Now search gcd with n-only expressions."
    )

    polys = n_only_polynomials(nv)

    for g, rv, sv, m3, m6 in hits[:20]:

        modulus = m3 * m6 // gcd(m3, m6)

        best = []

        for name, value in polys.items():

            gg = gcd(
                abs(value),
                modulus
            )

            if gg > 1:
                best.append(
                    (name, gg)
                )

        if best:

            print()
            print(
                f"r={rv}, s={sv}, combined modulus={modulus}"
            )

            for name, gg in best:
                print(
                    f"  gcd(modulus,{name})={gg}"
                )


# =============================================================================
# SECTION 9
# CRT ACCUMULATION OF u WHEN DPHI RESIDUES ARE KNOWN
# =============================================================================

def crt_pair(a1, m1, a2, m2):
    g = gcd(m1, m2)

    if (a2 - a1) % g != 0:
        return None

    m1g = m1 // g
    m2g = m2 // g

    if m2g == 1:
        return a1 % (m1 * m2g), m1 * m2g

    t = (
        ((a2 - a1) // g)
        * pow(m1g, -1, m2g)
    ) % m2g

    modulus = m1 * m2g

    residue = (
        a1 + m1 * t
    ) % modulus

    return residue, modulus


def section_crt_control():
    section("9. CONTROL: CRT GROWTH FROM DPHI RESIDUES")

    pv = 22612043
    qv = 26706517

    nv = pv * qv
    uv = pv + qv
    d = Dphi_true(pv, qv)

    print()
    print(
        f"n={nv}"
    )

    print(
        f"true u={uv}"
    )

    residue = 0
    modulus = 1

    print()
    print(
        "m       Dphi mod m       u residues       combined modulus bits"
    )
    line()

    for m in MODULI:

        target = d % m

        candidates = []

        for ur in range(m):

            if (
                2 * (nv + 1) * ur
            ) % m == target:

                candidates.append(ur)

        if len(candidates) == 1:

            local_u = candidates[0]

            result = crt_pair(
                residue,
                modulus,
                local_u,
                m
            )

            if result is not None:

                residue, modulus = result

        print(
            f"{m:3d}"
            f"{target:18d}"
            f"{str(candidates):>25s}"
            f"{bitlen(modulus):20d}"
        )

    print()
    print(
        f"Final CRT residue = {residue}"
    )

    print(
        f"Final CRT modulus = {modulus}"
    )

    print(
        f"Final modulus bits = {bitlen(modulus)}"
    )

    print(
        f"True u mod modulus = {uv % modulus}"
    )

    print(
        f"Exact after CRT? {residue == uv % modulus}"
    )


# =============================================================================
# SECTION 10
# FIXED-n COLLISION TEST FOR DPHI
# =============================================================================

def section_fixed_n_collisions():
    section("10. FIXED-n DPHI COLLISION SEARCH")

    print()
    print(
        "For composite n <= 1000, compare different factor pairs."
    )

    collisions = []

    for nv in range(4, 1001):

        pairs = []

        for pv in range(2, int(nv**0.5) + 1):

            if nv % pv != 0:
                continue

            qv = nv // pv

            if pv == qv:
                continue

            d = Dphi_true(
                pv,
                qv
            )

            pairs.append(
                ((pv, qv), d)
            )

        if len(pairs) < 2:
            continue

        grouped = {}

        for pair, value in pairs:

            grouped.setdefault(
                value,
                []
            ).append(pair)

        for value, pair_list in grouped.items():

            if len(pair_list) > 1:

                collisions.append(
                    (
                        nv,
                        value,
                        pair_list
                    )
                )

    print(
        f"DPhi collisions found = {len(collisions)}"
    )

    for nv, value, pairs in collisions[:50]:

        print()
        print(
            f"n={nv}"
        )
        print(
            f"DPhi={value}"
        )
        print(
            f"pairs={pairs}"
        )


# =============================================================================
# SECTION 11
# MULTI-TARGET DPHI / n RELATIONS
# =============================================================================

def section_multitarget():
    section("11. MULTI-TARGET EXACT STRUCTURE TEST")

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

    print()
    print(
        "p      q       n          u          Dphi/(2(n+1))"
    )
    line()

    failures = 0

    for pv, qv in examples:

        nv = pv * qv
        uv = pv + qv

        d = Dphi_true(
            pv,
            qv
        )

        quotient = d // (
            2 * (nv + 1)
        )

        ok = quotient == uv

        print(
            f"{pv:5d}"
            f"{qv:8d}"
            f"{nv:11d}"
            f"{uv:11d}"
            f"{quotient:18d}"
            f"   {ok}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"Failures={failures}"
    )

    if failures == 0:
        print(
            "PASS"
        )


# =============================================================================
# SECTION 12
# FINAL RESEARCH SUMMARY
# =============================================================================

def section_final():
    section("12. FINAL RESEARCH SUMMARY")

    print()
    print(
        "EXACT DUAL-CYCLOTOMIC IDENTITY"
    )

    print()
    print(
        "Dphi = A3-A6"
    )

    print(
        "      = 2(n+1)u"
    )

    print()
    print(
        "Therefore:"
    )

    print(
        "u = Dphi / [2(n+1)]"
    )

    print()
    print(
        "This changes the experimental target."
    )

    print()
    print(
        "We do NOT need A3."
    )

    print(
        "We do NOT need A6."
    )

    print(
        "We only need:"
    )

    print(
        "    Dphi mod useful moduli."
    )

    print()
    print(
        "The key open question is now:"
    )

    print()
    print(
        "Can n-only cyclotomic structure determine"
    )

    print(
        "a nontrivial residue of"
    )

    print(
        "    Dphi = 2(n+1)(p+q)"
    )

    print(
        "without knowing p,q?"
    )

    print()
    print(
        "If YES:"
    )

    print(
        "    Dphi residues"
    )

    print(
        "       -> u residues"
    )

    print(
        "       -> CRT"
    )

    print(
        "       -> u"
    )

    print(
        "       -> p,q"
    )

    print()
    print(
        "If NO:"
    )

    print(
        "the dual-cyclotomic structure still gives an exact"
    )

    print(
        "linear characterization of the hidden factor sum,"
    )

    print(
        "but not yet an n-only factoring mechanism."
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
        "DUAL CYCLOTOMIC DIFFERENCE / N-ONLY SEARCH"
    )
    print("=" * 78)

    section_exact_symbolic()

    section_exact_numeric()

    section_exact_recovery()

    section_modular_recovery()

    section_n_only_dphi()

    section_adaptive_moduli()

    section_auxiliary_dphi_constraints()

    section_dual_auxiliary_search()

    section_crt_control()

    section_fixed_n_collisions()

    section_multitarget()

    section_final()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

