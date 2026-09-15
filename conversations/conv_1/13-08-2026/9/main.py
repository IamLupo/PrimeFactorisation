#!/usr/bin/env python3

"""
==============================================================================
3-ADIC CUBIC-ROOT / GENERALIZED ROOT-OF-UNITY EXPERIMENT
==============================================================================

Purpose
-------

Investigate the strongest signal from the previous experiment:

    gcd(n, Phi_3(a^E))

where

    E = smooth exponent / 3^j

and compare this against:

    gcd(n, a^E - 1).

The experiment also extends the same idea to:

    Phi_3
    Phi_5
    Phi_7

No CSV files.
Everything is printed.

IMPORTANT
---------

This script is exploratory and does not assume the hidden factors during
the actual factor-search sections.

Known p,q are used only in controlled diagnostic sections.
"""

from __future__ import annotations

from math import gcd, lcm
from collections import Counter, defaultdict

import sympy as sp


# =============================================================================
# CONFIGURATION
# =============================================================================

BASES = list(range(2, 21))

SMOOTH_BOUNDS = [
    5,
    7,
    10,
    12,
    15,
    20,
    25,
    30,
]

TARGETS = [
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
    (211, 223),
    (227, 229),
    (233, 239),
    (241, 251),
]

# This is used only for controlled analysis.
REPRESENTATIVE_P = 157
REPRESENTATIVE_Q = 163
REPRESENTATIVE_N = REPRESENTATIVE_P * REPRESENTATIVE_Q


# =============================================================================
# BASIC HELPERS
# =============================================================================

def section(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def line() -> None:
    print("-" * 78)


def bitlen(x: int) -> int:
    return abs(int(x)).bit_length() if x else 1


def nontrivial(g: int, n: int) -> bool:
    return 1 < g < n


def valuation(x: int, p: int) -> int:
    x = abs(int(x))

    if x == 0:
        return 999999

    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def smooth_lcm(B: int) -> int:
    value = 1

    for k in range(2, B + 1):
        value = lcm(value, k)

    return value


def multiplicative_order(a: int, p: int) -> int:
    return int(sp.n_order(a % p, p))


# =============================================================================
# CYCLOTOMIC POLYNOMIALS
# =============================================================================

def cyclotomic_value(k: int, z: int) -> int:
    """
    Evaluate Phi_k(z) without symbolic substitution.

    The small k values used here are handled directly.
    """

    if k == 3:
        return z * z + z + 1

    if k == 5:
        return (
            z**4
            + z**3
            + z**2
            + z
            + 1
        )

    if k == 7:
        return (
            z**6
            + z**5
            + z**4
            + z**3
            + z**2
            + z
            + 1
        )

    raise ValueError(f"Unsupported k={k}")


def cyclotomic_expression(k: int) -> str:

    if k == 3:
        return "x^2 + x + 1"

    if k == 5:
        return "x^4 + x^3 + x^2 + x + 1"

    if k == 7:
        return "x^6 + x^5 + x^4 + x^3 + x^2 + x + 1"

    return str(sp.cyclotomic_poly(k, sp.Symbol("x")))


# =============================================================================
# 1. EXACT STRUCTURE
# =============================================================================

def section_exact_structure():

    section("1. EXACT CYCLOTOMIC STRUCTURE")

    print()
    print("Phi_3(x) =", cyclotomic_expression(3))
    print("Phi_5(x) =", cyclotomic_expression(5))
    print("Phi_7(x) =", cyclotomic_expression(7))

    print()
    print("For prime ell:")
    print()
    print("    Phi_ell(x) = 0 mod p")
    print()
    print("means that x has exact order ell modulo p, provided p != ell.")
    print()
    print("The cubic case is therefore:")
    print()
    print("    Phi_3(x)=0 mod p")
    print("    <=> x has order 3 mod p.")


# =============================================================================
# 2. REPRESENTATIVE TARGET
# =============================================================================

def section_representative():

    section("2. REPRESENTATIVE TARGET")

    p = REPRESENTATIVE_P
    q = REPRESENTATIVE_Q
    n = REPRESENTATIVE_N

    print()
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"n = {n}")

    print()
    print(f"p-1 = {p-1} = {sp.factorint(p-1)}")
    print(f"q-1 = {q-1} = {sp.factorint(q-1)}")

    print()
    print(f"v3(p-1) = {valuation(p-1, 3)}")
    print(f"v3(q-1) = {valuation(q-1, 3)}")

    print()
    print("This target is retained because the previous experiment")
    print("showed strong cubic separation for 157 * 163.")


# =============================================================================
# 3. MULTIPLICATIVE ORDER TABLE
# =============================================================================

def section_order_table():

    section("3. MULTIPLICATIVE ORDER TABLE")

    p = REPRESENTATIVE_P
    q = REPRESENTATIVE_Q
    n = p * q

    print()
    print(
        "a     ord_p(a)   v3(ord_p)    "
        "ord_q(a)   v3(ord_q)"
    )

    line()

    for a in BASES:

        if gcd(a, n) != 1:
            continue

        op = multiplicative_order(a, p)
        oq = multiplicative_order(a, q)

        print(
            f"{a:2d}"
            f"{op:12d}"
            f"{valuation(op,3):13d}"
            f"{oq:12d}"
            f"{valuation(oq,3):13d}"
        )


# =============================================================================
# 4. SMOOTH EXPONENT STRUCTURE
# =============================================================================

def section_smooth_exponents():

    section("4. SMOOTH EXPONENT STRUCTURE")

    print()
    print(
        " B       E=lcm(1,...,B)       bits(E)       v3(E)"
    )

    line()

    for B in SMOOTH_BOUNDS:

        E = smooth_lcm(B)

        print(
            f"{B:2d}"
            f"{E:24d}"
            f"{bitlen(E):14d}"
            f"{valuation(E,3):13d}"
        )


# =============================================================================
# 5. CUBIC LADDER — COMPACT VERSION
# =============================================================================

def section_cubic_ladder():

    section("5. CUBIC 3-ADIC EXPONENT LADDER")

    p = REPRESENTATIVE_P
    q = REPRESENTATIVE_Q
    n = p * q

    print()
    print(
        "Only nontrivial cubic gcds are printed."
    )

    print()
    print(
        " B   j   v3(E)   a    gcd(n,Phi3(a^E))"
    )

    line()

    hit_counter = Counter()

    hits = []

    for B in SMOOTH_BOUNDS:

        L = smooth_lcm(B)

        max_j = valuation(L, 3)

        for j in range(max_j + 1):

            E = L // (3 ** j)
            e3 = valuation(E, 3)

            for a in BASES:

                if gcd(a, n) != 1:
                    continue

                x = pow(a, E, n)

                g = gcd(
                    n,
                    cyclotomic_value(3, x)
                )

                if nontrivial(g, n):

                    hit_counter[(e3, g)] += 1

                    hits.append(
                        (B, j, E, e3, a, g)
                    )

                    print(
                        f"{B:2d}"
                        f"{j:4d}"
                        f"{e3:8d}"
                        f"{a:5d}"
                        f"{g:20d}"
                    )

    print()

    print(
        f"Total nontrivial cubic hits = {len(hits)}"
    )

    print()

    print("Compressed hit structure:")
    print()

    for (e3, g), count in sorted(hit_counter.items()):

        print(
            f"v3(E)={e3}: factor={g}, occurrences={count}"
        )


# =============================================================================
# 6. HIT MECHANISM
# =============================================================================

def section_hit_mechanism():

    section("6. HIT MECHANISM")

    p = REPRESENTATIVE_P
    q = REPRESENTATIVE_Q
    n = p * q

    print()
    print(
        "Representative successful cases:"
    )

    line()

    interesting = [
        (5, 1, 12),    # E=20
        (10, 0, 5),    # E=2520
        (10, 2, 12),   # E=280
        (15, 0, 5),
        (15, 2, 3),
        (30, 0, 2),
        (30, 3, 3),
    ]

    for B, j, a in interesting:

        L = smooth_lcm(B)
        E = L // (3 ** j)

        x = pow(a, E, n)
        g = gcd(n, cyclotomic_value(3, x))

        if not nontrivial(g, n):
            continue

        op = multiplicative_order(a, p)
        oq = multiplicative_order(a, q)

        xp = pow(a, E, p)
        xq = pow(a, E, q)

        print()
        print(
            f"B={B}, j={j}, E={E}, v3(E)={valuation(E,3)}, a={a}"
        )

        print(
            f"ord_p(a)={op}, v3={valuation(op,3)}"
        )

        print(
            f"ord_q(a)={oq}, v3={valuation(oq,3)}"
        )

        print(
            f"a^E mod p = {xp}"
        )

        print(
            f"a^E mod q = {xq}"
        )

        print(
            f"Phi3(a^E) mod p = {cyclotomic_value(3, xp) % p}"
        )

        print(
            f"Phi3(a^E) mod q = {cyclotomic_value(3, xq) % q}"
        )

        print(
            f"gcd = {g}"
        )


# =============================================================================
# 7. MULTI-TARGET SUMMARY
# =============================================================================

def section_multitarget():

    section("7. MULTI-TARGET CUBIC SEARCH")

    print()
    print(
        "For each semiprime we ask whether the compact 3-adic ladder"
    )
    print(
        "finds any nontrivial cubic gcd."
    )

    print()
    print(
        "p       q       v3(p-1)   v3(q-1)"
        "     cubic?       factor"
    )

    line()

    total = 0
    success = 0

    for p, q in TARGETS:

        n = p * q

        vp = valuation(p - 1, 3)
        vq = valuation(q - 1, 3)

        found = None

        for B in SMOOTH_BOUNDS:

            L = smooth_lcm(B)

            max_j = valuation(L, 3)

            for j in range(max_j + 1):

                E = L // (3 ** j)

                for a in BASES:

                    if gcd(a, n) != 1:
                        continue

                    x = pow(a, E, n)

                    g = gcd(
                        n,
                        cyclotomic_value(3, x)
                    )

                    if nontrivial(g, n):

                        found = (B, j, a, g)
                        break

                if found:
                    break

            if found:
                break

        total += 1

        if found:
            success += 1

            B, j, a, g = found

            print(
                f"{p:7d}"
                f"{q:8d}"
                f"{vp:11d}"
                f"{vq:10d}"
                f"{'YES':9s}"
                f"{g:14d}"
                f"   B={B},j={j},a={a}"
            )

        else:

            print(
                f"{p:7d}"
                f"{q:8d}"
                f"{vp:11d}"
                f"{vq:10d}"
                f"{'NO':9s}"
                f"{'-':14s}"
            )

    print()
    print(
        f"Cubic success rate = {success}/{total}"
    )


# =============================================================================
# 8. CUBIC VS ORDINARY ORDER
# =============================================================================

def section_compare_order():

    section("8. CUBIC VS ORDINARY ORDER TEST")

    cubic_total = 0
    order_total = 0
    cubic_only = 0
    order_only = 0

    for p, q in TARGETS:

        n = p * q

        got_cubic = False
        got_order = False

        for B in SMOOTH_BOUNDS:

            L = smooth_lcm(B)

            max_j = valuation(L, 3)

            for j in range(max_j + 1):

                E = L // (3 ** j)

                for a in BASES:

                    if gcd(a, n) != 1:
                        continue

                    x = pow(a, E, n)

                    g_order = gcd(n, x - 1)

                    g_cubic = gcd(
                        n,
                        cyclotomic_value(3, x)
                    )

                    if nontrivial(g_order, n):
                        got_order = True

                    if nontrivial(g_cubic, n):
                        got_cubic = True

        if got_cubic:
            cubic_total += 1

        if got_order:
            order_total += 1

        if got_cubic and not got_order:
            cubic_only += 1

        if got_order and not got_cubic:
            order_only += 1

        print(
            f"p={p:7d}"
            f" q={q:7d}"
            f" cubic={str(got_cubic):5s}"
            f" order={str(got_order):5s}"
        )

    print()

    print(
        f"Cubic successes = {cubic_total}/{len(TARGETS)}"
    )

    print(
        f"Order successes = {order_total}/{len(TARGETS)}"
    )

    print(
        f"Cubic-only = {cubic_only}"
    )

    print(
        f"Order-only = {order_only}"
    )


# =============================================================================
# 9. v3(E) SWEEP
# =============================================================================

def section_v3_sweep():

    section("9. v3(E) SWEEP")

    p = REPRESENTATIVE_P
    q = REPRESENTATIVE_Q
    n = p * q

    print()
    print(
        "e = v3(E)       cubic-success bases/factors"
    )

    line()

    for e in range(0, 7):

        E = (
            3 ** e
            * 2
            * 5
            * 7
            * 11
            * 13
            * 17
            * 19
        )

        successes = []

        for a in BASES:

            if gcd(a, n) != 1:
                continue

            x = pow(a, E, n)

            g = gcd(
                n,
                cyclotomic_value(3, x)
            )

            if nontrivial(g, n):
                successes.append((a, g))

        print(
            f"{e:10d}       {successes}"
        )


# =============================================================================
# 10. GENERALIZED ROOT-OF-UNITY TEST
# =============================================================================

def section_generalized_roots():

    section("10. GENERALIZED ROOT-OF-UNITY TEST")

    print()
    print(
        "The cubic case is extended to Phi_5 and Phi_7."
    )

    print()
    print(
        "For each ell:"
    )

    print(
        "    x = a^E mod n"
    )

    print(
        "    g = gcd(n, Phi_ell(x))"
    )

    print()
    print(
        "Only nontrivial hits are printed."
    )

    for ell in [3, 5, 7]:

        print()
        print(
            f"Phi_{ell}(x) = {cyclotomic_expression(ell)}"
        )

        line()

        total_hits = 0
        examples = []

        for p, q in TARGETS:

            n = p * q

            for B in [5, 10, 15, 20, 25, 30]:

                E = smooth_lcm(B)

                # Also test nearby divisions by ell where possible.
                exponents = [E]

                if E % ell == 0:
                    exponents.append(E // ell)

                for E_test in exponents:

                    for a in BASES:

                        if gcd(a, n) != 1:
                            continue

                        x = pow(a, E_test, n)

                        g = gcd(
                            n,
                            cyclotomic_value(ell, x)
                        )

                        if nontrivial(g, n):

                            total_hits += 1

                            record = (
                                p,
                                q,
                                B,
                                E_test,
                                a,
                                g,
                            )

                            if len(examples) < 20:
                                examples.append(record)

        print(
            f"Total hits = {total_hits}"
        )

        if examples:

            print()
            print("First examples:")

            for (
                p,
                q,
                B,
                E,
                a,
                g,
            ) in examples:

                print(
                    f"p={p:7d}"
                    f" q={q:7d}"
                    f" B={B:2d}"
                    f" E={E}"
                    f" a={a:2d}"
                    f" factor={g}"
                )

        else:
            print("No nontrivial hits.")


# =============================================================================
# 11. ROOT-ORDER CORRELATION
# =============================================================================

def section_order_correlation():

    section("11. ROOT-HIT / ORDER CORRELATION")

    p = REPRESENTATIVE_P
    q = REPRESENTATIVE_Q
    n = p * q

    print()
    print(
        "For the representative target, summarize cubic hits"
    )

    print(
        "according to the pair"
    )

    print(
        "    (v3(ord_p(a)), v3(ord_q(a)))"
    )

    print()

    correlation = defaultdict(int)

    for B in SMOOTH_BOUNDS:

        L = smooth_lcm(B)

        max_j = valuation(L, 3)

        for j in range(max_j + 1):

            E = L // (3 ** j)

            for a in BASES:

                if gcd(a, n) != 1:
                    continue

                x = pow(a, E, n)

                g = gcd(
                    n,
                    cyclotomic_value(3, x)
                )

                if not nontrivial(g, n):
                    continue

                op = multiplicative_order(a, p)
                oq = multiplicative_order(a, q)

                key = (
                    valuation(op,3),
                    valuation(oq,3),
                    valuation(E,3),
                    g,
                )

                correlation[key] += 1

    print(
        "v3(ord_p)  v3(ord_q)  v3(E)  factor   count"
    )

    line()

    for key, count in sorted(correlation.items()):

        vp, vq, ve, g = key

        print(
            f"{vp:10d}"
            f"{vq:11d}"
            f"{ve:8d}"
            f"{g:9d}"
            f"{count:7d}"
        )


# =============================================================================
# 12. FINAL SUMMARY
# =============================================================================

def final_summary():

    section("12. FINAL RESEARCH TARGET")

    print()
    print(
        "The strongest signal from the previous experiment was:"
    )

    print()
    print(
        "    gcd(n, Phi_3(a^E))"
    )

    print(
        "can separate the two hidden factors when their"
    )

    print(
        "3-adic order structures interact differently with E."
    )

    print()
    print(
        "This experiment now tests:"
    )

    print(
        "    1. exact v3(order) correlation"
    )

    print(
        "    2. exponent valuation dependence"
    )

    print(
        "    3. cubic versus ordinary order-finding"
    )

    print(
        "    4. generalized Phi_ell roots"
    )

    print(
        "    5. repeated behavior across targets"
    )

    print()
    print(
        "A particularly strong result would be a small, predictable"
    )

    print(
        "3-adic exponent strategy that gives a factor for a large"
    )

    print(
        "fraction of semiprimes without knowing p or q."
    )

    print()
    print(
        "NO CSV FILES WERE CREATED."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print()
    print("=" * 78)
    print("3-ADIC CUBIC-ROOT / GENERALIZED ROOT-OF-UNITY EXPERIMENT")
    print("=" * 78)

    section_exact_structure()
    section_representative()
    section_order_table()
    section_smooth_exponents()
    section_cubic_ladder()
    section_hit_mechanism()
    section_multitarget()
    section_compare_order()
    section_v3_sweep()
    section_generalized_roots()
    section_order_correlation()
    final_summary()

    print()
    print("=" * 78)
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()