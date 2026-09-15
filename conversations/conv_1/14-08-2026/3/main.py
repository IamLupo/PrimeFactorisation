#!/usr/bin/env python3

import math
import random
import time
from collections import Counter, defaultdict

import sympy


# =============================================================================
# KAPPA EXPERIMENT 46
# DIRECT INVESTIGATION OF F(x) = x^2 + x + 1
# NO CSV OUTPUT
# =============================================================================

SEED = 20260814

R_VALUES = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

FAMILY_SIZE = len(R_VALUES)


def F(r: int) -> int:
    return r * r + r + 1


def factor_dict(n: int) -> dict[int, int]:
    return dict(sympy.factorint(n))


def multiplicative_order_if_unit(a: int, m: int):
    """
    Return multiplicative order of a mod m if gcd(a,m)=1.
    Return None otherwise.
    """
    if math.gcd(a, m) != 1:
        return None
    return sympy.n_order(a, m)


def print_rule_block(title: str):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def main():
    random.seed(SEED)

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 46")
    print("DIRECT CYCLOTOMIC INVESTIGATION OF F(x) = x^2 + x + 1")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print(f"random seed      = {SEED}")
    print(f"R values         = {R_VALUES}")
    print()

    # -------------------------------------------------------------------------
    # 1. POLYNOMIAL IDENTITY
    # -------------------------------------------------------------------------
    print_rule_block("1. POLYNOMIAL IDENTITY")

    print("F(x) = x^2 + x + 1")
    print("x^3 - 1 = (x - 1)(x^2 + x + 1)")
    print()
    print("Therefore:")
    print("    F(r) | (r^3 - 1)")
    print()
    print("The central question is whether this structure is visible")
    print("inside the KAPPA modulus family.")

    # -------------------------------------------------------------------------
    # 2. MODULUS TABLE
    # -------------------------------------------------------------------------
    print_rule_block("2. MODULUS / CYCLOTOMIC TABLE")

    print(
        f"{'r':>5} {'F(r)':>12} {'r^3-1':>16} "
        f"{'(r^3-1)/F(r)':>16} {'F factorization'}"
    )
    print("-" * 78)

    moduli = {}
    factorizations = {}

    for r in R_VALUES:
        m = F(r)
        numerator = r**3 - 1
        quotient = numerator // m

        moduli[r] = m
        fac = factor_dict(m)
        factorizations[r] = fac

        fac_text = " * ".join(
            f"{p}^{e}" if e != 1 else str(p)
            for p, e in fac.items()
        )

        print(
            f"{r:5d} {m:12d} {numerator:16d} "
            f"{quotient:16d} {fac_text}"
        )

    # -------------------------------------------------------------------------
    # 3. EXACT CUBIC-ROOT TEST
    # -------------------------------------------------------------------------
    print_rule_block("3. CUBIC-ROOT-OF-UNITY TEST")

    all_ok = True

    print(
        f"{'r':>5} {'m=F(r)':>12} "
        f"{'r^3 mod m':>14} {'gcd(r,m)':>10} {'status'}"
    )
    print("-" * 78)

    for r in R_VALUES:
        m = moduli[r]
        residue = pow(r, 3, m)
        g = math.gcd(r, m)
        ok = residue == 1

        if not ok:
            all_ok = False

        print(
            f"{r:5d} {m:12d} {residue:14d} "
            f"{g:10d} {'PASS' if ok else 'FAIL'}"
        )

    print()
    print(f"global cubic identity check = {'PASS' if all_ok else 'FAIL'}")

    # -------------------------------------------------------------------------
    # 4. FACTOR-LEVEL ORDER TEST
    # -------------------------------------------------------------------------
    print_rule_block("4. PRIME-FACTOR MULTIPLICATIVE ORDER")

    print(
        "For every prime ell dividing F(r), test the order of r modulo ell."
    )
    print()

    order_histogram = Counter()
    factor_mod_3 = Counter()
    exceptional_factors = []

    print(
        f"{'r':>5} {'ell':>10} {'exp':>5} "
        f"{'ell mod 3':>10} {'ord_ell(r)':>14} {'primitive cube root?'}"
    )
    print("-" * 78)

    for r in R_VALUES:
        fac = factorizations[r]

        for ell in sorted(fac):
            exponent = fac[ell]

            if math.gcd(r, ell) != 1:
                order = None
                primitive = False
            else:
                order = multiplicative_order_if_unit(r, ell)
                primitive = order == 3

            factor_mod_3[ell % 3] += 1

            if order is not None:
                order_histogram[order] += 1

            if order != 3:
                exceptional_factors.append((r, ell, order))

            print(
                f"{r:5d} {ell:10d} {exponent:5d} "
                f"{ell % 3:10d} {str(order):>14} "
                f"{'YES' if primitive else 'NO'}"
            )

    print()
    print("ORDER HISTOGRAM")
    for order, count in sorted(order_histogram.items()):
        print(f"  order {order:>4}: {count}")

    print()
    print("PRIME-FACTOR RESIDUE CLASS MOD 3")
    for cls in sorted(factor_mod_3):
        print(f"  ell ≡ {cls} (mod 3): {factor_mod_3[cls]}")

    print()
    if exceptional_factors:
        print("EXCEPTIONS TO ORDER-3 BEHAVIOUR")
        for r, ell, order in exceptional_factors:
            print(f"  r={r}, ell={ell}, order={order}")
    else:
        print("ALL PRIME FACTORS HAVE ORDER 3")

    # -------------------------------------------------------------------------
    # 5. DISTINCT PRIME FACTOR INVENTORY
    # -------------------------------------------------------------------------
    print_rule_block("5. DISTINCT PRIME FACTOR INVENTORY")

    all_factor_multiplicities = Counter()

    for r in R_VALUES:
        for ell, exponent in factorizations[r].items():
            all_factor_multiplicities[ell] += exponent

    distinct_factors = sorted(all_factor_multiplicities)

    print(f"distinct prime factors across family = {len(distinct_factors)}")
    print()

    print(
        f"{'ell':>10} {'ell mod 3':>10} "
        f"{'total exponent':>15} {'appears for r values'}"
    )
    print("-" * 78)

    appearances = defaultdict(list)

    for r in R_VALUES:
        for ell in factorizations[r]:
            appearances[ell].append(r)

    for ell in distinct_factors:
        print(
            f"{ell:10d} {ell % 3:10d} "
            f"{all_factor_multiplicities[ell]:15d} "
            f"{appearances[ell]}"
        )

    # -------------------------------------------------------------------------
    # 6. PAIRWISE GCD STRUCTURE OF F(r)
    # -------------------------------------------------------------------------
    print_rule_block("6. PAIRWISE GCD STRUCTURE OF F(r)")

    print(
        "Testing gcd(F(r), F(s)) for all distinct r,s."
    )
    print()

    nontrivial_gcds = []

    for i, r in enumerate(R_VALUES):
        for s in R_VALUES[i + 1:]:
            m1 = moduli[r]
            m2 = moduli[s]
            g = math.gcd(m1, m2)

            if g > 1:
                nontrivial_gcds.append((r, s, m1, m2, g))

    if not nontrivial_gcds:
        print("No nontrivial pairwise gcds found.")
    else:
        print(
            f"{'r':>5} {'s':>5} {'F(r)':>10} "
            f"{'F(s)':>10} {'gcd':>10}"
        )
        print("-" * 60)

        for r, s, m1, m2, g in nontrivial_gcds:
            print(f"{r:5d} {s:5d} {m1:10d} {m2:10d} {g:10d}")

    # -------------------------------------------------------------------------
    # 7. PRIME-FACTOR REUSE
    # -------------------------------------------------------------------------
    print_rule_block("7. PRIME-FACTOR REUSE ACROSS THE FAMILY")

    reused = [
        (ell, rs)
        for ell, rs in appearances.items()
        if len(rs) > 1
    ]

    if not reused:
        print("No prime factor is shared by two different F(r).")
    else:
        for ell, rs in reused:
            print(
                f"prime factor {ell} occurs in F(r) for r values {rs}"
            )

    # -------------------------------------------------------------------------
    # 8. LCM GROWTH
    # -------------------------------------------------------------------------
    print_rule_block("8. INCREMENTAL LCM GROWTH")

    running_lcm = 1

    print(
        f"{'step':>5} {'r':>5} {'F(r)':>12} "
        f"{'gcd(prev,F)':>14} {'LCM':>20} {'bits':>8}"
    )
    print("-" * 78)

    for step, r in enumerate(R_VALUES, start=1):
        m = moduli[r]
        g = math.gcd(running_lcm, m)
        running_lcm = math.lcm(running_lcm, m)

        print(
            f"{step:5d} {r:5d} {m:12d} "
            f"{g:14d} {running_lcm:20d} {running_lcm.bit_length():8d}"
        )

    # -------------------------------------------------------------------------
    # 9. ROOTS OF UNITY AS POLYNOMIAL REPRESENTATION
    # -------------------------------------------------------------------------
    print_rule_block("9. REDUCTION TO THE {1,r} BASIS")

    print(
        "Because r^2 + r + 1 = 0 mod F(r),"
    )
    print()
    print("    r^2 ≡ -r - 1")
    print()
    print("and therefore every polynomial in r can be reduced to")
    print()
    print("    A + B*r.")
    print()
    print("The experiment verifies several powers explicitly.")

    print()
    print(
        f"{'r':>5} {'F(r)':>12} {'r^2 mod F':>14} "
        f"{'-r-1 mod F':>14} {'r^3 mod F':>14}"
    )
    print("-" * 78)

    for r in R_VALUES:
        m = moduli[r]
        left = pow(r, 2, m)
        right = (-r - 1) % m
        cube = pow(r, 3, m)

        print(
            f"{r:5d} {m:12d} {left:14d} "
            f"{right:14d} {cube:14d}"
        )

    # -------------------------------------------------------------------------
    # 10. TEST HIGHER POWERS
    # -------------------------------------------------------------------------
    print_rule_block("10. PERIOD-3 POWER COLLAPSE")

    print(
        "Testing whether powers reduce with period 3 modulo F(r)."
    )
    print()

    print(
        f"{'r':>5} {'F(r)':>12} "
        f"{'r^4 mod m':>14} {'r mod m':>14} "
        f"{'r^5 mod m':>14} {'r^2 mod m':>14} "
        f"{'r^6 mod m':>14} {'1'}"
    )
    print("-" * 78)

    for r in R_VALUES:
        m = moduli[r]

        r4 = pow(r, 4, m)
        r1 = r % m
        r5 = pow(r, 5, m)
        r2 = pow(r, 2, m)
        r6 = pow(r, 6, m)

        print(
            f"{r:5d} {m:12d} "
            f"{r4:14d} {r1:14d} "
            f"{r5:14d} {r2:14d} "
            f"{r6:14d} "
            f"{'PASS' if r6 == 1 else 'FAIL'}"
        )

    # -------------------------------------------------------------------------
    # 11. SEARCH FOR A CLOSED-FORM GCD RELATION
    # -------------------------------------------------------------------------
    print_rule_block("11. GCD RELATION SEARCH")

    print(
        "A useful algebraic question is whether gcd(F(r), F(s))"
    )
    print(
        "is controlled by a simple relation between r and s."
    )
    print()

    print("Testing candidate relations:")
    print("  r == s mod ell")
    print("  r^2 + r + 1 == 0 mod ell")
    print("  ratio r * s^(-1) has small order when ell is shared")
    print()

    relation_hits = []

    for ell, rs in reused:
        for r in rs:
            for s in rs:
                if r >= s:
                    continue

                if math.gcd(s, ell) == 1:
                    ratio = (r * pow(s, -1, ell)) % ell
                    ratio_order = sympy.n_order(ratio, ell)
                else:
                    ratio = None
                    ratio_order = None

                relation_hits.append(
                    (ell, r, s, r % ell, s % ell, ratio, ratio_order)
                )

    if not relation_hits:
        print("No shared prime factors, so no cross-r relation was available.")
    else:
        for row in relation_hits:
            ell, r, s, rr, ss, ratio, ratio_order = row
            print(
                f"ell={ell} r={r} s={s} "
                f"r mod ell={rr} s mod ell={ss} "
                f"ratio={ratio} ratio_order={ratio_order}"
            )

    # -------------------------------------------------------------------------
    # 12. CONNECTION TO THE OBSERVED ~4x CANDIDATE DENSITY
    # -------------------------------------------------------------------------
    print_rule_block("12. CONNECTION TO THE OBSERVED ~4x FACTOR")

    print(
        "Earlier experiments observed approximately:"
    )
    print()
    print("    candidate fraction ≈ C / M")
    print()
    print("with C approaching roughly 4 for deeper prefixes.")
    print()
    print("This experiment does NOT claim that C=4 is explained.")
    print("Instead it prints the algebraic quantities needed to test")
    print("whether the cyclotomic structure is responsible.")
    print()

    # Simple family-level diagnostics.
    total_distinct_factors = len(distinct_factors)

    count_ell_1_mod3 = sum(
        1 for ell in distinct_factors
        if ell % 3 == 1
    )

    count_ell_2_mod3 = sum(
        1 for ell in distinct_factors
        if ell % 3 == 2
    )

    count_shared = len(reused)

    print(f"distinct prime factors                 = {total_distinct_factors}")
    print(f"prime factors ≡ 1 mod 3                = {count_ell_1_mod3}")
    print(f"prime factors ≡ 2 mod 3                = {count_ell_2_mod3}")
    print(f"shared prime factors across F(r)       = {count_shared}")
    print(f"nontrivial gcd pairs                   = {len(nontrivial_gcds)}")
    print(f"all order-3 factor tests passed        = {not exceptional_factors}")

    # -------------------------------------------------------------------------
    # 13. FINAL INTERPRETATION
    # -------------------------------------------------------------------------
    print_rule_block("13. INTERPRETATION")

    print(
        "The experiment distinguishes three levels of structure:"
    )
    print()
    print("LEVEL 1")
    print("    F(r) = r^2 + r + 1")
    print()
    print("LEVEL 2")
    print("    F(r) | r^3 - 1")
    print("    so r is a cube root of unity modulo F(r).")
    print()
    print("LEVEL 3")
    print("    prime factors ell | F(r)")
    print("    should inherit order-3 behaviour when ell is non-exceptional.")
    print()
    print("The next mathematical question is therefore:")
    print()
    print("    Can the observed KAPPA residue-filter behaviour")
    print("    be rewritten in the algebra of third roots of unity?")
    print()
    print("A positive result would be substantially stronger than")
    print("simply saying that the moduli happen to be F(r).")
    print()
    print("It would give an algebraic explanation for:")
    print("    * the residue classes")
    print("    * the dependence between moduli")
    print("    * the repeated candidate-density multiplier")
    print("    * the exceptional CRT growth")
    print()
    print("No candidate-pair enumeration is performed.")
    print("No CSV files are produced.")

    total_time = time.perf_counter() - t0

    print()
    print("=" * 78)
    print("EXPERIMENT 46 COMPLETE")
    print("=" * 78)
    print(f"total runtime = {total_time:.3f}s")


if __name__ == "__main__":
    main()

