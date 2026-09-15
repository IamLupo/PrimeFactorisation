#!/usr/bin/env python3
"""
==============================================================================
CYCLOTOMIC ORDER-SPECTRUM / ORDINARY-ORDER EQUIVALENCE EXPERIMENT
==============================================================================

Purpose
-------
Test whether cyclotomic projections

    gcd(n, Phi_l(a^E))

provide information beyond ordinary order projections

    gcd(n, a^(lE) - 1).

This experiment is deliberately PRINT-ONLY.
No CSV files are created.

Main questions
--------------
1. Is every cyclotomic hit already contained in an ordinary order hit?
2. Does Phi_l ever reveal a factor when a^(lE)-1 does not?
3. Does a multi-l cyclotomic "order spectrum" distinguish factors
   better than ordinary order projections?
4. Can several small primes l combine to reveal a factor earlier?
5. Does the answer change when p-1 and q-1 have different
   l-adic structures?
6. What is the cheapest successful exponent ladder?

Core identity
-------------
If

    Phi_l(x) = 0 mod p

then

    x^l = 1 mod p

so

    p | x^l - 1.

For x=a^E:

    gcd(n, Phi_l(a^E))
        divides
    gcd(n, a^(lE)-1).

The experiment checks whether the cyclotomic projection gives
STRICTLY more useful factor information despite this containment.

==============================================================================
"""

from __future__ import annotations

import math
import random
from itertools import combinations
from collections import defaultdict


# =============================================================================
# CONFIGURATION
# =============================================================================

BASES = list(range(2, 21))

# Small prime cyclotomic indices.
ELLS = [3, 5, 7, 11]

# Smooth-exponent bounds.
B_VALUES = [5, 7, 10, 12, 15, 20, 25]

# Number of bases actually used in larger experiments.
MAX_BASES = 12

# Small deterministic semiprime targets.
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

# Larger prime targets generated deterministically.
RANDOM_SEED = 20260814
LARGE_TARGET_COUNT = 12

# Avoid giant output from every possible tuple.
MAX_SIGNATURE_PRINT = 25
MAX_HIT_PRINT = 40


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def is_prime(n: int) -> bool:
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


def next_prime(n: int) -> int:
    x = max(2, n)
    while not is_prime(x):
        x += 1
    return x


def phi_prime(l: int, x: int, mod: int | None = None) -> int:
    """
    Phi_l(x) for prime l:

        Phi_l(x) = 1 + x + ... + x^(l-1)

    If mod is supplied, evaluate modulo mod efficiently.
    """
    if mod is None:
        total = 0
        power = 1
        for _ in range(l):
            total += power
            power *= x
        return total

    total = 0
    power = 1
    xm = x % mod
    for _ in range(l):
        total = (total + power) % mod
        power = (power * xm) % mod
    return total


def valuation(n: int, p: int) -> int:
    if n == 0:
        return 0
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def integer_lcm_upto(B: int) -> int:
    result = 1
    for k in range(2, B + 1):
        result = math.lcm(result, k)
    return result


def factor_small(n: int) -> dict[int, int]:
    result: dict[int, int] = {}
    d = 2
    x = n
    while d * d <= x:
        while x % d == 0:
            result[d] = result.get(d, 0) + 1
            x //= d
        d = 3 if d == 2 else d + 2
    if x > 1:
        result[x] = result.get(x, 0) + 1
    return result


def multiplicative_order(a: int, p: int) -> int | None:
    """
    Compute ord_p(a) for prime p using factorization of p-1.
    Only used for analysis/control because the hidden factor is known
    in these experiments.
    """
    if math.gcd(a, p) != 1:
        return None

    order = p - 1
    factors = factor_small(order)

    for prime_factor in factors:
        while order % prime_factor == 0:
            candidate = order // prime_factor
            if pow(a, candidate, p) == 1:
                order = candidate
            else:
                break

    return order


# =============================================================================
# CYCLOTOMIC / ORDER TESTS
# =============================================================================

def ordinary_gcd(n: int, a: int, E: int, ell: int) -> int:
    """
    gcd(n, a^(ell*E)-1)
    """
    x = pow(a, ell * E, n)
    return math.gcd(n, x - 1)


def cyclotomic_gcd(n: int, a: int, E: int, ell: int) -> int:
    """
    gcd(n, Phi_ell(a^E))
    """
    x = pow(a, E, n)
    value = phi_prime(ell, x, mod=n)
    return math.gcd(n, value)


def prime_factor_signature(g: int, p: int, q: int) -> str:
    if g == 1:
        return "NONE"
    if g == p:
        return f"p={p}"
    if g == q:
        return f"q={q}"
    if g == p * q:
        return "BOTH"
    return f"g={g}"


def stripped_factor(g: int, n: int) -> int:
    """
    Normalize a gcd into 1, p, q, or n when n is known factored.
    """
    return math.gcd(g, n)


# =============================================================================
# TARGET GENERATION
# =============================================================================

def deterministic_targets() -> list[tuple[int, int]]:
    return list(TARGETS)


def generate_large_targets() -> list[tuple[int, int]]:
    rng = random.Random(RANDOM_SEED)
    result: list[tuple[int, int]] = []

    # Use primes around different bit sizes.
    bit_sizes = [18, 20, 22, 24, 26, 28]

    for bits in bit_sizes:
        for _ in range(2):
            low = 1 << (bits - 1)
            high = (1 << bits) - 1

            # Find p.
            p = rng.randrange(low, high)
            p = next_prime(p)

            # Independent q.
            q = rng.randrange(low, high)
            q = next_prime(q)

            if p == q:
                q = next_prime(q + 1)

            if p > q:
                p, q = q, p

            result.append((p, q))

    return result[:LARGE_TARGET_COUNT]


# =============================================================================
# SECTION 1
# EXACT THEORETICAL CONTROL
# =============================================================================

def section_exact_containment() -> None:
    print("=" * 78)
    print("1. EXACT CYCLOTOMIC -> ORDINARY ORDER CONTAINMENT")
    print("=" * 78)
    print()
    print("For prime ell:")
    print()
    print("    Phi_ell(x)=0 mod p")
    print("        => x^ell = 1 mod p")
    print()
    print("Therefore:")
    print()
    print("    gcd(n, Phi_ell(a^E))")
    print("        divides")
    print("    gcd(n, a^(ell*E)-1)")
    print()

    failures = 0

    for ell in ELLS:
        for x in range(2, 15):
            for mod in [5, 7, 11, 13, 17, 19, 23, 29]:
                cyclo = phi_prime(ell, x, mod=mod)
                ordinary = (pow(x, ell, mod) - 1) % mod

                if cyclo == 0 and ordinary != 0:
                    failures += 1
                    print(
                        f"FAIL ell={ell} x={x} mod={mod}: "
                        "Phi=0 but x^ell-1 != 0"
                    )

    print(f"Containment failures: {failures}")
    print("PASS" if failures == 0 else "FAIL")
    print()


# =============================================================================
# SECTION 2
# REPRESENTATIVE TARGET
# =============================================================================

def section_representative_target(p: int, q: int) -> None:
    n = p * q

    print("=" * 78)
    print("2. REPRESENTATIVE TARGET")
    print("=" * 78)
    print()
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"n = {n}")
    print(f"p-1 factors = {factor_small(p - 1)}")
    print(f"q-1 factors = {factor_small(q - 1)}")
    print()

    for ell in ELLS:
        print(f"ell={ell}:")
        print(f"  p mod ell = {p % ell}")
        print(f"  q mod ell = {q % ell}")
        print(
            f"  v_{ell}(p-1) = {valuation(p - 1, ell)}, "
            f"v_{ell}(q-1) = {valuation(q - 1, ell)}"
        )

        shown = 0
        for a in BASES[:MAX_BASES]:
            op = multiplicative_order(a, p)
            oq = multiplicative_order(a, q)

            if op is None or oq is None:
                continue

            vp = valuation(op, ell)
            vq = valuation(oq, ell)

            if vp != 0 or vq != 0:
                print(
                    f"    a={a:2d}: "
                    f"ord_p={op:<8} v={vp} | "
                    f"ord_q={oq:<8} v={vq}"
                )
                shown += 1

        if shown == 0:
            print("    no ell-adic order structure for selected bases")
        print()


# =============================================================================
# SECTION 3
# ONE TARGET: DETAILED CYCLIC VS ORDINARY COMPARISON
# =============================================================================

def section_single_target_comparison(p: int, q: int) -> None:
    n = p * q

    print("=" * 78)
    print("3. SINGLE-TARGET CYCLOTOMIC VS ORDINARY ORDER")
    print("=" * 78)
    print()
    print(f"Target: {p} * {q} = {n}")
    print()

    for B in [5, 7, 10, 12, 15]:
        E = integer_lcm_upto(B)

        print(f"B={B}, E=lcm(1..B)={E}, bits(E)={E.bit_length()}")
        print(
            "ell  base   cyclo gcd   ordinary gcd   "
            "cyclo-factor   ordinary-factor   equal?"
        )
        print("-" * 78)

        local_hits = 0

        for ell in ELLS:
            for a in BASES[:MAX_BASES]:
                g_c = cyclotomic_gcd(n, a, E, ell)
                g_o = ordinary_gcd(n, a, E, ell)

                if g_c != 1 or g_o != 1:
                    local_hits += 1
                    print(
                        f"{ell:3d} {a:5d} "
                        f"{g_c:12d} {g_o:14d} "
                        f"{prime_factor_signature(g_c,p,q):>14} "
                        f"{prime_factor_signature(g_o,p,q):>16} "
                        f"{g_c == g_o}"
                    )

        print(f"nontrivial rows: {local_hits}")
        print()


# =============================================================================
# SECTION 4
# SEARCH FOR CYCLICOTOMIC-ONLY ADVANTAGE
# =============================================================================

def section_cyclotomic_only_search(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("4. SEARCH FOR CYCLOTOMIC-ONLY FACTOR INFORMATION")
    print("=" * 78)
    print()

    found: list[tuple] = []
    order_only: list[tuple] = []
    both: list[tuple] = []

    for p, q in targets:
        n = p * q

        for B in [5, 7, 10, 12, 15, 20]:
            E = integer_lcm_upto(B)

            for ell in ELLS:
                for a in BASES[:MAX_BASES]:
                    g_c = cyclotomic_gcd(n, a, E, ell)
                    g_o = ordinary_gcd(n, a, E, ell)

                    c_hit = g_c not in (1, n)
                    o_hit = g_o not in (1, n)

                    if c_hit and not o_hit:
                        found.append((p, q, B, ell, a, g_c, g_o))

                    elif o_hit and not c_hit:
                        order_only.append((p, q, B, ell, a, g_c, g_o))

                    elif c_hit and o_hit:
                        both.append((p, q, B, ell, a, g_c, g_o))

    print(f"cyclotomic-only hits : {len(found)}")
    print(f"ordinary-only hits   : {len(order_only)}")
    print(f"both                 : {len(both)}")
    print()

    if found:
        print("CYClOTOMIC-ONLY EXAMPLES:")
        for row in found[:MAX_HIT_PRINT]:
            print(row)
    else:
        print("No cyclotomic-only hits found.")

    print()

    if order_only:
        print("ORDINARY-ONLY EXAMPLES:")
        for row in order_only[:MAX_HIT_PRINT]:
            print(row)

    print()


# =============================================================================
# SECTION 5
# FULL ORDER SPECTRUM
# =============================================================================

def spectrum_for(
    n: int,
    p: int,
    q: int,
    a: int,
    E: int,
) -> dict:
    """
    Produce the factor results for all tested ell.

    spectrum entries:
        ordinary[ell] = gcd(n, a^(ell E)-1)
        cyclo[ell]    = gcd(n, Phi_ell(a^E))
    """
    ordinary = {}
    cyclo = {}

    for ell in ELLS:
        ordinary[ell] = ordinary_gcd(n, a, E, ell)
        cyclo[ell] = cyclotomic_gcd(n, a, E, ell)

    return {
        "ordinary": ordinary,
        "cyclo": cyclo,
    }


def factor_label(g: int, p: int, q: int) -> str:
    if g == 1:
        return "NONE"
    if g == p:
        return "p"
    if g == q:
        return "q"
    if g == p * q:
        return "BOTH"
    return "OTHER"


def section_order_spectrum(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("5. MULTI-ELL ORDER SPECTRUM")
    print("=" * 78)
    print()

    stats = defaultdict(int)

    representative_count = 0

    for p, q in targets:
        n = p * q

        for B in [5, 10, 15]:
            E = integer_lcm_upto(B)

            for a in BASES[:MAX_BASES]:
                spec = spectrum_for(n, p, q, a, E)

                ordinary_sig = tuple(
                    factor_label(spec["ordinary"][ell], p, q)
                    for ell in ELLS
                )

                cyclo_sig = tuple(
                    factor_label(spec["cyclo"][ell], p, q)
                    for ell in ELLS
                )

                stats[("ordinary", ordinary_sig)] += 1
                stats[("cyclo", cyclo_sig)] += 1

                if representative_count < MAX_SIGNATURE_PRINT:
                    print(
                        f"p={p:<8} q={q:<8} "
                        f"B={B:<2} a={a:<2}  "
                        f"ordinary={ordinary_sig}  "
                        f"cyclo={cyclo_sig}"
                    )
                    representative_count += 1

    print()
    print("Distinct ordinary signatures:", len(
        {k for k in stats if k[0] == "ordinary"}
    ))
    print("Distinct cyclotomic signatures:", len(
        {k for k in stats if k[0] == "cyclo"}
    ))
    print()


# =============================================================================
# SECTION 6
# INFORMATION CONTENT OF SPECTRA
# =============================================================================

def section_information_content(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("6. FACTOR-DISCRIMINATION POWER OF THE SPECTRA")
    print("=" * 78)
    print()

    ordinary_success = 0
    cyclo_success = 0
    combined_success = 0

    ordinary_factor_first: list[tuple] = []
    cyclo_factor_first: list[tuple] = []
    combined_factor_first: list[tuple] = []

    for p, q in targets:
        n = p * q

        found_o = None
        found_c = None
        found_combined = None

        attempt = 0

        for B in [5, 7, 10, 12, 15, 20, 25]:
            E = integer_lcm_upto(B)

            for ell in ELLS:
                for a in BASES[:MAX_BASES]:
                    attempt += 1

                    g_o = ordinary_gcd(n, a, E, ell)
                    g_c = cyclotomic_gcd(n, a, E, ell)

                    if found_o is None and g_o not in (1, n):
                        found_o = (attempt, B, ell, a, g_o)

                    if found_c is None and g_c not in (1, n):
                        found_c = (attempt, B, ell, a, g_c)

                    # Combined spectrum: either observation produces a factor.
                    if found_combined is None:
                        if g_o not in (1, n):
                            found_combined = ("ordinary", attempt, B, ell, a, g_o)
                        elif g_c not in (1, n):
                            found_combined = ("cyclo", attempt, B, ell, a, g_c)

                if found_o is not None and found_c is not None:
                    # Enough for this target.
                    pass

            if found_o is not None and found_c is not None:
                break

        if found_o is not None:
            ordinary_success += 1
            ordinary_factor_first.append((p, q, found_o))

        if found_c is not None:
            cyclo_success += 1
            cyclo_factor_first.append((p, q, found_c))

        if found_combined is not None:
            combined_success += 1
            combined_factor_first.append((p, q, found_combined))

    total = len(targets)

    print(f"Targets                         : {total}")
    print(f"Ordinary-order successes        : {ordinary_success}/{total}")
    print(f"Cyclotomic successes            : {cyclo_success}/{total}")
    print(f"Union of both                   : {combined_success}/{total}")
    print()

    print("FIRST ORDINARY FACTOR HITS:")
    for row in ordinary_factor_first[:MAX_SIGNATURE_PRINT]:
        print(row)

    print()
    print("FIRST CYCLOTOMIC FACTOR HITS:")
    for row in cyclo_factor_first[:MAX_SIGNATURE_PRINT]:
        print(row)

    print()


# =============================================================================
# SECTION 7
# DOES MULTIPLE ELL HELP?
# =============================================================================

def section_multi_ell_union(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("7. DOES MULTIPLE ell OUTPERFORM THE BEST SINGLE ell?")
    print("=" * 78)
    print()

    for ell_subset_size in range(1, len(ELLS) + 1):
        subset_results = {}

        for subset in combinations(ELLS, ell_subset_size):
            successes = 0

            for p, q in targets:
                n = p * q
                got = False

                for B in [5, 7, 10, 12, 15, 20]:
                    E = integer_lcm_upto(B)

                    for a in BASES[:MAX_BASES]:
                        for ell in subset:
                            g = cyclotomic_gcd(n, a, E, ell)

                            if g not in (1, n):
                                got = True
                                break

                        if got:
                            break

                    if got:
                        break

                if got:
                    successes += 1

            subset_results[subset] = successes

        best = max(subset_results.values())
        best_subsets = [
            subset for subset, score in subset_results.items()
            if score == best
        ]

        print(
            f"subset size={ell_subset_size}: "
            f"best={best}/{len(targets)} "
            f"subsets={best_subsets}"
        )

    print()


# =============================================================================
# SECTION 8
# 3-ADIC GENERALIZATION
# =============================================================================

def section_ell_adic_explanation(p: int, q: int) -> None:
    print("=" * 78)
    print("8. ELL-ADIC ORDER PROJECTION")
    print("=" * 78)
    print()

    print(f"Target: p={p}, q={q}")
    print()

    for ell in ELLS:
        print(f"ell={ell}")
        print(
            "base   ord_p   v_ell(ord_p)   "
            "ord_q   v_ell(ord_q)"
        )
        print("-" * 60)

        for a in BASES[:MAX_BASES]:
            op = multiplicative_order(a, p)
            oq = multiplicative_order(a, q)

            if op is None or oq is None:
                continue

            vp = valuation(op, ell)
            vq = valuation(oq, ell)

            if vp or vq:
                print(
                    f"{a:4d} {op:8d} {vp:14d} "
                    f"{oq:8d} {vq:14d}"
                )

        print()


# =============================================================================
# SECTION 9
# NEW TEST:
# CAN WE RECONSTRUCT ORDINARY INFORMATION FROM CYCLOTOMIC DATA?
# =============================================================================

def section_reconstruction(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("9. CAN CYCLOTOMIC PROJECTIONS RECONSTRUCT ORDINARY ORDER HITS?")
    print("=" * 78)
    print()

    exact = 0
    total = 0
    examples = []

    for p, q in targets:
        n = p * q

        for B in [5, 7, 10, 15]:
            E = integer_lcm_upto(B)

            for a in BASES[:MAX_BASES]:
                for ell in ELLS:
                    g_c = cyclotomic_gcd(n, a, E, ell)
                    g_oe = ordinary_gcd(n, a, E, ell)

                    # The cyclotomic factor must divide the ordinary factor.
                    # Check divisibility in the actual gcd values.
                    total += 1

                    if g_oe % g_c == 0:
                        exact += 1
                    elif len(examples) < 10:
                        examples.append(
                            (p, q, B, a, ell, g_c, g_oe)
                        )

    print(f"Containment checks = {exact}/{total}")
    if examples:
        print("Unexpected examples:")
        for row in examples:
            print(row)
    else:
        print("PASS: every cyclotomic gcd divides the corresponding ordinary gcd.")

    print()


# =============================================================================
# SECTION 10
# SEARCH FOR A GENUINE ADVANTAGE THROUGH COMBINATIONS
# =============================================================================

def section_combined_signatures(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("10. COMBINED SPECTRUM SEARCH")
    print("=" * 78)
    print()

    print(
        "We now compare:"
        "\n"
        "  ordinary vector = [gcd(n,a^(ellE)-1)]"
        "\n"
        "  cyclotomic vector = [gcd(n,Phi_ell(a^E))]"
        "\n"
        "  mixed vector = both together"
    )
    print()

    ordinary_factor_hits = 0
    cyclo_factor_hits = 0
    mixed_factor_hits = 0

    for p, q in targets:
        n = p * q

        ordinary_vector = []
        cyclo_vector = []

        for B in [5, 10, 15]:
            E = integer_lcm_upto(B)

            for a in BASES[:MAX_BASES]:
                for ell in ELLS:
                    ordinary_vector.append(
                        ordinary_gcd(n, a, E, ell)
                    )
                    cyclo_vector.append(
                        cyclotomic_gcd(n, a, E, ell)
                    )

        ordinary_factor_hits += any(
            g not in (1, n) for g in ordinary_vector
        )

        cyclo_factor_hits += any(
            g not in (1, n) for g in cyclo_vector
        )

        mixed_factor_hits += any(
            g not in (1, n)
            for g in ordinary_vector + cyclo_vector
        )

    total = len(targets)

    print(f"ordinary vector factor coverage = {ordinary_factor_hits}/{total}")
    print(f"cyclotomic vector factor coverage = {cyclo_factor_hits}/{total}")
    print(f"mixed vector factor coverage = {mixed_factor_hits}/{total}")
    print()


# =============================================================================
# SECTION 11
# SCALING
# =============================================================================

def section_scaling(targets: list[tuple[int, int]]) -> None:
    print("=" * 78)
    print("11. TARGET-SIZE SCALING")
    print("=" * 78)
    print()

    print(
        "bits(n)   p_bits q_bits   "
        "ordinary_hit   cyclo_hit   first_method"
    )
    print("-" * 78)

    for p, q in targets:
        n = p * q

        ordinary = None
        cyclo = None

        attempt = 0

        for B in [5, 7, 10, 12, 15, 20, 25]:
            E = integer_lcm_upto(B)

            for ell in ELLS:
                for a in BASES[:MAX_BASES]:
                    attempt += 1

                    if ordinary is None:
                        g = ordinary_gcd(n, a, E, ell)
                        if g not in (1, n):
                            ordinary = (attempt, B, ell, a, g)

                    if cyclo is None:
                        g = cyclotomic_gcd(n, a, E, ell)
                        if g not in (1, n):
                            cyclo = (attempt, B, ell, a, g)

                    if ordinary is not None and cyclo is not None:
                        break

                if ordinary is not None and cyclo is not None:
                    break

            if ordinary is not None and cyclo is not None:
                break

        if ordinary is None and cyclo is None:
            method = "NONE"
        elif ordinary is None:
            method = "CYCLO"
        elif cyclo is None:
            method = "ORDER"
        else:
            method = "CYCLO" if cyclo[0] < ordinary[0] else "ORDER"

        print(
            f"{n.bit_length():7d} "
            f"{p.bit_length():6d} {q.bit_length():6d}   "
            f"{str(ordinary):>20} "
            f"{str(cyclo):>20}   "
            f"{method}"
        )

    print()


# =============================================================================
# SECTION 12
# KEY CONTROL:
# IS CYCLOTOMIC TEST COMPUTATIONALLY CHEAPER?
# =============================================================================

def section_cost_comparison() -> None:
    print("=" * 78)
    print("12. COMPUTATIONAL FORM COMPARISON")
    print("=" * 78)
    print()

    print("Cyclotomic test:")
    print("    x = a^E mod n")
    print("    Phi_ell(x) mod n")
    print("    gcd(n, Phi_ell(x))")
    print()

    print("Ordinary order test:")
    print("    y = a^(ell E) mod n")
    print("    gcd(n, y-1)")
    print()

    print(
        "Both require modular exponentiation.\n"
        "The cyclotomic test does not eliminate the need for\n"
        "order-related exponentiation; it only post-processes a^E."
    )
    print()

    print("The following experiment therefore focuses on INFORMATION,\n"
          "not merely arithmetic operation count.")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print()
    print("=" * 78)
    print("CYCLOTOMIC ORDER-SPECTRUM / ORDINARY-ORDER EQUIVALENCE EXPERIMENT")
    print("=" * 78)
    print()
    print("NO CSV FILES WILL BE CREATED.")
    print()

    targets = deterministic_targets()
    targets += generate_large_targets()

    representative = (157, 163)

    section_exact_containment()
    section_representative_target(*representative)
    section_single_target_comparison(*representative)
    section_cyclotomic_only_search(targets)
    section_order_spectrum(targets)
    section_information_content(targets)
    section_multi_ell_union(targets)
    section_ell_adic_explanation(*representative)
    section_reconstruction(targets)
    section_combined_signatures(targets)
    section_scaling(targets)
    section_cost_comparison()

    print("=" * 78)
    print("FINAL QUESTIONS FOR THE OUTPUT")
    print("=" * 78)
    print()
    print("1. Are there any cyclotomic-only factor hits?")
    print("2. Does adding ell values improve coverage?")
    print("3. Does the cyclotomic spectrum distinguish any target")
    print("   that the ordinary-order spectrum cannot?")
    print("4. Does cyclotomic information ever give a first hit")
    print("   substantially earlier than ordinary order information?")
    print("5. Does this remain true as n grows?")
    print()
    print("DONE")
    print("=" * 78)


if __name__ == "__main__":
    main()

