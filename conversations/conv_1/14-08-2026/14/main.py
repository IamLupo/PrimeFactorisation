"""
==============================================================================
KAPPA EXPERIMENT 56
CONJUGATE DISCRIMINANT PRODUCT / N-ONLY ELIMINATION
NO CSV OUTPUT
==============================================================================

GOAL
----
Test the new algebraic identity obtained by restricting the discriminant

    D(s) = s^2 - 4n

to the cyclotomic root locus

    F(s) = s^2 + s + 1 = 0.

If w and w^2 are the two roots of F modulo ell, then

    d1 = w^2 - 4n
    d2 = w - 4n

because w^2 = -w-1.

Their product is independent of w:

    d1*d2 = 16*n^2 + 4*n + 1.

Define

    G(n) = 16*n^2 + 4*n + 1.

The main experimental question is whether

    Legendre(d1, ell) * Legendre(d2, ell)
      =
    Legendre(G(n), ell)

holds exactly, and whether this N-only invariant gives any
useful information about the two conjugate discriminant signs.

SECONDARY QUESTIONS
-------------------
1. Does G(n) correctly predict whether the two discriminants have
   equal or opposite Legendre symbols?

2. Can G(n) distinguish:
       (+1,+1)
       (-1,-1)
       (+1,-1)/(-1,+1)

3. What happens when one of d1,d2,G(n) is zero mod ell?

4. Does the phenomenon occur specifically for cyclotomic ell
   and disappear for random control primes?

5. Among the existing discriminant-sieve survivor sums s, does
   the conjugate-root classification produce any additional
   rejection beyond the ordinary QR condition?

IMPORTANT
---------
This experiment does NOT assume that the new identity gives a
factorization algorithm. It tests whether the eliminated N-only
quantity actually carries information.

No prime-pair enumeration beyond the supplied target controls.
No CSV files.
==============================================================================
"""

from __future__ import annotations

import math
import random
import time
from collections import Counter
from dataclasses import dataclass


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47,
]

CONTROL_COUNT = 13

TARGETS = [
    (3318013, 4042603),
    (2129167, 3402323),
    (2224517, 3978749),
    (3685051, 4020281),
    (2399627, 2452649),
    (2593039, 2996527),
    (2149859, 2772097),
    (2060543, 2514401),
    (2675423, 2883973),
    (2828887, 3960137),
    (3497381, 3793241),
    (2193509, 4011353),
]

# Existing cyclotomic prime family from the previous experiments.
CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723,
]

# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def sieve_primes(limit: int) -> list[int]:
    """Return all primes < limit."""
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * limit
    sieve[0:2] = b"\x00\x00"

    root = int(math.isqrt(limit - 1))
    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit:p] = b"\x00" * (
                ((limit - 1 - start) // p) + 1
            )

    return [i for i, flag in enumerate(sieve) if flag]


def factor_small(n: int) -> dict[int, int]:
    """Small trial-division factorization for F(r), used only here."""
    out: dict[int, int] = {}
    x = n
    d = 2

    while d * d <= x:
        while x % d == 0:
            out[d] = out.get(d, 0) + 1
            x //= d

        d = 3 if d == 2 else d + 2

    if x > 1:
        out[x] = out.get(x, 0) + 1

    return out


def legendre_symbol(a: int, p: int) -> int:
    """
    Legendre symbol (a/p):
        0 if p | a
        1 if a is a nonzero square mod p
       -1 otherwise
    """
    a %= p

    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1

    if value == p - 1:
        return -1

    raise ArithmeticError(
        f"Unexpected Legendre result: a={a}, p={p}, value={value}"
    )


def is_quadratic_residue(a: int, p: int) -> bool:
    return legendre_symbol(a, p) >= 0


def roots_of_F(ell: int) -> list[int]:
    """
    Return roots of x^2+x+1 modulo ell.

    For ell != 3 and ell == 1 mod 3 there are exactly two:
        w, w^2
    """
    if ell == 3:
        # x^2+x+1 = (x-1)^2 mod 3
        return [1]

    roots = []
    for x in range(ell):
        if (x * x + x + 1) % ell == 0:
            roots.append(x)

    return roots


def find_control_primes(
    excluded: set[int],
    count: int,
    upper: int = 10_000,
) -> list[int]:
    """Choose reproducible random control primes."""
    primes = [
        p
        for p in sieve_primes(upper)
        if p not in excluded and p >= 100
    ]

    rng = random.Random(SEED)
    rng.shuffle(primes)

    return sorted(primes[:count])


# ============================================================================
# ALGEBRA
# ============================================================================

def F(x: int) -> int:
    return x * x + x + 1


def G(n: int) -> int:
    """
    N-only elimination polynomial:

        G(n) = 16*n^2 + 4*n + 1
    """
    return 16 * n * n + 4 * n + 1


def discriminant_from_root(w: int, n: int, ell: int) -> int:
    """
    d(w) = w^2 - 4n mod ell.

    On F(w)=0 we could also use:
        d(w) = -w - 1 - 4n.
    """
    return (w * w - 4 * n) % ell


def reduced_discriminant_from_root(w: int, n: int, ell: int) -> int:
    """
    Same discriminant using w^2 = -w-1.
    This is an independent implementation of the same quantity.
    """
    return (-w - 1 - 4 * n) % ell


def conjugate_product_formula(
    w: int,
    n: int,
    ell: int,
) -> tuple[int, int, int]:
    """
    Return:

        d1 = d(w)
        d2 = d(w^2)
        d1*d2 mod ell

    using the two actual conjugate roots.
    """
    w2 = (w * w) % ell

    d1 = discriminant_from_root(w, n, ell)
    d2 = discriminant_from_root(w2, n, ell)
    product = (d1 * d2) % ell

    return d1, d2, product


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class PairStats:
    both_zero: int = 0
    one_zero: int = 0
    both_positive: int = 0
    both_negative: int = 0
    opposite_sign: int = 0
    product_residue_matches: int = 0
    product_residue_failures: int = 0


# ============================================================================
# CONTROL FAMILY
# ============================================================================

def build_control_family() -> list[int]:
    return find_control_primes(
        excluded=set(CYCLOTOMIC),
        count=CONTROL_COUNT,
        upper=10_000,
    )


# ============================================================================
# LOCAL ALGEBRA CHECK
# ============================================================================

def local_algebra_check(ell: int) -> tuple[bool, dict]:
    roots = roots_of_F(ell)

    info = {
        "roots": roots,
        "root_count": len(roots),
        "identity_ok": True,
    }

    for w in roots:
        if (w * w + w + 1) % ell != 0:
            info["identity_ok"] = False

        if pow(w, 3, ell) != 1:
            info["identity_ok"] = False

        if ell != 3 and w != 0 and pow(w, 2, ell) == 1:
            info["identity_ok"] = False

    return (
        info["identity_ok"],
        info,
    )


# ============================================================================
# MAIN IDENTITY TEST
# ============================================================================

def test_identity_for_modulus(
    ell: int,
    targets: list[tuple[int, int]],
) -> dict:
    roots = roots_of_F(ell)

    stats = {
        "ell": ell,
        "roots": roots,
        "targets": len(targets),
        "identity_pass": 0,
        "identity_fail": 0,
        "symbol_product_pass": 0,
        "symbol_product_fail": 0,
        "zero_cases": 0,
        "nonzero_cases": 0,
    }

    for p, q in targets:
        n = p * q

        # With two distinct roots, evaluate both.
        if len(roots) == 2:
            w1, w2 = roots
        elif len(roots) == 1:
            # Exceptional ell=3 case.
            w1 = roots[0]
            w2 = roots[0]
        else:
            stats["identity_fail"] += 1
            continue

        d1 = discriminant_from_root(w1, n, ell)
        d2 = discriminant_from_root(w2, n, ell)

        product_direct = (d1 * d2) % ell
        product_expected = G(n) % ell

        if product_direct == product_expected:
            stats["identity_pass"] += 1
        else:
            stats["identity_fail"] += 1

        s1 = legendre_symbol(d1, ell)
        s2 = legendre_symbol(d2, ell)
        sg = legendre_symbol(product_expected, ell)

        # For Legendre symbols:
        #
        #     (d1/ell)*(d2/ell) = (d1*d2/ell)
        #
        if s1 * s2 == sg:
            stats["symbol_product_pass"] += 1
        else:
            stats["symbol_product_fail"] += 1

        if d1 == 0 or d2 == 0:
            stats["zero_cases"] += 1
        else:
            stats["nonzero_cases"] += 1

    return stats


# ============================================================================
# DETAILED TARGET REPORT
# ============================================================================

def print_target_root_report(
    target_index: int,
    p: int,
    q: int,
    ell: int,
) -> PairStats:

    n = p * q
    roots = roots_of_F(ell)

    stats = PairStats()

    if len(roots) != 2:
        return stats

    w1, w2 = roots

    d1, d2, direct_product = conjugate_product_formula(w1, n, ell)
    expected = G(n) % ell

    # Independent reduced calculation.
    rd1 = reduced_discriminant_from_root(w1, n, ell)
    rd2 = reduced_discriminant_from_root(w2, n, ell)

    assert d1 == rd1
    assert d2 == rd2

    s1 = legendre_symbol(d1, ell)
    s2 = legendre_symbol(d2, ell)
    sg = legendre_symbol(expected, ell)

    print(
        f"  ell={ell:5d}"
        f" roots=({w1},{w2})"
        f" d1={d1:5d}"
        f" d2={d2:5d}"
        f" product={direct_product:5d}"
        f" G={expected:5d}"
        f" symbols=({s1:+d},{s2:+d})"
        f" G_symbol={sg:+d}"
        f" product_match={s1*s2 == sg}"
    )

    if d1 == 0 and d2 == 0:
        stats.both_zero += 1
    elif d1 == 0 or d2 == 0:
        stats.one_zero += 1
    elif s1 == 1 and s2 == 1:
        stats.both_positive += 1
    elif s1 == -1 and s2 == -1:
        stats.both_negative += 1
    elif s1 * s2 == -1:
        stats.opposite_sign += 1

    if s1 * s2 == sg:
        stats.product_residue_matches += 1
    else:
        stats.product_residue_failures += 1

    return stats


# ============================================================================
# N-ONLY CLASSIFICATION
# ============================================================================

def classify_pair_from_G(s1: int, s2: int, sg: int) -> str:
    """
    Interpret the pair of Legendre signs and the eliminated N-only sign.
    """
    if s1 == 0 or s2 == 0:
        return "ZERO_CASE"

    if s1 == 1 and s2 == 1:
        actual = "++"
    elif s1 == -1 and s2 == -1:
        actual = "--"
    else:
        actual = "+-"

    predicted = "same_sign" if sg == 1 else "opposite_sign"

    return f"{actual}/{predicted}"


# ============================================================================
# CANDIDATE SUM TEST
# ============================================================================

def candidate_sums_for_target(
    n: int,
    primes: list[int],
) -> list[int]:
    """
    Existing discriminant-sum domain:

        s = p+q

    with p,q in the configured prime interval.

    This deliberately scans p only and reconstructs q=n/p.
    It is a control/check path rather than the main discovery engine.
    """
    out = []

    prime_set = set(primes)

    for p in primes:
        if n % p != 0:
            continue

        q = n // p

        if q not in prime_set:
            continue

        if p > q:
            continue

        s = p + q

        if s % 2 == 0:
            out.append(s)

    return out


def discriminant_survivor_sums(
    n: int,
    ell_family: list[int],
    s_min: int,
    s_max: int,
) -> list[int]:
    """
    Apply only the QR discriminant condition:

        s^2 - 4n is QR mod every ell.

    This is a deliberately direct reference implementation.
    """
    survivors = []

    qr_sets: dict[int, set[int]] = {}

    for ell in ell_family:
        qr_sets[ell] = {
            x
            for x in range(ell)
            if legendre_symbol(x, ell) >= 0
        }

    for s in range(s_min, s_max + 1, 2):
        good = True
        for ell in ell_family:
            d = (s * s - 4 * n) % ell
            if d not in qr_sets[ell]:
                good = False
                break

        if good:
            survivors.append(s)

    return survivors


# ============================================================================
# EXPERIMENT
# ============================================================================

def main() -> None:
    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 56")
    print("CONJUGATE DISCRIMINANT PRODUCT / N-ONLY ELIMINATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # ------------------------------------------------------------------------
    # 1. TARGETS
    # ------------------------------------------------------------------------

    print("\n1. TARGETS")
    print("-" * 78)

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q
        s = p + q
        print(
            f"target {i:2d}: p={p} q={q} n={n} s={s}"
        )

    # ------------------------------------------------------------------------
    # 2. CYCLOTOMIC FAMILY
    # ------------------------------------------------------------------------

    print("\n2. CYCLOTOMIC FAMILY")
    print("-" * 78)

    for ell in CYCLOTOMIC:
        roots = roots_of_F(ell)
        print(
            f"ell={ell:5d}"
            f" roots={roots}"
            f" root_count={len(roots)}"
            f" ell%3={ell % 3}"
        )

    # ------------------------------------------------------------------------
    # 3. CONTROL FAMILY
    # ------------------------------------------------------------------------

    control = build_control_family()

    print("\n3. RANDOM CONTROL FAMILY")
    print("-" * 78)
    print(control)

    # ------------------------------------------------------------------------
    # 4. LOCAL ALGEBRA
    # ------------------------------------------------------------------------

    print("\n4. LOCAL ROOT ALGEBRA CHECK")
    print("-" * 78)

    local_failures = 0

    for ell in CYCLOTOMIC:
        ok, info = local_algebra_check(ell)
        print(
            f"ell={ell:5d}"
            f" roots={info['roots']}"
            f" identity_ok={info['identity_ok']}"
        )
        if not ok:
            local_failures += 1

    print(f"local failures = {local_failures}")
    print(f"status = {'PASS' if local_failures == 0 else 'FAIL'}")

    # ------------------------------------------------------------------------
    # 5. CENTRAL IDENTITY
    # ------------------------------------------------------------------------

    print("\n5. CENTRAL IDENTITY TEST")
    print("-" * 78)
    print("For F(w)=0:")
    print("    d(w) = w^2 - 4n = -w - 1 - 4n")
    print("")
    print("For conjugate roots w,w^2:")
    print("    d1*d2 = 16*n^2 + 4*n + 1 = G(n)")
    print("")

    identity_total = 0
    identity_pass = 0
    symbol_total = 0
    symbol_pass = 0

    for family_name, family in [
        ("CYCLOTOMIC", CYCLOTOMIC),
        ("CONTROL", control),
    ]:
        print(f"\n{family_name}")

        for ell in family:
            stats = test_identity_for_modulus(ell, TARGETS)

            identity_total += stats["identity_pass"] + stats["identity_fail"]
            identity_pass += stats["identity_pass"]

            symbol_total += (
                stats["symbol_product_pass"]
                + stats["symbol_product_fail"]
            )
            symbol_pass += stats["symbol_product_pass"]

            print(
                f"ell={ell:5d}"
                f" roots={stats['roots']}"
                f" identity={stats['identity_pass']}/"
                f"{stats['targets']}"
                f" symbol_product={stats['symbol_product_pass']}/"
                f"{stats['targets']}"
            )

    print(
        f"\nglobal identity checks = {identity_pass}/{identity_total}"
    )
    print(
        f"global symbol checks   = {symbol_pass}/{symbol_total}"
    )

    # ------------------------------------------------------------------------
    # 6. DETAILED TRUE-FACTOR TEST
    # ------------------------------------------------------------------------

    print("\n6. TRUE-FACTOR CONJUGATE DISCRIMINANTS")
    print("-" * 78)

    global_stats = PairStats()

    for target_index, (p, q) in enumerate(TARGETS, 1):
        print(
            f"\nTARGET {target_index}"
            f" p={p} q={q}"
            f" n={p*q}"
        )

        target_stats = PairStats()

        for ell in CYCLOTOMIC:
            s = print_target_root_report(
                target_index,
                p,
                q,
                ell,
            )

            target_stats.both_zero += s.both_zero
            target_stats.one_zero += s.one_zero
            target_stats.both_positive += s.both_positive
            target_stats.both_negative += s.both_negative
            target_stats.opposite_sign += s.opposite_sign
            target_stats.product_residue_matches += (
                s.product_residue_matches
            )
            target_stats.product_residue_failures += (
                s.product_residue_failures
            )

        global_stats.both_zero += target_stats.both_zero
        global_stats.one_zero += target_stats.one_zero
        global_stats.both_positive += target_stats.both_positive
        global_stats.both_negative += target_stats.both_negative
        global_stats.opposite_sign += target_stats.opposite_sign
        global_stats.product_residue_matches += (
            target_stats.product_residue_matches
        )
        global_stats.product_residue_failures += (
            target_stats.product_residue_failures
        )

        print(
            "  target summary:"
            f" both_zero={target_stats.both_zero}"
            f" one_zero={target_stats.one_zero}"
            f" ++={target_stats.both_positive}"
            f" --={target_stats.both_negative}"
            f" +-={target_stats.opposite_sign}"
        )

    print("\nGLOBAL TRUE-FACTOR SIGN SUMMARY")
    print(
        f"both_zero          = {global_stats.both_zero}"
    )
    print(
        f"one_zero           = {global_stats.one_zero}"
    )
    print(
        f"both positive      = {global_stats.both_positive}"
    )
    print(
        f"both negative      = {global_stats.both_negative}"
    )
    print(
        f"opposite signs     = {global_stats.opposite_sign}"
    )

    # ------------------------------------------------------------------------
    # 7. DIRECT N-ONLY TEST OF G(n)
    # ------------------------------------------------------------------------

    print("\n7. N-ONLY LEGENDRE SIGN OF G(n)")
    print("-" * 78)

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q

        print(f"\nTARGET {i} n={n}")
        for ell in CYCLOTOMIC:
            g_mod = G(n) % ell
            sg = legendre_symbol(g_mod, ell)

            print(
                f"  ell={ell:5d}"
                f" G(n) mod ell={g_mod:5d}"
                f" Legendre(G)={sg:+d}"
            )

    # ------------------------------------------------------------------------
    # 8. DOES G(n) PREDICT SAME/OPPOSITE SIGNS?
    # ------------------------------------------------------------------------

    print("\n8. G(n) SAME-SIGN / OPPOSITE-SIGN PREDICTION")
    print("-" * 78)

    prediction_total = 0
    prediction_correct = 0
    zero_cases = 0

    for target_index, (p, q) in enumerate(TARGETS, 1):
        n = p * q

        for ell in CYCLOTOMIC:
            roots = roots_of_F(ell)

            if len(roots) != 2:
                continue

            w1, w2 = roots

            d1 = discriminant_from_root(w1, n, ell)
            d2 = discriminant_from_root(w2, n, ell)
            sg = legendre_symbol(G(n), ell)

            s1 = legendre_symbol(d1, ell)
            s2 = legendre_symbol(d2, ell)

            if s1 == 0 or s2 == 0 or sg == 0:
                zero_cases += 1
                continue

            predicted_same = (sg == 1)
            actual_same = (s1 == s2)

            prediction_total += 1

            if predicted_same == actual_same:
                prediction_correct += 1

    prediction_rate = (
        prediction_correct / prediction_total
        if prediction_total
        else 0.0
    )

    print(
        f"nonzero prediction cases = {prediction_total}"
    )
    print(
        f"correct predictions       = {prediction_correct}"
    )
    print(
        f"prediction accuracy       = {prediction_rate:.6f}"
    )
    print(
        f"zero/degenerate cases     = {zero_cases}"
    )

    # ------------------------------------------------------------------------
    # 9. CONTROL COMPARISON
    # ------------------------------------------------------------------------

    print("\n9. CONTROL COMPARISON")
    print("-" * 78)

    def family_sign_summary(
        name: str,
        family: list[int],
    ) -> None:
        total = 0
        same = 0
        opposite = 0
        zeros = 0

        for p, q in TARGETS:
            n = p * q

            for ell in family:
                roots = roots_of_F(ell)

                if len(roots) != 2:
                    continue

                w1, w2 = roots

                d1 = discriminant_from_root(w1, n, ell)
                d2 = discriminant_from_root(w2, n, ell)

                s1 = legendre_symbol(d1, ell)
                s2 = legendre_symbol(d2, ell)

                if s1 == 0 or s2 == 0:
                    zeros += 1
                    continue

                total += 1

                if s1 == s2:
                    same += 1
                else:
                    opposite += 1

        print(
            f"{name}:"
            f" nonzero_pairs={total}"
            f" same={same}"
            f" opposite={opposite}"
            f" zero={zeros}"
        )

    family_sign_summary("CYCLOTOMIC", CYCLOTOMIC)

    # Random controls do not necessarily contain roots of F.
    # We therefore report only primes ell == 1 mod 3, where the
    # same two-root construction exists.
    control_root_primes = [
        ell for ell in control
        if ell % 3 == 1
    ]

    print(
        f"control primes ell=1 mod 3 = {control_root_primes}"
    )

    family_sign_summary(
        "CONTROL ell=1 mod 3",
        control_root_primes,
    )

    # ------------------------------------------------------------------------
    # 10. TRUE SUM VS FALSE SUM TEST
    # ------------------------------------------------------------------------

    print("\n10. TRUE SUM VS NEARBY FALSE SUMS")
    print("-" * 78)

    """
    Test whether the G(n)-predicted sign relation is something special
    about the true sum s=p+q.

    For each target:
      - true s
      - several nearby even sums
      - compare actual discriminants and the G(n) prediction.

    This avoids doing a huge candidate scan while checking whether
    the algebra is merely an identity valid for every integer s on
    the root locus, versus something correlated with the true factor.
    """

    OFFSETS = [-20, -18, -16, -14, -12, -10, -8, -6, -4, -2,
               2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

    for target_index, (p, q) in enumerate(TARGETS, 1):
        true_s = p + q
        n = p * q

        print(f"\nTARGET {target_index} true_s={true_s}")

        for offset in OFFSETS:
            s = true_s + offset

            if s <= 0:
                continue

            good = 0
            tested = 0

            for ell in CYCLOTOMIC:
                roots = roots_of_F(ell)

                if len(roots) != 2:
                    continue

                w1, w2 = roots

                # The actual candidate sum does not enter d1,d2 because
                # this section intentionally tests the root-locus algebra.
                d1 = (w1 * w1 - 4 * n) % ell
                d2 = (w2 * w2 - 4 * n) % ell

                sg = legendre_symbol(G(n), ell)

                s1 = legendre_symbol(d1, ell)
                s2 = legendre_symbol(d2, ell)

                if s1 == 0 or s2 == 0 or sg == 0:
                    continue

                tested += 1

                if (s1 == s2) == (sg == 1):
                    good += 1

            accuracy = good / tested if tested else 0.0

            print(
                f"  s={s:8d}"
                f" offset={offset:+4d}"
                f" relation_accuracy={accuracy:.6f}"
            )

    # ------------------------------------------------------------------------
    # 11. OPTIONAL EXISTING-SIEVE SURVIVOR ANALYSIS
    # ------------------------------------------------------------------------

    print("\n11. EXISTING DISCRIMINANT-SURVIVOR CROSS-CHECK")
    print("-" * 78)

    # We use a smaller range around the known true sum to keep runtime low.
    # This is NOT intended to reproduce the full earlier 2.2M-sum sieve.
    WINDOW = 2000
    survivor_stats = []

    for target_index, (p, q) in enumerate(TARGETS, 1):
        n = p * q
        true_s = p + q

        local_survivors = []

        for s in range(
            max(4_000_000, true_s - WINDOW),
            min(8_399_998, true_s + WINDOW) + 1,
            2,
        ):
            good = True

            for ell in CYCLOTOMIC:
                d = (s * s - 4 * n) % ell
                if legendre_symbol(d, ell) < 0:
                    good = False
                    break

            if good:
                local_survivors.append(s)

        survivor_stats.append((target_index, true_s, local_survivors))

        print(
            f"target {target_index:2d}:"
            f" true_s={true_s}"
            f" local_survivors={len(local_survivors)}"
            f" true_survives={true_s in local_survivors}"
        )

    # ------------------------------------------------------------------------
    # 12. FINAL SUMMARY
    # ------------------------------------------------------------------------

    print("\n12. FINAL SUMMARY")
    print("-" * 78)

    print(
        f"central identity checks = "
        f"{identity_pass}/{identity_total}"
    )

    print(
        f"Legendre product checks = "
        f"{symbol_pass}/{symbol_total}"
    )

    print(
        f"G(n) sign prediction     = "
        f"{prediction_correct}/{prediction_total}"
        f" ({prediction_rate:.6f})"
    )

    print(
        "\nCore identity:"
    )
    print(
        "  (w^2 - 4n)(w - 4n) = 16n^2 + 4n + 1"
    )

    print(
        "\nDecision interpretation:"
    )
    print(
        "  PASS identity + PASS Legendre relation:"
    )
    print(
        "      The N-only elimination identity is exact."
    )
    print(
        "  High same/opposite prediction accuracy:"
    )
    print(
        "      G(n) genuinely determines the product of the two"
    )
    print(
        "      conjugate discriminant Legendre symbols."
    )
    print(
        "  But this still does NOT determine which conjugate root"
    )
    print(
        "  is QR/non-QR individually."
    )
    print(
        "  That remaining asymmetry is the next possible source"
    )
    print(
        "  of additional information about s=p+q."
    )

    print("\nNo CSV files produced.")
    print("=" * 78)
    print("EXPERIMENT 56 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

