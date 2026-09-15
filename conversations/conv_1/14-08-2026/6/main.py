#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 49
CUBIC-ROOT-OF-UNITY TWO-COORDINATE FACTOR EQUATION
NO CSV OUTPUT
==============================================================================

Goal
----
Test the algebraic representation induced by

    F(r) = r^2 + r + 1 = 0

and therefore

    omega^2 = -omega - 1
    omega^3 = 1.

For a prime ell | F(r), represent residues in the basis {1, omega}:

    x = a + b*omega

For

    p = a + b*omega
    q = c + d*omega

we have

    pq =
        (a*c - b*d)
        + (a*d + b*c - b*d) * omega.

Therefore pq == n (mod ell) is equivalent to:

    a*c - b*d             == n (mod ell)
    a*d + b*c - b*d       == 0 (mod ell)

The experiment compares this two-coordinate algebraic formulation
with the ordinary scalar residue relation

    p*q == n (mod ell).

It also tests the algebraic inverse explicitly.

No candidate-pair enumeration is performed.
No CSV files are produced.
"""

from __future__ import annotations

import math
import random
import time
from collections import Counter, defaultdict
from typing import Dict, List, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

TARGETS = 12

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29,
    31, 37, 41, 43, 47,
]

# Use the same control targets used repeatedly in the preceding experiments.
CONTROL_PAIRS = [
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


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    """Return all primes in [lo, hi)."""
    if hi <= 2 or hi <= lo:
        return []

    limit = hi
    sieve = bytearray(b"\x01") * limit
    sieve[0:2] = b"\x00\x00"

    root = math.isqrt(limit - 1)

    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit:p] = b"\x00" * (
                ((limit - 1 - start) // p) + 1
            )

    return [p for p in range(lo, hi) if sieve[p]]


# ============================================================================
# CYCLOTOMIC FAMILY
# ============================================================================

def F(r: int) -> int:
    return r * r + r + 1


def factor_small(n: int) -> Dict[int, int]:
    """Trial-factor the small cyclotomic moduli."""
    out: Dict[int, int] = {}
    x = n
    d = 2

    while d * d <= x:
        while x % d == 0:
            out[d] = out.get(d, 0) + 1
            x //= d
        d += 1

    if x > 1:
        out[x] = out.get(x, 0) + 1

    return out


# ============================================================================
# CUBIC ROOT ALGEBRA
# ============================================================================

def basis_mul(
    a: int,
    b: int,
    c: int,
    d: int,
    ell: int,
) -> Tuple[int, int]:
    """
    Multiply

        (a + b*w)(c + d*w)

    under w^2 = -w - 1.

    Result:

        real = ac - bd
        omega = ad + bc - bd
    """
    real = (a * c - b * d) % ell
    omega = (a * d + b * c - b * d) % ell
    return real, omega


def basis_pow(a: int, b: int, exponent: int, ell: int) -> Tuple[int, int]:
    """Exponentiate a+b*w in the two-coordinate algebra."""
    result = (1, 0)
    base = (a % ell, b % ell)

    e = exponent

    while e:
        if e & 1:
            result = basis_mul(
                result[0], result[1],
                base[0], base[1],
                ell,
            )
        base = basis_mul(
            base[0], base[1],
            base[0], base[1],
            ell,
        )
        e >>= 1

    return result


def basis_inverse(
    a: int,
    b: int,
    ell: int,
) -> Tuple[int, int]:
    """
    Invert a+b*w using the norm

        N(a+b*w) = a^2 - a*b + b^2.

    For the quadratic Eisenstein relation w^2+w+1=0,

        (a+b*w)(a+b*w^2)
    = a^2 - ab + b^2.

    Also

        w^2 = -1-w.

    The conjugate can therefore be written

        a + b*w^2 = (a-b) - b*w.

    Hence

        (a+b*w)^(-1)
        = ((a-b) - b*w) / N.
    """
    a %= ell
    b %= ell

    norm = (a * a - a * b + b * b) % ell

    if norm == 0:
        raise ZeroDivisionError(
            f"non-unit element ({a},{b}) modulo {ell}"
        )

    inv_norm = pow(norm, -1, ell)

    inv_a = ((a - b) * inv_norm) % ell
    inv_b = (-b * inv_norm) % ell

    # Independent verification.
    chk_a, chk_b = basis_mul(
        a, b,
        inv_a, inv_b,
        ell,
    )

    if chk_a != 1 or chk_b != 0:
        raise AssertionError(
            f"inverse verification failed for ({a},{b}) mod {ell}"
        )

    return inv_a, inv_b


def scalar_to_basis(x: int, ell: int) -> Tuple[int, int]:
    """Ordinary integer x embedded as x + 0*w."""
    return x % ell, 0


# ============================================================================
# COORDINATE EXTRACTION
# ============================================================================

def r_root_coordinates(r: int, ell: int) -> Tuple[int, int]:
    """
    In the quotient defined by ell | F(r), the distinguished root is
    simply the residue r.

    Therefore omega can be represented by r:

        omega = r mod ell.

    The scalar residue x can be viewed as x + 0*omega.

    This routine returns the coordinate representation of r itself
    in the ambient residue field, which is (0,1) conceptually,
    while also retaining the concrete scalar realization r mod ell.
    """
    if (r * r + r + 1) % ell != 0:
        raise ValueError(
            f"{ell} does not divide F({r})"
        )

    return 0, 1


# ============================================================================
# LOCAL ALGEBRA TESTS
# ============================================================================

def test_cyclotomic_algebra() -> None:
    print()
    print("=" * 78)
    print("1. LOCAL CUBIC-ROOT ALGEBRA CHECK")
    print("=" * 78)

    failures = 0

    for r in R_VALUES:
        fr = F(r)
        factors = factor_small(fr)

        print(f"r={r:3d} F(r)={fr:6d} factors={factors}")

        for ell in factors:
            if ell == 3:
                # ell=3 is exceptional: r == 1 (mod 3) for these r.
                print(
                    f"    ell={ell:4d} exceptional factor "
                    f"(order-1 case)"
                )
                continue

            omega = r % ell

            if (omega * omega + omega + 1) % ell != 0:
                print("    FAIL: omega^2 + omega + 1 != 0")
                failures += 1
                continue

            if pow(omega, 3, ell) != 1:
                print("    FAIL: omega^3 != 1")
                failures += 1
                continue

            # Test period-3 powers.
            ok = True
            for k in range(1, 10):
                expected = pow(omega, k % 3 or 3, ell)
                if pow(omega, k, ell) != expected:
                    ok = False
                    break

            print(
                f"    ell={ell:4d} omega={omega:4d} "
                f"order_candidate=3 powers={'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(f"local algebra failures = {failures}")
    print(f"status = {'PASS' if failures == 0 else 'FAIL'}")


# ============================================================================
# TWO-COORDINATE FACTOR EQUATION
# ============================================================================

def algebraic_factor_equation(
    p: int,
    q: int,
    n: int,
    r: int,
    ell: int,
) -> Tuple[Tuple[int, int], Tuple[int, int], bool]:
    """
    Represent p and q as scalar elements:

        p -> (p,0)
        q -> (q,0)

    and multiply in the cubic-root algebra.

    Because both are ordinary integers, the omega coordinate should
    vanish automatically.

    This is a deliberately simple control: it confirms the algebraic
    multiplication law is consistent with ordinary modular multiplication.
    """
    pa, pb = scalar_to_basis(p, ell)
    qa, qb = scalar_to_basis(q, ell)

    prod = basis_mul(pa, pb, qa, qb, ell)

    expected = (n % ell, 0)

    return (
        (pa, pb),
        (qa, qb),
        prod == expected,
    )


def test_two_coordinate_factor_equation(
    targets: List[Tuple[int, int, int]],
) -> None:
    print()
    print("=" * 78)
    print("2. TWO-COORDINATE FACTOR EQUATION")
    print("=" * 78)

    checked = 0
    failures = 0

    for target_id, (p, q, n) in enumerate(targets, start=1):
        print()
        print(
            f"TARGET {target_id:2d} "
            f"p={p} q={q} n={n}"
        )

        for r in R_VALUES:
            factors = factor_small(F(r))

            for ell in factors:
                if ell == 3:
                    continue

                pa, pb = scalar_to_basis(p, ell)
                qa, qb = scalar_to_basis(q, ell)

                prod_a, prod_b = basis_mul(
                    pa, pb,
                    qa, qb,
                    ell,
                )

                scalar_ok = (p * q - n) % ell == 0
                algebra_ok = (
                    prod_a == n % ell
                    and prod_b == 0
                )

                checked += 1

                if scalar_ok != algebra_ok:
                    failures += 1

                print(
                    f"  ell={ell:4d} "
                    f"p=({pa:4d},{pb:2d}) "
                    f"q=({qa:4d},{qb:2d}) "
                    f"product=({prod_a:4d},{prod_b:2d}) "
                    f"scalar={scalar_ok} "
                    f"algebra={algebra_ok}"
                )

    print()
    print(f"checked equations = {checked}")
    print(f"equation failures = {failures}")
    print(
        "status = "
        + ("PASS" if failures == 0 else "FAIL")
    )


# ============================================================================
# NONTRIVIAL BASIS ELEMENT EXPERIMENT
# ============================================================================

def test_nontrivial_basis_elements() -> None:
    print()
    print("=" * 78)
    print("3. NONTRIVIAL {1,omega} MULTIPLICATION")
    print("=" * 78)

    failures = 0
    cases = 0

    samples = [
        (2, 3, 5, 7),
        (4, 1, 6, 2),
        (7, 5, 3, 9),
        (11, 4, 8, 6),
    ]

    for r in R_VALUES:
        factors = factor_small(F(r))

        for ell in factors:
            if ell == 3:
                continue

            print()
            print(f"r={r:3d} ell={ell:4d}")

            for a, b, c, d in samples:
                left = basis_mul(a, b, c, d, ell)

                # Explicit expansion:
                # (a+bw)(c+dw)
                #
                # = ac + (ad+bc)w + bd*w^2
                # = ac-bd + (ad+bc-bd)w
                right = (
                    (a * c - b * d) % ell,
                    (a * d + b * c - b * d) % ell,
                )

                cases += 1

                ok = left == right
                if not ok:
                    failures += 1

                print(
                    f"  ({a}+{b}w)*({c}+{d}w)"
                    f" -> {left}"
                    f" expected={right}"
                    f" {'PASS' if ok else 'FAIL'}"
                )

    print()
    print(f"nontrivial basis cases = {cases}")
    print(f"failures = {failures}")
    print(
        "status = "
        + ("PASS" if failures == 0 else "FAIL")
    )


# ============================================================================
# INVERSE STRUCTURE
# ============================================================================

def test_inverse_structure() -> None:
    print()
    print("=" * 78)
    print("4. ALGEBRAIC INVERSE / NORM STRUCTURE")
    print("=" * 78)

    failures = 0
    cases = 0

    samples = [
        (1, 1),
        (2, 1),
        (3, 2),
        (5, 4),
        (7, 3),
    ]

    for r in R_VALUES:
        factors = factor_small(F(r))

        for ell in factors:
            if ell == 3:
                continue

            print()
            print(f"r={r:3d} ell={ell:4d}")

            for a, b in samples:
                try:
                    inv = basis_inverse(a, b, ell)
                except ZeroDivisionError:
                    print(
                        f"  ({a},{b}) non-unit mod {ell} "
                        f"(skipped)"
                    )
                    continue

                norm = (
                    a * a
                    - a * b
                    + b * b
                ) % ell

                product = basis_mul(
                    a, b,
                    inv[0], inv[1],
                    ell,
                )

                cases += 1

                ok = product == (1, 0)

                if not ok:
                    failures += 1

                print(
                    f"  x=({a},{b}) "
                    f"norm={norm:5d} "
                    f"inv=({inv[0]:5d},{inv[1]:5d}) "
                    f"x*inv={product} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

    print()
    print(f"inverse cases = {cases}")
    print(f"failures = {failures}")
    print(
        "status = "
        + ("PASS" if failures == 0 else "FAIL")
    )


# ============================================================================
# SIGNATURE COMPARISON
# ============================================================================

def cyclotomic_prime_factors() -> List[int]:
    factors = set()

    for r in R_VALUES:
        for ell in factor_small(F(r)):
            if ell != 3:
                factors.add(ell)

    return sorted(factors)


def signature_scalar(
    p: int,
    n: int,
    ell: int,
) -> int:
    """
    Ordinary scalar signature:

        n * p^{-1} mod ell.
    """
    p_mod = p % ell

    if p_mod == 0:
        return -1

    return (n % ell) * pow(p_mod, -1, ell) % ell


def signature_basis_embedding(
    p: int,
    n: int,
    ell: int,
) -> Tuple[int, int]:
    """
    Embed the scalar signature into the basis:

        q = q_scalar + 0*w.
    """
    q = signature_scalar(p, n, ell)

    if q < 0:
        return (-1, -1)

    return (q, 0)


def compare_signature_forms(
    targets: List[Tuple[int, int, int]],
) -> None:
    print()
    print("=" * 78)
    print("5. SCALAR SIGNATURE VS TWO-COORDINATE SIGNATURE")
    print("=" * 78)

    ells = cyclotomic_prime_factors()

    print()
    print("non-exceptional cyclotomic prime factors:")
    print(ells)

    total = 0
    failures = 0

    for target_id, (p, q, n) in enumerate(targets, start=1):
        print()
        print(
            f"TARGET {target_id:2d}: "
            f"true p={p} q={q}"
        )

        for ell in ells:
            scalar_q = signature_scalar(p, n, ell)

            if scalar_q < 0:
                print(
                    f"  ell={ell:4d} p divisible by ell "
                    f"(non-unit branch)"
                )
                continue

            basis_q = signature_basis_embedding(
                p,
                n,
                ell,
            )

            direct = q % ell

            ok = (
                scalar_q == direct
                and basis_q == (direct, 0)
            )

            total += 1

            if not ok:
                failures += 1

            print(
                f"  ell={ell:4d} "
                f"q_direct={direct:6d} "
                f"q_scalar={scalar_q:6d} "
                f"q_basis={basis_q} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(f"signature comparisons = {total}")
    print(f"failures = {failures}")
    print(
        "status = "
        + ("PASS" if failures == 0 else "FAIL")
    )


# ============================================================================
# NEW ALGEBRAIC QUESTION
# ============================================================================

def test_norm_relation(
    targets: List[Tuple[int, int, int]],
) -> None:
    """
    Examine whether the norm expression for the basis element
    corresponding to p/q produces a useful relation.

    For x = a+b*w:

        N(x)=a^2-ab+b^2.

    The experiment records the norm of:
        1) p as scalar (p,0)
        2) q as scalar (q,0)
        3) np^{-1} as scalar

    Since this is initially an embedding into the scalar axis,
    the norm is simply the square of the scalar residue.

    This section is intentionally diagnostic: we are looking for
    nontrivial reductions rather than assuming they exist.
    """
    print()
    print("=" * 78)
    print("6. NORM DIAGNOSTIC")
    print("=" * 78)

    ells = cyclotomic_prime_factors()

    for target_id, (p, q, n) in enumerate(targets, start=1):
        print()
        print(
            f"TARGET {target_id:2d} "
            f"n={n}"
        )

        interesting = 0

        for ell in ells:
            if p % ell == 0:
                continue

            qp = signature_scalar(p, n, ell)

            norm_p = (p * p) % ell
            norm_q = (q * q) % ell
            norm_sig = (qp * qp) % ell

            relation = (
                norm_sig == norm_q
                and (
                    (norm_p * norm_q - n * n)
                    % ell
                ) == 0
            )

            if relation:
                interesting += 1

            print(
                f"  ell={ell:4d} "
                f"N(p)={norm_p:6d} "
                f"N(q)={norm_q:6d} "
                f"N(sig)={norm_sig:6d} "
                f"relation={relation}"
            )

        print(
            f"  norm relations holding = "
            f"{interesting}/{len(ells)}"
        )


# ============================================================================
# TARGET GENERATION
# ============================================================================

def build_targets(primes: List[int]) -> List[Tuple[int, int, int]]:
    rng = random.Random(SEED)

    targets: List[Tuple[int, int, int]] = []

    available = set(primes)

    while len(targets) < TARGETS:
        p, q = rng.sample(primes, 2)

        if p == q:
            continue

        if p not in available or q not in available:
            continue

        n = p * q
        targets.append((min(p, q), max(p, q), n))

    return targets


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 78)
    print("KAPPA EXPERIMENT 49")
    print("CUBIC-ROOT-OF-UNITY TWO-COORDINATE FACTOR EQUATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(
        f"prime interval   = "
        f"[{PRIME_LO:,}, {PRIME_HI:,})"
    )
    print(f"R values         = {R_VALUES}")

    # ------------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PRIME POPULATION")
    print("=" * 78)

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    prime_time = time.perf_counter() - t0

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.4f}s")

    if len(primes) < 2:
        raise RuntimeError("prime population is too small")

    # ------------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. TARGETS")
    print("=" * 78)

    targets = build_targets(primes)

    for i, (p, q, n) in enumerate(targets, start=1):
        print(
            f"target {i:2d}: "
            f"p={p} q={q} n={n}"
        )

    # ------------------------------------------------------------------------
    # Algebra checks
    # ------------------------------------------------------------------------

    test_cyclotomic_algebra()

    test_two_coordinate_factor_equation(targets)

    test_nontrivial_basis_elements()

    test_inverse_structure()

    compare_signature_forms(targets)

    test_norm_relation(targets)

    # ------------------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. EXPERIMENT 49 DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment deliberately does NOT ask whether the CRT filter is fast.

It asks a more fundamental algebraic question:

    Does the factor equation acquire additional structure
    when written in Z[w]/(w^2+w+1)?

The tested identities are:

    w^2 = -w - 1
    w^3 = 1

and

    (a+bw)(c+dw)
      = (ac-bd)
        + (ad+bc-bd)w.

For the actual integer factorization

    pq = n

the omega-coordinate must vanish.

The important next question is whether nontrivial
two-coordinate representations can turn the KAPPA
signature into equations that determine p or q,
rather than merely identifying them inside a prime table.

Interpretation:

    PASS:
        algebraic representation is internally consistent.

    IMPORTANT:
        a nontrivial reduction of the factor equation appears.

    NEGATIVE:
        the two-coordinate language is only a re-expression
        of ordinary modular arithmetic and gives no new constraint.

This experiment is therefore a bridge from empirical CRT filtering
to a genuinely algebraic formulation.
"""
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 49 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

