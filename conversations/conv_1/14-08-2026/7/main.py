#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 50
CONJUGATE-ROOT INVARIANCE / TWISTED SIGNATURE INFORMATION TEST
NO CSV OUTPUT
==============================================================================

Purpose
-------
Experiment 49 established that embedding ordinary integers as (p, 0)
in Z[w]/(w^2+w+1) does not create a new constraint.

Experiment 50 moves one step further:

    1. Explicitly construct the conjugate roots
           w, w^2
       for every non-exceptional cyclotomic prime ell.

    2. Verify that the N-only equation
           p*q = n
       remains identical under both conjugate evaluations.

    3. Compare ordinary CRT compatibility against simultaneous
       conjugate-root compatibility.

    4. Construct genuinely nontrivial "twisted" signatures
           T1 = p + w*q
           T2 = p + w^2*q
       and the conjugate norm
           N(T) = T1*T2 = p^2 - p*q + q^2.

    5. Test whether these twisted signatures provide any
       N-only rejection power, or merely re-label the same
       candidate p values.

Important:
    This experiment does NOT claim a new factorization algorithm.
    It is designed to distinguish:
        - duplicated information
        - genuinely independent information
        - candidate labelling without candidate elimination.

No candidate-pair enumeration is performed.
No CSV files are produced.
"""

from __future__ import annotations

import math
import random
import time
from collections import defaultdict


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

TARGETS = 12
PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47
]

# The non-exceptional cyclotomic prime factors already established.
# These are reconstructed automatically below from F(r), but this list
# is useful as a stable expected family.
EXPECTED_ELL = [
    7, 13, 19, 31, 37, 61, 67, 79,
    127, 307, 331, 631, 1723
]


# ============================================================================
# BASIC NUMBER THEORY
# ============================================================================

def sieve_primes(low: int, high: int) -> list[int]:
    """Return all primes in [low, high)."""
    if high <= 2 or high <= low:
        return []

    limit = math.isqrt(high - 1) + 1
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    base_primes = [p for p in range(2, limit + 1) if sieve[p]]

    segment = bytearray(b"\x01") * (high - low)

    for p in base_primes:
        start = max(p * p, ((low + p - 1) // p) * p)
        for x in range(start, high, p):
            segment[x - low] = 0

    return [low + i for i, flag in enumerate(segment) if flag]


def factor_small(n: int) -> dict[int, int]:
    """Trial-factor a small positive integer."""
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


def primitive_cube_order(omega: int, ell: int) -> int:
    """Return multiplicative order of omega modulo prime ell."""
    if omega % ell == 0:
        return 0

    x = 1
    for order in range(1, ell):
        x = (x * omega) % ell
        if x == 1:
            return order

    return -1


def inverse_mod(a: int, m: int) -> int:
    """Modular inverse with explicit failure."""
    a %= m
    if math.gcd(a, m) != 1:
        raise ValueError(f"{a} is not invertible modulo {m}")
    return pow(a, -1, m)


# ============================================================================
# CYCLOTOMIC FAMILY
# ============================================================================

def cyclotomic_family() -> tuple[
    list[int],
    dict[int, list[int]],
    dict[int, list[tuple[int, int]]]
]:
    """
    Build:

        ell_values
        ell_to_r
        ell_to_roots

    where ell_to_roots[ell] contains all primitive cube roots
    seen as r mod ell.

    For a shared ell, several r values can correspond to one of the
    two conjugate roots.
    """
    ell_to_r: dict[int, list[int]] = defaultdict(list)

    for r in R_VALUES:
        m = r * r + r + 1
        factors = factor_small(m)

        for ell in factors:
            if ell == 3:
                continue
            if ell % 3 != 1:
                continue

            omega = r % ell
            if primitive_cube_order(omega, ell) == 3:
                ell_to_r[ell].append(r)

    ell_values = sorted(ell_to_r)

    ell_to_roots: dict[int, list[tuple[int, int]]] = {}

    for ell in ell_values:
        roots = []
        seen = set()

        for r in ell_to_r[ell]:
            omega = r % ell
            omega2 = (omega * omega) % ell
            pair = (omega, omega2)

            if pair not in seen:
                seen.add(pair)
                roots.append(pair)

        ell_to_roots[ell] = roots

    return ell_values, dict(ell_to_r), ell_to_roots


# ============================================================================
# TARGET GENERATION
# ============================================================================

def build_targets(primes: list[int], count: int, seed: int):
    rng = random.Random(seed)

    targets = []
    seen = set()

    while len(targets) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        a, b = sorted((p, q))
        if (a, b) in seen:
            continue

        seen.add((a, b))
        targets.append((a, b, a * b))

    return targets


# ============================================================================
# CONJUGATE ROOT OPERATIONS
# ============================================================================

def verify_root(omega: int, ell: int) -> bool:
    return (
        (omega * omega + omega + 1) % ell == 0
        and pow(omega, 3, ell) == 1
        and primitive_cube_order(omega, ell) == 3
    )


def twisted_signature(
    p: int,
    q: int,
    omega: int,
    ell: int
) -> tuple[int, int]:
    """
    T1 = p + omega*q
    T2 = p + omega^2*q

    The pair is the explicit conjugate signature.
    """
    w2 = (omega * omega) % ell

    t1 = (p + omega * q) % ell
    t2 = (p + w2 * q) % ell

    return t1, t2


def twisted_norm(
    p: int,
    q: int,
    ell: int
) -> int:
    """
    N(T) = (p + wq)(p + w^2 q)
         = p^2 - pq + q^2 mod ell.
    """
    return (p * p - p * q + q * q) % ell


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def main() -> None:
    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 50")
    print("CONJUGATE-ROOT INVARIANCE / TWISTED SIGNATURE INFORMATION TEST")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime interval   = [{PRIME_LOW:,}, {PRIME_HIGH:,})")
    print(f"R values         = {R_VALUES}")
    print()

    # ------------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("1. PRIME POPULATION")
    print("-" * 78)

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LOW, PRIME_HIGH)
    prime_time = time.perf_counter() - t0

    prime_set = set(primes)

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.4f}s")
    print()

    # ------------------------------------------------------------------------
    # 2. TARGETS
    # ------------------------------------------------------------------------
    targets = build_targets(primes, TARGETS, SEED)

    print("-" * 78)
    print("2. TARGETS")
    print("-" * 78)

    for i, (p, q, n) in enumerate(targets, 1):
        print(
            f"target {i:2d}: p={p} q={q} n={n}"
        )

    print()

    # ------------------------------------------------------------------------
    # 3. CYCLOTOMIC FAMILY
    # ------------------------------------------------------------------------
    ell_values, ell_to_r, ell_to_roots = cyclotomic_family()

    print("-" * 78)
    print("3. CYCLOTOMIC CONJUGATE-ROOT FAMILY")
    print("-" * 78)
    print(f"non-exceptional ell = {ell_values}")
    print()

    for ell in ell_values:
        roots = ell_to_roots[ell]
        root_text = ", ".join(
            f"(w={a}, w²={b})" for a, b in roots
        )
        print(
            f"ell={ell:4d} "
            f"source_r={ell_to_r[ell]} "
            f"roots={root_text}"
        )

    print()

    expected_ok = ell_values == EXPECTED_ELL
    print(f"expected family check = {'PASS' if expected_ok else 'FAIL'}")
    print()

    # ------------------------------------------------------------------------
    # 4. ROOT / CONJUGATION IDENTITY
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("4. CONJUGATE-ROOT IDENTITY CHECK")
    print("-" * 78)

    root_checks = 0
    root_failures = 0

    for ell in ell_values:
        for omega, omega2 in ell_to_roots[ell]:
            root_checks += 1

            ok1 = (omega * omega + omega + 1) % ell == 0
            ok2 = (omega2 * omega2 + omega2 + 1) % ell == 0
            ok3 = (omega * omega) % ell == omega2
            ok4 = (omega * omega2) % ell == 1

            if not (ok1 and ok2 and ok3 and ok4):
                root_failures += 1

            print(
                f"ell={ell:4d} "
                f"w={omega:4d} "
                f"w²={omega2:4d} "
                f"w³={pow(omega,3,ell):1d} "
                f"conjugate_ok={ok1 and ok2 and ok3 and ok4}"
            )

    print()
    print(f"root checks = {root_checks}")
    print(f"root failures = {root_failures}")
    print(f"status = {'PASS' if root_failures == 0 else 'FAIL'}")
    print()

    # ------------------------------------------------------------------------
    # 5. N-ONLY CONJUGATE DUPLICATION
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("5. N-ONLY EQUATION UNDER BOTH CONJUGATE ROOTS")
    print("-" * 78)
    print(
        "The scalar equation is:"
    )
    print(
        "    E = p*q - n"
    )
    print(
        "Its two root evaluations are:"
    )
    print(
        "    E(w)  = (p*q-n)*w"
    )
    print(
        "    E(w²) = (p*q-n)*w²"
    )
    print(
        "For the true factorization both are zero, so conjugation"
    )
    print(
        "does not create a second independent scalar constraint."
    )
    print()

    duplication_checks = 0
    duplication_failures = 0

    for ti, (p, q, n) in enumerate(targets, 1):
        for ell in ell_values:
            for omega, omega2 in ell_to_roots[ell]:
                duplication_checks += 1

                e = (p * q - n) % ell
                e1 = (e * omega) % ell
                e2 = (e * omega2) % ell

                if not (e == 0 and e1 == 0 and e2 == 0):
                    duplication_failures += 1

    print(f"duplication checks = {duplication_checks}")
    print(f"failures           = {duplication_failures}")
    print(
        f"status             = "
        f"{'PASS' if duplication_failures == 0 else 'FAIL'}"
    )
    print()

    # ------------------------------------------------------------------------
    # 6. FILTER EQUIVALENCE
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("6. ORDINARY VS CONJUGATE FILTER")
    print("-" * 78)

    """
    For each ell and each prime candidate p:

        q = n * p^-1 mod ell

    Ordinary condition:
        q exists as a residue.

    Conjugate condition:
        (p + w*q, p + w²*q)

    is retained.

    Since q is already derived from pq=n, the conjugate pair is always
    constructible. We therefore count whether the simultaneous
    conjugate construction rejects anything that the scalar condition
    accepts.
    """

    for ti, (p_true, q_true, n) in enumerate(targets, 1):
        print(f"TARGET {ti:2d}")

        for ell in ell_values:
            roots = ell_to_roots[ell]
            omega, omega2 = roots[0]

            ordinary = 0
            conjugate = 0
            signature_set = set()

            for p in primes:
                if p == 0 or math.gcd(p, ell) != 1:
                    continue

                q_res = (n % ell) * inverse_mod(p, ell) % ell

                # Ordinary scalar compatibility exists.
                ordinary += 1

                t1, t2 = twisted_signature(
                    p % ell, q_res, omega, ell
                )

                # Conjugate signature is always defined.
                if (
                    (t1 + t2) % ell
                    == (2 * p + (omega + omega2) * q_res) % ell
                ):
                    conjugate += 1

                signature_set.add((t1, t2))

            same = ordinary == conjugate

            print(
                f"  ell={ell:4d} "
                f"ordinary={ordinary:7d} "
                f"conjugate={conjugate:7d} "
                f"same={same}"
            )

        print()

    # ------------------------------------------------------------------------
    # 7. TWISTED SIGNATURE INFORMATION CONTENT
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("7. TWISTED SIGNATURE INFORMATION CONTENT")
    print("-" * 78)

    print(
        "For each ell we measure:"
    )
    print(
        "  - number of scalar-compatible p residues"
    )
    print(
        "  - number of distinct conjugate signatures"
    )
    print(
        "  - whether multiple p values collapse to the same signature"
    )
    print()

    for ell in ell_values:
        omega, omega2 = ell_to_roots[ell][0]

        p_to_sig: dict[int, tuple[int, int]] = {}
        sig_to_p: dict[tuple[int, int], list[int]] = defaultdict(list)

        for p in range(1, ell):
            if math.gcd(p, ell) != 1:
                continue

            # Use a fixed synthetic target n=1 here to isolate the
            # signature mapping itself:
            q = inverse_mod(p, ell)

            sig = twisted_signature(p, q, omega, ell)

            p_to_sig[p] = sig
            sig_to_p[sig].append(p)

        collisions = sum(
            1
            for values in sig_to_p.values()
            if len(values) > 1
        )

        max_bucket = max(
            (len(values) for values in sig_to_p.values()),
            default=0
        )

        print(
            f"ell={ell:4d} "
            f"units={len(p_to_sig):5d} "
            f"signatures={len(sig_to_p):5d} "
            f"colliding_signatures={collisions:4d} "
            f"max_bucket={max_bucket:2d}"
        )

    print()

    # ------------------------------------------------------------------------
    # 8. CONJUGATE NORM DIAGNOSTIC
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("8. CONJUGATE NORM DIAGNOSTIC")
    print("-" * 78)

    print(
        "For T = p + w*q:"
    )
    print(
        "    N(T) = (p+wq)(p+w²q)"
    )
    print(
        "         = p² - p*q + q²"
    )
    print(
        "         = (p+q)² - 3n"
    )
    print()
    print(
        "The final equality is important: N(T) is not a function"
    )
    print(
        "of n alone unless additional information about p+q is known."
    )
    print()

    norm_checks = 0
    norm_failures = 0

    for ti, (p, q, n) in enumerate(targets, 1):
        for ell in ell_values:
            omega, omega2 = ell_to_roots[ell][0]

            t1, t2 = twisted_signature(
                p % ell, q % ell, omega, ell
            )

            lhs = (t1 * t2) % ell
            rhs = twisted_norm(p % ell, q % ell, ell)

            direct = ((p + omega * q) * (p + omega2 * q)) % ell

            norm_checks += 1

            if not (lhs == rhs == direct):
                norm_failures += 1

        print(
            f"target {ti:2d}: "
            f"norm checks={len(ell_values):2d} "
            f"failures_so_far={norm_failures}"
        )

    print()
    print(f"norm checks = {norm_checks}")
    print(f"norm failures = {norm_failures}")
    print(
        f"status = {'PASS' if norm_failures == 0 else 'FAIL'}"
    )
    print()

    # ------------------------------------------------------------------------
    # 9. IS THE TWIST N-ONLY?
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("9. N-ONLY TEST FOR THE TWISTED SIGNATURE")
    print("-" * 78)

    print(
        "We now compare candidates having the same n mod ell."
    )
    print(
        "If T=(p+wq) were determined by n alone, all candidates"
    )
    print(
        "with the same n residue would have the same T."
    )
    print()

    for ell in ell_values:
        omega, omega2 = ell_to_roots[ell][0]

        n_res = targets[0][2] % ell

        signatures = set()
        sample_count = 0

        for p in primes:
            if math.gcd(p, ell) != 1:
                continue

            q = n_res * inverse_mod(p, ell) % ell

            t1, t2 = twisted_signature(p % ell, q, omega, ell)

            signatures.add((t1, t2))
            sample_count += 1

            if sample_count >= min(1000, len(primes)):
                break

        print(
            f"ell={ell:4d} "
            f"sample_candidates={sample_count:4d} "
            f"distinct_twisted_signatures={len(signatures):4d}"
        )

    print()

    # ------------------------------------------------------------------------
    # 10. TARGET FACTOR SIGNATURES
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("10. TRUE-FACTOR TWISTED SIGNATURES")
    print("-" * 78)

    for ti, (p, q, n) in enumerate(targets, 1):
        print(f"TARGET {ti:2d}")

        for ell in ell_values:
            omega, omega2 = ell_to_roots[ell][0]

            t1, t2 = twisted_signature(
                p % ell, q % ell, omega, ell
            )

            norm = (t1 * t2) % ell

            print(
                f"  ell={ell:4d} "
                f"T1={t1:5d} "
                f"T2={t2:5d} "
                f"norm={norm:5d}"
            )

        print()

    # ------------------------------------------------------------------------
    # 11. GLOBAL DIAGNOSTIC
    # ------------------------------------------------------------------------
    print("-" * 78)
    print("11. GLOBAL DIAGNOSTIC")
    print("-" * 78)

    print(
        "The experiment distinguishes three possibilities:"
    )
    print()
    print(
        "A. CONJUGATE DUPLICATION"
    )
    print(
        "   E = pq-n gives E(w)=E(w²)=0."
    )
    print(
        "   This adds no independent N-only constraint."
    )
    print()
    print(
        "B. TWISTED LABELLING"
    )
    print(
        "   T=(p+wq) creates a richer signature for a candidate p,"
    )
    print(
        "   but the signature depends on p and therefore does not"
    )
    print(
        "   automatically reject candidates from n alone."
    )
    print()
    print(
        "C. GENUINE NEW INFORMATION"
    )
    print(
        "   would require a twisted relation whose right-hand side"
    )
    print(
        "   is computable from n alone and which rejects candidates."
    )
    print()
    print(
        "This experiment is primarily designed to determine whether"
    )
    print(
        "case C exists for the natural conjugate-root constructions."
    )
    print()

    total_time = time.perf_counter() - total_start

    print("=" * 78)
    print("EXPERIMENT 50 COMPLETE")
    print("=" * 78)
    print(f"total runtime = {total_time:.4f}s")
    print()


if __name__ == "__main__":
    main()

