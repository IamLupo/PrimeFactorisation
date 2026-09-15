#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 73R
CORRECTED E1 ORACLE CRT RECONSTRUCTION
DIRECT CUBIC ROOT ENUMERATION
CORRECT EVEN-S DOMAIN
CYCLOTOMIC VS CONTROL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Core identity:

    E1 = (p^2 - 1)(q^2 - 1)(p + q)

with

    n = p*q
    s = p+q

Therefore

    E1 = s * ((n+1)^2 - s^2)

and hence

    f_n(s) = s^3 - (n+1)^2*s + E1

satisfies

    f_n(p+q) = 0.

IMPORTANT:
E1 is an ORACLE here.  The experiment tests what information E1 supplies.
It does NOT claim E1 is computable cheaply from n alone.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass


# ============================================================================
# CONFIG
# ============================================================================

SEED = 73073
NUM_TARGETS = 40

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

CYCLO_L3 = [7, 13, 19]
CYCLO_L5 = [7, 13, 19, 31, 37]
CYCLO_L7 = [7, 13, 19, 31, 37, 61, 67]

CONTROL_L3 = [673, 4561, 4759]
CONTROL_L5 = [673, 4561, 4759, 6211, 7879]
CONTROL_L7 = [673, 4561, 4759, 6211, 7879, 7951, 8689]


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Target:
    index: int
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> list[int]:
    if hi < 2 or lo > hi:
        return []

    root = math.isqrt(hi)

    base = bytearray(b"\x01") * (root + 1)
    base[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(root) + 1):
        if base[p]:
            start = p * p
            base[start:root + 1:p] = b"\x00" * (
                ((root - start) // p) + 1
            )

    base_primes = [
        p for p in range(2, root + 1)
        if base[p]
    ]

    size = hi - lo + 1
    flags = bytearray(b"\x01") * size

    if lo == 0:
        if size >= 1:
            flags[0] = 0
        if size >= 2:
            flags[1] = 0
    elif lo == 1:
        flags[0] = 0

    for p in base_primes:
        start = max(p * p, ((lo + p - 1) // p) * p)
        for x in range(start, hi + 1, p):
            flags[x - lo] = 0

    return [
        lo + i
        for i, flag in enumerate(flags)
        if flag
    ]


def generate_targets(
    primes: list[int],
    count: int,
    seed: int,
) -> list[Target]:
    rng = random.Random(seed)

    seen: set[tuple[int, int]] = set()
    targets: list[Target] = []

    while len(targets) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        n = p * q
        s = p + q

        # The experiment intentionally works in the EVEN s-domain.
        if s < S_MIN or s > S_MAX:
            continue

        if s % 2 != 0:
            continue

        seen.add((p, q))

        targets.append(
            Target(
                index=len(targets) + 1,
                p=p,
                q=q,
                n=n,
                s=s,
            )
        )

    return targets


# ============================================================================
# E1
# ============================================================================

def e1_oracle(t: Target) -> int:
    """
    Exact oracle value.
    """
    return (
        (t.p * t.p - 1)
        * (t.q * t.q - 1)
        * (t.p + t.q)
    )


def e1_from_ns_oracle(t: Target) -> int:
    """
    Same E1 expressed using n and the true s.

        E1 = s * ((n+1)^2 - s^2)
    """
    return (
        t.s
        * ((t.n + 1) ** 2 - t.s ** 2)
    )


def cubic_value(
    s: int,
    n: int,
    e1: int,
) -> int:
    return (
        s ** 3
        - (n + 1) ** 2 * s
        + e1
    )


# ============================================================================
# ROOTS MOD PRIME
# ============================================================================

def cubic_roots_mod_prime(
    n: int,
    e1: int,
    ell: int,
) -> list[int]:
    """
    Direct enumeration is completely safe here because all moduli
    are at most 1723.
    """
    a = (n + 1) % ell
    c = e1 % ell

    roots: list[int] = []

    for r in range(ell):
        value = (
            r ** 3
            - (a * a % ell) * r
            + c
        ) % ell

        if value == 0:
            roots.append(r)

    return roots


# ============================================================================
# CRT
# ============================================================================

def crt_pair(
    a: int,
    m: int,
    b: int,
    p: int,
) -> tuple[int, int]:
    if math.gcd(m, p) != 1:
        raise ValueError(
            f"CRT moduli are not coprime: {m}, {p}"
        )

    k = ((b - a) * pow(m, -1, p)) % p
    modulus = m * p
    x = (a + m * k) % modulus

    return x, modulus


def combine_classes(
    classes: list[int],
    modulus: int,
    roots: list[int],
    ell: int,
) -> tuple[list[int], int]:

    out: list[int] = []

    for a in classes:
        for b in roots:
            x, _ = crt_pair(
                a,
                modulus,
                b,
                ell,
            )
            out.append(x)

    return out, modulus * ell


def crt_classes(
    n: int,
    e1: int,
    moduli: list[int],
) -> tuple[list[int], int, list[int]]:

    classes = [0]
    modulus = 1
    counts: list[int] = []

    for ell in moduli:
        roots = cubic_roots_mod_prime(
            n,
            e1,
            ell,
        )

        counts.append(len(roots))

        if not roots:
            return [], modulus * ell, counts

        classes, modulus = combine_classes(
            classes,
            modulus,
            roots,
            ell,
        )

    return classes, modulus, counts


# ============================================================================
# EVEN DOMAIN
# ============================================================================

def domain_values(
    classes: list[int],
    modulus: int,
) -> set[int]:

    hits: set[int] = set()

    for r in classes:
        if r >= S_MIN:
            s = r
        else:
            k = (
                S_MIN - r + modulus - 1
            ) // modulus
            s = r + k * modulus

        # Explicitly align to even s.
        if s % 2 != 0:
            s += modulus

        step = 2 * modulus

        while s <= S_MAX:
            hits.add(s)
            s += step

    return hits


# ============================================================================
# FAMILY ANALYSIS
# ============================================================================

def analyze_family(
    targets: list[Target],
    name: str,
    moduli: list[int],
) -> None:

    domain_size = (
        (S_MAX - S_MIN) // 2 + 1
    )

    all_hits: list[int] = []
    crt_class_counts: list[int] = []
    root_means: list[float] = []

    recovered = 0
    unique = 0

    for t in targets:

        e1 = e1_oracle(t)

        classes, M, root_counts = crt_classes(
            t.n,
            e1,
            moduli,
        )

        hits = domain_values(
            classes,
            M,
        )

        all_hits.append(len(hits))
        crt_class_counts.append(len(classes))

        if root_counts:
            root_means.append(
                statistics.fmean(root_counts)
            )
        else:
            root_means.append(0.0)

        if t.s in hits:
            recovered += 1

        if len(hits) == 1 and t.s in hits:
            unique += 1

    print()
    print(name)
    print(f"  moduli              = {moduli}")
    print(f"  M                   = {math.prod(moduli)}")
    print(
        f"  mean local roots   = "
        f"{statistics.fmean(root_means):.4f}"
    )
    print(
        f"  median local roots = "
        f"{statistics.median(root_means):.4f}"
    )
    print(
        f"  mean CRT classes   = "
        f"{statistics.fmean(crt_class_counts):.3f}"
    )
    print(
        f"  median CRT classes = "
        f"{statistics.median(crt_class_counts):.3f}"
    )
    print(
        f"  mean domain hits   = "
        f"{statistics.fmean(all_hits):.3f}"
    )
    print(
        f"  median domain hits = "
        f"{statistics.median(all_hits):.3f}"
    )
    print(
        f"  min/max hits       = "
        f"{min(all_hits)}/{max(all_hits)}"
    )
    print(
        f"  true recovered     = "
        f"{recovered}/{len(targets)}"
    )
    print(
        f"  UNIQUE recovery    = "
        f"{unique}/{len(targets)}"
    )
    print(
        f"  mean domain frac   = "
        f"{statistics.fmean(all_hits) / domain_size:.12f}"
    )


# ============================================================================
# PREFIX ANALYSIS
# ============================================================================

def prefix_analysis(
    targets: list[Target],
    name: str,
    moduli: list[int],
) -> None:

    print()
    print(f"{name} PREFIX ANALYSIS")
    print("-" * 78)

    for k in range(1, len(moduli) + 1):

        prefix = moduli[:k]

        hit_counts: list[int] = []
        unique = 0
        recovered = 0

        for t in targets:

            e1 = e1_oracle(t)

            classes, M, _ = crt_classes(
                t.n,
                e1,
                prefix,
            )

            hits = domain_values(
                classes,
                M,
            )

            hit_counts.append(len(hits))

            if t.s in hits:
                recovered += 1

            if len(hits) == 1 and t.s in hits:
                unique += 1

        print(
            f"  first {k:2d}: "
            f"M={math.prod(prefix):>14d} "
            f"mean_hits={statistics.fmean(hit_counts):>10.3f} "
            f"median={statistics.median(hit_counts):>8.1f} "
            f"recovered={recovered:2d}/{len(targets)} "
            f"unique={unique:2d}/{len(targets)}"
        )


# ============================================================================
# LOCAL ROOT DIAGNOSTIC
# ============================================================================

def root_diagnostics(
    targets: list[Target],
    moduli: list[int],
) -> None:

    print()
    print("LOCAL ROOT DIAGNOSTICS")
    print("-" * 78)

    for ell in moduli:

        counts: list[int] = []

        for t in targets:

            e1 = e1_oracle(t)

            roots = cubic_roots_mod_prime(
                t.n,
                e1,
                ell,
            )

            counts.append(len(roots))

        print(
            f"ell={ell:5d} "
            f"mean={statistics.fmean(counts):.4f} "
            f"min={min(counts)} "
            f"max={max(counts)} "
            f"zero={sum(x == 0 for x in counts)} "
            f"one={sum(x == 1 for x in counts)} "
            f"two={sum(x == 2 for x in counts)} "
            f"three={sum(x == 3 for x in counts)}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 73R")
    print("CORRECTED E1 ORACLE CRT RECONSTRUCTION")
    print("DIRECT CUBIC ROOT ENUMERATION")
    print("CORRECT EVEN-S DOMAIN")
    print("CYCLOTOMIC VS CONTROL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # ----------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # ----------------------------------------------------------------------
    # 2. TARGETS
    # ----------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
        SEED,
    )

    print()
    print("2. TARGET SUMMARY")
    print("-" * 78)

    for t in targets:
        print(
            f"target {t.index:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    # ----------------------------------------------------------------------
    # 3. FAMILIES
    # ----------------------------------------------------------------------

    print()
    print("3. MODULUS FAMILIES")
    print("-" * 78)

    families = [
        ("E1_L3", CYCLO_L3),
        ("E1_L5", CYCLO_L5),
        ("E1_L7", CYCLO_L7),

        ("CONTROL_L3", CONTROL_L3),
        ("CONTROL_L5", CONTROL_L5),
        ("CONTROL_L7", CONTROL_L7),
    ]

    for name, mods in families:
        print(
            f"{name:<12} = {mods}"
        )

    # ----------------------------------------------------------------------
    # 4. CORRECT E1 VALIDATION
    # ----------------------------------------------------------------------

    print()
    print("4. E1 IDENTITY VALIDATION")
    print("-" * 78)

    failures = 0

    for t in targets:

        e1a = e1_oracle(t)
        e1b = e1_from_ns_oracle(t)

        cubic_zero = (
            cubic_value(
                t.s,
                t.n,
                e1a,
            )
            == 0
        )

        ok = (
            e1a == e1b
            and cubic_zero
        )

        if not ok:
            failures += 1

        print(
            f"target {t.index:3d}: "
            f"E1_factor={e1a == e1b} "
            f"cubic_zero={cubic_zero}"
        )

    print(
        f"E1 identity failures = {failures}"
    )
    print(
        f"status = "
        f"{'PASS' if failures == 0 else 'FAIL'}"
    )

    if failures:
        raise RuntimeError(
            "Correct E1 identity validation failed"
        )

    # ----------------------------------------------------------------------
    # 5. DOMAIN
    # ----------------------------------------------------------------------

    domain_size = (
        (S_MAX - S_MIN) // 2 + 1
    )

    print()
    print("5. EVEN S-DOMAIN")
    print("-" * 78)
    print(
        f"s minimum = {S_MIN}"
    )
    print(
        f"s maximum = {S_MAX}"
    )
    print(
        f"even candidate count = {domain_size}"
    )

    # ----------------------------------------------------------------------
    # 6. ROOT DIAGNOSTICS
    # ----------------------------------------------------------------------

    print()
    print("6. CYCLOTOMIC ROOT DIAGNOSTICS")
    print("-" * 78)

    root_diagnostics(
        targets,
        CYCLO_L7,
    )

    # ----------------------------------------------------------------------
    # 7. MAIN ANALYSIS
    # ----------------------------------------------------------------------

    print()
    print("7. CRT ANALYSIS")
    print("-" * 78)

    for name, mods in families:
        analyze_family(
            targets,
            name,
            mods,
        )

    # ----------------------------------------------------------------------
    # 8. PREFIX ANALYSIS
    # ----------------------------------------------------------------------

    prefix_analysis(
        targets,
        "CYCLOTOMIC",
        CYCLO_L7,
    )

    prefix_analysis(
        targets,
        "CONTROL",
        CONTROL_L7,
    )

    # ----------------------------------------------------------------------
    # 9. TARGET-LEVEL FULL C7 CHECK
    # ----------------------------------------------------------------------

    print()
    print("9. FULL C7 TARGET CHECK")
    print("-" * 78)

    for t in targets:

        e1 = e1_oracle(t)

        classes, M, root_counts = crt_classes(
            t.n,
            e1,
            CYCLO_L7,
        )

        hits = domain_values(
            classes,
            M,
        )

        print(
            f"target={t.index:3d} "
            f"roots={root_counts} "
            f"CRT={len(classes):4d} "
            f"domain_hits={len(hits):4d} "
            f"true_in={t.s in hits}"
        )

    # ----------------------------------------------------------------------
    # 10. FINAL DIAGNOSTIC
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The previous Experiment 73 failure was a script identity error."
    )
    print(
        "The corrected oracle identity is:"
    )
    print(
        "    E1 = (p^2 - 1)(q^2 - 1)(p + q)"
    )
    print(
        "and equivalently:"
    )
    print(
        "    E1 = s * ((n+1)^2 - s^2)"
    )
    print()
    print(
        "The cubic is:"
    )
    print(
        "    s^3 - (n+1)^2*s + E1 = 0"
    )
    print()
    print(
        "The key experimental quantity is UNIQUE recovery from the "
        "CRT classes inside the even s-domain."
    )
    print()
    print(
        "If CYCLOTOMIC_L7 repeatedly gives unique recovery while the "
        "controls do not, E1 contains extremely strong information "
        "about the factor-sum."
    )
    print()
    print(
        "The remaining hard problem is then not CRT reconstruction."
    )
    print(
        "It is obtaining E1 from N alone."
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - total_start:.6f}s"
    )
    print(
        "=" * 78
    )
    print(
        "EXPERIMENT 73R COMPLETE"
    )
    print(
        "=" * 78
    )


if __name__ == "__main__":
    main()

