#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 73
E1 ORACLE CRT RECONSTRUCTION / CORRECTED EVEN-S DOMAIN
DIRECT CUBIC ROOT ENUMERATION VS BRUTE-FORCE SANITY
CYCLOTOMIC VS GENERIC CONTROL MODULI
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Purpose
-------
Experiment 72 established the exact oracle identity

    E1(n) = (n^2 - n + 1) * sigma_1(n) - sigma_3(n)

and, for n = p*q,

    E1(n) = (p^2 - 1)(q^2 - 1)(p + q)

with

    f_n(s) = s^3 - (n+1)^2*s + E1(n)
           = 0

at the true factor sum s = p+q.

Experiment 72 also contained a domain/scanning inconsistency:
the program claimed an EVEN s-domain but could start a step-2 scan
from an odd lower bound. This experiment removes that issue completely.

NEW METHOD
----------
Do NOT scan millions of s values.

For each modulus ell:

    f_n(s) == 0 (mod ell)

is solved directly by enumerating the roots modulo ell.

The local root sets are then combined with CRT.

We measure:

    1. number of local roots
    2. number of CRT residue classes
    3. number of CRT classes that intersect the VALID EVEN s-domain
    4. whether the true s is recovered
    5. how many moduli are needed for unique recovery
    6. cyclotomic vs generic-control comparison
    7. a small corrected brute-force sanity check

IMPORTANT
---------
E1 is still an ORACLE in this experiment.

This experiment tests the INFORMATION CONTENT of E1.

It does NOT claim E1 can currently be computed cheaply from n alone.
==============================================================================
"""
from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable


# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

SEED = 73001

NUM_TARGETS = 40

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

# Main cyclotomic families.
CYCLO = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

C3 = [7, 13, 19]
C5 = [7, 13, 19, 31, 37]
C7 = [7, 13, 19, 31, 37, 61, 67]

# Controls used previously and all satisfying ell == 1 mod 3.
CONTROLS = [
    673, 4561, 4759, 6211,
    7879, 7951, 8689, 9781
]

CONTROL_C3 = CONTROLS[:3]
CONTROL_C5 = CONTROLS[:5]
CONTROL_C7 = CONTROLS[:7]

# Small direct scan sanity check only.
SANITY_SCAN_TARGETS = 3


# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    index: int
    p: int
    q: int
    n: int
    s: int


@dataclass
class FamilyResult:
    name: str
    moduli: list[int]
    modulus_product: int

    mean_local_roots: float
    median_local_roots: float

    mean_crt_classes: float
    median_crt_classes: float

    mean_domain_hits: float
    median_domain_hits: float

    true_recovered: int
    unique_recovered: int

    mean_candidate_reduction: float
    median_candidate_reduction: float

    min_domain_hits: int
    max_domain_hits: int


# ---------------------------------------------------------------------------
# PRIME GENERATION
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> list[int]:
    """
    Return primes in [lo, hi].

    Odd-only segmented sieve keeps memory modest.
    """
    if hi < 2 or hi < lo:
        return []

    root = math.isqrt(hi)

    base = bytearray(b"\x01") * (root + 1)
    base[:2] = b"\x00\x00"

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
        if size > 0:
            flags[0] = 0
        if size > 1:
            flags[1] = 0
    elif lo == 1:
        flags[0] = 0

    for p in base_primes:
        start = max(p * p, ((lo + p - 1) // p) * p)
        for x in range(start, hi + 1, p):
            flags[x - lo] = 0

    return [
        lo + i for i, flag in enumerate(flags)
        if flag
    ]


def generate_targets(primes: list[int], count: int, seed: int) -> list[Target]:
    rng = random.Random(seed)

    out: list[Target] = []
    seen: set[tuple[int, int]] = set()

    while len(out) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)
        if key in seen:
            continue

        seen.add(key)

        n = p * q
        s = p + q

        if not (S_MIN <= s <= S_MAX):
            continue

        # The intended domain is EVEN.
        if s % 2 != 0:
            continue

        out.append(
            Target(
                index=len(out) + 1,
                p=p,
                q=q,
                n=n,
                s=s,
            )
        )

    return out


# ---------------------------------------------------------------------------
# E1 IDENTITY
# ---------------------------------------------------------------------------

def e1_oracle(t: Target) -> int:
    """
    Exact E1 value using p,q.

    This is deliberately an ORACLE.
    """
    return (t.p * t.p - 1) * (t.q * t.q - 1) * (t.p + t.q)


def cubic_value(s: int, n: int, e1: int) -> int:
    return s * s * s - (n + 1) * (n + 1) * s + e1


def validate_identity(t: Target, e1: int) -> tuple[bool, bool]:
    direct = (
        e1
        == (t.n * t.n - t.n + 1) * (t.p + 1) * (t.q + 1)
        - (t.p ** 3 + t.q ** 3 + 2 * t.n * (t.p + t.q))
    )

    cubic_zero = cubic_value(t.s, t.n, e1) == 0

    return direct, cubic_zero


# ---------------------------------------------------------------------------
# ROOT SOLVING MODULO PRIME
# ---------------------------------------------------------------------------

def cubic_roots_mod_prime(
    n: int,
    e1: int,
    ell: int,
) -> list[int]:
    """
    Directly enumerate roots of

        s^3 - (n+1)^2 s + E1 == 0 mod ell

    Because ell <= 1723 in this experiment, direct enumeration is cheap,
    deterministic, and easy to audit.
    """
    a = (n + 1) % ell
    e = e1 % ell

    roots: list[int] = []

    for r in range(ell):
        if (
            (r * r * r)
            - (a * a % ell) * r
            + e
        ) % ell == 0:
            roots.append(r)

    return roots


# ---------------------------------------------------------------------------
# CRT
# ---------------------------------------------------------------------------

def crt_pair(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> tuple[int, int]:
    """
    CRT for coprime moduli.
    Returns x in [0, m1*m2) and the combined modulus.
    """
    if math.gcd(m1, m2) != 1:
        raise ValueError("CRT moduli must be coprime")

    # x = a1 + m1*k
    # m1*k == a2-a1 (mod m2)
    inv = pow(m1, -1, m2)
    k = ((a2 - a1) * inv) % m2

    mod = m1 * m2
    x = (a1 + m1 * k) % mod

    return x, mod


def combine_residue_sets(
    current: list[int],
    current_modulus: int,
    new_roots: list[int],
    new_modulus: int,
) -> tuple[list[int], int]:
    """
    Combine a residue-class set with one new modulus.
    """
    new_values: list[int] = []

    for a in current:
        for b in new_roots:
            x, _ = crt_pair(
                a,
                current_modulus,
                b,
                new_modulus,
            )
            new_values.append(x)

    return new_values, current_modulus * new_modulus


def all_crt_classes(
    n: int,
    e1: int,
    moduli: list[int],
) -> tuple[list[int], int, list[int]]:
    """
    Return:

        CRT classes in [0,M)
        M
        local root counts

    The root sets are combined sequentially.
    """
    if not moduli:
        raise ValueError("Empty modulus family")

    classes = [0]
    M = 1
    root_counts: list[int] = []

    for ell in moduli:
        roots = cubic_roots_mod_prime(n, e1, ell)
        root_counts.append(len(roots))

        if not roots:
            return [], M * ell, root_counts

        if M == 1:
            classes = roots[:]
            M = ell
        else:
            classes, M = combine_residue_sets(
                classes,
                M,
                roots,
                ell,
            )

    return classes, M, root_counts


# ---------------------------------------------------------------------------
# EVEN DOMAIN HANDLING
# ---------------------------------------------------------------------------

def even_domain_values_for_class(
    residue: int,
    modulus: int,
    lo: int = S_MIN,
    hi: int = S_MAX,
) -> list[int]:
    """
    Find all s in [lo, hi] satisfying:

        s == residue (mod modulus)
        s even
    """
    if lo > hi:
        return []

    # First value >= lo in this residue class.
    if residue >= lo:
        s = residue
    else:
        k = (lo - residue + modulus - 1) // modulus
        s = residue + k * modulus

    # Enforce evenness.
    if s & 1:
        s += modulus

    if s > hi:
        return []

    step = modulus
    # All experimental moduli are odd, therefore adding modulus flips parity.
    # So after enforcing evenness we need 2*M for subsequent even values.
    step = 2 * modulus

    out = []
    while s <= hi:
        out.append(s)
        s += step

    return out


def domain_hits(
    classes: list[int],
    modulus: int,
    lo: int = S_MIN,
    hi: int = S_MAX,
) -> set[int]:
    """
    Set of distinct even s values in the valid domain represented by CRT
    classes.
    """
    hits: set[int] = set()

    for r in classes:
        hits.update(
            even_domain_values_for_class(
                r,
                modulus,
                lo,
                hi,
            )
        )

    return hits


# ---------------------------------------------------------------------------
# SIMPLE CORRECTED SCAN SANITY CHECK
# ---------------------------------------------------------------------------

def corrected_even_scan(
    n: int,
    e1: int,
    ell: int,
    limit: int = 500_000,
) -> set[int]:
    """
    Direct scan using an explicitly-even starting point.

    This is only a sanity checker, not the main method.
    """
    lo = math.isqrt(4 * n)

    if lo * lo < 4 * n:
        lo += 1

    if lo < S_MIN:
        lo = S_MIN

    if lo & 1:
        lo += 1

    hi = min(S_MAX, lo + 2 * limit)

    survivors = set()

    for s in range(lo, hi + 1, 2):
        if cubic_value(s, n, e1) % ell == 0:
            survivors.add(s)

    return survivors


# ---------------------------------------------------------------------------
# ANALYSIS
# ---------------------------------------------------------------------------

def analyze_family(
    targets: list[Target],
    family_name: str,
    moduli: list[int],
) -> FamilyResult:
    local_root_counts: list[int] = []
    crt_counts: list[int] = []
    domain_counts: list[int] = []
    reductions: list[float] = []

    true_recovered = 0
    unique_recovered = 0

    M = math.prod(moduli)

    for t in targets:
        e1 = e1_oracle(t)

        classes, modulus_product, root_counts = all_crt_classes(
            t.n,
            e1,
            moduli,
        )

        hits = domain_hits(
            classes,
            modulus_product,
            S_MIN,
            S_MAX,
        )

        local_root_counts.append(
            sum(root_counts) / len(root_counts)
            if root_counts
            else 0.0
        )

        crt_counts.append(len(classes))
        domain_counts.append(len(hits))

        if t.s in hits:
            true_recovered += 1

        if len(hits) == 1 and t.s in hits:
            unique_recovered += 1

        domain_size = (
            (S_MAX - S_MIN) // 2 + 1
        )

        if domain_size > 0:
            reductions.append(
                len(hits) / domain_size
            )

    return FamilyResult(
        name=family_name,
        moduli=moduli,
        modulus_product=M,

        mean_local_roots=statistics.fmean(local_root_counts),
        median_local_roots=statistics.median(local_root_counts),

        mean_crt_classes=statistics.fmean(crt_counts),
        median_crt_classes=statistics.median(crt_counts),

        mean_domain_hits=statistics.fmean(domain_counts),
        median_domain_hits=statistics.median(domain_counts),

        true_recovered=true_recovered,
        unique_recovered=unique_recovered,

        mean_candidate_reduction=statistics.fmean(reductions),
        median_candidate_reduction=statistics.median(reductions),

        min_domain_hits=min(domain_counts),
        max_domain_hits=max(domain_counts),
    )


def print_family_result(r: FamilyResult) -> None:
    domain_size = (S_MAX - S_MIN) // 2 + 1

    print(f"\n{r.name} mods={r.moduli}")
    print(f"  M = {r.modulus_product}")
    print(f"  mean local roots/mod = {r.mean_local_roots:.4f}")
    print(f"  median local roots/mod = {r.median_local_roots:.4f}")
    print(f"  mean CRT classes = {r.mean_crt_classes:.2f}")
    print(f"  median CRT classes = {r.median_crt_classes:.2f}")
    print(f"  mean even-domain hits = {r.mean_domain_hits:.2f}")
    print(f"  median even-domain hits = {r.median_domain_hits:.2f}")
    print(f"  min/max domain hits = {r.min_domain_hits}/{r.max_domain_hits}")
    print(
        f"  true recovered = "
        f"{r.true_recovered}/{NUM_TARGETS}"
    )
    print(
        f"  unique true recovery = "
        f"{r.unique_recovered}/{NUM_TARGETS}"
    )
    print(
        f"  mean domain fraction = "
        f"{r.mean_candidate_reduction:.9f}"
    )
    print(
        f"  median domain fraction = "
        f"{r.median_candidate_reduction:.9f}"
    )
    print(
        f"  raw domain size = {domain_size}"
    )


# ---------------------------------------------------------------------------
# ROOT COUNT DIAGNOSTICS
# ---------------------------------------------------------------------------

def root_count_table(
    targets: list[Target],
    moduli: list[int],
) -> dict[int, list[int]]:
    table: dict[int, list[int]] = {
        ell: [] for ell in moduli
    }

    for t in targets:
        e1 = e1_oracle(t)

        for ell in moduli:
            roots = cubic_roots_mod_prime(
                t.n,
                e1,
                ell,
            )
            table[ell].append(len(roots))

    return table


def print_root_count_table(
    targets: list[Target],
    moduli: list[int],
) -> None:
    print()
    print("LOCAL CUBIC ROOT COUNTS")
    print("-" * 78)

    table = root_count_table(targets, moduli)

    for ell in moduli:
        vals = table[ell]
        print(
            f"ell={ell:5d} "
            f"mean={statistics.fmean(vals):.4f} "
            f"min={min(vals)} "
            f"max={max(vals)} "
            f"zero={sum(v == 0 for v in vals)} "
            f"three={sum(v == 3 for v in vals)}"
        )


# ---------------------------------------------------------------------------
# MINIMAL PREFIX SEARCH
# ---------------------------------------------------------------------------

def prefix_recovery_analysis(
    targets: list[Target],
    families: list[tuple[str, list[int]]],
) -> None:
    print()
    print("PREFIX CRT RECOVERY")
    print("-" * 78)

    for name, moduli in families:
        print(f"\n{name}")

        for k in range(1, len(moduli) + 1):
            prefix = moduli[:k]

            results: list[int] = []
            unique = 0

            for t in targets:
                e1 = e1_oracle(t)

                classes, M, _ = all_crt_classes(
                    t.n,
                    e1,
                    prefix,
                )

                hits = domain_hits(
                    classes,
                    M,
                    S_MIN,
                    S_MAX,
                )

                results.append(len(hits))

                if len(hits) == 1 and t.s in hits:
                    unique += 1

            print(
                f"  first {k} moduli "
                f"M={math.prod(prefix):>12d} "
                f"mean_hits={statistics.fmean(results):>10.3f} "
                f"median={statistics.median(results):>8.1f} "
                f"unique={unique:2d}/{len(targets)}"
            )


# ---------------------------------------------------------------------------
# BRUTE-FORCE / CRT EQUIVALENCE
# ---------------------------------------------------------------------------

def sanity_check_scan_vs_crt(
    t: Target,
    family: list[int],
) -> None:
    """
    Compare corrected direct scanning against CRT for a manageable lower
    segment of the domain.
    """
    ell = family[-1]

    e1 = e1_oracle(t)

    classes, M, _ = all_crt_classes(
        t.n,
        e1,
        [ell],
    )

    crt_hits = domain_hits(
        classes,
        M,
        S_MIN,
        min(S_MAX, S_MIN + 500_000),
    )

    scan_hits = corrected_even_scan(
        t.n,
        e1,
        ell,
        limit=250_000,
    )

    # The scan and CRT ranges are slightly different because the scan begins
    # at ceil(sqrt(4n)) while the CRT sanity interval begins at S_MIN.
    # Compare only values from the intersection of the two ranges.
    shared = scan_hits & crt_hits

    print(
        f"target={t.index:3d} ell={ell:4d} "
        f"scan_hits={len(scan_hits):5d} "
        f"crt_hits_in_window={len(crt_hits):5d} "
        f"common={len(shared):5d}"
    )


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 73")
    print("E1 ORACLE CRT RECONSTRUCTION / CORRECTED EVEN-S DOMAIN")
    print("DIRECT CUBIC ROOT ENUMERATION")
    print("CYCLOTOMIC VS GENERIC CONTROL MODULI")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # -----------------------------------------------------------------------
    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    prime_time = time.perf_counter() - t0

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time = {prime_time:.6f}s")

    # -----------------------------------------------------------------------
    # 2. TARGETS
    # -----------------------------------------------------------------------
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
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    # -----------------------------------------------------------------------
    # 3. FAMILIES
    # -----------------------------------------------------------------------
    families = [
        ("E1_L3", C3),
        ("E1_L5", C5),
        ("E1_L7", C7),

        ("CONTROL_L3", CONTROL_C3),
        ("CONTROL_L5", CONTROL_C5),
        ("CONTROL_L7", CONTROL_C7),
    ]

    print()
    print("3. MODULUS FAMILIES")
    print("-" * 78)

    for name, mods in families:
        print(f"{name:<12} = {mods}")

    # -----------------------------------------------------------------------
    # 4. IDENTITY VALIDATION
    # -----------------------------------------------------------------------
    print()
    print("4. PAPER-DRIVEN E1 IDENTITY VALIDATION")
    print("-" * 78)

    failures = 0

    for t in targets:
        e1 = e1_oracle(t)
        identity_ok, cubic_zero = validate_identity(t, e1)

        if not identity_ok or not cubic_zero:
            failures += 1

        print(
            f"target {t.index:3d}: "
            f"identity={identity_ok} "
            f"cubic_zero={cubic_zero}"
        )

    print(f"identity failures = {failures}")
    print(f"status = {'PASS' if failures == 0 else 'FAIL'}")

    if failures:
        raise RuntimeError("E1 identity validation failed")

    # -----------------------------------------------------------------------
    # 5. CORRECTED DOMAIN
    # -----------------------------------------------------------------------
    print()
    print("5. CORRECTED EVEN S-DOMAIN")
    print("-" * 78)

    domain_size = (S_MAX - S_MIN) // 2 + 1

    print(
        f"s minimum = {S_MIN}"
    )
    print(
        f"s maximum = {S_MAX}"
    )
    print(
        f"even candidate count = {domain_size}"
    )
    print(
        "The scan now explicitly aligns the lower bound to an even integer."
    )

    # -----------------------------------------------------------------------
    # 6. ROOT TABLES
    # -----------------------------------------------------------------------
    print()
    print("6. LOCAL CUBIC ROOT TABLES")
    print("-" * 78)

    print_root_family = C7 + CONTROL_C7

    root_count_table(
        targets,
        print_root_family,
    )

    print_root_count_table(
        targets,
        print_root_family,
    )

    # -----------------------------------------------------------------------
    # 7. MAIN CRT ANALYSIS
    # -----------------------------------------------------------------------
    print()
    print("7. CRT CANDIDATE COMPRESSION")
    print("-" * 78)

    results: list[FamilyResult] = []

    for name, mods in families:
        r = analyze_family(
            targets,
            name,
            mods,
        )

        results.append(r)
        print_family_result(r)

    # -----------------------------------------------------------------------
    # 8. CYCLO VS CONTROL COMPARISON
    # -----------------------------------------------------------------------
    print()
    print("8. CYCLOTOMIC VS CONTROL")
    print("-" * 78)

    by_name = {r.name: r for r in results}

    for k in ("L3", "L5", "L7"):
        c = by_name[f"E1_{k}"]
        r = by_name[f"CONTROL_{k}"]

        print(
            f"{k}: "
            f"cyclo mean_hits={c.mean_domain_hits:.3f} "
            f"control mean_hits={r.mean_domain_hits:.3f} "
            f"ratio="
            f"{(
                c.mean_domain_hits / r.mean_domain_hits
                if r.mean_domain_hits
                else float('inf')
            ):.4f}"
        )

    # -----------------------------------------------------------------------
    # 9. PREFIX RECOVERY
    # -----------------------------------------------------------------------
    prefix_recovery_analysis(
        targets,
        [
            ("CYCLOTOMIC", C7),
            ("CONTROL", CONTROL_C7),
        ],
    )

    # -----------------------------------------------------------------------
    # 10. TRUE-S RECOVERY DETAILS
    # -----------------------------------------------------------------------
    print()
    print("10. FULL C7 RECOVERY BY TARGET")
    print("-" * 78)

    for family_name in ("E1_L7", "CONTROL_L7"):
        mods = C7 if family_name == "E1_L7" else CONTROL_C7

        print(f"\n{family_name}")

        for t in targets:
            e1 = e1_oracle(t)

            classes, M, roots = all_crt_classes(
                t.n,
                e1,
                mods,
            )

            hits = domain_hits(
                classes,
                M,
                S_MIN,
                S_MAX,
            )

            status = "TRUE" if t.s in hits else "MISSING"

            print(
                f"target={t.index:3d} "
                f"roots={roots} "
                f"CRT_classes={len(classes):4d} "
                f"domain_hits={len(hits):4d} "
                f"true={status}"
            )

    # -----------------------------------------------------------------------
    # 11. CORRECTED SCAN SANITY
    # -----------------------------------------------------------------------
    print()
    print("11. CORRECTED BRUTE-FORCE SANITY CHECK")
    print("-" * 78)

    for t in targets[:SANITY_SCAN_TARGETS]:
        sanity_check_scan_vs_crt(
            t,
            C7,
        )

    # -----------------------------------------------------------------------
    # 12. FINAL DIAGNOSTIC
    # -----------------------------------------------------------------------
    e1_l7 = by_name["E1_L7"]
    e1_l5 = by_name["E1_L5"]
    e1_l3 = by_name["E1_L3"]

    print()
    print("=" * 78)
    print("12. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        f"E1_L3 mean domain hits = "
        f"{e1_l3.mean_domain_hits:.3f}"
    )
    print(
        f"E1_L5 mean domain hits = "
        f"{e1_l5.mean_domain_hits:.3f}"
    )
    print(
        f"E1_L7 mean domain hits = "
        f"{e1_l7.mean_domain_hits:.3f}"
    )

    print(
        f"E1_L3 unique recovery = "
        f"{e1_l3.unique_recovered}/{NUM_TARGETS}"
    )
    print(
        f"E1_L5 unique recovery = "
        f"{e1_l5.unique_recovered}/{NUM_TARGETS}"
    )
    print(
        f"E1_L7 unique recovery = "
        f"{e1_l7.unique_recovered}/{NUM_TARGETS}"
    )

    print()
    print(
        "INTERPRETATION"
    )
    print(
        "--------------"
    )
    print(
        "1. This experiment removes the Experiment 72 parity/scanning "
        "ambiguity."
    )
    print(
        "2. CRT classes are generated directly from the cubic roots; "
        "no million-scale s scan is required."
    )
    print(
        "3. If E1_L7 consistently gives one even-domain candidate and "
        "recovers the true s, then the oracle E1 information is "
        "sufficient for exact factor-sum recovery on this domain."
    )
    print(
        "4. Control families tell us whether the compression is merely "
        "a generic low-degree modular phenomenon."
    )
    print(
        "5. The result still does NOT show how to compute E1 from n alone."
    )
    print(
        "6. The next decisive problem is therefore E1-computation "
        "without p,q."
    )

    elapsed = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 73 COMPLETE")
    print(f"total runtime = {elapsed:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

