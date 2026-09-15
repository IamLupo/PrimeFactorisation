#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 72
MACMAHON PRIME-DETECTOR CUBIC RESOLVENT
ORACLE E1 MODULAR CANDIDATE COMPRESSION TEST
COMPARE E1-CUBIC VS D-QR VS INTERSECTION
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Research question
-----------------
Craig–van Ittersum–Ono derive the first prime-detecting expression

    E1(n) = (n^2 - n + 1) * sigma_1(n) - sigma_3(n)

because their MacMahon expression
    (n^2 - 3n + 2) M1(n) - 8 M2(n)
reduces to it.

For a semiprime n = p*q, this experiment derives the exact identity

    E1(n) = (p^2 - 1)(q^2 - 1)(p + q)
          = s * ((n + 1)^2 - s^2),

where s = p + q.

Therefore the true s satisfies

    s^3 - (n + 1)^2 s + E1(n) = 0.

This experiment deliberately computes E1 using the known factors as an ORACLE.
It does NOT claim E1 is efficiently computable from n alone.

The point is to measure the theoretical information content of the
MacMahon/Craig resolvent:
    1. How many s-candidates survive E1 mod each modulus?
    2. How many survive combined E1 congruences?
    3. How does that compare with D = s^2 - 4n quadratic-residue filters?
    4. How much does the intersection improve the candidate space?

A strong result here is evidence that this resolvent is structurally useful.
A practical result requires a separate way to compute E1 (or its residues)
without already knowing p and q.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 72072
TARGETS = 120

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

# Use prime moduli that were central in the earlier Kappa experiments.
MODULI = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723,
    673, 4561, 4759, 6211, 7879, 7951, 8689, 9781,
]

# Prefixes let us measure information growth.
E1_PREFIXES = {
    "E1_L3": [7, 13, 19],
    "E1_L5": [7, 13, 19, 31, 37],
    "E1_L7": [7, 13, 19, 31, 37, 61, 67],
}

DQR_PREFIXES = {
    "DQR_L3": [7, 13, 19],
    "DQR_L5": [7, 13, 19, 31, 37],
    "DQR_L7": [7, 13, 19, 31, 37, 61, 67],
}


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(math.isqrt(n)) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * (((n - start) // p) + 1)
    return [i for i, flag in enumerate(sieve) if flag]


def random_semiprime_targets(rng: random.Random, primes: list[int]) -> list[tuple[int, int, int, int]]:
    """
    Return (p, q, n, s). Distinct primes are used.
    """
    targets: list[tuple[int, int, int, int]] = []
    seen: set[int] = set()

    while len(targets) < TARGETS:
        p = rng.choice(primes)
        q = rng.choice(primes)
        if p == q:
            continue
        n = p * q
        if n in seen:
            continue
        seen.add(n)
        s = p + q
        if S_MIN <= s <= S_MAX:
            targets.append((p, q, n, s))

    return targets


# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

def sigma1_semiprime(p: int, q: int) -> int:
    return (p + 1) * (q + 1)


def sigma3_semiprime(p: int, q: int) -> int:
    return (p ** 3 + 1) * (q ** 3 + 1)


def macmahon_E1(n: int, sigma1: int, sigma3: int) -> int:
    """
    Craig et al.'s first prime-detector after eliminating M2:
        E1 = (n^2 - n + 1) sigma_1(n) - sigma_3(n)
    """
    return (n * n - n + 1) * sigma1 - sigma3


def direct_E1_identity(p: int, q: int) -> int:
    s = p + q
    n = p * q
    return s * ((n + 1) * (n + 1) - s * s)


def cubic_value(n: int, s: int, e1: int) -> int:
    return s ** 3 - (n + 1) ** 2 * s + e1


def legendre_symbol(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    v = pow(a, (p - 1) // 2, p)
    if v == 1:
        return 1
    if v == p - 1:
        return -1
    raise RuntimeError(f"unexpected Legendre result a={a}, p={p}, v={v}")


# ---------------------------------------------------------------------------
# Candidate-domain helpers
# ---------------------------------------------------------------------------

def domain_bounds(n: int) -> tuple[int, int]:
    lo = max(S_MIN, math.isqrt(4 * n))
    while lo * lo < 4 * n:
        lo += 1

    hi = min(S_MAX, int(math.isqrt(4 * n)) + 2_000_000)
    if hi < lo:
        raise RuntimeError(f"empty s-domain for n={n}")
    return lo, hi


def qr_d_filter(n: int, s: int, ell: int) -> bool:
    """
    Basic discriminant QR condition.
    For D = s^2 - 4n, require D to be a quadratic residue mod ell.
    Zero counts as admissible.
    """
    d = (s * s - 4 * n) % ell
    return d == 0 or legendre_symbol(d, ell) == 1


def e1_cubic_filter(n: int, e1_mod: int, s: int, ell: int) -> bool:
    """
    E1 cubic congruence:
        s^3 - (n+1)^2 s + E1 == 0 (mod ell)
    """
    return (
        (pow(s % ell, 3, ell)
         - ((n + 1) % ell) ** 2 % ell * (s % ell)
         + e1_mod)
        % ell
        == 0
    )


def combined_e1_filter(n: int, e1: int, s: int, moduli: list[int]) -> bool:
    return all(e1_cubic_filter(n, e1 % ell, s, ell) for ell in moduli)


def combined_dqr_filter(n: int, s: int, moduli: list[int]) -> bool:
    return all(qr_d_filter(n, s, ell) for ell in moduli)


# ---------------------------------------------------------------------------
# One-target experiment
# ---------------------------------------------------------------------------

@dataclass
class PrefixResult:
    name: str
    e1_candidates: int
    dqr_candidates: int
    intersection_candidates: int
    true_pass_e1: bool
    true_pass_dqr: bool
    true_pass_both: bool


def scan_prefix(
    n: int,
    s_true: int,
    e1: int,
    name: str,
    moduli: list[int],
    lo: int,
    hi: int,
) -> PrefixResult:
    e1_candidates = 0
    dqr_candidates = 0
    intersection_candidates = 0

    for s in range(lo, hi + 1, 2):
        # All target s are even in our generated population.
        e1_ok = combined_e1_filter(n, e1, s, moduli)
        dqr_ok = combined_dqr_filter(n, s, moduli)

        if e1_ok:
            e1_candidates += 1
        if dqr_ok:
            dqr_candidates += 1
        if e1_ok and dqr_ok:
            intersection_candidates += 1

    return PrefixResult(
        name=name,
        e1_candidates=e1_candidates,
        dqr_candidates=dqr_candidates,
        intersection_candidates=intersection_candidates,
        true_pass_e1=combined_e1_filter(n, e1, s_true, moduli),
        true_pass_dqr=combined_dqr_filter(n, s_true, moduli),
        true_pass_both=combined_e1_filter(n, e1, s_true, moduli)
        and combined_dqr_filter(n, s_true, moduli),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 72")
    print("MACMAHON PRIME-DETECTOR CUBIC RESOLVENT")
    print("ORACLE E1 MODULAR CANDIDATE COMPRESSION TEST")
    print("COMPARE E1-CUBIC VS D-QR VS INTERSECTION")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    t0 = time.perf_counter()
    primes = primes_upto(PRIME_MAX)
    primes = [p for p in primes if p >= PRIME_MIN]
    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time = {time.perf_counter() - t0:.6f}s")

    targets = random_semiprime_targets(rng, primes)

    print("\n2. TARGET SUMMARY")
    print("-" * 78)
    for i, (p, q, n, s) in enumerate(targets, 1):
        print(f"target {i:3d}: p={p} q={q} n={n} s={s}")

    print("\n3. PAPER-DRIVEN RESOLVENT IDENTITY")
    print("-" * 78)
    print("E1(n) = (n^2 - n + 1) sigma_1(n) - sigma_3(n)")
    print("For n=pq:")
    print("E1(n) = (p^2-1)(q^2-1)(p+q)")
    print("      = s * ((n+1)^2 - s^2)")
    print("Hence:")
    print("s^3 - (n+1)^2 s + E1(n) = 0")

    print("\n4. IDENTITY VALIDATION")
    print("-" * 78)
    identity_failures = 0

    e1_values: list[int] = []
    for i, (p, q, n, s) in enumerate(targets, 1):
        s1 = sigma1_semiprime(p, q)
        s3 = sigma3_semiprime(p, q)
        e1_a = macmahon_E1(n, s1, s3)
        e1_b = direct_E1_identity(p, q)
        cubic = cubic_value(n, s, e1_a)

        ok = (e1_a == e1_b) and (cubic == 0)
        if not ok:
            identity_failures += 1

        e1_values.append(e1_a)

        if i <= 10:
            print(
                f"target {i:3d}: E1_identity={e1_a == e1_b} "
                f"cubic_zero={cubic == 0} "
                f"log10(E1)={math.log10(e1_a):.4f}"
            )

    print(f"identity failures = {identity_failures}")
    print("status =", "PASS" if identity_failures == 0 else "FAIL")

    print("\n5. CANDIDATE-DOMAIN PARAMETERS")
    print("-" * 78)
    print(f"s domain = [{S_MIN}, {S_MAX}] even values")
    print("Only the positive even s-domain is scanned.")
    print("E1 residues are ORACLE values derived from the known p,q pair.")

    print("\n6. MODULAR INFORMATION TEST")
    print("-" * 78)

    all_prefixes = list(E1_PREFIXES.items())

    aggregate: dict[str, list[PrefixResult]] = {name: [] for name, _ in all_prefixes}

    scan_start = time.perf_counter()

    for idx, ((p, q, n, s), e1) in enumerate(zip(targets, e1_values), 1):
        lo, hi = domain_bounds(n)

        # Full scans are intentionally limited to first 40 targets to keep
        # runtime bounded. Algebraic identities are validated for all targets.
        if idx > 40:
            break

        for name, mods in all_prefixes:
            result = scan_prefix(n, s, e1, name, mods, lo, hi)
            aggregate[name].append(result)

        if idx % 5 == 0:
            print(
                f"processed {idx:3d}/40 targets "
                f"time={time.perf_counter() - scan_start:.3f}s"
            )

    print("\n7. AGGREGATE CANDIDATE COMPRESSION")
    print("-" * 78)

    for name, results in aggregate.items():
        if not results:
            continue

        e1_counts = [r.e1_candidates for r in results]
        dqr_counts = [r.dqr_candidates for r in results]
        both_counts = [r.intersection_candidates for r in results]

        e1_true = sum(r.true_pass_e1 for r in results)
        dqr_true = sum(r.true_pass_dqr for r in results)
        both_true = sum(r.true_pass_both for r in results)

        print(f"\n{name} mods={E1_PREFIXES[name]}")
        print(
            f"  mean E1 candidates       = {statistics.fmean(e1_counts):.2f}"
        )
        print(
            f"  median E1 candidates     = {statistics.median(e1_counts):.2f}"
        )
        print(
            f"  mean D-QR candidates     = {statistics.fmean(dqr_counts):.2f}"
        )
        print(
            f"  median D-QR candidates   = {statistics.median(dqr_counts):.2f}"
        )
        print(
            f"  mean intersection        = {statistics.fmean(both_counts):.2f}"
        )
        print(
            f"  median intersection      = {statistics.median(both_counts):.2f}"
        )
        print(f"  true passes E1           = {e1_true}/{len(results)}")
        print(f"  true passes D-QR         = {dqr_true}/{len(results)}")
        print(f"  true passes intersection = {both_true}/{len(results)}")

    print("\n8. INFORMATION GAIN RATIOS")
    print("-" * 78)

    domain_sizes = []
    for (p, q, n, s) in targets[:40]:
        lo, hi = domain_bounds(n)
        domain_sizes.append(((hi - lo) // 2) + 1)

    domain_mean = statistics.fmean(domain_sizes)
    print(f"mean raw even-domain size = {domain_mean:.2f}")

    for name, results in aggregate.items():
        if not results:
            continue
        mean_e1 = statistics.fmean(r.e1_candidates for r in results)
        mean_dqr = statistics.fmean(r.dqr_candidates for r in results)
        mean_both = statistics.fmean(r.intersection_candidates for r in results)

        print(f"{name}")
        print(
            f"  E1 survivor fraction       = {mean_e1 / domain_mean:.8f}"
        )
        print(
            f"  D-QR survivor fraction     = {mean_dqr / domain_mean:.8f}"
        )
        print(
            f"  intersection fraction     = {mean_both / domain_mean:.8f}"
        )
        if mean_both > 0:
            print(
                f"  E1->intersection retention = {mean_both / mean_e1:.6f}"
                if mean_e1 > 0 else
                "  E1->intersection retention = 0.000000"
            )

    print("\n9. LOCAL CUBIC-ROOT COUNT")
    print("-" * 78)
    for name, mods in all_prefixes:
        counts = []
        for ell in mods:
            c = 0
            for (p, q, n, s), e1 in zip(targets[:40], e1_values[:40]):
                lo = s % ell
                roots = [
                    r for r in range(ell)
                    if (r ** 3 - ((n + 1) % ell) ** 2 * r + e1) % ell == 0
                ]
                c += len(roots)
            counts.append(c / min(40, len(targets)))
        print(
            f"{name}: mean cubic roots/modulus = "
            + ", ".join(
                f"{ell}:{count:.3f}" for ell, count in zip(mods, counts)
            )
        )

    print("\n10. WHAT THIS EXPERIMENT CAN / CANNOT SHOW")
    print("-" * 78)
    print("CAN:")
    print("  * establish an exact semiprime cubic resolvent from the paper identity")
    print("  * quantify its modular candidate-compression potential")
    print("  * compare E1 information against D-QR information")
    print("  * test whether the two constraints intersect unusually strongly")
    print("CANNOT:")
    print("  * claim an N-only algorithm, because E1 is supplied as an ORACLE")
    print("  * prove that sigma_1 or sigma_3 is cheap to evaluate on unknown semiprimes")
    print("  * claim factorization from the resolvent without a practical E1 computation")

    print("\n11. FINAL DIAGNOSTIC")
    print("-" * 78)
    print(
        "The key research question is now computationally precise:\n"
        "  If E1 residues were cheaply available from n, would the induced\n"
        "  cubic congruences collapse the p+q search space substantially?\n\n"
        "A strong next result would require a second experiment that computes\n"
        "the same E1 residues WITHOUT using p and q, using the paper's\n"
        "MacMahon/q-series machinery or a provably cheaper recurrence."
    )

    print("\n===============================================================================")
    print("EXPERIMENT 72 COMPLETE")
    print(f"total runtime = {time.perf_counter() - t0:.6f}s")
    print("===============================================================================")


if __name__ == "__main__":
    main()
