#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 55
F(s) = s^2 + s + 1 MODULO CYClOTOMIC PRIME FACTORS
DIRECT SUM-POLYNOMIAL RESIDUE TEST
NO CSV OUTPUT
==============================================================================

Question:
    Does evaluating F(s) = s^2 + s + 1 directly on the unknown factor sum
    s = p + q reveal structure that is not explained by generic modular
    behaviour?

We compare:

    1. TRUE FACTOR SUMS
         s_true = p + q

    2. BASE SURVIVOR SUMS
         sums surviving the existing quadratic-discriminant sieve

    3. RANDOM CONTROL MODULI

For each modulus ell we measure:

    F(s) mod ell
    F(s) == 0 mod ell
    whether F(s) is a quadratic residue mod ell
    residue entropy / distinct residue count
    rank of the true sum within the candidate residue distribution

We also test the stronger condition:

    F(s) mod ell == F(r)^k mod ell

for small k, rather than assuming the true sum is itself a root.

The goal is to determine whether F(s) contains information beyond the
already-known condition

    d^2 = s^2 - 4n.

No prime-pair enumeration is required.
No CSV files are produced.
"""

from __future__ import annotations

import math
import random
import time
from collections import Counter
from statistics import mean, median

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

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

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47,
]

# Previously observed cyclotomic prime factors of F(r),
# excluding the exceptional ell=3 case.
BASE_CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723,
]

# Control primes deliberately unrelated to F(r).
CONTROL_PRIMES = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781,
]

# How many candidate sums to inspect at most per target.
# The complete even-sum range has 2,200,000 values, but limiting this
# keeps repeated diagnostics manageable.
MAX_CANDIDATE_SAMPLES = 250_000

# Small exponents for the optional cyclotomic-power comparison.
POWER_RANGE = range(0, 7)

# ---------------------------------------------------------------------------
# BASIC NUMBER THEORY
# ---------------------------------------------------------------------------


def F(x: int) -> int:
    """F(x) = x^2 + x + 1."""
    return x * x + x + 1


def factor_interval_sieve(lo: int, hi: int) -> list[int]:
    """Return all primes in [lo, hi)."""
    size = hi - lo
    is_prime = bytearray(b"\x01") * size

    if size > 0:
        is_prime[0] = 0

    limit = math.isqrt(hi - 1)

    # Small primes up to sqrt(hi).
    small = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        small[0] = 0
    if limit >= 1:
        small[1] = 0

    p = 2
    while p * p <= limit:
        if small[p]:
            start = p * p
            small[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )
        p += 1

    for p in range(2, limit + 1):
        if small[p]:
            start = max(p * p, ((lo + p - 1) // p) * p)
            if start < hi:
                is_prime[start - lo : size : p] = b"\x00" * (
                    ((hi - 1 - start) // p) + 1
                )

    return [lo + i for i, flag in enumerate(is_prime) if flag]


def is_quadratic_residue(a: int, ell: int) -> bool:
    """Legendre-symbol test for odd prime ell."""
    a %= ell
    if a == 0:
        return True
    return pow(a, (ell - 1) // 2, ell) == 1


def legendre_symbol(a: int, ell: int) -> int:
    """Return -1, 0, +1."""
    a %= ell
    if a == 0:
        return 0
    return 1 if pow(a, (ell - 1) // 2, ell) == 1 else -1


# ---------------------------------------------------------------------------
# DISCRIMINANT / SUM DOMAIN
# ---------------------------------------------------------------------------


def true_sum(p: int, q: int) -> int:
    return p + q


def even_sum_domain() -> tuple[int, int]:
    return 4_000_000, 8_399_998


def discriminant_square(n: int, s: int) -> tuple[bool, int]:
    """Check d^2 = s^2 - 4n."""
    d2 = s * s - 4 * n
    if d2 < 0:
        return False, -1

    d = math.isqrt(d2)
    return d * d == d2, d


def build_qr_tables(moduli: list[int]) -> dict[int, set[int]]:
    tables: dict[int, set[int]] = {}

    for ell in moduli:
        tables[ell] = {
            (x * x) % ell
            for x in range(ell)
        }

    return tables


def discriminant_sieve_candidates(
    n: int,
    moduli: list[int],
    qr_tables: dict[int, set[int]],
) -> list[int]:
    """
    Explicit candidate-sum sieve over the full even-sum domain.

    This is deliberately straightforward so the F(s) experiment is
    interpretable and independent of complicated indexing machinery.
    """
    s_min, s_max = even_sum_domain()

    values = range(s_min, s_max + 1, 2)

    survivors = []

    # First pass: use the smallest modulus first.
    ordered = sorted(moduli)

    for s in values:
        ok = True

        for ell in ordered:
            d2_mod = (s % ell) ** 2 - (4 * (n % ell))
            d2_mod %= ell

            if d2_mod not in qr_tables[ell]:
                ok = False
                break

        if ok:
            survivors.append(s)

        if len(survivors) > MAX_CANDIDATE_SAMPLES:
            # Keep deterministic prefix if pathological.
            break

    return survivors


# ---------------------------------------------------------------------------
# F(s) DIAGNOSTICS
# ---------------------------------------------------------------------------


def residue_profile(
    sums: list[int],
    ell: int,
) -> dict[str, object]:
    """
    Compute the complete residue profile of F(s) mod ell.
    """
    residues = [F(s) % ell for s in sums]
    counts = Counter(residues)

    zero_count = counts.get(0, 0)

    qr_count = sum(
        1
        for r in residues
        if is_quadratic_residue(r, ell)
    )

    nonzero = [r for r in residues if r != 0]

    entropy = 0.0
    total = len(residues)

    if total:
        for c in counts.values():
            p = c / total
            entropy -= p * math.log2(p)

    return {
        "count": total,
        "distinct": len(counts),
        "zero_count": zero_count,
        "zero_fraction": zero_count / total if total else 0.0,
        "qr_count": qr_count,
        "qr_fraction": qr_count / total if total else 0.0,
        "entropy": entropy,
        "counts": counts,
    }


def true_sum_statistics(
    target_s: int,
    candidate_sums: list[int],
    ell: int,
) -> dict[str, object]:
    """
    Determine where the true sum sits inside the candidate residue distribution.
    """
    target_residue = F(target_s) % ell

    profile = residue_profile(candidate_sums, ell)
    counts: Counter = profile["counts"]  # type: ignore[assignment]

    target_bucket = counts.get(target_residue, 0)

    # Rank by descending bucket frequency.
    frequencies = sorted(counts.values(), reverse=True)

    try:
        rank = frequencies.index(target_bucket) + 1
    except ValueError:
        rank = len(frequencies)

    return {
        "target_residue": target_residue,
        "target_bucket": target_bucket,
        "rank": rank,
        "zero": target_residue == 0,
        "qr": is_quadratic_residue(target_residue, ell),
        "profile": profile,
    }


def root_test(target_s: int, ell: int) -> bool:
    return F(target_s) % ell == 0


# ---------------------------------------------------------------------------
# CYClOTOMIC POWER TEST
# ---------------------------------------------------------------------------


def cyclotomic_power_values(ell: int, r_values: list[int]) -> set[int]:
    """
    Values generated by r^k mod ell for k in a small range, for source roots r.
    """
    values = set()

    for r in r_values:
        w = r % ell

        for k in POWER_RANGE:
            values.add(pow(w, k, ell))

    return values


def power_signature_test(
    target_s: int,
    ell: int,
    source_r_values: list[int],
) -> dict[str, object]:
    """
    Test whether F(s) mod ell falls into the small multiplicative
    signature set generated by the known cube-root source values.
    """
    value = F(target_s) % ell
    signature_values = cyclotomic_power_values(ell, source_r_values)

    return {
        "value": value,
        "matches_power_signature": value in signature_values,
        "signature_size": len(signature_values),
    }


# ---------------------------------------------------------------------------
# RANDOM CONTROL / EMPIRICAL NULL MODEL
# ---------------------------------------------------------------------------


def random_prime_family(
    count: int,
    lo: int = 43,
    hi: int = 12_000,
    excluded: set[int] | None = None,
) -> list[int]:
    """
    Deterministically select random odd primes from a small range.

    We generate candidates by sieve rather than relying on external libraries.
    """
    excluded = excluded or set()

    primes = [
        p
        for p in factor_interval_sieve(lo, hi)
        if p not in excluded and p > 3
    ]

    rng = random.Random(SEED + 55_001)
    rng.shuffle(primes)

    return sorted(primes[:count])


# ---------------------------------------------------------------------------
# REPORTING
# ---------------------------------------------------------------------------


def print_target_header(i: int, p: int, q: int) -> None:
    n = p * q
    s = p + q

    print()
    print("-" * 78)
    print(f"TARGET {i:2d}")
    print(f"p={p} q={q}")
    print(f"n={n}")
    print(f"s={s}")


def main() -> None:
    total_start = time.perf_counter()

    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 55")
    print("F(s) = s^2 + s + 1 MODULO CYCLOTOMIC PRIME FACTORS")
    print("DIRECT SUM-POLYNOMIAL RESIDUE TEST")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # PRIME POPULATION
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = factor_interval_sieve(PRIME_LO, PRIME_HI)
    generation_time = time.perf_counter() - t0

    prime_set = set(primes)

    print()
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {generation_time:.6f}s")

    # -----------------------------------------------------------------------
    # TARGETS
    # -----------------------------------------------------------------------

    print()
    print("2. TARGETS")
    print("-" * 78)

    for i, (p, q) in enumerate(TARGETS, 1):
        print(
            f"target {i:2d}: "
            f"p={p} q={q} n={p*q} s={p+q}"
        )

    # -----------------------------------------------------------------------
    # MODULI
    # -----------------------------------------------------------------------

    print()
    print("3. MODULUS FAMILIES")
    print("-" * 78)
    print(f"cyclotomic = {BASE_CYCLOTOMIC}")
    print(f"control    = {CONTROL_PRIMES}")

    # Check every cyclotomic modulus really comes from F(r).
    source_map: dict[int, list[int]] = {}

    for r in R_VALUES:
        m = F(r)

        for ell in BASE_CYCLOTOMIC:
            if m % ell == 0:
                source_map.setdefault(ell, []).append(r)

    print()
    print("cyclotomic source map:")
    for ell in BASE_CYCLOTOMIC:
        print(
            f"ell={ell:5d} "
            f"source_r={source_map.get(ell, [])}"
        )

    # -----------------------------------------------------------------------
    # QR TABLES
    # -----------------------------------------------------------------------

    all_moduli = BASE_CYCLOTOMIC + CONTROL_PRIMES

    t0 = time.perf_counter()
    qr_tables = build_qr_tables(all_moduli)
    qr_time = time.perf_counter() - t0

    print()
    print("4. QR TABLE PREPARATION")
    print("-" * 78)
    print(f"tables = {len(qr_tables)}")
    print(f"time   = {qr_time:.6f}s")

    # -----------------------------------------------------------------------
    # TRUE SUM DIRECT F(s) VALUES
    # -----------------------------------------------------------------------

    print()
    print("5. TRUE-SUM F(s) RESIDUES")
    print("-" * 78)

    true_values_by_target: dict[int, dict[int, int]] = {}

    for i, (p, q) in enumerate(TARGETS, 1):
        s = p + q

        values = {}

        print()
        print(f"TARGET {i:2d} s={s}")

        for ell in BASE_CYCLOTOMIC:
            v = F(s) % ell
            values[ell] = v

            print(
                f"  ell={ell:5d} "
                f"F(s) mod ell={v:5d} "
                f"zero={v == 0!s:5s} "
                f"QR={is_quadratic_residue(v, ell)!s}"
            )

        true_values_by_target[i] = values

    # -----------------------------------------------------------------------
    # ROOT TEST
    # -----------------------------------------------------------------------

    print()
    print("6. DIRECT ROOT TEST F(s) == 0")
    print("-" * 78)

    root_hits = 0
    root_tests = 0

    for i, (p, q) in enumerate(TARGETS, 1):
        s = p + q

        zero_ells = []

        for ell in BASE_CYCLOTOMIC:
            root_tests += 1
            if root_test(s, ell):
                root_hits += 1
                zero_ells.append(ell)

        print(
            f"target {i:2d}: "
            f"zero_moduli={zero_ells if zero_ells else 'none'}"
        )

    print()
    print(
        f"total zero hits = {root_hits}/{root_tests}"
    )

    # -----------------------------------------------------------------------
    # BASE DISCRIMINANT SURVIVORS
    # -----------------------------------------------------------------------

    print()
    print("7. EXISTING DISCRIMINANT BASE SURVIVORS")
    print("-" * 78)

    candidate_sets: dict[int, list[int]] = {}

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q

        t0 = time.perf_counter()
        candidates = discriminant_sieve_candidates(
            n,
            BASE_CYCLOTOMIC,
            qr_tables,
        )
        dt = time.perf_counter() - t0

        candidate_sets[i] = candidates

        print(
            f"target {i:2d}: "
            f"survivors={len(candidates):6d} "
            f"time={dt:.6f}s "
            f"true_survives={p+q in set(candidates)}"
        )

    # -----------------------------------------------------------------------
    # F(s) PROFILE ON SURVIVOR SET
    # -----------------------------------------------------------------------

    print()
    print("8. F(s) RESIDUE PROFILE ON BASE SURVIVORS")
    print("-" * 78)

    global_zero_hits = Counter()
    global_qr_advantage = Counter()

    for i, (p, q) in enumerate(TARGETS, 1):
        s_true = p + q
        candidates = candidate_sets[i]

        print()
        print(f"TARGET {i:2d} candidates={len(candidates)}")

        for ell in BASE_CYCLOTOMIC:
            stats = true_sum_statistics(
                s_true,
                candidates,
                ell,
            )

            profile = stats["profile"]

            global_zero_hits[ell] += profile["zero_count"]
            global_qr_advantage[ell] += profile["qr_count"]

            print(
                f"  ell={ell:5d} "
                f"true_F={stats['target_residue']:5d} "
                f"bucket={stats['target_bucket']:4d} "
                f"rank={stats['rank']:3d} "
                f"zero_frac={profile['zero_fraction']:.4f} "
                f"QR_frac={profile['qr_fraction']:.4f} "
                f"distinct={profile['distinct']:4d}"
            )

    # -----------------------------------------------------------------------
    # CONTROL FAMILY COMPARISON
    # -----------------------------------------------------------------------

    print()
    print("9. CYCLOTOMIC VS RANDOM-CONTROL F(s)")
    print("-" * 78)

    cyclo_zero_rates = []
    control_zero_rates = []

    cyclo_qr_rates = []
    control_qr_rates = []

    for i, (p, q) in enumerate(TARGETS, 1):
        candidates = candidate_sets[i]

        for family_name, family, zero_rates, qr_rates in [
            (
                "cyclotomic",
                BASE_CYCLOTOMIC,
                cyclo_zero_rates,
                cyclo_qr_rates,
            ),
            (
                "control",
                CONTROL_PRIMES,
                control_zero_rates,
                control_qr_rates,
            ),
        ]:
            for ell in family:
                profile = residue_profile(candidates, ell)

                zero_rates.append(profile["zero_fraction"])
                qr_rates.append(profile["qr_fraction"])

            print(
                f"target {i:2d} {family_name:12s}: "
                f"mean_zero={mean(zero_rates[-len(family):]):.6f} "
                f"mean_QR={mean(qr_rates[-len(family):]):.6f}"
            )

    print()
    print("global averages")
    print(
        f"cyclotomic zero rate = {mean(cyclo_zero_rates):.8f}"
    )
    print(
        f"control zero rate    = {mean(control_zero_rates):.8f}"
    )
    print(
        f"cyclotomic QR rate   = {mean(cyclo_qr_rates):.8f}"
    )
    print(
        f"control QR rate      = {mean(control_qr_rates):.8f}"
    )

    # -----------------------------------------------------------------------
    # POWER-SIGNATURE TEST
    # -----------------------------------------------------------------------

    print()
    print("10. CYClOTOMIC POWER-SIGNATURE TEST")
    print("-" * 78)

    total_power_hits = 0
    total_power_tests = 0

    for i, (p, q) in enumerate(TARGETS, 1):
        s = p + q

        print(f"target {i:2d}:")

        for ell in BASE_CYCLOTOMIC:
            source_r_values = source_map.get(ell, [0])

            if not source_r_values:
                continue

            result = power_signature_test(
                s,
                ell,
                source_r_values,
            )

            total_power_tests += 1
            if result["matches_power_signature"]:
                total_power_hits += 1

            print(
                f"  ell={ell:5d} "
                f"F(s)={result['value']:5d} "
                f"signature_size={result['signature_size']:3d} "
                f"match={result['matches_power_signature']}"
            )

    print()
    print(
        f"power-signature hits = "
        f"{total_power_hits}/{total_power_tests}"
    )

    # -----------------------------------------------------------------------
    # STRONGER TEST: DOES F(s) HAVE A SPECIAL VALUE?
    # -----------------------------------------------------------------------

    print()
    print("11. SPECIAL-VALUE TEST")
    print("-" * 78)

    for i, (p, q) in enumerate(TARGETS, 1):
        s = p + q

        print(f"target {i:2d}:")

        for ell in BASE_CYCLOTOMIC:
            v = F(s) % ell

            root = v == 0
            qr = is_quadratic_residue(v, ell)

            # Does F(s) itself satisfy the same cubic polynomial?
            cubic_root_test = (
                (v * v + v + 1) % ell == 0
            )

            print(
                f"  ell={ell:5d} "
                f"F={v:5d} "
                f"F=0:{root!s:5s} "
                f"QR:{qr!s:5s} "
                f"F(F)=0:{cubic_root_test!s:5s}"
            )

    # -----------------------------------------------------------------------
    # TRUE SUM VS SURVIVOR DISTRIBUTION
    # -----------------------------------------------------------------------

    print()
    print("12. TRUE-SUM POSITION TEST")
    print("-" * 78)

    true_rank_values = []

    for i, (p, q) in enumerate(TARGETS, 1):
        s_true = p + q
        candidates = candidate_sets[i]

        for ell in BASE_CYCLOTOMIC:
            stats = true_sum_statistics(
                s_true,
                candidates,
                ell,
            )

            true_rank_values.append(stats["rank"])

    print(
        f"median residue-frequency rank = "
        f"{median(true_rank_values):.2f}"
    )
    print(
        f"mean residue-frequency rank   = "
        f"{mean(true_rank_values):.2f}"
    )

    # -----------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # -----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The hypotheses being tested are:

A. NO SPECIAL STRUCTURE
   F(s) behaves like an ordinary polynomial evaluated at a
   modularly distributed sum.

B. ROOT STRUCTURE
   F(s) is unusually often 0 modulo the cyclotomic prime factors.

C. CYCLOTOMIC SIGNATURE
   F(s) modulo ell is systematically related to the root
   structure ell | F(r).

D. CANDIDATE-DISTRIBUTION EFFECT
   Among sums already surviving the discriminant sieve,
   F(s) produces a non-generic residue distribution.

The most important evidence is not one unusually small residue.
We need a repeatable difference between:

    true sums
    base survivor sums
    random/control moduli

A strong positive result would require the cyclotomic family to
show a systematic effect that is absent from the control family.

In particular:

    mean(F(s) == 0)       should differ materially
    QR(F(s)) frequency    should differ materially
    residue entropy       should differ materially
    power-signature hits  should be unusually frequent

If all of those remain close to the control family, then applying
F directly to s is mathematically interesting but does not provide
a new factorization constraint.

No prime-pair enumeration is performed.
No CSV files are produced.
"""
    )

    total_time = time.perf_counter() - total_start

    print("=" * 78)
    print("EXPERIMENT 55 COMPLETE")
    print(f"total runtime = {total_time:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

