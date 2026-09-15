#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 60
CRT-COMBINED MODULUS / JOINT CONGRUENCE TEST
TRAINED D-QR NULL VS HELD-OUT CRT INFORMATION
NO CSV OUTPUT
==============================================================================

QUESTION
--------
Does multiplying several coprime moduli into

    M = ell_1 * ell_2 * ... * ell_k

provide useful information beyond the individual quadratic-residue tests?

We compare:

A. TRAIN QR
   D(s) is a quadratic residue modulo every training ell.

B. TRAIN CRT-EXACT
   D(s) == D(true_s) modulo the product of the training ell's.

C. HOLDOUT QR
   Same local QR test on held-out ell's.

D. HOLDOUT CRT-EXACT
   D(s) == D(true_s) modulo the product of held-out ell's.

E. FULL CRT-EXACT
   D(s) == D(true_s) modulo the product of all ell's.

The matched false null is constructed from sums satisfying exactly the
TRAIN QR conditions, excluding the true sum.

This is deliberately an information experiment, not a factorization
algorithm claim.

IMPORTANT:
-----------
An exact D-residue test is mathematically equivalent to knowing all the
individual residues D mod ell_i. CRT does not create new information.
The experiment asks whether the combined congruence is nevertheless a
much stronger practical discriminator on the finite scan domain.

We also test whether CRT-combination can be used to enumerate candidate
s values directly through modular square-root congruences.

NO CSV FILES.
==============================================================================

"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple


# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

SEED = 20260815

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998
S_STEP = 2

TARGET_COUNT = 24

# Keep the null manageable.
MAX_FALSE_PER_TARGET = 160

# Use a split so that the "held-out" block is genuinely unseen.
TRAIN_COUNT = 8
HOLDOUT_COUNT = 5

# Cyclotomic quadratic-residue moduli.
CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127,
    307, 331, 631, 1723
]

# Random controls with ell == 1 (mod 3), matching the arithmetic family.
CONTROL = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]


# ---------------------------------------------------------------------------
# TARGETS
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q


# First 24 targets from the established family.
TARGETS = [
    Target(3318013, 4042603),
    Target(2129167, 3402323),
    Target(2224517, 3978749),
    Target(3685051, 4020281),
    Target(2399627, 2452649),
    Target(2593039, 2996527),
    Target(2149859, 2772097),
    Target(2060543, 2514401),
    Target(2675423, 2883973),
    Target(2828887, 3960137),
    Target(3497381, 3793241),
    Target(2193509, 4011353),
    Target(3429689, 3983927),
    Target(2515871, 3050581),
    Target(2261297, 3374827),
    Target(2345537, 3673349),
    Target(2212039, 3452809),
    Target(2175373, 2476921),
    Target(2438509, 4188577),
    Target(2367553, 2399407),
    Target(2072897, 2087077),
    Target(3552023, 3707453),
    Target(3112909, 3210167),
    Target(2965819, 3239963),
]


# ---------------------------------------------------------------------------
# BASIC NUMBER THEORY
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> List[int]:
    """Return primes in [lo, hi)."""
    if hi <= 2:
        return []

    size = hi
    is_prime = bytearray(b"\x01") * size
    is_prime[:2] = b"\x00\x00"

    limit = int(math.isqrt(hi - 1))
    for p in range(2, limit + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:hi:p] = b"\x00" * (((hi - 1 - start) // p) + 1)

    return [p for p in range(lo, hi) if is_prime[p]]


def legendre_symbol(a: int, p: int) -> int:
    """
    Legendre symbol (a/p), with p an odd prime.
    Returns +1, -1, or 0.
    """
    a %= p
    if a == 0:
        return 0

    x = pow(a, (p - 1) // 2, p)
    if x == 1:
        return 1
    if x == p - 1:
        return -1
    raise RuntimeError(f"Invalid Legendre result for a={a}, p={p}")


def is_qr(a: int, p: int) -> bool:
    return legendre_symbol(a, p) >= 0


def D_of(n: int, s: int) -> int:
    return s * s - 4 * n


def F_of(x: int) -> int:
    return x * x + x + 1


# ---------------------------------------------------------------------------
# MODULAR HELPERS
# ---------------------------------------------------------------------------

def product(values: Sequence[int]) -> int:
    out = 1
    for x in values:
        out *= x
    return out


def crt_pairwise(residues: Sequence[int], moduli: Sequence[int]) -> int:
    """
    CRT for pairwise-coprime moduli.
    Returns x in [0, M).
    """
    if len(residues) != len(moduli):
        raise ValueError("CRT residue/modulus length mismatch")

    M = product(moduli)
    x = 0

    for r, m in zip(residues, moduli):
        Mi = M // m
        inv = pow(Mi, -1, m)
        x += r * Mi * inv

    return x % M


def crt_exact_D_match(
    n: int,
    s: int,
    true_s: int,
    moduli: Sequence[int],
) -> bool:
    """
    Exact combined congruence:

        D(n,s) == D(n,true_s) (mod M)

    where M = product(moduli).
    """
    M = product(moduli)
    return (D_of(n, s) - D_of(n, true_s)) % M == 0


def quadratic_residue_train(
    n: int,
    s: int,
    moduli: Sequence[int],
) -> bool:
    D = D_of(n, s)
    return all(is_qr(D, ell) for ell in moduli)


# ---------------------------------------------------------------------------
# SCAN DOMAIN
# ---------------------------------------------------------------------------

def iter_even_sums() -> Iterable[int]:
    for s in range(S_MIN, S_MAX + 1, S_STEP):
        yield s


# ---------------------------------------------------------------------------
# MATCHED FALSE NULL
# ---------------------------------------------------------------------------

def collect_false_train_qr(
    target: Target,
    train_moduli: Sequence[int],
    max_count: int,
    rng: random.Random,
) -> List[int]:
    """
    Complete training-null pool first; then sample up to max_count.

    False sums:
      * even
      * in scan domain
      * satisfy exactly the same TRAIN QR conditions
      * are not the true sum
    """
    true_s = target.s
    pool: List[int] = []

    for s in iter_even_sums():
        if s == true_s:
            continue
        if quadratic_residue_train(target.n, s, train_moduli):
            pool.append(s)

    rng.shuffle(pool)
    return pool[:max_count]


# ---------------------------------------------------------------------------
# DIRECT CRT-ROOT ENUMERATION
# ---------------------------------------------------------------------------

def square_roots_prime(a: int, p: int) -> List[int]:
    """
    Brute-force root enumeration for these relatively small ell.
    This is only used on the small prime factors of M.
    """
    a %= p
    roots = []
    for x in range(p):
        if (x * x - a) % p == 0:
            roots.append(x)
    return roots


def crt_square_roots_of_D(
    n: int,
    true_s: int,
    moduli: Sequence[int],
) -> List[int]:
    """
    Enumerate s classes modulo M satisfying

        s^2 - 4n == true_D (mod M)

    equivalently

        s^2 == true_s^2 (mod M).

    This deliberately exposes the CRT structure.

    For odd pairwise-coprime primes, each factor has at most two roots.
    """
    roots_per_modulus: List[List[int]] = []

    true_D = D_of(n, true_s)

    for ell in moduli:
        rhs = (4 * n + true_D) % ell
        roots = square_roots_prime(rhs, ell)

        if not roots:
            return []

        roots_per_modulus.append(roots)

    classes = [0]

    current_modulus = 1

    for ell, roots in zip(moduli, roots_per_modulus):
        new_classes: List[int] = []

        for base in classes:
            for r in roots:
                # Solve:
                # x = base (mod current_modulus)
                # x = r    (mod ell)
                #
                # x = base + current_modulus * k
                # k ≡ (r-base) * current_modulus^{-1} (mod ell)
                inv = pow(current_modulus, -1, ell)
                k = ((r - base) * inv) % ell
                x = base + current_modulus * k
                new_classes.append(x)

        classes = new_classes
        current_modulus *= ell

    return sorted(set(x % current_modulus for x in classes))


def count_domain_hits(classes: Sequence[int], M: int) -> List[int]:
    """
    Expand CRT classes into actual even s values in the scan domain.

    This is practical because M becomes very large after combining even
    a modest number of these moduli.
    """
    hits: List[int] = []

    for r in classes:
        # smallest k such that s = r + kM >= S_MIN
        if r < S_MIN:
            k = (S_MIN - r + M - 1) // M
        else:
            k = 0

        s = r + k * M

        # Need exact even-sum domain.
        if s % 2:
            # Since M is odd, shifting by M toggles parity.
            # Add another M if needed.
            s += M

        while s <= S_MAX:
            if s >= S_MIN and s % 2 == 0:
                hits.append(s)
            s += 2 * M

    return sorted(set(hits))


# ---------------------------------------------------------------------------
# STATISTICS
# ---------------------------------------------------------------------------

def safe_fraction(a: int, b: int) -> float:
    return a / b if b else 0.0


def entropy_from_counts(counts: Sequence[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0

    h = 0.0
    for c in counts:
        if c:
            p = c / total
            h -= p * math.log2(p)
    return h


def summarize_hit_rate(
    true_values: Sequence[int],
    false_values: Sequence[int],
) -> Tuple[float, float, float]:
    t = safe_fraction(sum(true_values), len(true_values))
    f = safe_fraction(sum(false_values), len(false_values))
    return t, f, t - f


# ---------------------------------------------------------------------------
# MAIN EXPERIMENT
# ---------------------------------------------------------------------------

def main() -> None:
    rng = random.Random(SEED)

    start_total = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 60")
    print("CRT-COMBINED MODULUS / JOINT CONGRUENCE TEST")
    print("TRAINED D-QR NULL VS HELD-OUT CRT INFORMATION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # -----------------------------------------------------------------------
    # 1. PRIME POPULATION
    # -----------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LO, PRIME_HI)
    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {time.perf_counter() - t0:.6f}s")

    # Verify targets are indeed primes.
    prime_set = set(primes)
    for i, t in enumerate(TARGETS, 1):
        if t.p not in prime_set or t.q not in prime_set:
            raise RuntimeError(f"Target {i} contains non-prime factor.")

    # -----------------------------------------------------------------------
    # 2. TARGETS
    # -----------------------------------------------------------------------

    print("\n2. TARGETS")
    print("-" * 78)

    for i, t in enumerate(TARGETS, 1):
        print(
            f"target {i:2d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    # -----------------------------------------------------------------------
    # 3. MODULUS SPLIT
    # -----------------------------------------------------------------------

    if TRAIN_COUNT + HOLDOUT_COUNT > len(CYCLOTOMIC):
        raise RuntimeError("Not enough cyclotomic moduli for requested split.")

    shuffled_cyclo = list(CYCLOTOMIC)
    rng.shuffle(shuffled_cyclo)

    train = sorted(shuffled_cyclo[:TRAIN_COUNT])
    holdout = sorted(
        shuffled_cyclo[TRAIN_COUNT:TRAIN_COUNT + HOLDOUT_COUNT]
    )

    # Remaining cyclotomic moduli are kept as an optional extra block.
    remainder = sorted(
        shuffled_cyclo[TRAIN_COUNT + HOLDOUT_COUNT:]
    )

    all_used = sorted(train + holdout)

    train_control = CONTROL[: min(TRAIN_COUNT, len(CONTROL))]
    holdout_control = CONTROL[
        min(TRAIN_COUNT, len(CONTROL)):
        min(TRAIN_COUNT + HOLDOUT_COUNT, len(CONTROL))
    ]

    M_train = product(train)
    M_hold = product(holdout)
    M_all = product(all_used)

    print("\n3. MODULUS SPLIT")
    print("-" * 78)
    print(f"cyclotomic TRAIN    = {train}")
    print(f"cyclotomic HOLDOUT  = {holdout}")
    print(f"cyclotomic REMAINDER= {remainder}")
    print(f"control TRAIN       = {train_control}")
    print(f"control HOLDOUT     = {holdout_control}")

    print("\ncombined moduli")
    print(f"M_train    = {M_train:,}")
    print(f"M_holdout  = {M_hold:,}")
    print(f"M_all      = {M_all:,}")

    print(
        f"digits: M_train={len(str(M_train))}, "
        f"M_holdout={len(str(M_hold))}, "
        f"M_all={len(str(M_all))}"
    )

    # -----------------------------------------------------------------------
    # 4. TARGET TRUE RESIDUES
    # -----------------------------------------------------------------------

    print("\n4. TRUE RESIDUE CHECK")
    print("-" * 78)

    true_records = []

    for idx, t in enumerate(TARGETS, 1):
        D_true = D_of(t.n, t.s)

        train_sig = tuple(D_true % ell for ell in train)
        hold_sig = tuple(D_true % ell for ell in holdout)

        train_crt = crt_pairwise(train_sig, train)
        hold_crt = crt_pairwise(hold_sig, holdout)

        full_sig = tuple(D_true % ell for ell in all_used)
        full_crt = crt_pairwise(full_sig, all_used)

        true_records.append(
            {
                "target": idx,
                "n": t.n,
                "s": t.s,
                "D": D_true,
                "train_crt": train_crt,
                "hold_crt": hold_crt,
                "full_crt": full_crt,
            }
        )

        print(
            f"target {idx:2d}: "
            f"s={t.s} "
            f"D mod M_train={train_crt} "
            f"D mod M_hold={hold_crt} "
            f"D mod M_all={full_crt}"
        )

    # -----------------------------------------------------------------------
    # 5. MATCHED FALSE POOLS
    # -----------------------------------------------------------------------

    print("\n5. MATCHED FALSE D-QR POOLS")
    print("-" * 78)
    print(
        "False sums satisfy EXACTLY the TRAIN cyclotomic D-QR conditions."
    )

    false_pools: List[List[int]] = []

    pool_t0 = time.perf_counter()

    for idx, t in enumerate(TARGETS, 1):
        false_pool = collect_false_train_qr(
            t,
            train,
            MAX_FALSE_PER_TARGET,
            rng,
        )

        if not false_pool:
            raise RuntimeError(
                f"Target {idx}: no false TRAIN-D-QR sums found."
            )

        false_pools.append(false_pool)

        print(
            f"target {idx:2d}: "
            f"false_count={len(false_pool):4d} "
            f"true_s_in_pool={t.s in false_pool}"
        )

    print(f"pool construction time = {time.perf_counter() - pool_t0:.6f}s")

    # -----------------------------------------------------------------------
    # 6. INDIVIDUAL QR VS CRT-EXACT
    # -----------------------------------------------------------------------

    print("\n6. TRUE VS MATCHED-FALSE FILTER RATES")
    print("-" * 78)

    # Store per-target results for global summaries.
    per_target = []

    for idx, (t, false_pool) in enumerate(zip(TARGETS, false_pools), 1):
        true_s = t.s

        # These are always true for the true sum.
        train_true_qr = quadratic_residue_train(t.n, true_s, train)
        train_true_exact = crt_exact_D_match(
            t.n,
            true_s,
            true_s,
            train,
        )

        hold_true_exact = crt_exact_D_match(
            t.n,
            true_s,
            true_s,
            holdout,
        )

        full_true_exact = crt_exact_D_match(
            t.n,
            true_s,
            true_s,
            all_used,
        )

        train_false_qr = [
            1
            for s in false_pool
            if quadratic_residue_train(t.n, s, train)
        ]

        train_false_exact = [
            1
            for s in false_pool
            if crt_exact_D_match(t.n, s, true_s, train)
        ]

        hold_false_qr = [
            1
            for s in false_pool
            if quadratic_residue_train(t.n, s, holdout)
        ]

        hold_false_exact = [
            1
            for s in false_pool
            if crt_exact_D_match(t.n, s, true_s, holdout)
        ]

        full_false_exact = [
            1
            for s in false_pool
            if crt_exact_D_match(t.n, s, true_s, all_used)
        ]

        # Cross-check by direct CRT arithmetic.
        train_exact_direct = 0
        hold_exact_direct = 0

        train_true_sig = tuple(D_of(t.n, true_s) % ell for ell in train)
        hold_true_sig = tuple(D_of(t.n, true_s) % ell for ell in holdout)

        train_true_crt = crt_pairwise(train_true_sig, train)
        hold_true_crt = crt_pairwise(hold_true_sig, holdout)

        for s in false_pool:
            D = D_of(t.n, s)

            train_sig = tuple(D % ell for ell in train)
            hold_sig = tuple(D % ell for ell in holdout)

            train_crt = crt_pairwise(train_sig, train)
            hold_crt = crt_pairwise(hold_sig, holdout)

            if train_crt == train_true_crt:
                train_exact_direct += 1

            if hold_crt == hold_true_crt:
                hold_exact_direct += 1

        if train_exact_direct != sum(train_false_exact):
            raise RuntimeError("Train CRT cross-check failed.")

        if hold_exact_direct != sum(hold_false_exact):
            raise RuntimeError("Holdout CRT cross-check failed.")

        train_exact_rate = safe_fraction(
            sum(train_false_exact),
            len(false_pool),
        )
        hold_exact_rate = safe_fraction(
            sum(hold_false_exact),
            len(false_pool),
        )
        full_exact_rate = safe_fraction(
            sum(full_false_exact),
            len(false_pool),
        )

        hold_qr_rate = safe_fraction(
            sum(hold_false_qr),
            len(false_pool),
        )

        per_target.append(
            {
                "train_false_count": len(false_pool),
                "train_exact_rate": train_exact_rate,
                "hold_qr_rate": hold_qr_rate,
                "hold_exact_rate": hold_exact_rate,
                "full_exact_rate": full_exact_rate,
            }
        )

        print(f"\nTARGET {idx:2d} true_s={true_s}")
        print(
            f"  train true QR       = {train_true_qr}"
        )
        print(
            f"  train false QR      = {len(train_false_qr)}/{len(false_pool)} "
            f"= {safe_fraction(len(train_false_qr), len(false_pool)):.6f}"
        )
        print(
            f"  train false CRT     = {len(train_false_exact)}/{len(false_pool)} "
            f"= {train_exact_rate:.6f}"
        )
        print(
            f"  holdout false QR    = {len(hold_false_qr)}/{len(false_pool)} "
            f"= {hold_qr_rate:.6f}"
        )
        print(
            f"  holdout false CRT   = {len(hold_false_exact)}/{len(false_pool)} "
            f"= {hold_exact_rate:.6f}"
        )
        print(
            f"  full false CRT      = {len(full_false_exact)}/{len(false_pool)} "
            f"= {full_exact_rate:.6f}"
        )
        print(
            f"  exact CRT true      = "
            f"train={train_true_exact} "
            f"holdout={hold_true_exact} "
            f"full={full_true_exact}"
        )

    # -----------------------------------------------------------------------
    # 7. GLOBAL FALSE-NULL RATES
    # -----------------------------------------------------------------------

    print("\n7. GLOBAL NULL RATES")
    print("-" * 78)

    total_false = sum(len(x) for x in false_pools)

    train_exact_total = sum(
        round(r["train_exact_rate"] * r["train_false_count"])
        for r in per_target
    )

    hold_qr_total = 0
    hold_exact_total = 0
    full_exact_total = 0

    for idx, (t, pool) in enumerate(zip(TARGETS, false_pools), 1):
        for s in pool:
            hold_qr_total += int(
                quadratic_residue_train(t.n, s, holdout)
            )
            hold_exact_total += int(
                crt_exact_D_match(t.n, s, t.s, holdout)
            )
            full_exact_total += int(
                crt_exact_match if False else
                crt_exact_D_match(t.n, s, t.s, all_used)
            )

    print(
        f"total matched false sums = {total_false:,}"
    )
    print(
        f"TRAIN exact-CRT false rate = "
        f"{train_exact_total / total_false:.8f}"
    )
    print(
        f"HOLDOUT QR false rate      = "
        f"{hold_qr_total / total_false:.8f}"
    )
    print(
        f"HOLDOUT exact-CRT rate     = "
        f"{hold_exact_total / total_false:.8f}"
    )
    print(
        f"FULL exact-CRT rate        = "
        f"{full_exact_total / total_false:.8f}"
    )

    # -----------------------------------------------------------------------
    # 8. DIRECT DOMAIN SCAN
    # -----------------------------------------------------------------------

    print("\n8. CRT CLASS ENUMERATION")
    print("-" * 78)
    print(
        "For each target, enumerate the classes satisfying "
        "s^2 == true_s^2 (mod M)."
    )

    crt_domain_stats = []

    for idx, t in enumerate(TARGETS, 1):
        t0 = time.perf_counter()

        train_classes = crt_square_roots_of_D(
            t.n,
            t.s,
            train,
        )

        hold_classes = crt_square_roots_of_D(
            t.n,
            t.s,
            holdout,
        )

        full_classes = crt_square_roots_of_D(
            t.n,
            t.s,
            all_used,
        )

        train_hits = count_domain_hits(train_classes, M_train)
        hold_hits = count_domain_hits(hold_classes, M_hold)
        full_hits = count_domain_hits(full_classes, M_all)

        # The true sum must be recovered by the full CRT classes.
        full_true_recovered = t.s in full_hits

        crt_domain_stats.append(
            {
                "train_classes": len(train_classes),
                "hold_classes": len(hold_classes),
                "full_classes": len(full_classes),
                "train_hits": len(train_hits),
                "hold_hits": len(hold_hits),
                "full_hits": len(full_hits),
                "true_recovered": full_true_recovered,
            }
        )

        print(f"\nTARGET {idx:2d}")
        print(
            f"  train CRT classes = {len(train_classes):4d} "
            f"domain hits = {len(train_hits):4d}"
        )
        print(
            f"  holdout CRT classes = {len(hold_classes):4d} "
            f"domain hits = {len(hold_hits):4d}"
        )
        print(
            f"  full CRT classes = {len(full_classes):4d} "
            f"domain hits = {len(full_hits):4d}"
        )
        print(
            f"  true s recovered = {full_true_recovered}"
        )

        if not full_true_recovered:
            raise RuntimeError(
                f"Target {idx}: full CRT enumeration lost true sum."
            )

        print(
            f"  enumeration time = {time.perf_counter() - t0:.6f}s"
        )

    # -----------------------------------------------------------------------
    # 9. CRITICAL INFORMATION CHECK
    # -----------------------------------------------------------------------

    print("\n9. CRT INFORMATION EQUIVALENCE CHECK")
    print("-" * 78)

    print(
        "For each test sum we verify:"
    )
    print(
        "  all local residues equal"
        "  <=>  CRT-combined residue equal."
    )

    equivalence_failures = 0

    for idx, (t, pool) in enumerate(zip(TARGETS, false_pools), 1):
        true_D = D_of(t.n, t.s)

        for s in pool:
            D = D_of(t.n, s)

            local_same = all(
                (D - true_D) % ell == 0
                for ell in all_used
            )

            combined_same = (
                (D - true_D) % M_all == 0
            )

            if local_same != combined_same:
                equivalence_failures += 1

    print(f"equivalence failures = {equivalence_failures}")

    if equivalence_failures != 0:
        raise RuntimeError("CRT equivalence check failed.")

    print("status = PASS")

    # -----------------------------------------------------------------------
    # 10. PRACTICAL ADVANTAGE TEST
    # -----------------------------------------------------------------------

    print("\n10. PRACTICAL SCANNING TEST")
    print("-" * 78)

    print(
        "The question here is whether a giant M permits a scan shortcut."
    )

    for idx, t in enumerate(TARGETS, 1):
        full_stats = crt_domain_stats[idx - 1]

        domain_size = (
            ((S_MAX - S_MIN) // S_STEP) + 1
        )

        reduction = safe_fraction(
            full_stats["full_hits"],
            domain_size,
        )

        print(
            f"target {idx:2d}: "
            f"domain={domain_size:,} "
            f"full_CRT_hits={full_stats['full_hits']:,} "
            f"fraction={reduction:.12f}"
        )

    # -----------------------------------------------------------------------
    # 11. SUMMARY
    # -----------------------------------------------------------------------

    print("\n11. FINAL SUMMARY")
    print("-" * 78)

    med_train_classes = statistics.median(
        x["train_classes"] for x in crt_domain_stats
    )
    med_hold_classes = statistics.median(
        x["hold_classes"] for x in crt_domain_stats
    )
    med_full_classes = statistics.median(
        x["full_classes"] for x in crt_domain_stats
    )

    med_train_hits = statistics.median(
        x["train_hits"] for x in crt_domain_stats
    )
    med_hold_hits = statistics.median(
        x["hold_hits"] for x in crt_domain_stats
    )
    med_full_hits = statistics.median(
        x["full_hits"] for x in crt_domain_stats
    )

    print(f"median train CRT classes = {med_train_classes}")
    print(f"median holdout CRT classes = {med_hold_classes}")
    print(f"median full CRT classes = {med_full_classes}")

    print(f"median train domain hits = {med_train_hits}")
    print(f"median holdout domain hits = {med_hold_hits}")
    print(f"median full domain hits = {med_full_hits}")

    print("\nInterpretation:")
    print(
        "  1. CRT does NOT create new mathematical information beyond"
        "     the complete set of local congruences."
    )
    print(
        "  2. If full CRT classes are extremely sparse, it can still be"
        "     computationally useful as a compact representation of a"
        "     joint congruence condition."
    )
    print(
        "  3. If the full CRT class set is essentially identical to what"
        "     the individual residues already imply, there is no new"
        "     mathematical phenomenon."
    )
    print(
        "  4. A genuine advance would require a new N-only way to generate"
        "     or restrict these CRT classes without knowing true_s."
    )
    print(
        "  5. The most important follow-up is therefore to replace the"
        "     unknown true residue D(true_s) with information computable"
        "     from n alone, then test whether the CRT modulus yields a"
        "     stronger-than-local scan."
    )

    print("\n" + "=" * 78)
    print("EXPERIMENT 60 COMPLETE")
    print("=" * 78)

    print(
        f"total runtime = {time.perf_counter() - start_total:.6f}s"
    )


if __name__ == "__main__":
    main()

