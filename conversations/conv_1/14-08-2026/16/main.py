#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 58R
TRUE-SUM RESULTANT SIGNATURE VS EXACT MATCHED FALSE-SUM NULL
FULL LEGENDRE SIGNATURE TEST
NO CSV OUTPUT
==============================================================================

FIXES EXPERIMENT 58 FAILURE
----------------------------
The previous version tried to discover 240 false D-QR sums by random sampling.
For some targets there were too few hits in the random walk.

This version:

    1. Enumerates the full feasible EVEN sum domain.
    2. Applies the cyclotomic discriminant-QR condition exactly.
    3. Removes the true sum.
    4. Samples false sums DIRECTLY from that exact survivor pool.

Therefore:
    false_sums are genuinely matched D-QR survivors
    no random-search failure is possible
    no false sum is sampled twice

IMPORTANT
---------
We do NOT use the resultant as a rejection filter.

We compare:

    TRUE SUM:
        s = p + q

against

    FALSE NULL:
        arbitrary other s from the same exact cyclotomic D-QR
        survivor population for the same n.

This tests information content rather than sieve validity.
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

TARGET_COUNT = 120

# Number of false sums to sample from the exact survivor pool.
FALSE_SAMPLES_PER_TARGET = 240

# We also report the full pool size.
PRINT_POOL_SIZE = True

R_SAMPLE_MOD = 1000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

CONTROL = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

ALL_MODULI = CYCLOTOMIC + CONTROL


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


@dataclass
class SignatureStats:
    total: int = 0
    qr_r: int = 0
    nonqr_r: int = 0
    zero_r: int = 0

    def add(self, chi_r: int) -> None:
        self.total += 1

        if chi_r > 0:
            self.qr_r += 1
        elif chi_r < 0:
            self.nonqr_r += 1
        else:
            self.zero_r += 1

    def qr_fraction(self) -> float:
        nonzero = self.total - self.zero_r

        if nonzero == 0:
            return float("nan")

        return self.qr_r / nonzero


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi <= lo:
        return []

    limit = int(math.isqrt(hi - 1))

    small = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        small[0] = 0

    if limit >= 1:
        small[1] = 0

    for p in range(2, int(math.isqrt(limit)) + 1):
        if small[p]:
            start = p * p
            small[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    is_prime = bytearray(b"\x01") * (hi - lo)

    for p in range(2, limit + 1):
        if not small[p]:
            continue

        first = max(
            p * p,
            ((lo + p - 1) // p) * p
        )

        for x in range(first, hi, p):
            is_prime[x - lo] = 0

    return [
        lo + i
        for i, flag in enumerate(is_prime)
        if flag
    ]


# ============================================================================
# LEGENDRE / QR
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    v = pow(a, (p - 1) // 2, p)

    return 1 if v == 1 else -1


def quadratic_residue_table(p: int) -> bytearray:
    table = bytearray(p)

    for x in range(p):
        table[(x * x) % p] = 1

    return table


# ============================================================================
# ALGEBRA
# ============================================================================

def F(x: int) -> int:
    return x * x + x + 1


def discriminant(n: int, s: int) -> int:
    return s * s - 4 * n


def resultant(n: int, s: int) -> int:
    """
    Res_x(x^2-sx+n, x^2+x+1)

    = n^2 + n*s - n + s^2 + s + 1
    = F(p)F(q) when n=pq and s=p+q.
    """
    return (
        n * n
        + n * s
        - n
        + s * s
        + s
        + 1
    )


def verify_resultant_identity(t: Target) -> bool:
    return resultant(t.n, t.s) == F(t.p) * F(t.q)


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    rng: random.Random,
) -> List[Target]:

    out: List[Target] = []
    seen_n = set()

    attempts = 0
    max_attempts = count * 500

    while len(out) < count and attempts < max_attempts:

        attempts += 1

        i = rng.randrange(len(primes))
        j = rng.randrange(len(primes))

        if i == j:
            continue

        p = primes[min(i, j)]
        q = primes[max(i, j)]

        n = p * q

        if n in seen_n:
            continue

        out.append(
            Target(
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

        seen_n.add(n)

    if len(out) != count:
        raise RuntimeError(
            f"Could only generate {len(out)} targets "
            f"out of requested {count}."
        )

    return out


# ============================================================================
# SUM DOMAIN
# ============================================================================

def feasible_sum_bounds() -> Tuple[int, int]:
    s_min = 2 * PRIME_LO
    s_max = 2 * (PRIME_HI - 1)

    if s_min % 2:
        s_min += 1

    if s_max % 2:
        s_max -= 1

    return s_min, s_max


# ============================================================================
# EXACT MATCHED FALSE NULL
# ============================================================================

def build_exact_false_pool(
    n: int,
    true_s: int,
    s_min: int,
    s_max: int,
    qr_tables: Dict[int, bytearray],
) -> List[int]:
    """
    Build the COMPLETE population of false sums satisfying:

        D = s^2 - 4n >= 0

    and

        D is QR modulo every cyclotomic ell.

    The true sum itself is removed.

    This is the critical correction to the previous experiment.
    """

    survivors: List[int] = []

    # Since n = pq and p,q >= PRIME_LO, any candidate sum below
    # 2*sqrt(n) has D < 0. We still test explicitly.
    start = max(
        s_min,
        2 * math.isqrt(n)
    )

    # Preserve parity.
    if start % 2 != 0:
        start += 1

    for s in range(start, s_max + 1, 2):

        if s == true_s:
            continue

        d = s * s - 4 * n

        if d < 0:
            continue

        ok = True

        for ell in CYCLOTOMIC:
            if qr_tables[ell][d % ell] == 0:
                ok = False
                break

        if ok:
            survivors.append(s)

    return survivors


# ============================================================================
# SIGNATURE
# ============================================================================

def signature(
    n: int,
    s: int,
    ell: int,
) -> Tuple[int, int]:

    d = discriminant(n, s)
    r = resultant(n, s)

    return (
        legendre_symbol(d, ell),
        legendre_symbol(r, ell),
    )


# ============================================================================
# INFORMATION MEASURES
# ============================================================================

def entropy_from_counts(counts: Counter) -> float:
    total = sum(counts.values())

    if total == 0:
        return 0.0

    h = 0.0

    for c in counts.values():

        if c <= 0:
            continue

        p = c / total

        h -= p * math.log2(p)

    return h


def mutual_information(
    pairs: Sequence[Tuple[int, int]]
) -> float:

    if not pairs:
        return 0.0

    joint = Counter(pairs)

    x_counts = Counter(x for x, _ in pairs)
    y_counts = Counter(y for _, y in pairs)

    total = len(pairs)

    mi = 0.0

    for (x, y), c in joint.items():

        pxy = c / total
        px = x_counts[x] / total
        py = y_counts[y] / total

        mi += pxy * math.log2(
            pxy / (px * py)
        )

    return mi


def r_stats(
    pairs: Sequence[Tuple[int, int]]
) -> SignatureStats:

    out = SignatureStats()

    for _d, r in pairs:
        out.add(r)

    return out


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 58R")
    print("TRUE-SUM RESULTANT SIGNATURE VS EXACT MATCHED FALSE-SUM NULL")
    print("FULL LEGENDRE SIGNATURE TEST")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # ------------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI
    )

    prime_time = time.perf_counter() - t0

    print("\n" + "-" * 78)
    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.6f}s")

    # ------------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------------

    targets = generate_targets(
        primes,
        TARGET_COUNT,
        rng,
    )

    print("\n" + "-" * 78)
    print("2. GENERATED TARGETS")
    print("-" * 78)

    for i, t in enumerate(targets, 1):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    # ------------------------------------------------------------------------
    # RESULTANT IDENTITY
    # ------------------------------------------------------------------------

    failures = 0

    print("\n" + "-" * 78)
    print("3. RESULTANT IDENTITY CHECK")
    print("-" * 78)

    for i, t in enumerate(targets, 1):

        ok = verify_resultant_identity(t)

        if not ok:
            failures += 1

        print(
            f"target {i:3d}: "
            f"R(n,s)=F(p)F(q)={ok} "
            f"R mod 1000="
            f"{resultant(t.n, t.s) % 1000:03d}"
        )

    print(f"identity failures = {failures}")
    print(
        f"status = "
        f"{'PASS' if failures == 0 else 'FAIL'}"
    )

    # ------------------------------------------------------------------------
    # SUM DOMAIN
    # ------------------------------------------------------------------------

    s_min, s_max = feasible_sum_bounds()

    print("\n" + "-" * 78)
    print("4. EVEN-SUM DOMAIN")
    print("-" * 78)

    print(f"s minimum       = {s_min:,}")
    print(f"s maximum       = {s_max:,}")
    print(
        f"candidate sums  = "
        f"{(s_max - s_min) // 2 + 1:,}"
    )

    # ------------------------------------------------------------------------
    # QR TABLES
    # ------------------------------------------------------------------------

    t0 = time.perf_counter()

    qr_tables = {
        ell: quadratic_residue_table(ell)
        for ell in ALL_MODULI
    }

    qr_time = time.perf_counter() - t0

    print("\n" + "-" * 78)
    print("5. QR TABLE PREPARATION")
    print("-" * 78)

    print(f"tables = {len(qr_tables)}")
    print(f"time   = {qr_time:.6f}s")

    # ------------------------------------------------------------------------
    # ACCUMULATORS
    # ------------------------------------------------------------------------

    true_c = defaultdict(list)
    false_c = defaultdict(list)

    true_r = defaultdict(list)
    false_r = defaultdict(list)

    aggregate_true_c = Counter()
    aggregate_false_c = Counter()

    aggregate_true_r = Counter()
    aggregate_false_r = Counter()

    pool_sizes = []

    # ------------------------------------------------------------------------
    # EXACT NULL POPULATIONS
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("6. EXACT MATCHED FALSE-SUM POOLS")
    print("-" * 78)

    all_false_samples = []

    for i, t in enumerate(targets, 1):

        pool = build_exact_false_pool(
            n=t.n,
            true_s=t.s,
            s_min=s_min,
            s_max=s_max,
            qr_tables=qr_tables,
        )

        pool_sizes.append(len(pool))

        if len(pool) == 0:
            raise RuntimeError(
                f"Target {i} has no false cyclotomic "
                f"D-QR survivors."
            )

        sample_count = min(
            FALSE_SAMPLES_PER_TARGET,
            len(pool)
        )

        false_sums = rng.sample(
            pool,
            sample_count,
        )

        all_false_samples.append(
            (t, false_sums)
        )

        print(
            f"target {i:3d}: "
            f"full_false_pool={len(pool):8d} "
            f"sampled={len(false_sums):4d} "
            f"true_s_in_pool=False"
        )

        # True signature.
        for ell in CYCLOTOMIC:

            pair = signature(
                t.n,
                t.s,
                ell,
            )

            true_c[ell].append(pair)
            aggregate_true_c[pair] += 1

        for ell in CONTROL:

            pair = signature(
                t.n,
                t.s,
                ell,
            )

            true_r[ell].append(pair)
            aggregate_true_r[pair] += 1

        # False signatures.
        for sf in false_sums:

            for ell in CYCLOTOMIC:

                pair = signature(
                    t.n,
                    sf,
                    ell,
                )

                false_c[ell].append(pair)
                aggregate_false_c[pair] += 1

            for ell in CONTROL:

                pair = signature(
                    t.n,
                    sf,
                    ell,
                )

                false_r[ell].append(pair)
                aggregate_false_r[pair] += 1

    # ------------------------------------------------------------------------
    # POOL STATISTICS
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("7. FALSE-POOL STATISTICS")
    print("-" * 78)

    print(
        f"minimum pool = {min(pool_sizes):,}"
    )

    print(
        f"median pool  = "
        f"{statistics.median(pool_sizes):,.1f}"
    )

    print(
        f"maximum pool = {max(pool_sizes):,}"
    )

    print(
        f"mean pool    = "
        f"{statistics.mean(pool_sizes):,.1f}"
    )

    # ------------------------------------------------------------------------
    # SIGNATURE DISTRIBUTIONS
    # ------------------------------------------------------------------------

    signatures = [
        (-1, -1),
        (-1, 0),
        (-1, +1),
        (0, -1),
        (0, 0),
        (0, +1),
        (+1, -1),
        (+1, 0),
        (+1, +1),
    ]

    print("\n" + "-" * 78)
    print("8. SIGNATURE DISTRIBUTIONS")
    print("-" * 78)

    for name, counts in [
        ("CYCLOTOMIC TRUE", aggregate_true_c),
        ("CYCLOTOMIC FALSE", aggregate_false_c),
        ("CONTROL TRUE", aggregate_true_r),
        ("CONTROL FALSE", aggregate_false_r),
    ]:

        print(f"\n{name}")

        for sig in signatures:

            print(
                f"  {sig}: "
                f"{counts[sig]:8d}"
            )

    # ------------------------------------------------------------------------
    # RESULTANT QR FRACTIONS
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("9. RESULTANT QR FRACTIONS")
    print("-" * 78)

    def family_qr_report(
        name: str,
        true_map,
        false_map,
    ) -> None:

        print(f"\n{name}")

        trues = []
        falses = []

        for ell in sorted(true_map):

            ts = r_stats(true_map[ell])
            fs = r_stats(false_map[ell])

            tq = ts.qr_fraction()
            fq = fs.qr_fraction()

            trues.append(tq)
            falses.append(fq)

            print(
                f"  ell={ell:5d} "
                f"true={tq:.6f} "
                f"false={fq:.6f} "
                f"delta={tq-fq:+.6f}"
            )

        print(
            f"  mean true  = "
            f"{statistics.mean(trues):.6f}"
        )

        print(
            f"  mean false = "
            f"{statistics.mean(falses):.6f}"
        )

        print(
            f"  mean delta = "
            f"{statistics.mean(
                t-f for t,f in zip(trues,falses)
            ):+.6f}"
        )

    family_qr_report(
        "CYCLOTOMIC",
        true_c,
        false_c,
    )

    family_qr_report(
        "CONTROL",
        true_r,
        false_r,
    )

    # ------------------------------------------------------------------------
    # MI
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("10. CHI(D)-CHI(R) MUTUAL INFORMATION")
    print("-" * 78)

    def family_mi_report(
        name: str,
        true_map,
        false_map,
    ) -> None:

        print(f"\n{name}")

        true_mi = []
        false_mi = []

        for ell in sorted(true_map):

            a = mutual_information(
                true_map[ell]
            )

            b = mutual_information(
                false_map[ell]
            )

            true_mi.append(a)
            false_mi.append(b)

            print(
                f"  ell={ell:5d} "
                f"true={a:.6f} "
                f"false={b:.6f} "
                f"delta={a-b:+.6f}"
            )

        print(
            f"  mean true  = "
            f"{statistics.mean(true_mi):.6f}"
        )

        print(
            f"  mean false = "
            f"{statistics.mean(false_mi):.6f}"
        )

        print(
            f"  mean delta = "
            f"{statistics.mean(
                a-b for a,b in zip(true_mi,false_mi)
            ):+.6f}"
        )

    family_mi_report(
        "CYCLOTOMIC",
        true_c,
        false_c,
    )

    family_mi_report(
        "CONTROL",
        true_r,
        false_r,
    )

    # ------------------------------------------------------------------------
    # ENTROPY
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("11. SIGNATURE ENTROPY")
    print("-" * 78)

    tc = entropy_from_counts(
        aggregate_true_c
    )

    fc = entropy_from_counts(
        aggregate_false_c
    )

    tr = entropy_from_counts(
        aggregate_true_r
    )

    fr = entropy_from_counts(
        aggregate_false_r
    )

    print(
        f"cyclotomic true  = {tc:.6f} bits"
    )

    print(
        f"cyclotomic false = {fc:.6f} bits"
    )

    print(
        f"cyclotomic delta = {tc-fc:+.6f}"
    )

    print(
        f"control true     = {tr:.6f} bits"
    )

    print(
        f"control false    = {fr:.6f} bits"
    )

    print(
        f"control delta    = {tr-fr:+.6f}"
    )

    # ------------------------------------------------------------------------
    # GLOBAL R SIGN TEST
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("12. GLOBAL R-SIGN BALANCE")
    print("-" * 78)

    for name, tmap, fmap in [
        ("CYCLOTOMIC", true_c, false_c),
        ("CONTROL", true_r, false_r),
    ]:

        tp = tn = tz = 0
        fp = fn = fz = 0

        for ell in tmap:

            for _d, r in tmap[ell]:

                if r > 0:
                    tp += 1
                elif r < 0:
                    tn += 1
                else:
                    tz += 1

            for _d, r in fmap[ell]:

                if r > 0:
                    fp += 1
                elif r < 0:
                    fn += 1
                else:
                    fz += 1

        tnonzero = tp + tn
        fnonzero = fp + fn

        true_pos_frac = (
            tp / tnonzero
            if tnonzero else float("nan")
        )

        false_pos_frac = (
            fp / fnonzero
            if fnonzero else float("nan")
        )

        print(f"\n{name}")

        print(
            f"  true  +={tp} -={tn} 0={tz}"
        )

        print(
            f"  false +={fp} -={fn} 0={fz}"
        )

        print(
            f"  true positive fraction  = "
            f"{true_pos_frac:.6f}"
        )

        print(
            f"  false positive fraction = "
            f"{false_pos_frac:.6f}"
        )

        print(
            f"  difference              = "
            f"{true_pos_frac-false_pos_frac:+.6f}"
        )

    # ------------------------------------------------------------------------
    # TRUE SUM RESIDUE RANKING
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("13. TRUE-SUM R-RANK AMONG FALSE POOL")
    print("-" * 78)

    for i, (t, false_sums) in enumerate(
        all_false_samples,
        1
    ):

        ranks = []

        for ell in CYCLOTOMIC:

            true_r_mod = resultant(
                t.n,
                t.s
            ) % ell

            false_values = [
                resultant(t.n, sf) % ell
                for sf in false_sums
            ]

            less = sum(
                x < true_r_mod
                for x in false_values
            )

            rank = 1 + less

            ranks.append(
                rank / (len(false_values) + 1)
            )

        print(
            f"target {i:3d}: "
            f"median percentile="
            f"{statistics.median(ranks):.6f} "
            f"mean percentile="
            f"{statistics.mean(ranks):.6f}"
        )

    # ------------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # ------------------------------------------------------------------------

    print("\n" + "-" * 78)
    print("14. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        """
The null is now exact:

    FALSE = sums drawn directly from the COMPLETE set of
            sums satisfying the same 13 cyclotomic
            discriminant-QR conditions.

Therefore the comparison is not contaminated by random
search failure.

A positive result requires a persistent difference:

    TRUE sums
        versus
    matched false D-QR sums

and preferably:

    CYCLOTOMIC separation
        >
    CONTROL separation.

Interpretation:

    NO EFFECT
        True and false resultant signatures look alike.

    GENERIC EFFECT
        True/false separation also occurs in controls.

    CYCLOTOMIC EFFECT
        True/false separation is stronger in cyclotomic moduli
        than in the random controls.

Most important:
    This experiment does not claim that R is a valid sieve.

It asks a narrower information-theoretic question:

    does knowing that s is the actual factor sum leave a detectable
    resultant signature that is absent from other discriminant-
    compatible sums?

If the answer remains no, the resultant branch is effectively
exhausted for this line of attack.

==============================================================================
EXPERIMENT 58R COMPLETE
==============================================================================
"""
    )


if __name__ == "__main__":
    main()
