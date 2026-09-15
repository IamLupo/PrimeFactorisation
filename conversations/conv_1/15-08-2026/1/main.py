#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 60R
OUT-OF-SAMPLE CUBIC-CHARACTER HOLDOUT TEST
TRUE FACTOR SUM VS EXACT MATCHED FALSE D-QR NULL
PATCHED COMMON-POOL SAMPLING
NO CSV OUTPUT
==============================================================================

PATCH
-----
The previous version required every target to contain 240 false D-QR sums.

That assumption is false: some targets have smaller complete false pools.

This version:

  1. constructs the COMPLETE false D-QR pool for every target;
  2. records every pool size;
  3. chooses one common matched sample size

         COMMON_FALSE_COUNT =
             min(FALSE_CAP, minimum complete pool size)

  4. samples exactly COMMON_FALSE_COUNT false sums from EVERY target.

Therefore every target contributes the same number of false observations,
and no target is discarded merely because its null pool is smaller than 240.

The holdout cubic-character moduli are still completely excluded from
construction of the false pools.

No prime-pair enumeration is performed.
No CSV files are produced.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260815

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

CHUNK_SIZE = 200_000

TARGET_COUNT = 60

# Maximum desired false observations per target.
FALSE_CAP = 240

# Out-of-sample holdout families.
HOLDOUT_PRIME_MAX = 100_000
HOLDOUT_CYCLO_COUNT = 200
HOLDOUT_CONTROL_COUNT = 200

# New source-r window. This is disjoint from the small training source
# values used in the earlier experiments.
HOLDOUT_R_LO = 53
HOLDOUT_R_HI = 2_000


# ============================================================================
# TRAINING CYCLOTOMIC FAMILY
# ============================================================================

BASE_CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]


# ============================================================================
# TARGETS
# ============================================================================

TARGETS_FIXED = [
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

    (3429689, 3983927),
    (2515871, 3050581),
    (2261297, 3374827),
    (2345537, 3673349),
    (2212039, 3452809),
    (2175373, 2476921),
    (2438509, 4188577),
    (2367553, 2399407),
    (2072897, 2087077),
    (3552023, 3707453),
    (3112909, 3210167),
    (2965819, 3239963),
    (2053903, 2270987),
    (2857333, 3816073),
    (2258339, 2296909),
    (2096599, 2331377),
    (3549901, 3807889),
    (3667751, 4092703),
    (2350687, 3923893),
    (2689007, 3164827),
    (3214879, 3580091),
    (2082061, 3078997),
    (2870279, 3340847),
    (3394891, 3714769),
    (2644981, 3942373),
    (2285779, 2866837),
    (2873657, 3778679),
    (2809567, 3367097),
    (2170813, 3133399),
    (2294993, 4169213),
    (2804327, 3136657),
    (2774143, 3605869),
    (2147009, 3495413),
    (2484827, 3048511),
    (2094361, 2480909),
    (2455589, 3779437),
    (2095109, 2528489),
    (2372761, 2983441),
    (3045323, 3605293),
    (3412159, 3774709),
    (2497951, 2730811),
    (2179097, 2382389),
    (3592637, 3863173),
    (2294771, 2496409),
    (2415443, 4095331),
    (2549123, 3646637),
    (2065897, 2486041),
    (2881271, 4165229),
]


# ============================================================================
# DATA TYPES
# ============================================================================

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


@dataclass(frozen=True)
class CubicModulus:
    ell: int
    zeta: int
    source_r: int | None


# ============================================================================
# PRIME / QR UTILITIES
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi <= 2:
        return []

    sieve = np.ones(hi, dtype=np.bool_)
    sieve[:2] = False

    limit = int(math.isqrt(hi - 1))

    for p in range(2, limit + 1):
        if sieve[p]:
            sieve[p * p:hi:p] = False

    arr = np.flatnonzero(sieve)
    arr = arr[arr >= lo]
    return [int(x) for x in arr]


def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    x = pow(a, (p - 1) // 2, p)

    if x == 1:
        return 1

    return -1


def qr_table(p: int) -> np.ndarray:
    table = np.zeros(p, dtype=np.bool_)
    x = np.arange(p, dtype=np.int64)
    table[(x * x) % p] = True
    return table


# ============================================================================
# CUBIC CHARACTER
# ============================================================================

def primitive_cubic_root(ell: int) -> int:
    if ell % 3 != 1:
        raise ValueError(f"ell={ell} is not 1 mod 3")

    exponent = (ell - 1) // 3

    for g in range(2, min(100, ell)):
        z = pow(g, exponent, ell)
        if z != 1 and pow(z, 3, ell) == 1:
            return z

    for g in range(2, ell):
        z = pow(g, exponent, ell)
        if z != 1 and pow(z, 3, ell) == 1:
            return z

    raise RuntimeError(f"No nontrivial cubic root found for ell={ell}")


def cubic_roots_of_F(ell: int, zeta: int) -> Tuple[int, int]:
    r1 = zeta
    r2 = pow(zeta, 2, ell)

    if (r1 * r1 + r1 + 1) % ell != 0:
        raise AssertionError(f"Bad root r1 for ell={ell}")

    if (r2 * r2 + r2 + 1) % ell != 0:
        raise AssertionError(f"Bad root r2 for ell={ell}")

    return r1, r2


def build_holdout_families(
    used: set[int],
) -> Tuple[List[CubicModulus], List[CubicModulus]]:

    primes = sieve_primes(2, HOLDOUT_PRIME_MAX)

    cyclo: List[CubicModulus] = []
    control_candidates: List[CubicModulus] = []

    holdout_seen: set[int] = set()

    for ell in primes:
        if ell % 3 != 1:
            continue

        if ell in used:
            continue

        zeta = primitive_cubic_root(ell)
        r1, r2 = cubic_roots_of_F(ell, zeta)

        source_r = None

        for r in (r1, r2):
            if HOLDOUT_R_LO <= r <= HOLDOUT_R_HI:
                source_r = r
                break

        cm = CubicModulus(
            ell=ell,
            zeta=zeta,
            source_r=source_r,
        )

        if source_r is not None and len(cyclo) < HOLDOUT_CYCLO_COUNT:
            cyclo.append(cm)
            holdout_seen.add(ell)
        else:
            control_candidates.append(cm)

        if (
            len(cyclo) >= HOLDOUT_CYCLO_COUNT
            and len(control_candidates) >= HOLDOUT_CONTROL_COUNT
        ):
            break

    if len(cyclo) < HOLDOUT_CYCLO_COUNT:
        raise RuntimeError(
            f"Only found {len(cyclo)} cyclotomic holdout moduli."
        )

    rng = random.Random(SEED + 6000)
    rng.shuffle(control_candidates)

    control: List[CubicModulus] = []

    for cm in control_candidates:
        if cm.ell in holdout_seen:
            continue

        control.append(cm)

        if len(control) >= HOLDOUT_CONTROL_COUNT:
            break

    if len(control) < HOLDOUT_CONTROL_COUNT:
        raise RuntimeError(
            f"Only found {len(control)} control holdout moduli."
        )

    return cyclo, control


def cubic_ratio_class(
    a: int,
    b: int,
    cm: CubicModulus,
) -> int:

    ell = cm.ell

    a %= ell
    b %= ell

    if a == 0 or b == 0:
        raise ValueError(
            f"Zero denominator/numerator in cubic ratio mod {ell}"
        )

    rho = a * pow(b, -1, ell) % ell

    value = pow(rho, (ell - 1) // 3, ell)

    zeta = cm.zeta
    zeta2 = pow(zeta, 2, ell)

    if value == 1:
        return 0

    if value == zeta:
        return 1

    if value == zeta2:
        return 2

    raise AssertionError(
        f"Unexpected cubic value ell={ell}, rho={rho}, value={value}"
    )


# ============================================================================
# MODULAR SQUARE ROOT
# ============================================================================

def mod_sqrt(a: int, p: int) -> int | None:
    a %= p

    if a == 0:
        return 0

    if p == 2:
        return a

    if pow(a, (p - 1) // 2, p) != 1:
        return None

    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    q = p - 1
    s = 0

    while q % 2 == 0:
        q //= 2
        s += 1

    z = 2

    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1

    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    r = pow(a, (q + 1) // 2, p)

    while t != 1:

        i = 1
        t2 = t * t % p

        while i < m and t2 != 1:
            t2 = t2 * t2 % p
            i += 1

        if i == m:
            return None

        b = pow(c, 1 << (m - i - 1), p)

        r = r * b % p
        c = b * b % p
        t = t * c % p
        m = i

    return r


# ============================================================================
# TRAINING DISCRIMINANT FILTER
# ============================================================================

def discriminant_qr_survives(
    s: int,
    n: int,
    moduli: Sequence[int],
    tables: Dict[int, np.ndarray],
) -> bool:

    d = s * s - 4 * n

    for ell in moduli:
        if not tables[ell][d % ell]:
            return False

    return True


# ============================================================================
# IMPORTANT PATCH:
# CONSTRUCT THE COMPLETE FALSE POOL WITHOUT REQUIRING 240
# ============================================================================

def build_complete_false_pool(
    target: Target,
    moduli: Sequence[int],
    tables: Dict[int, np.ndarray],
) -> List[int]:

    """
    Return ALL false sums in the feasible even-sum domain satisfying
    the original 13 cyclotomic D-QR constraints.

    The true sum is excluded.

    This function deliberately has NO "needed" parameter.
    That is the core patch.
    """

    true_s = target.s
    out: List[int] = []

    for block_start in range(
        S_MIN,
        S_MAX + 1,
        2 * CHUNK_SIZE,
    ):
        block_end = min(
            S_MAX,
            block_start + 2 * CHUNK_SIZE - 2,
        )

        sums = np.arange(
            block_start,
            block_end + 1,
            2,
            dtype=np.int64,
        )

        d = sums * sums - np.int64(4) * np.int64(target.n)

        mask = np.ones(
            sums.shape,
            dtype=np.bool_,
        )

        for ell in moduli:

            residues = np.mod(
                d,
                ell,
            )

            mask &= tables[ell][residues]

            if not np.any(mask):
                break

        if np.any(mask):

            vals = sums[mask]

            for value in vals:

                s = int(value)

                if s != true_s:
                    out.append(s)

    return out


# ============================================================================
# FORMAL ROOT-RATIO CLASS FOR A CANDIDATE SUM
# ============================================================================

def formal_root_ratio_class(
    s: int,
    n: int,
    cm: CubicModulus,
) -> int | None:

    ell = cm.ell

    D = (s * s - 4 * n) % ell

    if D == 0:
        return 0

    if legendre_symbol(D, ell) != 1:
        return None

    root_D = mod_sqrt(D, ell)

    if root_D is None:
        return None

    inv2 = pow(2, -1, ell)

    x1 = (s + root_D) * inv2 % ell
    x2 = (s - root_D) * inv2 % ell

    if x1 == 0 or x2 == 0:
        return None

    return cubic_ratio_class(
        x1,
        x2,
        cm,
    )


# ============================================================================
# STATISTICS
# ============================================================================

def class_counts(values: Iterable[int]) -> List[int]:
    counts = [0, 0, 0]

    for c in values:
        counts[c] += 1

    return counts


def entropy(counts: Sequence[int]) -> float:

    total = sum(counts)

    if total == 0:
        return 0.0

    h = 0.0

    for c in counts:

        if c == 0:
            continue

        p = c / total
        h -= p * math.log2(p)

    return h


def nontrivial_fraction(values: Sequence[int]) -> float:

    if not values:
        return 0.0

    return sum(
        c in (1, 2)
        for c in values
    ) / len(values)


# ============================================================================
# EVALUATE ONE TARGET
# ============================================================================

def evaluate_target(
    target: Target,
    false_pool: Sequence[int],
    moduli: Sequence[CubicModulus],
) -> Tuple[List[int], List[int], int, int]:

    true_classes: List[int] = []
    false_classes: List[int] = []

    true_skips = 0
    false_skips = 0

    for cm in moduli:

        c_true = formal_root_ratio_class(
            target.s,
            target.n,
            cm,
        )

        if c_true is None:
            true_skips += 1
        else:
            true_classes.append(c_true)

        for s in false_pool:

            c_false = formal_root_ratio_class(
                s,
                target.n,
                cm,
            )

            if c_false is None:
                false_skips += 1
            else:
                false_classes.append(c_false)

    return (
        true_classes,
        false_classes,
        true_skips,
        false_skips,
    )


# ============================================================================
# SWAP TEST
# ============================================================================

def swap_test(
    targets: Sequence[Target],
    moduli: Sequence[CubicModulus],
) -> int:

    failures = 0

    for target in targets:

        for cm in moduli:

            c1 = cubic_ratio_class(
                target.p,
                target.q,
                cm,
            )

            c2 = cubic_ratio_class(
                target.q,
                target.p,
                cm,
            )

            if c2 != (-c1) % 3:
                failures += 1

    return failures


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    overall_start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 60R")
    print("OUT-OF-SAMPLE CUBIC-CHARACTER HOLDOUT TEST")
    print("PATCHED COMMON-POOL FALSE NULL")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1. TARGETS
    # ----------------------------------------------------------------------

    targets = [
        Target(p, q)
        for p, q in TARGETS_FIXED[:TARGET_COUNT]
    ]

    print()
    print("1. TARGETS")
    print("-" * 78)

    for i, t in enumerate(targets, 1):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    # ----------------------------------------------------------------------
    # 2. TRAINING MODULI
    # ----------------------------------------------------------------------

    print()
    print("2. TRAINING CYCLOTOMIC FAMILY")
    print("-" * 78)
    print(BASE_CYCLOTOMIC)

    tables: Dict[int, np.ndarray] = {}

    t0 = time.perf_counter()

    for ell in BASE_CYCLOTOMIC:
        tables[ell] = qr_table(ell)

    print(
        f"QR table preparation = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # ----------------------------------------------------------------------
    # 3. HOLDOUT MODULI
    # ----------------------------------------------------------------------

    print()
    print("3. OUT-OF-SAMPLE HOLDOUT FAMILIES")
    print("-" * 78)

    t0 = time.perf_counter()

    cyclo_holdout, control_holdout = build_holdout_families(
        used=set(BASE_CYCLOTOMIC),
    )

    print(
        f"cyclotomic holdout = "
        f"{len(cyclo_holdout)}"
    )

    print(
        f"control holdout    = "
        f"{len(control_holdout)}"
    )

    print(
        f"construction time  = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    print()
    print("first cyclotomic holdouts:")
    print(
        [
            (cm.ell, cm.source_r)
            for cm in cyclo_holdout[:20]
        ]
    )

    print()
    print("first control holdouts:")
    print(
        [
            cm.ell
            for cm in control_holdout[:20]
        ]
    )

    # ----------------------------------------------------------------------
    # 4. BUILD COMPLETE FALSE POOLS
    # ----------------------------------------------------------------------

    print()
    print("4. COMPLETE FALSE D-QR POOLS")
    print("-" * 78)

    pool_start = time.perf_counter()

    complete_pools: List[List[int]] = []

    for i, target in enumerate(targets, 1):

        pool = build_complete_false_pool(
            target=target,
            moduli=BASE_CYCLOTOMIC,
            tables=tables,
        )

        if target.s in pool:
            raise AssertionError(
                f"True sum entered false pool for target {i}"
            )

        complete_pools.append(pool)

        print(
            f"target {i:3d}: "
            f"complete_false_pool={len(pool):4d}"
        )

    pool_time = time.perf_counter() - pool_start

    print()
    print(
        f"pool construction time = "
        f"{pool_time:.6f}s"
    )

    # ----------------------------------------------------------------------
    # 5. PATCHED COMMON SAMPLE SIZE
    # ----------------------------------------------------------------------

    min_pool = min(
        len(pool)
        for pool in complete_pools
    )

    common_false_count = min(
        FALSE_CAP,
        min_pool,
    )

    if common_false_count < 20:
        raise RuntimeError(
            "Common false-pool size is too small "
            f"({common_false_count}); "
            "experiment should not proceed."
        )

    print()
    print("5. COMMON MATCHED FALSE SAMPLE")
    print("-" * 78)

    print(
        f"requested cap        = {FALSE_CAP}"
    )

    print(
        f"minimum complete pool = {min_pool}"
    )

    print(
        f"COMMON_FALSE_COUNT    = "
        f"{common_false_count}"
    )

    print(
        "All targets will contribute exactly "
        f"{common_false_count} false sums."
    )

    false_pools: List[List[int]] = []

    for i, pool in enumerate(complete_pools, 1):

        if len(pool) == common_false_count:

            selected = list(pool)

        else:

            selected = rng.sample(
                pool,
                common_false_count,
            )

        if target.s if False else False:
            pass

        false_pools.append(selected)

        print(
            f"target {i:3d}: "
            f"complete={len(pool):4d} "
            f"sampled={len(selected):4d}"
        )

    # ----------------------------------------------------------------------
    # 6. BASIC NULL VALIDATION
    # ----------------------------------------------------------------------

    print()
    print("6. NULL VALIDATION")
    print("-" * 78)

    validation_failures = 0

    for i, (target, pool) in enumerate(
        zip(targets, false_pools),
        1,
    ):

        for s in pool:

            if s == target.s:
                validation_failures += 1
                continue

            if not discriminant_qr_survives(
                s,
                target.n,
                BASE_CYCLOTOMIC,
                tables,
            ):
                validation_failures += 1

    print(
        f"training-null validation failures = "
        f"{validation_failures}"
    )

    if validation_failures:
        raise RuntimeError(
            "False-pool validation failed."
        )

    print("status = PASS")

    # ----------------------------------------------------------------------
    # 7. HOLDOUT CUBIC TEST
    # ----------------------------------------------------------------------

    print()
    print("7. HOLDOUT CUBIC-CHARACTER TEST")
    print("-" * 78)

    global_true_c: List[int] = []
    global_false_c: List[int] = []

    global_true_r: List[int] = []
    global_false_r: List[int] = []

    cyclo_deltas: List[float] = []
    control_deltas: List[float] = []

    cyclo_h_deltas: List[float] = []
    control_h_deltas: List[float] = []

    for i, (target, pool) in enumerate(
        zip(targets, false_pools),
        1,
    ):

        (
            true_c,
            false_c,
            true_skip_c,
            false_skip_c,
        ) = evaluate_target(
            target,
            pool,
            cyclo_holdout,
        )

        (
            true_r,
            false_r,
            true_skip_r,
            false_skip_r,
        ) = evaluate_target(
            target,
            pool,
            control_holdout,
        )

        global_true_c.extend(true_c)
        global_false_c.extend(false_c)

        global_true_r.extend(true_r)
        global_false_r.extend(false_r)

        c_true_non = nontrivial_fraction(true_c)
        c_false_non = nontrivial_fraction(false_c)

        r_true_non = nontrivial_fraction(true_r)
        r_false_non = nontrivial_fraction(false_r)

        c_delta = c_true_non - c_false_non
        r_delta = r_true_non - r_false_non

        c_counts_true = class_counts(true_c)
        c_counts_false = class_counts(false_c)

        r_counts_true = class_counts(true_r)
        r_counts_false = class_counts(false_r)

        c_hd = (
            entropy(c_counts_true)
            - entropy(c_counts_false)
        )

        r_hd = (
            entropy(r_counts_true)
            - entropy(r_counts_false)
        )

        cyclo_deltas.append(c_delta)
        control_deltas.append(r_delta)

        cyclo_h_deltas.append(c_hd)
        control_h_deltas.append(r_hd)

        print(
            f"target {i:3d}: "
            f"C true={c_true_non:.5f} "
            f"false={c_false_non:.5f} "
            f"delta={c_delta:+.5f} "
            f"skip=({true_skip_c},{false_skip_c}) | "
            f"R true={r_true_non:.5f} "
            f"false={r_false_non:.5f} "
            f"delta={r_delta:+.5f} "
            f"skip=({true_skip_r},{false_skip_r})"
        )

    # ----------------------------------------------------------------------
    # 8. GLOBAL RESULTS
    # ----------------------------------------------------------------------

    print()
    print("8. GLOBAL HOLDOUT DISTRIBUTIONS")
    print("-" * 78)

    ctc = class_counts(global_true_c)
    cfc = class_counts(global_false_c)

    rtc = class_counts(global_true_r)
    rfc = class_counts(global_false_r)

    c_true_non = nontrivial_fraction(global_true_c)
    c_false_non = nontrivial_fraction(global_false_c)

    r_true_non = nontrivial_fraction(global_true_r)
    r_false_non = nontrivial_fraction(global_false_r)

    print("CYCLOTOMIC")
    print(
        f"  TRUE  classes={ctc} "
        f"entropy={entropy(ctc):.6f}"
    )
    print(
        f"  FALSE classes={cfc} "
        f"entropy={entropy(cfc):.6f}"
    )
    print(
        f"  TRUE nontrivial  = {c_true_non:.8f}"
    )
    print(
        f"  FALSE nontrivial = {c_false_non:.8f}"
    )
    print(
        f"  delta            = "
        f"{c_true_non - c_false_non:+.8f}"
    )

    print()
    print("CONTROL")
    print(
        f"  TRUE  classes={rtc} "
        f"entropy={entropy(rtc):.6f}"
    )
    print(
        f"  FALSE classes={rfc} "
        f"entropy={entropy(rfc):.6f}"
    )
    print(
        f"  TRUE nontrivial  = {r_true_non:.8f}"
    )
    print(
        f"  FALSE nontrivial = {r_false_non:.8f}"
    )
    print(
        f"  delta            = "
        f"{r_true_non - r_false_non:+.8f}"
    )

    # ----------------------------------------------------------------------
    # 9. TARGET-LEVEL SUMMARY
    # ----------------------------------------------------------------------

    print()
    print("9. TARGET-LEVEL SEPARATION SUMMARY")
    print("-" * 78)

    cyclo_mean = statistics.mean(cyclo_deltas)
    cyclo_median = statistics.median(cyclo_deltas)

    control_mean = statistics.mean(control_deltas)
    control_median = statistics.median(control_deltas)

    cyclo_entropy_mean = statistics.mean(cyclo_h_deltas)
    control_entropy_mean = statistics.mean(control_h_deltas)

    print(
        f"cyclotomic mean delta   = "
        f"{cyclo_mean:+.8f}"
    )

    print(
        f"cyclotomic median delta = "
        f"{cyclo_median:+.8f}"
    )

    print(
        f"control mean delta      = "
        f"{control_mean:+.8f}"
    )

    print(
        f"control median delta    = "
        f"{control_median:+.8f}"
    )

    print(
        f"cyclo-control mean gap  = "
        f"{cyclo_mean - control_mean:+.8f}"
    )

    print(
        f"cyclo-control median gap = "
        f"{cyclo_median - control_median:+.8f}"
    )

    print()
    print(
        f"mean entropy delta cyclo   = "
        f"{cyclo_entropy_mean:+.8f}"
    )

    print(
        f"mean entropy delta control = "
        f"{control_entropy_mean:+.8f}"
    )

    # ----------------------------------------------------------------------
    # 10. TARGET-LEVEL EFFECT COUNTS
    # ----------------------------------------------------------------------

    print()
    print("10. TARGET-LEVEL EFFECT COUNTS")
    print("-" * 78)

    cyclo_positive = sum(
        d > 0
        for d in cyclo_deltas
    )

    control_positive = sum(
        d > 0
        for d in control_deltas
    )

    print(
        f"cyclotomic targets delta > 0 = "
        f"{cyclo_positive}/{TARGET_COUNT}"
    )

    print(
        f"control targets delta > 0     = "
        f"{control_positive}/{TARGET_COUNT}"
    )

    print()
    print("largest cyclotomic deltas:")

    for idx in sorted(
        range(TARGET_COUNT),
        key=lambda i: abs(cyclo_deltas[i]),
        reverse=True,
    )[:10]:

        print(
            f"  target {idx+1:3d}: "
            f"C={cyclo_deltas[idx]:+.6f} "
            f"R={control_deltas[idx]:+.6f}"
        )

    # ----------------------------------------------------------------------
    # 11. ROOT SWAP TEST
    # ----------------------------------------------------------------------

    print()
    print("11. ROOT-SWAP / CONJUGATION TEST")
    print("-" * 78)

    swap_failures = swap_test(
        targets,
        cyclo_holdout[:100] + control_holdout[:100],
    )

    print(
        f"swap failures = {swap_failures}"
    )

    print(
        "status =",
        "PASS" if swap_failures == 0 else "FAIL",
    )

    # ----------------------------------------------------------------------
    # 12. FINAL DIAGNOSTIC
    # ----------------------------------------------------------------------

    print()
    print("12. FINAL DIAGNOSTIC")
    print("-" * 78)

    global_c_gap = (
        c_true_non
        - c_false_non
    )

    global_r_gap = (
        r_true_non
        - r_false_non
    )

    cyclo_control_gap = (
        cyclo_mean
        - control_mean
    )

    print(
        f"cyclotomic global delta = "
        f"{global_c_gap:+.8f}"
    )

    print(
        f"control global delta    = "
        f"{global_r_gap:+.8f}"
    )

    print(
        f"target-level C-R gap    = "
        f"{cyclo_control_gap:+.8f}"
    )

    if abs(cyclo_control_gap) < 0.01:
        verdict = (
            "NO PERSISTENT CYCLOTOMIC HOLDOUT EFFECT"
        )

    elif abs(cyclo_control_gap) < 0.03:
        verdict = (
            "WEAK / INCONCLUSIVE HOLDOUT EFFECT"
        )

    else:
        verdict = (
            "POTENTIAL CYCLOTOMIC HOLDOUT EFFECT"
        )

    print()
    print("VERDICT")
    print(verdict)

    print(
        """
Interpretation
--------------
The false null is now properly matched.

Every false sum:

    * lies in the same feasible even-sum domain;
    * satisfies exactly the original 13 training D-QR conditions;
    * is excluded from the true sum;
    * is sampled in equal quantity across every target.

The cubic-character moduli are held out from the construction of those
false pools.

Therefore:

    training cyclotomic constraints
        -> build false null

and only afterwards:

    holdout cubic character
        -> test TRUE versus FALSE.

This removes the main leakage problem from the previous experiment.

A strong positive result should show:

    cyclotomic holdout separation
        >
    random-control separation

across many targets.

A weak or inconsistent result means the earlier cubic-character
differences were probably sampling/training-family effects rather than
a new invariant of the true factor sum.

Even a positive result is a statistical signature, not yet a factorization
algorithm. The next step would then be to identify the algebraic identity
behind the holdout effect.
        """
    )

    # ----------------------------------------------------------------------
    # 13. RUNTIME
    # ----------------------------------------------------------------------

    print("=" * 78)
    print(
        "EXPERIMENT 60R COMPLETE"
    )
    print(
        f"total runtime = "
        f"{time.perf_counter() - overall_start:.6f}s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()