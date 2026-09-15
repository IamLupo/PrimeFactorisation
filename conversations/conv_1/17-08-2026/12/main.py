#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 82R
FAST DIRECT PARTITION-COEFFICIENT SIEVE FOR s = p + q

CORRECTED EXPERIMENT 82
NO FACTOR-PAIR FILTER DURING MODULAR SIEVE
GENERAL SYMMETRIC COEFFICIENT EVALUATION
BYTEARRAY / RESIDUE-CLASS SIEVE
INCREMENTAL MODULAR FINGERPRINTS
EXACT BIG-INTEGER CHECK ONLY AFTER MODULAR COLLAPSE
TARGET HOLDOUT
NO CSV OUTPUT
NO SKLEARN
==============================================================================

CORE IDEA
---------
For the paper-derived coefficient

    B_(k,l)(p,q)
      = (1+q)^l p^k - (1+q)^k p^l
        + (1+p)^l q^k - (1+p)^k q^l

and n = p*q, s = p+q, the expression is symmetric in p,q and can
be evaluated directly from n,s.

This script does NOT test whether s^2 - 4n is a square while sieving.

Instead:

    unrestricted even s-domain
              |
              v
       B(n,s) mod 7
              |
              v
       B(n,s) mod 13
              |
              v
             ...
              |
              v
       tiny surviving s-set
              |
              v
       exact integer fingerprint
              |
              v
       ONLY NOW check whether s actually factors n

This fixes the methodological error in Experiment 81.

SPEED OPTIMIZATION
------------------
For each odd modulus m, the allowed residue classes of the candidate
index t are converted into bytearray slice masks:

    survivors[t] = 1

and disallowed residue classes are cleared in bulk using

    mask[r::m] = b'\\x00' * ...

This avoids millions of Python-level modular evaluations.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 820082

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

# Even candidate-s domain.
S_MIN = 4_000_000
S_MAX = 8_400_000

# Paper coefficient family.
# Keep the family moderate initially; the fast sieve makes expansion cheap.
PAIRS: List[Tuple[int, int]] = [
    (1, 3),
    (1, 5),
    (1, 7),
    (1, 9),
    (3, 5),
    (3, 7),
]

# Incremental modular fingerprint.
MODULI = [
    7,
    13,
    19,
    31,
    37,
    61,
    67,
]

# Prefixes reported explicitly.
MODULUS_PREFIXES = [
    [7],
    [7, 13],
    [7, 13, 19],
    [7, 13, 19, 31, 37],
    [7, 13, 19, 31, 37, 61, 67],
]

PRINT_TEST_TARGETS = 10

# Stop modular sieve early when the candidate set becomes this small.
# Exact fingerprint evaluation then finishes the job.
EARLY_STOP_SURVIVORS = 256


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


@dataclass
class SieveResult:
    domain_size: int
    survivors: int
    exact_matches: int
    factor_pair_candidates: int
    true_s_survives: bool
    unique_exact: bool
    unique_factor_candidate: bool
    surviving_s: List[int]


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi < 2 or lo > hi:
        return []

    lo = max(2, lo)
    root = math.isqrt(hi)

    base = bytearray(b"\x01") * (root + 1)
    base[:2] = b"\x00\x00"

    limit = math.isqrt(root)

    for p in range(2, limit + 1):
        if base[p]:
            start = p * p
            base[start:root + 1:p] = b"\x00" * (
                ((root - start) // p) + 1
            )

    small_primes = [
        p
        for p in range(2, root + 1)
        if base[p]
    ]

    out: List[int] = []

    segment_size = 1_000_000
    start = lo

    while start <= hi:
        end = min(start + segment_size - 1, hi)
        mark = bytearray(b"\x01") * (end - start + 1)

        for p in small_primes:
            if p * p > end:
                break

            first = max(
                p * p,
                ((start + p - 1) // p) * p,
            )

            for value in range(first, end + 1, p):
                mark[value - start] = 0

        out.extend(
            start + i
            for i, flag in enumerate(mark)
            if flag
        )

        start = end + 1

    return out


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    seed: int,
) -> List[Target]:
    rng = random.Random(seed)

    result: List[Target] = []
    seen = set()

    while len(result) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if q - p < 50_000:
            continue

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        result.append(
            Target(
                idx=len(result) + 1,
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return result


# ============================================================================
# MODULAR ARITHMETIC
# ============================================================================

def mod_inverse(a: int, m: int) -> int:
    """
    Extended Euclidean inverse.
    """
    old_r, r = a, m
    old_s, s = 1, 0

    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if old_r != 1:
        raise ValueError(
            f"{a} has no inverse modulo {m}"
        )

    return old_s % m


# ============================================================================
# POWER SUMS p^d + q^d FROM n,s
# ============================================================================

def power_sums_mod(
    n: int,
    s: int,
    max_degree: int,
    modulus: int,
) -> List[int]:
    """
    Returns P[d] = p^d + q^d mod modulus, where

        p + q = s
        p q = n.

    Recurrence:

        P_0 = 2
        P_1 = s
        P_d = s P_{d-1} - n P_{d-2}.
    """
    if max_degree < 0:
        return []

    n %= modulus
    s %= modulus

    pows = [0] * (max_degree + 1)

    pows[0] = 2 % modulus

    if max_degree == 0:
        return pows

    pows[1] = s

    for d in range(2, max_degree + 1):
        pows[d] = (
            s * pows[d - 1]
            - n * pows[d - 2]
        ) % modulus

    return pows


# ============================================================================
# SYMMETRIC MONOMIAL
# ============================================================================

def symmetric_monomial(
    p_power: int,
    q_power: int,
    n: int,
    s: int,
    modulus: int | None = None,
    power_sums: Sequence[int] | None = None,
) -> int:
    """
    Computes

        p^a q^b + q^a p^b

    using only n=pq and s=p+q.

    If a >= b:

        n^b * (p^(a-b) + q^(a-b)).
    """
    a = p_power
    b = q_power

    if a == b:
        value = 2 * (n ** a)
        return value if modulus is None else value % modulus

    if a < b:
        a, b = b, a

    diff = a - b

    if modulus is None:
        if diff == 0:
            power_sum = 2
        elif diff == 1:
            power_sum = s
        else:
            prev2 = 2
            prev1 = s

            if diff == 2:
                power_sum = s * s - 2 * n
            else:
                for d in range(2, diff + 1):
                    current = s * prev1 - n * prev2
                    prev2, prev1 = prev1, current

                power_sum = prev1

        return (n ** b) * power_sum

    # Modular version.
    if power_sums is None:
        power_sums = power_sums_mod(
            n,
            s,
            diff,
            modulus,
        )

    return (
        pow(n % modulus, b, modulus)
        * power_sums[diff]
    ) % modulus


# ============================================================================
# GENERAL PAPER COEFFICIENT
# ============================================================================

def paper_coeff_mod(
    n: int,
    s: int,
    k: int,
    ell: int,
    modulus: int,
) -> int:
    """
    General symmetric evaluation of

        (1+q)^ell p^k - (1+q)^k p^ell
        +(1+p)^ell q^k - (1+p)^k q^ell.

    No factorization is needed.
    """

    max_degree = max(k, ell)

    ps = power_sums_mod(
        n,
        s,
        max_degree,
        modulus,
    )

    total = 0

    # Expand:
    #
    # (1+q)^ell p^k
    # = sum_i C(ell,i) p^k q^i.
    #
    # Adding the swapped term gives the symmetric monomial.

    for i in range(ell + 1):
        coeff = math.comb(ell, i)

        term = symmetric_monomial(
            k,
            i,
            n,
            s,
            modulus,
            ps,
        )

        total += coeff * term

    for i in range(k + 1):
        coeff = math.comb(k, i)

        term = symmetric_monomial(
            ell,
            i,
            n,
            s,
            modulus,
            ps,
        )

        total -= coeff * term

    return total % modulus


def paper_coeff_exact(
    n: int,
    s: int,
    k: int,
    ell: int,
) -> int:
    """
    Exact integer version of the same symmetric expression.
    """
    max_degree = max(k, ell)

    ps = [0] * (max_degree + 1)
    ps[0] = 2

    if max_degree >= 1:
        ps[1] = s

    for d in range(2, max_degree + 1):
        ps[d] = (
            s * ps[d - 1]
            - n * ps[d - 2]
        )

    total = 0

    for i in range(ell + 1):
        total += (
            math.comb(ell, i)
            * symmetric_monomial(
                k,
                i,
                n,
                s,
                None,
                ps,
            )
        )

    for i in range(k + 1):
        total -= (
            math.comb(k, i)
            * symmetric_monomial(
                ell,
                i,
                n,
                s,
                None,
                ps,
            )
        )

    return total


# ============================================================================
# DIRECT p,q COEFFICIENT — VALIDATION ONLY
# ============================================================================

def paper_coeff_pq(
    p: int,
    q: int,
    k: int,
    ell: int,
) -> int:
    return (
        (1 + q) ** ell * p ** k
        - (1 + q) ** k * p ** ell
        + (1 + p) ** ell * q ** k
        - (1 + p) ** k * q ** ell
    )


# ============================================================================
# FINGERPRINTS
# ============================================================================

def exact_fingerprint(
    n: int,
    s: int,
    pairs: Sequence[Tuple[int, int]],
) -> Tuple[int, ...]:
    return tuple(
        paper_coeff_exact(n, s, k, ell)
        for k, ell in pairs
    )


def modular_fingerprint(
    n: int,
    s: int,
    pairs: Sequence[Tuple[int, int]],
    modulus: int,
) -> Tuple[int, ...]:
    return tuple(
        paper_coeff_mod(
            n,
            s,
            k,
            ell,
            modulus,
        )
        for k, ell in pairs
    )


# ============================================================================
# EVEN S DOMAIN
# ============================================================================

def domain_parameters() -> Tuple[int, int, int]:
    lo = S_MIN

    if lo & 1:
        lo += 1

    hi = S_MAX

    if hi & 1:
        hi -= 1

    count = ((hi - lo) // 2) + 1

    return lo, hi, count


# ============================================================================
# FACTOR PAIR CHECK
# ONLY AFTER FINGERPRINT MATCHING
# ============================================================================

def factor_pair_from_s(
    n: int,
    s: int,
) -> Tuple[int, int] | None:
    D = s * s - 4 * n

    if D < 0:
        return None

    r = math.isqrt(D)

    if r * r != D:
        return None

    if (s - r) & 1:
        return None

    p = (s - r) // 2
    q = (s + r) // 2

    if p <= 1 or q <= 1:
        return None

    if p * q != n:
        return None

    return p, q


# ============================================================================
# RESIDUE-CLASS PREPARATION
# ============================================================================

def allowed_s_residues(
    n: int,
    true_s: int,
    pairs: Sequence[Tuple[int, int]],
    modulus: int,
) -> bytearray:
    """
    allowed[r] = 1 iff s == r (mod modulus) can have the same modular
    fingerprint as the true s.
    """

    oracle = modular_fingerprint(
        n,
        true_s,
        pairs,
        modulus,
    )

    allowed = bytearray(modulus)

    for r in range(modulus):
        fp = modular_fingerprint(
            n,
            r,
            pairs,
            modulus,
        )

        if fp == oracle:
            allowed[r] = 1

    return allowed


# ============================================================================
# FAST BYTEARRAY SIEVE
# ============================================================================

def fast_modular_sieve(
    n: int,
    true_s: int,
    pairs: Sequence[Tuple[int, int]],
    moduli: Sequence[int],
    early_stop: int = EARLY_STOP_SURVIVORS,
    report: bool = False,
) -> Tuple[List[int], List[Tuple[int, int]]]:
    """
    Fast sieve over the complete even-s domain.

    Candidate indexing:

        s = S_MIN + 2*t

    For an odd modulus m, t has one residue class for each possible
    s residue because 2 is invertible mod m.

    Disallowed residue classes are removed by bytearray slicing.
    """

    lo, hi, count = domain_parameters()

    # Start with all candidates alive.
    mask = bytearray(b"\x01") * count

    # Track survivors after every modulus.
    history: List[Tuple[int, int]] = []

    for modulus in moduli:

        allowed_s = allowed_s_residues(
            n,
            true_s,
            pairs,
            modulus,
        )

        # Map allowed s residues to allowed t residues.
        inv2 = mod_inverse(2, modulus)

        allowed_t = bytearray(modulus)

        for s_residue in range(modulus):
            if allowed_s[s_residue]:
                t_residue = (
                    (s_residue - (lo % modulus))
                    * inv2
                ) % modulus

                allowed_t[t_residue] = 1

        # Remove every disallowed residue class in t.
        zero_block_cache: Dict[int, bytes] = {}

        for r in range(modulus):
            if allowed_t[r]:
                continue

            length = (
                (count - 1 - r) // modulus + 1
                if r < count
                else 0
            )

            if length <= 0:
                continue

            zeros = zero_block_cache.get(length)

            if zeros is None:
                zeros = b"\x00" * length
                zero_block_cache[length] = zeros

            mask[r:count:modulus] = zeros

        survivors = mask.count(1)

        history.append(
            (modulus, survivors)
        )

        if report:
            allowed_count = allowed_t.count(1)

            print(
                f"    m={modulus:4d} "
                f"allowed_residues={allowed_count:3d} "
                f"survivors={survivors:9d}"
            )

        if survivors <= early_stop:
            break

    # Efficiently extract surviving s values.
    surviving_s: List[int] = []

    pos = mask.find(1)

    while pos != -1:
        surviving_s.append(lo + 2 * pos)
        pos = mask.find(1, pos + 1)

    return surviving_s, history


# ============================================================================
# EXACT FINAL CHECK
# ============================================================================

def exact_filter_survivors(
    n: int,
    oracle_fp: Tuple[int, ...],
    survivors: Sequence[int],
    pairs: Sequence[Tuple[int, int]],
) -> List[int]:
    matches: List[int] = []

    for s in survivors:
        fp = exact_fingerprint(
            n,
            s,
            pairs,
        )

        if fp == oracle_fp:
            matches.append(s)

    return matches


# ============================================================================
# FULL TARGET RUN
# ============================================================================

def run_target(
    target: Target,
    pairs: Sequence[Tuple[int, int]],
    moduli: Sequence[int],
    report: bool = False,
) -> SieveResult:

    _, _, domain_size = domain_parameters()

    oracle_exact = exact_fingerprint(
        target.n,
        target.s,
        pairs,
    )

    survivors, history = fast_modular_sieve(
        target.n,
        target.s,
        pairs,
        moduli,
        report=report,
    )

    exact_matches_list = exact_filter_survivors(
        target.n,
        oracle_exact,
        survivors,
        pairs,
    )

    factor_candidates = 0

    for s in exact_matches_list:
        if factor_pair_from_s(
            target.n,
            s,
        ) is not None:
            factor_candidates += 1

    return SieveResult(
        domain_size=domain_size,
        survivors=len(survivors),
        exact_matches=len(exact_matches_list),
        factor_pair_candidates=factor_candidates,
        true_s_survives=(
            target.s in exact_matches_list
        ),
        unique_exact=(
            len(exact_matches_list) == 1
        ),
        unique_factor_candidate=(
            factor_candidates == 1
            and target.s in exact_matches_list
        ),
        surviving_s=exact_matches_list,
    )


# ============================================================================
# VALIDATE FORMULA FAMILY
# ============================================================================

def validate_coefficients(
    targets: Sequence[Target],
    pairs: Sequence[Tuple[int, int]],
) -> int:

    failures = 0

    print("\n3. SYMMETRIC POLYNOMIAL VALIDATION")
    print("-" * 78)

    for t in targets:

        ok = True

        for k, ell in pairs:
            direct = paper_coeff_pq(
                t.p,
                t.q,
                k,
                ell,
            )

            symmetric = paper_coeff_exact(
                t.n,
                t.s,
                k,
                ell,
            )

            if direct != symmetric:
                ok = False
                failures += 1

        print(
            f"target {t.idx:3d}: "
            f"coefficients_ok={ok}"
        )

    print(
        f"identity failures = {failures}"
    )

    if failures:
        raise RuntimeError(
            "General symmetric coefficient evaluation failed."
        )

    print("status = PASS")

    return failures


# ============================================================================
# PREFIX EXPERIMENT
# ============================================================================

def run_prefix_experiment(
    target: Target,
    prefixes: Sequence[Sequence[int]],
    pairs: Sequence[Tuple[int, int]],
) -> List[Tuple[List[int], int, int]]:
    """
    Returns prefix, final survivor count, exact-match count.
    """

    output = []

    for prefix in prefixes:

        survivors, _ = fast_modular_sieve(
            target.n,
            target.s,
            pairs,
            prefix,
            early_stop=EARLY_STOP_SURVIVORS,
            report=False,
        )

        oracle = exact_fingerprint(
            target.n,
            target.s,
            pairs,
        )

        exact_matches = exact_filter_survivors(
            target.n,
            oracle,
            survivors,
            pairs,
        )

        output.append(
            (
                list(prefix),
                len(survivors),
                len(exact_matches),
            )
        )

    return output


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    overall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 82R")
    print("FAST DIRECT PARTITION-COEFFICIENT SIEVE FOR s = p + q")
    print("BYTEARRAY RESIDUE-CLASS SIEVE")
    print("EXACT CHECK ONLY AFTER MODULAR COLLAPSE")
    print("TARGET HOLDOUT")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------
    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------
    targets = generate_targets(
        primes,
        NUM_TARGETS,
        SEED,
    )

    print("\n2. TARGET SUMMARY")
    print("-" * 78)

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print(
            "... remaining generated targets omitted"
        )

    # ------------------------------------------------------------------
    # Validate formulas
    # ------------------------------------------------------------------
    validate_coefficients(
        targets,
        PAIRS,
    )

    # ------------------------------------------------------------------
    # Holdout
    # ------------------------------------------------------------------
    train = targets[:TRAIN_TARGETS]
    test = targets[TRAIN_TARGETS:]

    print("\n4. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train)}"
    )
    print(
        f"test targets     = {len(test)}"
    )

    # ------------------------------------------------------------------
    # Domain
    # ------------------------------------------------------------------
    lo, hi, domain_size = domain_parameters()

    print("\n5. S-CANDIDATE DOMAIN")
    print("-" * 78)
    print(
        f"s minimum         = {lo}"
    )
    print(
        f"s maximum         = {hi}"
    )
    print(
        f"even candidates   = {domain_size}"
    )
    print(
        "factor-pair test during modular sieve = DISABLED"
    )

    # ------------------------------------------------------------------
    # Main full modular sieve
    # ------------------------------------------------------------------
    print("\n6. FAST MODULAR PARTITION SIEVE")
    print("-" * 78)

    all_results: List[Tuple[Target, SieveResult]] = []

    for i, target in enumerate(
        targets,
        start=1,
    ):

        t0 = time.perf_counter()

        result = run_target(
            target,
            PAIRS,
            MODULI,
            report=(i == TRAIN_TARGETS + 1),
        )

        elapsed = (
            time.perf_counter()
            - t0
        )

        all_results.append(
            (target, result)
        )

        if (
            i <= 3
            or i % 5 == 0
            or i == NUM_TARGETS
        ):
            print(
                f"processed {i:3d}/{NUM_TARGETS} "
                f"target={target.idx:3d} "
                f"modular_survivors="
                f"{result.survivors:5d} "
                f"exact_matches="
                f"{result.exact_matches:3d} "
                f"time={elapsed:.4f}s"
            )

    # ------------------------------------------------------------------
    # Train/test metrics
    # ------------------------------------------------------------------
    train_rows = [
        row
        for row in all_results
        if row[0].idx <= TRAIN_TARGETS
    ]

    test_rows = [
        row
        for row in all_results
        if row[0].idx > TRAIN_TARGETS
    ]

    print("\n7. TRAINING RESULTS")
    print("-" * 78)

    print(
        f"mean modular survivors = "
        f"{statistics.fmean(r.survivors for _, r in train_rows):.3f}"
    )

    print(
        f"median modular survivors = "
        f"{statistics.median(r.survivors for _, r in train_rows):.3f}"
    )

    print(
        f"mean exact matches = "
        f"{statistics.fmean(r.exact_matches for _, r in train_rows):.3f}"
    )

    print(
        f"true-s survival = "
        f"{sum(r.true_s_survives for _, r in train_rows)}/{len(train_rows)}"
    )

    print(
        f"unique exact = "
        f"{sum(r.unique_exact for _, r in train_rows)}/{len(train_rows)}"
    )

    print("\n8. OUT-OF-SAMPLE RESULTS")
    print("-" * 78)

    print(
        f"mean modular survivors = "
        f"{statistics.fmean(r.survivors for _, r in test_rows):.3f}"
    )

    print(
        f"median modular survivors = "
        f"{statistics.median(r.survivors for _, r in test_rows):.3f}"
    )

    print(
        f"mean exact matches = "
        f"{statistics.fmean(r.exact_matches for _, r in test_rows):.3f}"
    )

    print(
        f"true-s survival = "
        f"{sum(r.true_s_survives for _, r in test_rows)}/{len(test_rows)}"
    )

    print(
        f"unique exact = "
        f"{sum(r.unique_exact for _, r in test_rows)}/{len(test_rows)}"
    )

    print(
        f"unique factor candidate = "
        f"{sum(r.unique_factor_candidate for _, r in test_rows)}/{len(test_rows)}"
    )

    # ------------------------------------------------------------------
    # Test target details
    # ------------------------------------------------------------------
    print("\n9. TEST TARGET DETAILS")
    print("-" * 78)

    for target, result in test_rows[:PRINT_TEST_TARGETS]:

        print(
            f"target {target.idx:3d}: "
            f"modular={result.survivors:6d} "
            f"exact={result.exact_matches:4d} "
            f"factor={result.factor_pair_candidates:2d} "
            f"true={result.true_s_survives} "
            f"unique_exact={result.unique_exact}"
        )

        if result.exact_matches <= 20:
            print(
                f"    exact s candidates = "
                f"{result.surviving_s}"
            )

    # ------------------------------------------------------------------
    # Prefix study
    # ------------------------------------------------------------------
    print("\n10. MODULUS PREFIX STUDY")
    print("-" * 78)

    prefix_test_rows = test[:]

    for prefix in MODULUS_PREFIXES:

        survivor_counts = []
        exact_counts = []
        survival_count = 0

        for target in prefix_test_rows:

            survivors, _ = fast_modular_sieve(
                target.n,
                target.s,
                PAIRS,
                prefix,
                early_stop=EARLY_STOP_SURVIVORS,
                report=False,
            )

            oracle = exact_fingerprint(
                target.n,
                target.s,
                PAIRS,
            )

            exact = exact_filter_survivors(
                target.n,
                oracle,
                survivors,
                PAIRS,
            )

            survivor_counts.append(
                len(survivors)
            )

            exact_counts.append(
                len(exact)
            )

            if target.s in exact:
                survival_count += 1

        print(
            f"prefix={prefix}"
        )

        print(
            f"    mean modular survivors = "
            f"{statistics.fmean(survivor_counts):.3f}"
        )

        print(
            f"    median modular survivors = "
            f"{statistics.median(survivor_counts):.3f}"
        )

        print(
            f"    mean exact matches = "
            f"{statistics.fmean(exact_counts):.3f}"
        )

        print(
            f"    true-s survival = "
            f"{survival_count}/{len(prefix_test_rows)}"
        )

    # ------------------------------------------------------------------
    # IMPORTANT baseline
    # ------------------------------------------------------------------
    print("\n11. CONTROL: NO FINGERPRINT")
    print("-" * 78)

    print(
        f"unrestricted even-s candidates = "
        f"{domain_size}"
    )

    print(
        "This is the correct baseline."
    )

    print(
        "The factor-pair condition is NOT used to produce this baseline."
    )

    # ------------------------------------------------------------------
    # Final interpretation
    # ------------------------------------------------------------------
    test_mean_modular = statistics.fmean(
        r.survivors
        for _, r in test_rows
    )

    test_mean_exact = statistics.fmean(
        r.exact_matches
        for _, r in test_rows
    )

    test_true_survival = sum(
        r.true_s_survives
        for _, r in test_rows
    )

    test_unique = sum(
        r.unique_exact
        for _, r in test_rows
    )

    test_unique_factor = sum(
        r.unique_factor_candidate
        for _, r in test_rows
    )

    print("\n12. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        f"unrestricted even-s domain = "
        f"{domain_size}"
    )

    print(
        f"test mean modular survivors = "
        f"{test_mean_modular:.6f}"
    )

    print(
        f"test mean exact matches = "
        f"{test_mean_exact:.6f}"
    )

    print(
        f"test true-s survival = "
        f"{test_true_survival}/{len(test_rows)}"
    )

    print(
        f"test unique exact recovery = "
        f"{test_unique}/{len(test_rows)}"
    )

    print(
        f"test unique factor recovery = "
        f"{test_unique_factor}/{len(test_rows)}"
    )

    reduction_ratio = (
        domain_size / test_mean_modular
        if test_mean_modular > 0
        else float("inf")
    )

    print(
        f"mean modular search-space reduction = "
        f"{reduction_ratio:.3f}x"
    )

    print()
    print(
        "INTERPRETATION"
    )
    print(
        "--------------"
    )
    print(
        "A genuine positive result requires candidate-space collapse "
        "BEFORE any factor-pair test."
    )
    print()
    print(
        "The strongest result would be:"
    )
    print(
        "    millions of even s candidates"
    )
    print(
        "          -> small modular survivor set"
    )
    print(
        "          -> exact fingerprint leaves very few s"
    )
    print(
        "          -> true s survives on unseen targets"
    )
    print()
    print(
        "If the modular sieve reduces the domain but the exact "
        "fingerprint still leaves many candidates, we have useful "
        "partial information but not reconstruction."
    )

    if (
        test_true_survival == len(test_rows)
        and test_unique == len(test_rows)
        and test_mean_modular < 1000
    ):
        status = (
            "STRONG RESULT: the partition fingerprint sharply "
            "identifies s on unseen targets."
        )
    elif (
        test_true_survival == len(test_rows)
        and test_mean_modular < domain_size / 100
    ):
        status = (
            "INTERESTING RESULT: substantial N-conditioned "
            "candidate-space collapse survives target holdout."
        )
    elif test_true_survival == len(test_rows):
        status = (
            "INFORMATIONAL RESULT: the fingerprint retains the "
            "true s but leaves a large candidate set."
        )
    else:
        status = (
            "NO ROBUST RESULT: the fingerprint does not reliably "
            "retain the true s."
        )

    print(
        f"\nSTATUS = {status}"
    )

    runtime = (
        time.perf_counter()
        - overall_start
    )

    print(
        f"\ntotal runtime = {runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 82R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()