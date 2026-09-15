#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 85R
PATCHED EISENSTEIN / QUASIMODULAR TOWER RESOLVENT
H6 + H8 + H10 + H12

FIX:
    Exact candidate comparison now uses the fixed-denominator NUMERATORS,
    so arbitrary false candidate s values cannot trigger an integrality
    exception.

PAPER-DERIVED H_k COEFFICIENTS
SEMIPRIME REDUCTION TO (n, s)
MULTI-COEFFICIENT MODULAR SIEVE
STRICT TARGET HOLDOUT
CYCLOTOMIC VS CONTROL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

IMPORTANT:
    This is an oracle information experiment.

    The true H_k values are constructed from the known p,q target pair.
    Candidate reconstruction uses only n and candidate s.

    Positive results establish identifiability / information content,
    not yet an N-only factorization algorithm.
==============================================================================
"""

from __future__ import annotations

import math
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

NUM_TARGETS = 40
TRAIN_TARGETS = 30

P_MIN = 2_000_000
P_MAX = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_400_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37,
    61, 67,
]

CONTROL = [
    673, 4561, 4759,
    6211, 7879,
    7951, 8689,
]

RNG_SEED = 85085

PRINT_LIMIT = 20


# ============================================================================
# TARGET
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


TowerNumerators = Tuple[int, int, int, int]
TowerValues = Tuple[int, int, int, int]


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

    if hi < 2 or lo > hi:
        return []

    lo = max(lo, 2)

    sieve = bytearray(
        b"\x01" * (hi + 1)
    )

    sieve[0:2] = b"\x00\x00"

    for p in range(
        2,
        math.isqrt(hi) + 1,
    ):
        if sieve[p]:
            start = p * p
            count = (
                (hi - start) // p
            ) + 1

            sieve[start:hi + 1:p] = (
                b"\x00" * count
            )

    return [
        x
        for x in range(lo, hi + 1)
        if sieve[x]
    ]


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
) -> List[Target]:

    state = RNG_SEED
    mask = (1 << 64) - 1

    out: List[Target] = []
    seen = set()

    N = len(primes)

    while len(out) < count:

        state = (
            6364136223846793005 * state
            + 1442695040888963407
        ) & mask

        i = state % N

        state = (
            6364136223846793005 * state
            + 1442695040888963407
        ) & mask

        j = state % N

        if i == j:
            continue

        p = primes[i]
        q = primes[j]

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)

        out.append(
            Target(
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return out


# ============================================================================
# POWER SUMS
# ============================================================================

def power_sums_from_ns(
    n: int,
    s: int,
    max_power: int,
) -> List[int]:
    """
    R_j = p^j + q^j.

    R_0 = 2
    R_1 = s
    R_j = s R_(j-1) - n R_(j-2)
    """

    R = [0] * (max_power + 1)

    R[0] = 2

    if max_power >= 1:
        R[1] = s

    for j in range(
        2,
        max_power + 1,
    ):
        R[j] = (
            s * R[j - 1]
            - n * R[j - 2]
        )

    return R


def sigma_prime_product(
    n: int,
    s: int,
    j: int,
) -> int:
    """
    For n=pq:

        sigma_j(n)
          = (1+p^j)(1+q^j)
          = 1 + p^j + q^j + n^j.
    """

    R = power_sums_from_ns(
        n,
        s,
        j,
    )

    return (
        1
        + R[j]
        + n ** j
    )


# ============================================================================
# PAPER-DERIVED H_k NUMERATORS
# ============================================================================

def H6_numerator_from_ns(
    n: int,
    s: int,
) -> int:
    """
    H6 = N6 / 6
    """

    sig1 = sigma_prime_product(
        n,
        s,
        1,
    )

    sig3 = sigma_prime_product(
        n,
        s,
        3,
    )

    return (
        (n * n - n + 1) * sig1
        - sig3
    )


def H_even_numerator_from_ns(
    n: int,
    s: int,
    k: int,
) -> int:
    """
    For k >= 8 even:

        H_k = N_k / 24

        N_k =
            -n^2 sigma_(k-7)
            +(n^2+1) sigma_(k-5)
            -sigma_(k-3)
    """

    if k < 8 or k % 2:
        raise ValueError(
            "k must be even and >= 8"
        )

    a = k - 7
    b = k - 5
    c = k - 3

    sa = sigma_prime_product(
        n,
        s,
        a,
    )

    sb = sigma_prime_product(
        n,
        s,
        b,
    )

    sc = sigma_prime_product(
        n,
        s,
        c,
    )

    return (
        -n * n * sa
        + (n * n + 1) * sb
        - sc
    )


def tower_numerators_from_ns(
    n: int,
    s: int,
) -> TowerNumerators:

    return (
        H6_numerator_from_ns(n, s),
        H_even_numerator_from_ns(n, s, 8),
        H_even_numerator_from_ns(n, s, 10),
        H_even_numerator_from_ns(n, s, 12),
    )


def tower_values_from_ns(
    n: int,
    s: int,
) -> TowerValues:
    """
    Integer-valued tower.

    Only valid when the corresponding divisibility conditions hold.
    Intended for the true target, not arbitrary candidate s.
    """

    nums = tower_numerators_from_ns(
        n,
        s,
    )

    if nums[0] % 6 != 0:
        raise ArithmeticError(
            "H6 numerator is not divisible by 6"
        )

    if nums[1] % 24 != 0:
        raise ArithmeticError(
            "H8 numerator is not divisible by 24"
        )

    if nums[2] % 24 != 0:
        raise ArithmeticError(
            "H10 numerator is not divisible by 24"
        )

    if nums[3] % 24 != 0:
        raise ArithmeticError(
            "H12 numerator is not divisible by 24"
        )

    return (
        nums[0] // 6,
        nums[1] // 24,
        nums[2] // 24,
        nums[3] // 24,
    )


# ============================================================================
# INDEPENDENT TARGET VALIDATION
# ============================================================================

def direct_sigma_for_target(
    t: Target,
    j: int,
) -> int:
    return (
        (1 + t.p ** j)
        * (1 + t.q ** j)
    )


def validate_target(
    t: Target,
) -> Tuple[bool, TowerValues]:
    """
    Validate:
        n=pq
        s=p+q
        H6/H8/H10/H12
    using direct p,q formulas.
    """

    n = t.n
    s = t.s

    nums = tower_numerators_from_ns(
        n,
        s,
    )

    if nums[0] % 6 != 0:
        return (
            False,
            (0, 0, 0, 0),
        )

    if (
        nums[1] % 24
        or nums[2] % 24
        or nums[3] % 24
    ):
        return (
            False,
            (0, 0, 0, 0),
        )

    values = (
        nums[0] // 6,
        nums[1] // 24,
        nums[2] // 24,
        nums[3] // 24,
    )

    # Direct p,q verification.
    sig1 = direct_sigma_for_target(t, 1)
    sig3 = direct_sigma_for_target(t, 3)
    sig5 = direct_sigma_for_target(t, 5)
    sig7 = direct_sigma_for_target(t, 7)
    sig9 = direct_sigma_for_target(t, 9)

    direct_h6_num = (
        (n * n - n + 1) * sig1
        - sig3
    )

    direct_h8_num = (
        -n * n * sig1
        + (n * n + 1) * sig3
        - sig5
    )

    direct_h10_num = (
        -n * n * sig3
        + (n * n + 1) * sig5
        - sig7
    )

    direct_h12_num = (
        -n * n * sig5
        + (n * n + 1) * sig7
        - sig9
    )

    direct = (
        direct_h6_num // 6,
        direct_h8_num // 24,
        direct_h10_num // 24,
        direct_h12_num // 24,
    )

    return (
        direct == values,
        values,
    )


# ============================================================================
# MODULAR TOWER
# ============================================================================

def tower_mod(
    n_mod: int,
    s_mod: int,
    ell: int,
) -> TowerValues:
    """
    Evaluate H6/H8/H10/H12 modulo ell.

    This is safe for the tested primes because none divide 6 or 24.
    """

    n = n_mod % ell
    s = s_mod % ell

    R = [0] * 10

    R[0] = 2 % ell
    R[1] = s

    for j in range(
        2,
        10,
    ):
        R[j] = (
            s * R[j - 1]
            - n * R[j - 2]
        ) % ell

    def sigma(j: int) -> int:
        return (
            1
            + R[j]
            + pow(n, j, ell)
        ) % ell

    inv6 = pow(
        6,
        ell - 2,
        ell,
    )

    inv24 = pow(
        24,
        ell - 2,
        ell,
    )

    h6 = (
        (
            (n * n - n + 1)
            * sigma(1)
            - sigma(3)
        )
        * inv6
    ) % ell

    h8 = (
        (
            -n * n * sigma(1)
            + (n * n + 1) * sigma(3)
            - sigma(5)
        )
        * inv24
    ) % ell

    h10 = (
        (
            -n * n * sigma(3)
            + (n * n + 1) * sigma(5)
            - sigma(7)
        )
        * inv24
    ) % ell

    h12 = (
        (
            -n * n * sigma(5)
            + (n * n + 1) * sigma(7)
            - sigma(9)
        )
        * inv24
    ) % ell

    return (
        h6,
        h8,
        h10,
        h12,
    )


# ============================================================================
# RESIDUE TABLE
# ============================================================================

def build_allowed_table(
    t: Target,
    ell: int,
    tower: bool,
) -> bytearray:

    true_tower = tower_values_from_ns(
        t.n,
        t.s,
    )

    wanted = tuple(
        value % ell
        for value in true_tower
    )

    table = bytearray(
        b"\x00" * ell
    )

    n_mod = t.n % ell

    for r in range(ell):

        signature = tower_mod(
            n_mod,
            r,
            ell,
        )

        if tower:
            ok = (
                signature
                == wanted
            )
        else:
            ok = (
                signature[0]
                == wanted[0]
            )

        if ok:
            table[r] = 1

    return table


# ============================================================================
# CANDIDATE DOMAIN
# ============================================================================

def even_candidates(
    lo: int,
    hi: int,
) -> Iterable[int]:

    start = (
        lo
        if lo % 2 == 0
        else lo + 1
    )

    for s in range(
        start,
        hi + 1,
        2,
    ):
        yield s


# ============================================================================
# MODULAR SIEVE
# ============================================================================

def modular_survivors(
    t: Target,
    moduli: Sequence[int],
    tower: bool,
) -> List[int]:
    """
    Return candidate s values surviving all modular signatures.
    """

    tables = [
        (
            ell,
            build_allowed_table(
                t,
                ell,
                tower,
            )
        )
        for ell in moduli
    ]

    survivors: List[int] = []

    for s in even_candidates(
        S_MIN,
        S_MAX,
    ):

        ok = True

        for ell, table in tables:

            if table[s % ell] == 0:
                ok = False
                break

        if ok:
            survivors.append(s)

    return survivors


# ============================================================================
# EXACT NUMERATOR COMPARISON
# ============================================================================

def exact_tower_match(
    t: Target,
    candidate_s: int,
    true_numerators: TowerNumerators,
) -> bool:
    """
    PATCHED:

    DO NOT compute integer H_k values for arbitrary candidate s.

    False candidates are allowed to have numerators that are not divisible
    by 6 or 24.

    Since each H_k has a fixed denominator, equality of H_k values is
    equivalent to equality of the numerators.

        H6  = N6  / 6
        H8  = N8  / 24
        H10 = N10 / 24
        H12 = N12 / 24

    Therefore compare the integer numerator tuple directly.
    """

    candidate_numerators = (
        tower_numerators_from_ns(
            t.n,
            candidate_s,
        )
    )

    return (
        candidate_numerators
        == true_numerators
    )


def exact_h6_match(
    t: Target,
    candidate_s: int,
    true_h6_numerator: int,
) -> bool:

    return (
        H6_numerator_from_ns(
            t.n,
            candidate_s,
        )
        == true_h6_numerator
    )


# ============================================================================
# RESULT SUMMARY
# ============================================================================

@dataclass
class FamilyStats:
    name: str
    train_counts: List[int]
    test_counts: List[int]
    train_unique: int
    test_unique: int
    train_survival: int
    test_survival: int


# ============================================================================
# FAMILY RUN
# ============================================================================

def run_family(
    name: str,
    moduli: Sequence[int],
    tower: bool,
    train: Sequence[Target],
    test: Sequence[Target],
) -> FamilyStats:

    print("\n")
    print("=" * 78)
    print(name)
    print(
        f"moduli = {list(moduli)}"
    )
    print(
        f"features = "
        f"{'H6+H8+H10+H12' if tower else 'H6 only'}"
    )
    print("=" * 78)

    train_counts: List[int] = []
    test_counts: List[int] = []

    train_unique = 0
    test_unique = 0

    train_survival = 0
    test_survival = 0

    t0 = time.perf_counter()

    for target in train:

        survivors = modular_survivors(
            target,
            moduli,
            tower,
        )

        train_counts.append(
            len(survivors)
        )

        if target.s in survivors:
            train_survival += 1

        if len(survivors) == 1:
            train_unique += 1

    for target in test:

        survivors = modular_survivors(
            target,
            moduli,
            tower,
        )

        test_counts.append(
            len(survivors)
        )

        if target.s in survivors:
            test_survival += 1

        if len(survivors) == 1:
            test_unique += 1

    elapsed = (
        time.perf_counter()
        - t0
    )

    print(
        f"runtime = {elapsed:.6f}s"
    )

    print("\nTRAIN")
    print(
        f"mean survivors = "
        f"{statistics.fmean(train_counts):.3f}"
    )
    print(
        f"median survivors = "
        f"{statistics.median(train_counts):.3f}"
    )
    print(
        f"unique recovery = "
        f"{train_unique}/{len(train)}"
    )
    print(
        f"true survival = "
        f"{train_survival}/{len(train)}"
    )

    print("\nTEST")
    print(
        f"mean survivors = "
        f"{statistics.fmean(test_counts):.3f}"
    )
    print(
        f"median survivors = "
        f"{statistics.median(test_counts):.3f}"
    )
    print(
        f"unique recovery = "
        f"{test_unique}/{len(test)}"
    )
    print(
        f"true survival = "
        f"{test_survival}/{len(test)}"
    )

    return FamilyStats(
        name=name,
        train_counts=train_counts,
        test_counts=test_counts,
        train_unique=train_unique,
        test_unique=test_unique,
        train_survival=train_survival,
        test_survival=test_survival,
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    overall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 85R")
    print("PATCHED EISENSTEIN / QUASIMODULAR TOWER RESOLVENT")
    print("H6 + H8 + H10 + H12")
    print("NUMERATOR-SAFE EXACT RECONSTRUCTION")
    print("STRICT TARGET HOLDOUT")
    print("CYCLOTOMIC VS CONTROL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    prime_time = (
        time.perf_counter()
        - t0
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = "
        f"{len(primes)}"
    )
    print(
        f"generation time = "
        f"{prime_time:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    for i, target in enumerate(
        targets[:24],
        1,
    ):
        print(
            f"target {i:3d}: "
            f"p={target.p} "
            f"q={target.q} "
            f"n={target.n} "
            f"s={target.s}"
        )

    if NUM_TARGETS > 24:
        print(
            "... remaining generated targets omitted"
        )

    # ------------------------------------------------------------------
    # Identity validation
    # ------------------------------------------------------------------

    print("\n2. TOWER IDENTITY VALIDATION")
    print("-" * 78)

    failures = 0
    true_numerators: Dict[int, TowerNumerators] = {}
    true_values: Dict[int, TowerValues] = {}

    for i, target in enumerate(
        targets,
        1,
    ):

        ok, values = validate_target(
            target,
        )

        if not ok:
            failures += 1

        nums = tower_numerators_from_ns(
            target.n,
            target.s,
        )

        true_numerators[i] = nums
        true_values[i] = values

        if i <= 24:
            print(
                f"target {i:3d}: "
                f"tower_identity={ok} "
                f"H6={values[0]} "
                f"H8={values[1]}"
            )

    print(
        f"identity failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Tower identity validation failed"
        )

    print("status = PASS")

    # ------------------------------------------------------------------
    # Structural statement
    # ------------------------------------------------------------------

    print("\n3. STRUCTURAL REDUCTION")
    print("-" * 78)
    print(
        "For n=pq and s=p+q:"
    )
    print(
        "  p^j + q^j = R_j"
    )
    print(
        "  R_0 = 2"
    )
    print(
        "  R_1 = s"
    )
    print(
        "  R_j = s*R_(j-1) - n*R_(j-2)"
    )
    print(
        "Therefore every H_k in this experiment reduces "
        "to a polynomial in n and s."
    )

    print()
    print(
        "H6 = E1 / 6"
    )
    print(
        "E1 = s*((n+1)^2 - s^2)"
    )

    # ------------------------------------------------------------------
    # Holdout
    # ------------------------------------------------------------------

    train = targets[
        :TRAIN_TARGETS
    ]

    test = targets[
        TRAIN_TARGETS:
    ]

    print("\n4. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = "
        f"{len(train)}"
    )
    print(
        f"test targets = "
        f"{len(test)}"
    )

    # ------------------------------------------------------------------
    # Families
    # ------------------------------------------------------------------

    family_specs = [
        (
            "H6_ONLY_C3",
            CYCLOTOMIC[:3],
            False,
        ),
        (
            "TOWER_C3",
            CYCLOTOMIC[:3],
            True,
        ),
        (
            "H6_ONLY_C5",
            CYCLOTOMIC[:5],
            False,
        ),
        (
            "TOWER_C5",
            CYCLOTOMIC[:5],
            True,
        ),
        (
            "H6_ONLY_C7",
            CYCLOTOMIC,
            False,
        ),
        (
            "TOWER_C7",
            CYCLOTOMIC,
            True,
        ),
        (
            "TOWER_R3",
            CONTROL[:3],
            True,
        ),
        (
            "TOWER_R5",
            CONTROL[:5],
            True,
        ),
        (
            "TOWER_R7",
            CONTROL,
            True,
        ),
    ]

    family_results: Dict[
        str,
        FamilyStats,
    ] = {}

    for name, moduli, tower in family_specs:

        stats = run_family(
            name,
            moduli,
            tower,
            train,
            test,
        )

        family_results[name] = stats

    # ------------------------------------------------------------------
    # Exact tower reconstruction
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("5. EXACT TOWER RECONSTRUCTION")
    print("=" * 78)

    for family_name, moduli in [
        (
            "C3_TOWER",
            CYCLOTOMIC[:3],
        ),
        (
            "C5_TOWER",
            CYCLOTOMIC[:5],
        ),
        (
            "C7_TOWER",
            CYCLOTOMIC,
        ),
    ]:

        print("\n")
        print(family_name)
        print(
            f"moduli = {list(moduli)}"
        )
        print("-" * 78)

        exact_counts = []

        unique = 0
        true_mod_survival = 0
        true_exact_survival = 0

        for local_index, target in enumerate(
            test,
            1,
        ):

            survivors = modular_survivors(
                target,
                moduli,
                True,
            )

            if target.s in survivors:
                true_mod_survival += 1

            nums = true_numerators[
                TRAIN_TARGETS + local_index
            ]

            exact = [
                s
                for s in survivors
                if exact_tower_match(
                    target,
                    s,
                    nums,
                )
            ]

            exact_counts.append(
                len(exact)
            )

            if target.s in exact:
                true_exact_survival += 1

            if len(exact) == 1:
                unique += 1

        print(
            f"mean exact matches = "
            f"{statistics.fmean(exact_counts):.3f}"
        )

        print(
            f"median exact matches = "
            f"{statistics.median(exact_counts):.3f}"
        )

        print(
            f"unique exact recovery = "
            f"{unique}/{len(test)}"
        )

        print(
            f"true modular survival = "
            f"{true_mod_survival}/{len(test)}"
        )

        print(
            f"true exact survival = "
            f"{true_exact_survival}/{len(test)}"
        )

    # ------------------------------------------------------------------
    # Direct comparison H6 vs tower
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("6. TOWER GAIN OVER H6 / E1")
    print("=" * 78)

    for label, moduli in [
        (
            "C3",
            CYCLOTOMIC[:3],
        ),
        (
            "C5",
            CYCLOTOMIC[:5],
        ),
        (
            "C7",
            CYCLOTOMIC,
        ),
    ]:

        h6_counts = []
        tower_counts = []

        for target in test:

            h6 = modular_survivors(
                target,
                moduli,
                False,
            )

            tower = modular_survivors(
                target,
                moduli,
                True,
            )

            h6_counts.append(
                len(h6)
            )

            tower_counts.append(
                len(tower)
            )

        h6_mean = statistics.fmean(
            h6_counts
        )

        tower_mean = statistics.fmean(
            tower_counts
        )

        print(
            f"{label}: "
            f"H6 mean={h6_mean:.3f} "
            f"tower mean={tower_mean:.3f} "
            f"compression="
            f"{(
                h6_mean / tower_mean
                if tower_mean
                else float('inf')
            ):.3f}x"
        )

        print(
            f"{label}: "
            f"H6 unique={sum(x == 1 for x in h6_counts)}/{len(test)} "
            f"tower unique={sum(x == 1 for x in tower_counts)}/{len(test)}"
        )

    # ------------------------------------------------------------------
    # Sample exact diagnostics
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. SAMPLE TEST TARGET DIAGNOSTICS")
    print("=" * 78)

    for idx, target in enumerate(
        test,
        TRAIN_TARGETS + 1,
    ):

        survivors = modular_survivors(
            target,
            CYCLOTOMIC,
            True,
        )

        nums = true_numerators[idx]

        exact = [
            s
            for s in survivors
            if exact_tower_match(
                target,
                s,
                nums,
            )
        ]

        print(
            f"target {idx:3d}: "
            f"tower survivors={len(survivors):4d} "
            f"exact={len(exact):3d} "
            f"true_s={target.s} "
            f"true_in_mod={target.s in survivors} "
            f"true_in_exact={target.s in exact}"
        )

        if len(exact) <= PRINT_LIMIT:
            print(
                f"    exact candidates = {exact}"
            )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    h6_c7 = family_results[
        "H6_ONLY_C7"
    ]

    tower_c7 = family_results[
        "TOWER_C7"
    ]

    tower_c5 = family_results[
        "TOWER_C5"
    ]

    print(
        "HELD-OUT C7 COMPARISON"
    )

    print(
        f"  H6 mean survivors    = "
        f"{statistics.fmean(h6_c7.test_counts):.3f}"
    )

    print(
        f"  Tower mean survivors = "
        f"{statistics.fmean(tower_c7.test_counts):.3f}"
    )

    print(
        f"  H6 unique            = "
        f"{h6_c7.test_unique}/{len(test)}"
    )

    print(
        f"  Tower unique         = "
        f"{tower_c7.test_unique}/{len(test)}"
    )

    print()
    print(
        "HELD-OUT C5 COMPARISON"
    )

    print(
        f"  H6 mean survivors    = "
        f"{statistics.fmean(
            family_results['H6_ONLY_C5'].test_counts
        ):.3f}"
    )

    print(
        f"  Tower mean survivors = "
        f"{statistics.fmean(tower_c5.test_counts):.3f}"
    )

    print(
        f"  H6 unique            = "
        f"{family_results['H6_ONLY_C5'].test_unique}/{len(test)}"
    )

    print(
        f"  Tower unique         = "
        f"{tower_c5.test_unique}/{len(test)}"
    )

    print()
    print(
        "INTERPRETATION"
    )
    print(
        "--------------"
    )

    c7_tower_mean = statistics.fmean(
        tower_c7.test_counts
    )

    c7_h6_mean = statistics.fmean(
        h6_c7.test_counts
    )

    if (
        tower_c7.test_unique == len(test)
        and c7_tower_mean <= 1.0
    ):
        print(
            "STRONG TOWER IDENTIFIABILITY:"
        )
        print(
            "The full canonical tower uniquely identifies s "
            "on every held-out target."
        )

    elif (
        c7_tower_mean
        < c7_h6_mean
    ):
        print(
            "TOWER IMPROVEMENT:"
        )
        print(
            "Adding H8/H10/H12 gives a measurable reduction "
            "in the remaining s candidates."
        )

    else:
        print(
            "NO MATERIAL TOWER IMPROVEMENT:"
        )
        print(
            "The higher H_k levels add little beyond H6/E1 "
            "for this modulus family."
        )

    print()
    print(
        "CONTROL WARNING:"
    )
    print(
        "If generic control moduli also produce one survivor, "
        "uniqueness alone is not evidence of a cyclotomic advantage."
    )

    print()
    print(
        "CRITICAL NEXT QUESTION:"
    )
    print(
        "Can the H_k tower be computed from N alone without first "
        "knowing p and q?"
    )

    total_runtime = (
        time.perf_counter()
        - overall_start
    )

    print(
        f"\ntotal runtime = "
        f"{total_runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 85R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()