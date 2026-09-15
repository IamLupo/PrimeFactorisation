#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 86
INCREMENTAL EISENSTEIN TOWER INFORMATION TEST

H6
H6 + H8
H6 + H8 + H10
H6 + H8 + H10 + H12

MARGINAL CANDIDATE COLLAPSE
NESTED MODULAR SIGNATURES
STRICT TARGET HOLDOUT
CYCLOTOMIC VS CONTROL
ORACLE INFORMATION ANALYSIS

NO CSV OUTPUT
NO SKLEARN
==============================================================================

MAIN QUESTION
-------------
Experiment 85R showed:

    H6-only C5     mean survivors ~= 58.3
    full tower C5  mean survivors ~= 9.7

This experiment asks WHICH LEVEL supplies the information.

For each modulus family we measure:

    H6
    H6+H8
    H6+H8+H10
    H6+H8+H10+H12

and compute:

    - mean survivors
    - median survivors
    - unique recovery
    - true-s survival
    - marginal reduction at each level
    - candidate-elimination fraction
    - exact oracle recovery after modular filtering

The true p,q are used only to construct the oracle tower values.
Candidate reconstruction uses n and candidate s only.

==============================================================================
"""

from __future__ import annotations

import math
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

NUM_TARGETS = 40
TRAIN_TARGETS = 30

P_MIN = 2_000_000
P_MAX = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_400_000

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67
]

CONTROL = [
    673, 4561, 4759, 6211,
    7879, 7951, 8689
]

RNG_SEED = 86086

LEVELS = [
    "H6",
    "H6+H8",
    "H6+H8+H10",
    "H6+H8+H10+H12",
]

# Maximum detailed candidate list per target.
PRINT_LIMIT = 20


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


TowerNums = Tuple[int, int, int, int]


@dataclass
class LevelStats:
    level: str
    counts: List[int]
    unique: int
    true_survival: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(
    lo: int,
    hi: int,
) -> List[int]:

    lo = max(2, lo)

    if hi < lo:
        return []

    sieve = bytearray(
        b"\x01" * (hi + 1)
    )

    sieve[0:2] = b"\x00\x00"

    for p in range(
        2,
        math.isqrt(hi) + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        sieve[
            start:
            hi + 1:
            p
        ] = b"\x00" * (
            ((hi - start) // p) + 1
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

        if (p, q) in seen:
            continue

        seen.add((p, q))

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
# POWER SUMS / DIVISOR SUMS
# ============================================================================

def power_sums(
    n: int,
    s: int,
    max_power: int,
) -> List[int]:

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


def sigma_j_semiprime(
    n: int,
    s: int,
    j: int,
) -> int:

    R = power_sums(
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
# PAPER TOWER NUMERATORS
# ============================================================================

def h6_num(
    n: int,
    s: int,
) -> int:

    sig1 = sigma_j_semiprime(
        n,
        s,
        1,
    )

    sig3 = sigma_j_semiprime(
        n,
        s,
        3,
    )

    return (
        (n * n - n + 1) * sig1
        - sig3
    )


def h_even_num(
    n: int,
    s: int,
    k: int,
) -> int:

    if k < 8 or k % 2:
        raise ValueError(
            "k must be even and >= 8"
        )

    a = k - 7
    b = k - 5
    c = k - 3

    sa = sigma_j_semiprime(
        n,
        s,
        a,
    )

    sb = sigma_j_semiprime(
        n,
        s,
        b,
    )

    sc = sigma_j_semiprime(
        n,
        s,
        c,
    )

    return (
        -n * n * sa
        + (n * n + 1) * sb
        - sc
    )


def tower_nums(
    n: int,
    s: int,
) -> TowerNums:

    return (
        h6_num(n, s),
        h_even_num(n, s, 8),
        h_even_num(n, s, 10),
        h_even_num(n, s, 12),
    )


# ============================================================================
# TARGET VALIDATION
# ============================================================================

def validate_target(
    t: Target,
) -> bool:

    nums = tower_nums(
        t.n,
        t.s,
    )

    if nums[0] % 6 != 0:
        return False

    if nums[1] % 24 != 0:
        return False

    if nums[2] % 24 != 0:
        return False

    if nums[3] % 24 != 0:
        return False

    # Direct sigma checks against p,q.
    for j in [
        1, 3, 5, 7, 9
    ]:
        indirect = sigma_j_semiprime(
            t.n,
            t.s,
            j,
        )

        direct = (
            (1 + t.p ** j)
            * (1 + t.q ** j)
        )

        if indirect != direct:
            return False

    return True


# ============================================================================
# MODULAR TOWER SIGNATURE
# ============================================================================

def tower_mod(
    n_mod: int,
    s_mod: int,
    ell: int,
) -> TowerNums:

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

    H6 = (
        (
            (n * n - n + 1)
            * sigma(1)
            - sigma(3)
        )
        * inv6
    ) % ell

    H8 = (
        (
            -n * n * sigma(1)
            + (n * n + 1) * sigma(3)
            - sigma(5)
        )
        * inv24
    ) % ell

    H10 = (
        (
            -n * n * sigma(3)
            + (n * n + 1) * sigma(5)
            - sigma(7)
        )
        * inv24
    ) % ell

    H12 = (
        (
            -n * n * sigma(5)
            + (n * n + 1) * sigma(7)
            - sigma(9)
        )
        * inv24
    ) % ell

    return (
        H6,
        H8,
        H10,
        H12,
    )


# ============================================================================
# SIGNATURE PREFIX
# ============================================================================

def signature_prefix(
    signature: TowerNums,
    level: str,
) -> Tuple[int, ...]:

    if level == "H6":
        return (
            signature[0],
        )

    if level == "H6+H8":
        return (
            signature[0],
            signature[1],
        )

    if level == "H6+H8+H10":
        return (
            signature[0],
            signature[1],
            signature[2],
        )

    if level == "H6+H8+H10+H12":
        return signature

    raise ValueError(
        f"Unknown level: {level}"
    )


# ============================================================================
# ALLOWED RESIDUES
# ============================================================================

def allowed_residues(
    t: Target,
    ell: int,
    level: str,
) -> bytearray:

    true_signature = tower_mod(
        t.n,
        t.s,
        ell,
    )

    wanted = signature_prefix(
        true_signature,
        level,
    )

    table = bytearray(
        b"\x00" * ell
    )

    n_mod = t.n % ell

    for r in range(ell):

        candidate = tower_mod(
            n_mod,
            r,
            ell,
        )

        if (
            signature_prefix(
                candidate,
                level,
            )
            == wanted
        ):
            table[r] = 1

    return table


# ============================================================================
# DOMAIN
# ============================================================================

def even_domain() -> Iterable[int]:

    start = (
        S_MIN
        if S_MIN % 2 == 0
        else S_MIN + 1
    )

    for s in range(
        start,
        S_MAX + 1,
        2,
    ):
        yield s


# ============================================================================
# NESTED CANDIDATE COLLAPSE
# ============================================================================

def nested_survivors(
    t: Target,
    moduli: Sequence[int],
) -> Dict[str, List[int]]:

    tables = {
        level: [
            (
                ell,
                allowed_residues(
                    t,
                    ell,
                    level,
                )
            )
            for ell in moduli
        ]
        for level in LEVELS
    }

    survivors = {
        level: []
        for level in LEVELS
    }

    for s in even_domain():

        # Evaluate each level only until it fails.
        active = {
            level: True
            for level in LEVELS
        }

        for ell, _dummy in []:
            pass

        for level in LEVELS:

            if not active[level]:
                continue

            for ell, table in tables[level]:

                if table[s % ell] == 0:

                    active[level] = False
                    break

            if active[level]:
                survivors[level].append(s)

    return survivors


# ============================================================================
# MORE EFFICIENT NESTED VERSION
# ============================================================================

def nested_survivors_fast(
    t: Target,
    moduli: Sequence[int],
) -> Dict[str, List[int]]:

    """
    Build one table per modulus containing the four signatures.

    Then scan s once, updating nested masks.

    This avoids reconstructing the residue tables for every level.
    """

    true_signatures = {
        ell: tower_mod(
            t.n,
            t.s,
            ell,
        )
        for ell in moduli
    }

    tables = {}

    for ell in moduli:

        residue_map = []

        n_mod = t.n % ell

        wanted = true_signatures[ell]

        for r in range(ell):

            residue_map.append(
                tower_mod(
                    n_mod,
                    r,
                    ell,
                )
            )

        tables[ell] = (
            residue_map,
            wanted,
        )

    out = {
        level: []
        for level in LEVELS
    }

    for s in even_domain():

        states = [True] * len(
            LEVELS
        )

        for ell in moduli:

            residue_map, wanted = tables[ell]

            cand = residue_map[
                s % ell
            ]

            for i, level in enumerate(
                LEVELS
            ):

                if not states[i]:
                    continue

                prefix_len = i + 1

                if (
                    cand[:prefix_len]
                    != wanted[:prefix_len]
                ):
                    states[i] = False

        for i, level in enumerate(
            LEVELS
        ):

            if states[i]:
                out[level].append(s)

    return out


# ============================================================================
# STATS
# ============================================================================

def stats_for(
    counts: Sequence[int],
) -> Tuple[float, float]:

    return (
        statistics.fmean(counts),
        statistics.median(counts),
    )


# ============================================================================
# MARGINAL INFORMATION
# ============================================================================

def marginal_fraction(
    previous: Sequence[int],
    current: Sequence[int],
) -> float:

    """
    Fraction of candidates removed by adding the new tower level.
    """

    if not previous:
        return float("nan")

    total_removed = sum(
        max(
            0,
            a - b,
        )
        for a, b in zip(
            previous,
            current,
        )
    )

    total_previous = sum(
        previous
    )

    if total_previous == 0:
        return 0.0

    return (
        total_removed
        / total_previous
    )


def geometric_compression(
    counts: Sequence[int],
) -> float:

    if not counts:
        return 0.0

    logs = [
        math.log(
            max(1, x)
        )
        for x in counts
    ]

    return math.exp(
        statistics.fmean(logs)
    )


# ============================================================================
# TARGET ANALYSIS
# ============================================================================

def analyze_family(
    name: str,
    moduli: Sequence[int],
    train: Sequence[Target],
    test: Sequence[Target],
) -> Dict[str, LevelStats]:

    print("\n")
    print("=" * 78)
    print(name)
    print(
        f"moduli = {list(moduli)}"
    )
    print("=" * 78)

    train_data = {
        level: []
        for level in LEVELS
    }

    test_data = {
        level: []
        for level in LEVELS
    }

    train_unique = {
        level: 0
        for level in LEVELS
    }

    test_unique = {
        level: 0
        for level in LEVELS
    }

    train_survive = {
        level: 0
        for level in LEVELS
    }

    test_survive = {
        level: 0
        for level in LEVELS
    }

    t0 = time.perf_counter()

    for target in train:

        result = nested_survivors_fast(
            target,
            moduli,
        )

        for level in LEVELS:

            c = len(
                result[level]
            )

            train_data[level].append(
                c
            )

            if c == 1:
                train_unique[level] += 1

            if target.s in result[level]:
                train_survive[level] += 1

    for target in test:

        result = nested_survivors_fast(
            target,
            moduli,
        )

        for level in LEVELS:

            c = len(
                result[level]
            )

            test_data[level].append(
                c
            )

            if c == 1:
                test_unique[level] += 1

            if target.s in result[level]:
                test_survive[level] += 1

    elapsed = (
        time.perf_counter()
        - t0
    )

    print(
        f"runtime = {elapsed:.6f}s"
    )

    stats = {}

    for level in LEVELS:

        train_mean, train_median = stats_for(
            train_data[level]
        )

        test_mean, test_median = stats_for(
            test_data[level]
        )

        print("\n")
        print(level)
        print("-" * 78)

        print(
            f"TRAIN mean={train_mean:.3f} "
            f"median={train_median:.3f} "
            f"unique={train_unique[level]}/{len(train)} "
            f"survival={train_survive[level]}/{len(train)}"
        )

        print(
            f"TEST  mean={test_mean:.3f} "
            f"median={test_median:.3f} "
            f"unique={test_unique[level]}/{len(test)} "
            f"survival={test_survive[level]}/{len(test)}"
        )

        print(
            f"TEST geometric mean={geometric_compression(test_data[level]):.3f}"
        )

        stats[level] = LevelStats(
            level=level,
            counts=test_data[level],
            unique=test_unique[level],
            true_survival=test_survive[level],
        )

    # --------------------------------------------------------------
    # Marginal information
    # --------------------------------------------------------------

    print("\n")
    print("MARGINAL COLLAPSE")
    print("-" * 78)

    for i in range(
        1,
        len(LEVELS),
    ):

        previous = (
            stats[LEVELS[i - 1]]
            .counts
        )

        current = (
            stats[LEVELS[i]]
            .counts
        )

        fraction = marginal_fraction(
            previous,
            current,
        )

        print(
            f"{LEVELS[i-1]} "
            f"-> {LEVELS[i]}: "
            f"candidate reduction="
            f"{fraction:.6f}"
        )

    return stats


# ============================================================================
# EXACT ORACLE FILTER
# ============================================================================

def exact_tower_match(
    t: Target,
    candidate_s: int,
    true_nums: TowerNums,
) -> bool:

    """
    Numerator-safe exact comparison.

    Never require a false candidate to have integral H_k values.
    """

    return (
        tower_nums(
            t.n,
            candidate_s,
        )
        == true_nums
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 86")
    print("INCREMENTAL EISENSTEIN TOWER INFORMATION TEST")
    print("H6 -> H6+H8 -> H6+H8+H10 -> H6+H8+H10+H12")
    print("MARGINAL CANDIDATE COLLAPSE")
    print("STRICT TARGET HOLDOUT")
    print("CYCLOTOMIC VS CONTROL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Population
    # ------------------------------------------------------------------

    t0 = time.perf_counter()

    primes = sieve_primes(
        P_MIN,
        P_MAX,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = "
        f"{len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    for i, t in enumerate(
        targets[:24],
        1,
    ):

        print(
            f"target {i:3d}: "
            f"p={t.p} "
            f"q={t.q} "
            f"n={t.n} "
            f"s={t.s}"
        )

    if NUM_TARGETS > 24:
        print(
            "... remaining generated targets omitted"
        )

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    print("\n2. TOWER VALIDATION")
    print("-" * 78)

    failures = 0

    true_nums_by_id = {}

    for i, t in enumerate(
        targets,
        1,
    ):

        ok = validate_target(t)

        nums = tower_nums(
            t.n,
            t.s,
        )

        true_nums_by_id[i] = nums

        if not ok:
            failures += 1

        if i <= 24:
            print(
                f"target {i:3d}: "
                f"identity_ok={ok}"
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
    # Holdout
    # ------------------------------------------------------------------

    train = targets[:TRAIN_TARGETS]
    test = targets[TRAIN_TARGETS:]

    print("\n3. TARGET HOLDOUT")
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
    # Main experiments
    # ------------------------------------------------------------------

    cyclo_results = analyze_family(
        "CYCLOTOMIC C7",
        CYCLOTOMIC,
        train,
        test,
    )

    control_results = analyze_family(
        "CONTROL R7",
        CONTROL,
        train,
        test,
    )

    # ------------------------------------------------------------------
    # Print direct marginal comparison
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("4. CYCLOTOMIC VS CONTROL MARGINAL INFORMATION")
    print("=" * 78)

    print(
        "level                        "
        "C-test mean    R-test mean    C/R"
    )

    for level in LEVELS:

        c_mean = statistics.fmean(
            cyclo_results[level].counts
        )

        r_mean = statistics.fmean(
            control_results[level].counts
        )

        ratio = (
            c_mean / r_mean
            if r_mean
            else float("inf")
        )

        print(
            f"{level:28s} "
            f"{c_mean:12.3f} "
            f"{r_mean:12.3f} "
            f"{ratio:8.3f}"
        )

    # ------------------------------------------------------------------
    # Exact oracle check on test targets
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("5. EXACT ORACLE CHECK ON TEST TARGETS")
    print("=" * 78)

    for level in [
        "H6",
        "H6+H8",
        "H6+H8+H10",
        "H6+H8+H10+H12",
    ]:

        survivors_total = []
        exact_total = []

        unique_exact = 0
        true_survival = 0

        for i, target in enumerate(
            test,
            TRAIN_TARGETS + 1,
        ):

            result = nested_survivors_fast(
                target,
                CYCLOTOMIC,
            )

            survivors = result[level]

            true_num = (
                true_nums_by_id[i]
            )

            exact = [
                s
                for s in survivors
                if exact_tower_match(
                    target,
                    s,
                    true_num,
                )
            ]

            survivors_total.append(
                len(survivors)
            )

            exact_total.append(
                len(exact)
            )

            if target.s in survivors:
                true_survival += 1

            if len(exact) == 1:
                unique_exact += 1

        print("\n")
        print(level)
        print(
            f"mean modular survivors = "
            f"{statistics.fmean(survivors_total):.3f}"
        )
        print(
            f"mean exact matches = "
            f"{statistics.fmean(exact_total):.3f}"
        )
        print(
            f"unique exact recovery = "
            f"{unique_exact}/{len(test)}"
        )
        print(
            f"true modular survival = "
            f"{true_survival}/{len(test)}"
        )

    # ------------------------------------------------------------------
    # Sample target
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("6. SAMPLE HELD-OUT TARGETS")
    print("=" * 78)

    for i, target in enumerate(
        test,
        TRAIN_TARGETS + 1,
    ):

        result = nested_survivors_fast(
            target,
            CYCLOTOMIC,
        )

        print(
            f"target {i:3d} "
            f"true_s={target.s}"
        )

        for level in LEVELS:

            vals = result[level]

            print(
                f"  {level:24s} "
                f"survivors={len(vals):5d}"
            )

            if len(vals) <= PRINT_LIMIT:

                print(
                    f"    candidates={vals}"
                )

    # ------------------------------------------------------------------
    # Interpretation
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The experiment asks WHICH H_k actually contributes "
        "the candidate-space collapse."
    )

    c_h6 = statistics.fmean(
        cyclo_results["H6"].counts
    )

    c_h8 = statistics.fmean(
        cyclo_results["H6+H8"].counts
    )

    c_h10 = statistics.fmean(
        cyclo_results["H6+H8+H10"].counts
    )

    c_h12 = statistics.fmean(
        cyclo_results["H6+H8+H10+H12"].counts
    )

    print()
    print(
        f"CYCLOTOMIC TEST MEANS:"
    )
    print(
        f"  H6                  = {c_h6:.3f}"
    )
    print(
        f"  H6 + H8             = {c_h8:.3f}"
    )
    print(
        f"  H6 + H8 + H10       = {c_h10:.3f}"
    )
    print(
        f"  H6 + H8 + H10 + H12 = {c_h12:.3f}"
    )

    print()
    print(
        "MARGINAL REDUCTIONS:"
    )

    print(
        f"  H8 contribution  = "
        f"{(
            1 - c_h8 / c_h6
            if c_h6
            else 0
        ):.6f}"
    )

    print(
        f"  H10 contribution = "
        f"{(
            1 - c_h10 / c_h8
            if c_h8
            else 0
        ):.6f}"
    )

    print(
        f"  H12 contribution = "
        f"{(
            1 - c_h12 / c_h10
            if c_h10
            else 0
        ):.6f}"
    )

    print()

    if c_h12 < c_h10 < c_h8 < c_h6:
        print(
            "NESTED INFORMATION CONFIRMED:"
        )
        print(
            "Every added canonical H_k level contributes "
            "additional candidate-space compression."
        )

    elif c_h12 < c_h6:
        print(
            "TOWER INFORMATION CONFIRMED:"
        )
        print(
            "The higher tower improves identifiability, "
            "but the marginal contribution is not monotone."
        )

    else:
        print(
            "NO MATERIAL NESTED TOWER GAIN:"
        )
        print(
            "The higher H_k coefficients do not significantly "
            "improve the H6/E1 fingerprint."
        )

    print()
    print(
        "MOST IMPORTANT NEXT QUESTION:"
    )
    print(
        "If H8/H10/H12 add independent information, can their "
        "divisor-sum expressions be algebraically eliminated or "
        "compressed into an N-only relation?"
    )

    runtime = (
        time.perf_counter()
        - start
    )

    print(
        f"\ntotal runtime = "
        f"{runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 86 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

