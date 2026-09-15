#!/usr/bin/env python3

"""
================================================================================
EXPERIMENT 697
HIGH-C SCANNER: DISCOVER c IN [1, 1_000_000] FROM TRAINING CASES
                  THEN VALIDATE ON LARGER INDEPENDENT CASES
================================================================================

FRAME A
    n = p*q
    C = q + 3

SHIFT CHANNEL
    H_c = gcd(C, n+c)

IDENTITY
    n+c = p*q+c
         = p*(q+3) + (c-3p)
         = p*C + (c-3p)

therefore
    H_c | C
    H_c | (c-3p)

DISCOVERY IDEA
--------------
Do NOT hard-code c = 3, 81, 137.

Instead scan every odd c from 1 through 1,000,000 and look for
c-values that repeatedly expose nontrivial components of C across
multiple independent training semiprimes.

VALIDATION
----------
The discovered c-values are then tested on independent, larger
semiprimes.

IMPORTANT
---------
This experiment does NOT factor n+c.

It uses the known training factorization only to discover/validate
the observable H-channel.

That keeps the experiment focused on the structure of c itself.
================================================================================
"""

from __future__ import annotations

import itertools
import math
import random
from collections import Counter
from dataclasses import dataclass
from typing import Iterable


# =============================================================================
# CONFIGURATION
# =============================================================================

C_MIN = 1
C_MAX = 1_000_000

# Only odd c values are scanned.
ODD_ONLY = True

# Training cases:
TRAINING_CASES = 3

# Training primes:
TRAIN_P_MIN = 1_000_001
TRAIN_P_MAX = 5_000_000

# Validation primes are deliberately larger.
VALIDATION_CASES = 5
VALID_P_MIN = 10_000_001
VALID_P_MAX = 50_000_000

RANDOM_SEED = 697_000_001

# Discovery rule.
#
# A c gets into the candidate list when:
#
#   * at least MIN_NONTRIVIAL_CASES training cases have H_c > 2
#   * the combined LCM of the H-values has an odd component > 1
#
# This removes the universal trivial factor 2.
MIN_NONTRIVIAL_CASES = 2

# Rank candidates by this score.
#
# score =
#     number of cases with H > 2
#     +
#     total odd-component bit contribution
#     +
#     number of distinct odd prime factors exposed
#
# The exact score is only for ranking; candidate inclusion uses
# MIN_NONTRIVIAL_CASES.
#
# Number of strongest candidates to display.
TOP_K = 100

# Combination search.
#
# Searching all combinations of millions of candidates would explode,
# so we take the strongest TOP_COMBO_POOL individual c-values first.
TOP_COMBO_POOL = 30

# Search pairs and triples.
COMBINATION_SIZES = (2, 3)

# A combination is considered interesting on a case when its LCM
# has a nontrivial odd component.
REQUIRE_ODD_LCM_COMPONENT = True


# =============================================================================
# PRIME GENERATION
# =============================================================================

# Deterministic Miller-Rabin for 64-bit integers.
MR_BASES = (
    2, 325, 9375, 28178,
    450775, 9780504, 1795265022
)


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    for a in MR_BASES:
        a %= n
        if a == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False

    return True


def random_prime(rng: random.Random, lo: int, hi: int) -> int:
    if lo < 2:
        lo = 2

    if lo > hi:
        raise ValueError("invalid prime range")

    while True:
        x = rng.randrange(lo, hi + 1)

        if x <= 2:
            x = 3

        if x % 2 == 0:
            x += 1

        while x <= hi:
            if is_probable_prime(x):
                return x

            x += 2


def make_semiprime(
    rng: random.Random,
    lo: int,
    hi: int,
) -> tuple[int, int, int]:
    while True:
        p = random_prime(rng, lo, hi)
        q = random_prime(rng, lo, hi)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        return p, q, n


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class Case:
    p: int
    q: int
    n: int
    C: int


@dataclass
class Candidate:
    c: int
    score: float
    lcm_H: int
    H_values: tuple[int, ...]
    odd_parts: tuple[int, ...]
    nontrivial_cases: int
    distinct_odd_primes: int


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def odd_part(x: int) -> int:
    while x % 2 == 0 and x > 0:
        x //= 2
    return x


def distinct_prime_factors(n: int) -> set[int]:
    """
    Trial-factorization used only on H-values.

    H divides C and C is only a few million in the training range,
    so this is cheap.
    """
    result: set[int] = set()

    while n % 2 == 0 and n > 1:
        result.add(2)
        n //= 2

    d = 3

    while d * d <= n:
        if n % d == 0:
            result.add(d)

            while n % d == 0:
                n //= d

        d += 2

    if n > 1:
        result.add(n)

    return result


def lcm(a: int, b: int) -> int:
    return a // math.gcd(a, b) * b


def lcm_many(values: Iterable[int]) -> int:
    out = 1

    for x in values:
        out = lcm(out, x)

    return out


# =============================================================================
# CASE GENERATION
# =============================================================================

def generate_training_cases(
    rng: random.Random,
) -> list[Case]:

    cases: list[Case] = []

    while len(cases) < TRAINING_CASES:
        p, q, n = make_semiprime(
            rng,
            TRAIN_P_MIN,
            TRAIN_P_MAX,
        )

        C = q + 3

        cases.append(
            Case(
                p=p,
                q=q,
                n=n,
                C=C,
            )
        )

    return cases


def generate_validation_cases(
    rng: random.Random,
) -> list[Case]:

    cases: list[Case] = []

    while len(cases) < VALIDATION_CASES:
        p, q, n = make_semiprime(
            rng,
            VALID_P_MIN,
            VALID_P_MAX,
        )

        C = q + 3

        cases.append(
            Case(
                p=p,
                q=q,
                n=n,
                C=C,
            )
        )

    return cases


# =============================================================================
# DISCOVERY
# =============================================================================

def scan_candidate_c_values(
    training: list[Case],
) -> list[Candidate]:

    candidates: list[Candidate] = []

    values = range(
        C_MIN,
        C_MAX + 1,
        2 if ODD_ONLY else 1,
    )

    total = C_MAX // 2 if ODD_ONLY else C_MAX

    print()
    print("=" * 80)
    print("SCANNING C VALUES")
    print("=" * 80)
    print(f"range={C_MIN}..{C_MAX}")
    print(f"odd_only={ODD_ONLY}")
    print(f"scan_count~={total}")

    progress_points = {
        C_MAX // 10 * i
        for i in range(1, 11)
    }

    for c in values:

        H_values = tuple(
            math.gcd(case.C, case.n + c)
            for case in training
        )

        odd_parts = tuple(
            odd_part(h)
            for h in H_values
        )

        nontrivial_cases = sum(
            1
            for h in H_values
            if h > 2
        )

        combined_lcm = lcm_many(H_values)

        combined_odd = odd_part(combined_lcm)

        if combined_odd <= 1:
            continue

        if nontrivial_cases < MIN_NONTRIVIAL_CASES:
            continue

        odd_prime_set: set[int] = set()

        for op in odd_parts:
            if op > 1:
                odd_prime_set.update(
                    p
                    for p in distinct_prime_factors(op)
                    if p != 2
                )

        # Ranking score.
        score = 0.0

        # Number of cases with nontrivial information.
        score += 10.0 * nontrivial_cases

        # Prefer larger combined odd component.
        score += math.log2(combined_odd)

        # Prefer multiple different odd primes.
        score += 5.0 * len(odd_prime_set)

        candidates.append(
            Candidate(
                c=c,
                score=score,
                lcm_H=combined_lcm,
                H_values=H_values,
                odd_parts=odd_parts,
                nontrivial_cases=nontrivial_cases,
                distinct_odd_primes=len(odd_prime_set),
            )
        )

        if c in progress_points:
            print(f"progress c={c:,}")

    candidates.sort(
        key=lambda x: (
            x.score,
            x.lcm_H,
            x.distinct_odd_primes,
        ),
        reverse=True,
    )

    return candidates


# =============================================================================
# VALIDATION OF INDIVIDUAL C
# =============================================================================

def validate_candidate(
    candidate: Candidate,
    validation_cases: list[Case],
) -> dict:

    H_values = []

    nontrivial = 0
    odd_lcm = 1

    for case in validation_cases:
        H = math.gcd(case.C, case.n + candidate.c)

        H_values.append(H)

        if H > 2:
            nontrivial += 1

        odd_lcm = lcm(
            odd_lcm,
            odd_part(H),
        )

    return {
        "c": candidate.c,
        "nontrivial_cases": nontrivial,
        "odd_lcm": odd_lcm,
        "H_values": tuple(H_values),
        "survives_all": (
            odd_lcm > 1
            and nontrivial > 0
        ),
    }


# =============================================================================
# COMBINATION TEST
# =============================================================================

def evaluate_combination(
    c_values: tuple[int, ...],
    cases: list[Case],
) -> dict:

    case_results = []

    for case in cases:

        H_values = tuple(
            math.gcd(
                case.C,
                case.n + c,
            )
            for c in c_values
        )

        combined = lcm_many(H_values)
        combined_odd = odd_part(combined)

        case_results.append(
            {
                "H_values": H_values,
                "LCM": combined,
                "odd_LCM": combined_odd,
                "nontrivial": combined_odd > 1,
            }
        )

    successful = sum(
        1
        for result in case_results
        if result["nontrivial"]
    )

    return {
        "c_values": c_values,
        "successes": successful,
        "total": len(cases),
        "ratio": successful / len(cases),
        "case_results": case_results,
    }


# =============================================================================
# OUTPUT
# =============================================================================

def print_training_cases(cases: list[Case]) -> None:

    print()
    print("=" * 80)
    print("TRAINING CASES")
    print("=" * 80)

    for i, case in enumerate(cases, 1):
        print(
            f"case={i} "
            f"p={case.p} "
            f"q={case.q} "
            f"n={case.n} "
            f"C={case.C}"
        )


def print_validation_cases(cases: list[Case]) -> None:

    print()
    print("=" * 80)
    print("VALIDATION CASES")
    print("=" * 80)

    for i, case in enumerate(cases, 1):
        print(
            f"case={i} "
            f"p={case.p} "
            f"q={case.q} "
            f"n={case.n} "
            f"C={case.C}"
        )


def print_top_candidates(
    candidates: list[Candidate],
) -> None:

    print()
    print("=" * 80)
    print("TOP DISCOVERED C VALUES")
    print("=" * 80)

    print(
        f"{'rank':>5} "
        f"{'c':>10} "
        f"{'score':>10} "
        f"{'H-vector':>35} "
        f"{'LCM(H)':>15} "
        f"{'odd LCM':>15}"
    )

    for rank, candidate in enumerate(
        candidates[:TOP_K],
        1,
    ):
        print(
            f"{rank:5d} "
            f"{candidate.c:10d} "
            f"{candidate.score:10.3f} "
            f"{str(candidate.H_values):>35} "
            f"{candidate.lcm_H:15d} "
            f"{odd_part(candidate.lcm_H):15d}"
        )


def print_validation_results(
    candidates: list[Candidate],
    validation_cases: list[Case],
) -> list[dict]:

    print()
    print("=" * 80)
    print("VALIDATION OF DISCOVERED C VALUES")
    print("=" * 80)

    results = []

    for candidate in candidates:
        result = validate_candidate(
            candidate,
            validation_cases,
        )

        results.append(result)

    results.sort(
        key=lambda x: (
            x["survives_all"],
            x["nontrivial_cases"],
            x["odd_lcm"],
        ),
        reverse=True,
    )

    print(
        f"{'c':>10} "
        f"{'nontrivial':>12} "
        f"{'odd_LCM':>15} "
        f"{'survives':>10}"
    )

    for result in results[:TOP_K]:

        print(
            f"{result['c']:10d} "
            f"{result['nontrivial_cases']:12d} "
            f"{result['odd_lcm']:15d} "
            f"{str(result['survives_all']):>10}"
        )

    return results


def print_survivors(
    results: list[dict],
) -> list[int]:

    survivors = [
        result["c"]
        for result in results
        if result["survives_all"]
    ]

    print()
    print("=" * 80)
    print("GENERALIZING C VALUES")
    print("=" * 80)

    print(
        f"survivors={len(survivors)}"
    )

    if survivors:
        print(survivors)

    return survivors


# =============================================================================
# COMBINATION SEARCH
# =============================================================================

def search_combinations(
    candidates: list[Candidate],
    training: list[Case],
    validation: list[Case],
) -> None:

    pool = [
        candidate.c
        for candidate in candidates[:TOP_COMBO_POOL]
    ]

    print()
    print("=" * 80)
    print("COMBINATION SEARCH")
    print("=" * 80)

    print(
        f"combination_pool={pool}"
    )

    all_results = []

    for size in COMBINATION_SIZES:

        print()
        print(f"COMBINATIONS OF SIZE {size}")
        print("-" * 80)

        for combo in itertools.combinations(pool, size):

            train_result = evaluate_combination(
                combo,
                training,
            )

            validation_result = evaluate_combination(
                combo,
                validation,
            )

            all_results.append(
                (
                    combo,
                    train_result,
                    validation_result,
                )
            )

    # Sort by validation performance first.
    all_results.sort(
        key=lambda item: (
            item[2]["ratio"],
            max(
                r["odd_LCM"]
                for r in item[2]["case_results"]
            ),
        ),
        reverse=True,
    )

    print()
    print("TOP COMBINATIONS")
    print("-" * 80)

    for combo, train_result, validation_result in all_results[:50]:

        print(
            f"c={combo} "
            f"training="
            f"{train_result['successes']}/"
            f"{train_result['total']} "
            f"validation="
            f"{validation_result['successes']}/"
            f"{validation_result['total']}"
        )

        train_lcms = [
            r["odd_LCM"]
            for r in train_result["case_results"]
        ]

        validation_lcms = [
            r["odd_LCM"]
            for r in validation_result["case_results"]
        ]

        print(
            f"    training odd LCMs={train_lcms}"
        )

        print(
            f"    validation odd LCMs={validation_lcms}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    rng = random.Random(RANDOM_SEED)

    print("=" * 80)
    print("EXPERIMENT 697 START")
    print("=" * 80)

    print()
    print(f"C range       = [{C_MIN}, {C_MAX}]")
    print(f"odd only      = {ODD_ONLY}")
    print(f"training      = {TRAINING_CASES}")
    print(f"validation    = {VALIDATION_CASES}")
    print(
        f"training p,q  = "
        f"{TRAIN_P_MIN:,} .. {TRAIN_P_MAX:,}"
    )
    print(
        f"validation p,q = "
        f"{VALID_P_MIN:,} .. {VALID_P_MAX:,}"
    )

    # -------------------------------------------------------------------------
    # TEST 0
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("TEST 0: GENERATE TRAINING CASES")
    print("=" * 80)

    training = generate_training_cases(rng)

    print_training_cases(training)

    # -------------------------------------------------------------------------
    # TEST 1
    # -------------------------------------------------------------------------

    candidates = scan_candidate_c_values(training)

    print()
    print("=" * 80)
    print("DISCOVERY SUMMARY")
    print("=" * 80)

    print(
        f"discovered candidates={len(candidates)}"
    )

    print_top_candidates(candidates)

    # -------------------------------------------------------------------------
    # TEST 2
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("TEST 2: GENERATE LARGER VALIDATION CASES")
    print("=" * 80)

    validation = generate_validation_cases(rng)

    print_validation_cases(validation)

    # -------------------------------------------------------------------------
    # TEST 3
    # -------------------------------------------------------------------------

    results = print_validation_results(
        candidates,
        validation,
    )

    survivors = print_survivors(results)

    # -------------------------------------------------------------------------
    # TEST 4
    # -------------------------------------------------------------------------

    search_combinations(
        candidates,
        training,
        validation,
    )

    # -------------------------------------------------------------------------
    # REPRESENTATIVE SURVIVORS
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("REPRESENTATIVE GENERALIZING C VALUES")
    print("=" * 80)

    for c in survivors[:20]:

        candidate = next(
            x for x in candidates
            if x.c == c
        )

        print()
        print(f"c={c}")

        print(
            f"    training H={candidate.H_values}"
        )

        validation_H = tuple(
            math.gcd(
                case.C,
                case.n + c,
            )
            for case in validation
        )

        print(
            f"    validation H={validation_H}"
        )

        print(
            f"    training LCM={candidate.lcm_H}"
        )

        print(
            f"    training odd LCM="
            f"{odd_part(candidate.lcm_H)}"
        )

        validation_lcm = lcm_many(
            validation_H
        )

        print(
            f"    validation LCM={validation_lcm}"
        )

        print(
            f"    validation odd LCM="
            f"{odd_part(validation_lcm)}"
        )

    # -------------------------------------------------------------------------
    # FINAL SUMMARY
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 80)

    print(
        """
The experiment searches for c-values rather than assuming
c = 3, 81, 137.

For every training case:

    C = q + 3
    H_c = gcd(C, n+c)

The scan covers every odd:

    1 <= c <= 1,000,000

A c is retained when it exposes nontrivial odd structure
in at least the configured number of training cases.

The independent validation stage then asks:

    Does the same c continue to expose nontrivial structure
    for completely new, larger semiprimes?

This gives three levels:

    DISCOVERY
        c repeatedly looks interesting on training cases

    GENERALIZATION
        the same c remains interesting on unseen cases

    COMBINATION
        several c-values jointly expose more structure than
        the individual channels

The experiment does NOT assume a special high-c formula.

It searches the interval directly and lets the training cases
select the candidate shifts.

The most interesting output is therefore not simply a large list
of c-values, but the subset that survives independent validation.
"""
    )

    print()
    print(
        f"training candidates = {len(candidates)}"
    )

    print(
        f"generalizing c      = {len(survivors)}"
    )

    print("=" * 80)
    print("EXPERIMENT 697 FINISHED")
    print("=" * 80)


if __name__ == "__main__":
    main()