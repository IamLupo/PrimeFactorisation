#!/usr/bin/env python3

"""
EXPERIMENT 188
Candidate k-dependent convolution reconstruction.

Purpose
-------
Experiment 187 established that the current exact summand is missing
effective k-dependence.

This experiment compares several structurally plausible k-dependent
alternating convolutions.

IMPORTANT
---------
This is a STRUCTURAL AUDIT, not a claim that any candidate is the true
mathematical definition.

Rules:
    - No results.txt or previous experiment output is read.
    - All test data are generated here.
    - Every function is defined here.
    - Each candidate is evaluated independently.
    - Support, parity, symmetry, zero regions, and boundary behavior
      are measured separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from typing import Callable


# ============================================================================
# BASIC ARITHMETIC
# ============================================================================

def C(n: int, r: int) -> int:
    if n < 0:
        return 0
    if r < 0 or r > n:
        return 0
    return comb(n, r)


def sgn(n: int) -> int:
    return -1 if n % 2 else 1


# ============================================================================
# SUPPORT
# ============================================================================

def support_lower(k: int, ell: int) -> int:
    return (ell - k) // 2


def support_upper(k: int, ell: int) -> int:
    return (ell + k) // 2


def in_support(k: int, ell: int, s: int) -> bool:
    return (
        support_lower(k, ell)
        <= s
        <= support_upper(k, ell)
    )


# ============================================================================
# CANDIDATE DEFINITION
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    name: str
    summand: Callable[[int, int, int, int], int]
    upper_limit: Callable[[int, int, int], int]


# ============================================================================
# COMMON UPPER LIMITS
# ============================================================================

def upper_ell(
    k: int,
    ell: int,
    s: int,
) -> int:
    return ell


def upper_k(
    k: int,
    ell: int,
    s: int,
) -> int:
    return min(k, ell)


def upper_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return min(s, ell)


# ============================================================================
# CANDIDATE A
# ============================================================================

def summand_A(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k enters the second upper binomial argument.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell + k - j, s - j)
    )


# ============================================================================
# CANDIDATE B
# ============================================================================

def summand_B(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k shifts the target index.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - k - j)
    )


# ============================================================================
# CANDIDATE C
# ============================================================================

def summand_C(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k enters by truncating the alternating sum.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - j)
    )


def value_C(
    k: int,
    ell: int,
    s: int,
) -> int:
    if not in_support(k, ell, s):
        return 0

    total = 0

    for j in range(
        upper_k(k, ell, s) + 1
    ):
        total += summand_C(
            k,
            ell,
            s,
            j,
        )

    return total


# ============================================================================
# CANDIDATE D
# ============================================================================

def summand_D(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k enters both the target index and second upper argument.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell + k - j, s - k - j)
    )


# ============================================================================
# CANDIDATE E
# ============================================================================

def summand_E(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    Symmetric-looking variant with ell-k.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - k - j, s - j)
    )


# ============================================================================
# CANDIDATE F
# ============================================================================

def summand_F(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    Variant where k changes the first upper parameter.
    """
    return (
        sgn(j)
        * C(ell + k, j)
        * C(ell + k - j, s - j)
    )


# ============================================================================
# GENERIC EVALUATION
# ============================================================================

def evaluate_candidate(
    candidate: Candidate,
    k: int,
    ell: int,
    s: int,
) -> int:
    if not in_support(k, ell, s):
        return 0

    total = 0

    upper = candidate.upper_limit(
        k,
        ell,
        s,
    )

    for j in range(upper + 1):
        total += candidate.summand(
            k,
            ell,
            s,
            j,
        )

    return total


def make_candidates() -> list[Candidate]:
    return [
        Candidate(
            "A: second upper +k",
            summand_A,
            upper_ell,
        ),
        Candidate(
            "B: target index -k",
            summand_B,
            upper_ell,
        ),
        Candidate(
            "C: truncated j<=k",
            summand_C,
            upper_k,
        ),
        Candidate(
            "D: second upper +k, target -k",
            summand_D,
            upper_ell,
        ),
        Candidate(
            "E: first upper -k",
            summand_E,
            upper_ell,
        ),
        Candidate(
            "F: first upper +k",
            summand_F,
            upper_ell,
        ),
    ]


# ============================================================================
# GENERATED CASES
# ============================================================================

@dataclass(frozen=True)
class Case:
    k: int
    ell: int
    s: int
    location: str


def generate_cases() -> list[Case]:
    cases: list[Case] = []

    for k in [1, 3, 5, 7, 9]:
        for ell in range(
            8,
            41,
        ):
            lo = support_lower(
                k,
                ell,
            )

            hi = support_upper(
                k,
                ell,
            )

            # Full support.
            for s in range(
                lo,
                hi + 1,
            ):
                if s == lo:
                    location = "LOWER"
                elif s == hi:
                    location = "UPPER"
                elif s == lo + 1:
                    location = "LOWER+1"
                elif s == hi - 1:
                    location = "UPPER-1"
                else:
                    location = "INTERIOR"

                cases.append(
                    Case(
                        k=k,
                        ell=ell,
                        s=s,
                        location=location,
                    )
                )

    return cases


# ============================================================================
# BASIC STATISTICS
# ============================================================================

def count_nonzero(
    candidate: Candidate,
    cases: list[Case],
) -> int:
    return sum(
        evaluate_candidate(
            candidate,
            case.k,
            case.ell,
            case.s,
        ) != 0
        for case in cases
    )


def count_k_variation(
    candidate: Candidate,
) -> tuple[int, int]:
    """
    For fixed (ell,s), determine whether candidate changes with k.
    """
    groups: dict[
        tuple[int, int],
        list[int],
    ] = {}

    for ell in range(8, 41):
        for s in range(0, ell + 10):
            values = []

            for k in [1, 3, 5, 7, 9]:
                if in_support(
                    k,
                    ell,
                    s,
                ):
                    values.append(
                        evaluate_candidate(
                            candidate,
                            k,
                            ell,
                            s,
                        )
                    )

            if values:
                groups[(ell, s)] = values

    varying = sum(
        len(set(values)) > 1
        for values in groups.values()
    )

    return (
        varying,
        len(groups),
    )


# ============================================================================
# SUPPORT TEST
# ============================================================================

def count_outside_support_nonzero(
    candidate: Candidate,
) -> tuple[int, int]:
    """
    Generate points just outside the proven support and check whether the
    candidate incorrectly produces nonzero values.
    """
    failures = 0
    checked = 0

    for k in [1, 3, 5, 7, 9]:
        for ell in range(
            8,
            41,
        ):
            lo = support_lower(
                k,
                ell,
            )

            hi = support_upper(
                k,
                ell,
            )

            for s in [
                lo - 2,
                lo - 1,
                hi + 1,
                hi + 2,
            ]:
                checked += 1

                value = evaluate_candidate(
                    candidate,
                    k,
                    ell,
                    s,
                )

                if value != 0:
                    failures += 1

    return failures, checked


# ============================================================================
# PARITY TEST
# ============================================================================

def parity_signature(
    value: int,
) -> int:
    if value == 0:
        return 0

    return 1 if value > 0 else -1


def print_parity_behavior(
    candidate: Candidate,
) -> None:
    print()
    print(
        f"{candidate.name}: parity behavior"
    )
    print("-" * 72)

    for ell_parity in [0, 1]:
        values = []

        for k in [1, 3, 5, 7]:
            for ell in range(
                10 + ell_parity,
                31,
                2,
            ):
                lo = support_lower(
                    k,
                    ell,
                )

                hi = support_upper(
                    k,
                    ell,
                )

                if lo <= hi:
                    s = (lo + hi) // 2

                    value = evaluate_candidate(
                        candidate,
                        k,
                        ell,
                        s,
                    )

                    values.append(
                        parity_signature(
                            value
                        )
                    )

        positive = values.count(1)
        negative = values.count(-1)
        zero = values.count(0)

        print(
            f"ell%2={ell_parity}: "
            f"positive={positive:4d} "
            f"negative={negative:4d} "
            f"zero={zero:4d}"
        )


# ============================================================================
# BOUNDARY BEHAVIOR
# ============================================================================

def print_boundary_statistics(
    candidate: Candidate,
    cases: list[Case],
) -> None:
    print()
    print(
        f"{candidate.name}: boundary statistics"
    )
    print("-" * 72)

    for location in [
        "LOWER",
        "LOWER+1",
        "INTERIOR",
        "UPPER-1",
        "UPPER",
    ]:
        subset = [
            case
            for case in cases
            if case.location == location
        ]

        if not subset:
            continue

        nonzero = sum(
            evaluate_candidate(
                candidate,
                case.k,
                case.ell,
                case.s,
            ) != 0
            for case in subset
        )

        print(
            f"{location:10s}: "
            f"{nonzero:5d}/{len(subset):5d} nonzero"
        )


# ============================================================================
# CONCRETE TABLE
# ============================================================================

def print_concrete_table(
    candidates: list[Candidate],
) -> None:
    """
    Print values for a compact collection of cases where the support
    changes visibly with k.
    """
    cases = [
        (1, 12, 5),
        (3, 12, 5),
        (5, 12, 5),
        (7, 12, 5),
        (9, 12, 5),

        (1, 12, 6),
        (3, 12, 6),
        (5, 12, 6),
        (7, 12, 6),
        (9, 12, 6),

        (1, 20, 9),
        (3, 20, 9),
        (5, 20, 9),
        (7, 20, 9),
        (9, 20, 9),
    ]

    print()
    print("=" * 72)
    print("CONCRETE CANDIDATE VALUES")
    print("=" * 72)

    for k, ell, s in cases:
        print()
        print(
            f"k={k}, ell={ell}, s={s}"
        )

        for candidate in candidates:
            value = evaluate_candidate(
                candidate,
                k,
                ell,
                s,
            )

            print(
                f"  {candidate.name:35s}"
                f" {value:12d}"
            )


# ============================================================================
# RANKING
# ============================================================================

def rank_candidates(
    candidates: list[Candidate],
    cases: list[Case],
) -> None:
    """
    Rank candidates by structural behavior.

    This is NOT a correctness score because no true reference values
    are supplied here.
    """
    rows = []

    for candidate in candidates:
        varying, groups = count_k_variation(
            candidate
        )

        nonzero = count_nonzero(
            candidate,
            cases,
        )

        support_failures, support_checked = (
            count_outside_support_nonzero(
                candidate
            )
        )

        rows.append(
            (
                candidate.name,
                varying,
                groups,
                nonzero,
                support_failures,
                support_checked,
            )
        )

    print()
    print("=" * 72)
    print("STRUCTURAL RANKING")
    print("=" * 72)

    for row in rows:
        (
            name,
            varying,
            groups,
            nonzero,
            support_failures,
            support_checked,
        ) = row

        print()
        print(name)
        print(
            f"  k-varying groups : "
            f"{varying}/{groups}"
        )
        print(
            f"  nonzero support : "
            f"{nonzero}/{len(cases)}"
        )
        print(
            f"  outside-support "
            f"nonzero         : "
            f"{support_failures}/{support_checked}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 188")
    print("Candidate k-dependent convolution reconstruction")
    print("=" * 72)

    print()
    print(
        "All data are generated locally."
    )
    print(
        "No previous result file is read."
    )

    cases = generate_cases()
    candidates = make_candidates()

    print()
    print(
        f"Generated cases: {len(cases)}"
    )

    rank_candidates(
        candidates,
        cases,
    )

    print_concrete_table(
        candidates
    )

    for candidate in candidates:
        print_parity_behavior(
            candidate
        )

        print_boundary_statistics(
            candidate,
            cases,
        )

    print()
    print("=" * 72)
    print("END EXPERIMENT 188")
    print("=" * 72)


if __name__ == "__main__":
    main()

