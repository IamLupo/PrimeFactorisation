#!/usr/bin/env python3

"""
EXPERIMENT 186
Exact-convolution sanity check and bulk reconstruction.

This experiment is deliberately designed to catch a wrong transcription
of the exact summand before attempting any sign/boundary diagnosis.

Rules:
    - No previous output is read.
    - All data are generated here.
    - Every function is defined here.
    - exact_value() and bulk_value() are independent.
    - The script detects whether the current exact summand collapses
      identically by a binomial identity.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb


# ============================================================================
# BASIC FUNCTIONS
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
# EXACT DEFINITION
# ============================================================================

def exact_summand(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    CURRENT EXACT SUMMATION FORM.

    Replace ONLY this body with the actual Experiment 183 summand.

    This function is intentionally independent of bulk_value().
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - j)
    )


def exact_terms(
    k: int,
    ell: int,
    s: int,
) -> list[tuple[int, int, int]]:
    """
    Return (j, term, cumulative_sum) for every nonzero term.
    """
    result = []
    cumulative = 0

    for j in range(ell + 1):
        term = exact_summand(
            k,
            ell,
            s,
            j,
        )

        if term != 0:
            cumulative += term

            result.append(
                (
                    j,
                    term,
                    cumulative,
                )
            )

    return result


def exact_value(
    k: int,
    ell: int,
    s: int,
) -> int:
    if not in_support(k, ell, s):
        return 0

    return sum(
        term
        for _, term, _ in exact_terms(
            k,
            ell,
            s,
        )
    )


# ============================================================================
# INDEPENDENT BULK FORM
# ============================================================================

def bulk_value(
    k: int,
    ell: int,
    s: int,
) -> int:
    """
    CURRENT BULK TEMPLATE.

    Replace with the actual Experiment 183 bulk formula.
    """
    if not in_support(k, ell, s):
        return 0

    return sgn(s) * C(ell, s)


# ============================================================================
# FRESH DATA
# ============================================================================

@dataclass(frozen=True)
class Case:
    k: int
    ell: int
    s: int


def generate_cases() -> list[Case]:
    cases = []

    for k in [1, 3, 5, 7]:
        for ell in range(
            max(k, 5),
            31,
        ):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            for s in range(lo, hi + 1):
                cases.append(
                    Case(
                        k=k,
                        ell=ell,
                        s=s,
                    )
                )

    return cases


# ============================================================================
# CANCELLATION ANALYSIS
# ============================================================================

def count_zero_exact_values(
    cases: list[Case],
) -> tuple[int, int]:
    zero = 0
    nonzero = 0

    for case in cases:
        value = exact_value(
            case.k,
            case.ell,
            case.s,
        )

        if value == 0:
            zero += 1
        else:
            nonzero += 1

    return zero, nonzero


def print_cancellation_summary(
    cases: list[Case],
) -> None:
    zero, nonzero = count_zero_exact_values(
        cases
    )

    print()
    print("=" * 72)
    print("EXACT CANCELLATION SUMMARY")
    print("=" * 72)
    print(f"cases          : {len(cases)}")
    print(f"exact == 0     : {zero}")
    print(f"exact != 0     : {nonzero}")

    if nonzero == 0:
        print()
        print(
            "WARNING: the current exact summand produces zero "
            "for every generated support point."
        )
        print(
            "Therefore this summand cannot represent a nontrivial "
            "Experiment 183 exact quantity."
        )


# ============================================================================
# SMALL HAND-CHECK TABLE
# ============================================================================

def print_small_cases() -> None:
    """
    Print small cases where the full summation can be inspected manually.
    """
    selected = [
        (3, 8, 3),
        (3, 8, 4),
        (3, 9, 4),
        (3, 9, 5),
        (5, 10, 4),
        (5, 10, 5),
    ]

    print()
    print("=" * 72)
    print("SMALL CASE TERM-BY-TERM CHECK")
    print("=" * 72)

    for k, ell, s in selected:
        print()
        print(
            f"k={k}, ell={ell}, s={s}"
        )
        print("-" * 72)

        terms = exact_terms(
            k,
            ell,
            s,
        )

        for j, term, cumulative in terms:
            print(
                f"j={j:2d} "
                f"term={term:12d} "
                f"cumulative={cumulative:12d}"
            )

        exact = exact_value(
            k,
            ell,
            s,
        )

        bulk = bulk_value(
            k,
            ell,
            s,
        )

        print()
        print(f"exact = {exact}")
        print(f"bulk  = {bulk}")


# ============================================================================
# DIRECT CONVOLUTION IDENTITY TEST
# ============================================================================

def test_current_sum_identity(
    cases: list[Case],
) -> None:
    """
    Test whether the current exact summand has the expected
    binomial-convolution collapse.

    For the currently written summand,

        C(ell,j) C(ell-j,s-j)
        = C(ell,s) C(s,j),

    so the sum becomes

        C(ell,s) sum_j (-1)^j C(s,j),

    which is zero for s > 0.

    This diagnostic confirms exactly what happened in Experiment 185.
    """
    failures = 0

    for case in cases:
        k = case.k
        ell = case.ell
        s = case.s

        direct = exact_value(
            k,
            ell,
            s,
        )

        predicted = (
            C(ell, s)
            * (
                1
                if s == 0
                else 0
            )
        )

        if direct != predicted:
            failures += 1

    print()
    print("=" * 72)
    print("CURRENT CONVOLUTION IDENTITY TEST")
    print("=" * 72)
    print(f"checked : {len(cases)}")
    print(f"failures: {failures}")

    if failures == 0:
        print()
        print(
            "The current exact summand collapses exactly to "
            "C(ell,s) * sum_j (-1)^j C(s,j)."
        )
        print(
            "Hence it is zero for every positive s."
        )


# ============================================================================
# BULK COMPARISON
# ============================================================================

def print_bulk_comparison(
    cases: list[Case],
    limit: int = 40,
) -> None:
    rows = []

    for case in cases:
        exact = exact_value(
            case.k,
            case.ell,
            case.s,
        )

        bulk = bulk_value(
            case.k,
            case.ell,
            case.s,
        )

        if exact != bulk:
            rows.append(
                (
                    case,
                    exact,
                    bulk,
                )
            )

    print()
    print("=" * 72)
    print("EXACT VS BULK")
    print("=" * 72)

    print(
        f"failures: {len(rows)}"
    )

    for case, exact, bulk in rows[:limit]:
        print(
            f"k={case.k:2d} "
            f"ell={case.ell:2d} "
            f"s={case.s:2d} "
            f"exact={exact:12d} "
            f"bulk={bulk:12d}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 186")
    print("Exact-convolution sanity check")
    print("=" * 72)

    print()
    print(
        "All test data are generated locally."
    )
    print(
        "No results.txt or previous experiment output is read."
    )

    cases = generate_cases()

    print(
        f"Generated cases: {len(cases)}"
    )

    print_cancellation_summary(
        cases
    )

    print_small_cases()

    test_current_sum_identity(
        cases
    )

    print_bulk_comparison(
        cases
    )

    print()
    print("=" * 72)
    print("END EXPERIMENT 186")
    print("=" * 72)


if __name__ == "__main__":
    main()

