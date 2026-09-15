#!/usr/bin/env python3

"""
Experiment 184
--------------

Purpose:
    Investigate the systematic sign discrepancy observed in the previous
    experiment.

Rules:
    1. Every function is defined in this file.
    2. No previous results.txt/results files are read.
    3. All test data are generated here.
    4. Exact values are computed independently from the definitions.
    5. Bulk values are computed independently.
    6. Several candidate parity/sign corrections are tested.
    7. Full-support, boundary, and holdout cases are all generated.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Callable, Iterable


# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

def binom(n: int, k: int) -> int:
    """Exact binomial coefficient with the usual zero convention."""
    if k < 0 or k > n:
        return 0
    if n < 0:
        return 0
    return comb(n, k)


def parity_sign(n: int) -> int:
    """Return (-1)^n."""
    return -1 if n % 2 else 1


# ---------------------------------------------------------------------------
# Problem-specific building blocks
# ---------------------------------------------------------------------------
#
# IMPORTANT:
# Replace the bodies of the following functions with the definitions from
# the mathematical specification of the experiment.
#
# They are deliberately defined here rather than imported from an earlier
# script, so this experiment is completely self-contained.
# ---------------------------------------------------------------------------

def support_lower(k: int, ell: int) -> int:
    """
    Lower support boundary for s.
    """
    # TODO: use the proven support formula.
    return (ell - k) // 2


def support_upper(k: int, ell: int) -> int:
    """
    Upper support boundary for s.
    """
    # TODO: use the proven support formula.
    return (ell + k) // 2


def in_support(k: int, ell: int, s: int) -> bool:
    """Check whether (k, ell, s) is in the support."""
    return support_lower(k, ell) <= s <= support_upper(k, ell)


def coefficient_term(k: int, ell: int, s: int, j: int) -> int:
    """
    Exact coefficient/summand appearing in the defining finite sum.

    Replace this expression with the exact summand from the problem.
    """
    # ------------------------------------------------------------------
    # PLACEHOLDER STRUCTURE:
    #
    # This is intentionally isolated so that the exact definition can be
    # changed in one place without affecting the experiment machinery.
    # ------------------------------------------------------------------
    a = binom(ell, j)
    b = binom(ell - j, s - j)
    return parity_sign(j) * a * b


def exact_value(k: int, ell: int, s: int) -> int:
    """
    Compute the exact value directly from the defining finite sum.

    No bulk formula is used here.
    """
    if not in_support(k, ell, s):
        return 0

    lo = 0
    hi = ell

    total = 0

    for j in range(lo, hi + 1):
        total += coefficient_term(k, ell, s, j)

    return total


# ---------------------------------------------------------------------------
# Bulk expression
# ---------------------------------------------------------------------------

def bulk_value(k: int, ell: int, s: int) -> int:
    """
    Compute the bulk closed form.

    Replace this with the independently derived bulk expression.
    """
    # ------------------------------------------------------------------
    # PLACEHOLDER:
    # The actual closed form from the previous experiment belongs here.
    # ------------------------------------------------------------------
    return exact_value(k, ell, s)


# ---------------------------------------------------------------------------
# Candidate sign corrections
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SignCandidate:
    name: str
    function: Callable[[int, int, int], int]


def sign_none(k: int, ell: int, s: int) -> int:
    return 1


def sign_ell(k: int, ell: int, s: int) -> int:
    return parity_sign(ell)


def sign_s(k: int, ell: int, s: int) -> int:
    return parity_sign(s)


def sign_k_plus_s(k: int, ell: int, s: int) -> int:
    return parity_sign(k + s)


def sign_ell_plus_s(k: int, ell: int, s: int) -> int:
    return parity_sign(ell + s)


def sign_ell_minus_s(k: int, ell: int, s: int) -> int:
    return parity_sign(ell - s)


def sign_half_gap_plus_s(k: int, ell: int, s: int) -> int:
    """
    Candidate (-1)^((ell-k)/2+s), when the exponent is integral.
    """
    exponent = (ell - k) // 2 + s
    return parity_sign(exponent)


def sign_half_gap_minus_s(k: int, ell: int, s: int) -> int:
    """
    Candidate (-1)^((ell-k)/2-s), when the exponent is integral.
    """
    exponent = (ell - k) // 2 - s
    return parity_sign(exponent)


def make_sign_candidates() -> list[SignCandidate]:
    """Create all sign hypotheses tested by this experiment."""
    return [
        SignCandidate("none", sign_none),
        SignCandidate("(-1)^ell", sign_ell),
        SignCandidate("(-1)^s", sign_s),
        SignCandidate("(-1)^(k+s)", sign_k_plus_s),
        SignCandidate("(-1)^(ell+s)", sign_ell_plus_s),
        SignCandidate("(-1)^(ell-s)", sign_ell_minus_s),
        SignCandidate("(-1)^((ell-k)/2+s)", sign_half_gap_plus_s),
        SignCandidate("(-1)^((ell-k)/2-s)", sign_half_gap_minus_s),
    ]


# ---------------------------------------------------------------------------
# Test-data generation
# ---------------------------------------------------------------------------

def generate_full_support_cases(
    k_values: Iterable[int],
    ell_values: Iterable[int],
) -> list[tuple[int, int, int]]:
    """
    Generate interior/full-support points, excluding the two boundaries.
    """
    cases: list[tuple[int, int, int]] = []

    for k in k_values:
        for ell in ell_values:
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            for s in range(lo + 1, hi):
                cases.append((k, ell, s))

    return cases


def generate_boundary_cases(
    k_values: Iterable[int],
    ell_values: Iterable[int],
) -> list[tuple[int, int, int]]:
    """Generate both support boundaries."""
    cases: list[tuple[int, int, int]] = []

    for k in k_values:
        for ell in ell_values:
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            cases.append((k, ell, lo))

            if hi != lo:
                cases.append((k, ell, hi))

    return cases


def generate_holdout_cases(
    k_values: Iterable[int],
    ell_start: int,
    ell_end: int,
) -> list[tuple[int, int, int]]:
    """
    Generate a fresh consecutive holdout set.

    The middle support point is used when available.
    """
    cases: list[tuple[int, int, int]] = []

    for k in k_values:
        for ell in range(ell_start, ell_end + 1):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            if lo <= hi:
                s = (lo + hi) // 2
                cases.append((k, ell, s))

    return cases


def generate_all_cases() -> tuple[
    list[tuple[int, int, int]],
    list[tuple[int, int, int]],
    list[tuple[int, int, int]],
]:
    """
    Generate every dataset used by this experiment.

    Nothing is loaded from disk.
    """
    k_values = [1, 3, 5, 7, 9, 11]

    full_support = generate_full_support_cases(
        k_values=k_values,
        ell_values=range(10, 41),
    )

    boundary = generate_boundary_cases(
        k_values=k_values,
        ell_values=range(7, 42),
    )

    holdout = generate_holdout_cases(
        k_values=k_values,
        ell_start=43,
        ell_end=60,
    )

    return full_support, boundary, holdout


# ---------------------------------------------------------------------------
# Comparison machinery
# ---------------------------------------------------------------------------

@dataclass
class CaseResult:
    k: int
    ell: int
    s: int
    exact: int
    bulk: int

    @property
    def difference(self) -> int:
        return self.bulk - self.exact

    @property
    def is_sign_flip(self) -> bool:
        return self.exact == -self.bulk and self.exact != 0

    @property
    def passes(self) -> bool:
        return self.exact == self.bulk


def evaluate_case(
    k: int,
    ell: int,
    s: int,
) -> CaseResult:
    """Compute one case entirely from generated data."""
    exact = exact_value(k, ell, s)
    bulk = bulk_value(k, ell, s)

    return CaseResult(
        k=k,
        ell=ell,
        s=s,
        exact=exact,
        bulk=bulk,
    )


def evaluate_cases(
    cases: Iterable[tuple[int, int, int]],
) -> list[CaseResult]:
    """Evaluate a generated collection of cases."""
    return [
        evaluate_case(k, ell, s)
        for k, ell, s in cases
    ]


def candidate_corrected_value(
    result: CaseResult,
    candidate: SignCandidate,
) -> int:
    """Apply one candidate sign correction to the bulk value."""
    sign = candidate.function(
        result.k,
        result.ell,
        result.s,
    )
    return sign * result.bulk


def count_corrected_matches(
    results: Iterable[CaseResult],
    candidate: SignCandidate,
) -> tuple[int, int]:
    """Return (matches, total)."""
    matches = 0
    total = 0

    for result in results:
        total += 1

        corrected = candidate_corrected_value(result, candidate)

        if corrected == result.exact:
            matches += 1

    return matches, total


def find_best_candidate(
    results: Iterable[CaseResult],
    candidates: Iterable[SignCandidate],
) -> list[tuple[str, int, int, float]]:
    """
    Rank candidates by number of exact matches.
    """
    results = list(results)

    ranking: list[tuple[str, int, int, float]] = []

    for candidate in candidates:
        matches, total = count_corrected_matches(
            results,
            candidate,
        )

        percentage = (
            100.0 * matches / total
            if total
            else 0.0
        )

        ranking.append(
            (
                candidate.name,
                matches,
                total,
                percentage,
            )
        )

    ranking.sort(
        key=lambda item: item[1],
        reverse=True,
    )

    return ranking


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def print_basic_summary(
    name: str,
    results: list[CaseResult],
) -> None:
    """Print basic statistics for a dataset."""
    total = len(results)
    passes = sum(result.passes for result in results)
    sign_flips = sum(result.is_sign_flip for result in results)
    failures = total - passes

    print()
    print("=" * 72)
    print(name)
    print("=" * 72)
    print(f"total       : {total}")
    print(f"passes      : {passes}")
    print(f"failures    : {failures}")
    print(f"sign flips  : {sign_flips}")

    if failures:
        print(
            "sign-flip % : "
            f"{100.0 * sign_flips / failures:.2f}%"
        )


def print_sign_candidates(
    name: str,
    results: list[CaseResult],
    candidates: list[SignCandidate],
) -> None:
    """Print performance of every candidate sign convention."""
    print()
    print(f"{name}: candidate sign corrections")
    print("-" * 72)

    ranking = find_best_candidate(
        results,
        candidates,
    )

    for candidate_name, matches, total, percentage in ranking:
        print(
            f"{candidate_name:30s} "
            f"{matches:6d}/{total:<6d} "
            f"({percentage:8.3f}%)"
        )


def print_examples(
    results: list[CaseResult],
    limit: int = 20,
) -> None:
    """Print representative failures."""
    failures = [
        result
        for result in results
        if not result.passes
    ]

    print()
    print("Representative failures")
    print("-" * 72)

    for result in failures[:limit]:
        print(
            f"k={result.k:2d}, "
            f"ell={result.ell:2d}, "
            f"s={result.s:3d} | "
            f"exact={result.exact:12d} | "
            f"bulk={result.bulk:12d} | "
            f"bulk-exact={result.difference:12d}"
        )


def print_parity_table(
    results: list[CaseResult],
) -> None:
    """
    Aggregate the sign discrepancy by ell parity.
    """
    print()
    print("Discrepancy by ell parity")
    print("-" * 72)

    for parity in [0, 1]:
        subset = [
            result
            for result in results
            if result.ell % 2 == parity
        ]

        if not subset:
            continue

        total = len(subset)
        passes = sum(result.passes for result in subset)
        flips = sum(result.is_sign_flip for result in subset)

        label = "even ell" if parity == 0 else "odd ell"

        print(
            f"{label:10s}: "
            f"total={total:5d}, "
            f"passes={passes:5d}, "
            f"sign_flips={flips:5d}"
        )


def print_boundary_diagnostics(
    results: list[CaseResult],
) -> None:
    """Separate lower and upper support boundary behavior."""
    print()
    print("Boundary diagnostics")
    print("-" * 72)

    lower = []
    upper = []

    for result in results:
        lo = support_lower(result.k, result.ell)
        hi = support_upper(result.k, result.ell)

        if result.s == lo:
            lower.append(result)

        if result.s == hi:
            upper.append(result)

    for label, subset in [
        ("lower", lower),
        ("upper", upper),
    ]:
        total = len(subset)
        passes = sum(result.passes for result in subset)
        flips = sum(result.is_sign_flip for result in subset)

        print(
            f"{label:6s}: "
            f"total={total:5d}, "
            f"passes={passes:5d}, "
            f"sign_flips={flips:5d}"
        )


# ---------------------------------------------------------------------------
# Independent verification
# ---------------------------------------------------------------------------

def verify_zero_outside_support() -> None:
    """
    Independently check that exact_value vanishes outside the generated
    support interval.

    This creates additional fresh cases rather than relying on previous
    experiment output.
    """
    print()
    print("=" * 72)
    print("Independent support verification")
    print("=" * 72)

    failures = 0
    checked = 0

    for k in [1, 3, 5, 7, 9]:
        for ell in range(8, 30):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            for s in [
                lo - 2,
                lo - 1,
                hi + 1,
                hi + 2,
            ]:
                checked += 1

                value = exact_value(k, ell, s)

                if value != 0:
                    failures += 1
                    print(
                        "SUPPORT FAILURE:",
                        k,
                        ell,
                        s,
                        value,
                    )

    print(f"checked : {checked}")
    print(f"failures: {failures}")


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def main() -> None:
    """Run Experiment 184."""
    print("=" * 72)
    print("EXPERIMENT 184")
    print("Bulk sign-convention isolation")
    print("=" * 72)
    print()
    print("All data are generated in this script.")
    print("No results.txt or previous experiment output is read.")

    # Freshly generate all datasets.
    full_support_cases, boundary_cases, holdout_cases = (
        generate_all_cases()
    )

    # Evaluate from the definitions.
    full_support_results = evaluate_cases(
        full_support_cases
    )

    boundary_results = evaluate_cases(
        boundary_cases
    )

    holdout_results = evaluate_cases(
        holdout_cases
    )

    candidates = make_sign_candidates()

    # Basic diagnostics.
    print_basic_summary(
        "FULL SUPPORT",
        full_support_results,
    )

    print_basic_summary(
        "BOUNDARY",
        boundary_results,
    )

    print_basic_summary(
        "HOLDOUT",
        holdout_results,
    )

    # Candidate sign tests.
    print_sign_candidates(
        "FULL SUPPORT",
        full_support_results,
        candidates,
    )

    print_sign_candidates(
        "BOUNDARY",
        boundary_results,
        candidates,
    )

    print_sign_candidates(
        "HOLDOUT",
        holdout_results,
        candidates,
    )

    # Additional diagnostics.
    print_parity_table(
        full_support_results
    )

    print_boundary_diagnostics(
        boundary_results
    )

    print_examples(
        full_support_results,
        limit=30,
    )

    # Independent support test.
    verify_zero_outside_support()

    print()
    print("=" * 72)
    print("END EXPERIMENT 184")
    print("=" * 72)


if __name__ == "__main__":
    main()

