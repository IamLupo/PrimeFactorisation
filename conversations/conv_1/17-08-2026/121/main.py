#!/usr/bin/env python3

"""
EXPERIMENT 185
Independent exact-vs-bulk sign diagnosis.

IMPORTANT:
- This script generates all test data itself.
- It does NOT read results.txt or any previous output.
- Every function is defined locally.
- exact_value() and bulk_value() are completely separate.
- bulk_value() must never call exact_value().
- The script tests raw equality, sign reversal, and parity-dependent signs.

Replace ONLY the mathematical definitions in the clearly marked sections
with the definitions used by Experiment 183.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb
from typing import Callable


# ============================================================================
# BASIC ARITHMETIC
# ============================================================================

def C(n: int, r: int) -> int:
    """Binomial coefficient with zero outside the natural range."""
    if n < 0:
        return 0
    if r < 0 or r > n:
        return 0
    return comb(n, r)


def sign(n: int) -> int:
    """Return (-1)^n."""
    return -1 if n % 2 else 1


def parity(n: int) -> int:
    """Return n mod 2."""
    return n & 1


# ============================================================================
# SUPPORT
# ============================================================================

def support_lower(k: int, ell: int) -> int:
    """
    Proven lower support boundary.

    Replace only if Experiment 183 uses a different proven formula.
    """
    return (ell - k) // 2


def support_upper(k: int, ell: int) -> int:
    """
    Proven upper support boundary.

    Replace only if Experiment 183 uses a different proven formula.
    """
    return (ell + k) // 2


def in_support(k: int, ell: int, s: int) -> bool:
    """Return whether s lies in the proven support interval."""
    return (
        support_lower(k, ell)
        <= s
        <= support_upper(k, ell)
    )


# ============================================================================
# EXACT DEFINITION
# ============================================================================
#
# THIS SECTION MUST CONTAIN THE ACTUAL EXACT FORMULA FROM EXPERIMENT 183.
#
# It is intentionally independent of the bulk formula below.
# ============================================================================

def exact_summand(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    One term of the exact defining sum.

    IMPORTANT:
    Replace this body with the actual exact summand from Experiment 183.
    Do NOT call bulk_value() here.
    """

    # ------------------------------------------------------------------------
    # CURRENT TEMPLATE.
    #
    # Replace with the real summand.
    # ------------------------------------------------------------------------
    return (
        sign(j)
        * C(ell, j)
        * C(ell - j, s - j)
    )


def exact_value(
    k: int,
    ell: int,
    s: int,
) -> int:
    """
    Exact value from the defining sum.

    This function knows nothing about the bulk closed form.
    """
    if not in_support(k, ell, s):
        return 0

    total = 0

    for j in range(ell + 1):
        total += exact_summand(
            k,
            ell,
            s,
            j,
        )

    return total


# ============================================================================
# BULK CLOSED FORM
# ============================================================================
#
# THIS SECTION MUST CONTAIN THE ACTUAL BULK FORMULA FROM EXPERIMENT 183.
#
# CRITICAL:
# bulk_value() MUST NOT call exact_value().
# ============================================================================

def bulk_core(
    k: int,
    ell: int,
    s: int,
) -> int:
    """
    Core bulk closed form.

    Replace this with the actual independently derived bulk expression
    from Experiment 183.

    DO NOT compute this by evaluating the exact finite sum.
    """

    # ------------------------------------------------------------------------
    # CURRENT TEMPLATE ONLY.
    #
    # Replace this entire expression with the Experiment 183 bulk formula.
    # ------------------------------------------------------------------------

    return (
        C(ell, s)
        * sign(s)
    )


def bulk_value(
    k: int,
    ell: int,
    s: int,
) -> int:
    """
    Bulk value.

    This is deliberately a separate computation.
    """
    if not in_support(k, ell, s):
        return 0

    return bulk_core(
        k,
        ell,
        s,
    )


# ============================================================================
# INDEPENDENCE GUARD
# ============================================================================

def assert_independent_implementations() -> None:
    """
    Static sanity check that bulk_value's source does not contain
    'exact_value'.

    This is deliberately simple: the important rule is that the bulk
    implementation cannot silently delegate to the exact implementation.
    """
    source = bulk_value.__code__

    names = set(source.co_names)

    if "exact_value" in names:
        raise RuntimeError(
            "INVALID EXPERIMENT: bulk_value() calls exact_value()."
        )

    if "exact_summand" in names:
        raise RuntimeError(
            "INVALID EXPERIMENT: bulk_value() depends on exact_summand()."
        )


# ============================================================================
# GENERATED TEST DATA
# ============================================================================

@dataclass(frozen=True)
class Case:
    k: int
    ell: int
    s: int
    label: str


def generate_full_support_cases() -> list[Case]:
    """
    Fresh interior/full-support cases.
    """
    cases: list[Case] = []

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(8, 51):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            for s in range(lo + 1, hi):
                cases.append(
                    Case(
                        k=k,
                        ell=ell,
                        s=s,
                        label="FULL",
                    )
                )

    return cases


def generate_boundary_cases() -> list[Case]:
    """
    Fresh lower and upper boundary cases.
    """
    cases: list[Case] = []

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(7, 51):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            cases.append(
                Case(
                    k=k,
                    ell=ell,
                    s=lo,
                    label="LOWER",
                )
            )

            if hi != lo:
                cases.append(
                    Case(
                        k=k,
                        ell=ell,
                        s=hi,
                        label="UPPER",
                    )
                )

    return cases


def generate_near_boundary_cases() -> list[Case]:
    """
    Fresh cases one and two positions inside each boundary.
    """
    cases: list[Case] = []

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(8, 51):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            for d in [0, 1, 2]:
                left = lo + d
                right = hi - d

                if lo <= left <= hi:
                    cases.append(
                        Case(
                            k=k,
                            ell=ell,
                            s=left,
                            label=f"LOWER+d{d}",
                        )
                    )

                if lo <= right <= hi:
                    cases.append(
                        Case(
                            k=k,
                            ell=ell,
                            s=right,
                            label=f"UPPER-d{d}",
                        )
                    )

    return cases


def generate_holdout_cases() -> list[Case]:
    """
    Completely fresh holdout region.

    The holdout ell values do not overlap the earlier main range.
    """
    cases: list[Case] = []

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(53, 73):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            # Several independently selected interior points.
            candidates = [
                lo,
                lo + 1,
                (lo + hi) // 2,
                hi - 1,
                hi,
            ]

            for s in candidates:
                if lo <= s <= hi:
                    cases.append(
                        Case(
                            k=k,
                            ell=ell,
                            s=s,
                            label="HOLDOUT",
                        )
                    )

    return cases


def generate_all_cases() -> list[Case]:
    """Generate the complete fresh dataset."""
    cases = []

    cases.extend(generate_full_support_cases())
    cases.extend(generate_boundary_cases())
    cases.extend(generate_near_boundary_cases())
    cases.extend(generate_holdout_cases())

    return cases


# ============================================================================
# RESULT RECORD
# ============================================================================

@dataclass(frozen=True)
class Result:
    case: Case
    exact: int
    bulk: int

    @property
    def difference(self) -> int:
        return self.bulk - self.exact

    @property
    def equal(self) -> bool:
        return self.exact == self.bulk

    @property
    def sign_flip(self) -> bool:
        return (
            self.exact != 0
            and self.bulk == -self.exact
        )

    @property
    def ratio(self) -> Fraction | None:
        if self.bulk == 0:
            return None

        return Fraction(
            self.exact,
            self.bulk,
        )


def evaluate_case(case: Case) -> Result:
    """Evaluate exact and bulk independently."""
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

    return Result(
        case=case,
        exact=exact,
        bulk=bulk,
    )


def evaluate_cases(
    cases: list[Case],
) -> list[Result]:
    """Evaluate every generated case."""
    return [
        evaluate_case(case)
        for case in cases
    ]


# ============================================================================
# SIGN HYPOTHESES
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    name: str
    function: Callable[[int, int, int], int]


def candidate_none(k: int, ell: int, s: int) -> int:
    return 1


def candidate_ell(k: int, ell: int, s: int) -> int:
    return sign(ell)


def candidate_s(k: int, ell: int, s: int) -> int:
    return sign(s)


def candidate_k(k: int, ell: int, s: int) -> int:
    return sign(k)


def candidate_k_plus_s(k: int, ell: int, s: int) -> int:
    return sign(k + s)


def candidate_ell_plus_s(k: int, ell: int, s: int) -> int:
    return sign(ell + s)


def candidate_ell_minus_s(k: int, ell: int, s: int) -> int:
    return sign(ell - s)


def candidate_k_plus_ell(k: int, ell: int, s: int) -> int:
    return sign(k + ell)


def candidate_k_plus_ell_plus_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return sign(k + ell + s)


def candidate_gap_plus_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return sign((ell - k) // 2 + s)


def candidate_gap_minus_s(
    k: int,
    ell: int,
    s: int,
) -> int:
    return sign((ell - k) // 2 - s)


def make_candidates() -> list[Candidate]:
    """Construct every tested parity hypothesis."""
    return [
        Candidate("none", candidate_none),
        Candidate("(-1)^k", candidate_k),
        Candidate("(-1)^ell", candidate_ell),
        Candidate("(-1)^s", candidate_s),
        Candidate("(-1)^(k+s)", candidate_k_plus_s),
        Candidate("(-1)^(k+ell)", candidate_k_plus_ell),
        Candidate("(-1)^(ell+s)", candidate_ell_plus_s),
        Candidate("(-1)^(ell-s)", candidate_ell_minus_s),
        Candidate(
            "(-1)^(k+ell+s)",
            candidate_k_plus_ell_plus_s,
        ),
        Candidate(
            "(-1)^((ell-k)/2+s)",
            candidate_gap_plus_s,
        ),
        Candidate(
            "(-1)^((ell-k)/2-s)",
            candidate_gap_minus_s,
        ),
    ]


def candidate_match_count(
    results: list[Result],
    candidate: Candidate,
) -> tuple[int, int]:
    """Count exact matches after applying one candidate sign."""
    matches = 0

    for result in results:
        c = result.case

        corrected = (
            candidate.function(
                c.k,
                c.ell,
                c.s,
            )
            * result.bulk
        )

        if corrected == result.exact:
            matches += 1

    return matches, len(results)


# ============================================================================
# TERM-BY-TERM DIAGNOSTICS
# ============================================================================

def exact_term_vector(
    k: int,
    ell: int,
    s: int,
) -> list[int]:
    """
    Return every exact summand separately.

    This is useful if the global discrepancy is caused by an alternating
    sign introduced during summation.
    """
    return [
        exact_summand(k, ell, s, j)
        for j in range(ell + 1)
    ]


def alternating_term_vector(
    k: int,
    ell: int,
    s: int,
) -> list[int]:
    """
    Return the exact summands with their signs flipped.

    This does NOT alter exact_value(); it is only a diagnostic.
    """
    return [
        -term
        for term in exact_term_vector(k, ell, s)
    ]


def print_term_diagnostic(
    case: Case,
) -> None:
    """Print term-level information for one selected failure."""
    terms = exact_term_vector(
        case.k,
        case.ell,
        case.s,
    )

    print()
    print("=" * 72)
    print("TERM DIAGNOSTIC")
    print("=" * 72)

    print(
        f"k={case.k}, ell={case.ell}, "
        f"s={case.s}, label={case.label}"
    )

    print()
    print("j       exact summand")
    print("-" * 40)

    for j, value in enumerate(terms):
        if value != 0:
            print(
                f"{j:3d} {value:20d}"
            )

    print()
    print(
        "sum(exact terms) =",
        sum(terms),
    )

    print(
        "sum(negated terms) =",
        sum(alternating_term_vector(
            case.k,
            case.ell,
            case.s,
        )),
    )

    print(
        "independent bulk =",
        bulk_value(
            case.k,
            case.ell,
            case.s,
        ),
    )


# ============================================================================
# REPORTING
# ============================================================================

def print_summary(
    results: list[Result],
) -> None:
    """Print global comparison statistics."""
    total = len(results)

    equal = sum(
        result.equal
        for result in results
    )

    failures = total - equal

    flips = sum(
        result.sign_flip
        for result in results
    )

    print()
    print("=" * 72)
    print("GLOBAL SUMMARY")
    print("=" * 72)

    print(f"total       : {total}")
    print(f"equal       : {equal}")
    print(f"failures    : {failures}")
    print(f"sign flips  : {flips}")

    if failures:
        print(
            f"sign-flip % : "
            f"{100.0 * flips / failures:.3f}%"
        )


def print_by_label(
    results: list[Result],
) -> None:
    """Print statistics for each generated dataset."""
    labels = sorted(
        set(result.case.label for result in results)
    )

    print()
    print("=" * 72)
    print("BY DATASET")
    print("=" * 72)

    for label in labels:
        subset = [
            result
            for result in results
            if result.case.label == label
        ]

        total = len(subset)
        equal = sum(
            result.equal
            for result in subset
        )
        flips = sum(
            result.sign_flip
            for result in subset
        )

        print(
            f"{label:15s} "
            f"total={total:5d} "
            f"equal={equal:5d} "
            f"fail={total-equal:5d} "
            f"flip={flips:5d}"
        )


def print_by_ell_parity(
    results: list[Result],
) -> None:
    """Test the suspected even/odd ell pattern."""
    print()
    print("=" * 72)
    print("BY ELL PARITY")
    print("=" * 72)

    for p in [0, 1]:
        subset = [
            result
            for result in results
            if parity(result.case.ell) == p
        ]

        total = len(subset)
        equal = sum(
            result.equal
            for result in subset
        )
        flips = sum(
            result.sign_flip
            for result in subset
        )

        print(
            f"ell % 2 = {p}: "
            f"total={total:5d} "
            f"equal={equal:5d} "
            f"fail={total-equal:5d} "
            f"flip={flips:5d}"
        )


def print_by_k_parity(
    results: list[Result],
) -> None:
    """Test whether k parity participates in the discrepancy."""
    print()
    print("=" * 72)
    print("BY K PARITY")
    print("=" * 72)

    for p in [0, 1]:
        subset = [
            result
            for result in results
            if parity(result.case.k) == p
        ]

        if not subset:
            continue

        total = len(subset)
        equal = sum(
            result.equal
            for result in subset
        )
        flips = sum(
            result.sign_flip
            for result in subset
        )

        print(
            f"k % 2 = {p}: "
            f"total={total:5d} "
            f"equal={equal:5d} "
            f"fail={total-equal:5d} "
            f"flip={flips:5d}"
        )


def print_candidate_table(
    results: list[Result],
) -> None:
    """Rank all candidate sign corrections."""
    candidates = make_candidates()

    rows = []

    for candidate in candidates:
        matches, total = candidate_match_count(
            results,
            candidate,
        )

        rows.append(
            (
                candidate.name,
                matches,
                total,
                100.0 * matches / total,
            )
        )

    rows.sort(
        key=lambda row: row[1],
        reverse=True,
    )

    print()
    print("=" * 72)
    print("SIGN-CANDIDATE TEST")
    print("=" * 72)

    for name, matches, total, percentage in rows:
        print(
            f"{name:30s} "
            f"{matches:6d}/{total:<6d} "
            f"({percentage:8.3f}%)"
        )


def print_failures(
    results: list[Result],
    limit: int = 30,
) -> list[Result]:
    """Print the first failures and return them."""
    failures = [
        result
        for result in results
        if not result.equal
    ]

    print()
    print("=" * 72)
    print("FIRST FAILURES")
    print("=" * 72)

    if not failures:
        print("NONE")
        return failures

    for result in failures[:limit]:
        c = result.case

        print(
            f"k={c.k:2d} "
            f"ell={c.ell:3d} "
            f"s={c.s:4d} "
            f"{c.label:12s} "
            f"exact={result.exact:15d} "
            f"bulk={result.bulk:15d} "
            f"diff={result.difference:15d} "
            f"ratio={result.ratio}"
        )

    return failures


# ============================================================================
# RATIO ANALYSIS
# ============================================================================

def print_ratio_classes(
    results: list[Result],
) -> None:
    """
    Count exact/bulk ratios.

    If the previous sign hypothesis is correct, -1 should dominate failures.
    """
    counts: dict[Fraction | None, int] = {}

    for result in results:
        if result.equal:
            continue

        ratio = result.ratio

        counts[ratio] = counts.get(ratio, 0) + 1

    print()
    print("=" * 72)
    print("EXACT / BULK RATIO CLASSES FOR FAILURES")
    print("=" * 72)

    if not counts:
        print("No failures.")
        return

    for ratio, count in sorted(
        counts.items(),
        key=lambda item: item[1],
        reverse=True,
    ):
        print(
            f"ratio={str(ratio):>12s} "
            f"count={count:6d}"
        )


# ============================================================================
# SUPPORT CHECK
# ============================================================================

def verify_support_independently() -> None:
    """
    Generate points outside the support and verify exact_value is zero.
    """
    failures = 0
    checked = 0

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(8, 41):
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            for s in [
                lo - 3,
                lo - 2,
                lo - 1,
                hi + 1,
                hi + 2,
                hi + 3,
            ]:
                checked += 1

                value = exact_value(
                    k,
                    ell,
                    s,
                )

                if value != 0:
                    failures += 1

                    print(
                        "SUPPORT FAILURE:",
                        f"k={k}",
                        f"ell={ell}",
                        f"s={s}",
                        f"value={value}",
                    )

    print()
    print("=" * 72)
    print("SUPPORT CHECK")
    print("=" * 72)
    print(f"checked : {checked}")
    print(f"failures: {failures}")


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    """Run Experiment 185."""
    print("=" * 72)
    print("EXPERIMENT 185")
    print("Independent exact-vs-bulk sign diagnosis")
    print("=" * 72)

    print()
    print("Generating all data locally.")
    print("No previous result file is read.")

    # First make sure the implementations are genuinely independent.
    assert_independent_implementations()

    # Generate fresh data.
    cases = generate_all_cases()

    print()
    print(f"Generated cases: {len(cases)}")

    # Evaluate exact and bulk independently.
    results = evaluate_cases(cases)

    # Core diagnostics.
    print_summary(results)
    print_by_label(results)
    print_by_ell_parity(results)
    print_by_k_parity(results)

    # Sign diagnosis.
    print_candidate_table(results)
    print_ratio_classes(results)

    # Show concrete failures.
    failures = print_failures(
        results,
        limit=40,
    )

    # If there are failures, inspect the first one term by term.
    if failures:
        print_term_diagnostic(
            failures[0].case
        )

    # Independent support verification.
    verify_support_independently()

    print()
    print("=" * 72)
    print("END EXPERIMENT 185")
    print("=" * 72)


if __name__ == "__main__":
    main()

