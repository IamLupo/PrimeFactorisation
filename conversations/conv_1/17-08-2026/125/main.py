#!/usr/bin/env python3

"""
========================================================================
EXPERIMENT 189
Generating-function reconstruction from support geometry
========================================================================

Purpose
-------
Experiment 188 showed that arbitrary k-dependent shifts do not reproduce
the observed structure well.

The proven support is

    (ell-k)/2 <= s <= (ell+k)/2

for odd k.

That geometry strongly suggests a coefficient representation involving

    x^((ell-k)/2) (1+x)^k

and a polynomial in x^2.

This experiment systematically tests generating-function families of the
form

    x^shift * (1+x)^k * (1-x^2)^b

and sign/reflection variants.

Rules
-----
1. All data are generated locally.
2. No results.txt is read.
3. No previous experiment output is read.
4. Every function is defined in this script.
5. No candidate is declared correct merely because it has the right
   support.
6. The experiment measures:
      - exact support,
      - zero region,
      - k-dependence,
      - parity,
      - boundary activity,
      - symmetry.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from typing import Callable


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
# POLYNOMIAL COEFFICIENT HELPERS
# ============================================================================

def coefficient_one_plus_x_power(
    k: int,
    r: int,
) -> int:
    return C(k, r)


def coefficient_one_minus_x2_power(
    b: int,
    q: int,
) -> int:
    if q < 0:
        return 0

    if q % 1 != 0:
        return 0

    return (
        sgn(q)
        * C(b, q)
    )


def coefficient_one_plus_x2_power(
    b: int,
    q: int,
) -> int:
    if q < 0:
        return 0

    return C(b, q)


# ============================================================================
# GENERIC CONVOLUTION
# ============================================================================

def coeff_shifted_product(
    shift: int,
    k: int,
    b: int,
    s: int,
    sign_x2: int,
    sign_x: int,
) -> int:
    """
    Coefficient of

        x^shift * (1 + sign_x*x)^k
                  * (1 + sign_x2*x^2)^b

    at x^s.
    """
    target = s - shift

    if target < 0:
        return 0

    total = 0

    for r in range(k + 1):
        remaining = target - r

        if remaining < 0:
            continue

        if remaining % 2 != 0:
            continue

        q = remaining // 2

        if q > b:
            continue

        first = (
            sign_x ** r
            * C(k, r)
        )

        second = (
            sign_x2 ** q
            * C(b, q)
        )

        total += first * second

    return total


# ============================================================================
# CANDIDATE DESCRIPTIONS
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    name: str
    shift_fn: Callable[[int, int], int]
    b_fn: Callable[[int, int], int]
    sign_x: int
    sign_x2: int


def shift_left(
    k: int,
    ell: int,
) -> int:
    return (ell - k) // 2


def shift_right(
    k: int,
    ell: int,
) -> int:
    return (ell + k) // 2


def shift_zero(
    k: int,
    ell: int,
) -> int:
    return 0


def b_half_difference(
    k: int,
    ell: int,
) -> int:
    return (ell - k) // 2


def b_half_sum(
    k: int,
    ell: int,
) -> int:
    return (ell + k) // 2


def b_ell(
    k: int,
    ell: int,
) -> int:
    return ell


def b_difference_plus_one(
    k: int,
    ell: int,
) -> int:
    return (ell - k) // 2 + 1


def make_candidates() -> list[Candidate]:
    candidates = []

    shift_fns = [
        ("left", shift_left),
        ("right", shift_right),
        ("zero", shift_zero),
    ]

    b_fns = [
        ("half-difference", b_half_difference),
        ("half-sum", b_half_sum),
        ("ell", b_ell),
        ("difference+1", b_difference_plus_one),
    ]

    for shift_name, shift_fn in shift_fns:
        for b_name, b_fn in b_fns:
            for sign_x in [1, -1]:
                for sign_x2 in [1, -1]:
                    candidates.append(
                        Candidate(
                            name=(
                                f"shift={shift_name}, "
                                f"b={b_name}, "
                                f"x={sign_x:+d}, "
                                f"x2={sign_x2:+d}"
                            ),
                            shift_fn=shift_fn,
                            b_fn=b_fn,
                            sign_x=sign_x,
                            sign_x2=sign_x2,
                        )
                    )

    return candidates


# ============================================================================
# CANDIDATE VALUE
# ============================================================================

def candidate_value(
    candidate: Candidate,
    k: int,
    ell: int,
    s: int,
) -> int:
    shift = candidate.shift_fn(
        k,
        ell,
    )

    b = candidate.b_fn(
        k,
        ell,
    )

    if b < 0:
        return 0

    return coeff_shifted_product(
        shift=shift,
        k=k,
        b=b,
        s=s,
        sign_x=candidate.sign_x,
        sign_x2=candidate.sign_x2,
    )


# ============================================================================
# GENERATED DATA
# ============================================================================

@dataclass(frozen=True)
class Case:
    k: int
    ell: int
    s: int
    location: str


def generate_cases() -> list[Case]:
    cases = []

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(
            k,
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

            for s in range(
                lo - 2,
                hi + 3,
            ):
                if s < lo:
                    location = "BELOW"
                elif s > hi:
                    location = "ABOVE"
                elif s == lo:
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
# SUPPORT SCORE
# ============================================================================

def support_score(
    candidate: Candidate,
    cases: list[Case],
) -> tuple[int, int]:
    """
    Returns:
        inside failures,
        outside failures
    """
    inside_failures = 0
    outside_failures = 0

    for case in cases:
        value = candidate_value(
            candidate,
            case.k,
            case.ell,
            case.s,
        )

        inside = in_support(
            case.k,
            case.ell,
            case.s,
        )

        if inside and value == 0:
            inside_failures += 1

        if not inside and value != 0:
            outside_failures += 1

    return (
        inside_failures,
        outside_failures,
    )


# ============================================================================
# K-DEPENDENCE
# ============================================================================

def k_dependence_score(
    candidate: Candidate,
) -> tuple[int, int]:
    groups = {}

    for ell in range(
        8,
        41,
    ):
        for s in range(
            0,
            ell + 10,
        ):
            values = []

            for k in [1, 3, 5, 7, 9, 11]:
                if in_support(
                    k,
                    ell,
                    s,
                ):
                    values.append(
                        candidate_value(
                            candidate,
                            k,
                            ell,
                            s,
                        )
                    )

            if len(values) >= 2:
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
# PARITY STRUCTURE
# ============================================================================

def parity_score(
    candidate: Candidate,
) -> tuple[int, int, int]:
    positive = 0
    negative = 0
    zero = 0

    for k in [1, 3, 5, 7, 9]:
        for ell in range(
            8,
            31,
        ):
            lo = support_lower(
                k,
                ell,
            )

            hi = support_upper(
                k,
                ell,
            )

            for s in range(
                lo,
                hi + 1,
            ):
                value = candidate_value(
                    candidate,
                    k,
                    ell,
                    s,
                )

                if value > 0:
                    positive += 1
                elif value < 0:
                    negative += 1
                else:
                    zero += 1

    return (
        positive,
        negative,
        zero,
    )


# ============================================================================
# BOUNDARY ACTIVITY
# ============================================================================

def boundary_score(
    candidate: Candidate,
) -> dict[str, tuple[int, int]]:
    result = {}

    for location in [
        "LOWER",
        "LOWER+1",
        "INTERIOR",
        "UPPER-1",
        "UPPER",
    ]:
        total = 0
        nonzero = 0

        for k in [1, 3, 5, 7, 9]:
            for ell in range(
                max(k, 8),
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

                if location == "LOWER":
                    s = lo
                elif location == "UPPER":
                    s = hi
                elif location == "LOWER+1":
                    s = lo + 1
                elif location == "UPPER-1":
                    s = hi - 1
                else:
                    if hi - lo < 4:
                        continue

                    s = (lo + hi) // 2

                total += 1

                if candidate_value(
                    candidate,
                    k,
                    ell,
                    s,
                ) != 0:
                    nonzero += 1

        result[location] = (
            nonzero,
            total,
        )

    return result


# ============================================================================
# SUPPORT WIDTH TEST
# ============================================================================

def width_identity_score(
    candidate: Candidate,
) -> tuple[int, int]:
    """
    Test whether the candidate's first/last nonzero coefficients agree
    with the expected support endpoints.

    This is stronger than simply checking a few outside points.
    """
    failures = 0
    checked = 0

    for k in [1, 3, 5, 7, 9, 11]:
        for ell in range(
            k,
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

            nonzero_indices = []

            for s in range(
                0,
                ell + k + 3,
            ):
                if candidate_value(
                    candidate,
                    k,
                    ell,
                    s,
                ) != 0:
                    nonzero_indices.append(s)

            checked += 1

            if not nonzero_indices:
                failures += 1
                continue

            first = min(
                nonzero_indices
            )
            last = max(
                nonzero_indices
            )

            if first != lo or last != hi:
                failures += 1

    return (
        failures,
        checked,
    )


# ============================================================================
# SYMMETRY TEST
# ============================================================================

def symmetry_score(
    candidate: Candidate,
) -> tuple[int, int]:
    """
    Test the natural reflection

        s -> ell-s

    and record how often the magnitudes agree.

    No sign convention is assumed.
    """
    matches = 0
    checked = 0

    for k in [1, 3, 5, 7, 9]:
        for ell in range(
            8,
            31,
        ):
            lo = support_lower(
                k,
                ell,
            )

            hi = support_upper(
                k,
                ell,
            )

            for s in range(
                lo,
                hi + 1,
            ):
                reflected = ell - s

                if not in_support(
                    k,
                    ell,
                    reflected,
                ):
                    continue

                a = candidate_value(
                    candidate,
                    k,
                    ell,
                    s,
                )

                b = candidate_value(
                    candidate,
                    k,
                    ell,
                    reflected,
                )

                checked += 1

                if abs(a) == abs(b):
                    matches += 1

    return (
        matches,
        checked,
    )


# ============================================================================
# RANKING
# ============================================================================

@dataclass(frozen=True)
class Score:
    candidate: Candidate
    inside_support_failures: int
    outside_support_failures: int
    width_failures: int
    width_checked: int
    varying_groups: int
    total_groups: int
    positive: int
    negative: int
    zero: int
    symmetry_matches: int
    symmetry_checked: int


def score_candidate(
    candidate: Candidate,
    cases: list[Case],
) -> Score:
    inside_failures, outside_failures = (
        support_score(
            candidate,
            cases,
        )
    )

    width_failures, width_checked = (
        width_identity_score(
            candidate
        )
    )

    varying_groups, total_groups = (
        k_dependence_score(
            candidate
        )
    )

    positive, negative, zero = (
        parity_score(
            candidate
        )
    )

    symmetry_matches, symmetry_checked = (
        symmetry_score(
            candidate
        )
    )

    return Score(
        candidate=candidate,
        inside_support_failures=inside_failures,
        outside_support_failures=outside_failures,
        width_failures=width_failures,
        width_checked=width_checked,
        varying_groups=varying_groups,
        total_groups=total_groups,
        positive=positive,
        negative=negative,
        zero=zero,
        symmetry_matches=symmetry_matches,
        symmetry_checked=symmetry_checked,
    )


def print_ranking(
    scores: list[Score],
) -> None:
    print()
    print("=" * 72)
    print("GENERATING-FUNCTION STRUCTURAL RANKING")
    print("=" * 72)

    ordered = sorted(
        scores,
        key=lambda x: (
            x.width_failures,
            x.outside_support_failures,
            x.inside_support_failures,
        ),
    )

    for score in ordered:
        print()
        print(score.candidate.name)
        print(
            f"  inside support zeros : "
            f"{score.inside_support_failures}"
        )
        print(
            f"  outside support     : "
            f"{score.outside_support_failures}"
        )
        print(
            f"  endpoint failures   : "
            f"{score.width_failures}/"
            f"{score.width_checked}"
        )
        print(
            f"  k-varying groups    : "
            f"{score.varying_groups}/"
            f"{score.total_groups}"
        )
        print(
            f"  signs               : "
            f"+{score.positive} "
            f"-{score.negative} "
            f"0={score.zero}"
        )
        print(
            f"  abs symmetry        : "
            f"{score.symmetry_matches}/"
            f"{score.symmetry_checked}"
        )


# ============================================================================
# REPRESENTATIVE COEFFICIENT TABLE
# ============================================================================

def print_representatives(
    candidates: list[Candidate],
) -> None:
    representatives = [
        (3, 8),
        (3, 10),
        (3, 12),
        (5, 10),
        (5, 12),
        (7, 14),
        (9, 18),
    ]

    print()
    print("=" * 72)
    print("REPRESENTATIVE COEFFICIENTS")
    print("=" * 72)

    for k, ell in representatives:
        lo = support_lower(
            k,
            ell,
        )

        hi = support_upper(
            k,
            ell,
        )

        print()
        print(
            f"k={k}, ell={ell}, "
            f"support=[{lo},{hi}]"
        )

        for s in range(
            lo,
            hi + 1,
        ):
            values = []

            for candidate in candidates:
                value = candidate_value(
                    candidate,
                    k,
                    ell,
                    s,
                )

                if value != 0:
                    values.append(
                        (
                            candidate.name,
                            value,
                        )
                    )

            print(
                f"  s={s:2d}:"
            )

            for name, value in values:
                print(
                    f"    {name:42s}"
                    f"{value:12d}"
                )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 189")
    print("Generating-function reconstruction from support geometry")
    print("=" * 72)

    print()
    print(
        "All data are generated locally."
    )
    print(
        "No results.txt or previous experiment output is read."
    )

    cases = generate_cases()
    candidates = make_candidates()

    print()
    print(
        f"Generated cases: {len(cases)}"
    )

    scores = [
        score_candidate(
            candidate,
            cases,
        )
        for candidate in candidates
    ]

    print_ranking(
        scores
    )

    print_representatives(
        candidates
    )

    print()
    print("=" * 72)
    print("END EXPERIMENT 189")
    print("=" * 72)


if __name__ == "__main__":
    main()

