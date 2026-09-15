#!/usr/bin/env python3

"""
EXPERIMENT 190
Degree-k generating-polynomial reconstruction.

Key observation
---------------
The proven support is

    L = (ell-k)/2
    U = (ell+k)/2

so the support width is exactly

    U-L = k.

Therefore, after factoring x^L, the remaining polynomial must have
degree exactly k.

For odd k, write

    m = (k-1)/2.

This experiment tests the natural degree-k families

    (1+x)^k
    (1-x)^k
    (1+x)(1-x^2)^m
    (1-x)(1-x^2)^m

together with reflections/sign normalizations.

Rules
-----
- No results.txt is read.
- No previous experiment output is read.
- All data are generated locally.
- Every function is defined in this file.
- No candidate uses an exact-value oracle.
- Output is compact: no enormous coefficient dump.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
from math import comb


# ============================================================================
# BASIC FUNCTIONS
# ============================================================================

def C(n: int, r: int) -> int:
    """Binomial coefficient with zero convention."""
    if n < 0:
        return 0
    if r < 0 or r > n:
        return 0
    return comb(n, r)


def sgn(n: int) -> int:
    """Return (-1)^n."""
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
# POLYNOMIAL REPRESENTATION
# ============================================================================

def polynomial_multiply(
    a: list[int],
    b: list[int],
) -> list[int]:
    """Multiply two integer coefficient polynomials."""
    result = [0] * (len(a) + len(b) - 1)

    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] += ai * bj

    return result


def polynomial_power(
    base: list[int],
    n: int,
) -> list[int]:
    """Raise a polynomial to a nonnegative integer power."""
    result = [1]

    for _ in range(n):
        result = polynomial_multiply(
            result,
            base,
        )

    return result


def polynomial_reflect(
    p: list[int],
) -> list[int]:
    """Reverse coefficient order."""
    return list(reversed(p))


def polynomial_negate(
    p: list[int],
) -> list[int]:
    """Multiply every coefficient by -1."""
    return [-x for x in p]


# ============================================================================
# NATURAL DEGREE-k FAMILIES
# ============================================================================

def poly_binomial_plus(
    k: int,
) -> list[int]:
    """
    (1+x)^k
    """
    return [
        C(k, r)
        for r in range(k + 1)
    ]


def poly_binomial_minus(
    k: int,
) -> list[int]:
    """
    (1-x)^k
    """
    return [
        sgn(r) * C(k, r)
        for r in range(k + 1)
    ]


def poly_odd_plus(
    k: int,
) -> list[int]:
    """
    (1+x)(1-x^2)^m,
    m=(k-1)/2.
    """
    if k % 2 != 1:
        return []

    m = (k - 1) // 2

    even_part = [
        sgn(q) * C(m, q)
        for q in range(m + 1)
    ]

    return polynomial_multiply(
        [1, 1],
        even_part,
    )


def poly_odd_minus(
    k: int,
) -> list[int]:
    """
    (1-x)(1-x^2)^m,
    m=(k-1)/2.
    """
    if k % 2 != 1:
        return []

    m = (k - 1) // 2

    even_part = [
        sgn(q) * C(m, q)
        for q in range(m + 1)
    ]

    return polynomial_multiply(
        [1, -1],
        even_part,
    )


# ============================================================================
# CANDIDATES
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    name: str
    polynomial_fn: Callable[[int], list[int]]
    reflected: bool
    global_sign: int


def make_candidates() -> list[Candidate]:
    base = [
        ("(1+x)^k", poly_binomial_plus),
        ("(1-x)^k", poly_binomial_minus),
        (
            "(1+x)(1-x^2)^((k-1)/2)",
            poly_odd_plus,
        ),
        (
            "(1-x)(1-x^2)^((k-1)/2)",
            poly_odd_minus,
        ),
    ]

    candidates: list[Candidate] = []

    for name, fn in base:
        for reflected in [False, True]:
            for global_sign in [1, -1]:
                suffix = []

                if reflected:
                    suffix.append("reflected")

                if global_sign == -1:
                    suffix.append("global -")

                label = name

                if suffix:
                    label += " [" + ", ".join(suffix) + "]"

                candidates.append(
                    Candidate(
                        name=label,
                        polynomial_fn=fn,
                        reflected=reflected,
                        global_sign=global_sign,
                    )
                )

    return candidates


# ============================================================================
# CANDIDATE COEFFICIENT
# ============================================================================

def candidate_coefficients(
    candidate: Candidate,
    k: int,
) -> list[int]:
    """Return the degree-k polynomial coefficient vector."""
    coeffs = candidate.polynomial_fn(k)

    if candidate.reflected:
        coeffs = polynomial_reflect(
            coeffs
        )

    if candidate.global_sign == -1:
        coeffs = polynomial_negate(
            coeffs
        )

    return coeffs


def candidate_value(
    candidate: Candidate,
    k: int,
    ell: int,
    s: int,
) -> int:
    """
    Coefficient of x^s in

        x^L P_k(x),

    where L=(ell-k)/2.
    """
    if not in_support(
        k,
        ell,
        s,
    ):
        return 0

    L = support_lower(
        k,
        ell,
    )

    r = s - L

    coeffs = candidate_coefficients(
        candidate,
        k,
    )

    if r < 0 or r >= len(coeffs):
        return 0

    return coeffs[r]


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
    cases: list[Case] = []

    for k in [1, 3, 5, 7, 9, 11, 13]:
        for ell in range(
            k,
            51,
        ):
            L = support_lower(
                k,
                ell,
            )

            U = support_upper(
                k,
                ell,
            )

            for s in range(
                L - 1,
                U + 2,
            ):
                if s < L:
                    location = "BELOW"
                elif s > U:
                    location = "ABOVE"
                elif s == L:
                    location = "LOWER"
                elif s == U:
                    location = "UPPER"
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
# SUPPORT TEST
# ============================================================================

def support_failures(
    candidate: Candidate,
    cases: list[Case],
) -> tuple[int, int]:
    inside_zero = 0
    outside_nonzero = 0

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
            inside_zero += 1

        if not inside and value != 0:
            outside_nonzero += 1

    return (
        inside_zero,
        outside_nonzero,
    )


# ============================================================================
# ENDPOINT TEST
# ============================================================================

def endpoint_failures(
    candidate: Candidate,
) -> tuple[int, int]:
    failures = 0
    checked = 0

    for k in [1, 3, 5, 7, 9, 11, 13]:
        for ell in range(
            k,
            51,
        ):
            L = support_lower(
                k,
                ell,
            )

            U = support_upper(
                k,
                ell,
            )

            first = min(
                (
                    s
                    for s in range(
                        0,
                        U + 3,
                    )
                    if candidate_value(
                        candidate,
                        k,
                        ell,
                        s,
                    ) != 0
                ),
                default=None,
            )

            last = max(
                (
                    s
                    for s in range(
                        0,
                        U + 3,
                    )
                    if candidate_value(
                        candidate,
                        k,
                        ell,
                        s,
                    ) != 0
                ),
                default=None,
            )

            checked += 1

            if first != L or last != U:
                failures += 1

    return (
        failures,
        checked,
    )


# ============================================================================
# k-DEpendence
# ============================================================================

def k_variation_score(
    candidate: Candidate,
) -> tuple[int, int]:
    groups = {}

    for ell in range(
        10,
        41,
    ):
        for s in range(
            0,
            ell + 10,
        ):
            values = []

            for k in [1, 3, 5, 7, 9, 11, 13]:
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
# SHAPE TESTS
# ============================================================================

def coefficient_vector(
    candidate: Candidate,
    k: int,
) -> list[int]:
    return candidate_coefficients(
        candidate,
        k,
    )


def endpoint_normalization_score(
    candidate: Candidate,
) -> tuple[int, int]:
    """
    Check that both endpoint magnitudes equal 1.
    """
    matches = 0
    total = 0

    for k in [1, 3, 5, 7, 9, 11, 13]:
        coeffs = coefficient_vector(
            candidate,
            k,
        )

        total += 1

        if (
            len(coeffs) == k + 1
            and abs(coeffs[0]) == 1
            and abs(coeffs[-1]) == 1
        ):
            matches += 1

    return (
        matches,
        total,
    )


def central_coefficient_score(
    candidate: Candidate,
) -> tuple[int, int]:
    """
    Record whether the polynomial has a nonzero central coefficient.
    """
    nonzero = 0
    total = 0

    for k in [1, 3, 5, 7, 9, 11, 13]:
        coeffs = coefficient_vector(
            candidate,
            k,
        )

        center = k // 2

        total += 1

        if coeffs[center] != 0:
            nonzero += 1

    return (
        nonzero,
        total,
    )


# ============================================================================
# COMPACT VECTOR OUTPUT
# ============================================================================

def print_vectors(
    candidates: list[Candidate],
) -> None:
    print()
    print("=" * 72)
    print("POLYNOMIAL VECTORS FOR SMALL ODD k")
    print("=" * 72)

    selected_k = [1, 3, 5, 7]

    # Print only the four unreversed/base candidates.
    seen = set()

    for candidate in candidates:
        base_name = candidate.polynomial_fn.__name__

        if candidate.reflected:
            continue

        if candidate.global_sign == -1:
            continue

        if base_name in seen:
            continue

        seen.add(base_name)

        print()
        print(candidate.name)

        for k in selected_k:
            print(
                f"  k={k:2d}: "
                f"{candidate_coefficients(candidate, k)}"
            )


# ============================================================================
# RANKING
# ============================================================================

@dataclass(frozen=True)
class Score:
    candidate: Candidate
    inside_zero: int
    outside_nonzero: int
    endpoint_failures: int
    endpoint_total: int
    varying_groups: int
    total_groups: int
    endpoint_normalized: int
    endpoint_normalized_total: int
    central_nonzero: int
    central_total: int


def score_candidate(
    candidate: Candidate,
    cases: list[Case],
) -> Score:
    inside_zero, outside_nonzero = (
        support_failures(
            candidate,
            cases,
        )
    )

    endpoint_fail, endpoint_total = (
        endpoint_failures(
            candidate
        )
    )

    varying, total_groups = (
        k_variation_score(
            candidate
        )
    )

    normalized, normalized_total = (
        endpoint_normalization_score(
            candidate
        )
    )

    central, central_total = (
        central_coefficient_score(
            candidate
        )
    )

    return Score(
        candidate=candidate,
        inside_zero=inside_zero,
        outside_nonzero=outside_nonzero,
        endpoint_failures=endpoint_fail,
        endpoint_total=endpoint_total,
        varying_groups=varying,
        total_groups=total_groups,
        endpoint_normalized=normalized,
        endpoint_normalized_total=normalized_total,
        central_nonzero=central,
        central_total=central_total,
    )


def print_ranking(
    scores: list[Score],
) -> None:
    print()
    print("=" * 72)
    print("DEGREE-k STRUCTURAL RANKING")
    print("=" * 72)

    ordered = sorted(
        scores,
        key=lambda score: (
            score.inside_zero,
            score.outside_nonzero,
            score.endpoint_failures,
        ),
    )

    for score in ordered:
        print()
        print(
            score.candidate.name
        )

        print(
            f"  inside zeros       : "
            f"{score.inside_zero}"
        )

        print(
            f"  outside nonzero   : "
            f"{score.outside_nonzero}"
        )

        print(
            f"  endpoint failures : "
            f"{score.endpoint_failures}/"
            f"{score.endpoint_total}"
        )

        print(
            f"  k-varying groups  : "
            f"{score.varying_groups}/"
            f"{score.total_groups}"
        )

        print(
            f"  endpoint |coef|=1 : "
            f"{score.endpoint_normalized}/"
            f"{score.endpoint_normalized_total}"
        )

        print(
            f"  central nonzero   : "
            f"{score.central_nonzero}/"
            f"{score.central_total}"
        )


# ============================================================================
# DIRECT SUPPORT DEMONSTRATION
# ============================================================================

def print_support_demo() -> None:
    print()
    print("=" * 72)
    print("SUPPORT WIDTH DEMONSTRATION")
    print("=" * 72)

    for k, ell in [
        (3, 8),
        (3, 12),
        (5, 10),
        (7, 15),
        (9, 21),
    ]:
        L = support_lower(
            k,
            ell,
        )

        U = support_upper(
            k,
            ell,
        )

        print(
            f"k={k:2d}, ell={ell:2d}: "
            f"L={L:2d}, U={U:2d}, "
            f"width={U-L}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 190")
    print("Degree-k generating-polynomial reconstruction")
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

    print_support_demo()

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

    print_vectors(
        candidates
    )

    print()
    print("=" * 72)
    print("END EXPERIMENT 190")
    print("=" * 72)


if __name__ == "__main__":
    main()

