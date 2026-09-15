#!/usr/bin/env python3

"""
========================================================================
EXPERIMENT 191
Degree-k coefficient recovery from exact anchor constraints
========================================================================

Purpose
-------
Experiment 190 established that:

    support width = k

is compatible with degree-k polynomial families.

However, it did NOT establish that the true object is simply

    x^L (1 +/- x)^k.

The current experiment therefore asks a narrower question:

    Can the known exact anchor values be represented as

        exact(k, ell, s)
          = A(k, ell) * P_k[r],

    where

        r = s - L(k, ell),

    and P_k is one of several degree-k candidate coefficient vectors?

The scalar A(k,ell) is allowed to depend on k and ell, but NOT on s.

This is a strong test:
    - if the ratios vary with s, the candidate polynomial shape is wrong;
    - if ratios are constant, that candidate remains plausible.

IMPORTANT
---------
1. No results.txt is read.
2. No previous script output is read.
3. All test cases are constructed in this file.
4. All candidate functions are defined here.
5. The small set of exact anchors below is copied explicitly from the
   conversation, not loaded from any external file.
6. Those anchors are used only as constraints; all candidate grids are
   freshly generated.
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
# SUPPORT MODEL USED FOR THE RECOVERY TEST
# ============================================================================

def support_lower(k: int, ell: int) -> int:
    return (ell - k) // 2


def support_upper(k: int, ell: int) -> int:
    return (ell + k) // 2


def support_index(
    k: int,
    ell: int,
    s: int,
) -> int:
    """
    Index inside the degree-k coefficient vector.
    """
    return (
        s
        - support_lower(k, ell)
    )


# ============================================================================
# KNOWN EXACT ANCHORS
# ============================================================================
#
# These are explicit values that appeared in the conversation.
#
# They are NOT read from results.txt.
#
# More anchors can be added here as they become available.
# ============================================================================

@dataclass(frozen=True)
class Anchor:
    k: int
    ell: int
    s: int
    exact: int
    label: str


def generate_anchor_data() -> list[Anchor]:
    """
    Generate the explicit anchor dataset.

    The values below are taken from the previously reported exact values.
    """
    return [
        Anchor(
            k=1,
            ell=10,
            s=2,
            exact=27,
            label="known-183-a",
        ),
        Anchor(
            k=3,
            ell=20,
            s=6,
            exact=935,
            label="known-183-b",
        ),
        Anchor(
            k=7,
            ell=40,
            s=15,
            exact=-1797818,
            label="known-183-c",
        ),
    ]


# ============================================================================
# DEGREE-k POLYNOMIAL FAMILIES
# ============================================================================

@dataclass(frozen=True)
class PolynomialCandidate:
    name: str
    coefficient_fn: Callable[
        [int, int],
        int,
    ]


def plus_binomial(
    k: int,
    r: int,
) -> int:
    return C(k, r)


def minus_binomial(
    k: int,
    r: int,
) -> int:
    return (
        sgn(r)
        * C(k, r)
    )


def reversed_plus_binomial(
    k: int,
    r: int,
) -> int:
    return C(k, k - r)


def reversed_minus_binomial(
    k: int,
    r: int,
) -> int:
    return (
        sgn(k - r)
        * C(k, k - r)
    )


def central_difference_binomial(
    k: int,
    r: int,
) -> int:
    """
    First simple deformation:
        C(k,r) - C(k,r-1)
    """
    return (
        C(k, r)
        - C(k, r - 1)
    )


def central_sum_binomial(
    k: int,
    r: int,
) -> int:
    """
    First simple convolution deformation:
        C(k,r) + C(k,r-1)
    """
    return (
        C(k, r)
        + C(k, r - 1)
    )


def alternating_difference(
    k: int,
    r: int,
) -> int:
    return (
        sgn(r)
        * (
            C(k, r)
            - C(k, r - 1)
        )
    )


def alternating_sum(
    k: int,
    r: int,
) -> int:
    return (
        sgn(r)
        * (
            C(k, r)
            + C(k, r - 1)
        )
    )


def quadratic_binomial(
    k: int,
    r: int,
) -> int:
    """
    Second-neighbor binomial deformation:
        C(k,r) - C(k,r-2)
    """
    return (
        C(k, r)
        - C(k, r - 2)
    )


def quadratic_sum_binomial(
    k: int,
    r: int,
) -> int:
    return (
        C(k, r)
        + C(k, r - 2)
    )


def make_polynomial_candidates() -> list[PolynomialCandidate]:
    return [
        PolynomialCandidate(
            "(1+x)^k",
            plus_binomial,
        ),
        PolynomialCandidate(
            "(1-x)^k",
            minus_binomial,
        ),
        PolynomialCandidate(
            "reversed (1+x)^k",
            reversed_plus_binomial,
        ),
        PolynomialCandidate(
            "reversed (1-x)^k",
            reversed_minus_binomial,
        ),
        PolynomialCandidate(
            "C(k,r)-C(k,r-1)",
            central_difference_binomial,
        ),
        PolynomialCandidate(
            "C(k,r)+C(k,r-1)",
            central_sum_binomial,
        ),
        PolynomialCandidate(
            "(-1)^r[C(k,r)-C(k,r-1)]",
            alternating_difference,
        ),
        PolynomialCandidate(
            "(-1)^r[C(k,r)+C(k,r-1)]",
            alternating_sum,
        ),
        PolynomialCandidate(
            "C(k,r)-C(k,r-2)",
            quadratic_binomial,
        ),
        PolynomialCandidate(
            "C(k,r)+C(k,r-2)",
            quadratic_sum_binomial,
        ),
    ]


# ============================================================================
# SCALAR-COMPATIBILITY TEST
# ============================================================================

@dataclass(frozen=True)
class AnchorFit:
    candidate: PolynomialCandidate
    compatible: bool
    reason: str
    ratios: tuple[tuple[int, int, Fraction | None], ...]


def fit_single_anchor(
    candidate: PolynomialCandidate,
    anchor: Anchor,
) -> Fraction | None:
    """
    Compute the normalization A implied by one anchor.

    exact = A * P_k[r]
    """
    r = support_index(
        anchor.k,
        anchor.ell,
        anchor.s,
    )

    coefficient = candidate.coefficient_fn(
        anchor.k,
        r,
    )

    if coefficient == 0:
        return None

    return Fraction(
        anchor.exact,
        coefficient,
    )


def test_anchor_compatibility(
    candidate: PolynomialCandidate,
    anchors: list[Anchor],
) -> AnchorFit:
    """
    Test whether each anchor has a finite scalar A(k,ell).

    Since all three anchors have different (k,ell), this first test checks
    only whether the proposed polynomial can produce a nonzero coefficient
    at the observed support index.

    A stronger grouped test is performed later.
    """
    ratios = []

    for anchor in anchors:
        r = support_index(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

        coefficient = candidate.coefficient_fn(
            anchor.k,
            r,
        )

        ratio = fit_single_anchor(
            candidate,
            anchor,
        )

        ratios.append(
            (
                anchor.k,
                r,
                ratio,
            )
        )

        if coefficient == 0:
            return AnchorFit(
                candidate=candidate,
                compatible=False,
                reason=(
                    f"zero coefficient at "
                    f"k={anchor.k}, r={r}"
                ),
                ratios=tuple(ratios),
            )

    return AnchorFit(
        candidate=candidate,
        compatible=True,
        reason="all anchor coefficients nonzero",
        ratios=tuple(ratios),
    )


# ============================================================================
# GENERATED MULTI-s TESTS
# ============================================================================

def generate_same_keL_tests() -> list[tuple[int, int]]:
    """
    Generate fresh (k,ell) pairs on which multiple s values are tested.

    These do not rely on previous outputs.
    """
    return [
        (1, 10),
        (3, 20),
        (5, 20),
        (7, 40),
        (9, 30),
        (11, 42),
    ]


def generate_s_grid(
    k: int,
    ell: int,
) -> list[int]:
    """
    Generate every integer s in the assumed support.
    """
    lo = support_lower(
        k,
        ell,
    )

    hi = support_upper(
        k,
        ell,
    )

    return list(
        range(
            lo,
            hi + 1,
        )
    )


# ============================================================================
# SHAPE / RATIO DIAGNOSTIC
# ============================================================================

def print_candidate_vectors(
    candidates: list[PolynomialCandidate],
) -> None:
    print()
    print("=" * 72)
    print("SMALL DEGREE-K CANDIDATE VECTORS")
    print("=" * 72)

    for candidate in candidates:
        print()
        print(candidate.name)

        for k in [1, 3, 5, 7]:
            vector = [
                candidate.coefficient_fn(
                    k,
                    r,
                )
                for r in range(k + 1)
            ]

            print(
                f"  k={k}: {vector}"
            )


def print_anchor_ratios(
    candidates: list[PolynomialCandidate],
    anchors: list[Anchor],
) -> None:
    print()
    print("=" * 72)
    print("KNOWN ANCHOR COMPATIBILITY")
    print("=" * 72)

    for candidate in candidates:
        fit = test_anchor_compatibility(
            candidate,
            anchors,
        )

        print()
        print(candidate.name)
        print(
            f"  compatible: {fit.compatible}"
        )
        print(
            f"  reason    : {fit.reason}"
        )

        for (
            k,
            r,
            ratio,
        ) in fit.ratios:
            print(
                f"  k={k:2d}, r={r:4d}, "
                f"A={ratio}"
            )


# ============================================================================
# PARAMETER-FITTING TEST
# ============================================================================

def test_single_keL_shape(
    candidate: PolynomialCandidate,
    k: int,
    ell: int,
    normalization: Fraction,
) -> list[tuple[int, int, int]]:
    """
    Generate all s for a (k,ell) pair and produce predicted values.
    """
    rows = []

    for s in generate_s_grid(
        k,
        ell,
    ):
        r = support_index(
            k,
            ell,
            s,
        )

        coefficient = candidate.coefficient_fn(
            k,
            r,
        )

        predicted = (
            normalization
            * coefficient
        )

        rows.append(
            (
                s,
                r,
                int(predicted),
            )
        )

    return rows


def print_generated_shapes(
    candidates: list[PolynomialCandidate],
) -> None:
    """
    Display freshly generated coefficient shapes for several (k,ell)
    pairs.

    No previous exact data are used here.
    """
    pairs = generate_same_keL_tests()

    print()
    print("=" * 72)
    print("FRESHLY GENERATED SHAPES")
    print("=" * 72)

    for candidate in candidates[:4]:
        print()
        print(candidate.name)

        for k, ell in pairs:
            coeffs = [
                candidate.coefficient_fn(
                    k,
                    r,
                )
                for r in range(k + 1)
            ]

            print(
                f"  (k,ell)=({k},{ell}) "
                f"coeff={coeffs}"
            )


# ============================================================================
# SUPPORT INDEX AUDIT
# ============================================================================

def print_anchor_indices(
    anchors: list[Anchor],
) -> None:
    print()
    print("=" * 72)
    print("ANCHOR SUPPORT-INDEX AUDIT")
    print("=" * 72)

    for anchor in anchors:
        L = support_lower(
            anchor.k,
            anchor.ell,
        )

        U = support_upper(
            anchor.k,
            anchor.ell,
        )

        r = support_index(
            anchor.k,
            anchor.ell,
            anchor.s,
        )

        print(
            f"{anchor.label:12s} "
            f"k={anchor.k:2d} "
            f"ell={anchor.ell:2d} "
            f"s={anchor.s:3d} "
            f"L={L:3d} "
            f"U={U:3d} "
            f"r={r:3d} "
            f"exact={anchor.exact}"
        )


# ============================================================================
# CONSISTENCY CHECK
# ============================================================================

def verify_degree_k_property(
    candidates: list[PolynomialCandidate],
) -> None:
    """
    Verify that each polynomial really has degree <= k and endpoint
    coefficients are available.
    """
    print()
    print("=" * 72)
    print("DEGREE CHECK")
    print("=" * 72)

    for candidate in candidates:
        failures = 0

        for k in range(1, 14, 2):
            vector = [
                candidate.coefficient_fn(
                    k,
                    r,
                )
                for r in range(k + 1)
            ]

            if len(vector) != k + 1:
                failures += 1

        print(
            f"{candidate.name:40s} "
            f"degree-shape failures={failures}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 191")
    print("Degree-k coefficient recovery from exact anchors")
    print("=" * 72)

    print()
    print(
        "All test cases are generated locally."
    )
    print(
        "No results.txt or previous experiment output is read."
    )

    anchors = generate_anchor_data()
    candidates = make_polynomial_candidates()

    print()
    print(
        f"Anchor constraints: {len(anchors)}"
    )

    print(
        "Fresh (k,ell) test pairs:",
        len(generate_same_keL_tests()),
    )

    print_anchor_indices(
        anchors
    )

    verify_degree_k_property(
        candidates
    )

    print_candidate_vectors(
        candidates
    )

    print_anchor_ratios(
        candidates,
        anchors,
    )

    print_generated_shapes(
        candidates
    )

    print()
    print("=" * 72)
    print("END EXPERIMENT 191")
    print("=" * 72)


if __name__ == "__main__":
    main()

