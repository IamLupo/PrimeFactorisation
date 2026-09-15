#!/usr/bin/env python3

"""
EXPERIMENT 187
K-dependence / summand-definition audit.

Purpose
-------
Experiment 186 proved that the currently implemented exact summand is
independent of k and therefore collapses to zero for every positive s.

Before attempting another closed-form/sign derivation, this experiment
audits the exact definition itself.

Rules
-----
1. No results.txt or previous output is read.
2. Every function is defined locally.
3. Every test case is generated locally.
4. exact_value() is independent of bulk_value().
5. The experiment explicitly measures k-dependence.
6. Several plausible locations for k-dependence are tested.
7. The script does NOT claim any candidate is the true formula.
   It only reports which structural variants reproduce nontrivial
   k-dependent behavior.

IMPORTANT
---------
The ACTUAL Experiment-183 summand should eventually replace the marked
REFERENCE SUMMAND section. Until then this script is an audit of the
current transcription.
"""

from __future__ import annotations

from dataclasses import dataclass
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
# REFERENCE SUMMAND
# ============================================================================
#
# This is the summand currently responsible for the zero collapse in 186.
#
# DO NOT use this as evidence for the real formula.
# It is retained here so the audit can demonstrate exactly what is wrong.
# ============================================================================

def reference_summand(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    Current transcription from Experiment 186.

    Notice deliberately:
        k does not occur.

    This is the defect being audited.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - j)
    )


def reference_exact_value(
    k: int,
    ell: int,
    s: int,
) -> int:
    """Exact value using the currently transcribed summand."""
    if not in_support(k, ell, s):
        return 0

    total = 0

    for j in range(ell + 1):
        total += reference_summand(
            k,
            ell,
            s,
            j,
        )

    return total


# ============================================================================
# K-DEPENDENCE AUDIT
# ============================================================================

@dataclass(frozen=True)
class Case:
    k: int
    ell: int
    s: int


def generate_cases() -> list[Case]:
    """Generate fresh cases."""
    cases: list[Case] = []

    for ell in range(8, 31):
        for k in [1, 3, 5, 7, 9]:
            lo = support_lower(k, ell)
            hi = support_upper(k, ell)

            if lo > hi:
                continue

            for s in range(lo, hi + 1):
                cases.append(
                    Case(
                        k=k,
                        ell=ell,
                        s=s,
                    )
                )

    return cases


def compare_same_ell_s_different_k(
    ell: int,
    s: int,
) -> list[tuple[int, int]]:
    """
    Evaluate all admissible k for fixed ell,s.

    Returns (k,value).
    """
    values = []

    for k in [1, 3, 5, 7, 9]:
        if in_support(k, ell, s):
            values.append(
                (
                    k,
                    reference_exact_value(
                        k,
                        ell,
                        s,
                    ),
                )
            )

    return values


def count_k_variation(
    cases: list[Case],
) -> tuple[int, int]:
    """
    Count fixed-(ell,s) groups where the value changes with k.

    Returns:
        varying_groups,
        total_groups
    """
    groups: dict[tuple[int, int], list[int]] = {}

    for case in cases:
        key = (case.ell, case.s)

        groups.setdefault(
            key,
            [],
        ).append(
            reference_exact_value(
                case.k,
                case.ell,
                case.s,
            )
        )

    varying = 0

    for values in groups.values():
        if len(set(values)) > 1:
            varying += 1

    return varying, len(groups)


def print_k_audit(
    cases: list[Case],
) -> None:
    varying, total = count_k_variation(
        cases
    )

    print()
    print("=" * 72)
    print("K-DEPENDENCE AUDIT")
    print("=" * 72)

    print(f"groups checked     : {total}")
    print(f"k-varying groups   : {varying}")
    print(
        f"k-independent     : {total - varying}"
    )

    if varying == 0:
        print()
        print(
            "CONFIRMED: the current exact implementation "
            "contains no effective k-dependence."
        )


# ============================================================================
# SYMBOLIC-STRUCTURE TESTS
# ============================================================================
#
# These are NOT claims about the true formula.
# They test common ways a parameter k could enter an alternating binomial
# convolution.
# ============================================================================

@dataclass(frozen=True)
class Variant:
    name: str
    summand: Callable[[int, int, int, int], int]


def variant_shift_first_binomial(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k enters the first upper binomial parameter.
    """
    return (
        sgn(j)
        * C(ell + k, j)
        * C(ell + k - j, s - j)
    )


def variant_shift_second_binomial(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k enters the second upper binomial parameter.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell + k - j, s - j)
    )


def variant_shift_lower_argument(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k shifts the target coefficient index.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - k - j)
    )


def variant_shift_j_by_k(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k shifts the summation index in the second factor.
    """
    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - k - j)
    )


def variant_k_dependent_sign(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k enters only through an alternating sign.
    """
    return (
        sgn(j + k)
        * C(ell, j)
        * C(ell - j, s - j)
    )


def variant_k_dependent_truncation(
    k: int,
    ell: int,
    s: int,
    j: int,
) -> int:
    """
    k changes the summation range rather than the summand.
    This is handled separately below, but the function is included
    so every candidate remains explicit.
    """
    if j > k:
        return 0

    return (
        sgn(j)
        * C(ell, j)
        * C(ell - j, s - j)
    )


def make_variants() -> list[Variant]:
    return [
        Variant(
            "k in first upper binomial",
            variant_shift_first_binomial,
        ),
        Variant(
            "k in second upper binomial",
            variant_shift_second_binomial,
        ),
        Variant(
            "k shifts target index",
            variant_shift_lower_argument,
        ),
        Variant(
            "k shifts j/index",
            variant_shift_j_by_k,
        ),
        Variant(
            "k in alternating sign",
            variant_k_dependent_sign,
        ),
        Variant(
            "k truncates j-range",
            variant_k_dependent_truncation,
        ),
    ]


def variant_value(
    variant: Variant,
    k: int,
    ell: int,
    s: int,
) -> int:
    """Evaluate one candidate structural variant."""
    if not in_support(k, ell, s):
        return 0

    if variant.name == "k truncates j-range":
        hi = min(ell, k)
    else:
        hi = ell

    total = 0

    for j in range(hi + 1):
        total += variant.summand(
            k,
            ell,
            s,
            j,
        )

    return total


# ============================================================================
# VARIANT AUDIT
# ============================================================================

def variant_k_variation_score(
    variant: Variant,
    cases: list[Case],
) -> tuple[int, int, int]:
    """
    Return:
        varying groups,
        nonzero cases,
        total groups
    """
    groups: dict[tuple[int, int], list[int]] = {}

    for case in cases:
        key = (
            case.ell,
            case.s,
        )

        value = variant_value(
            variant,
            case.k,
            case.ell,
            case.s,
        )

        groups.setdefault(
            key,
            [],
        ).append(value)

    varying = 0

    for values in groups.values():
        if len(set(values)) > 1:
            varying += 1

    nonzero = sum(
        variant_value(
            variant,
            case.k,
            case.ell,
            case.s,
        ) != 0
        for case in cases
    )

    return (
        varying,
        nonzero,
        len(groups),
    )


def print_variant_audit(
    cases: list[Case],
) -> None:
    print()
    print("=" * 72)
    print("K-DEPENDENCE OF STRUCTURAL VARIANTS")
    print("=" * 72)

    for variant in make_variants():
        varying, nonzero, groups = (
            variant_k_variation_score(
                variant,
                cases,
            )
        )

        print()
        print(
            variant.name
        )
        print(
            f"  varying groups : {varying}/{groups}"
        )
        print(
            f"  nonzero cases  : {nonzero}/{len(cases)}"
        )


# ============================================================================
# DIRECT EXAMPLES
# ============================================================================

def print_fixed_ell_s_examples() -> None:
    """
    Show how the current reference implementation behaves when k changes.
    """
    examples = [
        (12, 5),
        (12, 6),
        (16, 7),
        (16, 8),
        (20, 9),
        (20, 10),
    ]

    print()
    print("=" * 72)
    print("FIXED (ELL,S), VARY K")
    print("=" * 72)

    for ell, s in examples:
        values = compare_same_ell_s_different_k(
            ell,
            s,
        )

        print()
        print(
            f"ell={ell}, s={s}"
        )

        for k, value in values:
            print(
                f"  k={k:2d} -> {value:12d}"
            )


# ============================================================================
# SUPPORT / NONZERO COMPATIBILITY
# ============================================================================

def print_support_nonzero_audit(
    cases: list[Case],
) -> None:
    """
    Compare support membership with whether the current exact definition
    is nonzero.
    """
    inside = 0
    inside_nonzero = 0

    for case in cases:
        if in_support(
            case.k,
            case.ell,
            case.s,
        ):
            inside += 1

            if reference_exact_value(
                case.k,
                case.ell,
                case.s,
            ) != 0:
                inside_nonzero += 1

    print()
    print("=" * 72)
    print("SUPPORT / NONZERO COMPATIBILITY")
    print("=" * 72)

    print(
        f"inside support : {inside}"
    )

    print(
        f"nonzero exact  : {inside_nonzero}"
    )

    if inside:
        print(
            f"nonzero rate   : "
            f"{100.0 * inside_nonzero / inside:.3f}%"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 72)
    print("EXPERIMENT 187")
    print("K-dependence / exact-definition audit")
    print("=" * 72)

    print()
    print(
        "All data are generated locally."
    )
    print(
        "No previous result file is read."
    )

    cases = generate_cases()

    print()
    print(
        f"Generated cases: {len(cases)}"
    )

    print_k_audit(
        cases
    )

    print_support_nonzero_audit(
        cases
    )

    print_fixed_ell_s_examples()

    print_variant_audit(
        cases
    )

    print()
    print("=" * 72)
    print("END EXPERIMENT 187")
    print("=" * 72)


if __name__ == "__main__":
    main()

