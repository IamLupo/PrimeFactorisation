#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 196
STAGED PRE-CONVOLUTION k-PLACEMENT SEARCH
==============================================================================

Purpose
-------
Experiment 195 established that the current reference convolution is
k-independent after collapse. Therefore this experiment searches for
candidate PRE-CONVOLUTION structures in which k enters before the final
binomial cancellation.

This is deliberately a STAGED search.

It does NOT attempt a massive unrestricted affine-family enumeration.
Candidates are filtered in the following order:

    1. exact anchor agreement
    2. fixed-(ell,s), varying-k behavior
    3. support compatibility
    4. broader diagnostic scan

Only survivors reach the more expensive stages.

Important
---------
This experiment is DIAGNOSTIC ONLY.

A surviving family is not accepted as the theorem unless it can later be
derived explicitly from the n = p*q construction.

==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from typing import Callable, Iterable


# ============================================================================
# 1. BASIC DATA
# ============================================================================

@dataclass(frozen=True)
class Anchor:
    label: str
    k: int
    ell: int
    s: int
    exact: int


ANCHORS = [
    Anchor("A1", 1, 10, 2, 27),
    Anchor("A2", 3, 20, 6, 935),
    Anchor("A3", 7, 40, 15, -1797818),
    Anchor("B1", 3, 7, 3, -3),
    Anchor("B2", 3, 9, 4, 9),
    Anchor("B3", 3, 11, 5, -16),
    Anchor("B4", 5, 11, 5, -5),
    Anchor("B5", 5, 19, 9, -196),
    Anchor("B6", 9, 21, 11, -84),
]


FIXED_ELL_S_GROUPS = [
    (8, 3),
    (8, 4),
    (9, 4),
    (10, 4),
    (10, 5),
    (12, 5),
    (12, 6),
    (20, 6),
    (20, 7),
    (40, 15),
]


# ============================================================================
# 2. SAFE BINOMIAL
# ============================================================================

def C(n: int, r: int) -> int:
    """
    Integer binomial coefficient with zero outside the natural range.
    """
    if n < 0 or r < 0 or r > n:
        return 0
    return comb(n, r)


# ============================================================================
# 3. EXACT REFERENCE CONVOLUTION
# ============================================================================

def exact_reference(k: int, ell: int, s: int) -> int:
    """
    Reference convolution from the preceding experiments.

        sum_j (-1)^j C(ell,j) C(ell-j, s-j)

    The identity

        C(ell,j) C(ell-j,s-j)
        = C(ell,s) C(s,j)

    causes the usual alternating cancellation.

    This function is retained as a control and is NOT expected to recover
    the observed anchor values.
    """
    total = 0

    for j in range(0, s + 1):
        total += (-1) ** j * C(ell, j) * C(ell - j, s - j)

    return total


# ============================================================================
# 4. PRE-CONVOLUTION CANDIDATE DEFINITIONS
# ============================================================================

@dataclass(frozen=True)
class Candidate:
    name: str
    fn: Callable[[int, int, int], int]
    description: str


def convolution_sum(
    k: int,
    ell: int,
    s: int,
    term_fn: Callable[[int, int, int, int], int],
    j_min: int = 0,
    j_max: int | None = None,
) -> int:
    """
    Generic pre-convolution summation.

    term_fn(k, ell, s, j) returns the j-th summand.

    The summation range is conservatively chosen to include all potentially
    nonzero terms while allowing candidates to implement their own zeros.
    """
    if j_max is None:
        j_max = max(ell + abs(k) + 5, s + abs(k) + 5)

    total = 0

    for j in range(j_min, j_max + 1):
        total += term_fn(k, ell, s, j)

    return total


# ----------------------------------------------------------------------------
# Candidate A
# ----------------------------------------------------------------------------

def cand_A_term(k: int, ell: int, s: int, j: int) -> int:
    """
    A:
        C(ell, j+k) C(ell-j-k, s-j-k)
    """
    return (
        (-1) ** j
        * C(ell, j + k)
        * C(ell - j - k, s - j - k)
    )


def cand_A(k: int, ell: int, s: int) -> int:
    return convolution_sum(k, ell, s, cand_A_term)


# ----------------------------------------------------------------------------
# Candidate B
# ----------------------------------------------------------------------------

def cand_B_term(k: int, ell: int, s: int, j: int) -> int:
    """
    B:
        C(ell+k, j) C(ell-j, s-j)
    """
    return (
        (-1) ** j
        * C(ell + k, j)
        * C(ell - j, s - j)
    )


def cand_B(k: int, ell: int, s: int) -> int:
    return convolution_sum(k, ell, s, cand_B_term)


# ----------------------------------------------------------------------------
# Candidate C
# ----------------------------------------------------------------------------

def cand_C_term(k: int, ell: int, s: int, j: int) -> int:
    """
    C:
        C(ell,j) C(ell-j+k, s-j)
    """
    return (
        (-1) ** j
        * C(ell, j)
        * C(ell - j + k, s - j)
    )


def cand_C(k: int, ell: int, s: int) -> int:
    return convolution_sum(k, ell, s, cand_C_term)


# ----------------------------------------------------------------------------
# Candidate D
# ----------------------------------------------------------------------------

def cand_D_term(k: int, ell: int, s: int, j: int) -> int:
    """
    D:
        C(ell,j) C(ell-j, s-j+k)
    """
    return (
        (-1) ** j
        * C(ell, j)
        * C(ell - j, s - j + k)
    )


def cand_D(k: int, ell: int, s: int) -> int:
    return convolution_sum(k, ell, s, cand_D_term)


# ----------------------------------------------------------------------------
# Candidate E
# ----------------------------------------------------------------------------

def cand_E_term(k: int, ell: int, s: int, j: int) -> int:
    """
    E:
        C(ell,j+k) C(ell-j, s-j)
    """
    return (
        (-1) ** j
        * C(ell, j + k)
        * C(ell - j, s - j)
    )


def cand_E(k: int, ell: int, s: int) -> int:
    return convolution_sum(k, ell, s, cand_E_term)


# ----------------------------------------------------------------------------
# Candidate F
# ----------------------------------------------------------------------------

def cand_F_term(k: int, ell: int, s: int, j: int) -> int:
    """
    F:
        C(ell+k,j+k) C(ell-j, s-j)
    """
    return (
        (-1) ** j
        * C(ell + k, j + k)
        * C(ell - j, s - j)
    )


def cand_F(k: int, ell: int, s: int) -> int:
    return convolution_sum(k, ell, s, cand_F_term)


BASE_CANDIDATES = [
    Candidate(
        "A_second_index_shift",
        cand_A,
        "C(ell,j+k) C(ell-j-k,s-j-k)",
    ),
    Candidate(
        "B_first_upper_plus_k",
        cand_B,
        "C(ell+k,j) C(ell-j,s-j)",
    ),
    Candidate(
        "C_second_upper_plus_k",
        cand_C,
        "C(ell,j) C(ell-j+k,s-j)",
    ),
    Candidate(
        "D_lower_index_plus_k",
        cand_D,
        "C(ell,j) C(ell-j,s-j+k)",
    ),
    Candidate(
        "E_first_index_shift",
        cand_E,
        "C(ell,j+k) C(ell-j,s-j)",
    ),
    Candidate(
        "F_both_first_upper_and_index",
        cand_F,
        "C(ell+k,j+k) C(ell-j,s-j)",
    ),
]


# ============================================================================
# 5. SMALL SHIFT EXTENSIONS
# ============================================================================

@dataclass(frozen=True)
class ShiftCandidate:
    name: str
    fn: Callable[[int, int, int], int]
    description: str


def make_shifted_candidate(
    base_name: str,
    description: str,
    builder: Callable[[int, int, int, int, int], int],
    delta_k: int,
    delta_j: int,
) -> ShiftCandidate:
    """
    Generic small-shift wrapper.

    This is intentionally restricted to small shifts:
        delta_k, delta_j in {-1, 0, +1}
    """

    def fn(k: int, ell: int, s: int) -> int:
        total = 0
        j_max = max(ell + abs(k) + 4, s + abs(k) + 4)

        for j in range(0, j_max + 1):
            jj = j + delta_j
            kk = k + delta_k
            total += builder(kk, ell, s, jj, j)

        return total

    return ShiftCandidate(
        name=(
            f"{base_name}"
            f"[dk={delta_k:+d},dj={delta_j:+d}]"
        ),
        fn=fn,
        description=description,
    )


# ============================================================================
# 6. SUPPORT MODEL
# ============================================================================

def support_bounds(ell: int, k: int) -> tuple[int, int]:
    """
    Current experimentally favored support geometry from Experiments 192-193:

        L = (ell - 3k)/4
        U = (ell + k)/2

    We use a conservative integer support enclosure.

    This is a DIAGNOSTIC support model, not yet a theorem.
    """
    L_real = (ell - 3 * k) / 4.0
    U_real = (ell + k) / 2.0

    L = int(L_real.__ceil__())
    U = int(U_real.__floor__())

    return L, U


def in_support(ell: int, k: int, s: int) -> bool:
    L, U = support_bounds(ell, k)
    return L <= s <= U


# ============================================================================
# 7. STAGE 1: ANCHOR FILTER
# ============================================================================

def anchor_filter(candidate: Candidate) -> bool:
    """
    Candidate must reproduce ALL nine anchor values.
    """
    for a in ANCHORS:
        value = candidate.fn(a.k, a.ell, a.s)

        if value != a.exact:
            return False

    return True


def anchor_report(candidate: Candidate) -> tuple[int, int]:
    matches = 0

    for a in ANCHORS:
        value = candidate.fn(a.k, a.ell, a.s)
        if value == a.exact:
            matches += 1

    return matches, len(ANCHORS)


# ============================================================================
# 8. STAGE 2: FIXED-(ell,s), VARY-k TEST
# ============================================================================

def k_dependence_filter(candidate: Candidate) -> bool:
    """
    The current exact implementation was k-independent for fixed (ell,s),
    which is precisely the behavior we are trying to escape.

    A structurally useful candidate should produce at least one genuine
    k-dependent group among the test families.
    """

    for ell, s in FIXED_ELL_S_GROUPS:
        values = []

        # odd k only, matching the branch under investigation
        for k in range(1, min(ell, 11) + 1, 2):
            values.append(candidate.fn(k, ell, s))

        if len(set(values)) > 1:
            return True

    return False


def k_dependence_score(candidate: Candidate) -> int:
    score = 0

    for ell, s in FIXED_ELL_S_GROUPS:
        values = [
            candidate.fn(k, ell, s)
            for k in range(1, min(ell, 11) + 1, 2)
        ]

        if len(set(values)) > 1:
            score += 1

    return score


# ============================================================================
# 9. STAGE 3: SUPPORT FILTER
# ============================================================================

def support_filter(candidate: Candidate) -> bool:
    """
    A candidate passes if it does not produce nonzero values far outside the
    currently observed support.

    We do NOT demand exact values here; this is only a structural screen.
    """

    checks = 0

    for k in range(1, 12, 2):
        for ell in range(max(k, 3), 25):
            if ell < k:
                continue

            L, U = support_bounds(ell, k)

            # Test a small band around support.
            for s in range(max(0, L - 2), U + 3):
                value = candidate.fn(k, ell, s)

                # Outside-support points should be zero.
                if s < L or s > U:
                    if value != 0:
                        return False

                checks += 1

                # Keep this diagnostic bounded.
                if checks >= 500:
                    return True

    return True


# ============================================================================
# 10. STAGE 4: BROADER ANCHOR-NEIGHBOR CHECK
# ============================================================================

def neighborhood_score(candidate: Candidate) -> tuple[int, int, int]:
    """
    Compare candidate against the anchor neighborhood.

    Returns:
        matches, total, absolute_error
    """

    matches = 0
    total = 0
    error = 0

    for a in ANCHORS:
        for dk in (-2, 0, 2):
            kk = a.k + dk

            if kk <= 0 or kk % 2 == 0:
                continue

            if kk > a.ell:
                continue

            # Keep (ell,s) fixed while varying k.
            target = candidate.fn(kk, a.ell, a.s)

            # Only exact anchor itself has a known target.
            if kk == a.k:
                total += 1
                if target == a.exact:
                    matches += 1
                error += abs(target - a.exact)

    return matches, total, error


# ============================================================================
# 11. RESULT RECORD
# ============================================================================

@dataclass
class Survivor:
    name: str
    description: str
    anchor_matches: int
    anchor_total: int
    k_dependence_score: int
    neighborhood_matches: int
    neighborhood_total: int
    neighborhood_error: int


# ============================================================================
# 12. STAGED SEARCH
# ============================================================================

def run_staged_search() -> list[Survivor]:
    print()
    print("=" * 78)
    print("STAGED SEARCH")
    print("=" * 78)

    survivors: list[Survivor] = []

    # ------------------------------------------------------------------
    # Stage 1
    # ------------------------------------------------------------------
    print()
    print("-" * 78)
    print("STAGE 1: ANCHOR FILTER")
    print("-" * 78)

    stage1: list[Candidate] = []

    for candidate in BASE_CANDIDATES:
        matches, total = anchor_report(candidate)

        print(
            f"{candidate.name:32s} "
            f"{matches:2d}/{total}"
        )

        if matches == total:
            stage1.append(candidate)

    print(f"\nStage 1 survivors = {len(stage1)}")

    # ------------------------------------------------------------------
    # Stage 2
    # ------------------------------------------------------------------
    print()
    print("-" * 78)
    print("STAGE 2: k-DEPENDENCE FILTER")
    print("-" * 78)

    stage2: list[Candidate] = []

    for candidate in stage1:
        score = k_dependence_score(candidate)

        print(
            f"{candidate.name:32s} "
            f"k-dependent groups={score}/{len(FIXED_ELL_S_GROUPS)}"
        )

        if score > 0:
            stage2.append(candidate)

    print(f"\nStage 2 survivors = {len(stage2)}")

    # ------------------------------------------------------------------
    # Stage 3
    # ------------------------------------------------------------------
    print()
    print("-" * 78)
    print("STAGE 3: SUPPORT FILTER")
    print("-" * 78)

    stage3: list[Candidate] = []

    for candidate in stage2:
        passed = support_filter(candidate)

        print(
            f"{candidate.name:32s} "
            f"{'PASS' if passed else 'FAIL'}"
        )

        if passed:
            stage3.append(candidate)

    print(f"\nStage 3 survivors = {len(stage3)}")

    # ------------------------------------------------------------------
    # Stage 4
    # ------------------------------------------------------------------
    print()
    print("-" * 78)
    print("STAGE 4: NEIGHBORHOOD SCORE")
    print("-" * 78)

    for candidate in stage3:
        matches, total, error = neighborhood_score(candidate)
        kd = k_dependence_score(candidate)

        print(
            f"{candidate.name:32s} "
            f"anchors={matches}/{total} "
            f"kdep={kd} "
            f"abs_error={error}"
        )

        survivors.append(
            Survivor(
                name=candidate.name,
                description=candidate.description,
                anchor_matches=matches,
                anchor_total=total,
                k_dependence_score=kd,
                neighborhood_matches=matches,
                neighborhood_total=total,
                neighborhood_error=error,
            )
        )

    return survivors


# ============================================================================
# 13. SEPARATE SHIFT SEARCH
# ============================================================================

def run_small_shift_search() -> list[ShiftCandidate]:
    """
    Only reached after the base structural families have been tested.

    We deliberately keep this search tiny:
        dk in {-1,0,+1}
        dj in {-1,0,+1}

    The shifts are applied to the basic structural families.
    """

    print()
    print("=" * 78)
    print("SMALL SHIFT SEARCH")
    print("=" * 78)

    shifted: list[ShiftCandidate] = []

    # Builder functions.
    #
    # Each builder receives:
    #   kk, ell, s, jj, original_j
    #
    # and constructs the summand.

    def A_builder(
        kk: int, ell: int, s: int, jj: int, original_j: int
    ) -> int:
        return (
            (-1) ** original_j
            * C(ell, jj + kk)
            * C(ell - jj - kk, s - jj - kk)
        )

    def B_builder(
        kk: int, ell: int, s: int, jj: int, original_j: int
    ) -> int:
        return (
            (-1) ** original_j
            * C(ell + kk, jj)
            * C(ell - jj, s - jj)
        )

    def C_builder(
        kk: int, ell: int, s: int, jj: int, original_j: int
    ) -> int:
        return (
            (-1) ** original_j
            * C(ell, jj)
            * C(ell - jj + kk, s - jj)
        )

    def D_builder(
        kk: int, ell: int, s: int, jj: int, original_j: int
    ) -> int:
        return (
            (-1) ** original_j
            * C(ell, jj)
            * C(ell - jj, s - jj + kk)
        )

    builders = [
        ("A", A_builder),
        ("B", B_builder),
        ("C", C_builder),
        ("D", D_builder),
    ]

    for label, builder in builders:
        for dk in (-1, 0, 1):
            for dj in (-1, 0, 1):
                cand = make_shifted_candidate(
                    label,
                    f"{label} with dk={dk}, dj={dj}",
                    builder,
                    dk,
                    dj,
                )

                matches = sum(
                    cand.fn(a.k, a.ell, a.s) == a.exact
                    for a in ANCHORS
                )

                if matches > 0:
                    print(
                        f"{cand.name:32s} "
                        f"anchors={matches}/{len(ANCHORS)}"
                    )

                    shifted.append(cand)

    return shifted


# ============================================================================
# 14. DIRECT GENERATING-FUNCTION CONTROL
# ============================================================================

def generating_control(k: int, ell: int, s: int) -> int:
    """
    Control family:
        coefficient extraction from (1+x)^k.

    Its degree is k, so it provides a useful sign/magnitude control, but is
    NOT assumed to be the answer.
    """
    r = s
    return C(k, r)


# ============================================================================
# 15. REPORT CONTROLS
# ============================================================================

def print_reference_sanity() -> None:
    print()
    print("=" * 78)
    print("REFERENCE CONVOLUTION SANITY CHECK")
    print("=" * 78)

    tests = [
        (1, 5, 2),
        (1, 6, 3),
        (3, 8, 3),
        (3, 9, 4),
        (5, 10, 5),
    ]

    failures = 0

    for k, ell, s in tests:
        value = exact_reference(k, ell, s)

        # Collapsed expression.
        collapsed = (
            C(ell, s)
            * sum((-1) ** j * C(s, j) for j in range(0, s + 1))
        )

        status = "PASS" if value == collapsed else "FAIL"
        if value != collapsed:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} s={s:2d} "
            f"sum={value:6d} collapsed={collapsed:6d} {status}"
        )

    print(f"\ncollapse failures = {failures}")


def print_anchor_baseline() -> None:
    print()
    print("=" * 78)
    print("REFERENCE ANCHOR BASELINE")
    print("=" * 78)

    failures = 0

    for a in ANCHORS:
        value = exact_reference(a.k, a.ell, a.s)

        status = "PASS" if value == a.exact else "FAIL"

        if value != a.exact:
            failures += 1

        print(
            f"{a.label}: "
            f"k={a.k:2d} ell={a.ell:2d} s={a.s:2d} "
            f"exact_anchor={a.exact:10d} "
            f"reference={value:10d} {status}"
        )

    print(f"\nreference anchor failures = {failures}/{len(ANCHORS)}")


# ============================================================================
# 16. MAIN
# ============================================================================

def main() -> None:
    print("=" * 78)
    print("EXPERIMENT 196")
    print("STAGED PRE-CONVOLUTION k-PLACEMENT SEARCH")
    print("=" * 78)

    print()
    print("This experiment is intentionally bounded.")
    print("No unrestricted affine candidate enumeration is performed.")
    print("Candidates are filtered before expensive support diagnostics.")

    print_reference_sanity()
    print_anchor_baseline()

    survivors = run_staged_search()

    # Only invoke the small shift search if the base search did not already
    # produce a convincing family.
    #
    # Even if it does produce one, the shift search remains diagnostic.
    if len(survivors) == 0:
        shifted = run_small_shift_search()
    else:
        shifted = []

    # ------------------------------------------------------------------
    # FINAL SUMMARY
    # ------------------------------------------------------------------
    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(f"base candidates tested      = {len(BASE_CANDIDATES)}")
    print(f"final base survivors        = {len(survivors)}")
    print(f"shift candidates reported   = {len(shifted)}")

    if survivors:
        print()
        print("SURVIVING BASE FAMILIES")
        for s in survivors:
            print(
                f"  {s.name:32s} "
                f"anchor={s.anchor_matches}/{s.anchor_total} "
                f"kdep={s.k_dependence_score} "
                f"error={s.neighborhood_error}"
            )
    else:
        print()
        print("No base family survived all staged filters.")

    if shifted:
        print()
        print("SHIFT CANDIDATES WITH NONZERO ANCHOR AGREEMENT")
        for cand in shifted:
            matches = sum(
                cand.fn(a.k, a.ell, a.s) == a.exact
                for a in ANCHORS
            )

            print(
                f"  {cand.name:32s} "
                f"{matches}/{len(ANCHORS)}"
            )

    print()
    print("-" * 78)
    print("INTERPRETATION")
    print("-" * 78)
    print(
        "A survivor is evidence only that the tested pre-convolution "
        "structure is compatible with the current diagnostics."
    )
    print(
        "It is NOT evidence that the candidate is the theoretical n=pq "
        "formula until the term is derived directly from the pq construction."
    )
    print()
    print(
        "The next mathematical step, after any survivors are identified, "
        "is to derive that summand from the exact n=pq map."
    )

    print()
    print("=" * 78)
    print("END EXPERIMENT 196")
    print("=" * 78)


if __name__ == "__main__":
    main()