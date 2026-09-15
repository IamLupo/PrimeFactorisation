#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 195
PQ -> (k, ell, s) REDUCTION AUDIT
==============================================================================
Purpose
-------
The previous experiments indicate:

    exact convolution
        -> k-independent cancellation
        -> contradiction with known k-dependent anchors.

Experiment 195 therefore does NOT search for another arbitrary closed form.

Instead it audits the reduction itself.

Main questions
--------------
1. At what stage can k enter before the binomial cancellation?
2. Does the current exact convolution erase that k-dependence?
3. Can the n=pq labels be carried through the reduction without losing k?
4. Which structural placement of k survives the audit?

Important
---------
The exact map

    (p, q) -> (k, ell, s)

must come from the original n=pq derivation.

This script therefore exposes `pq_to_reduced_state()` as an explicit hook.
Do not silently guess this mapping.
==============================================================================

"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
from itertools import combinations
from typing import Callable, Iterable, Optional


# ============================================================================
# 1. DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class PQCase:
    p: int
    q: int
    n: int


@dataclass(frozen=True)
class ReducedState:
    k: int
    ell: int
    s: int


@dataclass
class Comparison:
    name: str
    value: int


# ============================================================================
# 2. PRIME GENERATION
# ============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def odd_primes(limit: int) -> list[int]:
    return [n for n in range(3, limit + 1, 2) if is_prime(n)]


def generate_pq_cases(limit: int = 31) -> list[PQCase]:
    ps = odd_primes(limit)
    out: list[PQCase] = []

    for p, q in combinations(ps, 2):
        out.append(PQCase(p=p, q=q, n=p * q))

    return out


# ============================================================================
# 3. PQ -> REDUCED-STATE MAPPING
# ============================================================================

def pq_to_reduced_state(p: int, q: int) -> Optional[ReducedState]:
    """
    --------------------------------------------------------------------------
    REQUIRED ORIGINAL DERIVATION HOOK
    --------------------------------------------------------------------------

    Replace this function with the exact map from the n = p*q derivation.

    The experiment intentionally refuses to fabricate that map.

    Return:
        ReducedState(k=k, ell=ell, s=s)

    or:
        None

    for a pq-case that is outside the relevant branch/support.

    --------------------------------------------------------------------------
    TEMPORARY SAFETY DEFAULT
    --------------------------------------------------------------------------
    We use a simple diagnostic embedding only so the script can execute:

        k   = p - 2
        ell = q + p
        s   = q // 2

    THIS IS NOT CLAIMED TO BE THE THEORETICAL MAP.
    It only exercises the audit machinery.

    Replace immediately with the actual formula.
    """

    # ------------------------------------------------------------
    # TEMPORARY PLACEHOLDER
    # ------------------------------------------------------------
    k = p - 2
    ell = p + q
    s = q // 2

    if k <= 0 or ell < 0 or s < 0:
        return None

    return ReducedState(k=k, ell=ell, s=s)


# ============================================================================
# 4. CURRENT REFERENCE EXACT CONVOLUTION
# ============================================================================

def exact_reference(k: int, ell: int, s: int) -> int:
    """
    This is the convolution diagnosed in Experiments 185-186:

        sum_j (-1)^j C(ell,j) C(ell-j, s-j)

    The k parameter is present only as an argument to make k-independence
    explicit; it does not occur in the summand.
    """

    total = 0

    lo = max(0, s - ell)
    hi = min(ell, s)

    for j in range(lo, hi + 1):
        term = (
            (-1) ** j
            * comb(ell, j)
            * comb(ell - j, s - j)
        )
        total += term

    return total


def exact_collapsed(ell: int, s: int) -> int:
    """
    By

        C(ell,j) C(ell-j,s-j)
        = C(ell,s) C(s,j),

    this becomes

        C(ell,s) * sum_j (-1)^j C(s,j)

    and hence:
        s = 0 -> 1
        s > 0 -> 0
    """

    if s < 0 or s > ell:
        return 0

    if s == 0:
        return 1

    return 0


# ============================================================================
# 5. K-DEPENDENCE DIAGNOSTICS
# ============================================================================

def current_k_sweep(
    ell: int,
    s: int,
    k_values: Iterable[int],
) -> list[Comparison]:

    out = []

    for k in k_values:
        out.append(
            Comparison(
                name=f"k={k}",
                value=exact_reference(k, ell, s),
            )
        )

    return out


def values_equal(comparisons: list[Comparison]) -> bool:
    if not comparisons:
        return True

    first = comparisons[0].value
    return all(c.value == first for c in comparisons)


def print_k_sweep(ell: int, s: int, k_values: Iterable[int]) -> None:
    values = current_k_sweep(ell, s, k_values)

    print(f"ell={ell:2d} s={s:2d}")

    for item in values:
        print(f"  {item.name:>5s} -> {item.value}")

    print(
        "  k-independent = "
        + ("YES" if values_equal(values) else "NO")
    )


# ============================================================================
# 6. PRE-COLLAPSE K-DEPENDENT CANDIDATES
# ============================================================================
#
# These are not claimed to be the answer.
#
# They are the structural placements already identified in Experiment 188.
#
# A = second upper + k
# B = target index - k
# C = truncated j <= k
# D = second upper + k, target - k
# E = first upper - k
# F = first upper + k
#
# We test them here only as "where can k survive?" diagnostics.
# ============================================================================

def candidate_A(k: int, ell: int, s: int) -> int:
    """
    Second upper binomial shifted by +k:

        sum_j (-1)^j C(ell,j) C(ell-j+k, s-j)

    """
    total = 0

    for j in range(0, min(ell, s) + 1):
        a = ell - j + k
        b = s - j

        if a < 0 or b < 0 or b > a:
            continue

        total += (-1) ** j * comb(ell, j) * comb(a, b)

    return total


def candidate_B(k: int, ell: int, s: int) -> int:
    """
    Target index shifted by -k.
    """

    target = s - k
    if target < 0:
        return 0

    total = 0

    for j in range(0, ell + 1):
        if j > target:
            break

        a = ell - j
        b = target - j

        if 0 <= b <= a:
            total += (-1) ** j * comb(ell, j) * comb(a, b)

    return total


def candidate_C(k: int, ell: int, s: int) -> int:
    """
    Truncate the j-range at j <= k.
    """

    total = 0

    hi = min(ell, s, k)

    for j in range(0, hi + 1):
        total += (
            (-1) ** j
            * comb(ell, j)
            * comb(ell - j, s - j)
        )

    return total


def candidate_D(k: int, ell: int, s: int) -> int:
    """
    Combine A and B.
    """

    target = s - k
    if target < 0:
        return 0

    total = 0

    for j in range(0, min(ell, target) + 1):
        a = ell - j + k
        b = target - j

        if a < 0 or b < 0 or b > a:
            continue

        total += (-1) ** j * comb(ell, j) * comb(a, b)

    return total


def candidate_E(k: int, ell: int, s: int) -> int:
    """
    First upper binomial shifted by -k.
    """

    total = 0

    for j in range(0, min(ell, s) + 1):
        a = ell - k
        if a < j:
            continue

        b = j
        if b > a:
            continue

        rhs_top = a - j
        rhs_bottom = s - j

        if rhs_bottom < 0 or rhs_bottom > rhs_top:
            continue

        total += (
            (-1) ** j
            * comb(a, j)
            * comb(rhs_top, rhs_bottom)
        )

    return total


def candidate_F(k: int, ell: int, s: int) -> int:
    """
    First upper + k.
    """

    total = 0

    a0 = ell + k

    for j in range(0, min(a0, s) + 1):
        if j > a0:
            break

        b = s - j
        top2 = a0 - j

        if 0 <= b <= top2:
            total += (
                (-1) ** j
                * comb(a0, j)
                * comb(top2, b)
            )

    return total


CANDIDATES: dict[str, Callable[[int, int, int], int]] = {
    "A_second_upper_plus_k": candidate_A,
    "B_target_minus_k": candidate_B,
    "C_truncate_j_at_k": candidate_C,
    "D_A_plus_B": candidate_D,
    "E_first_upper_minus_k": candidate_E,
    "F_first_upper_plus_k": candidate_F,
}


# ============================================================================
# 7. KNOWN ANCHORS FROM EXPERIMENT 192
# ============================================================================

ANCHORS = [
    ("A1", 1, 10, 2, 27),
    ("A2", 3, 20, 6, 935),
    ("A3", 7, 40, 15, -1797818),
    ("B1", 3, 7, 3, -3),
    ("B2", 3, 9, 4, 9),
    ("B3", 3, 11, 5, -16),
    ("B4", 5, 11, 5, -5),
    ("B5", 5, 19, 9, -196),
    ("B6", 9, 21, 11, -84),
]


# ============================================================================
# 8. ANCHOR AUDIT
# ============================================================================

def anchor_audit() -> None:
    print("\n" + "=" * 78)
    print("ANCHOR AUDIT")
    print("=" * 78)

    current_failures = 0

    for label, k, ell, s, exact in ANCHORS:
        got = exact_reference(k, ell, s)

        ok = (got == exact)

        if not ok:
            current_failures += 1

        print(
            f"{label:>2s}: "
            f"k={k:2d} ell={ell:2d} s={s:2d} "
            f"exact_anchor={exact:>9d} "
            f"current={got:>9d} "
            f"{'PASS' if ok else 'FAIL'}"
        )

    print(f"\nCurrent reference anchor failures = {current_failures}/{len(ANCHORS)}")


# ============================================================================
# 9. STRUCTURAL CANDIDATE AUDIT
# ============================================================================

def candidate_audit() -> None:
    print("\n" + "=" * 78)
    print("PRE-COLLAPSE K-DEPENDENCE AUDIT")
    print("=" * 78)

    for name, fn in CANDIDATES.items():
        matches = 0

        for _, k, ell, s, exact in ANCHORS:
            got = fn(k, ell, s)

            if got == exact:
                matches += 1

        print(f"{name:28s}: {matches}/{len(ANCHORS)} anchor matches")


# ============================================================================
# 10. PQ PROPAGATION AUDIT
# ============================================================================

def pq_propagation_audit(
    pq_cases: list[PQCase],
    state_map: Callable[[int, int], Optional[ReducedState]],
) -> None:

    print("\n" + "=" * 78)
    print("PQ -> REDUCED STATE -> CONVOLUTION AUDIT")
    print("=" * 78)

    mapped = 0
    undefined = 0
    current_zero = 0
    current_nonzero = 0

    for case in pq_cases:
        state = state_map(case.p, case.q)

        if state is None:
            undefined += 1
            continue

        mapped += 1

        value = exact_reference(
            state.k,
            state.ell,
            state.s,
        )

        if value == 0:
            current_zero += 1
        else:
            current_nonzero += 1

        print(
            f"p={case.p:2d} q={case.q:2d} n={case.n:4d}"
            f" -> k={state.k:2d} ell={state.ell:3d} s={state.s:2d}"
            f" -> current={value:>8d}"
        )

    print("\nPQ summary")
    print(f"  total cases     = {len(pq_cases)}")
    print(f"  mapped          = {mapped}")
    print(f"  undefined       = {undefined}")
    print(f"  current zero    = {current_zero}")
    print(f"  current nonzero = {current_nonzero}")


# ============================================================================
# 11. FIXED-(ell,s), VARY-k AUDIT
# ============================================================================

def fixed_ell_s_k_audit() -> None:
    print("\n" + "=" * 78)
    print("FIXED-(ell,s), VARY-k AUDIT")
    print("=" * 78)

    test_pairs = [
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

    for ell, s in test_pairs:
        max_k = min(11, ell)

        print_k_sweep(
            ell=ell,
            s=s,
            k_values=range(1, max_k + 1, 2),
        )


# ============================================================================
# 12. SIGN / ODD-BRANCH AUDIT
# ============================================================================

def odd_branch_sign_audit() -> None:
    print("\n" + "=" * 78)
    print("ODD-BRANCH SIGN AUDIT")
    print("=" * 78)

    examples = [
        (1, 10, 2),
        (3, 20, 6),
        (7, 40, 15),
        (3, 7, 3),
        (3, 9, 4),
        (5, 11, 5),
        (9, 21, 11),
    ]

    for k, ell, s in examples:
        ref = exact_reference(k, ell, s)

        signs = {
            "none": 1,
            "(-1)^s": (-1) ** s,
            "(-1)^k": (-1) ** k,
            "(-1)^ell": (-1) ** ell,
            "(-1)^(k+s)": (-1) ** (k + s),
            "(-1)^((ell-k)/2+s)": (
                -1
                if ((ell - k) // 2 + s) % 2
                else 1
            ),
        }

        print(
            f"k={k:2d} ell={ell:2d} s={s:2d}"
            f" current={ref:>8d}"
        )

        for name, sign in signs.items():
            print(f"    {name:25s} -> {ref * sign:>8d}")


# ============================================================================
# 13. MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("EXPERIMENT 195")
    print("PQ -> (k, ell, s) REDUCTION AUDIT")
    print("=" * 78)

    print(
        """
This experiment is a reduction audit.

It is intentionally NOT claiming that the temporary
pq_to_reduced_state() map is the theoretical n=pq map.

Replace that function with the exact derivation before using the
pq-level results as mathematical evidence.
"""
    )

    # ------------------------------------------------------------
    # A. Reference convolution audit
    # ------------------------------------------------------------
    print("\n" + "=" * 78)
    print("1. REFERENCE CONVOLUTION SANITY CHECK")
    print("=" * 78)

    sanity_cases = [
        (1, 5, 2),
        (1, 6, 3),
        (3, 8, 3),
        (3, 9, 4),
        (5, 10, 5),
    ]

    sanity_failures = 0

    for k, ell, s in sanity_cases:
        a = exact_reference(k, ell, s)
        b = exact_collapsed(ell, s)

        ok = (a == b)

        if not ok:
            sanity_failures += 1

        print(
            f"k={k:2d} ell={ell:2d} s={s:2d}"
            f" sum={a:>6d}"
            f" collapsed={b:>6d}"
            f" {'PASS' if ok else 'FAIL'}"
        )

    print(f"\ncollapse failures = {sanity_failures}")

    # ------------------------------------------------------------
    # B. Existing anchor data
    # ------------------------------------------------------------
    anchor_audit()

    # ------------------------------------------------------------
    # C. Structural candidate audit
    # ------------------------------------------------------------
    candidate_audit()

    # ------------------------------------------------------------
    # D. Explicit k-independence audit
    # ------------------------------------------------------------
    fixed_ell_s_k_audit()

    # ------------------------------------------------------------
    # E. Odd-branch sign audit
    # ------------------------------------------------------------
    odd_branch_sign_audit()

    # ------------------------------------------------------------
    # F. pq propagation
    # ------------------------------------------------------------
    pq_cases = generate_pq_cases(limit=31)

    pq_propagation_audit(
        pq_cases=pq_cases,
        state_map=pq_to_reduced_state,
    )

    # ------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # ------------------------------------------------------------
    print("\n" + "=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The intended interpretation is:

A. If the exact convolution is k-independent for fixed (ell,s),
   then the k-dependence has already been lost before the final
   closed-form stage.

B. If the pq -> (k,ell,s) map varies k but the convolution does not,
   then specializing to n=pq cannot restore the lost dependence.

C. Therefore the next derivation target is the PRE-CONVOLUTION
   summand, truncation, or index transformation.

D. Candidate families A-F are diagnostic only. A surviving candidate
   must be traced back to an explicit pq derivation.

NEXT REQUIRED STEP:
Replace pq_to_reduced_state() with the exact n=pq map and rerun.
"""
    )


if __name__ == "__main__":
    main()

