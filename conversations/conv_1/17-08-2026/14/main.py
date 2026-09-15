#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 84
PRIME-POWER / FINITE-FIELD CONGRUENCE DIAGNOSTIC

FOLLOW-UP TO EXPERIMENT 83R2

FOCUS:
    The E8 mod 67 sequence showed only 12 mismatches at shift 67^2.

TESTS:
    E6(n + T) == E6(n) mod ell
    E8(n + T) == E8(n) mod ell

for:
    T = ell
    T = 2*ell
    T = ell^2
    T = 2*ell^2
    T = ell^3
    T = ell*(ell-1)
    T = ell*(ell+1)

Also:
    - mismatch counts
    - mismatch positions
    - mismatch residue classes mod ell
    - difference-value distributions
    - warmup sensitivity
    - direct comparison of ell=67 against all other moduli
    - optional gcd / offset diagnostics

NO FACTOR INPUT
NO CSV OUTPUT
NO SKLEARN
==============================================================================

The goal is NOT to fit another predictor.

The goal is to determine whether the near-67^2 phenomenon is:

    1. a genuine prime-power congruence,
    2. a finite-window coincidence,
    3. a boundary/warmup artifact,
    4. or a phenomenon occurring across several ell.

==============================================================================
"""

from __future__ import annotations

import math
import statistics
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


# ============================================================================
# CONFIGURATION
# ============================================================================

LIMIT = 12000

# Ignore an initial prefix when testing congruences.
WARMUPS = [0, 100, 500, 1000, 2000]

MODULI = [
    7,
    13,
    19,
    31,
    37,
    61,
    67,
]

# The previous anomaly.
FOCUS_ELL = 67

# Print at most this many mismatch positions per test.
MAX_PRINT_MISMATCHES = 30

# Also inspect mismatch neighborhoods around the first few failures.
NEIGHBORHOOD = 3


# ============================================================================
# MACMAHON q-SERIES
# ============================================================================

def compute_M123(
    limit: int,
    mod: int,
) -> Tuple[List[int], List[int], List[int]]:
    """
    Compute M1, M2, M3 modulo mod.

    U_a(q) =
        sum_{s1<...<sa}
            q^(s1+...+sa)
            / prod_i (1-q^si)^2

    and

        q^s/(1-q^s)^2 = sum_{m>=1} m q^(ms).

    dp[j][n] = coefficient using exactly j distinct part sizes.
    """

    dp = [
        [0] * (limit + 1)
        for _ in range(4)
    ]

    dp[0][0] = 1

    for s in range(1, limit + 1):

        max_j = min(
            3,
            limit // s,
        )

        for j in range(max_j, 0, -1):

            old = dp[j - 1]
            target = dp[j]

            for r in range(s):

                if r > limit:
                    break

                A = 0
                B = 0
                n = r

                while n + s <= limit:

                    A = (
                        A
                        + old[n]
                    ) % mod

                    B = (
                        B
                        + A
                    ) % mod

                    n += s

                    target[n] = (
                        target[n]
                        + B
                    ) % mod

    return dp[1], dp[2], dp[3]


# ============================================================================
# PAPER-DERIVED COEFFICIENTS
# ============================================================================

def E6_sequence(
    M1: Sequence[int],
    M2: Sequence[int],
    mod: int,
) -> List[int]:

    out = [0] * len(M1)

    for n in range(1, len(M1)):

        r = n % mod

        c1 = (
            r * r
            - 3 * r
            + 2
        ) % mod

        out[n] = (
            c1 * M1[n]
            - 8 * M2[n]
        ) % mod

    return out


def E8_sequence(
    M1: Sequence[int],
    M2: Sequence[int],
    M3: Sequence[int],
    mod: int,
) -> List[int]:

    out = [0] * len(M1)

    for n in range(1, len(M1)):

        r = n % mod

        c1 = (
            3 * r**3
            - 13 * r**2
            + 18 * r
            - 8
        ) % mod

        c2 = (
            12 * r**2
            - 120 * r
            + 212
        ) % mod

        out[n] = (
            c1 * M1[n]
            + c2 * M2[n]
            - 960 * M3[n]
        ) % mod

    return out


# ============================================================================
# CONGRUENCE TESTS
# ============================================================================

@dataclass
class CongruenceResult:
    ell: int
    shift: int
    warmup: int
    comparisons: int
    mismatches: int
    mismatch_rate: float
    match_rate: float
    difference_counts: Dict[int, int]
    mismatch_positions: List[int]
    mismatch_residues: Counter


def compare_shift(
    seq: Sequence[int],
    ell: int,
    shift: int,
    warmup: int,
) -> CongruenceResult:

    start = max(
        1,
        warmup,
    )

    end = len(seq) - 1 - shift

    if end < start:
        return CongruenceResult(
            ell=ell,
            shift=shift,
            warmup=warmup,
            comparisons=0,
            mismatches=0,
            mismatch_rate=1.0,
            match_rate=0.0,
            difference_counts={},
            mismatch_positions=[],
            mismatch_residues=Counter(),
        )

    mismatches = 0
    differences = Counter()
    mismatch_positions: List[int] = []
    mismatch_residues = Counter()

    for n in range(start, end + 1):

        a = seq[n]
        b = seq[n + shift]

        if a != b:
            mismatches += 1

            d = (
                b - a
            ) % ell

            differences[d] += 1

            mismatch_residues[
                n % ell
            ] += 1

            if len(mismatch_positions) < MAX_PRINT_MISMATCHES:
                mismatch_positions.append(n)

    comparisons = (
        end - start + 1
    )

    return CongruenceResult(
        ell=ell,
        shift=shift,
        warmup=warmup,
        comparisons=comparisons,
        mismatches=mismatches,
        mismatch_rate=(
            mismatches / comparisons
            if comparisons
            else 1.0
        ),
        match_rate=(
            1.0 - mismatches / comparisons
            if comparisons
            else 0.0
        ),
        difference_counts=dict(differences),
        mismatch_positions=mismatch_positions,
        mismatch_residues=mismatch_residues,
    )


# ============================================================================
# SHIFT FAMILY
# ============================================================================

def candidate_shifts(
    ell: int,
) -> List[Tuple[str, int]]:

    candidates = [
        ("ell", ell),
        ("2ell", 2 * ell),
        ("ell^2", ell**2),
        ("2ell^2", 2 * ell**2),
        ("ell^3", ell**3),
        ("ell(ell-1)", ell * (ell - 1)),
        ("ell(ell+1)", ell * (ell + 1)),
    ]

    # Remove duplicate shifts while retaining labels.
    seen = set()
    out = []

    for label, shift in candidates:

        if shift <= 0:
            continue

        if shift in seen:
            continue

        seen.add(shift)
        out.append((label, shift))

    return out


# ============================================================================
# MISMATCH RESIDUE SUMMARY
# ============================================================================

def summarize_residue_concentration(
    result: CongruenceResult,
) -> str:

    if not result.mismatch_residues:
        return "none"

    top = (
        result.mismatch_residues
        .most_common(10)
    )

    return str(top)


# ============================================================================
# NEIGHBORHOOD ANALYSIS
# ============================================================================

def print_mismatch_neighborhoods(
    seq: Sequence[int],
    result: CongruenceResult,
) -> None:

    if not result.mismatch_positions:
        print("  no mismatches")

        return

    print("  mismatch neighborhoods:")

    for pos in result.mismatch_positions[:10]:

        left = max(
            0,
            pos - NEIGHBORHOOD,
        )

        right = min(
            len(seq) - result.shift - 1,
            pos + NEIGHBORHOOD,
        )

        values = []

        for n in range(
            left,
            right + 1,
        ):

            a = seq[n]
            b = seq[n + result.shift]

            values.append(
                (
                    n,
                    a,
                    b,
                    (b - a) % result.ell,
                )
            )

        print(
            f"    n={pos}: "
            f"{values}"
        )


# ============================================================================
# Q-SERIES VALIDATION
# ============================================================================

def direct_Ma(
    n: int,
    depth: int,
    mod: int,
) -> int:

    total = 0

    def rec(
        next_part: int,
        remaining: int,
        chosen: int,
        multiplicity_product: int,
    ) -> None:

        nonlocal total

        if chosen == depth:
            if remaining == 0:
                total = (
                    total
                    + multiplicity_product
                ) % mod

            return

        if remaining <= 0:
            return

        for s in range(
            next_part,
            remaining + 1,
        ):

            max_mult = (
                remaining // s
            )

            for m in range(
                1,
                max_mult + 1,
            ):

                rec(
                    s + 1,
                    remaining - m * s,
                    chosen + 1,
                    (
                        multiplicity_product
                        * m
                    ) % mod,
                )

    rec(
        1,
        n,
        0,
        1,
    )

    return total % mod


def validate_small(
    M1: Sequence[int],
    M2: Sequence[int],
    M3: Sequence[int],
    mod: int,
    check_limit: int = 25,
) -> bool:

    for n in range(
        1,
        check_limit + 1,
    ):

        if (
            M1[n]
            != direct_Ma(n, 1, mod)
        ):
            return False

        if (
            M2[n]
            != direct_Ma(n, 2, mod)
        ):
            return False

        if (
            M3[n]
            != direct_Ma(n, 3, mod)
        ):
            return False

    return True


# ============================================================================
# SINGLE MODULUS REPORT
# ============================================================================

def analyze_modulus(
    ell: int,
) -> Tuple[List[int], List[CongruenceResult]]:

    print("\n")
    print("=" * 78)
    print(f"MODULUS ell = {ell}")
    print("=" * 78)

    t0 = time.perf_counter()

    M1, M2, M3 = compute_M123(
        LIMIT,
        ell,
    )

    elapsed = (
        time.perf_counter()
        - t0
    )

    print(
        f"M1/M2/M3 computation = "
        f"{elapsed:.6f}s"
    )

    valid = validate_small(
        M1,
        M2,
        M3,
        ell,
    )

    print(
        f"independent small-n validation = "
        f"{valid}"
    )

    if not valid:
        raise RuntimeError(
            f"q-series validation failed for ell={ell}"
        )

    E6 = E6_sequence(
        M1,
        M2,
        ell,
    )

    E8 = E8_sequence(
        M1,
        M2,
        M3,
        ell,
    )

    sequences = [
        ("E6", E6),
        ("E8", E8),
    ]

    all_results: List[CongruenceResult] = []

    for seq_name, seq in sequences:

        print("\n")
        print(
            f"{seq_name} PRIME-POWER / SHIFT TEST"
        )
        print("-" * 78)

        for shift_name, shift in candidate_shifts(ell):

            # Use warmup=500 for the main result.
            result = compare_shift(
                seq,
                ell,
                shift,
                500,
            )

            all_results.append(result)

            print(
                f"{shift_name:12s} "
                f"T={shift:8d} "
                f"match={result.match_rate:.8f} "
                f"mismatch="
                f"{result.mismatches:6d}/"
                f"{result.comparisons:<6d}"
            )

            if (
                ell == FOCUS_ELL
                and shift_name == "ell^2"
                and seq_name == "E8"
            ):

                print(
                    "\nFOCUS DIAGNOSTIC: "
                    "E8 mod 67 at shift 67^2"
                )

                print(
                    f"mismatch rate = "
                    f"{result.mismatch_rate:.8f}"
                )

                print(
                    f"difference distribution = "
                    f"{result.difference_counts}"
                )

                print(
                    "mismatch residue classes mod 67:"
                )

                print(
                    summarize_residue_concentration(
                        result
                    )
                )

                print_mismatch_neighborhoods(
                    seq,
                    result,
                )

    return E8, all_results


# ============================================================================
# WARMUP SENSITIVITY
# ============================================================================

def warmup_sensitivity(
    seq: Sequence[int],
    ell: int,
    shift: int,
) -> None:

    print("\nWARMUP SENSITIVITY")
    print("-" * 78)

    for warmup in WARMUPS:

        result = compare_shift(
            seq,
            ell,
            shift,
            warmup,
        )

        print(
            f"warmup={warmup:5d} "
            f"match={result.match_rate:.8f} "
            f"mismatches="
            f"{result.mismatches}/"
            f"{result.comparisons}"
        )


# ============================================================================
# CROSS-MODULUS PRIME-SQUARE TABLE
# ============================================================================

def print_cross_modulus_square_table(
    stored_E8: Dict[int, List[int]],
) -> None:

    print("\n")
    print("=" * 78)
    print("CROSS-MODULUS E8 PRIME-SQUARE TABLE")
    print("=" * 78)

    print(
        "ell | ell^2 | match rate | mismatch rate"
    )

    for ell in MODULI:

        seq = stored_E8[ell]

        result = compare_shift(
            seq,
            ell,
            ell * ell,
            500,
        )

        print(
            f"{ell:3d} | "
            f"{ell * ell:5d} | "
            f"{result.match_rate:10.8f} | "
            f"{result.mismatch_rate:13.8f}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    overall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 84")
    print("PRIME-POWER / FINITE-FIELD CONGRUENCE DIAGNOSTIC")
    print("FOLLOW-UP TO EXPERIMENT 83R2")
    print("FOCUS: THE 67^2 E8 ANOMALY")
    print("NO FACTOR INPUT")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    print("\n1. PARAMETERS")
    print("-" * 78)
    print(f"q-series limit = {LIMIT}")
    print(f"moduli         = {MODULI}")
    print(f"warmups        = {WARMUPS}")
    print(f"focus ell      = {FOCUS_ELL}")
    print(
        "candidate shifts = ell, 2ell, ell^2, 2ell^2, "
        "ell^3, ell(ell-1), ell(ell+1)"
    )

    stored_E8: Dict[int, List[int]] = {}
    all_results: Dict[int, List[CongruenceResult]] = {}

    # ----------------------------------------------------------------------
    # Analyze every modulus.
    # ----------------------------------------------------------------------

    for ell in MODULI:

        E8, results = analyze_modulus(
            ell,
        )

        stored_E8[ell] = E8
        all_results[ell] = results

    # ----------------------------------------------------------------------
    # Warmup sensitivity for the key anomaly.
    # ----------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("2. FOCUS: E8 mod 67, SHIFT = 67^2")
    print("=" * 78)

    focus_seq = stored_E8[FOCUS_ELL]

    warmup_sensitivity(
        focus_seq,
        FOCUS_ELL,
        FOCUS_ELL**2,
    )

    # ----------------------------------------------------------------------
    # Cross-modulus square comparison.
    # ----------------------------------------------------------------------

    print_cross_modulus_square_table(
        stored_E8,
    )

    # ----------------------------------------------------------------------
    # Search for the best prime-power shift for each modulus.
    # ----------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("3. BEST SHIFT PER MODULUS")
    print("=" * 78)

    print(
        "ell | best shift | T | match rate | mismatches"
    )

    for ell in MODULI:

        seq = stored_E8[ell]

        best = None

        for label, shift in candidate_shifts(ell):

            result = compare_shift(
                seq,
                ell,
                shift,
                500,
            )

            key = (
                result.match_rate,
                -result.mismatches,
            )

            if best is None or key > best[0]:
                best = (
                    key,
                    label,
                    shift,
                    result,
                )

        _, label, shift, result = best

        print(
            f"{ell:3d} | "
            f"{label:12s} | "
            f"{shift:7d} | "
            f"{result.match_rate:.8f} | "
            f"{result.mismatches}"
        )

    # ----------------------------------------------------------------------
    # Final focus summary.
    # ----------------------------------------------------------------------

    focus_result = compare_shift(
        focus_seq,
        FOCUS_ELL,
        FOCUS_ELL**2,
        500,
    )

    print("\n")
    print("=" * 78)
    print("4. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        f"E8 mod 67, shift 67^2:"
    )

    print(
        f"  comparisons = "
        f"{focus_result.comparisons}"
    )

    print(
        f"  mismatches   = "
        f"{focus_result.mismatches}"
    )

    print(
        f"  match rate   = "
        f"{focus_result.match_rate:.10f}"
    )

    print(
        f"  mismatch rate = "
        f"{focus_result.mismatch_rate:.10f}"
    )

    print(
        f"  difference distribution = "
        f"{focus_result.difference_counts}"
    )

    print(
        f"  mismatch residue classes = "
        f"{focus_result.mismatch_residues.most_common()}"
    )

    print()
    print("INTERPRETATION RULES")
    print("--------------------")

    if (
        focus_result.mismatch_rate < 0.01
    ):
        print(
            "VERY STRONG CONGRUENCE CANDIDATE:"
        )
        print(
            "The 67^2 shift is unusually close to an exact "
            "E8 congruence and deserves a dedicated proof/search."
        )

    elif (
        focus_result.mismatch_rate < 0.05
    ):
        print(
            "INTERESTING NEAR-CONGRUENCE:"
        )
        print(
            "The 67^2 phenomenon is unusually strong but not exact."
        )

    else:
        print(
            "NO SPECIAL 67^2 EFFECT:"
        )
        print(
            "The previous anomaly does not survive this wider test."
        )

    print()
    print(
        "A genuine prime-power congruence should ideally:"
    )

    print(
        "  1. survive larger warmup windows,"
    )
    print(
        "  2. remain strong on a longer q-series window,"
    )
    print(
        "  3. show a structured mismatch pattern rather than random noise,"
    )
    print(
        "  4. have a corresponding explanation in the finite-field/"
        "quasimodular structure,"
    )
    print(
        "  5. preferably recur for more than one modulus."
    )

    runtime = (
        time.perf_counter()
        - overall_start
    )

    print(
        f"\ntotal runtime = {runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 84 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

