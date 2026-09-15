#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 83R2
QUASIMODULAR FINITE-STATE / LINEAR-RECURRENCE SEARCH

FIXED PERIOD CONFIGURATION
HIGH-ORDER RECURRENCE DIAGNOSTIC
HELD-OUT VALIDATION
NO INVALID LARGE-N EXTRAPOLATION

NO FACTOR INPUT
NO CSV OUTPUT
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import math
import time
from typing import Dict, List, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

LIMIT = 5000
TRAIN_END = 3500
VALIDATION_END = LIMIT

MODULI = [7, 13, 19, 31, 37, 61, 67]

# Restored configuration variable from the previous script.
PERIOD_WARMUP = 500

# We regard recurrences above this order as "high complexity".
MAX_USEFUL_ORDER = 50

LARGE_N_TESTS = [
    6_536_304,
    6_067_312,
    6_647_970,
    7_548_800,
    5_111_112,
    6_584_270,
    7_021_700,
    6_494_300,
    8_107_520,
    6_733_028,
]


# ============================================================================
# MACMAHON q-SERIES
# ============================================================================

def compute_M123(
    limit: int,
    mod: int,
) -> Tuple[List[int], List[int], List[int]]:
    """
    Compute M1, M2, M3 modulo mod.

    dp[j][n] is the coefficient using exactly j distinct part sizes.
    """

    dp = [
        [0] * (limit + 1)
        for _ in range(4)
    ]

    dp[0][0] = 1

    for s in range(1, limit + 1):

        max_j = min(3, limit // s)

        for j in range(max_j, 0, -1):

            old = dp[j - 1]
            target = dp[j]

            for r in range(s):

                prev_a = 0
                prev_b = 0
                n = r

                while n + s <= limit:

                    prev_a = (
                        prev_a
                        + old[n]
                    ) % mod

                    prev_b = (
                        prev_b
                        + prev_a
                    ) % mod

                    n += s

                    target[n] = (
                        target[n]
                        + prev_b
                    ) % mod

    return dp[1], dp[2], dp[3]


# ============================================================================
# INDEPENDENT SMALL-N VALIDATION
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

            max_mult = remaining // s

            for m in range(1, max_mult + 1):

                rec(
                    s + 1,
                    remaining - m * s,
                    chosen + 1,
                    (
                        multiplicity_product * m
                    ) % mod,
                )

    rec(
        1,
        n,
        0,
        1,
    )

    return total % mod


def validate_q_series(
    M1: Sequence[int],
    M2: Sequence[int],
    M3: Sequence[int],
    mod: int,
    limit: int = 30,
) -> bool:

    for n in range(1, limit + 1):

        d1 = direct_Ma(n, 1, mod)
        d2 = direct_Ma(n, 2, mod)
        d3 = direct_Ma(n, 3, mod)

        if M1[n] != d1:
            print(
                f"FAIL M1 n={n}: "
                f"computed={M1[n]} direct={d1}"
            )
            return False

        if M2[n] != d2:
            print(
                f"FAIL M2 n={n}: "
                f"computed={M2[n]} direct={d2}"
            )
            return False

        if M3[n] != d3:
            print(
                f"FAIL M3 n={n}: "
                f"computed={M3[n]} direct={d3}"
            )
            return False

    return True


# ============================================================================
# PRIME-DETECTING SEQUENCES
# ============================================================================

def E6_sequence(
    M1: Sequence[int],
    M2: Sequence[int],
    mod: int,
) -> List[int]:

    seq = [0] * len(M1)

    for n in range(1, len(M1)):

        nn = n % mod

        c1 = (
            nn * nn
            - 3 * nn
            + 2
        ) % mod

        seq[n] = (
            c1 * M1[n]
            - 8 * M2[n]
        ) % mod

    return seq


def E8_sequence(
    M1: Sequence[int],
    M2: Sequence[int],
    M3: Sequence[int],
    mod: int,
) -> List[int]:

    seq = [0] * len(M1)

    for n in range(1, len(M1)):

        nn = n % mod

        c1 = (
            3 * nn**3
            - 13 * nn**2
            + 18 * nn
            - 8
        ) % mod

        c2 = (
            12 * nn**2
            - 120 * nn
            + 212
        ) % mod

        seq[n] = (
            c1 * M1[n]
            + c2 * M2[n]
            - 960 * M3[n]
        ) % mod

    return seq


# ============================================================================
# BERLEKAMP-MASSEY
# ============================================================================

def bm(
    sequence: Sequence[int],
    mod: int,
) -> List[int]:

    C = [1]
    B = [1]

    L = 0
    m = 1
    b = 1

    for n in range(len(sequence)):

        discrepancy = sequence[n] % mod

        for i in range(1, L + 1):
            discrepancy = (
                discrepancy
                + C[i] * sequence[n - i]
            ) % mod

        if discrepancy == 0:
            m += 1
            continue

        inv_b = pow(
            b,
            mod - 2,
            mod,
        )

        coef = (
            discrepancy
            * inv_b
        ) % mod

        T = C[:]

        needed = (
            len(B) + m
        )

        if len(C) < needed:
            C.extend(
                [0] * (
                    needed
                    - len(C)
                )
            )

        for j in range(len(B)):
            C[j + m] = (
                C[j + m]
                - coef * B[j]
            ) % mod

        if 2 * L <= n:
            L = n + 1 - L
            B = T
            b = discrepancy
            m = 1
        else:
            m += 1

    return C


def recurrence_order(C: Sequence[int]) -> int:
    return max(0, len(C) - 1)


def recurrence_errors(
    sequence: Sequence[int],
    C: Sequence[int],
    mod: int,
    start: int,
    end: int,
) -> int:

    L = recurrence_order(C)

    if L == 0:
        return 0

    errors = 0

    for n in range(
        max(start, L),
        min(end, len(sequence) - 1) + 1,
    ):

        predicted = 0

        for i in range(1, L + 1):
            predicted -= (
                C[i]
                * sequence[n - i]
            )

        predicted %= mod

        if predicted != sequence[n] % mod:
            errors += 1

    return errors


# ============================================================================
# PERIOD TEST
# ============================================================================

def period_errors(
    sequence: Sequence[int],
    period: int,
    start: int,
) -> int:

    if period <= 0:
        return 10**18

    errors = 0

    begin = max(
        start + period,
        period + 1,
    )

    for n in range(begin, len(sequence)):
        if sequence[n] != sequence[n - period]:
            errors += 1

    return errors


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    overall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 83R2")
    print("QUASIMODULAR FINITE-STATE / LINEAR-RECURRENCE SEARCH")
    print("FIXED PERIOD CONFIGURATION")
    print("HIGH-ORDER RECURRENCE DIAGNOSTIC")
    print("HELD-OUT VALIDATION")
    print("NO FACTOR INPUT")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    print("\n1. PARAMETERS")
    print("-" * 78)
    print(f"q-series limit       = {LIMIT}")
    print(f"training endpoint    = {TRAIN_END}")
    print(f"validation endpoint  = {VALIDATION_END}")
    print(f"period warmup        = {PERIOD_WARMUP}")
    print(f"moduli               = {MODULI}")
    print(f"max useful order     = {MAX_USEFUL_ORDER}")

    summaries = []

    for mod in MODULI:

        print("\n")
        print("=" * 78)
        print(f"MODULUS ell = {mod}")
        print("=" * 78)

        # --------------------------------------------------------------
        # q-series
        # --------------------------------------------------------------

        t0 = time.perf_counter()

        M1, M2, M3 = compute_M123(
            LIMIT,
            mod,
        )

        q_time = (
            time.perf_counter()
            - t0
        )

        print(
            f"M1/M2/M3 computation = "
            f"{q_time:.6f}s"
        )

        # --------------------------------------------------------------
        # Validate
        # --------------------------------------------------------------

        t0 = time.perf_counter()

        valid = validate_q_series(
            M1,
            M2,
            M3,
            mod,
        )

        validation_time = (
            time.perf_counter()
            - t0
        )

        print(
            f"independent validation = "
            f"{valid}"
        )

        print(
            f"validation time = "
            f"{validation_time:.6f}s"
        )

        if not valid:
            raise RuntimeError(
                f"q-series validation failed at ell={mod}"
            )

        # --------------------------------------------------------------
        # Sequences
        # --------------------------------------------------------------

        E6 = E6_sequence(
            M1,
            M2,
            mod,
        )

        E8 = E8_sequence(
            M1,
            M2,
            M3,
            mod,
        )

        # --------------------------------------------------------------
        # BM
        # --------------------------------------------------------------

        C6 = bm(
            E6[:TRAIN_END],
            mod,
        )

        C8 = bm(
            E8[:TRAIN_END],
            mod,
        )

        L6 = recurrence_order(C6)
        L8 = recurrence_order(C8)

        print("\n2. BERLEKAMP-MASSEY")
        print("-" * 78)
        print(f"E6 recurrence order = {L6}")
        print(f"E8 recurrence order = {L8}")

        # --------------------------------------------------------------
        # Held-out validation
        # --------------------------------------------------------------

        e6_errors = recurrence_errors(
            E6,
            C6,
            mod,
            TRAIN_END,
            VALIDATION_END,
        )

        e8_errors = recurrence_errors(
            E8,
            C8,
            mod,
            TRAIN_END,
            VALIDATION_END,
        )

        validation_length = (
            VALIDATION_END
            - TRAIN_END
            + 1
        )

        e6_error_rate = (
            e6_errors / validation_length
        )

        e8_error_rate = (
            e8_errors / validation_length
        )

        print("\n3. HELD-OUT RECURRENCE VALIDATION")
        print("-" * 78)
        print(
            f"E6 errors = "
            f"{e6_errors}/{validation_length} "
            f"rate={e6_error_rate:.6f}"
        )
        print(
            f"E8 errors = "
            f"{e8_errors}/{validation_length} "
            f"rate={e8_error_rate:.6f}"
        )

        # --------------------------------------------------------------
        # Complexity classification
        # --------------------------------------------------------------

        if L8 <= MAX_USEFUL_ORDER:
            complexity = "LOW"
        elif L8 <= TRAIN_END // 4:
            complexity = "MEDIUM"
        else:
            complexity = "HIGH"

        print(
            f"E8 recurrence complexity = "
            f"{complexity}"
        )

        # --------------------------------------------------------------
        # Period tests
        # --------------------------------------------------------------

        print("\n4. SIMPLE PERIOD TEST")
        print("-" * 78)

        period_candidates = [
            mod,
            2 * mod,
            mod * mod,
            2 * mod * mod,
            mod * (mod - 1),
        ]

        for period in period_candidates:

            if period >= LIMIT:
                continue

            errors = period_errors(
                E8,
                period,
                PERIOD_WARMUP,
            )

            print(
                f"period={period:6d} "
                f"errors={errors:6d}"
            )

        # --------------------------------------------------------------
        # Large-N extrapolation
        # --------------------------------------------------------------

        print("\n5. LARGE-N EXTRAPOLATION")
        print("-" * 78)

        if (
            L8 > 0
            and L8 <= MAX_USEFUL_ORDER
            and e8_errors == 0
        ):
            print(
                "VALIDATED LOW-ORDER RECURRENCE."
            )

            # This branch is intentionally conservative.
            print(
                "A recurrence exponentiation routine should be "
                "added only after this condition is observed."
            )

            for n in LARGE_N_TESTS:
                print(
                    f"N={n:16d}: "
                    "eligible for recurrence extrapolation"
                )

            large_n_status = True

        else:
            print(
                "NOT ELIGIBLE."
            )

            print(
                "Reason:"
            )

            if L8 == 0:
                print(
                    "  recurrence order is zero"
                )

            if L8 > MAX_USEFUL_ORDER:
                print(
                    f"  order {L8} > useful threshold "
                    f"{MAX_USEFUL_ORDER}"
                )

            if e8_errors != 0:
                print(
                    f"  held-out errors = {e8_errors}"
                )

            large_n_status = False

        # --------------------------------------------------------------
        # Occupancy
        # --------------------------------------------------------------

        counts: Dict[int, int] = {}

        for value in E8[
            PERIOD_WARMUP:
        ]:
            counts[value] = (
                counts.get(value, 0)
                + 1
            )

        top = sorted(
            counts.items(),
            key=lambda z: -z[1],
        )[:10]

        print("\n6. RESIDUE OCCUPANCY")
        print("-" * 78)
        print(
            f"distinct E8 residues = "
            f"{len(counts)}"
        )
        print(
            f"top residues = {top}"
        )

        summaries.append(
            {
                "mod": mod,
                "q_time": q_time,
                "L6": L6,
                "L8": L8,
                "E6_errors": e6_errors,
                "E8_errors": e8_errors,
                "E8_error_rate": e8_error_rate,
                "complexity": complexity,
                "large_n": large_n_status,
            }
        )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print("\n")
    print("=" * 78)
    print("7. CROSS-MODULUS SUMMARY")
    print("=" * 78)

    print(
        "ell | L(E6) | L(E8) | "
        "E8 valid errors | "
        "E8 error rate | complexity"
    )

    for row in summaries:

        print(
            f"{row['mod']:3d} | "
            f"{row['L6']:6d} | "
            f"{row['L8']:6d} | "
            f"{row['E8_errors']:15d} | "
            f"{row['E8_error_rate']:13.6f} | "
            f"{row['complexity']}"
        )

    strong = [
        row
        for row in summaries
        if (
            row["L8"] <= MAX_USEFUL_ORDER
            and row["E8_errors"] == 0
        )
    ]

    print("\n8. FINAL DIAGNOSTIC")
    print("-" * 78)

    if len(strong) >= 3:

        print(
            "STRONG RESULT:"
        )
        print(
            "Several finite fields have short exact recurrences "
            "that survive the held-out suffix."
        )
        print(
            "This would justify the next step: recurrence-based "
            "large-N modular evaluation."
        )

    elif len(strong) == 1 or len(strong) == 2:

        print(
            "PARTIAL RESULT:"
        )
        print(
            "A recurrence appears in some moduli but not robustly "
            "across the family."
        )

    else:

        print(
            "NO LOW-ORDER RECURRENCE SIGNAL."
        )
        print(
            "The ell=7 result already illustrates the issue:"
        )

        for row in summaries:
            if row["mod"] == 7:
                print(
                    f"  E8 order = {row['L8']}"
                )
                print(
                    f"  held-out errors = "
                    f"{row['E8_errors']}"
                )
                print(
                    f"  held-out error rate = "
                    f"{row['E8_error_rate']:.6f}"
                )

        print(
            "The sequence behaves as a high-complexity modular "
            "sequence rather than a short finite-state recurrence."
        )

    print()
    print(
        "RESEARCH IMPLICATION"
    )
    print(
        "A high BM order close to half the training length is not "
        "evidence for hidden structure. It is exactly what we expect "
        "when the observed sequence does not admit a short recurrence "
        "over the tested field."
    )

    print(
        "\nNo recurrence extrapolation is attempted unless the recurrence "
        "is both low-order and perfectly validated on held-out data."
    )

    total_runtime = (
        time.perf_counter()
        - overall_start
    )

    print(
        f"\ntotal runtime = {total_runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 83R2 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()