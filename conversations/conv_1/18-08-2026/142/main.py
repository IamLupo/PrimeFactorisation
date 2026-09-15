#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 142 — EXACT NORMALIZED COFACTOR / ADJUGATE AUDIT
==============================================================================

Continuation of Experiment 141.

The previous script failed because Fraction was not imported.

This replacement:

  * imports Fraction explicitly;
  * reconstructs the exact 12x12 F-core;
  * computes all 11x11 maximal minors;
  * divides them exactly by Delta_11;
  * builds the normalized cofactor matrix R;
  * checks the adjugate/determinant identity exactly;
  * audits row/column gcd structure;
  * audits prime valuations for 2, 3, 7;
  * checks selected modular ranks.

No floating point.
No SymPy.
No extrapolation.
No connection to the original (p,q)-kernel is assumed.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
import sys


# ============================================================================
# EXACT DATA
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}

D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}

TRANSITIONS = [
    (1, 3),
    (3, 5),
    (5, 7),
]

DELTA_11 = 9512681472

TERMINAL_D12 = (
    22698304379332609168018772364495783974334437093962503157858693669865148581087339234303140035372612096
)

TEST_PRIMES = [2, 3, 5, 7, 11, 13]


# ============================================================================
# BASIC HELPERS
# ============================================================================

def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


def valuation(n: int, p: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def gcd_many(values) -> int:
    g = 0

    for x in values:
        g = gcd(g, abs(int(x)))

    return abs(g)


def q_value(p: int, r: int) -> int:
    if r < 0:
        return 0

    if r >= len(Q[p]):
        return 0

    return Q[p][r]


def monomial_basis(p: int, d: int):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


# ============================================================================
# EXACT BAREISS DETERMINANT
# ============================================================================

def bareiss_det(A):
    if not A:
        return 1

    M = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(M)

    if n == 1:
        return M[0][0]

    previous = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):

            if M[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            M[k], M[pivot_row] = M[pivot_row], M[k]
            sign *= -1

        pivot = M[k][k]

        for i in range(k + 1, n):

            for j in range(k + 1, n):

                value = (
                    M[i][j] * pivot
                    - M[i][k] * M[k][j]
                )

                if k > 0:

                    if value % previous != 0:
                        raise ArithmeticError(
                            "Bareiss division was not exact."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD FULL 14x14 SYSTEM
# ============================================================================

def build_full_system():

    M = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            basis = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = basis[j] * q0

            for j in range(6):
                row[6 + j] = basis[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# BUILD 12x12 F-CORE
# ============================================================================

def build_core(M):

    boundary0 = M[5]
    boundary1 = M[10]

    H = [
        [boundary0[12], boundary0[13]],
        [boundary1[12], boundary1[13]],
    ]

    H_det = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if H_det == 0:
        raise ArithmeticError(
            "Boundary block is singular."
        )

    H_inv = [
        [
            Fraction(H[1][1], H_det),
            Fraction(-H[0][1], H_det),
        ],
        [
            Fraction(-H[1][0], H_det),
            Fraction(H[0][0], H_det),
        ],
    ]

    rhs0 = q_value(3, 5)
    rhs1 = q_value(5, 4)

    h_const = [
        H_inv[i][0] * rhs0
        + H_inv[i][1] * rhs1
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * boundary0[j]
                + H_inv[i][1] * boundary1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        reduced = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral core entry."
                )

            reduced.append(value.numerator)

        core.append(reduced)

    return core, H_det


# ============================================================================
# MAXIMAL MINORS
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [
        i for i in range(len(A))
        if i != omit_row
    ]

    cols = [
        j for j in range(len(A[0]))
        if j != omit_col
    ]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def normalized_cofactors(A):

    n = len(A)
    R = []

    for i in range(n):

        row = []

        for j in range(n):

            det_minor = bareiss_det(
                maximal_minor(A, i, j)
            )

            if det_minor % DELTA_11 != 0:
                raise ArithmeticError(
                    f"Delta_11 does not divide minor ({i},{j})."
                )

            row.append(
                det_minor // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MODULAR RANK
# ============================================================================

def rank_mod(A, p):

    M = [
        [x % p for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])
    rank = 0

    for col in range(cols):

        pivot = None

        for r in range(rank, rows):

            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = M[pivot], M[rank]

        inv = pow(
            M[rank][col],
            -1,
            p,
        )

        for j in range(col, cols):

            M[rank][j] = (
                M[rank][j] * inv
            ) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):

                M[r][j] = (
                    M[r][j]
                    - factor * M[rank][j]
                ) % p

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# ROW/COLUMN GCDS
# ============================================================================

def row_gcds(A):

    return [
        gcd_many(row)
        for row in A
    ]


def col_gcds(A):

    n = len(A[0])

    result = []

    for j in range(n):

        result.append(
            gcd_many(
                A[i][j]
                for i in range(len(A))
            )
        )

    return result


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 142 — EXACT NORMALIZED COFACTOR / "
        "ADJUGATE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in (1, 3, 5, 7):

        entries = len(Q[p])
        degree = -1

        for i in range(len(Q[p]) - 1, -1, -1):

            if Q[p][i] != 0:
                degree = i
                break

        expected = D[p]

        exact = (
            entries == expected + 1
            and degree == expected
        )

        print(
            f"  p={p}: entries={entries} "
            f"degree={degree} expected={expected} "
            f"exact={exact}"
        )

        data_exact = data_exact and exact

    print()
    print(
        f"  data_exact={data_exact}"
    )

    if not data_exact:
        print("\nDATA VALIDATION FAILED")
        return

    # ------------------------------------------------------------------
    # 2. CORE
    # ------------------------------------------------------------------

    M = build_full_system()
    core, H_det = build_core(M)

    det_core = bareiss_det(core)

    print()
    print("=" * 78)
    print("2. EXACT CORE")
    print("=" * 78)

    print(
        f"  core_shape=({len(core)},{len(core[0])})"
    )

    print(
        f"  H_det={H_det}"
    )

    print(
        f"  det_core_digits={digits(abs(det_core))}"
    )

    print(
        f"  det_core_sign="
        f"{'+' if det_core >= 0 else '-'}"
    )

    # ------------------------------------------------------------------
    # 3. ALL NORMALIZED COFACTORS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. NORMALIZED COFACTOR MATRIX")
    print("=" * 78)

    R = normalized_cofactors(core)

    print(
        "  R_ij = minor_ij / Delta_11"
    )

    print(
        f"  entry_count={len(R) * len(R[0])}"
    )

    max_digits = max(
        digits(x)
        for row in R
        for x in row
    )

    print(
        f"  max_entry_digits={max_digits}"
    )

    # ------------------------------------------------------------------
    # 4. COFACTOR CONTENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COFACTOR CONTENT")
    print("=" * 78)

    cofactor_content = gcd_many(
        x
        for row in R
        for x in row
    )

    print(
        f"  gcd_all_R={cofactor_content}"
    )

    print(
        f"  gcd_digits={digits(cofactor_content)}"
    )

    for p in TEST_PRIMES:

        v = valuation(
            cofactor_content,
            p,
        )

        if v:
            print(
                f"  prime={p}: valuation={v}"
            )

    # ------------------------------------------------------------------
    # 5. ROW / COLUMN PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. ROW / COLUMN GCD PROFILE")
    print("=" * 78)

    rg = row_gcds(R)
    cg = col_gcds(R)

    print()
    print("  row_gcds:")

    for i, value in enumerate(rg):

        print(
            f"    row={i}: gcd={value}"
        )

    print()
    print("  col_gcds:")

    for j, value in enumerate(cg):

        print(
            f"    col={j}: gcd={value}"
        )

    # ------------------------------------------------------------------
    # 6. PRIME VALUATION EXTREMES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. PRIME VALUATION EXTREMES")
    print("=" * 78)

    for p in (2, 3, 7):

        records = []

        for i in range(12):

            for j in range(12):

                value = R[i][j]

                if value == 0:
                    continue

                records.append(
                    (
                        valuation(value, p),
                        i,
                        j,
                    )
                )

        minimum = min(
            item[0]
            for item in records
        )

        minimum_positions = [
            (i, j)
            for v, i, j in records
            if v == minimum
        ]

        maximum = max(
            item[0]
            for item in records
        )

        print()
        print(
            f"  prime={p}"
        )

        print(
            f"    min_valuation={minimum}"
        )

        print(
            f"    min_positions={minimum_positions}"
        )

        print(
            f"    max_valuation={maximum}"
        )

    # ------------------------------------------------------------------
    # 7. MODULAR RANK OF R
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. NORMALIZED COFACTOR MODULAR RANK")
    print("=" * 78)

    for p in (2, 3, 7):

        rank = rank_mod(R, p)

        print(
            f"  prime={p}: rank_mod={rank}"
        )

    # ------------------------------------------------------------------
    # 8. EXACT ADJUGATE DETERMINANT IDENTITY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. EXACT ADJUGATE DETERMINANT IDENTITY")
    print("=" * 78)

    det_R = bareiss_det(R)

    # R = adj(A) / Delta_11.
    #
    # det(adj(A)) = det(A)^(11) for a 12x12 matrix.
    #
    # Therefore
    #
    # det(R)
    # =
    # det(A)^11 / Delta_11^12.

    expected_R = Fraction(
        abs(det_core) ** 11,
        DELTA_11 ** 12,
    )

    determinant_identity = (
        Fraction(det_R, 1)
        == expected_R
    )

    print(
        f"  det_R_digits={digits(abs(det_R))}"
    )

    print(
        f"  expected_R_is_integer="
        f"{expected_R.denominator == 1}"
    )

    print(
        f"  adjugate_identity_exact="
        f"{determinant_identity}"
    )

    # ------------------------------------------------------------------
    # 9. TERMINAL FACTOR RELATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. TERMINAL-SMITH RELATION")
    print("=" * 78)

    terminal_from_det = (
        abs(det_core) // DELTA_11
    )

    terminal_match = (
        terminal_from_det == TERMINAL_D12
    )

    print(
        f"  terminal_from_det={terminal_from_det}"
    )

    print(
        f"  terminal_reference={TERMINAL_D12}"
    )

    print(
        f"  terminal_match={terminal_match}"
    )

    for p in (2, 3, 7):

        v_delta = valuation(
            DELTA_11,
            p,
        )

        v_terminal = valuation(
            terminal_from_det,
            p,
        )

        v_det = valuation(
            abs(det_core),
            p,
        )

        exact = (
            v_delta + v_terminal
            == v_det
        )

        print(
            f"  prime={p}: "
            f"v_delta={v_delta} "
            f"v_terminal={v_terminal} "
            f"v_det={v_det} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The normalized cofactor matrix is an exact object derived from the
12x12 F-core:

    R_ij = minor_ij / Delta_11.

Because Delta_11 is the gcd of all maximal minors, R has integer
entries.

This experiment does not attempt a full Smith reduction.

Instead it asks whether the maximal-minor lattice has a secondary
structure visible directly in R:

    * common integer content;
    * row/column divisibility;
    * prime valuation localization;
    * modular rank defects;
    * exact adjugate determinant identity.

These are descriptive arithmetic properties.

No claim is made that any one of them is the hidden mathematical
mechanism producing the original kernel.

No connection to the original (p,q)-kernel is asserted.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    shape_ok = (
        len(R) == 12
        and len(R[0]) == 12
    )

    terminal_identity = (
        DELTA_11 * terminal_from_det
        == abs(det_core)
    )

    final_ok = (
        data_exact
        and shape_ok
        and DELTA_11 != 0
        and det_core != 0
        and terminal_identity
        and determinant_identity
        and terminal_match
    )

    print(
        f"  data_exact={data_exact}"
    )

    print(
        f"  normalized_cofactor_exact="
        f"{shape_ok}"
    )

    print(
        f"  Delta_11_divides_all_minors=True"
    )

    print(
        f"  adjugate_identity_exact="
        f"{determinant_identity}"
    )

    print(
        f"  terminal_identity_exact="
        f"{terminal_identity}"
    )

    print(
        f"  terminal_reference_match="
        f"{terminal_match}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 142 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise