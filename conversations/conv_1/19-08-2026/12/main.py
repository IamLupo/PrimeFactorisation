#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 181 — EXACT 17-ADIC LOCALIZATION / SCHUR-HIERARCHY AUDIT
==============================================================================

Experiment 180 established:

    gcd(q1(0), q3(0)) = 17

and

    gcd(s0,s1,s2) = 17

for the terminal Schur coefficients

    s0 = q1 - 2 q3
    s1 = q1 - 6 q3
    s2 = q1 - 18 q3.

The next question is whether the factor 17 is localized to the terminal
source pair or appears more broadly in the exact matrix hierarchy.

This experiment audits 17-adic content of:

    M        : full 14x14 system
    F_ret    : retained 12x12 block
    B        : 12x2 boundary-column block
    H        : 2x2 boundary block
    C        : 2x12 terminal-row block
    K        : exact Schur correction
    K_row_11 : distinguished terminal Schur row

and references:

    Delta_11
    d_12

For each object we compute:

    * gcd of all nonzero entries;
    * v17(gcd);
    * minimum finite v17 among nonzero entries;
    * maximum finite v17;
    * number of entries divisible by 17;
    * total number of nonzero entries.

The purpose is purely local:

    Does the factor 17 already occur in the bulk lattice?

or

    Is it concentrated in the terminal boundary/source layer?

No causal interpretation is attached.
No general theorem is inferred.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# DATA
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

FULL_N = 14
CORE_N = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)

P17 = 17


# ============================================================================
# HELPERS
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p, d):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def assert_rectangular(A, name="matrix"):
    if not A:
        return

    width = len(A[0])

    for i, row in enumerate(A):
        if len(row) != width:
            raise ValueError(
                f"{name}: row {i} width {len(row)}, expected {width}."
            )


def v_p(x, p):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % p == 0:
        x //= p
        e += 1

    return e


def gcd_nonzero(values):
    g = 0

    for x in values:
        if x != 0:
            g = gcd(g, abs(x))

    return g


def flatten(A):
    return [
        x
        for row in A
        for x in row
    ]


def matrix_stats(A, p, name):
    assert_rectangular(A, name)

    values = flatten(A)
    nonzero = [
        x
        for x in values
        if x != 0
    ]

    if not nonzero:
        return {
            "name": name,
            "entries": len(values),
            "nonzero": 0,
            "gcd": 0,
            "gcd_v": None,
            "min_v": None,
            "max_v": None,
            "divisible_count": 0,
            "divisible_fraction": 0.0,
        }

    valuations = [
        v_p(x, p)
        for x in nonzero
    ]

    divisible_count = sum(
        1
        for x in nonzero
        if x % p == 0
    )

    g = gcd_nonzero(values)

    return {
        "name": name,
        "entries": len(values),
        "nonzero": len(nonzero),
        "gcd": g,
        "gcd_v": v_p(g, p),
        "min_v": min(valuations),
        "max_v": max(valuations),
        "divisible_count": divisible_count,
        "divisible_fraction": (
            divisible_count / len(nonzero)
        ),
    }


def print_stats(stats):
    print(
        f"  {stats['name']}: "
        f"entries={stats['entries']} "
        f"nonzero={stats['nonzero']} "
        f"gcd={stats['gcd']} "
        f"v17(gcd)={stats['gcd_v']} "
        f"min_v17={stats['min_v']} "
        f"max_v17={stats['max_v']} "
        f"divisible_by_17={stats['divisible_count']}"
    )


# ============================================================================
# MATRIX OPERATIONS
# ============================================================================

def matmul(A, B):
    assert_rectangular(A, "left matrix")
    assert_rectangular(B, "right matrix")

    if not A or not B:
        return []

    if len(A[0]) != len(B):
        raise ValueError(
            f"Dimension mismatch: "
            f"{len(A)}x{len(A[0])} times "
            f"{len(B)}x{len(B[0])}."
        )

    out = [
        [0 for _ in range(len(B[0]))]
        for _ in range(len(A))
    ]

    for i in range(len(A)):
        for k in range(len(B)):

            if A[i][k] == 0:
                continue

            for j in range(len(B[0])):
                out[i][j] += (
                    A[i][k] * B[k][j]
                )

    return out


def inverse_2x2(H):
    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if abs(det_h) != 2:
        raise ArithmeticError(
            f"Expected |det(H)|=2, got {det_h}."
        )

    adj = [
        [H[1][1], -H[0][1]],
        [-H[1][0], H[0][0]],
    ]

    return det_h, adj


def exact_divide_matrix(A, d, name):
    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, x in enumerate(row):

            if x % d != 0:
                raise ArithmeticError(
                    f"{name}[{i}][{j}]={x} "
                    f"not divisible by {d}."
                )

            new_row.append(
                x // d
            )

        out.append(new_row)

    return out


# ============================================================================
# FULL SYSTEM
# ============================================================================

def build_full_system():
    M = []
    metadata = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * FULL_N

            for j in range(6):
                row[j] = (
                    b[j] * q0
                )

            for j in range(6):
                row[6 + j] = (
                    b[j] * q1
                )

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

            metadata.append(
                {
                    "full_row": len(M) - 1,
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "d": d,
                    "q0": q0,
                    "q1": q1,
                }
            )

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} rows, got {len(M)}."
        )

    return M, metadata


# ============================================================================
# SCHUR BLOCKS
# ============================================================================

def build_blocks(M):

    retained_indices = [
        i
        for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    retained = [
        M[i]
        for i in retained_indices
    ]

    boundary = [
        M[i]
        for i in BOUNDARY_ROWS
    ]

    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained
    ]

    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained
    ]

    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary
    ]

    C = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in boundary
    ]

    return F, B, H, C


def compute_schur(M):

    F, B, H, C = build_blocks(M)

    det_h, adj_h = inverse_2x2(H)

    numerator = matmul(
        B,
        matmul(
            adj_h,
            C,
        ),
    )

    if any(
        x % det_h != 0
        for row in numerator
        for x in row
    ):
        raise ArithmeticError(
            "Schur numerator not divisible by det(H)."
        )

    K = [
        [
            -(x // det_h)
            for x in row
        ]
        for row in numerator
    ]

    return {
        "F": F,
        "B": B,
        "H": H,
        "C": C,
        "adj_H": adj_h,
        "det_H": det_h,
        "K": K,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 181 — EXACT 17-ADIC LOCALIZATION / "
        "SCHUR-HIERARCHY AUDIT"
    )
    print("=" * 78)

    M, metadata = build_full_system()

    schur = compute_schur(M)

    F = schur["F"]
    B = schur["B"]
    H = schur["H"]
    C = schur["C"]
    adj_H = schur["adj_H"]
    det_H = schur["det_H"]
    K = schur["K"]

    # ------------------------------------------------------------------
    # 1. SOURCE TERMINAL VALUES
    # ------------------------------------------------------------------

    q1 = q_value(1, 5)
    q3 = q_value(3, 4)

    print()
    print("=" * 78)
    print("1. TERMINAL SOURCE VALUES")
    print("=" * 78)

    print(
        f"  q1=q_1(0)={q1} "
        f"v17={v_p(q1, 17)}"
    )

    print(
        f"  q3=q_3(0)={q3} "
        f"v17={v_p(q3, 17)}"
    )

    source_gcd = gcd(
        abs(q1),
        abs(q3),
    )

    print(
        f"  gcd(q1,q3)={source_gcd} "
        f"v17={v_p(source_gcd,17)}"
    )

    # ------------------------------------------------------------------
    # 2. BLOCK CONTENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. 17-ADIC BLOCK CONTENT")
    print("=" * 78)

    blocks = [
        ("M_full", M),
        ("F_ret", F),
        ("B", B),
        ("H", H),
        ("C", C),
        ("K", K),
    ]

    stats = {}

    for name, block in blocks:

        stats[name] = matrix_stats(
            block,
            17,
            name
        )

        print_stats(
            stats[name]
        )

    # ------------------------------------------------------------------
    # 3. DISTINGUISHED K ROW
    # ------------------------------------------------------------------

    K_row = K[11]

    K_row_stats = matrix_stats(
        [K_row],
        17,
        "K_row_11"
    )

    print()
    print("=" * 78)
    print("3. DISTINGUISHED TERMINAL SCHUR ROW")
    print("=" * 78)

    print_stats(
        K_row_stats
    )

    print(
        f"  K_row_11={K_row}"
    )

    # ------------------------------------------------------------------
    # 4. TERMINAL C ROWS
    # ------------------------------------------------------------------

    C0 = C[0]
    C1 = C[1]

    print()
    print("=" * 78)
    print("4. TERMINAL C-ROW 17-ADIC CONTENT")
    print("=" * 78)

    print_stats(
        matrix_stats(
            [C0],
            17,
            "C_row_0"
        )
    )

    print_stats(
        matrix_stats(
            [C1],
            17,
            "C_row_1"
        )
    )

    print(
        f"  gcd(C0)={gcd_nonzero(C0)}"
    )

    print(
        f"  gcd(C1)={gcd_nonzero(C1)}"
    )

    # ------------------------------------------------------------------
    # 5. BLOCKWISE 17 DIVISIBILITY FRACTIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. 17-DIVISIBILITY FRACTIONS")
    print("=" * 78)

    for name in [
        "M_full",
        "F_ret",
        "B",
        "H",
        "C",
        "K",
    ]:

        s = stats[name]

        print(
            f"  {name}: "
            f"{s['divisible_count']}/{s['nonzero']}"
        )

    # ------------------------------------------------------------------
    # 6. DELTA_11 / d12 REFERENCES
    # ------------------------------------------------------------------

    # Exact values from the established audit.
    Delta_11 = 9512681472

    d12 = (
        22698304379332609168018772364495783974334437093962503157858693669865148581087339234303140035372612096
    )

    print()
    print("=" * 78)
    print("6. SMITH TERMINAL REFERENCES")
    print("=" * 78)

    print(
        f"  Delta_11={Delta_11}"
    )

    print(
        f"  v17(Delta_11)={v_p(Delta_11,17)}"
    )

    print(
        f"  d12={d12}"
    )

    print(
        f"  v17(d12)={v_p(d12,17)}"
    )

    # ------------------------------------------------------------------
    # 7. SOURCE VS SCHUR CONTENT
    # ------------------------------------------------------------------

    schur_row_gcd = gcd_nonzero(
        K_row
    )

    print()
    print("=" * 78)
    print("7. SOURCE / SCHUR CONTENT COMPARISON")
    print("=" * 78)

    print(
        f"  source_gcd={source_gcd}"
    )

    print(
        f"  K_row_gcd={schur_row_gcd}"
    )

    print(
        f"  source_gcd_equals_K_row_gcd="
        f"{source_gcd == schur_row_gcd}"
    )

    print(
        f"  v17(source_gcd)="
        f"{v_p(source_gcd,17)}"
    )

    print(
        f"  v17(K_row_gcd)="
        f"{v_p(schur_row_gcd,17)}"
    )

    # ------------------------------------------------------------------
    # 8. DOES 17 APPEAR IN BULK?
    # ------------------------------------------------------------------

    full_bulk_has_17 = (
        stats["M_full"]["divisible_count"] > 0
    )

    full_bulk_all_17 = (
        stats["M_full"]["divisible_count"]
        == stats["M_full"]["nonzero"]
    )

    F_has_17 = (
        stats["F_ret"]["divisible_count"] > 0
    )

    F_all_17 = (
        stats["F_ret"]["divisible_count"]
        == stats["F_ret"]["nonzero"]
    )

    print()
    print("=" * 78)
    print("8. BULK-LOCALIZATION CONTROL")
    print("=" * 78)

    print(
        f"  full_system_has_17_entries="
        f"{full_bulk_has_17}"
    )

    print(
        f"  all_full_system_nonzero_entries_divisible_by_17="
        f"{full_bulk_all_17}"
    )

    print(
        f"  F_core_has_17_entries="
        f"{F_has_17}"
    )

    print(
        f"  all_F_core_nonzero_entries_divisible_by_17="
        f"{F_all_17}"
    )

    # ------------------------------------------------------------------
    # 9. 17-REMOVED TERMINAL SOURCE
    # ------------------------------------------------------------------

    q1_17 = q1 // 17
    q3_17 = q3 // 17

    terminal_s = [
        q1 - 2 * q3,
        q1 - 6 * q3,
        q1 - 18 * q3,
    ]

    terminal_s_17 = [
        x // 17
        for x in terminal_s
    ]

    terminal_17_exact = all(
        x % 17 == 0
        for x in terminal_s
    )

    print()
    print("=" * 78)
    print("9. TERMINAL SOURCE FACTOR REMOVAL")
    print("=" * 78)

    print(
        f"  q1/17={q1_17}"
    )

    print(
        f"  q3/17={q3_17}"
    )

    print(
        f"  terminal_coefficients={terminal_s}"
    )

    print(
        f"  all_terminal_coefficients_divisible_by_17="
        f"{terminal_17_exact}"
    )

    print(
        f"  terminal_coefficients_divided_by_17="
        f"{terminal_s_17}"
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
The terminal experiments established an exact common factor 17 in the
pair

    (q_1(0), q_3(0))

and consequently in the three terminal Schur coefficients.

Experiment 181 asks a different question:

    Is 17 a bulk lattice factor, or is it localized to the terminal
    source layer?

If the full matrix or F-core has gcd 1 while the terminal source pair
and distinguished Schur row have gcd 17, then the evidence supports a
localized interpretation.

If 17 divides a substantial fraction of the bulk entries, then the
terminal factor is less exceptional.

The Smith references Delta_11 and d12 are included to determine whether
17 survives into the higher determinantal invariants.

This experiment does not identify why 17 occurs. It only localizes
where it occurs in the exact hierarchy.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    source_gcd_ok = (
        source_gcd == 17
    )

    K_row_gcd_ok = (
        schur_row_gcd == 17
    )

    terminal_ok = (
        terminal_17_exact
    )

    full_content_one = (
        stats["M_full"]["gcd"] == 1
    )

    F_content_one = (
        stats["F_ret"]["gcd"] == 1
    )

    final_ok = (
        source_gcd_ok
        and K_row_gcd_ok
        and terminal_ok
        and full_content_one
        and F_content_one
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_gcd_is_17={source_gcd_ok}"
    )

    print(
        f"  K_row_gcd_is_17={K_row_gcd_ok}"
    )

    print(
        f"  terminal_factor_17_exact={terminal_ok}"
    )

    print(
        f"  full_matrix_content_is_1={full_content_one}"
    )

    print(
        f"  F_core_content_is_1={F_content_one}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 181 COMPLETE")


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

