#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 183 — EXACT 17-ADIC SOURCE-TABLE / FULL-SYSTEM PROVENANCE AUDIT
==============================================================================

Experiment 182 established the exact transfer

    C = 17 C'
    K = 17 K'

with C' and K' primitive, while

    B, H, Delta_11, d12

carry no factor 17.

The next question is:

    Where does the factor 17 first enter the original data?

The full matrix contains terms of the form

    basis(p,d) * q_p(r)

and

    basis(p,d) * q_p(r+1).

For the present dataset,

    p <= 5,
    d <= 5,

so none of

    1, d, d^2, p, p*d, p^2

is divisible by 17 unless the factor were already in the q-data.

Experiment 183 therefore audits the source table itself.

For every q_p(r):

    * exact value;
    * residue mod 17;
    * v17;
    * odd/content information.

It then groups all q entries divisible by 17 and counts how many full-system
entries they generate.

The final provenance test asks whether every nonzero full-system entry
divisible by 17 is explained by a q-entry divisible by 17, with no
additional 17 introduced by the basis.

This is a source-localization experiment only.

No theorem is inferred from the finite table.
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

P17 = 17


# ============================================================================
# HELPERS
# ============================================================================

def v_p(x, p):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % p == 0:
        x //= p
        e += 1

    return e


def basis(p, d):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def gcd_list(values):
    g = 0

    for x in values:
        g = gcd(g, abs(x))

    return g


def assert_rectangular(A, name="matrix"):
    if not A:
        return

    w = len(A[0])

    for i, row in enumerate(A):

        if len(row) != w:
            raise ValueError(
                f"{name}: row {i} has width {len(row)}, "
                f"expected {w}."
            )


def matrix_support(A):
    return [
        (i, j)
        for i in range(len(A))
        for j in range(len(A[0]))
        if A[i][j] != 0
    ]


# ============================================================================
# FULL SYSTEM
# ============================================================================

def build_full_system():

    M = []
    provenance = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * FULL_N
            source = [None] * FULL_N

            # First six coordinates come from q_p(r).
            for j in range(6):

                row[j] = b[j] * q0

                source[j] = {
                    "q_p": p,
                    "q_r": r,
                    "q_value": q0,
                    "basis_index": j,
                    "basis_value": b[j],
                    "kind": "q(r)",
                }

            # Second six coordinates come from q_p(r+1).
            for j in range(6):

                row[6 + j] = b[j] * q1

                source[6 + j] = {
                    "q_p": p,
                    "q_r": r + 1,
                    "q_value": q1,
                    "basis_index": j,
                    "basis_value": b[j],
                    "kind": "q(r+1)",
                }

            if d == 0:

                row[12] = 1
                row[13] = p

                source[12] = {
                    "kind": "boundary_1",
                }

                source[13] = {
                    "kind": "boundary_p",
                }

            M.append(row)

            provenance.append(
                {
                    "full_row": len(M) - 1,
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "d": d,
                    "sources": source,
                }
            )

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} rows, got {len(M)}."
        )

    return M, provenance


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 183 — EXACT 17-ADIC SOURCE-TABLE / "
        "FULL-SYSTEM PROVENANCE AUDIT"
    )
    print("=" * 78)

    M, provenance = build_full_system()

    # ------------------------------------------------------------------
    # 1. COMPLETE q TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. COMPLETE q-TABLE 17-ADIC PROFILE")
    print("=" * 78)

    q_entries = []

    for p in sorted(Q):

        row = Q[p]

        row_gcd = gcd_list(row)

        print(
            f"  p={p}: "
            f"count={len(row)} "
            f"gcd={row_gcd} "
            f"v17(gcd)={v_p(row_gcd,17)}"
        )

        for r, value in enumerate(row):

            residue = value % P17
            valuation = v_p(value, P17)

            item = (
                p,
                r,
                value,
                residue,
                valuation,
            )

            q_entries.append(item)

            print(
                f"    r={r}: "
                f"value={value} "
                f"mod17={residue} "
                f"v17={valuation}"
            )

    # ------------------------------------------------------------------
    # 2. q ENTRIES DIVISIBLE BY 17
    # ------------------------------------------------------------------

    divisible_q = [
        item
        for item in q_entries
        if item[4] is not None and item[4] >= 1
    ]

    print()
    print("=" * 78)
    print("2. q-ENTRIES DIVISIBLE BY 17")
    print("=" * 78)

    print(
        f"  count={len(divisible_q)}"
    )

    for p, r, value, residue, valuation in divisible_q:

        print(
            f"  (p={p}, r={r}): "
            f"value={value} "
            f"v17={valuation}"
        )

    # ------------------------------------------------------------------
    # 3. NONZERO SOURCE CONTENT
    # ------------------------------------------------------------------

    all_q_values = [
        value
        for _, _, value, _, _ in q_entries
        if value != 0
    ]

    global_q_gcd = gcd_list(
        all_q_values
    )

    print()
    print("=" * 78)
    print("3. GLOBAL q-TABLE CONTENT")
    print("=" * 78)

    print(
        f"  global_q_gcd={global_q_gcd}"
    )

    print(
        f"  v17(global_q_gcd)="
        f"{v_p(global_q_gcd,17)}"
    )

    # ------------------------------------------------------------------
    # 4. BASIS 17-ADIC PROFILE
    # ------------------------------------------------------------------

    basis_values = []

    for p in [1, 3, 5, 7]:

        for d in range(
            D[p] + 1
        ):

            b = basis(p, d)

            for j, value in enumerate(b):

                basis_values.append(
                    (p, d, j, value)
                )

    max_basis_v17 = max(
        v_p(value, 17) or 0
        for _, _, _, value in basis_values
    )

    basis_divisible = [
        item
        for item in basis_values
        if item[3] % 17 == 0
    ]

    print()
    print("=" * 78)
    print("4. BASIS 17-ADIC CONTROL")
    print("=" * 78)

    print(
        f"  total_basis_values="
        f"{len(basis_values)}"
    )

    print(
        f"  max_v17_basis={max_basis_v17}"
    )

    print(
        f"  basis_values_divisible_by_17="
        f"{len(basis_divisible)}"
    )

    for p, d, j, value in basis_divisible:

        print(
            f"  p={p} d={d} basis_index={j} "
            f"value={value}"
        )

    # ------------------------------------------------------------------
    # 5. FULL SYSTEM 17-DIVISIBLE ENTRIES
    # ------------------------------------------------------------------

    divisible_matrix_entries = []

    for i, row in enumerate(M):

        for j, value in enumerate(row):

            if value != 0 and value % 17 == 0:

                src = provenance[i]["sources"][j]

                divisible_matrix_entries.append(
                    (
                        i,
                        j,
                        value,
                        src,
                    )
                )

    print()
    print("=" * 78)
    print("5. FULL-SYSTEM ENTRIES DIVISIBLE BY 17")
    print("=" * 78)

    print(
        f"  count={len(divisible_matrix_entries)}"
    )

    for i, j, value, src in divisible_matrix_entries:

        print(
            f"  matrix=({i},{j}) "
            f"value={value} "
            f"source={src}"
        )

    # ------------------------------------------------------------------
    # 6. SOURCE GROUPING
    # ------------------------------------------------------------------

    source_counts = {}

    for i, j, value, src in divisible_matrix_entries:

        if src["kind"] not in (
            "q(r)",
            "q(r+1)",
        ):
            key = (
                src["kind"],
            )

        else:
            key = (
                src["q_p"],
                src["q_r"],
            )

        source_counts[key] = (
            source_counts.get(key, 0)
            + 1
        )

    print()
    print("=" * 78)
    print("6. 17-DIVISIBLE MATRIX ENTRIES BY SOURCE q")
    print("=" * 78)

    for key in sorted(
        source_counts,
        key=str
    ):

        print(
            f"  source={key} "
            f"matrix_entry_count="
            f"{source_counts[key]}"
        )

    # ------------------------------------------------------------------
    # 7. SOURCE-EXPLANATION TEST
    # ------------------------------------------------------------------

    unexplained = []

    for i, j, value, src in divisible_matrix_entries:

        if src["kind"] in (
            "q(r)",
            "q(r+1)",
        ):

            qv = src["q_value"]
            basisv = src["basis_value"]

            if qv % 17 != 0:
                unexplained.append(
                    (
                        i,
                        j,
                        value,
                        qv,
                        basisv,
                    )
                )

        else:

            unexplained.append(
                (
                    i,
                    j,
                    value,
                    src["kind"],
                )
            )

    explained_exactly = (
        len(unexplained) == 0
    )

    print()
    print("=" * 78)
    print("7. SOURCE PROVENANCE TEST")
    print("=" * 78)

    print(
        f"  unexplained_entries="
        f"{len(unexplained)}"
    )

    print(
        f"  every_17_entry_explained_by_q="
        f"{explained_exactly}"
    )

    if unexplained:

        for item in unexplained:
            print(
                f"  unexplained={item}"
            )

    # ------------------------------------------------------------------
    # 8. TERMINAL q SOURCE
    # ------------------------------------------------------------------

    terminal_sources = [
        (1, 5, q_value(1, 5)),
        (3, 4, q_value(3, 4)),
    ]

    print()
    print("=" * 78)
    print("8. TERMINAL q SOURCES")
    print("=" * 78)

    for p, r, value in terminal_sources:

        print(
            f"  p={p} r={r} "
            f"value={value} "
            f"v17={v_p(value,17)}"
        )

    # ------------------------------------------------------------------
    # 9. TERMINAL VS NONTERMINAL 17-DIVISIBILITY
    # ------------------------------------------------------------------

    terminal_divisible_q = [
        item
        for item in divisible_q
        if (
            item[0] == 1 and item[1] == 5
        )
        or (
            item[0] == 3 and item[1] == 4
        )
    ]

    nonterminal_divisible_q = [
        item
        for item in divisible_q
        if item not in terminal_divisible_q
    ]

    print()
    print("=" * 78)
    print("9. TERMINAL / NONTERMINAL SPLIT")
    print("=" * 78)

    print(
        f"  terminal_divisible_q_count="
        f"{len(terminal_divisible_q)}"
    )

    print(
        f"  nonterminal_divisible_q_count="
        f"{len(nonterminal_divisible_q)}"
    )

    for item in nonterminal_divisible_q:

        print(
            f"  nonterminal_17_source="
            f"(p={item[0]}, r={item[1]}) "
            f"value={item[2]} "
            f"v17={item[4]}"
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
Experiment 182 showed exact 17-content transfer

    C = 17 C'
    K = 17 K'.

Experiment 183 moves upstream to the q-table.

Because all basis values in the present finite system have
v17 = 0, a 17-divisible matrix coefficient can only arise from
a 17-divisible q-value.

The source-provenance audit therefore distinguishes:

    terminal-only source divisibility

from

    broader q-table divisibility.

If only the terminal q-values q_1(0) and q_3(0) are divisible by 17,
then the factor 17 is highly localized at the terminal source.

If additional q_p(r) values are divisible by 17, the factor is still
source-level but not terminal-only.

No attempt is made to infer a formula for future values of q.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    basis_clean = (
        len(basis_divisible) == 0
        and max_basis_v17 == 0
    )

    q_source_explains = (
        explained_exactly
    )

    final_ok = (
        basis_clean
        and q_source_explains
        and len(divisible_q) == 3
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  basis_has_no_17_factor="
        f"{basis_clean}"
    )

    print(
        f"  q-source_explains_all_matrix_17_entries="
        f"{q_source_explains}"
    )

    print(
        f"  divisible_q_count_is_3="
        f"{len(divisible_q) == 3}"
    )

    print(
        f"  terminal_sources_are_17_divisible="
        f"{len(terminal_divisible_q) == 2}"
    )

    print(
        f"  nonterminal_divisible_q_count="
        f"{len(nonterminal_divisible_q)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 183 COMPLETE")


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

