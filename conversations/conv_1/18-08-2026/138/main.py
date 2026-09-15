#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 138 — EXACT F-CORE SMITH / PRIME-SPECTRUM / CONTENT AUDIT
==============================================================================

Evidence-first continuation of Experiments 133–137.

Experiment 137 established an exact decomposition

    full 14x14 system
        ->
    12x12 boundary-eliminated F-core.

Experiment 138 studies that F-core independently.

Questions:

    1. What are its exact determinantal divisors?
    2. What is its terminal Smith factor?
    3. Which primes cause rank loss modulo p?
    4. Does the core have a nontrivial common integer content?
    5. Are the core invariants preserved by unimodular changes?

No assumption is made about the origin of the large integers.

Everything is exact integer/Fraction arithmetic.
No floating point.
No SymPy.
No extrapolation.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd
import sys


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


# ============================================================================
# BASIC HELPERS
# ============================================================================

def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


def lcm(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def poly_degree(values: list[int]) -> int:
    for i in range(len(values) - 1, -1, -1):
        if values[i] != 0:
            return i
    return -1


def valuation(n: int, p: int) -> int:
    n = abs(n)
    value = 0

    while n and n % p == 0:
        n //= p
        value += 1

    return value


# ============================================================================
# EXACT MATRIX RANK
# ============================================================================

def matrix_rank(A) -> int:
    if not A:
        return 0

    M = [
        [Fraction(x) for x in row]
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

        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]

        pv = M[rank][col]

        for j in range(col, cols):
            M[rank][j] /= pv

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):
                M[r][j] -= factor * M[rank][j]

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# EXACT MODULAR RANK
# ============================================================================

def mod_rank(A: list[list[int]], p: int) -> int:
    if not A:
        return 0

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
            if M[r][col] % p != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = M[pivot], M[rank]

        inv = pow(M[rank][col], -1, p)

        for j in range(col, cols):
            M[rank][j] = (M[rank][j] * inv) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col] % p

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
# BAREISS DETERMINANT
# ============================================================================

def bareiss_det(A: list[list[int]]) -> int:
    if not A:
        return 1

    B = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(B)

    if n == 1:
        return B[0][0]

    previous = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):
            if B[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            B[k], B[pivot_row] = B[pivot_row], B[k]
            sign *= -1

        pivot = B[k][k]

        for i in range(k + 1, n):

            for j in range(k + 1, n):

                value = (
                    B[i][j] * pivot
                    - B[i][k] * B[k][j]
                )

                if k > 0:

                    if value % previous != 0:
                        raise ArithmeticError(
                            "Bareiss division was not exact."
                        )

                    value //= previous

                B[i][j] = value

        for i in range(k + 1, n):
            B[i][k] = 0

        previous = pivot

    return sign * B[n - 1][n - 1]


# ============================================================================
# MINOR
# ============================================================================

def minor_matrix(
    A: list[list[int]],
    rows,
    cols,
) -> list[list[int]]:
    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


# ============================================================================
# q ACCESS
# ============================================================================

def q_value(p: int, r: int) -> int:
    if r < 0:
        return 0

    row = Q[p]

    if r >= len(row):
        return 0

    return row[r]


# ============================================================================
# BASIS
# ============================================================================

def monomial_basis(p: int, d: int) -> list[int]:
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


# ============================================================================
# ROW METADATA
# ============================================================================

def row_metadata(index: int) -> tuple[int, int, int]:

    if 0 <= index < 6:
        return 1, 3, index

    if 6 <= index < 11:
        return 3, 5, index - 6

    if 11 <= index < 14:
        return 5, 7, index - 11

    raise IndexError(index)


# ============================================================================
# BUILD FULL SYSTEM
# ============================================================================

def build_full_system() -> list[list[int]]:

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
# BOUNDARY ELIMINATION
# ============================================================================

def build_core(
    M: list[list[int]],
) -> tuple[list[list[int]], int]:

    # Boundary equations:
    # p=1,d=0 -> row 5
    # p=3,d=0 -> row 10

    row_a = M[5]
    row_b = M[10]

    H = [
        [row_a[12], row_a[13]],
        [row_b[12], row_b[13]],
    ]

    H_det = bareiss_det(H)

    if H_det == 0:
        raise ValueError(
            "Boundary H block is singular."
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

    p_a, p_a_next, r_a = row_metadata(5)
    p_b, p_b_next, r_b = row_metadata(10)

    rhs = [
        q_value(p_a_next, r_a),
        q_value(p_b_next, r_b),
    ]

    h_const = [
        H_inv[i][0] * rhs[0]
        + H_inv[i][1] * rhs[1]
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    A_boundary = [
        row_a[:12],
        row_b[:12],
    ]

    for i in range(2):
        for j in range(12):
            h_coeff[i][j] = -(
                H_inv[i][0] * A_boundary[0][j]
                + H_inv[i][1] * A_boundary[1][j]
            )

    core = []

    removed = {5, 10}

    for idx in range(14):

        if idx in removed:
            continue

        row = M[idx]

        p, p_next, r = row_metadata(idx)

        target = q_value(p_next, r)

        reduced = []

        for j in range(12):
            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Reduced F-core is not integral."
                )

            reduced.append(value.numerator)

        reduced.append(
            Fraction(target)
            - Fraction(row[12]) * h_const[0]
            - Fraction(row[13]) * h_const[1]
        )

        core.append(reduced)

    A = [
        row[:12]
        for row in core
    ]

    # Verify exact integral RHS independently.
    for row in core:
        if row[12].denominator != 1:
            raise ArithmeticError(
                "Reduced core RHS is not integral."
            )

    A_integer = A
    return A_integer, H_det


# ============================================================================
# DETERMINANTAL DIVISOR
# ============================================================================

def delta_11_core(A: list[list[int]]) -> int:

    if len(A) != 12 or len(A[0]) != 12:
        raise ValueError("Expected 12x12 matrix.")

    g = 0

    rows = range(12)

    for omitted_row in range(12):

        kept_rows = [
            i for i in rows
            if i != omitted_row
        ]

        for omitted_col in range(12):

            kept_cols = [
                j for j in range(12)
                if j != omitted_col
            ]

            value = bareiss_det(
                minor_matrix(
                    A,
                    kept_rows,
                    kept_cols,
                )
            )

            g = gcd(g, abs(value))

            if g == 1:
                return 1

    return abs(g)


# ============================================================================
# ENTRY CONTENT
# ============================================================================

def entry_content(A: list[list[int]]) -> int:

    g = 0

    for row in A:
        for x in row:
            g = gcd(g, abs(int(x)))

    return abs(g)


# ============================================================================
# UNIMODULAR TRANSFORMS
# ============================================================================

def row_add(
    A: list[list[int]],
    target: int,
    source: int,
    multiplier: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    for j in range(len(B[0])):
        B[target][j] += multiplier * B[source][j]

    return B


def col_add(
    A: list[list[int]],
    target: int,
    source: int,
    multiplier: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    for i in range(len(B)):
        B[i][target] += multiplier * B[i][source]

    return B


def row_swap(
    A: list[list[int]],
    i: int,
    j: int,
) -> list[list[int]]:

    B = [row[:] for row in A]
    B[i], B[j] = B[j], B[i]
    return B


def col_swap(
    A: list[list[int]],
    i: int,
    j: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    for row in B:
        row[i], row[j] = row[j], row[i]

    return B


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print(
        "EXPERIMENT 138 — EXACT F-CORE SMITH / "
        "PRIME-SPECTRUM / CONTENT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. DATA VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_exact = True

    for p in (1, 3, 5, 7):

        entries = len(Q[p])
        deg = poly_degree(Q[p])
        expected = D[p]

        exact = (
            entries == expected + 1
            and deg == expected
        )

        print(
            f"  p={p}: entries={entries} "
            f"degree={deg} expected={expected} "
            f"exact={exact}"
        )

        data_exact = data_exact and exact

    print()
    print(f"  data_exact={data_exact}")

    if not data_exact:
        print("\nDATA VALIDATION FAILED")
        return

    # ------------------------------------------------------------------
    # 2. FULL SYSTEM / CORE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FULL SYSTEM AND EXACT F-CORE")
    print("=" * 78)

    M = build_full_system()

    full_rank = matrix_rank(M)
    full_det = bareiss_det(M)

    core_with_rhs, H_det = build_core(M)

    # core_with_rhs is 12 rows x 13 columns.
    core = [
        row[:12]
        for row in core_with_rhs
    ]

    core_rank = matrix_rank(core)
    core_det = bareiss_det(core)

    determinant_identity = (
        abs(full_det) == abs(H_det * core_det)
    )

    print(f"  full_rank={full_rank}")
    print(
        f"  full_det_digits={digits(abs(full_det))}"
    )
    print(f"  H_det={H_det}")
    print(
        f"  core_shape=({len(core)},{len(core[0])})"
    )
    print(f"  core_rank={core_rank}")
    print(
        f"  core_det_digits={digits(abs(core_det))}"
    )
    print(
        f"  determinant_identity={determinant_identity}"
    )

    # ------------------------------------------------------------------
    # 3. CORE SMITH DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT CORE SMITH DATA")
    print("=" * 78)

    delta11 = delta_11_core(core)
    d12 = abs(core_det) // delta11

    print(f"  Delta_11={delta11}")
    print(
        f"  Delta_11_digits={digits(delta11)}"
    )
    print(f"  d_12={d12}")
    print(
        f"  d_12_digits={digits(d12)}"
    )

    print()
    print("  prime valuations:")

    for prime in (
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31,
    ):
        v_delta = valuation(delta11, prime)
        v_last = valuation(d12, prime)

        if v_delta or v_last:
            print(
                f"    prime={prime}: "
                f"v(Delta_11)={v_delta}, "
                f"v(d_12)={v_last}"
            )

    # ------------------------------------------------------------------
    # 4. CORE ENTRY CONTENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CORE INTEGER-CONTENT AUDIT")
    print("=" * 78)

    content = entry_content(core)

    print(
        f"  gcd_of_all_core_entries={content}"
    )
    print(
        f"  content_digits={digits(content)}"
    )

    content_removal_exact = True

    if content != 0:

        primitive = [
            [
                x // content
                for x in row
            ]
            for row in core
        ]

        primitive_det = bareiss_det(primitive)

        expected_det = (
            abs(primitive_det)
            * content ** 12
        )

        content_removal_exact = (
            abs(core_det) == expected_det
        )

        print(
            f"  primitive_core_det_digits="
            f"{digits(abs(primitive_det))}"
        )

        print(
            f"  det_content_scaling_exact="
            f"{content_removal_exact}"
        )

    # ------------------------------------------------------------------
    # 5. MODULAR RANK SPECTRUM
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. CORE MODULAR RANK SPECTRUM")
    print("=" * 78)

    primes = [
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31,
        37, 41, 43, 47,
        53, 59, 61, 67, 71,
        73, 79, 83, 89, 97,
    ]

    rank_defect_primes = []

    for prime in primes:

        rank_mod = mod_rank(core, prime)

        defect = rank_mod < 12

        if defect:
            rank_defect_primes.append(prime)

        print(
            f"  prime={prime}: "
            f"rank_mod_p={rank_mod} "
            f"defect={defect}"
        )

    print()
    print(
        f"  rank_defect_primes={rank_defect_primes}"
    )

    # ------------------------------------------------------------------
    # 6. UNIMODULAR INVARIANCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. UNIMODULAR CORE INVARIANCE")
    print("=" * 78)

    transforms = [
        (
            "row_add(0,1,+3)",
            lambda A: row_add(A, 0, 1, 3),
        ),
        (
            "row_add(5,2,-7)",
            lambda A: row_add(A, 5, 2, -7),
        ),
        (
            "col_add(3,7,+5)",
            lambda A: col_add(A, 3, 7, 5),
        ),
        (
            "col_add(10,4,-4)",
            lambda A: col_add(A, 10, 4, -4),
        ),
        (
            "row_swap(0,11)",
            lambda A: row_swap(A, 0, 11),
        ),
        (
            "col_swap(1,10)",
            lambda A: col_swap(A, 1, 10),
        ),
    ]

    invariant_ok = True

    for name, transform in transforms:

        T = transform(core)

        det_T = bareiss_det(T)
        delta_T = delta_11_core(T)
        d_T = abs(det_T) // delta_T

        same_det = abs(det_T) == abs(core_det)
        same_delta = delta_T == delta11
        same_last = d_T == d12

        exact = (
            same_det
            and same_delta
            and same_last
        )

        print()
        print(f"  {name}")
        print(f"    same_abs_det={same_det}")
        print(f"    same_Delta_11={same_delta}")
        print(f"    same_d_12={same_last}")
        print(f"    invariant_exact={exact}")

        invariant_ok = invariant_ok and exact

    # ------------------------------------------------------------------
    # 7. FULL VS CORE COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. FULL VS CORE TERMINAL FACTOR")
    print("=" * 78)

    # Full Delta_13
    full_delta13 = 0

    for i in range(14):

        for j in range(14):

            kept_rows = [
                r for r in range(14)
                if r != i
            ]

            kept_cols = [
                c for c in range(14)
                if c != j
            ]

            value = bareiss_det(
                minor_matrix(
                    M,
                    kept_rows,
                    kept_cols,
                )
            )

            full_delta13 = gcd(
                full_delta13,
                abs(value),
            )

            if full_delta13 == 1:
                break

        if full_delta13 == 1:
            break

    full_d14 = abs(full_det) // full_delta13

    ratio_full_core = Fraction(
        full_d14,
        d12,
    )

    print(f"  full_Delta_13={full_delta13}")
    print(
        f"  full_d_14_digits={digits(full_d14)}"
    )
    print(f"  core_Delta_11={delta11}")
    print(
        f"  core_d_12_digits={digits(d12)}"
    )
    print(
        f"  full_d14_over_core_d12="
        f"{ratio_full_core}"
    )

    # ------------------------------------------------------------------
    # 8. INTERIOR REFERENCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. INTERIOR REFERENCE")
    print("=" * 78)

    boundary_rows = {5, 10, 13}

    interior_rows = [
        i
        for i in range(14)
        if i not in boundary_rows
    ]

    interior = [
        [
            M[i][j]
            for j in range(12)
        ]
        for i in interior_rows
    ]

    interior_rank = matrix_rank(interior)

    interior_minor_gcd = 0

    # 11 x 12 -> twelve maximal 11x11 minors.
    for omitted_col in range(12):

        cols = [
            j
            for j in range(12)
            if j != omitted_col
        ]

        sub = [
            [interior[i][j] for j in cols]
            for i in range(11)
        ]

        value = bareiss_det(sub)

        interior_minor_gcd = gcd(
            interior_minor_gcd,
            abs(value),
        )

    print(
        f"  interior_shape="
        f"({len(interior)},{len(interior[0])})"
    )
    print(f"  interior_rank={interior_rank}")
    print(
        f"  interior_maximal_minor_gcd="
        f"{interior_minor_gcd}"
    )
    print(
        f"  interior_gcd_digits="
        f"{digits(interior_minor_gcd)}"
    )

    # ------------------------------------------------------------------
    # 9. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
This experiment intentionally makes no causal claim about the large
coefficients.

The exact hierarchy is:

    full 14x14 system
        |
        | eliminate the two H variables
        v
    12x12 F-core
        |
        | exact Smith analysis
        v
    Delta_11, d_12.

The F-core is therefore examined as an independent arithmetic object.

Three different quantities are separated:

    entry content:
        gcd of all 144 core entries;

    Delta_11:
        gcd of all 11x11 core minors;

    d_12:
        |det(core)| / Delta_11.

These measure increasingly global lattice structure.

The modular rank spectrum gives an independent prime-level view:
a prime is structurally visible when the core loses rank modulo that
prime.

Unimodular transformations test that the observed determinant and
Smith terminal invariant are coordinate invariants.

The comparison with the full-system d_14 is descriptive only.
No hypothesis is imposed that one must explain the other in a
particular way.

No connection to the original (p,q)-kernel is asserted.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    full_exact = (
        full_rank == 14
        and full_det != 0
    )

    core_exact = (
        core_rank == 12
        and core_det != 0
    )

    smith_exact = (
        delta11 != 0
        and d12 != 0
        and abs(core_det) == delta11 * d12
    )

    determinant_exact = determinant_identity

    final_ok = (
        data_exact
        and full_exact
        and core_exact
        and smith_exact
        and determinant_exact
        and content_removal_exact
        and invariant_ok
        and interior_rank == 11
        and interior_minor_gcd != 0
    )

    print(f"  data_exact={data_exact}")
    print(f"  full_system_exact={full_exact}")
    print(f"  core_exact={core_exact}")
    print(f"  core_smith_exact={smith_exact}")
    print(
        f"  determinant_identity="
        f"{determinant_exact}"
    )
    print(
        f"  content_audit_exact="
        f"{content_removal_exact}"
    )
    print(
        f"  modular_audit_completed=True"
    )
    print(
        f"  unimodular_invariance="
        f"{invariant_ok}"
    )
    print(
        f"  interior_reference_exact="
        f"{interior_rank == 11 and interior_minor_gcd != 0}"
    )
    print(
        f"  failures={0 if final_ok else 1}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 138 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:
        print(
            f"\nFATAL ERROR: {type(exc).__name__}: {exc}"
        )
        raise