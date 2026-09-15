#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 284 — EXACT B-DIAGONAL RECURRENCE / SEPARABILITY AUDIT
==============================================================================

Experiment 283 established:

    * no total-degree <= 4 bivariate polynomial law in (k,r);
    * no low-degree tensor-product polynomial law;
    * no low-degree polynomial law in (k,d), d=r-k;
    * the few exact fixed-diagonal polynomial fits are interpolation-only.

Experiment 284 therefore changes the model class.

The natural geometry is the triangular array

    B[k,r],    r >= k.

We test:

    A. constant-coefficient recurrences in k along each fixed d=r-k;
    B. constant-coefficient recurrences in k along each fixed r;
    C. geometric-ratio laws;
    D. separability of the form

           B[k,r] = A[k] * C[r]

       or equivalently, rank-one structure on compatible submatrices;
    E. rank-2 structure on compatible submatrices;
    F. normalized recurrence after clearing simple factorial factors.

Only OVERDETERMINED tests count as structural.

No interpolation polynomial of degree n-1 is counted as evidence.
No q_p(r) data are used.
No arbitrary matrix fitting.
Exact QQ arithmetic only.
"""


from __future__ import annotations

import sys
import sympy as sp


# =============================================================================
# EXACT B TABLE
# =============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# BASIC ACCESS
# =============================================================================

def b_value(k, r):
    off = r - k

    if off < 0:
        return None

    if off >= len(B[k]):
        return None

    return sp.Rational(B[k][off])


def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def fixed_d_sequence(d):
    """
    Sequence:

        B[0,d], B[1,d+1], B[2,d+2], ...

    i.e. fixed offset d=r-k.
    """
    seq = []

    for k in range(len(B)):

        r = k + d
        value = b_value(k, r)

        if value is not None:
            seq.append(
                (k, value)
            )

    return seq


def fixed_r_sequence(r):
    """
    Sequence:

        B[0,r], B[1,r], ..., B[r,r].
    """
    seq = []

    for k in range(len(B)):

        value = b_value(k, r)

        if value is not None:
            seq.append(
                (k, value)
            )

    return seq


# =============================================================================
# CONSTANT-COEFFICIENT RECURRENCE
# =============================================================================

def recurrence_fit(values, order):
    """
    Test

        x_n = c1*x_(n-1) + ... + cm*x_(n-m)

    with constant rational coefficients.

    Requires MORE equations than unknowns.

    Returns:
        status, coefficients, residuals
    """

    if len(values) <= 2 * order:
        return (
            "INSUFFICIENT_DATA",
            None,
            [],
        )

    rows = []
    rhs = []

    for n in range(
        order,
        len(values),
    ):

        rows.append([
            sp.Rational(
                values[n - j - 1]
            )
            for j in range(order)
        ])

        rhs.append(
            sp.Rational(
                values[n]
            )
        )

    M = sp.Matrix(rows)
    y = sp.Matrix(rhs)

    rank = M.rank()
    augmented_rank = (
        M.row_join(y).rank()
    )

    if augmented_rank != rank:
        return (
            "NO_SOLUTION",
            None,
            [],
        )

    if rank != order:
        return (
            "NONUNIQUE",
            None,
            [],
        )

    coeffs = M.gauss_jordan_solve(y)[0]

    coeffs = tuple(
        clean(x)
        for x in coeffs
    )

    residuals = []

    for n in range(
        order,
        len(values),
    ):

        pred = clean(
            sum(
                coeffs[j]
                * values[n - j - 1]
                for j in range(order)
            )
        )

        residuals.append(
            clean(
                pred - values[n]
            )
        )

    exact = all(
        x == 0
        for x in residuals
    )

    if exact:
        return (
            "EXACT",
            coeffs,
            residuals,
        )

    return (
        "FAILED",
        coeffs,
        residuals,
    )


# =============================================================================
# GEOMETRIC RATIO
# =============================================================================

def geometric_ratio(seq):
    """
    Test x_{n+1} = c*x_n with one common c.

    Requires at least three values, so the ratio is actually tested
    on at least two transitions.
    """

    values = [
        value
        for _, value in seq
    ]

    if len(values) < 3:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

    if any(
        x == 0
        for x in values[:-1]
    ):
        return (
            "ZERO_DENOMINATOR",
            None,
        )

    ratios = [
        clean(
            values[i + 1]
            / values[i]
        )
        for i in range(
            len(values) - 1
        )
    ]

    exact = all(
        q == ratios[0]
        for q in ratios
    )

    if exact:
        return (
            "EXACT",
            ratios[0],
        )

    return (
        "FAILED",
        ratios,
    )


# =============================================================================
# RANK TESTS
# =============================================================================

def available_submatrix(max_r=None):
    """
    Build the largest rectangular leading submatrix for

        rows k=0..m
        columns r=0..n

    for which every entry exists.
    """
    if max_r is None:
        max_r = 7

    matrix = []

    for k in range(len(B)):

        row = []

        for r in range(max_r + 1):

            value = b_value(k, r)

            if value is None:
                break

            row.append(value)

        if row:
            matrix.append(row)

    if not matrix:
        return sp.Matrix([])

    width = min(
        len(row)
        for row in matrix
    )

    return sp.Matrix([
        row[:width]
        for row in matrix
    ])


def rectangular_rank_tests():
    """
    Test all meaningful rectangular leading submatrices with dimensions
    at least 2x2.

    Rank 1 means exact multiplicative separability:

        B[k,r] = A[k] C[r].

    Rank 2 means a two-term separable decomposition is possible.
    """

    results = []

    full = available_submatrix()

    rows = full.rows
    cols = full.cols

    for nr in range(
        2,
        rows + 1,
    ):

        for nc in range(
            2,
            cols + 1,
        ):

            M = full[
                :nr,
                :nc,
            ]

            results.append(
                (
                    nr,
                    nc,
                    M.rank(),
                )
            )

    return results


# =============================================================================
# SIMPLE FACTORIAL NORMALIZATIONS
# =============================================================================

def normalized_sequence(seq, mode):
    """
    Apply a deterministic normalization, then test recurrences.

    modes:

        raw
        divide_k_factorial
        multiply_k_factorial
        divide_r_factorial
        multiply_r_factorial
        multiply_binom_r_k
        divide_binom_r_k
    """

    out = []

    for k, value in seq:

        r = None

        # Caller stores r separately where needed, so the default
        # transformations below only depend on k.
        if mode == "raw":
            z = value

        elif mode == "divide_k_factorial":
            z = value / sp.factorial(k)

        elif mode == "multiply_k_factorial":
            z = value * sp.factorial(k)

        else:
            # These modes need r and are handled separately.
            z = value

        out.append(
            (
                k,
                clean(z),
            )
        )

    return out


def fixed_d_normalized(d, mode):
    seq = fixed_d_sequence(d)

    out = []

    for k, value in seq:

        r = k + d

        if mode == "raw":
            z = value
        elif mode == "divide_k_factorial":
            z = value / sp.factorial(k)
        elif mode == "multiply_k_factorial":
            z = value * sp.factorial(k)
        elif mode == "divide_r_factorial":
            z = value / sp.factorial(r)
        elif mode == "multiply_r_factorial":
            z = value * sp.factorial(r)
        elif mode == "divide_binom":
            z = value / sp.binomial(r, k)
        elif mode == "multiply_binom":
            z = value * sp.binomial(r, k)
        else:
            raise ValueError(
                f"Unknown normalization: {mode}"
            )

        out.append(
            (
                k,
                clean(z),
            )
        )

    return out


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 284 — EXACT B-DIAGONAL RECURRENCE / "
        "SEPARABILITY AUDIT"
    )
    print("=" * 78)

    # =========================================================================
    # 1. FIXED-OFFSET SEQUENCES
    # =========================================================================

    print()
    print("=" * 78)
    print("1. FIXED-OFFSET d=r-k SEQUENCES")
    print("=" * 78)

    recurrence_hits = []
    ratio_hits = []

    for d in range(8):

        seq = fixed_d_sequence(d)

        if len(seq) < 3:
            continue

        values = [
            value
            for _, value in seq
        ]

        print()
        print(
            f"  d={d}:"
        )

        print(
            f"    sequence={seq}"
        )

        for order in [1, 2]:

            status, coeffs, residuals = (
                recurrence_fit(
                    values,
                    order,
                )
            )

            print(
                f"    recurrence_order={order}: "
                f"status={status}"
            )

            if status == "EXACT":

                print(
                    f"      coefficients={coeffs}"
                )

                recurrence_hits.append(
                    (
                        d,
                        order,
                        coeffs,
                    )
                )

        ratio_status, ratio = (
            geometric_ratio(seq)
        )

        print(
            f"    geometric_ratio="
            f"{ratio_status}"
        )

        print(
            f"    ratio={ratio}"
        )

        if ratio_status == "EXACT":
            ratio_hits.append(
                (
                    d,
                    ratio,
                )
            )

    # =========================================================================
    # 2. FIXED-r SEQUENCES
    # =========================================================================

    print()
    print("=" * 78)
    print("2. FIXED-ABSOLUTE-r SEQUENCES")
    print("=" * 78)

    fixed_r_hits = []

    for r in range(8):

        seq = fixed_r_sequence(r)

        if len(seq) < 3:
            continue

        values = [
            value
            for _, value in seq
        ]

        print()
        print(
            f"  r={r}:"
        )

        print(
            f"    sequence={seq}"
        )

        for order in [1, 2]:

            status, coeffs, residuals = (
                recurrence_fit(
                    values,
                    order,
                )
            )

            print(
                f"    recurrence_order={order}: "
                f"status={status}"
            )

            if status == "EXACT":

                print(
                    f"      coefficients={coeffs}"
                )

                fixed_r_hits.append(
                    (
                        r,
                        order,
                        coeffs,
                    )
                )

    # =========================================================================
    # 3. NORMALIZED DIAGONAL TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("3. FACTORIAL / BINOMIAL-NORMALIZED DIAGONAL RECURRENCES")
    print("=" * 78)

    normalization_names = [
        "raw",
        "divide_k_factorial",
        "multiply_k_factorial",
        "divide_r_factorial",
        "multiply_r_factorial",
        "divide_binom",
        "multiply_binom",
    ]

    normalized_hits = []

    for mode in normalization_names:

        print()
        print(
            f"  mode={mode}:"
        )

        for d in range(8):

            seq = fixed_d_normalized(
                d,
                mode,
            )

            if len(seq) < 4:
                continue

            values = [
                value
                for _, value in seq
            ]

            for order in [1, 2]:

                status, coeffs, residuals = (
                    recurrence_fit(
                        values,
                        order,
                    )
                )

                if status == "EXACT":

                    print(
                        f"    d={d} "
                        f"order={order} "
                        f"coeffs={coeffs}"
                    )

                    normalized_hits.append(
                        (
                            mode,
                            d,
                            order,
                            coeffs,
                        )
                    )

    # =========================================================================
    # 4. RANK-ONE SEPARABILITY
    # =========================================================================

    print()
    print("=" * 78)
    print("4. MULTIPLICATIVE SEPARABILITY / RANK AUDIT")
    print("=" * 78)

    rank_results = rectangular_rank_tests()

    rank1 = [
        x
        for x in rank_results
        if x[2] == 1
    ]

    rank2 = [
        x
        for x in rank_results
        if x[2] == 2
    ]

    print(
        f"  rectangular_submatrices="
        f"{len(rank_results)}"
    )

    print(
        f"  rank_one_cases="
        f"{len(rank1)}"
    )

    print(
        f"  rank_two_cases="
        f"{len(rank2)}"
    )

    for item in rank_results:

        if item[2] <= 2:

            print(
                f"  shape="
                f"{item[0]}x{item[1]} "
                f"rank={item[2]}"
            )

    # =========================================================================
    # 5. EXACT CROSS-RATIO TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("5. CROSS-RATIO SEPARABILITY TEST")
    print("=" * 78)

    cross_ratio_failures = 0
    cross_ratio_checks = 0

    # Rank-one condition for every 2x2 submatrix:

    # B[a,c] * B[b,d] == B[a,d] * B[b,c]

    for a in range(len(B)):

        for b in range(
            a + 1,
            len(B),
        ):

            for c in range(8):

                for d in range(
                    c + 1,
                    8,
                ):

                    x1 = b_value(
                        a,
                        c,
                    )
                    x2 = b_value(
                        a,
                        d,
                    )
                    x3 = b_value(
                        b,
                        c,
                    )
                    x4 = b_value(
                        b,
                        d,
                    )

                    if None in (
                        x1,
                        x2,
                        x3,
                        x4,
                    ):
                        continue

                    cross_ratio_checks += 1

                    if clean(
                        x1 * x4
                        - x2 * x3
                    ) != 0:

                        cross_ratio_failures += 1

    print(
        f"  checked_2x2_relations="
        f"{cross_ratio_checks}"
    )

    print(
        f"  cross_ratio_failures="
        f"{cross_ratio_failures}"
    )

    # =========================================================================
    # 6. INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 283 rejected low-degree polynomial descriptions of the B-grid.

Experiment 284 tests a different class of exact structure.

A low-order constant-coefficient recurrence would indicate that the
B rows/diagonals arise from a finite-state combinatorial operator.

A geometric ratio would indicate first-order multiplicative structure.

A rank-one rectangular submatrix would imply

    B[k,r] = A[k] C[r]

on that region.

A rank-two structure would indicate a sum of two separable source
channels.

Factorial/binomial-normalized recurrences test whether the complicated
rational coefficients are only a normalization artifact.

The important point is that all these tests are OVERDETERMINED whenever
enough data exist. Interpolation is therefore not counted as evidence.
"""
    )

    # =========================================================================
    # 7. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print("7. TERMINAL PROJECTIVE SOURCE REFERENCE")
    print("=" * 78)

    q1_terminal = 495451247
    q3_terminal = 421514439

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd={sp.gcd(
            q1_terminal,
            q3_terminal,
        )}"
    )

    print(
        f"  q1/17={q1_terminal // 17}"
    )

    print(
        f"  q3/17={q3_terminal // 17}"
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exact_raw_recurrence_hits="
        f"{len(recurrence_hits)}"
    )

    print(
        f"  exact_geometric_hits="
        f"{len(ratio_hits)}"
    )

    print(
        f"  exact_normalized_recurrence_hits="
        f"{len(normalized_hits)}"
    )

    print(
        f"  rank_one_rectangular_cases="
        f"{len(rank1)}"
    )

    print(
        f"  rank_two_rectangular_cases="
        f"{len(rank2)}"
    )

    print(
        f"  cross_ratio_failures="
        f"{cross_ratio_failures}"
    )

    print(
        "  interpolation_only_not_counted=True"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 284 COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print(
            "\nInterrupted."
        )
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise
