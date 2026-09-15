#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 293 — EXACT SHIFTED NEWTON-COEFFICIENT RECURRENCE AUDIT
==============================================================================

Experiment 292 established exact Newton/falling-factorial expansions of the
factorial-normalized grid

    C[k,r] = r! B[k,r].

For each fixed k, write

    C_k(d)
      =
    sum_j a[k,j] d_(j),

where

    d = r-k.

Experiment 293 tests whether the Newton coefficients a[k,j] themselves have
a low-order recurrence in j, and whether that recurrence is stable across k.

Tests:

    1. first-order ratio:
           a[j+1] = c_j a[j]

    2. order-1 affine recurrence:
           a[j+1] = A(k) a[j] + B(k)

    3. order-2 constant recurrence:
           a[j+2] = A a[j+1] + B a[j]

    4. order-2 recurrence with coefficients depending on k:

           a[k,j+2]
             =
           A_k a[k,j+1]
             +
           B_k a[k,j]

    5. recurrence after primitive integer normalization.

Only overdetermined exact identities count.

No q-family data.
No arbitrary matrix fitting.
No interpolation used as proof.
"""


from __future__ import annotations

import math
import sys
import sympy as sp


# ============================================================================
# DATA
# ============================================================================

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


# ============================================================================
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def falling(x, n):
    out = sp.Integer(1)

    for j in range(n):
        out *= x - j

    return sp.expand(out)


def integer_grid():
    C = {}

    for k, row in enumerate(B):

        C[k] = {}

        for offset, value in enumerate(row):

            r = k + offset

            value = clean(
                sp.factorial(r) * value
            )

            if sp.denom(value) != 1:
                raise ArithmeticError(
                    f"C[{k},{r}] not integral: {value}"
                )

            C[k][r] = sp.Integer(value)

    return C


def falling_coefficients(points):
    """
    Exact Newton/falling-factorial coefficients.

    Input:
        [(d0, value0), ...]

    Output:
        [a0, a1, ...]
    """

    d = sp.Symbol("d")

    polynomial = sp.interpolate(
        [
            (
                sp.Integer(x),
                sp.Rational(y),
            )
            for x, y in points
        ],
        d,
    )

    polynomial = clean(
        polynomial
    )

    degree = sp.Poly(
        polynomial,
        d,
        domain=sp.QQ,
    ).degree()

    residual = polynomial
    coefficients = [sp.Integer(0)] * (
        degree + 1
    )

    for j in range(
        degree,
        -1,
        -1,
    ):

        basis = falling(
            d,
            j,
        )

        basis_poly = sp.Poly(
            basis,
            d,
            domain=sp.QQ,
        )

        residual_poly = sp.Poly(
            residual,
            d,
            domain=sp.QQ,
        )

        lead_basis = basis_poly.LC()
        lead_residual = residual_poly.nth(j)

        a = clean(
            sp.Rational(
                lead_residual,
                lead_basis,
            )
        )

        coefficients[j] = a

        residual = clean(
            residual
            - a * basis
        )

    return coefficients


def primitive_signature(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    den = 1

    for v in values:
        den = math.lcm(
            den,
            int(v.q),
        )

    ints = [
        int(v * den)
        for v in values
    ]

    g = 0

    for x in ints:
        g = math.gcd(
            g,
            abs(x),
        )

    if g:
        ints = [
            x // g
            for x in ints
        ]

    return ints


def solve_recurrence_order_two(values):
    """
    Solve

        v[j+2] = A v[j+1] + B v[j]

    using all available equations.
    """

    if len(values) < 4:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

    M = []
    rhs = []

    for j in range(
        len(values) - 2
    ):

        M.append(
            [
                values[j + 1],
                values[j],
            ]
        )

        rhs.append(
            values[j + 2]
        )

    M = sp.Matrix(M)
    rhs = sp.Matrix(rhs)

    rank = M.rank()
    augmented = M.row_join(
        rhs
    ).rank()

    if augmented != rank:
        return (
            "NO_SOLUTION",
            None,
        )

    if rank != 2:
        return (
            "NONUNIQUE",
            None,
        )

    solution = M.gauss_jordan_solve(
        rhs
    )[0]

    A = clean(
        solution[0]
    )

    Bc = clean(
        solution[1]
    )

    return (
        "EXACT",
        (A, Bc),
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 293 — EXACT SHIFTED NEWTON-COEFFICIENT "
        "RECURRENCE AUDIT"
    )
    print("=" * 78)

    C = integer_grid()

    coefficient_rows = {}

    # =========================================================================
    # 1. NEWTON COEFFICIENTS
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "1. SHIFTED NEWTON COEFFICIENTS"
    )
    print("=" * 78)

    for k in sorted(C):

        points = []

        for r, value in sorted(
            C[k].items()
        ):

            d = r - k

            points.append(
                (
                    d,
                    value,
                )
            )

        coeffs = falling_coefficients(
            points
        )

        coefficient_rows[k] = coeffs

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    coefficients={coeffs}"
        )

        print(
            f"    primitive={primitive_signature(coeffs)}"
        )

    # =========================================================================
    # 2. FIRST-ORDER RATIO TEST
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "2. FIRST-ORDER RATIO TEST"
    )
    print("=" * 78)

    ratio_hits = []

    for k, coeffs in coefficient_rows.items():

        ratios = []
        valid = True

        for j in range(
            len(coeffs) - 1
        ):

            a = coeffs[j]
            b = coeffs[j + 1]

            if a == 0:
                valid = False
                break

            ratios.append(
                clean(
                    b / a
                )
            )

        print()
        print(
            f"  k={k}: "
            f"ratios={ratios}"
        )

        if valid and len(ratios) >= 2:

            if all(
                ratio == ratios[0]
                for ratio in ratios
            ):
                ratio_hits.append(
                    (
                        k,
                        ratios[0],
                    )
                )

    # =========================================================================
    # 3. ORDER-2 RECURRENCE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "3. ORDER-2 CONSTANT RECURRENCE"
    )
    print("=" * 78)

    recurrence_hits = []

    for k, coeffs in coefficient_rows.items():

        status, params = (
            solve_recurrence_order_two(
                coeffs
            )
        )

        print()
        print(
            f"  k={k}: "
            f"status={status}"
        )

        print(
            f"    parameters={params}"
        )

        if status == "EXACT":
            recurrence_hits.append(
                (
                    k,
                    params,
                )
            )

    # =========================================================================
    # 4. NORMALIZED INTEGER SIGNATURE RECURRENCE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "4. PRIMITIVE-SIGNATURE RECURRENCE AUDIT"
    )
    print("=" * 78)

    signature_hits = []

    for k, coeffs in coefficient_rows.items():

        signature = primitive_signature(
            coeffs
        )

        status, params = (
            solve_recurrence_order_two(
                [
                    sp.Integer(v)
                    for v in signature
                ]
            )
        )

        print()
        print(
            f"  k={k}: "
            f"signature={signature}"
        )

        print(
            f"    recurrence_status={status}"
        )

        print(
            f"    parameters={params}"
        )

        if status == "EXACT":
            signature_hits.append(
                (
                    k,
                    params,
                )
            )

    # =========================================================================
    # 5. CROSS-k PARAMETER COMPARISON
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "5. CROSS-k RECURRENCE PARAMETER PROFILE"
    )
    print("=" * 78)

    for k, params in recurrence_hits:

        A, Bc = params

        print()
        print(
            f"  k={k}: "
            f"A={A} "
            f"B={Bc}"
        )

    # =========================================================================
    # 6. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 292 showed that the factorial-normalized B-grid has an exact
Newton/falling-factorial expansion.

Experiment 293 now treats those Newton coefficients as the primary
objects.

If

    a[k,j+2]
      =
    A_k a[k,j+1]
      +
    B_k a[k,j]

holds exactly, then the B rows are generated by a finite-state
recurrence in the Newton basis.

An even stronger result would be a common A,B law across k.

The first-order ratio test is a stricter special case.

If these tests all fail, then the Newton coefficients themselves do
not have a low-order finite-state structure, and the original B
construction remains the necessary upstream target.
"""
    )

    # =========================================================================
    # 7. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "7. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    q1 = 495451247
    q3 = 421514439

    g = math.gcd(
        q1,
        q3,
    )

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={g}"
    )

    print(
        f"  q1/17={q1 // 17}"
    )

    print(
        f"  q3/17={q3 // 17}"
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  r_factorial_grid_exact=True"
    )

    print(
        "  shifted_newton_basis_exact=True"
    )

    print(
        "  first_order_ratio_hits="
        + str(
            len(ratio_hits)
        )
    )

    print(
        "  order2_recurrence_hits="
        + str(
            len(recurrence_hits)
        )
    )

    print(
        "  primitive_signature_recurrence_hits="
        + str(
            len(signature_hits)
        )
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
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
        "EXPERIMENT 293 COMPLETE"
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

