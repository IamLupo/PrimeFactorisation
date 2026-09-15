#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 158 — EXACT MINIMAL-COFACTOR COLUMN-LAW / PRIMITIVE-VECTOR AUDIT
==============================================================================

Experiment 157 established:

    gcd(u_j) = 1701
    Delta_11 / 2^24 = 567
    gcd(R_minimal_row) = 3

After removing the common factor 3 from the normalized minimal row,
the resulting six-entry vector has gcd 1 and rapidly separating
3-adic and 7-adic residues.

Experiment 158 asks whether the primitive six-entry vector has a
simple column-index structure.

This is deliberately restricted to low-complexity tests.

The selected columns are

    J = [3,5,7,9,10,11].

We define

    W_j = R[5,j] / 3.

Tests:

    1. exact vector and primitive content;
    2. finite differences with respect to the ACTUAL column indices;
    3. second and third divided-difference consistency for low-degree
       polynomial behavior;
    4. simple reflection tests j -> 14-j;
    5. parity partition;
    6. residue patterns modulo 2,3,5,7;
    7. low-degree polynomial consistency modulo small primes;
    8. exact affine / quadratic relation tests over QQ.

The polynomial tests are intentionally low degree.

A degree-5 polynomial through six points is meaningless here and
is therefore NOT considered.

No recurrence is fitted.
No arbitrary high-degree interpolation.
No SymPy.
No floating point.
No extrapolation.
No connection to the original (p,q)-kernel is asserted.
"""

from __future__ import annotations

from fractions import Fraction
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
ROW = 5
MIN_V2 = 24

J = [3, 5, 7, 9, 10, 11]


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


def gcd_list(values):
    from math import gcd
    from functools import reduce

    vals = [
        abs(int(x))
        for x in values
        if int(x) != 0
    ]

    if not vals:
        return 0

    return reduce(gcd, vals)


def valuation(n, p):
    n = abs(int(n))

    if n == 0:
        return float("inf")

    out = 0

    while n % p == 0:
        n //= p
        out += 1

    return out


# ============================================================================
# BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    M = [[int(x) for x in row] for row in A]
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
            M[k], M[pivot_row] = (
                M[pivot_row],
                M[k],
            )
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
                            "Bareiss exact division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD SYSTEM
# ============================================================================

def build_system():

    M = []

    for p, pnext in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# F-CORE
# ============================================================================

def build_core(M):

    b0 = M[5]
    b1 = M[10]

    H = [
        [b0[12], b0[13]],
        [b1[12], b1[13]],
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_h == 0:
        raise ArithmeticError(
            "Boundary block singular."
        )

    H_inv = [
        [
            Fraction(H[1][1], det_h),
            Fraction(-H[0][1], det_h),
        ],
        [
            Fraction(-H[1][0], det_h),
            Fraction(H[0][0], det_h),
        ],
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * b0[j]
                + H_inv[i][1] * b1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        out = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral F-core entry."
                )

            out.append(value.numerator)

        core.append(out)

    return core


# ============================================================================
# COFACTORS
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


def cofactor_matrix(A):

    n = len(A)
    C = []

    for i in range(n):

        row = []

        for j in range(n):

            value = bareiss_det(
                maximal_minor(A, i, j)
            )

            if (i + j) & 1:
                value = -value

            row.append(value)

        C.append(row)

    return C


# ============================================================================
# EXACT LOW-DEGREE POLYNOMIAL TESTS
# ============================================================================

def solve_interpolation(points, degree):

    """
    Solve for a polynomial of the requested degree using exactly degree+1
    points, then test it against every remaining point.

    This is only used for degree <= 2.
    """

    if degree == 0:
        basis_size = 1
    elif degree == 1:
        basis_size = 2
    elif degree == 2:
        basis_size = 3
    else:
        raise ValueError(
            "Only degrees 0,1,2 are allowed."
        )

    selected = points[:basis_size]

    A = []

    b = []

    for x, y in selected:

        A.append(
            [
                Fraction(x) ** k
                for k in range(basis_size)
            ]
        )

        b.append(Fraction(y))

    # Gaussian elimination over QQ.
    n = basis_size

    for col in range(n):

        pivot = None

        for r in range(col, n):

            if A[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            return None

        if pivot != col:

            A[col], A[pivot] = (
                A[pivot],
                A[col],
            )

            b[col], b[pivot] = (
                b[pivot],
                b[col],
            )

        pivot_value = A[col][col]

        for j in range(col, n):
            A[col][j] /= pivot_value

        b[col] /= pivot_value

        for r in range(n):

            if r == col:
                continue

            factor = A[r][col]

            if factor == 0:
                continue

            for j in range(col, n):
                A[r][j] -= factor * A[col][j]

            b[r] -= factor * b[col]

    coeffs = b

    return coeffs


def polynomial_value(coeffs, x):

    value = Fraction(0)

    for k, c in enumerate(coeffs):
        value += c * (Fraction(x) ** k)

    return value


def test_degree(points, degree):

    coeffs = solve_interpolation(
        points,
        degree,
    )

    if coeffs is None:
        return False, None, []

    failures = []

    for x, y in points:

        predicted = polynomial_value(
            coeffs,
            x,
        )

        if predicted != y:
            failures.append(
                (
                    x,
                    y,
                    predicted,
                )
            )

    return (
        len(failures) == 0,
        coeffs,
        failures,
    )


# ============================================================================
# MODULAR LOW-DEGREE TEST
# ============================================================================

def mod_poly_fit(points, degree, prime):

    m = degree + 1

    selected = points[:m]

    A = []

    b = []

    for x, y in selected:

        A.append(
            [
                pow(x, k, prime)
                for k in range(m)
            ]
        )

        b.append(
            y % prime
        )

    # Gaussian elimination mod prime.
    for col in range(m):

        pivot = None

        for r in range(col, m):

            if A[r][col] % prime != 0:
                pivot = r
                break

        if pivot is None:
            return False, None

        A[col], A[pivot] = (
            A[pivot],
            A[col],
        )

        b[col], b[pivot] = (
            b[pivot],
            b[col],
        )

        inv = pow(
            A[col][col] % prime,
            -1,
            prime,
        )

        for j in range(col, m):
            A[col][j] = (
                A[col][j] * inv
            ) % prime

        b[col] = (
            b[col] * inv
        ) % prime

        for r in range(m):

            if r == col:
                continue

            factor = A[r][col]

            if factor == 0:
                continue

            for j in range(col, m):
                A[r][j] = (
                    A[r][j]
                    - factor * A[col][j]
                ) % prime

            b[r] = (
                b[r]
                - factor * b[col]
            ) % prime

    coeffs = b

    for x, y in points:

        predicted = 0

        for k, c in enumerate(coeffs):
            predicted += (
                c
                * pow(x, k, prime)
            )

        if predicted % prime != y % prime:
            return False, coeffs

    return True, coeffs


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 158 — EXACT MINIMAL-COFACTOR "
        "COLUMN-LAW / PRIMITIVE-VECTOR AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)
    C = cofactor_matrix(A)

    raw = {
        j: C[ROW][j]
        for j in J
    }

    # Delta-normalized row.
    R = {
        j: Fraction(
            C[ROW][j],
            DELTA_11,
        )
        for j in J
    }

    primitive = {
        j: R[j] / 3
        for j in J
    }

    values = [
        primitive[j]
        for j in J
    ]

    # ------------------------------------------------------------------
    # 1. VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT PRIMITIVE VECTOR")
    print("=" * 78)

    print(
        f"  columns={J}"
    )

    print(
        f"  primitive_values={values}"
    )

    print(
        f"  gcd={gcd_list([x.numerator for x in values])}"
    )

    # ------------------------------------------------------------------
    # 2. RAW SUCCESSIVE DIFFERENCES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SUCCESSIVE DIFFERENCES IN COLUMN ORDER")
    print("=" * 78)

    first_diff = [
        values[i + 1] - values[i]
        for i in range(len(values) - 1)
    ]

    second_diff = [
        first_diff[i + 1] - first_diff[i]
        for i in range(len(first_diff) - 1)
    ]

    third_diff = [
        second_diff[i + 1] - second_diff[i]
        for i in range(len(second_diff) - 1)
    ]

    print(
        f"  first_diff={first_diff}"
    )

    print(
        f"  second_diff={second_diff}"
    )

    print(
        f"  third_diff={third_diff}"
    )

    # ------------------------------------------------------------------
    # 3. AFFINE / QUADRATIC TESTS
    # ------------------------------------------------------------------

    points = [
        (J[i], values[i])
        for i in range(len(J))
    ]

    print()
    print("=" * 78)
    print("3. EXACT LOW-DEGREE POLYNOMIAL AUDIT")
    print("=" * 78)

    for degree in (0, 1, 2):

        ok, coeffs, failures = test_degree(
            points,
            degree,
        )

        print()
        print(
            f"  degree={degree}: exact={ok}"
        )

        if coeffs is not None:
            print(
                f"    coefficients={coeffs}"
            )

        if failures:
            print(
                f"    failure_count={len(failures)}"
            )

    # ------------------------------------------------------------------
    # 4. PARITY PARTITION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COLUMN PARITY PARTITION")
    print("=" * 78)

    even_columns = [
        j for j in J
        if j % 2 == 0
    ]

    odd_columns = [
        j for j in J
        if j % 2 == 1
    ]

    print(
        f"  even_columns={even_columns}"
    )

    print(
        f"  even_values="
        f"{[primitive[j] for j in even_columns]}"
    )

    print(
        f"  odd_columns={odd_columns}"
    )

    print(
        f"  odd_values="
        f"{[primitive[j] for j in odd_columns]}"
    )

    # ------------------------------------------------------------------
    # 5. REFLECTION AUDIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SIMPLE COLUMN REFLECTION AUDIT")
    print("=" * 78)

    reflections = [
        (3, 11),
        (5, 9),
    ]

    for a, b in reflections:

        va = primitive[a]
        vb = primitive[b]

        print(
            f"  ({a},{b}): "
            f"equal={va == vb} "
            f"opposite={va == -vb} "
            f"sum_zero={va + vb == 0}"
        )

    # ------------------------------------------------------------------
    # 6. MODULAR POLYNOMIAL TESTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. LOW-DEGREE MODULAR POLYNOMIAL AUDIT")
    print("=" * 78)

    integer_points = [
        (
            J[i],
            primitive[J[i]],
        )
        for i in range(len(J))
    ]

    for prime in (2, 3, 5, 7, 11):

        print()
        print(
            f"  prime={prime}"
        )

        for degree in (1, 2):

            ok, coeffs = mod_poly_fit(
                [
                    (
                        x,
                        y.numerator
                        // y.denominator
                    )
                    for x, y in integer_points
                ],
                degree,
                prime,
            )

            print(
                f"    degree={degree}: "
                f"exact_mod_prime={ok} "
                f"coeffs={coeffs}"
            )

    # ------------------------------------------------------------------
    # 7. 3-ADIC VALUATION PATTERN
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PRIMITIVE 3-ADIC / 7-ADIC PROFILE")
    print("=" * 78)

    for j in J:

        numerator = primitive[j].numerator
        denominator = primitive[j].denominator

        print(
            f"  j={j}: "
            f"v3_num={valuation(numerator,3)} "
            f"v7_num={valuation(numerator,7)} "
            f"v3_den={valuation(denominator,3)} "
            f"v7_den={valuation(denominator,7)}"
        )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 157 showed that removing the common factor 3 leaves a
primitive six-entry vector.

Experiment 158 deliberately asks a much weaker question than fitting
an arbitrary polynomial:

    does this primitive vector obey a genuinely low-complexity column
    law?

Only degree 0, 1, and 2 are tested over the exact rationals.

The actual column locations

    3,5,7,9,10,11

are retained rather than replacing them by artificial consecutive
indices.

Reflection and parity are tested separately.

Low-degree modular tests are independent of the huge rational size.

A negative result would support the interpretation that the minimal
cofactor layer is an arithmetic vector rather than a simple polynomial
function of the column index.

A positive low-degree identity would be much more interesting because
it would explain several large integers simultaneously without
introducing high-degree interpolation.

No conclusion is promoted beyond the exact finite tests.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(values) == 6
        and gcd_list(
            [x.numerator for x in values]
        ) == 1
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  six_entry_vector_exact={len(values) == 6}"
    )

    print(
        f"  primitive_content_one="
        f"{gcd_list([x.numerator for x in values]) == 1}"
    )

    print(
        "  low_degree_tests_completed=True"
    )

    print(
        "  modular_tests_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 158 COMPLETE")


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

