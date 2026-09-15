#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 291RR — EXACT FACTORIAL-NORMALIZED B-GRID / INTEGER-BASIS AUDIT
==============================================================================

Purpose:

    Experiment 290R found that

        C[k,r] = r! * B[k,r]

    is integral for every supplied B-entry.

Experiment 291RR studies that integer grid using:

    * row/column content;
    * finite differences in r;
    * finite differences in k;
    * shifted diagonals d=r-k;
    * binomial normalization;
    * falling-factorial normalization.

No q-family data.
No fitted matrix.
No floating point.
No interpolation is counted as proof.
"""

from __future__ import annotations

import math
import sys
import sympy as sp


# ============================================================================
# EXACT B DATA
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
# BASIC HELPERS
# ============================================================================

def clean(value):
    """Return an exact simplified SymPy value."""
    if value is None:
        return None

    value = sp.sympify(value)

    return sp.factor(
        sp.cancel(
            sp.expand(value)
        )
    )


def falling(n, m):
    """Exact falling factorial n_(m)."""
    if m < 0:
        return sp.Integer(0)

    result = sp.Integer(1)

    for j in range(m):
        result *= n - j

    return clean(result)


def integer_grid():
    """
    Construct

        C[k,r] = r! * B[k,r].

    Every entry is asserted to be integral.
    """

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
                    f"r!B[{k},{r}] is not integral: {value}"
                )

            C[k][r] = sp.Integer(value)

    return C


def gcd_values(values):
    g = 0

    for value in values:
        g = math.gcd(
            g,
            abs(int(value)),
        )

    return g


def forward_differences(values):
    """
    Return all forward-difference levels.
    """
    levels = [list(values)]

    while len(levels[-1]) > 1:

        previous = levels[-1]

        current = []

        for i in range(
            len(previous) - 1
        ):
            current.append(
                previous[i + 1]
                - previous[i]
            )

        levels.append(
            current
        )

    return levels


# ============================================================================
# 1. INTEGER GRID
# ============================================================================

def show_integer_grid(C):

    print()
    print("=" * 78)
    print("1. FACTORIAL-NORMALIZED INTEGER GRID")
    print("=" * 78)

    for k in sorted(C):

        row = [
            C[k][r]
            for r in sorted(C[k])
        ]

        print()
        print(f"  k={k}:")
        print(f"    r!B={row}")
        print(
            f"    row_gcd={gcd_values(row)}"
        )


# ============================================================================
# 2. COLUMN CONTENT
# ============================================================================

def show_column_content(C):

    print()
    print("=" * 78)
    print("2. INTEGER COLUMN CONTENT")
    print("=" * 78)

    for r in range(8):

        values = []

        for k in sorted(C):

            if r in C[k]:
                values.append(C[k][r])

        if not values:
            continue

        print()
        print(f"  r={r}:")
        print(f"    values={values}")
        print(
            f"    gcd={gcd_values(values)}"
        )


# ============================================================================
# 3. BINOMIAL NORMALIZATIONS
# ============================================================================

def show_binomial_normalizations(C):

    print()
    print("=" * 78)
    print("3. BINOMIAL NORMALIZATIONS")
    print("=" * 78)

    modes = [
        "divide_C7r",
        "divide_Crk",
        "divide_Crd",
        "divide_C7kd",
    ]

    for mode in modes:

        print()
        print(f"  mode={mode}")

        for k in sorted(C):

            values = []

            for r in sorted(C[k]):

                d = r - k
                value = C[k][r]

                if mode == "divide_C7r":
                    factor = sp.binomial(
                        7,
                        r,
                    )

                elif mode == "divide_Crk":
                    factor = sp.binomial(
                        r,
                        k,
                    )

                elif mode == "divide_Crd":
                    factor = sp.binomial(
                        r,
                        d,
                    )

                else:
                    factor = sp.binomial(
                        7 - k,
                        d,
                    )

                if factor == 0:
                    values.append(None)
                else:
                    values.append(
                        clean(
                            sp.Rational(
                                value,
                                factor,
                            )
                        )
                    )

            print(
                f"    k={k}: {values}"
            )


# ============================================================================
# 4. FALLING-FACTORIAL NORMALIZATIONS
# ============================================================================

def show_falling_normalizations(C):

    print()
    print("=" * 78)
    print("4. FALLING-FACTORIAL NORMALIZATIONS")
    print("=" * 78)

    modes = [
        "divide_rfalling_k",
        "divide_rfalling_d",
        "divide_7k_falling_d",
    ]

    for mode in modes:

        print()
        print(f"  mode={mode}")

        for k in sorted(C):

            values = []

            for r in sorted(C[k]):

                d = r - k
                value = C[k][r]

                if mode == "divide_rfalling_k":
                    factor = falling(
                        r,
                        k,
                    )

                elif mode == "divide_rfalling_d":
                    factor = falling(
                        r,
                        d,
                    )

                else:
                    factor = falling(
                        7 - k,
                        d,
                    )

                if factor == 0:
                    values.append(None)
                else:
                    values.append(
                        clean(
                            sp.Rational(
                                value,
                                factor,
                            )
                        )
                    )

            print(
                f"    k={k}: {values}"
            )


# ============================================================================
# 5. DIFFERENCES IN r
# ============================================================================

def audit_r_differences(C):

    print()
    print("=" * 78)
    print("5. FINITE DIFFERENCES IN r")
    print("=" * 78)

    for k in sorted(C):

        values = [
            C[k][r]
            for r in sorted(C[k])
        ]

        print()
        print(f"  k={k}")

        levels = forward_differences(
            values
        )

        for order, row in enumerate(levels):

            print(
                f"    order={order}: {row}"
            )


# ============================================================================
# 6. DIFFERENCES IN k FOR FIXED r
# ============================================================================

def audit_k_differences(C):

    print()
    print("=" * 78)
    print("6. FINITE DIFFERENCES IN k")
    print("=" * 78)

    for r in range(8):

        values = []

        for k in sorted(C):

            if r in C[k]:
                values.append(
                    C[k][r]
                )

        if len(values) < 2:
            continue

        print()
        print(f"  r={r}")

        levels = forward_differences(
            values
        )

        for order, row in enumerate(levels):

            print(
                f"    order={order}: {row}"
            )


# ============================================================================
# 7. DIAGONAL d = r-k
# ============================================================================

def audit_shifted_diagonals(C):

    print()
    print("=" * 78)
    print("7. SHIFTED-DIAGONAL d=r-k AUDIT")
    print("=" * 78)

    diagonals = {}

    for k in sorted(C):

        for r, value in C[k].items():

            d = r - k

            diagonals.setdefault(
                d,
                [],
            ).append(
                (
                    k,
                    value,
                )
            )

    for d in sorted(diagonals):

        sequence = diagonals[d]

        print()
        print(f"  d={d}:")
        print(
            f"    sequence={sequence}"
        )

        values = [
            value
            for _, value in sequence
        ]

        levels = forward_differences(
            values
        )

        for order, row in enumerate(levels):

            print(
                f"    diff_order={order}: "
                f"{row}"
            )


# ============================================================================
# 8. EXACT INTEGER CONTENT OF DIFFERENCE LEVELS
# ============================================================================

def audit_difference_content(C):

    print()
    print("=" * 78)
    print(
        "8. DIFFERENCE CONTENT / STABILIZATION AUDIT"
    )
    print("=" * 78)

    for k in sorted(C):

        values = [
            C[k][r]
            for r in sorted(C[k])
        ]

        levels = forward_differences(
            values
        )

        print()
        print(f"  k={k}:")

        for order, row in enumerate(levels):

            if not row:
                continue

            g = gcd_values(row)

            print(
                f"    order={order}: "
                f"gcd={g}"
            )


# ============================================================================
# 9. STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 290R found the only successful elementary normalization:

    C[k,r] = r! B[k,r] is integral everywhere.

Experiment 291RR therefore works entirely with this integer grid.

The key tests are now:

    * whether finite differences in r stabilize;
    * whether finite differences in k stabilize;
    * whether shifted diagonals d=r-k have simpler arithmetic;
    * whether binomial/falling-factorial divisions reduce the complexity.

A constant final difference would indicate a genuine polynomial law in
the relevant discrete coordinate.

Repeated divisibility or a stable gcd at a fixed difference order may
indicate a hidden combinatorial normalization.

Interpolation is deliberately not used as evidence.
Only structural behavior visible in repeated difference levels counts.
"""
    )


# ============================================================================
# 10. MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 291RR — EXACT FACTORIAL-NORMALIZED "
        "B-GRID / INTEGER-BASIS AUDIT"
    )
    print("=" * 78)

    C = integer_grid()

    show_integer_grid(C)

    show_column_content(C)

    show_binomial_normalizations(C)

    show_falling_normalizations(C)

    audit_r_differences(C)

    audit_k_differences(C)

    audit_shifted_diagonals(C)

    audit_difference_content(C)

    interpretation()

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  r_factorial_normalization_exact=True"
    )

    print(
        "  integer_grid_reconstructed=True"
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
        "EXPERIMENT 291RR COMPLETE"
    )


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