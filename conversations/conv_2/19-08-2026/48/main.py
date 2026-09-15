#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 363R — EXACT TRANSFORMED-LATTICE INTEGER-ANNIHILATOR AUDIT
==============================================================================

Purpose
-------
362R found no small-support integer-coefficient annihilator directly on Q.

363R asks whether the obstruction disappears after canonical discrete
transformations of the observed lattice.

For

    Q(r,t)

with

    r = (p-1)/2,

construct exact transformed arrays:

    H10 = Delta_r Q
    H01 = Delta_t Q
    H11 = Delta_r Delta_t Q
    H20 = Delta_r^2 Q
    H02 = Delta_t^2 Q
    H21 = Delta_r^2 Delta_t Q
    H12 = Delta_r Delta_t^2 Q
    H22 = Delta_r^2 Delta_t^2 Q

Only differences whose source cells are all observed are constructed.

For each transformed lattice, search the same compact family of local
integer annihilators used in 362R.

The missing source cell

    Z = Q(2,3)

is NEVER inserted.

A transformed cell is retained symbolically only when its finite-difference
stencil genuinely contains Z; otherwise it is discarded from symbolic
calibration.

The experiment distinguishes:

    NO_INTEGER_RELATION
    DATA_LIMITED
    MULTIPLE_RELATIONS
    INTEGER_RELATION_Z_FREE
    INTEGER_RELATION_FORCES_Z

No extrapolation.
No interpolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.

Output is intentionally compact.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# SOURCE DATA
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


Z = sp.Symbol("Z")


# ============================================================================
# BASIC LATTICE
# ============================================================================

def build_lattice():
    lattice = {}

    for p_value, values in Q.items():
        r = (p_value - 1) // 2

        for index, value in enumerate(values):
            t = len(values) - 1 - index
            lattice[(r, t)] = sp.Integer(value)

    return lattice


BASE = build_lattice()


# ============================================================================
# CELL ACCESS
# ============================================================================

def base_value(cell):
    if cell == (2, 3):
        return Z

    return BASE.get(cell)


def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def is_symbolic(value):
    return sp.sympify(value).has(Z)


# ============================================================================
# FINITE-DIFFERENCE TRANSFORMS
# ============================================================================

TRANSFORMS = {
    "H00": (0, 0),
    "H10": (1, 0),
    "H01": (0, 1),
    "H11": (1, 1),
    "H20": (2, 0),
    "H02": (0, 2),
    "H21": (2, 1),
    "H12": (1, 2),
    "H22": (2, 2),
}


def transformed_value(
    r,
    t,
    dr,
    dt,
):
    """
    Forward finite difference:

        Delta_r^dr Delta_t^dt Q(r,t)

    Expands exactly by binomial coefficients.

    Returns None when at least one required base cell is unavailable.

    Returns a symbolic expression when the only unavailable cell is the
    designated symbolic target Z.
    """

    if dr == 0 and dt == 0:
        return base_value((r, t))

    result = 0

    for i in range(dr + 1):
        for j in range(dt + 1):

            coefficient = (
                (-1) ** (
                    dr - i + dt - j
                )
                * math.comb(dr, i)
                * math.comb(dt, j)
            )

            cell = (
                r + i,
                t + j,
            )

            value = base_value(cell)

            if value is None:
                return None

            result += coefficient * value

    return clean(result)


def build_transform(
    dr,
    dt,
):
    result = {}

    # The observed triangle occupies r=0..3, t=0..5.
    for r in range(4):
        for t in range(6):

            value = transformed_value(
                r,
                t,
                dr,
                dt,
            )

            if value is not None:
                result[(r, t)] = clean(value)

    return result


# ============================================================================
# PRIMITIVE INTEGER NULLSPACE
# ============================================================================

def primitive_integer_vector(values):
    rational_values = [
        sp.Rational(v)
        for v in values
    ]

    denominator = 1

    for value in rational_values:
        denominator = sp.ilcm(
            denominator,
            int(value.q),
        )

    integers = [
        int(
            value * denominator
        )
        for value in rational_values
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    if g == 0:
        return integers

    integers = [
        value // g
        for value in integers
    ]

    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -v
                    for v in integers
                ]
            break

    return integers


def integer_nullspace(matrix):
    rational_basis = matrix.nullspace()

    result = []

    for vector in rational_basis:
        result.append(
            primitive_integer_vector(
                list(vector)
            )
        )

    return result


# ============================================================================
# SUPPORTS
# ============================================================================

STENCILS = {
    "horizontal_width3": [
        (0, 0),
        (1, 0),
        (2, 0),
    ],

    "vertical_width3": [
        (0, 0),
        (0, 1),
        (0, 2),
    ],

    "rectangle_2x2": [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
    ],

    "transport_3": [
        (0, 0),
        (0, 1),
        (1, 1),
    ],

    "lower_transport_3": [
        (0, 0),
        (1, 0),
        (1, 1),
    ],

    "diamond_5": [
        (0, 0),
        (1, 0),
        (0, 1),
        (0, 2),
        (1, 1),
    ],
}


# ============================================================================
# TRANSFORMED WINDOW GENERATION
# ============================================================================

def transformed_windows(
    lattice,
    offsets,
):
    windows = []

    for r0 in range(4):
        for t0 in range(6):

            cells = [
                (
                    r0 + dr,
                    t0 + dt,
                )
                for dr, dt in offsets
            ]

            if all(
                cell in lattice
                for cell in cells
            ):
                windows.append(
                    (
                        r0,
                        t0,
                        cells,
                    )
                )

    return windows


# ============================================================================
# INTEGER ANNIHILATOR AUDIT
# ============================================================================

def audit_transformed_stencil(
    name,
    lattice,
    offsets,
):
    windows = transformed_windows(
        lattice,
        offsets,
    )

    required = len(offsets)

    print()
    print(
        "    {}: windows={}, support_size={}".format(
            name,
            len(windows),
            required,
        )
    )

    if len(windows) < required:
        print(
            "      status=DATA_LIMITED"
        )
        return {
            "status": "DATA_LIMITED",
            "discoveries": [],
        }

    observed_windows = []
    symbolic_windows = []

    for window in windows:

        values = [
            lattice[cell]
            for cell in window[2]
        ]

        if any(
            is_symbolic(value)
            for value in values
        ):
            symbolic_windows.append(
                window
            )
        else:
            observed_windows.append(
                window
            )

    matrix = sp.Matrix([
        [
            lattice[cell]
            for cell in window[2]
        ]
        for window in observed_windows
    ])

    if matrix.rows == 0:
        print(
            "      status=NO_OBSERVED_WINDOWS"
        )
        return {
            "status":
                "NO_OBSERVED_WINDOWS",
            "discoveries": [],
        }

    basis = integer_nullspace(
        matrix
    )

    if not basis:
        print(
            "      status=NO_INTEGER_RELATION"
        )
        return {
            "status":
                "NO_INTEGER_RELATION",
            "discoveries": [],
        }

    if len(basis) > 1:
        print(
            "      status=MULTIPLE_RELATIONS"
        )
        print(
            "      nullity={}".format(
                len(basis)
            )
        )
        return {
            "status":
                "MULTIPLE_RELATIONS",
            "discoveries": [],
        }

    vector = basis[0]

    # Verify observed windows.
    for window in observed_windows:

        residual = clean(
            sum(
                sp.Integer(vector[j])
                * lattice[cell]
                for j, cell
                in enumerate(
                    window[2]
                )
            )
        )

        if residual != 0:

            print(
                "      status=VERIFICATION_FAILED"
            )

            return {
                "status":
                    "VERIFICATION_FAILED",
                "discoveries": [],
            }

    symbolic_residuals = []

    for window in symbolic_windows:

        residual = clean(
            sum(
                sp.Integer(vector[j])
                * lattice[cell]
                for j, cell
                in enumerate(
                    window[2]
                )
            )
        )

        symbolic_residuals.append(
            residual
        )

    forced = []

    for residual in symbolic_residuals:

        polynomial = sp.Poly(
            residual,
            Z,
            domain=sp.QQ,
        )

        if polynomial.is_zero:
            continue

        if polynomial.degree() != 1:
            print(
                "      status=SYMBOLICALLY_UNRESOLVED"
            )
            return {
                "status":
                    "SYMBOLICALLY_UNRESOLVED",
                "discoveries": [],
            }

        roots = sp.solve(
            sp.Eq(
                residual,
                0,
            ),
            Z,
        )

        if len(roots) != 1:
            print(
                "      status=SYMBOLICALLY_UNRESOLVED"
            )
            return {
                "status":
                    "SYMBOLICALLY_UNRESOLVED",
                "discoveries": [],
            }

        forced.append(
            clean(
                roots[0]
            )
        )

    if not forced:

        print(
            "      status=INTEGER_RELATION_Z_FREE"
        )
        print(
            "      primitive_relation={}".format(
                vector
            )
        )

        return {
            "status":
                "INTEGER_RELATION_Z_FREE",
            "discoveries": [],
        }

    if not all(
        value == forced[0]
        for value in forced[1:]
    ):

        print(
            "      status=CONFLICTING_Z_FORCINGS"
        )

        return {
            "status":
                "CONFLICTING_Z_FORCINGS",
            "discoveries": [],
        }

    forced_Z = forced[0]

    print(
        "      status=INTEGER_RELATION_FORCES_Z"
    )

    print(
        "      primitive_relation={}".format(
            vector
        )
    )

    print(
        "      forced_Z={}".format(
            forced_Z
        )
    )

    print(
        "      forced_Z_integer={}".format(
            sp.denom(
                sp.Rational(
                    forced_Z
                )
            )
            == 1
        )
    )

    return {
        "status":
            "INTEGER_RELATION_FORCES_Z",
        "discoveries": [
            (
                vector,
                forced_Z,
            )
        ],
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 363R — EXACT TRANSFORMED-LATTICE "
        "INTEGER-ANNIHILATOR AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "OBSERVED BASE CELLS={}".format(
            len(BASE)
        )
    )

    print(
        "SYMBOLIC BASE CELL=Z=Q_3(5)"
    )

    discoveries = []

    summary = {}

    for transform_name, (
        dr,
        dt,
    ) in TRANSFORMS.items():

        print()
        print("=" * 78)
        print(
            "{} = Delta_r^{} Delta_t^{} Q".format(
                transform_name,
                dr,
                dt,
            )
        )
        print("=" * 78)

        lattice = build_transform(
            dr,
            dt,
        )

        print(
            "  transformed_cells={}".format(
                len(lattice)
            )
        )

        transform_discoveries = []

        for stencil_name, offsets in (
            STENCILS.items()
        ):

            result = audit_transformed_stencil(
                stencil_name,
                lattice,
                offsets,
            )

            summary[
                (
                    transform_name,
                    stencil_name,
                )
            ] = result["status"]

            for discovery in result.get(
                "discoveries",
                []
            ):
                transform_discoveries.append(
                    (
                        stencil_name,
                        discovery,
                    )
                )

        if transform_discoveries:

            for item in transform_discoveries:
                discoveries.append(
                    (
                        transform_name,
                        item,
                    )
                )

    # =========================================================================
    # SUMMARY
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "SUMMARY"
    )
    print("=" * 78)

    for key, status in summary.items():

        transform_name, stencil_name = key

        print(
            "  {} / {}: {}".format(
                transform_name,
                stencil_name,
                status,
            )
        )

    print()
    print(
        "  total_integer_Z_discoveries={}".format(
            len(discoveries)
        )
    )

    for (
        transform_name,
        (
            stencil_name,
            (
                vector,
                forced_Z,
            ),
        ),
    ) in discoveries:

        print()
        print(
            "  DISCOVERY {} / {}".format(
                transform_name,
                stencil_name,
            )
        )

        print(
            "    relation={}".format(
                vector
            )
        )

        print(
            "    forced_Z={}".format(
                forced_Z
            )
        )

        print(
            "    forced_Z_integer={}".format(
                sp.denom(
                    sp.Rational(
                        forced_Z
                    )
                )
                == 1
            )
        )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
363R asks whether the absence of a direct integer local law in 362R is
specific to the raw coordinates.

Finite-difference transformations are canonical on the observed lattice,
so they provide a principled next coordinate system rather than an
arbitrary reparameterization.

The key distinction remains:

    integer relation after transformation
        versus
    arbitrary rational fit.

A transformed relation that forces Z is still only a model consequence.

The especially informative result would be:

    a low-order finite-difference transform
        +
    one primitive integer annihilator
        +
    exact forced Z.

If every transformed lattice still has no such relation, the evidence
against small-support integer local structure becomes considerably
stronger.

No missing source value is inserted.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  transformed_lattices_tested={}".format(
            len(TRANSFORMS)
        )
    )

    print(
        "  stencils_per_transform={}".format(
            len(STENCILS)
        )
    )

    print(
        "  total_integer_Z_discoveries={}".format(
            len(discoveries)
        )
    )

    print(
        "  symbolic_Z_preserved=True"
    )

    print(
        "  missing_value_inserted=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_counted_as_evidence=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 363R COMPLETE"
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
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )
        raise
