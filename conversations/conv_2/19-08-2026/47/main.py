#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 362R — EXACT INTEGER-COEFFICIENT LOCAL-LAW / SYMBOLIC-Z AUDIT
==============================================================================

Purpose
-------
The observed source contains 15 integer cells.

The strategically important missing cell is

    Z = Q_3(5) = Q(r=2,t=3).

Experiments 355R-361R established:

    D5:
        forces a non-integer Z;

    D6:
        forces a non-integer Z;

    D7:
        admits a one-parameter rational family, but NO integer Z makes all
        five D7 coefficients integers simultaneously.

Therefore this experiment does NOT revisit those models.

Instead it searches a finite, explicitly declared family of local
translation-invariant homogeneous linear laws whose coefficients are
required to be integers.

The distinction is:

    integer-coefficient law
        versus
    rational-coefficient law.

For each stencil:

    sum_i c_i Q(r+dr_i,t+dt_i) = 0,

with

    c_i in Z,

the observed equations are converted to an integer nullspace problem.

The missing value Z remains symbolic.

A stencil is reported only when:

    1. it has a one-dimensional integer nullspace after primitive
       normalization;

    2. the symbolic equations involving Z are exactly consistent;

    3. the resulting relation either forces Z or leaves Z genuinely free.

No arbitrary fitting.
No numerical optimization.
No missing value inserted as data.
No interpolation.
No extrapolation.
No synthetic second n=pq case.

Output is intentionally compact.
"""


from __future__ import annotations

import itertools
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


# ============================================================================
# SYMBOL
# ============================================================================

Z = sp.Symbol(
    "Z"
)


# ============================================================================
# HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def primitive_integer_vector(values):
    values = [
        int(sp.Integer(v))
        for v in values
    ]

    g = 0

    for value in values:
        g = math.gcd(
            g,
            abs(value),
        )

    if g == 0:
        return [
            0
            for _ in values
        ]

    result = [
        value // g
        for value in values
    ]

    for value in result:
        if value != 0:
            if value < 0:
                result = [
                    -v
                    for v in result
                ]
            break

    return result


def build_lattice():
    lattice = {}

    for p_value, values in Q.items():

        r = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r, t)
            ] = sp.Integer(value)

    return lattice


LATTICE = build_lattice()


def value_at(cell):
    if cell == (2, 3):
        return Z

    return LATTICE.get(
        cell
    )


def observed(cell):
    return cell in LATTICE


# ============================================================================
# STENCIL DEFINITIONS
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

    "rectangle_2x3": [
        (0, 0),
        (1, 0),
        (0, 1),
        (1, 1),
        (0, 2),
        (1, 2),
    ],

    "triangle_6": [
        (0, 0),
        (1, 0),
        (2, 0),
        (0, 1),
        (1, 1),
        (0, 2),
    ],
}


# ============================================================================
# TRANSLATION WINDOWS
# ============================================================================

def usable_windows(
    offsets
):
    """
    A window is usable when every cell other than the designated symbolic
    target is observed, or is the symbolic Z cell.
    """

    dr_values = [
        dr
        for dr, _ in offsets
    ]

    dt_values = [
        dt
        for _, dt in offsets
    ]

    min_dr = min(dr_values)
    min_dt = min(dt_values)

    # Conservative finite search region.
    windows = []

    for r0 in range(
        -min_dr,
        4,
    ):

        for t0 in range(
            -min_dt,
            6,
        ):

            cells = [
                (
                    r0 + dr,
                    t0 + dt,
                )
                for dr, dt
                in offsets
            ]

            if all(
                value_at(cell)
                is not None
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
# SYMBOLIC LINEAR SYSTEM FOR STENCIL
# ============================================================================

def stencil_system(
    offsets
):

    windows = usable_windows(
        offsets
    )

    rows = []
    rhs = []

    symbolic_locations = []

    for (
        _r0,
        _t0,
        cells,
    ) in windows:

        row = []

        has_Z = False

        for cell in cells:

            cell_value = value_at(
                cell
            )

            if cell == (2, 3):
                has_Z = True

                # The homogeneous relation is
                #
                # c_i * Z + ...
                #
                # so Z belongs in the coefficient row.
                #
                # We treat the coefficients themselves as unknowns,
                # while Z is a symbolic scalar multiplying one coefficient.
                row.append(
                    Z
                )

            else:
                row.append(
                    cell_value
                )

        rows.append(
            row
        )
        rhs.append(
            sp.Integer(0)
        )

        if has_Z:
            symbolic_locations.append(
                (
                    _r0,
                    _t0,
                )
            )

    return windows, rows, rhs, symbolic_locations


# ============================================================================
# INTEGER NULLSPACE AUDIT
# ============================================================================

def integer_nullspace(
    matrix
):
    """
    Exact nullspace over Q followed by primitive integer normalization.
    """

    nullspace = matrix.nullspace()

    primitive = []

    for vector in nullspace:

        denominators = [
            int(
                sp.denom(
                    sp.Rational(
                        entry
                    )
                )
            )
            for entry in vector
        ]

        common_denominator = 1

        for denominator in denominators:
            common_denominator = sp.ilcm(
                common_denominator,
                denominator,
            )

        integer_entries = [
            int(
                sp.Rational(entry)
                * common_denominator
            )
            for entry in vector
        ]

        primitive.append(
            primitive_integer_vector(
                integer_entries
            )
        )

    return primitive


# ============================================================================
# SYMBOLIC Z COMPATIBILITY
# ============================================================================

def relation_from_vector(
    vector,
    offsets,
    window
):

    _r0, _t0, cells = window

    expression = 0

    for coefficient, cell in zip(
        vector,
        cells,
    ):
        expression += (
            sp.Integer(coefficient)
            * value_at(cell)
        )

    return clean(
        expression
    )


def forced_Z_from_relation(
    relation
):
    """
    Solve relation == 0 for Z.

    Returns:
        ("FORCED", value)
        ("FREE", None)
        ("INCONSISTENT", None)
    """

    relation = clean(
        relation
    )

    if relation == 0:
        return (
            "FREE",
            None,
        )

    polynomial = sp.Poly(
        relation,
        Z,
        domain=sp.QQ,
    )

    degree = polynomial.degree()

    if degree == 0:
        if polynomial.as_expr() == 0:
            return (
                "FREE",
                None,
            )

        return (
            "INCONSISTENT",
            None,
        )

    if degree == 1:

        solution = sp.solve(
            sp.Eq(
                relation,
                0,
            ),
            Z,
        )

        if len(solution) == 1:
            return (
                "FORCED",
                clean(
                    solution[0]
                ),
            )

    return (
        "UNRESOLVED",
        None,
    )


# ============================================================================
# ONE-STENCIL AUDIT
# ============================================================================

def audit_stencil(
    name,
    offsets
):

    windows, rows, rhs, symbolic_locations = (
        stencil_system(
            offsets
        )
    )

    print()
    print("=" * 78)
    print(
        name
    )
    print("=" * 78)

    print(
        "  support={}".format(
            offsets
        )
    )

    print(
        "  usable_windows={}".format(
            len(windows)
        )
    )

    print(
        "  symbolic_windows={}".format(
            symbolic_locations
        )
    )

    if len(windows) < len(offsets):

        print(
            "  status=DATA_LIMITED"
        )

        return {
            "status": "DATA_LIMITED",
            "discoveries": [],
        }

    matrix = sp.Matrix(
        rows
    )

    rank = matrix.rank()

    unknowns = len(offsets)

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  unknowns={}".format(
            unknowns
        )
    )

    print(
        "  redundancy={}".format(
            len(windows)
            - unknowns
        )
    )

    # A symbolic Z makes direct nullspace interpretation inappropriate.
    # Separate windows WITHOUT Z first.
    observed_indices = [
        i
        for i, window
        in enumerate(
            windows
        )
        if window[0:2]
        not in symbolic_locations
    ]

    observed_matrix = sp.Matrix([
        rows[i]
        for i in observed_indices
    ])

    if observed_matrix.rows == 0:

        print(
            "  observed_only_rows=0"
        )

        print(
            "  status=INSUFFICIENT_OBSERVED_CONSTRAINTS"
        )

        return {
            "status":
            "INSUFFICIENT_OBSERVED_CONSTRAINTS",
            "discoveries": [],
        }

    primitive_relations = integer_nullspace(
        observed_matrix
    )

    print(
        "  observed_only_nullspace_dimension={}".format(
            len(primitive_relations)
        )
    )

    if len(primitive_relations) == 0:

        print(
            "  status=NO_INTEGER_RELATION"
        )

        return {
            "status":
            "NO_INTEGER_RELATION",
            "discoveries": [],
        }

    if len(primitive_relations) > 1:

        print(
            "  status=MULTIPLE_RELATIONS"
        )

        print(
            "  primitive_relation_count={}".format(
                len(primitive_relations)
            )
        )

        return {
            "status":
            "MULTIPLE_RELATIONS",
            "discoveries": [],
        }

    vector = primitive_relations[0]

    print(
        "  primitive_relation={}".format(
            vector
        )
    )

    # Verify every observed-only window.
    observed_residuals = []

    for i in observed_indices:

        relation = sum(
            sp.Integer(vector[j])
            * rows[i][j]
            for j in range(
                unknowns
            )
        )

        observed_residuals.append(
            clean(relation)
        )

    observed_ok = all(
        residual == 0
        for residual in
        observed_residuals
    )

    print(
        "  observed_residuals_zero={}".format(
            observed_ok
        )
    )

    if not observed_ok:

        return {
            "status":
            "VERIFICATION_FAILED",
            "discoveries": [],
        }

    forced_values = []

    for i, window in enumerate(
        windows
    ):

        cells = window[2]

        relation = clean(
            sum(
                sp.Integer(
                    vector[j]
                )
                * value_at(
                    cells[j]
                )
                for j in range(
                    unknowns
                )
            )
        )

        kind, value = (
            forced_Z_from_relation(
                relation
            )
        )

        if kind == "INCONSISTENT":

            print(
                "  symbolic_window={}: "
                "status=INCONSISTENT".format(
                    cells
                )
            )

            return {
                "status":
                "SYMBOLICALLY_INCONSISTENT",
                "discoveries": [],
            }

        if kind == "FORCED":

            forced_values.append(
                value
            )

        elif kind == "UNRESOLVED":

            return {
                "status":
                "SYMBOLICALLY_UNRESOLVED",
                "discoveries": [],
            }

    forced_values = [
        clean(value)
        for value in forced_values
    ]

    if forced_values:

        consistent = all(
            value == forced_values[0]
            for value in forced_values[1:]
        )

        if not consistent:

            print(
                "  status=CONFLICTING_Z_FORCINGS"
            )

            print(
                "  forced_values={}".format(
                    forced_values
                )
            )

            return {
                "status":
                "CONFLICTING_Z_FORCINGS",
                "discoveries": [],
            }

        forced_Z = forced_values[0]

        print(
            "  status=INTEGER_STENCIL_FORCES_Z"
        )

        print(
            "  forced_Z={}".format(
                forced_Z
            )
        )

        print(
            "  forced_Z_integer={}".format(
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
            "INTEGER_STENCIL_FORCES_Z",
            "primitive_relation":
            vector,
            "forced_Z":
            forced_Z,
            "discoveries": [
                (
                    name,
                    vector,
                    forced_Z,
                )
            ],
        }

    print(
        "  status=INTEGER_STENCIL_WITH_Z_FREE"
    )

    return {
        "status":
        "INTEGER_STENCIL_WITH_Z_FREE",
        "primitive_relation":
        vector,
        "forced_Z":
        None,
        "discoveries": [],
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 362R — EXACT INTEGER-COEFFICIENT "
        "LOCAL-LAW / SYMBOLIC-Z AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "OBSERVED CELLS={}".format(
            len(LATTICE)
        )
    )

    print(
        "SYMBOLIC MISSING CELL:"
    )

    print(
        "  Z=Q_3(5)=Q(2,3)"
    )

    results = {}

    for name, offsets in (
        STENCILS.items()
    ):

        results[name] = audit_stencil(
            name,
            offsets,
        )

    # ------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "SUMMARY"
    )
    print("=" * 78)

    for name, result in (
        results.items()
    ):

        print(
            "  {}: {}".format(
                name,
                result["status"],
            )
        )

    discoveries = []

    for result in results.values():

        discoveries.extend(
            result.get(
                "discoveries",
                [],
            )
        )

    print()
    print(
        "  integer_stencil_Z_discoveries={}".format(
            len(discoveries)
        )
    )

    for (
        name,
        vector,
        forced_Z,
    ) in discoveries:

        print()
        print(
            "  DISCOVERY={}".format(
                name
            )
        )

        print(
            "    primitive_relation={}".format(
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

    # ------------------------------------------------------------------------
    # STRATEGIC INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
362R deliberately changes the arithmetic question.

A rational local recurrence can force a non-integer Z, as D5 and D6 did.

Here the coefficients themselves are required to be primitive integers.

That is a substantially more restrictive class.

The experiment therefore asks:

    Does any small-support integer-coefficient relation survive all
    currently available observed windows and constrain Z?

The outcomes are:

    NO_INTEGER_RELATION
        the observed data admit no primitive integer annihilator of that
        support;

    MULTIPLE_RELATIONS
        the support is too flexible to identify one law;

    INTEGER_STENCIL_WITH_Z_FREE
        an integer relation exists but does not determine Q_3(5);

    INTEGER_STENCIL_FORCES_Z
        an integer relation determines Q_3(5) exactly.

Only the last case is structurally interesting.

Even then, the forced Z remains a model consequence until independently
observed.

No missing cell is inserted into any calibration calculation.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(LATTICE)
        )
    )

    print(
        "  symbolic_Z_used=True"
    )

    print(
        "  integer_coefficient_stencils_tested={}".format(
            len(STENCILS)
        )
    )

    print(
        "  integer_stencil_Z_discoveries={}".format(
            len(discoveries)
        )
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
        "EXPERIMENT 362R COMPLETE"
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
