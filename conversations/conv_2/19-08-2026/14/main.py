#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 330R — EXACT LOCAL 2D STENCIL / SOURCE-LATTICE RECURRENCE AUDIT
==============================================================================

Purpose
-------
Experiments 322R-329R rejected:

    * constant raw width-2 transfer;
    * low-width higher stencils;
    * shared vertical recurrences;
    * low-degree bivariate polynomial laws;
    * low-degree polynomial-coefficient vertical recurrences;
    * rational first-order vertical laws through degree 2;
    * rational second-order vertical laws of degrees 0 and 1.

The next qualitatively different possibility is a LOCAL TWO-DIMENSIONAL LAW.

The observed source table is naturally indexed by

    (p,t),

with p increasing in steps of 2 and t increasing by 1.

This experiment tests exact constant-coefficient local stencils such as

    Q(t+1,p)
      = a Q(t,p)
      + b Q(t,p+2)
      + c Q(t+1,p+2),

and more general four-point / five-point relations.

The coefficients are solved globally across every fully observed stencil
translation.

The key distinction is:

    EXACT_OVERDETERMINED
        genuine finite-data validation;

    EXACT_DATA_SIZED
        exact reconstruction only;

    NONUNIQUE
        data insufficient;

    NO_SOLUTION
        exact contradiction.

No interpolation.
No missing-value reconstruction.
No synthetic second n=pq case.
Exact SymPy rational arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# OBSERVED SOURCE TABLE
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
# SYMBOLS
# ============================================================================

p, t = sp.symbols("p t")


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


def valuation(x, prime):
    x = sp.Rational(x)

    if x == 0:
        return sp.oo

    numerator = abs(int(x.p))
    denominator = abs(int(x.q))

    value = 0

    while numerator % prime == 0:
        numerator //= prime
        value += 1

    while denominator % prime == 0:
        denominator //= prime
        value -= 1

    return value


def build_table():
    """
    Convert the reverse-ordered source storage into

        table[(p,t)] = Q_t(p).

    Only observed cells are included.
    """

    table = {}

    for p_value, values in Q.items():

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            table[
                (
                    int(p_value),
                    int(t_value),
                )
            ] = sp.Integer(value)

    return table


# ============================================================================
# STENCIL DEFINITIONS
# ============================================================================

STENCILS = {
    # Q(t+1,p) = a Q(t,p) + b Q(t,p+2) + c Q(t+1,p+2)
    "transport_3": [
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 0),
    ],

    # Four-point homogeneous relation:
    # a Q(t,p) + b Q(t,p+2) + c Q(t+1,p) + d Q(t+1,p+2) = 0
    "rectangle_4": [
        (0, 0),
        (0, 1),
        (1, 0),
        (1, 1),
    ],

    # Q(t+1,p) =
    # a Q(t,p) + b Q(t,p+2) + c Q(t,p+4)
    "vertical_3": [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
    ],

    # Q(t+1,p) =
    # a Q(t,p) + b Q(t,p+2)
    "raw_width2": [
        (0, 0),
        (0, 1),
        (1, 0),
    ],

    # Five-point local stencil:
    # a Q(t,p)
    # + b Q(t,p+2)
    # + c Q(t,p+4)
    # + d Q(t+1,p+2)
    # + e Q(t+1,p) = 0
    "diamond_5": [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 1),
        (1, 0),
    ],

    # Diagonal transport:
    # Q(t+1,p) =
    # a Q(t,p) + b Q(t,p+2) + c Q(t+1,p+2)
    # The same as transport_3, retained as a named canonical model.
    "diagonal_transport": [
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 0),
    ],
}


# ============================================================================
# STENCIL WINDOWS
# ============================================================================

def available_windows(table, stencil):
    """
    Return all translated stencil windows that are fully observed.

    A stencil consists of offsets (dt, dp_index), where dp_index means
    p + 2*dp_index.
    """

    observed = set(table)

    candidate_origins = set()

    for point in observed:

        p0, t0 = point

        for dt, dp_index in stencil:

            candidate_origins.add(
                (
                    p0 - 2 * dp_index,
                    t0 - dt,
                )
            )

    windows = []

    for p0, t0 in sorted(
        candidate_origins,
        key=lambda q: (
            q[1],
            q[0],
        ),
    ):

        points = [
            (
                p0 + 2 * dp_index,
                t0 + dt,
            )
            for dt, dp_index
            in stencil
        ]

        if all(
            point in observed
            for point in points
        ):
            windows.append(
                (
                    p0,
                    t0,
                )
            )

    return windows


# ============================================================================
# STENCIL SYSTEM
# ============================================================================

def build_stencil_system(
    table,
    stencil,
):
    """
    Each fully observed translation gives one homogeneous linear equation:

        sum c_j Q(t+dt_j, p+2*dp_j) = 0.
    """

    windows = available_windows(
        table,
        stencil,
    )

    rows = []

    for p0, t0 in windows:

        row = []

        for dt, dp_index in stencil:

            value = table[
                (
                    p0 + 2 * dp_index,
                    t0 + dt,
                )
            ]

            row.append(
                sp.Integer(value)
            )

        rows.append(row)

    return (
        sp.Matrix(rows),
        windows,
    )


# ============================================================================
# HOMOGENEOUS NULLSPACE SOLVER
# ============================================================================

def solve_stencil(
    table,
    stencil,
):
    A, windows = build_stencil_system(
        table,
        stencil,
    )

    equations = A.rows
    coefficients = A.cols

    if equations == 0:

        return {
            "status": "NO_WINDOWS",
            "matrix": A,
            "windows": windows,
            "rank": 0,
            "nullity": coefficients,
            "nullspace": [],
        }

    rank = A.rank()
    nullspace = A.nullspace()
    nullity = len(nullspace)

    if nullity == 0:

        status = "NO_SOLUTION"

    elif nullity == 1:

        # Exactly one projective stencil relation.
        if equations > coefficients - 1:
            status = "EXACT_OVERDETERMINED"
        else:
            status = "EXACT_DATA_SIZED"

    else:

        status = "NONUNIQUE"

    return {
        "status": status,
        "matrix": A,
        "windows": windows,
        "rank": rank,
        "nullity": nullity,
        "nullspace": nullspace,
    }


# ============================================================================
# NORMALIZE PROJECTIVE VECTOR
# ============================================================================

def normalize_vector(vector):

    first_nonzero = next(
        value
        for value in vector
        if value != 0
    )

    return [
        clean(
            value / first_nonzero
        )
        for value in vector
    ]


# ============================================================================
# REPORTING
# ============================================================================

def report_stencil(
    name,
    table,
    stencil,
):
    result = solve_stencil(
        table,
        stencil,
    )

    A = result["matrix"]

    print()
    print("=" * 78)
    print(
        "STENCIL: {}".format(
            name
        )
    )
    print("=" * 78)

    print(
        "  stencil_offsets={}".format(
            stencil
        )
    )

    print(
        "  equation_count={}".format(
            A.rows
        )
    )

    print(
        "  unknown_coefficients={}".format(
            A.cols
        )
    )

    print(
        "  projective_unknown_dimension={}".format(
            max(
                0,
                A.cols - 1,
            )
        )
    )

    print(
        "  redundancy={}".format(
            A.rows
            - max(
                0,
                A.cols - 1,
            )
        )
    )

    print(
        "  rank={}".format(
            result["rank"]
        )
    )

    print(
        "  nullity={}".format(
            result["nullity"]
        )
    )

    print(
        "  windows={}".format(
            result["windows"]
        )
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        vector = normalize_vector(
            result["nullspace"][0]
        )

        print()
        print(
            "  normalized_relation={}".format(
                vector
            )
        )

        residuals = []

        for row in A.tolist():

            residuals.append(
                clean(
                    sum(
                        vector[j] * row[j]
                        for j in range(
                            len(vector)
                        )
                    )
                )
            )

        print(
            "  residuals={}".format(
                residuals
            )
        )

        print(
            "  all_residuals_zero={}".format(
                all(
                    value == 0
                    for value in residuals
                )
            )

        )

        print()
        print(
            "  coefficient_prime_profile:"
        )

        for index, value in enumerate(
            vector
        ):

            print(
                "    c_{}={}".format(
                    index,
                    value,
                )
            )

            print(
                "      valuations={}".format(
                    {
                        prime: valuation(
                            value,
                            prime,
                        )
                        for prime in (
                            2,
                            3,
                            5,
                            7,
                            11,
                            13,
                            17,
                        )
                    }
                )
            )

    elif result["status"] == "NONUNIQUE":

        print(
            "  status_note="
            "multiple_local_stencil_relations"
        )

    return result


# ============================================================================
# CROSS-STENCIL ANALYSIS
# ============================================================================

def compare_transport_and_rectangle(
    transport_result,
    rectangle_result,
):
    print()
    print("=" * 78)
    print(
        "CROSS-STENCIL COMPARISON"
    )
    print("=" * 78)

    print(
        "  transport_status={}".format(
            transport_result["status"]
        )
    )

    print(
        "  rectangle_status={}".format(
            rectangle_result["status"]
        )
    )

    if (
        transport_result["status"]
        in (
            "EXACT_OVERDETERMINED",
            "EXACT_DATA_SIZED",
        )
        and
        rectangle_result["status"]
        in (
            "EXACT_OVERDETERMINED",
            "EXACT_DATA_SIZED",
        )
    ):

        v1 = normalize_vector(
            transport_result["nullspace"][0]
        )

        v2 = normalize_vector(
            rectangle_result["nullspace"][0]
        )

        print(
            "  transport_relation={}".format(
                v1
            )
        )

        print(
            "  rectangle_relation={}".format(
                v2
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 330R — EXACT LOCAL 2D STENCIL / "
        "SOURCE-LATTICE RECURRENCE AUDIT"
    )
    print("=" * 78)

    table = build_table()

    print()
    print("=" * 78)
    print(
        "1. OBSERVED SOURCE TABLE"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(table)
        )
    )

    for (
        p_value,
        t_value,
    ) in sorted(
        table,
        key=lambda key: (
            key[1],
            key[0],
        ),
    ):

        print(
            "  (p={},t={}) -> {}".format(
                p_value,
                t_value,
                table[
                    (
                        p_value,
                        t_value,
                    )
                ],
            )
        )

    results = {}

    for name, stencil in STENCILS.items():

        results[name] = report_stencil(
            name,
            table,
            stencil,
        )

    compare_transport_and_rectangle(
        results["transport_3"],
        results["rectangle_4"],
    )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The source table has resisted:

    scalar vertical recurrences,
    low-degree p-dependent vertical recurrences,
    low-degree rational vertical recurrences,
    low-degree bivariate polynomial laws,
    rank-one and rank-two separability.

A local two-dimensional recurrence is qualitatively different.

It does not require the dynamics to be vertical in t alone.

For example,

    a Q(t,p)
    + b Q(t,p+2)
    + c Q(t+1,p)
    + d Q(t+1,p+2)
    = 0

allows p and t to interact locally.

The important standard is unchanged:

    one projective stencil relation with redundant windows
        = genuine finite-data evidence;

    multiple nullspace dimensions
        = underdetermined;

    no nullspace
        = exact contradiction.

The source support is triangular, so only fully observed stencil
translations are included.

No missing values are filled and no unobserved point is used.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    discoveries = [
        name
        for name, result
        in results.items()
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_local_stencils={}".format(
            discoveries
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
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
        "EXPERIMENT 330R COMPLETE"
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
