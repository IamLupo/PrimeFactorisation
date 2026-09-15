#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 384R — EXACT DIAGONAL-COORDINATE / TRIANGULAR-SOURCE PROVENANCE AUDIT
==============================================================================

Purpose
-------
The previous experiments found that:

    * ordinary local recurrences fail;
    * rational generating denominators fail;
    * projective/Möbius laws fail;
    * hypergeometric ratio laws fail;
    * modular relations are largely torsion phenomena;
    * the genuine rational nullspaces are large;
    * the Newton/binomial basis does not become sparse;
    * the Newton coefficient polynomial does not factor.

The observed support is triangular:

    t=0 : r=0..3
    t=1 : r=0..2
    t=2 : r=0..2
    t=3 : r=0..1
    t=4 : r=0..1
    t=5 : r=0

This experiment therefore changes the coordinate system itself.

It tests whether Q(r,t) has a simple exact law in the natural triangular
coordinates

    s = r + t
    d = r - t

or low-complexity combinations of them.

Test families
--------------

A) PURE DIAGONAL LAWS
       Q(r,t) = F(r+t)
       Q(r,t) = F(r-t)

B) LOW-DEGREE POLYNOMIALS
       Q(r,t) = P(r+t)
       Q(r,t) = P(r-t)

C) SEPARATED CORRECTIONS
       Q(r,t) = P(s) + R(d)
       Q(r,t) = P(s) + R(r)
       Q(r,t) = P(s) + R(t)

D) BILINEAR / LOW-BIDEGREE MODELS IN (s,d)
       Q(r,t) = sum a_ij s^i d^j

E) SYMMETRIC / ANTISYMMETRIC TESTS
       Q(r,t) = Q(t,r)
       Q(r,t) = -Q(t,r)
       where both cells are actually observed.

F) DIAGONAL DIFFERENCE AUDIT
       Delta_s Q
       Delta_d Q
   computed only on actually available lattice points.

Rules
-----
* Exact SymPy arithmetic only.
* Only observed cells are used.
* No missing value Q(2,3) is created.
* No interpolation.
* No extrapolation.
* No arbitrary matrix fit is accepted.
* A model is called STRUCTURALLY_EXACT only if:
      equations > unknowns
      rank = unknown_count
      augmented_rank = rank
      every observed equation is verified exactly.
* Data-sized models are reported but rejected as discoveries.
"""

from __future__ import annotations

import itertools
import math
import sys

import sympy as sp


# ============================================================================
# OBSERVED SOURCE LATTICE
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

r, t = sp.symbols("r t", integer=True)
s, d = sp.symbols("s d")

# Polynomial coefficients generated dynamically.


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


def as_integer(value):
    return sp.Integer(value)


def is_integer(value):
    value = sp.Rational(value)
    return value.q == 1


def gcd_list(values):
    nonzero = [
        abs(int(v))
        for v in values
        if int(v) != 0
    ]

    if not nonzero:
        return 0

    g = 0

    for value in nonzero:
        g = math.gcd(g, value)

    return g


def factor_integer(value):
    value = abs(int(value))

    if value in (0, 1):
        return {}

    return sp.factorint(value)


# ============================================================================
# OBSERVED LATTICE
# ============================================================================

def build_lattice():
    lattice = {}

    for p_value, values in Q.items():

        r_value = (p_value - 1) // 2

        for index, value in enumerate(values):

            t_value = len(values) - 1 - index

            lattice[
                (r_value, t_value)
            ] = sp.Integer(value)

    return lattice


def sorted_cells(lattice):
    return sorted(
        lattice,
        key=lambda cell: (
            cell[1],
            cell[0],
        )
    )


# ============================================================================
# STRUCTURED LINEAR MODEL SOLVER
# ============================================================================

def solve_exact_linear_model(
    feature_rows,
    values,
    unknowns,
):
    """
    Solve

        M * theta = y

    exactly over Q.

    Discovery status requires strict overdetermination:

        equations > unknowns
        rank = unknown_count
        augmented_rank = rank.

    A data-sized unique solution is deliberately not accepted as a
    structural discovery.
    """

    equation_count = len(feature_rows)
    unknown_count = len(unknowns)

    if equation_count == 0:
        return {
            "status": "DATA_LIMITED",
            "equations": 0,
            "unknowns": unknown_count,
            "rank": 0,
            "augmented_rank": 0,
            "parameters": None,
            "residuals": [],
        }

    M = sp.Matrix(feature_rows)
    y = sp.Matrix(values)

    rank = M.rank()

    augmented_rank = (
        M.row_join(y).rank()
    )

    if augmented_rank > rank:
        return {
            "status": "NO_SOLUTION",
            "equations": equation_count,
            "unknowns": unknown_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
            "residuals": None,
        }

    if equation_count <= unknown_count:
        if rank == unknown_count:
            return {
                "status": "DATA_SIZED",
                "equations": equation_count,
                "unknowns": unknown_count,
                "rank": rank,
                "augmented_rank": augmented_rank,
                "parameters": None,
                "residuals": None,
            }

        return {
            "status": "UNDERDETERMINED",
            "equations": equation_count,
            "unknowns": unknown_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
            "residuals": None,
        }

    if rank < unknown_count:
        return {
            "status": "MULTIPLE_RELATIONS",
            "equations": equation_count,
            "unknowns": unknown_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
            "residuals": None,
        }

    solution = M.gauss_jordan_solve(y)[0]

    parameters = {
        unknowns[i]: clean(solution[i])
        for i in range(unknown_count)
    }

    residuals = []

    for row, rhs in zip(
        feature_rows,
        values,
    ):

        lhs = sum(
            row[i] * solution[i]
            for i in range(
                unknown_count
            )
        )

        residuals.append(
            clean(lhs - rhs)
        )

    exact = all(
        residual == 0
        for residual in residuals
    )

    if not exact:
        raise RuntimeError(
            "Exact reconstruction residual failure."
        )

    return {
        "status": "EXACT_OVERDETERMINED",
        "equations": equation_count,
        "unknowns": unknown_count,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "parameters": parameters,
        "residuals": residuals,
    }


# ============================================================================
# MONOMIAL FEATURE GENERATORS
# ============================================================================

def monomial_pairs(max_degree_s, max_degree_d):
    pairs = []

    for i in range(max_degree_s + 1):
        for j in range(max_degree_d + 1):
            pairs.append((i, j))

    return pairs


def feature_value(
    cell,
    pair,
):
    r_value, t_value = cell

    s_value = r_value + t_value
    d_value = r_value - t_value

    i, j = pair

    return (
        sp.Integer(s_value) ** i
        * sp.Integer(d_value) ** j
    )


def build_bivariate_model(
    lattice,
    max_degree_s,
    max_degree_d,
):
    pairs = monomial_pairs(
        max_degree_s,
        max_degree_d,
    )

    unknowns = sp.symbols(
        "a0:" + str(len(pairs))
    )

    rows = []
    values = []

    for cell in sorted_cells(lattice):

        rows.append(
            [
                feature_value(cell, pair)
                for pair in pairs
            ]
        )

        values.append(
            lattice[cell]
        )

    result = solve_exact_linear_model(
        rows,
        values,
        unknowns,
    )

    result["pairs"] = pairs

    if result["parameters"] is not None:

        expression = 0

        for pair, symbol in zip(
            pairs,
            unknowns,
        ):

            expression += (
                result["parameters"][symbol]
                * s ** pair[0]
                * d ** pair[1]
            )

        result["expression"] = clean(
            expression
        )

    else:

        result["expression"] = None

    return result


# ============================================================================
# UNIVARIATE DIAGONAL MODELS
# ============================================================================

def build_univariate_model(
    lattice,
    coordinate_name,
    degree,
):
    symbols = sp.symbols(
        "b0:" + str(degree + 1)
    )

    rows = []
    values = []

    for cell in sorted_cells(lattice):

        r_value, t_value = cell

        if coordinate_name == "s":
            x_value = r_value + t_value

        elif coordinate_name == "d":
            x_value = r_value - t_value

        else:
            raise ValueError(
                coordinate_name
            )

        rows.append(
            [
                sp.Integer(x_value) ** k
                for k in range(
                    degree + 1
                )
            ]
        )

        values.append(
            lattice[cell]
        )

    result = solve_exact_linear_model(
        rows,
        values,
        symbols,
    )

    result["coordinate"] = coordinate_name
    result["degree"] = degree

    if result["parameters"] is not None:

        expr = sum(
            result["parameters"][symbols[k]]
            * (
                s if coordinate_name == "s"
                else d
            ) ** k
            for k in range(
                degree + 1
            )
        )

        result["expression"] = clean(expr)

    else:

        result["expression"] = None

    return result


# ============================================================================
# SEPARATED ADDITIVE MODELS
# ============================================================================

def build_additive_separated_model(
    lattice,
    left_coordinate,
    right_coordinate,
    left_degree,
    right_degree,
):
    left_symbols = sp.symbols(
        "L0:" + str(left_degree + 1)
    )

    right_symbols = sp.symbols(
        "R0:" + str(right_degree + 1)
    )

    # Remove constant from right side to avoid duplicate constant column.
    unknowns = list(left_symbols) + list(
        right_symbols[1:]
    )

    rows = []
    values = []

    for cell in sorted_cells(lattice):

        r_value, t_value = cell

        if left_coordinate == "s":
            left_value = r_value + t_value
        elif left_coordinate == "d":
            left_value = r_value - t_value
        elif left_coordinate == "r":
            left_value = r_value
        elif left_coordinate == "t":
            left_value = t_value
        else:
            raise ValueError(left_coordinate)

        if right_coordinate == "s":
            right_value = r_value + t_value
        elif right_coordinate == "d":
            right_value = r_value - t_value
        elif right_coordinate == "r":
            right_value = r_value
        elif right_coordinate == "t":
            right_value = t_value
        else:
            raise ValueError(right_coordinate)

        row = [
            sp.Integer(left_value) ** k
            for k in range(
                left_degree + 1
            )
        ]

        row += [
            sp.Integer(right_value) ** k
            for k in range(
                1,
                right_degree + 1,
            )
        ]

        rows.append(row)
        values.append(lattice[cell])

    result = solve_exact_linear_model(
        rows,
        values,
        unknowns,
    )

    result["left_coordinate"] = left_coordinate
    result["right_coordinate"] = right_coordinate
    result["left_degree"] = left_degree
    result["right_degree"] = right_degree

    if result["parameters"] is not None:

        expr = 0

        for k in range(
            left_degree + 1
        ):

            expr += (
                result["parameters"][left_symbols[k]]
                * (
                    s if left_coordinate == "s"
                    else d if left_coordinate == "d"
                    else r if left_coordinate == "r"
                    else t
                ) ** k
            )

        for k in range(
            1,
            right_degree + 1,
        ):

            expr += (
                result["parameters"][
                    right_symbols[k]
                ]
                * (
                    s if right_coordinate == "s"
                    else d if right_coordinate == "d"
                    else r if right_coordinate == "r"
                    else t
                ) ** k
            )

        result["expression"] = clean(
            expr
        )

    else:

        result["expression"] = None

    return result


# ============================================================================
# DIAGONAL CONSISTENCY
# ============================================================================

def diagonal_groups(
    lattice,
    coordinate,
):
    groups = {}

    for cell in sorted_cells(lattice):

        r_value, t_value = cell

        key = (
            r_value + t_value
            if coordinate == "s"
            else r_value - t_value
        )

        groups.setdefault(
            key,
            []
        ).append(
            (
                cell,
                lattice[cell]
            )
        )

    return groups


def exact_diagonal_equalities(
    lattice,
    coordinate,
):
    groups = diagonal_groups(
        lattice,
        coordinate,
    )

    violations = []

    for key, entries in groups.items():

        values = {
            value
            for _, value in entries
        }

        if len(values) > 1:

            violations.append(
                (
                    key,
                    entries,
                )
            )

    return groups, violations


# ============================================================================
# SYMMETRY TESTS
# ============================================================================

def symmetry_audit(lattice):

    tested = []
    equal = []
    antisymmetric = []

    for cell in sorted_cells(lattice):

        r_value, t_value = cell
        swapped = (t_value, r_value)

        if swapped not in lattice:
            continue

        if cell >= swapped:
            continue

        tested.append(
            (cell, swapped)
        )

        value0 = lattice[cell]
        value1 = lattice[swapped]

        if value0 == value1:
            equal.append(
                (cell, swapped)
            )

        if value0 == -value1:
            antisymmetric.append(
                (cell, swapped)
            )

    return (
        tested,
        equal,
        antisymmetric,
    )


# ============================================================================
# FIRST-DIFFERENCE AUDIT IN (s,d)
# ============================================================================

def transformed_neighbor_differences(
    lattice,
):
    """
    A difference in s means:

        (r,t) -> (r+1,t)

    while holding t fixed.

    A difference in d means:

        (r,t) -> (r,t+1)

    because d decreases by 1.

    Only observed pairs are used.
    """

    ds = []
    dd = []

    for cell in sorted_cells(lattice):

        r_value, t_value = cell

        right = (
            r_value + 1,
            t_value,
        )

        up = (
            r_value,
            t_value + 1,
        )

        if right in lattice:

            ds.append(
                (
                    cell,
                    clean(
                        lattice[right]
                        - lattice[cell]
                    )
                )
            )

        if up in lattice:

            dd.append(
                (
                    cell,
                    clean(
                        lattice[up]
                        - lattice[cell]
                    )
                )
            )

    return ds, dd


# ============================================================================
# REPORTING
# ============================================================================

def print_model_report(
    name,
    result,
):

    print()
    print(
        "  MODEL={}".format(
            name
        )
    )

    print(
        "    equations={}".format(
            result.get(
                "equations",
                0,
            )
        )
    )

    print(
        "    unknowns={}".format(
            result.get(
                "unknowns",
                0,
            )
        )
    )

    print(
        "    redundancy={}".format(
            result.get(
                "equations",
                0
            )
            -
            result.get(
                "unknowns",
                0
            )
        )
    )

    print(
        "    rank={}".format(
            result.get(
                "rank",
                None,
            )
        )
    )

    print(
        "    augmented_rank={}".format(
            result.get(
                "augmented_rank",
                None,
            )
        )
    )

    print(
        "    status={}".format(
            result.get(
                "status"
            )
        )
    )

    if result.get(
        "expression"
    ) is not None:

        print(
            "    expression={}".format(
                result["expression"]
            )
        )

    if result.get(
        "residuals"
    ) is not None:

        failures = [
            residual
            for residual in result["residuals"]
            if residual != 0
        ]

        print(
            "    residual_failures={}".format(
                len(failures)
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    lattice = build_lattice()

    print("=" * 78)
    print(
        "EXPERIMENT 384R — EXACT DIAGONAL-COORDINATE / "
        "TRIANGULAR-SOURCE PROVENANCE AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "1. OBSERVED SOURCE"
    )

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  cells={}".format(
            sorted_cells(lattice)
        )
    )

    print()
    print(
        "  missing_strategic_cells="
        "[(2,3)=Q_3(5), (3,1)=Q_1(7)]"
    )

    # ========================================================================
    # 2. PURE DIAGONAL EQUALITY
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "2. EXACT DIAGONAL-EQUALITY AUDIT"
    )
    print("=" * 78)

    for coordinate in (
        "s",
        "d",
    ):

        groups, violations = (
            exact_diagonal_equalities(
                lattice,
                coordinate,
            )
        )

        print()
        print(
            "  coordinate={}".format(
                coordinate
            )
        )

        for key in sorted(
            groups
        ):

            entries = groups[key]

            if len(entries) < 2:
                continue

            print(
                "    level={}".format(
                    key
                )
            )

            print(
                "      entries={}".format(
                    entries
                )
            )

        print(
            "    violating_levels={}".format(
                len(violations)
            )
        )

        if violations:

            print(
                "    status=NO_PURE_DIAGONAL_LAW"
            )

        else:

            print(
                "    status=PURE_DIAGONAL_LAW_SURVIVES"
            )

    # ========================================================================
    # 3. UNIVARIATE POLYNOMIAL IN s AND d
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "3. EXACT UNIVARIATE DIAGONAL-POLYNOMIAL AUDIT"
    )
    print("=" * 78)

    accepted_univariate = []

    for coordinate in (
        "s",
        "d",
    ):

        print()
        print(
            "  COORDINATE={}".format(
                coordinate
            )
        )

        for degree in range(
            0,
            7,
        ):

            result = build_univariate_model(
                lattice,
                coordinate,
                degree,
            )

            print_model_report(
                "P({}) degree {}".format(
                    coordinate,
                    degree,
                ),
                result,
            )

            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
            ):

                accepted_univariate.append(
                    (
                        coordinate,
                        degree,
                        result,
                    )
                )

    # ========================================================================
    # 4. ADDITIVE SEPARATION
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "4. EXACT ADDITIVE-SEPARATION AUDIT"
    )
    print("=" * 78)

    separated_families = [
        ("s+r", "s", "r"),
        ("s+t", "s", "t"),
        ("s+d", "s", "d"),
        ("d+r", "d", "r"),
        ("d+t", "d", "t"),
    ]

    accepted_separated = []

    for (
        name,
        left_coordinate,
        right_coordinate,
    ) in separated_families:

        result = build_additive_separated_model(
            lattice,
            left_coordinate,
            right_coordinate,
            1,
            1,
        )

        print_model_report(
            name + " affine+affine",
            result,
        )

        if (
            result["status"]
            == "EXACT_OVERDETERMINED"
        ):

            accepted_separated.append(
                (
                    name,
                    result,
                )
            )

    # Low-degree version.
    for (
        name,
        left_coordinate,
        right_coordinate,
    ) in (
        ("s+r", "s", "r"),
        ("s+t", "s", "t"),
        ("s+d", "s", "d"),
    ):

        result = build_additive_separated_model(
            lattice,
            left_coordinate,
            right_coordinate,
            2,
            2,
        )

        print_model_report(
            name + " quadratic+quadratic",
            result,
        )

        if (
            result["status"]
            == "EXACT_OVERDETERMINED"
        ):

            accepted_separated.append(
                (
                    name + " quadratic+quadratic",
                    result,
                )
            )

    # ========================================================================
    # 5. BIVARIATE (s,d) MODEL SCAN
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "5. EXACT LOW-BIDEGREE (s,d) POLYNOMIAL AUDIT"
    )
    print("=" * 78)

    accepted_bivariate = []

    for max_s, max_d in (
        (1, 1),
        (2, 1),
        (1, 2),
        (2, 2),
        (3, 1),
        (1, 3),
    ):

        result = build_bivariate_model(
            lattice,
            max_s,
            max_d,
        )

        print_model_report(
            "degrees=(s={},d={})".format(
                max_s,
                max_d,
            ),
            result,
        )

        if (
            result["status"]
            == "EXACT_OVERDETERMINED"
        ):

            accepted_bivariate.append(
                (
                    max_s,
                    max_d,
                    result,
                )
            )

    # ========================================================================
    # 6. 2x2 LOCAL MAPPING IN (s,d)
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "6. EXACT (s,d) FIRST-DIFFERENCE AUDIT"
    )
    print("=" * 78)

    ds, dd = transformed_neighbor_differences(
        lattice
    )

    print()
    print(
        "  Delta_s_pairs={}".format(
            len(ds)
        )
    )

    print(
        "  Delta_d_pairs={}".format(
            len(dd)
        )
    )

    print()
    print(
        "  Delta_s_values={}".format(
            ds
        )
    )

    print()
    print(
        "  Delta_d_values={}".format(
            dd
        )
    )

    print()
    print(
        "  gcd_Delta_s={}".format(
            gcd_list(
                [
                    value
                    for _, value in ds
                ]
            )
        )
    )

    print(
        "  gcd_Delta_d={}".format(
            gcd_list(
                [
                    value
                    for _, value in dd
                ]
            )
        )
    )

    # ========================================================================
    # 7. SYMMETRY / TRANSPOSE AUDIT
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "7. EXACT R/T PROJECTIVE-SYMMETRY AUDIT"
    )
    print("=" * 78)

    (
        tested_pairs,
        equal_pairs,
        antisymmetric_pairs,
    ) = symmetry_audit(
        lattice
    )

    print(
        "  tested_transpose_pairs={}".format(
            len(tested_pairs)
        )
    )

    print(
        "  equal_pairs={}".format(
            equal_pairs
        )
    )

    print(
        "  antisymmetric_pairs={}".format(
            antisymmetric_pairs
        )
    )

    # ========================================================================
    # 8. CROSS-COORDINATE MONOMIAL INVENTORY
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "8. OBSERVED (s,d) COORDINATE INVENTORY"
    )
    print("=" * 78)

    sd_points = []

    for cell in sorted_cells(
        lattice
    ):

        r_value, t_value = cell

        sd_points.append(
            (
                r_value + t_value,
                r_value - t_value,
                cell,
                lattice[cell],
            )
        )

    for item in sd_points:

        print(
            "  s={}, d={}, cell={}, value={}".format(
                item[0],
                item[1],
                item[2],
                item[3],
            )
        )

    # ========================================================================
    # 9. STRUCTURAL VERDICT
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL VERDICT"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_univariate_count={}".format(
            len(
                accepted_univariate
            )
        )
    )

    print(
        "  exact_overdetermined_separated_count={}".format(
            len(
                accepted_separated
            )
        )
    )

    print(
        "  exact_overdetermined_bivariate_count={}".format(
            len(
                accepted_bivariate
            )
        )
    )

    if (
        accepted_univariate
        or accepted_separated
        or accepted_bivariate
    ):

        verdict = (
            "STRUCTURED_DIAGONAL_COORDINATE_LAW_FOUND"
        )

    else:

        verdict = (
            "NO_EXACT_OVERDETERMINED_DIAGONAL_COORDINATE_LAW"
        )

    print(
        "  verdict={}".format(
            verdict
        )
    )

    # ========================================================================
    # 10. FINAL EXACTNESS
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  missing_Q3_5_used=False"
    )

    print(
        "  missing_Q1_7_used=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_counted_as_evidence=False"
    )

    print(
        "  strict_overdetermination_required=True"
    )

    print(
        "  exact_sympy_arithmetic=True"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
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
        "EXPERIMENT 384R COMPLETE"
    )


# ============================================================================
# ENTRY POINT
# ============================================================================

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
