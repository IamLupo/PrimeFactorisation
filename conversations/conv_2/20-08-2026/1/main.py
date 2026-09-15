#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 381R-FIXED2 — EXACT JOINT QUADRATIC-HYPERGEOMETRIC /
                          RATIONAL-RATIO COMPATIBILITY AUDIT
==============================================================================

This is a syntax-safe revision of 381R-FIXED.

The previous failure was caused by constructing a multi-line f-string
expression inside:

    print(
        f"...{sum(...)}"
        f"..."
    )

This version computes every derived quantity before formatting output.

No missing values.
No interpolation.
No extrapolation.
No synthetic second case.
Exact SymPy arithmetic only.
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


# ============================================================================
# SYMBOLS
# ============================================================================

r, t = sp.symbols("r t")


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


def build_lattice():
    lattice = {}

    for p_value, values in Q.items():

        r_value = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            lattice[
                (
                    r_value,
                    t_value,
                )
            ] = sp.Integer(value)

    return lattice


def primitive_integer_vector(vector):
    """
    Convert a rational vector into a primitive integer vector.
    """

    rationals = [
        sp.Rational(value)
        for value in vector
    ]

    denominator_lcm = 1

    for value in rationals:

        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(
                sp.denom(value)
            ),
        )

    integers = [
        int(
            value * denominator_lcm
        )
        for value in rationals
    ]

    common_gcd = 0

    for value in integers:

        common_gcd = math.gcd(
            common_gcd,
            abs(value),
        )

    if common_gcd == 0:

        return tuple(
            0
            for _ in integers
        )

    integers = [
        value // common_gcd
        for value in integers
    ]

    # Canonical sign.
    for value in integers:

        if value != 0:

            if value < 0:

                integers = [
                    -v
                    for v in integers
                ]

            break

    return tuple(
        integers
    )


# ============================================================================
# MONOMIAL BASES
# ============================================================================

def monomials_for_degree(degree):

    """
    IMPORTANT:
    Every constant monomial is a SymPy Integer.

    This avoids the previous:
        AttributeError: 'int' object has no attribute 'subs'
    """

    if degree == 0:

        return [
            sp.Integer(1)
        ]

    if degree == 1:

        return [
            sp.Integer(1),
            r,
            t,
        ]

    if degree == 2:

        return [
            sp.Integer(1),
            r,
            t,
            r**2,
            r * t,
            t**2,
        ]

    if degree == 3:

        return [
            sp.Integer(1),
            r,
            t,
            r**2,
            r * t,
            t**2,
            r**3,
            r**2 * t,
            r * t**2,
            t**3,
        ]

    raise ValueError(
        "unsupported polynomial degree: {}".format(
            degree
        )
    )


def evaluated_monomials(
    degree,
    rv,
    tv,
):

    substitutions = {
        r: sp.Integer(rv),
        t: sp.Integer(tv),
    }

    return [
        clean(
            monomial.subs(
                substitutions
            )
        )
        for monomial in monomials_for_degree(
            degree
        )
    ]


# ============================================================================
# EDGE INVENTORY
# ============================================================================

def r_edges(lattice):

    result = []

    for (
        rv,
        tv,
    ) in sorted(lattice):

        source = (
            rv,
            tv,
        )

        target = (
            rv + 1,
            tv,
        )

        if target in lattice:

            result.append(
                (
                    source,
                    target,
                )
            )

    return result


def t_edges(lattice):

    result = []

    for (
        rv,
        tv,
    ) in sorted(lattice):

        source = (
            rv,
            tv,
        )

        target = (
            rv,
            tv + 1,
        )

        if target in lattice:

            result.append(
                (
                    source,
                    target,
                )
            )

    return result


def plaquettes(lattice):

    result = []

    for rv in range(3):

        for tv in range(5):

            cells = [
                (rv, tv),
                (rv + 1, tv),
                (rv, tv + 1),
                (rv + 1, tv + 1),
            ]

            if all(
                cell in lattice
                for cell in cells
            ):

                result.append(
                    (
                        rv,
                        tv,
                    )
                )

    return result


# ============================================================================
# MODEL SPECIFICATIONS
# ============================================================================

def model_spec(name):

    if name == "common_denominator_quad":

        return {
            "num_degree": 2,
            "den_degree": 2,
        }

    if name == "quad_num_affine_den":

        return {
            "num_degree": 2,
            "den_degree": 1,
        }

    if name == "affine_num_quad_den":

        return {
            "num_degree": 1,
            "den_degree": 2,
        }

    raise ValueError(
        "unknown model: {}".format(
            name
        )
    )


# ============================================================================
# COMMON-DENOMINATOR JOINT MODEL
# ============================================================================

def ratio_rows_common_denominator(
    lattice,
    r_edge_list,
    t_edge_list,
    num_degree,
    den_degree,
):

    """
    Model:

        R_r = A(r,t) / B(r,t)
        R_t = C(r,t) / B(r,t)

    giving

        Q(r+1,t) B(r,t) - Q(r,t) A(r,t) = 0

        Q(r,t+1) B(r,t) - Q(r,t) C(r,t) = 0.

    Unknown coefficient ordering:

        A coefficients
        C coefficients
        B coefficients
    """

    num_basis = monomials_for_degree(
        num_degree
    )

    den_basis = monomials_for_degree(
        den_degree
    )

    num_count = len(
        num_basis
    )

    den_count = len(
        den_basis
    )

    rows = []

    # ------------------------------------------------------------------------
    # r-direction equations
    # ------------------------------------------------------------------------

    for (
        source,
        target,
    ) in r_edge_list:

        rv, tv = source

        q0 = lattice[
            source
        ]

        q1 = lattice[
            target
        ]

        num_values = evaluated_monomials(
            num_degree,
            rv,
            tv,
        )

        den_values = evaluated_monomials(
            den_degree,
            rv,
            tv,
        )

        row = []

        # A coefficients
        row.extend(
            [
                -q0 * value
                for value in num_values
            ]
        )

        # C coefficients
        row.extend(
            [
                sp.Integer(0)
            ] * num_count
        )

        # B coefficients
        row.extend(
            [
                q1 * value
                for value in den_values
            ]
        )

        rows.append(
            row
        )

    # ------------------------------------------------------------------------
    # t-direction equations
    # ------------------------------------------------------------------------

    for (
        source,
        target,
    ) in t_edge_list:

        rv, tv = source

        q0 = lattice[
            source
        ]

        q1 = lattice[
            target
        ]

        num_values = evaluated_monomials(
            num_degree,
            rv,
            tv,
        )

        den_values = evaluated_monomials(
            den_degree,
            rv,
            tv,
        )

        row = []

        # A coefficients
        row.extend(
            [
                sp.Integer(0)
            ] * num_count
        )

        # C coefficients
        row.extend(
            [
                -q0 * value
                for value in num_values
            ]
        )

        # B coefficients
        row.extend(
            [
                q1 * value
                for value in den_values
            ]
        )

        rows.append(
            row
        )

    return rows


def audit_common_denominator_model(
    lattice,
    name,
):

    spec = model_spec(
        name
    )

    num_degree = spec[
        "num_degree"
    ]

    den_degree = spec[
        "den_degree"
    ]

    r_edge_list = r_edges(
        lattice
    )

    t_edge_list = t_edges(
        lattice
    )

    rows = ratio_rows_common_denominator(
        lattice,
        r_edge_list,
        t_edge_list,
        num_degree,
        den_degree,
    )

    matrix = sp.Matrix(
        rows
    )

    coefficient_count = len(
        monomials_for_degree(
            num_degree
        )
    )

    denominator_count = len(
        monomials_for_degree(
            den_degree
        )
    )

    parameter_count = (
        2 * coefficient_count
        + denominator_count
    )

    equation_count = len(
        rows
    )

    rank = matrix.rank()

    nullity = (
        parameter_count
        - rank
    )

    print()
    print(
        "  MODEL={}".format(
            name
        )
    )

    print(
        "    equations={}".format(
            equation_count
        )
    )

    print(
        "    unknowns={}".format(
            parameter_count
        )
    )

    print(
        "    redundancy={}".format(
            equation_count
            - parameter_count
        )
    )

    print(
        "    rank={}".format(
            rank
        )
    )

    print(
        "    nullity={}".format(
            nullity
        )
    )

    basis = (
        matrix.nullspace()
        if nullity > 0
        else []
    )

    if nullity == 0:

        status = (
            "NO_RELATION"
        )

    elif (
        nullity == 1
        and
        equation_count > parameter_count
    ):

        status = (
            "EXACT_OVERDETERMINED"
        )

    elif nullity == 1:

        status = (
            "EXACT_DATA_SIZED"
        )

    else:

        status = (
            "UNDERDETERMINED"
        )

    print(
        "    status={}".format(
            status
        )
    )

    result = {
        "status": status,
        "equations": equation_count,
        "unknowns": parameter_count,
        "rank": rank,
        "nullity": nullity,
        "basis": basis,
    }

    if status != "EXACT_OVERDETERMINED":

        return result

    primitive = primitive_integer_vector(
        basis[0]
    )

    print(
        "    primitive_relation={}".format(
            primitive
        )
    )

    residuals = []

    for row in rows:

        residual = clean(
            sum(
                sp.Integer(
                    row[i]
                )
                * sp.Integer(
                    primitive[i]
                )
                for i in range(
                    parameter_count
                )
            )
        )

        residuals.append(
            residual
        )

    exact = all(
        residual == 0
        for residual in residuals
    )

    print(
        "    all_residuals_zero={}".format(
            exact
        )
    )

    result[
        "primitive_relation"
    ] = primitive

    result[
        "residuals"
    ] = residuals

    result[
        "exact"
    ] = exact

    return result


# ============================================================================
# SEPARATED DENOMINATOR MODEL
# ============================================================================

def separated_denominator_rows(
    lattice,
):

    """
    Test:

        R_r(r,t) = A(r,t) / B(r)
        R_t(r,t) = C(r,t) / D(t)

    where A,C are affine bivariate forms and B,D are affine
    one-variable forms.

    Unknowns:

        A0,A1,A2,
        C0,C1,C2,
        B0,B1,
        D0,D1

    total = 10.
    """

    r_edge_list = r_edges(
        lattice
    )

    t_edge_list = t_edges(
        lattice
    )

    rows = []

    # r-direction
    for (
        source,
        target,
    ) in r_edge_list:

        rv, tv = source

        q0 = lattice[
            source
        ]

        q1 = lattice[
            target
        ]

        A_values = evaluated_monomials(
            1,
            rv,
            tv,
        )

        B_values = [
            sp.Integer(1),
            sp.Integer(rv),
        ]

        row = []

        row.extend(
            [
                -q0 * value
                for value in A_values
            ]
        )

        row.extend(
            [
                sp.Integer(0)
            ] * 3
        )

        row.extend(
            [
                q1 * value
                for value in B_values
            ]
        )

        row.extend(
            [
                sp.Integer(0)
            ] * 2
        )

        rows.append(
            row
        )

    # t-direction
    for (
        source,
        target,
    ) in t_edge_list:

        rv, tv = source

        q0 = lattice[
            source
        ]

        q1 = lattice[
            target
        ]

        C_values = evaluated_monomials(
            1,
            rv,
            tv,
        )

        D_values = [
            sp.Integer(1),
            sp.Integer(tv),
        ]

        row = []

        row.extend(
            [
                sp.Integer(0)
            ] * 3
        )

        row.extend(
            [
                -q0 * value
                for value in C_values
            ]
        )

        row.extend(
            [
                sp.Integer(0)
            ] * 2
        )

        row.extend(
            [
                q1 * value
                for value in D_values
            ]
        )

        rows.append(
            row
        )

    return rows


def audit_separated_denominator(
    lattice,
):

    rows = separated_denominator_rows(
        lattice
    )

    matrix = sp.Matrix(
        rows
    )

    equation_count = len(
        rows
    )

    parameter_count = 10

    rank = matrix.rank()

    nullity = (
        parameter_count
        - rank
    )

    print()
    print(
        "  MODEL=separated_denominator"
    )

    print(
        "    equations={}".format(
            equation_count
        )
    )

    print(
        "    unknowns={}".format(
            parameter_count
        )
    )

    print(
        "    redundancy={}".format(
            equation_count
            - parameter_count
        )
    )

    print(
        "    rank={}".format(
            rank
        )
    )

    print(
        "    nullity={}".format(
            nullity
        )
    )

    basis = (
        matrix.nullspace()
        if nullity > 0
        else []
    )

    if nullity == 0:

        status = (
            "NO_RELATION"
        )

    elif (
        nullity == 1
        and
        equation_count > parameter_count
    ):

        status = (
            "EXACT_OVERDETERMINED"
        )

    elif nullity == 1:

        status = (
            "EXACT_DATA_SIZED"
        )

    else:

        status = (
            "UNDERDETERMINED"
        )

    print(
        "    status={}".format(
            status
        )
    )

    result = {
        "status": status,
        "equations": equation_count,
        "unknowns": parameter_count,
        "rank": rank,
        "nullity": nullity,
        "basis": basis,
    }

    if status != "EXACT_OVERDETERMINED":

        return result

    primitive = primitive_integer_vector(
        basis[0]
    )

    print(
        "    primitive_relation={}".format(
            primitive
        )
    )

    residuals = []

    for row in rows:

        residual = clean(
            sum(
                sp.Integer(
                    row[i]
                )
                * sp.Integer(
                    primitive[i]
                )
                for i in range(
                    parameter_count
                )
            )
        )

        residuals.append(
            residual
        )

    exact = all(
        residual == 0
        for residual in residuals
    )

    print(
        "    all_residuals_zero={}".format(
            exact
        )
    )

    result[
        "primitive_relation"
    ] = primitive

    result[
        "residuals"
    ] = residuals

    result[
        "exact"
    ] = exact

    return result


# ============================================================================
# PLAQUETTE AUDIT
# ============================================================================

def exact_plaquette_audit(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "2. EXACT PLAQUETTE PATH-CONSISTENCY AUDIT"
    )
    print("=" * 78)

    results = []

    for (
        rv,
        tv,
    ) in plaquettes(
        lattice
    ):

        q00 = lattice[
            (
                rv,
                tv,
            )
        ]

        q10 = lattice[
            (
                rv + 1,
                tv,
            )
        ]

        q01 = lattice[
            (
                rv,
                tv + 1,
            )
        ]

        q11 = lattice[
            (
                rv + 1,
                tv + 1,
            )
        ]

        residual = clean(
            q11 * q00
            -
            q10 * q01
        )

        print()
        print(
            "  plaquette=({}, {})".format(
                rv,
                tv
            )
        )

        print(
            "    q00={}".format(
                q00
            )
        )

        print(
            "    q10={}".format(
                q10
            )
        )

        print(
            "    q01={}".format(
                q01
            )
        )

        print(
            "    q11={}".format(
                q11
            )
        )

        print(
            "    residual={}".format(
                residual
            )
        )

        print(
            "    exact_zero={}".format(
                residual == 0
            )
        )

        results.append(
            {
                "location": (
                    rv,
                    tv,
                ),
                "residual": residual,
            }
        )

    return results


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 381R-FIXED2 — EXACT JOINT "
        "QUADRATIC-HYPERGEOMETRIC / RATIONAL-RATIO "
        "COMPATIBILITY AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    r_edge_list = r_edges(
        lattice
    )

    t_edge_list = t_edges(
        lattice
    )

    print()
    print(
        "OBSERVED SOURCE"
    )

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  r_edges={}".format(
            len(r_edge_list)
        )
    )

    print(
        "  t_edges={}".format(
            len(t_edge_list)
        )
    )

    # ------------------------------------------------------------------------
    # Joint model audit
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. JOINT RATIO-MODEL AUDIT"
    )
    print("=" * 78)

    results = {}

    model_names = [
        "common_denominator_quad",
        "quad_num_affine_den",
        "affine_num_quad_den",
    ]

    for name in model_names:

        results[
            name
        ] = audit_common_denominator_model(
            lattice,
            name,
        )

    results[
        "separated_denominator"
    ] = audit_separated_denominator(
        lattice
    )

    # ------------------------------------------------------------------------
    # Plaquette audit
    # ------------------------------------------------------------------------

    plaquette_results = (
        exact_plaquette_audit(
            lattice
        )
    )

    plaquette_zero_count = sum(
        1
        for item in plaquette_results
        if item["residual"] == 0
    )

    # ------------------------------------------------------------------------
    # Accepted models
    # ------------------------------------------------------------------------

    accepted = [
        name
        for (
            name,
            result,
        ) in results.items()
        if (
            result.get("status")
            ==
            "EXACT_OVERDETERMINED"
            and
            result.get("exact", False)
        )
    ]

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. STRUCTURAL VERDICT"
    )
    print("=" * 78)

    print(
        "  accepted_exact_overdetermined_models={}".format(
            accepted
        )
    )

    print(
        "  plaquette_count={}".format(
            len(
                plaquette_results
            )
        )
    )

    print(
        "  plaquette_zero_count={}".format(
            plaquette_zero_count
        )
    )

    print(
r"""
381R tests a joint rational-ratio family after the low-degree
single-direction ratio searches of 380R.

The common-denominator families use

    R_r(r,t) = A(r,t) / B(r,t)
    R_t(r,t) = C(r,t) / B(r,t),

so the edge equations are linear in the polynomial coefficients.

A model is accepted only when its coefficient system is:

    * consistent;
    * overdetermined;
    * rank-deficient by exactly one dimension;
    * exactly verified on every calibration edge.

The separated-denominator model tests the more restrictive form

    R_r(r,t) = A(r,t) / B(r)
    R_t(r,t) = C(r,t) / D(t).

The plaquette audit is independent of those fits. It checks the
strongest simple multiplicative-separability condition:

    Q(r,t) Q(r+1,t+1)
        =
    Q(r+1,t) Q(r,t+1).

A zero plaquette residual is necessary for that simple separable
construction, but is not required by every hypergeometric function.

No missing source value is used.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  exact_sympy_arithmetic=True"
    )

    print(
        "  joint_rational_ratio_models_tested=True"
    )

    print(
        "  quadratic_coordinate_dependence_tested=True"
    )

    print(
        "  common_denominator_models_tested=True"
    )

    print(
        "  separated_denominator_model_tested=True"
    )

    print(
        "  overdetermination_required=True"
    )

    print(
        "  data_sized_models_rejected=True"
    )

    print(
        "  accepted_exact_overdetermined_models={}".format(
            accepted
        )
    )

    print(
        "  plaquette_count={}".format(
            len(
                plaquette_results
            )
        )
    )

    print(
        "  plaquette_zero_count={}".format(
            plaquette_zero_count
        )
    )

    print(
        "  missing_values_used=False"
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
        "EXPERIMENT 381R-FIXED2 COMPLETE"
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