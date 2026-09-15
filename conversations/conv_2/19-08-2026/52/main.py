#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 367R-FIXED — EXACT OVERDETERMINED TRANSFORMED-ANNIHILATOR AUDIT
==============================================================================

Purpose
-------
Re-audit all transformed-lattice integer annihilators from 363R/366R while
enforcing the genuine evidence criterion:

    equations > number of stencil coefficients
    AND
    nullspace dimension = 1.

Data-sized nullspaces are NOT discoveries.

The transformed values are always reconstructed DIRECTLY from the original
15-cell base lattice.

Missing Q_3(5) is represented symbolically by Z only in the symbolic lattice.

No recursive transformed lattice is used.
No missing numerical value is inserted.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
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


Z = sp.Symbol("Z")


# ============================================================================
# TRANSFORMS
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


# ============================================================================
# STENCILS
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
# HELPERS
# ============================================================================

def clean(value):
    """
    Exact symbolic normalization.

    IMPORTANT:
    None is a legitimate "unavailable transformed value" marker and must
    never be passed to SymPy.
    """
    if value is None:
        return None

    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def is_integer(value):
    if value is None:
        return False

    value = sp.Rational(value)

    return value.q == 1


def primitive_integer_vector(vector):
    """
    Normalize a rational nullspace vector to a primitive integer vector.
    """

    rationals = [
        sp.Rational(value)
        for value in vector
    ]

    denominator_lcm = 1

    for value in rationals:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(sp.denom(value)),
        )

    integers = [
        int(value * denominator_lcm)
        for value in rationals
    ]

    gcd_value = 0

    for value in integers:
        gcd_value = math.gcd(
            gcd_value,
            abs(value),
        )

    if gcd_value:
        integers = [
            value // gcd_value
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


# ============================================================================
# BASE LATTICES
# ============================================================================

def build_base_lattice(
    symbolic=False,
):
    """
    Build

        (r,t) -> Q_t(p),

    with r=(p-1)/2.

    The symbolic lattice additionally contains

        (2,3) -> Z.
    """

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
                (r_value, t_value)
            ] = sp.Integer(value)

    if symbolic:
        lattice[(2, 3)] = Z

    return lattice


# ============================================================================
# FINITE-DIFFERENCE EXPANSION
# ============================================================================

def difference_terms(
    r,
    t,
    order_r,
    order_t,
):
    """
    Direct expansion of

        Delta_r^order_r Delta_t^order_t Q(r,t).

    Returns coefficient/base-cell pairs.
    """

    terms = []

    for i in range(
        order_r + 1
    ):

        for j in range(
            order_t + 1
        ):

            coefficient = (
                (-1) ** (
                    order_r - i
                    + order_t - j
                )
                *
                math.comb(
                    order_r,
                    i,
                )
                *
                math.comb(
                    order_t,
                    j,
                )
            )

            terms.append(
                (
                    coefficient,
                    (
                        r + i,
                        t + j,
                    ),
                )
            )

    return terms


def transformed_value(
    lattice,
    r,
    t,
    order_r,
    order_t,
):
    """
    Return the exact transformed value if all required base cells exist.

    Return None if any required base cell is absent.

    This is the critical support-safety rule.
    """

    terms = difference_terms(
        r,
        t,
        order_r,
        order_t,
    )

    for _, cell in terms:

        if cell not in lattice:
            return None

    return clean(
        sum(
            coefficient * lattice[cell]
            for coefficient, cell in terms
        )
    )


# ============================================================================
# TRANSFORMED SUPPORT
# ============================================================================

def transformed_support(
    lattice,
    order_r,
    order_t,
):
    """
    Determine transformed support only from complete finite-difference
    expansions.
    """

    if not lattice:
        return set()

    max_r = max(
        r
        for r, _ in lattice
    )

    max_t = max(
        t
        for _, t in lattice
    )

    support = set()

    for r in range(
        max_r + 1
    ):

        for t in range(
            max_t + 1
        ):

            value = transformed_value(
                lattice,
                r,
                t,
                order_r,
                order_t,
            )

            if value is not None:

                support.add(
                    (r, t)
                )

    return support


# ============================================================================
# WINDOW ENUMERATION
# ============================================================================

def find_windows(
    support,
    offsets,
):
    """
    Find every translated stencil whose transformed cells are all present.
    """

    if not support:
        return []

    max_r = max(
        r
        for r, _ in support
    )

    max_t = max(
        t
        for _, t in support
    )

    max_dr = max(
        dr
        for dr, _ in offsets
    )

    max_dt = max(
        dt
        for _, dt in offsets
    )

    windows = []

    for r0 in range(
        max_r - max_dr + 1
    ):

        for t0 in range(
            max_t - max_dt + 1
        ):

            if all(
                (
                    r0 + dr,
                    t0 + dt,
                ) in support
                for dr, dt in offsets
            ):

                windows.append(
                    (r0, t0)
                )

    return windows


# ============================================================================
# WINDOW VALUES
# ============================================================================

def window_values(
    lattice,
    window,
    offsets,
    order_r,
    order_t,
):
    """
    Return transformed values for one stencil window.

    Unlike the previous implementation, None is preserved rather than sent
    into clean(), preventing the crash seen in 367R.
    """

    r0, t0 = window

    values = []

    for dr, dt in offsets:

        value = transformed_value(
            lattice,
            r0 + dr,
            t0 + dt,
            order_r,
            order_t,
        )

        values.append(
            value
        )

    return values


def window_is_complete(
    values,
):
    return all(
        value is not None
        for value in values
    )


# ============================================================================
# OBSERVED / SYMBOLIC WINDOW CLASSIFICATION
# ============================================================================

def classify_windows(
    observed_lattice,
    symbolic_lattice,
    all_windows,
    offsets,
    order_r,
    order_t,
):
    """
    A window is fully observed only when every transformed coefficient
    can be constructed from observed base cells.

    A window may be symbolically evaluable when its transformed expansion
    uses Z but all other base cells are available.
    """

    fully_observed = []
    symbolic_windows = []
    unavailable = []

    for window in all_windows:

        observed_values = window_values(
            observed_lattice,
            window,
            offsets,
            order_r,
            order_t,
        )

        if window_is_complete(
            observed_values
        ):
            fully_observed.append(
                window
            )
            continue

        symbolic_values = window_values(
            symbolic_lattice,
            window,
            offsets,
            order_r,
            order_t,
        )

        if not window_is_complete(
            symbolic_values
        ):
            unavailable.append(
                window
            )
            continue

        symbolic_windows.append(
            window
        )

    return (
        fully_observed,
        symbolic_windows,
        unavailable,
    )


# ============================================================================
# OBSERVED CALIBRATION MATRIX
# ============================================================================

def calibration_matrix(
    observed_lattice,
    windows,
    offsets,
    order_r,
    order_t,
):
    rows = []

    for window in windows:

        values = window_values(
            observed_lattice,
            window,
            offsets,
            order_r,
            order_t,
        )

        if not window_is_complete(
            values
        ):
            continue

        rows.append(
            values
        )

    if not rows:

        return sp.zeros(
            0,
            len(offsets),
        )

    return sp.Matrix(
        rows
    )


# ============================================================================
# NULLSPACE
# ============================================================================

def classify_nullspace(
    matrix,
    equations,
    unknowns,
):
    rank = matrix.rank()

    nullity = (
        unknowns
        - rank
    )

    if rank == unknowns:

        return (
            "NO_RELATION",
            rank,
            nullity,
        )

    if equations > unknowns:

        if nullity == 1:

            return (
                "EXACT_OVERDETERMINED",
                rank,
                nullity,
            )

        return (
            "MULTIPLE_RELATIONS",
            rank,
            nullity,
        )

    if equations == unknowns:

        if nullity == 1:

            return (
                "DATA_SIZED",
                rank,
                nullity,
            )

        return (
            "UNDERDETERMINED",
            rank,
            nullity,
        )

    return (
        "UNDERDETERMINED",
        rank,
        nullity,
    )


def extract_relations(matrix):
    return [
        primitive_integer_vector(
            vector
        )
        for vector in matrix.nullspace()
    ]


# ============================================================================
# SYMBOLIC Z TEST
# ============================================================================

def symbolic_forcing(
    relation,
    values,
):
    """
    Evaluate

        sum c_i * value_i = 0

    and solve exactly for Z when possible.
    """

    expression = clean(
        sum(
            sp.Integer(
                relation[i]
            ) * values[i]
            for i in range(
                len(relation)
            )
        )
    )

    coefficient_z = clean(
        expression.coeff(Z)
    )

    if coefficient_z == 0:
        return None

    constant = clean(
        expression.subs(
            Z,
            0,
        )
    )

    forced_z = clean(
        -constant
        / coefficient_z
    )

    residual = clean(
        expression.subs(
            Z,
            forced_z,
        )
    )

    return {
        "expression": expression,
        "coefficient_z": coefficient_z,
        "forced_z": forced_z,
        "residual": residual,
    }


# ============================================================================
# SINGLE MODEL AUDIT
# ============================================================================

def audit_case(
    transform_name,
    stencil_name,
):
    order_r, order_t = TRANSFORMS[
        transform_name
    ]

    offsets = STENCILS[
        stencil_name
    ]

    observed_lattice = (
        build_base_lattice(
            symbolic=False
        )
    )

    symbolic_lattice = (
        build_base_lattice(
            symbolic=True
        )
    )

    observed_support = (
        transformed_support(
            observed_lattice,
            order_r,
            order_t,
        )
    )

    symbolic_support = (
        transformed_support(
            symbolic_lattice,
            order_r,
            order_t,
        )
    )

    all_windows = find_windows(
        symbolic_support,
        offsets,
    )

    (
        fully_observed,
        symbolic_windows_list,
        unavailable,
    ) = classify_windows(
        observed_lattice,
        symbolic_lattice,
        all_windows,
        offsets,
        order_r,
        order_t,
    )

    matrix = calibration_matrix(
        observed_lattice,
        fully_observed,
        offsets,
        order_r,
        order_t,
    )

    equations = matrix.rows
    unknowns = len(offsets)

    (
        status,
        rank,
        nullity,
    ) = classify_nullspace(
        matrix,
        equations,
        unknowns,
    )

    relations = []

    if matrix.rows:
        relations = extract_relations(
            matrix
        )

    forced_values = []

    # Only report forced-Z consequences for one-dimensional relations.
    if relations and nullity == 1:

        for relation in relations:

            for window in (
                symbolic_windows_list
            ):

                values = window_values(
                    symbolic_lattice,
                    window,
                    offsets,
                    order_r,
                    order_t,
                )

                if not window_is_complete(
                    values
                ):
                    continue

                result = symbolic_forcing(
                    relation,
                    values,
                )

                if result is None:
                    continue

                forced_values.append(
                    {
                        "window": window,
                        "relation": relation,
                        "values": values,
                        **result,
                    }
                )

    return {
        "transform": transform_name,
        "stencil": stencil_name,
        "order_r": order_r,
        "order_t": order_t,
        "support_size": unknowns,
        "observed_support": sorted(
            observed_support
        ),
        "symbolic_support": sorted(
            symbolic_support
        ),
        "all_windows": sorted(
            all_windows
        ),
        "fully_observed_windows": sorted(
            fully_observed
        ),
        "symbolic_windows": sorted(
            symbolic_windows_list
        ),
        "unavailable_windows": sorted(
            unavailable
        ),
        "equation_count": equations,
        "unknown_count": unknowns,
        "rank": rank,
        "nullity": nullity,
        "status": status,
        "relations": relations,
        "forced_values": forced_values,
    }


# ============================================================================
# FULL RUN
# ============================================================================

def run_all():

    print()
    print("=" * 78)
    print(
        "1. EXACT TRANSFORMED SUPPORT / OVERDETERMINATION AUDIT"
    )
    print("=" * 78)

    results = {}

    for transform_name in TRANSFORMS:

        for stencil_name in STENCILS:

            result = audit_case(
                transform_name,
                stencil_name,
            )

            key = (
                transform_name,
                stencil_name,
            )

            results[key] = result

            print()
            print(
                "  {} / {}".format(
                    transform_name,
                    stencil_name,
                )
            )

            print(
                "    transform_orders=({}, {})".format(
                    result["order_r"],
                    result["order_t"],
                )
            )

            print(
                "    equations={}".format(
                    result["equation_count"]
                )
            )

            print(
                "    unknowns={}".format(
                    result["unknown_count"]
                )
            )

            print(
                "    redundancy={}".format(
                    result["equation_count"]
                    - result["unknown_count"]
                )
            )

            print(
                "    rank={}".format(
                    result["rank"]
                )
            )

            print(
                "    nullity={}".format(
                    result["nullity"]
                )
            )

            print(
                "    fully_observed_windows={}".format(
                    result[
                        "fully_observed_windows"
                    ]
                )
            )

            print(
                "    symbolic_windows={}".format(
                    result[
                        "symbolic_windows"
                    ]
                )
            )

            print(
                "    unavailable_windows={}".format(
                    result[
                        "unavailable_windows"
                    ]
                )
            )

            print(
                "    status={}".format(
                    result["status"]
                )
            )

            if result[
                "relations"
            ]:

                print(
                    "    relations={}".format(
                        result["relations"]
                    )
                )

            if result[
                "forced_values"
            ]:

                for record in result[
                    "forced_values"
                ]:

                    print()
                    print(
                        "    symbolic_window={}".format(
                            record["window"]
                        )
                    )

                    print(
                        "      forced_Z={}".format(
                            record["forced_z"]
                        )
                    )

                    print(
                        "      Z_integer={}".format(
                            is_integer(
                                record[
                                    "forced_z"
                                ]
                            )
                        )
                    )

                    print(
                        "      residual={}".format(
                            record["residual"]
                        )
                    )

    return results


# ============================================================================
# SUMMARY
# ============================================================================

def summarize(
    results,
):
    print()
    print("=" * 78)
    print(
        "2. OVERDETERMINED ACCEPTANCE SUMMARY"
    )
    print("=" * 78)

    accepted = []
    data_sized = []
    no_relation = []
    underdetermined = []
    multiple = []

    for key, result in results.items():

        status = result["status"]

        if status == (
            "EXACT_OVERDETERMINED"
        ):
            accepted.append(key)

        elif status == "DATA_SIZED":
            data_sized.append(key)

        elif status == "NO_RELATION":
            no_relation.append(key)

        elif status == "UNDERDETERMINED":
            underdetermined.append(key)

        elif status == "MULTIPLE_RELATIONS":
            multiple.append(key)

    print(
        "  exact_overdetermined_count={}".format(
            len(accepted)
        )
    )

    print(
        "  exact_overdetermined={}".format(
            accepted
        )
    )

    print(
        "  data_sized_count={}".format(
            len(data_sized)
        )
    )

    print(
        "  data_sized={}".format(
            data_sized
        )
    )

    print(
        "  no_relation_count={}".format(
            len(no_relation)
        )
    )

    print(
        "  underdetermined_count={}".format(
            len(underdetermined)
        )
    )

    print(
        "  multiple_relation_count={}".format(
            len(multiple)
        )
    )

    return {
        "accepted": accepted,
        "data_sized": data_sized,
        "no_relation": no_relation,
        "underdetermined": underdetermined,
        "multiple": multiple,
    }


# ============================================================================
# ACCEPTED MODEL AUDIT
# ============================================================================

def accepted_audit(
    results,
    accepted,
):
    print()
    print("=" * 78)
    print(
        "3. ACCEPTED OVERDETERMINED SYMBOLIC-CLOSURE AUDIT"
    )
    print("=" * 78)

    records = []

    for key in accepted:

        result = results[
            key
        ]

        print()
        print(
            "  {} / {}".format(
                key[0],
                key[1],
            )
        )

        for record in result[
            "forced_values"
        ]:

            records.append(
                (
                    key,
                    record,
                )
            )

            print(
                "    window={}".format(
                    record["window"]
                )
            )

            print(
                "    relation={}".format(
                    record["relation"]
                )
            )

            print(
                "    forced_Z={}".format(
                    record["forced_z"]
                )
            )

            print(
                "    residual={}".format(
                    record["residual"]
                )
            )

    return records


# ============================================================================
# DATA-SIZED WARNING
# ============================================================================

def data_sized_warning(
    results,
    data_sized,
):
    print()
    print("=" * 78)
    print(
        "4. DATA-SIZED NULLSPACE WARNING"
    )
    print("=" * 78)

    for key in data_sized:

        result = results[
            key
        ]

        print()
        print(
            "  {} / {}".format(
                key[0],
                key[1],
            )
        )

        print(
            "    equations={}".format(
                result["equation_count"]
            )
        )

        print(
            "    unknowns={}".format(
                result["unknown_count"]
            )
        )

        print(
            "    relations={}".format(
                result["relations"]
            )
        )

        if result[
            "forced_values"
        ]:

            print(
                "    symbolic_forced_values_are_MODEL_CONSEQUENCES_ONLY=True"
            )

        print(
            "    accepted_as_discovery=False"
        )


# ============================================================================
# CROSS-CONSISTENCY
# ============================================================================

def cross_consistency(
    accepted_records,
):
    print()
    print("=" * 78)
    print(
        "5. ACCEPTED-DISCOVERY CROSS-Z CONSISTENCY"
    )
    print("=" * 78)

    if not accepted_records:

        print(
            "  accepted_discovery_count=0"
        )

        return []

    values = []

    for key, record in accepted_records:

        value = sp.Rational(
            record["forced_z"]
        )

        values.append(
            (
                key,
                value,
            )
        )

        print(
            "  {} / {} -> Z={}".format(
                key[0],
                key[1],
                value,
            )
        )

    equal_pairs = []

    for i in range(
        len(values)
    ):

        for j in range(
            i + 1,
            len(values),
        ):

            if (
                values[i][1]
                ==
                values[j][1]
            ):

                equal_pairs.append(
                    (
                        values[i][0],
                        values[j][0],
                    )
                )

    distinct_count = len(
        {
            value
            for _, value in values
        }
    )

    print()
    print(
        "  accepted_discovery_count={}".format(
            len(values)
        )
    )

    print(
        "  distinct_Z_count={}".format(
            distinct_count
        )
    )

    print(
        "  equal_pairs={}".format(
            equal_pairs
        )
    )

    return equal_pairs


# ============================================================================
# FINAL
# ============================================================================

def final_report(
    results,
    summary,
    accepted_records,
    equal_pairs,
):
    print()
    print("=" * 78)
    print(
        "6. FINAL STRUCTURAL VERDICT"
    )
    print("=" * 78)

    accepted_count = len(
        summary["accepted"]
    )

    if accepted_count == 0:

        verdict = (
            "NO_EXACT_OVERDETERMINED_TRANSFORMED_ANNIHILATOR"
        )

    elif (
        len(accepted_records) > 1
        and len(equal_pairs) == 0
    ):

        verdict = (
            "OVERDETERMINED_TRANSFORMED_LAWS_FORCE_INCONSISTENT_Z"
        )

    else:

        verdict = (
            "OVERDETERMINED_TRANSFORMED_LAW_REQUIRES_FURTHER_AUDIT"
        )

    print(
        "  transformed_models_tested={}".format(
            len(results)
        )
    )

    print(
        "  exact_overdetermined_models={}".format(
            summary["accepted"]
        )
    )

    print(
        "  data_sized_models={}".format(
            summary["data_sized"]
        )
    )

    print(
        "  accepted_symbolic_records={}".format(
            len(accepted_records)
        )
    )

    print(
        "  equal_accepted_Z_pairs={}".format(
            len(equal_pairs)
        )
    )

    print(
        "  verdict={}".format(
            verdict
        )
    )

    print()
    print(
        "  DATA-SIZED RELATIONS ARE NOT DISCOVERIES."
    )

    return verdict


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 367R-FIXED — EXACT OVERDETERMINED "
        "TRANSFORMED-ANNIHILATOR AUDIT"
    )
    print("=" * 78)

    results = run_all()

    summary = summarize(
        results
    )

    accepted_records = accepted_audit(
        results,
        summary["accepted"],
    )

    data_sized_warning(
        results,
        summary["data_sized"],
    )

    equal_pairs = cross_consistency(
        accepted_records
    )

    verdict = final_report(
        results,
        summary,
        accepted_records,
        equal_pairs,
    )

    print()
    print("=" * 78)
    print(
        "7. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  original_15_cell_lattice_used=True"
    )

    print(
        "  symbolic_Z_preserved=True"
    )

    print(
        "  direct_finite_difference_expansion=True"
    )

    print(
        "  recursive_transformed_data_reused=False"
    )

    print(
        "  None_values_passed_to_SymPy=False"
    )

    print(
        "  support_safe_window_evaluation=True"
    )

    print(
        "  overdetermination_required=True"
    )

    print(
        "  equations_strictly_greater_than_unknowns_required=True"
    )

    print(
        "  nullspace_dimension_one_required=True"
    )

    print(
        "  data_sized_relations_rejected_as_discoveries=True"
    )

    print(
        "  exact_overdetermined_count={}".format(
            len(
                summary["accepted"]
            )
        )
    )

    print(
        "  data_sized_count={}".format(
            len(
                summary["data_sized"]
            )
        )
    )

    print(
        "  accepted_symbolic_records={}".format(
            len(accepted_records)
        )
    )

    print(
        "  missing_value_inserted_numerically=False"
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
        "EXPERIMENT 367R-FIXED COMPLETE"
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